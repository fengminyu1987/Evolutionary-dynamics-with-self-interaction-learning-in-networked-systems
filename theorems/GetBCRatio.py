# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 18:04:52 2024
input: the adjacency matrix
output: coalescence time \tau_{ij}, in matrix form, with size n^2*n^2
@author: zziya
"""
import networkx as nx
import numpy as np
import scipy.sparse as sp
import time
fpfunc='e^-k'#(failure probability for replacement)
realset=['dolphins','americanfootball','contiguous-usa','zebra','ca-sandi_auths','rt-retweet']
fpfuncs=['(k+1)^-1','e^-k','(k+1)^-2','1-k^-1']

def adjacency(G):
    A = []
    for i in G.nodes:
        Atemp = []
        for j in G.nodes:
            if G.has_edge(i, j):
                Atemp.append(G[i][j]['weight'])
            else:
                Atemp.append(0)
        A.append(Atemp)
    return A

def Failure_Vector(G):
    n=G.number_of_nodes()
    FV=[]
    for i in range(n):
        FV.append(1-Failure_Probability(G.degree(i)))
    return np.array(FV)

def Failure_Probability(k, fpfunc):
    if fpfunc=='(k+1)^-1':
        return 1/(k+1)
    if fpfunc=='e^-k':
        return np.exp(-k)
    if fpfunc=='(k+1)^-2':
        return 1/(k+1)**2
    if fpfunc=='1-k^-1':
        return 1-1/k
    if fpfunc=='1-e^-k':
        return 1-np.exp(-k)
    if fpfunc=='1-k^-2':
        return 1-1/(k*k)
    if fpfunc != 'num':
        return float(fpfunc)
    
def normalized_laplacian(adj_matrix):
    n = len(adj_matrix)
    P = np.zeros((n, n))
    for row in range(n):
        sumOverRow = np.sum(adj_matrix[row])
        for col in range(n):
            P[row, col] = adj_matrix[row][col] / sumOverRow
    return P

def laplacian_with_failure(G):#p
    adj_matrix=adjacency(G)
    n = len(adj_matrix)
    X = sp.lil_matrix((n, n))
    for i in range(n):
        for j in range(n):
            ki=G.degree(i)
            delta=Failure_Probability(ki)
            if i==j:
                X[i,j]=delta
            else:
                ineis=G.neighbors(i)
                if j in ineis:
                    X[i,j]=(1-delta)/ki
    return X.todense()

def random_walk_with_failure(G):#pi
    n = G.number_of_nodes()
    transition_matrix=laplacian_with_failure(G)
    transition_matrix_transp = transition_matrix.T
    eigenvals, eigenvects = np.linalg.eig(transition_matrix_transp)
    close_to_1_idx = np.isclose(eigenvals,1)
    target_eigenvect = eigenvects[:,close_to_1_idx]
    target_eigenvect = target_eigenvect[:,0]
    stationary_distribution = target_eigenvect / sum(target_eigenvect)
    result = np.array(np.real(stationary_distribution))
    return np.array([i for arr in result for i in arr])

def Find_CT_with_failure(G):#tau
    adj = adjacency(G)
    #L = normalized_laplacian(adj)
    L = laplacian_with_failure(G)
    N = len(adj)
    X = sp.lil_matrix((N * N, N * N))
    Y = np.ones((N * N, 1))
    for i in range(N):
        for j in range(N):
            if i == j:
                linenum = i * N + j
                X[linenum, linenum] = 1
                Y[linenum] = 0
            else:
                linenum = i * N + j
                X[linenum, linenum] = 1 - 0.5 * L[j,j]
                neis_i = G.neighbors(i)
                neis_j = G.neighbors(j)
                for k in neis_i:
                    rownum = j * N + k
                    X[linenum, rownum] = -0.5 * L[i, k]
                for k in neis_j:
                    rownum = i * N + k
                    X[linenum, rownum] = -0.5 * L[j, k]
                k = i
                rownum = j * N + k
                X[linenum, rownum] = -0.5 * L[i, i]
    x = X.todense()
    tau = np.linalg.solve(x, Y)
    Tau = tau.reshape(N, N)
    return Tau

def Find_CT(G):
    adj = adjacency(G)
    L = normalized_laplacian(adj)
    N = len(adj)
    X = sp.lil_matrix((N * N, N * N))
    Y = np.ones((N * N, 1))
    for i in range(N):
        for j in range(N):
            if i == j:
                linenum = i * N + j
                X[linenum, linenum] = 1
                Y[linenum] = 0
            else:
                linenum = i * N + j
                X[linenum, linenum] = 1
                neis_i = G.neighbors(i)
                neis_j = G.neighbors(j)
                for k in neis_i:
                    rownum = j * N + k
                    X[linenum, rownum] = -0.5 * L[i, k]
                for k in neis_j:
                    if k !=j:
                        rownum = i * N + k
                        X[linenum, rownum] = -0.5 * L[j, k]
                    if k==j:
                        rownum = i * N + k
                        X[linenum, rownum] = 1-0.5 * L[j, k]
    x = X.todense()
    tau = np.linalg.solve(x, Y)
    Tau = tau.reshape(N, N)
    return Tau

def random_walk_stationary(G):
    n = G.number_of_nodes()
    pi = np.zeros(n)
    adj = adjacency(G)
    W = np.sum(adj)
    for i in range(n):
        pi[i] = np.sum(adj[i]) / W
    return pi


def Get_BCRatio(G):
    adj = adjacency(G)
    pi = random_walk_stationary(G)
    p = normalized_laplacian(adj)
    tau = Find_CT(G)
    t1 = np.sum(pi.T @ (p * tau))
    t2 = np.sum(pi.T @ (p @ p * tau))
    t3 = np.sum(pi.T @ (p @ p @ p * tau))
    return (t2) / (t3 - t1)

def Trans_to_Selfloop(g,fpfunc):
    for i in g.nodes:
        g.add_edge(i,i,weight=Failure_Probability(g.degree(i), fpfunc))
    return g
network_type = 'dolphins'
filename = '.\\datas\\' + network_type + '.txt'
# Note: the function nx.read_edgelist('dolphins.txt') is alternative
G = nx.Graph()
with open(filename) as file:
    for line in file:
        head, tail = [str(x) for x in line.split()]
        G.add_edge(int(head), int(tail), weight=1)
# relabel
if 0 not in G.nodes():
    j = 0
    new_mapping = {}
    for i in G.nodes:
        new_mapping[i] = j
        j += 1
    G = nx.relabel_nodes(G, new_mapping)
G=Trans_to_Selfloop(G, fpfunc)
print(Get_BCRatio(G))