# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 00:46:45 2024

@author: zziya
"""

import networkx as nx
import numpy as np
import scipy.sparse as sp
fpfunc='0'
def adjacency(G):
    A = []
    for i in G.nodes:
        Atemp = []
        for j in G.nodes:
            if G.has_edge(i, j):
                Atemp.append(1)
            else:
                Atemp.append(0)
        A.append(Atemp)
    return A
def Failure_Probability(k):
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
    if fpfunc=='0':
        return 0
    
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

def random_walk_stationary(G):
    n = G.number_of_nodes()
    pi = np.zeros(n)
    adj = adjacency(G)
    W = np.sum(adj)
    for i in range(n):
        pi[i] = np.sum(adj[i]) / W
    return pi

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
                    rownum = i * N + k
                    X[linenum, rownum] = -0.5 * L[j, k]
    x = X.todense()
    tau = np.linalg.solve(x, Y)
    Tau = tau.reshape(N, N)
    return Tau

network_type = 'dolphins'
filename = '.\\datas\\' + network_type + '.txt'

# Note: the function nx.read_edgelist('dolphins.txt') is alternative

G = nx.Graph()
with open(filename) as file:
    for line in file:
        head, tail = [str(x) for x in line.split()]
        G.add_edge(int(head), int(tail))

# relabel
j = 0
new_mapping = {}
for i in G.nodes:
    new_mapping[i] = j
    j += 1
G = nx.relabel_nodes(G, new_mapping)
p = laplacian_with_failure(G)-normalized_laplacian(adjacency(G))
pi=random_walk_with_failure(G)-random_walk_stationary(G)
tau=Find_CT(G)-Find_CT_with_failure(G)