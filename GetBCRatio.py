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
fpfunc='lnk'#(failure probability for replacement)
realset=['dolphins','americanfootball','contiguous-usa','zebra','ca-sandi_auths','rt-retweet']
fpfuncs=['(k+1)^-1','e^-k','(k+1)^-2','1-k^-1']
def adjacency(G):
    A = []
    N = len(G.nodes())
    for i in range(N):
        Atemp = []
        for j in range(N):
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
    if fpfunc == '0':
        return 0
    if fpfunc=='lnk':
        return np.log(k)
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
def star(N,alpha, beta):
    numerator = 3 * (1 + alpha) * (-1 + N + beta) * (-(8/3) + N**3 * (4/3 + (11 * alpha)/3 + alpha**2) + (14 * beta)/3 - (4 * beta**2)/3 + alpha**2 * (-2 + (8 * beta)/3) + alpha * (-6 + 10 * beta - (8 * beta**2)/3) + N**2 * (-(16/3) + 3 * beta + alpha**2 * (-4 + (4 * beta)/3) + alpha * (-(40/3) + 7 * beta)) + N * (20/3 + alpha**2 * (5 - 4 * beta) - (23 * beta)/3 + (4 * beta**2)/3 + alpha * (47/3 - 17 * beta + (8 * beta**2)/3)))
    denominator = 2 * (beta * (-N**3 + N * (-1 - beta/2) + N**2 * (2 - beta/2) + beta) + alpha**3 * (2 + N**4 - 5 * beta + 4 * beta**2 + N**3 * (-5 + (5 * beta)/2) + N * (-7 + (25 * beta)/2 - 6 * beta**2) + N**2 * (9 - 10 * beta + 2 * beta**2)) + alpha * (4 + N**4 - 11 * beta + 12 * beta**2 - 2 * beta**3 + N**3 * (-7 + (5 * beta)/2) + N**2 * (15 - 16 * beta + 5 * beta**2) + N * (-13 + (49 * beta)/2 - 17 * beta**2 + 2 * beta**3)) + alpha**2 * (8 + 4 * N**4 - 22 * beta + 21 * beta**2 - 4 * beta**3 + N**3 * (-20 + 12 * beta) + N**2 * (36 - 46 * beta + (27 * beta**2)/2) + N * (-28 + 56 * beta - (69 * beta**2)/2 + 4 * beta**3)))
    return numerator / denominator
def star_HH(N,alpha,gamma):
    numerator = (
        3 * (1 + alpha) * (N + gamma) * (
            N**4 * (alpha**3 + (16/3) * alpha**2 + (25/3) * alpha + (10/3)) +
            N**3 * (
                alpha**3 * ((4/3) * gamma - 2) +
                alpha**2 * (10 * gamma - (25/3)) +
                alpha * (20 * gamma - (34/3)) +
                (29/3) * gamma - (13/3)
            ) +
            N**2 * (
                alpha**3 * (-(7/3) * gamma + (1/3)) +
                alpha**2 * ((11/3) * gamma**2 - (38/3) * gamma + 1) +
                alpha * ((40/3) * gamma**2 - (64/3) * gamma + 1) +
                (26/3) * gamma**2 - 9 * gamma + (1/3)
            ) +
            N * (
                alpha**3 * (-(1/3) * gamma + (2/3)) +
                alpha**2 * (-(10/3) * gamma**2 - gamma + 2) +
                alpha * ((8/3) * gamma**3 - (32/3) * gamma**2 - gamma + 2) +
                3 * gamma**3 - (16/3) * gamma**2 - (1/3) * gamma + (2/3)
            ) +
            gamma * (
                (4/3) * alpha**3 + 4 * alpha**2 +
                alpha * (-(5/3) * gamma**2 + 4) +
                (1/3) * gamma**3 - gamma**2 + (4/3)
            )
        )
    )
    denominator = (
        2 * (
            N**5 * (alpha**4 + (11/2) * alpha**3 + 9 * alpha**2 + 5 * alpha + 2) +
            N**4 * (
                alpha**4 * ((5/2) * gamma - (5/2)) +
                alpha**3 * ((33/2) * gamma - (23/2)) +
                alpha**2 * ((61/2) * gamma - 19) +
                alpha * ((33/2) * gamma - (25/2)) +
                6 * gamma - 4
            ) +
            N**3 * (
                alpha**4 * (2 * gamma**2 - (11/2) * gamma + 2) +
                alpha**3 * ((37/2) * gamma**2 - 29 * gamma + (15/2)) +
                alpha**2 * (44 * gamma**2 - (107/2) * gamma + 11) +
                alpha * (29 * gamma**2 - 36 * gamma + 7) +
                (21/2) * gamma**2 - 12 * gamma + (3/2)
            ) +
            N**2 * (
                alpha**4 * (-(7/2) * gamma**2 + 2 * gamma - (1/2)) +
                alpha**3 * gamma * ((11/2) * gamma**2 - (49/2) * gamma + (15/2)) +
                alpha**2 * ((47/2) * gamma**3 - (111/2) * gamma**2 + (23/2) * gamma + (7/2)) +
                alpha * (21 * gamma**3 - (79/2) * gamma**2 + 7 * gamma + 5) +
                9 * gamma**3 - 14 * gamma**2 + gamma + 2
            ) +
            N * (
                alpha**4 * (-(1/2) * gamma**2 + gamma - 1) +
                alpha**3 * (-5 * gamma**3 - 2 * gamma**2 + 6 * gamma - (11/2)) +
                alpha**2 * (4 * gamma**4 - (41/2) * gamma**3 - (5/2) * gamma**2 + 13 * gamma - (21/2)) +
                alpha * (6 * gamma**4 - 16 * gamma**3 - (5/2) * gamma**2 + 12 * gamma - (17/2)) +
                (7/2) * gamma**4 - (13/2) * gamma**3 - (3/2) * gamma**2 + 4 * gamma - (5/2)
            ) +
            alpha**4 * (2 * gamma**2 + 1) +
            alpha**3 * (8 * gamma**2 - (3/2) * gamma + 4) +
            alpha**2 * (-(5/2) * gamma**4 + (25/2) * gamma**2 - (9/2) * gamma + 6) +
            alpha * ((1/2) * gamma**5 - 2 * gamma**4 - (1/2) * gamma**3 + 9 * gamma**2 - (9/2) * gamma + 4) +
            (1/2) * gamma**5 - gamma**4 - (1/2) * gamma**3 + (5/2) * gamma**2 - (3/2) * gamma + 1
        )
    )
    return numerator / denominator
def ceiling_fan(N,epsilon,beta):
    numerator = (2 + epsilon) * (-1 + N + beta) * (
    78 + 4 * N**2 * (9 + 11 * epsilon + 2 * epsilon**2) +
    epsilon**2 * (17 - 22 * beta) - 122 * beta + 24 * beta**2 +
    epsilon * (83 - 125 * beta + 22 * beta**2) + 
    N * (-114 + 68 * beta + 5 * epsilon**2 * (-5 + 2 * beta) + epsilon * (-127 + 74 * beta))
)
    denominator = (
    -57 + N**3 * (9 + 25 * epsilon + 34 * epsilon**2 + 6 * epsilon**3) +
    157 * beta - 155 * beta**2 + 15 * beta**3 +
    epsilon**3 * (-12 + 29 * beta - 22 * beta**2) +
    epsilon**2 * (-76 + 201 * beta - 177 * beta**2 + 22 * beta**3) +
    epsilon * (-109 + 294 * beta - 273 * beta**2 + 28 * beta**3) +
    N**2 * (25 * (-3 + beta) + 48 * epsilon**2 * (-3 + 2 * beta) + 2 * epsilon**3 * (-12 + 7 * beta) +
            3 * epsilon * (-53 + 26 * beta)) +
    N * (123 - 182 * beta + 47 * beta**2 + epsilon**3 * (30 - 43 * beta + 10 * beta**2) +
         3 * epsilon**2 * (62 - 99 * beta + 32 * beta**2) + 3 * epsilon * (81 - 124 * beta + 35 * beta**2))
)
    return numerator / denominator
'''
network_type = 'ca-sandi_auths'
filename = '.\\saves\\nets\\' + network_type + '.txt'
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
eta=Find_CT(G)
theorem=Get_BCRatio(G)
G=Trans_to_Selfloop(G, fpfunc)
print(Get_BCRatio(G))'''