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


def normalized_laplacian(adj_matrix):
    n = len(adj_matrix)
    P = np.zeros((n, n))
    for row in range(n):
        sumOverRow = np.sum(adj_matrix[row])
        for col in range(n):
            P[row, col] = adj_matrix[row][col] / sumOverRow
    return P


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
    return x,Tau


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
    X,tau = Find_CT(G)
    t1 = np.sum(pi.T @ (p * tau))
    t2 = np.sum(pi.T @ (p @ p * tau))
    t3 = np.sum(pi.T @ (p @ p @ p * tau))
    print(t1, t2, t3)
    return (t2) / (t3 - t1)
'''
network_type = 'dolphins'
filename = '.\\saves\\nets\\' + network_type + '.txt'

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
N = G.number_of_nodes()
print(Get_BCRatio(G))
'''
G=nx.generators.random_regular_graph(4,20)
adj = adjacency(G)
pi = random_walk_stationary(G)
p = normalized_laplacian(adj)
X,tau = Find_CT(G)
print(Get_BCRatio(G))