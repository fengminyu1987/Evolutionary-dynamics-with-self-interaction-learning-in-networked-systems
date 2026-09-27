# calculate coalescence times for any connected graph
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import time
fpfuncs=['(k+1)^-1','e^-k','(k+1)^-2','1-k^-1']
fpfunc='1-k^-1'
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
    if fpfunc=='(k+1)^-3':
        return 1/(k+1)**3
def adjacency(G):
    A=[]
    for i in G.nodes:
        Atemp=[]
        for j in G.nodes:
            if G.has_edge(i,j):
                Atemp.append(1)
            else:
                Atemp.append(0)
        A.append(Atemp)
    return A
# some examples for graphs of size five
# df = pd.DataFrame({'from': ['A', 'A', 'A', 'A'], 'to': ['B', 'C', 'D', 'E']})
# df = pd.DataFrame({'from': ['A', 'A', 'B', 'B', 'C', 'C'], 'to': ['B', 'C', 'D', 'E', 'D', 'E']})
# df = pd.DataFrame({'from': ['A', 'A', 'B', 'B', 'C', 'D'], 'to': ['B', 'C', 'C', 'D', 'E', 'E']})
# df = pd.DataFrame({'from': ['A', 'B', 'B', 'C', 'D'], 'to': ['B', 'C', 'D', 'D', 'E']})
# df = pd.DataFrame({'from': ['A', 'B', 'B', 'C', 'E'], 'to': ['B', 'C', 'D', 'E', 'A']})
# df = pd.DataFrame({'from': ['A', 'B', 'C', 'C', 'D'], 'to': ['B', 'C', 'D', 'E', 'E']})
# df = pd.DataFrame({'from': ['A', 'B', 'C', 'D', 'E'], 'to': ['B', 'C', 'D', 'E', 'A']})
# df = pd.DataFrame({'from': ['A', 'B', 'C', 'D'], 'to': ['B', 'C', 'D', 'E']})
#df = pd.DataFrame({'from': ['A', 'B', 'C', 'C'], 'to': ['B', 'C', 'D', 'E']})
#G = nx.from_pandas_edgelist(df, 'from', 'to')
#nx.draw(G, with_labels=True)
#plt.show()
st=time.time()
network_type='dolphins'
filename='.\\datas\\'+network_type+'.txt'
G=nx.Graph()
with open(filename) as file:
    for line in file:
        head, tail=[str(x) for x in line.split()]
        G.add_edge(int(head),int(tail))
if 0 not in G.nodes():
    new_mapping={}
    for i in G.nodes:
        new_mapping[i]=i-1
    G = nx.relabel_nodes(G, new_mapping)
# regular graphs
# G = nx.random_regular_graph(4, 100)
adjacencyMatrix = adjacency(G)
n = len(G.nodes)
# One step transition probability matrix
P = np.zeros((n, n))

for row in range(n):
    sumOverRow = np.sum(adjacencyMatrix[row])
    for col in range(n):
        P[row, col] = adjacencyMatrix[row][col] / sumOverRow

A = np.zeros((n * n, n * n))
B = np.ones((n * n, 1))
lists = [i * n + i for i in range(n)]
for row in range(n * n):
    if row in lists:
        A[row, row] = 1
        B[row, 0] = 0
    else:
        i = int(row / n)
        j = row % n
        for k in range(i * n, (i + 1) * n):
            if k % n == j:
                A[row, k] = 1 - 0.5 * P[j, k % n]
            else:
                A[row, k] = -0.5 * P[j, k % n]
        for k in range(j * n, (j + 1) * n):
            A[row, k] = -0.5 * P[i, k % n]

Pi = np.zeros(n)
W = np.sum(adjacencyMatrix)
for i in range(n):
    Pi[i] = np.sum(adjacencyMatrix[i]) / W

X = np.linalg.solve(A, B)
Tau = X.reshape(n, n)
P2 = np.linalg.matrix_power(P, 2)
P3 = np.linalg.matrix_power(P, 3)
tau1 = 0
for i in range(n):
    for j in range(n):
        tau1 += Pi[i] * P[i, j] * Tau[i, j]
tau2 = 0
for i in range(n):
    for j in range(n):
        tau2 += Pi[i] * P2[i, j] * Tau[i, j]
tau3 = 0
for i in range(n):
    for j in range(n):
        tau3 += Pi[i] * P3[i, j] * Tau[i, j]
if tau3 == tau1:
    print(np.inf)
else:
    threshold = tau2 / (tau3 - tau1)
    print(threshold)
en=time.time()
print(en-st)