# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 13:24:22 2023

@author: zziya
"""
import numpy as np
import random
import networkx as nx
import json
import GetBCRatio
fpfunc='e^-k'#(failure probability for replacement)
realset=['dolphins','americanfootball','contiguous-usa','zebra','ca-sandi_auths','rt-retweet']
fpfuncs=['(k+1)^-1','e^-k','1-k^-1','lnk']
network_type='contiguous-usa'
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
def normalized_laplacian(adj_matrix):
    n = len(adj_matrix)
    P = np.zeros((n, n))
    for row in range(n):
        sumOverRow = np.sum(adj_matrix[row])
        for col in range(n):
            P[row, col] = adj_matrix[row][col] / sumOverRow
    return P
def get_fitness(G, node_state, node, b, nl):
    c=1
    result=-c*(1-node_state[node])
    delta=0.01
    neis=G.neighbors(node)
    for nei in neis:
        result+=b*nl[node][nei]*(1-node_state[nei])
    return np.exp(delta*result)
def DB_update(changer, G, node_state,b, nl):
    kcfc=0
    kdfd=0
    neis=G.neighbors(changer)
    for nei in neis:
        if node_state[nei]==0:
            kcfc+=G[changer][nei]['weight']*get_fitness(G, node_state, nei, b, nl)
        if node_state[nei]==1:
            kdfd+=G[changer][nei]['weight']*get_fitness(G, node_state, nei, b, nl)
    if random.uniform(0,1)<(kcfc/(kcfc+kdfd)):
        return 0
    else:
        return 1
def fc(G, node_state):
    return 1-sum(node_state.values())/G.number_of_nodes()
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
    if fpfunc=='lnk':
        return np.log(k)
    if fpfunc=='e^k':
        return np.exp(k)
def Trans_to_Selfloop(g):
    for i in g.nodes:
        g.add_edge(i,i,weight=Failure_Probability(g.degree(i)))
    return g
N=100
k=4
avgs=500000
B=np.linspace(2,6,20)
c=1
#Cinvade
#0:C, 1:D
resultc={}
resultd={}
if network_type=='sln':
    L=10
    k=4
    N=L*L
    H=nx.grid_2d_graph(L,L, periodic=True)
    #change the label to accelerate
    new_mapping={}
    for i in H.nodes:
        new_mapping[i]=i[0]*L+i[1]
    G = nx.relabel_nodes(H, new_mapping)
if network_type=='ern':
    G=nx.random_graphs.erdos_renyi_graph(N, 0.2)
if network_type=='nwn':
    G=nx.newman_watts_strogatz_graph(N, k, 0.2)
if network_type=='rrg':
    G=nx.random_regular_graph(k,N)
if network_type=='wsn':
    G=nx.generators.watts_strogatz_graph(N, k, 0.05)
if network_type=='ban':
    G=nx.generators.barabasi_albert_graph(N, int(k/2))
failure_type='df'
'''
filename='.\\datas\\'+network_type+'.txt'
G=nx.Graph()
with open(filename) as file:
    for line in file:
        head, tail=[str(x) for x in line.split()]
        G.add_edge(int(head),int(tail))'''
if network_type in realset:
    k=0
    filename='.\\saves\\nets\\'+network_type+'.txt'
    G=nx.Graph()
    with open(filename) as file:
        for line in file:
            head, tail=[str(x) for x in line.split()]
            G.add_edge(int(head),int(tail),weight=1)
# relabel
if 0 not in G.nodes():
    j = 0
    new_mapping = {}
    for i in G.nodes:
        new_mapping[i] = j
        j += 1
    G = nx.relabel_nodes(G, new_mapping)
G=Trans_to_Selfloop(G)
nl=normalized_laplacian(adjacency(G))
THEOREM=GetBCRatio.Get_BCRatio(G)
if THEOREM-2>1:
    B=np.linspace(THEOREM-1,THEOREM+1,7)
else:
    B=np.linspace(1, THEOREM+1,7)
N=G.number_of_nodes()
#initialization

for b in B:
    c_absorb_number=0
    print('\r DF, C, b : {:.8f}, net : {}, fp : {}'.format(b, network_type,fpfunc), end='    ')
    for avg in range(avgs):
        node_state={}
        for i in G.nodes:
            node_state[i]=1
        node_state[random.choice(list(G.nodes))]=0
        t=0
        while (True):
            t+=1
            updater=random.choice(list(G.nodes))
            node_state[updater]=DB_update(updater, G, node_state,b,nl)
            fc_=fc(G, node_state)
            if (fc_==1):
                c_absorb_number+=1
                break
            elif (fc_==0):
                break
    resultc['{:.5f}'.format(b)]=c_absorb_number/avgs
json_str=json.dumps(resultc)
with open('.\\saves\\{}\\C_{}_{:d}_{:d}_{}.json'.format(failure_type,network_type,k,N,fpfunc), 'w') as json_file:
    json_file.write(json_str)
for b in B:
    d_absorb_number=0
    print('\r DF, D, b : {:.8f}, net : {}, fp : {}'.format(b, network_type,fpfunc), end='    ')
    for avg in range(avgs):
        node_state={}
        for i in G.nodes:
            node_state[i]=0
        node_state[random.choice(list(G.nodes))]=1
        t=0
        while (True):
            t+=1
            updater=random.choice(list(G.nodes))
            node_state[updater]=DB_update(updater, G, node_state,b,nl)
            fc_=fc(G, node_state)
            if (fc_==0):
                d_absorb_number+=1
                break
            elif (fc_==1):
                break
    resultd['{:.5f}'.format(b)]=d_absorb_number/avgs
json_str=json.dumps(resultd)
with open('.\\saves\\{}\\D_{}_{:d}_{:d}_{}.json'.format(failure_type,network_type,k,N,fpfunc), 'w') as json_file:
    json_file.write(json_str)