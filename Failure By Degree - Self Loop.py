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
realset=['dolphins','americanfootball','contiguous-usa','ca-sandi_auths','iceland','rt-retweet','soc-karate','mouse-visual-cortex','chesapeake']
regset=['rrg','lat','wsn']
fpfuncs=['(k+1)^-1','e^-k','1-k^-1','lnk']
network_type='iceland'
k=6
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
def bcrreg(N,k,fpfunc):
    return (N*(k*k+3*k*Failure_Probability(k)+2*Failure_Probability(k)*Failure_Probability(k))-2*(k+Failure_Probability(k))*(k+Failure_Probability(k)))/(N*(2*Failure_Probability(k)*Failure_Probability(k)+k*Failure_Probability(k)+k)-2*(k+Failure_Probability(k))*(k+Failure_Probability(k)))
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
    if fpfunc=='0':
        return 0
def Trans_to_Selfloop(g):
    for i in g.nodes:
        g.add_edge(i,i,weight=Failure_Probability(g.degree(i)))
    return g
N=100
avgs=500000
c=1
#Cinvade
#0:C, 1:D
resultc={}
resultd={}
failure_type='df'
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
if network_type in regset:
    if network_type=='rrg':
        G=nx.generators.random_regular_graph(k,N)
    if network_type=='lat':
        if k==4:
            G=nx.generators.lattice.grid_2d_graph(10,10,periodic=True)
            N=100
        if k==3:
            N=72
            G=nx.generators.lattice.hexagonal_lattice_graph(6,6, periodic=True, with_positions=False)
        if k==6:
            N=98
            G=nx.generators.lattice.triangular_lattice_graph(14,14, periodic=True, with_positions=False)
    if network_type=='wsn':
        G=nx.generators.connected_watts_strogatz_graph(100, k, 0.2, tries=100, seed=None)
    for i,j in G.edges:
        G.add_edge(i, j, weight=1)
    # relabel
    if 0 not in G.nodes():
        j = 0
        new_mapping = {}
        for i in G.nodes:
            new_mapping[i] = j
            j += 1
        G = nx.relabel_nodes(G, new_mapping)
    G=Trans_to_Selfloop(G)
    THEOREM=bcrreg(N,k,fpfunc)
    nl=normalized_laplacian(adjacency(G))
if THEOREM>2:
    B=np.linspace(THEOREM-1,THEOREM+1,7)
else:
    B=np.linspace(1, THEOREM+1,7)
N=G.number_of_nodes()
#initialization
for b in B:
    c_absorb_number=0
    print('\r {:d}, C, b : {:.8f}/{:.8f}, net : {}, fp : {}'.format(k, b,THEOREM, network_type,fpfunc), end='    ')
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
    resultc['{:.8f}'.format(b)]=c_absorb_number/avgs
json_str=json.dumps(resultc)
with open('.\\saves\\{}\\C_{}_{:d}_{:d}_{}.json'.format(failure_type,network_type,k,N,fpfunc), 'w') as json_file:
    json_file.write(json_str)
for b in B:
    d_absorb_number=0
    print('\r {:.2f}, D, b : {:.8f}/{:.8f}, net : {}, fp : {}'.format(k,b,THEOREM, network_type,fpfunc), end='    ')
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
    resultd['{:.8f}'.format(b)]=d_absorb_number/avgs
json_str=json.dumps(resultd)
with open('.\\saves\\{}\\D_{}_{:d}_{:d}_{}.json'.format(failure_type,network_type,k,N,fpfunc), 'w') as json_file:
    json_file.write(json_str)