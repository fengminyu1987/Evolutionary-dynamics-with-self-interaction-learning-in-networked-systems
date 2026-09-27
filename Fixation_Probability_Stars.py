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
p=1
q=2
network_type='joined_star'
N=6
avgs=1000000
c=1
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
def fp(k,p,q):
    if k==1 or k==2:
        return p
    if k==N-1 or k==N:
        return q
def Trans_to_Selfloop(g):
    for i in g.nodes:
        g.add_edge(i,i,weight=fp(G.degree(i),p,q))
    return g
def Generate_Star_HH(N):
    g1=nx.star_graph(N-1)
    g2=nx.star_graph(N-1)
    G=nx.disjoint_union(g1,g2)
    G.add_edge(0,N)
    return G
#Cinvade
#0:C, 1:D
resultc={}
resultd={}
failure_type='df'
if network_type=='star':
    G=nx.star_graph(N-1)
if network_type=='joined_star':
    G=Generate_Star_HH(N)
for u, v in G.edges():
    G.edges[u, v]['weight']=1
G=Trans_to_Selfloop(G)
THEOREM=GetBCRatio.star_HH(N,fp(1,p,q),fp(N,p,q))
nl=normalized_laplacian(adjacency(G))
print(GetBCRatio.Get_BCRatio(G),THEOREM)
if THEOREM>2:
    B=np.linspace(THEOREM-1,THEOREM+1,7)
else:
    B=np.linspace(1, THEOREM+1,7)
N=G.number_of_nodes()
#initialization
for b in B:
    c_absorb_number=0
    print('\r C, b : {:.8f}, net : {}, p: {:.2f}, q={:.2f}'.format(b, network_type,p,q), end='    ')
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
with open('.\\saves\\{}\\C_{}_{}_{:.2f}_{:.2f}.json'.format(failure_type,network_type,N,p,q), 'w') as json_file:
    json_file.write(json_str)
'''
for b in B:
    d_absorb_number=0
    print('\r D, b : {:.8f}, net : {}, p: {:.2f}, q={:.2f}'.format(b, network_type,p,q), end='    ')
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
with open('.\\saves\\{}\\D_{}_{}_{:.2f}_{:.2f}.json'.format(failure_type,network_type,N,p,q), 'w') as json_file:
    json_file.write(json_str)'''