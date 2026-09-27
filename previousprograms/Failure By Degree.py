# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 13:24:22 2023

@author: zziya
"""
import numpy as np
import random
import networkx as nx
import json
fpfunc='k^-1'#(failure probability for replacement)
realset=['dolphins','americanfootball','contiguous-usa','zebra','ca-sandi_auths','rt-retweet']
fpfuncs=['k^-1','e^-k','k^-2','1-k^-1','1-e^-k','1-k^-2']
def get_fitness(G, node_state, node, payoffmat):
    result=0
    delta=0.01
    neis=G.neighbors(node)
    for nei in neis:
        result+=payoffmat[node_state[node]][node_state[nei]]
    return np.exp(delta*result)
def DB_update(changer, G, node_state,payoffmat):
    kcfc=0
    kdfd=0
    neis=G.neighbors(changer)
    for nei in neis:
        if node_state[nei]==0:
            kcfc+=get_fitness(G, node_state, nei, payoffmat)
        if node_state[nei]==1:
            kdfd+=get_fitness(G, node_state, nei, payoffmat)
    if random.uniform(0,1)<(kcfc/(kcfc+kdfd)):
        return 0
    else:
        return 1
def fc(G, node_state):
    return 1-sum(node_state.values())/G.number_of_nodes()
def Failure_Probability(k):
    if fpfunc=='k^-1':
        return 1/k
    if fpfunc=='e^-k':
        return np.exp(-k)
    if fpfunc=='k^-2':
        return 1/(k*k)
    if fpfunc=='1-k^-1':
        return 1-1/k
    if fpfunc=='1-e^-k':
        return 1-np.exp(-k)
    if fpfunc=='1-k^-2':
        return 1-1/(k*k)
N=100
k=4
avgs=50000
B=np.linspace(2,6,20)
c=1
#Cinvade
#0:C, 1:D
resultc={}
resultd={}
network_type='dolphins'
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
failure_type='uf'
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
            G.add_edge(int(head),int(tail))
    if network_type=='dolphins':
        B=np.linspace(4.5,8.5,20)
    if network_type=='americanfootball':
        B=np.linspace(10.8,14.8,20)
    if network_type=='contiguous-usa':
        B=np.linspace(3.5,7.5,20)
    if network_type=='zebra':
        B=np.linspace(6,10,20)
    if network_type=='ca-sandi_auths':
        B=np.linspace(1.24,5.24,20)
    if network_type=='rt-retweet':
        B=np.linspace(1,4.65,20)
N=G.number_of_nodes()
#initialization
for b in B:
    c_absorb_number=0
    payoffmat=[[b-c, -c], [b, 0]]
    print('\r DF, C, b : {:.2f}, net : {}, fp : {:.2f}'.format(b, network_type,fpfunc), end='    ')
    for avg in range(avgs):
        node_state={}
        for i in G.nodes:
            node_state[i]=1
        node_state[random.choice(list(G.nodes))]=0
        t=0
        while (True):
            t+=1
            updater=random.choice(list(G.nodes))
            fp=Failure_Probability(G.degree(updater))
            if random.uniform(0,1)<1-fp:
                node_state[updater]=DB_update(updater, G, node_state,payoffmat)
            fc_=fc(G, node_state)
            if (fc_==1):
                c_absorb_number+=1
                break
            elif (fc_==0):
                break
    resultc['{:.2f}'.format(b)]=c_absorb_number/avgs
json_str=json.dumps(resultc)
with open('.\\saves\\{}\\C_{}_{:d}_{:d}_{:.2f}.json'.format(failure_type,network_type,k,N,fpfunc), 'w') as json_file:
    json_file.write(json_str)
for b in B:
    d_absorb_number=0
    payoffmat=[[b-c, -c], [b, 0]]
    print('\r UF, D, b : {:.2f}, net : {}, fp : {:.2f}'.format(b, network_type,fp), end='    ')
    for avg in range(avgs):
        node_state={}
        for i in G.nodes:
            node_state[i]=0
        node_state[random.choice(list(G.nodes))]=1
        t=0
        while (True):
            t+=1
            updater=random.choice(list(G.nodes))
            fp=Failure_Probability(G.degree(updater))
            if random.uniform(0,1)<1-fp:
                node_state[updater]=DB_update(updater, G, node_state,payoffmat)
            fc_=fc(G, node_state)
            if (fc_==0):
                d_absorb_number+=1
                break
            elif (fc_==1):
                break
    resultd['{:.2f}'.format(b)]=d_absorb_number/avgs
json_str=json.dumps(resultd)
with open('.\\saves\\{}\\D_{}_{:d}_{:d}_{:.2f}.json'.format(failure_type,network_type,k,N,fp), 'w') as json_file:
    json_file.write(json_str)