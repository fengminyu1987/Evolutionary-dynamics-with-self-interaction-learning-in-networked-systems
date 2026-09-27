# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 21:35:40 2024
Transition Diagrams of (b/c)_L^* and (b/c)^*, k, knn, density
Network Numbers for One Figure: 1.32 \times 10^5
@author: zziya
"""
import networkx as nx
import GetBCRatio
import pandas as pd
import csv
import numpy as np
import time
#Parameters
NetSet=['ban','wsn','nwn','ern','hkn','ddn']
network_type='nwn'
Ns=np.arange(40,141,1)
fpfunc='lnk'
fpfuncs=['(k+1)^-1','e^-k','1-k^-1','lnk']
def Failure_Probability(k):
    if fpfunc=='(k+1)^-1':
        return 1/(k+1)
    if fpfunc=='1-k^-1':
        return 1-1/k
    if fpfunc=='lnk':
        return np.log(k)
    if fpfunc=='e^-k':
        return np.exp(-k)
def Trans_to_Selfloop(g):
    for i in g.nodes:
        g.add_edge(i,i,weight=Failure_Probability(g.degree(i)))
    return g
def write_csv(row):
    with open('.\\saves\\comparison\\{}_{}.csv'.format(network_type,fpfunc),'a',newline='') as result:
        csv_write=csv.writer(result)
        csv_write.writerow(row)
#Preparation for Network Parameters
if network_type=='ban':#2e+4
    ms=np.arange(1,2,1)#10
    avgs=1
if network_type=='wsn':#2.4e+4
    ks=np.arange(2,9,2)#4
    ps=np.arange(0.1,0.251,0.03)#6
    avgs=10
if network_type=='ern':#1.2e+4
    ps=np.arange(0.1, 0.31, 0.025)#8
    avgs=15
if network_type=='nwn':#2.4e+4
    ks=np.arange(2,9,2)#4
    ps=np.arange(0.1,0.251,0.03)#6
    avgs=10
if network_type=='hkn':#2.4e+4
    ms=np.arange(1,8,2)#4
    ps=np.arange(0.1,0.251,0.03)#6
    avgs=10
if network_type=='rpt':#2e+4
    gammas=np.arange(2.5,3.5,0.1)#10
    avgs=20
if network_type=='ddn':#2.2e+4
    ps=np.arange(0.1,0.32,0.02)#11
    avgs=20
#Repair Each Line
Modified_Number=0
with open('.\\saves\\comparison\\{}_{}.csv'.format(network_type,fpfunc), mode='r', newline='') as infile, open('.\\saves\\comparison\\{}_{}_repair.csv'.format(network_type,fpfunc), mode='w', newline='') as outfile:
    reader = csv.DictReader(infile)
    fieldnames = reader.fieldnames
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in reader:
        if float(row['bcl'])<1 or float(row['bcl'])>200 or float(row['bc'])<1 or float(row['bc'])>200 or row['bcl']=='nan' or row['bc']=='nan':#Then Modify
            if network_type=='ban':
                temp=row['net']
                temp=temp.strip('()')
                temp=temp.split(", ",2)
                N=int(temp[0])
                m=int(temp[1])
                netparam=(N,m)
                G=nx.barabasi_albert_graph(N,m)
                while not nx.is_connected(G):
                    G=nx.barabasi_albert_graph(N,m)
            if network_type=='wsn':
                temp=row['net']
                temp=temp.strip("()")
                temp=temp.split(", ",2)
                N=int(temp[0])
                k=int(temp[1])
                p=float(temp[2])
                G=nx.watts_strogatz_graph(N, k, p)
                while not nx.is_connected(G):
                    G=nx.watts_strogatz_graph(N, k, p)
                netparam=(N,k,p)
            if network_type=='nwn':
                temp=row['net']
                temp=temp.strip("()")
                temp=temp.split(", ",2)
                N=int(temp[0])
                k=int(temp[1])
                p=float(temp[2])
                G=nx.newman_watts_strogatz_graph(N, k, p)
                while not nx.is_connected(G):
                    G=nx.newman_watts_strogatz_graph(N, k, p)
                netparam=(N,k,p)
            if network_type=='ern':
                temp=row['net']
                temp=temp.strip("()")
                temp=temp.split(", ",2)
                N=int(temp[0])
                p=float(temp[1])
                G=nx.erdos_renyi_graph(N, p)
                while not nx.is_connected(G):
                    G=nx.erdos_renyi_graph(N, p)
                netparam=(N,p)
            if network_type=='hkn':
                temp=row['net']
                temp=temp.strip("()")
                temp=temp.split(", ",2)
                N=int(temp[0])
                m=int(temp[1])
                p=float(temp[2])
                G=nx.powerlaw_cluster_graph(N, m, p)
                while not nx.is_connected(G):
                    G=nx.powerlaw_cluster_graph(N, m, p)
                netparam=(N,m,p)
            if network_type=='ddn':
                temp=row['net']
                temp=temp.strip("()")
                temp=temp.split(", ",2)
                N=int(temp[0])
                p=float(temp[1])
                G=nx.duplication_divergence_graph(N, p)
                while not nx.is_connected(G):
                    G=nx.duplication_divergence_graph(N, p)
                netparam=(N,p)
            #print("\r NET : {}, fpfunc : {}, PARAM : {}, Modified: {:d}".format(network_type, fpfunc, netparam,Modified_Number),end='        ')
            for u, v in G.edges():
                G.edges[u, v]['weight']=1
            row['bc']=GetBCRatio.Get_BCRatio(G)
            row['k']=2*G.number_of_edges()/G.number_of_nodes()
            knn=nx.average_neighbor_degree(G)
            row['knn']=np.mean(list(knn.values()))
            row['density']=nx.density(G)
            G=Trans_to_Selfloop(G)
            row['bcl']=GetBCRatio.Get_BCRatio(G)
            Modified_Number+=1
            print("\r NET : {}, fpfunc : {}, PARAM : {}, bc : {}, bcl : {}, Modified: {:d}".format(network_type, fpfunc, netparam,row['bc'],row['bcl'],Modified_Number),end='        ')
        writer.writerow(row)