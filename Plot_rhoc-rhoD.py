# -*- coding: utf-8 -*-
"""
Created on Fri May 17 16:32:26 2024

@author: zziya
"""
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import random
import json
import GetBCRatio
#fpfunc='(k+1)^-1'#(failure probability for replacement)
realset=['dolphins','americanfootball','contiguous-usa','rt-retweet','ca-sandi_auths','soc-karate', 'chesapeake','iceland']
regset=['lat','rrg','joined_star']
fpfuncs=['(k+1)^-1','e^-k','lnk','1-k^-1']
#fpfuncs=['lnk']
labs=['$(k+1)^{-1}$','$e^{-k}$','$\ln k$', '$1-k^{-1}$']
ns=[62,115,49,96,86,34,39,75]
network_type='chesapeake'
k=0
n=96
Cs=['#EAB080','#0001A1','#037F77','#C5272D']#Cs=['red', 'blue', 'green', 'black']
Ms=['^', '>', '<', 'v']
if network_type in realset:
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
    arrow=GetBCRatio.Get_BCRatio(G)
plt.figure(figsize=(8,8))
for fpfunc in fpfuncs:
    if network_type in realset:
        filename1='.\\saves\\df\\C_{}_{}.json'.format(network_type,fpfunc)
        filename2='.\\saves\\df\\D_{}_{}.json'.format(network_type,fpfunc)
        n=ns[realset.index(network_type)]
    if network_type in regset:
        filename1='.\\saves\\df\\C_{}_{:d}_{:d}_{}.json'.format(network_type,k,n,fpfunc)
        filename2='.\\saves\\df\\D_{}_{:d}_{:d}_{}.json'.format(network_type,k,n,fpfunc)
    f1=open(filename1)
    f2=open(filename2)
    data1=json.load(f1)
    data2=json.load(f2)
    B=list(data1.keys())
    Bf=[round(float(B[i]),5) for i in range(len(B))]
    THEOREM=(max(Bf)+min(Bf))/2
    rhoc_d=[]
    #rhod=[]
    for b in B:
        rhoc_d.append(data1[b]-data2[b])
    plt.axvline(THEOREM,color=Cs[fpfuncs.index(fpfunc)],linestyle='--')
    plt.scatter(Bf, n*np.array(rhoc_d), color=Cs[fpfuncs.index(fpfunc)],marker=Ms[fpfuncs.index(fpfunc)],s=300, label=labs[fpfuncs.index(fpfunc)])
    z1=np.polyfit(Bf, n*np.array(rhoc_d), 1)
    p1=np.poly1d(z1)
    plt.plot(np.linspace(min(Bf),max(Bf),15), p1(np.linspace(min(Bf),max(Bf),15)),color=Cs[fpfuncs.index(fpfunc)], lw=2)
plt.axhline(0, color = 'purple', alpha=1, linestyle='--',lw=3)
if network_type=='dolphins':
    plt.legend(fontsize=25,loc='best')
    plt.xlim([1,8])
    plt.ylabel('$N\\times(\\rho_C-\\rho_D)$',fontsize=35)
if network_type=='ca-sandi_auths':
    plt.ylabel('$N\\times(\\rho_C-\\rho_D)$',fontsize=35)
if network_type=='lat' and (k==3 or k==4):
    plt.ylabel('$N\\times(\\rho_C-\\rho_D)$',fontsize=35)
    #plt.legend(fontsize=25)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('$b/c$',fontsize=35)
plt.grid()
if network_type in realset:
    plt.arrow(arrow, -0.1, 0, 0.05,width=1.5, head_width=3, head_length=0.05, fc='grey', ec='grey',alpha=0.7)
    plt.savefig(".\\saves\\imgs\\{}.pdf".format(network_type), dpi=500, bbox_inches='tight')
if network_type in regset:
    plt.arrow((n-2)/(n/k-2), -0.2, 0, 0.1,width=0.1, head_width=0.3, head_length=0.1, fc='grey', ec='grey',alpha=0.7)
    plt.ylim([-0.5,0.5])
    plt.savefig(".\\saves\\imgs\\{}_{:d}_{:d}.pdf".format(network_type,n,k), dpi=500, bbox_inches='tight')
plt.show()