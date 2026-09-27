# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 23:39:58 2024

parameter pairs: N=50; 100
k=4; 30

@author: zziya
"""

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import GetBCRatio
markers=['^','>','<','v']
colors=['#A51C36','#7ABBDB','#84BA42','#DBB428']
def fp(x,alpha,gamma):
    if x==1:
        return alpha
    if x!=1:
        return gamma
def Trans_to_Selfloop(g,alpha,gamma):
    for i in g.nodes:
        g.add_edge(i,i,weight=fp(g.degree(i),alpha,gamma))
    return g
def Generate_Star_HH(N):
    g1=nx.star_graph(N-1)
    g2=nx.star_graph(N-1)
    G=nx.disjoint_union(g1,g2)
    G.add_edge(0,N)
    return G

N=10
alpha=0.2
gamma=0.8
'''
G=Generate_Star_HH(N)
for u, v in G.edges():
    G.edges[u, v]['weight']=1
G=Trans_to_Selfloop(G,alpha,gamma)
print(GetBCRatio.star_HH(N,alpha,gamma))
print(GetBCRatio.Get_BCRatio(G))
'''
#alpha=0
plt.figure(figsize=(8,8))
Ns=[5,7,9,11]
gammas=np.linspace(0,6,50)
alpha=0
plt.figure(figsize=(8,8))
for N in Ns:
    temp_result=[]
    for gamma in gammas:
        temp_result.append(GetBCRatio.star_HH(N,alpha, gamma))
    plt.scatter(gammas, temp_result, s=180,label='$N={:d}$'.format(N),marker=markers[Ns.index(N)],color=colors[Ns.index(N)])
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('$\\gamma$',fontsize=25)
plt.ylabel('$(\\frac{b}{c})^*_{\\mathbf{\\mathcal{L}}}$',fontsize=35,rotation=0,labelpad=30)
plt.legend(fontsize=25)
#plt.text(2,-30,'Spite',fontsize=35)
#plt.ylim([-60,-1])
plt.grid(ls='--')
plt.savefig('.\\imgs\\4(a).pdf',dpi=400,bbox_inches='tight')
plt.show()
#gamma=0
plt.figure(figsize=(8,8))
Ns=[5,7,9,11]
Alphas=np.linspace(0,3,50)
gamma=0
plt.figure(figsize=(8,8))
for N in Ns:
    temp_result=[]
    for alpha in Alphas:
        temp_result.append(GetBCRatio.star_HH(N,alpha, gamma))
    plt.scatter(Alphas, temp_result, s=180,label='$N={:d}$'.format(N),marker=markers[Ns.index(N)],color=colors[Ns.index(N)])
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('$\\alpha$',fontsize=25)
#plt.text(1.5,10,'Cooperation',fontsize=35)
#plt.ylabel('$(\\frac{b}{c})^*_{\\mathbf{\\mathcal{L}}}$',fontsize=35,rotation=0,labelpad=20)
#plt.legend(fontsize=25)
#plt.ylim([1,15])
plt.grid(ls='--')
plt.savefig('.\\imgs\\4(b).pdf',dpi=400,bbox_inches='tight')
plt.show()