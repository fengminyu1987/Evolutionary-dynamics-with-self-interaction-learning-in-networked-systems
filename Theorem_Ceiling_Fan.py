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
def fp(x,epsilon,beta):
    if x==1:
        return epsilon
    if x!=1:
        return beta
def Trans_to_Selfloop(g,epsilon,beta):
    for i in g.nodes:
        g.add_edge(i,i,weight=fp(g.degree(i),epsilon,beta))
    return g
def Generate_Celing_Fan(N):
    G=nx.star_graph(N-1)
    for i in range(1,N):
        if i>N-i:
            break
        G.add_edge(i,N-i)
    return G
'''
N=10
epsilon=1
beta=3
G=Generate_Star_HH(N)
for u, v in G.edges():
    G.edges[u, v]['weight']=1
G=Trans_to_Selfloop(G,epsilon,beta)
print(GetBCRatio.star_HH(N,epsilon,beta))
print(GetBCRatio.Get_BCRatio(G))
'''
#N>6
plt.figure(figsize=(8,8))
Ns=[7,9,11,13]
betas=np.linspace(0,6,50)
epsilon=0
plt.figure(figsize=(8,8))
for N in Ns:
    temp_result=[]
    for beta in betas:
        temp_result.append(GetBCRatio.ceiling_fan(N,epsilon, beta))
    plt.scatter(betas, temp_result, s=180,label='$N={:d}$'.format(N),marker=markers[Ns.index(N)],color=colors[Ns.index(N)])
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('$\\beta$',fontsize=25)
plt.ylabel('$(\\frac{b}{c})^*_{\\mathbf{\\mathcal{L}}}$',fontsize=35,rotation=0,labelpad=30)
plt.legend(fontsize=25)
#plt.text(2,-30,'Spite',fontsize=35)
#plt.ylim([-60,-1])
plt.grid(ls='--')
plt.savefig('.\\imgs\\6(a).pdf',dpi=400,bbox_inches='tight')
plt.show()
#beta=0
plt.figure(figsize=(8,8))
Ns=[7,9,11,13]
epsilons=np.linspace(0,2,50)
beta=0
plt.figure(figsize=(8,8))
for N in Ns:
    temp_result=[]
    for epsilon in epsilons:
        temp_result.append(GetBCRatio.ceiling_fan(N,epsilon, beta))
    plt.scatter(epsilons, temp_result, s=180,label='$N={:d}$'.format(N),marker=markers[Ns.index(N)],color=colors[Ns.index(N)])
plt.yticks(fontsize=25)
plt.xlabel('$\\epsilon$',fontsize=25)
plt.xticks(np.arange(0,2.1,step=0.5),fontsize=25)
#plt.text(1.5,10,'Cooperation',fontsize=35)
#plt.ylabel('$(\\frac{b}{c})^*_{\\mathbf{\\mathcal{L}}}$',fontsize=35,rotation=0,labelpad=20)
#plt.legend(fontsize=25)
#plt.ylim([1,15])
plt.grid(ls='--')
plt.savefig('.\\imgs\\6(b).pdf',dpi=400,bbox_inches='tight')
plt.show()
colors=['#CE9FCA','#F38F3C','#415C42','#FAE93B']
#N=3 and N=5, evolution of spite
plt.figure(figsize=(8,8))
N=3
epsilons=np.linspace(0,2,50)
betas=[0,0.7,1.1,1.5]
plt.figure(figsize=(8,8))
for beta in betas:
    temp_result=[]
    for epsilon in epsilons:
        temp_result.append(GetBCRatio.ceiling_fan(N,epsilon, beta))
    plt.scatter(epsilons, temp_result, s=180,label='$\\beta={:.2f}$'.format(beta),marker=markers[betas.index(beta)],color=colors[betas.index(beta)])
plt.yticks(fontsize=25)
plt.xlabel('$\\epsilon$',fontsize=25)
plt.xticks(np.arange(0,2.1,step=0.5),fontsize=25)
plt.text(0,50,'C',fontsize=25)
plt.text(0,-200,'S',fontsize=25)
plt.arrow(1.53, -50, 0, 40,width=0.05, head_width=0.13, head_length=10, fc='grey', ec='grey',alpha=0.7)
#plt.text(1.5,10,'Cooperation',fontsize=35)
#plt.ylabel('$(\\frac{b}{c})^*_{\\mathbf{\\mathcal{L}}}$',fontsize=35,rotation=0,labelpad=20)
plt.legend(fontsize=25)
plt.axhline(0,c='purple',lw=3)
plt.ylim([-250,100])
plt.ylabel('$(\\frac{b}{c})^*_{\\mathbf{\\mathcal{L}}}$',fontsize=35,rotation=0,labelpad=30)
plt.grid(ls='--')
plt.savefig('.\\imgs\\6(c).pdf',dpi=400,bbox_inches='tight')
plt.show()
plt.figure(figsize=(8,8))
N=5
epsilons=np.linspace(0,1,50)
betas=[0,0.7,1.1,1.5]
plt.figure(figsize=(8,8))
for beta in betas:
    temp_result=[]
    for epsilon in epsilons:
        temp_result.append(GetBCRatio.ceiling_fan(N,epsilon, beta))
    plt.scatter(epsilons, temp_result, s=180,label='$\\beta={:.2f}$'.format(beta),marker=markers[betas.index(beta)],color=colors[betas.index(beta)])
plt.yticks(fontsize=25)
plt.axhline(0,c='purple',lw=3)
plt.xlabel('$\\epsilon$',fontsize=25)
plt.xticks(np.arange(0,1.1,step=0.25),fontsize=25)
plt.text(0.5,50,'C',fontsize=25)
plt.text(0.5,-50,'S',fontsize=25)
plt.arrow(0.277, -50, 0, 40,width=0.025, head_width=0.08, head_length=10, fc='grey', ec='grey',alpha=0.7)
#plt.text(1.5,10,'Cooperation',fontsize=35)
#plt.ylabel('$(\\frac{b}{c})^*_{\\mathbf{\\mathcal{L}}}$',fontsize=35,rotation=0,labelpad=20)
plt.legend(fontsize=25)
plt.ylim([-250,100])
plt.grid(ls='--')
plt.savefig('.\\imgs\\6(d).pdf',dpi=400,bbox_inches='tight')
plt.show()