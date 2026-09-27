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
def fp(x,alpha,beta):
    if x==1:
        return alpha
    if x!=1:
        return beta
def Trans_to_Selfloop(g,alpha,beta):
    for i in g.nodes:
        g.add_edge(i,i,weight=fp(g.degree(i),alpha,beta))
    return g
def star(N,alpha, beta):
    # 分子部分
    numerator = 3 * (1 + alpha) * (-1 + N + beta) * (
        -(8/3) 
        + N**3 * (4/3 + (11 * alpha)/3 + alpha**2) 
        + (14 * beta)/3 
        - (4 * beta**2)/3 
        + alpha**2 * (-2 + (8 * beta)/3) 
        + alpha * (-6 + 10 * beta - (8 * beta**2)/3) 
        + N**2 * (-(16/3) + 3 * beta + alpha**2 * (-4 + (4 * beta)/3) + alpha * (-(40/3) + 7 * beta)) 
        + N * (20/3 + alpha**2 * (5 - 4 * beta) - (23 * beta)/3 + (4 * beta**2)/3 + alpha * (47/3 - 17 * beta + (8 * beta**2)/3))
    )

    # 分母部分
    denominator = 2 * (
        beta * (-N**3 + N * (-1 - beta/2) + N**2 * (2 - beta/2) + beta) 
        + alpha**3 * (2 + N**4 - 5 * beta + 4 * beta**2 + N**3 * (-5 + (5 * beta)/2) + N * (-7 + (25 * beta)/2 - 6 * beta**2) + N**2 * (9 - 10 * beta + 2 * beta**2)) 
        + alpha * (4 + N**4 - 11 * beta + 12 * beta**2 - 2 * beta**3 + N**3 * (-7 + (5 * beta)/2) + N**2 * (15 - 16 * beta + 5 * beta**2) + N * (-13 + (49 * beta)/2 - 17 * beta**2 + 2 * beta**3)) 
        + alpha**2 * (8 + 4 * N**4 - 22 * beta + 21 * beta**2 - 4 * beta**3 + N**3 * (-20 + 12 * beta) + N**2 * (36 - 46 * beta + (27 * beta**2)/2) + N * (-28 + 56 * beta - (69 * beta**2)/2 + 4 * beta**3))
    )
    # 返回分子/分母的值
    return numerator / denominator
'''N=10
alpha=1
beta=2
G=nx.star_graph(N-1)
for u, v in G.edges():
    G.edges[u, v]['weight']=1
G=Trans_to_Selfloop(G,alpha,beta)
print(star(N,alpha,beta))
print(GetBCRatio.Get_BCRatio(G))'''
#Spite, alpha=0
plt.figure(figsize=(8,8))
Ns=[5,7,9,11]
Betas=np.linspace(0,3,50)
alpha=0
plt.figure(figsize=(8,8))
for N in Ns:
    temp_result=[]
    for beta in Betas:
        temp_result.append(star(N,alpha, beta))
    plt.scatter(Betas, temp_result, s=180,label='$N={:d}$'.format(N),marker=markers[Ns.index(N)],color=colors[Ns.index(N)])
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('$\\beta$',fontsize=25)
plt.ylabel('$(\\frac{b}{c})^*_{\\mathbf{\\mathcal{L}}}$',fontsize=35,rotation=0,labelpad=30)
plt.legend(fontsize=25)
plt.text(2,-30,'S',fontsize=35)
plt.ylim([-60,-1])
plt.grid(ls='--')
plt.savefig('.\\imgs\\3(a).pdf',dpi=400,bbox_inches='tight')
plt.show()
#Cooperation, beta=0
plt.figure(figsize=(8,8))
Ns=[5,7,9,11]
Alphas=np.linspace(0,3,50)
beta=0
plt.figure(figsize=(8,8))
for N in Ns:
    temp_result=[]
    for alpha in Alphas:
        temp_result.append(star(N,alpha, beta))
    plt.scatter(Alphas, temp_result, s=180,label='$N={:d}$'.format(N),marker=markers[Ns.index(N)],color=colors[Ns.index(N)])
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('$\\alpha$',fontsize=25)
plt.text(1.5,10,'C',fontsize=35)
#plt.ylabel('$(\\frac{b}{c})^*_{\\mathbf{\\mathcal{L}}}$',fontsize=35,rotation=0,labelpad=20)
#plt.legend(fontsize=25)
plt.ylim([1,15])
plt.grid(ls='--')
plt.savefig('.\\imgs\\3(b).pdf',dpi=400,bbox_inches='tight')
plt.show()