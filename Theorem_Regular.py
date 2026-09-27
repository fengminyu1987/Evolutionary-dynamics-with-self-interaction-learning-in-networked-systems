# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 23:39:58 2024

parameter pairs: N=50; 100
k=4; 30

@author: zziya
"""

import matplotlib.pyplot as plt
import numpy as np
import math
markers=['^','>','<','v']
colors=['#A51C36','#7ABBDB','#84BA42','#DBB428']
def theorem_regular(N,k,l):
    num=N*(k*k+3*k*l+2*l*l)-2*(k+l)*(k+l)
    den=N*(2*l*l+k*l+k)-2*(k+l)*(k+l)
    return num/den
def transition(N,k):
    return (N*(4-k)+math.sqrt(N*N*k*k+8*N*k+16*(k-N)*(k*N-N-k)))/(4*(N-1))
#Sparse
L=np.linspace(0,6,50)
Pairs=[(50,2),(100,2),(50,3),(100,8)]
plt.figure(figsize=(8,8))
for (N,k) in Pairs:
    temp_result=[]
    for l in L:
        temp_result.append(theorem_regular(N, k, l))
    plt.scatter(L, temp_result, s=180,label='$N={:d}, k={:d}$'.format(N,k),marker=markers[Pairs.index((N,k))],color=colors[Pairs.index((N,k))])
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('$\\ell(k)$',fontsize=25)
plt.ylabel('$(\\frac{b}{c})^*_{\\mathbf{\\mathcal{L}}}$',fontsize=35,rotation=0,labelpad=20)
plt.legend(fontsize=25)
plt.ylim([1,10])
plt.grid(ls='--')
plt.savefig('.\\imgs\\1(a).pdf',dpi=400,bbox_inches='tight')
plt.show()
#Dense
Pairs=[(50,30),(100,60),(50,40),(100,80)]
L=np.linspace(0,6,100)
plt.figure(figsize=(8,8))
sign=40
for (N,k) in Pairs:
    temp_result=[]
    for l in L:
        temp_result.append(theorem_regular(N, k, l))
    plt.scatter(L, temp_result, s=180,label='$N={:d}, k={:d}$'.format(N,k),marker=markers[Pairs.index((N,k))],color=colors[Pairs.index((N,k))])
    sign=-sign
    plt.arrow(transition(N,k), sign, 0, -sign*0.75,width=0.1, head_width=0.3, head_length=10, fc=colors[Pairs.index((N,k))], ec=colors[Pairs.index((N,k))])
plt.axhline(0,c='purple',lw=3)
plt.legend(fontsize=25)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.text(4,150,'C',fontsize=25)
plt.text(1,-100,'S',fontsize=25)
plt.xlabel('$\\ell(k)$',fontsize=25)
#plt.ylabel('$(\\frac{b}{c})^*_{\\mathbf{\\mathcal{L}}}$',fontsize=25,rotation=0,labelpad=30)
plt.legend(fontsize=25)
plt.grid(ls='--')
plt.ylim([-200,200])
plt.savefig('.\\imgs\\1(b).pdf',dpi=400,bbox_inches='tight')
plt.show()