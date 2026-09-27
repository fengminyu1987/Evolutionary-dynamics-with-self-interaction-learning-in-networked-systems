# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 22:07:07 2024

@author: zziya
"""

import networkx as nx
import numpy as np
import scipy.sparse as sp
import time
import GetBCRatio
fpfunc='lnk'#(failure probability for replacement)
realset=['dolphins','americanfootball','contiguous-usa','zebra','ca-sandi_auths','rt-retweet']
fpfuncs=['(k+1)^-1','e^-k','1-k^-1','lnk']
def theorem_regular(N,k,l):
    num=N*(k*k+3*k*l+2*l*l)-2*(k+l)*(k+l)
    den=N*(2*l*l+k*l+k)-2*(k+l)*(k+l)
    return num/den
def Trans_to_Selfloop(g):
    for i in g.nodes:
        g.add_edge(i,i,weight=Failure_Probability(g.degree(i)))
    return g
def Failure_Probability(k):
    if fpfunc=='(k+1)^-1':
        return 1/(k+1)
    if fpfunc=='e^-k':
        return np.exp(-k)
    if fpfunc=='1-k^-1':
        return 1-1/k
    if fpfunc=='lnk':
        return np.log(k)
time1=time.time()
N=120
k=4
l=np.log(k)
testg=nx.random_regular_graph(k,N)
for u, v in testg.edges():
    testg.edges[u, v]['weight']=1
print(GetBCRatio.Get_BCRatio(testg))
testg=Trans_to_Selfloop(testg)
print(GetBCRatio.Get_BCRatio(testg))
time2=time.time()
print(theorem_regular(N,k,l))
print(time2-time1)