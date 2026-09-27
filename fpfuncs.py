# -*- coding: utf-8 -*-
"""
Created on Fri May 17 16:48:09 2024

@author: zziya
"""
import numpy as np
import matplotlib.pyplot as plt
#fpfunc='(k+1)^-1'#(failure probability for replacement)
fpfuncs=['(k+1)^-1','e^-k','(k+1)^-2','1-k^-1']
Cs=['red', 'blue', 'green', 'black']
Ms=['^', '>', '<', 'v']
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
    if fpfunc=='(k+1)^-3':
        return 1/(k+1)**3
Ks=np.arange(1,30)
plt.figure(figsize=(8,8))
for fpfunc in fpfuncs:
    res=[]
    for k in Ks:
        res.append(Failure_Probability(k))
    plt.scatter(Ks,res,color=Cs[fpfuncs.index(fpfunc)],marker=Ms[fpfuncs.index(fpfunc)],s=90,label='{}'.format(fpfunc))
plt.legend(fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel('Failure Probability',fontsize=25)
plt.xlabel('$k$',fontsize=25)
plt.show()