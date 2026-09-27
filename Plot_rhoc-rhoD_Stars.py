# -*- coding: utf-8 -*-
"""
Created on Fri May 17 16:32:26 2024

@author: zziya
"""
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import random
import GetBCRatio
import json
network_type='ceiling_fan'
k=4
n=15
Cs=['#EAB080','#0001A1','#037F77','#C5272D']#Cs=['red', 'blue', 'green', 'black']
Ms=['^', '>', '<', 'v']
Pairs=[(1,2),(0.5,1),(0.2,1),(1.5,4)]
plt.figure(figsize=(8,8))
for (p,q) in Pairs:
    filename1='.\\saves\\df\\C_{}_{:d}_{:.2f}_{:.2f}.json'.format(network_type,n,p,q)
    filename2='.\\saves\\df\\D_{}_{:d}_{:.2f}_{:.2f}.json'.format(network_type,n,p,q)
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
    plt.axvline(THEOREM,color=Cs[Pairs.index((p,q))],linestyle='--',lw=2)
    z1=np.polyfit(Bf, n*np.array(rhoc_d), 1)
    p1=np.poly1d(z1)
    plt.scatter(Bf, n*np.array(rhoc_d), color=Cs[Pairs.index((p,q))],marker=Ms[Pairs.index((p,q))],s=300, label='$\\epsilon={:.1f},\\beta={:.1f}$'.format(p,q))
    plt.plot(np.linspace(min(Bf),max(Bf),15), p1(np.linspace(min(Bf),max(Bf),15)),color=Cs[Pairs.index((p,q))], lw=2)
plt.axhline(0, color = 'purple', alpha=1, linestyle='--',lw=3)
plt.legend(fontsize=20)
plt.ylabel('$N\\times(\\rho_C-\\rho_D)$',fontsize=35)
#plt.legend(fontsize=25)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('$b/c$',fontsize=35)
plt.grid()
plt.arrow(GetBCRatio.ceiling_fan(int(n),0,0), -0.02, 0, 0.01,width=0.4, head_width=0.8, head_length=0.01, fc='grey', ec='grey',alpha=0.7)
plt.ylim([-0.05,0.05])
plt.savefig(".\\saves\\imgs\\{}_{:d}.pdf".format(network_type,n), dpi=500, bbox_inches='tight')
plt.show()