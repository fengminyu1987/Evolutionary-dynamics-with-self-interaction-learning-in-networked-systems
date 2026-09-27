# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 13:06:50 2024

@author: zziya
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
xkeys=['bc','k','knn','density']
xlabels=['$b/c$','$\\langle k \\rangle$','$\\langle k_{nn} \\rangle$', 'Density']
ykey=['bcl']
nets=['ban','ern','hkn','ddn','wsn','nwn']
labels=['BAN','ERN','HKN','DDN','WSN','NWN']
funcs=['e^-k','lnk','1-k^-1','(k+1)^-1']
func='(k+1)^-1'
colors=['#7ABBDB','#84BA42','#682487','#DBB428','#D4562E','#9D9EA3']
flag=1
for xkey in xkeys:
    plt.figure(figsize=(8,8))
    for net in nets:
        data=pd.read_csv(".\\saves\\comparison\\{}_{}.csv".format(net,func))
        plt.scatter(data[xkey],data[ykey],s=10,color=colors[nets.index(net)],alpha=0.6)
        plt.scatter(data[xkey].iloc[0],data[ykey].iloc[0],s=45,color=colors[nets.index(net)],label=labels[nets.index(net)])
    if xkey=='bc':
        if func=='e^-k':
            plt.legend(fontsize=25,ncol=2)
        plt.ylabel('$(b/c)^*_{\\mathbf{\\mathcal{L}}}$',fontsize=35)#,rotation=0,labelpad=45)
        plt.plot([-100,1000],[-100,1000],color='red',alpha=0.7,lw=3)
        plt.xlim([0,140])
    if xkey=='k':
        plt.plot([-100,1000],[-100,1000],color='red',alpha=0.7,lw=3)
        plt.xlim([0,45])
    if xkey=='knn':
        plt.plot([-100,1000],[-100,1000],color='red',alpha=0.7,lw=3)
        plt.xlim([0,50])
    if xkey=='density':
        plt.xlim([0,0.4])
        plt.xticks(np.arange(0,0.41,step=0.1),fontsize=25)
    plt.ylim([0,120])
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.xlabel('{}'.format(xlabels[xkeys.index(xkey)]),fontsize=35)
    plt.savefig(".\\saves\\imgs\\comparison_{}_{}.png".format(func,xkey),dpi=300,bbox_inches='tight')
    plt.show()
if flag==1:
    plt.figure(figsize=(8,8))
    xkey='bc'
    ykey='bcl'
    for net in nets:
        data=pd.read_csv(".\\saves\\comparison\\{}_{}.csv".format(net,func))
        plt.scatter((data[xkey]+1)/(data[xkey]-1),(data[ykey]+1)/(data[ykey]-1),s=10,color=colors[nets.index(net)],alpha=0.6)
        plt.xlim([1,3.5])
        plt.ylim([1,3.5])
        plt.plot([-100,1000],[-100,1000],color='red',alpha=0.7,lw=3)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.xlabel('$\\sigma$',fontsize=35)
        plt.ylabel('$\\sigma_{\\mathbf{\\mathcal{L}}}$',fontsize=35)
        #plt.scatter(data[xkey].iloc[0],data[ykey].iloc[0],s=45,color=colors[nets.index(net)],label=labels[nets.index(net)])
    plt.savefig(".\\saves\\imgs\\comparison_{}_sigma.png".format(func),dpi=300,bbox_inches='tight')
    plt.show()