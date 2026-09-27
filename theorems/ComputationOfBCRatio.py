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
fpfunc='e^-k'#(failure probability for replacement)
realset=['dolphins','americanfootball','contiguous-usa','zebra','ca-sandi_auths','rt-retweet']
fpfuncs=['(k+1)^-1','e^-k','(k+1)^-2','1-k^-1']
network_type = 'dolphins'
filename = '.\\datas\\' + network_type + '.txt'
# Note: the function nx.read_edgelist('dolphins.txt') is alternative
G = nx.Graph()
with open(filename) as file:
    for line in file:
        head, tail = [str(x) for x in line.split()]
        G.add_edge(int(head), int(tail), weight=1)
# relabel
if 0 not in G.nodes():
    j = 0
    new_mapping = {}
    for i in G.nodes:
        new_mapping[i] = j
        j += 1
    G = nx.relabel_nodes(G, new_mapping)
G=GetBCRatio.Trans_to_Selfloop(G, fpfunc)
print(GetBCRatio.Get_BCRatio(G))