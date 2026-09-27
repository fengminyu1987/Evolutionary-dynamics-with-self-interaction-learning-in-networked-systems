# -*- coding: utf-8 -*-
"""
Created on Fri Jun  6 13:40:18 2025

@author: ziyan
"""
import numpy as np
import matplotlib.pyplot as plt
import math
plt.rcParams['axes.linewidth'] = 3
c1="#000080"
c2="#FFD700"
def bcrreg(N,k,lk):
    return (N*(k*k+3*k*lk+2*lk*lk)-2*(k+lk)*(k+lk))/(N*(2*lk*lk+k*lk+k)-2*(k+lk)*(k+lk))
def star(N,alpha, beta):
    numerator = 3 * (1 + alpha) * (-1 + N + beta) * (-(8/3) + N**3 * (4/3 + (11 * alpha)/3 + alpha**2) + (14 * beta)/3 - (4 * beta**2)/3 + alpha**2 * (-2 + (8 * beta)/3) + alpha * (-6 + 10 * beta - (8 * beta**2)/3) + N**2 * (-(16/3) + 3 * beta + alpha**2 * (-4 + (4 * beta)/3) + alpha * (-(40/3) + 7 * beta)) + N * (20/3 + alpha**2 * (5 - 4 * beta) - (23 * beta)/3 + (4 * beta**2)/3 + alpha * (47/3 - 17 * beta + (8 * beta**2)/3)))
    denominator = 2 * (beta * (-N**3 + N * (-1 - beta/2) + N**2 * (2 - beta/2) + beta) + alpha**3 * (2 + N**4 - 5 * beta + 4 * beta**2 + N**3 * (-5 + (5 * beta)/2) + N * (-7 + (25 * beta)/2 - 6 * beta**2) + N**2 * (9 - 10 * beta + 2 * beta**2)) + alpha * (4 + N**4 - 11 * beta + 12 * beta**2 - 2 * beta**3 + N**3 * (-7 + (5 * beta)/2) + N**2 * (15 - 16 * beta + 5 * beta**2) + N * (-13 + (49 * beta)/2 - 17 * beta**2 + 2 * beta**3)) + alpha**2 * (8 + 4 * N**4 - 22 * beta + 21 * beta**2 - 4 * beta**3 + N**3 * (-20 + 12 * beta) + N**2 * (36 - 46 * beta + (27 * beta**2)/2) + N * (-28 + 56 * beta - (69 * beta**2)/2 + 4 * beta**3)))
    return numerator / denominator
def star_HH(N,alpha,gamma):
    numerator = (
        3 * (1 + alpha) * (N + gamma) * (
            N**4 * (alpha**3 + (16/3) * alpha**2 + (25/3) * alpha + (10/3)) +
            N**3 * (
                alpha**3 * ((4/3) * gamma - 2) +
                alpha**2 * (10 * gamma - (25/3)) +
                alpha * (20 * gamma - (34/3)) +
                (29/3) * gamma - (13/3)
            ) +
            N**2 * (
                alpha**3 * (-(7/3) * gamma + (1/3)) +
                alpha**2 * ((11/3) * gamma**2 - (38/3) * gamma + 1) +
                alpha * ((40/3) * gamma**2 - (64/3) * gamma + 1) +
                (26/3) * gamma**2 - 9 * gamma + (1/3)
            ) +
            N * (
                alpha**3 * (-(1/3) * gamma + (2/3)) +
                alpha**2 * (-(10/3) * gamma**2 - gamma + 2) +
                alpha * ((8/3) * gamma**3 - (32/3) * gamma**2 - gamma + 2) +
                3 * gamma**3 - (16/3) * gamma**2 - (1/3) * gamma + (2/3)
            ) +
            gamma * (
                (4/3) * alpha**3 + 4 * alpha**2 +
                alpha * (-(5/3) * gamma**2 + 4) +
                (1/3) * gamma**3 - gamma**2 + (4/3)
            )
        )
    )
    denominator = (
        2 * (
            N**5 * (alpha**4 + (11/2) * alpha**3 + 9 * alpha**2 + 5 * alpha + 2) +
            N**4 * (
                alpha**4 * ((5/2) * gamma - (5/2)) +
                alpha**3 * ((33/2) * gamma - (23/2)) +
                alpha**2 * ((61/2) * gamma - 19) +
                alpha * ((33/2) * gamma - (25/2)) +
                6 * gamma - 4
            ) +
            N**3 * (
                alpha**4 * (2 * gamma**2 - (11/2) * gamma + 2) +
                alpha**3 * ((37/2) * gamma**2 - 29 * gamma + (15/2)) +
                alpha**2 * (44 * gamma**2 - (107/2) * gamma + 11) +
                alpha * (29 * gamma**2 - 36 * gamma + 7) +
                (21/2) * gamma**2 - 12 * gamma + (3/2)
            ) +
            N**2 * (
                alpha**4 * (-(7/2) * gamma**2 + 2 * gamma - (1/2)) +
                alpha**3 * gamma * ((11/2) * gamma**2 - (49/2) * gamma + (15/2)) +
                alpha**2 * ((47/2) * gamma**3 - (111/2) * gamma**2 + (23/2) * gamma + (7/2)) +
                alpha * (21 * gamma**3 - (79/2) * gamma**2 + 7 * gamma + 5) +
                9 * gamma**3 - 14 * gamma**2 + gamma + 2
            ) +
            N * (
                alpha**4 * (-(1/2) * gamma**2 + gamma - 1) +
                alpha**3 * (-5 * gamma**3 - 2 * gamma**2 + 6 * gamma - (11/2)) +
                alpha**2 * (4 * gamma**4 - (41/2) * gamma**3 - (5/2) * gamma**2 + 13 * gamma - (21/2)) +
                alpha * (6 * gamma**4 - 16 * gamma**3 - (5/2) * gamma**2 + 12 * gamma - (17/2)) +
                (7/2) * gamma**4 - (13/2) * gamma**3 - (3/2) * gamma**2 + 4 * gamma - (5/2)
            ) +
            alpha**4 * (2 * gamma**2 + 1) +
            alpha**3 * (8 * gamma**2 - (3/2) * gamma + 4) +
            alpha**2 * (-(5/2) * gamma**4 + (25/2) * gamma**2 - (9/2) * gamma + 6) +
            alpha * ((1/2) * gamma**5 - 2 * gamma**4 - (1/2) * gamma**3 + 9 * gamma**2 - (9/2) * gamma + 4) +
            (1/2) * gamma**5 - gamma**4 - (1/2) * gamma**3 + (5/2) * gamma**2 - (3/2) * gamma + 1
        )
    )
    return numerator / denominator
def ceiling_fan(N,epsilon,beta):
    numerator = (2 + epsilon) * (-1 + N + beta) * (
    78 + 4 * N**2 * (9 + 11 * epsilon + 2 * epsilon**2) +
    epsilon**2 * (17 - 22 * beta) - 122 * beta + 24 * beta**2 +
    epsilon * (83 - 125 * beta + 22 * beta**2) + 
    N * (-114 + 68 * beta + 5 * epsilon**2 * (-5 + 2 * beta) + epsilon * (-127 + 74 * beta))
)
    denominator = (
    -57 + N**3 * (9 + 25 * epsilon + 34 * epsilon**2 + 6 * epsilon**3) +
    157 * beta - 155 * beta**2 + 15 * beta**3 +
    epsilon**3 * (-12 + 29 * beta - 22 * beta**2) +
    epsilon**2 * (-76 + 201 * beta - 177 * beta**2 + 22 * beta**3) +
    epsilon * (-109 + 294 * beta - 273 * beta**2 + 28 * beta**3) +
    N**2 * (25 * (-3 + beta) + 48 * epsilon**2 * (-3 + 2 * beta) + 2 * epsilon**3 * (-12 + 7 * beta) +
            3 * epsilon * (-53 + 26 * beta)) +
    N * (123 - 182 * beta + 47 * beta**2 + epsilon**3 * (30 - 43 * beta + 10 * beta**2) +
         3 * epsilon**2 * (62 - 99 * beta + 32 * beta**2) + 3 * epsilon * (81 - 124 * beta + 35 * beta**2))
)
    return numerator / denominator
def Regular_Phase_Diagram(N, krange,lkrange):
    K, LK = np.meshgrid(krange, lkrange)
    Z = bcrreg(N, K, LK)
    plt.figure(figsize=(8, 6))
    plt.contourf(K, LK, Z > 0, levels=1, colors=[c1,c2])
    plt.contour(K, LK, Z, levels=[0], colors="black", linewidths=4)
    #plt.plot(krange[int(N/2+20):], [(N*(4-k)+math.sqrt(N*N*k*k+8*k*N+16*(k-N)*(k*N-k-N)))/(4*(k-1)) for k in krange[int(N/2):]],lw=2,color='black')
    plt.xlabel("$k$",fontsize=25)
    plt.ylabel("$\\ell(k)$",fontsize=25)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.tick_params(axis='both', # 应用于x轴和y轴
    which='major', # 应用于主刻度线
    direction='in', # 刻度线向内外都显示
    length=6, # 刻度线长度
    width=2, # 刻度线宽度
    colors='black', # 刻度线颜色
    labelcolor='black', # 标签颜色
    bottom=True, # 显示x轴的底部刻度线
    top=True, # 不显示x轴的顶部刻度线
    left=True, # 显示y轴的左侧刻度线
    right=True)
    plt.text(10,1,'Cooperation',fontsize=35,color='grey')
    plt.text(40.5,0.05,'Spite',fontsize=35,color='grey')
    plt.savefig(".\\saves\\imgs\\Regular_Phase_{:d}.pdf".format(N),dpi=300,bbox_inches='tight')
    plt.show()
def Star_Phase_Diagram(N, alpharange, betarange):
    A,B=np.meshgrid(alpharange,betarange)
    Z=star(N,A,B)
    plt.figure(figsize=(8,8))
    plt.contourf(A, B, Z > 0, levels=1, colors=[c1,c2])
    #plt.imshow(Z,cmap=plt.cm.jet,extent=(0,0.4,0,3),aspect='auto')
    plt.contour(A, B, Z, levels=[0], colors="black", linewidths=2)
    plt.xlabel("$\\alpha$",fontsize=25)
    plt.ylabel("$\\beta$",fontsize=25)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.text(0.20,15,'Cooperation',fontsize=35,color='grey')
    plt.text(0.02,15,'Spite',fontsize=35,color='grey')
    plt.tick_params(axis='both', # 应用于x轴和y轴
    which='major', # 应用于主刻度线
    direction='in', # 刻度线向内外都显示
    length=6, # 刻度线长度
    width=2, # 刻度线宽度
    colors='black', # 刻度线颜色
    labelcolor='black', # 标签颜色
    bottom=True, # 显示x轴的底部刻度线
    top=True, # 不显示x轴的顶部刻度线
    left=True, # 显示y轴的左侧刻度线
    right=True) 
    plt.savefig(".\\saves\\imgs\\Star_Phase_{:d}.pdf".format(N),dpi=300,bbox_inches='tight')
    plt.show()
def Star_HH_Phase_Diagram(N, alpharange, gammarange):
    A,C=np.meshgrid(alpharange,gammarange)
    Z=star_HH(N,A,C)
    plt.figure(figsize=(8,8))
    #plt.contourf(A, C, Z > 0, levels=1, colors=["lightcoral", "lightblue"])
    plt.imshow(Z,cmap=plt.cm.jet,extent=(0,1,0,8),aspect='auto')
    #plt.contour(A, B, Z, levels=[0], colors="black", linewidths=1)
    cb=plt.colorbar()
    cb.set_label(label='$(\\frac{b}{c})^*_{\\mathbf{\\mathcal{L}}}$',rotation=0,fontsize=25)
    cb.ax.tick_params(labelsize=25)
    plt.xlabel("$\\alpha$",fontsize=25)
    plt.ylabel("$\\gamma$",fontsize=25)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.tick_params(axis='both', # 应用于x轴和y轴
    which='major', # 应用于主刻度线
    direction='in', # 刻度线向内外都显示
    length=6, # 刻度线长度
    width=2, # 刻度线宽度
    colors='black', # 刻度线颜色
    labelcolor='black', # 标签颜色
    bottom=True, # 显示x轴的底部刻度线
    top=True, # 不显示x轴的顶部刻度线
    left=True, # 显示y轴的左侧刻度线
    right=True)
    plt.savefig(".\\saves\\imgs\\Star_HH_Phase_{:d}.pdf".format(N),dpi=300,bbox_inches='tight')
    plt.show()
def Ceiling_Fan_Phase_Diagram(N, epsilonrange, betarange):
    E,B=np.meshgrid(epsilonrange, betarange)
    Z=ceiling_fan(N,E,B)
    plt.figure(figsize=(8,8))
    plt.contourf(E, B, Z > 0, levels=1, colors=[c1,c2])
    #plt.imshow(Z,cmap=plt.cm.jet,extent=(0,1,0,8),aspect='auto')
    plt.contour(E, B, Z, levels=[0], colors="black", linewidths=2)
    plt.xlabel("$\\epsilon$",fontsize=25)
    plt.ylabel("$\\beta$",fontsize=25)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.tick_params(axis='both', # 应用于x轴和y轴
    which='major', # 应用于主刻度线
    direction='in', # 刻度线向内外都显示
    length=6, # 刻度线长度
    width=2, # 刻度线宽度
    colors='black', # 刻度线颜色
    labelcolor='black', # 标签颜色
    bottom=True, # 显示x轴的底部刻度线
    top=True, # 不显示x轴的顶部刻度线
    left=True, # 显示y轴的左侧刻度线
    right=True)
    if N==3:
        plt.text(0.8,1.5,'Cooperation',fontsize=35,color='grey')
        plt.text(0.25,0.5,'Spite',fontsize=35,color='grey')
    if N==5:
        plt.text(0.75,1.5,'Cooperation',fontsize=35,color='grey')
        plt.text(0.08,0.25,'Spite',fontsize=35,color='grey',rotation=90)
    plt.savefig(".\\saves\\imgs\\Ceiling_Fan_Phase_{:d}.pdf".format(N),dpi=300,bbox_inches='tight')
    plt.show()
#Regular_Phase_Diagram(50, np.linspace(3,50,5000),np.linspace(0,1.5,5000))
#Star_Phase_Diagram(25, np.linspace(0,0.4,5000),np.linspace(0,30,5000))
Star_HH_Phase_Diagram(15, np.linspace(0,1,1000),np.linspace(0,8,1000))
#Ceiling_Fan_Phase_Diagram(3, np.linspace(0,2,5000),np.linspace(0,2,5000))