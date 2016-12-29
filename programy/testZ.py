# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 13:58:35 2016

@author: Kinia
"""

import numpy as np
from math import sqrt, floor
from scipy.stats import norm, hypergeom as hg #cdf - dystybuanta
import matplotlib.pyplot as plt

def countZ(n1, n2, k1, k2, variance):
    if ((k1 == 0 and k2 == 0) or variance == 0):
        Z = 0
    else:
        Z = (k1/n1 - k2/n2) / sqrt(variance)
    return Z
    
#def pvalueE(x1, n1, M1, N1, x2, n2, M2, N2, Z_k, Z_x):
    

counter = 1000 # proby monte carlo
alpha = 0.05

N1 = 100 # cała populacja
M1 = 10 # ilość z daną cechą
n1 = 60 # probka
p1 = M1/N1
population1 = np.append(np.ones(M1), np.zeros(N1-M1))
np.random.shuffle(population1)
    
N2 = 100 # cała populacja
M2 = 10        # ilość z daną cechą
n2 = 60 # probka
p2 = M2/N2
population2 = np.append(np.ones(M2), np.zeros(N2-M2))
np.random.shuffle(population1)

vectorExactSizeZ = np.array([])    
vectorExactSizeE = np.array([])
    
for n in range(1, n2+1):

    vectorZ_k = np.array([])
    vectorpvalueE = np.array([])
        
    for i in range(counter):
    
        sample1 = np.random.choice(population1, n)  
        k1 = np.sum(sample1)    
            
        sample2 = np.random.choice(population2, n)   
        k2 = np.sum(sample2)            
        
        variance_k = ((N1 - n)/(n*(N1 - 1)) + (N2 - n)/(n*(N2 - 1))) * ((k1 + k2)/(n + n)) * (1 - (k1 + k2)/(n + n))  
        
        Z_k = countZ(n, n, k1, k2, variance_k)        
        vectorZ_k = np.append(vectorZ_k, Z_k)
        
        estp = (k1 + k2)/(n + n)
        estM1 = floor(N1 * estp)
        estM2 = floor(N2 * estp) 
        L1 = max(0, estM1 - N1 + n)
        L2 = max(0, estM2 - N2 + n)
        U1 = min(n, estM1)
        U2 = min(n, estM2)
        
#        vectorX1 = np.arange(L1, U1 + 1)
#        vectorX2 = np.arange(L2, U2 + 1)
        hg1 = hg(N1, estM1, n) 
        hg2 = hg(N2, estM2, n)
        
        vectorZ_x = np.array([])
        
        pvalueE = 0
        
        for x1 in range(L1, U1 + 1):
            h1 = hg1.pmf(x1)
            for x2 in range(L2, U2 + 1):                     
                variance_x = ((N1 - n)/(n*(N1 - 1)) + (N2 - n)/(n*(N2 - 1))) * ((x1 + x2)/(n + n)) * (1 - (x1 + x1)/(n + n))
                Z_x = countZ(n, n, x1, x2, variance_x)
                if Z_x >= Z_k:
                    h2 = hg2.pmf(x2)
                    pvalueE += h1 * h2
        
        vectorpvalueE = np.append(vectorpvalueE, pvalueE)           
    
    vectorpvalueZ = 2 * (1 - norm.cdf(np.abs(vectorZ_k)))
    exactSizeZ = np.sum(vectorpvalueZ < alpha) / counter
    vectorExactSizeZ = np.append(vectorExactSizeZ, exactSizeZ)
    
    exactSizeE = np.sum(vectorpvalueE < alpha) / counter
    vectorExactSizeE = np.append(vectorExactSizeE, exactSizeE)    
        

vectorX = np.arange(1, n2+1)

fig = plt.figure()
plt.plot(vectorX, vectorExactSizeZ)
plt.plot(vectorX, vectorExactSizeE)
#plt.vlines(vectorX, 0, vectorP1)
#plt.axis([0, n, 0, maxY])
plt.grid(True)
plt.xlabel('$n$', fontsize=14)
plt.ylabel('Exact size')
title = "$p=" + str(p2) + "$"
plt.title(title, fontsize=14)
#plt.xticks(vectorX)
