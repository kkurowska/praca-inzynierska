# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 13:58:35 2016

@author: Kinia
"""

import numpy as np
from math import sqrt, floor
from scipy.stats import norm, hypergeom as hg
import matplotlib.pyplot as plt

def countZ(n1, n2, k1, k2, variance):
    if (variance == 0): # czy mogę tak to robić?
        Z = 0
    else:
        Z = (k1/n1 - k2/n2) / sqrt(variance)
    return Z
    
def countarrayZ(n1, n2, x1, x2, variance):
    Z = np.zeros(np.shape(x1))
    b = (variance > 0) # bo inaczej wywali błąd, czy może zostać Z=0?
    b0 = (variance == 0)
    Z[b] = (x1[b]/n1 - x2[b]/n2) / np.sqrt(variance[b])
    Z[b0] = x1[b0]/n1 - x2[b0]/n2
    return Z
      

counter = 1000 # proby monte carlo
alpha = 0.05

N1 = 100 # cała populacja
M1 = 5 # ilość z daną cechą
n1 = 10 # probka
p1 = M1/N1
population1 = np.append(np.ones(M1), np.zeros(N1-M1))
np.random.shuffle(population1)
    
N2 = 100 # cała populacja
M2 = 5        # ilość z daną cechą
n2 = 60 # probka
p2 = M2/N2
population2 = np.append(np.ones(M2), np.zeros(N2-M2))
np.random.shuffle(population2)

vectorExactSizeZ = np.array([])    
vectorExactSizeE = np.array([])
    
for n in range(1, n2+1):
    
    n1 = n

    exactSizeZ = 0
    exactSizeE = 0
    vectorpvalueE = np.array([])
        
    for i in range(counter):
    
        sample1 = np.random.choice(population1, n1)  
        k1 = np.sum(sample1)    
            
        sample2 = np.random.choice(population2, n)   
        k2 = np.sum(sample2)            
        
        variance_k = ((N1 - n1)/(n1*(N1 - 1)) + (N2 - n)/(n*(N2 - 1))) * ((k1 + k2)/(n1 + n)) * (1 - (k1 + k2)/(n1 + n))  
        if (variance_k != 0):
            
            # test Z
            Z_k = countZ(n1, n, k1, k2, variance_k)    
            pvalueZ = 2 * (1 - norm.cdf(abs(Z_k)))
            if (pvalueZ < alpha):
                exactSizeZ += 1        

            # test E       
        
            estp = (k1 + k2)/(n1 + n)
            estM1 = floor(N1 * estp)
            estM2 = floor(N2 * estp) 
            L1 = max(0, estM1 - N1 + n1)
            L2 = max(0, estM2 - N2 + n)
            U1 = min(n1, estM1)
            U2 = min(n, estM2)
            
            size1 = U1 + 1 - L1
            size2 = U2 + 1 - L2 
            vectorX1 = np.arange(L1, U1 + 1).reshape((size1, 1)) # pionowy wektor
            vectorX2 = np.arange(L2, U2 + 1)
            arrayX1 = np.tile(vectorX1, (1, size2)) # wektory X1 w pionie
            arrayX2 = np.tile(vectorX2, (size1, 1)) # wektory X2 w poziomie
            
            h1 = hg.pmf(vectorX1, N1, estM1, n1) 
            h2 = hg.pmf(vectorX2, N2, estM2, n)
            arrayH1 = np.tile(h1, (1, size2)) # wektory w pionie
            arrayH2 = np.tile(h2, (size1, 1)) # wektory w poziomie  
            arrayH = arrayH1 * arrayH2        
            
            vectorZ_x = np.array([])
            
            variance_x = ((N1 - n1)/(n1*(N1 - 1)) + (N2 - n)/(n*(N2 - 1))) * ((arrayX1 + arrayX2)/(n1 + n)) * (1 - (arrayX1 + arrayX2)/(n1 + n))
            Z_x = countarrayZ(n1, n, arrayX1, arrayX2, variance_x)
            bZxZk = np.abs(Z_x) >= np.abs(Z_k) # indykator
            b0 = variance_x == 0
            bZxZk[b0] = np.abs(Z_x[b0]) >= np.abs(k1/n1 - k2/n2)
            pvalueE = np.sum(arrayH[bZxZk])
            if (pvalueE < alpha):
                exactSizeE += 1
        
        else:
            if (k1/n1 != k2/n):
                exactSizeZ += 1
                exactSizeE += 1
    
   # vectorpvalueZ = 2 * (1 - norm.cdf(np.abs(vectorZ_k)))
   # exactSizeZ = np.sum(vectorpvalueZ < alpha) / counter
    vectorExactSizeZ = np.append(vectorExactSizeZ, exactSizeZ/counter)
    vectorExactSizeE = np.append(vectorExactSizeE, exactSizeE/counter)    
        

vectorX = np.arange(1, n2+1)
vectorAplha = np.tile(alpha, np.size(vectorX))

fig = plt.figure()
plt.plot(vectorX, vectorExactSizeZ)
plt.plot(vectorX, vectorExactSizeE)
plt.plot(vectorX, vectorAplha, '.r')
#plt.axis([0, n2, 0, 0.1])
plt.grid(True)
plt.xlabel('$n_2$', fontsize=14)
plt.ylabel('p-stwo bledu I rodzaju')
title = "$n_1=" + str(n1) + "," + "p=" + str(p2) + "$"
plt.title(title, fontsize=14)
plt.legend(["test Z", "test E"], loc=4)
