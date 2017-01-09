# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 13:25:11 2017

@author: Kinia
"""


import numpy as np
from math import sqrt, floor
from scipy.stats import norm, hypergeom as hg, binom
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
M1 = 10 # ilość z daną cechą
n1 = 10 # probka
p1 = M1/N1
population1 = np.append(np.ones(M1), np.zeros(N1-M1))
np.random.shuffle(population1)
    
N2 = 100 # cała populacja
M2 = 10        # ilość z daną cechą
n2 = 60 # probka
p2 = M2/N2
population2 = np.append(np.ones(M2), np.zeros(N2-M2))
np.random.shuffle(population2)

vectorExactSizeZ = np.array([])    
vectorExactSizeZb = np.array([])
    
for n in range(1, n2+1):
    
    n1 = n

    exactSizeZ = 0
    exactSizeZb = 0
    vectorpvalueZb = np.array([])
        
    for i in range(counter):
    
        sample1 = np.random.choice(population1, n1)  
        k1 = np.sum(sample1)    
            
        sample2 = np.random.choice(population2, n)   
        k2 = np.sum(sample2)            
        
        # test Z
        variance_k = ((N1 - n1)/(n1*(N1 - 1)) + (N2 - n)/(n*(N2 - 1))) * ((k1 + k2)/(n1 + n)) * (1 - (k1 + k2)/(n1 + n))  
        if (variance_k != 0):
            Z_k = countZ(n1, n, k1, k2, variance_k)    
            pvalueZ = 2 * (1 - norm.cdf(abs(Z_k)))
            if (pvalueZ < alpha):
                exactSizeZ += 1  
        elif (k1/n1 != k2/n):
                exactSizeZ += 1

        # test Zb       
    
        estp = (k1 + k2)/(n1 + n)
        varianceb_k = estp*(1 - estp)*(1/n1 + 1/n)  
    
        if (varianceb_k != 0):
            Zb_k = countZ(n1, n, k1, k2, varianceb_k)    
            
            size1 = n1 + 1
            size2 = n + 1
            
            vectorX1 = np.arange(size1).reshape((size1, 1))
            vectorX2 = np.arange(size2)
            arrayX1 = np.tile(vectorX1, (1, size2))
            arrayX2 = np.tile(vectorX2, (size1, 1))
            
            b1 = binom.pmf(vectorX1, n1, p1) 
            b2 = binom.pmf(vectorX2, n, p2)
            arrayB1 = np.tile(b1, (1, size2))
            arrayB2 = np.tile(b2, (size1, 1)) 
            arrayB = arrayB1 * arrayB2
    
            estp_x = (arrayX1 + arrayX2)/(n1 + n)
            varianceb_x = estp_x*(1 - estp_x)*(1/n1 + 1/n)  
       
            Zb_x = countarrayZ(n1, n, arrayX1, arrayX2, varianceb_x)       
            
            bZxZk = np.abs(Zb_x) >= np.abs(Zb_k) # indykator
            b0 = varianceb_x == 0
            bZxZk[b0] = np.abs(Zb_x[b0]) >= abs(k1/n1 - k2/n2)
            pvalueZb = np.sum(arrayB[bZxZk])
            if (pvalueZb < alpha):
                exactSizeZb += 1
                
        elif (k1/n1 != k2/n):
                exactSizeZb += 1
    
   # vectorpvalueZ = 2 * (1 - norm.cdf(np.abs(vectorZ_k)))
   # exactSizeZ = np.sum(vectorpvalueZ < alpha) / counter
    vectorExactSizeZ = np.append(vectorExactSizeZ, exactSizeZ/counter)
    vectorExactSizeZb = np.append(vectorExactSizeZb, exactSizeZb/counter)    
        

vectorX = np.arange(1, n2+1)
vectorAplha = np.tile(alpha, np.size(vectorX))

fig = plt.figure()
plt.plot(vectorX, vectorExactSizeZ)
plt.plot(vectorX, vectorExactSizeZb)
plt.plot(vectorX, vectorAplha, '.r')
#plt.axis([0, n2, 0, 0.1])
plt.grid(True)
plt.xlabel('$n_2$', fontsize=14)
plt.ylabel('p-stwo bledu I rodzaju')
title = "$n_1=" + str(n1) + "," + "p=" + str(p2) + "$"
plt.title(title, fontsize=14)
plt.legend(["test Z", "test Zb"], loc=4)
