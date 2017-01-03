# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 16:22:20 2016

@author: Kinia
"""

import numpy as np
from math import sqrt, floor
from scipy.stats import norm, hypergeom as hg
import matplotlib.pyplot as plt

def countZ(p1, p2, variance):
    if (variance == 0): # czy mogę tak to robić?
        Z = 0
    elif(variance < 0):
        print("wariancja<0!!")  
    else:
        Z = (p1 - p2) / sqrt(variance)
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
M2 = 10  # ilość z daną cechą
n2 = 60 # probka
p2 = M2/N2
population2 = np.append(np.ones(M2), np.zeros(N2-M2))
np.random.shuffle(population2)

vectorExactSizeZ = np.array([])    
vectorExactSizeZb = np.array([])
    
for n in range(1, n2+1):
    
   # n1=n

    vectorZ = np.array([])
    vectorZb = np.array([])
        
    for i in range(counter):
    
        sample1 = np.random.choice(population1, n1)  
        k1 = np.sum(sample1)    
            
        sample2 = np.random.choice(population2, n)   
        k2 = np.sum(sample2)    

        estp1 = k1/n1
        estp2 = k2/n  
        estp = (n1*estp1 + n*estp2)/(n1 + n)
        
        # test Z ze skończoną poprawką                 
        variance = ((N1 - n1)/(n1*(N1 - 1)) + (N2 - n)/(n*(N2 - 1))) * ((k1 + k2)/(n1 + n)) * (1 - (k1 + k2)/(n1 + n))  
        
        Z = countZ(estp1, estp2, variance)        
        vectorZ = np.append(vectorZ, Z)

        # test Z bez skończonej poprawki
        varianceb = estp*(1 - estp)*(1/n1 + 1/n)  
        
        Zb = countZ(estp1, estp2, varianceb)        
        vectorZb = np.append(vectorZb, Zb)  
              
    
    vectorpvalueZ = 2 * (1 - norm.cdf(np.abs(vectorZ)))
    exactSizeZ = np.sum(vectorpvalueZ < alpha) / counter
    vectorExactSizeZ = np.append(vectorExactSizeZ, exactSizeZ)

    vectorpvalueZb = 2 * (1 - norm.cdf(np.abs(vectorZb)))
    exactSizeZb = np.sum(vectorpvalueZb < alpha) / counter
    vectorExactSizeZb = np.append(vectorExactSizeZb, exactSizeZb)


vectorX = np.arange(1, n2+1)
vectorAplha = np.tile(alpha, np.size(vectorX))

fig = plt.figure()
plt.plot(vectorX, vectorExactSizeZ)
plt.plot(vectorX, vectorExactSizeZb)
plt.plot(vectorX, vectorAplha, '.r')
#plt.axis([0, n2, 0, 0.1])
plt.grid(True)
plt.xlabel('$n_2$', fontsize=14)
plt.ylabel('Exact size')
title = "$n_1=" + str(n1) + "," + "p=" + str(p2) + "$"
plt.title(title, fontsize=14)
plt.legend(["test Z", "test Z bez skonczonej poprawki"])
