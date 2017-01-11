# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 23:58:00 2017

@author: Kinia
"""

def countarrayZ(n1, n2, x1, x2, variance):
    Z = np.zeros(np.shape(x1))
    b = (variance > 0) # bo inaczej wywali błąd, czy może zostać Z=0?
    b0 = (variance == 0)
    Z[b] = (x1[b]/n1 - x2[b]/n2) / np.sqrt(variance[b])
    Z[b0] = x1[b0]/n1 - x2[b0]/n2
    return Z
def countZ(n1, n2, k1, k2, variance):
    if (variance == 0): # czy mogę tak to robić?
        Z = 0
    else:
        Z = (k1/n1 - k2/n2) / sqrt(variance)
    return Z
    

import numpy as np
from scipy.stats import norm, hypergeom as hg, binom
import matplotlib.pyplot as plt
from math import floor, sqrt

alpha = 0.05

N1 = 100
M1 = 10
p1 = M1/N1
    
N2 = 100
M2 = 10
p2 = M2/N2

n=23

vectorK = np.array([1,0,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0])

K = np.sum(vectorK)

K1 = np.array([0,0,1,0,0,0,1,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0])
K2 = np.array([0,0,1,0,0,0,1,0,1,1,0,1,0,0,0,0,0,1,0,1,0,0])
K3 = np.array([0,1,0,0,1,1,0,1,1,1,0,1,0,1,1,1,1,0,0,0,1,0])
K4 = np.array([0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0])
K5 = np.array([0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0])

KK = np.array([np.sum(K1),np.sum(K2),np.sum(K3),np.sum(K4),np.sum(K5)])
invKK = n-KK

print(invKK)

counter = 5
bE = np.zeros(counter)

for i in range(counter):
    k1 = K
    k2 = invKK[i]    
    
    variance_k = ((N1 - n)/(n*(N1 - 1)) + (N2 - n)/(n*(N2 - 1))) * ((k1 + k2)/(n + n)) * (1 - (k1 + k2)/(n + n))  
    if (variance_k != 0):
        Z_k = countZ(n, n, k1, k2, variance_k)    
    
        estp = (k1 + k2)/(n + n)
        estM1 = floor(N1 * estp)
        estM2 = floor(N2 * estp) 
        Lx1 = max(0, estM1 - N1 + n)
        Lx2 = max(0, estM2 - N2 + n)
        Ux1 = min(n, estM1)
        Ux2 = min(n, estM2)
        
        sizex1 = Ux1 + 1 - Lx1
        sizex2 = Ux2 + 1 - Lx2 
        vectorX1 = np.arange(Lx1, Ux1 + 1).reshape((sizex1, 1)) # pionowy wektor
        vectorX2 = np.arange(Lx2, Ux2 + 1)
        arrayX1 = np.tile(vectorX1, (1, sizex2)) # wektory X1 w pionie
        arrayX2 = np.tile(vectorX2, (sizex1, 1)) # wektory X2 w poziomie
        
        hx1 = hg.pmf(vectorX1, N1, estM1, n) 
        hx2 = hg.pmf(vectorX2, N2, estM2, n)
        arrayHx1 = np.tile(hx1, (1, sizex2)) # wektory w pionie
        arrayHx2 = np.tile(hx2, (sizex1, 1)) # wektory w poziomie  
        arrayHx = arrayHx1 * arrayHx2
        
        variance_x = ((N1 - n)/(n*(N1 - 1)) + (N2 - n)/(n*(N2 - 1))) * ((arrayX1 + arrayX2)/(n + n)) * (1 - (arrayX1 + arrayX2)/(n + n))
        Z_x = countarrayZ(n, n, arrayX1, arrayX2, variance_x)  
        bZx = np.abs(Z_x) >= np.abs(Z_k) # indykator
        bx0 = variance_x == 0
        bZx[bx0] = np.abs(Z_x[bx0]) >= abs(k1/n - k2/n)
    #                bZx0 = (variance_x == 0) #w tym przypadku porownujemy estp1 i estp2
    #                bZx[bZx0] = (arrayX1[bZx0]/n1 != arrayX2[bZx0]/n)
        pvalueE = np.sum(arrayHx[bZx])
        bE[i] = pvalueE <= alpha # zapisuje się 0 lub 1
    elif (k1/n != k2/n):
        bE[i] = 1
            
            
print(bE)
        
