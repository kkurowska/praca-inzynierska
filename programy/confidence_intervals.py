# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 00:34:54 2017

@author: Kinia
"""

import numpy as np
from scipy.stats import norm, hypergeom as hg, binom
import matplotlib.pyplot as plt
from math import floor
   
def countarrayZ(n1, n2, x1, x2, variance):
    Z = np.zeros(np.shape(x1))
    b = (variance > 0) # bo inaczej wywali błąd, czy może zostać Z=0?
    if (np.any(variance<0)):
        print("wariancja<0!!")
    b0 = (variance == 0)
    Z[b] = (x1[b]/n1 - x2[b]/n2) / np.sqrt(variance[b])
    Z[b0] = x1[b0]/n1 - x2[b0]/n2
    return Z
    
alpha = 0.05

N1 = 100
M1 = 10
n1 = 30 # probka
p1 = M1/N1
    
N2 = 100
M2 = 10
n2 = 10 # probka
p2 = M2/N2

vectorSizeE = np.zeros(n2)    
vectorSizeZb = np.zeros(n2)  

L1 = max(0, M1 - N1 + n1)
L2 = max(0, M2 - N2 + n2)
U1 = min(n1, M1)
U2 = min(n2, M2)

print(L1,U1,L2,U2)

size1 = U1 + 1 - L1
size2 = U2 + 1 - L2 
vectorK1 = np.arange(L1, U1 + 1).reshape((size1, 1)) # pionowy wektor
vectorK2 = np.arange(L2, U2 + 1)
arrayK1 = np.tile(vectorK1, (1, size2)) # wektory K1 w pionie
arrayK2 = np.tile(vectorK2, (size1, 1)) # wektory K2 w poziomie

h1 = hg.pmf(vectorK1, N1, M1, n1) 
h2 = hg.pmf(vectorK2, N2, M2, n2)
arrayH1 = np.tile(h1, (1, size2)) # wektory w pionie
arrayH2 = np.tile(h2, (size1, 1)) # wektory w poziomie  
arrayH = arrayH1 * arrayH2

variance_k = ((N1 - n1)/(n1*(N1 - 1)) + (N2 - n2)/(n2*(N2 - 1))) * ((arrayK1 + arrayK2)/(n1 + n2)) * (1 - (arrayK1 + arrayK2)/(n1 + n2))
Z_k = countarrayZ(n1, n2, arrayK1, arrayK2, variance_k)


bE = np.zeros((size1, size2))

for r in range(size1):
    for c in range(size2):
        
        if (variance_k[r,c] != 0):
            
             # test E

            estp = (arrayK1[r,c] + arrayK2[r,c])/(n1 + n2)
            estM1 = floor(N1 * estp)
            estM2 = floor(N2 * estp) 
            Lx1 = max(0, estM1 - N1 + n1)
            Lx2 = max(0, estM2 - N2 + n2)
            Ux1 = min(n1, estM1)
            Ux2 = min(n2, estM2)
            
            sizex1 = Ux1 + 1 - Lx1
            sizex2 = Ux2 + 1 - Lx2 
            vectorX1 = np.arange(Lx1, Ux1 + 1).reshape((sizex1, 1)) # pionowy wektor
            vectorX2 = np.arange(Lx2, Ux2 + 1)
            arrayX1 = np.tile(vectorX1, (1, sizex2)) # wektory X1 w pionie
            arrayX2 = np.tile(vectorX2, (sizex1, 1)) # wektory X2 w poziomie
            
            hx1 = hg.pmf(vectorX1, N1, estM1, n1) 
            hx2 = hg.pmf(vectorX2, N2, estM2, n2)
            arrayHx1 = np.tile(hx1, (1, sizex2)) # wektory w pionie
            arrayHx2 = np.tile(hx2, (sizex1, 1)) # wektory w poziomie  
            arrayHx = arrayHx1 * arrayHx2
            
            variance_x = ((N1 - n1)/(n1*(N1 - 1)) + (N2 - n2)/(n2*(N2 - 1))) * ((arrayX1 + arrayX2)/(n1 + n2)) * (1 - (arrayX1 + arrayX2)/(n1 + n2))
            Z_x = countarrayZ(n1, n2, arrayX1, arrayX2, variance_x)  
            bZx = np.abs(Z_x) >= np.abs(Z_k[r,c]) # indykator
            bx0 = variance_x == 0
            bZx[bx0] = np.abs(Z_x[bx0]) >= abs(arrayK1[r,c]/n1 - arrayK2[r,c]/n2)
#                bZx0 = (variance_x == 0) #w tym przypadku porownujemy estp1 i estp2
#                bZx[bZx0] = (arrayX1[bZx0]/n1 != arrayX2[bZx0]/n)
            pvalueE = np.sum(arrayHx[bZx])
            bE[r,c] = pvalueE <= alpha # zapisuje się 0 lub 1
            
        elif (arrayK1[r,c]/n1 != arrayK2[r,c]/n2):                
            bE[r,c] = 1
            

bE = bE == 1
P = (arrayK1/n1 - arrayK2/n2)
KP = P[bE]
KPD = KP[KP>0]
KPU = KP[KP<0]

a = np.max(KPU)
b = np.min(KPD)

print(a, b)