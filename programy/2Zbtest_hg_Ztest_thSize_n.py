# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 18:41:22 2017

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
M1 = 5
#n1 = 5 # probka
p1 = M1/N1
    
N2 = 100
M2 = 5
n2 = 60 # probka
p2 = M2/N2

vectorSizeZ = np.zeros(n2)    
vectorSizeZb = np.zeros(n2)  
    
for n in range(1, n2+1):
    
    n1=n # n1 = n2
                    
    L1 = max(0, M1 - N1 + n1)
    L2 = max(0, M2 - N2 + n)
    U1 = min(n1, M1)
    U2 = min(n, M2)
    
    size1 = U1 + 1 - L1
    size2 = U2 + 1 - L2 
    vectorK1 = np.arange(L1, U1 + 1).reshape((size1, 1)) # pionowy wektor
    vectorK2 = np.arange(L2, U2 + 1)
    arrayK1 = np.tile(vectorK1, (1, size2)) # wektory K1 w pionie
    arrayK2 = np.tile(vectorK2, (size1, 1)) # wektory K2 w poziomie
    
    h1 = hg.pmf(vectorK1, N1, M1, n1) 
    h2 = hg.pmf(vectorK2, N2, M2, n)
    arrayH1 = np.tile(h1, (1, size2)) # wektory w pionie
    arrayH2 = np.tile(h2, (size1, 1)) # wektory w poziomie  
    arrayH = arrayH1 * arrayH2
    
    variance_k = ((N1 - n1)/(n1*(N1 - 1)) + (N2 - n)/(n*(N2 - 1))) * ((arrayK1 + arrayK2)/(n1 + n)) * (1 - (arrayK1 + arrayK2)/(n1 + n))
    Z_k = countarrayZ(n1, n, arrayK1, arrayK2, variance_k)
  
    # test Z (ze skończoną poprawką)
    quantile = norm.ppf(1-alpha/2)
    bZ = np.abs(Z_k) > quantile # indykator
    b0 = (variance_k == 0) #w tym przypadku porownujemy estp1 i estp2
    bZ[b0] = (arrayK1[b0]/n1 != arrayK2[b0]/n)
    sizeZ = np.sum(arrayH[bZ])
    vectorSizeZ[n-1] = sizeZ
    
    # test bez skończonej poprawki
    
    bZb = np.zeros((size1, size2))
    
    for r in range(size1):
        for c in range(size2):
            
            if (variance_k[r,c] != 0):

                estpK = (arrayK1[r,c] + arrayK2[r,c])/(n1 + n)
                
                sizex1 = n1 + 1
                sizex2 = n + 1 
                vectorX1 = np.arange(sizex1).reshape((sizex1, 1)) # pionowy wektor
                vectorX2 = np.arange(sizex2)
                arrayX1 = np.tile(vectorX1, (1, sizex2)) # wektory X1 w pionie
                arrayX2 = np.tile(vectorX2, (sizex1, 1)) # wektory X2 w poziomie
                
                bx1 = binom.pmf(vectorX1, n1, estpK) 
                bx2 = binom.pmf(vectorX2, n, estpK)
                arrayBx1 = np.tile(bx1, (1, sizex2))
                arrayBx2 = np.tile(bx2, (sizex1, 1)) 
                arrayBx = arrayBx1 * arrayBx2
                
                estp = (arrayX1 + arrayX2)/(n1 + n)
                
                varianceb_x = estp*(1 - estp)*(1/n1 + 1/n)
                Zb_x = countarrayZ(n1, n, arrayX1, arrayX2, varianceb_x)  
                bZx = np.abs(Zb_x) >= np.abs(Z_k[r,c]) # indykator
                bx0 = varianceb_x == 0
                bZx[bx0] = np.abs(Zb_x[bx0]) >= abs(arrayK1[r,c]/n1 - arrayK2[r,c]/n)
#                bZx0 = (variance_x == 0) #w tym przypadku porownujemy estp1 i estp2
#                bZx[bZx0] = (arrayX1[bZx0]/n1 != arrayX2[bZx0]/n)
                pvalueZb = np.sum(arrayBx[bZx])
                bZb[r,c] = pvalueZb <= alpha # zapisuje się 0 lub 1
            elif (arrayK1[r,c]/n1 != arrayK2[r,c]/n):                
                bZb[r,c] = 1
                print("coś")
      
    bZb = (bZb == 1)
    sizeZb = np.sum(arrayH[bZb])
    vectorSizeZb[n-1] = sizeZb
                

vectorX = np.arange(1, n2+1)
vectorAplha = np.tile(alpha, np.size(vectorX))

fig = plt.figure()
plt.plot(vectorX, vectorSizeZ)
plt.plot(vectorX, vectorSizeZb)
plt.plot(vectorX, vectorAplha, '.r')
#plt.axis([0, n2, 0, 0.1])
plt.grid(True)
plt.xlabel('$n$', fontsize=14)
plt.ylabel('p-stwo bledu I rodzaju')
#title = "$n_1=" + str(n1) + "," + "p=" + str(p2) + "$"
title = "$p=" + str(p2) + "$"
plt.title(title, fontsize=14)
plt.legend(["test Z", "test Zb", "alpha=0.05"], loc=4, fancybox=True, framealpha=0.5)
#plt.ylim([0,0.1])