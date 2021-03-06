# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 19:47:36 2017

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

licznik=25

#N1 = 30 # cała populacja
#vectorM1 = np.array([21, 24, 27]) # ilość z daną cechą
N1 = 1000
vectorM1 = np.array([500, 400, 300])
n1 = 25 # probka
vectorp1 = vectorM1/N1
    
#N2 = 50 # cała populacja
#M2 = 30  # ilość z daną cechą
N2 = 3000
M2 = 1800
n2 = 50 # probka
p2 = M2/N2

arrayPowerE = np.zeros((3, licznik))    
arrayPowerZb = np.zeros((3, licznik))  
    
for j in range(1, licznik+1):
    
    n=j #n2
    
    n1=j # n1 = n2 ?
            
    for i in range(3):
        
        # test Z (ze skończoną poprawką)
        
        M1 = vectorM1[i]
        p1 = vectorp1[i]        
                    
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
        
               
         # test E

        bE = np.zeros((size1, size2))
        bZb = np.zeros((size1, size2))
        
        for r in range(size1):
            for c in range(size2):
                
                if (variance_k[r,c] != 0):
                    
                    #test E

                    estp = (arrayK1[r,c] + arrayK2[r,c])/(n1 + n)
                    estM1 = floor(N1 * estp)
                    estM2 = floor(N2 * estp) 
                    Lx1 = max(0, estM1 - N1 + n1)
                    Lx2 = max(0, estM2 - N2 + n)
                    Ux1 = min(n1, estM1)
                    Ux2 = min(n, estM2)
                    
                    sizex1 = Ux1 + 1 - Lx1
                    sizex2 = Ux2 + 1 - Lx2 
                    vectorX1 = np.arange(Lx1, Ux1 + 1).reshape((sizex1, 1)) # pionowy wektor
                    vectorX2 = np.arange(Lx2, Ux2 + 1)
                    arrayX1 = np.tile(vectorX1, (1, sizex2)) # wektory X1 w pionie
                    arrayX2 = np.tile(vectorX2, (sizex1, 1)) # wektory X2 w poziomie
                    
                    hx1 = hg.pmf(vectorX1, N1, estM1, n1) 
                    hx2 = hg.pmf(vectorX2, N2, estM2, n)
                    arrayHx1 = np.tile(hx1, (1, sizex2)) # wektory w pionie
                    arrayHx2 = np.tile(hx2, (sizex1, 1)) # wektory w poziomie  
                    arrayHx = arrayHx1 * arrayHx2
                    
                    variance_x = ((N1 - n1)/(n1*(N1 - 1)) + (N2 - n)/(n*(N2 - 1))) * ((arrayX1 + arrayX2)/(n1 + n)) * (1 - (arrayX1 + arrayX2)/(n1 + n))
                    Z_x = countarrayZ(n1, n, arrayX1, arrayX2, variance_x)  
                    bZx = np.abs(Z_x) >= np.abs(Z_k[r,c]) # indykator
                    bx0 = variance_x == 0
                    bZx[bx0] = np.abs(Z_x[bx0]) >= abs(arrayK1[r,c]/n1 - arrayK2[r,c]/n)
    #                bZx0 = (variance_x == 0) #w tym przypadku porownujemy estp1 i estp2
    #                bZx[bZx0] = (arrayX1[bZx0]/n1 != arrayX2[bZx0]/n)
                    pvalueE = np.sum(arrayHx[bZx])
                    bE[r,c] = pvalueE <= alpha # zapisuje się 0 lub 1
                

                    # test Zb                    
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
                    pvalueZb = np.sum(arrayBx[bZx])
                    bZb[r,c] = pvalueZb <= alpha # zapisuje się 0 lub 1
                elif (arrayK1[r,c]/n1 != arrayK2[r,c]/n):                
                    bZb[r,c] = 1
                    print("coś")
                    
        bE = (bE == 1)
        powerE = np.sum(arrayH[bE])
        arrayPowerE[i, j-1] = powerE     
          
        bZb = (bZb == 1)
        powerZb = np.sum(arrayH[bZb])
        arrayPowerZb[i, j-1] = powerZb
                

vectorX = np.arange(1, licznik+1)

fig = plt.figure()
color = ['b', 'g', 'r']
legE1 = "test E: $p_1=" + str(vectorp1[0]) + "$"
legE2 = "$p_1=" + str(vectorp1[1]) + "$"
legE3 = "$p_1=" + str(vectorp1[2]) + "$"
legZb1 = "test Zb: $p_1=" + str(vectorp1[0]) + "$"
legZb2 = "$p_1=" + str(vectorp1[1]) + "$"
legZb3 = "$p_1=" + str(vectorp1[2]) + "$"
legZb = [legZb1, legZb2, legZb3]
legE = [legE1, legE2, legE3]
for i in range (3):
    plt.plot(vectorX, arrayPowerE[i], color[i], label = legE[i])
for i in range (3):
    plt.plot(vectorX, arrayPowerZb[i], color[i], ls = '--', label = legZb[i])
plt.grid(True)
plt.xlabel("$n$", fontsize=14)
#plt.xlabel("$n_1(n_2=2n_1)$", fontsize=14)
plt.ylabel("moc testu")
title = "$N_1=" + str(N1) + "$, $" + "N_2=" + str(N2) + "$, $" + "p_2=" + str(p2) + "$"
#title = "$N_1=N_2=" + str(N2) + "$, $" + "p_2=" + str(p2) + "$"
plt.title(title, fontsize=14)
#plt.legend([legZ1, legZb1, legZ2, legZb2, legZ3, legZb3], loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True)
