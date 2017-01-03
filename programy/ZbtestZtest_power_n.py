# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 23:23:39 2016

@author: Kinia
"""

import numpy as np
from scipy.stats import norm, hypergeom as hg, binom
import matplotlib.pyplot as plt
   
def countarrayZ(n1, n2, x1, x2, variance):
    Z = np.zeros(np.shape(x1))
    b = (variance > 0) # bo inaczej wywali błąd, czy może zostać Z=0?
    if (np.any(variance<0)):
        print("wariancja<0!!")
    Z[b] = (x1[b]/n1 - x2[b]/n2) / np.sqrt(variance[b])
    return Z
    
alpha = 0.05

N1 = 100 # cała populacja
vectorM1 = np.array([20, 30, 40]) # ilość z daną cechą
n1 = 20 # probka
vectorp1 = vectorM1/N1
    
N2 = 100 # cała populacja
M2 = 10  # ilość z daną cechą
n2 = 50 # probka
p2 = M2/N2
population2 = np.append(np.ones(M2), np.zeros(N2-M2))
np.random.shuffle(population2)

arrayPowerZ = np.zeros((3, n2))    
arrayPowerZb = np.zeros((3, n2))  
    
for n in range(1, n2+1):
    
    n1=n # n1 = n2
            
    for i in range(3):
        
        quantile = norm.ppf(1-alpha/2)
        
        M1 = vectorM1[i]
        p1 = vectorp1[i]        
        
        population1 = np.append(np.ones(vectorM1), np.zeros(N1-M1))
        np.random.shuffle(population1)
        
        # test Z ze skończoną poprawką                 
             
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
        
        variance = ((N1 - n1)/(n1*(N1 - 1)) + (N2 - n)/(n*(N2 - 1))) * ((arrayK1 + arrayK2)/(n1 + n)) * (1 - (arrayK1 + arrayK1)/(n1 + n))
      #  Z = countarrayZ(n1, n, arrayK1, arrayK2, variance)
        Z = countarrayZ(n1, n, arrayK1, arrayK2, variance)
        bZ = np.abs(Z) > quantile # indykator
        b0 = (variance == 0) #w tym przypadku porownujemy estp1 i estp2
        bZ[b0] = (arrayK1[b0]/n1 != arrayK2[b0]/n)
        powerZ = np.sum(arrayH[bZ])
        arrayPowerZ[i, n-1] = powerZ


        # test Z bez skończonej poprawki

        size1b = n1 + 1 # nie jestem pewna jaki powinien być rozmiar
        size2b = n + 1 # nie jestem pewna jaki powinien być rozmiar
        vectorX1 = np.arange(size1b).reshape((size1b, 1))
        vectorX2 = np.arange(size2b)
        arrayX1 = np.tile(vectorX1, (1, size2b))
        arrayX2 = np.tile(vectorX2, (size1b, 1))
        
        b1 = binom.pmf(vectorX1, n1, p1) 
        b2 = binom.pmf(vectorX2, n, p2)
        arrayB1 = np.tile(b1, (1, size2b))
        arrayB2 = np.tile(b2, (size1b, 1)) 
        arrayB = arrayB1 * arrayB2

        estp = (arrayX1 + arrayX2)/(n1 + n)
        varianceb = estp*(1 - estp)*(1/n1 + 1/n)  
   
        Zb = countarrayZ(n1, n, arrayX1, arrayX2, varianceb)
        bZb = np.abs(Zb) > quantile # indykator
        bb0 = (varianceb == 0) #w tym przypadku porownujemy estp1 i estp2
        bZb[bb0] = (arrayX1[bb0]/n1 != arrayX2[bb0]/n)
        powerZb = np.sum(arrayB[bZb])
        arrayPowerZb[i, n-1] = powerZb

vectorX = np.arange(1, n2+1)

fig = plt.figure()
color = ['b', 'g', 'r']
legZ1 = "test Z: $p_1=" + str(vectorp1[0]) + "$"
legZ2 = "$p_1=" + str(vectorp1[1]) + "$"
legZ3 = "$p_1=" + str(vectorp1[2]) + "$"
legZb1 = "test Zb: $p_1=" + str(vectorp1[0]) + "$"
legZb2 = "$p_1=" + str(vectorp1[1]) + "$"
legZb3 = "$p_1=" + str(vectorp1[2]) + "$"
legZ = [legZ1, legZ2, legZ3]
legZb = [legZb1, legZb2, legZb3]
for i in range (3):
    plt.plot(vectorX, arrayPowerZ[i], color[i], label = legZ[i])
for i in range (3):
    plt.plot(vectorX, arrayPowerZb[i], color[i], ls = '--', label = legZb[i])
plt.grid(True)
plt.xlabel("$n$", fontsize=14)
plt.ylabel("Power")
title = "$N_1=N_2=" + str(N1) + "$, $" + "p_2=" + str(p2) + "$"
plt.title(title, fontsize=14)
#plt.legend([legZ1, legZb1, legZ2, legZb2, legZ3, legZb3], loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))