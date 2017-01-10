# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 22:21:11 2016

@author: Kinia
"""

from scipy.stats import binom, hypergeom
import numpy as np
import matplotlib.pyplot as plt

N = 20 #N
M = 5 #M
n = 20 #n
p = M/N


vectorV1 = np.array([])
vectorV2 = np.array([])

for i in range(n+1):

    hg = hypergeom(N, M, i)
    b = binom(i, p)
    vectorV1 = np.append(vectorV1, [hg.var()])
    vectorV2 = np.append(vectorV2, [b.var()])


vectorX = np.arange(0, n+1)

maxY = 4

fig = plt.figure()
plt.plot(vectorX, vectorV1, 'o')
plt.vlines(vectorX, 0, vectorV1)
plt.axis([0, n, 0, maxY])
plt.grid(True)
plt.xlabel('$n$', fontsize=14)
plt.ylabel('$var(X)$', fontsize=14)
plt.title('$X\sim \mathcal{H}(n,5,20)$', fontsize=14)
plt.xticks(vectorX)

fig = plt.figure()
plt.plot(vectorX, vectorV2, 'ro')
plt.vlines(vectorX, 0, vectorV2)
plt.axis([0, n, 0, maxY])
plt.grid(True)
plt.xlabel('$n$', fontsize=14)
plt.ylabel('$var(X)$', fontsize=14)
plt.title('$X\sim \mathcal{B}(n,0.25)$', fontsize=14)
plt.xticks(vectorX)