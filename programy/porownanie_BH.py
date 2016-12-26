# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 18:06:20 2016

@author: Kinia
"""

from scipy.stats import binom, hypergeom
import numpy as np
import matplotlib.pyplot as plt

N = 20 #N
M = 5 #M
n = 10 #n

hg = hypergeom(N, M, n)
p = M/N
b = binom(n, p)

m1 = hg.mean()
v1 = hg.var()
print("hg.mean =", m1, "hg.var =", v1)
m2 = b.mean()
v2 = b.var()
print("b.mean =", m2, "b.var =", v2)

vectorX = np.arange(0, n+1)
vectorP1 = hg.pmf(vectorX)
vectorP2 = b.pmf(vectorX)

fig = plt.figure()
plt.plot(vectorX, vectorP1, 'o')
plt.vlines(vectorX, 0, vectorP1)
plt.axis([0, n, 0, 1])
plt.grid(True)
plt.xlabel('$k$', fontsize=14)
plt.ylabel('$P(X=k)$', fontsize=14)

fig = plt.figure()
plt.plot(vectorX, vectorP2, 'ro')
plt.vlines(vectorX, 0, vectorP2)
plt.axis([0, n, 0, 1])
plt.grid(True)
plt.xlabel('$k$', fontsize=14)
plt.ylabel('$P(X=k)$', fontsize=14)








