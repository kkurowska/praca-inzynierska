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
n = 17 #n

hg = hypergeom(N, M, n)
p = M/N
b = binom(n, p)

m1 = hg.mean()
v1 = hg.var()
print("hg.mean =", m1, "hg.var =", v1)
m2 = b.mean()
v2 = b.var()
print("b.mean =", m2, "b.var =", v2)

L = max(0, M - N + n)
U = min(n, M)

vectorX1 = np.arange(L, U+1)
vectorX2 = np.arange(0, n+1)

vectorP1 = hg.pmf(vectorX1)
vectorP2 = b.pmf(vectorX2)

maxy = 0.5

fig = plt.figure()
plt.plot(vectorX1, vectorP1, 'o')
plt.vlines(vectorX1, 0, vectorP1)
plt.axis([0, n, 0, maxy])
plt.grid(True)
plt.xlabel('$k$', fontsize=14)
plt.ylabel('$P(X=k)$', fontsize=14)
plt.title('$X\sim h(17,5,20)$', fontsize=14)

fig = plt.figure()
plt.plot(vectorX2, vectorP2, 'ro')
plt.vlines(vectorX2, 0, vectorP2)
plt.axis([0, n, 0, maxy])
plt.grid(True)
plt.xlabel('$k$', fontsize=14)
plt.ylabel('$P(X=k)$', fontsize=14)
plt.title('$X\sim B(10,0.25)$', fontsize=14)








