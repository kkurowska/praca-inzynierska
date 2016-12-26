# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from math import *
# funkcja błędu erf(cc)
from scipy.stats import norm, hypergeom as hg
#ppf - gęstość, cdf - dystybuanta
#norm.ppf(0.95, loc=10, scale=2)
#hypergeom: pmf - gęstość
import matplotlib.pyplot as plt
import pickle
import random

M = 100 # cała populacja
n = 20 # ilość z daną cechą
N = 15 # probka
p = n/M

with open('populacja100.pickle', 'rb') as f:
    populacja = pickle.load(f)
    
m = np.mean(populacja)
s = np.std(populacja)

print(m,s)
 
probka = np.array([])

for i in populacja:
    a = random.random()
    if a <= p:
        probka = np.append(probka, populacja[i])
        
print(probka)


        
    

# wybrać n osobnikow z populacji
#
#dane = np.array(hg.rvs(M, n, N, size=10))
#print(dane)
#
#Z = norm.cdf(abs(dane))
#print (Z)
#
#E = hg.cdf(dane, M, n, N)
#print (E)