# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 13:58:35 2016

@author: Kinia
"""

import numpy as np
from math import sqrt
from scipy.stats import norm #cdf - dystybuanta
import pickle
import random

counter = 10000 # proby monte carlo
alpha = 0.05

N1 = 100 # cała populacja
M1 = 20 # ilość z daną cechą
n1 = 20 # probka
p1 = M1/N1
with open('populacja1_N100_M20.pickle', 'rb') as f:
    population1 = pickle.load(f)
    
N2 = 100 # cała populacja
M2 = 20 # ilość z daną cechą
n2 = 20 # probka
p2 = M2/N2
with open('populacja2_N100_M20.pickle', 'rb') as f:
    population2 = pickle.load(f)
    
arrayZ = np.array([])
    
for i in range(counter):

    sample1 = np.zeros(n1)
    m = N1 - 1
    population = population1
    
    for i in range(n1):
        x = round(random.random()*m) #losujemy elem z populacji
        sample1[i] = population[x]
        m-=1
        population=np.delete(population, x) #usuwa wylosowany element
    
    k1 = np.sum(sample1)    
        
    sample2 = np.zeros(n2)
    m = N2 - 1
    population = population2
    
    for i in range(n2):
        x = round(random.random()*m) #losujemy elem z populacji
        sample2[i] = population[x]
        m-=1
        population=np.delete(population, x) #usuwa wylosowany element
    
    k2 = np.sum(sample2)            
    
    variance = ((N1 - n1)/(n1*(N1 - 1)) + (N2 - n2)/(n2*(N2 - 1))) * ((k1 + k2)/(n1 + n2)) * (1 - (k1 + k2)/(n1 + n2))  
    
    Z = (k1/n1 - k2/n2) / sqrt(variance)
    
    arrayZ = np.append(arrayZ, Z)

pvalue = 2 * (1 - norm.cdf(np.abs(arrayZ)))

exactSize = np.sum(pvalue < alpha) / counter

print("exact size =", exactSize)