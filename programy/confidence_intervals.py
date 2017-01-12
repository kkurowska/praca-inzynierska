# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 19:08:38 2017

@author: Kinia
"""

import numpy as np
from scipy.stats import norm, hypergeom as hg, binom
import matplotlib.pyplot as plt
from math import floor, sqrt
   

alpha = 0.05

q = norm.ppf(1-alpha/2)

N = 20
M = 5
n = 10
p = M/N

k = 3

estp = k/n
print(estp)

ba = estp-q*sqrt(estp*(1-estp)/n)
bb = estp+q*sqrt(estp*(1-estp)/n)
print(ba,bb)

ha = estp-q*sqrt(estp*(1-estp)*(N-n)/(n*(N-1)))
hb = estp+q*sqrt(estp*(1-estp)*(N-n)/(n*(N-1)))
print(ha,hb)