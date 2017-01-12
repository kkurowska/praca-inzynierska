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
n = 6
p = M/N

ba = -q*sqrt(p*(1-p)/n)
bb = q*sqrt(p*(1-p)/n)
print(ba,bb)

ha = -q*sqrt(p*(1-p)*(N-n)/(n*(N-1)))
hb = q*sqrt(p*(1-p)*(N-n)/(n*(N-1)))
print(ha,hb)