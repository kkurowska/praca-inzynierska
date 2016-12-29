# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 18:25:59 2016

@author: Kinia
"""
from math import sqrt

def countZ(n1, n2, k1, k2, variance):
    if ((k1 == 0 and k2 == 0) or variance == 0):
        Z = 0
    else:
        Z = (k1/n1 - k2/n2) / sqrt(variance)
    return Z