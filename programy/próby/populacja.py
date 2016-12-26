# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 14:01:44 2016

@author: Kinia
"""

import pickle
import numpy as np
import random


N = 100 # cała populacja
M = 20 # ilość z daną cechą
p = M/N

populacja1 = np.zeros(N)

licznik = M

# dokładnie 20 z daną cechą
while licznik > 0:
    x = round(random.random()*(N-1)) #losujemy od 0 do N-1
    if populacja1[x] == 0:
        populacja1[x] = 1
        licznik -= 1
        
# losujemy czy osobnik ma daną cechę z prawdop p
#==============================================================================
# for i in range(N):
#     a = random.random()
#     if a <= p:
#         populacja1 = np.append(populacja1, [1])
#     else:
#         populacja1 = np.append(populacja1, [0])
#==============================================================================



with open('populacja1_N100_M20.pickle', 'wb') as f:
    pickle.dump(populacja1, f, pickle.HIGHEST_PROTOCOL)
    
    
populacja2 = np.zeros(N)

licznik = M

while licznik > 0:
    x = round(random.random()*(N-1)) #losujemy od 0 do N-1
    if populacja2[x] == 0:
        populacja2[x] = 1
        licznik -= 1
    
with open('populacja2_N100_M20.pickle', 'wb') as f:
    pickle.dump(populacja2, f, pickle.HIGHEST_PROTOCOL)
    
