# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 22:48:14 2016

@author: Kinia
"""

import pickle
from scipy.stats import norm
import numpy as np


M=100
populacja1 = np.array(norm.rvs(size=M))

with open('populacja100_norm.pickle', 'wb') as f:
    pickle.dump(populacja1, f, pickle.HIGHEST_PROTOCOL)
    
    
M=20
populacja2 = np.array(norm.rvs(size=M))

with open('populacja20_norm.pickle', 'wb') as f:
    # Pickle the 'data' dictionary using the highest protocol available.
    pickle.dump(populacja2, f, pickle.HIGHEST_PROTOCOL)