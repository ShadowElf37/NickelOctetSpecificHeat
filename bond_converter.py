"""
Test file
"""

import pickle
from chainlib import *
import numpy as np

with open('structures', 'rb') as f:
    chains = pickle.load(f)
    ss = []
    for chain in sorted(chains, key=lambda c: c.N):
        #print(f's{chain.N}.make_coupled_H{chain.make_H_pairs()}')
        ss.append(chain.structstr())
    print(ss)