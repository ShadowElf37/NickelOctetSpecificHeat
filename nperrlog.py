"""
Import this to make numpy shut up about overflows
"""

import numpy as np

WARNS = 0

def __error_hit(*e, **kwe):
    global WARNS
    WARNS += 1

np.seterrcall(__error_hit)
np.seterr('call')