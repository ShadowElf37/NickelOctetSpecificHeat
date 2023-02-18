"""
Creates a hundred doped samples
"""

import multiprocessing as mp
import pickle
import numpy as np
from dopelib import *


LENGTH = 100
#COUNTS = np.linspace(0, (LENGTH-1) ** 3, 101, dtype=int)

def produce_sampleL(cnt):
    x = cnt/LENGTH**3
    print('Producing sample at x=%.2f' % x)
    raw_sample = produce_sample(LENGTH, cnt)
    return raw_sample

    sample = to_vertices(raw_sample)
    #if cnt == COUNTS[-1]:print(sample)
    return cnt, count_polymer_pct(sample)

if __name__ == '__main__':
    import time
    ACTUAL_SAMPLES = mp.Queue()

    COUNTS = np.linspace(0, (LENGTH-1) ** 3, 101, dtype=int)

    pool = mp.Pool(16)
    t0 = time.time()
    results = pool.map(produce_sampleL, COUNTS)
    t = time.time() - t0

    print('Finished in', t)

    #ACTUAL_SAMPLES.put('STOP')

    #ACTUAL_SAMPLES_LIST = list(iter(ACTUAL_SAMPLES.get, 'STOP'))

    print('Saving...')

    #with open('doping_data', 'wb') as f:
    #    pickle.dump(results, f)
    with open('samples', 'wb') as f:
        pickle.dump(results, f)

    print('Done!')