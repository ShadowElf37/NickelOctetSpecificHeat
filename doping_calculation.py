import multiprocessing as mp
import pickle
import numpy as np
from doping import *


LENGTH = 100
#COUNTS = np.linspace(0, (LENGTH-1) ** 3, 101, dtype=int)

def produce_sampleL(cnt):
    x = cnt/LENGTH**3
    print('Producing sample at x=%.2f' % x)
    sample = produce_sample(LENGTH, cnt)
    #if cnt == COUNTS[-1]:print(sample)
    return cnt, count_polymer_pct(sample)

if __name__ == '__main__':
    COUNTS = np.linspace(0, (LENGTH-1) ** 3, 101, dtype=int)

    pool = mp.Pool(16)
    results = pool.map(produce_sampleL, COUNTS)

    print('Saving...')

    with open('doping_data.npy', 'wb') as f:
        pickle.dump(results, f)

    print('Done!')