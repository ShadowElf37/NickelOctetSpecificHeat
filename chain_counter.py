"""
Counts up all the types of chains in a raw sample
"""

import pickle
from collections import defaultdict

import numpy as np
import dopelib
import chainlib
import random
import multiprocessing as mp

#print(test)
def manhattan_dist(p1, p2):
    return np.sum(np.abs(np.array(p1) - np.array(p2)))
def manhattan_dist3(p1, p2):
    return manhattan_dist(p1, p2) <= 3
def manhattan_within_radius(p1, p2, radius):
    return all(np.less_equal(np.abs(np.array(p1) - np.array(p2)), np.ones_like(p1)*radius))

def branch(sample, point): # return all vertices that the point is bonded to
    new_vertices = []
    for (edge_index,) in chainlib.where(sample[point] == 0): # 1 marks a bond, so this gets all the bond indices
        shifted_point = chainlib.index_add(point, chainlib.Chain.SHIFTS[edge_index]) # shift to get the next vertex
        if all(sample.shape[0] > i > -1 for i in shifted_point) and not any(sample[shifted_point]): # check it's in the sample's bounds so we dont get index errors or looping around
            new_vertices.append(shifted_point)
    return new_vertices

#print(branch(test, (0,0,0)))

def probe_point(sample, origin, n=3, vertex_limit=27, octamer_if_limit_reached_at_n=3, radius_limit=1):
    chain = chainlib.Chain(default_base_point=origin)
    points = [origin]

    #print(sample[origin])

    if any(sample[origin]):
        return None

    for i in range(n):
        new_points = []
        STOP = False

        for p in points:
            #print(p)
            for new_p in branch(sample, p):
                if not manhattan_within_radius(new_p, origin, radius_limit):
                    continue
                if not chain.e_exists(p, new_p) and chain.bond(p, new_p, update_structure=False):
                    new_points.append(new_p)
                if chain.N >= vertex_limit:
                    #if i <= octamer_if_limit_reached_at_n: # if it's growing really fast its prob octamer
                    #    chain.force_octamer = True
                    STOP = True
                    break
            if STOP:
                break

        points = new_points

        if STOP or not points:
            break

    chain.structure = chain.analyze()
    return chain


def analyze_sample(sample, sample_i=None, npoints=1000000//27):
    if sample_i:
        print('Analyzing sample', sample_i)
    to_hit = dopelib.shuffled_point_box(sample.shape[0])[:npoints]
    structure_counts = defaultdict(int)
    examples = {}
    for _ in range(len(to_hit)):
        new_chain = probe_point(sample, to_hit.pop())
        ss = '0' if new_chain is None else new_chain.structstr()
        structure_counts[ss] += 1
        #if ss not in examples.keys():
        #    examples[ss] = new_chain
    return structure_counts#, examples

#print(analyze_sample(test))



if __name__ == "__main__":
    pool = mp.Pool(16)

    with open('samples', 'rb') as f:
        samples = pickle.load(f)

    structures = pool.starmap(analyze_sample, zip(samples, range(len(samples))))

    with open('real_structures0', 'wb') as f:
        pickle.dump(structures, f)

    print('Done!')