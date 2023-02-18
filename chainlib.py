"""
The basic Chain builder and analyzer
"""

import numpy as np
from collections import defaultdict
from copy import deepcopy
from typing import Iterable, List, Dict
import multiprocessing as mp
import pickle
import time

def flatten(l):
    return [i for a in l for i in a]

#print(diffs)

def index_add(p, dp):
    return tuple(np.array(p)+np.array(dp))
def where(condition):  # this returns actual points you can iterate over
    points_bad = np.array(np.where(condition))
    return points_bad.transpose()


def bool_setadd(s: set, obj):
    if obj not in s:
        s.add(obj)
        return True
    return False

class Chain:
    SHIFTS = tuple(map(np.array, ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))))
    SHIFTS_INVERSE = tuple(map(np.array, ((-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1))))

    def __init__(self, *bonds, default_base_point=(0, 0, 0)):
        self.force_octamer = False

        self.edges = set()
        self.vertices = set()
        self.origin = default_base_point

        if bonds:
            for bond in bonds:
                self.bond(*bond, update_structure=False)
        else:
            self.vertices.add(default_base_point)
        self.structure = self.analyze()

    def __repr__(self):
        return f'<{self.N}-Chain {self.structstr()}>'

    def structstr(self):
        return "".join(map(str, self.structure.values()))

    #@property
    #def vertices(self):
    #    return set(flatten(self.edges))

    @property
    def N(self):
        return len(self.vertices)

    def make_H_pairs(self):
        pairs = []
        vertices = list(self.vertices)
        for v1,v2 in self.edges:
            pairs.append((vertices.index(v1), vertices.index(v2)))
        return tuple(pairs)

    def copy(self):
        return deepcopy(self)

    def e_exists(self, v1: tuple, v2: tuple) -> bool:
        return (v1,v2) in self.edges or (v2,v1) in self.edges
    def v_exists(self, v: tuple) -> bool:
        return v in self.vertices

    def analyze(self) -> dict:  # sums all bond counts at each vertex to acquire the relevant structure
        counts = {k: 0 for k in range(min(self.N, 7))}
        #if self.force_octamer:
        #    counts[3] = 8
        #    return counts
        for v in self.vertices:
            counts[sum(1 for e in self.edges if v in e)] += 1
        return counts

    def bond(self, p1: tuple, p2: tuple, update_structure=True):  # bonds two vertices, or doesn't if they're already bonded
        # NOTE: This returns true if the SECOND POINT was added
        if self.e_exists(p1, p2):
            return False
        set.add(self.edges, (p1, p2))
        self.vertices.add(p1)
        bonded = bool_setadd(self.vertices, p2)
        if update_structure: self.structure = self.analyze()
        return bonded


    def branch(self, available_shifts=SHIFTS) -> list:  # return a list of all possible chains with 1 additional connected vertex compared to the current structure
        new_chains = []
        # add new vertices
        for v1 in self.vertices:
            for dv in available_shifts:
                v2 = index_add(v1, dv)
                if not self.e_exists(v1, v2): # if the bond doesn't already exist, this is a new structure we add to the list
                    new_chains.append(c := self.copy())
                    c.bond(v1, v2)

                    # bond possible cycles by checking for already-existing vertices around the new one
                    for dv in available_shifts:
                        v3 = index_add(v2, dv)
                        if c.v_exists(v3) and not c.e_exists(v2, v3):
                            new_chains.append(d := c.copy())
                            d.bond(v2, v3)
        return new_chains

    @staticmethod
    def find_unique_structures(chains: Iterable) -> list:  # finds 1 example chain for every unique structure in the list
        uniques = []
        for chain in chains:
            if not any(chain.structure == c.structure for c in uniques):
                uniques.append(chain)
        return uniques

    @staticmethod
    def count_unique_structures(chains: Iterable) -> defaultdict:  # counts instances of each unique structure type
        counts = defaultdict(int)
        for chain in chains:
            counts[chain.structstr()] += 1
        return counts

    @staticmethod
    def any_with_structure(chains: Iterable, structurestring: str):  # gets an example from the list that has the specified structure
        s = {i:int(v) for i,v in enumerate(structurestring)}
        for chain in chains:
            if chain.structure == s:
                return chain
        return None

    @staticmethod
    def all_with_structure(chains: Iterable, structurestring: str):  # returns every chain in list with specified structure
        s = {i: int(v) for i, v in enumerate(structurestring)}
        new = []
        for chain in chains:
            if chain.structure == s:
                new.append(chain)
        return new

    @staticmethod
    def nbranch(base_chains: list, ntimes: int) -> list:  # branches n times and returns all chains in between, including base chain
        #t0 = time.time_ns()
        branched_chains = map(Chain.branch, base_chains.copy())
        #t1 = time.time_ns() - t0
        #if t1: print('Branching performance', t1/1000000)

        if ntimes == 1:
            return base_chains + flatten(branched_chains)

        print(f'Running nbranch at n={ntimes}...')
        for chain in branched_chains:
            base_chains += Chain.nbranch(chain, ntimes-1)
        return base_chains

    @staticmethod
    def mp_nbranch(pool, base_chains: list, ntimes: int) -> list:  # nbranch with mp support
        #t0 = time.time_ns()
        branched_chains = pool.map(Chain.branch, base_chains.copy())
        #t1 = time.time_ns()-t0
        #if t1: print('Branching performance', t1/1000000)

        if ntimes == 1:
            return base_chains + flatten(branched_chains)

        #print(f'Running nbranch at n={ntimes}...')
        for chain in branched_chains:
            base_chains += Chain.mp_nbranch(pool, chain, ntimes - 1)
        return base_chains


base_chains = [Chain()]
if __name__ == "__main__":

    print('Starting up pool...')
    pool = mp.Pool(16)
    print('Here we go!')
    built = Chain.mp_nbranch(pool, base_chains, 6)
    #built = Chain.nbranch(base_chains, 5)
    print('Done!')
    us = Chain.find_unique_structures(built)

    print(us)

    counts = []
    for i in range(max(us, key=lambda s: s.N).N):
        counts.append(0)
        for s in us:
            if s.N == i+1:
                counts[-1] += 1
    print(counts)

    with open('structures', 'wb') as f:
        pickle.dump(us, f)





#for i, c in enumerate(built):
