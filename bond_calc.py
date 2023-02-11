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
def where(condition):
    points_bad = np.array(np.where(condition))
    return points_bad.transpose()


class Chain:
    SHIFTS = tuple(map(np.array, ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))))

    def __init__(self, *bonds):
        self.edges = set()
        self.vertices = set()

        if bonds:
            for bond in bonds:
                self.bond(*bond)
        else:
            self.vertices.add((0,0,0))
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

    def analyze(self) -> dict:
        counts = {k:0 for k in range(self.N)}
        for v in self.vertices:
            counts[len(tuple(filter(lambda e: v in e, self.edges)))] += 1
        return counts

    def bond(self, p1: tuple, p2: tuple, update_structure=True):
        self.edges.add((p1, p2))
        self.vertices.add(p1)
        self.vertices.add(p2)
        if update_structure: self.structure = self.analyze()

    def branch(self) -> list:
        new_chains = []
        # add new vertices
        for v1 in self.vertices:
            for dv in Chain.SHIFTS:
                v2 = index_add(v1, dv)
                if not self.e_exists(v1, v2):
                    new_chains.append(c := self.copy())
                    c.bond(v1, v2)
                    # bond possible cycles
                    for dv in Chain.SHIFTS:
                        v3 = index_add(v2, dv)
                        if c.v_exists(v3) and not c.e_exists(v2, v3):
                            new_chains.append(d := c.copy())
                            d.bond(v2, v3)
                                #print('Bonded a cycle!')
        return new_chains

    @staticmethod
    def find_unique_structures(chains: Iterable) -> list:
        uniques = []
        for chain in chains:
            if not any(chain.structure == c.structure for c in uniques):
                uniques.append(chain)
        return uniques

    @staticmethod
    def any_with_structure(chains: Iterable, structurestring: str):
        s = {i:int(v) for i,v in enumerate(structurestring)}
        for chain in chains:
            if chain.structure == s:
                return chain
        return None

    @staticmethod
    def all_with_structure(chains: Iterable, structurestring: str):
        s = {i: int(v) for i, v in enumerate(structurestring)}
        new = []
        for chain in chains:
            if chain.structure == s:
                new.append(chain)
        return new

    @staticmethod
    def nbranch(base_chains: list, ntimes: int) -> list:
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
    def mp_nbranch(pool, base_chains: list, ntimes: int) -> list:
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
