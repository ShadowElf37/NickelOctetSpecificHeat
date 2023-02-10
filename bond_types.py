import numpy as np
from collections import defaultdict

diffs = tuple(map(np.array, ((1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1))))
#print(diffs)

def safe_add(p, dp):
    return tuple(np.array(p)+np.array(dp))
def where(condition):
    points_bad = np.array(np.where(condition))
    return points_bad.transpose()

def analyze(w):
    counts = defaultdict(int)
    #print('#', w)
    for p, ion in np.ndenumerate(w):
        if ion:
            #print(p, ion)
            current_count = 0
            for dp in diffs:
                # bounds check
                if any(np.greater(np.array(p)+np.array(dp), L-1)) or any(np.less(np.array(p)+np.array(dp), 0)):
                    continue
                #print(safe_add(p, dp), w[safe_add(p, dp)])
                current_count += w[safe_add(p, dp)]
            counts[current_count] += 1

    return counts

def branch(w):
    structures = []
    unique_worlds = []

    for p in where(w == 1):
        #print(p)
        for dp in diffs:
            shifted = safe_add(p, dp)
            #print(shifted)
            # bounds check
            if any(np.greater(np.array(p) + np.array(dp), L - 1)) or any(np.less(np.array(p) + np.array(dp), 0)):
                continue
            # already was a vertex there check
            if w[shifted]:
                continue

            w[shifted] = 1

            structure = analyze(w)
            if structure not in structures:
                structures.append(structure)
                unique_worlds.append(w.copy())

            w[shifted] = 0

    return structures, unique_worlds


def iterate(base_structures, base_worlds, ntimes):
    #print(ntimes)
    if ntimes == 0:
        return base_structures, base_worlds
    #print('Made it through')

    total_structures, total_worlds = base_structures.copy(), base_worlds.copy()
    for world in base_worlds:
        #print('Doing world')
        branched = branch(world)
        #print(branched, world)
        new_structures, new_worlds = iterate(*branch(world), ntimes-1)
        #print('Branched!')
        total_structures += new_structures
        total_worlds += new_worlds
    return total_structures, total_worlds


L = 3
base_world = np.zeros((L, L, L), dtype=np.int8)
base_world[L // 2, L // 2, L // 2] = 1
base_structure = analyze(base_world)

#print(base_structure, base_world)

#print(branch(base_world))

for s, w in zip(*iterate([base_structure], [base_world], 3)):
    print(dict(s))
    print(w.tolist())