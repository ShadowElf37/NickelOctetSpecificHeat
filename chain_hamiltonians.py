"""
Breaks up large chains, constructs Hamiltonians, and computes Cv and Xm
"""

from collections import defaultdict

import numpy as np
import pickle
import chainlib

def agnostic_ge(x, y):
    return int(x) >= int(y)
def agnostic_le(x, y):
    return int(x) <= int(y)

with open('real_structures', 'rb') as f:
    data = pickle.load(f)
    counts = [pair[0] for pair in data]
    examples_by_x = [pair[1] for pair in data]
    X = np.linspace(0, 1, len(data))


print('Preparing data...')
examples: {str: chainlib.Chain} = {}
for x in examples_by_x:
    examples.update(x)
#print(examples)

#total_chains = [sum(sample.values()) for sample in counts]
large_chains = [chain for chain in examples.values() if chain.N > 10]
large_structures = [c.structstr() for c in large_chains]

decomposed = defaultdict(list)

print('Breaking up large chains...')
for chain in large_chains:
    origin = chain.origin#tuple(map(int, chain.origin))
    for cmp0 in (agnostic_ge, agnostic_le):
        for cmp1 in (agnostic_ge, agnostic_le):
            for cmp2 in (agnostic_ge, agnostic_le):
                octv = [v for v in chain.vertices if cmp0(v[0], origin[0]) and cmp1(v[1], origin[1]) and cmp2(v[2], origin[2])]
                #print(len(octv), octv)
                new_chain = chainlib.Chain(*(e for e in chain.edges if e[0] in octv and e[1] in octv))
                decomposed[chain.structstr()].append(new_chain)

    #print(chain)
    #print(decomposed)
    #exit()

for _,small_chains in decomposed.items():
    for sc in small_chains:
        examples[sc.structstr()] = sc

for s in large_structures:
    examples.pop(s)

print('Recounting structures...')
new_counts = []
for sample in counts:
    new_counts.append(defaultdict(int))
    for s,count in sample.items():
        if s not in large_structures:
            new_counts[-1][s] = count
        else:
            for decomp in decomposed[s]:
                new_counts[-1][decomp.structstr()] += 1

with open('final_structures', 'wb') as f:
    pickle.dump((examples, new_counts), f)

#print(new_counts)
print(len(examples))