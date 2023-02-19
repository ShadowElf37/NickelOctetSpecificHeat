"""
Graphs distribution of chain types in the sample
"""

import matplotlib.pyplot as plt
import numpy as np
import pickle
from collections import defaultdict
from chainlib import Chain


with open('real_structures', 'rb') as f:
    data = [pair[0] for pair in pickle.load(f)]

X = np.linspace(0, 1, len(data))

#print(data)

total_chains = [sum(sample.values()) for sample in data]
by_structure = [{k: v/total_chains[i] for (k,v) in (sample.items())} for (i,sample) in enumerate(data)]
by_ion_count = [defaultdict(int) for sample in data]
for i,sample in enumerate(data):
    for k,v in sample.items():
        by_ion_count[i][Chain.N_from_ss(k)] += v/total_chains[i]


nice_structures = [sum([frac for structure,frac in by_structure[i].items() if Chain.N_from_ss(structure) >= 8 and Chain.N_from_ss(structure, 3) > Chain.N_from_ss(structure, 0, 3)]) for i in range(len(X))]

fig, (ax1) = plt.subplots(1,1)


ax1.plot(X, [by_ion_count[i][1] for i in range(len(X))], label="Monomers", linewidth=1.5)
ax1.plot(X, [by_ion_count[i][2] for i in range(len(X))], label="Dimers", linewidth=1.5)
ax1.plot(X, [by_ion_count[i][3] for i in range(len(X))], label="Trimers", linewidth=1.5)
ax1.plot(X, [by_ion_count[i][4] for i in range(len(X))], label="Tetramers", linewidth=1.5)
ax1.plot(X, [by_ion_count[i][5] for i in range(len(X))], label="Pentamers", linewidth=1.5)
ax1.plot(X, [by_ion_count[i][6] for i in range(len(X))], label="Hexamers", linewidth=1.5)
ax1.plot(X, [by_ion_count[i][7] for i in range(len(X))], label="Septamers", linewidth=1.5)
ax1.plot(X, [sum(by_ion_count[i][j] for j in range(8, 28)) for i in range(len(X))], label="Octamers+", linewidth=1.5)
#ax1.plot(X, [by_ion_count[i][9] for i in range(len(X))], label="Nonamers", linewidth=1.5)
#ax1.plot(X, [by_ion_count[i][10] for i in range(len(X))], label="Decamer+", linewidth=1.5)
#ax1.plot(X, [by_ion_count[i][11] for i in range(len(X))], label="??+", linewidth=1.5)
ax1.plot(X, nice_structures, label='Nice Octamers', color='red')





ax1.set_xlabel('Doping Concentration (x)')
ax1.set_ylabel('Fraction of total structures')

ax1.legend()

plt.suptitle('Doped Structure Distribution')

plt.show()