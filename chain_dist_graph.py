"""
Graphs distribution of chain types in the sample
"""

import matplotlib.pyplot as plt
import numpy as np
import pickle
from collections import defaultdict
from chainlib import Chain


total_data = [defaultdict(int) for _ in range(101)]

for i in range(9):
    with open(f'real_structures{i}', 'rb') as f:
        #data = [pair[0] for pair in pickle.load(f)]
        data = pickle.load(f)
        for i in range(101):
            for key in data[i].keys():
                total_data[i][key] += data[i][key]

#print(data[1])

data = total_data

X = np.linspace(0, 1, len(data))

#print(data)
N = Chain.N_from_ss
Chain.N_from_ss = lambda k: 1

total_chains = [sum(v * (Chain.N_from_ss(k) or 1) for (k,v) in sample.items()) for sample in data]
by_structure = [{k: v * (Chain.N_from_ss(k))/total_chains[i] for (k,v) in (sample.items())} for (i,sample) in enumerate(data)]
by_ion_count = [[0]*28 for sample in data]
for i,sample in enumerate(data):
    for k,v in sample.items():
        by_ion_count[i][N(k)] += v*(Chain.N_from_ss(k) or 1)/total_chains[i]


#nice_structures = [sum([frac for structure,frac in by_structure[i].items() if Chain.N_from_ss(structure) >= 8 and Chain.N_from_ss(structure, 3) > Chain.N_from_ss(structure, 0, 3)]) for i in range(len(X))]

fig, (ax1) = plt.subplots(1,1)





from scipy.ndimage.filters import gaussian_filter1d

gaussian_filter1d = lambda x, sigma: x

data = [gaussian_filter1d([(by_ion_count[i][j]) for i in range(len(X))], sigma=0.5) for j in range(0, 28)]

#ax1.plot(X, [by_ion_count[i][0] for i in range(len(X))], label="Still Conducting", linewidth=1.5)
ax1.set_xlim([0,0.3])
ax1.set_ylim([0,0.03])
ax1.plot(X, data[1], label="Isolated", linewidth=1.5)
ax1.plot(X, data[2], label="Dimers", linewidth=1.5)
ax1.plot(X, data[3], label="Trimers", linewidth=1.5)
ax1.plot(X, data[4], label="Tetramers", linewidth=1.5)
#ax1.plot(X, data[4], label="Pentamers", linewidth=1.5)

##ax1.plot(X, 2*X**2)
"""
ax1.plot(X, [by_ion_count[i][0] for i in range(len(X))], label="Still Conducting", linewidth=1.5)
ax1.plot(X, [sum(by_ion_count[i][1:]) for i in range(len(X))], label="Isolated+", linewidth=1.5)
ax1.plot(X, [sum(by_ion_count[i][2:]) for i in range(len(X))], label="Dimers+", linewidth=1.5)
ax1.plot(X, [sum(by_ion_count[i][3:]) for i in range(len(X))], label="Trimers+", linewidth=1.5)
ax1.plot(X, [sum(by_ion_count[i][4:]) for i in range(len(X))], label="Tetramers+", linewidth=1.5)
ax1.plot(X, [sum(by_ion_count[i][5:]) for i in range(len(X))], label="Pentamers+", linewidth=1.5)
ax1.plot(X, [sum(by_ion_count[i][6:]) for i in range(len(X))], label="Hexamers+", linewidth=1.5)
ax1.plot(X, [sum(by_ion_count[i][7:]) for i in range(len(X))], label="Septamers+", linewidth=1.5)
ax1.plot(X, [sum(by_ion_count[i][8:]) for i in range(len(X))], label="Octamers+", linewidth=1.5)"""

#for j in range(1,28):
#    ax1.plot(X, [sum(by_ion_count[i][j:]) for i in range(len(X))], label=f"{j}+", linewidth=1.5)





ax1.set_xlabel('Doping Concentration (x)')
ax1.set_ylabel('Fraction of total structures')

ax1.legend()

plt.suptitle('Doped Structure Distribution')

plt.show()