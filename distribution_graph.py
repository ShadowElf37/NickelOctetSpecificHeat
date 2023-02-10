import matplotlib.pyplot as plt
import numpy as np
import pickle
from collections import defaultdict


with open('doping_data', 'rb') as f:
    data = pickle.load(f)

D = [list() for _ in range(7)]

X = np.linspace(0, 1, len(data))

for x in data:
    ddict = defaultdict(int)
    ddict.update(dict(x[1]))
    for i in range(7):
        D[i].append(ddict[i])

D = np.array(list(map(np.array, D)))


fig, (ax1) = plt.subplots(1,1)
#fig.set_size_inches(12, 5)

#print('$', C2, A2, M2, sep='%%\n')

ax1.plot(X, D[0], label="0 bonds (monomer)", linewidth=1.5)
ax1.plot(X, D[1], label="1 bond (dimer)", linewidth=1.5)
ax1.plot(X, D[2], label="2 bonds (tetramer)", linewidth=1.5)
ax1.plot(X, sum(D[3:]), label="3+ bonds (otcamer)", linewidth=1.5)

ax1.set_xlabel('Doping Concentration (x)')
ax1.set_ylabel('Fraction of total ions')

ax1.legend()

plt.suptitle('Doped Bond Distribution')

plt.show()