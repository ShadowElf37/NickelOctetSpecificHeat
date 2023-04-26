"""
Graphs distribution of chain types in the sample
"""

import matplotlib.pyplot as plt
import numpy as np
import pickle
from collections import defaultdict
from chainlib import Chain
import scipy

ncr = scipy.special.binom

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




"""
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
"""


#fac8 =
fac12 = 479001600

exact_data1 = [0, 0, 8, 192, 1536, 6720, 20160, 40320, 40320]
exact_pct1 = [e/np.prod(range(8, 8-i,-1)) for i,e in enumerate(exact_data1)]
exact_data2 =[0, 0, 0, 24, 1920, 40320, 463680, 3548160, 19514880,fac12/1/2/3,fac12/1/2,fac12/1,fac12]
exact_pct2 = [e/(np.prod(range(12, 12-i,-1) or 1)) for i,e in enumerate(exact_data2)]

#exact_data2p =[0,16,552,8448,84480,642240,3951360,19918080,fac12/1/2/3,fac12/1/2,fac12/1,fac12]
exact_data1p =[0,8,288,5184,62400,552960,3749760,19716480,fac12/1/2/3-4,fac12/1/2,fac12/1,fac12]
#exact_pct2p = [e/np.prod(range(12, 12-i-1,-1)) for i,e in enumerate(exact_data2p)]
exact_pct1p = [e/np.prod(range(12, 12-i,-1)) for i,e in enumerate(exact_data1p)]

print(*exact_pct1, sep=' ')
print(*exact_pct2, sep=' ')
print(*exact_pct1p, sep=' ')
print(*(np.array(exact_pct1p)**2), sep=' ')

#ax1.plot(np.array(list(range(1,9)))/8, exact_pct1, label="Exact Isolated+ (1)", linewidth=1.5)
#ax1.plot(np.array(list(range(1,13)))/12, exact_pct2p, label="Exact Isolated+ (2)", linewidth=1.5)
#ax1.plot(np.array(list(range(1,13)))/12, exact_pct2, label="Exact Dimer+ (2)", linewidth=1.5)
#ax1.plot(np.array(list(range(1,13)))/12, exact_pct1p, label="Exact Isolated+ (2)", linewidth=1.5)

#N = 10000
#attempt = np.array([ncr(36*N, k)*exact_pct2**k * (1-exact_pct2)**(36*N-k) for k in range(0, N, N//100)])

exact_isolated = sum(p*ncr(8, n)*np.power(X, n)*np.power(1-X, 8-n) for n,p in enumerate(exact_pct1))
exact_dimer = sum(p*ncr(12, n)*np.power(X, n)*np.power(1-X, 12-n) for n,p in enumerate(exact_pct2))
exact_isolated1 = sum(p*ncr(12, n)*np.power(X, n)*np.power(1-X, 12-n) for n,p in enumerate(exact_pct1p))


#ax1.plot(X, [by_ion_count[i][0] for i in range(len(X))], label="Simulated Conducting", linewidth=1.5)
ax1.plot(X, [sum(by_ion_count[i][2:]) for i in range(len(X))], label="Simulated Dimers+", linewidth=1.5)
ax1.plot(X, [sum(by_ion_count[i][1:]) for i in range(len(X))], label="Simulated Isolated+", linewidth=1.5)

ax1.plot(X, exact_dimer, label="Exact Dimer+", linewidth=1.5, linestyle='--')
ax1.plot(X, exact_isolated, label="Exact Isolated+", linewidth=1.5, linestyle='--')
#ax1.plot(X, exact_isolated1, label="Funny Isolated+", linewidth=1.5, linestyle='--')



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