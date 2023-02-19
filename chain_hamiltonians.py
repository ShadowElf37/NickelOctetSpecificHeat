"""
Breaks up large chains, constructs Hamiltonians, and computes Cv and Xm
"""

from collections import defaultdict
import numpy as np
import pickle
import chainlib
import spinlib
import nperrlog
import multiprocessing as mp

N = chainlib.Chain.N_from_ss

def agnostic_ge(x, y):
    return int(x) >= int(y)
def agnostic_le(x, y):
    return int(x) <= int(y)


def pool_c(H, T, dT, s, i):
    print(N(s), i, sep='/', end='\n')
    return s,spinlib.compute_C(H[s], T, dT) / N(s)

def pool_x(M, H, T, s, i):
    print(N(s), i, sep='/', end='\n')
    return s,spinlib.compute_X(M[s], H[s], T) / N(s) * 2

if __name__ == "__main__":
    BREAK_UP_CHAINS = False
    CALC_HAMILTONIANS = True

    if BREAK_UP_CHAINS:
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
    if CALC_HAMILTONIANS:
        examples, counts = pickle.load(open('final_structures', 'rb'))
        #totals = [sum(x.values()) for x in counts]

        spin_spaces = [spinlib.SpinSpace(i) for i in range(1, 11)]
        dT = 0.01
        T = np.arange(dT, 10, dT, dtype=np.float64)

        print('Generating %s Hamiltonians.' % len(examples))
        H = {}
        M = {}
        print('Finished: ', end='')
        for i, (s,chain) in enumerate(examples.items()):
            ss = spin_spaces[chain.N-1]
            H[s] = ss.make_hamiltonian(*chain.make_H_pairs())
            M[s] = ss.make_magnetization()**2
            print(i+1, end=',')
            if (i+1) % 20 == 0: print()
        print('Done.\nSpinning up a pool...')

        pool = mp.Pool(16)

        print('Calculating specific heat...')
        C_raw = pool.starmap(pool_c, ((H,T,dT,s,i) for i,s in enumerate(H.keys())))

        with open('sample_cvx_data', 'rb') as file:
            Tx, _, X = pickle.load(file)
        with open('sample_cvx_data_1', 'wb') as file:
            pickle.dump((T, Tx, dict(C_raw), X), file)
        exit()

        print('Calculating magnetic susceptibility...')
        X_raw = pool.starmap(pool_x, ((M,H,T,s,i) for i,s in enumerate(H.keys())))

        C = dict(C_raw)
        X = dict(X_raw) # [mu ** 2 / (kb * T)] +

        with open('sample_cvx_data', 'wb') as file:
            pickle.dump((T, C, X), file)

        import matplotlib.pyplot as plt
        plt.plot(T, C[0])
        plt.plot(T, X[0])
        plt.show()

