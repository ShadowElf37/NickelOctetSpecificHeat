"""
Tools to construct spin Hamiltonians and partition functions and stuff
"""


import numpy as np
from scipy.linalg import expm

J = -1
kb = 1
mu = -1

def tensor(*ops):
    total = ops[0]
    for op in ops[1:]:
        total = np.kron(total, op)
    return total

def check_symmetric(a, tol=1e-8):
    return np.all(np.abs(a-a.T) < tol)

def dim(H):
    return int(np.log2(H.shape[0]))


class SpinSpace:
    def __init__(self, n_particles, debug=False):
        # THIS IS SPIN 1/2!!!

        if debug: print('Preparing some arrays.')

        n = 2**n_particles
        self.n = n_particles

        self.zero = np.zeros(n)
        self.states = np.array([np.zeros(n) for _ in range(n)])
        for i, state in enumerate(self.states):
            state[i] = 1

        self.plus = np.array([np.identity(n) for _ in range(n_particles)])
        self.minus = np.array([np.identity(n) for _ in range(n_particles)])
        self.z = np.array([np.identity(n) for _ in range(n_particles)])

        canonical_id = np.identity(2)
        canonical_plus = np.array([[0,1],[0,0]])
        canonical_minus = np.array([[0,0],[1,0]])
        canonical_z = np.array([[1,0],[0,-1]])*0.5

        if debug: print('Generating S+...')

        for i, op in enumerate(self.plus):
            composition = [canonical_id.copy() for _ in range(n_particles)]
            composition[i] = canonical_plus
            self.plus[i] = tensor(*composition)

        if debug: print('Generating S-...')

        for i, op in enumerate(self.minus):
            composition = [canonical_id.copy() for _ in range(n_particles)]
            composition[i] = canonical_minus
            self.minus[i] = tensor(*composition)

        if debug: print('Generating Sz...')

        for i, op in enumerate(self.z):
            composition = [canonical_id.copy() for _ in range(n_particles)]
            composition[i] = canonical_z
            self.z[i] = tensor(*composition)


        self.x = (self.plus + self.minus) * 0.5
        self.y = (self.plus - self.minus) * 0.5

    def makeH(self, i, j):
        return -2 * J * (0.5 * (np.matmul(self.plus[i], self.minus[j]) + np.matmul(self.minus[i], self.plus[j])) + np.matmul(self.z[i], self.z[j]))

    def make_coupled_H(self, *pairs):
        return sum([self.makeH(*pair) for pair in pairs])

    def make_magnetization(self):
        return sum(self.z)*mu


def get_Eg(H):
    if H is 0: return (0,),(0,)

    energies_raw = np.linalg.eigvals(H).tolist()
    #print(np.linalg.eig(H))
    energies = []
    [energies.append(i) for i in energies_raw if i not in energies]
    degeneracies = [energies_raw.count(i) for i in energies]
    #print(H.shape, energies, degeneracies, sep='\n')
    return energies, degeneracies

#I = np.matrix([[-1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]])
#print(H - 1.5*I)
#print(np.linalg.eigvals(H - 1.5*I))

def Z(energies, degeneracies, T):
    total = 0
    for g, E in zip(degeneracies, energies):
        total += g*np.exp(-E/(kb*T))
    return total

def E_avg(energies, degeneracies, T):
    total = 0
    for g, E in zip(degeneracies, energies):
        total += g * E*np.exp(-E / (kb*T))
    return total / Z(energies, degeneracies, T)


