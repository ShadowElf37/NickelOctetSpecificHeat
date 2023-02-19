"""
Tools to construct spin Hamiltonians and partition functions and stuff

Xm computation from https://bingweb.binghamton.edu/~suzuki/ThermoStatFIles/15.5%20%20%20QM%20Magnetization%20susceptibility.pdf
"""


import numpy as np
import scipy.linalg
from scipy.linalg import expm, fractional_matrix_power
import time
from itertools import starmap

J = -1
kb = 1
mu = -1

def tensor(*ops):
    total = ops[0]
    if len(ops) > 1:
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

        self.zero = np.zeros((n,n))
        self.identity = np.identity(n)

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

    def __makeH(self, i, j):
        return -2 * J * (0.5 * (np.matmul(self.plus[i], self.minus[j]) + np.matmul(self.minus[i], self.plus[j])) + np.matmul(self.z[i], self.z[j]))

    def make_hamiltonian(self, *pairs):
        return sum(self.__makeH(*pair) for pair in pairs) + self.zero

    def make_magnetization(self):
        return sum(self.z)*mu


def get_Eg(H):
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


def compute_C(H, T, dT):
    #print(dim(H))
    E, g = get_Eg(H)
    return np.gradient(E_avg(E, g, T), dT)
def compute_X(MM, H, T):
    #print(dim(H))
    E, g = get_Eg(H)
    z = Z(E, g, T)
    kbt = T*kb
    N = len(T)


    #t0 = time.time()
    #th = -np.tile(H, (N, 1, 1)) / kbt[:,None,None]
    #np.tile(MM, (N, 1, 1))
    #th = np.tile(expm(-H), (N, 1, 1))
    #t4 = time.time()
    #rho = np.array(list(starmap(fractional_matrix_power, zip(th, (1/kbt)))))
    #t5 = time.time()
    #chi = np.array(list(map(lambda m: np.trace(np.matmul(MM, m)), rho))) / (z * kbt)
    #chi = np.array(list(map(np.trace, np.tile(MM, (N, 1, 1)) * np.array(list(map(lambda m: expm(m) ** (1 / kbt[:,None,None]), np.tile(-H, (N,1,1)))))))) / (z * kbt)

    #he = np.linalg.eigvalsh(-H)

    #chi = np.array([(np.dot(np.diagonal(MM), np.diagonal(np.exp(-H/t))))/(z[i]*t) for i,t in enumerate(T*kb)])


    #print(MM)

    #t1 = time.time()
    chi = np.array([np.trace(np.matmul(MM, expm(-H/t)))/(z[i]*t) for i,t in enumerate(T*kb)])
    #print(t1-t0, time.time()-t1, t5-t4, sep='/')
    return chi
