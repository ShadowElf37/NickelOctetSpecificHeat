import numpy as np
import decimal
from numpy import exp
from scipy.linalg import expm

def tensor(*ops):
    total = ops[0]
    for op in ops[1:]:
        total = np.kron(total, op)
    return total

def check_symmetric(a, tol=1e-8):
    return np.all(np.abs(a-a.T) < tol)


J = -1
kb = 1
mu = -1/2

#decimal.Decimal = lambda x: np.ndarray.astype(x, np.float64)

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

print('Creating spin spaces...')
s2 = SpinSpace(2)
s3 = SpinSpace(3)
s4 = SpinSpace(4)
s8 = SpinSpace(8)
#s16 = SpinSpace(16, True)

# 2 SPIN: s.makeH(0, 1)
# 4 SPIN: s.makeH(0, 1) + s.makeH(0, 2) + s.makeH(2, 3) + s.makeH(1, 3)
# 8 SPIN: s.makeH(0, 1) + s.makeH(0, 2) + s.makeH(2, 3) + s.makeH(1, 3)
#       + s.makeH(4, 5) + s.makeH(4, 6) + s.makeH(5, 7) + s.makeH(6, 7)
#       + s.makeH(4, 0) + s.makeH(5, 2) + s.makeH(6, 3) + s.makeH(7, 4)


print('Generating Hamiltonians...')
H2 = s2.make_coupled_H((0, 1))
H3 = s3.make_coupled_H((0, 1), (0, 2))
H4 = s4.make_coupled_H((0, 1), (0, 2), (2, 3), (1, 3))
H8 = s8.makeH(0, 1) + s8.makeH(0, 2) + s8.makeH(2, 3) + s8.makeH(1, 3) + s8.makeH(4, 5) + s8.makeH(4, 6) + s8.makeH(6, 7) + s8.makeH(5, 7) + s8.makeH(4, 0) + s8.makeH(5, 1) + s8.makeH(6, 2) + s8.makeH(7, 3)
#H16 = s16.makeH(0, 1) + s16.makeH(0, 2) + s16.makeH(2, 3) + s16.makeH(1, 3) + s16.makeH(4, 5) + s16.makeH(4, 6) + s16.makeH(5, 7) + s16.makeH(6, 7) + s16.makeH(4, 0) + s16.makeH(5, 1) + s16.makeH(6, 2) + s16.makeH(7, 3)\
#    + s16.makeH(8, 9) + s16.makeH(8, 10) + s16.makeH(10, 11) + s16.makeH(9, 11) + s16.makeH(12, 13) + s16.makeH(12, 14) + s16.makeH(13, 15) + s16.makeH(14, 15) + s16.makeH(12, 0) + s16.makeH(13, 9) + s16.makeH(14, 10) + s16.makeH(15, 11)\
 #   + s16.makeH(0, 8) + s16.makeH(1, 9) + s16.makeH(2, 10) + s16.makeH(3, 11) + s16.makeH(4, 12) + s16.makeH(5, 13) + s16.makeH(6, 14) + s16.makeH(7, 15)

M2 = s2.make_magnetization()
M3 = s3.make_magnetization()
M4 = s4.make_magnetization()
M8 = s8.make_magnetization()

print('Hermitian check:')
for H in (H2, H3, H4, H8):
    if not check_symmetric(H):
        print('H%d isn\'t Hermitian!' % np.log2(H.shape[0]))
        raise StopIteration
print('All passed.')


"""print(s2.states)
print(np.matmul(s2.plus[0], s2.states[0]))
print(np.matmul(s2.plus[0], s2.states[1]))
print(np.matmul(s2.plus[0], s2.states[2]))
print(np.matmul(s2.plus[0], s2.states[3]))
print(np.matmul(s2.plus[1], s2.states[0]))
print(np.matmul(s2.plus[1], s2.states[1]))
print(np.matmul(s2.plus[1], s2.states[2]))
print(np.matmul(s2.plus[1], s2.states[3]))"""

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
        try:
            total += g * E*np.exp(-E / (kb*T))
        except RuntimeWarning as e:
            print(e)

    return total / Z(energies, degeneracies, T)



dx = 0.01
T = np.arange(dx, 10, dx, dtype=np.float64)
#B = np.arange(dx, 5, dx, dtype=np.float64)


def compute_C(H):
    E, g = get_Eg(H)
    return np.gradient(E_avg(E, g, T), dx)


def compute_chi(M, H):
    E, g = get_Eg(H)
    z = Z(E, g, T)
    chi = np.zeros_like(T)
    MM = M * M
    for i,t in enumerate(T*kb):
        rho = expm(-H/t)/z[i]
        chi[i] = np.trace(MM * rho)/t
    return chi

#print(H4)
#e = list(map(np.round, np.linalg.eigvals(H4)))
#s = list(set(e))
#print(s, [e.count(i) for i in s])

print('Calculating H2...')
C2 = compute_C(H2)/2
print('Calculating H3...')
C3 = compute_C(H3)/3
print('Calculating H4...')
C4 = compute_C(H4)/4
#C4 = np.zeros_like(T, dtype=np.float64)
print('Calculating H8...')
C8 = compute_C(H8)/8
#C8 = np.zeros_like(T, dtype=np.float64)
print('Calculating X2...')
X2 = compute_chi(M2, H2)
print('Calculating X3...')
X3 = compute_chi(M3, H3)
print('Calculating X4...')
X4 = compute_chi(M4, H4)
print('Calculating X8...')
X8 = compute_chi(M8, H8)#np.zeros_like(X4)#


#print(H2)
#print(np.matrix([[-0.5, 0, 0, 0], [0, 0.5, -1, 0], [0, -1, 0.5, 0], [0, 0, 0, -0.5]]))
#M2 = compute_C(np.matrix([[-0.5, 0, 0, 0], [0, 0.5, -1, 0], [0, -1, 0.5, 0], [0, 0, 0, -0.5]]))

kbt = T*kb

G2 =  0.5 * 4/3 * J**2 * np.exp(2*J/(kbt)) / (kb* T**2 * (1/3 + np.exp(2*J/(kbt)))**2)
#THIS IS WRONG HAMILTONIAN: A2 = 2*J**2*exp(2*J/(kbt))*(25*exp(3*J/(kbt)) + 9*exp(5*J/(kbt)) + 2)/(T**2*kb*(exp(2*J/(kbt)) + 2*exp(5*J/(kbt)) + 1)**2)

XA2 = 8/(kbt * (3 + np.exp(J / kbt)))

print("Saving to file.")

with open('data.npy', 'wb') as file:
    np.save(file, (C2, C3, C4, C8, G2, T), allow_pickle=True)
    np.save(file, (X2,X3,X4,X8, XA2), allow_pickle=True)

print('Done.')

import graph