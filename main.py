import numpy as np
import decimal
from numpy import exp

def tensor(*ops):
    total = ops[0]
    for op in ops[1:]:
        total = np.kron(total, op)
    return total



class SpinSpace:
    def __init__(self, n_particles):
        # THIS IS SPIN 1/2!!!

        n = 2**n_particles
        self.n = n_particles

        self.zero = np.zeros(n)
        self.states = tuple(np.zeros(n) for _ in range(n))
        for i, state in enumerate(self.states):
            state[i] = 1

        self.plus = np.array([np.identity(n) for _ in range(n_particles)])
        self.minus = np.array([np.identity(n) for _ in range(n_particles)])
        self.z = np.array([np.identity(n) for _ in range(n_particles)])

        canonical_id = np.identity(2)
        canonical_plus = np.array([[0,1],[0,0]])
        canonical_minus = np.array([[0,0],[1,0]])
        canonical_z = np.array([[1,0],[0,-1]])*0.5

        for i, op in enumerate(self.plus):
            composition = [canonical_id.copy() for _ in range(n_particles)]
            composition[i] = canonical_plus
            self.plus[i] = tensor(*composition)

        for i, op in enumerate(self.minus):
            composition = [canonical_id.copy() for _ in range(n_particles)]
            composition[i] = canonical_minus
            self.minus[i] = tensor(*composition)

        for i, op in enumerate(self.z):
            composition = [canonical_id.copy() for _ in range(n_particles)]
            composition[i] = canonical_z
            self.z[i] = tensor(*composition)


    def state(self, *spins):
        return self.states[int(''.join(map(str, spins)), 2)]

    def makeH(self, i, j):
        return -2 * J * (0.5 * (np.matmul(self.plus[i], self.minus[j]) + np.matmul(self.minus[i], self.plus[j])) + np.matmul(self.z[i], self.z[j]))


s2 = SpinSpace(2)
s4 = SpinSpace(4)
s8 = SpinSpace(8)

J = 1

# 2 SPIN: s.makeH(0, 1)
# 4 SPIN: s.makeH(0, 1) + s.makeH(0, 2) + s.makeH(2, 3) + s.makeH(1, 3)
# 8 SPIN: s.makeH(0, 1) + s.makeH(0, 2) + s.makeH(2, 3) + s.makeH(1, 3)
#       + s.makeH(4, 5) + s.makeH(4, 6) + s.makeH(5, 7) + s.makeH(6, 7)
#       + s.makeH(4, 0) + s.makeH(5, 2) + s.makeH(6, 3) + s.makeH(7, 4)

H2 = s2.makeH(0, 1)
print('H=', H2)
H4 = (s4.makeH(0, 1) + s4.makeH(0, 2) + s4.makeH(2, 3) + s4.makeH(1, 3))
H8 = s8.makeH(0, 1) + s8.makeH(0, 2) + s8.makeH(2, 3) + s8.makeH(1, 3) + s8.makeH(4, 5) + s8.makeH(4, 6) + s8.makeH(5, 7) + s8.makeH(6, 7) + s8.makeH(4, 0) + s8.makeH(5, 2) + s8.makeH(6, 3) + s8.makeH(7, 4)

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
    return energies, degeneracies

#I = np.matrix([[-1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]])
#print(H - 1.5*I)
#print(np.linalg.eigvals(H - 1.5*I))

kb = 1

print(decimal.getcontext())

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
        except decimal.Overflow as e:
            print(e)
            #print(-E / (kb*T))
            raise KeyboardInterrupt
    return total / Z(energies, degeneracies, T)



dx = 0.01
T = np.arange(dx, 5, dx, dtype=np.float64)

decimal.getcontext().prec=200
conv = lambda arr: np.asarray([decimal.Decimal(t) for t in arr], dtype=object)
T = conv(T)

def compute_C(H, n=1):
    E, g = get_Eg(H)
    print(E)
    print(g)
    #try:
    #    print(max(E))
    #    print(min(E))
    #except: pass

    return np.gradient(np.array(list(map(float, E_avg(conv(E)/decimal.Decimal(n), conv(g), T)))), dx)

print(H4)
e = list(map(round, np.linalg.eigvals(H4)))
s = list(set(e))
print(s, [e.count(i) for i in s])

print('CALCULATING H2...')
C2 = compute_C(H2,1)
print('CALCULATING H4...')
C4 = compute_C(H4,2)
#C4 = np.zeros_like(T, dtype=np.float64)
print('CALCULATING H8...')
C8 = compute_C(H8,4)
#C8 = np.zeros_like(T, dtype=np.float64)



#print(H2)
#print(np.matrix([[-0.5, 0, 0, 0], [0, 0.5, -1, 0], [0, -1, 0.5, 0], [0, 0, 0, -0.5]]))
#M2 = compute_C(np.matrix([[-0.5, 0, 0, 0], [0, 0.5, -1, 0], [0, -1, 0.5, 0], [0, 0, 0, -0.5]]))

T = np.array(list(map(float, T)))

kbt = T*kb

G2 = 4/3 * J**2 * np.exp(2*J/(kb*T)) / (kb* T**2 * (1/3 + np.exp(2*J/(kb*T)))**2)
#THIS IS WRONG HAMILTONIAN: A2 = 2*J**2*exp(2*J/(kbt))*(25*exp(3*J/(kbt)) + 9*exp(5*J/(kbt)) + 2)/(T**2*kb*(exp(2*J/(kbt)) + 2*exp(5*J/(kbt)) + 1)**2)


print("Saving to file.")

with open('data.npy', 'wb') as file:
    np.save(file, C2)
    np.save(file, C4)
    np.save(file, C8)
    np.save(file, G2)
    np.save(file, T)

print('Done.')

import graph