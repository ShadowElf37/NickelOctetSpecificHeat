import numpy as np


# TODO: Make chain kron() for higher N
# TODO: fix factor of 4 in diagonal

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
        self.id = np.identity(n_particles)

        self.states = tuple(np.zeros(n) for _ in range(n))
        for i, state in enumerate(self.states):
            state[i] = 1

        self.plus = np.array([np.identity(n) for _ in range(n_particles)])
        self.minus = np.array([np.identity(n) for _ in range(n_particles)])
        self.z = np.array([np.identity(n) for _ in range(n_particles)])

        canonical_id = np.identity(2)
        canonical_plus = np.array([[0,1],[0,0]])
        canonical_minus = np.array([[0,0],[1,0]])
        canonical_z = np.array([[1,0],[0,-1]])

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
H4 = s4.makeH(0, 1) + s4.makeH(0, 2) + s4.makeH(2, 3) + s4.makeH(1, 3)
H8 = s8.makeH(0, 1) + s8.makeH(0, 2) + s8.makeH(2, 3) + s8.makeH(1, 3) + s8.makeH(4, 5) + s8.makeH(4, 6) + s8.makeH(5, 7) + s8.makeH(6, 7) + s8.makeH(4, 0) + s8.makeH(5, 2) + s8.makeH(6, 3) + s8.makeH(7, 4)



def get_Eg(H):
    energies_raw = np.linalg.eigvals(H).tolist()
    energies = []
    [energies.append(i) for i in energies_raw if i not in energies]
    degeneracies = [energies_raw.count(i) for i in energies]
    return energies, degeneracies

#I = np.matrix([[-1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]])
#print(H - 1.5*I)
#print(np.linalg.eigvals(H - 1.5*I))

kb = 1.3*10**-23

def Z(energies, degeneracies, T):
    total = 0
    for g, E in zip(degeneracies, energies):
        total += g*np.exp(-E/T)
    return total

def E_avg(energies, degeneracies, T):
    total = 0
    for g, E in zip(degeneracies, energies):
        try:
            total += g * E*np.exp(-E / T)
        except RuntimeWarning as e:
            print(e)
    return total / Z(energies, degeneracies, T)



dx = 0.0001
T = np.linspace(dx, 10, int(1/dx), dtype=np.float64)

def compute_C(H):
    E, g = get_Eg(H)
    print(E)
    print(g)
    try:
        print(max(E))
        print(min(E))
    except: pass

    return np.gradient(E_avg(E, g, T), dx)


C2 = compute_C(H2)
C4 = compute_C(H4)
C8 = compute_C(H8)



A2 = 4/3 * J**2 * np.exp(-2*J/T) / (T**2 * (1/3 + np.exp(-2*J/T))**2)


print(H2)
print(H2 - np.matrix([[-1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,-1]]))
M2 = compute_C(H2 - np.matrix([[-1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,-1]]))

import matplotlib.pyplot as plt

plt.plot(T, C2, color="blue", label="Dimer")
plt.plot(T, C4, color="green", label="Quartet")
plt.plot(T, C8, color="orange", label="Octet")

plt.plot(T, A2, color="lightcoral", label="Gregorio Dimer (analytic)")
plt.plot(T, M2, color="lightblue", label="Gregorio Dimer (computed)")

plt.legend()

plt.xlabel('Temperature')
plt.ylabel('Specific Heat')



plt.show()