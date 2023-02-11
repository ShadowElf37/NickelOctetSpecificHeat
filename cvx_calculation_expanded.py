import numpy as np
import pickle
from scipy.linalg import expm
from spinlib import *
import nperrlog


print('Creating spin spaces...')
s2 = SpinSpace(2)
s3 = SpinSpace(3)
s4 = SpinSpace(4)
s5 = SpinSpace(5)
s6 = SpinSpace(6)
s8 = SpinSpace(8)


print('Generating Hamiltonians...')
H = (
        #s2.make_coupled_H((1, 0), ),
        #s3.make_coupled_H((1, 0), (0, 2)),
        s4.make_coupled_H((0, 2), (1, 0), (0, 3)),
        s4.make_coupled_H((1, 0), (1, 2), (0, 3)),
        s4.make_coupled_H((0, 2), (1, 3), (1, 0), (3, 2)),
        s5.make_coupled_H((2, 0), (1, 2), (2, 4), (2, 3)),
        s5.make_coupled_H((2, 0), (1, 2), (1, 3), (2, 4)),
        s5.make_coupled_H((1, 0), (3, 0), (2, 1), (2, 3), (3, 4)),
        s5.make_coupled_H((1, 4), (0, 1), (3, 2), (0, 3)),
        s6.make_coupled_H((3, 0), (2, 3), (0, 1), (3, 4), (3, 5)),
        s6.make_coupled_H((3, 0), (1, 4), (2, 3), (0, 1), (3, 4), (3, 5)),
        s6.make_coupled_H((3, 1), (3, 0), (2, 3), (3, 4), (3, 5)),
        s6.make_coupled_H((3, 0), (2, 3), (0, 1), (3, 5), (2, 4)),
        s6.make_coupled_H((3, 0), (1, 5), (2, 3), (0, 1), (3, 5), (2, 4)),
        s6.make_coupled_H((3, 0), (2, 3), (3, 5), (1, 2), (0, 1), (2, 4)),
        s6.make_coupled_H((3, 0), (2, 1), (2, 3), (3, 5), (2, 4)),
        s6.make_coupled_H((1, 0), (4, 0), (3, 1), (3, 4), (0, 2), (2, 5), (4, 5)),
        s6.make_coupled_H((4, 1), (3, 4), (0, 2), (2, 5), (0, 3)),
        s6.make_coupled_H((1, 0), (2, 5), (3, 1), (3, 4), (0, 2), (4, 5)),
        #s8.make_coupled_H((0, 1), (0, 2), (2, 3), (1, 3), (4, 5), (4, 6), (6, 7), (5, 7), (4, 0), (5, 1), (6, 2), (7, 3)),
)

print('Generating magnetizations...')
M = (
    #s2.make_magnetization(),
    #s3.make_magnetization(),
    *(s4.make_magnetization() for _ in range(3)),
    *(s5.make_magnetization() for _ in range(4)),
    *(s6.make_magnetization() for _ in range(10)),
    #s8.make_magnetization(),
)

print('Hermitian check...')
for h in H:
    if not check_symmetric(h):
        print('H%d isn\'t Hermitian!' % dim(H)/2)
        raise StopIteration
print('All passed.')


dt = 0.01
T = np.arange(dt, 10, dt, dtype=np.float64)

def compute_C(H):
    print(dim(H))
    E, g = get_Eg(H)
    return np.gradient(E_avg(E, g, T), dt)
def compute_X(M, H):
    print(dim(H))
    E, g = get_Eg(H)
    z = Z(E, g, T)
    MM = M * M
    chi = np.array([np.trace(MM * expm(-H/t)/z[i])/t for i,t in enumerate(T*kb)])
    return chi

print('Calculating C...')
C = [compute_C(H[i])/dim(H[i]) for i in range(len(H))]

print('Calculating X...')
X = [compute_X(M[i], H[i])/dim(H[i]) for i in range(len(H))]#[mu**2/(kb*T)] +


print("Saving to file...")

with open('cvx_data_expanded', 'wb') as file:
    pickle.dump((*C, T, *X), file)

print('Done.')


print(nperrlog.WARNS, 'warnings were ignored.')




#print(H2)
#print(np.matrix([[-0.5, 0, 0, 0], [0, 0.5, -1, 0], [0, -1, 0.5, 0], [0, 0, 0, -0.5]]))
#M2 = compute_C(np.matrix([[-0.5, 0, 0, 0], [0, 0.5, -1, 0], [0, -1, 0.5, 0], [0, 0, 0, -0.5]]))

#kbt = T*kb

#G2 =  0.5 * 4/3 * J**2 * np.exp(2*J/(kbt)) / (kb* T**2 * (1/3 + np.exp(2*J/(kbt)))**2)
#THIS IS WRONG HAMILTONIAN: A2 = 2*J**2*exp(2*J/(kbt))*(25*exp(3*J/(kbt)) + 9*exp(5*J/(kbt)) + 2)/(T**2*kb*(exp(2*J/(kbt)) + 2*exp(5*J/(kbt)) + 1)**2)

#XA2 = 8/(kbt * (3 + np.exp(J / kbt)))

#import graph