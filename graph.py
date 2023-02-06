import matplotlib.pyplot as plt
import numpy as np


SPECIFICHEAT = False
SUSCEPTIBILITY = True


with open('data.npy', 'rb') as f:
    C2, C3, C4, C8, G2, T = np.load(f, allow_pickle=True)
    X2,X3,X4,X8,XA2 = np.load(f, allow_pickle=True)



#print('$', C2, A2, M2, sep='%%\n')
if SPECIFICHEAT:
    plt.plot(T, C2, color="lightblue", label="Dimer", linewidth=2)
    plt.plot(T, C3, color="blue", label="Trimer", linewidth=1.5)
    plt.plot(T, C4, color="green", label="Tetramer", linewidth=1.5)
    plt.plot(T, C8, color="orange", label="Octamer", linewidth=1.5)
    #plt.plot(T, C16, color="red", label="Hexadecamer", linewidth=1.5)

    #plt.plot(T, A2, color="red", label="Analytic Dimer", linestyle="dotted", linewidth=2)
    plt.plot(T, G2, color="purple", label="Gregorio Dimer", linestyle="dotted", linewidth=2)


    plt.legend()

    plt.title('J=-1 Specific Heat per Ion')

    plt.xlabel('Temperature')
    plt.ylabel('Specific Heat')

    plt.show()

if SUSCEPTIBILITY:
    plt.plot(T, XA2, label="Dimer XA", linewidth=2, linestyle="dotted")
    plt.plot(T, X2, label="Dimer X", linewidth=2)
    plt.plot(T, X3, label="Trimer X", linewidth=1.5)
    plt.plot(T, X4, label="Tetramer X", linewidth=1.5)
    plt.plot(T, X8, label="Octamer X", linewidth=1.5)
    #plt.plot(B, X2B, color="blue", label="B = 1", linewidth=1.5)
    #plt.plot(B, X2C, color="green", label="B = 2", linewidth=1.5)

    plt.title('J=-1 Magnetic Susceptibility')

    plt.xlabel('Temperature')
    plt.ylabel('Susceptibility')

    plt.legend()

    plt.ylim(0, 10)

    plt.show()