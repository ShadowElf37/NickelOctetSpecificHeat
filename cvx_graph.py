"""
Classic Cv and Xm graph
"""

import matplotlib.pyplot as plt
import numpy as np
import pickle


with open('cvx_data', 'rb') as f:
    C2, C3, C4, C8, T, X1,X2,X3,X4,X8 = pickle.load(f)


fig, (ax1,ax2) = plt.subplots(1,2)
fig.set_size_inches(12, 5)

#print('$', C2, A2, M2, sep='%%\n')

ax1.plot(T, C2, label="Dimer", linewidth=1.5)
ax1.plot(T, C3, label="Trimer", linewidth=1.5)
ax1.plot(T, C4, label="Tetramer", linewidth=1.5)
ax1.plot(T, C8, label="Octamer", linewidth=1.5)
#plt.plot(T, C16, color="red", label="Hexadecamer", linewidth=1.5)

#plt.plot(T, A2, color="red", label="Analytic Dimer", linestyle="dotted", linewidth=2)
#plt.plot(T, G2, color="purple", label="Gregorio Dimer", linestyle="dotted", linewidth=2)

ax1.set_xlim([0,5])

ax1.set_title('Cv per ion')

ax1.set_xlabel('Temperature')
ax1.set_ylabel('Specific Heat')

ax1.legend()

#plt.plot(T, XA2, label="Dimer XA", linewidth=2, linestyle="dotted")
ax2.plot(T, X1, color='purple', label="Monomer", linewidth=1.5)
ax2.plot(T, X2, label="Dimer", linewidth=1.5)
ax2.plot(T, X3, label="Trimer", linewidth=1.5)
ax2.plot(T, X4, label="Tetramer", linewidth=1.5)
ax2.plot(T, X8, label="Octamer", linewidth=1.5)
#plt.plot(B, X2B, color="blue", label="B = 1", linewidth=1.5)
#plt.plot(B, X2C, color="green", label="B = 2", linewidth=1.5)

ax2.set_title('χ per ion')

ax2.set_xlabel('Temperature')
ax2.set_ylabel('Susceptibility')

ax2.legend()

plt.ylim(0, 1)

plt.suptitle('J=-1 Specific Heat and Susceptibility')

plt.show()