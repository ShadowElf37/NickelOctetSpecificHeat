import matplotlib.pyplot as plt
import numpy as np


with open('data.npy', 'rb') as f:
    C2 = (np.load(f))
    C4 = (np.load(f))
    C8 = np.load(f)
    A2 = np.load(f)
    T = np.load(f)



#print('$', C2, A2, M2, sep='%%\n')

plt.plot(T, C2, color="lightblue", label="Dimer", linewidth=2)
plt.plot(T, C4, color="green", label="Quadramer", linewidth=1.5)
plt.plot(T, C8, color="orange", label="Octamer", linewidth=1.5)

#plt.plot(T, A2, color="red", label="Analytic Dimer", linestyle="dotted", linewidth=2)
plt.plot(T, A2, color="purple", label="Gregorio Dimer", linestyle="dotted", linewidth=2)


plt.legend()

plt.xlabel('Temperature')
plt.ylabel('Specific Heat')



plt.show()