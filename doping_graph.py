import matplotlib.pyplot as plt
import numpy as np
import pickle
from collections import defaultdict
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

with open('data.npy', 'rb') as f:
    C2, C3, C4, C8, T, X1,X2,X3,X4,X8 = pickle.load(f)

with open('doping_data.npy', 'rb') as f:
    data = pickle.load(f)


def compute_C(i):
    ddict = defaultdict(int)
    ddict.update(dict(data[i][1]))
    C = (ddict[3] + ddict[4] + ddict[5] + ddict[6]) * C8 + ddict[2] * C4 + ddict[1] * C2 + ddict[0] * 0
    return np.real(C)

def compute_X(i):
    ddict = defaultdict(int)
    ddict.update(dict(data[i][1]))
    C = (ddict[3] + ddict[4] + ddict[5] + ddict[6]) * X8 + ddict[2] * X4 + ddict[1] * X2 + ddict[0] * X1
    return np.real(C)


# SPECIFIC HEAT ======================
specific_heat = []
for i in range(len(data)):
    specific_heat.append(compute_C(i))

susc = []
for i in range(len(data)):
    susc.append(compute_X(i))

fig, (ax1,ax2) = plt.subplots(1, 2)
fig.set_size_inches(12, 5)

def animate(i):
    x = round(data[i][0] / 1000000, 2)

    C = specific_heat[i]
    ax1.clear()
    ax1.set_ylim((0, 0.5))
    ax1.set_xlabel('Temperature')
    ax1.set_ylabel('Specific Heat')
    ax1.set_title(f'Cv per ion (from spin)')
    ax1.plot(T, C, color='#32d')

    X = susc[i]
    ax2.clear()
    ax2.set_ylim((0, 1))
    ax2.set_xlabel('Temperature')
    ax2.set_ylabel('Susceptibility')
    ax2.set_title(f'χ per ion')
    ax2.plot(T, X, color='#d32')

    plt.suptitle(f'J=-1, x={x}')

ani = FuncAnimation(fig, animate, frames=len(specific_heat), interval=50, repeat=True)

writervideo = animation.FFMpegWriter(fps=20)
ani.save('doped_Cv_X.mp4', writer=writervideo)

plt.show()
