"""
Animated graph of the doped sample
"""

import matplotlib.pyplot as plt
import numpy as np
import pickle
from collections import defaultdict
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

Tc, Tx, Cdict, Xdict = pickle.load(open('sample_cvx_data_1', 'rb'))
_, counts = pickle.load(open('final_structures', 'rb'))
totals = [sum(x.values()) for x in counts]

C = [sum(Cdict[s] * count / totals[i] for s,count in x.items()) for i,x in enumerate(counts)]
X = [sum(Xdict[s] * count / totals[i] for s,count in x.items()) for i,x in enumerate(counts)]

DOPE_X_MAX = len(C)

fig, (ax1,ax2) = plt.subplots(1, 2)
fig.set_size_inches(12, 5)

def animate(i):
    x = round(i / DOPE_X_MAX, 2)

    ax1.clear()
    ax1.set_ylim((0, 0.5))
    ax1.set_xlabel('Temperature')
    ax1.set_ylabel('Specific Heat')
    ax1.set_title(f'Cv per ion (from spin)')
    ax1.plot(Tc, C[i], color='#32d')

    ax2.clear()
    ax2.set_ylim((0, 1))
    ax2.set_xlabel('Temperature')
    ax2.set_ylabel('Susceptibility')
    ax2.set_title(f'χ per ion')
    ax2.plot(Tx, X[i], color='#d32')

    plt.suptitle(f'J=-1, x={x}')

ani = FuncAnimation(fig, animate, frames=DOPE_X_MAX, interval=50, repeat=True)

writervideo = animation.FFMpegWriter(fps=30)
ani.save('doped_Cv_X.mp4', writer=writervideo, dpi=150)
