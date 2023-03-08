"""
Animated graph of the doped sample
"""

import matplotlib.pyplot as plt
import numpy as np
import pickle
from collections import defaultdict
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import chainlib

N = chainlib.Chain.N_from_ss

Tc, Cdict = pickle.load(open('sample_cv_data_2', 'rb'))
Tx, Xdict = pickle.load(open('sample_xm_data_2', 'rb'))


# UNCOMMENT BOTTOM TO USE UNBROKEN CHAINS
#_, counts = pickle.load(open('final_structures', 'rb'))

counts = [defaultdict(int) for _ in range(101)]
for i in range(9):
    with open(f'real_structures{i}', 'rb') as f:
        #data = [pair[0] for pair in pickle.load(f)]
        data = pickle.load(f)
        for i in range(101):
            for key in data[i].keys():
                counts[i][key] += data[i][key]

#counts =  pickle.load(open('real_structures', 'rb'))#[data[0] for data in]
totals = [sum(v * (N(k) or 1) for (k,v) in x.items()) for x in counts]

C = [sum(Cdict.get(s, Cdict['0.0.0.8.0.0.0']*N(s)/8) * count * (N(s) or 1) / totals[i] for s,count in x.items()) for i,x in enumerate(counts)]
X = [sum(Xdict.get(s, Xdict['0.0.0.8.0.0.0']*N(s)/8) * count * (N(s) or 1)/ totals[i] for s,count in x.items()) for i,x in enumerate(counts)]

DOPE_X_MAX = len(C)

fig, (ax1,ax2) = plt.subplots(1, 2)
fig.set_size_inches(12, 5)

def animate(i):
    x = round(i / DOPE_X_MAX, 2)

    ax1.clear()
    ax1.set_ylim((0, 1.2))
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

print('Animating...')
ani = FuncAnimation(fig, animate, frames=DOPE_X_MAX, interval=50, repeat=True)

#plt.show()

print('Saving...')
writervideo = animation.FFMpegWriter(fps=30)
ani.save('doped_CvXm.mp4', writer=writervideo, dpi=150)
