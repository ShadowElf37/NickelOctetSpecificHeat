import matplotlib.pyplot as plt
import numpy as np
import pickle


with open('cvx_data_expanded', 'rb') as f:
    data = pickle.load(f)

C4 = data[0:3]
C5 = data[3:7]
C6 = data[7:17]
T = data[17]
#X1 = data[21]
X4 = data[18:21]
X5 = data[21:25]
X6 = data[25:35]

print(len(data))

structures = list(reversed(['0301', '0220', '0040', '04001', '03110', '01310', '02300', '041010', '023010', '050001', '032100', '014100', '022200', '040200', '004200', '024000', '006000']))
structures *= 2

fig, ((pC4,pX4),(pC5,pX5),(pC6,pX6),) = plt.subplots(3,2)
#fig.set_size_inches(12, 5)

pC =  pC4, pC5, pC6,
pX =  pX4, pX5, pX6,


for c in C4:
    pC4.plot(T, c, label="Tetramer (" + structures.pop() + ")", linewidth=1.5)
for c in C5:
    pC5.plot(T, c, label="Pentamer (" + structures.pop() + ")", linewidth=1.5)
for c in C6:
    pC6.plot(T, c, label="Hexamer (" + structures.pop() + ")", linewidth=1.5)

for p in pC:
    p.set_xlim([0,5])
    p.set_ylabel('Specific Heat')
    p.legend(loc='upper right')

pC4.set_title('Cv per ion')
pC6.set_xlabel('Temperature')

for x in X4:
    pX4.plot(T, x, label="Tetramer (" + structures.pop() + ")", linewidth=1.5)
for x in X5:
    pX5.plot(T, x, label="Pentamer (" + structures.pop() + ")", linewidth=1.5)
for x in X6:
    pX6.plot(T, x, label="Hexamer (" + structures.pop() + ")", linewidth=1.5)

for p in pX:
    p.set_ylabel('Susceptibility')
    p.legend(loc='upper right')
    p.set_ylim(0, 0.5)

pX4.set_title('χ per ion')
pX6.set_xlabel('Temperature')

plt.suptitle('J=-1 Specific Heat and Susceptibility')

plt.show()