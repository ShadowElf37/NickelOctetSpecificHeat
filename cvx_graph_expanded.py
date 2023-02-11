import matplotlib.pyplot as plt
import numpy as np
import pickle


with open('cvx_data_expanded', 'rb') as f:
    data = pickle.load(f)

C2 = data[0]
C3 = data[1]
C4 = data[2:5]
C5 = data[5:9]
C6 = data[9:19]
C8 = data[19]
T = data[20]
#X1 = data[21]
X2 = data[21]
X3 = data[22]
X4 = data[23:26]
X5 = data[26:30]
X6 = data[30:40]
X8 = data[40]

print(len(data))

structures = list(reversed(['0301', '0220', '0040', '04001', '03110', '01310', '02300', '041010', '023010', '050001', '032100', '014100', '022200', '040200', '004200', '024000', '006000']))
structures *= 2

fig, ((pC2,pX2),(pC3,pX3),(pC4,pX4),(pC5,pX5),(pC6,pX6),(pC8,pX8),) = plt.subplots(6,2)
#fig.set_size_inches(12, 5)

pC = pC2, pC3, pC4, pC5, pC6, pC8
pX = pX2, pX3, pX4, pX5, pX6, pX8

#print('$', C2, A2, M2, sep='%%\n')

pC2.plot(T, C2, label="Dimer (02)", linewidth=1.5)
pC3.plot(T, C3, label="Trimer (021)", linewidth=1.5)
for c in C4:
    pC4.plot(T, c, label="Tetramer (" + structures.pop() + ")", linewidth=1.5)
for c in C5:
    pC5.plot(T, c, label="Pentamer (" + structures.pop() + ")", linewidth=1.5)
for c in C6:
    pC6.plot(T, c, label="Hexamer (" + structures.pop() + ")", linewidth=1.5)
pC8.plot(T, C8, label="Otcamer (006)", linewidth=1.5)

#plt.plot(T, C16, color="red", label="Hexadecamer", linewidth=1.5)

#plt.plot(T, A2, color="red", label="Analytic Dimer", linestyle="dotted", linewidth=2)
#plt.plot(T, G2, color="purple", label="Gregorio Dimer", linestyle="dotted", linewidth=2)

for p in pC:
    p.set_xlim([0,5])

    p.set_title('Cv per ion')

    p.set_xlabel('Temperature')
    p.set_ylabel('Specific Heat')

    #p.legend()

#plt.plot(T, XA2, label="Dimer XA", linewidth=2, linestyle="dotted")
pX2.plot(T, C2, label="Dimer (02)", linewidth=1.5)
pX3.plot(T, C3, label="Trimer (021)", linewidth=1.5)
for x in X4:
    pX4.plot(T, x, label="Tetramer (" + structures.pop() + ")", linewidth=1.5)
for x in X5:
    pX5.plot(T, x, label="Pentamer (" + structures.pop() + ")", linewidth=1.5)
for x in X6:
    pX6.plot(T, x, label="Hexamer (" + structures.pop() + ")", linewidth=1.5)
pX8.plot(T, X8, label="Otcamer (006)", linewidth=1.5)
#plt.plot(B, X2B, color="blue", label="B = 1", linewidth=1.5)
#plt.plot(B, X2C, color="green", label="B = 2", linewidth=1.5)

for p in pX:
    p.set_title('χ per ion')

    p.set_xlabel('Temperature')
    p.set_ylabel('Susceptibility')

    #p.legend()

    p.set_ylim(0, 0.5)

plt.suptitle('J=-1 Specific Heat and Susceptibility')

plt.show()