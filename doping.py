import numpy as np

length = 3
bonds = np.ones((length, length, length, 3), np.int8)


def boolv(x,y,z):
    return np.array((x,y,z), dtype=bool)

def dope(bonds, x,y,z):
    X = np.array((1, 0, 0))
    Y = np.array((0, 1, 0))
    Z = np.array((0, 0, 1))

    center = np.array((x,y,z))

    for point, val in np.ndenumerate(bonds[x:x+2,y:y+2,z:z+2]):
        #print(point)
        bonds[point] = np.multiply(val, (point[:-1]-center)[point[-1]])

    #bonds[point] = np.multiply(bonds[point], (0,0,0))
    #bonds[point+X] = np.multiply(bonds[point+X], (1,0,0))
    #bonds[point+Y] = np.multiply(bonds[point+Y], (0,1,0))
    #bonds[point+Z] = np.multiply(bonds[point+Z], (0,0,1))
    #bonds[point+X+Z] = np.multiply(bonds[point+X+Z], (1,0,1))
    #bonds[point+Y+Z] = np.multiply(bonds[point+Y+Z], (0,1,1))
    #bonds[point+X+Y] = np.multiply(bonds[point+X+Y], (1,1,0))


dope(bonds, 0,0,0)

print(bonds)