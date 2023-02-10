import numpy as np
from random import randint
import random


def make_bonds_array(dimx, dimy, dimz):
    bonds = np.ones((dimx, dimy, dimz, 6), np.int8)
    bonds[0,:,:] *= np.full((dimy, dimz, 6), [1, 0, 1, 1, 1, 1])
    bonds[-1, :, :] *= np.full((dimy, dimz, 6), [0, 1, 1, 1, 1, 1])
    bonds[:, 0, :] *= np.full((dimx, dimz, 6), [1, 1, 1, 0, 1, 1])
    bonds[:, -1, :] *= np.full((dimx, dimz, 6), [1, 1, 0, 1, 1, 1])
    bonds[:, :, 0] *= np.full((dimx, dimy, 6), [1, 1, 1, 1, 1, 0])
    bonds[:, :, -1] *= np.full((dimx, dimy, 6), [1, 1, 1, 1, 0, 1])
    return bonds

def dope(arr, x,y,z):
    doper = np.array([[[[0, 1, 0, 1, 0, 1],
                        [0, 1, 0, 1, 1, 0]],
                       [[0, 1, 1, 0, 0, 1],
                        [0, 1, 1, 0, 1, 0]]],
                      [[[1, 0, 0, 1, 0, 1],
                        [1, 0, 0, 1, 1, 0]],
                       [[1, 0, 1, 0, 0, 1],
                        [1, 0, 1, 0, 1, 0]]]], dtype=np.int8)
    dsm = arr[x:x + 2, y:y + 2, z:z + 2].shape # "doper shape mandate" - required to chop this matrix off at the corners
    arr[x:x + 2, y:y + 2, z:z + 2] *= doper[:dsm[0],:dsm[1],:dsm[2]]


def to_vertices(edges):
    vertices = np.zeros(edges.shape[:3], dtype=np.int8)
    for p, val in np.ndenumerate(edges):
        vertices[p[:3]] += val
    return vertices


def random_dope(arr, count, length, hit_unique=True):
    if hit_unique:
        to_hit = [(i, j, k) for i in range(length - 1) for j in range(length - 1) for k in range(length - 1)]
        random.shuffle(to_hit)
        for _ in range(count):
            x, y, z = to_hit.pop()
            dope(arr, x,y,z)
    else:
        for _ in range(count):
            x, y, z = randint(0, length-1), randint(0, length-1), randint(0, length-1)
            dope(arr, x,y,z)


# You probably want these 3 functions
def produce_sample(length, doping_count):
    length = length
    bonds = make_bonds_array(length, length, length)
    random_dope(bonds, doping_count, length, hit_unique=True)
    return to_vertices(bonds)

def count_polymers(vertex_arr):
    return list(zip(*np.unique(vertex_arr, return_counts=True)))

def count_polymer_pct(vertex_arr):
    total = np.prod(vertex_arr.shape)
    counts = count_polymers(vertex_arr)
    for i in range(len(counts)):
        counts[i] = (counts[i][0], counts[i][1]/total)
    return counts


if __name__=="__main__":
    import pickle as p
    print(*p.load(open('doping_data', 'rb')), sep='\n')