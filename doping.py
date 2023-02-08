import numpy as np
from random import randint


def make_bonds_array(dimx, dimy, dimz):
    return np.ones((dimx, dimy, dimz, 3), np.int8)

def dope(arr, x,y,z):
    doper = np.array([[[[0., 0., 0.],
                        [0., 0., 1.]],
                       [[0., 1., 0.],
                        [0., 1., 1.]]],
                      [[[1., 0., 0.],
                        [1., 0., 1.]],
                       [[1., 1., 0.],
                        [1., 1., 1.]]]], dtype=np.int8)
    arr[x:x + 2, y:y + 2, z:z + 2] *= doper


def to_vertices(edges):
    disp_picker = np.identity(3, dtype=np.int8)

    vertices = np.zeros(edges.shape[:3], dtype=np.int8)
    for p, val in np.ndenumerate(edges):
        coord = p[:3]
        #print(p, val, coord, disp_picker[p[3]])
        vertices[coord] += val
        try:
            vertices[tuple(coord+disp_picker[p[3]])] += val
        except IndexError:
            pass

    return vertices


def random_dope(arr, count):
    for _ in range(count):
        dope(arr, randint(1, arr.shape[0]-2),randint(1, arr.shape[1]-2),randint(1, arr.shape[2]-2))


def produce_sample(length, doping_count):
    length = length+1
    bonds = make_bonds_array(length, length, length)
    random_dope(bonds, doping_count)
    return to_vertices(bonds)[1:,1:,1:]

def count_polymers(vertex_arr):
    return list(zip(*np.unique(vertex_arr, return_counts=True)))

def count_polymer_pct(vertex_arr):
    total = np.prod(vertex_arr.shape)
    counts = count_polymers(vertex_arr)
    for i in range(len(counts)):
        counts[i] = (counts[i][0], counts[i][1]/total)
    return counts


sample: np.ndarray = produce_sample(100, 1000000//10)
print(count_polymer_pct(sample))