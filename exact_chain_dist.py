import dopelib
import chainlib
import numpy as np
import time

L = 5

#s1 = dopelib.make_bonds_array(L, L, L)
s2 = dopelib.make_bonds_array(L, L, L)

#to_hit1 = [(i, j, k) for i in range(L - 1) for j in range(L - 1) for k in range(L - 1)]
to_hit2 = [(i, j, k) for i in range(L - 1) for j in range(L - 1) for k in range(L - 1)]
for point in to_hit2.copy():
    if sum(point[j] == 0 or point[j] == 3 for j in (0,1,2)) > 1:
        to_hit2.remove(point)




def reverse_fac(n, i):
    return np.prod(range(n, n-i,-1))
class ReverseFactorialMachine:
    COMPUTED = {}
    @staticmethod
    def register(n, i):
        p = np.prod(range(n, n-i,-1))
        ReverseFactorialMachine.COMPUTED[(n,i)] = p
        return p
    @staticmethod
    def get(n, i):
        return ReverseFactorialMachine.COMPUTED.get((n,i), None) or ReverseFactorialMachine.register(n, i)



def new_without(list, item):
    new = list.copy()
    new.remove(item)
    return new
def check_isolated(sample, p):
    return not any(sample[p])
def isolated_checker(sample):
    """returns a wrapped function that can then be mapped onto points"""
    def check(p):
        return check_isolated(sample, p)
    return check

def walk_points(sample, points, points_to_check_all=((1,1,1),), points_to_check_any=(), n=0):
    isolated = 0
    #no_dimer_samples
    if n == 0:
        return 0

    for p in points:
        s = sample.copy()
        dopelib.dope(s, *p)
        if n > 1:
            new_points = new_without(points, p)
            result = walk_points(s, new_points, points_to_check_all, points_to_check_any, n=n-1)
            #if type(result) is int:
            isolated += result
            #else:

        elif all(map(isolated_checker(s), points_to_check_all)) and any(map(isolated_checker(s), points_to_check_any)):
            isolated += 1
    return isolated# or sample

def walk_points2(sample, points, points_to_check_all=((1,1,1),), points_to_check_any=(), n=0, N=32, already_doped=0):
    isolated = 0
    #no_dimer_samples

    if all(map(isolated_checker(sample), points_to_check_all)) and any(map(isolated_checker(sample), points_to_check_any)):
        isolated += reverse_fac(N-already_doped, n)
    elif n > 0:
        for p in points:
            s = sample.copy()
            dopelib.dope(s, *p)
            new_points = new_without(points, p)
            #if type(result) is int:
            isolated += walk_points2(s, new_points, points_to_check_all, points_to_check_any, n=n-1, already_doped=already_doped+1)

    return isolated# or sample

def walk_points3(sample, points, isolated_arr, points_to_check_all=((1,1,1),), points_to_check_any=(), n=0, N=32, already_doped=0):
    if all(map(isolated_checker(sample), points_to_check_all)) and any(map(isolated_checker(sample), points_to_check_any)):
        for i in range(n+1):
            isolated_arr[already_doped+i] += reverse_fac(N-already_doped, i)

    elif n > 0:
        for p in points:
            s = sample.copy()
            dopelib.dope(s, *p)
            new_points = new_without(points, p)
            #if type(result) is int:
            walk_points3(s, new_points, isolated_arr, points_to_check_all, points_to_check_any, n=n-1, already_doped=already_doped+1)

#sample_1_table = [walk_points(s1, to_hit1, n=i) for i in range(0, 8)]
#print(sample_1_table)
#print('[', end='')
n = 9
isolated = [0]*(n+1)
print('Running...')
t = time.time()
walk_points3(s2, to_hit2, isolated, points_to_check_all=((2, 2, 2),), points_to_check_any=((1,2,2),(2,1,2),(2,2,1),(3,2,2),(2,3,2),(2,2,3),), n=n)
print('Finished in', time.time()-t)
print(isolated/np.array([reverse_fac(32,i) for i in range(n+1)]))
#print(']', end='')

