import dopelib
import chainlib

L = 3

s1 = dopelib.make_bonds_array(L, L, L)
s2 = dopelib.make_bonds_array(L, L+1, L)

to_hit1 = [(i, j, k) for i in range(L - 1) for j in range(L - 1) for k in range(L - 1)]
to_hit2 = [(i, j, k) for i in range(L - 1) for j in range(L+1 - 1) for k in range(L - 1)]


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

def walk_points(sample, points, points_to_check=((1,1,1),), n=0):
    isolated = 0
    for p in points:
        s = sample.copy()
        dopelib.dope(s, *p)
        if n > 0:
            new_points = new_without(points, p)
            isolated += walk_points(s, new_points, points_to_check, n=n-1)
        elif any(map(isolated_checker(s), points_to_check)):
            isolated += 1
    #print(isolated)
    return isolated

#sample_1_table = [walk_points(s1, to_hit1, n=i) for i in range(0, 8)]
#print(sample_1_table)
print('[', end='')
for i in range(0,8):
    print(walk_points(s2, to_hit2, points_to_check=((1, 1, 1),), n=i), end=',')
print(']', end='')
sample_2_table = []
print(sample_2_table)
