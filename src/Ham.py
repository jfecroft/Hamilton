import sympy as sp
from itertools import groupby, combinations, product
from scipy.misc import comb

dim = 2 # spatial dimension
num = 3 # particle number
n = dim*num
pairs = int(comb(num, 2, exact=True))

def lj(r2, C12, C6):
    """
    lennard jones potentail r2 is r**2
    """
    return C12/r2**6 - C6/r2**3

def is_q(var):
    """
    test if a variable is a q
    """
    return var.name[0] == 'q'

def is_p(var):
    """
    test if a variable is a p
    """
    return var.name[0] == 'p'

def get_dim(var):
    """
    return the dimension of var
    """
    return int(var.name.split('_')[1])

def get_num(var):
    """
    return the paritcle number of var
    """
    return int(var.name.split('_')[2])


def T(coords, mass, **system_vars):
    """
    generic kinetic energy term
    """
    return sum([i**2/(2*mass[get_num(i)]) for i in coords.itervalues() if is_p(i)])

def V(coords, num, dim, C12, C6, **system_vars):
    """
    pairwise potential energy term
    """
    pot = 0
    for pair in combinations(range(num), 2):
        # loop over pairs of particles
        c6 = C6['C6_{}_{}'.format(pair[0], pair[1])]
        c12 = C12['C12_{}_{}'.format(pair[0], pair[1])]
        r2 = sum([(coords['p_{}_{}'.format(i, pair[0])]-coords['p_{}_{}'.format(i, pair[1])])**2 for i in range(dim)])
        pot += lj(r2,c12,c6)
    return pot

# define canonical ps and qs
ps = sp.var('p_:{}_:{}'.format(dim, num))
qs = sp.var('q_:{}_:{}'.format(dim, num))
coords = {i.name:i for i in ps+qs}
#define system varibales
mass = sp.var('m_:{}'.format(num))
C12 = {k.name:k for k in [sp.var('C12_{}_{}'.format(i,j)) for i,j in combinations(range(num), 2)]}
C6 = {k.name:k for k in [sp.var('C6_{}_{}'.format(i,j)) for i,j in combinations(range(num), 2)]}
system_vars = {'mass': mass, 'dim': dim, 'num': num, 'C12':C12, 'C6':C6}
KE = T(coords, **system_vars)
POT = V(coords, **system_vars)
H = KE+POT
import pdb
pdb.set_trace()



qs_dim = [list(j) for i,j  in groupby(qs, lambda x:  x.name[2])] #qs grouped by dim
all_pairs = [j for i in qs_dim for j in combinations(i,2)]


#define H
T = sum([p**2/(2*m) for p,m in zip(ps,mass)])

sp.Matrix([T]).jacobian(ps)


