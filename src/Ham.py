import sympy as sp
from itertools import groupby, combinations, product
from scipy.misc import comb
from scipy.integrate import solve_ivp
from yaml import load

dim = 2 # spatial dimension
num = 3 # particle number
n = dim*num
pairs = int(comb(num, 2, exact=True))

def reduce_output(func, item, *args, **kwargs):
    """
    simple function to reduce output from existing functions

    if func returns an iterable - just return item
    """
    def inner_func(*args, **kwargs):
        return func(*args, **kwargs)[item]
    return inner_func

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


def T(coords, system_vars):
    """
    generic kinetic energy term
    """
    ke = 0
    for particle in range(num):
        for dimension in range(dim):
            ke += coords['p_{}_{}'.format(dimension,particle)]**2/(2*system_vars['mass_{}'.format(particle)])
    return ke

def V(coords, system_vars):
    """
    pairwise potential energy term
    """
    pot = 0
    for pair in combinations(range(num), 2):
        # loop over pairs of particles
        c6 = system_vars['C6_{}_{}'.format(pair[0], pair[1])]
        c12 = system_vars['C12_{}_{}'.format(pair[0], pair[1])]
        r2 = sum([(coords['q_{}_{}'.format(i, pair[0])]-coords['p_{}_{}'.format(i, pair[1])])**2 for i in range(dim)])
        pot += lj(r2,c12,c6)
    return pot

def get_coords(coords):
    return [i for i in coords.itervalues()]

def load_yaml(filen):
    """
    load a yaml file and return the json object
    """
    with open('{}.yml'.format(filen), 'r') as open_file:
        return_dict = load(open_file)
    return return_dict

# define canonical ps and qs
ps = sp.var('p_:{}_:{}'.format(dim, num))
qs = sp.var('q_:{}_:{}'.format(dim, num))
t = sp.var('t')
coords = {i.name:i for i in ps+qs}
system_vars = load_yaml('input')
#define system varibales in yaml file
# mass = sp.var('m_:{}'.format(num))
# C12 = {k.name:k for k in [sp.var('C12_{}_{}'.format(i,j)) for i,j in combinations(range(num), 2)]}
# C6 = {k.name:k for k in [sp.var('C6_{}_{}'.format(i,j)) for i,j in combinations(range(num), 2)]}
# system_vars = {'mass': mass, 'dim': dim, 'num': num, 'C12':C12, 'C6':C6}
print system_vars
KE = T(coords, system_vars)
POT = V(coords, system_vars)
H = KE+POT
HJacp = sp.Matrix([H]).jacobian(ps)
HJacq = -sp.Matrix([H]).jacobian(qs)
RHS = sp.Matrix(sp.BlockMatrix([[HJacp, HJacq]]))
func = reduce_output(sp.lambdify((t,(ps+qs)), RHS.subs(system_vars)), 0)

import pdb
pdb.set_trace()
