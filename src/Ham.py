import sympy as sp
from itertools import groupby, combinations, product
from scipy.misc import comb
from scipy.integrate import solve_ivp
from yaml import load
from generic import reduce_output

dim = 2 # spatial dimension
num = 3 # particle number
n = dim*num
pairs = int(comb(num, 2, exact=True))


def lj(r2, C12, C6):
    """
    lennard jones potentail r2 is r**2
    """
    return C12/r2**6 - C6/r2**3

def T(coords, system_def, system_vars):
    """
    generic kinetic energy term
    """
    ke = 0
    for particle in range(system_def['num']):
        for dimension in system_def['dim']:
            ke += coords['p_{}_{}'.format(dimension,particle)]**2/(2*system_vars['mass_{}'.format(particle)])
    return ke

def V(coords, system_def, system_vars):
    """
    pairwise potential energy term
    """
    pot = 0
    for pair in combinations(range(system_def['num']), 2):
        # loop over pairs of particles
        c6 = system_vars['C6_{}_{}'.format(pair[0], pair[1])]
        c12 = system_vars['C12_{}_{}'.format(pair[0], pair[1])]
        r2 = sum([(coords['q_{}_{}'.format(i, pair[0])]-coords['p_{}_{}'.format(i, pair[1])])**2 for i in system_def['dim']])
        pot += lj(r2,c12,c6)
    return pot

def load_yaml(filen):
    """
    load a yaml file and return the json object
    """
    with open('{}.yml'.format(filen), 'r') as open_file:
        return_dict = load(open_file)
    return return_dict


class Hamilton:
    def __init__(self, T, V, input_file):
        self.data = load_yaml('input')
        self.system_def = self.data['system_def']
        self.dim = self.system_def['dim']
        self.num = self.system_def['num']
        self.system_vars = self.data['system_vars']
        self.coords = self.create_coords(self.dim, self.num)
        self.coord_dict = {i.name:i for i in self.coords}
        self.T = T(self.coord_dict, self.system_def, self.system_vars)
        self.V = V(self.coord_dict, self.system_def, self.system_vars)
        self.H = self.T + self.V
        self.ps = [i for i in self.coords if i.name[0] == 'p']
        self.qs = [i for i in self.coords if i.name[0] == 'q']

    def create_coords(self, dim, num):
        coords = []
        for i, j, k in product(dim, range(num), ('p', 'q')):
            coords.append(sp.var('{}_{}_{}'.format(k,i,j)))
        return tuple(coords)

    def prop(self, time):
        HJacp = sp.Matrix([self.H]).jacobian(self.ps)
        HJacq = -sp.Matrix([self.H]).jacobian(self.qs)
        RHS = sp.Matrix(sp.BlockMatrix([[HJacp, HJacq]]))
        t = sp.var('t')
        func = reduce_output(sp.lambdify((t,(self.ps+self.qs)), RHS), 0)
        # func is defind as ps+qs therefore we must pass qs + ps - see hamilton equations
        print solve_ivp(func, (0,100), range(12))
        import pdb
        pdb.set_trace()




Ham = Hamilton(T, V, 'input')
Ham.prop(10)
# define canonical ps and qs
data = load_yaml('input')
system_vars = data['system_vars']
system_def = data['system_def']
coords = create_coords(system_def['dim'], system_def['num'])
ps = [i for i in coords if i.name[0] == 'p']
qs = [i for i in coords if i.name[0] == 'q']
coord_dict = {i.name:i for i in coords}
#define system varibales in yaml file
# mass = sp.var('m_:{}'.format(num))
# C12 = {k.name:k for k in [sp.var('C12_{}_{}'.format(i,j)) for i,j in combinations(range(num), 2)]}
# C6 = {k.name:k for k in [sp.var('C6_{}_{}'.format(i,j)) for i,j in combinations(range(num), 2)]}
# system_vars = {'mass': mass, 'dim': dim, 'num': num, 'C12':C12, 'C6':C6}
print system_vars
KE = T(coord_dict,system_def, system_vars)
POT = V(coord_dict, system_def, system_vars)
H = KE+POT

