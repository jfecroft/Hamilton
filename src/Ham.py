import sympy as sp
from itertools import groupby, combinations, product
from scipy.misc import comb
from scipy.integrate import solve_ivp
from generic import reduce_output, load_yaml
import numpy as np

dim = 2 # spatial dimension
num = 3 # particle number
n = dim*num
pairs = int(comb(num, 2, exact=True))


def lj(r2, C12, C6):
    """
    lennard jones potentail r2 is r**2
    """
    return C12/r2**6 - C6/r2**3

def T(coords, num, dim, **kwargs):
    """
    generic kinetic energy term
    """
    ke = 0
    for particle in range(num):
        for dimension in dim:
            ke += coords['p_{}_{}'.format(dimension,particle)]**2/(2*kwargs['mass_{}'.format(particle)])
    return ke

def V(coords, num, dim, **kwargs):
    """
    pairwise potential energy term
    """
    pot = 0
    for pair in combinations(range(num), 2):
        # loop over pairs of particles
        c6 = kwargs['C6_{}_{}'.format(pair[0], pair[1])]
        c12 = kwargs['C12_{}_{}'.format(pair[0], pair[1])]
        r2 = sum([(coords['q_{}_{}'.format(i, pair[0])]-coords['q_{}_{}'.format(i, pair[1])])**2 for i in dim])
        pot += lj(r2,c12,c6)
    return pot


class Hamilton:
    def __init__(self, num, dim, T, V, **kwargs):
        self.dim = dim
        self.num = num
        self.ps, self.qs = self.create_coords(self.dim, self.num)
        self.coords = tuple(self.ps+self.qs)
        self.coord_dict = {i.name:i for i in self.coords}
        self.T = T(self.coord_dict, num, dim, **kwargs)
        self.V = V(self.coord_dict, num, dim, **kwargs)
        self.H = self.T + self.V

    def create_coords(self, dim, num):
        ps = []
        for i, j in product(dim, range(num)):
            ps.append(sp.var('p_{}_{}'.format(i,j)))
        qs = []
        for i, j in product(dim, range(num)):
            qs.append(sp.var('q_{}_{}'.format(i,j)))
        return tuple(ps), tuple(qs )

    def create_initial_condition(self, initial_condition, default=0):
        """
        create initial condition vector from input file and ps and qs
        consistant with the RHS definition


        if not in dict use default - avoids lots of zeros in intial conditions
        """
        return [initial_condition.get(i.name, default) for i in self.coords]

    def prop(self, time, initial_condition, rtol=1.0e-6):
        H = sp.Matrix([self.H])
        RHS = sp.Matrix(sp.BlockMatrix([[H.jacobian(self.qs), H.jacobian(self.ps)]]))
        t = sp.var('t')
        dydt_func = reduce_output(sp.lambdify((t,(self.coords)), RHS), 0)
        nrg_func = sp.lambdify((t,(self.coords)), self.H, initial_condition)

        # create some events
        y0 = self.create_initial_condition(initial_condition)
        events = []
        #want to set rtol = nrg_tol*0.1
        inital_energy = nrg_func(0, y0)
        nrg_condition = rtol_func(nrg_func, inital_energy, 1.0e-5)
        nrg_condition.terminal = False
        events.append(nrg_condition)
        sol = solve_ivp(dydt_func, (0,time), y0, rtol=rtol, events=events)
        final = sol['y'][:,-1]
        final_energy = nrg_func(0, final)
        energy_conservation = (inital_energy-final_energy)/inital_energy
        print 'change in total energy {}%'.format(energy_conservation*100)
        print (self.coords)
        traj = np.vstack((sol['t'],sol['y']))
        np.savetxt('traj.dat', traj.T)


def rtol_func(func, val, rtol, *args, **kwargs):
    """
    simple function - given a function which returns a value
    return a function which return the fractional change
    """
    def inner_func(*args, **kwargs):
        return abs(100*(func(*args, **kwargs) - val)/val) - rtol
    return inner_func


data = load_yaml('input')
num = data['system_def']['num']
dim = data['system_def']['dim']
initial_condition = data['initial_condition']
system_vars = data['system_vars']
Ham = Hamilton(num, dim, T, V, **system_vars)
Ham.prop(10, initial_condition, 1.0e-6)
# # define canonical ps and qs
# data = load_yaml('input')
# system_vars = data['system_vars']
# system_def = data['system_def']
# coords = create_coords(system_def['dim'], system_def['num'])
# ps = [i for i in coords if i.name[0] == 'p']
# qs = [i for i in coords if i.name[0] == 'q']
# coord_dict = {i.name:i for i in coords}
# #define system varibales in yaml file
# # mass = sp.var('m_:{}'.format(num))
# # C12 = {k.name:k for k in [sp.var('C12_{}_{}'.format(i,j)) for i,j in combinations(range(num), 2)]}
# # C6 = {k.name:k for k in [sp.var('C6_{}_{}'.format(i,j)) for i,j in combinations(range(num), 2)]}
# # system_vars = {'mass': mass, 'dim': dim, 'num': num, 'C12':C12, 'C6':C6}

