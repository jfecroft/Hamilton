"""
Generic hamiltonian propagator - as general as possible
NO units or system dependant stuff here
"""

import sympy as sp
from itertools import groupby, combinations, product
from scipy.misc import comb
from scipy.integrate import solve_ivp
from generic import reduce_output, load_yaml
import numpy as np
from pprint import pprint


def lj(r2, C12, C6):
    """
    lennard jones potential r2 is r**2
    """
    return C12/r2**6 - C6/r2**3


def T(coords, num, dim, **kwargs):
    """
    generic kinetic energy term
    """
    ke = 0
    for particle in range(num):
        for dimension in dim:
            ke += coords['p_{}_{}'.format(dimension, particle)]**2/(
                2*kwargs['mass_{}'.format(particle)])
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
        r2 = sum([(coords['q_{}_{}'.format(i, pair[0])] -
                   coords['q_{}_{}'.format(i, pair[1])])**2 for i in dim])
        pot += lj(r2, c12, c6)
    return pot


class Hamilton:
    def __init__(self, num, dim, T, V, **kwargs):
        """
        initialise the hamiltonian and the canonical coordiantes
        """
        self.dim = dim
        self.num = num
        self.ps, self.qs = self.create_coords(self.dim, self.num)
        self.coords = tuple(self.ps+self.qs)
        coord_dict = {i.name: i for i in self.coords}
        self.H = T(coord_dict, num, dim, **kwargs) + V(
                   coord_dict, num, dim, **kwargs)

    def create_coords(self, dim, num):
        ps = []
        for i, j in product(dim, range(num)):
            ps.append(sp.var('p_{}_{}'.format(i, j)))
        qs = []
        for i, j in product(dim, range(num)):
            qs.append(sp.var('q_{}_{}'.format(i, j)))
        return tuple(ps), tuple(qs)

    def create_initial_condition(self, initial_condition, default=0):
        """
        create initial condition vector from input file and ps and qs
        consistant with the order of the coords in self.coords

        if not in dict use default - avoids lots of zeros in intial conditions
        """
        return [initial_condition.get(i.name, default) for i in self.coords]

    def prop(self, time, initial_condition, nrgtol=1.0e-3, rtol=1.0e-5):
        """
        propagate the solution as a function of time
        Variables:
            time - tmax
            initial_condition - dictionary of coord:vals
            nrg_tol - desired energy tolerance
            rtol - relative tolerance passed to solve_ivp -
                   recommeded to be 0.1*nrgtol
        """
        t = sp.var('t')
        H = sp.Matrix([self.H])
        dydt = sp.Matrix(sp.BlockMatrix([[-H.jacobian(self.qs),
                                          H.jacobian(self.ps)]]))
        dydt_func = reduce_output(sp.lambdify((t, (self.coords)), dydt), 0)
        nrg_func = sp.lambdify((t, (self.coords)), self.H, initial_condition)
        y0 = self.create_initial_condition(initial_condition)
        inital_energy = nrg_func(0, y0)
        print 'inital_energy', inital_energy
        print self.coords

        # create some events
        events = []
        if nrgtol:
            # want to set rtol = nrg_tol*0.1
            nrg_condition = rtol_func(nrg_func, inital_energy, nrgtol)
            nrg_condition.terminal = False
            events.append(nrg_condition)
        sol = solve_ivp(dydt_func, (0, time), y0, rtol=rtol, events=events)
        final_y = sol['y'][:, -1]
        final_energy = nrg_func(0, final_y)
        energy_conservation = (inital_energy-final_energy)/inital_energy
        sol.energy_conservation = energy_conservation
        traj = np.vstack((sol['t'], sol['y']))
        np.savetxt('traj.dat', traj.T)
        return sol


def rtol_func(func, val, rtol, *args, **kwargs):
    """
    simple function - given a function which returns a value
    return a function which return the fractional change
    """
    def inner_func(*args, **kwargs):
        return abs((func(*args, **kwargs) - val)/val) - rtol
    return inner_func
