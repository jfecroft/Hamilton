from generic import reduce_output, load_yaml
from Ham import Hamilton, T, V
from pprint import pprint
from math import sin, cos, sqrt, pi
from scipy.constants import physical_constants
au2wavenumbers = physical_constants["hartree-inverse meter relationship"][0]/100.0
kelvin2hartree = physical_constants["kelvin-hartree relationship"][0]

def gen_initial_condition(theta, nrg, R, r, mass_0, mass_1, mass_2, **kwargs):
    """
    generate initial conditions for atom+dimer run
    not quite right want equal and opposite momentum not velocity
    """
    ic = {}
    dy = r*sin(theta)
    dx = r*cos(theta)

    ic['q_x_0'] = -R/2.0
    ic['q_y_0'] = 0

    ic['q_x_1'] = R/2.0 + dx
    ic['q_y_1'] = dy

    ic['q_x_2'] = R/2.0 - dx
    ic['q_y_2'] = -dy

    total_mass = mass_0 + mass_1 + mass_2
    p0 = sqrt(2.0*nrg*mass_0**2/total_mass)
    p1 = sqrt(2.0*nrg*mass_1**2/total_mass)
    p2 = sqrt(2.0*nrg*mass_2**2/total_mass)
    ic['p_x_0'] = p0
    ic['p_x_1'] = -p1
    ic['p_x_2'] = -p2
    return ic

data = load_yaml('input')
num = data['system_def']['num']
dim = data['system_def']['dim']
# initial_condition = data['initial_condition']
system_vars = data['system_vars']
Ham = Hamilton(num, dim, T, V, **system_vars)
r = 8.77967500274974/2
R = 200
nrg = 10.0*kelvin2hartree
theta = pi/1.345
ic = gen_initial_condition(theta, nrg, R, r, **system_vars)
pprint(ic)
sol = Ham.prop(10000000, ic)
pprint(sol)
