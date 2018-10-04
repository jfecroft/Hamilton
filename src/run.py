from generic import reduce_output, load_yaml
from Ham import Hamilton, T, V
from pprint import pprint
from math import sin, cos, sqrt, pi
from scipy.constants import physical_constants
au2wavenumbers = physical_constants["hartree-inverse meter relationship"][0]/100.0
kelvin2hartree = physical_constants["kelvin-hartree relationship"][0]

def gen_initial_condition():
    ic = {}
    R = 2000
    theta = pi/2
    r = 8.77967500274974/2
    r = 2000
    dy = r*sin(theta)
    dx = r*cos(theta)
    ic['q_x_0'] = -R/2.0
    ic['q_y_0'] = 0

    ic['q_x_1'] = R/2.0 + dx
    ic['q_y_1'] = dy

    ic['q_x_2'] = R/2.0 - dx
    ic['q_y_2'] = -dy
    ecol = 1.0*kelvin2hartree
    p = ecol*2.0*12652.668976691984/3.0
    p = sqrt(p)
    ic['p_x_0'] = p
    ic['p_x_1'] = -p
    ic['p_x_2'] = -p
    return ic

data = load_yaml('input')
num = data['system_def']['num']
dim = data['system_def']['dim']
initial_condition = data['initial_condition']
system_vars = data['system_vars']
Ham = Hamilton(num, dim, T, V, **system_vars)
ic = gen_initial_condition()
initial_condition.update(ic)
pprint(initial_condition)
sol = Ham.prop(10000000, initial_condition)
pprint(sol)
