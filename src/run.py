from generic import reduce_output, load_yaml
from Ham import Hamilton, T, V
from pprint import pprint
data = load_yaml('input')
num = data['system_def']['num']
dim = data['system_def']['dim']
initial_condition = data['initial_condition']
system_vars = data['system_vars']
Ham = Hamilton(num, dim, T, V, **system_vars)
sol = Ham.prop(1000, initial_condition, rtol=1.0e-9, nrgtol=1.0e-3)
pprint(sol)
