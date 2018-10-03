"""
generic functions
"""
from yaml import load

def reduce_output(func, item, *args, **kwargs):
    """
    simple function to reduce output from existing functions

    if func returns an iterable - just return item
    """
    def inner_func(*args, **kwargs):
        return func(*args, **kwargs)[item]
    return inner_func

def load_yaml(filen):
    """
    load a yaml file and return the json object
    """
    with open('{}.yml'.format(filen), 'r') as open_file:
        return_dict = load(open_file)
    return return_dict


