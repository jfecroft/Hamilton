"""
generic functions
"""

def reduce_output(func, item, *args, **kwargs):
    """
    simple function to reduce output from existing functions

    if func returns an iterable - just return item
    """
    def inner_func(*args, **kwargs):
        return func(*args, **kwargs)[item]
    return inner_func

