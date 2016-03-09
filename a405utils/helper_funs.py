import numpy as np
from collections import namedtuple

def test_scalar(*args):
    """
      return true if every argument is a scalar
    """
    isscalar=True
    for item in args:
        isscalar = isscalar & np.isscalar(item)
    return isscalar

    
def make_tuple(in_dict,tupname='values'):
    """
     make a named tuple from a dictionary
    """
    the_tup = namedtuple(tupname, in_dict.keys())
    the_tup = the_tup(**in_dict)
    return the_tup

def find_centers(x):
    """
    return a vector of bin centers given the bin edges
    """
    center = (x[1:] + x[:-1])/2.
    return center

