#!/usr/bin/env python
"""
  

"""
import numpy as np
from scipy import optimize

def find_interval(func, x, *args):
    """
      starting from a 2% difference, move out from a 
      point until func changes sign

      input func(x,*args): function to zero over x, with optional arguments
            x:  variable to search over for root of func(x,*args)
      output
            left,right  brackets for root 
    """
    if x == 0.:
        dx = 1./50.
    else:
        dx = x/50.
        
    maxiter = 40
    twosqrt = np.sqrt(2)

    failed=True
    for i in range(maxiter):
        dx = dx*twosqrt
        a = x - dx
        fa = func(a, *args)
        b = x + dx
        fb = func(b, *args)
        if (fa*fb < 0.):
            failed=False
            break
    if failed:
        raise ValueError("Couldn't find a suitable range.")
    return (a, b)
        

def fzero(the_func, root_bracket, *args, **parms):
    """
        simple wrapper for optimize.zeros.brenth

        the_func is the function we wish to find the zeros of
        root_bracket is an initial guess of the zero location 
        root_bracket must be a sequence of two floats specifying a range 
        (the_func must differ in sign when evaluated at these points)
        use find_interval to spot a sign change and set the bracket

        *args contains any other arguments needed for the_func
        **parms is passed to brenth and can be xtol (allowable error) or maxiter (max number of iterations.)
        see module tests for usage
    """
    answer=optimize.zeros.brenth(the_func, root_bracket[0], root_bracket[1], *args, **parms)
    return answer
    
def runtests():
    the_zero=fzero(np.sin, [12,13])*180./np.pi  #expecting 720 degrees
    np.testing.assert_almost_equal(the_zero,720.)
    the_zero=fzero(np.sin,[18,20], xtol=1.e-300, maxiter=80)*180./np.pi
    np.testing.assert_almost_equal(the_zero,1080.)
    brackets=find_interval(np.sin,25)
    the_zero=fzero(np.sin,brackets, xtol=1.e-300, maxiter=80)*180./np.pi
    np.testing.assert_almost_equal(the_zero,1440.)
    
 
if __name__=="__main__":
    runtests()
    
