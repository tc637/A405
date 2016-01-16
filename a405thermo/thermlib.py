from collections import Iterable
import numpy as np

class constants:
   """ 
   A list of constants relevant to atmospheric science.

   References
   - - - - - -
   Emanuel appendix 2
   
   """
   Tc = 273.15 #0 deg C in Kelvin
   eps = 0.622
   p0 = 1.e5
   eps = 0.622
   lv0 = 2.501e6 # J/kg
   Rv = 461.50 # J/kg/K
   Rd = 287.04 # J/kg/K
   cpv = 1870. # Heat capacity of water vapor (J/kg/K)
   cl = 4190.  # Heat capacity of liquid water (J/kg/K)
   cpd = 1005.7 # Heat capacity of dry air (J/kg/K)
   g0 = 9.8      # m/s^2
   D = 2.36e-5# Diffusivity m^2/s^1 
             # Note: fairly strong function of temperature
             #       and pressure -- this is at 100kPa, 10degC
   rhol = 1000.

def esat(temp):
    """
    esat(temp)

    Calculates the saturation water vapor pressure over a flat
    surface of water at temperature 'temp'.

    Parameters
    - - - - - -
    temp : float or array_like
           Temperature of parcel (K).

    Returns
    - - - -
    esatOut : float or list
        Saturation water vapour pressure (Pa).

    Examples
    - - - - -
    >>> esat(300.)
    3534.5196668891358
    >>> esat([300., 310.])
    [3534.5196668891358, 6235.5321818976754]

    References
    - - - - - -
    Emanuel 4.4.14 p. 117
      
    """
    # determine if temp has been input as a vector
    is_scalar=True
    if isinstance(temp,Iterable):
        is_scalar = False
    temp=np.atleast_1d(temp)
    Tc = temp - 273.15
    esatOut = 611.2 * np.exp(17.67 * Tc / (Tc + 243.5))
    # if temp is a vector
    if is_scalar:
        esatOut = esatOut[0]
    return esatOut
   
