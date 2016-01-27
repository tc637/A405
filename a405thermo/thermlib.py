from collections import Iterable
import numpy as np
from a405thermo.constants import constants as c

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

def find_esat(temp):
    """
    find_esat(temp)

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


def zero_rs(temp,rsat,press):
    """
      find the saturation temperature for 
      a given rsat,press, by rootfinding this zero
      
      input: temp (guess) (K)
             rsat (kg/kg)
             press (hPa)
      output: residual
      
      see thompkins 2.20
     
    """
    esat=find_esat(temp)*0.01  #convert to hPa
    residual=rsat - c.eps*esat/(press - esat)
    return residual

def find_rsat(temp,press):
    """
       calculate the saturation mixing ratio at (temp,press)

       input: temp (K)
              press (Pa)
        output: rsat (kg/kg)
    """
    esat = find_esat(temp)*0.01
    rsat=c.eps*esat/(press - esat)
    return rsat

def zero_find_rs(tstart,rsat,press):
    """
       rootfind the temp that produces rsat at press.


       input: temp (K)
              rsat (kg/kg)
              press (hPa)
        output: temp (K)
    """
    brackets=rf.find_interval(zero_rs,tstart,rsat,press)
    temp = rf.fzero(zero_rs,brackets,rsat,press)
    return temp
 
def find_theta(temp,press,wv=0):
    """
    theta(*args)

    Computes potential temperature.
    Allows for either temp,p or T,p,wv as inputs.
    

    Parameters
    - - - - - -
    temp : float
        Temperature (K).
    press : float
        Pressure (Pa).


    Returns
    - - - -
    thetaOut : float
        Potential temperature (K).


    Optional Parameters
    - - - - - - - - -
    wv : float, optional
        Vapour mixing ratio (kg,kg). Can be appended as an argument
        in order to increase precision of returned 'theta' value.
    
    
    Raises
    - - - -
    NameError
        If an incorrect number of arguments is provided.
    
    
    References
    - - - - - -
    Emanuel p. 111 4.2.11


    Examples
    - - - - -
    >>> theta(300., 8.e4) # Only 'temp' and 'p' are input.
    319.72798180767984
    >>> theta(300., 8.e4, wv=0.001) # 'temp', 'p', and 'wv' all input.
    319.72309475657323
    
    """

    power = c.Rd / c.cpd * (1. - 0.24 * wv);
    thetaOut = temp * (c.p0 / press) ** power;
    return thetaOut

 
def convertSkewToTemp(xcoord, press, skew):
    """
    convertSkewToTemp(xcoord, press, skew)

    Determines temperature from knowledge of a plotting coordinate
    system and corresponding plot skew.
    
    Parameters
    - - - - - -
    xcoord : int
        X coordinate in temperature plotting coordinates.
    press : float
        Pressure (hPa).
    skew : int
        Skew of a given coordinate system.

    Returns
    - - - -
    Temp : float
        Converted temperature in degC.

    Examples
    - - - - -
    >>> convertSkewToTemp(300, 8.e4, 30)
    638.6934574096806
    
    """
    Temp = xcoord  + skew * np.log(press);
    return Temp

def convertTempToSkew(Temp, press, skew):
    """
    convertTempToSkew(Temp, press, skew)

    Determines the transformed temperature in plotting coordinates.
    
    Parameters
    - - - - - -
    Temp : float
        Temperature (degC)
    press : float
        Pressure (hPa).
    skew : int
        Designated skew factor of temperature.

    Returns
    - - - -
    tempOut : float
        Converted temperature (degC).

    Examples
    - - - - -
    >>> convertTempToSkew(30., 8.e4, 30)
    -308.69345740968055
    
    """
    
    tempOut = Temp - skew * np.log(press);
    return tempOut
