from collections import Iterable
import numpy as np
import os,site

from a405thermo.constants import constants as c
from a405thermo import rootfinder as rf
import numpy.testing as ntest
import nose

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
    >>> find_esat(300.)
    3534.5196668891358
    >>> find_esat([300., 310.])
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
    esat = find_esat(temp)
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

def find_thetaes(Temp, press):
    """
    find_thetaes(Temp, press)

    Calculates the pseudo equivalent potential temperature of an air
    parcel.

    Parameters
    - - - - - -
    Temp : float
        Temperature (K).
    press : float
        Pressure (Pa).


    Returns
    - - - -
    thetaep : float
        Pseudo equivalent potential temperature (K).


    Notes
    - - -
    It should be noted that the pseudo equivalent potential
    temperature (thetaep) of an air parcel is not a conserved
    variable.


    References
    - - - - - -
    Emanuel 4.7.9 p. 132


    Examples
    - - - - -
    >>> test.assert_almost_equal(find_thetaes(300., 8.e4),412.9736,decimal=4)
    """
    # The parcel is saturated - prohibit supersaturation with Td > T.
    Tlcl = Temp;
    rv = find_rsat(Temp, press);
    thetaval = find_theta(Temp, press, rv);
    thetaep = thetaval * np.exp(rv * (1 + 0.81 * rv) * \
                                (3376. / Tlcl - 2.54))
    #
    # peg this at 450 so rootfinder won't blow up
    #
    if thetaep > 450.:
        thetaep = 450;
    return thetaep


def find_thetaep(Td, T, p):
    """
    thetaep(Td, T, p)

    Calculates the pseudo equivalent potential temperature of a
    parcel. 


    Parameters
    - - - - - -
    Td : float
        Dewpoint temperature (K).
    T : float
        Temperature (K).
    p : float
        Pressure (Pa).


    Returns
    - - - -
    thetaepOut : float
        Pseudo equivalent potential temperature (K).


    Notes
    - - -
    Note that the pseudo equivalent potential temperature of an air
    parcel is not a conserved variable.


    References
    - - - - - -
    Emanuel 4.7.9 p. 132


    Examples
    - - - - -
    >>> test.assert_almost_equal(thetaep(280., 300., 8.e4),344.998307,decimal=5) # Parcel is unsaturated.
    >>> test.assert_almost_equal(thetaep(300., 280., 8.e4),321.53029,decimal=5) # Parcel is saturated.
    """
    if Td < T:
        #parcel is unsaturated
        [Tlcl, plcl] = LCLfind(Td, T, p);
        rv = find_rsat(Td, p);
    else:
        #parcel is saturated -- prohibit supersaturation with Td > T
        Tlcl = T;
        rv = find_rsat(T, p);
    
    # $$$   disp('inside theate')
    # $$$   [Td,T,wv]
    thetaval = find_theta(T, p, rv);
    thetaepOut = thetaval * np.exp(rv * (1 + 0.81 * rv) \
                                   * (3376. / Tlcl - 2.54));
    #
    # peg this at 450 so rootfinder won't blow up
    #

    if (thetaepOut > 450.):
        thetaepOut = 450;

    return thetaepOut

 
def LCLfind(Td, T, p):
    """
    LCLfind(Td, T, p)

    Finds the temperature and pressure at the lifting condensation
    level (LCL) of an air parcel.

    Parameters
    - - - - - -
    Td : float
        Dewpoint temperature (K).
    T : float
        Temperature (K).
    p : float
        Pressure (Pa)

    Returns
    - - - -
    Tlcl : float
        Temperature at the LCL (K).
    plcl : float
        Pressure at the LCL (Pa).

    Raises
    - - - -
    NameError
        If the air is saturated at a given Td and T (ie. Td >= T)
    
    Examples
    - - - - -
    >>> [Tlcl, plcl] =  LCLfind(280., 300., 8.e4)
    >>> print [Tlcl, plcl]
    [275.76250387361404, 59518.928699453245]
    >>> LCLfind(300., 280., 8.e4)
    Traceback (most recent call last):
        ...
    NameError: parcel is saturated at this pressure

    References
    - - - - - -
    Emanuel 4.6.24 p. 130 and 4.6.22 p. 129
    
    """
    hit = Td >= T;
    if hit is True:
        raise NameError('parcel is saturated at this pressure');

    e = find_esat(Td);
    ehPa = e * 0.01; #Bolton's formula requires hPa.
    # This is is an empircal fit from for LCL temp from Bolton, 1980 MWR.
    Tlcl = (2840. / (3.5 * np.log(T) - np.log(ehPa) - 4.805)) + 55.;

    r = c.eps * e / (p - e);
    #disp(sprintf('r=%0.5g',r'))
    cp = c.cpd + r * c.cpv;
    logplcl = np.log(p) + cp / (c.Rd * (1 + r / c.eps)) * \
              np.log(Tlcl / T);
    plcl = np.exp(logplcl);
    #disp(sprintf('plcl=%0.5g',plcl))

    return Tlcl, plcl
 
def test_therm():
    Tlcl, plcl =  LCLfind(280., 300., 8.e4)
    ntest.assert_allclose(Tlcl,275.7625)
    ntest.assert_allclose(plcl,59518.928)
    ntest.assert_almost_equal(find_thetaep(280., 300., 8.e4),344.998307,decimal=5) # Parcel is unsaturated.
    ntest.assert_almost_equal(find_thetaep(300., 280., 8.e4),321.53029,decimal=5) # Parcel is saturated.
    ntest.assert_allclose(find_esat(300.),3534.51966)
    ntest.assert_allclose(find_esat([300., 310.]),[3534.51966, 6235.53218])
    ntest.assert_almost_equal(find_thetaes(300., 8.e4),412.9736,decimal=4)

if __name__ == "__main__":
   test_therm()
   
