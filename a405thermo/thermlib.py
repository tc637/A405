"""
  Module for functions to calculate vapor pressure,
  thetae, etc.
"""

from collections import Iterable
import numpy as np
import doctest


from a405thermo.constants import constants as c
from a405thermo import rootfinder as rf
from a405utils.helper_funs import test_scalar
import numpy.testing as ntest
np.set_printoptions(precision=4)

class constants:
   """  A list of constants relevant to atmospheric science.

   No constructor, class variables only

   References
   ----------


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
    ----------

    temp : float or array_like
           Temperature of parcel (K).

    Returns
    -------

    esatOut : float or list
        Saturation water vapour pressure (Pa).

    Examples
    --------

    >>> find_esat(300.)
    3534.5196668891358
    >>> find_esat([300., 310.])
    array([ 3534.5197,  6235.5322])

    References
    ----------
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

def zero_find_rs(Tstart,rsat,press):
    """
       rootfind the temp that produces rsat at press.

       input: temp (K)
              rsat (kg/kg)
              press (hPa)

        output: temp (K)
    """
    brackets=rf.find_interval(zero_rs,Tstart,rsat,press)
    temp = rf.fzero(zero_rs,brackets,rsat,press)
    return temp
 
def find_theta(temp,press,wv=0):
    """
    Computes potential temperature.
    Allows for either temp,p or T,p,wv as inputs.

    Parameters
    ----------

    temp : float
        Temperature (K).

    press : float
        Pressure (Pa).

    wv : float, optional
        Vapour mixing ratio (kg,kg). Can be appended as an argument
        in order to increase precision of returned 'theta' value.


    Returns
    -------

    thetaOut : float
        Potential temperature (K).

    
    Raises
    ------

    NameError
        If an incorrect number of arguments is provided.
    
    
    References
    ----------
    Emanuel p. 111 4.2.11


    Examples
    --------
    >>> find_theta(300., 8.e4) # Only 'temp' and 'p' are input.
    319.72798180767984
    >>> find_theta(300., 8.e4, wv=0.001) # 'temp', 'p', and 'wv' all input.
    319.72309475657323
    
    """

    power = c.Rd / c.cpd * (1. -0.24 * wv);
    thetaOut = temp * (c.p0 / press) ** power;
    return thetaOut

 
def convertSkewToTemp(xcoord, press, skew):
    """
    convertSkewToTemp(xcoord, press, skew)

    Determines temperature from knowledge of a plotting coordinate
    system and corresponding plot skew.
    
    Parameters
    ----------
    xcoord : int
        X coordinate in temperature plotting coordinates.
    press : float
        Pressure (hPa).
    skew : int
        Skew of a given coordinate system.

    Returns
    ----
    Temp : float
        Converted temperature in degC.

    Examples
    --------
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
    ----------
    Temp : float
        Temperature (degC)
    press : float
        Pressure (hPa).
    skew : int
        Designated skew factor of temperature.

    Returns
    ----
    tempOut : float
        Converted temperature (degC).

    Examples
    --------
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
    ----------
    Temp : float
        Temperature (K).
    press : float
        Pressure (Pa).


    Returns
    ----
    thetaep : float
        Pseudo equivalent potential temperature (K).


    Notes
    ---
    It should be noted that the pseudo equivalent potential
    temperature (thetaep) of an air parcel is not a conserved
    variable.


    References
    ----------
    Emanuel 4.7.9 p. 132


    Examples
    --------

    >>> find_thetaes(300., 8.e4)
    412.97362667593831

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
    Calculates the pseudo equivalent potential temperature of a
    parcel. 


    Parameters
    ----------

    Td : float
        Dewpoint temperature (K).

    T : float
        Temperature (K).

    p : float
        Pressure (Pa).


    Returns
    -------

    thetaepOut : float
        Pseudo equivalent potential temperature (K).


    Notes
    -----
    Note that the pseudo equivalent potential temperature of an air
    parcel is not a conserved variable.


    References
    ----------
    Emanuel 4.7.9 p. 132


    Examples
    --------
    >>> find_thetaep(280., 300., 8.e4) # Parcel is unsaturated.
    344.99830738253371

    >>> find_thetaep(300., 280., 8.e4) # Parcel is saturated.
    321.5302927767795

    """
    if Td < T:
        #parcel is unsaturated
        [Tlcl, plcl] = find_lcl(Td, T, p);
        rv = find_rsat(Td, p);
    else:
        #parcel is saturated -- prohibit supersaturation with Td > T
        Tlcl = T;
        rv = find_rsat(T, p);
    
    thetaval = find_theta(T, p, rv);
    thetaepOut = thetaval * np.exp(rv * (1 + 0.81 * rv) \
                                   * (3376. / Tlcl - 2.54));
    #
    # peg this at 450 so rootfinder won't blow up
    #

    if (thetaepOut > 450.):
        thetaepOut = 450;

    return thetaepOut

 
def find_lcl(Td, T, p):
    """
    find_lcl(Td, T, p)

    Finds the temperature and pressure at the lifting condensation
    level (LCL) of an air parcel.

    Parameters
    ----------
    Td : float
        Dewpoint temperature (K).

    T : float
        Temperature (K).

    p : float
        Pressure (Pa)

    Returns
    -------

    Tlcl : float
        Temperature at the LCL (K).
    plcl : float
        Pressure at the LCL (Pa).

    Raises
    ------

    NameError
        If the air is saturated at a given Td and T (ie. Td >= T)
    
    Examples
    --------

    >>> [Tlcl, plcl] =  find_lcl(280., 300., 8.e4)
    >>> print(np.array([Tlcl, plcl]))
    [   275.7625  59518.9287]
    >>> find_lcl(300., 280., 8.e4)
    Traceback (most recent call last):
        ...
    NameError: parcel is saturated at this pressure

    References
    ----------
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
    cp = c.cpd + r * c.cpv;
    logplcl = np.log(p) + cp / (c.Rd * (1 + r / c.eps)) * \
              np.log(Tlcl / T);
    plcl = np.exp(logplcl);

    return Tlcl, plcl

def find_rvrl(Temp, rT, press):
    """
    Computes the vapour and liquid water mixing ratios.

    Parameters
    ----------
    Temp : float
        Temperature (K).
    rT : float
        Total water mixing ratio (kg/kg).
    press : float
        Pressure (Pa).


    Returns
    ----
    rv : float
        Water vapour mixing ratio (kg/kg).
    rl : float
        Liquid water mixing ratio (kg/kg).


    Raises
    ----
    AssertionError
        If any of the inputs are in vector form.

    Examples
    --------

    >>> print(np.array(find_rvrl(250., 0.01, 8.e4))*1.e3)
    [ 0.7433  9.2567]

    >>> print(np.array(find_rvrl(305., 0.01, 8.e4))*1.e3)
    [ 10.   0.]

    """
    rsVal = find_rsat(Temp, press)
    if rsVal > rT: #unsaturated
        rv = rT
        rl = 0
    else:  #saturated
        rv = rsVal
        rl = rT - rv
    return rv, rl

def find_Tmoist(thetaE0, press):
    """
    Calculates the temperatures along a moist adiabat.
    
    Parameters
    ----------

    thetaE0 : float
        Initial equivalent potential temperature (K).
    press : float or array_like
        Pressure (Pa).

    Returns
    -------
    Temp : float or array_like
        Temperature (K) of thetaE0 adiabat at 'press'.

    Examples
    --------
    >>> find_Tmoist(300., 8.e4)
    270.59590841970277
    
    >>> find_Tmoist(330., 800)
    290.98965468303

    """
    Tstart=c.Tc
    brackets=rf.find_interval(thetaes_diff,Tstart,thetaE0,press)
    Temp = rf.fzero(thetaes_diff,brackets,thetaE0,press)
    return Temp
    

def thetaes_diff(Tguess, thetaE0, press):
    """
    Evaluates the equation and passes it back to brenth.

    Parameters
    ----------
    Tguess : float
        Trial temperature value (K).
    ws0 : float
        Initial saturated mixing ratio (kg/kg).
    press : float
        Pressure (Pa).

    Returns
    -------
    theDiff : float
        The difference between the values of 'thetaEguess' and
        'thetaE0'. This difference is then compared to the tolerance
        allowed by brenth.
        
    """
    thetaes_guess = find_thetaes(Tguess, press);
    
    #when this result is small enough we're done
    the_diff = thetaes_guess - thetaE0;
    return the_diff

 
 
def tinvert_thetae(thetaeVal, rT, press):
    """
    temp,rv,rl=tinvert_thetae(thetaeVal, rT, press)

    Uses a rootfinder to determine the temperature for which the
    pseudo equivilant potential temperature (thetaepress) is equal to the
    equivilant potential temperature (thetae) of the parcel.

    Parameters
    ----------
    thetaeVal : float
        Thetae of parcel (K).
    rT : float
        Total water mixing ratio (kg/kg).
    press : float
        Pressure of parcel in (Pa).

    Returns
    -------

    theTemp : float
        Temperature for which thetaep equals the parcel thetae (K).
    rv : float
        Vapor mixing ratio of the parcel (kg/kg).
    rl : float
        liquid water mixing ratio of the parcel (kg/kg) at 'press'.

    Raises
    ------
    IOError
        If 'press' is larger than 100000 Pa.

    Examples
    --------

    >>> tinvert_thetae(300., 0.001, 8.e4)
    (278.4050485684102, 0.001, 0)
    
    """
    if press > 1.e5:
        raise IOError('expecting pressure level less than 100000 Pa')
    # The temperature has to be somewhere between thetae
    # (T at surface) and -40 deg. C (no ice).
    Tstart=c.Tc
    brackets=rf.find_interval(find_Tchange,Tstart,thetaeVal,rT,press)
    theTemp = rf.fzero(find_Tchange,brackets,thetaeVal,rT,press)
    rv,rl = find_rvrl(theTemp, rT, press);
    return theTemp,rv,rl

 
def find_Tchange(Tguess, thetaeVal, rT, press):
    rv, rl = find_rvrl(Tguess, rT, press);
    tdGuess = find_Td(rv, press);
    # Iterate on Tguess until this function is
    # zero to within tolerance.
    return thetaeVal - find_thetaep(tdGuess,Tguess,press);

def find_Td(wv, press):
    """
    Calculates the due point temperature of an air parcel.

    Parameters
    ----------

    wv : float
        Mixing ratio (kg/kg).
    press : float
        Pressure (Pa).

    Returns
    -------

    Td : float
        Dew point temperature (K).

    Examples
    --------

    >>> find_Td(0.001, 8.e4)
    253.39429263963504

    References
    ----------

    Emanuel 4.4.14 p. 117
    
    """
    e = wv * press / (c.eps + wv);
    denom = (17.67 / np.log(e / 611.2)) - 1.;
    Td = 243.5 / denom;
    Td = Td + 273.15;
    return Td

 
def test_therm():
    Tlcl, plcl =  find_lcl(280., 300., 8.e4)
    ntest.assert_almost_equal(Tlcl,275.7625,decimal=3)
    ntest.assert_almost_equal(plcl,59518.928,decimal=2)
    ntest.assert_almost_equal(find_thetaep(280., 300., 8.e4),344.998307,decimal=5) # Parcel is unsaturated.
    ntest.assert_almost_equal(find_thetaep(300., 280., 8.e4),321.53029,decimal=5) # Parcel is saturated.
    ntest.assert_almost_equal(find_esat(300.),3534.51966,decimal=2)
    ntest.assert_almost_equal(find_thetaes(300., 8.e4),412.9736,decimal=4)
    ntest.assert_allclose(find_esat([300., 310.]),[3534.51966, 6235.53218])
    ntest.assert_almost_equal(find_Tmoist(300., 8.e4),270.59590841970277)
    ntest.assert_almost_equal(find_Tmoist(330., 8.e4), 282.92999,decimal=4)

    
if __name__ == "__main__":
   test_therm()
   doctest.testmod()
   
