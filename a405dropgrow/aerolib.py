import numpy as np
from a405thermo.constants import constants as c


def lognormal(x,mu,sigma):
    """
    parameters
    ----------
    x: vector (float)  
      aerosol masses (kg)  (for example)
      
    mu: float
       log(mean mass)
       
    sigma: float
       log(standard deviation)
       
    returns
    -------
    
    out: vector (float)
        lognormal pdf, normalized to 1
        units
    
    """
    out=(1/(x*sigma*np.sqrt(2*np.pi)))*np.exp(-(np.log(x) - mu)**2./(2*sigma**2.))
    return out


def create_koehler(aero,parcel):
    """
    generate a kohler function for specific
    values of the a,b coefficients

    Parameters
    ----------

    a: float (m)
      kelvin term in kohler equation
    b: float (m^3/kg)
      solution term in kohler equation, needs to 
      be multiplied by aerosol mass

    Returns
    -------

    Koheler function(r,m)
    """
    def find_S(r,m):
        """
        Parameters
        ----------
        r: float (m)
           drop radius
        m: float (kg)
           aerosol mass (kg)

        Returns
        -------

        Function to calculate saturation over the curved solution drop
        """
        #
        # reset negative radii to 0.001 micron
        #
        if r < 0.:
            r = 1.e-9
        # use exact koehler equation
        #
        a=(2.*aero.Sigma)/(c.Rv*parcel.Tinit*c.rhol)  #curvature term
        ns = m*aero.vanHoff/aero.Ms
        nw = 4/3.*np.pi*r**3.*c.rhol/aero.Mw
        S = (nw/(ns + nw))*np.exp(a/r)
        return S
    return find_S

def find_koehler_coeffs(aero,parcel):
    """
    given named tuples with the aerosol and parcel
    variables, find the a and b coefficients for the approximate koehler equation
    """
    a=(2.*aero.Sigma)/(c.Rv*parcel.Tinit*c.rhol)  #curvature term
    b=(aero.vanHoff*aero.Mw)/((4./3.)*np.pi*c.rhol*aero.Ms)  #solution term, no mass
    return a,b
