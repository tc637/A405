"""
calculate the grow history of droplet categories in a constant velocity updraft
"""
import ruamel.yaml as ry
import a405dropgrow.aerolib
from importlib import reload
reload(a405dropgrow.aerolib)
from a405dropgrow.aerolib import lognormal,create_koehler,find_koehler_coeffs
import numpy as np
import a405utils.helper_funs
reload(a405utils.helper_funs)
from a405utils.helper_funs import make_tuple, find_centers
from a405thermo.rootfinder import find_interval, fzero
from a405thermo.constants import constants as c
from a405thermo.thermlib import find_lv,find_esat
from collections import OrderedDict as od
from scipy.integrate import odeint
import pandas as pd
from matplotlib import pyplot as plt

def find_diff(logr,S_target,m):
    """
    zero function for rootfinder

    Parameters
    ----------

    logr: float
       log of radius
    S_target: float
          Satutation ratio to match
    m: float
      aerosol mass (kg)

    Returns
    -------

    Sdiff: float
       difference between target and guess
    """
    r = np.exp(logr)
    return S_target - koehler_fun(r,m)


def wlcalc(var_vec,cloud_tup):
    """
    calculate the liquid water for the distribution

    Parameters
    ----------

    var_vec: vector(float)
           vector of values to be integrated
    cloud_top: namedtuple
           tuple of necessary coefficients
    """
    wl=cloud_tup.ndist*(var_vec[:-3]**3.)
    wl=np.sum(wl)
    wl=wl*4./3.*np.pi*c.rhol
    return wl

def Scalc(var_vec,cloud_tup):
    """
    calculate the environmental saturation

    Parameters
    ----------

    var_vec: vector(float)
        vector of values to be integrated

    cloud_top: namedtuple
           tuple of necessary coefficients

    Returns
    -------

    Sout: float
       environmental saturation

    """
    temp,press,height = var_vec[-3:]
    wl=wlcalc(var_vec,cloud_tup)
    wv=cloud_tup.wt - wl
    e=wv*press/(c.eps + wv)
    Sout=e/find_esat(temp)
    return Sout

def find_derivs(var_vec,the_time,cloud_tup):
    """
    calcuate derivatives of var_vec 

    Parameters
    ----------

    var_vec: vector(float)
        vector of values to be integrated

    the_time: float
       timestep 

    cloud_tup: namedtuple
           tuple of necessary coefficients
    

    Returns
    -------

    deriv_vec: vector(float)
         derivatives of each of var_vec
    
    """
    #print('inside: ',var_vec)
    temp,press,height = var_vec[-3:]
    numrads = len(var_vec) - 3
    dry_radius = cloud_tup.dry_radius
    rho=press/(c.Rd*temp)
    S=Scalc(var_vec,cloud_tup)
    deriv_vec=np.zeros_like(var_vec)
    #dropgrow notes equaton 18 (W&H p. 170)
    for i in range(numrads):
        m=cloud_tup.masses[i]
        if var_vec[i] < dry_radius[i]:
            var_vec[i] = dry_radius[i]
        Seq=koehler_fun(var_vec[i],m)
        rhovr=(Seq*find_esat(temp))/(c.Rv*temp)
        rhovinf=S*find_esat(temp)/(c.Rv*temp)
        deriv_vec[i]=(c.D/(var_vec[i]*c.rhol))*(rhovinf - rhovr)
    deriv_vec[-3]=find_lv(temp)/c.cpd*wlderiv(var_vec,deriv_vec,cloud_tup) - c.g0/c.cpd*parcel.wvel
    deriv_vec[-2]= -1.*rho*c.g0*parcel.wvel
    deriv_vec[-1] = parcel.wvel
    return deriv_vec

def wlderiv(var_vec,deriv_vec,cloud_tup):
    """
    calculate the time derivative of the liquid water content
    
    Parameters
    ----------

    var_vec: vector(float)
        vector of values to be integrated

    deriv_vec: vector(float)
         derivatives of each of var_vec members

    cloud_tup: namedtuple
           tuple of input coefficients

    Returns
    -------

    dwldt: float
         rate of change of wl
    """
    wlderiv=(var_vec[:-3])**2.
    wlderiv=cloud_tup.ndist*wlderiv
    wlderiv=wlderiv*deriv_vec[:-3]
    dwldt = np.sum(wlderiv)*4.*np.pi*c.rhol
    return dwldt


if __name__ == "__main__":
    
    infile = 'dropgrow.yaml'
    with open(infile,'r') as f:
        input_dict=ry.load(f,Loader=ry.RoundTripLoader)
    #
    #set the edges of the mass bins
    #31 edges means we have 30 droplet bins
    #
    numrads = 30
    mass_vals = np.linspace(-20,-16,numrads+1)
    mass_vals = 10**mass_vals
    mu=input_dict['aerosol']['themean']
    sigma = input_dict['aerosol']['sd']
    totmass = input_dict['aerosol']['totmass']
    mdist = totmass*lognormal(mass_vals,np.log(mu),np.log(sigma))
    mdist = find_centers(mdist)*np.diff(mass_vals)
    center_mass = find_centers(mass_vals)
    ndist = mdist/center_mass

    cloud_vars = od()
    cloud_vars['mdist'] = mdist
    cloud_vars['ndist'] = ndist
    cloud_vars['center_mass'] = center_mass

    aero=make_tuple(input_dict['aerosol'])
    parcel=make_tuple(input_dict['initial_conditions'])

    a, b = find_koehler_coeffs(aero,parcel)

    #
    # sanity check
    #
    m=1.e-18
    Scrit=(4.*a**3./(27.*b*m))**0.5;
    rcrit = (3.*m*b/a)**0.5
    print("for aerosol with mass = {} kg, Scrit,rcrit are {}, {} microns".format(m,Scrit,rcrit))

    koehler_fun = create_koehler(aero,parcel)

    S_target = parcel.Sinit
    logr_start = np.log(0.1e-6)

    initial_radius = []
    dry_radius = []
    for mass in center_mass:
        brackets = np.array(find_interval(find_diff,logr_start,S_target,mass))
        left_bracket, right_bracket = np.exp(brackets)*1.e6  #get brackets in microns for printing
        equil_rad = np.exp(fzero(find_diff,brackets,S_target,mass))

        Scrit=(4.*a**3./(27.*b*mass))**0.5
        
        initial_radius.append(equil_rad)
        dry_rad = (mass/(4./3.*np.pi*aero.rhoaero))**(1./3.)
        dry_radius.append(dry_rad)

        print(('mass = {mass:6.3g} kg\n'
               'left bracket = {left_bracket:8.3e} microns\n'
               'right bracket={right_bracket:8.3e} microns\n'
               'critical supersaturation: {Scrit:6.3g}')
               .format_map(locals()))
        print('equlibrium radius at S={} is {:5.3f} microns\n'.format(S_target,equil_rad*1.e6))

    cloud_vars['initial_radiius'] = initial_radius
    cloud_vars['dry_radius'] = dry_radius
    cloud_vars['masses'] = center_mass
    numrads = len(initial_radius)
    var_vec = np.empty(numrads + 3)
    for i in range(numrads):
        var_vec[i] = initial_radius[i]
        
    var_vec[-3] = parcel.Tinit
    var_vec[-2] = parcel.Pinit
    var_vec[-1] = parcel.Zinit

    cloud_tup = make_tuple(cloud_vars)
    #calculate the total water (kg/kg)
    wl=wlcalc(var_vec,cloud_tup);
    e=parcel.Sinit*find_esat(parcel.Tinit);
    wv=c.eps*e/(parcel.Pinit - e);
    #save total water
    cloud_vars['wt'] = wv + wl;
    
    cloud_tup= make_tuple(cloud_vars)

    tinit=0
    # r = ode(find_derivs).set_integrator('dopri5')
    # r.set_f_params(cloud_tup)
    # r.set_initial_value(var_vec, tinit)

    var_out = []
    time_out =[]
    dt = 10
    tfin = 300
    t = np.arange(0,tfin,dt)
    sol = odeint(find_derivs,var_vec, t, args=(cloud_tup,))
    colnames = ["r{}".format(item) for item in range(30)]
    colnames.extend(['temp','press','z'])
    output = pd.DataFrame.from_records(sol,columns = colnames)

    if input_dict['dump_output']:
        with pd.HDFStore(input_dict['output_file'],'w') as store:
            store.put(input_dict['frame_name'],output,format='table')
            
    plt.close('all')
    fig, ax = plt.subplots(1,1)
    for i in colnames[:-3]:
        ax.plot(output[i]*1.e6,output['z']*1.e-3,label=i)
        #ax.plot(i,'z',data=output)
    plt.show()
    
