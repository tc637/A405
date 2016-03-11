"""
calculate the grow history of droplet categories in a constant velocity updraft
"""
import ruamel.yaml as ry
import a405dropgrow.aerolib
from importlib import reload
reload(a405dropgrow.aerolib)
from a405dropgrow.aerolib import lognormal,create_koehler
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
import ruamel.yaml, h5py, datetime

def find_diff(logr,S_target,m, koehler_fun):
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
    
    koehler_fun: function
       function that returns S as function of radius (m) and aerosol mass (kg)

    Returns
    -------

    Sdiff: float
       difference between target and guess
    """
    r = np.exp(logr)
    return S_target - koehler_fun(r,m)


def wlcalc(var_vec,cloud_tup):
    """
    calculate the liquid water by integrating n(r)r**3

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
    calculate the environmental saturation using conservation
    of total water mixing ratio cloud_top.wt and the current
    value of the liquid water mixing ratio wl

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
    #
    # find the evironmental S by water balance
    #
    S=Scalc(var_vec,cloud_tup)
    deriv_vec=np.zeros_like(var_vec)
    #dropgrow notes equaton 18 (W&H p. 170)
    for i in range(numrads):
        m=cloud_tup.masses[i]
        if var_vec[i] < dry_radius[i]:
            var_vec[i] = dry_radius[i]
        Seq=cloud_tup.koehler_fun(var_vec[i],m)  
        rhovr=(Seq*find_esat(temp))/(c.Rv*temp)
        rhovinf=S*find_esat(temp)/(c.Rv*temp)
        #day 25 drop_grow.pdf eqn. 18
        deriv_vec[i]=(c.D/(var_vec[i]*c.rhol))*(rhovinf - rhovr)
    #
    # moist adiabat day 25 equation 21a
    #
    deriv_vec[-3]=find_lv(temp)/c.cpd*wlderiv(var_vec,deriv_vec,cloud_tup) - c.g0/c.cpd*cloud_tup.wvel
    #
    # hydrostatic balance  dp/dt = -rho g dz/dt
    #
    deriv_vec[-2]= -1.*rho*c.g0*cloud_tup.wvel
    #
    # how far up have we traveled?
    #
    deriv_vec[-1] = cloud_tup.wvel
    return deriv_vec

def wlderiv(var_vec,deriv_vec,cloud_tup):
    """
    calculate the time derivative of the liquid water content
    using drop_grow.pdf eqn 21b
    
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

    from pathlib import Path
    util_dir, = a405utils.__path__._path
    data_dir = Path(util_dir).joinpath('../data')
    yaml_file = data_dir.joinpath('dropgrow.yaml')
    with yaml_file.open('r') as f:
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

    koehler_fun = create_koehler(aero,parcel)

    S_target = parcel.Sinit
    logr_start = np.log(0.1e-6)

    initial_radius = []
    dry_radius = []
    for mass in center_mass:
        brackets = np.array(find_interval(find_diff,logr_start,S_target,mass,koehler_fun))
        left_bracket, right_bracket = np.exp(brackets)*1.e6  #get brackets in microns for printing
        equil_rad = np.exp(fzero(find_diff,brackets,S_target,mass,koehler_fun))

        initial_radius.append(equil_rad)
        dry_rad = (mass/(4./3.*np.pi*aero.rhoaero))**(1./3.)
        dry_radius.append(dry_rad)

        print('mass = {mass:6.3g} kg\n'.format_map(locals()))
        print('equlibrium radius at S={} is {:5.3f} microns\n'.format(S_target,equil_rad*1.e6))

    cloud_vars['initial_radiius'] = initial_radius
    cloud_vars['dry_radius'] = dry_radius
    cloud_vars['masses'] = center_mass
    cloud_vars['koehler_fun'] = koehler_fun
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
    wv=c.eps*e/(parcel.Pinit - e)
    #save total water
    cloud_vars['wt'] = wv + wl
    cloud_vars['wvel'] = parcel.wvel
    
    cloud_tup= make_tuple(cloud_vars)

    var_out = []
    time_out =[]
    
    tinit=input_dict['integration']['dt']
    dt = input_dict['integration']['dt']
    tfin = input_dict['integration']['tend']
    
    t = np.arange(0,tfin,dt)
    sol = odeint(find_derivs,var_vec, t, args=(cloud_tup,))
    colnames = ["r{}".format(item) for item in range(30)]
    colnames.extend(['temp','press','z'])
    output = pd.DataFrame.from_records(sol,columns = colnames)

    plt.close('all')
    fig, ax = plt.subplots(1,1)
    for i in colnames[:-3]:
        ax.plot(output[i]*1.e6,output['z']*1.e-3,label=i)
        #ax.plot(i,'z',data=output)
    plt.show()
    
    if input_dict['dump_output']:
        with pd.HDFStore(input_dict['output_file'],'w') as store:
            store.put(input_dict['frame_name'],output,format='table')
            
        yaml_string = ruamel.yaml.dump(input_dict,None,Dumper=ruamel.yaml.RoundTripDumper,
                               default_flow_style=False)

        with h5py.File(input_dict['output_file'],'a') as f:
            f.attrs['yaml_string']=yaml_string
            date=datetime.datetime.now().strftime('%Y-%M-%d')
            history ="file produced by drop_grow.py on {}".format(date)
            print('history: ',history)
            f.attrs['history']=history
