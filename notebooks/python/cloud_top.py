
# coding: utf-8

# # Calculate CAPE for a Littlerock sounding

# In[1]:

import numpy as np
import pandas as pd
import h5py
from pprint import pformat
from importlib import reload
import a405thermo.thermlib
reload(a405thermo.thermlib)
from a405thermo.constants import constants as c
from a405thermo.thermlib import convertSkewToTemp, convertTempToSkew
import a405skewT.makeSkewII
reload(a405skewT.makeSkewII)
from a405skewT.makeSkewII import makeSkewWet


# In[2]:

filename='single_littlerock.h5';
print('reading file: %s\n' %filename)

attributes={}
with h5py.File(filename,'r') as f:
    keys=f.attrs.keys()
    for key in keys:
        try:
            attributes[key]=f.attrs[key]
        except IOError:
            print('empty key: ',key)
print('\nread in these attributes: \n\n',pformat(attributes))

separator= '\n' + '+'*30 + '\n'
sounding_dict={}
with pd.HDFStore(filename,'r') as store:
    times=store.keys()
    for the_time in times:
        sounding_dict[the_time]=store[the_time]
sounding=sounding_dict[times[3]]
# print('{0:}sounding times:{1:}{0:} '.format(separator,times))
# print('{0:}columns: {1:}{0:} '.format(separator,pformat(sounding.columns)))
title_string=attributes['header']
index=title_string.find(' Observations at')
location=title_string[:index]
title='{} at {}'.format(location,times[0][2:])
print('title: :',title)
units=attributes['units'].split(';')
units_dict={}
for count,var in enumerate(sounding.columns):
    units_dict[var]=units[count]


# In[3]:

skew=30.
triplets=zip(sounding['temp'],sounding['dwpt'],sounding['pres'])
xcoord_T=[]
xcoord_Td=[]
for a_temp,a_dew,a_pres in triplets:
    xcoord_T.append(convertTempToSkew(a_temp,a_pres,skew))
    xcoord_Td.append(convertTempToSkew(a_dew,a_pres,skew))


# In[4]:

get_ipython().magic('matplotlib inline')
skew=30
plt.close('all')
fig,ax =plt.subplots(1,1,figsize=(8,8))
corners=[-15,35]
ax,skew = makeSkewWet(ax,corners=corners,skew=skew)
ax.plot(xcoord_T,sounding['pres'],color='k',label='temp')
ax.plot(xcoord_Td,sounding['pres'],color='g',label='dew')
[line.set(linewidth=3) for line in ax.lines[-2:]]
out=ax.set(title=title)


# In[5]:

from a405thermo.thermlib import find_Tmoist,find_thetaep,find_rsat,find_Tv
#
# find thetae of the surface air
#
sfc_press,sfc_temp,sfc_td =[sounding[key][0] for key in ['pres','temp','dwpt']]
sfc_press,sfc_temp,sfc_td = sfc_press*100.,sfc_temp+c.Tc,sfc_td+c.Tc
sfc_rvap = find_rsat(sfc_temp,sfc_press)
sfc_thetae=find_thetaep(sfc_td,sfc_temp,sfc_press)
press=sounding['pres'].values*100.
#
# find the index for 200 hPa pressure -- searchsorted requires
# the pressure array to be increasing, so flip it for the search,
# then flip the index.  Above 100 hPa thetae goes bananas, so
# so trim so we only have 
#
toplim=len(press) - np.searchsorted(press[::-1],1.5e4)
clipped_press=press[:toplim]
#
# find temps along that adiabat
#
adia_temps= np.array([find_Tmoist(sfc_thetae,the_press) for the_press in clipped_press])
adia_rvaps = find_rsat(adia_temps,clipped_press)
adia_rls = sfc_rvap - adia_rvaps
env_temps = (sounding['temp'].values + c.Tc)[:toplim]
env_Td = (sounding['dwpt'].values + c.Tc)[:toplim]
height = sounding['hght'].values[:toplim]
pairs = zip(env_Td,clipped_press)
env_rvaps= np.array([find_rsat(td,the_press) for td,the_press in pairs])
env_Tv = find_Tv(env_temps,env_rvaps)
adia_Tv = find_Tv(adia_temps,adia_rvaps,adia_rls)
xcoord_thetae=[]
clipped_press_hPa = clipped_press*1.e-2
for a_temp,a_press in zip(adia_temps - c.Tc,clipped_press_hPa):
    out=convertTempToSkew(a_temp,a_press,skew)
    xcoord_thetae.append(out)
ax.plot(xcoord_thetae,clipped_press_hPa,color='r',label='rsat',linewidth=3.)
ax.set(ylim=[1000.,200.])
display(fig)


# In[6]:

def find_buoy(adia_Tv,env_Tv):
  buoy=c.g0*(adia_Tv - env_Tv)/env_Tv
  return buoy

#
# moved find_buoy to library
#
from a405thermo.thermlib import find_buoy
    
fig,ax =plt.subplots(1,1)
buoy=find_buoy(adia_Tv,env_Tv)
ax.plot(buoy,clipped_press*1.e-2)
ax.invert_yaxis()
out=ax.set(ylabel='press (hPa)',xlabel='buoyancy (m/s^2)')


# ### Find the bottom and top buoyancy zero crossings

# In[7]:

#np.searchsorted finds the first crossing
zerobot=np.searchsorted(buoy,0)
#flip the array and search backwards for top crossing
zerotop=len(buoy) - np.searchsorted(buoy[::-1],0)
print('pressure levels for crossing: {} hPa {} hPa'      .format(press[zerobot]*0.01,press[zerotop]*0.01))


# ### find the cape

# In[8]:

clipped_buoy = buoy[zerobot:zerotop]
clipped_height = height[zerobot:zerotop]
#average the levels to get layer buoyancy
layer_buoy = (clipped_buoy[:-1] + clipped_buoy[1:])/2.
cape = np.sum(layer_buoy*np.diff(clipped_height))
print('cape is {:6.2f} J/kg'.format(cape))


# ### Find the cumulative buoyancy

# In[9]:

clipped_buoy = buoy[zerobot:]
clipped_height = height[zerobot:]
new_press= clipped_press[zerobot:]
#average the levels to get layer buoyancy
layer_buoy = (clipped_buoy[:-1] + clipped_buoy[1:])/2.
layer_height = (clipped_height[:-1] + clipped_height[:-1])/2.
layer_press = (new_press[:-1] + new_press[1:])/2.
cape = np.cumsum(layer_buoy*np.diff(clipped_height))
fig,ax=plt.subplots(1,1)
ax.plot(cape,layer_press*1.e-2)
ax.invert_yaxis()


# ### Create an interpolator for the environmental Tv using interp1d

# In[10]:

from scipy.interpolate import interp1d
env_Tv_interp=interp1d(clipped_press,env_Tv)
fig,ax = plt.subplots(1,1)
ax.plot(env_Tv,clipped_press*0.01,label='orig sounding')
ax.invert_yaxis()
reg_press=np.linspace(9.8e4,6.e4,25)
ax.plot(env_Tv_interp(reg_press),reg_press*0.01,'ro',label='interp')
ax.set(ylim=(1000,650))
out=ax.legend(loc='best')


# In[11]:

len(env_temps),len(env_rvaps)


# In[12]:

from a405thermo.constants import constants as c


# In[13]:

dir(c)


# In[ ]:



