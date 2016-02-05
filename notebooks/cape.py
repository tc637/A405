
# coding: utf-8

# # Add a constant $r_s$ line to the Skew-T - ln P plot
# 
# Demonstrate how to construct dry adiabats and isotherms for
# a thermodynamic diagram using the functions in
# [makeSkew.py](https://github.com/phaustin/A405/blob/master/a405skewT/makeSkew.py)

# In[9]:

import numpy as np
import pandas as pd
import h5py
from pprint import pformat
from a405thermo.constants import constants as c
from a405skewT.makeSkew import makeSkewDry
from importlib import reload
import a405thermo.thermlib
reload(a405thermo.thermlib)
from a405thermo.thermlib import convertSkewToTemp, convertTempToSkew
from a405thermo.thermlib import find_Tmoist,find_thetaep,find_rsat,find_Tv


# In[6]:

filename='littlerock.h5';
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


# In[22]:

try:
    del fig
except NameError:
    pass
fig,ax =plt.subplots(1,1,figsize=(8,8))
ax,skew = makeSkewDry(ax)
triplets=zip(sounding['temp'],sounding['dwpt'],sounding['pres'])
xcoord_T=[]
xcoord_Td=[]
for a_temp,a_dew,a_pres in triplets:
    xcoord_T.append(convertTempToSkew(a_temp,a_pres,skew))
    xcoord_Td.append(convertTempToSkew(a_dew,a_pres,skew))
ax.plot(xcoord_T,sounding['pres'],color='k',label='temp')
ax.plot(xcoord_Td,sounding['pres'],color='g',label='dew')

#
#  make the last two lines added to the plot (temp and dewpoint)
#  thicker
#
[line.set(linewidth=3) for line in ax.lines[-2:]]
out=ax.set(title=title)


# ## Add a moist adiabat

# In[23]:

import a405skewT.makeSkewII
reload(a405skewT.makeSkewII)
from a405skewT.makeSkewII import makeSkewWet
plt.close('all')
fig,ax =plt.subplots(1,1,figsize=(8,8))
corners=[-15,35]
ax,skew = makeSkewWet(ax,corners=corners,skew=skew)
ax.plot(xcoord_T,sounding['pres'],color='k',label='temp')
ax.plot(xcoord_Td,sounding['pres'],color='g',label='dew')
[line.set(linewidth=3) for line in ax.lines[-2:]]
out=ax.set(title=title)
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
# then flip the index
#
toplim=len(press) - np.searchsorted(press[::-1],2.e4)
press=press[:toplim]
adia_temps= np.array([find_Tmoist(sfc_thetae,the_press) for the_press in press])
adia_rvaps = find_rsat(adia_temps,press)
adia_rls = sfc_rvap - adia_rvaps
env_temps = (sounding['temp'].values + c.Tc)[:toplim]
env_Td = (sounding['dwpt'].values + c.Tc)[:toplim]
pairs = zip(env_Td,press)
env_rvaps= np.array([find_rsat(td,the_press) for td,the_press in pairs])
env_Tv = find_Tv(env_temps,env_rvaps)
adia_Tv = find_Tv(adia_temps,adia_rvaps,adia_rls)
xcoord_thetae=[]
press_hPa = press*1.e-2
for a_temp,a_press in zip(adia_temps - c.Tc,press_hPa):
    out=convertTempToSkew(a_temp,a_press,skew)
    xcoord_thetae.append(out)
ax.plot(xcoord_thetae,press_hPa,color='r',label='rsat',linewidth=3.)
ax.set(ylim=[1000.,200.])
display(fig)                                   

fig,ax =plt.subplots(1,1)
ax.plot(adia_Tv - env_Tv,press)
ax.invert_yaxis()
display(fig)


# In[ ]:



