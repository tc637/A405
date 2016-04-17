
# coding: utf-8

# ## Working with a 3D large eddy simulation of shallow convection

# ### The simulation
# 
# * Objective: compare a single column of a GCM with large eddy simlations for three different cloud types (stratus, stratocumulus, trade cumulus)
# 
# [GCM paper](http://ezproxy.library.ubc.ca/login?url=http://doi.wiley.com/10.1002/jame.v5.4)
# 
# [LES paper](http://ezproxy.library.ubc.ca/login?url=http://doi.wiley.com/10.1002/jame.20025)
# 
# * We started with the trade cumulus simulation, then perturbed it by raising the temperature to 300 K and 301 K.
# 
# http://clouds.eos.ubc.ca/~phil/courses/atsc405/docs/cgils_ctl_s6_synthetic_albedo.mp4
# 
# http://clouds.eos.ubc.ca/~phil/courses/atsc405/docs/cgils_sst_300K_synthetic_albedo.mp4
# 
# http://clouds.eos.ubc.ca/~phil/courses/atsc405/docs/cgils_sst_301K_synthetic_albedo.mp4
# 

# ###  The dataset  -- netccdf
# 
# [An example of reading a netCDF4 file ](http://schubert.atmos.colostate.edu/~cslocum/netcdf_example.html)
# 

# In[2]:

import glob
from netCDF4 import Dataset
from a405utils.ncdump import ncdump
from a405utils.download import download_file

download = False
if download:
    #
    #  satelite data for day 127 of 2014  Modis Aqua level 3 cloud data
    #
    url = 'http://clouds.eos.ubc.ca/~phil/Downloads/a405/ENT_CGILS_CTL_S6_3D_384x384x194_25m_1s_96_0000014160.nc'
    local_file = download_file(url)
    print('downloaded {}'.format(local_file))
    
the_file = glob.glob("*CTL*")[0]
with Dataset(the_file,'r') as ncin:
    ncdump(ncin)


# ### liquid water cross section at 1 km
# 
# 

# In[3]:

def get_var(the_file,varname):
    with Dataset(the_file,'r') as ncin:
         out=ncin.variables[varname][...]
         x = ncin.variables['x'][...]
         y = ncin.variables['y'][...]
         z = ncin.variables['z'][...]
         out = out.squeeze()  #remove the time dimension, since we only have one timestep
    return x,y,z,out
x,y,z,qn = get_var(the_file, 'QN')
print(the_file)
print(qn.max())


# In[4]:

#
#  find the index for z = 1000 meters
#

level = np.searchsorted(z, 1000)
print(level)


# In[5]:

get_ipython().magic('matplotlib inline')

# get the cloud liquid water at 1000 m
#
horiz_cross_sec = qn[level,:,:]
#
# find the cross section cloud fraction
#
cloud_frac=np.sum(horiz_cross_sec > 0)/horiz_cross_sec.size
print('cloud fraction: {:5.3f}'.format(cloud_frac))


# In[6]:

get_ipython().magic('matplotlib inline')
plt.style.use('ggplot')
from matplotlib import pyplot as plt
plt.close('all')
fig,ax =plt.subplots(1,1,figsize=(10,10))
whole_scene=ax.imshow(horiz_cross_sec)
cax=plt.colorbar(whole_scene,ax=ax)
cax.set_label('liquid water content (g/kg)')
title = 'horizontal qn cross section at z=1000 m'
ax.set_title(title)


# ### zoom in on  the top left corner
# 
# Switch from [imshow](http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.imshow) to 
# [pcolormesh](http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.pcolormesh) so we can orient the axes along model x,y, and z coordinates.  Note that if y is north/south (north up), then imshow plots the image upside down.

# In[7]:

#
# it helps in checking your orientation to make the selection
# have different numbers of rows and columns
#
end_col = 200
end_row = 180
fig,ax =plt.subplots(1,1,figsize=(10,10))
image=ax.pcolormesh(x[:end_col],y[:end_row],horiz_cross_sec[:end_row,:end_col])
ax.set(xlabel='distance east',ylabel='distance north')
cax = plt.colorbar(image,ax=ax)
cax.set_label('liquid water content (g/kg)')
ax.set_title('zoomed horiz qn cross section at z=1000 m')


# ### Get a vertical cross section along y = 2km

# In[8]:

row_number = np.searchsorted(y,2000)  #(y index of 80)


# In[9]:

row_number


# In[10]:

vert_cross_sec = qn[:,row_number,:end_col]
print(vert_cross_sec.shape)
print(len(z),len(x[:end_col]))
fig,ax = plt.subplots(1,1,figsize=(10,10))
image=ax.pcolormesh(x[:end_col],z,vert_cross_sec[:,:end_col])
cax = plt.colorbar(image,ax=ax)
cax.set_label('liquid water mixing ratio qn (g/kg)')
ax.set_title('vertical qn cross section along y=2 km')
ax.set(xlabel='distance east (m)',ylabel='height (m)')


# ### Find the vapor mixing ratio along this cross section

# In[11]:

plt.close('all')
x,y,z,qv = get_var(the_file, 'QV')
vert_cross_sec = qv[:,row_number,:end_col]
fig,ax = plt.subplots(1,1,figsize=(10,10))
image=ax.pcolormesh(x[:end_col],z,vert_cross_sec[:,:end_col])
cax = plt.colorbar(image,ax=ax)
cax.set_label('water vapor mixing ratio qv (g/kg)')
ax.set_title('vertical qv cross section along y=2 km')


# ### For Monday
# 
# 1\.  Read Thompkins chapter 4 parameterization notes through section 4.7.1
# 
# 2\.  Read [Zhu and Randall, 1996](http://kiwi.atmos.colostate.edu/pubs/XuandRandall-semiempirical-1996.pdf)
# 
# 3\.  Hand in a notebook that adds cells to cgilsI.ipynb to:
# 
#      * use pcolormesh to plot a vertical cross section of the relative humidity for along y=2 km, x= 0-5 km
#      
#      * use plot to plot a vertical profile of the horizontal mean RH in for this cross section as a function of height
#      
#      * use plot to plot a vertical profile of the horizontal standard deviation of RH as a function of height

# ### Calculate the relative humidity from qv and rsat
# 
# Use newaxis to make the shapes work for [numpy broadcasting](http://www.scipy-lectures.org/intro/numpy/operations.html#broadcasting)

# In[12]:

from a405thermo.thermlib import find_rsat, find_esat
from a405thermo.constants import constants as c
x,y,z,temp = get_var(the_file,'TABS')
x,y,z,press = get_var(the_file,'p')
vert_cross_sec = temp[:,row_number,:end_col]
plt.close('all')
fig,ax = plt.subplots(1,1,figsize=(10,10))
image=ax.pcolormesh(x[:end_col],z,vert_cross_sec[:,:end_col])
cax = plt.colorbar(image,ax=ax)
cax.set_label('temperature (K)')
ax.set_title('vertical temp cross section along y=2 km')

esat = find_esat(temp)
print('esat is 3D: ',esat.shape)
print('press is 1D: ',press.shape)
#
# add dimensions to the 1D pressure vectory so
# that it can be used in the denominator 
#
press = press[:,np.newaxis,np.newaxis]*100.  #convert to Pa, make pressure 3D
rsat = c.eps*esat/(press - esat)*1.e3  #convert to g/kg, rsat has dimensions of esat
rh = qv/rsat

fig,ax = plt.subplots(1,1,figsize=(10,10))
vert_cross_sec = rh[:,row_number,:end_col]
image=ax.pcolormesh(x[:end_col],z,vert_cross_sec[:,:end_col])
cax = plt.colorbar(image,ax=ax)
cax.set_label('relative humidity')
ax.set_title('relative humidity cross section along y=2 km')


# ### screen out cells with RH > 1 using a masked array
# 
# Also, contract the range of the colormap so that all the 256 colors
# are spread between RH 0.8 to 0.98

# In[13]:

from matplotlib.colors import Normalize
fig,ax=plt.subplots(1,1,figsize=(10,10))
pal = plt.get_cmap('plasma')
pal.set_bad('0.75') #75% grey
pal.set_over('r')  #color cells > 0.98 red
pal.set_under('k')  #color cells < 0.75 black
vmin= 0.8
vmax= 0.98
#
#mask relative humidities > 1
#using a masked array
#
import numpy.ma as ma
mask = vert_cross_sec > 1
ma_rh_vert = ma.array(vert_cross_sec, mask = mask)
the_norm=Normalize(vmin=vmin,vmax=vmax,clip=False)
image=ax.pcolormesh(x[:end_col],z,ma_rh_vert[:,:end_col],cmap=pal,norm=the_norm)
cax=plt.colorbar(image,ax = ax,extend='both')
ax.set(ylim=[0,2500])



    



# ### find the mean and variance for the subset

# In[14]:

end_z = np.searchsorted(z,2500.)
rh_subset = rh[:end_z,:end_row,:end_col]
mean_rh_xy = rh_subset.mean(axis=(1,2))
sd_rh_xy = rh_subset.var(axis=(1,2))**0.5
fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10,10))
ax1.plot(mean_rh_xy,z[:end_z])
ax2.plot(sd_rh_xy*100.,z[:end_z])
ax1.set(ylabel='height (m)',xlabel='RH')
out=ax2.set(xlabel='RH standard dev x 100')


# ### RH histogram outside of cloud

# In[15]:

rh_boundary = rh[:end_z,:,:]
rh_environ = rh_boundary[rh_boundary < 1.]
fig,ax = plt.subplots(1,1,figsize=(10,10))
out=ax.hist(rh_environ.flatten(),bins=50)
ax.set_title('histogram of RH for boundary layer environment')


# ### Problem for Wednesday, April 6
# * check in a notebook to
# 
#     - find theate and qT (qv + qn) for our LES volume below z = 2500 m.
#     
#     - plot histograms of thetae and qT, and a scatterplot that looks like
#       Fig. 1 of Zhu and Randall II, 1996 (you don' have to put the
#       saturation line on the figure, but you should be able to tell me how to
#       use a rootfinder to find that line at a particular pressure.
# 

# In[16]:

import a405thermo.thermlib
reload(a405thermo.thermlib)
from a405thermo.thermlib import find_Td,find_thetaet, find_thetal
from a405utils.progressbar import ProgressBar
import numpy.random as nr
plt.close('all')


# ### 1\.  plot the average soundings

# In[26]:

end_z = np.searchsorted(z,2600.)
bot_z = 0
# bot_z = np.searchsorted(z,500.)

end_z = np.searchsorted(z,1200.)
bot_z = np.searchsorted(z,1000.)


qv_sound = qv.mean(axis=(1,2))[bot_z:end_z]*1.e-3 #kg/kg
temp_sound=temp.mean(axis=(1,2))[bot_z:end_z]
press_sound=np.squeeze(press)[bot_z:end_z]
z_sound = z[bot_z:end_z]
#print(qv_sound.shape,temp_sound.shape,z_sound.shape,press_sound.shape)
thetal_sound = [find_thetal(the_press,the_temp,the_qv) for 
                the_press,the_temp,the_qv in zip(press_sound,temp_sound,qv_sound)]
Td_sound = find_Td(qv_sound,press_sound)
thetae_sound = [find_thetaet(the_Td,the_qv,the_temp,the_press) for
               the_Td,the_press,the_temp,the_qv in zip(Td_sound,press_sound,temp_sound,qv_sound)]
fig,(ax1,ax2,ax3,ax4) = plt.subplots(1,4,figsize=(16,8))
ax1.plot(qv_sound*1.e3,z_sound)
ax2.plot(temp_sound,press_sound*1.e-2)
ax3.plot(thetal_sound,press_sound*1.e-2)
ax4.plot(thetae_sound,press_sound*1.e-2)
axes = [ax2,ax3,ax4]
[ax.invert_yaxis() for ax in axes]
p0=press_sound[0]*1.e-2
[ax.set(ylim=(p0,750)) for ax in axes]
axes = [ax1, ax2,ax3,ax4]
ax1.set(ylabel='heigth (m)')
ax2.set(ylabel='pressure (hPa)')
titles = ['$q_v\ (g/kg)$', '$temperature\ (K)$', r'$\theta_l\ (K)$',r'$\theta_e\ (K)$']
out=[ax.set_title(the_title) for ax,the_title in zip(axes,titles)]


# ### 2\.  select a subregion from bot_z to top_z

# In[27]:

qv_bl,qn_bl,temp_bl,press_bl = qv[bot_z:end_z,...],qn[bot_z:end_z,...],temp[bot_z:end_z,...],press[bot_z:end_z]


# In[28]:

#
# note that these are views into the original data, not copies
# so changing the sliced regions will change the original data
#
qv_bl.flags.owndata


# ### 3\. Calculate thetae and thetal for a random sample in the subregion
# 
# Make a scatterplot with cloud pixels in blue and clear in red
# Add the constant rsat(T,p) line for press=900 hPa
# 
# Include a progressbar for the impatient

# In[29]:



Td_bl = find_Td(qv_bl*1.e-3,press_bl)
qt_bl = (qv_bl + qn_bl)*1.e-3

nr.seed(10)
new_shape = [20,40,50]
old_shape = temp_bl.shape  #100 x 384 x 384
thetae_sample=np.empty(new_shape) 
thetal_sample = np.empty(new_shape)
qn_sample = np.empty(new_shape)
levels, rows, cols =[nr.randint(the_shape,size=the_size) 
                         for the_shape,the_size in zip(old_shape,new_shape)]
qt_sample =np.empty(new_shape)
pbar = ProgressBar(thetae_sample.size)
new_lev=0
for old_lev in levels:
    new_row = 0
    for old_row in rows:
        new_col = 0
        for old_col in cols:
            the_Td = Td_bl[old_lev,old_row,old_col]
            the_qt = qt_bl[old_lev,old_row,old_col]
            the_press = press_bl[old_lev]
            the_temp = temp_bl[old_lev,old_row,old_col]
            thetae_sample[new_lev,new_row,new_col] = find_thetaet(the_Td,the_qt,the_temp,the_press)
            thetal_sample[new_lev,new_row,new_col] = find_thetal(the_press,the_temp,the_qt)
            qt_sample[new_lev,new_row,new_col] = the_qt
            qn_sample[new_lev,new_row,new_col] = qn_bl[old_lev,old_row,old_col]
            new_col += 1
            pbar.increment()
        new_row += 1
    new_lev += 1
pbar.finish()


# In[31]:

cloud_points = qn_sample > 0.
clear_points = np.logical_not(cloud_points)
print(np.sum(cloud_points),np.sum(clear_points))
thetal_clear=thetal_sample[clear_points]
thetal_cloud=thetal_sample[cloud_points]
thetae_clear=thetae_sample[clear_points]
thetae_cloud=thetae_sample[cloud_points]
qt_cloud=qt_sample[cloud_points]
qt_clear=qt_sample[clear_points]
plt.close('all')
fig,(ax1,ax2) = plt.subplots(1,2,figsize=(12,8))
print(qt_cloud.size,thetae_cloud.size)
ax1.plot(qt_cloud.ravel()*1.e3,thetae_cloud.ravel(),'b+',label='cloud')
ax1.plot(qt_clear.ravel()*1.e3,thetae_clear.ravel(),'r+',label='clear')
press_lev=9.e4
qv_vals=np.linspace(7,15,30)*1.e-3
td_vals = find_Td(qv_vals,press_lev)
#
#  loop over qv and find the thetae and thetal for a saturated mixture at
#  that qv and press= 900 hPa
#
thetae_line = [find_thetaet(the_td,the_qv,the_td,press_lev) for the_td,the_qv in zip(td_vals,qv_vals)]
thetal_line = [find_thetal(press_lev,the_td,the_qv) for the_td,the_qv in zip(td_vals,qv_vals) ]
ax1.plot(qv_vals*1.e3,thetae_line,'g-',label='sat at 900hPa')
ax1.set(ylim=[325,340])
ax1.set(xlabel='qt (g/kg)',ylabel='thetae (K)', )
ax2.set(xlabel='qt (g/kg)',ylabel='thetal (K)', )
ax2.plot(qt_cloud.ravel()*1.e3,thetal_cloud.ravel(),'b+',label='cloud')
ax2.plot(qt_clear.ravel()*1.e3,thetal_clear.ravel(),'r+',label='clear')
ax2.plot(qv_vals*1.e3,thetal_line,'g-',label='sat at 900hPa')
ax1.legend(loc='best')
ax2.legend(loc='best')


# In[ ]:



