
# coding: utf-8

# 1.  Buoyancy and stability
# Questions
# 
# i)  (3pt) Our by-now familiar air parcel of moist (but non-cloudy)
# 
# air has a temperature of 24.9 deg C. It is neutrally buoyant sitting in an
# environment that has a temperature of 25.0 oC and a mixing ratio of 10
# g kg−1. What is the mixing ratio in the air parcel?

# In[28]:

from importlib import reload
import a405thermo.constants
reload(a405thermo.constants)
from a405thermo.constants import constants as c
# c.delta = 0.608
rvenv = 1.e-2
Tenv = c.Tc + 25
Tparcel = c.Tc + 24.9
Tvenv = Tenv*(1. + c.delta*rvenv)
#Tvparcel = Tparcel*(1 + c.delta*rvparcel)
#want Tvparcel = Tvenv
rvparcel =(Tvenv/Tparcel - 1)/c.delta
print('rvparcel: {:5.3f} g/kg'.format(rvparcel*1.e3))


# ii (2pt) the temperature of the air parcel is increased by 1 oC, so
# that it is no longer neutrally buoyant. Assuming the resulting
# acceleration is maintained (i.e. constant) for 10 seconds, what is the
# final parcel velocity?

# In[30]:

Tparcel = Tparcel + 1
Tvparcel = Tparcel*(1 + c.delta*rvparcel)
buoy = c.g0*(Tvparcel - Tvenv)/Tvenv
print('parcel vel at 10 seconds={:5.3f} m/s'.format(buoy*10))


# iii (2 pt) Dry air at 1000 hPa is measured to have a temperature of
# 27oC, while at 900 hPa the temperature is 22oC, is this layer of the
# atmosphere dry stable, neutral or unstable?

# In[33]:

temp1 = c.Tc + 27
press1 = 1.e3
p0 = 1.e3
theta1 = temp1*(p0/press1)**(c.Rd/c.cpd)
temp2 = c.Tc + 22
press2 = 9.e2
theta2 = temp2*(p0/press2)**(c.Rd/c.cpd)
print('theta1, theta2 stable',theta1,theta2)


# 2. Fog
# 
# In this question, for simplicity, assume the definition of relative humidity 
# defined in terms of mixing ratio: RH = rv
# 
# i (2pt) The environmental air with temperature of 25.0 oC and a mixing
# ratio of 10 g kg−1 is located near the surface at a pressure of
# 1000hPa. What is the relative humidity?
# 
#  ii (1pt) The airstarts to cool isobarically at nightnight until even-
#  tually a fog forms. Assuming the pressure is 1000hPa. Which is the
#  temperature at which the fogs occurs? (Use the tephigram to answer
#  this, please mark the answer on the tephigram)
# 
# iii (2pt) Derive this same temperature by inverting Teton’s formula,
# showing your working. (Out of interest note how close the answer is to
# your tephigram result).
# 
# iv (1pt) What is the common name given to this temperature?

# In[46]:

get_ipython().magic('matplotlib inline')
import a405skewT.makeSkewII
reload(a405skewT.makeSkewII)
from a405skewT.makeSkewII import makeSkewWet
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
ax, skew = makeSkewWet(ax,corners=[10,30])
ax.set(ylim=[1000,700])
#invert Teton
c1 = 273.16
c2=32.19
press=1.e5
rsat=1.e-2
c3=np.log(press*rsat/380.)/17.5
Tdew = (c3*c2 - c1)/(c3 - 1)
print("dewpoint is: {:5.3f} deg C".format(Tdew - c.Tc))


# 3.  Relative humidity
# 
# In this question, for simplicity, assume the definition of relative
# hu- midity defined in terms of mixing ratio: RH = rv .  s
# 
# i (2pt) A parcel of air at time t=0 has a temperature T0 of 25oC and a
# relative humidity RH0 of 0.1. What is the mixing ratio rv0?
# 
# ii (4pt) The parcel is cooled and moistened isobarically by the evap-
# oration of precipitation to reach a final temperature of tempera- ture
# of T1 and a final mixing ratio of rv1. If we ignore the humid-
# ity/precipitation in the parcel’s heat capacity (i.e. heat capacity
# is cp as for dry air) and treat the latent heat of vaporization as
# a constant Lv that is independent of temperature, then we can
# relate the change in temperature to the change in mixing ratio as
# follows:
# 
# dT
# calculated at temperature T0.)
# iii (1pt) If RH1 were instead 1.0, what is the temperature T1 com-
# monly known as?

# In this question we use the attached tephigram with the thermody-
# namic profiles
# 4.i (1pt) What is the lifting condensation level (LCL) of air originat-
# ing from 700 hPa?
# 4.ii (2pt) Does the wet bulb potential temperature $\theta_w$ increase or
# decrease between 700 and 600 hPa?
# 4.iii (1 pt) Given the answer to the previous question, you can now
# state whether the equivalent potential temperature is higher at
# 600 or 700 hPa, why?
# 4.iv (3 pts) Air at 700 hPa is brought to saturation by the evaporation
# of rainfall. It then decends to 1000 hPa, and during this descent
# is subject to further rainfall evaporation such that it arrives at
# 1000hPa still saturated. What is the nal temperature?

# In[47]:

fig, ax = plt.subplots(1, 1, figsize=(10, 10))
ax, skew = makeSkewWet(ax,corners=[5,30])


# In[ ]:



