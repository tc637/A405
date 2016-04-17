
# coding: utf-8

# ## The cell below shows how to use a rootfinder to solve cos(x) = 0.75

# In[3]:

get_ipython().magic('matplotlib inline')
from scipy import optimize
from matplotlib import pyplot as plt
import numpy as np
plt.style.use('ggplot')


xvals=np.linspace(0,10.)
fig,ax = plt.subplots(1,1,figsize=(10,8))
ax.plot(xvals,np.cos(xvals))
straight_line=np.ones_like(xvals)
ax.plot(xvals,straight_line*0.75,'b-')

def root_function(x):
    """Function we want to find the root of
       input: x value
       output: y value -- should be zero when x is a root
    """
    return np.cos(x) - 0.75

root1 = optimize.zeros.brentq(root_function,0,2)
root2 = optimize.zeros.brentq(root_function,4,6)
root3 = optimize.zeros.brentq(root_function,6,8)
xvals=np.array([root1,root2,root3])
yvals=np.cos(xvals)
ax.plot(xvals,yvals,'ro',markersize=20)

print(optimize.zeros.__file__)


# ## Problem for Wednesday:  9am
# 
# Write a function:
#     
#     def temp_from_theta(theta,press)
#     
# that uses brentq to find the temperature for a theta of 280 K and pressures of press=[1.e5, 7.e4, 4.e4] Pascals
# 
# Email a working notebook as an attachment, or send a url to a dropbox-like repository

# ### Solution

# In[4]:

def root_function(Tguess,theta,press):
    """Function we want to find the root of
       input: Tguess (K), target theta (K), press (Pa)
       output: difference between guess and target -- should be zero when x is a root
    """
    p0=1.e5  #Pa
    Rd=287. #J/kg/K
    cp=1004.  #J/kg/K     
    theta_guess=Tguess*(p0/press)**(Rd/cp)
    return theta - theta_guess


# In[5]:

def temp_from_theta(theta,press):
    """
       input: theta (K), press (Pa)
       output: temp (K) found by rootfinder
    """     
    left=10 #K
    right=1000 #K
    temp = optimize.zeros.brentq(root_function,left,right,args=(theta,press))
    return temp

for press in [1.e5,7.e4,4.e4]:
    print('Temp {:5.2f} (K) at pressure of {:5.2f} kPa'.format(temp_from_theta(280.,press),press*1e-2))


# ## Bracket finding
# 
# I've written a couple of convenience functions called rootfinder.find_interval and
# rootfinder.fzero to make rootfinding a little easier.   The new module is 
# [rootfinder.py](https://github.com/phaustin/A405/blob/master/a405thermo/rootfinder.py)
# 

# In[6]:

from a405thermo import rootfinder as rf


# In[7]:

print(help(rf.find_interval))


# ### example -- find a bracket for sin(x)=0 near x=12 radians ~ 700 degrees

# In[8]:

brackets=rf.find_interval(np.sin,12)
brackets


# ## now use the fzero wrapper to find the root of sin(x)=0  (720 degrees)

# In[9]:

print(rf.fzero(np.sin,brackets)*180/np.pi)


# ## Redo theta example with find_interval

# In[39]:

import a405thermo.rootfinder as rf
from importlib import reload
reload(rf)

def theta_zero(Tguess,theta,press):
    """Function we want to find the root of
       input: Tguess (K), target theta (K), press (Pa)
       output: difference between guess and target -- should be zero when x is a root
    """
    p0=1.e5  #Pa
    Rd=287. #J/kg/K
    cp=1004.  #J/kg/K     
    theta_guess=Tguess*(p0/press)**(Rd/cp)
    return theta - theta_guess


# In[43]:

def temp_from_theta(theta,press):
    """
       input: theta (K), press (Pa)
       output: temp (K) found by rootfinder
    """     
    #
    #  use theta as guess for bracket and pass theta,press to theta_zero
    #
    brackets=rf.find_interval(theta_zero,theta,theta,press)
    the_temp = rf.fzero(theta_zero,brackets,theta,press)
    return the_temp

for press in [1.e5,7.e4,4.e4]:
    print('Temp {:5.2f} (K) at pressure of {:5.2f} kPa'.format(temp_from_theta(280.,press),press*1e-2))


# In[ ]:




# In[ ]:



