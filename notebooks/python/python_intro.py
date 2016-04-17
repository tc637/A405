
# coding: utf-8

# ## Matlab vs. Python
# 
# This notebook uses the python matlab bridge to do some side by side examples for python and matlab

# In[4]:

import pymatbridge as pymat


# In[5]:

ip = get_ipython()
pymat.load_ipython_extension(ip)


# ### Try a plot in matlab

# In[6]:

get_ipython().run_cell_magic('matlab', '--size 400,350', "t = linspace(0,6*pi,100);\nplot(sin(t))\ngrid on\nhold on\nplot(cos(t), 'r')")


# In[7]:

get_ipython().magic('matplotlib inline')
from matplotlib import pyplot as plt
plt.style.use('ggplot')
#plt.style.use('dark_background')
import numpy as np
fig,ax = plt.subplots(1,1,figsize=(7,6))
t=np.linspace(0,6*np.pi,100)
ax.plot(np.sin(t))
out=ax.plot(np.cos(t))


# Why is python more verbose?  It uses *namespaces* to separate different code modules, so that you can make your own version of plot, or sin, etc.  These modules need to be imported.  See
# [Pine Chapter 2: modules](http://clouds.eos.ubc.ca/~phil/djpine_python/chap2/chap2_basics.html#python-modules)

# ### Python vs. Matlab functions

# ### A matlab function stored in the file average.m:
# 
# ```matlab
#     function y = average(x)
#     if ~isvector(x)
#         error('Input must be a vector')
#     end
#     y = sum(x)/length(x); 
#     end
# ```
# 
# 

# ### The python equivalent

# In[2]:

def average(x):
    x = np.atleast_1d(x)
    y = np.sum(x)/len(x)
    return y


# In[3]:

x=np.linspace(0,100,20)
print(average(x))


# ## Next steps

# 1. Read [Eric Firing's Python Intro Section](http://currents.soest.hawaii.edu/ocn_data_analysis/getting_started.html#how-do-i-learn-a-new-computer-language)
# 2. Start working on the notebooks in his [Python tutorial](http://currents.soest.hawaii.edu/ocn_data_analysis/python_tutorial.html).  Give it 60 minutes and let me know how far you get
# 
# Useful links:
# 
# [Numpy for Matlab users](https://docs.scipy.org/doc/numpy-dev/user/numpy-for-matlab-users.html)
# 
# [David Pine's Python Intro](http://clouds.eos.ubc.ca/~phil/djpine_python/index.html)
# 
# [Stack Overflow](http://stackoverflow.com/questions/tagged/python)
# 

# In[1]:

import test


# In[ ]:



