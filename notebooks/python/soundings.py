
# coding: utf-8

# In[3]:

from importlib import reload
import soundings.readsoundings
reload(soundings.readsoundings)
from soundings.readsoundings import readsound
import pandas as pd
from matplotlib import pyplot as plt


# In[4]:

import glob
out=glob.glob('../*soundings/*july*txt')
sounding_dict=readsound(out[0])
dict_keys=list(sounding_dict.keys())
print(dict_keys[:3])
print(sounding_dict[dict_keys[0]].head())


# In[5]:

outfile='july.h5'
with pd.HDFStore(outfile,'w') as store:
    for key,value in sounding_dict.items():
        thetime=key.strftime("Y%Y_%b_%d_%HZ")
        store.put(thetime,value,format='table')
        


# In[6]:

with pd.HDFStore(outfile,'r') as store:
    sounding=store['Y2006_Jul_12_00Z']
    print(sounding.loc[:5])


# In[55]:

get_ipython().magic('matplotlib inline')
m2km=1.e-3  #convert meters to km
plt.style.use('ggplot')
fig,ax=plt.subplots(1,1,figsize=(8,10))
ax.plot(sounding['temp'],sounding['hght']*m2km,label='temp')
ax.plot(sounding['dwpt'],sounding['hght']*m2km,label='dewpoint')
ax.legend()
out=ax.set(xlabel="temperature (K)",ylabel="height (km)",
      title ="July 2006 soundings, Port Hardy")


# In[ ]:



