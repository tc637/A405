
# coding: utf-8

# ### Retrieve soundings directly using python
# 
# This notebook uses the [wyominglib](https://github.com/phaustin/A405/blob/master/soundings/wyominglib.py) module, which parses  sounding data from [U. Wyoming](http://weather.uwyo.edu/upperair/sounding.html)
# with [BeautifulSoup](http://www.crummy.com/software/BeautifulSoup/)

# In[6]:

from importlib import reload
import soundings.wyominglib
reload(soundings.wyominglib)
from soundings.wyominglib import make_frames
import pandas as pd
from matplotlib import pyplot as plt
import requests
import datetime
from datetime import timezone as tz
import tzlocal   #pip install tzlocal
import h5py
from pprint import pformat
print(help(make_frames))


# ### 1.  Form the url using a template with variable substitution from a dictionary

# In[7]:

#this template holds the text that will always be the same,
# plus variables that will be substituted in based on dictionary values

url_template=("http://weather.uwyo.edu/cgi-bin/sounding?"
              "region={region:s}"
              "&TYPE=TEXT%3ALIST"
              "&YEAR={year:s}"
              "&MONTH={month:s}"
              "&FROM={start:s}"
              "&TO={stop:s}"
              "&STNM={station:s}")


# #### Uncomment one of these dictionaries to get a set of soundings.  Dictionary values are inserted into the url_tempate to form the query url.  See wikipedia on [query strings](https://en.wikipedia.org/wiki/Query_string)

# In[8]:

values=dict(region='samer',year='2013',month='2',start='0100',stop='2800',station='82965')
#values=dict(region='nz',year='2013',month='2',start='0100',stop='2800',station='93417')
#values=dict(region='naconf',year='2013',month='2',start='0100',stop='2800',station='71802')
#values=dict(region='ant',year='2013',month='07',start='0100',stop='2800',station='89009')
url=url_template.format_map(values)
print('here is the url we call to get the soundings: \n',url)


# ### 2. Use the requests module to grab the web page

# In[9]:

#
# Make do_web False to reuse a page for debugging
#
do_web = True
backup_file='backup.txt'
if do_web:
    #
    # grab the web page that is loaded with this url
    #
    html_doc = requests.get(url).text
    print('read {} bytes'.format(len(html_doc)))
    with open(backup_file,'w') as f:
        f.write(html_doc)
    if len(html_doc) < 2000:
        print('debug: short html_doc, something went wrong:',html_doc)
        sys.exit(1)
else:
    with open(backup_file,'r') as f:
        html_doc=f.read()


# ### 3.  Parse the sounding page into two dictionaries
# 
# attr_dict holds attributes ['header', 'site_id','longitude','latitude', 'elevation', 'units']
# 
# sounding_dict holds the soundings indexed by datetime

# In[10]:

attr_dict,sounding_dict = make_frames(html_doc)


# ### Create a timestamp and convert to UTC

# In[11]:

mytz=tzlocal.get_localzone()
now=datetime.datetime.now(tz=mytz)
now=now.astimezone(tz.utc)
timestamp=now.strftime('%Y-%m-%d %H:%M:%S UTC')
print(timestamp)


# ### Add history, query and timestamp attributes to the dictionary, and order the keys

# In[12]:

attr_dict['timestamp']=timestamp
attr_dict['history']="written by test_requests.py"
attr_dict['query']= url
#
# write the keys out in this order
#
key_list=['header', 'site_id','longitude','latitude', 'elevation', 'units','history','query','timestamp']


# ### Use HDFStore to write the soundings out keyed by date

# In[13]:

print(sounding_dict.keys())
name = 'out.h5'    
with pd.HDFStore(name,'w') as store:
    for key,value in sounding_dict.items():
        #need to insert Y in front of the year because
        #h5py groups need to be legal python variables
        #(remove Y to see error message)
        thetime=key.strftime("Y%Y_%b_%d_%HZ")  
        store.put(thetime,value,format='table')


# ### Use h5py to store the attributes as metadata

# In[14]:

with h5py.File(name,'a') as f:
        for key in key_list:
            print('writing key, value: ',key,attr_dict[key])
            f.attrs[key]=attr_dict[key]
        f.close()


# ### Show how to read the data and metadata back in
# 
# Need to trap IOError because there are 4 empty variables in the hdf file
# that can't be read.  See [this discussion](https://github.com/h5py/h5py/issues/279)
# about how to fix this in future h5py releases.

# In[15]:

attributes={}
with h5py.File(name,'r') as f:
    keys=f.attrs.keys()
    for key in keys:
        try:
            attributes[key]=f.attrs[key]
        except IOError:
            print('empty key: ',key)
print('\nread in these attributes: \n\n',pformat(attributes))


# In[16]:

name = 'out.h5'    
separator= '\n' + '+'*30 + '\n'
sounding_dict={}
with pd.HDFStore(name,'r') as store:
    times=store.keys()
    for the_time in times:
        sounding_dict[the_time]=store[the_time]
sounding=sounding_dict[times[0]]
print('{0:}sounding times:{1:}{0:} '.format(separator,times))
print('{0:}columns: {1:}{0:} '.format(separator,pformat(sounding.columns)))
title_string=attributes['header']
index=title_string.find(' Observations at')
location=title_string[:index]
title='{} at {}'.format(location,times[0][2:])
print('title: :',title)
units=attributes['units'].split(';')
units_dict={}
for count,var in enumerate(sounding.columns):
    units_dict[var]=units[count]
print('variables with units: \n',pformat(units_dict))


# In[17]:

get_ipython().magic('matplotlib inline')
key=times[3]  #pick the third sounding
the_sounding=sounding_dict[key]
m2km=1.e-3  #convert meters to km
plt.style.use('ggplot')
fig,ax=plt.subplots(1,1,figsize=(8,10))
ax.plot(the_sounding['temp'],the_sounding['hght']*m2km,label='temp')
ax.plot(the_sounding['dwpt'],the_sounding['hght']*m2km,label='dewpoint')
ax.legend()
out=ax.set(xlabel="temperature (deg C)",ylabel="height (km)",
      title =title)
out=ax.set(ylim=[0,2],xlim=[-30,30])


# In[18]:

fig,ax=plt.subplots(1,1,figsize=(8,10))
for the_time in times:
    the_sound=sounding_dict[the_time]
    ax.plot(the_sound['dwpt'],the_sound['hght']*m2km)
out=ax.set(xlabel="dew point temperature (deg C)",ylabel="height (km)",
      title =location)
out=ax.set(ylim=[0,2],xlim=[-30,30])
print('Dew point temperature between {} to {}'.format(times[0],times[-1]))


# In[ ]:



