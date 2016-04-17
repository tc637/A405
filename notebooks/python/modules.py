
# coding: utf-8

# # Using modules in python
# 
# 

# ## Putting our code repository in your import path
# 
# All the functionality in Python is organized in packages and modules, which are described in more detail in [Pine Section 2.4](http://clouds.eos.ubc.ca/~phil/djpine_python/chap2/chap2_basics.html#python-modules).  To see where Python is looking for modules, print the sys.path variable:

# In[8]:

import sys
print(sys.path)


# As in matlab, Python looks through each of these folders in order until if matches the name of the module you want to import.   When you install a module with conda, it goes by default into one of the site-packages folder:

# In[9]:

import site
print(site.getsitepackages())


# In order to import your own modules they need to be either in the folder you started ipython in, or in another folder you've added to sys.path.  We want our A405 git repository to be searched automatically for modules, so we need to add that directory to sys.path every time we start ipython.  To accomplish this, do the following:
# 
# 1.  From a cmd or bash shell, create a default IPython profile:
# 
#     ```
#         > ipython profile create
#     ```
#     
# 2.  Find out where ipython is keeping your profile:
# 
#     ```
#         > ipython locate profile
#         /Users/phil/.ipython/profile_default
#     ```
#     
# 3.  Using an editor, open a file in profile_default/startup  called **00first.py** and add two lines
#     to put your repository in the path (for example) here's a file for Windows machines in the M203 lab:
#     
#     ```
#         import site
#         site.addsitedir('z:\\repos\\A405')
#         print("executed startup file")
#     ```
#     
#     and here's a setting for a Mac:
#     
#     ```
#         import site
#         site.addsitedir('/Users/phil/repos/A405')
#         print("executed startup file")
#     ```
# 
# 4.  Restart ipython and do an  "import sys;print(sys.path)" to confirm that the directory has
#     been added to your path
#     
# 
# 
# 

# In[10]:

import sys
print(sys.path)


# ## Using packages

# In Python, a group of modules in the same directory is called a package.  Python recognizes a directory as a package because the directory contains a file called *__init__.py*.  I have created a package called *a405thermo*, which contains the module *intro.py* that defines the function *hello*.   You should be able to import and run the function like this:

# In[1]:

from a405thermo import intro
intro.hello()


# ## Reloading modules

# One quirk of python is that it ignores any changes in a module unless you explicitly reload it.  That's because checking dozens of modules every time you issue a command would cause python to be frustratingly slow.  So instead you need to import and call the reload function.   To see this, edit the intro.py file (making sure you are off of the master branch) so that hello() prints something different.  Then reload and run the new function like this:

# In[12]:

from importlib import reload
reload(intro)
intro.hello()


# In[ ]:



