
# coding: utf-8

# #### ATSC 405 exercise: eqstate.py
# 
# Hadleigh Thompson
# 
# January 2016
# 

# #### Instructions
# 
# Email a ipython notebook that defines a function called “eqstate” and uses it to calculate the density of dry air at a pressure of 80 kPa and temperatures of temp=[280, 290, 300] K.

# #### Answer

# In[63]:

import numpy as np
import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')



def eqstate(listoftemps):
    """
    list of floats --> float
    
    uses the ideal gas law to calculate density of an air parcel at fixed pressure of 80 kPa
    
    P = rho * Rd * T    # Stull, pg 14, eqn 1.18
    
    in: list of temperatures (K)
    
    out: list of densities of dry air (kg.m^-3)
    
    """
    P = 80e+3           # Pressure, Pa
    Rd = 287.05        # gas constant for dry air, J mole^-1 K^-1
    
    list_of_rho = [(P / (Rd * temp)) for temp in listoftemps]
    
    return list_of_rho


temps = [280, 290, 300]

rho = eqstate(temps)

results = zip(temps,rho)

[print('At temp: %2.1f K, density of dry air = %2.8f kg.m^-3' %(t,r)) for (t,r) in (results)]

plt.plot(temp, rho, label='density at 80kPa')
plt.xlabel('Temperature $K$')
plt.ylabel('Density $(kg.m^{-3})$')
plt.legend()

