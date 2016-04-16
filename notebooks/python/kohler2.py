
# coding: utf-8

#  Show using two Taylor series expansions how to get  Lohmann 6.24:
# 
# $S= a_w \exp \left [ \frac{2\sigma}{\rho_l R_v T r} \right ] \approx \left ( 1 + \frac{a}{r} - \frac{b}{r^3} \right )$
# 
# 
# **First taylor series:**
# 
# $a_w = \left ( \frac{n_w}{n_w + n_s} \right ) = \left ( \frac{1}{1 + n_s/n_w} \right ) = (1 + n_s/n_w)^{-1}$
# 
# If $n_s/n_w  = x \ll 1$  then expand $(1 + x)^{-1}$ in a Taylor series about x=0:
# 
# $f^\prime (0) = -(1 + 0 )^{-2} = -1$, $f^{\prime \prime} (0) = 2 ( 1 + 0 )^{-3} = 2$
# 
# so $(1 + x)^{-1} \approx 1 - x + 2x^2/2 + \ldots$
# 
# **Second taylor series:**
# 
# $\exp \left [ \frac{2\sigma}{\rho_l R_v T r} \right ] = \exp \left ( \frac{a}{r} \right ) = \exp(y)$
# 
# if $y \ll 1$ then expand exp in a taylor series about y=0:
# 
# $f^\prime (0) =\exp(0) = 1$
# 
# $f^{\prime \prime}(0) = 1$
# 
# so $\exp(y) \approx 1 + y  + y^2/2 + \ldots$
# 
# **Combining: **
# 
# $a_w \exp \left [ \frac{2\sigma}{\rho_l R_v T r} \right ] \approx (1 - x + x^2 )(1 + y + y^2/2) = 1 + y - x + y^2/2 + x^2 + \ldots$
# 
# and keeping only first order terms:
# 
# $S = 1 + y - x = 1 + \frac{a}{r} - n_s/n_w$
# 
# and since $n_w \propto r^3$: $S = 1 + y - x = 1 + \frac{a}{r} - \frac{b}/{r^3}$
# 
# 
# 

# 2\.  For the aerosol defined in the kohler.ipynb notebook ($10^{-19}$ kg, ammonium sulphate), inside a droplet of radius $r=1\ \mu m$
#       find the size of the smallest term you've kept (either $\frac{a}{r}$ or $\frac{b}{r^3}$ and compare it to
#       the size of the largest term you've dropped.   Repeat this for a droplet of radius   $r=0.1\ \mu m$
# 
# 

# In[1]:

from a405thermo.constants import constants as c
def find_terms(r):
    Tinit=c.Tc + 15 #Temperature K
    Ms=132 #ammonium sulphae kg/Kmole
    Mw=18. #water kg/Kmole
    Sigma=0.075  #N/m
    vanHoff=3. #van hoff for ammonium bisulphate
    #calculate kohler coefficients:
    a=(2.*Sigma)/(c.Rv*Tinit*c.rhol)  #curvature term
    mass_aero = 1.e-19  #kg
    b=(vanHoff*Mw)/((4./3.)*np.pi*c.rhol*Ms)*mass_aero   #solution term
    return (a/r,b/r**3.)

kelvin_term,raoult_term=find_terms(1.e-6)
print('kelvin: {:8.3e} vs raoult {:8.3e}'.format(kelvin_term,raoult_term))

kelvin_term,raoult_term=find_terms(1.e-7)
print('kelvin: {:8.3e} vs raoult {:8.3e}'.format(kelvin_term,raoult_term))


#  Suppose you have $r_l$ =1 g/kg of liquid water spread among N spherical drops of radius 10 $\mu m$.  Compare the surface energy of this
#       population (which we've been neglecting) with the enthalpy required to vaporize it:  $l_v r_l$.  Is it negligible in comparison?
# 
# 

# In[11]:

sigma = 0.075
rhol=1000
rl=1.e-3
r=1.e-5
N = rl/(4./3.*np.pi*r**3.*rhol)
surface_energy = N*sigma*4*np.pi*r**2.
vapor_energy = c.lv0*rl
print('there are N={:6.3g} drops/kg'.format(N))
print('surface_energy {:6.3f} J/kg, vapor_energy {:6.3f} J/kg'.format(surface_energy,vapor_energy))


# In[ ]:



