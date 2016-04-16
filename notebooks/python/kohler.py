
# coding: utf-8

# In[ ]:

### Kohler problem

1) Plot the Kohler curve for the aerosol on page 2 of the Kohler notes


# In[19]:

from a405thermo.constants import constants as c
Tinit=c.Tc + 15 #Temperature K
Sinit=0.9 #e/es
Ms=132 #ammonium sulphae kg/Kmole
Mw=18. #water kg/Kmole
Sigma=0.075  #N/m
vanHoff=3. #van hoff for ammonium bisulphate
#calculate kohler coefficients:
a=(2.*Sigma)/(c.Rv*Tinit*c.rhol)  #curvature term
mass_aero = 1.e-19  #kg
b=(vanHoff*Mw)/((4./3.)*np.pi*c.rhol*Ms)*mass_aero   #solution term
rcrit = (3*b/a)**0.5  #critical radius
print('critical radius = {:5.3f} microns'.format(rcrit*1.e6))


# In[20]:

get_ipython().magic('matplotlib inline')
def find_S(r):
    S= 1 + a/r - b/r**3
    return S

    
plt.close('all')
import matplotlib
matplotlib.style.use('ggplot')
fig,ax = plt.subplots(1,1)
rvals=np.linspace(0.05,0.3,50)*1.e-6
Svals = find_S(rvals)
ax.plot(rvals*1.e6,(Svals -1)*1.e2)
out=ax.set(ylim=[-1,1],xlim = [0.05,0.3],
       ylabel='Supersaturation (e/es - 1) %',xlabel='radius (microns)')


# ### 2.  Find haze particle equilibrium radius at S=0.9

# In[21]:

def find_diff(r,S_target):
    return S_target - find_S(r)

S_target = 0.90
r_start = 0.1e-6
from a405thermo.rootfinder import find_interval, fzero
brackets = np.array(find_interval(find_diff,r_start,S_target))
print('left bracket = {:8.3e} microns, right bracket={:8.3e} microns'.format(*(brackets*1.e6)))

equil_rad = fzero(find_diff,brackets,S_target)
print('equlibrium radius at S={} is {:5.3f} microns'.format(S_target,equil_rad*1.e6))


# In[ ]:



