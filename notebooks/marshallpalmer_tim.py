import numpy as np
from matplotlib import pyplot as plt

plt.style.use('ggplot')
plt.close("all")

# changed D steps to be much smaller
D=np.arange(0,5,0.00001)
Dmm=D
Dcm=D*0.1
N0=0.08*1.e6*1.e-1 #m**{-3} mm^{-1}
rhol = 1000.
rhoa = 1.2
grav= 9.81
n0 = 0.08 # cm^-4

def calc_fall(curve):
    """Calculates precipitation flux using Villermaux and Bossa (2009)
        Inputs:
            curve (m^-3, mm^-1), Marshall-Palmer distribution
        Outputs:
            rfall (mm/h), fall speed
    """
    wfall = np.sqrt(rhol/rhoa*grav*Dcm)
    
    rfall = np.sum(curve[:-1]/1000.*(np.pi*Dcm[:-1]**3/6)*wfall[:-1]*np.diff(Dcm))
    
    # return in mm/h
    return (rfall*10000/2)
      

            
if __name__ == "__main__":
    Rvec = [1., 5., 15., 25.] # mm/h
    fig, ax = plt.subplots(1,1,figsize=(10,8))
    
    for R in Rvec:
        theLambda=41*R**(-0.21)
        curve=N0*np.exp(-theLambda*Dcm)
        # numerical mean of distribution (weighted integral of the PSD)
        num_mean = np.sum(Dcm[:-1]*curve[:-1]*np.diff(Dcm))/np.sum(curve[:-1]*np.diff(Dcm))
        # analytical mean of distribution
        ana_mean = 1/theLambda
        fallspeed = calc_fall(curve)
        print("R = {} mm/h".format(R))
        print("Numerical mean diameter= {} cm, analytical mean diameter = {} cm".format(num_mean, ana_mean))
        print("Fallspeed = {} mm/hr".format(fallspeed))
        print("========================")
        
        label_str = "{} mm/hr".format(R)

        ax.semilogy(D,curve,label=label_str)
        ax.set_xlabel('Drop diameter (mm)')
        ax.set_ylabel('$n(D) m^{-3} mm^{-1}$')
        ax.set_title('Marshall Palmer distribution for three rain rates')
        ax.legend(loc='best')
