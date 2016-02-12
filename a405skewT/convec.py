#calculate thermodynamic variables
#for a convectively unstable layer
from a405thermo.constants import constants as c
from a405skewT.makeSkewII import makeSkewWet
from a405thermo.thermlib import find_thetaep, find_rsat, tinvert_thetae, find_Td, find_lcl
from a405thermo.thermlib import convertTempToSkew
from a405thermo.thermlib import find_Tmoist
import numpy as np
import matplotlib.pyplot as plt

pa2hPa = 1.e-2
hPa2pa = 100.

def makePlot(ax,Temp=None,Tdew=None,Tpseudo=None,press=None,Tlcl=None,pLCL=None,botLabel='LCL bot (835 hPa)',topLabel='LCL top (768 hPa)'):
  """
    Draw a tephigram of a convectively unstable layer with Temp, Tdew, Tpseudo soundings, as well as the LCL of the top and the bottom of
    the layer
  """
  xplot = convertTempToSkew(Tlcl[0] - c.Tc, pLCL[0] * pa2hPa, skew)
  bot, = ax.plot(xplot, pLCL[0] * pa2hPa, 'ro', markersize=12, markerfacecolor='r',label=botLabel)
  xplot = convertTempToSkew(Tlcl[-1] - c.Tc, pLCL[-1] * pa2hPa, skew)
  top, = ax.plot(xplot, pLCL[-1] * pa2hPa,'bd', markersize=12,
                 markerfacecolor='b',label = topLabel)
  xplot = convertTempToSkew(Tpseudo - c.Tc, press, skew)
  thetaEhandle, = ax.plot(xplot, press, 'c-', linewidth=2.5,label=r'$\theta_e$')
  xplot1 = convertTempToSkew(Temp - c.Tc, press, skew)
  Thandle, = ax.plot(xplot1, press, 'k-', linewidth=2.5,label='Temp (deg C)')
  xplot2 = convertTempToSkew(Tdew - c.Tc, press, skew)
  TdHandle, = ax.plot(xplot2, press, 'b--', linewidth=2.5,label='Dewpoint (deg C)')
  ax.legend(loc='best',numpoints=1)
  base=press[0]
  title_string = 'convectively unstable sounding: base at {} hPa'.format(base)
  ax.set_title(title_string)
  return ax

def lift_sounding(rTotal,theThetae,press,):
  #
  # return temp,dewpoint and Tspeudo soundings for an rT,thetae sounding at pressure levels press
  #
  numPoints, = press.shape
  for i in range(0, numPoints):
    #find the actual temperature and dewpoint
    Temp[i], rv, rl = tinvert_thetae(theThetae[i], rTotal[i], press[i] * hPa2pa)
    Tdew[i] = find_Td(rv, press[i] * hPa2pa)
    Tpseudo[i] = find_Tmoist(theThetae[i], press[i] * hPa2pa)
  return Tdew,Temp,Tpseudo


if __name__ == "__main__":
  
  plt.close('all')

  #
  # construct a 100 hPa thick layer with a convectively unstable sounding
  #
  Tbot = 20.
  Ttop = 17.5
  Tdewbot = 15
  Tdewtop = 5.5

  Pbot = 900.
  Ptop = 800.

  #calculate the temperature and dewpoint sounding
  #assuming linear profiles
  slope = (Ttop - Tbot) / (Ptop - Pbot)
  press = np.arange(Pbot, Ptop - 10, -10)
  Temp = (Tbot + slope * (press - Pbot))
  slope = (Tdewtop - Tdewbot) / (Ptop - Pbot)
  Tdew = (Tdewbot + slope * (press - Pbot))

  #
  #figure 1: plot the T,Tdew profile only
  #
  
  plt.close('all')

  fig, ax = plt.subplots(1, 1, figsize=(10, 10))
  ax, skew = makeSkewWet(ax,corners=[5,25])

  #zoom the axis to focus on layer
  skewLimits = convertTempToSkew([5, 25], 1.e3, skew)
  ax.axis([skewLimits[0], skewLimits[1], 1000, 600])
  xplot1 = convertTempToSkew(Temp, press, skew)
  #plot() returns a list of handles for each line plotted
  Thandle, = ax.plot(xplot1, press, 'k-', linewidth=2.5,label='Temp (deg C)')
  xplot2 = convertTempToSkew(Tdew, press, skew)
  TdHandle, = ax.plot(xplot2, press, 'b--', linewidth=2.5,label='Dewpoint (deg C)')
  ax.set_title('convectively unstable sounding: base at 900 hPa')
  ax.legend(numpoints=1,loc='best')
  fig.savefig('initial_sound.png')
  fig.savefig('initial_sound.pdf')

  #
  # figure 2, add thetae sounding and LCL for layer top and bottom
  #

  #put on the top and bottom LCLs and the thetae sounding
  Tlcl = np.zeros_like(press)
  pLCL = np.zeros_like(press)
  theTheta = np.zeros_like(press)
  theThetae = np.zeros_like(press)
  Tpseudo = np.zeros_like(press)
  rTotal = np.zeros_like(press)

  #
  # calculate the rTotal,thetae sounding from the original dewpoint and temperture
  #

  numPoints, = press.shape
  for i in range(0, numPoints):
      rTotal[i] = find_rsat(Tdew[i] + c.Tc, press[i] * hPa2pa)
      Tlcl[i], pLCL[i] = find_lcl(Tdew[i] + c.Tc, Temp[i] + c.Tc,
                                  press[i] * hPa2pa)
      theThetae[i] = find_thetaep(Tdew[i] + c.Tc, Temp[i] + c.Tc,
                                  press[i] * hPa2pa)
      #find the temperature along the pseudo adiabat at press[i]
      Tpseudo[i] = find_Tmoist(theThetae[i], press[i] * hPa2pa)


  #
  # given the total water and thetae, calcultate temp,dewpoint for pressure vector press
  #
  Tdew, Temp, Tpseudo = lift_sounding(rTotal,theThetae,press)

  fig, ax = plt.subplots(1, 1, figsize=(10, 10))
  ax, skew = makeSkewWet(ax,corners=[5,25])

  fig_dict=dict(press=press,Tpseudo=Tpseudo,Temp=Temp,Tdew=Tdew,Tlcl=Tlcl,pLCL=pLCL,botLabel='LCL bot (835 hPa)',
                topLabel='LCL top (768 hPa)')
  ax = makePlot(ax,**fig_dict)   

  fig.savefig('base900_thetae.png')
  fig.savefig('base900_thetae.pdf')


  # #figure 3: lift cloud base by 50 hPa to 850 hPa

  press = press - 50.
  Tdew, Temp, Tpseudo = lift_sounding(rTotal,theThetae,press)
  fig_dict['press']=press

  fig, ax = plt.subplots(1, 1, figsize=(10, 10))
  ax, skew = makeSkewWet(ax,corners=[5,25])
  ax = makePlot(ax,**fig_dict)   
  fig.savefig('base850_thetae.png')

  #
  # figure 4 -- lift by 14.7 hPa to 835.3
  #

  press = press - 14.7
  Tdew, Temp, Tpseudo = lift_sounding(rTotal,theThetae,press)
  fig_dict['press']=press

  fig, ax = plt.subplots(1, 1, figsize=(10, 10))
  ax, skew = makeSkewWet(ax,corners=[5,25])
  ax = makePlot(ax,**fig_dict)   

  fig.savefig('base835_thetae.png')
  fig.savefig('base835_thetae.pdf')
  
  #
  # figure 5 -- lift by 10.3 hPa to 825
  #

  press = press - 10.3
  Tdew, Temp, Tpseudo = lift_sounding(rTotal,theThetae,press)
  fig_dict['press']=press

  fig, ax = plt.subplots(1, 1, figsize=(10, 10))
  ax, skew = makeSkewWet(ax,corners=[5,25])
  ax = makePlot(ax,**fig_dict)   
  fig.savefig('base825_thetae.png')
  fig.savefig('base825_thetae.pdf')

  #
  # figure 6 -- lift by 25 hPa to 800 hPa
  #

  press = press - 25
  Tdew, Temp, Tpseudo = lift_sounding(rTotal,theThetae,press)
  fig_dict['press']=press

  fig, ax = plt.subplots(1, 1, figsize=(10, 10))
  ax, skew = makeSkewWet(ax,corners=[5,25])
  ax = makePlot(ax,**fig_dict)   
  fig.savefig('base800_thetae.png')
  fig.savefig('base800_thetae.pdf')

  #
  # figure 7 -- lift by 32.25 to 768 hPa
  #

  press = press - 32.25
  Tdew, Temp, Tpseudo = lift_sounding(rTotal,theThetae,press)
  fig_dict['press']=press

  fig, ax = plt.subplots(1, 1, figsize=(10, 10))
  ax, skew = makeSkewWet(ax,corners=[5,25])
  ax = makePlot(ax,**fig_dict)   

  fig.savefig('base768_thetae.png')
  fig.savefig('base768_thetae.pdf')
