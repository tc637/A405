import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from constants import constants as c
from new_thermo import convertSkewToTemp,convertTempToSkew,theta,wsat,thetaes

def convecSkew(ax):
      """       
      Usage:  convecSkew(ax)
      Input:  figNum = integer
       Takes any integer, creates figure(figNum), and plots a
       skewT logp thermodiagram.
      Output: skew=30 and the handle for the plot
      """
      yplot = range(1000,190,-10)
      xplot = range(-300,-139)
      pvals = np.size(yplot)
      tvals = np.size(xplot)
      temp = np.zeros([pvals, tvals])
      theTheta = np.zeros([pvals, tvals])
      ws = np.zeros([pvals, tvals])
      theThetae = np.zeros([pvals, tvals])      
      skew = 30 #skewness factor (deg C)

      # lay down a reference grid that labels xplot,yplot points 
      # in the new (skewT-lnP) coordinate system .
      # Each value of the temp matrix holds the actual (data)
      # temperature label (in deg C)  of the xplot, yplot coordinate.
      # pairs. The transformation is given by W&H 3.56, p. 78.  Note
      # that there is a sign difference, because rather than
      # taking y= -log(P) like W&H, I take y= +log(P) and
      # then reverse the y axis         
      
      for i in yplot:
            for j in xplot:
                  # Note that we don't have to transform the y
                  # coordinate, as it is still pressure.
                  iInd = yplot.index(i)
                  jInd = xplot.index(j)
                  temp[iInd, jInd] = convertSkewToTemp(j, i, skew)
                  Tk = c.Tc + temp[iInd, jInd]
                  pressPa = i * 100.
                  theTheta[iInd, jInd] = theta(Tk, pressPa)
                  ws[iInd, jInd] = wsat(Tk, pressPa)
                  theThetae[iInd, jInd] = thetaes(Tk, pressPa)
                  
      #
      # Contour the temperature matrix.
      #

      # First, make sure that all plotted lines are solid.
      mpl.rcParams['contour.negative_linestyle'] = 'solid'
      tempLabels = range(-40, 50, 10)
      tempLevs = ax.contour(xplot, yplot, temp, tempLabels, \
                            colors='k')
      
      #
      # Customize the plot
      #
      ax.set_yscale('log')
      locs = np.array(range(100, 1100, 100))
      labels = locs
      ax.set_yticks(locs)
      ax.set_yticklabels(labels) # Conventionally labels semilog graph.
      ax.set_ybound((200, 1000))
      plt.setp(ax.get_xticklabels(), weight='bold')
      plt.setp(ax.get_yticklabels(), weight='bold')
      ax.yaxis.grid(True)

      
      thetaLabels = list(range(200, 390, 10))
      thetaLevs = ax.contour(xplot, yplot, theTheta, thetaLabels, \
                        colors='b')


      wsLabels =[0.1,0.25,0.5,1,2,3] + list(range(4, 20, 2)) + [20,24,28]

      wsLevs = ax.contour(xplot, yplot, (ws * 1.e3), wsLabels, \
                        colors='g')

      thetaeLabels = np.arange(250, 410, 10)
      thetaeLevs = ax.contour(xplot, yplot, theThetae, thetaeLabels, \
                        colors='r') 
      
      # Transform the temperature,dewpoint from data coords to
      # plotting coords.
      ax.set_title('skew T - lnp chart')
      ax.set_ylabel('pressure (hPa)')
      ax.set_xlabel('temperature (deg C)')

      #
      # Crop image to a more usable size
      #    
      

      TempTickLabels = range(-15, 40, 5)

      TempTickCoords = TempTickLabels
      skewTickCoords = convertTempToSkew(TempTickCoords, 1.e3, skew)
      ax.set_xticks(skewTickCoords)
      ax.set_xticklabels(TempTickLabels)

      skewLimits = convertTempToSkew([-15, 35], 1.e3, skew)

      ax.axis([skewLimits[0], skewLimits[1], 300, 1.e3])
      
      #
      # Create line labels
      #
      fntsz = 9 # Handle for 'fontsize' of the line label.
      ovrlp = True # Handle for 'inline'. Any integer other than 0
                # creates a white space around the label.
                
      thetaeLevs.clabel(thetaeLabels, inline_spacing=0, inline=ovrlp, fmt='%5d', fontsize=fntsz,use_clabeltext=True)
      tempLevs.clabel(inline=ovrlp, inline_spacing=0,fmt='%2d', fontsize=fntsz,use_clabeltext=True)
      thetaLevs.clabel(inline=ovrlp, inline_spacing=0,fmt='%5d', fontsize=fntsz,use_clabeltext=True)
      wsLevs.clabel(inline=ovrlp, inline_spacing=0, fmt='%2d', fontsize=fntsz,use_clabeltext=True)
      #print thetaeLabels
      #
      # Flip the y axis
      #
      
      ax.invert_yaxis()
      ax.figure.canvas.draw()
      
      return skew, ax

if __name__== "__main__":
      skew, ax =convecSkew(1)
      plt.show()
      
