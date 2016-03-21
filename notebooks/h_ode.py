# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 09:59:58 2016

@author: Timothy
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

plt.style.use("ggplot")
plt.close("all")

def calc_derivs(y, t, G, SS, w, beta):
    rad, height = y
    dydt = [G*SS/rad, w - beta*rad**2]
    return dydt
    
G = 0.1*1.2/1000
SS = 0.001
w = 20. # m/s
beta = 1.2e8

y0 = [0.1e-6, 0] # m
t = np.linspace(0, 10, 101)
numsol = odeint(calc_derivs, y0, t, args=(G,SS,w,beta))

anlsol = y0[1] + t[-1]*(w - beta*y0[0]**2 - beta*G*SS*t[-1])

print("Numerical solution = {} m , analytic solution = {} m".format(numsol[-1][1], anlsol))
    
anlsol = y0[1] + t*(w - beta*y0[0]**2 - beta*G*SS*t)
fig,ax = plt.subplots(1,1,figsize=(6,6))

ax.plot(t, numsol[:,1], "o", label="Numerical solution")
ax.plot(t, anlsol, label="Analytic solution")
ax.set_xlabel("Time (s)")
ax.set_ylabel("h (m)")
ax.set_title("Comparison of Numerical and Analytic Solutions for Cloud Droplet Height")
ax.legend(loc="best")


#############################################

def calc_derivs1(y, t, G, SS, w, beta):
    rad, height = y
    dydt = [G*SS/rad, w - beta*rad]
    return dydt

G = 0.7*1.2/1000
SS = 0. # RH = 60%
w = 0. # m/s
beta = 6000.

y0 = [1000e-6, 0] # m
t = np.linspace(0, 1000, 1000)
numsol = odeint(calc_derivs1, y0, t, args=(G,SS,w,beta))

fig,ax1 = plt.subplots(1,1,figsize=(6,6))

ax1.plot(t, numsol[:,1], label="Numerical solution")
ax1.set_xlabel("Time (s)")
ax1.set_ylabel("h (m)")
ax1.set_title("Day 31 Problem")
ax1.legend(loc="best")

ind = len(numsol[:,1]) - np.searchsorted(numsol[:,1][::-1], -5000, side="left") 
print("radius at ground = {} m, time to reach ground = {} s".format(numsol[ind,0], t[ind]))



