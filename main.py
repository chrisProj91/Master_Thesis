# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 22:55:19 2020

@author: Χρήστος
"""

import numpy as np
from matplotlib import pyplot as plt
from math import pi, exp, cos, sin, sqrt, atan2
Nx = 500
nsteps = 1200
ex = np.zeros(Nx)
Jxm1 = np.zeros(Nx)
jx = np.zeros(Nx)
Jxm2 = np.zeros(Nx)
hy = np.zeros(Nx)
eym1 = np.zeros(Nx)        
eym2 = np.zeros(Nx) 

boundary_low = [0, 0]
boundary_high = [0, 0]
# FDTD parameter and incident wave parameter

eps0 = 8.854 * pow(10,-12)
m0 = 4 * pi * pow(10, -7)
c = 1/sqrt(eps0*m0)
ddx = pow(10, -4) #10^-4
dt = ddx/(2*c)
tc = 300
t_0 = 90

#ddx = 0.01 # Cell size
#dt = ddx / 6e8 # Time step size

# Parameter medium Lorentz 
eps_inf = 1.0  
Deps_p  = 0.5  #  eps-eps_inf
sigma = np.zeros(Nx)
omega_0 = np.zeros(Nx)
a1 = np.zeros(Nx)
a2 = np.zeros(Nx)
a3 = np.zeros(Nx)
c1 = np.zeros(Nx)
c2 = np.ones(Nx)

omega_p = np.zeros(Nx)
delta_p = np.zeros(Nx)
f = np.zeros(Nx)
p = np.zeros(Nx)

for k in range(1, Nx):
    c3 = dt / (eps0)
    
    
for k in range(100, 150):
    sigma[k] = 1.5
    omega_p[k] = 40*pi*pow(10,9)  ## 10^9
    delta_p[k] = (0.01*omega_p[k])
    
    f[k]=1
    a1[k]    = (2- pow((omega_p[k] * dt), 2)/(1+dt*delta_p[k])
    a2[k]    = (dt*delta_p[k]-1)/(1+dt*delta_p[k])
    a3[k]    = Deps_p*eps0*f[k]* pow((omega_p[k]*dt),2) /(delta_p[k]*dt+1)
    c1[k]    = (a3[k]/2)/(2*eps_inf*eps0+sigma[k]*dt+(a3[k]/2))         
    c2[k]    = (2*eps_inf*eps0-sigma[k]*dt)/(2*eps_inf*eps0+sigma[k]*dt+(a3[k]/2))      
    c3[k]    = 2*dt/(2*eps_inf*eps0+sigma[k]*dt+(a3[k]/2))
    
    
#Main FDTD loop
for time in range(1, nsteps+1):

 
    source = exp(-0.5*( pow((tc-time)/t_0), 2) )   #time*dt*3*10^12
    ex[2] = ex[2]+source

    for k in range(1, Nx):
        ex[k]=c1[k]*eym2[k]+c2[k]*ex[k]+(c3[k])*(hy[k-1] - hy[k])/(ddx)-0.5*c3[k]*((1+a1[k])*Jxm1[k]+a2[k]*Jxm2[k])

        eym2[k]=eym1[k]  
        eym1[k]=ex[k]
      