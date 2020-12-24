# -*- coding: utf-8 -*-
"""
Created on Thu Dec 24 19:19:57 2020

@author: Χρήστος
"""

import numpy as np
from matplotlib import pyplot as plt
from math import pi, exp, cos, sin, sqrt, atan2
from scipy.signal import hilbert, chirp

class Gain_media:
    """
        Class for computing electric field
    """
    Nx=5000
    # FDTD parameter and incident wave parameter
    eps0 = 1# float(8.854 * pow(10, -12))
    m0 = 1#float(4 * pi * pow(10,-7))
    c = float(1/sqrt(eps0*m0))
    ddx = pow(10, -3)#pow(10, -4)
    dt = float(ddx/(2*c))
    ex = np.zeros(Nx)
    Jxm1 = np.zeros(Nx)
    jx = np.zeros(Nx)
    Jxm2 = np.zeros(Nx)
    hy = np.zeros(Nx)
    eym1 = np.zeros(Nx)        
    eym2 = np.zeros(Nx)
    # Initial boundary conditions
    ey_low_m1=0
    ey_low_m2=0
    ey_low_m3=0
    ######
    ex_high_m1=0
    ex_high_m2=0
    ex_high_m3=0
    # Parameter medium Lorentz
    eps_inf = 1.0  
    Deps_p  =0.5
    sigma = np.zeros(Nx)
    omega_0 = np.zeros(Nx)
    omega_p = np.zeros(Nx)
    delta_p = np.zeros(Nx)
    f = np.zeros(Nx)
    p = np.zeros(Nx)
    a1 = np.zeros(Nx)
    a2 = np.zeros(Nx)
    a3 = np.zeros(Nx)
    c1 = np.zeros(Nx)
    c2 = np.ones(Nx)
    c3 = np.zeros(Nx)
    
    def __init__(self, nsteps): 
        """
            Constructor
        """
        self.nsteps = nsteps
        # self.clc_c3()
        # self.media()
        # self.Main_loop()
        # self.envel()
        # self.plot()
        
    def clc_c3(self):
        """
            Calculate parameter c3 in full space
        """    
        for k in range(1, self.Nx):
            self.c3[k] = self.dt / (self.eps0)
            
    def media(self):
        """
        Create the media
        """
        for k in range(3200, 3500):
            self.sigma[k] = -5      #-16.2#-6.2#-0.62
            self.omega_p[k] = 5     #40*pi*pow(10, 9)
            self.delta_p[k] = (0.01 / 2)  #(0.01*omega_p[k])
    
            self.f[k] = -1
            self.a1[k] = (2-pow((self.omega_p[k] * self.dt), 2))/(1+self.dt * self.delta_p[k])
            self.a2[k] = (self.dt*self.delta_p[k]-1)/(1+self.dt*self.delta_p[k])
            self.a3[k] = self.Deps_p * self.eps0 * self.f[k] * pow((self.omega_p[k] * self.dt), 2)/(self.delta_p[k] * self.dt+1)
            self.c1[k] = (self.a3[k]/2)/(2 * self.eps_inf * self.eps0 + self.sigma[k] * self.dt + (self.a3[k]/2))         
            self.c2[k] = (2*self.eps_inf * self.eps0 - self.sigma[k] * self.dt)/(2 * self.eps_inf * self.eps0 + self.sigma[k] * self.dt + (self.a3[k]/2))      
            self.c3[k] = (2*self.dt)/(2 * self.eps_inf * self.eps0 + self.sigma[k] * self.dt + (self.a3[k]/2))
    
    def Main_loop(self):
        """
         Main loop calculating electric field
        """
        for time in range(1, self.nsteps+1):
            time=time+1
    #source = sin(50*time*dt)*exp(-0.5 * ((tc-time) / t_0) **2)  #time*dt*3*10^12
            source = sin(50 * time * self.dt) * exp (-0.5 * ((1000-time) / 232) **2)
                 
            self.ex[2000] = self.ex[2000] + source
            
            #Mur absorbing condition

            self.hy[0] = self.ey_low_m2
            self.ey_low_m2 = self.ey_low_m1
            self.ey_low_m1 = self.hy[1]
 
 
 
            self.ex_high_m3 = self.ex_high_m1
            self.ex_high_m2 = self.ex_high_m1
            self.ex_high_m1 = self.ex[-1]
     
            self.ex[0] = self.ex[0] + (self.dt/(self.ddx * self.eps0)) * (self.ey_low_m3 - self.hy[0])
            
            #Electric field update
            for k in range(0, self.Nx):
                self.ex[k]=self.c1[k] * self.eym2[k] + self.c2[k] * self.ex[k] + (self.c3[k]) * (self.hy[k-1] - self.hy[k]) / (self.ddx)-0.5*self.c3[k] * ((1+self.a1[k]) * self.Jxm1[k] + self.a2[k] * self.Jxm2[k])

                self.eym2[k]=self.eym1[k]  
                self.eym1[k]=self.ex[k]
                
            for k in range(0, self.Nx):
                self.jx[k]=self.a1[k] * self.Jxm1[k] + self.a2[k] * self.Jxm2[k] + (self.a3[k]/(2 * self.dt))*(self.ex[k]-self.eym2[k])
        
        
        
            for k in range(1, self.Nx-1):
                self.hy[k] = self.hy[k] + (self.dt/(self.ddx * self.m0)) * (self.ex[k] - self.ex[k+1])
        
            self.hy[self.Nx-1] = self.hy[self.Nx-1] + (self.dt/(self.ddx * self.m0))*( self.ex[self.Nx-1] - self.ex_high_m3)     
            
            
    def envel(self):
        """
        Calculate the envelope of ex and maximum value 
        """
        envelope=np.abs(hilbert(self.ex))  
        maximum = max(envelope)
        print(maximum) 
        
    def plot(self):
        """
        Plot ex with Nx for different nsteps
        """
        plt.rcParams['font.size'] = 12
        fig = plt.figure(figsize=(8, 7))   
        plt.plot(self.envelope, color='k', linewidth=1)
        #plt.xticks(np.arange(0, 1000, step=200))
        plt.xlim(0, self.Nx)
        #plt.yticks(np.arange(-1, 5, step=1))
        plt.ylim(0, 3)
        plt.axvline(3200, color='b')
        plt.axvline(3500, color='b')
        
def main():
    """
        main def()
    """
    time_steps = 1000
    visual = Gain_media(time_steps)  
    visual.clc_c3()
    visual.media()
    visual.Main_loop()
    visual.envel()
    visual.plot()

if __name__ == "__main__":
    main()    
        
        

        
        
        
        
            
        