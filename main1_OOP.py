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
    # make envelope also global to call it from plot() def
    envelope = 0
    
    def __init__(self, nsteps): 
        """
            Constructor
        """
        self.nsteps = nsteps
        
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
        self.envelope=np.abs(hilbert(self.ex))  
        maximum = max(self.envelope)
        print("Maximum value is : %.5f" % maximum)

    def getMaximum(self):
        """
            return max(envelope), 
            needed to plot later from main def()
        """
        return max(self.envelope)
        
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

    def run_procedure(self):
        """
            A bundle to run all needed functions at once,
            in order not to call them one by one from main 
            function
        """
        self.clc_c3()
        self.media()
        self.Main_loop()
        self.envel()
        self.plot()

class Ploter:
    """
        A second class to plot based on steps and the value returned by
        the first Gain_media class
    """
    def __init__(self, values, steps):
        """
            Plotter class constructor
            @steps  : the used steps from main function
            @values : the values (maximum) returned when Gain_media class is called
        """
        self.x = values
        self.y = steps

    def plot(self):
        """
            plot
        """
        plt.plot(self.y, self.x)
        plt.show()

        
def main():
    """
        main def()
    """
    time_steps = [10, 20, 30, 40, 50]
    maximum_values = []
    i = 1
    for step in time_steps:
        print("##################################")
        print("Creating object %d : " % i)
        print("Calculating value for Steps: %d" % step)
        visual_obj = Gain_media(step)
        visual_obj.run_procedure()
        # save max values in the list 
        maximum_values.append(visual_obj.getMaximum())
        i += 1
        print("##################################\n")

    # print(maximum_values)

    print("Ploting ...")
    plt.plot(maximum_values, time_steps)
    plt.show()

    # visual = Gain_media(time_steps)  
    # visual.run_procedure()

if __name__ == "__main__":
    main()    
        
        

        
        
        
        
            
        