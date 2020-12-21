#!/usr/bin/python
# -*- coding: utf-8 -*-

#######################################################
# Created on : December 20, 2020
#
# @author: Christo Vagena
#
#######################################################

import math
# import numpy as np
# from matplotlib import pyplot as plt
# from math import pi, exp, cos, sin, sqrt, atan2

# INFO: Create the OOP style program template and vag will 
#       expand it later on ...


class Diplo:
    """ 
        Global Variables and lists to be used in this program
        Use these variables only inside the class's scope/functions !
    """
    Nx = 5  # 500
    nsteps = 1200
    eps0 = float( 8.854 * pow(10,-12) )
    m0 = float( 4 * math.pi * pow(10, -7) )
    c = float( 1/math.sqrt(eps0*m0) )
    ddx = float( pow(10, -4) )
    dt = float( ddx/(2*c) )
    tc = 300
    t_0 = 90
    # Initialize all lists/matrices with Nx*zeros or Nx*ones respectively
    #
    ex = jx = Jxm1 = hy = eym1 = eym2 = omega_p = sigma = omega_zero = delta_p = f = p = [[0] * Nx]
    a1 = a2 = a3 = [[0] * Nx]
    c1 = [[0] * Nx]
    c2 = [[1] * Nx]
    # Cell size
    ddx = 0.01
    # Time step size | dt = dx/(6*10^8)
    dt = float(ddx/(6*pow(10, 8)))
    # Parameter medium Lorentz
    eps_inf = 1.0
    # eps-eps_inf | Deps_p = 0.5
    deps_p = 0.5

    def __init__(self):
        """
            Contructor for the Diplo class
        """
        self.init_all_lists()

    def init_all_lists(self):
        """
            A simple (void) bundle function that demonstrates how matrices work in Python
            so you can see and use later in your program ...
        """
        all_matrices = [self.ex, self.a1, self.jx, self.Jxm1, self.hy, self.eym1, self.eym2, self.omega_p, \
                        self.sigma, self.omega_zero, self.delta_p, self.f, self.p, self.a2, self.a3]
        all_matrices_names = ["ex", "a1", "jx", "Jxm1", "hy", "eym1", "eym2", "omega_p", \
                               "sigma", "omega_zero", "delta_p", "f", "p", "a2", "a3"]
        index = 0
        for item in all_matrices:
            print("Now printing matrix : %s " %all_matrices_names[index])
            index = index + 1
            for sub_item in item:
                print("%s\n" %sub_item)

    def fill_c3(self):
        for k in range(1, Nx):
            c3 = float(dt/eps0)

    def step_une(self):
        """
            The 1st step of the procedure to generate the matrices, calculate equations, etc...
        """
        for i in range(100, 151):
            self.a1[i]= (2- pow((omega_p[i] * dt), 2)/(1+dt*delta_p[i]) )

    def ret_ex_size(self):
        """ 
            Test function for vag to check how everything works out in Python OOP
        """
        return len(self.ex)


def main():
    """
        main function
    """
    print("main function called - starting main function ... \n")
    DiploObject = Diplo()
    # print("Length of 'ex' list is: " + str(DiploObject.ret_ex_size()))

if __name__ == "__main__":
    main()
