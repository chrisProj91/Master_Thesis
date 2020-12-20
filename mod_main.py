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
        BeWare these variables belong only to the class's scope !
    """
    Nx = 500
    nsteps = 1200
    pi = 3.141593
    eps0 = float( 8.854 * pow(10,-12) )
    m0 = float( 4 * pi * pow(10, -7) )
    c = float( 1/math.sqrt(eps0*m0) )
    ddx = float( pow(10, -4) )
    dt = float( ddx/(2*c) )
    tc = 300
    t_0 = 90
    ex = a1 = jx = Jxm1 = hy = eym1 = eym2 = omega_p = list()
    def __init__(self):
        """
            Contructor for the Diplo class
        """
        self.init_all_lists()

    def init_matrix(self, my_matrix):
        """
            A template simple funtion to initialize a matrix with zeros and size=500
        """
        for i in range(1, self.Nx+1):
            my_matrix.append(0)

    def init_all_lists(self):
        """
            A simple (void) bundle function that calls the template init_matrix() function to initialize
            all needed lists that are declared in the beginning of this simple program
            Create a list of all lists
            and use this to initialize every list of this list
        """
        all_matrices = [self.ex, self.a1, self.jx, self.Jxm1, self.hy, self.eym1, self.eym2, self.omega_p]
        index = 0
        for item in all_matrices:
            print("Initializing %d item from all_matrices list of lists ...\n" %(index))
            self.init_matrix(item)
            index = index + 1

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
    print("main function called - starting main function ... ")
    DiploObject = Diplo()
    print("Length of 'ex' list is: " + str(DiploObject.ret_ex_size()))

if __name__ == "__main__":
    main()
