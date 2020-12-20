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
    Nx = 5 #500
    nsteps = 1200
    pi = 3.141593
    eps0 = float( 8.854 * pow(10,-12) )
    m0 = float( 4 * pi * pow(10, -7) )
    c = float( 1/math.sqrt(eps0*m0) )
    ddx = float( pow(10, -4) )
    dt = float( ddx/(2*c) )
    tc = 300
    t_0 = 90
    ex = a1 = jx = Jxm1 = hy = eym1 = eym2 = omega_p = [[0] * Nx]
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
        all_matrices = [self.ex, self.a1, self.jx, self.Jxm1, self.hy, self.eym1, self.eym2, self.omega_p]
        all_matrices_names = ["ex", "a1", "jx", "Jxm1", "hy", "eym1", "eym2", "omega_p"]
        index = 0
        for item in all_matrices:
            print("Now printing matrix : %s " %all_matrices_names[index])
            index = index + 1
            for sub_item in item:
                print("%s\n" %sub_item)

    def step_une(self):
        """
            The 1st step of the procedure to generate the matrices, calculate equations, etc...
        """
        for i in range(0,1): #(100, 151):
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
