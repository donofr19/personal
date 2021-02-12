""" 
########################################################################
# Vincenzo Donofrio
# AST 304, Fall 2020
# Michigan State University
########################################################################
"""

'''
Routine that calculates specific heating rate from pp chain
hydrogen burning; default pp_factor is 1.0
'''

import numpy as np

def pp_rate(T,rho,XH,pp_factor=1.0):
    """
    Specific heating rate from pp chain hydrogen burning. Approximate rate 
    taken from Hansen, Kawaler, & Trimble.
    
    Arguments
        T, rho
            temperature [K] and density [kg/m**3]
        XH
            mass fraction hydrogen
        pp_factor
            multiplicative factor for rate
    Returns
        heating rate from the pp-reaction chain [W/kg]
    """
    
    Tn = T / 10**9 #Converts T from Kelvin to GigaKelvin
    
    
    rate = pp_factor*(2.4e-3)*rho*XH**2/Tn**(2/3)*np.exp(-3.380/Tn**(1/3))
    return rate
