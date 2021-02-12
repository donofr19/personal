""" 
########################################################################
# Vincenzo Donofrio
# AST 304, Fall 2020
# Michigan State University
########################################################################
"""

'''
This routine computes the specific radius that satisfies
Lnuc(R) - 4piR**2sigmaTeff**4. We include the *for_main*
functions to allow for an easier time converting solar and 
scaled units.
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

from eos import get_rho_and_T, mean_molecular_weight
from ode import rk4
from astro_const import G, Msun, Rsun, Lsun, kB, m_u, fourpi, sigmaSB
from reactions import pp_rate
from zams import Teff, surface_luminosity
from structure import central_thermal, integrate
from structure_for_main import central_thermal_for_main, Teff_for_main

# Function to compute the luminosity difference: Lnuc(R) - 4piR**2sigmaTeff**4 
def lum_difference(radius,mass,delta_m,eta,xi,mue,pp_factor):
    """
    Computes the difference between nuclear luminosity
    and stellar luminosity.
    
    Arguments
       radius (scaled units)
       
       mass (scaled units)
       
       delta_m, eta, xi
           convergence paramters
       mue
           mean molecular weight 
       pp_factor
           multiplicative factor for rate
    Returns
       Lnuc(R) - 4piR**2sigmaTeff**4
    """
    m,r,p,Lnuc = integrate(mass,radius,delta_m,eta,xi,mue,pp_factor,max_steps=10000)
    return Lnuc[-1]-surface_luminosity(Teff_for_main(m[-1]),r[-1]) 
    
    
def find_radius(mass,delta_m,eta,xi,mue,pp_factor):
    """
    For a given mass calls rootfind over some range of radii,
    integrates over the function until the difference in luminosity is zero
    (nuclear luminosity = surface luminosity)
    
    Arguments
        mass (scaled units)

        delta_m, eta, xi
            convergence paramters  
        mue
            mean molecular weight  
        pp_factor
           multiplicative factor for rate
    Returns
        radius
            radius that satisfies the luminosity difference (scaled units)           
    """

    #range of radii; reason in detail under step 9 of report
    r_low = 0.01*Rsun # MKS
    r_high = 3*Rsun # MKS
    
    radius = brentq(lum_difference, r_low, r_high, xtol=1.0e-4, args = (mass,delta_m,eta,xi,mue,pp_factor))
    return radius

