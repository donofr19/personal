""" 
########################################################################
# Vincenzo Donofrio
# AST 304, Fall 2020
# Michigan State University
########################################################################
"""

'''
This script is a convergence test of our main file to make sure 
the convergence parameters are not too small but also occur
after covergence. 

Leaving one of the last three lines uncommented will test
a paramter's convergence; it is default set to delta_m
'''


# Necessary imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

from eos import get_rho_and_T, mean_molecular_weight
from ode import rk4
from astro_const import G, Msun, Rsun, Lsun, kB, m_u, fourpi, sigmaSB
from reactions import pp_rate
from zams import surface_luminosity
from structure import integrate
from rootfind import lum_difference, find_radius
from structure_for_main import central_thermal_for_main, Teff_for_main

# Using given pp_factor
pp_factor = 1.0

# Input constants
# Each array utilizes the assumption of current interstellar gas; so, arrays follow as
# [H,He,N]
Z = [1.0,2.0,7.0] 
A = [1.0,4.0,14.0] 
X = [0.706,.275,0.019] # Adds to 1

# Find mue
mue = mean_molecular_weight(Z, A, X)

# Setting mass *using just one should be fine*
mass_star = 0.2*Msun

# Initialize lists for L, Teff, rhoc, Tc for mass range
# For 10a
Luminosity_list = []
Teff_list = []

# For 10b
rhoc_list = []
Tc_list = []

# For convergence tests
# For delta_test
delta_test_array = np.array([1.0e-3,1.0e-4,1.0e-5,1.0e-6,1.0e-7,1.0e-8,1.0e-9])
eta = 1.0e-10
xi = 0.05

# For eta_test
delta_m = 1.0e-8*Msun
eta_test_array = np.array([1.0e-5,1.0e-6,1.0e-7,1.0e-8,1.0e-9,1.0e-10,1.0e-11])

# For xi_test
xi_test_array = np.array([1.0,0.5,0.1,0.05,0.01])

# convergence lists
mass_convergence_test = [] 
radius_convergence_test = [] 
pressure_convergence_test = [] 
luminosity_convergence_test = []
mass_keep = []

# Loop to find L, Teff, rhoc, and Tc for mass range
def delta_test(delta_test_array):
    """
    Computes convergence test for delta_m. Also includes a mass keep
    to see immediately which delta_m values correspond to convergence of mass
    
    Arguments
        delta_test_array
            array that contains testing delta_m values (input into function as solar units
            but later converted to MKS)
    Returns
        mass_keep
            array of delta_m to see immediately which values correspond to convergence of mass
        mass_convergence_test, radius_convergence_test, pressure_convergence_test, luminosity_convergence_test
            arrays of integrated values for parameter array
    """
    
    for i in delta_test_array:

        radius = find_radius(mass_star,i*Msun,eta,xi,mue,pp_factor) # converting to mks 
        Pc, rhoc, Tc = central_thermal_for_main(mass_star, radius, mue) 
        T_eff = Teff_for_main(mass_star)

        # Find integrated values; arrays
        int_mass, int_radius, int_pressure, int_Luminosity = integrate(mass_star,radius,i*Msun,eta,xi,mue,pp_factor,max_steps=10000)

        # Store in lists
        Luminosity_list.append(int_Luminosity[-1])
        Teff_list.append(T_eff)
        Tc_list.append(Tc)
        rhoc_list.append(rhoc) 

        # Store in convergence lists
        mass_convergence_test.append(int_mass[-1])
        radius_convergence_test.append(int_radius[-1])
        pressure_convergence_test.append(int_pressure[-1])
        luminosity_convergence_test.append(int_Luminosity[-1])
        
        # Save extremely senstive values
        if int_mass[-1]/mass_star >= .9999 and int_mass[-1]/mass_star <= 1.0001:
            mass_keep.append(i)
        
    return print('Converging values for delta_m:',mass_keep, 'Mass:',mass_convergence_test, 
                'Radius:',radius_convergence_test,'Pressure:',pressure_convergence_test, 'Luminosity:',luminosity_convergence_test)

def eta_test(eta_test_array):
    """
    Computes convergence test for eta. Also includes a mass keep
    to see immediately which eta values correspond to convergence of mass
    
    Arguments
        eta_test_array
            array that contains testing eta values
    Returns
        mass_keep
            array of eta to see immediately which values correspond to convergence of mass
        mass_convergence_test, radius_convergence_test, pressure_convergence_test, luminosity_convergence_test
            arrays of integrated values for parameter array
    """
    
    for i in eta_test_array:

        radius = find_radius(mass_star,delta_m,i,xi,mue,pp_factor) # converting to mks 
        Pc, rhoc, Tc = central_thermal_for_main(mass_star, radius, mue) 
        T_eff = Teff_for_main(mass_star)

        # Find integrated values; arrays
        int_mass, int_radius, int_pressure, int_Luminosity = integrate(mass_star,radius,delta_m,i,xi,mue,pp_factor,max_steps=10000)

        # Store in lists
        Luminosity_list.append(int_Luminosity[-1])
        Teff_list.append(T_eff)
        Tc_list.append(Tc)
        rhoc_list.append(rhoc) 

        # Store in convergence lists
        mass_convergence_test.append(int_mass[-1])
        radius_convergence_test.append(int_radius[-1])
        pressure_convergence_test.append(int_pressure[-1])
        luminosity_convergence_test.append(int_Luminosity[-1])
        
        # Save extremely senstive values       
        if int_mass[-1]/mass_star >= .9999 and int_mass[-1]/mass_star <= 1.0001:
            mass_keep.append(i)
        
    return print('Converging values for eta:',mass_keep, 'Mass:',mass_convergence_test, 
                'Radius:',radius_convergence_test,'Pressure:',pressure_convergence_test, 'Luminosity:',luminosity_convergence_test)
    
    
def xi_test(xi_test_array):
    """
    Computes convergence test for xi. Also includes a mass keep
    to see immediately which xi values correspond to convergence of mass
    
    Arguments
        xi_test_array
            array that contains testing xi values
    Returns
        mass_keep
            array of xi to see immediately which values correspond to convergence of mass
        mass_convergence_test, radius_convergence_test, pressure_convergence_test, luminosity_convergence_test
            arrays of integrated values for parameter array
    """
    
    for i in xi_test_array:

        radius = find_radius(mass_star,delta_m,eta,i,mue,pp_factor) # converting to mks 
        Pc, rhoc, Tc = central_thermal_for_main(mass_star, radius, mue) 
        T_eff = Teff_for_main(mass_star)

        # Find integrated values; arrays
        int_mass, int_radius, int_pressure, int_Luminosity = integrate(mass_star,radius,delta_m,eta,i,mue,pp_factor,max_steps=10000)

        # Store in lists
        Luminosity_list.append(int_Luminosity[-1])
        Teff_list.append(T_eff)
        Tc_list.append(Tc)
        rhoc_list.append(rhoc) 

        # Store in convergence lists
        mass_convergence_test.append(int_mass[-1])
        radius_convergence_test.append(int_radius[-1])
        pressure_convergence_test.append(int_pressure[-1])
        luminosity_convergence_test.append(int_Luminosity[-1])
        
        # Save extremely senstive values
        if int_mass[-1]/mass_star >= .9999 and int_mass[-1]/mass_star <= 1.0001:
            mass_keep.append(i)
        
    return print('Converging values for xi:',mass_keep, 'Mass:',mass_convergence_test, 
                'Radius:',radius_convergence_test,'Pressure:',pressure_convergence_test, 'Luminosity:',luminosity_convergence_test)

# *** Uncomment line for constant you want to test ***

delta_test(delta_test_array)
#eta_test(eta_test_array)
#xi_test(xi_test_array)

