""" 
########################################################################
# Vincenzo Donofrio
# AST 304, Fall 2020
# Michigan State University
########################################################################
"""

""" 
Routines for computing structure of fully convective star of a given mass and 
radius.

Through some trial and error, we have found that changing the unit in which Teff and central_thermal handles 
the star's mass allows for a smoother process to compute the missing radius in the main scripts. These versions
of Teff and central_thermal will be used in the main scripts, while the original versions used in structure.py
are fitted so the grader can check them with the given tests.

We have also included pp_factor as a parameter in the non-testable routines
to allow for easier time changing pp_factors.
"""

import numpy as np
from eos import get_rho_and_T, mean_molecular_weight
from ode import rk4
from astro_const import G, Msun, Rsun, Lsun, kB, m_u, fourpi
from reactions import pp_rate


def central_thermal_for_main(m,r,mue):
    """
    Computes the central pressure, density, and temperature from the polytropic
    relations for n = 3/2.

    Arguments
        m
            mass in scaled units
        r
            radius is scaled units
        mue
            mean molecular weight
    Returns
        Pc, rhoc, Tc
            central pressure, density, and temperature in solar units
    """
    
    # Not converting to scaled units as done in central_thermal under structure.py
    # Computing Pc, rhoc, and Tc
    Pc = 0.77*G*m**2/r**4
    rhoc = 5.99*(3*m/(fourpi*r**3))
    Tc = 0.54*(mue*m_u/kB)*(G*m/r)

    return Pc, rhoc, Tc

def Teff_for_main(Mwant):
    """
    Interpolates effective temperatures given mass from tabulated [1] values 
    for low-mass stars.
    [1] Chabrier, Baraffe, Allard, and Hauschildt. Evolutionary Models for Very 
    Low-Mass Stars and Brown Dwarfs with Dusty Atmospheres. Astrophys. Jour. 
    542:464--472, Oct. 2000.
    Parameters
        Mwant (float, scalar or array)
            Mass of the star(s) in units of solar masses
    Returns
       Teff (float, same type/shape as Mwant)
            Interpolated effective temperatures of the stars in Kelvin.
    """

    # tabulated values from Chabrier et al. (2000)
    masses = np.array([0.1*Msun,0.15*Msun,0.2*Msun,0.3*Msun]) # scaled units
    Teffs = np.array([2800.0,3150.0,3300.0,3400.0]) # [K]
    
    # fill this out to perform interpolation to find Teff for Mwant
    T_eff =  np.interp(Mwant,masses,Teffs)
    return T_eff

