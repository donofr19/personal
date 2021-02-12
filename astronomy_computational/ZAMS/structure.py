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

We have also included pp_factor as a parameter in the non-testable routines
to allow for easier time changing pp_factors.
"""

import numpy as np
from eos import get_rho_and_T, mean_molecular_weight
from ode import rk4
from astro_const import G, Msun, Rsun, Lsun, kB, m_u, fourpi
from reactions import pp_rate
from structure_for_main import central_thermal_for_main


def central_thermal(m,r,mue):
    """
    Computes the central pressure, density, and temperature from the polytropic
    relations for n = 3/2.

    Arguments
        m
            mass in scaled units
        r
            radius is scaled units
        mu
            mean molecular weight
    Returns
        Pc, rhoc, Tc
            central pressure, density, and temperature in solar units
    """
    
    # Converting to MKS
    m = m*Msun
    r = r*Rsun
    
    # Computing Pc, rhoc, and Tc
    Pc = 0.77*G*m**2/r**4
    rhoc = 5.99*(3*m/(fourpi*r**3))
    Tc = 0.54*(mue*m_u/kB)*(G*m/r)

    return Pc, rhoc, Tc


def stellar_derivatives(m,z,mue,Pc,rhoc,Tc,pp_factor):
    """
    RHS of Lagrangian differential equations for radius and pressure
    
    Arguments
        m
            current value of the mass in scaled units
        z (array)
            current values of (radius, pressure) in scaled units
        mue
            ratio, nucleons to electrons.  For a carbon-oxygen white dwarf, 
            mue = 2.
        Pc
            central pressure in scaled units
        rhoc 
            central density in scaled units
        Tc
            central Temperature in scaled units
        pp_factor
            multiplicative factor for rate    
    Returns
        dzdm (array)
            Lagrangian derivatives dr/dm, dP/dm in scaled units
    """
    # Init. dzdm
    dzdm = np.zeros(3)
    
    # setting radius and pressure to arrays indices
    r = z[0]
    P = z[1]
    
    # getting rho and T
    rho, T = get_rho_and_T(P,Pc,rhoc,Tc) 
    
    # updating dzdm
    dzdm[0] = 1.0/(fourpi*r**2*rho)
    dzdm[1] = -G*m/fourpi/r**4
    dzdm[2] = pp_rate(T,rho,XH=0.7,pp_factor=pp_factor)
    
    return dzdm


def central_values(m,r,delta_m,mue,pp_factor):
    """
    Constructs the boundary conditions at the edge of a small, constant density 
    core of mass delta_m with central pressure P_c
    
    Arguments
        m 
            mass in scaled units
        r 
            radius in scaled units
        Pc
            central pressure in scaled units
        delta_m
            core mass in scaled units
        mue
            nucleon/electron ratio
        pp_factor
            multiplicative factor for rate    
    Returns
        z = array([ r, p, l ])
            central values of radius, pressure, and luminosity in scaled units
    """
    # Init. dzdm    
    z = np.zeros(3)
    
    # Gettinf Pc, rhoc, and Tc
    Pc, rhoc, Tc = central_thermal_for_main(m,r,mue)
    
    # Updating z
    z[0] = (3.0*delta_m/fourpi/rhoc)**(1/3) 
    z[1] = Pc
    z[2] = pp_rate(Tc,rhoc,XH=0.7,pp_factor=pp_factor)*delta_m # boundary condition

    return z

def lengthscales(m,z,mue,Pc,rhoc,Tc,pp_factor):
    """
    Computes the radial length scale H_r and the pressure length H_P
    
    Arguments
        m
           current mass coordinate in scaled units
        z (array)
           [ r, p, l ] in scaled units
        mue
            mean electron weight
        Pc, rhoc, Tc
            central pressure, density, and temperature in solar units         
        pp_factor
            multiplicative factor for rate 
    Returns
        z/|dzdm| in units of solar masses
    """

    # add small number to dz/dm to prevent division by zero
    EPS = 1.0e-30
    return z/(np.abs(stellar_derivatives(m,z,mue,Pc,rhoc,Tc,pp_factor))+EPS) # z includes HL for L/dLdm


def integrate(mass,r,delta_m,eta,xi,mue,pp_factor,max_steps=10000):
    """
    Integrates the scaled stellar structure equations
    Arguments
        mass
            mass in scaled units; *using instead of m to differ from variable delta_m is later stored in*
        r
            radius in scaled units
        delta_m
            initial offset from center in units of solar mass
        eta
            The integration stops when P < eta * Pc
        xi
            The stepsize is set to be xi*min(p/|dp/dm|, r/|dr/dm|)
        mue
            mean electron mass
        pp_factor
            multiplicative factor for rate 
        max_steps
            solver will quit and throw error if this more than max_steps are 
            required (default is 10000)
                       
    Returns
        m_step, r_step, p_step, l_step
            arrays containing mass coordinates, radii, pressures, and luminosies during 
            integration in scaled units
    """
    
    # allocate storage
    m_step = np.zeros(max_steps)
    r_step = np.zeros(max_steps)
    p_step = np.zeros(max_steps)
    l_step = np.zeros(max_steps)
    
    # set starting conditions using central values
    z = central_values(mass,r,delta_m,mue,pp_factor)
    Pc, rhoc, Tc = central_thermal_for_main(mass,r,mue)

    m = delta_m 
    Nsteps = 0
    for step in range(max_steps):
        radius = z[0]
        pressure = z[1]
        luminosity = z[2]
        # are we at the surface?
        if (pressure < eta*Pc):
            break
        # store the step
        m_step[step] = m
        r_step[step] = radius
        p_step[step] = pressure
        l_step[step] = luminosity 
         
        # set the stepsize
        h = xi*min(lengthscales(m,z,mue,Pc,rhoc,Tc,pp_factor))

        # take a step
        z = rk4(stellar_derivatives,m,z,h,args=(mue,Pc,rhoc,Tc,pp_factor))
        m += h
        # increment the counter
        Nsteps += 1
    # if the loop runs to max_steps, then signal an error
    else:
        raise Exception('too many iterations')
        
    return m_step[0:Nsteps],r_step[0:Nsteps],p_step[0:Nsteps],l_step[0:Nsteps]
