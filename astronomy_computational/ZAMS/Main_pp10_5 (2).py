""" 
########################################################################
# Vincenzo Donofrio
# AST 304, Fall 2020
# Michigan State University
########################################################################
"""

'''
This script involves the main execution of the code. Included are the 
calculations for log(L/Lsun) against log(Teff/K) and log(Tc/K) 
against log(rhoc/g cm^-3). In addition, are calculations for 
T(r), T(m), L(r), and L(m) for a 0.3Msun star. We also 
find the radius at which L(r) is at 90% and the fraction of
stellar mass that is enclosed by this radius. Finally,
we include code to plot under namely titles. This script is
for a pp_factor of 10^5
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


pp_factor = 10e+5

# Input constants
eta = 1.0e-10
delta_m = 1.0e-8*Msun # convert to MKS
xi = 0.05
# Each array utilizes the assumption of current interstellar gas; so, arrays follow as
# [H,He,N]
Z = [1.0,2.0,7.0] 
A = [1.0,4.0,14.0] 
X = [0.706,.275,0.019] # Adds to 1

# Find mue
mue = mean_molecular_weight(Z, A, X)

# Setting mass range (0.1Msun <= M <= 0.3Msun)
mass_range = np.linspace(0.1, 0.3, num=25)

# Initialize lists for L, Teff, rhoc, Tc for mass range
# For 10a
Luminosity_list = []
Teff_list = []

# For 10b
rhoc_list = []
Tc_list = []

# Initialize list for T (for T(r) and T(m)):
# For 11
adiabatic_T = []

# Loop to find L, Teff, rhoc, and Tc for mass range
for i in mass_range:
    radius = find_radius(i*Msun,delta_m,eta,xi,mue,pp_factor) # converting to mks 
    Pc, rhoc, Tc = central_thermal_for_main(i*Msun, radius, mue) 
    T_eff = Teff_for_main(i*Msun)

    # Find integrated values; arrays
    int_mass, int_radius, int_pressure, int_Luminosity = integrate(i*Msun,radius,delta_m,eta,xi,mue,pp_factor,max_steps=10000)

    # Store in lists
    Luminosity_list.append(int_Luminosity[-1])
    Teff_list.append(T_eff)
    Tc_list.append(Tc)
    rhoc_list.append(rhoc) 



# Finding T(r), T(m), L(r), and L(m) for M = 0.3Msun; using eq. 13 of instruction to compute adiabatic 
# relation for temperature, or T = Tc(P/Pc)**(1-1/gamma)    
radius_dot3 = find_radius(0.3*Msun,delta_m,eta,xi,mue,pp_factor) # converting to mks
Pc_dot3, rhoc_dot3, Tc_dot3 = central_thermal_for_main(0.3*Msun, radius, mue) 

# Find integrated values; arrays
int_mass_dot3, int_radius_dot3, int_pressure_dot3, int_Luminosity_dot3 = integrate(0.3*Msun,radius,delta_m,eta,xi,mue,pp_factor,max_steps=10000)

# Loop to find adiabatic_T's array of values
for i in int_pressure_dot3:
    gamma = 5/3
    T = Tc_dot3*(i/Pc_dot3)**(1-(1/gamma))
    adiabatic_T.append(T)

# Find what radius L(r) reaches 90% of final value and what fraction of star's mass is enclosed by the radius
for i in range(len(int_Luminosity_dot3)):
            if int_Luminosity_dot3[i]/int_Luminosity_dot3[-1] >= 90/100:

                radius_90 = int_radius_dot3[i] # 90% radius
                mass_fraction_90 = int_mass_dot3[i]/int_mass_dot3[-1] # mass fraction
                break # so loop stops at the 90% value (first value), 
                #instead of printing out the 100% value (end value)


print('Radius when L(r) reaches 90%:',radius_90,'meters')
print(mass_fraction_90,'of the stars mass is enclosed by this radius')

# Plots for Part 10 a and b
fig, ax = plt.subplots(1,2,figsize=(10,5))
fig.subplots_adjust(wspace=0.6)
fig.tight_layout(pad=3.0)

# Plots log(L/Lsun) against log(Teff/K)
ax[0].scatter(np.log10(Teff_list),np.log10([x / Lsun for x in Luminosity_list]),color='darkcyan',marker='^') # converting to solar units
ax[0].invert_xaxis() # HR
ax[0].set_xlabel('log(Teff/K)')
ax[0].set_ylabel('log(L/Lsun)')
ax[0].set_title('Mass Range of 0.1MSun-0.3MSun; pp_factor=10^5')

# Plots log(Tc/K) against log(rhoc/g cm^-3)
ax[1].scatter(np.log10([(x*1000/100**3) for x in rhoc_list]),np.log10(Tc_list),color='darkcyan',marker='^') # converting to CGS
ax[1].set_xlabel('log(rhoc/g cm^-3)')
ax[1].set_ylabel('log(Tc/K)')
ax[1].set_title('Mass Range of 0.1MSun-0.3MSun; pp_factor=10^5')

fig.savefig('Part10_pp10_5.pdf',format='pdf',bbox_inches='tight',facecolor='white')
    
# Plots for Part 11

fig, ax = plt.subplots(2,2,figsize=(10,5))
fig.subplots_adjust(wspace=0.6)
fig.tight_layout(pad=3.0)

# Plots T(r) 
ax[1,0].plot(int_radius,adiabatic_T,color='darkcyan', linewidth=4) 
ax[1,0].set_xlabel('radius (m)')
ax[1,0].set_ylabel('T(r) (K)')
ax[1,0].set_title('T(r); 0.3Msun; pp_factor=10^5')

# Plots L(r) 
ax[1,1].plot(int_radius,int_Luminosity,color='orange', linewidth=4) 
ax[1,1].set_xlabel('radius (m)')
ax[1,1].set_ylabel('L(r) (W)')
ax[1,1].set_title('L(r); 0.3Msun; pp_factor=10^5')

# Plots T(m) 
ax[0,0].plot(int_mass,adiabatic_T,color='black', linewidth=4) 
ax[0,0].set_xlabel('mass (kg)')
ax[0,0].set_ylabel('T(m) (K)')
ax[0,0].set_title('T(m); 0.3Msun; pp_factor=10^5')

# Plots L(m)
ax[0,1].plot(int_mass,int_Luminosity,color='blue', linewidth=4) 
ax[0,1].set_xlabel('mass (kg)')
ax[0,1].set_ylabel('L(m) (W)')
ax[0,1].set_title('L(m); 0.3Msun; pp_factor=10^5')

fig.savefig('Part11_pp10_5.pdf',format='pdf',bbox_inches='tight',facecolor='white')



