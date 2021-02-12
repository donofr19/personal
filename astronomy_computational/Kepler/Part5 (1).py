########################################################################
# Vincenzo Donofrio
# AST 304, Fall 2020
# Michigan State University
########################################################################

"""
This script calls functions necessary from ode.py and kepler.py to set the initial
conditions that can be used to compute the trajectory of an orbit and its subsequent
total, kinetic, and potential energy as a function of time for Forward Euler, 
Second-Order Runge-Kutta, and Fourth-Order Runge-Kutta numerical 
integration algorithms.
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm
# calling necessary functions
from ode import fEuler, rk2, rk4
from kepler import integrate_orbit, set_initial_conditions

# The integration methods from ode.py
integration_methods = {
    'Euler':fEuler,
    'RK2':rk2,
    'RK4':rk4
    }

# Use set_initial_conditions to get the initial conditions (position, velocity,
# energy-per-reduced mass, and the period for one full orbit)
init, energy_per_rm, period = set_initial_conditions(1, 1, .5)

h0 = 0.1*period   # Base stepsize

# Array of stepsizes (here we use only the smallest and largest from part 4:h and h/1024)
h_array = np.array([h0, h0/1024]) 

Euler_orbit = [] #initialize Euler list
RK2_orbit = [] #initialize RK2 list
RK4_orbit = [] #initialize RK4 list

# This loop uses integrate_orbit to solve the orbital equations of motion for 
# each of the two (largest and smallest) stepsizes in h_array, for all three of the integration methods
# in ode.py. Two arrays are collected for each integration method:
for i in range(len(h_array)):
    Euler_orbit.append(integrate_orbit(init, 1, 3*period, h_array[i], method = "Euler")) #Euler orbit
    RK2_orbit.append(integrate_orbit(init, 1, 3*period, h_array[i], method = "RK2")) #RK2 orbit
    RK4_orbit.append(integrate_orbit(init, 1, 3*period, h_array[i], method = "RK4")) #RK4 orbit
    

# Next, we want to assign all lists necessary for plotting with sensibly callable variables
# We need to first take care of the particle trajectory of each step size for each integration method (6 total).
# Assigning each component to a single variable is necessary to show correct interpretation of the return statement
# of the integrate_orbit function as well as the created arrays that contain a list of orbit information
# e.g. x_Euler_hlarge = Euler_orbit[0][1] corresponds to the x positions for the euler integration method with
# the large step size (h*1) used.


# Large step size - Euler
x_Euler_hlarge = Euler_orbit[0][1] # x positions for euler with large stepsize
y_Euler_hlarge = Euler_orbit[0][2] # y positions for euler with large stepsize

# Small step size - Euler
x_Euler_hsmall = Euler_orbit[1][1] # x positions for euler with small stepsize
y_Euler_hsmall = Euler_orbit[1][2] # y positions for euler with small stepsize

# Large step size - RK2
x_RK2_hlarge = RK2_orbit[0][1] # x positions for RK2 with large stepsize
y_RK2_hlarge = RK2_orbit[0][2] # y positions for RK2 with large stepsize
    
# Small step size - RK2
x_RK2_hsmall = RK2_orbit[1][1] # x positions for RK2 with small stepsize
y_RK2_hsmall = RK2_orbit[1][2] # y positions for RK2 with small stepsize

# Large step size - RK4
x_RK4_hlarge = RK4_orbit[0][1] # x positions for RK4 with large stepsize
y_RK4_hlarge = RK4_orbit[0][2] # y positions for RK4 with large stepsize

# Small step size - RK4
x_RK4_hsmall = RK4_orbit[1][1] # x positions for RK4 with small stepsize
y_RK4_hsmall = RK4_orbit[1][2] # y positions for RK4 with small stepsize



# We also need to collect the energies (KE, PE, TE) as a function of time (t)


# Time arrays *only need two (one for each stepsize) - we will utilize the 
#ones associated with Euler - could be RK2 or RK4, though*
time_large = Euler_orbit[0][0]
time_small = Euler_orbit[1][0]

# Assigning each component to a single variable is necessary to show correct interpretation of the return statement
# of the integrate_orbit function as well as the created arrays that contain a list of orbit information
# e.g. KE_Euler_hlarge = Euler_orbit[0][3] corresponds to the kinetic energy for the euler integration method with
# the large step size (h*1) used.

# Large step size - Euler
KE_Euler_hlarge = Euler_orbit[0][3] # kinetic energy for euler with large stepsize
PE_Euler_hlarge = Euler_orbit[0][4] # potential energy for euler with large stepsize
TE_Euler_hlarge = Euler_orbit[0][5] # total energy for euler with large stepsize

# Small step size - Euler
KE_Euler_hsmall = Euler_orbit[1][3] # kinetic energy for euler with small stepsize
PE_Euler_hsmall = Euler_orbit[1][4] # potential energy for euler with small stepsize
TE_Euler_hsmall = Euler_orbit[1][5] # total energy for euler with small stepsize

# Large step size - RK2
KE_RK2_hlarge = RK2_orbit[0][3] # kinetic energy for RK2 with large stepsize
PE_RK2_hlarge = RK2_orbit[0][4] # potential energy for RK2 with large stepsize
TE_RK2_hlarge = RK2_orbit[0][5] # Total energy for RK2 with large stepsize

# Small step size - RK2
KE_RK2_hsmall = RK2_orbit[1][3] # kinetic energy for RK2 with small stepsize
PE_RK2_hsmall = RK2_orbit[1][4] # potential energy for RK2 with small stepsize
TE_RK2_hsmall = RK2_orbit[1][5] # Total energy for RK2 with small stepsize

# Large step size - RK4
KE_RK4_hlarge = RK4_orbit[0][3] # kinetic energy for RK4 with large stepsize
PE_RK4_hlarge = RK4_orbit[0][4] # potential energy for RK4 with large stepsize
TE_RK4_hlarge = RK4_orbit[0][5] # Total energy for RK4 with large stepsize

# Small step size - RK4
KE_RK4_hsmall = RK4_orbit[1][3] # kinetic energy for RK4 with small stepsize
PE_RK4_hsmall = RK4_orbit[1][4] # potential energy for RK4 with small stepsize
TE_RK4_hsmall = RK4_orbit[1][5] # Total energy for RK4 with small stepsize


# Now we can plot each orbit and their corresponding energies


#This subplot represents the particle trajectory as the y positions
# against the x positions for each integration method and each stepsize.
# A 3x2 subplot is necessary to compare stepsize and integration method conveniently
# against themselves.

fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(20, 10))
fig.suptitle('Particle Trajectory', size = 20)
ax1.plot(x_Euler_hlarge, y_Euler_hlarge, label = 'Euler; h')
ax2.plot(x_Euler_hsmall, y_Euler_hsmall, 'tab:orange', label = 'Euler; h/1024')
ax3.plot(x_RK2_hlarge, y_RK2_hlarge, 'tab:green', label = 'RK2; h')
ax4.plot(x_RK2_hsmall, y_RK2_hsmall, 'tab:red', label = 'RK2; h/1024')
ax5.plot(x_RK4_hlarge, y_RK4_hlarge, 'tab:pink', label = 'RK4; h')
ax6.plot(x_RK4_hsmall, y_RK4_hsmall, 'tab:brown', label = 'RK4; h/1024')

fig.legend(prop={'size': 20})

ax1.set_ylabel('y (AU)', size = 15)
ax3.set_ylabel('y (AU)', size = 15)
ax5.set_ylabel('y (AU)', size = 15)
ax5.set_xlabel('x (AU)', size = 15)
ax6.set_xlabel('x (AU)', size = 15)

#This subplot represents the kinetic, potential, and total energy as a function 
#of time for each integration method and each stepsize. A 3x2 subplot is necessary 
#to compare stepsize and integration method conveniently against themselves.

fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(20, 10))
fig.suptitle('Energy as function of Time', size = 20)

ax1.plot(time_large, KE_Euler_hlarge, 'r', label = 'KE') 
ax1.plot(time_large, PE_Euler_hlarge, 'b', label = 'PE')
ax1.plot(time_large, TE_Euler_hlarge, 'g', label = 'TE')  

ax2.plot(time_small, KE_Euler_hsmall, 'r') 
ax2.plot(time_small, PE_Euler_hsmall, 'b')
ax2.plot(time_small, TE_Euler_hsmall, 'g')

ax3.plot(time_large, KE_RK2_hlarge, 'r') 
ax3.plot(time_large, PE_RK2_hlarge, 'b')
ax3.plot(time_large, TE_RK2_hlarge, 'g')

ax4.plot(time_small, KE_RK2_hsmall, 'r') 
ax4.plot(time_small, PE_RK2_hsmall, 'b')
ax4.plot(time_small, TE_RK2_hsmall, 'g')

ax5.plot(time_large, KE_RK4_hlarge, 'r') 
ax5.plot(time_large, PE_RK4_hlarge, 'b')
ax5.plot(time_large, TE_RK4_hlarge, 'g')

ax6.plot(time_small, KE_RK4_hsmall, 'r') 
ax6.plot(time_small, PE_RK4_hsmall, 'b')
ax6.plot(time_small, TE_RK4_hsmall, 'g')

ax1.set_title('Euler; h', size = 15)
ax2.set_title('Euler; h/1024', size = 15)
ax3.set_title('RK2; h', size = 15)
ax4.set_title('RK2; h/1024', size = 15)
ax5.set_title('RK4; h', size = 15)
ax6.set_title('RK4; h/1024', size = 15)

ax1.set_ylabel('Energy (in terms of GM/a)', size = 15)
ax3.set_ylabel('Energy (in terms of GM/a)', size = 15)
ax5.set_ylabel('Energy (in terms of GM/a)', size = 15)
ax5.set_xlabel('Time (2\u03C0 * yr)', size = 15)
ax6.set_xlabel('Time (2\u03C0 * yr)', size = 15)

fig.tight_layout()
fig.legend(prop={'size': 20})

