########################################################################
# Vincenzo Donofrio
# AST 304, Fall 2020
# Michigan State University
########################################################################

"""
This script calls functions necessary from ode.py and kepler.py to set the initial
conditions that can be used to compute the total energy of an orbit which is then 
used to calculate the relative error in the Forward Euler, Second-Order Runge-Kutta, 
and Fourth-Order Runge-Kutta numerical integration algorithms as a function of error.
"""

# Necessary imports, including our routines from ode.py and kepler.py
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm
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

# Array of stepsizes for each run of integrations
h_array = np.array([h0, h0/2, h0/4, h0/8, h0/16, h0/32, h0/64, h0/128, h0/256, h0/512, h0/1024])   

# Initialize arrays to catch relative error calculations
relative_errors_euler = np.zeros(len(h_array))
relative_errors_rk2 = np.zeros(len(h_array))
relative_errors_rk4 = np.zeros(len(h_array))

# This loop uses integrate_orbit to solve the orbital equations of motion for 
# each of our stepsize in h_array, for all three of the intefration methods
# in ode.py. The relative errors in the total energy are then calculated and 
# collected in the relative_errors* arrays
for i in range(len(h_array)):
    te_euler = integrate_orbit(init, 1, 3*period, h_array[i], method = "Euler")[5]   # total energy of Euler method
    te_rk2 = integrate_orbit(init, 1, 3*period, h_array[i], method = "RK2")[5]   # total energy of RK2 method
    te_rk4 = integrate_orbit(init, 1, 3*period, h_array[i], method = "RK4")[5]   # total energy of RK4 method
    relative_errors_euler[i] = np.abs((te_euler[-1] - te_euler[0])/te_euler[0]) * 100   # relative error in TE of Euler
    relative_errors_rk2[i] = np.abs((te_rk2[-1] - te_rk2[0])/te_rk2[0]) * 100   # relative error in TE of RK2
    relative_errors_rk4[i] = np.abs((te_rk4[-1] - te_rk4[0])/te_rk4[0]) * 100   # relative error in TE of RK4

# Plots for relative error vs. stepsize for each integration method
plt.plot(h_array, relative_errors_euler, color = "blue", label = "Euler")
plt.plot(h_array, relative_errors_rk2, color = "orange", label = "RK2")
plt.plot(h_array, relative_errors_rk4, color = "red", label = "RK4")

plt.title("Relative Error vs. Stepsize", fontsize = 20, pad = 8)
plt.xlim(h0 + .01, h_array[-1])
plt.xlabel("Stepsize", fontsize = 12, labelpad = 8)
plt.ylabel("Relative Error", fontsize = 12, labelpad = 8)
plt.yscale("log")   # log scale on the y-axis to make the magnitude of the errors more apparent
plt.legend()
