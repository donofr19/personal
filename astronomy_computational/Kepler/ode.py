########################################################################
# Vincenzo Donofrio
# AST 304, Fall 2020
# Michigan State University
########################################################################

"""
A module that containing functions that run the next step in sequence for the 
Forward Euler, Second-Order Runge-Kutta, and Fourth-Order Runge-Kutta numerical 
integration algorithms.
"""

def fEuler(f,t,z,h,args=()):
    """
    Forward Euler integration method. First computes the slope of f(t,z) and then
    use Taylor series to expand the solution to t+h. Repeats this process until it 
    comes to a solution.
    
    To call this function use the function name: fEuler
    
    Arguments
        f(t,z,...)
            function that contains the RHS of the equation dz/dt = f(t,z,...)
        t 
            independent variable
        z 
            dependent variable, function of t
        h 
            step size
        args (tuple, optional)
            additional arguments to pass to f
    
    Returns
        znew = z(t+h)
            new value of z from which the next step will be taken
    """
    
    # The following trick allows us to pass additional parameters to f
    # first we make sure that args is of type tuple; if not, we make it into
    # that form
    if not isinstance(args,tuple):
        args = (args,)
    
    # when we call f, we use *args to pass it as a list of parameters.
    # for example, if elsewhere we define f like
    # def f(t,z,x,y):
    #    ...
    # then we would call this routine as
    # znew = fEuler(f,t,z,h,args=(x,y))
    #
    return z + h*f(t,z,*args)


def rk2(f,t,z,h,args=()):
    """
    Second-order Runge-Kutta integration method. First the slope f(t,z) is computed
    and the solution to the midpoint t+h is calculated. Then the slope for the midpoint
    is calcuated and used to extend the solution.
    
    To call this function use the function name: rk2
    
    Arguments
        f(t,z,...)
            function that contains the RHS of the equation dz/dt = f(t,z,...)
        t 
            independent variable
        z 
            dependent variable, function of t
        h 
            step size
        args (tuple, optional)
            additional arguments to pass to f
    
    Returns
        znew = z(t+h)
            new value of z from which the next step will be taken
    """
    
    if not isinstance(args,tuple):
        args = (args,)
    
    
    z_mp = z + (h/2)*f(t,z,*args)
    t_half = t + (h/2)
    
    return z + h*f(t_half, z_mp, *args)

def rk4(f,t,z,h,args=()):
    """
    Fourth-Order Runge-Kutta integration method. Integrates from z to t to t+h in four
    steps. Forward step is taken to the midpoint and the solution is estimated. This solution
    is used to make an estimate of the new slope k2 and then k2 is used to take a step from t
    to the midpoint. Another value of the slope is computed k3 which is used to step across the
    entire interval. The slope at the endpoint t+h is computed. The final solution is the
    sum of the slopes. 
    
    To call this function use the function name: rk4
    
    Arguments
        f(t,z,...)
            function that contains the RHS of the equation dz/dt = f(t,z,...)
        t 
            independent variable
        z 
            dependent variable, function of t
        h 
            step size
        args (tuple, optional)
            additional arguments to pass to f
    
    Returns
        znew = z(t+h)
            new value of z from which the next step will be taken
    """
   
    if not isinstance(args,tuple):
        args = (args,)
    
    k1 = f(t,z,*args)
    k2 = f(t+(h/2), z+(h/2)*k1, *args)
    k3 = f(t+(h/2), z+(h/2)*k2, *args)
    k4 = f(t+h, z + h*k3, *args)
    
    return z + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
