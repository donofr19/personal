########################################################################
# Vincenzo Donofrio
# AST 304, Fall 2020
# Michigan State University
########################################################################

"""
  routines for advancing the solution of z'(t) = f(z,t) from t to t+h
  1. Forward Euler
  2. Second-order Runge-Kutta
  3. Fourth-order Runge-Kutta

  All routines have the same calling sequence:
  znew = <stepper>(f,t,z,h,args=())
  Arguments:
     f     := user-defined function; RHS of z'(t) = f(t,z,...)
     t     := current value of independent variable
     z     := current value of dependent variable
     h     := desired step
     args  := tuple of optional arguments to pass to f
"""

def fEuler(f,t,z,h,args=()):
    """
    Takes one forward Euler step.
    znew = fEuler(f,t,z,h,args=())
    
    Arguments:
         f     := user-defined function; RHS of z'(t) = f(t,z,...)
         t     := current value of independent variable
         z     := current value of dependent variable
         h     := desired step
         args  := tuple of optional arguments to pass to f
    Returns:
        znew   := z(t+h)
    """
    
    if not isinstance(args,tuple):
        args = (args,)

    return z + h*f(t,z,*args)

def rk2(f,t,z,h,args=()):
    """
    Takes one second-order Runge-Kutta step.
      znew = rk2(f,t,z,h,args=())
      Arguments:
         f     := user-defined function; RHS of z'(t) = f(t,z,...)
         t     := current value of indendent variable
         z     := current value of dependent variable
         h     := desired step
         args  := tuple of optional arguments to pass to f
      Returns:
        znew   := value of z at t+h
    """
    
    if not isinstance(args,tuple):
        args = (args,)
    
    zp = z + 0.5*h*f(t,z,*args)
    return z + h*f(t+0.5*h,zp,*args)

def rk4(f,t,z,h,args=()):
    """
    Takes one fourth-order Runge-Kutta step.
     znew = rk4(f,t,z,h,args=())
     Arguments:
        f     := user-defined function; RHS of z'(t) = f(t,z,...)
        t     := current value of indendent variable
        z     := current value of dependent variable
        h     := desired step
        args  := tuple of optional arguments to pass to f
      Returns:
        znew  := value of z at t+h
    """
   
    if not isinstance(args,tuple):
        args = (args,)
    
    k1 = f(t,z,*args)
    k2 = f(t+0.5*h,z+0.5*h*k1,*args)
    k3 = f(t+0.5*h,z+0.5*h*k2,*args)
    k4 = f(t+h,z+h*k3,*args)
    return z + h*(k1+2.0*k2+2.0*k3+k4)/6.0
