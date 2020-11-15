import numpy as np

Neq = 1
a = 1 # advection velocity
def exact_solution(x,t):
    U = np.zeros( (Neq, x.shape[0], x.shape[1]) )
    #U[0,:,:] = 1 + (x-np.pi -t)**3
    #U[0,:,:] = x-a*t
    U[0,:,:] = np.sin(np.pi*(x-a*t))
    return U
def flux(U): # linear advection
    return a*U
def numerical_flux(Uleft, Uright): # upwinding
    return flux(Uleft)
def wavespeed(U):
    return a
