# this python script demonstrates the SRD approach to stabilizing the DG method in one dimension
import sys
import numpy as np
from rk import RK 
from legendre import LegendreVandermonde 
from advection import exact_solution, flux, numerical_flux, wavespeed 
#from burgers import exact_solution, flux, numerical_flux, wavespeed 
#from euler import exact_solution, flux, numerical_flux, wavespeed 
from plot import plot
from srd import srd_basis_functions, srd

import ipdb

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-P", "--POLYNOMIAL", type=int, help="polynomial degree of approximation")
parser.add_argument("-T", "--FINALTIME", type=float, help="final time")
parser.add_argument("-G", "--GRIDDNAME", type=str, help="grid name (without file extensions)")
parser.add_argument("-PLOT", "--PLOT", action="store_true", default=False)
args = parser.parse_args()




p = args.POLYNOMIAL
T = args.FINALTIME
grid = args.GRIDDNAME
plotfinal = args.PLOT

grid_name =           grid + ".dat"
preprocessing_name =  grid + ".pdat"
merging_meta_data  =  grid + ".mdat"
#ipdb.set_trace(context=21)

grid_x = np.loadtxt(grid_name) # load the element grid coordinates




### load quadrature information ###
gauss_legendre_quad = np.polynomial.legendre.leggauss(p+1)
x_gl, w_gl = gauss_legendre_quad[0], gauss_legendre_quad[1]
basis_vol,dbasis_vol = LegendreVandermonde(x_gl, p)
basis_edge,_         = LegendreVandermonde(np.array([-1, 1]), p)
basis_edge = np.transpose(basis_edge)
xval = grid_x[:-1, None] * (1.-x_gl[None, :])/2. + grid_x[1:, None] * (1. + x_gl[None, :])/2.
h = grid_x[1:] - grid_x[:-1]
N = h.size          # number of elements in the grid


## load merge data ##
with open(merging_meta_data, 'r') as file:
    indata = file.read().replace('\n', ' ').split()
#ipdb.set_trace(context=21)
dx = float(indata[0])
merge_type = str(indata[1])
TOL = float(indata[2])
grid_type = str(indata[3])


### load neighborhood data for SRD ###
preprocessing_data = np.loadtxt(preprocessing_name).astype(int)
m = preprocessing_data[:,0]  # starting index of merging neighborhood
M = preprocessing_data[:,1] # ending index of merging neighborhood
overlaps = preprocessing_data[:,2]        # overlap counts

mbasis_vol, mw_gl = srd_basis_functions(grid_x, h, m, M, overlaps, p)

def mass(c, h):
    return np.sum(c[:,:,0] * h[None,:], 1)

def compute_initial_condition(x): # return the modal coefficients on each element
    fvals = exact_solution(x, 0.)
    c = 0.5*np.matmul(fvals * w_gl[None,None,:], basis_vol)
    #ipdb.set_trace(context=21)
    #c = np.matmul(fvals * w_gl, basis_vol)
    return c
def compute_error(c,t):
#    ipdb.set_trace(context=21)
    gauss_legendre_quad_fine = np.polynomial.legendre.leggauss(p+2)
    x_gl_fine, w_gl_fine = gauss_legendre_quad_fine[0], gauss_legendre_quad_fine[1]
    basis_vol_fine,_ = LegendreVandermonde(x_gl_fine, p)
    x = grid_x[:-1, None] * (1.-x_gl_fine[None, :])/2. + grid_x[1:, None] * (1. + x_gl_fine[None, :])/2.

    U = np.matmul(c, np.transpose(basis_vol_fine))
    exact = exact_solution(x, t)
    diff = np.abs(U-exact)
    print("\n")
    print("######## Error at T = %.4e ######## "%t)
    for eq in range(c.shape[0]):
        L1 = np.sum( 0.5 * h * np.sum(diff[eq]*w_gl_fine, 1) )
        Linfinity = np.max(diff[eq])
        print("N = %i, L1[%i] error        %.16e"        % (c.shape[1],eq, L1))
        print("N = %i, Linfinity[%i] error %.16e" % (c.shape[1],eq, Linfinity))
    print("\n")

def surface(c): # surface integral
    # elem0 | elem1 | elem2
    #      f1      f2     
    Neq = c.shape[0]    
    U = np.matmul(c, basis_edge )
    Uminus = np.zeros( (Neq, N+1) )
    Uplus = np.zeros( (Neq, N+1) )
    
    # periodic boundary conditions
    Uminus[:,0]  = U[:,-1, 1] # last cell, right value
    Uplus[:,-1]  = U[:, 0, 0]  # first cell, left value
    
    Uminus[:, 1:] = U[:,:,1]
    Uplus[ :,:-1] = U[:,:,0]

    # numerical flux
    Fn = numerical_flux(Uminus, Uplus)
    surf = Fn[:,1:, None]*basis_edge[:,1] - Fn[:,:-1, None]*basis_edge[:,0]

    return surf 


def volume(c): # volume integral
    U = np.matmul(c, np.transpose(basis_vol))
    F = flux(U)
    vol = np.matmul( 0.5*F*w_gl , dbasis_vol)
    #vol = np.matmul( F*w_gl , dbasis_vol)
    return vol

def L(c, rk_time):
    rhs = -(surface(c) - volume(c))/h[None,:, None]
    return rhs

def compute_dt(c):
    U = np.matmul(c, np.transpose(basis_vol))
    speed = wavespeed(U)
    return  0.9*( dx / np.max(speed) ) / ( 2*p + 1 )


c0 = compute_initial_condition(xval)

m0 = mass(c0,h)
c0 = srd(c0, basis_vol, mbasis_vol, m, M, overlaps, w_gl, mw_gl)
m1 = mass(c0,h)

print("\n###### MERGING META DATA ######\n")
print("min volume fraction %.16e" % np.min(h/dx) )
if merge_type == "LRNP":
    print("merging type = left and right merging (not periodic)")
elif merge_type == "LRP":
    print("merging type = left and right merging (periodic)")
elif merge_type == "LP":
    print("merging type = left merging (periodic)")
elif merge_type == "RP":
    print("merging type = right merging (periodic)")
else:
    print("merging type not recognized...\n")
    quit()
print("merging tolerance %lf\n" % (TOL/dx) )
print("###### ----------------- #######\n\n")
print("###### MASS DIFFERENCE %.16e ######"%np.fabs(m0-m1))

max_step = -1
nstep = 0
time = 0


rk_b,rk_c,rk_d = RK(p) # load the coefficients from the Butcher tableau for the RK-(p+1) scheme
# storage of the RK stages
num_stages = rk_c.size
k = np.zeros( (c0.shape[0], c0.shape[1], c0.shape[2], num_stages) )

while (max_step == -1 and time < T) or (max_step > -1 and nstep < max_step) :
    dt = compute_dt(c0)
    if time + dt > T:
        dt = T - time
    for s in range(num_stages):
        ctemp = np.array(c0)

        # don't do this (or stabilize) for the first intermediate solution, U(1) = Un
        if s > 0: 
	        for ss in range(num_stages):
	            ctemp = ctemp + rk_b[s,ss] * k[:,:,:,ss]
        	ctemp =   srd(ctemp, basis_vol, mbasis_vol, m, M, overlaps, w_gl, mw_gl) 

        k[:,:,:,s] = dt * L(ctemp, time + dt * rk_d[s])

    for s in range(num_stages):
        c0 = c0 + rk_c[s] * k[:,:,:,s]
    c0 = srd(c0, basis_vol, mbasis_vol, m, M, overlaps, w_gl, mw_gl)

    time = time + dt
    nstep = nstep + 1

print("\n")
print("nstep = %i"%nstep)
compute_error(c0, time)

c0 = srd(c0, basis_vol, mbasis_vol, m, M, overlaps, w_gl, mw_gl)
m1 = mass(c0,h)

print("###### FINAL MASS DIFFERENCE %.16e ######"%np.fabs(m0-m1))

if plotfinal:
    plot(c0, grid_x, p)

