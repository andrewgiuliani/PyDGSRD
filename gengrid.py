import numpy as np
import ipdb
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-N", "--NUMELEM", type=int, help="number of elements", default = 0)
parser.add_argument("-L", "--LEFT", type=float, help="left endpoint", default = 0)
parser.add_argument("-R", "--RIGHT", type=float, help="right endpoint", default = 0)
parser.add_argument("-MESHTYPE", "--MESHTYPE", type=str, choices = ["uniform", "rand","perturb", "power", "paper", "bdry1", "bdry2"], default="load", help="mesh type")
parser.add_argument("-MERGETYPE", "--MERGETYPE", type=str, choices = ["LRP", "LRNP", "LP", "RP"], help="LRP = merge to left and right. LRNP = left, right, nonperiodic LP = left periodic. RP = right periodic")
parser.add_argument("-LOAD", "--LOAD", type = str, help="load a grid", default = "null")
args = parser.parse_args()
#ipdb.set_trace(context=21)

num_elem = args.NUMELEM 
left = args.LEFT
right = args.RIGHT
mesh_type = args.MESHTYPE
merge_type = args.MERGETYPE
loadgridname = args.LOAD

if mesh_type == "uniform":
    x = dx * np.linspace(left, right, N)
elif mesh_type == "rand": # grid points satisfy a uniform distribution
    temp = np.random.rand(num_elem+1)
    dx = (right-left)/num_elem
    temp[0] = 0
    x = np.cumsum(temp)
    x = left + (right-left) * x / x[-1]
elif mesh_type == "perturb":
    temp = np.linspace(0,1, num_elem+1)
    var = 1.
    perturb =  (np.random.rand(num_elem+1)-0.5)*2 * var /num_elem
    temp = temp + perturb
    temp = np.random.rand(num_elem+1)
    dx = (right-left)/num_elem
    temp[0] = 0
    x = np.cumsum(temp)
    x = left + (right-left) * x / x[-1]
elif mesh_type == "power": # grid points satisfy a power law distribution
    a = 0.25
    temp = np.random.power(a, num_elem + 1)
    dx = (right-left)/num_elem
    temp[0] = 0
    x = np.cumsum(temp)
    x = left + (right-left) * x / x[-1]
elif mesh_type == "paper":
    N = num_elem
    num_elem = 2*N + 1 + 2
    alpha = 1e-5 # volume fraction
    dx = (right-left) / (2*N + 1 + alpha * 2)
    reg = np.arange(N+1) * dx
    irreg = np.array( [-0.5 * dx, 0.5 * dx] )
    x = np.hstack( (reg + left, irreg, (alpha + 0.5) * dx + reg) )
elif mesh_type == "bdry1":
    alpha = 1e-5 # volume fraction
    dx = (right - left) / (num_elem-1 + alpha)
    reg = np.arange(num_elem) * dx + alpha*dx
    x = np.hstack( (np.array([0]), reg) ) + left
elif mesh_type == "bdry2":
    alpha = 1e-5 # volume fraction
    dx = (right - left) / (num_elem-2 + 2*alpha)
    reg = np.arange(num_elem-1) * dx + alpha*dx
    x = np.hstack( (np.array([left]), reg+left, np.array([right]) ) ) 
#    ipdb.set_trace(context=21)
elif mesh_type == "load":
    
    x = np.loadtxt(loadgridname + ".dat")
    with open(loadgridname + ".mdat", 'r') as file:
        indata = file.read().replace('\n', ' ').split()
    #ipdb.set_trace(context=21)
    dx = float(indata[0])
    merge_type = str(indata[1])
    #TOL = float(indata[2])
    grid_type = str(indata[3])
    num_elem = x.size-1



h = x[1:] - x[:-1]
print("Minimum volume fraction: " , np.min(h/dx) )
#print(h/dx)


TOL = dx/2

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











m = np.zeros(num_elem)
M = np.zeros(num_elem)
overlaps = np.ones(num_elem)
for elem in range(num_elem):
    left_length = 0
    right_length = 0

    m[elem] = elem
    M[elem] = elem

    if h[elem] > dx:
        continue

#    ipdb.set_trace(context=21)

    if merge_type == "LRNP": # merge to left and right (nonperiodic)
        curr = elem-1
        while left_length < TOL and curr > -1:
            left_length = left_length + h[curr]
            m[elem] = curr
            overlaps[curr] = overlaps[curr] + 1
            curr = curr - 1
        
        curr = elem+1
        while right_length < TOL and curr < num_elem:
            right_length = right_length + h[curr]
            M[elem] = curr
            overlaps[curr] = overlaps[curr] + 1
            curr = curr + 1

        # merge the left neighborhoods until a tolerance of TOL*2 now
        if left_length < TOL:
            curr = int(M[elem]+1)
            while right_length < 2*TOL and curr < num_elem:
                right_length = right_length + h[curr]
                M[elem] = curr
                overlaps[curr] = overlaps[curr] + 1
                curr = curr + 1
        
        # merge the right neighborhoods until a tolerance of TOL*2 now
        if right_length < TOL:
            curr = int(m[elem]-1)
            while left_length < 2*TOL and curr > -1:
                left_length = left_length + h[curr]
                m[elem] = curr
                overlaps[curr] = overlaps[curr] + 1
                curr = curr - 1

#        print(left_length/TOL, right_length/TOL)

    elif merge_type == "LRP": # merge to left and right (periodic)
        curr = (elem-1)%num_elem
        while left_length < TOL:
            left_length = left_length + h[curr]
            m[elem] = curr
            overlaps[curr] = overlaps[curr] + 1
            curr = (curr - 1)%num_elem
        
        curr = (elem+1)%num_elem
        while right_length < TOL:
            right_length = right_length + h[curr]
            M[elem] = curr
            overlaps[curr] = overlaps[curr] + 1
            curr = (curr + 1)%num_elem
    elif merge_type == "LP": # merge to left (periodic)
        curr = (elem-1)%num_elem
        while left_length < TOL:
            left_length = left_length + h[curr]
            m[elem] = curr
            overlaps[curr] = overlaps[curr] + 1
            curr = (curr - 1)%num_elem

    elif merge_type == "RP": # merge to right (periodic)
        curr = (elem+1)%num_elem
        while right_length < TOL:
            right_length = right_length + h[curr]
            M[elem] = curr
            overlaps[curr] = overlaps[curr] + 1
            curr = (curr + 1)%num_elem
       





f=open("grid_"+str(num_elem)+".mdat", 'w')
f.write(str(dx)+"\n")
f.write(str(merge_type)+"\n")
f.write(str(TOL)+"\n")
f.write(str(mesh_type)+"\n")
f.close()
np.savetxt("grid_"+str(num_elem)+".dat", x)

pdata = np.hstack( (m.reshape( (-1,1) ), M.reshape( (-1,1) ), overlaps.reshape( (-1,1) ) )  ).astype(int)
np.savetxt("grid_"+str(num_elem)+".pdat", pdata, fmt="%i")


