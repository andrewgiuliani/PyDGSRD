import numpy as np
import ipdb
import sys

N = int(sys.argv[1]) # N large cells, 2 small cells
num_elem = 2*N + 1 + 2
alpha = 1e-5 # volume fraction
left = -1
right = 1
dx = (right-left) / (2*N + 1 + alpha * 2)

reg = np.arange(N+1) * dx
irreg = np.array( [-0.5 * dx, 0.5 * dx] )
x = np.hstack( (reg + left, irreg, (alpha + 0.5) * dx + reg) )


h = x[1:]-x[:-1]
overlaps = np.ones(num_elem)
m = np.arange(num_elem)
M = np.arange(num_elem)

neigh_idx1 = N
neigh_idx2 = N+2

m[neigh_idx1] = neigh_idx1-1
M[neigh_idx1] = neigh_idx1+1

m[neigh_idx2] = neigh_idx2-1
M[neigh_idx2] = neigh_idx2+1

overlaps[neigh_idx1-1] = 2
overlaps[N+1] = 3
overlaps[neigh_idx2+1] = 2

pdata = np.hstack( (m.reshape( (-1,1) ), M.reshape( (-1,1) ), overlaps.reshape( (-1,1) ) )  ).astype(int)




f=open("grid_"+str(num_elem)+".dat", 'w')
f.write(str(dx)+"\n")
np.savetxt(f, x)
f.close()
np.savetxt("grid_"+str(N)+".pdat", pdata, fmt="%i")

#ipdb.set_trace(context=21)
