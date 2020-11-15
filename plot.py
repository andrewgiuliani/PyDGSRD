import matplotlib.pyplot as plt
import numpy as np
from legendre import LegendreVandermonde
from matplotlib.collections import LineCollection
import ipdb

def plot(c, grid_x, p):



    res = 10
    xi = np.linspace(-1, 1, res)
    basis_vol,_ = LegendreVandermonde(xi, p)
    xvals = grid_x[:-1, None] * (1.-xi[None, :])/2. + grid_x[1:, None] * (1. + xi[None, :])/2. 





    # (Neq, N, res, 2)
    Ux = np.zeros( (c.shape[0], c.shape[1], res, 2) )
    Ux[:,:,:,0] = xvals
    Ux[:,:,:,1] = np.matmul(c, np.transpose(basis_vol)) 

    #ipdb.set_trace(context=21)


    fig, ax = plt.subplots(1, c.shape[0])
    if c.shape[0] > 1:
        for eq in range(c.shape[0]):
            ax[eq].set_xlim(grid_x.min(), grid_x.max())
            ax[eq].set_ylim(Ux[eq,:,:,:].min(), Ux[eq,:,:,:].max())



            line_segments = LineCollection(Ux[eq,:,:,:])
            ax[eq].add_collection(line_segments)
    else:
        ax.set_xlim(grid_x.min(), grid_x.max())
        ax.set_ylim(Ux[0,:,:,:].min(), Ux[0,:,:,:].max())

        line_segments = LineCollection(Ux[0,:,:,:])
        ax.add_collection(line_segments)









    plt.show()
