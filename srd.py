from legendre import LegendreVandermonde
import numpy as np
import ipdb


def srd_basis_functions(grid_x, h, m, M, overlaps, p):
    gauss_legendre_quad = np.polynomial.legendre.leggauss(p+1)
    x_gl, w_gl = gauss_legendre_quad[0], gauss_legendre_quad[1]
    basis_vol,_ = LegendreVandermonde(x_gl, p+1)

    mbasis_vol = []
    mw_gl_list = []

    
    xval = grid_x[:-1, None] * (1.-x_gl[None, :])/2. + grid_x[1:, None] * (1. + x_gl[None, :])/2.
    
    for idx in range(h.size):
        
        if m[idx] > M[idx]:
            midxs_left = np.arange(m[idx], h.size)
            midxs_right = np.arange(0, M[idx]+1)
            
            mx_gl_left = xval[midxs_left, :].flatten()-grid_x[m[idx]]
            mx_gl_right = xval[midxs_right,:].flatten() -grid_x[0] + grid_x[-1]-grid_x[m[idx]]
            midxs = np.hstack( (midxs_left, midxs_right) )
            mx_gl = 2.*(np.hstack( (mx_gl_left, mx_gl_right) ) / np.sum( h[midxs] ) - 0.5)
            
        else:
            midxs = np.arange(m[idx], M[idx]+1)
            mx_gl = 2*(xval[midxs, :].flatten()-grid_x[m[idx]])/(np.sum(h[midxs])) - 1.
        
        mw_gl = np.outer(h[midxs]/overlaps[midxs], w_gl)
        vmerge = np.sum(np.sum(mw_gl,0))
        mw_gl = mw_gl / vmerge
        
        mbasis_vol_init,_ = LegendreVandermonde(mx_gl, p)
        wsqrt = np.sqrt(np.diagflat(mw_gl))
        q, r = np.linalg.qr(np.matmul(wsqrt,mbasis_vol_init))
        mbv = np.matmul( np.diagflat(mw_gl**(-0.5)) , q )
        if mbv[0,0] < 0 :
            mbv[:,0] = -mbv[:,0]
        mbasis_vol.append( mbv )
        mw_gl_list.append( mw_gl )
       
    return mbasis_vol, mw_gl_list
        
def srd(c_in, bv, mbv, m, M, overlaps, w_gl, mw_gl):
    qhat = np.zeros( c_in.shape )
    num_elem = c_in.shape[1]
    for idx in range(num_elem):
        
        if m[idx] > M[idx]:
            #ipdb.set_trace(context=21)
            midxs_left = np.arange(m[idx], num_elem)
            midxs_right = np.arange(0, M[idx]+1)
            midxs = np.hstack( (midxs_left, midxs_right) )
        else:
            midxs = np.arange(m[idx], M[idx]+1)
        
        # c_in * bv.T * 
        if m[idx] == M[idx]:
            qhat[:,idx, :] = np.matmul(c_in[:,idx,:], np.transpose(bv) )
        else:
            U = np.matmul(c_in[:,midxs,:], np.transpose(bv))
            Uw = (U * mw_gl[idx][None, :,:]).reshape( (c_in.shape[0], -1) ) #multiply each Uval by weights
            qhat[:,idx,:] = np.matmul(Uw, mbv[idx])
    
    c_out = np.zeros(c_in.shape)
    # project back onto the base grid
    for idx in range(num_elem):
        if m[idx] == M[idx]:
            c_out[:,idx,:] = c_out[:,idx,:]+ c_in[:,idx,:]/overlaps[None, idx, None]
        else:
            #ipdb.set_trace(context=21)
            if m[idx] > M[idx]:
                #ipdb.set_trace(context=21)
                midxs_left = np.arange(m[idx], num_elem)
                midxs_right = np.arange(0, M[idx]+1)
                midxs = np.hstack( (midxs_left, midxs_right) )
            else:
                midxs = np.arange(m[idx], M[idx]+1)

            Qhat = np.matmul(qhat[:,idx,:] , np.transpose(mbv[idx]) )           
            Qhat = Qhat.reshape((c_out.shape[0], midxs.size, -1))
            proj = 0.5*np.matmul(Qhat * w_gl[None,None,:], bv)
            c_out[:,midxs,:] = c_out[:,midxs,:] +  proj / overlaps[None, midxs, None]

    return c_out
