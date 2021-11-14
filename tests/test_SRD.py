import pytest
import numpy as np
import pydgsrd1d as pdg



@pytest.mark.parametrize("p", [1, 2, 3, 4, 5])
def test_srd(p):
    T = 1
    plotfinal = False
    
    grid_list = ['data/grid_78', 'data/grid_88', 'data/grid_98', 'data/grid_108', 'data/grid_118', 'data/grid_128']
    l1_list = []
    linf_list =[]

    for grid in grid_list:
        l1, linf = pdg.PyDGSRD1D(p, grid, T, plotfinal)
        l1_list.append(l1)
        linf_list.append(linf)
    
    h = np.array([78, 88, 98, 108, 118, 128])
    
    A = np.concatenate( (np.ones((h.size,1)), np.log(1/h).reshape((-1,1))), axis=1)
    b_l1 = np.log(l1_list)
    b_linf = np.log(linf_list)
    
    ls_l1 = np.linalg.solve(A.T @ A, A.T @ b_l1)
    ls_linf = np.linalg.solve(A.T @ A, A.T @ b_linf) 
    
    min_rate = min([ls_l1[1], ls_linf[1]])
    print(min_rate)
    assert min_rate/(p+1) > 0.9
