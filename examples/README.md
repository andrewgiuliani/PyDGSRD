## 🏗&nbsp; Grid generation and preprocessing 
Nonuniform grids on which state redistribution can be applied are generated using `gengrid.py`.  For example, the call

```
python gengrid.py -L -1.0 -R 1.0 -N 100 -MESHTYPE rand -MERGETYPE LRP
```
randomly generates a nonuniform grid on the interval [-1.0,1.0] with 100 elements.  The cell sizes are generated using a uniform distribution, and merging neighborhoods are generated by merging to the left and right periodically.  The different arguments that `gengrid.py` accepts are explained below:


-N 
number of cells on the grid

-L
left endpoint

-R
right endpoint

-MESHTYPE
* `uniform`: uniform grid.
* `rand`: cell sizes follow a uniform distribution.
* `perturb`: cells are generated by perturbing the endpoints of a uniform grid.
* `power`: cell sizes follow a power law distribution.
* `paper`: this is the grid used to demonstrate one-dimensional SRD in our paper.
* `bdry1`: one small cell on the left boundary, otherwise uniform.
* `brdy2`: one small cell on the left and right boundary, otherwise uniform.
* `brdy3`: one small cell on the left, and two small cells at the center, otherwise uniform.

-MERGETYPE

* `LRNP`: merge to the left, right, non-periodically.
<p align="center">
  <img src="https://github.com/andrewgiuliani/PyDGSRD/blob/main/images/LRPNP.png" alt="mergetype"  width="700">
</p>
<p align="center"> <i>The LRP option merges a small cell until the neighborhood has size TOL on the left and the right of the small cell, with size alpha h, only if it can.</i> <p align="center">


* `LRP`:  merge to the left, right, periodically.
<p align="center">
  <img src="https://github.com/andrewgiuliani/PyDGSRD/blob/main/images/LRP.png" alt="mergetype"  width="700" >
</p>
<p align="center"> <i>The LRP option merges a small cell until the neighborhood has size TOL on the left and the right of the small cell, with size alpha h.</i> <p align="center">


* `LP`: merge only to the left, periodically.
<p align="center">
  <img src="https://github.com/andrewgiuliani/PyDGSRD/blob/main/images/LP.png" alt="mergetype"  width="700" >
</p>
<p align="center"> <i>The LP option merges a small cell until the neighborhood has size TOL only to the left of the small cell, with size alpha h.</i> <p align="center">


* `RP`: merge only to the right, periodically.
<p align="center">
  <img src="https://github.com/andrewgiuliani/PyDGSRD/blob/main/images/RP.png" alt="mergetype"  width="700" >
</p>
<p align="center"> <i>The LRP option merges a small cell until the neighborhood has size TOL only to the right of the small cell, with size alpha h.</i> <p align="center">

All of the above options create neighbourhoods using a specified tolerance (TOL) in the code, whereby cells are merged to the left or to the right until the neighborhood satisfies a size constraint:



After the grid generator finishes, it output three files.  
- The file with extension `.dat`, contains the `N+1` grid endpoints.  
- The file with extension `.pdat` contains the preprocessing information that specifies the merging neighborhoods.  The file contains three columns corresponding to `m`, `M`, and `overlaps` in the code.  `m` and `M` specify the indices of the first and last cell in the merging neighbourhoods.  For example the merging neighborhood associated to cell 5 is made up of cells with indices `m[5], m[5]+1, ..., M[5]`.  `overlaps` contains the number of neighborhoods that overlap each cell in the grid.
- The file with extension `.mdat` contains metadata about the preprocessing stage: cell size to be used in the CFL condition, which merging algorithm was chosen, the merging tolerance, and type of random grid that was generated.

## Running the code 
`PyDGSRD` can be called after grid generation with, for example,
```
python PyDGSRD.py -P 5 -T 1 -G grid_100
```
which computes a sixth order (`p = 5`) approximation to the solution at the final time (`T = 1`) on `grid_100`.
The different arguments that `PyDGSRD.py` accepts are explained below:

-P
polynomial degree of approximation

-T
final time

-G
grid filename (without any file extensions!)

-PLOT
plot the numerical solution at the final time using matplotlib




## 🧪 &nbsp; Examples

1. Generates the grid in the figures of this readme.  The small cell volume fraction, `alpha`, is set to 1e-5 in the code, but this can be modified.
```
python gengrid.py -L -1.0 -R 1.0 -N 100 -MESHTYPE bdry3 -MERGETYPE LRP
python PyDGSRD.py -P 5 -T 1.0 -G grid_100
```
2. Reproduces the one-dimensional convergence test in [2], here, I've chosen a sixth order accurate numerical solution, but this can be changed.  Merging neighbourhoods are created by merging to the left and right of each small cell.
```
python gengrid.py -L -1.0 -R 1.0 -N 100 -MESHTYPE paper -MERGETYPE LRNP 
python PyDGSRD.py -P 5 -T 1.0 -G grid_100
```
3. Generates a grid where the cell sizes follow a power law distribution.  This means that the grid is composed of vastly different cell sizes.
```
python gengrid.py -L -1.0 -R 1.0 -N 100 -MESHTYPE power -MERGETYPE LRP
python PyDGSRD.py -P 5 -T 1.0 -G grid_100
```

