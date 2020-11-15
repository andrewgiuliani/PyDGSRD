# 🛫&nbsp; PyDGSRD 
This is a Python code that solves 1D hyperbolic conservation laws on nonuniform grids using the state redistribution method.  State redistribution is an algorithm that solves the small cell problem on cut cell grids.  That is, arbitrarily small cells on embedded boundary grids result in overly restrictive maximum stable time steps when using explicit time stepping algorithms. Similar in spirit to flux redistribution by Collela [1], state redistribution relaxes this time step restriction using a simple postprocessing operation.  Of course, this algorithm is most interesting in two and three dimensions, but this one-dimensinal code illustrates the important aspects of the algorithm.

<p align="center">
  <img src="https://github.com/andrewgiuliani/PyDGSRD/blob/main/srd.png" alt="SRD" width="300" >
</p>
<p align="center"> <i>High order approximation (p = 5) of an advecting pulse on a highly nonunform grid.  The DG solution on each element is plotted with a different colour.</i> <p align="center">

## 🏗&nbsp; Grid generation and preprocessing 
Nonuniform grids on which state redistribution can be applied are generated using `gengrid.py`.  For example, the call

```
python gengrid.py -L -1.0 -R 1.0 -N 100 -MESHTYPE uniform -MERGETYPE LRP
```
randomly generates a nonuniform grid on the interval [-1.0,1.0] with 100 elements.  The cell sizes are generated using a uniform distribution, and merging neighborhoods are generated by merging to the left and right periodically.  The different arguments that `gengrid.py` accepts are explained below:


-N 
number of cells on the grid

-L
left endpoint

-R
right endpoint

-MESHTYPE
* `uniform`: cell sizes follow a uniform distribution.
* `perturb`: cells are generated by perturbing the endpoints of a uniform grid.
* `power`: cell sizes follow a power law distribution.
* `paper`: this is the grid used to demonstrate one-dimensional SRD in our paper.
* `bdry1`: one small cell on the left boundary, otherwise uniform.
* `brdy2`: one small cell on the left and right boundary, otherwise uniform.

-MERGETYPE
* `LRNP`: merge to the left, right, non-periodically.
* `LRP`:  merge to the left, right, periodically.
* `LP`: merge only to the left, periodically.
* `RP`: merge only to the right, periodically.
All of the above options use a specified tolerance (TOL) in the code.

After the grid generator finishes, it output three files.  The file with extension `.dat`, contains the grid endpoints.  The file with extension `.pdat` contains the preprocessing information that specifies the merging neighborhoods.

## 🏃🏻‍♀️&nbsp; Running the code 
[1] Colella, Phillip, et al. "A Cartesian grid embedded boundary method for hyperbolic conservation laws." Journal of Computational Physics 211.1 (2006): 347-366.
