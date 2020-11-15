


# PyDGSRD
This is a Python code that solves 1D hyperbolic conservation laws on nonuniform grids using the state redistribution method.  State redistribution is an algorithm that solves the small cell problem on cut cell grids.  That is, arbitrarily small cells on embedded boundary grids result in overly restrictive maximum stable time steps when using explicit time stepping algorithms. Similar in spirit to flux redistribution by Collela, state redistribution relaxes this time step restriction using a simple postprocessing operation.  Of course, this algorithm is most interesting in two and three dimensions, but this one-dimensinal code illustrates the important aspects of the algorithm.

<p align="center">
  <img src="https://github.com/andrewgiuliani/PyDGSRD/blob/main/srd.png" alt="SRD" width="300" >
</p>
<p align="center"> <i>High order approximation of an advecting pulse on a highly nonunform grid.  The solutions on different elements is plotted with different colours.</i> <p align="center">



### Grid generation
Nonuniform grids on which state redistribution can be applied are generated using `gengrid.py`.

### Mesh preprocessing
Preprocessing of the generated grid is done automatically in `gengrid.py`.

### Running the code
