# PyDGSRD-1D
This is a Python code that solves 1D hyperbolic conservation laws on nonuniform grids using the state redistribution method.  State redistribution is an algorithm that solves the small cell problem on cut cell grids.  That is, arbitrarily small cells on embedded boundary grids result in overly restrictive maximum stable time steps when using explicit time stepping algorithms. Similar in spirit to flux redistribution by Collela [1], state redistribution relaxes this time step restriction using a simple postprocessing operation.  Of course, this algorithm is most interesting in two and three dimensions, but this one-dimensional code illustrates the important aspects of the algorithm.

<p align="center">
  <img src="https://github.com/andrewgiuliani/PyDGSRD/blob/main/images/srd.png" alt="SRD" width="300" >
</p>
<p align="center"> <i>High order approximation (p = 5) of an advecting pulse on a highly nonunform grid.  The DG solution on each element is plotted with a different colour.</i> <p align="center">
  
## Goal
The goal of this work is to use explicit Runge Kutta time steppers on highly nonuniform grids using the time step restriction 

<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\Delta&space;t&space;\leq&space;\frac{h}{a}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\Delta&space;t&space;\leq&space;\frac{1}{2p+1}\frac{h}{a}," title="  \Delta t \leq \frac{1}{2p+1} \frac{h}{a}" /></a>
</p>

where <a href="https://www.codecogs.com/eqnedit.php?latex=a" target="_blank"><img src="https://latex.codecogs.com/gif.latex?h" title="h" /></a> is the cell size, and <a href="https://www.codecogs.com/eqnedit.php?latex=a" target="_blank"><img src="https://latex.codecogs.com/gif.latex?a" title="a" /></a> is the maximum wavespeed in the numerical solution.   In the following grid, the majority of cells on the grid have size <a href="https://www.codecogs.com/eqnedit.php?latex=a" target="_blank"><img src="https://latex.codecogs.com/gif.latex?h" title="h" /></a> and only three have size <a href="https://www.codecogs.com/eqnedit.php?latex=\alpha&space;h" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha&space;h" title="\alpha h" /></a> with <a href="https://www.codecogs.com/eqnedit.php?latex=1&space;\leq&space;\alpha&space;<&space;1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?1&space;\leq&space;\alpha&space;<&space;1" title="1 \leq \alpha < 1" /></a>.  That is, <a href="https://www.codecogs.com/eqnedit.php?latex=\alpha" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha" title="\alpha" /></a> could potentially be quite small.
<p align="center">
  <img src="https://github.com/andrewgiuliani/PyDGSRD/blob/main/images/example.png" alt="example"  width="700">
</p>

State redistribution will allow the use of a time step that is proportional <a href="https://www.codecogs.com/eqnedit.php?latex=a" target="_blank"><img src="https://latex.codecogs.com/gif.latex?h" title="h" /></a> _not_ <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\alpha&space;h" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\alpha&space;h" title="\alpha h" /></a>.  The main idea of state redistribution is to temporarily merge, or coarsen, the numerical solution into neighborhoods located on the small cells.
Then, the solution on these neighborhoods is refined back onto the base, nonuniform, grid.  In the next section, we describe the different ways of determining these merging neighborhoods in one dimension (the preprocessing stage).

## Contact
For help running the code, or any other questions, send me an email at
`giuliani AT cims DOT nyu DOT edu`

## 📓&nbsp; License
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



## References

[1] Colella, Phillip, et al. "A Cartesian grid embedded boundary method for hyperbolic conservation laws." Journal of Computational Physics 211.1 (2006): 347-366.

[2] Giuliani, Andrew. "A two-dimensional stabilized discontinuous Galerkin method on curvilinear embedded boundary grids". [https://arxiv.org/abs/2102.01857](https://arxiv.org/abs/2102.01857)


