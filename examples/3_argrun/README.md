## 🧪 &nbsp; 3_argrun: launch the DG code with arguments from command line
This is a python code to solve a simple linear advection problem with periodic BCs with command line arguments.  Before trying this example, the user should generate a non-uniform grid using the `gengrid.py` utility in the `examples/1_gengrid/` directory.

First, a grid must be generated.  This can be done using 
```
./gengrid.py -L -1.0 -R 1.0 -N 100 -MESHTYPE bdry3 -MERGETYPE LRP
```
After the above is executed, `driver.py` can be called with command line arguments. For example,
```
./driver.py -P 5 -T 1 -G grid_204
```
which computes a sixth order (`p = 5`) approximation to the solution at the final time (`T = 1`) on `grid_100`.
The different arguments that `driver.py` accepts are explained below:

`-P`
polynomial degree of approximation

`-T`
final time

`-G`
grid filename (without any file extensions!)

`-PLOT`
plot the numerical solution at the final time using matplotlib

Some possible runs are provided below:
1. Generates the grid in the figures of this readme.  The small cell volume fraction, `alpha`, is set to 1e-5 in the code, but this can be modified.
```
./gengrid.py -L -1.0 -R 1.0 -N 100 -MESHTYPE bdry3 -MERGETYPE LRP
./driver.py -P 5 -T 1.0 -G grid_204
```
2. Reproduces the one-dimensional convergence test in [2], here, I've chosen a sixth order accurate numerical solution, but this can be changed.  Merging neighbourhoods are created by merging to the left and right of each small cell.
```
./gengrid.py -L -1.0 -R 1.0 -N 100 -MESHTYPE paper -MERGETYPE LRNP 
./driver.py -P 5 -T 1.0 -G grid_203
```
3. Generates a grid where the cell sizes follow a power law distribution.  This means that the grid is composed of vastly different cell sizes.
```
./gengrid.py -L -1.0 -R 1.0 -N 100 -MESHTYPE power -MERGETYPE LRP
./driver.py -P 5 -T 1.0 -G grid_100
```



