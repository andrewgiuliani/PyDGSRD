## 🧪 &nbsp; 3_argrun: launch the DG code with arguments from command line
This is a python code to solve a simple linear advection problem with periodic BCs with command line arguments.

1. Generates the grid in the figures of this readme.  The small cell volume fraction, `alpha`, is set to 1e-5 in the code, but this can be modified.
```
python gengrid.py -L -1.0 -R 1.0 -N 100 -MESHTYPE bdry3 -MERGETYPE LRP
python driver.py -P 5 -T 1.0 -G grid_100
```
2. Reproduces the one-dimensional convergence test in [2], here, I've chosen a sixth order accurate numerical solution, but this can be changed.  Merging neighbourhoods are created by merging to the left and right of each small cell.
```
python gengrid.py -L -1.0 -R 1.0 -N 100 -MESHTYPE paper -MERGETYPE LRNP 
python driver.py -P 5 -T 1.0 -G grid_100
```
3. Generates a grid where the cell sizes follow a power law distribution.  This means that the grid is composed of vastly different cell sizes.
```
python gengrid.py -L -1.0 -R 1.0 -N 100 -MESHTYPE power -MERGETYPE LRP
python driver.py -P 5 -T 1.0 -G grid_100
```
