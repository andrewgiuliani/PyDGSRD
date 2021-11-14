# arg run
This is a python code to solve a simple linear advection problem with periodic BCs with command line arguments.  The grid on which it is being solved was generated from
python gengrid.py -N 100 -L -1 -R 1 -MESHTYPE paper -MERGETYPE LRP

The example can be run by calling
python driver.py -P 3 -T 1 -G grid_23
