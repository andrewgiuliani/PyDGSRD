## 🧪 &nbsp; 2_run: launch the DG code with inputs in the driver.py 
This is a python code to solve a simple linear advection problem with periodic BCs.  The grid on which it is being solved was generated from

```
./gengrid.py -N 10 -L -1 -R 1 -MESHTYPE paper -MERGETYPE LRP
```

Run the `./driver.py` code to solve the linear advection problem on a nonuniform grid, stabilized by state redistribution.
