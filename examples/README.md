## 🧪 &nbsp; Examples

This folder contains three examples to learn how to use the code:
- `1_gengrid` shows how to generate nonuniform grids in 1D on which the stabilized DG method can be used.
- `2_run` shows how to run the DG code from a driver file, where inputs to the DG code are provided in the driver code.
- `3_argrun` shows how to run the DG code, where inputs to the DG code are provided at command line.

In example `1_gengrid`, we generate nonuniform grids. Then in examples `2_run` and `3_argrun`, we use a stabilized DG method to solve the following simple linear advection problem on nonuniform grids generated using the tool in `1_gengrid`:

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{align*}&space;\begin{aligned}&space;u_t&space;&plus;&space;u_x&space;&=&space;0&space;\\&space;u(x,0)&space;&=&space;\cos&space;\left(\pi&space;x&plus;\frac{\pi}{3}&space;\right)\\&space;u(-1,t)&space;&=&space;u(1,t)&space;\end{aligned}&space;\end{align*}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{align*}&space;\begin{aligned}&space;u_t&space;&plus;&space;u_x&space;&=&space;0&space;\\&space;u(x,0)&space;&=&space;\cos&space;\left(\pi&space;x&plus;\frac{\pi}{3}&space;\right)\\&space;u(-1,t)&space;&=&space;u(1,t)&space;\end{aligned}&space;\end{align*}" title="\begin{align*} \begin{aligned} u_t + u_x &= 0 \\ u(x,0) &= \cos \left(\pi x+\frac{\pi}{3} \right)\\ u(-1,t) &= u(1,t) \end{aligned} \end{align*}" /></a>

The initial condition and advection speed can be set in `pydgsrd1d/advection.py`.
