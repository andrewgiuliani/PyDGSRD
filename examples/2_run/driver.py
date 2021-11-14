#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import pydgsrd1d as pdg

p = 1
T = 1
grid = 'grid_23'
plotfinal = True

l1, linf = pdg.PyDGSRD1D(p, grid, T, plotfinal)
