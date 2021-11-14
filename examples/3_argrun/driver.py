#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import pydgsrd1d as pdg

parser = argparse.ArgumentParser()
parser.add_argument("-P", "--POLYNOMIAL", type=int, help="polynomial degree of approximation")
parser.add_argument("-T", "--FINALTIME", type=float, help="final time")
parser.add_argument("-G", "--GRIDDNAME", type=str, help="grid name (without file extensions)")
parser.add_argument("-PLOT", "--PLOT", action="store_true", default=False)
args = parser.parse_args()


if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

p = args.POLYNOMIAL
T = args.FINALTIME
grid = args.GRIDDNAME
plotfinal = args.PLOT

l1, linf = pdg.PyDGSRD1D(p, grid, T, plotfinal)

