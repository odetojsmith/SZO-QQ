# SZO-LP (Safe Zeroth-Order Optimization Using Linear Programs)
Author: Baiwei Guo (baiwei.guo@epfl.ch) & Yang Wang (y.wang-40@tudelft.nl)

## Description
SZO-LP is a safe zeroth-order algorithm aimed at solving unmodelled optimization problems. In this repository, we use SZO-LP to solve an Optimal Power Flow problem.
Regarding the datails, we refer the readers to https://arxiv.org/abs/2304.01797 for our paper "Safe Zeroth-Order Optimization Using Linear Programs"

## To run this code
One need to load matpower in code/packages/matpower7.1 (If there is a new matpower version, please also rename the filefolder.)
The benchmark case is pglib_opf_case30_as.m, copied from "A Library of IEEE PES Power Grid Benchmarks", https://github.com/power-grid-lib/pglib-opf.
Then run the code in SZO-LP/code/SZO-LP/SZO_LP_OPF.m By doing this, one can derive a stationary point of the OPF problem.

