# Computing-the-Exact-LogConcave-MLE
This repository contains computations for the paper "Exact solutions in low-rank approximation with zeros" by Kaie Kubjas, Luca Sodomaco and Elias Tsigaridas. 

## Macaulay2 code: orbits.m2
This file contains the code that computes the orbit representatives and orbit sizes of zero patterns in Table 7.

## Julia code: crit_points_nxn_corank1_one_zero.jl
This file contains the code that computes the critical points of constrained corank 1 approximation (in the square case) with zero pattern S={(1,1)}. The number of solutions appear in Table 2.


## Maple code: Degree-matrix-rank2-diag-constraints
A maple code to estimate the degree (and the dimension, but it is 0)
of the critical points of the approximation of a data matrix
by a $m \times n$ matrix of rank 2, when a number of diagonal entries equal to zero.
We based the code on the dual formulation of the problem,
that is we express the rank constraints by parametrizing the kernel.
