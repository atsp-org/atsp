#!/bin/bash
#
# Generate an active set expansion for the 3s(2)3p(1) 2P state
# of aluminum.
#
# This example is given in Table 4.4 of CAS [1].
#
# References:
#   1. Computational Atomic Structure: An MCHF Approach,
#      C.F. Fischer, T. Brage, and P. Jonsson, IoP Publishing
#      (2000).

BIN=../../bin

rm -f cfg.inp
$BIN/GENCL <<ST1

Al
 1s 2s 2p
3s(2)3p(1)

3s,3p,3d
0

2P

ST1

# Display the couplings found by GENCL 
cat cfg.inp



