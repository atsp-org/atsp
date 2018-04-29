#!/bin/bash
#
# GENCL example for non-orthogonal orbitals.
#
# This example is given in Table 4.3 of CAS [1].
#
# References:
#   1. Computational Atomic Structure: An MCHF Approach,
#      C.F. Fischer, T. Brage, and P. Jonsson, IoP Publishing
#      (2000).

BIN=../../bin

rm -f cfg.inp
$BIN/GENCL <<ST1

Boron
 1s
2s(2)2p1(1)
2p2(3)



2P

ST1

# Display the couplings found by GENCL 
cat cfg.inp



