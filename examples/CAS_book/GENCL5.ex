#!/bin/bash
#
# Generate configurations to allow for core-polarization in Ge III.
#
# This example is given in Table 6.6 of CAS [1].
#
# References:
#   1. Computational Atomic Structure: An MCHF Approach,
#      C.F. Fischer, T. Brage, and P. Jonsson, IoP Publishing
#      (2000).

BIN=../../bin

rm -f cfg.inp
$BIN/GENCL <<ST1

Ge III
 1s 2s 2p 3s 3p
3d(10)4s(2)

3d,4s,4p
1

1S

ST1

# Display the couplings found by GENCL 
cat cfg.inp



