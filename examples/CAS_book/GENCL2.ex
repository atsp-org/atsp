#!/bin/bash
#
# Generate an SD expansion for the 1s(2)2p(1) 2P term.
#
# This example is given in Table 4.2 of CAS [1].
#
# References:
#   1. Computational Atomic Structure: An MCHF Approach,
#      C.F. Fischer, T. Brage, and P. Jonsson, IoP Publishing
#      (2000).

BIN=../../bin

rm -f cfg.inp
$BIN/GENCL <<ST1

Lithium

1s(2)2p(1)


sd
2s,2p,3s,3p,3d
1
2
2P

ST1

# Display the couplings found by GENCL 
cat cfg.inp



