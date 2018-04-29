#!/bin/bash
#
# Example for generating configurations with single orbital replacements.
#
# This example is given in Table 7.6 of CAS [1].
#
# References:
#   1. Computational Atomic Structure: An MCHF Approach,
#      C.F. Fischer, T. Brage, and P. Jonsson, IoP Publishing
#      (2000).

BIN=../../bin

rm -f cfg.inp
$BIN/GENCL <<ST1

Mg I
 1s 2s
2p(6)3s(1)3p(1)


2p=3p
2p=4p

3P

ST1

# Display the couplings found by GENCL 
cat cfg.inp



