#!/bin/bash
#
# Generate the possible angular momentum couplings for a
# 3P term, with the electron configuration:
#
#    1s(2)2s(1)2p(2)3s(1)
#
# Coupling within the subshells leads to the four LS terms
# (see the program TERMS and CAS Table 2.1):
#
#   1S 2S 3P 2S
#
# Each of these may be coupled together in the following ways
# to produce a 3P term:
#
#   1) 1S + 2S --> 2S + 3P --> 2P + 2S --> 3P
#
#   2) 1S + 2S --> 2S + 3P --> 4P + 2S --> 3P
#
# This example is given in Table 2.2 of CAS [1].
#
# References:
#   1. Computational Atomic Structure: An MCHF Approach,
#      C.F. Fischer, T. Brage, and P. Jonsson, IoP Publishing
#      (2000).

BIN=../../bin

rm -f cfg.inp
$BIN/GENCL <<ST1

Example
1s
2s(1)2p(2)3s(1)



3P

ST1

# Display the couplings found by GENCL 
cat cfg.inp



