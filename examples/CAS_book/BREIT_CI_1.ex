#!/bin/bash
#
# Compute the fine structure of the He 1s2p 3P term, including
# only the spin-orbit interaction.
#
# See section 7.6 and Table 7.2 of CAS
#
BIN=../../bin
# ---------------------------------------
#   Step 1: Compute the non-relativistic HF wavefunction for the
#           He 1s2p 3P state.
# ---------------------------------------
$BIN/HF <<S0
He,3P,2.

1s(1)2p(1)
all
y
y
n
n
S0

cp -f wfn.out He_3P.w

# -------------------------------------
#   Step 2: Generate the configuration file for this state
# -------------------------------------
$BIN/GENCL <<S1

He 3P

1s(1)2p(1)



3P

S1

cp -f cfg.inp He_3P.c

# --------------------------------------
#   Step 3: Compute matrix elements of the Breit-Pauli
#           Hamiltonian
# --------------------------------------
$BIN/BREIT <<S2
2
y
n
y
n
n
y
S2

# ---------------------------------------
#   Step 4: Construct the full Hamiltonian and find
#           the eigenvalues
# ---------------------------------------
$BIN/CI <<S3
He_3P
y
n
1
4,0
n
S3

cat He_3P.j

# Note there are three fine structure levels for the 3P
# term: J = 2, 1, 0. The printed output displays 2J.

 
