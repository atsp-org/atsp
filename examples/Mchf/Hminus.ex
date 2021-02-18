#! /usr/bin/sh
#
# MCHF calculation for the H- 1s(2) 1S ground state
# with a wave function expansion over the set of 
# electron configurations:
#
#     { 1s(2), 2s(2), 2p(2), 3s(2), 3p(2), 3d(2) }
#
# References:
#  1. A. R. P. Rau, "The Negative Ion of Hydrogen",
#     J. Astrophys. Astr. v. 17, 113--145 (1996).
#     Link:
#  http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?bibcode=1996JApA...17..113R
#
# K. Myneni, 2021-02-14

set -e
BIN=../../bin

#----------------------------------------
# 1) Use HF to obtain radial function for He 1s
#----------------------------------------

# In case a previous wfn.inp exists, remove it
# or it will be used by HF.
rm -f wfn.inp

$BIN/HF <<EOF
He,1S,2.

1s(2)
all
y
y
n
n
EOF

mv wfn.out wfn.inp

#----------------------------------------
# 2) Use He 1s wavefunction as input to solve
#    radial function for H- 1s, using HF
#----------------------------------------

$BIN/HF <<EOF
H-,1S,1.

1s(2)
all
y
y
n
n
EOF

mv wfn.out wfn.inp

#----------------------------------------
# 3) Specify the set of configuration states
#----------------------------------------
$BIN/GENCL <<EOF

H-

1s(2)
2s(2)
2p(2)
3s(2)
3p(2)
3d(2)



1S

EOF

# The following is needed to prevent MCHF
# from reading past the end of cfg.inp
echo "*" >> cfg.inp

#----------------------------------------
# 4) Compute energy expression
#----------------------------------------
$BIN/NONH <<EOF
n
y
EOF


#----------------------------------------
# 5) Run MCHF
#----------------------------------------
$BIN/MCHF <<EOF
H-,1S,1.
all
y
y
y
EOF

