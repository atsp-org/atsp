#! /bin/sh
#
# MCHF calculation for the He 1S ground state
# with a wave function expansion over the set of
# following electron configurations:
#
#   { 1s(2), 2s(2), 2p(2), 3s(2), 3p(2), 3d(2), 
#     4s(2), 4p(2), 4d(2), 5s(2), 5p(2), 
#     6s(2), 6p(2) }
#
set -e
BIN=../../bin

#----------------------------------------
# 1) Use HF to obtain radial function for 1s
#----------------------------------------

# In case a previous wfn.inp exists, remove it
# or it will be used by HF.
rm -f wfn.inp

$BIN/HF <<EOF
He,1S,2

1s(2)
all
y
y
n
n
EOF

mv wfn.out wfn.inp

#----------------------------------------
# 2) Specify starting set of configuration states
#----------------------------------------
$BIN/GENCL <<EOF

He

1s(2)
2s(2)
3s(2)
4s(2)



1S

EOF

# The following is needed to prevent MCHF
# from reading past the end of cfg.inp
echo "*" >> cfg.inp

#----------------------------------------
# 3) Compute energy expression
#----------------------------------------
$BIN/NONH <<EOF
n
y
EOF


#----------------------------------------
# 4) Run MCHF
#----------------------------------------
$BIN/MCHF <<EOF
He,1S,2
all
y
y
y
n

EOF

mv wfn.out wfn.inp

#
# 5) Run Gencl to add more electron configurations
#
$BIN/GENCL <<EOF

He

1s(2)
2s(2)
3s(2)
4s(2)
5s(2)
2p(2)
3p(2)
4p(2)




1S

EOF
echo "*" >> cfg.inp

#
# 6) Run Nonh for energy expression
#
$BIN/NONH <<EOF
n
y
EOF

#
# 7) Run MCHF
#
$BIN/MCHF <<EOF
He,1S,2
all
y
y
y
n

EOF

mv wfn.out wfn.inp

#
# 8) Add more configurations
#
$BIN/GENCL <<EOF

He

1s(2)
2s(2)
3s(2)
4s(2)
5s(2)
2p(2)
3p(2)
4p(2)
5p(2)



1S

EOF
echo "*" >> cfg.inp

#
# 9) Run Nonh
#
$BIN/NONH <<EOF
n
y
EOF

#
# 10) Run Mchf again
#
$BIN/MCHF <<EOF
He,1S,2
all
y
y
y
n

EOF

mv wfn.out wfn.inp

#
# 11) Add remaining configurations
#
$BIN/GENCL <<EOF

He

1s(2)
2s(2)
3s(2)
4s(2)
5s(2)
6s(2)
2p(2)
3p(2)
4p(2)
5p(2)
6p(2)
3d(2)
4d(2)



1S

EOF
echo "*" >> cfg.inp

#
# 12) Obtain energy expression
#
$BIN/NONH <<EOF
n
y
EOF

#
# 13) Run MCHF final time for all configs
#
$BIN/MCHF <<EOF
He,1S,2
all
y
y
y
n

EOF

