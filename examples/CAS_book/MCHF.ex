#! /bin/sh

set -e

#----------------------------------------

../../src/GENCL <<EOF

Al
 1s  2s  2p
3s(2)3p(1)

3s,3p,3d
0

2P

EOF

#----------------------------------------

../../src/NONH <<EOF
n
y
EOF

#----------------------------------------

../../src/HF <<EOF
Al,2P,13.
  1s  2s  2p  3s
3p(1)
all
y
y
n
n
EOF

mv wfn.out wfn.inp

#----------------------------------------

# This is needed, because the code in MCHF first read all records from
# `cfg.inp` up until "*" and then later tries to read any additional records as
# "DETERMINE ADDITIONAL ORTHOGONALITY PAIRS". Without adding the explicit "*"
# at the end, the second read fails with "reading past EOF", even though the
# EOF specifier is there --- it might be that some Fortran compilers didn't
# fail. In any case, this fixes it.
echo "*" >> cfg.inp

../../src/MCHF <<EOF
Al,2P,13.
=3
y
y
y
EOF
