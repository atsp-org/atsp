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

echo "" >> cfg.inp
echo "" >> cfg.inp

../../src/MCHF <<EOF
Al,2P,13.
=3
y
y
y
EOF
