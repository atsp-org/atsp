#! /bin/sh

set -e
BIN=../../bin
#----------------------------------------

$BIN/GENCL <<EOF

C
 1s  2s
2p(2)




EOF

#----------------------------------------

$BIN/NONH <<EOF
y
y
EOF
