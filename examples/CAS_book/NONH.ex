#!/bin/bash
#
# Example of program NONH execution from CAS Table 2.3
#

set -e
BIN=../../bin
#----------------------------------------
rm -f cfg.inp
$BIN/GENCL <<EOF

Example

1s(1)2p(1)
2p(1)3d(1)



3P

EOF

cat cfg.inp

#----------------------------------------

$BIN/NONH <<EOF
y
y
EOF
