#! /bin/sh

set -e
rm -f wfn.inp

../../bin/HF <<EOF
Li,2S,3.
 1s
2s(1)
all
y
y
n
n
EOF
