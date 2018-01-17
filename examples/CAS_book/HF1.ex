#! /bin/sh

set -e
rm -f wfn.inp

../../bin/HF <<EOF
He,3P,2.

1s(1)2p(1)
all
y
y
n
n
EOF
