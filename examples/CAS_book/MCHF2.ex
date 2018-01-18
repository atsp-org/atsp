#! /bin/sh
#
# Example from CAS, ch. 5, Table 5.3:
#
# MCHF calculation for the He 1s(2) 1S ground state
# with a wave function expansion over the set of 
# configuration states, {1s(2), 2s(2), 2p(2)} 1S.
#
set -e

#----------------------------------------
# 1) Use HF to obtain radial function for 1s
#----------------------------------------

# In case a previous wfn.inp exists, remove it
# or it will be used by HF.
rm -f wfn.inp

../../src/HF <<EOF
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
# 2) Specify the set of configuration states
#----------------------------------------
../../src/GENCL <<EOF

He

1s(2)
2s(2)
2p(2)



1S

EOF

# The following is needed to prevent MCHF
# from reading past the end of cfg.inp
echo "*" >> cfg.inp

#----------------------------------------
# 3) Compute energy expression
#----------------------------------------
../../src/NONH <<EOF
n
y
EOF


#----------------------------------------
# 4) Run MCHF
#----------------------------------------
../../src/MCHF <<EOF
He,1S,2.
all
y
y
y
EOF

