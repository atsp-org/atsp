#! /bin/sh

set -e
#
#   MCHF Script file showing the importance of the method
#            used in solving the Mchf SCF equations for a core
#            excited state in Sodium, namely 2p(5)3s3p 2S
#
#  Step 1.  Obtain the Hartree-Fock orbitals for the core from 
#           a calculation for 2p(5)3s 3P.
#
#
rm -f Na.out
time ../../src/HF >Na.out <<S1
Na+,3P,11.
  1s  2s
2p(5)3s1(1)
all
y
y
n
n
S1
#
echo '  '
echo Move wfn.out wfn.inp
mv -f wfn.out wfn.inp
#
#  Step 2. Perform a fixed core calculation for the 2p(5)3s1(3P)3p1 2S 
#          state, first creating the cfg.inp file, then using Nonh, 
#          followed by Mchf.
#
rm -f cfg.inp
cat  >cfg.inp <<S2a
Na 2S
  1s  2s
  2p( 5) 3s1( 1) 3p1( 1)
     2P1     2S1     2P1     3P0     2S0
  2p( 5) 3p2( 1) 3d2( 1)
     2P1     2P1     2D1     3D0     2S0
*
S2a
time ../../src/NONH >>Na.out <<S2b
n
y
S2b
time ../../src/MCHF >>Na.out <<S2c
Na,2S,11.
=4
y
y
y
n
S2c
# Save results
mv -f cfg.out Na2S.c
mv -f wfn.out Na2S.w
