#
#  Example of hyperfine structure calculations for
#  2P excited state in Na.
#
#  See Tables 8.5--8.7 in Ref. 1.
#
#  References:
#
#  1. C.F. Fischer, T. Brage, and P. Jonsson, Computational
#     Atomic Structure,
#
#  2. M.S. Safronova, W.R. Johnson, and A. Derevianko,
#     Phys. Rev. A 60, 4476 (1999).
#
#  3. D. Das and V. Natarajan, J. Phys. B 40, 035001 (2008).
#

BIN=../../bin

# ----------------------------------------
# Step 1. Generate 1s, 2s, 2p, and 3p wavefunctions using HF
# ---------------------------------------- 
rm -f wfn.inp
$BIN/HF <<S0
Na,2P,11.
  1s  2s
2p(6)3p(1)
all
y
y
n
n
S0

# Use output wavefunctions from HF as inputs to MCHF
mv -f wfn.out wfn.inp

# -----------------------------------------
#  Step 2. Generate configuration list
#          (14 configurations) 
# -----------------------------------------

rm -f cfg.inp
cat >cfg.inp <<S1
 Na 2P                                                       
                                                                         
  1s( 2)  2s( 2)  2p( 6)  3p( 1)
     1S0     1S0     1S0     2P1     1S0     1S0     2P0
  1s( 1) 1s1( 1)  2s( 2)  2p( 6)  3p( 1)
     2S1     2S1     1S0     1S0     2P1     3S0     3S0     3S0     2P0
  1s( 1) 3d1( 1)  2s( 2)  2p( 6)  3p( 1)
     2S1     2D1     1S0     1S0     2P1     1D0     1D0     1D0     2P0
  1s( 1) 3d2( 1)  2s( 2)  2p( 6)  3p( 1)
     2S1     2D1     1S0     1S0     2P1     3D0     3D0     3D0     2P0
  1s( 2)  2s( 1) 2s2( 1)  2p( 6)  3p( 1)
     1S0     2S1     2S1     1S0     2P1     2S0     3S0     3S0     2P0
  1s( 2)  2s( 1) 3d3( 1)  2p( 6)  3p( 1)
     1S0     2S1     2D1     1S0     2P1     2S0     1D0     1D0     2P0
  1s( 2)  2s( 1) 3d4( 1)  2p( 6)  3p( 1)
     1S0     2S1     2D1     1S0     2P1     2S0     3D0     3D0     2P0
  1s( 2)  2s( 2)  2p( 5) 2p1( 1)  3p( 1)
     1S0     1S0     2P1     2P1     2P1     1S0     2P0     1P0     2P0
  1s( 2)  2s( 2)  2p( 5) 2p2( 1)  3p( 1)
     1S0     1S0     2P1     2P1     2P1     1S0     2P0     1D0     2P0
  1s( 2)  2s( 2)  2p( 5) 2p3( 1)  3p( 1)
     1S0     1S0     2P1     2P1     2P1     1S0     2P0     3S0     2P0
  1s( 2)  2s( 2)  2p( 5) 2p4( 1)  3p( 1)
     1S0     1S0     2P1     2P1     2P1     1S0     2P0     3P0     2P0
  1s( 2)  2s( 2)  2p( 5) 2p5( 1)  3p( 1)
     1S0     1S0     2P1     2P1     2P1     1S0     2P0     3D0     2P0
  1s( 2)  2s( 2)  2p( 5) 4f1( 1)  3p( 1)
     1S0     1S0     2P1     2F1     2P1     1S0     2P0     1D0     2P0
  1s( 2)  2s( 2)  2p( 5) 4f2( 1)  3p( 1)
     1S0     1S0     2P1     2F1     2P1     1S0     2P0     3D0     2P0
*
S1

# -----------------------------------------
#  Step 3. Obtain Energy Expression
# -----------------------------------------
$BIN/NONH <<S2
n
y
S2

# -----------------------------------------
#  Step 4. Obtain MCHF Wavefunctions
# -----------------------------------------
$BIN/MCHF <<S3
Na,2P,11.
1s1,3d1,3d2,2s2,3d3,3d4,2p1,2p2,2p3,2p4,2p5,4f1,4f2
y
y
y
S3
#
mv -f cfg.out Na_2P.c
mv -f wfn.out Na_2P.w

# -----------------------------------------
#  Step 5. Compute hyperfine constants A and B
#          for 2P_1/2 and 2P_3/2 states.
#          Nuclear data from Ref. 3
# -----------------------------------------
$BIN/HFS <<S4
Na_2P
0
y
m
y
0.0
3,2.2175219,0.1006
S4
#
#  Display results 
#
#  Compare with Table 8.7 in Ref. 1; Note that Ref. 1 
#  uses a dummy argument of 1 barns for the electric
#  quadrupole moment for the Na nucleus. In the calc
#  above, we used the known value for the Na (atomic
#  wt=23) nucleus. Our value for B is different from
#  Ref. 1, and is within 10% of the experimental value.
cat Na_2P.h

