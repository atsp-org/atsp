#!/bin/bash
#
# Isotope shift calculation of the 1s(1)2p(1) 3P state
# between 7Li II and 6Li II, with a wave function
# expansion over the set of configuration states:
#
#    { 1s(1)2p1(1), 2s(1)3p1(1), 2p2(1)3d(1) }
#
# See section 8.5.1 and Table 8.1 in ATSP book.
#
BIN=../../bin

# ----------------------------------------
# Step 1. Generate 1s, 2p wavefunctions using HF
# ---------------------------------------- 
rm -f wfn.inp
$BIN/HF <<ST0
LiII,3P,3.

1s(1)2p(1)
all
y
y
n
n
ST0

# Use output wavefunctions from HF as inputs to MCHF
mv -f wfn.out wfn.inp

# -----------------------------------------
#  Step 2. Generate configuration list
#          (single configuration) 
# -----------------------------------------

rm -f cfg.inp
cat >cfg.inp <<ST1
Li II 3P

  1s( 1) 2p1( 1)
     2S1     2P1     3P0
  2s( 1) 3p1( 1)
     2S1     2P1     3P0
 2p2( 1)  3d( 1)
     2P1     2D1     3P0
*    
ST1
cat cfg.inp

# -----------------------------------------
#  Step 3. Obtain Energy Expression
# -----------------------------------------
$BIN/NONH <<ST2
n
y
ST2

# -----------------------------------------
#  Step 4. Obtain MCHF Wavefunctions
# -----------------------------------------

$BIN/MCHF <<ST3
LiII,3P,3.
all
y
y
y
ST3
#
mv -f cfg.out LiII_3P.c
mv -f wfn.out LiII_3P.w

# -----------------------------------------
#  Step 5. Compute isotope shift corrections
# -----------------------------------------
$BIN/ISO <<ST4
LiII_3P
1
n
1
y
y
7.0160030,6.0151214
ST4

