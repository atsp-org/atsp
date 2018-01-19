#
#  Script file for hyperfine structure calculations of
#  1s(2)2p(1) 2P excited state in Li. See ref. 1.
#
#  References:
#
#  1. P. Jonsson, C.G. Wahlstrom, and C.F. Fischer,
#     Computer Physics Comm. 74, 399--414 (1993).
#

BIN=../../bin

# ----------------------------------------
# Step 1. Generate 1s, 2p wavefunctions using HF
# ---------------------------------------- 
rm -f wfn.inp
$BIN/HF <<S0
Li,2P,3.

1s(2)2p1(1)
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
#          (single configuration) 
# -----------------------------------------

rm -f cfg.inp Li1.out
cat >cfg.inp <<S1
Li  2P

  1s( 2) 2p1( 1)
     1S0     2P1     2P0
S1
cat cfg.inp

# -----------------------------------------
#  Step 3. Obtain Energy Expression
# -----------------------------------------
$BIN/NONH >Li1.out <<S2
n
y
S2

# -----------------------------------------
#  Step 4. Obtain MCHF Wavefunctions
# -----------------------------------------
echo "*" >> cfg.inp

$BIN/MCHF >>Li1.out <<S3
Li,2P,3.
all
y
y
y
S3
#
mv -f cfg.out li1.c
mv -f wfn.out li1.w

# -----------------------------------------
#  Step 5. Compute hyperfine parameters
#          for single configuration.
# -----------------------------------------
$BIN/HFS Li1.out <<S4
li1
0
y
m
y
0.0
3,3.256,-0.040
S4
#
#  Display results
#
cat li1.h >>Li1.out

# ------------------------------------------
#   Step 6: 
#   Case 2. Create configuration list (3 config)
# ------------------------------------------
rm -f cfg.inp  Li2.out
cat >cfg.inp <<S1
Li  2P

  1s( 2) 2p1( 1)
     1S0     2P1     2P0
  1s( 1) 2s2( 1) 2p2( 1)
     2S1     2S1     2P1     1S0     2P0
  1s( 1) 2s3( 1) 2p3( 1)
     2S1     2S1     2P1     3S0     2P0
S1
cat cfg.inp

# -------------------------------------------
#  Step 7: Get new energy expression
# -------------------------------------------
$BIN/NONH > Li2.out <<S2
n
y
S2

# --------------------------------------------  
#  Step 8: Get MCHF wave function (3 config)
# --------------------------------------------
echo "*" >> cfg.inp

cp li1.w wfn.inp
$BIN/MCHF >>Li2.out <<S3
Li,2P,3.
2p1,2s2,2p2,2s3,2p3
y
y
y
S3

#  Save results
mv -f cfg.out li2.c
mv -f wfn.out li2.w

# ---------------------------------------------
#  Step 9: Perform hyperfine calculation (3 config)
# --------------------------------------------- 
$BIN/HFS >>Li2.out <<S4
li2
1
y
m
y
1.0
3,3.256,-0.040
S4

#  Display results
cat li2.h >>Li2.out
