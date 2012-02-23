#
#  Script file for hyperfine program
# 
#  Generate configuration list
#
rm -f cfg.inp Li1.out
cat >cfg.inp <<S1
Li  2P

  1s( 2) 2p1( 1)
     1S0     2P1     2P0
S1
cat cfg.inp
#
#  Obtain Energy Expression
#
Nonh >Li1.out <<S2
n
y
S2
#
#  Obtain Wavefunction
#
Mchf >>Li1.out <<S3
Li,2P,3.
all
y
y
y
n
S3
#
mv -f cfg.out li1.c
mv -f wfn.out li1.w
#
#  Compute Hyperfine parameters
#
time Hfs Li1.out <<S4
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
#
#                         Case 2.
#
#   Create configuration list
#
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
#
#  Get energy expression
#
Nonh > Li2.out <<S2
n
y
S2
#  
#  Get MCHF wave function
#
cp li1.w wfn.inp
Mchf >>Li2.out <<S3
Li,2P,3.
2p1,2s2,2p2,2s3,2p3
y
y
y
n
S3
#  Save results
mv -f cfg.out li2.c
mv -f wfn.out li2.w
#
#  Perform Hyperfine calculation
time Hfs >>Li2.out <<S4
li2
1
y
m
y
1.0
3,3.256,-0.04
S4
#
#  Display results
cat li2.h >>Li2.out
