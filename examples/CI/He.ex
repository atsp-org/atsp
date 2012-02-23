#
#   Script (command) file illustrating the use of the MCHF_CI program
#   for study of the mass-polarization correction.
#
#   Case 1.  Mass-Polarization correction for Helium 1s2p 1P.
#
#   Step 1.  Clear cfg.inp, wfn.inp, and create cfg.inp
#
rm -f cfg.inp wfn.inp he1P.c he1P.w  He.out
#
cat >cfg.inp <<S1
  He  1P

  1s( 1) 2p1( 1)
     2S1     2P1     1P0
 2p2( 1) 3d2( 1)
     2P1     2D1     1P0
S1
#
#   Step 2. Obtain the energy expression
#
time Nonh >He.out <<S2
n
y
S2
#   Step 3.  Obtain the radial functions
#
time Mchf >>He.out  <<S3
He,1P,2.
all
y
y
y
n
S3
#
#   Step 4. Save cfg.out and wfn.out as he1P.c and he1P.w
mv cfg.out he1P.c
mv wfn.out he1P.w
#
#   Step 4.  Perform CI mass-polarization calculation
#
time Ci >>He.out <<S4
he1P
n
y
g
n
4.022
1
S4
echo ' '
echo Display the file he1P.l that was produced
cat he1p.l
