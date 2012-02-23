#
#   Script (command) file illustrating the use of the LSTR and LSJTR
#   programs in the MCHF Atomic Structure Package.
#
#   Case 1.  A non-relativistic transition probability calculation for
#            3s(2)3p -- 3s(2)3d in Al restricted to the complex
#  
#   Step 1.  Generate the configuration state list for the Al 2P
#            state
#
rm Al.out tr.lsj
#
time Gencl >Al.out <<S1

 AL 2P
 1s  2s  2p
3s(2)3p(1)

3s,3p,3d
0

2P

S1
echo ' '
echo Move the cfg.inp file produced to al2P.c
mv -f cfg.inp al2P.c
#
#   Step 2.  Determine the configuration state list for the Al 2D atomic
#            state
#
time Gencl >>Al.out <<S2

 AL 2D
 1s  2s  2p
3s(2)3d(1)

3s,3p,3d
0

2D

S2
echo ' '
echo Move the cfg.inp file produced to al2D.c
mv -f cfg.inp al2D.c
#
#
#   Step 3. Determine the expressions for the dipole operator for E1 
#            transitions 2P Al -- 2D Al
#
time Mltpol >>Al.out  <<S3
al2P
al2D
n

e1
*
S3
echo '  '
echo Display mltpol.lst file produced
cat mltpol.lst >>Al.out
