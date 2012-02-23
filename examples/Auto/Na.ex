#
#  Case 2.  Script (command) file illustrating the use of the MCHF 
#  Atomic Structure Package for the study of core excited quartets
#  2p(5)3s3p in Sodium.
#
#  Step 1.  Obtain the Hartree-Fock orbitals for the core from 
#           a calculation for 2p(5)3s 3P.
#
# Remove both cfg.inp and wfn.inp 
rm -f cfg.inp wfn.inp
# Remove Naquart.c Naquart.w
rm -f Naquart.c Naquart.w delete
#
time HF <<S1
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
#          state, first creating the cfg.inp file, then using NONH, 
#          followed by MCHF.
#
cat  >cfg.inp <<S2a
Na 2S
  1s  2s
  2p( 5) 3s1( 1) 3p1( 1)
     2P1     2S1     2P1     3P0     2S0
  2p( 5) 3p2( 1) 3d2( 1)
     2P1     2P1     2D1     3D0     2S0
S2a
rm -f Na.out
time Nonh >Na.out <<S2b
n
y
S2b
time Mchf >>Na.out <<S2c
Na,2S,11.
=4
y
y
y
n
S2c
echo ' '
echo Move wfn.out to wfn.inp and remove cfg.inp
rm -f cfg.inp wfn.inp
mv -f wfn.out wfn.inp
#
#  Step 3. Perform a fixed core calculation for the 2p(5)3s1(3P)3p1 4S
#          state by generating the configuration state list, then
#          using NONH followed by MCHF
#
time Gencl >>Na.out <<S3a

Na core-excited 4S
  1s  2s
2p(5)3s1(1)4s1(0)3p1(1)
2p(5)3p2(1)3d2(1)
2p(5)3s1(1)4p1(1)
2p(5)4s1(1)3p1(1)
2p(5)4s1(1)4p1(1)
2p(5)3p2(1)4d2(1)



4S

S3a
#
time Nonh >>Na.out <<S3b
n
y
S3b
#
time Mchf >>Na.out <<S3c
Na,4S,11.
=3
y
y
y
n
S3c
echo '  '
echo Move the wfn.out file produced to Naquart.w
mv -f wfn.out Naquart.w
echo Display the cfg.out file produced
cat cfg.out
#
#  Step 4. Produce the LSJ configuration state list for the
#          terms of 2p(5)3s3p quartet states
#
time Gencl >>Na.out  <<S4

Na quartet
  1s  2s
2p(5)3s1(1)4s1(0)3p1(1)
2p(5)3p2(1)3d2(1)
2p(5)3s1(1)4p1(1)
2p(5)4s1(1)3p1(1)
2p(5)4s1(1)4p1(1)
2p(5)3p2(1)4d2(1)




S4
#
#  Step 5. Determine the non-relativistic interactions
#          allowing for non-orthogonal orbitals
#
time Nonh >>Na.out <<S5
n
y
S5
echo ' '
echo Move int.lst file to int1.lst 
mv -f int.lst int1.lst
#
#  Step 6. Determine the relativistic corrections in the
#          restricted two-body mode
#
time Breit >>Na.out <<S6
1
n
y
n
0,54
y
y
y
S6
echo ' '
echo Append the relativistic int.lst to int1.lst and delete int.lst
cat int.lst >>int1.lst
rm -f int.lst
echo Delete lines  3264 to 3283 from int1.lst and copy to int.lst
#  
# ... create a file containing the delete commands
cat > delete <<D1
3264,3283d
D1
# ... apply the delete commands
sed -f delete int1.lst >int.lst
echo Move cfg.inp to Naquart.c
mv -f cfg.inp Naquart.c
#
#  Step 6.  Perform the CI calculations for a range of J-values
#
time Ci >> Na.out <<S6
Naquart
y
n
9
7,1
n
S6
echo ' '
echo Display the beginning of the Naquart.l file produced
head Naquart.l
echo Display the end of the Naquart.l file produced
tail Naquart.l
echo Display the beginning of the Naquart.j file produced
head Naquart.j
echo Display the end of the Naquart.j file produced
tail Naquart.j
#
#  Step 7.  Produce an energy level table
#
time Levels >>Na.out <<S7
Naquart.j
S7
#
#  Step 8.  Show the major components of the eigenvectors
#
time Comp >>Na.out <<S8
Naquart
0.1
3
S8
#
# Step 9. Add the energy of the target and the continuum state to
# the .c file.
#
rm -f ecore
cat > ecore << S9a
   Na    2S      -161.6769459
S9a
cat >cont << S9b
  2p( 6) kd1( 1)                         1.0
     1S0     2D1     2D0
S9b
rm -f temp.c
cat ecore Naquart.c cont >temp.c
rm -f delete
cat >delete << S9c
2d
S9c
# Apply deleting commands.
sed -f delete temp.c >cfg.inp
#
#  Display the cfg.inp file
#
cat cfg.inp
#
#  Step 10. Run nonh for this case.
#	    Observe that we only need to calculate the last row.
#
time Nonh >> Na.out << S10
n
n
1,0
S10
rm -f Naquart.c
mv -f cfg.inp Naquart.c
#
#  Step 11. Run auto for this case.
#
time Auto >>Na.out << S11
Naquart
3
Na,2D,11.
y
S11
#
# Display the auto.dat file.
#
cat auto.dat >>Na.out
