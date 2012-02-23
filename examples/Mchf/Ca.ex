#
#   Script (command) file for  Mchf example  1
#
#..... Step 0.  Obtain initial estimates for the Ca- calculation from a
#            Hartree-Fock calculation for Ca 4s3p 3P
#
rm -f Ca.out
time HF >Ca.out <<S0
Ca,3P,20.
  1s  2s  2p  3s  3p
4s(1)4p(1)
all
y
y
n
n
S0
echo ' '
echo Move wfn.out wfn.inp
mv -f wfn.out wfn.inp
#
#..... Step 1.  Obtain the configuration state expansion for Ca-
#
time Gencl >>Ca.out <<S1

 Ca- 2P
 1s  2s  2p  3s  3p
4s(2)4p(1)

4s,3d,4p
0

2P

S1
echo '  '
echo Display the cfg.inp file produced
cat cfg.inp >>Ca.out
#
#..... Step 2.  Obtain the energy expression for the non-relativistic 
#            hamiltonian for this atomic state expansion
#
time Nonh >>Ca.out  <<S2
n
y
S2
echo  ' '
echo Display the int.lst file produced
cat int.lst
#
#..... Step 3.  Determine the radial functions for the Mchf expansion
#
time Mchf >>Ca.out  <<S3
Ca-,2P,20.
all
y
y
y
n
S3
echo  ' '
echo Display the cfg.out file produced
cat cfg.out  >>Ca.out
echo  ' '
# Save results
cp cfg.out Ca2P1.c
cp wfn.out Ca2P1.w
rm wfn.inp
echo Move wfn.out to wfn.inp
mv -f wfn.out wfn.inp
#
#..... Step 4.  Obtain a larger configuration state expansion for Ca-
#
time Gencl >> Ca.out <<S4

 Ca- 2P
 1s  2s  2p  3s  3p
4s(2)4p(1)

4s,5s,3d,4d,4p,5p
0

2P

S4
echo '  '
#
#   Apply Brillouin's Theorem to the 4s3d(3D)4p configuration state
#
#   ... create a file containing the delete commands
rm delete cfg.out
cat  > delete <<D1
21,22d
25,26d
37,38d
D1
#   ... apply the delete commands to cfg.inp, creating cfg.out
sed -f delete cfg.inp >cfg.out
echo ' '
echo Apply Brillouins Theorem by deleting the following from cfg.inp
#   ... display the deleted lines
diff cfg.inp cfg.out
#   ... delete cfg.inp and replace it with cfg.out
rm cfg.inp
mv -f cfg.out cfg.inp
echo Display the cfg.inp file produced
cat cfg.inp >>Ca.out
#
#..... Step 5.  Obtain the energy expression for the non-relativistic 
#            hamiltonian
#
time Nonh >> Ca.out  <<S5
n
y
S5
echo  ' '
echo Display the int.lst file produced
cat int.lst >>Ca.out
#
#..... Step 6.  Determine the radial functions for the Mchf expansion
#
time Mchf  >>Ca.out  <<S6
Ca-,2P,20.
4s,4p,5s,3d,4d,5p
y
y
y
n
S6
echo  ' '
echo Display the cfg.out file produced
cat cfg.out >>Ca.out
mv -f cfg.out Ca2P2.c
mv -f wfn.out Ca2P2.w
