#
#   Script (command) file illustrating the use of the LSTR and LSJTR
#   programs in the MCHF Atomic Structure Package.
#
#   Case 1.  A non-relativistic transition probability calculation for
#            3s(2)3p -- 3s(2)3d in Al, using the MCHF approximation.
#  
#...To reproduce these results, there should be no wfn.inp file!
#
rm -f Al.out tr.lsj
#
#   Step 0.  Obtain intial estimates for the core radial functions
#
time HF >Al.out  <<S0
Al,1S,13.
  1s  2s  2p  3s

all
y
y
n
n
S0
# echo ' '
# echo Move wfn.out wfn.inp
mv -f wfn.out wfn.inp
#   
#   Step 1.  Generate the configuration state list for the Al 2P
#            state
#
time Gencl >>Al.out <<S1

 AL 2P
 1s  2s  2p
3s(2)3p(1)

3s,3p,3d
0

2P

S1
# echo '  '
# echo Display the cfg.inp file produced
cat cfg.inp 
#
#   Step 2.  Determine the energy expression for the non-relativistic
#            hamiltonian for this atomic state
#
time Nonh >>Al.out  <<S2
n
y
S2
# echo  ' '
# echo Display the int.lst file produced
cat int.lst   >>Al.out
#
#   Step 3.  Determine the radial function for the MCHF approximation
#
time Mchf >>Al.out  <<S3
Al,2P,13.
all
y
y
y
n
S3
# echo  ' '
# echo Move cfg.out to al2P.c
mv -f cfg.out al2P.c
# echo  ' '
# echo Move wfn.out to al2P.w
mv -f wfn.out al2P.w
#
#   Step 4.  Derive the energy expression for the Breit-Pauli 
#            interaction matrix
#
time Breit >>Al.out <<S4
2
n
y
y
S4
# echo  ' '
#   
#   Step 5.  Determine eigenvalues and eigenvectors of the Breit-Pauli 
#            interaction matrix for a range of J values
#
time Ci >>Al.out  <<S5
al2P
y
n
1
3,1
n
S5
# echo '  '
# echo Display the al2P.l file produced
cat al2P.l 
# echo '  '
# echo Display the al2P.j file produced
cat al2P.j
#
#   Step 6.  Determine the configuration state list for the Al 2D atomic
#            state
#
time Gencl >>Al.out <<S6

 AL 2D
 1s  2s  2p
3s(2)3d(1)

3s,3p,3d
0

2D

S6
#
#   Step 7.   Determine the non-relativistic interaction matrix for the 
#             Al 2D atomic state
#
time Nonh >>Al.out  <<S7
n
y
S7
# echo '  '
# echo Display the int.lst file produced
cat int.lst >>Al.out
#   Step 8.   Determine radial functions for the Al 2D state with the 
#             same core are the Al 2P state
#

time Mchf >>Al.out<<S8
Al,2D,13.
=3
y
y
y
n
S8
#              Recycle to improve convergence
mv -f wfn.out wfn.inp
mv -f cfg.out cfg.inp
time Mchf >>Al.out  <<S8b
Al,2D,13.
=3
y
y
y
n
S8b
# echo ' '
# echo Display the cfg.out file produced
cat cfg.out 
rm al2D.c al2D.w
# echo ' '
# echo Move the cfg.out file to al2D.c
mv -f cfg.out al2D.c
# echo ' '
# echo Move the wfn.out file to al2D.w
mv -f wfn.out al2D.w
#
#   Step 9.  Determine the expressions for the Breit-Pauli interaction 
#             matrix for the Al 2D atomic state
#
time Breit >>Al.out <<S9
2
n
y
y
S9
# echo ' '
# echo Display the int.lst produced
cat int.lst >>Al.out
#
#   Step 10.  Determine eigenvalues and eigenvectors for the Al 2D 
#            states for a range of J values
#
time Ci >>Al.out <<S10
al2D
y
n
1
5,3
n
S10
# echo '  '
# echo Display the al2D.l file produced
cat al2D.l 
# echo '  '
# echo Display the al2D.j file produced
cat al2D.j 
#
#   Step 11. Determine the expressions for the dipole operator for E1 
#            transitions 2P Al -- 2D Al
#
time Mltpol >>Al.out <<S11
al2P
al2D
n

e1
*
S11
# echo '  '
# echo Display the mltpol.lst file produced
cat mltpol.lst >> Al.out
#  
#   Step 12.  Determine the non-relativistic transition data
#
time Lstr >>Al.out  <<S12
al2P
al2D
y
0.001
y
S12
#
#   Step 13.  Determine the Breit-Pauli LSJ transition data
#
time Lsjtr >>Al.out <<S13
al2P
al2D
y
0.001
y
S13
# echo '  '
# echo Display the tr.lsj file produced
cat tr.lsj >>Al.out
#
#   Step 14.  Display a list of lines from the transition data in tr.lsj
#
Lines >>Al.out <<S14
0.0
2
S14
# echo '  '
# echo Concatenate the al2P.j and al2D.j files to form file al.j
cat al2P.j al2D.j >al.j
#
#   Step 15.  Display a list of levels from the combined list of 
#             eigenvalues and eigenvectors of the two states.
#
#
Levels >>Al.out <<S15
al.j
S15
