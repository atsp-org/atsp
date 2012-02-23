#
#  Case 2.  M1 nd E2 LSJTR transitions in Carbon-like Iron.
#
#   Script (command) file illustrating the use of the MCHF
#   atomic structure package for the study of the E2 and M1
#   transitions in Carbon.  This is the case where it may be
#   necessary to have two copies of some files.
#  
rm -f wfn.inp cfg.inp Lsjtr2.out
rm -f iron.c iron.w tr.lsj
#   
#   Step 1.  Generate the configuration state list for the 2p(2) 1S 
#            complex.
#
time Gencl >Lsjtr2.out <<S1

 Carbon-like Iron 
 1s  
2s(2)2p(2)

2s,2p
0

1S

S1
#
#   Step 2.  Determine the energy expression for the non-relativistic
#            hamiltonian
#
time Nonh >>Lsjtr2.out  <<S2
n
y
S2
echo  ' '
echo Display the int.lst file produced
cat int.lst >>Lsjtr2.out
#
#   Step 3.  Determine the radial function for the MCHF approximation
#
time Mchf >>Lsjtr2.out   <<S3
Fe+20,1S,26.
all
y
y
y
n
S3
echo  ' '
echo Display the cfg.out file produced
cat cfg.out  >>Lsjtr2.out
echo  ' '
echo Move wfn.out to iron.w
mv -f wfn.out iron.w
#
#   Step 4.  Generate the configuration state list for the 2p(2) LSJ 
#            terms that include single and double replacements
#
time Gencl >>Lsjtr2.out <<S4

 Carbon-like Iron
 1s  
2s(2)2p(2)
2p(4)





S4
echo '  '
echo Display the cfg.inp file produced
cat cfg.inp >>Lsjtr2.out
#
#   Step 5.  Derive the energy expression for the Breit-Pauli 
#            interaction matrix
#
time Breit >>Lsjtr2.out  <<S5
2
n
y
y
S5
echo  ' '
echo Display the int.lst file produced
cat int.lst >> Lsjtr2.out
echo ' '
echo Move cfg.inp to iron.c
mv -f cfg.inp iron.c
#   
#   Step 6.  Determine eigenvalues and eigenvectors of the Breit-Pauli 
#            interaction matrix for a range of J values
#
time Ci >>Lsjtr2.out <<S6
iron
y
n
6
4,0
n
S6
echo '  '
echo Display the iron.l file produced
cat iron.l
echo '  '
echo Display the iron.j file produced
cat iron.j
#
#   Step 7. Determine the expressions for the E2 and M1
#            transitions operators 
#
#   Make second copy of iron.c, iron.w, and iron.j
cp -f iron.c iron2.c
cp -f iron.w iron2.w
cp -f iron.j iron2.j
time Mltpoli >>Lsjtr2.out
iron
iron2
y
E2
M1
*
S7
echo '  '
echo Display the mltpol.lst file produced
cat mltpol.lst >>Lsjtr2.out
#
#   Step 8.  Determine the Breit-Pauli LSJ transition data
#
time Lsjtr >>Lsjtr2.out  <<S8
iron
iron2
y
0.0000001
y
S8
echo '  '
echo Display the tr.lsj file produced
cat tr.lsj  >>Lsjtr2.out
#
#   Step 9.  Display a list of lines from the transition data in tr.lsj
#
Lines >>Lsjtr2.out  <<S9
0.0
2
S9
#
#   Step 10.  Display a list of levels 
#
Levels >>Lsjtr2.out  <<S10
iron.j
S10
