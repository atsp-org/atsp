#
#  Case 1. Script (command) file illustrating the use of the MCHF
#  Atomic Structure Package for the study of autoionization of
#  1s2p(2) 2D in Litium
#
#  Step 1.  Generate cfg.inp for localized state.
#
rm -f cfg.inp
cat >cfg.inp <<S1
 Li                                                          
 
  1s( 1)  2p( 2)
     2S1     1D2     2D0
  1s( 1)  3p( 2)
     2S1     1D2     2D0
  1s( 1) 3d2( 2)
     2S1     1D2     2D0
  1s( 1)  2s( 1) 3d1( 1)
     2S1     2S1     2D1     1S0     2D0
  1s( 1)  2s( 1) 3d3( 1)
     2S1     2S1     2D1     3S0     2D0
S1
#
#  Step 2. Compute Angular integrals for this case.
#
rm -f Li.out
time Nonh >Li.out << S2
n
y
S2
#
#  Step 3. Run MCHF to obtain radial functions.
#
time Mchf >>Li.out << S3
Li,2D,3.
all
y
y
y
n
S3
echo ' '
echo Move wfn.out to Li2D.w
mv -f wfn.out Li2D.w
#
#  Step 4. Include continuum state in cfile
#
rm -f cont.c
cat >cont.c << S4a
  1s( 2) kdc( 1)                         1.0
     1S0     2D1     2D0
S4a
rm -f cfg.inp 
cat cfg.out cont.c > temp.c
echo ' '
echo Delete extra line in cfile and move it to cfg.inp
rm -f delete
cat >delete << S4b
13d
S4b
sed -f delete temp.c >cfg.inp
#
#  display cfg.inp
#
cat cfg.inp >>Li.out
#
#  Step 5. Run nonh to obtain the angular part of the interaction
#          matrix element.
#
time Nonh >>Li.out << S5
n
n
1,0
S5
echo ' '
echo Move cfg.inp to Li2D.c
mv -f cfg.inp Li2D.c
#
#  Step 6 Run auto to compute continuum orbitals and autoionization
#         rate for 1s2p(2) 2D.
#
time Auto >>Li.out << S6
Li2D
1
Li,2D,3.
-7.2364152
y
S6
#
#  Display the results
#
cat auto.dat >>Li.out
#
rm -f cfg.inp
cat >cfg.inp <<S5
  Li              -7.23694760
                                                          
  1s( 1)  2p( 2)
     2S1     1S0     2S0
  1s( 1)  2p( 2)
     2S1     1D2     2D0
  1s( 1)  2p( 2)
     2S1     3P2     2P0
  1s( 1)  2p( 2)
     2S1     3P2     4P0
  1s( 1)  2p( 1)  3p( 1)
     2S1     2P1     2P1     1P0     2S0
  1s( 1)  2p( 1)  3p( 1)
     2S1     2P1     2P1     1P0     2P0
  1s( 1)  2p( 1)  3p( 1)
     2S1     2P1     2P1     1P0     2D0
  1s( 1)  2p( 1)  3p( 1)
     2S1     2P1     2P1     3P0     2S0
  1s( 1)  2p( 1)  3p( 1)
     2S1     2P1     2P1     3P0     2P0
  1s( 1)  2p( 1)  3p( 1)
     2S1     2P1     2P1     3P0     2D0
  1s( 1)  2p( 1)  3p( 1)
     2S1     2P1     2P1     3P0     4S0
  1s( 1)  2p( 1)  3p( 1)
     2S1     2P1     2P1     3P0     4P0
  1s( 1)  2p( 1)  3p( 1)
     2S1     2P1     2P1     3P0     4D0
  1s( 1)  3p( 2)
     2S1     1S0     2S0
  1s( 1)  3p( 2)
     2S1     1D2     2D0
  1s( 1)  3p( 2)
     2S1     3P2     2P0
  1s( 1)  3p( 2)
     2S1     3P2     4P0
  1s( 1)  2s( 1) 3d1( 1)
     2S1     2S1     2D1     1S0     2D0
  1s( 1)  2s( 1) 3d3( 1)
     2S1     2S1     2D1     3S0     2D0
  1s( 1)  2s( 1) 3d3( 1)
     2S1     2S1     2D1     3S0     4D0
  1s( 1) 3d2( 2)
     2S1     1S0     2S0
  1s( 1) 3d2( 2)
     2S1     1D2     2D0
  1s( 1) 3d2( 2)
     2S1     1G2     2G0
  1s( 1) 3d2( 2)
     2S1     3P2     2P0
  1s( 1) 3d2( 2)
     2S1     3P2     4P0
  1s( 1) 3d2( 2)
     2S1     3F2     2F0
  1s( 1) 3d2( 2)
     2S1     3F2     4F0
  1s( 1)  2s( 2)
     2S1     1S0     2S0
S5
Breit  >>Li.out <<S6
2
n
y
y
S6
rm -f LiBP.c LiBP.w
cp cfg.inp LiBP.c
cp Li2D.w LiBP.w
Ci >>Li.out <<S7
LiBP
y
y
g
y
3
5,1
n
S7
Levels <<S8
LiBP.j
S8
Comp <<S9
LiBP
0.0001
3
S9
cat >> cfg.inp <<S10
  1s( 2) kdc( 1)                            1.0
     1S0     2D1     2D0
S10
Nonh >>Li.out <<S11 
n
y
S11
rm LiBP.c
mv -f cfg.inp LiBP.c
echo Starting Auto
Auto >>Li.out <<S12
LiBP
3
Li,LSJ,3.
y
S12
cat auto.dat >>Li.out
