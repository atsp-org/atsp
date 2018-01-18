#! /bin/sh

set -e
#
#   Script (command) file illustrating the use of the Nonh program
#
#   Case 1.  
#
#..... Step 1. Remove cfg.inp and then create cfg.inp
#
rm -f cfg.inp Nonh1.out
cat >cfg.inp <<S1
 Test Number 1

  1s( 2) 5s2( 2)
     1S0     1S0     1S0
 3s1( 2) 4s1( 2)
     1S0     1S0     1S0
S1
#
#..... Step 2. Obtain the energy expression
#
time ../../src/NONH >Nonh1.out <<S2
y
y
S2
#
# echo Display the int.lst file produced
#
cat int.lst >>Nonh1.out
#
#   Case 2.
#
#..... Step 1.  Remove cfg.inp and then create cfg.inp
#
rm cfg.inp
cat >cfg.inp <<S1
 Test 2

  1s( 1) 2p2( 1) 3d2( 1)
     2S1     2P1     2D1     1P0     2P0
  1s( 1) 2p2( 1) 3d2( 1)
     2S1     2P1     2D1     3P0     2P0
 2si( 1) 3pi( 1) 3d1( 1)
     2S1     2P1     2D1     1P0     2P0
 2si( 1) 3pi( 1) 3d1( 1)
     2S1     2P1     2D1     3P0     2P0
S1
#
#..... Step 2. Obtain the energy expression
#
time ../../src/NONH >Nonh2.out<<S2
y
y
S2
#
# echo Display the int.lst file produced
#
cat int.lst  >>Nonh2.out
