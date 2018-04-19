#
#  Script (command) file for the BREIT test run data
#
#  Case 1.   Breit-Pauli interactions for the odd parity configuration
#            states of the Al-like sequence
#   
#   Step 1.  Create a configuration list
#
BIN=../../bin

rm -f cfg.inp Breit.out
cat >cfg.inp <<ST1
   al                                                        
  1s  2s  2p
  3s( 2)  3p( 1)
     1S0     2P1     2P0
  3p( 3)
     2P1
  3p( 3)
     2D3
  3p( 3)
     4S3
ST1
#
#   Step 2.  Derive the energy expression for the Breit-Pauli 
#            interaction matrix
time $BIN/BREIT >Breit.out  <<ST2
1
y
y
y
ST2
echo  ' '
echo Display the int.lst file produced
