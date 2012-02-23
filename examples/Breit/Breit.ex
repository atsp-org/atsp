#
#  Script (command) file for the BREIT test run data
#
#  Case 1.   Breit-Pauli interactions for the odd parity configuration
#            states of the Al-like sequence
#   
#   Step 1.  Create a configuration list
#           
rm -f cfg.inp Breit.out
cat >cfg.inp <<S1
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
S1
#
#   Step 2.  Derive the energy expression for the Breit-Pauli 
#            interaction matrix
time Breit >Breit.out  <<S2
1
y
y
y
S2
echo  ' '
echo Display the int.lst file produced
