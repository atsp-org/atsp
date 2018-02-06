#
# CASE II. Prepare <name>.c, <name>.j, <name>.w and int.lst
#          files to run the MCHF_ISOTOPE program for Ni II.
#
BIN=../../bin

# Step 1.  Display a set of configuration state functions. 
#
rm -f cfg.inp Ni.out
cat > cfg.inp << STOP5
 Ni II, 4G
  1s  2s  2p  3s  3p
  3d( 8)  4p( 1)
     3F2     2P1     4G0
*
STOP5
#
# Step 2. Compute angular integrals using the MCHF_NONH program.
#
$BIN/NONH >Ni.out << STOP6
n
y
STOP6
#
# Step 3. Run the MCHF_88 program 
#

$BIN/MCHF >>Ni.out << STOP7
Ni,4G,28.0
all
y
y
y
n
STOP7
#
# Save the results into Ni.c and Ni.w files.
#
echo ' '
echo move wfn.out into Ni.w
mv -f wfn.out Ni.w
echo ' '
echo remove cfg.inp
rm -f cfg.inp
cat > cfg.inp <<STOP1
  Ni                                                         
  1s  2s  2p  3s  3p                                                     
  3d( 8)  4p( 1)
     1D2     2P1     2D0
  3d( 8)  4p( 1)
     1D2     2P1     2F0
  3d( 8)  4p( 1)
     1G2     2P1     2F0
  3d( 8)  4p( 1)
     3P2     2P1     2D0
  3d( 8)  4p( 1)
     3P2     2P1     4P0
  3d( 8)  4p( 1)
     3P2     2P1     4D0
  3d( 8)  4p( 1)
     3F2     2P1     2D0
  3d( 8)  4p( 1)
     3F2     2P1     2F0
  3d( 8)  4p( 1)
     3F2     2P1     4D0
  3d( 8)  4p( 1)
     3F2     2P1     4F0
  3d( 8)  4p( 1)
     3F2     2P1     4G0
*
STOP1
#
# Step 4. Compute angular integrals including relativistic 
#         shift operators using the MCHF_BREIT program.
#
$BIN/BREIT >>Ni.out << STOP8
0
n
y
STOP8
mv -f cfg.inp Ni.c
#
# Step 5. Run the CI program to generate expansion sets. 
#
$BIN/CI >>Ni.out << STOP9
Ni
y
n
11
5,5
n
STOP9
#
# Step 4. Run the ISOTOPE program
#
$BIN/ISO >> Ni.out << STOP10
Ni
3
y
0.01
3
5,5
y
y
56.0,54.0
y
1360,0.16
STOP10
#
# End of Case
#
