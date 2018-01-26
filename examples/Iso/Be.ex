#  
#  Script (command) file illustrating the use of the MCHF
#  Atomic Structure Package for the study of the isotope 
#  shift of 1s(2)2s(1)3s(1)_1S in Be I and 3d(8)4p(1)_4G
#  in Ni II.  
#
#  CASE I. Prepare <name>.w, <name>.c and int.lst files 
#          to run the MCHF_ISOTOPE program for Be I.
#
BIN=../../bin

# ------------------------------------------
#  Step 1. Run the HF program to get initial wavefunctions.
# ------------------------------------------
rm -f wfn.inp
$BIN/HF << S1
Be,1S,4.
  1s
2s(1)3s(1)
all
y
y
n
n
S1

# -------------------------------------------------
#  Step 2. Specify set of configuration state functions.
# -------------------------------------------------

rm -f cfg.inp Be.out
cat > cfg.inp << STOP1
   Be    1S 
  1s
  2s( 1)  3s( 1)
     2S1     2S1     1S0
  2s( 2)       
     1S0
  2p( 2)      
     1S0
  3p( 2)     
     1S0
  3d( 2)    
     1S0
STOP1

# -------------------------------------------
# Step 3. Compute angular integrals using the MCHF_nonh program. 
# --------------------------------------------
$BIN/NONH > Be.out << STOP2
n
y
STOP2

# --------------------------------------------
# Step 4. Run the MCHF_88 program to obtain radial functions
#         and mixing coefficients.
# --------------------------------------------
mv -f wfn.out wfn.inp
echo "*" >> cfg.inp

$BIN/MCHF >>Be.out << STOP3
Be,1S,4.0
all
y
y
y
STOP3
#
# Save the results into Be.c and Be.w files. 
#
echo ' '
echo move cfg.out into Be.c
echo move wfn.out into Be.w
mv -f cfg.out Be.c
mv -f wfn.out Be.w
echo ' '

# -----------------------------------------
# Step 5. Run the MCHF_ISOTOPE program
# -----------------------------------------
$BIN/ISO >> Be.out << STOP4
Be
1
y
0.001
3
y
y
9.0,10.0
n
STOP4
