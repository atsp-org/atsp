#
#  Case 2.  M1 nd E2 LSJTR transitions in Carbon-like Argon.
#           This illustrates the case where the same file is
#           used for initial states and final states
#
#   Step 7.  Generate the configuration state list for the 2p(2) LSJ 
#            terms that include single and double replacements
#
rm -f Argon.out
time Gencl >Argon.out <<S7

 Carbon-like Argon
 1s  
2s(2)2p(2)
2p(4)


s
3s,3p,3d
1
2

S7
echo '  '
echo 
cat cfg.inp
echo '  '
echo Move cfg.inp to argon.c
# When the system does not allow two files to be connected
# to two units, we need to make two copies
cp -f cfg.inp argon1.c
mv -f cfg.inp argon2.c
#
#   Step 8. Determine the expressions for the E2 and M1
#            transitions operators 
#
time Mltpol >>Argon.out <<S8
argon1
argon2
y
E2
M1
*
S8
echo '  '
echo Display the mltpol.lst file produced
cat mltpol.lst >>Argon.out
