#
#   Script (command) file illustrating the use of the CMCHF
#   program in the MCHF Atomic Structure Package.
#
#   1skp wave function including correltion
#  
#   Step 0.  Get expansion for initial state
#
rm -f wfn.inp  cfg.inp H-.out
cat >cfg.inp <<S0
   H-    
        
  1s( 2)
     1S0
  2s( 2)
     1S0
  2p( 2)
     1S0
  3s( 2)
     1S0
  3p( 2)
     1S0
  3d( 2)
     1S0
S0
#
#   Step 1.  Determine energy expression
#
time Nonh >H-.out  <<S1
n
y
S1
#
#   Step 2.  Determine the radial function for the MCHF approximation
#            Note the two-step process
#
time Mchf >>H-.out  <<S2
H-,1S,1.
1s,2s,2p
n 
0.3,0,1,0
-1,0,3,0
-1,0,3,0
y
y
n
S2
mv -f cfg.out cfg.inp
mv -f wfn.out wfn.inp
time Mchf >>H-.out  <<S2
H-,1S,1.
all
y 
y
y
n
S2
# echo  ' '
# echo Move cfg.out to H1S.c
mv -f cfg.out H1S.c
# echo  ' '
# echo Move wfn.out to H1S.w
mv -f wfn.out H1S.w
#
#   Step 3.  Get expansion for the final state
#            -- perturber and continuum
#
rm -f cfg.inp
cat >cfg.inp <<S3
   H-    1P        -0.5000000 

  2s( 1) 9p1( 1)             
     2S1     2P1     1P0    
 2ph( 1) 9s1( 1)           
     2P1     2S1     1P0  
 2ph( 1) 9d1( 1)         
     2P1     2D1     1P0                                                
 3sh( 1) 9p2( 1)        
     2S1     2P1     1P0                                                
 3ph( 1) 9s2( 1)       
     2P1     2S1     1P0                                                
 3ph( 1) 9d2( 1)      
     2P1     2D1     1P0                                                
 3dh( 1) 9p3( 1)     
     2D1     2P1     1P0                                                
 3dh( 1) 9f1( 1)    
     2D1     2F1     1P0                                                
  1s( 1) kpc( 1)                          1.000000
     2S1     2P1     1P0                                                
*****
 2ph 9p2
 2ph 9p3
 3sh 9s2
 3dh 9d2
S3
# echo  ' '
#   
#   Step 4.  Get energy expression for the new expansion
#
time Nonh >>H-.out  <<S4
n
y
S4
#
#   Step 5.  Determine continuum wavefunctions for an 
#            Energy Range
#
time Cmchf >>H-.out <<S5
H-,1P,1.
9p1,9s1,9d1,9p2,9s2,9d2,9p3,9f1,kpc
y
y
y
.014,-.001,.002
n
S5
cat phase.out >>H-.out
