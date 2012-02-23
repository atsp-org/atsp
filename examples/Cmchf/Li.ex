#
#   Script (command) file illustrating the use of the CMCHF
#   program in the MCHF Atomic Structure Package.
#
#   Case 1. Li wave function near the 1s.2p(2) resonance
#  
#...To reproduce these results, there should be no wfn.inp file!
#
rm  wfn.inp  
#
#   Step 0.  Obtain a HF 1s for the target
#
time HF >Li.out  <<S0
Li,1S,3.
  1s  

all
y
y
n
n
S0
# echo ' '
# echo Move wfn.out wfn.inp
mv wfn.out wfn.inp
#   
#   Step 1.  Obtain expansion for perturber
#
rm cfg.inp
cat >cfg.inp <<S1
   Li    2D    

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
#   Step 2.  Determine energy expression
#
time Nonh >>Li.out  <<S2
n
y
S2
#
#   Step 3.  Determine the radial function for the MCHF approximation
#
time Mchf >>Li.out  <<S3
Li,2D,3.
=6
y
y
y
n
S3
# echo  ' '
# echo Move cfg.out to al2P.c
rm li2Dp.c li2Dp.w
mv cfg.out li2Dp.c
# echo  ' '
# echo Move wfn.out to al2P.w
mv wfn.out li2Dp.w
#
#   Step 4.  Determine expansion  -- perturber and continuum
#
rm cfg.inp
cat >cfg.inp <<S4
   Li    2D        -7.2364152

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
  1s( 2) kd4( 1)                           1.000000
     1S0     2D1     2D1
S4
# echo  ' '
#   
#   Step 5.  Get energy expression for the new expansion
#
time Nonh >>Li.out  <<S5
n
y
S5
#
#   Step 6.  Determine continuum wavefunctions for an 
#            Energy Range
#
cp li2Dp.w wfn.inp
time Cmchf >>Li.out <<S6
Li,2D,3.
=1
y
y
y
4.0, .01, 4.2
n
S6
cat phase.out >>Li.out
