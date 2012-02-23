#
# Makefile for MCHF Atomic Structure Package
# Created by C. Froese Fischer, December 1988
#

#  ....   Libraries to be searched
LIB_ANG_COM = -lang -lnjgraf -lcom
LIB_COM_RAD3 = -lcom -lrad3 -lcom
LIB_ANG_COM_RAD3 = -lang -lnjgraf -lcom -lrad3 -lcom
LIB_COM_RAD6 = -lcom -lrad6 -lcom
LIB_ANG_COM_RAD6 = -lang -lnjgraf -lcom -lrad6 -lcom

#  ....  Location of libraries and binaries
LIB = ../bin

#  ....  Compiler for SUN, DEC, SGI, ..
FC = f77

#  ....  Compiler for IBM AIX
#FC = xlf

#  ....  Flags to be added for compilation SGI
#FFLAGS =  -O  -static 
#PFLAGS = -O  -static 

# ....  Flags to be added for compilation on SUN, IBM AIX
FFLAGS = -O
PFLAGS = -O2 
#PFLAGS = -g 

#  .... Flags to be added for compilaton on DEC ALpha
#FFLAGS = -O4 -r4 -non_shared -static -align dcommon -assume noacc -convert big_endian
#PFLAGS = -O2 -r4 -non_shared -static -align dcommon -assume noacc -convert big_endian

#  ....  One of the features of most UNIX is that library routines that 
#  ....  are compiled together are linked together.  For that reason
#  ....  the ANG.f, COM.f, RADn files are split into individual
#  ....  program units and then merged into the corresponding library.

#  ....  Make the ANGular library
 
#    ... Definitions for the ANGular library

ANGSRC  = analy1.f analy2.f cfgin2.f cfgn1.f cfgo1.f cfgtst.f \
          cfp.f cfpd.f cfpf.f cfpp.f find.f gstate.f initt.f \
	  mumdad.f ntab1.f orth.f qsort.f save.f \
	  tensor.f tritst.f vijout.f

ANGOBJS = analy1.o analy2.o cfgin2.o cfgn1.o cfgo1.o cfgtst.o \
          cfp.o cfpd.o cfpf.o cfpp.o find.o gstate.o initt.o \
	  mumdad.o ntab1.o orth.o qsort.o save.o \
	  tensor.o tritst.o vijout.o 

#ANGSRC = ANALY1.f ANALY2.f CFGIN2.f CFGN1.f CFGO1.f CFGTST.f \
#         CFP.f    CFPD.f   CFPF.f    CFPP.f  FIND.f GSTATE.f \
#         INITT.f  MUMDAD.f NTAB1.f   ORTH.f  SAVE.f TENSOR.f \
#         TRITST.f VIJOUT.f qsort.f 
 
#ANGOBJS= ANALY1.o ANALY2.o CFGIN2.o CFGN1.o CFGO1.o CFGTST.o \
#         CFP.o    CFPD.o   CFPF.o   CFPP.o  FIND.o  GSTATE.o \
#         INITT.o  MUMDAD.o NTAB1.o  ORTH.o  SAVE.o  TENSOR.o \
#         TRITST.o  VIJOUT.o  qsort.o 

# libang.a -- make the archive library

libang.a : ANG.f 
	fsplit ANG.f
	$(FC) $(FFLAGS) -c $(ANGSRC)
	ar r libang.a $(ANGOBJS) 
	mv libang.a $(LIB)
	rm $(ANGOBJS) $(ANGSRC) 

# libnjgraf.a -- make the archive and then the random library

libnjgraf.a : NJGRAF.f
	$(FC) $(FFLAGS) -c NJGRAF.f
	ar r libnjgraf.a NJGRAF.o
	mv libnjgraf.a $(LIB)

#  ....  Make the COMmon library
 
#  ...  Definitions for the COMmon library
#

COMSRC  = acurat.f bisect.f bwint.f ca.f cb.f eptr.f factrl.f fk.f gk.f\
	  gracah.f hnorm.f hwf.f inita.f initr.f intact.f lineqn.f\
	  lval.f pack.f reform.f rk.f rme.f sn.f symmeq.f\
          tinvit.f tql2.f trbak1.f tred1.f tred2.f vk.f

COMOBJS = acurat.o bisect.o bwint.o ca.o cb.o eptr.o factrl.o fk.o gk.o\
          gracah.o hnorm.o hwf.o inita.o initr.o intact.o lineqn.o\
          lval.o pack.o reform.o rk.o rme.o sn.o symmeq.o\
          tinvit.o tql2.o trbak1.o tred1.o tred2.o vk.o


#COMSRC = ACURAT.f BISECT.f BWINT.f  CA.f     CB.f     EPTR.f   \
#         FACTRL.f FK.f     GK.f     GRACAH.f HNORM.f  HWF.f    \
#         INITA.f  INITR.f  INTACT.f LINEQN.f LVAL.f   PACK.f   \
#         REFORM.f RK.f     RME.f    SN.f     SYMMEQ.f TINVIT.f \
#         TQL2.f   TRBAK1.f TRED1.f  TRED2.f  VK.f
 
 
#COMOBJS= ACURAT.o BISECT.o BWINT.o  CA.o     CB.o     EPTR.o   \
#         FACTRL.o FK.o     GK.o     GRACAH.o HNORM.o  HWF.o    \
#         INITA.o  INITR.o  INTACT.o LINEQN.o LVAL.o   PACK.o   \
#         REFORM.o RK.o     RME.o    SN.o     SYMMEQ.o TINVIT.o \
#         TQL2.o   TRBAK1.o TRED1.o  TRED2.o  VK.o

# make the archive library

libcom.a : COM.f
	fsplit COM.f
	$(FC) $(FFLAGS) -c $(COMSRC)
	ar r libcom.a $(COMOBJS)
	mv libcom.a $(LIB)
	rm $(COMOBJS) $(COMSRC)


#  ....  Make the RADial library

#   ...  Definitions for the RAD Libraries

RADSRC =  dyk.f ecore.f grad.f hl.f hlc.f quadr.f quads.f\
	  rlshft.f ykf.f zeta.f zk.f

RADOBJS = dyk.o ecore.o grad.o hl.o hlc.o quadr.o quads.o\
	  rlshft.o ykf.o zeta.o zk.o


#RADSRC = DYK.f   ECORE.f  GRAD.f HL.f   HLC.f QUADR.f \
#         QUADS.f RLSHFT.f YKF.f  ZETA.f ZK.f
 
 
#RADOBJS= DYK.o   ECORE.o  GRAD.o HL.o   HLC.o QUADR.o \
#         QUADS.o RLSHFT.o YKF.o  ZETA.o ZK.o
 
#   .... Make the RAD3 library
#
librad3.a : RAD3.f  
	 fsplit RAD3.f
	 $(FC) $(FFLAGS) -c $(RADSRC)
	 ar r librad3.a $(RADOBJS)
	 mv librad3.a $(LIB)
	 rm $(RADSRC) $(RADOBJS)

#   ...  Make the RAD6 library

librad6.a : RAD6.f  
	 fsplit RAD6.f
	 $(FC) $(FFLAGS) -c $(RADSRC)
	 ar r librad6.a $(RADOBJS)
	 mv librad6.a $(LIB)
	 rm $(RADSRC) $(RADOBJS)

#
#  ...  Make the UTILITIES
#

Comp : COMP.o
	$(FC) -o Comp $(FFLAGS)  COMP.o -L$(LIB) $(LIB_ANG_COM)
	mv Comp $(LIB)

Condens : CONDENS.o
	$(FC) -o Condens $(FFLAGS) CONDENS.o 
	mv Condens $(LIB)

Levels : LEVELS.o
	$(FC) -o Levels $(FFLAGS) LEVELS.o 
	mv Levels $(LIB)

Lines : LINES.o
	$(FC) -o Lines $(FFLAGS) LINES.o 
	mv Lines $(LIB)

Plotw:  PLOTW.o
	$(FC) -o Plotw PLOTW.o
	mv Plotw $(LIB)

Printw: PRINTW.o
	$(FC) -o Printw PRINTW.o
	mv Printw $(LIB)
	
Relabel : RELABEL.o
	$(FC) -o Relabel $(FFLAGS) RELABEL.o 
	mv Relabel $(LIB)

Terms : TERMS.o
	$(FC) -o Terms $(FFLAGS) TERMS.o 
	mv Terms $(LIB)

#
#  .... Make the programs. Because global optimizaiton on long
#  .... programs causes swap space problems, these will be 
#  .... compiled here with a lower level of optimization
#


Auto :  AUTO.f
	$(FC) -o Auto $(PFLAGS) AUTO.f -L$(LIB) $(LIB_COM_RAD3)
	mv Auto $(LIB)

Breit : BREIT.f
	$(FC) -o Breit $(PFLAGS) BREIT.f -L$(LIB) $(LIB_ANG_COM)
	mv Breit $(LIB)

Cmchf : CMCHF.f
	$(FC) -o Cmchf $(PFLAGS) CMCHF.f -L$(LIB) $(LIB_COM_RAD3)
	mv Cmchf $(LIB)
	
Ci :    CI.f
	$(FC) -o Ci $(PFLAGS) CI.f -L$(LIB) $(LIB_COM_RAD3)
	mv Ci $(LIB)

Gencl:  GENCL.f
	$(FC) -o Gencl $(FFLAGS) GENCL.f
	mv Gencl $(LIB)

HF :    HF.f
	$(FC) -o HF $(PFLAGS) HF.f
	mv HF $(LIB)

Hfs :   HFS.f
	$(FC) -o Hfs $(PFLAGS) HFS.f -L$(LIB) $(LIB_ANG_COM_RAD3)
	mv Hfs $(LIB)

Iso :   ISO.f
	$(FC) -o Iso $(PFLAGS) ISO.f -L$(LIB) $(LIB_COM_RAD3)
	mv Iso $(LIB)

Lsjtr : LSJTR.f
	$(FC) -o Lsjtr $(PFLAGS) LSJTR.f -L$(LIB) $(LIB_ANG_COM_RAD6)
	mv Lsjtr $(LIB)

Lstr :  LSTR.f
	$(FC) -o Lstr $(FFLAGS) LSTR.f -L$(LIB) $(LIB_COM_RAD6)
	mv Lstr $(LIB)

Mchf :  MCHF.f
	$(FC) -o Mchf $(PFLAGS) MCHF.f -L$(LIB) $(LIB_COM_RAD3)
	mv Mchf $(LIB)

Mltpol: MLTPOL.f
	$(FC) -o Mltpol $(PFLAGS) MLTPOL.f -L$(LIB) $(LIB_ANG_COM)
	mv Mltpol $(LIB)

Nonh :  NONH.f
	$(FC) -o Nonh $(PFLAGS) NONH.f -L$(LIB) $(LIB_ANG_COM)
	mv Nonh $(LIB)

Photo:  PHOTO.f
	$(FC) -o Photo $(PFLAGS) PHOTO.f -L$(LIB) $(LIB_COM_RAD6)
	mv Photo $(LIB)

