*    SUBROUTINE PACKAGE NJGRAF
*
*           C O P Y R I G H T -- 1994
*
*    A. Bar-Shalom and Marcel Klapisch
*    
*    Computer Physics Communications, Vol. 50, 375 (1988)
*
*    This version contains the corrections to errors found
*    in the published version of this package.  CFF
*
*    Gensum correct December, 1996
************************************************************************
      SUBROUTINE BUBBLE(JPOL,FAIL)
*
*    REDUCES A CIRCUIT OF ORDER 2,GIVING DELTA FUNCTION AND PHASE
*    FACTORS.
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,KFL2=12,KFL2A=2*KFL2,KFL2B=4*KFL2,
     +          KFL2C=2*KFL2+2,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20)
*
      LOGICAL FAIL
*
      INTEGER ARR,TAB1
      COMMON/GRAPH/JDIAG(KFL2B,3),ARR(KFL2B,3),TAB1(KFL1,2),IL(KFL2B),
     + IH(KFL2B),NPOINT(KFL2C),NBNODE,IFIRST,ILAST,IPARTS,IPARTL,NPART,
     + ICROSS,NFREE,ITFREE(KFL2A),NFIN,NC
      LOGICAL SUMVAR
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),KW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
*
      CHARACTER*6 NAME,NAMSUB
      COMMON/NAM/NAMSUB
      DATA NAME/'BUBBLE'/
*
      NAMSUB=NAME
      K2=2
      K23=3
      I1=1
      I2=1
      IT1=NPOINT(1)
      IT2=NPOINT(2)
      IF(IT2.NE.ILAST)GO TO 2
      IF(IT1.EQ.IFIRST)GO TO 1
      IT2=IT1
      IT1=ILAST
    1 I1=-1
      K23=2
      I2=2
    2 CALL PHASE(IT1,JDIAG,KFL2B)
      K=IABS((3*ARR(IT2,1)+2*ARR(IT2,2)+ARR(IT2,3))/2)+1
      IF(K.NE.4)CALL PHASE2(JDIAG(IT2,K))
      IF(NBNODE.EQ.2)GO TO 7
      IL1=IL(IT2)+I1
      IT=IH(IL1)
      ARR(IT,K23)=ARR(IT1,K23)
      L=JDIAG(IT1,K23)
      L1=JDIAG(IT,K23)
      JDIAG(IT,K23)=L
      IF(JPOL.NE.1)GO TO 3
      MP=MP-1
      KW(2,JWC)=L
      J6(J6C-1)=L
      J6(J6C)=L
      IF(K.EQ.2)J8(J8C)=L
      GO TO 4
    3 CALL DELTA(L,L1,FAIL)
      IF(FAIL)GO TO 7
    4 TAB1(L,I2)=IT
      IF(IT1.EQ.ILAST)GO TO 6
      IF(IT2.NE.ILAST)GO TO 9
      TAB1(L,1)=IH(2)
      IL1=2
      K2=1
    9 DO 5 I=IL1,NBNODE
      IT=IH(I)
      IL(IT)=I-K2
      IH(I-K2)=IT
    5 CONTINUE
    6 J9(J9C+1)=L
      J9C=J9C+2
      J9(J9C)=L
    7 RETURN
      END
************************************************************************
      SUBROUTINE CHANGE(L,K)
*
*     EXCHANGES THE FREE ENDS IN EITHER FIRST OR LAST TRIAD OF JDIAG.
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,KFL2=12,KFL2A=2*KFL2,KFL2B=4*KFL2,
     +          KFL2C=2*KFL2+2)
*
      INTEGER ARR,TAB1
      COMMON/GRAPH/JDIAG(KFL2B,3),ARR(KFL2B,3),TAB1(KFL1,2),IL(KFL2B),
     + IH(KFL2B),NPOINT(KFL2C),NBNODE,IFIRST,ILAST,IPARTS,IPARTL,NPART,
     + ICROSS,NFREE,ITFREE(KFL2A),NFIN,NC
*
      CALL PHASE(L,JDIAG,KFL2B)
      JP=JDIAG(L,K)
      JDIAG(L,K)=JDIAG(L,1)
      JDIAG(L,1)=JP
      JAR=ARR(L,K)
      ARR(L,K)=ARR(L,1)
      ARR(L,1)=JAR
      RETURN
      END
************************************************************************
c
c     ******************************
      subroutine chklp1(fail)
      implicit real*8(a-h,o-z)
c     ******************************
c
c
c  ***this routine checks if there are active triads with two identical
c  ***arguments.-This is a loop of order 1 ("lollypop").
c  *** The other argument must then be zero.i.e.j1(j)=1 in 2j+1 notation
c  ***Suppression of the loop introduces  factors and phases.
c  ***Two triads become inactive. all this is performed by invoking sub.
c  ***zero with first arg. = 1.
c  ***written by MK 4/92 for correcting an error detected by CFF
c
*
      PARAMETER(KFL1=60,KFL2=12,KFL2A=2*KFL2,KFL2B=4*KFL2,
     +          KFL2C=2*KFL2+2,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20)
c
      logical fail
c
      character*6 name,namsub 
c
      common/nam/namsub
*
      LOGICAL FREE
      COMMON/COUPLE/M,N,J1(KFL1),J2(KFL2,3),J3(KFL2,3),FREE(KFL1)
*
      LOGICAL TABS
      INTEGER ARROW
      COMMON/TREE/J23(KFL2A,3),ARROW(KFL2A,3),LINE(KFL1,2),LCOL(KFL1,2),
     +  TABS(KFL2A),NBTR
c
      data name/'chklp1'/
c
      namsub= name
*     the following line corrected June 12, 1992
      nbtr1=2*(n-1)
      do 10 l=1,nbtr1
      if(.not.tabs(l)) then
         jdif=0
         if(j23(l,1).eq.j23(l,2)) then
            jdif=j23(l,3)
         else if(j23(l,1).eq.j23(l,3))then
            jdif=j23(l,2)
         else if(j23(l,2).eq.j23(l,3))then
            jdif=j23(l,1)
         endif
         if(jdif.ne.0)then
c..........putting the link to 0.sub zero changes nbtr
          fail=.false.
          if (j1(jdif).ne.1.and..not.free(jdif))then
            fail=.true.
c           write(*,*)jdif,' should be 0.Is :',j1(jdif),'recup ->0'
            return
            else
             call zero(1,jdif,fail)
             if (fail)return
           endif
          endif
      endif
  10  continue 
      if(jdif.ne.0)call printj(name,4)
      return
      end
************************************************************************
      SUBROUTINE CHVAR(JP,NBC,KBC,JT,JINV,NSUM)
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
*     CHANGE THE ORDER OF SUMMATION VARIABLE TO BE ABLE TO PERFORM
*     SEPARATELY THE SUMMATIONS IN GENSUM.
*
      LOGICAL JT(NSUM)
      DIMENSION JP(NBC),JINV(NSUM)
      KB=KBC+1
      IF(KB.GT.NBC)GO TO 2
      DO 1 I=KB,NBC
      JK=JP(I)
      IF(.NOT.JT(JK))GO TO 1
      KBC=KBC+1
      JP(I)=JP(KBC)
      JP(KBC)=JINV(JK)
    1 CONTINUE
    2 RETURN
      END
************************************************************************
      SUBROUTINE CUT1L(FAIL)
*

*     CUT ON ONE LINE,THAT WAS LEFT AS A FREE END IN JDIAG.PUTS
*     CORRESPONDING DELTA IN J23.
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,KFL2=12,KFL2A=2*KFL2,KFL2B=4*KFL2,
     +          KFL2C=2*KFL2+2,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20)
*
      LOGICAL FAIL
*
      INTEGER ARR,TAB1
      COMMON/GRAPH/JDIAG(KFL2B,3),ARR(KFL2B,3),TAB1(KFL1,2),IL(KFL2B),
     + IH(KFL2B),NPOINT(KFL2C),NBNODE,IFIRST,ILAST,IPARTS,IPARTL,NPART,
     + ICROSS,NFREE,ITFREE(KFL2A),NFIN,NC
*
      LOGICAL SUMVAR
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),KW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
      LOGICAL FREE
      COMMON/COUPLE/M,N,J1(KFL1),J2(KFL2,3),J3(KFL2,3),FREE(KFL1)
*
      CHARACTER*6 NAME
      DATA NAME/'CUT1L '/
*
      IT=ITFREE(1)
      J0=JDIAG(IT,1)
      CALL DELTA (J0,M,FAIL)
      IF(FAIL)GO TO 2
      CALL DELTA(JDIAG(IT,3),JDIAG(IT,2),FAIL)
      IF(FAIL)GO TO 2
      JDIAG(IT+1,3)=JDIAG(IT,3)
      IF(ARR(IT,2)-ARR(IT,3))4,3,5
    3 ARR(IT+1,3)=1
      ARR(IT-1,2)=-1
      GO TO 5
    4 ARR(IT+1,3)=-1
      ARR(IT-1,2)=1
    5 J9C=J9C+1
      J9(J9C)=JDIAG(IT,3)
      J=2
      CALL ZERO(J,J0,FAIL)
      IF(FAIL)GO TO 2
      IL1=IL(IT+1)
      DO 1 I=IL1,NBNODE
      IT=IH(I)
      ILP=I-1
      IL(IT)=ILP
      IH(ILP)=IT
    1 CONTINUE
      NBNODE=NBNODE-1
    2 CALL PRINTJ(NAME,12)
      RETURN
      END
************************************************************************
      SUBROUTINE CUT2L(FAIL)
*
*     CUT ON TWO LINES THAT WERE LEFT AS FREE ENDS IN JDIAG.PUTS
*     CORRESPONDING DELTA IN J23.
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,KFL2=12,KFL2A=2*KFL2,KFL2B=4*KFL2,
     +          KFL2C=2*KFL2+2,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20)
*
      LOGICAL FAIL
*
      INTEGER ARR,TAB1
      COMMON/GRAPH/JDIAG(KFL2B,3),ARR(KFL2B,3),TAB1(KFL1,2),IL(KFL2B),
     + IH(KFL2B),NPOINT(KFL2C),NBNODE,IFIRST,ILAST,IPARTS,IPARTL,NPART,
     + ICROSS,NFREE,ITFREE(KFL2A),NFIN,NC
*
      LOGICAL TABS
      INTEGER ARROW
      COMMON/TREE/J23(KFL2A,3),ARROW(KFL2A,3),LINE(KFL1,2),LCOL(KFL1,2),
     +  TABS(KFL2A),NBTR
*
      LOGICAL SUMVAR
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),KW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
*
      CHARACTER*6 NAME
      DATA NAME/'CUT2L '/
*
      IT1=ITFREE(1)
      IT2=ITFREE(2)
      JT1=JDIAG(IT1,1)
      JT2=JDIAG(IT2,1)
      CALL DELTA(JT1,JT2,FAIL)
      IF(FAIL) GO TO 1
      IF(ARR(IT1,1).EQ.ARR(IT2,1))CALL PHASE2(JT1)
      ARR(IT2,1)=-ARR(IT1,1)
      JDIAG(IT2,1)=JT1
      TAB1(JT1,2)=IT2
      J9(J9C+1)=JT1
      J9C=J9C+2
      J9(J9C)=JT1
      CALL OTHERJ(0,JT1,L1,LC1,K1)
      CALL OTHERJ(0,JT2,L2,LC2,K2)
      J23(L2,LC2)=JT1
      LINE(JT1,K1)=L2
      LCOL(JT1,K1)=LC2
      ARROW(L2,LC2)=-ARROW(L1,LC1)
    1 CALL PRINTJ(NAME,12)
      RETURN
      END
************************************************************************
      SUBROUTINE CUTNL(FAIL)
*
*     THIS SUBROUTINE  EXAMINES THE CASE WHERE THERE ARE MORE THAN
*     TWO FREE ENDS,BUT THEY ARE CONTIGUOUS,SO THAT THE GRAPH CAN
*     BE CUT WITHOUT DESTROYING THE FLAT STRUCTURE.
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,KFL2=12,KFL2A=2*KFL2,KFL2B=4*KFL2,
     +          KFL2C=2*KFL2+2,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20)
*
      INTEGER ARROW
      LOGICAL TABS
      COMMON/TREE/J23(KFL2A,3),ARROW(KFL2A,3),LINE(KFL1,2),LCOL(KFL1,2),
     +  TABS(KFL2A),NBTR
*
      INTEGER ARR,TAB1
      COMMON/GRAPH/JDIAG(KFL2B,3),ARR(KFL2B,3),TAB1(KFL1,2),IL(KFL2B),
     + IH(KFL2B),NPOINT(KFL2C),NBNODE,IFIRST,ILAST,IPARTS,IPARTL,NPART,
     + ICROSS,NFREE,ITFREE(KFL2A),NFIN,NC
*
      COMMON/KEEP/JKP(2,3),JARR(2,3),IT2,IT3,IT5
*
      LOGICAL SUMVAR
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),KW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
*
      LOGICAL FAIL
*
      CHARACTER*6 NAME
      DATA NAME/'CUTNL '/
*
      NTF=ITFREE(NFREE)-ITFREE(1)
      IF(NTF.GT.NFREE)GO TO 8
      IT2=ITFREE(1)
      IT3=ITFREE(NFREE)
      IT1=IT2-1
      IT4=IT3+1
      IF(NTF.EQ.NFREE)GO TO 2
      JT=JDIAG(IT2,3)
      CALL DELTA(JT,JDIAG(IT3,2),FAIL)
      IF(FAIL) GO TO 9
      IF(ARR(IT2,3).NE.ARR(IT3,2))GO TO 1
      CALL PHASE2(JT)
      ARR(IT2,3)=-ARR(IT2,3)
      ARR(IT1,2)=-ARR(IT1,2)
    1 JDIAG(IT3,2)=JT
      JDIAG(IT4,3)=JT
      J9(J9C+1)=JT
      J9C=J9C+2
      J9(J9C)=JT
      NBTR=NBTR+NFREE
      IT5=0
      GO TO 6
    2 NFR=0
      DO 3 IT5=IT2,IT3
      NFR=NFR+1
      IF(ITFREE(NFR).GT.IT5)GO TO 4
    3 CONTINUE
    4 JKP(1,1)=JDIAG(IT5,1)
      JARR(1,1)=-ARR(IT5,1)
      JKP(1,2)=JDIAG(IT2,3)
      JARR(1,2)=-ARR(IT2,3)
      JKP(1,3)=JDIAG(IT3,2)
      JARR(1,3)=-ARR(IT3,2)
      DO 5 J=1,3
      JKP(2,J)=JDIAG(IT5,J)
    5 JARR(2,J)=ARR(IT5,J)
      JDIAG(IT5,2)=JDIAG(IT3,2)
      ARR(IT5,2)=ARR(IT3,2)
      JDIAG(IT5,3)=JDIAG(IT2,3)
      ARR(IT5,3)=ARR(IT2,3)
      ILP=IL(IT2)
      IL(IT5)=ILP
      IH(ILP)=IT5
      NBTR=NBTR+NFREE+2
      CALL PHASE(IT5,JDIAG,KFL2B)
*Changed according to message from M. Klapisch, Nov 5, 1993
*     K=IABS((3*ARR(IT5,1)+2*ARR(IT5,2)+ARR(IT5,3))/2)+1
      K=IABS(3*ARR(IT5,1)+2*ARR(IT5,2)+ARR(IT5,3))/2+1
      IF(K.NE.4) CALL PHASE2(JDIAG(IT5,K))
    6 IL1=IL(IT4)
      DO 7 I=IL1,NBNODE
      IT=IH(I)
      ILP=I-NFREE
      IL(IT)=ILP
      IH(ILP)=IT
    7 CONTINUE
      NBNODE=NBNODE-NFREE
      NFIN=0
      GO TO 8
    9 FAIL=.TRUE.
    8 CALL PRINTJ(NAME,8)
      RETURN
      END
************************************************************************
      SUBROUTINE DELTA(JA,JB,FAIL)
*
*     TEST FOR DELTA(JA,JB).IF THEY ARE SUMMATION VARIABLES,THE SECOND
*     IS CHANGED INTO THE FIRST EVERYWHERE.IF THEY ARE FIXED,THEIR
*     VALUE IS CHECKED,AND FAIL PUT TO .TRUE. IF THEY DIFFER.
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,KFL2=12,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20)
*
      LOGICAL FAIL
*
      LOGICAL CUT
      COMMON/CUTDIG/CUT
      LOGICAL SUMVAR
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),KW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
*
      COMMON/DIM/J6CC,J7CC,J8CC,J9CC,JWCC,JDELC
*
*
      LOGICAL FREE
      COMMON/COUPLE/M,N,J1(KFL1),J2(KFL2,3),J3(KFL2,3),FREE(KFL1)
*
      IF(IBUG3.EQ.1)PRINT 1000,JA,SUMVAR(JA),JB,SUMVAR(JB)
 1000 FORMAT(/2X,'FROM DELTA',2X,'JA=',I2,L2,5X,'JB=',I2,L2)
      IF(SUMVAR(JA).AND.SUMVAR(JB))GO TO 2
      IF(FREE(JA).OR.FREE(JB))GO TO 1
      IF(J1(JA).NE.J1(JB))FAIL=.TRUE.
      CUT=.TRUE.
      GO TO 18
    1 JDEL=JDEL+1
      LDEL(JDEL,1)=JA
      LDEL(JDEL,2)=JB
      SUMVAR(JA)=.FALSE.
      SUMVAR(JB)=.FALSE.
      GO TO 18
    2 IF(J6C.EQ.J6CC)GO TO 4
      J61=J6CC+1
      DO 3 I=J61,J6C
    3 IF(J6(I).EQ.JB)J6(I)=JA
    4 IF(J7C.EQ.J7CC)GO TO 6
      J71=J7CC+1
      DO 5 I=J71,J7C
    5 IF(J7(I).EQ.JB)J7(I)=JA
    6 IF(J8C.EQ.J8CC)GO TO 8
      J81=J8CC+1
      DO 7 I=J81,J8C
    7 IF(J8(I).EQ.JB)J8(I)=JA
    8 IF(J9C.EQ.J9CC)GO TO 10
      J91=J9CC+1
      DO 9 I=J91,J9C
    9 IF(J9(I).EQ.JB)J9(I)=JA
   10 IF(JWC.EQ.JWCC)GO TO 15
      JW1=JWCC+1
      DO 14 I=JW1,JWC
      DO 13 J=1,6
   13 IF(KW(J,I).EQ.JB)KW(J,I)=JA
   14 CONTINUE
   15 IF(JDEL.EQ.JDELC)GO TO 18
      JDEL1=JDELC+1
      DO 17 I=JDEL1,JDEL
      DO 16 J=1,2
   16 IF(LDEL(I,J).EQ.JB)LDEL(I,J)=JA
   17 CONTINUE
      SUMVAR(JB)=.FALSE.
   18 RETURN
      END
************************************************************************
      subroutine diagrm(jump)
************************************************************************
*
*     this subroutine builds up a flat diagram from the triads j23 and
*  ***places them in jdiag.arrows are in arr (integer).the diagram is
*  ***built so as to maximize the number of triads involved,within  a
*  ***one-step-forward-check process.if the diagram does not
*  ***include all the nbtr triads,it will have 'free ends'.jdiag has
*  ***dimension double that of j23,because the path may proceed either
*  ***way.
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,KFL2=12,KFL2A=2*KFL2,KFL2B=4*KFL2,
     +          KFL2C=2*KFL2+2,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20)
*
      LOGICAL SUMVAR
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),KW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
      LOGICAL TABS
      INTEGER ARROW
      COMMON/TREE/J23(KFL2A,3),ARROW(KFL2A,3),LINE(KFL1,2),LCOL(KFL1,2),
     +  TABS(KFL2A),NBTR
      INTEGER ARR,TAB1
      COMMON/GRAPH/JDIAG(KFL2B,3),ARR(KFL2B,3),TAB1(KFL1,2),IL(KFL2B),
     + IH(KFL2B),NPOINT(KFL2C),NBNODE,IFIRST,ILAST,IPARTS,IPARTL,NPART,
     + ICROSS,NFREE,ITFREE(KFL2A),NFIN,NC
*
      COMMON/BUILD/IAL(KFL2B),IF1,IF2,NODE
*
      LOGICAL FREE
      COMMON/COUPLE/M,N,J1(KFL1),J2(KFL2,3),J3(KFL2,3),FREE(KFL1)
*
      character*6 name
      data name/'diagrm'/
*
*
*
*  ***initialization
*
      if(jump .gt. 2) go to 17
      if(jump .lt. 2) nb=0
    1 nb=nb+1
      if(tabs(nb))go to 1
      node=nbtr
      ilast=nbtr
*
*      write(6,*) 'nb = ', nb
      do 2 j=1,3
        jdiag(node,j)=j23(nb,j)
        arr(node,j)=arrow(nb,j)
    2 continue
*
      tabs(nb)=.true.
*
      do 15 i=1,mp
        ial(i)=0
   15 continue
*
      if1=jdiag(node,1)
      if2=jdiag(node,3)
      ial(if1)=1
      ial(if2)=1
   17 ntime=0
      i1=1
      k1=1
      k2=2
      k3=3
    3 jb=jdiag(node,k2)
      call otherj(0,jb,l,lc,kp)
      call neibor(lc,l1,l2)
*  november 22 1989 check consistency of triads
      if(tabs(l))stop' building diagram impossible '
      call way(l,l1,l2,ich,nd)
      node=node+i1
      tabs(l)=.true.
      jdiag(node,k3)=j23(l,lc)
      arr(node,k3)=arrow(l,lc)
      ict=ich*i1
*
      if (ich .le. 0) then
        lp=l1
        l1=l2
        l2=lp
      endif
*
      if (ict .le. 0) call phase(l,j23,kfl2a)
      jdiag(node,k1)=j23(l,l1)
      arr(node,k1)=arrow(l,l1)
      jdiag(node,k2)=j23(l,l2)
      arr(node,k2)=arrow(l,l2)
      j=j23(l,l1)
      ial(j)=ial(j)+1
      j=j23(l,l2)
      ial(j)=ial(j)+1
      if(nd.lt.1)go to 3
      ntime=ntime+1
      ilast=max0(node,ilast)
      ifirst=min0(node,nbtr)
      nbp=ial(if1)+ial(if2)
      if (nbp .gt. 3 .or. ntime .gt. 1) then
        nbnode=ilast-ifirst+1
        nbtr=nbtr-nbnode
*
*  ***definition of free ends and other quantities.
*
        call intab
        call printj(name,mtriad)
        go to 50
      endif
*
      if (nbp .gt. 2) then
        if (ial(if1) .le. ial(if2)) then
          jt=jdiag(nbtr,1)
          jar=arr(nbtr,1)
          jdiag(nbtr,1)=jdiag(nbtr,3)
          arr(nbtr,1)=arr(nbtr,3)
          jdiag(nbtr,3)=jt
          arr(nbtr,3)=jar
          call phase(nbtr,jdiag,kfl2b)
      endif
      endif
*
      node=nbtr
      i1=-1
      k2=3
      k3=2
      go to 3
*
   50 return
      end
************************************************************************
      SUBROUTINE GENSUM(J6C,J7C,J8C,J9C,JWC,J6,J7,J8,J9,JW,JDEL,
     +            LDEL,SUMVAR,MP,
     +            J6P,J7P,J8P,J9P,JWORD,NLSUM,NBJ,NB6J,
     +            K6CP,K7CP,K8CP,K9CP,JSUM4,JSUM5,JSUM6,INV6J,
     +            RECUP)
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
*  CARRIES OUT THE SUMMATION OVER COEFFICIENTS DEFINED BY THE ARRAYS
*  J6,J7,J8,LDEL AND JW TO GIVE RECUP
*  THE ENTRY IS EITHER MADE FROM NJGRAF OR DIRECTLY ASSUMING THAT THE
*  ARRAYS J6,...,JW HAVE ALREADY BEEN DETERMINED BY A PREVIOUS
*  ENTRY TO NJGRAF AND THAT THE SUMMATION IS REQUIRED FOR ANOTHER SET
*  OF J VALUES DEFINED BY THE ARRAY J1
*
*  RECUP IS THE RECOUPLING COEFFICIENT
*
*  SUBROUTINE CALLED: GRACAH
*
      PARAMETER(KFL1=60,KFL2=12,KFL2A=2*KFL2,KFL2B=4*KFL2,
     +          KFL2C=2*KFL2+2,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20,
     +          KFLS=12,KFLN=10,KFLV=40)
*
      LOGICAL SUMVAR
      DIMENSION J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),JW(6,KFLW),LDEL(KFLW,2),SUMVAR(KFL1)
*
      DIMENSION J6P(KFLV),J7P(KFLV),J8P(KFLV),J9P(KFLV),JWORD(6,KFLW),
     + NBJ(KFLN),NB6J(KFLN),K6CP(KFLN),K7CP(KFLN),K8CP(KFLN),K9CP(KFLN),
     + JSUM6(KFLS),JSUM4(KFLS,KFLW),JSUM5(KFLS,KFLW),INV6J(KFLW)
*
      LOGICAL LDIAG,NOEL
      DIMENSION MAT(KFLS,KFLS),JMNP(5),JMXP(5),NOEL(KFLS),MAXLP(KFLS),
     + IST(6),JSUM2(KFLS),JSUM3(KFLS),JSUM8(KFLS),JSUM(2,KFLW),
     + JWTEST(KFLW),WSTOR(KFLW),IPAIR(2,2),LDIAG(KFLS)
      DIMENSION J12(4,KFLS,KFLS)
      DIMENSION XJ1(KFL1)
*
      LOGICAL FREE
      COMMON/COUPLE/M,N,J1(KFL1),J2(KFL2,3),J3(KFL2,3),FREE(KFL1)
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON /FACTS / GAM(100)
*
      parameter(mxcsvr=6)
c
      common/inform/ ird , ipd
c
      dimension jlow(2),jhig(2)
      dimension jlimit(2,2,6)
c
      logical already
c
      data jlimit/2,3,5,6,
     +1,4,5,6,
     +1,4,6,5,
     +2,3,6,5,
     +1,3,2,4,
     +1,2,3,4/
      data zero,one,two /0.D0,1.D0,2.D0/
      data epsil/1.D-12/
      data maxincr/30/
c
c  ***format statements used in gensum
c
  302 format('   sum nr.',i3)
  303 format(' no summation. recoupling coefficient=',g15.8)
  304 format(' recoupling coefficient=',g15.8)
  305 format(i6,2x,6f5.1,10x,g15.8)
  306 format(' number of independent sums:',i3)
  307 format(' sum nr.',i2,' sum value=',g15.8,' recup=',g15.8)
  308 format('recoupling coefficient zero because of unsatisfied delta',
     +'j1',i3,'=',i3,'j2',i3,'=',i3)
  309 format('   not involving summation variable')
  310 format(' summation variable reached maximum increment',2i6)
  400 format(//' print out from subroutine gensum'//' values of angular
     + momenta in *real* format'/,(14f5.1))
  401 format(/' racah w functions(6j)'/' arguments in *real* format'
     +, 18x,'value')
c
c
c
c   ***evaluates all terms in j6,j7,j8,j9,ldel,and jw which do not involve
c   ***a summation.the result is stored in recup and iastor
c
      if(ibug3.eq.1) then
        do  i=1,m
          xj1(i)=(j1(i)-one)/two
        enddo
        write(ipd ,400) (xj1(i),i=1,m)
        write(ipd ,306) nlsum
        write(ipd ,401)
      endif
c this puts the last "fixed" angular momentum to zero
      mm=m+1
      j1(mm)=1
c
c  ***test delta functions
c
      if(jdel.gt.0) then
          do  i=1,jdel
            i1=ldel(i,1)
            i2=ldel(i,2)
            if(i1.gt.mm)then
                j1(i1)=j1(i2)
            else if(i2.gt.mm)then
                j1(i2)=j1(i1)
            else if(j1(i1).ne.j1(i2)) then
                recup=zero
                if(ibug3.eq.1)write(ipd,308) i1,i2,j1(i1),j1(i2)
                return
            endif
          enddo
      endif !!end jdel.gt.0
c
      recup=one
c
c  ***multiply recup by all racah coefficients which do not involve a
c  ***summation
c
      if(jwc.ne.0)then
        if(ibug3.eq.1) write(ipd ,309)
c
        do  i=1,jwc
        if(inv6j(i).le.0)then
           do  j=1,6
           i1=jw(j,i)
           ist(j) = j1(i1) - 1
           enddo !! j
c
           call gracah(ist(1),ist(2),ist(3),ist(4),ist(5),ist(6),x1)
           if(ibug3.eq.1) write(ipd ,305) i,(xj1(jw(k,i)),k=1,6),x1
           recup = recup*x1
c
        endif
       enddo !! i
c
      endif !! end if(jwc.ne.0)
c
      sqr=1.0
c
      if(j6c.ne.0) then
        do  i=1,j6c
          i1=j6(i)
          sqr=sqr*j1(i1)
        enddo
      endif
c
      spr=1.0
c
      if(j9c.ne.0) then
        do  i=1,j9c
          i1=j9(i)
          spr=spr*j1(i1)
        enddo
      endif
c
      recup=recup*sqrt(sqr/spr)
      if(abs(recup).lt.epsil)then
          recup=zero
          return
       endif
      iastor = 0
c
      if(j7c.ne.0) then
        do   i=1,j7c
          i1=j7(i)
          iastor = iastor + j1(i1) -1
        enddo
      endif
c
      if(j8c.ne.0) then
        do  i=1,j8c
          i1=j8(i)
          iastor = iastor +2*(j1(i1)-1)
         enddo
      endif
c
      if(nlsum.le.0) then
        iastor=iastor/2
c
c  ***no summation involved. end of computation
c
        stor1=one
        stor=one
        if(mod(iastor,2).eq.1)recup=-recup
        if(ibug3.eq.1) write(ipd ,303) recup
        return
c
      endif
c
c
c  ***evaluation of the part involving summations.
c
c
      nfs=0
      jwr=0
      j6f=0
      j7f=0
      j8f=0
      j9f=0
      do 25 nps=1,nlsum
      if(ibug3.eq.1)write(ipd ,302) nps
c
c
c   *** loop on the disconnected summations
c
c
      ias=0
      nsum=nbj(nps)-nfs
      jwrd=nb6j(nps)-jwr
      j6cp=k6cp(nps)
      j7cp=k7cp(nps)
      j8cp=k8cp(nps)
      j9cp=k9cp(nps)
c
c
c     ***the range of values of each summation variable is
c     ***defined by establishing a matrix of the links between
c     ***variables.mat(i,j) contains:
c        i=j    number of possible values of i due to triangular
c               relations with non-variables,i.e. constants.
c       i.gt.j  number of links between i and j through constants
c       i.lt.j  value of the constant,if the above is 1.if not,
c               these values are srored in j12(l,i,j) where there
c               is room for mxcsvr such values (l.le.4)
c
c
      do 52 i=1,nsum
       do 152 j=1,nsum
         mat(i,j)=0
  152  continue
   52 continue
c
      do 66 i1=1,nsum
        i1t=i1+nfs
        i2=jsum6(i1t)
       do 65 i3=1,i2
         i=jsum5(i1t,i3)
         j=jsum4(i1t,i3)
c
c  ***the rows of the ipair arrays give limits of summation imposed
c
         do i5=1,2
           do i4=1,2
           ipair(i4,i5)=jword(jlimit(i4,i5,j),i)
           enddo
         enddo
        do 63 i4=1,2
              km=0
             do  i5=1,2
             if(ipair(i4,i5).gt.mp)km=km+1
             enddo
c
            jj1=ipair(i4,1)
            jj2=ipair(i4,2)
c
           if(km .lt.1)then
c  ***one variable linked to two constants.fix the diagonal mat(i,i)
c
          jt1=j1(jj1)-1
          jt2=j1(jj2)-1
          jmin=iabs(jt1-jt2)
          jmax=jt1+jt2
c
          if(mat(i1,i1) .gt. 1) then
c
c  ***if there are several couples of constants ,take the more
c  ***stringent combination
c
            jmin=max0(jmin,jsum(1,i1))
            jmax=min0(jmax,jsum(2,i1))
            if(jmax.ge.jmin)then
             jsum(1,i1)=jmin
             jsum(2,i1)=jmax
             mat(i1,i1)=(jmax-jmin)/2+1
             go to 63
            else
             recup=zero !! no possibility of value for this variable
             return
            endif
          else
          if(mat(i1,i1) .lt. 1) then
c
c  ***first time
c
            mat(i1,i1)=(jmax-jmin)/2+1
            jsum(1,i1)=jmin
            jsum(2,i1)=jmax
c
          endif
          endif
c
          elseif(km.eq.1)then
c
c
c  ***one variable linked to one constant and one variable  non diagonal
c  ***element
c
   67     jt1=min0(jj1,jj2)
          jt2=max0(jj1,jj2)-mp
          if(jt2.gt.i1)go to 63
          jt4=j1(jt1)-1
          k=mat(i1,jt2)
          if(k.ne.0)then
c
           do  ll=1,k
            if(jt4.eq.j12(ll,jt2,i1))go to 63
            enddo
c
       endif! k.ne.0
       k=k+1
      if(k.gt.mxcsvr)go to 63
      mat(i1,jt2)=k
      j12(k,jt2,i1)=jt4
c
      endif! endif km
c
   63 continue
   65 continue
   66 continue
c
c  ***reduce the diagonal elements by taking into account the non
c  ***diagonal elements,and keep the latter only if needed
c
  150 ichan=0
c
      do 74 i=1,nsum
        noel(i)=.true.
        i1=i-1
        jlow(1) = 1
        jhig(1) = i1
        jlow(2) = i+1
        jhig(2) = nsum
      do 720 ireduc = 1,2
        if(i1.eq.0)go to 170
       do 72  j=jlow(ireduc),jhig(ireduc)
       if(ireduc.eq.1)then
         ik1=i
         ik2=j
       else
         ik1 = j
         ik2 = i
       endif
         if(mat(ik1,ik2).eq.0 .or. mat(j,j) .eq. 0) go to 72
      jmin1=0
      jmax1=1000
      k=mat(ik1,ik2)
c
      do 203 l1=1,k
c
        l3=mat(j,j)
        jj1=jsum(1,j)
        jnd=j12(l1,ik2,ik1)
        jmin=1000
        jmax=0
        jmnp(l1)=0
        jmxp(l1)=1000
c
      do 204 l2=1,l3
c
        jmn=iabs(jnd-jj1)
        jmx=jnd+jj1
        jmin=min0(jmn,jmin)
        jmax=max0(jmx,jmax)
        jmnp(l1)=max0(jmn,jmnp(l1))
        jmxp(l1)=min0(jmx,jmxp(l1))
        jj1=jj1+2
c
  204 continue
c
      jmin1=max0(jmin1,jmin)
      jmax1=min0(jmax1,jmax)
c
  203 continue
c
      if(mat(i,i).eq.0) then
        jsum(1,i)=jmin1
        jsum(2,i)=jmax1
        mat(i,i)=(jmax1-jmin1)/2+1
        ichan=ichan+1
        go to 206
      endif
c
      if(jsum(1,i).lt.jmin1) then
        jsum(1,i)=jmin1
        ichan=ichan+1
      endif
c
      if(jsum(2,i).gt.jmax1) then
        jsum(2,i)=jmax1
        ichan=ichan+1
      endif
c
  206 k1=0
c
      do 207 l1=1,k
        if(jmnp(l1).le.jsum(1,i).and.jmxp(l1).ge.jsum(2,i))go to 207
        k1=k1+1
        j12(k1,ik2,ik1)=j12(l1,ik2,ik1)
  207 continue
c
      if(k1.ne.k) then
        mat(ik1,ik2)=k1
        ichan=ichan+1
      endif
c
      mat(ik2,ik1)=j12(1,ik2,ik1)
      if(ireduc.eq.1)noel(i)=.false.
   72  continue
c
c
  170 if(i.eq.nsum)go to 74
720   continue
c
c
   74 continue
c
      if(ichan.ne.0)go to 150
c
c
c  ***carry out the summations.
c
  220 do 230 i=1,nsum
        jsum3(i)=1
        ldiag(i)=.false.
        if(mat(i,i).eq.1)ldiag(i)=.true.
  230 continue
c
      do 231 i=1,jwrd
        jwtest(i)=1
  231 continue
c
      stor=zero
      stor1=one
      nolp=0
      ip=1
      ipold=ip
c   nolp is the index  of the summation variable
c  jsum2 is the value (format 2*J)
c  jsum3 is 1 if this variable is being changed from previous value
c
  240  nolp=nolp+1
c     find the range of values of summation variables, context dependent
          do  nsv=nolp,nsum
              jmin=jsum(1,nsv)
              jmax=jsum(2,nsv)
              if(.not.noel(nsv))then
c
              do   nj=1,nsv-1
                if(mat(nsv,nj) .eq. 1) then
                  jj1=mat(nj,nsv)
                  jj2=jsum2(nj)
                  jmin=max0(jmin,iabs(jj2-jj1))
                  jmax=min0(jmax,jj1+jj2)
                else if(mat(nsv,nj) .gt. 1) then
                  k=mat(nsv,nj)
                  jj2=jsum2(nj)
                  do  i=1,k
                  jj1=j12(i,nj,nsv)
                  jmin=max0(jmin,iabs(jj2-jj1))
                  jmax=min0(jmax,jj1+jj2)
                  enddo
                endif
               enddo
c
                      endif !! endif(noel)
                      jsum2(nsv)=jmin
                      maxlp(nsv)=jmax
                      if(ldiag(nsv))jsum3(nsv)=0
          enddo
c
          already=.false.
          ipold=ip
c loop on the last variable values
      do 260 jj=jmin,jmax,2
        jsum2(nsum)=jj
c
c  ***determine which racah coefficients need re-evaluating and
c  ***set jwtest appropriately
c
      do  j=ip,nsum
        if(jsum3(j).gt.0) then
           i2=jsum6(j)
           do i1=1,i2
           i3=jsum5(j,i1)
           jwtest(i3)=1
           enddo
        endif
      enddo
c
      do 98 j=1,jwrd
        if(jwtest(j).ne.0)then
        jwj=j+jwr
            do  i=1,6
            if(jword(i,jwj).le.mp) then
              i1=jword(i,jwj)
              ist(i) = j1(i1) - 1
            else
              i1=jword(i,jwj)-mp-nfs
              ist(i) = jsum2(i1)
            endif
           enddo
c
         call gracah(ist(1),ist(2),ist(3),ist(4),ist(5),ist(6),x1)
         wstor(j)=x1
          if(ibug3.eq.1) then
            do   i=1,6
            xj1(i)=ist(i)/two
            enddo
            write (ipd,305) jwj,(xj1(i), i=1,6),x1
          endif
       endif ! end jwtest(j).ne.0
   98 continue
c
c
c  ***form product of racah coefficients,(2j+1) factors and (-1)
c  ***factors in stor1
c
      do  i=1,jwrd
        stor1 = stor1*wstor(i)
      enddo
c
c  ***iastor contains the power of (-1)which is common to all terms
c
      ix2 = 0
      xij6cp=1.
      if(j6cp.ne.j6f) then
        jb=j6f+1
c [2*j+1]**1/2 factors
        do  i=jb,j6cp
        i1=j6p(i)-nfs
        xij6cp=xij6cp*(jsum2(i1)+1.)
        enddo
      endif
c 1/(2*j+1)**1/2 factors
      if(j9cp.ne.j9f) then
        jb=j9f+1
        do   i=jb,j9cp
        i1=j9p(i)-nfs
        xij6cp=xij6cp/(jsum2(i1)+1.)
         enddo
      endif
c
      stor1 = stor1*sqrt(xij6cp)
c
c simple phase factors: (-1)**j
      if(j7cp.ne.j7f) then
        jb=j7f+1
        do   i=jb,j7cp
          i1=j7p(i)-nfs
          ix2 = ix2 + jsum2(i1)
        enddo
      endif
c double phase factors: (-1)**2j
      if(j8cp.ne.j8f) then
        jb=j8f+1
        do   i=jb,j8cp
          i1=j8p(i)-nfs
          ix2 = ix2 + 2*(jsum2(i1))
         enddo
      endif
c
      if(mod(ix2,2).eq.1) then
        ias=-1
        ix2=ix2+1
      endif
c
      ix2 = ix2/2
c
c
c  ***add term into stor and reset stor1 to 1 ready for next term
c
      if (mod(ix2,2) .eq. 1) stor1 = -stor1
      stor = stor + stor1
      stor1=one
      nsum1 =nsum-1
      already=.true.
          if(nsum1.ne.0)then
              do  ik=1,nsum1
               jsum3(ik)=0
              enddo
              do  ik=1,jwrd
              jwtest(ik)=0
              enddo
           endif
  260 continue ! end sum on last variable values
c
c   control change of variable and
c  increment current value of jsum2(nolp)
c
       nolp=nsum
  250  nolp=nolp-1
       if(nolp.ne.0)then
        if(ldiag(nolp))go to 250
            jsum3(nolp)=1
            jsum2(nolp)=jsum2(nolp)+2
            if(jsum2(nolp).gt.maxlp(nolp))go to 250
                ip=nolp
                if(.not.already)ip=min(ipold,ip)
c....... proceed  to next variable
         go to 240
      endif
c end loop on summation variables
c
      recup=recup*stor
      if(ibug3.eq.1) write(ipd ,307) nps,stor,recup
      if(abs(recup).lt.epsil)then
        recup=zero
        return
      endif
      jwr=jwrd+jwr
      nfs=nsum+nfs
      j6f=j6cp
      j7f=j7cp
      j8f=j8cp
      j9f=j9cp
      iastor=iastor+ias
c
c  ***proceed to next sum
  25  continue
c
      iastor=iastor/2
      if(mod(iastor,2).ne.0)recup=-recup
      if(ibug3.eq.1) write(ipd ,304) recup
       return
c
      end
************************************************************************
      SUBROUTINE INTAB
*
*     THIS SUBROUTINE CALLED AT THE END OF DIAGRM,FIXES THE ARRAYS IH
*     AND IL-SO TO SPEAK HARDWARE AND LOGICAL ADDRESSES OF TRIADS IN
*     JDIAG.ALSO DETERMINES THE NUMBER OF FREE ENDS NFREE AND THEIR
*     LOCATION ITFREE.
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,KFL2=12,KFL2A=2*KFL2,KFL2B=4*KFL2,
     +          KFL2C=2*KFL2+2,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20)
*
      INTEGER ARR,TAB1
      COMMON/GRAPH/JDIAG(KFL2B,3),ARR(KFL2B,3),TAB1(KFL1,2),IL(KFL2B),
     + IH(KFL2B),NPOINT(KFL2C),NBNODE,IFIRST,ILAST,IPARTS,IPARTL,NPART,
     + ICROSS,NFREE,ITFREE(KFL2A),NFIN,NC
*
      COMMON/BUILD/IAL(KFL2B),IF1,IF2,NODE
*
*
      LOGICAL FREE
      COMMON/COUPLE/M,N,J1(KFL1),J2(KFL2,3),J3(KFL2,3),FREE(KFL1)
*
      DO 1 I=1,M
    1 IAL(I)=1
      DO 3 I=IFIRST,ILAST
      J=JDIAG(I,1)
      K=IAL(J)
      TAB1(J,K)=I
      IAL(J)=K+1
    3 CONTINUE
      IFR=IFIRST-1
      DO 4 I=IFIRST,ILAST
      IT=I-IFR
      IL(I)=IT
      IH(IT)=I
    4 CONTINUE
      J=JDIAG(IFIRST,3)
      K=IAL(J)
      IF(K-1)6,6,5
    5 TAB1(J,2)=TAB1(J,1)
    6 TAB1(J,1)=IFIRST
      IAL(J)=3
      J=JDIAG(ILAST,2)
      TAB1(J,2)=ILAST
      IAL(J)=3
      NFREE=0
      DO 7 I=IFIRST,ILAST
      J=JDIAG(I,1)
      IF(IAL(J).EQ.3)GO TO 7
      NFREE=NFREE+1
      ITT=ILAST+NFREE
      TAB1(J,2)=ITT
      IL(ITT)=NFREE*1000
      ITFREE(NFREE)=I
    7 CONTINUE
      RETURN
      END
************************************************************************
      SUBROUTINE KNJ(JD6C,JD7C,JD8C,JD9C,JDWC,JD6,JD7,JD8,JD9,KDW,JDDEL,
     +            LDDEL,DSUMVR,MDP,
     +            JD6P,JD7P,JD8P,JD9P,JDWORD,NDLSUM,NDBJ,NDB6J,
     +            KD6CP,KD7CP,KD8CP,KD9CP,JDSUM4,JDSUM5,JDSUM6,INVD6J)
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20,
     +          KFLS=12,KFLN=10,KFLV=40)
*
      LOGICAL DSUMVR
      DIMENSION JD6(KFL6),JD7(KFL7),JD8(KFL8),
     + JD9(KFL9),KDW(6,KFLW),LDDEL(KFLW,2),DSUMVR(KFL1)
*
      DIMENSION JD6P(KFLV),JD7P(KFLV),JD8P(KFLV),JD9P(KFLV),
     + JDWORD(6,KFLW),
     + NDBJ(KFLN),NDB6J(KFLN),KD6CP(KFLN),KD7CP(KFLN),KD8CP(KFLN),
     + KD9CP(KFLN),JDSUM6(KFLS),JDSUM4(KFLS,KFLW),JDSUM5(KFLS,KFLW),
     + INVD6J(KFLW)
*
      LOGICAL SUMVAR
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),KW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
*
      COMMON/SUMARG/J6P(KFLV),J7P(KFLV),J8P(KFLV),J9P(KFLV),
     + JWORD(6,KFLW),NLSUM,
     + NBJ(KFLN),NB6J(KFLN),K6CP(KFLN),K7CP(KFLN),K8CP(KFLN),K9CP(KFLN),
     + JSUM6(KFLS),JSUM4(KFLS,KFLW),JSUM5(KFLS,KFLW),INV6J(KFLW)
*
*
*     THIS ROUTINE STORES DATA FOR FUTURE CALLS TO GENSUM
*
      JD6C=J6C
      JD7C=J7C
      JD8C=J8C
      JD9C=J9C
      JDWC=JWC
      JDDEL=JDEL
      MDP=MP
      NDLSUM=NLSUM
      IF(J6C.NE.0) THEN
          DO 800 I=1,J6C
              JD6(I)=J6(I)
  800     CONTINUE
      ENDIF
      IF(J7C.NE.0) THEN
          DO 801 I=1,J7C
              JD7(I)=J7(I)
  801     CONTINUE
      ENDIF
      IF(J8C.NE.0) THEN
          DO 802 I=1,J8C
              JD8(I)=J8(I)
  802     CONTINUE
      ENDIF
      IF(J9C.NE.0) THEN
          DO 803 I=1,J9C
              JD9(I)=J9(I)
  803     CONTINUE
      ENDIF
      IF(JWC.NE.0) THEN
          DO 804 I=1,6
              DO 804 J=1,JWC
                  KDW(I,J)=KW(I,J)
  804     CONTINUE
          DO 805 I=1,JWC
              INVD6J(I)=INV6J(I)
  805     CONTINUE
      ENDIF
      IF(JDEL.NE.0) THEN
          DO 806 I=1,2
              DO 806 J=1,JDEL
                  LDDEL(J,I)=LDEL(J,I)
  806     CONTINUE
      ENDIF
      IF(MP.NE.0) THEN
          DO 807 I=1,MP
              DSUMVR(I)=SUMVAR(I)
  807     CONTINUE
      ENDIF
      IF(NLSUM.NE.0) THEN
          DO 808 I=1,NLSUM
              NDBJ(I)=NBJ(I)
              NDB6J(I)=NB6J(I)
              KD6CP(I)=K6CP(I)
              KD7CP(I)=K7CP(I)
              KD8CP(I)=K8CP(I)
              KD9CP(I)=K9CP(I)
  808     CONTINUE
      ENDIF
      DO 809 I=1,KFLV
          JD6P(I)=J6P(I)
          JD7P(I)=J7P(I)
          JD8P(I)=J8P(I)
          JD9P(I)=J9P(I)
  809 CONTINUE
      DO 810 I=1,KFLS
          JDSUM6(I)=JSUM6(I)
          DO 810 J=1,KFLW
              JDSUM4(I,J)=JSUM4(I,J)
              JDSUM5(I,J)=JSUM5(I,J)
  810 CONTINUE
      DO 811 I=1,6
          DO 811 J=1,KFLW
              JDWORD(I,J)=JWORD(I,J)
  811 CONTINUE
*
*
      RETURN
      END
************************************************************************
      subroutine lolpop(fail)
************************************************************************
*
*  ***reduces a loop with one line and one node in the flat graph.
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,KFL2=12,KFL2A=2*KFL2,KFL2B=4*KFL2,
     +          KFL2C=2*KFL2+2,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20)
*
      LOGICAL FAIL
*
      INTEGER ARR,TAB1
      COMMON/GRAPH/JDIAG(KFL2B,3),ARR(KFL2B,3),TAB1(KFL1,2),IL(KFL2B),
     + IH(KFL2B),NPOINT(KFL2C),NBNODE,IFIRST,ILAST,IPARTS,IPARTL,NPART,
     + ICROSS,NFREE,ITFREE(KFL2A),NFIN,NC
*
      LOGICAL SUMVAR
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),KW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
      LOGICAL FREE
      COMMON/COUPLE/M,N,J1(KFL1),J2(KFL2,3),J3(KFL2,3),FREE(KFL1)
*
      dimension kp(3), ks(3)
*
      character*6 name,namsub
      common/nam/namsub
*
      data name/'lolpop'/
      data kp/2,3,1/
      data ks/0,1,-1/
*
*
*
      namsub=name
      i1=npoint(1)
      k3=2
      if(i1.eq.ilast)k3=3
      l=jdiag(i1,k3)
      call delta(l,mp,fail)
      if(fail)return
      k=kp(k3)
      if(arr(i1,k).lt.0)call phase2(jdiag(i1,k))
      k1=ks(k3)
      il1=il(i1)+k1
      i2=ih(il1)
      l1=jdiag(i2,1)
      call delta(l1,jdiag(i2,k3),fail)
      if(fail)return
      if(arr(i2,k3).eq.k1)call phase2(l1)
      il2=il(i2)+k1
      i3=ih(il2)
      k2=k3+k1
      jdiag(i3,k2)=l1
      arr(i3,k2)=arr(i2,1)
      j9c=j9c+1
      j9(j9c)=l1
      j6c=j6c+1
      j6(j6c)=jdiag(i1,1)
      if(k3.eq.3)return
*
      do 1 i=3,nbnode
        it=ih(i)
        ilp=i-2
        il(it)=ilp
        ih(ilp)=it
    1 continue
*
      return
      end
************************************************************************
      SUBROUTINE NEIBOR(LC,L1,L2)
*
*    GIVES THE POSITIONS OF THE OTHER TWO ARGUMENTS IN THE TRIAD.
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      IF(LC-2)1,2,3
    1 L1=2
      L2=3
      GO TO 4
    2 L1=3
      L2=1
      GO TO 4
    3 L1=1
      L2=2
    4 RETURN
      END
************************************************************************

      SUBROUTINE NJGRAF(RECUP,FAIL)

************************************************************************
*
*     WRITTEN BY:
*                     A. BAR-SHALOM AND M. KLAPISCH
*                     RACAH INSTITUTE OF PHYSICS
*                     HEBREW UNIVERSITY
*                     91904 JERUSALEM
*                     ISRAEL.
*
************************************************************************
*
*-----THIS IS THE MAIN PROGRAM.IT HANDLES  ALL THE ANALYSIS OF THE
*     RECOUPLING COEFFICIENT WITHOUT REFERRING EXPLICITLY TO THE VALUES
*     OF ANGULAR MOMENTA WHICH ARE IN J1(J),EXCEPT FOR ZERO IN CASE FREE
*     =.FALSE. .LIKE NJSYM IT PREPARES ARRAYS OF ARGUMENTS FOR PHASE
*     FACTORS,(2*J+1) FACTORS AND 6J COEFFICIENTS TO BE COMPUTED IN
*     GENSUM,WHICH CAN ALSO BE CALLED SEPARATELY WHEN ONLY THE NUMERICAL
*     VALUES OF ANGULAR MOMENTA CHANGE.THESE VARIABLE ANGULAR MOMENTA
*     SHOULD BE DECLARED FREE(J)=.TRUE.,SO THAT THE FORMULA PREPARED FOR
*     GENSUM SHOULD BE CORRECT WHEN  J1 IS NOT ZERO.
*     FAIL WILL BE TRUE WHEN THE RECOUPLING COEFFICIENT IS ZERO BECAUSE
*     OF UNSATISFIED DELTA OR OTHER SIMILAR CAUSES.
*
************************************************************************
*
*     THIS VERSION HAS BEEN MODIFIED TO HOLD THE ARRAY DIMENSIONS IN
*     PARAMETER STATEMENTS. THE DIMENSIONS ARE LABELLED:
*
*     KFL1 - DIMENSION OF THE J1 AND FREE ARRAYS IN /COUPLE/, AND THE
*            FIRST DIMENSION OF THE LINE AND LCOL ARRAYS IN /TREE/.
*            ALSO THE DIMENSION OF THE SUMVAR ARRAY IN /ARGU/, AND
*            OF THE INVER ARRAY IN ROUTINE SPRATE. IT IS TESTED FOR M
*            ON ENTRY TO NJGRAF, AND FOR MP IN ROUTINE SPRATE.
*
*     KFL2 - DIMENSION OF THE J2 AND J3 ARRAYS IN /COUPLE/. THE
*            DIMENSIONS OF THESE ARRAYS ARE CHECKED ON ENTRY TO NJGRAF.
*
*     KFL2A = (2*KFL2) - DIMENSION OF THE J23, ARROW AND TABS ARRAYS IN
*                      /TREE/.
*     KFL2B = (4*KFL2) - DIMENSION OF THE JDIAG,ARR,IL AND IH ARRAYS
*                            IN /GRAPH/, AND OF THE IAL ARRAY IN /BUILD/
*     KFL2C = (2*KFL2+2) - DIMENSION OF THE NPOINT ARRAY IN /GRAPH/.
*
*     KFL6 - DIMENSION OF THE J6 ARRAY IN /ARGU/. TESTED IN SPRATE.
*
*     KFL7 - DIMENSION OF THE J7 ARRAY IN /ARGU/. TESTED IN SPRATE.
*
*     KFL8 - DIMENSION OF THE J8 ARRAY IN /ARGU/. TESTED IN SPRATE.
*
*     KFL9 - DIMENSION OF THE J9 ARRAY IN /ARGU/. TESTED IN SPRATE.
*
*     KFLW - DIMENSION OF THE JW(OR KW) AND LDEL ARRAYS IN /ARGU/, AND
*            OF THE JWORD AND INV6J ARRAYS IN /SUMARG/. ALSO THE SECOND
*            DIMENSION OF THE JSUM4 AND JSUM5 ARRAYS IN /SUMARG/.
*            IN ADDITION IT GIVES THE DIMENSIONS OF A NUMBER OF
*            TEMPORARY WORKING ARRAYS IN ROUTINES SPRATE AND GENSUM.
*            KFLW IS TESTED IN SPRATE.
*
*     KFLS - DIMENSION OF THE JSUM6 ARRAY AND FIRST DIMENSION OF THE
*            JSUM4 AND JSUM5 ARRAYS IN /SUMARG/. ALSO GIVES THE
*            DIMENSIONS OF SOME TEMPORARY WORKING ARRAYS IN SPRATE
*            AND GENSUM. KFLS IS THE MAXIMUM NUMBER OF SUMMATION
*            VARIABLES IN A PARTICULAR SUM, AND IS TESTED IN SPRATE.
*
*     KFLN - DIMENSION OF THE NBJ, NB6J, K6CP, K7CP, K8CP AND K9CP
*            ARRAYS IN /SUMARG/. KFLN IS THE MAXIMUM NUMBER OF
*            SUMS ALLOWED, AND IS TESTED IN ROUTINE SPRATE.
*
*     KFLV - DIMENSION OF THE J6P, J7P, J8P AND J9P ARRAYS IN
*            /SUMARG/, AND OF THE JNS ARRAY IN ROUTINE VAR.
*            KFLV IS TESTED IN VAR.
*
*     KFLZ - DIMENSION OF THE JZERO ARRAY IN /ZER/. KFLZ IS TESTED IN
*            ROUTINE ZERO.
*
*******************************************************************
*
*     OTHER CHANGES:
*
*  (1) THIS VERSION HAS BEEN ALTERED TO USE THE RACAH
*     COEFFICIENT ROUTINE WRITTEN BY STAN SCOTT AND ALAN HIBBERT
*     (CPC 28 189-200 1982). SUBROUTINE FACTT MUST BE CALLED TO
*     SET UP THE ARRAY GAM WHICH HOLDS LN(N!) BEFORE NJGRAF IS
*     CALLED.
*  (2) THE SUBROUTINES ORDER AND SETDIM HAVE BEEN RENAMED AS
*     ORDTRI AND SETDM, TO AVOID CONFLICT IN THE QUB CODES.
*  (3) COMMON BLOCKS /ARGU/ AND /SUMARG/ HAVE BEEN REMOVED FROM
*     SUBROUTINE GENSUM, AND THEIR CONTENTS ARE TRANSFERRED VIA
*     THE ARGUMENT LIST INSTEAD. /SUMARG/ HAS THUS BEEN INCLUDED
*     IN SUBROUTINE NJGRAF TO ACHIEVE THIS.
*  (4) THE TIMING ROUTINES HAVE BEEN REMOVED FROM THE NJGRAF
*     PACKAGE AS THEY WERE CDC SPECIFIC AND SEEMED UNNECESSARY.
*
*****************************************************************
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,KFL2=12,KFL2A=2*KFL2,KFL2B=4*KFL2,
     +          KFL2C=2*KFL2+2,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20,
     +          KFLS=12,KFLN=10,KFLV=40)
*
      LOGICAL FAIL,FIND
*
      INTEGER ARROW
      LOGICAL TABS
      COMMON/TREE/J23(KFL2A,3),ARROW(KFL2A,3),LINE(KFL1,2),LCOL(KFL1,2),
     +  TABS(KFL2A),NBTR
      LOGICAL CUT
      COMMON/CUTDIG/CUT
      INTEGER ARR,TAB1
      COMMON/GRAPH/JDIAG(KFL2B,3),ARR(KFL2B,3),TAB1(KFL1,2),IL(KFL2B),
     + IH(KFL2B),NPOINT(KFL2C),NBNODE,IFIRST,ILAST,IPARTS,IPARTL,NPART,
     + ICROSS,NFREE,ITFREE(KFL2A),NFIN,NC
*
      LOGICAL FREE
      COMMON/COUPLE/M,N,J1(KFL1),J2(KFL2,3),J3(KFL2,3),FREE(KFL1)
*
      LOGICAL SUMVAR
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),KW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
*
      COMMON/SUMARG/J6P(KFLV),J7P(KFLV),J8P(KFLV),J9P(KFLV),
     + JWORD(6,KFLW),NLSUM,
     + NBJ(KFLN),NB6J(KFLN),K6CP(KFLN),K7CP(KFLN),K8CP(KFLN),K9CP(KFLN),
     + JSUM6(KFLS),JSUM4(KFLS,KFLW),JSUM5(KFLS,KFLW),INV6J(KFLW)
*
      CHARACTER*6 NAME,NAMSUB
      COMMON/NAM/NAMSUB
      DATA NAME/'NJGRAF'/
*
*     TEST THE DIMENSION OF THE J1 ARRAY
*
      IF(M+1.GT.KFL1) THEN
          PRINT 100,M+1,KFL1
  100     FORMAT(1X,'DIMENSION ERROR IN NJGRAF. M+1=',I5,' KFL1=',I5)
          STOP
      ENDIF
*
*     TEST THE DIMENSIONS OF THE J2 AND J3 ARRAYS
*
      IF(N-1.GT.KFL2) THEN
          PRINT 101,N-1,KFL2
  101     FORMAT(1X,'DIMENSION ERROR IN NJGRAF. N-1=',I5,' KFL2=',I5)
          STOP
      ENDIF
*
*    1.0
*    BUILDING UP OF THE UNSTRUCTURED GRAPH
*
      do 501 i = N+1,kfl2
      do 501 j = 1,3
      j2(i,j) = 0
      j3(i,j) = 0
  501 continue
      FAIL=.FALSE.
      J6C=0
      J7C=0
      J8C=0
      J9C=0
      JWC=0
      JDEL=0
      CALL SETDM
      NFIN=0
      CUT=.FALSE.
      CALL SETTAB(FAIL)
      M=M+1
      IF(FAIL)GO TO 7
      M=M-1
      JF=0
      JF1=0
*
*   1.1
*      LOCATING AND HANDLING OF ZEROS
*
      CALL ZERO(JF1,JF,FAIL)
      IF(FAIL)GO TO 7
      MP=M
      IF(NBTR.EQ.0)GO TO 6
      JUMP=1
c
c   ***locating and handling order 1 loops - ususally consequences of
c      zeros or cut in preceding diagram
c
    1 call chklp1(fail)
      if (fail) go to 7

*
*    2.
*      BUILDING OF A FLAT DIAGRAM OUT OF THE UNSTRUCTURED GRAPH.
*      THERE MAY BE SEVERAL FLAT DIAGRAMS OUT OF THE ORIGINAL
*      GRAPH,IN CASE OF POSSIBLE CUTS.THEN THE FLAT DIAGRAMS
*      WILL HAVE FREE ENDS.
*
      CALL DIAGRM(JUMP)
      NFIN=MAX0(0,NFREE-2)
      IF(NFIN.EQ.0)GO TO 2
      JUMP=3
*
*    2.1
*     HANDLING OF FREE ENDS IF A CUT WAS FOUND
*
      CALL CUTNL(FAIL)
      IF(FAIL)GO TO 7
      GO TO 5
    2 JUMP=2
      IF(NFREE-1)5,3,4
    3 CALL CUT1L(FAIL)
      IF(FAIL)GO TO 7
      GO TO 5
    4 CALL CUT2L(FAIL)
      IF(FAIL)GO TO 7
    5 NBTR=NBTR+NFIN
      IF(NBTR.NE.0)CUT=.TRUE.
*
*    3.
*    ANALYSIS OF THE FLAT DIAGRAM.
*     CLOSED CIRCUITS OF INCREASING ORDER NC ARE SEARCHED,ANALYSED,AND
*     TAKEN OUT OF THE FLAT DIAGRAM,THUS REDUCING THE NUMBER OF NODES,
*     NBNODE.
*
      NC=0
   10 NC=NC+1
      CALL SEARCH(FIND)
      IF(.NOT.FIND)GO TO 10
      NCP=NC-2
      JPOL=0
      IF(M.EQ.MP.AND.NC.GT.3)CALL SETDM
      IF(IPARTL.GT.2)CALL POLYGN(JPOL)
      GO TO (11,12,13,14),NC
   11 CALL LOLPOP(FAIL)
      IF(FAIL)GO TO 7
      GO TO 15
   12 CALL BUBBLE(JPOL,FAIL)
      IF(FAIL)GO TO 7
      GO TO 15
   13 CALL TRIANG(FAIL)
      IF(FAIL)GO TO 7
      GO TO 15
   14 CALL SQUARE
   15 NBNODE=NBNODE-2
      IF(NBNODE.EQ.0)GO TO 9
      IFIRST=IH(1)
      ILAST=IH(NBNODE)
*      PRINTJ IS AN ALL PURPOSE PRINTING SUBROUTINE CALLED FROM MANY PLA
      CALL PRINTJ(NAMSUB,8)
      IF(NBNODE.EQ.NFIN)GO TO 9
*     .. correction for avoiding nc=0..4/13/94
      NC=max(NCP,0)
*
*    PROCEED TO OTHER CIRCUITS OF ORDER NC-1
*
      GO TO 10
    9 IF(NBTR.EQ.0)GO TO 6
      IF(JUMP.EQ.3)CALL ORDTRI
*
*    AT THIS STAGE,THE FLAT DIAGRAM HAS BEEN REDUCED TO NODES
*    INVOLVING FREE ENDS.PROCEED TO BUILD OTHER FLAT DIAGRAMS
*    IF NECESSARY.
*
      GO TO 1
*
*   ALL PARTS OF THE ORIGINAL GRAPH HAVE BEEN REDUCED.
*
    7 RECUP=0.
      M=M-1
      RETURN
    6 CALL PRINTJ(NAME,0)
*
*    4.PREPARATION OF THE RESULTS,AND SEPARATION IN SEVERAL SUMS
*      IF CUTS HAVE BEEN DETECTED,ALSO IN THE FLAT DIAGRAM ITSELF
*
      CALL SPRATE(M)
      M=M-1
*
*    5. GENSUM COMPUTES THE NUMERICAL VALUE OF THE RECOUPLING COEFFICIEN
*
      CALL GENSUM(J6C,J7C,J8C,J9C,JWC,J6,J7,J8,J9,KW,JDEL,
     +            LDEL,SUMVAR,MP,
     +            J6P,J7P,J8P,J9P,JWORD,NLSUM,NBJ,NB6J,
     +            K6CP,K7CP,K8CP,K9CP,JSUM4,JSUM5,JSUM6,INV6J,
     +            RECUP)
      RETURN
      END
************************************************************************
      subroutine ordtri
************************************************************************
*
*  ***this subroutine orders the triads which were left with free ends
*  ***as consequence of cutting,so that the new graph will start there.
*
* ... modified by Bar Shalom, July, 1989; November 1989
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,KFL2=12,KFL2A=2*KFL2,KFL2B=4*KFL2,
     +          KFL2C=2*KFL2+2,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20)
*
      LOGICAL SUMVAR
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),KW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
      INTEGER ARROW
      LOGICAL TABS
      COMMON/TREE/J23(KFL2A,3),ARROW(KFL2A,3),LINE(KFL1,2),LCOL(KFL1,2),
     +  TABS(KFL2A),NBTR
*
      INTEGER ARR,TAB1
      COMMON/GRAPH/JDIAG(KFL2B,3),ARR(KFL2B,3),TAB1(KFL1,2),IL(KFL2B),
     + IH(KFL2B),NPOINT(KFL2C),NBNODE,IFIRST,ILAST,IPARTS,IPARTL,NPART,
     + ICROSS,NFREE,ITFREE(KFL2A),NFIN,NC
*
      COMMON/BUILD/IAL(KFL2B),IF1,IF2,NODE
*
*
      COMMON/KEEP/JKP(2,3),JARR(2,3),IT2,IT3,IT5
*
      character*6 name
      data name/'ordtri'/
*
*
      do 10 i=1,mp
        ial(i)=0
   10 continue
*
      if(nfin.ne.0)then
        nbt1=nbtr-1
        nbt=nbt1+nfin
        nbtt=nbt+1
        nb=0
        go to 31
      endif
*
      nf=nbtr-itfree(1)
*
      if(it5.eq.0) then
        nbt1=nbtr-1
        n0=0
        nft=nfree
        isw=2
        go to 100
      endif
*
      nft=it5-it2
      nm=nft+nbtr+1
      nbt1=nbtr
*
      do 21 j=1,3
        jdiag(nbtr,j)=jkp(1,j)
        arr(nbtr,j)=jarr(1,j)
*  June 25  1989  (avi)
*        jdiag(nm,j)=jkp(2,j)
*        arr(nm,j)=jarr(2,j)
   21 continue
*
      jt=jdiag(nm,1)
      n0=0
      isw= 1
      go to 100
*
   22 n0=nft
* June 25 1989 (avi)
      do 211 j=1,3
        jdiag(nm,j)=jkp(2,j)
        arr(nm,j)=jarr(2,j)
 211  continue
* June 25 1989 (avi)
*      nbt1=nbt1+n0
      nbt1=nbt1+1
      nft=it3-it5
      isw= 3
      go to 100
*
* June 25 1989 (avi)
*   24 nft=nft+1
   24 nbt1=k-nft
*
   23 node=nbt1+nft
      call change(node,2)
      go to 40
*
*
   31 do 35 i=1,nbnode
        i1=ih(i)
        if(il(i1).gt.ilast)go to 35
        i2=nbt1+i
        if(i1.gt.nbtt)go to 33
        if(i1.eq.i2)go to 32
        if(il(i2).le.nbnode)go to 35
*
   33 do 34 j=1,3
        jdiag(i2,j)=jdiag(i1,j)
        arr(i2,j)=arr(i1,j)
   34 continue
*
      il(i1)=ilast+i
   32 nb=nb+1
      il(i2)=0
*
   35 continue
*
      if(nb.ne.nfin)go to 31
      node=nbt
   40 if1=jdiag(nbtr,1)
      if2=jdiag(nbtr,3)
*
      do 51 i=nbtr,node
       do 50 k=1,3
         j=jdiag(i,k)
         ial(j)=ial(j)+1
   50  continue
   51 continue
*
      ilast=node
      call printj(name,8)
*
      return
*
  100 if(nf.le.0)then
        nfr=n0
        i1=1
      else
        nfr=nft+1
        i1=-1
      endif
*
      do 4 i=1,nft
        ik=nfr+i1*i
        it=itfree(ik)
        k=nbt1+ik
*
      do 3 j=1,3
        jdiag(k,j)=jdiag(it,j)
        arr(k,j)=arr(it,j)
    3 continue
*
    4 continue
*
      go to (22,23,24),isw
      end
************************************************************************
      SUBROUTINE OTHERJ(LIN,J,LO,LCO,K)
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
*     GIVES THE OTHER TRIAD WHERE A GIVEN J OCCURS AND ITS POSITION.
*
*
      PARAMETER(KFL1=60,KFL2=12,KFL2A=2*KFL2)
*
      LOGICAL TABS
      INTEGER ARROW
      COMMON/TREE/J23(KFL2A,3),ARROW(KFL2A,3),LINE(KFL1,2),LCOL(KFL1,2),
     +  TABS(KFL2A),NBTR
*
      LO=LINE(J,1)
      IF(LO.EQ.LIN.OR.TABS(LO))GO TO 1
      K=2
      LCO=LCOL(J,1)
      GO TO 2
    1 K=1
      LO=LINE(J,2)
      LCO=LCOL(J,2)
    2 RETURN
      END
************************************************************************
      SUBROUTINE PHASE(L,JM,NDIM)
*
*    PHASE FACTOR ARISING FROM NON-CYCLIC PERMUTATION OF ARGUMENTS IN
*    TRIAD L.JM MAY BE EITHER J23 OR JDIAG.
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20)
*
      DIMENSION JM(NDIM,3)
      LOGICAL SUMVAR
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),KW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
*
      J7(J7C+1)=JM(L,1)
      J7(J7C+2)=JM(L,2)
      J7C=J7C+3
      J7(J7C)=JM(L,3)
      RETURN
      END
************************************************************************
      SUBROUTINE PHASE2(J)
*
*     ADDS A PHASE FACTOR (-1)**2J
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20)
*
      LOGICAL SUMVAR
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),KW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
*
      J8C=J8C+1
      J8(J8C)=J
      RETURN
      END
************************************************************************
      SUBROUTINE POLYGN(JPOL)
*
*    THIS ROUTINE REDUCES A CIRCUIT OF ARBITRARY ORDER NC.IT EXCHANGES
*    NODES ON THE FLAT DIAGRAM UNTIL THE DISTANCE ON THE AXIS BETWEEN
*    NODES EQUEALS ONE.EACH EXCHANGE INTRODUCES A SUMMATION VARIABLE
*    AND A 6J SYMBOL.THE CIRCUIT HAS A MAXIMUM OF NPART=2 DISCONNECTED
*    PARTS ON THE AXIS.
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,KFL2=12,KFL2A=2*KFL2,KFL2B=4*KFL2,
     +          KFL2C=2*KFL2+2,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20)
*
      INTEGER ARR,TAB1
      COMMON/GRAPH/JDIAG(KFL2B,3),ARR(KFL2B,3),TAB1(KFL1,2),IL(KFL2B),
     + IH(KFL2B),NPOINT(KFL2C),NBNODE,IFIRST,ILAST,IPARTS,IPARTL,NPART,
     + ICROSS,NFREE,ITFREE(KFL2A),NFIN,NC
      LOGICAL SUMVAR
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),KW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
*
      CHARACTER*6 NAME
      DATA NAME/'POLYGN'/
*
      NC1=NC+1
      NC2=NC
      NBC=IPARTL-2
   10 DO 8 I=1,NBC
      IT2=NPOINT(NC1-I)
      IT1=NPOINT(NC2-I)
      JB=JDIAG(IT1,1)
      JC=JDIAG(IT2,1)
      JDIAG(IT1,1)=JC
      JDIAG(IT2,1)=JB
      JAR=ARR(IT1,1)
      ARR(IT1,1)=ARR(IT2,1)
      ARR(IT2,1)=JAR
      JE=JDIAG(IT1,2)
      MP=MP+1
      SUMVAR(MP)=.TRUE.
      JDIAG(IT1,2)=MP
      JDIAG(IT2,3)=MP
      IF(TAB1(JB,1)-IT1)2,1,2
    1 TAB1(JB,1)=IT2
      GO TO 3
    2 TAB1(JB,2)=IT2
    3 IF(TAB1(JC,1)-IT2)5,4,5
    4 TAB1(JC,1)=IT1
      GO TO 6
    5 TAB1(JC,2)=IT1
    6 IF(ARR(IT1,2).GT.0)GO TO 7
      CALL PHASE2(JE)
      ARR(IT1,2)=1
      ARR(IT2,3)=-1
    7 JWC=JWC+1
      KW(1,JWC)=JB
      KW(2,JWC)=MP
      KW(3,JWC)=JE
      KW(4,JWC)=JC
      KW(5,JWC)=JDIAG(IT2,2)
      KW(6,JWC)=JDIAG(IT1,3)
      J6(J6C+1)=MP
      J6C=J6C+2
      J6(J6C)=MP
    8 CONTINUE
      NC=NC-NBC
      IF(NC.LE.4)GO TO 11
      NBC=IPARTS-2
      NC1=IPARTS+1
      NC2=IPARTS
      GO TO 10
   11 IF(NPART.EQ.1)GO TO 12
      NPOINT(3)=NPOINT(NC1)
      NPOINT(4)=NPOINT(NC1+1)
   12 IF(NC.EQ.2)JPOL=1
      CALL PRINTJ(NAME,10)
      RETURN
      END
************************************************************************
      SUBROUTINE PRINTJ(NAMES,JP)
*
*    THIS SUBROUTINE PRINTS INTERMEDIATE RESULTS IN STANDARD FORM FROM
*    WHEREVER IT IS CALLED.
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,KFL2=12,KFL2A=2*KFL2,KFL2B=4*KFL2,
     +          KFL2C=2*KFL2+2,
     +          KFLZ=30,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20)
*
      CHARACTER IM,IP,IS(3)
      CHARACTER*4 I6,I7,I8,I9,IJ1
      CHARACTER*6 NAMES,NSETTB
      CHARACTER*8 IBLANK,IFREE,IFR
*
      DIMENSION IX(6),JTAB(KFL1,3)
      INTEGER ARROW
      LOGICAL TABS
      COMMON/TREE/J23(KFL2A,3),ARROW(KFL2A,3),LINE(KFL1,2),LCOL(KFL1,2),
     +  TABS(KFL2A),NBTR
      LOGICAL SUMVAR
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),KW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
      INTEGER ARR,TAB1
      LOGICAL FREE
      COMMON/COUPLE/M,N,J1(KFL1),J2(KFL2,3),J3(KFL2,3),FREE(KFL1)
      COMMON/GRAPH/JDIAG(KFL2B,3),ARR(KFL2B,3),TAB1(KFL1,2),IL(KFL2B),
     + IH(KFL2B),NPOINT(KFL2C),NBNODE,IFIRST,ILAST,IPARTS,IPARTL,NPART,
     + ICROSS,NFREE,ITFREE(KFL2A),NFIN,NC
      COMMON/CONST/I6C,I7C,I8C,I9C,IDEL,IWC
      COMMON/ZER/NZERO,JZERO(KFLZ)
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      EQUIVALENCE(I6C,IX(1))
      DATA IBLANK,IFREE,IP,IM/'        ','FREE END','+','-'/
      DATA NSETTB/'SETTAB'/
      DATA I6,I7,I8,I9,IJ1/'I6= ','I7= ','I8= ','I9= ','J1= '/
      DATA I6C,I7C,I8C,I9C,IDEL,IWC/6*1/
 1000 FORMAT(/10X,'NBNODE=',I3,10X,'NBTR=',I3,10X,'NFIN=',I3,
     + /10X,'IFIRST=',I3,10X,'ILAST=',I3,9X,'NFREE=',I3)
 1001 FORMAT(//7X,'IL',3X,'IH',14X,'JDIAG'//)
 1002 FORMAT(28X,3(A1,2X))
 1003 FORMAT(7X,I2,3X,I2,2X,A8,2X,3I3/)
 1004 FORMAT(/5X,'TAB1'/)
 1005 FORMAT(4(I3,1H),2X,I3,I5,5X))
 1006 FORMAT(/2X,'SUMVAR=',15(I3,L1))
 1010 FORMAT(//10X,'J23',10X,'NBTR1=',I3//)
 1012 FORMAT(18X,3(A1,2X))
 1013 FORMAT(I9,I5,2X,3I3/)
 1014 FORMAT(/3X,'J  L1 K1  L2 K2')
 1015 FORMAT(4(I4,1H),I3,I3,I4,I3))
 1020 FORMAT(/3X,A4,3X,3(20I3/))
 1021 FORMAT(/3X,'DELTA=',7(I5,I3))
 1022 FORMAT(/3X,'KW(ARG. OF 6J)',6I3)
 1030 FORMAT(//2X,'NC=',I2,4X,'NPART=',I2,4X,'IPARTL=',I2,4X,
     +  'IPARTS=',I2,4X,'ICROSS=',I2,4X,/2X,'NPOINT=',20I3)
 1040 FORMAT(//2X,'NZERO=',I2,5X,12(I4,1H),I3))
 1050 FORMAT(///3X,'PRINT OUT AFTER CALLING SUBROUTINE ',A7)
*
      IF(IBUG3.NE.1)RETURN
      PRINT 1050,NAMES
      JUMP=JP
      IF(JUMP.NE.0)GO TO 10
      DO 9 I=1,6
    9 IX(I)=1
      PRINT 1020,IJ1,(J1(I),I=1,M)
   10 IF(JUMP.LT.8)GO TO 20
      PRINT 1000,NBNODE,NBTR,NFIN,IFIRST,ILAST,NFREE
      JUMP=JUMP-8
      PRINT 1001
      K=0
      DO 1 I=1,NBNODE
      IT=IH(I)
      IFR=IBLANK
      JT=JDIAG(IT,1)
      IF(TAB1(JT,2).EQ.IT.AND.JT.NE.JDIAG(IFIRST,3))GO TO 6
      K=K+1
      JTAB(K,1)=JT
      JTAB(K,2)=TAB1(JT,1)
      JTAB(K,3)=TAB1(JT,2)
    6 IF(TAB1(JT,2).GT.ILAST)IFR=IFREE
      DO 2 J=1,3
      IS(J)=IP
      IF(ARR(IT,J).LT.1)IS(J)=IM
    2 CONTINUE
      PRINT 1002,(IS(J),J=1,3)
      PRINT 1003,IL(IT),IT,IFR,(JDIAG(IT,J),J=1,3)
    1 CONTINUE
      PRINT 1004
      NTIME=0
      JT=JDIAG(IFIRST,3)
      IF(JT.EQ.JDIAG(ILAST,2))GO TO 4
      IF(TAB1(JT,2).GE.1000)GO TO 4
      GO TO 5
    4 K=K+1
      JTAB(K,1)=JT
      JTAB(K,2)=TAB1(JT,1)
      JTAB(K,3)=TAB1(JT,2)
    5 NTIME=NTIME+1
      IF(NTIME.EQ.2)GO TO 3
      JT=JDIAG(ILAST,2)
      IF(TAB1(JT,2).EQ.1000)GO TO 4
    3 PRINT 1005,((JTAB(I,J),J=1,3),I=1,K)
      PRINT 1006,(I,SUMVAR(I),I=1,MP)
   20 IF(JUMP.LT.4)GO TO 30
      JUMP=JUMP-4
      NBTR1=2*N-2
      PRINT 1010,NBTR1
      K=0
      DO 11 I=1,NBTR1
      IF(TABS(I))GO TO 11
      K=K+1
      DO 12 J=1,3
      IS(J)=IP
      IF(ARROW(I,J).LT.1)IS(J)=IM
   12 CONTINUE
      PRINT 1012,(IS(J),J=1,3)
      PRINT 1013,K,I,(J23(I,J),J=1,3)
   11 CONTINUE
      PRINT 1014
      MM=M
      IF(NAMES.NE.NSETTB) MM=M-1
      PRINT 1015,(I,(LINE(I,J),LCOL(I,J),J=1,2),I=1,MM)
   30 IF(JUMP.LT.2)GO TO 40
      JUMP=JUMP-2
      PRINT 1030,NC,NPART,IPARTL,IPARTS,ICROSS,(NPOINT(I),I=1,NC)
   40 IF(JUMP.LT.1)GO TO 50
      PRINT 1040,NZERO,(I,JZERO(I),I=1,NZERO)
   50 IF(J6C.GE.I6C)PRINT 1020,I6,(J6(I),I=I6C,J6C)
      IF(J7C.GE.I7C)PRINT 1020,I7,(J7(I),I=I7C,J7C)
      IF(J8C.GE.I8C)PRINT 1020,I8,(J8(I),I=I8C,J8C)
      IF(J9C.GE.I9C)PRINT 1020,I9,(J9(I),I=I9C,J9C)
      IF(JDEL.GE.IDEL)PRINT 1021,((LDEL(I,J),J=1,2),I=IDEL,JDEL)
      IF(JWC.GE.IWC)PRINT 1022,((KW(J,I),J=1,6),I=IWC,JWC)
      I6C=J6C+1
      I7C=J7C+1
      I8C=J8C+1
      I9C=J9C+1
      IDEL=JDEL+1
      IWC=JWC+1
      RETURN
      END
cc
c     *********************** 
      subroutine search(find) 
      implicit real*8(a-h,o-z)
c     *********************** 
c  ***this subroutine locates circuits or loops of order nc.npoint(nc)
c  ***are the indices of the points(triads) pertaining to the first
c  ***such loop found.
c  ***npart is the number of separate parts(groups of contiguous points)
c  ***on the axis of the flat graph.iparts is the number of points in 
c  ***the smallest part.ipartl is the number of points in the largest 
c  ***part.
c  ***this subroutine finds all the possible loops of order 3 and 4.for
c  ***nc.ge.5,it looks for only those who are partitionned in npart.le.2
c  ***which can eventually reduce tp a loop of order 4 without breaking
c  ***the basic structure of the flat graph. icross=-1,if lines cross 
c-------------------------------------------------------------------- 
c
c........corrected for free ends at last triad 11/23/93
*
      PARAMETER(KFL1=60,KFL2=12,KFL2A=2*KFL2,KFL2B=4*KFL2,
     +          KFL2C=2*KFL2+2)
      parameter(msum=10)
c
      logical find
      integer arr,tab1
      character*6 name
c
*
      COMMON/GRAPH/JDIAG(KFL2B,3),ARR(KFL2B,3),TAB1(KFL1,2),IL(KFL2B),
     + IH(KFL2B),NPOINT(KFL2C),NBNODE,IFIRST,ILAST,IPARTS,IPARTL,NPART,
     + ICROSS,NFREE,ITFREE(KFL2A),NFIN,NC
c
      data name/'search'/
c
 1000 format(' error in search.i,i1,i2,i3,npart,ipart,nc=',7i5)
c
c
c
c  ***initialization
c
      find=.false.
      ncm1=nc-1
      ncm=nc-2
      icross=0
c
c  ***first are treated two cases that do not involve do loops
c  ***1.one isolated point,either the first or the last
c
      npart=1
      ipartl=nc-1
      iparts=1
c
c  ***a.first
c
      i1=ifirst
      ja=jdiag(i1,1)
      jc=jdiag(i1,3)
c
      if(ja.eq.jc) then
        if(nc.gt.1) go to 800 
        npoint(1)=i1
        go to 900
      endif
c
      i2=tab1(ja,2)
      i3=tab1(jc,2)
      idist=iabs(il(i3)-il(i2))
c
      if(idist.lt.ncm) go to 800
c
      if(idist.gt.ncm.or.(i2.gt.ilast.and.i3.gt.ilast))then
c
c  ***b.last
c
        i1 = ilast
        ja=jdiag(i1,1)
        jc=jdiag(i1,2)
c
      if(ja.eq.jc) then
        if(nc.gt.1) go to 800 
        npoint(1)=i1
        go to 900
      endif
        i2=tab1(ja,1)
        i3=tab1(jc,1)
c.......11/23/93 correction for case of free ends
          if(i2.eq.i1) i2=tab1(ja,2)
          if(i3.eq.i1) i3=tab1(jc,2)
c.......4/13/94 case of inactive triads before first
          if(i2.lt.ifirst.or.i3.lt.ifirst)go to 250
      idist=iabs(il(i3)-il(i2))
      if(idist.lt.ncm) go to 800
c.........inactive triads beyond ilast
      if(idist.gt.ncm.or.(i2.gt.i1.or.i3.gt.i1)) go to 250
      endif
c  ***..found a loop
      ic=1
      npoint(ic)= i1
      i20=min0(i2,i3)
      i21=il(i20)
      i31=i21+ncm1
c
      do 203 ii=i21,i31
        ic=ic+1
        npoint(ic)=ih(ii)
  203 continue
c
      if(nc.le.2) then
        if(jdiag(ifirst,1).ne.jdiag(ilast,1))call phase(i1,jdiag,KFL2B)
        go to 900
      endif
c
      if(i1.ne.ilast) then
        it=i2
        jt=jdiag(ilast,2)
        k4=2
        i4=ilast
      else
        it=i3
        jt=jdiag(ifirst,3)
        k4=3
        i4=ifirst
      endif
c
      if(it.eq.i20)call phase(i1,jdiag,KFL2B)
      if(jt.eq.ja.or.jt.eq.jc)call change(i4,k4)
      go to 900
c
c  ***2.two isolated points,first and last.
c
  250 if(nc.eq.1)return
      if(nc.le.3) go to 100
      ipartl=nc-2
      iparts=1
      i1=ifirst
      i2=ilast
      ja=jdiag(i1,1)
      jb=jdiag(i1,3)
c
      if(tab1(ja,2).ne.i2) then
        ja=jdiag(i1,3)
        jb=jdiag(i1,1)
        if(tab1(ja,2).ne.i2) go to 100
      endif
c
      if(ja.eq.jdiag(i2,1)) then
        jc=jdiag(i2,2)
      else
        jc=jdiag(ilast,1)
      endif
c
      i3=tab1(jb,2) 
      i4=tab1(jc,1) 
      idist=il(i4)-il(i3)
c
      if(iabs(idist)-(ncm-1) .lt. 0) go to 800
      if(iabs(idist)-(ncm-1) .eq. 0) then
        npoint(1)= ilast
        npoint(2) = ifirst
        icross=isign(1,idist) 
        ic=2
        i20=min0(i3,i4)
        i21=il(i20) 
        i31=i21+ncm 
c
        do 261 ii=i21,i31
          ic=ic+1
          npoint(ic) = ih(ii) 
  261   continue
c
        if(ja.eq.jdiag(ifirst,1))call change(ifirst,3)
        if(ja.eq.jdiag(ilast,1))call change(ilast,2)
        go to 900
      endif
c
c  ***first general case:all points in one group
c
  100 npart=1
      iparts=0
      ipartl=nc
      k3=1
c
      do 101 in=1,nbnode
        i=ih(in)
  108   ja=jdiag(i,k3)
        if(i.ne.tab1(ja,2))then
          i2=tab1(ja,2)
c
          if(il(i2)-in-ncm1 .lt. 0) go to 800
          if(il(i2)-in-ncm1 .eq. 0)then 
            i21=il(i2)
            ic=0
c
            do 103 ii=in,i21
             ic=ic+1
             npoint(ic)=ih(ii)
  103       continue
c
            if(ja.eq.jdiag(ifirst,3))call change(ifirst,3)
            if(ja.eq.jdiag(ilast,2))call change(ilast,2)
            go to 900
        endif
      endif
c
        if(in.eq.1) then
          if(k3.ne.3) then
            k3=3
            go to 108
          else
            k3=1
          endif
        endif
c
  101 continue
c
c  ***search did not find loop nc.le.3
c
      if(nc.le.3) return
c
c  ***general case of loop partitionned in 2 groups.do loop 
c  ***on iparts
c
      npart=2
      nc2=nc/2
      k3=1
      k2=1
c
      do 400 ips=2,nc2
        jps=ips-1
        nbn=nbnode-jps
c
      do 301 i1=1,nbn
        i=ih(i1)
        i2=ih(i1+jps)
  302   ja=jdiag(i,k3)
        jd=jdiag(i2,k2)
c
        if(i.eq.tab1(ja,1)) then
          ii2=tab1(jd,2)
          ii1=tab1(ja,2)
        else
          ii1=tab1(ja,1)
          ii2=tab1(jd,1)
        endif
c
        idist=il(ii1)-il(ii2) 
c
        if(iabs(idist)-(ncm-jps) .lt. 0) go to 800
        if(iabs(idist)-(ncm-jps) .gt. 0) go to 320
  306   icross=isign(1,idist) 
        ic=0
        i21=il(i2)
c
      do 310 ii=i1,i21
        ic=ic+1
        npoint(ic)=ih(ii)
  310 continue
c
      i20=min0(ii1,ii2)
      i30=max0(ii1,ii2)
      i21=il(i20)
      i31=il(i30)
c
      do 311 ii=i21,i31
        ic=ic+1
        npoint(ic)=ih(ii)
  311 continue
c
      iparts=ips
      ipartl=nc-ips 
      if(jdiag(ifirst,3).eq.ja.or.jdiag(ifirst,3).eq.jd)call
     +change(ifirst,3)
      if(jdiag(ilast,2).eq.ja.or.jdiag(ilast,2).eq.jd)call
     +change(ilast,2)
      go to 900
c
  320 if(i1.eq.1) then
        if(k3.eq.3) then
          k3=1
          go to 301 
        else
          k3=3
          go to 302 
      endif
      endif
c
      if(i2.eq.ilast)then
       if(k2.ne.2) then
        k2=2
        go to 302
       endif
      endif
c
  301 continue
  400 continue
c
c  ***search did not find circuit of order nc
c
      return
c
c  ***loop found
c
  900 find=.true.
      call printj(name,msum)
c
      return
c
c  ***error printout
c
  800 print 1000,i,i1,i2,i3,npart,iparts,nc
      call printj(name,8)
      stop
c
      end 

************************************************************************
      SUBROUTINE SETDM
*
*     SET DIMENSIONS OF ARRAYS.
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20)
*
      LOGICAL SUMVAR
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),KW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
*
      COMMON/DIM/J6CC,J7CC,J8CC,J9CC,JWCC,JDELC
*
      JWCC=JWC
      JDELC=JDEL
      J6CC=J6C
      J7CC=J7C
      J8CC=J8C
      J9CC=J9C
      RETURN
      END
************************************************************************
      SUBROUTINE SETTAB(FAIL)
*
*     BUILDS UP THE UNSTRUCTURED GRAPH
*     SETS THE ARRAY J23,CONTAINING THE TWO LISTS OF ORIGINAL TRIADS
*     J2 AND J3,AND THE CORRESPONDING ARROWS ON THE ANGULAR MOMENTA
*     LINES.ALSO ESTABLISHES THE NUMERICAL AND PHASE FACTORS CONNECTING
*     RECOUPLING COEFFICIENT AND GRAPHS,ACCORDING TO YUTSIS,LEVINSON AND
*     VANAGAS.FOR THIS PURPOSE DETERMINES THE TOTAL J
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,KFL2=12,KFL2A=2*KFL2,KFL2B=4*KFL2,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20)
*
      LOGICAL FAIL
*
      LOGICAL TABS
      INTEGER ARROW
      COMMON/TREE/J23(KFL2A,3),ARROW(KFL2A,3),LINE(KFL1,2),LCOL(KFL1,2),
     +  TABS(KFL2A),NBTR
*
      LOGICAL FREE
      COMMON/COUPLE/M,N,J1(KFL1),J2(KFL2,3),J3(KFL2,3),FREE(KFL1)
*
      LOGICAL SUMVAR
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),KW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
*
      COMMON/BUILD/IAL(KFL2B),IF1,IF2,NODE
*
      CHARACTER*6 NAME
      DATA NAME/'SETTAB'/
*
      do 501 i = 1,kfl2a
      do 501 j = 1,3
  501 j23(i,j) = 0
      IPR=N-1
      NBTR=IPR+IPR
      DO 4 I=1,IPR
      DO 5 J=1,2
      J23(I,J)=J2(I,J)
      ARROW(I,J)=1
    5 CONTINUE
      TABS(I)=.FALSE.
      J23(I,3)=J2(I,3)
      ARROW(I,3)=-1
    4 CONTINUE
      IPR1=IPR+1
      DO 7 I=IPR1,NBTR
      II=I-IPR
      DO 6 J=1,2
      J23(I,J)=J3(II,J)
      ARROW(I,J)=-1
    6 CONTINUE
      TABS(I)=.FALSE.
      J23(I,3)=J3(II,3)
      ARROW(I,3)=1
    7 CONTINUE
      DO 11 J=1,NBTR
   11 J8(J)=J23(J,1)
      J8C=NBTR+IPR
      NB1=NBTR+1
      DO 12 J=NB1,J8C
      I=J-IPR
   12 J8(J)=J23(I,3)
      J6C=NBTR
      DO 13 J=1,J6C
   13 J6(J)=J23(J,3)
      DO 10 I=1,M
      SUMVAR(I)=.FALSE.
   10 IAL(I)=1
      DO 9 I=1,NBTR
      DO 8 J=1,3
      JI=J23(I,J)
      K=IAL(JI)
      LINE(JI,K)=I
      LCOL(JI,K)=J
      IAL(JI)=K+1
    8 CONTINUE
    9 CONTINUE
      IT=0
      DO 18 I=1,NBTR
      JT=J23(I,3)
      IF(IAL(JT).EQ.3)GO TO 17
      IF(IT.EQ.1)GO TO 16
      JT1=JT
      IT=1
      GO TO 18
   16 CALL DELTA(JT1,JT,FAIL)
      IF(FAIL)GO TO 20
      K=LINE(JT,1)
      KC=LCOL(JT,1)
      LINE(JT1,2)=K
      LCOL(JT1,2)=KC
      J23(K,KC)=JT1
      IAL(JT)=1
      GO TO 19
   17 CALL OTHERJ(I,JT,L,LC,K)
      IF(LC.EQ.3)GO TO 19
   18 CONTINUE
   19 J9(J9C+1)=JT
      J9C=J9C+2
      J9(J9C)=JT
   20 CALL PRINTJ(NAME,4)
      RETURN
      END
************************************************************************
      SUBROUTINE SPRATE(M)
*
*    THIS SUBROUTINE PREPARES THE INFORMATION TO BE TRANSFERED TO
*    GENSUM FOR NUMERICAL EVALUATION.
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20,
     +          KFLS=12,KFLN=10,KFLV=40)
*
      LOGICAL SUM6J,T6J,JT,JS
      DIMENSION JTEM4(KFLS,KFLW),JTEM5(KFLS,KFLW),JTEM6(KFLS),
     +  NSUM6J(KFLW),J6SUM(KFLW)
      DIMENSION SUM6J(KFLW),T6J(KFLW),JT(KFLS),JS(KFLS),INVER(KFL1),
     +  JNSUM(KFLS),JINV(KFLS),N6JN(KFLW),IN6J(KFLW),JSUMT(KFLW,6)
*
      LOGICAL SUMVAR
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),JW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
*
      COMMON/DIM/J6CC,J7CC,J8CC,J9CC,JWCC,JDELC
*
      COMMON/SUMARG/J6P(KFLV),J7P(KFLV),J8P(KFLV),J9P(KFLV),
     + JWORD(6,KFLW),NLSUM,
     + NBJ(KFLN),NB6J(KFLN),K6CP(KFLN),K7CP(KFLN),K8CP(KFLN),K9CP(KFLN),
     + JSUM6(KFLS),JSUM4(KFLS,KFLW),JSUM5(KFLS,KFLW),INV6J(KFLW)
*
      LOGICAL CUT
      COMMON/CUTDIG/CUT
*
      CHARACTER*4 NAME
*
*

*    1. TEST THAT ARRAY DIMENSIONS HAVE NOT BEEN EXCEEDED.
*
      IF(MP.LE.KFL1)GO TO 49
      NMX=KFL1
      NPX=MP
      NAME='KFL1'
      GO TO 60
   49 IF(JWC.LE.KFLW)GO TO 50
      NMX=KFLW
      NPX=JWC
      NAME='KFLW'
      GO TO 60
   50 IF(J6C.LE.KFL6)GO TO 51
      NMX=KFL6
      NPX=J6C
      NAME='KFL6'
      GO TO 60
   51 IF(J7C.LE.KFL7)GO TO 52
      NMX=KFL7
      NPX=J7C
      NAME='KFL7'
      GO TO 60
   52 IF(J8C.LE.KFL8)GO TO 53
      NMX=KFL8
      NPX=J8C
      NAME='KFL8'
      GO TO 60
   53 IF(J9C.LE.KFL9)GO TO 54
      NMX=KFL9
      NPX=J9C
      NAME='KFL9'
   60 PRINT 1000,NAME,NPX,NMX
 1000 FORMAT(2X,'DIMENSION ERROR FOR  ',A4,I5,' IS OUT OF ALLOWED RAN
     +GE',I3)
      STOP
*
*     2. DETERMINATION OF EFFECTIVE SUMMATION VARIABLES AND THEIR
*     RELATIONSHIPS WITH 6J COEFFICIENTS.
*
   54 DO 2 I=1,JWC
      INV6J(I)=0
    2 SUM6J(I)=.FALSE.
      NSUM=0
      NLSUM=0
      IF(MP.EQ.M)RETURN
      M1=M+1
      DO 1 I=M1,MP
      IF(.NOT.SUMVAR(I))GO TO 1
      NSUM=NSUM+1
      JSUM6(NSUM)=0
      INVER(I)=NSUM
    1 CONTINUE
      IF(NSUM.EQ.0)RETURN
      IF(NSUM.LE.KFLS)GO TO 70
      NMX=KFLS
      NPX=NSUM
      NAME='NSUM'
      GO TO 60
   70 KT=0
      DO 4 I=1,JWC
      DO 5 J=1,6
      IK=JW(J,I)
      IF(.NOT.SUMVAR(IK))GO TO 5
      IF(SUM6J(I))GO TO 6
      SUM6J(I)=.TRUE.
      KT=KT+1
      J6SUM(KT)=0
      NSUM6J(KT)=I
      INV6J(I)=KT
    6 ISK=INVER(IK)
      I2=JSUM6(ISK)+1
      JSUM6(ISK)=I2
      JSUM4(ISK,I2)=J
      JSUM5(ISK,I2)=KT
      I3=J6SUM(KT)+1
      J6SUM(KT)=I3
      JSUMT(KT,I3)=ISK
    5 CONTINUE
    4 CONTINUE
      CALL VAR(J6,J6P,J6C,J6CP,J6CC,SUMVAR,MP,M,INVER)
      CALL VAR(J7,J7P,J7C,J7CP,J7CC,SUMVAR,MP,M,INVER)
      CALL VAR(J8,J8P,J8C,J8CP,J8CC,SUMVAR,MP,M,INVER)
      CALL VAR(J9,J9P,J9C,J9CP,J9CC,SUMVAR,MP,M,INVER)
      IF(CUT)GO TO 17
      NLSUM=1
      NBJ(1)=NSUM
      NB6J(1)=KT
      K6CP(1)=J6CP
      K7CP(1)=J7CP
      K8CP(1)=J8CP
      K9CP(1)=J9CP
      DO 21 I=1,KT
      I1=NSUM6J(I)
      DO 22 J=1,6
      JWORD(J,I)=JW(J,I1)
   22 CONTINUE
   21 CONTINUE
      DO 80 I=1,NSUM
      ISU=JSUM6(I)
      DO 81 J=1,ISU
      I1=JSUM5(I,J)
      J1=JSUM4(I,J)
      JWORD(J1,I1)=MP+I
   81 CONTINUE
   80 CONTINUE
      GO TO 25
*
*     3.SEPARATION OF VARIABLES AND SUMS IN CASE A CUT WAS DETECTED.
*
   17 K6C=0
      K7C=0
      K8C=0
      K9C=0
      NJ=0
      N6J=0
      DO 9 I=1,KT
    9 T6J(I)=.FALSE.
      DO 7 I=1,NSUM
      JT(I)=.FALSE.
    7 JS(I)=.FALSE.
      J=1
   10 NJ=NJ+1
      JNSUM(NJ)=J
      JINV(J)=NJ
      JT(J)=.TRUE.
   18 JS(J)=.TRUE.
      JS6=JSUM6(J)
      DO 11 I=1,JS6
      I6J=JSUM5(J,I)
      IF(T6J(I6J))GO TO 14
      T6J(I6J)=.TRUE.
      N6J=N6J+1
      N6JN(N6J)=NSUM6J(I6J)
      IN6J(I6J)=N6J
   14 J6J=J6SUM(I6J)
      DO 12 K=1,J6J
      JK=JSUMT(I6J,K)
      IF(JT(JK))GO TO 12
      NJ=NJ+1
      JNSUM(NJ)=JK
      JINV(JK)=NJ
      JT(JK)=.TRUE.
   12 CONTINUE
   11 CONTINUE
      DO 13 JJ=1,NSUM
      J=JJ
      IF(JS(JJ).OR..NOT.JT(JJ))GO TO 13
      GO TO 18
   13 CONTINUE
   15 NLSUM=NLSUM+1
      IF(NLSUM.LE.KFLN)GO TO 71
      NMX=KFLN
      NPX=NLSUM
      NAME='KFLN'
      GO TO 60
   71 NBJ(NLSUM)=NJ
      NB6J(NLSUM)=N6J
      IF(J6CP.EQ.0)GO TO 30
      CALL CHVAR(J6P,J6CP,K6C,JT,JINV,NSUM)
   30 K6CP(NLSUM)=K6C
      IF(J7CP.EQ.0)GO TO 31
      CALL CHVAR(J7P,J7CP,K7C,JT,JINV,NSUM)
   31 K7CP(NLSUM)=K7C
      IF(J8CP.EQ.0)GO TO 32
      CALL CHVAR(J8P,J8CP,K8C,JT,JINV,NSUM)
   32 K8CP(NLSUM)=K8C
      IF(J9CP.EQ.0)GO TO 33
      CALL CHVAR(J9P,J9CP,K9C,JT,JINV,NSUM)
   33 K9CP(NLSUM)=K9C
      IF(NJ.EQ.NSUM)GO TO 20
      DO 16 JJ=1,NSUM
      J=JJ
      IF(.NOT.JT(JJ))GO TO 10
   16 CONTINUE
   20 DO 26 I=1,KT
      I1=N6JN(I)
      DO 27 J=1,6
   27 JWORD(J,I)=JW(J,I1)
   26 CONTINUE
      DO 28 I=1,NSUM
      IK=JNSUM(I)
      I2=JSUM6(IK)
      JTEM6(I)=I2
      DO 29 J=1,I2
      JTEM4(I,J)=JSUM4(IK,J)
      K=JSUM5(IK,J)
      JTEM5(I,J)=IN6J(K)
   29 CONTINUE
   28 CONTINUE
      DO 40 I=1,NSUM
      I2=JTEM6(I)
      JSUM6(I)=I2
      DO 41 J=1,I2
      I1=JTEM5(I,J)
      J1=JTEM4(I,J)
      JSUM4(I,J)=J1
      JSUM5(I,J)=I1
      JWORD(J1,I1)=I+MP
   41 CONTINUE
   40 CONTINUE
   25 RETURN
      END
************************************************************************
      SUBROUTINE SQUARE
*
*    REDUCES A CIRCUIT OF ORDER 4 IN THE TWO CASES WHICH ARE LEFT
*    OVER BY POLYGN,NAMELY TWO DISCONNECTED GROUPS OF TWO POINTS
*    AND ONE GROUP OF TWO POINTS PLUS THE TWO ENDS OF THE AXIS.IN
*    THE LATTER, THE END OF THE AXIS IS TRANSFERRED TO THE BEGINNING.
*    IN THIS PROCESS,ONE SUMMATION VARIABLE AND TWO 6J SYMBOLS ARE
*    INTRODUCED.
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,KFL2=12,KFL2A=2*KFL2,KFL2B=4*KFL2,
     +          KFL2C=2*KFL2+2,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20)
*
      INTEGER ARR,TAB1
      COMMON/GRAPH/JDIAG(KFL2B,3),ARR(KFL2B,3),TAB1(KFL1,2),IL(KFL2B),
     + IH(KFL2B),NPOINT(KFL2C),NBNODE,IFIRST,ILAST,IPARTS,IPARTL,NPART,
     + ICROSS,NFREE,ITFREE(KFL2A),NFIN,NC
      LOGICAL SUMVAR
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),KW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
*
      CHARACTER*6 NAME,NAMSUB
      COMMON/NAM/NAMSUB
      DATA NAME/'SQUARE'/
*
      NAMSUB = NAME
      MP = MP+1
      SUMVAR(MP) = .TRUE.
      IT1 = NPOINT(1)
      IT2 = NPOINT(2)
*
      IF (ICROSS .EQ. 1) THEN
         IT3 = NPOINT(3)
         IT4 = NPOINT(4)
         K23 = 3
         K32 = 2
      ELSE
         IT3 = NPOINT(4)
         IT4 = NPOINT(3)
         K23 = 2
         K32 = 3
      ENDIF
*
      L4 = JDIAG(IT2,1)
*
      IF (ARR(IT2,1) .LE. 0) THEN
         CALL PHASE2 (L4)
         ARR(IT2,1) = 1
         ARR(IT3,1) = -1
      ENDIF
*
      L2 = JDIAG(IT1,1)
      IF (ARR(IT1,1) .GT. 0) CALL PHASE2(L2)
      JWC = JWC+1
      KW(1,JWC) = L4
      KW(2,JWC) = L2
      KW(3,JWC) = JDIAG(IT2,2)
      JJ1 = JDIAG(IT1,3)
      KW(4,JWC) = JJ1
      KW(5,JWC) = MP
      KW(6,JWC) = JDIAG(IT1,2)
      IF (ARR(IT1,2) .LT. 0) CALL PHASE2 (JDIAG(IT1,2))
      JWC = JWC+1
      KW(1,JWC) = L4
      KW(2,JWC) = L2
      JJ3 = JDIAG(IT3,K23)
      JJ2 = JDIAG(IT4,K32)
      KW(3,JWC) = JJ3
      KW(4,JWC) = JJ2
      KW(5,JWC) = MP
      KW(6,JWC) = JDIAG(IT3,K32)
      IF (ARR(IT3,K32) .LT. 0) CALL PHASE2 (JDIAG(IT3,K32))
      J6(J6C+1) = MP
      J6C = J6C+2
      J6(J6C) = MP
*
      IF (NPART .EQ. 1) THEN
         ITMIN = IT2
         ITMAX = IT3
      ELSE
         ITMIN = MIN (IT2,IT3)
         ITMAX = MAX (IT2,IT3)
      ENDIF
      ITMN = MIN (IT1,IT4)
      ITMX = MAX (IT1,IT4)
*
      TAB1(MP,1) = ITMIN
      TAB1(MP,2) = ITMAX
      JDIAG(IT2,1) = MP
      JDIAG(IT3,1) = MP
      JDIAG(IT2,3) = JJ1
      ARR(IT2,3) = ARR(IT1,3)
      JDIAG(IT3,K32) = JJ2
      ARR(IT3,K32) = ARR(IT4,K32)
*
      IF (ICROSS .EQ. 1) THEN
         J7(J7C+1) = L2
         J7(J7C+2) = L4
         CALL PHASE2 (L4)
         J7C = J7C+3
         J7(J7C) = MP
      ELSE
         CALL PHASE2 (JJ2)
      ENDIF
*
      ITLL = IL(ITMN)
      ITHL = IL(ITMX)
*
      DO 5 I = ITLL+1,ITHL-1
         IT = IH(I)
         ILP = I-1
         IL(IT) = ILP
         IH(ILP) = IT
    5 CONTINUE
      IF (ITHL .NE. NBNODE) THEN
         DO 6 I = ITHL+1,NBNODE
            IT = IH(I)
            ILP = I-2
            IL(IT) = ILP
            IH(ILP) = IT
    6    CONTINUE
      ENDIF
*
      IF (NPART .NE. 2) THEN
         TAB1(JJ1,1) = IH(1)
         TAB1(JJ1,2) = IH(NBNODE-2)
      ENDIF
*
      RETURN
      END
************************************************************************
      SUBROUTINE TRDEL(JJ1,JJ2,JJ3,NBN,FAIL)
*
*     TEST FOR TRIANGULAR DELTA.IF NOT SATISFIED FAIL=.TRUE.
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,KFL2=12,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20)
*
      LOGICAL FAIL
      LOGICAL SUMVAR
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),KW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
      LOGICAL CUT
      COMMON/CUTDIG/CUT
      LOGICAL FREE
      COMMON/COUPLE/M,N,J1(KFL1),J2(KFL2,3),J3(KFL2,3),FREE(KFL1)
*
      IF(SUMVAR(JJ1).OR.SUMVAR(JJ2).OR.SUMVAR(JJ3))GO TO 1
      IF(NBN.GT.4)CUT=.TRUE.
      IF(FREE(JJ1).OR.FREE(JJ2).OR.FREE(JJ3))GO TO 1
      I1=J1(JJ1)
      I2=J1(JJ2)
      I3=J1(JJ3)
      IF(I1.GE.(IABS(I2-I3)+1).AND.I1.LE.(I2+I3-1))GO TO 1
      FAIL=.TRUE.
    1 RETURN
      END
************************************************************************
      SUBROUTINE TRIANG(FAIL)
*
*    REDUCES A TRIANGLE HAVING ONE APEX AT EITHER END OF THE AXIS OF
*    THE FLAT DIAGRAM.
*    THIS INTRODUCES ONE 6J SYMBOL AND SOME PHASE FACTORS .
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,KFL2=12,KFL2A=2*KFL2,KFL2B=4*KFL2,
     +          KFL2C=2*KFL2+2,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20)
*
      LOGICAL FAIL
*
      INTEGER ARR,TAB1
      COMMON/GRAPH/JDIAG(KFL2B,3),ARR(KFL2B,3),TAB1(KFL1,2),IL(KFL2B),
     + IH(KFL2B),NPOINT(KFL2C),NBNODE,IFIRST,ILAST,IPARTS,IPARTL,NPART,
     + ICROSS,NFREE,ITFREE(KFL2A),NFIN,NC
      LOGICAL SUMVAR
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),KW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
*
      CHARACTER*6 NAME,NAMSUB
      COMMON/NAM/NAMSUB
      DATA NAME/'TRIANG'/
*
      NAMSUB=NAME
      IT1=NPOINT(1)
      IT2=NPOINT(2)
      IT3=NPOINT(3)
      JWC=JWC+1
      KW(1,JWC)=JDIAG(IT3,2)
      KW(2,JWC)=JDIAG(IT2,3)
      KW(3,JWC)=JDIAG(IT3,1)
      IF(ARR(IT3,1).GT.0)CALL PHASE2(KW(3,JWC))
      KW(4,JWC)=JDIAG(IT2,1)
      IF(ARR(IT2,1).LT.0)CALL PHASE2(KW(4,JWC))
      K23=3
      IF(IT1.EQ.IFIRST)K23=2
      KW(5,JWC)=JDIAG(IT1,K23)
      KW(6,JWC)=JDIAG(IT3,3)
      CALL TRDEL(KW(1,JWC),KW(2,JWC),KW(5,JWC),NBNODE,FAIL)
      IF(FAIL)GO TO 15
      IF(ARR(IT3,3).GT.0)CALL PHASE2(KW(6,JWC))
      JT1=KW(5,JWC)
      JDIAG(IT3,1)=JT1
      JDIAG(IT3,3)=KW(2,JWC)
      ARR(IT3,1)=ARR(IT1,K23)
      ARR(IT3,3)=ARR(IT2,3)
      IF(IT1.EQ.IFIRST)GO TO 10
      TAB1(JT1,1)=IT3
      TAB1(JT1,2)=IH(NBNODE-1)
      K12=1
      GO TO 11
   10 TAB1(JT1,1)=IH(2)
      TAB1(JT1,2)=IT3
      K12=2
   11 IL3=IL(IT3)
      IF(IT1.EQ.ILAST)GO TO 13
      IL2=IL(IT2)-1
      DO 2 I=2,IL2
      IT=IH(I)
      ILP=I-1
      IL(IT)=ILP
      IH(ILP)=IT
    2 CONTINUE
   13 DO 1 I=IL3,NBNODE
      IT=IH(I)
      ILP=I-K12
      IL(IT)=ILP
      IH(ILP)=IT
    1 CONTINUE
   15 RETURN
      END
************************************************************************
      SUBROUTINE VAR(JN,JNS,JNC,JNSC,JBC,SUMVAR,MP,M,INVER)
*
*    TEST FOR VARIABLE CHARACTER AND PUT IN JNS IF YES,AND JN NOW
*    CONTAINS 0.
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFLV=40)
*
      LOGICAL SUMVAR(MP)
      DIMENSION JN(JNC),JNS(KFLV),INVER(MP)
      JNSC=0
      IF(JBC.EQ.JNC)GO TO 2
      JBBC=JBC+1
      DO 1 I=JBBC,JNC
      I1=JN(I)
      IF(.NOT.SUMVAR(I1))GO TO 1
      JNSC=JNSC+1
      IF(JNSC.GT.KFLV) THEN
          PRINT 100,JNSC,KFLV
  100     FORMAT(1X,'DIMENSION ERROR IN VAR. JNSC=',I5,' KFLV=',I5)
          STOP
      ENDIF
      J=INVER(I1)
      JNS(JNSC)=J
      JN(I)=M
    1 CONTINUE
    2 RETURN
      END
************************************************************************
      SUBROUTINE WAY(L,KA,KB,ICH,NB)
*
*
*     TESTS ONE STEP FORWARD  IF THE WAY IS FREE.FIRST AND SECOND
*     ARGUMENTS ARE INTERCHANGED OR NOT ACCORDING TO ICH=-1,OR +1
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,KFL2=12,KFL2A=2*KFL2,KFL2B=4*KFL2)
*
      LOGICAL TABS
      INTEGER ARROW
      COMMON/TREE/J23(KFL2A,3),ARROW(KFL2A,3),LINE(KFL1,2),LCOL(KFL1,2),
     +  TABS(KFL2A),NBTR
*
      COMMON/BUILD/IAL(KFL2B),IF1,IF2,NODE
*
      K1=J23(L,KA)
      K2=J23(L,KB)
      NB=IAL(K1)+IAL(K2)-1
      IF(NB)3,2,8
    2 NB1=IAL(K1)-IAL(K2)
      IF(NB1)9,8,8
    3 CALL OTHERJ(L,K1,L1,LC1,LA)
      CALL OTHERJ(L,K2,L2,LC2,LB)
      CALL NEIBOR(LC1,I1,I2)
      CALL NEIBOR(LC2,I3,I4)
      JI1=J23(L1,I1)
      JI2=J23(L1,I2)
      JI3=J23(L2,I3)
      JI4=J23(L2,I4)
      IA=IAL(JI1)+IAL(JI2)
      IB=IAL(JI3)+IAL(JI4)
      NBP=IB+IA+1
      NBM=IB-IA
      GO TO (8,4,5,4,6),NBP
    4 IF(NBM)9,8,8
    5 IF(NBM)9,6,8
    6 IF(JI3.EQ.IF1.OR.JI3.EQ.IF2.OR.JI4.EQ.IF1.OR.JI4.EQ.IF2)GO TO 9
    8 ICH=1
      GO TO 10
    9 ICH=-1
   10 RETURN
      END
************************************************************************
      subroutine zero(j,jz,fail)
************************************************************************
*
*  ***suppresses one line and two nodes of the unstructured graph
*  ***introduces  zeros in the triads j23.as a consequence the other
*  ***two arguments of the triad are put equal.if there was already
*  ***a zero in the triad which is changed,it is a special case.
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(KFL1=60,KFL2=12,KFL2A=2*KFL2,KFL2B=4*KFL2,
     +          KFLZ=30,
     +          KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20)
*
      LOGICAL FAIL
      COMMON/ZER/NZERO,JZERO(KFLZ)
      COMMON/BUILD/IAL(KFL2B),IF1,IF2,NODE
*
      LOGICAL TABS
      INTEGER ARROW
      COMMON/TREE/J23(KFL2A,3),ARROW(KFL2A,3),LINE(KFL1,2),LCOL(KFL1,2),
     +  TABS(KFL2A),NBTR
*
      LOGICAL FREE
      COMMON/COUPLE/M,N,J1(KFL1),J2(KFL2,3),J3(KFL2,3),FREE(KFL1)
*
      LOGICAL SUMVAR
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     + J9(KFL9),KW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
*
      LOGICAL CUT,NOCUT
      COMMON/CUTDIG/CUT
*
      character*6 name
      data name/'zero  '/
*
*
*
      nocut=.false.
      nzero=0
*
      if(j.ge.1)then
       call otherj(0,jz,lin,lc,k1)
        i=nzero
        go to 8
      endif
*
      do 11 i=1,m
        if(j1(i).ne.1.or.free(i).or.ial(i).le.1)go to 11
        nzero=nzero+1
        if (nzero .gt. kflz) then
           print 100, nzero, kflz
  100      format(1x,'Dimension error in ZERO: nzero=',I5,' KFLZ=',I5)
           stop
        endif
        jzero(nzero)=i
   11 continue
*
      nocut=.true.
      m=m+1
      j1(m)=1
      sumvar(m)=.false.
      free(m)=.false.
      if(nzero.eq.0)go to 7
      call printj(name,1)
      i=0
    1 i=i+1
      jz=jzero(i)
      j=0
   13 j=j+1
      lin=line(jz,j)
      if(tabs(lin))go to 2
      lc=lcol(jz,j)
    8 call neibor(lc,l1,l2)
      jj1=j23(lin,l1)
      jj2=j23(lin,l2)
*
      if(jj1.eq.jj2)then
        j6c=j6c+1
        j6(j6c)=jj1 
        lo1=lin
        lo2=lin
        lco1=l1
        lco2=l2
        go to 10
      endif
*
      call delta(jj1,jj2,fail)
      if(fail)go to 7
*
      if(j1(jj1).ne.1.and.j1(jj2).ne.1)go to 15
      if(j1(jj1) .lt. j1(jj2))go to 15
      if(j1(jj1) .gt. j1(jj2))go to 19
*
      if(nzero.ne.0)then
        do 17 jjx=i,nzero
          jjz=jzero(jjx)
          if(jj1 .eq. jjz)go to 15
   18     if(jj2 .eq. jjz)go to 19
   17   continue
      endif
*
      go to 15
*
   19 jjz=jj2
      jj2=jj1
      jj1=jjz
*
   15 call otherj(lin,jj1,lo1,lco1,k1)
      call otherj(lin,jj2,lo2,lco2,k2)
      j9c=j9c+1
      j9(j9c)=jj1
      j23(lo2,lco2)=jj1
      line(jj1,k1)=lo2
      lcol(jj1,k1)=lco2
*
   10 if(arrow(lin,l1) .lt. arrow(lin,l2)) then
        call phase2(jj1)
      else
      if(arrow(lin,l1) .eq. arrow(lin,l2)) then
        arrow(lo1,lco1)=1
        arrow(lo2,lco2)=-1
      endif
      endif
*
      tabs(lin)=.true.
      nbtr=nbtr-1
      if(nbtr.eq.0)go to 7
*      if(lo1.ne.lo2)go to 2
*      l=6-lco1-lco2
*      jt=j23(lo1,l)
*      if(j1(jt).eq.1.and..not.free(jt))go to 2
*      call delta(jt,m,fail)
*      if(fail)go to 7
*      call neibor(l,l1,l2)
*      jtf=j23(lo1,l1)
*      if(arrow(lo1,l1) .lt. arrow(lo1,l2))call phase2(jtf)
*      j6c=j6c+1
*      j6(j6c)=jtf
*      nbtr=nbtr-1
*      tabs(lo1)=.true.
*      call otherj(lo1,jt,lin,lc,k)
*      go to 8
*     november 22 1989
      if(lo1.eq.lo2)then
        l=6-lco1-lco2
        jt=j23(lo1,l)
        if(j1(jt).eq.1.and..not.free(jt))go to 2
        call delta(jt,m,fail)
        if(fail)go to 7
        nzero=nzero+1
        jzero(nzero)=jt
      end if
    2 if(j.eq.1)go to 13
*
      if (nbtr .ne. 0) then
        if(i.lt.nzero) go to 1
      endif
*
    7 call printj(name,4)
      if(nocut)cut=.false.
*
      return
      end
