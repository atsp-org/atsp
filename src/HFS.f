*=======================================================================
*
*      A General Program for Computing Magnetic Dipole and
*      Electric Quadrupole Hyperfine Constants 
*
*               C O P Y R I G H T -- 1994
*
*      P. Jonsson and C.G. Wahlstrom, Department of Physics,  
*                                     Lund Institute of Technology,
*                                     P.O. Box 118,S-221 00 Lund
*                                     Sweden 
*      C. Froese Fischer,             Department of Computer Science
*                                     Vanderbilt University,
*                                     Nashville, TN 37235
*                                     USA
*
*     Computer Physics Communication, Vol. 74, 399--414 (1993)
*======================================================================
*
*      The program allows for a limited degree of non-orthogonality
*      between the left and right configuration states.
*
*     The current limits are:
*     max NWD=30 different orbitals
*     max NCD=100 configurations (change the parameter NCD
*     in this routine and in the subroutines in this program unit)

      PROGRAM HFS 

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220,NWD=30,NCD=100,NORD=NWD*(NWD-1)/2)
      PARAMETER (EPSILO=1.D-9)
      CHARACTER*24 NAME(6)
      CHARACTER*1 ANS
      INTEGER SS,SS1,SS2
CSUN  REAL TIMES(2),DTIME

*                        Dimension statements

      DIMENSION ORBA(20),DIPA(20),CONTA(20),QUADB(20),WT(NCD,20)
      DIMENSION ORBF1(20),DIPF1(20),CONTF1(20),QUADF1(20),WTJIJF(20)

*                        Common blocks

      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/TENSOR/VVSHELL(50),IIHSH,NNOVLP(50),IIRHO(50),IISIG(50),
     : MMU(50),MMUP(50),NNU(50),NNUP(50)
      COMMON/COEF/ORBF(20),DIPF(20),CONTF(20),QUADF(20),J(20),JP(20),
     : JJMAX,JJMIN,JJ1MAX,JJ1MIN,JJ2MAX,JJ2MIN,LL1,SS1,LL2,SS2
      COMMON/PRINTA/ORB,DIP,CONT,QUAD,RADINT,AZSQR,CONTRI,INTGR,
     : IREPRI,ICOPRI,IDOPRI,OVLINT1,OVLINT2
      COMMON/FACT/GAM(100)
      COMMON/HYPER/NHY,MHY,IRHY(20),ISHY(20),JRHY(20),JSHY(20),VHY(20),
     1 IMUHY(20),IMUPHY(20),JMUHY(20),JMUPHY(20),
     2 INUHY(20),INUPHY(20),JNUHY(20),JNUPHY(20),NNNOVLP(20)
      COMMON/INFORM/IREADC,IWRITE,IOUT,IREADJ,IREADW,ISCW,IREAD,ISC(4)
      COMMON/MEDEFN/IHSH,NJ(10),LJ(10),NOSH1(10),NOSH2(10),J1QN1(19,3),
     :J1QN2(19,3),IJFUL(10)
      COMMON/OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
     1 ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
     2 IORTH(NORD)
      COMMON/OVRINT/IOVEL(2),IOVER(2),IOVEP(2)
      COMMON/STATES/NCFG,MAXORB,IAJCMP(NWD),LJCOMP(NWD),NJCOMP(NWD),
     :NOCCSH(NCD),NELCSH(5,NCD),NOCORB(5,NCD),J1QNRD(9,NCD)
      COMMON/EMS/IEM(4),JI,JF,IFL

*                       Data statements

      DATA ORBA,DIPA,CONTA,QUADB/80*0.D0/
      DATA AMHZ,BMHZ,GS/47.70534D0,234.9647D0,2.002319D0/
      DATA IEM/2HL1,2HCS,2HS1,2HC2/

*                       Output formats

610   FORMAT(///20X,'====================='/
     :          20X,'  H Y P E R F I N E  '/
     :          20X,'====================='//)
615   FORMAT(/20X,' Hyperfine structure calculation '/)
620   FORMAT(/20X,' The configuration set ')
630   FORMAT(/30X,'**',//20X,17HA factors in MHz ,
     : //,2X,8HJ     J'
     : ,8X,52HOrbital        Spin-dip        Cont           Total )
640   FORMAT(/,A4,2X,A4,4(2X,F13.5))
650   FORMAT(/20X,17HB factors in MHz ,
     : //,2X,8HJ     J'
     : ,8X,10HQuadrupole )
660   FORMAT(/,A4,2X,A4,2X,F13.5)
670   FORMAT(/30X,'**',//20X,28HHyperfine parameters in a.u.,
     : //20X,47Hal             ad             ac             bq ) 
680   FORMAT(/,9X,4(F15.5)/)

*     IBUG3   debug in recoupling package
*     NBUG6   debug in tensor package

*     Debug information
      INTGR = 0
      IBUG1 = 0
      IBUG2 = 0
      IBUG3 = 0
      NBUG6 = 0
      NBUG7 = 0
      IFULL = 0
      NFACTS = 32

CDBG  WRITE(ISCW,*) ' Input IBUG3, NBUG6 (0/1) '
CDBG  READ(IREAD,*) IBUG3,NBUG6

CSUN  RTIME=DTIME(TIMES)

*     The following section concerns input/output and may be 
*     system dependent. Check allowed unit numbers and file name
*     conventions - modify, if necessary

*     IWRITE  debug file
*     IREADC  configuration file
*     IREADJ  lsj or ls file
*     IREADW  wavefunction file
*     IOUT    output file
*     ISCW    unit number of standard error (screen)
*     IREAD   unit number of standard input

      IWRITE=6
      IREADC=4
      IREADJ=2
      IREADW=8
      IOUT=3
CDBG  ISC(1)=7
CTSS  ISCW=5
      ISCW=0
      IREAD=5

*     Write heading

      WRITE(ISCW,610)
10    WRITE(ISCW,*) 'Name of state'
      READ(IREAD,'(A)') NAME(1)
      K=INDEX(NAME(1),' ')
      IF (K.EQ.1) THEN
         WRITE(ISCW,*) 'Names may not start with a blank'
         GO TO 10
      ELSE
         NAME(1)(K:K+1)='.c'
         NAME(2)=NAME(1)(1:K-1)//'.j'
         NAME(3)=NAME(1)(1:K-1)//'.w'
         NAME(4)=NAME(1)(1:K-1)//'.h'
CDBG     NAME(5)=NAME(1)(1:K-1)//'.m'
         NAME(6)=NAME(1)(1:K-1)//'.l'

      ENDIF
      OPEN(UNIT=IREADC, FILE=NAME(1),STATUS='OLD')
      OPEN(UNIT=IOUT, FILE=NAME(4),STATUS='UNKNOWN')
CDBG  OPEN(UNIT=ISC(1), FILE=NAME(5),STATUS='UNKNOWN')

      WRITE(IOUT,615) 

20    WRITE(ISCW,'(A/A/A/A)') ' Indicate the type of calculation ',
     : ' 0 => diagonal A and B factors only;',
     : ' 1 => diagonal and off-diagonal A and B factors;',
     : ' 2 => coefficients of the radial matrix elements only:'
      READ(IREAD,*) ICALC
      IF (ICALC.EQ.2) THEN
         WRITE(IOUT,*) NAME(1)
         ANS='M'
      ELSE
33       WRITE(ISCW,*)'Hyperfine parameters al,ad,ac and bq ? (Y/N)'
         READ(IREAD,'(A)') ANS
         IF (ANS.EQ.'N'.OR.ANS.EQ.'n') THEN
            HFSPARAM=0
         ELSEIF (ANS.EQ.'Y'.OR.ANS.EQ.'y') THEN
            HFSPARAM=1 
         ELSE 
            GOTO 33
         ENDIF
30       WRITE(ISCW,*)'Input from an MCHF (M) or CI (C) calculation ?'
         READ(IREAD,'(A)') ANS
      ENDIF
      IF (ANS.EQ.'M'.OR.ANS.EQ.'m') THEN
         WRITE(IOUT,*) NAME(1)
         IMCHF=1
         NR=1
      ELSEIF (ANS.EQ.'C'.OR.ANS.EQ.'c') THEN
35       WRITE(ISCW,*)'Is the CI calculation J dependant ? (Y/N)'
         READ(IREAD,'(A)') ANS
         IF (ANS.EQ.'N'.OR.ANS.EQ.'n') THEN
            OPEN(UNIT=IREADJ, FILE=NAME(6),STATUS='OLD')
            WRITE(IOUT,*) NAME(6)
            IMCHF=2
         ELSEIF (ANS.EQ.'Y'.OR.ANS.EQ.'y') THEN
            OPEN(UNIT=IREADJ, FILE=NAME(2),STATUS='OLD')
            WRITE(IOUT,*) NAME(2)
            IMCHF=3
         ELSE
            GO TO 35
         ENDIF
         WRITE(ISCW,*) 'Give the number of the CSF in the CI expansion, 
     :for which the hfs. is to be calculated. '
         READ(IREAD,*) NR
      ELSE
         GO TO 30 
      ENDIF
      IF (ICALC.LT.2) THEN
         IF (ICALC.EQ.0) IDIAG=1
         IF (ICALC.EQ.1) IDIAG=0
         INTGR=1
         OPEN(UNIT=IREADW,FILE=NAME(3),STATUS='OLD',FORM='UNFORMATTED')
40       WRITE(ISCW,*) 'Full print-out ? (Y/N) '
         READ(IREAD,'(A)') ANS
         IF (ANS.EQ.'Y'.OR.ANS.EQ.'y') THEN 
            IPRINT=1
            WRITE(ISCW,*) 'Tolerance for printing (in MHz)'
            READ(IREAD,'(F11.6)') PRIMIN
            WRITE(IOUT,'(A,F11.6,A)') ' Tolerance for printing    '
     :      ,PRIMIN,' MHz '  
         ELSEIF (ANS.EQ.'N'.OR.ANS.EQ.'n') THEN
            IPRINT=0
         ELSE
            GO TO 40 
         ENDIF
50       WRITE(ISCW,*) 'Give 2*I and nuclear dipole and quadrupole momen
     :ts (in n.m. and barn) '
         READ(IREAD,*) II,DIPM,Q
         IF (II.EQ.0) GO TO 50 
         IF (DABS(DIPM).LT.EPSILO) THEN
            IF (DABS(Q).LT.EPSILO) GO TO 50 
         ENDIF
         GI=2.D0*DIPM/DFLOAT(II)
         WRITE(IOUT,'(A,F11.6,A)') ' Nuclear dipole moment     ',DIPM,
     :' n.m.'
         WRITE(IOUT,'(A,F11.6,A)') ' Nuclear quadrupole moment ',Q,
     :' barns'
         WRITE(IOUT,'(A,I4)') ' 2*I=',II
      ELSEIF (ICALC.EQ.2) THEN
         INTGR=0
         IPRINT=1
         IDIAG=0
         IMCHF=3
      ELSE
         GO TO 20 
      ENDIF

*     Set factorials

      CALL FACTRL(NFACTS)

*     Start to perform the calculations.
*     Read the configurationlist and print it.
     
      WRITE(ISCW,620)
      CALL CFGN1(NAME(1))

*     Determine L,S,Jmax,Jmin and the LSJ dependent factors for the main 
*     configuration. Save the LSJ dependent factors for later use.

      CALL LSJ(NR,LL,SS,JJMAX,JJMIN)
      CALL LSJ(NR,LL1,SS1,JJ1MAX,JJ1MIN)
      CALL LSJ(NR,LL2,SS2,JJ2MAX,JJ2MIN)
      CALL LSJFACT(NR,NR,NJQ,IDIAG)
      DO 60 K=1,NJQ
         ORBF1(K)=ORBF(K)
         DIPF1(K)=DIPF(K)
         CONTF1(K)=CONTF(K)
         QUADF1(K)=QUADF(K)
60    CONTINUE 

*     Read the J dependent weights of the CSFs and store them in a
*     vector WT(JI,K1). Read the radial wavefunctions and sort
*     them according to the order in IAJCMP.

      IF (INTGR.EQ.1) THEN
         CALL READWT(IMCHF,NCFG,NR,JJMAX,JJMIN,WT)
         CALL READWFN 
      ENDIF

*     Calculate and print the contributions to the A and B factors from
*     every pair of CSFs.

      INCFG=NCFG
      DO 70 JI=1,NCFG
         WRITE(ISCW,*) '   Ja=',JI
         IF (IMCHF.LT.3.OR.IDIAG.EQ.1) INCFG=JI
         DO 80 JF=1,INCFG
            ID=0
            IF ((IMCHF.LT.3.OR.IDIAG.EQ.1).AND.(JI.NE.JF)) ID=1
            NOVLPS=0
            JMUP=0
            JNUP=0
            JMU=0
            JNU=0
            IF (NORTH.NE.0) CALL NORTBP(JI,JF)
            CALL SETUP(JI,JF)

*  Calculate the LSJ dependent factors. In the J independant expansions
*  the factors are identical for every pair of CSFs and should not be 
*  recalculated. In J dependant expansions there is no need for  
*  recalculation when the pair of CSFs has the same LS terms as the 
*  main CSF.

            IF (IMCHF.EQ.3) THEN
               CALL LSJ(JI,LL1,SS1,JJ1MAX,JJ1MIN)
               CALL LSJ(JF,LL2,SS2,JJ2MAX,JJ2MIN)
               IF (LL1.EQ.LL.AND.LL2.EQ.LL.AND.SS1.EQ.SS.AND.SS2.EQ.SS) 
     :            THEN
                  DO 90 K=1,NJQ
                     ORBF(K)=ORBF1(K)
                     DIPF(K)=DIPF1(K)
                     CONTF(K)=CONTF1(K)
                     QUADF(K)=QUADF1(K)
90                CONTINUE
               ELSE
                  CALL LSJFACT(JI,JF,NJQ,IDIAG)
               ENDIF
            ENDIF

*           Multiply WT(JI,K1) and WT(JF,K2) and save it in a vector
*           WTJIJF(K)

            IF (INTGR.EQ.1) CALL MULTWT(JI,JF,IMCHF,IDIAG,WT,WTJIJF) 
            
            IREPRI=0

            ICOPRI=0
            IF (SS1.NE.SS2) GOTO 101 
            CALL ORBITAL(JI,JF)
            DO 100 N=1,NHY
               IDOPRI=0
               IF (DABS(VHY(N)).GT.EPSILO) THEN
                  IF (INTGR.EQ.1) CALL RADIAL1(N,RADINT,OVLINT1,
     :            OVLINT2)
                  DO 110  K=1,NJQ
                     IF (DABS(ORBF(K)).GT.EPSILO) THEN
                        ORB=ORBF(K)*VHY(N)
                        CONTRI=WTJIJF(K)*AMHZ*GI*ORB*RADINT
                        IF ((DABS(CONTRI).GT.PRIMIN.AND.INTGR.EQ.1)
     1                  .AND.(IPRINT.EQ.1))
     2                  CALL PRINT(JI,JF,N,K,WTJIJF(K),1)
                        IF (INTGR.EQ.0)
     1                  CALL PRINT(JI,JF,N,K,0.D0,1)
                        IF (ID.EQ.1) CONTRI=2.D0*CONTRI
                        ORBA(K)=ORBA(K)+CONTRI
                     ENDIF
110               CONTINUE
               ENDIF
100         CONTINUE

101         CONTINUE

            ICOPRI=0
            CALL DIPOLE(JI,JF)
            DO 120 N=1,NHY
               IDOPRI=0
               IF (DABS(VHY(N)).GT.EPSILO) THEN
                  IF (INTGR.EQ.1) CALL RADIAL1(N,RADINT,OVLINT1,
     :            OVLINT2)
                  DO 130 K=1,NJQ
                     IF (DABS(DIPF(K)).GT.EPSILO) THEN
                        DIP=DIPF(K)*VHY(N)
                        CONTRI=WTJIJF(K)*AMHZ*GI*DIP*RADINT
                        IF ((DABS(CONTRI).GT.PRIMIN.AND.INTGR.EQ.1).
     1                  AND.(IPRINT.EQ.1))
     2                  CALL PRINT(JI,JF,N,K,WTJIJF(K),2)
                        IF (INTGR.EQ.0)
     1                  CALL PRINT(JI,JF,N,K,0.D0,2)
                        IF (ID.EQ.1) CONTRI=2.D0*CONTRI
                        DIPA(K)=DIPA(K)+CONTRI
                     ENDIF
130               CONTINUE
               ENDIF
120         CONTINUE

            ICOPRI=0
            IF (LL1.NE.LL2) GOTO 141
            CALL CONTACT(JI,JF)
            DO 140  N=1,NHY
               IDOPRI=0
               IF (DABS(VHY(N)).GT.EPSILO) THEN
                  IF (INTGR.EQ.1) CALL RADIAL2(N,
     :            AZSQR,OVLINT1,OVLINT2)
                  DO 150  K=1,NJQ
                     IF (DABS(CONTF(K)).GT.EPSILO) THEN
                        CONT=CONTF(K)*VHY(N)
                        CONTRI=WTJIJF(K)*AMHZ*GI*CONT*AZSQR
                        IF ((DABS(CONTRI).GT.PRIMIN.AND.INTGR.EQ.1).
     1                  AND.(IPRINT.EQ.1)) 
     2                  CALL PRINT(JI,JF,N,K,WTJIJF(K),3)
                        IF (INTGR.EQ.0)
     1                  CALL PRINT(JI,JF,N,K,0.D0,3)
                        IF (ID.EQ.1) CONTRI=2.D0*CONTRI
                        CONTA(K)=CONTA(K)+CONTRI
                     ENDIF
150               CONTINUE
               ENDIF
140         CONTINUE

141         CONTINUE

            ICOPRI=0
            IF (SS1.NE.SS2) GOTO 161
            CALL QDRPOLE(JI,JF)
            DO 160 N=1,NHY
               IDOPRI=0
               IF (DABS(VHY(N)).GT.EPSILO) THEN
                  IF (INTGR.EQ.1) CALL RADIAL1(N,RADINT,OVLINT1,
     :            OVLINT2)
                  DO 170 K=1,NJQ
                     IF (DABS(QUADF(K)).GT.EPSILO) THEN
                        QUAD=QUADF(K)*VHY(N)
                        CONTRI=WTJIJF(K)*BMHZ*Q*QUAD*RADINT
                        IF ((DABS(CONTRI).GT.PRIMIN.AND.INTGR.EQ.1).
     1                  AND.(IPRINT.EQ.1))  
     2                  CALL PRINT(JI,JF,N,K,WTJIJF(K),4)
                        IF (INTGR.EQ.0)
     1                  CALL PRINT(JI,JF,N,K,0.D0,4)
                        IF (ID.EQ.1) CONTRI=2.D0*CONTRI
                        QUADB(K)=QUADB(K)+CONTRI
                     ENDIF
170               CONTINUE
               ENDIF 
160         CONTINUE

161         CONTINUE

80       CONTINUE
70    CONTINUE

*     Print the final results if radial integrations were performed.

      IF (INTGR.EQ.0) GO TO 999
      WRITE(IOUT,630)
      DO 180 K=1,NJQ
180   WRITE(IOUT,640) J(K),JP(K),ORBA(K),
     :   DIPA(K),CONTA(K),ORBA(K)+DIPA(K)+CONTA(K)
      WRITE(IOUT,650)
      DO 190 K=1,NJQ
190   WRITE(IOUT,660) J(K),JP(K),QUADB(K)
      IF (HFSPARAM.EQ.1) THEN
          IF (ICALC.EQ.1) THEN
             IF (JJMAX.GT.(JJMIN+2)) K=NJQ-2
             IF (JJMAX.EQ.(JJMIN+2)) K=NJQ-1
             IF (JJMAX.EQ.JJMIN) K=NJQ
          ELSE
             K=NJQ
          ENDIF
          BL=DFLOAT(LL)/2.D0
          BS=DFLOAT(SS)/2.D0
          BJ=DFLOAT(JJMAX)/2.D0
          BLJ=(BJ*(BJ+1.D0)+BL*(BL+1.D0)-BS*(BS+1.D0))/2.D0
          BSJ=(BJ*(BJ+1.D0)-BL*(BL+1.D0)+BS*(BS+1.D0))/2.D0
          BSL=(BJ*(BJ+1.D0)-BL*(BL+1.D0)-BS*(BS+1.D0))/2.D0
          IF (ABS(ORBA(K)).LT.EPSILO) THEN
             AL=0.D0
          ELSE
             AL=ORBA(K)*BL*BJ*(BJ+1.D0)/(2.D0*AMHZ*GI*BLJ)
          ENDIF
          IF (ABS(DIPA(K)).LT.EPSILO) THEN
             AD=0.D0
          ELSE
             AD=DIPA(K)*BS*BL*(2.D0*BL-1.D0)*BJ*(BJ+1.D0)/
     :       (AMHZ*GS*GI*(3.D0*BSL*BLJ-BL*(BL+1.D0)*BSJ))
          ENDIF
          IF (ABS(CONTA(K)).LT.EPSILO) THEN
             AC=0.D0
          ELSE
             AC=3.D0*CONTA(K)*BS*BJ*(BJ+1.D0)/(AMHZ*GS*GI*BSJ)
          ENDIF
          IF (ABS(QUADB(K)).LT.EPSILO) THEN
             BQ=0.D0
          ELSE
             BQ=-QUADB(K)*BL*(2.D0*BL-1.D0)*(BJ+1.D0)*(2.D0*BJ+3.D0)/
     :  (BMHZ*Q*(6.D0*BLJ*BLJ-3.D0*BLJ-2.D0*BL*(BL+1.D0)*BJ*(BJ+1.D0)))
          ENDIF 
         WRITE(IOUT,670)
         WRITE(IOUT,680) AL,AD,AC,BQ
      ENDIF
999   WRITE(IOUT,'(30X,A)') '**'
CSUN  RTIME=DTIME(TIMES)
CSUN  WRITE(ISCW,'(//A/A//A/3F10.3//)') ' END OF CASE',' ==========',
CSUN :  '     Real       User      System  Time  (in minutes)',
CSUN : RTIME/60.,TIMES(1)/60., TIMES(2)/60.
      END

*-----------------------------------------------------------------------
*     R E A D W F N
*-----------------------------------------------------------------------
*     This subroutine reads the radial wavefunctions and sorts them
*     according to the order in IAJCMP.
*
      SUBROUTINE READWFN 
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWD=30,NCD=100)
      INTEGER SMAX(NWD)
      DIMENSION SAZ(NWD),SP(NOD,NWD)
      CHARACTER*3 EL(NWD),SEL(NWD),AT*6,TT*6
      COMMON /PARAT/D0,D1,D2,D3,D4,D5,D6,D10,H,H1,NO,ND
      COMMON/ADATA/MNCFG,MNNN,P(NOD,NWD),R(NOD),RR(NOD),R2(NOD),
     :       ZED,W(NCD),AZ(NWD),L(NWD),MAX(NWD)
      COMMON/INFORM/IREADC,IWRITE,IOUT,IREADJ,IREADW,ISCW,IREAD,ISC(4)
      COMMON/STATES/NCFG,MAXORB,IAJCMP(NWD),LJCOMP(NWD),NJCOMP(NWD),
     :NOCCSH(NCD),NELCSH(5,NCD),NOCORB(5,NCD),J1QNRD(9,NCD)

*     Read radial wavefunctions 

      I=1
5     READ(IREADW,END=999) AT,TT,SEL(I),MR,Z,ETI,EKI,SAZ(I),
     :   (SP(JJ,I),JJ=1,MR)

*     Check if the orbital is one of those required in the present run.

      READ(SEL(I),'(A3)')ITEMP
      DO 10 J=1,MAXORB
         IF (IAJCMP(J).EQ.ITEMP) THEN
            DO 15 JJJ = MR,NO
15          SP(JJJ,I) = D0
            SMAX(I) = MR
            I = I+1
            IF (I.GT.(60)) STOP ' Too many electrons: max=(30)'
            GO TO 5   
         ENDIF
10    CONTINUE
      GO TO 5 
999   CONTINUE
      ZED=DBLE(Z)

*     Sort the orbitals according to the order in IAJCMP

      DO 20 II=1,MAXORB
         READ(SEL(II),'(A3)')ITEMP
         DO 25 I= 1,MAXORB
            IF (IAJCMP(I).EQ.ITEMP) THEN
               EL(I)=SEL(II)
               AZ(I)=SAZ(II)
               MAX(I)=SMAX(II)
               DO 30 JJ=1,NO
                  P(JJ,I)=SP(JJ,II)
30             CONTINUE
            ENDIF
25       CONTINUE
20    CONTINUE
      DO 35 I=1,MAXORB
         IF (EL(I).EQ.'   ') STOP 'Radial orbital missing.'
35    CONTINUE

      DO 40 I=1,MAXORB
         IF (EL(I)(1:1) .NE. ' ') THEN
            L(I) = LVAL(EL(I)(2:2))
         ELSE
            L(I) = LVAL(EL(I)(3:3))
         END IF
40    CONTINUE
      RHO=-4.D0
      DO 45 J=1,NO
         R(J)=DEXP(RHO)/ZED
         RR(J)=R(J)*R(J)
         R2(J)=DSQRT(R(J))
         RHO=RHO+H
45    CONTINUE
      RETURN
      END
*
*-----------------------------------------------------------------------
*     B L O C K  D A T A
*-----------------------------------------------------------------------
*
      BLOCK DATA
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /PARAT/D0,D1,D2,D3,D4,D5,D6,D10,H,H1,NO,ND
      DATA D0,D1,D2,D3,D4,D5/0.D0,1.D0,2.D0,3.D0,4.D0,.5D0/
      DATA D6,D10/6.D0,10.D0/
      DATA H,H1/.0625D0,.4166666666666667D-01/
      DATA NO,ND/220,218/
      END
*-----------------------------------------------------------------------
*     Q D R P O L E
*-----------------------------------------------------------------------
*
*     This subroutine calculates the quantities
*     VVSHELL(m)<l(rho)||C2||l(sigma)>
*     and saves them in the vector VHY.
*     The reduced one particle matrix element is calculated as
*     <l||C2||l'>=((-1)**(l+(l+l'+2)/2))*RME(l,l',2)
*     The phase of <l||C2||l'> agrees with the conventions in R. Cowan; 
*     The Theory of Atomic Structure and Spectra. The phase of RME 
*     agrees with the conventions in Fano and Racah; Irreducible 
*     Tensorial Sets.

      SUBROUTINE QDRPOLE(JI,JF)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWD=30,NCD=100,EPSILO=1.D-9)
      COMMON/TENSOR/VVSHELL(50),IIHSH,NNOVLP(50),IIRHO(50),IISIG(50),
     : MMU(50),MMUP(50),NNU(50),NNUP(50)
      COMMON/HYPER/NHY,MHY,IRHY(20),ISHY(20),JRHY(20),JSHY(20),VHY(20),
     1 IMUHY(20),IMUPHY(20),JMUHY(20),JMUPHY(20),
     2 INUHY(20),INUPHY(20),JNUHY(20),JNUPHY(20),NNNOVLP(20)
      COMMON/MEDEFN/IHSH,NJ(10),LJ(10),NOSH1(10),NOSH2(10),J1QN1(19,3),
     : J1QN2(19,3),IJFUL(10)
      COMMON/STATES/NCFG,MAXORB,IAJCMP(NWD),LJCOMP(NWD),NJCOMP(NWD),
     : NOCCSH(NCD),NELCSH(5,NCD),NOCORB(5,NCD),J1QNRD(9,NCD)

      KA=2
      KB=0
      ISPIN=0
      NHY=0
      CALL TENSOR2(KA,KB,ISPIN,IRHO,ISIG)
      DO 10 IS=1,IIHSH
      IF (DABS(VVSHELL(IS)).LE.EPSILO) GO TO 10
         LA=LJ(IIRHO(IS))
         LB=LJ(IISIG(IS))
         IF ((LA+LB).LT.2) GO TO 10
         REDMAT=RME(LA,LB,2)
         IF (REDMAT.LE.EPSILO) GO TO 10
         NHY=NHY+1
         JRHO=IJFUL(IIRHO(IS))
         JSIG=IJFUL(IISIG(IS))
         JRHY(NHY)=IAJCMP(JRHO)
         JSHY(NHY)=IAJCMP(JSIG)
         IRHY(NHY)=JRHO
         ISHY(NHY)=JSIG
         NNNOVLP(NHY)=NNOVLP(IS)
         IF (NNOVLP(IS).GE.1) THEN
            JJMU=IJFUL(MMU(IS))
            JJMUP=IJFUL(MMUP(IS))
            IMUHY(NHY)=JJMU
            IMUPHY(NHY)=JJMUP
            JMUHY(NHY)=IAJCMP(JJMU)
            JMUPHY(NHY)=IAJCMP(JJMUP)
            IF (NNOVLP(IS).GT.1) THEN
               JJNU=IJFUL(NNU(IS))
               JJNUP=IJFUL(NNUP(IS))
               INUHY(NHY)=JJNU
               INUPHY(NHY)=JJNUP
               JNUHY(NHY)=IAJCMP(JJNU)
               JNUPHY(NHY)=IAJCMP(JJNUP)
            ENDIF
         ENDIF 
         NNN=LA+(LA+LB+2)/2
         VHY(NHY)=VVSHELL(IS)*((-1)**NNN)*REDMAT
10    CONTINUE
      RETURN
      END
*-----------------------------------------------------------------------
*     R E A D W T
*-----------------------------------------------------------------------
*
*     This subroutine reads the weights of the CSFs in the wavefunction
*     expansion. The J independant weights are read from <name>.c or
*     <name>.l. and the J dependant from <name>.j
*     In case of a CI calculation NR is the ordernumber of the main CSF
*     component in the CI expansion, and determines where in the file
*     the weights are read. The weights are saved in an array W(K,J)
*     where K specifies the CSF and J defines the J quantum number
*

      SUBROUTINE  READWT(IMCHF,NCFG,NR,JJMAX,JJMIN,W)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NCD=100)
      INTEGER CFGNR
      DIMENSION W(NCD,20)
      CHARACTER*26 LINE
      COMMON/INFORM/IREADC,IWRITE,IOUT,IREADJ,IREADW,ISCW,IREAD,ISC(4)

*     J independent weights. Read from <name>.c

      MAX=JJMAX
      MIN=JJMIN

      IF (IMCHF.EQ.1) THEN
         REWIND(UNIT=IREADC)
         READ(IREADC,'(A26)') LINE
         READ(IREADC,*)
         DO 10 K=1,NCFG
            READ(IREADC,100) W(K,1)
            READ(IREADC,*)
10       CONTINUE
100      FORMAT(T41,F10.7)
         RETURN
      ENDIF

*     Weights from CI calc. If IMCHF=2 read from <name>.l 
*     else if IMCHF=3 read from <name>.j

      J=0
      READ(IREADJ,200) NNCFG
200   FORMAT(T38,I4)
      IF (NNCFG.NE.NCFG) THEN
         WRITE(ISCW,*) 'Unequal number of CSFs in <name>.c and <name>.j(
     :1)'
         STOP
      ENDIF

*     Determine how many lines the array of weights occupy
*     and the number of weights on the last line.

      IF (MOD(NNCFG,7).EQ.0) THEN
         NRLINES=NNCFG/7
         NRLL=0
      ELSE
         NRLINES=NNCFG/7+1
         NRLL=MOD(NNCFG,7)
      ENDIF

*     Find the right array of weights in the file.
*     Firstly the right J quantum number should be found.

20    READ(IREADJ,*)
      READ(IREADJ,*)
      READ(IREADJ,300) JJ,NUMBER
300   FORMAT(T8,I6,T23,I4)

      IF (IMCHF.EQ.2) THEN
         MAX=0
         MIN=0 
      ENDIF
      IF (JJ.GT.MAX) THEN
         DO 30 K=1,NUMBER*(NRLINES+2)
            READ(IREADJ,*)
30       CONTINUE
         GO TO 20
      ENDIF    

      J=J+1
      N=0

*     Secondly the right CSF should be found.

      DO 40 I=1,NUMBER
      READ(IREADJ,*)
      READ(IREADJ,400) CFGNR
      IF (CFGNR.EQ.NR) THEN 
         N=N+1 
         DO 50 K=1,NNCFG/7
            READ(IREADJ,500) (W(L,J),L=7*(K-1)+1,7*K)
50       CONTINUE
         IF (NRLL.NE.0) THEN
            READ(IREADJ,500) (W(L,J),L=NNCFG-NRLL+1,NNCFG)
         ENDIF
      ELSE  
         DO 60 K=1,NRLINES
            READ(IREADJ,*)
60       CONTINUE 
      ENDIF 
40    CONTINUE
400   FORMAT(I6)
500   FORMAT(7F10.7)
6     FORMAT(A,I4,A/A/A/A)
7     FORMAT(A,I4,A/A,I4,A/A/A)
8     FORMAT(A,/I4,A/A/A)
9     FORMAT(A,I4,A/A,I4,A/A/A)
      IF (N.EQ.0) THEN
         IF (IMCHF.EQ.2) THEN
            WRITE(ISCW,6)  'There is no eigenvector in the file <name>.l  
     : that has CSF',NR,' as the','dominating component. Change the numb
     :er of the dominating component in','the wanted eigenvector to the 
     :right one according to your','identification and rerun the program  
     :.'
            STOP
         ELSE 
            WRITE(ISCW,7)  'There is no eigenvector with 2*J=',JJ,' in t
     :he file <name>.j that has','CSF',NR,' as the dominating component.
     : Change the number of the','dominating component in the wanted eig
     :envector to the right one','according to your identification and r
     :erun the program.'
            STOP
         ENDIF
      ENDIF
      IF (N.GT.1) THEN
         IF (IMCHF.EQ.2) THEN
            WRITE(ISCW,8)  'There is more than one eigenvector in the fi
     :le <name>.l that has CSF',NR,' as the dominating component. Discar
     :d the unwanted eigenvectors and','change the digit giving the numb
     :er of eigenvektors to the right one','and rerun the program.'
            STOP
         ELSE
            WRITE(ISCW,9) 'There is more than one eigenvector with 2*J='
     :,JJ,' in the file <name>.j','that has CSF',NR,' as the dominating
     :component. Discard the unwanted','eigenvectors and change the digi
     :t giving the number of eigenvectors','to the right one and rerun t
     :he program.'
            STOP
         ENDIF
      ENDIF
      IF (JJ.EQ.MIN) THEN
         GO TO 999
      ENDIF
      GO TO 20
999   RETURN
      END


*-----------------------------------------------------------------------
*     C O N T A C T
*-----------------------------------------------------------------------
*
*     This subroutine calculates the quantities
*     VVSHELL(m)<s||s1||s>
*     and saves them in the vector VHY

      SUBROUTINE CONTACT(JI,JF)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWD=30,NCD=100,EPSILO=1.D-9)

      COMMON/TENSOR/VVSHELL(50),IIHSH,NNOVLP(50),IIRHO(50),IISIG(50),
     : MMU(50),MMUP(50),NNU(50),NNUP(50)
      COMMON/HYPER/NHY,MHY,IRHY(20),ISHY(20),JRHY(20),JSHY(20),VHY(20),
     1 IMUHY(20),IMUPHY(20),JMUHY(20),JMUPHY(20),
     2 INUHY(20),INUPHY(20),JNUHY(20),JNUPHY(20),NNNOVLP(20)
      COMMON/MEDEFN/IHSH,NJ(10),LJ(10),NOSH1(10),NOSH2(10),J1QN1(19,3),
     : J1QN2(19,3),IJFUL(10)
      COMMON/STATES/NCFG,MAXORB,IAJCMP(NWD),LJCOMP(NWD),NJCOMP(NWD),
     : NOCCSH(NCD),NELCSH(5,NCD),NOCORB(5,NCD),J1QNRD(9,NCD)

      KA=0
      KB=1
      ISPIN=1
      NHY=0
      CALL TENSOR2(KA,KB,ISPIN,IRHO,ISIG)
      DO 10 IS=1,IIHSH
         IF (DABS(VVSHELL(IS)).LE.EPSILO) GO TO 10
         LA=LJ(IIRHO(IS))
         LB=LJ(IISIG(IS))
         IF (LA.NE.0.OR.LB.NE.0) GO TO 10
         NHY=NHY+1
         JRHO=IJFUL(IIRHO(IS))
         JSIG=IJFUL(IISIG(IS))
         JRHY(NHY)=IAJCMP(JRHO)
         JSHY(NHY)=IAJCMP(JSIG)
         IRHY(NHY)=JRHO
         ISHY(NHY)=JSIG
         NNNOVLP(NHY)=NNOVLP(IS)
         IF (NNOVLP(IS).GE.1) THEN
            JJMU=IJFUL(MMU(IS))
            JJMUP=IJFUL(MMUP(IS))
            IMUHY(NHY)=JJMU
            IMUPHY(NHY)=JJMUP
            JMUHY(NHY)=IAJCMP(JJMU)
            JMUPHY(NHY)=IAJCMP(JJMUP)
            IF (NNOVLP(IS).GT.1) THEN
               JJNU=IJFUL(NNU(IS))
               JJNUP=IJFUL(NNUP(IS))
               INUHY(NHY)=JJNU
               INUPHY(NHY)=JJNUP
               JNUHY(NHY)=IAJCMP(JJNU)
               JNUPHY(NHY)=IAJCMP(JJNUP)
            ENDIF
         ENDIF
         VHY(NHY)=VVSHELL(IS)*DSQRT(3.D0/2.D0)
10    CONTINUE
      RETURN
      END
*-----------------------------------------------------------------------
*     L S J F A C T
*-----------------------------------------------------------------------
*
*     This subroutine calculates the LSJ dependent angular factors
*     for the orbital,spindipolar,contact and quadrupol contribution and
*     saves them in vectors ORBF(K),DIPF(K),CONTF(K),QUADF(K).
*     K is a variable to order the different combinations of J and J'.
*     Observe that in the Breit-Pauli approximation, the intermediate
*     coupling wave functions produced by the MCHF_ASP package
*     assumes the S + L coupling order while the decoupling formulas
*     correspond to the L + S coupling order of Cowan; The Theory of
*     Atomic Structure and Spectra. A correcting factor (-1)**(L+S-J)
*     for both the left hand side and the right hand side has to be
*     included.

      SUBROUTINE LSJFACT(JI,JF,NJQ,IDIAG)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NCD=100)
      INTEGER SS1,SS2,ALPHJ(24)

      COMMON/COEF/ORBF(20),DIPF(20),CONTF(20),QUADF(20),J(20),JP(20),
     :JJMAX,JJMIN,JJ1MAX,JJ1MIN,JJ2MAX,JJ2MIN,LL1,SS1,LL2,SS2

*                       DATA STATEMENTS

      DATA ALPHJ/4H  0 ,4H 1/2,4H  1 ,4H 3/2,4H  2 ,4H 5/2,4H  3 ,4H 7/2
     :          ,4H  4 ,4H 9/2,4H  5 ,4H11/2,4H  6 ,4H13/2,4H  7 ,4H15/2
     :          ,4H  8 ,4H17/2,4H  9 ,4H19/2,4H  10,4H21/2,4H  11,4H23/2
     :/
      DATA GS/2.002319D0/

      NJQ=0

*     The phase factor

      NPH2=(LL2+SS2-LL1-SS1)/2
      PHASE=DFLOAT((-1)**NPH2)
      
      DO 10 JJ=JJMIN,JJMAX,2
         IF (JJ.EQ.0) GO TO 10
         JJP=JJ
         NJQ=NJQ+1

*        J'=J

         IF (((JJ.GE.JJ1MIN).AND.(JJ.LE.JJ1MAX)).AND.
     :      ((JJ.GE.JJ2MIN).AND.(JJ.LE.JJ2MAX))) THEN
            CALL GRACAH(LL1,SS1,2,JJ,JJ,LL2,W61)
            CALL NINEJ(LL1,SS1,JJ,LL2,SS2,JJ,4,2,2,S9J)
            CALL GRACAH(LL1,SS1,JJ,2,JJ,SS2,W62)
            CALL GRACAH(LL1,SS1,4,JJ,JJ,LL2,W63)
            SQRT1=DSQRT(DFLOAT((JJ+1)*4)/DFLOAT(JJ*(JJ+2)))
            SQRT2=DSQRT(DFLOAT(JJ*(JJ+1)*(JJ-1))/DFLOAT((JJ+2)*(JJ+3)))
            NPH1=(SS2-SS1)/2
            ORBF(NJQ)=SQRT1*W61*2.D0*PHASE
            DIPF(NJQ)=(-1.D0)*DSQRT(30.D0)*SQRT1*S9J*GS*PHASE
            CONTF(NJQ)=DFLOAT((-1)**(NPH1))*SQRT1*W62*2.D0*GS
     :      *PHASE/3.D0
            QUADF(NJQ)=(-1.D0)*SQRT2*W63*2.D0*PHASE
         ELSE
            ORBF(NJQ)=0.D0
            DIPF(NJQ)=0.D0
            CONTF(NJQ)=0.D0
            QUADF(NJQ)=0.D0
         ENDIF

         J(NJQ)=ALPHJ(JJ+1)
         JP(NJQ)=ALPHJ(JJP+1)

*        J'=J-1

*        Observe that A(LS;J,J-1)=(-1)A(SL;J,J-1) and that the
*        phase factor therefore is PHASE instead of (-1)*PHASE


         IF (JJP.EQ.JJMIN.OR.IDIAG.EQ.1) GO TO 10
         JJP=JJ-2
         NJQ=NJQ+1
         IF (((JJ.GE.JJ1MIN).AND.(JJ.LE.JJ1MAX)).AND.
     :      ((JJP.GE.JJ2MIN).AND.(JJP.LE.JJ2MAX))) THEN
            CALL GRACAH(LL1,SS1,2,JJP,JJ,LL2,W61)
            CALL NINEJ(LL1,SS1,JJ,LL2,SS2,JJP,4,2,2,S9J1)
            CALL GRACAH(LL1,SS1,JJP,2,JJ,SS2,W62)
            CALL GRACAH(LL1,SS1,4,JJP,JJ,LL2,W63)
            SQRT1=DSQRT(DFLOAT(2)/DFLOAT(JJ))
            SQRT2=DSQRT(DFLOAT(JJ*(JJ-2))/DFLOAT((JJ+2)))
            NPH1=(SS2-SS1+2)/2
            ORBF(NJQ)=SQRT1*W61*2.D0*PHASE
            DIPF(NJQ)=(-1.D0)*DSQRT(30.D0)*SQRT1*S9J1*GS*PHASE
            CONTF(NJQ)=DFLOAT((-1)**(NPH1))*SQRT1*W62*2.D0*GS
     :      *PHASE/3.D0
            QUADF(NJQ)=(-1.D0)*SQRT2*W63*PHASE/(DSQRT(2.D0)*2.D0)
         ELSE
             ORBF(NJQ)=0.D0
             DIPF(NJQ)=0.D0
             CONTF(NJQ)=0.D0
             QUADF(NJQ)=0.D0
         ENDIF

         J(NJQ)=ALPHJ(JJ+1)
         JP(NJQ)=ALPHJ(JJP+1)

*        J'=J-2

         IF (JJP.EQ.JJMIN) GO TO 10
         JJP=JJ-4
         NJQ=NJQ+1
         IF (((JJ.GE.JJ1MIN).AND.(JJ.LE.JJ1MAX)).AND.
     :      ((JJP.GE.JJ2MIN).AND.(JJP.LE.JJ2MAX))) THEN
            CALL GRACAH(LL1,SS1,4,JJP,JJ,LL2,W61)
            SQRT2=DSQRT(DFLOAT((JJ-2)*JJ*(JJ-1)))
            QUADF(NJQ)=(-1.D0)*SQRT2*W61*PHASE/8.D0
         ELSE
            ORBF(NJQ)=0.D0
            DIPF(NJQ)=0.D0
            CONTF(NJQ)=0.D0
         ENDIF

40       J(NJQ)=ALPHJ(JJ+1)
         JP(NJQ)=ALPHJ(JJP+1)
10    CONTINUE
      RETURN
      END
*-----------------------------------------------------------------------
*     L S J
*-----------------------------------------------------------------------
*
*     This subroutine determines 2L, 2S, 2Jmax, and 2Jmin for the
*     CSF with the ordernumber NR.

      SUBROUTINE LSJ(NR,LL,SS,JJMAX,JJMIN)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NCD=100)
      INTEGER SS

      COMMON/STATES/NCFG,MAXORB,IAJCMP(NWD),LJCOMP(NWD),NJCOMP(NWD),
     :NOCCSH(NCD),NELCSH(5,NCD),NOCORB(5,NCD),J1QNRD(9,NCD)
      
      MI=2*NOCCSH(NR)-1
      MQ=J1QNRD(MI,NR)/64
      LL=MOD(MQ,64)-1
      SS=MQ/64-1
      JJMIN=IABS(LL-SS)
      JJMAX=LL+SS
      RETURN
      END
*
*     -----------------------------------------------------------------
*        N I N E J
*     -----------------------------------------------------------------
*
*
      SUBROUTINE NINEJ(I1,I2,I3,I4,I5,I6,I7,I8,I9,XD) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      PARAMETER (KFL1=60,KFL2=12)
*
      LOGICAL FAIL,FREE
      COMMON/COUPLE/M,N,J1(KFL1),J2(KFL2,3),J3(KFL2,3),FREE(KFL1)
* 
*     EVALUATES A 9-J SYMBOL.  THE FIRST NINE ELEMENTS OF THE 
*     PARAMETER LIST ARE DOUBLE THE ACTUAL J-VALUES 
* 
      J1(1)=I1+1 
      J1(2)=I2+1 
      J1(3)=I3+1 
      J1(4)=I4+1 
      J1(5)=I5+1 
      J1(6)=I6+1 
      J1(7)=I7+1 
      J1(8)=I8+1 
      J1(9)=I9+1 
* 
      J2(1,1)=1 
      J2(1,2)=2 
      J2(1,3)=3 
      J2(2,1)=4 
      J2(2,2)=5 
      J2(2,3)=6 
      J2(3,1)=3 
      J2(3,2)=6 
      J2(3,3)=9 
* 
      J3(1,1)=1 
      J3(1,2)=4 
      J3(1,3)=7 
      J3(2,1)=2 
      J3(2,2)=5 
      J3(2,3)=8 
      J3(3,1)=7 
      J3(3,2)=8 
      J3(3,3)=9 
* 
      M=9 
      N=4 
      DO 10 I=1,M
         FREE(I) = .FALSE.
10    CONTINUE
* 
      CALL NJGRAF(RECUP,FAIL) 
* 
      XD=(I3+1)*(I6+1)*(I7+1)*(I8+1) 
      XD=RECUP/DSQRT(XD) 
      RETURN 
      END 
*-----------------------------------------------------------------------
*     O R B I T A L
*-----------------------------------------------------------------------
*
*     This subroutine calculates the quantities
*     VVSHELL(m)<l(rho)||l1||l(sigma)> 
*     and saves them in the vector VHY

      SUBROUTINE ORBITAL(JI,JF)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWD=30,NCD=100,EPSILO=1.D-9)

      COMMON/TENSOR/VVSHELL(50),IIHSH,NNOVLP(50),IIRHO(50),IISIG(50),
     : MMU(50),MMUP(50),NNU(50),NNUP(50)
      COMMON/HYPER/NHY,MHY,IRHY(20),ISHY(20),JRHY(20),JSHY(20),VHY(20),
     1 IMUHY(20),IMUPHY(20),JMUHY(20),JMUPHY(20),
     2 INUHY(20),INUPHY(20),JNUHY(20),JNUPHY(20),NNNOVLP(20)
      COMMON/MEDEFN/IHSH,NJ(10),LJ(10),NOSH1(10),NOSH2(10),J1QN1(19,3),
     : J1QN2(19,3),IJFUL(10)
      COMMON/STATES/NCFG,MAXORB,IAJCMP(NWD),LJCOMP(NWD),NJCOMP(NWD),
     : NOCCSH(NCD),NELCSH(5,NCD),NOCORB(5,NCD),J1QNRD(9,NCD)

      KA=1 
      KB=0
      ISPIN=0
      NHY=0
      CALL TENSOR2(KA,KB,ISPIN,IRHO,ISIG)
      IF (IIHSH.EQ.0) RETURN
      DO 10 IS=1,IIHSH
         IF (DABS(VVSHELL(IS)).LE.EPSILO) GO TO 10
         LA=LJ(IIRHO(IS))
         LB=LJ(IISIG(IS))
         IF ((LA+LB).LT.1.OR.LA.NE.LB) GO TO 10
         NHY=NHY+1
         JRHO=IJFUL(IIRHO(IS))
         JSIG=IJFUL(IISIG(IS))
         IRHY(NHY)=JRHO
         ISHY(NHY)=JSIG
         JRHY(NHY)=IAJCMP(JRHO)
         JSHY(NHY)=IAJCMP(JSIG)
         NNNOVLP(NHY)=NNOVLP(IS)
         IF (NNOVLP(IS).GE.1) THEN
            JJMU=IJFUL(MMU(IS))
            JJMUP=IJFUL(MMUP(IS))
            IMUHY(NHY)=JJMU
            IMUPHY(NHY)=JJMUP
            JMUHY(NHY)=IAJCMP(JJMU)
            JMUPHY(NHY)=IAJCMP(JJMUP)
            IF (NNOVLP(IS).GT.1) THEN
               JJNU=IJFUL(NNU(IS))
               JJNUP=IJFUL(NNUP(IS))
               INUHY(NHY)=JJNU
               INUPHY(NHY)=JJNUP
               JNUHY(NHY)=IAJCMP(JJNU)
               JNUPHY(NHY)=IAJCMP(JJNUP)
            ENDIF
         ENDIF
         VHY(NHY)=VVSHELL(IS)*DSQRT(DFLOAT(LA*(LA+1)*(LA+LA+1)))
10     CONTINUE
      RETURN
      END
*
*-----------------------------------------------------------------------
*     P R I N T
*-----------------------------------------------------------------------
*
*     This subroutine prints different parameters obtained from the
*     calculation of the reduced matrix element between CSF JI and JF.

      SUBROUTINE PRINT(JI,JF,N,K,WTJIJF,ICASE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220,NWD=30,NCD=100)
      INTEGER SS1,SS2

      COMMON/PRINTA/ORB,DIP,CONT,QUAD,RADINT,AZSQR,CONTRI,INTGR,
     : IREPRI,ICOPRI,IDOPRI,OVLINT1,OVLINT2
      COMMON/COEF/ORBF(20),DIPF(20),CONTF(20),QUADF(20),J(20),JP(20),
     : JJMAX,JJMIN,JJ1MAX,JJ1MIN,JJ2MAX,JJ2MIN,LL1,SS1,LL2,SS2
      COMMON/TENSOR/VVSHELL(50),IIHSH,NNOVLP(50),IIRHO(50),IISIG(50),
     : MMU(50),MMUP(50),NNU(50),NNUP(50)
      COMMON/HYPER/NHY,MHY,IRHY(20),ISHY(20),JRHY(20),JSHY(20),VHY(20),
     1 IMUHY(20),IMUPHY(20),JMUHY(20),JMUPHY(20),
     2 INUHY(20),INUPHY(20),JNUHY(20),JNUPHY(20),NNNOVLP(20)
      COMMON/INFORM/IREADC,IWRITE,IOUT,IREADJ,IREADW,ISCW,IREAD,ISC(4)
      COMMON/OVRINT/IOVEL(2),IOVER(2),IOVEP(2)

*                       OUTPUT FORMATS


660   FORMAT(/8H (config,I4,11H|hfs|config,I4,2H ))
670   FORMAT(/13H Orbital term)
680   FORMAT(/17H Spin-dipole term)
690   FORMAT(/19H Fermi-contact term)
700   FORMAT(/16H Quadrupole term)

*     Orbital terms

710   FORMAT(15H matrix element,5X,1H<,A3,9H|R**(-3)|,A3,1H>/
     :2X,69HJ    J'   weight    coeff       matrix element values               
     :Aorb(MHz))
720   FORMAT(15H matrix element,5X,1H<,A3,9H|R**(-3)|,A3,1H>
     :,1X,1H<,A3,1H|,A3,3H>**,I2/
     :2X,69HJ    J'   weight    coeff       matrix element values               
     :Aorb(MHz))
730   FORMAT(15H matrix element,
     :5X,1H<,A3,9H|R**(-3)|,A3,1H>,1X,1H<,A3,1H|,A3,3H>**,I2,
     :1X,1H<,A3,1H|,A3,3H>**,I2/
     :2X,69HJ    J'   weight    coeff       matrix element values               
     :Aorb(MHz))

*     Spin-dipole terms

740   FORMAT(15H matrix element,5X,1H<,A3,9H|R**(-3)|,A3,1H>/
     :2X,69HJ    J'   weight    coeff       matrix element values               
     :Adip(MHz))
750   FORMAT(15H matrix element,5X,1H<,A3,9H|R**(-3)|,A3,1H>
     :,1X,1H<,A3,1H|,A3,3H>**,I2/
     :2X,69HJ    J'   weight    coeff       matrix element values               
     :Adip(MHz))
760   FORMAT(15H matrix element,
     :5X,1H<,A3,9H|R**(-3)|,A3,1H>,1X,1H<,A3,1H|,A3,3H>**,I2,
     :1X,1H<,A3,1H|,A3,3H>**,I2/
     :2X,69HJ    J'   weight    coeff       matrix element values               
     :Adip(MHz))

*     Contact terms

770   FORMAT(15H matrix element,5X,3HAZ(,A3,4H)AZ(,A3,1H)/
     :2X,69HJ    J'   weight    coeff       matrix element values      A        
     :cont(MHz))
780   FORMAT(15H matrix element,5X,3HAZ(,A3,4H)AZ(,A3,1H)
     :,1X,1H<,A3,1H|,A3,3H>**,I2/
     :2X,69HJ    J'   weight    coeff       matrix element values      A        
     :cont(MHz))
790   FORMAT(15H matrix element,
     :5X,3HAZ(,A3,4H)AZ(,A3,1H),1X,1H<,A3,1H|,A3,3H>**,I2,
     :1X,1H<,A3,1H|,A3,3H>**,I2/
     :2X,69HJ    J'   weight    coeff       matrix element values      A      
     :cont(MHz))
800   FORMAT(15H matrix element,5X,1H<,A3,9H|R**(-3)|,A3,1H>/
     :2X,69HJ    J'   weight    coeff       matrix element values               
     :   B(MHz))                                                     

*     Quadrupole terms

810   FORMAT(15H matrix element,5X,1H<,A3,9H|R**(-3)|,A3,1H>
     :,1X,1H<,A3,1H|,A3,3H>**,I2/
     :2X,69HJ    J'   weight    coeff       matrix element values               
     :   B(MHz))
820   FORMAT(15H matrix element,
     :5X,1H<,A3,9H|R**(-3)|,A3,1H>,1X,1H<,A3,1H|,A3,3H>**,I2,
     :1X,1H<,A3,1H|,A3,3H>**,I2/
     :2X,69HJ    J'   weight    coeff       matrix element values               
     :   B(MHz))
830   FORMAT( A4,X,A4,X,F9.6,X,F9.6,X,F9.4,17X,F12.6,F12.6)
840   FORMAT( A4,X,A4,X,F9.6,X,F9.6,X,F9.4,X,F7.3,9X,
     :F12.6,F12.6)
850   FORMAT( A4,X,A4,X,F9.6,X,F9.6,X,F9.4,X,F7.3,1X,F7.3,1X,
     :F12.6,F12.6)
860   FORMAT( A4,X,A4,12X,F10.7)

      IF (IREPRI.NE.1) THEN      
         IREPRI=1
         WRITE(IOUT,660) JI,JF
      ENDIF

      IF (ICASE.EQ.1) THEN

*     Orbital term

         IF (ICOPRI.NE.1) THEN     
            WRITE(IOUT,670)
            ICOPRI=1
         ENDIF
         IF (IDOPRI.NE.1) THEN 
            IDOPRI=1
            IF (NNNOVLP(N).EQ.0) THEN
               WRITE(IOUT,710) JRHY(N),JSHY(N)
            ELSEIF (NNNOVLP(N).EQ.1) THEN
               WRITE(IOUT,720) JRHY(N),JSHY(N),JMUHY(N),JMUPHY(N),
     :         IOVEP(1)
            ELSEIF (NNNOVLP(N).GT.1) THEN
               WRITE(IOUT,730) JRHY(N),JSHY(N),JMUHY(N),JMUPHY(N),
     :         IOVEP(1),JNUHY(N),JNUPHY(N),IOVEP(2)
            ENDIF
         ENDIF 
         IF (INTGR.EQ.1) THEN
            IF (NNNOVLP(N).EQ.0) THEN
                WRITE(IOUT,830) J(K),JP(K),WTJIJF
     :         ,ORB,RADINT,CONTRI
            ELSEIF (NNNOVLP(N).EQ.1) THEN
                WRITE(IOUT,840) J(K),JP(K),WTJIJF
     :         ,ORB,RADINT/OVLINT1,OVLINT1,CONTRI
            ELSEIF (NNNOVLP(N).EQ.2) THEN
                WRITE(IOUT,850) J(K),JP(K),WTJIJF
     :         ,ORB,RADINT/OVLINT1*OVLINT2,OVLINT1,OVLINT2,CONTRI
            ENDIF
         ELSE
            WRITE(IOUT,860) J(K),JP(K),ORB
         ENDIF
      ELSEIF (ICASE.EQ.2) THEN
 
*     Spin-dipole term

         IF (ICOPRI.NE.1) THEN     
            WRITE(IOUT,680)
            ICOPRI=1
         ENDIF
         IF (IDOPRI.NE.1) THEN 
            IDOPRI=1
            IF (NNNOVLP(N).EQ.0) THEN
               WRITE(IOUT,740) JRHY(N),JSHY(N)
            ELSEIF (NNNOVLP(N).EQ.1) THEN
               WRITE(IOUT,750) JRHY(N),JSHY(N),JMUHY(N),JMUPHY(N),
     :         IOVEP(1)
            ELSEIF (NNNOVLP(N).GT.1) THEN
               WRITE(IOUT,760) JRHY(N),JSHY(N),JMUHY(N),JMUPHY(N),
     :         IOVEP(1),JNUHY(N),JNUPHY(N),IOVEP(2)
            ENDIF
         ENDIF 
         IF (INTGR.EQ.1) THEN
            IF (NNNOVLP(N).EQ.0) THEN
                WRITE(IOUT,830) J(K),JP(K),WTJIJF
     :         ,DIP,RADINT,CONTRI
            ELSEIF (NNNOVLP(N).EQ.1) THEN
                WRITE(IOUT,840) J(K),JP(K),WTJIJF
     :         ,DIP,RADINT/OVLINT1,OVLINT1,CONTRI
            ELSEIF (NNNOVLP(N).EQ.2) THEN
                WRITE(IOUT,850) J(K),JP(K),WTJIJF
     :         ,DIP,RADINT/OVLINT1*OVLINT2,OVLINT1,OVLINT2,CONTRI
            ENDIF
         ELSE
            WRITE(IOUT,860) J(K),JP(K),DIP
         ENDIF
      ELSEIF (ICASE.EQ.3) THEN

*     Contact term

         IF (ICOPRI.NE.1) THEN     
            WRITE(IOUT,690)
            ICOPRI=1
         ENDIF
         IF (IDOPRI.NE.1) THEN 
            IDOPRI=1
            IF (NNNOVLP(N).EQ.0) THEN
               WRITE(IOUT,770) JRHY(N),JSHY(N)
            ELSEIF (NNNOVLP(N).EQ.1) THEN
               WRITE(IOUT,780) JRHY(N),JSHY(N),JMUHY(N),JMUPHY(N),
     :         IOVEP(1) 
            ELSEIF (NNNOVLP(N).GT.1) THEN
               WRITE(IOUT,790) JRHY(N),JSHY(N),JMUHY(N),JMUPHY(N),
     :         IOVEP(1),JNUHY(N),JNUPHY(N),IOVEP(2)
            ENDIF
         ENDIF 
         IF (INTGR.EQ.1) THEN
            IF (NNNOVLP(N).EQ.0) THEN
               WRITE(IOUT,830) J(K),JP(K),WTJIJF,CONT,AZSQR,CONTRI
            ELSEIF (NNNOVLP(N).EQ.1) THEN
               WRITE(IOUT,840) J(K),JP(K),WTJIJF
     :        ,CONT,AZSQR/OVLINT1,OVLINT1,CONTRI
            ELSEIF (NNNOVLP(N).EQ.2) THEN
               WRITE(IOUT,850) J(K),JP(K),WTJIJF
     :        ,CONT,AZSQR/OVLINT1*OVLINT2,OVLINT1,OVLINT2,CONTRI
            ENDIF
         ELSE
            WRITE(IOUT,860) J(K),JP(K),CONT
         ENDIF
      ELSEIF (ICASE.EQ.4) THEN

*     Quadrupole term

         IF (ICOPRI.NE.1) THEN     
            WRITE(IOUT,700)
            ICOPRI=1
         ENDIF
         IF (IDOPRI.NE.1) THEN     
            IDOPRI=1
            IF (NNNOVLP(N).EQ.0) THEN
               WRITE(IOUT,800) JRHY(N),JSHY(N)
            ELSEIF (NNNOVLP(N).EQ.1) THEN
               WRITE(IOUT,810) JRHY(N),JSHY(N),JMUHY(N),JMUPHY(N),
     :         IOVEP(1)
            ELSEIF (NNNOVLP(N).GT.1) THEN
               WRITE(IOUT,820) JRHY(N),JSHY(N),JMUHY(N),JMUPHY(N),
     :         IOVEP(1),JNUHY(N),JNUPHY(N),IOVEP(2)
            ENDIF
         ENDIF
         IF (INTGR.EQ.1) THEN
            IF (NNNOVLP(N).EQ.0) THEN
                WRITE(IOUT,830) J(K),JP(K),WTJIJF
     :         ,QUAD,RADINT,CONTRI
            ELSEIF (NNNOVLP(N).EQ.1) THEN
                WRITE(IOUT,840) J(K),JP(K),WTJIJF
     :         ,QUAD,RADINT/OVLINT1,OVLINT1,CONTRI
            ELSEIF (NNNOVLP(N).EQ.2) THEN
                WRITE(IOUT,850) J(K),JP(K),WTJIJF
     :         ,QUAD,RADINT/OVLINT1*OVLINT2,OVLINT1,OVLINT2,CONTRI
            ENDIF
         ELSE
            WRITE(IOUT,860) J(K),JP(K),QUAD
         ENDIF
      ENDIF
      RETURN
      END
*
*-----------------------------------------------------------------------
*     Q U A D R 
*-----------------------------------------------------------------------
*
*     quadr integrates P(I)*P(J)*R**KK by Simpson's rule 
*

      DOUBLE PRECISION FUNCTION QUADR(I,J,KK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWD=30,NCD=100)
      COMMON /PARAT/D0,D1,D2,D3,D4,D5,D6,D10,H,H1,NO,ND
      COMMON/ADATA/MNCFG,MNNN,P(NOD,NWD),R(NOD),RR(NOD),R2(NOD),
     :       ZED,W(NCD),AZ(NWD),L(NWD),MAX(NWD)
      
      Z=ZED
      K = KK + 2
      LI = L(I)
      LJ = L(J)
      DEN = LI + LJ + 1 + K
      ZR = Z*R(4)
      BI = (P(4,I)/(AZ(I)*R2(4)*R(4)**LI) - D1+ZR/(LI+1))/ZR**2
      BJ = (P(4,J)/(AZ(J)*R2(4)*R(4)**LJ) - D1+ZR/(LJ+1))/ZR**2
      ALPHA= (D1/(LI + 1) + D1/(LJ + 1))/(DEN + D1)
      ZR = Z*R(1)
      BETA = (DEN+D1)*ALPHA**2 - D2*(BI+BJ+D1/((LI+1)*(LJ+1)))/(DEN+D2)
      D = P(1,I)*P(1,J)*R(1)**K*(((BETA*ZR+ALPHA)*ZR+D1)/(DEN*H1)+D5)
      M = MIN0(MAX(I),MAX(J)) - 1
      DO 1 JJ = 2,M,2
      JP = JJ + 1
    1 D=D+D2*P(JJ,I)*P(JJ,J)*R(JJ)**K+P(JP,I)*P(JP,J)*R(JP)**K
      QUADR = D*H1
      RETURN
      END
*-----------------------------------------------------------------------
*        D I P O L E
*-----------------------------------------------------------------------
*
*    This subroutine calculates the quantities
*    VSHELL(m)<l(rho)||C2||l(sigma)><s||s1||s>
*    and saves them in the vector VHY
*    The reduced one particle matrix element is calculated as
*    <l||C2||l'>=((-1)**(l+(l+l'+2)/2))*RME(l,l',2)
*    The phase of <l||C2||l'> agrees with the conventions of R. Cowan;
*    The Theory of Atomic Structure and Spectra.The phase of RME(l,l',2)
*    agrees with the convention of Fano and Racah; Irreducible Tensorial
*    Sets.

      SUBROUTINE DIPOLE(JI,JF)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWD=30,NCD=100,EPSILO=1.D-9)

      COMMON/TENSOR/VVSHELL(50),IIHSH,NNOVLP(50),IIRHO(50),IISIG(50),
     : MMU(50),MMUP(50),NNU(50),NNUP(50)
      COMMON/HYPER/NHY,MHY,IRHY(20),ISHY(20),JRHY(20),JSHY(20),VHY(20),
     1 IMUHY(20),IMUPHY(20),JMUHY(20),JMUPHY(20),
     2 INUHY(20),INUPHY(20),JNUHY(20),JNUPHY(20),NNNOVLP(20)
      COMMON/MEDEFN/IHSH,NJ(10),LJ(10),NOSH1(10),NOSH2(10),J1QN1(19,3),
     : J1QN2(19,3),IJFUL(10)
      COMMON/STATES/NCFG,MAXORB,IAJCMP(NWD),LJCOMP(NWD),NJCOMP(NWD),
     : NOCCSH(NCD),NELCSH(5,NCD),NOCORB(5,NCD),J1QNRD(9,NCD)

      KA=2
      KB=1
      ISPIN=2
      NHY=0
      CALL TENSOR2(KA,KB,ISPIN,IRHO,ISIG)
      IF (IIHSH.EQ.0) RETURN
      DO 10 IS=1,IIHSH
         IF (DABS(VVSHELL(IS)).LE.EPSILO) GO TO 10
         LA=LJ(IIRHO(IS))
         LB=LJ(IISIG(IS))
         IF ((LA+LB).LT.2) GO TO 10
         REDMAT=RME(LA,LB,2)
         IF (REDMAT.LE.EPSILO) GO TO 10
         NHY=NHY+1
         JRHO=IJFUL(IIRHO(IS))
         JSIG=IJFUL(IISIG(IS))
         IRHY(NHY)=JRHO
         ISHY(NHY)=JSIG
         JRHY(NHY)=IAJCMP(JRHO)
         JSHY(NHY)=IAJCMP(JSIG)
         NNNOVLP(NHY)=NNOVLP(IS)
         IF (NNOVLP(IS).GE.1) THEN
            JJMU=IJFUL(MMU(IS))
            JJMUP=IJFUL(MMUP(IS))
            IMUHY(NHY)=JJMU
            IMUPHY(NHY)=JJMUP
            JMUHY(NHY)=IAJCMP(JJMU)
            JMUPHY(NHY)=IAJCMP(JJMUP)
            IF (NNOVLP(IS).GT.1) THEN
               JJNU=IJFUL(NNU(IS))
               JJNUP=IJFUL(NNUP(IS))
               INUHY(NHY)=JJNU
               INUPHY(NHY)=JJNUP
               JNUHY(NHY)=IAJCMP(JJNU)
               JNUPHY(NHY)=IAJCMP(JJNUP)
            ENDIF
         ENDIF
         NNN=LA+(LA+LB+2)/2
         VHY(NHY)=VVSHELL(IS)*((-1)**NNN)*REDMAT*DSQRT(3.D0/2.D0)
10    CONTINUE
      RETURN
      END
*-----------------------------------------------------------------------
*     M U L T W T
*-----------------------------------------------------------------------
*
*     The J dependent weights for CSF JI and JF are  multiplied and 
*     saved in WTJIJF(NJQ). NJQ is a parameter that orders the 
*     different J and J' combinations. This order agrees with the order 
*     in LSJFACT.

      SUBROUTINE MULTWT(JI,JF,IMCHF,IDIAG,WT,WTJIJF)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NCD=100)
      DIMENSION WT(NCD,20),WTJIJF(20)
      INTEGER SS1,SS2

      COMMON/COEF/ORBF(20),DIPF(20),CONTF(20),QUADF(20),J(20),JP(20),
     :JJMAX,JJMIN,JJ1MAX,JJ1MIN,JJ2MAX,JJ2MIN,LL1,SS1,LL2,SS2

      NJQ=0

      DO 10 JJ=JJMIN,JJMAX,2
         IF (JJ.EQ.0) GO TO 10

*        J'=J

         JJP=JJ
         NJQ=NJQ+1
         IF (((JJ.GE.JJ1MIN).AND.(JJ.LE.JJ1MAX)).AND.
     :      ((JJ.GE.JJ2MIN).AND.(JJ.LE.JJ2MAX))) THEN
            IF (IMCHF.LT.3) THEN
               K1=1
            ELSE
               K1=(JJMAX-JJ)/2+1
            ENDIF
            WTJIJF(NJQ)=WT(JI,K1)*WT(JF,K1)
         ELSE
            WTJIJF(NJQ)=0.D0
         ENDIF

*        J'=J-1

         IF (JJP.EQ.JJMIN.OR.IDIAG.EQ.1) GO TO 10
         JJP=JJ-2
         NJQ=NJQ+1
         IF (((JJ.GE.JJ1MIN).AND.(JJ.LE.JJ1MAX)).AND.
     :      ((JJP.GE.JJ2MIN).AND.(JJP.LE.JJ2MAX))) THEN
            IF (IMCHF.LT.3) THEN
               K1=1
               K2=1
            ELSE
               K1=(JJMAX-JJ)/2+1
               K2=(JJMAX-JJP)/2+1
            ENDIF
            WTJIJF(NJQ)=WT(JI,K1)*WT(JF,K2)
         ELSE
            WTJIJF(NJQ)=0.D0
         ENDIF

*        J'=J-2

         IF (JJP.EQ.JJMIN) GO TO 10
         JJP=JJ-4
         NJQ=NJQ+1
         IF (((JJ.GE.JJ1MIN).AND.(JJ.LE.JJ1MAX)).AND.
     :      ((JJP.GE.JJ2MIN).AND.(JJP.LE.JJ2MAX))) THEN
            IF (IMCHF.LT.3) THEN
               K1=1
               K2=1
            ELSE
               K1=(JJMAX-JJ)/2+1
               K2=(JJMAX-JJP)/2+1
            ENDIF
            WTJIJF(NJQ)=WT(JI,K1)*WT(JF,K2)
         ELSE
            WTJIJF(NJQ)=0.D0
         ENDIF
10    CONTINUE
      RETURN
      END

*
*-----------------------------------------------------------------------
*     B L O C K    D A T A
*-----------------------------------------------------------------------
*
      BLOCK DATA RCONST
*
      IMPLICIT REAL *8(A-H,O-Z)
*
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
*
*     SET GLOBAL REAL CONSTANTS
*
      DATA ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS/
     :     0.0D 00,0.1D 00,0.5D 00,1.0D 00,2.0D 00,3.0D 00,4.0D 00,
     :     7.0D 00,1.1D 01,1.0D-08/
*
      END
*
*-----------------------------------------------------------------------
*        T E N S O R 2
*-----------------------------------------------------------------------
*
      SUBROUTINE TENSOR2(KA,KB,ISPIN,IRHO,ISIG)
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(KFL1=60,KFL2=12)
      PARAMETER (NWD=30,NCD=100,NORD=NWD*(NWD-1)/2)
*
*     W. D. ROBB   -   NOVEMBER 1971
*     Modified by M. GODEFROID   --  1981
*     Modified by P. JONSSON     --  1991
* **********************************************************************
*
*     A ROUTINE FOR THE EVALUATION OF ANGULAR AND SPIN FACTORS IN THE
*     REDUCED MATRIX ELEMENT OF ANY ONE-ELECTRON TENSOR OPERATOR BETWEEN
*     ARBITRARILY COUPLED L-S CONFIGURATIONS (Modified for non-
*     orthogonal orbitals)
*
************************************************************************
*
*     **  NOTE THAT THE DEFINITIONS OF TENSOR OPERATORS USED ARE THOSE
*     OF FANO AND RACAH, IRREDUCIBLE TENSORIAL SETS, ACADEMIC PRESS 1959
*
************************************************************************
*
*     MODIFIED TO FIT THE STRUCTURE OF THE HFS PROGRAM - 1991
*
************************************************************************

*                       DIMENSION STATEMENTS
*
      DIMENSION J2STO(KFL2,3),J3STO(KFL2,3),NBAR(KFL2),
     :          JBAR(KFL2,3),JPBAR(KFL2,3)
  

*
*                       COMMON BLOCKS
*
      LOGICAL FAIL,FREE

      COMMON/TENSOR/VVSHELL(50),IIHSH,NNOVLP(50),IIRHO(50),IISIG(50),
     : MMU(50),MMUP(50),NNU(50),NNUP(50)
      COMMON/COUPLE/NJ1S,NJ23S,J1(KFL1),J2(KFL2,3),J3(KFL2,3),FREE(KFL1)
      COMMON/INFORM/IREADC,IWRITE,IOUT,IREADJ,IREADW,ISCW,IREAD,ISC(4)
      COMMON/KRON/IDEL(10,10)
      COMMON/MEDEFN/IHSH,NJ(10),LJ(10),NOSH1(10),NOSH2(10),J1QN1(19,3)
     :    ,J1QN2(19,3),IJFUL(10)
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/TERMS/NROWS,ITAB(24),JTAB(24),NTAB(333)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
     1 ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
     2 IORTH(NORD)
      COMMON/OVRINT/IOVEL(2),IOVER(2),IOVEP(2)
      COMMON/MACHOR/IMATCH(2)
      COMMON/FACT/GAM(100)
      COMMON/EMS/IEM(4),JI,JF,IFL
*
  203 FORMAT(//7H NJ,LJ ,10(I6,I3))
  204 FORMAT(//6H NOSH ,10I4)
  205 FORMAT(//6H J1QN ,30I3)
  207 FORMAT(8F15.8)
  208 FORMAT(// 23H PARENT TERMS NOT FOUND//)
  209 FORMAT(//3H J1)
  210 FORMAT(24I5)
  211 FORMAT(24H J2                   J3)
  212 FORMAT(3I5,I10,2I5)
  213 FORMAT(///26H ORBITAL RECOUPLING COEFF //)
  214 FORMAT(///23H SPIN RECOUPLING COEFF //)
  215 FORMAT(/29H THE CONTRIBUTION FROM SHELLS,2I2,3H IS,F15.8,2X,F15.8)
  216  FORMAT(//21H THIS IS NOT A PARENT)
  220 FORMAT(//6H SHELL,I2)
  221 FORMAT(//)
  302 FORMAT(/5X,89H NO CONTRIBUTION FROM TENSOR SINCE MORE THAN ONE ELE
     1CTRON DIFFERENT IN THE CONFIGURATIONS/)
  303 FORMAT(/5X,114H NO CONTRIBUTION FROM TENSOR SINCE THE TRIANGLE REL
     1ATION BETWEEN KA AND THE TOTAL ANGULAR MOMENTA IS NOT SATISFIED/)
  313 FORMAT(19H SPECTATOR SUBSHELL,I3,69H HAS DIFFERENT QUANTUM NUMBERS
     1 ON THE TWO SIDES OF THE MATRIX ELEMENT/)
  601 FORMAT(//5H KA =,I2,
     1        /5H KB =,I2,
     2        /8H ISPIN =,I2)
*
      IF(KA.EQ.2.AND.KB.EQ.0) IFL=4
      IF(KA.EQ.2.AND.KB.EQ.1) IFL=2
      IF(KA.EQ.0.AND.KB.EQ.1) IFL=3
      IF(KA.EQ.1.AND.KB.EQ.0) IFL=1

      DO 7836 I=1,20
         VVSHELL(I)=0.D0
7836  CONTINUE
      IIHSH=0




      IF(NOVLPS.EQ.0) THEN
         IOVEP(1)=0
         IOVEP(2)=0
      ELSE IF(NOVLPS.EQ.1) THEN
         IOVEP(2)=0
      ENDIF
      AJF=ONE
      RML=ZERO
      RML2=ZERO
      IONE=1
      IZERO=0
      NOVSTO=NOVLPS
      NTOT=0
      ISUM=0
      MUSTO=0
      NUSTO=0
      MUPSTO=0
      NUPSTO=0
      NOVSTO=0
      LMUSTO=0
      LNUSTO=0
      LMUPST=0
      LNUPST=0
      IF(NBUG6.NE.0) WRITE(IWRITE,601) KA,KB,ISPIN,JI,JF
*     DO 100 IS=1,IHSH
      VSHELL=ZERO
* 100 CONTINUE
      IHSHP1=IHSH+1
      I2HSH=IHSH*2-1
*
*     PRINT OUT THE OCCUPATION AND COUPLING ARRAYS
*
      IF(NBUG6.NE.1) GO TO 103
    2 WRITE(IWRITE,203) (NJ(I),LJ(I),I=1,IHSH)
      WRITE(IWRITE,204)(NOSH1(J),J=1,IHSH)
      WRITE(IWRITE,204)(NOSH2(J),J=1,IHSH)
      WRITE(IWRITE,205) ((J1QN1(J,K),K=1,3),J=1,I2HSH)
      WRITE(IWRITE,205) ((J1QN2(J,K),K=1,3),J=1,I2HSH)
*
*      TEST FOR AT MOST ONE ELECTRON DIFFERENCE IN CONFIGURATIONS
*
*
*     TEST FOR TRIANGLE RELATION BETWEEN KA AND TOTAL ANGULAR MOMENTA
*
  103 IF(ISPIN.EQ.0) GO TO 198
      K=3
      KC=KB
      IF(ISPIN.EQ.2) GO TO 199
      IF(J1QN1(I2HSH,2).NE.J1QN2(I2HSH,2)) GO TO 183
      GO TO 199
  198 K=2
      KC=KA
      IF(J1QN1(I2HSH,3).NE.J1QN2(I2HSH,3)) GO TO 183
  199 LB=J1QN1(I2HSH,K)-1
      NB=J1QN2(I2HSH,K)-1
      MB=KC+KC
      BTST=TRITST(MB,LB,NB)
      IF(DABS(BTST).GT.EPS) GO TO 301
      IF(K.EQ.2.OR.ISPIN.LT.2) GO TO 104
      K=2
      KC=KA
      GO TO 199
*
*      DETERMINE IRHO AND ISIGMA, THE NUMBERS OF THE OCCUPIED SHELLS
*
  104 IBACK3=0
      IF (NBUG6.NE.0) WRITE(IWRITE,1000)
 1000 FORMAT(11H GOT TO 104)
      IX=0
      IY=0
      DO 102 I=1,IHSH
      IF(I.EQ.MU.OR.I.EQ.NU.OR.I.EQ.MUP.OR.I.EQ.NUP) GO TO 102
      N=NOSH1(I)-NOSH2(I)
      IF(IABS(N).GT.1) GO TO 300
      IF(N.EQ.0) THEN
         GO TO 102
      ELSE IF(N.EQ.1) THEN
         IRHO=I
         IX=IX+1
      ELSE
         ISIG=I
         IY=IY+1
      ENDIF
  102 CONTINUE
      IF(IX.GT.1.OR.IY.GT.1) GO TO 300
  101 IF(NOVLPS.NE.0) THEN
         MUSTO=MU
         NUSTO=NU
         MUPSTO=MUP
         NUPSTO=NUP
         NOVSTO=NOVLPS
         LMUSTO=LMU
         LNUSTO=LNU
         LMUPST=LMUP
         LNUPST=LNUP
      ENDIF
*
*     INTERACTING SUBSHELLS FOUND
*
      IF(IX.EQ.1.AND.IY.EQ.1) THEN
         IF(NOVLPS.EQ.0) GO TO 108
         MGO=MATCH(IZERO,IZERO,IZERO,IZERO)
         IF(MGO.EQ.2) GO TO 108
         IRHO=0
         RETURN
      ENDIF
*
*     INTERACTING SUBSHELL ON R.H.S. FOUND. THE ONE ON L.H.S. MUST
*     BE EITHER MU OR NU
*
  506 IF(IX.EQ.0.AND.IY.EQ.1) THEN
         IF(IBACK3.EQ.0) THEN
            IRHO=MU
            IF(MU.NE.NU) IBACK3=1
            MGO=MATCH(IONE,IZERO,IZERO,IZERO)
            IF(MGO.EQ.0) THEN
               IRHO=0
               RETURN
            ELSEIF(MGO.EQ.1) THEN
               IF(IBACK3.EQ.1) THEN
                  MU=MUSTO
                  NU=NUSTO
                  MUP=MUPSTO
                  NUP=NUPSTO
                  NOVLPS=NOVSTO
                  LMU=LMUSTO
                  LNU=LNUSTO
                  LMUP=LMUPST
                  LNUP=LNUPST
                  GO TO 105
               ELSE
                  IRHO=0
                  RETURN
               ENDIF
            ELSE
               GO TO 108
            ENDIF
         ELSE
  105       IRHO=NU
            RML=ZERO
            RML2=ZERO
            IBACK3=0
            MGO=MATCH(IZERO,IONE,IZERO,IZERO)
            IF(MGO.EQ.2) GO TO 108
            IRHO=0
            RETURN
         ENDIF
*
*     INTERACTING SUBSHELL ON L.H.S. FOUND. THE ONE ON R.H.S. MUST
*     BE EITHER MUP OR NUP
*
      ELSE IF(IX.EQ.1.AND.IY.EQ.0) THEN
         IF(IBACK3.EQ.0) THEN
            ISIG=MUP
            IF(MUP.NE.NUP) IBACK3=1
            MGO=MATCH(IZERO,IZERO,IONE,IZERO)
            IF(MGO.EQ.0) THEN
               IRHO=0
               RETURN
            ELSEIF(MGO.EQ.1) THEN
               IF(IBACK3.EQ.1) THEN
                  MU=MUSTO
                  NU=NUSTO
                  MUP=MUPSTO
                  NUP=NUPSTO
                  NOVLPS=NOVSTO
                  LMU=LMUSTO
                  LNU=LNUSTO
                  LMUP=LMUPST
                  LNUP=LNUPST
                  GO TO 106
               ELSE
                  IRHO=0
                  RETURN
               ENDIF
            ELSE
               GO TO 108
            ENDIF
         ELSE
  106       ISIG=NUP
            RML=ZERO
            RML2=ZERO
            IBACK3=0
            MGO=MATCH(IZERO,IZERO,IZERO,IONE)
            IF(MGO.LT.2) RETURN
            GO TO 108
         ENDIF
      ELSE
         ISUM=1
         GO TO 107
      ENDIF
*
*     SUM OVER SHELLS
*
  107 IRHO=1
      ISIG=1
 1509 IF (NOSH1(IRHO).EQ.0) GO TO 189
      IF(NOVLPS.EQ.0) GO TO 108
      IF(NOVLPS.EQ.1) GO TO 508
  509 IF(IRHO.NE.MU.AND.IRHO.NE.NU) THEN
         MGO=MATCH(IZERO,IZERO,IZERO,IZERO)
         IF(MGO.EQ.0) THEN
            IRHO=0
            RETURN
         ELSEIF(MGO.EQ.1) THEN
            GO TO 189
         ELSE
            GO TO 108
         ENDIF
      ELSE
         IF(IBACK3.EQ.0) THEN
            ISIG=MUP
            IF(MUP.NE.NUP) IBACK3=1
            IF(IRHO.EQ.MU) THEN
               MGO=MATCH(IONE,IZERO,IONE,IZERO)
            ELSE
               MGO=MATCH(IZERO,IONE,IONE,IZERO)
            ENDIF
            IF(MGO.EQ.0) THEN
               IRHO=0
               RETURN
            ELSEIF(MGO.EQ.1) THEN
               IF(IBACK3.EQ.1) THEN
                  MU=MUSTO
                  NU=NUSTO
                  MUP=MUPSTO
                  NUP=NUPSTO
                  NOVLPS=NOVSTO
                  LMU=LMUSTO
                  LNU=LNUSTO
                  LMUP=LMUPST
                  LNUP=LNUPST
                  GO TO 507
               ELSE
                  GO TO 189
               ENDIF
            ELSE
               GO TO 108
            ENDIF
         ELSE
  507       ISIG=NUP
            IF(NBUG6.NE.0) WRITE(IWRITE,1215) IRHO,ISIG
 1215 FORMAT(15H AT 507, IRHO =,I3,7H ISIG =,I3)
            IBACK3=0
            RML=ZERO
            RML2=ZERO
            IF(IRHO.EQ.MU) THEN
               MGO=MATCH(IONE,IZERO,IZERO,IONE)
            ELSE
               MGO=MATCH(IZERO,IONE,IZERO,IONE)
            ENDIF
            IF(MGO.EQ.0) THEN
               IRHO=0
               RETURN
            ELSEIF(MGO.EQ.1) THEN
               GO TO 189
            ELSE
               GO TO 108
            ENDIF
         ENDIF
      ENDIF
  508 IF(IRHO.NE.MU) THEN
         MGO=MATCH(IZERO,IZERO,IZERO,IZERO)
      ELSE
         ISIG=MUP
         MGO=MATCH(IONE,IZERO,IONE,IZERO)
      ENDIF
      IF(MGO.EQ.0) THEN
         IRHO=0
         RETURN
      ELSEIF(MGO.EQ.1) THEN
         GO TO 189
      ELSE
         GO TO 108
      ENDIF
  108 MEMR=IRHO
      IF (NBUG6.NE.0) WRITE(IWRITE,1001) IRHO,ISIG
 1001 FORMAT(' AT 108, IRHO =',I3,2X,'ISIG =',I3)
      IF (NBUG6.NE.0) WRITE(IWRITE,1010) MU,MUP,NU,NUP,NOVLPS
 1010 FORMAT(' AT 108, MU =',I3,2X,'MUP =',I3,/,
     1       '         NU =',I3,2X,'NUP =',I3,/,
     2       '         NOVLPS =',I3)
*
*     THE BEGINNING OF THE LOOP OVER ALL SHELLS
*
  109 IF(ISUM.EQ.0) GO TO 309
      IF(NBUG6-1) 309,4,309
    4 WRITE(IWRITE,220) IRHO
  309 NTOT=NTOT+1
      LRHO=LJ(IRHO)
      LSIG=LJ(ISIG)
      L1=LRHO+1
      L2=LSIG+1
      AJF=DFLOAT(J1QN1(I2HSH,2))/DFLOAT(2*LRHO+1)
      IF(ISPIN.EQ.1) AJF=DFLOAT(J1QN1(I2HSH,3))*HALF
      IF(ISPIN.EQ.2) AJF=AJF*DFLOAT(J1QN1(I2HSH,3))*HALF
*
*     CHECK THE DIAGONAL CHARACTER OF QUANTUM NUMBERS OF SPECTATOR
*     SHELLS
*
      DO 255 J=1,IHSH
      IF(J.NE.IRHO) THEN
         DO 253 KK=1,3
         JBAR(J,KK)=J1QN1(J,KK)
  253    CONTINUE
      ENDIF
      IF(J.NE.ISIG) THEN
         DO 254 KK=1,3
         JPBAR(J,KK)=J1QN2(J,KK)
  254    CONTINUE
      ENDIF
      IF(J.EQ.MU.OR.J.EQ.NU.OR.J.EQ.MUP.OR.J.EQ.NUP) GO TO 255
      IF(J.EQ.IRHO.OR.J.EQ.ISIG) GO TO 255
      DO 256 KK=1,3
      IF(JBAR(J,KK).NE.JPBAR(J,KK)) GO TO 257
  256 CONTINUE
  255 CONTINUE
      DO 405 J = 1,IHSH
      NBAR(J) = NOSH1(J)-IDEL(J,IRHO)
  405 CONTINUE
      IF(MUP.NE.0) NBAR(MUP)=NOSH2(MUP)-IDEL(MUP,ISIG)
      IF(NUP.NE.0) NBAR(NUP)=NOSH2(NUP)-IDEL(NUP,ISIG)
      IF(NOVLPS.LT.2) GO TO 175
      IF(MU.EQ.NU) THEN
         NKAP=NBAR(MUP)*NBAR(NUP)
      ELSE
         NKAP=NBAR(MU)*NBAR(NU)
      ENDIF
  175 IDELP=2
      IF(IRHO.EQ.IHSH) GO TO 177
      JRHO=IRHO+1
      DO 178 J=JRHO,IHSH
      IF(J.EQ.MUP.OR.J.EQ.NUP) GO TO 178
      IDELP=IDELP+NBAR(J)
  178 CONTINUE
  177 IF(ISIG.EQ.IHSH) GO TO 481
      JSIG=ISIG+1
      DO 180 J=JSIG,IHSH
      IF(J.EQ.MU.OR.J.EQ.NU) GO TO 180
      IDELP=IDELP+NBAR(J)
  180 CONTINUE
  481 IF(NOVLPS.EQ.0) GO TO 181
      IF(NOVLPS.EQ.1) GO TO 925
      IF (MU.EQ.NU) GO TO 904
      IF (MUP.EQ.NUP) GO TO 905
      IF (NBAR(MU).EQ.NBAR(NUP)) GO TO 906
      MU1N = MIN0(MU,MUP) + 1
      MU1X = MAX0(MU,MUP) - 1
      MU2N = MIN0(NU,NUP) + 1
      MU2X = MAX0(NU,NUP) - 1
      MU1=MU
      NU1=NU
      GO TO 933
  906 MU1N = MIN0(MU,NUP) + 1
      MU1X = MAX0(MU,NUP) - 1
      MU2N = MIN0(NU,MUP) + 1
      MU2X = MAX0(NU,MUP) - 1
      MU1=MU
      NU1=NU
      GO TO 933
  905 MU1N = MIN0(MU,MUP) + 1
      MU1X = MAX0(MU,MUP) - 1
      MU2N = MIN0(MUP,NU) + 1
      MU2X = MAX0(MUP,NU) - 1
      MU1=MU
      NU1=NU
      GO TO 933
  904 MU1N = MIN0(MU,MUP) + 1
      MU1X = MAX0(MU,MUP) - 1
      MU2N = MIN0(MU,NUP) + 1
      MU2X = MAX0(MU,NUP) - 1
      MU1=MUP
      NU1=NUP
  933 IF (MU1N .GT. MU1X) GO TO 934
      DO 962 J = MU1N,MU1X
      IF (J.EQ.MU .OR. J.EQ.NU .OR. J.EQ.MUP .OR. J.EQ.NUP) GO TO 962
      IDELP = IDELP + NBAR(J)*NBAR(MU1)
  962 CONTINUE
  934 IF (MU2N .GT. MU2X) GO TO 181
      DO 963 J = MU2N,MU2X
      IF (J.EQ.MU .OR. J.EQ.NU .OR. J.EQ.MUP .OR. J.EQ.NUP) GO TO 963
      IDELP = IDELP + NBAR(J)*NBAR(NU1)
  963 CONTINUE
      GO TO 181
  925 MUMIN1=MIN0(MU,MUP)+1
      MUMAX1=MAX0(MU,MUP)-1
      IF(MUMIN1.GT.MUMAX1) GO TO 181
      DO 927 J=MUMIN1,MUMAX1
      IDELP=IDELP+NBAR(J)*NBAR(MU)
  927 CONTINUE
  181 MINUS=(-1)**IDELP
      FACTR=ONE
      IF(NOVLPS.LT.2) GO TO 127
      IF(MU.NE.NU.AND.MUP.NE.NUP) GO TO 127
*
*     MU=NU OR MUP=NUP
*
      N1=IOVEP(1)
      N2=IOVEP(2)
      N3=N1+N2
      FACTR=FACTR*DEXP((GAM(N3+1)-GAM(N1+1)-GAM(N2+1))/2)
      GO TO 127
  257 IF(NBUG6.EQ.1) WRITE(IWRITE,313) J
      GO TO 189
*
*     SET  J2  AND  J3 .  SAME FOR  L  AND  S
*
  127 M1=IHSH-2
      M2=2*M1+1
      M3=3*IHSH-1
      M4=M3+1
      M5=M3+2
      M10=M5+1
      MN1=M10+1
      NJ1S=MN1
      J2(1,1)=M10
      J2(1,2)=MN1
      J2(1,3)=M5
      J2(2,1)=IRHO
      J2(2,2)=M5
      J2(2,3)=M3
      J3(1,1)=ISIG
      J3(1,2)=M10
      J3(1,3)=M4
      IF(IRHO-1) 128,129,128
  129 J2(3,1)=M3
      GO TO 130
  128 J2(3,1)=1
  130 IF(IRHO-2) 131,132,131
  132 J2(3,2)=M3
      GO TO 133
  131 J2(3,2)=2
  133 J2(3,3)=IHSHP1
      IF(ISIG-1) 134,135,134
  135 J3(2,1)=M4
      GO TO 136
  134 J3(2,1) = 1
  136 IF(ISIG-2) 137,138,137
  138 J3(2,2)=M4
      GO TO 139
  137 J3(2,2)=2
  139 J3(2,3)=2*IHSH
      IF(IHSH-3) 149,140,140
  140 DO 148 J=4,IHSHP1
      L=J-1
      J2(J,1)=M1+L
      J2(J,3)=M1+J
      J3(L,1)=M2+L
      J3(L,3)=M2+J
  141 IF(IRHO-L) 142,143,142
  143 J2(J,2)=M3
      GO TO 144
  142 J2(J,2)=L
  144 IF(ISIG-L) 145,146,145
  146 J3(L,2)=M4
       GO TO 148
  145 J3(L,2)=L
  148 CONTINUE
  149 M6=IHSHP1
      J3(M6,1)=M3-1
      J3(M6,2)=MN1
      J3(M6,3)=I2HSH
      IF(IHSH-1) 450,451,450
  451 J3(M6,1) = M4
      J3(M6,3) = M3
  450 NJ23S=M6+1
      IL=2
      IF(NOVLPS.NE.0) CALL CNDENS(IL)
      NJ23=NJ23S-1
      DO 150 J=1,NJ23
      DO 151 K=1,3
      J2STO(J,K)=J2(J,K)
      J3STO(J,K)=J3(J,K)
  151 CONTINUE
  150 CONTINUE
*
*   SUM OVER PARENTS OF IRHO AND ISIG
*
      N1=NOSH1(IRHO)
      K1=NTAB1(N1,L1)
      KK1=ITAB(K1)
      N2=NOSH2(ISIG)
      K2=NTAB1(N2,L2)
      KK2=ITAB(K2)
      DO 111 JJ1=1,KK1
*
*   CHECK ON TRIANGULAR CONDITIONS
*
      IN3=2*LRHO
      IJK1=3*(JJ1-1)+JTAB(K1)
      DO 113 KK=2,3
      IN1=NTAB(IJK1+KK)-1
      IN2=J1QN1(IRHO,KK)-1
      IF(IN1.GT.(IN2+IN3)) GO TO 111
      IF(IN1.LT.IABS(IN2-IN3)) GO TO 111
      IN3=1
  113 CONTINUE
      DO 112 JJ2=1,KK2
      IN3=2*LSIG
      IJK2=3*(JJ2-1)+JTAB(K2)
      DO 114 KK=2,3
      IN1=NTAB(IJK2+KK)-1
      IN2=J1QN2(ISIG,KK)-1
      IF(IN1.GT.(IN2+IN3)) GO TO 112
      IF(IN1.LT.IABS(IN2-IN3)) GO TO 112
      IN3=1
  114 CONTINUE
      DO 115 KK=1,3
      JBAR(IRHO,KK)=NTAB(IJK1+KK)
      JPBAR(ISIG,KK)=NTAB(IJK2+KK)
  115 CONTINUE
      IF(IRHO.EQ.MU.OR.IRHO.EQ.NU) GO TO 116
      DO 117 KK=1,3
      IF(JBAR(IRHO,KK).NE.JPBAR(IRHO,KK)) GO TO 112
  117 CONTINUE
  116 IF(ISIG.EQ.MUP.OR.ISIG.EQ.NUP) GO TO 118
      DO 119 KK=1,3
      IF(JBAR(ISIG,KK).NE.JPBAR(ISIG,KK)) GO TO 112
  119 CONTINUE
  118 IF(NOVLPS.EQ.1) THEN
         DO 120 KK=1,3
         IF(JBAR(MU,KK).NE.JPBAR(MUP,KK)) GO TO 112
  120    CONTINUE
      ENDIF
  154 K=2
      M7=M3-IHSH
      M9=M7+1
      M11=M3-1
      M12=IHSH-1
      RECUPS=ONE
      VSHEL2=ZERO
      RECUPT=ZERO
      IF(NOVLPS.EQ.2) RECUPT=ONE
*
*     FIRST FRACTIONAL PARENTAGE COEFFICIENT
*
      IF(NBUG6.NE.0) WRITE(IWRITE,221)
      LIJ=LRHO
      COEFP=ONE
      IF(LIJ) 171,173,171
  171 N=NOSH1(IRHO)
      IV1=J1QN1(IRHO,1)
      IL1=(J1QN1(IRHO,2)-1)/2
      IS1=J1QN1(IRHO,3)
      IV2=JBAR(IRHO,1)
      IL2=(JBAR(IRHO,2)-1 )/2
      IS2=JBAR(IRHO,3)
      CALL CFP(LIJ,N,IV1,IL1,IS1,IV2,IL2,IS2,COEFP)
      RECUPS=RECUPS*COEFP
      IF(NOVLPS.EQ.2) RECUPT=RECUPT*COEFP
      IF(NBUG6.NE.0) WRITE(IWRITE,207) RECUPS,RECUPT
*
*     SECOND FRACTIONAL PARENTAGE COEFFICIENT
*
  173 LIJ=LSIG
      COEFP=ONE
      IF(LIJ) 176,176,174
  174 N=NOSH2(ISIG)
      IV1=J1QN2(ISIG,1)
      IL1=(J1QN2(ISIG,2)-1)/2
      IS1=J1QN2(ISIG,3)
      IV2=JPBAR(ISIG,1)
      IL2=(JPBAR(ISIG,2)-1)/2
      IS2=JPBAR(ISIG,3)
      CALL CFP(LIJ,N,IV1,IL1,IS1,IV2,IL2,IS2,COEFP)
  176 RECUPS=RECUPS*COEFP
      IF(NOVLPS.EQ.2) RECUPT=RECUPT*COEFP
      IF(NBUG6.NE.0) WRITE(IWRITE,207) RECUPS,RECUPT
      IF(DABS(RECUPS).LT.1.D-14.AND.DABS(RECUPT).LT.1.D-14) GO TO 112
*
*     SET UP THE J1 ARRAY FOR THE ANGULAR AND SPIN RECOUPLING
*     COEFFICIENTS
*
  155 IF(K-3) 156,157,157
  156 J1(M5)=2*LRHO+1
      J1(M10)=2*LSIG+1
      J1(MN1)=2*KA+1
      IF(ISPIN.EQ.1) J1(MN1)=1
      J1(M3)=J1QN1(IRHO,K)
      J1(M4)=J1QN2(ISIG,K)
      GO TO 158
  157 J1(M5)=2
      J1(M10)=2
      J1(MN1)=KB+KB+1
      IF(ISPIN.EQ.0) J1(MN1)=1
      J1(M3)=J1QN1(IRHO,K)
      J1(M4)=J1QN2(ISIG,K)
  158 DO 161 J=1,IHSH
      J1(J)=JBAR(J,K)
      IF(J.EQ.MUP.OR.J.EQ.NUP) J1(J)=JPBAR(J,K)
  161 CONTINUE
      IF(IHSH.EQ.1) GO TO 197
      DO 162 J=M6,M7
      J1(J)=J1QN1(J,K)
  162 CONTINUE
      DO 163 J=M9,M11
      JM12=J-M12
      J1(J)=J1QN2(JM12,K)
  163 CONTINUE
*
*     PRINT OUT THE J1,J2 AND J3 ARRAYS
*
  197 IF(NBUG6.EQ.0) GO TO 304
   17 WRITE(IWRITE,209)
      WRITE(IWRITE,210) (J1(J),J=1,NJ1S)
      WRITE(IWRITE,211)
      DO 166 I=1,NJ23
      WRITE(IWRITE,212) (J2(I,J),J=1,3),(J3(I,J),J=1,3)
  166 CONTINUE
      IF(K.EQ.2) WRITE(IWRITE,213)
  304 CONTINUE
      IF(DABS(RECUPS).LT.1.D-14) GO TO 1164
*
*     EVALUATE ORBITAL AND SPIN RECOUPLING COEFFICIENTS
*
  164 IF(K.EQ.3.AND.NBUG6.NE.0) WRITE(IWRITE,214)
      IF(NOVLPS.LT.2.OR.MU.EQ.NU.OR.MUP.EQ.NUP) GO TO 77
      IF(IMATCH(1).EQ.0) GO TO 79
      IF(LMU.NE.LMUP.OR.LNU.NE.LNUP) GO TO 79
      DO 937 KK=1,3
      IF(JBAR(MU,KK).NE.JPBAR(MUP,KK)) GO TO 79
      IF(JBAR(NU,KK).NE.JPBAR(NUP,KK)) GO TO 79
  937 CONTINUE
   77 DO 78 I = 1,NJ1S
        FREE(I) = .FALSE.
   78 CONTINUE
      CALL NJGRAF(RECUP,FAIL)
      GO TO 81
   79 RECUP=ZERO
   81 RECUPS=RECUPS*RECUP
      IF(NBUG6.NE.0) WRITE(IWRITE,207) RECUPS,RECUP
 1164 IF(NOVLPS.LT.2) GO TO 170
      IF(IMATCH(2).EQ.0) GO TO 84
      IF(LMU.NE.LNUP.OR.LNU.NE.LMUP) GO TO 84
      IF(MU.EQ.NU.OR.MUP.EQ.NUP) GO TO 84
      DO 938 KK=1,3
      IF(JBAR(MU,KK).NE.JPBAR(NUP,KK)) GO TO 84
      IF(JBAR(NU,KK).NE.JPBAR(MUP,KK)) GO TO 84
  938 CONTINUE
   68 DO 82 J=1,NJ23
      DO 83 KK=1,3
      J2(J,KK)=J2STO(J,KK)
      J3(J,KK)=J3STO(J,KK)
   83 CONTINUE
   82 CONTINUE
      JSTO=J3(IROWMU,ICOLMU)
      J3(IROWMU,ICOLMU)=J3(IROWNU,ICOLNU)
      J3(IROWNU,ICOLNU)=JSTO
*
*     EVALUATE ORBITAL AND SPIN RECOUPLING COEFFICIENTS
*
      DO 500 I = 1,NJ1S
         FREE(I) = .FALSE.
  500 CONTINUE
*
      CALL NJGRAF(RECUP,FAIL)
*
   86 RECUPT=RECUPT*RECUP
      IF(NBUG6.NE.0) WRITE(IWRITE,207) RECUPT,RECUP
      GO TO 170
   84 RECUP=ZERO
      GO TO 86
  170 K=K+1
      DO 168 J=1,NJ23
      DO 169 KK=1,3
      J2(J,KK)=J2STO(J,KK)
      J3(J,KK)=J3STO(J,KK)
  169 CONTINUE
  168 CONTINUE
      IF(K.EQ.3) GO TO 155
*
*     PERMUTATION FACTOR
*
  182 VALML=RECUPS
      IF(NOVLPS.EQ.2) VALML2=RECUPT
      RML=RML+VALML
      IF(NOVLPS.EQ.2) RML2=RML2+VALML2
      IF(NBUG6.NE.0) WRITE(IWRITE,2005) VALML,VALML2,RML,RML2
 2005 FORMAT(30X,4F10.6)
  112 CONTINUE
  111 CONTINUE
      GO TO 184
  183 RML=ZERO
      IF(NOVLPS.EQ.2) RML2=ZERO
  184 SQRN=DSQRT(DFLOAT(NOSH1(IRHO)*NOSH2(ISIG)))
      FACTR=FACTR*SQRN*DFLOAT(MINUS)*DSQRT(AJF)
      VSHELL=RML*FACTR
      IF(NOVLPS.EQ.2) VSHEL2=RML2*FACTR*(-ONE)**NKAP

      CALL PRNTML(IRHO,ISIG,VSHELL,VSHEL2,KA,KB,ISPIN)

  189 IF(NOVSTO.EQ.0) GO TO 289
      MU=MUSTO
      NU=NUSTO
      MUP=MUPSTO
      NUP=NUPSTO
      NOVLPS=NOVSTO
      LMU=LMUSTO
      LNU=LNUSTO
      LMUP=LMUPST
      LNUP=LNUPST
  289 RML=ZERO
      RML2=ZERO
      IF(ISUM.EQ.0) GO TO 190
      IF(IBACK3.EQ.1) GO TO 509
      IRHO=IRHO+1
      ISIG=IRHO
      IF(NBUG6.NE.0) WRITE(IWRITE,1216) IRHO,ISIG
 1216 FORMAT(15H AT 289, IRHO =,I3,7H ISIG =,I3)
      IF(IRHO.LE.IHSH) GO TO 1509
      RETURN
*
*     NO SUM OVER SHELLS IBACK3=1 IMPLIES ISIG=NUP OR IRHO=NU
*     TO BE CONSIDERED ALSO
*
  190 IF(IBACK3.EQ.1) GO TO 506
      RETURN
  300 IF(NBUG6.NE.0) WRITE(IWRITE,302)
      RETURN
  301 IF(NBUG6.NE.0) WRITE(IWRITE,303)
      RETURN
      END
*
*     ------------------------------------------------------------------
*        C N D E N S
*     ------------------------------------------------------------------
*
      SUBROUTINE CNDENS(I1L)
*
      PARAMETER (NWD=30,NCD=100,NORD=NWD*(NWD-1)/2)
      PARAMETER(KFL1=60,KFL2=12)
      LOGICAL FREE
      COMMON/COUPLE/NJ1S,NJ23S,J1(KFL1),J2(KFL2,3),J3(KFL2,3),FREE(KFL1)
      COMMON/OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
     : ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
     :     IORTH(NORD)
*
* --- Links the co-ordinates of electrons in non-orthogonal subshells 
*     on the two sides of the martix element, and condenses the angular 
*     momentum coupling schemes accordingly.
*
      NJ23ST=NJ23S
      NJ23=NJ23S-1
      CALL PIKUP1(J2,NJ23,MUP,NUP,IROWMU,IROWNU,ICOLMU,ICOLNU)
      IF(IROWMU.NE.IROWNU) GO TO 1
      IROWM1=IROWMU+1
      J2(IROWM1+1,1)=J2(IROWM1,2)
      CALL SCRAP(J2,NJ23,IROWMU)
      CALL SCRAP(J2,NJ23,IROWMU)
      GO TO 3
    1 IMU=J2(IROWMU,3)
      CALL PIKUP2(J2,NJ23,IROWMU,IMU,IROW,ICOL)
      IF(IROWMU.EQ.NJ23) GO TO 2
      J2(IROW,ICOL)=J2(IROWMU,3-ICOLMU)
    2 IF(NOVLPS.EQ.1.OR.MUP.EQ.NUP) GO TO 4
      IMU=J2(IROWNU,3)
      CALL PIKUP2(J2,NJ23,IROWNU,IMU,IROW,ICOL)
      IF(IROWNU.EQ.NJ23) GO TO 4
      J2(IROW,ICOL)=J2(IROWNU,3-ICOLNU)
    4 CALL SCRAP(J2,NJ23,IROWMU)
      IF(NOVLPS.LT.2.OR.MUP.EQ.NUP) GO TO 3
      IROWNU=IROWNU-1
      CALL SCRAP(J2,NJ23,IROWNU)
    3 NJ23=NJ23ST-1
      CALL PIKUP1(J3,NJ23,MU,NU,IROWMU,IROWNU,ICOLMU,ICOLNU)
      IF(IROWMU.NE.IROWNU) GO TO 11
      IROWM1=IROWMU+1
      J3(IROWM1+1,1)=J3(IROWM1,2)
      CALL SCRAP(J3,NJ23,IROWMU)
      CALL SCRAP(J3,NJ23,IROWMU)
      GO TO 13
   11 IMU=J3(IROWMU,3)
      CALL PIKUP2(J3,NJ23,IROWMU,IMU,IROW,ICOL)
      IF(IROWMU.EQ.NJ23) GO TO 12
      J3(IROW,ICOL)=J3(IROWMU,3-ICOLMU)
   12 IF(NOVLPS.EQ.1.OR.MU.EQ.NU) GO TO 14
      IMU=J3(IROWNU,3)
      CALL PIKUP2(J3,NJ23,IROWNU,IMU,IROW,ICOL)
      IF(IROWNU.EQ.NJ23) GO TO 14
      J3(IROW,ICOL)=J3(IROWNU,3-ICOLNU)
   14 CALL SCRAP(J3,NJ23,IROWMU)
      IF(NOVLPS.LT.2.OR.MU.EQ.NU) GO TO 13
      IROWNU=IROWNU-1
      CALL SCRAP(J3,NJ23,IROWNU)
   13 IF(MU.EQ.NU.AND.MUP.NE.NUP) GO TO 21
      IF(MUP.EQ.NUP.AND.MU.NE.NU) GO TO 22
      CALL PIKUP1(J3,NJ23,MUP,NUP,IROWMU,IROWNU,ICOLMU,ICOLNU)
      J3(IROWMU,ICOLMU)=MU
      IF(IROWNU.EQ.0) GO TO 23
      J3(IROWNU,ICOLNU)=NU
      GO TO 23
   21 NJ23=NJ23ST-3
      IROW=3-I1L
      CALL NSCRAP(J2,NJ23,IROW)
      J2(IROW,1)=MUP
      J2(IROW,2)=NUP
      J2(IROW,3)=MU
      GO TO 23
   22 NJ23=NJ23ST-3
      IROW=3-I1L
      CALL NSCRAP(J3,NJ23,IROW)
      J3(IROW,1)=MU
      J3(IROW,2)=NU
      J3(IROW,3)=MUP
   23 NJ23S=NJ23+1
      RETURN
      END
*
*     ------------------------------------------------------------------
*        M A T C H
*     ------------------------------------------------------------------
*
      INTEGER FUNCTION MATCH(IA, IB, IC, ID)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NCD=100,NORD=NWD*(NWD-1)/2)
*
      COMMON/OVRINT/IOVEL(2),IOVER(2),IOVEP(2)
      COMMON/OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
     : ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
     :     IORTH(NORD)
      COMMON/MEDEFN/IHSH,NJ(10),LJ(10),NOSH1(10),NOSH2(10),J1QN1(19,3),
     :     J1QN2(19,3),IJFUL(10)
      COMMON/INFORM/IREADC,IWRITE,IOUT,IREADJ,IREADW,ISCW,IREAD,ISC(4)
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
*     COMMON/DIAGNL/IDIAG,JA,JB
      COMMON/EMS/IEM(4),JA,JB,IFL

      COMMON/MACHOR/IMATCH(2)
*
* --- Tries to match the angular momenta of the spectator
*     non-orthogonal orbitals.
*
      IMATCH(1)=0
      IMATCH(2)=0
      IM = 0
      NA = 0
      NB = 0
      NC = 0
      ND = 0
  200 IF (MU .NE. 0) NA = NOSH1(MU)
      IF (NU .NE. 0 .AND. MU .NE. NU) NB = NOSH1(NU)
      IF (MUP .NE. 0) NC = NOSH2(MUP)
      IF (NUP .NE. 0 .AND. MUP .NE. NUP) ND = NOSH2(NUP)
      NA = NA - IA
      NB = NB - IB
      NC = NC - IC
      ND = ND - ID
      IF (NA.LT.0 .OR. NB.LT.0 .OR. NC.LT.0 .OR. ND.LT.0) GO TO 20
      IF ((NA+NB) .EQ. (NC+ND)) GO TO 14
   20 MATCH = 1
      RETURN
   14 IF (NA.EQ.0 .AND. NB.EQ.0 .AND. NC.EQ.0 .AND. ND.EQ.0) GO TO 17
      IF ((NA.EQ.0 .OR. NB.EQ.0).AND.(NC.EQ.0 .OR. ND.EQ.0)) GO TO 15
      NOVLPS=2
      GO TO 16
   15 IF (NA .EQ. 0) MU = NU
      IF (NC .EQ. 0) MUP = NUP
      NU = 0
      NUP = 0
      NOVLPS = 1
      LMU = LJ(MU)
      LMUP = LJ(MUP)
      IF (LMU .NE. LMUP) GO TO 20
      IOVEL(1) = MU
      IOVER(1) = MUP
      IOVEP(1) = NA + NB
      MATCH = 2
      GO TO 30
   17 MU = 0
      NU = 0
      MUP = 0
      NUP = 0
      NOVLPS = 0
      MATCH = 2
      RETURN
   16 IF(NA.EQ.0) MU=NU
      IF(NB.EQ.0) NU=MU
      IF(NC.EQ.0) MUP=NUP
      IF(ND.EQ.0) NUP=MUP
      LA = LJ(MU)
      LB = LJ(NU)
      LC = LJ(MUP)
      LD = LJ(NUP)
      IF (LA .EQ. LB) GO TO 21
      IF (LA .EQ. LC) GO TO 22
      IF (LA .NE. LD .OR. NA .NE. ND) GO TO 20
      IF (LB .NE. LC) GO TO 20
      IOVEL(1) = MU
      IOVER(1) = NUP
      IOVEP(1) = NA
      IOVEL(2) = NU
      IOVER(2) = MUP
      IOVEP(2) = NB
      MATCH = 2
      IM = 1
      GO TO 30
   22 IF (LB.NE.LD .OR. NB.NE.ND) GO TO 20
      IOVEL(1) = MU
      IOVER(1) = MUP
      IOVEP(1) = NA
      IOVEL(2) = NU
      IOVER(2) = NUP
      IOVEP(2) = NB
      MATCH = 2
      GO TO 30
   21 IF (LA.NE.LC .OR. LA.NE.LD) GO TO 20
      IF(MU .EQ. NU) GO TO 23
      IF (MUP .EQ. NUP) GO TO 24
      IF (NA.GT.1 .OR. NB.GT.1 .OR. NC.GT.1 .OR. ND.GT.1) GO TO 25
      IOVEL(1) = MU
      IOVER(1) = MUP
      IOVEP(1) = 1
      IOVEL(2) = NU
      IOVER(2) = NUP
      IOVEP(2) = 1
      MATCH = 2
      GO TO 30
   25 WRITE(IWRITE,300) MU,NU,MUP,NUP,JA,JB
  300 FORMAT('THE FOLLOWING SUBSHELLS HAVE A COMMON L-VALUE BUT',
     :  ' CONTAIN TOO MANY ELECTRONS FOR THIS CODE',5X,4I3/
     :    5X,'THE MATRIX ELEMENT IS (',I2,1H,,I2,1H))
      MATCH = 0
      RETURN
   23 IF ((NC+ND) .GT. 2) GO TO 25
      IOVEL(1) = MU
      IOVER(1) = MUP
      IOVEP(1) = NC
      IOVEL(2) = MU
      IOVER(2) = NUP
      IOVEP(2) = ND
      MATCH = 2
      GO TO 30
   24 IF ((NA+NB) .GT. 2) GO TO 25
      IOVEL(1) = MU
      IOVER(1) = MUP
      IOVEP(1) = NA
      IOVEL(2) = NU
      IOVER(2) = MUP
      IOVEP(2) = NB
      MATCH = 2
*
*      FINAL CHECK ON ALLOWED NON-ORTHOGONALITY
*
   30 K1=1
      K2=1
      IGO=1
   31 I1=IOVEL(K1)
      I2=IOVER(K2)
      I5=IJFUL(I1)
      I6=IJFUL(I2)
      I3=MIN0(I5,I6)
      I4=MAX0(I5,I6)
      I1=(I4-2)*(I4-1)/2+I3
      I2=IORTH(I1)
      GO TO (32,33,35,36),IGO
   32 IF(I2.NE.1) GO TO 34
      IF(NOVLPS.EQ.1) GO TO 33
      K1=2
      K2=2
      IGO=2
      GO TO 31
   33 IF(I2.EQ.1) IMATCH(IM+1)=1
   34 IF(NOVLPS.EQ.1) THEN
         IF(IMATCH(1).EQ.0) GO TO 20
         RETURN
      ENDIF
      K1=1
      K2=2
      IGO=3
      GO TO 31
   35 IF(I2.NE.1) GO TO 37
      K1=2
      K2=1
      IGO=4
      GO TO 31
   36 IF(I2.EQ.1) IMATCH(2-IM)=1
   37 IF(IMATCH(1).EQ.0.AND.IMATCH(2).EQ.0) GO TO 20
      RETURN
      END
*
*     ------------------------------------------------------------------
*        N O R T B P
*     ------------------------------------------------------------------
*
      SUBROUTINE NORTBP(JA,JB)
      PARAMETER (NWD=30,NCD=100,NORD=NWD*(NWD-1)/2)
*
      DIMENSION ILNO(2),IRNO(2)
      COMMON/INFORM/IREADC,IWRITE,IOUT,IREADJ,IREADW,ISCW,IREAD,ISC(4)
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
     : ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
     :     IORTH(NORD)
      COMMON/STATES/NCFG,MAXORB,IAJCMP(NWD),LJCOMP(NWD),
     :NJCOMP(NWD),NOCCSH(NCD),NELCSH(5,NCD),NOCORB(5,NCD),J1QNRD(9,NCD)
*
* === Determines, in each of the two configurations JA and JB in a 
*     matrix element, which subshells are non-orthogonal.
*
  101 FORMAT(/63H INCORRECT NON-ORTHOGONALITY SET UP IN THE MATRIX ELEME
     :NT  -  (,I3,3H/V/,I3,1H))
*
      N1=NOCCSH(JA)
      N2=NOCCSH(JB)
      JMU=0
      JNU=0
      JMUP=0
      JNUP=0
      IL=0
      IR=0
*
* --- BEGIN SEARCH FOR NON-ORTHOGONAL SUBSHELLS IN THIS MATRIX ELEMENT
*
      DO 1 I=1,N1
      NI=NOCORB(I,JA)
      DO 2 J=1,N2
      NJ=NOCORB(J,JB)
      IF(NI.EQ.NJ) GO TO 2
      NA=MIN0(NI,NJ)
      NB=MAX0(NI,NJ)
      NC=(NB-2)*(NB-1)/2+NA
      IF(IORTH(NC).NE.1) GO TO 2
      IF (IL .EQ. 0) GO TO 4
      IF (ILNO(IL) .EQ. I) GO TO 14
      IF (IL .EQ. 2) GO TO 100
    4 IL = IL+1
      ILNO(IL) = I
   14 IF ( IR .EQ. 0) GO TO 7
      DO 15 K = 1,IR
      IF (IRNO(K) .EQ. J) GO TO 2
   15 CONTINUE
      IF (IR .EQ. 2) GO TO 100
    7 IR = IR+1
      IRNO(IR) = J
    2 CONTINUE
    1 CONTINUE
      IF(IL.EQ.0) GO TO 8
      IF(IR.EQ.1) GO TO 11
      IF(IRNO(1).LE.IRNO(2)) GO TO 11
      ISTO=IRNO(1)
      IRNO(1)=IRNO(2)
      IRNO(2)=ISTO
   11 JMU=ILNO(1)
    3 IF(IL.EQ.1) GO TO 5
      JNU=ILNO(2)
    5 JMUP=IRNO(1)
    6 IF (IR .EQ. 1) GO TO 10
      JNUP=IRNO(2)
      GO TO 10
    8 NOVLPS=0
      RETURN
  100 WRITE(IWRITE,101) JA,JB
      STOP
   10 NMU=NOCORB(JMU,JA)
      NMUP=NOCORB(JMUP,JB)
      LMU=LJCOMP(NMU)
      LMUP=LJCOMP(NMUP)
      IF (JNU .EQ. 0 .AND. JNUP .EQ. 0) GO TO 9
      IF (JNU .EQ. 0) JNU = JMU
      IF (JNUP .EQ. 0) JNUP = JMUP
      NNU=NOCORB(JNU,JA)
      NNUP=NOCORB(JNUP,JB)
      LNU=LJCOMP(NNU)
      LNUP=LJCOMP(NNUP)
      NOVLPS=2
      RETURN
    9 NOVLPS=1
      RETURN
      END
*
*     ------------------------------------------------------------------
*        N S C R A P
*     ------------------------------------------------------------------
*
      SUBROUTINE NSCRAP(IX,IRS,IR1)
      PARAMETER(KFL2=12)
      DIMENSION IX(KFL2,3)
*
* === Create a blank row in array IX after row (IR1 -1)
*
      IR2=IRS+1-IR1
      DO 1 I=1,IR2
      II=IRS+1-I
      I1=II+1
      DO 2 K=1,3
      IX(I1,K)=IX(II,K)
    2 CONTINUE
    1 CONTINUE
      IRS=IRS+1
      RETURN
      END
*
*     ------------------------------------------------------------------
*        P I K U P 1
*     ------------------------------------------------------------------
*
      SUBROUTINE PIKUP1(IX,IRS,MU,NU,IR1,IR2,IC1,IC2)
      PARAMETER(KFL2=12)
      DIMENSION IX(KFL2,3)
*
* === Locates the position in array IX of MU and NU
*
      IR1=0
      IR2=0
      DO 1 I=1,IRS
      DO 2 K=1,2
      IA=IX(I,K)
      IF(IA.EQ.MU) GO TO 3
      IF(IA.EQ.NU) GO TO 4
      GO TO 2
    3 IR1=I
      IC1=K
      IF(MU.EQ.NU) RETURN
      GO TO 2
    4 IR2=I
      IC2=K
    2 CONTINUE
    1 CONTINUE
      RETURN
      END
*
*     ------------------------------------------------------------------
*        P I K U P 2
*     ------------------------------------------------------------------
*
      SUBROUTINE PIKUP2(IX,IRS,IR1,IMU,IR2,IC2)
      PARAMETER(KFL2=12)
      DIMENSION IX(KFL2,3)
*
* === Locates IMU in IX
*
      IR2=0
      IR3=IR1+1
      DO 1 I=IR3,IRS
      DO 2 K=1,2
      IA=IX(I,K)
      IF(IA.NE.IMU) GO TO 2
      IR2=I
      IC2=K
      RETURN
    2 CONTINUE
    1 CONTINUE
      RETURN
      END
*
*
*     ------------------------------------------------------------------
*        S C R A P
*     ------------------------------------------------------------------
*
      SUBROUTINE SCRAP(IX,IRS,IR1)
      PARAMETER(KFL2=12)
      DIMENSION IX(KFL2,3)
*
* === Deletes row IR1 from the coupling array I1
*
      IF(IR1.EQ.IRS) GO TO 3
      IR2=IR1+1
      DO 1 I=IR2,IRS
      I1=I-1
      DO 2 K=1,3
      IX(I1,K)=IX(I,K)
    2 CONTINUE
    1 CONTINUE
      GO TO 4
    3 IX(IRS-1,3)=IX(IRS,3)
    4 IRS=IRS-1
      RETURN
      END
*
*
*     ------------------------------------------------------------------
*        S E T U P
*     ------------------------------------------------------------------
*
      SUBROUTINE SETUP(JA,JB)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NCD=100,NORD=NWD*(NWD-1)/2)
*
      COMMON/INFORM/IREADC,IWRITE,IOUT,IREADJ,IREADW,ISCW,IREAD,ISC(4)
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/STATES/NCFG,MAXORB,IAJCMP(NWD),LJCOMP(NWD),
     :NJCOMP(NWD),NOCCSH(NCD),NELCSH(5,NCD),NOCORB(5,NCD),J1QNRD(9,NCD)
      COMMON/MEDEFN/IHSH,NJ(10),LJ(10),NOSH(10,2),J1QN(19,3,2),IJFUL(10)
      COMMON/OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
     : ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
     :     IORTH(NORD)
*
*     NOTICE THE DIFFERENT NAMES IN THE COMMON BLOCK MEDEFN  -  WE
*      STORE NOSH1(I=1,10) AS NOSH((I=1,10),1) AND NOSH2(I=1,10) AS
*     NOSH((I=1,10),2)   AND USE THE FACT THAT NOSH1 AND NOSH2 WILL THEN
*     BE EQUIVALENT TO THE SINGLE 2-DIMENSIONAL ARRAY NOSH.  SIMILARLY
*     FOR J1QN
*
* === GENERATES THE ARRAYS  NJ,LJ - DEFINING THE QUANTUM NUMBERS OF THE
*     SHELLS,   NOSH - DEFINING THE OCCUPATION OF THE SHELLS,  J1QN -
*     DEFINING THE COUPLING OF THE SHELLS,   FOR EACH OF THE TWO
*     CONFIGURATIONS CONSIDERED.    ONLY THOSE SHELLS OCCURRING IN AT
*     LEAST ONE CONFIGURATION ARE INCLUDED.
*                   AT LEAST TWO SHELLS MUST BE CONSIDERED OCCUPIED.
*     THUS (1S)**2    HELIUM  MUST BE TREATED AS ,E.G., (1S)**2(2S)**0
*     THE SIZE OF THE ARRAYS HERE CALCULATED IS ARRANGED TO BE NO
*     GREATER THAN IS NECESSARY TO INCLUDE ALL ORBITALS WHICH ARE
*     DEEMED TO BE OCCUPIED IN EITHER OR BOTH OF THE CONFIGURATIONS
*     JA,JB
*
* --- INITIALIZE BASIC QUANTITIES - (I1+1) RUNS OVER 1,MAXORB,  IHSH IS
*     THE CURRENT VALUE OF THE HIGHEST OCCUPIED SHELL YET CONSIDERED,
*     WHILE I2HSH=2*IHSH-1
*
      MU=0
      NU=0
      MUP=0
      NUP=0
      I1=0
      IHSH=0
      I2HSH=-1
      IA=NOCCSH(JA)
      IB=NOCCSH(JB)
*
* --- TEST ON WHETHER LIMIT OF I1 HAS BEEN REACHED
*
    1 IF(I1-MAXORB) 101,100,100
*
* --- INCREASE BASIC QUANTITIES
*
  101 I1=I1+1
      I3=IHSH+1
      I5=I2HSH+I3
*
* --- IS THE SHELL I1 OCCUPIED IN JA
*
      DO 2 J=1,IA
      IF(I1-NOCORB(J,JA)) 2,3,2
    2 CONTINUE
      NA=1
      GO TO 4
    3 NA=2
      J1=J
*
* --- IS THE SHELL I1 OCCUPIED IN JB
*
    4 DO 5 J=1,IB
      IF(I1-NOCORB(J,JB)) 5,6,5
    5 CONTINUE
      NB=1
      GO TO 7
    6 NB=2
      J2=J
*
*     IF THE SHELL I1 IS NOT OCCUPIED IN EITHER JA OR JB, IGNORE THE
*     SHELL, DO NOT INCREASE IHSH, AND CONSIDER NEXT SHELL BY INCREASING
*     I1
*
    7 IF(NA-1) 8,8,9
    8 IF(NB-1) 1,1,9
*
* --- IF THE SHELL I1 IS OCCUPIED IN EITHER JA OR JB -
*     (1)   IF IHSH.GT.1, THEN ALREADY AT LEAST TWO SHELLS AND THE
*     RESULTING COUPLINGS HAVE BEEN STORED. WE MUST THUS MAKE ROOM FOR
*     THE QUANTUM NUMBERS OF THIS NEW SHELL BETWEEN THE QUANTUM NUMBERS
*     OF THE PREVIOUS SHELLS AND THE QUANTUM NUMBERS OF THE INTERMEDIATE
*     COUPLINGS OF THE CONFIGURATIONS.  THUS THE LATTER SET ARE =MOVED
*     ALONG= TO MAKE ROOM FOR THE NEW SHELL
*     (2)   IF IHSH.LE.1, THERE ARE NO INTERMEDIATE COUPLING QUANTUM
*     NUMBERS, AND SO THERE IS NOTHING TO MOVE
*
    9 IF(IHSH-1) 11,11,10
   10 DO 12 I=1,2
      DO 13 J=I3,I2HSH
      I4=I5-J
      DO 14 K=1,3
      J1QN(I4+1,K,I)=J1QN(I4,K,I)
   14 CONTINUE
   13 CONTINUE
   12 CONTINUE
   11 IHSH=I3
      I2HSH=I2HSH+2
      NC=NA
      I=1
      IC=J1
      JC=JA
*
* --- FIRST CONSIDER THE L.H.S. (I=1) OF THE MATRIX ELEMENT.  NC=1 MEANS
*     UNOCCUPIED, REPRESENTED BY A DUMMY SINGLET S SHELL, AND THE
*     ADDITIONAL SET OF COUPLING QUANTUM NUMBERS WILL BE THE SAME AS THE
*     LAST SET OF COUPLING QUANTUM NUMBERS ALREADY OBTAINED.
*     NC=2 MEANS OCCUPIED.  THEN ALL THE NEW QUANTUM NUMBERS (BOTH FOR
*     THE SHELL AND FOR THE COUPLING OF THIS SHELL TO THE RESULTANT OF
*     THE PREVIOUS ONES) ARE DEFINED IN THE CORRESPONDING J1QNRD ARRAY.
*     NOSH - THE NUMBER OF ELECTRONS IN THIS SHELL, IS DEFINED BY THE
*     APPROPRIATE ENTRY IN NELCSH .  THE R.H.S. IS THEN CONSIDERED
*     SIMILARLY (I=2)
*
   25 GO TO (15,16),NC
   15 NOSH(IHSH,I)=0
      J1QN(IHSH,1,I)=0
      J1QN(IHSH,2,I)=1
      J1QN(IHSH,3,I)=1
      IF(IHSH-2) 22,18,19
   18 J1QN(3,1,I)=0
      J1QN(3,2,I)=J1QN(1,2,I)
      J1QN(3,3,I)=J1QN(1,3,I)
      GO TO 22
   19 DO 27 K=1,3
      J1QN(I2HSH,K,I)=J1QN(I2HSH-1,K,I)
   27 CONTINUE
      GO TO 22
   16 NOSH(IHSH,I)=NELCSH(IC,JC)
      IF(NOVLPS.EQ.0) GO TO 33
      GO TO (31,32),I
   31 IF(IC.EQ.JMU) MU=IHSH
      IF(IC.EQ.JNU) NU=IHSH
      GO TO 33
   32 IF(IC.EQ.JMUP) MUP=IHSH
      IF(IC.EQ.JNUP) NUP=IHSH
   33 JD = J1QNRD(IC,JC)
      J1QN(IHSH,1,I)=MOD(JD,64)
      JD = JD/64
      J1QN(IHSH,2,I) = MOD(JD,64)
      J1QN(IHSH,3,I) = JD/64
*
*     IS THIS THE FIRST OCCUPIED SHELL OF EITHER CONFIGURATION. IF SO,
*     THEN THERE ARE NO INTERMEDIATE COUPLINGS TO CONSIDER AT THIS STAGE
*
      IF(IHSH .GT. 1) THEN
*
*     IS THIS THE FIRST OCCUPIED SHELL OF THIS CONFIGURATION, THOUGH NOT
*     THE FIRST OF THE OTHER CONFIGURATION.  IF SO, THE INTERMEDIATE
*     COUPLING FORMED HAS THE SAME  L,S  VALUES AS THIS OCCUPIED SHELL,
*     SINCE WE COUPLE THE SHELL TO A DUMMY SINGLET S.
*
         IF(IC .LE.1) THEN 
            I2 = 1
         ELSE
            I2 = NOCCSH(JC)+IC-1
         END IF
         JD = J1QNRD(I2,JC)
         IF (IC .LE. 1) THEN
            J1QN(I2HSH,1,I) = 0
         ELSE
            J1QN(I2HSH,1,I) = MOD(JD,64)
         END IF
         JD = JD/64
         J1QN(I2HSH,2,I) = MOD(JD, 64)
         J1QN(I2HSH,3,I) = JD/64
      END IF
*
*     SENIORITY SET (ARBITRARILY) ZERO FOR INTERMEDIATE COUPLING
*
   22 IF(I-2) 23,24,24
   23 NC=NB
      I=2
      IC=J2
      JC=JB
      GO TO 25
*
* --- SET THE NJ AND LJ VALUES OF THE OCCUPIED SHELLS
*
   24 NJ(IHSH)=NJCOMP(I1)
      IJFUL(IHSH)=I1
      LJ(IHSH)=LJCOMP(I1)
*
* --- RETURN TO 1  TO SEE IF MAXORB HAS BEEN REACHED
*
      GO TO 1
  100 CONTINUE
      END
     

*--------------------------------------------------------------------
*     R A D I A L 1
*--------------------------------------------------------------------
*
*     This subroutine calculates the radial matrix elements
*     for the orbital, spin dipolar and quadrupole terms.
*

      SUBROUTINE RADIAL1(N,RADINT,OVLINT1,OVLINT2)
  
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON/HYPER/NHY,MHY,IRHY(20),ISHY(20),JRHY(20),JSHY(20),VHY(20),
     1 IMUHY(20),IMUPHY(20),JMUHY(20),JMUPHY(20),
     2 INUHY(20),INUPHY(20),JNUHY(20),JNUPHY(20),NNNOVLP(20)
      COMMON/OVRINT/IOVEL(2),IOVER(2),IOVEP(2)

      RADINT= QUADR(IRHY(N),ISHY(N),-3)
      IF (NNNOVLP(N).EQ.1) THEN
         OVLINT1=(QUADR(IMUHY(N),IMUPHY(N),0))**(IOVEP(1))
         RADINT=RADINT*OVLINT1
      ELSEIF (NNNOVLP(N).GT.1) THEN
         OVLINT1=(QUADR(IMUHY(N),IMUPHY(N),0))**(IOVEP(1))
         OVLINT2=(QUADR(INUHY(N),INUPHY(N),0))**(IOVEP(1))
         RADINT=RADINT*OVLINT1*OVLINT2
      ENDIF
      RETURN
      END

*--------------------------------------------------------------------
*     R A D I A L 2
*--------------------------------------------------------------------
*
*     This subroutine calculates the radial matrix element
*     for the Fermi contact term.
* 

      SUBROUTINE RADIAL2(N,AZSQR,OVLINT1,OVLINT2)
  
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220,NWD=30,NCD=100)
       
      COMMON/ADATA/MNCFG,MNNN,P(NOD,NWD),R(NOD),RR(NOD),R2(NOD),
     :       ZED,W(NCD),AZ(NWD),L(NWD),MAX(NWD)
      COMMON/HYPER/NHY,MHY,IRHY(20),ISHY(20),JRHY(20),JSHY(20),VHY(20),
     1 IMUHY(20),IMUPHY(20),JMUHY(20),JMUPHY(20),
     2 INUHY(20),INUPHY(20),JNUHY(20),JNUPHY(20),NNNOVLP(20)
      COMMON/OVRINT/IOVEL(2),IOVER(2),IOVEP(2)

      AZSQR= AZ(IRHY(N))*AZ(ISHY(N))
      IF (NNNOVLP(N).EQ.1) THEN
         OVLINT1=(QUADR(IMUHY(N),IMUPHY(N),0))**(IOVEP(1))
         AZSQR=AZSQR*OVLINT1
      ELSEIF (NNNOVLP(N).GT.1) THEN
         OVLINT1=(QUADR(IMUHY(N),IMUPHY(N),0))**(IOVEP(1))
         OVLINT2=(QUADR(INUHY(N),INUPHY(N),0))**(IOVEP(2))
         AZSQR=AZSQR*OVLINT1*OVLINT2
      ENDIF
      RETURN
      END
*     ------------------------------------------------------------------
*        P R N T M L
*     ------------------------------------------------------------------
*
*     This subroutine prints out the output from tensor2 on a
*     file <name>.m and has been added as a debug help.
*

      SUBROUTINE PRNTML(IRHO,ISIG,VSHELL,VSHEL2,KA,KB,ISPIN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NCD=100,NORD=NWD*(NWD-1)/2)
      LOGICAL UPDATE
*
      COMMON/INFORM/IREADC,IWRITE,IOUT,IREADJ,IREADW,ISCW,IREAD,ISC(4)
      COMMON/STATES/NCFG,MAXORB,IAJCMP(NWD),LJCOMP(NWD),NJCOMP(NWD),
     : NOCCSH(NCD),NELCSH(5,NCD),NOCORB(5,NCD),J1QNRD(9,NCD)
      COMMON/TENSOR/VVSHELL(50),IIHSH,NNOVLP(50),IIRHO(50),IISIG(50),
     : MMU(50),MMUP(50),NNU(50),NNUP(50)
      COMMON/EMS/IEM(4),JI,JF,IFL
      COMMON/OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
     : ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
     : IORTH(NORD)
      COMMON/OVRINT/IOVEL(2),IOVER(2),IOVEP(2)
      COMMON/MEDEFN/IHSH,NJ(10),LJ(10),NOSH1(10),NOSH2(10),J1QN1(19,3),
     : J1QN2(19,3),IJFUL(10)
    3 FORMAT(2X,F12.8,I5,I5,5X,1H<,A3,2H||,A2,2H||,A3,2H >)
    4 FORMAT(2X,F12.8,I5,I5,5X,1H<,A3,2H||,A2,2H||,A3,2H >,
     : 3X,1H<,A3,1H|,A3,3H>**,I2)
    5 FORMAT(2X,F12.8,I5,I5,5X,1H<,A3,2H||,A2,2H||,A3,2H >,
     : 3X,1H<,A3,1H|,A3,3H>**,I2,
     : 1X,1H<,A3,1H|,A3,3H>**,I2)
    6 FORMAT(F14.8,A2,1H(,A3,I3,1H,,A3,I3,1H):2(1H(,A3,1H|,A3,1H),I2:
     : ))


      UPDATE = .TRUE.
      IF (NOVLPS.EQ.0) THEN
         IF (DABS(VSHELL).LT.1.D-14) RETURN
         IIHSH=IIHSH+1
         VVSHELL(IIHSH)=VSHELL
         IIRHO(IIHSH)=IRHO
         IISIG(IIHSH)=ISIG
         NNOVLP(IIHSH)=0
      ELSEIF (NOVLPS.EQ.1) THEN
         IF (DABS(VSHELL).LT.1.D-14) RETURN
         IIHSH=IIHSH+1
         VVSHELL(IIHSH)=VSHELL
         IIRHO(IIHSH)=IRHO
         IISIG(IIHSH)=ISIG
         MMU(IIHSH)= MU
         MMUP(IIHSH)=MUP 
         NNOVLP(IIHSH)=1
      ELSE
         IF (DABS(VSHELL).LT.1.D-14) GO TO 15
         UPDATE = .FALSE.
         IIHSH=IIHSH+1
         VVSHELL(IIHSH)=VSHELL
         IIRHO(IIHSH)=IRHO
         IISIG(IIHSH)=ISIG
         MMU(IIHSH)= MU
         MMUP(IIHSH)=MUP
         NNU(IIHSH)=NU
         NNUP(IIHSH)=NUP
         NNOVLP(IIHSH)=2
   15    IF (DABS(VSHEL2).GT.1.D-14) THEN
         IF UPDATE IIHSH=IIHSH+1
         VVSHELL(IIHSH)=VSHEL2
         IIRHO(IIHSH)=IRHO
         IISIG(IIHSH)=ISIG
         MMU(IIHSH)= MU
         MMUP(IIHSH)=NUP
         NNU(IIHSH)=NU
         NNUP(IIHSH)=MUP
         NNOVLP(IIHSH)=2
         ENDIF
      ENDIF
*
*     PRINT OUT THE RESULTS AND OUTPUT THEM ON CHANNEL ISC(1) 
*
      RETURN
      END

