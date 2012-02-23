*     -------------------------------------------------------------
*     ISOTOPE SHIFT
*
*                C O P Y R I G H T -- 1994
*
*     Written by C. Froese Fischer,
*                L. Smentek - Mielczarek
*                Vanderbilt University
*                Nashville, TN 37235 USA
*
*                N. Vaeck
*                Van der Waals-Zeeman Laboratory
*                University of Amsterdam
*                The Netherlands
*
*                G. Miecznik
*                University of Lund,
*                S-22 362 Lund, Sweden
*
*     Computer Physics Communication, Vol. 74, 415--431 (1993).
*
*     ----------------------------------------------------------------
*                M A I N  P R O G R A M
*     -----------------------------------------------------------------
*
*     Open the input files, initialize and control the overall
*     calculations.
*
      PROGRAM ISOTOPE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=100,IDIM=550,NCDIM=3000)
      INTEGER OUT,ERR,PRI
      COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,IUL,IUJ,IU(2)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,
     :ID,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
      INTEGER KVAL,IEL,CPTR,IH,JH,OPTR
      COMMON /STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM),
     :VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
      COMMON /RYDBERG/ RYDBRG,ZMU,ENERGY
      LOGICAL PRINT,CORE,YM,YF
*
      COMMON /COR/COREF,COREG,CORER,CORES,ZM1,ZM2,FZE,FZA,DR2,CORE,YM,YF
*
      CHARACTER YES*1,NAME(5)*24,LABEL*72,HEADER*50
      DATA HEADER/' 2*J-Value     Energy       State'/
*
*     Define unit numbers and open files
*
*     UNIT NUMBERS AND FILE NAMES MAY BE MACHINE
*     DEPENDENT. CHECK THE FOLLOWING SECTION.
*
*     IN - Standard input unit, normally the terminal
*     OUT- Standard output unit, normally the terminal
*     ERR- Prompts and Error messages, always the terminal
*     PRI- Printer output unit or file
*     IUC- <name>.c file
*     IUF- <name>.w file
*     IUL- <name>.l file
*     IUJ- <name>.j file
*     IUD- int.lst file
*
      IN  = 5
      OUT = 6
      ERR = 0
      PRI = 3
      IUC = 7
      IUF = 8
      IUD = 9
      IUL =10
      IUJ =11
*
* ****   WRITE OUT HEADER
*
      WRITE(OUT,1)
1     FORMAT(/20X,'==================='/
     :        20X,'   ISOTOPE SHIFT   '/
     :        20X,'==================='/)

*
* ****   OPEN INPUT FILES AND DETERMINE TOLERACE FOR PRINTOUT
*
3     WRITE(ERR,'(/1X,A)') 'Name of State'
      READ(IN,'(A)') NAME(1)
      WRITE(ERR,'(/A/A/A/A)')
     :' Select type of the input file:',
     :'   1  name.c', '   2  name.l ','   3  name.j'
      READ(IN,*) ICASE
      J=INDEX(NAME(1),' ')
      NAME(1)=NAME(1)(1:J-1)//'.c'
      NAME(3)=NAME(1)(1:J-1)//'.w'
      NAME(4)=NAME(1)(1:J-1)//'.l'
      NAME(5)=NAME(1)(1:J-1)//'.j'
      NAME(2)='int.lst'
      OPEN(UNIT=IUC,FILE=NAME(1),STATUS='OLD',ERR=200)
      OPEN(UNIT=IUD,FILE=NAME(2),STATUS='OLD',ERR=200)
      OPEN(UNIT=IUF,FILE=NAME(3),STATUS='OLD',
     :                     FORM='UNFORMATTED',ERR=200)
      IF (ICASE.EQ.2) THEN
         OPEN(UNIT=IUL,FILE=NAME(4),STATUS='OLD',ERR=200)
         IFILE=IUL
      ELSEIF (ICASE.EQ.3) THEN
         OPEN(UNIT=IUJ,FILE=NAME(5),STATUS='OLD',ERR=200)
         IFILE=IUJ
      ENDIF
      GOTO 10
200   WRITE(ERR,'(A)')'Error in OPEN '
      GOTO 3
10    WRITE(ERR,'(/1X,A)') 'Intermediate printing (y/n)'
      READ(IN,'(A1)') YES
      IF (YES.EQ.'Y'.OR.YES.EQ.'y') THEN
         PRINT=.TRUE.
         WRITE(ERR,'(A)') ' Tolerance for printout'
         READ(IN,'(F6.4)')TOL
         WRITE(OUT,21) TOL
21       FORMAT(/38X,'Tolerance for printing = ',F7.5,/)
      ELSE
         PRINT=.FALSE.
      ENDIF
*
* ****  DETERMINE PROBLEM
*
      WRITE(ERR,'(/A)')
     :' Mass shift (1), field shift (2) or both (3) ?'
      READ(IN,*)IPROB
      IF (ICASE.EQ.1) THEN
         ISTATE=1
         JSTATE=1
      ELSEIF (ICASE.EQ.2) THEN
         JSTATE=1
         READ(IFILE,'(A)') SKIP
      ELSE
         WRITE(ERR,'(/A)')' Maximum and minimum value of 2*J'
         READ(IN,*)MAXJ,MINJ
         JSTATE=((MAXJ-MINJ)/2) + 1
         READ(IFILE,'(A)') SKIP
      ENDIF
*
* ****   INITIALIZE
*
      CALL ISDATA(ICASE)
      CALL INITA
      CALL INITR
      CALL INITM(IPROB)
      CALL INTGRL
*
* ****   CALL SMASS AND/OR FIELDSH SUBROUTNIES
*
      CORE=.TRUE.
      DO 750 JJ=1,JSTATE
         IF (ICASE.NE.1) THEN
            READ(IFILE,'(//8X,I4,10X,I4)',END=751) JV,NJ
            ISTATE=NJ
         ENDIF
      DO 750 II=1,ISTATE
         IF (ICASE.NE.1) THEN
             READ(IFILE,'(/6X,F16.8,A72/(7F10.6))')ENERGY,LABEL,
     :       (WT(NC),NC=1,NCFG)
             WRITE(OUT,101)HEADER,JV,ENERGY,LABEL
         ENDIF
      IF (IPROB.NE.2) CALL SMASS(PRINT)
      IF (IPROB.NE.1) CALL FIELDSH(PRINT)
      CORE=.FALSE.
      WRITE(OUT,102)
*
750   CONTINUE
751   CLOSE(IUC)
      CLOSE(IUF)
      CLOSE(IUD)
      CLOSE(IUL)
      CLOSE(IUJ)
*
101   FORMAT(//A//3X,I3,2X,F16.8,A72,/35('-'))
102   FORMAT(//73('=')//)
      STOP
      END
*
*----------------------------------------------------------------------
*       C O E F
*----------------------------------------------------------------------
*
*       This function returns the coefficient of a given Slater
*  integral in the expansion of the energy
*
        DOUBLE PRECISION FUNCTION COEF(INT)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
        INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
        PARAMETER (NCD=100,IDIM=550,NCDIM=3000)
        COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :         ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
        COEF = 0.D0
        IBEGIN = 1
        IF (INT .GT. 1) IBEGIN = CPTR(INT-1)+1
        IEND = CPTR(INT)
        DO 1 II = IBEGIN,IEND
           T = WT(IH(II))*WT(JH(II))*COEFF(II)
           IF (OPTR(II).NE.0) T = T*VALUE(OPTR(II))
           IF (IH(II) .NE. JH(II)) T = T+T
           COEF = COEF+T
  1     CONTINUE
        END
*-----------------------------------------------------------------------
*       D O P R I N T
*-----------------------------------------------------------------------
*
*     Controls the printing
*
      LOGICAL FUNCTION DOPRINT(A,B,C,PRINT,TOL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL PRINT
*
      DOPRINT=.FALSE.
      IF (PRINT.AND.(DABS(A).GT.TOL
     :           .OR.DABS(B).GT.TOL
     :           .OR.DABS(C).GT.TOL)) THEN
          DOPRINT=.TRUE.
      ENDIF
      RETURN
      END
*--------------------------------------------------------------
*     F I E L D S H
*--------------------------------------------------------------
*
*     This subroutine computes the field shift contribution
*
*
      SUBROUTINE FIELDSH(PRINT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=100,IDIM=550,NCDIM=3000)
      PARAMETER (NELEMENT=50)
*
      CHARACTER CONFIG*40,EL*3,ATOM*6,TERM*6,COUPLE*3
      COMMON /LABEL/CONFIG(NCD),EL(NWD),ATOM,TERM,COUPLE(NCD,9)
*
      INTEGER IN,OUT,ERR,PRI
      COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,IUL,IUJ,IU(2)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,
     :ID,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
      INTEGER KVAL,IEL,CPTR,IH,JH,OPTR
      COMMON /STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM),
     :VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD),
     :X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
      LOGICAL PRINT,CORE,DOPRINT,YM,YF
      COMMON /COR/COREF,COREG,CORER,CORES,ZM1,ZM2,FZE,FZA,DR2,CORE,YM,YF
*
      COMMON /RYDBERG/ RYDBRG,ZMU,ENERGY
      DATA FVP/5.003461D-6/,FSP/8.339102D-6/
      CHARACTER YES
*
*
* ****  WRITE OUT HEADER
*
      WRITE (OUT,101)
*
* ****   1) CLOSED SHELL CONTRIBUTION
*
      IF (.NOT.CORE) GOTO 30
      IF (PRINT) WRITE(OUT,'(/A/)')'  i)   Core'
*
      DO 20 I=1,NCLOSD
         IF (L(I).EQ.0) THEN
            CFC=2*AZ(I)**2
            COREF=COREF+CFC
            IF(DOPRINT(CFC,D0,D0,PRINT,TOL))
     :      WRITE(OUT,102) EL(I),EL(I),CFC
         ENDIF
 20   CONTINUE
 30   WRITE(OUT,103)'Total for core',COREF
*
* ****   c) CALCULATE THE MATRIX ELEMENTS BETEEN OUTER
*           ELECTRONS
*
      IF (PRINT) WRITE(OUT,'(/A/)')' ii)   Outer Electrons'
*
* ****      Evaluate the overlap integrals
*           needed for determining coef(i)
*
      IBEGIN=INTPTR(2)+1
      IEND=INTPTR(3)
      DO 60 I=IBEGIN,IEND
         VALUE(I)=QUADR(IEL(I,1),IEL(I,2),0)**KVAL(I)
 60   CONTINUE
*
      IBEGIN=IEND+1
      IEND=INTPTR(4)
      DO 70 I=IBEGIN,IEND
         K1=KVAL(I)/64
         K2=KVAL(I)-64*K1
         VALUE(I)=QUADR(IEL(I,1),IEL(I,2),0)**K1
     :           *QUADR(IEL(I,3),IEL(I,4),0)**K2
 70   CONTINUE
*
      IBEGIN=INTPTR(5)+1
      IEND = INTPTR(6)
      CFO=D0
      OUTER=D0
      DO 80 I = IBEGIN,IEND
         C = COEF(I)
         IF(L(IEL(I,1)).EQ.0) THEN
            C2=-D2*C
            CFO=C2*AZ(IEL(I,1))*AZ(IEL(I,2))
            OUTER=OUTER+CFO
            IF(DOPRINT(CFO,D0,D0,PRINT,TOL))
     :      WRITE(OUT,102) EL(IEL(i,1)),EL(IEL(I,2)),CFO
         ENDIF
 80   CONTINUE
      WRITE(OUT,103)'Total for outer electrons ',OUTER
      WRITE(OUT,*)
      WRITE(OUT,103)'Total contribution (au)',COREF+OUTER
*
*
* **** CALCULATE THE FINITE VOLUME CORRECTION
*
*
      A1 = 1.115*(ZMU**(1./3.))
      A2 = 2.151*(ZMU**(-1./3.))
      A3 = 1.742/ZMU
      REQ = A1 + A2 - A3
      REQ2 = REQ*REQ
      DENS = COREF + OUTER
      IF (INT(Z).LT.10) THEN
          FVC = FVP*DENS*FZE*REQ2*Z
      ELSE
          FVC = FVP*DENS*FZE*REQ2/Z
      ENDIF
      WRITE(OUT,108) ZMU
      WRITE(OUT,110) FVC
      WRITE(OUT,111) REQ
*
*
* **** CALCULATE THE FIELD SHIFT CONTRIBUTION
*
*
      IF (CORE) THEN
         IF(INT(Z).LT.10) WRITE(ERR,'(A)') ' Warning Z below 10 !'
         WRITE(ERR,'(A)')
     : ' Contribution to the isotope field shift (y/n)'
         READ(IN,'(A)') YES
         IF (YES.EQ.'Y'.OR.YES.EQ.'y') THEN
            YF=.TRUE.
            WRITE(ERR,'(A)')
     :           ' Enter f(Z) (MHz/fm^2) and delta <r^2> (fm^2)'
            READ(IN,*) FZA,DR2
         ENDIF
      ENDIF
      IF (YF) THEN
         FSC = FSP*DENS*FZA*DR2/Z
         WRITE(OUT,112) FSC
      ENDIF

 101  FORMAT(///10X,46('-'),/10X,
     :'  ELECTRONIC CONTRIBUTION TO THE FIELD SHIFT  ',/10X,46('-'))
 102  FORMAT(4X,2(A3,1X),34X,F20.8)
 103  FORMAT(/4X,A,T47,F20.8)
 108  FORMAT(//1X,'1  FINITE VOLUME CORRECTION FOR MASS'
     :,23X,F6.2)
 110  FORMAT(/6X,'Volume correction',27X,F16.8,'  cm-1')
 111  FORMAT(6X,'Req was',37X,F16.8,'  fm'/ )
 112  FORMAT(//1X,'2  FIELD SHIFT CONTRIBUTION',22X,F16.8,'  cm-1')
      RETURN
      END
*     ------------------------------------------------------------
*     I N I T M
*     ------------------------------------------------------------
*
*     Determine the Rydberg constant and
*     initialize the radial grid
*
      SUBROUTINE INITM(IPROB)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=100,IDIM=550,NCDIM=3000)
      PARAMETER (NELEMENT=50)
*
      INTEGER OUT,ERR,PRI
      COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,IUL,IUJ,IU(2)
*
      COMMON /PARAM/ H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,
     :ID,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD),
     :X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
      COMMON /RYDBERG/ RYDBRG,ZMU,ENERGY
*
      LOGICAL CORE,YM,YF
      COMMON /COR/COREF,COREG,CORER,CORES,ZM1,ZM2,FZE,FZA,DR2,CORE,YM,YF
*
      CHARACTER YES
      REAL ME,MN(NELEMENT),FZ(NELEMENT-9)
*
      DATA ME/548.579903D-6/,RY/109737.31534/ MN/1,4,7,9,11,12,14,16,
     : 19,20,23,24,27,28,31,32,35,40,39,40,45,48,51,52,55,56,59,58,
     : 63,64,69,74,75,80,79,84,85,88,89,90,93,98,99,102,103,108,
     : 107,114,115,120/
*
      DATA FZ/163.5,198.4,239.6,283.9,332.7,386.1,444.2,507.8,
     : 576.4,650.7,731.2,817.9,911.2,1011.4,1119.3,1234.9,1359.3,
     : 1491.9,1634.4,1786.6,1948.7,2122.7,2308.2,2507.4,2718.5,
     : 2944.3,3187.5,3442.2,3720.3,4015.4,4326.6,4664.8,5022.4,
     : 5405.9,5811.9,6245.0,6706.8,7201.3,7728.7,8296.5,8889.3/
*

*
* ****   DETERMINE THE RYDBERG CONSTANT AND F(Z) PARAMETER 
*
      WRITE(ERR,'(/1X,A)') 'Default Rydberg constant (y/n)'
      READ(IN,'(A1)') YES
      IF (YES.EQ.'Y'.OR.YES.EQ.'y'.AND.(INT(Z).LE.NELEMENT)) THEN
         ZMU=MN(INT(Z))
      ELSE
         WRITE(ERR,'(1X,A)') 'Enter the mass of the atom'
         READ(IN,*) ZMU
      END IF
      RYDBRG = RY/(1.+ME/ZMU)
      IF (IPROB.NE.1) THEN
         IF (INT(Z).LT.10) THEN
            FZE=1.566432
         ELSEIF (INT(Z).LE.NELEMENT) THEN
            FZE=FZ(INT(Z-9))
         ELSE
            WRITE(ERR,'(1X,A)') 'Enter the f(Z) value in MHz'
            READ(IN,*) FZE
         ENDIF
       ENDIF
*
* ****  INITIALIZE THE RADIAL GRID
*
      RHO=-4.0
      DO 10 J=1,NO
         R(J)=DEXP(RHO)/Z
         RR(J)=R(J)*R(J)
         R2(J)=DSQRT(R(J))
         RHO=RHO+H
10    CONTINUE
      RETURN
      END
*
*-----------------------------------------------------------------------
*               I N T G R L
*-----------------------------------------------------------------------
*
*    Read the integrals that define the energy expression
*
      SUBROUTINE INTGRL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=100)
*
      INTEGER IN,OUT,ERR,PRI,OUC,OUD,OUF,OUH
      COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,OUC,OUD,OUF,OUH
*
      CHARACTER CONFIG*40,EL*3,ATOM*6,TERM*6,COUPLE*3
      COMMON /LABEL/CONFIG(NCD),EL(NWD),ATOM,TERM,COUPLE(NCD,9)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
      INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
      PARAMETER (IDIM=550,NCDIM=3000)
      COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
      CHARACTER END*1, EL1*3, EL2*3, EL3*3, EL4*3
*
    1 FORMAT(1X,A1,I2,1X,A3,1X,A3,1X,I5)
    2 FORMAT(1X,A1,I2,1X,2A3,1X,2A3,1X,I5)
    3 FORMAT(1X,A1,I2,1X,A3,1X,A3,2X,I2,1X,A3,1X,A3,I5)
    4 FORMAT(F14.8,A1,3I3)
*
* ***** READ  THE LIST OF INTEGRALS
*
        LAST = 0
        IC = 1
        I = 1
        READ(IUD,'()')
        DO 10 INT = 1,6
          IF (INT.NE.4 .AND. INT.NE.5) THEN
*
*            ...F, G, L, or O1 integrals....
*
   12        READ(IUD,1) END, KVAL(I), EL1, EL2, ICPTR
             IF (END .EQ. '*') GO TO 16
             IF (ICPTR+LAST .LE. (NCDIM)) THEN
                   CPTR(I) = ICPTR + LAST
             ELSE
                WRITE(ERR,*) ' Too much data-current dimensions:' ,NCDIM
                STOP
             END IF
             CALL EPTR(EL, EL1,IEL(I,1),*999)
             CALL EPTR(EL, EL2,IEL(I,2),*999)
             I = I + 1
             IF (I .LE. (IDIM) ) GO TO 12
             WRITE(ERR,*) ' Too many integrals - MAX =',IDIM
             STOP
           ELSE
   14         IF (INT.EQ.5) THEN
*
*               ... R integrals ...
*
                READ(IUD,2) END, KVAL(I), EL1, EL2, EL3, EL4, ICPTR
*
              ELSE
*
*                ... O2 integrals ...
*
                READ(IUD, 3) END, K1, EL1, EL2, K2, EL3, EL4
                KVAL(I) = 64*K1 + K2
              END IF
              IF (ICPTR+LAST .LE. (NCDIM)) THEN
                   CPTR(I) = ICPTR + LAST
              ELSE
                   STOP ' Too much data - current dimensions = (NCDIM)'
              END IF
*
             IF ( END .EQ. '*') GO TO 16
             CALL EPTR(EL, EL1, IEL(I,1), *999)
             CALL EPTR(EL, EL2, IEL(I,2), *999)
             CALL EPTR(EL, EL3, IEL(I,3), *999)
             CALL EPTR(EL, EL4, IEL(I,4), *999)
             I = I + 1
             IF (I .LE. (IDIM) ) GO TO 14
             STOP ' Too many integrals - MAX = (IDIM)'
          END IF
 16       IF (INT .EQ. 3 .OR. INT .EQ. 4) GO TO 18
*
*       ... Read the data ...
*
   20     READ(IUD,4) COEFF(IC), END, IH(IC), JH(IC), OPTR(IC)
          IF ( END .NE. '*') THEN
            IF (INT .LE. 2) THEN
              COEFF(IC) = ACURAT(COEFF(IC))
            ELSE
*
*         ... Shift origin for overlap integrals
*
              IF (OPTR(IC).LT.0) THEN
                OPTR(IC) = INTPTR(3) - OPTR(IC)
              ELSE IF (OPTR(IC).GT.0) THEN
                OPTR(IC) = INTPTR(2) + OPTR(IC)
              END IF
            END IF
            IC = IC + 1
            GO TO 20
          END IF
*
*       ... Initialize for next set ..
*
   18     INTPTR(INT) = I-1
          LAST = IC-1
   10   CONTINUE
        RETURN
*
  999   WRITE(ERR,*)' Electron in ',END,'-data not found in ',
     :          'configuration list data'
        STOP
        END
*-----------------------------------------------------------------
*     I S D A T A
*-----------------------------------------------------------------
*
*     Read the wavefunctions and configurations
*     from the <name>.w and <name>.c files.
*
      SUBROUTINE ISDATA(ICASE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=100,IDIM=550,NCDIM=3000)
*
      CHARACTER CONFIG*40,EL*3,ATOM*6,TERM*6,COUPLE*3
      COMMON /LABEL/CONFIG(NCD),EL(NWD),ATOM,TERM,COUPLE(NCD,9)
*
      INTEGER OUT,ERR,PRI
      COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,IUL,IUJ,IU(2)
*
      COMMON /PARAM/ H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,
     :ID,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD),
     :X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
      INTEGER KVAL,IEL,CPTR,IH,JH,OPTR
      COMMON /STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM),
     :VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
      COMMON /RYDBERG/ RYDBRG,ZMU,ENERGY
*
      INTEGER QC(5)
      CHARACTER*3 ELC(5),ELCSD(18),COUPL(9)
      CHARACTER LABEL*36,PCONFIG*66
      DATA LABEL/'   Atom  State      Energy (in a.u.)'/
*
*  *****  READ THE RADIAL FUNCTIONS
*
      I=1
      M=1
12    READ(IUF,END=13)ATOM,TERM,EL(I),MR,Z,ETI,EKI,AZ(I),(P(J,I),J=1,MR)
      IF (EL(I)(1:1).NE.' ') THEN
         N(I)=ICHAR(EL(I)(1:1))-ICHAR('1')
         L(I)=LVAL(EL(I)(2:2))
      ELSE
         N(I)=ICHAR(EL(I)(2:2))-ICHAR('1')
         L(I)=LVAL(EL(I)(3:3))
      ENDIF
      MM=MR+1
      DO 24 J=MM,NOD
         P(J,I)=D0
24    CONTINUE
8     MAX(I)=MR
      I=I+1
      IF (I.GT.NWD) STOP 'TOO MANY ELECTRONS:MAX=(NWD)'
      GOTO 12
13    NWF=I-1
*
*  *****  READ THE CONFIGURATIONS
*
      READ(IUC,108,END=1000) ATOM,TERM,ENERGY
      IF (ICASE.EQ.1) THEN
         WRITE(OUT,102) LABEL,ATOM,TERM,ENERGY
         WRITE(OUT,107)
      ELSE
         WRITE(OUT,106)
      ENDIF
      READ(IUC,'(18(1X,A3))')(ELCSD(II),II=1,18)
      NCLOSD=0
6     IF (ELCSD(NCLOSD+1).EQ.EL(NCLOSD+1)) THEN
         NCLOSD=NCLOSD+1
         IF (NCLOSD.LT.18) GOTO 6
      ELSE IF (ELCSD(NCLOSD+1).NE.'   ') THEN
         WRITE(ERR,'(1X,A3,A)')  ELCSD(NCLOSD+1),
     :     ' Does not match with the electron list'
         STOP
      END IF
15    READ(IUC,103,END=999) (ELC(K),QC(K),K=1,5),WT(M)
      READ(IUC,104,END=999) (COUPL(J),J=1,9)
      NOCC=0
16    IF (ELC(NOCC+1) .NE. '   ' ) THEN
         NOCC=NOCC+1
         IF (NOCC.LT.5) GOTO 16
      ENDIF
      IF (NOCC.EQ.0) GOTO 18
      CALL PACK(NOCC,ELC,QC,COUPL,PCONFIG)
      K=66
18    IF (PCONFIG(K:K).EQ.' ') THEN
         K=K-1
         GOTO 18
      ENDIF
      IF (ICASE.EQ.1) THEN
         WRITE(OUT,105) M,PCONFIG(1:K),WT(M)
      ELSE
         WRITE(OUT,105) M,PCONFIG(1:K)
      ENDIF
      M=M+1
      IF (M.GT.(NCD)) THEN
         WRITE(ERR,*)
     :       ' TOO MANY CONFIGURATIONS IN THE INITIAL AND FINAL STATE',
     :       ' MAXIMUM FOR COMBINED SUM = ',(NCD)-1
      ENDIF
      GOTO 15
999   NCFG=M-1
102   FORMAT(/A/3X,31('-')/3X,2A6,F14.7/)
103   FORMAT(5(1X,A3,1X,I2,1X),F10.7)
104   FORMAT(9(5X,A3))
105   FORMAT(I4,3X,A,T60,F10.7)
106   FORMAT(///7X,'CONFIGURATIONS',/)
107   FORMAT(///7X,'CONFIGURATION',T61,'WEIGHT',/)
108   FORMAT(3X,2A6,F14.7)
1000  RETURN
      END
*-------------------------------------------------------------
*     S M A S S
*--------------------------------------------------------------
*
*     This subroutine computes the normal and specific mass
*     shift for a given set of configuration state functions.
*
      SUBROUTINE SMASS(PRINT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=100,IDIM=550,NCDIM=3000)
*
      CHARACTER CONFIG*40,EL*3,ATOM*6,TERM*6,COUPLE*3,ELC(4)*3,YES
      COMMON /LABEL/CONFIG(NCD),EL(NWD),ATOM,TERM,COUPLE(NCD,9)
*
      INTEGER IN,OUT,ERR,PRI
      COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,IUL,IUJ,IU(2)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,
     :ID,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
      INTEGER KVAL,IEL,CPTR,IH,JH,OPTR
      COMMON /STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM),
     :VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD),
     :X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
      LOGICAL PRINT,REL,DOPRINT,CORE,YM,YF
      COMMON /COR/COREF,COREG,CORER,CORES,ZM1,ZM2,FZE,FZA,DR2,CORE,YM,YF
*
      COMMON /RYDBERG/ RYDBRG,ZMU,ENERGY
      LOGICAL YI
      REAL ME
      DATA ME/548.579903D-6/,RY/109737.31534/,REL/.FALSE./
*
* **** INITIALIZE
*
      GRADM2=D0
      RMASS2=D0
      SLATM2=D0
      GRADM3=D0
      RMASS3=D0
      SLATM3=D0
*
* ****  WRITE OUT HEADER (OUTPUT FILE)
*
      WRITE (OUT,101)
      WRITE (OUT,1)
 1    FORMAT(//41X,'CONTRIBUTIONS TO SPECIFIC MASS SHIFT',/)
      WRITE(OUT,2)
 2    FORMAT(38X,'Gradient Form',7X,'Slater Integral Form',/38X
     :,13('-'),7X,20('-'))
      WRITE(OUT,3)
 3    FORMAT(41X,'TOTAL',12X,'PART1',11X,'PART2')
*
* ****   1) CLOSED SHELL CONTRIBUTION
*
      IF (PRINT.AND.CORE) WRITE(OUT,'(/A/)')' i)   Core'
*
* ****   a) CALCULATE THE MATRIX ELEMENTS BETWEEN
* ****      ELECTRONS BELONGING TO A COMMON CORE
*
      IF (.NOT.CORE) GOTO 50
      DO 20 I=1,NCLOSD
         SLATM2=D0
         GRADM2=D0
         SUMI=4*L(I)+2
         RMASS=FK(I,I,0,REL)
         DO 30 K=2,2*L(I),2
            RMASS=RMASS-CA(L(I),K)*FK(I,I,K,REL)
  30     CONTINUE
         RMASS2=SUMI*(SUMI-D1)*RMASS/D4
         CORER=CORER+RMASS2
         IF (DOPRINT(GRADM2,RMASS2,SLATM2,PRINT,TOL)) THEN
            WRITE(OUT,102) EL(I),EL(I),GRADM2,RMASS2,SLATM2
         ENDIF
*
* ****   b) CALCULATE THE MATRIX ELEMENTS BETWEEN
*           ELECTRONS BELONGING TO DIFFERENT CORES
*
         DO 20 J=1,I-1
            RMASS=FK(I,J,0,REL)
            SLATM1=D0
            GRADM1=D0
            SUMJ=4*L(J)+2
            DO 40 K=IABS(L(I)-L(J)),L(I)+L(J),2
               C=CB(L(I),L(J),K)
               RMASS=RMASS-C*GK(I,J,K,REL)
               IF (K.EQ.1) THEN
                  SLATM1=-C*Z*QUADR(I,J,-2)*QUADR(I,J,1)
                  GRADM1=-C*GRAD(I,J)**2
               ENDIF
  40        CONTINUE
            RMASS2=SUMI*SUMJ*RMASS/D2
            SLATM2=SUMI*SUMJ*SLATM1
            GRADM2=SUMI*SUMJ*GRADM1
            CORER=CORER+RMASS2
            CORES=CORES+SLATM2
            COREG=COREG+GRADM2
            IF (DOPRINT(GRADM2,RMASS2,SLATM2,PRINT,TOL)) THEN
               WRITE(OUT,102) EL(I),EL(J),GRADM2,RMASS2,SLATM2
            ENDIF
  20  CONTINUE
  50  WRITE(6,103)'Total for core',COREG,CORER,CORES
      RMASS3=CORER
      SLATM3=CORES
      GRADM3=COREG
*
* ****   c) CALCULATE THE MATRIX ELEMENTS BETEEN OUTER
*           ELECTRONS BELONGING TO A COMMON CONFIGURATION
*
      INTG=1
      IBEGIN=1
 80   IF (PRINT.AND.INTG.EQ.1) THEN
         WRITE(OUT,'(//A)') ' ii)  Fk-integrals'
      ELSE IF (PRINT.AND.INTG.EQ.2) THEN
         WRITE(OUT,'(//A)') ' iii) Gk-integrals'
      ENDIF
      IEND = INTPTR(INTG)
      DO 70 I=IBEGIN,IEND
         SLATM1=D0
         GRADM1=D0
         C=COEF(I)
         ELC(1)=EL(IEL(I,1))
         ELC(2)=EL(IEL(I,2))
         IF (INTG.EQ.1) THEN
            ELC(3)=ELC(1)
            ELC(4)=ELC(2)
            RMASS=FK(IEL(I,1),IEL(I,2),KVAL(I),REL)/D2
         ELSE
            RMASS=GK(IEL(I,1),IEL(I,2),KVAL(I),REL)/D2
            ELC(3)=ELC(2)
            ELC(4)=ELC(1)
            IF (KVAL(I).EQ.1) THEN
               SLATM1=Z*QUADR(IEL(I,1),IEL(I,2),1)*
     :                 QUADR(IEL(I,2),IEL(I,1),-2)
            ENDIF
         ENDIF
         IF (KVAL(I).EQ.1) THEN
            GRADM1=GRAD(IEL(I,1),IEL(I,2))**2
         ENDIF
         RMASS2=C*RMASS
         SLATM2=C*SLATM1
         GRADM2=C*GRADM1
         RMASS3=RMASS3+RMASS2
         SLATM3=SLATM3+SLATM2
         GRADM3=GRADM3+GRADM2
         IF (DOPRINT(GRADM2,RMASS2,SLATM2,PRINT,TOL)) THEN
            WRITE(OUT,104)ELC(1),ELC(2),ELC(3),ELC(4),
     :                           GRADM2,RMASS2,SLATM2
         ENDIF
 70   CONTINUE
      IBEGIN=IEND+1
      INTG=INTG+1
      IF (INTG.LE.2) GOTO 80
*
* ****   d) CALCULATE THE MATRIX ELEMENTS BETWEEN THE OUTER
*           ELECTRONS BELONGING TO DIFFERENT CONFIGURATIONS
*           Rk integrals
*
*           Evaluate the overlapp integrals
*           needed for determining coef(i)
*
      IBEGIN=INTPTR(2)+1
      IEND=INTPTR(3)
      DO 90 I=IBEGIN,IEND
         VALUE(I)=QUADR(IEL(I,1),IEL(I,2),0)**KVAL(I)
 90   CONTINUE
*
      IBEGIN=IEND+1
      IEND=INTPTR(4)
      DO 100 I=IBEGIN,IEND
         K1=KVAL(I)/64
         K2=KVAL(I)-64*K1
         VALUE(I)=QUADR(IEL(I,1),IEL(I,2),0)**K1
     :           *QUADR(IEL(I,3),IEL(I,4),0)**K2
 100  CONTINUE
*
      IF (PRINT) THEN
         WRITE(OUT,'(//A)') ' iv)  Rk-integrals'
      ENDIF
      IBEGIN=IEND+1
      IEND=INTPTR(5)
      DO 120 I=IBEGIN,IEND
         SLATM1=D0
         GRADM1=D0
         C=COEF(I)
         ELC(1) = EL(IEL(I,1))
         ELC(2) = EL(IEL(I,2))
         ELC(3) = EL(IEL(I,3))
         ELC(4) = EL(IEL(I,4))
         RMASS=D5*RK(IEL(I,1),IEL(I,2),IEL(I,3),IEL(I,4),KVAL(I),REL)
         IF (KVAL(I).EQ.1) THEN
            SLATM1=(Z/D2)*(QUADR(IEL(I,1),IEL(I,3),1)*
     :                    QUADR(IEL(I,2),IEL(I,4),-2)+
     :                    QUADR(IEL(I,1),IEL(I,3),-2)*
     :                    QUADR(IEL(I,2),IEL(I,4),1))
            GRADM1=-GRAD(IEL(I,1),IEL(I,3))*GRAD(IEL(I,2),IEL(I,4))
         ENDIF
         RMASS2=C*RMASS
         SLATM2=C*SLATM1
         GRADM2=C*GRADM1
         RMASS3=RMASS3+RMASS2
         SLATM3=SLATM3+SLATM2
         GRADM3=GRADM3+GRADM2
         IF (DOPRINT(GRADM2,RMASS2,SLATM2,PRINT,TOL)) THEN
            WRITE(OUT,104) ELC(1),ELC(2),ELC(3),ELC(4),
     :                            GRADM2,RMASS2,SLATM2
         ENDIF
 120  CONTINUE
      AM=GRADM3-COREG
      BM=RMASS3-CORER
      CM=SLATM3-CORES
      WRITE(OUT,103)'Total for integrals',AM,BM,CM
      AM=GRADM3
      BM=RMASS3
      CM=SLATM3
*
* ****   e) CALCULATE THE MATRIX ELEMENTS BETEEN THE
*           OUTER ELECTRONS AND THE CORE ELECTRONS
*
      IF (PRINT) WRITE(OUT,'(//A)') ' v)   Core-Outer'
      IBEGIN=IEND+1
      IEND=INTPTR(6)
      DO 130 I=IBEGIN,IEND
         SUMI=-D2*COEF(I)
         DO 130 J=1,NCLOSD
            RMASS=RK(IEL(I,1),J,IEL(I,2),J,0,REL)
            SLATM1=D0
            GRADM1=D0
            SUMJ=4*L(J)+2
            DO 140 K=IABS(L(IEL(I,1))-L(J)),L(IEL(I,1))+L(J),2
               C=CB(L(IEL(I,1)),L(J),K)
               RMASS=RMASS-C*RK(IEL(I,1),J,J,IEL(I,2),K,REL)
               IF (K.EQ.1) THEN
                  IF (IEL(I,1).EQ.IEL(I,2)) THEN
                     SLATM1=-C*Z*QUADR(IEL(I,1),J,1)*
     :                          QUADR(J,IEL(I,2),-2)
                     GRADM1=-C*GRAD(IEL(I,1),J)**2
                  ELSE
                     SLATM1=-(Z/D2)*C*(QUADR(IEL(I,1),J,1)*
     :                                QUADR(J,IEL(I,2),-2)+
     :                                 QUADR(IEL(I,2),J,1)*
     :                                QUADR(J,IEL(I,1),-2))
                     GRADM1=C*GRAD(IEL(I,1),J)*GRAD(J,IEL(I,2))
                  ENDIF
               ENDIF
 140        CONTINUE
            RMASS2=SUMI*SUMJ*RMASS/D2
            GRADM2=SUMI*SUMJ*GRADM1
            SLATM2=SUMI*SUMJ*SLATM1
            GRADM3=GRADM3+GRADM2
            RMASS3=RMASS3+RMASS2
            SLATM3=SLATM3+SLATM2
            IF (DOPRINT(GRADM2,RMASS2,SLATM2,PRINT,TOL)) THEN
                ELC(1)=EL(IEL(I,1))
                ELC(2)=EL(IEL(I,2))
                WRITE(OUT,104) ELC(1),EL(J),ELC(2),EL(J),
     :                               GRADM2,RMASS2,SLATM2
            ENDIF
 130  CONTINUE
      AM=GRADM3-AM
      BM=RMASS3-BM
      CM=SLATM3-CM
      WRITE(OUT,103)'Total for core-outer',AM,BM,CM
*
* ****  SUMMARY
*       CALCULATE THE TOTAL CONTRIBUTION TO THE ISOTOPE SHIFT
*
      RATIO=D2*RY*ME/ZMU
      SLATM3=SLATM3+RMASS3
      TOTSLAT=SLATM3*RATIO
      TOTGRAD=GRADM3*RATIO
      WRITE(OUT,106) RMASS3,SLATM3-RMASS3
      WRITE(OUT,107) GRADM3,SLATM3
      WRITE(OUT,108) ZMU
      WRITE(OUT,109) 'Gradient Form       ',TOTGRAD
      WRITE(OUT,109) 'Slater Integral Form',TOTSLAT
*
      WRITE(OUT,110) 'Rydberg constant was',RYDBRG
      DE1=-ENERGY*2*RYDBRG*ME/ZMU
      WRITE(OUT,112)DE1
      IF (CORE) THEN
         WRITE(ERR,'(/1X,A)')
     :  'Contribution to the isotope mass shift (y/n)'
         READ(IN,'(A)') YES
         IF (YES.EQ.'Y'.OR.YES.EQ.'y') THEN
             YI=.TRUE.
             WRITE(ERR,'(A)')' Enter two masses'
             READ(IN,*) ZM1,ZM2
         ENDIF
      ENDIF
      IF (YI) THEN
         WRITE(OUT,114)' 3  ISOTOPE SHIFT FOR MASSES',ZM1,ZM2
         WRITE(OUT,'(A)')'    a)  SPECIFIC MASS CONTRIBUTION'
         RATIO = (ZM2-ZM1)/((ZM1+ME)*(ZM2+ME))
         C=D2*RY*ME
         DE1=ENERGY*2*RY*RATIO*ME
         RATIO = (ZM2-ZM1)/(ZM1*ZM2)
         TOTGRAD=C*GRADM3*RATIO
         TOTSLAT=C*SLATM3*RATIO
         WRITE(OUT,109) 'Gradient Form       ',TOTGRAD
         WRITE(OUT,109) 'Slater Integral Form',TOTSLAT
         WRITE(OUT,116) DE1
      ENDIF
*
101   FORMAT(///20X,21('-'),/20X,'  MASS CONTRIBUTION  '/20X,21('-'))
102   FORMAT(4X,2(A3,1X),22X,F14.8,4X,2(F14.8,2X))
103   FORMAT(/4X,A,T35,F14.8,4X,2(F14.8,2X))
104   FORMAT(4X,4(A3,1X),14X,F14.8,4X,2(F14.8,2X))
106   FORMAT(/33X,49('-'),//4X,'Total contribution',30X,2(F14.8,2X))
107   FORMAT(4X,'in atomic units:',14X,F14.8,13X,F14.8)
108   FORMAT(//1X,'1  SPECIFIC MASS SHIFT CORRECTION FOR MASS'
     :,17X,F6.2)
109   FORMAT(6X,A20,26X,F14.8,'  cm-1')
110   FORMAT(6X,A20,26X,F14.4,'  cm-1')
112   FORMAT(/1X,'2  NORMAL MASS SHIFT CORRECTION  ',
     :18X,F14.8,'  cm-1')
114   FORMAT(//A,23X,F6.2,' - ',F6.2,/)
116   FORMAT('    b)  NORMAL MASS CONTRIBUTION  ',
     :18X,F14.8,'  cm-1')
      RETURN
      END
