*     ------------------------------------------------------------------
*        MULTICONFIGURATION HARTREE-FOCK PROGRAM
*
*                 C O P Y R I G H T -- 1994
*
*     Written and Revised by:  Charlotte Froese Fischer
*                              Department of Computer Science
*                              Vanderbilt University
*     September, 1988
*
*     Computer Physics Communication, Vol. 64, 431 (1991)
*     ------------------------------------------------------------------
*
*     All comments in the program assume the radial function P
*     is the solution of an equation of the form
*
*      P" + ( 2Z/R - Y - L(L+1)/R**2 - E)P = X + T
*
*     where Y is called a potential function
*           X is called an exchange function, and
*           T includes contributions from off-diagonal energy parameter,
*             interactions between configurations, etc.
*
*     The program uses LOG(Z*R) as independent variable and
*                      P/SQRT(R) as dependent variable.
*     As a result all equations must be transformed to these variables
*
*     For details of the numerical procedures, see
*     Computer Physics Reports, Vol. 3, 273--326 (1986).
*     ------------------------------------------------------------------
*           M A I N   P R O G R A M
*     ------------------------------------------------------------------
*
*       The MAIN program controls the overall calculation and  allows
*   a series of cases to be processed as one run.  Each  case  itself
*   may  consist  of  a  series of atoms or ions in an iso-electronic
*   sequence.  In each case, all but the initial  estimates  for  the
*   first  are  obtained  by  scaling  the previous results.
*   Mixing coefficients are left unchanged.
*
        PROGRAM MCHF
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220,NCD=100)
*
        INTEGER IN,OUT,ERR,PRI,OUC,OUD,OUF,OUH,IOU(7)
        COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,OUC,OUD,OUF,OUH
*
        CHARACTER CONFIG*40,EL*3,ATOM*6,TERM*6,COUPLE*3
        COMMON /LABEL/CONFIG(NCD),EL(NWD),ATOM,TERM,COUPLE(NCD,9)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
        LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :     ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
        LOGICAL PRINT,LD
        CHARACTER*24 ANS*1,NAME(7), FILE
CSUN  REAL TIMES(2),DTIME
      EQUIVALENCE (IUC,IOU(1)),(OUC,IOU(4))
      DATA NAME/'cfg.inp','int.lst','wfn.inp','cfg.out',' ',
     :          'wfn.out',' '/
*
CSUN     RTIME = DTIME(TIMES)
*  ***** Define unit numbers and open files *********************
*                                                               *
*        UNIT NUMBERS AND FILE NAMES MAY BE MACHINE             *
*        DEPENDENT. CHECK THE FOLLOWING SECTION.                *
*                                                               *
*        IN - Standard input unit, normally the terminal        *
*        OUT- Standard output unit, normally the terminal       *
*        ERR- Prompts and Error messages, always the terminal   *
*        PRI- Printer output unit or file.                      *
*                                                               *
        IN = 5
        OUT = 6
        ERR = 0
        PRI = 3
*
*  *****  WRITE OUT HEADER
*
      WRITE(OUT,9)
9     FORMAT(//20X,'==================='/
     :         20X,'  MCHF  ...  1988'/
     :         20X,'===================')
*
*  *****  WRITE OUT DIMENSION INFORMATION
*
      WRITE(OUT,99) 'NCFG',NCD,'NWF',NWD,'NO',NOD
99    FORMAT(//10X,'THE DIMENSIONS FOR THE CURRENT VERSION ARE:'/
     :       (10X,3(A6,'=',I3,4X)/)/)
*
*  *****  INITIALIZE COMMON DATA ARRAYS
*
      CALL INITA
      CALL INITR
*                                                               *
*  ***** IN THE OPEN STATEMENTS CHECK FOR VALID FILE NAMES ******
*                                                               *
1     WRITE(ERR,'(//A/A//)') ' START OF CASE',' ============='
CSUN       i = iargc()
CSUN       do 999 j = 1,i
CSUN          call getarg(j,FILE)
CSUN          jj = ichar(FILE(1:1)) - ichar('0')
CSUN          name(jj) = FILE(3:)
CSUN  999   continue
      OPEN(UNIT=PRI,FILE='summry',STATUS='UNKNOWN')
      DO 37 J = 1,7
         IOU(J) = 0
         IF (J.LE.3) THEN
            IF (NAME(J) .NE. '  ') THEN
               IF (J.NE. 3) THEN
                  IOU(J) = 20+J
                  OPEN(UNIT=IOU(J),FILE=NAME(J),STATUS='OLD')
               ELSE
                  INQUIRE(FILE=NAME(J),EXIST=LD)
                  IF (LD) THEN
                     IOU(J) = 20+J
                     OPEN(UNIT=IOU(J),FILE=NAME(J),STATUS='OLD',
     :                    FORM = 'UNFORMATTED')
                  END IF
               END IF
            END IF
         ELSE
            IF (NAME(J) .NE. '  ') THEN
               IOU(J) = 20+J
               IF (J.NE.6) THEN
                  OPEN(UNIT=IOU(J),FILE=NAME(J),STATUS='UNKNOWN')
               ELSE
                  OPEN(UNIT=IOU(J),FILE=NAME(J),STATUS='UNKNOWN',
     :                 FORM='UNFORMATTED')
               END IF
            END IF
         END IF
37    CONTINUE
      OPEN(UNIT=35,STATUS='SCRATCH',FORM='UNFORMATTED')
*
*  ***** END OF INPUT/OUTPUT INTERFACE **************************
*
      FAIL = .FALSE.
      DO 4 I=1,(NWD)
      DPM(I) = D10
      IEPTR(I) = 0

4     CONTINUE
      DO 5 I = 1,(98)
         IJE(I) = 0
5     CONTINUE
*
*  *****  DETERMINE DATA ABOUT THE PROBLEM
*
      CALL DATA
*
*  *****  SET PARAMETERS TO THEIR DEFAULT VALUE
*
13    PRINT = .FALSE.
      CFGTOL = 1.D-10
      SCFTOL = 1.D-7
      NSCF = 12
      IC = 0
      IF (NCFG .EQ. 1) IC = 3 + (NWF+1-IB)/4
      ACFG = D0
      LD = .TRUE.
      TRACE = .FALSE.
      WRITE(ERR,'(/A)') ' Default values for other parameters ? (Y/N) '
      READ (IN,'(A)') ANS
      IF (ANS .EQ. 'N' .OR. ANS .EQ. 'n') THEN
         WRITE(ERR,'(/A,A)') ' Default values for',
     :   ' PRINT, CFGTOL, SCFTOL ? (Y/N) '
         READ(IN,'(A)') ANS
         IF ( ANS .NE. 'Y' .AND. ANS .NE. 'y'  ) THEN
            WRITE(ERR,'(A)') ' Input free FORMAT(L, F, F) '
            READ(IN,*) PRINT, CFGTOL, SCFTOL
         END IF
         WRITE(ERR,'(/A)') ' Default values for NSCF, IC ? (Y/N) '
         READ(IN,'(A)') ANS
         IF (ANS .NE. 'Y' .AND. ANS .NE. 'y' ) THEN
            WRITE(ERR,'(A)') ' Input free FORMAT(I, I) '
            READ(IN,*) NSCF, IC
         END IF
         WRITE(ERR,'(/A)') ' Default values for ACFG,LD,TRACE ? (Y/N) '
         READ(IN,'(A)') ANS
         IF (ANS .NE. 'Y' .AND. ANS .NE. 'y') THEN
             WRITE(ERR,'(A)') ' Input free FORMAT( F, L, L) '
             READ(IN,*) ACFG,LD,TRACE
         END IF
      END IF
      WRITE(OUT,2) PRINT,CFGTOL,SCFTOL,NSCF,IC,ACFG,ID,LD,TRACE
2     FORMAT(/L3,2D6.1,2I3,F3.1,I3,2L3)
*
*  *****  PERFORM THE MCHF ITERATION
*
      CALL SCF(ACFG,SCFTOL,CFGTOL,LD)
*
*  *****  OUTPUT RESULTS IF PRINT = .TRUE.
*
      CALL OUTPUT(PRINT)
15    IF (FAIL) GO TO 6
      CALL SUMMRY
*
*  *****  CHECK FOR ISOELECTRONIC SEQUENCE OR END OF CASE.
*
*      WRITE(ERR,'(/A)') ' Do you wish to continue along the sequence ? '
*      READ(IN,'(A1)') ANS
*      IF (ANS .EQ. 'Y' .OR. ANS .EQ.'y') THEN
*          WRITE(ERR,*) ' ATOM, ZZ, (ACC(I),I=1,NWF) in ',
*     ;             ' format(A6,F6.0,(20F3.1))'
*          READ(IN,'(A6,F6.0,(20F3.1))') ATOM, ZZ, (ACC(I),I=1,NWF)
*          WRITE(OUT,10) ATOM,ZZ,(ACC(I),I=1,NWF)
*10        FORMAT(1X,A6,F6.0,(20F3.1))
**
**  *****  SCALE RESULTS FOR ANOTHER MEMBER OF THE ISOELECTRONIC SEQUENCE
**
*          CALL SCALE(ZZ)
*          WRITE(PRI,14) ATOM,TERM
*14        FORMAT(1H1,9X,2A6)
*          CALL ORTHOG
*          GO TO 13
*      END IF
*
*  *****  DETERMINE END OF CASE
*
6     CONTINUE
CSUN  RTIME = DTIME(TIMES)
CSUN  WRITE(ERR,'(//A/A//A/3F10.3//)') ' END OF CASE',' ===========',
CSUN :  '    Real      User      System  Time (in minutes)',
CSUN : RTIME/60.,TIMES(1)/60., TIMES(2)/60.
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
*
*----------------------------------------------------------------------
*       C O N T C
*----------------------------------------------------------------------
*
*       Contributions to first and second-order energy corrections 
*   from the rotation of orbitals from a given integral
*
        SUBROUTINE CONTC(INT,DC,U,COEF)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION U(*)
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
           T = COEFF(II)
           IF (OPTR(II).NE.0) T = T*VALUE(OPTR(II))
           U(IH(II)) = U(IH(II)) + WT(JH(II))*T*DC
           IF (IH(II) .NE. JH(II)) THEN
                U(JH(II)) = U(JH(II)) + WT(IH(II))*T*DC
                T = T+T
           END IF
           COEF = COEF+T*WT(IH(II))*WT(JH(II))
  1     CONTINUE
        END
*
*-----------------------------------------------------------------------
*       C O N T O V
*-----------------------------------------------------------------------
*     Contribution to the first and second order variation of the
*   energy from the presence of an overlap integral.  The present
*   implementation is incomplete in that it excludes the case
*   where one orbital is in an integral and the other in an overlap.
*
        SUBROUTINE CONTOV(M,DC,U,COV)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION U(*)
*
        INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
        PARAMETER (NCD=100,IDIM=550,NCDIM=3000)
        COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :         ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
        COV = 0.D0
        IBEGIN = INTPTR(4)+1
        IEND =INTPTR(6)
        DO 10 I = IBEGIN,IEND
           JBEGIN = CPTR(I-1)+1
           JEND = CPTR(I)
           DO 12 J = JBEGIN,JEND
              IF ( OPTR(J) .EQ. M) THEN
                 CC =  COEFF(J)*VALUE(I)
                 U(IH(J)) = U(IH(J)) + WT(JH(J))*CC*DC
                 IF (IH(J) .NE. JH(J)) THEN
                    U(JH(J)) = U(JH(J)) + WT(IH(J))*CC*DC
                    CC = 2*CC
                 END IF
                 COV = COV + CC*WT(IH(J))*WT(JH(J))
              END IF
 12        CONTINUE
 10     CONTINUE
        END
*
*-----------------------------------------------------------------------
*       C O V
*-----------------------------------------------------------------------
*
*       The coefficient of a particular overlap integral in the
*   expression for the energy.
*
        DOUBLE PRECISION FUNCTION COV(M)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
        INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
        PARAMETER (NCD=100,IDIM=550,NCDIM=3000)
        COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :         ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
        COV = 0.D0
        IBEGIN = INTPTR(4)+1
        IEND =INTPTR(6)
        DO 10 I = IBEGIN,IEND
           JBEGIN = CPTR(I-1)+1
           JEND = CPTR(I)
           DO 12 J = JBEGIN,JEND
              IF ( OPTR(J) .EQ. M) THEN
                 CC =  COEFF(J)*WT(IH(J))*WT(JH(J))*VALUE(I)
                 IF (IH(J) .NE. JH(J)) CC = 2*CC
                 COV = COV + CC
              END IF
 12        CONTINUE
 10     CONTINUE
        END
*
*     ------------------------------------------------------------------
*              D A T A
*     ------------------------------------------------------------------
*
*       Data concerning the number of configurations (NCFG), the number
*   and type of electrons in each  configuration,  as  well  as  data
*   associated with the energy expression are read and stored.
*
*
      SUBROUTINE DATA
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220,NCD=100,IDIM=550,NCDIM=3000)
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
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
        INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
        COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :         ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
        LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
        COMMON ZZ(NWD),IND(NWD),IELI(5),NOCCSH(NCD)
*
        LOGICAL SETORT,STRONG
        CHARACTER*3 EL1,EL2,ELCLSD(18),ELORT(10,2),
     :              ELI(5),ANS*1,STRING*40,LIST*72
*
    1 FORMAT(18(1X,A3))
    7 FORMAT(A3,F6.0,I3,I3,F3.1)
*
*  *****  READ 'ATOM' CARD
*
5     WRITE(ERR,'(/A)') ' ATOM, TERM, Z in FORMAT(A,A,F) : '
      READ(IN,'(A)') STRING
      I = INDEX(STRING,',')
      IF ( I .EQ. 0) THEN
          WRITE(ERR,*) ' ATOM, TERM, and Z must be separated by commas '
          GO TO 5
      END IF
      ATOM = STRING(1:I-1)
      J = INDEX(STRING(I+1:),',')
      IF ( J .EQ. 0) THEN
          WRITE(ERR,*) ' ATOM, TERM, and Z must be separated by commas '
          GO TO 5
      END IF
      TERM = STRING(I+1:I+J-1)
*     Some compilers do not allow the list directe I/O here.  In such
*     cases use the following, but then Z must have an overriding 
*     decimal point or expressed as a three digit integer with
*     leading blanks, if necessary
*     READ(STRING(I+J+1:),'(F3.0)') Z
      READ(STRING(I+J+1:),*) Z
*
*  ***** READ CONFIGURATION CARDS AND NORMALIZE THE WEIGHTS
*
      IF (IUC .EQ. IN ) THEN
*
*  *****  INPUT COMMON CLOSED SHELLS
*
      WRITE(ERR,*) ' List the common CLOSED shells - FORMAT(18(1X,A3))'
      READ(IUC,1) (ELCLSD(I),I=1,18)
*
*  *****  INPUT THE CONFIGURATIONS
*
      NCFG = 0
      WRITE(ERR,'(/A/A/)') ' Enter configurations followed by weights',
     :   ' Example:  1s(2)2s(2)2p(1),1.0 '
    2 WRITE(ERR,'($,I6,A)') NCFG+1,'.  '
      READ(IUC,'(A)',END=10)  STRING
      IF (STRING(1:1) .NE. '*' .AND. STRING(1:3) .NE. '   ') THEN
         ICOMMA = INDEX(STRING,',')
         IF (ICOMMA .EQ. 0) THEN
            W = D0
            ICOMMA = 41
          ELSE
            W = D0
            IF (ICOMMA .LT. 40) READ(STRING(ICOMMA+1:),'(F10.8)') W
         END IF
         NCFG = NCFG+1
         IF (NCFG .LE. (NCD)) THEN
            WT(NCFG) = W
            CALL REFORM(STRING, CONFIG(NCFG))
            DO 4 I = 1,9
               COUPLE(NCFG,I) = '   '
   4        CONTINUE
            GO TO 2
         ELSE
            WRITE(ERR,*) ' TOO MANY CONFIGURATIONS: MAX =', NCD
         END IF
      END IF
      ELSE
*
*  *****  READ CONFIGURATIONS FROM A FILE
*
         READ(IUC,'(/18(1X,A3))') (ELCLSD(I),I=1,18)
         NCFG = 0
    3    READ(IUC,'(A40,F10.8)',END=10) STRING,W
         IF (STRING(1:1) .NE. '*' .AND. STRING(1:3) .NE. '   ') THEN
            NCFG = NCFG+1
            IF (NCFG .LE. (NCD) ) THEN
               CONFIG(NCFG) = STRING
               WT(NCFG) = W
               READ(IUC,'(9(5X,A3))') (COUPLE(NCFG,J),J=1,9)
               GO TO 3
              ELSE
               WRITE(ERR,*) ' TOO MANY CONFIGURATIONS: MAX =',NCD
            END IF
         END IF
      END IF
10    CONTINUE
*
*  --- The following feature was used at one time to  determine a basis 
*   for a Rydberg series as in 3snd, n = 3, 4, 5, etc. but complicated 
*   orthogonality constraints can invalidate this approach, in general.  
*   Setting ID=1 turns OFF the matrix diagonalization and the changing 
*   of the mixing-coefficients. It also sets off-diagonal energy
*   parameters to zero.
*
      ID = 0
C     IF ( NCFG .GT. 1) THEN
C       WRITE(ERR,'(/A,A)')' Is this an MCHF calculation (Y)',
C    :           ' or a basis calculation (N) ? (Y/N) '
C       READ(IN,'(A)') ANS
C       IF (ANS .EQ. 'N' .OR. ANS .EQ. 'n') ID = 1
C     END IF
*
*  *****  DETERMINE NCLOSD SHELLS
*
      I = 0
      SS = D0
   12 IF (ELCLSD(I+1) .NE. '   ') THEN
         I = I+1
         VARIED(I) = .TRUE.
         EL(I) = ELCLSD(I)
         J = 3
         IF (EL(I)(1:1) .NE. ' ') J = 2
         L(I) = LVAL(EL(I)(J:J))
         N(I) = ICHAR(EL(I)(J-1:J-1)) - ICHAR('1') + 1
         IFULL = 2*(2*L(I)+1)
         S(I) = SS + IFULL/2
         SS = SS + IFULL
         METH(I) = 1
         ACC(I) = D0
         IND(I) = 0
         SUM(I) = 4*L(I)+2
         IF (IUF .NE. 0)  IND(I) = -1
         IF( I .LT. 18) GO TO 12
         WRITE(ERR,*) ' TOO MANY CLOSED SHELLS: MAX = 18'
         STOP 1
      END IF
      NCLOSD = I
*
*  *****  DETERMINE THE OTHER ELECTRONS
*
      MAXORB = NCLOSD
      DO 15 NC = 1,NCFG
         STRING = CONFIG(NC)
         J = 2
         I = 0
 16      IF (STRING(J:J+2) .NE. '   ' ) THEN
*
*  *****     An electron has been found; is it a new one?
*
            I = I+1
            EL1 = STRING(J:J+2)
            K = NCLOSD + 1
 17         IF (K .LE. MAXORB) THEN
               IF ( EL(K) .NE. EL1 ) THEN
                  K = K+1
                  IF (K .GT. (NWD)) THEN
                    WRITE(ERR,*) ' TOO MANY ELECTRONS: MAX=',NWD
                    STOP 1
                  END IF
                  GO TO 17
               END IF
              ELSE
*
*  *****         A new electron has been found; add it to the list
*
               MAXORB = K
               EL(MAXORB) = EL1
            END IF
            J = J+8
            IF (J .LT. 40) GO TO 16
         END IF
         NOCCSH(NC) = I
   15 CONTINUE
      NWF = MAXORB
*
*  *****  The list of electrons has been determined

      WRITE(ERR,19) MAXORB,(EL(J),J=1,MAXORB)
   19 FORMAT(/' There are ',I3,' orbitals as follows:'/(1X,18(1X,A3)))
21    WRITE(ERR,'(A,A)') ' Enter orbitals to be varied:',
     :     ' (ALL,NONE,SOME,NIT=,comma delimited list)'
      READ '(A)', LIST
      IF (LIST(1:3) .EQ. 'ALL' .OR. LIST(1:3) .EQ. 'all') THEN
         NIT = NWF
      ELSE IF (LIST(1:4).EQ.'SOME' .OR. LIST(1:4).EQ.'some') THEN
         NIT = NWF - NCLOSD
      ELSE IF (LIST(1:4).EQ.'NONE' .OR. LIST(1:4).EQ.'none') THEN
         NIT = 0
      ELSE IF (INDEX(LIST,'=') .NE. 0) THEN
         J = INDEX(LIST,'=')
         READ(LIST(J+1:),'(I3)') NIT
221      IF (NIT .gt. NWF) THEN
           NIT = NIT/10
           GO TO 221
         END IF
      ELSE
         NIT = 0
         J = 1
22       NEXT = INDEX(LIST(J:),',')
*
*        ***  Search for last electron label which need not be followed
*             by a comma
*
         IF (NEXT .EQ. 0 .AND. LIST(J:J+2) .NE. '   ')
     :       NEXT = INDEX(LIST(J+1:),' ') + 1
         IF (NEXT .NE. 0) THEN
            IF (NEXT .EQ. 4) THEN
              EL1 = LIST(J:J+2)
            ELSE IF (NEXT .EQ. 3) THEN
              EL1 = ' '//LIST(J:J+1)
            ELSE
              WRITE(ERR,*)'Electron labels must be separated by commas;'
              WRITE(ERR,*)' each label must contain 2 or 3 characters'
              GO TO 21
            END IF
            CALL REORD(EL,EL1,NWF,IERR)
            IF (IERR .EQ. 0) THEN
               NIT = NIT + 1
               J = J + NEXT
               IF (J .LT. 72) GO TO 22
            ELSE
               WRITE(ERR,*) ' Case must match as well as position of',
     :                  ' imbedded blanks'
               WRITE(ERR,*) ' For 3rd character of label to be blank',
     :                 ' follow blank with comma'
 
               GO TO 21
            END IF
         END IF
      END IF
*
      IB = NWF - NIT + 1
      WRITE(ERR,'(/,A)') ' Default electron parameters ? (Y/N) '
      READ '(A)', ANS
      IF ( ANS.NE.'Y' .AND. ANS.NE.'y' .AND. NIT.NE.0) WRITE(ERR,*)
     :   ' S, IND, METH, ACC for electrons to be varied (free-format)'
      DO 20 I = NCLOSD+1,NWF
         IF (ANS.EQ.'Y' .OR. ANS.EQ.'y' .OR. I.LT.IB) THEN
                S(I) = SS
                        METH(I) = 3
                        IF (ID .EQ. 1) METH(I) = 1
                        ACC(I) = D0
                        IND(I) = 0
                        IF (IUF .NE. 0) IND(I) = -1
         ELSE
            WRITE(ERR,*) EL(I)
            READ(IN,*) S(I),IND(I),METH(I),ACC(I)
         END IF
         VARIED(I) = .TRUE.
         J = 2
         IF (EL(I)(1:1) .EQ. ' ') J = 3
         L(I) = LVAL(EL(I)(J:J))
         N(I) = ICHAR(EL(I)(J-1:J-1)) - ICHAR('1') + 1
 20   CONTINUE
*
*  *****  CHECK METHOD AND SET ORTHOGONALITY
*
        DO 35 NC = 1,NCFG
         STRING = CONFIG(NC)
         J = 2
         I = 0
 30      IF (STRING(J:J+2) .NE. '   ' ) THEN
*
*  *****     An electron has been found; find its index
*
            I = I+1
            ELI(I) = STRING(J:J+2)
            CALL EPTR(EL,ELI(I),IELI(I),*99)
            READ(STRING(J+4:J+5),'(I2)') IQ
            IF (NC.EQ.1 .AND. IQ.NE.D0 .AND. METH(IELI(I)).EQ.3
     :          .AND. (ANS.EQ.'Y' .OR. ANS.EQ.'y'))
     :        METH(IELI(I)) = 1
            J = J+8
            IF (J .LT. 40) GO TO 30
         END IF
*
*  *****  DEFINE ALL ORBITALS IN THE CONFIGURATION TO BE ORTHOGONAL
*
         DO 34 I1 = 2,I
            J1 = IELI(I1)
            DO 33 I2 = 1,I1-1
               J2 = IELI(I2)
               IF (L(J1) .EQ. L(J2) ) THEN
                    CALL EIJSET(J1,J2,1.D-5)
                    CALL EIJSET(J2,J1,1.D-5)
                 END IF
   33       CONTINUE
   34    CONTINUE
   35 CONTINUE
*
*  ***** SET THE FOLLOWING ORBITALS ORTHOGONAL
*
*        1) ORBITALS WITH DIFFERENT L'S
*        2) IN THE SAME ORTHOGONAL SET
*        3) SPECIFIED ORTHOGONALITY
*
*  *****
*
      DO 38 I = 2,NWF
         DO 39 J = 1,I-1
            IF (L(I) .EQ. L(J) .AND. SETORT(EL(I),EL(J)) ) THEN
                C = 1.D-5
                IF ((I.LE.NCLOSD .AND. J.LE.NCLOSD)
     :                          .OR. (ID .NE. 0))  C = 1.D-10
                CALL EIJSET(I,J,C)
                CALL EIJSET(J,I,C)
            END IF
   39    CONTINUE
   38  CONTINUE
*
*  *****  DETERMINE ADDITIONAL ORTHOGONALITY PAIRS
*
      I = 0
      IF ( IUC .NE. IN) THEN
   40    READ(IUC,1,END=50) EL1,EL2
         IF ( EL1 .NE. '*  ' .AND. EL2 .NE. '   ') THEN
            ELORT(I,1) = EL1
            ELORT(I,2) = EL2
            CALL EPTR(EL,EL1,I1,*99)
            CALL EPTR(EL,EL2,I2,*99)
            CALL EIJSET(I1,I2,1.D-5)
            CALL EIJSET(I2,I1,1.D-5)
            I = I +1
            IF (I .GT. (10)) STOP ' TOO MANY ORTHOGONALITIES: MAX=(10)'
            GO TO 40
         END IF
         NORT = I
      END IF
   50 CONTINUE
*
*  *****  ADDITIONAL PARAMETERS
*
      WRITE(ERR,'(/,A)') ' Default values (NO,REL,STRONG) ? (Y/N) '
      READ(IN,'(A)') ANS
      IF (ANS .EQ. 'Y' .OR. ANS .EQ. 'y') THEN
         NO = (NOD)
         REL = .FALSE.
         STRONG = .FALSE.
         IF (NCFG .GT. 1) STRONG = .TRUE.
       ELSE
         WRITE(ERR,*) ' Enter values in FORMAT(I,L,L) '
         READ(IN,*) NO, REL, STRONG
         IF (NO .GT. (NOD)) STOP
     :     ' TOO MANY POINTS FOR EACH FUNCTION: MAX=(NOD)'
      END IF
      ND = NO - 2
      WRITE(OUT,61) ATOM,TERM,Z,NO,NWF,NIT,NCFG,REL,STRONG
61    FORMAT(/1X,2A6,F6.0,I6,3I3,2L3)
      WRITE(PRI,62) ATOM,TERM,Z,(EL(I),4*L(I)+2,I=1,NCLOSD)
62    FORMAT(1H1///9X,33HHARTREE-FOCK WAVE FUNCTIONS FOR  ,2A6,4H Z =,
     1   F5.1//14X,'CORE = ',5(A3,'(',I2,')'))
      WRITE(PRI,'(//11X,A,37X,A//)') 'CONFIGURATION','WEIGHT'
      OMIT = .NOT. STRONG
*
* *****  WRITE 'CONFIGURATION' CARDS  AND NORMALIZE THE WEIGHTS
*
      W = D0
      DO 63 I=1,NCFG
   63 W = W + WT(I)**2
*
*  *****  IF INPUT WEIGHTS ARE ALL ZERO, SET WT(1)=1.
*
      IF (W .EQ. D0) THEN
         WT(1) = D1
         W = D1
      END IF
      W = DSQRT(W)
      DO 68 I = 1,NCFG
      WT(I) = WT(I)/W
      NOCC=NOCCSH(I)
      WRITE(PRI,70) I, CONFIG(I), WT(I),(COUPLE(I,J),J=1,NOCC)
70    FORMAT(/3X,I3,6X,A40,F19.8/12X,9(5X,A3))
      WRITE(PRI,73) (COUPLE(I,J),J=NOCC+1,2*NOCC-1)
73    FORMAT(23X,4(5X,A3))
68    CONTINUE
      WRITE(PRI,71)
71    FORMAT(//9X,10HINPUT DATA/9X,10H----- ----//13X,13HWAVE FUNCTION,
     1   11H  PROCEDURE/17X,22HNL  SIGMA METH ACC OPT///)
      DO 79 I = 1,NWF
      WRITE(PRI,78) I,EL(I),N(I),L(I),S(I),METH(I),ACC(I),IND(I)
78    FORMAT(I8, 2X,A3,2I3,F7.1,I4,F4.1,I4)
79    CONTINUE
*
*  *****  INITIALIZE ARRAYS, IF NECESSARY
*
      CALL WAVEFN
      DO 100 I=1,6
         INTPTR(I) = 0
100   CONTINUE
*
*        Initialize the VIJ data-structure
*
      IBEGIN = 1
      IEND = IEPTR(NWF)
      DO 45 I = IBEGIN,IEND
        VIJ(I) = D1
 45   CONTINUE
*
      IF (IUD .NE. IN) CALL INTGRL
*
*  *****  DEFINE SUM(I)
*
        IBEGIN = INTPTR(5)+1
        IEND = INTPTR(6)
        DO 80 I = IBEGIN,IEND
           IF (IEL(I,1).EQ.IEL(I,2)) SUM(IEL(I,1)) = -2*COEF(I)
 80     CONTINUE
      RETURN
99    STOP
      END
*
*     ------------------------------------------------------------------
*              D E
*     ------------------------------------------------------------------
*
*       This routine controls the solution of the differenttial equation
*   for the radial function P  .  One of three methods is selected -
*                            I1
*   M1, M2, or M3 -  for solving the equations,  the  initial  choice
*   being determined by an input paramter, METH(I), except when no
*   exchange is present, in which case M2 is selected. 
*
*        Value of METH(I)     Method
*        ---------------      ------
*        < or =1            M1 with search for an acceptable solution
*             =2            M2 with search for an acceptable solution
*             =3            M3 without any checking
*
*   If M1 fails to find an acceptable solution, the radial  functions
*   are  orthogonalized,  off-diagonal  energy parameters recomputed,
*   and the method tried again.   Should it continue to fail, METH(I)
*   is set to 2.
*
*
      SUBROUTINE DE(I1)
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
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
        LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
        COMMON P2(NOD),HQ(NOD),XX(NOD),AC(30,NWD),BC(NWD),JV(NWD),
     :     AZZ,PP,FN,EM,FM,EU,FU,DELTAE,M,NODE,MK,KK,NJ
*
        LOGICAL DIAG
        CHARACTER*2 ASTER(3)
        DATA ASTER/'  ','* ','**'/
*
      I = I1
      ED2 = E(I,I)
      KK= MAX0(1,METH(I))
      IF (NWF .EQ. 1) KK = 2
      NODE = N(I) - L(I) - 1
*
*  *****  CALL METHD1 TO SOLVE THE DIFFERENTIAL EQUATION
*
      CALL METHD1(I)
      IF ( FAIL ) GO TO 25
*
12    PN = DSQRT(QUAD(I,M,PDE,PDE))
      DO 9 J = 1,M
9     PDE(J) = PDE(J)/PN
      AZZ = AZZ/PN
*
*  *****  CHECK IF METHOD 2 SHOULD BE USED
*
      IF ( KK .EQ. 1 ) THEN
         IF (DABS(D1 -ED2/E(I,I)) .LT. 0.005D0  .AND.
     :       DMAX1(DABS(D1 - PN), DABS(D1/PN - D1)) .GT. 0.2D0 ) THEN
             METH(I) = 2
             KK = 2
             GO TO 25
         END IF
      ELSE IF (DABS(D1 -ED2/E(I,I)) .LT. 0.0001D0 .AND. IC .GT. 1) THEN
         IC = IC-1
      END IF
*
*  *****  SET THE ACCELERATING PARAMETER
*
*
13    IF (IPR .NE. I ) GO TO 33
      ED2 = ED2 - E(I,I)
      IF (ED1*ED2 .GT. D0) ACC(I) = .75*ACC(I)
      IF (ED1*ED2 .LT. D0) ACC(I) = (D1 + D3*ACC(I))/D4
33    C = ACC(I)
      CD = D1 - C
*
*   *****  IMPROVE THE ESTIMATES
*
      MAX(I) = M
      DP     = D0
      DO 21 J = 1,M
      DIFF = P(J,I)-PDE(J)
      DP     = DMAX1(DP    ,DABS(DIFF)*R2(J))
21     P(J,I) = PDE(J) + C*DIFF
      IF (M .EQ. NO) GO TO 26
      M = M + 1
      DO 24 J = M,NO
24     P(J,I) = D0
      AZ(I) = CD*AZZ + C*AZ(I)
      AZZ = AZ(I)
*
*  *****  CHECK THE ORTHOGONALIZATION
*
26    CONTINUE
      IBEGIN = 1
      IF (I .GT. 1) IBEGIN = IEPTR(I-1) + 1
      IP = IBEGIN
      IJ = 0
50    JI = IJE(IP)
      IF (JI .NE. I ) THEN
      IF (JI .GE. IB .AND. DPM(JI) .GE. DPM(I)) THEN
*               
*               The JI orbital should be orthogonalized
*
         C = QUADR(I,JI,0)
         MM = MAX0(MAX(JI),MAX(I))
         DO 51 J = 1,MM
            P(J,JI) = P(J,JI) - C*P(J,I)
51       CONTINUE
         C2 = SQRT(QUADR(JI,JI,0))
         DO 52 J = 1,MM
            P(J,JI) = P(J,JI)/C2
52       CONTINUE
         VARIED(JI) = .TRUE.
         MAX(JI) = MM
         AZ(JI) = (AZ(JI) - C*AZ(I))/C2
         WRITE(OUT,63) EL(I),EL(JI),C
      ELSE
*
*              The I'th orbital must be orthogonalized
*
         IJ = IJ + 1 
         IF (IJ .GT. 20) STOP ' TOO MANY ORTHOGONALITY CONDITIONS' 
         JV(IJ) = JI 
      END IF
      END IF
      IP = IP + 1
      IF (IP .LE. IEPTR(I)) GO TO 50
      IF (IJ .NE. 0 ) THEN
         DIAG = .TRUE.
         DO 61 J = 1,IJ
            BC(J) = QUADR(I,JV(J),0)
            AC(J,J) = D1
            DO 62 JJ = J+1, IJ
               IF (E(JV(J),JV(JJ)) .NE. D0 ) THEN
                  AC(J,JJ) = D0
                  AC(JJ,J) = D0
                ELSE
                  AC(J,JJ) = QUADR(JV(J),JV(JJ),0)
                  AC(JJ,J) = AC(J,JJ)
                  DIAG = .FALSE.
               END IF
62          CONTINUE
61       CONTINUE
         IF ( .NOT. DIAG .AND. IJ .GT. 1) CALL LINEQN(30,IJ,AC,BC)
         M = MAX(I)
         DO 65 JJ = 1,IJ
            C = BC(JJ)
            WRITE(OUT,63) EL(JV(JJ)),EL(I),C
63          FORMAT(6X,'<',A3,'|',A3,'>=',1PD8.1)
            M = MAX0(M,MAX(JV(JJ)))
            DO 64 J = 1,M
               P(J,I) = P(J,I) - C*P(J,JV(JJ))
64          CONTINUE
            AZZ = AZZ - C*AZ(JV(JJ))
65       CONTINUE
         PNN = DSQRT(QUADR(I,I,0))
         DO 66 J = 1,M
            P(J,I) = P(J,I)/PNN
66       CONTINUE
         AZZ = AZZ/PNN
      END IF
      M = NO
67    IF (DABS(P(M,I)) .LT. 1.D-15) THEN
         P(M,I) = D0
         M = M-1
         GO TO 67
      END IF
      MAX(I) = M
      IF (AZZ .GT. D0) AZ(I) = DMAX1(AZZ,D5*AZ(I))
      WRITE(OUT,17) EL(I),E(I,I),AZ(I),PN,ASTER(KK),DP
17    FORMAT(20X,A3,2F15.7,F12.7, A2,1PD10.2)
      DPM(I) = DP
      IF (IPR .EQ. I1) ED1 = ED2
      IF (IPR .NE. I1) ED1 = ED2 - E(I1,I1)
      IPR = I1
      VARIED(I) = .TRUE.
      RETURN
*
*  *****  IF METHD1 FAILED TO FIND AN ACCEPTABLE SOLUTION, ORTHOGONALIZE
*  *****  THE ESTIMATES AND TRY AGAIN
*
25    IF (I .EQ. IB) GO TO 27
      CALL ORTHOG
      CALL UPDATE
      CALL GRANGE
27    CALL METHD1(I)
      IF ( FAIL ) GO TO 23
      GO TO 12
*
*  *****  ERROR RETURN FROM SECOND TRY.  IF M1 WAS USED,SWITCH TO
*         M2 AND TRY ONCE MORE.
*
23    IF ( KK .EQ. 2) RETURN
      KK = 2
      GO TO 27
      END
*-----------------------------------------------------------------------
*               D I A G
*-----------------------------------------------------------------------
*     This routines sets up the interaction matrix and, given an 
*  approximate eigenvector, finds an eigenvalue and eigenvector.
*
      SUBROUTINE DIAG(ECONV,ACFG,CFGTOL,LAST)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER(NWD=30,NOD=220,NCD=100)
*
        INTEGER IN,OUT,ERR,PRI,OUC,OUD,OUF,OUH
        COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,OUC,OUD,OUF,OUH
*
        CHARACTER CONFIG*40,EL*3,ATOM*6,TERM*6,COUPLE*3
        COMMON /LABEL/CONFIG(NCD),EL(NWD),ATOM,TERM,COUPLE(NCD,9)
*
        COMMON /MATRIX/ETOTAL,W(NCD,NCD)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
        INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
        PARAMETER (IDIM=550,NCDIM=3000)
        COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :         ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
        LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
      COMMON WP(NCD)
*
        LOGICAL ECONV,LAST
*
        IF (LAST.AND.OUH.NE.0) WRITE(OUH,'(2X,A6,I6)') 'NCFG =',NCFG
        DO 1 I = 1,NCFG
           DO 2 J = 1,NCFG
              W(I,J) = D0
  2        CONTINUE
           WP(I) = WT(I)
           W(I,I) = EC
  1     CONTINUE
*
        IBEGIN = 1
        IEND = INTPTR(6)
        J = 0
        DO 10 I = IBEGIN,IEND
 11        IF (CPTR(I) .GT. J) THEN
              J = J + 1
              C = COEFF(J)*VALUE(I)
              IF (OPTR(J) .NE. 0) C = C*VALUE(OPTR(J))
              W(IH(J),JH(J)) = W(IH(J),JH(J)) + C
              GO TO 11
           END IF
 10     CONTINUE
*
*  ***** SYMMETRIZE THE MATRIX
*
        DO 12 I = 1,NCFG-1
           DO 13 J = I+1,NCFG
              W(I,J) = W(J,I)
 13        CONTINUE
 12     CONTINUE
      IF (NCFG .EQ. 1) GO TO 37
      IF (ID .GT. 0) GO TO 38
*
*  ***** Write the matrix to disk
*
      DO 14 I = 1,NCFG
         WRITE(35) (W(J,I),J=I,NCFG)
 14   CONTINUE
      ENDFILE(35)
      REWIND(35)
*
*  *****  COMPUTE an EIGENVALUE AND EIGENVECTOR.
*  *****  THIS METHOD MAY CAUSE DIFFICULTIES WHEN NEAR
*  *****  DEGENERACY EFFECTS ARE PRESENT SINCE IT MAY CONVERGE TO
*  *****  THE WRONG EIGENVALUE.  THE CODE UP TO, BUT NOT INCLUDING
*  *****  STATEMENT NUMBER 31, MAY BE REPLACED BY A CALL TO A MORE
*  *****  REFINED LIBRARY SUBROUTINE.
*
      ETPREV = D0
      DO 30 II=1,5
      IF (II .NE. 1) THEN
         DO 5 I = 1,NCFG
           READ(35) (W(J,I),J=I,NCFG)
           DO 8 J = I+1,NCFG
             W(I,J) = W(J,I)
   8       CONTINUE
   5     CONTINUE
         REWIND(35)
      END IF
      ETL = D0
*
*  *****  DETERMINE ESTIMATE OF THE EIGENVALUE
*
      DO 16 I=1,NCFG
      CONT = D0
      DO 17 J=1,NCFG
17    CONT = CONT + WT(J)*W(I,J)
16    ETL = ETL + WT(I)*CONT
*
*  *****  SOLVE SYSTEM OF EQUATIONS FOR EIGENVECTOR
*
      DO 18 I=1,NCFG
18    W(I,I) = W(I,I) - ETL
      WT(1) = D1
      IF (NCFG .NE. 2) GO TO 54
      WT(2) = -W(2,1)/W(2,2)
      GO TO 20
*
54    CALL SYMMEQ((NCD),NCFG,W,WT)
*
20    DO 26 I=2,NCFG
26    WT(1) = WT(1) + WT(I)**2
      WT(1) = D1/DSQRT(WT(1))
      DO 27 I=2,NCFG
27    WT(I) = WT(I)*WT(1)
*
*  *****  ITERATE, IF NECESSARY, OTHERWISE OUTPUT RESULTS
*
      IF (DABS((ETPREV-ETL)/ETL) .LT. 1.D-7) GO TO 31
30    ETPREV = ETL
      WRITE(OUT,40)
40    FORMAT(///10X,47HMATRIX DIAGONALIZATION PROCEDURE NOT CONVERGING )
      DELTAE = D0
      GO TO 33
31    DELTAE = ETL -ETOTAL
33    ETOTAL = ETL
      CC = D0
      DO 28 I = 1,NCFG
      WT(I) = WT(I) + ACFG*(WP(I) - WT(I))
28    CC = CC + WT(I)**2
      CC = D1/DSQRT(CC)
      DO 35 I=1,NCFG
35    WT(I) = WT(I)*CC
      IF ( LAST ) THEN
         WRITE(PRI,34)
         WRITE(PRI,36) (I,WT(I),I=1,NCFG)
 36      FORMAT(6(I5,F8.4))
         WRITE(PRI,32) ETOTAL
*
* *****  Read matrix back from disk and print
*
         DO 60 I = 1,NCFG
           READ(35) (W(J,I),J=I,NCFG)
           DO 68 J = I+1,NCFG
             W(I,J) = W(J,I)
  68       CONTINUE
           IF (OUH.NE.0) WRITE(OUH,'(5F14.7)') (W(I,J),J=1,I)
           WRITE(PRI,62) I, (W(I,J),J=1,I)
 62        FORMAT(I5,2X,5F15.7/(7X,5F15.7))
  60     CONTINUE
         REWIND(35)
      END IF
      WRITE(OUT,32) ETOTAL  
 32   FORMAT(//10X,15HTOTAL ENERGY = ,F18.9 )
 39   WRITE(OUT,34)
 34   FORMAT(/6X,'WEIGHTS')
      WRITE(OUT,'(6(I3,F8.4))') (I,WT(I), I=1,MIN0(12,NCFG))
*
*  *****  REDEFINE SUM(I)
*
        IBEGIN = INTPTR(5)+1
        IEND = INTPTR(6)
        DO 50 I = IBEGIN,IEND
           IF (IEL(I,1).EQ.IEL(I,2)) SUM(IEL(I,1)) = -2*COEF(I)
 50     CONTINUE
70    IF (.NOT. LAST .OR. OUC .EQ. 0 ) GO TO 49
*
*  *****  PUNCH CONFIGURATIONS AND WEIGHTS ON UNIT OUC
*
*
      WRITE(OUC,46) ATOM,TERM,ETOTAL
46    FORMAT(3X,2A6,F14.7)
      WRITE(OUC,'(18(1X,A3))') (EL(J),J=1,NCLOSD)
      DO 47 J = 1,NCFG
47    WRITE(OUC,48) CONFIG(J),WT(J),(COUPLE(J,JJ),JJ=1,9)
48    FORMAT(A40,F10.7/9(5X,A3))
      WRITE (OUC,'(A)') '*'
49    ECONV = .FALSE.
      IF (DABS(DELTAE/ETOTAL) .LE. CFGTOL) ECONV = .TRUE.
      RETURN
37    ETOTAL = W(1,1)
      DELTAE = D0
      WRITE(OUT,32) ETOTAL
      GO TO 70
38    DELTAE = D0
      ETOTAL = W(1,1)
      GO TO 39
      END
*
*     ------------------------------------------------------------------
*              D I F F
*     ------------------------------------------------------------------
*
*
*       Stores LP  in the array YK.  The difference approximation of
*                i
*   Eq. (6-14) is used.
*
*
      SUBROUTINE DIFF(I)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220,NCD=100)
*
        INTEGER IN,OUT,ERR,PRI,OUC,OUD,OUF,OUH
        COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,OUC,OUD,OUF,OUH
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
*  *****  FORM DD + 2Z/R -L(L+1)/RR|P(I)>
*
      MM = MAX(I) - 3
      FL = L(I)
      TWOZ = Z + Z
      C = (FL+D5)**2
      HH = 180.D0*H*H
      DO 11 K =  4,MM
11    YK(K) = (D2*(P(K+3,I)+P(K-3,I)) - 27.D0*(P(K+2,I)+P(K-2,I)) +
     1   270.D0*(P(K+1,I)+P(K-1,I)) - 490.D0*P(K,I))/HH +
     2   P(K,I)*(TWOZ*R(K) - C)
*
*  *****  BECAUSE OF THE POSSIBILITY OF EXTENSIVE CANCELLATION NEAR THE
*  *****  ORIGIN, SEARCH FOR THE POINT WHERE THE ASYMPTOTIC BEHAVIOUR
*  *****  BEGINS AND SMOOTH THE ORIGIN.
*
      LEXP = L(I) + 2
      Y1 = YK(4)/R2(4)/R(4)**LEXP
      Y2 = YK(5)/R2(5)/R(5)**LEXP
      DO 1 K = 4,100
      KP = K+2
      Y3 = YK(KP)/R2(KP)/R(KP)**LEXP
      IF (Y2 .EQ. D0) GO TO 1
      IF (DABS(Y1/Y2 - D1) .LT..1D0 .OR. DABS(Y3/Y2 - D1) .LT..1D0)
     1       GO TO 2
      Y1 = Y2
      Y2 = Y3
1     CONTINUE
      WRITE(OUT,3)  I
3     FORMAT(6X, 'ASYMPTOTIC REGION NOT FOUND FOR FUNCTION NUMBER',I3)
      STOP
*
*  *****  ASYMPTOTIC REGION HAS BEEN FOUND
*
2     KP = K
      KM = KP - 1
      DO 4 K = 1,KM
4     YK(K) = Y1*R2(K)*R(K)**LEXP
      MM = MM + 1
      YK(MM) = (-(P(MM+2,I)+P(MM-2,I)) + D16*(P(MM+1,I)+P(MM-1,I))
     1      -D30*P(MM,I))/(D12*H*H) + P(MM,I)*(TWOZ*R(MM) - C)
      MM = MM + 1
      YK(MM) = (P(MM+1,I) + P(MM-1,I) - D2*P(MM,I))/(H*H) +
     1   P(MM,I)*(TWOZ*R(MM) - C)
      MM = MM+1
      DO 5 K =MM,NO
5     YK(K) = D0
      RETURN
      END
*
*-----------------------------------------------------------------------
*               E
*-----------------------------------------------------------------------
*      Returns the value of the off-diagonal energy parameter
*  for the (i,j) pair from the data structure.
*
        DOUBLE PRECISION FUNCTION E(I,J)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220,NCD=100)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
        IBEGIN = 1
        IF (I .GT. 1) IBEGIN = IEPTR(I-1) + 1
        IEND = IEPTR(I)
        E = 0.D0
        DO 10 II = IBEGIN,IEND
           IF (IJE(II) .EQ. J) THEN
              E = EIJ(II)
              RETURN
           END IF
 10     CONTINUE
        END
*
*-----------------------------------------------------------------------
*               E I J S E T
*-----------------------------------------------------------------------
*
*      Stores the value of the off-diagonal energy parameter for the
*   pair (i,j) in the data structure
*
        SUBROUTINE EIJSET(I,J,E)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220,NCD=100)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
        IBEGIN = 1
        IF (I .GT. 1) IBEGIN = IEPTR(I-1)+1
        IEND = IEPTR(I)
        DO 10 II = IBEGIN,IEND
           IF (IJE(II) .EQ. J) THEN
              EIJ(II) = E
              RETURN
           END IF
 10     CONTINUE
*
* ***** J-value not found - enter into list
*
        IF (IJE(98) .NE. 0)
     :   STOP ' Too many off-diagonal energy parameters'
*
*  ***** Find point at which the insertion should be made
*
        IEND = IEPTR(I)
        IF (IEND .NE. 0) THEN
           IP = 1
           IF (I .GT. 1) IP = IEPTR(I-1)+1
 30        IF (IJE(IP) .LT. J .AND. IP .LE. IEND) THEN
              IP = IP + 1
              GO TO 30
           END IF
        ELSE
           IP = 1
        END IF
*
* *****  IP is the location in which EIJ should be stored
*        Move other data
*
        
        DO 40 JJ = (98)-1,IP,-1
           IJE(JJ+1) = IJE(JJ)
           EIJ(JJ+1) = EIJ(JJ)
 40     CONTINUE
*
* ***** Space has been made - insert data
*
        IJE(IP) = J
        EIJ(IP) = E
*
* ***** Update pointers
*
        DO 50 II = I,NWF
           IEPTR(II) = IEPTR(II) +1
 50     CONTINUE
        END
*
*-----------------------------------------------------------------------
*               E P S I L O N
*-----------------------------------------------------------------------
*
*       Analysis the effect of a rotation of an (i,j) pair of orbitals 
*   on the total energy.  If the energy is invariant, the off-diagonal 
*   energy parameter data structure is set to reflect this result.
*   Otherwise, rotation for a stationary energy is determined.
*
        SUBROUTINE EPSILON(IIN,JIN)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220,NCD=100)
*
        INTEGER IN,OUT,ERR,PRI,OUC,OUD,OUF,OUH
        COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,OUC,OUD,OUF,OUH
*
        CHARACTER CONFIG*40,EL*3,ATOM*6,TERM*6,COUPLE*3
        COMMON /LABEL/CONFIG(NCD),EL(NWD),ATOM,TERM,COUPLE(NCD,9)
*
        COMMON /MATRIX/ETOTAL,W(NCD,NCD)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
        INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
        PARAMETER (IDIM=550,NCDIM=3000)
        COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :         ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
        LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
        COMMON U(NCD),DWT(NCD)
*
        LOGICAL ICLOSD,FOUND

        IF (IIN .LT. JIN) THEN
           I = IIN
           J = JIN
        ELSE
           I = JIN
           J = IIN
        END IF
        SUMI = SUM(I)
        SUMJ = SUM(J)
        LI = L(I)
        LJ = L(J)
        FULL = 4*LI+2
        IF (ABS(SUMI-FULL).LT.1.D-8 .AND. ABS(SUMJ-FULL).LT.1.D-8) THEN
           CALL EIJSET(I,J,1.D-10)
           CALL EIJSET(J,I,1.D-10)
           RETURN
        ELSE IF ((ABS(SUMI-FULL).LT.1.D-8.AND.ABS(SUMJ-FULL).LT.1.D-2)
     :   .OR.(ABS(SUMJ-FULL).LT.1.D-8.AND.ABS(SUMI-FULL).LT.1.D-2)) THEN
           CALL VIJSET(J,I,1.D-5)
           RETURN
        END IF
*
        VPR = V(I,J)
        DVPR = V(J,I)
        IF (VPR .NE. D0 .AND. DABS(VPR) .LT. 1.D-4) RETURN

        DO 2 II = 1,NCFG
           U(II) = D0
2       CONTINUE

        ALL = .TRUE.
        IF (VPR .NE. D0) THEN
           IF (ABS(DVPR/VPR) .LT. 0.002) ALL = .FALSE.
        END IF
        G = D0
        DG = D0
        DIFF = SUMI - SUMJ
        IF (I .LE. NCLOSD) THEN
            ICLOSD = .TRUE.
            HLIJ = HL(EL,I,J,REL)
            IF (ALL) 
     :          DG=DIFF*(HL(EL,I,I,REL)-HL(EL,J,J,REL))
            VI =  -SUMI*HLIJ
            VJ =  HLIJ
            TI = D0
            TJ = D0
            TT = D0
            VV = D0
            DO 10 K = 0,2*LI,2
                IF (K .EQ. 0) THEN
                   CIJ = (D1 - CB(LI,LJ,0))
                   CII = D1
                ELSE
                   CIJ = -CB(LI,LJ,K)
                   CII = -CA(LI,K)
                END IF
                RKIJJJ = RK(I,J,J,J,K,REL)
                RKJIII = RK(J,I,I,I,K,REL)
                TI = TI + CIJ*RKIJJJ
                TT = TT + (SUMI-D1)*CII*RKJIII
                TJ = TJ + CIJ*RKJIII
                IF (ALL) THEN
                   FKIJ = FK(I,J,K,REL)
                   GKIJ = GK(I,J,K,REL)
                   FKII = FK(I,I,K,REL)
                   FKJJ = FK(J,J,K,REL)
                   VV = VV + SUMJ*CIJ*(FKJJ + FKII - 2*FKIJ - 4*GKIJ)
     :                   +(SUMI-D1)*CII*(FKIJ + 2*GKIJ - FKII)
                END IF
 10        CONTINUE
           DG = DG + D2*SUMI*VV
           VI = VI + D2*SUMI*(TT + SUMJ*TI)
           VJ = VJ - D2*SUMI*TJ
           VV = D0
           DO 20 IA = 1,NCLOSD
                IF (IA .EQ. I) GO TO 20
                T = RK(I,IA,J,IA,0,REL)
                IF (ALL) 
     :              VV = FK(J,IA,0,REL)-FK(I,IA,0,REL)
                DO 22 K = IABS(LI-L(IA)),LI+L(IA),2
                   CIJ = CB(LI,L(IA),K)
                   T = T - CIJ*RK(I,IA,IA,J,K,REL)
                   IF (ALL)
     :                   VV = VV - CIJ*(GK(J,IA,K,REL)-GK(I,IA,K,REL))
 22             CONTINUE
                DG = DG + D2*DIFF*SUM(IA)*VV
                VI = VI + D2*SUMI*SUM(IA)*T
                VJ = VJ - D2*SUM(IA)*T
 20        CONTINUE
           G = G + VI
           VJ = -D2*VJ
           DO 24 II = 1,NCFG
                U(II) = U(II) + WT(II)*VI
 24        CONTINUE
        ELSE
*
*     ... Both orbitals are outer orbitals
*
           ICLOSD = .FALSE.
           HLCIJ = HLC(EL,I,J,REL)
           IF (ALL) 
     :        DG=DIFF*(HLC(EL,I,I,REL)-HLC(EL,J,J,REL))
           VI =   D2*HLCIJ
           VJ =  -D2*HLCIJ
        END IF
*
*
*       ... Add contributions from integrals between outer electrons
*
        IEND = 0
        DO 40 INT = 1,6
           IBEGIN = IEND + 1
           IEND = INTPTR(INT)
           DO 30 II = IBEGIN,IEND
                FOUND = .FALSE.
                I1 = IEL(II,1)
                I2 = IEL(II,2)
                IR = (I1-I)*(I2-I)
                JR = (I1-J)*(I2-J)
                IF (IR .EQ. 0 .OR. JR .EQ. 0 .OR.
     :             (INT.EQ.6 .AND. ICLOSD)) THEN
*               ... We have found a contribution ...
                   FOUND = .TRUE.
                   IF (INT.EQ.5) THEN
                        CALL RPERT(I,J,II,INT,DC,DV)
                   ELSE IF (INT.EQ.4) THEN
                        CALL O2PERT(I,J,II,INT,DC,DV)
                   ELSE IF (INT.EQ.6 .AND. (IR.EQ.0 .OR. JR.EQ.0) 
     :                  .AND. I1.EQ.I2) THEN
                      IF (I1.EQ.I) THEN
                          DC = VI
                          DV = D0
                      ELSE IF (I1.EQ.J) THEN
                          DC = VJ
                          DV = D0
                      END IF
                   ELSE
                        CALL PERT(I,J,II,INT,DC,DV)
                   END IF
                ELSE IF (INT.EQ.5 .OR. INT.EQ.4) THEN
*                  ... Continue to search for a contribution
                   I1 = IEL(II,3)
                   I2 = IEL(II,4)
                   IR = (I1-I)*(I2-I)
                   JR = (I1-J)*(I2-J)
                   IF (IR .EQ. 0 .OR. JR .EQ. 0) THEN
*                  ... Another contribution has been found
                      FOUND = .TRUE.
                      IF (INT.EQ.5) CALL RPERT(I,J,II,INT,DC,DV)
                      IF (INT.EQ.4) CALL O2PERT(I,J,II,INT,DC,DV)
                   END IF
                END IF
                IF (FOUND) THEN
                   IF (INT .NE. 3 .AND. INT .NE. 4) THEN
                      CALL CONTC(II,DC,U,C)
                   ELSE
                      CALL CONTOV(II,DC,U,C)
                   END IF
                   G = G + C*DC
                   DG = DG + C*DV
                END IF
 30        CONTINUE
 40     CONTINUE
*
        IF (.NOT. ALL) THEN
           DG = VPR
           GO TO 38
        END IF
        DO 32 II = 1,NCFG
           U(II) = U(II) - G*WT(II)
           W(II,1) = U(II)
 32     CONTINUE
        W(1,1) = D0
*
*  *****  SOLVE SYSTEM OF EQUATIONS FOR PERTURBATIONS TO THE EIGENVECTOR
*
        
      IF (NCFG .EQ. 1) GO TO 38
      IF (NCFG .EQ. 2) THEN
         W(2,1) = -W(2,1)/W(2,2)
      ELSE
           CALL SYMMSL((NCD),NCFG,W,W(1,1))
      END IF
      C = D0
      DO 36 JJ = 1,NCFG
         C = C + WT(JJ)*W(JJ,1)
 36   CONTINUE
      DO 37 JJ = 1,NCFG
         W(JJ,1) = W(JJ,1) - C*WT(JJ)
         DG = DG + D2*W(JJ,1)*U(JJ)
 37   CONTINUE
 38   IF (DABS(G)+DABS(DG) .GT. 4.D-5*((SUMI*SUMJ)**0.25)
     :         .OR. DABS(E(I,J)) .GT. 2.D-5) THEN
         EPS = -G/DG
         IF (DABS(DG) .LT. 1.D-4) THEN
            EPS = 2.D-5
         ELSE
            EPS = DSIGN(DMAX1(DMIN1(DABS(EPS),0.025D0),1.1D-10),EPS)
         END IF
         IF (ALL .AND. E(I,J).NE.1.D-5) THEN
            DO 42 JJ = 1,NCFG
                DWT(JJ) = DWT(JJ) + EPS*W(JJ,1)
  42        CONTINUE
         END IF
         CALL EIJSET(J,I,EPS)
         WRITE(OUT,100) EL(I),EL(J),G,EL(I),EL(J),DG,EPS
100      FORMAT(10X,'C(',2A3,') =',F12.5,3X,'V(',2A3,') =',
     :                  F12.5,3X,'EPS =',F9.6)
       ELSE
*
*  *****  THE ENERGY IS STATIONARY WITH RESPECT TO ROTATIONS
*
        CALL EIJSET(I,J,1.D-10)
        CALL EIJSET(J,I,1.D-10)
      END IF
        IF (ALL) THEN
           CALL VIJSET(I,J,DG)
           DDG = DG - VPR
           CALL VIJSET(J,I,DDG)
        END IF
      END
*
*     ------------------------------------------------------------------
*              G R A N G E
*     ------------------------------------------------------------------
*
*       Controls the calculation of off-diagonal energy parameters.
*   It searches for all pairs (i,j) which are constrained through  an
*   orthogonality requirement.   When one of the pair , say P
*                                                            i
*   must be orthogonal not only to P  but also to P  where n = n ,
*                                   j              k        j   k
*   a system of equations must be solved, the exact form depending on
*   whether or not any of the functions are part of the frozen  core.
*
      SUBROUTINE GRANGE
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
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
        INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
        PARAMETER (IDIM=550,NCDIM=3000)
        COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :         ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
        LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
        COMMON U(NCD),DWT(NCD),AC(30,NWD),BC(NWD),JV(NWD),IV(NWD)
        LOGICAL DIAG, SETEQL, FIRST
*
*       CLEAR THE ARRAY FOR CHANGING THE WT ARRAY
 
*
        DO 3 I = 1,NCFG
           DWT(I) = D0
 3      CONTINUE
*
*  *****  ROTATE PAIRS CONNECTED BY ORTHOGONALITY BUT NOT WHEN ONE OF
*         THE ORBITALS IS SIMULTANEOUSLY ORTHOGONAL TO A NON-ORTHOGONAL
*         PAIR
*
        DO 1 I = IB,NWF-1
           DO 2 J = I+1,NWF
              IF (DABS(E(I,J)) .GT. 1.D-10 .AND. SETEQL(I,J))
     :            CALL EPSILON(I,J)
2          CONTINUE
1       CONTINUE
        CALL ROTATE
*
*       Adjust the WT array, renormalize, and recompute SUM(i)
*
        C = D0
        DO 4 I = 1,NCFG
           WT(I) = WT(I) + DWT(I)
           C = C + WT(I)**2
 4      CONTINUE
        C = SQRT(C)
        DO 5 I = 1,NCFG
           WT(I) = WT(I)/C
 5      CONTINUE
*
        IBEGIN = INTPTR(5)+1
        IEND = INTPTR(6)
        DO 6 I = IBEGIN,IEND
           IF (IEL(I,1).EQ.IEL(I,2)) SUM(IEL(I,1)) = -2*COEF(I)
 6      CONTINUE
*
*  *****  FOR EACH l COMPUTE OFF-DIAGONAL ENERGY PARAMETERS
*
        DO 10 IL = 0,5
           IJ = 0
           DO 11 I = IB,NWF
           IF ( L(I) .NE. IL ) GO TO 11
           DO 12 J = 1,I-1
              IF (DABS(E(I,J)) .GT. 1.D-10 ) THEN
                 IJ = IJ + 1
                 IF ( IJ.GT.(NWD)) STOP '  TOO MANY LAGRANGE MULTIPLIERS'
                 IV(IJ) = I
                 JV(IJ) = J
              END IF
12         CONTINUE
11         CONTINUE
*
*  ***** IJ IS THE NUMBER OF LAGRANGE MULTIPLIERS FOR l = IL
*
           IF (IJ .EQ. 0) GO TO 10
           DO 13 II = 1,IJ
              BC(II) = D0
              DO 14 III = 1,IJ
                 AC(II,III) = D0
14            CONTINUE
13         CONTINUE
           DO 16 I = IB,NWF
              IF ( L(I) .NE. IL ) GO TO 16
                 FIRST = .TRUE.
                 DO 18 II = 1,IJ
                    J = 0
                    IF ( IV(II) .EQ. I) THEN
                       J = JV(II)
                     ELSE IF ( JV(II) .EQ. I) THEN
                       J = IV(II)
                    END IF
                    IF ( J .NE. 0) THEN
                       IF (FIRST) THEN
                          CALL XCH(I,2)
                          CALL POTL(I)
                          DO 20 JJ = 1,NO
                             YK(JJ) = YR(JJ)
20                        CONTINUE
                          FIRST = .FALSE.
                       END IF
                       DO 22 JJ = 1,NO
                          YR(JJ) = P(JJ,J)
22                     CONTINUE
                       BC(II) = BC(II) +
     :                    HL(EL,I,J,REL)-D2*QUADS(I,J,1)-QUAD(J,NO,YR,X)
                    END IF
18              CONTINUE
16         CONTINUE
           DO 24 II = 1,IJ
              DO 26 III = 1,II
                 IF ( II .EQ. III) THEN
                    AC(II,II) = D1/SUM(IV(II))
                    IF (JV(II) .GE. IB) THEN
                       AC(II,II) = AC(II,II) + D1/SUM(JV(II))
                    END IF
                  ELSE IF (IV(II) .EQ. IV(III) .AND.
     :                    E(JV(II),JV(III)) .EQ. D0 ) THEN
                     AC(II,III) = QUADR(JV(II),JV(III),0)/SUM(IV(II))
                     AC(III,II) = AC(II,III)
                  ELSE IF (JV(II) .EQ. JV(III) .AND. JV(II) .GE. IB
     :              .AND. E(IV(II),IV(III)) .EQ. D0) THEN
                     AC(II,III) = QUADR(IV(II),IV(III),0)/SUM(JV(II))
                     AC(III,II) = AC(II,III)
                  END IF
26            CONTINUE
24         CONTINUE
           CALL LINEQN(30,IJ,AC,BC)
           DO 28 II = 1,IJ
              CALL EIJSET(IV(II),JV(II),BC(II)/SUM(IV(II)))
              IF ( JV(II) .GE. IB )
     :            CALL EIJSET(JV(II),IV(II),BC(II)/SUM(JV(II)))
28         CONTINUE
10    CONTINUE
*
*  *****  PRINT THE OFF-DIAGONAL ENERGY PARAMETERS
*
        DO 30 I = IB,NWF
           DO 32 J = 1,I-1
              IF (DABS(E(I,J)) .GT. 1.D-10) THEN
                 WRITE(OUT,35) EL(I),EL(J),E(I,J),EL(J),EL(I),E(J,I)
35               FORMAT(7X,2(3X,'E(',2A3,') =',F12.5))
              END IF
32         CONTINUE
30      CONTINUE
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
     :         ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
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
*
*     ------------------------------------------------------------------
*              M E T H O D
*     ------------------------------------------------------------------
*
*       Uses M1, M2, or M3 to solve the radial equation. If the input
*   data indicated METH(I) = 3, then this  solution  is  returned  to
*   DE.  Otherwise,  the routine searches for an acceptable  solution
*   which  is  both  positive  near  the  origin and has the required
*   number  of nodes.
*
*
      SUBROUTINE METHD1(I)
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
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
        LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
      COMMON P2(NOD),HQ(NOD),XX(NOD),AC(30,NWD),BC(NWD),JV(NWD),
     1     AZZ,PP,FN,EM,FM,EU,FU,DELTAE,M,NODE,MK,KK,NJ
*
      LOGICAL V2, FIRST
      DIMENSION P1(NOD)
      EQUIVALENCE (PDE(1),P1(1))
*
*  *****  'FIRST' MUST BE 'TRUE' THE FIRST TIME SOLVE IS CALLED FOR
*  *****  POTENTIAL AND EXCHANGE TO BE COMPUTED
*  *****  'EU' IS THE UPPER BOUND OF THE ENERGY PARAMETER
*  *****  'EM' IS THE MINIMUM VALUE OF THE ENERGY PARAMETER
*
      FIRST = .TRUE.
      FAIL = .FALSE.
      EM = D0
      EU = ((Z - DMIN1(D5*S(I),D2*S(I)))/N(I))**2
      FU = EU
      MK = 0
17    CALL SOLVE(I,FIRST)
*
*  *****  IF KK EQUALS 3, OMIT THE NODE CHECKING
*
      IF (KK .EQ. 3) GO TO 51
*
*  *****  COUNT THE NUMBER OF NODES
*
      MN = M
      NC = NODEC(MN)
      IF (TRACE) WRITE(OUT,99) EL(I),NC,MN,NJ,PDE(MN),ED,EU,EM,DELTAE
99    FORMAT(2X,A3,' NC =',I3,' MN =',I3,' NJ =',I3,' PDE(MN) =',
     1   D10.2,' ED =',D10.2,' EU =',D10.2,' EM =',D10.2,
     2   ' DELTAE =',D10.2)
*
*  *****  IF NODE COUNT IS OFF BY NO MORE THAN 1 AND DELTAE IS STILL
*  *****  QUITE LARGE, APPLY THE DELTAE CORRECTION
*
      IF (IABS(NC-NODE) .EQ. 1 .AND. DABS(DELTAE/ED) .GT. 0.02D0)
     1      GO TO 46
*
*  *****  BRANCH ACCORDING TO WHETHER THE NODE COUNT IS TOO SMALL,
*  *****  JUST RIGHT, OR TOO LARGE
*
12    IF (NC - NODE ) 8,9,10
*
*  *****  THE SOLUTION HAS THE CORRECT NUMBER OF NODES
*
9     V2 = DABS(DELTAE)/ED .LT. 1.D-5
      IF (PDE(MN) .LT. D0 .AND. .NOT. V2) GO TO 46
      IF (PDE(MN) .GT. D0) GO TO 51
      DO 52 J = 1,NO
52    PDE(J) = - PDE(J)
      PP = -D2 - PP
51    AZZ = AZD*(D1 + PP)
      CALL EIJSET(I,I,ED)
      RETURN
*
*  *****  THE SOLUTION HAS TOO FEW NODES
*
8     IF (PDE(MN) .LE. D0) GO TO 11
      DEL = D1 - ED/EU
      EU = ED
      IF ( DEL .LT. .05D0) FU = FU*((L(I)+1+NC)/FN)**2.5
       IF (DEL  .GE. .05D0) FU = ED*((L(I)+1+NC)/FN)**2.5
      IF (FU .LT. EM) FU = D5*(EU + EM)
      IF (DABS(FU - ED) .LT. 0.001D0) GO TO 27
      ED = FU
      GO TO 33
*
*  *****  TRY A NEW VALUE OF ED WHICH MUST LIE WITHIN THE UPPER AND
*  *****  LOWER BOUND
*
11    EDP = ED
                    ED = ED*((L(I)+1+NC)/FN)**2.5
      IF (ED .GE. EU ) ED = D5*(EU + EDP)
      IF (ED .LE. EM ) ED = D5*(EM + EDP)
33    MK = MK + 1
      IF ( EU .LE. EM ) WRITE(OUT,30) EM,EU,ED
30    FORMAT(6X,48HWARNING: DIFFICULTY WITH NODE COUNTING PROCEDURE/
     1   6X,42HLOWER BOUND ON ED GREATER THAN UPPER BOUND/
     2   6X,5HEL = ,F10.6,7H  EU = ,F10.6,7H  ED = ,F10.6)
      FIRST = .FALSE.
      IF ( MK .GT. 3*N(I) .OR. EU-EM .LT. FN**(-3)) GO TO 27
      GO TO 17
*
*  *****  THE SOLUTION HAS TOO MANY NODES
*
10    IF (PDE(MN) .LT. D0) GO TO 11
      DEL = D1 - EM/ED
      EM = ED
      IF (DEL .LT. 0.05D0) FM = FM*((L(I)+1+NC)/FN)**2.5
      IF (DEL .GE. 0.05D0) FM = ED*((L(I)+1+NC)/FN)**2.5
      IF (FM .GT. EU) FM = D5*(EU + EM)
      IF (DABS(FM - ED) .LT. 0.001D0) GO TO 27
      ED = FM
       GO TO 33
*
*  *****  ADJUST ENERGY TO LIE BETWEEN UPPER AND LOWER BOUND
*
46    ED = ED - DELTAE
      IF ( ED .GE. EM .AND. ED .LE. EU ) GO TO 33
      EDP = ED
      IF ( NC-NODE .NE. 0 ) ED = (ED+DELTAE)*((L(I)+1+NC)/FN)**2.5
      IF ( ED .GE. EM .AND. ED .LE. EU ) GO TO 33
      ED = EDP + DELTAE + DELTAE
      IF ( ED .GE. EM .AND. ED .LE. EU ) GO TO 33
      ED = ED -DELTAE
      DELTAE = D5*DELTAE
      GO TO 46
*
*  *****  METHOD WAS UNABLE TO FIND AN ACCEPTABLE SOLUTION
*
27    WRITE(OUT,28) KK,EL(I),NC,NJ,ED,EM,EU
28    FORMAT(10X,6HMETHOD,I2,38H UNABLE TO SOLVE EQUATION FOR ELECTRON,
     1   A3/10X,5HNC = ,I3,3X,5HNJ = ,I3,3X,5HED = ,F10.6,3X,5HEL = ,
     2   F10.6,3X,5HEU = ,F10.6)
      FAIL = .TRUE.
      RETURN
      END
*
*     ------------------------------------------------------------------
*              N M R V S
*     ------------------------------------------------------------------
*
*       Given two starting values, PDE(1) and PDE(2), values of PDE(j),
*   j=3,4,...,NJ+1 are obtained by outward integration of
*               Y" = YR y + F
*   using the discretization  of  Eq.  (6-27 )  with  the  difference
*   correction.  With PDE(NJ) given, the tail procedure is applied to
*   PDE(j),j=NJ+1,  NJ+2,...,MM, where MM is determined automatically
*   and DELTA is the difference between  PDE(NJ+1)  for  outward  and
*   inward integration. (See Eq 6-32, 6-33, and 6-37 for further
*   details.)
*
*
      SUBROUTINE NMRVS(NJ,DELTA,MM,PP,F)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220,NCD=100)
*
        INTEGER IN,OUT,ERR,PRI,OUC,OUD,OUF,OUH
        COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,OUC,OUD,OUF,OUH

      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
      DIMENSION PP(NOD),F(NOD),A(150),D(150)
      EQUIVALENCE (G,G3)
*
*  *****  INTEGRATE OUTWARD TO NJ+1
*
      Y1 = PP(1)
      Y2= PP(2)
      G1 = YR(1)
      G2 = YR(2)
      M = NJ + 1
      DO 1 I = 3,M
      G3 = YR(I)
      Y3 = (Y2+Y2-Y1 + (D10*G2*Y2 + G1*Y1) + F(I-1)) / (D1 - G3)
      PP(I) = Y3
      Y1 = Y2
      Y2 = Y3
      G1 = G2
1     G2 = G3
      DELTA = Y3
*
*  *****  APPLY THE TAIL PROCEDURE
*
      K = 1
      PP(M) = -(D1 - G1)*Y1 + F(M)
      A(1) = D1 - G
      D(1) = -(D2 + D10*G)
22    RATIO = A(K)/D(K)
*
*  *****  THE INTEGER 149 IN THE NEXT STATEMENT IS THE DIMENSION OF A
*  *****  MINUS 1
*
      IF (K .GE. (150)-1 .OR. M .EQ. ND) GO TO 23
      K = K +1
      M = M+1
      G = YR(M)
      A(K) = D1 - G
      D(K) = -(D2 + D10*G) - A(K)*RATIO
      PP(M) = -PP(M-1)*RATIO + F(M)
      IF (DABS(PP(M))+DABS(PP(M-1)) .GT. TOL .OR. K .LT. 9) GO TO 22
20    CONTINUE
      PP(M) = PP(M)/D(K)
      J = M+1
      DO 2 I= J,NO
2     PP(I) = D0
      DO 3 J = 2,K
      I = M-J+1
      II = K-J+1
3     PP(I) = (PP(I)-A(II+1)*PP(I+1))/D(II)
*
*  *****  SET DELTA = DIFFERENCE OF THE TWO SOLUTIONS AT NJ+1
*  *****         MM = NUMBER OF POINTS IN THE RANGE OF THE SOLUTION
*
      DELTA = DELTA - PP(I)
      MM = M
      RETURN
23    WRITE(OUT,24)
24    FORMAT(6X,52HWARNING: FUNCTIONS TRUNCATED BY NMRVS IN TAIL REGION)
      GO TO 20
      END
*
*     ------------------------------------------------------------------
*              N O D E C
*     ------------------------------------------------------------------
*
*      Counts the number of nodes of the function PDE(j) in the range
*   j = 40,...,M-10.   The node counting procedure counts the local max
*   and min values.   Only nodes between sufficiently large max and
*   min values are counted.
*
*
      FUNCTION NODEC(M)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=100)
*
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
*   ***** FIND MAX|PDE(J)|
*
       MM = M - 10
      DM = 0.D0
      DO 1 J = 40,MM
1     DM = DMAX1(DM, DABS(PDE(J)))
*
*   *****  COUNT THE NUMBER OF LOCAL MAX OR MIN'S
*
      NCC = 0
      SIGN = 0.D0
      DIFF1 = PDE(40) - PDE(39)
      DO 2 J = 40, MM
      DIFF2 = PDE(J+1) - PDE(J)
      IF (DIFF2*DIFF1 .GT. 0.D0 .OR. DIFF1 .EQ. 0.D0) GO TO 2
*
*   *****  A MAX OR MIN HAS BEEN FOUND.   TEST IF IT IS
*          SUFFICIENTLY LARGE
*
      IF ( DABS(PDE(J))/DM .LT. 0.05D0 ) GO TO 2
*
*   ***** CHECK IF THIS IS THE FIRST SIGNIFICANT MAXIMUM
*
      IF (SIGN .NE. 0.D0 ) GO TO 4
      M = J
      GO TO 3
*
*   ***** IF NOT THE FIRST, TEST WHETHER A SIGN CHANGE HAS
*         OCCURRED SINCE THE LAST SIGNIFICANT MAX OR MIN
*
4     IF (PDE(J)*SIGN .GT. 0.D0 ) GO TO 2
      NCC = NCC + 1
*
*   ***** RESET FOR THE NEXT NODE
*
3     SIGN = PDE(J)
2     DIFF1 = DIFF2
      NODEC = NCC
      RETURN
      END
*     MCHF_ MCHF (Part 2 of 2)
*-----------------------------------------------------------------------
*               O 2 P E R T
*-----------------------------------------------------------------------
*
*     Contribution to the energy variation from a pair of overlap
*   integrals
*
      SUBROUTINE O2PERT(IIN,JIN,II,INT,DC,DV)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=100)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
        INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
        PARAMETER (IDIM=550,NCDIM=3000)
        COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :         ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
        LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
      DIMENSION IND(6),KV(3),OV(3)
      EQUIVALENCE (I1,IND(1)),(I2,IND(2)),(I3,IND(3)),(I4,IND(4)),
     :            (I5,IND(5)),(I6,IND(6))
*
        I = JIN
        J = IIN
        I1 = IEL(II,1)
        I2 = IEL(II,2)
        I3 = IEL(II,3)
        I4 = IEL(II,4)
        KK = KVAL(II)
        KV(1) = KK/64
        KV(2) = KK - 64*KV(1)
        OV(1) = QUADR(I1,I2,0)
        OV(2) = QUADR(I3,I4,0)
        OVL = OV(1)**KV(1)*OV(2)**KV(2)
        IF (ABS(OV(1)).LE.1.E-8 .OR. ABS(OV(2)).LE.1.E-8) THEN
           WRITE(6,20) 
     :             'Overlap Integrals too small for this code:',
     :             'Electrons',I1,I2,' have overlap =',OV(1),
     :             'Electrons',I3,I4,' have overlap =',OV(2)
 20        FORMAT(/A/2(1X,A,2I3,A,F8.5/))
           STOP
        END IF
1     DO 10 KP = 1,2
      DI = D0
      DII = D0
      DIJ = D0
      DO 2 K = 1,4
      IF (IND(K) .NE. I) GO TO 2
        KI = (K-1)/2 + 1
        I5 = J
        I6 = IND(K+1 - 2*MOD(K+1,2))
        OV(3) = QUADR(I5,I6,0)
        KV(3) = 1
        CK = KV(KI)
        TI = CK*OVL/OV(KI)*OV(3)
        IF (.NOT. ALL) GO TO 3
        KV(KI) = KV(KI)-1
      DO 4 K2 = 1,6
      IF (IND(K2) .NE. J) GO TO 5
        KIJ = (K2-1)/2 + 1
        IF (KV(KIJ) .EQ. 0) GO TO 5
        IJPAIR = K2+1 - 2*MOD(K2+1,2)
        CKIJ = KV(KIJ)
        DIJ = DIJ + CKIJ*TI/OV(KIJ)*QUADR(I,IND(IJPAIR),0)
5     IF (IND(K2) .NE. I) GO TO 4
        KII = (K2-1)/2 + 1
        IF (KV(KII) .EQ. 0) GO TO 4
        IIPAIR = K2+1 - 2*MOD(K2+1,2)
        CKII = KV(KII)
        DII = DII + CKII*TI/OV(KII)*QUADR(J,IND(IIPAIR),0)
4     CONTINUE
      KV(KI) = KV(KI) + 1
3     DI = DI + TI
2     CONTINUE
6     IF (KP .EQ. 2) THEN
          DC = DI - DJ
          DV = DII + DJJ - DIJ - DJI
          RETURN
      ELSE
         DJ = DI
         DJJ = DII
         DJI = DIJ
         I = IIN
         J = JIN
      END IF
10    CONTINUE
      END
*
*     ------------------------------------------------------------------
*              O R T H O G
*     ------------------------------------------------------------------
*
*       This routine orthogonalizes the set of radial functions when an
*   orthogonality constraint applies.  A Gram-Schmidt type of  process
*   is used.  When more than one radial function with a given (nl) is
*   present, it may be necessary to solve a 2x2 system of equations.
*
*
      SUBROUTINE ORTHOG
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
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
        LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
        LOGICAL DIAG
        COMMON AC(30,NWD),BC(NWD)
*
      IF (NWF .EQ. 1 .OR. IB .GT. NWF) RETURN
      II = MAX0(2,IB)
      DO 2 I = II,NWF
         DIAG = .TRUE.
         IBEGIN = IEPTR(I-1)+1
         IP = IBEGIN
         IJ = 0
 60      JV = IJE(IP)
         IF (JV .LT. I .AND. IP .LE. IEPTR(I)) THEN
            IJ = IJ+1
            IF ( IJ .GT. (NWD)) STOP ' TOO MANY ORTHOGONALITY CONDITIONS'
            BC(IJ) = QUADR(I,JV,0)
            AC(IJ,IJ) = D1
            DO 62 JJ = IBEGIN,IP-1
                   IK = JJ - IBEGIN + 1
               IF (E(IJE(IP),IJE(JJ)) .NE. D0 ) THEN
                  AC(IJ,IK) = D0
                  AC(IK,IJ) = D0
                ELSE
                  AC(IJ,IK) = QUADR(IJE(IP),IJE(JJ),0)
                  AC(IK,IJ) = AC(IJ,IK)
                  DIAG = .FALSE.
               END IF
62          CONTINUE
            IP = IP+1
            GO TO 60
         END IF
      IF ( IJ .GT. 0) THEN
         IF ( .NOT. DIAG .AND. IJ.GT.1) CALL LINEQN(30,IJ,AC,BC)
         M = MAX(I)
         AZZ = AZ(I)
         IP = IBEGIN
         CTOTAL = D0
         DO 65 JJ = 1,IJ
            C = BC(JJ)
            IF (DABS(C) .GT. 1.D-10) THEN
               WRITE(OUT,63) EL(IJE(IP)),EL(I),C
63             FORMAT(6X,'<',A3,'|',A3,'>=',1PD8.1)
               M = MAX0(M,MAX(IJE(IP)))
               DO 64 J = 1,M
                  P(J,I) = P(J,I) - C*P(J,IJE(IP))
64             CONTINUE
               AZZ = AZZ - C*AZ(IJE(IP))
            END IF
            IP = IP + 1
            CTOTAL = CTOTAL + ABS(C)
65       CONTINUE
         IF (CTOTAL .GT. 1.D-10) THEN
            PNN = DSQRT(QUADR(I,I,0))
            DO 66 JJ = 1,M
                  P(JJ,I) = P(JJ,I)/PNN
66          CONTINUE
            AZZ = AZZ/PNN
            M = NO
67          IF (DABS(P(M,I)) .LT. 1.D-15) THEN
                  P(M,I) = D0
                  M = M-1
                  GO TO 67
            END IF
            MAX(I) = M
            AZ(I) = AZZ
            VARIED(I) = .TRUE.
         END IF
      END IF
2     CONTINUE
      END
*
*     ------------------------------------------------------------------
*              O U T P U T
*     ------------------------------------------------------------------
*
*       The radial functions and orthogonality integrals are printed,
*   if PRINT is .TRUE.   The  functions  will  also  be  punched  (or
*   stored) on unit OUF, if OUF .NE. 0.
*
*
      SUBROUTINE OUTPUT(PRINT)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220,NCD=100)
        INTEGER MMX
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
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
        LOGICAL PRINT
      DIMENSION POUT(8)
      IF ( .NOT. PRINT ) GO TO 31
*
*  *****  PRINT RADIAL FUNCTIONS, 7 PER PAGE
*
      ML = 1
2     MU = MIN0(ML+7,NWF)
      I = MU - ML + 1
      MX = 0
      DO 1 J = ML,MU
1     MX = MAX0(MX,MAX(J))
      WRITE(PRI,5) ATOM,TERM,(EL(J),J=ML,MU)
5     FORMAT(1H1,9X,19HWAVE FUNCTIONS FOR  ,2A6//10X,1HR,8(10X,A3))
      K= 0
      KK = 0
      DO 6 J = 1,MX
      DO 9 JJ = ML,MU
      IJ = JJ - ML + 1
9     POUT(IJ) = P(J,JJ)*R2(J)
      K = K+1
      IF (K .LE. 10) GO TO 6
      K = 1
      KK = KK+1
      IF (KK .LT. 5) GO TO 21
      KK = 0
      WRITE(PRI,23)
23    FORMAT(1H1//)
      GO TO 6
21    WRITE(PRI,8)
8     FORMAT(1X)
6     WRITE(PRI,10) R(J),(POUT(JJ),JJ=1,I)
10    FORMAT(F13.5,F15.6,7F13.6)
      DO 15 J = ML,MU
      IJ = J - ML + 1
15    POUT(IJ) = DPM(J)
      WRITE(PRI,16) (POUT(J),J=1,I)
16    FORMAT(4X,10HMAX. DIFF. ,F15.7,7F13.7)
      ML = ML+8
      IF (ML .LE. NWF) GO TO 2
31    IF ( NWF .LE. 1) GO TO 30
*
*  *****  PRINT OVERLAP INTEGRALS
*
      WRITE(PRI,11) ATOM,TERM
11    FORMAT(////10X,'OVERLAP INTEGRALS FOR ATOM ',A6,' TERM ',A6
     1   //20X, 4H(NL),3X,4H(NL),7X,8HINTEGRAL //)
      LM = IB
      ML = MAX0(2,LM)
      DO 12 I = ML,NWF
      JF = I - 1
      DO 13 J = 1,JF
      IF (L(I) .NE. L(J)) GO TO 13
      T = QUADR(I,J,0)
      IF (ABS(T) .GT. 1.D-5) WRITE(PRI,17) EL(I),EL(J),T
17     FORMAT(21X,A3,4X,A3,F15.8)
13    CONTINUE
12    CONTINUE
30    IF ( OUF .EQ. 0) GO TO 14
*
*  *****  OUTPUT FUNCTIONS ON UNIT OUF FOR FUTURE INPUT
*
*         EKI retained only for compatibility with MCHF format
*
      DO 3 I = 1,NWF
      MMX = MAX(I)
      WRITE (OUF) ATOM,TERM,EL(I),MMX,Z,E(I,I),EKI,AZ(I),
     :   (P(J,I),J=1,MMX)
3     CONTINUE
*
14    RETURN
      END
*
*-----------------------------------------------------------------------
*               P E R T
*-----------------------------------------------------------------------
*       Contribution to the perturbation of the energy under a rotation
*   from either an Fk, Gk, or L integral.

      SUBROUTINE PERT(IIN,JIN,II,INT,DC,DV)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220,NCD=100)
*
        CHARACTER CONFIG*40,EL*3,ATOM*6,TERM*6,COUPLE*3
        COMMON /LABEL/CONFIG(NCD),EL(NWD),ATOM,TERM,COUPLE(NCD,9)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
      INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
      PARAMETER (IDIM=550,NCDIM=3000)
      COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :         ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
      LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :       ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
      I = IIN
      J = JIN
      DC = D0
      DV = D0
      I1 = IEL(II,1)
      I2 = IEL(II,2)
      K = KVAL(II)
      IF (I1 .GT. I2) THEN
         IT = I1
         I1 = I2
         I2 = IT
      END IF
      IF (I1.EQ.I .AND. I2.EQ.J) THEN
         IF (INT .LE. 2) THEN
        DC = 2*(RK(I,J,J,J,K,REL)-RK(J,I,I,I,K,REL))
              IF (ALL)   DV = 2*(FK(J,J,K,REL)+FK(I,I,K,REL)
     :               -2*FK(I,J,K,REL)-4*GK(I,J,K,REL))
         ELSE
        DC = HLC(EL,J,J,REL)-HLC(EL,I,I,REL)
        IF (ALL) DV = -4*HLC(EL,I,J,REL)
         END IF
      ELSE IF ((I1-I)*(I2-I)*(I1-J)*(I2-J).NE.0.AND.INT.EQ.6) THEN
*        ... i is a closed subshell
        LI = L(I)
        DC =D2*RK(I,I1,J,I2,0,REL)
        IF (ALL) DV = D2*(RK(J,I1,J,I2,0,REL)-RK(I,I1,I,I2,0,REL))
        DO 18 K = IABS(LI-L(I1)),LI+L(I1),2
           CIJ = CB(LI,L(I1),k)
           IF (I1.EQ.I2) THEN
              DC = DC - D2*CIJ*RK(I,I1,I2,J,K,REL)
           ELSE
              DC = DC-CIJ*(RK(I,I1,I2,J,K,REL)+RK(J,I1,I2,I,K,REL))
           END IF
           IF (ALL) DV = DV -
     :        D2*CIJ*(RK(J,I1,I2,J,k,REL) - RK(I,I1,I2,I,k,REL))
18      CONTINUE
        DC = -D2*SUM(I)*DC
        IF (ALL) DV =  - D2*SUM(I)*DV
      ELSE
         PHASE = D1
         DO 10 KK = 1,2
           IF (I1.EQ.I .OR. I2.EQ.I) THEN
        IF (I1.NE.I) THEN
           IT = I1
           I1 = I2
           I2 = IT
        END IF
        IF (I1.EQ.I2) THEN
           DC = 4*PHASE*RK(J,I,I,I,K,REL)
           IF (ALL) DV = 4*(FK(J,I,K,REL)
     :                           -FK(I,I,K,REL)+2*GK(I,J,K,REL))
        ELSE IF (INT.EQ.1) THEN
           DC = 2*PHASE*RK(I,I2,J,I2,K,REL)
           IF (ALL) DV = 2*(FK(J,I2,K,REL)-FK(I,I2,K,REL))
        ELSE IF (INT.EQ.2) THEN
           DC = 2*PHASE*RK(I,I2,I2,J,K,REL)
           IF (ALL) DV = 2*(GK(J,I2,K,REL)-GK(I,I2,K,REL))
        ELSE IF (INT.EQ.6) THEN
           DC = PHASE*HLC(EL,J,I2,REL)
           IF (ALL) DV = -HLC(EL,I,I2,REL)
        ELSE
           OV1 = QUADR(J,I2,0)
           OV2 = QUADR(I,I2,0)
           DC = K*PHASE*OV1*OV2**(K-1)
           IF (ALL) THEN
              DV = -K*OV2**K
              IF (K.GT.1) DV = DV+K*(K-1)*OV2**(K-2)*OV1**2
           END IF
        END IF
           END IF
           I = JIN
           J = IIN
           PHASE = -D1
 10     CONTINUE
      END IF
      END
*
*     ------------------------------------------------------------------
*              Q U A D
*     ------------------------------------------------------------------
*
*       Evaluates the integral of F(r)G(r) with respect to r , where
*   F(r) and G(r) have the same asymptotic properties as P (r).   The
*                                                         i
*   composite Simpson's rule is used.   The integrand is zero for r >
*   r  .
*    M
*
      DOUBLE PRECISION FUNCTION QUAD(I,M,F,G)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
        DIMENSION F(NOD),G(NOD)
*
      D = (D1 + D5*Z*R(1))/(H1*(2*L(I) + 3))
      QUAD = RR(1)* F(1)*G(1)*( D -D5)
      QUAD2 = D0
      DO 1 J = 2,M,2
      QUAD = QUAD + RR(J-1)*F(J-1)*G(J-1)
      QUAD2 = QUAD2 + RR(J)*F(J)*G(J)
1     CONTINUE
      QUAD = H1*(QUAD + D2*QUAD2)
      RETURN
      END
*
*     ------------------------------------------------------------------
*              P O T L
*     ------------------------------------------------------------------
*
*       Computes and stores the potential function
*                              2(k-1)
*              YR = SUM  a    Y      (j,j;r)
*                   j,k   ijk
*
        SUBROUTINE POTL(I)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220,NCD=100)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
        INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
        PARAMETER (IDIM=550,NCDIM=3000)
        COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :         ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
        LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
      DO 1 J=1,NO
1     YR(J) = D0
        DO 2 J = 1,NWF
           IF (I.GT.NCLOSD .AND. J.GT.NCLOSD) GO TO 2
           C = SUM(J)
           IF ( I.EQ.J ) C = C - D1
           CALL YKF(J,J,0,REL)
           DO 3 JJ = 1,NO
              YR(JJ) = YR(JJ) + C*YK(JJ)
3          CONTINUE
           IF ( I.EQ.J .AND. L(I) .GT. 0) THEN
              DO 4 K = 2,2*L(I),2
                 CC = -C*CA(L(I),K)
                 CALL YKF(I,I,K,REL)
                 DO 5 JJ = 1,NO
                    YR(JJ) = YR(JJ) + CC*YK(JJ)
5                CONTINUE
4             CONTINUE
           END IF
2       CONTINUE
*
        SUMI = SUM(I)
        IBEGIN = 1
        IEND = INTPTR(1)
        DO 10 J = IBEGIN,IEND
           IE = 0
           IF (IEL(J,1) .EQ. I) THEN
              IE = IEL(J,2)
           ELSE IF (IEL(J,2) .EQ. I) THEN
              IE = IEL(J,1)
           END IF
           IF (IE .NE. 0) THEN
              C = COEF(J)/SUMI
              IF (IEL(J,1) .EQ. IEL(J,2)) C = 2*C
              CALL YKF(IE,IE,KVAL(J),REL)
              DO 12 JJ = 1,NO
                 YR(JJ) = YR(JJ) + C*YK(JJ)
 12           CONTINUE
           END IF
 10     CONTINUE
        END
*
*   --------------------------------------------------------------------
*               R E O R D
*   --------------------------------------------------------------------
*
*       Reorder the list of first appearance so that the functions to be
*   iterated appear last in the list.
*
        SUBROUTINE REORD(EL, ELC, NWF, IERR)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220,NCD=100)
        CHARACTER*3 EL(NWD), ELC
*
*
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
        COMMON ZZ(NWD),IND(NWD),IELI(5),NOCCSH(NCD)
*
        IERR = 1
        CALL EPTR(EL, ELC, I, *99)
*
*       The i'th orbital is to be placed at the end
*
        SI = S(I)
        METHI = METH(I)
        ACCI = ACC(I)
        INDI = IND(I)
        DO 10 J = I, NWF-1
           EL(J) = EL(J+1)
           S(J) = S(J+1)
           METH(J) = METH(J+1)
           ACC(J) = ACC(J+1)
           IND(J) = IND(J+1)
10      CONTINUE
        EL(NWF) = ELC
        S(NWF) = SI
        METH(NWF) = METHI
        ACC(NWF) = ACCI
        IND(NWF) = INDI
        IERR = 0
99      RETURN
        END
*
*-----------------------------------------------------------------------
*               R O T A T E
*-----------------------------------------------------------------------
*     Rotate orbitals connected through orthogonality.
*
        SUBROUTINE ROTATE
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
        LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
        LOGICAL SETEQL
*
*  *****  ROTATE PAIRS CONNECTED BY ORTHOGONALITY BUT NOT WHEN ONE OF
*         THE ORBITALS IS SIMULTANEOUSLY ORTHOGONAL TO A NON-ORTHOGONAL
*         PAIR
*
        DO 1 I = IB,NWF-1
           DO 2 J = I+1,NWF
              IF (DABS(E(I,J)) .GT. 1.D-10 .AND. SETEQL(I,J)
     :            .AND. DABS(V(I,J)) .GT. 1.D-4) THEN
                 EPS = E(J,I)
                 DD = DSQRT(D1 + EPS*EPS)
                 DO 41 JJ = 1,NO
                    PI = (P(JJ,I) + EPS*P(JJ,J))/DD
                    P(JJ,J) = (P(JJ,J) - EPS*P(JJ,I))/DD
41               P(JJ,I) = PI
                 VARIED(I) = .TRUE.
                 VARIED(J) = .TRUE.
              END IF
2          CONTINUE
1       CONTINUE
        END
*
*     ------------------------------------------------------------------
*              R P E R T
*     -----------------------------------------------------------------
*
*       RPERT determines the effect of a perturbation in the form of a
*   rotation of the i'th and j'th orbital both on the energy and on the
*   stationary condition, due to the presence of a given RK integral.
*
*
      SUBROUTINE RPERT(IIN,JIN,II,INT,DC,DV)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=100)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
        INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
        PARAMETER (IDIM=550,NCDIM=3000)
        COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :         ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
        LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
      DIMENSION IND(4)
      EQUIVALENCE (I1,IND(1)),(I2,IND(2)),(I3,IND(3)),(I4,IND(4))
*
        I = JIN
        J = IIN
        I1 = IEL(II,1)
        I2 = IEL(II,2)
        I3 = IEL(II,3)
        I4 = IEL(II,4)
        KK = KVAL(II)
        INC = 1
        IF (I1.EQ.I2 .AND. I3.EQ.I4) INC = 2
1     DO 10 KP = 1,2
      DI = D0
      DII = D0
      DIJ = D0
      DO 2 K = 1,4,INC
      IF (IND(K) .NE. I) GO TO 2
      IND(K) = J
      DI = DI + RK(I1,I2,I3,I4,KK,REL)
      IF (.NOT. ALL) GO TO 3
      DO 4 K2 = 1,4
      IF (IND(K2) .NE. J) GO TO 5
      IND(K2) = I
      DIJ = DIJ + RK(I1,I2,I3,I4,KK,REL)
      IND(K2) = J
5     IF (IND(K2) .NE. I) GO TO 4
      IND(K2) = J
      DII = DII + RK(I1,I2,I3,I4,KK,REL)
      IND(K2) = I
4     CONTINUE
3     IND(K) = I
2     CONTINUE
6     IF (KP .EQ. 2) THEN
          DC = INC*(DI - DJ)
          DV = INC*(DII + DJJ - DIJ - DJI)
          RETURN
      ELSE
         DJ = DI
         DJJ = DII
         DJI = DIJ
         I = IIN
         J = JIN
      END IF
10    CONTINUE
      END
*
*     ------------------------------------------------------------------
*              S C A L E
*     ------------------------------------------------------------------
*
*       The current radial functions are scaled according to 
*   Z-dependent screening.   Values of AZ and E(I,I), the starting
*   values and the diagonal energy parameters are also scaled.
*
      SUBROUTINE SCALE(ZZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
      COMMON RATIO,SR,SC,SS,F0,F1,F2,F3,PNORM,THETA,K
      DIMENSION RS(NOD),PS(NOD)
      EQUIVALENCE (RS(1),YR(1)),(PS(1),X(1))
*
*  *****  SCALE VALUES OF R=RS, P=PS AND ONE-ELECTRON PARAMETERS.
*  *****  GENERATE NEW VALUES OF R, R*R, AND SQRT(R)
*
      RATIO = Z/ZZ
      SR = DSQRT(RATIO)
      DO 1 J = 1,NO
      R(J) = R(J)*RATIO
      RR(J) = R(J)*R(J)
1     R2(J) = R2(J)*SR
      DO 2 I = 1,NWF
      SC = (ZZ-S(I))/(Z-S(I))
      SS = SC*RATIO
      ED = E(I,I)*SC**2
      CALL EIJSET(I,I,ED)
      DO 3 J = 1,NO
      RS(J) = R(J)/SS
3     PS(J) = P(J,I)*SC
      SC = (ZZ - D5*S(I))/(Z - D5*S(I))
      AZ(I) = AZ(I)*SC**(L(I)+1)*DSQRT(SC)
      K = 3
*
*  *****  INTERPOLATE THE (RS,PS) FUNCTIONS FOR VALUES OF P AT THE SET
*  *****  OF POINTS R
*
      DO 4 J = 1,NO
*
*  *****  SEARCH FOR THE NEAREST ENTRIES IN THE (RS,PS) TABLE
*
5     IF (K .EQ. ND) GO TO 7
      IF (RS(K) .GT. R(J)) GO TO 6
      K = K + 1
      GO TO 5
*
*  *****  INTERPOLATE
*
6     THETA = DLOG(R(J)/RS(K-1))/H
      F0 = PS(K-2)
      F1 = PS(K-1)
      F2 = PS(K)
      F3 = PS(K+1)
      P(J,I) = D5*(F1+F2) + (THETA -D5)*(F2 - F1) +
     1   THETA*(THETA - D1)*(F0 - F1 - F2 + F3)/D4
      GO TO 4
7     P(J,I) = D0
4     CONTINUE
      MAX(I) = NO
*
*  *****NORMALIZE THE INTERPOLATED FUNCTION
*
      PNORM = DSQRT(QUADR(I,I,0))
      DO 10 J = 1,NO
10    P(J,I) = P(J,I)/PNORM
2     CONTINUE
      Z = ZZ
      RETURN
      END
*
*     ------------------------------------------------------------------
*              S C F
*     -----------------------------------------------------------------
*
*       This routine controls the SCF procedure described in Chapter
*   7.  If certain input parameters are zero (or blank) they will  be
*   set to their default value.
*
*          Parameter       Default Value
*          --------        -------------
*          CFGTOL          1.D-10
*          SCFTOL          1.D-7
*          IC              (NWF + 1 - IB)/4 + 3
*          NSCF            12
*
*   The self-consistency convergence criterion is
*
*          Z2 = SQRT( SCFTOL*(Z*NWF/2) )
*
*   It is increased by a factor two at the end of each iteration whereas
*   CFGTOL is increased by SQRT(2).
*
*
      SUBROUTINE SCF(ACFG,SCFTOL,CFGTOL,LD)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220,NCD=100)
*
        INTEGER IN,OUT,ERR,PRI,OUC,OUD,OUF,OUH
        COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,OUC,OUD,OUF,OUH
*
        CHARACTER CONFIG*40,EL*3,ATOM*6,TERM*6,COUPLE*3
        COMMON /LABEL/CONFIG(NCD),EL(NWD),ATOM,TERM,COUPLE(NCD,9)
*
        COMMON /MATRIX/ETOTAL,W(NCD,NCD)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
        LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
        LOGICAL LAST,LD,CONV,ECONV
        CHARACTER ANS
*
*  *****  SET THE SCF CONVERGENCE PARAMETER TO AN OPTIMISTIC VALUE
*
*     TOL = DSQRT(Z)*1.D-7
      TOL = DSQRT(Z)*1.D-10
      Z2 = SCFTOL*DSQRT(Z*NWF)
      WRITE(OUT,15)
15    FORMAT(//)
      WRITE(OUT,16) OMIT,ACFG,SCFTOL,NO,REL
   16 FORMAT(10X,44HWEAK ORTHOGONALIZATION DURING THE SCF CYCLE=,L4/
     :       10X,44HACCELERATING PARAMETER FOR MCHF ITERATION  =,F5.2/
     :       10X,44HSCF CONVERGENCE TOLERANCE (FUNCTIONS)      =,1PD9.2
     :      /10X,44HNUMBER OF POINTS IN THE MAXIMUM RANGE      =,I4/
     :       10X,44HRELATIVISTIC DIAGONAL  ENERGY CORRECTIONS  =,L4)
*
*  *****  SET ITERATION PARAMETERS
*
      IPR = 0
      ECONV = .FALSE.
      LAST = .FALSE.
      DP1 = D0
      ETOTAL = D0
      EC = D0
      ICYCLE = 0
      CALL UPDATE
      IF ( IB .GT. NWF ) GO TO 17
      IF ( .NOT. LD ) GO TO 9
19    IF ( ID .EQ. 1 .OR. NCFG .EQ. 1) GO TO 9
      CALL DIAG(ECONV,ACFG,CFGTOL,LAST)
*
*  *****  PERFORM NSCF SELF-CONSISTENT FIELD ITERATIONS
*
9     DO 100 I = 1,NSCF
      ICYCLE = ICYCLE + 1
      WRITE(OUT,7) ICYCLE,CFGTOL,Z2
7     FORMAT(//10X,17HITERATION NUMBER ,I2/10X,16H----------------//
     1 10X,50HCONVERGENCE CRITERIA:ENERGY  (CFGTOL)            =,1PD9.1/
     2 11X,49H                   :FUNCTION(SCFTOL*SQRT(Z*NWF))=,1PD9.1/)
      DP1 = D0
      IF (IB .GT. NWF) GO TO 17
      CALL GRANGE
*
*  *****  SOLVE EACH DIFFERENTIAL EQUATION IN TURN
*
      WRITE(OUT,14)
14    FORMAT(/20X,' EL',9X,'ED',13X,'AZ',11X,'NORM',7X,'DPM')
      DO 2 J = IB,NWF
      CALL DE(J)
      IF ( FAIL ) RETURN
      DP = DPM(J)*DSQRT(SUM(J))
      IF ( DP1 .GE. DP ) GO TO 2
      DP1 = DP
      JJ = J
2     CONTINUE
      IF ((NCFG .EQ. 1 .OR. ID .EQ. 1) .AND. DP1 .LT. Z2) GO TO 6
      IF ( IC .LE. 0) GO TO 6
*
*  *****  SOLVE IC DIFFERENTIAL EQUATIONS EACH TIME SELECTING THE
*  *****  ONE WITH THE LARGEST DPM
*
      DO 4 II =1,IC
      CALL DE(JJ)
      IF ( FAIL ) RETURN
      DP1 = D0
      DO 5 J = IB,NWF
      DP = DSQRT(SUM(J))*DPM(J)
      IF ( DP1 .GT. DP ) GO TO 5
      JJ = J
      DP1 = DP
5     CONTINUE
      IF (DP1 .LT. Z2) GO TO 6
4     CONTINUE
6     CALL ORTHOG
      CALL UPDATE
      IF ( LAST ) GO TO 17
      IF ( I .EQ. NSCF ) GO TO 1
      IF (.NOT.(NCFG .EQ. 1 .OR. ID .EQ. 1)) GO TO 12
      IF (DP1 .LE. Z2 )  LAST = .TRUE.
      GO TO 1
12    CALL DIAG(ECONV,ACFG,CFGTOL,LAST)
*
*  *****  IF FUNCTIONS APPEAR TO HAVE CONVERGED,SOLVE EACH AGAIN, IN
*  *****  TURN, AND TEST AGAIN
*
      CONV = ECONV .AND. DP1 .LE. Z2
      IF (CONV) LAST =.TRUE.
*
*  *****  INCREASE THE CONVERGENCE CRITERION FOR SELF-CONSISTENCY
*
1     Z2 = D2*Z2
      WRITE(OUT,8) EL(JJ),DP1
8     FORMAT(/ 6X,34HLEAST SELF-CONSISTENT FUNCTION IS ,A3,
     1   27H :WEIGHTED MAXIMUM CHANGE =,1PD10.2)
100   CFGTOL = 1.4D0*CFGTOL
18    WRITE(ERR,13)
13    FORMAT(10X/' SCF ITERATIONS HAVE CONVERGED TO THE ABOVE ACCURACY')
      WRITE(PRI,13)
      WRITE(ERR,*) ' Do you wish to continue ? (Y/N) '
      READ(IN,'(A)') ANS
      IF (ANS .EQ. 'Y' .OR. ANS .EQ. 'y') THEN
         WRITE(ERR,*) ' Enter the additional iterations and new IC '
         READ(IN,*) NSCF,IC
      CALL UPDATE
         GO TO 19
      END IF
      FAIL = .TRUE.
*
*  *****  PERFORM FINAL CALCULATIONS
*
17    ACFG = D0
      CALL DIAG(ECONV,ACFG,CFGTOL,.TRUE.)
      NIT = NWF - IB + 1
      WRITE(PRI, 105) NIT, DP1, CFGTOL
105   FORMAT(//10X,'NUMBER OF FUNCTIONS ITERATED          =',I6/
     1         10X,'MAXIMUM WEIGHTED CHANGE IN FUNCTIONS  =',D10.2/
     2         10X,'TOLERANCE FOR THE MCHF ITERATION      =',D10.2)
      RETURN
      END
*
*     ------------------------------------------------------------------
*              S E A R C H
*     ------------------------------------------------------------------
*
*       This routine searches for the NJ>70 such that YR(j) > 0 for all
*   j > NJ.
*
*
      SUBROUTINE SEARCH(NJ,I)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
      IA = 70
      IL = NO
4     IF (YR(IA) .LT. D0) GO TO 3
      IA = IA + 2
      IF (IA .LT. IL ) GO TO 4
      NJ = MAX0(70,MAX(I)-100)
      RETURN
3     NK = (IA + IL)/2
      IF (YR(NK) .LT. D0) GO TO 1
      IL = NK
      GO TO 2
1     IA = NK
2     IF (IL - IA .GT. 1) GO TO 3
      NJ = IL - 7
      RETURN
      END
*
*-----------------------------------------------------------------------
*               S E T Q L
*-----------------------------------------------------------------------
*
*       Determine if the orbitals (i,j) are in the same orthogonal set.
*
        LOGICAL FUNCTION SETEQL(I,J)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
        IBEGIN = 1
        IF (I .GT. 1) IBEGIN = IEPTR(I-1) + 1
        IEND = IEPTR(I)
        JBEGIN = 1
        IF (J .GT. 1) JBEGIN = IEPTR(J-1) + 1
        JEND = IEPTR(J)
        SETEQL = .FALSE.
        DO 10 II = IBEGIN,IEND
           DO 11 JJ = JBEGIN,JEND
              IF (IJE(II) .EQ. IJE(JJ)) GO TO 10
 11        CONTINUE
           RETURN
 10     CONTINUE
        SETEQL = .TRUE.
        END
*
*---------------------------------------------------------------------
*               S E T O R T
*---------------------------------------------------------------------
*
*       Determine if orbitals for electrons (el1, el2) should be
*   orthogonal.
*
      LOGICAL FUNCTION SETORT(EL1,EL2)
      CHARACTER*3 EL1,EL2
      CHARACTER*1 S1, S2
*
      IF (EL1(1:1) .EQ. ' ') THEN
          S1 = ' '
        ELSE
          S1 = EL1(3:3)
      END IF
      IF (EL2(1:1) .EQ. ' ') THEN
          S2 = ' '
        ELSE
          S2 = EL2(3:3)
      END IF
*
      IF (S1 .EQ. ' ' .OR. S2 .EQ. ' ') THEN
         SETORT = .TRUE.
        ELSE IF (S1 .EQ. S2) THEN
         SETORT = .TRUE.
       ELSE
         SETORT = .FALSE.
      END IF
      RETURN
      END
*
*     ------------------------------------------------------------------
*              S O L V E
*     ------------------------------------------------------------------
*
*       When FIRST is .TRUE., SOLVE computes the potential and exchange
*   function and initializes variables for the i'th radial  equation.
*   The vector P1 is the solution of the radial equation and P2 the
*   variation of the solution with respect to the energy parameter
*   E(I,I).
*
*
      SUBROUTINE SOLVE(I,FIRST)
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
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
        COMMON P2(NOD),HQ(NOD),XX(NOD),AC(30,NWD),BC(NWD),JV(NWD),
     :     AZZ,PP,FN,EM,FM,EU,FU,DELTAE,M,NODE,MK,KK,NJ
*
        LOGICAL FIRST
      DIMENSION ZERO(NOD),P1(NOD)
      EQUIVALENCE (ZERO(1),XX(1)),(PDE(1),P1(1))
*
*  *****  IF FIRST IS 'TRUE', CALL POTL AND XCH AND SET UP ARRAYS
*
      IF (.NOT. FIRST) GO TO 17
      CALL POTL(I)
      CALL XCH(I,3)
      ZINF = DMAX1(0.05D0, Z-YR(ND))
      FN = N(I)
      FL = L(I)
      V = YR(1)/R(1)
      B4 = Z*(FL+D4/D3)/((FL+D1)*(FL+D2))
      CN = (D2*Z/FN)**(L(I) +1)
      C = D4*FL +D6
      CD = (FL+D5)**2
      XY = X(1)
      XP = X(2)
      ED = E(I,I)
      X1 = X(1)
      X2 = X(2)
      X3 = X(3)
      X4 = X(4)
      DO 1 J = 3,ND
      X5 = X(J+2)
      X(J) =CH*(-X5+24.D0*(X4+X2) + 194.D0*X3 - X1)/20.D0
      X1 = X2
      X2= X3
      X3 = X4
1     X4 = X5
      X(NO-1) = CH*(X4 + D10*X3 + X2)
      DO 4 J = 1,NO
4     YK(J) = -D2*(Z - YR(J))*R(J) + CD
      X1 =    CH*P(1,I)*(YK(1)+ED*RR(1))
      X2 =    CH*P(2,I)*(YK(2)+ED*RR(2))
      X3 =    CH*P(3,I)*(YK(3)+ED*RR(3))
      X4 =    CH*P(4,I)*(YK(4)+ED*RR(4))
      DO 7 J = 3,ND
      X5 =    CH* P(J+2,I)*(YK(J+2)+ED*RR(J+2))
      X(J) = X(J) - (X5 -D4*(X2 + X4) + D6*X3 +X1)/20.D0
      X1 = X2
      X2 = X3
      X3 = X4
7     X4 = X5
      RL = L(I) + 2.5
      X(2) = R(2)**RL*(X(5)/R(5)**RL - D3*(X(4)/R(4)**RL -
     1     X(3)/R(3)**RL))
*
*  *****  DETERMINE LOWER BOUND ON THE ENERGY PARAMETER
*
      IF (KK .NE. 3) GO TO 80
      DO 11 JJ = 15,ND
      J = NO - JJ
      IF (YK(J) .LT. D0 ) GO TO 63
11    CONTINUE
      WRITE(OUT,12)
12    FORMAT(10X,'POTENTIAL FUNCTION TOO SMALL - 2R*(Z-Y)<(L+.5)**2')
*     STOP
      GO TO 80
63    EM = -YK(J)/RR(J)
      GO TO 81
80    EM = (ZINF/(FN + D5))**2
81    FM = EM
*
*  *****  DETERMINE DIAGONAL ENERGY PARAMETER
*
      F1 = D0
      C11 = D0
      M = MIN0(MAX(I),NO-1)
      DO 5 J = 2,M
      FNUM = P(J+1,I) - P(J,I) - P(J,I) + P(J-1,I)
      FNUM = FNUM - CH*(YK(J+1)*
     1   P(J+1,I) + D10*YK(J)*P(J,I) + YK(J-1)*P(J-1,I))-X(J)
      DEL1 = RR(J+1)*P(J+1,I) + D10*RR(J)*P(J,I) + RR(J-1)*P(J-1,I)
      F1 = F1 +P(J,I)*FNUM
      C11 = C11 + P(J,I)*DEL1
5     CONTINUE
      ED = F1/(C11*CH)
      IF (ED .GT. EM) GO TO 19
*
*  *****  ERROR MESSAGE AND ENERGY ADJUSTMENT FOR AN ENERGY PARAMETER
*  *****  TOO SMALL FOR THE RANGE OF THE FUNCTION
*
      WRITE(OUT,24) ED
24    FORMAT(10X,5HED = ,F10.6,36H; ADJUSTED TO ALLOWED MINIMUM ENERGY )
      ED = EM
      IF ( DABS(FM - E(I,I)) .GT. 1.D-6 .OR. KK .EQ. 3 ) GO TO 19
*
*  ***** RETURN HYDROGENIC FUNCTION
*
      PN = HNORM(N(I),L(I),ZINF)
      DO 65 J = 1,NO
65    PDE(J) = PN*HWF(N(I),L(I),ZINF,R(J))/R2(J)
      AZD = PN*(D2*ZINF/N(I))**(L(I)+1)
      PP = D0
      WRITE(OUT,66) EL(I), ZINF
66    FORMAT(//10X, 'RETURN HYDROGENIC FUNCTION FOR ',A3,
     1   ' WITH EFFECTIVE CHARGE ',F10.3)
      RETURN
*
*  *****  CHECK IF UPPER BOUND IS CORRECT
*
19    IF ( D10*ED .LT. EU) GO TO 18
      EU = D10*ED
      FU = EU
18    AZD = AZ(I)
17    DO 26 J=1,NO
      YR(J) = (YK(J) + ED*RR(J))*CH
26    ZERO(J) = D0
*
*  *****  SEARCH FOR THE POINT AT WHICH YR BECOMES POSITIVE
*
      CALL SEARCH(NJ,I)
*
*  *****  COMPUTE STARTING VALUES FROM SERIES EXPANSION
*
      B3 = (V + V + ED - (Z/FN)**2)/C
      DO 6 J = 1,2
      HW  = HWF(N(I),L(I),Z,R(J))/CN
6     HQ(J)   = AZD*(HW + R(J)**(L(I)+3)*B3*(D1-R(J)*B4))/R2(J)
*
*  *****  OBTAIN HOMOGENEOUS SOLUTION
*
      CALL NMRVS(NJ,DELH,MH,HQ,ZERO)
      P1(1) = HQ(1) + XY/C
      P1(2) = HQ(2) + XP/C
*
*  *****  OBTAIN PARTICULAR SOLUTION
*
      CALL NMRVS(NJ,DEL1,M1,P1,X)
*
*  *****  DETERMINE THE ENERGY ADJUSTMENT REQUIRED FOR A SOLUTION WITH
*  *****  GIVEN A0
*
      M = MAX0(M1,MH)
      PNORM = D0
      DO 50 J = 1,M
50    PNORM = PNORM + RR(J)*HQ(J)*P1(J)
      Y1 = P1(NJ-1)
      Y2 = P1(NJ)
      Y3 = P1(NJ+1)
      DELTA = Y2 - Y1 + Y2 - Y3 +YR(NJ-1)*Y1 + D10*YR(NJ)*Y2
     1   + YR(NJ+1)*Y3 + X(NJ)
      DELTAE = HQ(NJ)*DELTA/(H*H*PNORM)
      PP = -DEL1/DELH
*
*  *****  MATCH AT THE JOIN FOR A SOLUTION OF THE DIFFERENTIAL EQUATION
*
      DO 13 J = 1,NO
13    P1(J)   = P1(J) + PP*HQ(J)
*
*  *****  IF  THE EQUATIONS APPEAR TO BE NEARLY
*  ****  SINGULAR, SOLVE THE VARIATIONAL EQUATIONS
*
      IF (KK .NE. 2) RETURN
      X1 = P(1,I)*RR(1)
      X2 = P(2,I)*RR(2)
      P2(1) = X1/C
      P2(2) = X2/C
      DO 8 J = 3,NO
      X3 = P(J,I)*RR(J)
      XX(J-1) = (D10*X2 + X1 + X3)*CH
      X1 = X2
8     X2 = X3
      CALL NMRVS(NJ,DEL2,M2,P2,XX)
      AA = -DEL2/DELH
      M = MAX0(M,M2)
      DO 9 J = 1,NO
9     P2(J) = P2(J) + AA*HQ(J)
      A11 = QUAD(I,M,P2,P2)
      B11 = QUAD(I,M,P1,P2)
      C11 = QUAD(I,M,P1,P1) - D1
      DISC = B11*B11 - A11*C11
      IF ( DISC .LT. D0 ) GO TO 70
      DE1 = -(B11+DSQRT(DISC))/A11
      DE2 = C11/A11/DE1
      IF ( P1(3)+DE1*P2(3) .LT. D0) DE1 = DE2
      GO TO 71
70    DE1 = C11/A11
71    DO 301 J = 1,NO
      P1(J) = P1(J) + DE1*P2(J)
301   CONTINUE
      PP = PP + DE1*AA
      RETURN
      END
*
*     ------------------------------------------------------------------
*              S U M M R Y
*     ------------------------------------------------------------------
*
*       The results of a calculation are summarized.   These include
*   the following for each electron:
*
*          E(NL)   - diagonal energy parameter
*          AZ(NL)  - starting parameter, P(r)/r**(l+1) as r -> 0.
*          SIGMA   - screening parameter as defined by Eq. (6-  ).
*          1/R**3  - expected value of <1/r**3>
*          1/R     - expected value of <1/r>
*          R       - expected mean radius
*          R**2    - expected value of <r**2>
*          I(NL)   - -(1/2)<nl|L|nl>
*          KE      - I(NL) + Z <r>
*          REL     - Relativistic shift (mass-velocity, Darwin term,
*                    spin-spin contact term)
*
*   These results are followed by:
*
*          TOTAL ENERGY--RELATIVISTIC OR NON-RELATIVISTIC (ET)
*          KINETIC ENERGY-- NON-RELATIVISTIC (EN)
*          POTENTIAL ENERGY (EP) = ET - EN
*          RATIO                 = - EP/EN
*                      k   k   k
*   The values of all F , G , R  and <nl|L|n'l> integrals which enter
*   into the calculation are printed, but only if OUD > 0.
*
*
      SUBROUTINE SUMMRY
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220,NCD=100)
*
        INTEGER IN,OUT,ERR,PRI,OUC,OUD,OUF,OUH
        COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,OUC,OUD,OUF,OUH
*
        CHARACTER CONFIG*40,EL*3,ATOM*6,TERM*6,COUPLE*3
        COMMON /LABEL/CONFIG(NCD),EL(NWD),ATOM,TERM,COUPLE(NCD,9)
*
        COMMON /MATRIX/ETOTAL,W(NCD,NCD)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
        INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
        PARAMETER (IDIM=550,NCDIM=3000)
        COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :         ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
        LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
      COMMON R1(NWD),RM1(NWD),EK(NWD),SS(3)
      CHARACTER*1 SYMBOL
*
      PI = ACOS(-D1)
      WRITE (3,9) ATOM,TERM
9     FORMAT(/// 24X,'ATOM ',A6,3X,'TERM ',A6//
     :   2X,'nl',7X,'E(nl)',
     :   8X,'I(nl)',7X,'KE(nl)',8X,'Rel(nl)',3X,'S(nl)',5X,'Az(nl)')
      EN = D0
      REL = .FALSE.
*
*  *****  COMPUTE AND PRINT ONE-ELECTRON PARAMETERS
*
      DO 10 I = 1,NWF
      R1(I) = QUADR(I,I,1)
      EK(I) = -D5*HL(EL,I,I,REL)
      RM1(I) = QUADR(I,I,-1)
      EKINP = EK(I) + Z*RM1(I)
      EN = EN+ SUM(I)*EKINP
      RH = 3*N(I)*N(I) - L(I)*(L(I) + 1)
      SC = Z - D5*RH/R1(I)
      S(I) = SC
      RELS = RLSHFT(I,I)
      WRITE (3,15)EL(I),E(I,I),EK(I),EKINP,RELS,S(I),AZ(I)
15    FORMAT(1X,A3,F14.7,3F13.6,F7.2,F13.5)
10    CONTINUE
*
*  *****  Compute Moments
*
      WRITE(PRI,8) 'Delta(R)'
 8    FORMAT(//2X,'nl',6X,A8,5X,'1/R**3',7X,'1/R',9X,'R',8X,'R**2')
      DO 11 I = 1,NWF
      RM3 = 0
      IF (L(I) .NE. 0) RM3 = QUADR(I,I,-3)
      RP2 = QUADR(I,I,2)
      RZ = 0.
      IF ( L(I) .EQ. 0) RZ = AZ(I)**2/(4.*PI)
      WRITE(PRI,16) EL(I),RZ,RM3,RM1(I),R1(I),RP2
16    FORMAT(1X,A3,F14.3,F13.4,F11.5,F10.5,F11.5)
11    CONTINUE
*
*  *****  ADD CONTRIBUTION FROM THE 'L' INTEGRALS TO EN
*
      IBEGIN = INTPTR(5) + 1
      IEND = INTPTR(6)
      DO 32 I = IBEGIN,IEND
        IF (IEL(I,1) .NE. IEL(I,2)) THEN
          CONT = HL(EL,IEL(I,1),IEL(I,2),REL) 
     :           -D2*Z*QUADR(IEL(I,1),IEL(I,2),-1)
          EN = EN + CONT*COEF(I)
            END IF
32    CONTINUE
      EPOTL = ETOTAL - EN
      RATIO =-EPOTL/EN
      WRITE(OUT,26) ETOTAL,EPOTL,EN,RATIO
      WRITE(PRI,26) ETOTAL,EPOTL,EN,RATIO
26    FORMAT(//5X,'ENERGY (a.u.)'/5X,'------'/
     : 10X,' Total              ',F16.9/
     : 10X,' Potential          ',F16.9/
     : 10X,' Kinetic            ',F16.9/
     : 10X,' Ratio              ',F16.9)
*
*  *****  PRINT TABLES OF 'FK' AND 'GK' INTEGRALS WHICH WERE USED IN
*  *****  DETERMINING THE ENERGY
*
      IF ( OUD .EQ. 0 ) GO TO 13
      WRITE (OUD,126)
126   FORMAT(//2X,27HVALUES OF F AND G INTEGRALS        //)
      IBEGIN = 1
      IEND = INTPTR(2)
      DO 17 I = IBEGIN,IEND
         SYMBOL = 'F'
         IF (I .GT. INTPTR(1)) SYMBOL = 'G'
17       WRITE(OUD,19) SYMBOL,KVAL(I),EL(IEL(I,1)),EL(IEL(I,2)),VALUE(I)
19       FORMAT( 2X,A1,I2,1H(,A3,1H,,A3,4H ) =, F10.7)
*
*  *****  PRINT TABLES OF 'RK' INTEGRALS
*
      WRITE (OUD,21)
21    FORMAT(//2X,21HVALUES OF R INTEGRALS  //)
      IBEGIN = INTPTR(4) + 1
      IEND = INTPTR(5)
      DO 22 I = IBEGIN,IEND
      I1 = IEL(I,1)
      I2 = IEL(I,2)
      J1 = IEL(I,3)
      J2 = IEL(I,4)
22    WRITE (OUD,23) KVAL(I),EL(I1),EL(I2),EL(J1),EL(J2),VALUE(I)
23    FORMAT(2X,1HR,I2,1H(,2A3,1H,, 2A3,3H) =, F11.7 )
*
*  *****  PRINT TABLES OF 'L' INTEGRALS
*
      WRITE (OUD,28)
28    FORMAT(//2X,21HVALUES OF L INTEGRALS //)
      IBEGIN = IEND + 1
      IEND = INTPTR(6)
      DO 29 I = IBEGIN,IEND
29    WRITE(OUD,30) EL(IEL(I,1)),EL(IEL(I,2)),VALUE(I)
30    FORMAT(2X,2HL(,A3,1H,,A3,4H) = ,F12.7)
13    RETURN
      END
*
*-----------------------------------------------------------------------
*               U P D A T E
*-----------------------------------------------------------------------
*
*       Evaluate all integrals where at least on orbital has changed.
*
        SUBROUTINE UPDATE
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220,NCD=100)
*
        INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
        PARAMETER (IDIM=550,NCDIM=3000)
*
        CHARACTER CONFIG*40,EL*3,ATOM*6,TERM*6,COUPLE*3
        COMMON /LABEL/CONFIG(NCD),EL(NWD),ATOM,TERM,COUPLE(NCD,9)
        COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :         ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
        LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
        LOGICAL CHANGE
        IBEGIN =  1
        IEND = INTPTR(3)
        DO 1 I = IBEGIN,IEND
           IF (VARIED(IEL(I,1)) .OR. VARIED(IEL(I,2))) THEN
              IF (I .LE. INTPTR(1)) THEN
                 VALUE(I) = FK(IEL(I,1),IEL(I,2),KVAL(I),REL)
              ELSE IF (I .LE. INTPTR(2)) THEN
                 VALUE(I) = GK(IEL(I,1),IEL(I,2),KVAL(I),REL)
              ELSE
                 VALUE(I) = QUADR(IEL(I,1),IEL(I,2),0)**KVAL(I)
              END IF
           END IF
  1     CONTINUE
*
        IBEGIN = IEND + 1
        IEND = INTPTR(4)
        DO 30 I = IBEGIN,IEND
           CHANGE = .FALSE.
           DO 31 J = 1,4
             CHANGE = CHANGE .OR. VARIED(IEL(I,J))
 31        CONTINUE
           IF (CHANGE) THEN
              K1 = KVAL(I)/64
              K2 = KVAL(I) - 64*K1
              VALUE(I) = QUADR(IEL(I,1),IEL(I,2),0)**K1
     :                  *QUADR(IEL(I,3),IEL(I,4),0)**K2
           END IF
 30     CONTINUE
        IBEGIN = IEND + 1
        IEND = INTPTR(5)
        DO 10 I = IBEGIN,IEND
          CHANGE = .FALSE.
          DO 11 J = 1,4
             CHANGE = CHANGE .OR. VARIED(IEL(I,J))
 11       CONTINUE
          IF (CHANGE) VALUE(I)
     :        = RK(IEL(I,1),IEL(I,2),IEL(I,3),IEL(I,4),KVAL(I),REL)
 10     CONTINUE
*
        IBEGIN = IEND + 1
        IEND = INTPTR(6)
        DO 20 I = IBEGIN,IEND
           IF (VARIED(IEL(I,1)) .OR. VARIED(IEL(I,2))) 
     :                VALUE(I) = HLC(EL,IEL(I,1),IEL(I,2),REL)
 20     CONTINUE
*
*      ... Test if any of the core functions have changed
*
        CHANGE = .FALSE.
        DO 35 I = 1,NCLOSD
           CHANGE = CHANGE .OR. VARIED(I)
  35    CONTINUE
        IF (CHANGE .OR. EC.EQ.D0) CALL ECORE(EL,EC,REL)
*
        DO 40 I = 1,NWF
           VARIED(I) = .FALSE.
 40     CONTINUE
        END
*
*-----------------------------------------------------------------------
*               V
*-----------------------------------------------------------------------
*
*       Data structure for storing first- and second-order variations 
*   with respect to rotations.
*
        DOUBLE PRECISION FUNCTION V(I,J)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
        IBEGIN = 1
        IF (I .GT. 1) IBEGIN = IEPTR(I-1) + 1
        IEND = IEPTR(I)
        V = 0.D0
        DO 10 II = IBEGIN,IEND
           IF (IJE(II) .EQ. J) THEN
              V = VIJ(II)
              RETURN
           END IF
 10     CONTINUE
        END
*
*-----------------------------------------------------------------------
*               V I J S E T
*-----------------------------------------------------------------------
*
*       Enter the value V into the V(i,j) data structure
*
        SUBROUTINE VIJSET(I,J,V)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
        IBEGIN = 1
        IF (I .GT. 1) IBEGIN = IEPTR(I-1)+1
        IEND = IEPTR(I)
        DO 10 II = IBEGIN,IEND
           IF (IJE(II) .EQ. J) THEN
              VIJ(II) = V
              RETURN
           END IF
 10     CONTINUE
        END
*
*     ------------------------------------------------------------------
*              W A V E F N
*     ------------------------------------------------------------------
*
*       This routine initializes radial functions by the procedure
*   indicated by IND(I).
*
*         Value of IND(I)     Method
*         ---------------     ------
*             -1           Functions read from unit IU2
*              0           Screened hydrogenic functions with ZZ=Z-S(I)
*              1           Functions in memory left unchanged
*                                                  0
*   The set of functions are then orthogonalized, Y (i, i;r) and the
*   diagonal energy parameters computed, when necessary.
*
*
      SUBROUTINE WAVEFN
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
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
        LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
      COMMON ZZ(NWD),IND(NWD),PN,Z2,FN,M,K,ZT,
     1   ETI,EKI,AZI,PT(NOD),MT
*
        CHARACTER EL1*3,AT*6,TT*6,ATM(NWD)*6,TRM(NWD)*6,TITLE*24
*
*  *****  GENERATE ARRAYS FOR R,R*R AND SQRT(R) WITH A CONSTANT MESH
*  *****  SIZE IN THE LOG(Z*R) VARIABLE
*
      DO 1 I=1,NO
      R(I)= DEXP(RHO)/Z
      RR(I) = R(I)*R(I)
      R2(I) = DSQRT(R(I))
1     RHO = RHO + H
      RHO = RHO - NO*H
*
*  ***** READ THE WAVEFUNCTIONS
*
      IF (IUF .EQ. 0) GO TO 5
2     READ(IUF,END=5) AT,TT,EL1,MM,ZT,ETI,EKI,AZI,(PT(J),J=1,MM)
      M = MIN0(NO,MM)
      CALL EPTR(EL,EL1,I,*2)
      IF ( I .GT. 0 .AND. IND(I) .EQ. -1) THEN
         ATM(I) = AT
         TRM(I) = TT
         MAX(I) = M
         ZZ(I)  = ZT
         C = D1
         IF ( Z .NE. ZT ) C = Z/ZT
*
*  *****  SCALE RESULTS IF CARDS ARE FOR AN ATOM WITH A DIFFERENT Z
*
         CALL EIJSET(I,I,C*C*ETI)
         AZ(I)  = AZI*C**(L(I)+1)*DSQRT(C)
         DO 11 J = 1,M
            P(J,I) = C*PT(J)
11       CONTINUE
*
*  *****  SET REMAINING VALUES IN THE RANGE = 0.
*
         IF ( M .EQ. NO ) GO TO 12
         M = M +1
         DO 13  J=M,NO
13       P(J,I) = D0
12       IND(I) = -2
      ENDIF
      GO TO 2
*
*  *****  SET PARAMTERS FOR ELECTRONS AND INITIALIZE FUNCTIONS
*
5     DO 9 I = 1,NWF
      IF (IND(I)) 7,8,9
*
*  ***** WAVE FUNCTIONS NOT FOUND IN THE INPUT DATA, SET IND = 0
*
7     IF ( IND(I) .EQ. -2 ) GO TO 4
      IF (METH(I) .EQ. 4) GO TO 9
      IND(I) = 0
      WRITE(OUT,27) EL(I)
27    FORMAT(8X,'WAVE FUNCTIONS NOT FOUND FOR ',A3)
*
*  *****  DETERMINE ESTIMATES OF THE WAVE FUNCTIONS BY THE SCREENED
*  *****  HYDROGENIC APPROXIMATION
*
8     PN = HNORM(N(I),L(I),Z-S(I))
      DO 3 J=1,NO
      P(J,I) = PN*HWF(N(I),L(I),Z-S(I),R(J))/R2(J)
3     CONTINUE
      M = NO
30    IF ( DABS(P(M,I)) .GT. 1.D-15 ) GO TO 31
      P(M,I) = D0
      M = M-1
      GO TO 30
31    MAX(I) = M
*
*  ***** SET THE AZ(I) VALUE
*
      AZ(I) = PN*(D2*(Z - D5*S(I))/N(I))**(L(I) + 1)
      CALL EIJSET(I,I,D0)
*
*  *****  ORTHOGONALIZE TO INNER FUNCTIONS
*
4      IF (I .EQ. 1 ) GO TO 9
      IM = I - 1
      DO 6 II =1,IM
      IF (E(I,II) .EQ. D0) GO TO 6
      PN = QUADR(I,II,0)
      IF ( DABS(PN) .GT. 1.D-8 ) THEN
         PNN = DSQRT(D1 - PN*PN)
         IF (P(50,I) - PN*P(50,II) .LT. D0) PNN = -PNN
         M = MAX0(MAX(I),MAX(II))
         DO 25 J = 1,M
25          P(J,I) =(P(J,I) - PN*P(J,II))/PNN
      END IF
6     CONTINUE
9     CONTINUE
      WRITE(PRI,14)
14    FORMAT(/// 8X,18HINITIAL ESTIMATES  //10X,2HNL,
     1   4X,5HSIGMA,6X,5HE(NL),4X,6HAZ(NL),4X,9HFUNCTIONS//)
*
*  *****  COMPUTE ONE-ELECTRON ENERGY PARAMETERS IF THEY WERE NOT
*  *****  SPECIFIED ON INPUT.
*
      DO 15 I = 1,NWF
*     IF (E(I,I) .EQ. D0) E(I,I) = HL(EL,I,I,REL) - EKIN(I,I)
      K = IND(I) + 2
      IF ( IND(I) .EQ. -2 ) THEN
           TITLE = ' SCALED '//ATM(I)//TRM(I)
        ELSE IF (IND(I) .EQ. 0) THEN
           TITLE = ' SCREENED HYDROGENIC'
        ELSE
           TITLE = ' UNCHANGED'
      END IF
17    WRITE(PRI,19) EL(I),S(I),E(I,I),AZ(I),TITLE
19    FORMAT(9X,A3,F9.2,F11.3,F10.3,3X,A24)
15    CONTINUE
      IF ( IUF .NE. 0) REWIND(UNIT=IUF)
      RETURN
      END
*
*----------------------------------------------------------------------
*               X C H
*----------------------------------------------------------------------
*
*       Compute the function X(r) that includes contributions to the 
*   equations for orbital(i) arising from Gk (exchange), Rk (electro-
*   static interactions), L (one-electron part of hamiltonian), terms 
*   in the energy expression, and optionally also contributions from 
*   the off-diagonal energy parameters. The exact form depends on the 
*   value of IOPT.
*
*          
*          Value of IOPT      Function
*          -------------      --------
*              1             SQRT(r) X(r)
*              2             X(r)/SQRT(r)
*              3             r SQRT(r) ( X(r) + SUM e   P )
*                                                    ij  j
*          
        SUBROUTINE XCH(I,IOPT)
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
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
        INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
        PARAMETER (IDIM=550,NCDIM=3000)
        COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :         ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
        LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
        LOGICAL SAME,EXIT
*
        DO 1 J=1,NO
  1     X(J) = D0
        DO 2 J = 1,NWF
           IF ((I.LE.NCLOSD .AND. I.NE.J) .OR.
     :         (I.GT.NCLOSD .AND. J.LE.NCLOSD))  THEN
              DO 4 K = IABS(L(I)-L(J)),L(I)+L(J),2
                 C = - D2*CB(L(I),L(J),K)*SUM(J)
                 CALL YKF(J,I,K,REL)
                 DO 6 JJ = 1,NO
                    X(JJ) = X(JJ) + C*YK(JJ)*P(JJ,J)
  6              CONTINUE
  4           CONTINUE
           END IF
  2     CONTINUE
        SUMI = SUM(I)
        IF (I .LE. NCLOSD) GO TO 51
*
        IBEGIN = INTPTR(1)+1
        IEND = INTPTR(2)
        DO 7 INT = IBEGIN,IEND
           IE1 = 0
           IF (IEL(INT,1) .EQ. I) THEN
              IE1 = IEL(INT,1)
              IE2 = IEL(INT,2)
           ELSE IF (IEL(INT,2) .EQ. I) THEN
              IE1 = IEL(INT,2)
              IE2 = IEL(INT,1)
           END IF
           IF (IE1 .NE. 0) THEN
              C = D2*COEF(INT)/SUMI
              CALL YKF(IE1,IE2,KVAL(INT),REL)
              DO 8 JJ = 1,NO
                 X(JJ) = X(JJ) + C*YK(JJ)*P(JJ,IE2)
 8            CONTINUE
           END IF
 7      CONTINUE
*
        IBEGIN = INTPTR(4) + 1
        IEND = INTPTR(5)
        DO 50 INT = IBEGIN,IEND
           I1 = IEL(INT,1)
           I2 = IEL(INT,2)
           J1 = IEL(INT,3)
           J2 = IEL(INT,4)
           KK = KVAL(INT)
           IF ((I1-I)*(I2-I) .EQ. 0 .OR. (J1-I)*(J2-I) .EQ. 0) THEN
              C = COEF(INT)/SUMI
              CC = C
*
*  ***** COUNT THE NUMBER OF OCCURRENCES OF I
*
              IK = 0
              IF (I1 .EQ. I) IK = IK + 1
              IF (I2 .EQ. I) IK = IK + 1
              IF (J1 .EQ. I) IK = IK + 1
              IF (J2 .EQ. I) IK = IK + 1
              EXIT = .FALSE.
              DO 11 II2=1,2
              DO 12 II1=1,2
              GO TO (10, 20, 30, 40) IK
10          CONTINUE
*
*  ***** I OCCURS JUST ONCE IN RK
*
              IF (I1 .NE. I) GO TO 13
              GO TO 16
20            CONTINUE
*
*  ***** I OCCURS TWICE IN THE RK INTEGRAL
*
              IF (I1 .NE. I) GO TO 13
              IF (J1 .EQ. I) GO TO 17
*
*  ***** TEST IF THE PAIR (I1,J1) = PAIR (I2,J2)
*
              ICODE1 = 100*I1 + J1
              ICODE2 = 100*I2 + J2
              ICODE3 = 100*J2 + I2
              SAME = ICODE1 .EQ. ICODE2 .OR. ICODE1 .EQ. ICODE3
              IF ( .NOT. SAME ) GO TO 15
              GO TO 17
30            CONTINUE
*
*  ***** I OCCURS THREE TIMES IN THE RK INTEGRAL
*
*
              IF (I1 .EQ. I) GO TO 13
              CALL YKF(I2, J2, KK, REL)
              DO 33 J = 1,NO
33              X(J) = X(J) + CC*P(J,I1)*YK(J)
              CALL YKF(I1, J1, KK, REL)
              CC = D2*CC
              DO 34 J = 1,NO
34               X(J) = X(J) + CC*P(J,I2)*YK(J)
              GO TO 50
*
*  ***** I OCCURS FOUR TIMES IN RK INTEGRAL
*
40            CC = D4*CC
              GO TO 16
17            CC = D2*CC
16            EXIT = .TRUE.
15            CALL YKF(I2,J2,KK,REL)
              DO 14 J=1,NO
14               X(J) = X(J) +CC*P(J,J1)*YK(J)
              IF (EXIT) GO TO 50
13            III = I1
              I1= I2
              I2= III
              III = J1
              J1 = J2
12            J2 = III
              III = I1
              I1 = J1
              J1 = III
              III = I2
              I2= J2
11            J2= III
           END IF
50      CONTINUE
*
51      IBEGIN = INTPTR(5) + 1
        IEND = INTPTR(6)
        DO 60 INT = IBEGIN,IEND
*         ... Include only if off-diagonal ...
          IF (IEL(INT,1).NE.IEL(INT,2)) THEN
           I1 = IEL(INT,1)
           I2 = IEL(INT,2)
           IF (I1 .NE. I) THEN
              ITEMP = I1
              I1 = I2
              I2 = ITEMP
           END IF
           IF (I1 .EQ. I) THEN
              C = COEF(INT)/SUMI
              CALL DIFF(I2)
              DO 62 J = 1,NO
                 X(J) = X(J) + C*YK(J)/R(J)
 62           CONTINUE
              DO 64 II = 1,NCLOSD
                 CC = -D2*(4*L(II)+2)*C
                 CALL YKF(II,II,0,REL)
                 DO 65 J = 1,NO
                    X(J) = X(J) + CC*YK(J)*P(J,I2)
 65              CONTINUE
                 DO 66 K = IABS(L(I)-L(II)),L(I)+L(II),2
                    CCC = CC*CB(L(I),L(II),K)
                    CALL YKF(I2,II,K,REL)
                    DO 67 J = 1,NO
                       X(J) = X(J) - CCC*YK(J)*P(J,II)
 67                 CONTINUE
 66              CONTINUE
 64           CONTINUE
           END IF
           IF (I .LE. NCLOSD) THEN
              C = -D2*COEF(INT)
              CALL YKF(I1,I2,0,REL)
              CC = D2*C
              DO 61 J = 1,NO
                X(J) = X(J) + CC*YK(J)*P(J,I)
 61           CONTINUE
              DO 63 K = IABS(L(I)-L(I1)),L(I)+L(I1),2
                CC = C*CB(L(I),L(I1),K)
                CALL YKF(I2,I,K,REL)
                DO 68 J = 1,NO
                   X(J) = X(J) - CC*YK(J)*P(J,I1)
 68             CONTINUE
                CALL YKF(I1,I,K,REL)
                DO 69 J = 1,NO
                   X(J) = X(J) - CC*YK(J)*P(J,I2)
 69             CONTINUE
 63           CONTINUE
           END IF
          END IF
 60     CONTINUE
        IF (I .LE. NCLOSD) GO TO 71
*
        IBEGIN = INTPTR(2) + 1
        IEND = INTPTR(3)
        DO 70 INT = IBEGIN,IEND
           I1 = IEL(INT,1)
           I2 = IEL(INT,2)
           K1 = KVAL(INT)
           IF (I1 .NE. I) THEN
              ITEMP = I1
              I1 = I2
              I2 = ITEMP
           END IF
           IF (I1 .EQ. I) THEN
              C = COV(INT)/SUMI
              IF (K1 .GT. 1) C = C*K1*QUADR(I1,I2,0)**(K1-1)
              DO 72 J = 1,NO
                 X(J) = X(J) + C*P(J,I2)*R(J)
 72           CONTINUE
           END IF
 70     CONTINUE
*
        IBEGIN = IEND + 1
        IEND = INTPTR(4)
        DO 80 INT = IBEGIN,IEND
           I1 = IEL(INT,1)
           I2 = IEL(INT,2)
           I3 = IEL(INT,3)
           I4 = IEL(INT,4)
           K1 = KVAL(INT)/64
           K2 = KVAL(INT) - 64*K1
           OV1 = D0
           OV2 = D0
           DO 82 II = 1,2
              IF (I1 .NE. I) THEN
                 ITEMP = I1
                 I1 = I2
                 I2 = ITEMP
              END IF
              IF (I1 .EQ. I) THEN
                 IF (OV2 .EQ. D0) OV2 = QUADR(I3,I4,0)
                 C = OV2**K2*COV(INT)/SUMI
                 IF (OV1 .EQ. D0 .AND. K1 .GT. 1)
     :              OV1 = QUADR(I1,I2,0)
                 IF (K1 .GT. 1) C = K1*C*OV1**(K1-1)
                 DO 84 J = 1,NO
                    X(J) = X(J) + C*P(J,I2)*R(J)
 84              CONTINUE
              END IF
              ITEMP = I1
              I1 = I3
              I3 = ITEMP
              ITEMP = I2
              I2 = I4
              I4 = ITEMP
              ITEMP = K1
              K1 = K2
              K2 = ITEMP
              OTEMP = OV1
              OV1 = OV2
              OV2 = OTEMP
 82        CONTINUE
 80     CONTINUE
 71     GO TO (75,76,77),IOPT
 76     DO 78 J = 1,NO
 78        X(J) = X(J)/R(J)
        GO TO 75
 77     DO 79 J =1,NO
 79        X(J) = R(J)*X(J)
        DO 74 J = 1,NWF
           IF (J .NE. I) THEN
           C = E(I,J)
           IF (DABS(C) .LE. 1.D-20 ) GO TO 74
           DO 73 JJ = 1,NO
 73        X(JJ) = X(JJ) + C*P(JJ,J)*RR(JJ)
           END IF
 74     CONTINUE
*
*  *****  CHECK IF EXCHANGE IS ZERO: IF SO, METHOD 2 SHOULD BE USED.
*
 75     IF (METH(I) .EQ. 2 .OR. METH(I) .GT. 3) RETURN
        IF ( DABS(X(1)) + DABS(X(2)) + DABS(X(3)) .EQ. D0 ) METH(I) = 2
        END
