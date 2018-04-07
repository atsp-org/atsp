*     ------------------------------------------------------------------
*     		AUTOIONIZATION PROGRAM
*
*	Written by :	Charlotte Froese Fischer
*			and
*			Tomas Brage
*			Department of Computer Science
*			Vanderbilt University
*       Modified by:    Jinhua Xi, Vanderbilt University
*                       December, 1994
*     ------------------------------------------------------------------
*
*
*     All comments in the program listing assume the radial function P
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
*
*
*     ------------------------------------------------------------------
*              M A I N   P R O G R A M
*     ------------------------------------------------------------------
*	The main program: 
*		- sets unit numbers and opens files.
*		- prompts for type of calculation.
*		- initializes variabels. (by calling ADATA).
*		- determines data about problem.
*		- performs calculations. (by calling CALCAUTO).
*		
*
      PROGRAM AUTO
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500)
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
      LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
      COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
      COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :     ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
      COMMON /INPUT/ ECB,ICASE
*
      CHARACTER NAME(6)*24
      LOGICAL PRINT
CSUN  REAL TIMES(2),DTIME
*
CSUN  RTIME = DTIME(TIMES)
*
*  ***** Define unit numbers and open files *********************
*                                                               *
*      UNIT NUMBERS AND FILE NAMES MAY BE MACHINE               *
*      DEPENDENT. CHECK THE FOLLOWING SECTION.                  *
*                                                               *
*      IN - Standard input unit, normally the terminal          *
*      OUT- Standard output unit, normally the terminal         *
*      ERR- Prompts and Error messages, always the terminal     *
*      PRI- Printer output unit or file.                        *
*       10- Data file (auto.dat).				*
*                                                               *
      IN = 5
      OUT = 6
      ERR = 0
      PRI = 3
*								*
*      IUC - Input unit for configuration list.			*
*      IUD - Input unit for integral list.			*
*      IUF - Input unit for wave functions.			*
*     
*      OUC - Output unit for configuration list. (Not used). 	*
*      OUD - Output unit for integral list. (Not used).		*
*      OUF - Output unit for wave functions.			*
*     
      IUC = 21
      IUD = 22
      IUF = 23
      OUC = 0
      OUD = 0
      OUF = 26
*
*  *****  Write out header
*
      WRITE(OUT,9)
9     FORMAT(//10X,'============================================'/
     :         10X,'        A U T O _ I O N I Z A T I O N       '/
     :         10X,'============================================')
*
*  *****  Write out dimension information
*
      WRITE(OUT,99) 'NCD',NCD,'NWD',NWD,'NO',NOD
99    FORMAT(//10X,'THE DIMENSIONS FOR THE CURRENT VERSION ARE:'/
     :       (10X,3(A6,'=',I3,4X)/)/)
*
*  *****  Initialize common data arrays
*
      CALL INITA
      CALL INITR
*
1     WRITE(ERR,'(//A/A//)') ' START OF CASE',' ============='
*  
*  ***** Get parameters for the case
*
999   NAME(2) = 'int.lst'
      NAME(6) = 'wfn.out'
      write(ERR,*) ' Name of State'
      read(IN,'(A)') NAME(1)
      j = index(NAME(1),' ')
      if (j .eq. 1) then
         WRITE(ERR,*) ' Names may not start with a blank'
         GO TO 999
      else
         NAME(1) = NAME(1)(1:j-1)//'.c'
         NAME(3) = NAME(1)(1:j-1)//'.w'
         NAME(4) = NAME(1)(1:j-1)//'.l'
         NAME(5) = NAME(1)(1:j-1)//'.j'
      end if
*
      OPEN(UNIT=IUC,FILE=NAME(1),STATUS='OLD')
      OPEN(UNIT=IUD,FILE=NAME(2),STATUS='OLD')
      OPEN(UNIT=IUF,FILE=NAME(3),STATUS='OLD',
     :              FORM='UNFORMATTED')
      OPEN(UNIT=OUF,FILE=NAME(6),STATUS='UNKNOWN',
     :              FORM='UNFORMATTED')
      OPEN(UNIT=3, FILE='auto.log',STATUS='UNKNOWN',FORM='formatted')
*
*  ************************************************************	
*       The input for the discrete state can be made in four
*       different ways:
*       ICASE	weights from	ECORE from	EBOUND from
*	====================================================
*         1	 .c file	 terminal	 .c file
*	  2	 .l file	 .c file	 .l file
*	  3	 .j file	 .c file	 .j file
*  ************************************************************
*
      WRITE(ERR,'(//A/A/A/A)')
     :'Select input file for discrete state(s):',
     :'   1  name.c', '   2  name.l ','   3  name.j'
      READ(IN,*) ICASE
      IF(ICASE.EQ.2) THEN
        OPEN(UNIT=9,FILE=NAME(4),STATUS='OLD')
      ELSE IF(ICASE.EQ.3) THEN
        OPEN(UNIT=9,FILE=NAME(5),STATUS='OLD')
      END IF
      OPEN(UNIT=10,FILE='auto.dat',STATUS='UNKNOWN')
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
*  *****  Determine data about the problem
*
      CALL ADATA(ECORE)
*
*
*  *****  Set parameters to their default value
*
13    PRINT = .FALSE.
      SCFTOL = 1.D-6
      NSCF = 20
      IC = 0
      ACFG = D0
      TRACE = .FALSE.
      WRITE(OUT,2) PRINT,CFGTOL,SCFTOL,NSCF,IC,ACFG,ID,TRACE
2     FORMAT(/L3,2D6.1,2I3,F3.1,I3,L3)
*
* ***** Start actual calculations.
*
      CALL CALCAUTO(ECORE,ACFG,SCFTOL,PRINT)
*
*  *****  Determine end of case
*
6     CONTINUE
CSUN  RTIME = DTIME(TIMES)
CSUN  WRITE(ERR,'(//A/A//A/3F10.3//)') ' END OF CASE',' ===========',
CSUN :  '    Real      User      System  Time (in minutes)',
CSUN : RTIME/60.,TIMES(1)/60., TIMES(2)/60.
      END
*
*     ------------------------------------------------------------------
*           A D A T A
*     ------------------------------------------------------------------
*
*	ADATA performs the following tasks:
*	
*		1) Prompts user for ATOM, TERM, Z.
*		2) Reads the cfile.
*		3) interprets the cfile.
*		4) determine continuum data.
*		5) interprets configuration data.
*		6) sets orthogonality constraints.
*		7) first print outs.
*		8) Initialize arrays. (wavefunctions and integrals)
*	
*	The following routines are called during different stages:
*
*		5) EPTR.
*		6) EIJSET and EPTR.
*		7) WAVEFN and ANTGRL.
*		
*
      SUBROUTINE ADATA(ECORE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500,IDIM=1000,NCDIM=32768)
*
      INTEGER IN,OUT,ERR,PRI,OUC,OUD,OUF,OUH
      COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,OUC,OUD,OUF,OUH
*
      CHARACTER CONFIG*40,EL*3,ATOM*6,TERM*6,COUPLE*3, ans
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
     :       ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
      LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
      COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
      COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :     ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
      COMMON ZZ(NWD),IND(NWD),IELI(5),NOCCSH(NCD)
*
      COMMON /INPUT/ ECB,ICASE
*
      LOGICAL SETORT,STRONG
      CHARACTER*3 EL1,EL2,ELCLSD(18),ELORT(10,2),ELI(5),STRING*40
*
    1 FORMAT(18(1X,A3))
    7 FORMAT(A3,F6.0,I3,I3,F3.1)
*
*  *****  Read 'ATOM' card
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
*  ***** Read configuration cards and normalize the weights
*  *****  Read configurations from a file
*
      READ(IUC,'(15X,F14.7/18(1X,A3))') ECB,(ELCLSD(I),I=1,18)
      NCFG = 0
    3 READ(IUC,'(A40,F10.8)',END=10) STRING,W
      
      IF (STRING(1:1) .NE. '*' .AND. STRING(1:3) .NE. '   ') THEN
        NCFG = NCFG+1
        IF (NCFG .LE. (NCD) ) THEN
          CONFIG(NCFG) = STRING
          WT(NCFG) = W
          READ(IUC,'(9(5X,A3))') (COUPLE(NCFG,J),J=1,9)
          GO TO 3
        ELSE
          WRITE(ERR,*) ' TOO MANY CONFIGURATIONS: MAX =', NCD
          STOP
        END IF
      END IF
10    CONTINUE
      ID = -1
*
*  *****  Determine NCLOSD shells
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
        STOP ' TOO MANY CLOSED SHELLS: MAX = 18'
      END IF
      NCLOSD = I
*
*  *****  Determine the other electrons
*
      MAXORB = NCLOSD
      DO 15 NC = 1,NCFG
        STRING = CONFIG(NC)
        J = 2
        I = 0
 16     IF (STRING(J:J+2) .NE. '   ' ) THEN
*
*  *****     An electron has been found; is it a new one?
*
          I = I+1
          EL1 = STRING(J:J+2)
          K = NCLOSD + 1
 17       IF (K .LE. MAXORB) THEN
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
*  *****  A new electron has been found; add it to the list
*
            MAXORB = K
            EL(MAXORB) = EL1
            IF ((EL1(1:1) .EQ. 'k' .OR. EL1(2:2) .EQ. 'k' .OR.
     :          EL1(1:1) .EQ. 'n' .OR. EL1(2:2) .EQ. 'n') .AND.
     :          ID .EQ. -1) ID = NC-1
          END IF
          J = J+8
          IF (J .LT. 40) GO TO 16
        END IF
        NOCCSH(NC) = I
   15 CONTINUE
      IF (ID .EQ. -1) THEN
        WRITE(ERR,*) ' STOP in ADATA: No continuum function found'
        CALL EXIT(1)
      ELSE IF (ID.EQ.0) THEN
        ICASE = 4
      END IF
      IF (ICASE.EQ.1) THEN
        WRITE(ERR,*) ' Give energy of core state (a.u.):'
        READ(IN,*) ECORE
      ELSE IF (ICASE.EQ.2 .OR. ICASE.EQ.3) THEN
        ECORE = ECB
      ELSE IF (ICASE.EQ.4) THEN
        WRITE(ERR,*) 'Give energy of continuum electron (a.u.):'
        READ(IN,*) ECB
      END IF
*
*  *****  The list of electrons has been determined
*
      NWF = MAXORB
      WRITE(ERR,19) MAXORB,(EL(J),J=1,MAXORB)
   19 FORMAT(/' There are ',I3,' orbitals as follows:'/(1X,18(1X,A3)))
*
      write(err,*) ' Is continuum function to be computed? (y/n)'
      read (in,'(A)') ans
      if (ans .eq. 'y  ' .or. ans .eq. 'Y  ') then
        NIT = 1
      else 
        NIT = 0
      end if
      IB = NWF - NIT + 1
      DO 20 I = NCLOSD+1,NWF
        S(I) = SS
        METH(I) = 3
        ACC(I) = D0
        IND(I) =-1 
        VARIED(I) = .TRUE.
        J = 2
        IF (EL(I)(1:1) .EQ. ' ') J = 3
        L(I) = LVAL(EL(I)(J:J))
        N(I) = ICHAR(EL(I)(J-1:J-1)) - ICHAR('1') + 1
        IF (EL(I)(J-1:J-1) .EQ. 'n') METH(I) = 5
 20   CONTINUE
*
*  *****  Check that last electron is a continuum orbital.
*
      EL1 = EL(NWF)
      IF (.NOT. (EL1(1:1).EQ.'k' .OR. EL1(2:2).EQ.'k'))
     :     STOP ' Last orbital not a continuum orbital'
*
*  *****  Check method and initialize arrays for the continuum 
*  *****  electron.
*
      IF (METH(NWF) .NE. 4) THEN
        METH(NWF) = 4
*       IND(NWF) = 1
        DO 95 J = 1,NO
          P(J,NWF) = D0
95      CONTINUE
        AZ(NWF) = D1
        MAX(NWF) = NO
      END IF
      DO 35 NC = 1,NCFG
        STRING = CONFIG(NC)
        J = 2
        I = 0
 30     IF (STRING(J:J+2) .NE. '   ' ) THEN
*
*  *****     An electron has been found; find its index
*
          I = I+1
          ELI(I) = STRING(J:J+2)
          CALL EPTR(EL,ELI(I),IELI(I),*99)
          READ(STRING(J+4:J+5),'(I2)') IQ
          J = J+8
          IF (J .LT. 40) GO TO 30
        END IF
*
*  *****  Define all orbitals in the configuration to be orthogonal
*
        DO 34 I1 = 2,I
          J1 = IELI(I1)
          DO 33 I2 = 1,I1-1
            J2 = IELI(I2)
            IF (L(J1) .EQ. L(J2) ) THEN
              CALL EIJSET(J1,J2,1.D-5)
              CALL EIJSET(J2,J1,1.D-5)
            END IF
   33     CONTINUE
   34   CONTINUE
   35 CONTINUE
*
*  ***** Set the following orbitals orthogonal
*
*        1) Orbitals with different l's
*        2) In the same orthogonal set
*        3) Specified orthogonality
*
*  *****
*
      DO 38 I = 2,NWF
        DO 39 J = 1,I-1
          IF (L(I) .EQ. L(J) .AND. SETORT(EL(I),EL(J)) ) THEN
            C = 1.D-5
            IF (I.LE.NCLOSD .AND. J.LE.NCLOSD) C = 1.D-10
            CALL EIJSET(I,J,C)
            CALL EIJSET(J,I,C)
          END IF
   39   CONTINUE
   38 CONTINUE
*
*  *****  Determine additional orthogonality pairs
*
      I = 0
      IF ( IUC .NE. IN) THEN
   40   READ(IUC,1,END=50) EL1,EL2
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
*  *****  Additional parameters
*
      NO = (NOD)
      REL = .FALSE.
      STRONG = .FALSE.
      IF (NCFG .GT. 1) STRONG = .TRUE.
      ND = NO - 2
      WRITE(OUT,61) ATOM,TERM,Z,NO,NWF,NIT,NCFG,REL,STRONG
61    FORMAT(/1X,2A6,F6.0,I6,3I3,2L3)
      WRITE(PRI,62) ATOM,TERM,Z,(EL(I),4*L(I)+2,I=1,NCLOSD)
62    FORMAT(1H1///9X,33HHARTREE-FOCK WAVE FUNCTIONS FOR  ,2A6,4H Z =,
     1   F5.1//14X,'CORE = ',5(A3,'(',I2,')'))
      WRITE(PRI,'(//11X,A,37X,A//)') 'CONFIGURATION','WEIGHT'
      OMIT = .NOT. STRONG
*
* *****  Write 'CONFIGURATION' cards  and check the weights
*
      DO 68 I = 1,NCFG
        NOCC=NOCCSH(I)
        WRITE(PRI,70) I, CONFIG(I), WT(I),(COUPLE(I,J),J=1,NOCC)
70      FORMAT(/3X,I3,6X,A40,F19.8/12X,9(5X,A3))
        WRITE(PRI,73) (COUPLE(I,J),J=NOCC+1,2*NOCC-1)
73      FORMAT(23X,4(5X,A3))
68    CONTINUE
      WRITE(PRI,71)
71    FORMAT(//9X,10HINPUT DATA/9X,10H----- ----//13X,13HWAVE FUNCTION,
     1   11H  PROCEDURE/17X,22HNL  SIGMA METH ACC OPT///)
      DO 79 I = 1,NWF
        WRITE(PRI,78) I,EL(I),N(I),L(I),S(I),METH(I),ACC(I),IND(I)
78      FORMAT(I8, 2X,A3,2I3,F7.1,I4,F4.1,I4)
79    CONTINUE
*
*  *****  Initialize arrays, if necessary
*
      CALL WAVEFN
      DO 100 I=1,6
        INTPTR(I) = 0
100   CONTINUE
*
      IF (IUD .NE. IN) CALL ANTGRL
*
*  *****  Define SUM(I)
*
      IBEGIN = INTPTR(5)+1
      IEND = INTPTR(6)
      DO 80 I = IBEGIN,IEND
        IF (IEL(I,1).EQ.IEL(I,2)) SUM(IEL(I,1)) = -2*COEF(I)
 80   CONTINUE
      RETURN
99    STOP
      END
*
*     ------------------------------------------------------------------
*          A G R A N G E
*     ------------------------------------------------------------------
*
*       Controls the calculation of off-diagonal energy parameters.
*   It searches for all localized orbitals, i, to which the continuum
*   orbital, P   , is constrained to be orthogonal. If P    must be
*             NWF                                       NWF
*   orthogonal not only to P  but also to P , a system of equations
*                           i              j
*   must be solved, the exact form depending on whether or not any of 
*   the functions are part of the frozen  core.
*
      SUBROUTINE AGRANGE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500,IDIM=1000,NCDIM=32768)
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
     :     ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
      COMMON W(NCD,NCD),U(NCD),DWT(NCD),AC(30,NWD),BC(NWD),
     :     JV(NWD),IV(NWD)
      LOGICAL DIAG, FIRST
*
*
*  *****  For l=L(NWF) compute off-diagonal energy parameters
*
      IL = L(NWF)
      IJ = 0
      DO 12 J = 1,NWF-1
        IF (DABS(E(NWF,J)) .GT. 1.D-10 ) THEN
        IJ = IJ + 1
        IF ( IJ .GT. NWD) STOP '  TOO MANY LAGRANGE MULTIPLIERS'
          IV(IJ) = NWF
          JV(IJ) = J
        END IF
12    CONTINUE
*
*  ***** IJ is the number of lagrange multipliers for l = L(NWF)
*
      IF (IJ .EQ. 0) GO TO 10
      DO 13 II = 1,IJ
        BC(II) = D0
        DO 14 III = 1,IJ
          AC(II,III) = D0
14      CONTINUE
13    CONTINUE
      FIRST = .TRUE.
      DO 18 II = 1,IJ
        J = 0
        IF ( IV(II) .EQ. NWF) THEN
          J = JV(II)
          IF (FIRST) THEN
            CALL CXCH(NWF,2)
            CALL POTL(NWF)
            DO 20 JJ = 1,NO
              YK(JJ) = YR(JJ)
20          CONTINUE
            FIRST = .FALSE.
          END IF
          DO 22 JJ = 1,NO
            YR(JJ) = P(JJ,J)
22        CONTINUE
          BC(II) = BC(II) +
     :    HL(EL,NWF,J,REL)-D2*QUADS(NWF,J,1)-QUAD(J,NO,YR,X)
        END IF
18    CONTINUE
      DO 24 II = 1,IJ
        DO 26 III = 1,II
          IF ( II .EQ. III) THEN
            AC(II,II) = D1/SUM(IV(II))
          ELSE IF (IV(II) .EQ. IV(III) .AND.
     :               E(JV(II),JV(III)) .EQ. D0 ) THEN
            AC(II,III) = QUADR(JV(II),JV(III),0)/SUM(IV(II))
            AC(III,II) = AC(II,III)
            DIAG = .FALSE.
          END IF
26      CONTINUE
24    CONTINUE
      IF ( .NOT. DIAG ) CALL LINEQN(30,IJ,AC,BC)
      DO 28 II = 1,IJ
        CALL EIJSET(IV(II),JV(II),BC(II)/SUM(IV(II)))
28    CONTINUE
10    CONTINUE
*
*  *****  Print the off-diagonal energy parameters
*
      DO 32 J = 1,NWF-1
        IF (DABS(E(NWF,J)) .GT. 1.D-10) THEN
          WRITE(OUT,35) EL(NWF),EL(J),E(NWF,J),EL(J),EL(NWF),E(J,NWF)
35        FORMAT(7X,2(3X,'E(',2A3,') =',1PE12.5))
        END IF
32    CONTINUE
      RETURN
      END
*
*-----------------------------------------------------------------------
*		A N T G R L
*-----------------------------------------------------------------------
*
*    Read the integrals that define the energy expression
*
      SUBROUTINE ANTGRL
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220,NCD=500)
*
        INTEGER IN,OUT,ERR,PRI,OUC,OUD,OUF,OUH
        COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,OUC,OUD,OUF,OUH
*
        CHARACTER CONFIG*40,EL*3,ATOM*6,TERM*6,COUPLE*3
        character*3 eltemp(2)
        COMMON /LABEL/CONFIG(NCD),EL(NWD),ATOM,TERM,COUPLE(NCD,9)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
        INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
        PARAMETER (IDIM=1000,NCDIM=32768)
        COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :       ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
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
             CPTR(I) = ICPTR + LAST

             CALL EPTR(EL, EL1,IEL(I,1),*999)
             CALL EPTR(EL, EL2,IEL(I,2),*999)
             I = I + 1
             IF (I .LE. (IDIM) ) GO TO 12
             WRITE(ERR,*) ' Too many integrals - MAX =',IDIM
             STOP
           ELSE
   14        IF (INT.EQ.5) THEN
*
*	        ... R integrals ...
*
               READ(IUD,2) END, KVAL(I), EL1, EL2, EL3, EL4, ICPTR
               if((el1(1:1).eq.' ' .and. el1(2:2).eq.'k') .or.
     :           (el3(1:1).eq.' '. and .el3(2:2).eq.'k') .or.
     :           el1(1:1) .eq. 'k' .or. el3(1:1).eq.'k') then
                 eltemp(1) = el1
                 eltemp(2) = el3
                 el1 = el2
                 el3 = el4
                 el2 = eltemp(1)
                 el4 = eltemp(2)
               end if
*
             ELSE
*
*	         ... O2 integrals ...
*
                READ(IUD, 3) END, K1, EL1, EL2, K2, EL3, EL4
                KVAL(I) = 64*K1 + K2
             END IF
             CPTR(I) = ICPTR + LAST
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
*	... Read the data ...
*
   20     READ(IUD,4) COEFF(IC), END, IH(IC), JH(IC), OPTR(IC)
          IF ( END .NE. '*') THEN
            IF (INT .LE. 2) THEN
              COEFF(IC) = ACURAT(COEFF(IC))
            ELSE
*
*	  ... Shift origin for overlap integrals
*
              IF (OPTR(IC).GT.512) THEN
                OPTR(IC) = INTPTR(3) + OPTR(IC) - 512
              ELSE IF (OPTR(IC).GT.0) THEN
                OPTR(IC) = INTPTR(2) + OPTR(IC)
              END IF
            END IF
            IC = IC + 1
            IF (IC .LE. NCDIM) GO TO 20
            STOP ' Too much data - current dimensions = (NCDIM)'
          END IF
*
*	... Initialize for next set ..
*
   18     INTPTR(INT) = I-1
          LAST = IC-1
   10 	CONTINUE
        RETURN
*
  999   WRITE(ERR,*)' Electron in ',END,'-data not found in ',
     :          'configuration list data'
        STOP
        END
*
*     ------------------------------------------------------------------
*          A O R T H O 
*     ------------------------------------------------------------------
*
*       This routine orthogonalizes the set of radial functions when an
*   orthogonality constraint applies.  A Gram-Schmidt type of  process
*   is used.  When more than one radial function with a given (nl) is
*   present, it may be necessary to solve a 2x2 system of equations.
*
*
      SUBROUTINE AORTHO
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500)
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
     :     ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
      LOGICAL DIAG
      COMMON AC(30,NWD),BC(NWD)
*
      IF (NWF .EQ. 1 .OR. IB .GT. NWF) RETURN
      I=NWF
      DIAG = .TRUE.
      IBEGIN = IEPTR(NWF-1)+1
      IP = IBEGIN
      IJ = 0
 60   JV = IJE(IP)
      IF (JV .LT. NWF .AND. IP .LE. IEPTR(I)) THEN
        IJ = IJ+1
        IF ( IJ .GT. (NWD)) STOP ' TOO MANY ORTHOGONALITY CONDITIONS'
        BC(IJ) = QUADR(NWF,JV,0)
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
62      CONTINUE
        IP = IP+1
        GO TO 60
      END IF
      IF ( IJ .GT. 0) THEN
        IF ( .NOT. DIAG .AND. IJ.GT.1) CALL LINEQN(30,IJ,AC,BC)
        M = MAX(NWF)
        AZZ = AZ(NWF)
        IP = IBEGIN
        CTOTAL = D0
        DO 65 JJ = 1,IJ
          C = BC(JJ)
          IF (DABS(C) .GT. 1.D-10) THEN
            WRITE(OUT,63) EL(IJE(IP)),EL(NWF),C
63          FORMAT(6X,'<',A3,'|',A3,'>=',1PD8.1)
            DO 64 J = 1,M
              P(J,NWF) = P(J,NWF) - C*P(J,IJE(IP))
64          CONTINUE
            AZZ = AZZ - C*AZ(IJE(IP))
          END IF
          IP = IP + 1
          CTOTAL = CTOTAL + ABS(C)
65      CONTINUE
        IF (CTOTAL .GT. 1.D-10 ) THEN
          AZ(NWF) = AZZ
          VARIED(NWF) = .TRUE.
        END IF
      END IF
      END
*
*     ------------------------------------------------------------------
*          A O U T P U T
*     ------------------------------------------------------------------
*
*       The radial functions and orthogonality integrals are printed,
*   if PRINT is .TRUE.   The  functions  will  also  be  punched  (or
*   stored) on unit OUF, if OUF .NE. 0.
*
*
      SUBROUTINE AOUTPUT(PRINT,DELTA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500)
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
     :     ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
      LOGICAL PRINT, REL/.false./
      DIMENSION POUT(8)
      IF ( .NOT. PRINT ) GO TO 31
C
C  *****  PRINT RADIAL FUNCTIONS, 7 PER PAGE
C
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
          POUT(IJ) = P(J,JJ)*R2(J)
9       CONTINUE
        K = K+1
        IF (K .LE. 10) GO TO 6
        K = 1
        KK = KK+1
        IF (KK .LT. 5) GO TO 21
        KK = 0
        WRITE(PRI,23)
23      FORMAT(1H1//)
        GO TO 6
21      WRITE(PRI,8)
8       FORMAT(1X)
        WRITE(PRI,10) R(J),(POUT(JJ),JJ=1,I)
6     CONTINUE
10    FORMAT(F13.5,F15.6,7F13.6)
      DO 15 J = ML,MU
        IJ = J - ML + 1
        POUT(IJ) = DPM(J)
15    CONTINUE
      WRITE(PRI,16) (POUT(J),J=1,I)
16    FORMAT(4X,10HMAX. DIFF. ,F15.7,7F13.7)
      ML = ML+8
      IF (ML .LE. NWF) GO TO 2
31    IF ( NWF .LE. 1) GO TO 30
*
*  *****  PRINT ORTHOGONALITY INTEGRALS
*
      WRITE(PRI,11) ATOM,TERM
11    FORMAT(////10X,33HORTHOGONALITY INTEGRALS FOR ATOM ,A6,6H TERM ,A6
     :   //20X, 4H(NL),3X,4H(NL),7X,8HINTEGRAL //)
      LM = IB
      ML = MAX0(2,LM)
      DO 12 I = ML,NWF
        JF = I - 1
        DO 13 J = 1,JF
          IF (L(I) .NE. L(J)) GO TO 13
          T = QUADR(I,J,0)
          WRITE(PRI,17) EL(I),EL(J),T
17        FORMAT(21X,A3,4X,A3,F15.8)
13      CONTINUE
12    CONTINUE
30    IF ( OUF .EQ. 0) GO TO 14
*
*  *****  Output functions on unit OUF for future input
*
*       EKI retained only for compatibility with MCHF format
*
      DO 3 I = NCLOSD+1,NWF
        IF (METH(I) .NE. 4) THEN
          EKI = -D5*HL(EL,I,I,REL)
        ELSE
          EKI = DELTA
        END IF
        MMX = MAX(I)
        WRITE (OUF) ATOM,TERM,EL(I),MMX,Z,E(I,I),EKI,AZ(I),
     :   (P(J,I),J=1,MMX)
3     CONTINUE
      WRITE (OUF) ATOM,TERM,EL(NWF),0,Z,E(NWF,NWF),EKI,AZ(NWF)
*
14    RETURN
      END
*     ------------------------------------------------------------------
*           A R A T E
*     ------------------------------------------------------------------
*
*	The ARATE subroutine computes an energy matrix, H-E.
*	From the Golden Rule, the autoionisation rate is calculated as
*	square of the interaction element:
*
*		Sum C * H  * A
*               i,j  i   ij   j
*
*	Where C and A is the coefficient vectors for the discrete and the
*	continuum states respectively. The sum over i runs from 1 to ID
*	(number of discrete configuration states) and the sum over j runs
*	from ID+1 to NCFG.
*	Strictly speeking, the whole matrix is not needed. It is still
*	calculated for future purposes. The time wasted is minor in most
*	cases.
*		
*
      SUBROUTINE ARATE(ETOTAL,ACFG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500,IDIM=1000,NCDIM=32768)
*
      INTEGER IN,OUT,ERR,PRI,OUC,OUD,OUF,OUH
      COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,OUC,OUD,OUF,OUH
*
      CHARACTER CONFIG*40,EL*3,ATOM*6,TERM*6,COUPLE*3
      COMMON /LABEL/CONFIG(NCD),EL(NWD),ATOM,TERM,COUPLE(NCD,9)
*
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :  ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
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
     :     ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
      COMMON WP(NCD),W(NCD,NCD)
*
      DO 1 I = 1,NCFG
        DO 2 J = 1,NCFG
          W(I,J) = D0
  2     CONTINUE
        W(I,I) = -ETOTAL
  1   CONTINUE
*
      IBEGIN = 1
      IEND = INTPTR(6)
      J = 0
      DO 10 I = IBEGIN,IEND
 11     IF (CPTR(I) .GT. J) THEN
          J = J + 1
          C = COEFF(J)*VALUE(I)
          IF (OPTR(J) .NE. 0) C = C*VALUE(OPTR(J))
          W(IH(J),JH(J)) = W(IH(J),JH(J)) + C
          GO TO 11
        END IF
 10   CONTINUE
*
*  ***** Symmetrize the matrix
*
      DO 12 I = 1,NCFG-1
        DO 13 J = I+1,NCFG
          W(I,J) = W(J,I)
 13     CONTINUE
 12   CONTINUE
      IF (TRACE) THEN
        WRITE(OUT,'(/2(10X,A,F16.8)/)')
     :          'ETOTAL =',ETOTAL,'K**2=',-ED
        WRITE(OUT,'(4X,6F12.7/(4X,6F12.7))')
     :                 (W(I,NCFG),I=1,NCFG-1)
      END IF
      WRITE(PRI,'(/2(10X,A,F16.8)//10X,A/)')
     :  'ETOTAL =',ETOTAL,'K**2 =',-ED,
     :  'Interaction elements (H)'
      DO 17 K = ID+1,NCFG
        WRITE(PRI,'((4X,6F12.7))') (W(I,K),I=1,ID)
 17   CONTINUE 
      WRITE(PRI,'(/10X,A/(4X,6F12.7))')
     :   'Coefficients - discrete',(WT(I),I=1,ID)
      WRITE(PRI,'(/10X,A/(4X,6F12.7))')
     :   'Coefficients - continuum',(WT(I),I=ID+1,NCFG)
      IF (ID .GE. 1) THEN
        V = D0
        DO 20 I = 1,ID
          DV = D0
          DO 21 J = ID+1,NCFG
            DV  = DV  + WT(J)*W(I,J)
21        CONTINUE
          V = V + WT(I)*DV
20      CONTINUE
        IF (V .NE. D0) THEN
          VSQ = V*V
          PI = ACOS(-D1)
          WIDTH = PI*VSQ
          A = 2.5976E+17*VSQ
          TAU = D1/A
          WRITE (10,19) V, WIDTH,27.21*WIDTH,219474.*WIDTH,A,TAU
          WRITE (PRI,19) V,WIDTH,27.21*WIDTH,219474.*WIDTH,A,TAU
19        FORMAT(/10X,'Auto-ionization Data'/
     :         10X,'--------------------'//
     :         20X,' Interaction, V       =',1PE13.4,' au'/
     :         20X,' Half-Line Width      =',1PE13.4,' au'/
     :         20X,'                       ',1PE13.4,' ev'/
     :         20X,'                       ',1PE13.4,' cm-1'/
     :         20X,' Auto-ionization rate =',1PE13.4,' s-1'/
     :         20X,' Half-Life            =',1PE13.4,' s'//)
        END IF
      END IF
*
      RETURN
      END
*
*     ------------------------------------------------------------------
*          A S C F
*     -----------------------------------------------------------------
*
*       This routine controls the 'SCF' procedure for the continuum
*	function through the following steps.

*    1. Set parameters.
*       If certain input parameters are zero (or blank) 
*       they will be set to their default value.
*
*         Parameter       Default Value (set in MAIN).
*         --------        -------------
*         SCFTOL          1.D-6
*         NSCF            20
*
*   	The self-consistency convergence criterion is
*
*          Z2 =  SCFTOL
*
*   	It is increased by a factor 1.3 at the end of each iteration
*
*    2. Call COULOM to calculate direct function in outer region.
*
*    3. Find initial estimate of the continuum function (by calling
*       CSOLVE) and orthogonalize (AORTHO).
*
*    4. Do SCF-iterations in the following steps (from 1 to NSCF unless
*	converged):
*            i) Calculates off-diagonal Energy parameters (AGRANGE).
*            ii) Solve the Differential Equation (CSOLVE).
*	     iii) Orthogonalize (AORTHO).
*            iv) Update the integrals (UPDATE).
*
*    5. Calculate the autoionization properties (ARATE) and do final
*       printouts (AOUTPUT and ASUMMRY).
*
      SUBROUTINE ASCF(ECORE,EKK,ACFG,SCFTOL,PRINT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500)
*
      INTEGER IN,OUT,ERR,PRI,OUC,OUD,OUF,OUH
      COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,OUC,OUD,OUF,OUH
*
      CHARACTER CONFIG*40,EL*3,ATOM*6,TERM*6,COUPLE*3
      COMMON /LABEL/CONFIG(NCD),EL(NWD),ATOM,TERM,COUPLE(NCD,9)
*
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
     :     ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
      COMMON /INPUT/ ECB,ICASE
*
      LOGICAL LAST,PRINT
      CHARACTER ANS*1
*
*  *****  Set the SCF convergence parameter to an optimistic value
*
      ETOTAL = ECORE - D5*EKK
      TOL = DSQRT(Z)*1.D-10
      Z2 = SCFTOL
      WRITE(OUT,15)
15    FORMAT(//)
      WRITE(OUT,16) OMIT,ACFG,Z2,NO,REL
   16 FORMAT(10X,44HWEAK ORTHOGONALIZATION DURING THE SCF CYCLE=,L4/
     :       10X,44HACCELERATING PARAMETER FOR MCHF ITERATION  =,F5.2/
     :       10X,44HSCF CONVERGENCE TOLERANCE (FUNCTIONS)      =,1PD9.2
     :      /10X,44HNUMBER OF POINTS IN THE MAXIMUM RANGE      =,I4/
     :       10X,44HRELATIVISTIC DIAGONAL  ENERGY CORRECTIONS  =,L4//)
      IPR = 0
*
* *****  Calculate the direct function in the asymptotic region.
*
      CALL COULOM(NWF,EKK)
*
* *****  Get initial estimates of the continuum function.
*
      IF (P(1,NWF) .EQ. D0) THEN
        CALL CSOLVE(NWF,DELTA)
        CALL AORTHO
      END IF
*
*  *****  Set iteration parameters
*
      LAST = .FALSE.
      ICYCLE = 0
      CALL UPDATE
      IF ( IB .GT. NWF ) GO TO 17
19    IF ( ID .EQ. 0 .OR. NCFG .EQ. 1) GO TO 9
*
*  *****  Perform NSCF self-consistent field iterations
*
9     DO 100 I = 1,NSCF
        ICYCLE = ICYCLE + 1
        WRITE(OUT,7) ICYCLE,Z2
7       FORMAT(//10X,17HITERATION NUMBER ,I2/10X,16H----------------/
     :     10X,'CONVERGENCE CRITERIA =',1PD9.1/)
        DP1 = D0
        IF (IB .GT. NWF) GO TO 17
        CALL AGRANGE
*
*  *****  Solve the differential equation
*
        WRITE(OUT,14)
14      FORMAT(/20X,' EL',6X,'ED/DELTA',10X,'AZ',11X,'NORM',7X,'DPM')
        CALL CSOLVE(NWF,DELTA)
        DP1 = DPM(NWF)*DSQRT(SUM(NWF))
        CALL AORTHO
        CALL UPDATE
        IF ( LAST ) GO TO 17
        IF ( I .EQ. NSCF ) GO TO 1
*
*  *****  If function appear to have converged,solve and test again
*
        IF (DP1 .LE. Z2) LAST =.TRUE.
1       CONTINUE
        WRITE(OUT,8) EL(NWF),DP1
8       FORMAT(/ 6X,34HLEAST SELF-CONSISTENT FUNCTION IS ,A3,
     :   27H :WEIGHTED MAXIMUM CHANGE =,1PD10.2)
        Z2=1.3*Z2
100   CONTINUE
18    WRITE(ERR,13)
13    FORMAT(10X/' SCF ITERATIONS HAVE CONVERGED TO THE ABOVE ACCURACY')
      WRITE(PRI,13)
      WRITE(ERR,*) ' Do you wish to continue ? (Y/N) '
      READ(IN,'(A)') ANS
      IF (ANS .EQ. 'Y' .OR. ANS .EQ. 'y') THEN
        WRITE(ERR,*) ' Enter the additional iterations '
        READ(IN,*) NSCF
        CALL UPDATE
        GO TO 19
      END IF
      FAIL = .TRUE.
*
*  *****  Perform final calculations
*
17    ACFG = D0
      IF(ICASE.NE.4) THEN
        CALL ARATE(ETOTAL,ACFG)
      END IF
      CALL AOUTPUT(PRINT,DELTA)
      CALL ASUMMRY(DELTA)
      NIT = NWF - IB + 1
      WRITE(PRI, 105) DELTA,DP1
105   FORMAT(//10X,'PHASE SHIFT FOR CONTINUUM FUNCTION =',F10.6/
     :         10X,'FINAL WEIGHTED CHANGE IN FUNCTION  =',D10.2)
      RETURN
      END
*
*     ------------------------------------------------------------------
*          A S U M M R Y
*     ------------------------------------------------------------------
*
*       The results of a calculation are summarized.   These include
*   the following for each electron:
*
*          E(NL)   - diagonal energy parameter
*          AZ(NL)  - starting parameter, P(r)/r**(l+1) as r -> 0.
*          SIGMA   - screening parameter as defined by Eq. (7-2).
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
*                      k   k   k
*   The values of all F , G , R  and <nl|L|n'l> integrals which enter
*   into the calculation are printed, but only if OUD > 0.
*
*
      SUBROUTINE ASUMMRY(DELTA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500,IDIM=1000,NCDIM=32768)
*
      INTEGER IN,OUT,ERR,PRI,OUC,OUD,OUF,OUH
      COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,OUC,OUD,OUF,OUH
*
      CHARACTER CONFIG*40,EL*3,ATOM*6,TERM*6,COUPLE*3
      COMMON /LABEL/CONFIG(NCD),EL(NWD),ATOM,TERM,COUPLE(NCD,9)
*
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
     :     ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
      COMMON R3(NWD),SS(3),R1,RM,RMM,RH,SC,QI,QJ,SP,C,CC,
     :   EKINP,EN,EPOT,RATIO,LI,LJ,K,KF,I1,I2,J1,J2,I,J,MIN
      CHARACTER*1 SYMBOL
*
      WRITE(PRI,'(///9X,A)') 'FINAL RESULTS:'
      WRITE(PRI,9) ATOM,TERM
9     FORMAT(/ 24X,'ATOM',1X,A6,3X,'TERM',1X,A6//
     : 2X,'nl',7X,'E(nl)',
     : 5X,'I(nl)/Delta',5X,'Rel(nl)',3X,'S(nl)',5X,'Az(nl)')
*
*  *****  Compute and print one-electron parameters
*
      DO 10 I = 1,NWF
        EK = DELTA
        S(I) = D0
        IF (METH(I) .NE. 4) THEN
          RM = QUADR(I,I,1)
          EK = -D5*HL(EL,I,I,REL)
          RH = 3*N(I)*N(I) - L(I)*(L(I) + 1)
          SC = Z - D5*RH/RM
          S(I) = SC
        END IF
        RELS = RLSHFT(I,I)
        WRITE (PRI,15)EL(I),E(I,I),EK,RELS,S(I),AZ(I)
15      FORMAT(1X,A3,F14.7,2F13.6,F7.2,F13.5)
10    CONTINUE
*
*  *****  Compute and print moments.
*
      WRITE(PRI,8)
 8    FORMAT(//2X,'nl',7X,'1/R**3',8X,'1/R',10X,'R',11X,'R**2')
      DO 11 I = 1,NWF
        R1 = D0
        RM = D0
        RMM = D0
        IF (METH(I) .NE. 4) THEN
          R1 = QUADR(I,I,-1)
          RM = QUADR(I,I,1)
          RMM = QUADR(I,I,2)
          R3(I) = D0
        END IF
        IF (L(I) .NE. 0) R3(I) = QUADR(I,I,-3)
        WRITE(PRI,16) EL(I),R3(I),R1,RM,RMM
16      FORMAT(1X,A3,F14.3,F13.4,F13.5,F13.5,F13.5)
11    CONTINUE
*
*  *****  Print tables of 'FK' and 'GK' integrals which were used in
*  *****  determining the energy
*
      IF ( OUD .EQ. 0 ) GO TO 13
      WRITE (OUD,126)
126   FORMAT(//2X,27HVALUES OF F AND G INTEGRALS        //)
      IBEGIN = 1
      IEND = INTPTR(2)
      DO 17 I = IBEGIN,IEND
        SYMBOL = 'F'
        IF (I .GT. INTPTR(1)) SYMBOL = 'G'
17    WRITE(OUD,19) SYMBOL,KVAL(I),EL(IEL(I,1)),EL(IEL(I,2)),VALUE(I)
19    FORMAT( 2X,A1,I2,1H(,A3,1H,,A3,4H ) =, F10.7)
*
*  *****  Print tables of 'RK' integrals
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
*  *****  Print tables of 'L' integrals
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
*     ------------------------------------------------------------------
*       C A L C A U T O
*     ------------------------------------------------------------------
*
*	Calcauto performs the following tasks:
*		1) some printouts.
*		2) determines maximum range of the discrete orbitals 
*		   (this will define the maximum R) by searching for an
*		    MX;	ABS(PJ(R)) < 10**(-5) for all J; 0 < J < NWF
*		    if R > R(MX).
*		3) it can follow two branches, depending on the value 
*		   of ICASE. 
*	           * If ICASE = 1 or 4:
*		     all data is known from the array WT and the
*	             variables ECB and ECORE and we can procede quickly
*		     to the scf-routine.
*      	 	   * If ICASE = 2 or 3: 
*		     it reads in data from unit 9 (.l or .j).
*		     For each J-value (if ICASE = 2) we do calculation
*                    for all states found (1 to MFOUND) and
*		     read in weights and discrete energy.
*
*
      SUBROUTINE CALCAUTO(ECORE,ACFG,SCFTOL,PRINT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500,IDIM=1000,NCDIM=32768)
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
     :     ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
      COMMON ZZ(NWD),IND(NWD),IELI(5),NOCCSH(NCD)
*
      COMMON /INPUT/ ECB,ICASE
*
      CHARACTER*72 HEADER,LSI*2,LSJ*2
      LOGICAL PRINT
* 
*  *****  Some printouts.
*
      IF(ICASE.NE.4) THEN
      WRITE(10,'(/10X,A/10X,A/10X,A//)')
     :                        '===================================='
     :                      , '    AUTOIONIZATION CALCULATION       '
     :                      , '===================================='
      WRITE(10,3) ECORE
3     FORMAT(//10X,'Energy of Target:',F14.8//
     :         10X,'Continuum State:'/10X,'===============')
      END IF
      DO 5 I = ID+1,NCFG
        WRITE(10,4) CONFIG(I),(COUPLE(I,J),J=1,9)
4       FORMAT(/10X,A40/10X,9(5X,A3)/)
5     CONTINUE
*
*  ***** Determine extent of outer region of bound functions
*
      MX = 0.
      DO 1 I = 1,NWF-1
         J = MAX(I)
 2       IF (ABS(P(J,I)) .LT. 1.D-5) THEN
            J = J-1
            GO TO 2
         END IF
         MX = MAX0(MX,J)
 1    CONTINUE
      SUM(NWF) = 1.
*
* ***** ICASE = 1:
* ***** Eigen vector (WT) and Discrete energy (EBOUND = ECB) from cfile.
*
      IF(ICASE.EQ.1) THEN
        ICFG = 1
        EBOUND = ECB
        V = D0

        DO 40 I = 1,ID
          DV = D0
          CALL LSTERM(COUPLE,I,LSI)
          DO 42 J = ID+1,NCFG
            CALL LSTERM(COUPLE,J,LSJ)
            IF (LSI .EQ. LSJ) DV = DV + ABS(WT(J))
 42       CONTINUE
          V = V + ABS(WT(I))*DV
 40     CONTINUE
        IF ( V .NE. D0) THEN
          EKK = 2*(ECORE-EBOUND)
          WRITE(10,6) -EKK/2
          WRITE(10,'(10X,A)') 
     :      '======================================='
          WRITE(10,'(//10X,A)') 'Discrete State, main component:'
          WRITE(10,'(10X,A)')   '------------------------------'
          WRITE(OUT,'(//A,F10.5/)') ' *** CALCULATIONS FOR k^2 =',-EKK
          WRITE(10,103) ICFG,CONFIG(ICFG),
     :                  (COUPLE(ICFG,J),J=1,9)
          WRITE(PRI,103) ICFG,CONFIG(ICFG),
     :                  (COUPLE(ICFG,J),J=1,9)
          IF(EKK.GT.D0) THEN
            WRITE(OUT,'(/A/)') 'THIS IS A BOUND STATE: E<ECORE'
            WRITE(10,'(/10X,A/)') 'This is a bound state: E<ECORE'
            GO TO 999
          END IF
          CALL EIJSET(NWF,NWF,EKK)
*
*  *****  Perform the MCHF iteration
*
          MAX(NWF) = MX
          CALL ASCF(ECORE,EKK,ACFG,SCFTOL,PRINT)
        END IF
*
* ***** ICASE = 4:
* ***** Hartree-Fock case - no discrete state defined.
*
      ELSE IF(ICASE.EQ.4) THEN
        EKK = -2*ECB
        WRITE(10,6) -EKK/2
        WRITE(10,'(10X,A)') 
     :    '======================================='
        WRITE(OUT,'(//A,F10.5/)') ' *** CALCULATIONS FOR k^2 =',-EKK
        CALL EIJSET(NWF,NWF,EKK)
*
*  *****  Perform the MCHF iteration
*  *****  ECORE is only used in ARATE, which will be bypassed.
*
        MAX(NWF) = MX
        CALL ASCF(ECORE,EKK,ACFG,SCFTOL,PRINT)
*
* ***** ICASE = 2(3) 
* ***** Read in weights and energies from .l (.j) file. 
* ***** Looping over J-values (goto 10) and states (do 100).
*
      ELSE
        READ(9,'(A)') HEADER
*
 10     READ(9,'(//8X,I4,10X,I4)',END=999) JV, MFOUND
        DO 100 M = 1,MFOUND
*         READ(9,101) ICFG,EBOUND,(WT(I),I=1,NCFG-1)
          READ(9,101) ICFG,EBOUND,(WT(I),I=1,ID)
101       FORMAT(/I6,F16.8/(7F10.6))
          V = D0
          DO 50 I = 1,ID
            DV = D0
            CALL LSTERM(COUPLE,I,LSI)
            DO 52 J = ID+1,NCFG
              CALL LSTERM(COUPLE,J,LSJ)
              IF (LSI .EQ. LSJ) DV = DV + ABS(WT(J))
 52         CONTINUE
            V = V + ABS(WT(I))*DV
 50       CONTINUE
          IF ( V .NE. D0) THEN
            EKK = 2*(ECORE-EBOUND)
            WRITE(10,6) -EKK/2
            WRITE(10,'(10X,A)') 
     :      '======================================='
            WRITE(10,'(//10X,A)') 'Discrete State, main component:'
            WRITE(10,'(10X,A)')   '------------------------------'
            WRITE(OUT,'(//A,F10.5/)') ' *** CALCULATIONS FOR k^2 =',-EKK
            IF (ICASE.EQ.2) THEN
              WRITE(10,103) ICFG,CONFIG(ICFG),
     :                  (COUPLE(ICFG,J),J=1,9)
              WRITE(PRI,103) ICFG,CONFIG(ICFG),
     :                  (COUPLE(ICFG,J),J=1,9)
            ELSE
              WRITE(10,102) ICFG,CONFIG(ICFG),JV,
     :                  (COUPLE(ICFG,J),J=1,9)
              WRITE(PRI,102) ICFG,CONFIG(ICFG),JV,
     :                  (COUPLE(ICFG,J),J=1,9)
            END IF
            IF (EKK.GT.D0) THEN
              WRITE(OUT,'(/A/)') 'THIS IS A BOUND STATE: E<ECORE'
              WRITE(10,'(/A/)') 'This is a bound state: E<ECORE'
              GO TO 100
            END IF
 6          FORMAT(//10X,'Continuum electron: k^2/2=',F10.5,' au.')
102         FORMAT(/3X,I3,6X,A40,'2*J = ',I2/12X,9(5X,A3))
103         FORMAT(/3X,I3,6X,A40/12X,9(5X,A3))
            CALL EIJSET(NWF,NWF,EKK)
*
*  *****  Perform the MCHF iteration
*
            MAX(NWF) = MX
            CALL ASCF(ECORE,EKK,ACFG,SCFTOL,PRINT)
          END IF
100     CONTINUE
        GO TO 10
      END IF
999   RETURN
      END
*
*     ------------------------------------------------------------------
*       C N M R V
*     ------------------------------------------------------------------
*
* 	CNMRV Solves the differential equation
*
*               y" = FK y + X
*	
*	In two different regions - inner (1,NJ+1) and outer (NJ+2,M) - by
*	outwards integration.
*
*
      SUBROUTINE CNMRV(NJ,M,M0,AZZ,FK,X,FH,XH,PH,PDE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500)
      DIMENSION PDE(*),FK(*),X(*),FH(*),XH(*),PH(*)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :  ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
c
c	starting from JH=5, JH corresponding to J=NJ+2, 
c	two points each J, JH=2*(J-NJ) + 1
c	R(J) --> rh(JH)
cxi	comment by Jinhua Xi

*
*
*    ....INNER REGION .....
*
*  *****  Integrate outward to NJ+1
*
      Y1 = PDE(M0)
      Y2 = PDE(M0+1)
      G1 = FK(M0)
      G2 = FK(M0+1)
      DO 11 J = M0+2,NJ+1
        G3 = FK(J)
        Y3 = (Y2+Y2-Y1 + (D10*G2*Y2 + G1*Y1) + X(J-1)) / (D1 - G3)
        PDE(J) = Y3
        Y1 = Y2
        Y2 = Y3
        IF (ABS(Y3) .GT. 1.E5) THEN
*
*         ....Scale down to avoid overflow
*
          AZZ = AZZ/100.D0
          DO 12 JJ = 1,J
            PDE(JJ) = PDE(JJ)/100.D0
12        CONTINUE
          Y1 = Y1/100.D0
          Y2 = Y2/100.D0
        END IF
        G1 = G2
      G2 = G3
11    CONTINUE
*
*    ....OUTER REGION  ....
*
* ***** Redefine Y2 for half the step-size
*
      G1 = FH(1)
      G2 = FH(2)
      G3 = FH(3)
      Y2 = (-XH(2) + Y1+Y3 - G1*Y1 - G3*Y3)/(D2 + D10*G2)
*
*     ...Integrate outwards from NJ+2 to M.
*
      Y1 = Y2
      Y2 = Y3
      G1 = G2
      G2 = G3
      JH = 4
      DO 20 J = NJ+2,M
        DO 22 K = 1,2
          G3 = FH(JH)
          Y3 = (Y2+Y2-Y1 + (D10*G2*Y2 + G1*Y1) + XH(JH-1)) / (D1 - G3)
          ph(jh) = y3
          Y1 = Y2
          Y2 = Y3
          G1 = G2
          G2 = G3
          JH = JH+1
22      CONTINUE
        PDE(J) = Y3
20    CONTINUE
      END
*
*----------------------------------------------------------------------
*	C O E F
*----------------------------------------------------------------------
*
*	This function returns the coefficient of a given Slater
*  integral in the expansion of the energy
*
        DOUBLE PRECISION FUNCTION COEF(INT)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
        INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
        PARAMETER (NCD=500,IDIM=1000,NCDIM=32768)
        COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :       ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
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
*     ------------------------------------------------------------------
*       		C O U L O M
*     ------------------------------------------------------------------
*
*	This subroutine handles three different tasks:
*
*	  1) Computes and stores the FK function;
*
* 	  	FK(r) =[-2(Z-YR(r))*r + (l+.5)**2 + (k*r)**2]CH
*
*	     It will remain constant for the rest of the calculations.
*
*	  2) Determines the onset of the asymptotic, Coulomb region, as
*	     the classical turning point, R(NJ).
*
*         3) Determines the point, R(MP)>R(NJ), outside which the potential
*            could be described as Zeff/r, ie the change in YR is negligible.
*
*	  4) Interpolates FK, to obtain and store FH, in the asymptotic
*	     region. (That is between R(NJ) and R(MM), where MM =
*	     MIN(NO-2,NJ+99)) The FK functions is determined in the normal
*	     grid points, while the FH is defined for Half the stepsize.
*
      SUBROUTINE COULOM(I,EKK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :  ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
      COMMON /CONTIN/FK(NOD),FH(260),XH(260),rh(260),r2h(260),ph(260),
     :               CD,FL,ZL,ZF,V,NJ,MJ,MP,IX
      COMMON /COULFG/FKC(NOD),XC(NOD),PDC(NOD),FHC(260),PHC(260),NJM
*
      LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
      COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)

*
*  *****  Set arrays for the asymptotic region
*
*
      JJ = 0
      M = NOD
      MAX(I) = M
      CALL POTL(I)
      FL = L(I)
      ZL = Z/(FL+D1)
      ZF = Z - YR(NO)
      V = YR(1)/R(1)
      CD = (FL+D5)**2
      yrno = yr(no)

*
*        1) Compute the direct function, FK.
*
      NJM=0
      DO 4 J = 1,M

cc      write(76, '(2f15.7)') r(j), z-yr(j)
      if( dabs( yrno-yr(j) ) .lt. 1.d-7 .and. NJM .eq. 0 ) NJM=J

      FKC(J) = (-D2*ZF*R(J) + CD + EKK*RR(J))*CH
4     FK(J) = (-D2*(Z - YR(J))*R(J) + CD + EKK*RR(J))*CH
*
*        2) Search for the point at which FK(J)<0 for J>NJ
*
      NJ = M
5     IF ( FK(NJ) .LT. D0 ) THEN
        NJ = NJ-1
        IF (NJ .GT. 90 ) GO TO 5
      END IF
      NJ = NJ+1

*
*        3) Search for the point, MP, outside which YR is constant
*           and we can define a Zeff.
*
      MP = M
6     IF ( YR(MP) - YR(MP-1) .LT. r(mP)*1.D-4 ) THEN
        MP = MP-1
        IF (MP .GE. NJ) GO TO 6
      END IF
      MP = MP+1
*
*        4) Interpolate FK to obtain FH in "half-grid-points".
*
      EXPH = EXP(H/D2)
      CHH = CH/D4
      JH = 1
      MM = MIN0(NO-2,NJ+129)
      DO 50 J = NJ,MM
        FH(JH) = FK(J)/D4
        FHC(JH) = FKC(J)/D4
        YRH = (9.*(YR(J)+YR(J+1))-YR(J-1)-YR(J+2))/D16
        RHH = R(J)*EXPH
        RRH = RHH*RHH
        FH(JH+1) = (-D2*(Z-YRH)*RHH +CD + EKK*RRH)*CHH
        FHC(JH+1) = (-D2*ZF*RHH +CD + EKK*RRH)*CHH
        JH = JH+2
50    CONTINUE
      IX = 0
      END
*
*     ------------------------------------------------------------------
*       C S O L V E
*     ------------------------------------------------------------------
*       Modified by Jinhua Xi, December 1994
*       CSOLVE performs the following tasks:
*
*         1) Computes the Exchange function, X(r), by calling CXCH.
*	  2) Redefines the outer region;
*            * Initialization stage (IX = 0) - first call:
*	       reduce outer region (MAX(NWF) -> M) if step size is too
*              big (R(I) - R(I-1) > 2/K) and where NJ <= M <= MAX(NWF).
*              M-NJ is also, by dimensions, bound to be less than 130.
*            * During SCF-cycle (IX = 1):
*              find MJ; MP < MJ < M and 
*			X(r)/FK(r) < 0.0025 if R > R(MJ).
*              This will be used to calculate phase shift. If MJ = M 
*              the X function might be truncated -> warning.
*	  3) Interpolates X(r) to XH(r), for half the step size, in 
*            the azymptotic region; r(NJ) <= r <= r(M).
*	  4) Prepares for the Numerov method of solving the 
*	     differential equation, both in azymptotic and inner region.
*	  5) Solves the differential equation, by calling CNMRV.
*	  6) Calculates the phase shift and renormalization, by using
*	     regular, FC, and irregular, GC, Coulomb functions; 
*	     Renormalizes the continuumfunction.
*	  7) Orthogonalize the continuum function.
*	The following subroutines are called in the different
*	steps:
*	  1) CXCH
*	  5) CNMRV
*	  6) FGCOUL 
*	
*
      SUBROUTINE CSOLVE(I,DELTA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500)
*
      INTEGER IN,OUT,ERR,PRI,OUC,OUD,OUF,OUH
      COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,OUC,OUD,OUF,OUH
*
      CHARACTER CONFIG*40,EL*3,ATOM*6,TERM*6,COUPLE*3
      COMMON /LABEL/CONFIG(NCD),EL(NWD),ATOM,TERM,COUPLE(NCD,9)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :  ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
      LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
      COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
      COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :     ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
      COMMON /CONTIN/FK(NOD),FH(260),XH(260),rh(260),r2h(260),ph(260),
     :               CD,FL,ZL,ZF,V,NJ,MJ,MP,IX
      COMMON /COULFG/FKC(NOD),XC(NOD),PDC(NOD),FHC(260),PHC(260),NJM
*
      PI = ACOS(-D1)
      ED = E(I,I)
      ekk = -ed
      eK = SQRT(ekk)

      icase = 0
      AZD = AZ(I)
      IPR = NWF
      M = MAX(NWF)
*
* ***** Computes the Exchange function.
*
        CALL CXCH(I,3)
cxi
cxi	save X(J)  for repeated use with different
cxi	E(I,J) values
cxi
        do 12 j=1,no
          p(j,nwf+1) = X(j)
12      continue

*
* ***** Redefine the outer region.
*

678     IF (IX .EQ. 0) THEN
5       IF (eK*(R(M)-R(M-1)) .GT. 2.0D0) THEN
          M = M-1
          IF (M .GT. MP) GO TO 5
        END IF
        IF (M .LT. MAX(NWF)) THEN
          WRITE(OUT,*) 'Outer region reduced by ',MAX(NWF)-M,'points',
     :                 ' in CSOLVE. '
        END IF
        MJ = M
        print *, ' End of range: ', r(m), ' au.'
      ELSE
2       IF (ABS(X(MJ)/FK(MJ)) .LT. 0.0025) THEN
          MJ = MJ-1
          IF (MJ .GT. MP) GO TO 2
        END IF
      END IF
      IF (M-NJ .GE. (130)) THEN
        WRITE(OUT,*) ' Outer region contains ',M-NJ+1,' points'
        WRITE(OUT,*) ' Maximum allowed value is (130)'
        STOP 'In CSOLVE - Too large outer region'
      END IF
      IF (MJ .EQ. M .AND. IX.NE.0) THEN
        if ( icase .eq. 0 ) then 
          WRITE(OUT,'(/1X,A,I4)')
     :      'WARNING: Outer region may be truncated. M = MJ =',M
          WRITE(OUT,'(A)') 
     :      '          Exchange function not small at M!'
        endif
      ELSE
        IX = 1
      END IF
*
* ***** Interpolate X in outer region to obtain XH.
*
      ehh = exp(h/2.d0)
      e2hh = exp(h/4.d0)
      JH = 1
      DO 3 J = NJ,M
        XH(JH) = X(J)
        rh(jh) = r(j)
        r2h(jh) = r2(j)
        XH(JH+1) = (9.*(X(J)+X(J+1)) -X(J-1)-X(J+2))/D16
        rh(jh+1) = r(j)*ehh
        r2h(jh+1) = r2(j)*e2hh
        JH = JH + 2
3     CONTINUE
*
* ***** Prepare for Numerovs method.
* ***** i) Compute the RHS of the Numerov equation for the outer region
*
      CHH = CH/D4
      JH = JH-2
      X1 = XH(1)
      X2 = XH(2)
      DO 4 J = 2,JH
        X3 = XH(J+1)
        XH(J) = CHH*(X1 + D10*X2 +X3)
        X1 = X2
        X2 = X3
4     CONTINUE
*
* ***** ii) Compute the RHS of the Numerov equation
*
      XY = X(1)
      XP = X(2)
      X1 = X(1)
      X2 = X(2)
      X3 = X(3)
      X4 = X(4)
      DO 1 J = 3,NJ
        X5 = X(J+2)
        X(J) =CH*(-X5+24.D0*(X4+X2) + 194.D0*X3 - X1)/20.D0
        X1 = X2
        X2 = X3
        X3 = X4
        X4 = X5
1     CONTINUE
*
* ***** iii) Add the deferred difference correction to the exchange
* *****      for the outward integration region
*
      X1 =    P(1,I)*FK(1)
      X2 =    P(2,I)*FK(2)
      X3 =    P(3,I)*FK(3)
      X4 =    P(4,I)*FK(4)
      DO 7 J = 3,NJ
        X5 =     P(J+2,I)*FK(J+2)
        X(J) = X(J) - (X5 -D4*(X2 + X4) + D6*X3 +X1)/20.D0
        X1 = X2
        X2 = X3
        X3 = X4
        X4 = X5
7     CONTINUE
      RL = L(I) + 2.5
      X(2) = R(2)**RL*(X(5)/R(5)**RL - D3*(X(4)/R(4)**RL -
     :     X(3)/R(3)**RL))
*
* ***** iv) Compute starting values from series expansion
*
      CC = D2*FL + D3
      A2 = (V + ED/D2 + ZL*Z)/CC
      A3 = -((V + ED/D2)*ZL + Z*A2)/(D3*(FL+D2))
      DO 6 J = 1,2
6     PDE(J) = AZ(I)*R(J)**L(I)*R2(J)*
     :                           (R(J)*(R(J)*(R(J)*A3 + A2) -ZL) +D1)

      PDE(1) = PDE(1) + XY/(D2*CC)
      PDE(2) = PDE(2) + XP/(D2*CC)

*
*  *****  Solve the differential equation.
* 
cxi
      CALL CNMRV(NJ,M,1,AZ(I),FK,X,FH,XH,PH,PDE)


      CALL FGCOUL(I,M,CN,DELTA)
cxi
cxi	note: the sign of the overlop depends on cn, 
cxi	because we need to adjust the off-diagonal parameters
cxi	repeatedly to find a zero overlap, we can not change the sign of
cxi	the wavefunction
cxi	so, save the sign in cnn
cxi
      if ( cn .lt. 0.d0 ) then
        cnn = -1.d0
        cn = - cn
      else
        cnn = 1.d0
      endif
*
* ***** Renormalise.
*
      DO 92 J = 1,M
        PDE(J) = PDE(J)*CN
92    CONTINUE
      AZ(I)=CN*AZ(I)

      DO 13 J = 1,M
        DIF = P(J,I) - PDE(J)
        P(J,I) = PDE(J) + acc(i)*dif
13    CONTINUE
      AZ(I) = (1.d0 - acc(i))*AZ(i) + acc(i)*AZD
cxi
cxi	if the acc(i) is used, the wavefunction needs to be
cxi	re-normalized
cxi
      MAX(I) = M
*
* *****  Orthogonalize the continuum function
*
      IBEGIN = 1
      IF (I .GT. 1) IBEGIN = IEPTR(I-1)+1
      IP = IBEGIN
50    JI = IJE(IP)
      IF ( JI .NE. I) THEN
        CC = QUADR(I,JI,0)

        if ( ip .eq. ibegin ) then
cxi
cxi the first orthogonal requirements, use new approach
cxi
          if ( e(i,ji) .eq. 0.d0  .or. dabs(cc) .lt. 1.d-7) goto 65

          azcc= CC*AZ(JI)/az(i)
          if ( dabs(azcc) .gt. 0.1d0) goto 66 
        endif

65      WRITE(OUT,63) EL(JI),EL(I),CC, e(i,ji)
63      FORMAT(1X,'<',A3,'|',A3,'>=',D8.2, '  eij=',d18.10)

        AZ(I) = AZ(I) - CC*AZ(JI)
        DO 51 J = 1,M
          P(J,I) = P(J,I) - CC*P(J,JI)
51      CONTINUE
c
cxi	re-normalize the wavefunction
cxi
        cnn=cnn*pde(M)/P(M,i)
        if (az(i) .lt. 0.d0) cnn = - cnn 
        Az(i) = cnn*Az(i)
        do 52 j = 1,m
          p(j,i) = cnn*p(j,i)
52      continue
cxi
cxi	normalize the off diagonal parameters
cxi
        eijv = cn*cnn*e(i,ji)
        call eijset(i,ji,eijv)

        IP = IP+1
        IF (IP .LE. IEPTR(I)) GO TO 50
      END IF

*
      VARIED(I) = .TRUE.
      DP = ABS((AZD - AZ(I))/AZ(I))
      DPM(I) = DP
      WRITE(OUT,17) EL(I),DELTA,AZ(I),CN,'c',DP
17    FORMAT(20X,A3,2F15.7,F12.7,A2,1PD10.2)

      return
cxi
cxi	adjusting the off-diagonal parameters E(I,JI)
cxi
66    eijv = e(i,ji)
      call geteij(icase,eijv,cc)
      deteij = eijv - e(i,ji) 
      call eijset(i,ji,eijv)
cxi
cxi	for new E(i,j) values, always keep updating 
cxi	the X function
cxi
      do 18 j=1,no
        p(j,nwf+1) = p(j,nwf+1) + deteij*P(j,ji)*rr(j)
        x(j) = p(j,nwf+1)
18    continue
      az(i)= azd

      goto 678
      END

*
*----------------------------------------------------------------------
*               C X C H
*----------------------------------------------------------------------
*
*	  This subroutine computes the function X(r) for the continuum
*	orbital. It includes contributions from Exchange, Configuration
*	Interaction, One-electron part and off-diagonal Lagrange
*	Multipiers.
*
      SUBROUTINE CXCH(I,IOPT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500,IDIM=1000,NCDIM=32768)
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
     :     ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
      DO 1 J=1,NO
  1   X(J) = D0
      DO 2 J = 1,NWF
        IF ((I.LE.NCLOSD .AND. I.NE.J) .OR.
     :         (I.GT.NCLOSD .AND. J.LE.NCLOSD))  THEN
          DO 4 K = IABS(L(I)-L(J)),L(I)+L(J),2
            C = - D2*CB(L(I),L(J),K)*SUM(J)
            CALL YKF(J,I,K,REL)
            DO 6 JJ = 1,NO
              X(JJ) = X(JJ) + C*YK(JJ)*P(JJ,J)
  6         CONTINUE
  4       CONTINUE
        END IF
  2   CONTINUE
      SUMI = SUM(I)
      IF (I .LE. NCLOSD) GO TO 71
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
 8        CONTINUE
        END IF
 7    CONTINUE
*
 71   GO TO (75,76,77),IOPT
 76   DO 78 J = 1,NO
        X(J) = X(J)/R(J)
 78   CONTINUE
      GO TO 75
 77   DO 79 J =1,NO
        X(J) = R(J)*X(J)
 79   CONTINUE
      DO 74 J = 1,NWF
        IF (J .NE. I) THEN
          C = E(I,J)
          IF (DABS(C) .LE. 1.D-20 ) GO TO 74
          DO 73 JJ = 1,NO
 73         X(JJ) = X(JJ) + C*P(JJ,J)*RR(JJ)
        END IF
 74   CONTINUE
C
C  *****  Check if exchange is zero: if so, method 2 should be used.
C
 75   IF (METH(I) .GE. 2) RETURN
      IF ( DABS(X(1)) + DABS(X(2)) + DABS(X(3)) .EQ. D0 ) METH(I) = 2
      END
*
*-----------------------------------------------------------------------
* 		E
*-----------------------------------------------------------------------
*      Returns the value of the off-diagonal energy parameter
*  for the (i,j) pair from the data structure.
*
        DOUBLE PRECISION FUNCTION E(I,J)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220,NCD=500)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :     ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
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
*		E I J S E T
*-----------------------------------------------------------------------
*
*      Stores the value of the off-diagonal energy parameter for the
*    pair (i,j) in the data structure
*
        SUBROUTINE EIJSET(I,J,E)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220,NCD=500)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
*
      COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :    ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
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
     :     STOP ' Too many off-diagonal energy parameters'
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
*     ------------------------------------------------------------------
*        L S T E R M
*     ------------------------------------------------------------------
*
*       Determine the LS term value from the COUPLE array
*
        SUBROUTINE LSTERM(COUPLE,I,LS)
        PARAMETER (NWD=30,NOD=220,NCD=500)
*
        INTEGER IN,OUT,ERR,PRI,OUC,OUD,OUF,OUH
        COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,OUC,OUD,OUF,OUH
        CHARACTER*3 COUPLE(NCD,9),LS*2
*
        J = 9
  10    IF (COUPLE(I,J) .EQ. '   ') THEN
          J = J-1
          IF ( J .GT. 0) THEN
            GO TO 10
          ELSE
            WRITE(ERR,*) 'COUPLING information missing'
            STOP
          END IF
        END IF
        LS = COUPLE(I,J)(1:2)
        IF (LS(2:2) .GE. 'a' .AND. LS(2:2) .LE. 'z')
     :      LS(2:2) = CHAR(ICHAR(LS(2:2)) + ICHAR('A') - ICHAR('a'))
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
        PARAMETER (NWD=30,NOD=220,NCD=500)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
        COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
        INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
        PARAMETER (IDIM=1000,NCDIM=32768)
        COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :       ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
        LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
        COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
        DO 1 J=1,NO
1       YR(J) = D0
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
*---------------------------------------------------------------------
*		S E T O R T
*---------------------------------------------------------------------
*
*	Determine if orbitals for electrons (el1, el2) should be
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
*-----------------------------------------------------------------------
*		U P D A T E
*-----------------------------------------------------------------------
*
*	Evaluate all integrals where at least on orbital has changed.
*
        SUBROUTINE UPDATE
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (NWD=30,NOD=220,NCD=500)
*
        INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
        PARAMETER (IDIM=1000,NCDIM=32768)
*
        CHARACTER CONFIG*40,EL*3,ATOM*6,TERM*6,COUPLE*3
        COMMON /LABEL/CONFIG(NCD),EL(NWD),ATOM,TERM,COUPLE(NCD,9)
        COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :       ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
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
        PARAMETER (NWD=30,NOD=220,NCD=500)
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

      write(78,'(2f15.7)') R(j), P(J,i)


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
*     ------------------------------------------------------------------
*      FGCOUL  
*     ------------------------------------------------------------------
*
      SUBROUTINE FGCOUL(I,M,cn,delta)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500)
      DOUBLE PRECISION K
*
      INTEGER IN,OUT,ERR,PRI,OUC,OUD,OUF,OUH
      COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,OUC,OUD,OUF,OUH
      CHARACTER CONFIG*40,EL*3,ATOM*6,TERM*6,COUPLE*3
      COMMON /LABEL/CONFIG(NCD),EL(NWD),ATOM,TERM,COUPLE(NCD,9)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :  ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
      LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
      COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
      COMMON /CONTIN/FK(NOD),FH(260),XH(260),rh(260),r2h(260),ph(260),
     :     CD,FL,ZL,ZF,V,NJ,MJ,MP,IX
      COMMON /COULFG/FKC(NOD),XC(NOD),PDC(NOD),FHC(260),PHC(260),NJM
      dimension xhc(260)
      data ED0,LL0 / 1.d0, -1 / 

      PI = ACOS(-D1)
      PII = PI + PI
      ED = E(I,I)
      LL = L(I)

      K = SQRT(-ED)

      ifirst = 0
      if( ed .eq. ed0 .and. ll .eq. ll0 ) ifirst = 1

      if( ifirst .eq. 1 ) goto 10 

      it =1 
      azd = d1
      do j=1,no
        pdc(j)= d0
        xc(j)=0.d0
      enddo

      ed0 = ed
      ll0 = ll

*
* ***** iv) Compute starting values from series expansion
*
      CC = D2*FL + D3
      ZLF= ZF/(FL+D1)
      A2 = ( ED/D2 + ZLF*ZF)/CC
      A3 = -(( ED/D2)*ZLF + ZF*A2)/(D3*(FL+D2))
      DO 6 J = 1,2
6     PDC(J) = AZD*R(J)**LL*R2(J)*
     :               (R(J)*(R(J)*(R(J)*A3 + A2) -ZLF) +D1)
*
*  *****  Solve the differential equation.
*

       MM = MIN0( NO-2, NJ+129)
       do j=1,260
       xhc(j)=0.d0
       enddo

2      CALL CNMRV(NJ,MM,1,AZD,FKC,XC,FHC,XHC,PHC,PDC) 

      if( it .eq. 3 ) goto 10 

*
*
* ***** iii) Add the deferred difference correction to the exchange
* *****      for the outward integration region
*
      it = it + 1

      X1 =    PDC(1)*FKC(1)
      X2 =    PDC(2)*FKC(2)
      X3 =    PDC(3)*FKC(3)
      X4 =    PDC(4)*FKC(4)
      DO 7 J = 3,NJ
        X5 =     PDC(J+2)*FKC(J+2)
        XC(J) =  - (X5 -D4*(X2 + X4) + D6*X3 +X1)/20.D0
        X1 = X2
        X2 = X3
        X3 = X4
        X4 = X5
7     CONTINUE
      RL = LL + 2.5
      XC(2) = R(2)**RL*(XC(5)/R(5)**RL - D3*(XC(4)/R(4)**RL -
     :     XC(3)/R(3)**RL))

      goto 2

*
* ***** From the Coulomb functions we calculate the phase shift and
* ***** the amplitude.
*
10      H60 = 2.d0/(60.d0*H)

        MJ = MAX0(NJ+2,NJM)
        MJ = MAX0(MJ,MP)
        if( MJ .gt. M - 4 ) MJ = M - 4


95      MJ = MJ + 2

        if( MJ .gt. M-2 ) then
          write(*,*) 
     :    ' Failed in Calculating Coulomb Functions, R(M) too small '
          write(*,*) ' R(M) =', R(M)
          stop
        endif

      MJH = 2*(MJ-NJ)+1
cc
      DO 90 J = MJH-1,MJH
        xval = RH(j)
        yval = PHC(J)*R2h(J)


c
c	y'(r)
c
        dyval = ( PHC(J+3) - PHC(J-3) - 9.d0*( PHC(J+2)-PHC(J-2) ) 
     :           + 45.d0*( PHC(J+1) - PHC(J-1) ) )*h60
        dyval = ( dyval + 0.5d0*PHC(J) )/R2h(J)

        CALL  fgwkb(k,zf,ll,xval,yval,dyval,f2,df,g2,dg,ierr)
        if( ierr .ne. 0 ) goto 95

        if( J .eq. MJH-1 ) then 
          F1 = f2 
          G1 = G2
        endif

90    CONTINUE

      DNUM = R2h(MJH)*PH(MJH)*F1 - R2H(MJH-1)*PH(MJH-1)*F2
      DENO = R2h(MJH-1)*PH(MJH-1)*G2 - R2H(MJH)*PH(MJH)*G1
      DELTA = ATAN( DNUM/DENO )
      AMP = R2H(MJH)*PH(MJH)/(F2 + G2*(DNUM/DENO) )
      CN = COS(DELTA) / AMP 
      END
*     ------------------------------------------------------------------
*           F G W K B
*     ------------------------------------------------------------------
*
*	This routine determines the energy normalized 
*	Coulomb functions F , and G,
*	and their derivatives F' , G' , with respect to r
*
*	Normalization obtained using the  WKB method as proposed by 
*	Liu, Xi, and Li, PRA48, 228(1993)
*
*       Written by Jinhua Xi, December, 1994
*     ------------------------------------------------------------------
*
*	F  =  sin( phi) * sqrt(2/pi k)
*
        subroutine fgwkb(ek,z,l,r,yr,dyr,f,df,g,dg,ierr)
        Implicit double precision (a-h,o-z)
        data pi,pii/3.141592653589793d0,6.283185307179586d0/

        ierr=0
        call  wkb(ek,z,l,r,zeta,dz,deltaz)


        if( dabs(deltaz) .gt. 1.d-4 ) ierr=1
        if( dabs(deltaz) .gt. 1.d-1) then
           ierr=2
           return
        endif
*
*	determine the phase function phi(r) and the
*	normalization constant 
*
        dzz = 0.5d0*dz/zeta
        pn= dsqrt(2.d0/pi/zeta)
*

        phi = atan( zeta/(dyr/yr + dzz) )
        if( sin(phi)*yr .lt. 0.d0 ) phi =phi + pi 

        if( phi .lt. 0.d0 ) phi=phi+pii
        if( phi .gt. pii ) phi = phi - pii  
c
c	for Coulomb potential , get the F, G, F',G'
c
        f  =   dsin(phi)*pn
        g  =   dcos(phi)*pn
        df =   zeta*g - dzz*f
        dg = - zeta*f - dzz*g
        return 
        end
*
*     ------------------------------------------------------------------
*           W K B
*     ------------------------------------------------------------------
*
*	This routine performs the WKB iteration proposed by 
*	Liu, Xi, and Li, PRA48, 228(1993)
*
*       The exact formulas for the iterative procedure were
*       derived by C. F. Fischer, using the MAPLE Symbol 
*       manipulation package.
*	
*	This routine is called by asympn and fgwkb.
*	asympn: determines phase and normalization,
*		as proposed by LXL paper
*	fgwkb:  computes only f,g, f',g' 
*
*       Written by C. F. Fischer, July, 1994
*       Modified by Jinhua Xi, December, 1994
*     ------------------------------------------------------------------

        subroutine wkb(ek,z,l,r,zeta,dz,deltaz)
        Implicit double precision (a-h,o-z)
        double precision w(0:8), u(0:8,0:4)
cxi
cxi	it is ABSOLUTELY necessary to initiate the arrays
cxi	if you would like to have a correct result. 
cxi	
        do 1 i=0,8
          w(i)=0.d0
          do 2 j=0,4
            u(i,j)=0.d0
  2       continue
1       continue
cxi
cxi	...................................................
cxi
        a= 2*z/r
        b= -l*(l+1)/r/r
        ekk=ek*ek
        w(0) = ekk + a + b
        if( w(0) .lt. 0.3d0*ekk ) then
          write(*,*) '  in WKB: r-value too small, r=', r
          deltaz=99.d0
          return
        endif

        u(0,0) = w(0)
        do 10 i = 1,8
          a = -i*a/r
          b = -(i+1)*b/r
          w(i) = a+b
          u(i,0) = w(i)/w(0)
10      continue
        do 20 j = 0,3
           u(0,j+1) = w(0) + (5*u(1,j)**2 -4*u(2,j))/16.d0
           do 30 i = 1, 6-2*j

             if (i .eq. 1) then
               u(1,j+1) = 7*u(1,j)*u(2,j)-5*u(1,j)**3 -2*u(3,j)

             else if (i .eq. 2) then
               u(2,j+1) = 7*u(2,j)**2 -29*u(1,j)**2*u(2,j) + 
     :                    9*u(1,j)*u(3,j) +15*u(1,j)**4 -2*u(4,j)

             else if (i .eq. 3) then
               u(3,j+1) = 23*u(2,j)*u(3,j) -72*u(2,j)**2*u(1,j) +
     :                    147*u(1,j)**3*u(2,j) - 47*u(1,j)**3*u(2,j) +
     :                    11*u(1,j)*u(4,j) - 15*u(1,j)**5 -2*u(5,j)
             else if (i .eq. 4) then
               u(4,j+1) = 23*u(3,j)**2 -284*u(2,j)*u(3,j)*u(1,j) +
     :                    34*u(2,j)*u(4,j) + 657*(u(2,j)*u(1,j))**2 -
     :                    72*u(2,j)**3 - 888*u(1,j)**4*u(2,j) +
     :                    288*u(1,j)**3*u(3,j) -69*u(1,j)**2*u(4,j) +
     :                    13*u(1,j)*u(5,j) + 300*u(1,j)**6 - 2*u(6,j)

             else if (i .eq. 5) then
cxi		do iii=0,8
cxi	     write(*,'(5d15.5)') (u(iii,jjj), jjj=0,4)
cxi		enddo


               u(5,j+1) = -2*u(7,j) + 3030*u(2,j)*u(3,j)*u(1,j)**2 -
     :                    490*u(2,j)*u(4,j)*u(1,j) -
     :                    500*u(2,j)**2*u(3,j) +47*u(2,j)*u(5,j) -
     :                    2040*u(1,j)**4*u(3,j) +495*u(1,3)**3*u(4,j) -
     :                    95*u(1,j)**2*u(5,j) +15*u(1,j)*u(6,j) -
     :                    1800*u(1,j)**7 +80*u(3,j)*u(4,j)- 
     :                    330*u(3,j)**2*u(1,j) -
     :                    6180*u(2,j)**2*u(1,j)**3 + 
     :                    1530*u(2,j)**3*u(1,j) + 6240*u(1,j)**5*u(2,j)

             else if (i .eq. 6) then
               u(6,j+1) = -2*u(8,j) + 62100*u(2,j)**2*u(1,j)**4 -
     :                    24660*u(2,j)**3*u(1,j)**2 +
     :                    4020*u(3,j)**2*u(1,j)**2 -
     :                    32640*u(2,j)*u(3,j)*u(1,j)**3 +
     :                    5985*u(2,j)*u(4,j)*u(1,j)**2 +
     :                    12150*u(2,j)**2*u(3,j)*u(1,j) -
     :                    1310*u(3,j)*u(4,j)*u(1,j) - 
     :                    774*u(2,j)*u(5,j)*u(1,j) -
     :                    990*u(2,j)**2*u(4,j) -1330*u(2,j)*u(3,j)**2 +
     :                    16440*u(1,j)**5*u(3,j) +17*u(7,j)*u(1,j) +
     :                    127*u(3,j)*u(5,j) +62*u(2,j)*u(6,j) -
     :                    4020*u(1,j)**4*u(1,4) + 780*u(1,j)**3*u(5,j) -
     :                    125*u(1,j)**2*u(6,j) -50040*u(1,j)**6*u(2,j) +
     :                    12600*u(1,j)**8 +80*u(4,j)**2 + 1530*u(2,j)**4
            end if
            u(i,j+1) = (w(i) + u(i,j+1)/8.d0)/u(0,j+1)

  30      continue
20      continue

        if( u(0,3) .le. 0.d0 .or. u(0,4) .le. 0.d0 ) then
                write(*,*) '  in WKB: r-value too small, r=', r
                deltaz=99.d0
                return
        endif
        zeta = sqrt(u(0,4))
        dz = (w(1)+(7*u(1,3)*u(2,3)-5*u(1,3)**3-2*u(3,2))/8)/(2*zeta)
        deltaz = zeta - sqrt(u(0,3))

        return 
        end
*
*     ------------------------------------------------------------------
*           G E T E I J
*     ------------------------------------------------------------------
*
*       This routine adjusts the lagrange multiplier so as to 
*       satisfy the orthogonality condition
*
        subroutine geteij(n0,eij,ov)
        implicit real*8(a-h,o-z)
        common /saveeij/e1,e2,o1,o2,n

        EF(e1,e2,o1,o2)  = (e1*o2-e2*o1)/(o2-o1)
        EFF(e1,e2,eij,k) = (k+1.)*(e1+e2) - (k+k+1.)*eij
        k0=3

        if( n0 .ne. 0 ) goto 5 
        n0=1
        if( n .eq. 2 ) then
           n=1
           eij0 = eij
           if( eij .gt. e2 ) then
                eij = e1
           else if( eij .lt. e1 ) then
                eij = e2
           else
                ediff = e2-e1
                eij = eij + ediff
           endif
           e1=eij0
           o1=ov
           return

        else
           n = 0
        endif
5       if( n .eq. 2 ) goto 10 
        n = n + 1

        if( n .eq. 1 ) then
               e1=eij
               o1=ov
               eij  = - 2.d0* dabs(eij)
               return
        else
               e2 = eij
               o2 = ov
cxi	let e1 < e2
cxi
               call order12(e1,e2,o1,o2)

               if( o1*o2 .lt. 0.d0 ) then
                  eij = EF(e1,e2,o1,o2)
               else
cxi
cxi	if e0 = e1: eij = e2 + k del, del = e2-e1
cxi	if e0 = e2: eij = e1 - k del
cxi	==>: 	    eij = (k+1)(e1+e2) - (2k+1)e
cxi
                  eij =EFF(e1,e2,eij,k0) 
               endif
               return
        endif
cxi	
cxi	already saved 2 eij values
cxi	the current one is the third one
cxi
cxi
cxi	if level 00000000000000000000000000000000000
cxi
10      if( o1*o2 .lt. 0.d0 ) then
cxi	
cxi	o1, o2 opposite sign
cxi
          if( ov*o1 .lt. 0.d0 ) then
            e2 = eij
            o2 = ov
          else if( ov*o2 .lt. 0.d0 ) then
            e1 = eij
            o1 = ov
          else
            stop ' impossible case in geteij '
          endif
          call order12(e1,e2,o1,o2)
          eij = EF(e1,e2,o1,o2) 
          return
cxi
cxi	o1, o2 same sign
cxi
cxi	else level 00000000000000000000000000000000000
cxi
        else
cxi..........................................
          if( ov*o1 .lt. 0.d0  ) then
cxi
cxi	case: e --- e1 --- e2
cxi
            if( eij .lt. e1 ) then
                e2 = e1
                o2 = o1
                e1 = eij
                o1 = ov
            else
cxi
cxi	case: e1 --- e2 --- e
cxi
                e1 = e2
                o1 = o2
                e2 = eij
                o2 = ov
            endif
            eij =  EF(e1,e2,o1,o2)
            return
          else
cxi
cxi	o,o1,o2 same sign
cxi
             if( eij .lt. e1 ) then
                e1 = eij
                o1 = ov
             else
                e2 = eij
                o2 = ov
             endif
             eij = EFF(e1,e2,eij,k0) 
             return
          endif
cxi............................................
cxi
cxi	endif level 00000000000000000000000000000000000
cxi
        endif
        end

        subroutine order12(e1,e2,o1,o2)
        implicit real*8(a-h,o-z)

        if( e1 .gt. e2 ) then
            et = e1
            ot = o1
            e1 = e2
            o1 = o2
            e2 = et
            o2 = ot
        endif
        end
