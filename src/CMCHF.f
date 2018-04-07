*     ------------------------------------------------------------------
*		CONTINUUM PROGRAM
*
*       Written by :  Charlotte Froese Fischer
*                     Department of Computer Science
*                     Vanderbilt University
*
*       1987 Version, Updated October 1993.
*
*       Modified by:  Jinhua Xi, December 1994
*                     for use of WKB method for normalization
*     ------------------------------------------------------------------
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
*     As a result all equations must be transformed as described in
*     Sec. 6-2 and 6-4.
*
*     ------------------------------------------------------------------
*            M A I N   P R O G R A M
*     ------------------------------------------------------------------
*
      PROGRAM CMCHF
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500,IDIM=1000,NCDIM=32768)
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
      INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
      COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :         ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
      LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
      COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
      COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR

*
      LOGICAL PRINT,LD
      CHARACTER*3 ANS*1,NAME(7)*24,ELORT(20,2)
      EQUIVALENCE (IUC,IOU(1)),(OUC,IOU(4))
      DATA NAME/'cfg.inp','int.lst','wfn.inp','cfg.out',
     :           '   ','wfn.out', '   '/
      REAL times(2),rtime,dtime
*
*  ***** Define unit numbers and open files *********************
*								*
*	 UNIT NUMBERS AND FILE NAMES MAY BE MACHINE		*
*	 DEPENDENT. CHECK THE FOLLOWING SECTION. 		*
*								*
*	 IN - Standard input unit, normally the terminal	*
*	 OUT- Standard output unit, normally the terminal	*
*	 PRI- Printer output unit or file.			*
*								*
	IN = 5
	OUT = 6
	PRI = 3
*
*  *****  WRITE OUT HEADER
*
      WRITE(OUT,9)
9     FORMAT(//20X,'======================='/
     :         20X,' UNIX  CMCHF  ...  1987'/
     :         20X,'=======================')
*
*  *****  WRITE OUT DIMENSION INFORMATION
*
      WRITE(OUT,99) 'NCFG',NCFG,'NWF',NWD,'NO',NOD
99    FORMAT(//10X,'THE DIMENSIONS FOR THE CURRENT VERSION ARE:'/
     :       (10X,3(A6,'=',I3,4X)/)/)
*
*  *****  INITIALIZE COMMON DATA ARRAYS
*
      CALL INITA
      CALL INITR
*								*
*  ***** IN THE OPEN STATEMENTS CHECK FOR VALID FILE NAMES ******
*								*
1     WRITE(ERR,'(//A/A//)') ' START OF CASE',' ============='
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
      OPEN(UNIT=PRI,FILE='summry',STATUS='UNKNOWN')
*     .. on systems where 'APPEND' is allowed, this is useful
*     OPEN(UNIT=31,FILE='phase.out',STATUS='UNKNOWN',ACCESS='APPEND')
      OPEN(UNIT=31,FILE='phase.out',STATUS='UNKNOWN')
      write(31,*) '   k^2         sum(c^2)      delta       delta/pi'
*
*  ***** END OF INPUT/OUTPUT INTERFACE **************************
*
*      The following is a non-standard procedure for timing a
*     calculation.  It may be deleted or replaced.
*
*	  RTIME = DTIME(TIMES)
*
      FAIL = .FALSE.
      DO 4 I=1,(30)
      DPM(I) = D10
      IEPTR(I) = 0
     
4     CONTINUE
      DO 5 I = 1,(98)
         IJE(I) = 0
5     CONTINUE
*
*  *****  DETERMINE DATA ABOUT THE PROBLEM
*
	iend= - 1

      CALL CDATA(Etarget,nort,ELORT,iend)

*
*  *****  SET PARAMETERS TO THEIR DEFAULT VALUE
*
13    PRINT = .FALSE.
      SCFTOL = 1.D-6
      NSCF = 20
      IC = 0
      ACFG = D0
      TRACE = .FALSE.
      WRITE(ERR,'(/A)') ' Default values for other parameters ? (Y/N) '
      READ (IN,'(A)') ANS
      IF (ANS .EQ. 'N' .OR. ANS .EQ. 'n') THEN
         WRITE(ERR,'(/A,A)') ' Default values for',
     :   ' PRINT, SCFTOL ? (Y/N) '
         READ(IN,'(A)') ANS
         IF ( ANS .NE. 'Y' .AND. ANS .NE. 'y'  ) THEN
            WRITE(ERR,'(A)') ' Input free FORMAT(L, F) '
            READ(IN,*) PRINT, SCFTOL
         END IF
         WRITE(ERR,'(/A)') ' Default values for NSCF, IC ? (Y/N) '
         READ(IN,'(A)') ANS
         IF (ANS .NE. 'Y' .AND. ANS .NE. 'y' ) THEN
            WRITE(ERR,*) ' Input free FORMAT(I, I) '
            READ(IN,*) NSCF, IC
         END IF
         WRITE(ERR,'(/A)') ' Default values for ACFG,TRACE ? (Y/N) '
         READ(IN,'(A)') ANS
         IF (ANS .NE. 'Y' .AND. ANS .NE. 'y') THEN
             WRITE(ERR,'(A)') ' Input free FORMAT( F, L) '
             READ(IN,*) ACFG,TRACE
         END IF
      END IF
      WRITE(OUT,2) PRINT,CFGTOL,SCFTOL,NSCF,IC,ACFG,ID,TRACE
2     FORMAT(/L3,2D6.1,2I3,F3.1,I3,L3)
*
*
11	ekk = e(nwf,nwf)
	if( ekk .ge. d0 ) goto 7
      WRITE(OUT,'(//A,F15.7/)') ' *** CALCULATIONS FOR k^2 =',-EKK
      WRITE(ERR,'(//A,F15.7/)') ' *** CALCULATIONS FOR k^2 =',-EKK

	call  CSCF(Etarget,EKK,ACFG,SCFTOL,PRINT,nort,elort,iend)
        CALL CDATA(Etarget,nort,ELORT,iend)

	if( iend .le. 0 ) goto 11
22	WRITE(ERR,'(//A)') ' continue for other energies? '
         READ(IN,'(A)') ANS
         IF (ANS .ne. 'Y' .and. ANS .ne. 'y') goto 6 

*     Determine the range of values for K*K
*
7     WRITE(ERR,'(//A)') ' Enter k^2: MIN, DELTA, MAX '
      READ(5,*) ELOW, EDELTA, EHIGH
      EK = ELOW
100   EKK = -EK
      WRITE(OUT,'(//A,F10.4/)') ' *** CALCULATIONS FOR k^2 =',-EKK
      CALL EIJSET(NWF,NWF,EKK)
*
*  *****  PERFORM THE MCHF ITERATION
*
      CALL CSCF(Etarget,EKK,ACFG,SCFTOL,PRINT,nort,elort,iend)
      EK = EK + EDELTA
      IF (EDELTA .GT. D0) THEN
         IF (EK .LE. EHIGH + 0.5*EDELTA) GO TO 100
      ELSE
         IF (EK .GE. EHIGH + 0.5*EDELTA) GO TO 100
      END IF
      goto 22
*
*  *****  DETERMINE END OF CASE
*
6     CONTINUE
*    
*     ON systems where dtime is implemented, these are useful
*     RTIME = DTIME(TIMES)
*     WRITE(ERR,'(//A/A//A/3F12.3//)') ' END OF CASE',' ===========',
*    :  '  Real        User        System  Time (in minutes)',
*    : RTIME/60.,TIMES(1)/60., TIMES(2)/60.
      END
*
*     ------------------------------------------------------------------
*             C D A T A
*     ------------------------------------------------------------------
*
*       Data concerning the number of configurations (NCFG), the number
*   and type of electrons in each  configuration,  as  well  as  data
*   associated with the energy expression are read and stored.
*
*
      SUBROUTINE CDATA(Etarget,nort,ELORT,nend)
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
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,30),YK(NOD),YR(NOD)
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
      COMMON /zzind/ZZ(NWD),IND(NwD),IELI(5),NOCCSH(NCD)
*
      LOGICAL SETORT,STRONG
      CHARACTER*3 ELORT(20,2)
      CHARACTER*3 EL1,EL2,ELCLSD(18),ELI(5),ANS*1,STRING*40,LIST*72
*
    1 FORMAT(18(1X,A3))
    7 FORMAT(A3,F6.0,I3,I3,F3.1)


      if( nend .eq. -1 ) goto 5
cxi
cxi	read configuration weight for another continuum orbital energy
cxi

	nend = 0
	read(IUC,'(A)', end=99) string
	read(IUC,'(A)', end=99) string
	do 77 nc = 1, ncfg
	read(IUC,'(A40,F10.8)', end=99) string, WT(NC)
	read(IUC,'(A)', end=99) string
77	continue
	read(IUC,'(A)', end=99) string
	goto 28

*
*  *****  First time to enter CDATA,  READ 'ATOM' CARD
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
      READ(STRING(I+J+1:),*) Z
*
*
*  *****  READ CONFIGURATIONS FROM A FILE
*
      READ(IUC,'(15X,F14.7/18(1X,A3))') Etarget,(ELCLSD(I),I=1,18)
      NCFG = 0
 3    READ(IUC,'(A40,F10.8)',END=10) STRING,W
      IF (STRING(1:1) .NE. '*' .AND. STRING(1:3) .NE. '   ') THEN
         NCFG = NCFG+1
	if( w .eq. 0.d0 ) w = 1.d0
         IF (NCFG .LE. (NCD) ) THEN
            CONFIG(NCFG) = STRING
            WT(NCFG) = W
            READ(IUC,'(9(5X,A3))') (COUPLE(NCFG,J),J=1,9)
            GO TO 3
           ELSE
            STOP ' TOO MANY CONFIGURATIONS: MAX =(NCD)'
         END IF
      END IF
10    CONTINUE
	  ID = -1
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
         STOP ' TOO MANY CLOSED SHELLS: MAX = 18'
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
                  IF (K .GT. (30)) STOP ' TOO MANY ELECTRONS: MAX= (30)'
                  GO TO 17
               END IF
              ELSE
*
*  *****         A new electron has been found; add it to the list
*
               MAXORB = K
               EL(MAXORB) = EL1
               IF ((EL1(1:1) .EQ. 'k' .OR. EL1(2:2) .EQ. 'k' .OR.
     :            EL1(1:1) .EQ. 'n' .OR. EL1(2:2) .EQ. 'n') ) then
               		korb = k
			if( ID .EQ. -1) ID = NC-1
		endif

            END IF
            J = J+8
            IF (J .LT. 40) GO TO 16
         END IF
         NOCCSH(NC) = I

   15 CONTINUE
	  IF (ID .EQ. -1) THEN
	    WRITE(ERR,*) ' STOP in DATA: No continuume function found'
	     CALL EXIT(1)
	  END IF
*
*  *****  The list of electrons has been determined
*
      NWF = MAXORB
cxi
cxi	move the continuum orbital to the last
cxi
	if( korb .ne. nwf ) then
	  el1 = el(korb)
	  do 82 j=korb,nwf-1
	  el(j) = el(j+1)
82	  continue
	  el(nwf) = el1
	endif

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
         READ(LIST(J+1:),*) NIT
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
      IF ( ANS.NE.'Y' .AND. ANS.NE.'y' .AND. NIT.NE.0) then 
	WRITE(ERR,*)
     :   ' S, IND, METH, ACC for electrons to be varied (free-format)'
      ENDIF

      DO 20 I = NCLOSD+1,NWF
     	 IF (ANS.EQ.'Y' .OR. ANS.EQ.'y' .OR. I.LT.IB) THEN
	   S(I) = SS
	   METH(I) = 3
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
	 IF (EL(I)(J-1:J-1) .EQ. 'n') METH(I) = 5
 20   CONTINUE
*
*  *****  CHECK METHOD AND SET ORTHOGONALITY
*
28    EL1 = EL(NWF)
      IF (.NOT. (EL1(1:1).EQ.'k' .OR. EL1(2:2).EQ.'k'))
     :     STOP ' Last orbital not a continuum orbital'
      IF (METH(NWF) .NE. 4) THEN
          METH(NWF) = 4
*         IND(NWF) = 1
          DO 95 J = 1,NO
             P(J,NWF) = D0
95        CONTINUE
          AZ(NWF) = D1
          MAX(NWF) = NO
      END IF
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
		IF (I.LE.NCLOSD .AND. J.LE.NCLOSD) C = 1.D-10
	        CALL EIJSET(I,J,C)
		CALL EIJSET(J,I,C)
            END IF
   39    CONTINUE
   38  CONTINUE
*
*  *****  DETERMINE ADDITIONAL ORTHOGONALITY PAIRS
*
      if( nend .eq. 0 ) then
	do 87 i = 1, NORT
	    EL1 = ELORT(I,1)
    	    EL2 = ELORT(I,2)
	    CALL EPTR(EL,EL1,I1,*99)
            CALL EPTR(EL,EL2,I2,*99)
	    CALL EIJSET(I1,I2,1.D-5)
	    CALL EIJSET(I2,I1,1.D-5)
87	continue
	goto 66
      endif

      I = 0 
      IF ( IUC .NE. IN) THEN
   40    READ(IUC,1,END=50) EL1,EL2

         IF ( EL1(1:1) .NE. '*' .AND. EL2 .NE. '   ') THEN
            I = I +1
            IF (I .GT. (20)) STOP ' TOO MANY ORTHOGONALITIES: MAX=(20)'

            ELORT(I,1) = EL1
            ELORT(I,2) = EL2

            CALL EPTR(EL,EL1,I1,*99)
            CALL EPTR(EL,EL2,I2,*99)
	    CALL EIJSET(I1,I2,1.D-5)
	    CALL EIJSET(I2,I1,1.D-5)

            GO TO 40
         END IF
      END IF
   50 CONTINUE
      NORT = I
*
*  *****  ADDITIONAL PARAMETERS
*
      WRITE(ERR,'(/,A)') ' Default values (NO,REL,STRONG) ? (Y/N) '
      READ(IN,'(A)') ANS
      IF (ANS .EQ. 'Y' .OR. ANS .EQ. 'y') THEN
         NO = (220)
         REL = .FALSE.
         STRONG = .FALSE.
         IF (NCFG .GT. 1) STRONG = .TRUE.
       ELSE
         WRITE(ERR,*) ' Enter values in FORMAT(I,L,L) '
         READ(IN,*) NO, REL, STRONG
         IF (NO .GT. (220)) STOP
     :     ' TOO MANY POINTS FOR EACH FUNCTION: MAX=(220)'
      END IF
      ND = NO - 2
      WRITE(OUT,61) ATOM,TERM,Z,NO,NWF,NIT,NCFG,REL,STRONG
61    FORMAT(/1X,2A6,F6.0,I6,3I3,2L3)
      WRITE(PRI,62) ATOM,TERM,Z,(EL(I),4*L(I)+2,I=1,NCLOSD)
62    FORMAT(1H1///9X,33HHARTREE-FOCK WAVE FUNCTIONS FOR  ,2A6,4H Z =,
     1   F5.1//14X,'CORE = ',5(A3,'(',I2,')'))

      OMIT = .NOT. STRONG
66    WRITE(PRI,'(//11X,A,37X,A//)') 'CONFIGURATION','WEIGHT'
*
* *****  WRITE 'CONFIGURATION' CARDS  AND CHECK THE WEIGHTS
*
      W = D0
      DO 63 I=1,NCFG
   63 W = W + WT(I)**2
      IF (W .EQ. D0) STOP ' WEIGHT information omitted'

      DO 68 I = 1, NCFG
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
	if( nend .eq. 0 ) then
cxi	
cxi	keep core and target orbitals unchanged
cxi
	   do 83 i = 1, nwf 
	     if( i .lt. ib ) then
         	   ind(I) = 1
	     else
		   ind(i) = - 1
	     endif
83	  continue
	endif

      CALL WAVEFN(nend)

      if( nend .ne. -1 ) return

      DO 100 I=1,6
	 INTPTR(I) = 0
100   CONTINUE
*
      IF (IUD .NE. IN) CALL INTGRL
*
*  *****  DEFINE SUM(I)
*
	IBEGIN = INTPTR(5)+1
	IEND = INTPTR(6)
	DO 80 I = IBEGIN,IEND
	   IF (IEL(I,1).EQ.IEL(I,2)) SUM(IEL(I,1)) = -2*COEF(I)
 80	CONTINUE

*
      RETURN
99    if( nend .eq. -1 ) then
	STOP
      else
	nend = 1
	return
      endif
      END
*
*     ------------------------------------------------------------------
*            C D E
*     ------------------------------------------------------------------
*
*       This routine controls the solution of the differenttial equation
*   for the radial function P  .  One of three methods is selected -
*                            I1
*   M1, M2, or M3 -  for solving the equations,  the  initial  choice
*   being determined by an input paramter, METH(I), except when no
*   exchange is present, in which case M2 is selected. (For further
*   information see Sec. 7-4)
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
      SUBROUTINE CDE(I1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500,IDIM=1000,NCDIM=32768)
*
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
      COMMON P2(NOD),HQ(NOD),XX(NOD),AC(20,20),BC(20),JV(20),
     :     AZZ,PP,FN,EM,FM,EU,FU,DELTAE,M,NODE,MK,KK,NJ
*
      LOGICAL DIAG
      CHARACTER*2 ASTER(5)
      DATA ASTER/'  ','* ','**',' c',' b'/
      DOUBLE PRECISION ED1
      SAVE ED1
      DATA ED1 /0/
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
      IF (METH(I) .LE. 3) THEN
         DO 9 J = 1,M
9        PDE(J) = PDE(J)/PN
         AZZ = AZZ/PN
      END IF
*
*  *****  CHECK IF METHOD 2 SHOULD BE USED
*
      IF ( KK .NE. 1 ) GO TO 13
      IF (DABS(D1 -ED2/E(I,I)) .GT. 0.01D0  .OR.
     1    DMAX1(DABS(D1 - PN), DABS(D1/PN - D1)) .LT. 0.10D0 )
     2    GO TO 13
       METH(I) = 2
       KK = 2
       GO TO 25
*
*  *****  SET THE ACCELERATING PARAMETER
*
*
13    IF (IPR .NE. I .OR. KK.GT.3) GO TO 33
      ED2 = ED2 - E(I,I)
      IF (ED1*ED2 .GT. D0) ACC(I) = .75*ACC(I)
      IF (ED1*ED2 .LT. D0) ACC(I) = (D1 + D3*ACC(I))/D4
33    C = ACC(I)
      CD = D1 - C
      DP = DP/PN
*
*   *****  IMPROVE THE ESTIMATES
*
      AZ(I) = CD*AZZ + C*AZ(I)
      AZZ = AZ(I)
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
*
*  *****  CHECK THE ORTHOGONALIZATION
* 
26    NN = NWF
      IF (OMIT) NN = IB - 1
      IBEGIN = 1
      IF (I .GT. 1) IBEGIN = IEPTR(I-1) + 1
      IP = IBEGIN
      IJ = 0
50    JI = IJE(IP)
      IF (JI .NE. I .AND. JI .NE. NWF) THEN
      IF (JI .GE. IB .AND. DPM(JI) .GE. DPM(I)) THEN
*
*       The JI orbital should be orthogonalized
*
         C = QUADR(I,JI,0)
         MM = MIN0(MAX(JI),MAX(I))
         DO 51 J = 1,MM
            P(J,JI) = P(J,JI) - C*P(J,I)
51       CONTINUE
         AZ(JI) = AZ(JI) - C*AZ(I)
         IF ( METH(JI) .LE. 3) THEN
            C2 = SQRT(QUADR(JI,JI,0))
            DO 52 J = 1,MM
               P(J,JI) = P(J,JI)/C2
52          CONTINUE
            AZ(JI) = AZ(JI)/C2
         END IF
         VARIED(JI) = .TRUE.
         MAX(JI) = MM
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
         IF ( .NOT. DIAG .AND. IJ .GT. 1) CALL LINEQN(20,IJ,AC,BC)
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
         IF (METH(I) .LE. 3) THEN
            PNN = DSQRT(QUADR(I,I,0))
            DO 66 J = 1,M
               P(J,I) = P(J,I)/PNN
66          CONTINUE
            AZZ = AZZ/PNN
         END IF
      END IF
      M = NO
67    IF (DABS(P(M,I)) .LT. 1.D-15) THEN
         P(M,I) = D0
         M = M-1
         GO TO 67
      END IF
      MAX(I) = M
*     IF (AZZ .GT. D0) AZ(I) = DMAX1(AZZ,D5*AZ(I))
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
      CALL CUPDATE
      CALL CGRANGE
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
*
*     ------------------------------------------------------------------
*    3-5       D I A G
*     ------------------------------------------------------------------
*
*       The CDIAG subroutine computes an energy matrix and if the number
*   of configurations is greater than  1,  finds  an  eigenvalue  and
*   eigenvector of this matrix.
*
*       Given the accelerating parameter for configuration mixing,
*   ACFG, the  new mixing coefficients are stored in COMMON.
*
*
      SUBROUTINE CDIAG(ETOTAL,ACFG,PCONV)
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
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
      LOGICAL PCONV
      COMMON WP(NCD),W(NCD,NCD)
*
      DO 1 I = 1,NCFG
         DO 2 J = 1,NCFG
            W(I,J) = D0
  2        CONTINUE
          WP(I) = WT(I)
         W(I,I) = EC-ETOTAL
  1     CONTINUE
*
      IBEGIN = 1
      IEND = INTPTR(6)
      J = 0
      DO 10 I = IBEGIN,IEND
 11      IF (CPTR(I) .GT. J) THEN
            J = J + 1
            C = COEFF(J)*VALUE(I)
            IF (OPTR(J) .NE. 0) C = C*VALUE(OPTR(J))
            W(IH(J),JH(J)) = W(IH(J),JH(J)) + C
            GO TO 11
         END IF
 10   CONTINUE
*
*  ***** SYMMETRIZE THE MATRIX
*
      DO 12 I = 1,NCFG-1
         DO 13 J = I+1,NCFG
            W(I,J) = W(J,I)
 13      CONTINUE
 12   CONTINUE
      IF (TRACE) THEN
          WRITE(OUT,'(/10X,A,F16.8,10X,A,F16.8/)')
     :          'EC =',EC,'ETOTAL =',ETOTAL
          DO 15 I=1,ID
        WRITE(OUT,'(I4,6F12.7/(4X,6F12.7))')
     :                 I,(W(I,J),J=1,NCFG)
15       CONTINUE
      END IF
*
14    IF (NCFG .EQ. 1) GO TO 37
*
*
*     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
*
      DO 35 K = 1,ID-1
         KP1= K+1
*
*        FIND PIVOT
*
         M = K
         DO 21 I = KP1,ID
            IF (DABS(W(I,K)) .GT. DABS(W(M,K))) M = I
   21    CONTINUE
         T = W(M,K)
         W(M,K) = W(K,K)
         W(K,K) = T
*
*        SKIP STEP IF PIVOT IS ZERO
*
         IF (T .EQ. 0.0D0) GO TO 35
*
*        COMPUTE MULTIPLIERS
*
         DO 20 I = KP1,ID
             W(I,K) = -W(I,K)/T
   20    CONTINUE
*
*        INTERCHANGE AND ELIMINATE BY COLUMNS
*
         DO 30 J = KP1,NCFG
             T = W(M,J)
             W(M,J) = W(K,J)
             W(K,J) = T
             IF (T .EQ. 0.0D0) GO TO 30
             DO 23 I = KP1,ID
                W(I,J) = W(I,J) + W(I,K)*T
   23        CONTINUE
   30    CONTINUE
   35 CONTINUE
      DO 24 J = ID,1,-1
      JP = J+1
      WT(J) = D0
      DO 25 K = JP,NCFG
25    WT(J) = WT(J) - W(J,K)*WT(K)
24    WT(J) = WT(J)/W(J,J)
      DO 28 I = 1,ID
      WT(I) = WT(I) + ACFG*(WP(I) - WT(I))
28    CONTINUE
*
 37   WRITE(OUT,636) (I,WT(I),I=1,ID)
636   FORMAT(/(5(I4,F12.6)))
*
*  *****  REDEFINE SUM(I)
*
      IBEGIN = INTPTR(5)+1
      IEND = INTPTR(6)
      DO 50 I = IBEGIN,IEND
         IF (IEL(I,1).EQ.IEL(I,2)) SUM(IEL(I,1)) = -2*COEF(I)
 50   CONTINUE
      END
*
*     ------------------------------------------------------------------
*             C G R A N G E
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
*   When  only one pair with a given set of principal quantum numbers
*   is present, ELAGR(I,J) is used to  determine  the  off-  diagonal
*   energy  parameters  as  long  as  |q  -q | > 0.05.  Otherwise Eq.
*                                       i   j
*   (7-10) modified for configuration interaction is used.
*
*
      SUBROUTINE CGRANGE
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
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,30),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
      INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
      COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :         ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
      COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
      COMMON W(NCD,NCD),U(NCD),DWT(NCD),AC(20,20),BC(20),JV(20),IV(20)
      LOGICAL DIAG, FIRST, REL
*
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
         IF ( IJ .GT. 20) STOP '  TOO MANY LAGRANGE MULTIPLIERS'
                 IV(IJ) = I
         JV(IJ) = J
            END IF
12         CONTINUE
11       CONTINUE
*
*  ***** IJ IS THE NUMBER OF LAGRANGE MULTIPLIERS FOR l = IL
*
           IF (IJ .EQ. 0) GO TO 10
         DO 13 II = 1,IJ
            BC(II) = D0
            DO 14 III = 1,IJ
         AC(II,III) = D0
14          CONTINUE
13       CONTINUE
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
20                CONTINUE
                  FIRST = .FALSE.
               END IF
               DO 22 JJ = 1,NO
                  YR(JJ) = P(JJ,J)
22             CONTINUE
               BC(II) = BC(II) +
     :                    HL(EL,I,J,REL)-D2*QUADS(I,J,1)-QUAD(J,NO,YR,X)
            END IF
18      CONTINUE
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
             DIAG = .FALSE.
          ELSE IF (JV(II) .EQ. JV(III) .AND. JV(II) .GE. IB
     :              .AND. E(IV(II),IV(III)) .EQ. D0) THEN
             AC(II,III) = QUADR(IV(II),IV(III),0)/SUM(JV(II))
             AC(III,II) = AC(II,III)
             DIAG = .FALSE.
          END IF
26          CONTINUE
24       CONTINUE
         IF ( .NOT. DIAG ) CALL LINEQN(20,IJ,AC,BC)
         DO 28 II = 1,IJ
              CALL EIJSET(IV(II),JV(II),BC(II)/SUM(IV(II)))
            IF ( JV(II) .GE. IB )
     :            CALL EIJSET(JV(II),IV(II),BC(II)/SUM(JV(II)))
28       CONTINUE
10    CONTINUE
*
*  *****  PRINT THE OFF-DIAGONAL ENERGY PARAMETERS
*
      DO 30 I = IB,NWF
         DO 32 J = 1,I-1
            IF (DABS(E(I,J)) .GT. 1.D-10) THEN
         WRITE(OUT,35) EL(I),EL(J),E(I,J),EL(J),EL(I),E(J,I)
35               FORMAT(7X,2(3X,'E(',2A3,') =',1PE12.5))
            END IF
32       CONTINUE
30    CONTINUE
      RETURN
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
*          C O U T P U T
*     ------------------------------------------------------------------
*
*       The radial functions and orthogonality integrals are printed,
*   if PRINT is .TRUE.   The  functions  will  also  be  punched  (or
*   stored) on unit OUF, if OUF .NE. 0.
*
*
      SUBROUTINE COUTPUT(PRINT,DELTA,Nstart)
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
      LOGICAL PRINT, REL
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
cxi      DO 3 I = NCLOSD+1,NWF
      DO 3 I = Nstart,NWF
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
*
*     ------------------------------------------------------------------
*    3-32     C S C F
*     -----------------------------------------------------------------
*
*       This routine controls the SCF procedure described in Chapter
*   7.  If certain input parameters are zero (or blank) they will  be
*   set to their default value.
*
*          Parameter       Default Value
*          --------        -------------
*          SCFTOL          1.D-5
*          IC              (NWF + 1 - IB)/4 + 3
*          NSCF            20
*
*   The self-consistency convergence criterion is
*
*          Z2 = SCFTOL
*
*   It is increased by a factor two at the end of each iteration
*
*
      SUBROUTINE CSCF(Etarget,EKK,ACFG,SCFTOL,print,nort,elort,nend)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500,IDIM=1000,NCDIM=32768)
*
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
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,30),YK(NOD),YR(NOD)
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
*
      LOGICAL PCONV,PRINT
       CHARACTER*3 ELORT(20,2)
*
*  *****  SET THE SCF CONVERGENCE PARAMETER TO AN OPTIMISTIC VALUE
*


cxi================================================================
cxi	Nstart is the first orbital position to be written to wfn.out.
cxi	We need only to write the target orbitals once for the 
cxi	first continuum orbitel energy.
cxi	BUT if so, we need to modify PHOTO, PHOTO reads all target
cxi	and channel orbitals for each energy
cxi
cxi	IF you want to use wfn.out as wfn.inp again, then you need
cxi	to set Nstart = 1 when nend = - 1
cxi
	Nstart = IB 
	if( nend .eq. -1 )  Nstart = 1
cxi
cxi===============================================================
cxi	
cxi	Nstart = NCLOSD + 1    !!! PHOTO requires so

	PI = dacos(-1.d0)
      ETOTAL = Etarget - D5*EKK
      TOL = DSQRT(Z)*1.D-5
      Z2 = SCFTOL
      WRITE(OUT,15)
15    FORMAT(//)
      WRITE(OUT,16) OMIT,ACFG,Z2,NO,REL
   16 FORMAT(10X,44HWEAK ORTHOGONALIZATION DURING THE SCF CYCLE=,L4/
     :       10X,44HACCELERATING PARAMETER FOR MCHF ITERATION  =,F5.2/
     :       10X,44HSCF CONVERGENCE TOLERANCE (FUNCTIONS)      =,1PD9.2
     :      /10X,44HNUMBER OF POINTS IN THE MAXIMUM RANGE      =,I4/
     :       10X,44HRELATIVISTIC DIAGONAL  ENERGY CORRECTIONS  =,L4//)
*
*  *****  SET ITERATION PARAMETERS
*

      ICYCLE = 0
      IPR = 0
      CALL COULOM(NWF,EKK)
      if (P(1,nwf) .eq. D0) then
        CALL CSOLVE(NWF,DELTA)
        call orthog
      end if
      PCONV = DPM(NWF) .LE. Z2
      CALL CUPDATE
      IF (ID .GE. 1)  CALL CDIAG(ETOTAL,ACFG,PCONV)
      IF ( IB .GT. NWF ) GO TO 17
*
*  *****  PERFORM NSCF SELF-CONSISTENT FIELD ITERATIONS
*
9     DO 100 I = 1,NSCF

      ICYCLE = ICYCLE + 1
      WRITE(OUT,7) ICYCLE,Z2
7     FORMAT(//10X,17HITERATION NUMBER ,I2/10X,16H----------------/
     1 10X,'CONVERGENCE CRITERIA =',1PD9.1/)
      IF (IB .GT. NWF) GO TO 17
        CALL CGRANGE
*
*  *****  SOLVE EACH DIFFERENTIAL EQUATION IN TURN
*
      WRITE(OUT,14)
14    FORMAT(/20X,' EL',6X,'ED/DELTA',10X,'AZ',11X,'NORM',7X,'DPM')
      CALL CSOLVE(NWF,DELTA)
      DP1 = DPM(NWF)
      JJ = NWF
      DO 2 J = IB,NWF-1
         CALL CDE(J)
         IF ( FAIL ) RETURN
         DP = DPM(J)*DSQRT(SUM(J))
         IF ( DP1 .GE. DP ) GO TO 2
         DP1 = DP
         JJ = J
2     CONTINUE
      IF (((NCFG .EQ. 1 .OR. ID .EQ. 0) .AND. DP1 .LT. Z2) .OR.
     :      IC .LE. 0) GO TO 6
*
*  *****  SOLVE IC DIFFERENTIAL EQUATIONS EACH TIME SELECTING THE
*  *****  ONE WITH THE LARGEST DPM
*
      DO 4 II =1,IC
      IF (JJ .EQ. NWF) THEN
        CALL CSOLVE(NWF,DELTA)
        FAIL = .FALSE.
      ELSE
        CALL CDE(JJ)
      END IF
      IF ( FAIL ) RETURN
      DP1 = DPM(NWF)
      JJ = NWF
      DO 5 J = IB,NWF-1
         DP = DSQRT(SUM(J))*DPM(J)
         IF ( DP .LT. DP1 ) THEN
            JJ = J
            DP1 = DP
         END IF
5     CONTINUE
      IF (DP1 .LT. Z2) GO TO 6
4     CONTINUE
6     CALL ORTHOG
      CALL CUPDATE
      PCONV = DPM(NWF) .LE. Z2 .and. DP1 .LE. Z2
12    IF (.NOT. (NCFG .EQ. 1 .OR. ID.EQ.0)) 
     :   CALL CDIAG(ETOTAL,ACFG,PCONV)
      IF (PCONV) GO TO 17
*
*  *****  IF FUNCTIONS APPEAR TO HAVE CONVERGED,SOLVE EACH AGAIN, IN
*  *****  TURN, AND TEST AGAIN
*
1     CONTINUE
      WRITE(OUT,8) EL(JJ),DP1
8     FORMAT(/ 6X,34HLEAST SELF-CONSISTENT FUNCTION IS ,A3,
     1   27H :WEIGHTED MAXIMUM CHANGE =,1PD10.2)
100   Z2=1.3*Z2
*
*  *****  OUTPUT FINAL CALCULATIONS
*
17    CONTINUE
          WRITE(PRI,'(/10X,A,F16.8,10X,A,F16.8/)')
     :          'EC =',EC,'ETOTAL =',ETOTAL
*
*  *****  PUNCH CONFIGURATIONS AND WEIGHTS ON UNIT OUC
*
         WRITE(3,'(//A/)') '     Final Mixing'
         WRITE(OUC,46) ATOM,TERM,Etarget,ETOTAL,delta,-EKK
46       FORMAT(3X,2A6,4F14.7)
         WRITE(OUC,'(18(1X,A3))') (EL(J),J=1,NCLOSD)
cxi	
cxi	get the sum of wt^2
cxi
	 wtol=0.d0
         DO 47 J = 1,NCFG

	 if( j .le. id ) wtol = wtol + wt(j)*wt(j)
cxi==========================================================
cxi	these two lines are used to confine the data for output
cxi	you can comment out these two lines, and modify the 
cxi	format statement 648 .
cxi
	if( wt(j) .gt. 999.999999 ) wt(j) = 999.999999
	if( wt(j) .lt. -99.999999 ) wt(j) = -99.999999
cxi===========================================================

         WRITE(3,648) J,CONFIG(J),WT(J),(COUPLE(J,JJ),JJ=1,9)
47       WRITE(OUC,48) CONFIG(J),WT(J),(COUPLE(J,JJ),JJ=1,9)
48       FORMAT(A40,F10.6/9(5X,A3))
648      FORMAT(I4,2X,A40,F10.6/(6X,9(5X,A3)))
         WRITE (OUC,'(A)') '****'


      if( nend .eq. -1 ) then
	nend = 0
cxi
cxi==========================================================
cxi
cxi	if you include the following lines, the output file
cxi	cfg.out can be used as input file cfg.inp again, 
cxi
cxi	The PHOTO program has been modified to fit this 
cxi	format requirement
cxi
	do 55 j=1,nort
	write(OUC,'(1x,A3,1x,A3)' ) ELORT(j,1),ELORT(J,2)
55	continue	
         WRITE (OUC,'(A)') '****'
cxi==========================================================
cxi
      endif

      CALL CSUMMRY(DELTA)
      CALL COUTPUT(PRINT,DELTA,Nstart)
      NIT = NWF - IB + 1

      WRITE(PRI, 105) NIT, DP1
105   FORMAT(//10X,'NUMBER OF FUNCTIONS ITERATED          =',I6/
     1         10X,'MAXIMUM WEIGHTED CHANGE IN FUNCTIONS  =',D10.2)
cxi
cxi	save useful informations seperately, in file phase.out
cxi	unit = 31
cxi
	write(31,'(F10.6,F13.7,2F12.6)') 
     :           -sngl(ekk),sngl(wtol),sngl(delta),sngl(delta/pi)

      RETURN
      END
*     Jinhua's version without adjusting lagrange multiplier
*     ------------------------------------------------------------------
*       C S O L V E
*     ------------------------------------------------------------------
*
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
*	     regular, FC, and irregular, GC, Coulomb functions, as
*	     obtained from RCWFN. Renormalizes the continuumfunction.
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


      IBEGIN = 1
      IF (I .GT. 1) IBEGIN = IEPTR(I-1)+1
      IP = IBEGIN

      PI = ACOS(-D1)

      ED = E(I,I)
      ekk = -ed
      eK = SQRT(ekk)

      AZD = AZ(I)
      DELTA0 = DELTA
      IPR = NWF
      M = MAX(NWF)
*
* ***** Computes the Exchange function.
*
        CALL XCH(I,3)

*
* ***** Redefine the outer region.
*

678      IF (IX .EQ. 0) THEN
5       IF (eK*(R(M)-R(M-1)) .GT. 1.0D0) THEN
          M = M-1
CORRECTION: ??
*         IF (M .GT. NJ) GO TO 5
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
        WRITE(OUT,'(/1X,A,I4)')
     :       'WARNING: Outer region may be truncated. M = MJ =',M
	WRITE(OUT,'(A)') 
     :       '          Exchange function not small at M!'
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
      CALL CNMRV(NJ,M,1,AZ(I),FK,X,FH,XH,PH,PDE)

       CALL FGCOUL(I,M,CN,DELTA)

cxi
cxi	note: the sign of the overlop depends on cn, 
cxi	because we need to adjust the off-diagonal parameters
cxi	repeatedly to find a zero overlap, we can not change the sign of
cxi	the wavefunction
cxi	so, save the sign in cnn
cxi
       if( cn .lt. 0.d0 ) then
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
50    JI = IJE(IP)
      IF ( JI .NE. I) THEN
        CC = QUADR(I,JI,0)

cxi
cxi the first orthogonal requirements, use new approach
cxi
	if( e(i,ji) .eq. 0.d0  .or. dabs(cc) .lt. 1.d-7) goto 65 

        AZ(I) = AZ(I) - CC*AZ(JI)
        DO 51 J = 1,M
          P(J,I) = P(J,I) - CC*P(J,JI)
51      CONTINUE

65        WRITE(OUT,63) EL(JI),EL(I),CC
63      FORMAT(6X,'<',A3,'|',A3,'>=',1PD8.1 )
c
cxi	re-normalize the wavefunction
cxi
	cnn= cnn*pde(M)/P(M,i)
	Az(i) = cnn*Az(i)
	do 52 j = 1,m
	    p(j,i) = cnn*p(j,i)
52      continue
cxi
cxi	normalize the off diagonal parameters
cxi
        IP = IP+1
        IF (IP .LE. IEPTR(I)) GO TO 50
      END IF

*
      VARIED(I) = .TRUE.
      DP = ABS((delta0-delta)/delta)

      DPM(I) = DP
      WRITE(OUT,17) EL(I),DELTA,AZ(I),CN,'c',DP
17    FORMAT(20X,A3,2F15.7,F12.7,A2,1PD10.2)
      return	
      END


C
C     ------------------------------------------------------------------
C    3-35      C S U M M R Y
C     ------------------------------------------------------------------
C
C       The results of a calculation are summarized.   These include
C   the following for each electron:
C
C          E(NL)   - diagonal energy parameter
C          AZ(NL)  - starting parameter, P(r)/r**(l+1) as r -> 0.
C          SIGMA   - screening parameter as defined by Eq. (6-  ).
C          1/R**3  - expected value of <1/r**3>
C          1/R     - expected value of <1/r>
C          R       - expected mean radius
C          R**2    - expected value of <r**2>
C          I(NL)   - -(1/2)<nl|L|nl>
C          KE      - I(NL) + Z <r>
C          REL     - Relativistic shift (mass-velocity, Darwin term,
C                    spin-spin contact term)
C
C   These results are followed by:
C
C          TOTAL ENERGY--RELATIVISTIC OR NON-RELATIVISTIC (ET)
C          KINETIC ENERGY-- NON-RELATIVISTIC (EN)
C          POTENTIAL ENERGY (EP) = ET - EN
C          RATIO                 = - EP/EN
C                      k   k   k
C   The values of all F , G , R  and <nl|L|n'l> integrals which enter
C   into the calculation are printed, but only if OUD > 0.
C
C
      SUBROUTINE CSUMMRY(DELTA)
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
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,30),YK(NOD),YR(NOD)
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
*
      COMMON R3(NWD),SS(3),R1,RM,RMM,RH,R0,QI,QJ,SP,C,CC,
     1   EKINP,EN,EPOT,RATIO,LI,LJ,K,KF,I1,I2,J1,J2,I,J,MIN
      CHARACTER*1 SYMBOL
*
      WRITE(PRI,9) ATOM,TERM
9     FORMAT(/// 24X,5HATOM ,A6,3X,5HTERM ,A6//45X,13HMEAN VALUE OF,
     1  /3X,2HNL,9X,5HE(NL),9X,6HAZ(NL),
     2   4X,5H R0  ,6X,3H1/R, 9X,1HR, 6X,
     3   11HI(NL)/DELTA,5X,3HREL)
      EN = D0
C
C  *****  COMPUTE AND PRINT ONE-ELECTRON PARAMETERS
C
      DO 10 I = 1,NWF
      R0 = D1
      R1 = D0
      RM = D0
      EK = DELTA
      S(I) = D0
      IF (METH(I) .NE. 4) THEN
         R0 = QUADR(I,I,0)
         R1 = QUADR(I,I,-1)/R0
         RM = QUADR(I,I,1)/R0
         EK = -D5*HL(EL,I,I,REL)/R0
      END IF
      RELS = RLSHFT(I,I)/R0
      WRITE(PRI,15)EL(I),E(I,I),AZ(I),R0,R1,RM,EK,RELS
15    FORMAT( 2X,A3,F14.7,F15.7,F10.5,F11.5,F10.5,F11.6,2F13.6)
10    CONTINUE
C
C  *****  PRINT TABLES OF 'FK' AND 'GK' INTEGRALS WHICH WERE USED IN
C  *****  DETERMINING THE ENERGY
C
      IF ( OUD .EQ. 0 ) GO TO 13
      WRITE (OUD,126)
126   FORMAT(//2X,27HVALUES OF F AND G INTEGRALS        //)
      IBEGIN = 1
      IEND = INTPTR(2)
      DO 17 I = IBEGIN,IEND
         SYMBOL = 'F'
         IF (I .GT. INTPTR(1)) SYMBOL = 'G'
17     WRITE(OUD,19) SYMBOL,KVAL(I),EL(IEL(I,1)),EL(IEL(I,2)),VALUE(I)
19       FORMAT( 2X,A1,I2,1H(,A3,1H,,A3,4H ) =, F10.7)
C
C  *****  PRINT TABLES OF 'RK' INTEGRALS
C
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
C
C  *****  PRINT TABLES OF 'L' INTEGRALS
C
      WRITE (OUD,28)
28    FORMAT(//2X,21HVALUES OF L INTEGRALS //)
      IBEGIN = IEND + 1
      IEND = INTPTR(6)
      DO 29 I = IBEGIN,IEND
29    WRITE(OUD,30) EL(IEL(I,1)),EL(IEL(I,2)),VALUE(I)
30    FORMAT(2X,2HL(,A3,1H,,A3,4H) = ,F12.7)
13    RETURN
      END
*     ------------------------------------------------------------------
*           C U P D A T E
*     ------------------------------------------------------------------
*
      SUBROUTINE CUPDATE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500,IDIM=1000,NCDIM=32768)
*
      CHARACTER CONFIG*40,EL*3,ATOM*6,TERM*6,COUPLE*3
      COMMON /LABEL/CONFIG(NCD),EL(NWD),ATOM,TERM,COUPLE(NCD,9)
*
      INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
      COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :         ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,30),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
      LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
      COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
      COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
      LOGICAL CHANGE
*
      IBEGIN =  1
      IEND = INTPTR(3)
      DO 1 I = IBEGIN,IEND
*        Omit if one of the orbitals is a continuum orbital.
*        Value(i) should remain zero in this case
*
	 IF (IEL(I,1) .EQ. NWF .OR. IEL(I,2) .EQ. NWF) GO TO 1
         IF (VARIED(IEL(I,1)) .OR. VARIED(IEL(I,2))) THEN
            IF (I .LE. INTPTR(1)) THEN
	         VALUE(I) = FK(IEL(I,1),IEL(I,2),KVAL(I),REL)
              ELSE IF (I .LE. INTPTR(2)) THEN
      		 VALUE(I) = GK(IEL(I,1),IEL(I,2),KVAL(I),REL)
            ELSE
       		 VALUE(I) = QUADR(IEL(I,1),IEL(I,2),0)**KVAL(I)
            END IF
         END IF
  1   CONTINUE
*
      IBEGIN = IEND + 1
      IEND = INTPTR(4)
      DO 30 I = IBEGIN,IEND
         CHANGE = .FALSE.
         DO 31 J = 1,4
           CHANGE = CHANGE .OR. VARIED(IEL(I,J))
 31      CONTINUE
         IF (CHANGE) THEN
            K1 = KVAL(I)/64
            K2 = KVAL(I) - 64*K1
              VALUE(I) = QUADR(IEL(I,1),IEL(I,2),0)**K1
     :                  *QUADR(IEL(I,3),IEL(I,4),0)**K2
         END IF
 30   CONTINUE
      IBEGIN = IEND + 1
      IEND = INTPTR(5)
      DO 10 I = IBEGIN,IEND
        CHANGE = .FALSE.
        DO 11 J = 1,4
           CHANGE = CHANGE .OR. VARIED(IEL(I,J))
 11     CONTINUE
        IF (CHANGE) VALUE(I)
     :        = RK(IEL(I,1),IEL(I,2),IEL(I,3),IEL(I,4),KVAL(I),REL)
 10   CONTINUE
*
      IBEGIN = IEND + 1
      IEND = INTPTR(6)
      DO 20 I = IBEGIN,IEND
	 IF (IEL(I,1) .EQ. NWF .OR. IEL(I,2) .EQ. NWF) GO TO 20
         IF (VARIED(IEL(I,1)) .OR. VARIED(IEL(I,2)))
     :                VALUE(I) = HLC(EL,IEL(I,1),IEL(I,2),REL)
 20   CONTINUE
*
*      ... Test if any of the core functions have changed
*
      CHANGE = .FALSE.
      DO 35 I = 1,NCLOSD
         CHANGE = CHANGE .OR. VARIED(I)
  35  CONTINUE
      IF (CHANGE .OR. EC.EQ.D0) CALL ECORE(EL,EC,REL)
*
      DO 40 I = 1,NWF
         VARIED(I) = .FALSE.
 40   CONTINUE
*      print *, ' Values of Integrals'
*      do 100 i = 1, intptr(6)
*	print *, i, (iel(i,j),j=1,4), value(i)
*100   continue
*      print *, ' Value of EC =', ec
      END
      DOUBLE PRECISION FUNCTION COEF(INT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500,IDIM=1000,NCDIM=32768)
*
      INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
      COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :         ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
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
  1   CONTINUE
      END
      DOUBLE PRECISION FUNCTION COV(M)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500,IDIM=1000,NCDIM=32768)
*
      INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
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
 12      CONTINUE
 10   CONTINUE
      END
C
C     ------------------------------------------------------------------
C    3-6       D I F F
C     ------------------------------------------------------------------
C
C
C       Stores LP  in the array YK.  The difference approximation of
C                i
C   Eq. (6-14) is used.
C
C
      SUBROUTINE DIFF(I)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500,IDIM=1000,NCDIM=32768)
*
      INTEGER IN,OUT,ERR,PRI,OUC,OUD,OUF,OUH
      COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,OUC,OUD,OUF,OUH
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,30),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
      COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
C  *****  FORM DD + 2Z/R -L(L+1)/RR|P(I)>
C
      MM = MAX(I) - 3
      FL = L(I)
      TWOZ = Z + Z
      C = (FL+D5)**2
      HH = 180.D0*H*H
      DO 11 K =  4,MM
11    YK(K) = (D2*(P(K+3,I)+P(K-3,I)) - 27.D0*(P(K+2,I)+P(K-2,I)) +
     1   270.D0*(P(K+1,I)+P(K-1,I)) - 490.D0*P(K,I))/HH +
     2   P(K,I)*(TWOZ*R(K) - C)
C
C  *****  BECAUSE OF THE POSSIBILITY OF EXTENSIVE CANCELLATION NEAR THE
C  *****  ORIGIN, SEARCH FOR THE POINT WHERE THE ASYMPTOTIC BEHAVIOUR
C  *****  BEGINS AND SMOOTH THE ORIGIN.
C
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
C
C  *****  ASYMPTOTIC REGION HAS BEEN FOUND
C
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
      DOUBLE PRECISION FUNCTION E(I,J)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500,IDIM=1000,NCDIM=32768)
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
 10   CONTINUE
      END
      SUBROUTINE EIJSET(I,J,E)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500,IDIM=1000,NCDIM=32768)
*
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
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
 30      IF (IJE(IP) .LT. J .AND. IP .LE. IEND) THEN
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
 40   CONTINUE
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
 50   CONTINUE
      END
      SUBROUTINE INTGRL
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
      INTEGER KVAL, IEL, CPTR, IH, JH, OPTR
      COMMON/STATE/WT(NCD),INTPTR(6),KVAL(IDIM),IEL(IDIM,4),CPTR(IDIM)
     :         ,VALUE(IDIM),COEFF(NCDIM),IH(NCDIM),JH(NCDIM),OPTR(NCDIM)
*
      CHARACTER END*1, EL1*3, EL2*3, EL3*3, EL4*3
*
    1 FORMAT(1X,A1,I2,1X,A3,1X,A3,1X,I5)
    2 FORMAT(1X,A1,I2,1X,2A3,1X,2A3,1X,I5)
    3 FORMAT(1X,A1,I2,1X,A3,1X,A3,2X,I2,1X,A3,1X,A3,I5)
    4 FORMAT(F14.8,A1,2I3,I3)
C
C ***** READ  THE LIST OF INTEGRALS
C
      LAST = 0
      IC = 1
      I = 1
      READ(IUD,'()')
      DO 10 INT = 1,6
        IF (INT.NE.4 .AND. INT.NE.5) THEN
*
*            ...F, G, L, or O1 integrals....
*
   12      READ(IUD,1) END, KVAL(I), EL1, EL2, ICPTR
           IF (END .EQ. '*') GO TO 16
           IF (ICPTR+LAST .LE. (NCDIM)) THEN
           CPTR(I) = ICPTR + LAST
           ELSE
           PRINT *,' Too much data - current dimensions =',NCDIM
           STOP
           END IF
           CALL EPTR(EL, EL1,IEL(I,1),*999)
           CALL EPTR(EL, EL2,IEL(I,2),*999)
           I = I + 1
           IF (I .LE. (IDIM) ) GO TO 12
           PRINT *, ' Too many integrals - MAX =',IDIM
           STOP
         ELSE
   14         IF (INT.EQ.5) THEN
*
*             ... R integrals ...
*
              READ(IUD,2) END, KVAL(I), EL1, EL2, EL3, EL4, ICPTR
*
            ELSE
*
*              ... O2 integrals ...
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
 16     IF (INT .EQ. 3 .OR. INT .EQ. 4) GO TO 18
*
*     ... Read the data ...
*
   20   READ(IUD,4) COEFF(IC), END, IH(IC), JH(IC), OPTR(IC)
        IF ( END .NE. '*') THEN
          IF (INT .LE. 2) THEN
            COEFF(IC) = ACURAT(COEFF(IC))
          ELSE
*
*       ... Shift origin for overlap integrals
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
*     ... Initialize for next set ..
*
   18   INTPTR(INT) = I-1
        LAST = IC-1
   10   CONTINUE
      RETURN
*
  999   PRINT *,' Electron in ',END,'-data not found in ',
     :          'configuration data'
        STOP
      END
*
*     ------------------------------------------------------------------
*    3-19      M E T H O D
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
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,30),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
      LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
      COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
      COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
*
      COMMON P2(NOD),HQ(NOD),XX(NOD),AC(20,20),BC(20),JV(20),
     :     AZZ,PP,FN,EM,FM,EU,FU,DELTAE,M,NODE,MK,KK,NJ
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
      IF (KK .GE. 3) GO TO 51
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
      IF (KK .LE. 3)  CALL EIJSET(I,I,ED)
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
C
C
C     ------------------------------------------------------------------
C    3-20      N M R V S
C     ------------------------------------------------------------------
C
C       Given two starting values, PDE(1) and PDE(2), values of PDE(j),
C   j=3,4,...,NJ+1 are obtained by outward integration of
C               Y" = YR y + F
C   using the discretization  of  Eq.  (6-27 )  with  the  difference
C   correction.  With PDE(NJ) given, the tail procedure is applied to
C   PDE(j),j=NJ+1,  NJ+2,...,MM, where MM is determined automatically
C   and DELTA is the difference between  PDE(NJ+1)  for  outward  and
C   inward integration. (See Eq 6-32, 6-33, and 6-37 for further
C   details.)
C
C
      SUBROUTINE NMRVS(NJ,DELTA,MM,PP,F)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500,IDIM=1000,NCDIM=32768)
*
      INTEGER IN,OUT,ERR,PRI,OUC,OUD,OUF,OUH
      COMMON /INOUT/ IN,OUT,ERR,PRI,IUC,IUD,IUF,OUC,OUD,OUF,OUH
*    
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,30),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
      DIMENSION PP(NOD),F(NOD),A(150),D(150)
      EQUIVALENCE (G,G3)
C
C  *****  INTEGRATE OUTWARD TO NJ+1
C
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
C
C  *****  APPLY THE TAIL PROCEDURE
C
      K = 1
      PP(M) = -(D1 - G1)*Y1 + F(M)
      A(1) = D1 - G
      D(1) = -(D2 + D10*G)
22    RATIO = A(K)/D(K)
C
C  *****  THE INTEGER 149 IN THE NEXT STATEMENT IS THE DIMENSION OF A
C  *****  MINUS 1
C
      IF (K .GE. (150)-1 .OR. M .EQ. ND) GO TO 23
      K = K +1
      M = M+1
      G = YR(M)
      A(K) = D1 - G
      D(K) = -(D2 + D10*G) - A(K)*RATIO
      PP(M) = -PP(M-1)*RATIO + F(M)
      IF (DABS(PP(M))+DABS(PP(M-1)) .GT. TOL .OR. K .LT. 9
     :     .OR. M .LT. 120) GO TO 22
20    CON =DSQRT(EH)*DEXP(-DSQRT(DABS(G/CH-.25)/RR(M))*(R(M+1)-R(M)))
      PP(M) = PP(M)/(D(K) + CON*(D1-  YR(M+1)))
      J = M+1
      DO 2 I= J,NO
2     PP(I) = D0
      DO 3 J = 2,K
      I = M-J+1
      II = K-J+1
3     PP(I) = (PP(I)-A(II+1)*PP(I+1))/D(II)
C
C  *****  SET DELTA = DIFFERENCE OF THE TWO SOLUTIONS AT NJ+1
C  *****         MM = NUMBER OF POINTS IN THE RANGE OF THE SOLUTION
C
      DELTA = DELTA - PP(I)
      MM = M
      RETURN
23    WRITE(OUT,24)
24    FORMAT(6X,52HWARNING: FUNCTIONS TRUNCATED BY NMRVS IN TAIL REGION)
      GO TO 20
      END
C
C     ------------------------------------------------------------------
C    3-21      N O D E C
C     ------------------------------------------------------------------
C
C      Counts the number of nodes of the function PDE(j) in the range
C   j = 40,...,M-10.   The node counting procedure counts the local max
C   and min values.   Only nodes between sufficiently large max and
C   min values are counted.
C
C
      FUNCTION NODEC(M)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500,IDIM=1000,NCDIM=32768)
*
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,30),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
      COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
C
C   ***** FIND MAX|PDE(J)|
C
       MM = M - 10
      DM = 0.D0
      DO 1 J = 40,MM
1     DM = DMAX1(DM, DABS(PDE(J)))
C
C   *****  COUNT THE NUMBER OF LOCAL MAX OR MIN'S
C
      NCC = 0
      SIGN = 0.D0
      DIFF1 = PDE(40) - PDE(39)
      DO 2 J = 40, MM
      DIFF2 = PDE(J+1) - PDE(J)
      IF (DIFF2*DIFF1 .GT. 0.D0 .OR. DIFF1 .EQ. 0.D0) GO TO 2
C
C   *****  A MAX OR MIN HAS BEEN FOUND.   TEST IF IT IS
C          SUFFICIENTLY LARGE
C
      IF ( DABS(PDE(J))/DM .LT. 0.05D0 ) GO TO 2
C
C   ***** CHECK IF THIS IS THE FIRST SIGNIFICANT MAXIMUM
C
      IF (SIGN .NE. 0.D0 ) GO TO 4
      M = J
      GO TO 3
C
C   ***** IF NOT THE FIRST, TEST WHETHER A SIGN CHANGE HAS
C         OCCURRED SINCE THE LAST SIGNIFICANT MAX OR MIN
C
4     IF (PDE(J)*SIGN .GT. 0.D0 ) GO TO 2
      NCC = NCC + 1
C
C   ***** RESET FOR THE NEXT NODE
C
3     SIGN = PDE(J)
2     DIFF1 = DIFF2
      NODEC = NCC
      RETURN
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
*
	LOGICAL DIAG
	COMMON AC(20,20),BC(20)
*
      IF (NWF .EQ. 1 .OR. IB .GT. NWF) RETURN
      II = MAX0(2,IB)
      DO 2 I = II,NWF
         DIAG = .TRUE.
	 IBEGIN = IEPTR(I-1)+1
	 IP = IBEGIN
	 IJ = 0
 60	 JV = IJE(IP)
	 IF (JV .LT. I .AND. IP .LE. IEPTR(I)) THEN
	    IJ = IJ+1
	    IF ( IJ .GT. (NWD)) 
     :          STOP ' TOO MANY ORTHOGONALITY CONDITIONS'
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
         IF ( .NOT. DIAG .AND. IJ.GT.1) CALL LINEQN(20,IJ,AC,BC)
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


	 IF (CTOTAL .GT. 1.D-10  ) THEN
c
c	only for bound orbitals
c
	    if( e(i,i) .gt. d0 ) then

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
	    endif

            AZ(I) = AZZ
	    VARIED(I) = .TRUE.
	 END IF
      END IF
2     CONTINUE
      END
C
C     ------------------------------------------------------------------
C    3-24      P O T L
C     ------------------------------------------------------------------
C
C       Computes and stores the potential function
C                              2(k-1)
C              YR = SUM  a    Y      (j,j;r)
C                   j,k   ijk
C
      SUBROUTINE POTL(I)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500,IDIM=1000,NCDIM=32768)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,30),YK(NOD),YR(NOD)
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
5        CONTINUE
4           CONTINUE
         END IF
2     CONTINUE
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
 12         CONTINUE
         END IF
 10   CONTINUE
      END
C
C     ------------------------------------------------------------------
C    3-25      Q U A D
C     ------------------------------------------------------------------
C
C       Evaluates the integral of F(r)G(r) with respect to r , where
C   F(r) and G(r) have the same asymptotic properties as P (r).   The
C                                                         i
C   composite Simpson's rule is used.   The integrand is zero for r >
C   r  .
C    M
C
      DOUBLE PRECISION FUNCTION QUAD(I,M,F,G)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NOD=220,NCD=500,IDIM=1000,NCDIM=32768)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,30),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
      DIMENSION F(220),G(220)
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
*   --------------------------------------------------------------------
*               R E O R D
*   --------------------------------------------------------------------
*
*       Reorder the list of first appearance so that the functions to be
*   iterated appear last in the list.
*
        SUBROUTINE REORD(OF, ELC, NWF, IERR)
	PARAMETER (NWD=30)
        CHARACTER*3 OF(NWD), ELC
*
        IERR = 1
        CALL EPTR(OF, ELC, I, *99)
        DO 10 J = I, NWF-1
           OF(J) = OF(J+1)
10      CONTINUE
        OF(NWF) = ELC
        IERR = 0
99      RETURN
        END
C
C     ------------------------------------------------------------------
C    3-32      S E A R C H
C     ------------------------------------------------------------------
C
C       This routine searches for the NJ>70 such that YR(j) > 0 for all
C   j > NJ.
C
C
      SUBROUTINE SEARCH(NJ,I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NWD=30,NOD=220)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,30),YK(NOD),YR(NOD)
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
*   ---------------------------------------------------------------------
*               S E T O R T
*   ---------------------------------------------------------------------
*
      LOGICAL FUNCTION SETORT(EL1,EL2)
      CHARACTER*3 EL1,EL2
      CHARACTER*1 S1, S2
C
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
C
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
*    3-34      S O L V E
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
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,30),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
      COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
      COMMON P2(NOD),HQ(NOD),XX(NOD),AC(20,20),BC(20),JV(20),
     :     AZZ,PP,FN,EM,FM,EU,FU,DELTAE,M,NODE,MK,KK,NJ
*
      LOGICAL FIRST
      DIMENSION ZERO(NOD),P1(NOD)
      EQUIVALENCE (ZERO(1),XX(1)),(PDE(1),P1(1))
      SAVE Zinf,fl,V,b4,CN,C,CD,XY,XP
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
      IF (KK .LT. 3) THEN
         GO TO 80
      ELSE IF (KK .GT. 3) THEN
         GO TO 18
      END IF
*      DO 11 JJ = 10,ND
*      J = NO - JJ
*      IF (YK(J) .LT. D0 ) GO TO 63
*11    CONTINUE
*      WRITE(OUT,12)
*12    FORMAT(10X,'POTENTIAL FUNCTION TOO SMALL - 2R*(Z-Y)<(L+.5)**2')
**     STOP
*      GO TO 80
*63    EM = -YK(J)/RR(J)
*      print * , 'EM computed from Yk: ',em,j,yk(j),rr(j)
*      GO TO 81
80    EM = (ZINF/(FN + .01d0))**2
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
      IF (KK .LE. 3) THEN
         PNORM = D0
         DO 50 J = 1,M
50          PNORM = PNORM + RR(J)*HQ(J)*P1(J)
         Y1 = P1(NJ-1)
         Y2 = P1(NJ)
         Y3 = P1(NJ+1)
         DELTA = Y2 - Y1 + Y2 - Y3 +YR(NJ-1)*Y1 + D10*YR(NJ)*Y2
     :      + YR(NJ+1)*Y3 + X(NJ)
         DELTAE = HQ(NJ)*DELTA/(H*H*PNORM)
      END IF

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
C
C     ------------------------------------------------------------------
C    3-36      W A V E F N
C     ------------------------------------------------------------------
C
C       This routine initializes radial functions by the procedure
C   indicated by IND(I).
C
C         Value of IND(I)     Method
C         ---------------     ------
C             -1           Functions read from unit IU2
C              0           Screened hydrogenic functions with ZZ=Z-S(I)
C              1           Functions in memory left unchanged
C                                                  0
C   The set of functions are then orthogonalized, Y (i, i;r) and the
C   diagonal energy parameters computed, when necessary.
C
C
      SUBROUTINE WAVEFN(nend)
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
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,30),YK(NOD),YR(NOD)
     :          ,X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
*
      LOGICAL FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
      COMMON /TEST/FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED(NWD)
*
      COMMON /WAVE/EC,ED,AZD,PDE(NOD),SUM(NWD),S(NWD),DPM(NWD),
     :         ACC(NWD),METH(NWD),IEPTR(NWD),IJE(98),EIJ(98),VIJ(98),IPR
      COMMON PN,Z2,FN,M,K,ZT,ETI,EKI,AZI,PT(NOD),MT
      COMMON /zzind/ZZ(NWD),IND(NwD),IELI(5),NOCCSH(NCD)
*
      CHARACTER EL1*3,AT*6,TT*6,ATM(250)*6,TRM(250)*6,TITLE*24
C
C  *****  GENERATE ARRAYS FOR R,R*R AND SQRT(R) WITH A CONSTANT MESH
C  *****  SIZE IN THE LOG(Z*R) VARIABLE
C
	if( nend .eq. -1 ) then
      DO 1 I=1,NO
      R(I)= DEXP(RHO)/Z
      RR(I) = R(I)*R(I)
      R2(I) = DSQRT(R(I))
1     RHO = RHO + H
      RHO = RHO - NO*H
	endif
C
C  ***** READ THE WAVEFUNCTIONS
C
      IF (IUF .EQ. 0) GO TO 5
2     READ(IUF,END=5) AT,TT,EL1,MM,ZT,ETI,EKI,AZI,(PT(J),J=1,MM)

cxi
cxi	wavefunctions for each continuum energy are ended by a line
cxi	with MM=0
cxi
	if( mm .eq. 0 ) goto 5
      M = MIN0(NO,MM)
      CALL EPTR(EL,EL1,I,*2)
      IF ( I .GT. 0 .AND. IND(I) .EQ. -1) THEN
         ATM(I) = AT
         TRM(I) = TT
         MAX(I) = M
         ZZ(I)  = ZT
         C = D1
         IF ( Z .NE. ZT ) C = Z/ZT
C
C  *****  SCALE RESULTS IF CARDS ARE FOR AN ATOM WITH A DIFFERENT Z
C
         CALL EIJSET(I,I,C*C*ETI)
         AZ(I)  = AZI*C**(L(I)+1)*DSQRT(C)
         DO 11 J = 1,M
            P(J,I) = C*PT(J)
11       CONTINUE
C
C  *****  SET REMAINING VALUES IN THE RANGE = 0.
C
         IF ( M .EQ. NO ) GO TO 12
         M = M +1
         DO 13  J=M,NO
13       P(J,I) = D0
12       IND(I) = -2
      ENDIF
      GO TO 2
C
C  *****  SET PARAMTERS FOR ELECTRONS AND INITIALIZE FUNCTIONS
C
5     DO 9 I = 1,NWF
      IF (IND(I)) 7,8,9
C
C  ***** WAVE FUNCTIONS NOT FOUND IN THE INPUT DATA, SET IND = 0
C
7     IF ( IND(I) .EQ. -2 ) GO TO 4
      IF ( METH(I) .EQ. 4) GO TO 9
      IND(I) = 0
      WRITE(OUT,27) EL(I)
27    FORMAT(8X,'WAVE FUNCTIONS NOT FOUND FOR ',A3)
C
C  *****  DETERMINE ESTIMATES OF THE WAVE FUNCTIONS BY THE SCREENED
C  *****  HYDROGENIC APPROXIMATION
C
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
C
C  ***** SET THE AZ(I) VALUE
C
      AZ(I) = PN*(D2*(Z - D5*S(I))/N(I))**(L(I) + 1)
      CALL EIJSET(I,I,D0)
C
C  *****  ORTHOGONALIZE TO INNER FUNCTIONS
C
4      IF (I .EQ. 1 ) GO TO 9
      IM = I - 1
      DO 6 II =1,IM
      IF (E(I,II) .EQ. D0) GO TO 6
      PN = QUADR(I,II,0)
      if (abs(pn) .gt. 1.d0) then
	write(out,*) ' Overlap greater than unity: ', i,ii,pn
	stop
      end if
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
C
C  *****  COMPUTE ONE-ELECTRON ENERGY PARAMETERS IF THEY WERE NOT
C  *****  SPECIFIED ON INPUT.
C
      DO 15 I = 1,NWF
C     IF (E(I,I) .EQ. D0) E(I,I) = HL(EL,I,I,REL) - EKIN(I,I)
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
cxi      IF ( IUF .NE. 0) REWIND(UNIT=IUF)
      RETURN
      END
      SUBROUTINE XCH(I,IOPT)
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
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,30),YK(NOD),YR(NOD)
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
  6      CONTINUE
  4         CONTINUE
         END IF
  2   CONTINUE
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
 8          CONTINUE
         END IF
 7    CONTINUE
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
C
C  ***** COUNT THE NUMBER OF OCCURRENCES OF I
C
            IK = 0
            IF (I1 .EQ. I) IK = IK + 1
            IF (I2 .EQ. I) IK = IK + 1
            IF (J1 .EQ. I) IK = IK + 1
            IF (J2 .EQ. I) IK = IK + 1
            EXIT = .FALSE.
            DO 11 II2=1,2
            DO 12 II1=1,2
            GO TO (10, 20, 30, 40) IK
10        CONTINUE
C
C  ***** I OCCURS JUST ONCE IN RK
C
            IF (I1 .NE. I) GO TO 13
            GO TO 16
20          CONTINUE
C
C  ***** I OCCURS TWICE IN THE RK INTEGRAL
C
            IF (I1 .NE. I) GO TO 13
            IF (J1 .EQ. I) GO TO 17
C
C  ***** TEST IF THE PAIR (I1,J1) = PAIR (I2,J2)
C
            ICODE1 = 100*I1 + J1
            ICODE2 = 100*I2 + J2
            ICODE3 = 100*J2 + I2
            SAME = ICODE1 .EQ. ICODE2 .OR. ICODE1 .EQ. ICODE3
            IF ( .NOT. SAME ) GO TO 15
            GO TO 17
30          CONTINUE
C
C  ***** I OCCURS THREE TIMES IN THE RK INTEGRAL
C
C
            IF (I1 .EQ. I) GO TO 13
            CALL YKF(I2, J2, KK, REL)
            DO 33 J = 1,NO
33            X(J) = X(J) + CC*P(J,I1)*YK(J)
            CALL YKF(I1, J1, KK, REL)
            CC = D2*CC
            DO 34 J = 1,NO
34             X(J) = X(J) + CC*P(J,I2)*YK(J)
            GO TO 50
C
C  ***** I OCCURS FOUR TIMES IN RK INTEGRAL
C
40          CC = D4*CC
            GO TO 16
17          CC = D2*CC
16          EXIT = .TRUE.
15          CALL YKF(I2,J2,KK,REL)
            DO 14 J=1,NO
14             X(J) = X(J) +CC*P(J,J1)*YK(J)
            IF (EXIT) GO TO 50
13            III = I1
            I1= I2
            I2= III
            III = J1
            J1 = J2
12          J2 = III
            III = I1
            I1 = J1
            J1 = III
            III = I2
            I2= J2
11          J2= III
         END IF
50    CONTINUE
*
51    IBEGIN = INTPTR(5) + 1
      IEND = INTPTR(6)
      DO 60 INT = IBEGIN,IEND
C       ... Include only if off-diagonal ...
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
 65      CONTINUE
         DO 66 K = IABS(L(I)-L(II)),L(I)+L(II),2
            CCC = CC*CB(L(I),L(II),K)
            CALL YKF(I2,II,K,REL)
            DO 67 J = 1,NO
               X(J) = X(J) - CCC*YK(J)*P(J,II)
 67         CONTINUE
 66      CONTINUE
 64         CONTINUE
         END IF
         IF (I .LE. NCLOSD) THEN
            C = -D2*COEF(INT)
            CALL YKF(I1,I2,0,REL)
            CC = D2*C
            DO 61 J = 1,NO
        X(J) = X(J) + CC*YK(J)*P(J,I)
 61         CONTINUE
            DO 63 K = IABS(L(I)-L(I1)),L(I)+L(I1),2
        CC = C*CB(L(I),L(I1),K)
        CALL YKF(I2,I,K,REL)
        DO 68 J = 1,NO
           X(J) = X(J) - CC*YK(J)*P(J,I1)
 68     CONTINUE
        CALL YKF(I1,I,K,REL)
        DO 69 J = 1,NO
           X(J) = X(J) - CC*YK(J)*P(J,I2)
 69     CONTINUE
 63         CONTINUE
         END IF
        END IF
 60   CONTINUE
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
 72         CONTINUE
         END IF
 70   CONTINUE
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
 84      CONTINUE
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
 80   CONTINUE
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
C
C  *****  CHECK IF EXCHANGE IS ZERO: IF SO, METHOD 2 SHOULD BE USED.
C
 75     IF (METH(I) .EQ. 2 .OR. METH(I) .GT. 3) RETURN
        IF ( DABS(X(1)) + DABS(X(2)) + DABS(X(3)) .EQ. D0 ) METH(I) = 2
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
     :		     CD,FL,ZL,ZF,V,NJ,MJ,MP,IX
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
	

95	MJ = MJ + 2

	if( MJ .gt. M-2 ) then
	write(ERR,*) 
     : 	' Failed in Calculating Coulomb Functions, R(M) too small '
	write(ERR,*) ' R(M) =', R(M)
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
     :	           + 45.d0*( PHC(J+1) - PHC(J-1) ) )*h60
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
*       Written my Jinhua Xi, December, 1994
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
10	continue
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
20 	continue

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
