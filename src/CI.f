*     ------------------------------------------------------------------ 
* 
	PROGRAM CI
*
*                   C O P Y R I G H T -- 1994
*
*	by C. Froese Fischer 
*	   Vanderbilt University 
*	   Nashville, TN 37235 USA 
* 
*	May, 1983 
* 
*      Computer Physics Communications, Vol. 64, 473--485 (1991).
*     ------------------------------------------------------------------ 
* 
*	A CONFIGURATION INTERACTION PROGRAM EITHER NON-RELATIVISTIC 
*	      OR IN THE BREIT-PAULI APPROXIMATION 
* 
*	The PARAMETER values in this program define the following:

*		IREAD - The unit number of standard input
*		IWRITE- The unit number of printed output
*		MD    - Maximum number of eigen-pairs
*		NCD2  - Maximum number of configuration state functions
*		NOD   - Maximum number of points in the range of a
*		      - function
*		NTERMD- Maximum number of terms
*		NWD  - Maximum number of functions (or electrons)
*		NZ    - Maximum number of configuration state functions
*		      - defining the zero-order set
* 
*     ------------------------------------------------------------------ 
* 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      PARAMETER (IREAD=5,IWRITE=6) 
      CHARACTER*1 PP, NAME(5)*24
      LOGICAL PRINT, LS , REL
*
      WRITE(6,9999)
9999  FORMAT(/20X,'==========================='/
     :        20X,' CONFIGURATION INTERACTION '/
     :        20X,'==========================='/)
CSUN  i = iargc()
999   NAME(2) = 'int.lst'
CSUN  if (i .eq. 0) then
	 WRITE(0,*) ' Name of State'
         read(5,'(A)') NAME(1)
CSUN  else
CSUN	 call getarg(1,NAME)
CSUN	 if (i .eq. 2) call getarg(2,NAME(2))
CSUN  end if
      j = index(NAME(1),' ')
      if (j .eq. 1) then
	 WRITE(0,*) ' Names may not start with a blank'
	 go to 999
      else
	 NAME(1) = NAME(1)(1:j-1)//'.c'
	 NAME(3) = NAME(1)(1:j-1)//'.w'
	 NAME(4) = NAME(1)(1:j-1)//'.l'
	 NAME(5) = NAME(1)(1:j-1)//'.j'
      end if
*
      CALL INITA
      CALL INITR
      REL = .TRUE.
      WRITE(0,*) ' Is this a relativistic calculation ? (Y/N) : ' 
      READ(IREAD,'(A1)') PP 
      IF (PP.EQ.'N' .OR. PP.EQ.'n') REL = .FALSE.
      WRITE(0,*) ' Is mass-polarization to be included ? (Y/N) : '
      READ(IREAD,'(A1)') PP
      MASS = 0
      IF (PP.EQ.'Y' .OR. PP.EQ.'y') THEN
         WRITE(0,*) '   Gradient or Slater integral form ? (G/S) : '
         READ(IREAD,'(A1)') PP
         MASS = 1
         IF (PP.EQ.'S' .OR. PP.EQ.'s') MASS = 2
      END IF
*
      OPEN(UNIT=1,FILE=NAME(1),STATUS='OLD') 
      OPEN(UNIT=2,FILE=NAME(2),STATUS='OLD') 
      OPEN(UNIT=4,FILE=NAME(3),STATUS='OLD',FORM='UNFORMATTED') 
      OPEN(UNIT=7,FILE=NAME(4),STATUS='UNKNOWN') 
      IF (REL) THEN
         OPEN(UNIT=8,FILE=NAME(5),STATUS='UNKNOWN') 
	 OPEN(UNIT=3,STATUS='SCRATCH',FORM='UNFORMATTED')
      END IF
      OPEN(UNIT=9,STATUS='SCRATCH',FORM='UNFORMATTED') 
* 
      CALL EVAL(N,REL,NZERO) 
      LS = .TRUE. 
      PRINT = .FALSE. 
      CALL LSJMAT(N, NZERO, 0, PRINT, LS) 
      IF ( .NOT. REL) GO TO 99
* 
*  *****  DETERMINE DATA ABOUT THE CASE 
* 
      WRITE(0,*) ' Maximum and minimum values of 2*J ? ' 
      READ( IREAD,*) MAXJ,MINJ 
      WRITE(0,*) ' Do you want the matrix printed? (Y or N) ' 
      READ (IREAD,'(A)') PP 
      IF (PP .EQ. 'Y' .OR. PP .EQ. 'y') THEN 
          PRINT = .TRUE. 
        ELSE 
          PRINT = .FALSE. 
      END IF 
      LS = .FALSE. 
* 
*  ***** PERFORM CALCULATION FOR EACH J VALUE 
* 
      DO 1 J = MAXJ, MINJ, -2 
         CALL LSJMAT(N, NZERO,  J,  PRINT, LS) 
    1 CONTINUE 
99    IF (REL) THEN
	 CLOSE(UNIT=3) 
         WRITE(8,'(A)') '***'
      END IF
      WRITE(7,'(A)') '***'
      CLOSE(UNIT=7) 
      CLOSE(UNIT=8) 
      CLOSE(UNIT=9) 
      END 
*
*     -----------------------------------------------------------------
*       E V A L
*     -----------------------------------------------------------------
*     Read the configuration state list and radial functions.  Call 
*   MATRIX to read the non-relativistic portion of the int.lst, then
*   read the fine-structure integrals, evaluating the integrals and 
*   storing the coefficient data.
*
*
      SUBROUTINE EVAL(NC,REL,NZERO) 
* 
      PARAMETER(IREAD=5,IWRITE=6)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      PARAMETER(NOD=220,NWD=30,NCD2=200,MD=20,NZ=200,NTERMD=10,
     :  IDIM2=1000)
      LOGICAL REL, PRINT/.false./
      CHARACTER CONFIG*66,EL*3,ATOM*6,END*1,TT*6,HEADER*72,CLSDEL*72 
      CHARACTER*3 EL1, EL2, EL3, EL4, COUPLE(9), ELC(5), ANS*1 
      INTEGER S, IQ(5) 
      COMMON /INTGRL/VALUE(IDIM2),LL(NCD2),S(NCD2),LENGTH(NCD2),EC, 
     :	   ICPTR(IDIM2),IOV(2),OVALUE(20),LSP(NCD2),INDEX(NTERMD),NTERM
      COMMON /LABEL/CONFIG(NCD2),EL(NWD) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
 
* 
*  ***** READ THE CONFIGURATIONS 
* 
      I = 1 
      READ(1,'(A/A)' ) HEADER,CLSDEL 
    1 READ(1,'(5(1X,A3,1X,I2,1X))', END= 2) (ELC(K),IQ(K),K=1,5) 
      NOCC = 0 
    8 IF (ELC(NOCC+1) .NE. '   ') THEN 
         NOCC = NOCC + 1 
         IF (NOCC .LT. 5) GO TO 8 
      END IF 
      IF (NOCC .EQ. 0) GO TO 2 
      IF (I .GT. NCD2) THEN
	 WRITE(0,*) ' Too many configurations: Max = ',NCD2
	 STOP 1
      END IF
      NEL = IQ(1)+IQ(2)+IQ(3)+IQ(4)+IQ(5)
      READ(1,'(9(5X,A3))') (COUPLE(J),J=1,9) 
      CALL PACK(NOCC,ELC,IQ,COUPLE,CONFIG(I)) 
      K = 66 
    9 IF (CONFIG(I)(K:K) .EQ. ' ') THEN 
         K = K-1 
         GO TO 9 
      END IF 
      LENGTH(I) = K 
      LAST = 2*NOCC - 1 
      LL(I) = 2*LVAL(COUPLE(LAST)(2:2)) 
      S(I) = ICHAR(COUPLE(LAST)(1:1)) - ICHAR('1') 
      I = I + 1 
      GO TO 1 
    2 NC = I - 1 
      NZERO = MIN0(NC,(NZ))
*
*  *****  Determine the list of TERMS
*
      DO 14 I = 1,NC
	 LSP(I) = 0
 14   CONTINUE
      NTERM = 0
      DO 15 I = 1,NC
	 IF (LSP(I) .EQ. 0) THEN
	    NTERM = NTERM + 1
	    LSP(I) = NTERM
	    INDEX(NTERM) = I
	    DO 16 J = I+1,NC
	       IF (LSP(J) .EQ. 0 .AND.
     :		   LL(I).EQ.LL(J) .AND. S(I).EQ.S(J)) LSP(J) = NTERM
 16	    CONTINUE
	 END IF
 15   CONTINUE
      IF (NTERM .GT. NTERMD) THEN
         WRITE(0,'(/A,I3/A,I3)') ' The number of terms is ', NTERM,
     :     ' The maximum allowed in current dimensions is ', NTERMD
        STOP 1
      END IF
*
*     WRITE(0,'(/A,A,I4)' 
*    :    ' The zero-order calculation will diagonalize',
*    :	    ' a matrix of order ', NZERO
*     WRITE(0,*) ' Enter a new value IF you wish to reduce NZERO: '
*     READ(IREAD,*) NEW
*     IF (NEW .NE. 0) NZERO = NEW
* 
*  *****  READ THE RADIAL FUNCTIONS 
* 
      I = 1 
      NCLOSD = 18 
  12  READ(4,END=13) ATOM,TT,EL(I),M,Z,ETI,EKI,AZ(I),(P(J,I),J=1,M) 
    7 FORMAT(24X,8X,A3,I6,F6.0/40X,D14.7/(5D14.7)) 
      IF (I.GT.NWD) THEN
         WRITE(0,*) 'Too many functions: MAX = ',NWD
         STOP
      END IF
      IF (I .LE. NCLOSD) THEN 
         II = 4*(I-1)+2 
         EL1 = CLSDEL(II:II+2) 
         IF (EL1 .NE. '   ' .AND. EL(I) .NE. EL1) THEN 
             STOP 'Radial functions do match list of closed shells' 
	 ELSE IF (EL1 .EQ. EL(I)) THEN
	     J = 3
	     IF (EL1(1:1) .NE. ' ') J = 2
	     NEL = 4*LVAL(EL1(J:J)) + 2 + NEL
         ELSE IF (EL1 .EQ. '   ') THEN 
             NCLOSD = I-1 
         END IF 
      END IF 
      IF (EL(I)(1:1) .NE. ' ') THEN 
         N(I) = ICHAR(EL(I)(1:1)) - ICHAR('1') + 1 
         L(I) = LVAL(EL(I)(2:2)) 
      ELSE 
         N(I) = ICHAR(EL(I)(2:2)) - ICHAR('1') + 1 
         L(I) = LVAL(EL(I)(3:3)) 
      ENDIF 
      MM = M+1 
      DO 24 J = MM,NO 
   24 P(J,I) = D0 
      MAX(I)=M 
      I = I+1 
      GO TO 12 
   13 NWF = I-1 
* 
*  *****  SET UP DATA FOR AN INTERACTION CALCULATION 
* 
      WRITE (6,3) ATOM,Z 
3     FORMAT(//3X,'ATOM = ',A6,3X,'Z = ',F3.0/) 
      IF (REL) 
     :WRITE (8,'(2X,A6,A,F5.1,A,I3,A,I3)' ) ATOM,'  Z = ',Z ,'  N = ',
     :	   NEL, '   NCFG = ',NC
      WRITE (7,'(2X,A6,A,F5.1,A,I3,A,I3)' ) ATOM,'  Z = ',Z ,'  N = ',
     :	   NEL, '   NCFG = ',NC
      DO 4 J=1,NO 
      R(J)=DEXP(RHO)/Z 
      RR(J)=R(J)*R(J) 
      R2(J)=DSQRT(R(J)) 
 4    RHO=RHO+H 
*
      IF (MASS .GT. 0) THEN
      WRITE(0,*) ' Default Rydberg constant (Y or N) ? ' 
      READ(IREAD,'(A1)') ANS 
      IF ( ANS .EQ. 'Y' .OR. ANS .EQ. 'y' ) THEN 
         IF ( Z .EQ. 1.) THEN 
            ZMU = 1. 
          ELSE IF ( Z .GT. 10.) THEN 
            ZMU = 2*Z+1 + (Z-11)/2 
          ELSE IF ( MOD(INT(Z),2) .EQ. 0 .OR. Z .EQ. 7. ) THEN 
            ZMU = 2*Z 
          ELSE 
            ZMU = 2*Z+1 
         END IF 
        ELSE 
         WRITE(0,*) ' Enter the mass of the atom' 
         READ(IREAD,*)  ZMU 
      END IF 
      RMASS = 548.579903E-6/ZMU
      END IF
* 
* ***** COMPUTE ENERGY OF THE CORE 
* 
      CALL ECORE(EL,EC,REL) 
      CALL MATRIX(NC, REL, MASS, NZERO, PRINT)
      IF (.NOT. REL) RETURN
*
      I = 1
   50 FORMAT(1X,A1,I2,1X,2A3,1X,2A3,1X,I5)
   60 FORMAT(1X,A1,3X,A3,1X,A3,1X,I5)
* 
*  ***** READ Z INTEGRALS 
* 
   71 READ(2,60,END=111) END, EL1, EL2, ICPTR(I) 
      IF ( END .EQ. '*' ) GO TO 72 
         CALL EPTR(EL, EL1, I1, *999) 
         CALL EPTR(EL, EL2, I2, *999) 
         VALUE(I) = ZETA(I1,I2) 
         I = I + 1 
         IF( I .LE. IDIM2) THEN
	    GO TO 71 
	 ELSE
	   WRITE(0,*)
     :    ' Too many integrals for current dimensions: MAX=',IDIM2
	   STOP 1
	 END IF
   72 CALL LSJPTR
* 
*  ***** READ Nk INTEGRALS 
* 
   81 READ(2,50,END=111) END, K, EL1, EL2, EL3, EL4, ICPTR(I) 
      IF ( END .EQ. '*' ) GO TO 82 
         CALL EPTR(EL, EL1, I1, *999) 
         CALL EPTR(EL, EL2, I2, *999) 
         CALL EPTR(EL, EL3, I3, *999) 
         CALL EPTR(EL, EL4, I4, *999) 
         VALUE(I) = SN(I1, I2, I3, I4, K) 
         I = I + 1 
         IF (I .LE. IDIM2) THEN
	    GO TO 81 
	 ELSE
  	    WRITE(0,*)
     :    ' Too many integrals for current dimensions: MAX=',IDIM2
	    STOP 1
	 END IF
   82 CALL LSJPTR
* 
*  ***** READ Vk INTEGRALS 
* 
   91 READ(2,50,END=111) END, K, EL1, EL2, EL3, EL4, ICPTR(I) 
      IF ( END .EQ. '*' ) GO TO 92
         CALL EPTR(EL, EL1, I1, *999) 
         CALL EPTR(EL, EL2, I2, *999) 
         CALL EPTR(EL, EL3, I3, *999) 
         CALL EPTR(EL, EL4, I4, *999) 
         VALUE(I) = VK(I1, I2, I3, I4, K) 
         I = I + 1 
         IF (I .LE. IDIM2) THEN
	    GO TO 91 
	 ELSE
	    WRITE(0,*)
     :    ' Too many integrals for current dimensions: MAX=',IDIM2
	    STOP 1
	 END IF
   92 CALL LSJPTR
* 
*  *****  READ SPIN-SPIN INTEGRALS 
* 
  101 READ(2,50,END=111) END, K, EL1, EL2, EL3, EL4, ICPTR(I) 
      IF ( END .EQ. '*' ) GO TO 102 
         CALL EPTR(EL, EL1, I1, *999) 
         CALL EPTR(EL, EL2, I2, *999) 
         CALL EPTR(EL, EL3, I3, *999) 
         CALL EPTR(EL, EL4, I4, *999) 
         VALUE(I) = SN(I1, I2, I3, I4, K) 
         I = I + 1 
         IF (I .LE. IDIM2) THEN
	    GO TO 101 
	 ELSE
	    WRITE(0,*)
     :    ' Too many integrals for current dimensions: MAX=',IDIM2
	    STOP 1
	 END IF
  102 CALL LSJPTR
  111 CLOSE(UNIT=1) 
      CLOSE(UNIT=2) 
      CLOSE(UNIT=4) 
      REWIND(UNIT=3)
* 
      RETURN 
  999 STOP 
      END 
* 
*     ------------------------------------------------------------------ 
*             F I R S T
*     ------------------------------------------------------------------ 
* 
*	   Compute the first-order corrections to the eigenvalue and
*      eigenvectors.
*
      SUBROUTINE FIRST(N, NZERO, MFOUND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      PARAMETER(NOD=220,NWD=30,NCD2=200,MD=20,NZ=200)
      COMMON H(NCD2,NZ),RLB,RUB,HD(NCD2),W(MD),U(NCD2,MD)
      DIMENSION EIGVAL(MD),EIGVEC(NCD2,MD) 
      EQUIVALENCE (W(1),EIGVAL(1)),(U(1,1),EIGVEC(1,1)) 
      CHARACTER*3 CONFIG*66, COUPLE, END*1, EL
      COMMON /LABEL/CONFIG(NCD2),EL(NWD) 
*
	DO 1 I = 1,MFOUND
	   E0 = W(I)
	   E1 = 0.D0
	   U2 = 0.D0
	   DO 2 J = NZERO+1, N
*
*	   ...  Clear the components
*
	      U(J,I) = 0.D0
  2	   CONTINUE
	   DO 4 J = NZERO+1, N
*
*	   ...  Compute inner product of zero-order vectors with 
*		first-order interactions
*
	      V = 0.D0
	      DO 6 JJ = 1,NZERO
		V = V - H(J,JJ)*U(JJ,I)
  6	      CONTINUE
	      U(J,I) = V/(HD(J)-E0)
	      E1 = E1 - U(J,I)*V
	      U2 = U2 + U(J,I)*U(J,I)
  4	   CONTINUE
	   W(I) = W(I) + E1/(1.D0 +U2)
	   SCALE = 1.D0/SQRT(1.D0+U2)
	   DO 8 J = 1,N
	      U(J,I) = SCALE*U(J,I)
  8	   CONTINUE
  1	CONTINUE
	END
* 
*     ------------------------------------------------------------------ 
*             L S J M A T 
*     ------------------------------------------------------------------ 
*     Read the non-fine structure matrix, and the fine structure 
*   corrections (if any), for the current case, find eigenvalues and 
*   corresponding eigenvectors.
* 
      SUBROUTINE LSJMAT(N, NZERO, JJ, PRINT ,LS) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INTEGER L,S
      LOGICAL PRINT,LS 
      PARAMETER(NOD=220,NWD=30,NCD2=200,MD=20,NZ=200,NTERMD=10,
     :          IDIM2=1000)
      COMMON /INTGRL/VALUE(IDIM2),L(NCD2),S(NCD2),LENGTH(NCD2),EC, 
     :	   ICPTR(IDIM2),IOV(2),OVALUE(20),LSP(NCD2),INDEX(NTERMD),NTERM 
      COMMON H(NCD2,NZ),RLB,RUB,HD(NCD2),W(MD),U(NCD2,MD),V1(NZ),V2(NZ),
     :    V3(NZ),V4(NZ),V5(NZ),V6(NZ),D(NZ),E(NZ),E2(NZ),IND(NZ) 
      DIMENSION EIGVAL(MD),EIGVEC(NCD2,MD),FLSJ(NTERMD,NTERMD,2) 
      EQUIVALENCE (W(1),EIGVAL(1)),(U(1,1),EIGVEC(1,1)) 
      CHARACTER*3 CONFIG*66, END*1, EL
      COMMON /LABEL/CONFIG(NCD2),EL(NWD) 
* 
*  ***** READ THE LS-INTERACTION MATRIX 
* 
      READ (9) ((H(I,J),I=J,N),J=1,NZERO), (HD(J),J=NZERO+1,N) 
       REWIND 9 
	IF (LS) GO TO 30 
* 
*  *****  INCLUDE ONLY THOSE INTERACTIONS FOR WHICH 
*         |L - S| <= JJ <= L + S 
* 
       DO 1 J = 1,N
          IF ( JJ .LT. ABS(L(J) - S(J)) .OR. JJ .GT. L(J) + S(J) ) THEN
	     IF (J .LE. NZERO) THEN
	       DO 2 I = J+1,N
		 H(I,J) = 0.
   2	       CONTINUE
	       H(J,J) = 1.D0
	     ELSE
		HD(J) = 1.D0
	     END IF
          END IF 
    1  CONTINUE 
*
*  *****  CLEAR TABLE OF J-DEPENDENT FACTORS
*
      DO 3 II = 1,NTERM
	 I = INDEX(II)
	 DO 4 JJJ = 1,NTERM
	    J = INDEX(JJJ)
            PHASE =(-1)**((L(I)+S(J)-JJ +L(I)+S(I)-JJ +L(J)+S(J)-JJ)/2) 
            CALL GRACAH(L(J),S(J),L(I),S(I),JJ,2,W1) 
            CALL GRACAH(L(J),S(J),L(I),S(I),JJ,4,W2) 
	    FLSJ(II,JJJ,1) = PHASE*W1
	    FLSJ(II,JJJ,2) = PHASE*W2
 4       CONTINUE
 3    CONTINUE
* 
*  ***** ADD THE LSJ CONTRIBUTIONS 
* 
C      READ(3,'()' ) 
      INT = 1
      V = VALUE(INT)
      LAST = ICPTR(INT)
      DO 10 II = 1,4 
      IC = 1
   11 READ(3) C,END,I,J 
*  12 FORMAT(F14.8,A1,2I3) 
      IF ( END .EQ. '*' ) GO TO 10 
      IF (J .GT. NZERO .AND. .NOT.(J .EQ. N) ) THEN
	 WRITE(0,'(A,I3,A,I3,A)') ' Data for (',I,',',J,') IGNORED!'
	 GO TO 13
      END IF
	 IF (II .EQ. 4) THEN
	    FACTOR = FLSJ(LSP(I),LSP(J),2)
	 ELSE
	    FACTOR = FLSJ(LSP(I),LSP(J),1)
	 END IF
         H(I,J) = H(I,J) + C*FACTOR*V
   13    IC = IC + 1
	 IF (IC .GT. LAST) THEN
	    INT = INT + 1
	    V = VALUE(INT)
	    LAST = ICPTR(INT)
	 END IF
         GO TO 11 
   10 CONTINUE 
      REWIND 3 
   30      IF ( .NOT. PRINT ) GO TO 32 
* 
*  *****  PRINT THE MATRIX 
* 
      WRITE(6,'(//A)') '   LSJ interaction matrix' 
      DO 70 J=1,N 
      IF (J .LE. NZERO) THEN
         WRITE(6,71) (H(J,K),K=1,J) 
      ELSE
	 WRITE(6,71) (H(J,K),K=1,NZERO),HD(J)
      END IF
 70   CONTINUE
 71   FORMAT(/(8F16.7)) 
* 
*  ***** COMPUTE THE EIGENVALUES AND EIGENVECTORS 
* 
* 
*     NOW THAT THE INPUT MATRIX H IS READIED, WE MAY CALL THE 
*     EISPACK ROUTINE TRED1 WHICH REDUCES H TO A SYMMETRIC 
*     TRIDIAGONAL MATRIX USING ORTHOGONAL SIMILARITY 
*     TRANSFORMATIONS 
* 
32    CALL TRED1(NCD2,NZERO,H,D,E,E2) 
* 
*     PREPARE TO CALL THE EISPACK ROUTINE BISECT 
* 
      EPS1=-1D0 
* 
*     CALL THE EISPACK ROUTINE BISECT WHICH WILL USE THE 
*     TRIDIAGONAL MATRIX FOUND BY TRED1 TO ZERO IN ON ALL 
*     OF THE EIGENVALUES IN THE RANGE FROM RLB TO RUB 
* 
      CALL BISECT(NZERO,EPS1,D,E,E2,RLB,RUB,MD,MFOUND,W,IND,IERR,V4,V5) 
      IF (IERR .NE. 0 ) THEN 
* 
*     PRINT OUT THE ERROR CONDITION INDICATOR SET BY BISECT 
* 
         WRITE (6,400) IERR 
400      FORMAT(/1H1,'  IERR (FROM BISECT) =',I5) 
         STOP 
      END IF 
* 
*    CONTINUE  IF NO EIGENVALUES WERE FOUND IN THE DESIRED RANGE, 
* 
      IF (MFOUND .EQ. 0) GO TO 50 
      WRITE (6,73) MFOUND 
73    FORMAT(// I5,' EIGENVALUES FOUND ') 
      IF (LS) THEN 
         IU = 7 
      ELSE 
         IU = 8 
      END IF 
      WRITE (IU, '(//A8,I4,2X,A8,I4)' ) '  2*J = ',JJ,'NUMBER =',MFOUND 
* 
*     THE EISPACK ROUTINE TINVIT WILL NEXT BE CALLED TO FIND 
*     THE EIGENVECTORS (OF THE TRIDIAGONAL MATRIX FOUND BY TRED1) 
*     CORRESPONDING TO THE EIGENVALUES FOUND BY ROUTINE BISECT 
* 
      CALL TINVIT(NCD2,NZERO,D,E,E2,MFOUND,W,IND,U,IERR,V1,V2,V3,V4,V6) 
      IF ( IERR .NE. 0 ) THEN 
* 
*     PRINT OUT THE ERROR CONDITION INDICATOR SET BY TINVIT 
* 
         WRITE(6,460)IERR 
460      FORMAT('-','  IERR (FROM TINVIT) =',I5) 
         STOP 
      END IF 
* 
*     THE EIGENVECTORS FOUND BY TINVIT WILL NOW BE BACKTRANSFORMED 
*     BY THE EISPACK ROUTINE TRBAK1 TO FORM THE DESIRED EIGENVECTORS 
*     OF THE ORIGINAL INPUT MATRIX H, USING THE INFORMATION ABOUT THE 
*     ORTHOGONAL TRANSFORMATIONS USED IN REDUCING H TO TRIDIAGONAL 
*     FORM (THE LOWER PART OF H NOW CONTAINS THIS INFORMATION) 
* 
      CALL TRBAK1(NCD2,NZERO,H,E,MFOUND,U) 
* 
*  *****  COMPUTE FIRST-ORDER CORRECTIONS
*
      IF (NZERO .LT. N) CALL FIRST(N,NZERO,MFOUND)
*
*  *****  PRINT OUT THE EIGENVALUES AND EIGENVECTORS 
* 
      DO 40 K = 1,MFOUND 
* 
*  *****  SEARCH FOR THE LARGEST COMPONENT IN THE EIGENVECTOR FOR 
*         LABELLING PURPOSES
*
      VMAX = 0.D0 
      JMAX = 0 
      DO 90 J = 1,N 
      ABSEIG = DABS(EIGVEC(J,K)) 
      IF ( ABSEIG .LT. VMAX ) GO TO 90 
      JMAX = J 
      VMAX = ABSEIG 
90    CONTINUE 
      IF (EIGVEC(JMAX,K) .LT. 0.D0 ) THEN 
         DO 91 J = 1,N 
            EIGVEC(J,K) = - EIGVEC(J,K) 
91       CONTINUE 
      END IF 
      WRITE (IU,204)  JMAX, EIGVAL(K), CONFIG(JMAX)(1:LENGTH(JMAX)) 
204   FORMAT(/I6,F16.8,2X,A) 
      WRITE (IU,'(7F10.7)') (EIGVEC(J,K),J=1,N) 
40    WRITE(6,41) EIGVAL(K),CONFIG(JMAX)(1:LENGTH(JMAX)), 
     :            (EIGVEC(J,K),J=1,N) 
41    FORMAT(/1X,F14.8,2X,A/(1X,7F10.6)) 
      RETURN 
50    WRITE (6,51) RLB,RUB,JJ
51    FORMAT(1X,'NO EIGENVALUES FOUND  IN (',F14.7,',',F14.7, 
     :    ') FOR J = ',I3) 
      RETURN 
      END 
*
*-----------------------------------------------------------------------
*               L S J P T R
*-----------------------------------------------------------------------
*       Read the J-dependent data and write onto a scratch file
*
        SUBROUTINE LSJPTR
        CHARACTER*1 END
	DOUBLE PRECISION C
*
  1     READ(2,12) C, END, I, J
 12     FORMAT(F14.8,A1,2I3)
        WRITE(3) C, END, I, J
        IF (END .NE. '*') GO TO 1
        END
*
*----------------------------------------------------------------------- 
*      		M A T R I X 
*----------------------------------------------------------------------- 
*	This subroutine computes and stores on disk the LS interaction 
*     matrix. 
* 
      SUBROUTINE MATRIX(N, REL, MASS, NZERO, PRINT) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      PARAMETER (IREAD=5,IWRITE=6) 
      LOGICAL REL, PRINT
      CHARACTER CONFIG*66,EL*3,END*1
      CHARACTER*3 EL1, EL2, EL3, EL4
      INTEGER S
      PARAMETER(NOD=220,NWD=30,NCD2=200,MD=20,NZ=200,NTERMD=10,
     :          IDIM2=1000)
      COMMON /INTGRL/VALUE(IDIM2),LL(NCD2),S(NCD2),LENGTH(NCD2),EC, 
     :     ICPTR(IDIM2),IOV(2),OVALUE(20),LSP(NCD2),INDEX(NTERMD),NTERM
      COMMON /LABEL/CONFIG(NCD2),EL(NWD) 
      COMMON H(NCD2,NZ),RLB,RUB,HD(NCD2),DIAG(NZ) 
*
    1 FORMAT(1X,A1,I2,1X,A3,1X,A3,1X,I5)
    2 FORMAT(1X,A1,I2,1X,2A3,1X,2A3,1X,I5)
    3 FORMAT(1X,A1,I2,1X,A3,1X,A3,2X,I2,1X,A3,1X,A3,I5)
    4 FORMAT(F14.8,A1,3I3)
* 
*  ***** INITIALIZE  THE INTERACTION MATRIX TO ZERO 
* 
      DO 5 J=1,N
	IF (J .LE. NZERO) THEN
          DO 6 I = 1,N 
             H(I,J) = 0. 
    6     CONTINUE 
	  H(J,J) = EC
        ELSE
          HD(J) = EC
        END IF
    5 CONTINUE 
C
C ***** READ AND ADD THE LIST OF INTEGRALS
C
	READ(2,'()')
	DO 10 INT = 1,6
	  I = 1
	  IF (INT.NE.4 .AND. INT.NE.5) THEN
*
*            ...F, G, L, or O1 integrals....
*
   12	     READ(2,1) END, KVAL, EL1, EL2, ICPTR(I)
	     IF (END .EQ. '*') GO TO 16
	     CALL EPTR( EL, EL1,IEL1,*999)
	     CALL EPTR( EL, EL2,IEL2,*999)
	     IF (INT .EQ. 1) THEN
		VALUE(I) = FK(IEL1,IEL2,KVAL,REL)
	     ELSE IF (INT .EQ. 2) THEN
		VALUE(I) = GK(IEL1,IEL2,KVAL,REL)
 	     ELSE IF (INT .EQ. 3) THEN
		OVALUE(I) = QUADR(IEL1,IEL2,0)**KVAL
	     ELSE
		VALUE(I) = HLC(EL, IEL1, IEL2, REL)
	     END IF
	     I = I + 1
	     IF (INT .NE. 3) THEN
	        IF (I .LE. (IDIM2) ) GO TO 12
	        WRITE(0,*) ' Too many integrals - MAX = ',(IDIM2)
	     ELSE
		IF (I .LE. (20) ) GO TO 12
		WRITE(0,*) ' Too many overlap integrals - MAX = ',(20)
	     END IF
	   ELSE
   14         IF (INT.EQ.5) THEN
*
*	        ... R integrals ...
*
     	        READ(2,2) END, KVAL, EL1, EL2, EL3, EL4, ICPTR(I)
*
	      ELSE
*
*	         ... O2 integrals ...
*
    	     	READ(2, 3) END, K1, EL1, EL2, K2, EL3, EL4
	      END IF
*
	     IF ( END .EQ. '*') GO TO 16
             CALL EPTR( EL, EL1, IEL1, *999)
             CALL EPTR( EL, EL2, IEL2, *999)
             CALL EPTR( EL, EL3, IEL3, *999)
             CALL EPTR( EL, EL4, IEL4, *999)
	     IF (INT .EQ. 5) THEN
		VALUE(I) = RK( IEL1, IEL2, IEL3, IEL4, KVAL, REL)
	     ELSE
		OVALUE(I+IOV(1))
     :			= QUADR(IEL1,IEL2,0)**K1*QUADR(IEL3,IEL4,0)**K2
	     END IF
	     I = I + 1
	     IF (INT .NE. 3) THEN
	        IF (I .LE. (IDIM2) ) GO TO 14
	        WRITE(0,*) ' Too many integrals - MAX = ',IDIM2
		STOP 1
	     ELSE
		IF (I .LE. (20) ) GO TO 14
		WRITE(0,*) ' Too many overlap integrals - MAX = ',(20)
		STOP 1
	     END IF
	  END IF
 16	  IF (INT .EQ. 3 .OR. INT .EQ. 4) THEN
	    IOV(INT-2) = I-1
	    GO TO 10
	  END IF
*
*	... Read the data ...
*
	  I = 1
	  IC = 1

   20	  READ(2,4) COEFF, END, IH, JH, IOVPTR
	  IF ( END .NE. '*') THEN
	    IF (IOVPTR .LT. 0) IOVPTR = IOV(1) - IOVPTR
	    C = COEFF*VALUE(I)
	    IF (IOVPTR .NE. 0) C = C*OVALUE(IOVPTR)
	    IF (JH .LE. NZERO) THEN
	      H(IH,JH) = H(IH,JH) + C
	    ELSE
	      IF(IH .EQ. JH) HD(JH) = HD(JH) + C
	    END IF
	    IC = IC + 1
	    IF (IC .GT. ICPTR(I)) I = I+1
	    GO TO 20
	  END IF
   10 	CONTINUE
 
* 
*  *****  WRITE THE LS MATRIX ONTO SCRATCH DISK 
* 
      WRITE (9) ((H(I,J),I=J,N),J=1,NZERO), (HD(J),J=NZERO+1,N) 
      ENDFILE 9 
      REWIND 9 
      IF (.NOT. PRINT) GO TO 32
* 
*  *****  PRINT THE LS MATRIX 
* 
      WRITE (6,'(//A)') '   LS interaction matrix' 
      DO 30 I = 1,N 
	IF (I .LE. NZERO ) THEN
         WRITE(6,'(/(8F16.7))') (H(I,J),J=1,I) 
	ELSE
	 WRITE(6,'(/(8F16.7))') (H(I,J),J=1,NZERO),HD(I)
	END IF
30    CONTINUE 
32    WRITE(0,'(/A,I5/A)') 'The size of the matrix is ', N,
     :   ' Enter the approximate number of eigenvalues required ' 
      READ( IREAD, *) MEIV 
      IF (MEIV .GT. (MD)) THEN 
	WRITE(0,*)' Maximum for current dimensions is',MD,' :Re-enter' 
        GO TO 32 
      END IF 
* 
*  *****  DETERMINE THE APPROXIMATE RANGE OF THE MEIV LOWEST 
*         EIGENVALUES 
* 
      DO 33 I = 1,NZERO 
         DIAG(I) = H(I,I) 
   33 CONTINUE 
      NUMBER = 0 
      DO 34 I = 1,MIN(MEIV+1,NZERO) 
         DL = DIAG(I) 
         K = I 
         DO 35 J = I+1,NZERO 
            IF (DIAG(J) .LT. DL) THEN 
               DL = DIAG(J) 
               K = J 
            END IF 
   35    CONTINUE 
         DIAG(K) = DIAG(I) 
         DIAG(I) = DL 
	 IF (DL .NE. 1.D0) NUMBER = NUMBER + 1 
   34 CONTINUE 
      IF (NUMBER .EQ. 0) STOP ' ERROR IN MATRIX' 
      RLB = 1.5*DIAG(1) 
      IF (MEIV .LT. NUMBER) THEN 
         RUB = 0.5*(DIAG(MEIV) + DIAG(MEIV+1)) 
      ELSE 
         RUB = 2*DIAG(NUMBER)/3 
      END IF 
      RETURN 
*
  999   WRITE(0,*) ' Electron in ',END,'-data not found in ',
     :          'configuration data'
        STOP
	END
