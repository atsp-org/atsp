*     ------------------------------------------------------------------
*
*       COMP -- A PROGRAM FOR PRINTING DOMINANT CONTRIBUTIONS
*
*                   C O P Y R I G H T -- 1994
*
*       by C. Froese Fischer
*          Vanderbilt University
*          Nashville, TN 37235 USA
*
*       July, 1984
*       
*       Computer Physics Communication, Vol. 64, 399-405 (1991)
*-----------------------------------------------------------------------
*
      PROGRAM COMP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NOD=220,NWD=30,NWD2=2*NWD,NCD=100,NCD2=2*NCD)
      INTEGER Q(5),JV(NCD2),JL(NCD2),IP(NCD2)
      CHARACTER*1 ASTER,CH1,CH2
      CHARACTER*3 FIN1
      CHARACTER*3 EL(NWD),COUPLE(9),FINT(NCD2)
      CHARACTER*6 ATOM,TERM,ALABEL
      CHARACTER*24 CFILE, NAME
      CHARACTER*66  CONFIG(NCD2), LAB1
      COMMON /STATE/LENGTH(NCD2),NCFG(2),NCF(2),IST(2)
      DIMENSION ET(NCD2),COEF(NCD2,NCD2)
      PARAMETER (IREAD=5,IWRITE=6)
      DATA ASTER/'*'/
 
*
      IU = 7
      DO 100 ICASE = 1,1
*
*  *****  DETERMINE INFORMATION ABOUT FILES
*
CSUN  iarg = iargc()
CSUN  if ( iarg .eq. 0) then
        WRITE(0,*) 'Enter name of .c file'
        READ(IREAD,'(A)') NAME
CSUN  else
CSUN    call getarg(1,NAME)
CSUN  end if
      jend = index(NAME,' ')
      CFILE = NAME(1:JEND-1)//'.c'
      OPEN(UNIT=7,FILE=CFILE,STATUS='OLD')
      M = 1
      WRITE(0,*) ' What tolerance? :'
      READ(IREAD,*) TOL
*
*  *****  READ THE CONFIGURATIONS
*
      READ(IU,10) ATOM,ALABEL,ET(1)
   10 FORMAT(3X,2A6,F16.8/)
   15 READ(IU,14,END=18) (EL(K),Q(K),K=1,5),COEF(M,1)
   14 FORMAT(5(1X,A3,1X,I2,1X),F10.8)
      READ(IU,'(9(5X,A3))',END=18) (COUPLE(J),J=1,9)
      NOCC = 0
16    IF (EL(NOCC+1) .NE. '   ' ) THEN
         NOCC = NOCC+1
         IF (NOCC .LT. (5)) GO TO 16
      END IF
      IF (NOCC .EQ. 0) GO TO 18
      CALL PACK(NOCC,EL,Q,COUPLE,LAB1)
 
*
*       1. Separate the final term
*
            J = INDEX(LAB1,' ')-1
           CH1 = LAB1(J:J)
           IF (CH1 .GE. '0' .AND. CH1 .LE. '9') THEN
                J = J-4
            ELSE
                J = J-3
            ENDIF
            FIN1 = LAB1(J+2:J+3)
            LAB1(J+1:J+4) = '   '
 
*
*       2. Delete set subscriP2
*
*           CH1 = LAB1(J:J)
*           IF (CH1.GE.'0' .AND. CH1.LE.'9') LAB1(J:J) = ' '
 
*
*       3. If after removing the final term, there are no other
*   intermediate couplings prefaced by '_' and the last coupling
*   is the same as the final term, then the coupling for the final
*   term is omitted .
*
            IF (INDEX(LAB1,'_') .EQ. 0) THEN
                CH2 = LAB1(J-2:J-1)
                IF (CH2 .EQ. FIN1) LAB1(J-2:J-1) = '  '
            ENDIF
 
      CONFIG(M) = LAB1
      K = 45
   17 IF (CONFIG(M)(K:K) .EQ. ' ') THEN
         K = K-1
         GO TO 17
      END IF
      LENGTH(M) = K
      FINT(M) = FIN1
      M = M+1
      IF (M.GT.(NCD2)) THEN
         WRITE(0,'(A,I4)')
     :       ' TOO MANY CONFIGURATIONS:  MAXIMUM IS ', (NCD2)-1
      END IF
      GO TO 15
   18 NCF(ICASE)= M-1
      CLOSE(UNIT=7)
100   CONTINUE
 
      WRITE(0,'(A/5X,A/5X,A/5X,A/A)') ' Compositions from:',
     :     '1  name.c','2  name.l','3  name.j',' Enter selection'
      READ(IREAD,*) ICASE
      IF (ICASE .EQ. 1) THEN
	 NCFG(1) = NCF(1)
	 IP(1) = 1
	 JL(1) = 1
	 IST(1) = 1
	 JV(1) = 0
         CALL OUTPUT(1,JV,JL,CONFIG,FINT,COEF,IP,TOL,ET)
      ELSE
      IF (ICASE .EQ. 2) THEN
	 CFILE = NAME(1:JEND-1)//'.l'
      ELSE 
	 CFILE = NAME(1:JEND-1)//'.j'
      END IF
*
*   ****  If LS format, read sets of coefficients
*
      IS = 1
      DO 200 ICASE = 1,1
      IL = NCF(ICASE)
      OPEN(UNIT=7,FILE=CFILE,STATUS='OLD')
      READ(7,'(A6,7X,F6.1,5X,I4,9X,I4)') ATOM,Z,NEL,NCFG(ICASE)
65    READ(7,'(//8X,I4,10X,I4)',END=190) JVV,MFOUND
      DO 63 III = 1,MFOUND
        READ(7,64) JL(IS),ET(IS),(COEF(I,IS),I=1,NCFG(ICASE))
64      FORMAT(/I6,F16.8/(7F10.7))
        IF (ICASE .EQ. 2) JL(IS) = JL(IS) + NCF(1)
        JV(IS) = JVV
        IS = IS+1
63    CONTINUE
      GO TO 65
190   IST(ICASE) = IS-1
200   CONTINUE
      IS = IS-1
      CALL SORT(IS,ET,IP)
      CALL OUTPUT(IS,JV,JL,CONFIG,FINT,COEF,IP,TOL,ET)
      END IF
      END
*
*     ------------------------------------------------------------------
*       OUTPUT
*     ------------------------------------------------------------------
*
      SUBROUTINE OUTPUT(IS,JV,JL,CONFIG,FINT,COEF,IP,TOL,ET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NOD=220,NWD=30,NCD=100,NCD2=2*NCD)
      PARAMETER (IREAD=5,IWRITE=6)
      DIMENSION ET(*)
      COMMON /STATE/LENGTH(NCD2),NCFG(2),NCF(2),IST(2)
      CHARACTER*66 CONFIG(*)
      CHARACTER*3 FINT(*)
      INTEGER INDEX(NCD2),JV(NCD2),JL(NCD2),IP(*)
      DIMENSION COEF(NCD2,NCD2)
      LOGICAL FIRST,SECOND
      DO 100  II = 1,IS
         JS = IP(II)
         JSL = JL(JS)
         ICASE = 1
         IF (JSL .GT. NCF(1)) ICASE = 2
         K = LENGTH(JSL)
           WRITE(IWRITE,'(//1X,A,2X,A,2X,F4.1,F14.8)')
     :            CONFIG(JSL)(1:K),FINT(JSL),JV(JS)/2.,ET(II)
      DO 1 I = 1,NCFG(ICASE)
         INDEX(I) = I
1     CONTINUE
      DO 2 I = 1,NCFG(ICASE)
         JP = I
         DO 3 J = I+1,NCFG(ICASE)
            IF (ABS(COEF(J,JS)) .GT. ABS(COEF(JP,JS))) JP = J
3        CONTINUE
         TEMP = COEF(I,JS)
         COEF(I,JS) = COEF(JP,JS)
         COEF(JP,JS) = TEMP
         ITEMP = INDEX(I)
         INDEX(I) = INDEX(JP)
         INDEX(JP) = ITEMP
2     CONTINUE
        IB = 0
        IF (JS .GT. IST(1)) IB = NCF(1)
          FIRST = .TRUE.
          SECOND = .FALSE.
      DO 10 I = 1,NCFG(ICASE)
         J = INDEX(I) + IB
         K = LENGTH(J)
         IF (ABS(COEF(I,JS)) .GT. TOL) THEN
            IF (FIRST) THEN
                   WRITE(IWRITE,'(1X,F12.7,2X,A,2X,A)')
     :                  COEF(I,JS),CONFIG(J)(1:K),FINT(J)
                FIRST = .FALSE.
                SECOND = .TRUE.
            ELSE IF (SECOND) THEN
                   WRITE(IWRITE,'(1X,F12.7,2X,A,2X,A)')
     :                  COEF(I,JS),CONFIG(J)(1:K),FINT(J)
                SECOND = .FALSE.
            ELSE
                   WRITE(IWRITE,'(1X,F12.7,2X,A,2X,A)')
     :                  COEF(I,JS),CONFIG(J)(1:K),FINT(J)
               SECOND = .TRUE.
            END IF
 
        END IF
10      CONTINUE
100    CONTINUE
      END
*
*     ------------------------------------------------------------------
*       SORT
*     ------------------------------------------------------------------
*
        SUBROUTINE SORT(M,ET,IP)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION ET(*),IP(*)
 
        DO 1 I = 1,M
           IP(I) = I
1       CONTINUE
 
        DO 10 I = 1,M-1
           JP = I
           DO 12 J = I+1,M
              IF (ET(J) .LT. ET(JP)) JP = J
12         CONTINUE
           TEMP = ET(I)
           ET(I) = ET(JP)
           ET(JP) = TEMP
           ITEMP = IP(I)
           IP(I) = IP(JP)
           IP(JP) = ITEMP
10      CONTINUE
        END
 
