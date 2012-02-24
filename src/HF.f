*     ------------------------------------------------------------------
*
*     A GENERAL HARTREE-FOCK PROGRAM
*
*                    C O P Y R I G H T -- 1996
*
*     by Gediminas GAIGALAS 
*        Institute of Theoretical Physics and Astronomy
*        A. Gostauto str. 12
*        Vilnius, 2600, Lithuania
*       
*     and
*             
*        C. Froese Fischer
*        Vanderbilt University
*        Nashville, TN 37235 USA
*
*     January, 1996
*     Computer Physics Communication, Vol. 98, 255-264 (1996)
*
*     This program is an extension of the earlier, HF86, published by
*     C. Froese Fischer in Computer Physics Communication, Vol. 43, 
*     355-365 (1987) to include LS terms for partially filled 
*     f-subshells
*
* -------------------------------------------------------------------
*
*     All comments in the program listing assume the radial function P
*     is the solution of an equation of the form
*
*      P" + ( 2Z/R - Y - L(L+1)/R**2 - E)P = X + T
*
*     where Y is called a potential function
*           X is called an exchange function, and
*           T includes contributions from off-diagonal energy parameter
*
*     The program uses LOG(Z*R) as independent variable and
*                      P/SQRT(R) as dependent variable.
*     As a result all equations must be transformed as described in
*     Sec. 6-2 and 6-4 of the book - ``The Hartree-Fock Method for
*     Atoms'',Wiley Interscience, 1977, by Charlotte FROESE FISCHER.
*     (All references to equations and sections refer to this book.)
*
*     Numerical procedures are the same as those for MCHF and are
*     described in Computer Physics Reports, Vol. 3, 273--326 (1986).
*     ------------------------------------------------------------------
*       M A I N    P R O G R A M
*     ------------------------------------------------------------------
*
*       The MAIN program controls the overall calculation  that
*   may  consist  of  a  series of atoms or ions in an iso-electronic
*   sequence.  Initial  estimates  for  the first are obtained either
*   from a file WFN.INP, if it exists, and scaled for the appropriate
*   Z, if necessary, or from a screened hydrogenic approximation.  All
*   others are  obtained  by  scaling  the previous results using the
*   scaling of Sec.  (7-2).
*
*
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      LOGICAL FAIL,OMIT,REL,ALL,TRACE,PRINT,STRONG
      CHARACTER CONFIG*50,EL*3,ATOM*6,TERM*6
      CHARACTER ANS*1, ASTER
      COMMON /TEST/FAIL,OMIT,REL,ALL,TRACE
      COMMON /LABEL/CONFIG,EL(NWFD),ATOM,TERM
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :  S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      INTEGER OUF
      COMMON /INOUT/ IUF,OUF
      LOGICAL OLD
      DATA ASTER/'*'/
*
*  *****  WRITE OUT HEADER
*
      WRITE(6,9)
9     FORMAT(///////22X,'=============================',
     :             /22X,' H A R T R E E - F O C K . 96',
     :             /22X,'=============================')
*
*  *****  WRITE OUT DIMENSION INFORMATION
*
      WRITE(6,99) 'NWF',NWFD,'NO',NOD
99    FORMAT(//15X,'THE DIMENSIONS FOR THE CURRENT VERSION ARE:'/13X,
     :       2(10X,2(A6,'=',I3,4X)/)/)
*
*  *****  INITIALIZE
*
      CALL INIT
*
*  ***** SET UNIT NUMBERS AND OPEN FILES
*
1     WRITE(6,'(//A/A//)') ' START OF CASE',' ============='
      INQUIRE(FILE='wfn.inp',EXIST=OLD)
      IF (OLD) THEN
         IUF =  21
         OPEN(UNIT=IUF,FILE='wfn.inp',STATUS='OLD',
     :        FORM='UNFORMATTED')
      ELSE
         IUF = 0
      END IF
      OUF = 31
      OPEN(UNIT=OUF,FILE='wfn.out',STATUS='UNKNOWN',FORM='UNFORMATTED')
      OPEN(UNIT=3,FILE='hf.log',STATUS='UNKNOWN')
*
      FAIL = .FALSE.
      DO 4 I=1,(NWFD)
      DPM(I) = D10
      DO 4 J=1,(NWFD)
      E(I,J) = D0
4     CONTINUE
*
*  *****  DETERMINE DATA ABOUT THE PROBLEM
*
      CALL DATA
*
*  *****  SET PARAMETERS TO THEIR DEFAULT VALUE
*
13    PRINT = .FALSE.
      SCFTOL = 1.D-8
      NSCF = 12
      IC = 2 + (NWF+1-IB)/4
      TRACE = .FALSE.
      IF (IB .LE. NWF) THEN
         WRITE(0,'(/A)')
     :            ' Default values for remaining parameters? (Y/N/H) '
         READ(5,'(A)') ANS
         IF (ANS .EQ. 'H' .OR. ANS .EQ. 'h') THEN
            CALL HELP(4)
            GO TO 13
         END IF
         IF (ANS .NE. 'Y' .AND. ANS .NE. 'y') THEN
*
*  *****  ADDITIONAL PARAMETERS
*
   50       WRITE(0,'(/A)') ' Default values (NO,STRONG) ? (Y/N/H) '
            READ(5,'(A)') ANS
            IF (ANS .EQ. 'H' .OR. ANS .EQ. 'h') THEN
               CALL HELP(3)
               GO TO 50
            END IF
            IF (ANS .NE. 'Y' .AND. ANS .NE. 'y') THEN
               WRITE(0,*) ' Enter values in FORMAT(I3,1X,L1) '
               READ(5,'(I3,1X,L1)') NO, STRONG
               IF (NO .GT. (NOD)) THEN
                  WRITE(0,'(A,A,I4)') ' TOO MANY POINTS: the allowed',
     :              ' MAXIMUM is ', NOD
                  GO TO 50
               END IF
               ND = NO - 2
               OMIT = .NOT. STRONG
            END IF
 16         WRITE(0,'(A)')' Default values for PRINT, SCFTOL ? (Y/N/H)'
            READ(5,'(A)') ANS
            IF ( ANS .EQ. 'H' .OR. ANS .EQ. 'h' ) THEN
               CALL HELP(5)
               GO TO 16
            ENDIF
            IF ( ANS .NE. 'Y' .AND. ANS .NE. 'y'  ) THEN
               WRITE(0,'(A)') ' Input FORMAT(L1, 1X, E6.1) '
               READ(5,'(L1,1X,E6.1)') PRINT, SCFTOL
            END IF
17          WRITE(0,'(A)') ' Default values for NSCF, IC ? (Y/N/H) '
            READ(5,'(A)') ANS
            IF (ANS .EQ. 'H' .OR. ANS .EQ. 'h' ) THEN
               CALL HELP(6)
               GO TO 17
            END IF
            IF (ANS .NE. 'Y' .AND. ANS .NE. 'y' ) THEN
               WRITE(0,'(A)') ' Input FORMAT(I2, 1X, I1) '
               READ(5,'(I2,1X,I1)') NSCF, IC
            END IF
18          WRITE(0,'(A)') ' Default values for TRACE ? (Y/N/H) '
            READ(5,'(A)') ANS
            IF (ANS .EQ. 'H' .OR. ANS .EQ. 'h' ) THEN
               CALL HELP(7)
               GO TO 18
            END IF
            IF (ANS .EQ. 'N' .OR. ANS .EQ. 'n') TRACE = .TRUE.
         END IF
      END IF
*
*
*  *****  PERFORM THE MCHF ITERATION
*
      CALL SCF(ETOTAL,SCFTOL,EREL)
*
*  *****  OUTPUT RESULTS IF PRINT = .TRUE.
*
      CALL OUTPUT(PRINT)
15    IF (FAIL) GO TO 6
      CALL SUMMRY(ETOTAL,EREL)
19    WRITE(0,'(/A)') ' Additional parameters ? (Y/N/H) '
      READ(5,'(A)') ANS
      IF (ANS .EQ. 'H' .OR. ANS .EQ. 'h') THEN
         CALL HELP(8)
         GO TO 19
      END IF
      IF (ANS .EQ. 'Y' .OR. ANS .EQ. 'y') CALL MENU
*
*  *****  CHECK FOR ISOELECTRONIC SEQUENCE OR END OF CASE.
*
20    WRITE(0,'(/A)')' Do you wish to continue along the sequence ? '
      READ(5,'(A1)') ANS
      IF (ANS .EQ. 'H' .OR. ANS .EQ. 'h') THEN
         CALL HELP(9)
         GO TO 20
      END IF
      IF (ANS .EQ. 'Y' .OR. ANS .EQ.'y') THEN
          WRITE(0,*) '   Enter: ATOM, ZZ, (ACC(I),I=1,NWF) in ',
     :             ' format(A6,F6.0,(20F3.1))'
          READ(5,'(A6,F6.0,(20F3.1))') ATOM, ZZ, (ACC(I),I=1,NWF)
*
*  *****  SCALE RESULTS FOR ANOTHER MEMBER OF THE ISOELECTRONIC SEQUENCE
*
          CALL SCALE(ZZ)
          WRITE (3,14) ATOM,TERM
14        FORMAT(1H1,9X,2A6)
          CALL ORTHOG
          GO TO 13
      END IF
*
*  *****  DETERMINE END OF CASE
*
6     WRITE(6,'(//A/A//)') ' END OF CASE',' ==========='
      STOP
      END
*
*     ----------------------------------------------------------------
*               A
*     ----------------------------------------------------------------
*
*       Determine the coefficient in the potential for electron i of
*       Y^k(j,j)

      DOUBLE PRECISION FUNCTION A(I,J,K)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      COMMON /COEFF/COEF(100),IJPTR(5,5)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :  S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
*
      IF (I.GT.NCLOSD .AND. J.GT.NCLOSD) THEN
         ISTART = IJPTR(I-NCLOSD,J-NCLOSD) + 1
         A = COEF(ISTART + K/2)
      ELSE IF (I.EQ.J) THEN
         C = SUM(I) - D1
         IF (K.EQ.0) THEN
            A = C
         ELSE
            A = -C*CA(L(I),K)
         END IF
      ELSE IF (K.EQ.0) THEN
         A = SUM(J)
      ELSE
         A = D0
      END IF
      END
*
*     ----------------------------------------------------------------
*               A D D
*     ----------------------------------------------------------------
*
*     Add a Slater integral to the data structure associated with the
*     energy expression
*
      SUBROUTINE ADD(C,K,I,J,FIRST)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      COMMON /COEFF/COEF(100),IJPTR(5,5)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :   S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
      LOGICAL FIRST
*
      IP = IJPTR(I-NCLOSD,J-NCLOSD)

      IF (FIRST) THEN
         COEF(IP+K/2+1) = C/SUM(I) + COEF(IP+K/2+1)
        ELSE
           IP = IP + MIN(L(I),L(J)) +1 + (K-ABS(L(I)-L(J)))/2 + 1
         COEF(IP) = COEF(IP) + C/SUM(I)
        END IF
      END
*
*     ----------------------------------------------------------------
*               A R R A Y
*     ----------------------------------------------------------------
*
*     Set up the data structure associated with the average energy
*
      SUBROUTINE ARRAY
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      COMMON /COEFF/COEF(100),IJPTR(5,5)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :   S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
*
      IP = 0
      DO 10 I = NCLOSD+1,NWF
         ISUMI = int(SUM(I))
         DSUMI = SUM(I) - ISUMI
         DO 12 J = NCLOSD+1,NWF
            ISUMJ = int(SUM(J))
            DSUMJ = SUM(J) - ISUMJ
            IF ( I .NE. J) THEN
         C = SUM(J)
         IF (DSUMI .NE. D0 .AND. DSUMJ .NE. D0)
     :              C = (DSUMI*(ISUMI+1)*ISUMJ +
     :                   DSUMJ*(ISUMJ+1)*ISUMI)/SUM(I)
            ELSE
         C = SUM(I) - D1
         IF (DSUMI .NE. D0)
     :              C = (ISUMI*(SUM(I)+DSUMI-1))/SUM(I)
            END IF
*
            IJPTR(I-NCLOSD,J-NCLOSD) = IP
*
*  *****        Direct contribution
*
            DO 14 K = 0,2*MIN0(L(I),L(J)),2
         IP = IP + 1
         IF (IP .GT. (100))
     :               STOP ' COEF array too small: MAX = (100)'
         COEF(IP) = D0
               IF (K .EQ. 0) THEN
                COEF(IP) = C
         ELSE IF (I .EQ. J) THEN
                COEF(IP) = -C*CA(L(I),K)
         END IF
 14         CONTINUE
*
*  *****        Exchange contribution
*
            IF (I .NE. J) THEN
            DO 16 K = ABS(L(I)-L(J)),L(I)+L(J),2
         IP = IP + 1
         IF (IP .GT. (100))
     :               STOP ' COEF array too small: MAX = (100)'
         COEF(IP) = -C*CB(L(I),L(J),K)
 16         CONTINUE
            END IF
 12      CONTINUE
 10   CONTINUE
      END
*
*     ----------------------------------------------------------------
*               B
*     ----------------------------------------------------------------
*
*     Determine the coefficient of the Y^k(i,j)P(j) term in the exchange
*     expression of electron i
*
      DOUBLE PRECISION FUNCTION B(I,J,K)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      COMMON /COEFF/COEF(100),IJPTR(5,5)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :  S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
*
      IF (I.EQ.J) THEN
         B = D0
      ELSE IF (I.GT.NCLOSD .AND. J.GT.NCLOSD) THEN
*
*   ..... LL is the number of direct terms
*         ISTART the beginning of the exchange terms
*
         LL = MIN(L(I),L(J)) + 1
         ISTART = IJPTR(I-NCLOSD,J-NCLOSD) + 1 + LL
         KK = (K - ABS(L(I)-L(J)))/2
         B = COEF(ISTART + KK)
      ELSE
         B = -SUM(J)*CB(L(I),L(J),K)
      END IF
      END
*
*     ------------------------------------------------------------------
*               B W I N T
*     ------------------------------------------------------------------
*
      SUBROUTINE BWINT(LC,LO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/BLUME/COEFN2(4),COEFNK(4),COEFVK(4)
*
* ... LC IS THE L-VALUE OF THE FILLED SUBSHELL, LO IS THE L-VALUE
*     OF THE PARTIALLY-FILLED SUBSHELL.
*
      IF(LC.LE.3.AND.LO.LE.4) GO TO 1
      WRITE(0,100) LC,LO
  100 FORMAT (37H INCORRECT CALLING OF BWINT WITH LC =,I2,6H, LO =,I2)
    1 LC1 = LC + 1
      GO TO (10,20,30,40), LC1
   10 GO TO (11,12,13,14), LO
*
* ... S-P
*
   11 COEFNK(1) = 1.D0
      COEFN2(1) = -2.D0
      COEFVK(1) = 1.D0
      RETURN
*
* ... S-D
*
   12 COEFNK(1) = 6.D0/5.D0
      COEFN2(1) = -9.D0/5.D0
      COEFVK(1) = 3.D0/5.D0
      RETURN
*
* ... S-F
*
   13 COEFNK(1) = 9.D0/7.D0
      COEFN2(1) = -12.D0/7.D0
      COEFVK(1) = 3.D0/7.D0
      RETURN
*
* ... S-G
*
   14 COEFNK(1) = 4.D0/3.D0
      COEFN2(1) = -5.D0/3.D0
      COEFVK(1) = 1.D0/3.D0
      RETURN
   20 GO TO (21,22,23,24), LO
*
* ... P-P
*
   21 COEFNK(1) = 0.D0
      COEFN2(1) = 3.D0
      COEFVK(1) = 9.D0/5.D0
      RETURN
*
* ... P-D
*
   22 COEFNK(1) = 3.D0/7.D0
      COEFNK(2) = 36.D0/35.D0
      COEFN2(1) = -12.D0/5.D0
      COEFN2(2) = 0.D0
      COEFVK(1) = 3.D0/5.D0
      COEFVK(2) = 36.D0/35.D0
      RETURN
*
* ... P-F
*
   23 COEFNK(1) = 1.D0/7.D0
      COEFNK(2) = 10.D0/7.D0
      COEFN2(1) = -18.D0/7.D0
      COEFN2(2) = 0.D0
      COEFVK(1) = 18.D0/35.D0
      COEFVK(2) = 5.D0/7.D0
      RETURN
*
* ... P-G
*
   24 COEFNK(1) = 5.D0/77.D0
      COEFNK(2) = 18.D0/11.D0
      COEFN2(1) = -18.D0/7.D0
      COEFN2(2) = 0.D0
      COEFVK(1) = 3.D0/7.D0
      COEFVK(2) = 6.D0/11.D0
      RETURN
   30 GO TO (31,32,33,34), LO
*
* ... D-P
*
   31 COEFNK(1) = 59.D0/7.D0
      COEFNK(2) = -18.D0/7.D0
      COEFN2(1) = -4.D0
      COEFN2(2) = 0.D0
      COEFVK(1) = -1.D0
      COEFVK(2) = 18.D0/7.D0
      RETURN
*
* ... D-D
*
   32 COEFNK(1) = 6.D0/7.D0
      COEFNK(2) = 0.D0
      COEFN2(1) = 3.D0
      COEFN2(2) = 0.D0
      COEFVK(1) = 3.D0/7.D0
      COEFVK(2) = 10.D0/7.D0
      RETURN
*
* ... D-F
*
   33 COEFNK(1) = 9.D0/7.D0
      COEFNK(2) = -13.D0/77.D0
      COEFNK(3) = 75.D0/77.D0
      COEFN2(1) = -18.D0/7.D0
      COEFN2(2) = 0.D0
      COEFN2(3) = 0.D0
      COEFVK(1) = 3.D0/7.D0
      COEFVK(2) = 3.D0/7.D0
      COEFVK(3) = 75.D0/77.D0
      RETURN
*
* ... D-G
*
   34 COEFNK(1) = 741.D0/693.D0
      COEFNK(2) = -215.D0/429.D0
      COEFNK(3) = 210.D0/143.D0
      COEFN2(1) = -3.D0
      COEFN2(2) = 0.D0
      COEFN2(3) = 0.D0
      COEFVK(1) = 3.D0/7.D0
      COEFVK(2) = 255.D0/693.D0
      COEFVK(3) = 105.D0/143.D0
      RETURN
   40 GO TO (41,42,43,44), LO
*
* ... F-P
*
   41 COEFNK(1) = 52.D0/3.D0
      COEFNK(2) = -20.D0/3.D0
      COEFN2(1) = -9.D0
      COEFN2(2) = 0.D0
      COEFVK(1) = -9.D0/5.D0
      COEFVK(2) = 10.D0/3.D0
      RETURN
*
* ... F-D
*
   42 COEFNK(1) = 5.D0
      COEFNK(2) = 142.D0/55.D0
      COEFNK(3) = -20.D0/11.D0
      COEFN2(1) = -18.D0/5.D0
      COEFN2(2) = 0.D0
      COEFN2(3) = 0.D0
      COEFVK(1) = -3.D0/5.D0
      COEFVK(2) = 2.D0/5.D0
      COEFVK(3) = 20.D0/11.D0
      RETURN
*
* ... F-F
*
   43 COEFNK(1) = 1.D0
      COEFNK(2) = 5.D0/11.D0
      COEFNK(3) = 0.D0
      COEFN2(1) = 3.D0
      COEFN2(2) = 0.D0
      COEFN2(3) = 0.D0
      COEFVK(1) = 1.D0/5.D0
      COEFVK(2) = 5.D0/11.D0
      COEFVK(3) = 175.D0/143.D0
      RETURN
*
* ... F-G
*
   44 COEFNK(1) = 53.D0/33.D0
      COEFNK(2) = 57.D0/143.D0
      COEFNK(3) = -115.D0/429.D0
      COEFNK(4) = 392.D0/429.D0
      COEFN2(1) = -8.D0/3.D0
      COEFN2(2) = 0.D0
      COEFN2(3) = 0.D0
      COEFN2(4) = 0.D0
      COEFVK(1) = 1.D0/3.D0
      COEFVK(2) = 3.D0/11.D0
      COEFVK(3) = 57.D0/143.D0
      COEFVK(4) = 392.D0/429.D0
      RETURN
      END
*
*     ----------------------------------------------------------------
*               B W Z E T A
*     ----------------------------------------------------------------
*
*
*  ***** COMPUTES THE NUCLEAR SPIN-ORBIT PARAMETER AND THE
*        CORRECTIONS FOR THE OTHER ELECTRONS
*      USING THE FORMULA DERIVED BY Blume and Watson.
*
      DOUBLE PRECISION FUNCTION BWZETA(I1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      COMMON/BLUME/COEFN2(4),COEFNK(4),COEFVK(4)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
*
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :   S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
      DIMENSION SS(3)
*
*
      ZETA = FINE*Z*QUADR(I1,I1,-3)
      LB = L(I1)
      DO 10 I = 1,NWF
       IF (I .EQ. I1) GO TO 10
         LA = L(I)
         ZETA = ZETA -SUM(I)*SN(I1, I, I1, I, 0)
       IF (SUM(I) .NE. 4*L(I)+2) GO TO 10
         CALL BWINT(LA,LB)
         KE1 = 2
         IF (LA .NE. LB) KE1 = IABS(LA-LB)
         IP = 0
         DO 20 K = KE1,LA+LB,2
            IP = IP+1
            ZETA = ZETA+COEFN2(IP)*SN(I1, I, I, I1, K-2)
     :                 +COEFNK(IP)*SN(I, I1, I1, I, K)
     :                 +COEFVK(IP)*(VK(I1,I,I,I1,K-1)-VK(I,I1,I1,I,K-1))
   20    CONTINUE
   10 CONTINUE
      ZETA = D2*ZETA
      C= SUM(I1)
      IF (C .NE. D1) THEN
         SS(1) = SN(I1,I1,I1,I1,0)
         C = C + C - D3
         ZETA = ZETA - C*SS(1)
         IF (LB .EQ. 2) THEN
            SS(2) = SN(I1,I1,I1,I1,2)
            ZETA = ZETA + SS(2)*6.D0/7.D0
         ELSE IF (LB .EQ. 3) THEN
            SS(2) = SN(I1,I1,I1,I1,2)
            SS(3) = SN(I1,I1,I1,I1,4)
            ZETA = ZETA + SS(2) + SS(3)/2.2D0
         END IF
      END IF
      BWZETA = ZETA
      END
*
*     ------------------------------------------------------------------
*               C A
*     ------------------------------------------------------------------
*
*
      DOUBLE PRECISION FUNCTION CA(L,K)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /EAV/CCA(10),CCB(35)
*
      IF (L .LE. 4) THEN
         CA = CCA((L*(L-1) + K)/2)
       ELSE
*        correctedx according to Prof. P. Bogdanovich 1996.03.18
         CA = RME(L,L,K)**2/((2*L+1)*(4*L+1))
      END IF
      END
*
*     -----------------------------------------------------------------
*                 C B
*     -----------------------------------------------------------------
*
*
      DOUBLE PRECISION FUNCTION CB(L,LP,K)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /EAV/CCA(10),CCB(35)
      INTEGER ICBPTR(0:4)
      DATA    ICBPTR/1,6,14,23,31/
*
             IF (L .LE. LP) THEN
                 L1 = L
                 L2 = LP
              ELSE
                 L1 = LP
                 L2 = L
             END IF
             IF ( L2 .LE. 4) THEN
                CB = CCB(ICBPTR(L1)+(K+L1-L2)/2+(L1+1)*(L2-L1))
               ELSE
                CB = RME(L,LP,K)**2/(2*(2*L+1)*(2*LP+1))
             END IF
      END
*
*     ------------------------------------------------------------------
*               D A T A
*     ------------------------------------------------------------------
*
*       Data concerning the number of configurations (NCFG), the number
*   and type of electrons in each  configuration,  as  well  as  data
*   associated with the energy expression are read and stored.
*
      SUBROUTINE DATA
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      CHARACTER CONFIG*50,ATOM*6,TERM*6,ANS*1,STRING*50
      CHARACTER*3 EL,EL1,EL2,ELCSD(18)
      LOGICAL FAIL,OMIT,REL,ALL,TRACE,FIRST,STRONG,DONE,ORDERD
      COMMON /TEST/FAIL,OMIT,REL,ALL,TRACE
      COMMON /LABEL/CONFIG,EL(NWFD),ATOM,TERM
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :   S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      COMMON ZZ(NWFD),IND(NWFD)
      INTEGER OUF
      COMMON /INOUT/ IUF,OUF
      CHARACTER*1 ASTER,W
      DATA ASTER/'*'/
*
    1 FORMAT(18(1X,A3))
    7 FORMAT(A3,F6.0,I3,I3,F3.1)
*
*  *****  READ 'ATOM' CARD
*
5     WRITE(0,'(A/A)') ' Enter ATOM,TERM,Z',
     :   ' Examples: O,3P,8. or Oxygen,AV,8.'
      READ(5,'(A50)') STRING
      I = INDEX(STRING,',')
      IF ( I .EQ. 0) THEN
          WRITE(0,*)' ATOM, TERM, and Z must be separated by commas '
          GO TO 5
      END IF
      ATOM = STRING(1:I-1)
      J = INDEX(STRING(I+1:),',')
      IF ( J .EQ. 0) THEN
          WRITE(0,*)' ATOM, TERM, and Z must be separated by commas '
          GO TO 5
      END IF
      TERM = STRING(I+1:I+J-1)
      READ(STRING(I+J+1:), '(F3.0)') Z
*
*  *****  INPUT COMMON CLOSED SHELLS
*
   2  WRITE(0,*)
      WRITE(0,'(A,A)')' List the CLOSED shells in the fields indicated',
     :        ' (blank line if none)'
      WRITE(0,'(A)') ' ... ... ... ... ... ... ... ... etc.'
      READ(5,1) (ELCSD(I),I=1,18)
*
*  *****  INPUT THE CONFIGURATION
*
      WRITE(0,'(/A,A/A)')' Enter electrons outside CLOSED shells ',
     : '(blank line if none)',' Example: 2s(1)2p(3)'
      READ(5,'(A)')  STRING
      CALL REFORM(STRING, CONFIG)
*
*      Determine the number of closed shells
*
      I = 0
      SS = D0
   12 IF (ELCSD(I+1) .NE. '   ') THEN
         I = I+1
         EL(I) = ELCSD(I)
         J = 3
         IF (EL(I)(1:1) .NE. ' ') J = 2
         L(I) = LVAL(EL(I)(J:J))
         N(I) = ICHAR(EL(I)(J-1:J-1)) - ICHAR('1') + 1
         IFULL = 2*(2*L(I)+1)
         SUM(I) = IFULL
         S(I) = SS + IFULL/2
         SS = SS + IFULL
         METH(I) = 1
         ACC(I) = D0
         IND(I) = 0
         IF (IUF .NE. 0)  IND(I) = -1
         IF( I .LT. 18) GO TO 12
         STOP ' TOO MANY CLOSED SHELLS: MAX = 18'
      END IF
      NCLOSD = I
*
*  *****  DETERMINE THE OTHER ELECTRONS
*
      MAXORB = NCLOSD
      STRING = CONFIG
      J = 2
      I = 0
 16   IF (STRING(J:J+2) .NE. '   ' ) THEN
*
*  --------- An electron has been found; is it a new one?
*
         I = I+1
         IF (I .GT. (5)) STOP ' TOO MANY SHELLS: MAX= (5)'
         EL1 = STRING(J:J+2)
         K = NCLOSD + 1
 17      IF (K .LE. MAXORB) THEN
            IF ( EL(K) .NE. EL1 ) THEN
               K = K+1
               IF (K .GT. (NWFD)) THEN
                  WRITE(0,'(A,I4)')' TOO MANY ELECTRONS: MAX =',NWFD
                  GO TO 2
               ELSE
                  GO TO 17
               END IF
            END IF
         ELSE
*
*  ------------  A new electron has been found; add it to the list
*
            MAXORB = K
            EL(MAXORB) = EL1
            READ(STRING(J+4:J+7),'(F4.0)') SUM(K)
         END IF
         J = J+10
         IF (J .LT. (50)) GO TO 16
      END IF
*
*  -----  The list of electrons has been determined
*
      WRITE(0,19) MAXORB,(EL(J),J=1,MAXORB)
   19 FORMAT(/' There are ',I3,' orbitals as follows:'/(1X,18(1X,A3)))
      NWF = MAXORB
      IF (NIT .LT. 0) NIT=NWF
21    WRITE(0,'(/A,A)') ' Orbitals to be varied: ',
     :  'ALL/NONE/=i (last i)/comma delimited list/H'
      READ(5, '(A)') STRING
      IF (STRING(1:1) .EQ. 'h' .OR. STRING(1:1) .EQ. 'H') THEN
       CALL HELP(1)
       GO TO 21
      ELSE IF (STRING(1:3) .EQ. 'ALL' .OR. STRING(1:3) .EQ. 'all') THEN
         NIT = NWF
      ELSE IF (STRING(1:4).EQ.'NONE' .OR. STRING(1:4).EQ.'none') THEN
         NIT = 0
      ELSE IF (INDEX(STRING,'=') .NE. 0) THEN
         J = INDEX(STRING,'=')
         JJ = INDEX(STRING,' ')
         READ(STRING(J+1:),'(I2)') NIT
         IF (JJ .EQ. J+2) NIT= MOD(NIT,10)
      ELSE
         NIT = 0
         J = 1
22       NEXT = INDEX(STRING(J:),',')
*
*        ***  Search for last electron label which need not be followed
*             by a comma
*
         IF (NEXT .EQ. 0 .AND. STRING(J:J+2) .NE. '   ')
     :       NEXT = INDEX(STRING(J+1:),' ') + 1
         IF (NEXT .GE. 1) THEN
            IF (NEXT .EQ. 4) THEN
               EL1 = STRING(J:J+2)
            ELSE IF (NEXT .EQ. 3) THEN
               EL1 = ' '//STRING(J:J+1)
            ELSE
               WRITE(0,*) 'Electron labels must be separated by commas;'
               WRITE(0,*) ' each label must contain 2 or 3 characters'
               GO TO 21
            END IF
            CALL REORD(EL,EL1,NWF,IERR)
            IF (IERR .EQ. 0) THEN
               NIT = NIT + 1
               J = J + NEXT
               IF (J .LT. 72) GO TO 22
            ELSE
               WRITE(0,*) ' Case must match as well as position of',
     :                  ' imbedded blanks'
               WRITE(0,*) ' For 3rd character of label to be blank',
     :                 ' follow blank with comma'

               GO TO 21
            END IF
         END IF
      END IF
*
      IB = NWF - NIT + 1
      IF (NIT .NE. 0) THEN
23       WRITE(0,'(/A)') ' Default electron parameters ? (Y/N/H) '
         READ(5,'(A)') ANS
         IF ( ANS .EQ. 'H' .OR. ANS .EQ. 'h' ) THEN
            CALL HELP(2)
            GO TO 23
         END IF
      ELSE
         ANS = 'Y'
      END IF
      IF ( ANS .NE. 'Y' .AND. ANS .NE. 'y') WRITE(0,'(A,A)')
     :   ' S, IND, METH, ACC for non-closed Shell electrons: '
      DO 20 I = NCLOSD+1,NWF
         IF ( ANS .EQ. 'Y' .OR. ANS .EQ. 'y') THEN
            S(I) = SS + (SUM(I)-D1)/D2
            SS = SS + SUM(I)
            METH(I) = 1
            ACC(I) = D0
            IND(I) = 0
            IF (IUF .NE. 0)  IND(I) = -1
         ELSE
            WRITE(0,'(A,A)') EL(I),':  '
            READ(5,*) S(I),IND(I),METH(I),ACC(I)
         END IF
         J = 2
         IF (EL(I)(1:1) .EQ. ' ') J = 3
         L(I) = LVAL(EL(I)(J:J))
         N(I) = ICHAR(EL(I)(J-1:J-1)) - ICHAR('1') + 1
         IF (IND(I) .NE. 1) THEN
            EK(I) = D0
            AZ(I) =  D0
         END IF
 20   CONTINUE
*
*  *****  DEFINE ALL ORBITALS IN THE CONFIGURATION TO BE ORTHOGONAL
*
       DO 34 I = 1,NWF
         E(I,I) = D0
         DO 33 J = 1,I-1
            E(I,J) = D0
            IF (L(I) .EQ. L(J)) E(I,J) = 1.D-5
            E(J,I) = E(I,J)
   33       CONTINUE
   34 CONTINUE
      IB = NWF - NIT + 1
      NO = NOD
      ND = NO - 2
      STRONG = .FALSE.
      WRITE (3,62) ATOM,TERM,Z,(EL(I),INT(SUM(I)),I=1,NCLOSD)
62    FORMAT(1H1///9X,33HHARTREE-FOCK WAVE FUNCTIONS FOR  ,2A6,4H Z =,
     : F5.1//14X,'Core =',5(1X,A3,'(',I4,')')/(20X,5(1X,A3,'(',I4,')')))
      WRITE (3,'(5X,A15,A50)')  'Configuration =',CONFIG
68    CONTINUE
      WRITE (3,71)
71    FORMAT(//9X,10HINPUT DATA/9X,10H----- ----//13X,13HWAVE FUNCTION,
     :   11H  PROCEDURE/17X,22HNL  SIGMA METH ACC OPT///)
      DO 79 I = 1,NWF
         WRITE (3,78) I,EL(I),N(I),L(I),S(I),METH(I),ACC(I),IND(I)
78       FORMAT(I8, 2X,A3,2I3,F7.1,I4,F4.1,I4)
79    CONTINUE
      OMIT = .NOT. STRONG
*
      CALL ARRAY
      CALL ENEXPR(TERM, DONE)
      IF (.NOT. DONE) THEN
*
*  ---  Case needs additional data
*
         WRITE(0,85)
85       FORMAT(/' The program could not derive the energy expression'/
     :           ' Select one of the following options and enter:'/
     :           '    1  Re-enter the term and configuration'/
     :           '    2  Enter the deviations from Eav as input'/
     :           '    3  STOP'/)
         READ(5,*) ISELEC
         GO TO (5,86,99) ISELEC
86       WRITE(0,83)
83       FORMAT(/' Input data for deviations from the average energy'/
     :   ' First FK integrals, then GK integrals in indicated format'/
     :   '  cc.ccccccccccFkk(el1,el2)  - terminate each list with an *',
     :   ' in the F column')
         FIRST = .TRUE.
*
*  *****  READ 'FK' AND 'GK' CARDS, OMITTING THE HEADER IF A FILE
*
82       READ (5,84) CFG,W,KFG,EL1,EL2
84       FORMAT(F14.8,A1,I2,1X,A3,1X,A3)
         IF ( W .NE. ASTER ) THEN
            CALL EPTR(EL,EL1, IFG, *99)
            CALL EPTR(EL,EL2, JFG, *99)
            CALL ADD(CFG,KFG,IFG,JFG,FIRST)
            CALL ADD(CFG,KFG,JFG,IFG,FIRST)
            GO TO 82
         ELSE IF (FIRST) THEN
            FIRST = .FALSE.
            GO TO 82
         END IF
      END IF
*
*  *****  COMPUTE THE INITIAL ARRAY AND INITIAL RADIAL FUNCTIONS
*
      CALL WAVEFN
*
*      ... Define an order for the functions to be iterated
*
      DO 90 JP = 1,NWF
         IORD(JP) = JP
90    CONTINUE
91    ORDERD = .TRUE.
      DO 92 JP = IB,NWF-1
         N1 = N(IORD(JP))
         L1 = L(IORD(JP))
         N2 = N(IORD(JP+1))
         L2 = L(IORD(JP+1))
         IF (N1.GT.N2 .OR. (N1.EQ.N2 .AND. L1.GT.L2)) THEN
            ITEMP = IORD(JP)
            IORD(JP) = IORD(JP+1)
            IORD(JP+1) = ITEMP
            ORDERD = .FALSE.
         END IF
92    CONTINUE
      IF (.NOT. ORDERD) GO TO 91
      RETURN
99    STOP
      END
*
*     ------------------------------------------------------------------
*                       D E
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
      SUBROUTINE DE(I1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      CHARACTER CONFIG*50,EL*3,ATOM*6,TERM*6
      LOGICAL FAIL,OMIT,REL,ALL,TRACE,CHANGE
      COMMON /TEST/FAIL,OMIT,REL,ALL,TRACE
      COMMON /COEFF/COEF(100),IJPTR(5,5)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :   S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
      COMMON /LABEL/CONFIG,EL(NWFD),ATOM,TERM
      COMMON P2(NOD),HQ(NOD),XX(NOD),V,B4,CN,C,XY,XP,
     :     AZZ,PP,FN,EM,FM,EU,FU,DELTAE,M,NODE,MK,KK,NJ
      CHARACTER*3 ASTER(3)
      DATA ASTER(1),ASTER(2),ASTER(3)/'  ','* ','**'/
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
*  *****  CHECK IF DIFFERENT METHOD SHOULD BE USED
*
      IF ( KK .EQ. 1 ) THEN
         IF (DABS(D1 -ED2/E(I,I)) .LT. 0.005D0  .AND.
     :       DMAX1(DABS(D1 - PN), DABS(D1/PN - D1)) .GT. 0.20D0 ) THEN
            METH(I) = 2
            KK = 2
            GO TO 25
         END IF
      ELSE
         IF (DABS(D1 - ED2/E(I,I)) .LT. 0.0001D0 .AND.
     :       IC .GT. 1) IC = IC -1
      END IF
*
*  *****  SET THE ACCELERATING PARAMETER
*
13    IF (IPR .NE. I) THEN
         ACC(I) = .75*ACC(I)
      ELSE
         ED2 = ED2 - E(I,I)
         print *, "This could uses uninitialized variable ED1."
         print *, "Aborting for now, until it gets fixed..."
         stop 1
         IF (ED1*ED2 .GT. D0 ) THEN
            ACC(I) = .75*ACC(I)
         ELSE
            ACC(I) = (D1 + D3*ACC(I))/D4
         END IF
      END IF
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
21    P(J,I) = PDE(J) + C*DIFF
      IF (M .EQ. NO) GO TO 26
      M = M + 1
      DO 24 J = M,NO
24    P(J,I) = D0
      AZ(I) = CD*AZZ + C*AZ(I)
      AZZ = AZ(I)
*
*  *****  CHECK THE ORTHOGONALIZATION
*
26    NN = NWF
      IF (OMIT) NN = IB - 1
      IJ = 0
      DPW = DP/DSQRT(SUM(I))
      M = MAX(I)
      CHANGE = .FALSE.
      DO 60 J = 1,NN
         IF (E(I,J) .NE. D0 .AND. I .NE. J ) THEN
          IF (DPM(J) .LT. DSQRT(SUM(J))*DPW .OR. J .LT. IB ) THEN
*
*        ORTHOGONALITY CONDITION APPLIES
*
            C = QUADR(I,J,0)
            WRITE(6,63) EL(J),EL(I),C
63          FORMAT(6X,'<',A3,'|',A3,'>=',1PD8.1)
            M = MAX0(M,MAX(J))
            DO 64 JJ = 1,M
               P(JJ,I) = P(JJ,I) - C*P(JJ,J)
64          CONTINUE
            AZZ = AZZ - C*AZ(J)
          CHANGE = .TRUE.
          END IF
       END IF
60    CONTINUE
      IF (CHANGE .OR. C.NE.D0) THEN
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
      WRITE (6,17) EL(I),E(I,I),AZ(I),PN,ASTER(KK),DP
17    FORMAT(20X,A3,2F15.7,F12.7, A2,1PD10.2)
      DPM(I) = DP
      IF (IPR .EQ. I1) THEN
         ED1 = ED2
      ELSE
         ED1 = ED2 - E(I1,I1)
      END IF
      IPR = I1
      RETURN
*
*  *****  IF METHD1 FAILED TO FIND AN ACCEPTABLE SOLUTION, ORTHOGONALIZE
*  *****  THE ESTIMATES AND TRY AGAIN
*
25    IF (I .EQ. IB) GO TO 27
      CALL ORTHOG
      CALL GRANGE
27    CALL METHD1(I)
      IF ( FAIL ) THEN
*
*  *****  ERROR RETURN FROM SECOND TRY.  IF M1 WAS USED,SWITCH TO
*         M2 AND TRY ONCE MORE.
*
         IF ( KK .EQ. 2) RETURN
         KK = 2
         GO TO 27
      ELSE
         GO TO 12
      END IF
      END
*
*     ------------------------------------------------------------------
*                       D E V
*     ------------------------------------------------------------------
*
*     Add the deviations to the average energy for a partially filled
*       p- or d- shell
*
      SUBROUTINE DEV(IEL, L, Q, I, DONE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      LOGICAL DONE
      INTEGER F2PP(6), F2DD(45), F4DD(45)
      DATA    F2PP/-3,3,12,-9,0,6/
      DATA    F2DD/-58,77,50,-13,140,
*             ... d3 coefficients
     :       -93,42,-12,-57,123,105,69,-12,
*             ... d4 coefficients
     :        -105,-69,-24,66,12,39,21,57,
     :        -51,30,48,84,219,111,210,138,
*             ... d5 coefficients
     :        -175,-85,23,-22,-112,-76,-58,167,
     :        23,-85,59,140,104,86,320,113/
      DATA    F4DD/5,-70,15,50,140,
*             ... d3 coefficients
     :        -30,-105,30,55,-45,105,-15,30,
*             ... d4 coefficients
     :        -105,15,-10,45,-30,-45,70,-55,
     :        75,135,20, 0,30,-15,210,-30,
*             ... d5 coefficients
     :        -175,-50,-40,-85,35,50,110,-15,
     :        -5,125,-25,140,20,-40,-100,-55/

      DONE = .TRUE.
      N = int(Q)
      IF (N .GT. 2*L+1) N = 4*L+2 - N
      IF (N .GT. 1) THEN
         IF (L .EQ. 1) THEN
            CALL ADD(2*F2PP(I)/25.D0,2,IEL,IEL,.TRUE.)
         ELSE IF (L .EQ. 2) THEN
            I = I-6
            CALL ADD(2*F2DD(I)/441.D0,2,IEL,IEL,.TRUE.)
            CALL ADD(2*F4DD(I)/441.D0,4,IEL,IEL,.TRUE.)
CGG
         ELSE IF (L .EQ. 3) THEN
            CALL DEVF(IEL, L, N, I, DONE)
CGG
         ELSE
            DONE = .FALSE.
         END IF
      END IF
      RETURN
      END
*
*     ------------------------------------------------------------------
*                              D Y K
*     ------------------------------------------------------------------
*
*       Stores in YK the values of the integral of
*              k
*       P (s/r) (dP /ds - P /s) integrated over the interval (0,r)
*        i         j       j
*
*   which enter into the spin-orbit calculation.
*
*
      SUBROUTINE DYK(I,J,K)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      DEN = L(I)+ L(J)+ 2 + K
      FACT = L(J)
      DO 1 JJ =1,2
1     YK(JJ) = FACT*P(JJ,I)*P(JJ,J)*R(JJ)/DEN
      A = EH**K
      AA = A*A
      A = D4*A
      F1 = FACT*P(1,I)*P(1,J)*R(1)
      F2 = FACT*P(2,I)*P(2,J)*R(2)
      DO 8 M =3,ND
         F3 = (-P(M+2,J) + D8*(P(M+1,J)-P(M-1,J)) + P(M-2,J))/(D6*H)
         F3 = D5*P(M,I)*(F3 - P(M,J))*R(M)
         YK(M) = YK(M-2)*AA + H3*(F3+ A*F2 + AA*F1)
         F1 = F2
         F2 = F3
8     CONTINUE
      A = A*(EH)**3
      AA = A*A/D16
      C = 2*K+3
      HH = C*H3
      YK(NO)= YK(ND)
      F1 = YK(NO)
      F2 = F1
      DO 9 MM = 3,NO
         M = NO -MM+1
         F3 =YK(M)
         YK(M) = YK(M+2)*AA + HH*(F3 +A*F2 + AA*F1)
         F1 = F2
         F2 = F3
9     CONTINUE
      RETURN
      END
*
*     ------------------------------------------------------------------
*                       E K I N
*     ------------------------------------------------------------------
*
*       Returns the value of the integral of
*
*         (2/r)P (Y P  + X )
*               j  i i    i
*
*   integrated with respect to r.
*
*
      DOUBLE PRECISION FUNCTION EKIN(I,II,REL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      LOGICAL REL
      CALL XCH(I,2)
      CALL POTL(I,REL)
      DO 1 J=1,NO
         YK(J) = YR(J)
         YR(J) = P(J,II)
1     CONTINUE
      EKIN = D2*QUADS(I,II,1) + QUAD(II,NO ,YR,X)
      RETURN
      END
*
*     ------------------------------------------------------------------
*                       E N E R G Y
*     ------------------------------------------------------------------
*
*       Determines the position of the electron in the electron list
*
      SUBROUTINE ENERGY(ETOTAL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      LOGICAL FAIL,OMIT,REL,ALL,TRACE
      COMMON /TEST/FAIL,OMIT,REL,ALL,TRACE
      CHARACTER CONFIG*50,EL*3,ATOM*6,TERM*6
      COMMON /LABEL/CONFIG,EL(NWFD),ATOM,TERM
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :   S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
*
*  *****  COMPUTE KINETIC ENERGY IF NECESSARY
*
      DO 50 I = 1,NWF
50    EK(I) = -D5*HL(EL,I,I,REL)
*
      ETOTAL = D0
      DO 10 I = 1,NWF
         ETOTAL = ETOTAL + SUM(I)*(EK(I))
         DO 12 J = 1,I
            DO 11 K = 0,2*MIN0(L(I),L(J)),2
               C = A(I,J,K)*SUM(I)
               IF (I .EQ. J) C = C/D2
               IF (ABS(C).NE.D0) ETOTAL = ETOTAL + C*FK(I,J,K,REL)
 11         CONTINUE
 12      CONTINUE
         DO 15 J = 1,I-1
            DO 14 K = ABS(L(I)-L(J)),L(I)+L(J),2
              C = B(I,J,K)*SUM(I)
              IF (ABS(C).NE.D0) ETOTAL=ETOTAL+C*GK(I,J,K,REL)
 14         CONTINUE
 15      CONTINUE
 10     CONTINUE
      END
*
*     ------------------------------------------------------------------
*                       E N E X P R
*     ------------------------------------------------------------------
*
*     Determine the deviations to the average energy for the following:
*        i) an open p- , d- or f-shell
*       ii) a single electron or hole, any l
*      iii) an s-electron and a single electron, any l
*       iv) an s-electron and an open p- , d- or f-shell
*        v) an open p-shell and a single electron, any l
*       vi) a single electron in f- and a single electron in d- shells
*
      SUBROUTINE ENEXPR(TERM, DONE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      CHARACTER TERM*6,SL*2,SENOR*1,PSL*2,SLM*2,SLP*2
      LOGICAL DONE
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :   S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      INTEGER SUMTAB(5),PARTAB(11),PTRTAB(11),LTAB(54),
     :   NOS(2),PLVAL(11),PACVAL,SP,PS1,PS2
*
*     ... FINT, GINT1, and GINT2 are coefficients of polynomials
*         in l, tabulated by Slater,
*
*
      INTEGER FINT(3,54),GINT1(3,54),GINT2(3,54)
      CHARACTER*1 PARCH(11)
*
*     ... coefficients of F2 integrals for p(n)l(1) configurations
*
      DATA FINT/2,-1,0,-4,-4,3,2,5,3,2,-1,0,-4,-4,3,2,5,3,
     :   -2,1,0,4,4,-3,-2,-5,-3,-2,1,0,4,4,-3,-2,-5,-3,
     :   4,-2,0,-2,-11,6,-4,-4,15,-2,7,15,4,10,6,0,0,0,
     :   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     :   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     :   2,-1,0,-4,-4,3,2,5,3,2,-1,0,-4,-4,3,2,5,3,
     :   -4,2,0,2,11,-6,4,4,-15,2,-7,-15,-4,-10,-6,0,0,0,
     :   -2,1,0,4,4,-3,-2,-5,-3,-2,1,0,4,4,-3,-2,-5,-3/
*
*     ... coefficients of G(l-1) integrals
*
      DATA GINT1/-10,5,0,2,11,-6,2,-1,-6,14,-7,0,2,-13,6,2,-1,6,
     :   -8,4,0,-8,4,0,4,10,0,10,-5,0,10,-5,0,4,-8,0,
     :   -8,4,0,-2,13,-6,2,11,-12,4,4,-12,4,-2,-12,0,0,0,-6,3,0,10,-5,0,
     :   -6,3,0,-6,3,0,-4,8,3,0,12,3,6,9,0,-6,3,0,6,21,12,12,12,-27,
     :   12,-6,27,6,-15,-12,-6,3,0,0,6,-3,0,0,-3,6,-3,0,0,-6,3,12,6,3,
     :   -4,2,0,-4,2,0,-4,2,0,-4,2,0,14,11,-9,14,-7,-9,-4,-2,0,-4,2,0,
     :   -2,7,3,2,11,3,8,8,0,0,0,0,-2,1,0,-2,1,0,-2,1,0,-2,1,0,-2,1,0,
     :   22,13,0/
*
*     ... coefficients of G(l+1) integrals
*
      DATA GINT2/2,5,-3,2,-7,-15,-10,-25,-15,2,5,9,2,17,21,14,35,21,
     :   4,-2,-6,-8,-20,-12,-8,-20,-12,4,16,12,10,25,15,10,25,15,
     :   4,10,0,4,4,-12,2,-7,-21,-2,-17,-21,-8,-20,-12,0,0,0,
     :   -6,-15,-9,10,25,15,6,3,-3,0,-12,-9,-4,-16,-9,-6,-15,-9,
     :   -6,-15,-9,6,27,9,12,30,-9,12,12,-27,6,-9,-27,-6,-15,-9,
     :   0,0,-3,0,-6,-9,-6,-15,-9,12,18,9,0,6,9,6,15,9,
     :   -4,-10,-6,-4,-10,-6,-4,-10,-6,14,35,12,14,17,-6,-4,-10,-6,
     :   8,8,0,2,-7,-6,-2,-11,-6,-4,-10,-6,-4,-10,-6,0,0,0,
     :   -2,-5,-3,-2,-5,-3,-2,-5,-3,22,31,9,-2,-5,-3,-2,-5,-3/
*
*     ... Encoded term value -- S = LTAB/10
*                         Lterm = L + (LTAB mod 10 - 5)
*         Example:  LTAB = 36 with L = 2  is 3F
*
      DATA LTAB/36,35,34,16,15,14,46,45,44,26,25,24,27,26,25,24,
     :   23,25,55,35,37,36,35,34,33,17,16,15,14,13,36,35,34,16,15,14,
     :   46,45,44,26,25,24,27,26,25,24,23,25,36,35,34,16,15,14/
      DATA SUMTAB/1,4,7,10,11/
      DATA PARTAB/2,3,1,1,4,2,2,3,1,1,2/
      DATA PTRTAB/6,12,17,18,20,30,36,42,47,48,54/
      DATA PLVAL/1,1,2,0,0,2,1,1,2,0,1/
      DATA PARCH/'P','P','D','S','S','D','P','P','D','S','P'/
*
      IP = 1
 1    IF (TERM(IP:IP) .EQ. ' ') THEN
         IP = IP+1
         GO TO 1
      END IF
      SL = TERM(IP:IP+1)
      SENOR = ' '
      IF (IP.LE.4) SENOR = TERM(IP+2:IP+2)
*
*   ---  convert lowercase L symbol to uppercase
*
      IF (SL(2:2).GT.'a' .AND. SL(2:2).LT.'z')
     :    SL(2:2) = CHAR(ICHAR(SL(2:2)) + ICHAR('A') - ICHAR('a'))
*
*  ---  determine if FK or GK data needs to be input
*
      IL = 0
      IS = 0
      J = 1
      DO 2 I = NCLOSD+1, NWF
         IF (SUM(I) .NE. 4*L(I)+2 .AND. SUM(I) .NE. 0.D0) THEN
            IF (J.GT.2) THEN
               IF (SL .NE. 'AV' .AND. SL .NE. 'aV') THEN
                  DONE=.FALSE.
               ELSE
                  DONE = .TRUE.
               END IF
               RETURN
            ENDIF
            NOS(J) = I
            J = J + 1
            IF (L(I) .EQ. 0 .AND. IS .EQ. 0) THEN
               IS = IS + 1
               IIS = I
            ELSE
               IL = IL + 1
               IIL = I
            END IF
         END IF
 2    CONTINUE
      IF (SL .NE. 'AV' .AND. SL .NE. 'aV' .AND. IS+IL.NE.0) THEN
         DONE = .FALSE.
         C = 0.D0
         IF (IS+IL .LE. 2 .AND. IL .LE. 1) THEN
            IF (IS .EQ. 0 .AND. IL .EQ. 1) THEN
3             CALL LOOKTM(L(IIL),SL,SENOR,SUM(IIL),IP,NSL)
              IF (NSL .GT. 1) THEN
                IF (L(IIL) .NE. 3) THEN
                  WRITE(0,*)' Ambiguous term: enter seniority'
                ELSE
                  WRITE(0,*)' Ambiguous term: enter Nr for f-subshells'
                END IF
                READ (5,'(A1)') SENOR
                GO TO 3
              END IF
              CALL DEV(IIL,L(IIL),SUM(IIL),IP,DONE)
            ELSE IF (IS .EQ. 1 .AND. IL .EQ. 1) THEN
              SLM = SL
              SLP = SL
              SLM(1:1) = CHAR(ICHAR(SLM(1:1)) - 1)
              SLP(1:1) = CHAR(ICHAR(SLP(1:1)) + 1)
              CALL LOOKTM(L(IIL),SLM,SENOR,SUM(IIL),IPM,NSLM)
              CALL LOOKTM(L(IIL),SLP,SENOR,SUM(IIL),IPP,NSLP)
              IF (NSLM+NSLP .EQ. 0) THEN
                DONE = .FALSE.
                RETURN
              ELSE IF (NSLM .EQ. 1 .AND. NSLP .EQ. 0) THEN
                SL = SLM
                IP = IPM
              ELSE IF (NSLM .EQ. 0 .AND. NSLP .EQ. 1) THEN
                SL = SLP
                IP = IPP
              ELSE IF (NSLM .EQ. 1 .AND. NSLP .EQ. 1) THEN
4               WRITE(0,'(A,A3,A,A3)')' Ambiguous l**n term: enter',
     :                   SLM, ' or ',SLP
                READ(5,'(A2)') SL
                IF (SL .EQ. SLM) THEN
                  IP = IPM
                ELSE IF (SL .EQ. SLP) THEN
                  IP = IPP
                ELSE
                  WRITE(0,*) ' Term not allowed: re-enter'
                  GO TO 4
                END IF
              ELSE
5             CONTINUE
              IF (L(IIL) .NE. 3) THEN
                WRITE(0,'(A,A)') ' Ambiguous l**n parent term:',
     :                   'Enter term and seniority'
              ELSE
                WRITE(0,'(A,A)') ' Ambiguous l**n parent term:',
     :                   'Enter term and Nr for f-subshells'
              END IF
              READ(5,'(A2,A1)') SL, SENOR
              IF (SENOR .EQ. ' ') THEN
                IF (L(IIL) .NE. 3) THEN
                  WRITE(0,*) 'Seniority is needed'
                ELSE
                  WRITE(0,*) 'Nr for f-subshells is needed'
                END IF
                GO TO 5
              END IF
              CALL LOOKTM(L(IIL),SL,SENOR,SUM(IIL),IP,NSL)
              IF (NSL .NE. 1) THEN
                     WRITE(0,'(A,A3,A,A3,A)') ' Allowed terms are ',
     :                    SLM, ' or ', SLP,' plus seniority'
                     GO TO 5
                  END IF
               END IF
               CALL DEV(IIL,L(IIL),SUM(IIL),IP,DONE)
               IF (DONE ) THEN
                  SP = ICHAR(SL(1:1)) - ICHAR('0')
                  CSP = (SP - 1)/2.
                  IF (SL .EQ. SLM) THEN
                     C = -CSP/(2*L(IIL)+1)
                  ELSE
                     C = (CSP + 1)/(2*L(IIL)+1)
                  END IF
                  CALL ADD(C,L(IIL),IIS,IIL,.FALSE.)
                  CALL ADD(C,L(IIL),IIL,IIS,.FALSE.)
               END IF
            ELSE IF (IS .EQ. 1 .AND. IL .EQ. 0) THEN
               DONE = .TRUE.
            END IF
CGG
CGG
         ELSEIF ((L(NOS(1)).EQ.2).AND.(L(NOS(2)).EQ.3).OR.
     :          (L(NOS(1)).EQ.3).AND.(L(NOS(2)).EQ.2)) THEN
           IF (SUM(NOS(1)).EQ.1.AND.SUM(NOS(2)).EQ.1.D0) THEN
   11          CALL LOOKTMDF(SL,IP,NSL)
               IF (NSL .EQ. 0) THEN
                  WRITE(0,*)' Re-enter the term of configuration'
                  READ (5,'(A2)') SL
                  GO TO 11
               END IF
               IF (L(NOS(1)).EQ.2) THEN
                 NGGD=NOS(1)
                 NGGF=NOS(2)
               ELSE
                 NGGD=NOS(2)
                 NGGF=NOS(1)
               ENDIF
               CALL DEVDF(NGGD, NGGF, IP, DONE)
           ENDIF
CGG
CGG
         ELSE
            IF (((L(NOS(1)).EQ.1).AND.(SUM(NOS(2)).EQ.1.D0)).OR.
     :          ((L(NOS(2)).EQ.1).AND.(SUM(NOS(1)).EQ.1.D0))) THEN
               IF (L(NOS(1)).EQ.1.AND.SUM(NOS(2)).EQ.1.D0) THEN
                  ISUMP=int(SUM(NOS(1)))
                  NP = NOS(1)
                  NL = NOS(2)
               ELSE
                  ISUMP = int(SUM(NOS(2)))
                  NP = NOS(2)
                  NL = NOS(1)
               ENDIF
               SP=ICHAR(SL(1:1))-ICHAR('0')
               LP=LVAL(SL(2:2))
               PS1=SP+1
               PS2=SP-1
               IF (ISUMP.EQ.1) THEN
                  IPTR1=1
               ELSE
                  IPTR1=SUMTAB(ISUMP-1)+1
               END IF
               IPTR2=SUMTAB(ISUMP)
               NOMACH=0
               CALL LOOKUP(PARTAB,IPTR1,IPTR2,IND,NOMACH,PS1)
               CALL LOOKUP(PARTAB,IPTR1,IPTR2,IND,NOMACH,PS2)
               PSL(1:1)=CHAR(PARTAB(IND)+ ICHAR('0'))
               PSL(2:2)=PARCH(IND)
               IF (NOMACH.GT.1) THEN
                  WRITE(0,*)' AMBIGUOUS PARENT CASE'
 10               WRITE(0,*)' ENTER THE SL TERM FOR p(n) SUBSHELL'
                  READ (5,'(A)')PSL
                  IF (PSL(2:2).GT.'a'.AND.PSL(2:2).LT.'z')
     :                   PSL(2:2)=CHAR(ICHAR(PSL(2:2))+ICHAR('A')
     :                                                -ICHAR('a'))
                  PS1=ICHAR(PSL(1:1))-ICHAR('0')
                  PS2=LVAL(PSL(2:2))
                  CALL LOOKUP(PLVAL,IPTR1,IPTR2,IND,NOMACH,PS2)
                  IF ((NOMACH.NE.1).AND.(PARTAB(IND).NE.PS1))
     :                    GO TO 10
               END IF
               IF (ISUMP.EQ.1) THEN
                  IPTR1=1
               ELSE
                  IPTR1=PTRTAB(IND-1)+1
               END IF
               IPTR2 = PTRTAB(IND)
               LV=L(NL)
               PACVAL=SP*10+LP-LV+5
               NOMACH=0
               CALL LOOKUP(LTAB,IPTR1,IPTR2,IND,NOMACH,PACVAL)
               IF (NOMACH.NE.1) THEN
                  DONE=.FALSE.
                  RETURN
               ENDIF
               VAL1=((FINT(1,IND)*LV+FINT(2,IND))*LV+FINT(3,IND))
     :                  /(5.D0*(2*LV-1)*(2*LV+3))
               VAL2=((GINT1(1,IND)*LV+GINT1(2,IND))*LV+GINT1(3,IND))
     :                  /(2.D0*(2*LV+1)*(2*LV-1)**2)
               VAL3=((GINT2(1,IND)*LV+GINT2(2,IND))*LV+GINT2(3,IND))
     :                  /(2.D0*(2*LV+1)*(2*LV+3)**2)
*
*     ...  Add contributions from between p-subshell and l-electron
*
               CALL ADD(VAL1,2,NP,NL,.TRUE.)
               CALL ADD(VAL1,2,NL,NP,.TRUE.)
               CALL ADD(VAL2,LV-1,NP,NL,.FALSE.)
               CALL ADD(VAL2,LV-1,NL,NP,.FALSE.)
               CALL ADD(VAL3,LV+1,NP,NL,.FALSE.)
               CALL ADD(VAL3,LV+1,NL,NP,.FALSE.)
*
*     ... Add deviations for p-subshell
*
               CALL LOOKTM(1,PSL,' ',SUM(NP),IP,NSL)
               CALL DEV(NP,1,SUM(NP),IP,DONE)
            ELSE
               DONE = .FALSE.
            END IF
         END IF
      ELSE
         DONE = .TRUE.
      END IF
      RETURN
      END
*
*     ------------------------------------------------------------------
*                       E P T R
*     ------------------------------------------------------------------
*
*       Determines the position of the electron in the electron list
*
      SUBROUTINE EPTR(EL,ELSYMB, IEL, *)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER EL(*)*3, ELSYMB*3, BL*3
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      DATA BL/'   '/
*
* ***** SEARCH ELECTRON LIST FOR LSYMB
*
      IF ( ELSYMB .EQ. BL ) THEN
         IEL = 0
         RETURN
      ENDIF
      DO 10 I=1,NWF
         IF (EL(I) .EQ. ELSYMB ) THEN
            IEL = I
            RETURN
         ENDIF
10    CONTINUE
      IEL = -1
      WRITE(0,20) ELSYMB
20    FORMAT(/10X,A3,' NOT FOUND IN ELECTRON LIST')
      RETURN 1
      END
*
*     -----------------------------------------------------------------
*           F A C T R L
*     -----------------------------------------------------------------
*
*
      SUBROUTINE FACTRL(NFACT)
*
*      GAM(I) = LOG( GAMMA(I-1) ), WHERE GAMMA(I) = FACTORIAL I-1
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      COMMON/FACT/GAM(100)
      DATA ZERO,ONE,TWO/0.D0,1.D0,2.D0/
*
      GAMMA=ONE
      GAM(1) = ZERO
      DO 1 I=1,NFACT-1
         GAMMA=I*GAMMA
         GAM(I+1) = DLOG(GAMMA)
    1 CONTINUE
      DO 20 I = NFACT+1,(100)
         X = I-1
         GAM(I) = GAM(I-1) + DLOG(X)
   20 CONTINUE
      RETURN
      END
*
*     ------------------------------------------------------------------
*                 F K
*     ------------------------------------------------------------------
*                             k
*       Returns the value of F (i,j)
*
*
      DOUBLE PRECISION FUNCTION FK(I,J,K,REL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL REL
      CALL YKF(I,I,K,REL)
      FK = QUADS(J,J,1)
      RETURN
      END
*
*     ------------------------------------------------------------------
*                 G K
*     ------------------------------------------------------------------
*                             k
*       Returns the value of G (i,j).
*
*
      DOUBLE PRECISION FUNCTION GK(I,J,K,REL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL REL
      CALL YKF(I,J,K,REL)
      GK = QUADS(I,J,1)
      RETURN
      END
*
*     ------------------------------------------------------------------
*               G R A N G E
*     ------------------------------------------------------------------
*
*       Controls the calculation of off-diagonal energy parameters.
*   It searches for all pairs (i,j) which are constrained through  an
*   orthogonality requirement.   Eq. (7-10) is used to calculate the
*   parameter.
*
      SUBROUTINE GRANGE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      CHARACTER CONFIG*50,EL*3,ATOM*6,TERM*6
      LOGICAL FAIL,OMIT,REL,ALL,TRACE
      COMMON /TEST/FAIL,OMIT,REL,ALL,TRACE
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :   S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
      COMMON /LABEL/CONFIG,EL(NWFD),ATOM,TERM
      COMMON /COEFF/COEF(100),IJPTR(5,5)
*
*  *****  ROTATE PAIRS CONNECTED BY ORTHOGONALITY BUT NOT WHEN ONE OF
*         THE ORBITALS IS SIMULTANEOUSLY ORTHOGONAL TO A NON-ORTHOGONAL
*         PAIR
*
      DO 1 I = IB,NWF-1
         DO 2 J = I+1,NWF
            IF (DABS(E(I,J)) .GT. 1.D-10) CALL ROTATE(I,J)
2        CONTINUE
1     CONTINUE
*
*  *****   COMPUTE OFF-DIAGONAL ENERGY PARAMETERS
*
      DO 10 I = MAX0(2,IB),NWF
         DO 12 J = 1,I-1
           IF (DABS(E(I,J)) .GT. 1.D-10) THEN
        IF ( J .LT. IB) THEN
           E(I,J) = HL(EL,I,J,REL) - EKIN(I,J,REL)
           E(J,I) = D0
        ELSE IF (SUM(I) .EQ. SUM(J)) THEN
           C=HL(EL,I,J,REL)-(EKIN(I,J,REL)+EKIN(J,I,REL))/D2
           E(I,J) = C
           E(J,I) = C
        ELSE
           RES = D0
           DO 14 II = 1,NWF
             IF (II.EQ.I .OR. II.EQ.J) THEN
                DO 22 K = 0,2*L(I),2
                  IF(II.EQ.I) THEN
                    C= A(I,I,K)-A(J,I,K)-B(J,I,K)
                    IF (DABS(C).GT.1.D-10)
     :                          RES = RES + C*RK(I,I,I,J,K,REL)
                  ELSE IF (II.EQ.J) THEN
                    C= A(J,J,K)-A(I,J,K)-B(I,J,K)
                    IF (DABS(C).GT.1.D-10)
     :                         RES = RES - C*RK(J,J,J,I,K,REL)
                  END IF
22              CONTINUE
             ELSE
                DO 24 K = 0,2*MIN0(L(I),L(II)),2
                    C = A(I,II,K) - A(J,II,K)
                    IF (DABS(C).GT.1.D-10)
     :                          RES = RES + C*RK(I,II,J,II,K,REL)
                    KK = ABS(L(I)-L(II)) + K
                    C = B(I,II,KK) - B(J,II,KK)
                    IF(DABS(C) .GT. 1.D-10)
     :                            RES = RES + C*RK(I,II,II,J,KK,REL)
24              CONTINUE
             END IF
 14        CONTINUE
           E(I,J) = D2*SUM(J)*RES/(SUM(I)-SUM(J))
           E(J,I) = SUM(I)*E(I,J)/SUM(J)
              END IF
           END IF
           IF (DABS(E(I,J)) .GT. 1.D-10) THEN
              WRITE(6,35) EL(I),EL(J),E(I,J),EL(J),EL(I),E(J,I)
35            FORMAT(7X,2(3X,'E(',2A3,') =',F12.5))
           END IF
 12      CONTINUE
 10   CONTINUE
      RETURN
      END

*
*     ------------------------------------------------------------------
*               H E L P
*     ------------------------------------------------------------------
*
*       Provide HELP information about the data requested
*
      SUBROUTINE HELP(CASE)
      INTEGER CASE
*
      GO TO (10,20,30,40,50,60,70,80,90,100) CASE
*
10    WRITE(0,11)
*  ***** Which orbitals varied?
11    FORMAT(//1X,'Response  ALL will vary all orbitals'/
     :         11X,'NONE will not vary any orbitals'/
     :         11X,'=n (integer n), will vary last n'/
     :         11X,'with comma delimited list will vary',
     :         1X,'only the orbitals in the list'//)
      RETURN

20    WRITE(0,21)
*  ***** Default electron parameters ?
21    FORMAT(//1X,
     :  'Response  N  will prompt the user for:'//
     :   7X,'S     : Screening parameter  (Real number) '/7X,
     :  'IND   : Indicator specifying the type of initial estimate'/
     :  15X,'0 - Screened hydrogenic functions'/
     :  14X,'-1 - Search for functions in wavefunction file; if'/
     :  20X,'not present use screened hydrogenic.'/
     :  7X,'METH  : Method for solving the differential equation'/
     :  15X,'1 - Method 1 solves the boundary value problem for an'
     :  /20X,'acceptable solution which need not be normalized.'/
     :  15X,'2 - Method 2 solves the boundary value problem for '/
     :  20X,'an acceptable solution which is normalized to first'/
     :  20X,'order. If the exchange function is identically zero the '/
     :  20X,'program will automatically select Method 2.'/
     :  15X,'3 - Method 3 is similar to Method 1 but omits all checks'/
     :  20X,'for acceptability.'/
     :  7X,'ACC   : Inital accelerating factor '/
     :  18X,'( Real number such that 0 .LE. ACC .LT. 1 )'//)
      RETURN

30    WRITE(0,31)
*  ***** Default values (NO,STRONG) ?
31    FORMAT(//1X,
     :  'Response  Y will set default values as follows'//
     :  10X,'NO - Maximum number of points in the range of the '/
     :  15X,'function is set to 200.'/
     :  10X,'Strong - is logically set to .FALSE.'//
     :  1X,'Response N  will prompt the user for'//
     :  10X,'NO - Maximum number of points in the range of '/
     :  15X,'the function which should be a positive integer'/
     :  15X,'from 160 for a small atom to 220 for a large atom.'/
     :  10X,'Strong - may be set to .TRUE. or .FALSE by user.'//)
      RETURN

40    WRITE(0,41)
*  ***** Default values for remaining parameters ?
41    FORMAT(//1X,
     :  1X,'Response  Y --Sets the following default values'//
     :  10X,'PRINT=.FALSE.'/10X,'SCFTOL=1.D-8'/10X,'NSCF=12'/
     :  10X,'IC=2 + (NWF + 1 - IB)/4'/10X,'TRACE=.FALSE.'//
     :  1X,'Response  N --prompts user for new parameter values.'//)
      RETURN

50    WRITE(0,51)
*  ***** Default values for PRINT, SCFTOL ?
51    FORMAT(//1X,
     :  1X,'Response  Y --Default value for PRINT is .FALSE. thus'/
     :  16X,'radial functions are NOT printed.'//
     :  16X,'SCFTOL -the initial value of the parameter defining the'/
     :  16X,'self-consistency tolerance for radial functions is set'/
     :  16X,'to a default value of  1.D-8 .'//
     :  1X,'Response  N --prompts user for new values of PRINT and '
     :  'SCFTOL .'//)
      RETURN

60    WRITE(0,61)
*  ***** Default values for NSCF,  IC ?
61    FORMAT(//1X,
     :  1X,'Response  Y --NSCF, the maximum number of cycles for'/
     :  16X,'the SCF process is set to default value of 12 and'/
     :  'IC is set to 2 + (NWF + 1 -IB)/4 .'//
     :  1X,'Response  N --user prompted for new NSCF and IC values'//)
      RETURN

70    WRITE(0,71)
*  ***** Default values for TRACE ?
71    FORMAT(//1X,
     :  'Response  Y --sets trace to default value of .FALSE. thus a'/
     :  16X,'trace of energy adjustment will  NOT  be printed.'//
     :  1X,'Response  N --prompts the user for the new value of trace.'/
     :  16X,'If .TRUE. a trace will be printed showing the energy ',
     :  'adjustment'/16X,'process used by METHD1 for finding an',
     :  ' acceptable solution'/16X,'with the correct number of',
     :  ' nodes.'//)
      RETURN

80    WRITE(0,81)
*  ***** Additional parameters ?
81    FORMAT(//1X,
     :  'Response  Y --additional values may be computed for '/
     :  16X,'SLATER OR MAGNETIC INTEGRALS'/
     :  16X,'EXPECTATION VALUES OF R**K'/
     :  /16X,'ELECTRON DENSITY AT THE NUCLEUS'/16X,'SPIN-ORBIT '
     :  'PARAMETER'/16X,'TRANSITION INTEGRALS'//
     :  1X,'Response  N --additional parameter computation is ',
     :  'skipped'//)
      RETURN

90    WRITE(0,91)
*  ***** Do you wish to continue along the sequence ?
91    FORMAT(//1X,
     :  'Response  Y --sequence is continued'//
     :  1X,'Response  N --current case is ended, but new case may be'/
     :  16X,'started'//)
      RETURN

100   WRITE(0,101)
*  ***** Do you wish to continue ?
101   FORMAT(//1X,
     :  'Response  Y --prompts user for additional iterations and'/
     :  16X,'new IC.  Then performs additional NSCF self-consistent'/
     :  16X,'field iterations.'//
     :  1X,'Response  N --terminates the calculation.'//)
      RETURN
      END
*
*     ------------------------------------------------------------------
*                 H L
*     ------------------------------------------------------------------
*
*       Returns the value of <i^L^j>, using a special formula to
*  preserve symmetry.
*
      DOUBLE PRECISION FUNCTION HL(EL,I,J,REL)
      PARAMETER (NOD=220,NWFD=20)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER EL(*)*3
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      LOGICAL REL
      IF (IABS(L(I)-L(J)) .EQ. 0) GO TO 3
      WRITE(0,4) EL(I),L(I),EL(J),L(J)
4     FORMAT(10X,'UNALLOWED L VALUES OCCURRED IN HL SUBROUTINE'/
     :   2(10X,A3,' HAS L = ',I3))
      STOP
3     LI = L(I)
      C = 2*LI + 1
      A1 = -D2/(C*(LI+1))
      A2 = A1/((C+D2)*(LI+1))
      A3 = A2/((LI+2)*(LI+1))
      ZR = Z*R(1)
      HL = H*C*P(1,I)*P(1,J)*(D1+ZR*(A1+ZR*(A2+ZR*A3)))
      MM = MIN0(MAX(I)+3,MAX(J)+3,ND-1)
      K = 2
      C = D4/D3
      DI1 = P(K+1,I) - P(K-1,I)
      DI2 = P(K+1,I) - D2*P(K,I) + P(K-1,I)
      DJ1 = P(K+1,J) - P(K-1,J)
      DJ2 = P(K+1,J) - D2*P(K,J) + P(K-1,J)
      HL = HL + DI1*DJ1 + C*DI2*DJ2
      DO 1 K = 4,MM,2
      DI1 = P(K+1,I) - P(K-1,I)
      DI2 = P(K+1,I) - D2*P(K,I) + P(K-1,I)
      DI4 = P(K+2,I) - D4*(P(K+1,I)+P(K-1,I)) + D6*P(K,I) +P(K-2,I)
      DI3 = P(K+2,I) - P(K-2,I) - D2*DI1
      DI5 = P(K+3,I)-P(K-3,I) - D4*(P(K+2,I)-P(K-2,I))
     :   + 5.D0*(P(K+1,I)-P(K-1,I))
      DI6 = P(K+3,I)+P(K-3,I) - D6*(P(K+2,I)+P(K-2,I))
     :   + 15.D0*(P(K+1,I)+P(K-1,I)) - 20.D0*P(K,I)
      DJ1 = P(K+1,J) - P(K-1,J)
      DJ2 = P(K+1,J) - D2*P(K,J) + P(K-1,J)
      DJ4 = P(K+2,J) - D4*(P(K+1,J)+P(K-1,J)) + D6*P(K,J) +P(K-2,J)
      DJ3 = P(K+2,J) - P(K-2,J) - D2*DJ1
      DJ5 = P(K+3,J)-P(K-3,J) - D4*(P(K+2,J)-P(K-2,J))
     :   + 5.D0*(P(K+1,J)-P(K-1,J))
      DJ6 = P(K+3,J)+P(K-3,J) - D6*(P(K+2,J)+P(K-2,J))
     :   + 15.D0*(P(K+1,J)+P(K-1,J)) - 20.D0*P(K,J)
1     HL = HL + DI1*DJ1 + C*DI2*DJ2 + (DI3*DJ3 + DI2*DJ4+DI4*DJ2)/45.D0
     :  -(DI3*DJ5+DI5*DJ3)/252.D0 - (DI2*DJ6+DI6*DJ2-1.1*DI4*DJ4)/378.D0
      TZ = Z + Z
      C = (LI + D5)**2
      HL2 = D5*(TZ*R(1) - C)*P(1,I)*P(1,J)
      DO 2 K = 2,MM,2
2     HL2 = HL2 + D2*(TZ*R(K) - C)*P(K,I)*P(K,J)
     :   + (TZ*R(K+1) - C)*P(K+1,I)*P(K+1,J)
      HL = -HL/(D2*H) + HL2*H1
      IF (REL) HL=HL-D2*RLSHFT(I,J)
      RETURN
      END
*
*     ------------------------------------------------------------------
*               H N O R M
*     ------------------------------------------------------------------
*
*       Returns the value of the normalization constant for an (nl)
*   hydrogenic function with nuclear charge ZZ.
*
*
      DOUBLE PRECISION FUNCTION HNORM(N,L,ZZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      M = L + L + 1
      A = N + L
      B = M
      T = A
      D = B
      M = M - 1
      IF (M .EQ. 0) GO TO 2
      DO 1 I = 1,M
      A = A - D1
      B = B - D1
      T = T*A
1     D = D*B
2     HNORM = DSQRT(ZZ*T)/( N*D)
      RETURN
      END
*    MCHF_HF (Part2 of 2)
*     ------------------------------------------------------------------
*               H W F
*     ------------------------------------------------------------------
*
*       Returns the value of an unnormalized (nl) hydrogenic function
*   with nuclear charge ZZ and radius r.
*
*
      DOUBLE PRECISION FUNCTION HWF(N,L,ZZ,R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      K = N-L-1
      P = D1
      A = D1
      B = K
      C = N+ L
      X = -D2*ZZ*R/N
*
*  *****  TEST IF UNDERFLOW MAY OCCUR, IF SO SET HWF = 0
*
      IF ( X .LT. -150.D0 ) GO TO 5
      IF (K) 1,2,3
3     DO 4 I = 1,K
      P = D1 + A/B*P/C*X
      A = A + D1
      B = B - D1
4     C = C - D1
2     HWF = P*DEXP(X/D2)*(-X)**(L+1)
      RETURN
1     WRITE(0,7) N,L,ZZ,R
7     FORMAT(51H FORBIDDEN COMBINATION OF N AND L IN HWF SUBPROGRAM/
     :    4H N = ,I4,6H   L = ,I4,6H   Z = ,F6.1,6H   R = ,F8.4)
      STOP
5     HWF = D0
      RETURN
      END
*
*     ------------------------------------------------------------------
*               I N I T
*     ------------------------------------------------------------------
*
*       Initializes basic constants of the program including those
*   which define the average energy of a configuration.
*
*
      SUBROUTINE INIT
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /EAV/CCA(10),CCB(35)
      COMMON /FACT/GAM(100)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
*
*  *****  SET THE COMMONLY USED DOUBLE PRECISION CONSTANTS
*
      D0  =  0.D0
      D1  =  1.D0
      D2  =  2.D0
      D3  =  3.D0
      D4  =  4.D0
      D5  =  1.D0/2.D0
      D6  =  6.D0
      D8  =  8.D0
      D10 = 10.D0
      D12 = 12.D0
      D16 = 16.D0
      D30 = 30.D0
*
*  ***** Set the factorial needed by RME
*
      CALL FACTRL(32)
*
*  *****  SET FINE STRUCTURE CONSTANT
*
      FINE = 0.25D0/(137.036)**2
*
*  *****  SET THE STARTING POINT, STEP SIZE, AND RELATED PARAMETERS
*
      RHO = -4.D0
      H   = 1./16.D0
      H1 = H/1.5
      H3 = H/3.
      CH = H*H/12.
      EH = DEXP(-H)
      NO=220
      ND = NO - 2
*
*  *****  AVERAGE INTERACTIONS FOR EQUIVALENT ELECTRONS
*
*  *****  P - P
*
      CCA(1) = 2.D0/25.D0
*
*  *****  D - D
*
      CCA(2) = 2.D0/63.D0
      CCA(3) = 2.D0/63.D0
*
*  *****  F - F
*
      CCA(4) =   4.D0/ 195.D0
      CCA(5) =   2.D0/ 143.D0
      CCA(6) = 100.D0/5577.D0
*
*  *****  G - G
*
      CCA(7) =   20.D0/  1309.D0
      CCA(8) =  162.D0/ 17017.D0
      CCA(9) =   20.D0/  2431.D0
      CCA(10) = 4410.D0/371943.D0
*
*
*  ***** AVERAGE INTERACTIONS FOR NON-EQUIVALENT ELECTRONS
*
*  *****  S - ( S, P, D, F, G )
*
      CCB(1) = 1.D0/ 2.D0
      CCB(2) = 1.D0/ 6.D0
      CCB(3) = 1.D0/10.D0
      CCB(4) = 1.D0/14.D0
      CCB(5) = 1.D0/18.D0
*
*  *****  P - ( P, D, F, G )
*
      CCB(6) = 1.D0/  6.D0
      CCB(7) = 1.D0/ 15.D0
      CCB(8) = 1.D0/ 15.D0
      CCB(9) = 3.D0/ 70.D0
      CCB(10) = 3.D0/ 70.D0
      CCB(11) = 2.D0/ 63.D0
      CCB(12) = 2.D0/ 63.D0
      CCB(13) = 5.D0/198.D0
*
*  *****  D - ( D, F, G )
*
      CCB(14) =  1.D0/ 10.D0
      CCB(15) =  1.D0/ 35.D0
      CCB(16) =  1.D0/ 35.D0
      CCB(17) =  3.D0/ 70.D0
      CCB(18) =  2.D0/105.D0
      CCB(19) =  5.D0/231.D0
      CCB(20) =  1.D0/ 35.D0
      CCB(21) = 10.D0/693.D0
      CCB(22) =  5.D0/286.D0
*
*  *****  F - ( F, G )
*
      CCB(23) =  1.D0/  14.D0
      CCB(24) =  2.D0/ 105.D0
      CCB(25) =  1.D0/  77.D0
      CCB(26) = 50.D0/3003.D0
      CCB(27) =  2.D0/  63.D0
      CCB(28) =  1.D0/  77.D0
      CCB(29) = 10.D0/1001.D0
      CCB(30) = 35.D0/2574.D0
*
*  *****  G - ( G )
*
      CCB(31) =   1.D0/   18.D0
      CCB(32) =  10.D0/  693.D0
      CCB(33) =   9.D0/ 1001.D0
      CCB(34) =  10.D0/ 1287.D0
      CCB(35) = 245.D0/21879.D0
      RETURN
      END
*
*     ------------------------------------------------------------------
*                       L O O K - T M
*     ------------------------------------------------------------------
*
*     Add the deviations to the average energy for a partially filled
*       p- or d- shell
*
      SUBROUTINE LOOKTM( L, SL, SEN, Q, IP, NSL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      CHARACTER SL*2, SEN*1
      INTEGER IPTR(5)
      CHARACTER*3 TERMS(51)
      DATA    IPTR/6,11,19,35,51/
      DATA        TERMS/'3P2','1D2','1S0','4S3','2D3','2P1',
*             .. d2 and d3 terms
     :      '3F2','3P2','1G2','1D2','1S0','4F3','4P3','2H3','2G3',
     :      '2F3','2D1','2D3','2P3',
*            ... d4 terms ...
     :       '5D4','3H4','3G4','3F2','3F4','3D4','3P2','3P4',
     :       '1I4','1G2','1G4','1F4','1D2','1D4','1S0','1S4',
*            ... d5 terms ...
     :       '6S5','4G5','4F3','4D5','4P3','2I5','2H3','2G3',
     :       '2G5','2F3','2F5','2D1','2D3','2D5','2P3','2S5'/

*
*  --- search for a partially unfilled p- or d-shell
*
      N = int(Q)
      IF (N .GT. 2*L+1) N = 4*L+2 - N
      IP = 0
      NSL = 0
      IF (N .GT. 1  .AND. L .LE. 2) THEN
         IF (L .EQ. 1) THEN
            IBEGIN = 1
            IEND = 6
         ELSE
            IBEGIN = IPTR(N-1) + 1
            IEND = IPTR(N)
         END IF
1        I = IBEGIN
10       IF (SL .EQ. TERMS(I)(1:2)) THEN
            IF (SEN .EQ. ' ' .OR. SEN .EQ. TERMS(I)(3:3)) THEN
               NSL = NSL + 1
               IP = I
            END IF
         END IF
         I = I+1
         IF (I .LE. IEND) GO TO 10
      ELSE IF ( N .EQ. 1 .AND. SL(1:1) .EQ. '2') THEN
         NSL = 1
CGG
      ELSE IF ( L .EQ. 3) THEN
         CALL LOOKF( L, SL, SEN, N, IP, NSL)
CGG
      END IF
      RETURN
      END
*
*     -----------------------------------------------------------------
*                L O O K - U P
*     -----------------------------------------------------------------
*
      SUBROUTINE LOOKUP(TAB,P1,P2,IND,NO,KEY)
      INTEGER TAB(*),P1,P2,IND,NO,KEY
      DO 40 I = P1,P2
         IF (TAB(I).EQ.KEY) THEN
            NO = NO + 1
            IND = I
         END IF
40    CONTINUE
      END
*
*     -----------------------------------------------------------------
*                 L V A L
*     -----------------------------------------------------------------
*
*
      INTEGER FUNCTION LVAL(SYMBOL)
      CHARACTER*1 SYMBOL
      CHARACTER*22 SET
      DATA         SET/'spdfghiklmnSPDFGHIKLMN'/
      LOCATE = INDEX(SET,SYMBOL)
      IF ( LOCATE .LE. 11) THEN
            LVAL = LOCATE - 1
         ELSE
            LVAL = LOCATE - 12
      ENDIF
      RETURN
      END
*
*     ------------------------------------------------------------------
*               M E N U
*     ------------------------------------------------------------------
*
*
*     This routine evaluates a variety of atomic parameters as
*     requested by the user.
*
*
      SUBROUTINE MENU

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      CHARACTER CONFIG*50,EL*3,ATOM*6,TERM*6
      COMMON /LABEL/CONFIG,EL(NWFD),ATOM,TERM
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      LOGICAL FAIL,OMIT,REL,ALL,TRACE
      COMMON /TEST/FAIL,OMIT,REL,ALL,TRACE
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :  S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR

      CHARACTER*3 EL1,EL2,EL3,EL4,FUNC*1

4     WRITE(0,5)
5     FORMAT(//5X,
     :  ' These various functions are available:',//
     :  10X,'1 - EXPECTATION VALUES OF R**K'/
     :  10X,'2 - SLATER OR MAGNETIC INTEGRALS'/
     :  10X,'3 - ELECTRON DENSITY AT THE NUCLEUS'/
     :  10X,'4 - SPIN-ORBIT PARAMETER'/
     :  10X,'5 - TRANSITION INTEGRALS'/
     :  10X,'6 - EXIT TO MAIN PROGRAM'/)
      WRITE(0,'(5X,A)')'Input number corresponding to your selection:'
      READ(5,'(I1)') IFUNC
      GO TO (10,20,30,40,50,60) IFUNC

*  ****  COMPUTE EXPECTATION VALUES

10    WRITE(0,'(/5X,A,/A,T22,A)')
     :  'INPUT LABEL FOR ELECTRON FOLLOWED BY k: Example',
     :  '   2p  3','FORMAT(1X,A3,I3) '
      READ(5,'(1X,A3,I3)') EL1,K
      CALL EPTR(EL,EL1,I,*10)
      RKEV = QUADR(I,I,K)
      WRITE(3,12) EL1,K,EL1,RKEV
      WRITE(0,12) EL1,K,EL1,RKEV
12    FORMAT(/15X,' VALUE OF <',A3,'|R**',I2,'|',A3,'> = ',1PD14.7,
     :         ' a.u.'/)
      GO TO 4

*  ****  DETERMINE SLATER INTEGRALS  FK, GK, RK, NK, MK, VK


20    WRITE(0,'(/5X,A/A,T22,A)')
     : 'INPUT PARAMETERS FOR  Fk,Gk,Rk,Nk,Mk or Vk INTEGRAL: Example',
     : ' F 0( 1s, 2s)','FORMAT:  (A1,I2,1X,4(A3,1X)) '
      READ(5,'(A1,I2,1X,4(A3,1X))') FUNC,K,EL1,EL2,EL3,EL4
      IF (FUNC .GE. 'a' .AND. FUNC .LE. 'z')
     :   FUNC = CHAR(ICHAR(FUNC) + ICHAR('A') - ICHAR('a'))
      CALL EPTR(EL,EL1,I1,*20)
      CALL EPTR(EL,EL2,I2,*20)
      IF ( EL3 .NE. ' ') CALL EPTR(EL,EL3,I3,*20)
      IF ( EL4 .NE. ' ') CALL EPTR(EL,EL4,I4,*20)
      IF (FUNC .EQ. 'F' ) THEN
        SI = FK(I1,I2,K,REL)
      ELSE IF (FUNC .EQ. 'G') THEN
        SI = GK(I1,I2,K,REL)
      ELSE IF (FUNC .EQ. 'R') THEN
        SI = RK(I1,I2,I3,I4,K,REL)
      ELSE IF (FUNC .EQ. 'N') THEN
        SI = SN(I1,I2,I2,I1,K)
      ELSE IF (FUNC .EQ. 'M') THEN
        SI = SN(I1,I2,I1,I2,K)
      ELSE IF (FUNC .EQ. 'V') then
        SI = VK(I1,I2,I2,I1,k) - VK(I2,I1,I1,I2,K)
      ELSE
        WRITE(0,41)
41      FORMAT(15X,'INTEGRAL UNKNOWN: RE-ENTER')
        GO TO 20
      END IF
      IF (FUNC .NE. 'R') THEN
         WRITE(3,25) FUNC,K,EL1,EL2,SI,219474.D0*SI
         WRITE(0,25) FUNC,K,EL1,EL2,SI,219474.D0*SI
25       FORMAT(/15X,
     :  'INTEGRAL  ',A1,I2,'(',A3,',',A3,') = ',1PD14.7,' a.u.'/
     :   40X,0PF14.3,' cm-1'/)
      ELSE
         WRITE(3,26) FUNC,K,EL1,EL2,EL3,EL4,SI,219474.D0*SI
         WRITE(0,26) FUNC,K,EL1,EL2,EL3,EL4,SI,219474.D0*SI
26       FORMAT(/15X,
     :  'INTEGRAL  ',A1,I2,'(',2A3,',',2A3,') = ',1PD14.7,' a.u.'/
     :   46X,0PF14.3,' cm-1'/)
      END IF
      GO TO 4

*  ****  COMPUTE ELECTRON DENSITY AT THE NUCLEUS

30    WRITE(0,'(/5X,A/A,T22,A)')
     :  'INPUT IDENTIFYING LABEL FOR ELECTRON: Example',
     :  '   1s','FORMAT(1X,A3) '
      READ(5,'(1X,A3)') EL1
      CALL EPTR (EL,EL1,I,*30)
      LL = L(I)
      IF (LL .EQ. 0) THEN
         D = AZ(I)**2
      ELSE
         D = 0
      END IF
      WRITE(3,32) EL1,D
      WRITE(0,32) EL1,D
32    FORMAT(/15X,'DENSITY AT THE NUCLEUS FOR ',A3,' = ',1PD14.7,
     :     ' a.u.'/)
      GO TO 4

*  ****  COMPUTE SPIN-ORBIT PARAMETER

40    WRITE(0,'(/,5X,A/A,T22,A)')
     :  'INPUT IDENTIFYING LABEL FOR ELECTRON: Example',
     :  '   2p','FORMAT(1X,A3) '
      READ(5,'(1X,A3)') EL1
      CALL EPTR (EL,EL1,I,*40)
      ZETA = 0.D0
      IF (L(I) .NE. 0)  ZETA = BWZETA(I)
      ZETACM = 219474*ZETA
      WRITE(3,43) EL1,ZETA,ZETACM
      WRITE(0,43) EL1,ZETA,ZETACM
43    FORMAT(/15X,'SPIN-ORBIT PARAMETER FOR ',A3,' = ',1PD14.7,
     :        ' a.u.'/46X,0PF14.3,' cm-1'/)
      GO TO 4

*  ****  COMPUTE TRANSITION INTEGRALS

50    WRITE(0,'(/5X,A/A,T22,A)')
     :  'INPUT IDENTIFYING LABELS AND POWER OF R: Example',
     :  ' T 1( 2s, 2p)','FORMAT:  (A1,I2,2(1X,A3)) '
      READ(5,'(A1,I2,1X,A3,1X,A3)') FUNC,K,EL1,EL2
      CALL EPTR (EL,EL1,I1,*50)
      CALL EPTR (EL,EL2,I2,*50)
      TI = QUADR(I1,I2,K)
      WRITE(3,52) FUNC,K,EL1,EL2,TI
      WRITE(0,52) FUNC,K,EL1,EL2,TI
52    FORMAT(/15X,
     :  'INTEGRAL  ',A1,I2,'(',A3,',',A3,') = ',1PD14.7,' a.u.'/)
      GO TO 4
60    RETURN
      END
*
*     ------------------------------------------------------------------
*               M E T H O D
*     ------------------------------------------------------------------
*
*       Uses M1, M2, or M3 to solve the radial equation. If the input
*   data indicated METH(I) = 3, then this  solution  is  returned  to
*   DE.  Otherwise,  the routine searches for an acceptable  solution
*   which  is  both  positive  near  the  origin and has the required
*   number  of nodes.
*
      SUBROUTINE METHD1(I)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      CHARACTER CONFIG*50,EL*3,ATOM*6,TERM*6
      LOGICAL V2,FIRST
      LOGICAL FAIL,OMIT,REL,ALL,TRACE
      COMMON /TEST/FAIL,OMIT,REL,ALL,TRACE
      COMMON /COEFF/COEF(100),IJPTR(5,5)
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :   S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /LABEL/CONFIG,EL(NWFD),ATOM,TERM
      COMMON P2(NOD),HQ(NOD),XX(NOD),V,B4,CN,C,XY,XP,
     :     AZZ,PP,FN,EM,FM,EU,FU,DELTAE,M,NODE,MK,KK,NJ
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
17    CALL SOLVE(I,FIRST,REL)
*
*  *****  IF KK EQUALS 3, OMIT THE NODE CHECKING
*
      IF (KK .EQ. 3) GO TO 51
*
*  *****  COUNT THE NUMBER OF NODES
*
      MN = M
      NC = NODEC(MN)
      IF (TRACE) WRITE (6,99) EL(I),NC,MN,NJ,PDE(MN),ED,EU,EM,DELTAE
99    FORMAT(2X,A3,' NC =',I3,' MN =',I3,' NJ =',I3,' PDE(MN) =',
     :   D10.2,' ED =',D10.2,' EU =',D10.2,' EM =',D10.2,
     :   ' DELTAE =',D10.2)
*
*  *****  IF NODE COUNT IS OFF BY NO MORE THAN 1 AND DELTAE IS STILL
*  *****  QUITE LARGE, APPLY THE DELTAE CORRECTION
*
      IF (IABS(NC-NODE) .EQ. 1 .AND. DABS(DELTAE/ED) .GT. 0.02D0)
     :      GO TO 46
*
*  *****  BRANCH ACCORDING TO WHETHER THE NODE COUNT IS TOO SMALL,
*  *****  JUST RIGHT, OR TOO LARGE
*
12    IF (NC - NODE ) 8,9,10
*
*  *****  THE SOLUTION HAS THE CORRECT NUMBER OF NODES
*
9     V2 = DABS(DELTAE) .LT. 1.D-3 .OR. DABS(DELTAE)/ED .LT. 1.D-5
      IF (PDE(MN) .LT. D0 .AND. .NOT. V2) GO TO 46
      IF (PDE(MN) .GT. D0) GO TO 51
      DO 52 J = 1,NO
52    PDE(J) = - PDE(J)
      PP = -D2 - PP
51    AZZ = AZD*(D1 + PP)
      E(I,I) = ED
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
      IF ( EU .LE. EM ) WRITE(6,30) EM,EU,ED
30    FORMAT(6X,48HWARNING: DIFFICULTY WITH NODE COUNTING PROCEDURE/
     :   6X,42HLOWER BOUND ON ED GREATER THAN UPPER BOUND/
     :   6X,5HEL = ,F10.6,7H  EU = ,F10.6,7H  ED = ,F10.6)
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
27    WRITE (6,28) KK,EL(I),NC,NJ,ED,EM,EU
28    FORMAT(10X,6HMETHOD,I2,38H UNABLE TO SOLVE EQUATION FOR ELECTRON,
     :   A3/10X,5HNC = ,I3,3X,5HNJ = ,I3,3X,5HED = ,F10.6,3X,5HEL = ,
     :   F10.6,3X,5HEU = ,F10.6)
      FAIL = .TRUE.
      RETURN
      END

*
*     ------------------------------------------------------------------
*               N M R V S
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
      SUBROUTINE NMRVS(NJ,DELTA,MM,PDE,F)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      DIMENSION PDE(NOD),F(NOD),A(150),D(150)
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      EQUIVALENCE (G,G3)
*
*  *****  INTEGRATE OUTWARD TO NJ+1
*
      Y1 = PDE(1)
      Y2= PDE(2)
      G1 = YR(1)
      G2 = YR(2)
      M = NJ + 1
      DO 1 I = 3,M
      G3 = YR(I)
      Y3 = (Y2+Y2-Y1 + (D10*G2*Y2 + G1*Y1) + F(I-1)) / (D1 - G3)
      PDE(I) = Y3
      Y1 = Y2
      Y2 = Y3
      G1 = G2
1     G2 = G3
      DELTA = Y3
*
*  *****  APPLY THE TAIL PROCEDURE
*
      K = 1
      PDE(M) = -(D1 - G1)*Y1 + F(M)
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
      PDE(M) = -PDE(M-1)*RATIO + F(M)
      IF (DABS(PDE(M))+DABS(PDE(M-1)) .GT. TOL .OR. K .LT. 9) GO TO 22
20    CON =DSQRT(EH)*DEXP(-DSQRT(DABS(G/CH-.25)/RR(M))*(R(M+1)-R(M)))
      PDE(M) = PDE(M)/(D(K) + CON*(D1-  YR(M+1)))
      J = M+1
      DO 2 I= J,NO
2     PDE(I) = D0
      DO 3 J = 2,K
      I = M-J+1
      II = K-J+1
3     PDE(I) = (PDE(I)-A(II+1)*PDE(I+1))/D(II)
*
*  *****  SET DELTA = DIFFERENCE OF THE TWO SOLUTIONS AT NJ+1
*  *****         MM = NUMBER OF POINTS IN THE RANGE OF THE SOLUTION
*
      DELTA = DELTA - PDE(I)
      MM = M
      RETURN
23    WRITE (6,24)
24    FORMAT(6X,52HWARNING: FUNCTIONS TRUNCATED BY NMRVS IN TAIL REGION)
      GO TO 20
      END
*
*     ------------------------------------------------------------------
*               N O D E C
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
      PARAMETER(NOD=220,NWFD=20)
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :   S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
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
*
*     ------------------------------------------------------------------
*               O R T H O G
*     ------------------------------------------------------------------
*
*       This routine orthogonalizes the set of radial functions when an
*   orthogonality constraint applies.  A Gram-Schmidt type of  process
*   is used.
*
      SUBROUTINE ORTHOG
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      CHARACTER CONFIG*50,EL*3,ATOM*6,TERM*6
      LOGICAL FAIL,OMIT,REL,ALL,TRACE,CHANGE
      COMMON /TEST/FAIL,OMIT,REL,ALL,TRACE
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      COMMON P2(NOD),HQ(NOD),XX(NOD),V,B4,CN,C,XY,XP,
     :     AZZ,PP,FN,EM,FM,EU,FU,DELTAE,M,NODE,MK,KK,NJ
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :   S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
      COMMON /LABEL/CONFIG,EL(NWFD),ATOM,TERM
      COMMON /COEFF/COEF(100),IJPTR(5,5)
*
      IF (NWF .EQ. 1 .OR. IB .GT. NWF) RETURN
      WRITE (6,26)
26    FORMAT(/)
      II = MAX0(2,IB)
      DO 2 I = II,NWF
         CHANGE = .FALSE.
         AZZ = AZ(I)
         DO 60 J = 1,I-1
         IF (E(I,J) .NE. D0 ) THEN
*
*        ORTHOGONALITY CONDITION APPLIES
*
            C = QUADR(I,J,0)
            IF (DABS(C) .GT. 1.D-10) THEN
                WRITE(6,63) EL(J),EL(I),C
63              FORMAT(6X,'<',A3,'|',A3,'>=',1PD8.1)
                M = MAX0(M,MAX(J))
                DO 64 JJ = 1,M
                   P(JJ,I) = P(JJ,I) - C*P(JJ,J)
64              CONTINUE
                AZZ = AZZ - C*AZ(J)
                CHANGE = .TRUE.
            END IF
         END IF
60       CONTINUE
         IF (CHANGE) THEN
            PNN = DSQRT(QUADR(I,I,0))
         IF (P(1,I) .LT. D0) PNN = - PNN
         DO 66 JJ = 1,M
            P(JJ,I) = P(JJ,I)/PNN
66       CONTINUE
         AZZ = AZZ/PNN
         M = NO
67       IF (DABS(P(M,I)) .LT. 1.D-15) THEN
            P(M,I) = D0
            M = M-1
            GO TO 67
         END IF
         MAX(I) = M
         AZ(I) = AZZ
      END IF
2     CONTINUE
      RETURN
      END
*
*     ------------------------------------------------------------------
*               O U T P U T
*     ------------------------------------------------------------------
*
*       The radial functions and orthogonality integrals are printed,
*   if PRINT is .TRUE.   The  functions  will  also  be  punched  (or
*   stored) on unit OUF, if OUF .NE. 0.
*
      SUBROUTINE OUTPUT(PRINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      CHARACTER CONFIG*50,EL*3,ATOM*6,TERM*6
      LOGICAL PRINT
      INTEGER OUF
      COMMON /INOUT/ IUF,OUF
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :   S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      COMMON /LABEL/CONFIG,EL(NWFD),ATOM,TERM
      DIMENSION OUT(8)
      IF ( .NOT. PRINT ) GO TO 30
*
*  *****  PRINT RADIAL FUNCTIONS, 7 PER PAGE
*
      ML = IB
2     MU = MIN0(ML+7,NWF)
      I = MU - ML + 1
      MX = 0
      DO 1 J = ML,MU
1     MX = MAX0(MX,MAX(J))
      WRITE (3,5) ATOM,TERM,(EL(J),J=ML,MU)
5     FORMAT(1H1,9X,19HWAVE FUNCTIONS FOR  ,2A6//10X,1HR,8(10X,A3))
      K= 0
      KK = 0
      DO 6 J = 1,MX
      DO 9 JJ = ML,MU
      IJ = JJ - ML + 1
9     OUT(IJ) = P(J,JJ)*R2(J)
      K = K+1
      IF (K .LE. 10) GO TO 6
      K = 1
      KK = KK+1
      IF (KK .LT. 5) GO TO 21
      KK = 0
      WRITE (3,23)
23    FORMAT(1H1//)
      GO TO 6
21    WRITE (3,8)
8     FORMAT(1X)
6     WRITE (3,10) R(J),(OUT(JJ),JJ=1,I)
10    FORMAT(F13.5,F15.6,7F13.6)
      DO 15 J = ML,MU
      IJ = J - ML + 1
15    OUT(IJ) = DPM(J)
      WRITE (3,16) (OUT(J),J=1,I)
16    FORMAT(4X,10HMAX. DIFF. ,F15.7,7F13.7)
      ML = ML+8
      IF (ML .LE. NWF) GO TO 2
30    IF ( OUF .EQ. 0) GO TO 14
*
*  *****  OUTPUT FUNCTIONS ON UNIT OUF FOR FUTURE INPUT
*
      DO 3 I = 1,NWF
      MX = MAX(I)
      WRITE (OUF) ATOM,TERM,EL(I),MX,Z,E(I,I),EK(I),AZ(I),
     :   (P(J,I),J=1,MX)
3     CONTINUE
*
14    RETURN
      END
*
*     ------------------------------------------------------------------
*               P O T L
*     ------------------------------------------------------------------
*
*       Computes and stores the potential function
*                              2(k-1)
*              YR = SUM  a    Y      (j,j;r)
*                   j,k   ijk
*
      SUBROUTINE POTL(I,REL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      LOGICAL REL
      DO 1 J=1,NO
1     YR(J) = D0
      DO 2 J=1,NWF
      DO 5 K = 0,2*MIN0(L(I),L(J)),2
      C = A(I,J,K)
      IF (DABS(C) .LE. 1.D-8) GO TO 5
      CALL YKF(J,J,K,REL)
      DO 7 JJ=1,NO
7     YR(JJ) = YR(JJ) + C*YK(JJ)
5     CONTINUE
2     CONTINUE
      RETURN
      END
*
*     ------------------------------------------------------------------
*               Q U A D
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      DIMENSION F(NOD),G(NOD)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
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
*                 Q U A D R
*     ------------------------------------------------------------------
*
*                                   kk
*       Evaluates the integral of  r   P (r) P (r) with respect to r
*                                       i     j
*
      DOUBLE PRECISION FUNCTION QUADR(I,J,KK)
      PARAMETER (NOD=220,NWFD=20)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      K = KK + 2
      LI = L(I)
      LJ = L(J)
      DEN = LI + LJ + 1 + K
      ZR = Z*R(4)
      BI = (P(4,I)/(AZ(I)*R2(4)*R(4)**LI) - D1+ZR/(LI+1) )/ZR**2
      BJ = (P(4,J)/(AZ(J)*R2(4)*R(4)**LJ) - D1+ZR/(LJ+1) )/ZR**2
      ALPHA= (D1/(LI + 1) + D1/(LJ + 1))/(DEN + D1)
      ZR = Z*R(1)
      BETA = (DEN+D1)*ALPHA**2 - D2*(BI+BJ+D1/((LI+1)*(LJ+1)))/(DEN+D2)
      D = P(1,I)*P(1,J)*R(1)**K*(((BETA*ZR+ALPHA)*ZR+D1)/(DEN*H1)+D5)
      DD = D0
      M = MIN0(MAX(I),MAX(J)) - 1
      DO 1 JJ = 2,M,2
      JP = JJ + 1
      D = D +P(JP,I)*P(JP,J)*R(JP)**K
      DD = DD + P(JJ,I)*P(JJ,J)*R(JJ)**K
1     CONTINUE
      QUADR = H1*(D + D2*DD)
      RETURN
      END
*
*     ------------------------------------------------------------------
*                 Q U A D S
*     ------------------------------------------------------------------
*
*                                       kk
*       Evaluates the integral of  (1/r)   YK(r) P (r) P (r)  with
*                                                 i     j
*   respect to r.
*
      DOUBLE PRECISION FUNCTION  QUADS(I,J,KK)
      PARAMETER (NOD=220,NWFD=20)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      DEN = L(I) + L(J) + 3
      K = 2 - KK
      CD = D1 + Z*R(1)*(DEN-D1)/((DEN+D1)*((L(I)+1)*(L(J)+1)))
      D = YK(1)*P(1,I)*P(1,J)*R(1)**K*( CD/(DEN*H1)+ D5)
      DD = D0
      MX = MIN0(MAX(I),MAX(J)) - 1
      DO 1 M = 2,MX,2
      DD = DD + YK(M)*P(M,I)*P(M,J)*R(M)**K
      D= D+  YK(M+1)*P(M+1,I)*P(M+1,J)*R(M+1)**K
1     CONTINUE
      QUADS = H1*(D + D2*DD)
      RETURN
      END
*
*     ------------------------------------------------------------------
*               R E F O R M
*     ------------------------------------------------------------------
*
*     Convert the free-format STR1 to the fixed 5(1X,A3,1X,I4,1X) format
*     for STR2
*
      SUBROUTINE REFORM(STR1,STR2)
      CHARACTER*50 STR1,STR2,BLANK
      DATA                   BLANK/'   '/
*
    1 I = 0
      STR2 = BLANK
      IS = 0
    2 JS = INDEX(STR1(IS+1:),'(')
      IF (JS .NE. 0) THEN
         IF (JS .GT. 5) GO TO 10
         I = I+5
         STR2(I-JS+1:I) = STR1(IS+1:IS+JS)
         IS = IS + JS
         JS = INDEX(STR1(IS+1:),')')
         IF (JS .EQ. 0 .OR. JS .GT. 5) GO TO 10
         I = I+5
         STR2(I-JS+1:I) = STR1(IS+1:IS+JS)
         IS = IS + JS
         GO TO 2
      END IF
      RETURN
   10 WRITE(0,*)' Error in ',STR1,': Re-enter'
      READ(5,'(A)') STR1
      GO TO 1
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
        CHARACTER*3 OF(20), ELC
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
*
*
*     ------------------------------------------------------------------
*                 R K
*     ------------------------------------------------------------------
*
*                   k
*       Evaluates  R (i, j; ii, jj)
*
      DOUBLE PRECISION FUNCTION RK(I,J,II,JJ,K,REL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL REL
      CALL YKF(I,II,K,REL)
      RK = QUADS(J,JJ,1)
      RETURN
      END
*
*     ------------------------------------------------------------------
*                R M E
*     ------------------------------------------------------------------
*
*
      DOUBLE PRECISION FUNCTION RME(L,LP,K)
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      COMMON/FACT/GAM(100)
*
*--- EVALUATES THE REDUCED MATRIX ELEMENT (L//C(K)//LP)  -  SEE FANO
*    AND RACAH, IRREDUCIBLE TENSORIAL SETS, CHAP. 14, P. 81
*
      IF (MIN0(L,LP) .EQ. 0) THEN
         RME = 1.D0
       ELSE IF ( K .EQ. 0) THEN
         RME = 2*L+1
         RME = DSQRT(RME)
       ELSE
         I2G=L+LP+K
         IG=I2G/2
         IF (I2G - 2*IG .NE. 0) THEN
             RME = 0.D0
         ELSE
            I1=IG-L
            I2=IG-LP
            I3=IG-K
            QUSQRT=(2*L+1)*(2*LP+1)
            RME=DSQRT(QUSQRT)*DEXP((GAM(2*I1+1)+GAM(2*I2+1)+GAM(2*I3+1)-
     :        GAM(I2G+2))/2.D0 +GAM(IG+1)-GAM(I1+1)-GAM(I2+1)-GAM(I3+1))
         END IF
      END IF
      RETURN
      END
*
*     ------------------------------------------------------------------
*               R O T A T E
*     ------------------------------------------------------------------
*
*        This routine analyses the energy expression to determine the
*   stationary condition with respect to rotation of orbials i and j.
*   If the condition is zero, the off-diagonal energy parameters may
*   be set to zero;  otherwise the orbials are rotated so as to satisfy
*   the stationay condition to first order in the perturbation.
*
*
      SUBROUTINE ROTATE(I,J)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      CHARACTER CONFIG*50,EL*3,ATOM*6,TERM*6
      LOGICAL FAIL,OMIT,REL,ALL,TRACE
      COMMON /TEST/FAIL,OMIT,REL,ALL,TRACE
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      COMMON QI,QJ,G,DI,DJ,DII,DJJ,DIJ,DJI
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :   S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
      COMMON /LABEL/CONFIG,EL(NWFD),ATOM,TERM
      COMMON /COEFF/COEF(100),IJPTR(5,5)
      ALL = .TRUE.
      G = D0
      DG = D0
      QI = SUM(I)
      QJ = SUM(J)
      IF (QI .EQ. D2*(2*L(I)+1) .AND. QJ .EQ. D2*(2*L(J)+1))
     :   GO TO 44
      IF (DABS(QI - QJ) .LT. 1.D-14) GO TO 16
      C = D5*(QI - QJ)
      G = G -C*HL(EL,I,J,REL)
      DG = DG -C*(HL(EL,I,I,REL) - HL(EL,J,J,REL))
*
16    DO 13 K = 0,2*L(I),2
      C = QI*(A(I,I,K) - A(I,J,K) - B(I,J,K))
      IF (DABS(C) .LT. 1.D-8) GO TO 21
      G = G  + C*RK(I,I,I,J,K,REL)
      FKII = FK(I,I,K,REL)
      FKIJ = FK(I,J,K,REL)
      GKIJ = GK(I,J,K,REL)
      DG = DG +C*(FKII - FKIJ - D2*GKIJ)
21    CJ = QJ*(A(J,J,K) - A(J,I,K) - B(J,I,K))
      IF (DABS(CJ) .LT. 1.D-8) GO TO 13
      FKJJ = FK(J,J,K,REL)
      IF (DABS(C) .GE. 1.D-8) GO TO 22
      FKIJ = FK(I,J,K,REL)
      GKIJ = GK(I,J,K,REL)
22    G = G - CJ*RK(J,J,J,I,K,REL)
      DG = DG + CJ*(FKJJ -FKIJ - D2*GKIJ)
13    CONTINUE
      DO 14 M = 1,NWF
      IF (M .EQ. I .OR. M.EQ. J) GO TO 14
      DO 15 K = 0,2*MIN0(L(I),L(M)),2
      C = A(I,M,K)*QI - A(J,M,K)*QJ
      IF (DABS(C) .LT. 1.D-8) GO TO 23
      G = G + C*RK(I,M,J,M,K,REL)
      DG = DG + C*(FK(I,M,K,REL) - FK(J,M,K,REL))
23    KK = IABS(L(I)-L(M)) + K
      C = B(I,M,KK)*QI - B(J,M,KK)*QJ
      IF  (DABS(C) .LT. 1.D-8) GO TO 15
      G = G + C*RK(I,J,M,M,KK,REL)
      DG = DG + C*(GK(I,M,KK,REL) - GK(J,M,KK,REL))
15    CONTINUE
14    CONTINUE
      IF (DABS(QI-QJ) + DABS(G) + DABS(DG)
     :   .LE. 2.D-8 ) GO TO 44
101   IF (DABS(G)+DABS(DG) .GT. 1.D-8 .OR. DABS(E(I,J)) .GT. 2.D-5) THEN
         EPS = G/DG
         EPS = DSIGN(DMIN1(DABS(EPS),0.2D0),EPS)
         DD = DSQRT(D1 + EPS*EPS)
         DO 41 JJ = 1,NO
            PI = (P(JJ,I) + EPS*P(JJ,J))/DD
            P(JJ,J) = (P(JJ,J) - EPS*P(JJ,I))/DD
41          P(JJ,I) = PI
      ELSE
         EPS = D0
      END IF
      WRITE (6,100) EL(I),EL(J),G,EL(I),EL(J),DG,EPS
100   FORMAT(10X,'C(',2A3,') =',F12.5,3X,'V(',2A3,') =',F12.5,
     :       3X,'EPS =',F9.6)
      RETURN
*
*  *****  THE ENERGY IS STATIONARY WITH RESPECT TO ROTATIONS
*
44    E(I,J) = 1.D-10
      E(J,I) = 1.D-10
      RETURN
      END
*
*     ------------------------------------------------------------------
*                 S H I F T
*     ------------------------------------------------------------------
*
*       Computes the mass velocity  and one-body   Darwin term
*   corrections for the relativistic shift in the energy of the electron
*   including non-diagonal corrections
*
*
      DOUBLE PRECISION FUNCTION RLSHFT(I1,I2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
*
*  *****  FORM DD  -L(L+1)/RR|P(I)>
*
      FL = L(I1)
      C = (FL+D5)**2
      LL = L(I1) + 1
      L2 = 2*L(I1) + 1
      L3 = 2*L(I1) + 3
      ZZ = Z*Z
      HH = 180.D0*H*H
      MX = MAX0(MAX(I1), MAX(I2))
      DO 1 J = 2,MX
1     YK(J) = -D1/RR(J)
*
*  *****  FORM THE INTEGRAND
*
      I = I1
      A1 = D0
      B1 = D0
      DO 2 KK = 1,2
      B2 = B1
      YY = (P(3,I)+P(1,I) - D2*P(2,I))/(H*H) - C*P(2,I)
      YK(2) = YY*YK(2)
      YK(3) = YK(3)*((-(P(5,I)+P(1,I)) + D16*(P(4,I)+P(2,I))
     :        -D30*P(3,I))/(D12*H*H) - C*P(3,I))
      MM = MAX(I) - 3
      DO 11 K =  4,MM
      YY = D2*(P(K+3,I)+P(K-3,I))
      YY = YY - 27.D0*(P(K+2,I)+P(K-2,I))
      YY = YY +  270.D0*(P(K+1,I)+P(K-1,I)) - 490.D0*P(K,I)
      YY = YY/HH - C*P(K,I)
      YK(K) = YY*YK(K)
      IF (K .EQ. 4) B1 = (YY/(D2*Z*P(4,I)*R(4)) + D1)/R(4)
11    CONTINUE
      MM = MM + 1
      YK(MM) = YK(MM)*((-(P(MM+2,I)+P(MM-2,I))+D16*(P(MM+1,I)+P(MM-1,I))
     :      -D30*P(MM,I))/(D12*H*H) - C*P(MM,I))
      MM = MM + 1
      YK(MM) = YK(MM)*((P(MM+1,I) + P(MM-1,I) - D2*P(MM,I))/(H*H)
     :         - C*P(MM,I))
      A2 = A1
      A1 = (P(1,I)/(AZ(I)*R(1)**L(I)*R2(1)) - D1 +Z*R(1)/LL)/RR(1)
      I = I2
2     CONTINUE
*
*  ***** DETERMINE CONTRIBUTION FROM NEAR THE NUCLEUS
*
      A = (Z/LL - L2*(B1 + B2)/D2)/LL
      B = (L2*B1*B2 - D2*(A1 + A2) + (Z/LL**2)*(D2*Z*(D1 + D1/LL)
     :     - L2*(B1 + B2)))/L3
      RELSH = -P(4,I1)*P(4,I2)*(D1 + A*R(4) + B*RR(4))*D4*ZZ/L2
      RELSH = RELSH/H1 - D5*YK(4)
      RELSH2 = D0
*
*  *****  INTEGRATE
*
      DO 50 J = 5,MX,2
      RELSH2 = RELSH2 + YK(J)
50    RELSH = RELSH  + YK(J-1)
      RELSH  = (RELSH + D2*RELSH2)*H1
      IF ( L(I1) .EQ. 0 ) RELSH = RELSH + Z*AZ(I1)*AZ(I2)
      RLSHFT = RELSH*D5*FINE
      RETURN
      END
*
*     ------------------------------------------------------------------
*               S C A L E
*     ------------------------------------------------------------------
*
*       The current radial functions are scaled according to the
*   procedures of Sec. 7-2  .   Values of AZ and E(I,I), the starting
*   values and the diagonal energy parameters are also scaled.
*
*
      SUBROUTINE SCALE(ZZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :   S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
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
      E(I,I) = E(I,I)*SC**2
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
     :   THETA*(THETA - D1)*(F0 - F1 - F2 + F3)/D4
      GO TO 4
7     P(J,I) = D0
4     CONTINUE
      MAX(I) = NO
*
*  ***** NORMALIZE THE INTERPOLATED FUNCTION
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
*               S C F
*     -----------------------------------------------------------------
*
*       This routine controls the SCF procedure described in Chapter
*   7.  If certain input parameters are zero (or blank) they will  be
*   set to their default value.
*
*          Parameter       Default Value
*          --------        -------------
*          SCFTOL          1.D-7
*          I*              (NWF + 1 - IB)/4 + 3
*          NSCF            12
*
*   The self-consistency convergence criterion is
*
*          Z2 = SQRT( SCFTOL*(Z*NWF/2) )
*
*   It is increased by a factor two at the end of each iteration.
*
*
      SUBROUTINE SCF(ETOTAL,SCFTOL,EREL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      CHARACTER CONFIG*50,EL*3,ATOM*6,TERM*6,ANS*1
      LOGICAL FAIL,OMIT,REL,ALL,TRACE,LAST
      COMMON /TEST/FAIL,OMIT,REL,ALL,TRACE
      COMMON /LABEL/CONFIG,EL(NWFD),ATOM,TERM
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :   S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
      INTEGER OUF
      COMMON /INOUT/ IUF,OUF
*
*  *****  SET THE SCF CONVERGENCE PARAMETER TO AN OPTIMISTIC VALUE
*
      REL = .FALSE.
      TOL = DSQRT(Z)*1.D-10
      Z2 = SCFTOL*DSQRT(Z*NWF)
      WRITE (6,15)
15    FORMAT(//)
      WRITE (6,16) OMIT,SCFTOL,NO
   16 FORMAT(10X,44HWEAK ORTHOGONALIZATION DURING THE SCF CYCLE=,L4/
     :       10X,44HSCF CONVERGENCE TOLERANCE (FUNCTIONS)      =,1PD9.2
     :      /10X,44HNUMBER OF POINTS IN THE MAXIMUM RANGE      =,I4)
*
*  *****  SET ITERATION PARAMETERS
*
      IPR = 0
      DP1 = D0
      ETOTAL = D0
      ICYCLE = 0
      IF ( IB .GT. NWF ) GO TO 17
*
*  *****  PERFORM NSCF SELF-CONSISTENT FIELD ITERATIONS
*
      LAST = .FALSE.
9     DO 100 I = 1,NSCF
      ICYCLE = ICYCLE + 1
      WRITE(6,7) ICYCLE,Z2
7     FORMAT(//10X,17HITERATION NUMBER ,I2/10X,16H----------------//
     : 10X,'SCF CONVERGENCE CRITERIA (SCFTOL*SQRT(Z*NWF)) = ',1PD9.1/)
      DP1 = D0
      CALL GRANGE
*
*  *****  SOLVE EACH DIFFERENTIAL EQUATION IN TURN
*
      WRITE(6,14)
14    FORMAT(/20X,' EL',9X,'ED',13X,'AZ',11X,'NORM',7X,'DPM')
      DO 2 JP = IB,NWF
      J = IORD(JP)
      CALL DE(J)
      IF ( FAIL ) RETURN
      DP = DPM(J)*DSQRT(SUM(J))
      IF ( DP1 .GE. DP ) GO TO 2
      DP1 = DP
      JJ = J
2     CONTINUE
      IF ( DP1 .LT. Z2) GO TO 6
*
*  *****  SOLVE IC DIFFERENTIAL EQUATIONS EACH TIME SELECTING THE
*  *****  ONE WITH THE LARGEST DPM
*
      DO 4 II =1,IC
      CALL DE(JJ)
      IF ( FAIL ) RETURN
      DP1 = D0
      DO 5 JP = IB,NWF
      J = IORD(JP)
      DP = DSQRT(SUM(J))*DPM(J)
      IF ( DP1 .GT. DP ) GO TO 5
      JJ = J
      DP1 = DP
5     CONTINUE
      IF (DP1 .LT. Z2) GO TO 6
4     CONTINUE
6     CALL ORTHOG
      IF (DP1 .LT. Z2 .AND. LAST) GO TO 17
      IF ( I .EQ. NSCF ) GO TO 1
*
*  *****  IF FUNCTIONS APPEAR TO HAVE CONVERGED,SOLVE EACH AGAIN, IN
*  *****  TURN, AND TEST AGAIN
*
      IF (DP1 .LE. Z2) LAST = .TRUE.
*
*  *****  INCREASE THE CONVERGENCE CRITERION FOR SELF-CONSISTENCY
*
1     Z2 = D2*Z2
      WRITE(3,8) EL(JJ),DP1
      WRITE(0,8) EL(JJ),DP1
8     FORMAT(/ 6X,34HLEAST SELF-CONSISTENT FUNCTION IS ,A3,
     :   27H :WEIGHTED MAXIMUM CHANGE =,1PD10.2)
100   CFGTOL = 1.4D0*CFGTOL
18    WRITE(0,13)
13    FORMAT(10X/' SCF ITERATIONS HAVE CONVERGED TO THE ABOVE ACCURACY')
      WRITE(3,13)
20    WRITE(0,'(A)')' Do you wish to continue ? (Y/N/H) '
      READ(5,'(A)') ANS
      IF (ANS .EQ. 'H' .OR. ANS .EQ. 'h') THEN
         CALL HELP(10)
         GO TO 20
         END IF
      IF (ANS .EQ. 'Y' .OR. ANS .EQ. 'y') THEN
         WRITE(0,'(A,A)')' Enter the additional iterations and new IC ',
     :                  'in FORMAT(I2, 1X, I2)'
         READ(5,'(I2, 1X, I2)') NSCF,IC
         GO TO 9
      END IF
      FAIL = .TRUE.
*
*  *****  PERFORM RELATIVISTIC AND NON-RELATIVISTIC CALCULATIONS
*
17    CONTINUE
      CALL ENERGY(ENONR)
      REL = .TRUE.
      CALL ENERGY(ETOTAL)
      REL = .FALSE.
      EREL = ETOTAL - ENONR
      NIT = NWF - IB + 1
      WRITE (3, 105) NIT, DP1
105   FORMAT(//10X,'NUMBER OF FUNCTIONS ITERATED          =',I6/
     :         10X,'MAXIMUM WEIGHTED CHANGE IN FUNCTIONS  =',D10.2/)
      RETURN
      END
*
*     ------------------------------------------------------------------
*               S E A R C H
*     ------------------------------------------------------------------
*
*       This routine searches for the NJ>70 such that YR(j) > 0 for all
*   j > NJ.
*
*
      SUBROUTINE SEARCH(NJ,I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
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
*     ------------------------------------------------------------------
*                 S N
*     ------------------------------------------------------------------
*
*                                      3              k
*       Evaluates the integral of (1/r)  P (r) P (r) Z (i, j; r)  with
*                                         i     j
*   respect to r.
*
      DOUBLE PRECISION FUNCTION SN(I,J,II,JJ,K)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      CALL ZK(J,JJ,K)
      SN = QUADS(I,II,3)*FINE
      RETURN
      END
*
*     ------------------------------------------------------------------
*               S O L V E
*     ------------------------------------------------------------------
*
*       When FIRST is .TRUE., SOLVE computes the potential and exchange
*   function and initializes variables for the i'th radial  equation.
*   The vector P1 is the solution of the radial equation and P2 the
*   variation of the solution with respect to the energy parameter
*   E(I,I).
*
*
      SUBROUTINE SOLVE(I,FIRST,REL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      CHARACTER CONFIG*50,EL*3,ATOM*6,TERM*6
      LOGICAL FIRST
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :   S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON P2(NOD),HQ(NOD),XX(NOD),V,B4,CN,C,XY,XP,
     :     AZZ,PP,FN,EM,FM,EU,FU,DELTAE,M,NODE,MK,KK,NJ
      COMMON /LABEL/CONFIG,EL(NWFD),ATOM,TERM
      DIMENSION ZERO(NOD),P1(NOD)
      EQUIVALENCE (ZERO(1),XX(1)),(PDE(1),P1(1))
      LOGICAL REL
*
*  *****  IF FIRST IS 'TRUE', CALL POTL AND XCH AND SET UP ARRAYS
*
      IF (.NOT. FIRST) GO TO 17
      CALL POTL(I,REL)
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
     :     X(3)/R(3)**RL))
*
*  *****  DETERMINE LOWER BOUND ON THE ENERGY PARAMETER
*
      IF (KK .NE. 3) GO TO 80
      DO 11 JJ = 15,ND
      J = NO - JJ
      IF (YK(J) .LT. D0 ) GO TO 63
11    CONTINUE
      WRITE (6,12)
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
     :   P(J+1,I) + D10*YK(J)*P(J,I) + YK(J-1)*P(J-1,I))-X(J)
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
      WRITE (6,24) ED
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
      DELTAE = D0
      ED = (ZINF/N(I))**2
      WRITE (6,66) EL(I), ZINF
66    FORMAT(//10X, 'RETURN HYDROGENIC FUNCTION FOR ',A3,
     :   ' WITH EFFECTIVE CHARGE ',F10.3)
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
     :   + YR(NJ+1)*Y3 + X(NJ)
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
*               S U M M R Y
*     ------------------------------------------------------------------
*
*       The results of a calculation are summarized.   These include
*   the following for each electron:
*
*          E(NL)   - diagonal energy parameter
*          I(NL)   - -(1/2)<nl|L|nl>
*          KE      - I(NL) + Z <r>
*          REL     - Relativistic shift (mass-velocity, Darwin term,
*                    spin-spin contact term)
*          SIGMA   - screening parameter as defined by Eq. (6-  ).
*          AZ(NL)  - starting parameter, P(r)/r**(l+1) as r -> 0.
*          1/R**3  - expected value of <1/r**3>
*          1/R     - expected value of <1/r>
*          R       - expected mean radius
*          R**2    - expected value of <r**2>
*
*   These results are followed by:
*
*          KINETIC ENERGY (EK)
*          POTENTIAL ENERGY (EP) = ET - EN
*          RATIO                 = EP/EN
*          NON- RELATIVISTIC ENERGY (ET - EREL)
*          RELATIVISTIC SHIFT (EREL) FOR THE STATE
*          TOTAL ENERGY (ET)
*
      SUBROUTINE SUMMRY(ET,EREL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      LOGICAL FAIL,OMIT,REL,ALL,TRACE
      COMMON /TEST/FAIL,OMIT,REL,ALL,TRACE
      CHARACTER CONFIG*50,EL*3,ATOM*6,TERM*6
      INTEGER OUF
      COMMON /INOUT/ IUF,OUF
      COMMON /LABEL/CONFIG,EL(NWFD),ATOM,TERM
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :   S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /COEFF/COEF(100),IJPTR(5,5)
      COMMON R1(NWFD),RM1(NWFD)
*
      PI = ACOS(-D1)
      WRITE (3,9) ATOM,TERM
9     FORMAT(/// 24X,'ATOM ',A6,3X,'TERM ',A6//
     :   2X,'nl',7X,'E(nl)',
     :   8X,'I(nl)',7X,'KE(nl)',8X,'Rel(nl)',3X,'S(nl)',7X,'Az(nl)')
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
15    FORMAT(1X,A3,F14.7,3F13.6,F8.3,F14.6)
10    CONTINUE
*
*  *****  Compute Moments
*
      WRITE(3,8) 'Delta(R)'
 8    FORMAT(//2X,'nl',6X,A8,5X,'1/R**3',7X,'1/R',9X,'R',8X,'R**2')
      DO 11 I = 1,NWF
      RM3 = 0
      IF (L(I) .NE. 0) RM3 = QUADR(I,I,-3)
      RP2 = QUADR(I,I,2)
      RZ = 0.
      IF ( L(I) .EQ. 0) RZ = AZ(I)**2/(4.*PI)
      WRITE(3,16) EL(I),RZ,RM3,RM1(I),R1(I),RP2
16    FORMAT(1X,A3,F14.3,F13.4,F11.5,F10.5,F11.5)
11    CONTINUE
31    ETN = ET - EREL
      EPOTL = ETN - EN
      RATIO = EPOTL/EN
      WRITE(0,26) ETN,EN,EREL,EPOTL,ET,RATIO
      WRITE(3,26) ETN,EN,EREL,EPOTL,ET,RATIO
26    FORMAT(//5X,'TOTAL ENERGY (a.u.)'/5X,'----- ------'/
     : 10X,' Non-Relativistic   ',F15.8,T50,'Kinetic   ',F15.8/
     : 10X,' Relativistic Shift ',F15.8,T50,'Potential ',F15.8/
     : 10X,' Relativistic       ',F15.8,T50,'Ratio     ',F15.9)
13    RETURN
      END
*
*     ------------------------------------------------------------------
*                 V
*     ------------------------------------------------------------------
*
*                  k
*       Evaluates V (i,j) as defined by Blume and Watson (1962).
*
      DOUBLE PRECISION FUNCTION VK(I,J,II,JJ,K)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      CALL DYK(I,II,K)
      VK = QUADS(J,JJ,2)*FINE
      RETURN
      END
*
*     ------------------------------------------------------------------
*               W A V E F N
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
*
*   The set of functions are then orthogonalized.
*
*
      SUBROUTINE WAVEFN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      CHARACTER CONFIG*50,EL*3,ATOM*6,TERM*6,
     :          EL1*3,AT*6,TT*6,ATM(NWFD)*6,
     :          TRM(NWFD)*6,TITLE*24
      LOGICAL FAIL,OMIT,REL,ALL,TRACE
      INTEGER OUF
      COMMON /INOUT/ IUF,OUF
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /TEST/FAIL,OMIT,REL,ALL,TRACE
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :   S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
      COMMON ZZ(NWFD),IND(NWFD),PN,Z2,FN,M,K,ZT,
     :   ETI,EKI,AZI,PT(NOD),MT
      COMMON /LABEL/CONFIG,EL(NWFD),ATOM,TERM
      COMMON /COEFF/COEF(100),IJPTR(5,5)
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
2     READ(IUF,END=5) AT,TT,EL1,M,ZT,ETI,EKI,AZI,(PT(J),J=1,M)
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
         E(I,I) = C*C*ETI
         EK(I)  = C*C*EKI
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
      IND(I) = 0
      WRITE(6,27) EL(I)
27    FORMAT(8X,'WAVE FUNCTIONS NOT FOUND FOR ',A3)
*
*  *****  DETERMINE ESTIMATES OF THE WAVE FUNCTIONS BY THE SCREENED
*  *****  HYDROGENIC APPROXIMATION
*
8     IF ( Z-S(I) .GT. D0 ) THEN
         PN = HNORM(N(I),L(I),Z-S(I))
      ELSE
       WRITE(0,'(A,A)') ' Effective nuclear charge zero for ',EL(I)
       STOP
      END IF
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
      WRITE (3,14)
14    FORMAT(/// 8X,18HINITIAL ESTIMATES  //10X,2HNL,
     :   4X,5HSIGMA,6X,5HE(NL),4X,6HAZ(NL),4X,9HFUNCTIONS//)
*
*  *****  COMPUTE ONE-ELECTRON ENERGY PARAMETERS IF THEY WERE NOT
*  *****  SPECIFIED ON INPUT.
*
      DO 15 I = 1,NWF
         K = IND(I) + 2
         IF ( IND(I) .EQ. -2 ) THEN
            TITLE = ' SCALED '//ATM(I)//TRM(I)
         ELSE IF (IND(I) .EQ. 0) THEN
            TITLE = ' SCREENED HYDROGENIC'
         ELSE
            TITLE = ' UNCHANGED'
         END IF
17       WRITE (3,19) EL(I),S(I),E(I,I),AZ(I),TITLE
19       FORMAT(9X,A3,F9.2,F11.3,F10.3,3X,A24)
15    CONTINUE
      RETURN
      END
*
*     ------------------------------------------------------------------
*               X C H
*     ------------------------------------------------------------------
*
*       This routine computes functions associated with the exchange
*   function for the i'th radial equation,  including   contributions
*   from   the  interactions.    The   exact   form  of  the function
*   depends on the value of IOPT.
*
*          Value of IOPT      Function
*          -------------      --------
*              1             SQRT(r) X(r)
*              2             X(r)/SQRT(r)
*              3             r SQRT(r) ( X(r) + SUM e   P )
*                                                    ij  j
*
      SUBROUTINE XCH(I,IOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NOD=220,NWFD=20)
      LOGICAL FAIL,OMIT,REL,ALL,TRACE
      COMMON /TEST/FAIL,OMIT,REL,ALL,TRACE
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /WAVE/PDE(NOD),EK(NWFD),E(NWFD,NWFD),ED,AZD,SUM(NWFD),
     :  S(NWFD),DPM(NWFD),ACC(NWFD),METH(NWFD),IORD(NWFD),IPR
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      DO 1 J=1,NO
1     X(J) = D0
      DO 2 J=1,NWF
      IF (J .NE. I) THEN
        DO 3 K=IABS(L(I)-L(J)),L(I)+L(J),2
          C = B(I,J,K)*D2
          IF (DABS(C) .LE. 1.D-10) GO TO 3
          CALL YKF(I,J,K,REL)
          DO 6 JJ =1,NO
6           X(JJ) = X(JJ)+ C*YK(JJ)*P(JJ,J)
3       CONTINUE
      END IF
2     CONTINUE
4     GO TO (75,76,77),IOPT
76    DO 78 J = 1,NO
78    X(J) = X(J)/R(J)
      GO TO 75
77    DO 5 J =1,NO
5     X(J) = R(J)*X(J)
      DO 7 J = 1,NWF
      C = E(I,J)
      IF (DABS(C) .LE. 1.D-10 .OR. (J .EQ. I)) GO TO 7
      DO 9 JJ = 1,NO
9     X(JJ) = X(JJ) + C*P(JJ,J)*RR(JJ)
7     CONTINUE
*
*  *****  CHECK IF EXCHANGE IS ZERO: IF SO, METHOD 2 SHOULD BE USED.
*
75    IF (METH(I) .EQ. 2) RETURN
      IF ( DABS(X(1)) + DABS(X(2)) + DABS(X(3)) .EQ. D0 ) METH(I) = 2
      RETURN
      END
*
*     ------------------------------------------------------------------
*                 Y K F
*     ------------------------------------------------------------------
*
*               k
*       Stores Y (i, j; r) in the array YK
*
*
      SUBROUTINE YKF(I,J,K,REL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220,NWFD=20)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      LOGICAL REL
      CALL ZK(I,J,K)
      A =    EH**(K+1)
      C = 2*K+1
      A2 = A*A
      H90 = C*H3/D30
      A3 = A2*A*H90
      AI = H90/A
      AN = 114.D0*A*H90
      A34 = 34.D0*H90
      MX = (MIN0(MAX(I),MAX(J))/2)*2
      F1 = YK(MX)*EH**K
      F2 = YK(MX)
      F3 = YK(MX-1)
      F4 = YK(MX-2)
      DO 9 M = MX-2,2,-1
      F5 = YK(M-1)
      YK(M) = YK(M+2)*A2 + ( AN*F3 + A34*(F4+A2*F2)-F5*AI-F1*A3)
      F1 = F2
      F2 = F3
      F3 = F4
9     F4 = F5
      YK(1) = YK(3)*A2+C*H3*(F4 + D4*A*F3 + A2*F2)
      IF (.NOT.REL) RETURN
      C = C*FINE
      DO 10 M = 1,MX
      YK(M) = YK(M) + C*P(M,I)*P(M,J)
10    CONTINUE
      RETURN
      END
*
*     ------------------------------------------------------------------
*                 Z K
*     ------------------------------------------------------------------
*
*               k
*       Stores Z (i, j; r) in the array YK.
*
*
      SUBROUTINE ZK(I,J,K)
      PARAMETER (NOD=220,NWFD=20)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWFD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWFD),L(NWFD),MAX(NWFD),N(NWFD)
      DEN = L(I) + L(J) + 3+ K
      FACT = (D1/(L(I)+1) + D1/(L(J)+1))/(DEN + D1)
      A = EH**K
      A2 = A*A
      H90 = H/90.D0
      A3 = A2*A*H90
      AI = H90/A
      AN = 114.D0*A*H90
      A34 = 34.D0*H90
      F1 = RR(1)*P(1,I)*P(1,J)
      F2 = RR(2)*P(2,I)*P(2,J)
      F3 = RR(3)*P(3,I)*P(3,J)
      F4 = RR(4)*P(4,I)*P(4,J)
      YK(1) = F1*(D1 + Z*R(1)*FACT)/DEN
      YK(2) = F2*(D1 + Z*R(2)*FACT)/DEN
      YK(3) = YK(1)*A2 + H3*(F3 + D4*A*F2 + A2*F1)
      MX = (MIN0(MAX(I),MAX(J))/2)*2
      DO 8 M = 5,MX
      F5 = (RR(M)*P(M,I))*P(M,J)
      YK(M-1) = YK(M-3)*A2 + ( AN*F3 + A34*(F4+A2*F2)-F5*AI-F1*A3)
      F1 = F2
      F2 = F3
      F3 = F4
8     F4 = F5
      M1 = MX - 1
      IF (IABS(I-J)  +  IABS(K) .NE. 0) GO TO 2
*
*  *****  FOR Y0(I,I) SET THE LIMIT TO 1 AND REMOVE OSCILLATIONS
*  *****  INTRODUCED BY THE USE OF SIMPSON'S RULE
*
      M2 = M1 - 1
      C1 = D1 - YK(M1)
      C2 = D1 - YK(M2)
      DO 3 M = 1,M1,2
      YK(M) = YK(M) + C1
3     YK(M+1) = YK(M+1) + C2
2     DO 1 M = M1+1,NO
         YK(M) = A*YK(M-1)
1     CONTINUE
      END
*
*     ------------------------------------------------------------------
*                       L O O K - T M - D F
*     ------------------------------------------------------------------
*
      SUBROUTINE LOOKTMDF(SL, IP, NSL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*2 SL, TERMS(10)
*             .. d1 f1 terms
      DATA   TERMS/'1P','1D','1F','1G','1H','3P','3D','3F','3G','3H'/
*
      IP = 0
      NSL = 0
1     I = 1
2     IF (SL .EQ. TERMS(I)) THEN
        NSL = NSL + 1
        IP = I
      END IF
      I = I+1
      IF (I .LE. 10) GO TO 2
      RETURN
      END
*     ------------------------------------------------------------------
*                       L O O K - F
*     ------------------------------------------------------------------
*
*     Add the deviations to the average energy for a partially filled
*      f- hell
*
      SUBROUTINE LOOKF( L, SL, SEN, N, IP, NSL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      CHARACTER SL*2, SEN*1
      INTEGER IPTR(4)
      CHARACTER*3 TER,TERMS(144),TERMS6(119),TERMS7(119)
      DATA    IPTR/7,24,71,144/
      DATA       TERMS/'3F1','3P1','3H1','1D1','1G1','1S1','1I1',
*             .. f3 terms
     :     '4S1','4D1','4F1','4G1','4I1','2P1','2D1','2D2','2F1',
     :     '2F2','2G1','2G2','2H1','2H2','2I1','2K1','2L1',
*            ... f4 terms ...
     :       '5D1','5F1','5G1','5S0','5I1','3F1','3F2','3D1',
     :       '3D2','3F3','3G1','3G2','3F4','3G3','3P1','3P2',
     :       '3P3','3H1','3H2','3H3','3H4','3I1','3I2','3K1',
     :       '3K2','3L1','3M1','1D1','1D2','1D3','1F1','1G1',
     :       '1G2','1G3','1D4','1G4','1H1','1H2','1S1','1I1',
     :       '1S2','1I2','1I3','1K1','1L1','1L2','1N1',
*            ... f5 terms ...
     :       '6P0','6F0','6H0','4S1','4P1','4P2','4D1','4D2',
     :       '4D3','4F1','4F2','4F3','4F4','4G1','4G2','4G3',
     :       '4G4','4H1','4H2','4H3','4I1','4I2','4I3','4K1',
     :       '4K2','4L1','4M0','2P1','2P2','2P3','2P4','2D1',
     :       '2D2','2D3','2D4','2D5','2F1','2F2','2F3','2F4',
     :       '2F5','2F6','2F7','2G1','2G2','2G3','2G4','2G5',
     :       '2G6','2H1','2H2','2H3','2H4','2H5','2H6','2H7',
     :       '2I1','2I2','2I3','2I4','2I5','2K1','2K2','2K3',
     :       '2K4','2K5','2L1','2L2','2L3','2M1','2M2','2N1',
     :       '2O0'/
*            ... f6 terms ...
      DATA TERMS6/'7F0','5D1','5D2','5D3','5F1','5F2','5G1','5G2',
     :       '5G3','5P0','5H1','5H2','5S0','5I1','5I2','5K0',
     :       '5L0','3F1','3F2','3F6','3F8','3D1','3D2','3D3',
     :       '3D4','3F3','3F5','3G1','3G2','3G4','3G5','3D5',
     :       '3F4','3F7','3F9','3G3','3G6','3G7','3P1','3P2',
     :       '3P3','3H1','3H2','3H3','3H4','3P4','3H5','3H6',
     :       '3P5','3P6','3H7','3H8','3H9','3I1','3I2','3I3',
     :       '3I4','3I5','3I6','3K1','3K2','3K3','3K4','3K5',
     :       '3K6','3L1','3L2','3L3','3M1','3M2','3M3','3N0',
     :       '3O0','1F2','1F3','1F4','1D1','1D2','1D3','1F1',
     :       '1G1','1G2','1G3','1D5','1G5','1D6','1G6','1G7',
     :       '1G8','1D4','1G4','1H1','1H2','1P0','1H3','1H4',
     :       '1S1','1I1','1S2','1I2','1I3','1S3','1I4','1I5',
     :       '1S4','1I6','1I7','1K1','1K2','1K3','1L1','1L2',
     :       '1L3','1L4','1M1','1M2','1N1','1N2','1Q0'/
*            ... f7 terms ...
      DATA TERMS7/'8S0','6P0','6D0','6F0','6G0','6H0','6I0','4S1',
     :       '4S2','4P1','4P2','4D1','4D2','4D3','4D4','4D5',
     :       '4D6','4F1','4F2','4F3','4F4','4F5','4G1','4G2',
     :       '4G3','4G4','4G5','4G6','4G7','4H1','4H2','4H3',
     :       '4H4','4H5','4I1','4I2','4I3','4I4','4I5','4K1',
     :       '4K2','4K3','4L1','4L2','4L3','4M0','4N0','2S1',
     :       '2S2','2P1','2P2','2P3','2P4','2P5','2D1','2D2',
     :       '2D3','2D4','2D5','2D6','2D7','2F1','2F2','2F3',
     :       '2F4','2F5','2F6','2F7','2F8','2F9','2FA','2G1',
     :       '2G2','2G3','2G4','2G5','2G6','2G7','2G8','2G9',
     :       '2GA','2H1','2H2','2H3','2H4','2H5','2H6','2H7',
     :       '2H8','2H9','2I1','2I2','2I3','2I4','2I5','2I6',
     :       '2I7','2I8','2I9','2K1','2K2','2K3','2K4','2K5',
     :       '2K6','2K7','2L1','2L2','2L3','2L4','2L5','2M1',
     :       '2M2','2M3','2M4','2N1','2N2','2O0','2Q0'/

*
*  --- search for a partially unfilled f- shell
*
      IF (L .EQ. 3) THEN
         IF (N .EQ. 2) THEN
            IBEGIN = 1
            IEND = 7
         ELSE IF (N .EQ. 6) THEN
            IBEGIN = 1
            IEND = 119
         ELSE IF (N .EQ. 7) THEN
            IBEGIN = 1
            IEND = 119
         ELSE
            IBEGIN = IPTR(N-2) + 1
            IEND = IPTR(N-1)
         END IF
1        I = IBEGIN
10       IF (N .EQ. 6) THEN
            TER = TERMS6(I)
         ELSEIF (N .EQ. 7) THEN
            TER = TERMS7(I)
         ELSE
            TER = TERMS(I)
         END IF
         IF (SL .EQ. TER(1:2)) THEN
            IF (SEN .EQ. ' ' .OR. SEN .EQ. TER(3:3)) THEN
               NSL = NSL + 1
               IP = I
               IF (N .EQ. 6) IP = I + 144
               IF (N .EQ. 7) IP = I + 263
            END IF
         END IF
         I = I+1
         IF (I .LE. IEND) GO TO 10
      ELSE IF ( N .EQ. 1 .AND. SL(1:1) .EQ. '2') THEN
         NSL = 0
      END IF
      RETURN
      END
*
*     ------------------------------------------------------------------
*                       D E V F
*     ------------------------------------------------------------------
*
*     Add the deviations to the average energy for a partially filled
*       f- shell
*
      SUBROUTINE DEVF(IEL, L, N, I, DONE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL DONE
      INTEGER F2FF2(7), F4FF2(7), F6FF2(7),F2FF3(17),F4FF3(17),
     :F6FF3(17),F2FF4(47),F4FF4(47),F6FF4(47),F2FF5(73),F4FF5(73),
     :F6FF5(73),F2VV5(73),F4VV5(73),F6VV5(73),F2FF6(119),F4FF6(119),
     :F6FF6(119),F2VV6(119),F4VV6(119),F6VV6(119),F2FF7(119),F4FF7(119),
     :F6FF7(119),F2VV7(119),F4VV7(119),F6VV7(119)
      DATA F2FF2/-2092090,19277115,-7920055,9175309,-9862710,25105080,
     :11506495/
      DATA F4FF2/-1426425,3871725,-2871375,-6724575,9009325,17117100,
     :1945125/
      DATA F6FF2/-1828750,-13715625,1413125,10058125,2493750,
     :21945000,1579375/
*             ... f3 coefficients
      DATA F2FF3/-6936930,16681665,-6936930,1651650,-21966945,
     :-4789785,-6139419,12019293,23123100,36005970,7172880,10731006,
     :5230225,-13092079,3798795,-11231220,5945940/
      DATA F4FF3/-4729725,1126125,-4729725,-2600325,-8456175,150150,
     :4890600,-2372175,15765750,2190825,-1497600,13953225,2634450,
     :-2777775,3521700,4142775,-1535625/
      DATA F6FF3/-6063750,-19201875,-6063750,-10841250,2296875,7074375,
     :3320625,3661875,20212500,-8452500,4515000,-10395000,1500625,
     :12464375,1194375,5696250,4961250/
*             ... f4 coefficients
      DATA F2FF4/74340,-105840,-40320,-105840,-220500,123480,-3920,
     :-54828,216612,68880,53280,-29088,-9016,71064,159516,-59976,
     :234864,113652,11368,11732,81984,-23940,107100,44856,-145152,
     :-138600,-89460,362628,23445,-75420,156240,110376,24264,307872,
     :274995,5544,-17388,284004,352800,146412,549360,71883,50589,
     :-59976,36120,-136416,8820/
      DATA F4FF4/-55440,-145530,-112770,-145530,-202860,169785,
     :-5390,73260,-57195,-10395,-71430,81125,-56350,98350,187803,
     :-33418,180250,164871,2254,127925,-74550,10080,18270,-120505,
     :43750,-18585,-80010,151767,33363,234135,243495,373443,71037,
     :-6105,202545,-45885,-84105,83790,485100,181251,67410,72324,
     :239400,93345,49665,7350,7350/
      DATA F6FF4/-1732500,-831600,-1159200,-831600,-258300,970200,
     :-30800,-188100,-492300,-231000,57600,-297000,186200,-667800,
     :790020,249480,-1486800,1019340,-107240,-605500,-16800,-31500,
     :-459900,277200,201600,428400,522900,1150380,564795,69300,
     :-667800,1035720,413280,-138600,-864675,932400,1152900,-422100,
     :2772000,1687140,-1159200,264285,-332325,579600,462000,1150800,
     :711900/
*             ... f5 coefficients
      DATA F2FF5/1,-28,-179,14,74,-436,14,256,443,14,
     :-70,-134,5452,14,-538,874,344,-506,5881,-436,
     :14,-4,-101,-1814,578,-38,-17,22,131,5939,
     :1789,1738,3064,74,4742,10487,224,37,196,2543,
     :5413,1676,929,128,1639,8,14251,1241,821,46,
     :2996,-37,-13426,3209,269,43627,2*14,16,-406,
     :2831,14,6566,422,174706,-3776,14,-14,-5672,-5,
     :197,-583,-88/
      DATA F2VV5/195,117,585,195,1755,8775,195,6825,20475,195,
     :1053,1755,26325,195,4095,61425,8775,5265,26325,8775,
     :195,65,585,8775,8775,195,65,65,702,35100,
     :11700,6825,6825,6825,45045,32175,585,117,5265,19305,
     :26325,3861,8775,455,20475,4095,61425,8775,2925,351,
     :8775,10530,289575,35100,70200,154440,2*195,585,14625,
     :14625,195,96525,8775,740025,67275,195,143,32175,39,
     :2925,2925,585/
      DATA F4FF5/-4,-70,-848,7,-71,-331,7,-269,942,7,
     :-175,-94,2339,7,223,-226,-2705,-2237,-2458,485,
     :7,-289,-79,584,-3953,-536,-31,206,241,3989,
     :5195,192,355,4,7177,188,112,335,98,5557,
     :3197,652,283,3242,163,3226,4049,1702,1954,1538,
     :2149,59,19957,1633,97,7057,168,140,89,553,
     :-16,2149,-22631,337,-107135,24194,7,2177,3334,-97,
     :397,-932,-40/
      DATA F4VV5/39,429,4719,143,3861,3861,143,3003,11011,143,
     :3861,4719,42471,143,99099,2457,42471,42471,42471,14157,
     :143,4719,1573,42471,42471,4719,429,1287,3861,84942,
     :28314,1001,3003,3003,121121,51909,429,4719,3861,51909,
     :42471,17303,4719,33033,819,99099,297297,42471,14157,14157,
     :14157,3861,467181,28314,1716,207636,1573,4719,3146,9438,
     :4719,14157,467181,42471,3581721,325611,143,51909,51909,9438,
     :9438,14157,1573/
      DATA F6FF5/-1925,-3500,-31675,350,-5950,-1400,350,-100,-37525,
     :350,
     :-8750,-42700,-396200,350,-14450,-200,-42700,-68950,-239225,
     :-74200,
     :350,-2800,-3325,-35000,-33250,2800,175,350,-5075,-68425,
     :-107975,3*350,-41300,-91525,5600,27475,4900,-435575,
     :-113225,-679000,-11725,2800,175,7000,-12425,-8575,-27475,
     :98350,
     :38150,16625,824600,-31675,6475,-2796325,12250,350,
     :6125,4585,
     :1015,9800,532700,24500,687400,86800,1400,42700,141400,
     :2450,1400,35525,25900/
      DATA F6VV5/5577,16731,184041,5577,50193,50193,5577,1859,184041,
     :5577,
     :150579,552123,1656369,5577,184041,4563,552123,1656369,1656369,
     :552123,
     :5577,61347,184041,552123,552123,61347,5577,5577,100386,2208492,
     :736164,3*5577,2024451,2024451,16731,184041,150579,6073353,
     :1656369,6073353,552123,20449,1521,184041,552123,552123,184041,
     :552123,
     :552123,301158,18220059,2208492,401544,48586824,61347,20449,
     :184041,184041,
     :184041,61347,6073353,552123,46562373,4232943,5577,674817,2024451,
     :20449,184041,184041,184041/
*             ... f6 coefficients
      DATA F2FF6/-14,-227,-587,-223,-28,-16,-136,-2174,332,-41,
     :-61,-1031,-28,7,-331,-16,-184,154,476,3184,
     :1546,157,2713,983,18869,488,-14,128,3524,142,
     :-188,3023,11558,3407,1751,2458,2042,614,209,8672,
     :118,809,664,-961,274,359,-187,21721,3797,1069,
     :659,1387,1013,-1,59,-41,53,1219,-47,3602,
     :976,-46048,-1741,-77111,6686,-58,-1288,674,-97,-571,
     :-193,-647,-158,98,1166,1238,401,94541,961,8,
     :274,13648,17236,1543,974,287,692,-28,358066,499,
     :5716,1511,209,887,-23,341,112,679,4,1037,
     :233,98,259,8,1166,4724,322903,646,484,854,
     :-58,2336,3382,388,-16,-452,-35,-121,-82/
      DATA F2VV6/39,1755,61425,2275,585,117,1755,12285,6825,585,
     :351,8775,585,1755,1755,117,585,585,5265,19305,
     :19305,975,8775,20475,135135,1755,1053,585,26325,4095,
     :184275,32175,26325,26325,8775,26325,26325,2925,975,26325,
     :26325,2925,26325,26325,975,5265,5265,289575,26325,2925,
     :8775,8775,19305,39,1755,585,3510,9750,4875,26325,
     :26325,289575,26325,740025,67275,1755,19305,10725,1755,3510,
     :1950,2925,585,1755,8775,2925,585,245700,4095,195,
     :975,61425,61425,8775,19305,975,2925,15795,868725,2340,
     :8775,2925,975,2925,2925,2925,195,2925,13,35100,
     :11700,1755,19305,195,8775,30225,997425,8775,2925,8775,
     :1755,8775,90675,3627,117,2925,351,8775,325/
      DATA F4FF6/-35,-68,-1087,-688,-14,-76,-566,-6133,-1412,-107,
     :-1541,-1216,-14,-280,-1775,-2230,-466,7,238,-400,
     :15601,64,97,71,12916,379,-35,422,179,-1297,
     :16064,5755,446,772,370,1514,95,-85,359,778,
     :1202,4313,878,4867,158,43,-619,-11788,1747,-631,
     :2*-85,9257,236,370,-89,-19,-379,467,-293,
     :166,16942,-86,261440,-17321,955,-11270,-4255,58,-1457,
     :-373,-397,-323,49,1501,5,137,6824,9335,769,
     :23027,67511,31,487,10453,161,227,30695,-21961,652,
     :2411,17,262,461,2*213,56,1351,-34,5191,
     :899,49,3017,314,3139,3522,-50785,2021,-787,301,
     :773,2242,2417,311,106,174,266,-193,-193/
      DATA F4VV6/143,1287,9009,33033,429,4719,14157,99099,9009,1287,
     :14157,14157,429,14157,14157,14157,4719,39,3861,17303,
     :155727,429,14157,3003,1090089,4719,3861,14157,1287,99099,
     :297297,51909,3861,42471,14157,14157,3861,14157,2145,6435,
     :14157,23595,19305,42471,4719,3861,42471,467181,42471,14157,
     :2*14157,155727,4719,14157,4719,2574,15730,7865,14157,
     :1287,467181,42471,3581721,325611,14157,155727,51909,14157,28314,
     :9438,14157,4719,1287,14157,1573,715,45045,33033,4719,
     :70785,495495,9009,14157,155727,1573,14157,127413,1401543,4719,
     :14157,143,4719,14157,2*1573,143,7865,1573,70785,
     :4719,1287,155727,4719,14157,48763,1609179,14157,14157,4719,
     :14157,14157,48763,48763,4719,1573,14157,4719,4719/
      DATA F6FF6/-1750,175,-8275,-4575,-700,-39550,-14000,-88750,
     :-850,-2975,
     :-66325,-88025,-700,-39025,-55475,-16450,-15050,350,11900,-358750,
     :-850150,-125,29875,-1025,-495575,-350,-1750,11800,-50,
     :-4450,
     :-118850,-312025,-12950,-10675,-63175,62650,-10850,-16450,1435,
     :-8120,
     :144550,40985,17360,123725,2450,-8575,1225,-491225,
     :-180775,-19775,
     :-24325,-39725,-574175,2625,71575,2975,-175,-5495,-1505,226100,
     :5600,-52150,-175,-193025,-43750,4900,379750,18550,
     :98875,6475,
     :5425,13825,23450,2450,-20650,-8050,3395,5185,-7375,8050,
     :13510,51910,9650,6475,150850,-5425,1400,-555800,-1272950,
     :70175,
     :-44450,175,10325,-5425,-2975,-10675,2800,68635,4200,474005,
     :83825,2450,483875,1750,-5950,-29050,1289575,76300,17500,
     :6650,
     :139300,32900,255850,340900,23450,10150,169225,78575,350/
      DATA F6VV6/5577,2*50193,20449,16731,184041,552123,552123,
     :5577,16731,
     :2*552123,16731,2*552123,2*184041,1521,150579,6073353,
     :6073353,5577,552123,16731,6073353,552123,150579,184041,11583,
     :184041,
     :1656369,2024451,150579,127413,552123,1656369,150579,184041,5577,
     :150579,
     :1656369,184041,150579,1656369,61347,150579,1656369,18220059,
     :1656369,184041,
     :2*552123,6073353,20449,552123,184041,7722,122694,20449,1656369,
     :150579,18220059,1656369,3581721,4232943,42471,6073353,674817,
     :552123,84942,
     :122694,2*184041,50193,552123,184041,16731,200772,184041,61347,
     :61347,552123,50193,552123,6073353,61347,184041,4969107,54660177,
     :736164,
     :552123,16731,61347,3*184041,5577,184041,20449,2208492,
     :736164,50193,6073353,61347,42471,1901757,62757981,552123,184041,
     :552123,
     :2*552123,2*5705271,2*184041,2*552123,1573/
*             ... f7 coefficients
      DATA F2FF7/-98,-17,-937,-112,-4,-73,-203,14,-34,1,
     :1793,67,5323,-3517,-2629,919,17,14,-98,418,
     :-7328,14,158,-82,6038,688,-109,-3322,-3,-137,
     :-4559,-991,-391,-157,301,-73,109,-503,661,-638,
     :-274,-74,2,-362,-926,-11,-59,28,4144,821,
     :2,241,557,959,22157,11407,9067,13073,15277,517,
     :779,28,58,448,542,3163,3220,4079,28,1804,
     :272,
     :4304,4282,776,3022,662,-274,1556,374,7246,
     :38924,55,
     :4391,526,73559,1957,13807,11801,623,779,29,-73,
     :179,8453,3661,-175,44,12488,-9881,236,1604,43,
     :102031,7504,-1342,272,44,-108,1136,3128,-638,-9,
     :-577,-4,272,467,-191,-164,-872/
      DATA F2VV7/195,39,2925,585,39,585,585,117,585,117,
     :2925,1755,20475,20475,122850,4095,130,117,5265,1755,
     :26325,117,1755,4095,20475,2925,2457,20475,325,5265,
     :26325,8775,2925,2925,1755,585,585,3510,5850,2925,
     :2925,975,65,1755,8775,65,195,1755,8775,1755,
     :65,650,5850,8775,61425,20475,20475,45045,32175,8775,
     :2925,65,195,5265,19305,26325,19305,8775,1755,8775,
     :2925,
     :12285,61425,4095,20475,2925,2925,19305,8775,15795,
     :868725,351,
     :8775,5265,289575,17550,35100,77220,2925,2925,351,585,
     :1170,29250,14625,3861,585,90675,332475,1755,10725,975,
     :740025,67275,8775,2925,585,715,10725,90675,18135,130,
     :1950,39,2925,2925,975,585,2925/
      DATA F4FF7/-49,-82,-10,-56,-3538,-538,-266,35,541,-23,
     :4219,79,-41,1928,-545,-4546,-205,35,-49,124,
     :-1891,-434,1051,5017,-44,241,-6046,514,284,-617,
     :-5752,557,34,73,1337,-265,281,-1993,-1193,1346,
     :-4751,-1588,-58,-1172,-1237,-45,-1237,14,-574,274,
     :172,4999,2689,476,2152,5477,388,45739,3182,778,
     :320,42,137,224,1767,289,-406,267,14,128,
     :736,
     :13480,1957,9892,4001,3478,934,-542,398,-12290,
     :224416,1922,
     :3001,2620,64453,1181,547,2347,-84,-4,1844,8,
     :757,7451,742,602,-31,4655,218216,2663,-3485,107,
     :-119963,44492,424,206,111,857,8044,-30,2076,79,
     :-257,-239,-149,622,284,-58,-304/
      DATA F4VV7/143,429,143,429,14157,4719,1573,429,4719,3861,
     :42471,1287,3003,33033,18018,33033,9438,429,3861,1573,
     :42471,4719,14157,99099,27027,42471,99099,9009,14157,42471,
     :42471,14157,4719,4719,14157,4719,4719,28314,9438,42471,
     :42471,14157,1573,14157,14157,1573,14157,1287,14157,1287,
     :3861,42471,14157,4719,9009,33033,3003,363363,51909,14157,
     :4719,143,1573,3861,17303,3861,51909,1573,1287,14157,
     :4719,
     :99099,9009,99099,297297,42471,14157,155727,4719,127413,
     :1401543,14157,
     :14157,42471,467181,14157,3146,103818,1573,429,14157,4719,
     :9438,47190,23595,155727,4719,146289,1609179,14157,467181,3861,
     :3581721,325611,4719,14157,1573,17303,51909,48763,48763,3146,
     :9438,4719,1573,14157,14157,1573,4719/
      DATA F6FF7/-2450,-175,-5075,-2800,-13300,-37625,-37975,1750,
     :-44450,-175,
     :-20125,7525,-1825,-5125,-12575,-20725,-7525,1750,-2450,-101500,
     :164500,-21700,66850,-13250,-200,-21700,-49825,-200,-10675,
     :-20125,
     :239575,-16975,-22225,-28525,41825,10675,-27475,-53725,-31325,700,
     :5950,-700,-700,-700,-25550,2975,1225,700,-51800,2975,
     :1750,-27475,-16975,-53725,3625,13075,-1825,-110975,-9275,-21875,
     :-22225,700,14000,11200,330400,7175,639800,-6475,700,
     :-26600,
     :-25900,92200,9500,400,1900,350,23450,30100,-33250,
     :-197050,
     :-7008400,
     :128275,27125,109550,-972475,47425,-291725,2017225,875,-875,
     :144725,35875,
     :11725,-2345,7805,247625,1750,-260050,-1702925,106750,
     :224350,175,
     :6418825,-1400,350,-2800,58450,108850,-4200,102200,
     :17150,19075,875,6650,15050,875,1575,51800,34300/
      DATA F6VV7/2*5577,2*16731,61347,2*184041,16731,
     :184041,5577,
     :61347,50193,16731,184041,100386,184041,40898,16731,150579,552123,
     :1656369,184041,552123,184041,1521,184041,552123,1287,61347,
     :1656369,
     :1656369,552123,2*184041,552123,2*184041,1104246,368082,61347,
     :61347,20449,61347,42471,552123,61347,61347,50193,552123,50193,
     :16731,2*368082,552123,50193,184041,16731,2024451,155727,552123,
     :184041,1859,61347,150579,6073353,150579,6073353,42471,50193,
     :552123,
     :184041,552123,50193,184041,20449,61347,184041,6073353,552123,
     :4969107,
     :54660177,
     :2*552123,1656369,18220059,1104246,2208492,24293412,184041,16731,
     :552123,184041,
     :368082,28314,184041,6073353,184041,5705271,20919327,552123,
     :2024451,1521,
     :46562373,325611,552123,2*184041,674817,224939,5705271,
     :5705271,122694,3146,61347,184041,14157,20449,2*184041/
      DONE = .TRUE.
      IF (L .EQ. 3) THEN
         IF (N .EQ. 2) THEN
            CALL ADD(2*F2FF2(I)/87419475.D0,2,IEL,IEL,.TRUE.)
            CALL ADD(2*F4FF2(I)/87419475.D0,4,IEL,IEL,.TRUE.)
            CALL ADD(2*F6FF2(I)/87419475.D0,6,IEL,IEL,.TRUE.)
         ELSE IF (N.EQ. 3) THEN
            J=I-7
            CALL ADD(2*F2FF3(J)/96621525.D0,2,IEL,IEL,.TRUE.)
            CALL ADD(2*F4FF3(J)/96621525.D0,4,IEL,IEL,.TRUE.)
            CALL ADD(2*F6FF3(J)/96621525.D0,6,IEL,IEL,.TRUE.)
         ELSE IF (N.EQ. 4) THEN
            J=I-24
            CALL ADD(2*F2FF4(J)/737100.D0,2,IEL,IEL,.TRUE.)
            CALL ADD(2*F4FF4(J)/1486485.D0,4,IEL,IEL,.TRUE.)
            CALL ADD(2*F6FF4(J)/6625476.D0,6,IEL,IEL,.TRUE.)
         ELSE IF (N.EQ. 5) THEN
            J=I-71
            VV=F2VV5(J)*1.D0
            CALL ADD(2*F2FF5(J)/VV,2,IEL,IEL,.TRUE.)
            VV=F4VV5(J)*1.D0
            CALL ADD(2*F4FF5(J)/VV,4,IEL,IEL,.TRUE.)
            VV=F6VV5(J)*1.D0
            CALL ADD(2*F6FF5(J)/VV,6,IEL,IEL,.TRUE.)
         ELSE IF (N.EQ. 6) THEN
            J=I-144
            VV=F2VV6(J)*1.D0
            CALL ADD(2*F2FF6(J)/VV,2,IEL,IEL,.TRUE.)
            VV=F4VV6(J)*1.D0
            CALL ADD(2*F4FF6(J)/VV,4,IEL,IEL,.TRUE.)
            VV=F6VV6(J)*1.D0
            CALL ADD(2*F6FF6(J)/VV,6,IEL,IEL,.TRUE.)
         ELSE IF (N.EQ. 7) THEN
            J=I-263
            VV=F2VV7(J)*1.D0
            CALL ADD(2*F2FF7(J)/VV,2,IEL,IEL,.TRUE.)
            VV=F4VV7(J)*1.D0
            CALL ADD(2*F4FF7(J)/VV,4,IEL,IEL,.TRUE.)
            VV=F6VV7(J)*1.D0
            CALL ADD(2*F6FF7(J)/VV,6,IEL,IEL,.TRUE.)
         ELSE
            DONE = .FALSE.
         END IF
      ELSE
         DONE = .FALSE.
      END IF
      RETURN
      END
*
*     ------------------------------------------------------------------
*                       D E V D F
*     ------------------------------------------------------------------
*
*     Add the deviations to the average energy for a partially filled
*       p- or d- shell
*
      SUBROUTINE DEVDF(IED, IEF, I, DONE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL DONE
      INTEGER F1DF(10), F2DF(10), F3DF(10), F4DF(10), F5DF(10)
      DATA    F1DF/5,-3,15,-17,33,1,9,-9,23,-27/
      DATA    F2DF/24,6,-11,-15,10,24,6,-11,-15,10/
      DATA    F3DF/30,-36,25,41,16,-18,48,-13,-29,-4/
      DATA    F4DF/66,-99,66,-22,3,66,-99,66,-22,3/
      DATA    F5DF/1815,990,440,220,170,-1485,-660,-110,110,160/
      DONE = .TRUE.
      CALL ADD(F2DF(I)/105.D0,2,IED,IEF,.TRUE.)
      CALL ADD(F2DF(I)/105.D0,2,IEF,IED,.TRUE.)
      CALL ADD(F4DF(I)/693.D0,4,IED,IEF,.TRUE.)
      CALL ADD(F4DF(I)/693.D0,4,IEF,IED,.TRUE.)
      CALL ADD(F1DF(I)/70.D0,1,IED,IEF,.FALSE.)
      CALL ADD(F1DF(I)/70.D0,1,IEF,IED,.FALSE.)
      CALL ADD(F3DF(I)/315.D0,3,IED,IEF,.FALSE.)
      CALL ADD(F3DF(I)/315.D0,3,IEF,IED,.FALSE.)
      CALL ADD(F5DF(I)/7623.D0,5,IED,IEF,.FALSE.)
      CALL ADD(F5DF(I)/7623.D0,5,IEF,IED,.FALSE.)
      RETURN
      END
