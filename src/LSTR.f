*     ------------------------------------------------------------------
*
*       LSTR -- A GENERAL PROGRAM FOR ALLOWED (E1,E2) TRANSITIONS
*
*                  C O P Y R I G H T -- 1994
*
*       by C. Froese Fischer
*          Vanderbilt University
*          Nashville, TN 37235 USA
*
*          K.M.S. Saxena, 
*          Department of Theoretical Chemistry
*          University of Alberta
*
*       and
*
*          M.R. Godefroid
*          Universite Libre de Bruxelles
*          Bruxelles, 1050 Belgium
*
*     Computer Physics Communications, Vol.  9, 381 (1975).
*                                      Vol. 14, 275 (1978).
*                                      Vol. 17, 427 (1979).
*     ------------------------------------------------------------------
*
      PROGRAM TRANS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER (IREAD=5,IWRITE=6)
      PARAMETER(NOD=220,NWD=30,NWD2=2*NWD,NCD=100,NCD2=2*NCD)
      INTEGER Q(5),NW(2)
      CHARACTER*1 END,ASTER,PP
      CHARACTER*3 EL(75),ELC(8),COUPLE(9)
      CHARACTER*6 ATOM,ALABEL(2),AT,TT
      CHARACTER*24 CFILE(2),WFILE(2),DFILE
      CHARACTER*66 CONFIG(NCD2)
      CHARACTER*13 LABEL(2)
      CHARACTER*72 HEADER
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD2),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWD2),L(NWD2),MAX(NWD2),N(NWD2)
      COMMON /STATE/WT(NCD2),ET(2),JV(2),NCF(2),NJ(2)
      COMMON /MULT/OVINT(10,10),RLINT(10,10),RVINT(10,10)
      DIMENSION IP(3),IPT(8), LMULT(2), NOCCSH(NCD2), ISMULT(2)
      LOGICAL PRINT
      DATA ASTER/'*'/, LABEL/'Initial State','Final   State'/
*
*  *****  DETERMINE DATA
*
999   WRITE(IWRITE,998)
998   format(//20x,'============================'/
     :         20X,'  LS - T R A N S I T I O N '/
     :         20X,'============================'/)
*
CSUN  i = iargc()
1     continue
CSUN  if (i .lt. 2) then
         WRITE(0,*) ' Name of Initial State'
         read(iread,'(A)') cfile(1)
         WRITE(0,*) ' Name of Final State'
         read(iread,'(A)') cfile(2)
CSUN  else
CSUN     call getarg(1,cfile(1))
CSUN     call getarg(2,cfile(2))
CSUN  end if
CSUN  if (i .eq. 3) then
CSUN     call getarg(3,dfile)
CSUN  else
         dfile = 'mltpol.lst'
CSUN  end if
      do 2 i = 1,2
         j = index(cfile(i),' ')
         if (j .eq. 1) then
            WRITE(0,*) ' Names may not start with blanks'
            go to 1
         end if
         cfile(i) = cfile(i)(1:j-1)//'.c'
         wfile(i) = cfile(i)(1:j-1)//'.w'
2     continue
*
*     Some systems have options for having more than one unit
*     number connected or for appending. Default is without.
*     OPEN(UNIT=1,FILE=WFILE(1),STATUS='OLD',FORM='UNFORMATTED',
*    :     READONLY)
*     OPEN(UNIT=2,FILE=WFILE(2),STATUS='OLD',FORM='UNFORMATTED',
*    :     READONLY)
*     OPEN(UNIT=7,FILE=CFILE(1),STATUS='OLD',READONLY)
*     OPEN(UNIT=8,FILE=CFILE(2),STATUS='OLD',READONLY)
      OPEN(UNIT=1,FILE=WFILE(1),STATUS='OLD',FORM='UNFORMATTED')
      OPEN(UNIT=2,FILE=WFILE(2),STATUS='OLD',FORM='UNFORMATTED')
      OPEN(UNIT=7,FILE=CFILE(1),STATUS='OLD')
      OPEN(UNIT=8,FILE=CFILE(2),STATUS='OLD')
      OPEN(UNIT=9,FILE=DFILE,STATUS='OLD')
      WRITE(0,*) ' INTERMEDIATE PRINTING (Y OR N) ?  '
      READ(IREAD,'(A1)') PP
      IF ( PP .EQ. 'Y' .OR. PP .EQ. 'y') THEN
         PRINT = .TRUE.
         WRITE(0,*) '  TOLERANCE FOR PRINTING ?  '
         READ(IREAD,'(F6.4)') TOL
      ELSE
         PRINT = .FALSE.
      ENDIF
      DO 1000 III = 1,2
      READ(9,'(A)') HEADER
      I = 1
      M = 1
      DO 6 II = 1,2
*
*  *****  READ THE RADIAL FUNCTIONS DURING THE FIRST CYCLE
*
      IF ( III .EQ. 1 ) THEN
12    READ(II,END=13) AT,TT,EL(I),MR,Z,ETI,EKI,AZ(I),
     :  (P(JJ,I),JJ=1,MR)
      IF (I.GT.(75)) STOP ' TOO MANY ELECTRONS: MAX=(75)'
      IF (EL(I)(1:1) .NE. ' ') THEN
         N(I) = ICHAR(EL(I)(1:1)) - ICHAR('1')
         L(I) = LVAL(EL(I)(2:2))
      ELSE
         N(I) = ICHAR(EL(I)(2:2)) - ICHAR('1')
         L(I) = LVAL(EL(I)(3:3))
      ENDIF
      MM = MR+1
      DO 24 J = MM,NO
24    P(J,I) = D0
 8    MAX(I)=MR
      I = I+1
      GO TO 12
   13 NW(II) = I-1
      ENDIF
*
*  *****  READ THE CONFIGURATIONS
*
      IU = II+6
      READ(IU,10,END=1000) ATOM,ALABEL(II),ET(II)
   10 FORMAT(3X,2A6,F14.7/)
      WRITE(IWRITE, 11) LABEL(II),ATOM,ALABEL(II),ET(II)
   11 FORMAT(/3X,A/3X,13('-')/3X,2A6,F16.7/)
   15 READ(IU,14,END=18) (ELC(K),Q(K),K=1,5),WT(M)
   14 FORMAT(5(1X,A3,1X,I2,1X),F10.7)
      READ(IU,'(9(5X,A3))',END=18) (COUPLE(J),J=1,9)
      NOCC = 0
16    IF (ELC(NOCC+1) .NE. '   ' ) THEN
         NOCC = NOCC+1
         IF (NOCC .LT. (5)) GO TO 16
      END IF
      IF (NOCC .EQ. 0) GO TO 18
      MM = M
      IF (II.EQ.2) MM = M - NCF(1)
      NOCCSH(M) = NOCC
      CALL PACK(NOCC,ELC,Q,COUPLE,CONFIG(M))
      K = 66
   17 IF (CONFIG(M)(K:K) .EQ. ' ') THEN
         K = K-1
         GO TO 17
      END IF
      WRITE(IWRITE,'(I4,3X,F10.7,3X,A)') MM,WT(M),CONFIG(M)(1:K)
      M = M+1
      IF (M.GT.(NCD2)) THEN
         WRITE(0,*)
     :       ' TOO MANY CONFIGURATIONS IN THE INITIAL AND FINAL STATE',
     :       ' MAXIMUM FOR COMBINED SUM = ',(NCD2)-1
      END IF
      GO TO 15
   18 NCF(II) = M-1
      DO 5 J = 6,1,-1
         IF (ALABEL(II)(J:J) .NE. ' ') THEN
             LMULT(II) = 2*LVAL(ALABEL(II)(J:J)) + 1
             ISMULT(II) = ICHAR(ALABEL(II)(J-1:J-1))-ICHAR('0')
             GO TO 6
         ENDIF
    5 CONTINUE
      WRITE(IWRITE,'(A)') 'No term value found in ',ALABEL(II)
      STOP
    6 CONTINUE
*
      IF( ISMULT(1). NE. ISMULT(2)) THEN
         WRITE(IWRITE) 'Delta S .NE. 0 strictly forbidden'
         STOP
      END IF
      IF (ET(1) .LT. ET(2)) THEN
         LMULTU = LMULT(2)
      ELSE
         LMULTU = LMULT(1)
      END IF
*
*  ***  DETERMINE THE RYDBERG CONSTANT
*
      IF (III.EQ.2) GO TO 66
      WRITE(0,'(/A)') ' Default Rydberg constant (Y or N) ? '
      READ(IREAD,'(A1)') PP
      IF ( PP .EQ. 'Y' .OR. PP .EQ. 'y' ) THEN
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
         WRITE(0,'(A)') ' Enter the mass of the atom'
         READ(IREAD,*) ZMU
      END IF
      RYDBRG = 109737.31534/(1. + 548.579903D-6/ZMU)
   66 CONTINUE
      CALL INITA
      CALL INITR
*
*  *****  SET UP DATA FOR AN OSCILLATOR STRENGTH CALCULATION
*
      IF (III .EQ. 1) THEN
      WRITE (6,3) ATOM,Z,LMULT,ISMULT
3     FORMAT(/3X,A6,3X,'Z = ',F3.0,3X,'Angular multiplicities =',2I3
     :       /22X,                    'Spin    multiplicities =',2I3//
     :       3X,'Orbitals')
      WRITE(IWRITE,21) NW(1),(EL(I),I=1,NW(1))
21    FORMAT(' Initial: ',I3,';',30(1X,A3))
      WRITE(IWRITE,22) NW(2)-NW(1),(EL(I),I=NW(1)+1,NW(2))
22    FORMAT(' Final  : ',I3,';',30(1X,A3))
      WRITE(IWRITE,'(/)')
      RHO=-4.0
      DO 4 J=1,NO
      R(J)=DEXP(RHO)/Z
      RR(J)=R(J)*R(J)
      R2(J)=DSQRT(R(J))
 4    RHO=RHO+H
      DO 9 I = 1,10
      DO 9 J = 1,10
      RLINT(I,J) = D0
      RVINT(I,J) = D0
      OVINT(I,J) = D0
9     CONTINUE
      END IF
*
*  *****  PERFORM THE CALCULATION FOR THE MULTIPLET
*
      DE = ET(2) - ET(1)
      IF (PRINT) WRITE (IWRITE,20)
20    FORMAT(//12X,'Transition Integrals',23X,'Length',5X,'Velocity'/)
      SL = D0
      SV = D0
   30 READ(9,31) C,END,KVAL,ELC(1),J1,ELC(2),J2,(ELC(I),ELC(I+1),
     :                       IP(I/2),I=3,7,2)
   31 FORMAT(F14.8,A1,1X,I1,1X,2(A3,I3,1X),3(1X,A3,1X,A3,1X,I2))
      IF ( END .EQ. ASTER ) GO TO 61
      KV = KVAL
      DO 32 I = 1,8
         IPT(I) = 0
         IF (ELC(I) .EQ. '   ') GO TO 32
         IF (MOD(I,2) .EQ. 1) THEN
            IB = 1
            IL = NW(1)
         ELSE
            IB = NW(1) + 1
            IL = NW(2)
         ENDIF
         DO 33 J = IB,IL
            IF (ELC(I) .EQ. EL(J)) THEN
               IPT(I) = J
               GO TO 32
            ENDIF
   33    CONTINUE
         WRITE(IWRITE,34) ELC(I),I,ELC
   34    FORMAT(3X,A3,' ELECTRON NOT FOUND - #',I3,'  IN THE LIST ',
     :   8(1X,A3))
         STOP
   32 CONTINUE
      J2 = J2 + NCF(1)
      L1 = L(IPT(1))
      L2 = L(IPT(2))
      LM = MAX0(L1,L2)
      IF (KVAL. EQ. 1) THEN
         K = LM - L1
         D = (-1)**K*DSQRT(LM*D1)
      ELSE
         IF (L1 .EQ. L2) THEN
            XX = D1*L1*(L1+D1)*(D2*L1+D1)
            XX = XX/((D2*L1+D3)*(D2*L1-D1))
            D = -1.D0 * DSQRT(XX)
         ELSE
            XX = D3*LM*(LM-D1)
            XX = XX/(D2*(D2*LM-D1))
            D = DSQRT(XX)
         END IF
      END IF
      D = D*C*WT(J1)*WT(J2)
      DO 50 I = 1,3
      IF (IP(I) .EQ. 0) GO TO 50
      II = 2*I +1
      I1 = IPT(II)
      I2 = IPT(II+1)-NW(1)
      IF (I1 .GT. 10 .OR. I2 .GT. 10) THEN
         OV =QUADR(IPT(II),IPT(II+1),0)
      ELSE
         IF (OVINT(I1,I2).EQ.0) OVINT(I1,I2)=
     :                  QUADR(IPT(II),IPT(II+1),0)
         OV = OVINT(I1,I2)
      ENDIF
      D = D*OV**IP(I)
   50 CONTINUE
      I1 = IPT(1)
      I2 = IPT(2) - NW(1)
      IF (I1 .GT. 10 .OR. I2 .GT. 10) THEN
         QR = QUADR(IPT(1),IPT(2),KVAL)
         IF (KVAL .EQ. 1) THEN
            QV = GRAD(IPT(1),IPT(2))
         ELSE
            QV = GRAD2(IPT(1),IPT(2))
         END IF
      ELSE
         IF (RLINT(I1,I2) .EQ. D0) RLINT(I1,I2)=
     :             QUADR(IPT(1),IPT(2),KVAL)
         IF (RVINT(I1,I2) .EQ. D0) THEN
            IF (KVAL .EQ. 1) THEN
               RVINT(I1,I2)=GRAD(IPT(1),IPT(2))
            ELSE
               RVINT(I1,I2)=GRAD2(IPT(1),IPT(2))
            END IF
         END IF
         QR = RLINT(I1,I2)
         QV = RVINT(I1,I2)
      ENDIF
*
42    TL = D*QR
      SL = SL + TL
      TV = D*QV
      SV = SV + TV
      IF ((.NOT. PRINT) .OR. DABS(TL) .LE. TOL) GO TO 76
      WRITE(IWRITE,35) J1, J2-NCF(1),
     :     ELC(1),ELC(2),(ELC(I),ELC(I+1),IP(I/2),I=3,5,2),TL,TV/DE
35    FORMAT(2X,I3,'->',I3,':',2X,A3,'->',A3,2X,
     :   2('<',A3,'|',A3,'>^',I2,2X),2F13.7)
76    CONTINUE
      GO TO 30
61    WRITE(IWRITE,36) SL,SV/DE
36    FORMAT(/10X,'SUM OF LENGTH FORM  ',16X,'= ',F16.8/
     :        10X,'SUM OF VELOCITY FORM',16X,'= ',F16.8)
      D = ET(2) - ET(1)
      KVAL = KV
      IF (KVAL .EQ. 1) THEN
         CL = D2*D*ISMULT(1)/D3
         CV = D2*ISMULT(1)/(D*D3)
         AFAC = 1.D0
      ELSE
         CL = D1*D**3*ISMULT(1)/D30
         CV = D2*D*ISMULT(1)/15.D0
         AFAC = 5.325132272D-05
      END IF
      GL = CL*SL**2*AFAC
      GV = CV*SV**2*AFAC
      WRITE (6,37) KVAL,GL,GV
37    FORMAT(7X,'E',I1,' FINAL gf-VALUES: LENGTH  = ',
     :     D14.6/,27X,'VELOCITY= ', D14.6)
*
*  *****  COMPUTE TRANSITION ENERGY  IN ANGSTROMS
*
82    DD = DABS(D)*D2*RYDBRG
      ANGS = D10**8/DD
      WRITE(IWRITE,38) DD,ANGS,D,RYDBRG
38    FORMAT(/10X,'ENERGY DIFFERENCE OF THE STATES:',6X,1PD16.8,' CM-1',
     :   /48X,1PD16.8,' ANGSTROMS(vac)'/48X,1PD16.8,' A. U.'//
     :   10X,'Rydberg constant for conversion:',0PF20.2)
      REWIND 9
      AKIL = 6.6702*1.D7*GL/(ANGS**2*LMULTU*ISMULT(1))
      IF (AKIL .LT. 0) AKIL = -AKIL
      WRITE (6,39) AKIL
39    FORMAT(/10X,'TRANSITION PROBABILITY: (final->initial) = ',
     :        1PD16.7,' 10^8 sec-1')
1000  CONTINUE
   99 CLOSE(UNIT=1)
      CLOSE(UNIT=2)
      CLOSE(UNIT=7)
      CLOSE(UNIT=8)
      CLOSE(UNIT=9)
      END
*     ------------------------------------------------------------------
* *** GRAD2
*     ------------------------------------------------------------------
*
* *** THE GRAD2 FUNCTION SUBPROGRAM COMPUTES THE FOLLOWING DIRECTLY
* ***        <P(I)[R.D + F(DELTA)[ P(J)>
*
      DOUBLE PRECISION FUNCTION GRAD2(I,J)
      PARAMETER (NOD=220,NWD=30,NWD2=2*NWD)
      PARAMETER (IWRITE=6)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD2),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWD2),L(NWD2),MAX(NWD2),N(NWD2)
      DIMENSION Q(NOD)
      JJ=I
      II=J
      DO1K=1,NO
1     Q(K)=P(K,JJ)*R(K)
      LI=L(I)
      LJ=L(J)
      IL = IABS(LI - LJ)
      IF (IL .NE. 0 .AND. IL .NE. 2) GO TO 100
      A1=(LI+LJ+D2)/((LI+D1)*(LJ+D1))
      A2=((LJ+D5)*(LJ+D1)+(LJ+D1+D5)*(LI+D1))/((LI+D1)*(LJ+D1))
      A=A1-A2*(LI+LJ+D3)/((LI+LJ+D4)*(LJ+D5))
      FACT=(LJ+D5)/(LI+LJ+D3)
      G=R(1)**2*P(1,I)*P(1,J)*FACT*(1.+A*Z*R(1))
      MM=MIN0(MAX(I)+1,MAX(J)+1,ND)
      K=2
      F1=D5*(P(K+1,II)-P(K-1,II))
      F2=P(K+1,II)-D2*P(K,II)+P(K-1,II)
      G0=Q(K)*R(K)
      G1=D5*(Q(K+1)*R(K+1)-Q(K-1)*R(K-1))
      G2=Q(K+1)*R(K+1)-D2*Q(K)*R(K)+Q(K-1)*R(K-1)
      G=G+D2*F1*G0+(D2*F2*G1+F1*G2)/D3
      DO2K=4,MM,2
      F1=D5*(P(K+1,II)-P(K-1,II))
      F2=P(K+1,II)-D2*P(K,II)+P(K-1,II)
      F3=D5*(P(K+1,II)-P(K-2,II))-D2*F1
      F4=P(K+2,II)+P(K-2,II)-D4*(P(K+1,II)+P(K-1,II))
     1   +D6*P(K,II)
      G0=Q(K)*R(K)
      G1=D5*(Q(K+1)*R(K+1)-Q(K-1)*R(K-1))
      G2=Q(K+1)*R(K+1)-D2*Q(K)*R(K)+Q(K-1)*R(K-1)
      G3=D5*(Q(K+2)*R(K+2)-Q(K-2)*R(K-2))-D2*G1
      G4=Q(K+2)*R(K+2)+Q(K-2)*R(K-2)-D4*(Q(K+1)*R(K+1)
     1   +Q(K-1)*R(K-1))+D6*Q(K)*R(K)
      G=G+D2*F1*G0+(D2*F2*G1+F1*G2)/D3
     1  -(F1*G4-F4*G1+D4*(F2*G3-F3*G2))/90.E0
2     CONTINUE
      U=QUADR(JJ,II,0)
      G=G-D5*U
      DELTA=LJ-LI
      IF(DELTA)10,11,12
10    GRAD2=G-(LI-D2)*U
      RETURN
11    GRAD2=G+(D1+D5)*U
      RETURN
12    GRAD2=G+(LI+D3)*U
      RETURN
100   WRITE(6,101) I,J
101   FORMAT(5X,'L(I)-L(J) NOT=0,2 FOR I = ',I2,' AND J = ',I2)
      STOP
      END

