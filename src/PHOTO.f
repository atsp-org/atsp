*     ------------------------------------------------------------------ 
*	PHOTO -- A GENERAL PROGRAM FOR PHOTOIONIZATION CROSS-SECTIONS
*
*                C O P Y R I G H T -- 1994
* 
*	by C. Froese Fischer 
*	   Vanderbilt University 
*	   Nashville, TN 37235 USA 
* 
*     Based on the LSTR program for bound transitions
*     Computer Physics Communications, Vol.  9, 381 (1975).
*                                      Vol. 14, 275 (1978).
*     ------------------------------------------------------------------
* 
      PROGRAM PHOTO
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      PARAMETER(NOD=220,NWD=30,NWD2=2*NWD,NCD=500,NCD2=2*NCD) 
      INTEGER Q(5)
      CHARACTER*1 END,ASTER,PP 
      CHARACTER*3 EL(NWD2),ELC(8),COUPLE(9),ELI 
      CHARACTER*6 ATOM,TERM,ALABEL(2),AT,TT 
      CHARACTER*13 LABEL(2)
      CHARACTER*24 CFILE(2),WFILE(2),DFILE,pfile*48
      CHARACTER*66 CONFIG(NCD2) 
      CHARACTER*72 HEADER 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID 
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD2),YK(NOD), 
     :   YR(NOD),X(NOD),AZ(NWD2),L(NWD2),MAX(NWD2),N(NWD2) 
      COMMON /STATE/WT(NCD2),ETarg(2),Etotal(2),JV(2),NCF(2),NJ(2) 
      COMMON /MULT/OVINT(10,10),RLINT(10,10),RVINT(10,10) 
      DIMENSION IP(3),IPT(8), LMULT(2), NOCCSH(NCD2), NW(2), convert(2)
      LOGICAL PRINT 
      PARAMETER (IREAD=5,IWRITE=6) 
      DATA ASTER/'*'/, LABEL/'Initial State','Final   State'/ 
* 
*  *****  DETERMINE DATA 
*
999   write(iwrite,998)
998   format(//20x,'================================'/
     :         20X,'  PHOTIONIZATION CROSS-SECTION '/
     :         20X,'================================'/)
*
1     write(0,'(A)') ' Name of Initial State'
      read(iread,'(A)') cfile(1)
      write(0,'(A)') ' Name of Final State'
      read(iread,'(A)') cfile(2)
      dfile = 'mltpol.lst'
      do 2 i = 1,2
        j = index(cfile(i),' ')
        if (j .eq. 1) then
          write(0,*) ' Names may not start with blanks'
      	  go to 1
        end if
        cfile(i) = cfile(i)(1:j-1)//'.c'
        wfile(i) = cfile(i)(1:j-1)//'.w'
	if (i .eq. 1) then
	  pfile = cfile(i)(1:j-1)
	else
	  jj = index(pfile,' ')
	  pfile(jj:) = '.'//cfile(i)(1:j-1)//'.ph'
	end if
2     continue
*
      call transunit(convert)

      OPEN(UNIT=1,FILE=WFILE(1),STATUS='OLD',FORM='UNFORMATTED')
      OPEN(UNIT=2,FILE=WFILE(2),STATUS='OLD',FORM='UNFORMATTED')
      OPEN(UNIT=7,FILE=CFILE(1),STATUS='OLD') 
      OPEN(UNIT=8,FILE=CFILE(2),STATUS='OLD') 
      OPEN(UNIT=9,FILE=DFILE,STATUS='OLD') 
      OPEN(UNIT=10,FILE=pfile,STATUS='UNKNOWN')
*
      WRITE(0,'(A)') ' INTERMEDIATE PRINTING (Y OR N) ?  ' 
      READ(IREAD,'(A1)') PP 
      IF ( PP .EQ. 'Y' .OR. PP .EQ. 'y') THEN 
         PRINT = .TRUE. 
         WRITE(0,'(A)') '  TOLERANCE FOR PRINTING ?  ' 
         READ(IREAD,'(F6.4)') TOL 
      ELSE 
         PRINT = .FALSE. 
      ENDIF 
*
*     The conversion factor is 4(Pi*a0)^2/c/3 where
*        a0 (Bohr radius) = 5.29177259 10^(-9) cm
*        c (fine-structure constant) = 137.0359895
      Factor = 2.689095
*
*     .. Write header for photo.out file
      i = index(cfile(1),'.')
      j = index(cfile(2),'.')
      WRITE(10,*) 'Files: ', cfile(1)(1:i-1),'.', cfile(2)(1:j-1)
      WRITE(10,70)
70    FORMAT(/2X,'k**2',5X,'Phen',7X,'CRL',7X,'CRV',7x,'WaveL',
     :       7X,'SL',8X,'SV',6x,'DPM'/)
      II = 1
      NCASE = 1
60    READ(9,'(A)') HEADER 
      I = 1 
      M = 1 
      IF (NCASE .GT. 1) THEN
      	 II = 2
cxi	 I = NW(1) +1
      	 M = NCF(1) + 1
	 NWF1 = NW(1)
	 NWF2 = NW(2) - NW(1)
      END IF

*     Set value of NWF in COMMON to be that for the second set.
      NWF = NWF2
* 
*  *****  READ THE RADIAL FUNCTIONS DURING THE FIRST CYCLE 
* 
12    if( Ncase .eq. 1 ) then 
	READ(II,END=13) AT,TT,EL(I),MR,Z,ETI,EKI,AZ(I), 
     :  (P(JJ,I),JJ=1,MR) 
      IF (I.GT.(NWD2)) STOP ' TOO MANY ELECTRONS: MAX=NWD2' 
      IF (MR .EQ. 0) GO TO 13
      IF (EL(I)(1:1) .NE. ' ') THEN 
         N(I) = ICHAR(EL(I)(1:1)) - ICHAR('1') 
         L(I) = LVAL(EL(I)(2:2)) 
      ELSE 
         N(I) = ICHAR(EL(I)(2:2)) - ICHAR('1') 
         L(I) = LVAL(EL(I)(3:3)) 
      ENDIF 

      else
	READ(II,END=13) AT,TT,ELI,MR,Z,ETI,EKI,AZI, 
     :  (X(JJ),JJ=1,MR) 
      IF (MR .EQ. 0) GO TO 13
	call eptr(EL(NWF1+1),ELI,I,*99)
	if( I .le. 0 ) goto 12
	I = I + NWF1
	AZ(I) = AZI
	do 71 jj=1,MR
71	P(JJ,I) = X(JJ)
      endif

      IF (II .gt. 1) then
	ekk = -eti
      END IF
      MM = MR+1 
      DO 24 J = MM,NO 
24    P(J,I) = D0 
 8    MAX(I)=MR 
      I = I+1 
      GO TO 12 
   13 if( Ncase .eq. 1 ) NW(II) = I-1 
* 
*  *****  READ THE CONFIGURATIONS 
* 
      IU = II+6 
      READ(IU,10,END=99) ATOM,ALABEL(II),ETarg(II),Etotal(II)
   10 FORMAT(3X,2A6,2F14.7/) 
      WRITE(IWRITE, 11) LABEL(II),ATOM,ALABEL(II),ETarg(II),Etotal(II)
   11 FORMAT(/3X,A/3X,13('-')/3X,2A6,2F16.7/) 
   15 READ(IU,14,END=18) (ELC(K),Q(K),K=1,5),WT(M) 
      print 14, (ELC(K),Q(K),K=1,5),WT(M) 
   14 FORMAT(5(1X,A3,1X,I2,1X),F12.5) 
      IF (ELC(1)(1:1) .EQ. '*') GO TO 18
      READ(IU,'(9(5X,A3))',END=18) (COUPLE(J),J=1,9) 
      print '(9(5X,A3))', (COUPLE(J),J=1,9) 
      NOCC = 0
16    IF (ELC(NOCC+1) .NE. '   ' ) THEN 
         NOCC = NOCC+1 
         IF (NOCC .LT. (5)) GO TO 16 
      END IF 
      MM = M 
      IF (II.EQ.2) MM = M - NCF(1) 
      NOCCSH(M) = NOCC 
      CALL PACK(NOCC,ELC,Q,COUPLE,CONFIG(M)) 
      K = 66 
   17 IF (CONFIG(M)(K:K) .EQ. ' ') THEN 
         K = K-1 
         GO TO 17 
      END IF 
      WRITE(IWRITE,'(I4,3X,F12.5,3X,A)') MM,WT(M),CONFIG(M)(1:K) 
      M = M+1 
      IF (M.GT.(NCD2)) THEN 
         WRITE(0,'(A/A,I5)') 
     :       ' TOO MANY CONFIGURATIONS IN THE INITIAL AND FINAL STATE', 
     :       ' MAXIMUM FOR COMBINED SUM = ',(NCD2)-1 
      END IF 
      GO TO 15 
*
   18 NCF(II) = M-1 
      DO 5 J = 6,1,-1 
         IF (ALABEL(II)(J:J) .NE. ' ') THEN 
             LMULT(II) = 2*LVAL(ALABEL(II)(J:J)) + 1 
             GO TO 6 
         ENDIF 
    5 CONTINUE 
      WRITE(IWRITE,'(A)') 'No term value found in ',ALABEL(II) 
      STOP 
    6 II = II + 1
      IF (II .LE. 2) GO TO 12
*
      IF (NCASE .EQ. 1) THEN
cxi
cxi	read the lines in file(2).c for orthogonalities
cxi
72	read(8,'(A3)') ELI
	if( ELI(1:1) .ne. '*' ) goto 72

        CALL INITA
        CALL INITR
        DO 4 J=1,NO 
        R(J)=DEXP(RHO)/Z 
        RR(J)=R(J)*R(J) 
        R2(J)=DSQRT(R(J)) 
  4     RHO=RHO+H 
*
        WRITE (6,3) ATOM,Z,LMULT 
  3     FORMAT(/3X,A6,3X,'Z = ',F3.0,3X,'Angular multiplicities =',2I3// 
     :       3X,'Orbitals') 
        WRITE(6,21) NW(1),(EL(I),I=1,NW(1)) 
  21    FORMAT(' Initial: ',I3,';',30(1X,A3)) 
        WRITE(6,22) NW(2)-NW(1),(EL(I),I=NW(1)+1,NW(2)) 
  22    FORMAT(' Final  : ',I3,';',30(1X,A3)) 
        WRITE(6,'(/)') 
      END IF
      PHEN = etotal(2) - etarg(1)
* 
*  *****  SET UP DATA FOR AN OSCILLATOR STRENGTH CALCULATION 
* 
      DO 9 I = 1,10 
      DO 9 J = 1,10 
      RLINT(I,J) = D0 
      RVINT(I,J) = D0 
      OVINT(I,J) = D0 
9     CONTINUE 
* 
*  *****  PERFORM THE CALCULATION FOR THE MULTIPLET 
* 
      IF (PRINT) WRITE (IWRITE,20) 
20    FORMAT(//12X,'Transition Integrals',24X,'Length',7X,'Velocity'/) 
      SL = D0 
      SV = D0 
   30 READ(9,31) C,END,ELC(1),J1,ELC(2),J2,(ELC(I),ELC(I+1),IP(I/2), 
     :                                       I=3,7,2) 
   31 FORMAT(F14.8,A1,3X,2(A3,I3,1X),3(1X,A3,1X,A3,1X,I2)) 
      IF ( END .EQ. ASTER ) GO TO 61 
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
      K = LM - L1 
      D = (-1)**K*DSQRT(LM*D1)*C*WT(J1)*WT(J2) 
      DO 50 I = 1,3 
      IF (IP(I) .EQ. 0) GO TO 50 
      II = 2*I +1 
      I1 = IPT(II) 
      I2 = IPT(II+1)-NW(1) 
      IF (I1 .GT. 10 .OR. I2 .GT. 10) THEN 
         OV =QUADR(IPT(II),IPT(II+1),0) 
      ELSE 
         IF (OVINT(I1,I2) .EQ. 0)OVINT(I1,I2)=QUADR(IPT(II),IPT(II+1),0) 
         OV = OVINT(I1,I2) 
      ENDIF 
      D = D*OV**IP(I) 
   50 CONTINUE 
      I1 = IPT(1) 
      I2 = IPT(2) - NW(1) 
      IF (I1 .GT. 10 .OR. I2 .GT. 10) THEN 
         QR = QUADR(IPT(1),IPT(2),1) 
         QV = GRAD(IPT(1),IPT(2)) 
      ELSE 
         IF (RLINT(I1,I2) .EQ. D0) RLINT(I1,I2)=QUADR(IPT(1),IPT(2),1) 
         IF (RVINT(I1,I2) .EQ. D0) RVINT(I1,I2)=GRAD(IPT(1),IPT(2)) 
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
     :     ELC(1),ELC(2),(ELC(I),ELC(I+1),IP(I/2),I=3,5,2),TL,TV/PHEN 
35    FORMAT(2X,I3,'->',I3,':',2X,A3,'->',A3,2X, 
     :   2('<',A3,'|',A3,'>^',I2,2X),2F13.7) 
76    CONTINUE 
      GO TO 30 
61    CONTINUE
      WRITE(IWRITE,36) SL, SV/PHEN
36    FORMAT(/10X,'SUM OF LENGTH FORM  ',16X,'= ',F16.8/ 
     :        10X,'SUM OF VELOCITY FORM',16X,'= ',F16.8) 
      CRL = (Factor * PHEN * SL**2 )/(LMULT(1))
      CRV = (Factor * SV**2 )/(PHEN*LMULT(1))
      WRITE (6,37)  CRL,CRV 
37    FORMAT(10X,'PHOTOIONIZATION CROSS-SECTION: LENGTH  = ',F13.8/ 
     :    38X,'VELOCITY= ', F16.8) 
* 
*  *****  COMPUTE TRANSITION ENERGY  IN   eV
* 
82    PHENEV = PHEN*convert(2)
      WL = 1.D8/(PHEN*convert(1))
      WRITE(IWRITE,38) PHEN, PHENEV
38    FORMAT(/28X,'PHOTON ENERGY:',6X,1PD16.8,' au', 
     :   /48X,1PD16.8,' eV'// )
      WRITE(10,103) EKK,PHEN,CRL,CRV,WL,SL,SV/PHEN
103   FORMAT(F7.4,3F10.5,F10.1,2F10.4,F9.4)
      REWIND 9 
      NCASE = NCASE + 1
      GO TO 60
99    continue
      END 

*     ------------------------------------------------------------------
*       T R A N S U N I T
*     ------------------------------------------------------------------
	subroutine transunit(convert)
	implicit real*8 (a-h,o-z)
	dimension convert(*)
*
	write(0,*) ' convert the energy unit for atomic weight A '
	write(0,*) ' input A '
	read(*,*) AWT
*        
*        Electron Mass
	 
	 EM = 548.579909070D-6
	 ratio = AWT/(AWT + EM)
*
*	Hartree --> cm^1
*
	 convert(1) = 2*109737.31534*ratio
*
*	Hartree --> eV 
*
	convert(2) = 27.2113961*ratio

100	format(1x, ' Atomic Weight = ', I3 )
101	format(1x, ' 1  Hartree    = ', F15.7, '  cm^-1 ' )
102	format(1x, '               = ', F15.7, '  eV    ' )
	
	write(*,100) AWT
	write(*,101) convert(1)
	write(*,102) convert(2)
	end

