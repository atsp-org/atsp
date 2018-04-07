*---------------------------------------------------------------- 
*
	PROGRAM LSJTR
*
*                 C O P Y R I G H T -- 1994
*
*  M.R.GODEFROID LABORATOIRE DE CHIMIE PHYSIQUE MOLECULAIRE 
*                UNIVERSITE LIBRE OF BRUSSELS, BELGIUM 
*  DECEMBER 1981, VANDERBILT UNIVERSITY, NASHVILLE. 
*
*  Modified by C. Froese Fischer 
*              Vanderbilt University, Nashville TN 
*  April 1983 
* 
*  Computer Physics Communications, Vol. 64, 501--519 (1991).
* ---------------------------------------------------------------- 
*
*  THIS PROGRAM EVALUATES THE GF VALUES, TRANSITION PROBABILITIES 
*  AND TRANSITION ENERGIES OF ELECTRIC AND MAGNETIC PROCESSES IN 
*  THE LSJ INTERMEDIATE COUPLING SCHEME. 
* 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      PARAMETER(NOD=220,NWD=30,NWD2=2*NWD,NCD=500,NCD4=4*NCD,
     :          MD=20,MD2=2*MD,NTD=500) 
      INTEGER IQ(5),NW(2) 
      CHARACTER*1 END,ENDST,ASTER,PP,MS
      CHARACTER*3 EL(NWD2),ELC(8),COUPLE(9) 
      CHARACTER*6 ATOM,TT,LABEL(2)*13
      CHARACTER*72 HEADER,BIDON,ECLOSD 
      CHARACTER*24 CFILE(2),WFILE(2),JFILE(2),DFILE 
      CHARACTER*66 CONFIG(NCD4) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID 
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS 
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD2),YK(NOD), 
     :   YR(NOD),X(NOD),AZ(NWD2),L(NWD2),MAX(NWD2),N(NWD2) 
      COMMON /INFORM/IREAD,IWRITE,IOUT,ISC(8) 
      COMMON /DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL 
      COMMON /STATE/WT(NCD4,MD2),ET(2,MD2),JV(2),NCF(2),NJ(2),IS(NCD4), 
     : ILABEL(2,NCD4),LL(NCD4) 
      COMMON /MULT/TTL(NTD),JJ1(NTD),JJ2(NTD),II1(NTD),II2(NTD),
     :             RLINT(NWD,NWD,3),FLINE(NTD) 
      DIMENSION IPT(8),IP(3),LENGTH(NCD4) 
      LOGICAL PRINT,FIRST
      DATA ASTER/'*'/, LABEL/'Initial State','Final   State'/ 
      IREAD = 5
      IWRITE = 6
      IBUG1 = 0 
* 
* --- SET FACTORIALS 
* 
      CALL INITA
      CALL INITR
* 
*  *****  DETERMINE DATA 
* 
999   write(IWRITE,998)
998   format(//20x,'============================='/
     :         20X,'  LSJ - T R A N S I T I O N '/
     :         20X,'============================='/)
*
  1   CONTINUE

CSUN  i = iargc()
CSUN  if (i .lt. 2) then
         WRITE(0,*)  ' Name of Initial State'
         READ(IREAD,'(A)') cfile(1)
         WRITE(0,*)  ' Name of Final State'
       	 READ(IREAD,'(A)') cfile(2)
CSUN  else
CSUN   	 call getarg(1,cfile(1))
CSUN   	 call getarg(2,cfile(2))
CSUN  end if
CSUN  if (i .eq. 3) then
CSUN   	 call getarg(3,dfile)
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
	 jfile(i) = cfile(i)(1:j-1)//'.j'
2     continue
*
      WRITE(0,*) ' INTERMEDIATE PRINTING (Y OR N) ?  ' 
      READ(IREAD,'(A1)') PP 
      IF ( PP .EQ. 'Y' .OR. PP .EQ. 'y' ) THEN 
         PRINT = .TRUE. 
         WRITE(0,*) '  TOLERANCE FOR PRINTING ?  ' 
         READ(IREAD,*) TOL 
      ELSE 
         PRINT = .FALSE. 
      ENDIF 
*
*     Not all compilers allow the "READONLY" option which is
*     sometimes needed when the initial and final states are
*     the same.  In MicroSoft FORTRAN, "READONLY" should be
*     replaced by "MODE='READ'".
*     OPEN(UNIT=1,FILE=WFILE(1),STATUS='OLD',FORM='UNFORMATTED',
*    :     READONLY)
*     OPEN(UNIT=2,FILE=WFILE(2),STATUS='OLD',FORM='UNFORMATTED', 
*    :     READONLY)
*     OPEN(UNIT=3,FILE=CFILE(1),STATUS='OLD',READONLY) 
*     OPEN(UNIT=4,FILE=CFILE(2),STATUS='OLD',READONLY) 
*     OPEN(UNIT=7,FILE=JFILE(1),STATUS='OLD',READONLY) 
*     OPEN(UNIT=8,FILE=JFILE(2),STATUS='OLD',READONLY) 
*     OPEN(UNIT=9,FILE=DFILE,STATUS='OLD') 
*     OPEN(UNIT=10,FILE='tr.lsj',STATUS='UNKNOWN',ACCESS='APPEND') 
CSUN  OPEN(UNIT=10,FILE='tr.lsj',STATUS='UNKNOWN',fileopt='eof') 
*
*  *****  The fileopt='eof' allows for information to be
*         appended to an existing file on SUN UNIX.
*     On other systems, it may simply not be possible to have
*     two units attached to one file. The default here.
      OPEN(UNIT=1,FILE=WFILE(1),STATUS='OLD',FORM='UNFORMATTED')
      OPEN(UNIT=2,FILE=WFILE(2),STATUS='OLD',FORM='UNFORMATTED')
      OPEN(UNIT=3,FILE=CFILE(1),STATUS='OLD') 
      OPEN(UNIT=4,FILE=CFILE(2),STATUS='OLD') 
      OPEN(UNIT=7,FILE=JFILE(1),STATUS='OLD') 
      OPEN(UNIT=8,FILE=JFILE(2),STATUS='OLD') 
      OPEN(UNIT=9,FILE=DFILE,STATUS='OLD') 
      OPEN(UNIT=10,FILE='tr.lsj',STATUS='UNKNOWN') 

* 
*  *****  READ FOR THE INITIAL (II=1) AND FINAL (II=2) STATES, 
* 
      IC = 1 
      I = 1 
      DO 6 II=1,2 
* 
*  *****  THE RADIAL FUNCTIONS ...(UNITS 1 AND 2) 
* 
   12 READ(II,END=13) ATOM,TT,EL(I),MR,Z,ETI,EKI,AZ(I), 
     :                (P(J,I),J=1,MR) 
      IF (EL(I)(1:1) .NE. ' ') THEN 
         N(I) = ICHAR(EL(I)(1:1)) - ICHAR('1') + 1 
         L(I) = LVAL(EL(I)(2:2)) 
      ELSE 
         N(I) = ICHAR(EL(I)(2:2)) - ICHAR('1') + 1 
         L(I) = LVAL(EL(I)(3:3)) 
      ENDIF 
      MM = MR + 1 
      DO 24 J = MM,NO 
   24 P(J,I) = D0 
      MAX(I) = MR 
      IF(IBUG1.NE.0) WRITE(IWRITE,512) II,I,EL(I),N(I),L(I) 
  512 FORMAT(/,1H ,' STATE ',I2,' I = ',I3,' ELECTRON ',A3,' N = ',I3, 
     :             ' L = ',I3) 
      I = I+1 
      GO TO 12 
   13 NW(II) = I-1 
* 
*  *****  AND THE LIST OF CONFIGURATIONS ...(UNITS 3 AND 4). 
* 
      IU = II+2 
      WRITE(IWRITE,'(/1X,A/1X,A)') LABEL(II),'-------------' 
      READ(IU,'(A/A)') HEADER,ECLOSD 
   15 READ(IU,14,END=56) (ELC(K),IQ(K),K=1,5) 
   14 FORMAT(5(1X,A3,1X,I2,1X)) 
      READ(IU,'(9(5X,A3))',END=56) (COUPLE(J),J=1,9) 
      NOCC = 0 
16    IF (ELC(NOCC+1) .NE. '   ' ) THEN 
         NOCC = NOCC+1 
         GO TO 16 
      END IF 
      IF (NOCC .EQ. 0) GO TO 56 
      NEL = IQ(1)+IQ(2)+IQ(3)+IQ(4)+IQ(5)
      IF (IC .GT. (NCD4)) STOP ' TOO MANY CONFIGURATIONS: MAX=(800)' 
      CALL PACK(NOCC,ELC,IQ,COUPLE,CONFIG(IC)) 
      K = 66 
   17 IF (CONFIG(IC)(K:K) .EQ. ' ') THEN 
         K = K-1 
         GO TO 17 
      END IF 
      LENGTH(IC) = K 
      MM = IC 
      IF (II .EQ. 2) MM = IC - NCF(1) 
      WRITE(IWRITE,'(I4,3X,A)') MM,CONFIG(IC)(1:K) 
      LAST = 2*NOCC - 1 
      LL(IC) = 2*LVAL(COUPLE(LAST)(2:2)) 
      IS(IC) = ICHAR(COUPLE(LAST)(1:1)) - ICHAR('1') 
      IC = IC+1 
      GO TO 15 
   56 IF (II .EQ. 1) NCF(II) = IC - 1 
      IF (II .EQ. 2) NCF(II) = IC - 1 - NCF(1) 
      IF(IBUG1.NE.0) WRITE(IWRITE,514) II,NW(II),NCF(II) 
  514 FORMAT(/,1H ,' STATE ',I2,' NWF = ',I3,' NCF = ',I3) 
    6 CONTINUE 
      J = 4
    7 IF (ECLOSD(J:J) .NE. ' ') THEN
	 NEL = NEL + 4*LVAL(ECLOSD(J:J)) + 2
	 J = J + 4
	 IF (J .LE. 72) GO TO 7
      END IF
      WRITE(IWRITE,3) ATOM,Z,NEL 
      WRITE(10,3) ATOM,Z,NEL 
    3 FORMAT(/10X,A6,3X,'Z = ',F3.0,3X,'N = ',I3) 
* 
* --- Determine suitable Rydberg constant 
* 
      WRITE(0,*) ' Default Rydberg constant (Y or N) ? ' 
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
         WRITE(0,*) ' Enter the mass of the atom' 
         READ(IREAD,*)  ZMU 
      END IF 
      RYDBRG = 109737.31568508/(1. + 548.579909070D-6/ZMU) 
* 
*  *****  SET UP DATA FOR AN OSCILLATOR STRENGTH CALCULATION 
* 
      RHO=-4.0 
      DO 4 J=1,NO 
      R(J)=DEXP(RHO)/Z 
      RR(J)=R(J)*R(J) 
      R2(J)=DSQRT(R(J)) 
    4 RHO=RHO+H 
      DO 9 I=1,30 
      DO 9 J=1,30 
      DO 9 K=1,3 
      RLINT(I,J,K) = D0 
    9 CONTINUE 
      READ(9,'(A)') HEADER 
      WRITE(IWRITE,'(/A)') HEADER 
* 
*  *****  SKIP THE HEADER CARD ON UNITS 7 AND 8 
* 
   71 DO 72 II = 1,2 
      IU = II + 6 
      READ(IU,'(A)') BIDON 
   72 CONTINUE 
* 
*  *****  PERFORM THE CALCULATION FOR THE MULTIPLET 
* 
      NMA = 0 
      NMB = 0 
      NTR = 0 
      FIRST = .TRUE. 
      K = 1 
   30 READ(9,31,END=99) C,END,MS,LAM,ELC(1),J1,ELC(2),J2, 
     :(ELC(I),ELC(I+1),IP(I/2),I=3,7,2) 
   31 FORMAT(F14.8,2A1,I1,1X,2(A3,I3,1X),3(1X,A3,1X,A3,1X,I2)) 
      IF( END .EQ. ASTER ) GO TO 61 
      IF (K .GT. (7500)) STOP ' TOO MANY TERMS IN THE SUM: MAX=(7500)' 
      IF (END .EQ. 'e' ) THEN 
         END = 'E' 
        ELSE IF (END .EQ. 'm') THEN 
         END = 'M' 
      END IF 
      IF(K .EQ. 1) THEN 
         ENDST = END 
         LAMST = LAM 
         IF( END .EQ. 'E' ) THEN 
            WRITE(IWRITE,401) LAM 
         ELSE 
            WRITE(IWRITE,402) LAM 
         ENDIF 
      ENDIF 
  401 FORMAT(//,' ELECTRIC TRANSITION OF ORDER ',I3,/) 
  402 FORMAT(//,' MAGNETIC TRANSITION OF ORDER ',I3,/) 
* 
*  ***** FOR EACH TRANSITION INTEGRAL (CURRENT INDICE K), SEARCH THE 
*  ***** ELECTRONS WITH MATCHING LABELS 
* 
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
         WRITE(IWRITE,34) ELC(I) 
   34    FORMAT(3X,A3,' ELECTRON NOT FOUND') 
         STOP 
   32 CONTINUE 
* 
*  *****  PICK UP THE L VALUES OF THE INTERACTING SHELLS AND 
*  *****  STORE THE CONFIGURATION INDICES 
* 
      J2 = J2 + NCF(1) 
      L1 = L(IPT(1)) 
      L2 = L(IPT(2)) 
      II1(K) = IPT(1) 
      II2(K) = IPT(2) 
      JJ1(K) = J1 
      JJ2(K) = J2 
      IF(IBUG1.NE.0) WRITE(IWRITE,516) K,EL(II1(K)),L1,JJ1(K), 
     :                                   EL(II2(K)),L2,JJ2(K) 
  516 FORMAT(/,1H ,' K = ',I3,/,10X,A3,' L = ',I3,' CFG1 = ',I3,' -> ', 
     :                              A3,' L = ',I3,' CFG2 = ',I3) 
* 
*  *****  CALCULATE THE ADEQUATE REDUCED MATRIX ELEMENT 
*  *****  ORIGINAL USE OF FANO AND RACAH PHASE CONVENTION 
*  *****  BUT IFAC2 CORRECTS TO GET THE CONDON AND SHORTLEY DEFINITION 
* 
      L12 = L2 - L1 
      IF (END .EQ. 'E') THEN 
         D = RME(L1,L2,LAM) 
         IFAC2 = IABS((L12 + LAM)/2) 
         D = D*(-1.D0)**IFAC2 
      ELSE 
         D = RME(L1,L2,LAM-1) 
         IFAC2 = IABS((L12 + LAM - 1)/2) 
         D = D*(-1.D0)**IFAC2 
         IF (MS .EQ. 'A') THEN 
            LL1 = L1 + L1 
            LL2 = L2 + L2 
            L3 = LAM + LAM 
            L4 = L3 - 2 
            L5 = 2 
            CALL GRACAH(LL2,L5,LL1,L4,LL2,L3,W) 
            FACTR = DFLOAT(L2*(L2+1)*(LL2+1)*(L3+1)) 
            D = D*DSQRT(FACTR)*W 
         ELSE 
            D = D*DSQRT(D3*D5) 
         ENDIF 
      ENDIF 
* 
*  ***** MULTIPLY BY THE VSHELL COEFFICIENT 
* 
      IF(IBUG1.NE.0) WRITE(IWRITE,517) C,D 
  517 FORMAT(10X,' VSHELL = ',F14.8,'  RME = ',F14.8) 
      D = D*C 
* 
*  ***** MULTIPLY BY MAXIMUM 3 OVERLAPS AND STORE THEM IN RLINT(I,J,1) 
*  ***** FOR I AND J LE 10 
* 
      DO 50 I =1,3 
      IF (IP(I) .EQ. 0) GO TO 50 
      II = 2*I + 1 
      I1 = IPT(II) 
      I2 = IPT(II+1) - NW(1) 
      IF (I1 .GT. 30 .OR. I2 .GT. 30) THEN 
         OV = QUADR(IPT(II),IPT(II+1),0) 
      ELSE 
         IF (RLINT(I1,I2,1) .EQ. 0) RLINT(I1,I2,1)=QUADR(IPT(II), 
     :                                             IPT(II+1),0) 
         OV = RLINT(I1,I2,1) 
      ENDIF 
      D = D*OV**IP(I) 
      IF(IBUG1.NE.0) WRITE(IWRITE,518) EL(IPT(II)),EL(IPT(II+1)),IP(I) 
     :                                 ,OV,D 
  518 FORMAT(10X,' OVRLP <',A3,'|R0|',A3,'>**',I3,' = ',F14.8,'  D = ', 
     :       D16.8) 
   50 CONTINUE 
* 
*  *****  MULTIPLY BY THE RADIAL INTEGRALS AND STORE THEM ONLY FOR 
*  *****  (M1), (E1,M2) AND (E2,M3) PROCESSES IN RLINT(I,J,K) : 
*  *****               R0   R1   R2 
*  *****          K     1    2    3 
*  *****          E     -    1    2 
*  *****          M     1    2    3 
* 
      I1 = IPT(1) 
      I2 = IPT(2) - NW(1) 
      LAM1 = LAM 
      IF (END .EQ. 'M') LAM1 = LAM - 1 
      LAM2 = LAM1 + 1 
      IF (I1 .GT. 10 .OR. I2 .GT. 10 .OR. LAM1 .GT. 2) THEN 
         QR = QUADR(IPT(1),IPT(2),LAM1) 
      ELSE 
         IF (RLINT(I1,I2,LAM2) .EQ. D0) RLINT(I1,I2,LAM2)=QUADR(IPT(1), 
     :                                                    IPT(2),LAM1) 
         QR = RLINT(I1,I2,LAM2) 
      ENDIF 
* 
*  *****  FOR EACH TRANSITION INTEGRAL, STORE THE RESULT IN TTL(K) 
* 
      TTL(K) = D*QR 
      IF(IBUG1.NE.0) WRITE(IWRITE,519) EL(IPT(1)),LAM1,EL(IPT(2)),QR, 
     :                                 TTL(K) 
  519 FORMAT(10X,' TRANS <',A3,'|R',I1,'|',A3,'> = ',F14.8,'  TTL = ', 
     :           D16.8) 
      K = K + 1 
      IF ( K .LE. NTD) THEN
	  GO TO 30 
      ELSE
	 WRITE(0,*) ' NTD dimension too small'
	 STOP
      END IF
* 
*  *****  END OF THE LIST OF TRANSITION INTEGRALS 
* 
   61 IF (.NOT. FIRST) THEN 
         NTR = K-1 
         NMB = NTR - NMA 
      ELSE 
         NMA = K-1 
         FIRST = .FALSE. 
         GO TO 30 
      ENDIF 
      IF ( NTR .EQ. 0 ) THEN 
         FIRST = .TRUE. 
         K = 1 
         GO TO 30 
      ENDIF 
      IF (IBUG1 .NE. 0) WRITE(IWRITE,501) NMA,NMB,NTR 
  501 FORMAT(/,1H ,' NMA =',I3,2X,' NMB =',I3,2X,' NTR =',I3,/) 
* 
*  *****  READ EIGENVECTORS FOR A PAIR OF J-VALUES 
* 
   67 IC = 1 
      II = 1 
   65 IU = II + 6 
      READ(IU,60,END=300) JV(II),NJ(II) 
   60 FORMAT(//8X,I4,10X,I4) 
      JC = IC + NCF(II) - 1 
      MFOUND = NJ(II) 
      DO 64 M = 1,MFOUND 
      READ(IU,63) ILABEL(II,M),ET(II,M),(WT(I,M),I=IC,JC) 
   63 FORMAT(/I6,F16.8/(7F10.6)) 
   64 CONTINUE 
      IC = JC + 1 
      IF (II .EQ. 1) ICST = IC 
      II = II + 1 
      IF (II .LE. 2) GO TO 65 
* 
*  ***** TEST ON TRIANGULAR CONDITION (J,J',K) 
* 
      JITOT = JV(1) 
      JFTOT = JV(2) 
      L3 = LAMST+LAMST 
      ATST = TRITST(JITOT,JFTOT,L3) 
      WRITE(IWRITE,507) JITOT,JFTOT 
  507 FORMAT(/,1H ,' 2*JI =',I3,' 2*JF =',I3, 
     :/,' --------------------') 
      IF (ATST .EQ. 1.D0) THEN 
         WRITE(IWRITE,560) 
  560 FORMAT(/,' TRIANGULAR CONDITION NOT SATISFIED ',/) 
         GO TO 536 
      ENDIF 
* 
*  *****  A PAIR OF J-VALUES HAS BEEN FOUND, 
*  *****  COMPUTE LINE FACTORS FOR THE PAIR 
* 
*  *****  THE SPIN-ORBIT COUPLING ORDER L+S IS USED ORIGINALLY, 
*  *****  BUT IFL AND IFR CORRECT TO GET S+L ORDER 
* 
*  *****  FOR E AND MA INTEGRALS 
* 
      DO 140 I = 1,NMA 
      J1 = JJ1(I) 
      J2 = JJ2(I) 
      LL1 = LL(J1) 
      LL2 = LL(J2) 
      IS1 = IS(J1) 
      IS2 = IS(J2) 
      CALL GRACAH(LL2,JFTOT,LL1,JITOT,IS1,L3,F) 
      IFAC = IABS((IS1+L3-LL2-JITOT)/2) + 2 
      PHAZ = (-1.D0)**IFAC 
      F = F*PHAZ 
      IF (ENDST .EQ. 'M') F = F/(LAMST + 1) 
      IFL = IABS((LL1+IS1-JITOT)/2) 
      IFR = IABS((LL2+IS2-JFTOT)/2) 
      F = F*(-1.D0)**(IFL+IFR) 
      IF (IBUG1.NE.0) WRITE(IWRITE,505) I,F 
  505 FORMAT(1H ,' I =',I3,' FLINE =',F10.6) 
  140 FLINE(I) = F 
* 
*  *****  FOR MB INTEGRALS 
* 
      NMA1 =  NMA + 1 
      DO 141 I = NMA1,NTR 
      J1 = JJ1(I) 
      J2 = JJ2(I) 
      LL1 = LL(J1) 
      LL2 = LL(J2) 
      IS1 = IS(J1) 
      IS2 = IS(J2) 
      L4 = L3 - 2 
      L5 = 2 
      CALL NINEJ(LL2,IS2,JFTOT,L4,L5,L3,LL1,IS1,JITOT,F) 
      FACTR = DFLOAT(L3 + 1) 
      F = F*DSQRT(FACTR) 
      IFL = IABS((LL1+IS1-JITOT)/2) 
      IFR = IABS((LL2+IS2-JFTOT)/2) 
      F = F*(-1.D0)**(IFL+IFR) 
      IF (IBUG1.NE.0) WRITE(IWRITE,505) I,F 
  141 FLINE(I) = F 
* 
*  *****  COMPUTE OSCILLATOR STRENGTHS FOR ALL POSSIBLE PAIRS 
* 
      NI = NJ(1) 
      NF = NJ(2) 
      DO 80 N1 = 1,NI 
      DO 81 N2 = 1,NF 
      JI = ILABEL(1,N1)
      JF = ILABEL(2,N2) + NCF(1)
      ILEN = LENGTH(JI) + LENGTH(JF) + 9
      WRITE(IWRITE,'(/2X,A,2X,A,3H - ,A)')'Line',
     :    CONFIG(JI)(1:LENGTH(JI)), CONFIG(JF)(1:LENGTH(JF))
      WRITE(IWRITE,*) ' ',('-',I=1,ILEN)
      IF (IBUG1.NE.0) WRITE(IWRITE,522) N1,N2 
  522 FORMAT(/,' N1 = ',I3,' N2 = ',I3) 
* 
*  ***** OMIT DUPLICATES IF INITIAL AND FINAL STATES ARE THE SAME 
* 
      IF (CFILE(1) .EQ. CFILE(2) .AND. 
     :    ET(2,N2) - ET(1,N1) .LE. D0 ) GO TO 81 
* 
*  *****  COMPUTE THE LINE STRENGTH FOR THE PAIR 
*  *****  NO EXTRA FACTORS AS IN EQ.(8.15),CH.10 OF SHORE/MENZEL 
* 
      SL = D0 
      DO 130 I = 1,NTR 
      J1 = JJ1(I) 
      J2 = JJ2(I) 
      I1 = II1(I) 
      I2 = II2(I) 
      D = WT(J1,N1)*WT(J2,N2) 
      IF (D .EQ. D0) GO TO 130 
      D = D*FLINE(I) 
      TL = D*TTL(I) 
      IF (IBUG1.NE.0) WRITE(IWRITE,520) I,WT(J1,N1),WT(J2,N2), 
     :                FLINE(I),TTL(I) 
  520 FORMAT(10X,' I = ',I3,' WT1 = ',F10.6,' WT2 = ',F10.6, 
     :       ' FLINE = ',F10.6,' TTL = ',D16.8) 
      IF (PRINT .AND. ABS(TL).GE.sqrt(TOL)) 
     :  WRITE(IWRITE,133) CONFIG(J1)(1:24), CONFIG(J2)(1:24), 
     :                    EL(I1),EL(I2),TL 
  133 FORMAT(1X,A24,' -> ',A24,'  : ',A3,' -> ', A3,F13.7) 
      SL = SL + TL 
      IF (IBUG1.NE.0) WRITE(IWRITE,508) SL 
  130 CONTINUE 
      SL = SL**2 
      SL = SL*(JV(1)+1)*(JV(2)+1) 
      IF (ENDST .EQ. 'M') SL = SL*D4*LAMST*(2*LAMST-1) 
  508 FORMAT(10X,' SL = ',D16.8,/) 
* 
*  *****  COMPUTE THE TRANSITION ENERGIES 
*  *****  DA > 0   1->2  ABSORPTION 
*  *****     < 0   1->2  EMISSION 
* 
      DA = ET(2,N2)-ET(1,N1) 
      IF (DA .EQ. D0) GO TO 536 
      D = DABS(DA) 
      DD = D*D2*RYDBRG 
      ANGS = D10**8/DD 
      ANGSA = ANGS 
      IF (ANGS .GT. 2000.D0) THEN 
         SIGMA = (1.D8/ANGS)**2 
         ANGSA = ANGS/(D1 +8342.13D-8 +2406030./(130.D+8 - SIGMA) 
     :                                 +15997./(38.9D+8 - SIGMA)) 
      END IF 
      IF (PRINT) WRITE(IWRITE,35) DD,ANGS,D 
   35 FORMAT(/10X,'ENERGY DIFFERENCE OF THE STATES:',9X,1PD15.7,' CM-1', 
     :   /51X,1PD15.7,' ANGSTROMS'/51X,1PD15.7,' A. U.') 
* 
*  *****  COMPUTE THE TRANSITION PROBABILITIES (SEC-1) IN EMISSION 
*  *****  AND THE GF(1->2) VALUE 
* 
      I = 2 
      IF (DA .LT. D0) I = 1 
      TRPT = TRP(ENDST,LAMST,DD) 
      AKI = SL*TRPT/(JV(I)+1) 
      IF(IBUG1.NE.0) WRITE(IWRITE,509) SL,TRPT,AKI 
  509 FORMAT(/,1H ,' SL =',1PE15.7,' TRPT =',1PE15.7,' AKI =',1PE15.7) 
      GL = (ANGS**2)*1.4992D-16 
      IF (IBUG1.NE.0) WRITE(IWRITE,525) I,GL,JV(I),D1 
  525 FORMAT(1H ,' I = ',I3,' GL = ',1PD15.7,' J = ',I3, ' D1 = ',F3.0) 
      GL = GL*AKI*(JV(I)+1)*(-D1)**I 
      IF (DABS(GL) .LE. TOL) GO TO 81 
      IF (.NOT.PRINT) GO TO 82 
* 
*  *****  PRINT HEADING FOR THE PAIR (J,J') 
* 
      IC = 1 
      NC = N1 
      DO 70 II = 1,2 
      JC = NCF(II) + IC - 1 
      III = ILABEL(II,NC) + IC - 1 
      WRITE(IWRITE,44) CONFIG(III)(1:LENGTH(III)),JV(II), 
     :   ET(II,NC),(WT(I,NC),I=IC,JC) 
   44 FORMAT(/3X,A,2X,'2J =',I4/3X,'TOTAL ENERGY = ', 
     :   F16.7/(3X,7F10.6)) 
      IC = JC + 1 
      NC = N2 
   70 CONTINUE 
      WRITE(IWRITE,36) SL 
   36 FORMAT(/10X,'SUM OF LENGTH FORM                       = ',1PD15.7) 
      WRITE(IWRITE,234) GL 
  234 FORMAT(10X,'FINAL OSCILLATOR STRENGTH (GF)           = ',1PD15.7) 
      WRITE(IWRITE,134) AKI 
  134 FORMAT(10X,'TRANSITION PROBABILITY IN EMISSION (Aki) =',1PD16.7//) 
* 
*  *****  OUTPUT ON TRANS.LSJ FILE 
* 
   82 JI = ILABEL(1,N1) 
      JF = ILABEL(2,N2) + NCF(1) 
      NOCCI = LENGTH(J1) 
      ITI = 2*NOCCI-1 
      NOCCF = LENGTH(JF) 
      ITF = 2*NOCCF-1 
      WRITE(10,'(/)') 
      WRITE(10,'(I4,F14.8,2X,A)')JV(1),ET(1,N1),CONFIG(JI)(1:LENGTH(JI)) 
      WRITE(10,'(I4,F14.8,2X,A)')JV(2),ET(2,N2),CONFIG(JF)(1:LENGTH(JF)) 
      WRITE(10,38) DD,ANGS,ANGSA,ENDST,LAMST,SL,GL,AKI 
   38 FORMAT(F11.2,' CM-1',2X,F11.2,' ANGS(VAC)',2X,F11.2,' ANGS(AIR)'/ 
     :    1X,A1,I1,2X,'S = ',1PD12.5,3X,'GF = ',D12.5,3X,'AKI = ',D12.5) 
   81 CONTINUE 
   80 CONTINUE 
* 
*  *****  DO IT AGAIN FOR ANOTHER PAIR (J,J') 
* 
  536 IC = ICST 
      II = 2 
      GO TO 65 
* 
*  *****  REWIND STATE2.LSJ FILE IF EOF ON UNIT 8 AND SKIP 
*  *****  THE HEADER CARD 
*  *****  REWIND STATE1.LSJ AND STATE2.LSJ IF EOF ON UNIT 7 AND 
*  *****  CONSIDER THE NEXT SET OF TRANSITION INTEGRALS 
* 
  300 IF (II .EQ. 2) THEN 
         REWIND 8 
         READ(8,'(A)') BIDON 
         GO TO 67 
      ELSE 
         REWIND 7 
         REWIND 8 
         GO TO 71 
      ENDIF 
   99 CLOSE(UNIT=1) 
      CLOSE(UNIT=2) 
      CLOSE(UNIT=7) 
      CLOSE(UNIT=8) 
      CLOSE(UNIT=9) 
      WRITE(10,'(A1)') '*'
      CLOSE(UNIT=10) 
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
*
*     ---------------------------------------------------------------
*        T R P 
*     ---------------------------------------------------------------
*
*
      DOUBLE PRECISION FUNCTION TRP(END,LAM,DD) 
* 
*  *****  THIS FUNCTION CALCULATES THE FACTOR CONNECTING THE TRANSITION 
*  *****  PROBABILITY (SEC-1) AND THE LINE STRENGTH (A.U.) 
* 
*  *****  WE ARE TAKING THE NUMERICAL FACTORS OF PAGES 437,439 OF SHORE 
*  *****  AND MENZEL FOR M AND E TRANSITIONS, AND NOT THOSE OF PAGE 445 
*  *****  DOING THAT, OUR DEFINITION OF LINE STRENGTHS DOES NOT CORRES- 
*  *****  POND TO EQ.(8-15),CHAP.10 BUT RATHER TO THE SQUARE OF THE RME 
*  *****  WITHOUT EXTRA FACTORS (DEF.OF SOBEL'MAN) 
* 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      CHARACTER*1 END 
      DIMENSION DBLFAC(10) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,NP,NCFG,IB,IC,ID 
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS 
      COMMON /INFORM/IREAD,IWRITE,IOUT,ISC(8) 
      COMMON /DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL 
      DATA       C    ,   ALPHA   ,    RYD      ,     PI        / 
     :    2.997925D 10,7.29720D-03,1.0973731D 05,3.141592654D 00/ 
* 
*  *****  DBLFAC(I) = (2I+1)!! 
* 
      DATA DBLFAC/      3.0000000000D 00,1.5000000000D 01 
     :,1.0500000000D 02,9.4500000000D 02,1.0395000000D 04 
     :,1.3513500000D 05,2.0270250000D 06,3.4459425000D 07 
     :,6.5472907500D 08,1.3749310575D 10 
     :/ 
      IF (IBUG1.NE.0) WRITE(IWRITE,501) END,LAM,DD 
  501 FORMAT(/,1H ,' IN TRP ',A1,I2,'  DD =',D16.8) 
      IF (LAM .GT. 10) THEN 
         WRITE(IWRITE,100) 
  100    FORMAT(1X,' THE VECTOR DBLFAC IN TRP HAS TO BE EXTENDED') 
         STOP 
      ENDIF 
      L3 = LAM + LAM 
      L4 = L3 + 1 
      F = DFLOAT(L4*(LAM+1))/LAM 
      F2 = DBLFAC(LAM) 
      IF(IBUG1.NE.0) WRITE(IWRITE,500) F,LAM,F2 
  500 FORMAT(8X,' F =',D16.8,' LAM =',I3,' F2 =',D16.8) 
      F2 = F2**2 
      F = F/F2 
      R = RYD**L3 
      F = F*PI*C/R 
      IF (END .EQ. 'E') THEN 
         F2 = D4/(D2**L3) 
         F2 = F2*ALPHA**L4 
      ELSE 
         F2 = D2**(-L3) 
         F2 = F2*ALPHA**(L3+3) 
      ENDIF 
      F = F*F2 
* 
*  *****  MULTIPLY BY THE WAVENUMBER (CM-1) FUNCTION 
* 
      TRP = F*DD**L4 
      RETURN 
      END 
