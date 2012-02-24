*      
*     Routines for MCHF_LIB_ANG
*
*                   C O P Y R I G H T -- 1994
*
*     Computer Physics Communication, Vol. 64, 399-405 (1991)
*
*     NOTE:  This file contains routines for programs using one atomic 
*            state and two atomic states as in transitions. The COMMONS
*            are different, thus there will be warning messages if this
*            file compiled without first using fsplit.
*     ------------------------------------------------------------------
*       A N A L Y S E 1
*     ------------------------------------------------------------------
*
      SUBROUTINE ANALY1(IREAD,IWRITE,NCLOSD,MAXORB,N,NCFG,NOCCSH,LIST,
     :                   NCD)
*
*        This routine analyzes the format of the configuration input
*        data and determines a consistent ordering of the electrons
*
      PARAMETER (NWD=30)
      INTEGER NOCCSH(NCD),AFTER(NWD,NWD),IEL(5)
      CHARACTER LIST(NWD)*3, LINE*72, OF(NWD)*3, EL(5)*3
*
  1   FORMAT(A72)
*
      DO 2 I = 1,(NWD)
         DO 3 J = 1,(NWD)
            AFTER(I,J) = 0
  3      CONTINUE
  2   CONTINUE
*
*  ---  Determine the number of common closed subshells
*
      READ(IREAD,'(/A72)' ) LINE
      NCLOSD = 0
      J = 2
 10   IF (LINE(J:J+2) .NE. '   ' ) THEN
         NCLOSD = NCLOSD + 1
         J = J+4
         IF (J .LT. 72) GO TO 10
      END IF
*
*  ---  Determine the number or configurations and electrons
*
      MAXORB = 0
      NCFG = N
 20   READ(IREAD,1,END=55) LINE
      IF (LINE(1:1) .NE. '*'  .AND. LINE(2:2) .NE. '*' ) THEN
*
*  ------  A new configuration has been read; find the electrons
*
         NCFG = NCFG + 1
         IF (NCFG .GT. NCD )
     :      WRITE(IWRITE,'(A,I5)') ' TOO MANY CONFIGURATIONS: MAX=',NCD
         J = 2
         I = 0
 30      IF (LINE(J:J+2) .NE. '   ' .AND. I.LT.(5)) THEN
*
*  --------- An electron has been found; is it a new one?
*
            I = I+1
            EL(I) = LINE(J:J+2)
            K = 1
 40         IF (K .LE. MAXORB) THEN
               IF ( OF(K) .NE. EL(I) ) THEN
                  K = K+1
                  IF (K .GT. NWD) THEN
                     WRITE(IWRITE,*) ' TOO MANY ELECTRONS: MAX=',NWD
                     STOP
                  END IF
                  GO TO 40
                 ELSE
                  IEL(I) = K
               END IF
              ELSE
*
*  ------------  A new electron has been found; add it to the list
*
               MAXORB = K
               OF(MAXORB) = EL(I)
               IEL(I) = K
            END IF
            J = J+8
            GO TO 30
         END IF
         NOCCSH(NCFG) = I
*
*  ------  Add data to the AFTER matrix
*
         DO 50 I1 = 2,I
            DO 51 I2 = 1,I1-1
               AFTER(IEL(I1),IEL(I2)) = 1
 51         CONTINUE
 50      CONTINUE
         READ(IREAD,*)
         IF (I .GT. 5) READ(IREAD,*)
         GO TO 20
      END IF
*
*  ---  Check if the ordering of the electrons is inconsistent
*
 55   DO 60 I = 1,MAXORB
         DO 61 J = 1,MAXORB
            IF (AFTER(I,J) .EQ. 1 .AND. AFTER(J,I) .EQ. 1) THEN
                WRITE(IWRITE,*) ' The order of ',OF(I),' and ',
     :                OF(J),' is inconsistent'
                STOP
            END IF
 61      CONTINUE
 60   CONTINUE
*
*  ---  Reorder the electrons to satisfy the after relations found
*         in the different configurations
*
      IORD = 1
 70   IF (IORD .LE. MAXORB ) THEN
*
*  ------  Search for a row with no 1's
*
         DO 71 I = 1,MAXORB
            DO 72 J = 1,MAXORB
               IF (AFTER(I,J) .EQ. 1 ) GO TO 71
 72         CONTINUE
*
*  ---------  The current row contains all 0's or 2's
*
            IF (AFTER(I,I) .NE. 2 ) THEN
*
*  ------------  We have the next electron; delete the corresponding
*                  rows and columns from the AFTER matrix
*
               LIST(IORD) = OF(I)
               IORD = IORD+1
               DO 73 J = 1,MAXORB
                  AFTER(I,J) = 2
                  AFTER(J,I) = 2
 73            CONTINUE
               GO TO 70
            END IF
 71      CONTINUE
      END IF
      RETURN
      END
*
*     ------------------------------------------------------------------
*       A N A L Y S 2
*     ------------------------------------------------------------------
*
      SUBROUTINE ANALY2(NCLOSI,NCLOSF,MCFG,KCFG,LIST,LORTH)
*
*        This routine analyzes the format of the configuration input
*        data, for two sets, not necessarily orthogonal.
*
      PARAMETER (NWD=30,NWD2=2*NWD,NCD=100,NCD4=4*NCD)
      INTEGER AFTER(3*NWD,3*NWD),IEL(5),NORB(2),NCLOS(2),ICFG(2)
      CHARACTER*3 LIST(*), LINE*72, OF(NWD,2), ELC(NWD), EL(5), FIND
      CHARACTER*7 LABEL(2)
      CHARACTER*6 ANS
      LOGICAL LORTH
      COMMON/INFORM/ IREADI,IWRITE,IOUT,IREADF,ISC(7)
      COMMON/STATES/NCFG,MAXORB,IAJCMP(2*NWD),LJCOMP(2*NWD),
     :NJCOMP(2*NWD),NOCCSH(NCD4),NELCSH(5,NCD4),NOCORB(5,NCD4),
     :J1QNRD(9,NCD4)
      COMMON/NOR/NCOM,NORBI,NORBF,IWAR
      DATA LABEL/'Initial','Final  '/
*
  1   FORMAT(A72)
    4 FORMAT(/10H THERE ARE,I3,' INITIAL STATE ORBITALS AS FOLLOWS: '/
     :      (1X,18(1X,A3)))
    5 FORMAT(/10H THERE ARE,I3,' FINAL STATE ORBITALS AS FOLLOWS: '/
     :      (1X,18(1X,A3)))
    6 FORMAT(' List common orbitals, terminating with a blank orbital.'/
     :       ' Upper and lower case characters must match.'/
     :       ' Fixed format (18(1X,A3)) as inicated below:'/
     :       ' AAA AAA AAA AAA AAA AAA AAA .... etc (up to 18/line)')
    7 FORMAT(18(1X,A3))
    8 FORMAT(/10H THERE ARE,I3,' COMMON ORBITALS AS FOLLOWS: '/
     :      (1X,18(1X,A3)))
*
      DO 2 I = 1,(3*NWD)
         DO 3 J = 1,(3*NWD)
            AFTER(I,J) = 0
  3      CONTINUE
  2   CONTINUE
*
      IREAD = IREADI
      NCFG = 0
      DO 100 ISTATE = 1,2
*
*  ---  Determine the number of common closed subshells
*
      READ(IREAD,'(/A72)' ) LINE
      NCLO = 0
      J = 2
 10   IF (LINE(J:J+2) .NE. '   ' ) THEN
         NCLO = NCLO + 1
         J = J+4
         IF (J .LT. 72) GO TO 10
      END IF
      NCLOS(ISTATE) = NCLO
*
*  ---  Determine the number or configurations and electrons
*
      IORB = 0
 20   READ(IREAD,1,END=55) LINE
      IF (LINE(1:1) .NE. '*'  .AND. LINE(2:2) .NE. '*' ) THEN
*
*  ------  A new configuration has been read; find the electrons
*
         NCFG = NCFG + 1
         IF (NCFG .GT. (NCD4) ) THEN
            WRITE(IWRITE,*) ' TOO MANY CONFIGURATIONS: MAX=',NCD4
            STOP
         END IF
         J = 2
         I = 0
 30      IF (LINE(J:J+2) .NE. '   ' .AND. I.LT.(5)) THEN
*
*  --------- An electron has been found; is it a new one?
*
            I = I+1
            EL(I) = LINE(J:J+2)
            K = 1
 40         IF (K .LE. IORB) THEN
               IF ( OF(K,ISTATE) .NE. EL(I) ) THEN
                  K = K+1
                  IF (K .GT. (NWD)) THEN
                     WRITE(IWRITE,*) ' TOO MANY ELECTRONS: MAX=',NWD
                     STOP
                  END IF
                  GO TO 40
                 ELSE
                  IEL(I) = K
               END IF
              ELSE
*
*  ------------  A new electron has been found; add it to the list
*
               IORB = K
               OF(IORB,ISTATE) = EL(I)
               IEL(I) = K
            END IF
            J = J+8
            GO TO 30
         END IF
         NOCCSH(NCFG) = I
*
*  ------  Add data to the AFTER matrix
*
         DO 50 I1 = 2,I
            DO 51 I2 = 1,I1-1
               J1 = (NWD)*ISTATE + IEL(I1)
               J2 = (NWD)*ISTATE + IEL(I2)
               AFTER(J1,J2) = 1
 51         CONTINUE
 50      CONTINUE
         READ(IREAD,*)
         GO TO 20
      END IF
 55   NORB(ISTATE) = IORB
      ICFG(ISTATE) = NCFG
      IREAD = IREADF
  100 CONTINUE
*
*  ---   SET PARAMETERS
*
      NORBI = NORB(1)
      NORBF = NORB(2)
*
*  ---  DETERMINE THE COMMON INITIAL/FINAL STATE ORBITALS
*
      WRITE(0,4) NORBI,(OF(I,1),I=1,NORBI)
      WRITE(0,5) NORBF,(OF(I,2),I=1,NORBF)
      ANS = 'Y'
      IF (.NOT. LORTH) THEN
         WRITE(0,'(/A)')
     :     ' Initial & final state orbitals an orthonormal set ? (Y/N) '
         READ(5,'(A1)') ANS
      END IF
      IF ( ANS .EQ. 'Y' .OR. ANS .EQ. 'y') THEN
         DO 53 I = 1,NORBI
            ELC(I) = OF(I,1)
   53    CONTINUE
         NCOM = NORBI
*
*  ---  ADD OTHERS FROM FINAL STATE
*
         DO 54 I = 1,NORBF
            DO 56 J = 1,NORBI
               IF (OF(I,2) .EQ. OF(J,1)) GO TO 54
   56       CONTINUE
            NCOM = NCOM + 1
            IF (NCOM .GT. (NWD))
     :         STOP ' Too many common electrons: MAX=(30)'
            ELC(NCOM) = OF(I,2)
   54    CONTINUE
        ELSE
         WRITE(0,6)
         READ(5,7) (ELC(I),I=1,18)
         IF (ELC(18) .NE. '   ') READ(5,7) (ELC(I),I=19,(NWD))
         NCOM = 0
   52    IF (ELC(NCOM+1) .NE. '   ') THEN
            NCOM = NCOM + 1
            IF (NCOM .LT. (NWD)) GO TO 52
         END IF
      END IF
      WRITE(0,'(//)')
*
*  ---  Transfer electrons to common orthogonal set
*
      DO 200 ISTATE = 1,2
         IORIG = (NWD)*ISTATE
         LAST = NORB(ISTATE)
         DO 201 I=1,NCOM
*
*        Find electron and transfer AFTER information
*
         J = 1
  202    IF (J .LE. LAST) THEN
            IF (ELC(I) .NE. OF(J,ISTATE) ) THEN
               J = J+1
               GO TO 202
              ELSE
               II = IORIG + J
               DO 210 K = 1, IORIG+LAST
                  IF (AFTER(I,K) .EQ. 0) AFTER(I,K) = AFTER(II,K)
                  IF (AFTER(K,I) .EQ. 0) AFTER(K,I) = AFTER(K,II)
                  AFTER(II,K) = 2
                  AFTER(K,II) = 2
  210          CONTINUE
               NORB(ISTATE) = NORB(ISTATE) - 1
            END IF
           ELSE
            WRITE(0,*) ' Common electron ',ELC(I),' not found in ',
     :             LABEL(ISTATE),' state'
         END IF
  201    CONTINUE
  200 CONTINUE
*
*  ---  Check if the ordering of the electrons is inconsistent
*
      DO 60 I = 1,(NWD)*3
         EL(1) = FIND(I,OF,ELC)
         DO 61 J = 1,(NWD)*3
            EL(2) = FIND(J,OF,ELC)
            IF (AFTER(I,J) .EQ. 1 .AND. AFTER(J,I) .EQ. 1) THEN
                WRITE(0,*) ' The order of ',EL(1),' and ',
     :                EL(2),' is inconsistent'
                STOP
            END IF
 61      CONTINUE
 60   CONTINUE
*
*  ---  Reorder the electrons to satisfy the after relations found
*         in the different configurations
*
      IORD = 1
 70   IF (IORD .LE. NCOM ) THEN
*
*  ------  Search for a row with no 1's in the NCOM rows
*
         DO 71 I = 1,NCOM
            DO 72 J = 1,(NWD)*2+NORBF
               IF (AFTER(I,J) .EQ. 1 ) GO TO 71
 72         CONTINUE
*
*  ---------  The current row contains all 0's or 2's
*
            IF (AFTER(I,I) .NE. 2 ) THEN
*
*  ------------  We have the next electron; delete the corresponding
*                  rows and columns from the AFTER matrix
*
               LIST(IORD) = ELC(I)
               IORD = IORD+1
               DO 74 J = 1,(NWD)*2+NORBF
                  AFTER(I,J) = 2
                  AFTER(J,I) = 2
   74          CONTINUE
               GO TO 70
            END IF
 71      CONTINUE
      END IF
      IF (IORD .NE. NCOM+1) THEN
*
*        SEARCH FOR THE ELECTRON NOT INCLUDED
*
      DO 73 I = 1,NCOM
         IF (AFTER(I,I) .NE. 2) THEN
         DO 75 J = (NWD)+1,(NWD)*2+NORBF
            IF (AFTER(I,J) .EQ. 1) THEN
               WRITE(0,*) ELC(I),' cannot be included in the common set'
               IL = 1
               IF ( J .GT. (NWD)*2 ) IL = 2
               WRITE(0,*) ' Occurs AFTER ',FIND(J,OF,ELC),' in ',
     :                     LABEL(IL),' state'
               STOP
            END IF
   75    CONTINUE
         END IF
   73 CONTINUE
      END IF
*
*  ---  ORDER THE REMAINING ELECTRONS FOR THE INITIAL AND FINAL STATE
*
      LAST = NCOM

      LASTEL = NORBI
      DO 300 ISTATE = 1,2
         LAST = LAST + NORB(ISTATE)
  304    IF (IORD .LE. LAST) THEN
         IORIG = (NWD)*ISTATE
         DO 301 I = IORIG+1, IORIG+LASTEL
            DO 302 J = 1,IORIG+LASTEL
               IF (AFTER(I,J) .EQ. 1) GO TO 301
  302       CONTINUE
*
*           The current row contains no 1's
*
            IF (AFTER(I,I) .NE. 2) THEN
*
*               We have the next electron
*
                IF (IORD.GT.(2*NWD)) THEN
                  WRITE(IWRITE,*) ' Too many electrons: MAX=',2*NWD
                  STOP
                END IF
                LIST(IORD) = OF(I-IORIG,ISTATE)
                IORD = IORD+1
                DO 303 J = 1,IORIG+LASTEL
                   AFTER(I,J) = 2
                   AFTER(J,I) = 2
  303           CONTINUE
                GO TO 304
             END IF
  301    CONTINUE
         END IF
         LASTEL = NORBF
  300 CONTINUE
*
      NORBI = NORB(1)
      NORBF = NORB(2)
      NCLOSI = NCLOS(1)
      NCLOSF = NCLOS(2)
      MCFG = ICFG(1)
      KCFG = ICFG(2) - MCFG
      IF (NCOM .GT. 0) WRITE(IWRITE,8) NCOM,(LIST(I),I=1,NCOM)
      WRITE(IWRITE,4) NORBI,(LIST(I),I=NCOM+1,NCOM+NORBI)
      NOR11 = NCOM + NORBI
      WRITE(IWRITE,5) NORBF,(LIST(I),I=NOR11+1,NOR11+NORBF)
      RETURN
      END
*
*     ------------------------------------------------------------------
*       C F G I N 2
*     ------------------------------------------------------------------
*
      SUBROUTINE CFGIN2(MCFG,KCFG,LORTH,INPUT)
*
*       Read two sets of configurations and determine the orthogonality
*       conditions between them
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NCD=100,NCD4=4*NCD)
      CHARACTER*3 EL(2*NWD), ELC(NWD), JAJCMP(2*NWD,3)*1
      CHARACTER INPUT(2)*24,HEADI*72,HEADF*72,HEADER*72
      LOGICAL LORTH
      COMMON/INFORM/ IREADI,IWRITE,IOUT,IREADF,ISC(7)
      COMMON/STATES/NCFG,MAXORB,IAJCMP(2*NWD),LJCOMP(2*NWD),
     :NJCOMP(2*NWD),NOCCSH(NCD4),NELCSH(5,NCD4),NOCORB(5,NCD4),
     :J1QNRD(9,NCD4)
      COMMON/NOR/NCOM,NORBI,NORBF,IWAR
*
    3 FORMAT(18(1X,A3))
    7 FORMAT(A72)
   22 FORMAT(// 7H STATE ,' (WITH',I3,' CONFIGURATIONS):'/1H ,31(1H-)/)
   23 FORMAT(/10H THERE ARE,I3,21H ORBITALS AS FOLLOWS://
     1 5X,21(1X,A3):/5X,21(1X,A3))
*
*     The "readonly" option is needed on some computers when the
*     two files are in fact the same.  Others ignore the option 
*     which is OK in most cases.
*     Microsoft Fortran requires "READ" instead of "READONLY"
*     OPEN(UNIT=1,FILE=INPUT(1),STATUS='OLD',READONLY)
*     OPEN(UNIT=2,FILE=INPUT(2),STATUS='OLD',READONLY)
      OPEN(UNIT=1,FILE=INPUT(1),STATUS='OLD')
      OPEN(UNIT=2,FILE=INPUT(2),STATUS='OLD')
*
* --- ANALYZE INITIAL AND FINAL STATE DATA
*
      CALL ANALY2(NCLOSI,NCLOSF,MCFG,KCFG,EL,LORTH)
      REWIND(UNIT=IREADI)
      REWIND(UNIT=IREADF)
*
      MAXORB = NCOM + NORBI + NORBF
*   SET UP THE ELECTRONS
*
      READ(EL,'(A3)') (IAJCMP(I),I=1,MAXORB)
      READ(EL,'(3A1)')((JAJCMP(I,J),J=1,3),I=1,MAXORB)
*
*   SET UP OF LJCOMP
*
      DO 60 I = 1,MAXORB
      IF (JAJCMP(I,1) .EQ. ' ') THEN
         JAJCMP(I,1) = JAJCMP(I,2)
         JAJCMP(I,2) = JAJCMP(I,3)
         JAJCMP(I,3) = ' '
      ENDIF
      LJCOMP(I) = LVAL(JAJCMP(I,2))
      NJCOMP(I) = ICHAR(JAJCMP(I,1)) - ICHAR('1') + 1
   60 CONTINUE
*
* ---- CHECK COMMON CLOSED SHELLS
*
      IF (NCLOSI .NE. NCLOSF)
     :   STOP ' Common closed shells not the same in the two states'
*
      READ(IREADI,7) HEADI
      READ(IREADF,7) HEADF
      HEADER = HEADI(1:34)//'=>'//HEADF(1:34)
      WRITE(IOUT,7) HEADER
*
* --- CHECK CLOSED SHELLS FURTHER
*
      READ(IREADI,3) (ELC(I),I=1,NCLOSI)
      READ(IREADF,3) (EL(I),I=1,NCLOSF)
      DO 1 I = 1,NCLOSF
         J = 1
    2    IF (EL(I) .NE. ELC(J) ) THEN
            J = J+1
            IF (J .LE. NCLOSI) THEN
               GO TO 2
              ELSE
               STOP ' Common closed sub-shells not the same'
            END IF
         END IF
    1 CONTINUE
*
*  MAXORB < 2*NWD+1... LINKED TO JAJCMP(2*NWD,3)
*                            IAJCMP(2*NWD)
*                            LJCOMP(2*NWD)
*                            NJCOMP(2*NWD)
*  THE DIMENSION OF IORTH IS GIVEN BY THE PRODUCT OF THE ALLOWED
*  NORBI AND NORBF, I.E. ACTUALLY NWD by NWD = NWD^2
*
*
*   GET INITIAL STATE CONFIGURATIONS
*
      CALL GSTATE(1,MCFG)
*
*
      MCFG1 = MCFG + 1
      NCFG = MCFG + KCFG
      CALL GSTATE(MCFG1,NCFG)
*
*  ---  CHECK THE DATA
*
      CALL CFGTST(NCFG,LJCOMP,NOCCSH,NELCSH,NOCORB,J1QNRD,NCD4)
      RETURN
      END
*
*     ------------------------------------------------------------------
*       C F G N 1
*     ------------------------------------------------------------------
*
*       Read the configurations for a state and determine the
*       non-orthogonal orbitals
*
      SUBROUTINE CFGN1(INPUT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NCD=100)
*
      CHARACTER BUFFER*3
      CHARACTER*1 JAJCMP(NWD,3), INPUT*24
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,JSC0,JSC1,
     :JSC2,JSC3
      COMMON/OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
     : ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
     :     IORTH(NWD*(NWD-1)/2)
      COMMON/STATES/NCFG,MAXORB,IAJCMP(NWD),LJCOMP(NWD),
     :NJCOMP(NWD),NOCCSH(NCD),NELCSH(5,NCD),NOCORB(5,NCD),J1QNRD(9,NCD)
      INTEGER IEL(2),IBUFF(2)
      DATA IASTER,IBLANK/3H*  ,3H   /
*
      CALL CFGO1(NCFG,MAXORB,IAJCMP,LJCOMP,NJCOMP,NOCCSH,NELCSH,
     :            NOCORB,J1QNRD,NCD,INPUT)
*
*  ---  SEPARATE THE ELECTRON LABEL CHARACTERS AND LEFT JUSTIFY
*
      DO 10 I = 1,MAXORB
         WRITE(BUFFER,'(A3)') IAJCMP(I)
         READ(BUFFER,'(3A1)') (JAJCMP(I,J),J=1,3)
         IF (JAJCMP(I,1) .EQ. ' ') THEN
            JAJCMP(I,1) = JAJCMP(I,2)
            JAJCMP(I,2) = JAJCMP(I,3)
            JAJCMP(I,3) = ' '
         END IF
10    CONTINUE
*
*  ---  INITIALIZE THE ORTHOGONALITY ARRAY
*
      M1 = (MAXORB*(MAXORB-1))/2
      DO 20 I = 1,M1
         IORTH(I) = -1
20    CONTINUE
*
*  ---  SET ORBITALS IN THE SAME CONFIGURATION TO BE ORTHOGONAL
*
      DO 63 I = 1,NCFG
      N = NOCCSH(I)
      DO 66 J = 1,N-1
      DO 66 JJ = J+1,N
         I1 = NOCORB(J,I)
         J1 = NOCORB(JJ,I)
         IF (J1 .GT. I1) THEN
            M = I1
            I1 = J1
            J1 = M
         ENDIF
         IORTH(J1+((I1-1)*(I1-2))/2) = 0
   66 CONTINUE
   63 CONTINUE
*
* --- DETERMINE THE NON-ORTHOGONAL ORBITALS
*
      NORTH = 0
      DO 18 J = 1,MAXORB-1
      DO 19 I = J+1,MAXORB
         IJ = J + ((I-1)*(I-2))/2
         IF (JAJCMP(I,2) .EQ. JAJCMP(J,2) .AND.
     :       JAJCMP(I,3) .NE. ' ' .AND.
     :       JAJCMP(J,3) .NE. ' ' .AND.
     :       JAJCMP(I,3) .NE. JAJCMP(J,3) .AND.
     :       IORTH(IJ) .NE. 0 ) THEN
             NORTH = NORTH + 1
             IORTH(IJ) = 1
         ENDIF
19    CONTINUE
18    CONTINUE
*
      READ(IREAD,*,END=90)
 79   READ(IREAD,'(2(1X,A3))',END=90) IBUFF(1),IBUFF(2)
      IF (IBUFF(1) .NE. IASTER .AND. IBUFF(1) .NE. IBLANK) THEN
         DO 80 I = 1,2
            DO 81 J = 1,MAXORB
               IF (IBUFF(I) .EQ. IAJCMP(J)) THEN
                  IEL(I) = J
                  GO TO 80
               ENDIF
 81         CONTINUE
            WRITE(*,'(A,A3,A)') ' ELECTRON ',IBUFF(I),' NOT FOUND'
            STOP
 80      CONTINUE
         IF (IEL(1) .GT. IEL(2) ) THEN
            I = IEL(1)
            IEL(1) = IEL(2)
            IEL(2) = I
         END IF
         IJ = IEL(1) + ((IEL(2)-1)*(IEL(2)-2))/2
         IF (IORTH(IJ) .EQ. 1) NORTH = NORTH - 1
         IORTH(IJ) = 0
         WRITE(IWRITE,'(1X,A3,A,A3)')
     :           IBUFF(1),' is orthogonal to ',IBUFF(2)
         GO TO 79
      END IF
  90  RETURN
      END
*
*     ------------------------------------------------------------------
*       C F G O 1
*     ------------------------------------------------------------------
*
*       Read configurations for one state, assuming orthogonality of
*       the orbitals
*
      SUBROUTINE CFGO1(NCFG,MAXORB,IAJCMP,LJCOMP,NJCOMP,NOCCSH,
     :                  NELCSH,NOCORB,J1QNRD,NCD,INPUT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30)
*
      CHARACTER BUFFER(NWD)*3, HEADER*72, INPUT*24
      CHARACTER*1 JAJCLD(NWD,3),JAJCMP(NWD,3),JCQN(9)
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,JSC0,JSC(3)
      COMMON /CLOSED/B1ELC(4),NCLOSD,IAJCLD(NWD),LJCLSD(NWD),IBK
*
      DIMENSION IAJCMP(NWD),LJCOMP(NWD),NJCOMP(NWD),NOCCSH(NCD),
     :          NELCSH(5,NCD),NOCORB(5,NCD),J1QNRD(9,NCD),
     :          J3QN(9),J2QN(9),J1QN(9)
*
    3 FORMAT(18(1X,A3))
    4 FORMAT(3A1)
*   5 FORMAT(9(1X,A3,'(',I2,')'))
    5 FORMAT(9(1X,A3,1X,I2,1X))
    6 FORMAT(9(1X,4X,I1,A1,I1))
    7 FORMAT(A72)
    8 FORMAT(A3)
   22 FORMAT(// 7H STATE ,' (WITH',I3,' CONFIGURATIONS):'/1H ,31(1H-)/)
   23 FORMAT(/10H THERE ARE,I3,21H ORBITALS AS FOLLOWS://
     : 5X,21(1X,A3):/5X,21(1X,A3))
   25 FORMAT(/14H CONFIGURATION,I3,' ( OCCUPIED ORBITALS=',I2,' ):'
     : ,5(1X,A3,1H(,I2,1H)))
   26 FORMAT(26X,17H COUPLING SCHEME:,5(1X,4X,I1,A1,I1))
   27 FORMAT(54X,4(1X,4X,I1,A1,I1))
   28 FORMAT(/10H THERE ARE ,I3,31H CLOSED SUBSHELLS COMMON TO ALL ,
     :  27H CONFIGURATIONS AS FOLLOWS: //
     :  5X, 21(1X,A3))
*
* --- ANALYZE INPUT DATA
*
      OPEN(UNIT=4,FILE=INPUT,STATUS='OLD')
      CALL ANALY1(IREAD,IWRITE,NCLOSD,MAXORB,0,NCFG,NOCCSH,BUFFER,NCD)
      REWIND(UNIT=IREAD)
*
* ---  Process the configuration data
*
      READ(BUFFER,8) (IAJCMP(I),I=1,MAXORB)
      WRITE(IWRITE,22) NCFG
      WRITE(IWRITE,23) MAXORB,(IAJCMP(I),I=1,MAXORB)
      READ(BUFFER,4)((JAJCMP(I,J),J=1,3),I=1,MAXORB)
      DO 60 I=1,MAXORB
      IF (JAJCMP(I,1) .EQ. ' ') THEN
         JAJCMP(I,1) = JAJCMP(I,2)
         JAJCMP(I,2) = JAJCMP(I,3)
         JAJCMP(I,3) = ' '
      ENDIF
      LJCOMP(I) = LVAL(JAJCMP(I,2))
      NJCOMP(I) = ICHAR(JAJCMP(I,1)) - ICHAR('1') + 1
   60 CONTINUE
*
* --- READ HEADER CARD FOR THE CASE
*
      READ(IREAD,7) HEADER
      WRITE(IOUT,7) HEADER
*
* --- READ IN THE COMMON SET OF CLOSED SUBSHELLS
*
      READ(IREAD,3) (BUFFER(I),I=1,NCLOSD)
      IF (NCLOSD .EQ. 0) GO TO 70
      READ(BUFFER,8) (IAJCLD(I),I=1,NCLOSD)
      WRITE(IWRITE,28) NCLOSD,(IAJCLD(I),I=1,NCLOSD)
      READ(BUFFER,4) ((JAJCLD(I,J),J=1,3),I=1,NCLOSD)
      DO 71 I=1,NCLOSD
      J = 3
      IF (JAJCLD(I,1) .NE. ' ') J = 2
      LJCLSD(I) = LVAL(JAJCLD(I,J))
71    CONTINUE
 70   CONTINUE
*
* --- READ IN (AND PRINT OUT) CONFIGURATIONS ETC. FOR THE STATE UNDER
* --- CONSIDERATION
*
      DO 63 I=1,NCFG
      N=NOCCSH(I)
      READ(IREAD,5)        (NOCORB(J,I),NELCSH(J,I),J=1,N)
      WRITE(IWRITE,25) I,N,(NOCORB(J,I),NELCSH(J,I),J=1,N)
      DO 61 J=1,N
      DO 61 JJ=1,MAXORB
   61 IF(NOCORB(J,I).EQ.IAJCMP(JJ)) NOCORB(J,I)=JJ
      M=2*N-1
      N1=N+1
      READ(IREAD,6)    (J3QN(J),JCQN(J),J1QN(J),J=1,M) 
      WRITE(IWRITE,26) (J3QN(J),JCQN(J),J1QN(J),J=1,N) 
      IF(N.GT.1) WRITE(IWRITE,27) (J3QN(J),JCQN(J),J1QN(J),J=N1,M)
      DO 62 J=1,M
      J2QN(J) = 2*LVAL(JCQN(J)) + 1
      J1QNRD(J,I) = (J3QN(J)*64 + J2QN(J))*64 + J1QN(J)
   62 CONTINUE
   63 CONTINUE
      CALL CFGTST(NCFG,LJCOMP,NOCCSH,NELCSH,NOCORB,J1QNRD,NCD)
      RETURN
      END
*
*     ------------------------------------------------------------------
*       C F G T S T
*     ------------------------------------------------------------------
*
      SUBROUTINE CFGTST(NCFG,LJCOMP,NOCCSH,NELCSH,NOCORB,J1QNRD,NCD)
*
*     THIS SUBROUTINE CHECKS ALL THE CONFIGURATION SET TO ENSURE THAT
*     IT SATISFIES ALL THE FOLLOWING CONDITIONS:
*        (1)  EACH CONFIGURATION HAS THE SAME NUMBER OF ELECTRONS
*        (2)  NO SUBSHELL HAS TOO MANY (.GT.2*(2*L+1))  ELECTRONS
*        (3)  THE ELECTRONS IN ANY ONE SUBSHELL ARE COUPLED TO FORM AN
*             ALLOWED TRIAD OF QUANTUM NUMBERS
*        (4)  THE TRIADS COUPLE TOGETHER IN AN ALLOWED WAY
*
*     IN THE EVENT OF AN ERROR, THE PROGRAM HALTS AT THE COMPLETION
*     OF THE CHECKING.  ANY NUMBER OF S, P, D  ELECTRONS ARE ALLOWED,
*     (BUT .LE.2*(2*L+1)), BUT ONLY UP TO TWO ELECTRONS, L >=3.
*     WHEN L>4, THE ONLY ALLOWED TERMS ARE THOSE FOR L=4.
*     A FILLED F-SHELL IS ALSO ALLOWED AS WELL AS A SINGLE ELECTRON
*     WITH L.GT.4
*
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,JSC0,JSC(3)
      COMMON/TERMS/NROWS,ITAB(24),JTAB(24),NTAB(333)
*
      DIMENSION LJCOMP(*),NOCCSH(NCD),NELCSH(5,NCD),NOCORB(5,NCD),
     :          J1QNRD(9,NCD)
*
    5 FORMAT(/38H THE TRIAD OF QUANTUM NUMBERS OF SHELL,I3,17H IN CONFIG
     :URATION,I3,24H IS NOT A RECOGNIZED SET)
    7 FORMAT(/22H THE COUPLING OF SHELL,I3,17H IN CONFIGURATION,I3,
     : 38H RESULTS IN AN ILLEGAL COUPLING SCHEME)
   12 FORMAT(//41H CONFIGURATION DATA WRONG, PROGRAM HALTED//)
   15 FORMAT(/17H IN CONFIGURATION,I3,7H, SHELL,I3,28H CONTAINS TOO MANY
     : ELECTRONS)
   17 FORMAT(/14H CONFIGURATION,I3,68H INCLUDES A SHELL OF ANGULAR MOMEN
     :TUM L.GE.3 WITH TOO MANY ELECTRONS)
   18 FORMAT(/14H CONFIGURATION,I3,28H HAS AN INCORRECT NUMBER OF ,
     :        9HELECTRONS)
*
      IALLOW=1
      DO 1 I=1,NCFG
         NELSUM = 0
         N=NOCCSH(I)
         DO 2 J=1,N
            NA=NOCORB(J,I)
            LQU=LJCOMP(NA)
            NC=NELCSH(J,I)
            NELSUM = NELSUM + NC
            JD = J1QNRD(J,I)
            JA = MOD(JD,64)
            JD = JD/64
            JB = MOD(JD,64)
            JC = JD/64
            LQUMAX = 4*LQU + 2
            IF (NC .GT. LQUMAX) THEN
               WRITE(IWRITE,15) I,J
               IALLOW = 0
               GO TO 2
            ELSE IF ((LQU.EQ.3 .AND. NC.GT.2 .AND. NC.LT.14) .OR.
     :          (LQU.GT.4.AND.NC.GT.2)) THEN
               WRITE(IWRITE,17) I
               IALLOW = 0
               GO TO 2
            ELSE IF (NC .EQ. 1) THEN
               IF (JA.EQ.1 .AND. JB.EQ.(2*LQU+1) .AND. JC.EQ.2) GO TO 21
            ELSE
               IF (LQU .GT. 4 .AND. NC .EQ. 2) LQU = 4
               IF (NC .EQ. LQUMAX) THEN
                  NROW = 2
               ELSE
                  NROW = NTAB1(NC+1,LQU+1)
               END IF
               I1 = ITAB(NROW)
               I2 = JTAB(NROW)
               DO 4 IA = 1,I1
                  I3 = I2+3*IA-1
                  IF (JB .EQ. NTAB(I3)) THEN
                     I3 = I3+1
                     IF (JC .EQ. NTAB(I3)) THEN
                        I3 = I3-2
                        IF (JA .EQ. NTAB(I3)) GO TO 21
                     END IF
                  END IF
    4          CONTINUE
            END IF
            IALLOW = 0
            WRITE(IWRITE,5) J,I
            GO TO 2
*
*     CHECK ON THE COUPLING OF THE TRIADS
*
   21       IF (N.GT.1 .AND. J.GT.1) THEN
               J2 = N+J-1
               J1 = J2-1
               IF (J.EQ.2) J1 = 1
               JE = J1QNRD(J1,I)/64
               JD = MOD(JE,64)
               JE = JE/64
               JG = J1QNRD(J2,I)/64
               JF = MOD(JG,64)
               JG = JG/64
               IF (JF.GE.(JB+JD) .OR. JF.LE.IABS(JB-JD) .OR.
     :             JG.GE.(JC+JE) .OR. JG.LE.IABS(JC-JE) .OR.
     :             MOD(JC+JE-JG,2).EQ.0 ) THEN
                   WRITE(IWRITE,7) J,I
                   IALLOW = 0
               END IF
            END IF
    2    CONTINUE
         IF (I .EQ. 1) THEN
            NELCS = NELSUM
         ELSE IF (NELSUM .NE. NELCS) THEN
            WRITE(IWRITE,18) I
            IALLOW = 0
         END IF
    1 CONTINUE
      IF (IALLOW .EQ. 0) THEN
         WRITE(IWRITE,12)
         STOP
      END IF
      END
*
*     ------------------------------------------------------------------
*       B L O C K    D A T A
*     ------------------------------------------------------------------
*
      BLOCK DATA INITT
*
      COMMON/KRON/IDEL(10,10)
      COMMON/TERMS/NROWS,I(24),J(24),N(333)
*
* --- SETS QUANTUM NUMBERS OF TERMS WHICH CAN BE FORMED FROM
*     CONFIGURATIONS  L**Q . ONLY THE FIRST HALF OF THAT PART OF THE
*     TABLE, CORRESPONDING TO A GIVEN  L, IS INCLUDED, BECAUSE OF THE
*     SYMMETRY OF THE TABLE.  E.G. D**7 FORMS THE SAME TERMS AS D**3
*
*       The tables are set for a maximum value of L=9; the terms
*     for L>3 are assumed to be the same as those for L=3
*
*     S - SHELLS (ROWS 1 AND 2)
*
*     P - SHELLS (ROWS 3 TO 5)
*
*     D - SHELLS (ROWS 6 TO 10)
*
*     F - SHELLS (ROWS 11 AND 12)
*
*     G - SHELLS (ROWS 13 AND 14)
*
*     H - SHELLS (ROWS 15 AND 16)
*
*     I - SHELLS (ROWS 17 AND 18)
*
*     K - SHELLS (ROWS 19 AND 20)
*
*     L - SHELLS (ROWS 21 AND 22)
*
*     M - SHELLS (ROWS 23 AND 24)
*
      DATA NROWS/24/
*
*     THE ARRAYS I,J,N CORRESPOND TO THE ARRAYS ITAB,JTAB,NTAB
*
      DATA I/1,1, 1,3,3, 1,5,8,16,16, 1,7, 1,7, 1,7, 1,7, 1,7, 1,7, 1,7/
      DATA J/0,3, 6,9,18, 27,30,45,69,117, 165,168, 189,192,
     : 213,216, 237,240, 261,264, 285,288, 309,312/
      DATA N/1,1,2,  0,1,1,  1,3,2,  0,1,1, 2,5,1, 2,3,3, 1,3,2,
     : 3,5,2, 3,1,4,  1,5,2,  0,1,1, 2,5,1, 2,9,1, 2,3,3, 2,7,3,
     : 1,5,2,  3,3,2, 3,5,2, 3,7,2, 3,9,2, 3,11,2, 3,3,4, 3,7,4,
     : 0,1,1, 2,5,1, 2,9,1, 2,3,3, 2,7,3, 4,1,1, 4,5,1, 4,7,1, 4,9,1,
     : 4,13,1, 4,3,3, 4,5,3, 4,7,3, 4,9,3, 4,11,3, 4,5,5,
     : 1,5,2, 3,3,2, 3,5,2, 3,7,2, 3,9,2, 3,11,2, 3,3,4, 3,7,4, 5,1,2,
     : 5,5,2, 5,7,2, 5,9,2, 5,13,2, 5,5,4, 5,9,4, 5,1,6,
     : 1,7,2,  2,3,3, 2,7,3, 2,11,3 ,0,1,1, 2,5,1, 2,9,1 ,2,13,1,
     : 1,9,2,  2,3,3, 2,7,3, 2,11,3 ,0,1,1, 2,5,1, 2,9,1 ,2,13,1,
     : 1,11,2,  2,3,3, 2,7,3, 2,11,3 ,0,1,1, 2,5,1, 2,9,1 ,2,13,1,
     : 1,13,2,  2,3,3, 2,7,3, 2,11,3 ,0,1,1, 2,5,1, 2,9,1 ,2,13,1,
     : 1,15,2,  2,3,3, 2,7,3, 2,11,3 ,0,1,1, 2,5,1, 2,9,1 ,2,13,1,
     : 1,17,2,  2,3,3, 2,7,3, 2,11,3 ,0,1,1, 2,5,1, 2,9,1 ,2,13,1,
     : 1,19,2,  2,3,3, 2,7,3, 2,11,3 ,0,1,1, 2,5,1, 2,9,1 ,2,13,1/
*
* --- READ IN OTHER INITIALIZATION DATA
*
      DATA IDEL/1,10*0,1,10*0,1,10*0,1,10*0,1,10*0,1,10*0,1,10*0,
     :          1,10*0,1,10*0,1/
*
      END
*
*     ------------------------------------------------------------------
*       C F P
*     ------------------------------------------------------------------
*
      SUBROUTINE CFP(LIJ,N,IVI,ILI,ISI,IVJ,ILJ,ISJ,COEFP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(8)
*     MicroSoft FORTRAN does not allow the EXTERNAL declaration here
      EXTERNAL INITT
*
* === CHOOSES APPROPRIATE FRACTIONAL PARENTAGE SUBROUTINE
*
    9 FORMAT(69H UNNECESSARY ATTEMPT TO FORM CFP OF AN S-ELECTRON - THER
     :E IS AN ERROR)
      K=LIJ+1
      IF (K .GT. 4) K = 4
*
*     IF F-SHELL OR G-SHELL COEFFICIENT-OF-FRACTIONAL-PARENTAGE ROUTINES
*     ARE INCLUDED, THIS COMPUTED GO TO NEEDS MODIFYING TO ACCOUNT FOR
*     THIS
*
        
      GO TO (1,2,3,4,4) K
*
* --- FALSE CALL FOR S-SHELLS
*
    1 WRITE(IWRITE,9)
      STOP
*
* --- P-SHELLS
*
    2 CALL CFPP(N,ILI,ISI,ILJ,ISJ,COEFP)
      RETURN
*
* --- D-SHELLS
*
    3 CALL CFPD(N,IVI,ILI,ISI,IVJ,ILJ,ISJ,COEFP)
      RETURN
*
* --- F-SHELLS, G-SHELLS ETC. WITH UP TO TWO ELECTRONS
*
    4 CALL CFPF(N,IVI,ILI,ISI,IVJ,ILJ,ISJ,COEFP)
      RETURN
      END
*
*     ------------------------------------------------------------------
*       C F P D
*     ------------------------------------------------------------------
*
      SUBROUTINE CFPD(N,IVI,LI,ISI,IVJ,LJ,ISJ,COEFP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
*
*     THIS SUBROUTINE EVALUATES THE COEFFICIENTS OF FRACTIONAL PARENTAGE
*     FOR EQUIVALENT D SHELL ELECTRONS FROM TABLES GIVEN IN J.C.SLATER
*     QUANTUM THEORY OF ATOMIC STRUCTURE,VOLUME2,P350(1960)
*     IN THE SUBROUTINE LIST N,THE NO.OF ELECTRONS,V THE SENIORITY QUAN
*     TUM NO.,L THE ANGULAR MOMENTUM QUANTUM NO.,(2S+1) THE SPIN QUANTUM
*     NO. OF BOTH THE STATE IN QUESTION AND ITS PARENT STATE ARE INPUT
*     PARAMETERS  THE RESULT IS OUTPUT AS COEFP
*
*
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(8)
      DIMENSION      IV(5,16),IL(5,16),IS(5,16),
     :          ITAB1(5,1),ITAB2(8,5),ITAB3(16,8),ITAB4(16,16),
     :          NORM1(5),NORM2(8),NORM3(16),NORM4(16)
      DATA IV/1,2,3,4,5,0,2,3,4,5,0,2,3,4,3,0,2,3,2,5,0,0,3,4,3,0,0,1,4,
     :   5,0,0,3,2,3,0,0,3,4,3,0,0,0,4,5,0,0,0,2,3,0,0,0,4,5,0,0,0,4,1,
     :0,0,0,2,3,0,0,0,4,5,0,0,0,0,3,0,0,0,4,5/
      DATA IL/2,3,3,2,0,0,1,1,5,4,0,4,5,4,3,0,2,4,3,2,0,0,3,3,1,0,0,2,2,
     :   6,0,0,2,1,5,0,0,1,1,4,0,0,0,6,4,0,0,0,4,3,0,0,0,4,3,0,0,0,3,2,
     :   0,0,0,2,2,0,0,0,2,2,0,0,0,0,1,0,0,0,0,0/
      DATA IS/2,3,4,5,6,0,3,4,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,0,2,3,
     :   2,0,0,2,3,2,0,0,2,3,2,0,0,0,1,2,0,0,0,1,2,0,0,0,1,2,0,0,0,1,2,
     :   0,0,0,1,2,0,0,0,1,2,0,0,0,1,2,0,0,0,1,2/
      DATA ITAB1/1,1,1,1,1/
      DATA ITAB2/4,-7,-1,21,7,-21,21,-8,-1,-8,0,0,28,-9,-49,7,0,0,1,11,
     :   -25,-9,-25,0,0,0,0,-10,-10,-5,45,15,0,0,0,0,0,16,0,0/
      DATA ITAB3/7,20,-560,224,-112,-21,-56,16,0,0,0,0,0,0,0,0,3,0,0,-56
     :   ,-448,49,-64,-14,0,0,0,0,0,0,0,0,0,26,308,110,220,0,0,0,7,-154,
     :   -28,-132,0,0,0,0,0,-9,297,90,-405,45,0,0,3,66,-507,-3,-60,15,
     :    0,0,0,5,315,-14,-175,-21,-56,-25,0,70,385,-105,28,63,0,0,0,0,
     :   0,315,0,0,135,0,0,189,0,0,105,0,1,0,0,0,200,15,120,60,-35,10,0,
     :   -25,88,200,45,20,0,1,0,0,0,16,-200,-14,-14,25,0,0,0,120,-42,42,
     :    0,0/
      DATA ITAB4/1,-105,-175,-175,-75,12*0,154,-110,0,0,231,286,924,-308
     :   ,220,-396,6*0,-66,-90,180,0,99,-99,891,-5577,-405,-9,0,45,45,0,
     :   0,0,0,224,0,-56,0,-220,1680,0,112,0,-21,21,0,-16,0,0,-70,14,-84
     :   ,56,0,55,945,4235,-175,-315,0,-21,189,-25,0,0,25,-15,-135,35,0,
     :   0,600,968,120,600,0,60,60,10,3,0,0,-56,0,-64,4*0,448,0,-9,-49,
     :   0,14,0,0,0,-16,126,14,4*0,-200,360,0,-14,126,25,0,5*0,-175,182,
     :   -728,-2184,7*0,6*0,220,880,0,-400,0,-9,-25,0,0,0,5*0,-45,-5,845
     :    ,-1215,275,495,0,-11,99,0,0,6*0,33,-7,-2541,105,-525,0,35,35,
     :   -15,0,7*0,-800 ,0,-160,0,-5,45,0,30,0,7*0,-100,1452,180,-100,0,
     :   -10,90,15,-2,11*0,6,16*0,-14,-56,0,0/
      DATA NORM1/1,1,1,1,1/
      DATA NORM2/5,15,2,42,70,60,140,30/
      DATA NORM3/10,60,1680,840,1680,210,360,90,10,504,1008,560,280,140,
     :   1,1/
      DATA NORM4/1,420,700,700,300,550,1100,8400,18480,2800,2800,50,350,
     :   700,150,5/
*
*     READ IN D SHELL PARAMETERS AND TABLES
*     PERIPHERAL 1 IS THE CARD READER
*
*     TEST IF N IS IN THE FIRST HALF OF SHELL
*
99    IF(N-6) 40,103,103
*
*     TEST IF STATE IN QUESTION IS ALLOWED
*     IF IT IS, IDENTIFY THE ROW OF THE TABLE BY J1
*
40    J = 0
101   J = J+1
      IF(J-17) 41,11,11
41    IF(IV(N,J)-IVI) 101,42,101
42    IF(IL(N,J)-LI) 101,43,101
43    IF(IS(N,J)-ISI) 101,44,101
44    J1=J
*
*     TEST IF PARENT STATE IS ALLOWED
*     IF IT IS, IDENTIFY THE COLUMN OF THE TABLE BY J2
*
      IF(N-1) 45,30,45
30    IF(IVJ) 11,31,11
31    IF(LJ) 11,32,11
32    IF(ISJ-1) 11,1,11
45    J = 0
102   J = J+1
      IF(J-17) 46,11,11
46    IF(IV(N-1,J)-IVJ) 102,47,102
47    IF(IL(N-1,J)-LJ)  102,48,102
48    IF(IS(N-1,J)-ISJ) 102,49,102
49    J2=J
      GO TO 100
*
*     SIMILAR SETTING OF J1 AND J2 IF N IS IN SECOND HALF OF SHELL
*
103   M = 10-N
      IF(M) 36,33,36
33    IF(IVI) 11,34,11
34    IF(LI) 11,35,11
35    IF(ISI-1) 11,37,11
36    J = 0
104   J = J+1
      IF(J-17) 50,11,11
50    IF(IV(M,J)-IVI) 104,51,104
51    IF(IL(M,J)-LI) 104,52,104
52    IF(IS(M,J)-ISI) 104,53,104
53    J1=J
37    J = 0
105   J = J+1
      IF(J-17) 54,11,11
54    IF(IV(M+1,J)-IVJ) 105,55,105
55    IF(IL(M+1,J)-LJ)  105,56,105
56    IF(IS(M+1,J)-ISJ) 105,57,105
57    J2=J
*
*     IDENTIFY THE F.P.C AS A UNIQUE ELEMENT OF ITABN(J1,J2)
*
100   GO TO (1,2,3,4,5,12,12,12,12,1),N
1     COEFP = 1.0D0
      GO TO 10
2     COEFP = ITAB1(J1,J2)
      IF(COEFP) 60,10,81
60    COEFP = - DSQRT(-COEFP/NORM1(J1))
      GO TO 10
81    COEFP = DSQRT(COEFP/NORM1(J1))
      GO TO 10
3     COEFP = ITAB2(J1,J2)
      IF(COEFP) 61,10,82
61    COEFP = -DSQRT(-COEFP/NORM2(J1))
      GO TO 10
82    COEFP = DSQRT(COEFP/NORM2(J1))
      GO TO 10
4     COEFP = ITAB3(J1,J2)
      IF(COEFP) 62,10,83
62    COEFP = -DSQRT(-COEFP/NORM3(J1))
      GO TO 10
83    COEFP = DSQRT(COEFP/NORM3(J1))
      GO TO 10
5     COEFP = ITAB4(J1,J2)
      IF(COEFP) 63,10,84
63    COEFP = -DSQRT(-COEFP/NORM4(J1))
      GO TO 10
84    COEFP = DSQRT(COEFP/NORM4(J1))
      GO TO 10
*
*     USE RECURRENCE RELATION EQUATION (19) OF RACAH FOR SECOND HALF OF
*     SHELL
*
12    ISIGN = (-1)**((ISI+ISJ-7)/2 +LI +LJ)
      FACTOR = DSQRT(DFLOAT((11-N)*ISJ*(2*LJ+1))/DFLOAT(N*ISI*(2*LI+1)))
      M1 =N-5
      GO TO(6,7,8,9),M1
6     COEFP = ITAB4(J2,J1)
      IF(COEFP) 64,10,85
64    COEFP = -DSQRT(-COEFP/NORM4(J2))
      GO TO 86
85    COEFP = DSQRT(COEFP/NORM4(J2))
86    COEFP = COEFP*ISIGN*FACTOR
      IF(MOD((IVJ-1)/2,2))  87,10,87
87    COEFP = -COEFP
      GO TO 10
7     COEFP = ITAB3(J2,J1)
      IF(COEFP) 65,10,88
65    COEFP = -DSQRT(-COEFP/NORM3(J2))
      GO TO 89
88    COEFP = DSQRT(COEFP/NORM3(J2))
89    COEFP = COEFP * ISIGN * FACTOR
      GO TO 10
8     COEFP = ITAB2(J2,J1)
      IF(COEFP) 66,10,90
66    COEFP = -DSQRT(-COEFP/NORM2(J2))
      GO TO 91
90    COEFP = DSQRT(COEFP/NORM2(J2))
91    COEFP = COEFP * ISIGN * FACTOR
      GO TO 10
9     COEFP = ITAB1(J2,J1)
      IF(COEFP) 67,10,92
67    COEFP = -DSQRT(-COEFP/NORM1(J2))
      GO TO 93
92    COEFP = DSQRT(COEFP/NORM1(J2))
93    COEFP = COEFP * ISIGN * FACTOR
      GO TO 10
*
*     AN UNALLOWED STATE OR AN UNALLOWED PARENT
*
11    WRITE(IWRITE,1111)
 1111 FORMAT(' ERROR IN SUBROUTINE CFPD - THE STATE OR IT''S PARENT IS N
     :OT ALLOWED')
      CALL EXIT
10    CONTINUE
      RETURN
      END
*
*     ------------------------------------------------------------------
*       C F P F
*     ------------------------------------------------------------------
*
      SUBROUTINE CFPF(N,IVI,ILI,ISI,IVJ,ILJ,ISJ,COEFP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
*     THIS IS A DUMMY SUBROUTINE TO CALCULATE CFP OF F-ELECTRONS.  IT IS
*     VALID ONLY FOR ONE OR TWO ELECTRONS IN THE F-SHELL UNDER
*     CONSIDERATION.
*
      COEFP=1.D0
      RETURN
      END
*
*     ------------------------------------------------------------------
*       C F P P
*     ------------------------------------------------------------------
*
      SUBROUTINE CFPP(N,LI,ISI,LJ,ISJ,COEFP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
*     THIS SUBROUTINE EVALUATES THE COEFFICIENTS OF FRACTIONAL PARENTAGE
*     FOR EQUIVALENT P SHELL ELECTRONS FROM TABLES GIVEN IN J.C.SLATER
*     QUANTUM THEORY OF ATOMIC STRUCTURE,VOLUME2,P350(1960)
*     IN THE SUBROUTINE LIST N,THE NO. OF ELECTRONS,L THE ANGULAR
*     MOMENTUM QUANTUM NO.,(2S+1) THE SPIN QUANTUM NO. OF BOTH THE STATE
*     IN QUESTION AND ITS PARENT STATE ARE INPUT PARAMETERS.THE RESULT
*     IS OUTPUT AS COEFP
*
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(8)
      DIMENSION IL(3,3),IS(3,3),ITAB1(3,1),ITAB2(3,3),NORM1(3),NORM2(3)
      DATA IS(1,2),IS(1,3)/0,0/
      DATA IL(1,2),IL(1,3)/0,0/
*
*
*     SET UP P SHELL PARAMETERS AND TABLES
*
      DATA IL(1,1),IL(2,1),IL(2,2),IL(2,3),IL(3,1),IL(3,2),IL(3,3)/1,1,2
     :     ,0,0,2,1/
      DATA IS(1,1),IS(2,1),IS(2,2),IS(2,3),IS(3,1),IS(3,2),IS(3,3)/2,3,1
     :     ,1,4,2,2/
      DATA ITAB1(1,1),ITAB1(2,1),ITAB1(3,1)/1,1,1/
      DATA ITAB2(1,1),ITAB2(1,2),ITAB2(1,3),ITAB2(2,1),ITAB2(2,2),ITAB2(
     :     2,3),ITAB2(3,1),ITAB2(3,2),ITAB2(3,3)/1,0,0,1,-1,0,-9,-5,4/
      DATA NORM1(1),NORM1(2),NORM1(3)/1,1,1/
      DATA NORM2(1),NORM2(2),NORM2(3)/1,2,18/
*
*     TEST IF N IS IN THE FIRST HALF OF SHELL
*
99    IF(N-4) 40,103,103
*
*     TEST IF STATE IN QUESTION IS ALLOWED
*     IF IT IS, IDENTIFY THE ROW OF THE TABLE BY J1
*
40    J = 0
101   J = J+1
      IF(J-4) 41,8,8
41    IF(IL(N,J)-LI) 101,42,101
42    IF(IS(N,J)-ISI) 101,43,101
43    J1 = J
*
*     TEST IF PARENT STATE IS ALLOWED
*     IF IT IS, IDENTIFY THE COLUMN OF THE TABLE BY J2
*
      IF(N-1) 44,70,44
70    IF(LJ) 8,71,8
71    IF(ISJ-1) 8,1,8
44    J = 0
102   J = J+1
      IF(J-4) 45,8,8
45    IF(IL(N-1,J)-LJ) 102,46,102
46    IF(IS(N-1,J)-ISJ) 102,47,102
47    J2 = J
      GO TO 100
*
*     SIMILAR SETTING OF J1 AND J2 IF N IS IN SECOND HALF OF SHELL
*
103   M =6-N
      IF(M) 72,73,72
73    IF(LI) 8,74,8
74    IF(ISI-1) 8,75,8
72    J = 0
104   J = J+1
      IF(J-4) 48,8,8
48    IF(IL(M,J)-LI) 104,49,104
49    IF(IS(M,J)-ISI) 104,50,104
50    J1 = J
75    J = 0
105   J = J+1
      IF(J-4) 51,8,8
51    IF(IL(M+1,J)-LJ) 105,52,105
52    IF(IS(M+1,J)-ISJ) 105,53,105
53    J2 = J
*
*
*     IDENTIFY THE F.P.C AS A UNIQUE ELEMENT OF ITABN(J1,J2)
*
100   GO TO (1,2,3,4,4,1),N
1     COEFP = 1.D0
      GO TO 10
2     COEFP = ITAB1(J1,J2)
      IF(COEFP) 54,10,31
54    COEFP = -DSQRT(-COEFP/NORM1(J1))
      GO TO 10
31    COEFP = DSQRT(COEFP/NORM1(J1))
      GO TO 10
3     COEFP = ITAB2(J1,J2)
      IF(COEFP) 55,10,32
55    COEFP = -DSQRT(-COEFP/NORM2(J1))
      GO TO 10
32    COEFP =DSQRT(COEFP/NORM2(J1))
      GO TO 10
*
*     USE RECURRENCE RELATION EQUATION (19) OF RACAH FOR SECOND HALF OF
*     SHELL
*
4     ISIGN = (-1)**((ISI+ISJ-5)/2+LI+LJ)
      FACTOR=DFLOAT((7-N)*ISJ*(2*LJ+1))/DFLOAT(N*ISI*(2*LI+1))
      IF(N-5) 56,5,8
56    COEFP = ITAB2(J2,J1)
      IF(COEFP) 57,10,33
57    COEFP = -DSQRT(-COEFP/NORM2(J2))
      GO TO 34
33    COEFP = DSQRT(COEFP/NORM2(J2))
34    COEFP = COEFP * ISIGN * DSQRT(FACTOR)
      IF(LJ-1) 35,10,35
35    COEFP = -COEFP
      GO TO 10
5     COEFP = ITAB1(J2,J1)
      IF(COEFP) 58,10,36
58    COEFP = -DSQRT(-COEFP/NORM1(J2))
      GO TO 37
36    COEFP = DSQRT(COEFP/NORM1(J2))
37    COEFP = COEFP * ISIGN * DSQRT(FACTOR)
      GO TO 10
*
*     AN UNALLOWED STATE OR AN UNALLOWED PARENT
*
8     WRITE(IWRITE,8888)
 8888 FORMAT(' ERROR IN SUBROUTINE CFPP - THE STATE OR IT''S PARENT IS N
     :OT ALLOWED')
      CALL EXIT
10    CONTINUE
      RETURN
      END
*
*     ------------------------------------------------------------------
*       N T A B 1
*     ------------------------------------------------------------------
*
      FUNCTION NTAB1(NELCTS,K)
      INTEGER IROW(0:9)
      DATA IROW/0,2,5,10,12,14,16,18,20,22/
*
*     THIS SUBROUTINE CALCULATES THE ROW OF NTAB CORRESPONDING TO THE
*     PARENTS WHICH MAY GIVE RISE TO THE TERM ASSOCIATED WITH SHELL
*     LAMBDA .  E.G. IF WE SEEK THE ROW OF NTAB CONTAINING THE PARENTS
*     OF ONE OF THE P**3 TERMS, THE ROW = VALUE OF NTAB1 IS THAT
*     CONTAINING THE P**2 TERMS
*
*     USE IS MADE OF THE FACT THAT THE LIST OF POSSIBLE PARENTS (SEE
*     WHITE - ATOMIC SPECTRA - APPENDIX)  IS SYMMETRICAL ABOUT THE
*     CONFIGURATION L**(2L+1)
*
*
* --- FOR ONE ELECTRON IN A TERM, THE PARENT IS ALWAYS A SINGLET S TERM
*
      IF (NELCTS .EQ. 1) THEN
         NTAB1 = 2
      ELSE
         NPAR = NELCTS - 1
         L = K-1
         LHALF = 2*L+1
         IF (NPAR .GT. LHALF) NPAR = 2*LHALF - NPAR
         NTAB1 = IROW(L) + NPAR
      END IF
      END
*
*     ------------------------------------------------------------------
*       M U M D A D
*     ------------------------------------------------------------------
*
      SUBROUTINE MUMDAD(II,IJ,IK,M,X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/MEDEFN/IHSH,NJ(10),LJ(10),NOSH(10,2),J1QN(19,3,2),IJFUL(10)
      COMMON/INTERM/J1B(10,3,2),J1T(3,2)
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,JSC0,JSC1,
     :JSC2,JSC3
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
*
*     NOTICE THE NAMES IN THE COMMON BLOCKS. SEE SETUP FOR DESCRIPTION
*
* --- CALLS AND EVALUATES FRACTIONAL PARENTAGE COEFFICIENTS
*
   10 FORMAT(8H COEFP =,F15.9)
      X=1.0D0
      LIJ=LJ(IJ)
      IF(LIJ) 12,12,11
   12 IF(M)4,5,4
   11 N=NOSH(IJ,II)
      IVI=J1QN(IJ,1,II)
      ILI=(J1QN(IJ,2,II)-1)/2
      ISI=J1QN(IJ,3,II)
*
*     IF M=0 THERE ARE QUANTUM NUMBERS WITH TILDES TO CONSIDER
*
      IF(M) 1,2,1
    1 IVJ=J1B(IJ,1,II)
      ILJ=(J1B(IJ,2,II)-1)/2
      ISJ= J1B(IJ,3,II)
      GO TO 3
    2 IVJ=J1T(1,II)
      ILJ=(J1T(2,II)-1)/2
      ISJ=J1T(3,II)
    3 CALL CFP(LIJ,N,IVI,ILI,ISI,IVJ,ILJ,ISJ,COEFP)
      IF(IBUG2.GT.0) WRITE(IWRITE,10) COEFP
      X=X*COEFP
      IF(DABS(X).LT.1.D-14) GO TO 5
    4 LIJ=LJ(IK)
      IF(LIJ) 5,5,14
   14 IF(M) 6,7,6
    6 N=NOSH(IK,II)
      IVI=J1QN(IK,1,II)
      ILI=(J1QN(IK,2,II)-1)/2
      ISI=J1QN(IK,3,II)
      IVJ = J1B(IK,1,II)
      ILJ =(J1B(IK,2,II)-1)/2
      ISJ = J1B(IK,3,II)
      GO TO 8
    7 N=NOSH(IJ,II)-1
      IVI=IVJ
      ILI=ILJ
      ISI=ISJ
      IVJ=J1B(IJ,1,II)
      ILJ=(J1B(IJ,2,II)-1)/2
      ISJ = J1B(IJ,3,II)
    8 CALL CFP(LIJ,N,IVI,ILI,ISI,IVJ,ILJ,ISJ,COEFP)
      IF(IBUG2.GT.0) WRITE(IWRITE,10) COEFP
      X=X*COEFP
    5 CONTINUE
      RETURN
      END
*
*     ------------------------------------------------------------------
*       F I N D
*     ------------------------------------------------------------------
*
      CHARACTER*3 FUNCTION FIND (I,OF,EL)
*
*  ---  THIS ROUTINE FINDS ELECTRONS IN ONE OF THREE LISTS
*
      PARAMETER (NWD=30)
      CHARACTER*3 OF(NWD,2),EL(NWD)
*
      IF ( I .LE. (NWD)) THEN
         FIND = EL(I)
        ELSE IF ( I .LE. (NWD)*2 ) THEN
         FIND = OF(I-(NWD),1)
        ELSE
         FIND = OF(I-2*(NWD),2)
      END IF
      RETURN
      END
*
*     ------------------------------------------------------------------
*       G S T A T E
*     ------------------------------------------------------------------
*
      SUBROUTINE GSTATE(NFIRST,NLAST)
      PARAMETER (NWD=30,NCD=100,NWD2=2*NWD,NCD4=4*NCD)
*
      COMMON/INFORM/ IREADI,IWRITE,IOUT,IREADF,ISC(7)
      COMMON/STATES/NCFG,MAXORB,IAJCMP(NWD2),LJCOMP(NWD2),NJCOMP(NWD2),
     :NOCCSH(NCD4),NELCSH(5,NCD4),NOCORB(5,NCD4),J1QNRD(9,NCD4)
      COMMON/NOR/NCOM,NORBI,NORBF,IWAR
      CHARACTER*1 JCQN(9)
      DIMENSION J1QN(9),J2QN(9),J3QN(9)
      CHARACTER*8 LABEL(2)
      DATA LABEL/'INITIAL','FINAL'/
*
*      DATA DEFINING THE STATE IS READ IN AND PRINTED OUT.
*
    5 FORMAT(5(1X,A3,1H(,I2,1H)))
    6 FORMAT(9(1X,4X,I1,A1,I1))
   24 FORMAT(//31H INITIAL STATE CONFIGURATIONS:-)
   25 FORMAT(/5H     ,I3,3H.  ,10(1X,A3,1H(,I2,1H)))
   26 FORMAT(11X,10(1X,4X,I1,A1,I1))
   27 FORMAT(22X,9(1X,4X,I1,A1,I1))
   28 FORMAT(  31H ----------------------------  /)
   29 FORMAT(//29H FINAL STATE CONFIGURATIONS:-)
   30 FORMAT(2X,'ELECTRON ',A3,' NOT FOUND IN THE LIST OF ELECTRONS',
     :   ' FOR THE ',A8,' STATE')
      IF (NFIRST .EQ. 1) THEN
         WRITE(IWRITE,24)
         IREAD =IREADI
        ELSE
         WRITE(IWRITE,29)
         IREAD = IREADF
      END IF
      WRITE(IWRITE,28)
      DO 2 I=NFIRST,NLAST
      N=NOCCSH(I)
      READ(IREAD,5)        (NOCORB(J,I),NELCSH(J,I),J=1,N)
      K=I
      IF(NFIRST.NE.1) K=I-NFIRST+1
      WRITE(IWRITE,25) K,(NOCORB(J,I),NELCSH(J,I),J=1,N)
      NCOM1 = NCOM + 1
      NOR11 = NCOM1 + NORBI
      DO 61 J=1,N
      DO 63 JJ = 1,MAXORB
      IF (NFIRST .EQ. 1 .AND. JJ .GE. NOR11) GO TO 65
      IF(NFIRST .NE. 1 .AND. JJ .GE. NCOM1 .AND. JJ .LT. NOR11) GO TO 63
      IF(NOCORB(J,I).EQ.IAJCMP(JJ)) THEN
         NOCORB(J,I) = JJ
         GO TO 61
      END IF
   63    CONTINUE
*
*        ELECTRON NOT FOUND IN THE LIST
*
   65 WRITE(IWRITE,30) NOCORB(J,I),LABEL(NFIRST)
      STOP
   61 CONTINUE
      M=2*N-1
      N1=N+1
      READ(IREAD,6)    (J3QN(J),JCQN(J),J1QN(J),J=1,M)
      WRITE(IWRITE,26) (J3QN(J),JCQN(J),J1QN(J),J=1,N)
      IF(N.EQ.1) GO TO 64
      WRITE(IWRITE,27) (J3QN(J),JCQN(J),J1QN(J),J=N1,M)
   64 CONTINUE
      DO 62 J=1,M
      J2QN(J) = 2*LVAL(JCQN(J)) + 1
   62 J1QNRD(J,I)= (J3QN(J)*64 + J2QN(J))*64 + J1QN(J)
    2 CONTINUE
      RETURN
      END
*
*     ------------------------------------------------------------------
*       O R T H
*
*     ------------------------------------------------------------------
*
      SUBROUTINE ORTH
*
*       Determine the orthogonality between initial and final state
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NCD=100,NCD4=4*NCD)
      COMMON/STATES/NCFG,MAXORB,IAJCMP(2*NWD),LJCOMP(2*NWD),
     :NJCOMP(2*NWD),NOCCSH(NCD4),NELCSH(5,NCD4),NOCORB(5,NCD4),
     :J1QNRD(9,NCD4)
      COMMON/OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
     1 ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
     2 IORTH(NWD*(NWD-1)/2)

      COMMON/NOR/NCOM,NORBI,NORBF,IWAR
*
*   SET UP OF IORTH VECTOR
*   THE COMMON SET (NCOM) IS ASSUMED TO BE ORTHOGONAL TO BOTH
*   NORBI AND NORBF SETS
*
      IF (NORBI .EQ. 0) RETURN
      M1 = NORBI*NORBF
      DO 70 I = 1,M1
      IORTH(I) = 0
   70 CONTINUE
      NORTH = 0
*
*  | 1....NCOM | NCOM1.....NOR1 | NOR11.....NOR2|
*  |   NCOM    |     NORBI      |     NORBF     |
*  |   <= NWD  |     <= NWD     |     <=NWD     |
*
*  THIS LIMITATION IS LINKED TO THE DIMENSION OF BUFFER(NWD) IN
*  ANALYSE SUBROUTINE, where <= stands for LESS THAN or EQUAL.
*
      NCOM1 = NCOM+1
      NOR1  = NCOM + NORBI
      NOR11 = NOR1 + 1
      NOR2  = NOR1 + NORBF
      DO 78 J = NCOM1,NOR1
      DO 79 I = NOR11,NOR2
         IJ = NORBF*(J-NCOM1) + I - NOR1
         IF (LJCOMP(I) .EQ. LJCOMP(J)) THEN
            NORTH = NORTH + 1
            IORTH(IJ) = 1
         ENDIF
   79 CONTINUE
   78 CONTINUE
      RETURN
      END
*     ------------------------------------------------------------------
*       T E N S O R
*     ------------------------------------------------------------------
*
      SUBROUTINE TENSOR(KA,KB,ISPIN,IRHO,ISIG,VSHELL)
*
      IMPLICIT REAL *8(A-H,O-Z)
      PARAMETER(KFL1=60,KFL2=12)
*
*
*     W. D. ROBB   -   NOVEMBER 1971
*
*    Modified by C. FROESE FISCHER for use with NJGRAF
*
************************************************************************
*
*     A ROUTINE FOR THE EVALUATION OF ANGULAR AND SPIN FACTORS IN THE
*     REDUCED MATRIX ELEMENT OF ANY ONE-ELECTRON TENSOR OPERATOR BETWEEN
*     ARBITRARILY COUPLED L-S CONFIGURATIONS
*
*
*     **  NOTE THAT THE DEFINITIONS OF TENSOR OPERATORS USED ARE THOSE
*     OF FANO AND RACAH, IRREDUCIBLE TENSORIAL SETS, ACADEMIC PRESS 1959
*
************************************************************************
*
*                       DIMENSION STATEMENTS
*
      DIMENSION J2STO(KFL2,3),J3STO(KFL2,3),JMEM(119),VSHELL(20)
      LOGICAL FAIL,FREE
*
*                       COMMON BLOCKS
*
      COMMON/COUPLE/MN1,M0,J1(KFL1),J2(KFL2,3),J3(KFL2,3),FREE(KFL1)
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(7),IALL
      COMMON/MEDEFN/IHSH,NJ(10),LJ(10),NOSH1(10),NOSH2(10),J1QN1(19,3)
     :    ,J1QN2(19,3),IJFUL(10)
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/TERMS/NROWS,ITAB(24),JTAB(24),NTAB(333)
*
  203 FORMAT(//7H NJ,LJ ,10(I6,I3))
  204 FORMAT(//6H NOSH ,10I4)
  205 FORMAT(//6H J1QN ,30I3)
  207 FORMAT(8F15.8)
  208 FORMAT(// 23H PARENT TERMS NOT FOUND//)
  209 FORMAT(//3H J1)
  210 FORMAT(24I5)
  211 FORMAT(24H J2                   J3)
  212 FORMAT(3I5,I10,2I5)
  213 FORMAT(///26H ORBITAL RECOUPLING COEFF=,D20.8)
  214 FORMAT(///23H SPIN RECOUPLING COEFF=,D20.8//)
  215 FORMAT(/28H THE CONTRIBUTION FROM SHELL,I2,3H IS,F15.8)
  216 FORMAT(//21H THIS IS NOT A PARENT)
  217 FORMAT(///8H VSHELL=,8F15.8)
  218 FORMAT(//24H FRACTIONAL PARENT TERMS,I2)
  219 FORMAT(//49H THE CONTRIBUTION FROM FRACTIONAL PARENTAGE TERMS,I2,
     : 3H IS,F15.8)
  220 FORMAT(//6H SHELL,I2)
  302 FORMAT(/5X,89H NO CONTRIBUTION FROM TENSOR SINCE MORE THAN ONE ELE
     :CTRON DIFFERENT IN THE CONFIGURATIONS/)
  303 FORMAT(/5X,114H NO CONTRIBUTION FROM TENSOR SINCE THE TRIANGLE REL
     :ATION BETWEEN KA AND THE TOTAL ANGULAR MOMENTA IS NOT SATISFIED/)
  313 FORMAT(19H SPECTATOR SUBSHELL,I3,69H HAS DIFFERENT QUANTUM NUMBERS
     : ON THE TWO SIDES OF THE MATRIX ELEMENT/)
*
      AJF=1.D0
      RML=0.D0
      RPL=0.D0
      NTOT=0
      DO 100 IS=1,IHSH
      VSHELL(IS)=0.D0
  100 CONTINUE
      IHSHP1=IHSH+1
      I2HSH=IHSH*2-1
*
*     PRINT OUT THE OCCUPATION AND COUPLING ARRAYS
*
      IF(NBUG6-1) 101,2,101
    2 WRITE(IWRITE,203) (NJ(I),LJ(I),I=1,IHSH)
      WRITE(IWRITE,204)(NOSH1(J),J=1,IHSH)
      WRITE(IWRITE,204)(NOSH2(J),J=1,IHSH)
      WRITE(IWRITE,205) ((J1QN1(J,K),K=1,3),J=1,I2HSH)
      WRITE(IWRITE,205) ((J1QN2(J,K),K=1,3),J=1,I2HSH)
*
*      TEST FOR AT MOST ONE ELECTRON DIFFERENCE IN CONFIGURATIONS
*
  101 NOSHUM=0
      DO 102 K=1,IHSH
      NOSHUM=NOSHUM+IABS(NOSH1(K)-NOSH2(K))
  102 CONTINUE
      IF(NOSHUM.GT.2) GO TO 300
*
*     TEST FOR TRIANGLE RELATION BETWEEN KA AND TOTAL ANGULAR MOMENTA
*
  103 IF(ISPIN.EQ.0) GO TO 198
      K=3
      KC=KB
      IF(ISPIN.EQ.2) GO TO 199
      IF(J1QN1(I2HSH,2).NE.J1QN2(I2HSH,2)) GO TO 183
      GO TO 199
  198 K=2
      KC=KA
      IF(J1QN1(I2HSH,3).NE.J1QN2(I2HSH,3)) GO TO 183
  199 LB=J1QN1(I2HSH,K)-1
      NB=J1QN2(I2HSH,K)-1
      MB=KC+KC
      BTST=TRITST(MB,LB,NB)
      IF(DABS(BTST).GT.1.D-14) GO TO 301
      IF(K.EQ.2.OR.ISPIN.LT.2) GO TO 104
      K=2
      KC=KA
      GO TO 199
*
*      DETERMINE IRHO AND ISIGMA, THE NUMBERS OF THE OCCUPIED SHELLS
*
  104 IRHO=0
      ISIG=0
      DO 105 J=1,IHSH
      NX=NOSH1(J)-NOSH2(J)+2
      GO TO (107,105,106),NX
  107 ISIG=J
      GO TO 105
  106 IRHO=J
  105 CONTINUE
      IF(IRHO.NE.0 ) GO TO 108
      IRHO=1
      ISIG=1
  108 MEMR = IRHO
*
*     THE BEGINNING OF THE LOOP OVER ALL SHELLS
*
  109 IF(IRHO.NE.ISIG) GO TO 309
      IF(NBUG6-1) 309,4,309
    4 WRITE(IWRITE,220) IRHO
  309 NTOT=NTOT+1
      L1=LJ(IRHO)+1
      L2=LJ(ISIG)+1
      AJF=DFLOAT(J1QN1(I2HSH,2))/DFLOAT(2*LJ(IRHO)+1)
      IF(ISPIN.EQ.1) AJF=J1QN1(I2HSH,3)/2.D0
      IF(ISPIN.EQ.2) AJF=AJF*J1QN1(I2HSH,3)/2.D0
*
*     CHECK THE DIAGONAL CHARACTER OF QUANTUM NUMBERS OF SPECTATOR
*     SHELLS
*
      DO 255 J=1,IHSH
      IF(J.EQ.IRHO.OR.J.EQ.ISIG) GO TO 255
      DO 256 KK=1,3
      IF(J1QN1(J,KK).NE.J1QN2(J,KK)) GO TO 257
  256 CONTINUE
  255 CONTINUE
      GO TO 258
  257 IF(NBUG7.EQ.1) WRITE(IWRITE,313) J
      IF(IRHO.NE.ISIG) GO TO 190
      GO TO 189
  258 IF(IRHO-ISIG) 120,111,120
*
*     FIND THE PARENT TERMS GIVEN BY ALLOWED J VALUES IN NTAB WITH IRHO
*
  111 NELCTS=NOSH1(IRHO)
      K1=NTAB1(NELCTS,L1)
      KK1=ITAB(K1)
      DO 112 JJ1=1,KK1
      IJK1=3*(JJ1-1)+JTAB(K1)
      DO 113 K=2,3
      IJKK=IJK1+K
      IF(K.EQ.3) GO TO 114
      LA=NTAB(IJKK)
      MA=2*LJ(IRHO)+1
      NA=J1QN1(IRHO,K)
      GO TO 115
  114 LA=NTAB(IJKK)-1
      MA=1
      NA=J1QN1(IRHO,K)-1
  115 ATST=TRITST(LA,MA,NA)
      IF(DABS(ATST).GT.1.D-14) GO TO 116
  117 IF(K-3) 113,118,113
  116 JMEM(JJ1)=0
      GO TO 112
  118 JMEM(JJ1)=1
  113 CONTINUE
  112 CONTINUE
*
*     PARENTAGE CHECK
*
  120 IF(IRHO-ISIG) 121,127,121
  121 NELCTS=NOSH1(IRHO)
      K1=NTAB1(NELCTS,L1)
      NELCTS=NOSH2(ISIG)
      K2=NTAB1(NELCTS,L2)
      KK1=ITAB(K1)
      KK2=ITAB(K2)
      DO 122 JJ1=1,KK1
      IJK1=3*(JJ1-1)+JTAB(K1)
      DO 123 K=2,3
      IJKK=IJK1+K
      MSAM1=NTAB(IJKK)-J1QN2(IRHO,K)
      IF(MSAM1.NE.0) GO TO 122
      IF(K.EQ.3) GO TO 124
  123 CONTINUE
  122 CONTINUE
      IF(NBUG6-1) 192,7,192
    7 WRITE(IWRITE,208)
      GO TO 192
  124 DO 125 JJ2=1,KK2
      IJK2=3*(JJ2-1)+JTAB(K2)
      DO 126 K=2,3
      IJKK=IJK2+K
      MSAM2=NTAB(IJKK)-J1QN1(ISIG,K)
      IF(MSAM2.NE.0) GO TO 125
      IF(K.EQ.3) GO TO 127
  126 CONTINUE
  125 CONTINUE
      IF(NBUG6-1) 192,8,192
    8 WRITE(IWRITE,208)
      GO TO 192
*
*     SET  J2  AND  J3 .  SAME FOR  L  AND  S
*
  127 M1=IHSH-2
      M2=2*M1+1
      M3=3*IHSH-1
      M4=M3+1
      M5=M3+2
      M10=M5+1
      MN1=M10+1
      J2(1,1)=M10
      J2(1,2)=MN1
      J2(1,3)=M5
      J2(2,1)=IRHO
      J2(2,2)=M5
      J2(2,3)=M3
      J3(1,1)=ISIG
      J3(1,2)=M10
      J3(1,3)=M4
      IF(IRHO-1) 128,129,128
  129 J2(3,1)=M3
      GO TO 130
  128 J2(3,1)=1
  130 IF(IRHO-2) 131,132,131
  132 J2(3,2)=M3
      GO TO 133
  131 J2(3,2)=2
  133 J2(3,3)=IHSHP1
      IF(ISIG-1) 134,135,134
  135 J3(2,1)=M4
      GO TO 136
  134 J3(2,1) = 1
  136 IF(ISIG-2) 137,138,137
  138 J3(2,2)=M4
      GO TO 139
  137 J3(2,2)=2
  139 J3(2,3)=2*IHSH
      IF(IHSH-3) 149,140,140
  140 DO 148 J=4,IHSHP1
      L=J-1
      J2(J,1)=M1+L
      J2(J,3)=M1+J
      J3(L,1)=M2+L
      J3(L,3)=M2+J
  141 IF(IRHO-L) 142,143,142
  143 J2(J,2)=M3
      GO TO 144
  142 J2(J,2)=L
  144 IF(ISIG-L) 145,146,145
  146 J3(L,2)=M4
       GO TO 148
  145 J3(L,2)=L
  148 CONTINUE
  149 M6=IHSHP1
      J3(M6,1)=M3-1
      J3(M6,2)=MN1
      J3(M6,3)=I2HSH
      IF(IHSH-1) 450,451,450
  451 J3(M6,1) = M4
      J3(M6,3) = M3
  450 DO 150 J=1,IHSHP1
      DO 151 K=1,3
      J2STO(J,K)=J2(J,K)
      J3STO(J,K)=J3(J,K)
  151 CONTINUE
  150 CONTINUE
*
*     RECOUPLING COEFFICIENTS
*
      JMEM1=J1QN1(IRHO,1)
      JMEM2=J1QN1(IRHO,2)
      JMEM3=J1QN1(IRHO,3)
      JMEM4=J1QN2(ISIG,1)
      JMEM5=J1QN2(ISIG,2)
      JMEM6=J1QN2(ISIG,3)
      IF(IRHO-ISIG) 154,152,154
*
*     BEGINNING OF LOOP OVER ALL PARENT TERMS
*
  152 JJ1=1
 1152 IF(NBUG6-1) 12,11,12
   11 WRITE(IWRITE,218) JJ1
   12 IF(JMEM(JJ1).EQ.1) GO TO 153
      IF(NBUG6-1) 186,16,186
   16 WRITE(IWRITE,216)
      GO TO 186
  153 IJK1=3*(JJ1-1)+JTAB(K1)
      NI1=NTAB(IJK1+1)
      NI2=NTAB(IJK1+2)
      NI3=NTAB(IJK1+3)
      J1QN2(IRHO,1)=NI1
      J1QN1(ISIG,1)=NI1
      J1QN2(IRHO,2)=NI2
      J1QN1(ISIG,2)=NI2
      J1QN2(IRHO,3)=NI3
      J1QN1(ISIG,3)=NI3
  154 K=2
      M7=M3-IHSH
      M9=M7+1
      M11=M3-1
      M12=IHSH-1
      RECUPS=1.D0
      M0=M6+1
*
*     SET UP THE J1 ARRAY FOR THE ANGULAR AND SPIN RECOUPLING
*     COEFFICIENTS
*
  155 IF(K-3) 156,157,157
  156 J1(M5)=2*LJ(IRHO)+1
      J1(M10)=2*LJ(ISIG)+1
      J1(MN1)=2*KA+1
      IF(ISPIN.EQ.1) J1(MN1)=1
      J1(M3)=JMEM2
      J1(M4)=JMEM5
      IF(IRHO.EQ.ISIG) GO TO 158
      J1(M3)=J1QN1(IRHO,K)
      J1(M4)=J1QN2(ISIG,K)
      GO TO 158
  157 J1(M5)=2
      J1(M10)=2
      J1(MN1)=KB+KB+1
      IF(ISPIN.EQ.0) J1(MN1)=1
      J1(M3)=JMEM3
      J1(M4)=JMEM6
      IF(IRHO.EQ.ISIG) GO TO 158
      J1(M3)=J1QN1(IRHO,K)
      J1(M4)=J1QN2(ISIG,K)
  158 DO 161 J=1,IHSH
      IF(IRHO-J) 160,159,160
  159 J1(J)=J1QN2(IRHO,K )
      GO TO 161
  160 J1(J)=J1QN1(J,K)
  161 CONTINUE
      IF(IHSH.EQ.1) GO TO 197
      DO 162 J=M6,M7
      J1(J)=J1QN1(J,K)
  162 CONTINUE
      DO 163 J=M9,M11
      JM12=J-M12
      J1(J)=J1QN2(JM12,K)
  163 CONTINUE
*
*     PRINT OUT THE J1,J2 AND J3 ARRAYS
*
  197 IF(NBUG6-1) 304,9,304
    9 IF(K-3) 165,164,164
  165 IF(NBUG6-1) 304,17,304
   17 WRITE(IWRITE,209)
      WRITE(IWRITE,210) (J1(J),J=1,MN1)
      WRITE(IWRITE,211)
      DO 166 I=1,IHSHP1
      WRITE(IWRITE,212) (J2(I,J),J=1,3),(J3(I,J),J=1,3)
  166 CONTINUE
  304 CONTINUE
*
*     EVALUATE ORBITAL AND SPIN RECOUPLING COEFFICIENTS
*
  164 DO 500 I = 1,MN1
         FREE(I) = .FALSE.
  500 CONTINUE
*
      CALL NJGRAF(RECUP,FAIL)
*
      RECUPS=RECUPS*RECUP
      IF(K-3) 167,170,170
  167 IF(NBUG6-1) 305,18,305
   18 WRITE(IWRITE,213) RECUP
  305 CONTINUE
  170 K=K+1
      DO 168 J=1,IHSHP1
      DO 169 KK=1,3
      J2(J,KK)=J2STO(J,KK)
      J3(J,KK)=J3STO(J,KK)
  169 CONTINUE
  168 CONTINUE
      IF(K.EQ.3) GO TO 155
      IF(NBUG6-1) 306,19,306
   19 WRITE(IWRITE,214) RECUP
*
*     FIRST FRACTIONAL PARENTAGE COEFFICIENT
*
  306 LIJ=LJ(IRHO)
      COEFP=1.D0
      IF(LIJ) 171,272,171
  171 N=NOSH1(IRHO)
      IV1=JMEM1
      IL1=(JMEM2-1)/2
      IS1= JMEM3
      IV2=J1QN2(IRHO,1)
      IL2=(J1QN2(IRHO,2)-1 )/2
      IS2=J1QN2(IRHO,3)
      CALL CFP(LIJ,N,IV1,IL1,IS1,IV2,IL2,IS2,COEFP)
      RECUPS=RECUPS*COEFP
  272 IF(IRHO-ISIG) 172,173,172
  172 IF(DABS(RECUPS).LT.1.D-14) GO TO 183
*
*     SECOND FRACTIONAL PARENTAGE COEFFICIENT
*
  173 LIJ=LJ(ISIG)
      COEFP=1.D0
      IF(LIJ) 176,176,174
  174 N=NOSH2(ISIG)
      IV1=JMEM4
      IL1=(JMEM5-1)/2
      IS1=JMEM6
      IV2=J1QN1(ISIG,1)
      IL2=(J1QN1(ISIG,2)-1)/2
      IS2=J1QN1(ISIG,3)
      CALL CFP(LIJ,N,IV1,IL1,IS1,IV2,IL2,IS2,COEFP)
  176 RECUPS=RECUPS*COEFP
      IF(DABS(RECUPS).LT.1.D-14.AND.IRHO.NE.ISIG) GO TO 183
*
*     PERMUTATION FACTOR
*
  175 IDELP=2
      IF(IRHO-ISIG) 177,181,179
  177 JRHO = IRHO+1
      DO 178 J=JRHO,ISIG
  178 IDELP=IDELP+NOSH1(J)
      GO TO 181
  179 JSIG = ISIG+1
      DO 180 J=JSIG,IRHO
  180 IDELP = IDELP+NOSH2(J)
  181 MINUS=(-1)**IDELP
*
*     MULTIPLICATIVE FACTOR
*
      IF(IRHO-ISIG) 182,185,182
  182 SQRN=DSQRT(DFLOAT(NOSH1(IRHO)*NOSH2(ISIG)))
      VALML=SQRN*RECUPS*DFLOAT(MINUS)
      GO TO 184
  183 VALML=0.D0
  184 RML = RML+VALML
*     RESULT STORED IN VSHELL
      IF(NTOT.EQ.0) NTOT=1
      VSHELL(NTOT)=RML*DSQRT(AJF)
      GO TO 190
  185 VALUML=RECUPS
      IF(NBUG6.NE.0) WRITE(IWRITE,219) JJ1,VALUML
      RPL = RPL+VALUML
  186 IF(IRHO.NE.ISIG)GO TO 1186
      JJ1=JJ1+1
      IF(JJ1.LE.KK1)GO TO 1152
 1186 J1QN1(IRHO,1)=JMEM1
      J1QN1(IRHO,2)=JMEM2
      J1QN1(IRHO,3)=JMEM3
      J1QN2(ISIG,1)=JMEM4
      J1QN2(ISIG,2)=JMEM5
      J1QN2(ISIG,3)=JMEM6
      ANL=DFLOAT(NOSH1(IRHO))*RPL
*
*      RESULTS STORED IN VSHELL
*
      IF(NTOT.EQ.0) NTOT=1
      VSHELL(NTOT)=ANL*DSQRT(AJF)
  194 IF(NBUG6-1) 189,196,189
  196 WRITE(IWRITE,215) IRHO,ANL
  189 IRHO=IRHO+1
      ISIG=ISIG+1
      RPL=0.D0
      IF(IRHO-IHSH) 109,109,190
  190 IF(NBUG6-1) 192,13,192
   13 WRITE(IWRITE,217) (VSHELL(N),N=1,NTOT)
  192 RETURN
  300 IF(NBUG6.NE.0) WRITE(IWRITE,302)
      RETURN
  301 IF(NBUG6.NE.0) WRITE(IWRITE,303)
      RETURN
      END
*-----------------------------------------------------------------------
*        Q S O R T
*-----------------------------------------------------------------------

*     The method use to sort the data is quick sort with a pivot value
*     that is the larger value of the first 2 different value from the
*     the sublist to be sorted.
*     This sorting method used a stack to maintain the unsorted section,
*     and sorting will be finished when the stack is empty.
      subroutine qsort(n,key,pt,stack,nstack,ierr)
      integer top, from, to, pivot, left, right
      integer key(*),pt(*),stack(*)
*
*     Set the initial pointer values
*
      do 10 i=1,n
10    pt(i)=i
*
*     Initialize the stack and error indicator
*
      top=1
      stack(top)=1
      ierr = 0
*
*     Repeat Until the Stack is empty
*
100   continue
*
*     Determine the next section
*
      if (top .ne. 0) then
        from=stack(top)
        if (top.ne.1) then
          to=stack(top-1) - 1
        else
          to=n
        endif
        top=top-1
*
*        Find the position k of the pivot value that partitions
*        the current section. Return a value of k=0 when there is
*        no disinct value.
*
        if (from .eq. to) then
          k=0
        else
          k=0
          ismall=key( pt(from) )
          do 210 i=from+1,to
            if (key( pt(i) ) .ne. ismall) then
              k=i
              goto 200
            endif
210       continue
        endif
200     continue
        if (k.ne.0) then
          if( ismall .gt. key(pt(k)) ) then
            k=from
          endif
        endif
        if (k .ne. 0)then
*
*         Rearange the section of the keys such that all values
*         smaller than the pivot value are stored on the left and
*         larger values on the right.
*
          pivot=key(pt(k))
          left=from
          right=to
300       continue
310       if ( key(pt(left)) .lt. pivot ) then
            left = left+1
            goto 310
          endif
320       if ( key(pt(right)) .ge. pivot ) then
            right = right-1
            goto 320
          endif
          if (left .lt. right) then
            i=pt(left)
            pt(left)=pt(right)
            pt(right)=i
            goto 300
          endif
          kp=left
          if (top+2 .le. nstack) then
            stack(top+1)=kp
            stack(top+2)=from
            top=top+2
          else
            ierr =1
            return
          endif
        endif
        goto 100
      endif
*
*     Keys are sorted
*
      return
      end
*-----------------------------------------------------------------------
*        S A V E
*-----------------------------------------------------------------------
*
      SUBROUTINE SAVE(ICASE,C,K,I1,I2,I3,I4,JA,JB,IPTR)
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /FOUT/NOV(2),IOVLAP(10,2),NCOUNT(8),IFLAG,NIJ
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(8)

      IF (ICASE .LE. 2 .or. ICASE .EQ. 4 .or. ICASE .EQ. 5) THEN
*
*     Fk, Gk, L, or Z data
*
         IF (I2 .GT. I4) THEN
            II2 = I4
            II4 = I2
         ELSE
            II2 = I2
            II4 = I4
         END IF
         IPACK = (K*64 + II2)*64 + II4
         IF (ICASE .NE. 4) THEN
           WRITE(ISC(ICASE)) C,IPACK,JA,JB
         ELSE
           WRITE(ISC(ICASE)) C,IPACK,JA,JB,IPTR
         END IF
      ELSE IF(ICASE .EQ. 3) THEN
         J = 1
         IMIN = I1
         IF (I2 .LT. IMIN) THEN
            IMIN=I2
            J = 2
         END IF
         IF (I3 .LT. IMIN) THEN
            IMIN = I3
            J = 3
         END IF
         IF (I4 .LT. IMIN) THEN
            IMIN = I4
            J = 4
         END IF
         GO TO (10,20,30,40) J
10       II1 = I1
         II2 = I2
         II3 = I3
         II4 = I4
         Go to 50
        
20       II1 = I2
         II2 = I1
         II3 = I4
         II4 = I3
         GO TO 50

30       II1 = I3
         II2 = I4
         II3 = I1
         II4 = I2
         GO TO 50

40       II1 = I4
         II2 = I3
         II3 = I2
         II4 = I1

50       IPACK = (((K*64+II1)*64+II2)*64+II3)*64+II4
         WRITE(ISC(3)) C,IPACK,JA,JB,IPTR
      ELSE 
         II1 = I1
         II3 = I3
         IF (I2 .GT. I4) THEN
            II2 = I4
            II4 = I2
         ELSE
            II2 = I2
            II4 = I4
         END IF
         IF (ICASE .NE. 7) THEN
            IF (I1 .GT. I3) THEN
               II1 = I3
               II3 = I1
            END IF
         END IF
*
* ... Because the k-value may be -1, for these integrals, a
*     value of k+1 is stored.
*
         KK = K + 1
         IPACK = (((KK*64+II1)*64+II2)*64+II3)*64+II4
         WRITE(ISC(ICASE)) C,IPACK,JA,JB
      END IF
      NCOUNT(ICASE) = NCOUNT(ICASE) + 1
      IFLAG = 1
      END
*
*     ------------------------------------------------------------------
*
*       T R I T S T
*     ------------------------------------------------------------------
*
      DOUBLE PRECISION FUNCTION TRITST(L,M,N)
*
*
*      IF  TRITST=1.0   THE TRIANGLE RELATION IS NOT SATISFIED
*      IF  TRITST=0.0   THE TRIANGLE RELATION IS SATISFIED
*
      LMN=IABS(L-M)
      LM=L+M
      IF(N-LMN) 1,2,2
    2 IF(LM-N) 1,3,3
    3 TRITST=0.D0
      RETURN
    1 TRITST=1.D0
      RETURN
      END
*
*     ------------------------------------------------------------------
*       V I J O U T
*     ------------------------------------------------------------------
*
      SUBROUTINE VIJOUT(JA,JB)
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,JSC0,JSC1,
     :JSC2,JSC3
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/MEDEFN/IHSH,NJ(10),LJ(10),NOSH1(10),NOSH2(10),J1QN1(19,3),
     :     J1QN2(19,3),IJFUL(10)
*
*     THIS SUBROUTINE IS ENTERED ONLY IF IBUG2 IS GREATER THAN ZERO
*
* --- PRINT OUT OF QUANTUM NUMBERS AND COUPLING SCHEMES FOR EACH
*     MATRIX ELEMENT AS DEFINED BY SETUP
*
    5 FORMAT(//48H L.H.S. OF HAMILTONIAN MATRIX ELEMENT DEFINED BY)
    6 FORMAT(//48H R.H.S. OF HAMILTONIAN MATRIX ELEMENT DEFINED BY)
    7 FORMAT(9H1(CONFIG ,I2,10H/V/CONFIG ,I2,1H))
    8 FORMAT(/7H NJ,LJ ,10(I6,I3))
    9 FORMAT(/6H NOSH ,10I4)
   10 FORMAT(6H J1QN ,10(I5,2I3))
      I2HSH=2*IHSH-1
      WRITE(IWRITE,7) JA,JB
      WRITE(IWRITE,8) (NJ(I),LJ(I),I=1,IHSH)
      WRITE(IWRITE,5)
      WRITE(IWRITE,9) (NOSH1(J),J=1,IHSH)
      WRITE(IWRITE,10) ((J1QN1(J,K),K=1,3),J=1,I2HSH)
      WRITE(IWRITE,6)
      WRITE(IWRITE,9) (NOSH2(J),J=1,IHSH)
      WRITE(IWRITE,10) ((J1QN2(J,K),K=1,3),J=1,I2HSH)
    1 RETURN
      END
