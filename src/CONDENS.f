*     ------------------------------------------------------------------
*
*       A GENERAL PROGRAM TO CONDENSE THE LIST OF CONFIGURATIONS
*
*       by C. Froese Fischer
*          Vanderbilt University
*          Nashville, TN 37235 USA
*
*       March, 1984
*
*     Computer Physics Communication, Vol. 64, 399-405 (1991)
*     ------------------------------------------------------------------
*
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        CHARACTER S1*1,S4*4,S7*7,S19*19,S21*21,S23*23,S24*24
        CHARACTER*24 INFILE,OUTFILE,NAME
        CHARACTER*40 CONFIG(800)
        CHARACTER*72 COUPLE(800),HEADER,CLSDSH
        DIMENSION U(800),V(800),IDEL(100),INDX(5),IP(800)
        LOGICAL KEEP(800)
        DATA KEEP/800*.TRUE./
        DATA V/800*0.D0/
 
 
        WRITE(0,'(//A/A/A/A/A)') ' Selection Process:','    1  name.c',
     :        '    2  name.l','    3  name.j','    4  User Defined'
        WRITE(0,'(//A/A/A/A)') ' List to be condensed: ',
     :        '    1  int.lst','   2  name.i','    3  NONE'
        WRITE(0,'(//A)') ' Enter both selections'
        READ(5, *) ICASE,ILIST
CSUN    iarg = iargc()
CSUN    if ( iarg .eq. 0) then
           WRITE(0,*) 'Enter name of atom'
           read(5,'(A)') NAME
CSUN    else
CSUN       call getarg(1,name)
CSUN    end if
	j = index(name,' ')
	infile = name(1:j-1)//'.c'
 
        OPEN(UNIT=1,FILE=infile,STATUS='OLD')
 
        OUTFILE = 'cfg.out'
        OPEN(UNIT=3,FILE=OUTFILE,STATUS='UNKNOWN')
 
        READ(1,'(A72/A72)') HEADER, CLSDSH
        WRITE(3,'(A72/A72)') HEADER, CLSDSH
 
        I = 1
 10     READ(1,'(A40,F10.8)',END=20) CONFIG(I),U(I)
        IF (CONFIG(I)(1:1) .EQ. '*') GO TO 20
        IF (ICASE .EQ. 1) V(I) = ABS(U(I))
        READ(1,'(A72)') COUPLE(I)
        I = I+1
        IF (I .LE. (800)) GO TO 10
 
 20     N = I-1
        IF (ICASE .EQ. 2 .OR. ICASE .EQ. 3) THEN
              IF (ICASE .EQ. 2) THEN
              INFILE =NAME(1:j-1)//'.l'
           ELSE
              INFILE =NAME(1:j-1)//'.j'
           END IF
           OPEN(UNIT=2,FILE=INFILE,STATUS='OLD')
	   READ(2,'(A)') HEADER
25         READ(2,'(//22X,I4)',END=45) NUMBER
           DO 30 J = 1,NUMBER
              READ(2,'(//(7F10.6))') (U(I),I=1,N)
              DO 40 I = 1,N
                 V(I) = MAX(V(I),ABS(U(I)))
 40           CONTINUE
 30        CONTINUE
           GO TO 25
 45        CONTINUE
        ELSE IF (ICASE .EQ. 4) THEN
           DO 32 I = 1,N
              V(I) = 1.D0
 32        CONTINUE
           WRITE(0,'(A/A)') ' Enter the index of the configurations to'
	   
     :         ,' be deleted  (five at a time, terminate with a zero)'
           M = 1
 34        WRITE(0,'(I3,A2)') M,': '
           READ(5,*) INDX
           DO 35 K = 1,5
              I = INDX(K)
              IF (I .NE. 0) THEN
                 V(I) = 0.
              ELSE
                 GO TO 36
              END IF
 35        CONTINUE
           GO TO 34
        END IF
 
 36     TOL = 1.D-8
        IF (.NOT.(ICASE .EQ. 4)) THEN
           WRITE(0,'(//A)')' Tolerance for acceptance -- FORMAT(F): '
           READ(5,*) TOL
        END IF
	CALL SORT(N,V,IP)
        DO 70 II = 1,N
	   I = IP(II)
           IF (V(I) .GE. TOL) THEN
                K = 72
72              IF (COUPLE(I)(K:K) .EQ. ' ') THEN
                   K = K-1
                   GO TO 72
                END IF
		IF (ICASE .EQ. 1 .OR. ICASE .EQ. 4) THEN
                  WRITE(3,71) CONFIG(I),U(I),COUPLE(I)
		ELSE
                  WRITE(3,71) CONFIG(I),V(I),COUPLE(I)
		END IF
           ELSE
                KEEP(I) = .FALSE.
           END IF
 71        FORMAT(A40,F10.7/A72)
 70     CONTINUE
 
        IF (ICASE .EQ. 2 .OR. ICASE .EQ. 3)  CLOSE(UNIT=2)
        CLOSE(UNIT=3)
        IF (ILIST .EQ. 3) STOP
 
        NINT = 8
           OUTFILE = 'int.out'
        IF (ILIST .EQ. 1) THEN
           INFILE = 'int.lst'
        ELSE
	   j = index(name,' ')
	   infile = name(1:j-1)//'.i'
        END IF
 
        OPEN(UNIT=2,FILE=INFILE,STATUS='OLD')
        OPEN(UNIT=3,FILE=OUTFILE,STATUS='UNKNOWN')
 
        READ(2,'(A72)') HEADER
        WRITE(3,'(A72)') HEADER
 
  1     FORMAT(A21,I3,A4,I3,A1)
  2     FORMAT(A24,I3,A7,I3,A)
  3     FORMAT(A19,I3,A4,I3,A)
  4     FORMAT(A19,I3,A4,I3,A1)
  5     FORMAT(A24,I3,A7,I3,A1)
 
 
        DO 100 INT = 1,NINT
 
           IF (INT .LE. 2) THEN
201             READ(2,1) S21,IH,S4,JH,S1
                IF (IH .EQ. 0) GO TO 99
                IF (KEEP(IH).AND.KEEP(JH))
     :             WRITE(3,1) S21,IP(IH),S4,IP(JH),S1
                GO TO 201
           ELSE IF (INT .EQ. 3) THEN
202             READ(2,2) S24,IH,S7,JH,S23
                IF (IH .EQ. 0) GO TO 99
                IF (KEEP(IH).AND.KEEP(JH)) THEN
                   K = 23
 90                IF (S23(K:K) .EQ. ' ') THEN
                        K = K-1
                        GO TO 90
                   END IF
                   WRITE(3,92) S24,IP(IH),S7,IP(JH),(S23(J:J),J=1,K)
 92                FORMAT(A24,I3,A7,I3,23A1)
                END IF
                GO TO 202
           ELSE IF (INT .EQ. 4) THEN
203             READ(2,3) S19,IH,S4,JH,S23
                IF (IH .EQ. 0) GO TO 99
                IF (KEEP(IH).AND.KEEP(JH)) THEN
                   K = 23
 94                IF (S23(K:K) .EQ. ' ') THEN
                        K = K-1
                        GO TO 94
                   END IF
                   WRITE(3,96) S19,IP(IH),S4,IP(JH),(S23(J:J),J=1,K)
 96                FORMAT(A19,I3,A4,I3,23A1)
                END IF
                GO TO 203
           ELSE IF (INT .EQ. 5) THEN
204             READ(2,4) S19,IH,S4,JH,S1
                IF (IH .EQ. 0) GO TO 99
                IF (KEEP(IH).AND.KEEP(JH))
     :             WRITE(3,4) S19,IP(IH),S4,IP(JH),S1
                GO TO 204
           ELSE
205             READ(2,5) S24,IH,S7,JH,S1
                IF (IH .EQ. 0) GO TO 99
                IF (KEEP(IH).AND.KEEP(JH))
     :             WRITE(3,5)S24,IP(IH),S7,IP(JH),S1
                GO TO 205
           END IF
 99        WRITE(3,'(14X,A1)') '*'
100     CONTINUE
        CLOSE(UNIT=2)
        CLOSE(UNIT=3)
        END
*--------------------------------------------------------------
*       S O R T
*--------------------------------------------------------------
*
	SUBROUTINE SORT(N,V,IP)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION V(*),IP(*)

	DO 1 I = 1,N
	   IP(I) = I
  1     CONTINUE

	DO 10 I = 1,N-1
          JP = I
          DO 12 J = I+1,N
             IF (V(IP(J)) .GT. V(IP(JP))) JP = J
 12       CONTINUE
	  ITEMP = IP(I)
	  IP(I) = IP(JP)
	  IP(JP) = ITEMP
 10     CONTINUE
	END
 
