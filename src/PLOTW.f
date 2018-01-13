*	Program to list wavefunctions
*	
*	Created by C. Froese Fischer June 30, 1987
*	Vanderbilt University

      PROGRAM PLOTW
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER AT*6,TT*6,EL1*3,EL*3,NEW*3,INPUT*24,OUTPUT*24,
     :          ATOM*6,TERM*6
      DIMENSION P(220,30),EL(30),R(220),R2(220)
*
      iarg = iargc()
*     if (iarg .gt. 0) then
*	 call getarg(1,INPUT)
*     else
         INPUT = 'wfn.inp'
*     end if
*     if (iarg .gt. 1) then
*	 call getarg(1,OUTPUT)
*     else
         OUTPUT = 'plot.dat'
*     end if
      inc = 1
*     if (iarg .gt. 3) call getarg(3,inc)
      OPEN(UNIT=3,FILE=INPUT,STATUS='OLD',
     :   FORM='UNFORMATTED')
      OPEN(UNIT=4,FILE=OUTPUT,STATUS='UNKNOWN')
      IUF=3
      nwf=1
      MM = 0
2     P(1,nwf) = 0.D0
      READ(IUF,END=5) AT,TT,EL1,M,ZT,ETI,EKI,AZI,(P(J,NWF),J=2,M+1)
      WRITE(6,'($,2x,A,A)') EL1,' = '
      READ(5,'(A)') NEW
      IF ( NEW .NE. 'd  ' .AND. NEW .NE. 'D  ' ) THEN
         IF ( NEW .NE. '   ') THEN
            EL1 = NEW
         ENDIF
         EL(NWF) = EL1
         NWF = NWF + 1
         Z = ZT
         ATOM = AT
         TERM = TT
         MM = MAX(M,MM)
      END IF
      GO TO 2
5     CLOSE(UNIT=3)
      RHO = -4.0
      H = 1/16.
      R(1) = 0.d0
      R2(1) = 0.d0
      DO 10 J = 2,220
        R(J) = EXP(RHO)/Z
        R2(J) = SQRT(R(J))
        RHO = RHO + H
10    CONTINUE
*
*	LIST TABLES OF RADIAL FUNCTIONS
*
      NWF = NWF -1
      write(4,*) ' sqrt(r)    P(nl;r)'
      DO 20 I = 1,NWF
         write(4,*) EL(I)
         DO 21 J = 1,MM,inc
            Pwave = P(J,I)*R2(J)
            IF (abs(Pwave) .gt. 0.0005 .OR. J.EQ.1)
     :         WRITE(4,'(F10.4,F12.3)') R2(J),Pwave
21       CONTINUE
20    CONTINUE
      END
