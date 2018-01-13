*	Program to print wavefunctions
*	
*	Created by C. Froese Fischer June 16, 1987
*	Vanderbilt University

      PROGRAM PRINT
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER AT*6,TT*6,EL1*3,NEW*3,NAME*24
      DIMENSION PT(220)
*
*     iarg = iargc()
*     if (iarg .gt. 0) then
*	 call getarg(1,NAME)
*      else 
         NAME = 'wfn.inp'
*     end if
      OPEN(UNIT=3,FILE=NAME,STATUS='OLD',FORM='UNFORMATTED')
      IUF=3
2     READ(IUF,END=5) AT,TT,EL1,M,ZT,ETI,EKI,AZI,(PT(J),J=1,M)
      WRITE(6,'(2x,A,A,$)') EL1,' = '
      READ(5,'(A)') NEW
      IF ( NEW .NE. 'd  ' .AND. NEW .NE. 'D  ' ) THEN
         IF ( NEW .NE. '   ') THEN
            EL1 = NEW
         ENDIF
         WRITE(6,'(5F14.6)') (PT(J),J=1,M)
      END IF
      GO TO 2
5     CLOSE(UNIT=3)
      CLOSE(UNIT=4)
      END
