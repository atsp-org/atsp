*-----------------------------------------------------------------------
*
*	RELABEL  -- A PROGRAM TO RELABEL ORBITAL LABELS AND DELELTE
*                   ORBITALS
*
*	by C. Froese Fischer
*	   Vanderbilt University
*	   Nashville, TN 37235 USA
*
*       May, 1983
*
*     Computer Physics Communication, Vol. 64, 399-405 (1991)
*-----------------------------------------------------------------------
*
      PROGRAM RELABEL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220)
      CHARACTER AT*6,TT*6,EL1*3,NEW*3, NAME*24
      DIMENSION PT(220)
*
      NAME =  'wfn.inp'
*
*    Opens input files
*
      OPEN(UNIT=3,FILE=NAME,STATUS='OLD',
     :   FORM='UNFORMATTED')
      OPEN(UNIT=4,FILE='wfn.out',STATUS='UNKNOWN',
     :   FORM='UNFORMATTED')
      IUF=3
2     READ(IUF,END=5) AT,TT,EL1,M,ZT,ETI,EKI,AZI,(PT(J),J=1,M)
      WRITE(6,'(2X,A,A)') EL1,' = '
      READ(5,'(A)') NEW
      IF ( NEW .NE. 'd  ' .AND. NEW .NE. 'D  ' ) THEN
         IF ( NEW .NE. '   ') THEN
            EL1 = NEW
         ENDIF
         WRITE(4) AT,TT,EL1,M,ZT,ETI,EKI,AZI,(PT(J),J=1,M)
      END IF
      GO TO 2
5     CLOSE(UNIT=3)
      CLOSE(UNIT=4)
      END
