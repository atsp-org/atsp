*     ------------------------------------------------------------------
     
*       PROGRAM LINES  -- PROGRAM TO PRINT THE LINES IN THE
*			  LSJ TRANSITION ARRAY
*
*                   C O P Y R I G H T --1994
*
*                By C. Froese Fischer
*                   Bin Liu
*                   Vanderbilt University
*                   Nashville, TN 37235 USA
*
*                May, 1983
*
*     Computer Physics Communication, Vol. 64, 399-405 (1991)
*
*     ------------------------------------------------------------------
*
*       Output Table of Transition Probabilities
*
        PARAMETER (LINED=1000)
        CHARACTER*2 TYPE(LINED),TERM1,TERM2,TERM1P,TERM2P
        CHARACTER*6 ATOM
        CHARACTER*25 LAB1,LAB2,LAB1P,LAB2P
        CHARACTER*66 CONFIG(LINED,2)
        INTEGER JV(LINED,2), INDEX(LINED)
        DOUBLE PRECISION  ET(LINED,2)
        REAL V(LINED,6)
*
CTSS   CALL LINK('READ5, PRINT5, UNIT5=TERMINAL,
CTSS :           UNIT6=(out,create,text)//')
CTSS   CALL MPROMPT('>',1)
*
        WRITE(0,*)'Enter tolerance on line strength'
        READ(5, *) TOL
        OPEN(UNIT=1,FILE='tr.lsj',STATUS='OLD')
        N = 1
        READ(1,'(/10X,A6,7X,F3.0,7X,I3)') ATOM,Z,NEL
10      READ(1,'(/A)',END=999) LAB1
        IF (LAB1(20:20) .EQ. 'Z') GO TO 10
        READ(1,11) JV(N,1),ET(N,1),CONFIG(N,1)
11      FORMAT(I4,F14.8,2X,A)
        READ(1,11) JV(N,2),ET(N,2),CONFIG(N,2)
        READ(1,12) V(N,1),V(N,2),V(N,3),TYPE(N),V(N,4),V(N,5),V(N,6)
12      FORMAT(F11.2,7X,F11.2,12X,F11.2/
     :         1X,A2,6X,D12.5,8X,D12.5,9X,D12.5)
        IF (V(N,4) .LE. TOL) GO TO 10
        N = N+1
        IF (N .LE. (LINED)) GO TO 10
999     N = N-1
        WRITE(0,*) 'Number of transitions =',N
 
        DO 20 I = 1,N
           INDEX(I) = I
20      CONTINUE
 
        WRITE(0,'(A,7(/A))') ' Select the line list order:'
     :      , ' 1:  Energy (cm-1)'
     :      , ' 2:  Wavelength (Angstroms) in Vacuum'
     :      , ' 3:  Wavelength (Angstroms) in Air'
     :      , ' 4:  Line Strength'
     :      , ' 5:  gf Value'
     :      , ' 6:  Transition Probability'
     :      , ' Enter your selection: '
        READ(5,*) I
 
        DO 30 IPOS = 1,N-1
           J = IPOS
           DO 32 JPOS = IPOS+1,N
              JL = INDEX(J)
              JJ = INDEX(JPOS)
              IF (V(JJ,I) .LT. V(JL,I)) J = JPOS
32         CONTINUE
           ITEMP = INDEX(IPOS)
           INDEX(IPOS) = INDEX(J)
           INDEX(J) = ITEMP
30      CONTINUE
*
*               Data has been sorted
*
        WRITE(6,1) ATOM,Z,NEL
1       FORMAT(/' Line List for ',A6,' ( Z = ',F3.0,') with',
     :         I3,' electrons'//)
        WRITE(6,2)
2       FORMAT(80('-')/)
        WRITE(6,3)
3       FORMAT(1X,'Transition Array'/8X,'Multiplet',2x,'Line',
     :          3X,'Type',3X,'E(cm-1)',5X,'L(air)',6x,'S',
     :          6X,'gf',4X,'Aki')
        WRITE(6,2)
        LAB1P = ' '
        LAB2P = ' '
        TERM1P = ' '
        TERM2P = ' '
        DO 40 K = 1,N
           II = INDEX(K)
           CALL SEPAR(CONFIG(II,1),LAB1,TERM1)
           CALL SEPAR(CONFIG(II,2),LAB2,TERM2)
           IF (LAB1 .NE. LAB1P) THEN
                LAB1P = LAB1
                LAB2P = LAB2
                TERM1P = TERM1
                TERM2P = TERM2
                WRITE(6,*)
           ELSE
                LAB1 = ' '
                IF (LAB2 .NE. LAB2P) THEN
                    LAB2P = LAB2
                    TERM1P = TERM1
                    TERM2P = TERM2
                    WRITE(6,*)
                ELSE
                    LAB2 = ' '
                    IF (TERM1P .NE. TERM1) THEN
                        TERM1P = TERM1
                        TERM2P = TERM2
                    ELSE
                        TERM1 = ' '
                        IF (TERM2P .NE. TERM2) THEN
                            TERM2P = TERM2
                        ELSE
                            TERM2 = ' '
                        ENDIF
                    ENDIF
                ENDIF
            ENDIF
           IF (LAB1//LAB2 .NE. ' ') WRITE(6,35) LAB1,LAB2
	   WRITE(6,36) TERM1,TERM2,JV(II,1)/2.,JV(II,2)/2.,
     :          TYPE(II),V(II,1),V(II,3),V(II,4),V(II,5),V(II,6)
35         FORMAT(1X,A25,2X,A25)
36         FORMAT(8X,A3,1X,A3,2X,F4.1,'-',F4.1,
     :            2X,A2,F11.1,F11.1,2X,F6.3,2X,F6.3,1PE9.2)
40      CONTINUE
        END
* ----------------------------------------------------------------------
*       SEPARATE
* ----------------------------------------------------------------------
*
        SUBROUTINE SEPAR(CONFIG,LAB,TERM)
*
        CHARACTER CONFIG*66, LAB*25,TERM*2,CH2*2,CH1*1,COPY*66
 
                COPY = CONFIG
*
*       Delete the extra characters
*
                IF (INDEX(COPY,'(0)') .GT. 0) THEN
                    J = INDEX(COPY,'.')
                    COPY = CONFIG(J+1:66)
                ENDIF
 
*
*       1. SEPARATE the final term
*
                N = INDEX(COPY,' ')-1
                IF (COPY(N:N) .GE. '0' .AND. COPY(N:N) .LE.  '9') THEN
                    N = N-4
		    K = 4
                ELSE
                    N = N-3
		    K = 3
                ENDIF
                TERM = COPY(N+2:N+3)
               COPY(N+1:N+K) = '   '
 
*
*       2. Delete set subscript
*
                CH1 = COPY(N:N)
                IF (CH1.GE.'0' .AND. CH1.LE.'9') COPY(N:N) = ' '
 
*
*       3. If after removing the final term, there are no other
*   intermediate couplings prefaced by '_' and the last coupling
*   is the same as the final term, then the coupling for the TERMal
*   term is omitted .
*
                J = INDEX(COPY,'_')
                IF (J .EQ. 0) THEN
                    CH2 = COPY(N-2:N-1)
                    IF (CH2 .EQ. TERM) COPY(N-2:N-1) = '  '
                ENDIF
                LAB = COPY(1:30)
        END
 
