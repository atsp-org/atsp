*     ------------------------------------------------------------------
     
        PROGRAM GENCL
*                
*                C O P Y R I G H T -- 1994
*
*                By C. Froese Fischer
*                   Bin Liu
*                   Vanderbilt University
*                   Nashville, TN 37235 USA
*
*                Written: May, 1983
*    
*     Computer Physics Communications, Vol. 64, 406--416 (1991).
*     ------------------------------------------------------------------
*
*
*       This program computes all possible couplings for each member of
*   the reference set, generates all unique, possible configurations
     
*   and their couplings  from the active set, and for each replacement,
*   generates configurations and their couplings, then outputs the
*   configurations and couplings with given final term(s).
*     Input : (Interactive)
*        i)   Header ;
*       ii)   List of closed shells ;
*      iii)   The "reference" set ;
*       iv)   The "active" set ;
*        v)   Replacements from the reference set ;
*             Virtual Set if Replacement is 's' or 'd' or 'sd'
*       vi)   Final Term
*     Output :
*       i )   Header ;
*       ii)   List of closed shells ;
*      iii)   Configurations (FORMAT(5(1X,A3,'(',I2,')'))
*               and their couplings (FORMAT(9(5X,A3))
     
     
*   I/O allocation
*    FILE1(.)  :  Configurations by Active Set  (internal file)
*    FILE2(.)  :  Configurations by Replacement (internal file)
*    FILE3(.)  :  couplings                     (internal file)
*           6  :  Terminal output
*           7  :  File CLIST
*  FBETA(1,.)  :  Information of Beta2          (internal file )
*  FBETA(2,.)  :  Information of Beta3          (internal file )
*  FBETA(3,.)  :  Information of Beta4          (internal file )
*  FBETA(4,.)  :  Information of Beta5          (internal file )
*
*
*
* ---------------------------------------------------------------------
*               M A I N         P R O G R A M
* ---------------------------------------------------------------------
*
        PARAMETER      (NELS=15,NSHEL=5,NCOUPL=2*NSHEL-1)
        PARAMETER      (NCFG=500,NSCOUP=500)
        CHARACTER*60   REF(NELS),REPL(NELS),FINAL(NELS),ACT
        CHARACTER*60   STRL,STRR
        CHARACTER*72   HEADER,SHELLS,VIRTUL,TEMP
        CHARACTER*3    EL(NELS),ELL(NELS),ELR(NELS),ELS(NELS),
     :                 ELA(NELS)
        CHARACTER      FBETA*8, FILE1*32, FILE2*40, FILE3*32, NAME*24
        CHARACTER*3    COUPLE(NCOUPL),CH3,ELB(NELS),ELC(NELS),ELV(NELS)
        CHARACTER      CH1,CH2*2,FTM(NELS)*2
        INTEGER        ORDLA,ORDLZ,ORDUA,ORDUZ,ORD0,ORD9
        INTEGER        Q(NELS),QL(NELS),QR(NELS),
     :                 QS(NELS,NELS),MS(NELS),RL(NELS),
     :                 Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,
     :                 Q10,Q11,Q12,Q13,Q14,Q15
        INTEGER        F(NELS),PL,QA(NELS),QB(NELS),VL(NELS)
	INTEGER        H(NELS)
        INTEGER        LEFT,RIGHT,PARITY,CONST,SFLAG,DFLAG
        COMMON         NF,NR,NFTM,MAX,MIN,PARITY,CONST,NQ
     :                 /BLK0/ORDLA,ORDLZ,ORDUA,ORDUZ,ORD0,ORD9
     :                 /BLK1/HEADER,SHELLS,ACT,VIRTUL
     :                 /BLK2/REF,REPL,FTM
     :                 /BLK3/EL,ELL,ELR,ELS,ELA
     :                 /BLK4/Q,QL,ML,QR,MR,M,QS,MS,MA,RL,NREF,
     :                       Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,
     :                       Q10,Q11,Q12,Q13,Q14,Q15
                COMMON /FILES/FBETA(4,NSCOUP),FILE1(NCFG),FILE2(NCFG),
     :                        FILE3(NSCOUP)
     
     
***************
*       Declaration of variable for the input data
*
*        HEADER  =  header of output
*        SHELLS  =  Closed Shells
*           REF  =  Reference set
*           ACT  =  Active Set
*          REPL  =  Replacement
*        VIRTUL  = Virtual Set
*           FTM  = Final Term (NELS)
*
     
***************
CSUN            i = iargc( )
CSUN            if (i .eq. 0) then
                   NAME = 'cfg.inp'
CSUN            else
CSUN 		   call getarg(1,NAME)
CSUN            end if

*       Obtain ASCII value for the bound of character and digit
*
        ORDLA = ICHAR('a')
        ORDLZ = ICHAR('z')
        ORDUA = ICHAR('A')
        ORDUZ = ICHAR('Z')
        ORD0 = ICHAR('0')
        ORD9 = ICHAR('9')
        NQ = 0
        IFLAG = 0
        SFLAG = 0
        DFLAG = 0
     
****************
*       Disply for the beginning of the program
*
     
010     WRITE(0,*)
        WRITE(0,*)
        WRITE(0,*)
     :'     -----------------------------------------------------------'
        WRITE(0,*)
        WRITE(0,*) '     You are under the program GENCL'
        WRITE(0,*)  '            which GENerates a Configuration List'
        WRITE(0,*)
        WRITE(0,*)
        WRITE(0,*) 
     :'     Type H (Help) to obtain information about the input format '
        WRITE(0,*)
        WRITE(0,*) '     Type <RETURN> if you already know '
        WRITE(0,*)
        WRITE(0,*) '     -----------------------------------------------
     :---------'
        WRITE(0,*)
     
***************
*       If the user needs to obtain information about input format,
*    call subroutine HELP .
*
        READ(5,'(A)') STRL
        I = INDEX(STRL,'H')
        J = INDEX(STRL,'h')
        IF (I.NE.0 .OR. J.NE.0) CALL HELP
     
***********************************************************************
*            READ AND ANALYSIS THE INPUT DATA
*
*
*                    Input Header
*
100     WRITE(0,101)
101     FORMAT(T13,'Header  ?  ' )
        READ(5,'(A)') HEADER
     
                CH1 = CH3(2:2)
     
*
***************
*       Input Closed Shell, check if it satisfies the FORMAT(18(1X,A3))
*
110     WRITE(0,113)
113     FORMAT('     Closed Shells  ?  ' )
        READ(5,'(A)') TEMP
        I = INDEX(TEMP,'B')
        J = INDEX(TEMP,'b')
        IF (I.NE.0 .OR. J.NE.0) GOTO 100
        N = 2
        SHELLS = ' '
115     IF (INDEX(TEMP,'     ') .GT. 2) THEN
            CALL DEL(TEMP)
            J = ICHAR(TEMP(1:1))
            IF (J.LT.ORD0 .OR. J.GT.ORD9) THEN
                WRITE(0,*)  '     INPUT ERROR !'
                GOTO 110
            ENDIF
            SHELLS(N:N+2) = TEMP(:3)
            N = N+4
            CALL STRSH(TEMP,4)
            GO TO 115
        END IF
     
     
***************
*          Input Reference Set, and check the input error
*
*        LEFT  =  position of the first '(' in the string
*       RIGHT  =  position of the first ')' in the string
*         SUM  =  total number of Q(i) for each configuration
*          NL  =  Total quantum number for each reference set
*       CONST  =  Total number of Qi for the first reference set
*      PARITY  =  Parity for the first reference set
*
*
120     WRITE(0,121)
121     FORMAT('     Reference Set  ?  ' )
        CALL INPUT (NELS,NREF,REF,IFLAG,SFLAG,DFLAG,*110,*120)
        DO 122 I=1,NREF
            SUM = 0
            NL = 0
            TEMP = REF(I)
123         IF (TEMP(1:5) .NE. '     ') THEN
                CALL DEL(TEMP)
*
*       Error if the input has unmatched parenthesis
*
                LEFT = INDEX(TEMP,'(')
                RIGHT = INDEX(TEMP(:7),')')
                IF (LEFT.EQ.0 .OR. RIGHT.EQ.0) THEN
                    WRITE(0,*) '       Unmatched parenthesis!',
     :                 ' Please input again .'
                    GOTO 120
                ENDIF
*
*       Error if number of electron is more than FULL
*
                CH3 = TEMP(LEFT+1:RIGHT-1)
                N = ICTOI(CH3)
                L = LVAL(TEMP(2:2))
                IF (N .GT. L*4+2) THEN
                    WRITE(0,*) '     Number of electrons in a shell is
     :more than FULL !'
                    GOTO 120
                ENDIF
                SUM = SUM+N
                NL = NL+N*L
                CALL STRSH(TEMP,RIGHT+1)
                GO TO 123
            END IF
     
***************
*       The first Reference Set is defined as the standard value
*
            IF (I .EQ. 1) THEN
                CONST = SUM
                PARITY = MOD(NL,2)
            ELSE
*
*       Error if members of the reference set have different electron
*   numbers or parity
*
                IF (SUM .NE. CONST) THEN
                   WRITE(0,*) '     Total number of electrons is wrong!'
                   GOTO 120
                ENDIF
                IF (MOD(NL,2) .NE. PARITY) THEN
                   WRITE(0,*) '       Parity is wrong!'
                   GOTO 120
                ENDIF
            ENDIF
122     CONTINUE
     
*******************
*                  Input Active Set
*
130     WRITE(0,131)
131     FORMAT('        Active Set  ?  ' )
        READ(5,'(A)') ACT
        I = INDEX(ACT,'B')
        J = INDEX(ACT,'b')
        IF (I.NE.0 .OR. J.NE.0) GOTO 120
        IF (ACT .ne. '   ') then
	   WRITE(ISCW,132)
132        FORMAT('Type of set generation ?')
           READ *, ITYPE
	END IF
********************
*                  Input Replacement
*
140     WRITE(0,141)
141     FORMAT('      Replacements  ?  ')
        IFLAG = 1
        CALL INPUT (NELS,NREPL,REPL,IFLAG,SFLAG,DFLAG,*130,*150)
     
********************
*       If Replacement = S or D or SD, input Virtual Set
*
*      ELV  =  Electron label for virtual set
*       MV  =  Number of electron for Virtual Set
*       VL  =  Parity for each Qi in Virtual Set
*
142     IF (SFLAG.EQ.0 .AND. DFLAG.EQ.0) GOTO 150
     
        WRITE(0,143)
143     FORMAT ('       Virtual Set  ?  ')
        READ(5,'(A)') VIRTUL
        TEMP = VIRTUL
        J = INDEX(TEMP,'     ')
        TEMP(J:J) = ','
        MV = 0
     
*
*       Decompose the input of Virtual Set
*
144     IF (TEMP(:5) .NE. '     ') THEN
            MV = MV+1
            IF (MV.GT.(NELS)) THEN
	       WRITE(0,*) ' Virtual set too large: MAX=',NELS
	       STOP
	    END IF
            CALL DEL(TEMP)
            J = INDEX(TEMP,',')
            VL(MV) = MOD(LVAL(TEMP(2:2)),2)
*
*       Convert the input of uppercase to lowercase and assign value
*    ELVi
*
            N = ICHAR(TEMP(2:2))
            IF (N.GE.ORDUA .AND. N.LE.ORDUZ)
     :             TEMP(2:2) = CHAR(N-ORDUA+ORDLA)
            ELV(MV) = TEMP(:J-1)
            CALL STRSH(TEMP,J+1)
            GO TO 144
        END IF
     
145     WRITE(0,146)
146     FORMAT ('  From which shell  ?  ')
        READ(5,'(A)') TEMP
        I = INDEX(TEMP,'B')
        J = INDEX(TEMP,'b')
        IF (I.NE.0 .OR. J.NE.0) GOTO 140
        CALL DEL(TEMP)
        NVIR = ICHAR(TEMP(1:1))-ORD0
        IF (NVIR.LT.1 .OR. NVIR.GT.9) THEN
            WRITE(0,*) '    Please input digit for the shell position!'
            GOTO 145
        ENDIF
     
147     WRITE(0,148)
148     FORMAT ('    To which shell  ?  ')
        READ(5,'(A)') TEMP
        I = INDEX(TEMP,'B')
        J = INDEX(TEMP,'b')
        IF (I.NE.0 .OR. J.NE.0) GOTO 140
        CALL DEL(TEMP)
        LVIR = ICHAR(TEMP(:1))-ORD0
        IF (LVIR.LT.0 .OR. LVIR.GT.9) THEN
           WRITE(0,*) '     Please input digit for the shell position!'
           GOTO 147
        ENDIF
     
*******************
*                    Input Final Terms
*
*      NFTM  =  Number of input for final term
*
150     WRITE(0,151)
151     FORMAT ('       Final Terms  ?  ')
	IFLAG = 0
        CALL INPUT (NELS,NFTM,FINAL,IFLAG,SFLAG,DFLAG,*140,*160)
        DO 152 I=1,NFTM
            CH2 = FINAL(I)(:2)
     
*
*       Conver the input of lowercase to uppercase
*
            N = ICHAR(CH2(2:2))
            IF (N.GE.ORDLA .AND. N.LE.ORDLZ)
     :           CH2(2:2) = CHAR(N-ORDLA+ORDUA)
            FTM(I) = CH2
152     CONTINUE
     
***********************************************************************
*       Open the file for CI.LST
*
160     OPEN (7,FILE=NAME,STATUS='UNKNOWN')
     
***********************************************************************
*            PRINT OUT ALL INPUT DATA OF THE USER
*
        WRITE (6,170)
170     FORMAT (///////T5,'***************          I N P U T',
     :  '  D A T A          **********'///)
        WRITE (6,171) HEADER
171     FORMAT (T5,'          Header  :  ',A60/)
        WRITE (6,172) SHELLS
172     FORMAT (T5,'    Closed shell  :  ',A60/)
        WRITE (6,173) REF(1)
173     FORMAT (T5,'   Reference Set  :  ',A60/)
174     FORMAT (T29,I2,'  :  ',A60/)
        DO 175 I=2,NREF
            WRITE (6,174) I,REF(I)
175     CONTINUE
        WRITE (6,176) ACT
176     FORMAT (T5,'      Active Set  :  ',A60/)
        WRITE (6,178) REPL(1)
178     FORMAT (T5,'    Replacements  :  ',A60/)
        IF (SFLAG.NE.0 .OR. DFLAG.NE.0) THEN
            WRITE (6,179) VIRTUL
179         FORMAT (T5,'     Virtual Set  :  ',A72/)
            WRITE (6,*)'             From which shell  :  ',
     : CHAR(NVIR+ORD0)
            WRITE (6,*)
            WRITE (6,*)'               To which shell  :  ',
     : CHAR(LVIR+ORD0)
            WRITE (6,*)
        ENDIF
        DO 180 I=2,NREPL
            WRITE (6,174) I,REPL(I)
180     CONTINUE
        WRITE (6,181) FINAL(1)
181     FORMAT (T5,'     Final Terms  :  ',A60/)
        DO 182 I=2,NFTM
            WRITE (6,174) I,FINAL(I)
182     CONTINUE
     
        WRITE (7,183) HEADER
183     FORMAT (' ',A60)
        WRITE (7,'(1X,A72)') SHELLS
     
***********************************************************************
*                    PROCESS THE REFERENCE SET
*
*       M  =  Maximin number of ELi for reference set
*     ELS  =  Storage of ELi for different Reference Set
*      QS  =  Storage of Qi for each Reference Set
*      MS  =  Storage of M for each Reference Set
*    STRL  =  Temporary string
*
*
200     DO 201 NF=1,NREF
            M = 0
            TEMP = REF(NF)
            DO 207 I=1,NELS
     
***************
*       Decompose the input for Reference Set
*
                IF (TEMP(:5) .EQ. '     ') THEN
                    Q(I) = 0
                    GOTO 204
                ENDIF
                CALL DEL(TEMP)
                LEFT = INDEX(TEMP,'(')
                RIGHT = INDEX(TEMP,')')
                M = M+1
     
***************
*       Convert the input of uppercase to lowercase, and assign initial
*   values for ELi and Qi
*
                CH3 = TEMP(:LEFT-1)
                N = ICHAR(CH3(2:2))
                IF (N.GE.ORDUA .AND. N.LE.ORDUZ)
     :              CH3(2:2) = CHAR(N-ORDUA+ORDLA)
                EL(M) = CH3
                CH3 = TEMP(LEFT+1:RIGHT-1)
                Q(M) = ICTOI(CH3)
                CALL STRSH(TEMP,RIGHT+1)
207         CONTINUE
204         DO 212 J=1,NELS
                QS(J,NF) = 0
212         CONTINUE
     
***************
*       Add new electrons to the first configuration
*
            IF (NF .EQ. 1) THEN
               DO 211 J = 1,M
                  ELS(J) = EL(J)
                  QS(J,NF) = Q(J)
211            CONTINUE
               MS(NF) = M
            ELSE
               DO 210 J = 1,M
                  N = MS(1)
                  DO 208 I = 1,N
                     IF (EL(J) .EQ. ELS(I)) THEN
                        QS(I,NF) = Q(J)
                        GO TO 210
                     END IF
208               CONTINUE
                  N = N + 1
                IF (N.GT.NELS) THEN
                   WRITE(0,*) ' Too many shells in reference set: MAX=',
     :                NELS
                   STOP
                END IF
                  MS(1) = N
                  ELS(N) = EL(J)
                  QS(N,NF) = Q(J)
210            CONTINUE
               MS(NF) = N
            END IF
     
****************
*          Generate all couplings for the reference set
*          MAX = -5  means the first time to call subroutine COUPLD
*
            K = 0
            DO 206 I=1,M
                IF (Q(I) .NE. 0) THEN
                    K = K+1
                    ELB(K) = EL(I)
                    QB(K) = Q(I)
                ENDIF
206         CONTINUE
     
            MAX = -5
            CALL COUPLD (ELB,QB,K,ELC,NC,*500)
            IF (NC .NE. 0) THEN
                WRITE (6,202)
202             FORMAT (//'      GENERATE ALL COUPLINGS FOR EACH MEMBER
     : OF THE REFERENCE SET'//)
                CALL PRINT (ELC,QB,K,NC,1)
            ENDIF
     
***************
*            Compute MAX and MIN value for the finals .
*            Rule :  MAX = 2*|S+L| ,  MIN = 2*|S-L|
*
            IF (NF .EQ. 1) THEN
                MAX = 0
                MIN = 100
                DO 215 I=1,NC
                    N = 2*K-1
                    READ (FILE3(I),205) (COUPLE(J),J=1,N)
205                 FORMAT (9(A3))
                    CH3 = COUPLE(N)(1:1)
                    CH1 = COUPLE(N)(2:2)
                    S = (ICTOI(CH3)-1)/2.
                    L = LVAL(CH1)
                    LMIN = 2*ABS(S-L)
                    LMAX = 2*ABS(S+L)
                    IF (LMIN .LT. MIN) MIN = LMIN
                    IF (LMAX .GT. MAX) MAX = LMAX
215             CONTINUE
            ENDIF
201     CONTINUE
        IF (ACT(:5) .EQ. '     ') GOTO 400
     
***********************************************************************
*                    PROCESS THE ACTIVE SET
*
*       The member of the input is limited to NELS .
*
*      NQ  =  Number of configurations
*     ELA  =  ELi for Active Set
*      QA  =  Qi for Active Set
*      MA  =  Number of ELi for Active Set
*      RL  =  L-value for each shell
*       F  =  Full value for each shell
*
        MA = 0
300     TEMP = ACT
     
***************
*       Decompose the input of Active Set, Conver the input of uppercase
*    to lowercase, assign values to ELA,QA,RL,F .
*
301     IF (MA .LT. 15) THEN
            IF (TEMP(1:5) .EQ. '     ') GOTO 302
            CALL DEL(TEMP)
            MA = MA+1
            N = ICHAR(TEMP(2:2))
            IF (N.GE.ORDUA .AND. N.LE.ORDUZ)
     :          TEMP(2:2) = CHAR(N-ORDUA+ORDLA)
            ELA(MA) = TEMP(:2)
            RL(MA) = LVAL(TEMP(2:2))
	    IF (RL(MA) .GT. 3) THEN
	       F(MA) = 2
	    ELSE
               F(MA) = MIN0(const,4*RL(MA)+2)
	    END IF
            CALL STRSH(TEMP,4)
            GO TO 301
	ELSE
	   WRITE(0,*) 
     :    'Too many electrons in the active set: MAX= 15'
	   STOP
        END IF
     
****************
*       NELS is maximun value of active electrons
*
302     DO 303 I=MA+1,15
            F(I) = 0
303     CONTINUE
     
***************
*       Generate other possible configurations from the given set
*       F(I) is the maximum allowed number of electrons 
*       H(I) is the number of unassigned electrons or ``holes''
*
305     IF (itype .eq. 0 ) then 
	   MINQ1 = 0
	else if (itype .eq. 1) then
	   MINQ1 = max0(MIN0(F(1),CONST)-1,0)
	else if (itype .eq. 2) then
 	   MINQ1 = max0(MIN0(F(1),CONST)-2,0)
	else if (itype .eq. 3) then
 	   MINQ1 = max0(MIN0(F(1),CONST)-3,0)
	else
	   WRITE(ISCW,*) ' Unknown type: Re-enter'
	   READ *, itype
	   GO TO 305
	end if
        DO 310 Q1=MIN0(F(1),CONST),MINQ1,-1
         H(1) = MAX0(0,CONST - Q1)
         DO 310 Q2=MIN0(F(2),H(1)),0,-1
          H(2) = MAX0(0,H(1)-Q2)
          DO 310 Q3=MIN0(F(3),H(2)),0,-1
           H(3) = MAX0(0,H(2)-Q3)
           DO 310 Q4=MIN0(F(4),H(3)),0,-1
            H(4) = MAX0(0,H(3)-Q4)
            DO 310 Q5=MIN0(F(5),H(4)),0,-1
             H(5) = MAX0(0,H(4)-Q5)
             DO 310 Q6=MIN0(F(6),H(5)),0,-1
              H(6) = MAX0(0,H(5)-Q6)
              DO 310 Q7=MIN0(F(7),H(6)),0,-1
               H(7) = MAX0(0,H(6)-Q7)
               DO 310 Q8=MIN0(F(8),H(7)),0,-1
                H(8) = MAX0(0,H(7)-Q8)
                DO 310 Q9=MIN0(F(9),H(8)),0,-1
                 H(9) = MAX0(0,H(8)-Q9)
                 DO 310 Q10=MIN0(F(10),H(9)),0,-1
                  H(10) = MAX0(0,H(9)-Q10)
                  DO 310 Q11=MIN0(F(11),H(10)),0,-1
                   H(11) = MAX0(0,H(10)-Q11)
                   DO 310 Q12=MIN0(F(12),H(11)),0,-1
                    H(12) = MAX0(0,H(11)-Q12)
                    DO 310 Q13=MIN0(F(13),H(12)),0,-1
                     H(13) = MAX0(0,H(12)-Q13)
                     DO 310 Q14=MIN0(F(14),H(13)),0,-1
                      H(14) = MAX0(0,H(13)-Q14)
                      IF ( H(14) .LE. F(15)) THEN
                         Q15 = H(14)
                         CALL CONFIG
                      END IF
310     CONTINUE
***************
*                 Print Header of the output file
*
340     IF (NQ .EQ. 0) GOTO 400
        WRITE (6,341)
341     FORMAT (//'      GENERATE ALL POSSIBLE CONFIGURATIONS FROM THE
     : ACTIVE SET'//)
        J = INDEX(HEADER,'         ')
        WRITE (6,321) HEADER(1:J)
321     FORMAT (' '/T10,'-------------       ',A,'  --------'/)
        WRITE (6,322) SHELLS
322     FORMAT ('          Closed Shells  :  ',A60/)
        WRITE (6,323) ACT
323     FORMAT ('             Active Set  :  ',A60/)
     
***************
*       Print all configurations generated from Active Set
*
        WRITE (6,*) '         Configurations  :'
        DO 325 I=1,NQ
            READ (FILE1(I),326) (QA(J), J=1,MA)
326         FORMAT (15(I2))
            IF (ELA(1)(1:1) .EQ. ' ') THEN
                N = 27
            ELSE
                N = 28
            ENDIF
     
            K = 0
            DO 370 J=1,MA
                IF (QA(J).NE.0 .OR. MA.EQ.2) THEN
                    K = K+1
                    ELB(K) = ELA(J)
                    QB(K) = QA(J)
                ENDIF
370         CONTINUE
            WRITE (6,327) (ELB(J),QB(J), J=1,K)
327         FORMAT (T28,5(1X,A3,'(',I2,')'))
325     CONTINUE
     
****************
*       For each new configuration, generate all possible couplings
*
        N1 = 0
        N2 = 0
        DO 328 I=1,NQ
            READ (FILE1(I),326) (QA(J), J=1,MA)
     
***************
*               Omit ELA(i) if corresponding QA(i)=0
*
            K = 0
            DO 350 J=1,MA
                IF (QA(J).NE.0 .OR. MA.EQ.2) THEN
                    K = K+1
                    ELB(K) = ELA(J)
                    QB(K) = QA(J)
                ENDIF
350         CONTINUE
     
****************
*       Omit configurations which have more than 5 shells
*
            IF (K .LE. 5) THEN
                CALL COUPLD (ELB,QB,K,ELC,NC,*500)
                IF (NC .GT. 0) THEN
                    WRITE (6,343)
343                 FORMAT (//'    LIST THE COUPLINGS FOR EACH',
     :               ' CONFIGURATION GENERATED BY THE ACTIVE SET'//)
                    CALL PRINT (ELC,QB,K,NC,2)
                ELSE
                    N2 = N2+1
                ENDIF
            ELSE
                N1 = N1+1
            ENDIF
328     CONTINUE
        IF (N1 .NE. 0) PRINT 344,N1
344     FORMAT (T5,'Too many occuplied shells --- ',I3,
     : ' configuration omitted!')
        IF (N2 .NE. 0) PRINT 345, N2
345     FORMAT (T5,'No final term as your selection for ',I3,
     :' Active set!')
     
349     IF (NREPL .EQ. 0) GOTO 500
     
***********************************************************************
*                    PROCESS THE REPACEMENTS
*
*       STRL  :  String of the left side of '='
*       STRR  =  String of the right side of '='
*
     
400     IF (SFLAG.NE.0 .OR. DFLAG.NE.0) GOTO 450
        N1 = 0
        DO 401 NR=1,NREPL
            STRR = REPL(NR)
            J = INDEX(STRR,'=')
            STRL = STRR(1:J-1)
            CALL STRSH(STRR,J+1)
     
***************
*               Decompose the substring on the left of ''=''
*
*     ELL  =  Old value of ELi to be replaced
*      QL  =  Old value of Qi to be replaced
*      ML  =  Number of Old value of ELi to be replaced
*
            CALL DEL(STRL)
            CALL DECOMP (STRL,ELL,QL,ML)
     
***************
*               Decompose the substring on the right of ''=''
*
*     ELR  =  New value of ELi
*      QR  =  New vaRue of Qi
*      MR  =  Number of the new value of ELi
*
            CALL DEL(STRR)
            CALL DECOMP (STRR,ELR,QR,MR)
     
***************
*            For each Replacement, replace all Reference Sets
*
            DO 401 NF=1,NREF
                M = MS(NF)
                IF (M .LT. ML) GOTO 401
                DO 402 I=1,M
                    EL(I) = ELS(I)
                    Q(I) = QS(I,NF)
402             CONTINUE
                CALL REPLAC(ELB,QB,K,*401)
                CALL COUPLD (ELB,QB,K,ELC,NC,*500)
                IF (NC .GT. 0) THEN
                    WRITE (6,403)
403                 FORMAT (//'      FOR EACH REPLACEMENT, GENERATE',
     : ' CONFIGURATIONS AND COUPLINGS FOR THE REFERENCE SET'//)
                    CALL PRINT (ELC,QB,K,NC,3)
                ELSE
                    N1 = N1+1
                ENDIF
401     CONTINUE
        IF (N1 .NE. 0) WRITE(0,405) N1
405     FORMAT (T5,'No Final Term as selected for ',I3,
     :' Replacements!')
        GOTO 500
     
***********************************************************************
*               PROCESS THE VIRTUAL SET
*
*     RL  =  Parity of ELi of Reference Set
*     VL  =  Parity of ELi of Virtual Set
*
*
450     DO 451 NF=1,NREF
     
***************
*       Set the initial values of Reference Set
*
            M = MS(NF)
            DO 452 I=1,NELS
                CH3 = ELS(I)
                EL(I) = CH3
                ELB(I) = CH3
                Q(I) = QS(I,NF)
                QB(I) = Q(I)
                RL(I) = MOD(LVAL(CH3(2:2)),2)
452         CONTINUE
     
            IF (SFLAG.EQ.0 .AND. DFLAG.NE.0) GOTO 460
     
***************
*       Preparation for the Single Replacement
*
            QL(1) = 1
            ML = 1
            QR(1) = 1
            MR = 1
            DO 454 I=NVIR,LVIR
                IF (Q(I).EQ.0) GOTO 454
                ELL(1) = EL(I)
                IF (RL(I).GT.0) THEN
                    N = 1
                ELSE
                    N = 0
                ENDIF
                DO 455 J=1,MV
     
***************
*       Replace Qi of Reference Set which has the same parity with Qj
*    of Virtual Set
*
                    IF (VL(J) .EQ. N) THEN
                        ELR(1) = ELV(J)
                        CALL REPLAC(ELB,QB,K,*455)
                        CALL COUPLD (ELB,QB,K,ELC,NC,*500)
                        IF (NC .GT. 0) THEN
                            WRITE (6,453)
453                         FORMAT (//'        FOR VIRTUAL SET, ',
     : 'GENERATE CONFIGURATION AND COUPLINGS FOR S-REPLACEMENT'//)
                            NR = NELS
                            TEMP = ELL(1)//' = '//ELR(1)
                              CALL DEL(TEMP)
                              REPL(NR) = TEMP
                            CALL PRINT (ELC,QB,K,NC,4)
                        ENDIF
                    ENDIF
455             CONTINUE
454         CONTINUE
     
     
**************
*       Rreplace pairs of Q(i) and Q(j) by Double Virtual Set
*
*     PL  =  Pairty for the pair Qi and Qj in Reference Set
*
460         IF (DFLAG .EQ. 0) GOTO 500
            ML = 2
            QL(1) = 1
            QL(2) = 1
            DO 461 I=NVIR,LVIR-1
                ELL(1) = EL(I)
		LL1 = LVAL(ELL(1)(2:2))
                DO 461 J=I+1,M
                    IF (Q(I).EQ.0 .OR. Q(J).EQ.0) GOTO 461
                    ELL(2) = EL(J)
		    LL2 = LVAL(ELL(2)(2:2))
                    TEMP = ELL(1)//'.'//ELL(2)//' = '
		    LLMIN = IABS(LL1 - LL2)
		    LLMAX = LL1 + LL2
                    PL = MOD(LLMAX, 2)
                    CALL VPAIR (ELV,MV,PL,LLMIN,LLMAX,TEMP,*500)
461         CONTINUE
     
***************
*       Replace pairs of (Qi)=2 by Double Virtual Set
*
            ML = 1
            QL(1) = 2
            DO 464 I=NVIR,LVIR
                IF (Q(I) .GT. 1) THEN
                    ELL(1) = EL(I)
		    LL1 = LVAL(ELL(1)(2:2))
		    LLMIN = 0
		    LLMAX = LL1 + LL1
                    TEMP = ELL(1)//'(2) = '
                    PL = MOD(LLMAX, 2)
                    CALL VPAIR (ELV,MV,PL,LLMIN,LLMAX,TEMP,*500)
                ENDIF
464         CONTINUE
451     CONTINUE
     
***********************************************************************
*                    THE END OF THE PROGRAM
*
500     CLOSE (7)
        WRITE(0,*)
        WRITE(0,*) '         OK !'
        WRITE(0,*) '         List of configurations and their couplings'
        WRITE(0,*) '         is in the file ', NAME
        END
*
* ----------------------------------------------------------------------
*               FUNCTION        D E L
* ----------------------------------------------------------------------
*
*       Delete the leading space of the string
*
        SUBROUTINE DEL(STR)
        CHARACTER*(*)   STR
        CHARACTER*72 TEMP
     
            LENGTH = LEN(STR)
            I = 0
10          IF (STR(I+1:I+1) .EQ. ' ') THEN
               I = I+1
                IF (I .LT. LENGTH) GO TO 10
        END IF
            TEMP = STR(I+1:)
            STR = TEMP
            RETURN
        END
*
* ----------------------------------------------------------------------
*               SUBROUTINE        S T R S H
* ----------------------------------------------------------------------
*
*       Shift the string left
*
        SUBROUTINE STRSH(STR,I)
        CHARACTER*(*)   STR
        CHARACTER*72 TEMP
     
            TEMP = STR(I:)
            STR = TEMP
            RETURN
        END
*
* ----------------------------------------------------------------------
*               FUNCTION        I C T O I
* ----------------------------------------------------------------------
*
*       Convert character string into its corresponding integer
*
        INTEGER FUNCTION ICTOI(STR)
        CHARACTER*(*)   STR
     
            N = ICHAR(STR(1:1))-ICHAR('0')
            IF (STR(2:2) .NE. ' ') N = N*10+ICHAR(STR(2:2))-ICHAR('0')
            ICTOI = N
            RETURN
        END
*
* ----------------------------------------------------------------------
*               FUNCTION        L V A L
* ----------------------------------------------------------------------
*
*       convert the symbol into its corresponding quantum number
*
        INTEGER FUNCTION LVAL(SYMBOL)
        CHARACTER      SYMBOL
        CHARACTER*22   SET
        DATA SET/'spdfghiklmnSPDFGHIKLMN'/
     
            LOCATE = INDEX(SET,SYMBOL)
            IF (LOCATE .LE. 11) THEN
                LVAL = LOCATE-1
            ELSE
                LVAL = LOCATE-12
            ENDIF
            RETURN
         END
*
* ----------------------------------------------------------------------
*               FUNCTION        S Y M B
* ----------------------------------------------------------------------
*
*       Convert the quantum number into its corresponding symbol
*
        CHARACTER FUNCTION SYMB(L)
            CHARACTER*11   SET
            DATA SET/'SPDFGHIKLMN'/
     
            SYMB = SET(L+1:L+1)
            RETURN
        END
*
* ----------------------------------------------------------------------
*               SUBROUTINE      I N P U T
* ----------------------------------------------------------------------
*
*       Process the input set and check the input error
*
        SUBROUTINE INPUT (MAXSET,NSET,SET,MARK,SFLAG,DFLAG,*,*)
            CHARACTER*60   SET(*),TEMP,CH1,ACT
            CHARACTER*72   HEADER,SHELLS,VIRTUL
            INTEGER        SFLAG,DFLAG,MARK
            COMMON         /BLK1/HEADER,SHELLS,ACT,VIRTUL
     
*   MAXSET  =  Maximum number of input elements
*     NSET  =  Number of members in the set
*      SET  =  Character array with NSET elements
*       *1  =  Return label if input is 'B'
*       *2  =  Return label if the set is empty
*     MARK  =  1 if input is Replacement, 0 otherwise
     
            NSET = 0
            SET(1) = '    '
12          IF (NSET .LT. MAXSET) THEN
                READ(5,'(A60)') TEMP
     
***************
*       If input is 'B' or 'b', go back one step
*
                I = INDEX(TEMP,'B')
                J = INDEX(TEMP,'b')
                IF (I.NE.0 .OR. J.NE.0) THEN
                    IF (NSET .EQ. 0) THEN
                        RETURN 1
                    ELSE
                        NSET = NSET-1
                        GOTO 13
                    ENDIF
                ENDIF
     
***************
*       Go to for next input if the input is empty
*       Return if the input is finished
*
                IF (TEMP(1:5) .EQ. '     ') THEN
                    IF (NSET .EQ. 0) THEN
                        RETURN 2
                    ELSE
                        RETURN
                    ENDIF
                ENDIF
                CALL DEL(TEMP)
     
***************
*       If Replacement is 's' or 'd' or 'sd', set single and
*    double flag for Virtual Set
*
                IF (MARK.NE.0) THEN
                    CH1 = TEMP(:1)
                    IF (CH1.EQ.'S' .OR. CH1.EQ.'s') SFLAG = 1
                    I = INDEX(TEMP,'SD')
                    J = INDEX(TEMP,'sd')
                    IF (CH1.EQ.'D' .OR. CH1.EQ.'d' .OR.
     :                  I.NE.0 .OR. J.NE.0 )
     :                 DFLAG = 1
                    IF (SFLAG.NE.0 .OR. DFLAG.NE.0) THEN
                        NSET = 1
                        SET(1) = TEMP
                        RETURN
                    ENDIF
                ENDIF
     
***************
*       READ the input and delete the repeated member
*
                DO 14 I=1,NSET
                    IF (SET(I) .EQ. TEMP) THEN
                        WRITE(0,*) '     You give a repeated input!'
                        GOTO 13
                    ENDIF
14              CONTINUE
                NSET = NSET+1
                SET(NSET) = TEMP
13              WRITE(0,16) NSET+1,'  ?  '
16              FORMAT(T7,I10,A)
                GO TO 12
            END IF
            RETURN
        END
*
* ----------------------------------------------------------------------
*               SUBROUTINE      C O U P L D
* ----------------------------------------------------------------------
*
*          This subroutine generates all possible couplings .
*       First, compute Alpha value from the given configuration,
*       then compute Beta from each value of Alpha .
     
        SUBROUTINE COUPLD (EL,Q,M,ELC,NC,*)
     
*   INPUT :
*       El  =  electron label
*              where EL(I)(1=1)  ---  blank
*                    EL(I)(2=2)  ---  n-symbol
*                    EL(I)(3=3)  ---  L-symbol
*       Q  =  Occupation number
*               0 (empty) <= Q(i)  <= 2(2Li+1) (full)
*       M  =  number of shells
*               0  <=  M  <=  5
*   OUTPUT :
*       NC  =  number of couplings
*        *  =  Return label if the maximun number of couplings > NSCOUP
*
            PARAMETER      (NELS=15,NSHEL=5,NCOUPL=2*NSHEL-1)
            PARAMETER      (NCFG=500,NSCOUP=500)
            CHARACTER*60   REF(NELS),REPL(NELS)
            CHARACTER*65   TERM(22)
            CHARACTER*3    EL(*),ELC(*),ALFA(NSHEL,NSCOUP)
            CHARACTER*3    COUPLE(NCOUPL,NSCOUP),CH3,CCH3
            CHARACTER*2    CH2,A2,FTM(NELS)
            CHARACTER      CH1,CCH1,CALFA*15,SYMB,B1*2
            CHARACTER      FBETA*8, FILE1*32, FILE2*40, FILE3*32
            INTEGER        ORDLA,ORDLZ,ORDUA,ORDUZ,ORD0,ORD9
            INTEGER        Q(*),NTERM(22),LPOSIT(9),POSIT(NELS)
            INTEGER        PTR,PARENT,CHILD,BETA(NELS),PARITY,CONST
            COMMON         NF,NR,NFTM,MAX,MIN,PARITY,CONST,NQ
     :                     /BLK0/ORDLA,ORDLZ,ORDUA,ORDUZ,ORD0,ORD9
     :                     /BLK2/REF,REPL,FTM
                COMMON /FILES/FBETA(4,NSCOUP),FILE1(NCFG),FILE2(NCFG),
     :                        FILE3(NSCOUP)
     
***************
*           Number of possible terms for configurations P(1-3),
*           D(1-5), F(1-2), G(1-2), ... M(1-2)
*
        DATA  (NTERM(I), I=1,22) / 1,3,3,1,5,8,16,16,1,7,1,9,
     :                             1,9,1,9,1,9,1,9,1,9 /
***************
*           Starting position in term table for given L
*
	DATA  (LPOSIT(I), I=1,9)/1,4,9,11,13,15,17,19,21/

***************
*        Possible terms for configurations P(1-3),D(1-5),F(1-2),G(1-2)
*
	DATA (TERM(I), I=1,22) /
     :'2P1',
     :'1S0 1D2 3P2',
     :'2P1 2D3 4S3',
     :'2D1',
     :'1S0 1D2 1G2 3P2 3F2',
     :'2D1 2P3 2D3 2F3 2G3 2H3 4P3 4F3',
     :'1S0 1D2 1G2 3P2 3F2 1S4 1D4 1F4 1G4 1I4 3P4 3D4 3F4 3G4 3H4 5D4',
     :'2D1 2P3 2D3 2F3 2G3 2H3 4P3 4F3 2S5 2D5 2F5 2G5 2I5 4D5 4G5 6S5',
     :'2F1', '1S0 1D2 1G2 1I2 3P2 3F2 3H2',
     :'2G1', '1S0 1D2 1G2 1I2 1L2 3P2 3F2 3H2 3K2',
     :'2H1', '1S0 1D2 1G2 1I2 1L2 3P2 3F2 3H2 3K2',
     :'2I1', '1S0 1D2 1G2 1I2 1L2 3P2 3F2 3H2 3K2',
     :'2K1', '1S0 1D2 1G2 1I2 1L2 3P2 3F2 3H2 3K2',
     :'2L1', '1S0 1D2 1G2 1I2 1L2 3P2 3F2 3H2 3K2',
     :'2M1', '1S0 1D2 1G2 1I2 1L2 3P2 3F2 3H2 3K2' /
     
***********************************************************************
*               COMPUTE THE POSSIBLE VALUES FOR ALPHA
*
*   EMPTY  =  1 when Q(i)=0
*    FULL  =   2*(2L+1)
*    HALF  =  2L+1
*    ALFA  =  matrix of (NELS,NSCOUP)
*   NALFA  =  Number of ALFA
*   POSIT  =  array of NELS, store the position in table NTERM and TERM
*             corresponding Q(i). Rule  :
*             position = (L-1)*2+Q(i)        if 1 <= Qi <= HALF
*             position = (L-1)*2+(FULL-Qi)   if Qi > HALF
*
     
            NALFA = 1
            DO 10 I=1,M
                CH3 = EL(I)
               CH1 = CH3(2:2)
                IF (CH3(3:3) .EQ. ' ') THEN
                    CCH3 = CH3
                    CH3 = ' '//CCH3(1:2)
                END IF
                ELC(I) = CH3
                CCH1 = CHAR(ICHAR(CH1)-ORDLA+ORDUA)
                 CH1 = CCH1
                L = LVAL(CH1)
                FULL = 4*L+2
                K = Q(I)
     
***************
*               If shell is full, ALFA(i) = 1S0
*
                IF (K.EQ.0 .OR. K.EQ.FULL) THEN
                    CH3 = '1S0'
     
***************
*               If Q(i) = 1,  then ALFA(i)=2<L-symbol>1
*
                ELSE IF (K .EQ. 1) THEN
                    CH3 = '2'//CH1//'1'
     
***************
*       Otherwise, get the possible value from array NTERM and TERM
*
                ELSE
                    HALF = FULL/2
                    IF (K .LE. HALF) THEN
                        POSIT(I) = LPOSIT(L)+K-1
                    ELSE
                        POSIT(I) = LPOSIT(L)+(FULL-K)-1
                    ENDIF
                    NALFA = NALFA*NTERM(POSIT(I))
                    CH3 = '   '
                ENDIF
                IF (NALFA .GT. NSCOUP) THEN
                   WRITE(0,*)  'Array ALFA in routine COUPLD exceeded'
                   STOP
                END IF

***************
*       CALFA is a string storing the elements in one ALFA
*
                CALFA(3*I-2:3*I) = CH3
10          CONTINUE
     
***************
*            Assign values to all elements of ALFA
*
*       NT  =
*     LOCT  =  Current position in the table TERM
*     LOCA  =  Current position in the matrix ALFA
*
            NT = 1
            DO 12 I=M,1,-1
                CH3 = CALFA(I*3-2:I*3)
                IF (CH3 .NE. '   ') THEN
                    DO 11 J=1,NALFA
                        ALFA(I,J) = CH3
11                  CONTINUE
                ELSE
                    LOCT = POSIT(I)
                    N = NTERM(LOCT)
                    LOCA = 1
15                  IF (LOCA .LE. NALFA) THEN
                        DO 13 J=1,N
                            CH3 = TERM(LOCT)(J*4-3:J*4-1)
                            DO 13 K=1,NT
                                ALFA(I,LOCA) = CH3
                                LOCA = LOCA+1
13                      CONTINUE
                        GO TO 15
                    END IF
                    NT = NT*N
                ENDIF
12          CONTINUE
     
***********************************************************************
*       GENERATE POSSIBLE VALUE OF BETA FROM ALFA
*
     
            NC = 0
            DO 20 NB=1,NALFA
     
***************
*       There is only one coupling if M = 1
*
                IF (M .EQ. 1) THEN
                    BETA(1) = 1
                    COUPLE(1,1) = ALFA(1,NB)
                    NBETA = 1
                    GOTO 37
                ENDIF
     
***************
*       Define BETA(1)=ALFA(1), then the next basic coupling steps is :
*       S1 = (BETA(1)(1:1)-1)/2 ,   S2 = (ALFA(2)(1:1)-1)/2 ;
*               | S1-S2 | <= BETA(2)(1:1) <= | S1+S2 |
*       L1 = L-number of BETA(1)(2:2),   L2 = L-number of ALFA(2)(2:2) ;
*       Symbol(| L1-L2 |) <= BETA(2)(2:2) <= Symbol(| L1+L2 |) .
*
                B1 = ALFA(1,NB)(1:2)
                DO 21 J=2,M
                    A2 = ALFA(J,NB)(1:2)
                    PARENT = 1
                    CHILD = 1
                    BETA(J) = 0
26                  S1 = (ICHAR(B1(1:1))-ORD0-1)/2.
                    S2 = (ICHAR(A2(1:1))-ORD0-1)/2.
                    S3 = ABS(S1-S2)
                    S4 = ABS(S1+S2)
                    L1 = LVAL(B1(2:2))
                    L2 = LVAL(A2(2:2))
                    L3 = ABS(L1-L2)
                    L4 = ABS(L1+L2)
                    MBETA = (S4-S3+1)*(L4-L3+1)
*
     
**************
*       Generate Beta from each alpha .
*       There are four scratch files for storing the information about
*   Beta(i), 1<i<6, shown as follows :
*           --------------------------------------------
*           |    Parent    |    Value    |    Child    |
*           --------------------------------------------
*       Define   Beta(i) is child of Beta(i-1) and parent of beta(i+1) ;
*   Parent is a pointer to the parent of Beta(i), that is Beta(i-1) ;
*   Value is one of the possible value for BETA(i), and Child is the
*   number of children for Beta(i) .
*
*       MBETA  =  Number of couplings generated from ALFA(i)
*     PARENT  =  Current pointer to the parent of Beta(i)
*
                    DO 22 S=S3,S4
                        KS = 2*S
                        CH1 = CHAR(ORD0+1+KS)
                        DO 22 L=L3,L4
                            WRITE (FBETA(J-1,CHILD),25)
     :                           PARENT,CH1//SYMB(L),1
25                          FORMAT(I3,A2,I3)
                            CHILD = CHILD+1
22                  CONTINUE
50                  BETA(J) = BETA(J)+MBETA
                    IF (J .EQ. 2) GOTO 23
     
***************
*       Correct the number of its children for Beta(j-1)
*
                    LOCB = PARENT
                    DO 24 K=J-1,2,-1
                        READ (FBETA(K-1,LOCB),25) PTR,CH2,N
                        WRITE (FBETA(K-1,LOCB),25) PTR,CH2,N+MBETA-1
                        LOCB = PTR
24                  CONTINUE
     
***************
*       If pointer to the parent is the end of file for Beta(j-1),
*    prepare to generate Beta(i+1) ; otherwise, generate next Beta(j)
*    according Alfa(j-1) and Beta(j-1) .
*
52                  PARENT = PARENT+1
                    IF (PARENT .GT. BETA(J-1)) THEN
                       GOTO 23
                    ELSE
                        READ (FBETA(J-2,PARENT),25) PTR,B1,N
                        GOTO 26
                    ENDIF
23                  READ (FBETA(J-1,1),25) PTR,B1,N
21              CONTINUE
     
***************
*       Assign values to the couplings forward ,
*       COUPLE(I) = Alpha(i) for coupling(j) if i <= M ;
*
                NBETA = BETA(M)
                IF (NBETA .GT. NSCOUP) THEN
                    WRITE(0,*)  'Array COUPLE in routine COUPL exceeded'
                    STOP
                END IF
55              DO 27 J=1,NBETA
                    DO 31 I=1,M
                        COUPLE(I,J) = ALFA(I,NB)
31                  CONTINUE
                    DO 32 K=2*M-1,9
                        COUPLE(K,J) = '   '
32                  CONTINUE
27              CONTINUE
     
***************
*       Assign values to the couplings backward ,
*       COUPLE(i) = Beta(i-m) for coupling(j) if i > m
*
53              DO 28 I=M,2,-1
                    N = 1
                    DO 28 J=1,BETA(I)
                        READ (FBETA(I-1,J),25) PTR,CH2,NT
                        DO 28 K=1,NT
                            COUPLE(M+I-1,N) = CH2//'0'
                            N = N+1
                            IF (N .GT. NSCOUP) THEN
                                WRITE(0,*)  'Array COUPLE exceeded'
                                STOP
                            END IF
28              CONTINUE
     
***************
*    Selection from generated couplings according the following rules
*
37              DO 29 I=1,NBETA
                    N = 2*M-1
     
***************
*       If the first time to call COUPLD, not compute MAX and MIN
*
                    IF (MAX .EQ. -5) GOTO 34
     
***************
*       Compute MAX and MIN value for each final term, keep it if
*   intersection is non-empty
*
                    CH2 = COUPLE(N,I)(1:2)
                    CH3 = CH2(1:1)
                    CH1 = CH2(2:2)
                    S = (ICTOI(CH3)-1)/2.
                    L = LVAL(CH1)
                    LMIN = 2*ABS(S-L)
                    LMAX = 2*ABS(S+L)
                    IF (LMIN.GT.MAX .OR. LMAX.LT.MIN) GOTO 29
     
***************
*       If Final Terms are given, do selection
*
34                  IF (NFTM .NE. 0) THEN
                        CH2 = COUPLE(N,I)(:2)
                        DO 40 K=1,NFTM
                            IF (CH2 .EQ. FTM(K)) GOTO 42
40                      CONTINUE
41                      GOTO 29
                    ENDIF
     
***************
*       Waining if the number of couplings > NSCOUP
*
                    IF (NC .EQ. NSCOUP) THEN
                        WRITE(0,*)  '          WARNING !'
                        WRITE(0,*)  '          The number of couplings', 
     :                   ' is greater than 500 .',
     :                   '          Please select the Final Term .'
                        RETURN 1
                    ENDIF
     
***************
*       Write configurations and couplings to CI.LST
*
42                  NC = NC+1
                    WRITE (FILE3(NC),30) (COUPLE(J,I),J=1,N)
30                  FORMAT(9(A3))
                    WRITE (7,35) (ELC(J),Q(J), J=1,M)
35                  FORMAT (5(1X,A3,'(',I2,')'))
                    WRITE (7,36) (COUPLE(J,I), J=1,N)
36                  FORMAT (9(5X,A3))
29              CONTINUE
20          CONTINUE
            RETURN
        END
*
* ----------------------------------------------------------------------
*               SUBROUTINE      C O N F I G
* ----------------------------------------------------------------------
*
*               Examine if the new configuration has the same
*       electron number and parity .
*
        SUBROUTINE CONFIG
        PARAMETER     (NELS=15,NSHEL=5,NCOUPL=2*NSHEL-1)
        PARAMETER      (NCFG=500,NSCOUP=500)
        CHARACTER*3    EL(NELS),ELL(NELS),ELR(NELS),ELS(NELS),
     :                 ELA(NELS)
        CHARACTER*3    CH3
        CHARACTER      FBETA*8, FILE1*32, FILE2*40, FILE3*32
        INTEGER        QA(NELS),PARITY,CONST
        INTEGER        Q(NELS),QL(NELS),QR(NELS),
     :                 QS(NELS,NELS),MS(NELS),RL(NELS),
     :                 Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,
     :                 Q10,Q11,Q12,Q13,Q14,Q15
        COMMON         NF,NR,NFTM,MAX,MIN,PARITY,CONST,NQ
     :                 /BLK3/EL,ELL,ELR,ELS,ELA
     :                 /BLK4/Q,QL,ML,QR,MR,M,QS,MS,MA,RL,NREF,
     :                       Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,
     :                       Q10,Q11,Q12,Q13,Q14,Q15
        COMMON         /FILES/FBETA(4,NSCOUP),FILE1(NCFG),FILE2(NCFG),
     :                        FILE3(NSCOUP)
        EQUIVALENCE    (Q1,QA(1)),(Q2,QA(2)),(Q3,QA(3)),(Q4,QA(4)),
     :                    (Q5,QA(5)),(Q6,QA(6)),(Q7,QA(7)),(Q8,QA(8)),
     :                    (Q9,QA(9)),(Q10,QA(10)),(Q11,QA(11)),
     :                    (Q12,QA(12)),(Q13,QA(13)),(Q14,QA(14)),
     :                    (Q15,QA(15))
     
*   INPUT :
*     Q,Qi  =  Occupation number
*        L  =  L-value corresponding ELi
*   OUTPUT :
*       NQ  =  Number of configurations
     
***************
*       Return if the new configuration has different pairty
*
            NEWP = 0
            DO 11 I=1,MA
                NEWP = NEWP+QA(I)*RL(I)
11          CONTINUE
     
            NEWP = MOD(NEWP,2)
            IF (NEWP .NE. PARITY) RETURN
     
***************
*   Return if the new condiguration is the same as the Reference Set
*
            DO 12 I=1,NREF
                M = MS(I)
                N = 0
                DO 14 J=1,M
                    CH3 = ELS(J)
                    IQ = QS(J,I)
                    DO 14 K=1,MA
                        IF (CH3.EQ.ELA(K) .AND. IQ.EQ.QA(K)) N = N+1
14              CONTINUE
                IF (N .EQ. M) RETURN
12          CONTINUE
     
***************
*       Otherwise, write them into the configuration file
*
            NQ = NQ+1
            WRITE (FILE1(NQ),10) (QA(J), J=1,MA)
10          FORMAT (15(I2))
            RETURN
        END
*
* ----------------------------------------------------------------------
*       SUBROUTINE      D E C O M P
* ----------------------------------------------------------------------
*
*       Decompose the string of Replacement
*
        SUBROUTINE DECOMP (STR,EL,Q,MR)
            PARAMETER      (NELS=15,NSHEL=5,NCOUPL=2*NSHEL-1)
            PARAMETER      (NCFG=500,NSCOUP=500)
            CHARACTER*60   STR
            CHARACTER*3    CH3,EL(NELS)
            INTEGER        Q(NELS),LEFT,RIGHT
            INTEGER        ORDLA,ORDLZ,ORDUA,ORDUZ,ORD0,ORD9
            COMMON         /BLK0/ORDLA,ORDLZ,ORDUA,ORDUZ,ORD0,ORD9
     
*
*       INPUT  :
*         STR  =  String to be decomposed
*         EL and Q
*       OUTPUT  :
*           MR  =  Number of EL to be replaced
*
            DO 10 I=1,5
                IF (STR(:5) .EQ. '     ') THEN
                    MR = I-1
                    RETURN
                ENDIF
                CALL DEL(STR)
     
                LEFT = INDEX(STR,'(')
                RIGHT = INDEX(STR,')')
     
***************
*       If the Replacement is like 2s.2p = 3s.3p
*
                IF (LEFT .EQ. 0) THEN
                    K = 1
                    N = ICHAR(STR(3:3))
                    IF (N.GE.ORD0 .AND. N.LE.ORD9) THEN
                        J = 3
                    ELSE
                        J = 2
                    ENDIF
                    CH3 = STR(:J)
     
***************
*       If the Replacement is like 2p(2) = 3p(2)
*
                ELSE
                    CH3 = STR(LEFT+1:RIGHT-1)
                    K = ICTOI(CH3)
                    CH3 = STR(:LEFT-1)
                ENDIF
     
***************
*       Convert uppercase to lowercase, and assign value to ELi,Qi
*
                N = ICHAR(CH3(2:2))
                IF (N.GE.ORDUA .AND. N.LE.ORDUZ)
     :              CH3(2:2) = CHAR(N-ORDUA+ORDLA)
                EL(I) = CH3
                Q(I) = K
                IF (LEFT .EQ. 0) THEN
                    CALL STRSH(STR,(J+2))
                ELSE
                    CALL STRSH(STR,(RIGHT+1))
                ENDIF
10          CONTINUE
            RETURN
        END
*
* ----------------------------------------------------------------------
*       SUBROUTINE      R E P L A C E
* ----------------------------------------------------------------------
*
        SUBROUTINE REPLAC(ELB,QB,MB,*)
     
*   OUTPUT :
*       ELB,QB,MB
     
        PARAMETER      (NELS=15,NSHEL=5,NCOUPL=2*NSHEL-1)
        PARAMETER      (NCFG=500,NSCOUP=500)
        CHARACTER*3    EL(NELS),ELL(NELS),ELR(NELS),ELS(NELS),
     :                 ELA(NELS)
        CHARACTER*3    ELB(NELS),ELC(NELS),CH3
        CHARACTER      FBETA*8, FILE1*32, FILE2*40, FILE3*32
        INTEGER        Q(NELS),QL(NELS),QR(NELS),
     :                 QS(NELS,NELS),MS(NELS),RL(NELS),
     :                 Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,
     :                 Q10,Q11,Q12,Q13,Q14,Q15
        INTEGER        QA(NELS),QB(NELS),QC(NELS),PARITY,CONST,PP(NCFG),
     :                 NP,QQ
        COMMON         NF,NR,NFTM,MAX,MIN,PARITY,CONST,NQ
     :                 /BLK3/EL,ELL,ELR,ELS,ELA
     :                 /BLK4/Q,QL,ML,QR,MR,M,QS,MS,MA,RL,NREF,
     :                       Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,
     :                       Q10,Q11,Q12,Q13,Q14,Q15
        COMMON         /FILES/FBETA(4,NSCOUP),FILE1(NCFG),FILE2(NCFG),
     :                        FILE3(NSCOUP)
     
     
            DO 10 I=1,NELS
                ELC(I) = EL(I)
                QC(I) = Q(I)
10          CONTINUE
     
*
**********     Correct Q(i) by subtraction
*
50          DO 16 I=1,ML
                CH3 = ELL(I)
                MARK = 1
                DO 16 J=1,M
                IF (CH3 .EQ. ELC(J)) THEN
                    QC(J) = QC(J)-QL(I)
                    IF (QC(J) .GE. 0) MARK = 0
                ENDIF
16          CONTINUE
51          IF (MARK.NE.0) RETURN 1
     
*
**********     Correct QC(i) by adding
*
            MC = M
            DO 15 I=1,MR
                CH3 = ELR(I)
                MARK = 0
                DO 17 J=1,MC
                    IF (CH3 .EQ. ELC(J)) THEN
                        QC(J) = QC(J)+QR(I)
                        MARK = 1
                    ENDIF
17              CONTINUE
52              IF (MARK .EQ. 0) THEN
                    MC = MC+1
                    ELC(MC) = CH3
                    QC(MC) = QR(I)
                ENDIF
15          CONTINUE
     
***************
*       Delete EL(i) if Q(i) = 0
*
65          MB = 0
            DO 30 I=1,MC
                IF (QC(I).NE.0 .OR. MC.EQ.2) THEN
                    MB = MB+1
                    ELB(MB) = ELC(I)
                    QB(MB) = QC(I)
                ENDIF
30          CONTINUE
     
*
**********     Check the input error after replacement
*
            J = 0
            K = 0
            DO 27 I=1,MB
                CH3 = ELB(I)
                L = LVAL(CH3(2:2))
                IF (QB(I) .GT. L*4+2) THEN
                    RETURN 1
                ENDIF
                J = J+QB(I)
                K = K+QB(I)*L
27          CONTINUE
     
53          IF (J .NE. CONST) RETURN 1
     
            IF (MOD(K,2) .NE. PARITY) RETURN 1
     
*
**********    If the replacement duplicates a configuration in the
*           active set, it should not be sent to CI.LST
*
            DO 14 I=1,NQ
                READ (FILE1(I),25) (QA(J), J=1,MA)
25              FORMAT (15(I2))
                N = 0
                DO 35 J=1,MB
                    CH3 = ELB(J)
                    DO 35 K=1,MA
                        IF (CH3.EQ.ELA(K) .AND. QB(K).EQ.QA(K)) N = N+1
35              CONTINUE
                IF (N .EQ. MB) RETURN 1
14          CONTINUE
     
***************
*       If the replacement duplicates a configuration in the
*   Reference Set, it should not be sent to CI.LST
*
            DO 31 I=1,NREF
                L = MS(I)
                N = 0
                DO 36 J=1,L
                    CH3 = ELS(J)
                    QQ = QS(J,I)
                    DO 36 K=1,MB
                        IF (CH3.EQ.ELB(K) .AND. QQ.EQ.QB(K)) N = N+1
36              CONTINUE
                IF (N .EQ. MB) RETURN 1
31          CONTINUE
     
***************
*       If the replacement duplicates a configuration in the previous
*   replacement, it should not sent to CI.LST
*
            DO 40 I=1,NP
                L = PP(I)
                READ (FILE2(I),42) (ELC(J),QC(J), J=1,L)
42              FORMAT (8(A3,I2))
                N = 0
                DO 41 J=1,L
                    CH3 = ELC(J)
                    DO 41 K=1,MB
                        IF (CH3.EQ.ELB(K) .AND. QC(K).EQ.QB(K)) N = N+1
41              CONTINUE
                IF (N .EQ. MB) RETURN 1
40          CONTINUE
     
            NP = NP+1
            PP(NP) = MB
            WRITE (FILE2(NP),37) (ELB(J),QB(J), J=1,MB)
37          FORMAT (8(A3,I2))
            RETURN
        END
*
* ----------------------------------------------------------------------
*       SUBROUTINE      V P A I R
* ----------------------------------------------------------------------
*
*       Generate occupied or virtual pairs for D-Replacement
*
        SUBROUTINE VPAIR (ELV,MV,PL,LLMIN,LLMAX,STR,*)
     
*   INPUT :
*       ELV  =  ELi for Virtual Set
*        MV  =  Number of ELi for Virtual Set
*        PL  =  Parity of Qi for Reference Set
*     LLMIN  =  Minimum angular coupling of the pair
*     LLMAX  =  Maximum angular coupling of the pair
*       STR  =  String to be packed as output for Replacement
*
     
     
        PARAMETER     (NELS=15,NSHEL=5,NCOUPL=2*NSHEL-1)
        PARAMETER     (NCFG=500,NSCOUP=500)
        CHARACTER*72   HEADER,SHELLS,VIRTUL,STR,SSTR
        CHARACTER*60   REF(NELS),ACT,REPL(NELS)
        CHARACTER*2    FTM(NELS)
        CHARACTER*3    EL(NELS),ELL(NELS),ELR(NELS),ELS(NELS),
     :                 ELA(NELS)
        CHARACTER*3    ELV(NELS),ELB(NELS),ELC(NELS)
        INTEGER        Q(NELS),QL(NELS),QR(NELS),QB(NELS),
     :                 QS(NELS,NELS),MS(NELS),RL(NELS),
     :                 Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,
     :                 Q10,Q11,Q12,Q13,Q14,Q15
        INTEGER        PL,PR,PARITY,CONST
        COMMON         NF,NR,NFTM,MAX,MIN,PARITY,CONST,NQ
     :                 /BLK1/HEADER,SHELLS,ACT,VIRTUL
     :                 /BLK2/REF,REPL,FTM
     :                 /BLK3/EL,ELL,ELR,ELS,ELA
     :                 /BLK4/Q,QL,ML,QR,MR,M,QS,MS,MA,RL,NREF,
     :                       Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,
     :                       Q10,Q11,Q12,Q13,Q14,Q15
     
            N = INDEX(STR,'=')
            NR = NELS
     
***************
*       D-Replacement for the pair of two single ELi
*
            MR = 2
            QR(1) = 1
            QR(2) = 1
            DO 10 I = 1,MV-1
                ELR(1) = ELV(I)
		LR1 = LVAL(ELR(1)(2:2))
                DO 10 J=I+1,MV
                    ELR(2) = ELV(J)
		    LR2 = LVAL(ELR(2)(2:2))
		    LRMIN = IABS(LR1 - LR2)
		    LRMAX = LR1 + LR2
                    PR = MOD(LRMAX, 2)
*
*       If the pair has the same parity with the left side, and the
*   angular coupling of the two pairs have values in common, replace 
*   them, then generate couplings for the new configuration
*
                    IF (PL .NE. PR .OR. LRMIN .GT. LLMAX .OR.
     :                                  LRMAX .LT. LLMIN ) GOTO 10
                    CALL REPLAC(ELB,QB,K,*10)
                    CALL COUPLD (ELB,QB,K,ELC,NC,*25)
                    IF (NC .GT. 0) THEN
                        WRITE (6,12)
12                      FORMAT (//'   FOR VIRTUAL SET, GENERATE',
     :  ' CONFIGURATION AND COUPLINGS FOR D-REPLACEMENT'//)
                        SSTR = STR(:N)//ELR(1)//'.'//ELR(2)
                        CALL DEL(SSTR)
                        REPL(NR) = SSTR
                        CALL PRINT (ELC,QB,K,NC,4)
                    ENDIF
10          CONTINUE
     
***************
*       D-Replacement for the pairs which has the value Qi=2
*
50          MR = 1
            QR(1) = 2
            DO 20 I=1,MV
                ELR(1) = ELV(I)
		LR1 = LVAL(ELR(1)(2:2))
		LRMAX = LR1 + LR1
                PR = MOD(LRMAX, 2)
*
*       If it has the same parity with the left side, replace them,
*   then generate couplings for the new configuration
*
                IF (PL .EQ. PR .AND. LRMAX .GE. LLMIN) THEN
                    CALL REPLAC(ELB,QB,K,*20)
                    CALL COUPLD (ELB,QB,K,ELC,NC,*25)
                    IF (NC .GT. 0) THEN
                        WRITE (6,22)
22                      FORMAT (//'   FOR VIRTUAL SET, GENERATE ',
     : ' CONFIGURATION AND COUPLINGS FOR D-REPLACEMENT'//)
                        SSTR = STR(:N)//ELR(1)//'(2)'
                       CALL DEL(SSTR)
                       REPL(NR) = SSTR
                        CALL PRINT (ELC,QB,K,NC,4)
                    ENDIF
                ENDIF
20          CONTINUE
     
            RETURN
25          RETURN 1
        END
*
* ----------------------------------------------------------------------
*               SUBROUTINE      P R I N T
* ----------------------------------------------------------------------
*
*       Print out the values of couplings
*
        SUBROUTINE PRINT (EL,Q,M,NC,MARK)
            PARAMETER      (NELS=15,NSHEL=5,NCOUPL=2*NSHEL-1)
            PARAMETER      (NCFG=500,NSCOUP=500)
            CHARACTER*60   REF(NELS),REPL(NELS),ACT
            CHARACTER*72   HEADER,SHELLS,VIRTUL
            CHARACTER*3    EL(NELS),COUPLE(NCOUPL)
            CHARACTER*2    FTM(NELS)
            CHARACTER      CH1
            CHARACTER      FBETA*8, FILE1*32, FILE2*40, FILE3*32
            INTEGER        Q(NELS),PARITY,CONST
            COMMON         NF,NR,NFTM,MAX,MIN,PARITY,CONST,NQ
     :                     /BLK1/HEADER,SHELLS,ACT,VIRTUL
     :                     /BLK2/REF,REPL,FTM
            COMMON       /FILES/FBETA(4,NSCOUP),FILE1(NCFG),FILE2(NCFG),
     :                        FILE3(NSCOUP)
     
     
*   MARK  =  1 for Reference Set
*         =  2 for Reference Set, Active Set
*         =  3 for Reference Set, Replacement
*         =  4 for Reference Set, Virtual Set
*
***************
*               Print the input of Header
*
            J = INDEX(HEADER,'          ')
            WRITE (6,11) HEADER(1:J)
11          FORMAT (' '/T10,'-------------       ',A,'  --------',/)
     
***************
*               Print the input of Closed Shells
*
            J = INDEX(SHELLS,'     ')
            WRITE (6,12) SHELLS(:J)
12          FORMAT ('          Closed Shells  :  ',A/)
            IF (MARK.EQ.1 .OR. MARK.EQ.2) GOTO 20
     
***************
*               Print the input of Reference Set
*
            WRITE (6,10) REF(NF)
10          FORMAT ('          Reference set  :  ',A60/)
            IF (MARK .EQ. 3) GOTO 21
     
***************
*               Print the input of Virtual Set
*
            IF (MARK .EQ. 4) CH1 = 'S'
            IF (MARK .EQ. 5) CH1 = 'D'
            IF (MARK.EQ.4 .OR. MARK.EQ.5) THEN
                WRITE (6,30) VIRTUL
30              FORMAT ('            Virtual Set  :  ',A60/)
                WRITE (6,31) CH1,REPL(NELS)
31              FORMAT ('          ',A1,'-Replacement  :  ',A60/)
                GOTO 20
            ENDIF
     
***************
*               Print the input of Active Set
*
            WRITE (6,19) ACT
19          FORMAT ('             Active Set  :  ',A60/)
            IF (MARK .EQ. 2) GOTO 20
     
***************
*               Print the input of Replacement
*
21          WRITE (6,22) REPL(NR)
22          FORMAT ('            Replacement  :  ',A60/)
     
***************
*       Print the new configuration by Replacement
*
20          IF (EL(1)(1:1) .EQ. ' ') THEN
                K = 1
            ELSE
                K = 2
            ENDIF
            WRITE (6,13) (EL(J),Q(J), J=1,M)
13          FORMAT ('          Configuration  :',
     :          5(2X,A3,'(',I2,')')/)
     
***************
*       Print couplings generated from the configuration
*
            K = 2*M-1
            IF (M .LE. 3) THEN
                READ (FILE3(1),17) (COUPLE(J), J=1,K)
                WRITE (6,23) (COUPLE(J), J=1,K)
23              FORMAT ('        Their couplings  :  ',9(5X,A3))
                DO 24 I=2,NC
                    READ (FILE3(I),17) (COUPLE(J), J=1,K)
                    WRITE (6,25) (COUPLE(J), J=1,K)
25                  FORMAT (T29,9(5X,A3))
24              CONTINUE
            ELSE
                WRITE (6,15)
15              FORMAT ('         Their couplings  :'/)
                DO 16 I=1,NC
                    READ (FILE3(I),17) (COUPLE(J), J=1,K)
17                  FORMAT (9(A3))
                    WRITE (6,18) (COUPLE(J), J=1,K)
18                  FORMAT (T5,9(5X,A3))
16              CONTINUE
            ENDIF
            RETURN
        END
*
* ----------------------------------------------------------------------
*       SUBROUTINE      H E L P
* ----------------------------------------------------------------------
*
*       Explanation of the input format
*
        SUBROUTINE HELP ()
            CHARACTER*10   STR
     
10          WRITE(0,11)
11          FORMAT(//,
     : 5X,'This program prompts for each required input.  The user'/
     : 5X,'should enter data or a RETURN after a question (?) mark')
            WRITE(0,*)
            WRITE(0,*) '     Example 1 :'
            WRITE(0,*) '    --------------------'
            WRITE(0,*) '                 Header  ?  S II ground state'
            WRITE(0,*) '          Closed shells  ?   1s 2s 2p'
            WRITE(0,*) '          Reference Set  ?  3s(2) 3p(3)'
            WRITE(0,*) '                      2  ?  RETURN'
            WRITE(0,*) '             Active Set  ?  3s,3p'
            WRITE(0,*) 'Type of set generation   ?  0'
            WRITE(0,*) '            Replacement  ?  3s(2) = 3d(2)'
            WRITE(0,*) '                      2  ?  3s = 3d'
            WRITE(0,*) '                      3  ?  3s.3p = 4s.3d'
            WRITE(0,*) '                      4  ?  <RETURN>'
            WRITE(0,*) '             Final Term  ?  4S'
            WRITE(0,*) '                      2  ?  RETURN'
            WRITE(0,12)
12          FORMAT (/
     : 5X,'Header and Closed Shells cannot exceed 72 characters and'/
     : 5X,'will be copied to the output file. The electrons are'/
     : 5X,'separated by blanks in the Closed Shells.')
	    WRITE(0,13)
13          FORMAT(5X,'Press RETURN for more... ')
            READ(5,'(A)') STR
     
20          WRITE(0,21)
21          FORMAT(///////
     : 5X,'Items are separated by a blank in the Reference Set, by a'/
     : 5X,'comma or a blank in the Active set, and by a period or a'/
     : 5X,'blank in Replacements.'//
     : 5X,'Reference Set, Replacement, and Final Term are three sets'/
     : 5X,'of input, each with 0 to 10 members.  Each member must be'/
     : 5X,'entered on a separate line.')

            WRITE(0,*) '       PRINT RETURN to terminate the input set.'
            WRITE(0,*) '       PRINT RETURN if the set is empty.'
            WRITE(0,*)
            WRITE(0,*) '     Example 2 :'
            WRITE(0,*) '    --------------------'
            WRITE(0,*)
            WRITE(0,*) '         Reference Set  ?  2s(1) 2p(2) 3s(1)'
            WRITE(0,*) '                     2  ?  RETURN'
            WRITE(0,*) '            Active Set  ?  2s,2p,3s'
            WRITE(0,*) 'Type of set generation  ?  0'
            WRITE(0,*) '           Replacement  ?  RETURN '
            WRITE(0,*) '            Final Term  ?  RETURN'
            WRITE(0,*)
            WRITE(0,*) 
     : '     Where the Replacement and the Final Term are empty.'
            WRITE(0,*)
            WRITE(0,*)
            WRITE(0,22)
22          FORMAT(5X,'Press ''b'' for BACK or RETURN for more... ')
            READ(5,'(A)') STR
            I = INDEX(STR,'B')
            J = INDEX(STR,'b')
            IF (I.NE.0 .OR. J.NE.0) GOTO 10
     
30          WRITE(0,31)
31          FORMAT(///////
     : 5X,'By inputing "s" or "d" or "sd" you can compute the config-'/
     : 5X,'urations from the Virtual Set, where  S means Single'/
     : 5X,'Replacement, D means Double Replacement, SD means Single'/
     : 5X,'and Double Replacement.'//
     : 5X,'GENCL will give you prompts for the Virtual Set automati-'/
     : 5X,'cally, then you need to specify the range of shells that'/
     : 5X,'are to participate in the replacements. For instance, a'/
     : 5X,'response of 2 to "From which shell" and of 3 to "To which'/
     : 5X,'shell" implies that shells 2 and 3 participates in the'/
     : 5X,'replacements, and shell 1 does not enter into any',
     : 1X,'replacements.')
            WRITE(0,*)
            WRITE(0,*) '     Example 3  :'
            WRITE(0,*) '    -------------------------'
            WRITE(0,*) '                ...'
            WRITE(0,*) '         Reference Set  ?  2s(1) 2p(1) 3s(1)'
            WRITE(0,*) '            Active Set  ?  RETURN'
            WRITE(0,*) '           Replacement  ?  sd'
            WRITE(0,*) '           Virtual Set  ?  3p,3d,4s'
            WRITE(0,*) '      From which shell  ?  2'
            WRITE(0,*) '        To which shell  ?  3'
            WRITE(0,*) '            Final Term  ?  RETURN'
            WRITE(0,*)
            WRITE(0,*)
            WRITE(0,*)
            WRITE(0,22)
            READ(5,'(A)') STR
            I = INDEX(STR,'B')
            J = INDEX(STR,'b')
            IF (I.NE.0 .OR. J.NE.0) GOTO 20
40          WRITE(0,41)
41          FORMAT(///////
     : 5X,'After terminating an input line, if you find the previous'/
     : 5X,'input to be wrong, type ''B'' or ''b'' to go back one step.'
     :/5X,'For example, before inputing Active Set, if you find the'/
     : 5X,'wrong spelling in the Header, type "B" and GENCL will '/
     : 5X,'prompt for the Header again.'//)
            WRITE(0,*) '     Example 4 :'
            WRITE(0,*) '    ----------------------'
            WRITE(0,*) '                Header  ?  OXYYEN'
            WRITE(0,*) '            Active Set  ?  B'
            WRITE(0,*) '                Header  ?  OXYGEN '
            WRITE(0,*) '            Active Set  ?  2s '
            WRITE(0,*)
            WRITE(0,*) '     Then the following prompts will continue.'
            WRITE(0,*)
            WRITE(0,*)
            WRITE(0,22)
            READ(5,'(A)') STR
            I = INDEX(STR,'B')
            J = INDEX(STR,'b')
            IF (I.NE.0 .OR. J.NE.0) GOTO 30
     
50          WRITE(0,51)
51          FORMAT(///////
     : 5X,'Example 4 shows the procedure for going back four steps'/
     : 5X,'to correct the Closed Shells from  5s  to 1s  2s.'//)
            WRITE(0,*) '     Example 5 :'
            WRITE(0,*) '    -----------------------'
            WRITE(0,*)
            WRITE(0,*) '            Active Set  ?  5s'
            WRITE(0,*) 'Type of set generation  ?  0'
            WRITE(0,*) '         Reference Set  ?  2s(1) 2p(2) 3s(1)'
            WRITE(0,*) '                     2  ?  2P(4)'
            WRITE(0,*) '                     3  ?  b '
            WRITE(0,*) '                     2  ?  b '
            WRITE(0,*) '         Reference Set  ?  b '
            WRITE(0,*) '         Closed Shells  ?    1s  2s '
            WRITE(0,*) '         Reference Set  ?  '
            WRITE(0,*)
            WRITE(0,*)
            WRITE(0,*) '   Then reenter the data for the Reference Set',
	    
     :                 ' and continue the input.'
            WRITE(0,*)
            WRITE(0,22)
            READ(5,'(A)') STR
            I = INDEX(STR,'B')
            J = INDEX(STR,'b')
            IF (I.NE.0 .OR. J.NE.0) GOTO 40
     
60          WRITE(0,61)
61          FORMAT(///////
     : 5X,'When the following error conditions are detected,',
     :    ' a message is given.'/
     :10X,' 1). The parentheses are not matched'/
     :10X,' 2). The number of electrons in a shell is more than FULL'//
     :15X,'For each member of the Reference Set ,'/
     :10X,' 3). The number of electrons is not the same '/
     :10X,' 4). The parity is not the same '//
     :10X,' 5). The number of couplings generated by a configuration'/
     :10X,'     is more than 500'///////)
            WRITE(0,65)
65          FORMAT (5X,'Press ''b'' for BACK or RETURN to begin the',
     :              1X,'program.')
            READ(5,'(A)') STR
            I = INDEX(STR,'B')
            J = INDEX(STR,'b')
            IF (I.NE.0 .OR. J.NE.0) GOTO 50
            DO 18 I=1,30
                WRITE(0,*)
18          CONTINUE
            FLAG = 1
            RETURN
            END
