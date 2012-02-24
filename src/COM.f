*    
*     Routines for MCHF_LIB_COM
*     Computer Physics Communication, Vol. 64, 399-405 (1991)
*
*     C O P Y R I G H T -- 1994
*
*-----------------------------------------------------------------------
*       A C U R A T
*-----------------------------------------------------------------------
*
*       Coefficients of Slater integrals are square roots of
*  rational numbers. To improve the accuracy, certain commonly
*  occurring coefficients are improved to machine precision.
*
        DOUBLE PRECISION FUNCTION ACURAT(C)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        INTEGER NUM
        INTEGER DEN(11)
        DATA DEN/2,3,7,9,15,35,49,175,189,315,441/
        DATA D1/1.D0/
*
        C2 = C*C
        ACURAT = C
        DO 1 I = 1,11
           PROD = DEN(I)*C2
           NUM = NINT(PROD)
           EPS = ABS(NUM-PROD)/DEN(I)
           IF (EPS .LE. 1.E-8) THEN
              IF (EPS .NE. 0.) THEN
                 ACURAT = SQRT((NUM*D1)/DEN(I))
                 IF (C .LT. 0.) ACURAT = -ACURAT
                 RETURN
              END IF
           END IF
  1     CONTINUE
        END
*
*     -------------------------------------------------------------
*        B I S E C T
*     -------------------------------------------------------------
*
*
      SUBROUTINE BISECT(N,EPS1,D,E,E2,LB,UB,MM,M,W,IND,IERR,RV4,RV5) 
* 
      INTEGER I,J,K,L,M,N,P,Q,R,S,II,MM,M1,M2,TAG,IERR,ISTURM 
      DOUBLE PRECISION D(N),E(N),E2(N),W(MM),RV4(N),RV5(N) 
      DOUBLE PRECISION U,V,LB,T1,T2,UB,XU,X0,X1,EPS1,MACHEP 
      DOUBLE PRECISION DABS,DMAX1,DMIN1,DFLOAT 
      INTEGER IND(MM) 
* 
*     THIS SUBROUTINE IS A TRANSLATION OF THE BISECTION TECHNIQUE 
*     IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON. 
*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971). 
* 
*     THIS SUBROUTINE FINDS THOSE EIGENVALUES OF A TRIDIAGONAL 
*     SYMMETRIC MATRIX WHICH LIE IN A SPECIFIED INTERVAL, 
*     USING BISECTION. 
* 
*     ON INPUT: 
* 
*        N IS THE ORDER OF THE MATRIX; 
* 
*        EPS1 IS AN ABSOLUTE ERROR TOLERANCE FOR THE COMPUTED 
*          EIGENVALUES.  IF THE INPUT EPS1 IS NON-POSITIVE, 
*          IT IS RESET FOR EACH SUBMATRIX TO A DEFAULT VALUE, 
*          NAMELY, MINUS THE PRODUCT OF THE RELATIVE MACHINE 
*          PRECISION AND THE 1-NORM OF THE SUBMATRIX; 
* 
*        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX; 
* 
*        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX 
*          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY; 
* 
*        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E. 
*          E2(1) IS ARBITRARY; 
* 
*        LB AND UB DEFINE THE INTERVAL TO BE SEARCHED FOR EIGENVALUES. 
*          IF LB IS NOT LESS THAN UB, NO EIGENVALUES WILL BE FOUND; 
* 
*        MM SHOULD BE SET TO AN UPPER BOUND FOR THE NUMBER OF 
*          EIGENVALUES IN THE INTERVAL.  WARNING: IF MORE THAN 
*          MM EIGENVALUES ARE DETERMINED TO LIE IN THE INTERVAL, 
*          AN ERROR RETURN IS MADE WITH NO EIGENVALUES FOUND. 
* 
*     ON OUTPUT: 
* 
*        EPS1 IS UNALTERED UNLESS IT HAS BEEN RESET TO ITS 
*          (LAST) DEFAULT VALUE; 
* 
*        D AND E ARE UNALTERED; 
* 
*        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED 
*          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE 
*          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES. 
*          E2(1) IS ALSO SET TO ZERO; 
* 
*        M IS THE NUMBER OF EIGENVALUES DETERMINED TO LIE IN (LB,UB); 
* 
*        W CONTAINS THE M EIGENVALUES IN ASCENDING ORDER; 
* 
*        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES 
*          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W -- 
*          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM 
*          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.; 
* 
*        IERR IS SET TO 
*          ZERO       FOR NORMAL RETURN, 
*          3*N+1      IF M EXCEEDS MM; 
* 
*        RV4 AND RV5 ARE TEMPORARY STORAGE ARRAYS. 
* 
*     THE ALGOL PROCEDURE STURMCNT CONTAINED IN TRISTURM 
*     APPEARS IN BISECT IN-LINE. 
* 
*     NOTE THAT SUBROUTINE TQL1 OR IMTQL1 IS GENERALLY FASTER THAN 
*     BISECT, IF MORE THAN N/4 EIGENVALUES ARE TO BE FOUND. 
* 
*     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW, 
*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY 
* 
*     ------------------------------------------------------------------ 
* 
*     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING 
*                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC. 
*                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC 
*                ON S360 :::::::::: 
      DATA MACHEP/1.D-12/ 
* 
      IERR = 0 
      TAG = 0 
      T1 = LB 
      T2 = UB 
*     :::::::::: LOOK FOR SMALL SUB-DIAGONAL ENTRIES :::::::::: 
      DO 40 I = 1, N 
         IF (I .EQ. 1) GO TO 20 
         IF (DABS(E(I)) .GT. MACHEP * (DABS(D(I)) + DABS(D(I-1)))) 
     :      GO TO 40 
   20    E2(I) = 0.0D0 
   40 CONTINUE 
*     :::::::::: DETERMINE THE NUMBER OF EIGENVALUES 
*                IN THE INTERVAL :::::::::: 
      P = 1 
      Q = N 
      X1 = UB 
      ISTURM = 1 
      GO TO 320 
   60 M = S 
      X1 = LB 
      ISTURM = 2 
      GO TO 320 
   80 M = M - S 
      IF (M .GT. MM) GO TO 980 
      Q = 0 
      R = 0 
*     :::::::::: ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING 
*                INTERVAL BY THE GERSCHGORIN BOUNDS :::::::::: 
  100 IF (R .EQ. M) GO TO 1001 
      TAG = TAG + 1 
      P = Q + 1 
      XU = D(P) 
      X0 = D(P) 
      U = 0.0D0 
* 
      DO 120 Q = P, N 
         X1 = U 
         U = 0.0D0 
         V = 0.0D0 
         IF (Q .EQ. N) GO TO 110 
         U = DABS(E(Q+1)) 
         V = E2(Q+1) 
  110    XU = DMIN1(D(Q)-(X1+U),XU) 
         X0 = DMAX1(D(Q)+(X1+U),X0) 
         IF (V .EQ. 0.0D0) GO TO 140 
  120 CONTINUE 
* 
  140 X1 = DMAX1(DABS(XU),DABS(X0)) * MACHEP 
      IF (EPS1 .LE. 0.0D0) EPS1 = -X1 
      IF (P .NE. Q) GO TO 180 
*     :::::::::: CHECK FOR ISOLATED ROOT WITHIN INTERVAL :::::::::: 
      IF (T1 .GT. D(P) .OR. D(P) .GE. T2) GO TO 940 
      M1 = P 
      M2 = P 
      RV5(P) = D(P) 
      GO TO 900 
  180 X1 = X1 * DFLOAT(Q-P+1) 
      LB = DMAX1(T1,XU-X1) 
      UB = DMIN1(T2,X0+X1) 
      X1 = LB 
      ISTURM = 3 
      GO TO 320 
  200 M1 = S + 1 
      X1 = UB 
      ISTURM = 4 
      GO TO 320 
  220 M2 = S 
      IF (M1 .GT. M2) GO TO 940 
*     :::::::::: FIND ROOTS BY BISECTION :::::::::: 
      X0 = UB 
      ISTURM = 5 
* 
      DO 240 I = M1, M2 
         RV5(I) = UB 
         RV4(I) = LB 
  240 CONTINUE 
*     :::::::::: LOOP FOR K-TH EIGENVALUE 
*                FOR K=M2 STEP -1 UNTIL M1 DO -- 
*                (-DO- NOT USED TO LEGALIZE COMPUTED-GO-TO) :::::::::: 
      K = M2 
  250    XU = LB 
*     :::::::::: FOR I=K STEP -1 UNTIL M1 DO -- :::::::::: 
         DO 260 II = M1, K 
            I = M1 + K - II 
            IF (XU .GE. RV4(I)) GO TO 260 
            XU = RV4(I) 
            GO TO 280 
  260    CONTINUE 
* 
  280    IF (X0 .GT. RV5(K)) X0 = RV5(K) 
*     :::::::::: NEXT BISECTION STEP :::::::::: 
  300    X1 = (XU + X0) * 0.5D0 
         IF ((X0 - XU) .LE. (2.0D0 * MACHEP * 
     :      (DABS(XU) + DABS(X0)) + DABS(EPS1))) GO TO 420 
*     :::::::::: IN-LINE PROCEDURE FOR STURM SEQUENCE :::::::::: 
  320    S = P - 1 
         U = 1.0D0 
* 
         DO 340 I = P, Q 
            IF (U .NE. 0.0D0) GO TO 325 
            V = DABS(E(I)) / MACHEP 
            GO TO 330 
  325       V = E2(I) / U 
  330       U = D(I) - X1 - V 
            IF (U .LT. 0.0D0) S = S + 1 
  340    CONTINUE 
* 
         GO TO (60,80,200,220,360), ISTURM 
*     :::::::::: REFINE INTERVALS :::::::::: 
  360    IF (S .GE. K) GO TO 400 
         XU = X1 
         IF (S .GE. M1) GO TO 380 
         RV4(M1) = X1 
         GO TO 300 
  380    RV4(S+1) = X1 
         IF (RV5(S) .GT. X1) RV5(S) = X1 
         GO TO 300 
  400    X0 = X1 
         GO TO 300 
*     :::::::::: K-TH EIGENVALUE FOUND :::::::::: 
  420    RV5(K) = X1 
      K = K - 1 
      IF (K .GE. M1) GO TO 250 
*     :::::::::: ORDER EIGENVALUES TAGGED WITH THEIR 
*                SUBMATRIX ASSOCIATIONS :::::::::: 
  900 S = R 
      R = R + M2 - M1 + 1 
      J = 1 
      K = M1 
* 
      DO 920 L = 1, R 
         IF (J .GT. S) GO TO 910 
         IF (K .GT. M2) GO TO 940 
         IF (RV5(K) .GE. W(L)) GO TO 915 
* 
         DO 905 II = J, S 
            I = L + S - II 
            W(I+1) = W(I) 
            IND(I+1) = IND(I) 
  905    CONTINUE 
* 
  910    W(L) = RV5(K) 
         IND(L) = TAG 
         K = K + 1 
         GO TO 920 
  915    J = J + 1 
  920 CONTINUE 
* 
  940 IF (Q .LT. N) GO TO 100 
      GO TO 1001 
*     :::::::::: SET ERROR -- UNDERESTIMATE OF NUMBER OF 
*                EIGENVALUES IN INTERVAL :::::::::: 
  980 IERR = 3 * N + 1 
 1001 LB = T1 
      UB = T2 
      RETURN 
*     :::::::::: LAST CARD OF BISECT :::::::::: 
      END 
*
*     ------------------------------------------------------------------
*       B W I N T
*     ------------------------------------------------------------------
*
      SUBROUTINE BWINT(LC,LO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/BLUME/COEFN2(4),COEFNK(4),COEFVK(4)
*
* ... LC IS THE L-VALUE OF THE FILLED SUBSHELL, LO IS THE L-VALUE
*     OF THE PARTIALLY-FILLED SUBSHELL.
*
      IF(LC.LE.3.AND.LO.LE.4) GO TO 1
      PRINT 100, LC,LO
  100 FORMAT (37H INCORRECT CALLING OF BWINT WITH LC =,I2,6H, LO =,I2)
    1 LC1 = LC + 1
      GO TO (10,20,30,40), LC1
   10 GO TO (11,12,13,14), LO
*
* ... S-P
*
   11 COEFNK(1) = 1.D0
      COEFN2(1) = -2.D0
      COEFVK(1) = 1.D0
      RETURN
*
* ... S-D
*
   12 COEFNK(1) = 6.D0/5.D0
      COEFN2(1) = -9.D0/5.D0
      COEFVK(1) = 3.D0/5.D0
      RETURN
*
* ... S-F
*
   13 COEFNK(1) = 9.D0/7.D0
      COEFN2(1) = -12.D0/7.D0
      COEFVK(1) = 3.D0/7.D0
      RETURN
*
* ... S-G
*
   14 COEFNK(1) = 4.D0/3.D0
      COEFN2(1) = -5.D0/3.D0
      COEFVK(1) = 1.D0/3.D0
      RETURN
   20 GO TO (21,22,23,24), LO
*
* ... P-P
*
   21 COEFNK(1) = 0.D0
      COEFN2(1) = 3.D0
      COEFVK(1) = 9.D0/5.D0
      RETURN
*
* ... P-D
*
   22 COEFNK(1) = 3.D0/7.D0
      COEFNK(2) = 36.D0/35.D0
      COEFN2(1) = -12.D0/5.D0
      COEFN2(2) = 0.D0
      COEFVK(1) = 3.D0/5.D0
      COEFVK(2) = 36.D0/35.D0
      RETURN
*
* ... P-F
*
   23 COEFNK(1) = 1.D0/7.D0
      COEFNK(2) = 10.D0/7.D0
      COEFN2(1) = -18.D0/7.D0
      COEFN2(2) = 0.D0
      COEFVK(1) = 18.D0/35.D0
      COEFVK(2) = 5.D0/7.D0
      RETURN
*
* ... P-G
*
   24 COEFNK(1) = 5.D0/77.D0
      COEFNK(2) = 18.D0/11.D0
      COEFN2(1) = -18.D0/7.D0
      COEFN2(2) = 0.D0
      COEFVK(1) = 3.D0/7.D0
      COEFVK(2) = 6.D0/11.D0
      RETURN
   30 GO TO (31,32,33,34), LO
*
* ... D-P
*
   31 COEFNK(1) = 59.D0/7.D0
      COEFNK(2) = -18.D0/7.D0
      COEFN2(1) = -4.D0
      COEFN2(2) = 0.D0
      COEFVK(1) = -1.D0
      COEFVK(2) = 18.D0/7.D0
      RETURN
*
* ... D-D
*
   32 COEFNK(1) = 6.D0/7.D0
      COEFNK(2) = 0.D0
      COEFN2(1) = 3.D0
      COEFN2(2) = 0.D0
      COEFVK(1) = 3.D0/7.D0
      COEFVK(2) = 10.D0/7.D0
      RETURN
*
* ... D-F
*
   33 COEFNK(1) = 9.D0/7.D0
      COEFNK(2) = -13.D0/77.D0
      COEFNK(3) = 75.D0/77.D0
      COEFN2(1) = -18.D0/7.D0
      COEFN2(2) = 0.D0
      COEFN2(3) = 0.D0
      COEFVK(1) = 3.D0/7.D0
      COEFVK(2) = 3.D0/7.D0
      COEFVK(3) = 75.D0/77.D0
      RETURN
*
* ... D-G
*
   34 COEFNK(1) = 741.D0/693.D0
      COEFNK(2) = -215.D0/429.D0
      COEFNK(3) = 210.D0/143.D0
      COEFN2(1) = -3.D0
      COEFN2(2) = 0.D0
      COEFN2(3) = 0.D0
      COEFVK(1) = 3.D0/7.D0
      COEFVK(2) = 255.D0/693.D0
      COEFVK(3) = 105.D0/143.D0
      RETURN
   40 GO TO (41,42,43,44), LO
*
* ... F-P
*
   41 COEFNK(1) = 52.D0/3.D0
      COEFNK(2) = -20.D0/3.D0
      COEFN2(1) = -9.D0
      COEFN2(2) = 0.D0
      COEFVK(1) = -9.D0/5.D0
      COEFVK(2) = 10.D0/3.D0
      RETURN
*
* ... F-D
*
   42 COEFNK(1) = 5.D0
      COEFNK(2) = 142.D0/55.D0
      COEFNK(3) = -20.D0/11.D0
      COEFN2(1) = -18.D0/5.D0
      COEFN2(2) = 0.D0
      COEFN2(3) = 0.D0
      COEFVK(1) = -3.D0/5.D0
      COEFVK(2) = 2.D0/5.D0
      COEFVK(3) = 20.D0/11.D0
      RETURN
*
* ... F-F
*
   43 COEFNK(1) = 1.D0
      COEFNK(2) = 5.D0/11.D0
      COEFNK(3) = 0.D0
      COEFN2(1) = 3.D0
      COEFN2(2) = 0.D0
      COEFN2(3) = 0.D0
      COEFVK(1) = 1.D0/5.D0
      COEFVK(2) = 5.D0/11.D0
      COEFVK(3) = 175.D0/143.D0
      RETURN
*
* ... F-G
*
   44 COEFNK(1) = 53.D0/33.D0
      COEFNK(2) = 57.D0/143.D0
      COEFNK(3) = -115.D0/429.D0
      COEFNK(4) = 392.D0/429.D0
      COEFN2(1) = -8.D0/3.D0
      COEFN2(2) = 0.D0
      COEFN2(3) = 0.D0
      COEFN2(4) = 0.D0
      COEFVK(1) = 1.D0/3.D0
      COEFVK(2) = 3.D0/11.D0
      COEFVK(3) = 57.D0/143.D0
      COEFVK(4) = 392.D0/429.D0
      RETURN
      END
*
*     ------------------------------------------------------------------
*          C A
*     ------------------------------------------------------------------
*
*
      DOUBLE PRECISION FUNCTION CA(L,K)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /EAV/CCA(10),CCB(35)
*
      IF (L .LE. 4) THEN
         CA = CCA((L*(L-1) + K)/2)
       ELSE
         CA = RME(L,L,K)**2
      END IF
      END
*
*     -----------------------------------------------------------------
*          C B
*     -----------------------------------------------------------------
*
*
      DOUBLE PRECISION FUNCTION CB(L,LP,K)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /EAV/CCA(10),CCB(35)
      INTEGER ICBPTR(0:4)
      DATA ICBPTR/1,6,14,23,31/
*
             IF (L .LE. LP) THEN
                 L1 = L
                 L2 = LP
              ELSE
                 L1 = LP
                 L2 = L
             END IF
             IF ( L2 .LE. 4) THEN
                CB = CCB(ICBPTR(L1)+(K+L1-L2)/2+(L1+1)*(L2-L1))
               ELSE
                CB = RME(L,LP,K)**2/(2*(2*L+1)*(2*LP+1))
             END IF
      END
* 
*     ------------------------------------------------------------------ 
*             E P T R
*     ------------------------------------------------------------------ 
* 
*       Determines the position of the electron in the electron list 
* 
      SUBROUTINE EPTR(EL,ELSYMB, IEL, *) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      PARAMETER(IWRITE=6)
      CHARACTER CONFIG*66, COUPLE*3, EL(*)*3, ELSYMB*3, BL*3 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      DATA BL/'   '/ 
* 
* ***** SEARCH ELECTRON LIST FOR LSYMB 
* 
      IF ( ELSYMB .EQ. BL ) THEN 
         IEL = 0 
         RETURN 
      ENDIF 
      DO 10 I=1,NWF 
         IF (EL(I) .EQ. ELSYMB ) THEN 
            IEL = I 
            RETURN 
         ENDIF 
10    CONTINUE 
      IEL = -1 
      WRITE (IWRITE,20) ELSYMB 
20    FORMAT(/10X,A3,' NOT FOUND IN ELECTRON LIST') 
      RETURN 1 
      END 
*
*     -----------------------------------------------------------------
*           F A C T R L
*     -----------------------------------------------------------------
*
*
      SUBROUTINE FACTRL(NFACT) 
* 
*      GAM(I) = LOG( GAMMA(I-1) ), WHERE GAMMA(I) = FACTORIAL I-1 
* 
      IMPLICIT REAL *8(A-H,O-Z) 
* 
      COMMON/FACT/GAM(100) 
      DATA ZERO,ONE,TWO/0.D0,1.D0,2.D0/ 
* 
      GAMMA=ONE 
      GAM(1) = ZERO 
      DO 1 I=1,NFACT-1 
         GAMMA=I*GAMMA 
         GAM(I+1) = DLOG(GAMMA) 
    1 CONTINUE 
      DO 20 I = NFACT+1,(100) 
         X = I-1 
         GAM(I) = GAM(I-1) + DLOG(X) 
   20 CONTINUE 
      RETURN 
      END 
* 
*     ------------------------------------------------------------------ 
*              F K 
*     ------------------------------------------------------------------ 
*                             k 
*       Returns the value of F (i,j) 
* 
* 
      DOUBLE PRECISION FUNCTION FK(I,J,K,REL) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      LOGICAL REL
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      CALL YKF(I,I,K,REL) 
      FK = QUADS(J,J,1) 
      IF (MASS .EQ. 2) FK = FK*(D1 + RMASS/D2)
      RETURN 
      END 
* 
*     ------------------------------------------------------------------ 
*              G K 
*     ------------------------------------------------------------------ 
*                             k 
*       Returns the value of G (i,j). 
* 
* 
      DOUBLE PRECISION FUNCTION GK(I,J,K,REL) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      LOGICAL REL
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      CALL YKF(I,J,K,REL) 
      GK = QUADS(I,J,1) 
      IF (MASS .GT. 0) THEN
         IF (MASS .EQ. 1) THEN
            IF (K .EQ. 1) GK = GK + RMASS*GRAD(I,J)**2
         ELSE
            GK = GK*(D1 + RMASS/D2)
            IF (K .EQ. 1) GK = GK + Z*RMASS*QUADR(I,J,1)*QUADR(J,I,-2)
         END IF
      END IF
      RETURN 
      END 
*
*     -----------------------------------------------------------------
*           G R A C A H
*     -----------------------------------------------------------------
*
*
      SUBROUTINE GRACAH(I,J,K,L,M,N,RAC) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
* 
*      SUBROUTINE TO CALCULATE RACAH COEFFICIENTS. 
*      THE ARGUMENTS I,J,K,L,M,N SHOULD BE TWICE THEIR ACTUAL VALUE. 
*      WRITTEN BY N. S. SCOTT 
*      Modified by C. Froese Fischer, March 11, 1988 to use
*         table look-up
*
      LOGICAL SAVE
      COMMON/FACT/GAM(100) 
      DIMENSION RACA(0:4,0:4,0:4,0:4,0:4,0:4)
      DATA RACA/15625*1.D-20/
      DATA ZERO,ONE,TWO,UNDEF/0.D0,1.D0,2.D0,1.D-20/ 
* 
      SAVE = .FALSE.
      IMAX = MAX(I,J,K,L,M,N)
      IF (IMAX .LE. 4) THEN
         RAC = RACA(I,J,K,L,M,N)
         IF (RAC .EQ. UNDEF) THEN
            SAVE = .TRUE.
         ELSE
            RETURN
         END IF
      END IF
      J1 = I+J+M 
      J2 = K+L+M 
      J3 = I+K+N 
      J4 = J+L+N 
      IF (MOD(J1,2) .EQ. 0  .AND.  MOD(J2,2) .EQ. 0   .AND. 
     :    MOD(J3,2) .EQ. 0  .AND.  MOD(J4,2) .EQ. 0 )  THEN 
          J1 = J1/2 
          J2 = J2/2 
          J3 = J3/2 
          J4 = J4/2 
          IF (MAX(I,J,M) .LE. J1 .AND.  MAX(K,L,M) .LE. J2  .AND. 
     :        MAX(I,K,N) .LE. J3 .AND.  MAX(J,L,N) .LE. J4  )  THEN 
              J5 = (I+J+K+L)/2 
              J6 = (I+L+M+N)/2 
              J7 = (J+K+M+N)/2 
              NUMIN = MAX(J1, J2, J3, J4) + 1 
              NUMAX = MIN(J5, J6, J7)     + 1 
              RAC = ONE 
              ICOUNT = 0 
              DO 10 KK = NUMIN+1,NUMAX 
                KI = NUMAX - ICOUNT 
                RAC = ONE - (RAC*(KI*(J5-KI+2)*(J6-KI+2)*(J7-KI+2)))/ 
     :                   ((KI-1-J1)*(KI-1-J2)*(KI-1-J3)*(KI-1-J4)) 
                ICOUNT = ICOUNT+1 
  10          CONTINUE 
              RAC = RAC*EXP(
     :              (GAM(NUMIN+1) - GAM(NUMIN-J1) - GAM(NUMIN-J2) - 
     :               GAM(NUMIN-J3) - GAM(NUMIN-J4) - GAM(J5+2-NUMIN)- 
     :               GAM(J6+2-NUMIN)-GAM(J7+2-NUMIN)) + 
     :              (GAM(J1+1-I)+GAM(J1+1-J)+GAM(J1+1-M)-GAM(J1+2) + 
     :               GAM(J2+1-K)+GAM(J2+1-L)+GAM(J2+1-M)-GAM(J2+2) + 
     :               GAM(J3+1-I)+GAM(J3+1-K)+GAM(J3+1-N)-GAM(J3+2) + 
     :               GAM(J4+1-J)+GAM(J4+1-L)+GAM(J4+1-N)-GAM(J4+2))/TWO) 
              IF (MOD(J5+NUMIN,2) .EQ. 0) RAC = -RAC
          ELSE
              RAC = ZERO
          END IF
      ELSE 
         RAC = ZERO 
      END IF 
      IF (SAVE) RACA(I,J,K,L,M,N) = RAC
      RETURN 
      END 
*     ------------------------------------------------------------------ 
*              H N O R M 
*     ------------------------------------------------------------------ 
* 
*       Returns the value of the normalization constant for an (nl) 
*   hydrogenic function with nuclear charge ZZ. 
* 
* 
      DOUBLE PRECISION FUNCTION HNORM(N,L,ZZ) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID 
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS 
      M = L + L + 1 
      A = N + L 
      B = M 
      T = A 
      D = B 
      M = M - 1 
      IF (M .EQ. 0) GO TO 2 
      DO 1 I = 1,M 
      A = A - D1 
      B = B - D1 
      T = T*A 
1     D = D*B 
2     HNORM = DSQRT(ZZ*T)/( N*D) 
      RETURN 
      END 
* 
*     ------------------------------------------------------------------ 
*              H W F 
*     ------------------------------------------------------------------ 
* 
*       Returns the value of an unnormalized (nl) hydrogenic function 
*   with nuclear charge ZZ and radius r. 
* 
* 
      DOUBLE PRECISION FUNCTION HWF(N,L,ZZ,R) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID 
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS 
      K = N-L-1 
      P = D1 
      A = D1 
      B = K 
      C = N+ L 
      X = -D2*ZZ*R/N 
* 
*  *****  TEST IF UNDERFLOW MAY OCCUR, IF SO SET HWF = 0 
* 
      IF ( X .LT. -150.D0 ) GO TO 5 
      IF (K) 1,2,3 
3     DO 4 I = 1,K 
      P = D1 + A/B*P/C*X 
      A = A + D1 
      B = B - D1 
4     C = C - D1 
2     HWF = P*DEXP(X/D2)*(-X)**(L+1) 
      RETURN 
1     WRITE(6,7) N,L,ZZ,R 
7     FORMAT(51H FORBIDDEN COMBINATION OF N AND L IN HWF SUBPROGRAM/ 
     :    4H N = ,I4,6H   L = ,I4,6H   Z = ,F6.1,6H   R = ,F8.4) 
      STOP 
5     HWF = D0 
      RETURN 
      END 
*
*    ------------------------------------------------------------------
*             I N I T A
*    ------------------------------------------------------------------
*         
*      Initializes basi*constants of the program including those
*  which define the average energy of a configuration.
*         
*         
      SUBROUTINE INITA
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /EAV/CCA(10),CCB(35)
*         
* *****  AVERAGE INTERACTIONS FOR EQUIVALENT ELECTRONS
*         
* *****  P - P
*         
      CCA(1) = 2.D0/25.D0
*         
* *****  D - D
*         
      CCA(2) = 2.D0/63.D0
      CCA(3) = 2.D0/63.D0
*         
* *****  F - F
*         
      CCA(4) =   4.D0/ 195.D0
      CCA(5) =   2.D0/ 143.D0
      CCA(6) = 100.D0/5577.D0
*         
* *****  G - G
*         
      CCA(7) =   20.D0/  1309.D0
      CCA(8) =  162.D0/ 17017.D0
      CCA(9) =   20.D0/  2431.D0
      CCA(10) = 4410.D0/371943.D0
*         
*         
* ***** AVERAGE INTERACTIONS FOR NON-EQUIVALENT ELECTRONS
*         
* *****  S - ( S, P, D, F, G )
*         
      CCB(1) = 1.D0/ 2.D0
      CCB(2) = 1.D0/ 6.D0
      CCB(3) = 1.D0/10.D0
      CCB(4) = 1.D0/14.D0
      CCB(5) = 1.D0/18.D0
*         
* *****  P - ( P, D, F, G )
*         
      CCB(6) = 1.D0/  6.D0
      CCB(7) = 1.D0/ 15.D0
      CCB(8) = 1.D0/ 15.D0
      CCB(9) = 3.D0/ 70.D0
      CCB(10) = 3.D0/ 70.D0
      CCB(11) = 2.D0/ 63.D0
      CCB(12) = 2.D0/ 63.D0
      CCB(13) = 5.D0/198.D0
*         
* *****  D - ( D, F, G )
*         
      CCB(14) =  1.D0/ 10.D0
      CCB(15) =  1.D0/ 35.D0
      CCB(16) =  1.D0/ 35.D0
      CCB(17) =  3.D0/ 70.D0
      CCB(18) =  2.D0/105.D0
      CCB(19) =  5.D0/231.D0
      CCB(20) =  1.D0/ 35.D0
      CCB(21) = 10.D0/693.D0
      CCB(22) =  5.D0/286.D0
*         
* *****  F - ( F, G )
*         
      CCB(23) =  1.D0/  14.D0
      CCB(24) =  2.D0/ 105.D0
      CCB(25) =  1.D0/  77.D0
      CCB(26) = 50.D0/3003.D0
      CCB(27) =  2.D0/  63.D0
      CCB(28) =  1.D0/  77.D0
      CCB(29) = 10.D0/1001.D0
      CCB(30) = 35.D0/2574.D0
*         
* *****  G - ( G )
*         
      CCB(31) =   1.D0/   18.D0
      CCB(32) =  10.D0/  693.D0
      CCB(33) =   9.D0/ 1001.D0
      CCB(34) =  10.D0/ 1287.D0
      CCB(35) = 245.D0/21879.D0
*
* *** Initialize /FACT/
*
      CALL FACTRL(32)
      RETURN
      END  
*
*     ---------------------------------------------------------------
*        I N I T R
*     ---------------------------------------------------------------
*
*
      SUBROUTINE INITR
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*         
* *****  SET THE COMMONLY USED DOUBLE PRECISION CONSTANTS
*         
      D0  =  0.D0
      D1  =  1.D0
      D2  =  2.D0
      D3  =  3.D0
      D4  =  4.D0
      D5  =  1.D0/2.D0
      D6  =  6.D0
      D8  =  8.D0
      D10 = 10.D0
      D12 = 12.D0
      D16 = 16.D0
      D30 = 30.D0
*         
* *****  SET THE STARTING POINT, STEP SIZE, AND RELATED PARAMETERS
*         
      RHO = -4.D0
      H   = 1./16.D0
      H1 = H/1.5
      H3 = H/3.
      CH = H*H/12.
      EH = DEXP(-H)
      NO=NOD
      ND = NO - 2
*
* *****  SET THE FINE-STRUCTURE CONSTANT
*
      FINE = 0.25D0/(137.036D0)**2
      END
*
*     ------------------------------------------------------------------
*       I N T A C T
*     ------------------------------------------------------------------
*
      SUBROUTINE INTACT(L,LP,IEQUIV)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/ENAV/NINTS,KVALUE(15),COEFCT(15)
      COMMON /EAV/CCA(10),CCB(35)
      INTEGER ICBPTR(0:4)
      DATA ICBPTR/1,6,14,23,31/
*
*     THIS SUBROUTINE GIVES THE INTERACTION ENERGY BETWEEN TWO SHELLS,
*     ONE WITH ORBITAL ANGULAR MOMENTUM  L , THE OTHER WITH ORBITAL
*     ANGULAR MOMENTUM  LP .   NOTICE THAT THE FIRST TERM OF THIS
*     INTERACTION ENERGY IS ALWAYS   F0(L,LP)   AND THIS IS NOT GIVEN
*     IN THIS SUBROUTINE.   THUS ONLY THE EXTRA TERMS ARE HERE PRODUCED.
*     FOR EQUIVALENT ELECTRONS (IEQUIV = 1) ,  THERE WILL BE  FK
*     INTEGRALS ONLY.   FOR NON-EQUIVALENT ELECTRONS (IEQUIV = 2) ,
*     THERE WILL BE  GK  INTEGRALS ONLY.
*
*     THE EXPRESSIONS FOR THE INTERACTION ENERGIES ARE GIVEN BY
*     R. D. COWAN, THE THEORY OF ATOMIC SPECTRA, EQUATIONS (6.38)
*      AND (6.39).
*
      I = 0
      IF (IEQUIV .EQ. 1) THEN
          DO 1 K = 2,2*L,2
             I = I+1
             KVALUE(I) = K
             IF (L .LE. 4) THEN
                COEFCT(I) = -CCA((L*(L-1) + K)/2)
             ELSE
               COEFCT(I) = -RME(L,L,K)**2/((2*L+1)*(4*L+1))
             END IF
    1     CONTINUE
      ELSE
          DO 2 K = IABS(L-LP),L+LP,2
             I = I+1
             KVALUE(I) = K
             IF (L .LE. LP) THEN
                L1 = L
                 L2 = LP
              ELSE
                 L1 = LP
                 L2 = L
             END IF
             IF ( L2 .LE. 4) THEN
                COEFCT(I) = -CCB(ICBPTR(L1)+(K+L1-L2)/2+(L1+1)*(L2-L1))
               ELSE
                COEFCT(I) = -RME(L,LP,K)**2/(2*(2*L+1)*(2*LP+1))
             END IF
    2     CONTINUE
      END IF
      NINTS = I
      RETURN
      END
* 
*     ------------------------------------------------------------------ 
*               L I N E Q N 
*     ------------------------------------------------------------------ 
* 
*        This routine is a modification of the one in "Computer Methods 
*     for Mathematical Computation" by Forsythe, Malcolm, and Moler 
*     (Prentice Hall, 1975) 
* 
      SUBROUTINE LINEQN(NDIM,N,A,B) 
* 
      INTEGER NDIM,N 
      DOUBLE PRECISION A(NDIM,N),B(NDIM) 
* 
*     SOLVE A SYSTEM OF LINEAR EQUATIONS A*X = B 
* 
*     INPUT.. 
* 
*        NDIM = DECLARED ROW DIMENSION OF THE ARRAY CONTAINING  A. 
*        N = ORDER OF THE MATRIX 
*        A =  COEFFICIENT MATRIX 
*        B = CONSTANT VECTOR 
* 
*     OUTPUT.. 
* 
*        B = SOLUTION VECTOR 
* 
      DOUBLE PRECISION T 
      INTEGER NM1, I, J, K, KP1, KB, KM1, M 
      NM1 = N-1 
* 
* 
*     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING 
* 
      DO 35 K = 1,NM1 
         KP1= K+1 
* 
*        FIND PIVOT 
* 
         M = K 
         DO 15 I = KP1,N 
            IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I 
   15    CONTINUE 
         T = A(M,K) 
         A(M,K) = A(K,K) 
         A(K,K) = T 
* 
*        SKIP STEP IF PIVOT IS ZERO 
* 
         IF (T .EQ. 0.0D0) GO TO 35 
* 
*        COMPUTE MULTIPLIERS 
* 
         DO 20 I = KP1,N 
             A(I,K) = -A(I,K)/T 
   20    CONTINUE 
* 
*        INTERCHANGE AND ELIMINATE BY COLUMNS 
* 
         DO 30 J = KP1,N 
             T = A(M,J) 
             A(M,J) = A(K,J) 
             A(K,J) = T 
             IF (T .EQ. 0.0D0) GO TO 30 
             DO 25 I = KP1,N 
                A(I,J) = A(I,J) + A(I,K)*T 
   25        CONTINUE 
   30    CONTINUE 
          T = B(M) 
          B(M) = B(K) 
          B(K) = T 
          IF (T .EQ. 0.0D0) GO TO 35 
          DO 32 I = KP1, N 
              B(I) = B(I) + A(I,K)*T 
   32     CONTINUE 
   35 CONTINUE 
* 
*     BACK SUBSTITUTION 
* 
      DO 40 K = N,2,-1 
        IF (A(K,K) .EQ. 0.D0) THEN 
             B(K) = 0.D0 
          ELSE 
             B(K) = B(K)/A(K,K) 
        END IF 
         T = -B(K) 
         DO 41 I = 1, K-1 
             B(I) = B(I) + A(I,K)*T 
   41    CONTINUE 
   40 CONTINUE 
   50 B(1) = B(1)/A(1,1) 
      RETURN 
      END 
*     -----------------------------------------------------------------
*       L V A L
*     -----------------------------------------------------------------
*
*
      INTEGER FUNCTION LVAL(SYMBOL) 
      CHARACTER*1 SYMBOL 
      CHARACTER*22 SET
      DATA SET/'spdfghiklmnSPDFGHIKLMN'/ 
*
      LOCATE = INDEX(SET,SYMBOL) 
      IF ( LOCATE .LE. 11) THEN 
            LVAL = LOCATE - 1 
         ELSE 
            LVAL = LOCATE - 12 
      ENDIF 
      RETURN 
      END 
* 
* 
*     ------------------------------------------------------------------ 
*       P A C K 
*     ------------------------------------------------------------------ 
*  Subroutine written by Bin LIU 
* 
*  Rules for encoding 
*  1. All blanks deleted 
*  2. If Qi=1, omit Qi 
*  3. If Qi=1 or Qi>=4l+1, omit ALFAi 
*  4. If i=1 or (Qi=4l+2 and i<>m), insert '.'; else _BETAi. 
* 
      SUBROUTINE PACK (M, EL, Q, COUPLE, STR) 
 
          INTEGER FULL,Q(5),CONST 
          CHARACTER*3 EL(5),COUPLE(9),CH3
          CHARACTER   CH1
          CHARACTER*66 STR
* 
*  FULL   :  4l+2 
*  CONST  :   constant for converting lowercase to uppercase 
*  CH*    :  temporary variables 
* 
          CONST = ICHAR('a') - ICHAR('A') 
          STR=' ' 
          J = 0
* 
*   -----  begin to encode  ----- 
* 
          DO 100 I=1,M 
            N = Q(I)
            IF (N .EQ. 0) GO TO 100
            K = 3
            IF (EL(I)(3:3) .EQ. ' ') K = 2
            IF (EL(I)(1:1) .EQ. ' ') THEN
                EL(I)=EL(I)(2:3)//' ' 
                K = 2
            END IF
            CH1=EL(I)(2:2) 
            IF ((CH1.GE.'A') .AND. (CH1.LE.'Z')) 
     :            EL(I)(2:2)=CHAR(ICHAR(CH1)+CONST) 
            FULL=4*LVAL(CH1)+2 
* 
*  -----  convert Qi into character  ----- 
* 
            WRITE(CH3,'(I2)') Q(I) 
            STR=STR(1:J)//EL(I)(1:K)
            J= J + K
* 
*  -----  If Qi<>1, add Qi 
*           If Qi<4l+1, add TERMi for the shell ----- 
* 
            IF (N .NE. 1) THEN 
                IF (N .GT. 9) THEN 
                  STR=STR(1:J)//'('//CH3(1:2)//')' 
                  J = J + 4
                ELSE 
                  STR=STR(1:J)//'('//CH3(2:2)//')' 
                  J = J + 3
                ENDIF 
                IF (N .LT. FULL-1 .AND. M .NE. 1) THEN 
                    CH3=COUPLE(I) 
                    CH1= CH3(2:2) 
                    IF (CH1.GE.'a' .AND. CH1.LE.'z') 
     :                  CH3(2:2) = CHAR(ICHAR(CH1)-CONST) 
                    STR=STR(1:J)//CH3 
                    J= J + 3
                ENDIF 
            ENDIF 
* 
*  -----  If i=1 or Qi=4l+2 and i<>m, 
*           insert '.'; else _RESULTANTi.  ----- 
* 
50          IF ((I.NE.1 .AND. N.NE.FULL )
     :           .OR. I.EQ.M  ) THEN 
                CH3=COUPLE(M+I-1) 
                CH1 = CH3(2:2) 
                IF (CH1.GE.'a' .AND. CH1.LE.'z') 
     :                  CH3(2:2) = CHAR(ICHAR(CH1)-CONST) 
                K = 2 
                IF (M .EQ. 1) K = 3 
                STR=STR(1:J)//'_'//CH3(1:K) 
                J = J + K + 1
            ENDIF 
            IF (I .NE. M .AND. N.NE.0 ) THEN 
                J = J + 1
                STR(J:J)='.' 
            ENDIF 
100       CONTINUE 
*
*>>>>>    Because of a compiler error on the SUN, the following is
*         needed to have the string printed correctly.
          STR = STR(1:J)
          RETURN 
        END 
*
*    --------------------------------------------------------------
*            R E F O R M
*    --------------------------------------------------------------
*
*
      SUBROUTINE REFORM(STR1,STR2)
      CHARACTER*40 STR1,STR2,BLANK
      DATA BLANK/'   '/ 
* 
    1 I = 0 
      STR2 = BLANK 
      IS = 0 
    2 JS = INDEX(STR1(IS+1:),'(') 
      IF (JS .NE. 0) THEN 
         IF (JS .GT. 5) GO TO 10 
         I = I+5 
         STR2(I-JS+1:I) = STR1(IS+1:IS+JS) 
         IS = IS + JS 
         JS = INDEX(STR1(IS+1:),')') 
         IF (JS .EQ. 0 .OR. JS .GT. 3) GO TO 10 
         I = I+3 
         STR2(I-JS+1:I) = STR1(IS+1:IS+JS) 
         IS = IS + JS 
         GO TO 2 
      END IF 
      RETURN 
   10 PRINT *,' Error in ',STR1,': Re-enter' 
      READ '(A)',STR1 
      GO TO 1 
      END 
* 
*     ------------------------------------------------------------------ 
*              R K 
*     ------------------------------------------------------------------ 
* 
*                   k 
*       Evaluates  R (i, j; ii, jj) 
* 
* 
      DOUBLE PRECISION FUNCTION RK(I,J,II,JJ,K,REL) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      LOGICAL REL
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      CALL YKF(I,II,K,REL) 
      RK = QUADS(J,JJ,1) 
      IF (MASS .GT. 0) THEN
         IF (MASS .EQ. 1) THEN
            IF (K .EQ. 1) RK = RK - RMASS*GRAD(I,II)*GRAD(J,JJ)
         ELSE
            RK = RK*(D1 + RMASS/D2)
            IF (K .EQ. 1) RK = RK + Z*RMASS/D2*(
     :        QUADR(I,II,1)*QUADR(J,JJ,-2)+QUADR(I,II,-2)*QUADR(J,JJ,1))
         END IF
      END IF
      RETURN 
      END 
*
*     ------------------------------------------------------------------
*         R M E
*     ------------------------------------------------------------------
*
*
      DOUBLE PRECISION FUNCTION RME(L,LP,K)
*
      IMPLICIT REAL *8(A-H,O-Z)
*
      COMMON/FACT/GAM(100)
*
*--- EVALUATES THE REDUCED MATRIX ELEMENT (L//C(K)//LP)  -  SEE FANO
*    AND RACAH, IRREDUCIBLE TENSORIAL SETS, CHAP. 14, P. 81
*
*
      IF (MIN0(L,LP) .EQ. 0) THEN
         RME = 1.D0
       ELSE IF ( K .EQ. 0) THEN
         RME = 2*L+1
         RME = DSQRT(RME)
       ELSE 
         I2G=L+LP+K
         IG=I2G/2
         IF (I2G - 2*IG .NE. 0) THEN
             RME = 0.D0
          ELSE
            I1=IG-L
            I2=IG-LP
            I3=IG-K
            QUSQRT=(2*L+1)*(2*LP+1)
            RME=DSQRT(QUSQRT)*DEXP((GAM(2*I1+1)+GAM(2*I2+1)+GAM(2*I3+1)-
     :        GAM(I2G+2))/2.D0 +GAM(IG+1)-GAM(I1+1)-GAM(I2+1)-GAM(I3+1))
         END IF
      END IF
      RETURN
      END
* 
*     ------------------------------------------------------------------ 
*              S N 
*     ------------------------------------------------------------------ 
* 
*                                      3              k 
*       Evaluates the integral of (1/r)  P (r) P (r) Z (i, j; r)  with 
*                                         i     j 
*   respect to r. 
* 
      DOUBLE PRECISION FUNCTION SN(I,J,II,JJ,K) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      CALL ZK(J,JJ,K) 
      SN = QUADS(I,II,3)*FINE
      RETURN 
      END 
* 
*     ------------------------------------------------------------------ 
*        S Y M M E Q   A N D  S Y M M S L
*     ------------------------------------------------------------------ 
* 
*        This routine is a modification of the one in "Computer Methods 
*     for Mathematical Computation" by Forsythe, Malcolm, and Moler 
*     (Prentice Hall, 1975) to solve a singular system of equations.
* 
      SUBROUTINE SYMMEQ(NDIM,N,A,X)
* 
      DOUBLE PRECISION A(NDIM,N),X(NDIM)
* 
*     SOLVE A SYSTEM OF LINEAR EQUATIONS A*X = 0 
* 
*     INPUT.. 
* 
*        NDIM = DECLARED ROW DIMENSION OF THE ARRAY CONTAINING  A. 
*        N = ORDER OF THE MATRIX 
*        A =  COEFFICIENT MATRIX 
*        B = CONSTANT VECTOR 
* 
*     OUTPUT.. 
* 
*        X = SOLUTION VECTOR with X(1) unchanged.
* 
      DOUBLE PRECISION T 
      INTEGER NM1, I, J, K, KP1, KB, KM1, M 
      NM1 = N-1 
* 
* 
*      L U FACTORIZATION WITHOUT PIVOTING
* 
      DO 35 K = N,3,-1
         KP1= K-1 
         T = A(K,K)
         IF (T .EQ. 0.0D0) GO TO 35 
* 
*        COMPUTE MULTIPLIERS 
* 
         DO 20 I = KP1,2,-1
             A(I,K) = -A(I,K)/T 
   20    CONTINUE 
* 
*        INTERCHANGE AND ELIMINATE BY COLUMNS 
* 
         DO 30 J = KP1,2,-1
             T = A(K,J) 
             IF (T .EQ. 0.0D0) GO TO 30 
             DO 25 I = KP1,2,-1
                A(I,J) = A(I,J) + A(I,K)*T 
   25        CONTINUE 
   30    CONTINUE 
   35 CONTINUE 
*
*     At this point it is assumed that the LU factorization 
*     has already been performed.
*
      ENTRY SYMMSL(NDIM,N,A,X)
*
      DO 36 I = 2,N
        X(I) = - A(I,1)
   36 CONTINUE
      DO 50 K=N,3,-1
        T=X(K)
        IF (T.EQ.0.0D0) GO TO 50
        DO 40 I=K-1,2,-1
          X(I)=X(I)+A(I,K)*T
   40   CONTINUE
   50   CONTINUE
* 
*     BACK SUBSTITUTION 
* 
      DO 80 K =2,N-1
        IF (A(K,K) .EQ. 0.D0) THEN 
             X(K) = 0.D0 
          ELSE 
             X(K) = X(K)/A(K,K) 
        END IF 
         T = -X(K) 
         DO 70 I =K+1,N
             X(I) = X(I) + A(I,K)*T 
   70    CONTINUE 
   80 CONTINUE 
      X(N) = X(N)/A(N,N)
      RETURN 
      END 
*
*     ---------------------------------------------------------------
*        T I N V I T
*     ---------------------------------------------------------------
*
*
      SUBROUTINE TINVIT(NM,N,D,E,E2,M,W,IND,Z, 
     :                  IERR,RV1,RV2,RV3,RV4,RV6) 
* 
      INTEGER I,J,M,N,P,Q,R,S,II,IP,JJ,NM,ITS,TAG,IERR,GROUP 
      DOUBLE PRECISION D(N),E(N),E2(N),W(M),Z(NM,M), 
     :       RV1(N),RV2(N),RV3(N),RV4(N),RV6(N) 
      DOUBLE PRECISION U,V,UK,XU,X0,X1,EPS2,EPS3,EPS4,NORM,ORDER,MACHEP 
      DOUBLE PRECISION DSQRT,DABS,DFLOAT 
      INTEGER IND(M) 
* 
*     THIS SUBROUTINE IS A TRANSLATION OF THE INVERSE ITERATION TECH- 
*     NIQUE IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON. 
*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971). 
* 
*     THIS SUBROUTINE FINDS THOSE EIGENVECTORS OF A TRIDIAGONAL 
*     SYMMETRIC MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES, 
*     USING INVERSE ITERATION. 
* 
*     ON INPUT: 
* 
*        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL 
*          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM 
*          DIMENSION STATEMENT; 
* 
*        N IS THE ORDER OF THE MATRIX; 
* 
*        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX; 
* 
*        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX 
*          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY; 
* 
*        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E, 
*          WITH ZEROS CORRESPONDING TO NEGLIGIBLE ELEMENTS OF E. 
*          E(I) IS CONSIDERED NEGLIGIBLE IF IT IS NOT LARGER THAN 
*          THE PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE SUM 
*          OF THE MAGNITUDES OF D(I) AND D(I-1).  E2(1) MUST CONTAIN 
*          0.0D0 IF THE EIGENVALUES ARE IN ASCENDING ORDER, OR 2.0D0 
*          IF THE EIGENVALUES ARE IN DESCENDING ORDER.  IF  BISECT, 
*          TRIDIB, OR  IMTQLV  HAS BEEN USED TO FIND THE EIGENVALUES, 
*          THEIR OUTPUT E2 ARRAY IS EXACTLY WHAT IS EXPECTED HERE; 
* 
*        M IS THE NUMBER OF SPECIFIED EIGENVALUES; 
* 
*        W CONTAINS THE M EIGENVALUES IN ASCENDING OR DESCENDING ORDER; 
* 
*        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES 
*          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W -- 
*          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM 
*          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC. 
* 
*     ON OUTPUT: 
* 
*        ALL INPUT ARRAYS ARE UNALTERED; 
* 
*        Z CONTAINS THE ASSOCIATED SET OF ORTHONORMAL EIGENVECTORS. 
*          ANY VECTOR WHICH FAILS TO CONVERGE IS SET TO ZERO; 
* 
*        IERR IS SET TO 
*          ZERO       FOR NORMAL RETURN, 
*          -R         IF THE EIGENVECTOR CORRESPONDING TO THE R-TH 
*                     EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS; 
* 
*        RV1, RV2, RV3, RV4, AND RV6 ARE TEMPORARY STORAGE ARRAYS. 
* 
*     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW, 
*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY 
* 
*     ------------------------------------------------------------------ 
* 
*     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING 
*                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC. 
*                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC 
*                ON S360 :::::::::: 
      DATA MACHEP/1.D-12/ 
* 
      IERR = 0 
      IF (M .EQ. 0) GO TO 1001 
      TAG = 0 
      ORDER = 1.0D0 - E2(1) 
      Q = 0 
*     :::::::::: ESTABLISH AND PROCESS NEXT SUBMATRIX :::::::::: 
  100 P = Q + 1 
* 
      DO 120 Q = P, N 
         IF (Q .EQ. N) GO TO 140 
         IF (E2(Q+1) .EQ. 0.0D0) GO TO 140 
  120 CONTINUE 
*     :::::::::: FIND VECTORS BY INVERSE ITERATION :::::::::: 
  140 TAG = TAG + 1 
      S = 0 
* 
      DO 920 R = 1, M 
         IF (IND(R) .NE. TAG) GO TO 920 
         ITS = 1 
         X1 = W(R) 
         IF (S .NE. 0) GO TO 510 
*     :::::::::: CHECK FOR ISOLATED ROOT :::::::::: 
         XU = 1.0D0 
         IF (P .NE. Q) GO TO 490 
         RV6(P) = 1.0D0 
         GO TO 870 
  490    NORM = DABS(D(P)) 
         IP = P + 1 
* 
         DO 500 I = IP, Q 
  500    NORM = NORM + DABS(D(I)) + DABS(E(I)) 
*     :::::::::: EPS2 IS THE CRITERION FOR GROUPING, 
*                EPS3 REPLACES ZERO PIVOTS AND EQUAL 
*                ROOTS ARE MODIFIED BY EPS3, 
*                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW :::::::::: 
         EPS2 = 1.0D-3 * NORM 
         EPS3 = MACHEP * NORM 
         UK = DFLOAT(Q-P+1) 
         EPS4 = UK * EPS3 
         UK = EPS4 / DSQRT(UK) 
         S = P 
  505    GROUP = 0 
         GO TO 520 
*     :::::::::: LOOK FOR CLOSE OR COINCIDENT ROOTS :::::::::: 
  510    IF (DABS(X1-X0) .GE. EPS2) GO TO 505 
         GROUP = GROUP + 1 
         IF (ORDER * (X1 - X0) .LE. 0.0D0) X1 = X0 + ORDER * EPS3 
*     :::::::::: ELIMINATION WITH INTERCHANGES AND 
*                INITIALIZATION OF VECTOR :::::::::: 
  520    V = 0.0D0 
* 
         DO 580 I = P, Q 
            RV6(I) = UK 
            IF (I .EQ. P) GO TO 560 
            IF (DABS(E(I)) .LT. DABS(U)) GO TO 540 
*     :::::::::: WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF 
*                E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY :::::::::: 
            XU = U / E(I) 
            RV4(I) = XU 
            RV1(I-1) = E(I) 
            RV2(I-1) = D(I) - X1 
            RV3(I-1) = 0.0D0 
            IF (I .NE. Q) RV3(I-1) = E(I+1) 
            U = V - XU * RV2(I-1) 
            V = -XU * RV3(I-1) 
            GO TO 580 
  540       XU = E(I) / U 
            RV4(I) = XU 
            RV1(I-1) = U 
            RV2(I-1) = V 
            RV3(I-1) = 0.0D0 
  560       U = D(I) - X1 - XU * V 
            IF (I .NE. Q) V = E(I+1) 
  580    CONTINUE 
* 
         IF (U .EQ. 0.0D0) U = EPS3 
         RV1(Q) = U 
         RV2(Q) = 0.0D0 
         RV3(Q) = 0.0D0 
*     :::::::::: BACK SUBSTITUTION 
*                FOR I=Q STEP -1 UNTIL P DO -- :::::::::: 
  600    DO 620 II = P, Q 
            I = P + Q - II 
            RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I) 
            V = U 
            U = RV6(I) 
  620    CONTINUE 
*     :::::::::: ORTHOGONALIZE WITH RESPECT TO PREVIOUS 
*                MEMBERS OF GROUP :::::::::: 
         IF (GROUP .EQ. 0) GO TO 700 
         J = R 
* 
         DO 680 JJ = 1, GROUP 
  630       J = J - 1 
            IF (IND(J) .NE. TAG) GO TO 630 
            XU = 0.0D0 
* 
            DO 640 I = P, Q 
  640       XU = XU + RV6(I) * Z(I,J) 
* 
            DO 660 I = P, Q 
  660       RV6(I) = RV6(I) - XU * Z(I,J) 
* 
  680    CONTINUE 
* 
  700    NORM = 0.0D0 
* 
         DO 720 I = P, Q 
  720    NORM = NORM + DABS(RV6(I)) 
* 
         IF (NORM .GE. 1.0D0) GO TO 840 
*     :::::::::: FORWARD SUBSTITUTION :::::::::: 
         IF (ITS .EQ. 5) GO TO 830 
         IF (NORM .NE. 0.0D0) GO TO 740 
         RV6(S) = EPS4 
         S = S + 1 
         IF (S .GT. Q) S = P 
         GO TO 780 
  740    XU = EPS4 / NORM 
* 
         DO 760 I = P, Q 
  760    RV6(I) = RV6(I) * XU 
*     :::::::::: ELIMINATION OPERATIONS ON NEXT VECTOR 
*                ITERATE :::::::::: 
  780    DO 820 I = IP, Q 
            U = RV6(I) 
*     :::::::::: IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE 
*                WAS PERFORMED EARLIER IN THE 
*                TRIANGULARIZATION PROCESS :::::::::: 
            IF (RV1(I-1) .NE. E(I)) GO TO 800 
            U = RV6(I-1) 
            RV6(I-1) = RV6(I) 
  800       RV6(I) = U - RV4(I) * RV6(I-1) 
  820    CONTINUE 
* 
         ITS = ITS + 1 
         GO TO 600 
*     :::::::::: SET ERROR -- NON-CONVERGED EIGENVECTOR :::::::::: 
  830    IERR = -R 
         XU = 0.0D0 
         GO TO 870 
*     :::::::::: NORMALIZE SO THAT SUM OF SQUARES IS 
*                1 AND EXPAND TO FULL ORDER :::::::::: 
  840    U = 0.0D0 
* 
         DO 860 I = P, Q 
  860    U = U + RV6(I)**2 
* 
         XU = 1.0D0 / DSQRT(U) 
* 
  870    DO 880 I = 1, N 
  880    Z(I,R) = 0.0D0 
* 
         DO 900 I = P, Q 
  900    Z(I,R) = RV6(I) * XU 
* 
         X0 = X1 
  920 CONTINUE 
* 
      IF (Q .LT. N) GO TO 100 
 1001 RETURN 
*     :::::::::: LAST CARD OF TINVIT :::::::::: 
      END 
*
*     ----------------------------------------------------------------
*        T Q L 2
*     ----------------------------------------------------------------
*
*
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      DOUBLE PRECISION MACHEP 
      DIMENSION D(N),E(N),Z(NM,N) 
      MACHEP=1.387878878078144568D-17 
      IERR = 0 
      IF (N .EQ. 1) GO TO 1001 
      DO 100 I = 2, N 
  100 E(I-1) = E(I) 
      F = 0.0 
      B = 0.0 
      E(N) = 0.0 
      DO 240 L = 1, N 
         J = 0 
         H =MACHEP*(ABS(D(L)) + ABS(E(L))) 
         IF (B .LT. H) B = H 
         DO 110 M = L, N 
            IF (ABS(E(M)) .LE. B) GO TO 120 
  110    CONTINUE 
  120    IF (M .EQ. L) GO TO 220 
  130    IF (J .EQ. 30) GO TO 1000 
         J = J + 1 
         P = (D(L+1) - D(L)) / (2.0 * E(L)) 
         R = SQRT(P*P+1.0) 
         H = D(L) - E(L) / (P + SIGN(R,P)) 
         DO 140 I = L, N 
  140    D(I) = D(I) - H 
         F = F + H 
         P = D(M) 
         C = 1.0 
         S = 0.0 
         MML = M - L 
         DO 200 II = 1, MML 
            I = M - II 
            G = C * E(I) 
            H = C * P 
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150 
            C = E(I) / P 
            R = SQRT(C*C+1.0) 
            E(I+1) = S * P * R 
            S = C / R 
            C = 1.0 / R 
            GO TO 160 
  150       C = P / E(I) 
            R = SQRT(C*C+1.0) 
            E(I+1) = S * E(I) * R 
            S = 1.0 / R 
            C = C * S 
  160       P = C * D(I) - S * G 
            D(I+1) = H + S * (C * G + S * D(I)) 
            DO 180 K = 1, N 
               H = Z(K,I+1) 
               Z(K,I+1) = S * Z(K,I) + C * H 
               Z(K,I) = C * Z(K,I) - S * H 
  180       CONTINUE 
  200    CONTINUE 
         E(L) = S * P 
         D(L) = C * P 
         IF (ABS(E(L)) .GT. B) GO TO 130 
  220    D(L) = D(L) + F 
  240 CONTINUE 
      GO TO 1001 
 1000 IERR = L 
 1001 RETURN 
      END 
*
*     --------------------------------------------------------------
*        T R B A K 1
*     --------------------------------------------------------------
*
*
      SUBROUTINE TRBAK1(NM,N,A,E,M,Z) 
* 
      INTEGER I,J,K,L,M,N,NM 
      DOUBLE PRECISION A(NM,N),E(N),Z(NM,M) 
      DOUBLE PRECISION S 
* 
*     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK1, 
*     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON. 
*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971). 
* 
*     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC 
*     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING 
*     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRED1. 
* 
*     ON INPUT: 
* 
*        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL 
*          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM 
*          DIMENSION STATEMENT; 
* 
*        N IS THE ORDER OF THE MATRIX; 
* 
*        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS- 
*          FORMATIONS USED IN THE REDUCTION BY  TRED1 
*          IN ITS STRICT LOWER TRIANGLE; 
* 
*        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL 
*          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY; 
* 
*        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED; 
* 
*        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED 
*          IN ITS FIRST M COLUMNS. 
* 
*     ON OUTPUT: 
* 
*        Z CONTAINS THE TRANSFORMED EIGENVECTORS 
*          IN ITS FIRST M COLUMNS. 
* 
*     NOTE THAT TRBAK1 PRESERVES VECTOR EUCLIDEAN NORMS. 
* 
*     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW, 
*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY 
* 
*     ------------------------------------------------------------------ 
* 
      IF (M .EQ. 0) GO TO 200 
      IF (N .EQ. 1) GO TO 200 
* 
      DO 140 I = 2, N 
         L = I - 1 
         IF (E(I) .EQ. 0.0D0) GO TO 140 
* 
         DO 130 J = 1, M 
            S = 0.0D0 
* 
            DO 110 K = 1, L 
  110       S = S + A(I,K) * Z(K,J) 
*     :::::::::: DIVISOR BELOW IS NEGATIVE OF H FORMED IN TRED1. 
*                DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW :::::::::: 
            S = (S / A(I,L)) / E(I) 
* 
            DO 120 K = 1, L 
  120       Z(K,J) = Z(K,J) + S * A(I,K) 
* 
  130    CONTINUE 
* 
  140 CONTINUE 
* 
  200 RETURN 
*     :::::::::: LAST CARD OF TRBAK1 :::::::::: 
      END 
*
*     --------------------------------------------------------------
*       T R E D 1
*     --------------------------------------------------------------
*
*
      SUBROUTINE TRED1(NM,N,A,D,E,E2) 
* 
      INTEGER I,J,K,L,N,II,NM,JP1 
      DOUBLE PRECISION A(NM,N),D(N),E(N),E2(N) 
      DOUBLE PRECISION F,G,H,SCALE 
      DOUBLE PRECISION DSQRT,DABS,DSIGN 
* 
*     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1, 
*     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON. 
*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971). 
* 
*     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX 
*     TO A SYMMETRIC TRIDIAGONAL MATRIX USING 
*     ORTHOGONAL SIMILARITY TRANSFORMATIONS. 
* 
*     ON INPUT: 
* 
*        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL 
*          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM 
*          DIMENSION STATEMENT; 
* 
*        N IS THE ORDER OF THE MATRIX; 
* 
*        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE 
*          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED. 
* 
*     ON OUTPUT: 
* 
*        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS- 
*          FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER 
*          TRIANGLE.  THE FULL UPPER TRIANGLE OF A IS UNALTERED; 
* 
*        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX; 
* 
*        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL 
*          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO; 
* 
*        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E. 
*          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED. 
* 
*     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW, 
*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY 
* 
*     ------------------------------------------------------------------ 
* 
      DO 100 I = 1, N 
  100 D(I) = A(I,I) 
*     :::::::::: FOR I=N STEP -1 UNTIL 1 DO -- :::::::::: 
      DO 300 II = 1, N 
         I = N + 1 - II 
         L = I - 1 
         H = 0.0D0 
         SCALE = 0.0D0 
         IF (L .LT. 1) GO TO 130 
*     :::::::::: SCALE ROW (ALGOL TOL THEN NOT NEEDED) :::::::::: 
         DO 120 K = 1, L 
  120    SCALE = SCALE + DABS(A(I,K)) 
* 
         IF (SCALE .NE. 0.0D0) GO TO 140 
  130    E(I) = 0.0D0 
         E2(I) = 0.0D0 
         GO TO 290 
* 
  140    DO 150 K = 1, L 
            A(I,K) = A(I,K) / SCALE 
            H = H + A(I,K) * A(I,K) 
  150    CONTINUE 
* 
         E2(I) = SCALE * SCALE * H 
         F = A(I,L) 
         G = -DSIGN(DSQRT(H),F) 
         E(I) = SCALE * G 
         H = H - F * G 
         A(I,L) = F - G 
         IF (L .EQ. 1) GO TO 270 
         F = 0.0D0 
* 
         DO 240 J = 1, L 
            G = 0.0D0 
*     :::::::::: FORM ELEMENT OF A*U :::::::::: 
            DO 180 K = 1, J 
  180       G = G + A(J,K) * A(I,K) 
* 
            JP1 = J + 1 
            IF (L .LT. JP1) GO TO 220 
* 
            DO 200 K = JP1, L 
  200       G = G + A(K,J) * A(I,K) 
*     :::::::::: FORM ELEMENT OF P :::::::::: 
  220       E(J) = G / H 
            F = F + E(J) * A(I,J) 
  240    CONTINUE 
* 
         H = F / (H + H) 
*     :::::::::: FORM REDUCED A :::::::::: 
         DO 260 J = 1, L 
            F = A(I,J) 
            G = E(J) - H * F 
            E(J) = G 
* 
            DO 260 K = 1, J 
               A(J,K) = A(J,K) - F * E(K) - G * A(I,K) 
  260    CONTINUE 
* 
  270    DO 280 K = 1, L 
  280    A(I,K) = SCALE * A(I,K) 
* 
  290    H = D(I) 
         D(I) = A(I,I) 
         A(I,I) = H 
  300 CONTINUE 
* 
      RETURN 
*     :::::::::: LAST CARD OF TRED1 :::::::::: 
      END 
*
*     ----------------------------------------------------------------
*        T R E D 2
*     ----------------------------------------------------------------
*
*
      SUBROUTINE TRED2(NM,N,A,D,E,Z) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      DIMENSION A(NM,N),D(N),E(N),Z(NM,N) 
      DO 100 I = 1, N 
         DO 100 J = 1, I 
            Z(I,J) = A(I,J) 
  100 CONTINUE 
      IF (N .EQ. 1) GO TO 320 
      DO 300 II = 2, N 
         I = N + 2 - II 
         L = I - 1 
         H = 0.0 
         SCALE = 0.0 
         IF (L .LT. 2) GO TO 130 
         DO 120 K = 1, L 
  120    SCALE = SCALE + ABS(Z(I,K)) 
         IF (SCALE .NE. 0.0) GO TO 140 
  130    E(I) = Z(I,L) 
         GO TO 290 
  140    DO 150 K = 1, L 
            Z(I,K) = Z(I,K) / SCALE 
            H = H + Z(I,K) * Z(I,K) 
  150    CONTINUE 
         F = Z(I,L) 
         G = -SIGN(SQRT(H),F) 
         E(I) = SCALE * G 
         H = H - F * G 
         Z(I,L) = F - G 
         F = 0.0 
         DO 240 J = 1, L 
            Z(J,I) = Z(I,J) / (SCALE * H) 
            G = 0.0 
            DO 180 K = 1, J 
  180       G = G + Z(J,K) * Z(I,K) 
            JP1 = J + 1 
            IF (L .LT. JP1) GO TO 220 
            DO 200 K = JP1, L 
  200       G = G + Z(K,J) * Z(I,K) 
  220       E(J) = G / H 
            F = F + E(J) * Z(I,J) 
  240    CONTINUE 
         HH = F / (H + H) 
         DO 260 J = 1, L 
            F = Z(I,J) 
            G = E(J) - HH * F 
            E(J) = G 
            DO 260 K = 1, J 
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K) 
  260    CONTINUE 
         DO 280 K = 1, L 
  280    Z(I,K) = SCALE * Z(I,K) 
  290    D(I) = H 
  300 CONTINUE 
  320 D(1) = 0.0 
      E(1) = 0.0 
      DO 500 I = 1, N 
         L = I - 1 
         IF (D(I) .EQ. 0.0) GO TO 380 
         DO 360 J = 1, L 
            G = 0.0 
            DO 340 K = 1, L 
  340       G = G + Z(I,K) * Z(K,J) 
            DO 360 K = 1, L 
               Z(K,J) = Z(K,J) - G * Z(K,I) 
  360    CONTINUE 
  380    D(I) = Z(I,I) 
         Z(I,I) = 1.0 
         IF (L .LT. 1) GO TO 500 
         DO 400 J = 1, L 
            Z(I,J) = 0.0 
            Z(J,I) = 0.0 
  400    CONTINUE 
  500 CONTINUE 
      RETURN 
      END 
* 
*     ------------------------------------------------------------------ 
*              V K
*     ------------------------------------------------------------------ 
* 
*                  k 
*       Evaluates V (i,j) as defined by Blume and Watson (1962). 
* 
      DOUBLE PRECISION FUNCTION VK(I,J,II,JJ,K) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      CALL DYK(I,II,K) 
      VK = QUADS(J,JJ,2)*FINE
      RETURN 
      END 
