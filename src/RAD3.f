*
*     Routines for MCHF_LIB_RAD3
*     Computer Physics Communication, Vol. 64, 399-405 (1991)
*
*     C O P Y R I G H T -- 1994
*
*     ------------------------------------------------------------------ 
*              D Y K 
*     ------------------------------------------------------------------ 
* 
*       Stores in YK the values of the integral of 
*              k 
*       P (s/r) (dP /ds - P /s) integrated over the interval (0,r) 
*        i         j       j 
* 
*   which enter into the spin-orbit calculation. 
* 
* 
      SUBROUTINE DYK(I,J,K) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      PARAMETER(NOD=220,NWD=30) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID 
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS 
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD), 
     :   YR(NOD),X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD) 
      DEN = L(I)+ L(J)+ 2 + K 
      FACT = L(J) 
      DO 1 JJ =1,2 
1     YK(JJ) = FACT*P(JJ,I)*P(JJ,J)*R(JJ)/DEN 
      A = EH**K 
      AA = A*A 
      A = D4*A 
      F1 = FACT*P(1,I)*P(1,J)*R(1) 
      F2 = FACT*P(2,I)*P(2,J)*R(2) 
      DO 8 M =3,ND 
      F3 = (-P(M+2,J) + D8*(P(M+1,J)-P(M-1,J)) + P(M-2,J))/(D6*H) 
      F3 = D5*P(M,I)*(F3 - P(M,J))*R(M) 
      YK(M) = YK(M-2)*AA + H3*(F3+ A*F2 + AA*F1) 
      F1 = F2 
8     F2 = F3 
      A = A*(EH)**3 
      AA = A*A/D16 
      C = 2*K+3 
      HH = C*H3 
      YK(NO)= YK(ND) 
      F1 = YK(NO) 
      F2 = F1 
      DO 9 MM = 3,NO 
      M = NO -MM+1 
      F3 =YK(M) 
      YK(M) = YK(M+2)*AA + HH*(F3 +A*F2 + AA*F1) 
      F1 = F2 
9     F2 = F3 
      RETURN 
      END 
*
*    ------------------------------------------------------------------
*      E C O R E
*    ------------------------------------------------------------------
*
*
      SUBROUTINE ECORE(EL,EC,REL) 
* 
* *** COMPUTES THE ENERGY OF THE COMMON CLOSED SHELLS 
* 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      PARAMETER (NOD=220,NWD=30)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID 
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS 
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD), 
     :   YR(NOD),X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD) 
      CHARACTER EL(*)*3
      LOGICAL REL
* 
      EC = D0 
      DO 10 I = 1,NCLOSD 
         SUMI = 4*L(I)+2 
         TI   = FK(I,I,0,REL) 
         DO 20 K = 2,2*L(I),2 
            TI = TI - CA(L(I),K)*FK(I,I,K,REL) 
   20    CONTINUE 
         EC = EC + SUMI*((SUMI-1)*TI - HL(EL,I,I,REL))/D2 
         DO 30 J = 1,I-1 
            SUMJ = 4*L(J)+2 
            TIJ = FK(I,J,0,REL) 
            DO 40 K=IABS(L(I)-L(J)),L(I)+L(J),2 
               TIJ = TIJ -CB(L(I),L(J),K)*GK(I,J,K,REL) 
   40       CONTINUE 
            EC = EC + SUMI*SUMJ*TIJ 
   30    CONTINUE 
   10 CONTINUE 
      END 
* 
*     ------------------------------------------------------------------ 
*              G R A D 
*     ------------------------------------------------------------------ 
* 
*  *****  THE GRAD FUNCTION SUBPROGRAM COMPUTES THE FOLLOWING DIRECTLY 
*  *****         <P(J)^D + L(I)/R ^P(I)> WITH L(I) > L(J) 
* 
      DOUBLE PRECISION FUNCTION GRAD(I,J) 
      PARAMETER (NOD=220,NWD=30) 
      PARAMETER (IWRITE=6)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
      IF ( IABS(L(I) - L(J)) .NE. 1) GO TO 100 
      LL = MAX0(L(I),L(J)) 
      II = I 
      JJ = J 
      IF ( L(I) .GT. L(J) ) GO TO 1 
      II = J 
      JJ = I 
1     A1 = (LL+D5)/(LL*(LL+1)*(2*LL+1)) 
      GRAD = R(1)*P(1,I)*P(1,J)*(D1 + A1*Z*R(1)) 
      DL = D5*P(1,I)*P(1,J)*R(1) 
      MM = MIN0(MAX(I)+1,MAX(J)+1,ND) 
      K = 2 
      F1 = D5*(P(K+1,II) - P(K-1,II)) 
      F2 = P(K+1,II) - D2*P(K,II) + P(K-1,II) 
      G0 = P(K,JJ)*R(K) 
      G1 = D5*(P(K+1,JJ)*R(K+1) - P(K-1,JJ)*R(K-1)) 
      G2 = P(K+1,JJ)*R(K+1) - D2*P(K,JJ)*R(K) + P(K-1,JJ)*R(K-1) 
      GRAD = GRAD + D2*F1*G0 +(D2*F2*G1 + F1*G2)/D3 
      DL = DL + D2*P(K,II)*P(K,JJ)*R(K) + P(K+1,II)*P(K+1,JJ)*R(K+1) 
      DO 2 K = 4,MM,2 
      F1 = D5*(P(K+1,II) - P(K-1,II)) 
      F2 = P(K+1,II) - D2*P(K,II) + P(K-1,II) 
      F3 = D5*(P(K+2,II) - P(K-2,II)) - D2*F1 
      F4 = P(K+2,II) + P(K-2,II) - D4*(P(K+1,II) + P(K-1,II)) 
     :   + D6*P(K,II) 
      G0 = P(K,JJ)*R(K) 
      G1 = D5*(P(K+1,JJ)*R(K+1) - P(K-1,JJ)*R(K-1)) 
      G2 = P(K+1,JJ)*R(K+1) - D2*P(K,JJ)*R(K) + P(K-1,JJ)*R(K-1) 
      G3 = D5*(P(K+2,JJ)*R(K+2) - P(K-2,JJ)*R(K-2)) -D2*G1 
      G4 = P(K+2,JJ)*R(K+2) + P(K-2,JJ)*R(K-2) - D4*(P(K+1,JJ)*R(K+1) 
     :   + P(K-1,JJ)*R(K-1)) + D6*P(K,JJ)*R(K) 
      GRAD = GRAD + D2*F1*G0 +(D2*F2*G1 + F1*G2)/D3 
     :   - (F1*G4-F4*G1 + D4*(F2*G3-F3*G2))/90.D0 
 2    DL = DL + D2*P(K,II)*P(K,JJ)*R(K) + P(K+1,II)*P(K+1,JJ)*R(K+1) 
      GRAD = GRAD + (LL+D5)*DL*H1 
      IF (II .EQ. I) GRAD = - GRAD 
      RETURN 
100   WRITE(IWRITE,101) I,J 
101   FORMAT(5X,'L(I)-L(J) NOT =1 FOR I = ',I2,' AND J = ',I2) 
      STOP 
      END 
* 
*     ------------------------------------------------------------------ 
*              H L 
*     ------------------------------------------------------------------ 
* 
*       Returns the value of <i|L|j>, using a special formula to 
*  preserve symmetry. 
* 
      DOUBLE PRECISION FUNCTION HL(EL,I,J,REL) 
      PARAMETER (NOD=220,NWD=30) 
      PARAMETER(IWRITE=6)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      CHARACTER EL(*)*3 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
      LOGICAL REL
      IF (IABS(L(I)-L(J)) .EQ. 0) GO TO 3 
      WRITE(IWRITE,4) EL(I),L(I),EL(J),L(J) 
4     FORMAT(10X,'UNALLOWED L VALUES OCCURRED IN HL SUBROUTINE'/ 
     :   2(10X,A3,' HAS L = ',I3)) 
3     LI = L(I) 
      C = 2*LI + 1 
      A1 = -D2/(C*(LI+1)) 
      A2 = A1/((C+D2)*(LI+1)) 
      A3 = A2/((LI+2)*(LI+1)) 
      ZR = Z*R(1) 
      HL = H*C*P(1,I)*P(1,J)*(D1+ZR*(A1+ZR*(A2+ZR*A3))) 
      MM = MIN0(MAX(I)+3,MAX(J)+3,ND-1) 
      K = 2 
      C = D4/D3 
      DI1 = P(K+1,I) - P(K-1,I) 
      DI2 = P(K+1,I) - D2*P(K,I) + P(K-1,I) 
      DJ1 = P(K+1,J) - P(K-1,J) 
      DJ2 = P(K+1,J) - D2*P(K,J) + P(K-1,J) 
      HL = HL + DI1*DJ1 + C*DI2*DJ2 
      DO 1 K = 4,MM,2 
      DI1 = P(K+1,I) - P(K-1,I) 
      DI2 = P(K+1,I) - D2*P(K,I) + P(K-1,I) 
      DI4 = P(K+2,I) - D4*(P(K+1,I)+P(K-1,I)) + D6*P(K,I) +P(K-2,I) 
      DI3 = P(K+2,I) - P(K-2,I) - D2*DI1 
      DI5 = P(K+3,I)-P(K-3,I) - D4*(P(K+2,I)-P(K-2,I)) 
     :   + 5.D0*(P(K+1,I)-P(K-1,I)) 
      DI6 = P(K+3,I)+P(K-3,I) - D6*(P(K+2,I)+P(K-2,I)) 
     :   + 15.D0*(P(K+1,I)+P(K-1,I)) - 20.D0*P(K,I) 
      DJ1 = P(K+1,J) - P(K-1,J) 
      DJ2 = P(K+1,J) - D2*P(K,J) + P(K-1,J) 
      DJ4 = P(K+2,J) - D4*(P(K+1,J)+P(K-1,J)) + D6*P(K,J) +P(K-2,J) 
      DJ3 = P(K+2,J) - P(K-2,J) - D2*DJ1 
      DJ5 = P(K+3,J)-P(K-3,J) - D4*(P(K+2,J)-P(K-2,J)) 
     :   + 5.D0*(P(K+1,J)-P(K-1,J)) 
      DJ6 = P(K+3,J)+P(K-3,J) - D6*(P(K+2,J)+P(K-2,J)) 
     :   + 15.D0*(P(K+1,J)+P(K-1,J)) - 20.D0*P(K,J) 
1     HL = HL + DI1*DJ1 + C*DI2*DJ2 + (DI3*DJ3 + DI2*DJ4+DI4*DJ2)/45.D0 
     :  -(DI3*DJ5+DI5*DJ3)/252.D0 - (DI2*DJ6+DI6*DJ2-1.1*DI4*DJ4)/378.D0 
      TZ = Z + Z 
      C = (LI + D5)**2 
      HL2 = D5*(TZ*R(1) - C)*P(1,I)*P(1,J) 
      DO 2 K = 2,MM,2 
2     HL2 = HL2 + D2*(TZ*R(K) - C)*P(K,I)*P(K,J) 
     :   + (TZ*R(K+1) - C)*P(K+1,I)*P(K+1,J) 
      HL = -HL/(D2*H) + HL2*H1 
      IF (REL) HL=HL-D2*RLSHFT(I,J)
      RETURN 
      END 
*
*     -----------------------------------------------------------------
*       H L C
*     -----------------------------------------------------------------
*
*
      DOUBLE PRECISION FUNCTION HLC(EL,I,J,REL) 
* 
* *** COMPUTES HL(I,J) MODIFIED BY THE INTERACTIONS WITH THE CLOSED 
*     SHELL 
* 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      PARAMETER(NOD=220,NWD=30)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID 
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS 
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD), 
     :   YR(NOD),X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD) 
* 
      CHARACTER EL(*)*3
      LOGICAL REL
      HLC = HL(EL,I,J,REL) 
      DO 10 IP = 1,NCLOSD 
         SUMIP = 4*L(IP)+2 
         T = RK(I,IP,J,IP,0,REL) 
         DO 20 K = IABS(L(I)-L(IP)),L(I)+L(IP),2 
            T = T - CB(L(I),L(IP),K)*RK(I,IP,IP,J,K,REL) 
   20    CONTINUE 
         HLC = HLC - D2*SUMIP*T 
   10 CONTINUE 
      END 
* 
*     ------------------------------------------------------------------ 
*              Q U A D R 
*     ------------------------------------------------------------------ 
* 
*                                   kk 
*       Evaluates the integral of  r   P (r) P (r) with respect to r 
*                                       i     j 
* 
      DOUBLE PRECISION FUNCTION QUADR(I,J,KK) 
      PARAMETER (NOD=220,NWD=30) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
      K = KK + 2 
      LI = L(I) 
      LJ = L(J) 
      DEN = LI + LJ + 1 + K 
      ZR = Z*R(4) 
      BI = (P(4,I)/(AZ(I)*R2(4)*R(4)**LI) - D1+ZR/(LI+1) )/ZR**2 
      BJ = (P(4,J)/(AZ(J)*R2(4)*R(4)**LJ) - D1+ZR/(LJ+1) )/ZR**2 
      ALPHA= (D1/(LI + 1) + D1/(LJ + 1))/(DEN + D1) 
      ZR = Z*R(1) 
      BETA = (DEN+D1)*ALPHA**2 - D2*(BI+BJ+D1/((LI+1)*(LJ+1)))/(DEN+D2) 
      D = P(1,I)*P(1,J)*R(1)**K*(((BETA*ZR+ALPHA)*ZR+D1)/(DEN*H1)+D5) 
      DD = D0
      M = MIN0(MAX(I),MAX(J)) - 1 
      DO 1 JJ = 2,M,2 
      JP = JJ + 1 
      D = D +P(JP,I)*P(JP,J)*R(JP)**K 
      DD = DD + P(JJ,I)*P(JJ,J)*R(JJ)**K
1     CONTINUE
      QUADR = H1*(D + D2*DD)
      RETURN 
      END 
* 
*     ------------------------------------------------------------------ 
*              Q U A D S 
*     ------------------------------------------------------------------ 
* 
* 
*                                       kk 
*       Evaluates the integral of  (1/r)   YK(r) P (r) P (r)  with 
*                                                 i     j 
*   respect to r. 
* 
* 
      DOUBLE PRECISION FUNCTION  QUADS(I,J,KK) 
      PARAMETER (NOD=220,NWD=30) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
      DEN = L(I) + L(J) + 3 
      K = 2 - KK 
      CD = D1 + Z*R(1)*(DEN-D1)/((DEN+D1)*((L(I)+1)*(L(J)+1))) 
      D = YK(1)*P(1,I)*P(1,J)*R(1)**K*( CD/(DEN*H1)+ D5) 
      DD = D0
      MX = MIN0(MAX(I),MAX(J)) - 1 
      DO 1 M = 2,MX,2 
      DD = DD + YK(M)*P(M,I)*P(M,J)*R(M)**K  
      D= D+  YK(M+1)*P(M+1,I)*P(M+1,J)*R(M+1)**K 
1     CONTINUE
      QUADS = H1*(D + D2*DD)
      RETURN 
      END 
* 
*     ------------------------------------------------------------------ 
*              S H I F T 
*     ------------------------------------------------------------------ 
* 
* 
*       Computes the mass velocity  and one-body   Darwin term 
*   corrections for the relativistic shift in the energy of the electron 
*   including non-diagonal corrections 
* 
* 
      DOUBLE PRECISION FUNCTION RLSHFT(I1,I2) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      PARAMETER(NOD=220,NWD=30) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
* 
*  *****  FORM DD  -L(L+1)/RR|P(I)> 
* 
      FL = L(I1) 
      C = (FL+D5)**2 
      LL = L(I1) + 1 
      L2 = 2*L(I1) + 1 
      L3 = 2*L(I1) + 3 
      ZZ = Z*Z 
      HH = 180.D0*H*H 
      MX = MAX0(MAX(I1), MAX(I2)) 
      DO 1 J = 2,MX 
1     YK(J) = -D1/RR(J) 
* 
*  *****  FORM THE INTEGRAND 
* 
      I = I1 
      A1 = D0 
      B1 = D0 
      DO 2 KK = 1,2 
      B2 = B1 
      YY = (P(3,I)+P(1,I) - D2*P(2,I))/(H*H) - C*P(2,I) 
      YK(2) = YY*YK(2) 
      YK(3) = YK(3)*((-(P(5,I)+P(1,I)) + D16*(P(4,I)+P(2,I)) 
     :        -D30*P(3,I))/(D12*H*H) - C*P(3,I)) 
      MM = MAX(I) - 3 
      DO 11 K =  4,MM 
      YY = D2*(P(K+3,I)+P(K-3,I)) 
      YY = YY - 27.D0*(P(K+2,I)+P(K-2,I)) 
      YY = YY +  270.D0*(P(K+1,I)+P(K-1,I)) - 490.D0*P(K,I) 
      YY = YY/HH - C*P(K,I) 
      YK(K) = YY*YK(K) 
      IF (K .EQ. 4) B1 = (YY/(D2*Z*P(4,I)*R(4)) + D1)/R(4) 
11    CONTINUE 
      MM = MM + 1 
      YK(MM) = YK(MM)*((-(P(MM+2,I)+P(MM-2,I))+D16*(P(MM+1,I)+P(MM-1,I)) 
     :      -D30*P(MM,I))/(D12*H*H) - C*P(MM,I)) 
      MM = MM + 1 
      YK(MM) = YK(MM)*((P(MM+1,I) + P(MM-1,I) - D2*P(MM,I))/(H*H) 
     :         - C*P(MM,I)) 
      A2 = A1 
      A1 = (P(1,I)/(AZ(I)*R(1)**L(I)*R2(1)) - D1 +Z*R(1)/LL)/RR(1) 
      I = I2 
2     CONTINUE 
* 
*  ***** DETERMINE CONTRIBUTION FROM NEAR THE NUCLEUS 
* 
      A = (Z/LL - L2*(B1 + B2)/D2)/LL 
      B = (L2*B1*B2 - D2*(A1 + A2) + (Z/LL**2)*(D2*Z*(D1 + D1/LL) 
     :     - L2*(B1 + B2)))/L3 
      RELSH = -P(4,I1)*P(4,I2)*(D1 + A*R(4) + B*RR(4))*D4*ZZ/L2 
      RELSH = RELSH/H1 - D5*YK(4) 
      RELSH2 = D0 
* 
*  *****  INTEGRATE 
* 
      DO 50 J = 5,MX,2 
      RELSH2 = RELSH2 + YK(J) 
50    RELSH = RELSH  + YK(J-1) 
      RELSH  = (RELSH + D2*RELSH2)*H1 
      IF ( L(I1) .EQ. 0 ) RELSH = RELSH + Z*AZ(I1)*AZ(I2) 
      RLSHFT = RELSH*D5*FINE 
      RETURN 
      END 
* 
*     ------------------------------------------------------------------ 
*              Y K F 
*     ------------------------------------------------------------------ 
* 
*               k 
*       Stores Y (i, j; r) in the array YK 
* 
* 
      SUBROUTINE YKF(I,J,K,REL) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220,NWD=30) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
      LOGICAL REL
      CALL ZK(I,J,K) 
      A =    EH**(K+1) 
      C = 2*K+1 
      A2 = A*A 
      H90 = C*H3/D30 
      A3 = A2*A*H90 
      AI = H90/A 
      AN = 114.D0*A*H90 
      A34 = 34.D0*H90 
      F1 = YK(NO)*EH**K 
      F2 = YK(NO) 
      F3 = YK(NO-1) 
      F4 = YK(ND) 
      DO 9 MM = 2,ND 
      M = NO -MM 
      F5 = YK(M-1) 
      YK(M) = YK(M+2)*A2 +     ( AN*F3 + A34*(F4+A2*F2)-F5*AI-F1*A3) 
      F1 = F2 
      F2 = F3 
      F3 = F4 
9     F4 = F5 
      YK(1) = YK(3)*A2+C*H3*(F4 + D4*A*F3 + A2*F2) 
      IF (.NOT.REL) RETURN
      MM = MAX0( MAX(I), MAX(J) ) 
      C = C*FINE 
      DO 10 M = 1,MM 
      YK(M) = YK(M) + C*P(M,I)*P(M,J) 
10    CONTINUE 
      RETURN 
      END 
*
*     -----------------------------------------------------------------
*       Z E T A
*     -----------------------------------------------------------------
*
*
      DOUBLE PRECISION FUNCTION ZETA(I1,I2) 
* 
*  ***** COMPUTES THE NUCLEAR SPIN-ORBIT PARAMETER AND THE 
*        CORRECTIONS FOR THE COMMON CLOSED SHELLS 
* 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      PARAMETER (NOD=220,NWD=30)
      COMMON/BLUME/COEFN2(4),COEFNK(4),COEFVK(4) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID 
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS 
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD), 
     :   YR(NOD),X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD) 
* 
      ZETA = FINE*Z*QUADR(I1,I2,-3) 
      LB = L(I1) 
      DO 10 I = 1,NCLOSD 
         LA = L(I) 
         ZETA = ZETA -(4*LA+2)*SN(I1, I, I2, I, 0) 
         CALL BWINT(LA,LB) 
         KE1 = 2 
         IF (LA .NE. LB) KE1 = IABS(LA-LB) 
         IP = 0 
         DO 20 K = KE1,LA+LB,2 
            IP = IP+1 
            ZETA = ZETA+COEFN2(IP)*SN(I1, I, I, I2, K-2) 
     :               +COEFNK(IP)*SN(I, I1, I2, I, K) 
     :               +COEFVK(IP)*(VK(I1,I,I,I2,K-1)-VK(I,I1,I2,I,K-1)) 
   20   CONTINUE 
   10 CONTINUE 
      END 
* 
*     ------------------------------------------------------------------ 
*              Z K 
*     ------------------------------------------------------------------ 
* 
*               k 
*       Stores Z (i, j; r) in the array YK. 
* 
* 
      SUBROUTINE ZK(I,J,K) 
      PARAMETER (NOD=220,NWD=30) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
      DEN = L(I) + L(J) + 3+ K 
      FACT = (D1/(L(I)+1) + D1/(L(J)+1))/(DEN + D1) 
      A = EH**K 
      A2 = A*A 
      H90 = H/90.D0 
      A3 = A2*A*H90 
      AI = H90/A 
      AN = 114.D0*A*H90 
      A34 = 34.D0*H90 
      F1 = RR(1)*P(1,I)*P(1,J) 
      F2 = RR(2)*P(2,I)*P(2,J) 
      F3 = RR(3)*P(3,I)*P(3,J) 
      F4 = RR(4)*P(4,I)*P(4,J) 
      YK(1) = F1*(D1 + Z*R(1)*FACT)/DEN 
      YK(2) = F2*(D1 + Z*R(2)*FACT)/DEN 
      YK(3) = YK(1)*A2 + H3*(F3 + D4*A*F2 + A2*F1) 
      DO 8 M = 5,NO 
      F5 = (RR(M)*P(M,I))*P(M,J) 
      YK(M-1) = YK(M-3)*A2 +    ( AN*F3 + A34*(F4+A2*F2)-F5*AI-F1*A3) 
      F1 = F2 
      F2 = F3 
      F3 = F4 
8     F4 = F5 
      YK(NO) = A*YK(NO-1) 
      IF (IABS(I-J)  +  IABS(K) .NE. 0) GO TO 2 
* 
*  *****  FOR Y0(I,I) SET THE LIMIT TO 1 AND REMOVE OSCILLATIONS 
*  *****  INTRODUCED BY THE USE OF SIMPSON'S RULE 
* 
      M1 = (NO/2)*2 - 1 
      M2 = M1 - 1 
      C1 = D1 - YK(M1) 
      C2 = D1 - YK(M2) 
      DO 3 M = 1,M1,2 
      YK(M) = YK(M) + C1 
3     YK(M+1) = YK(M+1) + C2 
      YK(NO) = D1 
      YK(NO-1) = D1 
2     RETURN 
      END 
