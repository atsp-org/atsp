*     -----------------------------------------------------------------
      PROGRAM LSTERMS
*
*                   C O P Y R I G H T -- 1994
*
*      A. Hibbert, Dep't of Applied Mathematics and Theoretical Physics
*                  Queen's University of Belfast
*
*     From the program WEIGHTS
*     Computer Physics Communications, Vol. 1, 359--377 (1969).
*     -----------------------------------------------------------------
*
*    This program prints out a table of terms and seniority for the
*    the partially occupied shells of different l symmetries.

      CALL TMSOUT
      END

      BLOCK DATA
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/TERMS/NROWS,I(12),J(12),N(171)
C
C --- READS IN QUANTUM NUMBERS OF TERMS WHICH CAN BE FORMED FROM
C     CONFIGURATIONS  L**Q . ONLY THE FIRST HALF OF THAT PART OF THE
C     TABLE, CORRESPONDING TO A GIVEN  L, IS INCLUDED, BECAUSE OF THE
C     SYMMETRY OF THE TABLE.  E.G. D**7 FORMS THE SAME TERMS AS D**3
C
C     THE ARRAYS I,J,N CORRESPOND TO THE ARRAYS ITAB,JTAB,NTAB
C
      DATA NROWS/12/
      DATA I( 1),I( 2),I( 3),I( 4),I( 5),I( 6)/ 1, 1, 1, 3, 3, 1/
      DATA I( 7),I( 8),I( 9),I(10),I(11),I(12)/ 5, 8,16,16, 1, 1/
      DATA J( 1),J( 2),J( 3),J( 4),J( 5),J( 6)/  0,  3,  6,  9, 18, 27/
      DATA J( 7),J( 8),J( 9),J(10),J(11),J(12)/ 30, 45, 69,117,165,168/
      DATA N(  1),N(  2),N(  3),N(  4),N(  5),N(  6)/ 1, 1, 2, 0, 1, 1/
      DATA N(  7),N(  8),N(  9),N( 10),N( 11),N( 12)/ 1, 3, 2, 0, 1, 1/
      DATA N( 13),N( 14),N( 15),N( 16),N( 17),N( 18)/ 2, 5, 1, 2, 3, 3/
      DATA N( 19),N( 20),N( 21),N( 22),N( 23),N( 24)/ 1, 3, 2, 3, 5, 2/
      DATA N( 25),N( 26),N( 27),N( 28),N( 29),N( 30)/ 3, 1, 4, 1, 5, 2/
      DATA N( 31),N( 32),N( 33),N( 34),N( 35),N( 36)/ 0, 1, 1, 2, 5, 1/
      DATA N( 37),N( 38),N( 39),N( 40),N( 41),N( 42)/ 2, 9, 1, 2, 3, 3/
      DATA N( 43),N( 44),N( 45),N( 46),N( 47),N( 48)/ 2, 7, 3, 1, 5, 2/
      DATA N( 49),N( 50),N( 51),N( 52),N( 53),N( 54)/ 3, 3, 2, 3, 5, 2/
      DATA N( 55),N( 56),N( 57),N( 58),N( 59),N( 60)/ 3, 7, 2, 3, 9, 2/
      DATA N( 61),N( 62),N( 63),N( 64),N( 65),N( 66)/ 3,11, 2, 3, 3, 4/
      DATA N( 67),N( 68),N( 69),N( 70),N( 71),N( 72)/ 3, 7, 4, 0, 1, 1/
      DATA N( 73),N( 74),N( 75),N( 76),N( 77),N( 78)/ 2, 5, 1, 2, 9, 1/
      DATA N( 79),N( 80),N( 81),N( 82),N( 83),N( 84)/ 2, 3, 3, 2, 7, 3/
      DATA N( 85),N( 86),N( 87),N( 88),N( 89),N( 90)/ 4, 1, 1, 4, 5, 1/
      DATA N( 91),N( 92),N( 93),N( 94),N( 95),N( 96)/ 4, 7, 1, 4, 9, 1/
      DATA N( 97),N( 98),N( 99),N(100),N(101),N(102)/ 4,13, 1, 4, 3, 3/
      DATA N(103),N(104),N(105),N(106),N(107),N(108)/ 4, 5, 3, 4, 7, 3/
      DATA N(109),N(110),N(111),N(112),N(113),N(114)/ 4, 9, 3, 4,11, 3/
      DATA N(115),N(116),N(117),N(118),N(119),N(120)/ 4, 5, 5, 1, 5, 2/
      DATA N(121),N(122),N(123),N(124),N(125),N(126)/ 3, 3, 2, 3, 5, 2/
      DATA N(127),N(128),N(129),N(130),N(131),N(132)/ 3, 7, 2, 3, 9, 2/
      DATA N(133),N(134),N(135),N(136),N(137),N(138)/ 3,11, 2, 3, 3, 4/
      DATA N(139),N(140),N(141),N(142),N(143),N(144)/ 3, 7, 4, 5, 1, 2/
      DATA N(145),N(146),N(147),N(148),N(149),N(150)/ 5, 5, 2, 5, 7, 2/
      DATA N(151),N(152),N(153),N(154),N(155),N(156)/ 5, 9, 2, 5,13, 2/
      DATA N(157),N(158),N(159),N(160),N(161),N(162)/ 5, 5, 4, 5, 9, 4/
      DATA N(163),N(164),N(165),N(166),N(167),N(168)/ 5, 1, 6, 1, 7, 2/
      DATA N(169),N(170),N(171)                     / 1, 9, 2/
      END
      SUBROUTINE TMSOUT
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/TERMS/NROWS,ITAB(12),JTAB(12),NTAB(171)
      DIMENSION ICFG(12),LSYM(10),LLSYM(46)
      DATA ICFG,LSYM/4Hs(1),4Hs(2),4Hp(1),4Hp(2),4Hp(3),4Hd(1),4Hd(2),
     1               4Hd(3),4Hd(4),4Hd(5),4Hf(1),4Hg(1),
     2               1HS,1HP,1HD,1HF,1HG,1HH,1HI,1HK,1HL,1HM/
C
C --- PRINT-OUT OF TABLE OF TERMS, SET IN BLOCK DATA
C
 13   FORMAT( 1H1,23X,
     :  'PACKAGE - WEIGHTS'/
     :1X,'CALCULATION OF MATRIX ELEMENTS OF THE ONE- AND TWO-ELECTRON'/
     :10X,'PARTS OF THE HAMILTONIAN'//////21X
     :,' TABLE OF POSSIBLE TERMS'//
     : ' CONFIGURATION  TERMS (MULTIPLICITY, SYMMETRY, AND SENIORITY)'/)
    4 FORMAT(/1X,A4,4X,16(I1,A1,I1,1X))
      WRITE(6,13)
      DO 5 I=1,NROWS
      JI=JTAB(I)
      JJ=3*ITAB(I)
      DO 2 J=1,JJ,3
      LP=(NTAB(JI+J+1)+1)/2
    2 LLSYM(J)=LSYM(LP)
      WRITE(6,4) ICFG(I),(NTAB(JI+J+2),LLSYM(J),NTAB(JI+J),
     1 J=1,JJ,3)
    5 CONTINUE
      RETURN
      END
