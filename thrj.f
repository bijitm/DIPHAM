      FUNCTION THRJ(F1,F2,F3,G1,G2,G3)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  THIS FUNCTION CALCULATES THE 3J SYMBOL (F1 F2 F3)
C                                         (G1 G2 G3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE MUNG,X,Y
      DIMENSION X(502),Y(502)
      DATA MUNG/0/,MXIX/502/
      IF (MUNG.EQ.21) GOTO 69
      MUNG = 21
      X(1) = 0.D0
C
C  SET UP ARRAYS X(I)=LOG[(I-1)!] AND Y(I)=LOG(I-1)
      DO I = 1, MXIX-1
        A = I
        X(I+1) = LOG(A) +X(I)
        Y(I+1) = LOG(A)
      ENDDO

C  CRLS 06-2022: ARITHMETIC IFS REPLACED
C  CHECK THAT ALL M-VALUES ARE IN VALID RANGE
C  69 IF (F1-ABS(G1)) 1,13,13
   69 IF (F1.LT.ABS(G1)) GOTO 1
C  13 IF (F2-ABS(G2)) 1,14,14
      IF (F2.LT.ABS(G2)) GOTO 1
C  14 IF (F3-ABS(G3)) 1,15,15
      IF (F3.LT.ABS(G3)) GOTO 1
C  CHECK THAT SUM IS >= 0
      TOTAL=F1+F2+F3
      NTOTAL=TOTAL+.001D0
C     IF (TOTAL-NTOTAL) 2,2,1
      IF (TOTAL.LE.NTOTAL) GOTO 2
    1 THRJ=0.D0
      RETURN
C  CHECK SUM OF M-VALUES IS ZERO
C   2 IF (ABS(G1+G2+G3)-1.D-8) 3,3,1
    2 IF (ABS(G1+G2+G3).GT.1.D-8) GOTO 1
C  CHECK ALL TRIANGLE INEQUALITIES
C   3 IF (F1+F2-F3) 1,4,4
      IF (F1+F2.LT.F3) GOTO 1
C   4 IF (F1+F3-F2) 1,5,5
      IF (F1+F3.LT.F2) GOTO 1
C   5 IF (F2+F3-F1) 1,6,6
      IF (F2+F3.LT.F1) GOTO 1
C  SET UP SOME INTEGER VARIABLES
    6 J1=2.D0*F3+2.001D0
      J2=F1+F2-F3+1.001D0
      J3=F1-F2+F3+1.001D0
      J4=-F1+F2+F3+1.001D0
      J5=F1+F2+F3+2.001D0
      J6=F1+G1+1.001D0
      J7=F1-G1+1.001D0
      J8=F2+G2+1.001D0
      J9=F2-G2+1.001D0
      J10=F3+G3+1.001D0
      J11=F3-G3+1.001D0
      IF (J5.GT.MXIX) THEN
        WRITE(6,601) J5,MXIX
  601   FORMAT(' *** DIMENSION ERROR IN THRJ - INDEX.GT.MXIX',2I5)
        STOP
      ENDIF
C  GENERATE THE CLEBSCH-GORDAN COEFFICIENT (SEE, FOR EXAMPLE,
C  'ANGULAR MOMENTUM' BY ZARE EQN 2.25)
C
C  THIS IS THE PREFACTOR
      R=0.5D0*(Y(J1)+X(J2)+X(J3)+X(J4)-X(J5)+X(J6)+X(J7)+X(J8)+
     1         X(J9)+X(J10)+X(J11))
      TOTAL=0.D0
      F=-1.D0
      KZ=-1
C  LOOP LABELLED BY 7 DOES THE SUM IN THAT EXPRESSION
    7 KZ=KZ+1
      F=-F
      J1=KZ+1
      J2=F1+F2-F3-KZ+1.001D0
C     IF (J2) 20,20,8
      IF (J2.LT.1) GOTO 20
    8 J3=F1-G1-KZ+1.001D0
C     IF (J3) 20,20,9
      IF (J3.LT.1) GOTO 20
    9 J4=F2+G2-KZ+1.001D0
C     IF (J4) 20,20,10
      IF (J4.LT.1) GOTO 20
   10 J5=F3-F2+G1+KZ+1.001D0
C     IF (J5) 7,7,11
      IF (J5.LT.1) GOTO 7
   11 J6=F3-F1-G2+KZ+1.001D0
C     IF (J6) 7,7,12
      IF (J6.LT.1) GOTO 7
   12 JMAX=MAX(J1,J2,J3,J4,J5,J6)
      IF (JMAX.GT.MXIX) THEN
        WRITE(6,601) JMAX,MXIX
        STOP
      ENDIF
      S=-(X(J1)+X(J2)+X(J3)+X(J4)+X(J5)+X(J6))
      TOTAL=TOTAL+F*EXP(R+S)
      GOTO 7
C
C  CONVERT THE VALUE OF THE CLEBSCH-GORDAN CONSTANT TO THE RELEVANT 3J
C  SYMBOL
   20 INTVAL=ABS(F1-F2-G3)+0.0001D0
      VAL=((-1.D0)**INTVAL)*TOTAL/SQRT(2.D0*F3+1.D0)
      IF (ABS(VAL).LE.1.D-6) VAL=0.D0
      THRJ=VAL

      RETURN
      END
