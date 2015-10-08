module bessel_mod

contains

  FUNCTION BESSJ0 (X)
    
    REAL *8 X,BESSJ0,AX,FR,FS,Z,FP,FQ,XX
    REAL *8 Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  &
         ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
    DATA P1,P2,P3,P4,P5 /1.D0,-.1098628627D-2,.2734510407D-4, &
         -.2073370639D-5,.2093887211D-6 /
    DATA Q1,Q2,Q3,Q4,Q5 /-.1562499995D-1,.1430488765D-3, &
         -.6911147651D-5,.7621095161D-6,-.9349451520D-7 /
    DATA R1,R2,R3,R4,R5,R6 /57568490574.D0,-13362590354.D0, &
         651619640.7D0,-11214424.18D0,77392.33017D0,-184.9052456D0 /
    DATA S1,S2,S3,S4,S5,S6 /57568490411.D0,1029532985.D0, &
         9494680.718D0,59272.64853D0,267.8532712D0,1.D0 /
    IF(X.EQ.0.D0) GO TO 1
    AX = ABS (X)
    IF (AX.LT.8.) THEN
       Y = X*X
       FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
       FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
       BESSJ0 = FR/FS
    ELSE
       Z = 8./AX
       Y = Z*Z
       XX = AX-.785398164
       FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
       FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
       BESSJ0 = SQRT(.636619772/AX)*(FP*COS(XX)-Z*FQ*SIN(XX))
    ENDIF
    RETURN
1   BESSJ0 = 1.D0
    RETURN
  END FUNCTION BESSJ0

!-----------------------------------------------------------------------

  FUNCTION BESSJ1 (X)

    REAL *8 X,BESSJ1,AX,FR,FS,Z,FP,FQ,XX
    REAL *8 Y,P1,P2,P3,P4,P5,P6,R1,R2,R3,R4,R5,R6  &
         ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
    DATA P1,P2,P3,P4,P5 /1.D0,.183105D-2,-.3516396496D-4,  &
         .2457520174D-5,-.240337019D-6 /,P6 /.636619772D0 /
    DATA Q1,Q2,Q3,Q4,Q5 /.04687499995D0,-.2002690873D-3,   &
         .8449199096D-5,-.88228987D-6,.105787412D-6 /
    DATA R1,R2,R3,R4,R5,R6 /72362614232.D0,-7895059235.D0, & 
         242396853.1D0,-2972611.439D0,15704.48260D0,-30.16036606D0 /
    DATA S1,S2,S3,S4,S5,S6 /144725228442.D0,2300535178.D0, &
         18583304.74D0,99447.43394D0,376.9991397D0,1.D0 /
    
    AX = ABS(X)
    IF (AX.LT.8.) THEN
       Y = X*X
       FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
       FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
       BESSJ1 = X*(FR/FS)
    ELSE
       Z = 8./AX
       Y = Z*Z
       XX = AX-2.35619491
       FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
       FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
       BESSJ1 = SQRT(P6/AX)*(COS(XX)*FP-Z*SIN(XX)*FQ)*SIGN(S6,X)
    ENDIF
    RETURN
  END FUNCTION BESSJ1

!-----------------------------------------------------------------------

  FUNCTION BESSJ (N,X)

    PARAMETER (IACC = 40,BIGNO = 1.D10, BIGNI = 1.D-10)
    REAL *8 X,BESSJ,TOX,BJM,BJ,BJP,SUM
    IF (N.EQ.0) THEN
       BESSJ = BESSJ0(X)
       RETURN
    ENDIF
    IF (N.EQ.1) THEN
       BESSJ = BESSJ1(X)
       RETURN
    ENDIF
    IF (X.EQ.0.) THEN
       BESSJ = 0.
       RETURN
    ENDIF
    TOX = 2./X
    IF (X.GT.FLOAT(N)) THEN
       BJM = BESSJ0(X)
       BJ  = BESSJ1(X)
       DO 11 J = 1,N-1
          BJP = J*TOX*BJ-BJM
          BJM = BJ
          BJ  = BJP
11        CONTINUE
       BESSJ = BJ
    ELSE
       M = 2*((N+INT(SQRT(FLOAT(IACC*N))))/2)
       BESSJ = 0.
       JSUM = 0
       SUM = 0.
       BJP = 0.
       BJ  = 1.
       DO 12 J = M,1,-1
          BJM = J*TOX*BJ-BJP
          BJP = BJ
          BJ  = BJM
          IF (ABS(BJ).GT.BIGNO) THEN
             BJ  = BJ*BIGNI
             BJP = BJP*BIGNI
             BESSJ = BESSJ*BIGNI
             SUM = SUM*BIGNI
          ENDIF
          IF (JSUM.NE.0) SUM = SUM+BJ
          JSUM = 1-JSUM
          IF (J.EQ.N) BESSJ = BJP
12        CONTINUE
       SUM = 2.*SUM-BJ
       BESSJ = BESSJ/SUM
    ENDIF
    RETURN
  END FUNCTION BESSJ

!-----------------------------------------------------------------------

  FUNCTION BESSY0 (X)
    REAL *8 X,BESSY0,FS,FR,Z,FP,FQ,XX
    REAL *8 Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  &
         ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
    DATA P1,P2,P3,P4,P5 /1.D0,-.1098628627D-2,.2734510407D-4, &
         -.2073370639D-5,.2093887211D-6 /
    DATA Q1,Q2,Q3,Q4,Q5 /-.1562499995D-1,.1430488765D-3,  &
         -.6911147651D-5,.7621095161D-6,-.9349451520D-7 /
    DATA R1,R2,R3,R4,R5,R6 /-2957821389.D0,7062834065.D0, &
         -512359803.6D0,10879881.29D0,-86327.92757D0,228.4622733D0 /
    DATA S1,S2,S3,S4,S5,S6 /40076544269.D0,745249964.8D0, &
         7189466.438D0,47447.26470D0,226.1030244D0,1.D0 /
    IF (X.EQ.0.D0) THEN
       BESSY0 = -1.E30
       RETURN
    ENDIF
    IF (X.LT.8.D0) THEN
       Y = X*X
       FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
       FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
       BESSY0 = FR/FS+.636619772D0*BESSJ0(X)*LOG(X)
    ELSE
       Z = 8.D0/X
       Y = Z*Z
       XX = X-.785398164D0
       FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
       FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
       BESSY0 = SQRT(.636619772D0/X)*(FP*SIN(XX)+Z*FQ*COS(XX))
    ENDIF
    RETURN
  END FUNCTION BESSY0

!-----------------------------------------------------------------------

  FUNCTION BESSY1 (X)
    REAL *8 X,BESSY1,FR,FS,Z,FP,FQ,XX
    REAL *8 Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  &
         ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6,S7
    DATA P1,P2,P3,P4,P5 /1.D0,.183105D-2,-.3516396496D-4, &
         .2457520174D-5,-.240337019D-6 /
    DATA Q1,Q2,Q3,Q4,Q5 /.04687499995D0,-.2002690873D-3,  &
         .8449199096D-5,-.88228987D-6,.105787412D-6 /
    DATA R1,R2,R3,R4,R5,R6 /-.4900604943D13,.1275274390D13,   &
         -.5153438139D11,.7349264551D9,-.4237922726D7,.8511937935D4 /
    DATA S1,S2,S3,S4,S5,S6,S7 /.2499580570D14,.4244419664D12, &
         .3733650367D10,.2245904002D8,.1020426050D6,.3549632885D3,1.D0 /
    IF (X.EQ.0.) THEN
       BESSY1 = -1.E30
       RETURN
    ENDIF
    IF (X.LT.8.) THEN
       Y = X*X
       FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
       FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*(S6+Y*S7)))))
       BESSY1 = X*(FR/FS)+.636619772*(BESSJ1(X)*LOG(X)-1./X)
    ELSE
       Z = 8./X
       Y = Z*Z
       XX = X-2.356194491
       FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
       FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
       BESSY1 = SQRT(.636619772/X)*(SIN(XX)*FP+Z*COS(XX)*FQ)
    ENDIF
    RETURN
  END FUNCTION BESSY1
  
!-----------------------------------------------------------------------

  FUNCTION BESSY (N,X)
      
    REAL *8 X,BESSY,TOX,BY,BYM,BYP
    IF (N.EQ.0) THEN
       BESSY = BESSY0(X)
       RETURN
    ENDIF
    IF (N.EQ.1) THEN
       BESSY = BESSY1(X)
       RETURN
    ENDIF
    IF (X.EQ.0.) THEN
       BESSY = -1.E30
       RETURN
    ENDIF
    TOX = 2./X
    BY  = BESSY1(X)
    BYM = BESSY0(X)
    DO 11 J = 1,N-1
       BYP = J*TOX*BY-BYM
       BYM = BY
       BY  = BYP
11     CONTINUE
    BESSY = BY
    RETURN
  END FUNCTION BESSY

!-----------------------------------------------------------------------

  REAL (KIND=8) FUNCTION BESSI0(X)
  
    REAL *8 X,Y,P1,P2,P3,P4,P5,P6,P7,  &
         Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
    DATA P1,P2,P3,P4,P5,P6,P7/1.D0,3.5156229D0,3.0899424D0,1.2067429D0,  &
         0.2659732D0,0.360768D-1,0.45813D-2/
    DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1, &
         0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,  &
         0.2635537D-1,-0.1647633D-1,0.392377D-2/
    IF(ABS(X).LT.3.75D0) THEN
       Y=(X/3.75D0)**2
       BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
    ELSE
       AX=ABS(X)
       Y=3.75D0/AX
       BX=EXP(AX)/SQRT(AX)
       AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
       BESSI0=AX*BX
    ENDIF
    RETURN
  END FUNCTION BESSI0

!-----------------------------------------------------------------------

  REAL (KIND=8) FUNCTION BESSI1(X)
  
    REAL *8 X,Y,P1,P2,P3,P4,P5,P6,P7,  &
         Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
    DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,  &
         0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
    DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1, &
         -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1, &
         -0.2895312D-1,0.1787654D-1,-0.420059D-2/
    IF(ABS(X).LT.3.75D0) THEN
       Y=(X/3.75D0)**2
       BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
    ELSE
       AX=ABS(X)
       Y=3.75D0/AX
       BX=EXP(AX)/SQRT(AX)
       AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
       BESSI1=AX*BX
    ENDIF
    RETURN
  END FUNCTION BESSI1

!-----------------------------------------------------------------------

  REAL (KIND=8) FUNCTION BESSI(N,X)

    REAL (KIND=8), PARAMETER :: BIGNO = 1.D10
    REAL (KIND=8), PARAMETER :: BIGNI = 1.D-10
    INTEGER (KIND=4), PARAMETER :: IACC = 40
    REAL *8 X,TOX,BIM,BI,BIP
    
    IF (N.EQ.0) THEN
       BESSI = BESSI0(X)
       RETURN
    ENDIF
    IF (N.EQ.1) THEN
       BESSI = BESSI1(X)
       RETURN
    ENDIF
    IF(X.EQ.0.D0) THEN
       BESSI=0.D0
       RETURN
    ENDIF
    TOX = 2.D0/X
    BIP = 0.D0
    BI  = 1.D0
    BESSI = 0.D0
    M = 2*((N+INT(SQRT(FLOAT(IACC*N)))))
    DO J = M,1,-1
       BIM = BIP+DFLOAT(J)*TOX*BI
       BIP = BI
       BI  = BIM
       IF (ABS(BI).GT.BIGNO) THEN
          BI  = BI*BIGNI
          BIP = BIP*BIGNI
          BESSI = BESSI*BIGNI
       ENDIF
       IF (J.EQ.N) BESSI = BIP
    END DO
    BESSI = BESSI*BESSI0(X)/BI
    RETURN
  END FUNCTION BESSI

!-----------------------------------------------------------------------

  REAL (KIND=8) FUNCTION BESSK0(X)

    REAL*8 X,Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7

    DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,0.23069756D0, &
         0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
    DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,0.2189568D-1, & 
         -0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/
    IF(X.EQ.0.D0) THEN
       BESSK0=1.D30
       RETURN
    ENDIF
    IF(X.LE.2.D0) THEN
       Y=X*X/4.D0
       AX=-LOG(X/2.D0)*BESSI0(X)
       BESSK0=AX+(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
    ELSE
       Y=(2.D0/X)
       AX=EXP(-X)/SQRT(X)
       BESSK0=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
    ENDIF
    RETURN
  END FUNCTION BESSK0

!-----------------------------------------------------------------------

  REAL (KIND=8) FUNCTION BESSK1(X)

    REAL*8 X,Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7

    DATA P1,P2,P3,P4,P5,P6,P7/1.D0,0.15443144D0,-0.67278579D0,  &
         -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
    DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,0.23498619D0,-0.3655620D-1, &
         0.1504268D-1,-0.780353D-2,0.325614D-2,-0.68245D-3/
    IF(X.EQ.0.D0) THEN
       BESSK1=1.D32
       RETURN
    ENDIF
    IF(X.LE.2.D0) THEN
       Y=X*X/4.D0
       AX=LOG(X/2.D0)*BESSI1(X)
       BESSK1=AX+(1.D0/X)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
    ELSE
       Y=(2.D0/X)
       AX=EXP(-X)/SQRT(X)
       BESSK1=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
    ENDIF
    RETURN
  END FUNCTION BESSK1

!-----------------------------------------------------------------------

  REAL (KIND=8) FUNCTION BESSK(N,X)
  
    REAL *8 X,TOX,BK,BKM,BKP
    
    IF (N.EQ.0) THEN
       BESSK = BESSK0(X)
       RETURN
    ENDIF
    IF (N.EQ.1) THEN
       BESSK = BESSK1(X)
       RETURN
    ENDIF
    IF (X.EQ.0.D0) THEN
       BESSK = 1.D30
       RETURN
    ENDIF
    TOX = 2.D0/X
    BK  = BESSK1(X)
    BKM = BESSK0(X)
    DO J=1,N-1
       BKP = BKM+DFLOAT(J)*TOX*BK
       BKM = BK
       BK  = BKP
    ENDDO
    BESSK = BK
    RETURN
  END FUNCTION BESSK

!-----------------------------------------------------------------------

end module bessel_mod
