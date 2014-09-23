C ***********************************************************************
C MODUL enthaelt Clebsch-Gordan-Koeffizienten
C  vor benutzen CALL fkltt im Hauptprogramm
C
C ***********************************************************************
      MODULE clebsch
       PRIVATE
     
       PARAMETER (NFAK=50)
       REAL *8 FFAK(0:NFAK)
       REAL *8 WFAK(0:NFAK)

       REAL,PUBLIC ::  DCH(0:NFAK),SDCH(0:NFAK),SD2CH(0:NFAK),
     X                 GHI(0:NFAK)

       PUBLIC C6J,C6J2,C9J,C9J2,CG,CG000,FKLTT,RM1H,FFAK,WFAK,
     X        DACH,GAMMAHI,NFAK


      CONTAINS



      FUNCTION CG(I,J,K,L,M,N)
C*******************************************************************************
C
C  CLEBSCH-GORDON   < I/2 J/2 L/2 M/2 / K/2 N/2 >
C
C*******************************************************************************
      INTEGER T
      IF (L+M .NE. N .OR. ABS(I-J) .GT. K .OR. K .GT. I+J .OR.
     .    ABS(L) .GT. I .OR. ABS(M) .GT. J .OR. ABS(N) .GT. K) THEN
       CG=0.
      RETURN
      ENDIF
      I1=(I+L)/2
      I2=(I-L)/2
      I3=(J+M)/2
      I4=(J-M)/2
      I5=(K-N)/2
      I6=(K+N)/2
      I7=(I+J-K)/2
      I8=(J+K-I)/2
      I9=(K+I-J)/2
      I10=(I+J+K+2)/2
      XX=WFAK(I7)*WFAK(I8)*WFAK(I9)*WFAK(I1)*WFAK(I2)*
     X     WFAK(I3)*WFAK(I4)*WFAK(I5)*
     .   WFAK(I6)/WFAK(I10)
      J1=(K-J+L)/2
      J2=(K-I-M)/2
      NT=MIN(I1,I2,I3,I4,I5,I6,I7,I8,I9)+1
      IT=0
      SUM=0.
      T=-1
   10   T=T+1
        L1=J1+T
        L2=J2+T
        L3=I7-T
        L4=I2-T
        L5=I3-T
        IF (MIN(L1,L2,L3,L4,L5) .LT. 0) GOTO 10
       SUM=SUM+RM1H(T)/(FFAK(T)*FFAK(L1)*FFAK(L2)*FFAK(L3)
     X       *FFAK(L4)*FFAK(L5))
       IT=IT+1
       IF (IT .LT. NT) GOTO 10
      CG=XX*SUM*SQRT(K+1.)
      END FUNCTION CG


      SUBROUTINE FKLTT
C*******************************************************************************
C
C  I!  SQRT(I!)   0 <= I <= NFAK
C
C*******************************************************************************
      FFAK(0)=1.0
      DO 10 I=1,NFAK
   10  FFAK(I)=I*FFAK(I-1)
      DO 20 I=0,NFAK
   20  WFAK(I)=SQRT(FFAK(I))
      END SUBROUTINE FKLTT



      FUNCTION RM1H (I)
C*******************************************************************************
C
C  RM1H = (-1) ** I
C
C*******************************************************************************
      RM1H=1-2*MOD(I,2)
      END FUNCTION RM1H




      END MODULE clebsch

