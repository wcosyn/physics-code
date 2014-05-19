CCC   PROGRAM GLOBAL.FOR  Written in FORTRAN - 77  18-AUG-1992
CCC   LAST CHANGE 10-SEP-2009
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C     This program calculates the Lorentz scalar and vector optical
C     potentials and/or the Schroedinger equivalent central and spin-orbit
C     potentials for Elastic Proton Scattering from Spin-0 Targets.
C
C     AUTHORS:  S. HAMA, E. D. COOPER, B. C. CLARK, R. L. MERCER
C
C     ADDRESS:  THE DEPARTMENT OF PHYSICS
C               THE OHIO STATE UNIVERSITY
C               174 WEST 18TH AVENUE
C               COLUMBUS, OHIO 43210-1106
C               U.S.A.
C
C     e-mail:  sn-hama@hue.ac.jp
C              tim.cooper@ufv.ca
C              bcc@mps.ohio-state.edu
C
C     REFERENCES:
C     (1) PHYS. REV. C 41, 2737 (1990)
C       " GLOBAL DIRAC OPTICAL POTENTIALS FOR ELASTIC PROTON SCATTERING
C         FROM HEAVY NUCLEI "
C
C     (2) PHYS. REV. C 47, 297 (1993)
C       " GLOBAL DIRAC PHENOMENOLOGY FOR PROTON-NUCLEUS ELASTIC SCATTERING "
C
C     (3) PHYS. REV. C 80, 034605 (2009)
C       " Global Dirac optical potential from helium to lead "
C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     The program provides tabulated scalar and vector 
C     and/or central and spin-orbit 
C     optical potentials for elastic proton scattering from
C     spin-0 nuclear targets (12C - 208Pb)
C     from 0.0 fm to RMAX fm with a step size given by DX fm's
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  **** Warning about extrapolation in energy and mass number ****
C            DO NOT DO IT !!!
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  **** Warning about Couloumb potentials and recoil corrections ****
C     The Coulomb potentials used in this Global Optical Potential fit
C     are calculated from the 2- / 3- parameter fermi forms of the
C     experimental charge density distributions, given in
C     ATOMIC DATA AND NUCLEAR DATA TABLES, VOL.36 NO.3 (1987) P.495
C     and VOL.14 NO.5-6 (1974) P.479.
C     To take care of recoil corrections an effctive Z is used given by
C                   ZEFF = Z * RECV
C     where Z is the charge of the target, RECV the recoil factor for vector
C     potentials. It will be given in the output.
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Recoil --- See E.D.Cooper and B.K.Jennings, Nucl.Phys. A483,(1988)601
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Glossary
C     -- Input ---
C        TPLAB -- Lab kinetic energy of the incident proton (MeV)
C        AA    -- The mass number of the target nucleus (amu)
C        DX    -- The step (mesh) size of tabulated potentials (fm)
C        RMAX  -- The maximum radius. (fm)
C     -- Output ---
C        EPCM  -- The total c.m. energy of incident proton (MeV)
C        E     -- Deffined as EPCM/1000 MeV
C        ETCM  -- The total c.m. energy of target nucleus (MeV)
C        WT    -- The mass of the target nucleus (MeV)
C        SR    -- The total c.m. energy of p-A system ( sqrt(s) ) (MeV)
C        RECV  -- The Cooper-Jennings recoil factor for vector potentials
C        RECS  -- The Cooper-Jennings recoil factor for scalar potentials
C        *****    NOTICE--The Cooper-Jennings recoil factors are included 
C                 in the tabulated scalar and vector potentials.
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(1000),VV(1000),WV(1000),VS(1000),WS(1000)
      DIMENSION VCEN(1000),WCEN(1000),VSO(1000),WSO(1000)
      CHARACTER*1 ANS
      CHARACTER*20 CPOTS
 5    WRITE(6,* ) ' ENTER 1 FOR DIRAC, 2 FOR SCHROEDINGER '
      READ(5,*)NDS
      IF(NDS.EQ.1) THEN
         WRITE(6,102)
 102     FORMAT(' ',/,
     $' **************************************************',/,
     $' GLOBAL DIRAC OPTICAL POTENTIALS FOR ELASTIC PROTON',/,
     $' SCATTERING FROM SPIN-0 TARGET NUCLEI  ',/,
     $' **************************************************',/)
      ELSE
         WRITE(6,600)
 600     FORMAT(' ',/,
     $' *******************************************************',/,
     $' SCHROEDINGER EQUIVALENT OPTICAL POTENTIALS CALCULATED ',/,
     $' FROM GLOBAL DIRAC OPTICAL POTENTIALS FOR ELASTIC PROTON',/,
     $' SCATTERING FROM SPIN-0 TARGET NUCLEI ',/,
     $' *******************************************************',/)
      WRITE(6,610)
 610  FORMAT(' ','  ** QUESTION ABOUT KINEMATICAL FACTORS ** ',/,
     &' THIS GIVES THE SCHROEDINGER EQUIVALENT POTENTIALS WHICH',
     &' SATISFY',/,
     &' ( DELSQ + K**2 - 2*F*( UCEN(R) + USO(R)*SIGDL))',
     &' PSI = 0 ',/,
     &' SIGDL   IS THE EXPECTATION VALUE OF SIGMA DOT L',
     &' AND = -(1+KAPPA).',/,
     &' K       IS C.M. MOMENTUM OF INCIDENT PROTON. ',/,
     &' USO(R)  IS THE SCHRODINGER EQUIVALENT SPIN-ORBIT POTENTIAL.',/,
     &' UCEN(R) IS THE SCHRODINGER EQUIVALENT CENTRAL POTENTIAL.',/,
     &' F       IS A KINEMATICAL FACTOR IN QUESTION. ',/,
     &' **************************************************',/)
      WRITE(6,620)
 620  FORMAT(' ','  ** NOTE ON COULOMB POTENTIALS ** ',/,
     &' VC      IS THE COULOMB POTENTIALS ',
     &' ( INCLUDING RECOIL FACTOR RECV ).',/,
     &' BY DEFAULT, UCEN(R) INCLUDES ALL CONTRIBUTIONS FROM VC ',/,
     &' EVEN FROM LINEAR PURE COULOMB TERM.',/,
     &' CHOICE FOR EXCLUDING LINEAR PURE COULOMB IS GIVEN LATER.',/,
     &' IN SCHROEDINGER EQUIVALENT SENSE, IT ENTERS AS (E/F)VC,',/,
     &' WHERE E IS C.M. TOTAL ENERGY OF INCIDENT PROTON AND ',/,
     &'       F IS A KINEMATICAL FACTOR IN QUESTION. ',/,
     &' **************************************************',/)
      END IF
CCC------------------------------------------------------------
CCC   CHOICES OF THE GLOBAL PARAMETERS
CCC------------------------------------------------------------
      WRITE(6,210)
 210  FORMAT(' ',' THE CHOICES ARE : ',/,
     &' ',' (1) CA40-PB208 (P,P) 65-1040 MEV  FIT.1 ',/,
     &' ',' (2) CA40-PB208 (P,P) 65-1040 MEV  FIT.2 ',/,
     &' ',' (3) C12 (P,P)   29-1040 MEV (EDAI C12)',/,
     &' ',' (4) O16 (P,P)   23-1040 MEV (EDAI O16)',/,
     &' ',' (5) CA40 (P,P)  21-1040 MEV (EDAI CA40)',/,
     &' ',' (6) ZR90 (P,P)  22-800  MEV (EDAI ZR90)',/,
     &' ',' (7) PB208 (P,P) 21-1040 MEV (EDAI PB208)',/,
     &' ',' (8) C12-PB208 (P,P) 21-1040 MEV (EDAD FIT.1)',/,
     &' ',' (9) C12-PB208 (P,P) 21-1040 MEV (EDAD FIT.2)',/,
     &' ','(10) C12-PB208 (P,P) 21-1040 MEV (EDAD FIT.3)',/,
     &' ','(11) HE4-PB208 (P,P) 21-1040 MEV (UNDEMOCRATIC-Fit2)',/,
     &' ','(12) HE4-PB208 (P,P) 21-1040 MEV (DEMOCRATIC-Fit1)')
      WRITE(6,107)
 107  FORMAT(' ',' YOUR CHOICE IS ( 1 - 12 ) ??')
      READ(5,*)IFIT
      WRITE(6,100)
 100  FORMAT(' ',' ENTER PROTON LAB K.E (MEV)')
      READ(5,*)TPLAB
      WRITE(6,103)
 103  FORMAT(' ',' ENTER AT. MASS NO. AND CHARGE OF TARGET ')
      READ(5,*)AA,ZT
 7    WRITE(6,105)
 105  FORMAT(' ',' ENTER MAX RADIUS-RMAX (FM) AND STEP SIZE-DX (FM) ')
      READ(5,*)RMAX,DX
C
      NSTEP= 1 + (RMAX+1.D-10)/DX
      IF (NSTEP.GT.1000) THEN
         WRITE(6,*)' (RMAX/DX+1) IS BIGGER THAN 1000, RE-ENTER '
         GO TO 7
      END IF
C
      WRITE(6,110)TPLAB
 110  FORMAT(' ',' TPLAB(MEV)=',F8.2)
      WRITE(6,112)AA,ZT
 112  FORMAT(' ',' AT. MASS (AMU)=',F8.3,'        CHARGE ZT =',F8.3)
      WRITE(6,115)RMAX,DX
 115  FORMAT(' ',' RMAX (FM) =',F8.2,'          DX (FM) =',F8.3)
      WRITE(6,117)IFIT
 117  FORMAT(' ',' FIT.',I2,' IS CHOSEN'/)
C
      DO 10 N=1,NSTEP
        R(N)=(N-1)*DX
 10   CONTINUE
      IF(NDS.EQ.1) THEN 
      CALL GLOBAL(IFIT,TPLAB,AA,R,VV,WV,VS,WS,NSTEP)
C
      WRITE(6,120)
      DO 15 N=1,NSTEP
        WRITE(6,125) R(N),VV(N),WV(N),VS(N),WS(N)
 15   CONTINUE
 120  FORMAT(' ',/,'    RADIUS,   RE. VECTOR   IM. VECTOR',
     & '   RE. SCALAR   IM. SCALAR')
 125  FORMAT(' ',F8.3,4F14.6)
C
      ELSE
      CALL SCHREQ(IFIT,TPLAB,AA,ZT,R,DX,VSO,WSO,VCEN,WCEN,NSTEP)
C
      WRITE(6,121)
      DO 16 N=2,NSTEP
        WRITE(6,125) R(N),VSO(N),WSO(N),VCEN(N),WCEN(N)
 16   CONTINUE
 121  FORMAT(' ',/,'    RADIUS,   RE. S-O      IM. S-O   ',
     & '   RE. CENT.    IM. CENT. ')
      END IF
C
      WRITE(6,128)
 128  FORMAT(' WOULD YOU LIKE NUMBERS WRITTEN TO A FILE? (Y/N)',/)
      READ(5,150) ANS
      IF(ANS.NE.'Y') GO TO 6
      WRITE(6,130)
 130  FORMAT(' ',' ON WHICH FILE NAME',
     & 'WOULD YOU LIKE THESE WRITTEN? (20 CHARACTERS MAX)')
      READ(5,*) CPOTS
        IUNIT=8
        OPEN(UNIT=IUNIT,FILE=CPOTS//'.TXT')
      IF(NDS.EQ.1)THEN
         DO 20 N=1,NSTEP
            WRITE(IUNIT,125) R(N),VV(N),WV(N),VS(N),WS(N) 
 20      CONTINUE
      ELSE
         DO 21 N=2,NSTEP
            WRITE(IUNIT,125) R(N),VSO(N),WSO(N),VCEN(N),WCEN(N)
 21      CONTINUE
      END IF
      CLOSE(UNIT=IUNIT)
   6  WRITE(6,140)
 140  FORMAT(' ','   DO YOU WISH TO CONTINUE ? (Y/N) ')
      READ(5,150) ANS
 150  FORMAT(1A)
      IF(ANS.EQ.'Y') GO TO 5
      IF(ANS.NE.'Y'.AND.ANS.NE.'N') GO TO 6
      STOP
      END
      SUBROUTINE GLOBAL(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(NSTEP),VVA(NSTEP),WVA(NSTEP),VSA(NSTEP),WSA(NSTEP)
CCC------------------------------------------------------------
CCC   CHOICES OF THE GLOBAL PARAMETERS
CCC------------------------------------------------------------
CCC   (1) CA40-PB208 (P,P) 65-1040 MEV  FIT.1 
CCC   (2) CA40-PB208 (P,P) 65-1040 MEV  FIT.2 
CCC   (3) C12 (P,P)   29-1040 MEV (EDAI C12)
CCC   (4) O16 (P,P)   23-1040 MEV (EDAI O16)
CCC   (5) CA40 (P,P)  21-1040 MEV (EDAI CA40)
CCC   (6) ZR90 (P,P)  22-800  MEV (EDAI ZR90)
CCC   (7) PB208 (P,P) 21-1040 MEV (EDAI PB208)
CCC   (8) C12-PB208 (P,P) 21-1040 MEV (EDAD FIT.1)
CCC   (9) C12-PB208 (P,P) 21-1040 MEV (EDAD FIT.2)
CCC  (10) C12-PB208 (P,P) 21-1040 MEV (EDAD FIT.3)
CCC  (11) 4HE-PB208 (P,P) 21-1040 MEV (UNDEMOCRATIC, TC)
CCC  (12) 4HE-PB208 (P,P) 21-1040 MEV (DEMOCRATIC, SH)
CCC------------------------------------------------------------
      GOTO (10,10,20,30,40,50,60,70,80,90,100,110),IFIT
 10   CALL GLOB12(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP)
      GOTO 99
 20   CALL EDAIC(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP)
      GOTO 99
 30   CALL EDAIO(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP)
      GOTO 99
 40   CALL EDAICA(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP) 
      GOTO 99
 50   CALL EDAIZR(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP)
      GOTO 99
 60   CALL EDAIPB(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP)
      GOTO 99
 70   CALL EDAD1(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP)
      GOTO 99
 80   CALL EDAD2(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP)
      GOTO 99
 90   CALL EDAD3(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP)
      GOTO 99
 100  CALL HEPB1(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP)
      GOTO 99
 110  CALL HEPB1(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP)      
 99   CONTINUE
      RETURN
      END
      SUBROUTINE GLOB12(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(NSTEP),VVA(NSTEP),WVA(NSTEP),VSA(NSTEP),WSA(NSTEP)
      IF(IFIT.NE.1.AND.IFIT.NE.2) THEN
           STOP ' ONLY IFIT = 1 AND 2 ARE CURRENTLY IMPLEMENTED'
      END IF
      IF(TPLAB.LT.65.0D0.OR.TPLAB.GT.1040.D0) WRITE(6,601)
 601  FORMAT(' ',' YOU ARE EXTRAPOLATING IN ENERGY !!?? '/)
      IF(AA.LT.40.D0.OR.AA.GT.208.D0) WRITE(6,602)
 602  FORMAT(' ',' YOU ARE EXTRAPOLATING IN MASS NUMBER !!?? '/)
C
CCC  
C      PI=3.141592653589793D0
      AMU=931.5016D0
      WT=AA*AMU
      EE=TPLAB
      ACB=AA**0.33333333333333333D0
CC
      IF(IFIT.EQ.1) THEN
          A=AA
      END IF
CC
      IF(IFIT.EQ.2) THEN
          A =AA
          A2=AA**(-2./3.)
          A3=AA**(1./3.)
      END IF
CC
CC
CC--------------------------------------------------------------
CC TRANSFORM TO 2 BODY CM FRAME AND CALCULATE RECOIL FACTORS
CC
        WP=1.0072545D0*AMU
C        CWT=AA**0.3333333333333333D0
        EL=EE+WP
        WP2=WP**2
        WT2=WT**2
        PCM=DSQRT(WT2*(EL**2-WP2)/(WP2+WT2+2.0D0*WT*EL))
        EPCM=DSQRT(WP2+PCM**2)
        ETCM=DSQRT(WT2+PCM**2)
        SR=EPCM+ETCM
C        EKINCM=SR-WP-WT
C        TCM=EKINCM
CC
        EPS=EPCM/1000.0
        E=EPS
        RECV=(ETCM /SR)
        RECS=(WT/SR)
        WRITE(6,111) E,ETCM,WT,SR,RECV,RECS
 111     FORMAT('  E,ETCM,WT,SR,RECV,RECS',/,1X,6(1PE10.3))
CC-------------------------------------------------------------------
      IF(IFIT.EQ.1) THEN
CC    Fit.1
CC    Real Vector Volume
       VV =  -0.593711E+02    +0.243991E+03/E     +0.422820E+03/E**2  
     &       -0.281139E+03/E**3
      RV1 =   0.164460E+01    -0.133264E+02/A     +0.334714E+03/A**2
     &       -0.163304E+01/E  +0.169182E+01/E**2  -0.505995E+00/E**3
      AV1 =   0.309170E+01    -0.942102E+01/A     +0.170789E+03/A**2
     &       -0.996093E+01/E  +0.140240E+02/E**2  -0.651832E+01/E**3
CC    Imag. vector Volume
       WV =  -0.187370E+03    +0.261183E+03/E     -0.247983E+03/E**2
     &       +0.127567E+03/E**3
      RV2 =   0.731599E+00    -0.378503E+01/A     +0.710760E+02/A**2
     &       +0.159654E+01/E  -0.223619E+01/E**2  +0.116266E+01/E**3
      AV2 =   0.100835E+01    -0.326173E+02/A     +0.968429E+03/A**2
     &       -0.870205E+00/E  +0.101693E+01/E**2  -0.375525E+00/E**3
CC    Real Scalar Volume
      VS =    0.144635E+03    -0.115015E+04/E     +0.864680E+03/E**2
     &       -0.281139E+03/E**3
      RS1 =   0.165499E+01    -0.133258E+02/A     +0.325457E+03/A**2
     &       -0.166097E+01/E  +0.174220E+01/E**2  -0.545308E+00/E**3
      AS1 =   0.211241E+01    -0.785078E+01/A     +0.133156E+03/A**2
     &       -0.603129E+01/E  +0.897855E+01/E**2  -0.440066E+01/E**3
CC    Imag. Scalar Volume
       WS =   0.767390E+02    +0.896685E+02/E     -0.247983E+03/E**2
     &       +0.127567E+03/E**3
      RS2 =   0.106401E+01    -0.233260E+01/A     +0.229048E+02/A**2
     &       +0.603736E-01/E  +0.537633E-01/E**2  +0.680720E-01/E**3
      AS2 =   0.646036E+00    -0.406376E+02/A     +0.132273E+04/A**2
     &       +0.160396E+00/E  +0.166333E+00/E**2  -0.278618E+00/E**3
CC
CC--------------------------------------------------------------------
CC    Imag. Vector Surface Peaked
      WVSP=    0.391326E+03 *EXP(  -0.265868E+01*E)
      RV3 =    0.109691E+01  +0.212733E+02/A  -0.100298E+04/A**2
      AV3 =    0.950649E+00  -0.531595E+02/A  +0.136272E+04/A**2 
CC    Imag. Scalar Surface Peaked
      WSSP=   -0.778311E+05 *EXP(  -0.726225E+01*E)
      RS3 =    0.117144E+01  +0.445545E+01/A  -0.418481E+03/A**2
      AS3 =    0.840840E+00  -0.434704E+02/A  +0.131545E+04/A**2 
CC
      END IF
CC
CC--------------------------------------------------------------------
CC
      IF(IFIT.EQ.2) THEN
CC    Fit.2
CC    Real Vector Volume 
       VV =   0.162663E+03    -0.428584E+03/E     +0.959310E+03/E**2  
     &       -0.374010E+03/E**3
      RV1 =   0.118178E+01    -0.123158E+00*A2    +0.234456E-01*A3
     &       -0.638310E+00/E  +0.514879E+00/E**2  -0.513080E-01/E**3
      AV1 =   0.367009E+01    +0.900712E-03*A     -0.125634E-05*A**2
     &       -0.125369E+02/E  +0.167826E+02/E**2  -0.742252E+01/E**3
CC    Imag. Vector Volume
       WV =  -0.949660E+03    +0.308267E+04/E     -0.374717E+04/E**2
     &       +0.156799E+04/E**3
       RV2 =  0.228599E+00    +0.332582E+01*A2    +0.528099E-01*A3
     &       +0.290166E+01/E  -0.453304E+01/E**2  +0.230291E+01/E**3
       AV2 =  0.189637E+01    -0.276306E-02*A     +0.131470E-04*A**2
     &       -0.513141E+01/E  +0.736379E+01/E**2  -0.346781E+01/E**3
CC    Real Scalar Volume
       VS =  -0.786320E+03    +0.215619E+04/E     -0.281098E+04/E**2
     &       +0.103025E+04/E**3
      RS1 =   0.968875E+00    -0.281597E+00*A2    +0.219828E-01*A3
     &       +0.227499E+00/E  -0.527801E+00/E**2  +0.346211E+00/E**3
      AS1 =   0.290644E+01    +0.881106E-03*A     -0.141390E-05*A**2
     &       -0.944671E+01/E  +0.128925E+02/E**2  -0.582424E+01/E**3
CC    Imag. Scalar Volume
       WS =   0.127600E+04    -0.438397E+04/E     +0.537798E+04/E**2
     &       -0.222344E+04/E**3
      RS2 =   0.205572E+01    +0.127931E+01*A2    +0.206096E-01*A3
     &       -0.454134E+01/E  +0.561357E+01/E**2  -0.206883E+01/E**3
      AS2 =   0.338803E+00    -0.388022E-02*A     +0.196852E-04*A**2
     &       +0.219169E+01/E  -0.274821E+01/E**2  +0.114420E+01/E**3
CC
CC-----------------------------------------------------------------------
CC    Imag. Vector Surface Peaked
      WVSP=   0.665325E+03  -0.122515E+04/E  +0.661619E+03/E**2
     &       -0.378568E+02/E**3
      RV3 =   0.817165E+00  +0.228100E+01*A2 +0.481430E-01*A3
      AV3 =   0.769104E+00  -0.259660E-02*A  +0.131676E-04*A**2 
CC    Imag. Scalar Surface Peaked
      WSSP=  -0.392815E+03  +0.887654E+03/E  -0.558651E+03/E**2
     &       +0.100805E+02/E**3
      RS3 =   0.820507E+00  +0.211267E+01*A2 +0.627332E-01*A3
      AS3 =   0.970102E+00  -0.337844E-02*A  +0.155953E-04*A**2 
CC
      END IF
CC
CC---------------------------------------------------------------------
CC
        DO 10 I=1,NSTEP
         TMP1 = R(I)
           E1=DEXP( (RV1*ACB-TMP1)/AV1)
           E2=DEXP(-(RV1*ACB+TMP1)/AV1)
         RVA1 = RECV*VV*E1/((1+E1)*(1+E2))
           E1=DEXP( (RV2*ACB-TMP1)/AV2)
           E2=DEXP(-(RV2*ACB+TMP1)/AV2)
           E3=DEXP( (RV3*ACB-TMP1)/AV3)
           E4=DEXP(-(RV3*ACB+TMP1)/AV3)
         RVSP=RECV*WVSP*E3*(1+2.*E4+E3*E4)/( ((1+E3)*(1+E4))**2 )
         RVA2 = RECV*WV*E1/((1+E1)*(1+E2)) + RVSP
           E1=DEXP( (RS1*ACB-TMP1)/AS1)
           E2=DEXP(-(RS1*ACB+TMP1)/AS1)
         RSA1 = RECS*VS*E1/((1+E1)*(1+E2))
           E1=DEXP( (RS2*ACB-TMP1)/AS2)
           E2=DEXP(-(RS2*ACB+TMP1)/AS2)
           E3=DEXP( (RS3*ACB-TMP1)/AS3)
           E4=DEXP(-(RS3*ACB+TMP1)/AS3)
         RSSP=RECS*WSSP*E3*(1+2.*E4+E3*E4)/( ((1+E3)*(1+E4))**2 )
         RSA2 = RECS*WS*E1/((1+E1)*(1+E2)) + RSSP
CC
         VVA(I)=RVA1
         WVA(I)=RVA2
         VSA(I)=RSA1
         WSA(I)=RSA2
 10    CONTINUE
      RETURN
      END
      SUBROUTINE SCHREQ(IFIT,TPLAB,AA,ZT,X,DX,USR,USI,UER,UEI,NSTEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(NSTEP),USR(NSTEP),USI(NSTEP),UER(NSTEP),UEI(NSTEP)
      DIMENSION VVA(1000),WVA(1000),VSA(1000),WSA(1000)
      DIMENSION AC(5000),DC(5000),DDC(5000)
      COMMON/RECOIL/RECV,RECS
CC
      HC1=197.32040D0
      HC2=HC1**2
CCC  
C      PI=3.141592653589793D0
      AMU=931.5016D0
      WT=AA*AMU
      EE=TPLAB
CC
CC--------------------------------------------------------------
CC TRANSFORM TO 2 BODY CM FRAME AND CALCULATE RECOIL FACTORS
CC
        WP=1.0072545D0*AMU
        AM=WP
        CWT=AA**0.3333333333333333D0
        EL=EE+WP
        WP2=WP**2
        WT2=WT**2
        PCM=DSQRT(WT2*(EL**2-WP2)/(WP2+WT2+2.0D0*WT*EL))
        EPCM=DSQRT(WP2+PCM**2)
        AE=EPCM
        ETCM=DSQRT(WT2+PCM**2)
        SR=EPCM+ETCM
CC
        RECV=(ETCM /SR)
        RECS=(WT/SR)
C--------------------------------------------------------------
C-- REDUCED ENERGY & REDUCED MASS , THEIR RATIOS, ETC.
        RENER = EPCM*ETCM/SR
        RMASS = WP*WT/(WP+WT)
        RAT1 = RENER/RMASS
        RAT2 = EPCM/RMASS
        RAT3 = EPCM/RENER
      WRITE(6,610)RENER,RMASS,RAT1,RAT2,RAT3
 610  FORMAT(' ',' REDUCED ENERGY (MEV)                 = ',D20.10,/,
     &'  REDUCED MASS (MEV)                   = ',D20.10,/,
     &'  RATIO (RED. ENERGY)/(RED. MASS)      = ',D20.10,/,
     &'  RATIO (C.M.PRO.ENERGY)/(RED. MASS)   = ',D20.10,/,
     &'  RATIO (C.M.PRO.ENERGY)/(RED. ENERGY) = ',D20.10,/)
C--------------------------------------------------------------
C     FORM OF SCHROEDINGER EQUATION. REDUCED MASS OR ENERGY OR WHAT?
CCCC
      ICHC=4
      WRITE(6,620)
 620  FORMAT(' ','PLEASE ENTER YOUR CHOICE OF KINEMATICAL FACTOR ?',/,
     &'  1 -- SET F=0.5 TO OBTAIN UCEN AND USO IN FM**(-2)  ',/,
     &'  2 -- REDUCED ENERGY ( FOR ECIS ) :  POT. GIVEN IN MEV ',/,
     &'  3 -- REDUCED MASS ( FOR RUNT ):  POT. GIVEN IN MEV  ',/,
     &'  4 -- C.M. PROJECTILE ENERGY:  POT. GIVEN IN MEV   ',/,
     &'  5 -- PROTON MASS:  POT. GIVEN IN MEV   ',/)
CCCC
      READ(5,510)ICHC
 510  FORMAT(I1)
        IF(ICHC.EQ.1) REDU=0.5*HC2
        IF(ICHC.EQ.2) REDU=EPCM*ETCM/SR
        IF(ICHC.EQ.3) REDU=WP*WT/(WP+WT)
        IF(ICHC.EQ.4) REDU=AE
        IF(ICHC.EQ.5) REDU=AM
CCCC
CXX    ICHC=1-- RECV*(2*AE)/HBARC           ON VC & ZT
CXX    ICHC=2-- RECV*AE/(REDUCED ENERGY)    ON VC & ZT
CXX    ICHC=3-- RECV*AE/(REDUCED MASS)      ON VC & ZT
CXX    ICHC=4-- RECV                        ON VC & ZT
CXX    ICHC=5-- RECV*AE/AM                  ON VC & ZT
CCCC
      IF(ICHC.EQ.1) FACOU=RECV*2.0*AE/HC1
      IF(ICHC.EQ.2) FACOU=RECV*AE/REDU
      IF(ICHC.EQ.3) FACOU=RECV*AE/REDU
      IF(ICHC.EQ.4) FACOU=RECV
      IF(ICHC.EQ.5) FACOU=RECV*AE/AM
CCCC
      IF(ICHC.EQ.1)WRITE(6,641) FACOU
      IF(ICHC.EQ.2)WRITE(6,642) FACOU
      IF(ICHC.EQ.3)WRITE(6,643) FACOU
      IF(ICHC.EQ.4)WRITE(6,644) FACOU
      IF(ICHC.EQ.5)WRITE(6,645) FACOU
 641  FORMAT(' ','SET F=0.5 TO OBTAIN UCEN AND USO IN FM**(-2)',/,
     &' SCALE FACTOR ON VC & ZT IS ',D20.10)
 642  FORMAT(' ','REDUCED ENERGY:  POT. GIVEN IN MEV ',/,
     &' SCALE FACTOR ON VC & ZT IS ',D20.10)
 643  FORMAT(' ','REDUCED MASS:  POT. GIVEN IN MEV   ',/,
     &' SCALE FACTOR ON VC & ZT IS ',D20.10)
 644  FORMAT(' ','C.M. PROJECTILE ENERGY:  POT. GIVEN IN MEV ',/,
     &' SCALE FACTOR ON VC & ZT IS ',D20.10)
 645  FORMAT(' ','PROTON MASS:  POT. GIVEN IN MEV ',/,
     &' SCALE FACTOR ON VC & ZT IS ',D20.10)
CCCC
      ICOU=1
      WRITE(6,630)
 630  FORMAT(' ',/,' PLEASE ENTER YOUR CHOICE OF COULOMB ?',/,
     &'  0 -- INCLUDE LINEAR + SQUARED COULOMB TERMS ',/,
     &'  1 -- EXCLUDE LINEAR + SQUARED COULOMB TERMS ',/,
     &'  2 -- EXCLUDE LINEAR BUT INCLUDE SQUARED COULOMB TERM ',/,
     &'  3 -- TURN OFF ALL COULOMB CONTRIBUTION ',/)
      READ(5,520)ICOU
 520  FORMAT(I1)
CCCC
      IF(ICOU.EQ.0)WRITE(6,650)
      IF(ICOU.EQ.1)WRITE(6,651)
      IF(ICOU.EQ.2)WRITE(6,652)
      IF(ICOU.EQ.3)WRITE(6,653)
 650  FORMAT(' ','INCLUDE LINEAR + SQUARED COULOMB TERMS ')
 651  FORMAT(' ','EXCLUDE LINEAR + SQUARED COULOMB TERMS ')
 652  FORMAT(' ','EXCLUDE LINEAR BUT INCLUDE SQUARED COULOMB TERM ')
 653  FORMAT(' ','TURN OFF ALL COULOMB CONTRIBUTION ',/)
C----------------------------------------------------------------- 
CCCCC   COULOMB
      COUF1=1.0
      COUF2=1.0
      IF(ICOU.EQ.1) THEN
         COUF1=0.0
         COUF2=0.0
      END IF
      IF(ICOU.EQ.2) THEN
         COUF1=0.0
         COUF2=1.0
      END IF
      DO 40 J=1,NSTEP
         AC(J)=0.
         DC(J)=0.
         DDC(J)=0.
  40  CONTINUE
      IF(ICOU.EQ.3) GO TO 42
      RV1=1.05
      AV1=0.67
      RCOU=RV1*CWT
      ACOU=AV1*2.0D0*0.56418958D0
      NC=20.0/DX+10
      IF (NC.GT.5000) STOP ' STOP DUE TO NC.GT.5000 '
      CALL COULP3
      CALL OSUCOU(AC,DC,DDC,DX,NSTEP,NC,ZT,WT,RCOU,ACOU)
  42  CONTINUE
CCC---RECOIL FACTOR RECV ON COULOMB IS TAKEN CARE IN OSUCOU
CCC---ERFCOU IS NOT AVAILABLE IN THIS
C----------------------------------------------------------------- 
      CALL GLOBAL(IFIT,TPLAB,AA,X,VVA,WVA,VSA,WSA,NSTEP)
CCC       DO I=1,NSTEP
CCC       WRITE(6,999)X(I),VVA(I),WVA(I),VSA(I),WSA(I)
CCC       END DO
CCC 999   FORMAT(' ',F8.3,4F14.6)
C ---------------------------------------------------------------
C     EXPRESSION OF US,UO,UT
C ---------------------------------------------------------------
      DO 50 I=1,NSTEP
      X(I)=DX*FLOAT(I-1)
      IF(I.EQ.1)X(I)=0.1D-9
      U1=VSA(I)
      W1=WSA(I)
      U2=VVA(I)
      W2=WVA(I)
      CALL DERI6(I,NSTEP,DX,VSA,UD1,UDD1)
      CALL DERI6(I,NSTEP,DX,WSA,WD1,WDD1)
      CALL DERI6(I,NSTEP,DX,VVA,UD2,UDD2)
      CALL DERI6(I,NSTEP,DX,WVA,WD2,WDD2)
CCC      WRITE(6,998)X(I),U1,W1,U2,W2
CCC      WRITE(6,998)X(I),UD1,WD1,UD2,WD2
CCC      WRITE(6,998)X(I),UDD1,WDD1,UDD2,WDD2
C 998  FORMAT(' ',F8.3,4E14.5)
C---------------------------------------------------------------
      AAB=AE+AM+U1-U2-AC(I)
      ADR=UD1-UD2-DC(I)
      ADI=WD1-WD2
      ADDR=UDD1-UDD2-DDC(I)
      ADDI=WDD1-WDD2
      AAI=W1-W2
C-------------------------------------------------------------
C     CALCULATION OF SPIN-ORBIT TERMS
      USR(I)=-0.5D0*HC2/REDU/X(I)*(ADR*AAB+ADI*AAI)/(AAB**2+AAI**2)
      USI(I)=-0.5D0*HC2/REDU/X(I)*(ADI*AAB-ADR*AAI)/(AAB**2+AAI**2)
C-------------------------------------------------------------
C     CALCULATION OF DARWIN TERMS
      UDR1=-1.0D0/X(I)*(ADR*AAB+ADI*AAI)/(AAB**2+AAI**2)
      UDR2=-0.5D0*(ADDR*AAB+ADDI*AAI)/(AAB**2+AAI**2)
      UDR3=0.75D0*((ADR**2-ADI**2)*(AAB**2-AAI**2)+4.0D0*ADR*ADI*
     1AAB*AAI)/((AAB**2-AAI**2)**2+4.0D0*AAB**2*AAI**2)
      UDRW=0.5D0*HC2/REDU*(UDR1+UDR2+UDR3)
      UDI1=-1.0D0/X(I)*(ADI*AAB-ADR*AAI)/(AAB**2+AAI**2)
      UDI2=-0.5D0*(ADDI*AAB-ADDR*AAI)/(AAB**2+AAI**2)
      UDI3=1.5D0*(ADR*ADI*(AAB**2-AAI**2)-(ADR**2-ADI**2)*AAB*AAI)
     1/((AAB**2-AAI**2)**2+4.0D0*AAB**2*AAI**2)
      UDIW=0.5D0*HC2/REDU*(UDI1+UDI2+UDI3)
C-------------------------------------------------------------
C     CALCULATION OF CENTRAL POTENTIAL TERMS (NO DARWIN)
      UCR1=2.0D0*AE*U2+2.0D0*AM*U1-U2**2+W2**2+U1**2
     1-W1**2-2.0D0*AC(I)*U2
      UCRW=0.5D0/REDU*UCR1
      UCI1=2.0D0*(AE*W2+AM*W1-U2*W2+U1*W1-AC(I)*W2)
      UCIW=0.5D0/REDU*UCI1
C   ------------------------------------------------------
C     CALCULATION OF CENTRAL POTENTIAL (WITH DARWIN)
      UER(I)=UCRW+UDRW-COUF2*0.5D0*AC(I)**2/REDU 
     &     + COUF1*(AE/REDU)*AC(I)
      UEI(I)=UCIW+UDIW
CCC      WRITE(6,333)X(I),USR(I),USI(I),UER(I),UEI(I)
C  333 FORMAT(F7.3,4(1PD16.8))
   50 CONTINUE
C    --------------------------------------------------------
      RETURN
      END
      SUBROUTINE COULP3
C     SUBROUTINE COULP3   ( 24-JUN-1992 : LAST CHANGE )
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 CLP(3,92,238)
      COMMON/PARA3/CLP
CC==  SOURCES OF 2- AND 3-PARAMETER FERMI CHARGE DENSITY DISTRIBUTION ---
CC    ATOMIC DATA AND NUCLEAR DATA TABLES, VOL.36 NO.3 (1987) P.495
CC    ATOMIC DATA AND NUCLEAR DATA TABLES, VOL.14 NO.5-6 (1974) P.479
CC    PHYS.REV. 174, (1968) 1380
CC-------------------------------------
C   --4HE------------------
      CLP(1,2,4)=0.445
      CLP(2,2,4)=1.008
      CLP(3,2,4)=0.327
C   --12C-----------------
      CLP(1,6,12)=-0.149
      CLP(2,6,12)=2.355
      CLP(3,6,12)=0.5224
C   --14N------------------
      CLP(1,7,14)=-0.180
      CLP(2,7,14)=2.570
      CLP(3,7,14)=0.5052
C   ---15N-------------------   
      CLP(1,7,15)=0.139
      CLP(2,7,15)=2.334
      CLP(3,7,15)=0.498
C   --16O---------------------   
      CLP(1,8,16)=-0.051
      CLP(2,8,16)=2.608
      CLP(3,8,16)=0.513
C   ---19F------------------   
      CLP(1,9,19)=0.0
      CLP(2,9,19)=2.58
      CLP(3,9,19)=0.567
C   ---20NE--------------------   
      CLP(1,10,20)=-0.168
      CLP(2,10,20)=2.791
      CLP(3,10,20)=0.698
C   ----22NE-------------------   
      CLP(1,10,22)=0.0
      CLP(2,10,22)=2.782
      CLP(3,10,22)=0.549
C   -----24MG------------------   
      CLP(1,12,24)=-0.249
      CLP(2,12,24)=3.192
      CLP(3,12,24)=0.604
C   ----25MG-------------------   
      CLP(1,12,25)=-0.2360
      CLP(2,12,25)=3.22
      CLP(3,12,25)=0.58
C   ----26MG-------------------   
      CLP(1,12,26)=0.0
      CLP(2,12,26)=3.05
      CLP(3,12,26)=0.523
C   ----27AL-------------------   
      CLP(1,13,27)=0.0
      CLP(2,13,27)=3.07
      CLP(3,13,27)=0.519
C   ----28SI-------------------   
      CLP(1,14,28)=-0.233
      CLP(2,14,28)=3.340
      CLP(3,14,28)=0.580
C   ----29SI-------------------   
      CLP(1,14,29)=-0.203
      CLP(2,14,29)=3.338
      CLP(3,14,29)=0.547
C   ----30SI-------------------   
      CLP(1,14,30)=-0.078
      CLP(2,14,30)=3.252
      CLP(3,14,30)=0.553
C   ----31P-------------------   
      CLP(1,15,31)=-0.173
      CLP(2,15,31)=3.369
      CLP(3,15,31)=0.582
C   ----35CL-------------------   
      CLP(1,17,35)=-0.10
      CLP(2,17,35)=3.476
      CLP(3,17,35)=0.599
C   ----37CL-------------------   
      CLP(1,17,37)=-0.13
      CLP(2,17,37)=3.554
      CLP(3,17,37)=0.588
C   ----36AR-------------------   
      CLP(1,18,36)=-0.0
      CLP(2,18,36)=3.54
      CLP(3,18,36)=0.507
C   ----40AR-------------------   
      CLP(1,18,40)=-0.19
      CLP(2,18,40)=3.73
      CLP(3,18,40)=0.62
C   ----39K-------------------   
      CLP(1,19,39)=-0.201
      CLP(2,19,39)=3.743
      CLP(3,19,39)=0.585
C   ----40CA-------------------   
      CLP(1,20,40)=-0.161
      CLP(2,20,40)=3.766
      CLP(3,20,40)=0.586
C   ----42CA-------------------   
      CLP(1,20,42)=-0.1158
      CLP(2,20,42)=3.7278
      CLP(3,20,42)=0.5911
C   ----44CA-------------------   
      CLP(1,20,44)=-0.0948
      CLP(2,20,44)=3.7481
      CLP(3,20,44)=0.5715
C   ----48CA-------------------   
      CLP(1,20,48)=-0.030
      CLP(2,20,48)=3.7369
      CLP(3,20,48)=0.5245
C   ----46TI--(from difference between isotopes )-----------------   
      CLP(1,22,46)=-0.0
      CLP(2,22,46)=3.788
      CLP(3,22,46)=0.613
C   ----48TI-------------------   
      CLP(1,22,48)=-0.0
      CLP(2,22,48)=3.843
      CLP(3,22,48)=0.588
C   ----50TI--(from difference between isotopes )-----------------   
      CLP(1,22,50)=-0.0
      CLP(2,22,50)=3.888
      CLP(3,22,50)=0.563
C   ----51V-------------------   
      CLP(1,23,51)=-0.0
      CLP(2,23,51)=3.94
      CLP(3,23,51)=0.505
C   ----50CR-------------------   
      CLP(1,24,50)=-0.0
      CLP(2,24,50)=3.941
      CLP(3,24,50)=0.566
C   ----52CR-------------------   
      CLP(1,24,52)=-0.0
      CLP(2,24,52)=3.984
      CLP(3,24,52)=0.542
C   ----53CR-------------------   
      CLP(1,24,53)=-0.0
      CLP(2,24,53)=4.000
      CLP(3,24,53)=0.557
C   ----54CR-------------------   
      CLP(1,24,54)=-0.0
      CLP(2,24,54)=4.010
      CLP(3,24,54)=0.578
C   ----55MN-------------------   
      CLP(1,25,55)=-0.0
      CLP(2,25,55)=3.89
      CLP(3,25,55)=0.567
C   ----54FE-------------------   
      CLP(1,26,54)=-0.0
      CLP(2,26,54)=4.074
      CLP(3,26,54)=0.536
C   ----56FE-------------------   
      CLP(1,26,56)=-0.0
      CLP(2,26,56)=4.111
      CLP(3,26,56)=0.558
C   ----58FE-------------------   
      CLP(1,26,58)=-0.0
      CLP(2,26,58)=4.027
      CLP(3,26,58)=0.576
C   ----59CO-------------------   
      CLP(1,27,59)=-0.0
      CLP(2,27,59)=4.158
      CLP(3,27,59)=0.575
C   ----58NI-------------------   
      CLP(1,28,58)=-0.1308
      CLP(2,28,58)=4.3092
      CLP(3,28,58)=0.5169
C   ----60NI-------------------   
      CLP(1,28,60)=-0.2668
      CLP(2,28,60)=4.4891
      CLP(3,28,60)=0.5369
C   ----61NI-------------------   
      CLP(1,28,61)=-0.1983
      CLP(2,28,61)=4.4024
      CLP(3,28,61)=0.5401
C   ----62NI-------------------   
      CLP(1,28,62)=-0.2090
      CLP(2,28,62)=4.4425
      CLP(3,28,62)=0.5386
C   ----64NI-------------------   
      CLP(1,28,64)=-0.2284
      CLP(2,28,64)=4.5211
      CLP(3,28,64)=0.5278
C   ----89Y-------------------   
      CLP(1,39,89)=-0.0
      CLP(2,39,89)=4.86
      CLP(3,39,89)=0.542
C   ----90ZR-------( from NP A179, 529 (1972))------------   
      CLP(1,40,90)=-0.086
      CLP(2,40,90)=4.86
      CLP(3,40,90)=0.57
C   ----124SN-------------------   
      CLP(1,50,124)=0.
      CLP(2,50,124)=5.490
      CLP(3,50,124)=0.534
C   ----154SM-------------------   
      CLP(1,62,154)=0.
      CLP(2,62,154)=5.9387
      CLP(3,62,154)=0.522
C   ----176YB-------------------   
      CLP(1,70,176)=0.
      CLP(2,70,176)=6.127
      CLP(3,70,176)=0.363
C   ----208PB--------( PRL 23, 1402 (1969))-----------   
      CLP(1,82,208)=0.32
      CLP(2,82,208)=6.4
      CLP(3,82,208)=0.54
C   ----209BI-------------------   
      CLP(1,83,209)=0.
      CLP(2,83,209)=6.75
      CLP(3,83,209)=0.468
C   ----238U-------------------   
      CLP(1,92,238)=-0.0
      CLP(2,92,238)=6.874
      CLP(3,92,238)=0.556
C===================================================
C-------------------------------------------------
      RETURN
      END
      SUBROUTINE OSUCOU(AC,DC,DDC,DMESH,NPTS,NC,ZT,WT,RCOU,ACOU)
C      SUBROUTINE OSUCOU(AC,DMESH,NPTS,NDAT)
C ----------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AC(NC),DC(NC),DDC(NC)
      REAL*8 CLP(3,92,238)
      COMMON/PARA3/CLP
      COMMON/RECOIL/RECV,RECS
      COMMON W,C,Z
      DATA AMU /931.5016D0/
      EXTERNAL FCI,FCJ
C------------------------------------------------
      IAA=INT(WT/AMU + 0.0002)
      W=CLP(1,INT(ZT),IAA)
      C=CLP(2,INT(ZT),IAA)
      Z=CLP(3,INT(ZT),IAA)
      IF(DABS(C).LT.1.0D-33) THEN
         WRITE(6,888)
C----- USE W=0.0, C=1.05*A**.3 AND Z=0.65 INSTEAD
         W=0.0D0
         C=1.05D0*IAA**0.33333333333333333D0
         Z=0.65D0
         WRITE(6,889)W,C,Z
      END IF
 888  FORMAT(' ',' *** OSUCOU: CHARGE DENSITY IS NOT AVAILABLE ****')
 889  FORMAT(' ',' *** CHARGE DENSITY ASSUMED: W, C, Z = ',3F8.3,/)
CCC      IF(C.EQ.0.0D0) GO TO 877
C--------------------------------------------------
      ATNO=RECV*ZT
      CP=1.0D0
C -------------------------------------------------------------
      HC1=197.32040D0
C      HC2=HC1**2
C -------------------------------------------------------------
C     CALCULATION OF COULOMB POTENTIAL AND ITS NORMALIZATION
CCC            CALCULATED UP TO 20 FM
      SUMI=0.0D0
      SUMJ=0.0D0
      I=2
      XL=0.0D0
      XU=20.0D0
      DX=DMESH
      XB=XL+DX
    2 CALL DGI16(XL,XB,FCI,Y)
      CALL DGI16(XL,XB,FCJ,Y2)
      SUMI=SUMI+Y
      SUMJ=SUMJ+Y2
      AC(I)=SUMI/XB-SUMJ
      DC(I)=-SUMI/XB**2+FCI(XB)/XB-FCJ(XB)
      DDC(I)=(2.0D0*SUMI/XB-FCI(XB))/XB**2
      I=I+1
      IF(XB.GE.XU)GO TO 90
      XL=XB
      XB=XB+DX
      GO TO 2
 90   CONTINUE
CH   90 WRITE(6,100)SUMI,SUMJ
CH  100 FORMAT(1H ,'NORMALIZATION INTEGRALS(COULOMB) I=',D18.10,5X,
CH     1'J=',D18.10/)
      GAM=0.00729729D0*ATNO*HC1/SUMI*CP
      DO 20 I=2,NPTS
      AC(I)=(AC(I)+SUMJ)*GAM
      DC(I)=DC(I)*GAM
      DDC(I)=DDC(I)*GAM
   20 CONTINUE
      I=1
      XL=0.0D0
      XB=0.1D-9
      CALL DGI16(XL,XB,FCI,Y0)
      CALL DGI16(XL,XB,FCJ,Y1)
      AC(I)=Y0/XB-Y1
      AC(I)=(AC(I)+SUMJ)*GAM
      DC(I)=(-Y0/XB**2+FCI(XB)/XB-FCJ(XB))*GAM
      DDC(I)=(2.0D0*Y0/XB-FCI(XB))/XB**2*GAM
CH-------- AC(I) CHANGE TO IN UNIT OF FM*-1
CH      DO 30 J=1,NPTS
CH      AC(J)=AC(J)/HC1
CH30    CONTINUE
         GO TO 876
cx877      CONTINUE
cx         STOP 
CHH      CALL ERFCOU(AC,ZTARG,RCOU,ACOU,DMESH,NPTS)
876      CONTINUE
      RETURN
      END
      FUNCTION FCI(R)
C     FUNCTION FCI(R)
C ----------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON W,C,Z
      E1=1.0D0+W*(R/C)**2
      E2=1.0D0+DEXP((R-C)/Z)
      FCI=(E1/E2)*R**2
      RETURN
      END
      FUNCTION FCJ(R)
C     FUNCTION FCJ(R)
C ----------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON W,C,Z
      E1=1.0D0+W*(R/C)**2
      E2=1.0D0+DEXP((R-C)/Z)
      FCJ=(E1/E2)*R
      RETURN
      END
      SUBROUTINE DGI16(XL,XU,FCT,Y)
C     SUBROUTINE DGI16(XL,XU,FCT,Y)
C   --------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
C
      A=.5D0*(XU+XL)
      B=XU-XL
      C=.49470046749582497D0*B
      Y=.13576229705877047D-1*(FCT(A+C)+FCT(A-C))
      C=.47228751153661629D0*B
      Y=Y+.31126761969323946D-1*(FCT(A+C)+FCT(A-C))
      C=.43281560119391587D0*B
      Y=Y+.47579255841246392D-1*(FCT(A+C)+FCT(A-C))
      C=.37770220417750152D0*B
      Y=Y+.62314485627766936D-1*(FCT(A+C)+FCT(A-C))
      C=.30893812220132187D0*B
      Y=Y+.7479799440828837D-1*(FCT(A+C)+FCT(A-C))
      C=.22900838882861369D0*B
      Y=Y+.8457825969750127D-1*(FCT(A+C)+FCT(A-C))
      C=.14080177538962946D0*B
      Y=Y+.9130170752246179D-1*(FCT(A+C)+FCT(A-C))
      C=.47506254918818720D-1*B
      Y=B*(Y+.9472530522753425D-1*(FCT(A+C)+FCT(A-C)))
      RETURN
      END
      SUBROUTINE VINT1(M,VSUM,V,R,DR)
C     SUBROUTINE VINT1(M,VSUM,V,R,DR)
C   ------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION V(M),R(M)
      PI=3.141592654D0
      SUM1=0.0D0
      SUM2=0.0D0
      N1=M-1
      N2=M-2
      DO 10 I=2,N1,2
   10 SUM1=SUM1+4.0D0*V(I)*R(I)**2
      DO 20 I=3,N2,2
   20 SUM2=SUM2+2.0D0*V(I)*R(I)**2
      SUM3=V(M)*R(M)**2
      VSUM=4.0D0*PI*DR/3.0D0*(SUM1+SUM2+SUM3)
      RETURN
      END
      SUBROUTINE DERI6(I,N,DX,V,UD,UDD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V(1000)
      IF(I.GE.3) GO TO 1110
      UD=(-274.0D0*V(I)+600.0D0*V(I+1)-600.0D0*V(I+2)+
     &400.0D0*V(I+3)-150.0D0*V(I+4)+24.0D0*V(I+5))/(120.D0*DX)
      UDD=(225.D0*V(I)-770.D0*V(I+1)+1070.D0*V(I+2)-
     &780.D0*V(I+3)+305.D0*V(I+4)-50.D0*V(I+5))/(60.D0*DX**2)
      GO TO 1130
 1110 IF(I.LE.N-3) GO TO 1120
      UD=(-24.0D0*V(I-5)+150.0D0*V(I-4)-400.0D0*V(I-3)+
     &600.0D0*V(I-2)-600.0D0*V(I-1)+274.0D0*V(I))/(120.0D0*DX)
      UDD=(-50.D0*V(I-5)+305.D0*V(I-4)-780.D0*V(I-3)+
     &1070.D0*V(I-2)-770.D0*V(I-1)+225.D0*V(I))/(60.D0*DX**2)
      GO TO 1130
 1120 CONTINUE
      UD=(6.0D0*V(I-2)-60.0D0*V(I-1)-40.0D0*V(I)+120.0D0*
     &V(I+1)-30.0D0*V(I+2)+4.0D0*V(I+3))/(120.D0*DX)
      UDD=(-5.0D0*V(I-2)+80.0D0*V(I-1)-150.0D0*V(I)+
     &80.0D0*V(I+1)-5.0D0*V(I+2))/(60.0D0*DX**2)
 1130 CONTINUE
      RETURN
      END
      SUBROUTINE EDAIC(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(NSTEP),VVA(NSTEP),WVA(NSTEP),VSA(NSTEP),WSA(NSTEP)
      DIMENSION P(22,8),PT1(8,8),PT2(8,8),PT3(8,6)
C23456789012345678901234567890123456789012345678901234567890123456789012
CC    COSH.SPV 
CC     VOLUME = (COSH(R/A)-1)/(COSH(R/A)+COSH(X/A)-2)
CC     S.P = (COSH(R/A)-1)*(COSH(X/A)-1)/((COSH(R/A)+COSH(X/A)-2)**2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   C12_K.UNIT12 
       DATA PT1/
     & -4.501206E+00,  2.445517E+01, -2.616971E+01, -3.655777E+00,
     &  1.063368E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  3.210640E+00, -1.324620E+01,  2.712683E+01, -2.330649E+01,
     &  7.152633E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.306250E-02,  4.569257E+00, -8.858008E+00,  8.123912E+00,
     & -3.202942E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -5.006626E+02,  2.371955E+03, -4.161205E+03,  3.186804E+03,
     & -8.980641E+02,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.258661E+01, -5.940515E+01,  1.082025E+02, -8.419066E+01,
     &  2.341066E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -3.216678E+00,  2.527721E+01, -5.200693E+01,  4.221788E+01,
     & -1.121544E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.840370E+01, -9.793797E+01,  2.015341E+02, -1.829884E+02,
     &  6.194047E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  4.572097E+00, -2.093001E+01,  4.208205E+01, -3.541818E+01,
     &  1.061756E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00/
       DATA PT2/
     &  2.274882E+00, -6.968236E+00,  1.381223E+01, -1.181611E+01,
     &  3.423874E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.425824E+01,  1.647004E+01, -2.044816E+02,  3.398014E+02,
     & -1.656334E+02,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -4.155769E+00,  2.006808E+01, -2.893194E+01,  1.852043E+01,
     & -4.867382E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -1.791053E+01,  1.143694E+02, -2.449500E+02,  2.217276E+02,
     & -7.232253E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -8.572854E+01,  4.852726E+02, -9.926830E+02,  8.809676E+02,
     & -2.847979E+02,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  7.000000E-01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  6.409709E+01, -2.371326E+02,  2.649673E+02, -6.775234E+01,
     & -2.120731E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00/
       DATA PT3/ 
     &  1.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  7.000000E-01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00/
CC
      DSECH(X) = 2.*DEXP(-DABS(X))/( 1.0D0 + DEXP(-DABS(X))**2 )
CCC
      IF(TPLAB.LT.29.0D0.OR.TPLAB.GT.1040.D0) WRITE(6,601)
 601  FORMAT(' ',' YOU ARE EXTRAPOLATING IN ENERGY !!?? '/)
      IF(DABS(AA-12.0D0).GT.1.0D-33) WRITE(6,602)
 602  FORMAT(' ',' YOU ARE EXTRAPOLATING IN MASS NUMBER !!?? '/)
CCC
      DO 20 I=1,8
         DO 30 J=1,8
            P(J,I)=PT1(I,J)
 30      CONTINUE
         DO 40 J=1,8
            P(8+J,I)=PT2(I,J)
 40      CONTINUE
         DO 50 J=1,6
            P(16+J,I)=PT3(I,J)
 50      CONTINUE
 20   CONTINUE
CCC
C      PI=3.141592653589793D0
      AMU=931.5016D0
      WT=AA*AMU
      EE=TPLAB
      ACB=AA**0.33333333333333333D0
CC
CC
CC--------------------------------------------------------------
CC TRANSFORM TO 2 BODY CM FRAME AND CALCULATE RECOIL FACTORS
CC
        WP=1.0072545D0*AMU
C        CWT=AA**0.3333333333333333D0
        EL=EE+WP
        WP2=WP**2
        WT2=WT**2
        PCM=DSQRT(WT2*(EL**2-WP2)/(WP2+WT2+2.0D0*WT*EL))
        EPCM=DSQRT(WP2+PCM**2)
        ETCM=DSQRT(WT2+PCM**2)
        SR=EPCM+ETCM
C        EKINCM=SR-WP-WT
C        TCM=EKINCM
CC
        E=1000.0/EPCM
        X=E
        X2=X**2
        X3=X**3
        X4=X**4
        RECV=(ETCM /SR)
        RECS=(WT/SR)
        WRITE(6,111) E,ETCM,WT,SR,RECV,RECS
 111     FORMAT('  E,ETCM,WT,SR,RECV,RECS',/,1X,6(1PE10.3))
CC-------------------------------------------------------------------
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUMR=-100.* (P(1, 1)+P(1, 2)*X+P(1 ,3)*X2+P(1, 4)*X3+P(1, 5)*X4)
      RV1=        (P(2, 1)+P(2, 2)*X+P(2 ,3)*X2+P(2, 4)*X3+P(2, 5)*X4)
      AV1=0.7*    (P(3, 1)+P(3, 2)*X+P(3 ,3)*X2+P(3, 4)*X3+P(3, 5)*X4)
      SUMI=-15.0* (P(4, 1)+P(4, 2)*X+P(4 ,3)*X2+P(4, 4)*X3+P(4, 5)*X4)
      RV2=        (P(5, 1)+P(5, 2)*X+P(5 ,3)*X2+P(5, 4)*X3+P(5, 5)*X4)
      AV2=0.7*    (P(6, 1)+P(6, 2)*X+P(6 ,3)*X2+P(6, 4)*X3+P(6, 5)*X4)
      DIFFR=700.* (P(7, 1)+P(7, 2)*X+P(7 ,3)*X2+P(7, 4)*X3+P(7, 5)*X4)
      RS1=        (P(8, 1)+P(8, 2)*X+P(8 ,3)*X2+P(8, 4)*X3+P(8, 5)*X4)
      AS1=0.7*    (P(9, 1)+P(9, 2)*X+P(9 ,3)*X2+P(9, 4)*X3+P(9, 5)*X4)
      DIFFI=-150.*(P(10,1)+P(10,2)*X+P(10,3)*X2+P(10,4)*X3+P(10,5)*X4)
      RS2=        (P(11,1)+P(11,2)*X+P(11,3)*X2+P(11,4)*X3+P(11,5)*X4)
      AS2=0.7*    (P(12,1)+P(12,2)*X+P(12,3)*X2+P(12,4)*X3+P(12,5)*X4)
CC-----------------------------------------------
      VV=0.5*(SUMR+DIFFR)
      VS=0.5*(SUMR-DIFFR)
      WV=0.5*(SUMI+DIFFI)
      WS=0.5*(SUMI-DIFFI)
CC-------------------------------------------------------------------
C      VV2 =0.0
C      RV12=  1.0
C      AV12=  0.65
      WV2 = -100.0* (P(13,1)+P(13,2)*X+P(13,3)*X2+P(13,4)*X3+P(13,5)*X4)
      RV22=         (P(14,1)+P(14,2)*X+P(14,3)*X2+P(14,4)*X3+P(14,5)*X4)
      AV22=0.7*     (P(15,1)+P(15,2)*X+P(15,3)*X2+P(15,4)*X3+P(15,5)*X4)
C
C      VS2 =0.0
C      RS12=  1.0
C      AS12=  0.65
      WS2 = 100.0*(P(16,1)+P(16,2)*X+P(16,3)*X2+P(16,4)*X3+P(16,5)*X4)
      RS22=       (P(17,1)+P(17,2)*X+P(17,3)*X2+P(17,4)*X3+P(17,5)*X4)
      AS22=0.7*   (P(18,1)+P(18,2)*X+P(18,3)*X2+P(18,4)*X3+P(18,5)*X4)
CC
CC--- S.P  GEOMETRIES ARE FIXED BY THEIR VOLUME PART.
      RV22=RV2
      AV22=AV2
      RS22=RS2
      AS22=AS2
CC
      DSMALL=0.0
cx      DR=DMESH
      DO 10 I=1,NSTEP
         TMP1 = R(I)
            S1 = DSECH(RV1*ACB/AV1)
            S2 = DSECH(TMP1/AV1)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RVA1=RECV*VV*B/(A+B)
            S1 = DSECH(RV2*ACB/AV2)
            S2 = DSECH(TMP1/AV2)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RVA2=RECV*WV*B/(A+B)
            S1 = DSECH(RV22*ACB/AV22)
            S2 = DSECH(TMP1/AV22)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECV*WV2*A*B/(A+B)**2
         RVA2=RVA2+SURF
            S1 = DSECH(RS1*ACB/AS1)
            S2 = DSECH(TMP1/AS1)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RSA1=RECS*VS*B/(A+B)
            S1 = DSECH(RS2*ACB/AS2)
            S2 = DSECH(TMP1/AS2)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RSA2=RECS*WS*B/(A+B)
            S1 = DSECH(RS22*ACB/AS22)
            S2 = DSECH(TMP1/AS22)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECS*WS2*A*B/(A+B)**2
         RSA2=RSA2+SURF
CC
         VVA(I)=RVA1
         WVA(I)=RVA2
         VSA(I)=RSA1
         WSA(I)=RSA2
 10    CONTINUE
CC
      RETURN
      END
      SUBROUTINE EDAIO(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(NSTEP),VVA(NSTEP),WVA(NSTEP),VSA(NSTEP),WSA(NSTEP)
      DIMENSION P(22,8),PT1(8,8),PT2(8,8),PT3(8,6)
C23456789012345678901234567890123456789012345678901234567890123456789012
CC    COSH.SPV 
CC     VOLUME = (COSH(R/A)-1)/(COSH(R/A)+COSH(X/A)-2)
CC     S.P = (COSH(R/A)-1)*(COSH(X/A)-1)/((COSH(R/A)+COSH(X/A)-2)**2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   O16_K.UNIT12 
       DATA PT1/
     & -1.347740E+01,  3.728751E+01, -1.012196E+00, -5.602503E+01,
     &  3.394493E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.532862E+01, -6.241333E+01,  9.886500E+01, -6.718622E+01,
     &  1.637505E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.918290E+01, -1.020230E+02,  2.062910E+02, -1.799009E+02,
     &  5.713442E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -7.524496E+01,  5.105817E+02, -1.147151E+03,  1.047417E+03,
     & -3.360008E+02,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.743018E+01, -9.406074E+01,  1.967567E+02, -1.790756E+02,
     &  5.973467E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -2.946938E+01,  1.648218E+02, -3.301414E+02,  2.867085E+02,
     & -9.082262E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  8.049073E+00, -6.647033E+01,  1.753108E+02, -1.825030E+02,
     &  6.650062E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.341266E+01, -5.318056E+01,  8.214481E+01, -5.364332E+01,
     &  1.222869E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00/
       DATA PT2/
     &  2.109042E+01, -1.098495E+02,  2.185050E+02, -1.886238E+02,
     &  5.964062E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -2.673338E+01,  1.515865E+02, -3.378537E+02,  3.589903E+02,
     & -1.456723E+02,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  4.039061E+00, -3.158135E+01,  8.939129E+01, -9.802777E+01,
     &  3.695629E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -1.988325E+01,  1.223279E+02, -2.619458E+02,  2.400034E+02,
     & -7.948324E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -7.726141E+01,  4.188432E+02, -8.357413E+02,  7.327201E+02,
     & -2.362788E+02,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  7.000000E-01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -1.327100E+02,  7.479223E+02, -1.548295E+03,  1.395019E+03,
     & -4.595383E+02,  0.000000E+00,  0.000000E+00,  0.000000E+00/
       DATA PT3/ 
     &  1.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  7.000000E-01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00/
CC
      DSECH(X) = 2.*DEXP(-DABS(X))/( 1.0D0 + DEXP(-DABS(X))**2 )
CCC
C
      IF(TPLAB.LT.21.0D0.OR.TPLAB.GT.1040.D0) WRITE(6,601)
 601  FORMAT(' ',' YOU ARE EXTRAPOLATING IN ENERGY !!?? '/)
      IF(DABS(AA-16.0D0).GT.1.D-33) WRITE(6,602)
 602  FORMAT(' ',' YOU ARE EXTRAPOLATING IN MASS NUMBER !!?? '/)
CCC
      DO 20 I=1,8
         DO 30 J=1,8
            P(J,I)=PT1(I,J)
 30      CONTINUE
         DO 40 J=1,8
            P(8+J,I)=PT2(I,J)
 40      CONTINUE
         DO 50 J=1,6
            P(16+J,I)=PT3(I,J)
 50      CONTINUE
 20   CONTINUE
CCC
C      PI=3.141592653589793D0
      AMU=931.5016D0
      WT=AA*AMU
      EE=TPLAB
      ACB=AA**0.33333333333333333D0
CC
CC
CC--------------------------------------------------------------
CC TRANSFORM TO 2 BODY CM FRAME AND CALCULATE RECOIL FACTORS
CC
        WP=1.0072545D0*AMU
C        CWT=AA**0.3333333333333333D0
        EL=EE+WP
        WP2=WP**2
        WT2=WT**2
        PCM=DSQRT(WT2*(EL**2-WP2)/(WP2+WT2+2.0D0*WT*EL))
        EPCM=DSQRT(WP2+PCM**2)
        ETCM=DSQRT(WT2+PCM**2)
        SR=EPCM+ETCM
C        EKINCM=SR-WP-WT
C        TCM=EKINCM
CC
        E=1000.0/EPCM
        X=E
        X2=X**2
        X3=X**3
        X4=X**4
        RECV=(ETCM /SR)
        RECS=(WT/SR)
        WRITE(6,111) E,ETCM,WT,SR,RECV,RECS
 111     FORMAT('  E,ETCM,WT,SR,RECV,RECS',/,1X,6(1PE10.3))
CC-------------------------------------------------------------------
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUMR=-100.* (P(1, 1)+P(1, 2)*X+P(1 ,3)*X2+P(1, 4)*X3+P(1, 5)*X4)
      RV1=        (P(2, 1)+P(2, 2)*X+P(2 ,3)*X2+P(2, 4)*X3+P(2, 5)*X4)
      AV1=0.7*    (P(3, 1)+P(3, 2)*X+P(3 ,3)*X2+P(3, 4)*X3+P(3, 5)*X4)
      SUMI=-15.0* (P(4, 1)+P(4, 2)*X+P(4 ,3)*X2+P(4, 4)*X3+P(4, 5)*X4)
      RV2=        (P(5, 1)+P(5, 2)*X+P(5 ,3)*X2+P(5, 4)*X3+P(5, 5)*X4)
      AV2=0.7*    (P(6, 1)+P(6, 2)*X+P(6 ,3)*X2+P(6, 4)*X3+P(6, 5)*X4)
      DIFFR=700.* (P(7, 1)+P(7, 2)*X+P(7 ,3)*X2+P(7, 4)*X3+P(7, 5)*X4)
      RS1=        (P(8, 1)+P(8, 2)*X+P(8 ,3)*X2+P(8, 4)*X3+P(8, 5)*X4)
      AS1=0.7*    (P(9, 1)+P(9, 2)*X+P(9 ,3)*X2+P(9, 4)*X3+P(9, 5)*X4)
      DIFFI=-150.*(P(10,1)+P(10,2)*X+P(10,3)*X2+P(10,4)*X3+P(10,5)*X4)
      RS2=        (P(11,1)+P(11,2)*X+P(11,3)*X2+P(11,4)*X3+P(11,5)*X4)
      AS2=0.7*    (P(12,1)+P(12,2)*X+P(12,3)*X2+P(12,4)*X3+P(12,5)*X4)
CC-----------------------------------------------
      VV=0.5*(SUMR+DIFFR)
      VS=0.5*(SUMR-DIFFR)
      WV=0.5*(SUMI+DIFFI)
      WS=0.5*(SUMI-DIFFI)
CC-------------------------------------------------------------------
C      VV2 =0.0
C      RV12=  1.0
C      AV12=  0.65
      WV2 = -100.0* (P(13,1)+P(13,2)*X+P(13,3)*X2+P(13,4)*X3+P(13,5)*X4)
      RV22=         (P(14,1)+P(14,2)*X+P(14,3)*X2+P(14,4)*X3+P(14,5)*X4)
      AV22=0.7*     (P(15,1)+P(15,2)*X+P(15,3)*X2+P(15,4)*X3+P(15,5)*X4)
C
C      VS2 =0.0
C      RS12=  1.0
C      AS12=  0.65
      WS2 = 100.0*(P(16,1)+P(16,2)*X+P(16,3)*X2+P(16,4)*X3+P(16,5)*X4)
      RS22=       (P(17,1)+P(17,2)*X+P(17,3)*X2+P(17,4)*X3+P(17,5)*X4)
      AS22=0.7*   (P(18,1)+P(18,2)*X+P(18,3)*X2+P(18,4)*X3+P(18,5)*X4)
CC
CC--- S.P  GEOMETRIES ARE FIXED BY THEIR VOLUME PART.
      RV22=RV2
      AV22=AV2
      RS22=RS2
      AS22=AS2
CC
      DSMALL=0.0
cx      DR=DMESH
      DO 10 I=1,NSTEP
         TMP1 = R(I)
            S1 = DSECH(RV1*ACB/AV1)
            S2 = DSECH(TMP1/AV1)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RVA1=RECV*VV*B/(A+B)
            S1 = DSECH(RV2*ACB/AV2)
            S2 = DSECH(TMP1/AV2)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RVA2=RECV*WV*B/(A+B)
            S1 = DSECH(RV22*ACB/AV22)
            S2 = DSECH(TMP1/AV22)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECV*WV2*A*B/(A+B)**2
         RVA2=RVA2+SURF
            S1 = DSECH(RS1*ACB/AS1)
            S2 = DSECH(TMP1/AS1)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RSA1=RECS*VS*B/(A+B)
            S1 = DSECH(RS2*ACB/AS2)
            S2 = DSECH(TMP1/AS2)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RSA2=RECS*WS*B/(A+B)
            S1 = DSECH(RS22*ACB/AS22)
            S2 = DSECH(TMP1/AS22)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECS*WS2*A*B/(A+B)**2
         RSA2=RSA2+SURF
CC
         VVA(I)=RVA1
         WVA(I)=RVA2
         VSA(I)=RSA1
         WSA(I)=RSA2
 10    CONTINUE
CC
      RETURN
      END
      SUBROUTINE EDAICA(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(NSTEP),VVA(NSTEP),WVA(NSTEP),VSA(NSTEP),WSA(NSTEP)
      DIMENSION P(22,8),PT1(8,8),PT2(8,8),PT3(8,6)
C23456789012345678901234567890123456789012345678901234567890123456789012
CC    COSH.SPV 
CC     VOLUME = (COSH(R/A)-1)/(COSH(R/A)+COSH(X/A)-2)
CC     S.P = (COSH(R/A)-1)*(COSH(X/A)-1)/((COSH(R/A)+COSH(X/A)-2)**2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   CA40_K.UNIT12 
       DATA PT1/
     & -7.318954E+00,  5.398397E+01, -1.167047E+02,  1.046894E+02,
     & -3.370271E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  3.630568E+00, -1.110169E+01,  1.726085E+01, -1.151544E+01,
     &  2.774694E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  2.133154E+01, -1.082708E+02,  2.135381E+02, -1.846674E+02,
     &  5.886944E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  2.948211E+01, -1.374328E+02,  2.374793E+02, -1.799036E+02,
     &  5.057989E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  9.605070E+00, -4.859204E+01,  1.012495E+02, -9.188551E+01,
     &  3.083748E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -1.207646E+01,  7.067181E+01, -1.437512E+02,  1.281344E+02,
     & -4.211021E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.580615E+00, -9.748488E+00,  2.581606E+01, -2.597218E+01,
     &  9.420992E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  2.219730E+00, -5.548538E+00,  9.643271E+00, -7.363312E+00,
     &  2.087623E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00/
       DATA PT2/
     &  1.988024E+01, -1.010234E+02,  2.003805E+02, -1.741487E+02,
     &  5.576321E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -6.443876E+01,  3.708830E+02, -7.661438E+02,  6.922514E+02,
     & -2.321457E+02,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.433837E+01, -7.087190E+01,  1.406530E+02, -1.224460E+02,
     &  3.957107E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -1.917052E+00,  8.549128E+00, -2.013102E+01,  1.824768E+01,
     & -5.599311E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -1.217677E+01,  7.681134E+01, -1.687796E+02,  1.541948E+02,
     & -5.009802E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  7.000000E-01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -3.616350E+01,  1.862461E+02, -3.602911E+02,  3.032989E+02,
     & -9.335186E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00/
       DATA PT3/ 
     &  1.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  7.000000E-01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00/
CC
      DSECH(X) = 2.*DEXP(-DABS(X))/( 1.0D0 + DEXP(-DABS(X))**2 )
CCC
      IF(TPLAB.LT.21.0D0.OR.TPLAB.GT.1040.D0) WRITE(6,601)
 601  FORMAT(' ',' YOU ARE EXTRAPOLATING IN ENERGY !!?? '/)
      IF(DABS(AA-40.0D0).GT.1.D-33) WRITE(6,602)
 602  FORMAT(' ',' YOU ARE EXTRAPOLATING IN MASS NUMBER !!?? '/)
CCC
      DO 20 I=1,8
         DO 30 J=1,8
            P(J,I)=PT1(I,J)
 30      CONTINUE
         DO 40 J=1,8
            P(8+J,I)=PT2(I,J)
 40      CONTINUE
         DO 50 J=1,6
            P(16+J,I)=PT3(I,J)
 50      CONTINUE
 20   CONTINUE
CCC
C      PI=3.141592653589793D0
      AMU=931.5016D0
      WT=AA*AMU
      EE=TPLAB
      ACB=AA**0.33333333333333333D0
CC
CC
CC--------------------------------------------------------------
CC TRANSFORM TO 2 BODY CM FRAME AND CALCULATE RECOIL FACTORS
CC
        WP=1.0072545D0*AMU
C        CWT=AA**0.3333333333333333D0
        EL=EE+WP
        WP2=WP**2
        WT2=WT**2
        PCM=DSQRT(WT2*(EL**2-WP2)/(WP2+WT2+2.0D0*WT*EL))
        EPCM=DSQRT(WP2+PCM**2)
        ETCM=DSQRT(WT2+PCM**2)
        SR=EPCM+ETCM
C        EKINCM=SR-WP-WT
C        TCM=EKINCM
CC
        E=1000.0/EPCM
        X=E
        X2=X**2
        X3=X**3
        X4=X**4
        RECV=(ETCM /SR)
        RECS=(WT/SR)
        WRITE(6,111) E,ETCM,WT,SR,RECV,RECS
 111     FORMAT('  E,ETCM,WT,SR,RECV,RECS',/,1X,6(1PE10.3))
CC-------------------------------------------------------------------
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUMR=-100.* (P(1, 1)+P(1, 2)*X+P(1 ,3)*X2+P(1, 4)*X3+P(1, 5)*X4)
      RV1=        (P(2, 1)+P(2, 2)*X+P(2 ,3)*X2+P(2, 4)*X3+P(2, 5)*X4)
      AV1=0.7*    (P(3, 1)+P(3, 2)*X+P(3 ,3)*X2+P(3, 4)*X3+P(3, 5)*X4)
      SUMI=-15.0* (P(4, 1)+P(4, 2)*X+P(4 ,3)*X2+P(4, 4)*X3+P(4, 5)*X4)
      RV2=        (P(5, 1)+P(5, 2)*X+P(5 ,3)*X2+P(5, 4)*X3+P(5, 5)*X4)
      AV2=0.7*    (P(6, 1)+P(6, 2)*X+P(6 ,3)*X2+P(6, 4)*X3+P(6, 5)*X4)
      DIFFR=700.* (P(7, 1)+P(7, 2)*X+P(7 ,3)*X2+P(7, 4)*X3+P(7, 5)*X4)
      RS1=        (P(8, 1)+P(8, 2)*X+P(8 ,3)*X2+P(8, 4)*X3+P(8, 5)*X4)
      AS1=0.7*    (P(9, 1)+P(9, 2)*X+P(9 ,3)*X2+P(9, 4)*X3+P(9, 5)*X4)
      DIFFI=-150.*(P(10,1)+P(10,2)*X+P(10,3)*X2+P(10,4)*X3+P(10,5)*X4)
      RS2=        (P(11,1)+P(11,2)*X+P(11,3)*X2+P(11,4)*X3+P(11,5)*X4)
      AS2=0.7*    (P(12,1)+P(12,2)*X+P(12,3)*X2+P(12,4)*X3+P(12,5)*X4)
CC-----------------------------------------------
      VV=0.5*(SUMR+DIFFR)
      VS=0.5*(SUMR-DIFFR)
      WV=0.5*(SUMI+DIFFI)
      WS=0.5*(SUMI-DIFFI)
CC-------------------------------------------------------------------
C      VV2 =0.0
C      RV12=  1.0
C      AV12=  0.65
      WV2 = -100.0* (P(13,1)+P(13,2)*X+P(13,3)*X2+P(13,4)*X3+P(13,5)*X4)
      RV22=         (P(14,1)+P(14,2)*X+P(14,3)*X2+P(14,4)*X3+P(14,5)*X4)
      AV22=0.7*     (P(15,1)+P(15,2)*X+P(15,3)*X2+P(15,4)*X3+P(15,5)*X4)
C
C      VS2 =0.0
C      RS12=  1.0
C      AS12=  0.65
      WS2 = 100.0*(P(16,1)+P(16,2)*X+P(16,3)*X2+P(16,4)*X3+P(16,5)*X4)
      RS22=       (P(17,1)+P(17,2)*X+P(17,3)*X2+P(17,4)*X3+P(17,5)*X4)
      AS22=0.7*   (P(18,1)+P(18,2)*X+P(18,3)*X2+P(18,4)*X3+P(18,5)*X4)
CC
CC--- S.P  GEOMETRIES ARE FIXED BY THEIR VOLUME PART.
      RV22=RV2
      AV22=AV2
      RS22=RS2
      AS22=AS2
CC
      DSMALL=0.0
cx      DR=DMESH
      DO 10 I=1,NSTEP
         TMP1 = R(I)
            S1 = DSECH(RV1*ACB/AV1)
            S2 = DSECH(TMP1/AV1)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RVA1=RECV*VV*B/(A+B)
            S1 = DSECH(RV2*ACB/AV2)
            S2 = DSECH(TMP1/AV2)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RVA2=RECV*WV*B/(A+B)
            S1 = DSECH(RV22*ACB/AV22)
            S2 = DSECH(TMP1/AV22)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECV*WV2*A*B/(A+B)**2
         RVA2=RVA2+SURF
            S1 = DSECH(RS1*ACB/AS1)
            S2 = DSECH(TMP1/AS1)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RSA1=RECS*VS*B/(A+B)
            S1 = DSECH(RS2*ACB/AS2)
            S2 = DSECH(TMP1/AS2)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RSA2=RECS*WS*B/(A+B)
            S1 = DSECH(RS22*ACB/AS22)
            S2 = DSECH(TMP1/AS22)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECS*WS2*A*B/(A+B)**2
         RSA2=RSA2+SURF
CC
         VVA(I)=RVA1
         WVA(I)=RVA2
         VSA(I)=RSA1
         WSA(I)=RSA2
 10    CONTINUE
CC
      RETURN
      END
      SUBROUTINE EDAIZR(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(NSTEP),VVA(NSTEP),WVA(NSTEP),VSA(NSTEP),WSA(NSTEP)
      DIMENSION P(22,8),PT1(8,8),PT2(8,8),PT3(8,6)
C23456789012345678901234567890123456789012345678901234567890123456789012
CC    COSH.SPV 
CC     VOLUME = (COSH(R/A)-1)/(COSH(R/A)+COSH(X/A)-2)
CC     S.P = (COSH(R/A)-1)*(COSH(X/A)-1)/((COSH(R/A)+COSH(X/A)-2)**2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   ZR90_K.UNIT12 
       DATA PT1/
     & -4.204949E+00,  1.582555E+01, -2.926990E+00, -2.085166E+01,
     &  1.308523E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  9.622712E-01,  2.567417E-01, -1.097522E-01, -4.247417E-01,
     &  4.181596E-01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  5.709230E-01,  2.071186E+00, -2.771482E+00,  1.164592E+00,
     & -2.527041E-01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -9.990783E+00,  1.018055E+01,  4.023136E+00, -2.936995E+00,
     & -1.358213E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  2.047547E+00, -2.642592E+00,  1.131523E+00,  2.025188E+00,
     & -1.517759E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.635045E+00, -3.403219E+00,  3.347312E+00,  1.588919E+00,
     & -2.347243E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -3.459307E+00,  1.275814E+01, -8.870889E+00, -3.485597E+00,
     &  4.107062E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.240151E+00, -4.515656E-01, -8.056938E-04,  3.614969E-01,
     & -5.332916E-02,  0.000000E+00,  0.000000E+00,  0.000000E+00/
       DATA PT2/
     &  5.431843E-02,  2.239243E+00,  1.609727E+00, -6.068567E+00,
     &  2.998854E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  5.592459E+00, -7.400733E+00,  2.308148E+00,  6.794854E+00,
     & -6.606420E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  2.160236E+00, -2.726284E+00, -4.136885E-01,  5.535261E+00,
     & -3.403082E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.803567E+00, -7.252405E+00,  2.431111E+00,  6.176691E+00,
     & -4.112217E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  2.619201E+00, -1.489489E+01,  1.192544E+01,  1.755134E+01,
     & -1.733801E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  7.000000E-01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -1.959226E+00,  9.617828E-02,  1.897640E+00,  2.330366E+00,
     & -3.184638E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00/
       DATA PT3/ 
     &  1.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  7.000000E-01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00/
CC
      DSECH(X) = 2.*DEXP(-DABS(X))/( 1.0D0 + DEXP(-DABS(X))**2 )
CCC
      IF(TPLAB.LT.21.0D0.OR.TPLAB.GT.1040.D0) WRITE(6,601)
 601  FORMAT(' ',' YOU ARE EXTRAPOLATING IN ENERGY !!?? '/)
      IF(DABS(AA-90.0D0).GT.1.D-33) WRITE(6,602)
 602  FORMAT(' ',' YOU ARE EXTRAPOLATING IN MASS NUMBER !!?? '/)
CCC
      DO 20 I=1,8
         DO 30 J=1,8
            P(J,I)=PT1(I,J)
 30      CONTINUE
         DO 40 J=1,8
            P(8+J,I)=PT2(I,J)
 40      CONTINUE
         DO 50 J=1,6
            P(16+J,I)=PT3(I,J)
 50      CONTINUE
 20   CONTINUE
CCC
C      PI=3.141592653589793D0
      AMU=931.5016D0
      WT=AA*AMU
      EE=TPLAB
      ACB=AA**0.33333333333333333D0
CC
CC
CC--------------------------------------------------------------
CC TRANSFORM TO 2 BODY CM FRAME AND CALCULATE RECOIL FACTORS
CC
        WP=1.0072545D0*AMU
C        CWT=AA**0.3333333333333333D0
        EL=EE+WP
        WP2=WP**2
        WT2=WT**2
        PCM=DSQRT(WT2*(EL**2-WP2)/(WP2+WT2+2.0D0*WT*EL))
        EPCM=DSQRT(WP2+PCM**2)
        ETCM=DSQRT(WT2+PCM**2)
        SR=EPCM+ETCM
C        EKINCM=SR-WP-WT
C        TCM=EKINCM
CC
        E=1000.0/EPCM
        X=E
        X2=X**2
        X3=X**3
        X4=X**4
        RECV=(ETCM /SR)
        RECS=(WT/SR)
        WRITE(6,111) E,ETCM,WT,SR,RECV,RECS
 111     FORMAT('  E,ETCM,WT,SR,RECV,RECS',/,1X,6(1PE10.3))
CC-------------------------------------------------------------------
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUMR=-100.* (P(1, 1)+P(1, 2)*X+P(1 ,3)*X2+P(1, 4)*X3+P(1, 5)*X4)
      RV1=        (P(2, 1)+P(2, 2)*X+P(2 ,3)*X2+P(2, 4)*X3+P(2, 5)*X4)
      AV1=0.7*    (P(3, 1)+P(3, 2)*X+P(3 ,3)*X2+P(3, 4)*X3+P(3, 5)*X4)
      SUMI=-15.0* (P(4, 1)+P(4, 2)*X+P(4 ,3)*X2+P(4, 4)*X3+P(4, 5)*X4)
      RV2=        (P(5, 1)+P(5, 2)*X+P(5 ,3)*X2+P(5, 4)*X3+P(5, 5)*X4)
      AV2=0.7*    (P(6, 1)+P(6, 2)*X+P(6 ,3)*X2+P(6, 4)*X3+P(6, 5)*X4)
      DIFFR=700.* (P(7, 1)+P(7, 2)*X+P(7 ,3)*X2+P(7, 4)*X3+P(7, 5)*X4)
      RS1=        (P(8, 1)+P(8, 2)*X+P(8 ,3)*X2+P(8, 4)*X3+P(8, 5)*X4)
      AS1=0.7*    (P(9, 1)+P(9, 2)*X+P(9 ,3)*X2+P(9, 4)*X3+P(9, 5)*X4)
      DIFFI=-150.*(P(10,1)+P(10,2)*X+P(10,3)*X2+P(10,4)*X3+P(10,5)*X4)
      RS2=        (P(11,1)+P(11,2)*X+P(11,3)*X2+P(11,4)*X3+P(11,5)*X4)
      AS2=0.7*    (P(12,1)+P(12,2)*X+P(12,3)*X2+P(12,4)*X3+P(12,5)*X4)
CC-----------------------------------------------
      VV=0.5*(SUMR+DIFFR)
      VS=0.5*(SUMR-DIFFR)
      WV=0.5*(SUMI+DIFFI)
      WS=0.5*(SUMI-DIFFI)
CC-------------------------------------------------------------------
C      VV2 =0.0
C      RV12=  1.0
C      AV12=  0.65
      WV2 = -100.0* (P(13,1)+P(13,2)*X+P(13,3)*X2+P(13,4)*X3+P(13,5)*X4)
      RV22=         (P(14,1)+P(14,2)*X+P(14,3)*X2+P(14,4)*X3+P(14,5)*X4)
      AV22=0.7*     (P(15,1)+P(15,2)*X+P(15,3)*X2+P(15,4)*X3+P(15,5)*X4)
C
C      VS2 =0.0
C      RS12=  1.0
C      AS12=  0.65
      WS2 = 100.0*(P(16,1)+P(16,2)*X+P(16,3)*X2+P(16,4)*X3+P(16,5)*X4)
      RS22=       (P(17,1)+P(17,2)*X+P(17,3)*X2+P(17,4)*X3+P(17,5)*X4)
      AS22=0.7*   (P(18,1)+P(18,2)*X+P(18,3)*X2+P(18,4)*X3+P(18,5)*X4)
CC
CC--- S.P  GEOMETRIES ARE FIXED BY THEIR VOLUME PART.
      RV22=RV2
      AV22=AV2
      RS22=RS2
      AS22=AS2
CC
      DSMALL=0.0
cx      DR=DMESH
      DO 10 I=1,NSTEP
         TMP1 = R(I)
            S1 = DSECH(RV1*ACB/AV1)
            S2 = DSECH(TMP1/AV1)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RVA1=RECV*VV*B/(A+B)
            S1 = DSECH(RV2*ACB/AV2)
            S2 = DSECH(TMP1/AV2)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RVA2=RECV*WV*B/(A+B)
            S1 = DSECH(RV22*ACB/AV22)
            S2 = DSECH(TMP1/AV22)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECV*WV2*A*B/(A+B)**2
         RVA2=RVA2+SURF
            S1 = DSECH(RS1*ACB/AS1)
            S2 = DSECH(TMP1/AS1)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RSA1=RECS*VS*B/(A+B)
            S1 = DSECH(RS2*ACB/AS2)
            S2 = DSECH(TMP1/AS2)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RSA2=RECS*WS*B/(A+B)
            S1 = DSECH(RS22*ACB/AS22)
            S2 = DSECH(TMP1/AS22)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECS*WS2*A*B/(A+B)**2
         RSA2=RSA2+SURF
CC
         VVA(I)=RVA1
         WVA(I)=RVA2
         VSA(I)=RSA1
         WSA(I)=RSA2
 10    CONTINUE
CC
      RETURN
      END
      SUBROUTINE EDAIPB(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(NSTEP),VVA(NSTEP),WVA(NSTEP),VSA(NSTEP),WSA(NSTEP)
      DIMENSION P(22,8),PT1(8,8),PT2(8,8),PT3(8,6)
C23456789012345678901234567890123456789012345678901234567890123456789012
CC    COSH.SPV 
CC     VOLUME = (COSH(R/A)-1)/(COSH(R/A)+COSH(X/A)-2)
CC     S.P = (COSH(R/A)-1)*(COSH(X/A)-1)/((COSH(R/A)+COSH(X/A)-2)**2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   PB208_K.UNIT12 
       DATA PT1/
     &  1.925516E+01, -8.188806E+01,  1.382767E+02, -1.042260E+02,
     &  2.953352E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.052039E+00, -1.034336E+00,  4.394117E+00, -5.739941E+00,
     &  2.453279E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  7.536453E+00, -2.986209E+01,  5.012585E+01, -3.655076E+01,
     &  9.623429E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -5.038134E+00,  6.653209E+01, -1.866829E+02,  1.991322E+02,
     & -7.387818E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -2.530138E+00,  1.507098E+01, -2.308380E+01,  1.646445E+01,
     & -4.809330E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.679953E+01, -7.868251E+01,  1.456568E+02, -1.202732E+02,
     &  3.753159E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  4.737970E+00, -2.099953E+01,  3.959807E+01, -3.233373E+01,
     &  1.006902E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  6.620497E-01,  2.943505E-01,  3.125690E+00, -5.683778E+00,
     &  2.722779E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00/
       DATA PT2/
     &  1.207582E+01, -5.106658E+01,  8.751429E+01, -6.575640E+01,
     &  1.814573E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.620897E+01, -6.981451E+01,  1.120843E+02, -7.170996E+01,
     &  1.381926E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -2.109859E+00,  1.371687E+01, -2.040366E+01,  1.281788E+01,
     & -2.793603E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     & -1.146930E+01,  8.554716E+01, -2.027845E+02,  2.013710E+02,
     & -7.183535E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.850639E+00,  1.464488E+01, -6.528241E+01,  7.398588E+01,
     & -2.479158E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  3.271647E+00, -8.370088E+00, -1.688656E+00,  1.260083E+01,
     & -6.145986E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00/
       DATA PT3/ 
     &  1.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  7.000000E-01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00/
CC
      DSECH(X) = 2.*DEXP(-DABS(X))/( 1.0D0 + DEXP(-DABS(X))**2 )
CCC
      IF(TPLAB.LT.21.0D0.OR.TPLAB.GT.1040.D0) WRITE(6,601)
 601  FORMAT(' ',' YOU ARE EXTRAPOLATING IN ENERGY !!?? '/)
      IF(DABS(AA-208.0D0).GT.1.0D-33) WRITE(6,602)
 602  FORMAT(' ',' YOU ARE EXTRAPOLATING IN MASS NUMBER !!?? '/)
CCC
      DO 20 I=1,8
         DO 30 J=1,8
            P(J,I)=PT1(I,J)
 30      CONTINUE
         DO 40 J=1,8
            P(8+J,I)=PT2(I,J)
 40      CONTINUE
         DO 50 J=1,6
            P(16+J,I)=PT3(I,J)
 50      CONTINUE
 20   CONTINUE
CCC
C      PI=3.141592653589793D0
      AMU=931.5016D0
      WT=AA*AMU
      EE=TPLAB
      ACB=AA**0.33333333333333333D0
CC
CC
CC--------------------------------------------------------------
CC TRANSFORM TO 2 BODY CM FRAME AND CALCULATE RECOIL FACTORS
CC
        WP=1.0072545D0*AMU
C        CWT=AA**0.3333333333333333D0
        EL=EE+WP
        WP2=WP**2
        WT2=WT**2
        PCM=DSQRT(WT2*(EL**2-WP2)/(WP2+WT2+2.0D0*WT*EL))
        EPCM=DSQRT(WP2+PCM**2)
        ETCM=DSQRT(WT2+PCM**2)
        SR=EPCM+ETCM
C        EKINCM=SR-WP-WT
C        TCM=EKINCM
CC
        E=1000.0/EPCM
        X=E
        X2=X**2
        X3=X**3
        X4=X**4
        RECV=(ETCM /SR)
        RECS=(WT/SR)
        WRITE(6,111) E,ETCM,WT,SR,RECV,RECS
 111     FORMAT('  E,ETCM,WT,SR,RECV,RECS',/,1X,6(1PE10.3))
CC-------------------------------------------------------------------
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUMR=-100.* (P(1, 1)+P(1, 2)*X+P(1 ,3)*X2+P(1, 4)*X3+P(1, 5)*X4)
      RV1=        (P(2, 1)+P(2, 2)*X+P(2 ,3)*X2+P(2, 4)*X3+P(2, 5)*X4)
      AV1=0.7*    (P(3, 1)+P(3, 2)*X+P(3 ,3)*X2+P(3, 4)*X3+P(3, 5)*X4)
      SUMI=-15.0* (P(4, 1)+P(4, 2)*X+P(4 ,3)*X2+P(4, 4)*X3+P(4, 5)*X4)
      RV2=        (P(5, 1)+P(5, 2)*X+P(5 ,3)*X2+P(5, 4)*X3+P(5, 5)*X4)
      AV2=0.7*    (P(6, 1)+P(6, 2)*X+P(6 ,3)*X2+P(6, 4)*X3+P(6, 5)*X4)
      DIFFR=700.* (P(7, 1)+P(7, 2)*X+P(7 ,3)*X2+P(7, 4)*X3+P(7, 5)*X4)
      RS1=        (P(8, 1)+P(8, 2)*X+P(8 ,3)*X2+P(8, 4)*X3+P(8, 5)*X4)
      AS1=0.7*    (P(9, 1)+P(9, 2)*X+P(9 ,3)*X2+P(9, 4)*X3+P(9, 5)*X4)
      DIFFI=-150.*(P(10,1)+P(10,2)*X+P(10,3)*X2+P(10,4)*X3+P(10,5)*X4)
      RS2=        (P(11,1)+P(11,2)*X+P(11,3)*X2+P(11,4)*X3+P(11,5)*X4)
      AS2=0.7*    (P(12,1)+P(12,2)*X+P(12,3)*X2+P(12,4)*X3+P(12,5)*X4)
CC-----------------------------------------------
      VV=0.5*(SUMR+DIFFR)
      VS=0.5*(SUMR-DIFFR)
      WV=0.5*(SUMI+DIFFI)
      WS=0.5*(SUMI-DIFFI)
CC-------------------------------------------------------------------
C      VV2 =0.0
C      RV12=  1.0
C      AV12=  0.65
      WV2 = -100.0* (P(13,1)+P(13,2)*X+P(13,3)*X2+P(13,4)*X3+P(13,5)*X4)
      RV22=         (P(14,1)+P(14,2)*X+P(14,3)*X2+P(14,4)*X3+P(14,5)*X4)
      AV22=0.7*     (P(15,1)+P(15,2)*X+P(15,3)*X2+P(15,4)*X3+P(15,5)*X4)
C
C      VS2 =0.0
C      RS12=  1.0
C      AS12=  0.65
      WS2 = 100.0*(P(16,1)+P(16,2)*X+P(16,3)*X2+P(16,4)*X3+P(16,5)*X4)
      RS22=       (P(17,1)+P(17,2)*X+P(17,3)*X2+P(17,4)*X3+P(17,5)*X4)
      AS22=0.7*   (P(18,1)+P(18,2)*X+P(18,3)*X2+P(18,4)*X3+P(18,5)*X4)
CC
CC--- S.P  GEOMETRIES ARE FIXED BY THEIR VOLUME PART.
      RV22=RV2
      AV22=AV2
      RS22=RS2
      AS22=AS2
CC
      DSMALL=0.0
cx      DR=DMESH
      DO 10 I=1,NSTEP
         TMP1 = R(I)
            S1 = DSECH(RV1*ACB/AV1)
            S2 = DSECH(TMP1/AV1)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RVA1=RECV*VV*B/(A+B)
            S1 = DSECH(RV2*ACB/AV2)
            S2 = DSECH(TMP1/AV2)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RVA2=RECV*WV*B/(A+B)
            S1 = DSECH(RV22*ACB/AV22)
            S2 = DSECH(TMP1/AV22)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECV*WV2*A*B/(A+B)**2
         RVA2=RVA2+SURF
            S1 = DSECH(RS1*ACB/AS1)
            S2 = DSECH(TMP1/AS1)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RSA1=RECS*VS*B/(A+B)
            S1 = DSECH(RS2*ACB/AS2)
            S2 = DSECH(TMP1/AS2)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RSA2=RECS*WS*B/(A+B)
            S1 = DSECH(RS22*ACB/AS22)
            S2 = DSECH(TMP1/AS22)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECS*WS2*A*B/(A+B)**2
         RSA2=RSA2+SURF
CC
         VVA(I)=RVA1
         WVA(I)=RVA2
         VSA(I)=RSA1
         WSA(I)=RSA2
 10    CONTINUE
CC
      RETURN
      END
      SUBROUTINE EDAD1(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(NSTEP),VVA(NSTEP),WVA(NSTEP),VSA(NSTEP),WSA(NSTEP)
      DIMENSION P(22,8),PT1(8,8),PT2(8,8),PT3(8,6)
C23456789012345678901234567890123456789012345678901234567890123456789012
CC    A_DEP_FIT4.SPV ------ COSH.SPV  
CC     VOLUME = (COSH(R/A)-1)/(COSH(R/A)+COSH(X/A)-2)
CC     S.P = (COSH(R/A)-1)*(COSH(X/A)-1)/((COSH(R/A)+COSH(X/A)-2)**2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     X = E-DEP (EPS)   Y = A-DEP (CAA)
C--------------------------------------------------------------------
CC   AADEP_BIG.UNIT12
       DATA PT1/
     &  2.313828E+01, -1.233380E+02,  2.432987E+02, -2.092616E+02,
     &  6.606549E+01,  5.645371E+00, -1.014396E+01,  5.733192E+00,
     & -6.841731E+00,  3.747771E+01, -6.858889E+01,  5.484194E+01,
     & -1.617735E+01,  9.973483E-01, -1.044247E+00,  4.707822E-01,
     &  1.688843E+01, -8.404292E+01,  1.550125E+02, -1.251468E+02,
     &  3.722361E+01,  2.912163E+00, -3.560279E+00,  1.645079E+00,
     &  2.143528E+02, -9.614898E+02,  1.576657E+03, -1.147987E+03,
     &  3.136021E+02,  2.226300E+01, -3.363301E+01,  1.704743E+01,
     &  7.852790E+00, -4.551995E+01,  9.525029E+01, -8.644521E+01,
     &  2.891801E+01,  5.105014E+00, -7.512490E+00,  3.681682E+00,
     & -8.332038E+00,  4.210331E+01, -7.878047E+01,  6.498050E+01,
     & -1.997590E+01,  4.603199E+00, -7.773754E+00,  4.281708E+00,
     &  1.771733E+01, -9.514249E+01,  1.826431E+02, -1.507498E+02,
     &  4.563918E+01,  5.566540E+00, -1.001441E+01,  5.551678E+00,
     & -6.557412E+00,  3.596972E+01, -6.571700E+01,  5.233473E+01,
     & -1.534723E+01,  1.007355E+00, -9.669271E-01,  4.066760E-01/
       DATA PT2/
     &  2.183389E+01, -1.085667E+02,  2.012175E+02, -1.637375E+02,
     &  4.926579E+01,  2.874114E+00, -3.559180E+00,  1.647730E+00,
     & -5.944315E+01,  3.517880E+02, -6.896441E+02,  6.080210E+02,
     & -2.026092E+02, -3.492561E+01,  5.197784E+01, -2.516421E+01,
     &  6.414280E-01, -9.685444E+00,  2.992500E+01, -3.446931E+01,
     &  1.361126E+01,  5.209321E+00, -7.563083E+00,  3.671130E+00,
     & -7.989428E+00,  4.318908E+01, -8.862678E+01,  7.860778E+01,
     & -2.582041E+01,  7.390220E+00, -1.207679E+01,  6.360783E+00,
     &  1.871682E+01, -1.144552E+02,  2.566304E+02, -2.509193E+02,
     &  9.003415E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  3.713738E+01, -1.952369E+02,  3.842865E+02, -3.379853E+02,
     &  1.116151E+02,  0.000000E+00,  0.000000E+00,  0.000000E+00/
       DATA PT3/ 
     &  1.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  1.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00/
CC
      DSECH(X) = 2.*DEXP(-DABS(X))/( 1.0D0 + DEXP(-DABS(X))**2 )
C
      IF(TPLAB.LT.21.0D0.OR.TPLAB.GT.1040.D0) WRITE(6,601)
 601  FORMAT(' ',' YOU ARE EXTRAPOLATING IN ENERGY !!?? '/)
      IF(AA.LT.12.D0.OR.AA.GT.208.D0) WRITE(6,602)
 602  FORMAT(' ',' YOU ARE EXTRAPOLATING IN MASS NUMBER !!?? '/)
C
CCC
      DO 20 I=1,8
         DO 30 J=1,8
            P(J,I)=PT1(I,J)
 30      CONTINUE
         DO 40 J=1,8
            P(8+J,I)=PT2(I,J)
 40      CONTINUE
         DO 50 J=1,6
            P(16+J,I)=PT3(I,J)
 50      CONTINUE
 20   CONTINUE
CCC
C      PI=3.141592653589793D0
      AMU=931.5016D0
      WT=AA*AMU
      EE=TPLAB
      ACB=AA**0.33333333333333333D0
CC
      CAA=AA/(AA+20.)
      Y =CAA
      Y2=Y**2
      Y3=Y**3
      Y4=Y**4
CC
CC--------------------------------------------------------------
CC TRANSFORM TO 2 BODY CM FRAME AND CALCULATE RECOIL FACTORS
CC
        WP=1.0072545D0*AMU
C        CWT=AA**0.3333333333333333D0
        EL=EE+WP
        WP2=WP**2
        WT2=WT**2
        PCM=DSQRT(WT2*(EL**2-WP2)/(WP2+WT2+2.0D0*WT*EL))
        EPCM=DSQRT(WP2+PCM**2)
        ETCM=DSQRT(WT2+PCM**2)
        SR=EPCM+ETCM
C        EKINCM=SR-WP-WT
C        TCM=EKINCM
CC
        E=1000.0/EPCM
        X=E
        X2=X**2
        X3=X**3
        X4=X**4
        RECV=(ETCM /SR)
        RECS=(WT/SR)
        WRITE(6,111) E,ETCM,WT,SR,RECV,RECS
 111     FORMAT('  E,ETCM,WT,SR,RECV,RECS',/,1X,6(1PE10.3))
CC-------------------------------------------------------------------
C23456789012345678901234567890123456789012345678901234567890123456789012
      CV1B=1.0
      AV1B=0.7
      CV2B=CV1B
      CS1B=CV1B
      CS2B=CV1B
      AV2B=AV1B
      AS1B=AV1B
      AS2B=AV1B
CC-------------------------------------------------------------------
      SUMR=-100.* (P(1, 1)+P(1, 2)*X+P(1, 3)*X2+P(1, 4)*X3+P(1, 5)*X4
     &                    +P(1, 6)*Y+P(1, 7)*Y2+P(1, 8)*Y3+P(20,1)*Y4)
      RV1=  CV1B* (P(2, 1)+P(2, 2)*X+P(2, 3)*X2+P(2, 4)*X3+P(2, 5)*X4
     &                    +P(2, 6)*Y+P(2, 7)*Y2+P(2, 8)*Y3+P(21,1)*Y4)
      AV1=  AV1B* (P(3, 1)+P(3, 2)*X+P(3, 3)*X2+P(3, 4)*X3+P(3, 5)*X4
     &                    +P(3, 6)*Y+P(3, 7)*Y2+P(3, 8)*Y3)
      SUMI=-15.0* (P(4, 1)+P(4, 2)*X+P(4, 3)*X2+P(4, 4)*X3+P(4, 5)*X4
     &                    +P(4, 6)*Y+P(4, 7)*Y2+P(4, 8)*Y3+P(20,2)*Y4)
      RV2=  CV2B* (P(5, 1)+P(5, 2)*X+P(5, 3)*X2+P(5, 4)*X3+P(5, 5)*X4
     &                    +P(5, 6)*Y+P(5, 7)*Y2+P(5, 8)*Y3+P(21,2)*Y4)
      AV2=  AV2B* (P(6, 1)+P(6, 2)*X+P(6, 3)*X2+P(6, 4)*X3+P(6, 5)*X4
     &                    +P(6, 6)*Y+P(6, 7)*Y2+P(6, 8)*Y3)
      DIFFR=700.* (P(7, 1)+P(7, 2)*X+P(7, 3)*X2+P(7, 4)*X3+P(7, 5)*X4
     &                    +P(7, 6)*Y+P(7, 7)*Y2+P(7, 8)*Y3+P(20,3)*Y4)
      RS1=  CS1B* (P(8, 1)+P(8, 2)*X+P(8, 3)*X2+P(8, 4)*X3+P(8, 5)*X4
     &                    +P(8, 6)*Y+P(8, 7)*Y2+P(8, 8)*Y3+P(21,3)*Y4)
      AS1=  AS1B* (P(9, 1)+P(9, 2)*X+P(9, 3)*X2+P(9, 4)*X3+P(9, 5)*X4
     &                    +P(9, 6)*Y+P(9, 7)*Y2+P(9, 8)*Y3)
      DIFFI=-150.*(P(10,1)+P(10,2)*X+P(10,3)*X2+P(10,4)*X3+P(10,5)*X4
     &                    +P(10,6)*Y+P(10,7)*Y2+P(10,8)*Y3+P(20,4)*Y4)
      RS2=  CS2B* (P(11,1)+P(11,2)*X+P(11,3)*X2+P(11,4)*X3+P(11,5)*X4
     &                    +P(11,6)*Y+P(11,7)*Y2+P(11,8)*Y3+P(21,4)*Y4)
      AS2=  AS2B* (P(12,1)+P(12,2)*X+P(12,3)*X2+P(12,4)*X3+P(12,5)*X4
     &                    +P(12,6)*Y+P(12,7)*Y2+P(12,8)*Y3)
CC-----------------------------------------------
      VV=0.5*(SUMR+DIFFR)
      VS=0.5*(SUMR-DIFFR)
      WV=0.5*(SUMI+DIFFI)
      WS=0.5*(SUMI-DIFFI)
CC-------------------------------------------------------------------
C      VV2 =0.0
C      RV12=  1.0
C      AV12=  0.65
      WV2 = -100.0*(P(13,1)+P(13,2)*X+P(13,3)*X2+P(13,4)*X3+P(13,5)*X4
     &                     +P(13,6)*Y+P(13,7)*Y2+P(13,8)*Y3+P(20,5)*Y4)
      RV22=         (P(14,1)+P(14,2)*X+P(14,3)*X2+P(14,4)*X3+P(14,5)*X4)
      AV22=0.7*     (P(15,1)+P(15,2)*X+P(15,3)*X2+P(15,4)*X3+P(15,5)*X4)
C
C      VS2 =0.0
C      RS12=  1.0
C      AS12=  0.65
      WS2 = 100.0*(P(16,1)+P(16,2)*X+P(16,3)*X2+P(16,4)*X3+P(16,5)*X4
     &                    +P(16,6)*Y+P(16,7)*Y2+P(16,8)*Y3+P(20,6)*Y4)
      RS22=       (P(17,1)+P(17,2)*X+P(17,3)*X2+P(17,4)*X3+P(17,5)*X4)
      AS22=0.7*   (P(18,1)+P(18,2)*X+P(18,3)*X2+P(18,4)*X3+P(18,5)*X4)
CC
CC--- S.P  GEOMETRIES ARE FIXED BY THEIR VOLUME PART.
      RV22=RV2
      AV22=AV2
      RS22=RS2
      AS22=AS2
CC
      DSMALL=0.0
cx      DR=DMESH
      DO 10 I=1,NSTEP
         TMP1 = R(I)
            S1 = DSECH(RV1*ACB/AV1)
            S2 = DSECH(TMP1/AV1)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RVA1=RECV*VV*B/(A+B)
            S1 = DSECH(RV2*ACB/AV2)
            S2 = DSECH(TMP1/AV2)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RVA2=RECV*WV*B/(A+B)
            S1 = DSECH(RV22*ACB/AV22)
            S2 = DSECH(TMP1/AV22)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECV*WV2*A*B/(A+B)**2
         RVA2=RVA2+SURF
            S1 = DSECH(RS1*ACB/AS1)
            S2 = DSECH(TMP1/AS1)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RSA1=RECS*VS*B/(A+B)
            S1 = DSECH(RS2*ACB/AS2)
            S2 = DSECH(TMP1/AS2)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RSA2=RECS*WS*B/(A+B)
            S1 = DSECH(RS22*ACB/AS22)
            S2 = DSECH(TMP1/AS22)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECS*WS2*A*B/(A+B)**2
         RSA2=RSA2+SURF
CC
         VVA(I)=RVA1
         WVA(I)=RVA2
         VSA(I)=RSA1
         WSA(I)=RSA2
 10    CONTINUE
CC
      RETURN
      END
      SUBROUTINE EDAD2(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(NSTEP),VVA(NSTEP),WVA(NSTEP),VSA(NSTEP),WSA(NSTEP)
      DIMENSION P(22,8),PT1(8,8),PT2(8,8),PT3(8,6)
C23456789012345678901234567890123456789012345678901234567890123456789012
CC    A_DEP_FIT4D4.SPV ------ COSH.SPV  
CC     VOLUME = (COSH(R/A)-1)/(COSH(R/A)+COSH(X/A)-2)
CC     S.P = (COSH(R/A)-1)*(COSH(X/A)-1)/((COSH(R/A)+COSH(X/A)-2)**2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     X = E-DEP (EPS)   Y = A-DEP (CAA)
C--------------------------------------------------------------------
CCC   AADEP_BIG4.UNIT12 
       DATA PT1/
     &  2.419272E+01, -1.347303E+02,  2.565835E+02, -2.130657E+02,
     &  6.640990E+01,  1.370940E+01, -2.041843E+01,  7.075551E+00,
     & -7.996201E+00,  4.447337E+01, -8.209251E+01,  6.581308E+01,
     & -1.958522E+01,  4.233910E-01, -4.699721E-01,  4.984338E-01,
     &  1.175859E+01, -5.709638E+01,  1.018926E+02, -8.029220E+01,
     &  2.348163E+01,  2.467907E+00, -4.579546E+00,  2.374239E+00,
     &  2.983074E+02, -1.504920E+03,  2.599143E+03, -1.941105E+03,
     &  5.424824E+02,  8.374511E+01, -9.919251E+01,  1.445072E+01,
     &  6.406925E+00, -3.597956E+01,  6.812012E+01, -5.464849E+01,
     &  1.708801E+01,  4.927107E+00, -6.184351E+00,  1.744822E+00,
     & -1.390377E+01,  6.935514E+01, -1.199980E+02,  9.329035E+01,
     & -2.795019E+01,  5.004301E+00, -1.251040E+00,  1.823093E+00,
     &  2.049078E+01, -1.113498E+02,  2.078705E+02, -1.672745E+02,
     &  4.995035E+01,  8.214217E+00, -1.343718E+01,  5.939270E+00,
     & -7.614336E+00,  4.246558E+01, -7.769977E+01,  6.164974E+01,
     & -1.820658E+01,  1.267682E-01,  3.590679E-02,  3.889028E-01/
       DATA PT2/
     &  1.594873E+01, -7.751876E+01,  1.397383E+02, -1.114398E+02,
     &  3.316194E+01,  2.398806E+00, -4.620640E+00,  2.264939E+00,
     & -6.784464E+01,  4.019749E+02, -7.220732E+02,  5.760043E+02,
     & -1.782302E+02, -5.040930E+01,  7.959955E+01, -2.771373E+01,
     & -4.752988E+00,  8.257322E+00, -4.803583E+00, -6.860379E-01,
     &  2.084371E+00,  1.283227E+01, -1.491639E+01,  4.350842E+00,
     & -1.485421E+01,  8.094352E+01, -1.543248E+02,  1.272100E+02,
     & -3.885699E+01,  7.132415E+00, -1.182152E+01,  6.334432E+00,
     &  2.144962E+01, -1.353171E+02,  2.755466E+02, -2.605909E+02,
     &  9.401126E+01,  1.111419E+01, -3.690369E+01,  1.774420E+01,
     & -2.062301E-02,  9.013484E-01, -7.911715E-01,  3.133864E+00,
     & -1.699538E+00, -4.994396E-01,  0.000000E+00,  0.000000E+00,
     &  1.658819E+00, -4.383188E+00,  2.324352E+00, -1.355365E+01,
     &  9.547320E+00, -1.238569E+00,  0.000000E+00,  0.000000E+00,
     &  4.877378E+01, -2.516897E+02,  4.693758E+02, -3.746622E+02,
     &  1.099748E+02, -7.058762E-01,  1.635300E+01, -1.092225E+01/
       DATA PT3/ 
     & -8.331280E-02,  1.215059E+00, -1.131819E+00,  3.421545E+00,
     & -2.214773E+00, -1.626730E-01,  0.000000E+00,  0.000000E+00,
     & -3.889444E+00, -3.422776E+00,  6.241782E+00, -4.801031E+00,
     &  2.252025E+00,  1.982958E+00,  0.000000E+00,  0.000000E+00,
     &  1.221290E+00, -7.093230E+00,  7.266262E+00, -2.368339E+00,
     & -5.509552E+01,  6.461826E+01,  0.000000E+00,  0.000000E+00,
     &  2.000496E+00, -3.189503E+00,  2.003387E+00, -2.493790E+01,
     &  3.281269E+01, -1.877461E+01,  0.000000E+00,  0.000000E+00,
     &  3.119507E+01, -1.996679E+01,  2.433381E+00, -1.527184E+01,
     &  5.443329E+00,  2.462658E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00/
CC
      DSECH(X) = 2.*DEXP(-DABS(X))/( 1.0D0 + DEXP(-DABS(X))**2 )
C
      IF(TPLAB.LT.21.0D0.OR.TPLAB.GT.1040.D0) WRITE(6,601)
 601  FORMAT(' ',' YOU ARE EXTRAPOLATING IN ENERGY !!?? '/)
      IF(AA.LT.12.D0.OR.AA.GT.208.D0) WRITE(6,602)
 602  FORMAT(' ',' YOU ARE EXTRAPOLATING IN MASS NUMBER !!?? '/)
C
CCC
      DO 20 I=1,8
         DO 30 J=1,8
            P(J,I)=PT1(I,J)
 30      CONTINUE
         DO 40 J=1,8
            P(8+J,I)=PT2(I,J)
 40      CONTINUE
         DO 50 J=1,6
            P(16+J,I)=PT3(I,J)
 50      CONTINUE
 20   CONTINUE
CCC
C      PI=3.141592653589793D0
      AMU=931.5016D0
      WT=AA*AMU
      EE=TPLAB
      ACB=AA**0.33333333333333333D0
CC
      CAA=AA/(AA+20.)
      Y =CAA
      Y2=Y**2
      Y3=Y**3
C      Y4=Y**4
CC
CC--------------------------------------------------------------
CC TRANSFORM TO 2 BODY CM FRAME AND CALCULATE RECOIL FACTORS
CC
        WP=1.0072545D0*AMU
C        CWT=AA**0.3333333333333333D0
        EL=EE+WP
        WP2=WP**2
        WT2=WT**2
        PCM=DSQRT(WT2*(EL**2-WP2)/(WP2+WT2+2.0D0*WT*EL))
        EPCM=DSQRT(WP2+PCM**2)
        ETCM=DSQRT(WT2+PCM**2)
        SR=EPCM+ETCM
C        EKINCM=SR-WP-WT
C        TCM=EKINCM
CC
        E=1000.0/EPCM
        X=E
        X2=X**2
        X3=X**3
        X4=X**4
        XY=X*Y
        X2Y=X*X*Y
        XY2=X*Y*Y
        RECV=(ETCM /SR)
        RECS=(WT/SR)
        WRITE(6,111) E,ETCM,WT,SR,RECV,RECS
 111     FORMAT('  E,ETCM,WT,SR,RECV,RECS',/,1X,6(1PE10.3))
CC-------------------------------------------------------------------
C23456789012345678901234567890123456789012345678901234567890123456789012
      CV1B=1.0
      AV1B=0.7
      CV2B=CV1B
      CS1B=CV1B
      CS2B=CV1B
      AV2B=AV1B
      AS1B=AV1B
      AS2B=AV1B
CC-------------------------------------------------------------------
      SUMR=-100.* (P(1, 1)+P(1, 2)*X+P(1, 3)*X2+P(1, 4)*X3+P(1, 5)*X4
     &                    +P(1, 6)*Y+P(1, 7)*Y2+P(1, 8)*Y3
     &                    +P(19,1)*XY+P(19,2)*X2Y+P(19,3)*XY2)
      RV1=  CV1B* (P(2, 1)+P(2, 2)*X+P(2, 3)*X2+P(2, 4)*X3+P(2, 5)*X4
     &                    +P(2, 6)*Y+P(2, 7)*Y2+P(2, 8)*Y3
     &                    +P(14,1)*XY+P(14,2)*X2Y+P(14,3)*XY2)
      AV1=  AV1B* (P(3, 1)+P(3, 2)*X+P(3, 3)*X2+P(3, 4)*X3+P(3, 5)*X4
     &                    +P(3, 6)*Y+P(3, 7)*Y2+P(3, 8)*Y3
     &                    +P(14,4)*XY+P(14,5)*X2Y+P(14,6)*XY2)
      SUMI=-15.0* (P(4, 1)+P(4, 2)*X+P(4, 3)*X2+P(4, 4)*X3+P(4, 5)*X4
     &                    +P(4, 6)*Y+P(4, 7)*Y2+P(4, 8)*Y3
     &                    +P(19,4)*XY+P(19,5)*X2Y+P(19,6)*XY2)
      RV2=  CV2B* (P(5, 1)+P(5, 2)*X+P(5, 3)*X2+P(5, 4)*X3+P(5, 5)*X4
     &                    +P(5, 6)*Y+P(5, 7)*Y2+P(5, 8)*Y3
     &                    +P(15,1)*XY+P(15,2)*X2Y+P(15,3)*XY2)
      AV2=  AV2B* (P(6, 1)+P(6, 2)*X+P(6, 3)*X2+P(6, 4)*X3+P(6, 5)*X4
     &                    +P(6, 6)*Y+P(6, 7)*Y2+P(6, 8)*Y3
     &                    +P(15,4)*XY+P(15,5)*X2Y+P(15,6)*XY2)
      DIFFR=700.* (P(7, 1)+P(7, 2)*X+P(7, 3)*X2+P(7, 4)*X3+P(7, 5)*X4
     &                    +P(7, 6)*Y+P(7, 7)*Y2+P(7, 8)*Y3
     &                    +P(20,1)*XY+P(20,2)*X2Y+P(20,3)*XY2)
      RS1=  CS1B* (P(8, 1)+P(8, 2)*X+P(8, 3)*X2+P(8, 4)*X3+P(8, 5)*X4
     &                    +P(8, 6)*Y+P(8, 7)*Y2+P(8, 8)*Y3
     &                    +P(17,1)*XY+P(17,2)*X2Y+P(17,3)*XY2)
      AS1=  AS1B* (P(9, 1)+P(9, 2)*X+P(9, 3)*X2+P(9, 4)*X3+P(9, 5)*X4
     &                    +P(9, 6)*Y+P(9, 7)*Y2+P(9, 8)*Y3
     &                    +P(17,4)*XY+P(17,5)*X2Y+P(17,6)*XY2)
      DIFFI=-150.*(P(10,1)+P(10,2)*X+P(10,3)*X2+P(10,4)*X3+P(10,5)*X4
     &                    +P(10,6)*Y+P(10,7)*Y2+P(10,8)*Y3
     &                    +P(20,4)*XY+P(20,5)*X2Y+P(20,6)*XY2)
      RS2=  CS2B* (P(11,1)+P(11,2)*X+P(11,3)*X2+P(11,4)*X3+P(11,5)*X4
     &                    +P(11,6)*Y+P(11,7)*Y2+P(11,8)*Y3
     &                    +P(18,1)*XY+P(18,2)*X2Y+P(18,3)*XY2)
      AS2=  AS2B* (P(12,1)+P(12,2)*X+P(12,3)*X2+P(12,4)*X3+P(12,5)*X4
     &                    +P(12,6)*Y+P(12,7)*Y2+P(12,8)*Y3
     &                    +P(18,4)*XY+P(18,5)*X2Y+P(18,6)*XY2)
CC-----------------------------------------------
      VV=0.5*(SUMR+DIFFR)
      VS=0.5*(SUMR-DIFFR)
      WV=0.5*(SUMI+DIFFI)
      WS=0.5*(SUMI-DIFFI)
CC-------------------------------------------------------------------
C      VV2 =0.0
C      RV12=  1.0
C      AV12=  0.65
      WV2 = -100.0*(P(13,1)+P(13,2)*X+P(13,3)*X2+P(13,4)*X3+P(13,5)*X4
     &                     +P(13,6)*Y+P(13,7)*Y2+P(13,8)*Y3
     &                    +P(21,1)*XY+P(21,2)*X2Y+P(21,3)*XY2)
      RV22=  1.0
      AV22=  0.65
C
C      VS2 =0.0
C      RS12=  1.0
C      AS12=  0.65
      WS2 = 100.0*(P(16,1)+P(16,2)*X+P(16,3)*X2+P(16,4)*X3+P(16,5)*X4
     &                    +P(16,6)*Y+P(16,7)*Y2+P(16,8)*Y3
     &                    +P(21,4)*XY+P(21,5)*X2Y+P(21,6)*XY2)
      RS22=  1.0
      AS22=  0.65
CC
CC--- S.P  GEOMETRIES ARE FIXED BY THEIR VOLUME PART.
      RV22=RV2
      AV22=AV2
      RS22=RS2
      AS22=AS2
CC
      DSMALL=0.0
cx      DR=DMESH
      DO 10 I=1,NSTEP
         TMP1 = R(I)
            S1 = DSECH(RV1*ACB/AV1)
            S2 = DSECH(TMP1/AV1)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RVA1=RECV*VV*B/(A+B)
            S1 = DSECH(RV2*ACB/AV2)
            S2 = DSECH(TMP1/AV2)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RVA2=RECV*WV*B/(A+B)
            S1 = DSECH(RV22*ACB/AV22)
            S2 = DSECH(TMP1/AV22)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECV*WV2*A*B/(A+B)**2
         RVA2=RVA2+SURF
            S1 = DSECH(RS1*ACB/AS1)
            S2 = DSECH(TMP1/AS1)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RSA1=RECS*VS*B/(A+B)
            S1 = DSECH(RS2*ACB/AS2)
            S2 = DSECH(TMP1/AS2)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RSA2=RECS*WS*B/(A+B)
            S1 = DSECH(RS22*ACB/AS22)
            S2 = DSECH(TMP1/AS22)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECS*WS2*A*B/(A+B)**2
         RSA2=RSA2+SURF
CC
         VVA(I)=RVA1
         WVA(I)=RVA2
         VSA(I)=RSA1
         WSA(I)=RSA2
 10    CONTINUE
CC
      RETURN
      END
      SUBROUTINE EDAD3(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(NSTEP),VVA(NSTEP),WVA(NSTEP),VSA(NSTEP),WSA(NSTEP)
      DIMENSION P(22,8),PT1(8,8),PT2(8,8),PT3(8,6)
C23456789012345678901234567890123456789012345678901234567890123456789012
CC    A_DEP_FIT4D5.SPV ------ COSH.SPV  
CC     VOLUME = (COSH(R/A)-1)/(COSH(R/A)+COSH(X/A)-2)
CC     S.P = (COSH(R/A)-1)*(COSH(X/A)-1)/((COSH(R/A)+COSH(X/A)-2)**2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     X = E-DEP (EPS)   Y = A-DEP (CAA)
C--------------------------------------------------------------------
CCC   AADEP_BIG6.UNIT12 
       DATA PT1/
     & -2.368907E+01,  1.013198E+02, -1.692467E+02,  1.251361E+02,
     & -3.314906E+01,  1.591492E+01, -1.817534E+01,  3.102698E+00,
     &  4.471023E+00, -1.467462E+01,  2.652704E+01, -2.159166E+01,
     &  6.048228E+00, -2.441693E+00,  3.357096E+00, -4.663679E-01,
     &  1.360393E+01, -6.917896E+01,  1.382276E+02, -1.200362E+02,
     &  3.807423E+01, -2.314294E+00,  3.356424E+00, -1.618471E+00,
     & -7.578083E+01,  3.845006E+02, -7.474217E+02,  6.430927E+02,
     & -2.061245E+02,  1.481217E+01,  1.670286E+01,  5.139371E+00,
     &  4.858070E+00, -2.483372E+01,  4.738376E+01, -3.798950E+01,
     &  1.241031E+01,  3.797412E+00, -4.079241E+00,  7.813657E-01,
     & -1.061282E+01,  4.921347E+01, -6.624069E+01,  3.407809E+01,
     & -7.532193E+00,  6.490836E+00, -4.285493E+00,  8.012679E+00,
     &  1.384674E+00, -1.911635E+01,  4.990400E+01, -4.841091E+01,
     &  1.663529E+01,  6.709990E+00, -7.597097E+00,  2.632476E+00,
     &  1.941580E+00, -2.193781E+00,  3.824921E+00, -4.261957E+00,
     &  1.500230E+00, -2.263988E+00,  5.146258E-01,  1.090186E+00/
       DATA PT2/
     &  8.353493E+00, -4.280869E+01,  9.318841E+01, -8.662498E+01,
     &  2.866569E+01, -1.905610E+00,  3.064385E+00, -2.802112E-01,
     & -3.811906E+01,  2.261281E+02, -4.654927E+02,  4.283301E+02,
     & -1.461205E+02,  7.697116E+00, -1.966631E+00, -1.357537E+01,
     &  2.125486E+00, -3.189404E+01,  6.931285E+01, -5.498451E+01,
     &  1.687566E+01,  1.836083E+01, -1.332512E+01,  1.437681E+00,
     & -1.047202E+01,  3.846231E+01, -2.626984E+01, -1.075669E+01,
     &  8.458622E+00,  3.787689E+00,  5.592298E+00,  5.617703E+00,
     & -1.423291E+00, -1.395066E+01,  2.258138E+01, -4.070188E+00,
     & -4.088940E+00,  1.035303E+01, -1.580248E+01,  1.119403E+01,
     &  5.562317E-01,  2.325463E+00, -2.954182E+00,  2.352114E+00,
     & -1.042135E+00, -4.944172E-01, -3.599670E+01,  1.765904E+02,
     & -2.695062E-02, -4.119186E+00,  2.927944E+00, -1.816234E+01,
     &  2.156224E+01, -1.128396E+01, -3.421408E+02,  2.777086E+02,
     & -7.121706E+00,  9.348271E+01, -2.178505E+02,  2.041844E+02,
     & -7.165217E+01, -4.048138E+01,  3.144254E+01,  3.474007E+00/
       DATA PT3/ 
     &  3.195839E+00,  4.199934E-01, -2.636037E+00, -2.993818E-01,
     &  1.999406E+00, -2.390207E+00, -7.564810E+01,  1.075273E+01,
     & -1.578821E+01, -1.756925E+00,  1.082204E+01, -2.416553E+01,
     &  2.778435E+01, -1.707945E+01, -2.354979E+01,  3.406717E+00,
     & -8.497654E+00, -4.355649E+00,  1.272074E+01, -4.166584E+01,
     &  3.522784E+01, -2.837601E+01,  1.820885E+01, -2.990817E+01,
     & -3.796977E+00,  4.023453E-01,  2.400880E+00, -2.177108E+01,
     & -5.928764E+00,  3.136848E+01,  1.917163E+01,  4.568120E+00,
     & -3.922251E+00, -4.954481E-02, -6.934859E-02,  3.201427E+01,
     &  5.588527E+00, -3.304645E+01, -3.255172E+00, -2.388756E+01,
     &  4.414031E+01, -2.018493E+01, -4.662180E+00,  4.483327E+01,
     & -2.783741E+01, -3.675735E+01,  1.486115E+01,  8.281107E+00/
CC
      DSECH(X) = 2.*DEXP(-DABS(X))/( 1.0D0 + DEXP(-DABS(X))**2 )
C
      IF(TPLAB.LT.21.0D0.OR.TPLAB.GT.1040.D0) WRITE(6,601)
 601  FORMAT(' ',' YOU ARE EXTRAPOLATING IN ENERGY !!?? '/)
      IF(AA.LT.12.D0.OR.AA.GT.208.D0) WRITE(6,602)
 602  FORMAT(' ',' YOU ARE EXTRAPOLATING IN MASS NUMBER !!?? '/)
C
CCC
      DO 20 I=1,8
         DO 30 J=1,8
            P(J,I)=PT1(I,J)
 30      CONTINUE
         DO 40 J=1,8
            P(8+J,I)=PT2(I,J)
 40      CONTINUE
         DO 50 J=1,6
            P(16+J,I)=PT3(I,J)
 50      CONTINUE
 20   CONTINUE
CCC
C      PI=3.141592653589793D0
      AMU=931.5016D0
      WT=AA*AMU
      EE=TPLAB
      ACB=AA**0.33333333333333333D0
CC
      CAA=AA/(AA+20.)
      Y =CAA
      Y2=Y**2
      Y3=Y**3
C      Y4=Y**4
CC
CC--------------------------------------------------------------
CC TRANSFORM TO 2 BODY CM FRAME AND CALCULATE RECOIL FACTORS
CC
        WP=1.0072545D0*AMU
C        CWT=AA**0.3333333333333333D0
        EL=EE+WP
        WP2=WP**2
        WT2=WT**2
        PCM=DSQRT(WT2*(EL**2-WP2)/(WP2+WT2+2.0D0*WT*EL))
        EPCM=DSQRT(WP2+PCM**2)
        ETCM=DSQRT(WT2+PCM**2)
        SR=EPCM+ETCM
C        EKINCM=SR-WP-WT
C        TCM=EKINCM
CC
        E=1000.0/EPCM
        X=E
        X2=X**2
        X3=X**3
        X4=X**4
        XY=X*Y
        X2Y=X*X*Y
        XY2=X*Y*Y
        RECV=(ETCM /SR)
        RECS=(WT/SR)
        WRITE(6,111) E,ETCM,WT,SR,RECV,RECS
 111     FORMAT('  E,ETCM,WT,SR,RECV,RECS',/,1X,6(1PE10.3))
CC-------------------------------------------------------------------
C23456789012345678901234567890123456789012345678901234567890123456789012
      CV1B=1.0
      AV1B=0.7
      CV2B=CV1B
      CS1B=CV1B
      CS2B=CV1B
      AV2B=AV1B
      AS1B=AV1B
      AS2B=AV1B
CC-------------------------------------------------------------------
      SUMR=-100.* (P(1, 1)+P(1, 2)*X+P(1, 3)*X2+P(1, 4)*X3+P(1, 5)*X4
     &                    +P(1, 6)*Y+P(1, 7)*Y2+P(1, 8)*Y3
     &                    +P(19,1)*XY+P(19,2)*X2Y+P(19,3)*XY2)
      RV1=  CV1B* (P(2, 1)+P(2, 2)*X+P(2, 3)*X2+P(2, 4)*X3+P(2, 5)*X4
     &                    +P(2, 6)*Y+P(2, 7)*Y2+P(2, 8)*Y3
     &                    +P(14,1)*XY+P(14,2)*X2Y+P(14,3)*XY2)
      AV1=  AV1B* (P(3, 1)+P(3, 2)*X+P(3, 3)*X2+P(3, 4)*X3+P(3, 5)*X4
     &                    +P(3, 6)*Y+P(3, 7)*Y2+P(3, 8)*Y3
     &                    +P(14,4)*XY+P(14,5)*X2Y+P(14,6)*XY2)
      SUMI=-15.0* (P(4, 1)+P(4, 2)*X+P(4, 3)*X2+P(4, 4)*X3+P(4, 5)*X4
     &                    +P(4, 6)*Y+P(4, 7)*Y2+P(4, 8)*Y3
     &                    +P(19,4)*XY+P(19,5)*X2Y+P(19,6)*XY2)
      RV2=  CV2B* (P(5, 1)+P(5, 2)*X+P(5, 3)*X2+P(5, 4)*X3+P(5, 5)*X4
     &                    +P(5, 6)*Y+P(5, 7)*Y2+P(5, 8)*Y3
     &                    +P(15,1)*XY+P(15,2)*X2Y+P(15,3)*XY2)
      AV2=  AV2B* (P(6, 1)+P(6, 2)*X+P(6, 3)*X2+P(6, 4)*X3+P(6, 5)*X4
     &                    +P(6, 6)*Y+P(6, 7)*Y2+P(6, 8)*Y3
     &                    +P(15,4)*XY+P(15,5)*X2Y+P(15,6)*XY2)
      DIFFR=700.* (P(7, 1)+P(7, 2)*X+P(7, 3)*X2+P(7, 4)*X3+P(7, 5)*X4
     &                    +P(7, 6)*Y+P(7, 7)*Y2+P(7, 8)*Y3
     &                    +P(20,1)*XY+P(20,2)*X2Y+P(20,3)*XY2)
      RS1=  CS1B* (P(8, 1)+P(8, 2)*X+P(8, 3)*X2+P(8, 4)*X3+P(8, 5)*X4
     &                    +P(8, 6)*Y+P(8, 7)*Y2+P(8, 8)*Y3
     &                    +P(17,1)*XY+P(17,2)*X2Y+P(17,3)*XY2)
      AS1=  AS1B* (P(9, 1)+P(9, 2)*X+P(9, 3)*X2+P(9, 4)*X3+P(9, 5)*X4
     &                    +P(9, 6)*Y+P(9, 7)*Y2+P(9, 8)*Y3
     &                    +P(17,4)*XY+P(17,5)*X2Y+P(17,6)*XY2)
      DIFFI=-150.*(P(10,1)+P(10,2)*X+P(10,3)*X2+P(10,4)*X3+P(10,5)*X4
     &                    +P(10,6)*Y+P(10,7)*Y2+P(10,8)*Y3
     &                    +P(20,4)*XY+P(20,5)*X2Y+P(20,6)*XY2)
      RS2=  CS2B* (P(11,1)+P(11,2)*X+P(11,3)*X2+P(11,4)*X3+P(11,5)*X4
     &                    +P(11,6)*Y+P(11,7)*Y2+P(11,8)*Y3
     &                    +P(18,1)*XY+P(18,2)*X2Y+P(18,3)*XY2)
      AS2=  AS2B* (P(12,1)+P(12,2)*X+P(12,3)*X2+P(12,4)*X3+P(12,5)*X4
     &                    +P(12,6)*Y+P(12,7)*Y2+P(12,8)*Y3
     &                    +P(18,4)*XY+P(18,5)*X2Y+P(18,6)*XY2)
CC-----------------------------------------------
      VV=0.5*(SUMR+DIFFR)
      VS=0.5*(SUMR-DIFFR)
      WV=0.5*(SUMI+DIFFI)
      WS=0.5*(SUMI-DIFFI)
CC-------------------------------------------------------------------
      VV2 =  100.0*(P(14,7)+P(14,8)*X+P(15,7)*X2+P(15,8)*X3+P(17,7)*X4
     &                     +P(17,8)*Y+P(18,7)*Y2+P(18,8)*Y3
     &                    +P(19,7)*XY+P(19,8)*X2Y+P(20,7)*XY2)
      RV12=  1.0
      AV12=  0.65
      WV2 = -100.0*(P(13,1)+P(13,2)*X+P(13,3)*X2+P(13,4)*X3+P(13,5)*X4
     &                     +P(13,6)*Y+P(13,7)*Y2+P(13,8)*Y3
     &                    +P(21,1)*XY+P(21,2)*X2Y+P(21,3)*XY2)
      RV22=  1.0
      AV22=  0.65
C
      VS2 = -100.0*(P(20,8)+P(21,7)*X+P(21,8)*X2+P(22,1)*X3+P(22,2)*X4
     &                     +P(22,3)*Y+P(22,4)*Y2+P(22,5)*Y3
     &                    +P(22,6)*XY+P(22,7)*X2Y+P(22,8)*XY2)
      RS12=  1.0
      AS12=  0.65
      WS2 = 100.0*(P(16,1)+P(16,2)*X+P(16,3)*X2+P(16,4)*X3+P(16,5)*X4
     &                    +P(16,6)*Y+P(16,7)*Y2+P(16,8)*Y3
     &                    +P(21,4)*XY+P(21,5)*X2Y+P(21,6)*XY2)
      RS22=  1.0
      AS22=  0.65
CC
CC--- S.P  GEOMETRIES ARE FIXED BY THEIR VOLUME PART.
      RV12=RV1
      AV12=AV1
      RS12=RS1
      AS12=AS1
      RV22=RV2
      AV22=AV2
      RS22=RS2
      AS22=AS2
CC
      DSMALL=0.0
cx      DR=DMESH
      DO 10 I=1,NSTEP
         TMP1 = R(I)
            S1 = DSECH(RV1*ACB/AV1)
            S2 = DSECH(TMP1/AV1)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RVA1=RECV*VV*B/(A+B)
            S1 = DSECH(RV12*ACB/AV12)
            S2 = DSECH(TMP1/AV12)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECV*VV2*A*B/(A+B)**2
         RVA1=RVA1+SURF
            S1 = DSECH(RV2*ACB/AV2)
            S2 = DSECH(TMP1/AV2)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RVA2=RECV*WV*B/(A+B)
            S1 = DSECH(RV22*ACB/AV22)
            S2 = DSECH(TMP1/AV22)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECV*WV2*A*B/(A+B)**2
         RVA2=RVA2+SURF
            S1 = DSECH(RS1*ACB/AS1)
            S2 = DSECH(TMP1/AS1)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RSA1=RECS*VS*B/(A+B)
            S1 = DSECH(RS12*ACB/AS12)
            S2 = DSECH(TMP1/AS12)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECV*VS2*A*B/(A+B)**2
         RSA1=RSA1+SURF
            S1 = DSECH(RS2*ACB/AS2)
            S2 = DSECH(TMP1/AS2)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RSA2=RECS*WS*B/(A+B)
            S1 = DSECH(RS22*ACB/AS22)
            S2 = DSECH(TMP1/AS22)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECS*WS2*A*B/(A+B)**2
         RSA2=RSA2+SURF
CC
         VVA(I)=RVA1
         WVA(I)=RVA2
         VSA(I)=RSA1
         WSA(I)=RSA2
 10    CONTINUE
CC
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCC HEPB1 FITS BELOW
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE HEPB1(IFIT,TPLAB,AA,R,VVA,WVA,VSA,WSA,NSTEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(NSTEP),VVA(NSTEP),WVA(NSTEP),VSA(NSTEP),WSA(NSTEP)
      DIMENSION P(32,16)
      DIMENSION Z1(16),Z2(16),Z3(16),Z4(16),Z5(16),Z6(16),Z7(16),Z8(16)
      DIMENSION Z9(16),Z10(16),Z11(16),Z12(16),Z13(16),Z14(16),Z15(16)
      DIMENSION Z16(16),Z17(16),Z18(16),Z19(16),Z20(16),Z32(16)
      DIMENSION SH1(16),SH2(16),SH3(16),SH4(16),SH5(16),SH6(16)
      DIMENSION SH7(16),SH8(16),SH9(16),SH10(16),SH11(16),SH12(16)
      DIMENSION SH13(16),SH14(16),SH15(16)
      DIMENSION SH16(16),SH17(16),SH18(16),SH19(16),SH20(16),SH32(16)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CC    TIM'S UNDEMOCRATIC FIT  
CC     VOLUME = (COSH(R/A)-1)/(COSH(R/A)+COSH(X/A)-2)
CC     S.P = (COSH(R/A)-1)*(COSH(X/A)-1)/((COSH(R/A)+COSH(X/A)-2)**2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     X = E-DEP (EPS)   Y = A-DEP (CAA)
C--------------------------------------------------------------------
C 
      DSECH(X) = 2.*DEXP(-DABS(X))/( 1.0D0 + DEXP(-DABS(X))**2 )
C
C THE PARAMETERS THAT SHINICHI TWEAKED TO FINISH OFF.
C
      DATA Z1/
     &1.210009265D+00,2.164440631D-01,-1.632087915D-02,5.066222234D-02,
     &9.966173735D-02,1.729965500D-01,-1.113903710D-01,3.522140699D-02,
     &-4.588215380D-02,8.959125073D-02,9.421916258D-02,1.589589431D-01,
     &-1.275691609D-01,3.483462811D-02,-1.817554851D-01,-1.551888059D-1/
      DATA Z2/
     &7.240002066D-01,7.646860239D-03,-5.880651938D-02,-5.496517698D-02,
     &-1.438780212D-03,3.938464456D-02,6.459432505D-02,-8.100095939D-02,
     &7.668355577D-02,-3.463832697D-02,-1.336492177D-01,-6.747536604D-2,
     &1.801139117D-01,3.441068233D-02,7.886862457D-02,-5.082788875D-02/
      DATA Z3/
     &6.705349424D-01,-4.948276028D-02,8.940624656D-02,1.487242371D-01,
     & 5.408660994D-02,2.156668384D-01,-1.590567064D-01,3.159123699D-02,
     &-2.455568022D-03,3.002819319D-02,4.215350421D-02,9.171161981D-02,
     &-1.077904172D-01,1.610423516D-02,-8.649927007D-02,-3.218865927D-2/
      DATA Z4/
     &-1.389439211D+00,-1.674766667D+00,1.086812781D+00,7.731116915D-01,
     &8.281628954D-02,-1.570539120D+00,-1.299727193D+00,6.986198000D-01,
     &-1.505613593D+0,-5.198446507D-01,-6.996185145D-01,1.110163365D+00,
     &-1.324743287D+00,1.360097624D+00,1.833586323D+00,8.712987357D-01/
      DATA Z5/
     &2.618760610D-01,-7.439706896D-01,1.848591694D-01,1.065010884D-01,
     &9.309755283D-02,7.717715817D-02,6.744522800D-02,1.558625722D-02,
     &-2.422608451D-01,3.693103613D-01,7.617981797D-01,-2.717065334D-02,
     &-6.197575003D-01,-5.913409262D-02,2.475044261D-01,1.397198725D-01/
      DATA Z6/
     &7.772742060D-01,-2.006192277D-02,-3.466697033D-03,1.042060034D-01,
     &4.272258913D-02,2.666106279D-01,-3.563810061D-01,3.397640492D-02,
     &2.258757651D-01,-1.187993136D-01,1.639472506D-01,9.492328349D-02,
     &-3.054804542D-01,8.088311832D-03,-8.760803717D-02,1.075867069D-01/
      DATA Z7/
     &7.598827896D-01,-2.041734807D-02,-1.427708304D-1,-9.953086766D-02,
     &-9.72674816D-03,-2.417106012D-02,-8.226193185D-2,-1.633278524D-01,
     &2.345022654D-01,-1.965641605D-02,-2.082703303D-1,-3.483206088D-01,
     &-3.040756751D-01,-2.418069139D-01,4.536684171D-2,-1.091336848D-01/
      DATA Z8/
     &-5.092162979D-02,-4.930163562D-01,-2.9392423D-01,-2.739511910D-01,
     &-1.476426007D-01,-8.568171415D-01,7.661317173D-01,1.451300153D-01,
     &-8.964460353D-01,5.957433293D-01,-3.536651462D-1,-4.621628819D-01,
     & 3.114207591D-01,-3.276819890D-01,1.232954220D-01,1.182752629D-01/
      DATA Z9/
     &1.057536013D+00,3.807972635D-01,-1.420442056D-01,1.756531133D-01,
     &6.253517906D-02,-5.862139577D-02,5.186685622D-02,-1.520533985D-01,
     &9.379220488D-02,4.439871569D-02,5.931581678D-02,1.353251355D-01,
     &-9.488799943D-02,-6.586252719D-02,-2.918810293D-1,-2.059297763D-1/
      DATA Z10/
     &7.268561779D-01,1.287101206D-02,-5.069605444D-02,-3.897316418D-02,
     &-5.409464864D-03,2.905530364D-02,7.002828545D-02,-8.081606404D-02,
     &7.037408381D-02,-3.169516749D-02,-1.269405753D-01,-6.456187994D-2,
     &1.668402114D-01,1.659313898D-02,6.667738366D-02,-5.036794620D-02/
      DATA Z11/
     &7.306733408D-01,-1.960273597D-02,7.117359780D-02,1.140564383D-01,
     &2.560041135D-02,2.197110310D-01,-1.384013040D-01,1.889366471D-02,
     &-6.248092298D-04,2.775593077D-02,4.486727999D-02,8.599987371D-02,
     &-7.219179535D-02,5.864124698D-03,-5.058772453D-02,-5.710487284D-2/
      DATA Z12/
     &-1.574117945D+00,-1.943241207D+00,1.079285745D+00,5.653756520D-01,
     &1.841524592D-01,-1.229207967D+00,-1.314543415D+00,1.060765277D+00,
     &-1.25583451D+00,-4.790777349D-01,-7.595398752D-01,7.756273716D-01,
     &-1.219690406D+00,1.256730881D+00,1.850250890D+00,1.409690501D+00/
      DATA Z13/ 
     &4.133068767D-01,-5.966002718D-01,8.718250230D-02,-4.187273977D-02,
     &9.781736018D-02,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00/
      DATA Z14/
     &8.667204897D-01,4.843814744D-02,1.943901230D-02,7.277630631D-02,
     &4.367025891D-03,2.421222236D-01,-3.243050729D-01,2.051539254D-02,
     &1.576666310D-01,-1.242266442D-01,1.317526294D-01,9.099633318D-02,
     &-1.308920875D-01,-5.438766318D-02,-1.74516131D-01,-1.28198623D-02/
      DATA Z15/
     &6.499862089D-01,-5.540368374D-02,-3.310148723D-3,-4.371206135D-02,
     &-4.002163791D-02,3.484448853D-01,-1.668995676D-01,-1.60938085D-01,
     &4.161664958D-01,-1.749206862D-01,4.026344996D-01,4.669397793D-02,
     &-7.363547319D-02,-1.746537817D-01,1.700586666D-01,-7.81345462D-02/
      DATA Z16/
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00/
      DATA Z17/
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00/
      DATA Z18/
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00/
      DATA Z19/
     &-1.054108004D+00,4.342367218D-01,3.220792900D-01,1.411454528D+00,
     &-4.049756013D-01,-1.567243470D+00,5.219996087D-01,7.775902613D-02,
     &-5.335183251D-01,1.749402180D+00,-3.168476073D+0,-8.720883944D-01,
     &-1.699834025D+00,-5.194345518D-01,2.137403944D+0,-5.383190860D-01/
      DATA Z20/
     &3.133973960D-01,1.923209701D+00,2.870597182D-01,-5.287319875D-01,
     &-5.150448597D-01,-7.392983844D-01,3.927585488D-01,4.346482813D-01,
     &-2.282140040D-01,-4.254069733D-01,-3.390948565D+0,-1.253381272D+0,
     &2.566089335D+00,5.642504009D-01,9.445967945D-01,-9.934121111D-01/
C 
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C 
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C 
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C 
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C 
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C 
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C 
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C 
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C 
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C 
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C 
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
      DATA Z32/ 
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00/

C   Above are final FIT1 parameters from re-start-log-4-log-fit1.TXT   12/02/2007
C 
C   Below are THE PARAMETERS BELOW ARE FOR sHINICHI'S DEMOCRATIC FIT.
C
      DATA SH1/
     &1.255823718D+00,2.144542531D-01,-7.982000662D-02,-7.519529022D-02,
     &-2.888255545D-02,4.098492271D-02,1.762742074D-01,-7.430210820D-02,
     &-1.651996035D-01,1.042334867D-01,-1.026634130D-01,9.624749493D-02,
     &-3.318421766D-02,2.000644894D-01,-1.634272474D-01,-8.697887844D-2/
       DATA SH2/
     &7.442249727D-01,7.743629946D-02,3.325018138D-02,3.855402908D-03,
     &1.137009217D-03,7.788037460D-02,-4.884314931D-02,3.702100571D-03,
     &1.775969555D-02,-7.613280931D-03,-1.205083703D-01,-1.533451674D-1,
     &1.857152958D-02,-5.877077508D-03,9.406579298D-02,4.172911406D-02/
       DATA SH3/
     &6.579792559D-01,-5.349300459D-02,7.728961454D-02,1.252019569D-01,
     &5.143735432D-02,1.583674196D-01,-1.142531043D-01,2.731127221D-02,
     &2.013747259D-02,2.085846203D-02,-1.262168992D-01,-8.271072075D-02,
     &-7.312651145D-03,-3.942923256D-02,1.315342766D-03,-4.116082876D-2/
       DATA SH4/
     &-2.363365236D+0,-3.461003447D+0,-3.285042774D-01,-2.511375997D-01,
     &5.221618618D-02,-3.152970018D+00,-1.971500562D+00,9.646501061D-02,
     &1.136628963D-01,-5.095659939D-01,-2.333218232D+00,1.057270697D+00,
     &-2.163961089D+00,8.508835406D-01,2.187741423D-02,-1.079140423D-01/
       DATA SH5/
     &3.313479491D-01,-4.138199429D-01,3.879415821D-01,-1.921806680D-01,
     &1.276727616D-01,-4.279560777D-01,7.327066166D-01,-3.677208520D-01,
     &5.633771106D-02,2.451031443D-01,4.100417986D-01,-2.067603551D-01,
     &6.299351168D-02,1.973632158D-01,2.822140693D-01,-4.599065349D-01/
       DATA SH6/
     &6.723008295D-01,-1.135292374D-01,7.229324617D-02,5.712671874D-02,
     &-9.563697413D-03,2.825092368D-1,-2.747358372D-01,-6.929816646D-02,
     &1.636619261D-01,-7.972082590D-2,1.856832446D-01,-6.992123545D-02,
     &-2.481392639D-01,-7.077252279D-2,-9.117186990D-2,-3.477867994D-02/
       DATA SH7/
     &7.483745497D-01,7.110504851D-02,-1.743802604D-01,-8.414817501D-02,
     &6.689752879D-03,9.044546136D-02,-9.050903727D-02,3.262809477D-02,
     &1.059922172D-01,-1.667714026D-02,-2.777223586D-1,-2.692514546D-01,
     &1.035681958D-01,-1.313363792D-01,1.965236740D-01,-3.521642580D-02/
       DATA SH8/
     &4.189583682D-01,-3.421484518D-02,-4.256645266D-1,-1.064505615D-01,
     &2.240790983D-01,-6.009298445D-01,5.343650438D-01,-1.377573766D-01,
     &-7.635281603D-01,5.173255917D-01,1.731085791D-01,1.154236088D-01,
     &-3.606352594D-02,-1.551175348D-01,1.248206589D-01,7.496019707D-02/
       DATA SH9/
     &1.103267187D+00,4.188630297D-01,-1.552475283D-01,9.169946919D-02,
     &-3.823230731D-02,-1.763588784D-01,2.559768187D-1,-2.290887612D-01,
     &-1.790197114D-02,5.594007706D-02,-1.476343327D-01,2.918742535D-02,
     &2.686460055D-02,9.781347678D-03,-2.326051439D-01,-2.028161019D-01/
       DATA SH10/
     &7.357449157D-01,6.931820895D-02,2.344200487D-02,9.275189367D-03,
     &-9.980694292D-04,8.842858860D-02,-4.999236591D-2,-3.067849266D-03,
     &2.624480528D-02,-1.261042989D-02,-9.509460793D-2,-1.252370655D-01,
     &9.543469707D-03,-1.242747958D-02,7.792240216D-02,3.560009178D-02/
       DATA SH11/
     &7.186018689D-01,-2.527718150D-02,5.902889233D-02,9.890917525D-02,
     &3.490272258D-02,1.734589738D-01,-9.638642142D-02,1.704954798D-02,
     &1.878821929D-02,1.978106214D-02,-9.875253512D-02,-5.106906966D-02,
     &1.801364547D-02,-2.461938519D-02,1.967975757D-02,-5.482663894D-02/
        DATA SH12/
     &-2.873955517D+00,-4.359800538D+0,-3.053883404D-1,-9.434892596D-02,
     &3.770792037D-01,-4.443034656D+0,-3.058089245D+00,-9.625683909D-02,
     &1.092442629D-01,-4.573974270D-1,-4.918509911D+00,-6.207772999D-01,
     &-3.926103715D+00,6.502868236D-1,-9.085656131D-01,-4.787328690D-01/
       DATA SH13/
     &4.477237643D-01,-6.836773934D-01,-6.360112835D-02,8.444229710D-02,
     &-5.958418854D-02,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00/
       DATA SH14/
     &8.950546568D-01,1.176020928D-01,4.456044915D-02,-4.898991323D-02,
     &-2.139964784D-02,1.171645472D-02,-1.971165751D-01,6.184775882D-02,
     &5.502256943D-02,-5.630624901D-2,-3.192625392D-01,-1.442878659D-01,
     &1.134012615D-01,2.356683023D-02,-7.006727372D-02,-3.707414694D-02/
       DATA SH15/
     &7.009442784D-01,2.041609312D-01,-6.148906355D-03,2.321980697D-02,
     &-2.090430629D-02,5.249030162D-01,-2.777147236D-01,5.187778178D-02,
     &1.771422601D-01,-5.682470810D-02,8.558162206D-02,2.910037266D-02,
     &4.631471406D-01,-1.201336908D-01,2.693829394D-01,-1.826971385D-01/
       DATA SH16/
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00/
       DATA SH17/
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00/
       DATA SH18/
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00/
       DATA SH19/
     &-1.660261984D+00,1.200081642D-01,2.355905301D+00,6.839658567D-01,
     &-5.575955777D-01,-2.414911764D+0,3.622083848D+00,-3.842793095D+00,
     &2.327142748D+00,2.550724181D-01,-9.809073643D-02,-4.138249354D+00,
     &-1.057222318D+00,-2.748924995D+0,4.629142583D+00,-1.222202263D+00/
       DATA SH20/
     &-7.620562250D-01,-1.032597535D-1,2.574145924D-01,-2.544282479D-01,
     &-3.887437967D-01,7.028102777D-01,-6.668010684D-01,7.489436945D-01,
     &2.515123994D-01,-6.622390750D-01,7.953060319D-01,7.199382077D-01,
     &-1.329122313D+00,9.297090651D-01,-3.938619999D-01,2.355345600D-01/
C       DATA SH1/
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C       DATA SH1/
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C       DATA SH1/
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C       DATA SH1/
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C       DATA SH1/
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C       DATA SH1/
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C       DATA SH1/
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
CC    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C       DATA SH1/
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C       DATA SH1/
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C        DATA SH1/    V
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C       DATA SH1/
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
C    0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00
       DATA SH32/
     &0.0D+00,0.0D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00,
     &0.000000000D+00,0.000000000D+00,0.000000000D+00,0.000000000D+00/
CC
Cfinal undemocratic fit 2,from  log-tim-undemo  8.txt    12/16/2007 
C    last 2,entries 8.097746149D-01,3.990023189D-01,set to zero,since they are not used
C
C==========END OF DATA!=============
C
C LEGENDRE POLYNOMIALS...
C
C      PN0(Z)=1
C      PN1(Z)=Z
      PN2(Z)=(3*Z**2-1)/2.0D0
      PN3(Z)=(5*Z**3 - 3*Z)/2.0D0
      PN4(Z)=(35*Z**4 - 30*Z**2 + 3)/8.0D0
      PN5(Z)=(63*Z**5-70*Z**3+15*Z)/8.0D0
C     PN6(Z)=(231*Z**6-315*Z**4+105*Z**2 -5)/16.0D0
C
      IF(TPLAB.LT.21.0D0.OR.TPLAB.GT.1040.D0) WRITE(6,601)
 601  FORMAT(' ',' DANGER: YOU ARE EXTRAPOLATING IN ENERGY !! '/)
      IF(AA.LT.4.D0.OR.AA.GT.208.D0) WRITE(6,602)
 602  FORMAT(' ',' DANGER: YOU ARE EXTRAPOLATING IN MASS NUMBER !! '/)
C
CCC NOW LOAD UP THE BIG MATRIX.
      IF(IFIT.EQ.11) THEN
       DO 5 I=1,16
        P(1,I)=Z1(I)
        P(2,I)=Z2(I)
        P(3,I)=Z3(I)
        P(4,I)=Z4(I)
        P(5,I)=Z5(I)
        P(6,I)=Z6(I)
        P(7,I)=Z7(I)
        P(8,I)=Z8(I)
        P(9,I)=Z9(I)
        P(10,I)=Z10(I)
        P(11,I)=Z11(I)
        P(12,I)=Z12(I)
        P(13,I)=Z13(I)
        P(14,I)=Z14(I)
        P(15,I)=Z15(I)
        P(16,I)=Z16(I)
        P(17,I)=Z17(I)
        P(18,I)=Z18(I)
        P(19,I)=Z19(I)
        P(20,I)=Z20(I)
        P(32,I)=Z32(I)
5      CONTINUE
      END IF
      IF(IFIT.EQ.12) THEN
        DO 6 I=1,16
        P(1,I)=SH1(I)
        P(2,I)=SH2(I)
        P(3,I)=SH3(I)
        P(4,I)=SH4(I)
        P(5,I)=SH5(I)
        P(6,I)=SH6(I)
        P(7,I)=SH7(I)
        P(8,I)=SH8(I)
        P(9,I)=SH9(I)
        P(10,I)=SH10(I)
        P(11,I)=SH11(I)
        P(12,I)=SH12(I)
        P(13,I)=SH13(I)
        P(14,I)=SH14(I)
        P(15,I)=SH15(I)
        P(16,I)=SH16(I)
        P(17,I)=SH17(I)
        P(18,I)=SH18(I)
        P(19,I)=SH19(I)
        P(20,I)=SH20(I)
        P(32,I)=SH32(I)
6      CONTINUE
      END IF
CCC
      PI=3.141592653589793D0
      AMU=931.5016D0
      WT=AA*AMU
      EE=TPLAB
      ACB=AA**0.33333333333333333D0
CC
C      CAA=AA/(AA+20.)
C      Y =CAA
C      Y2=Y**2
C      Y3=Y**3
C      Y4=Y**4
CC
CC--------------------------------------------------------------
CC TRANSFORM TO 2 BODY CM FRAME AND CALCULATE RECOIL FACTORS
CC
        WP=1.0072545D0*AMU
C        CWT=AA**0.3333333333333333D0
        EL=EE+WP
        WP2=WP**2
        WT2=WT**2
        PCM=DSQRT(WT2*(EL**2-WP2)/(WP2+WT2+2.0D0*WT*EL))
        EPCM=DSQRT(WP2+PCM**2)
        ETCM=DSQRT(WT2+PCM**2)
        SR=EPCM+ETCM
        EKINCM=SR-WP-WT
        TCM=EKINCM
CC
C        E=1000.0/EPCM
C        X=E
C        X2=X**2
C        X3=X**3
C        X4=X**4
C        XY=X*Y
C        X2Y=X*X*Y
C        XY2=X*Y*Y
        RECV=(ETCM /SR)
        RECS=(WT/SR)
C        WRITE(6,111) E,ETCM,WT,SR,RECV,RECS
C 111     FORMAT('  E,ETCM,WT,SR,RECV,RECS',/,1X,6(1PE10.3))
C

C
C      AMU=931.5016D0
C      E=EA(NDAT)
C      WT=WTA(NDAT)
C      WT=AA*AMU
C      AA=WT/AMU
C      FZZ=ZT(NDAT)
C      FNN=AA-FZZ
C      FAA=AA  
      EE=TPLAB
C      WRITE(6,9871) TPLAB,FAA,FZZ,FNN
C9871  FORMAT(' TPLAB,FAA,FZZ,FNN',4F15.3)      
C      ACB=AA**0.33333333333333333D0
      ACB=AA**0.38D0
C
      FI=21.0D0
      YP=AA/(AA+FI)
      YMIN=4.0D0/(4.0D0+FI)
      YMAX=208.0D0/(208.0D0+FI)
C
      YMID=0.5D0*(YMIN+YMAX)
      Y = 2.0D0*(YP-YMID)/(YMAX-YMIN)
C      Y1=Y
      Y2=PN2(Y)
      Y3=PN3(Y)
      Y4=PN4(Y)
      Y5=PN5(Y)
CC
CC--------------------------------------------------------------
CC TRANSFORM TO 2 BODY CM FRAME AND CALCULATE RECOIL FACTORS
CC
C        WP=1.0072545D0*AMU
C        CWT=AA**0.3333333333333333D0
C        EL=EE+WP
C        WP2=WP**2
C        WT2=WT**2
C        PCM=DSQRT(WT2*(EL**2-WP2)/(WP2+WT2+2.0D0*WT*EL))
C        EPCM=DSQRT(WP2+PCM**2)
C        ETCM=DSQRT(WT2+PCM**2)
C        SR=EPCM+ETCM
C        EKINCM=SR-WP-WT
C        TCM=EKINCM
CC
C      XP=DLOG(TCM)/DLOG(1044.0D0)
C      XMIN=DLOG(21.0D0)/DLOG(1044.0D0)
C      XMAX=DLOG(1044.0D0)/DLOG(1044.0D0)
CC      ESCALE=100.0D0*P(32,15)
CC      XP=  DASINH(TCM     /ESCALE)/DASINH(1044.0D0/ESCALE)
CC      XMIN=DASINH(21.0D0  /ESCALE)/DASINH(1044.0D0/ESCALE)
CC      XMAX=DASINH(1044.0D0/ESCALE)/DASINH(1044.0D0/ESCALE)
      XP=1.0D0/DLOG(TCM)
      XMIN=1.0D0/DLOG(1040.0D0)
      XMAX=1.0D0/DLOG(19.0D0)
C
C      TSCALE=300.0D0*P(32,1)
C      IF(TSCALE.LT.20.0D0.OR.TSCALE.GT.1.D5) STOP 'TSCALE ERROR'
C      XP=  DEXP(-TCM     /TSCALE)
C      XMIN=DEXP(-1040.0D0/TSCALE)
C      XMAX=DEXP(  -19.0D0/TSCALE)
C
      XMID=0.5D0*(XMAX+XMIN)
      X=2.0D0*(XP-XMID)/(XMAX-XMIN)
C      X1=X
      X2=PN2(X)
      X3=PN3(X)
      X4=PN4(X)
C      X5=PN5(X)
C      X6=PN6(X)
C
      RECV=(ETCM /SR)
      RECS=(WT/SR)
C      RECVTAB(NDAT)=RECV
C      RECSTAB(NDAT)=RECS
C        RECV=1.0D0
C        RECS=1.0D0
C        WRITE(6,111) E,ETCM,WT,SR,RECV,RECS
C 111     FORMAT('  E,ETCM,WT,SR,RECV,RECS',/,1X,6(1PE10.3))
CC-------------------------------------------------------------------
C23456789012345678901234567890123456789012345678901234567890123456789012
C      CV1B=1.0
C      AV1B=0.7
C      CV2B=CV1B
C      CS1B=CV1B
C      CS2B=CV1B
C      AV2B=AV1B
C      AS1B=AV1B
C      AS2B=AV1B
CC NOW THE WOODS-SAXON PARAMETERS...
CC-------------------------------------------------------------------
C
C SUM OF REAL STRENGHTS
C
C STRAIGHT X-DEPENDENCE
      A0=P(1,1)
      A1=P(1,2)
      A2=P(1,3)
      A3=P(1,4)
      A4=P(1,5)
C      A5=P(1,6)
C      A6=P(1,7)
C
      B1 = A0 + A1*X + A2*X2 + A3*X3 + A4*X4 
C STRAIGHT Y DEPENDENCE
      A0=0.0D0
      A1=P(1,6)
      A2=P(1,7)
      A3=P(1,8)
      A4=P(1,9)
      A5=P(1,10)
C
      B2 = A0 + A1*Y + A2*Y2 + A3*Y3 + A4*Y4 + A5*Y5
CROSS TERMS
      A0=P(1,11)
      A1=P(1,12)
      A2=P(1,13)
      A3=P(1,14)
      A4=P(1,15)
      A5=P(1,16)
C
      B3 = A0*X*Y + A1*X2*Y + A2*X*Y2 + A3*X3*Y + A4*X2*Y2 + A5*X*Y3
C
      SUMR = -75.0D0*(B1+B2+B3)
C
C REAL VECTOR RADIUS
C
C STRAIGHT X-DEPENDENCE
      A0=P(2,1)
      A1=P(2,2)
      A2=P(2,3)
      A3=P(2,4)
      A4=P(2,5)
C      A5=P(1,6)
C      A6=P(1,7)
C
      B1 = A0 + A1*X + A2*X2 + A3*X3 + A4*X4 
C STRAIGHT Y DEPENDENCE
      A0=0.0D0
      A1=P(2,6)
      A2=P(2,7)
      A3=P(2,8)
      A4=P(2,9)
      A5=P(2,10)
C
      B2 = A0 + A1*Y + A2*Y2 + A3*Y3 + A4*Y4 + A5*Y5
CROSS TERMS
      A0=P(2,11)
      A1=P(2,12)
      A2=P(2,13)
      A3=P(2,14)
      A4=P(2,15)
      A5=P(2,16)
C
      B3 = A0*X*Y + A1*X2*Y + A2*X*Y2 + A3*X3*Y + A4*X2*Y2 + A5*X*Y3
C
      RV1 = B1+B2+B3
C
C REAL VECTOR DIFFUSENESS
C
C STRAIGHT X-DEPENDENCE
      A0=P(3,1)
      A1=P(3,2)
      A2=P(3,3)
      A3=P(3,4)
      A4=P(3,5)
C      A5=P(3,6)
C      A6=P(3,7)
C
      B1 = A0 + A1*X + A2*X2 + A3*X3 + A4*X4 
C STRAIGHT Y DEPENDENCE
      A0=0.0D0
      A1=P(3,6)
      A2=P(3,7)
      A3=P(3,8)
      A4=P(3,9)
      A5=P(3,10)
C
      B2 = A0 + A1*Y + A2*Y2 + A3*Y3 + A4*Y4 + A5*Y5
CROSS TERMS
      A0=P(3,11)
      A1=P(3,12)
      A2=P(3,13)
      A3=P(3,14)
      A4=P(3,15)
      A5=P(3,16)
C
      B3 = A0*X*Y + A1*X2*Y + A2*X*Y2 + A3*X3*Y + A4*X2*Y2 + A5*X*Y3
C
      AV1 = 0.7D0*(B1+B2+B3)
C
C REAL VECTOR PARABOLIC FERMI PARAMETER
C
C STRAIGHT X-DEPENDENCE
      A0=P(4,1)
      A1=P(4,2)
      A2=P(4,3)
      A3=P(4,4)
      A4=P(4,5)
C      A5=P(1,6)
C      A6=P(1,7)
C
      B1 = A0 + A1*X + A2*X2 + A3*X3 + A4*X4 
C STRAIGHT Y DEPENDENCE
      A0=0.0D0
      A1=P(4,6)
      A2=P(4,7)
      A3=P(4,8)
      A4=P(4,9)
      A5=P(4,10)
C
      B2 = A0 + A1*Y + A2*Y2 + A3*Y3 + A4*Y4 + A5*Y5
CROSS TERMS
      A0=P(4,11)
      A1=P(4,12)
      A2=P(4,13)
      A3=P(4,14)
      A4=P(4,15)
      A5=P(4,16)
C
      B3 = A0*X*Y + A1*X2*Y + A2*X*Y2 + A3*X3*Y + A4*X2*Y2 + A5*X*Y3
C
      WGV1 =(B1+B2+B3)
C
C SUM OF IMAGINARY STRENGHTS
C
C STRAIGHT X-DEPENDENCE
      A0=P(5,1)
      A1=P(5,2)
      A2=P(5,3)
      A3=P(5,4)
      A4=P(5,5)
C      A5=P(1,6)
C      A6=P(1,7)
C
      B1 = A0 + A1*X + A2*X2 + A3*X3 + A4*X4 
C STRAIGHT Y DEPENDENCE
      A0=0.0D0
      A1=P(5,6)
      A2=P(5,7)
      A3=P(5,8)
      A4=P(5,9)
      A5=P(5,10)
C
      B2 = A0 + A1*Y + A2*Y2 + A3*Y3 + A4*Y4 + A5*Y5
CROSS TERMS
      A0=P(5,11)
      A1=P(5,12)
      A2=P(5,13)
      A3=P(5,14)
      A4=P(5,15)
      A5=P(5,16)
C
      B3 = A0*X*Y + A1*X2*Y + A2*X*Y2 + A3*X3*Y + A4*X2*Y2 + A5*X*Y3
C
      SUMI = -20.D0*(B1+B2+B3)
C
C IMIGINARY VECTOR RADIUS
C
C STRAIGHT X-DEPENDENCE
      A0=P(6,1)
      A1=P(6,2)
      A2=P(6,3)
      A3=P(6,4)
      A4=P(6,5)
C      A5=P(1,6)
C      A6=P(1,7)
C
      B1 = A0 + A1*X + A2*X2 + A3*X3 + A4*X4 
C STRAIGHT Y DEPENDENCE
      A0=0.0D0
      A1=P(6,6)
      A2=P(6,7)
      A3=P(6,8)
      A4=P(6,9)
      A5=P(6,10)
C
      B2 = A0 + A1*Y + A2*Y2 + A3*Y3 + A4*Y4 + A5*Y5
CROSS TERMS
      A0=P(6,11)
      A1=P(6,12)
      A2=P(6,13)
      A3=P(6,14)
      A4=P(6,15)
      A5=P(6,16)
C
      B3 = A0*X*Y + A1*X2*Y + A2*X*Y2 + A3*X3*Y + A4*X2*Y2 + A5*X*Y3
C
      RV2 = B1+B2+B3
C
C IMAGINEARY  VECTOR DIFFUSENESS
C
C STRAIGHT X-DEPENDENCE
      A0=P(7,1)
      A1=P(7,2)
      A2=P(7,3)
      A3=P(7,4)
      A4=P(7,5)
C      A5=P(1,6)
C      A6=P(1,7)
C
      B1 = A0 + A1*X + A2*X2 + A3*X3 + A4*X4 
C STRAIGHT Y DEPENDENCE
      A0=0.0D0
      A1=P(7,6)
      A2=P(7,7)
      A3=P(7,8)
      A4=P(7,9)
      A5=P(7,10)
C
      B2 = A0 + A1*Y + A2*Y2 + A3*Y3 + A4*Y4 + A5*Y5
CROSS TERMS
      A0=P(7,11)
      A1=P(7,12)
      A2=P(7,13)
      A3=P(7,14)
      A4=P(7,15)
      A5=P(7,16)
C
      B3 = A0*X*Y + A1*X2*Y + A2*X*Y2 + A3*X3*Y + A4*X2*Y2 + A5*X*Y3
C
      AV2 = 0.7D0*(B1+B2+B3)
C C
C IMAGINARY VECTOR PARABOLIC FERMI PARAMETER
C
C STRAIGHT X-DEPENDENCE
      A0=P(8,1)
      A1=P(8,2)
      A2=P(8,3)
      A3=P(8,4)
      A4=P(8,5)
C      A5=P(1,6)
C      A6=P(1,7)
C
      B1 = A0 + A1*X + A2*X2 + A3*X3 + A4*X4 
C STRAIGHT Y DEPENDENCE
      A0=0.0D0
      A1=P(8,6)
      A2=P(8,7)
      A3=P(8,8)
      A4=P(8,9)
      A5=P(8,10)
C
      B2 = A0 + A1*Y + A2*Y2 + A3*Y3 + A4*Y4 + A5*Y5
CROSS TERMS
      A0=P(8,11)
      A1=P(8,12)
      A2=P(8,13)
      A3=P(8,14)
      A4=P(8,15)
      A5=P(8,16)
C
      B3 = A0*X*Y + A1*X2*Y + A2*X*Y2 + A3*X3*Y + A4*X2*Y2 + A5*X*Y3
C
      WGV2 =(B1+B2+B3)
C
CC DIFF OF REAL STRENGHTS
C
C STRAIGHT X-DEPENDENCE
C STRAIGHT X-DEPENDENCE
      A0=P(9,1)
      A1=P(9,2)
      A2=P(9,3)
      A3=P(9,4)
      A4=P(9,5)
C      A5=P(1,6)
C      A6=P(1,7)
C
      B1 = A0 + A1*X + A2*X2 + A3*X3 + A4*X4 
C STRAIGHT Y DEPENDENCE
      A0=0.0D0
      A1=P(9,6)
      A2=P(9,7)
      A3=P(9,8)
      A4=P(9,9)
      A5=P(9,10)
C
      B2 = A0 + A1*Y + A2*Y2 + A3*Y3 + A4*Y4 + A5*Y5
CROSS TERMS
      A0=P(9,11)
      A1=P(9,12)
      A2=P(9,13)
      A3=P(9,14)
      A4=P(9,15)
      A5=P(9,16)
C
      B3 = A0*X*Y + A1*X2*Y + A2*X*Y2 + A3*X3*Y + A4*X2*Y2 + A5*X*Y3
C
      DIFFR = 700.D0*(B1+B2+B3)
C
C REAL SCALAR RADIUS
C
C STRAIGHT X-DEPENDENCE
C STRAIGHT X-DEPENDENCE
      A0=P(10,1)
      A1=P(10,2)
      A2=P(10,3)
      A3=P(10,4)
      A4=P(10,5)
C      A5=P(1,6)
C      A6=P(1,7)
C
      B1 = A0 + A1*X + A2*X2 + A3*X3 + A4*X4 
C STRAIGHT Y DEPENDENCE
      A0=0.0D0
      A1=P(10,6)
      A2=P(10,7)
      A3=P(10,8)
      A4=P(10,9)
      A5=P(10,10)
C
      B2 = A0 + A1*Y + A2*Y2 + A3*Y3 + A4*Y4 + A5*Y5
CROSS TERMS
      A0=P(10,11)
      A1=P(10,12)
      A2=P(10,13)
      A3=P(10,14)
      A4=P(10,15)
      A5=P(10,16)
C
      B3 = A0*X*Y + A1*X2*Y + A2*X*Y2 + A3*X3*Y + A4*X2*Y2 + A5*X*Y3
C
      RS1 = B1+B2+B3
C
C REAL SCALAR DIFFUSENESS
C
C STRAIGHT X-DEPENDENCE
C STRAIGHT X-DEPENDENCE
      A0=P(11,1)
      A1=P(11,2)
      A2=P(11,3)
      A3=P(11,4)
      A4=P(11,5)
C      A5=P(1,6)
C      A6=P(1,7)
C
      B1 = A0 + A1*X + A2*X2 + A3*X3 + A4*X4 
C STRAIGHT Y DEPENDENCE
      A0=0.0D0
      A1=P(11,6)
      A2=P(11,7)
      A3=P(11,8)
      A4=P(11,9)
      A5=P(11,10)
C
      B2 = A0 + A1*Y + A2*Y2 + A3*Y3 + A4*Y4 + A5*Y5
CROSS TERMS
      A0=P(11,11)
      A1=P(11,12)
      A2=P(11,13)
      A3=P(11,14)
      A4=P(11,15)
      A5=P(11,16)
C
      B3 = A0*X*Y + A1*X2*Y + A2*X*Y2 + A3*X3*Y + A4*X2*Y2 + A5*X*Y3
C
      AS1 = 0.7D0*(B1+B2+B3)
C
C
C REAL SCALAR PARABOLIC FERMI PARAMETER
C
C STRAIGHT X-DEPENDENCE
C STRAIGHT X-DEPENDENCE
      A0=P(12,1)
      A1=P(12,2)
      A2=P(12,3)
      A3=P(12,4)
      A4=P(12,5)
C      A5=P(1,6)
C      A6=P(1,7)
C
      B1 = A0 + A1*X + A2*X2 + A3*X3 + A4*X4 
C STRAIGHT Y DEPENDENCE
      A0=0.0D0
      A1=P(12,6)
      A2=P(12,7)
      A3=P(12,8)
      A4=P(12,9)
      A5=P(12,10)
C
      B2 = A0 + A1*Y + A2*Y2 + A3*Y3 + A4*Y4 + A5*Y5
CROSS TERMS
      A0=P(12,11)
      A1=P(12,12)
      A2=P(12,13)
      A3=P(12,14)
      A4=P(12,15)
      A5=P(12,16)
C
      B3 = A0*X*Y + A1*X2*Y + A2*X*Y2 + A3*X3*Y + A4*X2*Y2 + A5*X*Y3
C
      WGS1 =(B1+B2+B3)
C
CC DIFF OF IMGAGINARY STRENGHTS
C
C STRAIGHT X-DEPENDENCE
C STRAIGHT X-DEPENDENCE
      A0=P(13,1)
      A1=P(13,2)
      A2=P(13,3)
      A3=P(13,4)
      A4=P(13,5)
C      A5=P(1,6)
C      A6=P(1,7)
C
      B1 = A0 + A1*X + A2*X2 + A3*X3 + A4*X4 
C STRAIGHT Y DEPENDENCE
      A0=0.0D0
      A1=P(13,6)
      A2=P(13,7)
      A3=P(13,8)
      A4=P(13,9)
      A5=P(13,10)
C
      B2 = A0 + A1*Y + A2*Y2 + A3*Y3 + A4*Y4 + A5*Y5
CROSS TERMS
      A0=P(13,11)
      A1=P(13,12)
      A2=P(13,13)
      A3=P(13,14)
      A4=P(13,15)
      A5=P(13,16)
C
      B3 = A0*X*Y + A1*X2*Y + A2*X*Y2 + A3*X3*Y + A4*X2*Y2 + A5*X*Y3
C
C      TSCALE1=P(32,2)*200
      DIFFI=-160*(B1+B2+B3)
C
C IMAG SCALAR RADIUS
C
C STRAIGHT X-DEPENDENCE
C STRAIGHT X-DEPENDENCE
      A0=P(14,1)
      A1=P(14,2)
      A2=P(14,3)
      A3=P(14,4)
      A4=P(14,5)
C      A5=P(1,6)
C      A6=P(1,7)
C
      B1 = A0 + A1*X + A2*X2 + A3*X3 + A4*X4 
C STRAIGHT Y DEPENDENCE
      A0=0.0D0
      A1=P(14,6)
      A2=P(14,7)
      A3=P(14,8)
      A4=P(14,9)
      A5=P(14,10)
C
      B2 = A0 + A1*Y + A2*Y2 + A3*Y3 + A4*Y4 + A5*Y5
CROSS TERMS
      A0=P(14,11)
      A1=P(14,12)
      A2=P(14,13)
      A3=P(14,14)
      A4=P(14,15)
      A5=P(14,16)
C
      B3 = A0*X*Y + A1*X2*Y + A2*X*Y2 + A3*X3*Y + A4*X2*Y2 + A5*X*Y3
C
      RS2 = B1+B2+B3
C
C IMAG SCALAR DIFFUSENESS
C
C STRAIGHT X-DEPENDENCE
C STRAIGHT X-DEPENDENCE
      A0=P(15,1)
      A1=P(15,2)
      A2=P(15,3)
      A3=P(15,4)
      A4=P(15,5)
C      A5=P(1,6)
C      A6=P(1,7)
C
      B1 = A0 + A1*X + A2*X2 + A3*X3 + A4*X4 
C STRAIGHT Y DEPENDENCE
      A0=0.0D0
      A1=P(15,6)
      A2=P(15,7)
      A3=P(15,8)
      A4=P(15,9)
      A5=P(15,10)
C
      B2 = A0 + A1*Y + A2*Y2 + A3*Y3 + A4*Y4 + A5*Y5
CROSS TERMS
      A0=P(15,11)
      A1=P(15,12)
      A2=P(15,13)
      A3=P(15,14)
      A4=P(15,15)
      A5=P(15,16)
C
      B3 = A0*X*Y + A1*X2*Y + A2*X*Y2 + A3*X3*Y + A4*X2*Y2 + A5*X*Y3
      AS2 = 0.7D0*(B1+B2+B3)
C
C
C IMAGINARY SCALAR PARABOLIC FERMI PARAMETER
C
C STRAIGHT X-DEPENDENCE
C STRAIGHT X-DEPENDENCE
      A0=P(16,1)
      A1=P(16,2)
      A2=P(16,3)
      A3=P(16,4)
      A4=P(16,5)
C      A5=P(1,6)
C      A6=P(1,7)
C
      B1 = A0 + A1*X + A2*X2 + A3*X3 + A4*X4 
C STRAIGHT Y DEPENDENCE
      A0=0.0D0
      A1=P(16,6)
      A2=P(16,7)
      A3=P(16,8)
      A4=P(16,9)
      A5=P(16,10)
C
      B2 = A0 + A1*Y + A2*Y2 + A3*Y3 + A4*Y4 + A5*Y5
CROSS TERMS
      A0=P(16,11)
      A1=P(16,12)
      A2=P(16,13)
      A3=P(16,14)
      A4=P(16,15)
      A5=P(16,16)
C
      B3 = A0*X*Y + A1*X2*Y + A2*X*Y2 + A3*X3*Y + A4*X2*Y2 + A5*X*Y3
C
      WGS2 =(B1+B2+B3)
C
C
C NOW THE REAL PEAKING SUM OF STRENTHS
C
C
C STRAIGHT X-DEPENDENCE
C STRAIGHT X-DEPENDENCE
      A0=P(17,1)
      A1=P(17,2)
      A2=P(17,3)
      A3=P(17,4)
      A4=P(17,5)
C      A5=P(1,6)
C      A6=P(1,7)
C
      B1 = A0 + A1*X + A2*X2 + A3*X3 + A4*X4 
C STRAIGHT Y DEPENDENCE
      A0=0.0D0
      A1=P(17,6)
      A2=P(17,7)
      A3=P(17,8)
      A4=P(17,9)
      A5=P(17,10)
C
      B2 = A0 + A1*Y + A2*Y2 + A3*Y3 + A4*Y4 + A5*Y5
CROSS TERMS
      A0=P(17,11)
      A1=P(17,12)
      A2=P(17,13)
      A3=P(17,14)
      A4=P(17,15)
      A5=P(17,16)
C
      B3 = A0*X*Y + A1*X2*Y + A2*X*Y2 + A3*X3*Y + A4*X2*Y2 + A5*X*Y3
C
      SUMRIV = 16.0D0*(B1+B2+B3)
      
C NOW THE SURFACE PEAKING DIFFERENCE OF REAL STRENTHS
C
C STRAIGHT X-DEPENDENCE
      A0=P(18,1)
      A1=P(18,2)
      A2=P(18,3)
      A3=P(18,4)
      A4=P(18,5)
C      A5=P(1,6)
C      A6=P(1,7)
C
      B1 = A0 + A1*X + A2*X2 + A3*X3 + A4*X4 
C STRAIGHT Y DEPENDENCE
      A0=0.0D0
      A1=P(18,6)
      A2=P(18,7)
      A3=P(18,8)
      A4=P(18,9)
      A5=P(18,10)
C
      B2 = A0 + A1*Y + A2*Y2 + A3*Y3 + A4*Y4 + A5*Y5
CROSS TERMS
      A0=P(18,11)
      A1=P(18,12)
      A2=P(18,13)
      A3=P(18,14)
      A4=P(18,15)
      A5=P(18,16)
C
      B3 = A0*X*Y + A1*X2*Y + A2*X*Y2 + A3*X3*Y + A4*X2*Y2 + A5*X*Y3
C
      DIFRIV = -80.0D0*(B1+B2+B3)
C NOW THE IMAGINARY SURFACE PEAKING SUM OF STRENTHS
C
C STRAIGHT X-DEPENDENCE
      A0=P(19,1)
      A1=P(19,2)
      A2=P(19,3)
      A3=P(19,4)
      A4=P(19,5)
C      A5=P(1,6)
C      A6=P(1,7)
C
      B1 = A0 + A1*X + A2*X2 + A3*X3 + A4*X4 
C STRAIGHT Y DEPENDENCE
      A0=0.0D0
      A1=P(19,6)
      A2=P(19,7)
      A3=P(19,8)
      A4=P(19,9)
      A5=P(19,10)
C
      B2 = A0 + A1*Y + A2*Y2 + A3*Y3 + A4*Y4 + A5*Y5
CROSS TERMS
      A0=P(19,11)
      A1=P(19,12)
      A2=P(19,13)
      A3=P(19,14)
      A4=P(19,15)
      A5=P(19,16)
C
      B3 = A0*X*Y + A1*X2*Y + A2*X*Y2 + A3*X3*Y + A4*X2*Y2 + A5*X*Y3
C
      SUMIIV = 16.0D0*(B1+B2+B3)
      
C NOW THE SURFACE PEAKING IMAGINARY DIFFERENCE OF STRENTHS
C
C STRAIGHT X-DEPENDENCE
C STRAIGHT X-DEPENDENCE
      A0=P(20,1)
      A1=P(20,2)
      A2=P(20,3)
      A3=P(20,4)
      A4=P(20,5)
C      A5=P(1,6)
C      A6=P(1,7)
C
      B1 = A0 + A1*X + A2*X2 + A3*X3 + A4*X4 
C STRAIGHT Y DEPENDENCE
      A0=0.0D0
      A1=P(20,6)
      A2=P(20,7)
      A3=P(20,8)
      A4=P(20,9)
      A5=P(20,10)
C
      B2 = A0 + A1*Y + A2*Y2 + A3*Y3 + A4*Y4 + A5*Y5
CROSS TERMS
      A0=P(20,11)
      A1=P(20,12)
      A2=P(20,13)
      A3=P(20,14)
      A4=P(20,15)
      A5=P(20,16)
C
      B3 = A0*X*Y + A1*X2*Y + A2*X*Y2 + A3*X3*Y + A4*X2*Y2 + A5*X*Y3
C
      DIFIIV = -80.0D0*(B1+B2+B3)
C
CC-----------------------------------------------
C
      VV=0.5D0*(SUMR+DIFFR)
      VS=0.5D0*(SUMR-DIFFR)
      WV=0.5D0*(SUMI+DIFFI)
      WS=0.5D0*(SUMI-DIFFI)
C PLEASE EXCUSE USE OF OLD 'IV' LABELS!
      WVSP=0.5D0*(SUMIIV+DIFIIV)
      WSSP=0.5D0*(SUMIIV-DIFIIV)
      VVSP=0.5D0*(SUMRIV+DIFRIV)
      VSSP=0.5D0*(SUMRIV-DIFRIV)
C
C      VV = VV + VVIV*(FNN-FZZ)/FAA
C      VS = VS + VSIV*(FNN-FZZ)/FAA
C      WV = WV + WVIV*(FNN-FZZ)/FAA
C      WS = WS + WSIV*(FNN-FZZ)/FAA      
C
      RV1SP=  RV1
      AV1SP=  AV1
      RS1SP=  RS1
      AS1SP=  AS1
C
      RV22=RV2
      AV22=AV2
      RS22=RS2
      AS22=AS2
C TRY TO PREVENT DIFFUSENESS GETTING TOO SMALL
C
      AV1 = DSQRT(AV1**2 + 0.04D0)
      AV2 = DSQRT(AV2**2 + 0.04D0)
      AS1 = DSQRT(AS1**2 + 0.04D0)
      AS2 = DSQRT(AS2**2 + 0.04D0)
C
C TRY TO PREVENT RADII GETTING TOO SMALL
C
      RV1 = DSQRT(RV1**2 + 0.4D0**2)
      RV2 = DSQRT(RV2**2 + 0.4D0**2)
      RS1 = DSQRT(RS1**2 + 0.4D0**2)
      RS2 = DSQRT(RS2**2 + 0.4D0**2)
C
      RV22=RV2
      AV22=AV2
      RS22=RS2
      AS22=AS2      
C
C TRY TO PREVENT WEIRD POTENTIALS, MAKE SURE WGV1,S1 NEVER GO BELOW -0.2
C
      PI = DACOS(-1.0D0)
      WGV1= .3D0+(2.0D0/PI)*datan(WGV1)/2.0D0
      WGS1= .3D0+(2.0D0/PI)*datan(WGS1)/2.0D0
C CONSTRAIN THE IMAGINARY W PARAMTERS TO BE CLOSE.
C
C      WGP2= (WGV2+WGS2)/2.0D0
C      WGM2= (WGV2-WGS2)/2.0D0
C      WGM2=0.2D0*DATAN(WGM2)*(2.0D0/3.14159D0)
C      WGV2=(WGP2+WGM2)/2.0D0
C      WGS2=(WGP2-WGM2)/2.0D0
C---- TRY FIXING WGS2 TO BE THE SAME AS WGS1...
      WGS2=WGV2      
C      
CC
C FIX UP TRANSMITTED BACK VARIABLES TO BE SAME AS SHINICHI USED.
C
C      RS12=RS1SP
C      RV12=RV1SP
C      AS12=AS1SP
C      AV12=AV1SP
C      VVSP=0.0D0
C      VSSP=0.0D0
C      WVSP=0.0D0
C      WSSP=0.0D0
C      VV2=VVSP
C      VS2=VSSP
      WV2=WVSP
      WS2=WSSP
C      
C      WVV=0.0D0
C      WWV=0.0D0
C      WVS=0.0D0
C      WWS=0.0D0
C
C
C
CC
      DSMALL=1.D-160
      DO 10 I=1,NSTEP
         TMP1 = R(I)
C RE VEC
         WG=1.0D0+WGV1*(TMP1/(RV1*ACB))**2
            S1 = DSECH(RV1*ACB/AV1)
            S2 = DSECH(TMP1/AV1)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
            RVA1=RECV*VV*WG*B/(A+B)
            S1 = DSECH(RV1SP*ACB/AV1SP)
            S2 = DSECH(TMP1/AV1SP)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECV*VVSP*A*B/(A+B)**2
         RVA1=RVA1+SURF
C IM VEC          
         WG=1.0D0+WGV2*(TMP1/(RV2*ACB))**2
            S1 = DSECH(RV2*ACB/AV2)
            S2 = DSECH(TMP1/AV2)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
            RVA2=RECV*WV*WG*B/(A+B)
            S1 = DSECH(RV22*ACB/AV22)
            S2 = DSECH(TMP1/AV22)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECV*WV2*A*B/(A+B)**2
         RVA2=RVA2+SURF
C RE SCA         
         WG=1.0D0+WGS1*(TMP1/(RS1*ACB))**2
            S1 = DSECH(RS1*ACB/AS1)
            S2 = DSECH(TMP1/AS1)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RSA1=RECS*VS*WG*B/(A+B)
            S1 = DSECH(RS1SP*ACB/AS1SP)
            S2 = DSECH(TMP1/AS1SP)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECS*VSSP*A*B/(A+B)**2
         RSA1=RSA1+SURF
C IM SCA
         WG=1.0D0+WGS2*(TMP1/(RS2*ACB))**2
            S1 = DSECH(RS2*ACB/AS2)
            S2 = DSECH(TMP1/AS2)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
          RSA2=RECS*WS*WG*B/(A+B)
            S1 = DSECH(RS22*ACB/AS22)
            S2 = DSECH(TMP1/AS22)
            PROD=S1*S2
            A=S1-PROD + DSMALL
            B=S2-PROD + DSMALL
         SURF=RECS*WS2*A*B/(A+B)**2
         RSA2=RSA2+SURF
CC
         VVA(I)=RVA1
         WVA(I)=RVA2
         VSA(I)=RSA1
         WSA(I)=RSA2
 10    CONTINUE
C      DO 20 I=(NSTEP+2),1000
C        VVA(I)=0.0D0
C        WVA(I)=0.0D0
C        VSA(I)=0.0D0
C        WSA(I)=0.0D0
C 20   CONTINUE
CC
      RETURN
      END
