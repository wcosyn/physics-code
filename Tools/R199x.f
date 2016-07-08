C ---------------------------------------------------------------------
      REAL*8 Function R1990(xx,QQ) 
C ---------------------------------------------------------------------
C L. W. Withlow et al., Phys. Lett. B250 194(1990)
      REAL*8 xx,QQ
      REAL*8 Theta,LOGprt,QQthr 
      REAL*8 a1,a2,a3,R1, b1,b2,b3,R2, c1,c2,c3,R3

      DATA a1,a2,a3 /0.0672D0,0.4671D0,1.8979D0/
      DATA b1,b2,b3 /0.0635D0,0.5747D0,-0.3534D0/
      DATA c1,c2,c3 /0.0599D0,0.5088D0,2.1081D0/

      Theta=1.D0+12.D0*QQ/(QQ+1.D0)*0.125D0**2
     -      /(0.125D0**2+xx*xx) 

      LOGprt=Dlog(QQ/0.04D0)
      QQthr=5.D0*(1.D0-xx)**5 

      R1=a1*Theta/LOGprt+a2/(QQ**4+a3**4)**0.25D0
      R2=b1*Theta/LOGprt+b2/QQ+b3/(QQ*QQ+0.3D0**2)
      R3=c1*Theta/LOGprt+c2*((QQ-QQthr)**2+c3*c3)**(-0.5D0)

      R1990=(R1+R2+R3)/3.D0

      RETURN
      END
C ---------------------------------------------------------------------
      REAL*8 Function R1998(xx,QQ) 
C ---------------------------------------------------------------------
C  (E143) K. Abe et al., Phys. Lett B452 194 (1999)
      IMPLICIT REAL*8(A-H,O-Z)

      DATA a1,a2,a3,a4,a5,a6 / 0.0485D0,  0.5470D0,  2.0621D0
     >                      , -0.3804D0,  0.5090D0, -0.0285/
      DATA b1,b2,b3,b4,b5,b6 / 0.0481D0,  0.6114D0, -0.3509D0,
     >                        -0.4611D0,  0.7172D0, -0.0317/
      DATA c1,c2,c3,c4,c5,c6 / 0.0577D0,  0.4644D0,  1.8288D0,
     >                        12.3708D0,-43.1043D0, 41.7415/


      Theta=1.D0+12.D0*QQ*0.125D0*0.125D0
     -      /((QQ+1.D0)*(0.125D0*0.125D0+xx*xx)) 

      RLOGprt=Dlog(QQ/0.04D0)
      QQthr=c4*xx+c5*xx**2.D0+c6*xx**3

      R1=a1*Theta/RLOGprt
     > + a2*(1.D0+a4*xx+a5*xx*xx)*(xx**a6)/((QQ**4+a3**4)**0.25D0)

      R2=b1*Theta/RLOGprt
     > +(b2/QQ+b3/(QQ**2.D0+0.3D0*0.3D0))*(1.D0+b4*xx+b5*xx*xx)*(xx**b6)

      R3=c1*Theta/RLOGprt +c2/DSQRT((QQ-QQthr)**2+c3**2)

      R1998=(R1+R2+R3)/3.D0

      RETURN
      END
C ---------------------------------------------------------------------
