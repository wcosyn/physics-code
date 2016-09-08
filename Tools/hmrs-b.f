c------------------------------------------------------------------
c added by S. Kumano for F1,   2015-08-04, HMRS
      real*8 function f1pn_hmrs(x,q2)
      implicit real*8(a-h,o-z)
      dmp=0.938272046d0
      c1=1.d0/(2.d0*x)
      c2=4.d0*dmp*dmp*x*x/q2
      R_1990=R1990(x,q2) 
      f1pn_hmrs=c1/(1.d0+R_1990)*(1.d0+c2)*f2pn_hmrs(x,q2)
      return
      end
c---------------------------
      real*8 function f1p_hmrs(x,q2)
      implicit real*8(a-h,o-z)
      dmp=0.938272046d0
      c1=1.d0/(2.d0*x)
      c2=4.d0*dmp*dmp*x*x/q2
      R_1990=R1990(x,q2) 
      f1p_hmrs=c1/(1.d0+R_1990)*(1.d0+c2)*f2p_hmrs(x,q2)
      return
      end
c--------------------------------------      
      real*8 function f1n_hmrs(x,q2)
      implicit real*8(a-h,o-z)
      dmp=0.938272046d0
      c1=1.d0/(2.d0*x)
      c2=4.d0*dmp*dmp*x*x/q2
      R_1990=R1990(x,q2) 
      f1n_hmrs=c1/(1.d0+R_1990)*(1.d0+c2)*f2n_hmrs(x,q2)
      return
      end
! c------------------------------------------------------------------
! c------------------------------------------------------------------
! c The following subroutine is from AAC polarized PDF fit code.
! C ---------------------------------------------------------------------
!       REAL*8 Function R1990(xx,QQ) 
! C ---------------------------------------------------------------------
! C L. W. Withlow et al., Phys. Lett. B250 194(1990)
!       REAL*8 xx,QQ
!       REAL*8 Theta,LOGprt,QQthr 
!       REAL*8 a1,a2,a3,R1, b1,b2,b3,R2, c1,c2,c3,R3
! 
!       DATA a1,a2,a3 /0.0672D0,0.4671D0,1.8979D0/
!       DATA b1,b2,b3 /0.0635D0,0.5747D0,-0.3534D0/
!       DATA c1,c2,c3 /0.0599D0,0.5088D0,2.1081D0/
! 
!       Theta=1.D0+12.D0*QQ/(QQ+1.D0)*0.125D0**2
!      -      /(0.125D0**2+xx*xx) 
! 
!       LOGprt=Dlog(QQ/0.04D0)
!       QQthr=5.D0*(1.D0-xx)**5 
! 
!       R1=a1*Theta/LOGprt+a2/(QQ**4+a3**4)**0.25D0
!       R2=b1*Theta/LOGprt+b2/QQ+b3/(QQ*QQ+0.3D0**2)
!       R3=c1*Theta/LOGprt+c2*((QQ-QQthr)**2+c3*c3)**(-0.5D0)
! 
!       R1990=(R1+R2+R3)/3.D0
! 
!       RETURN
!       END
! C ---------------------------------------------------------------------


c=====================================================================
c---------------------------------------------------------------
c (proton+neutron)/2 sturcture function.
      real*8 function f2pn_hmrs(x,q2)
      implicit real*8(a-h,o-z)
      if(x.ge.1.d0)go to 5510
      f2pn_hmrs=(f2p_hmrs(x,q2)+f2n_hmrs(x,q2))/2.d0
      go to 5511
 5510 f2pn_hmrs=0.d0
 5511 return
      end
c------------------------------------------------------
c proton structure function F2p_hmrs(x,q2).
      real*8 function f2p_hmrs(x,q2)
      implicit real*8(a-h,o-z)
      f2p_hmrs=1.0d0/9.d0*(4.d0*xuplus_hmrs(x,q2)+xdplus_hmrs(x,q2)
     1                +xsplus_hmrs(x,q2)+4.d0*xcplus_hmrs(x,q2))
      return
      end
c------------------------------------------------------
c neutron structure function F2n_hmrs(x,q2).
      real*8 function f2n_hmrs(x,q2)
      implicit real*8(a-h,o-z)
      f2n_hmrs=1.0d0/9.d0*(4.d0*xdplus_hmrs(x,q2)+xuplus_hmrs(x,q2)
     1                +xsplus_hmrs(x,q2)+4.d0*xcplus_hmrs(x,q2))
      return
      end
C ---------------------------------------------------------------------
c modified on Jan. 21, 2002.
c===================================================================
c u+ubar quark distribution.
      real*8 function xuplus_hmrs(x,q2)
      implicit real*8(a-h,o-z)
      xuplus_hmrs=xuv_hmrs_ana(x,q2)+2.d0*xubar_hmrs_ana(x,q2)
      return
      end
c-----------------------------------------------------------
c d+dbar quark distribution.
      real*8 function xdplus_hmrs(x,q2)
      implicit real*8(a-h,o-z)
      xdplus_hmrs=xdv_hmrs_ana(x,q2)+2.d0*xdbar_hmrs_ana(x,q2)
      return
      end
c-----------------------------------------------------------
c s+sbar quark distribution.
      real*8 function xsplus_hmrs(x,q2)
      implicit real*8(a-h,o-z)
      xsplus_hmrs=2.d0*xsbar_hmrs_ana(x,q2)
      return
      end
c-----------------------------------------------------------
c c+cbar quark distribution.
      real*8 function xcplus_hmrs(x,q2)
      implicit real*8(a-h,o-z)
      xcplus_hmrs=2.d0*xcbar_hmrs_ana(x,q2)
      return
      end
c-----------------------------------------------------------

c=====================================================================
c The following part is taken from bfv2t.for on 2015-08-07
c HMRS-B analytic form at Q2=4 GeV2.
c=====================================================================
C ---------------------------------------------------------------------
      REAL*8 FUNCTION QQ(X,A,B,C)
      IMPLICIT REAL*8(A-H,O-Z)
      QQ=A*(X**B)*((1.0D0-X)**C)
      RETURN
      END
c=====================================================================
c uv quark distribution.
      real*8 function xuv_hmrs_ana(x,q2)
      implicit real*8(a-h,o-z)
C  HMRS-B U+UBAR, D+DBAR OR S+SBAR DISTRIBUTION
C  AT Q**2=4 GeV**2 (LAMBDA=0.19 GeV). 
        A1=0.5469D0
        B1=0.237D0
        C1=4.07D0
        A2=A1*23.8D0
        B2=B1+1.D0
        C2=C1
      QV=QQ(X,A1,B1,C1)+QQ(X,A2,B2,C2)
        AA1=0.6957D0
        BB1=0.426D0
        CC1=4.82D0
        AA2=AA1*6.32D0
        BB2=BB1+1.D0
        CC2=CC1
      DV=QQ(X,AA1,BB1,CC1)+QQ(X,AA2,BB2,CC2)
      UV=QV-DV
      xuv_hmrs_ana=UV
      return
      end
c-----------------------------------------------------------
c dv quark distribution.
      real*8 function xdv_hmrs_ana(x,q2)
      implicit real*8(a-h,o-z)
C  HMRS-B U+UBAR, D+DBAR OR S+SBAR DISTRIBUTION
C  AT Q**2=4 GeV**2 (LAMBDA=0.19 GeV). 
        AA1=0.6957D0
        BB1=0.426D0
        CC1=4.82D0
        AA2=AA1*6.32D0
        BB2=BB1+1.D0
        CC2=CC1
      DV=QQ(X,AA1,BB1,CC1)+QQ(X,AA2,BB2,CC2)
      xdv_hmrs_ana=DV
      return
      end
c-----------------------------------------------------------
c ubar distribution.
      real*8 function xubar_hmrs_ana(x,q2)
      implicit real*8(a-h,o-z)
C  HMRS-B U+UBAR, D+DBAR OR S+SBAR DISTRIBUTION
C  AT Q**2=4 GeV**2 (LAMBDA=0.19 GeV). 
        AAA1=5.25D0
        BBB1=0.401D0
        CCC1=9.75D0
      QSEA0=QQ(X,AAA1,BBB1,CCC1)        
      qsea=QSEA0/5.d0
      xubar_hmrs_ana=qsea
      return
      end
c-----------------------------------------------------------
c dbar distribution.
      real*8 function xdbar_hmrs_ana(x,q2)
      implicit real*8(a-h,o-z)
C  HMRS-B U+UBAR, D+DBAR OR S+SBAR DISTRIBUTION
C  AT Q**2=4 GeV**2 (LAMBDA=0.19 GeV). 
        AAA1=5.25D0
        BBB1=0.401D0
        CCC1=9.75D0
      QSEA0=QQ(X,AAA1,BBB1,CCC1)        
      qsea=QSEA0/5.d0
      xdbar_hmrs_ana=qsea
      return
      end
c-----------------------------------------------------------
c sbar distribution.
      real*8 function xsbar_hmrs_ana(x,q2)
      implicit real*8(a-h,o-z)
C  HMRS-B U+UBAR, D+DBAR OR S+SBAR DISTRIBUTION
C  AT Q**2=4 GeV**2 (LAMBDA=0.19 GeV). 
        AAA1=5.25D0
        BBB1=0.401D0
        CCC1=9.75D0
      QSEA0=QQ(X,AAA1,BBB1,CCC1)        
      ssea=QSEA0/10.d0
      xsbar_hmrs_ana=ssea
      return
      end
c-----------------------------------------------------------
c cbar distribution.
      real*8 function xcbar_hmrs_ana(x,q2)
      implicit real*8(a-h,o-z)
      chm=0.d0
      xcbar_hmrs_ana=chm
      return
      end
c-----------------------------------------------------------
c gluon distribution.
      real*8 function xg_hmrs_ana(x,q2)
      implicit real*8(a-h,o-z)
C  HMRS-B  GLUON DISTRIBUTION
C  AT Q**2=4 GeV**2 (LAMBDA=0.19 GeV).
        A1=2.855D0
        B1=0.D0
        C1=5.1D0
      G0=QQ(X,A1,B1,C1)
      xg_hmrs_ana=G0
      return
      end
c-----------------------------------------------------------------


