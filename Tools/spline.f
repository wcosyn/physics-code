      MODULE spline
       PRIVATE

       PUBLIC cubbasis,cubherm 

       CONTAINS


        SUBROUTINE cubbasis(xold,n,xnew,m,spl)

C *** Die Routine bereitet Basissplines mit Hilfe der SPLMOD und SPLFAK vor.
C *** Die Parameter sind wie bei cubherm definiert
C *** Die Dimension ist jetzt anders : spl(m,n) statt spl(m,4)
C *** Der Parameter index ist nicht noetig !

         IMPLICIT NONE
         INTEGER n,m

         REAL xold(n),xnew(m),spl(m,n)
         REAL pspl(n)
         INTEGER i,j

         CALL splfak(xold,n)

         DO i=1,m
          IF((xnew(i).GE.xold(1)).AND.(xnew(i).LE.xold(n))) THEN
            CALL splmod(xold,n,xnew(i),pspl)
            DO j=1,n
C***  hier kann Interpolation aendern
             spl(i,j)=pspl(j)
     X            *xold(j)/xnew(i)
            END DO
           ELSE
            DO j=1,n
             spl(i,j)=0.0
            END DO
            IF(xnew(i).LE.xold(1)) THEN
              spl(i,1)=(xold(2)-xnew(i))/(xold(2)-xold(1))
     X                      *xold(1)/xnew(i)
              spl(i,2)=(xnew(i)-xold(1))/(xold(2)-xold(1))
     X                      *xold(2)/xnew(i)
            END IF
          END IF   
         END DO

        END SUBROUTINE cubbasis




      subroutine cubherm(xold,n,xnew,m,spl,index)
c
c     This subroutine prepares the interpolation of a function
c     given at the n grid points xold to the m points xnew.
c     The interpolating functions are cubic hermitian splines.
c     The first derivatives at the grid points are taken from a
c     parabola through the actual grid point and its two neighbours.
c     For the end points the parabola is taken through the two
c     right or left neighbours, resp.
c
c     In the calling routine you still have to take the sum i=1,4
c     ynew(j)=sum_i spl(j,i)*yold(index(j,i))
c
c     by Dirk Hueber, 08.02.1996
c
cc    parameter (nmax=100,mmax=100)
      dimension xold(n),xnew(m),spl(m,4),index(m,4)
      logical enough
      enough=n.ge.3
cc      if (n.gt.nmax) stop'nmax to small in cubherm'
cc      if (m.gt.mmax) stop'mmax to small in cubherm'
c
c     evaluation of the indices
c
      do 10 j=1,m
c
c     If xnew(j) is less than xold(1), we extrapolate using the first
c     two grid points xold(1) and xold(2).
c
       index(j,2)=1
10    continue
      do 11 i=1,n
       do 11 j=1,m
        if (xnew(j).gt.xold(i)) index(j,2)=i
11    continue
      do 12 j=1,m
       index(j,2)=min(index(j,2),n-1)
       index(j,1)=index(j,2)-1
       index(j,3)=index(j,2)+1
       index(j,4)=index(j,2)+2
c
c     Indices 1 and 4 are used only for the ewertcation of the derivatives.
c     The following settings provide the derivatives at the borders of the
c     grid xold.
c
       if (index(j,1).eq.0) index(j,1)=3
       if (index(j,4).eq.n+1) index(j,4)=n-2
12    continue
      do 20 j=1,m
c
c     We don't extrapolate to the right!
c
       if (xnew(j).le.xold(n).and.enough) then
        i0=index(j,1)
        i1=index(j,2)
        i2=index(j,3)
        i3=index(j,4)
        x0=xold(i0)
        x1=xold(i1)
        x2=xold(i2)
        x3=xold(i3)
c
c      Factors for the derivatives
c
        d10=x1-x0
        d21=x2-x1
        d32=x3-x2
        d20=x2-x0
        d31=x3-x1
        dfak13=(d21/d10-d10/d21)/d20
        dfak14=-d32/(d21*d31)
        dfak23=d10/(d21*d20)
        dfak24=(d32/d21-d21/d32)/d31
        dfak03=-d21/(d10*d20)
        dfak34=d21/(d32*d31)
c
c     the cubic hermitian splines
c
        xn=xnew(j)
        dn1=xn-x1
        d2n=x2-xn
        phidiv=1./(d21*d21*d21)
        phi1=d2n*d2n*phidiv*(d21+2.*dn1)
        phi2=dn1*dn1*phidiv*(d21+2.*d2n)
        phidiv=phidiv*d21*dn1*d2n
        phi3=phidiv*d2n
        phi4=-phidiv*dn1
c
c     combining everything to the final factors
c
        spl(j,2)=phi1+phi3*dfak13+phi4*dfak14
        spl(j,3)=phi2+phi3*dfak23+phi4*dfak24
        spl(j,1)=phi3*dfak03
        spl(j,4)=phi4*dfak34
       else
        spl(j,2)=0.0  
        spl(j,3)=0.0
        spl(j,1)=0.0
        spl(j,4)=0.0
       endif

C *** AN inserted interpolation of q*f(q)
c$$$       spl(j,2)=spl(j,2)*x1/xn  
c$$$       spl(j,3)=spl(j,3)*x2/xn
c$$$       spl(j,1)=spl(j,1)*x0/xn
c$$$       spl(j,4)=spl(j,4)*x3/xn
C *** AN

20    continue
      end subroutine



      SUBROUTINE SPLFAK(X,N)
C     =================
C     THE 2ND DERIVATIVES AT THE BOUNDARIES MUST BE ZERO !
C     THE SPLINE IS BUILT UP USING ALL THE N MESHPOINTS .
C
      PARAMETER(NQ=500)
      DIMENSION X(N)
      DIMENSION HI(NQ),U(NQ)
      REAL PI
      COMMON/FAKTOR/FAK1(NQ,NQ),FAK2(NQ,NQ),FAK3(NQ,NQ),
     C              Q(NQ,NQ),C(NQ,NQ)

      IF(N.GT.NQ) STOP 'SPLFAK: NQ nicht gross genug'

      U(1)=0.
      HI(2)=X(2)-X(1)
      N1=N-1
      DO 10 I=2,N1
      AAX=X(I+1)-X(I)
      HI(I+1)=AAX
      BX=X(I+1)-X(I-1)
      CX=X(I)-X(I-1)
      AL=AAX/BX
      AM=1.-AL
      PI=1./(2.-AM*U(I-1))
      U(I)=AL*PI
      DO 20 J=1,N
      Q(1,J)=0.
      H1=0.
      H2=0.
      H3=0.
      IF(J.EQ.I-1) H1=1./(CX*BX)
      IF(J.EQ.I) H2=1./(CX*AAX)
      IF(J.EQ.I+1) H3=1./(AAX*BX)
      Q(I,J)=-PI*(AM*Q(I-1,J)-H1+H2-H3)
  20  CONTINUE
  10  CONTINUE
      N2=N+1
      N3=N+2
      DO 30 K=3,N2
      J1=N3-K
      H1=1./HI(J1+1)
      DO 40 L=1,N
      C(N,L)=0.
      C(J1,L)=Q(J1,L)-C(J1+1,L)*U(J1)
      FAK1(J1,L)=-HI(J1+1)*(2.*C(J1,L)+C(J1+1,L))
      IF(L.EQ.J1) FAK1(J1,L)=FAK1(J1,L)-H1
      IF(L.EQ.J1+1) FAK1(J1,L)=FAK1(J1,L)+H1
      FAK2(J1,L)=3.*C(J1,L)
      FAK3(J1,L)=(C(J1+1,L)-C(J1,L))*H1
      FAK1(N,L)=0.
      FAK2(N,L)=0.
      FAK3(N,L)=0.
  40  CONTINUE
  30  CONTINUE
      RETURN
      END SUBROUTINE
      
      SUBROUTINE SPLMOD(X,N,XA,SPL)
C     =================
C     X  : ARRAY OF N WHERE THE FUNCTION TO BE INTERPOLATED IS KNOWN
C          (USUALLY MESHPPOINTS)
C     XA : POINT OF INTERPOLATION
C
C     F(XA) = SUM(I=1,N) [ SPL(I)*F(X(I)) ]
C
C     (THIS SUMMATION HAS TO BE DONE IN THE CALLING PROGRAM ! )
C
      PARAMETER(NQ=500)
      DIMENSION X(N)
      DIMENSION SPL(N)
      COMMON/FAKTOR/FAK1(NQ,NQ),FAK2(NQ,NQ),FAK3(NQ,NQ),
     C              QDUM(NQ,NQ),CDUM(NQ,NQ)
   1  I=0
   2  I=I+1
      IF(I.GT.N) GO TO 3
      IF(XA.GE.X(I)) GO TO 2
   3  I1=I-1
      DX=XA-X(I1)
      DO 10 J=1,N
      SPL(J)=((FAK3(I1,J)*DX+FAK2(I1,J))*DX+FAK1(I1,J))*DX
      IF(J.EQ.I1) SPL(J)=SPL(J)+1.
  10  CONTINUE
      RETURN
      END SUBROUTINE
     

      END MODULE spline

