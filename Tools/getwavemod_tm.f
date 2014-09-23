
       MODULE getwavemod
       USE spline
       USE clebsch
       USE ISO_C_BINDING

       PRIVATE
       REAL,ALLOCATABLE :: PSI(:,:,:)
       REAL,ALLOCATABLE :: P12P(:),P34P(:),
     X                     P3P(:),
     X                     Q4P(:),QP(:)
       REAL,ALLOCATABLE :: coupfaktor(:)
       LOGICAL,ALLOCATABLE :: chanf(:)
       INTEGER alpha1amax,alpha1bmax,alpha1cmax,alpha2amax,
     X                    alpha2bmax,np12,np3,nq4,nq,np34,
     X                    jmaxNN,mbj,mbt
       INTEGER,ALLOCATABLE :: qnl12(:),qns12(:),qnj12(:),
     X                        qnt12(:),qnl3(:),qnI3(:),qntau(:)
       COMPLEX,ALLOCATABLE ::  spharm12(:)
       COMPLEX,ALLOCATABLE ::  spharm3(:)
       REAL,ALLOCATABLE    ::  sqfakt(:,:)
       REAL,ALLOCATABLE    ::  wavef(:)

       PUBLIC readwave,getwave,getwave_s,getwave_d,mbjreset
       character(kind=C_CHAR), bind(C, name="he3filename"):: datname*256
       PUBLIC :: datname

      CONTAINS

C     getwave calculates the wave function in 3D representation
C     in fm**3
C     p1x,p1y,p1z,ms1,mt1   momentum in fm-1, spin projection and
c                    isospin projection  of first nucleon 
C     p2x,p2y,p2z,ms2,mt2   momentum in fm-1, spin projection and
c                    isospin projection  of second nucleon 
c     ms3,mt3 spin projection and isospin projection  of third nucleon
c     the momentum of the third nucleon is calculated 
c     from p1+p2+p3=0
c     before calling getwave readwave is required
  
      FUNCTION getwave(p1x,p1y,p1z,ms1,mt1,
     X                 p2x,p2y,p2z,ms2,mt2,ms3,mt3)
     Y   		bind(C, name="gethe3wave")
       IMPLICIT NONE
       REAL(kind=C_DOUBLE) p1x,p1y,p1z,p2x,p2y,p2z
       INTEGER ms1,mt1,ms2,mt2,ms3,mt3
       DOUBLE COMPLEX getwave,csum,sph12,sph3
       REAL fakt,faktw
       REAL p3x,p3y,p3z,p12x,p12y,p12z,p12mag,p3mag
       REAL phi3,phi12,x3,x12
       REAL newp(1),sum12,sum3
       REAL spl12(4),spl3(4),spfakt12,spfakt3
       INTEGER ind12(4),ind3(4)
       INTEGER is12,is3,indx12,indx3
      
       INTEGER l12,l3,I3,j12,s12,t12,mu,nu,tau
       INTEGER alpha,mindx,msum,cindx,isp12,isp3
C      prepare Jacobi coordinates,magnitudes and angles

       p3x=-p1x-p2x
       p3y=-p1y-p2y
       p3z=-p1z-p2z
       p12x=0.5*(p1x-p2x)
       p12y=0.5*(p1y-p2y)
       p12z=0.5*(p1z-p2z)
       
       p12mag=sqrt(p12x*p12x+p12y*p12y+p12z*p12z)
       p3mag=sqrt(p3x*p3x+p3y*p3y+p3z*p3z)
       
       IF(p12mag.NE.0) THEN
         x12=p12z/p12mag
       ELSE
         x12=1.0
       END IF

       IF(p3mag.NE.0) THEN
         x3=p3z/p3mag
       ELSE
         x3=1.0
       END IF
       
       IF(p12x.GT.0.0) THEN
         phi12=atan(p12y/p12x)
       ELSE IF(p12x.LT.0.0) THEN
         IF(p12y.GT.0.0) THEN
           phi12=acos(-1.0)+atan(p12y/p12x)
         ELSE
           phi12=-acos(-1.0)+atan(p12y/p12x)
         END IF
       ELSE                     ! p12x.EQ.0
         IF(p12y.GT.0.0) THEN
           phi12=acos(-1.0)/2.0
         ELSE
           phi12=-acos(-1.0)/2.0
         END IF           
       END IF

       IF(p3x.GT.0.0) THEN
         phi3=atan(p3y/p3x)
       ELSE IF(p3x.LT.0.0) THEN
         IF(p3y.GT.0.0) THEN
           phi3=acos(-1.0)+atan(p3y/p3x)
         ELSE
           phi3=-acos(-1.0)+atan(p3y/p3x)
         END IF
       ELSE                     ! p3x.EQ.0
         IF(p3y.GT.0.0) THEN
           phi3=acos(-1.0)/2.0
         ELSE
           phi3=-acos(-1.0)/2.0
         END IF           
       END IF
  
 
       
C     prepare spherical harmonics

       DO l12=0,jmaxNN+1
        isp12=l12*(2*jmaxNN+3)+jmaxNN+2
        DO mu=0,l12
*       write(6,*)"numbers",mu+isp12
         spharm12(mu+isp12) 
     X     =sqfakt(l12,mu)*PLGNDR(l12,mu,x12)*CEXP((0,1.0)*mu*phi12)
         spharm12(-mu+isp12)= (-1.0)**mu*CONJG(spharm12(mu+isp12))
        END DO
       END DO


       DO l3=0,jmaxNN+1
        isp3=l3*(2*jmaxNN+3)+jmaxNN+2
        DO nu=0,l3
         spharm3(nu+isp3)
     X     =sqfakt(l3,nu)*PLGNDR(l3,nu,x3)*CEXP((0,1.0)*nu*phi3)
         spharm3(-nu+isp3)=(-1.0)**nu*CONJG(spharm3(nu+isp3))
        END DO
       END DO

*       getwave = (2,0)
*       return


C      prepare wave function at off-grid momenta

       newp(1)=p12mag
       CALL cubherm(P12P,NP12,newp,1,spl12,ind12) 
       newp(1)=p3mag
       CALL cubherm(P3P,NP3,newp,1,spl3,ind3)

       wavef(1:alpha1amax)=0.0
c       CALL fort_pat_start(2)

       DO alpha=1,alpha1amax
        wavef(alpha)=wavef(alpha)
     X       +spl12(1)*spl3(1)*psi(alpha,ind12(1),ind3(1))
     X       +spl12(2)*spl3(1)*psi(alpha,ind12(2),ind3(1))
     X       +spl12(3)*spl3(1)*psi(alpha,ind12(3),ind3(1))
     X       +spl12(4)*spl3(1)*psi(alpha,ind12(4),ind3(1))
       END DO
       DO alpha=1,alpha1amax
        wavef(alpha)=wavef(alpha)
     X       +spl12(1)*spl3(2)*psi(alpha,ind12(1),ind3(2))
     X       +spl12(2)*spl3(2)*psi(alpha,ind12(2),ind3(2))
     X       +spl12(3)*spl3(2)*psi(alpha,ind12(3),ind3(2))
     X       +spl12(4)*spl3(2)*psi(alpha,ind12(4),ind3(2))
       END DO
       DO alpha=1,alpha1amax
        wavef(alpha)=wavef(alpha)
     X       +spl12(1)*spl3(3)*psi(alpha,ind12(1),ind3(3))
     X       +spl12(2)*spl3(3)*psi(alpha,ind12(2),ind3(3))
     X       +spl12(3)*spl3(3)*psi(alpha,ind12(3),ind3(3))
     X       +spl12(4)*spl3(3)*psi(alpha,ind12(4),ind3(3))
       END DO
       DO alpha=1,alpha1amax
        wavef(alpha)=wavef(alpha)
     X       +spl12(1)*spl3(4)*psi(alpha,ind12(1),ind3(4))
     X       +spl12(2)*spl3(4)*psi(alpha,ind12(2),ind3(4))
     X       +spl12(3)*spl3(4)*psi(alpha,ind12(3),ind3(4))
     X       +spl12(4)*spl3(4)*psi(alpha,ind12(4),ind3(4))
       END DO
c       CALL fort_pat_stop(2)

C     recouple wavefunction
       
       csum=0.0
       
       mindx=(ms1+1)/2+2*(ms2+1)/2+4*(ms3+1)/2
     X      +8*(mt1+1)/2+16*(mt2+1)/2+32*(mt3+1)/2

       msum=(mbj-ms1-ms2-ms3)/2
     
c       CALL fort_pat_start(3)

       DO alpha=1,alpha1amax
        cindx=(alpha+mindx*alpha1amax-1)*(2*jmaxNN+3)+jmaxNN+2
        l12=qnl12(alpha)
        l3=qnl3(alpha)
        IF(chanf(alpha+alpha1amax*mindx)) THEN
          isp12=l12*(2*jmaxNN+3)+jmaxNN+2
          isp3=l3*(2*jmaxNN+3)+jmaxNN+2
          DO mu=-l12,l12
           nu=msum-mu
           IF(abs(nu).LE.l3) THEN
             csum=csum
     X          +wavef(alpha)*coupfaktor(mu+cindx)
     X            *spharm12(mu+isp12)*spharm3(nu+isp3)
           END IF
          END DO
        END IF
       END DO       
c       CALL fort_pat_stop(3)

       getwave=csum
!        write(6,*)'wf ', getwave

      END FUNCTION getwave


      FUNCTION getwave_s(nal,p1x,p1y,p1z,ms1,mt1,
     X                 p2x,p2y,p2z,ms2,mt2,ms3,mt3)

*****************************************************************
* getwave_s
*****************************************************************
*     getwave_s calculates the wave function in 3D representation
*     for given channel nal (the wave functions is in fm**3)
*     in difference from getwave, it calculates the partial 
*     wave for given channel
*     p1x,p1y,p1z,ms1,mt1   momentum in fm-1, spin projection and
*                    isospin projection  of first nucleon 
*     p2x,p2y,p2z,ms2,mt2   momentum in fm-1, spin projection and
*                    isospin projection  of second nucleon 
*     ms3,mt3 spin projection and isospin projection  of third nucleon
*     the momentum of the third nucleon is calculated 
*     from p1+p2+p3=0
*     before calling getwave readwave is required
*****************************************************************  
       IMPLICIT NONE
       REAL p1x,p1y,p1z,p2x,p2y,p2z
       INTEGER ms1,mt1,ms2,mt2,ms3,mt3
       COMPLEX getwave_s,csum,sph12,sph3
       REAL fakt,faktw
       REAL p3x,p3y,p3z,p12x,p12y,p12z,p12mag,p3mag
       REAL phi3,phi12,x3,x12
       REAL newp(1),sum12,sum3
       REAL spl12(4),spl3(4),spfakt12,spfakt3
       INTEGER ind12(4),ind3(4)
       INTEGER is12,is3,indx12,indx3
      
       INTEGER l12,l3,I3,j12,s12,t12,mu,nu,tau
       INTEGER alpha,mindx,msum,cindx,isp12,isp3
       integer nal
C      prepare Jacobi coordinates,magnitudes and angles
 
       p3x=-p1x-p2x
       p3y=-p1y-p2y
       p3z=-p1z-p2z
       p12x=0.5*(p1x-p2x)
       p12y=0.5*(p1y-p2y)
       p12z=0.5*(p1z-p2z)
       
       p12mag=sqrt(p12x*p12x+p12y*p12y+p12z*p12z)
       p3mag=sqrt(p3x*p3x+p3y*p3y+p3z*p3z)
       
************************************************
* Calculation of polar angles of Jacobi momenta
************************************************
       IF(p12mag.NE.0) THEN
         x12=p12z/p12mag
       ELSE
         x12=1.0
       END IF

       IF(p3mag.NE.0) THEN
         x3=p3z/p3mag
       ELSE
         x3=1.0
       END IF
       
       IF(p12x.GT.0.0) THEN
         phi12=atan(p12y/p12x)
       ELSE IF(p12x.LT.0.0) THEN
         IF(p12y.GT.0.0) THEN
           phi12=acos(-1.0)+atan(p12y/p12x)
         ELSE
           phi12=-acos(-1.0)+atan(p12y/p12x)
         END IF
       ELSE                     ! p12x.EQ.0
         IF(p12y.GT.0.0) THEN
           phi12=acos(-1.0)/2.0
         ELSE
           phi12=-acos(-1.0)/2.0
         END IF           
       END IF

       IF(p3x.GT.0.0) THEN
         phi3=atan(p3y/p3x)
       ELSE IF(p3x.LT.0.0) THEN
         IF(p3y.GT.0.0) THEN
           phi3=acos(-1.0)+atan(p3y/p3x)
         ELSE
           phi3=-acos(-1.0)+atan(p3y/p3x)
         END IF
       ELSE                     ! p3x.EQ.0
         IF(p3y.GT.0.0) THEN
           phi3=acos(-1.0)/2.0
         ELSE
           phi3=-acos(-1.0)/2.0
         END IF           
       END IF

**************************************************
*            prepare spherical harmonics
**************************************************
       DO l12=0,jmaxNN+1
        isp12=l12*(2*jmaxNN+3)+jmaxNN+2
        DO mu=0,l12
         spharm12(mu+isp12) 
     X     =sqfakt(l12,mu)*PLGNDR(l12,mu,x12)*CEXP((0,1.0)*mu*phi12)
         spharm12(-mu+isp12)= (-1.0)**mu*CONJG(spharm12(mu+isp12))
        END DO
       END DO


       DO l3=0,jmaxNN+1
        isp3=l3*(2*jmaxNN+3)+jmaxNN+2
        DO nu=0,l3
         spharm3(nu+isp3)
     X     =sqfakt(l3,nu)*PLGNDR(l3,nu,x3)*CEXP((0,1.0)*nu*phi3)
         spharm3(-nu+isp3)=(-1.0)**nu*CONJG(spharm3(nu+isp3))
        END DO
       END DO


*******************************************************
*      prepare wave function at off-grid momenta
*******************************************************
       newp(1)=p12mag
       CALL cubherm(P12P,NP12,newp,1,spl12,ind12) 
       newp(1)=p3mag
       CALL cubherm(P3P,NP3,newp,1,spl3,ind3)

       wavef(1:alpha1amax)=0.0
c       CALL fort_pat_start(2)

*       alpha = nal
       DO alpha=1,nal !,14 !,alpha1amax
        wavef(alpha)=wavef(alpha)
     X       +spl12(1)*spl3(1)*psi(alpha,ind12(1),ind3(1))
     X       +spl12(2)*spl3(1)*psi(alpha,ind12(2),ind3(1))
     X       +spl12(3)*spl3(1)*psi(alpha,ind12(3),ind3(1))
     X       +spl12(4)*spl3(1)*psi(alpha,ind12(4),ind3(1))
       END DO
       DO alpha=1,nal !15,14 !,alpha1amax
        wavef(alpha)=wavef(alpha)
     X       +spl12(1)*spl3(2)*psi(alpha,ind12(1),ind3(2))
     X       +spl12(2)*spl3(2)*psi(alpha,ind12(2),ind3(2))
     X       +spl12(3)*spl3(2)*psi(alpha,ind12(3),ind3(2))
     X       +spl12(4)*spl3(2)*psi(alpha,ind12(4),ind3(2))
       END DO
       DO alpha=1,nal !15,14 !,alpha1amax
        wavef(alpha)=wavef(alpha)
     X       +spl12(1)*spl3(3)*psi(alpha,ind12(1),ind3(3))
     X       +spl12(2)*spl3(3)*psi(alpha,ind12(2),ind3(3))
     X       +spl12(3)*spl3(3)*psi(alpha,ind12(3),ind3(3))
     X       +spl12(4)*spl3(3)*psi(alpha,ind12(4),ind3(3))
       END DO
       DO alpha=1,nal !15,14 !,alpha1amax
        wavef(alpha)=wavef(alpha)
     X       +spl12(1)*spl3(4)*psi(alpha,ind12(1),ind3(4))
     X       +spl12(2)*spl3(4)*psi(alpha,ind12(2),ind3(4))
     X       +spl12(3)*spl3(4)*psi(alpha,ind12(3),ind3(4))
     X       +spl12(4)*spl3(4)*psi(alpha,ind12(4),ind3(4))
       END DO
c       CALL fort_pat_stop(2)

******************************************
*     recouple wavefunction
******************************************
       
       csum=0.0       
       mindx=(ms1+1)/2+2*(ms2+1)/2+4*(ms3+1)/2
     X      +8*(mt1+1)/2+16*(mt2+1)/2+32*(mt3+1)/2
       msum=(mbj-ms1-ms2-ms3)/2
     
c       CALL fort_pat_start(3)

       DO alpha=1,nal !15,14 !,alpha1amax
        cindx=(alpha+mindx*alpha1amax-1)*(2*jmaxNN+3)+jmaxNN+2
        l12=qnl12(alpha)
        l3=qnl3(alpha)
        IF(chanf(alpha+alpha1amax*mindx)) THEN
          isp12=l12*(2*jmaxNN+3)+jmaxNN+2
          isp3=l3*(2*jmaxNN+3)+jmaxNN+2
          DO mu=-l12,l12
           nu=msum-mu
           IF(abs(nu).LE.l3) THEN
             csum=csum
     X          +wavef(alpha)*coupfaktor(mu+cindx)
     X            *spharm12(mu+isp12)*spharm3(nu+isp3)
           END IF
          END DO
        END IF
       END DO       
c       CALL fort_pat_stop(3)
       getwave_s=csum
      END FUNCTION getwave_s


*====================================================
      FUNCTION getwave_d(p1x,p1y,p1z,ms1,mt1,
     X                 p2x,p2y,p2z,ms2,mt2,ms3,mt3)

*****************************************************************
* getwave_d
*****************************************************************
*     getwave_d calculates the wave function in 3D representation
*     for channels 15,16,17,18, which correspond to deuteron state 
*     for 12 nucleon system
* (the wave functions is in fm**3)
*     in difference from getwave, it calculates the partial 
*     wave for given channel
*     p1x,p1y,p1z,ms1,mt1   momentum in fm-1, spin projection and
*                    isospin projection  of first nucleon 
*     p2x,p2y,p2z,ms2,mt2   momentum in fm-1, spin projection and
*                    isospin projection  of second nucleon 
*     ms3,mt3 spin projection and isospin projection  of third nucleon
*     the momentum of the third nucleon is calculated 
*     from p1+p2+p3=0
*     before calling getwave readwave is required
*****************************************************************  
       IMPLICIT NONE
       REAL p1x,p1y,p1z,p2x,p2y,p2z
       INTEGER ms1,mt1,ms2,mt2,ms3,mt3
       COMPLEX getwave_d,csum,sph12,sph3
       REAL fakt,faktw
       REAL p3x,p3y,p3z,p12x,p12y,p12z,p12mag,p3mag
       REAL phi3,phi12,x3,x12
       REAL newp(1),sum12,sum3
       REAL spl12(4),spl3(4),spfakt12,spfakt3
       INTEGER ind12(4),ind3(4)
       INTEGER is12,is3,indx12,indx3
      
       INTEGER l12,l3,I3,j12,s12,t12,mu,nu,tau
       INTEGER alpha,mindx,msum,cindx,isp12,isp3
       integer nal
C      prepare Jacobi coordinates,magnitudes and angles
 
       p3x=-p1x-p2x
       p3y=-p1y-p2y
       p3z=-p1z-p2z
       p12x=0.5*(p1x-p2x)
       p12y=0.5*(p1y-p2y)
       p12z=0.5*(p1z-p2z)
       
       p12mag=sqrt(p12x*p12x+p12y*p12y+p12z*p12z)
       p3mag=sqrt(p3x*p3x+p3y*p3y+p3z*p3z)
       
************************************************
* Calculation of polar angles of Jacobi momenta
************************************************
       IF(p12mag.NE.0) THEN
         x12=p12z/p12mag
       ELSE
         x12=1.0
       END IF

       IF(p3mag.NE.0) THEN
         x3=p3z/p3mag
       ELSE
         x3=1.0
       END IF
       
       IF(p12x.GT.0.0) THEN
         phi12=atan(p12y/p12x)
       ELSE IF(p12x.LT.0.0) THEN
         IF(p12y.GT.0.0) THEN
           phi12=acos(-1.0)+atan(p12y/p12x)
         ELSE
           phi12=-acos(-1.0)+atan(p12y/p12x)
         END IF
       ELSE                     ! p12x.EQ.0
         IF(p12y.GT.0.0) THEN
           phi12=acos(-1.0)/2.0
         ELSE
           phi12=-acos(-1.0)/2.0
         END IF           
       END IF

       IF(p3x.GT.0.0) THEN
         phi3=atan(p3y/p3x)
       ELSE IF(p3x.LT.0.0) THEN
         IF(p3y.GT.0.0) THEN
           phi3=acos(-1.0)+atan(p3y/p3x)
         ELSE
           phi3=-acos(-1.0)+atan(p3y/p3x)
         END IF
       ELSE                     ! p3x.EQ.0
         IF(p3y.GT.0.0) THEN
           phi3=acos(-1.0)/2.0
         ELSE
           phi3=-acos(-1.0)/2.0
         END IF           
       END IF

**************************************************
*            prepare spherical harmonics
**************************************************
       DO l12=0,jmaxNN+1
        isp12=l12*(2*jmaxNN+3)+jmaxNN+2
        DO mu=0,l12
         spharm12(mu+isp12) 
     X     =sqfakt(l12,mu)*PLGNDR(l12,mu,x12)*CEXP((0,1.0)*mu*phi12)
         spharm12(-mu+isp12)= (-1.0)**mu*CONJG(spharm12(mu+isp12))
        END DO
       END DO


       DO l3=0,jmaxNN+1
        isp3=l3*(2*jmaxNN+3)+jmaxNN+2
        DO nu=0,l3
         spharm3(nu+isp3)
     X     =sqfakt(l3,nu)*PLGNDR(l3,nu,x3)*CEXP((0,1.0)*nu*phi3)
         spharm3(-nu+isp3)=(-1.0)**nu*CONJG(spharm3(nu+isp3))
        END DO
       END DO


*******************************************************
*      prepare wave function at off-grid momenta
*******************************************************
       newp(1)=p12mag
       CALL cubherm(P12P,NP12,newp,1,spl12,ind12) 
       newp(1)=p3mag
       CALL cubherm(P3P,NP3,newp,1,spl3,ind3)

       wavef(1:alpha1amax)=0.0
c       CALL fort_pat_start(2)

*       alpha = nal
       DO alpha=15, 18 !nal !,14 !,alpha1amax
        wavef(alpha)=wavef(alpha)
     X       +spl12(1)*spl3(1)*psi(alpha,ind12(1),ind3(1))
     X       +spl12(2)*spl3(1)*psi(alpha,ind12(2),ind3(1))
     X       +spl12(3)*spl3(1)*psi(alpha,ind12(3),ind3(1))
     X       +spl12(4)*spl3(1)*psi(alpha,ind12(4),ind3(1))
       END DO
       DO alpha=15,18 !nal !15,14 !,alpha1amax
        wavef(alpha)=wavef(alpha)
     X       +spl12(1)*spl3(2)*psi(alpha,ind12(1),ind3(2))
     X       +spl12(2)*spl3(2)*psi(alpha,ind12(2),ind3(2))
     X       +spl12(3)*spl3(2)*psi(alpha,ind12(3),ind3(2))
     X       +spl12(4)*spl3(2)*psi(alpha,ind12(4),ind3(2))
       END DO
       DO alpha=15, 18 !nal !15,14 !,alpha1amax
        wavef(alpha)=wavef(alpha)
     X       +spl12(1)*spl3(3)*psi(alpha,ind12(1),ind3(3))
     X       +spl12(2)*spl3(3)*psi(alpha,ind12(2),ind3(3))
     X       +spl12(3)*spl3(3)*psi(alpha,ind12(3),ind3(3))
     X       +spl12(4)*spl3(3)*psi(alpha,ind12(4),ind3(3))
       END DO
       DO alpha=15,18 !nal !15,14 !,alpha1amax
        wavef(alpha)=wavef(alpha)
     X       +spl12(1)*spl3(4)*psi(alpha,ind12(1),ind3(4))
     X       +spl12(2)*spl3(4)*psi(alpha,ind12(2),ind3(4))
     X       +spl12(3)*spl3(4)*psi(alpha,ind12(3),ind3(4))
     X       +spl12(4)*spl3(4)*psi(alpha,ind12(4),ind3(4))
       END DO
c       CALL fort_pat_stop(2)

******************************************
*     recouple wavefunction
******************************************
       
       csum=0.0       
       mindx=(ms1+1)/2+2*(ms2+1)/2+4*(ms3+1)/2
     X      +8*(mt1+1)/2+16*(mt2+1)/2+32*(mt3+1)/2
       msum=(mbj-ms1-ms2-ms3)/2
     
c       CALL fort_pat_start(3)

       DO alpha=15, 18 !nal !15,14 !,alpha1amax
        cindx=(alpha+mindx*alpha1amax-1)*(2*jmaxNN+3)+jmaxNN+2
        l12=qnl12(alpha)
        l3=qnl3(alpha)
        IF(chanf(alpha+alpha1amax*mindx)) THEN
          isp12=l12*(2*jmaxNN+3)+jmaxNN+2
          isp3=l3*(2*jmaxNN+3)+jmaxNN+2
          DO mu=-l12,l12
           nu=msum-mu
           IF(abs(nu).LE.l3) THEN
             csum=csum
     X          +wavef(alpha)*coupfaktor(mu+cindx)
     X            *spharm12(mu+isp12)*spharm3(nu+isp3)
           END IF
          END DO
        END IF
       END DO       
c       CALL fort_pat_stop(3)
       getwave_d=csum
      END FUNCTION getwave_d













c     readwave reads in the wave function in file DATNAME
c     the third component of the isospin is fixed by the data file
c     the third component of the spin of the 3n system is given 
c     by mbjset= +/- 1, it can only be called once.
c     if usecdep = .true. then it uses the isospin violating parts 
c     of the wave function, otherwise it neglects them (good approximation,faster)

      SUBROUTINE readwave(mbjset,usecdep) bind(C, name="readhe3wave")
       IMPLICIT NONE

!        CHARACTER(LEN=256) DATNAME
       LOGICAL usecdep

       INTEGER j3max,lammax,l4max,potnr_yak,potynnr_yak,
     X           l3max,lsummax,btmax,taumax,mbjset 
       LOGICAL cdep,cdep_app       
       INTEGER alpha,beta,ip,iq
       INTEGER alpha1a,l12,s12,j12,t12,l3,I3,j3,tau,
     X        l4,I4,bj,t4,bt,alpha1b,t3,alpha1c,t2,alpha2a,l34,
     $        s34,j34,t34,lam,I,alpha2b

       CHARACTER(LEN=80) erkennung
       REAL,ALLOCATABLE :: AP12(:),AP3(:),AQ4(:),AQ(:),AP34(:)
       REAL,ALLOCATABLE :: PSIDUM1(:,:,:),PSIDUM2(:,:,:)
       REAL sum,sump,sumq
C     define polarization of 3n system
       mbj=mbjset
!        write(6,*)'bla', DATNAME
       
*21apr       WRITE(*,*) 'Read 3N wave function'

       OPEN(84,FILE=DATNAME,FORM='FORMATTED',STATUS='OLD')
       
       READ(84,*) erkennung
*21apr       write(4,*) erkennung
       READ(84,*) jmaxNN,j3max,lammax,l4max,potnr_yak,potynnr_yak,
     X      l3max,lsummax,mbt,btmax,taumax,cdep,cdep_app
*21apr       write(4,*) jmaxNN,j3max,lammax,l4max,potnr_yak,potynnr_yak,
*21apr     X      l3max,lsummax,mbt,btmax,taumax,cdep,cdep_app
       READ(84,*) NP12,NP3,NQ4,NP34,NQ
*21apr       write(4,*) NP12,NP3,NQ4,NP34,NQ
       READ(84,*) alpha1amax,alpha1bmax,alpha1cmax,
     X      alpha2amax,alpha2bmax
*21apr       write(4,*) alpha1amax,alpha1bmax,alpha1cmax,
*21apr     X      alpha2amax,alpha2bmax
       
       
*21apr       WRITE(*,*) 'Signature: ',erkennung
*21apr       WRITE(*,*) 'MESH: ',NP12,NP3
*21apr       WRITE(*,*) 'CHANNEL: ',jmaxNN,alpha1amax
*21apr       WRITE(*,*) 'CHARGE: ',mbt/2
       
       IF(NQ4.NE.1) STOP'inconsistent NQ4'
       
       ALLOCATE(PSI(alpha1amax,NP12,NP3))
       ALLOCATE(PSIDUM1(NP12,NP3,NQ4))
       ALLOCATE(PSIDUM2(NP12,NP34,NQ))
       
       ALLOCATE(P12P(NP12),P3P(NP3),P34P(NP34),Q4P(NQ4),QP(NQ))
       ALLOCATE(AP12(NP12),AP3(NP3),AP34(NP34),AQ4(NQ4),AQ(NQ))
       ALLOCATE(qnl12(alpha1amax))
       ALLOCATE(qns12(alpha1amax))
       ALLOCATE(qnj12(alpha1amax))
       ALLOCATE(qnt12(alpha1amax))
       ALLOCATE(qnl3(alpha1amax))
       ALLOCATE(qnI3(alpha1amax))
       ALLOCATE(qntau(alpha1amax))

       
       READ(84,*) P12P,P3P,Q4P,P34P,QP
*21apr       write(4,*) P12P,P3P,Q4P,P34P,QP
       READ(84,*) AP12,AP3,AQ4,AP34,AQ
*21apr       write(4,*) AP12,AP3,AQ4,AP34,AQ
       
       DO alpha1a=1,alpha1amax
        READ(84,*) l12,s12,j12,t12,l3,I3,j3,tau,l4,I4,
     X       bj,t4,bt
*21apr        write(4,13) l12,s12,j12,t12,l3,I3,j3,tau,l4,I4,
*21apr     X       bj,t4,bt
 13     format(13(2x,I2))
        qnl12(alpha1a)=l12
        qns12(alpha1a)=s12
        qnj12(alpha1a)=j12
        qnt12(alpha1a)=t12
        qnl3(alpha1a)=l3
        qnI3(alpha1a)=I3
        qntau(alpha1a)=tau
*        write(6,*)"tau",tau
        IF(t4.NE.1) STOP'inconsistent t4'
        IF(l4.NE.0) STOP'inconsistent l4'
        IF(I4.NE.1) STOP'inconsistent I4'
        IF(j3.NE.1) STOP'inconsistent j3'
        
       END DO
       
       DO alpha1b=1,alpha1bmax
        READ(84,*) l12,s12,j12,t12,l3,I3,j3,
     X       t3,tau,l4,I4,bj,bt
*21apr        write(4,*) l12,s12,j12,t12,l3,I3,j3,
*21apr     X       t3,tau,l4,I4,bj,bt
       END DO
       
       DO alpha1c=1,alpha1cmax
        READ(84,*) l12,s12,j12,t2,t12,l3,
     $       I3,j3,tau,l4,I4,bj,bt
*21apr        write(4,*) l12,s12,j12,t2,t12,l3,
*21apr     $       I3,j3,tau,l4,I4,bj,bt
       END DO
       
       DO alpha2a=1,alpha2amax
        READ(84,*) l12,s12,j12,t12,l34,
     $       s34,j34,t4,t34,lam,I,bt
*21apr        write(4,*) l12,s12,j12,t12,l34,
*21apr     $       s34,j34,t4,t34,lam,I,bt
       END DO
       
       DO alpha2b=1,alpha2bmax
        READ(84,*) l12,s12,j12,t2,t12,l34,
     $       s34,j34,t34,lam,I,bt
*21apr        write(4,*) l12,s12,j12,t2,t12,l34,
*21apr     $       s34,j34,t34,lam,I,bt
       END DO
       
       
       
       DO alpha=1,alpha1amax
        READ(84,*) PSIDUM1
*21apr        write(4,*) PSIDUM1

        DO iq=1,NP3
         DO ip=1,NP12
          PSI(alpha,ip,iq)=PSIDUM1(ip,iq,1)
         END DO
        END DO

       END DO
       
       DO alpha=1,alpha1bmax
        READ(84,*) PSIDUM1
*21apr        write(4,*) PSIDUM1
       END DO
       
       
       DO alpha=1,alpha1cmax
        READ(84,*) PSIDUM1
*21apr        write(4,*) PSIDUM1
       END DO
       
       IF(NQ4.ne.1) THEN
         
         DO beta=1,alpha2amax
          READ(84,*) PSIDUM2
*21apr          write(4,*) PSIDUM2
         END DO
         
         DO beta=1,alpha2bmax
          READ(84,*) PSIDUM2
*21apr          write(4,*) PSIDUM2
         END DO
         
       END IF                   ! NQ4
       
       CLOSE(84)
       
C        NORMCHECK

       sum=0.0
       DO alpha=1,alpha1amax
        sumq=0.0
        DO iq=1,NP3
         sump=0.0
         DO ip=1,NP12  
          sump=sump+(P12P(ip)*PSI(alpha,ip,iq))**2*AP12(ip)
         END DO
         sumq=sumq+sump*AP3(iq)*P3P(iq)**2
        END DO
        sum=sum+sumq
       END DO

*21par       WRITE(*,*) 'NORM: ',sum

       sum=0.0
       DO alpha=1,alpha1amax
        sum=sum+PSI(alpha,1,1)**2
       END DO

*21apr       WRITE(*,*) 'ZERO: ',sum
      
       DEALLOCATE(PSIDUM1,PSIDUM2,AP12,AP3,AP34,AQ4)

*21apr       WRITE(*,*) 'Preparation of Clebsch-Gordan coefficients'
     

       CALL fkltt    ! preparing I! and sqrt(I!) 
       CALL prepcoup(mbt/2,mbj,usecdep) ! preparing Clebsch- Gordon couplings
        
       END SUBROUTINE readwave


      SUBROUTINE mbjreset(mbjset,usecdep) bind(C, name="he3reset")

C     mbjreset changes the spin orientation of the 3N system to mbjset = -/+1
C     without rereading the file 
C     it can be called more than once
c     if usecdep = true then it uses the isospin violating parts 
c     of the wave function, otherwise it neglects them
       IMPLICIT NONE
       LOGICAL usecdep,flag
       INTEGER l12,l3,I3,j12,s12,mu,nu,
     X         ms1,ms2,ms3,mt1,mt2,mt3,tau,t12,mbjset
       INTEGER alpha,mindx,cindx
       REAL faktl12
C     define polarization of 3n system
       mbj=mbjset

       coupfaktor(1:alpha1amax*2**6*(2*jmaxNN+3))=0.0
       chanf(1:alpha1amax*2**6)=.false.

       DO alpha=1,alpha1amax
        l12=qnl12(alpha)
        s12=qns12(alpha)
        j12=qnj12(alpha)
        t12=qnt12(alpha)
        l3=qnl3(alpha)
        I3=qnI3(alpha)
        tau=qntau(alpha)

        IF(usecdep .OR. tau.EQ.1) THEN
          DO mt1=-1,1,2
           DO mt2=-1,1,2
            DO mt3=-1,1,2
             IF(mt1+mt2+mt3.EQ.mbt/2) THEN
               DO ms1=-1,1,2
                DO ms2=-1,1,2
                 DO ms3=-1,1,2
                  flag=.false.
                  DO mu=-l12,l12
                   nu=(mbj-ms1-ms2-ms3-2*mu)/2
                   mindx=(ms1+1)/2+2*(ms2+1)/2+4*(ms3+1)/2
     X                  +8*(mt1+1)/2+16*(mt2+1)/2+32*(mt3+1)/2
                   cindx=(alpha+mindx*alpha1amax-1)*(2*jmaxNN+3)
     X                    +jmaxNN+2
                   coupfaktor(mu+cindx)
     X                =CG(2*l3,1,I3,2*nu,ms3,2*nu+ms3)
     X                *CG(2*l12,2*s12,2*j12,2*mu,ms1+ms2,2*mu+ms1+ms2)
     X                *CG(1,1,2*s12,ms1,ms2,ms1+ms2)
     X                *CG(2*j12,I3,1,2*mu+ms1+ms2,2*nu+ms3,mbj)
     X                *CG(1,1,t12,mt1,mt2,mt1+mt2)
     X                *CG(t12,1,tau,mt1+mt2,mt3,mbt/2)       
                   IF(coupfaktor(mu+cindx).NE.0.0) flag=.true.
                  END DO
                  chanf(alpha+alpha1amax*mindx)=flag
                 END DO
                END DO
               END DO
             END IF
            END DO
           END DO
          END DO
        END IF
       END DO

      END SUBROUTINE  mbjreset


      SUBROUTINE prepcoup(mbt,mbj,usecdep)
       IMPLICIT NONE
       INTEGER mbt,mbj,l12,l3,I3,j12,s12,mu,nu,
     X         ms1,ms2,ms3,mt1,mt2,mt3,tau,t12
       INTEGER alpha,mindx,cindx
       REAL faktl12
       LOGICAL usecdep,flag

       ALLOCATE(coupfaktor(alpha1amax*2**6*(2*jmaxNN+3)))
       ALLOCATE(chanf(alpha1amax*2**6))
       coupfaktor(1:alpha1amax*2**6*(2*jmaxNN+3))=0.0
       chanf(1:alpha1amax*2**6)=.false.

       DO alpha=1,alpha1amax
        l12=qnl12(alpha)
        s12=qns12(alpha)
        j12=qnj12(alpha)
        t12=qnt12(alpha)
        l3=qnl3(alpha)
        I3=qnI3(alpha)
        tau=qntau(alpha)

        IF(usecdep .OR. tau.EQ.1) THEN
          DO mt1=-1,1,2
           DO mt2=-1,1,2
            DO mt3=-1,1,2
             IF(mt1+mt2+mt3.EQ.mbt) THEN
               DO ms1=-1,1,2
                DO ms2=-1,1,2
                 DO ms3=-1,1,2
                  flag=.false.
                  DO mu=-l12,l12
                   nu=(mbj-ms1-ms2-ms3-2*mu)/2
                   mindx=(ms1+1)/2+2*(ms2+1)/2+4*(ms3+1)/2
     X                  +8*(mt1+1)/2+16*(mt2+1)/2+32*(mt3+1)/2
                   cindx=(alpha+mindx*alpha1amax-1)*(2*jmaxNN+3)
     X                    +jmaxNN+2
                   coupfaktor(mu+cindx)
     X                =CG(2*l3,1,I3,2*nu,ms3,2*nu+ms3)
     X                *CG(2*l12,2*s12,2*j12,2*mu,ms1+ms2,2*mu+ms1+ms2)
     X                *CG(1,1,2*s12,ms1,ms2,ms1+ms2)
     X                *CG(2*j12,I3,1,2*mu+ms1+ms2,2*nu+ms3,mbj)
     X                *CG(1,1,t12,mt1,mt2,mt1+mt2)
     X                *CG(t12,1,tau,mt1+mt2,mt3,mbt)
                   IF(coupfaktor(mu+cindx).NE.0.0) flag=.true.
                  END DO
                  chanf(alpha+alpha1amax*mindx)=flag
                 END DO
                END DO
               END DO
             END IF
            END DO
           END DO
          END DO
        END IF
       END DO

C     preparation for spherical harmonics


       ALLOCATE(spharm12((jmaxNN+2)*(2*jmaxNN+3)),
     X          spharm3((jmaxNN+2)*(2*jmaxNN+3)))
       
       ALLOCATE(sqfakt(0:jmaxNN+1,0:jmaxNN+1))

       DO l12=0,jmaxNN+1
        faktl12=sqrt((2.0*l12+1.0)/(4.0*acos(-1.0)))
        sqfakt(l12,0)=faktl12
        DO mu=1,l12
         faktl12=faktl12/sqrt((l12-mu+1.0)*(l12+mu))
         sqfakt(l12,mu)=faktl12
        END DO
       END DO

C     for interpolated wave function
       ALLOCATE(wavef(alpha1amax))

      END SUBROUTINE prepcoup



      FUNCTION PLGNDR(L,M,X)
      IF(M.LT.0.OR.M.GT.L.OR.ABS(X).GT.1.)PAUSE 'bad arguments'
      PMM=1.
      IF(M.GT.0) THEN
        SOMX2=SQRT((1.-X)*(1.+X))
        FACT=1.
        DO 11 I=1,M
          PMM=-PMM*FACT*SOMX2
          FACT=FACT+2.
11      CONTINUE
      ENDIF
      IF(L.EQ.M) THEN
        PLGNDR=PMM
      ELSE
        PMMP1=X*(2*M+1)*PMM
        IF(L.EQ.M+1) THEN
          PLGNDR=PMMP1
        ELSE
          DO 12 LL=M+2,L
            PLL=(X*(2*LL-1)*PMMP1-(LL+M-1)*PMM)/(LL-M)
            PMM=PMMP1
            PMMP1=PLL
12        CONTINUE
          PLGNDR=PLL
        ENDIF
      ENDIF
      RETURN
      END FUNCTION plgndr

 



       END MODULE getwavemod

