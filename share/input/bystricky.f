
       subroutine f_M(s,t,is1,is2,is1p,is2p,
     &                      irs,f_Mre,f_Mim,B_ampl)
*--------------------------------------------------------------------------
*              Nucleon-nucleon elastic scattering amplitudes 
*
*                 NOTE !!!   IMPLICIT REAL(k-n)   !!!  NOTE
*
*         Parametrisation as in [1] Rev. Mod. Phys. 65, 47(1993)
*                         input parametrs
*
* s                  --  (p_1^mu+p_2^mu)^2
* t                  --  (p_3^mu-p_1^mu)^2
*
* is1,is2,is1p,is2p  -- initial and final spins of scattered nucleons
*
* irs                -- describing isospin configuration 
*
*                              ATTENTION!
* the definition of "irs" from a1f.f assumed here to be GENERALIZED
* (in order to account charge exchange proccesses) as follows:
*                            We still have
* irs = 22 (pp->pp rescattering), irs = 21 (pn->pn rescattering)    
* irs = 121 - np case - knocked-out neutron with rescattering off proton-r    
* irs = 122 - np case - knocked-out neutron with rescattering off proton-s 
*
*                            BUT NOW
*
* irs = -22 (pp->pp rescattering WITHOUT coulomb's effects)
* irs = -21  (pn -> np rescattering)
*                              and
* irs = -121 - np case - knocked-out neutron with rescattering off proton-r  
* irs = -122 - np case - knocked-out neutron with rescattering off proton-s 
*             (np -> pn rescattering)
*
*                         output parameters
* f_Mre,f_Mim         --  real and imaginary M amplitudes for given spin
* B                   --  Bystricky amplitudes B(1-5)`re, B(6-10)`im
*---------------------------------------------------------------------------
       IMPLICIT REAL(k-n)
       dimension n(3),m(3),l(3)
       dimension B(10), Bpn(10), Bpp(10),B_ampl(10)
       common/par/pi,pm,pmp,pmn,tm,eb  
       common/Bystr/ ibys
       common/fM/  B,plab_old,acm_old,ich_old ! common only for this subr.
       if(ibys.ne.0) ibys =1
*~~~~~~~~~~~~~~~~~~~~~~~~~KINEMATICS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    We assume m1=m2

       plab = sqrt(s**2-4.*s*pm**2)/2./pm
       argum = (2.*t+s-4.*pm**2)/(s-4.*pm**2)
      
       if(argum.gt.1.)   argum =  1.0
       if(argum.lt.-1.0) argum = -1.0
       acm_rad = acos(argum)
       acm  = (acm_rad/pi)*180.
       ich = 1+ibys+irs
*~~~~~~~Obtains Bystricky amplitudes for given irs, see (2.4) from[1]~~~~~~

*      B(1)-B(5) are a_re, b_re, c_re, d_re, e_re
*      B(6)-B(10) are a_im, b_im, c_im, d_im, e_im
*   a-e are Bystricky amplitudes (see [1]), with normalization as in T2.1
*
*   if ibys = 1, gets amplitudes through data files, 
*   if ibys = 0, gets amplitudes through direct calculation (said's "nnsola")
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       if(ibys.eq.0) then
*                          Direct calculation
         if(plab.ne.plab_old.or.acm.ne.acm_old.or.ich.ne.ich_old)then
*avoids to run Bystricky subroutine if only spins 've been changed
           if (abs(irs).eq.22) then    ! pp --> pp case
              if(irs.eq.22) then
                IR = 0      ! Corresponding to M(pp=>pp) = M1
                CALL Bystricky(plab,acm,IR,B)
              else          ! irs.eq.-22
                IR = 20     ! Corresponding to no_coulomb M(pp=>pp) = M1
                CALL Bystricky(plab,acm,IR,B)
              endif
           elseif (irs.gt.0)then   !pn case
             IR = 1     ! Corresponding to M(np=>np)=M(pn=>pn) = (M1+M0)/2
             CALL Bystricky(plab,acm,IR,B)
           elseif (irs.lt.0)then
* CHARGE EXCHANGE: M(np=>pn)=M(pn=>np) = M(pp=>pp)- M(pn=>pn) = (M1-M0)/2
             IR = 0
             CALL Bystricky(plab,acm,IR,Bpp)
             IR = 1
             CALL Bystricky(plab,acm,IR,Bpn)
             B(1) = -Bpp(1) + Bpn(1); B(6) = -Bpp(6) + Bpn(6)
             B(2) = -Bpp(3) + Bpn(3); B(7) = -Bpp(8) + Bpn(8)
             B(3) = -Bpp(2) + Bpn(2); B(8) = -Bpp(7) + Bpn(7)
             B(4) =  Bpp(4) - Bpn(4); B(9) =  Bpp(9) - Bpn(9)
             B(5) =  Bpp(5) - Bpn(5); B(10)= Bpp(10) - Bpn(10)
             acm_rad = pi - acm_rad
           endif 
*             sum = 0.0
*             sum=dot_product(B,B)
*             write(*,*) 'Output'
*             write(*,*) irs,plab*1000.,acm,'DSG(before polar.)=',sum*0.5 
         endif
       
       else ! ibys = 1
*                          Calculation using data files
         if(plab.ne.plab_old.or.acm.ne.acm_old.or.ich.ne.ich_old)then 
*avoids to run Bystricky_dat subroutine if only spins 've been changed 
           if (abs(irs).eq.22) then
              if(irs.eq.22) then
                IR = 0      ! Corresponding to M(pp=>pp) = M1
                CALL Bystricky_dat(plab,acm,IR,B)
              else          ! irs.eq.-22
                IR = 20     ! Corresponding to no_coulomb M(pp=>pp) = M1
                CALL Bystricky_dat(plab,acm,IR,B)
              endif
           elseif (irs.gt.0)then
             IR = 1     ! Corresponding to M(np=>np)=M(pn=>pn) = (M1+M0)/2
             CALL Bystricky_dat(plab,acm,IR,B)
           elseif (irs.lt.0)then
* CHARGE EXCHANGE: M(np=>pn)=M(pn=>np) = M(pp=>pp)- M(pn=>pn) = (M1-M0)/2
             IR = 20
             CALL Bystricky_dat(plab,acm,IR,Bpp)
             IR = 1
             CALL Bystricky_dat(plab,acm,IR,Bpn)
             B(1) = -Bpp(1) + Bpn(1); B(6) = -Bpp(6) + Bpn(6)
             B(2) = -Bpp(3) + Bpn(3); B(7) = -Bpp(8) + Bpn(8)
             B(3) = -Bpp(2) + Bpn(2); B(8) = -Bpp(7) + Bpn(7)
             B(4) =  Bpp(4) - Bpn(4); B(9) =  Bpp(9) - Bpn(9)
             B(5) =  Bpp(5) - Bpn(5); B(10)= Bpp(10) - Bpn(10)
             acm_rad = pi - acm_rad
           endif 

*             sum = 0.0
*             sum=dot_product(B,B)
*             write(*,*) 'Output'
*             write(*,*) irs,plab*1000.,acm,'DSG(before polar.)=',sum*0.5 
         endif
       
       endif
       plab_old = plab
       acm_old  = acm
       ich_old  = ich



*~~~~~~~~~~ Obtains basis C.M. l,m,n vectors  (see [1] (2.2) )~~~~~~~~~
       n = (/       0.,             1.,          0.        /)
       l = (/  sin(acm_rad/2.),     0.,   cos(acm_rad/2.)  /)
       m = (/  cos(acm_rad/2.),     0.,  -sin(acm_rad/2.)  /)
       
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SPINS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*
*      Prepare to find M(k_,k) for given spin configurations
*                There are 16 (f1)(f2)M(i1)(i2) possibilities
*udMdu means that FIRST nucleon's initial spin is down,final spin is up
*                     numerically "u" is 1, "d" is -1 
*          16 spin configurations can be calculated through 5 paterns
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       ifinal = is1p + is2p 
       initial= is1  + is2
       if(abs(ifinal+initial).eq.2)then  
*--------------------------Patern 1-----------------------------------
         ipatern = 1 !(8 elements falls under this pattern out of 16)

*        There are 4 cases(2 identical elements in each) in Patern 1

*        No adjusment needed for    uuMud and uuMdu (ifinal=2,initial=0)(1.1)

         if(abs(initial+is1p).eq.3)then !       udMuu and duMdd (1.3)  
           n(2)=-n(2)
           m(2)=-m(2)
           l(2)=-l(2)
         elseif(abs(initial+is2p).eq.3)then!   udMdd and duMuu (1.4)
           n(3)=-n(3)
           m(3)=-m(3)
           l(3)=-l(3)
         elseif(ifinal.eq.-2)then!    ddMud and ddMdu (1.2)
           n(2)=-n(2)
           m(2)=-m(2)
           l(2)=-l(2)
           n(3)=-n(3)
           m(3)=-m(3)
           l(3)=-l(3)
         endif !End Patern 1 
*--------------------------Patern 2-----------------------------------         
       elseif(abs(ifinal+initial).eq.4)then

         ipatern = 2 !
*      No adjusment needed for        uuMuu (ifinal=2,initial=2)
         if(initial.eq.-2)then!       ddMdd 
           n(3)=-n(3)
           m(3)=-m(3)
           l(3)=-l(3)
         endif !End Patern 2
*--------------------------Patern 3----------------------------------- 
         elseif(abs(ifinal-initial).eq.4)then

         ipatern = 3 !
*      No adjusment needed for        uuMdd (ifinal=2,initial=-2)
         if(initial.eq.2)then!        ddMuu 
           n(2)=-n(2)
           m(2)=-m(2)
           l(2)=-l(2)
         endif !End Patern 3
*--------------------------Paterns 4 and 5---------------------------- 
       elseif(ifinal.eq.0.and.initial.eq.0)then
         if(is1p.eq.1.and.is2p.eq.-1.and.is1.eq.1.and.is2.eq.-1)
     &   ipatern=4!                   udMud
         if(is1p.eq.-1.and.is2p.eq.1.and.is1.eq.-1.and.is2.eq.1)
     &   ipatern=4!                   duMdu
         if(is1p.eq.-1.and.is2p.eq.1.and.is1.eq.1.and.is2.eq.-1)
     &   ipatern=5!                   duMud
         if(is1p.eq.1.and.is2p.eq.-1.and.is1.eq.-1.and.is2.eq.1)
     &   ipatern=5!                   udMdu
       ! End Paterns 4 and 5
*---------------------------------------------------------------------
       else
         write(6,*)" Bystricky amplitudes: ERROR: "
         write(6,*)"fails to chose relevant patern for matrix elements"
         stop
       endif
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END SPINS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       
*                        BEGIN CALCULATIONS

       select case(ipatern)
!      Patern expressions for scattering amplitudes M
!      see T.'s "write-up"s as of 05-01-2003 

*------------------- Patern 1 (uuMud)---------------------------------
       case(1)
            Mre = (B(1)-B(2))*n(1)*n(3) + (B(6)-B(7))*n(2)*n(3)
     &           +(B(3)+B(4))*m(1)*m(3) + (B(8)+B(9))*m(2)*m(3)
     &           +(B(3)-B(4))*l(1)*l(3) + (B(8)-B(9))*l(2)*l(3)
     &           + B(5)      *n(1)      +  B(10)     *n(2)

            Mim =-(B(1)-B(2))*n(2)*n(3) + (B(6)-B(7))*n(1)*n(3)
     &           -(B(3)+B(4))*m(2)*m(3) + (B(8)+B(9))*m(1)*m(3)
     &           -(B(3)-B(4))*l(2)*l(3) + (B(8)-B(9))*l(1)*l(3)
     &           - B(5)      *n(2)      +  B(10)     *n(1) 

*------------------- Patern 2 (uuMuu)---------------------------------
       case(2)
            Mre = (B(1)+B(2))           + (B(1)-B(2))*n(3)**2
     &           +(B(3)+B(4))*m(3)**2   + (B(3)-B(4))*l(3)**2
     &           +   2.*B(5) *n(3)
       
            Mim = (B(6)+B(7))           + (B(6)-B(7))*n(3)**2
     &           +(B(8)+B(9))*m(3)**2   + (B(8)-B(9))*l(3)**2
     &           +   2.*B(10)*n(3)
 
*------------------- Patern 3 (uuMdd)---------------------------------
       case(3)
            Mre = (B(1)-B(2))*(n(1)**2-n(2)**2)+2.*(B(6)-B(7))*n(1)*n(2)
     &           +(B(3)+B(4))*(m(1)**2-m(2)**2)+2.*(B(8)+B(9))*m(1)*m(2)
     &           +(B(3)-B(4))*(l(1)**2-l(2)**2)+2.*(B(8)-B(9))*l(1)*l(2)


            Mim =-2.*(B(1)-B(2))*n(1)*n(2)+(B(6)-B(7))*(n(1)**2-n(2)**2)
     &           -2.*(B(3)+B(4))*m(1)*m(2)+(B(8)+B(9))*(m(1)**2-m(2)**2)
     &           -2.*(B(3)-B(4))*l(1)*l(2)+(B(8)-B(9))*(l(1)**2-l(2)**2)

*------------------- Patern 4 (udMud)---------------------------------
       case(4)
            Mre = (B(1)+B(2))           - (B(1)-B(2))*n(3)**2
     &           -(B(3)+B(4))*m(3)**2   - (B(3)-B(4))*l(3)**2

            Mim = (B(6)+B(7))           - (B(6)-B(7))*n(3)**2
     &           -(B(8)+B(9))*m(3)**2   - (B(8)-B(9))*l(3)**2

*------------------- Patern 5 (duMud)---------------------------------
       case(5)
            Mre = (B(1)-B(2))*(n(1)**2+n(2)**2)
     &           +(B(3)+B(4))*(m(1)**2+m(2)**2)
     &           +(B(3)-B(4))*(l(1)**2+l(2)**2)

            Mim = (B(6)-B(7))*(n(1)**2+n(2)**2)
     &           +(B(8)+B(9))*(m(1)**2+m(2)**2)
     &           +(B(8)-B(9))*(l(1)**2+l(2)**2)

*--------------------------case default----------------------
       case default
            write(6,*)" Bystricky amplitudes: ERROR: "
            write(6,*)ipatern,"is invalid value for ipatern"
            write(6,*) "ipatern = 1,2,3,4,5"
            stop

       end select
       f_Mre = Mre/2.
       f_Mim = Mim/2.
       B_ampl = B
       return
       end

      



       subroutine Bystricky_dat(p_lab,a_cm,IR,B)
*-------------------------------------------------------------
* Obtains Bystricky a,b,c,d,e amplitudes using      
* B_ampl_pp.dat and B_ampl_pn.dat data files
* For more information see main program for amplitudes 
*
*                          ATTENTION! 
*                        B_ampl_pp.dat 
*         data file takes into account Coulomb's effects
*
*input
*  p_lab  -- lab momentum in GeV
*  acm    -- cm angle in deg
*  IR     -- pp/pn flag IR=0 -> pp
*                       IR=1 -> pn
*                       IR<0 -> initialisation
*output 
*  B      -- Bystricky amplitudes in sqrt(mb)
*            B(1)-B(5) are REAL a-e
*            B(6)-B(10)are IMAGINARY a-e
*
* May 20, 2003
*-------------------------------------------------------------

      parameter(nplab = 400, nacm = 90)
      parameter(h_p= 10., h_a=2.)
      dimension B_pp(10,nplab,0:nacm),B_pn(10,nplab,0:nacm)
      dimension B_pp_nc(10,nplab,0:nacm)
      dimension TR(5,8),TI(5,8),H(10)
      dimension B00(10),B01(10),B10(10),B11(10),B(10)
      common/B_NN/ B_pp,B_pn,B_pp_nc

      B = 0.0
     
*     GeV to MeV
      plab = 1000.*p_lab 
      acm  = a_cm
*                          initialization
      if (IR.lt.0) then
*      print*, ' '
*      print*, '         Subroutine Bystricky_dat: Reading data files'
*      print*, ' '
      open(22,file='B_ampl_pp.dat', status='unknown')
      do J=1,10
      read(22,*) ((B_pp(J,k,i),i=0,nacm),k=1,nplab)
      enddo
      close(22)
      open(23,file='B_ampl_pn.dat', status='unknown')
      do J=1,10
      read(23,*) ((B_pn(J,k,i),i=0,nacm),k=1,nplab)
      enddo
      close(23)
      open(24,file='B_ampl_pp_nocoul.dat', status='unknown')
      do J=1,10
      read(24,*) ((B_pp_nc(J,k,i),i=0,nacm),k=1,nplab)
      enddo
      close(24)
      return
      endif
      
*                           gets amplitudes
*              by using first order interpolating polinom

*     Finding the the mesh points plab_k and acm_i closest to plab and acm
*     whose value does not exceed plab and acm 

      k = int(plab/10.)
      i = int(acm/2.)

      
*~~~~~~~~~~~~~~~~~~~~~~~~~~Boundary cases~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(k.lt.1) k = 1
      if(plab.ge.3995.) k=nplab-1  
      fk = float(k)
      fi = float(i)
    
      if(acm.ge.179.) then
      acm = 180.
      i  = nacm
      fi = float(i)
        if (IR.eq.0) then     ! PP
          do J =1,10
          B00(J) = B_pp(J,k,i)
          B10(J) = B_pp(J,k+1,i)
          enddo
        elseif(IR.eq.20) then ! PP nocoul
          do J =1,10
          B00(J) = B_pp_nc(J,k,i)
          B10(J) = B_pp_nc(J,k+1,i)
          enddo
        elseif(IR.eq.1) then  ! NP
          do J =1,10
          B00(J) = B_pn(J,k,i)
          B10(J) = B_pn(J,k+1,i)
          enddo
        endif
        B = B00*( (h_p*(fk+1.)-plab)/h_p * (h_a*(fi+1.)-acm)/h_a )
     &     +B10*( (plab-h_p*fk)/h_p      * (h_a*(fi+1.)-acm)/h_a )
      return
      endif

      if(k*i.lt.0) then
      write(6,*) " Bystricky amplitudes: ERROR: "
      write(6,*) " invalid plab or acm value "
      write(6,*) " 10<=plab<=4000 MeV ; 0<=acm<=180 deg "
      return
      stop
      endif
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*                       discriminates against pp/pn reactions 
      if (IR.eq.0) then     ! PP
      do J =1,10
      B00(J) = B_pp(J,k,i)
      B10(J) = B_pp(J,k+1,i)
      B01(J) = B_pp(J,k,i+1)
      B11(J) = B_pp(J,k+1,i+1)
      enddo
      elseif(IR.eq.20) then ! PP nocoul
      do J =1,10
      B00(J) = B_pp_nc(J,k,i)
      B10(J) = B_pp_nc(J,k+1,i)
      B01(J) = B_pp_nc(J,k,i+1)
      B11(J) = B_pp_nc(J,k+1,i+1)
      enddo
      elseif(IR.eq.1) then  ! NP
      do J =1,10
      B00(J) = B_pn(J,k,i)
      B10(J) = B_pn(J,k+1,i)
      B01(J) = B_pn(J,k,i+1)
      B11(J) = B_pn(J,k+1,i+1)
      enddo
      endif
*                 interpolating polinom in Lagrange representation

      B = B00*( (h_p*(fk+1.)-plab)/h_p * (h_a*(fi+1.)-acm)/h_a )
     &   +B10*( (plab-h_p*fk)/h_p      * (h_a*(fi+1.)-acm)/h_a )
     &   +B01*( (h_p*(fk+1.)-plab)/h_p * (acm-h_a*fi)/h_a      )
     &   +B11*( (plab-h_p*fk)/h_p      * (acm-h_a*fi)/h_a      )

      return
      end

*----------------- Bystricky amplitudes via NNSOLA ----------------------------

      subroutine Bystricky(plab,acm,IR,B) 
*     Bystricky amplitudes in sqrt(mb)
      common/par/pi,pm,pmp,pmn,tm,eb
      dimension H(10),B(10)      
      DIMENSION TR(5,8),TI(5,8)
      H = 0.0
      tlab = sqrt(plab**2+pm**2)-pm
* GeV to MeV
      pmass = 1000.*pm
      tlab = 1000.*tlab 
*      CALL NNSOLA(tlab,acm,IR,TR,TI,H,TTL)
      pcm = sqrt(pmass*tlab/20.)/197.327 ! C. M. momentum in sqrt(mb)
      z=cos(acm*3.1419/180.0)
      s=sin(acm*3.1419/180.0)
      B(1)=2.*H(4)*s/pcm + z*(H(3)+H(5))/pcm
      B(2)=(H(2)+H(1))/pcm
      B(3)=(H(2)-H(1))/pcm
      B(4)=(H(5)-H(3))/pcm
      B(5)=s*(H(8)+H(10))/pcm - 2.*H(9)*z/pcm
      B(6)=2.*H(9)*s/pcm + z*(H(8)+H(10))/pcm
      B(7)=(H(7)+H(6))/pcm
      B(8)=(H(7)-H(6))/pcm
      B(9)=(H(10)-H(8))/pcm
      B(10)=-s*(H(3)+H(5))/pcm + 2.*H(4)*z/pcm

      return
      END
      
      function Trace(A)
      dimension A(4,4)
      sum = 0.0
      do i=1,4
      sum = sum+A(i,i)
      enddo
      Trace = sum
      end 
