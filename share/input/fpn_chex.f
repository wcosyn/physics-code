*        program test
*        pi = acos(-1.0)
*        pm    = 0.938187
*        ins = 0
*    	call f_pn_chex(is1,is2,is1p,is2p,s,t,f_re,f_im,ins)
*        thcm  = pi
*        is1 =  1
*        is2 =  -1
*        is1p = 1
*        is2p = -1

*        do ip = 1,25
*        plab = float(ip)*0.4
*        elab = sqrt(pm**2 + plab**2)
*        s    = 2.0*pm**2 + 2.0*pm*elab
*        pcm  = sqrt(s-4.0*pm**2)/2.0
*        ecm  = sqrt(pm**2+pcm**2)



        

*        thcm  = pi
*        t    = -2.0*(pcm**2) * (1.0-cos(thcm))
*        ins = 21 
*    	call f_pn_chex(is1,is2,is1p,is2p,s,t,f_reb,f_imb,ins)

*        thcm  = 0.
*        t    = -2.0*(pcm**2) * (1.0-cos(thcm))
*        ins = 10021 
*    	call  f_pn_chex(is1,is2,is1p,is2p,s,t,f_ref,f_imf,ins)      

        
*         print 12,    plab,f_reb,f_imb,f_ref,f_imf
*        write(24,12) plab,f_reb,f_imb,f_ref,f_imf
*        enddo
* 12     format(5f12.3)
*        end





    	subroutine f_pn_chex(is1,is2,is1p,is2p,s,t,f_re,f_im,ins)
        common/par/pi,pm,pmp,pmn,tm,eb  
******************************************************************************************
*   Subroutine calculates charge exchange pn amplitudes
*
* is1, is1p - initial and final spin projection of particle 1  (1 - up -1 - down)
* is2, is2p - initial and final spin projection of particle 2  (1 - up -1 - down)
* s         - total invariant energy of the system in GeV^2
* t         - invariant transferred momentum as defined for elastic scattering in GeV^2
*
* ins = 0 - initialization
* ins = 21 - calculates using only SAID parameterization, with charge exchange 
*            interpreting as a backward angle scattering of the pn amplitude
*            so in this case t shoud be identified with u (correct till plan=1.4 extended till 3.68)
* ins = 10 - calculates using only Gibbs Loiseau parameterization, t should be identfied with u
*            as above (for all the range of plab
*
* ins = 1021 - combinaation of the above, SAID bellow plab=3.68 and GL above 
*
* ina=10021  - SAID till p0 and imaginary from SAID and real from GL 
*              after p0, withe real part calculated as  SAID*s0/s
*              Since GL is dereved with the assumption that the amplitude is fully real
*              we renormalized it to make room for the imaginary part. (look at ren2)
*
* ins = -21 - charge exchange as defined by the SAID, in this case one enters the same t
*             as defined by elastic scattering
*
********************************************************************************************
        f_re = 0.0
        f_im = 0.0
        if(ins.eq.0)then
        call Bystricky_dat(0.01,0.,-1,BB)
        pi  = acos(-1.0)       
        pm  = 0.938279         
        return
        endif 
       
        if(ins.eq.10)then
         call  f_NN_chex_gl(is1,is2,is1p,is2p,s,t,f_re,f_im)
        elseif(abs(ins).eq.21)then
         ini = 1
    	 call f_NN_sd(is1,is2,is1p,is2p,s,t,f_re,f_im,ins,ini)
        elseif(abs(ins).eq.1021)then
         ep = (s - 2.0*pm**2)/(2.0*pm)
         p = sqrt(ep**2-pm**2)
*         if(p.lt.3.68)then
         if(p.le.1.4)then
         irs = 21
         ini = 1
    	 call f_NN_sd(is1,is2,is1p,is2p,s,t,f_re,f_im,irs,ini)
         else
         call  f_NN_chex_gl(is1,is2,is1p,is2p,s,t,f_re,f_im)
         endif
        elseif(ins.eq.10021)then
         p0 = 3.6
         ep = (s - 2.0*pm**2)/(2.0*pm)
         p = sqrt(ep**2-pm**2)
         if(p.le.p0)then
         ini = 1
         irs = 21
    	  call f_NN_sd(is1,is2,is1p,is2p,s,t,f_re,f_im,irs,ini)
         else
         e0 = sqrt(pm**2 + p0**2)
         s0 = 2.0*pm**2 + 2.0*pm*e0
         ini = 1
         irs = 21
    	 call f_NN_sd(is1,is2,is1p,is2p,s0,t,f_re0,f_im0,irs,ini)
         f_im = f_im0*s0/s
         call  f_NN_chex_gl(is1,is2,is1p,is2p,s,t,f_re1,f_im0)
         ren2 = (f_re1**2-f_im**2)/f_re1**2 ! renormalization factor
         if(ren2.le.0.0)ren2 = 1.0
         f_re = f_re1 * sqrt(ren2) 
         endif
       endif
        return
        end


       subroutine f_NN_chex_gl(is1,is2,is1p,is2p,s,t,f_NN_re,f_NN_im)
******************************************************************************
*             Parameterization of charge exchange amplitudes based on        * 
*                 Gibbs, Loiseau Phys. Rev. C 50, 2742 paper                 *
******************************************************************************
*             spin polarized charge exchange amplitudes in GeV-2
*             Normalization is: Im(f)=sigma_tot (so, 4*pi/pcm sits in f)
*input
* s                  --  (p_1^mu+p_2^mu)^2
* t                  --  (p_1^mu-p_3^mu)^2 
*
* is1,is2,is1_,is2_  -- initial and final spins of scattered nucleons
*
*output 
*f_NN_re           -- real part of chex amplitude
*f_NN_im           -- imaginary part of chex amplitude
* ini                -- ini=0 initializes, everything else -- computes
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        dimension B(10), BB(10) 
        common/par/pi,pm,pmp,pmn,tm,eb  
        f_NN_re = 0.0
        f_NN_im = 0.0

*~~~~~~~~~~~~~~~~~~~~~~~~~~~ Kinematic variables ~~~~~~~~~~~~~~~~~~~    
        ep = (s - 2.0*pm**2)/(2.0*pm)
        p = sqrt(ep**2-pm**2)
        pcm = sqrt(s-4*pm**2)/2.0 ! in GeV
       

*        argum = (2.*t+s-4.*pm**2)/(s-4.*pm**2) 
*        if(argum.gt.1.)   argum =  1.0
*        if(argum.lt.-1.0) argum = -1.0
*        acm_rad = acos(argum)

        sn_thkes2 = -t/4.0/pcm**2
        if(sn_thkes2.gt.1.0)return
        acm_rad = 2.0*asin(sqrt(sn_thkes2))
        acm = acm_rad*180.0/pi


        acm  = (acm_rad/pi)*180.

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Gibbs chex ~~~~~~~~~~~~~~~~~~~~~~~~~
             call gibbs(p,acm,g_p,g_q,f_p,f_q,dN)
             B(1) = dN/3.*(g_p + 2.*g_q)                 ; B(6) =  0.
             B(2) = dN/3.*(-g_p - 6.*f_q*g_q + 4.*g_q)   ; B(7) =  0.
             B(3) = dN/3.*(-3.*f_p*g_p + 2.*g_p - 2.*g_q); B(8) =  0.
             B(4) = dN*(-f_p*g_p + 2.*f_q*g_q)           ; B(9) =  0.
             B(5) =  0.                                  ; B(10)=  0. 
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SPINS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
*        (8 elements falls under this pattern out of 16)
*        There are 4 cases(2 identical elements in each) in Patern 1
         if(abs(initial+is1p).eq.3)then    !   udMuu, duMdd   
         fMre_pn_gibbs = -B(4)*sin(acm_rad)
         elseif(abs(initial+is2p).eq.3)then!   udMdd, duMuu
         fMre_pn_gibbs =  B(4)*sin(acm_rad)
         elseif(ifinal.eq.-2)then!             ddMud, ddMdu 
         fMre_pn_gibbs =  B(4)*sin(acm_rad)
         else                  !               uuMud, uuMdu 
         fMre_pn_gibbs = -B(4)*sin(acm_rad)
         endif !End Patern 1 

*--------------------------Patern 2-----------------------------------         
       elseif(abs(ifinal+initial).eq.4)then! uuMuu, ddMdd

       fMre_pn_gibbs = B(1) + B(2) + B(3) - B(4)*cos(acm_rad)
             
*--------------------------Patern 3----------------------------------- 
       elseif(abs(ifinal-initial).eq.4)then ! uuMdd, ddMuu 

       fMre_pn_gibbs =-B(1) + B(2) + B(3) + B(4)*cos(acm_rad)

*--------------------------Paterns 4 and 5---------------------------- 
       elseif(ifinal.eq.0.and.initial.eq.0)then
         if(is1p.eq.is1)then                ! udMud, duMdu

         fMre_pn_gibbs = B(1) + B(2) - B(3) + B(4)*cos(acm_rad)

         else                               ! udMdu, duMud

         fMre_pn_gibbs = B(1) - B(2) + B(3) + B(4)*cos(acm_rad)

         endif  
       ! End Paterns 4 and 5
*---------------------------------------------------------------------
       else
         write(6,*)" fchecx.f: ERROR: "
         write(6,*)"fails to chose relevant patern for matrix elements"
         stop
       endif
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END SPINS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       
          fMre = (fMre_pn_gibbs)/2.  ! GeV-1
          fMim = 0.0
          f_NN_re = 4.*pi/pcm*fMre  ! in GeV-2
          f_NN_im = 4.*pi/pcm*fMim  ! in GeV-2
        return
        end

      subroutine gibbs(plab,acm,g_p,g_q,f_p,f_q,dN)
      common/par/pi,pm,pmp,pmn,tm,eb  
      dLambda = 0.748!0.785
      dmu = 0.140 ! pion's mass
      f2_pi_over_4pi = 0.079
       dmu2  = sqrt(4.0*pm**2*f2_pi_over_4pi/14.4)      
   
      acm_rad = acm*pi/180.
      Elab = sqrt(plab**2+pm**2)
      s = 2.*pm**2+2.*pm*Elab       
      pcm = pm*plab/sqrt(s)

      p2 = 2.*pcm**2*(1.-cos(acm_rad))
      q2 = 2.*pcm**2*(1.+cos(acm_rad))

      g_p = (dLambda**2/(p2+dLambda**2))**2
      g_q = (dLambda**2/(q2+dLambda**2))**2

      
      f_p = (p2/(p2+dmu2**2))
      f_q = (q2/(q2+dmu2**2))

      dN  = (pm/dmu)**2 *
     &   2./sqrt(s)* f2_pi_over_4pi*((dLambda**2-dmu2**2)/dLambda**2)**2 
    
      return
      end

