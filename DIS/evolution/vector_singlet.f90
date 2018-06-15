module vector_singlet

  use kinds
  use msbar
  use utils
  use vector_nonsinglet

  implicit none
  private

  public :: evolve_singlet_step

  contains



    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Evolve singlet by a small step (public routine)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    subroutine evolve_singlet_step(nx, xx, Sigma, Gluons, Q2, dt, xi)
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in) :: xx
        real(dp), dimension(nx), intent(inout) :: Sigma, Gluons
        real(dp), intent(in) :: Q2, dt, xi

        real(dp), dimension(:), allocatable :: dSigma, dGluons, dSigmaQ, dSigmaG, dGluonsQ, dGluonsG
        integer :: ib1, ib2, nFl

        allocate(dSigma(nx), dGluons(nx), dSigmaQ(nx), dSigmaG(nx), dGluonsQ(nx), dGluonsG(nx))

        dSigmaQ = 0.
        dSigmaG = 0.
        dGluonsQ = 0.
        dGluonsG = 0.

        ! Find boundaries of ERBL region
        ib1 = locate(xx, -xi)
        ib2 = locate(xx,  xi)

        ! How many flavors? (passed into GG splitting functions
        nFl = get_effective_flavors(Q2)

        ! anti-DLGAP region increments
        if(minval(xx) .lt. 0.0_dp) then
          call evolve_singlet_antiDGLAP_QinQ(ib1, xx(1:ib1), Sigma(1:ib1),  dSigmaQ(1:ib1),  xi)
          call evolve_singlet_antiDGLAP_QinG(ib1, xx(1:ib1), Gluons(1:ib1), dSigmaG(1:ib1),  xi)
          call evolve_singlet_antiDGLAP_GinQ(ib1, xx(1:ib1), Sigma(1:ib1),  dGluonsQ(1:ib1), xi)
          call evolve_singlet_antiDGLAP_GinG(ib1, xx(1:ib1), Gluons(1:ib1), dGluonsG(1:ib1), xi, nFl)
        endif

        ! ERBL region increments
        if(xi.ne.0.0_dp) then
          call evolve_singlet_ERBL_QinQ(nx, ib1, ib2, xx, Sigma,  dSigmaQ(ib1+1:ib2),  xi)
          call evolve_singlet_ERBL_QinG(nx, ib1, ib2, xx, Gluons, dSigmaG(ib1+1:ib2),  xi)
          call evolve_singlet_ERBL_GinQ(nx, ib1, ib2, xx, Sigma,  dGluonsQ(ib1+1:ib2), xi)
          call evolve_singlet_ERBL_GinG(nx, ib1, ib2, xx, Gluons, dGluonsG(ib1+1:ib2), xi, nFl)
        endif

        ! DLGAP region increments
        if(maxval(xx) .gt. 0.0_dp) then
          call evolve_singlet_DGLAP_QinQ(nx-ib2, xx(ib2+1:nx), Sigma(ib2+1:nx),  dSigmaQ(ib2+1:nx),  xi)
          call evolve_singlet_DGLAP_QinG(nx-ib2, xx(ib2+1:nx), Gluons(ib2+1:nx), dSigmaG(ib2+1:nx),  xi)
          call evolve_singlet_DGLAP_GinQ(nx-ib2, xx(ib2+1:nx), Sigma(ib2+1:nx),  dGluonsQ(ib2+1:nx), xi)
          call evolve_singlet_DGLAP_GinG(nx-ib2, xx(ib2+1:nx), Gluons(ib2+1:nx), dGluonsG(ib2+1:nx), xi, nFl)
        endif

        ! Add contributions
        dSigma = dSigmaQ + 2.*real(nFl)*dSigmaG
        dGluons = dGluonsQ + dGluonsG

        ! Weight by alpha and add the increments
        dSigma = dSigma * get_alpha_QCD(Q2) / (2.*pi) * dt
        dGluons = dGluons * get_alpha_QCD(Q2) / (2.*pi) * dt
        Sigma = Sigma + dSigma
        Gluons = Gluons + dGluons

        deallocate(dSigma, dGluons, dSigmaQ, dSigmaG, dGluonsQ, dGluonsG)
    end subroutine evolve_singlet_step



    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Regions: quark in quark (same as non-singlet)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    subroutine evolve_singlet_antiDGLAP_QinQ(nx, xx, Quark, dQuark, xi)
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in) :: xx, Quark
        real(dp), dimension(nx), intent(out) :: dQuark
        real(dp), intent(in) :: xi
        call evolve_NS_antiDGLAP(nx, xx, Quark, dQuark, xi)
    end subroutine evolve_singlet_antiDGLAP_QinQ



    subroutine evolve_singlet_ERBL_QinQ(nx, ib1, ib2, xx, Quark, dQuark, xi)
        integer, intent(in) :: nx, ib1, ib2
        real(dp), dimension(nx), intent(in) :: xx, Quark
        real(dp), dimension(ib1+1:ib2), intent(out) :: dQuark
        real(dp), intent(in) :: xi
        call evolve_NS_ERBL(nx, ib1, ib2, xx, Quark, dQuark, xi)
    end subroutine evolve_singlet_ERBL_QinQ



    subroutine evolve_singlet_DGLAP_QinQ(nx, xx, Quark, dQuark, xi)
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in) :: xx, Quark
        real(dp), dimension(nx), intent(out) :: dQuark
        real(dp), intent(in) :: xi
        call evolve_NS_DGLAP(nx, xx, Quark, dQuark, xi)
    end subroutine evolve_singlet_DGLAP_QinQ



    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Regions: quark in gluon
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    subroutine evolve_singlet_antiDGLAP_QinG(nx, xx, Gluon, dQuark, xi)
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in) :: xx, Gluon
        real(dp), dimension(nx), intent(out) :: dQuark
        real(dp), intent(in) :: xi

        real(dp), dimension(:,:), allocatable :: intd
        integer :: ix, iy

        allocate(intd(nx,nx))
        intd = 0.

        forall(ix=1:nx)
            forall(iy=1:ix)
                intd(ix,iy) = P_qG(xx(ix), xx(iy), xi)*Gluon(iy)
            end forall
            dQuark(ix) = trapezoid(nx, xx, intd(ix,:))
        end forall

        dQuark = dQuark * TF
        deallocate(intd)
    end subroutine evolve_singlet_antiDGLAP_QinG



    subroutine evolve_singlet_ERBL_QinG(nx, ib1, ib2, xx, Gluon, dQuark, xi)
        integer, intent(in) :: nx, ib1, ib2
        real(dp), dimension(nx), intent(in) :: xx, Gluon
        real(dp), dimension(ib1+1:ib2), intent(out) :: dQuark
        real(dp), intent(in) :: xi

        real(dp), dimension(:,:), allocatable :: intd
        integer :: ix, iy

        allocate(intd(nx,nx))
        intd = 0.

        forall(ix=ib1+1:ib2)
            forall(iy=1:nx)
                intd(ix,iy) = P_qG(xx(ix), xx(iy), xi)*Gluon(iy)
            end forall
            dQuark(ix) = trapezoid(nx, xx, intd(ix,:))
        end forall

        dQuark = dQuark * TF
        deallocate(intd)
    end subroutine evolve_singlet_ERBL_QinG



    subroutine evolve_singlet_DGLAP_QinG(nx, xx, Gluon, dQuark, xi)
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in) :: xx, Gluon
        real(dp), dimension(nx), intent(out) :: dQuark
        real(dp), intent(in) :: xi

        real(dp), dimension(:,:), allocatable :: intd
        integer :: ix, iy

        allocate(intd(nx,nx))
        intd = 0.

        forall(ix=1:nx)
            forall(iy=ix:nx)
                intd(ix,iy) = P_qG(xx(ix), xx(iy), xi)*Gluon(iy)
            end forall
            dQuark(ix) = trapezoid(nx, xx, intd(ix,:))
        end forall

        dQuark = dQuark * TF
        deallocate(intd)
    end subroutine evolve_singlet_DGLAP_QinG



    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Regions: gluon in quark
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    subroutine evolve_singlet_antiDGLAP_GinQ(nx, xx, Quark, dGluon, xi)
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in) :: xx, Quark
        real(dp), dimension(nx), intent(out) :: dGluon
        real(dp), intent(in) :: xi

        real(dp), dimension(:,:), allocatable :: intd
        integer :: ix, iy

        allocate(intd(nx,nx))
        intd = 0.

        forall(ix=1:nx)
            forall(iy=1:ix)
                !intd(ix,iy) = intdfunc(xx(ix), xx(iy), xi)*Quark(iy)
                intd(ix,iy) = P_Gq(xx(ix), xx(iy), xi)*Quark(iy)
            end forall
            dGluon(ix) = trapezoid(nx, xx, intd(ix,:))
        end forall

        dGluon = dGluon * CF
        deallocate(intd)
    end subroutine evolve_singlet_antiDGLAP_GinQ



    subroutine evolve_singlet_ERBL_GinQ(nx, ib1, ib2, xx, Quark, dGluon, xi)
        integer, intent(in) :: nx, ib1, ib2
        real(dp), dimension(nx), intent(in) :: xx, Quark
        real(dp), dimension(ib1+1:ib2), intent(out) :: dGluon
        real(dp), intent(in) :: xi

        real(dp), dimension(:,:), allocatable :: intd
        integer :: ix, iy

        allocate(intd(nx,nx))
        intd = 0.

        forall(ix=ib1+1:ib2)
            forall(iy=1:nx)
                intd(ix,iy) = P_Gq(xx(ix), xx(iy), xi)*Quark(iy)
            end forall
            dGluon(ix) = trapezoid(nx, xx, intd(ix,:))
        end forall

        dGluon = dGluon * CF
        deallocate(intd)
    end subroutine evolve_singlet_ERBL_GinQ



    subroutine evolve_singlet_DGLAP_GinQ(nx, xx, Quark, dGluon, xi)
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in) :: xx, Quark
        real(dp), dimension(nx), intent(out) :: dGluon
        real(dp), intent(in) :: xi

        real(dp), dimension(:,:), allocatable :: intd
        integer :: ix, iy

        allocate(intd(nx,nx))
        intd = 0.

        forall(ix=1:nx)
            forall(iy=ix:nx)
                intd(ix,iy) = P_Gq(xx(ix), xx(iy), xi)*Quark(iy)
            end forall
            dGluon(ix) = trapezoid(nx, xx, intd(ix,:))
        end forall

        dGluon = dGluon * CF
        deallocate(intd)
    end subroutine evolve_singlet_DGLAP_GinQ



    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Regions: gluon in gluon
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    subroutine evolve_singlet_antiDGLAP_GinG(nx, xx, Gluon, dGluon, xi, nFl)
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in) :: xx, Gluon
        real(dp), dimension(nx), intent(out) :: dGluon
        real(dp), intent(in) :: xi
        integer, intent(in) :: nFl

        real(dp), dimension(:,:), allocatable :: intdR, intdP
        integer :: ix, iy

        allocate(intdR(nx,nx), intdP(nx,nx))
        intdR = 0.
        intdP = 0.

        forall(ix=1:nx)
            forall(iy=1:ix)
                intdR(ix,iy) = P_GG_regular(xx(ix), xx(iy), xi)*Gluon(iy)
                intdP(ix,iy) = P_GG_plus(xx(ix), xx(iy), xi)*(Gluon(iy)-Gluon(ix))
            end forall
            dGluon(ix) = P_GG_delta(xx(ix),xi,nFl)*Gluon(ix)
            dGluon(ix) = dGluon(ix) + trapezoid(nx, xx, intdR(ix,:))
            dGluon(ix) = dGluon(ix) + trapezoid(nx, xx, intdP(ix,:))
        end forall

        dGluon = dGluon * CA
        deallocate(intdR, intdP)
    end subroutine evolve_singlet_antiDGLAP_GinG



    subroutine evolve_singlet_ERBL_GinG(nx, ib1, ib2, xx, Gluon, dGluon, xi, nFl)
        integer, intent(in) :: nx, ib1, ib2
        real(dp), dimension(nx), intent(in) :: xx, Gluon
        real(dp), dimension(ib1+1:ib2), intent(out) :: dGluon
        real(dp), intent(in) :: xi
        integer, intent(in) :: nFl

        real(dp), dimension(:,:), allocatable :: intdR, intdP
        integer :: ix, iy

        allocate(intdR(nx,nx), intdP(nx,nx))
        intdR = 0.
        intdP = 0.

        forall(ix=ib1+1:ib2)
            forall(iy=1:nx)
                intdR(ix,iy) = P_GG_regular(xx(ix), xx(iy), xi)*Gluon(iy)
                intdP(ix,iy) = P_GG_plus(xx(ix), xx(iy), xi)*(Gluon(iy)-Gluon(ix))
            end forall
            dGluon(ix) = P_GG_delta(xx(ix),xi,nFl)*Gluon(ix)
            dGluon(ix) = dGluon(ix) + trapezoid(nx, xx, intdR(ix,:))
            dGluon(ix) = dGluon(ix) + trapezoid(nx, xx, intdP(ix,:))
        end forall

        dGluon = dGluon * CA
        deallocate(intdR, intdP)
    end subroutine evolve_singlet_ERBL_GinG



    subroutine evolve_singlet_DGLAP_GinG(nx, xx, Gluon, dGluon, xi, nFl)
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in) :: xx, Gluon
        real(dp), dimension(nx), intent(out) :: dGluon
        real(dp), intent(in) :: xi
        integer, intent(in) :: nFl

        real(dp), dimension(:,:), allocatable :: intdR, intdP
        integer :: ix, iy

        allocate(intdR(nx,nx), intdP(nx,nx))
        intdR = 0.
        intdP = 0.

        forall(ix=1:nx)
            forall(iy=ix:nx)
                intdR(ix,iy) = P_GG_regular(xx(ix), xx(iy), xi)*Gluon(iy)
                intdP(ix,iy) = P_GG_plus(xx(ix), xx(iy), xi)*(Gluon(iy)-Gluon(ix))
            end forall
            dGluon(ix) = P_GG_delta(xx(ix),xi,nFl)*Gluon(ix)
            dGluon(ix) = dGluon(ix) + trapezoid(nx, xx, intdR(ix,:))
            dGluon(ix) = dGluon(ix) + trapezoid(nx, xx, intdP(ix,:))
        end forall

        dGluon = dGluon * CA
        deallocate(intdR, intdP)
    end subroutine evolve_singlet_DGLAP_GinG



    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Splitting functions
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    pure function P_qG(x, y, xi) result(res)
        real(dp), intent(in) :: x, y, xi
        real(dp) :: res

        ! Avoid singularities
        if(tooclose(y,xi) .or. tooclose(y,-xi)) then
          res = 0.
          return
        end if

        ! ERBL region
        if( abs(x).lt.xi ) then
          if(y.gt.0.0_dp) then
            res = (x+xi)*(y-2.*x+xi) / ( 2.*xi*(y+xi)*(y**2-xi**2) )
          else
            res = (x-xi)*(y-2.*x-xi) / ( 2.*xi*(y-xi)*(y**2-xi**2) )
          end if
          return
        end if

        ! DGLAP, safely away from singularities
        res = ( x**2 + (y-x)**2 - xi**2 ) / (2.*(y**2-xi**2)**2)
        if(x.lt.0.0_dp) res = -res ! accounts for x < 0 DGLAP
    end function P_qG



    pure function P_Gq(x, y, xi) result(res)
        real(dp), intent(in) :: x, y, xi
        real(dp) :: res

        ! Avoid singularities
        if(tooclose(y,xi) .or. tooclose(y,-xi)) then
          res = 0.
          return
        end if

        ! ERBL
        if( abs(x).lt.xi ) then
          if(y.gt.0.0_dp) then
            res = (x+xi)*(2.*y-x+xi) / (2.*xi*(y+xi))
          else
            res = (x-xi)*(2.*y-x-xi) / (2.*xi*(y-xi))
          end if
          return
        end if

        ! DGLAP, safely away from singularities
        res = ( y**2 + (y-x)**2 - xi**2 ) / (y**2-xi**2)
        if(x.lt.0.0_dp) res = -res ! account for x < 0 DGLAP
    end function P_Gq



    pure function P_GG_delta(x, xi, nFl) result(res)
        real(dp), intent(in) :: x, xi
        integer, intent(in) :: nFl
        real(dp) :: res

        ! Avoid pathological x values
        if(tooclose(x,1.0_dp) .or. tooclose(x,-1.0_dp)) then
          res = 0.
          return
        end if

        ! TODO: x-xi boundary...!
        ! Maybe need to take a hint from Vinnikov here.
        if(tooclose(x,xi) .or. tooclose(x,-xi)) then
          res = 0.
          return
        end if

        ! ERBL
        if(abs(x).lt.xi) then
          res = log(1.-x**2) - log(xi**2-x**2) + 11./6. - 2.*real(nFl)/(3.*CA)
          return
        end if

        ! x=0 is OK in ERBL, but is a singularity for DGLAP ... ?
        if(tooclose(x,0.0_dp)) then
          res = 0.
          return
        end if

        ! DGLAP, safely away from singular region
        res = 2.*log(1.+abs(x)) - log(x**2-xi**2) + 11./6. - 2.*real(nFl)/(3.*CA)
    end function P_GG_delta



    pure function P_GG_plus(x, y, xi) result(res)
        real(dp), intent(in) :: x, y, xi
        real(dp) :: res

        ! Avoid singularities
        if(tooclose(y,x)) then
          res = 0.
          return
        end if

        ! ERBL
        if(abs(x).lt.xi) then
          if(y.gt.0.0_dp) then
            res = 2. / (y-x)
          else
            res = -2. / (y-x)
          end if
          return
        end if

        ! DGLAP, away from singularities
        res = 2. / (y-x)
        if(x.lt.0.0_dp) res = -res
    end function P_GG_plus



    pure function P_GG_regular(x, y, xi) result(res)
        real(dp), intent(in) :: x, y, xi
        real(dp) :: res

        ! Avoid singularities
        if(tooclose(y,xi) .or. tooclose(y,-xi)) then
          res = 0.
          return
        end if

        ! TODO: x=xi boundaries
        if(tooclose(x,xi) .or. tooclose(x,-xi)) then
          res = 0.
          return
        end if

        ! ERBL
        if(abs(x).lt.xi) then
          if(y.gt.0.0_dp) then
            res = ( &
                (xi**2-x**2)*(1.-2.*(y**2+x**2)/((y+xi)*(x-xi)))/(2.*xi) - (y+x) &
                ) / (y**2-xi**2)
          else
            res = ( &
                (xi**2-x**2)*(1.-2.*(y**2+x**2)/((y-xi)*(x+xi)))/(2.*xi) + (y+x) &
                ) / (y**2-xi**2)
          end if
          return
        end if

        ! DGLAP, safely away from singularities
        res = 2.*(y-x)*(y**2+x**2)/(y**2-xi**2)**2 - 2.*(y+x)/(y**2-xi**2)
        if(x.lt.0.0_dp) res = -res ! account for x < 0 DGLAP
    end function P_GG_regular



end module vector_singlet
