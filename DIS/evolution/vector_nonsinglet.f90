module vector_nonsinglet

  use kinds
  use msbar
  use utils

  implicit none
  private

  public :: evolve_NS_step, &
      evolve_NS_DGLAP, evolve_NS_ERBL, evolve_NS_antiDGLAP

  contains



    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Evolve NS by a small step (public routine)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    subroutine evolve_NS_step(nx, xx, HNS, Q2, dt, xi)
        ! Expects positive xi ... one can pass it abs(xi) from a routine
        ! that deals with negative xi
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in) :: xx
        real(dp), dimension(nx), intent(inout) :: HNS
        real(dp), intent(in) :: Q2, dt, xi

        real(dp), dimension(:), allocatable :: dHNS
        integer :: ib1, ib2

        allocate(dHNS(nx))
        dHNS = 0.

        ! Find boundaries of ERBL region
        ib1 = locate(xx, -xi)
        ib2 = locate(xx,  xi)

        if(minval(xx) .lt. 0.0_dp) then
          call evolve_NS_antiDGLAP(ib1, xx(1:ib1), HNS(1:ib1), dHNS(1:ib1), xi)
        endif

        if(xi.ne.0.0_dp) then
          call evolve_NS_ERBL(nx, ib1, ib2, xx, HNS, dHNS(ib1+1:ib2), xi)
        end if

        if(maxval(xx) .gt. 0.0_dp) then
          call evolve_NS_DGLAP(nx-ib2, xx(ib2+1:nx), HNS(ib2+1:nx), dHNS(ib2+1:nx), xi)
        endif

        ! Weight by alpha and add the increments
        dHNS = dHNS * get_alpha_QCD(Q2) / (2.*pi) * dt
        HNS = HNS + dHNS

        deallocate(dHNS)
    end subroutine evolve_NS_step



    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Routines to evolve three regions
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    subroutine evolve_NS_antiDGLAP(nx, xx, HNS, dHNS, xi)
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in) :: xx, HNS
        real(dp), dimension(nx), intent(out) :: dHNS
        real(dp), intent(in) :: xi

        real(dp), dimension(:,:), allocatable :: intd
        integer :: ix, iy

        allocate(intd(nx,nx))
        intd = 0.

        forall(ix=1:nx)
            forall(iy=1:ix)
                intd(ix,iy) = P_NS_plus(xx(ix), xx(iy), xi)*(HNS(iy)-HNS(ix))
            end forall
            dHNS(ix) = trapezoid(nx, xx, intd(ix,:))
            dHNS(ix) = dHNS(ix) + P_NS_delta(xx(ix), xi)*HNS(ix)
        end forall

        dHNS = dHNS * CF
        deallocate(intd)
    end subroutine evolve_NS_antiDGLAP



    subroutine evolve_NS_ERBL(nx, ib1, ib2, xx, HNS, dHNS, xi)
        integer, intent(in) :: nx, ib1, ib2
        real(dp), dimension(nx), intent(in) :: xx, HNS
        real(dp), dimension(ib1+1:ib2), intent(out) :: dHNS
        real(dp), intent(in) :: xi

        real(dp), dimension(:,:), allocatable :: intd
        integer :: ix, iy

        allocate(intd(nx,nx))
        intd = 0.

        ! Outer loop is for generating increments; only applies to ERBL region.
        forall(ix=ib1+1:ib2)
            forall(iy=1:nx)
                intd(ix,iy) = P_NS_plus(xx(ix), xx(iy), xi)*(HNS(iy)-HNS(ix))
            end forall
            dHNS(ix) = P_NS_delta(xx(ix), xi) * HNS(ix)
            dHNS(ix) = dHNS(ix) + trapezoid(nx, xx, intd(ix,:))
        end forall

        dHNS = dHNS * CF
        deallocate(intd)
    end subroutine evolve_NS_ERBL



    subroutine evolve_NS_DGLAP(nx, xx, HNS, dHNS, xi)
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in) :: xx, HNS
        real(dp), dimension(nx), intent(out) :: dHNS
        real(dp), intent(in) :: xi

        real(dp), dimension(:,:), allocatable :: intd
        integer :: ix, iy

        allocate(intd(nx,nx))
        intd = 0.

        forall(ix=1:nx)
            forall(iy=ix:nx)
                intd(ix,iy) = P_NS_plus(xx(ix), xx(iy), xi)*(HNS(iy)-HNS(ix))
            end forall
            dHNS(ix) = trapezoid(nx, xx, intd(ix,:))
            dHNS(ix) = dHNS(ix) + P_NS_delta(xx(ix), xi)*HNS(ix)
        end forall

        dHNS = dHNS * CF
        deallocate(intd)
    end subroutine evolve_NS_DGLAP



    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Splitting functions
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    pure function P_NS_delta(x, xi) result(res)
        real(dp), intent(in) :: x, xi
        real(dp) :: res

        ! Return 0 for pathological x values
        if(tooclose(x,1.0_dp) .or. tooclose(x,-1.0_dp) .or. tooclose(x,0.0_dp)) then
          res = 0.
          return
        end if

        ! Pure DGLAP case
        if(tooclose(xi,0.0_dp)) then
          res = 0.5 + abs(x) + 2.*log(1.-abs(x)) - log(abs(x))
          return
        end if

        ! DGLAP-ERBL boundary limits
        if(tooclose(x,xi) .or. tooclose(x,-xi)) then
          res = 1.5 + log(1.-xi) - log(2.*xi)
          return
        end if

        ! ERBL region
        if(abs(x).lt.xi) then
          res = 1.5 + log( (1.-x**2)/(1.+xi) ) - 0.5*log(xi**2-x**2) &
              + 0.5*x/xi * log( (xi-x)/(xi+x) )
          return
        end if

        ! DGLAP, safely away from singular cases
        res = 1.5 + 2.*log(1.-abs(x)) - 0.5*log( (1.-xi**2)*(x**2-xi**2) ) &
            - 0.5*abs(x)/xi*log( ((1.-xi)*(abs(x)+xi))/((1.+xi)*(abs(x)-xi)) )
    end function P_NS_delta



    pure function P_NS_plus(x, y, xi) result(res)
        real(dp), intent(in) :: x, y, xi
        real(dp) :: res

        ! Avoid singularities
        if(tooclose(y,xi) .or. tooclose(y,-xi) .or. tooclose(y,x)) then
          res = 0.
          return
        end if

        ! ERBL (contributions from both DGLAP regions in one integrand)
        if( abs(x).lt.xi ) then
          if(y.gt.0.0_dp) then
            res = (x+xi)*(y-x+2.*xi) / ( 2.*xi*(y+xi)*(y-x) )
          else
            res = (x-xi)*(y-x-2.*xi) / ( 2.*xi*(y-xi)*(y-x) )
          end if
          return
        end if

        ! DGLAP, safely away from singular cases
        res = ( x**2 + y**2 - 2.*xi**2 ) / ( (y-x)*(y**2-xi**2) )
        if(x.lt.0.0_dp) res = -res ! accounts for x < 0 DGLAP
    end function P_NS_plus



end module vector_nonsinglet
