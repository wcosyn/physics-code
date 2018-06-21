module evolution_speedy
  ! Herein is the minimal amount of needed evolution code
  ! to evolve F2, with the gluons neglected.

  use kinds
  use msbar
  USE ISO_C_BINDING

  implicit none
  private

  public :: evolve_NS_speedy

  contains



    subroutine evolve_NS_speedy(nx, xx, HNS, Q2i, Q2f, oA, oQ2step) bind(C, name="evolve_ns_speedy_adam")
        ! Evolves an array of non-singlet GPDs from Q2=Q2i to Q2=Q2f
        ! Input:
        !   - nx  : integer
        !           number of points in array
        !   - xx  : real(dp), dimension(nx)
        !           x values
        !   - Q2i : real(dp)
        !           initial Q2 value, GeV**2
        !   - Q2f : real(dp)
        !           final (target) Q2 value, GeV**2
        !
        ! Optional input:
        !   - oA      : real(dp)
        !               nuclear mass number (maximum value of x)
        !   - oQ2step : real(dp)
        !               factor to multiply Q2 by in each step
        !
        ! In/out:
        !   - HNS : real(dp), dimension(nx)
        !           array of PDF values at the arguments xx

        ! I/O
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in) :: xx
        real(dp), dimension(nx), intent(inout) :: HNS
        real(dp), intent(in) :: Q2i, Q2f

        ! Optional I/O
        real(dp), intent(in), optional :: oA, oQ2step
        real(dp) :: A, Q2step

        ! Internal variables
        real(dp) :: Q2, Q2next, dt

        ! Defaults/processing of optional input
        A = 1.
        if(present(oA)) A = oA

        Q2step = 1.2
        if(present(oQ2step)) Q2step = oQ2step

        ! Ensure consistency between step size and evolution direction
        if(Q2step.eq.1.0_dp .or. Q2step.le.0.0_dp) stop ! don't feed the trolls
        if(Q2f.gt.Q2i .and. Q2step.lt.1.0_dp) Q2step = 1./Q2step
        if(Q2f.lt.Q2i .and. Q2step.gt.1.0_dp) Q2step = 1./Q2step

        Q2 = Q2i
        !write(*,*) Q2i, Q2f, A, Q2step
        do
          Q2next = Q2 * Q2step
          ! Avoid overshooting target
          if(comp(Q2next,Q2f,Q2step)) Q2next = Q2f
          dt = log(Q2next/Q2)
          call evolve_NS_step(nx, xx/A, HNS, Q2, dt)
          Q2 = Q2next
          if(compeq(Q2,Q2f,Q2step)) return
        end do

    end subroutine evolve_NS_speedy



    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Evolve NS by a small step
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    subroutine evolve_NS_step(nx, xx, HNS, Q2, dt)
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in) :: xx
        real(dp), dimension(nx), intent(inout) :: HNS
        real(dp), intent(in) :: Q2, dt

        real(dp), dimension(:), allocatable :: dHNS

        allocate(dHNS(nx))
        dHNS = 0.

        call evolve_NS_DGLAP(nx, xx(1:nx), HNS(1:nx), dHNS(1:nx))

        ! Weight by alpha and add the increments
        dHNS = dHNS * get_alpha_QCD(Q2) / (2.*pi) * dt
        HNS = HNS + dHNS

        deallocate(dHNS)
    end subroutine evolve_NS_step



    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Routine to evolve the only region for this code
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    subroutine evolve_NS_DGLAP(nx, xx, HNS, dHNS)
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in) :: xx, HNS
        real(dp), dimension(nx), intent(out) :: dHNS

        real(dp), dimension(:,:), allocatable :: intd
        integer :: ix, iy

        allocate(intd(nx,nx))
        intd = 0.

        forall(ix=1:nx)
            forall(iy=ix:nx)
                intd(ix,iy) = P_NS_plus(xx(ix), xx(iy))*(HNS(iy)-HNS(ix))
            end forall
            dHNS(ix) = trapezoid(nx, xx, intd(ix,:))
            dHNS(ix) = dHNS(ix) + P_NS_delta(xx(ix))*HNS(ix)
        end forall

        dHNS = dHNS * CF
        deallocate(intd)
    end subroutine evolve_NS_DGLAP



    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Splitting functions
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    pure function P_NS_delta(x) result(res)
        real(dp), intent(in) :: x
        real(dp) :: res

        ! Return 0 for pathological x values
        if(tooclose(x,1.0_dp) .or. tooclose(x,-1.0_dp) .or. tooclose(x,0.0_dp)) then
          res = 0.
          return
        end if

        ! Pure DGLAP case ... which is the only case here.
        res = 0.5 + abs(x) + 2.*log(1.-abs(x)) - log(abs(x))
    end function P_NS_delta



    pure function P_NS_plus(x, y) result(res)
        real(dp), intent(in) :: x, y
        real(dp) :: res

        ! Avoid singularities
        if(tooclose(y,x)) then
          res = 0.
          return
        end if

        ! DGLAP, safely away from singular cases
        res = ( x**2 + y**2 ) / ( (y-x)*(y**2) )
    end function P_NS_plus



    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Some utilities
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    pure function trapezoid(nx, xx, yy) result(integral)
        ! Trapezoidal integration of y(x)
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in) :: xx, yy
        real(dp) :: integral
        integer :: ix
        integral = 0.0_dp
        do ix=1, nx-1, 1
          integral = integral + 0.5_dp*( yy(ix+1) + yy(ix) )*( xx(ix+1) - xx(ix) )
        end do
    end function trapezoid



    pure function comp(x1, x2, dir) result(res)
        ! Compares x1 and x2, with comparison depending on dir
        real(dp), intent(in) :: x1, x2, dir
        logical :: res
        if(dir.gt.1.0_dp) then
          res = (x1.gt.x2)
        else
          res = (x1.lt.x2)
        end if
    end function comp



    pure function compeq(x1, x2, dir) result(res)
        ! Compares x1 and x2, with comparison depending on dir
        real(dp), intent(in) :: x1, x2, dir
        logical :: res
        if(dir.gt.1.0_dp) then
          res = (x1.ge.x2)
        else
          res = (x1.le.x2)
        end if
    end function compeq



    pure function tooclose(x, y) result(res)
        ! Tell us whether x and y are "too close."
        real(dp), intent(in) :: x, y
        logical :: res
        real(dp), parameter :: eps_ = 1e-6_dp
        res = .false.
        if( abs(x-y) < eps_) res = .true.
        return
    end function tooclose



end module evolution_speedy
