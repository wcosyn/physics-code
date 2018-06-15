module utils

  use kinds, only: dp

  implicit none
  public

  contains



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



    pure function locate(xx,x) result(iLoc)
        ! Finds an index iLoc such that xx(iLoc) <= x < xx(iLoc+1)
        ! Assumes the sequence is ordered in ascending order.
        real(dp), intent(in) :: xx(:), x
        integer :: iLoc
        integer :: N, iStart, iMid, iEnd

        ! First, determine the array size
        N = size(xx)

        ! If x is too small or too large, then return an endpoint
        if( x <= xx(1) ) then
          iLoc = 1
          return
        else if( x >= xx(n) ) then
          iLoc = n - 1
          return
        end if

        ! Define initial endpoints for our search
        iStart = 0
        iEnd = N+1

        ! Search until the start and midpoints come together
        do while( iEnd-iStart > 1 )
          ! Find the midpoint of current start and end
          iMid = (iStart + iEnd) / 2
          ! If x is above the midpoint, move the start there. Else move the end.
          if( x >= xx(iMid) ) then
            iStart = iMid
          else
            iEnd = iMid
          end if
        end do

        ! The desired location is now the narrowed start point
        iLoc = iStart
    end function locate



end module utils
