module evolvers

  use decouple
  use kinds
  use msbar
  use utils
  use vector_nonsinglet
  use vector_singlet

  implicit none
  private

  public :: evolve_GPDs, evolve_NS, evolve_singlet

  contains



    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Public routines
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    subroutine evolve_GPDs(nx, xx, H_mesh, Q2i, Q2f, xi, oA, oQ2step)
        ! Evolves a mesh of GPDs in the physical basis from Q2i to Q2f
        ! Input:
        !   - nx  : integer
        !           number of points in array
        !   - xx  : real(dp), dimension(nx)
        !           x values
        !   - Q2i : real(dp)
        !           initial Q2 value, GeV**2
        !   - Q2f : real(dp)
        !           final (target) Q2 value, GeV**2
        !   - xi  : real(dp)
        !           skewness
        !
        ! Optional input:
        !   - oA      : real(dp)
        !               nuclear mass number (maximum value of x and xi)
        !   - oQ2step : real(dp)
        !               factor to multiply Q2 by in each step
        !
        ! In/out:
        !   - H_mesh : real(dp), dimension(-6:6,nx)
        !           array of GPD values at the arguments xx

        ! I/O
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in) :: xx
        real(dp), dimension(-6:6,nx), intent(inout) :: H_mesh
        real(dp), intent(in) :: Q2i, Q2f, xi

        ! Optional I/O
        real(dp), intent(in), optional :: oA, oQ2step
        real(dp) :: A, Q2step

        ! Internal variables
        real(dp) :: Q2, eps_
        integer :: nFl, i
        real(dp), dimension(3) :: HQM2

        A = 1.
        if(present(oA)) A = oA

        Q2step = 1.2
        if(present(oQ2step)) Q2step = oQ2step

        ! Ensure consistency between step size and evolution direction
        if(Q2step.eq.1.0_dp .or. Q2step.le.0.0_dp) stop ! don't feed the trolls
        if(Q2f.gt.Q2i .and. Q2step.lt.1.0_dp) Q2step = 1./Q2step
        if(Q2f.lt.Q2i .and. Q2step.gt.1.0_dp) Q2step = 1./Q2step

        ! The heavy quark masses are ordered depending on if evolution is up or down
        if(Q2step.gt.1.0_dp) then
          HQM2 = [ cMass2, bMass2, tMass2 ]
          eps_ = tiny(1.0_dp)
        else
          HQM2 = [ tMass2, bMass2, cMass2 ]
          eps_ = -tiny(1.0_dp)
        end if

        Q2 = Q2i

        ! Loop to be careful around flavor thresholds
        do i=1, 3, 1
          nFl = get_effective_flavors(Q2)

          ! Case 1: target Q2 falls short of next flavor threshold
          if( comp(HQM2(i),Q2f,Q2step) ) then
            call evolve_GPDs_fixed_nFl(nFl, nx, xx, H_mesh, Q2, Q2f, xi, oA=A, oQ2step=Q2step)
            return

          ! Case 2: target Q2 is beyond next flavor threshold
          else
            call evolve_GPDs_fixed_nFl(nFl, nx, xx, H_mesh, Q2, HQM2(i)-eps_, xi, oA=A, oQ2step=Q2step)
            Q2 = HQM2(i) + eps_

          end if

        end do

        ! Case 3: passed all flavor thresholds, still need to reach target
        nFl = get_effective_flavors(Q2)
        call evolve_GPDs_fixed_nFl(nFl, nx, xx, H_mesh, Q2, Q2f, xi, oA=A, oQ2step=Q2step)

    end subroutine evolve_GPDs



    subroutine evolve_NS(nx, xx, HNS, Q2i, Q2f, xi, &
            ! Optional input
            oA, oQ2step)
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
        !   - xi  : real(dp)
        !           skewness
        !
        ! Optional input:
        !   - oA      : real(dp)
        !               nuclear mass number (maximum value of x and xi)
        !   - oQ2step : real(dp)
        !               factor to multiply Q2 by in each step
        !
        ! In/out:
        !   - HNS : real(dp), dimension(nx)
        !           array of GPD values at the arguments xx

        ! I/O
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in) :: xx
        real(dp), dimension(nx), intent(inout) :: HNS
        real(dp), intent(in) :: Q2i, Q2f, xi

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

        do
          Q2next = Q2 * Q2step
          ! Avoid overshooting target
          if(comp(Q2next,Q2f,Q2step)) Q2next = Q2f
          dt = log(Q2next/Q2)
          call evolve_NS_step(nx, xx/A, HNS, Q2, dt, abs(xi/A))
          Q2 = Q2next
          if(compeq(Q2,Q2f,Q2step)) return
        end do

    end subroutine evolve_NS



    subroutine evolve_singlet(nx, xx, Sigma, Gluons, Q2i, Q2f, xi, &
            ! Optional input
            oA, oQ2step)
        ! Evolves the singlet quark and gluon GPDs from Q2=Q2i to Q2=Q2f
        ! Input:
        !   - nx  : integer
        !           number of points in array
        !   - xx  : real(dp), dimension(nx)
        !           x values
        !   - Q2i : real(dp)
        !           initial Q2 value, GeV**2
        !   - Q2f : real(dp)
        !           final (target) Q2 value, GeV**2
        !   - xi  : real(dp)
        !           skewness
        !
        ! Optional input:
        !   - oA      : real(dp)
        !               nuclear mass number (maximum value of x and xi)
        !   - oQ2step : real(dp)
        !               factor to multiply Q2 by in each step
        !
        ! In/out:
        !   - Sigma  : real(dp), dimension(nx)
        !              array of quark singlet GPD values at the arguments xx
        !   - Gluons : real(dp), dimension(nx)
        !              array of gluon GPD values at the arguments xx
        !              IMPORTANT: gluon GPD is defined as in Diehl [2]

        ! I/O
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in) :: xx
        real(dp), dimension(nx), intent(inout) :: Sigma, Gluons
        real(dp), intent(in) :: Q2i, Q2f, xi

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

        do
          Q2next = Q2 * Q2step
          ! Avoid overshooting target
          if(comp(Q2next,Q2f,Q2step)) Q2next = Q2f
          dt = log(Q2next/Q2)
          call evolve_singlet_step(nx, xx/A, Sigma, Gluons, Q2, dt, abs(xi/A))
          Q2 = Q2next
          if(compeq(Q2,Q2f,Q2step)) return
        end do

    end subroutine evolve_singlet



    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Implementation details: fixed-flavor physical basis evolution
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    subroutine evolve_GPDs_fixed_nFl(nFl, nx, xx, H_mesh, Q2i, Q2f, xi, oA, oQ2step)
        ! I/O
        integer, intent(in) :: nFl, nx
        real(dp), dimension(nx), intent(in) :: xx
        real(dp), dimension(-6:6,nx), intent(inout) :: H_mesh
        real(dp), intent(in) :: Q2i, Q2f, xi

        ! Optional I/O
        real(dp), intent(in), optional :: oA, oQ2step
        real(dp) :: A, Q2step

        ! Internal variables: the evolution basis
        real(dp), dimension(:,:), allocatable :: q_minus, T_NS
        real(dp), dimension(:), allocatable :: Gluons, Sigma
        integer :: i

        ! Defaults/processing of optional input
        A = 1.
        if(present(oA)) A = oA
        Q2step = 1.2
        if(present(oQ2step)) Q2step = oQ2step

        ! Go to evolution basis
        allocate(q_minus(1:nFl,nx), T_NS(2:nFl,nx), Gluons(nx), Sigma(nx))
        call decouple_arrays(nx, H_mesh(-nFl:nFl,:), q_minus, T_NS, Sigma, Gluons, nFl)

        ! Evolve in the evolution basis
        !!!$OMP PARALLEL DO
        do i=1, nFl, 1
          call evolve_NS(nx, xx, q_minus(i,:), Q2i, Q2f, xi, oA=A, oQ2step=Q2step)
        end do
        !!!$OMP END PARALLEL DO

        !!!$OMP PARALLEL DO
        do i=2, nFl, 1
          call evolve_NS(nx, xx, T_NS(i,:), Q2i, Q2f, xi, oA=A, oQ2step=Q2step)
        end do
        !!!$OMP END PARALLEL DO

        call evolve_singlet(nx, xx, Sigma, Gluons, Q2i, Q2f, xi, oA=A, oQ2step=Q2step)

        ! Go back to physical basis
        call recouple_arrays(nx, H_mesh(-nFl:nFl,:), q_minus, T_NS, Sigma, Gluons, nFl)

        deallocate(q_minus, T_NS, Gluons, Sigma)
    end subroutine evolve_GPDs_fixed_nFl



end module evolvers
