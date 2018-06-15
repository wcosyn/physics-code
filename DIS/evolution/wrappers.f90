! f2py can't figure out how to use kind=dp unless this is included.
include "kinds.f90"

module evowrappers
  use evolvers
  use kinds

  implicit none
  public

  contains



    subroutine evolve_gpds_(nx, x, nxi, xi, nt, t, Q2i, Q2f, A, H_mesh)
        integer, intent(in) :: nx, nxi, nt
        real(dp), dimension(nx), intent(in) :: x
        real(dp), dimension(nxi), intent(in) :: xi
        real(dp), dimension(nt), intent(in) :: t
        real(dp), intent(in) :: Q2i, Q2f
        real(dp), intent(in) :: A
        real(dp), dimension(-6:6,nx,nxi,nt), intent(inout) :: H_mesh
        integer :: i, j

        !$OMP PARALLEL DO
        do i=1, nxi, 1
        do j=1, nt, 1
          call evolve_gpds(nx, x, H_mesh(:,:,i,j), Q2i, Q2f, xi(i), oA=A)
        end do
        end do
        !$OMP END PARALLEL DO

    end subroutine evolve_gpds_



    subroutine evolve_ns_(nx, xx, HNS, Q2i, Q2f, xi, A)
        ! A is required input, since f2py doesn't seem to work with optionals.
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in) :: xx
        real(dp), dimension(nx), intent(inout) :: HNS
        real(dp), intent(in) :: Q2i, Q2f, xi
        real(dp), intent(in) :: A
        call evolve_ns(nx, xx, HNS, Q2i, Q2f, xi, oA=A)
    end subroutine evolve_ns_



    subroutine evolve_singlet_(nx, xx, Sigma, Gluons, Q2i, Q2f, xi, A)
        ! A is required input, since f2py doesn't seem to work with optionals.
        integer, intent(in) :: nx
        real(dp), dimension(nx), intent(in) :: xx
        real(dp), dimension(nx), intent(inout) :: Sigma, Gluons
        real(dp), intent(in) :: Q2i, Q2f, xi
        real(dp), intent(in) :: A
        call evolve_singlet(nx, xx, Sigma, Gluons, Q2i, Q2f, xi, oA=A)
    end subroutine evolve_singlet_



end module evowrappers
