module decouple
  ! A module to decouple the PDFs/GPDs for evolution,
  ! or to recombine them into the "physical" PDFs/GPDs after evolution.

  use kinds, only : dp

  implicit none
  private

  public :: decouple_arrays, recouple_arrays

  contains



    subroutine decouple_arrays(Nx, pdf_mesh, q_minus, TNS, Sigma, Gluons, nFl)
        ! Required input:
        !   - Nx : number of x evaluations
        !   - pdf_mesh : array with dimension(-6:6,Nx) with PDFs (or GPDs)
        ! The output:
        !   - q_minus is dimension(1:6,Nx)
        !   - TNS is dimension(2:6,Nx)
        !   - Sigma, Gluons are dimension(Nx)

        ! I/O
        integer, intent(in) :: Nx, nFl
        real(dp), dimension(-nFl:nFl,Nx), intent(in) :: pdf_mesh
        real(dp), dimension(1:nFl,Nx), intent(out) :: q_minus
        real(dp), dimension(2:nFl,Nx), intent(out) :: TNS
        real(dp), dimension(Nx), intent(out) :: Sigma, Gluons

        ! Intermediates
        real(dp), dimension(:,:), allocatable :: q_plus
        integer :: i, k

        allocate(q_plus(1:nFl,Nx))

        Gluons(:) = pdf_mesh(0,:)

        forall(i=1:nFl)
            q_plus(i,:)  = pdf_mesh(i,:) + pdf_mesh(-i,:)
            q_minus(i,:) = pdf_mesh(i,:) - pdf_mesh(-i,:)
        end forall

        forall(i=1:Nx)
            Sigma(i) = sum( q_plus(:,i) )
        end forall

        forall(i=1:Nx)
            forall(k=2:nFl)
                TNS(k,i) = sum( q_plus(1:k,i) ) - k*q_plus(k,i)
            end forall
        end forall

        deallocate(q_plus)
    end subroutine decouple_arrays



    subroutine recouple_arrays(Nx, pdf_mesh, q_minus, TNS, Sigma, Gluons, nFl)
        ! Basically, decouple_arrays in reverse.
        ! Input:
        !   - Nx, q_minus, TNS, Sigma, Gluons
        ! Output:
        !   - pdf_mesh

        ! I/O
        integer, intent(in) :: Nx, nFl
        real(dp), dimension(1:nFl,Nx), intent(in) :: q_minus
        real(dp), dimension(2:nFl,Nx), intent(in) :: TNS
        real(dp), dimension(Nx), intent(in) :: Sigma, Gluons
        real(dp), dimension(-nFl:nFl,Nx), intent(out) :: pdf_mesh

        ! Intermediates
        real(dp), dimension(:,:), allocatable :: q_plus
        real(dp), dimension(:), allocatable :: dummy_sigma
        integer :: i, k

        allocate(q_plus(1:nFl,Nx), dummy_sigma(Nx))

        pdf_mesh(0,:) = Gluons
        dummy_sigma = Sigma

        do k=nFl, 2, -1
          q_plus(k,:) = ( dummy_sigma(:) - TNS(k,:) ) / real(k)
          dummy_sigma(:) = ( real(k-1)*dummy_sigma(:) + TNS(k,:) ) / real(k)
        end do
        q_plus(1,:) = dummy_sigma(:)

        !q_plus = 0. ! TEST

        forall(i=1:nFl)
            pdf_mesh( i,:) = 0.5_dp*( q_plus(i,:) + q_minus(i,:) )
            pdf_mesh(-i,:) = 0.5_dp*( q_plus(i,:) - q_minus(i,:) )
        end forall

        deallocate(q_plus, dummy_sigma)
    end subroutine recouple_arrays

  !subroutine debrief_arrays()
  !  ! Inversion of brief_arrays; using q_minus, TNS, Sigma, and G,
  !  ! we get back the pdf_mesh. To be used after evolution.
  !  implicit none

  !  real(dp), dimension(:,:), allocatable :: q_plus
  !  integer i, k

  !  ! Gluon is easy part
  !  pdf_mesh(0,:) = Gluons

  !  ! Reduction sequence for unpacking TNS
  !  allocate( q_plus(6, Nx) )
  !  do k=6, 2, -1
  !    q_plus(k,:) = ( Sigma(:) - TNS(k,:) ) / k
  !    Sigma(:) = ( (k-1)*Sigma(:) + TNS(k,:) ) / k
  !  end do
  !  q_plus(1,:) = Sigma(:)

  !  ! With q_plus and q_minus, the rest is easy
  !  forall(i=1:6)
  !      pdf_mesh( i,:) = 0.5_dp*( q_plus(i,:) + q_minus(i,:) )
  !      pdf_mesh(-i,:) = 0.5_dp*( q_plus(i,:) - q_minus(i,:) )
  !  end forall

  !  ! This step is partnered with brief_arrays, so clean up its garbage too.
  !  deallocate( Gluons, TNS, q_plus, q_minus, Sigma )

  !  return

  !end subroutine debrief_arrays




    !! Creates linear combinations of "physical" PDFs that (mostly) decouple
    !! in their evolution.
    !implicit none

    !! q_plus is a temporary array since it's not used in the evolution
    !real(dp), dimension(:,:), allocatable :: q_plus
    !integer i, k

    !! Gluons
    !allocate( Gluons(Nx) )
    !Gluons(:) = pdf_mesh(0,:)

    !! q+ and q- are easiest
    !allocate( q_plus(6, Nx) )
    !allocate( q_minus(6, Nx) )
    !forall(i=1:6)
    !    q_plus(i,:)  = pdf_mesh(i,:) + pdf_mesh(-i,:)
    !    q_minus(i,:) = pdf_mesh(i,:) - pdf_mesh(-i,:)
    !end forall

    !! Sigma for singlet evolution
    !allocate( Sigma(Nx) )
    !forall(i=1:Nx)
    !    Sigma(i) = sum( q_plus(:,i) )
    !end forall

    !! T, like q_minus, evovles by non-singlet evolution
    !allocate( TNS(2:6, Nx) )
    !forall(i=1:Nx)
    !    forall(k=2:6)
    !        TNS(k,i) = sum( q_plus(1:k,i) ) - k*q_plus(k,i)
    !    end forall
    !end forall

    !! Garbage collection and return
    !deallocate( q_plus )

    !return


end module decouple
