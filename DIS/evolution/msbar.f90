module msbar

  use kinds, only : dp

  implicit none
  public

  ! Casimir invariants
  real(dp), public, parameter :: CF = 4.0_dp/3.0_dp
  real(dp), public, parameter :: CA = 3.0_dp
  real(dp), public, parameter :: TF = 1.0_dp/2.0_dp

  ! Mathematical constants
  real(dp), public, parameter :: pi = acos(-1.0_dp)
  real(dp), public, parameter :: zeta2 = pi**2 / 6.0_dp
  real(dp), public, parameter :: zeta3 = 1.2020569031_dp
  real(dp), public, parameter :: zeta4 = pi**4 / 90.0_dp

  ! Current quark masses for effective quark number thresholds
  real(dp), private, parameter :: cMass = 1.29_dp
  real(dp), private, parameter :: bMass = 4.5_dp
  real(dp), private, parameter :: tMass = 172.44_dp
  real(dp), public, parameter :: cMass2 = cMass**2
  real(dp), public, parameter :: bMass2 = bMass**2
  real(dp), public, parameter :: tMass2 = tMass**2

  real(dp), public, parameter, dimension(3:6) :: &
       & LambdaQCD = [ 0.332_dp, 0.292_dp, 0.210_dp, 0.089_dp ]



   contains



     function get_effective_flavors(Q2) result(Neff)
         real(dp), intent(in) :: Q2
         integer :: Neff
         if(Q2 .gt. tMass2) then
           Neff = 6
         else if(Q2 .gt. bMass2) then
           Neff = 5
         else if(Q2 .gt. cMass2) then
           Neff = 4
         else
           Neff = 3
         end if
     end function get_effective_flavors



     function get_effective_lambda(Q2) result(Lambda_eff)
         real(dp), intent(in) :: Q2
         real(dp) :: Lambda_eff
         integer :: Neff
         Neff = get_effective_flavors(Q2)
         Lambda_eff = LambdaQCD(Neff)
     end function get_effective_Lambda



     function get_tau(Q2) result(tau)
         real(dp), intent(in) :: Q2
         real(dp) :: tau
         real(dp) :: Lambda
         Lambda = get_effective_lambda(Q2)
         tau = log(Q2/Lambda**2)
     end function get_tau



     function get_alpha_QCD(Q2, onOrder) result(alphaQCD)
         ! Leading order for now.
         real(dp), intent(in) :: Q2
         real(dp) :: alphaQCD

         integer, intent(in), optional :: onOrder
         integer :: nOrder

         real(dp), dimension(0:4) :: a, b
         real(dp) :: t, nFl

         ! Allow nOrder to be given ... assume 1 if not
         nOrder = 1
         if(present(onOrder)) nOrder = onOrder

         ! Effective, Q2-dependent parameters
         t = get_tau(Q2)
         nFl = real(get_effective_flavors(Q2))

         b(0) = ( 11.*CA - 4.*nFl*TF ) / (12.*pi)
         b(1) = ( 17.*CA**2 - nFl*TF*(10.*CA+6.*CF) ) / (24.*pi**2)
         b(2) = ( 2857.-5033./9.*nFl+325/27.*nFl**2 ) / (128.*pi**3)
         b(3) = 0. ! TODO
         b(4) = 0. ! TODO

         ! Get the alpha contributions from various orders
         ! NB these are only approximations anyway (except at LO).
         a(0) = 1.

         a(1) = -( b(1)*log(t) ) / ( b(0)**2*t )

         a(2) = &
             ( b(1)**2*( log(t)**2 - log(t) - 1. ) + b(0)*b(2) ) &
             / ( b(0)**4*t**2 )

         a(3) = - ( &
             b(1)**3*( log(t)**3 - 5./2.*log(t)**2 - 2.*log(t) + 1./2. ) &
             + 3.*b(0)*b(1)*b(2)*log(t) &
             - 1./2.*b(0)**2*b(3) &
             ) / ( b(0)**6*t**3 )

         a(4) = 0.0_dp ! TODO

         ! Sum up the contributions
         alphaQCD = sum( a(0:nOrder-1) )  / ( b(0)*t )
     end function get_alpha_QCD



end module msbar
