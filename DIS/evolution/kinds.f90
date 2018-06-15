! kinds.f90
!
! This just defines double precision. Needed by everything else.

module kinds
  implicit none
  integer, parameter, public :: dp = kind(1d0)
end module kinds
