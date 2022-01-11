module global_param
  use iso_fortran_env, only: REAL64
  implicit none

  logical,parameter :: debugON = .TRUE.

  integer,parameter :: rp = REAL64, n_sd = 2, buffer = 20

  real(kind=rp),parameter :: pi = acos(-1._rp)
end module global_param
