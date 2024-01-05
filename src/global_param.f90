module global_param
  use iso_fortran_env, only: REAL64
  implicit none

  logical,parameter :: debugON = .TRUE., Tscheme = .TRUE.
  integer,parameter :: rp = REAL64, n_sd = 2, buffer = 0, epslim2 = 1,  &
                       nen = 3
  real(kind=rp),parameter :: pi = acos(-1._rp), tol = 1.E-6
end module global_param
