module global_param
  use iso_fortran_env, only: REAL64
  implicit none

  logical,parameter ::                                           &
    debugON = .TRUE.

  integer,parameter ::                                           &
    rp = REAL64,                                                 &
    N_buffer = 0

  real(kind=rp),parameter ::                                     &
    n_sd = 2,                                                    &
    pi   = acos(-1._rp),                                         &
    h_shear = 0.01_rp,                                           &
    R2_rnk  = 0.5_rp*h_shear**2,                                 &
    eps_adjst = 1.2_rp
end module global_param
