module vortsol
  use dashboard
  implicit none

  ! STOCK ARRAY:
  real(kind=rp),private,save :: 
    r_vor (n_sd,n_vor),          & ! location of vortices
    r_layr(n_sd,nmax_contr),     & ! location of shear layer
    omega_vor(n_vor),            & ! strength of vortices
    omega_layr(nmax_contr)         ! strength of shear layer

  ! ALIASES:
  real(kind=rp),private,save,pointer ::
    rp_layr(:,:),                & ! location of shear layer
    rp_contr(:,:)                  ! location of solid boundary

contains
  subroutine update_layer(r_contr,n_contr)
    real(kind=rp),intent(in) :: r_contr(:,:)
    integer      ,intent(in) :: n_contr

    rp_layr  => r_layr (:,1:n_contr)
    rp_contr => r_contr(:,1:n_contr)

!    rp_layr = layer(rp_contr,CSHIFT(rp_contr,1))
!  contains
!    pure elemental function layer(vec,nxt_vec)
!      type(vec),intent(in) :: vec, nxt_vec
!      type(vec) :: layer, normal
!
!      normal = rot90(nxt_vec - vec)
!      layer  = 0.5_rp*(nxt_vec + vec) + (h/mag(normal))*normal
!    end function layer
  end subroutine update_layr

!  subroutine gen_vortices
!    real(kind=rp) :: v_slip(n_contr), omega_layr(n_contr)
!
!    v_contr = velocity(r_contr)
!    call fix_slip(omega_layr,v_contr,r_contr)
!    call insert_vortices(omega_layr,r_layr)
!  contains
!  end subroutine gen_vortices
!
!  pure elemental function velocity(r)
!    type(vec),intent(in) :: r
!    type(vec)            :: velocity
!
!    if(vel_exact) then
!        velocity = SUM_vec(exact_vel(r-r_vor,omega_vor))
!    end if
!  contains
!    ! exact velocity induced by a single vortex
!    pure elemental function exact_vel(r0,w)
!      type(vec)    ,intent(in) :: r0
!      real(kind=rp),intent(in) :: w
!      type(vec)     :: exact_vel
!      real(kind=rp) :: dist
!
!      dist      = mag(r0)
!      exact_vel = 0.5_rp*w*rot90(r0)
!
!      if(dist > r_rnk) exact_vel = exact_vel*(r_rnk/dist)**2
!    end function exact_vel 
!  end function velocity
end module vortsol
