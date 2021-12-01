module vortsol
  use dashboard
  use vortex
  use contour
  implicit none

  ! domain and shear layer vortices
  type(vortex_t),private,target  :: vortex(n_vor)                 &
    = vortex_t(0._rp,vector_t([huge_real,huge_real]))

contains
  subroutine evolve_flow(Uinfty,contr,dt)
    type(vector_t), intent(in) :: Uinfty
    type(contour_t),intent(in) :: contr
    real(kind=rp)  ,intent(in) :: dt
    integer                    :: n_contr

    call move_vortex(vortex,Uinfty,dt)
    v_slip = Vel(contr) ! should we move this to L1?
    call drop_vortices(vortex,n_contr)
!L1
    associate (shl_vortex => vortex(1:n_contr),                   &
               dom_vortex => vortex(n_contr+1:))
      call locate_shlayer(shl_vortex,contr)
      call calc_shl_gama(shl_vortex,)
    end associate

  end subroutine evolve_flow

  pure subroutine move_vortex(vor,U,dt)
    type(vortex_t),intent(inout) :: vor(:)
    type(vector_t),intent(in)    :: U
    real(kind=rp) ,intent(in)    :: dt

    do concurrent (i = 1:n_vor, j = 1:n_vor, i /= j)
      r_ij = vor
    end do
    vor%r = vor%r + dt*(U + Vel_t(vor%r))
  end subroutine move_vortex


  pure elemental function Vel(x)
    type(vector_t),intent(in) :: x
    type(vector_t)            :: Vel
    real(kind=rp)             :: r2
    integer                   :: i

    if(vel_exact) then
      do concurrent (i = 1:n_vor)
        r2 = vortex(i)%dist2_from(x)
        Vel = Vel + vortex(i)%Vel_ind(x) 
      end do
    end if
  end function Vel

  pure elemental function Vel_d(r,rc)
    type(vector_t),intent(in) :: r, rc(:)
    type(vector_t)            :: Vel_d, I2, I3, r_ij
    real(kind=rp)             :: I0, I1, tmp
    integer                   :: j, k

    eps = c*sqrt(MIN(vortex%dist2_from(r_i)))

    I0 = 0._rp; I1 = 0._rp
    I2 = zero ; I3 = zero
    do concurrent (j = 1:n_vor)
      associate (Gama_j => vortex(j)%Gama, r_j => vortex(j)%r)
        r_ij = mag(r_i - r_j)
        tmp  = Gama_j*exp(-r_ij/eps)
        I1   = I1 + tmp
        I2   = I2 + (tmp/(eps*r_ij))*(r_i - r_j)
      end associate
    end do
  end function Vel_d

  pure subroutine drop_vortices(vor,n)
    type(vortex_t),intent(inout) :: vor
    integer       ,integer(in)   :: n
    integer                      :: i, f(n)

    call find_farthest_vortices(f,n)
    do concurrent(i = 1:n)
      vor(f(i)) = vor(i)
    end do
  end subroutine drop_vortices

  pure subroutine find_farthest_vortices(f,n)
    integer,intent(in)    :: n
    integer,intent(inout) :: f(n)
    integer               :: i
    real(kind=rp)         :: r2, r2_f(n)

    f    = [(i, i = 1,n)]
    r2_f = vortex(1:n)%dist2_from()

    do i = n+1,n_vor
      r2 = vortex(i)%dist2_from()
      if(r2 > MAX(r2_f)) then
        k       = MINLOC(r2_f)
        f(k)    = i
        r2_f(k) = r2
      end if
    end do
  end subroutine find_farthest_vortices

  ! calc. the location of shear layer vortices
  subroutine locate_shlayer(layer,contr)
    type(vortex_t),intent(inout) :: layer(:)
    type(vector_t),intent(in)    :: contr(:)

    call locate(layer,r_contr,CSHIFT(r_contr,1))
  contains
    pure elemental function layer(shl,vec,nxt_vec)
      type(vortex_t) ,intent(inout) :: shl
      type(vector_t),intent(in)    :: vec, nxt_vec
      type(vector_t)               :: normal

      normal = rot90(nxt_vec - vec)
      shl    = 0.5_rp*(nxt_vec + vec) + (h/mag(normal)**2)*normal
    end function layer
  end subroutine locate_shlayer

  ! calc. the strength of shear layer votices
  subroutine gen_vortices(r_contr)
    type(vector_t),intent(in) :: r_contr(:)
    type(vector_t) :: v_slip(SIZE(r_contr))

    v_slip = velocity(r_contr)
    call fix_slip(omega_layr,v_contr,r_contr)
    call insert_vortices(omega_layr,r_layr)
  contains
  end subroutine gen_vortices

  pure elemental function velocity(r)
    type(vector_t),intent(in) :: r
    type(vector_t)            :: velocity

    if(vel_exact) then
        velocity = SUM_vec(exact_vel(r-r_vor,omega_vor))
    end if
  contains
    ! exact velocity induced by a single vortex
    pure elemental function exact_vel(r0,w)
      type(vector_t)    ,intent(in) :: r0
      real(kind=rp),intent(in) :: w
      type(vector_t)     :: exact_vel
      real(kind=rp) :: dist

      dist      = mag(r0)
      exact_vel = 0.5_rp*w*rot90(r0)

      if(dist > r_rnk) exact_vel = exact_vel*(r_rnk/dist)**2
    end function exact_vel 
  end function velocity
end module vortsol
