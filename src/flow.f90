module flow
  use global_param
  use utility
  use solid_boundary
  use lapack95, only: gesv
  implicit none

  type :: flow_field
    integer :: N = 0, Nmax
    real(kind=rp)             :: R2_rnk, c = 1.2_rp
    real(kind=rp),allocatable :: Gama(:), r(:,:)
    logical      ,allocatable :: has_vortex(:)

  contains
    procedure :: init
    procedure :: Vel
    procedure :: enforce_vel
    procedure :: purge_vor
    procedure :: import_vor
    procedure :: step
    procedure :: write_vor
  end type flow_field

contains
  pure subroutine init(ffield,Nvor,R_rnk,adjst_eps)
    class(flow_field),intent(inout) :: ffield
    integer          ,intent(in)    :: Nvor
    real(kind=rp)    ,intent(in)    :: R_rnk
    real(kind=rp)    ,intent(in),optional :: adjst_eps

    ffield%Nmax = Nvor
    ALLOCATE(ffield%Gama(Nvor))
    ALLOCATE(ffield%r(n_sd,Nvor))
    ALLOCATE(ffield%has_vortex(Nvor))

    ffield%has_vortex = .FALSE.

    ffield%R2_rnk = R_rnk**2
    if(PRESENT(adjst_eps)) ffield%c = adjst_eps
  end subroutine init

  pure function Vel(ffield,r,U,bndry)
    class(flow_field) ,intent(in)          :: ffield
    real(kind=rp)     ,intent(in)          :: r(:,:), U(n_sd)
    type(boundary)    ,intent(in),optional :: bndry
    real(kind=rp)                          :: Vel(n_sd,SIZE(r,2))
    integer                                :: i, j

    associate(r_v    => ffield%r,                                 &
              R2_rkn => ffield%R2_rnk,                            &
              gama   => ffield%Gama)

      do concurrent (i = 1:SIZE(r,2))
        if(PRESENT(bndry)) then
          Vel(:,i) = U + bndry%Vel_src(r(:,i))
        else
          Vel(:,i) = U
        end if

        do j = 1,ffield%Nmax
          if(ffield%has_vortex(j)) Vel(:,i) = Vel(:,i)            &
                 + gama(j)*rkn_vortex(r(:,i),r_v(:,j),R2_rkn)
        end do
      end do
    end associate
  end function Vel

  pure subroutine enforce_vel(ffield,bndry,U)
    class(flow_field),intent(in)    :: ffield
    type(boundary)   ,intent(inout) :: bndry
    real(kind=rp)    ,intent(in)    :: U(n_sd)

    associate(V_res => bndry%b(1:n_sd*bndry%N),                   &
              r_c   => bndry%r_col(:,1:bndry%N),                  &
              r_sh  => bndry%r_shl(:,1:bndry%N),                  &
              ds    => bndry%ds(1:bndry%N),                       &
              e_s   => bndry%e_s(:,1:bndry%N),                    &
              A     => bndry%A(1:n_sd*bndry%N,1:n_sd*bndry%N),    &
              R2    => ffield%R2_rnk)
      V_res = -ffield%Vel(r_c,U)
      calc_inflcoeff: block
        integer :: i, j
        do concurrent(i = 1:bndry%N, j = 1:bndry%N)
          A(2*i-1:2*i,2*j-1) = rkn_vortex(r_c(:,i),r_sh(:,j),R2)
          A(2*i-1:2*i,2*j  ) = src_sheet (r_c(:,i),r_c(:,j),ds(j),&
                                          e_s(:,j))
        end do
      end block calc_inflcoeff
      call gesv(A,V_res)
    end associate
  end subroutine enforce_vel

  pure subroutine step(ffield,bndry,U,Re,dt)
    class(flow_field),intent(inout) :: ffield
    type (boundary)  ,intent(in)    :: bndry
    real(kind=rp)    ,intent(in)    :: U(n_sd), Re, dt
    integer                         :: i_vor

    do concurrent (i_vor = 1:ffield%N)
      call move(i_vor)
    end do
    if(ffield%N + bndry%N > ffield%Nmax)                         &
      call ffield%purge_vor(ffield%N + bndry%N - ffield%Nmax)
    call ffield%enforce_vel(bndry,U)
    call ffield%import_vor(bndry)
  contains
    pure subroutine move(i)
      integer,intent(in) :: i
      real(kind=rp)      :: r, dr(n_sd), V_pot(n_sd), eps, tmp,   &
                            I1, I2(n_sd), I0, I3(n_sd)
      integer            :: j

      eps = 1.e5
      do concurrent (j = 1:ffield%Nmax,                           &
                     j /= i .AND. ffield%has_vortex(j))
        eps = MIN(mag2(ffield%r(:,i) - ffield%r(:,j)),eps)
      end do
      eps = ffield%c * sqrt(eps)

      V_pot = bndry%Vel_src(ffield%r(:,i)) + U
      I1    = 0._rp
      I2    = 0._rp
      do concurrent (j = 1:ffield%Nmax,                           &
                     j /= i .AND. ffield%has_vortex(j))
        dr    = ffield%r(:,i) - ffield%r(:,j)
        r     = NORM2(dr)
        V_pot = V_pot + rot90(dr)                                 &
                      * ffield%Gama(j)/MAX(r**2,ffield%R2_rnk)
        tmp   = ffield%Gama(j)*exp(-r/eps)
        I1    = I1 + tmp
        I2    = I2 + dr*tmp/(eps*r)
      end do

!      do concurrent (j = 1:bndry%N)
!        dr(:,j) = ffield%r(:,i) - bndry%r_col(:,j)
!        r (j) = NORM2(dr(:,j))
!      end do
!      eps = ffield%c * MINVAL(r(1:bndry%N))

      I0 = 2._rp*pi*eps**2
      I3 = 0._rp
      do concurrent (j = 1:bndry%N)
        r = NORM2(ffield%r(:,i) - bndry%r_col(:,j))
        tmp = exp(-r/eps)
        I0  = I0 + eps*bndry%rej(bndry%r_col(:,j),j)              &
                 * (r+eps)*tmp/r**2
        I3  = I3 + bndry%ds(j)*tmp*rot90(bndry%e_s(:,j))
      end do
      ffield%r(:,i) = ffield%r(:,i)                               &
                    + (dt/2._rp/pi)*V_pot                         &
                    + (dt/Re/I0)*I3                               &
                    - (dt/Re/I1)*I2
    end subroutine move
  end subroutine step

  pure subroutine purge_vor(ffield,num)
    class(flow_field),intent(inout) :: ffield
    integer          ,intent(in)    :: num

    ffield%has_vortex(                                            &
      largestloc(mag2(ffield%r,ffield%has_vortex),                &
                 num))                                            &
    = .FALSE.
    ffield%N = ffield%N - num
  end subroutine purge_vor

  pure subroutine import_vor(ffield,bndry)
    class(flow_field),intent(inout) :: ffield
    type(boundary)   ,intent(in)    :: bndry
    integer                         :: i, j

    j = 0
    do i = 1,ffield%Nmax
      if(.NOT.ffield%has_vortex(i)) then
        j = j + 1
        ffield%has_vortex(i) = .TRUE.
        ffield%Gama(i)       = bndry%b(2*j-1)
        ffield%r(:,i)        = bndry%r_shl(:,j)
      end if
      if(j == bndry%N) exit
    end do
    ffield%N = ffield%N + bndry%N
  end subroutine import_vor

  subroutine write_vor(ffield,flname,flnum)
    class(flow_field),intent(in) :: ffield
    character(len=*) ,intent(in) :: flname
    integer          ,intent(in),optional :: flnum
    character(len=7) :: stamp
    integer          :: i

    if(PRESENT(flnum)) then
      write(stamp,'(i3.3a4)') flnum,'.dat'
      open(unit=10,file=flname//stamp)
    else
      open(unit=10,file=flname//'.dat')
    end if
              
    do i = 1,ffield%Nmax
      if(ffield%has_vortex(i))                                   &
        write(10,*) ffield%r(:,i),ffield%Gama(i)
    end do
    close(10)
  end subroutine write_vor
end module flow
