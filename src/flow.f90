module flow
  use global_param
  implicit none

  type :: flow_field
    integer :: N = 0, Nmax
    real(kind=rp)             :: R2_rnk, c = 1.2_rp
    real(kind=rp),allocatable :: Gama(:), r(:,:)
    logical      ,allocatable :: has_vortex(:)

  contains
    procedure :: init
    procedure :: Vel
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
    ALLOCATE(ffield%r(Nvor))
    ALLOCATE(ffield%has_vortex(Nvor))

    has_vortex = .FALSE.

    ffield%R2_rnk = R_rnk**2
    if(PRESENT(adjst_eps)) ffield%c = adjst_eps
  end subroutine init

  pure function Vel(ffield,bndry,r,U)
    class(ffield) ,intent(in) :: ffield
    type(boundary),intent(in) :: bndry
    real(kind=rp ,intent(in)  :: r(:,:), U(n_sd)
    real(kind=rp)             :: Vel(n_sd,SIZE(r,2))
    integer                   :: i

    do concurrent (i = 1:SIZE(r))
      Vel(:,i) = U + bndry%Vel_src(r(:,i))
      do concurrent (j = 1:ffield%Nmax, ffield%has_vortex(j))
        Vel(:,i) = Vel(:,i) + Vel_ind(j,r(:,i))
      end do
    end do

  contains
    pure function Vel_ind(i_vor,r_p)
      integer      ,intent(in) :: i_vor
      real(kind=rp),intent(in) :: r_p(n_sd)
      real(kind=rp)            :: Vel_ind(n_sd)

      Vel_ind = rot90(r_p - ffield%r(:,i_vor))
      r2      = mag2(Vel_ind)

      if(r2 < R2_rnk) then
        Vel_ind = ffield%Gama(i_vor)/(2._rp*pi*ffield%R2_rnk)    &
                * Vel_ind
      else
        Vel_ind = (ffield%Gama(i_vor)/2._rp*pi*r2)*Vel_ind
      end if
    end function Vel_ind
  end function Vel

  pure subroutine step(ffield,bndry,U,Re,dt)
    class(flow_field),intent(inout) :: ffield
    type (boundary)  ,intent(in)    :: bndry
    real(kind=rp)    ,intent(in)    :: U(n_sd), Re, dt

    do concurrent (i_vor = 1:ffield%N)
      call move(i_vor)
    end do
    if(ffield%N + bndry%N > ffield%Nmax)                         &
      call ffield%purge_vor(ffield%N + bndry%N - ffield%Nmax)
    call bndry%enforce_vel(ffield,U)
    call ffield%import_vor(bndry)
  contains
    pure subroutine move(i)
      integer,intent(in) :: i

      do concurrent (j = 1:ffield%Nmax,                          &
                     j /= i .AND. ffield%has_vortex(j))
        dr(j) = ffield%r(:,i) - ffield%r(:,j)
        r (j) = NORM2(dr(:,j))
      end do
      eps = ffield%c * MIN(r)

      V_pot = bndry%Vel_src(ffield%r(i)) + U
      I1    = 0._rp
      I2    = 0._rp
      do concurrent (j = 1:ffield%Nmax,                          &
                     j /= i .AND. ffield%has_vortex(j))
        V_pot = V_pot + rot90(dr(j))                             &
                      * ffield%Gama(j)/MAX(r(j)**2,ffield%R2_rnk)
        tmp   = ffield%Gama(j)*exp(-r(j)/eps)
        I1    = I1 + tmp
        I2    = I2 + dr(j)*tmp/(eps*r(j))
      end do

      do concurrent (j = 1:bndry%N)
        dr(j) = ffield%r(:,i) - bndry%r_col(:,j)
        r (j) = NORM2(dr(:,j))
      end do
      eps = ffield%c * MIN(r(1:bndry%N))

      I0 = 2._rp*pi*eps**2
      I3 = 0._rp
      do concurrent (j = 1:bndry%N)
        tmp = exp(-r(j)/eps)
        I0  = I0 + eps*rej(bndry%r_col(:,j),j)                   &
                 * (r(j)+eps)*tmp/r(j)**2
        I3  = I3 + bndry%ds(j)*tmp*rot90(bndry%e_s(:,j))
      end do
      ffield%r(:,i) = ffield%r(:,i)
                    + (dt/2._rp/pi)*V_pot
                    + (dt/Re/I0)*I3
                    - (dt/Re/I1)*I2
    end subroutine move
  end subroutine step

  pure subroutine purge_vor(ffield,num)
    class(flow_field),intent(inout) :: ffield
    integer          ,intent(in)    :: num

    ffield%has_vortex(largestloc(mag2(ffield%r),num)) = .FALSE.
    ffield%N = ffield%N - num
  end subroutine purge_vor

  pure subroutine import_vor(ffield,bndry)
    class(flow_field),intent(inout) :: ffield
    type(boundry)    ,intent(in)    :: bndry
    integer                         :: i, j

    j = 0
    do i = 1,ffield%Nmax
      if(.NOT.ffield%has_vortex(i)) then
        j = j + 1
        ffield%has_vortex(i) = .TRUE.
        ffield%Gama(i)       = bndry%b(1,j)
        ffield%r(i)          = bndry%r(:,j)
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
