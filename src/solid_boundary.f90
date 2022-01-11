module solid_boundary
  use global_param
  use utility
  use lapack95, only: gesv
  implicit none

  type :: boundary
    integer :: N = 0, Nmax = 0
    real(kind=rp),allocatable :: r_nod(:,:), r_col(:,:),          &
      r_shl(:,:), e_s(:,:), ds(:), b(:,:), A(:,:)

  contains
    procedure :: set
    procedure :: enforce_vel
    procedure :: Vel_src
    procedure :: shapeout
  end type boundary

contains
  subroutine set(this,grid,h_shear)
    class(boundary),intent(inout),target :: this
    real(kind=rp)  ,intent(in)           :: grid(:,:), h_shear

    this%N     = SIZE(grid,2)
    if(this%N > this%Nmax) then
      DEALLOCATE(this%r_nod)
      DEALLOCATE(this%r_col)
      DEALLOCATE(this%r_shl)
      DEALLOCATE(this%e_s)
      DEALLOCATE(this%ds)
      DEALLOCATE(this%b)
      DEALLOCATE(this%A)
    end if
    if(.NOT.ALLOCATED(this%r_nod)) then
      this%Nmax  = (this%N*(100+buffer))/100
      ALLOCATE(this%r_nod(n_sd,this%Nmax))
      ALLOCATE(this%r_col(n_sd,this%Nmax))
      ALLOCATE(this%r_shl(n_sd,this%Nmax))
      ALLOCATE(this%e_s(n_sd,this%Nmax))
      ALLOCATE(this%ds(this%Nmax))
      ALLOCATE(this%b(n_sd,this%Nmax))
      ALLOCATE(this%A(2*this%Nmax,2*this%Nmax))
    end if

!    this%gama   => this%b(1,1:this%N)
!    this%sigma  => this%b(2,1:this%N)
!    this%V_res  => this%b(:,1:this%N)

    associate(r    => this%r_nod(:,1:this%N),                    &
              r_ph => this%r_col(:,1:this%N),                    &
              r_p1 => CSHIFT(this%r_col(:,1:this%N),1),          &
              e    => this%e_s  (:,1:this%N),                    &
              rl   => this%r_shl(:,1:this%N),                    &
              ds_2 => this%ds(1:this%N))
      r    = grid
      r_ph = (r_p1+r)/2._rp
      e    =  r_p1-r
      ds   = NORM2(e,1)
      e    = e/ds
      rl   = r_ph - h_shear*rot90(e)
    end associate
  end subroutine set

  subroutine enforce_vel(bndry,ffield,U)
    class(boundary) ,intent(inout) :: bndry
    type(flow_field),intent(in)    :: ffield
    real(kind=rp)   ,intent(in)    :: U(n_sd)

    associate(V_res => bndry%b(:,1:bndry%N),                     &
              r_c   => bndry%r_col(1:bndry%N),                   &
              A     => bndry%A(1:this%N,1:this%N))
      V_res = -ffield%Vel(r_c,U)
      calc_inflcoeff: block
        integer :: i, j
        do concurrent(i = 1:bndry%N, j = 1:bndry%N)
          A(2*i-1:2*i,2*j-1:2*j) = coeff(i,j)
        end do
      end block calc_inflcoeff
      call gesv(A,V_res)
    end associate

  contains
    pure function coeff(m,n)
      integer      ,intent(in) :: m, n
      real(kind=rp)            :: dr(n_sd), xip, xim, eta, C1,   &
        C2, C3, coeff(2,2)

      dr  = bndry%r_col(:,m) - bndry%r_col(:,n)
      xip = proj(dr,n) + bndry%ds(n)/2._rp
      xim = proj(dr,n) - bndry%ds(n)/2._rp
      eta = bndry%rej (dr,n)
      C1  = 2._rp*pi*mag2(dr)
      C2  = 0.5_rp*log((xim**2+eta**2)/(xip**2+eta**2))
      C3  = atan2(xim,eta) - atan2(xip,eta)

      coeff(:,1) = [-dr(2)/C1, C2*bndry%e_s(1,n)-C3*bndry%e_s(2,n)]
      coeff(:,2) = [ dr(1)/C1, C2*bndry%e_s(2,n)+C3*bndry%e_s(1,n)]
    end function coeff
  end subroutine enforce_vel

  pure function proj(bndry,u,n)
    class(boundary),intent(in) :: bndry
    real(kind=rp)  ,intent(in) :: u(n_sd)
    integer        ,intent(in) :: n
    real(kind=rp) :: proj

    proj = DOT_PRODUCT(u,bndry%e_s(:,n))
  end function proj

  pure function rej(bndry,u,n)
    class(boundary),intent(in) :: bndry
    real(kind=rp)  ,intent(in) :: u(n_sd)
    integer        ,intent(in) :: n
    real(kind=rp) :: rej

    rej = DOT_PRODUCT(u,rot90(bndry%e_s(:,n)))
  end function rej

  pure function Vel_src(bndry,r)
    class(boundary),intent(in) :: bndry
    real(kind=rp)  ,intent(in) :: r(n_sd)
    real(kind=rp)              :: Vel_src(n_sd), dr(n_sd), xip,    &
                                  xim, eta, C2, C3
    integer                    :: i

    Vel_src = 0._rp
    do concurrent (i = 1:bndry%N)
      dr  = r - bndry%r_col(:,i)
      xip = bndry%proj(dr,i) + bndry%ds(i)/2._rp
      xim = bndry%proj(dr,i) - bndry%ds(i)/2._rp
      eta = bndry%rej (dr,i)
      C2  = 0.5_rp*log((xim**2+eta**2)/(xip**2+eta**2))
      C3  = atan2(xim,eta) - atan2(xip,eta)

      Vel_src = [C2*bndry%e_s(1,i)-C3*bndry%e_s(2,i),             &
                [C2*bndry%e_s(2,i)+C3*bndry%e_s(1,i)]*bndry%b(2,i)
    end do
  end function Vel_src

  subroutine shapeout(bndry,flname)
    class(bndry)    ,intent(in) :: bndry
    character(len=*),intent(in) :: flname
    integer :: i

    open(10,file=flname)
    do i = 1,bndry%N
      write(10,*) bndry%r(i)
    end do
    write(10,*) bndry%r(1)
    close(10)
  end subroutine shapeout
end module solid_boundary
