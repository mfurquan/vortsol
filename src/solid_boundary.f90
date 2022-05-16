module solid_boundary
  use global_param
  use utility
  use flow_elements
  implicit none

  type :: boundary
    integer :: N = 0, Nmax = 0
    real(kind=rp),allocatable :: r_nod(:,:), r_col(:,:),         &
      r_shl(:,:), e_s(:,:), ds(:), b(:), A(:,:)

  contains
    procedure :: set
    procedure :: Vel_src
    procedure :: proj
    procedure :: rej
    procedure :: shapeout
  end type boundary

contains
  subroutine set(this,grid,h_shear)
    class(boundary),intent(inout) :: this
    real(kind=rp)  ,intent(in)    :: grid(:,:), h_shear

    this%N     = SIZE(grid,2)
    if(this%Nmax > 0 .AND. this%N > this%Nmax) then
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
      ALLOCATE(this%b(n_sd*this%Nmax))
      ALLOCATE(this%A(2*this%Nmax,2*this%Nmax))
    end if

    this%b(1:2*this%N)     = 0._rp
    this%r_nod(:,1:this%N) = grid
    associate(r    => this%r_nod(:,1:this%N),                    &
              r_ph => this%r_col(:,1:this%N),                    &
              r_p1 => CSHIFT(this%r_nod(:,1:this%N),1,2),        &
              e    => this%e_s  (:,1:this%N),                    &
              rl   => this%r_shl(:,1:this%N),                    &
              ds   => this%ds(1:this%N))
      r_ph = (r_p1+r)/2._rp
      e    =  r_p1-r
      ds   = NORM2(e,1)
      e    = e/ds
      rl   = r_ph - h_shear*rot90(e)
    end associate
  end subroutine set

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
    real(kind=rp)              :: Vel_src(n_sd)
    integer                    :: i

    associate(q   => bndry%b(2:2*bndry%N:2),                     &
              rm  => bndry%r_col(:,1:bndry%N),                   &
              ds  => bndry%ds(1:bndry%N),                        &
              e_s => bndry%e_s(:,1:bndry%N))
      Vel_src = 0._rp
      do concurrent (i = 1:bndry%N)
        Vel_src = Vel_src                                        &
                + q(i)*src_sheet(r,rm(:,i),ds(i),e_s(:,i))
      end do
    end associate
  end function Vel_src

  subroutine shapeout(bndry,flname)
    class(boundary) ,intent(in) :: bndry
    character(len=*),intent(in) :: flname
    integer :: i

    open(10,file=flname)
    do i = 1,bndry%N
      write(10,*) bndry%r_nod(:,i)
    end do
    write(10,*) bndry%r_nod(:,1)
    close(10)
  end subroutine shapeout
end module solid_boundary
