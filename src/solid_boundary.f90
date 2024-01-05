module solid_boundary
  use global_param
  use utility
  implicit none

  private :: decom

  type :: boundary
    integer :: Ncn = 0, Nmax = 0
    logical :: flexible, moving, var_coeff
    real(kind=rp) :: delta       ! distance to shear layer
    real(kind=rp),allocatable ::                                          &
      r_nod(:,:),               &! nodal coords. for the boundary
      r_col(:,:),               &! collocation pts. (centre of panels)
      tau(:,:),                 &! tangent at r_col
      ds(:),                    &! panel length
      gama(:),                  &! free vortex sheet
      gama_att(:),              &! attached vortex sheet
      q_att(:),                 &! attached source sheet
      A(:,:)                     ! Matrix of coeffs.

    integer,allocatable ::      &
      ipiv(:),                  &! order of partial pivoting
      marker(:)                  ! location of marker particles

  contains
    procedure :: set
    procedure :: n_1
    procedure :: n_n
    procedure :: tau_dot_1
    procedure :: n_dot_1
    procedure :: tau_dot_n
    procedure :: n_dot_n
    generic   :: tau_dot => tau_dot_1, tau_dot_n
    generic   :: n_dot => n_dot_1, n_dot_n
    generic   :: n => n_1, n_n
    procedure :: inside
    procedure :: shapeout
    procedure :: decom
    !final     :: decom
  end type boundary

contains
  subroutine set(this,grid,h,is_flexible,is_moving)
    class(boundary),intent(inout) :: this
    real(kind=rp)  ,intent(in)    :: grid(:,:), h
    logical        ,intent(in)    :: is_flexible
    logical        ,intent(in),optional :: is_moving

    this%delta  = h
    this%Ncn      = SIZE(grid,2)
    if(is_flexible) then
      this%flexible = .TRUE.
      this%moving   = .TRUE.
    else
      this%flexible = .FALSE.
      this%moving   = is_moving
    end if
    this%var_coeff = .TRUE.

    if(this%Nmax > 0 .AND. this%Ncn > this%Nmax) then
      DEALLOCATE(this%r_nod)
      DEALLOCATE(this%r_col)
      DEALLOCATE(this%tau)
      DEALLOCATE(this%ds)
      DEALLOCATE(this%gama)
      DEALLOCATE(this%gama_att)
      DEALLOCATE(this%q_att)
      DEALLOCATE(this%A)
      DEALLOCATE(this%ipiv)
      DEALLOCATE(this%marker)
    end if
    if(.NOT.ALLOCATED(this%r_nod)) then
      this%Nmax  = (this%Ncn*(100+buffer))/100
      ALLOCATE(this%r_nod(n_sd,this%Nmax))
      ALLOCATE(this%r_col(n_sd,this%Nmax))
      ALLOCATE(this%tau(n_sd,this%Nmax))
      ALLOCATE(this%ds(this%Nmax))
      ALLOCATE(this%gama(this%Nmax+1))
      ALLOCATE(this%gama_att(this%Nmax))
      ALLOCATE(this%q_att(this%Nmax))
      ALLOCATE(this%A(this%Nmax+1,this%Nmax+1))
      ALLOCATE(this%ipiv(this%Nmax+1))
      ALLOCATE(this%marker(this%Nmax+1))
    end if

    this%gama     = 0._rp
    this%gama_att = 0._rp
    this%q_att    = 0._rp
    this%r_nod(:,1:this%Ncn)  = grid
    associate(r    => this%r_nod(:,1:this%Ncn),                    &
              r_ph => this%r_col(:,1:this%Ncn),                    &
              r_p1 => CSHIFT(this%r_nod(:,1:this%Ncn),1,2),        &
              ds   => this%ds(1:this%Ncn),                         &
              tau  => this%tau(:,1:this%Ncn))
      r_ph = (r + r_p1)/2._rp
      tau  =  r_p1 - r
      ds   =  NORM2(tau,1)
      tau  =  tau/ds !NORM2(tau,1)
    end associate
  end subroutine set

  pure subroutine decom(this)
    class(boundary),intent(inout) :: this
      DEALLOCATE(this%r_nod)
      DEALLOCATE(this%r_col)
      DEALLOCATE(this%tau)
      DEALLOCATE(this%ds)
      DEALLOCATE(this%gama)
      DEALLOCATE(this%gama_att)
      DEALLOCATE(this%q_att)
      DEALLOCATE(this%A)
      DEALLOCATE(this%ipiv)
      DEALLOCATE(this%marker)
  end subroutine decom

  pure function n_1(bndry,k)
    class(boundary),intent(in) :: bndry
    integer        ,intent(in) :: k
    real(kind=rp) :: n_1(n_sd)

    n_1(1) = bndry%tau(2,k); n_1(2) = -bndry%tau(1,k)
  end function n_1

  pure function n_n(bndry)
    class(boundary),intent(in) :: bndry
    real(kind=rp) :: n_n(n_sd,bndry%Nmax)

    n_n(1,:) = bndry%tau(2,:); n_n(2,:) = -bndry%tau(1,:)
  end function n_n

  pure function tau_dot_n(bndry,u)
    class(boundary),intent(in) :: bndry
    real(kind=rp)  ,intent(in) :: u(:,:)
    real(kind=rp)              :: tau_dot_n(SIZE(u,2))
    integer                    :: i

    do concurrent (i = 1:SIZE(u,2))
      tau_dot_n(i) = bndry%tau(1,i)*u(1,i) + bndry%tau(2,i)*u(2,i)
    end do
  end function tau_dot_n

  pure function n_dot_n(bndry,u)
    class(boundary),intent(in) :: bndry
    real(kind=rp)  ,intent(in) :: u(:,:)
    real(kind=rp)              :: n_dot_n(SIZE(u,2))
    integer                    :: i

    do concurrent (i = 1:SIZE(u,2))
      n_dot_n(i) = bndry%tau(2,i)*u(1,i) - bndry%tau(1,i)*u(2,i)
    end do
  end function n_dot_n

  pure function tau_dot_1(bndry,u,k)
    class(boundary),intent(in) :: bndry
    real(kind=rp)  ,intent(in) :: u(n_sd)
    integer        ,intent(in) :: k
    real(kind=rp)              :: tau_dot_1

    tau_dot_1 = bndry%tau(1,k)*u(1) + bndry%tau(2,k)*u(2)
  end function tau_dot_1

  pure function n_dot_1(bndry,u,k)
    class(boundary),intent(in) :: bndry
    real(kind=rp)  ,intent(in) :: u(n_sd)
    integer        ,intent(in) :: k
    real(kind=rp)              :: n_dot_1

    n_dot_1 = bndry%tau(2,k)*u(1) - bndry%tau(1,k)*u(2)
  end function n_dot_1

  pure function inside(bndry,r)
  ! use ray casting to check if a point is inside the polygon
    class(boundary),intent(in) :: bndry
    real(kind=rp)  ,intent(in) :: r(n_sd)
    real(kind=rp) :: m
    integer       :: i, ip
    logical :: inside

    inside = .FALSE.
    do i = 1,bndry%Ncn
      ip = i + 1
      if(i == bndry%Ncn) ip = 1
      associate (x1 => bndry%r_nod(1,i) , y1 => bndry%r_nod(2,i),        &
                 x2 => bndry%r_nod(1,ip), y2 => bndry%r_nod(2,ip),       &
                 x  => r(1)             , y  => r(2))
        if(y > MIN(y1,y2) .AND. y < MAX(y1,y2)) then
          if(ABS(x2-x1) > tol) then
            m = (y2 - y1)/(x2 - x1)
            if((m >  0 .AND. y1 + m*(x-x1) > y).OR.(m <= 0 .AND. y1 + m*(x-x1) < y))         &
              inside = .NOT.inside
          else
            if(x > x1) inside = .NOT.inside
          end if
        end if
      end associate
    end do
  end function inside

  subroutine shapeout(bndry,flname)
    class(boundary) ,intent(in) :: bndry
    character(len=*),intent(in) :: flname
    integer :: i

    open(10,file=flname)
    do i = 1,bndry%Ncn
      write(10,'(4F16.9)') bndry%r_nod(:,i),bndry%tau(:,i)
    end do
    write(10,'(4F16.9)') bndry%r_nod(:,1),bndry%tau(:,1)
    close(10)
  end subroutine shapeout
end module solid_boundary
