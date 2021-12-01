module contours
  use iso_fortran_env, only: rp => REAL64
  use vectors
  implicit none

  integer      ,private,save :: buff
  real(kind=rp),private,save :: h

  type :: contour
    integer ::                                                   &
      n_pts = 0            ,& ! no. of points in the contours
      n_mem = 0               ! max. no. of points w/o reallocation
    
    type(vector) ,allocatable ::                                 &
      r(:)                 ,& ! nodal points
      rc(:)                ,& ! centre of segment
      n(:)                 ,& ! normal to segment
      t(:)                 ,& ! tangent to segment
      rl(:)                   ! location of shear layer
    real(kind=rp),allocatable :: ds(:)
    
  contains
    procedure :: set_via_array
    procedure :: set_via_func
    generic   :: set => set_via_array, set_via_func
!   procedure :: move_contour
    procedure :: shapeout
  end type contour

contains
  subroutine set_via_func(this,contr,n_pts,u_max,dn,buffer)
    class(contour),intent(inout)       :: this
    real(kind=rp) ,intent(in)          :: dn
    real(kind=rp) ,intent(in),optional :: u_max
    integer       ,intent(in)          :: n_pts
    integer       ,intent(in),optional :: buffer
    real(kind=rp) :: du
    integer       :: i
    type(vector)  :: grid(n_pts)

    interface
      elemental function contr(u)
        use iso_fortran_env, only: rp => REAL64
        use vectors
        real(kind=rp),intent(in) :: u
        type(vector) :: contr
      end function contr
    end interface

    du = 1._rp/n_pts
    if(PRESENT(u_max)) du = du*u_max
    grid = contr([(du*i, i = 0,n_pts-1)])
    if(PRESENT(buffer)) then
      call set_via_array(this,grid,dn,buffer)
    else
      call set_via_array(this,grid,dn)
    end if
  end subroutine set_via_func

  subroutine set_via_array(this,grid,dn,buffer)
    class(contour),intent(inout)       :: this
    type(vector)  ,intent(in)          :: grid(:)
    real(kind=rp) ,intent(in)          :: dn
    integer       ,intent(in),optional :: buffer

    h = dn
    if(PRESENT(buffer)) then
      buff = buffer
    else
      buff = 0
    end if
    this%n_pts = SIZE(grid)
    this%n_mem = (this%n_pts*(100+buff))/100

    if(this%n_pts > this%n_mem) then
      DEALLOCATE(this%r)
      DEALLOCATE(this%rc)
      DEALLOCATE(this%n)
      DEALLOCATE(this%rl)
      DEALLOCATE(this%t)
      DEALLOCATE(this%ds)
    end if
    if(.NOT.ALLOCATED(this%r)) then
      this%n_mem = (this%n_pts*(100+buff))/100
      ALLOCATE(this%r (this%n_mem))
      ALLOCATE(this%rc(this%n_mem))
      ALLOCATE(this%n (this%n_mem))
      ALLOCATE(this%t (this%n_mem))
      ALLOCATE(this%rl(this%n_mem))
      ALLOCATE(this%ds(this%n_mem))
    end if

    associate(r  => this%r (1:this%n_pts),                       &
              rc => this%rc(1:this%n_pts),                       &
              n  => this%n (1:this%n_pts),                       &
              t  => this%t (1:this%n_pts),                       &
              rl => this%rl(1:this%n_pts),                       &
              ds => this%ds(1:this%n_pts))
      r  = grid
      rc = (CSHIFT(r,1) + r)/2._rp
      t  =  CSHIFT(r,1) - r
      ds = mag(t)
      t  = t/ds
      n  = -rot90(t)
      rl = rc + h*n
    end associate
  end subroutine set_via_array

  subroutine shapeout(this,flname)
    class(contour)  ,intent(in) :: this
    character(len=*),intent(in) :: flname
    integer :: i

    open(10,file=flname)
    do i = 1,this%n_pts
      write(10,*) this%r(i)
    end do
    write(10,*) this%r(1)
    close(10)
  end subroutine shapeout
end module contours
