module smart_array
  use global_param
  implicit none

  integer,parameter :: predict_rule = 3
  real   ,parameter :: over_compen  = 1.2

  type :: real_array
    integer                   :: rank, N, Nmax
    integer,      allocatable :: dimen(:)
    real(kind=rp),allocatable :: val(:)
  contains
  end type real_array

contains
  pure subroutine init_with_array(this,dimen,array)
    class(real_array),intent(inout) :: this
    integer          ,intent(in)    :: dimen
    real(kind=rp)    ,intent(in)   ,optional :: array(:)

    this%dimen = dimen
    this%rank  = SIZE(this%dimen)
    this%N     = this%dimen(-1)
    this%Nmax  = over_compen*this%N
    if(PRESENT(array)) then
      this%val = array
    else
      ALLOCATE(this%val(this%dimen))
  end subroutine init_with_array
end module smart_array
