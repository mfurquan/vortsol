module vector
  use iso_fortran_env, only: rp => REAL64
  implicit none

  type :: vec(dimen)
    integer,len   :: dimen
    real(kind=rp) :: e(dimen)
  end type vec

  interface operator (+)
    module procedure :: add_vec
  end interface operator (+)

  interface operator (-)
    module procedure :: sub_vec
  end interface operator (-)

  interface operator (*)
    module procedure :: scal_vec
  end interface operator (*)

  interface assignment (=)
    module procedure :: array_to_vec, array_to_vecArray
  end interface assignment (=)
contains
  pure subroutine array_to_vec(v,a)
    real(kind=rp)     ,intent(in)    :: a(:)
    type(vec(SIZE(a))),intent(inout) :: v

    v%e = a
  end subroutine array_to_vec

  pure subroutine array_to_vecArray(v,a)
    real(kind=rp)                   ,intent(in)    :: a(:,:)
    type(vec(SIZE(a,1))),allocatable,intent(inout) :: v(:)
    integer :: i, n

    n = SIZE(a,2)

    if(ALLOCATED(v)) then
      if(SIZE(v) /= n) ERROR STOP                                &
        'ERROR: dimension mistmatch in array to vec assignment'
    end if
    if(.NOT.ALLOCATED(v)) ALLOCATE(v(n))

    do concurrent (i = 1:n)
      v(i)%e = a(:,i)
    end do
  end subroutine array_to_vecArray

  pure elemental function add_vec(a,b) result(c)
    type(vec(*)),intent(in) :: a, b
    type(vec(a%dimen))      :: c

    c%e = a%e + b%e
  end function add_vec

  pure elemental function sub_vec(a,b) result(c)
    type(vec(*)),intent(in) :: a, b
    type(vec(a%dimen))      :: c

    c%e = a%e - b%e
  end function sub_vec

  pure elemental function scal_vec(k,a) result(b)
    type(vec(*)) ,intent(in) :: a
    real(kind=rp),intent(in) :: k
    type(vec(a%dimen))       :: b

    b%e = k*a%e
  end function scal_vec

  pure elemental function mag(a)
    type(vec(*)),intent(in) :: a
    real(kind=rp)           :: mag

    mag = NORM2(a%e)
  end function mag

  pure elemental function mag2(a)
    type(vec(*)),intent(in) :: a
    real(kind=rp)           :: mag2

    mag2 = SUM(a%e**2)
  end function mag2

  pure elemental function rot90(a)
    type(vec(2)),intent(in) :: a
    type(vec(2))            :: rot90

    rot90%e(1) = -a%e(2); rot90%e(2) = a%e(1)
  end function rot90
end module vector
