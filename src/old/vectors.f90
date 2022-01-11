module vectors
  use iso_fortran_env, only: rp => REAL64
  implicit none

  integer,parameter,private :: dimen = 2

  type :: vector
    real(kind=rp) :: e(dimen) = 0._rp
  contains
    procedure :: add_to
    generic   :: incr_by => add_to
  end type vector

  interface vector
    module procedure :: const_vector
  end interface vector

  integer,private :: i

  type(vector),parameter :: O = vector([(0._rp, i = 1,dimen)])

  interface operator (+)
    module procedure :: add_vector
  end interface operator (+)

  interface operator (-)
    module procedure :: sub_vector, neg_vector
  end interface operator (-)

  interface operator (*)
    module procedure :: scal_mult
  end interface operator (*)

  interface operator (/)
    module procedure :: scal_div
  end interface operator (/)

  interface operator (.dot.)
    module procedure :: dot_prod
  end interface operator (.dot.)

  interface assignment (=)
    module procedure :: array_to_vector, array_to_vectorArray,   &
                        vectorArray_to_array
  end interface assignment (=)

contains
  pure subroutine array_to_vector(v,a)
    real(kind=rp),intent(in)    :: a(:)
    type(vector) ,intent(inout) :: v

    v%e = a
  end subroutine array_to_vector

  pure subroutine array_to_vectorArray(v,a)
    real(kind=rp),intent(in)    :: a(:,:)
    type(vector) ,intent(inout) :: v(:)
    integer :: i, n

    n = SIZE(a,2)

    if(SIZE(v) /= n) ERROR STOP                                  &
      'ERROR: dimension mistmatch in array to vec assignment'

    do concurrent (i = 1:n)
      v(i)%e = a(:,i)
    end do
  end subroutine array_to_vectorArray

  pure subroutine vectorArray_to_array(a,v)
    type(vector)             ,intent(in)    :: v(:)
    real(kind=rp),allocatable,intent(inout) :: a(:,:)
    integer :: i, n

    n = SIZE(v)

    if(ALLOCATED(a)) then
      if(SIZE(v) /= n) ERROR STOP                                &
        'ERROR: dimension mistmatch in array to vec assignment'
    else
      ALLOCATE(a(dimen,n))
    end if

    do concurrent (i = 1:n)
      a(:,i) = v(i)%e
    end do
  end subroutine vectorArray_to_array

  pure elemental function add_vector(a,b) result(c)
    type(vector),intent(in) :: a, b
    type(vector)            :: c

    c%e = a%e + b%e
  end function add_vector

  pure elemental function sub_vector(a,b) result(c)
    type(vector),intent(in) :: a, b
    type(vector)            :: c

    c%e = a%e - b%e
  end function sub_vector

  pure elemental function neg_vector(a) result(c)
    type(vector),intent(in) :: a
    type(vector)            :: c

    c%e = -a%e
  end function neg_vector

  pure elemental function scal_mult(k,a) result(b)
    type(vector) ,intent(in) :: a
    real(kind=rp),intent(in) :: k
    type(vector)             :: b

    b%e = k*a%e
  end function scal_mult

  pure elemental function scal_div(a,k) result(b)
    type(vector) ,intent(in) :: a
    real(kind=rp),intent(in) :: k
    type(vector)             :: b

    b%e = a%e/k
  end function scal_div

  pure elemental function dot_prod(a,b) result(c)
    type(vector),intent(in) :: a, b
    real(kind=rp) :: c

    c = DOT_PRODUCT(a%e,b%e)
  end function dot_prod

  pure elemental subroutine add_to(a,b)
    class(vector),intent(inout) :: a
    type(vector) ,intent(in)    :: b

    a%e = a%e + b%e
  end subroutine add_to

  pure function const_vector(a) result(v)
    real(kind=rp),intent(in) :: a
    type(vector) :: v

    v%e = a
  end function const_vector

  pure elemental function mag(a)
    type(vector),intent(in) :: a
    real(kind=rp)           :: mag

    mag = NORM2(a%e)
  end function mag

  pure elemental function mag2(a)
    type(vector),intent(in) :: a
    real(kind=rp)           :: mag2

    mag2 = SUM(a%e*a%e)
  end function mag2

  pure elemental function dist(a,b)
    type(vector),intent(in) :: a, b
    real(kind=rp)           :: dist

    dist = NORM2(a%e - b%e)
  end function dist

  pure elemental function dist2(a,b)
    type(vector),intent(in) :: a, b
    real(kind=rp)           :: dist2

    dist2 = NORM2(a%e - b%e)
  end function dist2

  pure elemental function rot90(a)
    type(vector),intent(in) :: a
    type(vector)            :: rot90

    rot90%e(1) = -a%e(2); rot90%e(2) = a%e(1)
  end function rot90
end module vectors
