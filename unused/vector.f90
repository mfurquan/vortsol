module vector
  ! Defines dynamically expandable vector array
  !
  use iso_fortran_env, only: rp => REAL64
  implicit none

  type :: vec
    integer :: n_sd, n_max, n
    real(kind=rp),allocatable :: e(:,:)
  contains
    generic :: operator(+) => add_vec
    generic :: operator(-) => sub_vec
    generic :: operator(*) => scal_vec
    generic :: assignment(=) => array_to_vec
  end type vec

contains
  pure subroutine array_to_vec(v,a)
    real(kind=rp) ,intent(in)    :: a(:,:)
    type(vec(*,*)),intent(inout) :: v
    integer :: n

    if(SIZE(a,1) /= v%n_sd) &
      error stop 'ERROR! Failed array to vec assign, unequal n_sd'
    if(SIZE(a,2) < v%n_max) DEALLOCATE(v%e)
    n = SIZE(a,2)
    if(ALLOCATED(v)) then
      if(SIZE(v) /= n) DEALLOCATE(v)
    end if
    if(.NOT.ALLOCATED(v)) ALLOCATE(v(n))

    do concurrent (i = 1:n)
      v(i)%e = a(:,i)
    end do
  end subroutine array_to_vec

  pure elemental function add_vec(a,b)
    type(vec),intent(in) :: a, b
    type(vec) :: add_vec

    add_vec%e = a%e + b%e
  end function add_vec

  pure elemental function sub_vec(a,b)
    type(vec),intent(in) :: a, b
    type(vec) :: sub_vec

    sub_vec%e = a%e - b%e
  end function sub_vec

  pure elemental function scal_vec(k,a)
    type(vec),    intent(in) :: a
    real(kind=rp),intent(in) :: k
    type(vec) :: scal_vec

    scal_vec%e = k*a%e
  end function scal_vec

  pure elemental function mag(a)
    type(vec),intent(in) :: a
    real(kind=rp) :: mag

    mag = NORM2(a%e)
  end function mag

  pure elemental function ang(a)
    type(vec(2)),intent(in) :: a
    real(kind=rp)           :: ang

    ang = atan(vec%e(2)/vec%e(1))
  end function ang

  pure function SUM_vec(a)
    type(vec),intent(in) :: a(:)
    type(vec) :: SUM_vec
    integer   :: i

    SUM_vec = zero_vec
    do concurrent (i = 1:SIZE(a))
      SUM_vec = SUM_vec + a(i)
    end do
  end function SUM_vec

  pure elemental function rot90(a)
    type(vec),intent(in) :: a
    type(vec) :: rot90

    if(n_sd == 2) then
      rot90%e(1) = -a%e(2); rot90%e(2) = a%e(1)
    else
      error stop 'vector::rot90 implemented only for n_sd=2'
    end if
  end function rot90
end module vector
