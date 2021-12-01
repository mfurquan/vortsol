module array_utils
  use iso_fortran_env, only: rp => REAL64
  implicit none
contains
  pure function map(func,ar)
    real(kind=rp),intent(in) :: ar(:,:)
    real(kind=rp) :: map(SIZE(ar,2))
    interface
      pure function func(a)
        use iso_fortran_env, only: rp => REAL64
        real(kind=rp),intent(in) :: a(:)
        real(kind=rp) :: func
      end function func
    end interface

    do concurrent (ia = 1:SIZE(ar,2))
      map(ia) = func(ar(:,ia))
    end do
  end function map
end module array_utils
