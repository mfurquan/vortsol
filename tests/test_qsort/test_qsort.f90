program test_qsort
  use utility, only: quicksort
  implicit none

  integer :: N, k
  real   ,allocatable :: U(:)
  integer,allocatable :: I(:)

  read(*,*) N
  ALLOCATE(U(N))
  ALLOCATE(I(N))
  call RANDOM_NUMBER(U)
  U = U*10.0
  !N = 3
  !U = [2.0, 1.0, 3.0]
  write(*,*) 'Unsorted:',U
  I = [(k, k = 1,N)]
  call quicksort(I,gt)
  write(*,*) 'Sorted:',U(I)
contains
  pure function gt(a,b)
    integer,intent(in) :: a, b
    logical :: gt

    gt = U(a) > U(b)
  end function gt
end program test_qsort
