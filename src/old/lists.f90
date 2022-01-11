module lists
  use iso_fortran_env, only: rp => REAL64
  implicit none

  integer,parameter :: n_sd = 2
  real   ,parameter :: buffer = 1.2

  type :: scalars
    integer :: N = 0, Nmax = 0
    real(kind=rp),allocatable :: val(:)

  contains
    procedure :: array_to_scalars
    generic   :: assignment(=) => array_to_scalars
  end type scalars

  type :: vectors
    integer :: N = 0, Nmax = 0
    real(kind=rp),allocatable :: val(:,:)

  contains
    procedure :: array_to_vctors
    procedure :: add_vectors
    generic   :: assignment(=) => array_to_vectors
    generic   :: operator(+)   => add_vectors
  end type vectors

contains
  subroutine array_to_scalars(s,a)
    class(scalars),intent(inout) :: s
    real(kind=rp) ,intent(in)    :: a(:)

    s%N = SIZE(a)
    if(this%N > this%Nmax) then
      DEALLOCATE(s%val)
      s%Nmax = buffer*s%N
    end if
    if(.NOT.ALLOCATED(s%val)) ALLOCATE(s%val(s%Nmax))
    s%val(1:s%N) = a
  end subroutine array_to_scalars

  subroutine array_to_vectors(v,a)
    class(vectors),intent(inout) :: v
    real(kind=rp) ,intent(in)    :: a(:)

    v%N = SIZE(a,2)
    if(this%N > this%Nmax) then
      DEALLOCATE(v%val)
      v%Nmax = buffer*v%N
    end if
    if(.NOT.ALLOCATED(v%val)) ALLOCATE(v%val(n_sd,s%Nmax))
    v%val(1:n_sd,1:v%N) = a
  end subroutine array_to_vectors

  pure function CSHIFT(u,step) result (v)
    type(vectors),intent(in) :: u
    integer      ,intent(in) :: step
    type(vector)             :: v
    real(kind=rp)            :: tmp(n_sd,step)

    v%N = u%N; v%Nmax = v%Nmax
    ALLOCATE(v%val(Nmax))

    tmp                  = u%val(:,N-step+1:N)
    u%val(:,step+1:N)    = u%val(:,       1:N-step)
    u%val(:,     1:step) = tmp
  end function CSHIFT

  pure function add_vectors(u,v) result(w)
    class(vectors),intent(in) :: u, v
    type(vectors)             :: w

    w%N = u%N; w%Nmax = u%Nmax
    ALLOCATE(w%val(Nmax))

    w%val(1:)
  end function add_vectors
end module lists
