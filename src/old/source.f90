module source
  use globalParam
  use vectors
  implicit none

  type :: sourceList
    integer                   :: maxnum, num = 0
    real(kind=rp),allocatable :: m(:)
    type(vector) ,allocatable :: r(:)
  contains
    procedure :: init
  end type sourceList
end module source
