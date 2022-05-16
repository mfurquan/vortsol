module utility
  use global_param, only: rp, n_sd
  implicit none

  interface rot90
    module procedure rot90_1, rot90_n
  end interface rot90

  interface mag2
    module procedure mag2_1, mag2_n
  end interface mag2

  interface operator (/)
    module procedure divide_comp
  end interface operator (/)

  interface assignment (=)
    module procedure array2d_to_1d
  end interface assignment(=)

contains
  pure function rot90_1(v)
    real(kind=rp),intent(in) :: v(:)
    real(kind=rp)            :: rot90_1(n_sd)

    rot90_1(1) = -v(2); rot90_1(2) = v(1)
  end function rot90_1

  pure function rot90_n(v)
    real(kind=rp),intent(in) :: v(:,:)
    real(kind=rp)            :: rot90_n(n_sd,SIZE(v,2))

    rot90_n(1,:) = -v(2,:); rot90_n(2,:) = v(1,:)
  end function rot90_n

  pure function mag2_1(v)
    real(kind=rp),intent(in) :: v(:)
    real(kind=rp) :: mag2_1

    mag2_1 = v(1)**2 + v(2)**2
  end function mag2_1

  pure function mag2_n(v,mask)
    real(kind=rp),intent(in)          :: v(:,:)
    logical      ,intent(in),optional :: mask(:)
    real(kind=rp),allocatable         :: mag2_n(:)
    integer                           :: i, n

    if(PRESENT(mask)) then
      ALLOCATE(mag2_n(COUNT(mask)))
      n = 0
      do i = 1,SIZE(v,2)
        if(mask(i)) then
          n = n + 1
          mag2_n(n) = v(1,i)**2 + v(2,i)**2
        end if
      end do
    else
      ALLOCATE(mag2_n(SIZE(v,2)))
      mag2_n = v(1,:)**2 + v(2,:)**2
    end if
  end function mag2_n

  pure function divide_comp(v,c)
    real(kind=rp),intent(in) :: v(:,:), c(:)
    real(kind=rp)            :: divide_comp(n_sd,SIZE(v,2))

    divide_comp(1,:) = v(1,:)/c
    divide_comp(2,:) = v(2,:)/c
  end function divide_comp

  pure subroutine array2d_to_1d(a1d,a2d)
    real(kind=rp),intent(inout) :: a1d(:)
    real(kind=rp),intent(in)    :: a2d(:,:)
    integer                     :: i, j

    do concurrent (i = 1:n_sd, j = 1:SIZE(a2d,2))
      a1d((j-1)*n_sd+i) = a2d(i,j)
    end do
  end subroutine array2d_to_1d

  pure function largestloc(a,k)
    real(kind=rp),intent(in) :: a(:)
    integer      ,intent(in) :: k
    integer                  :: largestloc(k), i, j

    largestloc = [(i, i = 1,k)]
    call quicksort(largestloc,larger)
    do i = k+1,SIZE(a)
      do j = 1,k
        if(a(i) > a(largestloc(j))) then
          largestloc(j) = i
          continue
        end if
      end do
    end do

    contains
      pure function larger(p,q)
        integer,intent(in) :: p, q
        logical            :: larger

        larger = a(p) > a(q)
      end function larger
  end function largestloc

  pure subroutine quicksort(a,crit)
    integer,intent(inout) :: a(:)
    interface
      pure function crit(p,q)
        integer,intent(in) :: p, q
        logical            :: crit
      end function crit
    end interface

    call qsort(a)

  contains
    recursive pure subroutine qsort(ar)
      integer,intent(inout) :: ar(:)
      integer               :: pIndex
      
      if(SIZE(ar)>0) then
        call partition(pIndex,ar)
        call qsort(ar(:pIndex-1))
        call qsort(ar(pIndex+1:))
      end if
    end subroutine qsort

    pure subroutine partition(pIndex,ar)
      integer,intent(inout) :: pIndex, ar(:)
      integer               :: pivot, last, tmp, i

      last   = SIZE(ar)
      pivot  = ar(last)
      pIndex = 1
      do i = 1,last-1
        if(crit(ar(i),pivot)) then
          ! swap
          tmp        = ar(i)
          ar(i)      = ar(pIndex)
          ar(pIndex) = tmp

          pIndex     = pIndex + 1
        end if
      end do
      ! swap
      tmp        = ar(pIndex)
      ar(pIndex) = ar(last)
      ar(last)   = tmp
    end subroutine partition
  end subroutine quicksort
end module utility
