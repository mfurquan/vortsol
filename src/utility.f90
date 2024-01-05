module utility
  use global_param, only: rp, n_sd
  implicit none

  interface rot90
    module procedure rot90_1, rot90_n
  end interface rot90

  interface mag2
    module procedure mag2_1, mag2_n
  end interface mag2

  interface mag
    module procedure mag_1, mag_n
  end interface mag

  interface operator (/)
    module procedure divide_comp
  end interface operator (/)

  interface operator (.cross.)
    module procedure cross_prod_1
  end interface operator (.cross.)

  interface cross_k
    module procedure cross_k_1, cross_k_n
  end interface cross_k

  interface operator (.dot.)
          module procedure dot_prod_1, dot_prod_n
  end interface operator (.dot.)

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

  pure function cross_k_1(v)
    real(kind=rp),intent(in) :: v(:)
    real(kind=rp)            :: cross_k_1(n_sd)

    cross_k_1(1) = v(2); cross_k_1(2) = -v(1)
  end function cross_k_1

  pure function cross_k_n(v)
    real(kind=rp),intent(in) :: v(:,:)
    real(kind=rp)            :: cross_k_n(n_sd,SIZE(v,2))

    cross_k_n(1,:) = v(2,:); cross_k_n(2,:) = -v(1,:)
  end function cross_k_n

  pure function k_cross(v)
    real(kind=rp),intent(in) :: v(:)
    real(kind=rp)            :: k_cross(n_sd)

    k_cross(1) = -v(2); k_cross(2) = v(1)
  end function k_cross

  pure function cross_prod_1(u,v)
    real(kind=rp),intent(in) :: u(n_sd), v(n_sd)
    real(kind=rp)            :: cross_prod_1

    cross_prod_1 = u(1)*v(2) - u(2)*v(1)
  end function cross_prod_1

  pure function dot_prod_1(u,v)
    real(kind=rp),intent(in) :: u(n_sd), v(n_sd)
    real(kind=rp)            :: dot_prod_1

    dot_prod_1 = u(1)*v(1) + u(2)*v(2)
  end function dot_prod_1

  pure function dot_prod_n(u,v)
    real(kind=rp),intent(in) :: u(:,:), v(:,:)
    real(kind=rp)            :: dot_prod_n(SIZE(u,2))
    integer                  :: i

    do concurrent (i = 1:SIZE(u,2))
      dot_prod_n(i) = u(1,i)*v(1,i) + u(2,i)*v(2,i)
    end do
  end function dot_prod_n

  pure function mag_1(v)
    real(kind=rp),intent(in) :: v(:)
    real(kind=rp) :: mag_1

    mag_1 = NORM2(v)
  end function mag_1

  pure function mag_n(v)
    real(kind=rp),intent(in) :: v(:,:)
    real(kind=rp) :: mag_n(SIZE(v,2))

    mag_n = NORM2(v,1)
  end function mag_n

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
    integer                  :: i

    !divide_comp(1,:) = v(1,:)/c
    !divide_comp(2,:) = v(2,:)/c
    do concurrent (i = 1:SIZE(v,2))
      divide_comp(:,i) = v(:,i)/c(i)
    end do
  end function divide_comp

  pure subroutine array2d_to_1d(a1d,a2d)
    real(kind=rp),intent(inout) :: a1d(:)
    real(kind=rp),intent(in)    :: a2d(:,:)
    integer                     :: i, j

    do concurrent (i = 1:n_sd, j = 1:SIZE(a2d,2))
      a1d((j-1)*n_sd+i) = a2d(i,j)
    end do
  end subroutine array2d_to_1d

  pure subroutine swap(a,b)
    real(kind=rp),intent(inout) :: a, b
    real(kind=rp) :: tmp

    tmp = a; a = b; b = tmp
  end subroutine swap

  pure function same_sign(a,b)
    real(kind=rp),intent(in) :: a, b
    logical :: same_sign

    same_sign = (a > 0 .AND. b > 0).OR.(a < 0 .AND. b < 0)
  end function same_sign

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
          exit
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
