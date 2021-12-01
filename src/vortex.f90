module vortex
  use iso_fortran_env, only: rp => REAL64
  use vectors
  use contours
  use lapack95, only: gesv
  implicit none

  integer      ,parameter :: n_sd = 2
  real(kind=rp),parameter :: pi = acos(-1._rp)

  type :: vortexList
    integer                                :: maxnum, num = 0
    real(kind=rp)                          :: R2_rnk, c = 1.2_rp
    real(kind=rp),allocatable,dimension(:) :: Gama, eps
    type(vector) ,allocatable,dimension(:) :: r

  contains
    procedure :: init
    procedure :: move
    procedure :: delete
    procedure :: farthest
    procedure :: swap
    procedure :: gen_shear
    procedure :: Vel
    procedure :: Vel_ind
    procedure :: write_vor
  end type vortexList

contains
  pure subroutine init(vor,n_vor,R_rnk,adjstmnt_eps)
    class(vortexList),intent(inout)          :: vor
    integer          ,intent(in)             :: n_vor
    real(kind=rp)    ,intent(in)             :: R_rnk
    real(kind=rp)    ,intent(in)   ,optional :: adjstmnt_eps

    vor%maxnum = n_vor
    allocate(vor%Gama(n_vor))
    allocate(vor%eps(n_vor))
    allocate(vor%r(n_vor))

    vor%R2_rnk = R_rnk*R_rnk
    if(PRESENT(adjstmnt_eps)) vor%c = adjstmnt_eps
  end subroutine init

  pure subroutine move(vor,cn,U,Re,dt)
    class(vortexList),intent(inout) :: vor
    type(contour)    ,intent(in)    :: cn
    type(vector)     ,intent(in)    :: U
    real(kind=rp)    ,intent(in)    :: Re, dt
    type(vector)  :: I2, I3, r_ij, V_ind
    real(kind=rp) :: eps_i, I0, I1, r, tmp
    integer       :: i, j

    eps_i = vor%eps(1)
    do concurrent (i = 1:vor%num)
      I0 = 2._rp*pi*vor%eps(i)**2
      I1 = 0._rp
      I2 = O
      I3 = O
      V_ind = U
      do concurrent(j = 1:vor%num, i /= j)
        r_ij = vor%r(i) - vor%r(j)
        r    = mag(r_ij)

        V_ind = V_ind + (vor%Gama(j)/MAX(r*r,vor%R2_rnk))        &
                      * rot90(r_ij)
        tmp   = vor%Gama(j)*exp(-r/vor%eps(i))
        I1    = I1 + tmp
        I2    = I2 + (tmp/(vor%eps(i)*r))*r_ij

        if(r < eps_i) eps_i = r
      end do
      do concurrent (j = 1:cn%n_pts)
        r_ij = vor%r(i) - cn%rc(j)
        r    = mag(r_ij)
        tmp  = exp(-r/vor%eps(i))
        I0   = I0 + (vor%eps(i)*(cn%n(j).dot.r_ij)               &
                  * (r+vor%eps(i))/r**2)*tmp
        I3   = I3 + (cn%ds(j)*tmp)*cn%n(j)
      end do
      V_ind = (dt/2._rp/pi)*V_ind + (dt/Re/I0)*I3 - (dt/Re/I1)*I2
      call vor%r(i)%incr_by(V_ind)
      vor%eps(i) = vor%c*eps_i
    end do
  end subroutine move

  pure subroutine delete(vor,ndel)
    class(vortexList),intent(inout) :: vor
    integer          ,intent(in)    :: ndel

    call vor%swap(vor%farthest(ndel))
  end subroutine delete

  pure function farthest(vor,n) result(f)
    class(vortexList),intent(in) :: vor
    integer          ,intent(in) :: n
    integer       :: f(n), i, j
    real(kind=rp) :: r2(n), r2_i

    r2 = 0._rp
    do i = 1,vor%num
      r2_i = mag2(vor%r(i))
      do j = 1,n
        if(r2_i > r2(j)) then
          call insert(f,r2,r2_i,i,j)
          exit
        end if
      end do
    end do

  contains
    pure subroutine insert(far_vor,dist,new_dist,new,pos)
      integer      ,intent(in)    :: new, pos
      real(kind=rp),intent(in)    :: new_dist
      integer      ,intent(inout) :: far_vor(n)
      real(kind=rp),intent(inout) :: dist(n)

      far_vor(pos+1:n) = far_vor(pos:n-1)
      dist(pos+1:n)    = dist(pos:n-1)

      far_vor(pos)     = new
      dist(pos)        = new_dist
    end subroutine insert
  end function farthest

  pure subroutine swap(vor,list)
    class(vortexList),intent(inout) :: vor
    integer          ,intent(in)    :: list(:)
    real(kind=rp) :: Gama_tmp, eps_tmp
    type(vector)  :: r_tmp
    integer       :: i

    do concurrent (i = 1:SIZE(list))
      Gama_tmp = vor%Gama(i)
      eps_tmp  = vor%eps(i)
      r_tmp    = vor%r(i)

      vor%Gama(i) = vor%Gama(list(i))
      vor%eps(i)  = vor%eps(list(i))
      vor%r(i)    = vor%r(list(i))

      vor%Gama(list(i)) = Gama_tmp
      vor%eps(list(i))  = eps_tmp
      vor%r(list(i))    = r_tmp
    end do
  end subroutine swap

  subroutine gen_shear(vor,cn,U)
    class(vortexList),intent(inout) :: vor
    type(contour)    ,intent(in)    :: cn
    type(vector)     ,intent(in)    :: U
    real(kind=rp) :: A(cn%n_pts,cn%n_pts), b(cn%n_pts)
    integer       :: i, j, m, n

    ! calc. strength of shear layer vortices via no-slip
    b = -vor%Vel(cn%r(1:cn%n_pts),U).dot.cn%t(1:cn%n_pts)
    do concurrent (i = 1:cn%n_pts, j = 1:cn%n_pts)
      A(i,j) = coeff(i,j)
    end do
    call gesv(A,b)

    !delete old vortices, if required
    m = 1
    if(vor%num + cn%n_pts <= vor%maxnum) then
      m       = vor%num + 1 
      vor%num = vor%num + cn%n_pts
    else
      call vor%delete(vor%num + cn%n_pts - vor%maxnum)
      vor%num = vor%maxnum
    end if
    n = m + cn%n_pts - 1

    ! generate shear layer vortices
    vor%Gama(m:n) = b
    vor%r   (m:n) = cn%rl(1:cn%n_pts)
    vor%eps (m:n) = cn%ds(1:cn%n_pts)*vor%c

  contains
    pure function coeff(p,q)
      integer,intent(in) :: p, q
      type(vector)  :: v
      real(kind=rp) :: coeff
      real(kind=rp) :: r2
  
      v     = cn%r(p) - cn%rl(q)
      r2    = mag2(v)
      coeff = rot90(v).dot.cn%n(p)
  
      if(r2 < vor%R2_rnk) then
        coeff = coeff/(2._rp*pi*vor%R2_rnk)
      else
        coeff = coeff/(2._rp*pi*r2)
      end if
    end function coeff
  end subroutine gen_shear

  pure function Vel(vor,coord,U)
    class(vortexList),intent(in) :: vor
    type(vector)     ,intent(in) :: coord(:), U
    type(vector) :: Vel(SIZE(coord))
    integer      :: i, j

    Vel = U
    do concurrent (i = 1:SIZE(coord), j = 1:vor%num)
      Vel(i) = Vel(i) + vor%Vel_ind(j,coord(i))
    end do
  end function Vel

  pure function Vel_ind(vor,i_vor,at_pos)
    class(vortexList),intent(in) :: vor
    integer          ,intent(in) :: i_vor
    type(vector)     ,intent(in) :: at_pos
    type(vector)                 :: Vel_ind
    real(kind=rp)                :: r2

    Vel_ind = rot90(at_pos - vor%r(i_vor))
    r2      = mag2(Vel_ind)

    if(r2 < vor%R2_rnk) then
      Vel_ind = (vor%Gama(i_vor)/(2._rp*pi*vor%R2_rnk))*Vel_ind
    else
      Vel_ind = (vor%Gama(i_vor)/(2._rp*pi*r2)        )*Vel_ind
    end if
  end function Vel_ind

  subroutine write_vor(vor,flname,flnum)
    class(vortexList),intent(in)          :: vor
    character(len=*) ,intent(in)          :: flname
    integer          ,intent(in),optional :: flnum
    character(len=7) :: stamp
    integer          :: i

    if(PRESENT(flnum)) then
      write(stamp,'(i3.3a4)') flnum,'.dat'
      open(unit=10,file=flname//stamp)
    else
      open(unit=10,file=flname//'.dat')
    end if

    do i = 1,vor%num
      write(10,*) vor%r(i),vor%Gama(i)
    end do
    close(10)
  end subroutine write_vor

end module vortex
