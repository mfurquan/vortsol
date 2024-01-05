module flow
  use global_param
  use utility
  use solid_boundary
  use flow_elements
  use lapack95, only: getrf, getrs
  implicit none

  private :: decom

  type :: flow_field
    integer :: N = 0, Nmax
    real(kind=rp)             :: R2_rnk, epsmax2
    real(kind=rp),allocatable :: Gama(:), r(:,:), V_vor(:,:)
    logical    ,allocatable   :: has_vortex(:), not_marker(:)

  contains
    procedure :: init
    procedure :: step
    procedure :: import_vor
    procedure :: calc_V_vor
    procedure :: move_blobs
    procedure :: unmark_vor
    procedure :: clean_vor
    procedure :: Vel
    procedure :: enforce_bc
    procedure :: purge_vor
    procedure :: same_orient
    procedure :: write_vor
    procedure :: decom
  end type flow_field
      
contains
  pure subroutine init(ffield,U,Nvor,R_rnk,bndry,V_K)
    class(flow_field)         ,intent(inout) :: ffield
    real(kind=rp)             ,intent(in)    :: U(n_sd)
    integer                   ,intent(in)    :: Nvor
    real(kind=rp)             ,intent(in)    :: R_rnk
    type(boundary)   ,optional,intent(inout) :: bndry
    real(kind=rp)    ,optional,intent(in)    :: V_K(:,:)

    ffield%Nmax = Nvor
    ALLOCATE(ffield%Gama(Nvor))
    ALLOCATE(ffield%r(n_sd,Nvor))
    ALLOCATE(ffield%V_vor(n_sd,Nvor))
    ALLOCATE(ffield%has_vortex(Nvor))
    ALLOCATE(ffield%not_marker(Nvor))

    ffield%has_vortex = .FALSE.
    ffield%not_marker = .TRUE.
    ffield%R2_rnk  = R_rnk**2
    ffield%epsmax2 = epslim2*ffield%R2_rnk

    if(PRESENT(bndry)) then
      if(PRESENT(V_K)) then
        call ffield%enforce_bc(bndry,U,V_K)
      else
        call ffield%enforce_bc(bndry,U)
      end if
    end if
  end subroutine init

  pure subroutine decom(ffield)
    class(flow_field),intent(inout) :: ffield
      DEALLOCATE(ffield%Gama)
      DEALLOCATE(ffield%r)
      DEALLOCATE(ffield%V_vor)
      DEALLOCATE(ffield%has_vortex)
      DEALLOCATE(ffield%not_marker)
  end subroutine decom

  subroutine step(ffield,U,Re,dt,bndry,V_K)
    class(flow_field)         ,intent(inout) :: ffield
    real(kind=rp)             ,intent(in)    :: U(n_sd), Re, dt
    type (boundary)  ,optional,intent(inout) :: bndry
    real(kind=rp)    ,optional,intent(in)    :: V_K(:,:)

    if(PRESENT(bndry)) then
      if(ffield%N + bndry%Ncn > ffield%Nmax)                                                  &
        call ffield%purge_vor(ffield%N + bndry%Ncn - ffield%Nmax)
      call ffield%import_vor(bndry)
      call ffield%calc_V_vor(Re,U,dt,bndry)
      call ffield%move_blobs(dt)
      call ffield%unmark_vor(bndry)
      call ffield%clean_vor(bndry)
      if(PRESENT(V_K)) then
        call ffield%enforce_bc(bndry,U,V_K)
      else
        call ffield%enforce_bc(bndry,U)
      end if
    else
      call ffield%calc_V_vor(Re,U,dt)
      call ffield%move_blobs(dt)
    end if
  end subroutine step

  subroutine calc_V_vor(ffield,Re,U,dt,bndry)
    class(flow_field)         ,intent(inout) :: ffield
    real(kind=rp)             ,intent(in)    :: Re, U(n_sd), dt
    type(boundary)   ,optional,intent(in)    :: bndry

    real(kind=rp),dimension(n_sd) :: dr, dr1, dr2, h, V_pot, I2, I3, d, s, s0
    real(kind=rp)                 :: r, r1, r2, eps, tmp, tmp1, tmp2, I1, I0, z0, alpha, beta
    integer                       :: i, j, jp

    do concurrent (i = 1:ffield%Nmax, ffield%has_vortex(i))
      eps = calc_eps(i)

      V_pot = U
      I1 = 0._rp
      if(ANY(.NOT.ffield%not_marker)) I1 = 1._rp
      I2 = 0._rp
      do j = 1,ffield%Nmax
        if(ffield%has_vortex(j).AND.ffield%not_marker(j).AND.i/=j) then
          dr = ffield%r(:,i) - ffield%r(:,j)
          r  = NORM2(dr)
          if(ffield%not_marker(i))                                                            &
            V_pot = V_pot + k_cross(dr)*ffield%Gama(j)/(2._rp*pi*MAX(r**2,ffield%R2_rnk))
          if(ffield%same_orient(i,j)) then
            tmp = ffield%Gama(j)*exp(-r/eps)
            I1  = I1 + tmp
            I2  = I2 + dr*tmp/(eps*r)
          end if
        end if
      end do
      
      I0 = 2._rp*pi*eps**2
      I3 = 0._rp
      if(PRESENT(bndry)) then
        if(ffield%N == bndry%Ncn) I1 = 1._rp ! I think this has been taken care of earlier
        do j = 1,bndry%Ncn
          dr  = ffield%r(:,i) - bndry%r_col(:,j)
          r   = NORM2(dr)
!          if(r < 20*eps) then
!            if(r > 3*bndry%ds(j)) then
              ! Mid-point rule
!              tmp = exp(-r/eps)
!              I3  = I3 + bndry%ds(j)*tmp*bndry%n(j)
!              I0  = I0 + eps*bndry%n_dot(dr,j)*(r+eps)*bndry%ds(j)*tmp/r**2
!          else !if(r > 10*bndry%ds(j)) then
            ! 3 point Gauss-Legendre rule
            h = sqrt(3._rp/5._rp)*bndry%tau(:,j)*bndry%ds(j)/2._rp
            dr1 = dr - h; r1 = NORM2(dr1); tmp1 = exp(-r1/eps)
            dr2 = dr + h; r2 = NORM2(dr2); tmp2 = exp(-r2/eps)
            tmp = exp(-r/eps)

            I3 = I3 + (5._rp*tmp1 + 8._rp*tmp + 5._rp*tmp2)*(bndry%ds(j)/18._rp)*bndry%n(j)
            I0 = I0 + bndry%n_dot(                                                   &
                        5._rp*dr1*(r1+eps)*tmp1/r1**2                                    &
                      + 8._rp*dr *(r +eps)*tmp /r**2                                     &
                      + 5._rp*dr2*(r2+eps)*tmp2/r2**2, j)*eps*bndry%ds(j)/18._rp
          !else
          !  I3 = 2._rp*eps*(1._rp - exp(-bndry%ds(j)/(2._rp*eps)))*bndry%n(j)
!            end if
!          end if
          if(ffield%not_marker(i)) then
            jp = j + 1
            if(j == bndry%Ncn) jp = 1
            d  = bndry%r_nod(:,jp) - bndry%r_nod(:,j)
            s  = ffield%r(:,i) - bndry%r_nod(:,jp)
            s0 = ffield%r(:,i) - bndry%r_nod(:,j)
            z0 = d.cross.s0
            alpha = atan((d.dot.s)/z0) - atan((d.dot.s0)/z0)
            beta  = 0.5_rp*log(mag2(s)/mag2(s0))
            V_pot = V_pot                                                            & 
                  + ((bndry%gama_att(j) + bndry%gama(j))                             &
                  * (alpha*d + beta*cross_k(d))                                      &
                  -  bndry%q_att(j)*(alpha*cross_k(d) - beta*d))/(2._rp*pi*mag(d))
          end if
        end do
      end if
      if(ffield%not_marker(i)) then
        ffield%V_vor(:,i) = V_pot + I2/(Re*I1) + I3/(Re*I0)
      else
        ffield%V_vor(:,i) = I2/(Re*I1) + I3/(Re*I0)
      end if
    end do

    !do concurrent (i = 1:bndry%Ncn)
    !  ffield%V_vor(:,bndry%marker(i)) = ffield%V_vor(:,bndry%marker(i))  &
    !                                  + bndry%gama(i)*bndry%tau(:,i)
    !end do

    contains
    pure function calc_eps(i)
      integer      ,intent(in) :: i
      real(kind=rp)            :: calc_eps, e(3), dr2
      integer                  :: j

      e = 1.e7

      do j = 1,ffield%Nmax
        if(j/=i .AND. ffield%has_vortex(j)) then
          dr2 = mag2(ffield%r(:,i) - ffield%r(:,j))
          if(dr2 < e(1)) then
            e(3) = e(2)
            e(2) = e(1)
            e(1) = dr2
          else if(dr2 < e(2)) then
            e(3) = e(2)
            e(2) = dr2
          else if(dr2 < e(3)) then
            e(3) = dr2
          end if
        end if
      end do

      !calc_eps = MAX(sqrt(e(2)),0.02)
      calc_eps = SQRT(SUM(e)/3._rp)
      !calc_eps = SQRT(MAX(SUM(e)/3._rp,2*dt/Re))
      !calc_eps = SQRT(MIN(MAX(SUM(e)/3._rp,0.02),1.0/4.0))
      !calc_eps = SQRT(MIN(MAX(SUM(e)/3._rp,ffield%epsmax2),1.0/4.0))
    end function calc_eps
  end subroutine calc_V_vor

  pure subroutine move_blobs(ffield,dt)
    class(flow_field),intent(inout) :: ffield
    real(kind=rp)    ,intent(in)    :: dt
    integer                         :: i

    do concurrent (i = 1:ffield%Nmax, ffield%has_vortex(i))
      ffield%r(:,i) = ffield%r(:,i) + dt*ffield%V_vor(:,i)
    end do
  end subroutine move_blobs

  pure subroutine unmark_vor(ffield,bndry)
    class(flow_field),intent(inout) :: ffield
    type(boundary)   ,intent(in)    :: bndry
    integer :: i

    do i = 1,bndry%Ncn
      !ffield%Gama(bndry%marker(i))       = bndry%gama(i)*bndry%ds(i)
      ffield%not_marker(bndry%marker(i)) = .TRUE.
    end do
  end subroutine unmark_vor

  pure subroutine clean_vor(ffield,bndry)
    class(flow_field) ,intent(inout) :: ffield
    type(boundary)    ,intent(in)    :: bndry
    integer :: i, j

    do concurrent (i = 1:ffield%Nmax, ffield%has_vortex(i))
      if(bndry%inside(ffield%r(:,i))) then
        ffield%has_vortex(i) = .FALSE.
        ffield%N = ffield%N - 1
      else
        do j = 1,ffield%Nmax
          !if(j/=i .AND. ffield%has_vortex(j) .AND. ffield%same_orient(i,j)) then
          if(j/=i .AND. ffield%has_vortex(j)) then
          !if(j/=i .AND. ffield%has_vortex(j) .AND. .NOT.ffield%same_orient(i,j)) then
            if(mag2(ffield%r(:,i)-ffield%r(:,j)) < ffield%R2_rnk/epslim2) then
              ffield%Gama(i) = ffield%Gama(i) + ffield%Gama(j)
              ffield%has_vortex(j) = .FALSE.
            end if
          end if
        end do
      end if
    end do
  end subroutine clean_vor

  pure function Vel(ffield,r,U,bndry)
    class(flow_field) ,intent(in) :: ffield
    real(kind=rp)     ,intent(in) :: r(:,:), U(n_sd)
    type(boundary)    ,intent(in),optional :: bndry
    real(kind=rp)                 :: Vel(n_sd,SIZE(r,2)), R2_rkn
    integer                       :: i

    R2_rkn    = ffield%R2_rnk
    
    do concurrent (i = 1:SIZE(r,2))
      Vel(:,i) = U + Vel_wake(r(:,i))
      if(PRESENT(bndry)) Vel(:,i) = Vel(:,i) + Vel_bndry(r(:,i))
    end do
  contains
    pure function Vel_wake(r)
      real(kind=rp),intent(in) :: r(n_sd)
      real(kind=rp)            :: Vel_wake(n_sd)
      integer                  :: i

      Vel_wake = 0._rp
      do concurrent(i = 1:ffield%N)
        Vel_wake = Vel_wake                                              &
                 + ffield%Gama(i)                                        &
                 * rkn_vortex(r,ffield%r(:,i),R2_rkn)
      end do
    end function Vel_wake

    pure function Vel_bndry(r)
      real(kind=rp),intent(in) :: r(n_sd)
      real(kind=rp)            :: Vel_bndry(n_sd)
      integer                  :: i, ip
      real(kind=rp)            :: d(n_sd), s(n_sd), s0(n_sd), z0, alpha, &
                                  beta

      Vel_bndry = 0._rp
      if(bndry%moving) then
        do concurrent (i = 1:bndry%Ncn)
          ip = i + 1
          if(i == bndry%Ncn) ip = 1
          d  = bndry%r_nod(:,ip) - bndry%r_nod(:,i)
          s  = bndry%r_nod(:,ip) - r
          s0 = bndry%r_nod(:,i)  - r
          z0 = d.cross.s0
          alpha = atan((d.dot.s)/z0) - atan((d.dot.s0)/z0)
          beta  = 0.5_rp*log(mag2(s)/mag2(s0))
          Vel_bndry = Vel_bndry                                          &
                    + ((bndry%gama_att(i) + bndry%gama(i))             &
                    * (alpha*d + beta*cross_k(d))                        &
                    + bndry%q_att(i)*(alpha*cross_k(d) - beta*d))        &
                    / (2._rp*pi*mag(d))
        end do
      end if
    end function Vel_bndry
  end function Vel

  pure subroutine enforce_bc(ffield,bndry,U,V_K)
    class(flow_field),intent(in)    :: ffield
    type(boundary)   ,intent(inout) :: bndry
    real(kind=rp)    ,intent(in)    :: U(n_sd)
    real(kind=rp)    ,intent(in),optional :: V_K(:,:)

    real(kind=rp)                   :: a_ij, b_ij
    integer                         :: i, ip, j, jp, w, info

    ! computing the strength of the attached vortex and source sheets
    if(bndry%moving) then
      bndry%gama_att = bndry%tau_dot(V_K)
      bndry%q_att    = bndry%  n_dot(V_K)
    end if

    ! computing the strength of the free vortex sheet
    do concurrent (i = 1:bndry%Ncn)
      ip = i + 1
      if(i == bndry%Ncn) ip = 1
      bndry%gama(i) = 0._rp
      do j = 1,bndry%Ncn
        if(i/=j) then
          jp = j + 1
          if(j == bndry%Ncn) jp = 1
          call calc_coeffs(a_ij,b_ij,                                                         &
                         ! ith panel:
                          bndry%r_nod(:,i),bndry%r_nod(:,ip),bndry%tau(:,i),                  &
                         ! jth panel:
                          bndry%r_nod(:,j),bndry%r_nod(:,jp),                                 &
                         ! neighborhood condition 1 (p1 = 0):
                          i==jp,                                                              &
                         ! neighborhood condition 2 (s2 = 0):
                          ip==j)

          if(bndry%var_coeff) bndry%A(i,j) = a_ij
          bndry%gama(i)  = bndry%gama(i) - a_ij*bndry%gama_att(j) - b_ij*bndry%q_att(j)
        else
          if(bndry%var_coeff) bndry%A(i,j) = -0.5_rp
        end if
      end do

      do w = 1,ffield%Nmax
        if(ffield%has_vortex(w).AND.ffield%not_marker(w)) then
          bndry%gama(i) = bndry%gama(i)                              &
                        - ffield%Gama(w)                               &
                        * bndry%tau_dot(V_iw_wake(bndry%r_nod(:,i),    &
                                                  bndry%r_nod(:,ip),   &
                                                  ffield%r(:,w)),i)
        end if
      end do
      
      if(bndry%moving) then
        bndry%gama(i) = bndry%gama(i) - bndry%tau_dot(U - V_K(:,i),i) - bndry%gama_att(i)/2._rp
      else
        bndry%gama(i) = bndry%gama(i) - bndry%tau_dot(U,i)
      end if
    end do

    if(bndry%var_coeff) then
      !bndry%A(:bndry%Ncn,bndry%Ncn+1) = bndry%ds
      bndry%A(bndry%Ncn+1,:bndry%Ncn) = bndry%ds
      ! Regularizing over-determined system:
      !bndry%A(bndry%Ncn+1,:bndry%Ncn)  = 1._rp
      bndry%A(:bndry%Ncn,bndry%Ncn+1)  = 1._rp
      bndry%A(bndry%Ncn+1,bndry%Ncn+1) = 0._rp
      ! LU factorization:
      call getrf(bndry%A(:bndry%Ncn+1,:bndry%Ncn+1),bndry%ipiv(:bndry%Ncn+1),info)

      if(.NOT.bndry%flexible) bndry%var_coeff = .FALSE.
    end if

    bndry%gama(bndry%Ncn+1) = 0._rp
    if(bndry%moving)  bndry%gama(bndry%Ncn+1) = SUM(bndry%tau_dot(V_K))

    call getrs(bndry%A(:bndry%Ncn+1,:bndry%Ncn+1),                       &
               bndry%ipiv(:bndry%Ncn+1),                                 &
               bndry%gama(:bndry%Ncn+1),'N',info)
  contains
    pure subroutine calc_coeffs(a_ij,b_ij,r_i,r_ip,tau_i,r_j,r_jp,zero_p1,zero_s2)
      real(kind=rp)                ,intent(inout) :: a_ij, b_ij
      real(kind=rp),dimension(n_sd),intent(in)    :: r_i, r_ip, tau_i, r_j, r_jp
      logical,                      intent(in)    :: zero_p1, zero_s2

      real(kind=rp),dimension(n_sd) :: d, d0, p1, p2, s1, s2, c1, c2, c3, u, v, P, Q
      real(kind=rp)                 :: z1, z2, z3, R

      if(zero_s2) then
        d  = r_i  - r_ip ; d0 = r_j - r_jp
        p1 = r_ip - r_j  ; p2 = r_i - r_j
        s1 = r_ip - r_jp ; s2 = r_i - r_jp
      else
        d  = r_ip - r_i  ; d0 = r_jp - r_j
        p1 = r_i  - r_jp ; p2 = r_ip - r_jp
        s1 = r_i  - r_j  ; s2 = r_ip - r_j
      end if
      
      z1 = p1 .cross. p2
      z2 = s1 .cross. s2
      z3 = s2 .cross. p2

      P = 0._rp
      c2 = (d0.dot.s1)*d + (d.dot.s1)*d0 - (d.dot.d0)*s1
      c3 = ( d.dot.d)*d0
      if(ABS(z2) >= tol) P  =     (atan(( d.dot.s2)/z2) - atan(( d.dot.s1)/z2))*c2
      if(ABS(z3) >= tol) P  = P + (atan((d0.dot.p2)/z3) - atan((d0.dot.s2)/z3))*c3

      Q = (log(mag2(s1)/mag2(s2))*c2 + log(mag2(p2)/mag2(s2))*c3)*0.5_rp

      if(.NOT.(zero_p1.OR.zero_s2)) then
        c1 = (d0.dot.p1)*d + (d.dot.p1)*d0 - (d.dot.d0)*p1
        if(ABS(z1) >= tol) P  = P + (atan(( d.dot.p1)/z1) - atan(( d.dot.p2)/z1))*c1
        Q  = Q + log(mag2(p2)/mag2(p1))*c1*0.5_rp
      end if

      R  = 2._rp*pi*mag(d0)*mag2(d)

      u = P + cross_k(Q)
      v = cross_k(P) - Q
      a_ij = (u.dot.tau_i)/R
      b_ij = (v.dot.tau_i)/R
    end subroutine calc_coeffs

    pure function V_iw_wake(r_i,r_ip,r_w)
      real(kind=rp),dimension(n_sd),intent(in) :: r_i, r_ip, r_w
      real(kind=rp),dimension(n_sd) :: V_iw_wake, d, s, s0
      real(kind=rp)                 :: z0, alpha, beta

      d     = r_ip - r_i
      s     = r_ip - r_w
      s0    = r_i  - r_w

      z0    = d.cross.s0
      alpha = atan((d.dot.s)/z0) - atan((d.dot.s0)/z0)
      beta  = 0.5_rp*log(mag2(s)/mag2(s0))

      V_iw_wake = -(alpha*d + beta*cross_k(d))/(2._rp*pi*mag2(d))
    end function V_iw_wake
  end subroutine enforce_bc

  pure subroutine purge_vor(ffield,num)
    class(flow_field),intent(inout) :: ffield
    integer          ,intent(in)    :: num

    ffield%has_vortex(                                            &
      largestloc(mag2(ffield%r,ffield%has_vortex),           &
                 num))                                            &
    = .FALSE.
    ffield%N = ffield%N - num
  end subroutine purge_vor

  pure subroutine import_vor(ffield,bndry)
    class(flow_field),intent(inout) :: ffield
    type(boundary)   ,intent(inout) :: bndry
    integer                         :: i, j

    j = 0
    do i = 1,ffield%Nmax
      if(.NOT.ffield%has_vortex(i)) then
        j = j + 1
        ffield%has_vortex(i) = .TRUE.
        ffield%not_marker(i) = .FALSE.
        ffield%Gama(i)       = bndry%gama(j)*bndry%ds(j)
        ffield%r(:,i)        = bndry%r_col(:,j) + bndry%delta*bndry%n(j)
        bndry%marker(j)      = i
      end if
      if(j == bndry%Ncn) exit
    end do
    ffield%N = ffield%N + bndry%Ncn
  end subroutine import_vor

  pure function same_orient(ffield,a,b)
    class(flow_field),intent(in) :: ffield
    integer          , intent(in) :: a, b
    logical :: same_orient

    same_orient =  (ffield%Gama(a) > 0 .AND. ffield%Gama(b) > 0)                    &
               .OR.(ffield%Gama(a) < 0 .AND. ffield%Gama(b) < 0)
  end function same_orient

  subroutine write_vor(ffield,flname,flnum)
    class(flow_field),intent(in) :: ffield
    character(len=*) ,intent(in) :: flname
    integer          ,intent(in),optional :: flnum
    character(len=7) :: stamp
    integer          :: i

    if(PRESENT(flnum)) then
      write(stamp,'(i3.3a4)') flnum,'.dat'
      open(unit=10,file=flname//stamp)
    else
      open(unit=10,file=flname//'.dat')
    end if
              
    do i = 1,ffield%Nmax
      if(ffield%has_vortex(i))                                           &
        write(10,*) ffield%r(:,i),ffield%Gama(i)
    end do
    close(10)
  end subroutine write_vor
end module flow
