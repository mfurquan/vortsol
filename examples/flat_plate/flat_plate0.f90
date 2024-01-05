program flat_plate
  use global_param
  use solid_boundary
  use flow
  implicit none

  integer           :: Ncn, Nvor, Nstep, istep, wfreq
  real(kind=rp)     :: Rcore, delta, AR, Rey, dt, Uinfty(n_sd)
  type(boundary)    :: plate
  type(flow_field)  :: plate_flow
  namelist /inputs/ Ncn, Nvor, Nstep, wfreq, Rcore, delta, AR, Rey, dt

  read(*,nml=inputs)
  Uinfty = [1._rp,0._rp]

  call plate%set(RESHAPE([(plate_contr(istep*(1._rp/Ncn)), istep = 0,Ncn-1)],         &
                         [n_sd,Ncn]),delta,.FALSE.,.FALSE.)
  call plate%shapeout('plate.dat')
  call plate_flow%init(Uinfty,Nvor,Rcore,plate)

  do istep = 1,Nstep
    call plate_flow%step(Uinfty,Rey,dt,plate)
    if(MOD(istep,wfreq)==0)                                      &
      call plate_flow%write_vor('vortex',istep/wfreq)
    write(*,*) 'Step:',istep,'complete'
!   call write_surfvel(istep)
  end do
!   call write_surfvel(istep)

contains
  ! parameteric representation of plate contour
  pure function plate_contr(t)
    real(kind=rp),intent(in) :: t
    real(kind=rp)            :: plate_contr(n_sd), a, b, s, theta

    a = 1._rp/(2._rp*AR)                           ! cap radius
    b = 0.5_rp - a                                 ! half length of the shank

    s = trans(t,a,b)*(4._rp*b + 2._rp*pi*a)        ! t -> arc length (s)
    
    if(s <= pi*a) then
      theta = s/a - pi/2._rp
      plate_contr = [b + a*cos(theta),a*sin(theta)]
    else if(s <= 2._rp*b + pi*a) then
      plate_contr = [b + pi*a - s,a]
    else if(s <= 2._rp*(b + pi*a)) then
      theta = (s - 2._rp*b)/a - pi/2._rp
      plate_contr = [-b + a*cos(theta),a*sin(theta)]
    else
      plate_contr = [s - 3._rp*b - 2._rp*pi*a,-a]
    end if
  end function plate_contr

  ! transformation used for clustering nodes at the edges
  pure function trans(t,a,b)
    real(kind=rp),intent(in) :: t, a, b
    real(kind=rp)            :: trans, phi_2pi, u

    phi_2pi = pi*a/(8._rp*b + 4._rp*pi*a)
    u = 2._rp*pi*(t - phi_2pi)
    if(t < phi_2pi) then
      trans = cos(u) - cos(2._rp*pi*phi_2pi)
    else if(t < 0.5 + 3._rp*phi_2pi) then
      trans = 2._rp - cos(u) - cos(2._rp*pi*phi_2pi)
    else
      trans = 4._rp + cos(u) - cos(2._rp*pi*phi_2pi)
    end if

    trans = trans/4._rp
  end function trans

  subroutine write_surfvel(step)
    integer,intent(in) :: step
    character(len=7)   :: stamp
    real(kind=rp)      :: v(n_sd,plate%Ncn)
    integer            :: i

    v = plate_flow%Vel(plate%r_col(:,1:plate%Ncn),Uinfty,plate)
    
    write(stamp,'(i3.3a4)') step,'.dat'
    open(10,file='vel_pol'//stamp)
    do i = 1,Ncn
      write(10,*) 2._rp*pi*(i-1)/Ncn, plate%n_dot(v(:,i),i),plate%tau_dot(v(:,i),i)
    end do
    write(10,*) 2._rp*pi, plate%n_dot(v(:,1),1),plate%tau_dot(v(:,1),1)
    close(10)
  end subroutine write_surfvel
end program flat_plate
