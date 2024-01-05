program flat_plate
  use global_param
  use solid_boundary
  use flow
  implicit none

  integer           :: Ncn, Nvor, Nstep, istep, wfreq
  real(kind=rp)     :: Rcore, delta, AR, Rey, dt, Uinfty(n_sd), a, b, Lc
  type(boundary)    :: plate
  type(flow_field)  :: plate_flow
  namelist /inputs/ Ncn, Nvor, Nstep, wfreq, Rcore, delta, AR, Rey, dt

  read(*,nml=inputs)
  Uinfty = [1._rp,0._rp]

  a = 1._rp/(2._rp*AR)                           ! cap radius
  b = 0.5_rp - a                                 ! half length of the shank
  Lc = 4._rp*b + 2._rp*pi*a                      ! contour length
  call plate%set(RESHAPE([(plate_contr(trans(istep*(Lc/Ncn),Lc),a,b), istep = 0,Ncn-1)],    &
                         [n_sd,Ncn]),delta,.FALSE.,.FALSE.)
  call plate%shapeout('plate.dat')
  call plate_flow%init(Uinfty,Nvor,Rcore,plate)
  !call plate_flow%write_vor('vortex',0)

  do istep = 1,Nstep
    call plate_flow%step(Uinfty,Rey,dt,plate)
    if(MOD(istep,wfreq)==0)  call plate_flow%write_vor('vortex',istep/wfreq)
    write(*,*) 'Step:',istep,'complete'
!   call write_surfvel(istep)
  end do
!   call write_surfvel(istep)

contains
  ! parameteric representation of plate contour
  pure function plate_contr(s,a,b)
    real(kind=rp),intent(in) :: s, a, b
    real(kind=rp)            :: plate_contr(n_sd), theta

    if(s <= pi*a/2._rp) then
      theta = s/a
      plate_contr = [b + a*cos(theta),a*sin(theta)]
    else if(s <= 2._rp*b + pi*a/2._rp) then
      plate_contr = [b + pi*a/2._rp - s,a]
    else if(s <= 2._rp*b + 3._rp*pi*a/2._rp) then
      theta = (s - 2._rp*b)/a
      plate_contr = [-b + a*cos(theta),a*sin(theta)]
    else if(s <= 4._rp*b + 3._rp*pi*a/2._rp) then
      plate_contr = [s - 3._rp*b - 3._rp*pi*a/2._rp,-a]
    else
      theta = (s - 4._rp*b)/a
      plate_contr = [b + a*cos(theta),a*sin(theta)]
    end if
  end function plate_contr

  pure function trans(x,L)
    real(kind=rp),intent(in) :: x, L
    real(kind=rp)            :: trans, k

    ! smooth step function used for clustering nodes at the edges
    k = 4._rp*pi/L
    !trans = (k*x - sin(sin(sin(sin(sin(k*x))))))/k ! more sin => more uniform distribution

     ! identity function for uniform distribution of nodes
     trans = x
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
