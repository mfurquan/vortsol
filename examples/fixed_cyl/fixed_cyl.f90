program fixed_cyl
  use global_param
  use solid_boundary
  use flow
  implicit none

  integer           :: Ncn, Nvor, Nstep, istep, wfreq
  real(kind=rp)     :: h, Rey, dt, Uinfty(n_sd)
  type(boundary)    :: cyl
  type(flow_field)  :: cyl_flow
  namelist /inputs/ Ncn, Nvor, Nstep, wfreq, h, Rey, dt

  read(*,nml=inputs)
  Uinfty = [1._rp,0._rp]

  call cyl%set(RESHAPE([(circle(istep*2._rp*pi/Ncn),             &
               istep = 0,Ncn-1)],[n_sd,Ncn]),h)
  call cyl%shapeout('cylshape.dat')
  call cyl_flow%init(Nvor,h/2._rp,1.2_rp,0.02_rp)

  do istep = 1,Nstep
    call cyl_flow%step(Uinfty,Rey,dt,cyl)
    if(MOD(istep,wfreq)==0)                                      &
      call cyl_flow%write_vor('vortex',istep/wfreq)
!    call write_surfvel(istep)
  end do

contains
  pure function circle(theta)
    real(kind=rp),intent(in) :: theta
    real(kind=rp)            :: circle(n_sd)

    circle = 0.5_rp*[cos(theta),sin(theta)]
  end function circle

  subroutine write_surfvel(step)
    integer,intent(in) :: step
    character(len=7)   :: stamp
    real(kind=rp)      :: v(n_sd,cyl%N)
    integer            :: i

    v = cyl_flow%Vel(cyl%r_col(:,1:cyl%N),Uinfty,cyl)
    
    write(stamp,'(i3.3a4)') step,'.dat'
    open(10,file='vel_pol'//stamp)
    do i = 1,Ncn
      write(10,*) 2._rp*pi*(i-1)/Ncn, cyl% rej(v(:,i),i),        &
                                      cyl%proj(v(:,i),i)
    end do
    write(10,*) 2._rp*pi, cyl%rej(v(:,1),1),cyl%proj(v(:,1),1)
    close(10)
  end subroutine write_surfvel
end program fixed_cyl
