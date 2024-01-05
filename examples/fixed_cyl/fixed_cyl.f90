program fixed_cyl
  use global_param
  use solid_boundary
  use flow
  implicit none

  integer           :: Ncn, Nvor, Nstep, istep, wfreq
  real(kind=rp)     :: Rcore, delta, Rey, dt, Uinfty(n_sd)
  type(boundary)    :: cyl
  type(flow_field)  :: cyl_flow
  namelist /inputs/ Ncn, Nvor, Nstep, wfreq, Rcore, delta, Rey, dt

  read(*,nml=inputs)
  Uinfty = [1._rp,0._rp]

  call cyl%set(RESHAPE([(circle(2._rp*pi*istep/Ncn),             &
               istep = 0,Ncn-1)],[n_sd,Ncn]),delta,.FALSE.,.FALSE.)
  call cyl%shapeout('cylshape.dat')
  call cyl_flow%init(Uinfty,Nvor,Rcore,cyl)
  call cyl_flow%write_vor('vortex',0)

  do istep = 1,Nstep
    call cyl_flow%step(Uinfty,Rey,dt,cyl)
    if(MOD(istep,wfreq)==0)                                     &
      call cyl_flow%write_vor('vortex',istep/wfreq)
    write(*,*) 'Step:',istep,'complete'
    call write_Vn(istep)
!    call write_surfvel(istep)
  end do
!  call write_surfvel(istep)

  call cyl_flow%decom
  call cyl%decom

contains
  pure function circle(theta)
    real(kind=rp),intent(in) :: theta
    real(kind=rp)            :: circle(n_sd)

    circle = 0.5_rp*[cos(theta),sin(theta)]
  end function circle

  subroutine write_surfvel(step)
    integer,intent(in) :: step
    character(len=7)   :: stamp
    real(kind=rp)      :: v(n_sd,cyl%Ncn)
    integer            :: i

    v = cyl_flow%Vel(cyl%r_col(:,1:cyl%Ncn),Uinfty,cyl)
    
    write(stamp,'(i3.3a4)') step-1,'.dat'
    open(10,file='vel_pol'//stamp)
    do i = 1,Ncn
      write(10,*) 2._rp*(i-1)/Ncn, cyl%n_dot  (v(:,i),i),           &
                                   cyl%tau_dot(v(:,i),i)
    end do
    write(10,*) 2._rp, cyl%n_dot(v(:,1),1),cyl%tau_dot(v(:,1),1)
    close(10)
  end subroutine write_surfvel

  subroutine write_Vn(step)
    integer,intent(in) :: step
    character(len=7)   :: stamp
    integer            :: i

    
    write(stamp,'(i3.3a4)') step,'.dat'
    open(10,file='vel_n'//stamp)
    do i = 1,cyl_flow%Nmax
      if(cyl_flow%has_vortex(i)) then
              write(10,'(4F16.9)') cyl_flow%r(:,i), (cyl_flow%V_vor(:,i).dot.cyl_flow%r(:,i))/NORM2(cyl_flow%r(:,i))
      end if
    end do
    close(10)
  end subroutine write_Vn
end program fixed_cyl
