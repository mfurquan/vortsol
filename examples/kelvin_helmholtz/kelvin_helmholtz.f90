program kelvin_helmholtz
  use global_param
  use flow
  implicit none

  integer           :: Nvor, Nstep, istep, wfreq
  real(kind=rp)     :: R_rnk, Rey, dt, R, del
  type(flow_field)  :: sheet
  namelist /inputs/ Nvor, R, Nstep, wfreq, R_rnk, Rey, dt

  read(*,nml=inputs)

  ! initialize oseen-lamb vortex
  call sheet%init(Nvor,R_rnk,1.2_rp,0.02_rp)
  sheet%N          = Nvor
  sheet%has_vortex = .TRUE.
  sheet%Gama       = 1._rp/Nvor

  ! make a square
  del              = 1._rp/(Nvor/4)
  sheet%r(1,1:Nvor/4)          = [(del*istep, istep = 0,Nvor/4-1)]
  sheet%r(2,1:Nvor/4)          = 0._rp
  sheet%r(1,Nvor/4+1:Nvor/2)   = 1._rp
  sheet%r(2,Nvor/4+1:Nvor/2)   = [(del*istep, istep = 0,Nvor/4-1)]
  sheet%r(1,Nvor/2+1:3*Nvor/4) = [(del*istep, istep = 1,Nvor/4)]
  sheet%r(2,Nvor/2+1:3*Nvor/4) = 1._rp
  sheet%r(1,3*Nvor/4+1:Nvor)   = 0._rp
  sheet%r(2,3*Nvor/4+1:Nvor)   = [(del*istep, istep = 1,Nvor/4)]

!  ! make a circle
!  sheet%r(1,:) = [(R*cos(2._rp*pi*istep/Nvor), istep = 0,Nvor-1)]
!  sheet%r(2,:) = [(R*sin(2._rp*pi*istep/Nvor), istep = 0,Nvor-1)]

  call sheet%write_vor('vortex',0)

  do istep = 1,Nstep
    call sheet%step([0._rp,0._rp],Rey,dt)
    if(MOD(istep,wfreq)==0) then
      call sheet%write_vor('vortex',istep/wfreq)
    end if
  end do
end program kelvin_helmholtz
