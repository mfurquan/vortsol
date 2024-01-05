program oseen_lamb
  use global_param
  use flow
  implicit none

  integer           :: Nvor, Nstep, istep, wfreq, Nr
  real(kind=rp)     :: R_rnk, Rey, dt, delta, Rmax
  type(flow_field)  :: vortex
  real(kind=rp),allocatable :: rth(:,:), r(:)
  namelist /inputs/ Nvor, delta, Nstep, wfreq, R_rnk, Rey, dt,  &
                    Nr, Rmax

  read(*,nml=inputs)

  ! initialize oseen-lamb vortex
  call vortex%init([0._rp,0._rp],Nvor,R_rnk)
  vortex%N          = Nvor
  vortex%has_vortex = .TRUE.
  vortex%Gama       = 1._rp/Nvor
  ALLOCATE(rth(n_sd,Nvor))
  call RANDOM_NUMBER(rth)
  vortex%r(1,:) = delta*rth(1,:)*cos(2._rp*pi*rth(2,:))
  vortex%r(2,:) = delta*rth(1,:)*sin(2._rp*pi*rth(2,:))
  DEALLOCATE(rth)

  r = [(istep*Rmax/(Nr-1), istep = 0,Nr-1)]
  call vortex%write_vor('vortex',0)
  call write_circulation(0)

  do istep = 1,Nstep
    call vortex%step([0._rp,0._rp],Rey,dt)
    if(MOD(istep,wfreq)==0) then
      call vortex%write_vor('vortex',istep/wfreq)
      call write_circulation(istep)
    end if
  end do

contains
  subroutine write_circulation(step)
    integer,intent(in) :: step
    character(len=7)   :: stamp
    real(kind=rp)      :: G(Nr), gama
    integer            :: i

    gama    = vortex%Gama(1)
    G(1) = 0._rp
    do concurrent(i = 2:Nr)
      G(i) = COUNT(NORM2(vortex%r,1) < r(i))*gama
    end do

    write(stamp,'(i3.3a4)') step,'.dat'
    open(10,file='circulation'//stamp)
    do i = 1,Nr
      write(10,*) r(i),G(i)
    end do
    close(10)
  end subroutine write_circulation
end program oseen_lamb
