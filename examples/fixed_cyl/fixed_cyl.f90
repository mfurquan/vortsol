program fixed_cyl
  use vectors
  use contours
  use vortex
  implicit none

  integer           :: Ncn, Nvor, Nstep, istep
  real(kind=rp)     :: h, Rey, dt
  type(vector)      :: Uinfty
  type(contour)     :: cyl
  type(vortexList)  :: vor
  namelist /inputs/ Ncn, Nvor, Nstep, h, Rey, dt

  read(*,nml=inputs)
  Uinfty = [1._rp,0._rp]

  call cyl%set(circle,Ncn,2._rp*pi,h)
  call cyl%shapeout('cylshape.dat')
  call vor%init(Nvor,h/2._rp)
  call vor%gen_shear(cyl,Uinfty)
  call vor%write_vor('vortex',0)
  call write_surfvel(0)

  do istep = 1,Nstep
    call vor%move(cyl,Uinfty,Rey,dt)
    call vor%gen_shear(cyl,Uinfty)
    call vor%write_vor('vortex',istep)
    call write_surfvel(istep)
  end do

contains
  pure elemental function circle(theta)
    real(kind=rp),intent(in) :: theta
    type(vector) :: circle

    circle = [cos(theta),sin(theta)]
  end function circle

  subroutine write_surfvel(step)
    integer,intent(in) :: step
    character(len=7)   :: stamp
    type(vector)       :: v(cyl%n_pts)
    integer            :: i

    v = vor%Vel(cyl%rc(1:cyl%n_pts),Uinfty)
    
    write(stamp,'(i3.3a4)') step,'.dat'
    open(10,file='vel_pol'//stamp)
    do i = 1,Ncn
      write(10,*) (i-1)/(2._rp*pi), v(i).dot.cyl%n(i),           &
                                    v(i).dot.cyl%t(i)
    end do
    close(10)
  end subroutine write_surfvel
end program fixed_cyl
