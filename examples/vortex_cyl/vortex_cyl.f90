program fixed_cyl
  use global_param
  use solid_boundary
  use flow
  implicit none

  integer           :: Ncn, Nvor, istep
  real(kind=rp)     :: R, R_rnk, d, Gama, theta, Uinfty(n_sd)
  type(boundary)    :: cyl
  type(flow_field)  :: cyl_flow
  namelist /inputs/ Ncn, Nvor, R, R_rnk, d, Gama, theta

  read(*,nml=inputs)

  Uinfty = [0._rp,0._rp]
  call cyl%set(RESHAPE([(circle(2._rp*pi*istep/Ncn),             &
               istep = 0,Ncn-1)],[n_sd,Ncn]),0.01_rp,.FALSE.,.FALSE.)
  call cyl%shapeout('cylshape.dat')
  call cyl_flow%init(Uinfty,Nvor,R_rnk,cyl)
  cyl_flow%N = 1
  cyl_flow%has_vortex(1) = .TRUE.
  cyl_flow%Gama(1) = Gama
  theta = theta*pi/Ncn
  cyl_flow%r(:,1) = (R+d)*[cos(theta),sin(theta)]
  call cyl_flow%enforce_bc(cyl,Uinfty)
  call write_gama
  call cyl_flow%decom
  call cyl%decom

contains
  pure function circle(th)
    real(kind=rp),intent(in) :: th
    real(kind=rp)            :: circle(n_sd)

    circle = R*[cos(th),sin(th)]
  end function circle

  subroutine write_gama
    integer :: i
    real(kind=rp) :: phi
    open(10,file='gama.dat')
    do i = 1,Ncn
      phi = (i-1)*2._rp*pi/Ncn
      write(10,*) phi,cyl%gama(i),gama_ex(phi)
    end do
  end subroutine write_gama

  pure function gama_ex(ang)
    real(kind=rp),intent(in) :: ang
    real(kind=rp) :: gama_ex

    gama_ex = (Gama/pi)*((R+d)*cos(ang-theta)-R)                         &
            / (2._rp*R*(R+d)*(cos(ang-theta)-1._rp)-d*d)
  end function gama_ex
end program fixed_cyl
