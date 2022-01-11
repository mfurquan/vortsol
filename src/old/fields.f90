module fields
  use global_param
  implicit none

  type :: flowField
    type(vector)     :: Uinfty
    type(contour)    :: SF_interface
    type(vortexList) :: vortex

  contains
    procedure :: init
  end type flowField

contains
  pure subroutine init(this,U0,surf_nods,h_shear,pc_buff,n_vor)
    class(flowField),intent(inout) :: this
    type(vector)    ,intent(in)    :: U0
    real(kind=rp)   ,intent(in)    :: surf_nods(:,:), h_shear
    integer         ,intent(in)    :: pc_buff, n_vor

    this%Uinfty = U0
    call this%SF_interface%init(surf_nods,h_shear,0)
    call this%vortex%init(n_vor,h_shear/2,1.2_rp)
  end subroutine init
end module flowField_m
