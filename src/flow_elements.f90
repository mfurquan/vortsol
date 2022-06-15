module flow_elements
  use global_param
  use utility

contains
  ! velocity induced by a Rankine vortex of unit strength
  pure function rkn_vortex(r_p,r_v,R2_rkn)
    real(kind=rp),intent(in) :: r_p(n_sd), r_v(n_sd), R2_rkn
    real(kind=rp)            :: rkn_vortex(n_sd), dr2

    rkn_vortex = rot90(r_p - r_v)
    dr2        = mag2(rkn_vortex)

    rkn_vortex = rkn_vortex/(2._rp*pi*MAX(dr2,R2_rkn))
  end function rkn_vortex

  ! velocity induced by a source sheet of unit strength
  pure function src_sheet(r_p,r_s,ds,e_s)
    real(kind=rp),intent(in) :: r_p(n_sd), r_s(n_sd), ds, e_s(n_sd)
    real(kind=rp)            :: src_sheet(n_sd), dr(n_sd), xi,  &
                                eta, xip, xim, Vr, Vth

    dr  = r_p - r_s
    xi  = DOT_PRODUCT(dr,e_s)
    eta = DOT_PRODUCT(dr,rot90(e_s))

    xip = xi + ds/2._rp
    xim = xi - ds/2._rp

    Vr  = 0.5_rp*log((xim**2+eta**2)/(xip**2+eta**2))
    !Vth = atan2(xim,eta) - atan2(xip,eta)
    Vth = atan(eta/xim) - atan(eta/xip)

    src_sheet(1) = Vr*e_s(1) - Vth*e_s(2)
    src_sheet(2) = Vr*e_s(2) + Vth*e_s(1)
  end function src_sheet
end module flow_elements

