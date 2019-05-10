module params_module

  use type_kinds, ONLY: Double, FP
  implicit none

  public

 !------------------------------------------------------------------------
 !    Profile Parameters
 !------------------------------------------------------------------------
  INTEGER, PARAMETER  :: mxLev=120
  integer, parameter  :: mxlay=mxlev-1
  
 !------------------------------------------------------------------------
 !    Multiple Scattering Parameters
 !------------------------------------------------------------------------
  integer, parameter :: mxcang=3         ! For nstr=4, this must be >= 2
  integer, parameter :: maxcmu=2*mxcang  ! must be >= nstr in ossdrv_ir
  integer, parameter :: mxang=mxcang+1
  integer, parameter :: mxang2=mxang*mxang
  real(fp), save :: accur
  
  data accur/1.e-04/

  real(fp), parameter :: one_fp=1., two_fp=2., ten_fp=10.
  real(fp), parameter :: small_fp=ten_fp*tiny(one_fp)
  real(Double), parameter :: one_dp=1., two_dp=2., ten_dp=10.
  real(Double), parameter :: small_dp=ten_dp*tiny(one_dp)

 !------------------------------------------------------------------------
 ! used in oss_ir_module and main
 !------------------------------------------------------------------------
  INTEGER, PARAMETER :: mxHmol=25
  integer, public, parameter :: MxmolS=20
  INTEGER, PARAMETER :: mxParG=(mxHmol+1)*mxlev
  

end module params_module
