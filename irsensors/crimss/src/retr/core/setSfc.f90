MODULE SurfaceTypingModule
  USE StateIndexModule
  IMPLICIT NONE
  PRIVATE
  !------------------------------------------------------
  !	Public items made available to calling program(s)   
  !------------------------------------------------------
  PUBLIC :: setSfc
  !---------------------------------------------------------
  !     Declarations of private data (private to the module)
  !---------------------------------------------------------
  integer, parameter :: mxLength=256
  logical :: dbg = .false.
CONTAINS

  subroutine setSfc(igeo,plandAvg,plandMaxOc)
    !--I/O variables
    real,               intent(in)    :: plandAvg
    real,               intent(in)    :: plandMaxOc
    integer,            intent(inout) :: igeo
    !Local
    character (len=mxLength) :: procName

    procName = '[SurfaceTypingModule::setSfc]:'
    if (dbg) print *, trim(procName)//' starting ...'
    if(plandAvg <= plandMaxOc)then
       igeo=1
    else
       igeo=2
    end if
    if (dbg) print *, trim(procName)//' ending ...'
    return
  end subroutine setSfc
  
end MODULE SurfaceTypingModule

