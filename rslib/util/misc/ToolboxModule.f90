MODULE ToolboxModule
  !
  ! Collection of utility tools
  !
  ! AER, Inc., 04/29/2016
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: &
       getUnit

  integer, parameter :: mxLength = 128

CONTAINS
  integer function getUnit(startLUN,endLUN)
    integer, intent(in), optional :: startLUN
    integer, intent(in), optional :: endLUN
    integer :: unit
    integer :: sLUN=11, eLUN=99
    logical :: opened

    opened = .true.
    if (present(startLUN)) sLUN = startLUN
    if (present(endLUN))   eLUN = endLUN

    do unit = sLUN, eLUN
       INQUIRE(UNIT=unit,OPENED=opened)
       if (.not. opened) then
          getUnit = unit
          exit
       end if
    end do
    if (opened) then
       print *, '[ToolBoxModule::getUnit]: warning - no valid logic unit found, returning -999'
       getUnit = -999
    end if
  end function getUnit

END MODULE ToolboxModule
