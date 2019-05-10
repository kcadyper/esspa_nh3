MODULE MeasErrStructure
  !
  ! Module that defines measurement error structure for MW and IR
  !
  ! USE:
  !
  !    None
  !
  ! yhe@aer.com, 08/09/2016
  !
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: &
       MeasErr_t

  TYPE MeasErr_t
     integer :: nChan
     real, dimension(:), allocatable :: device    !size=nChan
     real, dimension(:), allocatable :: fwdModel  !size=nChan
     real :: rms
  END TYPE MeasErr_t

END MODULE MeasErrStructure
