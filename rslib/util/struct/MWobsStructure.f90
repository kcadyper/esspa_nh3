MODULE MWobsStructure
  !
  ! Module that defines all the microwave (MW) observation I/O's
  !
  ! yhe@aer.com, 03/18/2016
  !
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: &
       MWgran_t, &
       MWdefin_t, &
       MWob_t

  TYPE MWdefin_t
     integer :: nChan
     integer :: radOrTb
     integer, dimension(6) :: startTime ! year,month,day,hour,minute,second
     integer, dimension(6) :: endTime
     real, dimension(:), allocatable :: frq  !size=nChan
     integer, dimension(:), allocatable :: pol  !size=nChan
     integer, dimension(:), allocatable :: chanID  !size=nChan
     real, dimension(:), allocatable :: planckAlpha  !size=nChan
     real, dimension(:), allocatable :: planckBeta   !size=nChan
  END TYPE MWdefin_t

  TYPE MWob_t
     real, dimension(:), allocatable :: rad     !size=nChan
     integer, dimension(:), allocatable :: QC   !size=nChan
     real :: lat
     real :: lon
     real :: EIA
     real :: EAA
     real :: scanAng
     integer :: scanPos
     real :: surf_alt
     real :: land_frac
     integer, dimension(8) :: obs_time_utc
  END TYPE MWob_t

  TYPE MWgran_t
     integer :: nFOV
     TYPE(MWdefin_t) :: def
     TYPE(MWob_t), dimension(:), allocatable :: ob  !size=nFOV
  END TYPE MWgran_t

END MODULE MWobsStructure
