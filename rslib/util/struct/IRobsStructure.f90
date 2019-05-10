MODULE IRobsStructure
  !
  ! Module that defines all the infrared observation I/O's
  !
  ! yhe@aer.com, 03/18/2016
  !
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: &
       IRdefin_t, &
       IRob_t,    &
       IRgran_t,  &
       IRFOR_t

  TYPE IRdefin_t
     integer :: nFOV
     integer :: nChan
     integer, dimension(6) :: startTime ! year,month,day,hour,minute,second
     integer, dimension(6) :: endTime
     real, dimension(:,:), allocatable :: NEdN !size=nChan
     real, dimension(:), allocatable :: wvn  !size=nChan
     character (len=12), dimension(:), allocatable :: chanID  !size=nChan
     integer, dimension(:), allocatable :: chanNum !size=nChan
  END TYPE IRdefin_t

  TYPE IRob_t
     real, dimension(:), allocatable :: rad  !size=nChan
     integer, dimension(:), allocatable :: QC   !size=nChan
     real :: lat
     real :: lon
     real :: EIA
     real :: EAA
     real :: SolIA
     real :: SolAA
     real :: scanAng
     real :: surf_alt
     real :: land_frac
  END TYPE IRob_t

  TYPE IRFOR_t
     TYPE(IRob_t), dimension(:), allocatable :: ob  !size=nFOV
     real :: lat
     real :: lon
     real :: EIA
     real :: EAA
     real :: SolIA
     real :: SolAA
     real :: scanAng
     integer, dimension(8) :: obs_time_utc
  END TYPE IRFOR_t

  TYPE IRgran_t
     integer :: nFOR
     TYPE(IRdefin_t) :: def
     TYPE(IRFOR_t), dimension(:), allocatable :: FOR  !size=nFOR
  END TYPE IRgran_t

END MODULE IRobsStructure
