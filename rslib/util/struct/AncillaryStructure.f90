MODULE AncillaryStructure
  !
  ! Module that defines all the structures related to ancillary I/O's
  !
  ! USE:
  !
  !    StateIndexModule
  !
  ! yhe@aer.com, 03/18/2016
  !
  USE StateIndexModule, Only: &
       StateIndex_t, &
       maxMol

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: &
       AncillaryGran_t, &
       AncillaryDefin_t, &
       AncillaryFOR_t, &
       AncillaryData_t

  TYPE AncillaryDefin_t
     integer :: nFOV
     integer :: nParG
     integer :: nEmMW
     integer :: nEmIR
     integer :: nEmRfIR
     integer, dimension(maxMol) :: MolID
     real, dimension(:), allocatable :: wvn  !size=nEmIR
     real, dimension(:), allocatable :: frq  !size=nEmMW
     real, dimension(:), allocatable :: pol  !size=nEmMW
     TYPE(StateIndex_t) :: IG
     TYPE(StateIndex_t) :: NG
  END TYPE AncillaryDefin_t

  TYPE AncillaryData_t
     real :: pSfc
     real :: elev
     real :: landFrac
     real, dimension(:), allocatable :: stateV !size=nParG
     real, dimension(:), allocatable :: emRfIR !size=nEmRfIR
  END TYPE AncillaryData_t

  TYPE AncillaryFOR_t
     TYPE(AncillaryData_t), dimension(:), allocatable :: data !size=nFOV
     integer :: nEmMW
     real :: elev
     real :: pSfc
     real :: landFrac
     real, dimension(:), allocatable :: emMW  !size=nEmMW
  END TYPE AncillaryFOR_t

  TYPE AncillaryGran_t
     integer :: nFOR
     TYPE(AncillaryDefin_t) :: def
     TYPE(AncillaryFOR_t), dimension(:), allocatable :: FOR  !size=nFOR
  END TYPE AncillaryGran_t

END MODULE AncillaryStructure
