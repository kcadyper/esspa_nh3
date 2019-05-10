Module OutputStructure
  USE StateIndexModule, Only: &
       StateIndex_t

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: &
       OutputFOR_t, &
       OutputDefin_t, &
       OutputGran_t, &
       LinearFOR_t, &
       LinearDefin_t, &
       LinearGran_t

  !OutputDefin_t, OutputFOR_t, OutputGran_t were introduced on
  !06/08/2016, to replace RetrOutput_t
  TYPE OutputDefin_t
     integer :: nFOV
     integer :: nParG
     integer :: nEmMW
     integer :: nEmRfIR
     integer :: nLev
     integer :: roUnits
     integer, dimension(6) :: startTime !granule start time
     TYPE(StateIndex_t) :: IG
     TYPE(StateIndex_t) :: NG
  END TYPE OutputDefin_t

  TYPE OutputFOR_t
     integer :: nIter
     integer :: landType
     real :: rms
     real :: chisq
     real :: DOFS
     real :: lat
     real :: lon
     real :: EIA
     real :: SolIA
     real :: landFrac
     real, dimension(:), allocatable :: state !size=(nParG)
     real, dimension(:), allocatable :: press !size=(nLev)
     real, dimension(:), allocatable :: emMW  !size=(nEmMW)
     real, dimension(:), allocatable :: emisIR !size=(nEmRfIR)
     real, dimension(:), allocatable :: DiagError   !size=(nParG)
     integer :: atmosClass
  END TYPE OutputFOR_t

  TYPE OutputGran_t
     integer :: nFOR
     TYPE(OutputDefin_t) :: def
     TYPE(OutputFOR_t), dimension(:), allocatable :: FOR !size=(nFOR)
     TYPE(OutputFOR_t), dimension(:,:), allocatable :: FOR2D!size=(nFOV,nFOR)
  END TYPE OutputGran_t

  TYPE LinearFOR_t
     real, dimension(:), allocatable :: xOffG
     real, dimension(:,:), allocatable :: xGainG
     real :: Tsfc
     real :: chisq
  END TYPE LinearFOR_t

  TYPE LinearDefin_t
     integer :: nParG
     integer :: nChan
     integer, dimension(6) :: startTime !granule start time
     integer, dimension(:), allocatable :: kchan  !size=nChanMW+nChanIR
     TYPE(StateIndex_t) :: IG
     TYPE(StateIndex_t) :: NG
  END TYPE LinearDefin_t

  TYPE LinearGran_t
     integer :: nFOR
     TYPE(LinearDefin_t) :: def
     TYPE(LinearFOR_t), dimension(:), allocatable :: FOR
  END TYPE LinearGran_t

END Module OutputStructure
