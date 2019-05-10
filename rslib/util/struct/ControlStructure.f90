MODULE ControlStructure
  !
  ! Module that defines the general control structure
  !
  ! Derived data type:
  !
  !     GenControl_t
  !
  ! USE:
  !
  !     StateIndexModule
  !
  ! Notes: Many modules, e.g., MWRTmodule and IRRTmodule, rely on the
  ! general control structure to be made available before they can be
  ! comppiled. This dependency criterion demands that the GenControl_t
  ! be located somewhere up in the compiling chain. The structure,
  ! StateSpaceDefin_t, is kept together with GenControl_t for now,
  ! since it might fall into the same situation, but can be relocated,
  ! if necessary, to a different directory.
  !
  ! yhe@aer.com, 04/05/2016
  !
  USE StateIndexModule, Only: &
       StateIndex_t, &
       StateFlags_t, &
       getNmol, &
       maxMol, &
       stateIndicesMatch

  USE VertCoord, ONLY: &
       mxCoordTyp, &
       Pcoord

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: &
       GenControl_t, &
       StateSpaceDefin_t, &
       RetrResult_t, &
       LinearResult_t, &
       State_t, &
       Jacobian_t

  INTEGER, PARAMETER :: mxLength=256
  real, parameter :: deltaP=0.05
  real, parameter :: deltaFrq=0.001
  logical :: dbg = .TRUE.

  PUBLIC :: stateSpacesMatch

  TYPE GenControl_t
     integer :: nFORmax
     integer :: apodizationType
     integer, dimension(:), allocatable :: kchan !mxchan - temporary
     real    :: LWP1stGuess
     real    :: IWP1stGuess
     logical :: enableScatter
     logical :: statRetrOn
     logical :: MWretrOn
     logical :: primaryRetrOn
     logical :: secondaryRetrOn
     logical :: MWdataOn
     logical :: IRdataOn
     logical :: imgDataOn
     logical :: logSize
     logical :: independentFOV
     character(len=mxLength) :: atmosClassConfig
     character(len=mxLength) :: MWsurfaceClassConfig
     character(len=mxLength) :: IRsurfaceClassConfig
     character(len=mxLength) :: limitsConfig
     character(len=mxLength) :: backgroundConfig
     character(len=mxLength) :: MWnoiseFile
     character(len=mxLength) :: MWnoiseAvgFactorFile
     character(len=mxLength) :: IRnoiseFile
     character(len=mxLength) :: MWobsFile
     character(len=mxLength) :: IRobsFile
     character(len=mxLength) :: imgObsFile
     character(len=mxLength) :: ancillaryAtmosFile
     character(len=mxLength) :: ancillarySfcFile
     character(len=mxLength) :: ancillaryMWemisFile
     character(len=mxLength) :: ancillaryIRemisFile
     character(len=mxLength) :: statRetrOutFile
     character(len=mxLength) :: MWretrOutFile
     character(len=mxLength) :: primaryRetrOutFile
     character(len=mxLength) :: secondaryRetrOutFile
     character(len=mxLength) :: backgroundOutFile
     character(len=mxLength) :: linearInvertOutFile
     character(len=mxLength), dimension(:), allocatable :: backgroundFiles !size=maxBkg
     character(len=mxLength), dimension(:), allocatable :: channelSelectFiles !size=maxChSet
     TYPE(StateFlags_t) :: externDataFlags
  END TYPE GenControl_t

  TYPE StateSpaceDefin_t
     integer :: nParG
     integer :: nParGatm
     integer, dimension(maxMol) :: MolID
     integer, dimension(maxMol) :: MolTran
     logical :: logSize
     TYPE(StateIndex_t) :: IG
     TYPE(StateIndex_t) :: NG
     integer :: nLev
     character(LEN=mxCoordTyp) :: vCoordTyp
     real, dimension(:), allocatable :: pRef   !size=nLev
     integer :: nChanMW
     real, dimension(:), allocatable :: frqMW  !size=nChanMW
     integer, dimension(:), allocatable :: polMW  !size=nChanMW
     integer :: nEmIR
     real, dimension(:), allocatable :: wvnIR  !size=nEmIR
  END TYPE StateSpaceDefin_t

  !Relocated from irsensors/generic/src/retr/RetrModule.f90
  TYPE RetrResult_t
     integer :: nParG
     integer :: nParGmax
     TYPE(StateIndex_t) :: IG
     TYPE(StateIndex_t) :: NG
     integer :: nEmMW
     integer :: nEmMWmax
     integer :: nEmRfIR
     integer :: nEmRfIRmax
     real, dimension(:), allocatable :: press  !size=nLev
     real, dimension(:), allocatable :: stateV !size=nParGmax
     real, dimension(:), allocatable :: diagError !size=nParGmax
     real, dimension(:), allocatable :: emMW   !size=nEmMWmax
     real, dimension(:), allocatable :: emRfIR !size=nEmRfIRmax
     integer :: nIter
     real :: rms
     real :: chiSq
     real :: dofs
     integer :: roUnits
  END TYPE RetrResult_t

  !Relocated from irsensors/generic/src/retr/RetrModule.f90
  TYPE LinearResult_t
     integer :: nParG
     integer :: nChan
     real :: chiSq
     real :: Tsfc
     real, dimension(:), allocatable :: xOffG  !size=nParG
     real, dimension(:,:), allocatable :: xGainG  !size=(nParG,nChan)
  END TYPE LinearResult_t

  !temporary location
  TYPE State_t
     real, dimension(:), allocatable :: temp
     real :: Tskin
     real, dimension(:,:), allocatable :: mol
     real :: pSfc
     real, dimension(:), allocatable :: wind
     real, dimension(:), allocatable :: cloudLiq
     real, dimension(:), allocatable :: cloudIce
     !TYPE(CloudScene_t) :: cldScene !type CloudScene (in cloud handling)
     real, dimension(:), allocatable :: emMW
     real, dimension(:), allocatable :: emIR
     real, dimension(:), allocatable :: rfIR
  END TYPE State_t

  !temporary location
  TYPE Jacobian_t
     real, dimension(:,:),   allocatable :: temp !(:,nChan)
     real, dimension(:),     allocatable :: Tskin !(nChan)
     real, dimension(:,:,:), allocatable :: mol !(:,:,nChan)
     real, dimension(:),     allocatable :: pSfc !(nChan)
     real, dimension(:,:),   allocatable :: wind !(:,nChan)
     real, dimension(:,:),   allocatable :: cloudLiq !(:,nChan)
     real, dimension(:,:),   allocatable :: cloudIce !(:,nChan)
     real, dimension(:),     allocatable :: cldScene !(nChan)
     real, dimension(:,:),   allocatable :: emMW !(:,nChan)
     real, dimension(:,:),   allocatable :: emIR !(:,nChan)
     real, dimension(:,:),   allocatable :: rfIR !(:,nChan)
  END TYPE Jacobian_t

CONTAINS
  logical function stateSpacesMatch(stateSpaceA,stateSpaceB,atmos,sfcMW,sfcIR)
    TYPE(StateSpaceDefin_t), intent(in) :: stateSpaceA
    TYPE(StateSpaceDefin_t), intent(in) :: stateSpaceB
    logical, intent(in), optional :: atmos
    logical, intent(in), optional :: sfcMW
    logical, intent(in), optional :: sfcIR
    !local
    integer :: nMolA
    integer :: nMolB
    logical :: NGstatus
    logical :: IGstatus
    character (len=mxLength) :: procName

    procName = ' [ControlStructure::stateSpacesMatch]: '
    if (dbg) print *, trim(procName)//'starting ...'

    stateSpacesMatch = .TRUE.
    if (present(atmos)) then
       if (atmos) then
          if (stateSpaceA%nParGatm /= stateSpaceB%nParGatm) then
             print*,trim(procName)//'error - nParGatm inconsistent', &
                  stateSpaceA%nParGatm,stateSpaceB%nParGatm
             stateSpacesMatch = .FALSE.
          endif
          if (trim(stateSpaceA%vCoordTyp) /= trim(stateSpaceB%vCoordTyp)) then
             print*,trim(procName)//'error - vCoordTyp inconsistent', &
                  trim(stateSpaceA%vCoordTyp),trim(stateSpaceB%vCoordTyp)
             stateSpacesMatch = .FALSE.
          endif
          if (stateSpaceA%nLev /= stateSpaceB%nLev) then
             print*,trim(procName)//'error - nLev inconsistent', &
                  stateSpaceA%nLev,stateSpaceB%nLev
             stateSpacesMatch = .FALSE.
          endif
          if (trim(stateSpaceA%vCoordTyp) == Pcoord .and. &
             trim(stateSpaceB%vCoordTyp) == Pcoord) then
             IF (ANY(abs(stateSpaceA%pRef(1:stateSpaceA%nLev)- &
                  stateSpaceB%pRef(1:stateSpaceB%nLev)) > deltaP )) then
                print*,trim(procName)//'error - pRef inconsistent'
                stateSpacesMatch = .FALSE.
             endif
          else if (trim(stateSpaceA%vCoordTyp) /= trim(stateSpaceB%vCoordTyp)) then
             print *, trim(procName)//'error - inconsistent V-coordinates'
             call errorHalt(1)
          endif
          nMolA = getNmol(stateSpaceA%MolID)
          nMolB = getNmol(stateSpaceB%MolID)
          if (nMolA /= nMolB) then
             print*,trim(procName)//'error - nMol inconsistent', nMolA, nMolB
             stateSpacesMatch = .FALSE.
          endif
          if (ANY((stateSpaceA%MolID(1:nMolA)-stateSpaceB%MolID(1:nMolB)) /= 0)) then
             print*,trim(procName)//'error - MolID inconsistent'
             stateSpacesMatch = .FALSE.
          endif
          if (ANY((stateSpaceA%MolTran(1:nMolA)-stateSpaceB%MolTran(1:nMolB)) /= 0)) then
             print*,trim(procName)//'error - MolTran inconsistent'
             stateSpacesMatch = .FALSE.
          endif
          if ( (.not. stateSpaceA%logSize) .and. stateSpaceB%logSize) then
             print*,trim(procName)//'error - logSize inconsistent', &
                  stateSpaceA%logSize,stateSpaceB%logSize
             stateSpacesMatch = .FALSE.
          endif
       endif
    endif
    if (present(sfcMW)) then
       if (sfcMW) then
          if (allocated(stateSpaceA%frqMW) .and. allocated(stateSpaceB%frqMW)) then
             if (stateSpaceA%nChanMW  /= stateSpaceB%nChanMW) then
                print*,trim(procName)//'error - nChanMW inconsistent', &
                     stateSpaceA%nChanMW,stateSpaceB%nChanMW
                stateSpacesMatch = .FALSE.
             endif
             if (ANY(abs(stateSpaceA%frqMW(1:stateSpaceA%nChanMW)- &
                  stateSpaceB%frqMW(1:stateSpaceB%nChanMW)) > deltaFrq )) then
                print*,trim(procName)//'error - MW Frequencies inconsistent'
                stateSpacesMatch = .FALSE.
             endif
             if (ANY(stateSpaceA%polMW(1:stateSpaceA%nChanMW) /= &
                     stateSpaceB%polMW(1:stateSpaceB%nChanMW))) then
                print*,trim(procName)//'error - MW polarizatoins inconsistent'
                stateSpacesMatch = .FALSE.
             endif
          endif
       endif
    endif
    if (present(sfcIR)) then
       if (sfcIR) then
          if ( allocated(stateSpaceA%wvnIR) .and. allocated(stateSpaceB%wvnIR)) then
             if (stateSpaceA%nEmIR  /= stateSpaceB%nEmIR) then
                print*,trim(procName)//'error - nEmIR inconsistent', &
                     stateSpaceA%nEmIR,stateSpaceB%nEmIR
                stateSpacesMatch = .FALSE.
             endif
             if ( ANY(stateSpaceA%wvnIR /= stateSpaceB%wvnIR) ) then
                print*,trim(procName)//'error - wvnIR inconsistent', &
                     stateSpaceA%wvnIR,stateSpaceB%wvnIR
                stateSpacesMatch = .FALSE.
             endif
          endif
       endif
    endif
    IGstatus = stateIndicesMatch(stateSpaceA%IG,stateSpaceB%IG,atmos=atmos,sfcMW=sfcMW,sfcIR=sfcIR)
    if (.not. IGstatus) stateSpacesMatch = .FALSE.
    NGstatus = stateIndicesMatch(stateSpaceA%NG,stateSpaceB%NG,atmos=atmos,sfcMW=sfcMW,sfcIR=sfcIR)
    if (.not. NGstatus) stateSpacesMatch = .FALSE.
    if (dbg) print *, trim(procName)//'ending ...'
  end function stateSpacesMatch

END MODULE ControlStructure
