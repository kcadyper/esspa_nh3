MODULE MWRTmodule
  !
  ! Module for executing MW RTM
  !
  ! Subroutine:
  !
  !     MWRTinit
  !     MWRTexec (calls ossdrv_mw)
  !     MWRTdestroy
  !
  ! Derived data type:
  !
  !     MWRTcontrol_t
  !
  ! USE:
  !
  !     DimParameters
  !     oss_mw_module
  !
  ! yhe@aer.com, 03/16/2016
  !
  USE ToolboxModule, Only: &
       getUnit

  USE StateIndexModule, Only: &
       maxMol, &
       whereH2O, &
       getAtmosVectorLength

  USE ControlStructure, Only: &
       StateSpaceDefin_t

  USE MWobsStructure, Only: &
       MWdefin_t

  USE oss_mw_module, Only: &
       ossinit_mw, &
       ossdrv_mw, &
       oss_destroy

  USE OSSMWmoduleSubs, Only: &
       getPlanckLinCtr, &
       set2ndOrderPlanck

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: &
       MWRTinit, &
       MWRTexec, &
       MWRTdestroy, &
       MWRTcontrol_t

  integer, parameter :: mxLength = 256
  logical :: dbg = .false.
  TYPE MWRTcontrol_t
     character(len=mxLength) :: coefFile
     character(len=mxLength) :: LUTfile
     integer :: nMol
     integer, dimension(maxMol) :: molID
     integer :: radOrTb
     logical :: Planck2ndOrderTaylor
  END TYPE MWRTcontrol_t

CONTAINS

  subroutine MWRTinit(MWRTcontrol,stateSpaceDefin,enableScatter, &
       MWdefin,nChanRT,frqRT,polRT,nLevRT,pRefRT)

    TYPE(MWRTcontrol_t), intent(inout) :: MWRTcontrol
    TYPE(StateSpaceDefin_t), intent(in) :: stateSpaceDefin
    TYPE(MWdefin_t), intent(inout) :: MWdefin
    logical, intent(in) :: enableScatter
    real, dimension(:), intent(out) :: frqRT
    real, dimension(:), allocatable, intent(out) :: pRefRT
    integer, dimension(:), intent(out) :: polRT
    integer, intent(out) :: nchanRT
    integer, intent(out) :: nLevRT

    !Local
    character (len=*), parameter :: procName= ' [MWRTModule::MWRTInit]: '
    integer :: ih2o
    integer :: NParGAtm

    if (dbg) print *, procName, 'starting ...'
    ih2o=whereH2O(MWRTcontrol%MolID)
    NParGAtm=getAtmosVectorLength(stateSpaceDefin%NG)
    MWRTcontrol%radOrTb=MWdefin%radOrTb

    call ossinit_mw(MWRTcontrol%coefFile, &
         MWRTcontrol%LUTfile, &
         ih2o, &
         MWRTcontrol%MolID(1:ih2o),&
         NParGAtm, &
         nchanRT, &
         frqRT, &
         nLevRT, &
         pRefRT, &
         polRT, &
         MWRTcontrol%radOrTb)

!!$       allocate(alphaCtr(MWdefin%nChan),betaCtr(MWdefin%nChan))
!!$       CALL getPlanckLinCtr(alphaCtr,betaCtr)

       ! already allocated in IOmodule::getMWsensorData
!!$       allocate(MWdefin%planckAlpha(MWdefin%nChan),MWdefin%planckBeta(MWdefin%nChan))
    call getPlanckLinCtr(MWdefin%planckAlpha,MWdefin%planckBeta)
    if ( MWRTcontrol%Planck2ndOrderTaylor .and. (MWRTcontrol%radOrTb==0)) &
                                                  call set2ndOrderPlanck()
    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine MWRTinit

  subroutine MWRTexec()
    !Local
    character (len=mxLength) :: procName
    procName = ' [MWRTModule::MWRTexec]: '
    if (dbg) print *, trim(procName)//'starting ...'
    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine MWRTexec

  subroutine MWRTdestroy()
    !Local
    character (len=mxLength) :: procName
    procName = ' [MWRTModule::MWRTdestroy]: '
    if (dbg) print *, trim(procName)//'starting ...'
    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine MWRTdestroy

END MODULE MWRTmodule
