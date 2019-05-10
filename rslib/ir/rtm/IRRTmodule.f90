MODULE IRRTmodule
  !
  ! Module for executing IR RTM
  !
  ! Subroutine:
  !
  !     IRRTinit
  !     IRRTexec (calls ossdrv_mw)
  !     IRRTdestroy
  !
  ! Derived data type:
  !
  !     IRRTcontrol_t
  !
  ! USE:
  !
  !     DimParameters
  !     oss_ir_module
  !
  ! yhe@aer.com, 03/16/2016
  !
  USE ControlStructure, Only: &
       StateSpaceDefin_t

  USE StateIndexModule, Only: &
       maxMol, &
       getAtmosVectorLength

  USE oss_ir_module_scat, Only: &
       ossdrv_ir, &
       ossinit_ir

  USE ToolboxModule, Only: &
       getUnit

  USE CloudDataStruct, only : &
       CloudScene

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: &
       IRRTinit, &
       IRRTexec, &
       IRRTdestroy, &
       IRRTcontrol_t, &
       IRRTsetActiveChannelSelectionSet

  integer, parameter :: mxLength = 256
  logical :: dbg = .false.
  TYPE IRRTcontrol_t
     character(len=mxLength) :: coefFile
     character(len=mxLength) :: LUTfile
     character(len=mxLength) :: auxFile
     character(len=mxLength) :: defaultProfFile
     character(len=mxLength) :: cloudOptFile
     character(len=mxLength) :: solarFile
     character(len=mxLength) :: scatterCoefFile
     integer :: nMol
     integer, dimension(maxMol) :: molID
     logical :: scatterApprox
     integer :: nHyd
     !-- size of number of hydrometeors
     character(len=mxLength), dimension(:), allocatable :: dataTypeCloud
     !-- cloud type IDs, size of number of hydrometeors
     integer, dimension(:), allocatable :: cloudModel
  END TYPE IRRTcontrol_t

CONTAINS

  subroutine IRRTinit( &
       IRRTcontrol, &
       stateSpaceDefin, &
       enableScatter, & ! (from genControl%enableScatter)
       ChannelSelections, &
       nChanRT, &
       wvnRT, &
       nLevRT, &
       pRefRT &
       )

    use CloudModule, only : &
        typeSlab
    use ChannelSelectionModule, only : &
        ChannelSelections_t
    USE OSSIRmoduleSubs, only: &
        addChanSelectMask

    TYPE(IRRTcontrol_t),             intent(in)   :: IRRTcontrol
    TYPE(StateSpaceDefin_t),         intent(in)   :: stateSpaceDefin
    logical,                         intent(in)   :: enableScatter
    TYPE(ChannelSelections_t),       intent(inout):: ChannelSelections
    integer,                         intent(out)  :: nChanRT
    integer,                         intent(out)  :: nLevRT
    real, dimension(:),              intent(out)  :: wvnRT
    real, dimension(:), allocatable, intent(out)  :: pRefRT

    !Local
    integer :: NParGAtm
    integer, allocatable         :: chanIndex(:)
    LOGICAL                      :: ossScatRegrIRon=.FALSE. !should be loaded in
    integer                      :: cldTyp=typeSlab !'slab' !should be loaded in
    character (len=*), parameter :: procName=' [IRRTmodule::IRRTinit]: '
     !-- cloud scene, size of number of hydrometeors
    type(CloudScene) , dimension(:), allocatable :: cldScene
    CHARACTER (LEN=7), DIMENSION(:), allocatable :: propertyNames
    integer, parameter                           :: npProperty = 3
    integer                                      :: iHydr, ipp
    integer, parameter                           :: nLevHydr = 2

    !RetrModule:genControl, but
    !currently not associated with
    !any data structure


    if (dbg) print *, procName, 'starting ...'

    allocate(cldScene(IRRTcontrol%nHyd))
    ! npProperty has to be consistent with the number of elements in the
    ! propertyNames

    allocate(propertyNames(npProperty))
    propertyNames(1:npProperty) = (/'SizP', 'SizS', 'Temp'/)

    do iHydr = 1, IRRTcontrol%nHyd
       cldScene(iHydr)%name = IRRTcontrol%dataTypeCloud(iHydr)
       cldScene(iHydr)%nProperty = npProperty
       DO ipp = 1, npProperty
          cldScene(iHydr)%physDescr(ipp)%name = trim(propertyNames(ipp))
       ENDDO
    end do


    !----Initialize IR OSS module and get back Info about LUTs
    NParGAtm=getAtmosVectorLength(stateSpaceDefin%NG)

    call ossinit_ir( &
         IRRTcontrol%coefFile, &
         IRRTcontrol%LUTfile, &
         IRRTcontrol%solarFile, &
         IRRTcontrol%cloudOptFile, &
         IRRTcontrol%scatterCoefFile, &
         cldTyp, &
         IRRTcontrol%nHyd, &
         cldScene, &
         enableScatter, &
         ossScatRegrIRon, &
         IRRTcontrol%defaultProfFile, &
         IRRTcontrol%nMol, &
         IRRTcontrol%molID(1:IRRTcontrol%nMol), &
         nchanRT, &
         chanIndex, &
         wvnRT, &
         nLevRT, &
         pRefRT)

    ! add channel selections
      if (ChannelSelections%CloudClearing%nIR>0) &
          call addChanSelectMask(ChannelSelections%CloudClearing%IR, &
                                  ChannelSelections%CloudClearing%nIR, &
                                  chSetIDout=ChannelSelections%CloudClearing%chanSetInd)

      if (ChannelSelections%PrimaryClear%nIR>0) &
          call addChanSelectMask(ChannelSelections%PrimaryClear%IR, &
                                  ChannelSelections%PrimaryClear%nIR, &
                                  chSetIDout=ChannelSelections%PrimaryClear%chanSetInd)

      if (ChannelSelections%PrimaryCloudy%nIR>0) &
          call addChanSelectMask(ChannelSelections%PrimaryCloudy%IR, &
                                  ChannelSelections%PrimaryCloudy%nIR, &
                                  chSetIDout=ChannelSelections%PrimaryCloudy%chanSetInd)

      if (ChannelSelections%Secondary%nIR>0) &
          call addChanSelectMask(ChannelSelections%Secondary%IR, &
                                  ChannelSelections%Secondary%nIR, &
                                  chSetIDout=ChannelSelections%Secondary%chanSetInd)


    if (allocated(chanIndex)) deallocate(chanIndex)
    if (allocated(propertyNames)) deallocate(propertyNames)
    if (allocated(cldScene)) deallocate(cldScene)
    if (dbg) print *, trim(procName)//'ending ...'

  end subroutine IRRTinit

  subroutine IRRTexec()
    character (len=mxLength) :: procName

    procName = ' [IRRTModule::IRRTexec]: '
    if (dbg) print *, trim(procName)//'starting ...'
    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine IRRTexec

!
! sets a channel selection set active. The set must be added
! during call of IRRTinit.
!
  subroutine IRRTsetActiveChannelSelectionSet(chanSet)
    use OSSIRmoduleSubs, only: &
        setChanSelect, getSelectionIndex
    use ChannelSelectionModule, only : &
        ChannelSet_t

    type(ChannelSet_t), intent(in) :: chanSet

    character (len=*), parameter :: procName = &
                        ' [IRRTModule::IRRTsetActiveChannelSelectionSet]: '

    if (dbg) print *, trim(procName)//'starting ...'

    if (chanSet%nIR<0) then
        print *, 'err: ', procName, 'IR channel Selection is not supplied'
        call exit(1)
    end if

    if (getSelectionIndex() /= chanSet%chanSetInd) &
        call setChanSelect(chanSet%chanSetInd)

    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine IRRTsetActiveChannelSelectionSet

  subroutine IRRTdestroy()
    character (len=mxLength) :: procName

    procName = ' [IRRTModule::IRRTdestroy]: '
    if (dbg) print *, trim(procName)//'starting ...'
    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine IRRTdestroy

END MODULE IRRTmodule
