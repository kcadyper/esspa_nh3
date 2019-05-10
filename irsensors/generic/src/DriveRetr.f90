!<f90File>**************************************************************
!
! CONTACT:
!
!   Atmospheric & Environmental Research, Inc
!   131 Hartwell Ave
!   Lexington ,MA 02421-3126 USA
!   Phone: 781.761.2288
!   E-mail: yhe@aer.com
!
! COPYRIGHT NOTICE:
!
!   Copyright AER, Inc 2016-, All Rights Reserved
!   See the file README-DATARIGHTS.txt included with this release
!   for additional details.
!
!*************************************************************</f90File>

program DriveRetr

  ! <f90Main>*************************************************************
  !
  ! NAME:
  !
  !   DriveRetr
  !
  ! PURPOSE:
  !
  !   Main driver for calling all the subprocesses for the retrieval
  !
  ! SUBROUTINES:
  !
  !   initDriveRetr
  !   destroyDriveRetr
  !
  ! USE
  !
  !   RetrModule
  !   OutputStructure
  !   IOmodule
  !   StateIndexModule
  !   StatRetrModule
  !   ControlStructure
  !   MWobsStructure
  !   IRobsStructure
  !   ImgObsStructure
  !   MWretrModule
  !   PrimaryRetrModule
  !   SecondaryRetrModule
  !   MWRTmodule
  !   IRRTmodule
  !   ImgRTModule
  !   AncillaryStructure
  !   ChannelSelectionModule
  !   GranuleMemModule
  !
  !
  !*************************************************************</f90Main>
  USE RetrModule, Only: &
       getControl, &
       init, &
       retr, &
       destroy

  USE OutputStructure, Only: &
       OutputGran_t, &
       OutputDefin_t, &
       LinearDefin_t, &
       LinearGran_t

  USE IOModule, Only: &
       getMWsensorData, &
       getIRsensorData, &
       getImgSensorData, &
       getAncillaryData, &
       putRetr, &
       putLinear

  USE StatRetrModule, Only: &
       statRetrControl_t

  USE ControlStructure, Only: &
       GenControl_t, &
       StateSpaceDefin_t, &
       RetrResult_t, &
       LinearResult_t

  USE MWobsStructure, Only: &
       MWgran_t

  USE IRobsStructure, Only: &
       IRgran_t

  USE ImgObsStructure, Only: &
       ImgDef_t, &
       ImgOb_t, &
       ImgGran_t

  USE MWretrModule, Only: &
       MWretrControl_t

  USE PrimaryRetrModule, Only: &
       PrimaryRetrControl_t

  USE SecondaryRetrModule, Only: &
       SecondaryRetrControl_t

  USE MWRTmodule, Only: &
       MWRTcontrol_t

  USE IRRTmodule, Only: &
       IRRTcontrol_t

  USE ImgRTmodule, Only: &
       imgRTcontrol_t

  USE AncillaryStructure, Only: &
       AncillaryGran_t

  USE ChannelSelectionModule, Only: &
       ChannelSelections_t

  USE GranuleMemModule, Only: &
       granAlloc, &
       granDestroy

  IMPLICIT NONE

  INTEGER, PARAMETER :: mxLength = 256
  TYPE(GenControl_t)           :: genControl
  TYPE(statRetrControl_t)      :: statRetrControl
  TYPE(MWretrControl_t)        :: MWretrControl
  TYPE(primaryRetrControl_t)   :: primaryRetrControl
  TYPE(secondaryRetrControl_t) :: secondaryRetrControl
  TYPE(MWRTcontrol_t)          :: MWRTcontrol
  TYPE(IRRTcontrol_t)          :: IRRTcontrol
  TYPE(imgRTcontrol_t)         :: imgRTcontrol
  TYPE(MWgran_t)               :: MWgran
  TYPE(IRgran_t)               :: IRgran
  TYPE(ImgGran_t)              :: ImgGran
  TYPE(AncillaryGran_t)        :: ancillaryGran
  TYPE(StateSpaceDefin_t)      :: stateSpaceDefinStat
  TYPE(StateSpaceDefin_t)      :: stateSpaceDefinMW
  TYPE(StateSpaceDefin_t)      :: stateSpaceDefinPrimary
  TYPE(StateSpaceDefin_t)      :: stateSpaceDefinSecondary
  TYPE(ChannelSelections_t)    :: ChannelSelections
  TYPE(OutputGran_t) :: retrOutputStat
  TYPE(OutputGran_t) :: retrOutputMW
  TYPE(OutputGran_t) :: retrOutputPrimary
  TYPE(OutputGran_t) :: retrOutputSecondary
  TYPE(OutputGran_t) :: backgroundOutput
  TYPE(LinearGran_t) :: linearOutput
  integer :: iFOR, nFOR
  integer :: iFOR_MW, iFOR_IR, iFOR_Img
  character (len=mxLength)     :: procName
  logical :: dbg=.true.

  procName = " [DriveRetr]: "
  if (dbg) print *, trim(procName)//"starting ..."

  call getControl(&
       genControl, &
       statRetrControl, &
       MWretrControl, &
       primaryRetrControl, &
       secondaryRetrControl, &
       MWRTcontrol, &
       IRRTcontrol, &
       imgRTcontrol)

  call getMWsensorData(genControl, MWgran)
  call getIRsensorData(genControl, IRgran)
  call getImgSensorData(genControl, ImgGran)
  call getAncillaryData(genControl, ancillaryGran)
  call init( &
       genControl, &
       statRetrControl, &
       MWretrControl, &
       primaryRetrControl, &
       secondaryRetrControl, &
       MWRTcontrol, &
       IRRTcontrol, &
       imgRTcontrol, &
       MWgran%def, &
       IRgran%def, &
       ancillaryGran%def, &
       stateSpaceDefinStat, &
       stateSpaceDefinMW, &
       stateSpaceDefinPrimary, &
       stateSpaceDefinSecondary, &
       ChannelSelections &
       )


  call granAlloc( &
       genControl, &
       MWRTcontrol, &
       primaryRetrControl%genLinInvert, &
       IRgran, &
       MWgran, &
       stateSpaceDefinStat, &
       stateSpaceDefinMW, &
       stateSpaceDefinPrimary, &
       stateSpaceDefinSecondary, &
       retrOutputStat, &
       retrOutputMW, &
       retrOutputPrimary, &
       LinearOutput, &
       backgroundOutput, &
       retrOutputSecondary &
       )

  nFOR = genControl%nFORmax
  if (genControl%statRetrOn)      nFOR = min(nFOR,retrOutputStat%nFOR)
  if (genControl%MWretrOn)        nFOR = min(nFOR,retrOutputMW%nFOR)
  if (genControl%primaryRetrOn) then
     nFOR = min(nFOR,retrOutputPrimary%nFOR)
     nFOR = min(nFOR,backgroundOutput%nFOR)
     if (primaryRetrControl%genLinInvert) nFOR = min(nFOR,linearOutput%nFOR)
  endif
  if (genControl%secondaryRetrOn) nFOR = min(nFOR,retrOutputSecondary%nFOR)
  nFOR = max(0,nFOR)

  FORLoop: do iFOR = 1, nFOR
  !FORLoop: do iFOR = 1,50
     !if (dbg) print *, trim(procName)//"iFOR = ", iFOR
     iFOR_MW = 1  ! Default for memory control
     iFOR_IR = 1  ! Default for memory control
     iFOR_Img = 1 ! Default for memory control
     if (genControl%MWdataOn) iFOR_MW = iFOR
     if (genControl%IRdataOn) iFOR_IR = iFOR
     if (genControl%ImgDataOn) iFOR_Img = iFOR
     call retr(&
          genControl, &
          statRetrControl, &
          MWretrControl, &
          primaryRetrControl, &
          secondaryRetrControl, &
          MWRTcontrol, &
          IRRTcontrol, &
          imgRTcontrol, &
          MWgran%def, &
          MWgran%ob(iFOR_MW), &
          IRgran%def, &
          IRgran%FOR(iFOR_IR), &
          ImgGran%def, &
          ImgGran%ob(:,iFOR_Img), &
          ancillaryGran%def, &
          ancillaryGran%FOR(iFOR), &
          channelSelections, &
          retrOutputStat%FOR2D(:,iFOR), &
          retrOutputMW%FOR(iFOR), &
          retrOutputPrimary%FOR2D(:,iFOR), &
          linearOutput%FOR(iFOR), &
          backgroundOutput%FOR(iFOR), &
          retrOutputSecondary%FOR2D(:,iFOR) &
          )
  end do FORLoop

  !Statistical
  if (genControl%statRetrOn) then
     call putRetr(&
          genControl%statRetrOutFile, &
          stateSpaceDefinStat, &
          retrOutputStat &
          )
  end if

  !MW
  if (genControl%MWretrOn) &
       call putRetr(&
       genControl%MWretrOutFile, &
       stateSpaceDefinMW, &
       retrOutputMW &
       )

  !Primary
  if (genControl%primaryRetrOn) &
       call putRetr(&
       genControl%primaryRetrOutFile, &
       stateSpaceDefinPrimary, &
       retrOutputPrimary &
       )

  !background
  if (genControl%primaryRetrOn) &
       call putRetr(&
       genControl%backgroundOutFile, &
       stateSpaceDefinPrimary, &
       backgroundOutput &
       )

  !Linear
  if (primaryRetrControl%genLinInvert) &
       call putLinear(&
       genControl%linearInvertOutFile, &
       stateSpaceDefinPrimary, &
       primaryRetrControl, &
       linearOutput &
       )

  !Secondary
  if (genControl%secondaryRetrOn) &
       call putRetr(&
       genControl%secondaryRetrOutFile, &
       stateSpaceDefinSecondary, &
       retrOutputSecondary &
       )

  call destroy(genControl, &
       MWgran, &
       IRgran, &
       imgGran, &
       ancillaryGran, &
       ChannelSelections &
       )

   call granDestroy(&
        retrOutputStat, &
        retrOutputMW, &
        retrOutputPrimary, &
        LinearOutput, &
        backgroundOutput, &
        retrOutputSecondary &
        )

  if (dbg) print *, trim(procName)//"ending ..."

END PROGRAM DriveRetr
