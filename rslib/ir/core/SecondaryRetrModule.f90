MODULE SecondaryRetrModule
  !
  ! Module for executing MW retrieval process
  !
  ! Subroutine:
  !
  !     secondaryRetrInit
  !     secondaryRetr
  !     secondaryRetrDestroy
  !
  ! Derived data type:
  !
  !     SecondaryRetrControl_t
  !
  ! USE:
  !
  !     ControlStructure
  !     StateIndexModule
  !     IRRTmodule
  !     IRobsStructure
  !     IRmeasErrModule
  !     AncillaryStructure
  !     ChannelSelectionModule
  !
  ! History
  !     initial, 03/16/2016, yhe@aer.com
  !
  USE ControlStructure, Only: &
       GenControl_t, &
       RetrResult_t, &
       StateSpaceDefin_t

  USE StateIndexModule, Only: &
       StateIndex_t

  USE IRRTmodule, Only: &
       IRRTcontrol_t

  USE IRobsStructure, Only: &
       IRdefin_t, &
       IRob_t

  USE MeasErrStructure, Only: &
       MeasErr_t

  USE AncillaryStructure, Only: &
       AncillaryDefin_t, &
       AncillaryData_t

  USE ChannelSelectionModule, Only: &
       ChannelSet_t

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: &
       secondaryRetrInit, &
       secondaryRetr, &
       secondaryRetrDestroy, &
       SecondaryRetrControl_t

  integer, parameter :: mxLength = 256
  logical :: dbg = .true.

  TYPE SecondaryRetrControl_t
     TYPE(StateIndex_t) :: nRetr
  END TYPE SecondaryRetrControl_t

CONTAINS

  subroutine secondaryRetrInit( &
       genControl, &
       secondaryRetrControl, &
       IRRTcontrol, &
       stateSpaceDefin &
       )

    TYPE(GenControl_t), intent(in) :: genControl
    TYPE(SecondaryRetrControl_t), intent(in) :: secondaryRetrControl
    TYPE(IRRTcontrol_t), intent(in) :: IRRTcontrol
    TYPE(StateSpaceDefin_t), intent(inout) :: stateSpaceDefin

    !Local
    character (len=mxLength) :: procName

    procName = ' [SecondaryRetrModule::secondaryRetrInit]: '
    if (dbg) print *, trim(procName)//' starting ...'
    if (dbg) print *, trim(procName)//' to be implemented...'
    if (dbg) print *, trim(procName)//' ending ...'

  end subroutine secondaryRetrInit

  subroutine secondaryRetr(&
       genControl, &
       secondaryRetrControl, &
       IRdefin, &
       IRob, &
       IRmeasErr, &
       ancillaryDefin, &
       ancillaryOb, &
       channelSet, &
       retrPrior, &
       retrResult)
    TYPE(GenControl_t), intent(in) :: genControl
    TYPE(SecondaryRetrControl_t), intent(in) :: secondaryRetrControl
    TYPE(IRdefin_t), intent(in) :: IRdefin
    TYPE(IRob_t), intent(in) :: IRob
    TYPE(MeasErr_t) :: IRmeasErr
    TYPE(AncillaryDefin_t), intent(in) :: ancillaryDefin
    TYPE(AncillaryData_t), intent(in) :: ancillaryOb
    TYPE(ChannelSet_t), intent(in) :: ChannelSet
    TYPE(RetrResult_t), intent(in) :: retrPrior
    TYPE(RetrResult_t), intent(inout) :: retrResult

    !Local
    character (len=mxLength) :: procName

    procName = ' [SecondaryRetrModule::secondaryRetr]: '
    if (dbg) print *, trim(procName)//' starting ...'

    if (channelSet%nIR < 0) then
        print *, 'err: ', procName, 'Expected IR channel Selection for ', &
           'for Secondary retrieval not supplied'
        call exit(1)
    end if

    if (dbg) print *, trim(procName)//' to be implemented...'
    if (dbg) print *, trim(procName)//' ending ...'

  end subroutine secondaryRetr

  subroutine secondaryRetrDestroy()
    !Local
    character (len=mxLength) :: procName

    procName = ' [SecondaryRetrModule::secondaryRetrDestroy]: '
    if (dbg) print *, trim(procName)//' starting ...'
    if (dbg) print *, trim(procName)//' to be implemented...'
    if (dbg) print *, trim(procName)//' ending ...'

  end subroutine secondaryRetrDestroy

END MODULE SecondaryRetrModule
