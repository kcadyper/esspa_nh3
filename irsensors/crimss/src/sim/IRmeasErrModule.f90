MODULE IRmeasErrModule
  !
  ! Module for setting IR sensor noise
  !
  ! Subroutines:
  !
  !     IRnoiseInit
  !     setIRmeasErr
  !     IRnoiseDestroy
  !
  ! USE:
  !
  !     IRNoiseModule
  !
  ! yhe@aer.com, 03/16/2016
  !
  USE ControlStructure, Only: &
       GenControl_t

  USE IRobsStructure, Only: &
       IRdefin_t, &
       IRob_t

  USE MeasErrStructure, Only: &
       MeasErr_t

  USE IRNoiseModule, Only: &
       initIRdeviceNoise, &
       IRdeviceNoise, &
       destroyIRdeviceNoise

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: &
       IRnoiseInit, &
       setIRmeasErr, &
       IRnoiseDestroy

  integer, parameter :: mxLength = 256
  TYPE(MeasErr_t) :: IRmeasErrLocal

  logical, dimension(:), allocatable :: kchan
  logical :: dbg = .false.

CONTAINS

  subroutine IRnoiseInit( &
       genControl, &
       IRdefin, &
       nchanRT &
       )
    TYPE(GenControl_t), intent(in) :: genControl
    TYPE(IRdefin_t), intent(in) :: IRdefin
    integer, intent(in) :: nChanRT
    !Local
    character (len=mxLength) :: procName

    procName = ' [IRmeasErrModule::IRnoiseInit]: '
    if (dbg) print *, trim(procName)//'starting ...'
    !nChanRT or IRdefin%nChan?
    !design: "construct IRmeasErr structure from size nchanIR" (nchanRT?)
    !What about DevNoise and AtmNoise above? Are the local?
    if (.not. allocated(IRmeasErrLocal%device)) allocate(IRmeasErrLocal%device(nChanRT))
    if (.not. allocated(IRmeasErrLocal%fwdModel)) allocate(IRmeasErrLocal%fwdModel(nChanRT))
    call initIRdeviceNoise(nChanRT,genControl%IRnoiseFile)
    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine IRnoiseInit

  subroutine setIRmeasErr(&
       IRdefin, &
       IRob, &
       IRmeasErr &
       )
    TYPE(IRdefin_t), intent(in) :: IRdefin
    TYPE(IROb_t), intent(inout) :: IRob
    TYPE(MeasErr_t), intent(out) :: IRmeasErr

    !Local
    character (len=mxLength) :: procName

    procName = ' [IRmeasErrModule::setIRmeasErr]: '
    if (dbg) print *, trim(procName)//'starting ...'

    if (.not. allocated(kchan)) allocate(kchan(IRdefin%nChan))
    kchan = .true. ! This assumes that channel selection can be disregared here

    if (.not. allocated(IRmeasErr%device)) allocate(IRmeasErr%device(IRdefin%nChan))
    if (.not. allocated(IRmeasErr%fwdModel)) allocate(IRmeasErr%fwdModel(IRdefin%nChan))
    call IRdeviceNoise( &
         IRdefin%wvn, &
         IRob%rad, &
         IRmeasErr%device, &
         IRmeasErr%fwdModel, &
         IRmeasErr%rms, &
         kchan, &
         nchan=IRdefin%nChan, &
         chanIDreq=IRdefin%chanID, &
         chanNum=IRdefin%chanNum &
         )
    if (allocated(IRdefin%NEdN)) then
      !print *,'use NEdN from L1B file'
      ! Overwrite device noise and rms. For now, retain fwdModel noise from
      ! IRdeviceNoise. Later, forward model error could be produced in a 
      ! separate subroutine, with compatible revisions in all software that
      ! calls IRdeviceNoice
      IRmeasErr%device(:) = sum(IRdefin%NEdN, dim=2)/real(IRdefin%nFOV)
      IRmeasErr%rms = sqrt(dot_product(IRmeasErr%device(:),IRmeasErr%device(:)) &
                           /real(IRdefin%nChan))
    else
      print *,'use NEdN from IRmeasErr%device file'
    end if

    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine setIRmeasErr

  subroutine IRnoiseDestroy()
    !Local
    character (len=mxLength) :: procName

    procName = ' [IRmeasErrModule::IRnoiseDestroy]: '
    if (dbg) print *, trim(procName)//'starting ...'
    call destroyIRdeviceNoise()
    if (allocated(IRmeasErrLocal%device)) deallocate(IRmeasErrLocal%device)
    if (allocated(IRmeasErrLocal%fwdModel)) deallocate(IRmeasErrLocal%fwdModel)
    if (allocated(kchan)) deallocate(kchan)
    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine IRnoiseDestroy

END MODULE IRmeasErrModule
