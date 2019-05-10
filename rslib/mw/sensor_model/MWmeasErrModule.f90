MODULE MWmeasErrModule
  !
  ! Module for setting MW sensor noise
  !
  ! Subroutine:
  !
  !     MWnoiseInit
  !     setMWmeasErr
  !     MWnoiseDestroy
  !
  ! USE
  !
  !     constants
  !     NoiseModule
  !     MeasErrStructure
  !     MWobsStructure
  !     ChannelSelectionModule
  !     ControlStructure
  !
  ! yhe@aer.com, 03/16/2016
  !
  USE constants, Only : &
       MISSING_REAL

  USE NoiseModule, Only: &
       deviceNoise

  USE MeasErrStructure, Only: &
       MeasErr_t

  USE MWobsStructure, Only: &
       MWdefin_t, &
       MWob_t

  USE ChannelSelectionModule, Only: &
       ChannelSet_t

  USE ControlStructure, Only: &
       GenControl_t

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: &
       MWnoiseInit, &
       setMWmeasErr, &
       MWnoiseDestroy

  integer, parameter :: mxLength = 256
  integer, parameter :: iCell = 1
  TYPE(MeasErr_t) :: MWmeasErrLocal
  real, dimension(:), allocatable :: tbEquiv
  logical, dimension(:), allocatable :: kchan
  logical :: dbg = .false.
CONTAINS

  subroutine MWnoiseInit(MWdefin)
    TYPE(MWdefin_t), intent(inout) :: MWdefin
    !Local
    character (len=mxLength) :: procName

    procName = ' [MWmeasErrModule::MWnoiseInit]: '
    if (dbg) print *, trim(procName)//'starting ...'
    if (.not. allocated(MWmeasErrLocal%device)) allocate(MWmeasErrLocal%device(MWdefin%nChan))
    if (.not. allocated(MWmeasErrLocal%fwdModel)) allocate(MWmeasErrLocal%fwdModel(MWdefin%nChan))
    if (.not. allocated(tbEquiv)) allocate(tbEquiv(MWdefin%nChan))
    if (dbg) print *, trim(procName)//'ending ...'

  end subroutine MWnoiseInit

  subroutine setMWmeasErr(&
       genControl, &
       MWdefin, &
       MWob, &
       chanSet, &
       MWmeasErr &
       )
    TYPE(GenControl_t), intent(in) :: genControl
    TYPE(MWdefin_t), intent(in) :: MWdefin
    TYPE(MWOb_t), intent(inout) :: MWob
    TYPE(ChannelSet_t), intent(in) :: chanSet
    TYPE(MeasErr_t), intent(out) :: MWmeasErr

    !Local
    character (len=mxLength) :: procName
    logical :: DynamicNoise = .FALSE.
    real :: sumSq,dummy
    integer :: i,nChOn

    procName = ' [MWmeasErrModule::setMWmeasErr]: '
    if (dbg) print *, trim(procName)//'starting ...'

    if (.not. allocated(kchan)) allocate(kchan(chanSet%nMW))
    kchan = chanSet%MW == 1

    if (.not. allocated(MWmeasErr%device)) allocate(MWmeasErr%device(MWdefin%nChan))
    if (.not. allocated(MWmeasErr%fwdModel)) allocate(MWmeasErr%fwdModel(MWdefin%nChan))
    if (MWdefin%radOrTb == 0) then
       tbEquiv(1:MWdefin%nChan) = MWob%rad(1:MWdefin%nChan)
    else               ! Convert from radiance to radiance temperature
       tbEquiv(1:MWdefin%nChan) = &
            MWob%rad(1:MWdefin%nChan)*MWdefin%planckAlpha(1:MWdefin%nChan)+ &
            MWdefin%planckBeta(1:MWdefin%nChan)
    end if

    CALL deviceNoise(icellres=iCell, &
         tb=tbEquiv(1:MWdefin%nChan), &
         rerr=MWmeasErr%device(1:MWdefin%nChan), &
         sum_rerr=dummy, &
         kchan=kchan(1:MWdefin%nChan), &
         nchan=MWdefin%nChan, &
         DynamicNoise=DynamicNoise, &                !add scan position indext to arg list
         Fnedt=genControl%MWnoiseFile, &
         Fnrf=genControl%MWnoiseAvgFactorFile &
         )

    MWmeasErr%device(1:MWdefin%nChan) = sqrt(MWmeasErr%device(1:MWdefin%nChan)**2+(0.05)**2)

    IF (MWdefin%radOrTb == 1) &
         MWmeasErr%device(1:MWdefin%nChan) = &
         MWmeasErr%device(1:MWdefin%nChan)/MWdefin%planckAlpha(1:MWdefin%nChan)

    nChOn=0
    sumSq=0.
    do i=1,MWdefin%nChan
       if (kchan(i)) then
          nChOn=nChOn+1
          sumSq=sumSq+MWmeasErr%device(i)**2
       endif
    enddo
    if (nChOn == 0) then
       MWmeasErr%rms = MISSING_REAL
    else
       MWmeasErr%rms=SQRT(sumSq/float(nChOn))
    endif

    MWmeasErr%fwdModel(1:MWdefin%nChan)=0.    !goes in MWmeasErr%fwdModel

    MWmeasErrLocal%device = MWmeasErr%device
    MWmeasErrLocal%fwdModel = MWmeasErr%fwdModel
    if (dbg) print *, trim(procName)//'ending ...'

  end subroutine setMWmeasErr

  subroutine MWnoiseDestroy()
    !Local
    character (len=mxLength) :: procName

    procName = ' [MWmeasErrModule::MWnoiseDestroy]: '
    if (dbg) print *, trim(procName)//'starting ...'
    if (allocated(MWmeasErrLocal%device)) deallocate(MWmeasErrLocal%device)
    if (allocated(MWmeasErrLocal%fwdModel)) deallocate(MWmeasErrLocal%fwdModel)
    if (allocated(tbEquiv)) deallocate(tbEquiv)
    if (allocated(kchan)) deallocate(kchan)
    if (dbg) print *, trim(procName)//'ending ...'

  end subroutine MWnoiseDestroy

END MODULE MWmeasErrModule
