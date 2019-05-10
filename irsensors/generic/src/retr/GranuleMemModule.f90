Module GranuleMemModule
  !
  ! Memory management module
  !
  ! History:
  !
  !   07/25/2016, initial, yhe@aer.com
  !
  USE ControlStructure, Only: &
       GenControl_t, &
       StateSpaceDefin_t

  USE MWRTmodule, Only: &
       MWRTcontrol_t

  USE StateIndexModule, Only: &
       StateIndex_t, &
       initLengths, &
       getVectorLength, &
       genIndices

  USE scene_io_module, ONLY: &
       TB_UNITS, &
       RADMW_UNITS, &
       RAD_UNITS

  USE OutputStructure, Only: &
       OutputGran_t, &
       OutputDefin_t, &
       LinearGran_t, &
       LinearDefin_t

  USE MWobsStructure, Only: &
       MWgran_t, &
       MWdefin_t

  USE IRobsStructure, Only: &
       IRgran_t, &
       IRdefin_t

  USE LinInvert, Only: &
       getLinInvData

  IMPLICIT None
  Private
  Public :: &
       granAlloc, &
       granDestroy

  integer, parameter :: mxLength=256
  logical :: dbg = .false.

Contains
  subroutine granAlloc( &
       genControl, &
       MWRTcontrol, &
       genLinear, &
       IRgran, &
       MWgran, &
       stateSpaceDefinStat, &
       stateSpaceDefinMW, &
       stateSpaceDefinPrimary, &
       stateSpaceDefinSecondary, &
       retrOutputStat, &
       retrOutputMW, &
       retrOutputPrimary, &
       linearOutput, &
       backgroundOutput, &
       retrOutputSecondary &
       )
    TYPE(GenControl_t), intent(in) :: genControl
    TYPE(MWRTcontrol_t), intent(in) :: MWRTcontrol
    logical,            intent(in) :: genLinear
    TYPE(IRgran_t), intent(in)     :: IRgran
    TYPE(MWgran_t), intent(in)     :: MWgran
    TYPE(StateSpaceDefin_t)        :: stateSpaceDefinStat
    TYPE(StateSpaceDefin_t)        :: stateSpaceDefinMW
    TYPE(StateSpaceDefin_t)        :: stateSpaceDefinPrimary
    TYPE(StateSpaceDefin_t)        :: stateSpaceDefinSecondary
    TYPE(OutputGran_t), intent(inout) :: retrOutputStat
    TYPE(OutputGran_t), intent(inout) :: retrOutputMW
    TYPE(OutputGran_t), intent(inout) :: retrOutputPrimary
    TYPE(OutputGran_t), intent(inout) :: backgroundOutput
    TYPE(LinearGran_t), intent(inout) :: linearOutput
    TYPE(OutputGran_t), intent(inout) :: retrOutputSecondary
    !Local
    integer :: nFOR
    integer :: nIRFOV
    integer :: nChanMW
    integer :: nChanIR
    TYPE(StateSpaceDefin_t)        :: unitStateSpaceDefin
    TYPE(StateSpaceDefin_t)        :: stateSpaceDefinLoc
    character (len=mxLength) :: procName

    procName = "[GranuleMemModule::memAlloc]:"
    if (dbg) print *, trim(procName)//" starting ..." 

    unitStateSpaceDefin%NG = initLengths(NemMW=1,NemRfIR=1)
    unitStateSpaceDefin%nLev = 1

    !-This is the way the ESDR is designed to do
    !Notes: How about imgObs%nFOR,ancillaryGran%nFOR?
!!$  nFOR = max(0,maxval((/genControl%nFORmax,MWgran%nFOV,IRgran%nFOR/)))

    !-This is the way for debugging when genControl%nFORmax is the smallest
    nFOR = genControl%nFORmax
    if (genControl%MWdataOn) nFOR = min(nFOR,MWgran%nFOV)
    if (genControl%IRdataOn) nFOR = min(nFOR,IRgran%nFOR)
    nFOR = max(0,nFOR)
    if (dbg) print *, trim(procName)//"nFOR = ", nFOR

    if (genControl%MWdataOn) then
       nChanMW = MWgran%def%nChan
    else
       nChanMW = 0
    endif

    if (genControl%IRdataOn) then
       nIRFOV = IRgran%def%nFOV
       nChanIR = IRgran%def%nChan
    else
       nIRFOV = 1
       nChanIR = 0
    endif

    retrOutputStat%nFOR = nFOR
    if (genControl%statRetrOn) then
       stateSpaceDefinLoc=stateSpaceDefinStat
    else
       stateSpaceDefinLoc=unitStateSpaceDefin
    endif
    call defGranule(retrOutputStat%def,nIRFOV,stateSpaceDefinLoc,"retrOutputStat%def")

    retrOutputMW%nFOR = nFOR
    if (genControl%MWretrOn) then
       stateSpaceDefinLoc=stateSpaceDefinMW
    else
       stateSpaceDefinLoc=unitStateSpaceDefin
    endif
    call defGranule(retrOutputMW%def,1,stateSpaceDefinLoc,"retrOutputMW%def")

    retrOutputPrimary%nFOR = nFOR
    backgroundOutput%nFOR = nFOR
    if (genControl%primaryRetrOn) then
       stateSpaceDefinLoc=stateSpaceDefinPrimary
    else
       stateSpaceDefinLoc=unitStateSpaceDefin
    endif
    call defGranule(retrOutputPrimary%def,nIRFOV,stateSpaceDefinLoc,"retrOutputPrimary%def")
    call defGranule(backgroundOutput%def,nIRFOV,stateSpaceDefinLoc,"backgroundOutput%def")

    retrOutputSecondary%nFOR = nFOR
    if (genControl%secondaryRetrOn) then
       stateSpaceDefinLoc=stateSpaceDefinSecondary
    else
       stateSpaceDefinLoc=unitStateSpaceDefin
    endif
    call defGranule(retrOutputSecondary%def,nIRFOV,stateSpaceDefinLoc,"retrOutputSecondary%def")

    linearOutput%nFOR = nFOR
    if (genControl%primaryRetrOn .and. genLinear) then
       call defGranuleLin(linearOutput%def,nChanMW,nChanIR,"linearOutput%def")
    else
       linearOutput%def%nParG = 1
       linearOutput%def%nChan = 1
    endif

    call granAlloc1D(retrOutputStat,aName="retrOutputStat")
    if (MWRTcontrol%radOrTb == 0) then
       retrOutputMW%def%roUnits = TB_UNITS
    else
       retrOutputMW%def%roUnits = RADMW_UNITS
    endif
    retrOutputPrimary%def%roUnits = RAD_UNITS
    backgroundOutput%def%roUnits = RAD_UNITS
    call granAlloc1D(retrOutputMW,aName="retrOutputMW")
    call granAlloc2D(retrOutputPrimary,aName="retrOutputPrimary")
    call granAlloc1D(backgroundOutput,aName="backgroundOutput")
    call granAllocLin(linearOutput,aName="linearOutput")
    call granAlloc2D(retrOutputSecondary,aName="retrOutputSecondary")
    if (dbg) print *, trim(procName)//" ending ..." 
  end subroutine granAlloc

  subroutine defGranule(aGranDefin,nFOV,stateSpaceDefin,aName)
    TYPE(OutputDefin_t), intent(inout) :: aGranDefin
    integer, intent(in) :: nFOV
    TYPE(StateSpaceDefin_t), intent(in) :: stateSpaceDefin
    character (len=*), intent(in) :: aName
    !Local
    integer :: iMW
    character (len=mxLength) :: procName

    procName = "[GranuleMemModule::defGranule]:"
    if (dbg) print *, trim(procName)//" starting ..." 
    if (dbg) print *, "...defining "//trim(aName)//" ..."
    iMW = index(aName,"MW")
    if (iMW == 0) then !under MW, these are undefined
       aGranDefin%nEmRfIR = stateSpaceDefin%NG%emRfIR
       aGranDefin%nFOV = nFOV
    end if
    aGranDefin%NG = stateSpaceDefin%NG
    aGranDefin%IG = genIndices(stateSpaceDefin%NG)
    aGranDefin%nEmMW = stateSpaceDefin%NG%emMW
    aGranDefin%nParG = getVectorLength(stateSpaceDefin%NG)
    aGranDefin%nLev = stateSpaceDefin%nLev
    if (dbg) print *, trim(procName)//" ending ..." 
  end subroutine defGranule

  subroutine defGranuleLin(aGranDefin,nChanMW,nChanIR,aName)
    USE LinInvert, Only: getLinInvData
    TYPE(LinearDefin_t), intent(inout) :: aGranDefin
    integer, intent(in) :: nChanMW
    integer, intent(in) :: nChanIR
    character (len=*), intent(in) :: aName
    !Local
    TYPE(StateIndex_t) :: NGout
    TYPE(StateIndex_t) :: IGout
    integer :: nParGon
    character (len=mxLength) :: procName

    procName = "[GranuleMemModule::defGranuleLin]:"
    if (dbg) print *, trim(procName)//" starting ..." 
    if (dbg) print *, "...defining "//trim(aName)//" ..."
    call getLinInvData(NGout,IGout,nParGon)
    aGranDefin%nParG = nParGon
    aGranDefin%nChan = nChanMW+nChanIR
    aGranDefin%NG = NGout
    aGranDefin%IG = IGout

    if (dbg) print *, trim(procName)//" ending ..." 
  end subroutine defGranuleLin

  subroutine granAllocLin(aGranule,aName)
    TYPE(LinearGran_t), intent(inout) :: aGranule
    character (len=*), intent(in) :: aName
    !Local
    integer :: iFOR
    character (len=mxLength) :: procName

    procName = "[GranuleMemModule::granAllocLin]:"
    if (dbg) print *, trim(procName)//" starting ..." 
    if (dbg) print *, "...allocating memory for "//trim(aName)//" ..."
    if (.not. allocated(aGranule%FOR)) allocate(aGranule%FOR(aGranule%nFOR))
    do iFOR = 1, aGranule%nFOR
       if (.not. allocated(aGranule%FOR(iFOR)%xOffG)) &
            allocate(aGranule%FOR(iFOR)%xOffG(aGranule%def%nParG))
       if (.not. allocated(aGranule%FOR(iFOR)%xGainG)) &
            allocate(aGranule%FOR(iFOR)%xGainG(aGranule%def%nParG,aGranule%def%nChan))
    end do
    if (dbg) print *, trim(procName)//" ending ..." 
  end subroutine granAllocLin

  subroutine granAlloc1D(aGranule,aName)
    TYPE(OutputGran_t), intent(inout) :: aGranule
    character (len=*), intent(in) :: aName
    !Local
    integer :: iFOR
    integer :: iMW
    integer :: nEmisIR
    character (len=mxLength) :: procName

    procName = "[GranuleMemModule::granAlloc1D]:"
    if (dbg) print *, trim(procName)//" starting ..." 
    if (dbg) print *, "...allocating memory for "//trim(aName)//" ..."
    if (.not. allocated(aGranule%FOR)) allocate(aGranule%FOR(aGranule%nFOR))
    iMW = index(aName,"MW")  !could use toLowerCase(aName) when available
    do iFOR = 1, aGranule%nFOR
       nEmisIR = 0.5*aGranule%def%nEmRfIR
       if (.not. allocated(aGranule%FOR(iFOR)%state)) &
            allocate(aGranule%FOR(iFOR)%state(aGranule%def%nParG))
       if (.not. allocated(aGranule%FOR(iFOR)%press)) &
            allocate(aGranule%FOR(iFOR)%press(aGranule%def%nLev))
       if (.not. allocated(aGranule%FOR(iFOR)%emMW)) &
            allocate(aGranule%FOR(iFOR)%emMW(aGranule%def%nEmMW))
       if (iMW >= 1) cycle
       if (.not. allocated(aGranule%FOR(iFOR)%emisIR)) &
            allocate(aGranule%FOR(iFOR)%emisIR(nEmisIR))
    end do
    if (dbg) print *, trim(procName)//" ending ..." 
  end subroutine granAlloc1D

  subroutine granAlloc2D(aGranule,aName)
    TYPE(OutputGran_t), intent(inout) :: aGranule
    character (len=*), intent(in) :: aName
    !Local
    integer :: iFOV
    integer :: iFOR
    integer :: iMW
    integer :: nEmisIR
    character (len=mxLength) :: procName

    procName = "[GranuleMemModule::granAlloc1D]:"
    if (dbg) print *, trim(procName)//" starting ..." 
    if (dbg) print *, "...allocating memory for "//trim(aName)//" ..."
    if (.not. allocated(aGranule%FOR2D)) &
         allocate(aGranule%FOR2D(aGranule%def%nFOV,aGranule%nFOR))
    iMW = index(aName,"MW")  !could use toLowerCase(aName) when available
    do iFOR = 1, aGranule%nFOR
       do iFOV = 1, aGranule%def%nFOV
          nEmisIR = 0.5*aGranule%def%nEmRfIR
          if (.not. allocated(aGranule%FOR2D(iFOV,iFOR)%state)) &
               allocate(aGranule%FOR2D(iFOV,iFOR)%state(aGranule%def%nParG))
          if (.not. allocated(aGranule%FOR2D(iFOV,iFOR)%press)) &
               allocate(aGranule%FOR2D(iFOV,iFOR)%press(aGranule%def%nLev))
          if (.not. allocated(aGranule%FOR2D(iFOV,iFOR)%emMW)) &
               allocate(aGranule%FOR2D(iFOV,iFOR)%emMW(aGranule%def%nEmMW))
          if (iMW >= 1) cycle
          if (.not. allocated(aGranule%FOR2D(iFOV,iFOR)%emisIR)) &
               allocate(aGranule%FOR2D(iFOV,iFOR)%emisIR(nEmisIR))
       end do
    end do
    if (dbg) print *, trim(procName)//" ending ..." 
  end subroutine granAlloc2D

  subroutine granDestroy(&
       retrOutputStat, &
       retrOutputMW, &
       retrOutputPrimary, &
       linearOutput, &
       backgroundOutput, &
       retrOutputSecondary &
       )
    TYPE(OutputGran_t), intent(inout) :: retrOutputStat
    TYPE(OutputGran_t), intent(inout) :: retrOutputMW
    TYPE(OutputGran_t), intent(inout) :: retrOutputPrimary
    TYPE(OutputGran_t), intent(inout) :: backgroundOutput
    TYPE(LinearGran_t), intent(inout) :: linearOutput
    TYPE(OutputGran_t), intent(inout) :: retrOutputSecondary
    !Local
    character (len=mxLength) :: procName

    procName = "[GranuleMemModule::granDestroy]:"
    if (dbg) print *, trim(procName)//" starting..."
    print *, trim(procName)//" starting..."
    if (allocated(retrOutputStat%FOR)) deallocate(retrOutputStat%FOR)
    print *,'retrOutputStat%FOR'
    if (allocated(retrOutputStat%FOR2D)) deallocate(retrOutputStat%FOR2D)
    print *,'retrOutputStat%FOR2D'
    if (allocated(retrOutputMW%FOR)) deallocate(retrOutputMW%FOR)
    print *,'retrOutputMW%FOR'
    if (allocated(retrOutputPrimary%FOR)) deallocate(retrOutputPrimary%FOR)
    print *,'retrOutputPrimary%FOR'
    if (allocated(retrOutputPrimary%FOR2D)) deallocate(retrOutputPrimary%FOR2D)
    print *,'retrOutputPrimary%FOR2D'
    if (allocated(linearOutput%FOR)) deallocate(linearOutput%FOR)
    print *,'linearOutput%FOR'
    if (allocated(backgroundOutput%FOR)) deallocate(backgroundOutput%FOR)
    print *,'backgroundOutput%FOR'
    if (allocated(backgroundOutput%FOR2D)) deallocate(backgroundOutput%FOR2D)
    print *,'backgroundOutput%FOR2D'
    if (allocated(retrOutputSecondary%FOR)) deallocate(retrOutputSecondary%FOR)
    print *,'retrdOutputSecondary%FOR'
    if (allocated(retrOutputSecondary%FOR2D)) deallocate(retrOutputSecondary%FOR2D)
    print *,'retrdOutputSecondary%FOR2D'
    if (dbg) print *, trim(procName)//" ending..."
    print *, trim(procName)//" ending..."
  end subroutine granDestroy

End Module GranuleMemModule
