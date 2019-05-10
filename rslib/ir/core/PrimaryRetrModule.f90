MODULE PrimaryRetrModule
  !
  ! Module for executing MW retrieval process
  !
  ! Subroutine:
  !
  !     primaryRetrInit
  !     primaryRetr
  !     primaryRetrDestroy
  !
  ! Derived data type:
  !
  !     PrimaryRetrControl_t
  !
  ! USE:
  !
  !     ControlStructure
  !     StateIndexModule
  !     LinInvert
  !     InvertModule
  !     LvlInterp
  !     IRMapInvert
  !     MapInvert
  !     oss_ir_module
  !     oss_mw_module
  !     FiniteDiffModule
  !     CHKValid
  !     MWRTmodule
  !     IRRTmodule
  !     MWobsStructure
  !     IRobsStructure
  !     BackgroundModule
  !
  ! yhe@aer.com, 03/16/2016
  !
  USE constants, Only: &
       MISSING_REAL

  USE ToolboxModule, Only: &
       getUnit

  USE ControlStructure, Only: &
       GenControl_t, &
       StateSpaceDefin_t, &
       RetrResult_t, &
       LinearResult_t

  USE StateIndexModule, Only: &
       initLengths, &
       genIndices, &
       StateIndex_t, &
       getVectorLength, &
       getAtmosVectorLength, &
       getNmol, &
       maxMol, &
       whereH2O, &
       nullMolID

  USE LinInvert, Only: &
       initBuildLin, &
       buildLinearCoef

  USE InvertModule, Only: &
       dminv, &
       invrt1, &
       DBL_INVT

  USE LvlInterp, Only: &
       lvl_int

  USE IRMapInvert, Only: &
       IRmap_retr2geo, &
       set_IRMW_invert

  USE MapInvert, Only: &
       map_geo2retr

  USE oss_ir_module, Only: &
       ossdrv_ir

  USE oss_mw_module, Only: &
       ossdrv_mw

  USE FiniteDiffModule, Only: &
       initFiniteDiff, &
       loadFiniteDiff, &
       finiteDiff, &
       destroyFiniteDiff

  USE CHKValid, Only: &
       getChiSq, &
       chkges

  USE scene_io_module, ONLY: &
       TB_UNITS, &
       RADMW_UNITS

  USE MWRTmodule, Only: &
       MWRTcontrol_t

  USE IRRTmodule, Only: &
       IRRTcontrol_t

  USE MWobsStructure, Only: &
       MWdefin_t, &
       MWob_t

  USE IRobsStructure, Only: &
       IRdefin_t, &
       IRob_t

  USE MeasErrStructure, Only: &
       MeasErr_t

  USE ChannelSelectionModule, Only: &
       ChannelSelections_t

  USE AncillaryStructure, Only: &
       AncillaryDefin_t, &
       AncillaryData_t

  USE VertCoord, ONLY: &
       mxCoordTyp, &
       Pcoord

  USE BackgroundModule

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: &
       primaryRetrInit, &
       primaryRetr, &
       primaryRetrDestroy, &
       PrimaryRetrControl_t

  integer, parameter :: mxiter=15

  TYPE PrimaryRetrControl_t
     TYPE(StateIndex_t) :: nRetr
     integer :: maxIter
     integer, dimension(maxMol) :: nGas
     real    :: drad
     real    :: chiSqConv
     real    :: chiSqRatioConv
     logical :: genLinInvert
     logical :: MWon
  END TYPE PrimaryRetrControl_t

  TYPE(StateIndex_t) :: NR_orig
  TYPE(StateIndex_t) :: IR_orig
  TYPE(StateIndex_t) :: NR
  TYPE(StateIndex_t) :: IR
  TYPE(StateIndex_t) :: IGRout
  TYPE(StateIndex_t) :: NGRout
  TYPE(StateIndex_t) :: IG
  TYPE(StateIndex_t) :: NG
  TYPE(background_t) :: bkgObj
  integer :: nChanIR
  integer :: nChanMW
  integer :: nChan
  integer :: ncldLiqOrig
  integer :: ncldIceOrig
  integer :: nEmIR_orig
  integer :: nParAtm
  integer :: nParSfcMW
  integer :: nParMax
  integer :: nPar
  integer :: iH2O
  integer :: nParG
  integer :: nParGAtm
  integer :: nParGSfcMW
  integer :: nParGSfcIR
  integer :: nParGon
  integer :: nLev
  integer :: IEmMwG
  integer :: nMol
  integer :: IEmIRG
  integer :: nEmIRG
  integer, dimension(:), allocatable :: MolTran
  logical, dimension(:), allocatable :: kchan
  real,    dimension(:), allocatable :: pRef
  real, dimension(:), allocatable :: press
  real, dimension(:), allocatable :: rad
  real, dimension(:), allocatable :: radObs
  real, dimension(:), allocatable :: Sy,delY,xkEmMw,emMW
  real, dimension(:), allocatable :: Ylast
  real, dimension(:), allocatable :: xGes,xGesG,xGesNew
  real, dimension(:), allocatable :: xGesLast
  real, dimension(:), allocatable :: background
  real, dimension(:), allocatable :: xBak
  real, dimension(:), allocatable :: xBakSav
  real, dimension(:), allocatable :: EmRf
  real, dimension(:), allocatable :: EmRfGrid
  real, dimension(:,:), allocatable :: EmRfdrv

  real, dimension(:), allocatable :: dradChan
  real, dimension(:), allocatable :: xOffG
  real, dimension(:), allocatable :: rerr
  real, dimension(:), allocatable :: atmNoise
  real, dimension(:), allocatable :: dummyArr
  real, dimension(:,:), allocatable :: xktG,xkt
  real, dimension(:,:), allocatable :: clw_cov,ciw_cov
  real, dimension(:,:), allocatable :: ut_cov_net,covRetr,hm1,aPostErrorCovG
  real, dimension(:), allocatable :: diagError
  real, dimension(:,:), allocatable :: umtx,vmtx,wk
  real, dimension(:,:), allocatable :: xGainG
  real, dimension(:,:,:), allocatable :: xkEmRf
  real, dimension(:,:,:), allocatable :: xkEmRfdrv
  real(kind=DBL_INVT), dimension(:,:), allocatable :: dcovRetr
  character (len=mxCoordTyp) :: vCoordTyp
  logical :: extAtmosOn
  logical :: logSize
  logical :: debug = .false.
  logical :: dbg = .false.

CONTAINS

  subroutine primaryRetrInit( &
       genControl, &
       primaryRetrControl, &
       MWRTcontrol, &
       IRRTcontrol, &
       MWdefin, &
       IRdefin, &
       ancillaryDefin, &
       stateSpaceDefin &
       )
    TYPE(GenControl_t), intent(in) :: genControl
    TYPE(PrimaryRetrControl_t), intent(in) :: primaryRetrControl
    TYPE(MWRTcontrol_t), intent(in) :: MWRTcontrol
    TYPE(IRRTcontrol_t), intent(in) :: IRRTcontrol
    TYPE(MWdefin_t), intent(inout) :: MWdefin
    TYPE(IRdefin_t), intent(inout) :: IRdefin
    TYPE(AncillaryDefin_t), intent(in) :: ancillaryDefin
    TYPE(StateSpaceDefin_t), intent(inout) :: stateSpaceDefin

    !Local
    integer :: iFOR
    integer :: nProfs
    character (len=*), parameter :: procName=&
                                       ' [PrimaryRetrModule::primaryRetrInit]: '

    if (dbg) print *, trim(procName)//'starting ...'

    call bkgObj%bkgInitDefault()

    call bkgObj%bkgInitAtmos( &
         stateSpaceDefin%nParGatm, &
         stateSpaceDefin%IG, &
         stateSpaceDefin%NG, &
         stateSpaceDefin%nLev, &
         stateSpaceDefin%vCoordTyp, &
         stateSpaceDefin%pRef, &
         stateSpaceDefin%MolID, &
         stateSpaceDefin%MolTran, &
         genControl%limitsConfig, &
         genControl%backgroundFiles(1))
    ! This handling of logSize is temporary until background files store it
    stateSpaceDefin%logSize=genControl%logSize

    if (primaryRetrControl%MWon) then
       call bkgObj%bkgInitSfcMW( &
            stateSpaceDefin%nChanMW, &
            stateSpaceDefin%frqMW, &
            stateSpaceDefin%polMW, &
            stateSpaceDefin%IG%emMW, &
            genControl%backgroundFiles(2))
    else
       stateSpaceDefin%nChanMW = 0
    endif

    call bkgObj%bkgInitSfcIR( &
         stateSpaceDefin%nEmIR, &
         stateSpaceDefin%wvnIR, &
         stateSpaceDefin%IG%emRfIR, &
         genControl%backgroundFiles(3))

    stateSpaceDefin%NG%emMW = stateSpaceDefin%nChanMW
    stateSpaceDefin%NG%emRfIR = 2*stateSpaceDefin%nEmIR
    stateSpaceDefin%IG = genIndices(stateSpaceDefin%NG)
    stateSpaceDefin%nParG = getVectorLength(stateSpaceDefin%NG)

    if (primaryRetrControl%maxiter > mxiter) then
       print *, trim(procName)//'Error: primaryRetrControl%maxiter > mxiter ', &
            primaryRetrControl%maxiter, mxiter
       call exit(1)
    endif

    !----Set retrieval-space pointers
    call setXPtrR( &
         primaryRetrControl%nRetr%Tskin, &
         primaryRetrControl%nRetr%temp,  &
         primaryRetrControl%nRetr%mol, &      !nGas
         primaryRetrControl%nRetr%cldLiq, &
         primaryRetrControl%nRetr%cldIce, &
         IRRTcontrol%molID, &
         IRRTcontrol%nMol, &
         stateSpaceDefin%NG, &
         IR_orig, NR_orig)

    nParAtm = getVectorLength(NR_orig)
    NR_orig%emMW = primaryRetrControl%nRetr%emMW
    if (.not. primaryRetrControl%MWon) NR_orig%emMW = 0
    nParSfcMW = NR_orig%emMW
    if (primaryRetrControl%nRetr%EmRfIR == 0) then
       NR_orig%emRfIR = 0
    else
       NR_orig%emRfIR = 2*stateSpaceDefin%nEmIR
    endif
    IR_orig = genIndices(NR_orig)

    !----Verify that the selected number of variables to retrieve are
    !----compatible with the transformation matrices
    call bkgObj%checkTransformConform(NR_orig)

    nParGsfcMW = stateSpaceDefin%nChanMW
    nParGsfcIR = stateSpaceDefin%NG%emRfIR

    !--- Save the original retrieval-space pointers
    !----------------------------------------------------------
    !  Save original versions of variables that are dynamically
    !  changed in the ProfLoop via call to setRetrVec()
    !  to fall back if needed.
    !----------------------------------------------------------
    ncldLiqOrig = primaryRetrControl%nRetr%cldLiq
    ncldIceOrig = primaryRetrControl%nRetr%cldIce
    NG = stateSpaceDefin%NG
    IG = genIndices(NG)
    nParMax = getVectorLength(NR_orig)
    nParG = stateSpaceDefin%nParG !=retrPrior%nParGmax from RetrModule::init
    call bkgObj%bkgAllocate(nParG)
    nParGAtm = stateSpaceDefin%nParGAtm
    nLev = stateSpaceDefin%nLev
    IEmMwG = stateSpaceDefin%IG%emMW
    IEmIRG = stateSpaceDefin%IG%emRfIR
    nEmIrG = stateSpaceDefin%nEmIR
    iH2O = whereH2O(stateSpaceDefin%MolID)
    nMol = getNmol(stateSpaceDefin%MolID)
    logSize = stateSpaceDefin%logSize
    if (primaryRetrControl%MWon) then
       nChanMW = MWdefin%nChan
    else
       nChanMW = 0
    endif
    nChanIR = IRdefin%nChan
    nChan = nChanMW+nChanIR
    nProfs = genControl%nFORmax
    vCoordTyp = trim(stateSpaceDefin%vCoordTyp)

    if (.not. allocated(pRef)) allocate(pRef(nLev))
    if (.not. allocated(press)) allocate(press(nLev))
    if (.not. allocated(MolTran)) allocate(MolTran(size(stateSpaceDefin%MolTran)))
    if (.not. allocated(EmRfGrid))  allocate(EmRfGrid (nEmIrG))
    if (.not. allocated(EmRfdrv))   allocate(EmRfdrv  (nEmIrG,2))
    if (.not. allocated(xkEmRfdrv)) allocate(xkEmRfdrv(nEmIrG,nChanIR,2))
    if (.not. allocated(kchan)) allocate(kchan(nChan))
    if (.not. allocated(rad)) allocate(rad(nChan))
    if (.not. allocated(radObs)) allocate(radObs(nChan))
    if (.not. allocated(Sy)) allocate(Sy(nChan))
    if (.not. allocated(delY)) allocate(delY(nChan))
    if (.not. allocated(Ylast)) allocate(Ylast(nChan))
    if (.not. allocated(xGesLast)) allocate(xGesLast(nParMax))
    if (.not. allocated(rerr)) allocate(rerr(nChan))
    if (.not. allocated(atmNoise)) allocate(atmNoise(nChan))
    if (.not. allocated(dummyArr)) allocate(dummyArr(nChan))
    if (.not. allocated(EmRf)) allocate(EmRf(2*nEmIrG))
    if (.not. allocated(xkEmRf)) allocate(xkEmRf(2,nEmIrG,nChanIR))
    if (.not. allocated(xBakSav)) allocate(xBakSav(nParG))
    if (.not. allocated(xGesNew)) allocate(xGesNew(nParMax))
    if (.not. allocated(xktG)) allocate(xktG(nParG,nChan))
    if (.not. allocated(xkt)) allocate(xkt(nParMax,nChan))
    if (.not. allocated(background)) allocate(background(nParG))
    if (.not. allocated(clw_cov)) allocate(clw_cov(nCldLiqOrig,nCldLiqOrig))
    if (.not. allocated(ciw_cov)) allocate(ciw_cov(nCldIceOrig,nCldIceOrig))
    if (.not. allocated(ut_cov_net)) allocate(ut_cov_net(nParG,nParG))
    if (.not. allocated(hm1)) allocate(hm1(nParMax,nParMax))
    if (.not. allocated(wk)) allocate(wk(nParg,nParMax))
    if (.not. allocated(covRetr)) allocate(covRetr(nParMax,nParMax))
    if (.not. allocated(aPostErrorCovG)) allocate(aPostErrorCovG(nParG,nParG))
    if (.not. allocated(diagError)) allocate(diagError(nParG))
    if (.not. allocated(dcovRetr)) allocate(dcovRetr(nParMax,nParMax))
    if (.not. allocated(umtx)) allocate(umtx(nParMax,nParG))
    if (.not. allocated(vmtx)) allocate(vmtx(nParG,nParMax))
    if (.not. allocated(xGesG)) allocate(xGesG(nParG))
    if (.not. allocated(xGes)) allocate(xGes(nParMax))
    if (.not. allocated(xkEmMw)) allocate(xkEmMw(max(nChanMW,1)))
    if (.not. allocated(emMW)) allocate(emMW(max(nChanMW,1)))
    EmRfGrid = stateSpaceDefin%wvnIR

    !initialization of kchan, necessary to avoid non-zero values for
    !tripping some logic conditions in IRmap_Invert::set_irmw_invert,
    !which would cause segmentation fault due to accesing to invalid
    !memory
    kchan = .false.

    pRef = stateSpaceDefin%pRef
    MolTran = stateSpaceDefin%MolTran

    ! Initialize linear coefficient building: should be modified per
    ! introduction of the hierarchical definitions for the linear
    ! retrieval output
    IF (primaryRetrControl%genLinInvert) then
       call initBuildLin(stateSpaceDefin%IG,stateSpaceDefin%NG, &
            stateSpaceDefin%nParG,NR_orig,IGRout,NGRout,nParGon)
       if (.not. allocated(xOffG)) allocate (xOffG(nParGon))
       if (.not. allocated(xGainG)) allocate(xGainG(nParGon,nChan))
    ENDIF

    !--------------------------------------------------------
    ! Allocate and initialize dradChan with drad
    !--------------------------------------------------------
    IF (primaryRetrControl%MWon) THEN
       if (.not. allocated(dradChan)) allocate(dradChan(nChanMW))
       ! In radiance model, drad must be converted from dTB to dradiance
       IF (MWdefin%radOrTb == 0) THEN
          dradChan=primaryRetrControl%drad
       ELSE
          dradChan(1:nChanMW)=primaryRetrControl%drad/MWdefin%planckAlpha(1:nChanMW)
       ENDIF
    ELSE
       if (.not. allocated(dradChan)) allocate(dradChan(1))
    ENDIF

    call bkgObj%validateExtAtmos(&
         AncillaryDefin%NG, &
         nMol, &
         IRRTcontrol%MolID, &
         ancillaryDefin%MolID, &
         genControl%externDataFlags%temp, &
         genControl%externDataFlags%mol(1:nMol), &
         genControl%externDataFlags%pSfc, &
         genControl%externDataFlags%wind, &
         genControl%externDataFlags%cloudLiq, &
         genControl%externDataFlags%cloudIce)

    if (genControl%externDataFlags%emRfIR) &
         call bkgObj%validateExtSfcIR(ancillaryDefin%nEmIR,ancillaryDefin%wvn)

    extAtmosOn=.FALSE.
    if (ANY((/genControl%externDataFlags%temp, &
         genControl%externDataFlags%mol, &
         genControl%externDataFlags%pSfc, &
         genControl%externDataFlags%wind, &
         genControl%externDataFlags%cloudLiq, &
         genControl%externDataFlags%cloudIce/))) extAtmosOn=.TRUE.

    !- This is temporary, until scattering model with Jacobians is
    !- integrated ... this block needs to move down, to come after
    !- NR_orig has been defined.
    if (genControl%enableScatter .AND. &
         ANY((/NR_orig%temp,NR_orig%mol(1:nMol), &
         NR_orig%emMW,NR_orig%emRfIR/) /= 0)) then
       print*,trim(procName),' err - Jacobians not available for all',&
            ' retrieved variables in scattering mode; '
       print*,'See config file entries for temp, Tskin, gas, emis'
       call exit(1)
    endif

    if (genControl%enableScatter) then
       ! nchan = MWdefin%nChan+IRdefin%nChan
       call loadFiniteDiff(genControl%backgroundConfig)
       !pass in individual fields instead of the structure is for
       !backward compatibility in retr_sig.f90
       call initFiniteDiff(stateSpaceDefin%IG,stateSpaceDefin%NG,&
            stateSpaceDefin%vCoordTyp,nChan)
    endif

    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine primaryRetrInit

  subroutine primaryRetr(&
       genControl, &
       primaryRetrControl, &
       MWob, &
       MWmeasErr, &
       IRdefin, &
       IRob, &
       IRmeasErr, &
       ancillaryDefin, &
       ancillaryOb, &
       channelSelections, &
       atmosClass, &
       cloudClass, &
       surfaceClassMW, &
       surfaceClassIR, &
       isLandMW, &
       isLandIR, &
       retrPrior, &
       retrResult, &
       linearResult, &
       backResult) !rename xGesMwSave; vector, not structure

    use IRRTmodule, only: &
          IRRTsetActiveChannelSelectionSet

    TYPE(GenControl_t),         intent(in)    :: genControl
    TYPE(PrimaryRetrControl_t), intent(in)    :: primaryRetrControl
    TYPE(MWOb_t),               intent(in)    :: MWob
    TYPE(MeasErr_t)                           :: MWmeasErr
    TYPE(IRdefin_t),            intent(in)    :: IRdefin
    TYPE(IRob_t),               intent(in)    :: IRob
    TYPE(MeasErr_t)                           :: IRmeasErr
    TYPE(AncillaryDefin_t),     intent(in)    :: ancillaryDefin
    TYPE(AncillaryData_t),      intent(in)    :: ancillaryOb
    TYPE(ChannelSelections_t),  intent(in)    :: ChannelSelections
    TYPE(RetrResult_t),         intent(in)    :: retrPrior
    integer,                    intent(in)    :: atmosClass
    integer, dimension(:),      intent(in)    :: cloudClass
    integer,                    intent(in)    :: surfaceClassMW
    integer,                    intent(in)    :: surfaceClassIR
    logical,                    intent(in)    :: isLandMW
    logical,                    intent(in)    :: isLandIR
    TYPE(RetrResult_t),         intent(inout) :: retrResult
    TYPE(LinearResult_t),       intent(inout) :: linearResult
    real, dimension(:),         intent(inout) :: backResult

    !Local
    integer :: nch
    integer :: nChanOn
    integer :: iFOV
    integer :: iter
    integer :: ifail
    integer :: nParPrior
    integer :: nParGverif
    integer :: U_algconfig
    real :: xSfc
    real :: Tsfc
    real :: TsfcLast
    real :: chiSq0
    real :: chiSqIR
    real :: chiSqMW
    real :: chiSqRatio
    real :: rmsIR
    real :: rmsMW
    real(kind=DBL_INVT) :: detmnt
    real, dimension(mxiter) :: rms_IR,chiSq_IR
    logical :: xOutRangeFlag
    logical :: doVerif
    TYPE(retFlags_t) :: retrievalFlags = RET_FLAGS_ALL_ON
    real    :: azAngle,DOFS
    integer :: jj, iChan, kk,ip

    logical :: chanDataCompressed = .true. ! optional mode for getChiSq stating that
!                         the computed radaince is compressed.

    ! obsLevel the level at which the instrumnet is located
    ! TOA corresponds to obsLevel=1
    integer :: obsLevel=1

    character (len=*), parameter :: procName = &
                                    ' [PrimaryRetrModule::primaryRetr]: '
    if (dbg) print *, procName, 'starting ...'
    ! Disable retrieval of cloud parameters per cloud mask
    ! assumption: use the first element of cloudClass for the cloudModeTest
    iFOV = 1
    if (genControl%imgDataOn) &
         call bkgObj%setCloudMode(cloudClass(iFOV),retrievalFlags,&
         genControl%imgDataOn)

    ! Retrieve parameters based on the values of the flags in retrievalFlags
    ! Currently it is all ON.  Need to change that.
    CALL bkgObj%setRetrVec(retrievalFlags,NR_orig,nPar,NR,IR)

    umtx = 0.0
    vmtx = 0.0
    ut_cov_net = 0.0
    call bkgObj%setBkgClimAtmos(atmosClass,ancillaryOb%pSfc,vcoordTyp,&
         press,umtx(1:nParAtm,1:nParGAtm),vmtx(1:nParGatm,1:nParAtm))
    if (primaryRetrControl%MWon) &
       call bkgObj%setBkgClimSfcMW(surfaceClassMW, &
            umtx(nParAtm+1:nParAtm+nParSfcMW,nParGatm+1:nParGatm+nParGSfcMW), &
            vmtx(nParGatm+1:nParGatm+nParGSfcMW,nParAtm+1:nParAtm+nParSfcMW))
    call bkgObj%setBkgClimSfcIR(surfaceClassIR, &
         umtx(nParAtm+nParSfcMW+1:nParMax,nParGatm+nParGsfcMW+1:nParG), &
         vmtx(nParGatm+nParGsfcMW+1:nParG,nParAtm+nParSfcMW+1:nParMax))
    if (extAtmosOn) call bkgObj%setBkgExt(ancillaryOb%stateV)
    if (genControl%externDataFlags%emRfIR) &
         call bkgObj%setBkgSfcIR(ancillaryOb%emRfIR(1:ancillaryDefin%nEmRfIR),&
              ancillaryDefin%nEmIR)
    call bkgObj%combineBkgAtmos(ancillaryDefin%IG,ancillaryDefin%NG,&
         ancillaryOb%pSfc,isLandIR,MolTran,vCoordTyp,press,background(1:nParGatm), &
         ut_cov_net(1:nParGatm,1:nParGatm))
    if (primaryRetrControl%MWon) &
       call bkgObj%combineBkgSfcMW(isLandMW,background(nParGatm+1:nParGatm+nParGSfcMW), &
            ut_cov_net(nParGatm+1:nParGatm+nParGSfcMW,nParGatm+1:nParGatm+nParGSfcMW))
    call bkgObj%combineBkgSfcIR(isLandIR,&
         background(nParGatm+nParGSfcMW+1:nParGatm+nParGSfcMW+nParGsfcIR), &
         ut_cov_net(nParGatm+nParGSfcMW+1:nParGatm+nParGSfcMW+nParGsfcIR,&
                    nParGatm+nParGSfcMW+1:nParGatm+nParGSfcMW+nParGsfcIR))
    call bkgObj%overrideBkg(background,ut_cov_net,isLandIR)
    call bkgObj%transformBkg(ut_cov_net,covRetr(1:nPar,1:nPar), &
         umtx(1:nPar,1:nParG))

    !Concatenate MW and IR entities before calling set_IRMW_Invert
    if (primaryRetrControl%MWon) then
       rerr(1:nChanMW) = MWmeasErr%device
       atmNoise(1:nChanMW) = MWmeasErr%fwdModel
       radObs(1:nChanMW) = MWob%rad
    end if
    rerr(nChanMW+1:nChan) = IRmeasErr%device
    atmNoise(nChanMW+1:nChan) = IRmeasErr%fwdModel
    radObs(nChanMW+1:nChan) = IRob%rad

    !---Save the background for later storage--
    backResult(1:nParG)=background(1:nParG)
    backResult(IG%mol(iH2O):IG%mol(iH2O)+NG%mol(iH2O)-1)=&
         EXP(background(IG%mol(iH2O):IG%mol(iH2O)+NG%mol(iH2O)-1))  !should check molTran

    call lvl_int(backResult,pRef,Nlev,backResult(IG%Psfc),xSfc)
    backResult(IG%Tskin)=background(IG%Tskin)+xSfc

    rms_IR   = MISSING_REAL
    chiSq_IR = MISSING_REAL

    !do this if retrPrior%nParG < stateSpaceDefin%nParG
    if (retrPrior%nParG < nParG) then
       xGes=0.
       if(NR%cldLiq > 0)xGes(IR%cldliq+2)=genControl%LWP1stGuess
       if(NR%cldIce > 0)xGes(IR%cldIce+2)=genControl%IWP1stGuess
       call IRmap_retr2geo(xGes(1:nPar),background,xGesG,emMW,Tsfc, &
            vmtx(1:nParG,1:nPar),nParG,nPar,IG,NG,vCoordTyp, &
            pRef,nLev,nChanMW,IEmMwG,iH2O,nMol,molTran,IEmIRG,nEmIrG, &
            EmRf,xOutRangeFlag,logSize)
       if (xOutRangeFlag) &
            print *, procName, ' warning: x out of range may affect convergence'
    end if

    if (retrPrior%nParG > 0) then
       !Generate the first guess in the retrieval space based on the
       !previous retrieval process
       xGesG(1:retrPrior%nParG) = retrPrior%stateV(1:retrPrior%nParG)

       doVerif = .TRUE.
       nParPrior=IR%emMW-1
       nParGverif=IG%emMW-1
       if (retrPrior%nEmMW > 0 .and. primaryRetrControl%MWon) then
          emMW(1:retrPrior%nEmMW) = retrPrior%emMW(1:retrPrior%nEmMW)
          nParPrior=IR%emMW+NR%emMW-1
          nParGverif=IG%emMW+NG%emMW-1
       elseif (retrPrior%nEmMW > 0 .and. .not. primaryRetrControl%MWon) then
          doVerif = .FALSE.  ! Verification does not work in this case
       endif
       if (retrPrior%nEmRfIR > 0) then
          EmRf(1:retrPrior%nEmRfIR) = retrPrior%emRfIR(1:retrPrior%nEmRfIR)
          if (retrPrior%nEmMW == 0) then
             print *, procName,' Error: IR emis/refl 1st guess from '// &
                                      'prior not allowed without MW emis '
             call exit(1)
          endif
          nParPrior=nPar
          nParGverif=nParG
       end if
       if ((nParGverif /= retrPrior%nParG) .and. doVerif) then
          print *, trim(procName)//'Error: 1st guess from prior stage '// &
               'is not a supported size;'
          print *,'    nParGverif, retrPrior%nParG:',nParGverif,retrPrior%nParG
          call exit(1)
       endif

       CALL lvl_int(xGesG(1:retrPrior%nParG),pref,Nlev,xGesG(IG%Psfc),TSfc)
       CALL map_geo2retr(xGesG(1:retrPrior%nParG), &
            background,Tsfc,xGes(1:nParPrior), &
            umtx(1:nParPrior,1:retrPrior%nParG),retrPrior%nParG, &
            nParPrior,IG,NG,nMol,MolTran,vCoordTyp,logSize)
    end if

    xktG=0.  ! To initialize portions not filled later
    if (primaryRetrControl%MWon) &
      kchan(1:nChanMW) = channelSelections%PrimaryCloudy%MW /= 0

    kchan(nChanMW+1:nChan) = channelSelections%PrimaryClear%IR /= 0
    call IRRTsetActiveChannelSelectionSet(channelSelections%PrimaryClear)

   ! ossdrv_mw has a flag lambertion controling the reflection from the surface
   ! if omitted it set to .false. specular reflection
    azAngle = IRob%EAA-IRob%SolAA
    do iChan=1,nEmIrG
      EmRfdrv(iChan,1) = EmRf(2*iChan-1)
      EmRfdrv(iChan,2) = EmRf(2*iChan)
    end do
    rad = 0.
    if (trim(vCoordTyp) == Pcoord) then
       if (primaryRetrControl%MWon) &
            CALL ossdrv_mw(xGesG,emMW,MWob%EIA,obsLevel,rad,xktG,xkEmMw,kchan, &
                 ancillaryOb%landfrac,IRob%lat,IG%temp,IG%Tskin,IG%Psfc,&
                 IG%mol(1:iH2O), ICldLiqG_in=IG%cldLiq)

       call ossdrv_ir(xGesG(1:nParGAtm),EmRfGrid, EmRfdrv,IRob%EIA,IRob%SolIA,&
            azAngle,obsLevel,y=rad(nChanMW+1:nChan),&
            xkt=xktG(1:nParGAtm,nChanMW+1:nChan),xkEmRf=xkEmRfdrv,lat=IRob%lat,&
            tempIndex=IG%temp,tSkinIndex=IG%Tskin,pSurfIndex=IG%Psfc,&
            varMolIndex=IG%mol(1:nMol))
    else
       if (primaryRetrControl%MWon) &
        CALL ossdrv_mw(xGesG,emMW,MWob%EIA,obsLevel,rad,xktG,xkEmMw,&
            kchan,ancillaryOb%landfrac,IRob%lat, IG%temp,IG%Tskin,IG%Psfc,&
            IG%mol(1:iH2O),ICldLiqG_in=IG%cldLiq, puser=press)

       call ossdrv_ir(xGesG(1:nParGAtm),EmRfGrid, EmRfdrv,IRob%EIA,IRob%SolIA,&
            azAngle,obsLevel,y=rad(nChanMW+1:nChan),&
            xkt=xktG(1:nParGAtm,nChanMW+1:nChan),xkEmRf=xkEmRfdrv,lat=IRob%lat,&
            tempIndex=IG%Temp,tSkinIndex=IG%Tskin,pSurfIndex=IG%Psfc,&
            varMolIndex=IG%mol(1:nMol),puser=press)
       do jj=1,npargatm,10
          print *, jj,xgesg(jj:jj+9)
        end do

    end if

    do jj=1,2
      do iChan=1,nChanIR
        xkEmRf(jj,1:nEmIrG,iChan)=xkEmRfdrv(1:nEmIrG,iChan, jj)
      end do
    end do

    if (genControl%enableScatter) &
         call finiteDiff(xGesG(1:nParGAtm),EmRfGrid,EmRfdrv,IRob%EIA,IRob%SolIA, &
         azAngle,IRob%lat,obsLevel,rad(nChanMW+1:nChan),&
         xktG(1:nParGAtm,nChanMW+1:nChan),xkEmRfdrv, ancillaryOb%landfrac,pRef,NR)


    if (chanDataCompressed) &
      call compressRadJacobian(kchan, nChanMW, nChanIR, rad, xktG, xkEmMw)

    chiSq0 = MISSING_REAL

    dcovRetr=real(covRetr,KIND=DBL_INVT)
    call dminv(dcovRetr(1:nPar,1:nPar),nPar,detmnt)
    covRetr=real(dcovRetr)

    if (NR%cldLiq > 0) then
       clw_cov(1:NR%cldLiq,1:NR%cldLiq) = &
            covRetr(IR%cldLiq:IR%cldLiq+NR%CldLiq-1,IR%cldLiq:IR%cldLiq+NR%CldLiq-1)
    endif
    if (NR%cldIce > 0) then
       ciw_cov(1:NR%cldIce,1:NR%cldIce) = &
            covRetr(IR%cldIce:IR%cldIce+NR%CldIce-1,IR%cldIce:IR%cldIce+NR%CldIce-1)
    endif

    if (dbg) print '('' IR iter.  nchan      CC-RMS         RMS   Q-chi_square IR    chi_sq MW '')'

    !to construct Ym, use IRob and IRdefin%chanNum
    IRIterLoop: do iter=1,primaryRetrControl%maxiter

       print *,'radobs ',radobs(1:10)
       print *,'rad ',rad(1:10)
       print *, size(radobs), size(rad)
       call set_IRMW_Invert(xGesG,radObs,rad,xktG,xkEmMw,rerr,atmNoise, &
            dradChan,xkt,xkEmRf,Sy,delY,nch,iter,nMol,MolTran, &
            kchan,vmtx(1:nParG,1:nPar),nParG,nPar,nChanMW, &
            nChan,IG,NG,IEmMwG,IEmIrG,nEmIrG,vCoordTyp,pRef,dummyArr,logSize, &
            chanDataCompressed)

       
       ! now xktG contains emissivity Jacobians
       call invrt1(xkt(1:nPar,1:nch),Sy(1:nch),covRetr(1:nPar,1:nPar), &
            delY(1:nch),xGes(1:nPar),xGesNew(1:nPar),nPar,nch,hm1(1:nPar,1:nPar),dofs)

       !Convert error covariance matrix from retrieved to geophysical space
       wk(1:nParG,1:nPar) = matmul(vmtx(1:nParG,1:nPar),hm1(1:nPar,1:nPar))
       aPostErrorCovG = matmul(wk,transpose(vmtx(1:nParG,1:nPar)))

       !Extract diagonal of the error covariance matrix and get DOFS
       do ip=1,nParg
          diagError(ip) = aPostErrorCovG(ip,ip)
       enddo

       Ylast=rad
       xGesLast=xGes

       xGes(1:nPar) = xGesNew(1:nPar)

       call IRmap_retr2geo(xGesNew(1:nPar),background,xGesG,emMW,Tsfc, &
            vmtx(1:nParG,1:nPar),nParG,nPar,IG,NG,vCoordTyp,pRef,nLev,nChanMW, &
            IEmMwG,iH2O,nMol,MolTran,IEmIRG,nEmIrG,EmRf,xOutRangeFlag,logSize)
       if (xOutRangeFlag) &
            print *, trim(procName)//' warning: x out of range may affect '// &
            'convergence in main loop'
       TsfcLast=Tsfc

       U_algconfig = getUnit()
       call chkges(U_algconfig,genControl%backgroundConfig, &
            xGesG,ifail,covRetr(1:nPar,1:nPar),clw_cov,ciw_cov,nLev,press, &
            IG,NG,IR,NR,iH2O,debug,invCldCov_in=.true.)

       !- The two MW related calls separated by the call to ossdrv_ir
       !- could be merged into one block - YHE, 07/20/2016
       do iChan=1,nEmIrG
        EmRfdrv(iChan,1) = EmRf(2*iChan-1)
        EmRfdrv(iChan,2) = EmRf(2*iChan)
       end do

       if (trim(vCoordTyp) == Pcoord) then
          if (primaryRetrControl%MWon) &
               CALL ossdrv_mw(xGesG,emMW,MWob%EIA,obsLevel,rad,xktG,xkEmMw,kchan,&
                    ancillaryOb%landfrac,IRob%lat,IG%temp,IG%Tskin,IG%Psfc,&
                    IG%mol(1:iH2O), ICldLiqG_in=IG%cldLiq)

          call ossdrv_ir(xGesG(1:nParGAtm),EmRfGrid, EmRfdrv,IRob%EIA,IRob%SolIA,&
               azAngle,obsLevel,y=rad(nChanMW+1:nChan),&
               xkt=xktG(1:nParGAtm,nChanMW+1:nChan),xkEmRf=xkEmRfdrv,lat=IRob%lat,&
               tempIndex=IG%temp,tSkinIndex=IG%Tskin,pSurfIndex=IG%Psfc,&
               varMolIndex=IG%mol(1:nMol))
       else
          if (primaryRetrControl%MWon) &
          CALL ossdrv_mw(xGesG,emMW,MWob%EIA,obsLevel,rad,xktG,xkEmMw,kchan,&
               ancillaryOb%landfrac,IRob%lat,IG%temp,IG%Tskin,IG%Psfc,&
               IG%mol(1:iH2O),ICldLiqG_in=IG%cldLiq,puser=press)

         call ossdrv_ir(xGesG(1:nParGAtm),EmRfGrid,EmRfdrv,IRob%EIA,IRob%SolIA,&
              azAngle,obsLevel,y=rad(nChanMW+1:nChan),&
              xkt=xktG(1:nParGAtm,nChanMW+1:nChan),xkEmRf=xkEmRfdrv,lat=IRob%lat,&
              tempIndex=IG%Temp,tSkinIndex=IG%Tskin,pSurfIndex=IG%Psfc, &
              varMolIndex=IG%mol(1:nMol),puser=press)
       end if

       do jj=1,2
         do iChan=1,nChanIR
            xkEmRf(jj,1:nEmIrG,iChan)=xkEmRfdrv(1:nEmIrG,iChan, jj)
         end do
       end do

       if (genControl%enableScatter) &
            call finiteDiff(xGesG(1:nParGatm),EmRfGrid,EmRfdrv,IRob%EIA,IRob%SolIA, &
            azAngle,IRob%lat,obsLevel,rad(nChanMW+1:nChan),&
            xktG(1:nParGatm,nChanMW+1:nChan),xkEmRf,ancillaryOb%landfrac,press,NR)

       if (chanDataCompressed) &
         call compressRadJacobian(kchan, nChanMW, nChanIR, rad, xktG, xkEmMw)

       if (primaryRetrControl%MWon) then
          call getChiSq(radObs,rad,rerr,chiSqMW,rmsMW,kchan,nChanMW,nChanOn,&
               chanDataCompressed=chanDataCompressed)
       else
          chiSqMW=MISSING_REAL
       endif
       call getChiSq(radObs,rad,rerr,chiSqIR,rmsIR,kchan,nChan,nChanOn,&
            chanDataCompressed=chanDataCompressed)

       rms_IR(iter)=rmsIR
       chiSq_IR(iter)=chisqIR

       if (dbg) print '(i9,i7,2f12.7,2f15.3)',iter,nChanOn,IRmeasErr%rms,rmsIR,&
                                              chiSqIR,chiSqMW

       chiSqRatio=abs((ChiSq0-chiSqIR)/chiSq0)
       chiSq0=chiSqIR

       if (chiSqIR < primaryRetrControl%chisqconv .and. &
            chiSqRatio < primaryRetrControl%chisqRatioconv) exit IRIterLoop

       CALL map_geo2retr(xGesG,background,Tsfc,xGes(1:nPar), &
            umtx(1:nPar,1:nParG),nParG,nPar,IG,NG,nMol,MolTran,vCoordTyp,logSize)

    enddo IRIterLoop

    iter=min(iter,primaryRetrControl%maxiter)
    retrResult%nParG = nParG
    retrResult%nEmRfIR = 2*nEmIRG
    retrResult%stateV(1:nParG) = xGesG(1:nParG)
    if (primaryRetrControl%MWon) then
       retrResult%nEmMW = nChanMW
       retrResult%emMW(1:retrResult%nEmMW) = emMW(1:retrResult%nEmMW)
    end if
    retrResult%diagError = diagError(1:NparG)
    retrResult%DOFS = DOFS
    retrResult%emRfIR(1:retrResult%nEmRfIR) = EmRf(1:retrResult%nEmRfIR)
    retrResult%nIter = iter
    retrResult%rms = rms_IR(iter)
    retrResult%chiSq = chiSq_IR(iter)
    retrResult%IG = IG
    retrResult%NG = NG
    retrResult%press = press

    if (primaryRetrControl%genLinInvert) then
       call buildLinearCoef(xkt(1:nPar,1:nch),Sy(1:nch),covRetr(1:nPar,1:nPar), &
            Ylast,kchan,xGesLast(1:nPar),vmtx(1:nParG,1:nPar), &
            nPar,nch,nParGon,background,xOffG,xGainG(:,1:nch))
       linearResult%xOffG(1:nParGon)=xOffG(1:nParGon)
       linearResult%xGainG(1:nParGon,1:nch)=xGainG(1:nParGon,1:nch)
       linearResult%chiSq=chiSqIR
       linearResult%Tsfc=TsfcLast
    endif
    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine primaryRetr

  subroutine primaryRetrDestroy()
    !Local
    character (len=*), parameter :: procName = &
                                   ' [PrimaryRetrModule::primaryRetrDestroy]: '
    if (dbg) print *, procName, 'starting ...'
    if (allocated(pRef)) deallocate(pRef)
    if (allocated(press)) deallocate(press)
    if (allocated(MolTran)) deallocate(MolTran)
    if (allocated(kchan)) deallocate(kchan)
    if (allocated(rad)) deallocate(rad)
    if (allocated(radObs)) deallocate(radObs)
    if (allocated(Sy)) deallocate(Sy)
    if (allocated(delY)) deallocate(delY)
    if (allocated(Ylast)) deallocate(Ylast)
    if (allocated(xGesLast)) deallocate(xGesLast)
    if (allocated(rerr)) deallocate(rerr)
    if (allocated(atmNoise)) deallocate(atmNoise)
    if (allocated(dummyArr)) deallocate(dummyArr)
    if (allocated(EmRf)) deallocate(EmRf)
    if (allocated(xkEmRf)) deallocate(xkEmRf)
    if (allocated(xBakSav)) deallocate(xBakSav)
    if (allocated(xGesNew)) deallocate(xGesNew)
    if (allocated(xktG)) deallocate(xktG)
    if (allocated(xkt)) deallocate(xkt)
    if (allocated(background)) deallocate(background)
    if (allocated(clw_cov)) deallocate(clw_cov)
    if (allocated(ciw_cov)) deallocate(ciw_cov)
    if (allocated(ut_cov_net)) deallocate(ut_cov_net)
    if (allocated(covRetr)) deallocate(covRetr)
    if (allocated(dcovRetr)) deallocate(dcovRetr)
    if (allocated(umtx)) deallocate(umtx)
    if (allocated(vmtx)) deallocate(vmtx)
    if (allocated(xGesG)) deallocate(xGesG)
    if (allocated(xGes)) deallocate(xGes)
    if (allocated(xkEmMw)) deallocate(xkEmMw)
    if (allocated(emMW)) deallocate(emMW)
    if (allocated(xOffG)) deallocate(xOffG)
    if (allocated(xGainG)) deallocate(xGainG)
    if (allocated(dradChan)) deallocate(dradChan)
    call bkgObj%bkgDestroy()
    if (dbg) print *, procName, 'ending ...'
  end subroutine primaryRetrDestroy

  !----------------------------------------------------------------------
  ! function performs compression of computed radiances, Jacobians
  !
  !
  subroutine compressRadJacobian(kchan, nChMW, nChIR, y, xkt, xkEmMw)
    logical, dimension(:),               intent(in)    :: kchan
    integer,                             intent(in)    :: nChMW
    integer,                             intent(in)    :: nChIR
    real,    dimension(:),               intent(inout) :: y
    real,    dimension(:,:),             intent(inout) :: xkt
    real,    dimension(:),               intent(inout) :: xkEmMw

    ! Local
    integer                         :: nChanTot
    integer                         :: k, jj
    integer                         :: idx
    character(len=*), parameter     :: procName=&
                               ' [PrimaryRetrModule::compressRadJacobian]: '

    nChanTot = nChMW + nChIR
    idx = 0

    ! process MW (they are uncompressed and located at kchan=.true.)
    do k=1,nChMW
      if (kchan(k)) then
        idx = idx + 1
        y(idx) = y(k)
        xkt(:,idx) = xkt(:,k)
        xkEmMw(idx) = xkEmMw(k)
      end if
    end do

    ! process IR (they are already compressed)
    jj=nChMW
    do k=nChMW+1,nChanTot
      if (kchan(k)) then
        jj=jj+1
        idx = idx + 1
        y(idx) = y(jj)
        xkt(:,idx) = xkt(:,jj)
      end if
    end do
  end subroutine compressRadJacobian

END MODULE PrimaryRetrModule
