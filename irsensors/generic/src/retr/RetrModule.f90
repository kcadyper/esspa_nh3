MODULE RetrModule
  !
  ! Module that contains the major components of retrieval
  !
  ! Public Subroutines:
  !
  !     getControl
  !     init
  !     retr
  !     destroy
  !     MWRTverify
  !     IRRTverify
  !
  ! Private Subroutines:
  !
  !     setPropertiesAlgConfig
  !     loadAlgconfig
  !
  ! Derived data types:
  !
  !     GenControl_t
  !     RetrResult_t
  !     RetrOutput_t
  !     LinearResult_t
  !     StateSpaceDefin_t
  !
  ! USE:
  !
  !     constants
  !     ToolboxModule
  !     StateIndexModule
  !     SurfaceClassModule
  !     AtmosClassModule
  !     ChannelSelectionModule
  !     MWRTmodule
  !     IRRTmodule
  !     IMGRTmodule
  !     StatRetrModule
  !     MWretrModule
  !     IRretrModule
  !     PrimaryRetrModule
  !     SecondaryRetrModule
  !     MWmeasErrModule
  !     IRmeasErrModule
  !     LimitsModule
  !     CloudMitigateModule
  !     MWobsStructure
  !     IRobsStructure
  !     ImgObsStructure
  !     oss_mw_module
  !
  ! yhe@aer.com, 03/16/2016
  !
  USE constants, Only: &
       MISSING_INT, &
       MISSING_REAL, &
       MISSING_CHAR

  USE StateIndexModule, Only: &
       StateIndex_t, &
       StateFlags_t, &
       maxMol, &
       initLengths, &
       genIndices, &
       charToMolID

  USE SurfaceClassModule, Only: &
       MWsurfaceClassInit, &
       MWsurfaceClassify, &
       IRsurfaceClassInit, &
       IRsurfaceClassify

  USE AtmosClassModule, Only: &
       atmosClassInit, &
       atmosClassify, &
       MWmode, &
       IRmode, &
       MWIRmode

  USE ChannelSelectionModule, Only: &
       ChannelSet_t, &
       ChannelSelections_t, &
       channelSelectionInit, &
       channelSetUnion, &
       channelSelectionDestroy

  USE MWRTmodule, Only: &
       MWRTinit, &
       MWRTexec, & !(calls ossdrv_mw)
       MWRTdestroy, &
       MWRTcontrol_t

  USE IRRTmodule, Only: &
       IRRTinit, &
       IRRTexec, & !(calls ossdrv_ir)
       IRRTdestroy, &
       IRRTcontrol_t

  USE ImgRTmodule, Only: &
       ImgRTinit, &
       ImgRTexec, &
       ImgRTdestroy, &
       imgRTcontrol_t

  USE StatRetrModule, Only: &
       statRetrInit, &
       statRetr, &
       statRetrDestroy, &
       statRetrControl_t

  USE MWretrModule, Only: &
       MWretrInit, &
       MWretr, &
       MWretrDestroy, &
       MWretrControl_t

  USE PrimaryRetrModule, Only: &
       primaryRetrInit, &
       primaryRetr, &
       primaryRetrDestroy, &
       PrimaryRetrControl_t

  USE SecondaryRetrModule, Only: &
       secondaryRetrInit, &
       secondaryRetr, &
       secondaryRetrDestroy, &
       SecondaryRetrControl_t

  USE MeasErrStructure, Only: &
       MeasErr_t

  USE MWmeasErrModule, Only: &
       MWnoiseInit, &
       setMWmeasErr, &
       MWnoiseDestroy

  USE IRmeasErrModule, Only: &
       IRnoiseInit, &
       setIRmeasErr, &
       IRnoiseDestroy

  USE LimitsModule, Only: &
       limitsInit, & !(calls loadGuessLimits)
       imposeLimits, & !(calls chkges)
       limitsDestroy

  USE CloudMitigateModule, Only: &
       cloudClassifyInitial

  USE MWobsStructure, Only: &
       MWgran_t, &
       MWdefin_t, &
       MWob_t

  USE IRobsStructure, Only: &
       IRgran_t, &
       IRdefin_t, &
       IRFOR_t

  USE ImgObsStructure, Only: &
       ImgDef_t, &
       ImgOb_t, &
       ImgGran_t

  USE AncillaryStructure, Only: &
       AncillaryGran_t, &
       AncillaryDefin_t, &
       AncillaryFOR_t

  USE ControlStructure, Only: &
       GenControl_t, &
       StateSpaceDefin_t, &
       RetrResult_t, &
       LinearResult_t, &
       stateSpacesMatch

  USE VertCoord, Only: &
       mxCoordTyp, Pcoord

  USE OutputStructure, Only: &
       OutputFOR_t, &
       OutputDefin_t, &
       OutputGran_t, &
       LinearGran_t, &
       LinearFOR_t

  USE IOmodule, Only: &
       destroyMWsensorData, &
       destroyIRsensorData, &
       destroyImgSensorData, &
       destroyAncillaryData

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: &
       getControl, &
       init, &
       retr, &
       destroy

  INTEGER, PARAMETER :: mxLength=256
  INTEGER, PARAMETER :: mxchan = 890
  INTEGER, PARAMETER :: NCDF_DUMMY=999
  INTEGER, PARAMETER :: maxBkg = 3
  INTEGER, PARAMETER :: maxChSet = 1
  REAL, PARAMETER :: deltaP=0.05    !units: hpa
  REAL, PARAMETER :: deltaFrq=0.01  !units: GHz
  REAL, PARAMETER :: deltaWvn=0.001 !units: cm-1

  logical :: DynamicNoise=.FALSE.
  logical :: loadAlgCfgDone=.FALSE.
  logical :: enableScatter=.FALSE.
  logical :: extFlgTemp=.FALSE.
  logical :: extFlgTskin=.FALSE.
  logical :: extFlgPsfc=.FALSE.
  logical :: extFlgWind=.FALSE.
  logical :: extFlgCldLiq=.FALSE.
  logical :: extFlgCldIce=.FALSE.
  logical :: statRetrOn=.FALSE.
  logical :: primaryRetrOn=.FALSE.
  logical :: secondaryRetrOn=.FALSE.
  logical :: MWretrOn=.FALSE.
  logical :: primaryMWon=.FALSE.
  logical :: PhysOn=.FALSE.
  logical :: imgDataOn=.FALSE.
  logical :: ossScatRegrIRon=.FALSE.
  logical :: extFlgEmIR=.FALSE.
  logical :: genLinInvert=.FALSE.
  logical :: logSize
  logical :: independentFOV=.TRUE.
  logical, DIMENSION(maxMol) :: extFlgMol=.FALSE.
  integer :: apodizationType=0
  integer :: nTskin
  integer :: nTemp
  integer :: nMol
  integer :: nEmMW
  integer :: nEmIR
  integer :: nCldLiq
  integer :: nCldIce
  integer :: nxIterMw
  integer :: nxIter
  integer :: nParMax_MW
  integer :: nParMax_IR
  integer :: nParGmax
  integer :: nEmMWmax
  integer :: nEmRfIRmax
  integer, dimension(maxMol) :: ngas
  integer, dimension(:), allocatable :: MolID
  integer, dimension(:), allocatable :: MolTran
  real :: drad
  real :: alpha1slope
  real :: alpha1int
  real :: chisqconvMw
  real :: chisqconv
  real :: chisqRatioconv
  real :: cldamt1stguess
  real :: iceamt1stguess
  real :: TightCLDcov
  real :: plandMaxOc
  TYPE(RetrResult_t) :: retrResult
  TYPE(RetrResult_t) :: retrPrior
  integer, dimension(:), allocatable :: cloudClass
  TYPE(LinearResult_t) :: linearResult
  integer :: nSav
  real, dimension(:), allocatable :: retrSavIRsfc
  real, dimension(:), allocatable :: backResult
  character (len=mxCoordTyp) :: vCoordTyp

  logical :: dbg=.false.

CONTAINS

  subroutine getControl(&
       genControl, &
       statRetrControl, &
       MWretrControl, &
       primaryRetrControl, &
       secondaryRetrControl, &
       MWRTcontrol, &
       IRRTcontrol, &
       imgRTcontrol)

    TYPE(GenControl_t), INTENT(OUT)           ::  genControl
    TYPE(statRetrControl_t), INTENT(OUT)      ::  statRetrControl
    TYPE(MWretrControl_t), INTENT(OUT)        ::  MWretrControl
    TYPE(primaryRetrControl_t), INTENT(OUT)   ::  primaryRetrControl
    TYPE(secondaryRetrControl_t), INTENT(OUT) ::  secondaryRetrControl
    TYPE(MWRTcontrol_t), INTENT(OUT)          ::  MWRTcontrol
    TYPE(IRRTcontrol_t), INTENT(OUT)          ::  IRRTcontrol
    TYPE(imgRTcontrol_t), INTENT(OUT)         ::  imgRTcontrol

    !---Local variables
    integer :: nprofs
    real, dimension(:), allocatable :: psfc
    real, dimension(:), allocatable :: plandAvg
    REAL                       :: pobs=0.0
    INTEGER                    :: scatterSizeOpt=0
    CHARACTER(LEN=12), DIMENSION(maxMol) :: molOnRT
    CHARACTER(LEN=mxLength) :: AuxFile,RadMWfile,cloudTabFile,ldEmissExtern,bkgrdExtern
    CHARACTER(LEN=mxLength) :: nedtFile,noiseReductFactFile='NULL', &
         osscoefs_mw,ossoptdpth_mw,classAtmosConfig,classSfcMWconfig,classSfcIRconfig, &
         algConfig, RetrMWfile,GuessFile
    CHARACTER(LEN=mxLength), DIMENSION(maxBkg) :: bkgDataFiles
    CHARACTER(LEN=mxLength)    :: osscoefs_ir,ossoptdpth_ir
    CHARACTER(LEN=mxLength)    :: osspred_ir
    CHARACTER(len=mxLength)    :: chanFileName,nednFile
    CHARACTER(LEN=mxLength)    :: RadIRfile,solarfluxFile,regrFile,RetrRegrFile,RetrIrFile
    CHARACTER(LEN=mxLength)    :: defProfFile
    CHARACTER(LEN=mxLength)    :: linInvertFile,IRemissExtFile
    integer :: imol, icov, ibkg
    logical                    :: Planck2ndOrderTaylor
    character (len=*), parameter  :: procName=' [RetrModule::getControl]: '

    NAMELIST /RTinput/ molOnRT,cloudTabFile,solarfluxFile,defProfFile,&
         Planck2ndOrderTaylor
    NAMELIST /externControl/  &
         ldEmissExtern,IRemissExtFile,bkgrdExtern,extFlgTemp,extFlgTskin, &
         extFlgPsfc,extFlgWind,extFlgCldLiq,extFlgCldIce,extFlgMol,extFlgEmIR
    NAMELIST /sensorSpecific/ nedtFile,noiseReductFactFile,osscoefs_mw,  &
         ossoptdpth_mw,osscoefs_ir,ossoptdpth_ir,osspred_ir, &
         classAtmosConfig,classSfcMWconfig,classSfcIRconfig,chanFileName,nednFile
    NAMELIST /generalControl/ nprofs,enableScatter,statRetrOn,MWretrOn, &
         primaryRetrOn,secondaryRetrOn, &
         imgDataOn,algConfig,bkgDataFiles,GuessFile,AuxFile,RadMWfile,RadIRfile, &
         RetrRegrFile, &
         RetrMWfile,RetrIRfile,linInvertFile,independentFOV,apodizationType
    NAMELIST /statControl/ regrFile
    NAMELIST /primaryControl/ primaryMWon,genLinInvert

    if (dbg) print *, procName, 'starting ...'

    molOnRT=MISSING_CHAR  ! all elements not given are marked as unused
    Planck2ndOrderTaylor = .false.
    ! The namelist may come from standard input (default) or a specified file (optional).
    READ(*,generalControl)
    READ(*,externControl)
    READ(*,sensorSpecific)
    READ(*,RTinput)
    READ(*,statControl)
    READ(*,primaryControl)

    !---Load retr tuning parameters
    genControl%backgroundConfig = algConfig
    genControl%limitsConfig = algConfig
    call loadAlgconfig(GenControl%limitsConfig)
    !----Derive nMol and molID from inputs
    nMol=COUNT(molOnRT /= MISSING_CHAR)
    if (nMol <= 0) then
       print *,'err:[retr] No molecules selected for RT'
       call errorHalt(1)
    endif

    !need to assign molID and molTran to StateSpaceDefin_t defined
    !variable(s)
    allocate (MolID(nMol),MolTran(nMol))
    call charToMolID(molOnRT(1:nMol),MolID)

    !-- genControl data structure assignment
    genControl%nFORmax = nprofs
    if (.not. allocated(genControl%kchan)) allocate(genControl%kchan(mxchan))
    genControl%enableScatter = enableScatter
    genControl%statRetrOn = statRetrOn
    genControl%MWRetrOn = MWretrOn
    genControl%primaryRetrOn = primaryRetrOn
    genControl%secondaryRetrOn = secondaryRetrOn
    if (MWretrOn .or. primaryRetrOn .or. secondaryRetrOn) PhysOn = .TRUE.
    if (MWretrOn .or. (primaryRetrOn .and. primaryMWon)) then
       genControl%MWdataOn=.TRUE.
    else
       genControl%MWdataOn=.FALSE.
    endif
    if (primaryRetrOn .or. secondaryRetrOn) then
       genControl%IRdataOn = .TRUE.
    else
       genControl%IRdataOn = .FALSE.
    endif
    genControl%imgDataOn = imgDataOn
    genControl%LWP1stGuess = cldamt1stguess
    genControl%IWP1stGuess = iceamt1stguess
    genControl%logSize = logSize
    genControl%atmosClassConfig = ''
    genControl%MWsurfaceClassConfig = ''
    genControl%IRsurfaceClassConfig = ''
    genControl%MWnoiseFile = trim(nedtFile)
    genControl%MWnoiseAvgFactorFile = trim(noiseReductFactFile)
    genControl%IRnoiseFile = trim(nednFile)
    genControl%MWobsFile = trim(RadMWfile)
    genControl%IRobsFile = trim(RadIRfile)
    genControl%imgObsFile = ''
    genControl%ancillaryAtmosFile = trim(bkgrdExtern)
    genControl%ancillarySfcFile = trim(AuxFile)
    genControl%ancillaryMWemisFile = trim(ldEmissExtern)
    genControl%ancillaryIRemisFile = trim(IRemissExtFile)
    genControl%statRetrOutFile = trim(RetrRegrFile)
    genControl%MWretrOutFile = trim(RetrMWfile)
    genControl%primaryRetrOutFile = trim(RetrIRfile)
    genControl%secondaryRetrOutFile = ''
    genControl%backgroundOutFile = trim(GuessFile)
    genControl%linearInvertOutFile = trim(linInvertFile)
    genControl%independentFOV = independentFOV
    genControl%apodizationType = apodizationType
    if (.not. allocated(genControl%backgroundFiles)) &
         allocate(genControl%backgroundFiles(maxBkg))
    do ibkg = 1, maxBkg
       genControl%backgroundFiles(ibkg) = trim(bkgDataFiles(ibkg))
    end do
    if (.not. allocated(genControl%channelSelectFiles)) &
         allocate(genControl%channelSelectFiles(1))
    genControl%channelSelectFiles(1) = trim(chanFileName)
    genControl%externDataFlags%temp = extFlgTemp
    genControl%externDataFlags%Tskin = extFlgTskin
    genControl%externDataFlags%mol = .FALSE.
    do imol = 1, nMol
       genControl%externDataFlags%mol(imol) = extFlgMol(imol)
    end do
    genControl%externDataFlags%pSfc = extFlgPsfc
    genControl%externDataFlags%wind = extFlgWind
    genControl%externDataFlags%cloudLiq = extFlgCldLiq
    genControl%externDataFlags%cloudIce = extFlgCldIce
    genControl%externDataFlags%emMW = .TRUE.
    genControl%externDataFlags%emRfIR = extFlgEmIR
    genControl%atmosClassConfig = classAtmosConfig
    genControl%MWsurfaceClassConfig = classSfcMWconfig
    genControl%IRsurfaceClassConfig = classSfcIRconfig

    !-- MWRTcontrol
    MWRTcontrol%coefFile = trim(osscoefs_mw)
    MWRTcontrol%LUTfile = trim(ossoptdpth_mw)
    MWRTcontrol%nMol = nMol
    MWRTcontrol%MolID = MISSING_INT
    MWRTcontrol%MolID(1:nMol) = MolID
    MWRTcontrol%radOrTb = MISSING_INT
    MWRTcontrol%Planck2ndOrderTaylor=Planck2ndOrderTaylor

    !-- IRRTcontrol
    IRRTcontrol%coefFile = trim(osscoefs_ir)
    IRRTcontrol%LUTfile = trim(ossoptdpth_ir)
    IRRTcontrol%cloudOptFile = trim(cloudTabFile)
    IRRTcontrol%solarFile = trim(solarfluxFile)
    IRRTcontrol%defaultProfFile = trim(defProfFile)
    IRRTcontrol%scatterCoefFile = trim(osspred_ir)
    IRRTcontrol%nMol = nMol
    IRRTcontrol%MolID = MISSING_INT
    IRRTcontrol%MolID(1:nMol) = MolID
    IRRTcontrol%scatterApprox = .FALSE.
    IRRTcontrol%nHyd = 0  !temporary - TBD by scene or cloud table file?
    if (IRRTcontrol%nHyd > 0) then
       if (.not. allocated(IRRTcontrol%dataTypeCloud)) &
            allocate(IRRTcontrol%dataTypeCloud(IRRTcontrol%nHyd))
       if (.not. allocated(IRRTcontrol%cloudModel)) &
            allocate(IRRTcontrol%cloudModel(IRRTcontrol%nHyd))
       IRRTcontrol%dataTypeCloud = MISSING_CHAR !temporary - to be replaced by the input in the namelist?
       IRRTcontrol%cloudModel = MISSING_INT  !temporary - to be replaced by the input in the namelist?
    end if

    !-- statRetrControl
    statRetrControl%statRetrCoefFile = trim(regrFile)

    !-- MWretrControl
    MWretrControl%nRetr%temp = nTemp
    MWretrControl%nRetr%Tskin = nTskin
    MWretrControl%nRetr%mol = nGas
    MWretrControl%nRetr%cldLiq = nCldLiq
    MWretrControl%nRetr%cldIce = nCldIce
    MWretrControl%nRetr%EmMW = nEmMW
    MWretrControl%maxIter = nxiterMw
    MWretrControl%drad = drad
    MWretrControl%alpha1slope = alpha1slope
    MWretrControl%alpha1int = alpha1int
    MWretrControl%chiSqConv = chiSqConvMW
    MWretrControl%chisqRatioConv = chisqRatioConv
    !-- primaryRetrControl
    primaryRetrControl%nRetr%temp = nTemp
    primaryRetrControl%nRetr%Tskin = nTskin
    primaryRetrControl%nRetr%mol = nGas
    primaryRetrControl%nRetr%cldLiq = nCldLiq
    primaryRetrControl%nRetr%cldIce = nCldIce
    primaryRetrControl%nRetr%EmMW = nEmMW
    primaryRetrControl%nRetr%EmRfIR = nEmIR
    primaryRetrControl%maxIter = nxiter
    primaryRetrControl%drad = drad
    primaryRetrControl%chiSqConv = chiSqConv
    primaryRetrControl%chisqRatioConv = chisqRatioConv
    primaryRetrControl%genLinInvert = genLinInvert  !on/off switch for linear inversion process
    primaryRetrControl%MWon = primaryMWon
    !-- secondaryRetrControl
    !-- imgRTcontrol

    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine getControl

  subroutine init( &
       genControl, &
       statRetrControl, &
       MWretrControl, &
       primaryRetrControl, &
       secondaryRetrControl, &
       MWRTcontrol, &
       IRRTcontrol, &
       imgRTcontrol, &
       MWdefin, &
       IRdefin, &
       ancillaryDefin, &
       stateSpaceDefinStat, &
       stateSpaceDefinMW, &
       stateSpaceDefinPrimary, &
       stateSpaceDefinSecondary, &
       ChannelSelections &
       )
    TYPE(GenControl_t), intent(in) :: genControl
    TYPE(StatRetrControl_t), intent(in) :: statRetrControl
    TYPE(MWretrControl_t), intent(in) ::  MWretrControl
    TYPE(PrimaryRetrControl_t), intent(in) ::  primaryRetrControl
    TYPE(secondaryRetrControl_t), intent(in) :: secondaryRetrControl
    TYPE(MWdefin_t), intent(inout) :: MWdefin
    TYPE(IRdefin_t), intent(inout) :: IRdefin
    TYPE(AncillaryDefin_t), intent(in) :: ancillaryDefin
    TYPE(StateSpaceDefin_t), intent(out) :: stateSpaceDefinStat
    TYPE(StateSpaceDefin_t), intent(out) :: stateSpaceDefinMW
    TYPE(StateSpaceDefin_t), intent(out) :: stateSpaceDefinPrimary
    TYPE(StateSpaceDefin_t), intent(out) :: stateSpaceDefinSecondary
    TYPE(ChannelSelections_t), intent(out) :: ChannelSelections
    TYPE(MWRTcontrol_t), intent(inout) :: MWRTcontrol
    TYPE(IRRTcontrol_t), intent(inout) :: IRRTcontrol
    TYPE(ImgRTcontrol_t), intent(in) :: imgRTcontrol
    real :: maxLandInOceanClassMW !land and ocean classification threshold?
    real :: maxLandInOceanClassIR !land and ocean classification threshold?
    ! local to init
    TYPE(StateSpaceDefin_t) :: stateSpaceDefinMWRT
    TYPE(StateSpaceDefin_t) :: stateSpaceDefinIRRT
    real, dimension(:), allocatable :: frqRT
    real, dimension(:), allocatable :: pRefRT
    integer, dimension(:), allocatable :: polRT
    real, dimension(:), allocatable :: wvnRT
    integer :: nchanRTmw
    integer :: nLevRTmw
    integer :: nchanRTir
    integer :: nLevRTir
    logical :: stateSpaceMatchStatus
    integer :: atmosClassMode
    character (len=*), parameter :: procName = ' [RetrModule::init]: '

    if (dbg) print *, procName, ' starting ...'
    if (genControl%MWdataOn) &
       call MWsurfaceClassInit(genControl)
    if (genControl%IRdataOn) &
       call IRsurfaceClassInit(genControl)
    if (genControl%MWdataOn .and. genControl%IRdataOn) then
       atmosClassMode=MWIRmode
    else
       if (genControl%MWdataOn) atmosClassMode=MWmode
       if (genControl%IRdataOn) atmosClassMode=IRmode
    endif

    call channelSelectionInit(genControl,MWdefin,IRdefin,channelSelections)

    if (genControl%MWdataOn) then
       if (.not. allocated(frqRT)) allocate(frqRT(MWdefin%nChan))
       if (.not. allocated(polRT)) allocate(polRT(MWdefin%nChan))
! These definitions can be removed when OSS v1.2 is integrated and stateSpaceDefin
! is no longer needed by MWRTinit for calling oss_init_mw. Already these definitions
! are not used by oss_init_mw, except setNparG temporarily is needed below.
! Also, the use statements
! for initLengths and genIndices should be removed
       stateSpaceDefinMWRT%NG=initLengths()
       stateSpaceDefinMWRT%IG=genIndices(stateSpaceDefinMWRT%NG)
       call MWRTinit(MWRTcontrol,stateSpaceDefinMWRT,genControl%enableScatter, &
            MWdefin,nChanRTmw,frqRT,polRT,nLevRTmw,pRefRT)
    end if

    nParGMax = 0
    nEmMWmax = 0
    if (genControl%MWretrOn) then
       stateSpaceDefinMW%MolID = MISSING_INT
       stateSpaceDefinMW%MolID(1:nMol) = MolID
       call MWretrInit( &
            genControl, &
            MWretrControl, &
            MWRTcontrol, &
            MWdefin, &
            stateSpaceDefinMW &
            )
       nParGMax = max(nParGMax, stateSpaceDefinMW%nParG)
       nEmMWmax = max(nEmMWmax, stateSpaceDefinMW%NG%emMW)
       stateSpaceDefinMWRT=stateSpaceDefinMW
    end if

    nEmRfIRmax = 0
    if (genControl%primaryRetrOn) then
       stateSpaceDefinPrimary%MolID = MISSING_INT
       stateSpaceDefinPrimary%MolID(1:nMol) = MolID
       call primaryRetrInit( &
            genControl, &
            primaryRetrControl, &
            MWRTcontrol, &
            IRRTcontrol, &
            MWdefin, &
            IRdefin, &
            ancillaryDefin, &
            stateSpaceDefinPrimary &
            )
       nParGMax = max(nParGMax, stateSpaceDefinPrimary%nParG)
       nEmMWmax = max(nEmMWmax, stateSpaceDefinPrimary%NG%emMW)
       nEmRfIRmax = max(nEmRfIRmax, stateSpaceDefinPrimary%NG%emRfIR)
       stateSpaceDefinIRRT=stateSpaceDefinPrimary

       if (genControl%MWretrOn) then
          if (stateSpaceDefinMW%nLev /= stateSpaceDefinPrimary%nLev) then
             print *, procName, ' err - inconsistent nLev b/w MW and primary', &
                  stateSpaceDefinMW%nLev, stateSpaceDefinPrimary%nLev
             call errorHalt(-1)
          end if

          if ( .not. stateSpacesMatch(stateSpaceDefinMW,stateSpaceDefinPrimary,&
                            atmos=.True.,sfcMW=primaryRetrControl%MWon) ) then
             print *, procName, ' err - Mismatch MW and primary state spaces'
             call errorHalt(-1)
          end if
       else
          if (genControl%MWdataOn) stateSpaceDefinMWRT=stateSpaceDefinPrimary
       end if
    end if

    if(genControl%secondaryRetrOn) then
       call secondaryRetrInit( &
            genControl, &
            secondaryRetrControl, &
            IRRTcontrol, &
            stateSpaceDefinSecondary &
            )
       if (genControl%primaryRetrOn) then
          stateSpaceMatchStatus = &
               stateSpacesMatch(stateSpaceDefinSecondary,stateSpaceDefinPrimary, &
                   atmos=.True.,sfcMW=.True.,sfcIR=.True.)
          if (.not. stateSpaceMatchStatus) then
             print *, procName, ' err - Mismatch primary and secondary state spaces'
             call errorHalt(-1)
          end if
       else
          stateSpaceDefinIRRT=stateSpaceDefinSecondary
       end if
    end if

    if (genControl%MWdataOn) then
       vCoordTyp = trim(stateSpaceDefinMWRT%vCoordTyp)
       call MWRTverify(MWRTcontrol,MWdefin,stateSpaceDefinMWRT,nChanRTmw, &
            frqRT,polRT,nLevRTmw,pRefRT)
    end if

    !note: nChanRT, wvnRT, nLevRT, pRefRT here are for IR set - may be
    !named differently from those of MW, if both sets are used further
    !down the line.
    !
    ! added condition of primaryRetrOn, to be consistent with the
    ! availability of the primary fields. (05/09/2017, YHE)
    if (genControl%IRdataOn) then
       if (.not. allocated(wvnRT)) allocate(wvnRT(IRdefin%nChan))
       vCoordTyp = trim(stateSpaceDefinIRRT%vCoordTyp)
       call IRRTinit(&
            IRRTcontrol, &
            stateSpaceDefinIRRT , &
            genControl%enableScatter, &
            ChannelSelections, &
            nChanRTir, &
            wvnRT, &
            nLevRTir, &
            pRefRT &
            )

       call IRRTverify( &
            IRdefin, &
            stateSpaceDefinIRRT, &
            nChanRTir, &
            wvnRT, &
            nLevRTir, &
            pRefRT &
            )

    !this is temporary – will be superseded by new channel selection;
    !GetIRsensorData will ensure expected channel order
    !note: IRgran%FOR%ob%rad is unavaible to this subroutine - hold off
!!$    IF (IRdefin%nChan /= nChanRT)then
!!$       call adjustChanSet(IRgran%FOR%ob%rad,wvnRad,wvnOSS(1:nChanIR), &
!!$            nchanIRrad,nChanIR,kchan,radAll,chanIDtmp,chanIDir)
!!$    else
!!$       radAll = radTmp
!!$       chanIDir = chanIDtmp
!!$    endif
!!$    deallocate(radTmp,chanIDtmp)
    end if

    call ImgRTinit(imgRTcontrol)
    call ImgRTverify()

    call atmosClassInit(genControl,MWdefin,IRdefin,atmosClassMode)
    print *,'Irdefin%frq'
    print *,Irdefin%wvn(501:511)
    !notes: use of either stateSpaceDefinMW or stateSpaceDefinPrimary
    !depends on the MWretrOn and/or PrimaryRetrOn - needs to be
    !addressed, YHE, 05/05/2017
    if (genControl%statRetrOn) then
       CALL statRetrInit( &
            genControl, &
            statRetrControl, &
            MWdefin, &
            IRdefin, &
            stateSpaceDefinStat &
            )
       if (genControl%MWretrOn) then
          stateSpaceMatchStatus = &
               stateSpacesMatch(stateSpaceDefinMW,stateSpaceDefinStat,atmos=.True.,sfcMW=.True.)
       end if
       if (genControl%primaryRetrOn) then
          stateSpaceMatchStatus = &
               stateSpacesMatch(stateSpaceDefinPrimary,stateSpaceDefinStat,atmos=.True.,sfcMW=.True.,sfcIR=.True.)
       end if
    end if

    if (genControl%MWdataOn) call MWnoiseInit(MWdefin)

    if (genControl%IRdataOn) then
       call IRnoiseInit( &
            genControl, &
            IRdefin, &
            nChanRTir &
            )
    end if

    call limitsInit( &
         genControl &
!         genControl, &
!         stateSpaceDefinPrimary &
         )

    ! Comments from the design document as of 04/27/2016: "this is
    ! temporary, until scattering model with Jacobians is integrated I
    ! think this block needs to move down, to come after NR_orig has
    ! been defined."
    !
    ! notes: NR_orig is assigned with NR as indicated in
    ! retr_sig.f90, which is defined within MWretrInit(...) by the
    ! call to setXPtrR(...) per design document. Therefore, this block
    ! can be revisited when MWretrInit(...) is addressed.

!!$    if (genControl%enableScatter .AND. ANY( &
!!$         (/NR_orig%temp,NR_orig%mol(1:nmol), &
!!$         nEmMw_orig,nEmIR_orig/) /= 0)) then
!!$       print*,'err:[retr] Jacobians not available for all retrieved '// &
!!$            'variables in scattering mode; '
!!$       print*,'See config file entries for temp, Tskin, gas, emis'
!!$       call errorHalt(1)
!!$    endif
!!$
!!$    nchan=nchanmw+nchanir
!!$    if (genControl%enableScatter) then
!!$       call loadFiniteDiff(F_algconfig)
!!$       call initFiniteDiff()
!!$    endif

    if (allocated(frqRT)) deallocate(frqRT)
    if (allocated(polRT)) deallocate(polRT)
    if (allocated(pRefRT)) deallocate(pRefRT)
    if (allocated(wvnRT)) deallocate(wvnRT)

    retrResult%nParGmax = nParGmax
    retrResult%nEmMWmax = nEmMWmax
    retrResult%nEmRfIRmax = nEmRfIRmax
    if (.not. allocated(retrResult%press)) then
       if (genControl%MWretrOn .and. (.not. genControl%primaryRetrOn) ) then
          allocate(retrResult%press(stateSpaceDefinMW%nLev))
       else
          allocate(retrResult%press(stateSpaceDefinPrimary%nLev))
       end if
    end if
    if (.not. allocated(retrResult%stateV)) allocate(retrResult%stateV(nParGmax))
    if (.not. allocated(retrResult%emMW))   allocate(retrResult%emMW(nEmMWmax))
    if (.not. allocated(retrResult%emRfIR)) allocate(retrResult%emRfIR(nEmRfIRmax))

    retrPrior%nParGmax = nParGmax
    retrPrior%nEmMWmax = nEmMWmax
    retrPrior%nEmRfIRmax = nEmRfIRmax
    if (.not. allocated(retrPrior%press)) then
       if (genControl%MWretrOn .and. (.not. genControl%primaryRetrOn) ) then
          allocate(retrPrior%press(stateSpaceDefinMW%nLev))
       else
          allocate(retrPrior%press(stateSpaceDefinPrimary%nLev))
       end if
    end if
    if (.not. allocated(retrPrior%stateV)) allocate(retrPrior%stateV(nParGmax))
    if (.not. allocated(retrPrior%emMW))   allocate(retrPrior%emMW(nEmMWmax))
    if (.not. allocated(retrPrior%emRfIR)) allocate(retrPrior%emRfIR(nEmRfIRmax))
    if (.not. allocated(retrSavIRsfc))     allocate(retrSavIRsfc(nEmRfIRmax))

   if (genControl%primaryRetrOn) then
       if (.not. allocated(backResult)) allocate(backResult(nParGmax))
       if (.not. allocated(linearResult%xOffG)) allocate(linearResult%xOffG(nParGmax))
       if (.not. allocated(linearResult%xGainG)) then
          if (genControl%MWdataOn) then
             allocate(linearResult%xGainG(nParGmax,MWdefin%nChan+IRdefin%nChan))
          else
             allocate(linearResult%xGainG(nParGmax,IRdefin%nChan))
          endif
       endif
    end if

    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine init

  subroutine retr( &
       genControl, &
       statRetrControl, &
       MWretrControl, &
       primaryRetrControl, &
       secondaryRetrControl, &
       MWRTcontrol, &
       IRRTcontrol, &
       imgRTcontrol, &
       MWdefin, &
       MWob, &
       IRdefin, &
       IRFOR, &
       ImgDef, &
       ImgOb, &
       ancillaryDefin, &
       ancillaryFOR, &
       channelSelections, &
       retrOutputStatFOR, &
       retrOutputMWFOR, &  !Technically, should be retrOutputMWFOV, for FOR is not associated with MW
       retrOutputPrimaryFOR, &
       linearOutputFOR, &
       backgroundOutputFOR, &
       retrOutputSecondaryFOR &
       )

    TYPE(GenControl_t), intent(in) :: genControl
    TYPE(StatRetrControl_t), intent(in) :: statRetrControl
    TYPE(MWretrControl_t), intent(in) ::  MWretrControl
    TYPE(PrimaryRetrControl_t), intent(in) ::  primaryRetrControl
    TYPE(secondaryRetrControl_t), intent(in) :: secondaryRetrControl
    TYPE(MWRTcontrol_t), intent(in) :: MWRTcontrol
    TYPE(IRRTcontrol_t), intent(in) :: IRRTcontrol
    TYPE(ImgRTcontrol_t), intent(in) :: imgRTcontrol
    TYPE(MWdefin_t), intent(in) :: MWdefin
    TYPE(MWOb_t), intent(inout) :: MWob
    TYPE(IRdefin_t), intent(in) :: IRdefin
    TYPE(IRFOR_t), intent(inout) :: IRFOR
    TYPE(ImgDef_t), intent(in) :: ImgDef
    TYPE(ImgOb_t), dimension(:), intent(in) :: ImgOb
    TYPE(AncillaryDefin_t), intent(in) :: ancillaryDefin
    TYPE(AncillaryFOR_t), intent(in) :: AncillaryFOR
    TYPE(ChannelSelections_t), intent(in) :: channelSelections
    TYPE(OutputFOR_t), dimension(:), intent(inout) :: retrOutputStatFOR !retrOutputStat%FOR2D(:,iFOR)
    TYPE(OutputFOR_t), intent(inout) :: retrOutputMWFOR !retrOutputMW%FOR(ifor)
    TYPE(LinearFOR_t), intent(inout) :: linearOutputFOR !linearOutput%FOR2D(:,iFOR)
    TYPE(OutputFOR_t), intent(inout) :: backgroundOutputFOR !backgroundOutput%FOR(iFOR)
    TYPE(OutputFOR_t), dimension(:), intent(inout) :: retrOutputPrimaryFOR !retrOutputPrimary%FOR2D(:,iFOR)
    TYPE(OutputFOR_t), dimension(:), intent(inout) :: retrOutputSecondaryFOR !retrOutputSecondary%FOR2D(:,iFOR)

    !Local
    TYPE(MeasErr_t) :: MWmeasErr
    TYPE(MeasErr_t) :: IRmeasErr
    TYPE(ChannelSet_t) :: chanSetUnion
    integer :: iFOV
    integer :: iFOVtmp
    integer :: surfaceClassMW
    integer :: surfaceClassIR
    integer :: atmosClass
    logical :: isLandMW
    logical :: isLandIR
    TYPE(StateIndex_t) :: NR
    TYPE(StateIndex_t) :: IR
    integer :: nPar
    character (len=mxLength) :: procName

    procName = ' [RetrModule::retr]: '
    if (dbg) print *, trim(procName)//'starting ...'

    chanSetUnion = channelSetUnion(genControl,channelSelections)

    if (genControl%MWdataOn) then
       call setMWmeasErr(&
            genControl, &
            MWdefin, &
            MWob, &
            chanSetUnion, &
            MWmeasErr &
            )
    end if

!!$    do iFOV = 1, IRdefin%nFOV
    iFOV = 1 !temporary
    retrPrior%nParG = 0
    retrPrior%nEmMW = 0
    retrPrior%nEmRfIR = 0

    if (genControl%IRdataOn) then
       call setIRmeasErr(&
            IRdefin, &
            IRFOR%ob(iFOV), &
            IRmeasErr &
            )
    end if

    if (genControl%MWdataOn) then
       call MWsurfaceClassify( &
            ancillaryFOR%data(iFOV)%landFrac, &
            MWdefin, &
            MWob, &
            surfaceClassMW, &
            isLandMW &
            )
    endif

    if (genControl%IRdataOn) then
      call IRsurfaceClassify( &
           ancillaryFOR%data(iFOV)%landFrac, &
           IRdefin, &
           IRFOR%ob(iFOV), &
           surfaceClassIR, &
           isLandIR &
           )
    endif

    iFOVtmp = 1  ! Default for memory control
    if (genControl%IRdataOn) iFOVtmp = iFOV
    call atmosClassify( &
         MWdefin, &
         MWob, &
         IRdefin, &
         IRFOR%ob(iFOVtmp), &
         surfaceClassMW, &
         surfaceClassIR, &
         atmosClass &
         )

    if(genControl%statRetrOn) then
       CALL statRetr( &
            genControl, &
            statRetrControl, &
            MWdefin, &
            MWob, &
            IRdefin, &
            IRFOR, &   !should it be IRFOR%ob(iFOV)?
            ancillaryFOR, &
            retrResult &
            )
       retrPrior = retrResult
       retrOutputStatFOR(iFOV)%state = retrResult%stateV(1:size(retrOutputStatFOR(iFOV)%state))
       retrOutputStatFOR(iFOV)%press = retrResult%press
       retrOutputStatFOR(iFOV)%emMW = retrResult%emMW(1:size(retrOutputStatFOR(iFOV)%emMW))
       retrOutputStatFOR(iFOV)%emisIR = retrResult%emRfIR(1:size(retrOutputStatFOR(iFOV)%emisIR))
       retrOutputStatFOR(iFOV)%niter = retrResult%niter
       retrOutputStatFOR(iFOV)%rms = retrResult%rms
       retrOutputStatFOR(iFOV)%chiSq = retrResult%chiSq
       retrOutputStatFOR(iFOV)%lat = IRFOR%ob(iFOV)%lat
       retrOutputStatFOR(iFOV)%lon = IRFOR%ob(iFOV)%lon
       retrOutputStatFOR(iFOV)%EIA = IRFOR%ob(iFOV)%EIA
       retrOutputStatFOR(iFOV)%landFrac = ancillaryFOR%data(iFOV)%landFrac
       retrOutputStatFOR(iFOV)%landType = surfaceClassIR
    end if

    if (genControl%MWretrOn .or. genControl%primaryRetrOn .or. genControl%secondaryRetrOn) then
       if (.not. allocated(cloudClass)) allocate(cloudClass(ImgDef%nFOV))
       call cloudClassifyInitial(&
            genControl, &
            imgOb, &
            cloudClass &
            )

       if (genControl%MWretrOn) then
          call MWretr(&
               genControl, &
               MWretrControl, &
               MWdefin, &
               MWob, &
               MWmeasErr, &
               ancillaryDefin, &
               ancillaryFOR%data(iFOV), &
               channelSelections%MW, &
               atmosClass, &
               cloudClass, &
               surfaceClassMW, &
               isLandMW, &
               retrPrior, &
               retrResult)
          nSav=retrPrior%nEmRfIR
          retrSavIRsfc=retrPrior%emRfIR !retain result not updated at this stage
          retrPrior = retrResult
          retrPrior%nEmRfIR=nSav
          retrPrior%emRfIR=retrSavIRsfc
          retrOutputMWFOR%state(1:retrResult%nParG) = retrResult%stateV(1:retrResult%nParG)
          retrOutputMWFOR%press = retrResult%press
          retrOutputMWFOR%emMW(1:retrResult%nEmMW) = retrResult%emMW(1:retrResult%nEmMW)
          retrOutputMWFOR%niter = retrResult%niter
          retrOutputMWFOR%rms = retrResult%rms
          retrOutputMWFOR%chiSq = retrResult%chiSq
          retrOutputMWFOR%lat = MWob%lat
          retrOutputMWFOR%lon = MWob%lon
          retrOutputMWFOR%EIA = MWob%EIA
          retrOutputMWFOR%landFrac = ancillaryFOR%data(iFOV)%landFrac
          retrOutputMWFOR%landType = surfaceClassMW
       end if

       if (genControl%primaryRetrOn) then
          if (IRdefin%nFOV /= 1) then
             print *, trim(procName)//' error - nFOV /=1'
             call errorHalt(1)
          end if
         call primaryRetr(&
               genControl, &
               primaryRetrControl, &
               MWob, &
               MWmeasErr, &
               IRdefin, &
               IRFOR%ob(iFOV), &
               IRmeasErr, &
               ancillaryDefin, &
               ancillaryFOR%data(iFOV), &
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
          retrPrior = retrResult
          retrOutputPrimaryFOR(iFOV)%state(1:retrResult%nParG) = retrResult%stateV(1:retrResult%nParG)
          retrOutputPrimaryFOR(iFOV)%press = retrResult%press
          !- is emMW retrieved in primaryRetr?
          retrOutputPrimaryFOR(iFOV)%emMW = retrResult%emMW(1:retrResult%nEmMW)
          retrOutputPrimaryFOR(iFOV)%emisIR = retrResult%emRfIR(1:size(retrOutputPrimaryFOR(iFOV)%emisIR))
          retrOutputPrimaryFOR(iFOV)%niter = retrResult%niter
          retrOutputPrimaryFOR(iFOV)%rms = retrResult%rms
          retrOutputPrimaryFOR(iFOV)%chiSq = retrResult%chiSq
          retrOutputPrimaryFOR(iFOV)%DiagError = retrResult%DiagError
          retrOutputPrimaryFOR(iFOV)%DOFS = retrResult%DOFS
          retrOutputPrimaryFOR(iFOV)%lat = IRFOR%ob(iFOV)%lat
          retrOutputPrimaryFOR(iFOV)%lon = IRFOR%ob(iFOV)%lon
          retrOutputPrimaryFOR(iFOV)%EIA = IRFOR%ob(iFOV)%EIA
          retrOutputPrimaryFOR(iFOV)%SolIA = IRFOR%ob(iFOV)%SolIA
          retrOutputPrimaryFOR(iFOV)%landFrac = ancillaryFOR%data(iFOV)%landFrac
          retrOutputPrimaryFOR(iFOV)%landType = surfaceClassIR
          retrOutputPrimaryFOR(iFOV)%atmosClass = atmosClass

          if (primaryRetrControl%genLinInvert) then
             linearOutputFOR%xOffG(1:linearResult%nParG) = linearResult%xOffG(1:linearResult%nParG)
             linearOutputFOR%xGainG(1:linearResult%nParG,1:linearResult%nChan) = &
                  linearResult%xGainG(1:linearResult%nParG,1:linearResult%nChan)
             linearOutputFOR%chiSq = linearResult%chiSq
             linearOutputFOR%Tsfc = linearResult%Tsfc
          end if
          backgroundOutputFOR%state(1:retrResult%nParG) = backResult(1:retrResult%nParG)
       end if

       if (genControl%secondaryRetrOn) then
          if (IRdefin%nFOV /= 1) then
             print *, trim(procName)//' error - nFOV /=1'
             call errorHalt(1)
          end if
          call secondaryRetr(&
               genControl, &
               secondaryRetrControl, &
               IRdefin, &
               IRFOR%ob(iFOV), &
               IRmeasErr, &
               ancillaryDefin, &
               ancillaryFOR%data(iFOV), &
               channelSelections%secondary, &
               retrPrior, &
               retrResult)
       end if
    end if
    if (allocated(MWmeasErr%device)) deallocate(MWmeasErr%device)
    if (allocated(MWmeasErr%fwdModel)) deallocate(MWmeasErr%fwdModel)
    if (allocated(IRmeasErr%device)) deallocate(IRmeasErr%device)
    if (allocated(IRmeasErr%fwdModel)) deallocate(IRmeasErr%fwdModel)
    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine retr

  subroutine destroy(genControl,MWgran,IRgran,imgGran,ancillaryGran,ChannelSelections)
    TYPE(GenControl_t), intent(inout) :: genControl
    TYPE(MWgran_t), intent(inout) :: MWgran
    TYPE(IRgran_t), intent(inout) :: IRgran
    TYPE(ImgGran_t), intent(inout) :: imgGran
    TYPE(AncillaryGran_t), intent(inout) :: ancillaryGran
    TYPE(ChannelSelections_t), intent(inout) :: ChannelSelections
    character (len=mxLength) :: procName

    procName = "[RetrModule::destroy]:"
    if (dbg) print *, trim(procName)//"starting ..."
    if (genControl%statRetrOn) then
       call statRetrDestroy()
    end if
    if (genControl%MWretrOn) then
       call MWretrDestroy()
    end if
    if (genControl%primaryRetrOn) then
       call primaryRetrDestroy()
    end if
    if (genControl%secondaryRetrOn) then
       call secondaryRetrDestroy()
    end if
    if (genControl%MWdataOn) then
       call destroyMWsensorData(MWgran)
       call MWRTdestroy()
       call MWnoiseDestroy()
    end if
    if (genControl%IRdataOn) then
       call destroyIRsensorData(IRgran)
       call IRRTdestroy()
       call IRnoiseDestroy()
    end if
    if (genControl%imgDataOn) then
       call destroyImgSensorData(imgGran)
       call ImgRTdestroy()
    end if
    call destroyAncillaryData(ancillaryGran)
    call channelSelectionDestroy(ChannelSelections)
    call limitsDestroy()
    if (dbg) print *, trim(procName)//"deallocating ... "
    if (allocated(genControl%kchan)) deallocate(genControl%kchan)
    if (allocated(genControl%backgroundFiles)) deallocate(genControl%backgroundFiles)
    if (allocated(genControl%channelSelectFiles)) deallocate(genControl%channelSelectFiles)
    if (allocated(MolID)) deallocate(MolID)
    if (allocated(MolTran)) deallocate(MolTran)
    if (allocated(retrResult%stateV)) deallocate(retrResult%stateV)
    if (allocated(retrResult%emMW)) deallocate(retrResult%emMW)
    if (allocated(retrResult%emRfIR)) deallocate(retrResult%emRfIR)
    if (allocated(retrPrior%stateV)) deallocate(retrPrior%stateV)
    if (allocated(retrPrior%emMW)) deallocate(retrPrior%emMW)
    if (allocated(retrPrior%emRfIR)) deallocate(retrPrior%emRfIR)
    if (allocated(retrSavIRsfc)) deallocate(retrSavIRsfc)
    if (allocated(cloudClass)) deallocate(cloudClass)
    if (allocated(backResult)) deallocate(backResult)
    if (allocated(linearResult%xOffG)) deallocate(linearResult%xOffG)
    if (allocated(linearResult%xGainG)) deallocate(linearResult%xGainG)
    if (dbg) print *, trim(procName)//"ending ..."
  end subroutine destroy

  subroutine MWRTverify(MWRTcontrol,MWdefin,stateSpaceDefin, &
       nChanRT,frqRT,polRT,nLevRT,pRefRT)

    TYPE(MWRTcontrol_t), intent(in) :: MWRTcontrol
    TYPE(MWdefin_t), intent(in) :: MWdefin
    TYPE(StateSpaceDefin_t), intent(in) :: stateSpaceDefin
    real, dimension(:), intent(in) :: frqRT
    real, dimension(:), intent(in) :: pRefRT
    integer, dimension(:), intent(in) :: polRT
    integer, intent(in) :: nChanRT
    integer, intent(in) :: nLevRT

    !Local
    character (len=*), parameter :: procName = ' [RetrModule::MWRTverify]: '
    if (dbg) print *, procName, 'starting ...'

    IF (nchanRT  /= MWdefin%nChan)then
       print*,procName,'error - nChan inconsistent in MW (OSS/Bkg) ', &
            nchanRT, MWdefin%nChan
       call errorHalt(1)
    endif
    IF (ANY(abs(MWdefin%frq(1:MWdefin%nChan)-frqRT(1:MWdefin%nChan)) > deltaFrq ))then
       print*,procName,'error - MW Frequencies inconsistent (Bkg/OSS)'
       call errorHalt(1)
    endif
    IF (ANY(MWdefin%pol(1:MWdefin%nChan)  /= polRT(1:MWdefin%nChan) ))then
       print*,procName,'error - Polarization inconsistent (Bkg/OSS)'
       call errorHalt(1)
    endif
    IF (TRIM(vCoordTyp) == Pcoord) THEN
       IF (nLevRT  /= stateSpaceDefin%nLev)then
          print*,procName,'error - nLev inconsistent in MW (OSS/Bkg) ', &
               nLevRT, stateSpaceDefin%nLev
          call errorHalt(1)
       endif
       IF ( ANY(abs(stateSpaceDefin%pRef(1:stateSpaceDefin%nLev) - &
            pRefRT(1:stateSpaceDefin%nLev)) > deltaP) ) then
          print*,procName,'error - Pressure Grid inconsistent in MW (Bkg/OSS)'
          call errorHalt(1)
       endif
    END IF

    if (MWdefin%radOrTb /= MWRTcontrol%radOrTb) then
       print*,procName,'error - inconsistent radOrTb flag in MW (Bkg/OSS)', &
            MWdefin%radOrTb, MWRTcontrol%radOrTb
       call errorHalt(1)
    end if

    if (dbg) print *, procName,'ending ...'
  end subroutine MWRTverify

  subroutine IRRTverify( &
       IRdefin , &
       stateSpaceDefin, &
       nChanRT, &
       wvnRT, &
       nLevRT, &
       pRefRT &
       )

    TYPE(IRdefin_t), intent(inout) :: IRdefin
    TYPE(StateSpaceDefin_t), intent(in) :: stateSpaceDefin
    integer, intent(in) :: nChanRT
    integer, intent(in) :: nLevRT
    real, dimension(:), intent(in) :: wvnRT
    real, dimension(:), intent(in) :: pRefRT

    !Local
    character (len=mxLength) :: procName

    procName = ' [RetrModule::IRRTverify]: '

    if (dbg) print *, trim(procName)//'starting ...'
    if (vCoordTyp == Pcoord) then
       IF (nLevRT /= stateSpaceDefin%nLev)then
          print *, trim(procName)//'Nlev inconsistent in IR (OSS/Bkg) ', &
               nLevRT,stateSpaceDefin%nLev
          call errorHalt(1)
       endif

       IF (ANY(abs(pRefRT(1:nlevRT)-stateSpaceDefin%pRef(1:nlevRT)) > deltaP )) then
          print*,trim(procName)//'Pressure Grid inconsistent in IR (OSS/Bkg)'
          call errorHalt(1)
       endif
    end if

    IF (nChanRT /= IRdefin%nChan)then
       print *, trim(procName)//'nChan inconsistent in IR (OSS/Bkg) ', &
            nChanRT,IRdefin%nChan
       call errorHalt(1)
    endif

    !add verification of IRobs%wvn against wvnRT
    IF (ANY(abs(wvnRT(1:nChanRT)-IRdefin%wvn(1:nChanRT)) > deltaWvn )) then
       print*,trim(procName)//'Freq inconsistent in IR (OSS/Bkg)'
       call errorHalt(1)
    endif

    !using IRdefin%chanNum to match to the right wvnRT
    !what should be the invalid wvnRT value and where should it be set?
    where (wvnRT == 0)
       IRdefin%chanNum = 0
    end where

    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine IRRTverify

  subroutine ImgRTverify()
    !Local
    character (len=mxLength) :: procName

    procName = ' [RetrModule:ImgRTverify]: '
    if (dbg) print *, trim(procName)//'starting ...'
    if (dbg) print *, trim(procName)//'ending ...'

  end subroutine ImgRTverify

  !--
  ! Set global variables from argument list
  !--
  subroutine setPropertiesAlgConfig(ngas_in, &
       nTskin_in,ntemp_in,nemmw_in,nemir_in, &
       nCldLiq_in,nCldIce_in, &
       nxiterMw_in,nxiter_in,drad_in,alpha1slope_in,alpha1int_in, &
       chisqconvMw_in,chisqconv_in,chisqRatioconv_in,cldamt1stguess_in, &
       iceamt1stguess_in,TightCLDcov_in,plandMaxOc_in,logSize_in)
    !<f90Subroutine>********************************************************
    !
    ! NAME:
    !
    !   setPropertiesAlgConfig
    !
    ! PURPOSE:
    !
    !   Set global variables from argument list, for retrieval tuning
    !   parameters.
    !
    ! SYNTAX:
    !
    !   CALL setPropertiesAlgConfig(ngas_in, nTskin_in, ntemp_in,
    !      nemmw_in, nemir_in, nCldLiq_in, nCldIce_in, nxiterMw_in,
    !      nxiter_in, drad_in, alpha1slope_in, alpha1int_in,
    !      chisqconvMw_in, chisqconv_in, chisqRatioconv_in,
    !      cldamt1stguess_in, iceamt1stguess_in, TightCLDcov_in,
    !      plandMaxOc_in, logSize_in)
    !
    ! ARGUMENTS:
    !
    !   INPUTS:
    !
    !   ngas_in            INTEGER  Number of elements for molecular
    !                               concentrations
    !   nTskin_in          INTEGER  Number of eigenvectors for Tskin
    !   ntemp_in           INTEGER  Number of elements for temperature
    !   nemmw_in           INTEGER  Number of retrieved variables for
    !                               MW emissivity
    !   nemir_in           INTEGER  Number of hinge points for
    !                               retrieving IR emissivity
    !   nCldLiq_in         INTEGER  Number of elements for liquid cloud
    !   nCldIce_in         INTEGER  Number of elements for ice cloud
    !   nxiterMw_in        INTEGER  maximum number of iterations in MW
    !   nxiter_in          INTEGER  maximum number of iterations
    !   drad_in            REAL     Noise adjustment to stabilize
    !                               convergence, for input
    !   alpha1slope_in     REAL     Slope of factor to model
    !                               linearization error, for input
    !   alpha1int_in       REAL     Intercept of factor to model
    !                               linearization error, for input
    !   chisqconvMw_in     REAL     value of chi-squared needed for
    !                               convergence in MW
    !   chisqconv_in       REAL     value of chi-squared needed for
    !                               convergence
    !   chisqRatioconv_in  REAL     value of chisq fractional change
    !                               needed for convergence
    !   cldamt1stguess_in  REAL     First guess cloud liquid water path
    !   iceamt1stguess_in  REAL     First guess cloud ice water path
    !   TightCLDcov_in     REAL     Cloud water variance to use when
    !                               holding cloud water near zero
    !   plandMaxOc_in      REAL     Upper limit on pland for
    !                               classifiying as water surface
    !   logSize_in         LOGICAL  Flag for type of cloud size
    !                               parameter retrieval: Logarithmic vs
    !                               Linear
    !
    !   * OPTIONAL
    !
    ! INCLUDES:
    !
    !   None
    !
    !*******************************************************</f90Subroutine>

    !--I/O variables
    integer, intent(in)        :: ngas_in(maxMol)
    integer, intent(in)        :: nTskin_in,ntemp_in,nemmw_in,nemir_in
    integer, intent(in)        :: nCldLiq_in,nCldIce_in
    integer, intent(in)        :: nxiterMw_in,nxiter_in
    real,    intent(in)        :: drad_in,alpha1slope_in,alpha1int_in
    real,    intent(in)        :: chisqconvMw_in,chisqconv_in,chisqRatioconv_in
    real,    intent(in)        :: cldamt1stguess_in,iceamt1stguess_in
    real,    intent(in)        :: TightCLDcov_in
    real,    intent(in)        :: plandMaxOc_in
    logical, intent(in)        :: logSize_in
    character (len=mxLength) :: procName

    procName = ' [RetrModule::setPropertiesAlgConfig]: '
    !------------------------------------------------------------------------
    ! Set global variables from input arguments
    !------------------------------------------------------------------------
    ngas           = ngas_in
    nTskin         = nTskin_in
    ntemp          = ntemp_in
    nemmw          = nemmw_in
    nemir          = nemir_in
    nCldLiq        = nCldLiq_in
    nCldIce        = nCldIce_in
    nxiterMw       = nxiterMw_in
    nxiter         = nxiter_in
    drad           = drad_in
    alpha1slope    = alpha1slope_in
    alpha1int      = alpha1int_in
    chisqconvMw    = chisqconvMw_in
    chisqconv      = chisqconv_in
    chisqRatioconv = chisqRatioconv_in
    cldamt1stguess = cldamt1stguess_in
    iceamt1stguess = iceamt1stguess_in
    TightCLDcov    = TightCLDcov_in
    plandMaxOc     = plandMaxOc_in
    logSize        = logSize_in

    loadAlgCfgDone=.true.
    return
  end subroutine setPropertiesAlgConfig

  !--
  ! Read items from namelist into local variables
  !--
  subroutine loadAlgconfig(file_algcfg)
    !<f90Subroutine>********************************************************
    !
    ! NAME:
    !
    !   loadAlgconfig
    !
    ! PURPOSE:
    !
    !   Load retrieval tuning parameters.
    !
    ! SYNTAX:
    !
    !   CALL loadAlgconfig(file_algcfg)
    !
    ! ARGUMENTS:
    !
    !   INPUTS:
    !
    !   file_algcfg  CHAR  File path for algorithm configuration
    !                      parameters
    !
    !   * OPTIONAL
    !
    ! INCLUDES:
    !
    !   None
    !
    !*******************************************************</f90Subroutine>

    !Local parameters
    INTEGER, PARAMETER :: U_algconfig   = 30

    !--I/O variables
    character(len=*), intent(in) :: file_algcfg

    !---retrTuning Namelist items
    integer            :: ngas(maxMol)
    integer            :: nTskin=0 ! Overridden if ntemp>0 in setXPtrR
    integer            :: ntemp,nemmw,nCldLiq,nCldIce,nxiterMw,nxiter
    real               :: drad,alpha1slope,alpha1int
    real               :: chisqconvMw,chisqconv,chisqRatioconv
    real               :: cldamt1stguess,iceamt1stguess
    real               :: TightCLDcov
    logical            :: logSize
    character (len=mxLength) :: procName

    namelist /retrtuning/ nxiterMw,nxiter,drad,alpha1slope,alpha1int, &
         chisqconvMw,chisqconv,chisqRatioconv,cldamt1stguess, &
         iceamt1stguess,nTskin,ntemp,ngas,nemmw,nemir,nCldLiq, &
         nCldIce,TightCLDcov,logSize

    procName = ' [RetrModule::loadAlgconfig]: '

    if(loadAlgCfgDone) return

    !------------------------------------------------------------------------
    ! Read namelist file items into local variables
    !------------------------------------------------------------------------
    open(U_algconfig,file=file_algcfg)
    read(U_algconfig,retrtuning)
    close(U_algconfig)

    !------------------------------------------------------------------------
    ! Set global variables
    !------------------------------------------------------------------------
    call setPropertiesAlgConfig(ngas,nTskin,ntemp,nemmw,nemir,nCldLiq,nCldIce, &
         nxiterMw,nxiter,drad,alpha1slope,alpha1int,chisqconvMw, &
         chisqconv,chisqRatioconv,cldamt1stguess,iceamt1stguess,TightCLDcov, &
         plandMaxOc,logSize)

    return
  end subroutine loadAlgconfig

END MODULE RetrModule
