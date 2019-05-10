!<f90File>**************************************************************
!
! CONTACT:
!
!   Atmospheric & Environmental Research, Inc
!   131 Hartwell Ave
!   Lexington ,MA 02421-3126 USA
!   Phone: 781.761.2288
!   E-mail: guymin@aer.com
!
! COPYRIGHT NOTICE:
!
!   Copyright AER, Inc 2001-2009, All Rights Reserved
!   See the file README-DATARIGHTS.txt included with this release
!   for additional details.
!
!*************************************************************</f90File>

!--------------------------------------------------------------
!
!  MODULE ReadStdInputs.
!
!--------------------------------------------------------------
MODULE ReadStdInputs

! <f90Module>***********************************************************
!
! NAME:
!
!   ReadStdInputs
!
! PURPOSE:
!
!   Load the primary control input parameters. Control all of the hard-coded
!   file unit numbers. Set all file names. Unit file names and numbers are
!   public.
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

  USE StateIndexModule, ONLY: maxMol
  USE constants, ONLY : MISSING_CHAR
  IMPLICIT NONE
  INTEGER, PRIVATE, PARAMETER :: mxchan = 890, NCDF_DUMMY=999
  !------------------------------------------------------------------------
  !	Unit Numbers                            File Names
  !------------------------------------------------------------------------
  !---Sensor-specific files: oss opt depth, oss coefs, nedt, noise map, nrf, etc.
  INTEGER,PARAMETER :: U_ossoptdpth  = 10;  CHARACTER(LEN=200) :: F_ossoptdpth
  INTEGER,PARAMETER :: U_osscoefs    = 11;  CHARACTER(LEN=200) :: F_osscoefs
  INTEGER,PARAMETER :: U_nedt        = 12;  CHARACTER(LEN=200) :: F_nedt
  INTEGER,PARAMETER :: U_nrf         = 14;  CHARACTER(LEN=200) :: F_nrf
  INTEGER,PARAMETER :: U_Classatm    = 15;  CHARACTER(LEN=200) :: F_Classatm
  INTEGER,PARAMETER :: U_Classsfc    = 16;  CHARACTER(LEN=200) :: F_Classsfc
  INTEGER,PARAMETER :: U_beam        = 17;  CHARACTER(LEN=200) :: F_beam
  INTEGER,PARAMETER :: U_spov        = 18;  CHARACTER(LEN=200) :: F_spov
  INTEGER,PARAMETER :: U_emChange    = 19;  CHARACTER(LEN=200) :: F_emChange
  !---Scene, Rad and auxill. files
  INTEGER :: U_scene         = NCDF_DUMMY;  CHARACTER(LEN=200) :: F_scene
  INTEGER :: U_Rad           = NCDF_DUMMY;  CHARACTER(LEN=200) :: F_Rad
  INTEGER :: U_Aux           = NCDF_DUMMY;  CHARACTER(LEN=200) :: F_Aux
  !---Retrieval files, used exclusively by Retr
  INTEGER :: U_retr          = NCDF_DUMMY;  CHARACTER(LEN=200) :: F_Retr
  INTEGER :: U_guess         = NCDF_DUMMY;  CHARACTER(LEN=200) :: F_Guess
  INTEGER :: U_cov           = NCDF_DUMMY;  CHARACTER(LEN=200),DIMENSION(2) :: F_covar
  INTEGER,PARAMETER :: U_algconfig   = 30;  CHARACTER(LEN=200) :: F_algconfig
  INTEGER,PARAMETER :: U_cascIn      = 31;  CHARACTER(LEN=200) :: F_cascin
  INTEGER,PARAMETER :: U_cascOut     = 32;  CHARACTER(LEN=200) :: F_cascout
  INTEGER,PARAMETER :: U_SpecEmLUT   = 33;  CHARACTER(LEN=200) :: F_SpecEmLUT
  INTEGER,PARAMETER :: U_QC          = 34;  CHARACTER(LEN=200) :: F_QC
  !---External files
  INTEGER,PARAMETER :: U_cloudTab    = 40;  CHARACTER(LEN=200) :: F_cloudTab
  INTEGER :: U_ldEmissExt    = NCDF_DUMMY;  CHARACTER(LEN=200) :: F_ldEmissExt
  INTEGER :: U_bkgrdExtern   = NCDF_DUMMY;  CHARACTER(LEN=200) :: F_bkgrdExtern
  INTEGER,PARAMETER :: U_namelist    = 70;  CHARACTER(LEN=200) :: F_namelist

CONTAINS

  SUBROUTINE readStd(nprofs,debugIn,mwCldIn,iCellIn,icascIn, &
       iembkgflgIn,iextbkgIn,kchanOut,DynamicNoiseIn,addNoiseIn, &
       nRetrAttemptsIn,spilloverNoiseIn,calibAmpIn,&
       enableScatterIn,scatterSizeOptIn,finiteDiffOnIn, &
       use2ndOrderPlanckIn,pobsIn, &
       extFlgTempOut,extFlgMolOut,extFlgPsfcOut, &
       extFlgCldLiqOut,extFlgCldIceOut,extFlgWindOut, &
       molOnRTin, &
       nlfile)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   readStd
!
! PURPOSE:
!
!   Load the primary control input parameters.
!
! SYNTAX:
!
!   CALL readStd(nprofs, debug, mwCld, iCellIn, icasc, iembkgflg,
!      iextbkg, kchanOut, DynamicNoise, addNoiseIn, nRetrAttemptsIn,
!      spilloverNoiseIn, calibAmpIn, enableScatterIn, scatterSizeOptIn,
!      finiteDiffOnIn,
!      use2ndOrderPlanckIn, pobsIn, extFlgTempOut, extFlgMolOut,
!      extFlgPsfcOut, extFlgCldLiqOut, extFlgCldIceOut,
!      extFlgWindOut, molOnRTin, nlfile)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   nlfile*              CHAR     File path from which to read the namelist
!
!   INPUTS/OUTPUTS:
!
!   nprofs               INTEGER  Number profiles
!   debug                LOGICAL  Flag for debugging mode
!   mwCld                INTEGER  Interger flag to include clouds in MW
!                                 (inactive)
!   iCellIn              INTEGER  Index for horizontal cell size (resolution
!                                 of processing)
!   icasc                INTEGER  Index for cascade level (when in cascade
!                                 mode)
!   iembkgflg            INTEGER  Flag to use external emissivity background
!                                 in MW
!   iextbkg              INTEGER  Flag to use external atmosphere background
!                                 in MW
!   kchanOut             logical  Channel on/off mask, for output
!   DynamicNoise         LOGICAL  Flag for using dynamic noise model
!   addNoiseIn*          LOGICAL  Flag to add noise, in simulation
!   nRetrAttemptsIn*     INTEGER  maximum number of retrieval attempts, with
!                                 different settings
!   spilloverNoiseIn*    LOGICAL  Include antenna spillover factor in
!                                 microwave noise computation
!   calibAmpIn*          LOGICAL  Include calibration amplification factor in
!                                 microwave noise computation
!   enableScatterIn*     LOGICAL  Flag for including multiple scattering in
!                                 calculations, for input
!   scatterSizeOptIn*    INTEGER  Index for how size parameter is specified
!                                 for multiple scattering, for input
!   finiteDiffOnIn*      LOGICAL  Finite differencing flag, for input
!   use2ndOrderPlanckIn* LOGICAL  flag to use MW Planck approximation from
!                                 2nd-order Taylor expansion
!   pobsIn*              REAL     Pressure at observer level
!   extFlgTempOut*       LOGICAL  Flag for using external atmospheric
!                                 temperature
!   extFlgMolOut*        LOGICAL  Flag for using external set of molecular
!                                 concentrations
!   extFlgPsfcOut*       LOGICAL  Flag for using external surface pressure
!   extFlgCldLiqOut*     LOGICAL  Flag for using external data on liquid
!                                 clouds
!   extFlgCldIceOut*     LOGICAL  Flag for using external data on ice clouds
!   extFlgWindOut*       LOGICAL  Flag for using external wind data
!   molOnRTin*           CHAR     List of molecular species turned on for
!                                 radiative transfer
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !---In/Out variables
    INTEGER,                        intent(inout)          :: nprofs
    LOGICAL,                        intent(inout), optional :: debugIn
    INTEGER,                        intent(inout), optional :: mwCldIn,iCellIn
    INTEGER,                        intent(inout), optional :: icascIn,iembkgflgIn
    INTEGER,                        intent(inout), optional :: iextbkgIn
    LOGICAL,          DIMENSION(:), intent(inout), optional :: kchanOut
    LOGICAL,                        intent(inout), optional :: DynamicNoiseIn
    LOGICAL,                        intent(inout),optional :: addNoiseIN
    LOGICAL,                        intent(inout),optional :: spilloverNoiseIn
    LOGICAL,                        intent(inout),optional :: calibAmpIn
    LOGICAL,                        intent(inout),optional :: enableScatterIn
    INTEGER,                        intent(inout),optional :: scatterSizeOptIn
    LOGICAL,                        intent(inout),optional :: finiteDiffOnIn
    LOGICAL,                        intent(inout),optional :: use2ndOrderPlanckIn
    LOGICAL,                        intent(inout),optional :: extFlgTempOut
    LOGICAL,                        intent(inout),optional :: extFlgPsfcOut
    LOGICAL,                        intent(inout),optional :: extFlgWindOut
    LOGICAL,                        intent(inout),optional :: extFlgCldLiqOut
    LOGICAL,                        intent(inout),optional :: extFlgCldIceOut
    LOGICAL,          DIMENSION(:), intent(inout),optional :: extFlgMolOut
    CHARACTER(LEN=*), DIMENSION(:), intent(inout),optional :: molOnRTin
    REAL,                           intent(inout),optional :: pobsIn
    INTEGER,                        intent(inout),optional :: nRetrAttemptsIN
    CHARACTER(LEN=*),               intent(in),   optional :: nlfile
    !---Local variables
    LOGICAL                    :: addnoise=.TRUE.
    LOGICAL                    :: spilloverNoise=.FALSE.
    LOGICAL                    :: calibAmp=.FALSE.
    LOGICAL                    :: enableScatter=.FALSE.
    INTEGER                    :: scatterSizeOpt=0
    LOGICAL                    :: finiteDiffOn=.FALSE.
    LOGICAL                    :: extFlgTemp=.FALSE.
    LOGICAL                    :: extFlgPsfc=.FALSE.
    LOGICAL                    :: extFlgWind=.FALSE.
    LOGICAL, DIMENSION(maxMol) :: extFlgMol=.FALSE.
    LOGICAL                    :: extFlgCldLiq=.FALSE.,extFlgCldIce=.FALSE.
    INTEGER                    :: nmol
    LOGICAL                    :: use2ndOrderPlanck=.FALSE.
    REAL                       :: pobs=0.0
    INTEGER                    :: nRetrAttempts=2
    INTEGER, DIMENSION(mxchan) :: kchan=9
    INTEGER                    :: nchmwIn,icell=1
    INTEGER                    :: mwCld, icasc, iembkgflg, iextbkg
    LOGICAL                    :: debug, DynamicNoise
    CHARACTER(LEN=12), DIMENSION(maxMol) :: molOnRT
    CHARACTER(LEN=200) :: AuxFile,RadFile,cloudTabfile,ldEmissExtern,bkgrdExtern
    CHARACTER(LEN=200) :: Scenefile,nedtFile,noiseReductFactFile='NULL', &
         osscoefs,ossoptdpth,coefsClassatm,coefsClasssfc,algconfig,QCFile, &
         RetrFile,GuessFile,cascinFIle,cascoutFile,SpecEmLUTfile
    CHARACTER(LEN=200), DIMENSION(2)  :: covarFile
    NAMELIST /stdContrl/ nprofs,debug,mwCld,iCell,icasc,iembkgflg,iextbkg,&
         DynamicNoise,enableScatter,scatterSizeOpt,pobs,kchan,AuxFile,RadFile,&
         use2ndOrderPlanck,molOnRT,finiteDiffOn
    NAMELIST /externFiles/ cloudTabfile,ldEmissExtern,bkgrdExtern,&
         extFlgTemp,extFlgPsfc,extFlgWind,extFlgCldLiq,extFlgCldIce,&
         extFlgMol
    NAMELIST /SimFiles/ Scenefile, addNoise
    NAMELIST /sensorSpecificFiles/ nedtFile,noiseReductFactFile,osscoefs,  &
         ossoptdpth,coefsClassatm,coefsClasssfc,F_beam,F_spov,spilloverNoise,calibAmp, &
         F_emChange
    NAMELIST /RetrFiles/ algconfig,covarFile,QCFile,RetrFile,         &
         GuessFile,cascinFIle,cascoutFile,SpecEmLUTfile,nRetrAttempts

    ! The namelist may come from standard input (default) or a specified file (optional).

    if (present(nlfile)) then
       open(U_namelist,file=F_namelist)
    endif

    molOnRT=MISSING_CHAR  ! all elements not given are marked as unused
    if (present(nlfile)) then
       read(U_namelist,stdContrl)
    else
       READ(*,stdContrl)
    endif
    if(present(pobsIn))pobsIn=pobs

    if (present(nlfile)) then
       read(U_namelist,SimFiles)
    else
       READ(*,SimFiles)
    endif
    if(present(addnoiseIn))addnoiseIn=addnoise
    if(present(enableScatterIn))enableScatterIn=enableScatter
    if(present(scatterSizeOptIn))scatterSizeOptIn=scatterSizeOpt
    if(present(finiteDiffOnIn))finiteDiffOnIn=finiteDiffOn

    if (present(nlfile)) then
       read(U_namelist,sensorSpecificFiles)
    else
       READ(*,sensorSpecificFiles)
    endif
    if(present(spilloverNoiseIn))spilloverNoiseIn=spilloverNoise
    if(present(calibAmpIn))calibAmpIn=calibAmp
    if(present(use2ndOrderPlanckIn))use2ndOrderPlanckIn=use2ndOrderPlanck

    if (present(nlfile)) then
       read(U_namelist,externFiles)
    else
       READ(*,externFiles)
    endif

    if (present(nlfile)) then
       read(U_namelist,RetrFiles)
    else
       READ(*,RetrFiles)
    endif

    if (present(nlfile)) then
       close(U_namelist)
    endif
    if (present(extFlgTempOut))extFlgTempOut=extFlgTemp
    if (present(extFlgPsfcOut))extFlgPsfcOut=extFlgPsfc
    if (present(extFlgCldLiqOut))extFlgCldLiqOut=extFlgCldLiq
    if (present(extFlgCldIceOut))extFlgCldIceOut=extFlgCldIce
    if (present(extFlgWindOut))extFlgWindOut=extFlgWind
    if (present(extFlgMolOut)) extFlgMolOut(:)=extFlgMol(1:size(extFlgMolOut))
    if (present(molOnRTin)) molOnRTin(:)=molOnRT(1:size(molOnRTin))
    if(present(nRetrAttemptsIn))nRetrAttemptsIn=nRetrAttempts
    nchmwIn=COUNT(kchan .NE. 9)
    if (present(kchanOut)) kchanOut(1:nchmwIn)=kchan(1:nchmwIn)/=0
    if (present(icellIn)) icellIn = icell
    if (present(debugIn)) debugIn = debug
    if (present(mwCldIn)) mwCldIn = mwCld
    if (present(iCellIn)) iCellIn = iCell
    if (present(icascIn)) icascIn = icasc
    if (present(iembkgflgIn)) iembkgflgIn = iembkgflg
    if (present(iextbkgIn)) iextbkgIn = iextbkg
    if (present(DynamicNoiseIn)) DynamicNoiseIn = DynamicNoise

    !----Patch could be removed if we change names in namelist file AND script
    F_Aux         = AuxFile
    F_Rad         = RadFile
    F_Scene       = Scenefile
    F_nedt        = nedtFile
    F_nrf         = noiseReductFactFile
    F_osscoefs    = osscoefs
    F_ossoptdpth  = ossoptdpth
    F_Classatm    = coefsClassatm
    F_Classsfc    = coefsClasssfc
    F_cloudTab    = cloudTabfile
    F_ldEmissExt  = ldEmissExtern
    F_bkgrdExtern = bkgrdExtern
    F_algconfig   = algconfig
    F_covar(1:2)  = covarFile(1:2)
    F_QC          = QCFile
    F_Retr        = RetrFile
    F_Guess       = GuessFile
    F_cascin      = cascinFIle
    F_cascout     = cascoutFile
    F_SpecEmLUT   = SpecEmLUTfile
  END SUBROUTINE readStd

  !Return the global path(s) upon request
  SUBROUTINE getGlobalPath(P_ossoptdpth,P_osscoefs,P_nedt,P_nrf, &
       P_Classatm,P_Classsfc,P_beam,P_spov,P_emChange,P_scene, &
       P_Rad,P_Aux,P_Retr,P_Guess,P_covar,P_algconfig,P_cascin, &
       P_cascout,P_SpecEmLUT,P_QC,P_cloudTab,P_ldEmissExt, &
       P_bkgrdExtern,P_namelist)
    character(len=*), intent(inout), optional :: P_ossoptdpth
    character(len=*), intent(inout), optional :: P_osscoefs
    character(len=*), intent(inout), optional :: P_nedt
    character(len=*), intent(inout), optional :: P_nrf
    character(len=*), intent(inout), optional :: P_Classatm
    character(len=*), intent(inout), optional :: P_Classsfc
    character(len=*), intent(inout), optional :: P_beam
    character(len=*), intent(inout), optional :: P_spov
    character(len=*), intent(inout), optional :: P_emChange
    character(len=*), intent(inout), optional :: P_scene
    character(len=*), intent(inout), optional :: P_Rad
    character(len=*), intent(inout), optional :: P_Aux
    character(len=*), intent(inout), optional :: P_Retr
    character(len=*), intent(inout), optional :: P_Guess
    character(len=*), dimension(2), intent(inout), optional :: P_covar
    character(len=*), intent(inout), optional :: P_algconfig
    character(len=*), intent(inout), optional :: P_cascin
    character(len=*), intent(inout), optional :: P_cascout
    character(len=*), intent(inout), optional :: P_SpecEmLUT
    character(len=*), intent(inout), optional :: P_QC
    character(len=*), intent(inout), optional :: P_cloudTab
    character(len=*), intent(inout), optional :: P_ldEmissExt
    character(len=*), intent(inout), optional :: P_bkgrdExtern
    character(len=*), intent(inout), optional :: P_namelist

    if (present(P_ossoptdpth)) P_ossoptdpth = F_ossoptdpth
    if (present(P_osscoefs)) P_osscoefs = F_osscoefs
    if (present(P_nedt)) P_nedt = F_nedt
    if (present(P_nrf)) P_nrf = F_nrf
    if (present(P_Classatm)) P_Classatm = F_Classatm
    if (present(P_Classsfc)) P_Classsfc = F_Classsfc
    if (present(P_beam)) P_beam = F_beam
    if (present(P_spov)) P_spov = F_spov
    if (present(P_emChange)) P_emChange = F_emChange
    if (present(P_scene)) P_scene = F_scene
    if (present(P_Rad)) P_Rad = F_Rad
    if (present(P_Aux)) P_Aux = F_Aux
    if (present(P_Retr)) P_Retr = F_Retr
    if (present(P_Guess)) P_Guess = F_Guess
    if (present(P_covar)) P_covar = F_covar
    if (present(P_algconfig)) P_algconfig = F_algconfig
    if (present(P_cascin)) P_cascin = F_cascin
    if (present(P_cascout)) P_cascout = F_cascout
    if (present(P_SpecEmLUT)) P_SpecEmLUT = F_SpecEmLUT
    if (present(P_QC)) P_QC = F_QC
    if (present(P_cloudTab)) P_cloudTab = F_cloudTab
    if (present(P_ldEmissExt)) P_ldEmissExt = F_ldEmissExt
    if (present(P_bkgrdExtern)) P_bkgrdExtern = F_bkgrdExtern
    if (present(P_namelist)) P_namelist = F_namelist
  END SUBROUTINE getGlobalPath

END MODULE ReadStdInputs



