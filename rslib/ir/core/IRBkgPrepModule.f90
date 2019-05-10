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
!  MODULE BkgPrepModule: contains procedures related to processing
!          retrieval background data (a priori and error covariance)
!
!--------------------------------------------------------------
MODULE IRBkgPrepModule

! <f90Module>***********************************************************
!
! NAME:
!
!   IRBkgPrepModule
!
! PURPOSE:
!
!   Processing retrieval background data (a priori and error covariance)
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

  USE StateIndexModule
  USE constants, ONLY: MISSING_REAL
  USE VertCoord, ONLY: mxCoordTyp,Pcoord,Scoord,Hcoord, &
     putSigmDefin,putHybrDefin,getSigmPres,getHybrPres
  USE ToolboxModule, Only: getUnit

  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------
  !  List of Public subroutines (accessible from outside module) 
  !---------------------------------------------------------------
  PUBLIC :: IRbkgPrepInit,loadBkgClimo, &
            IRsetFOVbkg,IRcheckTransformConform,IRcombBkg,IRtuneBkg, &
            IRtransBkg,IRsetRetrVec,&
            IRsetExtBkg,IRsetEmisBkg,validateExt,validateEmis,cloudModeTest
  !---------------------------------------------------------------
  !  List of Public structures (accessible from outside module) 
  !  The structure retFlags_t is used to indicate whether or not
  !  an element such as temp, cldLiq, etc. should be retrieved.
  !  This structure of flags is intended to be passed in as a
  !  parameter to the subroutine IRsetRetrVec.
  !---------------------------------------------------------------
  PUBLIC :: retFlags_t, RET_FLAGS_ALL_ON
  PUBLIC :: extFlags_t
  
  TYPE retFlags_t
    LOGICAL                     :: tempFlag, tskinFlag, cldLiqFlag, cldIceFlag
    LOGICAL                     :: emissMwFlag, emissIrFlag, psfcFlag
    LOGICAL                     :: iceCldThkFlag, liqCldThkFlag, windFlag    
    LOGICAL, Dimension(1:maxMol):: molFlags
  END TYPE retFlags_t

  TYPE(retFlags_t), PARAMETER  :: RET_FLAGS_ALL_ON = retFlags_t(.TRUE.,    &
                                   .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., &
                                   .TRUE., .TRUE., .TRUE., .TRUE., .TRUE.)

  TYPE extFlags_t
    LOGICAL                     :: xTemp, xTskin
    LOGICAL                     :: xCldLiq, xCldIce
    LOGICAL                     :: xPsfc, xWind    
    LOGICAL, Dimension(1:maxMol):: xMols
  END TYPE extFlags_t

  TYPE(extFlags_t), PARAMETER  :: EXT_FLAGS_ALL_OFF = extFlags_t(.FALSE.,    &
                                   .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &
                                   .FALSE.)

  !---------------------------------------------------------------
  !     Declarations of private data (private to the module)
  !---------------------------------------------------------------
  INTEGER, DIMENSION(:), allocatable      :: myPres
  integer                                 :: myNemIrG
  INTEGER, DIMENSION(:), allocatable      :: AtmTypes,SfcTypes
  INTEGER, DIMENSION(:), allocatable      :: extrvec
  INTEGER                                 :: nSfcTypes,nAtmTypes
  INTEGER                                 :: nPar,nparG
  INTEGER                                 :: nparGatm
  INTEGER                                 :: nparGsfc
  INTEGER                                 :: nparAtm
  INTEGER                                 :: nparSfc
  INTEGER                                 :: nemmw,iemmw,nemmwG,iemmwG
  INTEGER                                 :: nemir
  INTEGER                                 :: nmol,iemIRG,ih2o
  INTEGER                                 :: nmolExt
  INTEGER, DIMENSION(maxMol)              :: mapExtMols
  TYPE(extFlags_t)                        :: externalFlags
  LOGICAL                                 :: trnFixedAtm
  LOGICAL                                 :: trnFixedSfc
  !--- my... are inserted as intermediate variables to avoid confusing
  !--- the pgf90 optimizer and getting memory violations
  real,    dimension(:), allocatable      :: myFreq,myFrqEmIR
  integer, dimension(:), allocatable      :: myPol
  REAL, DIMENSION(:,:), ALLOCATABLE       :: tmpl

  !---BackgTuning Namelist items
  REAL :: backcldtopoc,backcldtopld,backcldthkoc,backcldthkld
  REAL :: backcldamtoc,backcldamtld,varcldtopoc,varcldtopld
  REAL :: varcldthkoc,varcldthkld,varcldamtoc,varcldamtld
  REAL :: backicetopoc,backicetopld,backicethkoc,backicethkld
  REAL :: backiceamtoc,backiceamtld,varicetopoc,varicetopld
  REAL :: varicethkoc,varicethkld,variceamtoc,variceamtld
  REAL :: backclddeffoc,backclddeffld
  REAL :: backicedeffoc,backicedeffld
  REAL :: varclddeffoc,varclddeffld
  REAL :: varicedeffoc,varicedeffld
  REAL :: varTskinOc,varTskinLd

  ! The transformation matrix    
  REAL, dimension(:,:), allocatable      :: umtx_clim
  REAL, dimension(:,:), allocatable      :: vmtx_clim


  !---Climo Background Matrix (atm and sfc)
  REAL, DIMENSION(:,:), allocatable       :: back_clim_atm_all,back_clim_sfc_all
  !---Climo Transformation Matrices (atm and sfc)
  REAL, DIMENSION(:,:,:), allocatable     :: u_clim_atm_all,u_clim_sfc_all
  REAL, DIMENSION(:,:,:), allocatable     :: v_clim_atm_all,v_clim_sfc_all
  !---Climo Covariance Matrix (atm and sfc)
  REAL, DIMENSION(:,:,:), allocatable     :: ut_cov_clim_atm_all,ut_cov_clim_sfc_all
  !---Total Background & Covariance Matrix / Selected Climo - External
  REAL, DIMENSION(:)  , allocatable       :: xbakg_clim,xbakg_emissdb
  REAL, DIMENSION(:,:), allocatable       :: ut_cov_clim, ut_cov_ext,ut_cov_emissdb
  !---External Data Background Vector:
  REAL, DIMENSION(:),   allocatable       :: xExt
  REAL, DIMENSION(:),   allocatable       :: emRfExt
  REAL, DIMENSION(:,:), allocatable       :: cvtIR
  !---Background Matrix / Selected Climo - External
  !---Indices
  TYPE(StateIndex_t)             :: IG,NG,IC,NC
  INTEGER                        :: nTranMW,nTranIR

  ! loadDone implementation provides the option to call load*() in advance
  logical                            :: loadBkgTuningDone=.false.
  logical                            :: loadBkgClimoDone=.false.
  logical :: extAtmosOn = .false.
  logical :: extEmRfOn = .false.

CONTAINS

  !--------------------------------------------------------------
  ! Get the necessary parameters and allocate memory for vectors 
  ! such as the background vectors and covariance matrices in
  ! geophysical space in preparation for retreival of IR emissivities.
  !--------------------------------------------------------------
  SUBROUTINE IRbkgPrepInit(nparGout,IGout,NGout,nlevel, &
       vCoordTyp,pref,nchmw,freq,polarity,nemIRG,frqEmIR, &
       iemmwGout,iEmIRGout,molID,molTran,limitCfg,bkgCfg)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   IRbkgPrepInit
!
! PURPOSE:
!
!   Load and initialize data, and allocate memory
!
! SYNTAX:
!
!   CALL IRbkgPrepInit(nparGout, IGout, NGout, nlevel, vCoordTyp, 
!      pref, 
!      nchmw, freq, polarity, nemIRG, frqEmIR, iemmwGout, iEmIRGout, 
!      molID,molTran)
!
! ARGUMENTS:
!
!   INPUTS/OUTPUTS:
!   
!   nparGout   INTEGER             Total number of elements in geophysical 
!                                  state vector for output 
!   IGout      TYPE(STATEINDEX_T)  Starting indices for sections of 
!                                  geophysical state vector for output 
!   NGout      TYPE(STATEINDEX_T)  Number of elements for sections of 
!                                  geophysical state vector for output 
!   nlevel     INTEGER             Number of atmospheric levels
!   pref       REAL                Standard pressure profile
!   nchmw      INTEGER             Number of MW channels
!   freq       REAL                Frequency
!   polarity   INTEGER             Polarization of MW channels
!   nemIRG     INTEGER             Number of IR emissivity hinge points
!   frqEmIR    REAL                Wavenumbers at IR emissivity hinge points
!   iemmwGout  INTEGER             Starting index for MW emissivity in 
!                                  geophysical state vector for output 
!   iEmIRGout  INTEGER             Starting index for IR emissivity hinge 
!                                  points for output 
!   molID      INTEGER             List of IDs of relevant molecular species
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !--I/O variables
    integer,               intent(inout) :: nparGout
    type(StateIndex_t),    intent(inout) :: IGout,NGout  !---Indices
    integer,               intent(inout) :: nLevel
    CHARACTER(LEN=mxCoordTyp),intent(inout) :: vCoordTyp
    real,    dimension(:), allocatable       :: pref
    integer,               intent(inout) :: nchmw
    real,    dimension(:), allocatable   :: freq
    integer, dimension(:), allocatable   :: polarity
    integer,               intent(inout) :: nemIRG 
    real,    dimension(:), allocatable   :: frqEmIR
    integer,               intent(inout) :: iemmwGout,iEmIRGout
    integer, dimension(:), intent(inout) :: molID
    integer, dimension(:), intent(inout) :: molTran
    character (len=*), intent(in) :: limitCfg
    character (len=*), dimension(:), intent(in) :: bkgCfg

    ih2o=whereH2O(molID)
    !------------------------------------------------------------------------
    ! Set global variables from file variables
    !------------------------------------------------------------------------

    ! Load bkgtuning parameters
    call loadBkgTuning(limitCfg)

    !  Load Atmospheric covariance (assumed to be the 1st in covarFile(1:2))
    !  Load Surface covariance     (assumed to be the 2nd in covarFile(1:2))
    call loadBkgClimo(bkgCfg(1),bkgCfg(2),molID,molTran, &
        nlevel,vCoordTyp,pref)

    !------------------------------------------------------------------------
    ! Compute global variables from file and input argument values
    !------------------------------------------------------------------------
    nParG = nParGatm + npargsfc
    iemmwG        = IG%wind+NG%wind
    IEmIrG        = iemmwG+NemmwG ! Starting index for IR emissivity

    !------------------------------------------------------------------------
    ! Copy to output arguments
    !------------------------------------------------------------------------
    allocate(freq(NemmwG),polarity(NemmwG),frqEmIR(myNemIrG))
    nparGout      = nparG
    IGout         = IG
    NGout         = NG
    nchmw         = NemmwG
    freq          = myFreq
    polarity      = myPol
    nemIRG        = myNemIrG
    frqEmIR       = myFrqEmIR
    iemmwGout     = iemmwG
    iEmIRGout     = IEmIrG
 
    !------------------------------------------------------------------------
    ! Allocate workspace. 
    !------------------------------------------------------------------------
    allocate(xbakg_clim(NParG),xbakg_emissdb(NParG))
    allocate(ut_cov_clim(NParG,NParG),ut_cov_ext(NParG,NParG),ut_cov_emissdb(NParG,NParG))
    ALLOCATE(umtx_clim(nParAtm+nParSfc,nParG))
    ALLOCATE(vmtx_clim(nParG,nParAtm+nParSfc))
    ALLOCATE(extrvec(nParG))
    ALLOCATE(tmpl(nParG, nParG))

    RETURN
  END SUBROUTINE IRbkgPrepInit

  !--
  ! Set global variables from argument list
  !--
  subroutine setPropertiesBkgTuning(backcldtopoc_in,backcldtopld_in, &
       backcldthkoc_in,backcldthkld_in,backcldamtoc_in,backcldamtld_in, &
       backclddeffoc_in,backclddeffld_in, &
       varcldtopoc_in,varcldtopld_in,varcldthkoc_in,varcldthkld_in, &
       varcldamtoc_in,varcldamtld_in,varclddeffoc_in,varclddeffld_in, &
       backicetopoc_in,backicetopld_in, &
       backicethkoc_in,backicethkld_in,backiceamtoc_in,backiceamtld_in, &
       backicedeffoc_in,backicedeffld_in, &
       varicetopoc_in,varicetopld_in,varicethkoc_in,varicethkld_in, &
       variceamtoc_in,variceamtld_in,varicedeffoc_in,varicedeffld_in, &
       varTskinOc_in,varTskinLd_in)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   setPropertiesBkgTuning
!
! PURPOSE:
!
!   Set global variables from argument list, for background tuning parameters
!
! SYNTAX:
!
!   CALL setPropertiesBkgTuning(backcldtopoc_in, backcldtopld_in, 
!      backcldthkoc_in, backcldthkld_in, backcldamtoc_in, 
!      backcldamtld_in, backclddeffoc_in, backclddeffld_in, 
!      varcldtopoc_in, varcldtopld_in, varcldthkoc_in, 
!      varcldthkld_in, varcldamtoc_in, varcldamtld_in, 
!      varclddeffoc_in, varclddeffld_in, backicetopoc_in, 
!      backicetopld_in, backicethkoc_in, backicethkld_in, 
!      backiceamtoc_in, backiceamtld_in, backicedeffoc_in, 
!      backicedeffld_in, varicetopoc_in, varicetopld_in, 
!      varicethkoc_in, varicethkld_in, variceamtoc_in, 
!      variceamtld_in, varicedeffoc_in, varicedeffld_in, 
!      varTskinOc_in, varTskinLd_in)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   backcldtopoc_in   REAL  Background override for liquid cloud top - ocean
!   backcldtopld_in   REAL  Background override for liquid cloud top - land
!   backcldthkoc_in   REAL  Background override for liquid cloud thickness - 
!                           ocean 
!   backcldthkld_in   REAL  Background override for liquid cloud thickness - 
!                           land 
!   backcldamtoc_in   REAL  Background override for liquid cloud water - 
!                           ocean 
!   backcldamtld_in   REAL  Background override for liquid cloud water - land
!   backclddeffoc_in  REAL  Background override for liquid cloud effective 
!                           diameter - ocean 
!   backclddeffld_in  REAL  Background override for liquid cloud effective 
!                           diameter - land 
!   varcldtopoc_in    REAL  Error variance override for liquid cloud top - 
!                           ocean 
!   varcldtopld_in    REAL  Error variance override for liquid cloud top - 
!                           land 
!   varcldthkoc_in    REAL  Error variance override for liquid cloud 
!                           thickness - ocean 
!   varcldthkld_in    REAL  Error variance override for liquid cloud 
!                           thickness - land 
!   varcldamtoc_in    REAL  Error variance override for liquid cloud water - 
!                           ocean 
!   varcldamtld_in    REAL  Error variance override for liquid cloud water - 
!                           land 
!   varclddeffoc_in   REAL  Error variance override for liquid cloud 
!                           effective diameter - ocean 
!   varclddeffld_in   REAL  Error variance override for liquid cloud 
!                           effective diameter - land 
!   backicetopoc_in   REAL  Background override for ice cloud top - ocean
!   backicetopld_in   REAL  Background override for ice cloud top - land
!   backicethkoc_in   REAL  Background override for ice cloud thickness - 
!                           ocean 
!   backicethkld_in   REAL  Background override for ice cloud thickness - 
!                           land 
!   backiceamtoc_in   REAL  Background override for ice cloud water - ocean
!   backiceamtld_in   REAL  Background override for ice cloud water - land
!   backicedeffoc_in  REAL  Background override for ice cloud effective 
!                           diameter - ocean 
!   backicedeffld_in  REAL  Background override for ice cloud effective 
!                           diameter - land 
!   varicetopoc_in    REAL  Error variance override for ice cloud top - ocean
!   varicetopld_in    REAL  Error variance override for ice cloud top - land
!   varicethkoc_in    REAL  Error variance override for ice cloud thickness - 
!                           ocean 
!   varicethkld_in    REAL  Error variance override for ice cloud thickness - 
!                           land 
!   variceamtoc_in    REAL  Error variance override for ice cloud water - 
!                           ocean 
!   variceamtld_in    REAL  Error variance override for ice cloud water - 
!                           land 
!   varicedeffoc_in   REAL  Error variance override for ice cloud effective 
!                           diameter - ocean 
!   varicedeffld_in   REAL  Error variance override for ice cloud effective 
!                           diameter - land 
!   varTskinOc_in     REAL  Error variance override for Tskin - ocean
!   varTskinLd_in     REAL  Error variance override for Tskin - land
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !--I/O variables
    real,  intent(in) ::  backcldtopoc_in
    real,  intent(in) ::  backcldtopld_in
    real,  intent(in) ::  backcldthkoc_in
    real,  intent(in) ::  backcldthkld_in
    real,  intent(in) ::  backcldamtoc_in
    real,  intent(in) ::  backcldamtld_in
    real,  intent(in) ::  backclddeffoc_in
    real,  intent(in) ::  backclddeffld_in
    real,  intent(in) ::  varcldtopoc_in
    real,  intent(in) ::  varcldtopld_in
    real,  intent(in) ::  varcldthkoc_in
    real,  intent(in) ::  varcldthkld_in
    real,  intent(in) ::  varcldamtoc_in
    real,  intent(in) ::  varcldamtld_in
    real,  intent(in) ::  varclddeffoc_in
    real,  intent(in) ::  varclddeffld_in
    real,  intent(in) ::  backicetopoc_in
    real,  intent(in) ::  backicetopld_in
    real,  intent(in) ::  backicethkoc_in
    real,  intent(in) ::  backicethkld_in
    real,  intent(in) ::  backiceamtoc_in
    real,  intent(in) ::  backiceamtld_in
    real,  intent(in) ::  backicedeffoc_in
    real,  intent(in) ::  backicedeffld_in
    real,  intent(in) ::  varicetopoc_in
    real,  intent(in) ::  varicetopld_in
    real,  intent(in) ::  varicethkoc_in
    real,  intent(in) ::  varicethkld_in
    real,  intent(in) ::  variceamtoc_in
    real,  intent(in) ::  variceamtld_in
    real,  intent(in) ::  varicedeffoc_in
    real,  intent(in) ::  varicedeffld_in
    real,  intent(in) ::  varTskinOc_in
    real,  intent(in) ::  varTskinLd_in

    !------------------------------------------------------------------------
    ! Set global variables from input arguments
    !------------------------------------------------------------------------
    backcldtopoc   = backcldtopoc_in
    backcldtopld   = backcldtopld_in
    backcldthkoc   = backcldthkoc_in
    backcldthkld   = backcldthkld_in
    backcldamtoc   = backcldamtoc_in
    backcldamtld   = backcldamtld_in
    backclddeffoc  = backclddeffoc_in
    backclddeffld  = backclddeffld_in
    varcldtopoc    = varcldtopoc_in
    varcldtopld    = varcldtopld_in
    varcldthkoc    = varcldthkoc_in
    varcldthkld    = varcldthkld_in
    varcldamtoc    = varcldamtoc_in
    varcldamtld    = varcldamtld_in
    varclddeffoc   = varclddeffoc_in
    varclddeffld   = varclddeffld_in
    backicetopoc   = backicetopoc_in
    backicetopld   = backicetopld_in
    backicethkoc   = backicethkoc_in
    backicethkld   = backicethkld_in
    backiceamtoc   = backiceamtoc_in
    backiceamtld   = backiceamtld_in
    backicedeffoc  = backicedeffoc_in
    backicedeffld  = backicedeffld_in
    varicetopoc    = varicetopoc_in
    varicetopld    = varicetopld_in
    varicethkoc    = varicethkoc_in
    varicethkld    = varicethkld_in
    variceamtoc    = variceamtoc_in
    variceamtld    = variceamtld_in
    varicedeffoc   = varicedeffoc_in
    varicedeffld   = varicedeffld_in
    varTskinOc     = varTskinOc_in
    varTskinLd     = varTskinLd_in

    loadBkgTuningDone=.true.
    return
  end subroutine setPropertiesBkgTuning

  !--
  ! Read items from namelist into local variables
  !--
  subroutine loadBkgTuning(file_algcfg)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   loadBkgTuning
!
! PURPOSE:
!
!   Load background tuning parameters from tuning data file.
!
! SYNTAX:
!
!   CALL loadBkgTuning(file_algcfg)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   file_algcfg  CHAR  Tuning data file
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !--I/O variables
    character(len=*),              intent(in)  :: file_algcfg

    !---Local variables
    REAL :: backcldtopoc,backcldtopld,backcldthkoc,backcldthkld
    REAL :: backcldamtoc,backcldamtld,varcldtopoc,varcldtopld
    REAL :: varcldthkoc,varcldthkld,varcldamtoc,varcldamtld
    REAL :: backicetopoc,backicetopld,backicethkoc,backicethkld
    REAL :: backiceamtoc,backiceamtld,varicetopoc,varicetopld
    REAL :: varicethkoc,varicethkld,variceamtoc,variceamtld
    REAL :: backclddeffoc,backclddeffld,backicedeffoc,backicedeffld
    REAL :: varclddeffoc,varclddeffld,varicedeffoc,varicedeffld
    REAL :: varTskinOc=MISSING_REAL,varTskinLd=MISSING_REAL
    integer :: U_algconfig

    NAMELIST /bkgtuning/                                                  &
         backcldtopoc,backcldtopld,backcldthkoc,backcldthkld,             &
         backcldamtoc,backcldamtld,backclddeffoc,backclddeffld,           &
         varcldtopoc,varcldtopld,                                         &
         varcldthkoc,varcldthkld,varcldamtoc,varcldamtld,                 &
         varclddeffoc,varclddeffld,                                       &
         backicetopoc,backicetopld,backicethkoc,backicethkld,             &
         backiceamtoc,backiceamtld,backicedeffoc,backicedeffld,           &
         varicetopoc,varicetopld,                                         &
         varicethkoc,varicethkld,variceamtoc,variceamtld,                 &
         varicedeffoc,varicedeffld,                                       &
         varTskinOc,varTskinLd

    if(loadBkgTuningDone) return

    !------------------------------------------------------------------------
    ! Read namelist file items into local variables
    !------------------------------------------------------------------------
    !---Load retr bkgtuning parameters
    U_algconfig = getUnit()
    open(U_algconfig,file=file_algcfg)
    read(U_algconfig,bkgtuning)
    close(U_algconfig)

    !------------------------------------------------------------------------
    ! Set global variables
    !------------------------------------------------------------------------
    call setPropertiesBkgTuning(backcldtopoc,backcldtopld,backcldthkoc, &
         backcldthkld,backcldamtoc,backcldamtld,backclddeffoc,backclddeffld, &
         varcldtopoc,varcldtopld, &
         varcldthkoc,varcldthkld,varcldamtoc,varcldamtld, &
         varclddeffoc,varclddeffld, &
         backicetopoc,backicetopld,backicethkoc,backicethkld, &
         backiceamtoc,backiceamtld,backicedeffoc,backicedeffld, &
         varicetopoc,varicetopld, &
         varicethkoc,varicethkld,variceamtoc,variceamtld, &
         varicedeffoc,varicedeffld, &
         varTskinOc,varTskinLd)

    return
  end subroutine loadBkgTuning

  subroutine loadBkgClimo(file_covar1,file_covar2,molID,molTran, &
     nlevel,vCoordTyp,pref)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   loadBkgClimo
!
! PURPOSE:
!
!   Load background and covariance data from file for atmospheric and surface 
!   quantities
!
! SYNTAX:
!
!   CALL loadBkgClimo(file_covar1, file_covar2, molID, molTran)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   file_covar1  CHAR     Input background/covariance file name
!   file_covar2  CHAR     Input background/covariance file name
!   
!   INPUTS/OUTPUTS:
!   
!   molID        INTEGER  List of IDs of relevant molecular species
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    use bkg_io_module

    !--I/O variables
    character(len=*),          intent(in)    :: file_covar1,file_covar2
    integer, dimension(:),     intent(inout) :: molID
    integer, dimension(:),     intent(inout) :: molTran
    integer,                   intent(inout) :: nlevel
    CHARACTER(LEN=mxCoordTyp), intent(inout) :: vCoordTyp
    real,    dimension(:),     allocatable   :: pref
    !---Local variables
    integer                              :: itype
    integer                              :: nchmw
    integer                              :: invNotTrnsp
    real, dimension(:,:), allocatable    :: uv_atm_tmp,uv_sfc_tmp
    real                                 :: sigPtop
    real, dimension(:),   allocatable    :: vCoordA,vCoordB
    integer :: U_cov

    if(loadBkgClimoDone) return

    !------------------------------------------------------------------------
    !  Atmospheric covariance (assumed to be the 1st in covarFile)
    !------------------------------------------------------------------------
    call queryCov(file_covar1,nLevel=nlevel,IG=IG,NG=NG, &
       invNotTrnsp=invNotTrnsp,vCoordTyp=vCoordTyp)
    allocate(pref(nlevel),vCoordA(nlevel),vCoordB(nlevel),myPres(nlevel))
    nmol=getnmol(NG)
    NParGAtm = getVectorLength(NG)
    if (invNotTrnsp > 0) then
       call queryCov(file_covar1,molid=molID,IC=IC,NC=NC)
       nParAtm=getVectorLength(NC)
       trnFixedAtm=.TRUE.  ! transformation cannot be truncated
    else
       IC=IG
       NC=NG
       nParAtm=NParGAtm
       trnFixedAtm=.FALSE.  ! transformation can be truncated (EOF)
    endif
    SELECT CASE (TRIM(vCoordTyp))
    CASE (Pcoord)
       CALL queryCov(file_covar1,pressure=pref)
       myPres=pref
    CASE (Scoord)
       CALL queryCov(file_covar1,sigCoord=vCoordA,sigPtop=sigPtop)
       CALL putSigmDefin(nlevel,sigPtop,vCoordA)
    CASE (Hcoord)
       CALL queryCov(file_covar1,hybCoordA=vCoordA,hybCoordB=vCoordB)
       CALL putHybrDefin(nlevel,vCoordA,vCoordB)
    CASE default
       print*,'err[IRBkgPrepModule::loadBkgClimo]: ',&
       'Unrecognized vertical coordinate type in cov file: ', &
       TRIM(vCoordTyp)
       CALL errorHalt(1)
    END SELECT

    U_cov = getUnit()
    CALL openCov(ncid=U_cov, file=file_covar1,nbkg=nAtmTypes, molid=molID,&
         MolTran=molTran,nParG=nParGatm,IG=IG,NG=NG)
    allocate(back_clim_atm_all(NParGAtm,nAtmTypes))
    allocate(u_clim_atm_all(NParAtm,NParGAtm,nAtmTypes))
    allocate(v_clim_atm_all(NParGAtm,NParAtm,nAtmTypes))
    allocate(ut_cov_clim_atm_all(NParGAtm,NParGAtm,nAtmTypes))
    allocate(uv_atm_tmp(NParGAtm,NParGAtm))
    allocate(AtmTypes(nAtmTypes))


    DO itype = 1, nAtmTypes
       CALL getCov(ncid=U_cov,irec=itype, &
            dmean=back_clim_atm_all(1:npargatm,itype),             &
            un_atm=uv_atm_tmp(1:npargatm,1:nparAtm), &
            cov=ut_cov_clim_atm_all(1:npargatm,1:npargatm,itype),&
            TypeFlag=AtmTypes(itype))
       u_clim_atm_all(1:nparAtm,1:npargatm,itype)= &
           TRANSPOSE(uv_atm_tmp(1:npargatm,1:nparAtm))
       if (invNotTrnsp > 0) then
          CALL getCov(ncid=U_cov,irec=itype, &
               vn_atm=uv_atm_tmp(1:nparAtm,1:npargatm))
          v_clim_atm_all(1:npargatm,1:nparAtm,itype)= &
             TRANSPOSE(uv_atm_tmp(1:nparAtm,1:npargatm))
       else
          v_clim_atm_all(1:npargatm,1:nparAtm,itype)= &
             TRANSPOSE(u_clim_atm_all(1:nparAtm,1:npargatm,itype))
       endif
       if (AtmTypes(itype)<0) then
          print*,'err[IRBkgPrepModule::loadBkgClimo]: ',&
               ' AtmTypes undefined in Atm Covariance'
          call errorHalt(1)
       endif
    ENDDO
    CALL closeCov(ncid=U_cov)

    !------------------------------------------------------------------------
    !      Surface covariance (assumed to be the 2nd in covarFile)
    !------------------------------------------------------------------------
    call queryCov(file=file_covar2,nchmw=nchmw,nemir=myNemirG, &
    invNotTrnsp=invNotTrnsp)
    allocate(myFreq(nchmw),myPol(nchmw),myFrqEmIR(myNemirG))

    if (invNotTrnsp > 0) then
       call queryCov(file_covar2,nTranMW=nTranMW,nTranIR=nTranIR)
       trnFixedSfc=.TRUE.  ! transformation cannot be truncated
    else
       nTranMW=nchmw
       nTranIR=2*myNemirG
       trnFixedSfc=.FALSE.  ! transformation can be truncated (EOF)
    endif
    npargsfc=nchmw+2*myNemirG
    nparSfc=nTranMW+nTranIR
    ! The below 3 lines were inserted to support validation of size of priorResult in
    ! primmaryyRetr. Insertion caused an error.
    ! Revisit when MW and primary bkg are separate
!    NG%emMW=nchmw
!    NG%emRfIR=2*myNemirG
!    IG=genIndices(NG)

    U_cov = getUnit()
    CALL openCov(ncid=U_cov, file=file_covar2,nbkg=nSfcTypes, &
         freq=myFreq(1:nchmw),frqemir=myFrqEmIR(1:myNemirG), &
         polarity=myPol(1:nchmw))

    allocate(back_clim_sfc_all(NParGSfc,nSfcTypes))
    allocate(u_clim_sfc_all(NParSfc,NParGSfc,nSfcTypes))
    allocate(v_clim_sfc_all(NParGSfc,NParSfc,nSfcTypes))
    allocate(ut_cov_clim_sfc_all(NParGSfc,NParGSfc,nSfcTypes))
    allocate(uv_sfc_tmp(NParGSfc,NParGSfc))
    allocate(SfcTypes(nSfcTypes))
    u_clim_sfc_all(:,:,:)=0.
    v_clim_sfc_all(:,:,:)=0.
    uv_sfc_tmp(:,:)=0.
    ut_cov_clim_sfc_all(:,:,:)=0.
    DO itype = 1, nSfcTypes
       CALL getCov(ncid=U_cov,irec=itype,&
            dmEmMw=back_clim_sfc_all(1:nchmw,itype),&
            un_emmw=uv_sfc_tmp(1:nchmw,1:nTranMW),&
            covEmMw=ut_cov_clim_sfc_all(1:nchmw,1:nchmw,itype),&
            dmEmIR=back_clim_sfc_all(nchmw+1:npargsfc,itype),&
            un_emIR=uv_sfc_tmp(nchmw+1:npargsfc,nTranMW+1:nparSfc),&
            covEmIR=ut_cov_clim_sfc_all(nchmw+1:npargsfc,nchmw+1:npargsfc,itype),&
            TypeFlag=SfcTypes(itype))
       u_clim_sfc_all(1:nparSfc,1:npargsfc,itype)= &
           TRANSPOSE(uv_sfc_tmp(1:npargsfc,1:nparSfc))
       if (invNotTrnsp > 0) then
          CALL getCov(ncid=U_cov,irec=itype,&
               vn_emmw=uv_sfc_tmp(1:nTranMW,1:nchmw),&
               vn_emIR=uv_sfc_tmp(nTranMW+1:nparSfc,nchmw+1:npargsfc))
          v_clim_sfc_all(1:npargsfc,1:nparSfc,itype)= &
              TRANSPOSE(uv_sfc_tmp(1:nparSfc,1:npargsfc))
       else
          v_clim_sfc_all(1:npargsfc,1:nparSfc,itype)= &
             TRANSPOSE(u_clim_sfc_all(1:nparSfc,1:npargsfc,itype))
       endif
       if (SfcTypes(itype) < 0) then
          print*,'err[IRBkgPrepModule::loadBkgClimo]: ',&
               ' SfcTypes undefined in Sfc Covariance'
          call errorHalt(1)
       endif
    ENDDO

    CALL closeCov(ncid=U_cov)

    !------------------------------------------------------------------------
    ! Set global variables from file variables
    !------------------------------------------------------------------------
    NemmwG        = nchmw

    deallocate(uv_atm_tmp,uv_sfc_tmp)

    loadBkgClimoDone=.true.
    return
  end subroutine loadBkgClimo

  SUBROUTINE IRsetFOVbkg(iclassatm,igeo,psfc,vCoordTyp,press,umtx,vmtx)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   IRsetFOVbkg
!
! PURPOSE:
!
!   Check consistency between the selected and available classes, for atmosphere 
!   and surface.
!
! SYNTAX:
!
!   CALL IRsetFOVbkg(iclassatm, igeo, psfc, vCoordTyp, press, umtx, vmtx)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   iclassatm  INTEGER  Atmospheric class index
!   igeo       INTEGER  Surface (geography) class index
!   psfc       REAL     Surface pressure
!   
!   INPUTS/OUTPUTS:
!   
!   umtx       REAL     Transformation matrix between geophysical and 
!                       retrieval spaces
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !---Input variables
    INTEGER,                 INTENT(IN)    :: iclassatm
    INTEGER,                 INTENT(IN)    :: igeo
    REAL,                    INTENT(IN)    :: psfc
    CHARACTER(LEN=mxCoordTyp),INTENT(IN)   :: vCoordTyp
    REAL,    DIMENSION(:),   INTENT(INOUT) :: press
    REAL,    DIMENSION(:,:), INTENT(INOUT) :: umtx
    REAL,    DIMENSION(:,:), INTENT(INOUT) :: vmtx

    !---Check the consistency between the igeo & SfcTypes
    if(ALL(SfcTypes(1:nSfcTypes)/=igeo)) then
       print *,'err[IRBkgPrepModule::IRsetFOVbkg]:  Sfctype not specified'
       call errorHalt(1)
    end if
    !---Check the consistency between the iclassatm and AtmTypes
    if(ALL(AtmTypes(1:nAtmTypes)/=iclassatm)) then
       print *,'err[IRBkgPrepModule::IRsetFOVbkg]:Atmtype not specified'
       call errorHalt(1)
    end if
    !---Select climatology information first
    CALL IRgetBkgClim(iclassatm,igeo,psfc,umtx,vmtx)
    !---Get scene pressure profile
    SELECT CASE (TRIM(vCoordTyp))
    CASE (Pcoord)
       press=myPres
    CASE (Scoord)
       press=getSigmPres(psfc)
    CASE (Hcoord)
       press=getHybrPres(psfc)
    END SELECT
    RETURN
  END SUBROUTINE IRsetFOVbkg

  SUBROUTINE IRcheckTransformConform(NR,nemMWsel,nemIRsel)

    !---Input variables
    TYPE(StateIndex_t), INTENT(IN)   :: NR
    INTEGER,            INTENT(IN)   :: nemMWsel
    INTEGER,            INTENT(IN)   :: nemIRsel
    !---Local variables
    INTEGER :: j

    !--- Verify that the selected number of variables to retrieve are
    !--- compatible with the transformation matrices
    IF (trnFixedAtm) THEN
       IF ((NR%temp /= 0) .AND. (NR%temp /= NC%temp)) THEN
          print*,'err[IRBkgPrepModule::IRcheckTransformConform]: ',&
               'Cannot truncate transform of temperature profile'
          call errorHalt(1)
       ENDIF
       IF ((NR%tskin /= 0) .AND. (NR%tskin /= NC%tskin)) THEN
          print*,'err[IRBkgPrepModule::IRcheckTransformConform]: ',&
               'Cannot truncate transform of surface temperature'
          call errorHalt(1)
       ENDIF
       IF ((NR%psfc /= 0) .AND. (NR%psfc /= NC%psfc)) THEN
          print*,'err[IRBkgPrepModule::IRcheckTransformConform]: ',&
               'Cannot truncate transform of surface temperature'
          call errorHalt(1)
       ENDIF
       DO j=1,nmol
          IF ((NR%mol(j) /= 0) .AND. (NR%mol(j) /= NC%mol(j))) THEN
             print*,'err[IRBkgPrepModule::IRcheckTransformConform]: ',&
                  'Cannot truncate transform of gas profile',j
             call errorHalt(1)
          ENDIF
       ENDDO
       IF ((NR%cldliq /= 0) .AND. (NR%cldliq /= NC%cldliq)) THEN
          print*,'err[IRBkgPrepModule::IRcheckTransformConform]: ',&
               'Cannot truncate transform of cloud liquid profile'
          call errorHalt(1)
       ENDIF
       IF ((NR%cldice /= 0) .AND. (NR%cldice /= NC%cldice)) THEN
          print*,'err[IRBkgPrepModule::IRcheckTransformConform]: ',&
               'Cannot truncate transform of cloud ice profile'
          call errorHalt(1)
       ENDIF
       IF ((NR%wind /= 0) .AND. (NR%wind /= NC%wind)) THEN
          print*,'err[IRBkgPrepModule::IRcheckTransformConform]: ',&
               'Cannot truncate transform of wind'
          call errorHalt(1)
       ENDIF
    ENDIF

    IF (trnFixedSfc) THEN
       IF ((nemMWsel /= 0) .AND. (nemMWsel /= nTranMW)) THEN
          print*,'err[IRBkgPrepModule::IRcheckTransformConform]: ',&
               'Cannot truncate transform of MW emissivity'
          call errorHalt(1)
       ENDIF
       IF ((nemIRsel /= 0) .AND. (2*nemIRsel /= nTranIR)) THEN
          print*,'err[IRBkgPrepModule::IRcheckTransformConform]: ',&
               'Cannot truncate transform of IR emissivity'
          call errorHalt(1)
       ENDIF
    ENDIF

  END SUBROUTINE IRcheckTransformConform

  SUBROUTINE IRgetBkgClim(iclassatm,igeo,psfc,umtx,vmtx)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   IRgetBkgClim
!
! PURPOSE:
!
!   Get climatology information on background, covariance and transformation 
!   matrix for atmosphere and surface, for selected classes.
!
! SYNTAX:
!
!   CALL IRgetBkgClim(iclassatm, igeo, psfc, umtx, vmtx)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   iclassatm  INTEGER  Atmospheric class index
!   igeo       INTEGER  Surface (geography) class index
!   psfc       REAL     Surface pressure
!   
!   INPUTS/OUTPUTS:
!   
!   umtx       REAL     Transformation matrix between geophysical and 
!                       retrieval spaces
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !---Input variables
    INTEGER,                 INTENT(IN)    :: iclassatm,igeo
    REAL,                    INTENT(IN)    :: psfc
    REAL,    DIMENSION(:,:), INTENT(INOUT) :: umtx
    REAL,    DIMENSION(:,:), INTENT(INOUT) :: vmtx
    !---Local variables
    INTEGER :: i,j
    
    !----Background selection
    xbakg_clim(1:npargatm)        = back_clim_atm_all(1:npargatm,iclassatm)
    xbakg_clim(npargatm+1:nparg)  = back_clim_sfc_all(1:npargsfc,igeo)

    !----Covariance matrix selection
    ut_cov_clim=0.
    !--atmosphere
    ut_cov_clim(1:nparGAtm,1:nparGAtm)= &
         ut_cov_clim_atm_all(1:nparGAtm,1:nparGATm,iclassatm)
    !--surface
    ut_cov_clim(nparGATm+1:nparG,nparGATm+1:nparG)= &
         ut_cov_clim_sfc_all(1:npargsfc,1:npargsfc,igeo)
    
    !----Transformation matrix selection
    umtx_clim=0.
    umtx_clim(1:NParAtm,1:NParGatm)                = &
         u_clim_atm_all(1:nparAtm,1:npargatm,iclassatm)
    umtx_clim(NParAtm+1:nparAtm+nparSfc,NParGatm+1:NparG)= &
         u_clim_sfc_all(1:nparSfc,1:npargsfc,igeo)
    vmtx_clim=0.
    vmtx_clim(1:NParGatm,1:NParAtm)                = &
         v_clim_atm_all(1:npargatm,1:nparAtm,iclassatm)
    vmtx_clim(NParGatm+1:NparG,NParAtm+1:nparAtm+nparSfc)= &
         v_clim_sfc_all(1:npargsfc,1:nparSfc,igeo)
    !----Truncate (compress) the transformation matrix
    !--This initialization is necessary on Linux
    !--because otherwise the implicit do loop will 
    !--begin at j=0 and cause a segmentation fault
    !--This needs to be further investigated!
    j=1

    umtx(1:NPar,:)=umtx_clim(extrvec(1:NPar),:)
    vmtx(:,1:NPar)=vmtx_clim(:,extrvec(1:NPar))
    RETURN
  END SUBROUTINE IRgetBkgClim


 !-----------------------------------------------------------------
 ! The following subroutine is intended to be used to set 
 ! the o/p parameters such as IR_out, NR_out, nPar_out 
 ! using the values in corresponding input parameters and
 ! conditional on the retrievalFlags.
 ! Dynamically selected flags from pre-defined sets are used
 ! to determine which variable is to be set. 
 !----------------------------------------------------------------- 
  SUBROUTINE IRsetRetrVec(retrievalFlags, NR_in, nEmMw_in, nemIR_in, &
                                           nPar_out, NR_out, IR_out)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   IRsetRetrVec
!
! PURPOSE:
!
!   Set the sizes and on/off status of components of retrieval vector using the 
!   values in corresponding input parameters and conditional on the 
!   retrievalFlags.
!
! SYNTAX:
!
!   CALL IRsetRetrVec(retrievalFlags, NR_in, nEmMw_in, nemIR_in, 
!      nPar_out, NR_out, IR_out)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   retrievalFlags  TYPE(RETFLAGS_T)    Flags telling which parameters to be 
!                                       retrieved 
!   NR_in           TYPE(STATEINDEX_T)  Number of retrieved elements for 
!                                       sections of retrieval state vector, 
!                                       for input 
!   nEmMw_in        INTEGER             Number of retrieved variables for MW 
!                                       emissivity 
!   nemIR_in        INTEGER             Number of hinge points for retrieving 
!                                       IR emissivity 
!   
!   INPUTS/OUTPUTS:
!   
!   nPar_out        INTEGER             Total number of elements in retrieval 
!                                       vector 
!   NR_out          TYPE(STATEINDEX_T)  Number of retrieved elements for 
!                                       sections of retrieval state vector, 
!                                       for output 
!   IR_out          TYPE(STATEINDEX_T)  Starting indexes for the groups in 
!                                       retrieval vector
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>


    TYPE(retFlags_t),     INTENT(IN)   :: retrievalFlags
    TYPE(StateIndex_t),   INTENT(IN)   :: NR_in
    INTEGER,              INTENT(IN)   :: nEmMw_in
    INTEGER,              INTENT(IN)   :: nemIR_in
    INTEGER,              INTENT(INOUT):: nPar_out
    TYPE(StateIndex_t),   INTENT(INOUT):: NR_out
    TYPE(StateIndex_t),   INTENT(INOUT):: IR_out

    ! some local variables
    INTEGER          :: i, j
    INTEGER          :: IEmMwC
    INTEGER          :: IEmIrC
    TYPE(StateIndex_t) :: IR,NR

    nemmw         = nEmMw_in
    nemir         = nemIR_in

    NR  = NR_in

    ! While the parent subroutines do not set the surface parts of NR_in:
    NR%emMW = nemmw
    NR%emRfIR = 2 * nemir

    IF (.NOT. retrievalFlags%tempFlag) THEN
        NR%temp = 0
    END IF

    WHERE (retrievalFlags%molFlags)
        NR%mol = NR_in%mol
    ELSEWHERE
        NR%mol = 0
    END WHERE

    IF (.NOT. retrievalFlags%cldLiqFlag) THEN
       NR%cldLiq = 0
    END IF
    
    IF (.NOT. retrievalFlags%cldIceFlag) THEN
       NR%cldIce = 0
    END IF
    
    IF (.NOT. retrievalFlags%psfcFlag ) THEN
       NR%psfc = 0
    END IF
    
    IF (.NOT. retrievalFlags%tskinFlag ) THEN
       NR%tskin = 0 
    END IF
    
    IF (.NOT. retrievalFlags%windFlag ) THEN
       NR%wind = 0
    END IF
    
    IF (.NOT. retrievalFlags%emissMwFlag ) THEN
       NR%emMW = 0
    END IF
    
    IF (.NOT. retrievalFlags%emissIrFlag ) THEN
       NR%emRfIR = 0
    END IF
    
    NPar = getVectorLength(NR)
    npar_out = NPar
    IR = genIndices(NR)
    iemmw         = IR%wind + NR%wind

    NR_out = NR
    IR_out = IR

    IEmMwC = IC%wind+NC%wind
    IEmIrC = IEmMwC+nTranMW

    j = 1
    extrvec(1:npar)= (/(i,i=IC%Temp,IC%Temp+NR%Temp-1), &
        (i,i=IC%Tskin,IC%Tskin+NR%Tskin-1), (i, i=IC%Psfc,IC%Psfc+NR%Psfc-1), &
        ((i,i=IC%mol(j),IC%mol(j)+NR%mol(j)-1),j=1,nmol),   &
        (i,i=IC%Cldliq,IC%Cldliq+NR%Cldliq-1),              &
        (i,i=IC%CldIce,IC%CldIce+NR%CldIce-1),              &
        (i,i=IEmMwC,IEmMwC+NR%emMW-1),                 &
        (i,i=IEmIrC,IEmIrC+NR%emRfIR-1)/)
    RETURN
    END SUBROUTINE IRsetRetrVec

  SUBROUTINE IRsetExtBkg(xbakg_ext)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   IRsetExtBkg
!
! PURPOSE:
!
!   Load external background data into global vector
!
! SYNTAX:
!
!   CALL IRsetExtBkg(xbakg_ext)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   xbakg_ext  REAL  External background state vector
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    REAL,DIMENSION(:),INTENT(in) :: xbakg_ext
    !---Load external background data into global vector 'xExt' which 
    !    will be accessed in IRcombBkg.  When accessing external data,
    !    this subroutine should be called immediately before calling 
    !    IRcombBkg.
    xExt = xbakg_ext
    extAtmosOn = .true.
    RETURN
  END SUBROUTINE IRsetExtBkg

  SUBROUTINE IRsetEmisBkg(xbakg_emis,nemirExt)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   IRsetEmisBkg
!
! PURPOSE:
!
!   Load external emissivity and reflectivity data into global vector.
!
! SYNTAX:
!
!   CALL IRsetEmisBkg(xbakg_emis, nemirExt)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   xbakg_emis  REAL     External background IR emissivity/reflectivity 
!                        vector 
!   nemirExt    INTEGER  Number of IR emissivity hinge points for external 
!                        data
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !---Load external emissivity and reflectivity data into global vector 
    !   'emRfExt' that will be accessed in IRcombBkg.
    !    When accessing external data, this subroutine should be 
    !    called immediately before calling IRcombBkg.

    REAL,    DIMENSION(:), INTENT(in) :: xbakg_emis
    INTEGER,               INTENT(in) :: nemirExt

    ! Constants
    real,    parameter             :: emIRcap=0.995
    real,    parameter             :: emIRfloor=0.6

    ! Local variables
    INTEGER                        :: i

    IF (ALL(xbakg_emis == MISSING_REAL)) THEN
      emRfExt=MISSING_REAL
      RETURN
    ENDIF

    ! Interpolate spectrally, as needed
    IF (nemirExt /= myNemIrG) THEN
      emRfExt(1:myNemIrG)=matmul(cvtIR,xbakg_emis(1:myNemIrG))
      emRfExt(myNemIrG+1:myNemIrG+myNemIrG)=&
           matmul(cvtIR,xbakg_emis(myNemIrG+1:myNemIrG+myNemIrG))
    ELSE
      emRfExt(1:myNemIrG)=xbakg_emis(1:myNemIrG)
      emRfExt(myNemIrG+1:myNemIrG+myNemIrG)=&
           xbakg_emis(myNemIrG+1:myNemIrG+myNemIrG)
    ENDIF
    DO i=1,myNemIrG
      emRfExt(i)=min(emRfExt(i),emIRcap) 
      emRfExt(i)=max(emRfExt(i),emIRfloor) 
    ENDDO

    DO i=myNemIrG+1,myNemIrG+myNemIrG
      emRfExt(i)=min(emRfExt(i),1.0-emIRcap) 
      emRfExt(i)=max(emRfExt(i),1.0-emIRfloor) 
    ENDDO

    extEmRfon = .true.
    RETURN
  END SUBROUTINE IRsetEmisBkg

  SUBROUTINE IRcombBkg(IGExt,NGExt,psfc,isLand,molTran, &
     vCoordTyp,press,dback0,ut_cov_net)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   IRcombBkg
!
! PURPOSE:
!
!   Load climatology background into the worskpace array. Background may be 
!   overwritten with external data.
!
! SYNTAX:
!
!   CALL IRcombBkg(IGExt, NGExt, psfc, pland, plandMaxOc, molTran, dback0, 
!      ut_cov_net)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   IGExt       TYPE(STATEINDEX_T)  Starting indices for sections of 
!                                   geophysical state vector for external 
!                                   data 
!   NGExt       TYPE(STATEINDEX_T)  Number of elements for sections of 
!                                   geophysical state vector for external 
!                                   data 
!   psfc        REAL                Surface pressure
!   isLand      logical             flag for land or not
!   molTran     INTEGER             Molecule transformation index
!   
!   INPUTS/OUTPUTS:
!   
!   dback0      REAL                background state vector for the current 
!                                   profile 
!   ut_cov_net  REAL                background error covariance for the 
!                                   current profile
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    USE LvlInterp
    USE MapInvert, ONLY: MOL_TRAN_LOG
    implicit none
    TYPE(StateIndex_t),   INTENT(IN)    :: IGExt,NGExt
    REAL,                 INTENT(IN)    :: psfc
    LOGICAL,              INTENT(IN)    :: isLand
    INTEGER, DIMENSION(:),INTENT(IN)    :: molTran
    CHARACTER(LEN=mxCoordTyp),INTENT(IN):: vCoordTyp
    REAL, DIMENSION(:),   INTENT(IN)    :: press
    REAL, DIMENSION(:),   INTENT(INOUT) :: dback0
    REAL, DIMENSION(:,:), INTENT(INOUT) :: ut_cov_net
    !---Local variables
    INTEGER :: i,ii,nsurf
    REAL    :: xSfc,dxlvldp,deltExtrap
    character (len=256) :: procName
    logical :: dbg=.false.
    logical :: noExt
    
    procName = '[IRBkgPrepModule::IRBcombBkg]:'
    !---By default load the climatology background into the worskpace array
    dback0(1:NParG)=xbakg_clim(1:NParG)
    ut_cov_net(1:nparG,1:nparG)    = ut_cov_clim(1:nparG,1:nparG)
    dback0(IG%Psfc)=psfc

    noExt = .TRUE.
    IF (getVectorLength(NGExt) > 0) THEN
       IF ( ANY(xExt /= MISSING_REAL) ) noExt=.FALSE.
    ENDIF
    IF ( noExt .AND. &
         ANY( (/ &
         externalFlags%xTemp  , externalFlags%xTskin  , &
         externalFlags%xPsfc  , externalFlags%xMols(:), &
         externalFlags%xCldLiq, externalFlags%xCldIce , &
         externalFlags%xWind /) ) ) THEN
       print*,'err[IRBkgPrepModule::IRcombBkg]:  ', &
            'External data requested, but external data vector does not ',&
            'contain any valid data.',&
            'Call IRsetExtBkg prior to calling IRcombBkg.'
       call errorHalt(1)
    ENDIF

    IF (externalFlags%xPsfc) then
       dback0(IG%Psfc)=xExt(IGExt%psfc)
    ENDIF

    IF (TRIM(vCoordTyp) == Pcoord) THEN
       !---background may be overwritten with external data
       !---LVL_INT anywhere below is called to extrapolate the values below surface
       !---and get rid of '-9999' in modification of BKG
       IF (externalFlags%xTemp) then
          call lvl_int(xExt(IGExt%temp:IGExt%temp+NGExt%temp-1),press, &
             NGExt%temp,xExt(IGExt%Psfc),xSfc,dxlvldp,Nsurf)
          deltExtrap=xExt(IGExt%temp+Nsurf-1)-dback0(IG%temp+Nsurf-1)
          dback0(IG%temp:IG%temp+Nsurf-1)=&
             xExt(IGExt%temp:IGExt%temp+Nsurf-1)
!    For all levels below surface (levels >= Nsurf) we set a constant deviation
!    from background for an external Temp profile = xExt(Nsurf)-dback0(Nsurf) 
          dback0(IG%temp+Nsurf:IG%temp+NG%temp-1)=&
             dback0(IG%temp+Nsurf:IG%temp+NG%temp-1)+deltExtrap
       ENDIF
       IF (externalFlags%xTskin) then
          dback0(IG%tskin)=xExt(IGExt%tskin)-xSfc ! Delta-Tskin with fixed grid
       ENDIF
       DO i=1,nmolExt
          ii=mapExtMols(i)
          IF (externalFlags%xMols(i)) then
             call lvl_int(xExt(IGExt%mol(i):IGExt%mol(i)+NGExt%mol(i)-1),press, &
                  NGExt%mol(i),xExt(IGExt%Psfc),xSfc,dxlvldp,Nsurf)
!    The same way as for Temp profile, we do for gases profiles.
             if (molTran(ii) == MOL_TRAN_LOG) then
                deltExtrap=alog(xExt(IGExt%mol(i)+Nsurf-1))-dback0(IG%mol(ii)+Nsurf-1)
                dback0(IG%mol(ii):IG%mol(ii)+Nsurf-1)=&
                     alog(xExt(IGExt%mol(i):IGExt%mol(i)+Nsurf-1))
             else
                deltExtrap=xExt(IGExt%mol(i)+Nsurf-1)-dback0(IG%mol(ii)+Nsurf-1)
                dback0(IG%mol(ii):IG%mol(ii)+Nsurf-1)=&
                     xExt(IGExt%mol(i):IGExt%mol(i)+Nsurf-1)
             endif
             dback0(IG%mol(ii)+Nsurf:IG%mol(ii)+NG%mol(ii)-1)=&
                  dback0(IG%mol(ii)+Nsurf:IG%mol(ii)+NG%mol(ii)-1)+&
                  deltExtrap
          ENDIF
       ENDDO
    ELSE
       IF (externalFlags%xTemp) &
          dback0(IG%temp:IG%temp+NG%temp-1)=&
             xExt(IGExt%temp:IGExt%temp+NGExt%temp-1)
       IF (externalFlags%xTskin) &
          dback0(IG%tskin)=xExt(IGExt%tskin)
       DO i=1,nmolExt
          ii=mapExtMols(i)
          IF (externalFlags%xMols(i)) then
             if (molTran(ii) == MOL_TRAN_LOG) then
                dback0(IG%mol(ii):IG%mol(ii)+NG%mol(ii)-1)= &
                     alog(xExt(IGExt%mol(i):IGExt%mol(i)+NGExt%mol(i)-1))
             else
                dback0(IG%mol(ii):IG%mol(ii)+NG%mol(ii)-1)= &
                          xExt(IGExt%mol(i):IGExt%mol(i)+NGExt%mol(i)-1)
             endif
          ENDIF
       ENDDO
    ENDIF
    IF (externalFlags%xCldLiq) then
       dback0(IG%cldliq:IG%cldliq+NG%cldliq-1)=&
            xExt(IGExt%cldliq:IGExt%cldliq+NGExt%cldliq-1)
    ENDIF
    IF (externalFlags%xCldIce) then
       dback0(IG%cldice:IG%cldice+NG%cldice-1)=&
            xExt(IGExt%cldice:IGExt%cldice+NGExt%cldice-1)
    ENDIF
    IF (externalFlags%xWind) then
       dback0(IG%wind:IG%wind+NG%wind-1)=&
            xExt(IGExt%wind:IGExt%wind+NGExt%wind-1)
    ENDIF
    if (dbg) then
       print *, trim(procName)//' nParGatm,NemmwG,,NParG= ', nParGatm,NemmwG,NParG
       print *, trim(procName)//' size(dback0) = ', size(dback0)
       print *, trim(procName)//' size(emRfExt) = ', size(emRfExt)
    end if
    IF (isLand .and. extAtmosOn .and. extEmRfon) THEN
       IF (.NOT. ANY(emRfExt == MISSING_REAL)) THEN
          dback0(nParGatm+NemmwG+1:NParG)=emRfExt
       ENDIF
    ELSE
       ! placeholder for use of model that depends on conditions (wind, etc.)
    ENDIF

    RETURN
  END SUBROUTINE IRcombBkg


   SUBROUTINE IRtuneBkg(dback0,ut_cov_net,iclassatm,igeo)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   IRtuneBkg
!
! PURPOSE:
!
!   Tune background (a priori and error covariance) data, overriding any prior 
!   settings.
!
! SYNTAX:
!
!   CALL IRtuneBkg(dback0, ut_cov_net, iclassatm, igeo)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   iclassatm   INTEGER  Atmospheric class index
!   igeo        INTEGER  Surface (geography) class index
!   
!   INPUTS/OUTPUTS:
!   
!   dback0      REAL     background state vector for the current profile
!   ut_cov_net  REAL     background error covariance for the current profile
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !---Input variables
    REAL, DIMENSION(:),   INTENT(INOUT) :: dback0
    REAL, DIMENSION(:,:), INTENT(INOUT) :: ut_cov_net
    INTEGER,              INTENT(IN)    :: iclassatm,igeo

    IF(igeo.eq.1.or.igeo.eq.3.or.igeo.eq.4)THEN     !---Ocean
       dback0(IG%cldliq)    = backcldtopoc
       dback0(IG%cldliq+1)  = backcldthkoc
       dback0(IG%cldliq+2)  = backcldamtoc
       if (NG%CldLiq .GT. 3) dback0(IG%cldliq+3) =backclddeffoc
       !Liq
       ut_cov_net(IG%cldliq,IG%cldliq)     = varcldtopoc
       ut_cov_net(IG%cldliq+1,IG%cldliq+1) = varcldthkoc
       ut_cov_net(IG%cldliq+2,IG%cldliq+2) = varcldamtoc
       if (NG%CldLiq .GT. 3) ut_cov_net(IG%cldliq+3,IG%cldliq+3) =varclddeffoc
       !Ice
       IF (NG%cldIce .GT. 0) THEN
          dback0(IG%cldIce)   = backicetopoc
          dback0(IG%cldIce+1) = backicethkoc
          dback0(IG%cldIce+2) = backiceamtoc
          dback0(IG%cldice+3) = backicedeffoc
          ut_cov_net(IG%cldice,IG%cldice)     = varicetopoc
          ut_cov_net(IG%cldice+1,IG%cldice+1) = varicethkoc
          ut_cov_net(IG%cldice+2,IG%cldice+2) = variceamtoc
          ut_cov_net(IG%cldIce+3,IG%cldIce+3) = varicedeffoc
       ENDIF
       IF (NG%Tskin .GT. 0 .AND. varTskinOc .NE. MISSING_REAL) &
          ut_cov_net(IG%Tskin,IG%Tskin) = varTskinOc
    ELSE                                   !---Land
       dback0(IG%cldliq)   = backcldtopld
       dback0(IG%cldliq+1) = backcldthkld
       dback0(IG%cldliq+2) = backcldamtld
       if (NG%CldLiq .GT. 3) dback0(IG%cldliq+3) = backclddeffld
       ut_cov_net(IG%cldliq,IG%cldliq)     = varcldtopld
       ut_cov_net(IG%cldliq+1,IG%cldliq+1) = varcldthkld
       ut_cov_net(IG%cldliq+2,IG%cldliq+2) = varcldamtld
       if (NG%CldLiq .GT. 3) ut_cov_net(IG%cldliq+3,IG%cldliq+3) = varclddeffld
       IF (NG%cldIce .GT. 0) THEN
          dback0(IG%cldice)   = backicetopld
          dback0(IG%cldice+1) = backicethkld
          dback0(IG%cldice+2) = backiceamtld
          dback0(IG%cldice+3) = backicedeffld
          ut_cov_net(IG%cldice,IG%cldice)     = varicetopld
          ut_cov_net(IG%cldice+1,IG%cldice+1) = varicethkld
          ut_cov_net(IG%cldice+2,IG%cldice+2) = variceamtld
          ut_cov_net(IG%cldIce+3,IG%cldIce+3) =varicedeffld
       ENDIF
       IF (NG%Tskin .GT. 0 .AND. varTskinLd .NE. MISSING_REAL) &
          ut_cov_net(IG%Tskin,IG%Tskin) = varTskinLd
    ENDIF

    RETURN
  END SUBROUTINE IRtuneBkg

  SUBROUTINE IRtransBkg(ut_cov_net,ut_cov,umtx)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   IRtransBkg
!
! PURPOSE:
!
!   Transform covariance from physical space to retrieval space.
!
! SYNTAX:
!
!   CALL IRtransBkg(ut_cov_net, ut_cov, umtx)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   ut_cov_net  REAL  background error covariance for the current profile
!   
!   INPUTS/OUTPUTS:
!   
!   ut_cov      REAL  background error covariance for the current profile, in 
!                     retrieval space 
!   umtx        REAL  Transformation matrix between geophysical and retrieval 
!                     spaces
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !---Input variables
    REAL, DIMENSION(:,:), INTENT(IN)    :: ut_cov_net
    REAL, DIMENSION(:,:), INTENT(INOUT) :: ut_cov,umtx

    tmpl(1:nPar, 1:nParG) = MATMUL(umtx(1:NPar,1:NParG), &
         ut_cov_net(1:NParG,1:NParG))
    ut_cov(1:NPar,1:NPar)=MATMUL(tmpl(1:nPar, 1:nParG), &
         TRANSPOSE(umtx(1:NPar,1:NParG)))

    RETURN
  END SUBROUTINE IRtransBkg

  SUBROUTINE validateExt(NGExt,nmol,molID,molIDext, &
                         xFlgTemp,xFlgMols,xFlgPsfc,xFlgWind, &
                         xFlgCldLiq,xFlgCldIce)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   validateExt
!
! PURPOSE:
!
!   Check external data structure against background data structure and 
!   retrieval switches.
!
! SYNTAX:
!
!   CALL validateExt(NGExt, nmol, molID, molIDext, xFlgTemp, 
!      xFlgMols, xFlgPsfc, xFlgWind, xFlgCldLiq, xFlgCldIce)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   NGExt       TYPE(STATEINDEX_T)  Number of elements for sections of 
!                                   geophysical state vector for external 
!                                   data 
!   nmol        INTEGER             Number of molecular species
!   molID       INTEGER             List of IDs of relevant molecular species
!   molIDext    INTEGER             List of IDs of relevant molecular species 
!                                   for external data 
!   xFlgTemp    LOGICAL             Flag for using external atmospheric 
!                                   temperature  
!   xFlgMols    LOGICAL             Flag for using external set of molecular 
!                                   concentrations 
!   xFlgPsfc    LOGICAL             Flag for using external surface pressure 
!   xFlgWind    LOGICAL             Flag for using external wind data 
!   xFlgCldLiq  LOGICAL             Flag for using external data on liquid 
!                                   clouds
!   xFlgCldIce  LOGICAL             Flag for using external data on ice  
!                                   clouds
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    TYPE(StateIndex_t),   INTENT(IN)    :: NGExt
    INTEGER,              INTENT(IN)    :: nmol
    INTEGER, DIMENSION(:),INTENT(IN)    :: molID
    INTEGER, DIMENSION(:),INTENT(IN)    :: molIDext
    LOGICAL,              INTENT(IN)    :: xFlgTemp
    LOGICAL, DIMENSION(:),INTENT(IN)    :: xFlgMols
    LOGICAL,              INTENT(IN)    :: xFlgPsfc,xFlgWind
    LOGICAL,              INTENT(IN)    :: xFlgCldLiq,xFlgCldIce
    !Local variables
    INTEGER                             :: i,j
    INTEGER                             :: lenNGExt

    externalFlags = EXT_FLAGS_ALL_OFF
    lenNGExt = getVectorLength(NGExt)
    if (lenNGExt == 0) return

    !---Allocate vector to contain external background data:
    allocate(xExt(lenNGExt))
    !---Initialize to tell downstream there are no valid data
    xExt = MISSING_REAL

    if (xFlgTemp) then
       if (NGExt%temp .NE. NG%temp .OR. NGExt%tskin .NE. NG%tskin) then
          print*,'err[IRBkgPrepModule::validateExt]:  ', &
              'TEMP/or/TSKIN external fields not matching BKG data'
          call errorHalt(1)
       endif
       if (NG%Temp .GT. 0) externalFlags%xTemp=xFlgTemp
       if (NG%Tskin .GT. 0) externalFlags%xTskin=xFlgTemp
    endif
    if (xFlgPsfc) then
       if (NGExt%psfc .NE. NG%psfc) then
       print*,'err[IRBkgPrepModule::validateExt]:  ', &
           'PSFC external fields not matching BKG data'
          call errorHalt(1)
       endif
       if (NG%psfc .GT. 0) externalFlags%xPsfc=xFlgPsfc
    endif

    nmolExt=getNMol(molIDext)
    mapExtMols=0
    do j=1,nmolExt
       do i=1,nmol
          if (molIDext(j) == molID(i) .AND. xFlgMols(i)) then
             mapExtMols(j)=i
             if (NGExt%mol(j) .NE. NG%mol(i)) then
                print*,'MolID #',molID(i),':'
                print*,'err[IRBkgPrepModule::validateExt]:  ', &
                    'External fields not matching BKG data'
                call errorHalt(1)
             endif
             exit
          endif
       enddo
       if (NGExt%mol(j) .GT. 0) externalFlags%xMols(j)=xFlgMols(i)
    enddo
    if (xFlgCldLiq) then
       if (NGExt%cldLiq .NE. NG%cldLiq) then
          print*,'err[IRBkgPrepModule::validateExt]:  ', &
              'CloudLiq external fields not matching BKG data'
          call errorHalt(1)
       endif
       if (NG%cldLiq .GT. 0) externalFlags%xCldLiq=xFlgCldLiq
       externalFlags%xCldLiq=xFlgCldLiq
    endif
    if (xFlgCldIce) then
       if (NGExt%cldIce .NE. NG%cldIce) then
          print*,'err[IRBkgPrepModule::validateExt]:  ', &
              'CloudIce external fields not matching BKG data'
          call errorHalt(1)
       endif
       if (NG%cldIce .GT. 0) externalFlags%xCldIce=xFlgCldIce
    endif
    if (xFlgWind) then
       if (NGExt%wind .NE. NG%wind) then
          print*,'err[IRBkgPrepModule::validateExt]:  ', &
              'WIND external fields not matching BKG data'
          call errorHalt(1)
       endif
       if (NG%wind .GT. 0) externalFlags%xWind=xFlgWind
    endif
    RETURN
  END SUBROUTINE validateExt


  SUBROUTINE validateEmis(nemirExt,wvnExt)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   validateEmis
!
! PURPOSE:
!
!   Check external emissivity against background parameters.
!
! SYNTAX:
!
!   CALL validateEmis(nemirExt, wvnExt)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   nemirExt  INTEGER  Number of IR emissivity hinge points for external data
!   wvnExt    REAL     Wavenumbers at external IR emis hinge points
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    ! Check external emissivity against background parameters

    USE SpectralOperations, ONLY: spectralConvert

    INTEGER,               INTENT(IN)    :: nemirExt
    REAL,    DIMENSION(:), INTENT(IN)    :: wvnExt

    ! Local variables
    INTEGER                        :: i

    ! Constants
    real,    parameter             :: dSpectLim=0.001 ! Spectral match threshold
                                                      ! (fractional) 

    !---Allocate vector to contain external emissivity and reflectivity data:
    allocate(emRfExt(2*myNemIrG))
    !---Initialize to tell downstream there are no valid data
    emRfExt = MISSING_REAL

    if (nemirExt == 0) return  ! No external emis data to process

    if (nemirExt /= myNemIrG) then
      if (allocated(cvtIR)) deallocate(cvtIR)
      allocate(cvtIR(myNemIrG,nemirExt))
      call spectralConvert(wvnExt,myFrqEmIR,cvtIR,extrapOK=.TRUE.)
    elseif (any((abs(wvnExt-myFrqEmIR(1:myNemIrG))/wvnExt) > dSpectLim)) then
      print *,'err[IRBkgPrepModule::validateEmis]:  Channel set mismatch;'
      do i=1,myNemIrG
        write(*,'(i4,2f10.4)')i,wvnExt(i),myFrqEmIR(i)
      enddo
      call errorHalt(1)    
    endif

    RETURN
  END SUBROUTINE validateEmis

  subroutine cloudModeTest(cloudType,retrievalFlags,flagDynCldMask)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   cloudModeTest
!
! PURPOSE:
!
!   Set controls to disable retrieval of cloud parameters, depending 
!   on dynamic cloud mask
!
! SYNTAX:
!
!   CALL cloudModeTest(cloudType, retrievalFlags)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   cloudType       INTEGER           ID of cloud type class
!   
!   INPUTS/OUTPUTS:
!   
!   retrievalFlags  TYPE(RETFLAGS_T)  Flags telling which 
!                                     parameters to be retrieved
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
! NOTES: This subroutine was relocated from the main routine,
!        retr_sig, to be shared by both ESDR and SIGMA versions.
!        (yhe@aer.com, 06/22/2016)
!
!*******************************************************</f90Subroutine>

  ! Set controls to disable retrieval of cloud parameters,
  ! depending on dynamic cloud mask

    integer,                   intent(in)     :: cloudType
    LOGICAL,                   intent(in)     :: flagDynCldMask
    TYPE(retFlags_t),          intent(inout)  :: retrievalFlags

    select case (cloudType)
    case (0)
       retrievalFlags%cldLiqFlag=.FALSE.
       retrievalFlags%cldIceFlag=.FALSE.
    case (1)
       retrievalFlags%cldLiqFlag=.TRUE.
       retrievalFlags%cldIceFlag=.FALSE.
    case (2)
       retrievalFlags%cldLiqFlag=.FALSE.
       retrievalFlags%cldIceFlag=.TRUE.
    case (3)
       retrievalFlags%cldLiqFlag=.TRUE.
       retrievalFlags%cldIceFlag=.TRUE.
    case default
       print*,'Err[retr::cloudModeTest]: ', &
          'Invalid flagDynCldMask:',flagDynCldMask
       call errorHalt(1)
    end select

    return

  end subroutine cloudModeTest

END MODULE IRBkgPrepModule
