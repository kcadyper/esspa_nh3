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
!   Copyright AER, Inc 2001-2009, All Right Reserved
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
MODULE BackgroundModule

  ! <f2003Class>***********************************************************
  !
  ! NAME:
  !
  !   BackgroundModule
  !
  ! PURPOSE:
  !
  !   class for processing IR and MW background data (a priori and error
  !   covariance) for retrieval
  !
  ! INCLUDES:
  !
  !   None
  !
  !***********************************************************</f90Module>

  USE StateIndexModule
  USE constants, ONLY: &
       MISSING_REAL
  USE VertCoord, ONLY: &
       mxCoordTyp, &
       Pcoord, &
       Scoord, &
       Hcoord, &
       putSigmDefin, &
       putHybrDefin, &
       getSigmPres, &
       getHybrPres
  USE ToolboxModule, Only: &
       getUnit

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: RET_FLAGS_ALL_ON

  integer, parameter :: mxLength = 256
  TYPE, PUBLIC ::retFlags_t
     LOGICAL                     :: tempFlag, tskinFlag, cldLiqFlag, cldIceFlag
     LOGICAL                     :: emissMwFlag, emissIrFlag, psfcFlag
     LOGICAL                     :: iceCldThkFlag, liqCldThkFlag, windFlag    
     LOGICAL, Dimension(1:maxMol):: molFlags
  END TYPE retFlags_t

  TYPE, PUBLIC :: extFlags_t
     LOGICAL                     :: xTemp, xTskin
     LOGICAL                     :: xCldLiq, xCldIce
     LOGICAL                     :: xPsfc, xWind    
     LOGICAL, Dimension(1:maxMol):: xMols
  END TYPE extFlags_t

  TYPE(retFlags_t), PARAMETER  :: &
       RET_FLAGS_ALL_ON = retFlags_t(.TRUE.,   &
       .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., &
       .TRUE., .TRUE., .TRUE., .TRUE., .TRUE.)

  TYPE(extFlags_t), PARAMETER  :: &
       EXT_FLAGS_ALL_OFF = extFlags_t(.FALSE.,      &
       .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &
       .FALSE.)


  TYPE, PUBLIC :: background_t
     !instance variables
     !---------------------------------------------------------------
     !     Declarations of private data (private to the module)
     !---------------------------------------------------------------
     INTEGER, DIMENSION(:), allocatable      :: Pres
     integer                                 :: NemIrG
     INTEGER, DIMENSION(:), allocatable      :: AtmTypes
     INTEGER, DIMENSION(:), allocatable      :: SfcTypesMW
     INTEGER, DIMENSION(:), allocatable      :: SfcTypesIR
     INTEGER, DIMENSION(:), allocatable      :: extrvec
     INTEGER                                 :: nAtmTypes
     INTEGER                                 :: nSfcTypesMW,nSfcTypesIR
     INTEGER                                 :: nParG,nParGatm,nParGsfcMW,nParGsfcIR
     INTEGER                                 :: nPar,nParAtm,nParSfcMW,nParSfcIR
     INTEGER                                 :: nParGsum,nParSum
     INTEGER                                 :: sGsfcMW,eGsfcMW,sGsfcIR,eGsfcIR
     INTEGER                                 :: sSfcMW,eSfcMW,sSfcIR,eSfcIR
     INTEGER                                 :: nParR,nParRatm
     INTEGER                                 :: sRsfcMW,eRsfcMW,sRsfcIR,eRsfcIR
     INTEGER                                 :: iemmw,nemmwG,iemmwG
     INTEGER                                 :: nmol,ih2o
     INTEGER                                 :: nmolExt
     INTEGER, DIMENSION(maxMol)              :: mapExtMols
     TYPE(extFlags_t)                        :: externalFlags
     LOGICAL                                 :: trnFixedAtm
     LOGICAL                                 :: trnFixedSfcMW,trnFixedSfcIR
     !--- my... are inserted as intermediate variables to avoid confusing
     !--- the pgf90 optimizer and getting memory violations
     real,    dimension(:), allocatable      :: Freq,FrqEmIR
     integer, dimension(:), allocatable      :: Pol
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
     REAL, DIMENSION(:,:), allocatable       :: back_clim_atm_all
     REAL, DIMENSION(:,:), allocatable       :: back_clim_sfc_mw
     REAL, DIMENSION(:,:), allocatable       :: back_clim_sfc_ir
     !---Climo Transformation Matrices (atm and sfc)
     REAL, DIMENSION(:,:,:), allocatable     :: u_clim_atm_all
     REAL, DIMENSION(:,:,:), allocatable     :: v_clim_atm_all
     REAL, DIMENSION(:,:,:), allocatable     :: u_clim_sfc_mw
     REAL, DIMENSION(:,:,:), allocatable     :: v_clim_sfc_mw
     REAL, DIMENSION(:,:,:), allocatable     :: u_clim_sfc_ir
     REAL, DIMENSION(:,:,:), allocatable     :: v_clim_sfc_ir
     !---Climo Covariance Matrix (atm and sfc)
     REAL, DIMENSION(:,:,:), allocatable     :: ut_cov_clim_atm_all
     REAL, DIMENSION(:,:,:), allocatable     :: ut_cov_clim_sfc_mw
     REAL, DIMENSION(:,:,:), allocatable     :: ut_cov_clim_sfc_ir
     !---Total Background & Covariance Matrix / Selected Climo - External
     REAL, DIMENSION(:)  , allocatable       :: xbakg_clim,xbakg_emissdb
     REAL, DIMENSION(:,:), allocatable       :: ut_cov_clim,ut_cov_ext,ut_cov_emissdb
     !---External Data Background Vector:
     REAL, DIMENSION(:),   allocatable       :: xExt
     REAL, DIMENSION(:),   allocatable       :: emRfExt
     REAL, DIMENSION(:,:), allocatable       :: cvtIR
     !---Background Matrix / Selected Climo - External
     !---Indices
     TYPE(StateIndex_t)             :: IG,NG,IC,NC,IR,NR
     INTEGER                        :: nTranMW,nTranIR

     ! loadDone implementation provides the option to call load*() in advance
     logical :: loadBkgTuningDone=.false.
     logical :: loadBkgAtmosDone=.false.
     logical :: loadBkgSfcMWdone=.false.
     logical :: loadBkgSfcIRdone=.false.
     logical :: extAtmosOn = .false.
     logical :: extEmRfOn = .false.

   CONTAINS
     !---------------------------------------------------------------
     !  List of Public subroutines (accessible from outside module) 
     !---------------------------------------------------------------
     PROCEDURE, PASS :: bkgInitDefault
     PROCEDURE, PASS :: bkgInitAtmos
     PROCEDURE, PASS :: bkgInitSfcMW
     PROCEDURE, PASS :: bkgInitSfcIR
     PROCEDURE, PASS :: bkgAllocate
     PROCEDURE, PASS :: loadBkgAtmos
     PROCEDURE, PASS :: loadBkgSfcMW
     PROCEDURE, PASS :: loadBkgSfcIR
     PROCEDURE, PASS :: loadBkgTuning
     PROCEDURE, PASS :: setBkgClimAtmos
     PROCEDURE, PASS :: setBkgClimSfcMW
     PROCEDURE, PASS :: setBkgClimSfcIR
     PROCEDURE, PASS :: setPropertiesBkgTuning
     PROCEDURE, PASS :: checkTransformConform
     PROCEDURE, PASS :: combineBkgAtmos
     PROCEDURE, PASS :: combineBkgSfcMW
     PROCEDURE, PASS :: combineBkgSfcIR
     PROCEDURE, PASS :: overrideBkg
     PROCEDURE, PASS :: transformBkg
     PROCEDURE, PASS :: setRetrVec
     PROCEDURE, PASS :: setBkgExt
     PROCEDURE, PASS :: setBkgSfcIR
     PROCEDURE, PASS :: validateExtAtmos
     PROCEDURE, PASS :: validateExtSfcIR
     PROCEDURE, PASS :: setCloudMode
     PROCEDURE, PASS :: bkgDestroy !to be populated depending on dynamic fields
  END TYPE background_t

  !---------------------------------------------------------------
  !  List of Public structures (accessible from outside module) 
  !  The structure retFlags_t is used to indicate whether or not
  !  an element such as temp, cldLiq, etc. should be retrieved.
  !  This structure of flags is intended to be passed in as a
  !  parameter to the subroutine IRsetRetrVec.
  !---------------------------------------------------------------
CONTAINS


  SUBROUTINE bkgInitDefault(this)

    !<f90Subroutine>********************************************************
    !
    ! NAME:
    !
    !   bkgInitDefault
    !
    ! PURPOSE:
    !
    !   Initialize variables that are needed but might otherwise not 
    !      get initialized
    !
    ! SYNTAX:
    !
    !   call bkgInitDefault(this)
    !
    ! ARGUMENTS:
    !
    !   INPUTS/OUTPUTS:
    !   
    !   None
    !
    !
    ! INCLUDES:
    !
    !   None
    !
    !*******************************************************</f90Subroutine>

    CLASS(background_t), INTENT(out) :: this

    !local
    character (len=mxLength) :: procName = "[BackgroundModule::bkgInitDefault]"

    this%nEmMwG = 0     ! needed by bkgInitSfcIR to set iEmIRG index
    this%iemmwG = 0     ! needed by bkgInitSfcIR to set iEmIRG index
    this%nTranMW = 0    ! needed by setRetrVec to set iEmIrC index
    this%nParGSfcMW = 0 ! needed by combineBkgSfcIR as an offset

  END SUBROUTINE bkgInitDefault

  !--------------------------------------------------------------
  ! Get the necessary parameters and allocate memory for vectors 
  ! such as the background vectors and covariance matrices in
  ! geophysical space in preparation for retreival of IR emissivities.
  !--------------------------------------------------------------
  SUBROUTINE bkgInitAtmos(this,nParGout,IGout,NGout,nlevel, &
       vCoordTyp,pref,molID,molTran,limitCfg,bkgCfg)

    !<f90Subroutine>********************************************************
    !
    ! NAME:
    !
    !   bkgInitAtmos
    !
    ! PURPOSE:
    !
    !   Load and initialize atmospheric data, and allocate memory
    !
    ! SYNTAX:
    !
    !   call bkgInitAtmos(this,nParGout,IGout,NGout,nlevel, &
    !                     vCoordTyp,pref,molID,molTran,limitCfg,bkgCfg)
    !
    ! ARGUMENTS:
    !
    !   INPUTS/OUTPUTS:
    !   
    !   nParGout   INTEGER             Total number of elements in geophysical 
    !                                  state vector for output 
    !   IGout      TYPE(STATEINDEX_T)  Starting indices for sections of 
    !                                  geophysical state vector for output 
    !   NGout      TYPE(STATEINDEX_T)  Number of elements for sections of 
    !                                  geophysical state vector for output 
    !   nlevel     INTEGER             Number of atmospheric levels
    !   pref       REAL                Standard pressure profile
    !   molID      INTEGER             List of IDs of relevant molecular species
    !   molTran
    !   limitCfg   character           tuning parameter file
    !   bkgCfg     character           atmospheric covariance file
    !
    !   * OPTIONAL
    !
    ! INCLUDES:
    !
    !   None
    !
    !*******************************************************</f90Subroutine>

    CLASS(background_t), INTENT(out) :: this

    !--I/O variables
    integer,               intent(inout) :: nParGout
    type(StateIndex_t),    intent(inout) :: IGout,NGout  !---Indices
    integer,               intent(inout) :: nLevel
    character(LEN=mxCoordTyp),intent(inout) :: vCoordTyp
    real,    dimension(:), allocatable       :: pref
    integer, dimension(:), intent(inout) :: molID
    integer, dimension(:), intent(inout) :: molTran
    character (len=*), intent(in) :: limitCfg
    character (len=*), intent(in) :: bkgCfg
    !local
    character (len=mxLength) :: procName = "[BackgroundModule::bkgInitAtmos]"

    this%iH2O=whereH2O(molID)
    !------------------------------------------------------------------------
    ! Set global variables from file variables
    !------------------------------------------------------------------------

    ! Load bkgtuning parameters
    call this%loadBkgTuning(limitCfg)

    ! Load Atmospheric covariance
    call this%loadBkgAtmos(bkgCfg,molID,molTran,nLevel,vCoordTyp,pref)

    !------------------------------------------------------------------------
    ! Copy just the atmospheric part to the output arguments
    !------------------------------------------------------------------------
    nParGout      = this%nParGatm
    IGout         = this%IG
    NGout         = this%NG

  END SUBROUTINE bkgInitAtmos

  SUBROUTINE bkgInitSfcMW(this,nchMW,freq,polarity,iEmMwGout,bkgCfg)

    !<f90Subroutine>********************************************************
    !
    ! NAME:
    !
    !   bkgInitSfcMW
    !
    ! PURPOSE:
    !
    !   Load MW surface fields
    !
    ! SYNTAX:
    !
    !   call bkgInitSfcMW(this,nchMW,freq,polarity,nEmIRG,frqEmIR,iEmMwGout, &
    !                     iEmIRGout,bkgCfg)
    !
    ! ARGUMENTS:
    !
    !   INPUTS/OUTPUTS:
    !   
    !   nchmw      INTEGER             Number of MW channels
    !   freq       REAL                Frequency
    !   polarity   INTEGER             Polarization of MW channels
    !   nemIRG     INTEGER             Number of IR emissivity hinge points
    !   frqEmIR    REAL                Wavenumbers at IR emissivity hinge points
    !   iemmwGout  INTEGER             Starting index for MW emissivity in 
    !                                  geophysical state vector for output 
    !   iEmIRGout  INTEGER             Starting index for IR emissivity hinge 
    !                                  points for output 
    !   bkgCfg     character           surface field covariance file
    !
    !   * OPTIONAL
    !
    ! INCLUDES:
    !
    !   None
    !
    !*******************************************************</f90Subroutine>

    CLASS(background_t), intent(inout) :: this
    !--I/O variables
    integer,               intent(inout) :: nchmw
    real,    dimension(:), allocatable, intent(inout) :: freq
    integer, dimension(:), allocatable, intent(inout) :: polarity
    integer,               intent(inout) :: iemmwGout
    character (len=*),        intent(in) :: bkgCfg

    !  Load Surface covariance
    call this%loadBkgSfcMW(bkgCfg)

    !------------------------------------------------------------------------
    ! Compute global variables from file and input argument values
    !------------------------------------------------------------------------
    this%iemmwG        = this%IG%wind+this%NG%wind

    !------------------------------------------------------------------------
    ! Copy to output arguments
    !------------------------------------------------------------------------
    allocate(freq(this%nEmMwG),polarity(this%nEmMwG))
    nchmw         = this%nEmMwG
    freq          = this%Freq
    polarity      = this%Pol
    iemmwGout     = this%iemmwG

  END SUBROUTINE bkgInitSfcMW

  SUBROUTINE bkgInitSfcIR(this,nemIRG,frqEmIR,iEmIRGout,bkgCfg)

    !<f90Subroutine>********************************************************
    !
    ! NAME:
    !
    !   bkgInitSfcIR
    !
    ! PURPOSE:
    !
    !   Load and initialize atmospheric background data, and allocate memory
    !
    ! SYNTAX:
    !
    !   CALL bkgInitSfcIR(nParG)
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

    CLASS(background_t) :: this
    !--I/O variables
    integer,               intent(inout) :: nemIRG 
    real, dimension(:), allocatable, intent(inout) :: frqEmIR
    integer,               intent(inout) :: iEmIRGout
    character (len=*),     intent(in) :: bkgCfg

    !------------------------------------------------------------------------
    ! Set global variables from file variables
    !------------------------------------------------------------------------

    !  Load Surface covariance     (assumed to be the 2nd in covarFile(1:2))
    call this%loadBkgSfcIR(bkgCfg)

    !------------------------------------------------------------------------
    ! Copy to output arguments
    !------------------------------------------------------------------------
    allocate(frqEmIR(this%NemIrG))
    nemIRG        = this%NemIrG
    frqEmIR       = this%FrqEmIR
    iEmIRGout     = this%iemmwG+this%nEmMwG ! Starting index for IR emissivity

  END SUBROUTINE bkgInitSfcIR

  SUBROUTINE bkgAllocate(this, nParG)

    !<f90Subroutine>********************************************************
    !
    ! NAME:
    !
    !   bkgAllocate
    !
    ! PURPOSE:
    !
    !   Allocate all the background fields
    !
    ! SYNTAX:
    !
    !   CALL bkgAllocate(nParG)
    !
    ! ARGUMENTS:
    !
    !   INPUTS/OUTPUTS:
    !   
    !   nParG   INTEGER    Total number of elements in geophysical state vector
    !
    ! INCLUDES:
    !
    !   None
    !
    !*******************************************************</f90Subroutine>

    CLASS(background_t), intent(inout) :: this
    integer, intent(in) :: nParG

    !------------------------------------------------------------------------
    ! Allocate workspace. 
    !------------------------------------------------------------------------
    allocate(this%xbakg_clim(nParG))
    allocate(this%xbakg_emissdb(nParG))
    allocate(this%ut_cov_clim(nParG,nParG))
    allocate(this%ut_cov_ext(nParG,nParG))
    allocate(this%ut_cov_emissdb(nParG,nParG))
    allocate(this%umtx_clim(nParG,nParG))
    allocate(this%vmtx_clim(nParG,nParG))
    allocate(this%extrvec(nParG))
    allocate(this%tmpl(nParG,nParG))
    this%nParG = nParG

  END SUBROUTINE bkgAllocate

  !--
  ! Set global variables from argument list
  !--
  subroutine setPropertiesBkgTuning(this,backcldtopoc_in,backcldtopld_in, &
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
    !   Set global variables from argument list, for  tuning parameters
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

    CLASS(background_t), intent(inout) :: this
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
    this%backcldtopoc   = backcldtopoc_in
    this%backcldtopld   = backcldtopld_in
    this%backcldthkoc   = backcldthkoc_in
    this%backcldthkld   = backcldthkld_in
    this%backcldamtoc   = backcldamtoc_in
    this%backcldamtld   = backcldamtld_in
    this%backclddeffoc  = backclddeffoc_in
    this%backclddeffld  = backclddeffld_in
    this%varcldtopoc    = varcldtopoc_in
    this%varcldtopld    = varcldtopld_in
    this%varcldthkoc    = varcldthkoc_in
    this%varcldthkld    = varcldthkld_in
    this%varcldamtoc    = varcldamtoc_in
    this%varcldamtld    = varcldamtld_in
    this%varclddeffoc   = varclddeffoc_in
    this%varclddeffld   = varclddeffld_in
    this%backicetopoc   = backicetopoc_in
    this%backicetopld   = backicetopld_in
    this%backicethkoc   = backicethkoc_in
    this%backicethkld   = backicethkld_in
    this%backiceamtoc   = backiceamtoc_in
    this%backiceamtld   = backiceamtld_in
    this%backicedeffoc  = backicedeffoc_in
    this%backicedeffld  = backicedeffld_in
    this%varicetopoc    = varicetopoc_in
    this%varicetopld    = varicetopld_in
    this%varicethkoc    = varicethkoc_in
    this%varicethkld    = varicethkld_in
    this%variceamtoc    = variceamtoc_in
    this%variceamtld    = variceamtld_in
    this%varicedeffoc   = varicedeffoc_in
    this%varicedeffld   = varicedeffld_in
    this%varTskinOc     = varTskinOc_in
    this%varTskinLd     = varTskinLd_in

    this%loadBkgTuningDone=.true.
    return
  end subroutine setPropertiesBkgTuning

  !--
  ! Read items from namelist into local variables
  !--
  subroutine loadBkgTuning(this,file_algcfg)

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

    CLASS(background_t), intent(inout) :: this
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

    if(this%loadBkgTuningDone) return

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
    call this%setPropertiesBkgTuning(backcldtopoc,backcldtopld,backcldthkoc, &
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

  subroutine loadBkgAtmos(this,file_covar,molID,molTran,nlevel,vCoordTyp,pref)

    !<f90Subroutine>********************************************************
    !
    ! NAME:
    !
    !   loadBkgAtmos
    !
    ! PURPOSE:
    !
    !   Load atmospheric background and covariance data
    !
    ! SYNTAX:
    !
    !   CALL loadBkgAtmos(file_covar,molID,molTran,nlevel,vCoordTyp,pref)
    !
    ! ARGUMENTS:
    !
    !   INPUTS:
    !   
    !   file_covar  CHAR     Input background/covariance file name
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

    CLASS(background_t), intent(inout) :: this
    !--I/O variables
    character(len=*),          intent(in)    :: file_covar
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
    character (len=mxLength) :: procName = '[BackgroundModule::loadBkgAtmos]:'

    if(this%loadBkgAtmosDone) return

    !------------------------------------------------------------------------
    !  Atmospheric covariance (assumed to be the 1st in covarFile)
    !------------------------------------------------------------------------
    call queryCov(file_covar,nLevel=nlevel,IG=this%IG,NG=this%NG, &
         invNotTrnsp=invNotTrnsp,vCoordTyp=vCoordTyp)
    allocate(pref(nlevel),vCoordA(nlevel),vCoordB(nlevel),this%Pres(nlevel))
    this%nMol=getnmol(this%NG)
    this%nParGAtm = getVectorLength(this%NG)
    this%nParGsum = this%nParGAtm
    if (invNotTrnsp > 0) then
       call queryCov(file_covar,molid=molID,IC=this%IC,NC=this%NC)
       this%nParAtm=getVectorLength(this%NC)
       this%trnFixedAtm=.TRUE.  ! transformation cannot be truncated
    else
       this%IC=this%IG
       this%NC=this%NG
       this%nParAtm=this%nParGAtm
       this%trnFixedAtm=.FALSE.  ! transformation can be truncated (EOF)
    endif
    this%nParSum=this%nParAtm

    SELECT CASE (TRIM(vCoordTyp))
    CASE (Pcoord)
       CALL queryCov(file_covar,pressure=pref)
       this%Pres=pref
    CASE (Scoord)
       CALL queryCov(file_covar,sigCoord=vCoordA,sigPtop=sigPtop)
       CALL putSigmDefin(nlevel,sigPtop,vCoordA)
    CASE (Hcoord)
       CALL queryCov(file_covar,hybCoordA=vCoordA,hybCoordB=vCoordB)
       CALL putHybrDefin(nlevel,vCoordA,vCoordB)
    CASE default
       print*,trim(procName)//'Error: ',&
            'Unrecognized vertical coordinate type in cov file: ', &
            TRIM(vCoordTyp)
       CALL errorHalt(1)
    END SELECT

    U_cov = getUnit()
    CALL openCov(ncid=U_cov, file=file_covar,nbkg=this%nAtmTypes, molid=molID,&
         MolTran=molTran,nParG=this%nParGatm,IG=this%IG,NG=this%NG)
    allocate(this%back_clim_atm_all(this%nParGAtm,this%nAtmTypes))
    allocate(this%u_clim_atm_all(this%nParAtm,this%nParGAtm,this%nAtmTypes))
    allocate(this%v_clim_atm_all(this%nParGAtm,this%nParAtm,this%nAtmTypes))
    allocate(this%ut_cov_clim_atm_all(this%nParGAtm,this%nParGAtm,this%nAtmTypes))
    allocate(uv_atm_tmp(this%nParGAtm,this%nParGAtm))
    allocate(this%AtmTypes(this%nAtmTypes))


    DO itype = 1, this%nAtmTypes
       CALL getCov(ncid=U_cov,irec=itype, &
            dmean=this%back_clim_atm_all(1:this%nParGatm,itype),             &
            un_atm=uv_atm_tmp(1:this%nParGatm,1:this%nParAtm), &
            cov=this%ut_cov_clim_atm_all(1:this%nParGatm,1:this%nParGatm,itype),&
            TypeFlag=this%AtmTypes(itype))
       this%u_clim_atm_all(1:this%nParAtm,1:this%nParGatm,itype)= &
            TRANSPOSE(uv_atm_tmp(1:this%nParGatm,1:this%nParAtm))
       if (invNotTrnsp > 0) then
          CALL getCov(ncid=U_cov,irec=itype, &
               vn_atm=uv_atm_tmp(1:this%nParAtm,1:this%nParGatm))
          this%v_clim_atm_all(1:this%nParGatm,1:this%nParAtm,itype)= &
               TRANSPOSE(uv_atm_tmp(1:this%nParAtm,1:this%nParGatm))
       else
          this%v_clim_atm_all(1:this%nParGatm,1:this%nParAtm,itype)= &
               TRANSPOSE(this%u_clim_atm_all(1:this%nParAtm,1:this%nParGatm,itype))
       endif
       if (this%AtmTypes(itype)<0) then
          print*,trim(procName)//'Error: ',&
               ' AtmTypes undefined in Atm Covariance'
          call errorHalt(1)
       endif
    ENDDO
    CALL closeCov(ncid=U_cov)

    deallocate(uv_atm_tmp)

    this%loadBkgAtmosDone=.true.
    return
  end subroutine loadBkgAtmos

  subroutine loadBkgSfcMW(this,file_covar)

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

    CLASS(background_t), intent(inout) :: this
    !--I/O variables
    character(len=*), intent(in)    :: file_covar

    !---Local variables
    integer                              :: itype
    integer                              :: nchmw
    integer                              :: invNotTrnsp
    real, dimension(:,:), allocatable    :: uv_sfc_tmp
    real                                 :: sigPtop
    integer                              :: U_cov
    character (len=mxLength) :: procName = '[BackgroundModule::loadBkgSfcMW]'

    if(this%loadBkgSfcMWdone) return

    !------------------------------------------------------------------------
    !      Surface covariance
    !------------------------------------------------------------------------
    call queryCov(file=file_covar,nchmw=nchmw,invNotTrnsp=invNotTrnsp)
    allocate(this%Freq(nchmw),this%Pol(nchmw))

    if (invNotTrnsp > 0) then
       call queryCov(file_covar,nTranMW=this%nTranMW)
       this%trnFixedSfcMW=.TRUE.  ! transformation cannot be truncated
    else
       this%nTranMW=nchmw
       this%trnFixedSfcMW=.FALSE.  ! transformation can be truncated (EOF)
    endif
    this%nParGsfcMW=nchmw
    this%sGsfcMW=this%nParGsum+1
    this%eGsfcMW=this%nParGsum+this%nParGsfcMW
    this%nParGsum=this%eGsfcMW
    this%nParSfcMW=this%nTranMW
    this%sSfcMW=this%nParSum+1
    this%eSfcMW=this%nParSum+this%nParSfcMW
    this%nParSum=this%eSfcMW

    U_cov = getUnit()
    CALL openCov(ncid=U_cov,file=file_covar,nbkg=this%nSfcTypesMW, &
         freq=this%Freq(1:nchmw),polarity=this%Pol(1:nchmw))

    allocate(this%back_clim_sfc_mw(this%nParGsfcMW,this%nSfcTypesMW))
    allocate(this%u_clim_sfc_mw(this%nParSfcMW,this%nParGsfcMW,this%nSfcTypesMW))
    allocate(this%v_clim_sfc_mw(this%nParGsfcMW,this%nParSfcMW,this%nSfcTypesMW))
    allocate(this%ut_cov_clim_sfc_mw(this%nParGsfcMW,this%nParGsfcMW,this%nSfcTypesMW))
    allocate(uv_sfc_tmp(this%nParGsfcMW,this%nParGsfcMW))
    allocate(this%SfcTypesMW(this%nSfcTypesMW))
    this%u_clim_sfc_mw(:,:,:)=0.
    this%v_clim_sfc_mw(:,:,:)=0.
    uv_sfc_tmp(:,:)=0.
    this%ut_cov_clim_sfc_mw(:,:,:)=0.
    DO itype = 1, this%nSfcTypesMW
       CALL getCov(ncid=U_cov,irec=itype,&
            dmEmMw=this%back_clim_sfc_mw(1:nchmw,itype),&
            un_emmw=uv_sfc_tmp(1:nchmw,1:this%nTranMW),&
            covEmMw=this%ut_cov_clim_sfc_mw(1:nchmw,1:nchmw,itype),&
            TypeFlag=this%SfcTypesMW(itype))
       this%u_clim_sfc_mw(1:this%nParSfcMW,1:this%nParGsfcMW,itype)= &
            TRANSPOSE(uv_sfc_tmp(1:this%nParGsfcMW,1:this%nParSfcMW))
       if (invNotTrnsp > 0) then
          CALL getCov(ncid=U_cov,irec=itype,&
               vn_emmw=uv_sfc_tmp(1:this%nTranMW,1:nchmw))
          this%v_clim_sfc_mw(1:this%nParGsfcMW,1:this%nParSfcMW,itype)= &
               TRANSPOSE(uv_sfc_tmp(1:this%nParSfcMW,1:this%nParGSfcMW))
       else
          this%v_clim_sfc_mw(1:this%nParGSfcMW,1:this%nParSfcMW,itype)= &
               TRANSPOSE(this%u_clim_sfc_mw(1:this%nParSfcMW,1:this%nParGSfcMW,itype))
       endif
       if (this%SfcTypesMW(itype) < 0) then
          print*,trim(procName)//'Error: ',&
               ' MW SfcTypes undefined in Sfc Covariance'
          call errorHalt(1)
       endif
    ENDDO

    CALL closeCov(ncid=U_cov)

    !------------------------------------------------------------------------
    ! Set global variables from file variables
    !------------------------------------------------------------------------
    this%nEmMwG        = nchmw

    deallocate(uv_sfc_tmp)

    this%loadBkgSfcMWdone=.true.
    return
  end subroutine loadBkgSfcMW

  subroutine loadBkgSfcIR(this,file_covar)

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

    CLASS(background_t), intent(inout) :: this
    !--I/O variables
    character(len=*),    intent(in)    :: file_covar
    !---Local variables
    integer                              :: itype
    integer                              :: nchmw
    integer                              :: invNotTrnsp
    real, dimension(:,:), allocatable    :: uv_sfc_tmp
    real                                 :: sigPtop
    real, dimension(:),   allocatable    :: vCoordA,vCoordB
    integer :: U_cov
    character (len=mxLength) :: procName = '[BackgroundModule::loadBkgSfcIR]'

    if(this%loadBkgSfcIRdone) return

    !------------------------------------------------------------------------
    !      Surface covariance (assumed to be the 2nd in covarFile)
    !------------------------------------------------------------------------
    call queryCov(file=file_covar,nemir=this%NemirG,invNotTrnsp=invNotTrnsp)
    allocate(this%FrqEmIR(this%NemirG))

    if (invNotTrnsp > 0) then
       call queryCov(file_covar,nTranIR=this%nTranIR)
       this%trnFixedSfcIR=.TRUE.  ! transformation cannot be truncated
    else
       this%nTranIR=2*this%NemirG
       this%trnFixedSfcIR=.FALSE.  ! transformation can be truncated (EOF)
    endif
    ! nParGsfc=nchmw+2*myNemirG
    this%nParGsfcIR=2*this%NemirG
    this%sGsfcIR=this%nParGsum+1
    this%eGsfcIR=this%nParGsum+this%nParGsfcIR
    this%nParGsum=this%eGsfcIR
    this%nParSfcIR=this%nTranIR
    this%sSfcIR=this%nParSum+1
    this%eSfcIR=this%nParSum+this%nParSfcIR
    this%nParSum=this%eSfcIR

    U_cov = getUnit()
    CALL openCov(ncid=U_cov,file=file_covar,nbkg=this%nSfcTypesIR, &
         frqemir=this%FrqEmIR(1:this%NemirG))

    allocate(this%back_clim_sfc_ir(this%nParGsfcIR,this%nSfcTypesIR))
    allocate(this%u_clim_sfc_ir(this%nParSfcIR,this%nParGsfcIR,this%nSfcTypesIR))
    allocate(this%v_clim_sfc_ir(this%nParGsfcIR,this%nParSfcIR,this%nSfcTypesIR))
    allocate(this%ut_cov_clim_sfc_ir(this%nParGsfcIR,this%nParGsfcIR,this%nSfcTypesIR))
    allocate(uv_sfc_tmp(this%nParGsfcIR,this%nParGsfcIR))
    allocate(this%SfcTypesIR(this%nSfcTypesIR))
    this%u_clim_sfc_ir(:,:,:)=0.
    this%v_clim_sfc_ir(:,:,:)=0.
    uv_sfc_tmp(:,:)=0.
    this%ut_cov_clim_sfc_ir(:,:,:)=0.
    DO itype = 1, this%nSfcTypesIR
       CALL getCov(ncid=U_cov,irec=itype,&
            dmEmIR=this%back_clim_sfc_ir(1:this%nParGsfcIR,itype),&
            un_emIR=uv_sfc_tmp(1:this%nParGsfcIR,1:this%nParSfcIR),&
            covEmIR=this%ut_cov_clim_sfc_ir(1:this%nParGsfcIR,1:this%nParGsfcIR,itype),&
            TypeFlag=this%SfcTypesIR(itype))
       this%u_clim_sfc_ir(1:this%nParSfcIR,1:this%nParGsfcIR,itype)= &
            TRANSPOSE(uv_sfc_tmp(1:this%nParGsfcIR,1:this%nParSfcIR))
       if (invNotTrnsp > 0) then
          CALL getCov(ncid=U_cov,irec=itype,&
               vn_emIR=uv_sfc_tmp(1:this%nParSfcIR,1:this%nParGsfcIR))
          this%v_clim_sfc_ir(1:this%nParGsfcIR,1:this%nParSfcIR,itype)= &
               TRANSPOSE(uv_sfc_tmp(1:this%nParSfcIR,1:this%nParGsfcIR))
       else
          this%v_clim_sfc_ir(1:this%nParGsfcIR,1:this%nParSfcIR,itype)= &
               TRANSPOSE(this%u_clim_sfc_ir(1:this%nParSfcIR,1:this%nParGsfcIR,itype))
       endif
       if (this%SfcTypesIR(itype) < 0) then
          print*,trim(procName)//'Error: ',&
               ' IR SfcTypes undefined in Sfc Covariance'
          call errorHalt(1)
       endif
    ENDDO

    CALL closeCov(ncid=U_cov)

    deallocate(uv_sfc_tmp)

    this%loadBkgSfcIRdone=.true.
    return
  end subroutine loadBkgSfcIR

  SUBROUTINE checkTransformConform(this,NR)

    CLASS(background_t), intent(in) :: this

    !---Input variables
    TYPE(StateIndex_t), INTENT(IN) :: NR
    !---Local variables
    INTEGER :: j
    character (len=mxLength) :: &
       procName = '[BackgroundModule::checkTransformConform]'

    !--- Verify that the selected number of variables to retrieve are
    !--- compatible with the transformation matrices
    IF (this%trnFixedAtm) THEN
       IF ((NR%temp /= 0) .AND. (NR%temp /= this%NC%temp)) THEN
          print*,trim(procName)//'Error: ',&
               'Cannot truncate transform of temperature profile'
          call errorHalt(1)
       ENDIF
       IF ((NR%tskin /= 0) .AND. (NR%tskin /= this%NC%tskin)) THEN
          print*,trim(procName)//'Error: ',&
               'Cannot truncate transform of surface temperature'
          call errorHalt(1)
       ENDIF
       IF ((NR%psfc /= 0) .AND. (NR%psfc /= this%NC%psfc)) THEN
          print*,trim(procName)//'Error: ',&
               'Cannot truncate transform of surface temperature'
          call errorHalt(1)
       ENDIF
       DO j=1,this%nMol
          IF ((NR%mol(j) /= 0) .AND. (NR%mol(j) /= this%NC%mol(j))) THEN
             print*,trim(procName)//'Error: ',&
                  'Cannot truncate transform of gas profile',j
             call errorHalt(1)
          ENDIF
       ENDDO
       IF ((NR%cldLiq /= 0) .AND. (NR%cldLiq /= this%NC%cldLiq)) THEN
          print*,trim(procName)//'Error: ',&
               'Cannot truncate transform of cloud liquid profile'
          call errorHalt(1)
       ENDIF
       IF ((NR%cldIce /= 0) .AND. (NR%cldIce /= this%NC%cldIce)) THEN
          print*,trim(procName)//'Error: ',&
               'Cannot truncate transform of cloud ice profile'
          call errorHalt(1)
       ENDIF
       IF ((NR%wind /= 0) .AND. (NR%wind /= this%NC%wind)) THEN
          print*,trim(procName)//'Error: ',&
               'Cannot truncate transform of wind'
          call errorHalt(1)
       ENDIF
    ENDIF

    IF (this%loadBkgSfcMWdone) THEN
       IF (this%trnFixedSfcMW) THEN
          IF ((NR%emMW /= 0) .AND. (NR%emMW /= this%nTranMW)) THEN
             print*,trim(procName)//'Error: ',&
                  'Cannot truncate transform of MW emissivity'
             call errorHalt(1)
          ENDIF
       ENDIF
    ENDIF

    IF (this%loadBkgSfcIRdone) THEN
       IF (this%trnFixedSfcIR) THEN
          IF ((NR%emRfIR /= 0) .AND. (NR%emRfIR /= this%nTranIR)) THEN
             print*,trim(procName)//'Error: ',&
                  'Cannot truncate transform of IR emissivity'
             call errorHalt(1)
          ENDIF
       ENDIF
    ENDIF
  END SUBROUTINE checkTransformConform

  SUBROUTINE setBkgClimAtmos(this,iclassatm,psfc,vCoordTyp,press,umtx,vmtx)

    !<f90Subroutine>********************************************************
    !
    ! NAME:
    !
    !   setBkgClimAtmos (IRgetBkgClim)
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

    CLASS(background_t), intent(inout) :: this

    !---Input variables
    INTEGER,                 INTENT(IN)    :: iclassatm
    REAL,                    INTENT(IN)    :: psfc
    REAL,    DIMENSION(:,:), INTENT(INOUT) :: umtx
    REAL,    DIMENSION(:,:), INTENT(INOUT) :: vmtx
    REAL,    DIMENSION(:),   INTENT(INOUT) :: press
    CHARACTER(LEN=mxCoordTyp),INTENT(IN)   :: vCoordTyp
    character (len=mxLength) :: procName = '[BackgroundModule::setBkgClimAtmos]:'

    if (iclassatm > this%nAtmTypes) then
       print *, trim(procName)//'Error: igeo > nAtmTypes ', &
            iclassatm,this%nAtmTypes
       call errorHalt(1)
    endif

    !----Background selection
    this%xbakg_clim(1:this%nParGAtm)        = this%back_clim_atm_all(1:this%nParGAtm,iclassatm)

    !----Covariance matrix selection
    this%ut_cov_clim(1:this%nParGAtm,1:this%nParGAtm)=0.
    !--atmosphere
    this%ut_cov_clim(1:this%nParGAtm,1:this%nParGAtm)= &
         this%ut_cov_clim_atm_all(1:this%nParGAtm,1:this%nParGAtm,iclassatm)

    !----Transformation matrix selection
    this%umtx_clim(1:this%nParAtm,1:this%nParGatm)=0.
    this%umtx_clim(1:this%nParAtm,1:this%nParGatm) = &
         this%u_clim_atm_all(1:this%nParAtm,1:this%nParGAtm,iclassatm)
    this%vmtx_clim(1:this%nParGatm,1:this%nParAtm)=0.
    this%vmtx_clim(1:this%nParGatm,1:this%nParAtm) = &
         this%v_clim_atm_all(1:this%nParGAtm,1:this%nParAtm,iclassatm)

    umtx = this%umtx_clim(this%extrvec(1:this%nParRatm),1:this%nParGatm)
    vmtx = this%vmtx_clim(1:this%nParGatm,this%extrvec(1:this%nParRatm))
    SELECT CASE (TRIM(vCoordTyp))
    CASE (Pcoord)
       press=this%Pres
    CASE (Scoord)
       press=getSigmPres(psfc)
    CASE (Hcoord)
       press=getHybrPres(psfc)
    END SELECT
    RETURN
  END SUBROUTINE setBkgClimAtmos

  SUBROUTINE setBkgClimSfcMW(this,igeo,umtx,vmtx)

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

    CLASS(background_t), intent(inout) :: this

    !---Input variables
    INTEGER,                 INTENT(IN)    :: igeo
    REAL,    DIMENSION(:,:), INTENT(INOUT) :: umtx
    REAL,    DIMENSION(:,:), INTENT(INOUT) :: vmtx

    !---Local
    integer :: j
    integer :: sParG
    integer :: eParG
    integer :: sPar
    integer :: ePar
    integer :: sParR
    integer :: eParR
    character (len=mxLength) :: procName = '[BackgroundModule::setBkgClimSfcMW]:'

    if (igeo > this%nSfcTypesMW) then
       print *, trim(procName)//'Error: igeo > nSfcTypesMW ', &
            igeo,this%nSfcTypesMW
       call errorHalt(1)
    endif

    sParG = this%sGsfcMW
    eParG = this%eGsfcMW
    !----Background selection
    this%xbakg_clim(sParG:eParG) = this%back_clim_sfc_mw(1:this%nParGsfcMW,igeo)
    !----Covariance matrix selection
    this%ut_cov_clim(sParG:eParG,sParG:eParG)=0.
    !--surface
    this%ut_cov_clim(sParG:eParG,sParG:eParG)= &
         this%ut_cov_clim_sfc_mw(1:this%nParGsfcMW,1:this%nParGsfcMW,igeo)
    !----transformation matrix selection
    sPar = this%sSfcMW
    ePar = this%eSfcMW
    this%umtx_clim(sPar:ePar,sParG:eParG)=0.
    this%umtx_clim(sPar:ePar,sParG:eParG)= &
         this%u_clim_sfc_mw(1:this%nParSfcMW,1:this%nParGsfcMW,igeo)
    this%vmtx_clim(sParG:eParG,sPar:ePar)=0.
    this%vmtx_clim(sParG:eParG,sPar:ePar)= &
         this%v_clim_sfc_mw(1:this%nParGsfcMW,1:this%nParSfcMW,igeo)
    sParR = this%sRsfcMW
    eParR = this%eRsfcMW
    umtx = this%umtx_clim(this%extrvec(sParR:eParR),sParG:eParG)
    vmtx = this%vmtx_clim(sParG:eParG,this%extrvec(sParR:eParR))
  END SUBROUTINE setBkgClimSfcMW

  SUBROUTINE setBkgClimSfcIR(this,igeo,umtx,vmtx)

    !<f90Subroutine>********************************************************
    !
    ! NAME:
    !
    !   setBkgClimSfcIR
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

    CLASS(background_t), intent(inout) :: this

    !---Input variables
    INTEGER,                 INTENT(IN)    :: igeo
    REAL,    DIMENSION(:,:), INTENT(INOUT) :: umtx
    REAL,    DIMENSION(:,:), INTENT(INOUT) :: vmtx
    !---Local variables
    ! INTEGER :: j
    INTEGER :: sParG
    INTEGER :: eParG
    INTEGER :: sPar
    INTEGER :: ePar
    INTEGER :: sParR
    INTEGER :: eParR
    character (len=mxLength) :: procName = '[BackgroundModule::setBkgClimSfcIR]:'

    if (igeo > this%nSfcTypesIR) then
       print *, trim(procName)//'Error: igeo > nSfcTypesIR ', &
            igeo,this%nSfcTypesIR
       call errorHalt(1)
    endif

    sParG = this%sGsfcIR
    eParG = this%eGsfcIR
    !----Background selection
    this%xbakg_clim(sParG:eParG) = this%back_clim_sfc_ir(1:this%nParGsfcIR,igeo)

    !----Covariance matrix selection
    this%ut_cov_clim(sParG:eParG,sParG:eParG)=0.
    !--surface
    this%ut_cov_clim(sParG:eParG,sParG:eParG)= &
         this%ut_cov_clim_sfc_ir(1:this%nParGsfcIR,1:this%nParGsfcIR,igeo)

    !----Transformation matrix selection
    sPar = this%sSfcIR
    ePar = this%eSfcIR
    this%umtx_clim=0.
    this%umtx_clim(sPar:ePar,sParG:eParG) = this%u_clim_sfc_ir(1:this%nParSfcIR,1:this%nParGsfcIR,igeo)
    this%vmtx_clim=0.
    this%vmtx_clim(sParG:eParG,sPar:ePar) = this%v_clim_sfc_ir(1:this%nParGsfcIR,1:this%nParSfcIR,igeo)
    !----Truncate (compress) the transformation matrix
    !--This initialization is necessary on Linux
    !--because otherwise the implicit do loop will 
    !--begin at j=0 and cause a segmentation fault
    !--This needs to be further investigated!
    ! j=1

    !assume umtx and vmtx are IR only
    sParR = this%sRsfcIR
    eParR = this%eRsfcIR
    umtx = this%umtx_clim(this%extrvec(sParR:eParR),sParG:eParG)
    vmtx = this%vmtx_clim(sParG:eParG,this%extrvec(sParR:eParR))
  END SUBROUTINE setBkgClimSfcIR


  !-----------------------------------------------------------------
  ! The following subroutine is intended to be used to set 
  ! the o/p parameters such as IR_out, NR_out, nPar_out 
  ! using the values in corresponding input parameters and
  ! conditional on the retrievalFlags.
  ! Dynamically selected flags from pre-defined sets are used
  ! to determine which variable is to be set. 
  !----------------------------------------------------------------- 
  SUBROUTINE setRetrVec(this,retrievalFlags, NR_in, nPar_out, NR_out, IR_out)

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

    CLASS(background_t), intent(inout) :: this

    TYPE(retFlags_t),     INTENT(IN)   :: retrievalFlags
    TYPE(StateIndex_t),   INTENT(IN)   :: NR_in
    INTEGER,              INTENT(INOUT):: nPar_out
    TYPE(StateIndex_t),   INTENT(INOUT):: NR_out
    TYPE(StateIndex_t),   INTENT(INOUT):: IR_out

    ! some local variables
    INTEGER          :: i, j
    INTEGER          :: IEmMwC
    INTEGER          :: IEmIrC
    character (len=mxLength) :: procName = '[BackgroundModule::setRetrVec]:'

    this%NR  = NR_in

    ! While the parent subroutines do not set the surface parts of NR_in:
    ! NR%emMW = NR_in%emMW
    ! NR%emRfIR = 2 * NR_in%emRfIR

    IF (.NOT. retrievalFlags%tempFlag) THEN
       this%NR%temp = 0
    END IF

    WHERE (retrievalFlags%molFlags)
       this%NR%mol = NR_in%mol
    ELSEWHERE
       this%NR%mol = 0
    END WHERE

    IF (.NOT. retrievalFlags%cldLiqFlag) THEN
       this%NR%cldLiq = 0
    END IF

    IF (.NOT. retrievalFlags%cldIceFlag) THEN
       this%NR%cldIce = 0
    END IF

    IF (.NOT. retrievalFlags%psfcFlag ) THEN
       this%NR%psfc = 0
    END IF

    IF (.NOT. retrievalFlags%tskinFlag ) THEN
       this%NR%tskin = 0 
    END IF

    IF (.NOT. retrievalFlags%windFlag ) THEN
       this%NR%wind = 0
    END IF

    IF (.NOT. retrievalFlags%emissMwFlag ) THEN
       this%NR%emMW = 0
    END IF

    IF (.NOT. retrievalFlags%emissIrFlag ) THEN
       this%NR%emRfIR = 0
    END IF

    this%nParR = getVectorLength(this%NR)
    this%nParRatm = this%nParR - this%NR%emMW - this%NR%emRfIR
    this%sRsfcMW = this%nParRatm+1
    this%eRsfcMW = this%nParRatm+this%NR%emMW
    this%sRsfcIR = this%eRsfcMW+1
    this%eRsfcIR = this%eRsfcMW+this%NR%emRfIR
    npar_out = this%nParR
    this%IR = genIndices(this%NR)
    this%iemmw = this%IR%wind + this%NR%wind

    NR_out = this%NR
    IR_out = this%IR

    IEmMwC = this%IC%wind+this%NC%wind
    IEmIrC = IEmMwC+this%nTranMW

    j = 1
    this%extrvec(1:this%nParR)= &
         (/(i,i=this%IC%Temp,this%IC%Temp+this%NR%Temp-1), &
         (i,i=this%IC%Tskin,this%IC%Tskin+this%NR%Tskin-1), &
         (i, i=this%IC%Psfc,this%IC%Psfc+this%NR%Psfc-1), &
         ((i,i=this%IC%mol(j),this%IC%mol(j)+this%NR%mol(j)-1),j=1,this%nMol),   &
         (i,i=this%IC%cldLiq,this%IC%cldLiq+this%NR%cldLiq-1),              &
         (i,i=this%IC%cldIce,this%IC%cldIce+this%NR%cldIce-1),              &
         (i,i=IEmMwC,IEmMwC+this%NR%emMW-1),                 &
         (i,i=IEmIrC,IEmIrC+this%NR%emRfIR-1)/)
    RETURN
  END SUBROUTINE setRetrVec

  SUBROUTINE setBkgExt(this,xbakg_ext)

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

    CLASS(background_t), intent(inout) :: this

    REAL,DIMENSION(:),INTENT(in) :: xbakg_ext
    !---Load external background data into global vector 'xExt' which 
    !    will be accessed in IRcombBkg.  When accessing external data,
    !    this subroutine should be called immediately before calling 
    !    IRcombBkg.
    this%xExt = xbakg_ext
    this%extAtmosOn = .true.
    RETURN
  END SUBROUTINE setBkgExt

  SUBROUTINE setBkgSfcIR(this,xbakg_emis,nemirExt)

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

    CLASS(background_t), intent(inout) :: this

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
       this%emRfExt=MISSING_REAL
       RETURN
    ENDIF
    ! Interpolate spectrally, as needed
    IF (nemirExt /= this%NemIrG) THEN
       this%emRfExt(1:this%NemIrG)=matmul(this%cvtIR,xbakg_emis(1:nemirExt))
       this%emRfExt(this%NemIrG+1:this%NemIrG+this%NemIrG)=&
            matmul(this%cvtIR,xbakg_emis(nemirExt+1:nemirExt+nemirExt))
    ELSE
       this%emRfExt(1:this%NemIrG)=xbakg_emis(1:this%NemIrG)
       this%emRfExt(this%NemIrG+1:this%NemIrG+this%NemIrG)=&
            xbakg_emis(this%NemIrG+1:this%NemIrG+this%NemIrG)
    ENDIF
    DO i=1,this%NemIrG
       this%emRfExt(i)=min(this%emRfExt(i),emIRcap) 
       this%emRfExt(i)=max(this%emRfExt(i),emIRfloor) 
    ENDDO

    DO i=this%NemIrG+1,this%NemIrG+this%NemIrG
       this%emRfExt(i)=max(this%emRfExt(i),1.0-emIRcap)
       this%emRfExt(i)=min(this%emRfExt(i),1.0-emIRfloor)
    ENDDO

    this%extEmRfOn = .true.
    RETURN
  END SUBROUTINE setBkgSfcIR

  SUBROUTINE combineBkgAtmos(this,IGExt,NGExt,psfc,isLand,molTran, &
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
    CLASS(background_t), intent(in) :: this
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
    logical :: noExt
    character (len=mxLength) :: procName = '[BackgroundModule::combineBkgAtmos]:'

    !---By default load the climatology background into the worskpace array
    dback0(1:this%nParGatm)=this%xbakg_clim(1:this%nParGatm)
    ut_cov_net(1:this%nParGatm,1:this%nParGatm) = this%ut_cov_clim(1:this%nParGatm,1:this%nParGatm)
    ! dback0(1:this%nParG)=this%xbakg_clim(1:this%nParG)
    ! ut_cov_net(1:this%nParG,1:this%nParG) = this%ut_cov_clim(1:this%nParG,1:this%nParG)
    dback0(this%IG%Psfc)=psfc

    noExt = .TRUE.
    IF (getVectorLength(NGExt) > 0) THEN
       IF ( ANY(this%xExt /= MISSING_REAL) ) noExt=.FALSE.
    ENDIF
    IF ( noExt .AND. &
         ANY( (/ &
         this%externalFlags%xTemp  , this%externalFlags%xTskin  , &
         this%externalFlags%xPsfc  , this%externalFlags%xMols(:), &
         this%externalFlags%xCldLiq, this%externalFlags%xCldIce , &
         this%externalFlags%xWind /) ) ) THEN
       print*,trim(procName)//'Error:  ', &
            'External data requested, but external data vector does not ',&
            'contain any valid data.',&
            'Call IRsetExtBkg prior to calling IRcombBkg.'
       call errorHalt(1)
    ENDIF

    IF (this%externalFlags%xPsfc) then
       dback0(this%IG%Psfc)=this%xExt(IGExt%psfc)
    ENDIF

    IF (TRIM(vCoordTyp) == Pcoord) THEN
       !---background may be overwritten with external data
       !---LVL_INT anywhere below is called to extrapolate the values below surface
       !---and get rid of '-9999' in modification of BKG
       IF (this%externalFlags%xTemp) then
          call lvl_int(this%xExt(IGExt%temp:IGExt%temp+NGExt%temp-1),press, &
               NGExt%temp,this%xExt(IGExt%Psfc),xSfc,dxlvldp,Nsurf)
          deltExtrap=this%xExt(IGExt%temp+Nsurf-1)-dback0(this%IG%temp+Nsurf-1)
          dback0(this%IG%temp:this%IG%temp+Nsurf-1)=&
               this%xExt(IGExt%temp:IGExt%temp+Nsurf-1)
          !    For all levels below surface (levels >= Nsurf) we set a constant deviation
          !    from background for an external Temp profile = xExt(Nsurf)-dback0(Nsurf) 
          dback0(this%IG%temp+Nsurf:this%IG%temp+this%NG%temp-1)=&
               dback0(this%IG%temp+Nsurf:this%IG%temp+this%NG%temp-1)+deltExtrap
       ENDIF
       IF (this%externalFlags%xTskin) then
          dback0(this%IG%tskin)=this%xExt(IGExt%tskin)-xSfc ! Delta-Tskin with fixed grid
       ENDIF
       DO i=1,this%nMolExt
          ii=this%mapExtMols(i)
          IF (this%externalFlags%xMols(i)) then
             call lvl_int(this%xExt(IGExt%mol(i):IGExt%mol(i)+NGExt%mol(i)-1),press, &
                  NGExt%mol(i),this%xExt(IGExt%Psfc),xSfc,dxlvldp,Nsurf)
             !    The same way as for Temp profile, we do for gases profiles.
             if (molTran(ii) == MOL_TRAN_LOG) then
                deltExtrap=alog(this%xExt(IGExt%mol(i)+Nsurf-1))-dback0(this%IG%mol(ii)+Nsurf-1)
                dback0(this%IG%mol(ii):this%IG%mol(ii)+Nsurf-1)=&
                     alog(this%xExt(IGExt%mol(i):IGExt%mol(i)+Nsurf-1))
             else
                deltExtrap=this%xExt(IGExt%mol(i)+Nsurf-1)-dback0(this%IG%mol(ii)+Nsurf-1)
                dback0(this%IG%mol(ii):this%IG%mol(ii)+Nsurf-1)=&
                     this%xExt(IGExt%mol(i):IGExt%mol(i)+Nsurf-1)
             endif
             dback0(this%IG%mol(ii)+Nsurf:this%IG%mol(ii)+this%NG%mol(ii)-1)=&
                  dback0(this%IG%mol(ii)+Nsurf:this%IG%mol(ii)+this%NG%mol(ii)-1)+&
                  deltExtrap
          ENDIF
       ENDDO
    ELSE
       IF (this%externalFlags%xTemp) &
            dback0(this%IG%temp:this%IG%temp+this%NG%temp-1)=&
            this%xExt(IGExt%temp:IGExt%temp+NGExt%temp-1)
       IF (this%externalFlags%xTskin) &
            dback0(this%IG%tskin)=this%xExt(IGExt%tskin)
       DO i=1,this%nMolExt
          ii=this%mapExtMols(i)
          IF (this%externalFlags%xMols(i)) then
             if (molTran(ii) == MOL_TRAN_LOG) then
                dback0(this%IG%mol(ii):this%IG%mol(ii)+this%NG%mol(ii)-1)= &
                     alog(this%xExt(IGExt%mol(i):IGExt%mol(i)+NGExt%mol(i)-1))
             else
                dback0(this%IG%mol(ii):this%IG%mol(ii)+this%NG%mol(ii)-1)= &
                     this%xExt(IGExt%mol(i):IGExt%mol(i)+NGExt%mol(i)-1)
             endif
          ENDIF
       ENDDO
    ENDIF
    IF (this%externalFlags%xCldLiq) then
       dback0(this%IG%cldLiq:this%IG%cldLiq+this%NG%cldLiq-1)=&
            this%xExt(IGExt%cldLiq:IGExt%cldLiq+NGExt%cldLiq-1)
    ENDIF
    IF (this%externalFlags%xCldIce) then
       dback0(this%IG%cldIce:this%IG%cldIce+this%NG%cldIce-1)=&
            this%xExt(IGExt%cldIce:IGExt%cldIce+NGExt%cldIce-1)
    ENDIF
    IF (this%externalFlags%xWind) then
       dback0(this%IG%wind:this%IG%wind+this%NG%wind-1)=&
            this%xExt(IGExt%wind:IGExt%wind+NGExt%wind-1)
    ENDIF
    RETURN
  END SUBROUTINE combineBkgAtmos

  SUBROUTINE combineBkgSfcMW(this,isLand,dback0,ut_cov_net)

    !<f90Subroutine>********************************************************
    !
    ! NAME:
    !
    !   combineBkgSfcMW
    !
    ! PURPOSE:
    !
    !   Load climatology background into the worskpace array. Background may be 
    !   overwritten with external data.
    !
    ! SYNTAX:
    !
    !   CALL combineBkgSfcMW(pland, plandMaxOc, dback0, ut_cov_net)
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

    implicit none
    CLASS(background_t), intent(in) :: this
    LOGICAL,              INTENT(IN)    :: isLand
    REAL, DIMENSION(:),   INTENT(INOUT) :: dback0
    REAL, DIMENSION(:,:), INTENT(INOUT) :: ut_cov_net
    !---Local variables
    integer :: sParG
    integer :: eParG
    character (len=mxLength) :: procName = '[BackgroundModule::combineBkgSfcMW]:'

    !---By default load the climatology background into the worskpace array
    sParG = this%nParGatm+1
    eParG = this%nParGatm+this%nParGSfcMW
    !assume that back0 and ut_cov_net are MW only 
    dback0(1:this%nParGSfcMW)=this%xbakg_clim(sParG:eParG)
    ut_cov_net(1:this%nParGsfcMW,1:this%nParGsfcMW) = this%ut_cov_clim(sParG:eParG,sParG:eParG)

  END SUBROUTINE combineBkgSfcMW

  SUBROUTINE combineBkgSfcIR(this,isLand,dback0,ut_cov_net)

    !<f90Subroutine>********************************************************
    !
    ! NAME:
    !
    !   combineBkgSfcIR
    !
    ! PURPOSE:
    !
    !   Load climatology background into the worskpace array. Background may be 
    !   overwritten with external data.
    !
    ! SYNTAX:
    !
    !   CALL combineBkgSfcIR(pland, plandMaxOc, dback0, ut_cov_net)
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

    implicit none
    CLASS(background_t), intent(in) :: this
    LOGICAL,              INTENT(IN)    :: isLand
    REAL, DIMENSION(:),   INTENT(INOUT) :: dback0
    REAL, DIMENSION(:,:), INTENT(INOUT) :: ut_cov_net
    !---Local variables
    integer :: sParG
    integer :: eParG
    character (len=mxLength) :: procName = '[BackgroundModule::combineBkgSfcIR]:'

    !---By default load the climatology background into the worskpace array
    sParG = this%nParGatm+this%nParGSfcMW+1
    eParG = this%nParGatm+this%nParGSfcMW+this%nParGsfcIR
    dback0(1:this%nParGSfcIR)=this%xbakg_clim(sParG:eParG)
    ut_cov_net(1:this%nParGsfcIR,1:this%nParGsfcIR) = this%ut_cov_clim(sParG:eParG,sParG:eParG)

    IF (isLand .and. this%extAtmosOn .and. this%extEmRfOn) THEN
       IF (.NOT. ANY(this%emRfExt == MISSING_REAL)) THEN
          dback0(1:this%nParGSfcIR)=this%emRfExt
       ENDIF
    ELSE
       ! placeholder for use of model that depends on conditions (wind, etc.)
    ENDIF

    RETURN
  END SUBROUTINE combineBkgSfcIR


  SUBROUTINE overrideBkg(this,dback0,ut_cov_net,isLand)

    !<f90Subroutine>********************************************************
    !
    ! NAME:
    !
    !   overrideBkg
    !
    ! PURPOSE:
    !
    !   Tune background (a priori and error covariance) data, overriding any prior 
    !   settings.
    !
    ! SYNTAX:
    !
    !   CALL overrideBkg(dback0, ut_cov_net, isLand)
    !
    ! ARGUMENTS:
    !
    !   INPUTS:
    !   
    !   isLand      LOGICAL  flag for land surface
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

    CLASS(background_t), intent(in) :: this

    !---Input variables
    REAL, DIMENSION(:),   INTENT(INOUT) :: dback0
    REAL, DIMENSION(:,:), INTENT(INOUT) :: ut_cov_net
    LOGICAL,              INTENT(IN)    :: isLand

    IF (.not. isLand) THEN     !---Ocean
       !Liq
       IF (this%NG%cldLiq .GT. 0) THEN
          dback0(this%IG%cldLiq)    = this%backcldtopoc
          dback0(this%IG%cldLiq+1)  = this%backcldthkoc
          dback0(this%IG%cldLiq+2)  = this%backcldamtoc
          if (this%NG%cldLiq .GT. 3) dback0(this%IG%cldLiq+3) =this%backclddeffoc
          ut_cov_net(this%IG%cldLiq,this%IG%cldLiq)     = this%varcldtopoc
          ut_cov_net(this%IG%cldLiq+1,this%IG%cldLiq+1) = this%varcldthkoc
          ut_cov_net(this%IG%cldLiq+2,this%IG%cldLiq+2) = this%varcldamtoc
          if (this%NG%cldLiq .GT. 3) &
             ut_cov_net(this%IG%cldLiq+3,this%IG%cldLiq+3) =this%varclddeffoc
       ENDIF
       !Ice
       IF (this%NG%cldIce .GT. 0) THEN
          dback0(this%IG%cldIce)   = this%backicetopoc
          dback0(this%IG%cldIce+1) = this%backicethkoc
          dback0(this%IG%cldIce+2) = this%backiceamtoc
          dback0(this%IG%cldIce+3) = this%backicedeffoc
          ut_cov_net(this%IG%cldIce,this%IG%cldIce)     = this%varicetopoc
          ut_cov_net(this%IG%cldIce+1,this%IG%cldIce+1) = this%varicethkoc
          ut_cov_net(this%IG%cldIce+2,this%IG%cldIce+2) = this%variceamtoc
          ut_cov_net(this%IG%cldIce+3,this%IG%cldIce+3) = this%varicedeffoc
       ENDIF
       IF (this%NG%Tskin .GT. 0 .AND. this%varTskinOc .NE. MISSING_REAL) &
            ut_cov_net(this%IG%Tskin,this%IG%Tskin) = this%varTskinOc
    ELSE                                   !---Land
       IF (this%NG%cldLiq .GT. 0) THEN
          dback0(this%IG%cldLiq)   = this%backcldtopld
          dback0(this%IG%cldLiq+1) = this%backcldthkld
          dback0(this%IG%cldLiq+2) = this%backcldamtld
          if (this%NG%cldLiq .GT. 3) dback0(this%IG%cldLiq+3) = this%backclddeffld
          ut_cov_net(this%IG%cldLiq,this%IG%cldLiq)     = this%varcldtopld
          ut_cov_net(this%IG%cldLiq+1,this%IG%cldLiq+1) = this%varcldthkld
          ut_cov_net(this%IG%cldLiq+2,this%IG%cldLiq+2) = this%varcldamtld
          if (this%NG%cldLiq .GT. 3) &
             ut_cov_net(this%IG%cldLiq+3,this%IG%cldLiq+3) = this%varclddeffld
       ENDIF
       IF (this%NG%cldIce .GT. 0) THEN
          dback0(this%IG%cldIce)   = this%backicetopld
          dback0(this%IG%cldIce+1) = this%backicethkld
          dback0(this%IG%cldIce+2) = this%backiceamtld
          dback0(this%IG%cldIce+3) = this%backicedeffld
          ut_cov_net(this%IG%cldIce,this%IG%cldIce)     = this%varicetopld
          ut_cov_net(this%IG%cldIce+1,this%IG%cldIce+1) = this%varicethkld
          ut_cov_net(this%IG%cldIce+2,this%IG%cldIce+2) = this%variceamtld
          ut_cov_net(this%IG%cldIce+3,this%IG%cldIce+3) =this%varicedeffld
       ENDIF
       IF (this%NG%Tskin .GT. 0 .AND. this%varTskinLd .NE. MISSING_REAL) &
            ut_cov_net(this%IG%Tskin,this%IG%Tskin) = this%varTskinLd
    ENDIF

    RETURN
  END SUBROUTINE overrideBkg

  SUBROUTINE transformBkg(this,ut_cov_net,ut_cov,umtx)

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

    CLASS(background_t), intent(inout) :: this

    !---Input variables
    REAL, DIMENSION(:,:), INTENT(IN)    :: ut_cov_net
    REAL, DIMENSION(:,:), INTENT(INOUT) :: ut_cov,umtx

    this%tmpl(1:this%nParR, 1:this%nParG) = MATMUL(umtx(1:this%nParR,1:this%nParG), &
         ut_cov_net(1:this%nParG,1:this%nParG))
    ut_cov(1:this%nParR,1:this%nParR)=MATMUL(this%tmpl(1:this%nParR, 1:this%nParG), &
         TRANSPOSE(umtx(1:this%nParR,1:this%nParG)))

    RETURN
  END SUBROUTINE transformBkg

  SUBROUTINE validateExtAtmos(this,NGExt,nmol,molID,molIDext, &
       xFlgTemp,xFlgMols,xFlgPsfc,xFlgWind,xFlgCldLiq,xFlgCldIce)

    !<f90Subroutine>********************************************************
    !
    ! NAME:
    !
    !   validateExtAtmos
    !
    ! PURPOSE:
    !
    !   Check external data structure against background data structure and 
    !   retrieval switches.
    !
    ! SYNTAX:
    !
    !   CALL validateExtAtmos(NGExt, nmol, molID, molIDext, xFlgTemp, 
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

    CLASS(background_t), intent(inout) :: this

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
    character (len=mxLength) :: procName = '[BackgroundModule::validateExtAtmos]:'

    this%externalFlags = EXT_FLAGS_ALL_OFF
    lenNGExt = getVectorLength(NGExt)
    if (lenNGExt == 0) return

    !---Allocate vector to contain external background data:
    allocate(this%xExt(lenNGExt))
    !---Initialize to tell downstream there are no valid data
    this%xExt = MISSING_REAL

    if (xFlgTemp) then
       if (NGExt%temp .NE. this%NG%temp .OR. NGExt%tskin .NE. this%NG%tskin) then
          print*,trim(procName)//'Error:  ', &
               'TEMP/or/TSKIN external fields not matching BKG data'
          call errorHalt(1)
       endif
       if (this%NG%Temp .GT. 0) this%externalFlags%xTemp=xFlgTemp
       if (this%NG%Tskin .GT. 0) this%externalFlags%xTskin=xFlgTemp
    endif
    if (xFlgPsfc) then
       if (NGExt%psfc .NE. this%NG%psfc) then
          print*,trim(procName)//'Error:  ', &
               'PSFC external fields not matching BKG data'
          call errorHalt(1)
       endif
       if (this%NG%psfc .GT. 0) this%externalFlags%xPsfc=xFlgPsfc
    endif

    this%nMolExt=getNMol(molIDext)
    this%mapExtMols=0
    do j=1,this%nMolExt
       do i=1,nMol
          if (molIDext(j) == molID(i) .AND. xFlgMols(i)) then
             this%mapExtMols(j)=i
             if (NGExt%mol(j) .NE. this%NG%mol(i)) then
                print*,'MolID #',molID(i),':'
                print*,trim(procName)//'Error:  ', &
                     'External fields not matching BKG data'
                call errorHalt(1)
             endif
             exit
          endif
       enddo
       if (NGExt%mol(j) .GT. 0) this%externalFlags%xMols(j)=xFlgMols(i)
    enddo
    if (xFlgCldLiq) then
       if (NGExt%cldLiq .NE. this%NG%cldLiq) then
          print*,trim(procName)//'Error:  ', &
               'CloudLiq external fields not matching BKG data'
          call errorHalt(1)
       endif
       if (this%NG%cldLiq .GT. 0) this%externalFlags%xCldLiq=xFlgCldLiq
       this%externalFlags%xCldLiq=xFlgCldLiq
    endif
    if (xFlgCldIce) then
       if (NGExt%cldIce .NE. this%NG%cldIce) then
          print*,trim(procName)//'Error:  ', &
               'CloudIce external fields not matching BKG data'
          call errorHalt(1)
       endif
       if (this%NG%cldIce .GT. 0) this%externalFlags%xCldIce=xFlgCldIce
    endif
    if (xFlgWind) then
       if (NGExt%wind .NE. this%NG%wind) then
          print*,trim(procName)//'Error:  ', &
               'WIND external fields not matching BKG data'
          call errorHalt(1)
       endif
       if (this%NG%wind .GT. 0) this%externalFlags%xWind=xFlgWind
    endif
    RETURN
  END SUBROUTINE validateExtAtmos


  SUBROUTINE validateExtSfcIR(this,nemirExt,wvnExt)

    !<f90Subroutine>********************************************************
    !
    ! NAME:
    !
    !   validateExtSfcIR
    !
    ! PURPOSE:
    !
    !   Check external emissivity against background parameters.
    !
    ! SYNTAX:
    !
    !   CALL validateExtSfcIR(nemirExt, wvnExt)
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

    CLASS(background_t),intent(inout) :: this

    INTEGER,               INTENT(IN)    :: nemirExt
    REAL,    DIMENSION(:), INTENT(IN)    :: wvnExt

    ! Local variables
    INTEGER                        :: i
    character (len=mxLength) :: procName = '[BackgroundModule::validateExtSfcIR]:'

    ! Constants
    real,    parameter             :: dSpectLim=0.001 ! Spectral match threshold
    ! (fractional) 

    !---Allocate vector to contain external emissivity and reflectivity data:
    allocate(this%emRfExt(2*this%NemIrG))
    !---Initialize to tell downstream there are no valid data
    this%emRfExt = MISSING_REAL

    if (nemirExt == 0) return  ! No external emis data to process

    if (nemirExt /= this%NemIrG) then
       if (allocated(this%cvtIR)) deallocate(this%cvtIR)
       allocate(this%cvtIR(this%NemIrG,nemirExt))
       call spectralConvert(wvnExt,this%FrqEmIR,this%cvtIR,extrapOK=.TRUE.)
    elseif (any((abs(wvnExt-this%FrqEmIR(1:this%NemIrG))/wvnExt) > dSpectLim)) then
       print *,trim(procName)//'Error:  Channel set mismatch;'
       do i=1,this%NemIrG
          write(*,'(i4,2f10.4)')i,wvnExt(i),this%FrqEmIR(i)
       enddo
       call errorHalt(1)    
    endif

    RETURN
  END SUBROUTINE validateExtSfcIR

  subroutine setCloudMode(this,cloudType,retrievalFlags,flagDynCldMask)

    !<f90Subroutine>********************************************************
    !
    ! NAME:
    !
    !   setCloudMode
    !
    ! PURPOSE:
    !
    !   Set controls to disable retrieval of cloud parameters, depending 
    !   on dynamic cloud mask
    !
    ! SYNTAX:
    !
    !   CALL setCloudMode(cloudType, retrievalFlags)
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

    CLASS(background_t) :: this

    ! Set controls to disable retrieval of cloud parameters,
    ! depending on dynamic cloud mask

    integer,                   intent(in)     :: cloudType
    LOGICAL,                   intent(in)     :: flagDynCldMask
    TYPE(retFlags_t),          intent(inout)  :: retrievalFlags
    character (len=mxLength) :: procName = '[BackgroundModule::setCloudMode]:'

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
       print*,trim(procName)//'Error: ', &
            'Invalid flagDynCldMask:',flagDynCldMask
       call errorHalt(1)
    end select

    return

  end subroutine setCloudMode

  SUBROUTINE bkgDestroy(this)
    CLASS(background_t), intent(inout) :: this
    if (allocated(this%Pres)) deallocate(this%Pres)
    if (allocated(this%AtmTypes)) deallocate(this%AtmTypes)
    if (allocated(this%SfcTypesMW)) deallocate(this%SfcTypesMW)
    if (allocated(this%SfcTypesIR)) deallocate(this%SfcTypesIR)
    if (allocated(this%extrvec)) deallocate(this%extrvec)
    if (allocated(this%Freq)) deallocate(this%Freq)
    if (allocated(this%FrqEmIR)) deallocate(this%FrqEmIR)
    if (allocated(this%Pol)) deallocate(this%Pol)
    if (allocated(this%tmpl)) deallocate(this%tmpl)
    if (allocated(this%umtx_clim)) deallocate(this%umtx_clim)
    if (allocated(this%vmtx_clim)) deallocate(this%vmtx_clim)
    if (allocated(this%back_clim_atm_all)) deallocate(this%back_clim_atm_all)
    if (allocated(this%back_clim_sfc_mw)) deallocate(this%back_clim_sfc_mw)
    if (allocated(this%back_clim_sfc_ir)) deallocate(this%back_clim_sfc_ir)
    if (allocated(this%u_clim_atm_all)) deallocate(this%u_clim_atm_all)
    if (allocated(this%v_clim_atm_all)) deallocate(this%v_clim_atm_all)
    if (allocated(this%u_clim_sfc_mw)) deallocate(this%u_clim_sfc_mw)
    if (allocated(this%v_clim_sfc_mw)) deallocate(this%v_clim_sfc_mw)
    if (allocated(this%u_clim_sfc_ir)) deallocate(this%u_clim_sfc_ir)
    if (allocated(this%v_clim_sfc_ir)) deallocate(this%v_clim_sfc_ir)
    if (allocated(this%ut_cov_clim_atm_all)) deallocate(this%ut_cov_clim_atm_all)
    if (allocated(this%ut_cov_clim_sfc_mw)) deallocate(this%ut_cov_clim_sfc_mw)
    if (allocated(this%ut_cov_clim_sfc_ir)) deallocate(this%ut_cov_clim_sfc_ir)
    if (allocated(this%xbakg_clim)) deallocate(this%xbakg_clim)
    if (allocated(this%xbakg_emissdb)) deallocate(this%xbakg_emissdb)
    if (allocated(this%ut_cov_clim)) deallocate(this%ut_cov_clim)
    if (allocated(this%ut_cov_ext)) deallocate(this%ut_cov_ext)
    if (allocated(this%ut_cov_emissdb)) deallocate(this%ut_cov_emissdb)
    if (allocated(this%xExt)) deallocate(this%xExt)
    if (allocated(this%emRfExt)) deallocate(this%emRfExt)
    if (allocated(this%cvtIR)) deallocate(this%cvtIR)

  END SUBROUTINE bkgDestroy

END MODULE BackgroundModule
