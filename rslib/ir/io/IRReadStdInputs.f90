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
!  MODULE IRReadStdInputs.
!
!--------------------------------------------------------------
MODULE IRReadStdInputs

! <f90Module>***********************************************************
!
! NAME:
!
!   IRReadStdInputs
!
! PURPOSE:
!
!   Load the primary IR-related control input parameters. Control all of the
!   IR-related hard-coded file unit numbers. Set all file names. Unit file names
!   and numbers are public.
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

  IMPLICIT NONE
  INTEGER, PRIVATE, PARAMETER :: NCDF_DUMMY=999
  !------------------------------------------------------------------------
  !	Unit Numbers                            File Names
  !------------------------------------------------------------------------
  !---IR-related files
  INTEGER,PARAMETER :: U_osscoefs_ir = 50;  CHARACTER(LEN=200) :: F_osscoefs_ir
  INTEGER,PARAMETER :: U_ossoptdpth_ir=51;  CHARACTER(LEN=200) :: F_ossoptdpth_ir
                                            CHARACTER(LEN=200) :: F_defProfFile

  INTEGER,PARAMETER :: U_chan        = 52;  CHARACTER(LEN=200) :: F_chan
  INTEGER,PARAMETER :: U_nedn        = 53;  CHARACTER(LEN=200) :: F_nedn
  INTEGER,PARAMETER :: U_solarflux   = 54;  CHARACTER(LEN=200) :: F_solarflux
  INTEGER,PARAMETER :: U_ossScatRegr_ir = 55;  CHARACTER(LEN=200) :: F_ossScatRegr_ir
  INTEGER,PARAMETER :: U_namelist    = 71;  CHARACTER(LEN=200) :: F_namelist
  INTEGER :: U_RadIR             = NCDF_DUMMY;  CHARACTER(LEN=200) :: F_RadIR
  INTEGER :: U_regr              = NCDF_DUMMY;  CHARACTER(LEN=200) :: F_regr
  INTEGER :: U_retrRegr          = NCDF_DUMMY;  CHARACTER(LEN=200) :: F_retrRegr
  INTEGER :: U_retrIR            = NCDF_DUMMY;  CHARACTER(LEN=200) :: F_retrIR
  INTEGER :: U_linInvert         = NCDF_DUMMY;  CHARACTER(LEN=200) :: F_linInvert
  INTEGER :: U_IRemissExt        = NCDF_DUMMY;  CHARACTER(LEN=200) :: F_IRemissExt
CONTAINS

  SUBROUTINE IRreadStd(REGRESon,PHYSon,MWon,iStrat,flagDynCldMask, &
     extFlgEmIRin,genLinInvertIn,ossScatRegrIrIn,nlfile)
    !---In/Out variables

!<f90Subroutine>********************************************************
!
! NAME:
!
!   IRreadStd
!
! PURPOSE:
!
!   Read IR input parameters/file names used in simulations and retrieval
!
! SYNTAX:
!
!   CALL IRreadStd(REGRESon, PHYSon, MWon, iStrat,
!      flagDynCldMask, extFlgEmIRin, genLinInvertIn, ossScatRegrIrIn,
!      nlfile)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   nlfile*          CHAR     File path from which to read the namelist
!
!   INPUTS/OUTPUTS:
!
!   REGRESon         LOGICAL  Flag to include regression mode in IR retrieval
!   PHYSon           LOGICAL  Flag to include physical mode in IR retrieval
!   MWon             LOGICAL  Flag to include microwave
!   iStrat           INTEGER  Stratification level flag
!   flagDynCldMask   LOGICAL  Use dynamic cloud mask?
!   extFlgEmIRin*    LOGICAL  Flag for using external IR emissivity data
!   genLinInvertIn*  LOGICAL  Flag to generate output for linear inversion
!   ossScatRegrIrIn* LOGICAL  Flag for using predictors in scattering
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    LOGICAL,          INTENT(inout)           :: REGRESon,PHYSon,MWon
    INTEGER,          INTENT(inout)           :: iStrat
    LOGICAL,          INTENT(inout)           :: flagDynCldMask
    LOGICAL,          INTENT(inout), OPTIONAL :: extFlgEmIRin
    LOGICAL,          INTENT(INOUT), OPTIONAL :: genLinInvertIn
    LOGICAL,          INTENT(INOUT), OPTIONAL :: ossScatRegrIrIn
    character(len=*), intent(in),    optional :: nlfile
    !---Local variables
    LOGICAL               :: ossScatRegrIRon=.FALSE.
    LOGICAL               :: extFlgEmIR=.FALSE.
    LOGICAL               :: genLinInvert=.FALSE.
    CHARACTER(LEN=200)    :: osscoefs_ir,ossoptdpth_ir
    CHARACTER(LEN=200)    :: osspred_ir
    character(len=200)    :: chanFileName,nednFile
    CHARACTER(LEN=200)    :: RadIRFile,solarfluxFile,regrFile
    CHARACTER(LEN=200)    :: RetrRegrFile,RetrIrFile
    CHARACTER(LEN=200)    :: linInvertFile,IRemissExtFile
    CHARACTER(LEN=200)    :: defProfFile
    NAMELIST /IRInput/ REGRESon,PHYSon,MWon,iStrat,ossScatRegrIRon, &
         flagDynCldMask,RadIRFile, &
         osscoefs_ir,ossoptdpth_ir,osspred_ir,chanFileName,nednFile,&
         solarfluxFile, defProfFile,regrFile,RetrRegrFile,RetrIRFile,&
         genLinInvert,linInvertFile,extFlgEmIR,IRemissExtFile

    if (present(nlfile)) then
       open(U_namelist,file=F_namelist)
       read(U_namelist,IRInput)
       close(U_namelist)
    else
    READ(*,IRInput)
    endif

    if (present(extFlgEmIRin)) extFlgEmIRin=extFlgEmIR
    if (present(genLinInvertIn)) genLinInvertIn=genLinInvert
    if (present(ossScatRegrIrIn)) ossScatRegrIrIn=ossScatRegrIRon

    F_RadIR         = RadIRFile
    F_osscoefs_ir   = osscoefs_ir
    F_ossoptdpth_ir = ossoptdpth_ir
    F_ossscatregr_ir= osspred_ir
    F_defProfFile   = defProfFile
    F_chan          = chanFileName
    F_nedn          = nednFile
    F_solarflux     = solarfluxFile
    F_regr          = regrFile
    F_RetrRegr      = RetrRegrFile
    F_RetrIR        = RetrIRFile
    F_linInvert     = linInvertFile
    F_IRemissExt    = IRemissExtFile
  END SUBROUTINE IRreadStd

  !Returns the global path(s) upon request
  SUBROUTINE IRgetGlobalPath(P_osscoefs_ir,P_ossoptdpth_ir,P_chan, &
       P_nedn,P_solarflux,P_ossScatRegr_ir,P_namelist,P_RadIR, &
       P_regr,P_retrRegr,P_retrIR,P_linInvert,P_IRemissExt)

    CHARACTER(LEN=*), intent(inout), optional :: P_osscoefs_ir
    CHARACTER(LEN=*), intent(inout), optional :: P_ossoptdpth_ir
    CHARACTER(LEN=*), intent(inout), optional :: P_chan
    CHARACTER(LEN=*), intent(inout), optional :: P_nedn
    CHARACTER(LEN=*), intent(inout), optional :: P_solarflux
    CHARACTER(LEN=*), intent(inout), optional :: P_ossScatRegr_ir
    CHARACTER(LEN=*), intent(inout), optional :: P_namelist
    CHARACTER(LEN=*), intent(inout), optional :: P_RadIR
    CHARACTER(LEN=*), intent(inout), optional :: P_regr
    CHARACTER(LEN=*), intent(inout), optional :: P_retrRegr
    CHARACTER(LEN=*), intent(inout), optional :: P_retrIR
    CHARACTER(LEN=*), intent(inout), optional :: P_linInvert
    CHARACTER(LEN=*), intent(inout), optional :: P_IRemissExt

    if (present(P_osscoefs_ir)) P_osscoefs_ir = F_osscoefs_ir
    if (present(P_ossoptdpth_ir)) P_ossoptdpth_ir = F_ossoptdpth_ir
    if (present(P_chan)) P_chan = F_chan
    if (present(P_nedn)) P_nedn = F_nedn
    if (present(P_solarflux)) P_solarflux = F_solarflux
    if (present(P_ossScatRegr_ir)) P_ossScatRegr_ir = F_ossScatRegr_ir
    if (present(P_namelist)) P_namelist = F_namelist
    if (present(P_RadIR)) P_RadIR = F_RadIR
    if (present(P_regr)) P_regr = F_regr
    if (present(P_retrRegr)) P_retrRegr = F_retrRegr
    if (present(P_retrIR)) P_retrIR = F_retrIR
    if (present(P_linInvert)) P_linInvert = F_linInvert
    if (present(P_IRemissExt)) P_IRemissExt = F_IRemissExt

  END SUBROUTINE IRgetGlobalPath

END MODULE IRReadStdInputs



