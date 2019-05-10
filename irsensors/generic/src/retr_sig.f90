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

program retr

! <f90Main>*************************************************************
!
! NAME:
!
!   retr
!
! PURPOSE:
!
!   Main retrieval program.
!
! INCLUDES:
!
!   None
!
!*************************************************************</f90Main>

  USE StateIndexModule
  USE LvlInterp
  USE ncdf_module, only: writeNcdfData, writeNcdfAttr,nf_open,nf_nowrite,nf_noerr,closencdffile
  USE scene_io_module, ONLY: openScene, getScene, putScene, closeScene, &
       TB_UNITS, RADMW_UNITS
  USE rad_io_module
  USE RadFile
  USE IRBkgPrepModule, ONLY:retFlags_t,RET_FLAGS_ALL_ON,extFlags_t,&
                       IRbkgPrepInit,IRsetFOVbkg,IRcombBkg,IRtuneBkg,&
                       IRsetRetrVec,IRtransBkg,validateExt,validateEmis,&
                       IRsetExtBkg,IRsetEmisBkg,IRcheckTransformConform,&
                       cloudModeTest
  USE RegrPrepModule
  USE MapInvert, ONLY: map_geo2retr,set_MW_Invert,MOL_TRAN_LOG
  USE IRMapInvert
  USE ChkValid  
  USE oss_mw_module
!  USE oss_ir_module_scat
  USE oss_ir_module
  USE NoiseModule
  USE IRNoiseModule
  USE IRChannels
  USE SurfaceTypingModule
  USE ReadStdInputs
  USE IRReadStdInputs
  USE InvertModule
  use FiniteDiffModule, Only: &
       initFiniteDiff, &
       loadFiniteDiff, &
       finiteDiff, &
       destroyFiniteDiff
  use constants, ONLY : MISSING_REAL,MISSING_CHAR
  USE VertCoord, ONLY: mxCoordTyp,Pcoord
  USE LinInvert, ONLY: initBuildLin,buildLinearCoef,destroyLinInvert
  implicit none
  
  !-----------------------------------------------------
  !                DECLARATION SECTION
  !-----------------------------------------------------
  INTEGER, PARAMETER                     :: mxiter=15,mxchan=2600,mxlev=101 
  !----Geophy and Retr Space Arrays
  real,    dimension(mxiter)             :: rms_mw,chisq_mw,rms_ir,chisq_ir
  real,    dimension(:), allocatable     :: Y,Ym,Sy,delY,xkEmMw,EmMw
  real,    dimension(:), allocatable     :: Ylast,xGesLast
  real,    dimension(:), allocatable     :: tbEquiv
  real,    dimension(:), allocatable     :: alphaCtr,betaCtr
  real,    dimension(:), allocatable     :: dradChan
  real,    dimension(:), allocatable     :: xBakG,xGesMWSav,xGesG
  real,    dimension(:), allocatable     :: xGes,xGesNew
  real,    dimension(:,:), allocatable   :: xktG,xkt
  real,    dimension(:), allocatable     :: PrefOSS ! From OSS file, has to be dimensioned here.
  real,    dimension(:), allocatable     :: pref
  integer, dimension(mxchan)             :: kchan
  integer, dimension(:), allocatable     :: molID,molTran
  character(len=12), dimension(maxMol)   :: molOnRT

  !---Radiance file data for all locations
  real, dimension(:,:), allocatable :: radAll,radTmp
  real, dimension(:  ), allocatable :: eiaAll,latAll,lonAll
  integer(selected_int_kind(1)), dimension(:), allocatable :: cldTypAll
  
  !----Covariance/Backg Arrays
  real,    dimension(:,:), allocatable   :: clw_cov,ciw_cov,ut_cov_net,Sx
  real,    dimension(:,:), allocatable   :: umtx,vmtx
  real(kind=DBL_INVT), dimension(:,:), allocatable   :: dSx
  real(kind=DBL_INVT)                                :: detmnt
  
  !----IR-related arrays
  REAL,    DIMENSION(:), allocatable     :: frqEmIR
  REAL,    DIMENSION(:), allocatable     :: EmRf
  REAL,    DIMENSION(:,:,:), allocatable :: xkEmRf
  REAL,    DIMENSION(mxchan)             :: wvnOSS ! From OSS file, has to be max-dimensioned here.
  real, dimension(:), allocatable        :: wvnRad,Airs_Err
  INTEGER, PARAMETER                     :: lenIDir=12
  CHARACTER(LEN=lenIDir),DIMENSION(:),ALLOCATABLE :: chanIDtmp,chanIDir
  integer, dimension(:), allocatable     :: chanNum
  !-- Noise related variables:
  real, dimension(:), allocatable        :: DevNoise,AtmNoise
  real                                   :: sumMwDevNoise,sumIRDevNoise
      
  !----MW Instrument Data (from OSS file, have to be max-dimensioned here).
  INTEGER, DIMENSION(mxchan)             :: polOSS
  REAL,    DIMENSION(mxchan)             :: FrqOSS

  !----MW Instrument Data (from SDR file)
  INTEGER, DIMENSION(:), allocatable     :: pol_rad
  REAL,    DIMENSION(:), allocatable     :: frq_rad

  !----MW Instrument Data (from cov file)
  integer, dimension(:),   allocatable   :: pol
  real, dimension(:),      allocatable   :: frq

  !----Regression
  INTEGER, PARAMETER                     :: nmolRegr=1
  INTEGER, PARAMETER                     :: mxRdEof=200
  INTEGER, DIMENSION(nmolRegr)           :: molIDRegr=(/idH2O/)
  INTEGER                                :: nParRegr,nparGRegr
  !nParMax is nPar when nothing is turned OFF but eig truncation applied.
  INTEGER                                :: nParMax
  INTEGER                                :: nlevRegr,iclassatmRegr 
  INTEGER                                :: nchanRegr,neofRegr
  REAL                                   :: tsfcRegr
  REAL, DIMENSION(:), pointer            :: prefRegr
  REAL, DIMENSION(:), allocatable        :: xbackG_regr,xregr
  REAL, DIMENSION(mxchan)                :: rback_regr
  REAL, DIMENSION(mxchan)                :: wght_regr
  REAL, DIMENSION(:,:), allocatable      :: regr_coef,umtx_regr
  REAL, DIMENSION(mxchan,MxRdEof)        :: umtx_rad
  TYPE(StateIndex_t)                     :: IGRegr,NGRegr

  !---External Data
  Logical,Dimension(maxMol):: extFlgMol
  Logical                  :: extFlgTemp,extFlgPsfc
  Logical                  :: extFlgCldLiq,extFlgCldIce,extFlgWind
  integer                  :: nforExt,nparGExt
  integer                  :: nforEmis,nemirExt=0
  integer,dimension(maxMol):: molIDext
  real, dimension(:,:), allocatable      :: xbakG_Ext
  real, dimension(:,:), allocatable      :: xbakG_Emis
  real, dimension(:),   allocatable      :: wvnExt
  TYPE(StateIndex_t)       :: IGExt,NGExt
  Logical                  :: extFlgEmIR

  !---Single variables
  character(len=8)   :: IdProf
  logical            :: debug,REGRESon,PHYSon,MWon,DynamicNoise
  LOGICAL            :: flagDynCldMask
  LOGICAL            :: use2ndOrderPlanck
  LOGICAL            :: genLinInvert
  logical            :: isLand
  integer            :: nfor,itermw,iter,iter_tot,iemmwG,IEmIrG
  integer            :: nmol
  INTEGER            :: nemIrG ! number of IR emissivities in Geospace
  INTEGER            :: nemIRsel
  integer            :: iexit
  integer            :: ifail,nch,nchanOSS,nlev,nlevOSS
  integer            :: ih2o,i,ii,ifor,time(6)
  integer            :: ich
  integer            :: nprofs,mwCld,iCell,icasc,iStrat,iembkgflg,iextbkg
  integer            :: nchanmw,nchanmwRad,nchanir,nchanirRad,nchan
  integer            :: igeo,iclassatm
  INTEGER            :: nPar,nparG,npargAtm
  INTEGER            :: iRadOrTb
  integer            :: roUnits
  integer            :: nParGon
  integer            :: cloudType
  real               :: deltaP=0.05,deltaFrq=0.001
  real               :: tsfc
  real               :: lat,lon
  real               :: rmsMw,rms,chiSqMw,chisq0,chisq,chisqratio
  real               :: xSfc
  character(len=7)   :: cldTyp='layer  '
  integer,parameter  :: nHyd=2
  LOGICAL            :: enableScatter,ossScatRegrIRon
  LOGICAL            :: xOutRangeFlag
  CHARACTER(LEN=mxCoordTyp)              :: vCoordTyp
  character(len=20),dimension(nHyd)      :: &
       dataType_cld=(/'CldLiq ','CldIce '/)
  TYPE(StateIndex_t) :: IG, NG, IR, NR
  TYPE(StateIndex_t) :: IGRout,NGRout
  integer                                :: nAng
  REAL, DIMENSION(:), ALLOCATABLE        :: press
  real, dimension(:), allocatable        :: eia
  real, dimension(:), allocatable        :: psfc
  real, dimension(:), allocatable        :: plandAvg

  real, dimension(:),     allocatable    :: xOffG
  real, dimension(:,:),   allocatable    :: xGainG
  real, dimension(:,:),   allocatable    :: xOffGall
  real, dimension(:,:,:), allocatable    :: xGainGall
  real, dimension(:),     allocatable    :: chiSqAll
  real, dimension(:),     allocatable    :: TsfcAll
  real                                   :: TsfcLast

  !---retrTuning Namelist items
  integer            :: ngas(maxMol)
  integer            :: nTskin,ntemp,nemmw,nemir,nCldLiq,nCldIce,nxiterMw,nxiter
  logical            :: logSize
  real               :: drad,alpha1slope,alpha1int
  real               :: chisqconvMw,chisqconv,chisqRatioconv
  real               :: cldamt1stguess,iceamt1stguess
  real               :: TightCLDcov
  REAL               :: plandMaxOc

  !--- items to save original values prior to calling IRsetRetrVec
  INTEGER            :: nCldLiqOrig, nCldIceOrig, nEmMw_orig, nEmIR_orig
  TYPE(StateIndex_t) :: IR_orig, NR_orig

  !--------------------------------------------------------------
  ! global variable to set On/Off all the retrieval flags in a
  ! call to IRsetRetrVec(). For now, it will be set to all ON.
  !--------------------------------------------------------------
  TYPE(retFlags_t)   :: retrievalFlags = RET_FLAGS_ALL_ON

  ! loadDone implementation provides the option to call load*() in advance
  logical            :: loadAlgCfgDone=.false.
  logical            :: loadAuxDone=.false.
  logical            :: loadExtDone=.false.
  logical            :: loadEmisDone=.false.

  !-----------------------------------------------------
  !    READ DATA THROUGH NAMELISTS/LUTs/Attributes
  !-----------------------------------------------------
  call readStd(nprofs,debug,mwCld,iCell,icasc, &
       iembkgflgIn=iembkgflg,iextbkgIn=iextbkg, &
       kchanOut=kchan,DynamicNoiseIn=DynamicNoise, &
       enableScatterIn=enableScatter,use2ndOrderPlanckIn=use2ndOrderPlanck, &
       extFlgTempOut=extFlgTemp,extFlgMolOut=extFlgMol, &
       extFlgPsfcOut=extFlgPsfc,extFlgCldLiqOut=extFlgCldLiq, &
       extFlgCldIceOut=extFlgCldIce,extFlgWindOut=extFlgWind, &
       molOnRTin=molOnRT)
  call IRreadStd(REGRESon,PHYSon,MWon,iStrat,flagDynCldMask, &
       extFlgEmIRin=extFlgEmIR,genLinInvertIn=genLinInvert, &
       ossScatRegrIrIn=ossScatRegrIRon)

IF (enableScatter) THEN
print *
print *,'SCATTERING IS DISABLED UNTIL oss_ir_module_scat SUPPORTS VARIABLE GRID'
print *
ENDIF

  !---Load retr tuning parameters
  call loadAlgconfig(F_algconfig)
  IF (nxiterMw > mxiter)then
     print*,' Error: NxiterMw > mxiter ',nxiterMw,mxiter
     call errorHalt(1)
  endif
  IF (nxiter > mxiter)then
     print*,' Error: Nxiter > mxiter ',nxiter,mxiter
     call errorHalt(1)
  endif
  
  !----Derive nmol and molID from inputs
  nmol=COUNT(molOnRT /= MISSING_CHAR)
  if (nmol <= 0) then
     print *,'err:[retr] No molecules selected for RT'
     call errorHalt(1)
  endif
  allocate (molID(nmol),molTran(nmol))
  call charToMolID(molOnRT(1:nmol),molID)
  ih2o=whereH2O(molID)
  
  if(REGRESon) then
     call getRegrParam(nparGRegr,nparRegr,nchanregr,neofregr,&
          IGRegr,NGRegr,nlevRegr,prefRegr,MolIdRegr) 
     allocate(xbackG_regr(NParGRegr),xregr(NParGRegr),umtx_regr(NParGRegr,NParRegr))
     allocate(regr_coef(MxRdEof,nparRegr)) 
  endif
  CALL IRbkgPrepInit(nparG,IG,NG,nlev,vCoordTyp,pref,nchanmw,Frq,pol, &
       nemIrG,frqEmIR,iemmwG,IEmIrG,MolId,molTran, &
       limitCfg=F_algconfig,bkgCfg=F_covar)  

  ALLOCATE (press(nlev))

  NParGATm=getVectorLength(NG)

  
  !----Set retrieval-space pointers
  call setXPtrR(nTskin, Ntemp, NGas, nCldLiq, nCldIce, molID, nmol, NG, IR, NR)

  if (nemir == 0) then
     nemIRsel  = 0
  else
     nemIRsel  = nEmIRG  ! if emIR is retrieved, must use all hinge points
  endif

  !----Verify that the selected number of variables to retrieve are
  !----compatible with the transformation matrices
  call IRcheckTransformConform(NR,nemmw,nemIRsel)

  !--- Save the original retrieval-space pointers
  !----------------------------------------------------------
  !  Save original versions of variables that are dynamically
  !  changed in the ProfLoop via call to IRsetRetrVec()
  !  to fall back if needed.
  !----------------------------------------------------------
  ncldLiqOrig = nCldLiq 
  ncldIceOrig = nCldIce
  nEmMw_orig  = nemmw
  nEmIR_orig  = nemIRsel
  IR_orig     = IR
  NR_orig     = NR
  nParMax = getVectorLength(NR) + nEmMw_orig + 2 * nEmIR_orig
  allocate(clw_cov(nCldLiqOrig, nCldLiqOrig),ciw_cov(nCldIceOrig, nCldIceOrig))
  !--------------------------------------------
  !         Initialize OSS models for MW and IR
  !         and Perform Consistency Checks
  !--------------------------------------------
  if (MWon) then
     IF (use2ndOrderPlanck) call set2ndOrderPlanck()  ! Instead of default
     call oss_init_mw(U_osscoefs,U_ossoptdpth,F_osscoefs,F_ossoptdpth, &
          molID(1:ih2o),&
          npargAtm,IG%temp,IG%tskin,IG%psfc,IG%mol(1:ih2o), &
          IG%cldliq,nchanOSS,FrqOSS,nlevOSS,PrefOSS,polOSS,iRadOrTb)
     allocate(alphaCtr(nchanmw),betaCtr(nchanmw))
     CALL getPlanckLinCtr(alphaCtr,betaCtr)
     IF (nchanOSS  /= nchanmw)THEN
        print*,' Nchan inconsistent in MW (Oss/Bkg) ',nchanOSS,nchanmw
        call errorHalt(1)
     ENDIF
     IF (ANY(abs(frq(1:nchanmw)-frqOSS(1:nchanmw)) > deltaFrq ))THEN
        print*,' MW Frequencies inconsistent (Bkg/OSS)'
        call errorHalt(1)
     ENDIF
     IF (ANY(pol(1:nchanmw)  /= polOSS(1:nchanmw) ))THEN
        print*,' Polarization inconsistent (Bkg/OSS)'
        call errorHalt(1)
     ENDIF
     IF (TRIM(vCoordTyp) == Pcoord) THEN
        IF (nlevOSS  /= nlev)THEN
           print*,' Nlev inconsistent in MW (Oss/Bkg) ',nlevOSS,nlev
           call errorHalt(1)
        ENDIF
        IF (ANY(abs(pref(1:nlev)-prefOSS(1:nlev)) > deltaP ))THEN
           print*,' Pressure Grid inconsistent in MW (Bkg/OSS)'
           call errorHalt(1)
        ENDIF
     ENDIF
  else
     nchanmw=0
  endif
  
  !----Initialize IR OSS module and get back Info about LUTs
!  call oss_init_ir(U_osscoefs_ir,U_ossoptdpth_ir,U_ossScatRegr_ir, &
!       F_osscoefs_ir,F_ossoptdpth_ir,U_solarflux,F_solarflux,F_cloudTab, &
!       F_ossScatRegr_ir,nmol,molid,npargAtm,IG%temp,IG%tskin,IG%psfc,&
!       IG%mol(1:nmol),IG%cldLiq,IG%cldIce,nemIrG,frqEmIR(1:nemIrG), &
!       ossScatRegrIRon,enableScatter,dataType_cld,cldTyp,nHyd, &
!       nchanOSS,wvnOSS,nlevOSS,PrefOSS)  
  call oss_init_ir(U_osscoefs_ir,U_ossoptdpth_ir, &
       F_osscoefs_ir,F_ossoptdpth_ir,U_solarflux,F_solarflux, &
       nmol,molid,npargAtm,IG%temp,IG%tskin,IG%psfc,&
       IG%mol(1:nmol),nemIrG,frqEmIR(1:nemIrG), &
       nchanOSS,wvnOSS,nlevOSS,PrefOSS)  
  nchanir=nchanOSS
  IF (TRIM(vCoordTyp) == Pcoord) THEN
     IF (nlevOSS /= nlev)THEN
        print*,' Nlev inconsistent in IR (Oss/Bkg) ',nlevOSS,nlev
        call errorHalt(1)
     ENDIF
     IF (ANY(abs(pref(1:nlev)-prefOSS(1:nlev)) > deltaP ))THEN
        print*,' Pressure Grid inconsistent in IR (Bkg/OSS)'
        call errorHalt(1)
     ENDIF
  ENDIF

  ! Initialize linear coefficient building
  IF (genLinInvert) call initBuildLin(IG,NG,nParG,NR_orig,IGRout,NGRout,nParGon)

  !--------------------------------------
  !    OPEN SDR files and read in headers
  !--------------------------------------
  if (MWon) then
     call queryRad(fname=F_Rad,nchan=nchanmwrad)

     IF (nchanmwRad /= nchanmw) then
        print*,' Nchanmw inconsistent (Rad/Bkg) ',nchanmwRad,nchanmw
        call errorHalt(1)
     endif

     allocate(pol_rad(nchanmw),frq_rad(nchanmw))
     call openRad(ncid=U_Rad,fname=F_Rad,pol=pol_rad, &
          frq=frq_rad,nchan=nchanmwrad,nrec=NFOR,status='old')

     IF (ANY(abs(frq(1:nchanmw)-frq_rad(1:nchanmw)) > deltaFrq ))then
        print*,' MW Frequencies inconsistent (Bkg/Rad)'
        call errorHalt(1)
     endif
     IF (ANY(pol(1:nchanmw)  /= pol_rad(1:nchanmw) ))then     
        print*,' Polarizations inconsistent (Bkg/Rad)'
        call errorHalt(1)
     endif

     allocate(xkEmMw(nchanmw),EmMw(nchanmw))

  endif
  
  !---Get radiance file dimensions and allocate arrays:
  call getDimRad(U_RadIR,F_RadIR,nChanIRrad,nfor,nAng)
  allocate(eia(nAng))
  IF(nprofs < 1 .or. nprofs > nfor) nprofs=nfor
  allocate( wvnRad(nChanIRrad) )
  allocate( radAll(nChanIR,nfor) )
  allocate( radTmp(nChanIRrad,nfor) )
  allocate(eiaAll(nfor),latAll(nfor),lonAll(nfor))
  allocate(chanIDir(nChanIR),chanIDtmp(nChanIRrad))
  !---Load radiance data:
  call getRadData(U_RadIR,F_RadIR,radTmp,eiaAll,latAll,lonAll,wvnRad, &
     chanID=chanIDtmp)
  IF (nChanIRrad /= nChanIR)then
     call adjustChanSet(radTmp,wvnRad,wvnOSS(1:nChanIR), &
        nchanIRrad,nChanIR,kchan,radAll,chanIDtmp,chanIDir)
  else
     radAll = radTmp
     chanIDir = chanIDtmp
  endif
  deallocate(radTmp,chanIDtmp)
  if (flagDynCldMask) then
     allocate(cldTypAll(nfor))
     call getCloudData(U_RadIR,F_RadIR,cldTypAll)
  endif

  call loadAux(F_Aux)

  if (ANY((/extFlgTemp,extFlgMol(1:nmol),extFlgPsfc,extFlgWind, &
       extFlgCldLiq,extFlgCldIce/))) then
     call loadExt(F_bkgrdExtern)
     if (nforExt < nprofs) then
        print *,'err:[retr] Not enough external data;  nforExt,nprofs:', &
           nforExt,nprofs
        call errorHalt(1)
     endif
  endif
  call validateExt(NGExt,nmol,molID,molIDext,&
       extFlgTemp,extFlgMol(1:nmol),extFlgPsfc,extFlgWind, &
       extFlgCldLiq,extFlgCldIce)
  if (extFlgEmIR) then
    call loadEmis(F_IRemissExt)
    if (nforEmis < nprofs) then
       print *,'err:[retr] Not enough external emissivity data; '// &
          'nforEmis,nprofs:',nforEmis,nprofs
       call errorHalt(1)
    endif
  endif
  call validateEmis(nemirExt,wvnExt)
  if (enableScatter .AND. ANY( &
       (/NR_orig%temp,NR_orig%mol(1:nmol), &
       nEmMw_orig,nEmIR_orig/) /= 0)) then
     print*,'err:[retr] Jacobians not available for all retrieved '// &
         'variables in scattering mode; ' 
     print*,'See config file entries for temp, Tskin, gas, emis'
     call errorHalt(1)
  endif

  nchan=nchanmw+nchanir
  if (enableScatter) then
     call loadFiniteDiff(F_algconfig)
     call initFiniteDiff(IG,NG,vCoordTyp,nchan)
  endif

  !--- Get IR, NR, nParMax etc when all flags are ON

  allocate(DevNoise(nchan),AtmNoise(nchan),Airs_Err(nchan))
  allocate(Y(nchan),Ym(nchan),Sy(nchan),delY(nchan))
  allocate(Ylast(nchan),xGesLast(NParMax))

  IF (genLinInvert) THEN
    allocate (xOffG(nParGon),xGainG(nParGon,nChan))
    allocate (xOffGall(nParGon,nprofs),xGainGall(nParGon,nChan,nprofs))
    allocate (chiSqAll(nprofs),TsfcAll(nprofs))
  ENDIF

  allocate(tbEquiv(nchan))
  allocate(EmRf(2*nemIrG),xkEmRf(2,nemIrG,nchanir)) 
  
  allocate(xBakG(NParG),xGesMWSav(NParG),xGesG(NParG))
  allocate(xGes(NParMax),xGesNew(NParMax))
  allocate(xktG(NParG,nchan),xkt(NParMax,nchan))
  allocate(ut_cov_net(NParG,NParG),Sx(nParMax,nParMax))
  allocate(dSx(nParMax,nParMax))
  allocate(umtx(nParMax,NParG),vmtx(NParG,nParMax))
  !---------------------------------------
  !             IR channel selection flags
  !---------------------------------------
  call initIRChannels(nchanmw,nchan,wvnOSS,kchan)
  !--------------------------------------------------------
  ! Allocate and initialize dradChan with drad
  !--------------------------------------------------------
  IF (MWon) THEN
     allocate(dradChan(nchanmw))
     ! In radiance model, drad must be converted from dTB to dradiance
     IF (iRadOrTb == 0) THEN
        dradChan=drad
        roUnits=TB_UNITS
     ELSE
        dradChan(1:nchanmw)=drad/alphaCtr(1:nchanmw)
        roUnits=RADMW_UNITS
     ENDIF
  ENDIF

  !-----------------------------------------------------
  !    OPEN Output Files
  !-----------------------------------------------------
  if (MWon) then
     call openScene(ncid=U_retr,file=F_Retr,MolId=molId, &
          polarity=pol(1:nchanmw),freq=frq(1:nchanmw),casename=F_Retr, &
          nLevel=nlev,nchmw=nchanmw,vCoordTyp=vCoordTyp,status='new', &
          IG=IG,NG=NG)
     IF (TRIM(vCoordTyp) == Pcoord) CALL writeNcdfAttr(U_retr, &
        attrName='standardPressureGrid',attr=pref)
  endif
  if (PHYSon) then
     call openScene(ncid=U_retrIR,file=F_RetrIR,MolId=molId, &
          casename=F_RetrIR, &
          sfcGrid=frqEmIR(1:nemIrG),nSfcGrid=nemIrG, &
          nLevel=nlev,nchmw=nchanmw,vCoordTyp=vCoordTyp,status='new', &
          IG=IG,NG=NG)
     IF (TRIM(vCoordTyp) == Pcoord) CALL writeNcdfAttr(U_retrIR, &
        attrName='standardPressureGrid',attr=pref)
     IF (MWon) CALL writeNcdfAttr(U_retrIR, &
        attrName='mwfrequencies',attr=frq(1:nchanmw))
     call openScene(ncid=U_guess,file=F_Guess,MolId=molId, &
          casename=F_Guess, &
          sfcGrid=frqEmIR(1:nemIrG),nSfcGrid=nemIrG, &
          nLevel=nlev,nchmw=nchanmw,IG=IG,NG=NG,vCoordTyp=vCoordTyp, &
          status='new')
     IF (TRIM(vCoordTyp) == Pcoord) CALL writeNcdfAttr(U_guess, &
        attrName='standardPressureGrid',attr=pref)
     IF (MWon) CALL writeNcdfAttr(U_guess, &
        attrName='mwfrequencies',attr=frq(1:nchanmw))
  endif
  if(REGRESon) then
     IF (TRIM(vCoordTyp) /= Pcoord) THEN
        print*,'err:[retr] Regression not supported for non-standard pressures'
        call errorHalt(1)
     ENDIF
     call openScene(ncid=U_retrRegr,file=F_RetrRegr,MolId=molIdRegr, &
          pressure=pref,freq=frq(1:nchan),casename=F_RetrRegr, &
          nLevel=nlevRegr,nchmw=nchan,status='new',IG=IGRegr,NG=NGRegr)
  endif

  if (.not. allocated(chanNum)) allocate(chanNum(nchanir))
  chanNum = (/(ich, ich=1, nchanir)/)
  !-----------------------------------------------------
  !    Main Loop over the profiles
  !-----------------------------------------------------
  ProfLoop: do ifor=1,nprofs
     WRITE(*,'(a13,i5,a9)') '-----Profile#',ifor,'---------'

     !---get SDR data
     
     if (MWon) then
        call getRad(ncid=U_Rad,irec=ifor,xid=Idprof,lat=lat,lon=lon,time=time, &
             eia=eia)
        IF (iRadOrTb == 0) THEN
           CALL getRad(ncid=U_Rad,irec=ifor,tb=Ym(1:nchanmw))
        ELSE
           CALL getRad(ncid=U_Rad,irec=ifor,rad=Ym(1:nchanmw))
        ENDIF
        
        IF (iRadOrTb == 0) THEN
           tbEquiv(1:nchanmw)=Ym(1:nchanmw)
        ELSE               ! Convert from radiance to radiance temperature
           tbEquiv(1:nchanmw)= &
              Ym(1:nchanmw)*alphaCtr(1:nchanmw)+betaCtr(1:nchanmw)
        ENDIF
     endif
     !---Access pre-loaded data for current profile:
     Ym(nChanMW+1:nChan) = radAll(:,ifor)
     eia(1)              = eiaAll(ifor)
     lat                 = latAll(ifor)
     lon                 = lonAll(ifor)
     if (flagDynCldMask) cloudType = INT(cldTypAll(ifor))
     tbEquiv(nchanmw+1:nchan)=Ym(nchanmw+1:nchan)  !for regr
     
     print *,'Profile No. and Scan ID=',ifor,Idprof
     print *,'EIA =', eia(1)
     print *,'lat/lon = ',lat,lon
     
     !---get noise amplitude
     
     if (MWon) then
        CALL deviceNoise(icellres=iCell,tb=tbEquiv(1:nchanmw), &
             rerr=DevNoise(1:nchanmw), &
             sum_rerr=sumMwDevNoise,kchan=kchan(1:nchanmw),nchan=nchanmw, &
             DynamicNoise=DynamicNoise)
        DevNoise(1:nchanmw)=sqrt(DevNoise(1:nchanmw)**2+(0.05)**2)
        IF (iRadOrTb == 1) &
             DevNoise(1:nchanmw)=DevNoise(1:nchanmw)/alphaCtr(1:nchanmw)
        AtmNoise(1:nchanmw)=0.
     endif
     
     call IRdeviceNoise(wvnOSS,Ym(nchanmw+1:nchan), &
          DevNoise(nchanmw+1:nchan),AtmNoise(nchanmw+1:nchan),sumIRDevNoise, &
          kchan(nchanmw+1:nchan),nchanir,chanIDir,chanNum)
     
     !---atmosphere/surface classification
     iexit=0

     call setSfc(igeo,plandAvg(ifor),plandMaxOc)

!cccc
igeo=2
print *,'FORCING igeo'
     print "(' iStrat = ',i2)",iStrat
     print "(' igeo = ',i2)",igeo
     print "(' psurface = ',f8.3)",psfc(ifor)
     print "(' pland = ',f8.5)",plandAvg(ifor)
     
100  continue

     
     ! Regression Retrieval
     if(REGRESon)then
        call classatmRegr(tbEquiv,lat,lon,time,eia(1),iclassatmRegr)
        call setFOVRegr(iclassatmRegr,psfc(ifor),ifor,&
             xbackG_regr,rback_regr,wght_regr,&
             umtx_regr,umtx_rad,regr_coef)
        call invrtRegr(tbEquiv(1:nchanregr),&
             xbackG_regr(1:npargregr),&
             rback_regr(1:npargregr),&
             wght_regr(1:nchanregr),&
             umtx_regr(1:npargregr,1:nparregr),&
             umtx_rad(1:nchanregr,1:neofregr),&
             regr_coef(1:neofregr,1:nparregr),&
             nchanregr,npargregr,&
             neofregr,nparregr,xregr)
        !--Some explicit assignment and special treatment
        tsfc=tsfcregr
        xregr(IGRegr%psfc)=psfc(ifor)
        if (molTran(ih2o) /= MOL_TRAN_LOG) then
           print *,'err:[retr] Regression implementation expects log '//&
              'transform of H2O but background is not log'
           call errorHalt(1)
        endif
        xregr(IGRegr%mol(ih2o):IGRegr%mol(ih2o)+NGRegr%mol(ih2o)-1)= &
             exp(xregr(IGRegr%mol(ih2o):IGRegr%mol(ih2o)+NGRegr%mol(ih2o)-1))

        xgesg(1:npargregr)=xregr(1:npargregr)
        call putScene(ncid=U_retrRegr,xid=IdProf,lat=lat,lon=lon,time=time, &
             x=xGesG(1:nparGRegr),pland=plandAvg(ifor),landType=igeo,eia=eia)
     endif

     if(PHYSon) then
        
        call classatm(Ym,kchan,lat,lon,eia(1),time,igeo,iclassatm)
 
        ! Disable retrieval of cloud parameters per cloud mask
        if (flagDynCldMask) call cloudModeTest(cloudType,retrievalFlags,flagDynCldMask)

        ! Retrieve parameters based on the values of the flags in retrievalFlags
        ! Currently it is all ON.  Need to change that.
        CALL IRsetRetrVec(retrievalFlags, NR_orig, nEmMw_orig, nEmIR_orig, nPar, &
             NR, IR)
        !----Select/combine/tune background cov
        call IRsetFOVbkg(iclassatm,igeo,psfc(ifor), &
                      vCoordTyp,press,umtx(1:nPar,1:nParG),vmtx(1:nParG,1:nPar))
        if (loadExtDone) call IRsetExtBkg(xbakg_Ext(:,ifor))
        if (loadEmisDone) call IRsetEmisBkg(xbakg_Emis(:,ifor),nemirExt)
        IF (plandAvg(ifor) > plandMaxOc) THEN
           isLand=.TRUE.
        ELSE
           isLand=.FALSE.
        ENDIF
        call IRcombBkg(IGExt,NGExt,psfc(ifor),isLand, &
             molTran,vCoordTyp,press,xbakg,ut_cov_net)
        call IRtuneBkg(xbakg,ut_cov_net,iclassatm,igeo)
        call IRtransBkg(ut_cov_net,sx(1:nParMax,1:nParMax), &
             umtx(1:nPar,1:nParG))
        
        if (NR%cldLiq > 0) then
           clw_cov(1:NR%cldLiq,1:NR%cldLiq) = &
                Sx(IR%cldLiq:IR%cldLiq+NR%cldLiq-1,IR%cldLiq:IR%cldLiq+NR%cldLiq-1)
        endif
        if (NR%cldIce > 0) then
           ciw_cov(1:NR%cldIce,1:NR%cldIce) = &
                Sx(IR%cldIce:IR%cldIce+NR%cldIce-1,IR%cldIce:IR%cldIce+NR%cldIce-1)
        endif
        
        xGes=0.
        if(NR%cldLiq > 0)xGes(IR%cldliq+2)=cldamt1stguess
        if(NR%cldIce > 0)xGes(IR%cldIce+2)=iceamt1stguess

        !---Fwd model on 1st guess
        if(REGRESon) then
           call map_regr2retr_geo(xregr,NParGRegr,xGes(1:nPar),xBakG,&
                xGesG,EmMw,TSfc,umtx(1:nPar,1:nParG),NParG,NPar, &
                IG,NG,pref,nlev,nchanmw,IEmMwG,nmol,molTran,IEmIRG,nemIrG,EmRf)
        else
           call IRmap_retr2geo(xGes(1:nPar),xBakG,xGesG,EmMw,TSfc, &
                     vmtx(1:nParG,1:nPar),nParG,nPar,IG,NG,vCoordTyp, &
                     pref,nlev,nchanmw,IEmMwG,ih2o,nmol,molTran,IEmIRG,nemIrG, &
                     EmRf,xOutRangeFlag,logSize)
           if (xOutRangeFlag) &
              print *,'warning:[retr] x out of range may affect convergence'
        endif
        rms_mw   = MISSING_REAL
        chisq_mw = MISSING_REAL
        rms_ir   = MISSING_REAL
        chisq_ir = MISSING_REAL

        xktG=0.  ! To initialize portions not filled later

        if (MWon) then
           IF (TRIM(vCoordTyp) == Pcoord) THEN
              call ossdrv_mw(xGesG,EmMw,eia(1),Y,xktG,xkEmMw,kchan, &
                plandAvg(ifor))
           ELSE
              call ossdrv_mw(xGesG,EmMw,eia(1),Y,xktG,xkEmMw,kchan, &
                plandAvg(ifor),puser=press,iTempG_in=IG%temp, &
                iTskinG_in=IG%tskin,iPsfcG_in=IG%psfc, &
                ImolG_in=IG%mol(1:ih2o),ICldLiqG_in=IG%cldliq)
           ENDIF
           call getChiSq(Ym,Y,DevNoise,chiSqMw,rmsMw,kchan,nchanmw,nch)
           itermw = 0
           print "(' MW iter.   RMS(',f11.7,')      NormChiSq')",sumMWdevNoise
           print '(i9,f19.7,f14.3)',itermw,rmsMw,chiSqMw

           !---Start mw iteration:
           MwIterLoop: do itermw=1,nxiterMw
              CALL set_MW_Invert(xGesG,Ym,Y,xktG,xkEmMw,DevNoise, &
                     xkt(1:nPar,1:nch),Sy,delY,nch,itermw,kchan, &
                     vmtx(1:nParG,1:nPar),nchanmw,dradChan,alpha1slope, &
                     alpha1int,nparG,nPar,IG,NG,iemmwG,nmol,molTran, &
                     vCoordTyp,pref,logSize)
              
              call invrt2(xkt(1:npar,1:nch),Sy(1:nch),Sx(1:npar,1:npar), &
                   delY(1:nch),xGes(1:npar),xGesNew(1:npar),NPar,nch)

              xGes = xGesNew 
              call IRmap_retr2geo(xGesNew(1:nPar),xBakG,xGesG,EmMw,TSfc, &
                    vmtx(1:nParG,1:nPar),nParG,nPar,IG,NG,vCoordTyp,pref,nlev, &
                    nchanmw,IEmMwG,ih2o,nmol,molTran,IEmIRG,nemIrG,EmRf, &
                    xOutRangeFlag,logSize) 
              if (xOutRangeFlag) &
                 print *,'warning:[retr] x out of range may affect '// &
                     'convergence in MW loop'
              call chkges(U_algconfig,F_algconfig,xGesG,ifail, &
                   Sx(1:nPar,1:nPar),clw_cov,ciw_cov,nlev,press, &
                   IG,NG,IR,NR,ih2o,debug)
              IF (TRIM(vCoordTyp) == Pcoord) THEN
                 call ossdrv_mw(xGesG,EmMw,eia(1),Y,xktG,xkEmMw,kchan, &
                   plandAvg(ifor))
              ELSE
                 call ossdrv_mw(xGesG,EmMw,eia(1),Y,xktG,xkEmMw,kchan, &
                   plandAvg(ifor),puser=press,iTempG_in=IG%temp, &
                   iTskinG_in=IG%tskin,iPsfcG_in=IG%psfc, &
                   ImolG_in=IG%mol(1:ih2o),ICldLiqG_in=IG%cldliq)
              ENDIF
              call getChiSq(ym,y,DevNoise,chiSqMw,rmsMw,kchan,nchanmw,nch)
              CALL map_geo2retr(xGesG,xBakG,Tsfc,xGes(1:nPar), &
                                umtx(1:nPar,1:nParG),NParG,NPar,IG,NG, &
                                nmol,molTran,vCoordTyp,logSize)
              print '(i9,f19.7,f14.3)',itermw,rmsMw,chiSqMw
              rms_mw(itermw)   = rmsMw
              chisq_mw(itermw) = chiSqMw
              
              if (chiSqMw < chisqconvMw .and. iStrat /= 3) then
                 exit MwIterLoop  
              elseif(chiSqMw < chisqconvMw .and. iStrat == 3)then
                 if(iexit == 0) then
                    iexit=1
                    
                    call setSfc(igeo,plandAvg(ifor),plandMaxOc)
                    
!cccc
     igeo=2
     print *,'HOLDING igeo=2'
!cccc                    print "('switching to stratified Cov,','igeo = ',i2)",igeo
                    
                    goto 100
                 else
                    exit MwIterLoop 
                 endif
              endif
        
           enddo MwIterLoop
           itermw=min(itermw,nxiterMw)
           iter_tot=itermw
           call putScene(ncid=U_retr,xid=IdProf,lat=lat,lon=lon,time=time, &
                x=xGesG(1:nparGAtm),pland=plandAvg(ifor),landType=igeo,eia=eia, &
                iter=iterMw,rms=rms_mw,chisq=chisq_mw,noise=sumMWdevNoise, &
                EmMW=EmMw,roUnits=roUnits)
           IF (TRIM(vCoordTyp) /= Pcoord) call putScene(ncid=U_retr,pressure=press)
           
        else
           iter_tot=0  ! to initialize for output
        endif  ! MWon
     
        !---Save the retrieval of the FOR profile for later storage--
        xGesMwSav(1:NParG)=xBakG(1:NParG)
        xGesMwSav(IG%mol(iH2o):IG%mol(iH2o)+NG%mol(iH2o)-1)=&
             EXP(xBakG(IG%mol(iH2o):IG%mol(iH2o)+NG%mol(iH2o)-1))
        !-- get surface temperature 
        IF (TRIM(vCoordTyp) == Pcoord) THEN
!           IF (ValidRad) THEN !no valid psfc loaded => interp fails; activate when QC implemented
              call lvl_int(xGesMwSav,pref,Nlev,xGesMwSav(IG%Psfc),xSfc)
!           ELSE
!              xsfc = xBakG(Nlev)
!           ENDIF
           xGesMwSav(IG%Tskin)=xBakG(IG%Tskin)+xsfc
        ELSE
           xGesMwSav(IG%Tskin)=xBakG(IG%Tskin)
        ENDIF

        IF (TRIM(vCoordTyp) == Pcoord) THEN
           call ossdrv_ir(xgesG(1:npargAtm),EmRf,eia(1),sunang=90.,&
             y=y(nchanmw+1:nchan),xkt=xktg(1:npargAtm,nchanmw+1:nchan),&
             xkEmRf=xkEmRf,pland=plandAvg(ifor))
        ELSE
           call ossdrv_ir(xgesG(1:npargAtm),EmRf,eia(1),sunang=90.,&
             y=y(nchanmw+1:nchan),xkt=xktg(1:npargAtm,nchanmw+1:nchan),&
             xkEmRf=xkEmRf,pland=plandAvg(ifor),puser=press, &
             iTempG_in=IG%temp,iTskinG_in=IG%tskin,iPsfcG_in=IG%psfc, &
             iMolInd_in=IG%mol(1:nmol))
        ENDIF
        if (enableScatter)  call finiteDiff(xgesG(1:npargAtm),EmRf,&
             eia(1),90.,y(nchanmw+1:nchan),xktg(1:npargAtm,nchanmw+1:nchan),&
             xkEmRf,plandAvg(ifor),press,NR)
        
        chiSq0 = MISSING_REAL
        
        dSx=real(Sx,KIND=DBL_INVT)
        call dminv(dSx(1:Npar,1:NPar),npar,detmnt)
        Sx=real(dSx)

        if (NR%cldLiq > 0) then
           clw_cov(1:NR%cldLiq,1:NR%cldLiq) = &
                Sx(IR%cldLiq:IR%cldLiq+NR%cldLiq-1,IR%cldLiq:IR%cldLiq+NR%cldLiq-1)
        endif
        if (NR%cldIce > 0) then
           ciw_cov(1:NR%cldIce,1:NR%cldIce) = &
                Sx(IR%cldIce:IR%cldIce+NR%cldIce-1,IR%cldIce:IR%cldIce+NR%cldIce-1)
        endif
        
        IRIterLoop: do iter=1,nxiter
           iter_tot=iter_tot+1

           call set_IRMW_Invert(xGesG,Ym,Y,xktG,xkEmMw, &
                DevNoise(1:nchan),AtmNoise,dradChan,xkt,xkEmRf,Sy,delY,nch, &
                iter,nmol,molTran, &
                kchan,vmtx(1:nParG,1:nPar),nparG,nPar,nchanmw,nchan, &
                IG,NG,IEmMwG,IEmIrG,nemIrG,vCoordTyp,pref,Airs_Err,logSize)
        
           call invrt1(xkt(1:Npar,1:nch),Sy(1:nch),Sx(1:Npar,1:NPar), &
                delY(1:nch),xGes(1:nPar),xGesNew(1:Npar),NPar,nch)
          
           Ylast=Y
           xGesLast=xGes

           xges=xgesNew
           
           call IRmap_retr2geo(xGesNew(1:nPar),xBakG,xGesG,EmMw,TSfc, &
                vmtx(1:nParG,1:nPar),nParG,nPar,IG,NG,vCoordTyp,pref, &
                nlev,nchanmw,IEmMwG,ih2o,nmol,molTran,IEmIRG,nemIrG,EmRf, &
                xOutRangeFlag,logSize)
           if (xOutRangeFlag) &
              print *,'warning:[retr] x out of range may affect '// &
                  'convergence in main loop'
           TsfcLast=TSfc

           call chkges(U_algconfig,F_algconfig, &
                xGesG,ifail,Sx(1:nPar, 1:nPar),clw_cov,ciw_cov,nlev,press, &
                IG,NG,IR,NR,ih2o,debug,invCldCov_in=.true.)
           
           IF (TRIM(vCoordTyp) == Pcoord) THEN
              if (MWon)call ossdrv_mw(xGesG,EmMw,eia(1),Y,xktG,xkEmMw,kchan, &
                plandAvg(ifor))
              call ossdrv_ir(xgesG(1:npargAtm),EmRf,eia(1),sunang=90., &
                y=y(nchanmw+1:nchan),xkt=xktg(1:npargAtm,nchanmw+1:nchan), &
                xkEmRf=xkEmRf,pland=plandAvg(ifor))
           ELSE
              if (MWon)call ossdrv_mw(xGesG,EmMw,eia(1),Y,xktG,xkEmMw,kchan, &
                plandAvg(ifor),puser=press,iTempG_in=IG%temp, &
                iTskinG_in=IG%tskin,iPsfcG_in=IG%psfc, &
                ImolG_in=IG%mol(1:ih2o),ICldLiqG_in=IG%cldliq)
              call ossdrv_ir(xgesG(1:npargAtm),EmRf,eia(1),sunang=90., &
                y=y(nchanmw+1:nchan),xkt=xktg(1:npargAtm,nchanmw+1:nchan), &
                xkEmRf=xkEmRf,pland=plandAvg(ifor),puser=press, &
                iTempG_in=IG%temp,iTskinG_in=IG%tskin,iPsfcG_in=IG%psfc, &
                iMolInd_in=IG%mol(1:nmol))
           ENDIF

           if (MWon) then
              call getChiSq(ym,y,DevNoise,chiSqMw,rmsMw,kchan,nchanmw,nch)
           else
              chisqMw=MISSING_REAL
           endif
           call getChiSq(Ym,Y,DevNoise,chiSq,rms,kchan,nchan,nch)
           
           rms_ir(iter)=rms
           chisq_ir(iter)=chisq
           
           
           chiSqRatio=abs((ChiSq0-chiSq)/chiSq0)
           chiSq0=chiSq
           
           if (chiSq < chisqconv .and. chiSqRatio < chisqRatioconv)exit IRIterLoop

           if (enableScatter)call finiteDiff(xgesG(1:npargAtm),EmRf,eia(1),90.,&
             y(nchanmw+1:nchan), xktg(1:npargAtm,nchanmw+1:nchan),&
             xkEmRf,plandAvg(ifor),press,NR)
           
           CALL map_geo2retr(xGesG,xBakG,Tsfc,xGes(1:nPar),umtx(1:nPar,1:nParG), &
                     nParG,nPar,IG,NG,nmol,molTran,vCoordTyp,logSize)

        enddo IRIterLoop
        iter=min(iter,nxiter)
        
        IF (genLinInvert) THEN
           call buildLinearCoef(xkt(1:Npar,1:nch),Sy(1:nch),Sx(1:Npar,1:NPar), &
                Ylast,kchan,xGesLast(1:nPar),vmtx(1:nParG,1:nPar), &
                nPar,nch,nParGon,xBakG,xOffG,xGainG(:,1:nch))
           xOffGall(:,ifor)=xOffG
           xGainGall(:,1:nch,ifor)=xGainG(:,1:nch)
           chiSqAll(ifor)=chiSq
           TsfcAll(ifor)= TsfcLast
        ENDIF

        call putScene(ncid=U_retrIR,xid=IdProf,lat=lat,lon=lon,time=time, &
             x=xGesG(1:nparGAtm),pland=plandAvg(ifor),landType=igeo,eia=eia, &
             iter=iter,rms=rms_ir,chisq=chisq_ir,noise=sumIRdevNoise, &
             EmRf=EmRf)
        IF (TRIM(vCoordTyp) /= Pcoord) call putScene(ncid=U_retrIR,pressure=press)
        if (MWon) call putScene(ncid=U_retrIR,EmMW=EmMw)
        
        call putScene(ncid=U_guess,xid=IdProf,lat=lat,lon=lon,time=time, &
             x=xGesMwSav(1:nparGAtm),pland=plandAvg(ifor),landType=igeo,eia=eia, &
             EmRf=xGesMwSav(IEmIRG:IEmIRG+2*nEmIRG-1))
        IF (TRIM(vCoordTyp) /= Pcoord) call putScene(ncid=U_guess,pressure=press)
        if (MWon) call putScene(ncid=U_guess, &
             EmMW=xGesMwSav(IEmMwG:IEmMwG+nchanmw-1))
        
     endif ! PHYSon
     
  enddo ProfLoop

  if(REGRESon) call closeScene(U_retrRegr)
  if(MWon)     call closeScene(U_retr)
  if(PHYSon)   call closeScene(U_retrIR)
  call closeScene(U_guess)
  if(MWon)     call closeRad(U_rad)
  if (enableScatter) call destroyFiniteDiff()
  IF (genLinInvert) THEN
    call writeLinInvert(F_linInvert)
    call destroyLinInvert()
    deallocate (xOffG,xGainG,xOffGall,xGainGall,chiSqAll,TsfcAll)
  ENDIF
  IF (loadExtDone) deallocate(xbakg_Ext)
  IF (loadEmisDone) deallocate(xbakg_Emis,wvnExt)
  deallocate(Ylast,xGesLast)
  DEALLOCATE(press)
  deallocate(dSx)
  call oss_destroy()
  call oss_destroy_ir()
  CALL destroyNoise()

  deallocate(psfc)
  deallocate(plandAvg)
  if (allocated(pref)) deallocate(pref)
  if (allocated(frq)) deallocate(frq)
  if (allocated(pol)) deallocate(pol)
  if (allocated(frqEmIR)) deallocate(frqEmIR)
  if (associated(prefRegr)) deallocate(prefRegr)
  if (allocated(chanNum)) deallocate(chanNum)
contains 

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

  !--I/O variables
  character(len=*),                  intent(in)     :: file_algcfg

  !---retrTuning Namelist items
  integer            :: ngas(maxMol)
  integer            :: nTskin=0 ! Overridden if ntemp>0 in setXPtrR
  integer            :: ntemp,nemmw,nCldLiq,nCldIce,nxiterMw,nxiter
  real               :: drad,alpha1slope,alpha1int
  real               :: chisqconvMw,chisqconv,chisqRatioconv
  real               :: cldamt1stguess,iceamt1stguess
  real               :: TightCLDcov
  real               :: plandMaxOc
  logical            :: logSize

  namelist /retrtuning/ nxiterMw,nxiter,drad,alpha1slope,alpha1int, &
       chisqconvMw,chisqconv,chisqRatioconv,cldamt1stguess, &
       iceamt1stguess, &
       nTskin,ntemp,ngas,nemmw,nemir,nCldLiq,nCldIce,TightCLDcov, &
       logSize
  namelist /IRsfcClass/pLandMaxOc

  if(loadAlgCfgDone) return

  !------------------------------------------------------------------------
  ! Read namelist file items into local variables
  !------------------------------------------------------------------------
  open(U_algconfig,file=file_algcfg)
  read(U_algconfig,retrtuning)
  close(U_algconfig)

  read(*,IRsfcClass) !separated out the pLandMaxOc due to the impact
                     !from the change in EDSR modularization

  !------------------------------------------------------------------------
  ! Set global variables
  !------------------------------------------------------------------------
  call setPropertiesAlgConfig(ngas,nTskin,ntemp,nemmw,nemir,nCldLiq,nCldIce, &
       nxiterMw,nxiter,drad,alpha1slope,alpha1int,chisqconvMw, &
       chisqconv,chisqRatioconv,cldamt1stguess,iceamt1stguess,TightCLDcov, &
       plandMaxOc,logSize)

  return
end subroutine loadAlgconfig

  !--
  ! Read items from file into global memory
  !--
subroutine loadAux(file_aux)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   loadAux
!
! PURPOSE:
!
!   Load auxiliary data from file.
!
! SYNTAX:
!
!   CALL loadAux(file_aux)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   file_aux  CHAR  File path for auxiliary data
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

  !--I/O variables
  character(len=*),intent(in)   :: file_aux
  !--Local variables
  integer                       :: ncidtmp,openStatus
  logical                       :: is_ncdf
  integer                       :: i,nscans
  character(len=8)              :: dummy
  real,dimension(:),allocatable :: psfcTmp
  type(StateIndex_t)            :: IGAux,NGAux

  if(loadAuxDone) return
  
  !---Test to determine if aux file is netCDF:
  openStatus = nf_open(file_aux,nf_nowrite,ncidtmp)
  is_ncdf = (openStatus == nf_noerr)
  if ( is_ncdf ) call closeNcdfFile(ncidtmp)
  
  if ( is_ncdf ) then
     !---Reading of netCDF auxiliary file:
     call openScene(ncid=U_Aux,file=file_aux,IG=IGAux,NG=NGAux,status='old')
     allocate(psfc(nprofs),plandAvg(nprofs))
     allocate( psfcTmp(getVectorLength(NGAux)) )
     do i=1,nprofs
        call getScene(ncid=U_Aux,irec=i,x=psfcTmp,pland=plandAvg(i))
        psfc(i) = psfcTmp(IGAux%psfc)
     enddo
     deallocate(pSfcTmp)
     call closeScene(U_Aux)
  else
     !---Reading of ASCII auxiliary file, for backward compatibility:
     open(U_Aux, file=file_aux, status='old')
     read(U_Aux,'(2i6)') nscans
     allocate(psfc(nprofs),plandAvg(nprofs))
     do i=1,nprofs
        read(U_Aux,'(a8,2f10.3)')dummy,psfc(i),plandAvg(i)
     enddo
     close(U_aux)
  endif

  loadAuxDone=.true.
  return
end subroutine loadAux

  SUBROUTINE adjustChanSet(radTmp,wvnRad,wvnOSS,nchanRad,nChanIR,kchan,radAll, &
     chanIDtmp,chanID)
    

!<f90Subroutine>********************************************************
!
! NAME:
!
!   adjustChanSet
!
! PURPOSE:
!
!   Make adjustments, if OSS channel set does not match the radiance 
!   file channel set.
!
! SYNTAX:
!
!   CALL adjustChanSet(radTmp, wvnRad, wvnOSS, nchanRad, nChanIR, 
!      kchan, radAll, chanIDtmp, chanID)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   radTmp     REAL     Input set of radiances 
!   wvnRad     REAL     Center wavenumbers of channels from rad 
!                       file 
!   wvnOSS     REAL     Center wavenumbers of all channels provided 
!                       by OSS 
!   nchanRad   INTEGER  Number of channels in rad data source
!   nChanIR    INTEGER  Number of IR channels
!   kchan      INTEGER  Channel on/off mask
!   chanIDtmp  CHAR     Channel identifier string temporary version
!   
!   INPUTS/OUTPUTS:
!   
!   radAll     REAL     Radiances for all locations and channels
!   chanID     CHAR     Channel identifier string
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    ! Subroutine to make adjustments if OSS channel set does not match
    !   the radiance file channel set.
    real,             dimension(:,:) ,intent(in)    :: radTmp
    real,             dimension(:  ) ,intent(in)    :: wvnOSS
    integer,          dimension(:),   intent(in)    :: kchan
    real,             dimension(:  ) ,intent(in)    :: wvnRad
    real,             dimension(:,:) ,intent(inout) :: radAll
    integer,                          intent(in)    :: nchanRad,nChanIR
    character(len=*), dimension(:),   intent(in)    :: chanIDtmp
    character(len=*), dimension(:),   intent(inout) :: chanID

    !---Local variables:
    integer                            :: nChanOSS,iOSS,iRad
    real   ,parameter                  :: wvThresh=5., missingRad=0.

    !---Loop over OSS and radfile channel sets.  Organize the radfile
    !    channels into the same sequence as the OSS channels.  Fill
    !    radiance data with a missing value if it is unavailable.
    ossLoop: do iOSS=1,nChanIR
       radLoop: do iRad=1,nChanRad
          if ( ABS(wvnOSS(iOSS)-wvnRad(iRad)) < wvThresh ) then
             radAll(iOSS,:) = radTmp(iRad,:)
             chanID(iOSS) = chanIDtmp(iRad)
             exit radLoop
          endif
          if ( iRad == nChanRad ) then
             if (kchan(iOSS) > 0) then
                print*,'Err[retr::adjustChanSet]: '
                print*,'Radiance data unavailable for channel # ',iOSS
                print*,'Adjust input kchan vector to turn off this channel.'
                call errorHalt(1)
             endif
             radAll(iOSS,:) = missingRad
             chanID(iOSS) = repeat(' ',len(chanID(iOSS)))
          endif
       enddo radLoop
    enddo ossLoop

    return

  END SUBROUTINE adjustChanSet

  subroutine writeLinInvert(fileLinInv)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   writeLinInvert
!
! PURPOSE:
!
!   Write data for linear inversion to a file
!
! SYNTAX:
!
!   CALL writeLinInvert(fileLinInv)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   fileLinInv  CHAR  Path of file of linear inversion coefficients
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

  !  Write data for linear inversion to a file

    character(len=*),  intent(in)    :: fileLinInv     

    call openScene(U_linInvert,fileLinInv,MolId=molId, &
      nLevel=nlev,status='new',IG=IGRout,NG=NGRout)
    CALL writeNcdfAttr(U_linInvert, &
       attrName='VertCoordType',attr=TRIM(vCoordTyp))
    IF (TRIM(vCoordTyp) == Pcoord) CALL writeNcdfAttr(U_linInvert, &
       attrName='standardPressureGrid',attr=pref)

    call writeNcdfAttr(U_linInvert,attrName='nchan',attr=nchan)
    call writeNcdfAttr(U_linInvert,attrName='kchan',attr=kchan(1:nchan))
    call writeNcdfAttr(U_linInvert,attrName='chiSqThres',attr=chisqconv)
    call writeNcdfAttr(U_linInvert,attrName='wavenumbers', &
         attr=wvnOSS(1:nChanIR))

    call writeNcdfData(U_linInvert,xOffGAll, &
         varname='xOff', &
         varlenName=(/'nPar  ','nGroup'/), &
         varLongName='Linear inversion offset state vector')

    call writeNcdfData(U_linInvert,xGainGAll(:,1:nch,:), &
         varname='xGain', &
         varlenName=(/'nPar  ','nChOn ','nGroup'/), &
         varLongName='Linear inversion gain matrix')

    call writeNcdfData(U_linInvert,chiSqAll, &
         varname='chiSqFinal', &
         varlenName=('nGroup'), &
         varLongName='Chi-Squared at last iteration', &
         varUnit='unitless')

    call writeNcdfData(U_linInvert,psfc, &
         varname='pSfc', &
         varlenName=('nGroup'), &
         varLongName='Surface pressure', &
         varUnit='mb')

    IF (NR_orig%Tskin > 0) THEN
      call writeNcdfData(U_linInvert,TsfcAll, &
           varname='Tsfc', &
           varlenName=('nGroup'), &
           varLongName='Air temperature at the surface', &
           varUnit='K')
    ENDIF

    call closeScene(U_linInvert)

    return

  end subroutine writeLinInvert

  subroutine loadExt(file_ext)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   loadExt
!
! PURPOSE:
!
!   Load external background data from file
!
! SYNTAX:
!
!   CALL loadExt(file_ext)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   file_ext  CHAR  Path of file of external background data
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !--I/O variables
    character(len=*),       intent(in)     :: file_ext
    integer                                :: n
    real, dimension (:), allocatable       :: xExtTmp
    if(loadExtDone) return
    molIDext=nullMolID
    call openScene(ncid=U_bkgrdExtern,file=file_ext,nprf=nforExt, &
         molID=molIDext,status='old',IG=IGExt,NG=NGExt)
  
    nparGExt = getVectorLength(NGExt)
    allocate (xbakg_Ext(nparGExt,nforExt),xExtTmp(nparGExt))
    do n=1,nforExt
       CALL getScene(ncid=U_bkgrdExtern,irec=n,x=xExtTmp)
       xbakg_Ext(:,n)=xExtTmp(:)
    enddo
    deallocate(xExtTmp)
    call closeScene(U_bkgrdExtern)
    loadExtDone=.true.
    return
  end subroutine loadExt

  subroutine loadEmis(file_emis)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   loadEmis
!
! PURPOSE:
!
!   Load external emissivity data for use in background
!
! SYNTAX:
!
!   CALL loadEmis(file_emis)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   file_emis  CHAR  Path of file of external emissivity data
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

  ! Load external emissivity data for use in background

    USE bkg_io_module, ONLY: queryCov, openCov, getCov, closeCov

    !--I/O variables
    character(len=*),       intent(in)     :: file_emis

    integer                                :: i,n,nqcEmis
    real,    dimension (:),   allocatable  :: xEmisTmp
    integer, dimension (:),   allocatable  :: qcEmis

    if(loadEmisDone) return

    call queryCov(file=file_emis,nbkg=nforEmis,nemir=nemirExt,nqc=nqcEmis)
    allocate (wvnExt(nemirExt))
    call openCov(ncid=U_IRemissExt,file=file_emis, &
         frqEmir=wvnExt,status='old')
    
    allocate (xbakg_Emis(2*nemirExt,nforEmis),xEmisTmp(nemirExt),qcEmis(nqcEmis))
    do n=1,nforEmis
       CALL getCov(ncid=U_IRemissExt,irec=n,dmEmIR=xEmisTmp,qc=qcEmis)
       xbakg_Emis(1:nemirExt,n)=xEmisTmp(:)
       if (btest(qcEmis(1),0)) xbakg_Emis(:,n)=MISSING_REAL
    enddo
    deallocate(xEmisTmp,qcEmis)
    call closeCov(U_IRemissExt)
    loadEmisDone=.true.
    return

  end subroutine loadEmis

!
!Notes:
!
!      cloudModeTest is made part of IRBkgPrepModule
!
!      operations of finite differences have been constructed as an
!      independent module, FiniteDiffModule, under rslib/ir/rtm
!
! 07/27/2016, yhe@aer.com
!
end program retr
