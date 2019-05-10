!--------------------------------------------------------------
!
! SIM.f90:  Driver program for performing retrievals
!
! Atmospheric and Environmental Research, Inc., 2004.
!
!--------------------------------------------------------------

PROGRAM simScat
  !---This version of sim is modified for use with the cloud module
  !---and scattering RTM.
  USE StateIndexModule, ONLY: StateIndex_t, maxMol, whereH2O, getVectorLength, &
       charToMolID
  USE ncdf_module, ONLY: readNcdfDim
  USE scene_io_module
  USE rad_io_module
  USE IOutilities, ONLY: protectNonStdRadfiles
  USE oss_mw_module, ONLY : ossinit_mw, ossdrv_mw, oss_destroy
  USE OSSMWmoduleSubs, ONLY :  getOSSselMW,&
                       set2ndOrderPlanck,getPlanckLinCtr
  USE NoiseModule
  USE IRNoiseModule
  USE IRChannels
  USE ReadStdInputs
  USE IRReadStdInputs
  USE oss_ir_module_scat, only : OSSdestroyIR, ossdrv_ir, ossinit_ir
  USE constants, ONLY: MISSING_CHAR,MISSING_INT

  USE CloudModule, ONLY : iCLW,iIce,destroyCloud,typeSlab

  USE CloudParameters
  USE CloudDataStruct

  IMPLICIT NONE

  !----Arrays
  INTEGER, PARAMETER                      :: mxchan=2600,mxlev=101
  REAL,    DIMENSION(mxchan)              :: Ym,xkEmMw,rerr
  REAL,    DIMENSION(mxchan)              :: tbForNoise
  REAL,    DIMENSION(mxchan)              :: alphaCtr,betaCtr
  REAL,    DIMENSION(mxlev)               :: Pref_scene
  REAL,    DIMENSION(:), allocatable      :: pref
  REAL,    DIMENSION(mxchan)              :: frq,frq_scene
  INTEGER, DIMENSION(6)                   :: time
  LOGICAL, DIMENSION(mxchan)              :: kchan
  INTEGER, DIMENSION(mxchan)              :: pol,pol_scene
  INTEGER, DIMENSION(:),     allocatable  :: molID

  !----IR-related arrays
  REAL,    DIMENSION(:),     allocatable  :: frqEmIR,EmRf
  REAL,    DIMENSION(:,:),   allocatable  :: xktG
  REAL,    DIMENSION(:,:,:), allocatable  :: xkEmRf
  REAL,    DIMENSION(mxchan)              :: frq_ir
  REAL,    DIMENSION(:,:),   allocatable  :: EmRfdrv
  REAL,    DIMENSION(:,:,:), allocatable  :: xkEmRfdrv

  integer                                 :: nparGAtm,nEmIR,nchanOSS
  integer                                 :: nlevOSS
  real, dimension(:), allocatable         :: prefOSS
  INTEGER, PARAMETER                      :: lenIDir=12
  CHARACTER(LEN=lenIDir),DIMENSION(:),ALLOCATABLE :: chanIDir

  !-- Noise related variables:
  real                                    :: sumIRDevNoise

  !----Single variables
  LOGICAL            :: debug,MWon,addNoise=.TRUE.,DynamicNoise
  LOGICAL            :: REGRESon,PHYSon,ossScatRegrIRon
  LOGICAL            :: flagDynCldMask
  LOGICAL            :: use2ndOrderPlanck
  REAL               :: lat,lon,pland
  REAL               :: totalNoise
  real               :: deltaP=0.05,deltaFrq=0.001
  CHARACTER*8        :: xid
  INTEGER            :: nmol
  INTEGER            :: ih2o,nfor,nprofs,mwCld
  INTEGER            :: icell,icasc,iembkgflg,iextbkg,ifor
  INTEGER            :: nlev,nlevel,nchan,nchanmw,nchanir
  INTEGER            :: iStrat
  INTEGER,parameter  :: maxhyd=5
  integer,parameter  :: nHyd=2
  INTEGER            :: iRadOrTb
  CHARACTER(LEN=12)  :: molOnRT(maxMol)
  LOGICAL            :: enableScatter
  TYPE(StateIndex_t) :: IG, NG
  integer                                :: nAng,xidLen,timeDim
  real, dimension(:), allocatable        :: eia,EmMW,x
  LOGICAL            :: copyTime=.false.,copyXid=.false.

  logical            :: lambertian = .false.
  integer            :: cldTyp=typeSlab !'slab'
  CHARACTER (LEN=mxNameLength), DIMENSION(:), allocatable :: propertyNames
  TYPE(CloudScene), DIMENSION(:), ALLOCATABLE :: cldScene
  integer                                 :: iHydr, iLyr, ipp, ippSizP
  integer                                 :: nLevHydr, npProperty, nHydr
  integer,allocatable,dimension(:)        :: chanIndex
  integer                                 :: jj, ichan
  ! obsLevel the level at which the instrumnet is located
  ! TOA corresponds to obsLevel=1
  integer                                 :: obsLevel
  obsLevel=1


  !------------------------------------------
  !      Read Global Data to be shared
  !------------------------------------------
  call readStd(nprofs,debug,mwCld,iCell,icasc, &
       iembkgflg,iextbkg, &
       kchanOut=kchan,DynamicNoiseIn=DynamicNoise,addNoiseIn=addNoise, &
       enableScatterIn=enableScatter, &
       use2ndOrderPlanckIn=use2ndOrderPlanck,molOnRTin=molOnRT)
  call IRreadStd(REGRESON,PHYSON, &
       MWon,iStrat,flagDynCldMask,ossScatRegrIrIn=ossScatRegrIRon)

  !----Derive nmol and molID from inputs
  nmol=COUNT(molOnRT /= MISSING_CHAR)
  if (nmol <= 0) then
     print *,'err:[retr] No molecules selected for RT'
     call errorHalt(1)
  endif
  allocate (molID(nmol))
  call charToMolID(molOnRT(1:nmol),molID)
  ih2o=whereH2O(molID)

  if (MWon) then
     CALL queryScene(file=F_scene,nchmw=nchanmw)
     CALL queryScene(file=F_scene,polarity=pol_scene(1:nchanmw),&
          freq=frq_scene(1:nchanmw))
     allocate(EmMW(nchanmw))
  endif
  CALL queryScene(file=F_scene,nLevel=nlevel,nAng=nAng, &
       nSfcGrid=nEmIR)
  allocate(EmRf(2*nEmIR),frqEmIR(nEmIR), EmRfdrv(nEmIR,2))
  CALL openScene(ncid=U_scene,file=F_scene,nprf=nfor,IG=IG,NG=NG,&
       pressure=pref_scene(1:nlevel),MolId=molId,sfcGrid=frqEmIR)

  nHydr    = 2
  nLevHydr = 0

  !!--Fill in the non-missing hydrometeor types:
  allocate(cldScene(nHydr))
  cldScene(iCLW)%name = 'CldLiq '
  cldScene(iIce)%name = 'CldIce '
! npProperty has to be consistent with the number of elements in the
  ! propertyNames
  npProperty = 3
  allocate(propertyNames(npProperty))
  propertyNames(1:npProperty) = (/'SizP', 'SizS', 'Temp'/)
  do ipp = 1, npProperty
     if (trim(propertyNames(ipp)) == 'SizP') then
        ippSizP = ipp
        exit
     end if
  end do

  do iHydr = 1, nHydr
     cldScene(iHydr)%nLevel = nLevHydr
     cldScene(iHydr)%nLayer = nLevHydr - 1
     cldScene(iHydr)%nProperty = npProperty
     allocate(cldScene(iHydr)%physProp(nLevHydr))
     allocate(cldScene(iHydr)%physDescr(npProperty))
     allocate(cldScene(iHydr)%hydrAmt(nLevHydr))
     DO iLyr = 1, nLevHydr
        cldScene(iHydr)%physProp(iLyr)%nProperty = npProperty
        allocate(cldScene(iHydr)%physProp(iLyr)%value(npProperty))
        cldScene(iHydr)%physProp(iLyr)%value  = 0.0
     end do
     DO ipp = 1, npProperty
        cldScene(iHydr)%physDescr(ipp)%name = trim(propertyNames(ipp))
     ENDDO
  end do

  nparGAtm = getVectorLength(NG)
  allocate(eia(nAng),x(nparGAtm))
  !------------------------------------------
  !         Initialize the OSS model
  !------------------------------------------
  if (MWon) then
     IF (use2ndOrderPlanck) call set2ndOrderPlanck()  ! Instead of default
     CALL ossinit_mw(F_osscoefs,F_ossoptdpth, ih2o, &
          molID(1:ih2o),IG%cldliq+NG%cldliq-1,&
          nchan,Frq,nlev,Pref,pol,iRadOrTb)
     CALL getPlanckLinCtr(alphaCtr,betaCtr)
     pol_scene=pol  !hes fix because of inappropriate scene file
     if (nchan /= nchanmw) then
        print*,'nchan for OSS MW inconsistent with scene file'
        call errorHalt(1)
     endif
  else
     nchanmw=0
  endif ! MWon

  !----Initialize IR OSS module and get back Info about LUTs
  call ossinit_ir(F_osscoefs_ir,F_ossoptdpth_ir,F_solarflux, &
       F_cloudTab, F_ossScatRegr_ir, cldTyp,nHyd,cldScene, &
       ossScatRegrIRon, enableScatter,F_defProfFile,nmol,molid, &
       nchanOSS, chanIndex, Frq_ir, nlevOSS,PrefOSS)

  nchanir=nchanOSS
  nchan=nchanmw+nchanir  ! total number of spectral channels

  allocate(xkEmRf(2,nemir,nchanir))
  allocate(xkEmRfdrv(nemir,nchanir,2))
  allocate(xktg(npargAtm,MAX(nchanmw,nchanir)))
  allocate(chanIDir(nchanir))
  chanIDir(:)=repeat(' ',lenIDir)  ! No IDs available currently
  !------------------------------------------
  !             IR channel selection flags
  !------------------------------------------
  call initIRChannels(nchanmw,nchan,frq_ir,kchan)
  frq(nchanmw+1:nchan)=frq_ir(1:nchanir)

  !------------------------------------------
  !             Consistency Checks
  !------------------------------------------
  if (MWon) then
     IF (nlevOSS  /= nlev)  then
        print*,' Nlev inconsistent in IR (Oss/Bkg)'
        call errorHalt(1)
     endif
  else
     nlev=nlevOSS
     pref(1:nlev)=prefOSS(1:nlev)
  endif

  IF (ANY(abs(pref(1:nlev)-prefOSS(1:nlev)) > deltaP ))then
     print*,' Pressure Grid inconsistent in IR (Bkg/OSS)'
     call errorHalt(1)
  endif

  if (MWon) then
     IF (ANY(pol(1:nchanmw) /= pol_scene(1:nchanmw))) THEN
        PRINT *,'err[sim]:  Polarization inconsistent'
        call errorHalt(1)
     END IF
     IF (ANY(ABS(frq(1:nchanmw)-frq_scene(1:nchanmw)) > deltaFrq )) THEN
        PRINT *,'err[sim]: MW Frequencies inconsistent'
        call errorHalt(1)
     END IF
  endif
  IF((nprofs.LT.1).OR.(nprofs.GT.nfor)) nprofs=nfor
  !------------------------------------------
  !  Open SDR & Auxill files and write hdrs
  !------------------------------------------
  if (Mwon) then
     CALL protectNonStdRadfiles(TRIM(F_Rad))
     CALL openRad(ncid=U_Rad,fname=TRIM(F_Rad),pol=pol(1:nchanmw), &
          frq=frq(1:nchanmw),nchan=nchanmw,status='new')
  endif
  CALL protectNonStdRadfiles(TRIM(F_RadIR))
  CALL openRad(ncid=U_RadIR,fname=TRIM(F_RadIR), &
       wvn=frq(nchanmw+1:nchan),nchan=nchanir,status='new')

  !---Set flags to determine whether optional variables
  !    will be included in I/O:
  timeDim=readNcdfDim(U_scene,'nTime',silent=.true.)
  if (timeDim /= 0) copyTime=.true.
  call queryScene(file=F_scene,xidLen=xidLen)
  if (xidLen  /= 0) copyXid =.true.

  OPEN(U_Aux,file=F_Aux,status='unknown')
  WRITE(U_Aux,'(2i6)') nprofs,1
  !------------------------------------------
  !             MAIN PROFILE LOOP
  !------------------------------------------

  xid=MISSING_CHAR
  time=MISSING_INT
  DO ifor=1,nprofs
     !---READ THE GEOPHYSICAL DATA

     if (MWon) then
        CALL getScene(ncid=U_scene,irec=ifor,EmMw=EmMw(1:nchanmw))
     endif
     CALL getScene(ncid=U_scene,irec=ifor,lat=lat,lon=lon,&
          x=x,pland=pland,eia=eia,EmRf=EmRf)
     if (copyTime) CALL getScene(ncid=U_scene,irec=ifor,time=time)
     if (copyXid ) CALL getScene(ncid=U_scene,irec=ifor,xid=xid)

     print*, 'IFOR = ',ifor

     !---RADIATIVE TRANSFER MODELING
     do ichan=1,nEmIR
        EmRfdrv(ichan,1) = EmRf(2*ichan-1)
        EmRfdrv(ichan,2) = EmRf(2*ichan)
     end do
     if (MWon)CALL ossdrv_mw(x,EmMw,eia(1),obsLevel,Ym,xktG,xkEmMw,kchan,pland,&
                    lat,tempIndex=IG%temp,tSkinIndex=IG%tskin,&
                    pSurfIndex=IG%psfc,varMolIndex=IG%mol(1:ih2o),&
                    ICldLiqG_in=IG%cldliq, lambertian=lambertian)

     call  ossdrv_ir(x,frqEmIR(1:nEmIR),EmRfdrv,eia(1), 90.,0.,&
              obsLevel,y=ym(nchanmw+1:nchan),xkt=xktg,xkEmRf=xkEmRfdrv,&
              lat=lat,tempIndex=IG%temp,tSkinIndex=IG%tskin,pSurfIndex=IG%psfc,&
              varMolIndex=IG%mol(1:nmol))

     do jj=1,2
      do ichan=1,nchanir
        xkEmRf(jj,1:nemir,ichan)=xkEmRfdrv(1:nemir,ichan, jj)
      end do
     end do

     !---ADD INSTRUMENT NOISE
     IF (addNoise) then
        if (MWon) THEN
           IF (iRadOrTb == 0) THEN
              tbForNoise=Ym
           ELSE               ! Convert from radiance to radiance temperature
              tbForNoise(1:nchanmw)= &
                  Ym(1:nchanmw)*alphaCtr(1:nchanmw)+betaCtr(1:nchanmw)
           ENDIF
           CALL deviceNoise(icellres=icell,tb=tbForNoise,rerr=rerr, &
                sum_rerr=totalNoise,kchan=kchan,nchan=nchanmw, &
                DynamicNoise=DynamicNoise)
           IF (iRadOrTb == 1) rerr(1:nchanmw)=rerr(1:nchanmw)/alphaCtr(1:nchanmw)
           CALL addDeviceNoise(Y=Ym,devNoise=rerr,nchan=nchanmw)
        ENDIF
        CALL IRaddDeviceNoise(frq_ir,Ym(nchanmw+1:nchan), &
             rerr(nchanmw+1:nchan),sumIRDevNoise, &
             kchan(nchanmw+1:nchan),nchanir,chanIDir)
     endif

     !---OUTPUT RADIANCES AND AUXILL INFO
     if (Mwon) then
        CALL putRad(ncid=U_Rad,xid=xId,lat=lat,lon=lon,&
             time=Time,eia=eia)
        IF (iRadOrTb == 0) THEN
           CALL putRad(ncid=U_Rad,tb=Ym(1:nchanmw))
        ELSE
           CALL putRad(ncid=U_Rad,radMW=Ym(1:nchanmw))
        ENDIF
     endif
     CALL putRad(ncid=U_RadIR,xid=xId,lat=lat,lon=lon,&
          time=Time,eia=eia,rad=Ym(nchanmw+1:nchan))

     WRITE(U_Aux,'(a8,f10.3,9f10.3)') xid,x(IG%psfc),pland
  ENDDO
  !------------------------------------------
  !    CLOSE SCENE AND AUXILL FILES
  !------------------------------------------
  CALL oss_destroy()
  CALL OSSdestroyIR()
  CALL destroyCloud()
  CALL destroyNoise()
  CALL closeScene(ncid=U_scene)
  if(MWon)CALL closeRad(ncid=U_Rad)
  CALL closeRad(ncid=U_RadIR)
  CLOSE(U_Aux)

  IF (ALLOCATED(molID))   DEALLOCATE(molID)
  IF (ALLOCATED(EmMw))    DEALLOCATE(EmMw)
  IF (ALLOCATED(frqEmIR)) DEALLOCATE(frqEmIR)
  IF (ALLOCATED(eia))     DEALLOCATE(eia)
  IF (ALLOCATED(x))       DEALLOCATE(x)
  IF (ALLOCATED(EmRf))    DEALLOCATE(EmRf)
  IF (ALLOCATED(xkEmRf))  DEALLOCATE(xkEmRf)
  IF (ALLOCATED(xktg))    DEALLOCATE(xktg)
  IF (ALLOCATED(chanIDir))DEALLOCATE(chanIDir)

END PROGRAM simScat


