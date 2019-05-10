MODULE MWretrModule
  !
  ! Module for executing MW retrieval process
  !
  ! Subroutine:
  !
  !     MWretrInit
  !     MWretr
  !     MWretrDestroy
  !
  ! Derived data type:
  !
  !     MWretrControl_t
  !
  ! USE:
  !
  !     constants
  !     ToolboxModule
  !     ControlStructure
  !     StateIndexModule
  !     IRMapInvert
  !     oss_mw_module
  !     ChkValid
  !     MapInvert
  !     InvertModule
  !     MWRTmodule
  !     MWobsStructure
  !     MeasErrStructure
  !     AncillaryStructure
  !     ChannelSelectionModule
  !     VertCoord
  !     BackgroundModule
  !
  ! initial, 03/16/2016, yhe@aer.com
  ! revised using background object, 06/15/2017, yhe@aer.com
  !
  USE constants, Only: &
       MISSING_REAL

  USE ToolboxModule, Only: &
       getUnit

  USE ControlStructure, Only: &
       GenControl_t, &
       RetrResult_t, &
       StateSpaceDefin_t

  USE StateIndexModule, Only: &
       initLengths, &
       genIndices, &
       getVectorLength, &
       StateIndex_t, &
       maxMol, &
       whereH2O, &
       getNmol

  USE MapInvert, Only: &
       map_retr2geo, &
       map_geo2retr, &
       set_MW_Invert

  USE oss_mw_module, Only: &
       ossdrv_mw

  USE ChkValid, Only: &
       getChiSq, &
       chkges

  USE InvertModule, Only: &
       invrt2, &
       DBL_INVT

  USE MWRTmodule, Only: &
       MWRTcontrol_t

  USE MWobsStructure, Only: &
       MWdefin_t, &
       MWob_t

  USE MeasErrStructure, Only: &
       MeasErr_t

  USE AncillaryStructure, Only: &
       AncillaryDefin_t, &
       AncillaryData_t

  USE ChannelSelectionModule, Only: &
       ChannelSet_t

  USE VertCoord, Only: &
       mxCoordTyp, &
       Pcoord

  USE BackgroundModule

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: &
       MWretrInit, &
       MWretr, &
       MWretrDestroy, &
       MWretrControl_t

  integer, parameter :: mxLength=256
  integer, parameter :: mxiter=15
  real, parameter :: deltaFrq=0.01

  TYPE MWretrControl_t
     TYPE(StateIndex_t) :: nRetr
     integer :: maxIter
     real :: drad
     real :: alpha1slope
     real :: alpha1int
     real :: chiSqConv
     real :: chiSqRatioConv
  END TYPE MWretrControl_t

  TYPE(StateIndex_t) :: NR_orig
  TYPE(StateIndex_t) :: IR_orig
  TYPE(StateIndex_t) :: NR
  TYPE(StateIndex_t) :: IR
  TYPE(StateIndex_t) :: NG
  TYPE(StateIndex_t) :: IG
  TYPE(background_t) :: bkgObj
  integer :: ncldLiqOrig
  integer :: ncldIceOrig
  integer :: nEmMw_orig
  integer :: nPar
  integer :: nParAtm
  integer :: nParSfcMW
  integer :: nParMax
  integer :: iH2O
  integer :: nChan
  integer :: nParG
  integer :: nParGatm
  integer :: nParGsfcMW
  integer :: nLev
  integer :: IEmMwG
  integer :: nEmMwG
  integer :: nMol
  integer, dimension(maxMol) :: MolTran
  logical, dimension(:), allocatable :: kchan
  real, dimension(:), allocatable :: pRef
  real, dimension(:), allocatable :: press
  real, dimension(:), allocatable :: Y,Sy,delY,xkEmMw,emMW
  real, dimension(:), allocatable :: xGes,xGesNew
  real, dimension(:), allocatable :: background,xGesG
  real, dimension(:), allocatable :: dradChan
  real, dimension(:), allocatable :: rad
  real, dimension(:,:), allocatable :: xktG,xkt
  real, dimension(:,:), allocatable :: covRetr
  real, dimension(:,:), allocatable :: ut_cov_net !internal (nParGatm,nParGatm)
  real, dimension(:,:), allocatable :: umtx,vmtx
  real, dimension(:,:), allocatable :: clw_cov,ciw_cov
  real(kind=DBL_INVT), dimension(:,:), allocatable :: dSx
  logical :: extAtmosOn
  logical :: logSize
  logical :: debug=.false.  !temporary solution to the input to chkges, 7/11/2016
  logical :: dbg = .true.
  character (len=mxCoordTyp) :: vCoordTyp

CONTAINS

  subroutine MWretrInit( &
       genControl, &
       MWretrControl, &
       MWRTcontrol, &
       MWdefin, &
       stateSpaceDefin &
       )
    TYPE(GenControl_t), intent(in) :: genControl
    TYPE(MWretrControl_t), intent(in) :: MWretrControl
    TYPE(MWRTcontrol_t), intent(in) :: MWRTcontrol
    TYPE(MWdefin_t), intent(inout) :: MWdefin
    TYPE(StateSpaceDefin_t), intent(inout) :: stateSpaceDefin

    !Local
    integer kk
    character (len=*), parameter :: procName = &
                           ' [MWretrModule::MWretrInit]: '
    if (dbg) print *, trim(procName)//'starting ...'

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

    nParGatm = stateSpaceDefin%nParGatm

    CALL bkgObj%bkgInitSfcMW( &
         stateSpaceDefin%nChanMW, &
         stateSpaceDefin%frqMW, &
         stateSpaceDefin%polMW, &
         stateSpaceDefin%IG%emMW, &
         genControl%backgroundFiles(2))

    stateSpaceDefin%NG%emMW = stateSpaceDefin%nChanMW
    stateSpaceDefin%IG = genIndices(stateSpaceDefin%NG)
    stateSpaceDefin%nParG = getVectorLength(stateSpaceDefin%NG)

    if (MWretrControl%maxiter > mxiter) then
       print *, trim(procName)//'Error: MWretrControl%maxiter > mxiter ', &
            MWretrControl%maxiter, mxiter
       call errorHalt(1)
    endif

    if (stateSpaceDefin%nChanMW /= MWdefin%nChan) then
       print *, trim(procName)//' error - nChanMW inconsistent (Bkg/Rad) ', &
            stateSpaceDefin%nChanMW, MWdefin%nChan
       call errorHalt(1)
    endif
    if (any(abs(stateSpaceDefin%frqMW(1:stateSpaceDefin%nChanMW) - &
         MWdefin%frq(1:stateSpaceDefin%nChanMW)) > deltaFrq)) then
       print *, trim(procName)//' error - MW Frequencies inconsistent (Bkg/Rad)'
       call errorHalt(1)
    endif
    if (any(stateSpaceDefin%polMW(1:stateSpaceDefin%nChanMW) /=  &
         MWdefin%pol(1:stateSpaceDefin%nChanMW))) then
       print *, trim(procName)//' error - Polarizations inconsistent (Bkg/Rad)'
       call errorHalt(1)
    endif

    !----Set retrieval-space pointers
    call setXPtrR( &
         MWretrControl%nRetr%Tskin, &
         MWretrControl%nRetr%temp,  &
         MWretrControl%nRetr%mol,   &
         MWretrControl%nRetr%cldLiq, &
         MWretrControl%nRetr%cldIce, &
         MWRTcontrol%molID, &
         MWRTcontrol%nMol, &
         stateSpaceDefin%NG, &
         IR, NR)

    !----Verify that the selected number of variables to retrieve are
    !----compatible with the transformation matrices
    !
    nParAtm = getVectorLength(NR)
    NR%emMW = MWretrControl%nRetr%emMW
    IR = genIndices(NR)
    call bkgObj%checkTransformConform(NR)

    !--- Save the original retrieval-space pointers: *may not be needed - YHE
    !----------------------------------------------------------
    !  Save original versions of variables that are dynamically
    !  changed in the ProfLoop via call to IRsetRetrVec()
    !  to fall back if needed.
    !----------------------------------------------------------
    NG = stateSpaceDefin%NG
    IG = genIndices(NG)
    ncldLiqOrig = MWretrControl%nRetr%cldLiq
    ncldIceOrig = MWretrControl%nRetr%cldIce
    IR_orig = IR
    NR_orig = NR
    nParSfcMW = NR_orig%emMW
    nParMax = getVectorLength(NR) !%emMW included
    nChan = stateSpaceDefin%nchanMW
    nParG = getVectorLength(NG)
    call bkgObj%bkgAllocate(nParG)
    nLev = stateSpaceDefin%nLev
    IEmMwG = stateSpaceDefin%IG%emMW
    nEmMwG = nChan
    nParGsfcMW = nChan
    iH2O = whereH2O(stateSpaceDefin%MolID)
    nMol = getNmol(stateSpaceDefin%MolID)
    MolTran = stateSpaceDefin%MolTran
    logSize = stateSpaceDefin%logSize
    vCoordTyp = trim(stateSpaceDefin%vCoordTyp)

    if (.not. allocated(pRef)) allocate(pRef(nLev))
    pRef = stateSpaceDefin%pRef
    if (.not. allocated(press)) allocate(press(nLev))
    if (.not. allocated(kchan)) allocate(kchan(nChan))
    if (.not. allocated(rad)) allocate(rad(nChan))
    if (.not. allocated(Sy)) allocate(Sy(nChan))
    if (.not. allocated(delY)) allocate(delY(nChan))
    if (.not. allocated(xGesNew)) allocate(xGesNew(nParMax))
    if (.not. allocated(xktG)) allocate(xktG(nParG,nChan))
    if (.not. allocated(xkt)) allocate(xkt(nParMax,nChan))
    if (.not. allocated(background)) allocate(background(nParG))
    if (.not. allocated(ut_cov_net)) allocate(ut_cov_net(nParG,nParG))
    if (.not. allocated(covRetr)) allocate(covRetr(nParMax,nParMax))
    if (.not. allocated(dSx)) allocate(dSx(nParMax,nParMax))
    if (.not. allocated(umtx)) allocate(umtx(nParMax,nParG))
    if (.not. allocated(vmtx)) allocate(vmtx(nParG,nParMax))
    if (.not. allocated(xGesG)) allocate(xGesG(nParG))
    if (.not. allocated(xGes)) allocate(xGes(nParMax))
    if (.not. allocated(xkEmMw)) allocate(xkEmMw(nChan))
    if (.not. allocated(emMW)) allocate(emMW(nChan))
    if (.not. allocated(clw_cov)) allocate(clw_cov(nCldLiqOrig,nCldLiqOrig))
    if (.not. allocated(ciw_cov)) allocate(ciw_cov(nCldIceOrig,nCldIceOrig))

    !--------------------------------------------------------
    ! Allocate and initialize dradChan with drad
    !--------------------------------------------------------
    if (.not. allocated(dradChan)) allocate(dradChan(nChan))
    ! In radiance model, drad must be converted from dTB to dradiance
    IF (MWdefin%radOrTb == 0) THEN
       dradChan = MWretrControl%drad
    ELSE
       dradChan(1:nChan) = MWretrControl%drad/MWdefin%planckAlpha(1:nChan)
    ENDIF

    extAtmosOn=.FALSE.
    if (ANY((/genControl%externDataFlags%temp, &
         genControl%externDataFlags%mol, &
         genControl%externDataFlags%pSfc, &
         genControl%externDataFlags%wind, &
         genControl%externDataFlags%cloudLiq, &
         genControl%externDataFlags%cloudIce/))) extAtmosOn=.TRUE.

    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine MWretrInit

  subroutine MWretr(&
       genControl, &
       MWretrControl, &
       MWdefin, &
       MWob, &
       MWmeasErr, &
       ancillaryDefin, &
       ancillaryOb, &
       channelSetMW, &
       atmosClass, &
       cloudClass, &
       surfaceClassMW, &
       isLandMW, &
       retrPrior, &
       retrResult)

    TYPE(GenControl_t), intent(in) :: genControl
    TYPE(MWretrControl_t), intent(in) :: MWretrControl
    TYPE(MWdefin_t), intent(in) :: MWdefin
    TYPE(MWOb_t), intent(in) :: MWob
    TYPE(MeasErr_t) :: MWmeasErr
    TYPE(AncillaryDefin_t), intent(in) :: ancillaryDefin
    TYPE(AncillaryData_t), intent(in) :: ancillaryOb
    TYPE(ChannelSet_t), intent(in) :: ChannelSetMW
    TYPE(RetrResult_t), intent(in) :: retrPrior
    integer, intent(in) :: atmosClass
    integer, dimension(:), intent(in) :: cloudClass
    integer, intent(in) :: surfaceClassMW
    logical, intent(in) :: isLandMW
    TYPE(RetrResult_t), intent(inout) :: retrResult

    !Local
    integer :: nchanOn
    integer :: iter
    integer :: ifail
    integer :: U_algconfig
    real :: Tsfc
    real :: chiSqMw
    real :: rmsMw
    real, dimension(mxiter) :: rms_mw,chisq_mw
    logical :: xOutRangeFlag
    character (len=mxLength) :: procName
    integer :: ix
    integer :: nCopy
    TYPE(retFlags_t) :: retrievalFlags = RET_FLAGS_ALL_ON
    integer :: iFOV
    integer :: j
    ! obsLevel the level at which the instrumnet is located
    ! TOA corresponds to obsLevel=1
    integer, parameter :: obsLevel=1

    procName = ' [MWretrModule::MWretr]: '
    if (dbg) print *, trim(procName)//'starting ...'
    print *, trim(procName)//' warning - verification on dimensions needed before assignments...'

    ! Disable retrieval of cloud parameters per cloud mask
    ! assumption: use the first element of cloudClass for the cloudModeTest
    iFOV = 1
    if (genControl%imgDataOn) &
         call bkgObj%setCloudMode(cloudClass(iFOV),retrievalFlags,genControl%imgDataOn)

    ! Retrieve parameters based on the values of the flags in retrievalFlags
    ! Currently it is all ON.  Need to change that.
    CALL bkgObj%setRetrVec(retrievalFlags,NR_orig,nPar,NR,IR)

    umtx = 0.0
    vmtx = 0.0
    ut_cov_net = 0.0
    call bkgObj%setBkgClimAtmos(atmosClass,ancillaryOb%pSfc,vCoordTyp,press, &
         umtx(1:nParAtm,1:nParGAtm),vmtx(1:nParGatm,1:nParAtm))
    call bkgObj%setBkgClimSfcMW(surfaceClassMW, &
         umtx(nParAtm+1:nPar,nParGAtm+1:nParG),vmtx(nParGAtm+1:nParG,nParAtm+1:nPar))
    if (extAtmosOn) call bkgObj%setBkgExt(ancillaryOb%stateV)
    call bkgObj%combineBkgAtmos(ancillaryDefin%IG,ancillaryDefin%NG,ancillaryOb%pSfc,isLandMW, &
         MolTran,vCoordTyp,press,background(1:nParG),ut_cov_net(1:nParG,1:nParG))
    call bkgObj%combineBkgSfcMW(isLandMW,background(nParGatm+1:nParGatm+nParGsfcMW), &
         ut_cov_net(nParGatm+1:nParGatm+nParGsfcMW,nParGatm+1:nParGatm+nParGsfcMW))
    call bkgObj%overrideBkg(background,ut_cov_net,isLandMW)
    call bkgObj%transformBkg(ut_cov_net,covRetr(1:nParMax,1:nParMax), &
         umtx(1:nPar,1:nParG))

    rms_mw   = MISSING_REAL !local vector, needed outside?
    chisq_mw = MISSING_REAL

    !do this if retrPrior%nParG < stateSpaceDefin%nParG
    if (retrPrior%nParG < nParG) then
       xGes=0.
       if(NR%cldLiq > 0)xGes(IR%cldliq+2)=genControl%LWP1stGuess
       if(NR%cldIce > 0)xGes(IR%cldIce+2)=genControl%IWP1stGuess
       !nParG,pRef,nLev, etc are taken from module copy of stateSpaceDefin set by MWretrInit
       call map_retr2geo(xGes(1:nPar),background,xGesG,emMW,TSfc, &
            vmtx(1:nParG,1:nPar),nParG,nPar,IG,NG,vCoordTyp, &
            pRef,nLev,nChan,IEmMwG,iH2O,nMol,MolTran, &
            xOutRangeFlag,logSize)
       if (xOutRangeFlag) &
            print *,trim(procName)//' warning - out-of-range background may affect convergence'
    end if
    nCopy=min(nParG,retrPrior%nParG)
    if (retrPrior%nParG > 0) xGesG(1:nCopy) = retrPrior%stateV(1:nCopy)
    !"and also transfer MW emis if retrPrior sizes /=0"
    if (retrPrior%nEmMW > 0) then
       emMW(1:retrPrior%nEmMW) = retrPrior%emMW(1:retrPrior%nEmMW)
    end if

    if (NR%cldLiq > 0) then
       clw_cov(1:NR%cldLiq,1:NR%cldLiq) = &
            covRetr(IR%cldLiq:IR%cldLiq+NR%CldLiq-1,IR%cldLiq:IR%cldLiq+NR%CldLiq-1)
    endif
    if (NR%cldIce > 0) then
       ciw_cov(1:NR%cldIce,1:NR%cldIce) = &
            covRetr(IR%cldIce:IR%cldIce+NR%CldIce-1,IR%cldIce:IR%cldIce+NR%CldIce-1)
    endif

    xktG=0.  ! To initialize portions not filled later

    if (channelSetMW%nMW<0) then
      print *, 'err: ', procName, 'Expected MW channel Selection for ', &
         'for MW retrieval only not supplied'
      call exit(1)
    end if

    kchan = channelSetMW%MW /= 0
   ! ossdrv_mw has a flag lambertion controling the reflection from the surface
   ! if omitted it set to .false. specular reflection


    if (trim(vCoordTyp) == Pcoord) then
      CALL ossdrv_mw(xGesG,emMW,MWob%EIA,obsLevel,rad,xktG,xkEmMw,kchan,ancillaryOb%landfrac,&
            MWob%lat, IG%temp,IG%Tskin,IG%pSfc,IG%mol(1:iH2O), ICldLiqG_in=IG%cldLiq)
    else
      CALL ossdrv_mw(xGesG,emMW,MWob%EIA,obsLevel,rad,xktG,xkEmMw,kchan,ancillaryOb%landfrac,&
          MWob%lat, IG%temp,IG%Tskin, IG%pSfc, IG%mol(1:iH2O),&
          ICldLiqG_in=IG%cldLiq, puser=press)
    end if
    call getChiSq(MWob%rad,rad,MWmeasErr%device,chiSqMw,rmsMw,kchan,nChan,nChanOn)
    iter = 0
    if (dbg) print "(' MW iter.   RMS(',f11.7,')      NormChiSq')",MWmeasErr%rms
    if (dbg) print '(i9,f19.7,f14.3)',iter,rmsMw,chiSqMw

    !---Start mw iteration:
    MwIterLoop: do iter=1,MWretrControl%maxIter
       CALL set_MW_Invert(xGesG,MWob%rad,rad,xktG,xkEmMw,MWmeasErr%device, &
            xkt(1:nPar,1:nChanOn),Sy,delY,nChanOn,iter,kchan, &
            vmtx(1:nParG,1:nPar),nChan,dradChan,MWretrControl%alpha1slope, &
            MWretrControl%alpha1int,nParG,nPar,IG,NG,IEmMwG,nMol,MolTran,vCoordTyp, &
            pRef,logSize)

       call invrt2(xkt(1:nPar,1:nChanOn),Sy(1:nChanOn),covRetr(1:nPar,1:nPar), &
            delY(1:nChanOn),xGes(1:nPar),xGesNew(1:nPar),nPar,nChanOn)

       xGes(1:nPar) = xGesNew(1:nPar)

       call map_retr2geo(xGesNew(1:nPar),background,xGesG,emMW,TSfc, &
            vmtx(1:nParG,1:nPar),nParG,nPar,IG,NG,vCoordTyp,pref,nlev, &
            nChan,IEmMwG,iH2O,nMol,MolTran, &
            xOutRangeFlag,logSize)

       if (xOutRangeFlag) &
            print *,'warning:[retr] x out of range may affect '// &
            'convergence in MW loop'

       U_algconfig = getUnit()
       call chkges(U_algconfig,genControl%backgroundConfig,xGesG,ifail, &
            covRetr(1:nPar,1:nPar),clw_cov,ciw_cov,nLev,press, &
            IG,NG,IR,NR,iH2O,debug)  !put-off: remove debug (make it parameter in ChkValid)
       if (trim(vCoordTyp) == Pcoord) then
          CALL ossdrv_mw(xGesG,emMW,MWob%EIA,obsLevel,rad,xktG,xkEmMw,kchan,ancillaryOb%landfrac,&
               MWob%lat, IG%temp,IG%Tskin,IG%pSfc,IG%mol(1:iH2O), ICldLiqG_in=IG%cldLiq)
       else
          CALL ossdrv_mw(xGesG,emMW,MWob%EIA,obsLevel,rad,xktG,xkEmMw,kchan,ancillaryOb%landfrac,&
               MWob%lat, IG%temp,IG%Tskin, IG%pSfc, IG%mol(1:iH2O),&
               ICldLiqG_in=IG%cldLiq, puser=press)
       end if
       call getChiSq(MWob%rad,rad,MWmeasErr%device,chiSqMw,rmsMw,kchan,nChan,nChanOn)

       rms_mw(iter)   = rmsMw
       chisq_mw(iter) = chiSqMw

       if (dbg) print '(i9,f19.7,f14.3)',iter,rmsMw,chiSqMw

       if (chiSqMw < MWretrControl%chiSqConv) exit MwIterLoop

       CALL map_geo2retr(xGesG,background,Tsfc,xGes(1:nPar), &
            umtx(1:nPar,1:nParG),nParG,nPar,IG,NG, &
            nMol,MolTran,vCoordTyp,logSize)

    enddo MwIterLoop

    iter=min(iter,MWretrControl%maxIter)
    retrResult%nParG = nParG
    retrResult%nEmMW = nEmMwG
    retrResult%stateV(1:nParG) = xGesG(1:nParG)
    retrResult%emMW(1:nEmMwG) = emMW(1:nEmMwG)
    retrResult%nIter = iter
    retrResult%rms = rms_mw(iter)
    retrResult%chiSq = chisq_mw(iter)
    retrResult%IG = IG
    retrResult%NG = NG
    retrResult%press = press
    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine MWretr

  subroutine MWretrDestroy()
    !Local
    character (len=mxLength) :: procName

    procName = ' [MWretrModule::MWretrDestroy]: '
    if (dbg) print *, trim(procName)//' starting ...'
    if (allocated(pRef)) deallocate(pRef)
    if (allocated(press)) deallocate(press)
    if (allocated(kchan)) deallocate(kchan)
    if (allocated(rad)) deallocate(rad)
    if (allocated(Sy)) deallocate(Sy)
    if (allocated(delY)) deallocate(delY)
    if (allocated(xGesNew)) deallocate(xGesNew)
    if (allocated(xktG)) deallocate(xktG)
    if (allocated(xkt)) deallocate(xkt)
    if (allocated(background)) deallocate(background)
    if (allocated(covRetr)) deallocate(covRetr)
    if (allocated(dSx)) deallocate(dSx)
    if (allocated(umtx)) deallocate(umtx)
    if (allocated(vmtx)) deallocate(vmtx)
    if (allocated(xGesG)) deallocate(xGesG)
    if (allocated(xGes)) deallocate(xGes)
    if (allocated(xkEmMw)) deallocate(xkEmMw)
    if (allocated(emMW)) deallocate(emMW)
    if (allocated(dradChan)) deallocate(dradChan)
    call bkgObj%bkgDestroy()
    if (dbg) print *, trim(procName)//' ending ...'
  end subroutine MWretrDestroy

END MODULE MWretrModule
