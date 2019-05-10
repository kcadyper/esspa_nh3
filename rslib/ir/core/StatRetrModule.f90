MODULE StatRetrModule
  !
  ! Module for executing statistical retrieval process
  !
  ! Subroutine:
  !
  !     StatRetrInit
  !     StatRetr
  !     StatRetrDestroy
  !
  ! Derived data types:
  !
  !    StatRetrControl_t
  !
  ! USE:
  !
  !    StateIndexModule
  !    AncillaryStructure
  !    MWobsStructure
  !    IRobsStructure
  !    ControlStructure
  !    RegrPrepModule
  !    InvertModule
  !
  ! yhe@aer.com, 03/16/2016
  !
  USE StateIndexModule, Only: &
       StateIndex_t, &
       maxMol, &
       idH2O, &
       whereH2O, &
       genIndices, &
       getVectorLength

  USE AncillaryStructure, Only: &
       AncillaryFOR_t

  USE MWobsStructure, Only: &
       MWdefin_t, &
       MWOb_t

  USE IRobsStructure, Only: &
       IRdefin_t, &
       IRFOR_t

  USE ControlStructure, Only: &
       GenControl_t, &
       StateSpaceDefin_t, &
       RetrResult_t

  USE RegrPrepModule, Only: &
       getRegrParam, &
       setFOVRegr

  USE InvertModule, Only: &
       invrtRegr

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: &
       statRetrInit, &
       statRetr, &
       statRetrDestroy, &
       statRetrControl_t

  integer, parameter :: mxLength = 256
  integer, parameter :: nmolRegr = 1
  integer, parameter :: mxRdEof  = 200

  TYPE StatRetrControl_t
     character(len=mxLength) :: statRetrCoefFile
  END TYPE StatRetrControl_t

  TYPE(StateIndex_t) :: NGregr
  TYPE(StateIndex_t) :: IGregr
  real, dimension(:), pointer :: prefRegr
  integer :: nChanRegr
  integer :: nChanMW
  integer :: nChanIR
  integer :: nChan
  integer :: iH2O
  integer :: nParGRegr
  integer :: nParRegr
  integer :: nEofRegr
  integer, dimension(nmolRegr) :: molIDRegr=(/idH2O/) !different from statRetrControl%molID?
  real, dimension(:), allocatable   :: xbackG_regr
  real, dimension(:), allocatable   :: xregr
  real, dimension(:), allocatable   :: tbEquiv
  real, dimension(:), allocatable   :: emMW
  real, dimension(:), allocatable   :: emRfIR
  real, dimension(:), allocatable   :: rback_regr
  real, dimension(:), allocatable   :: wght_regr
  real, dimension(:,:), allocatable :: regr_coef
  real, dimension(:,:), allocatable :: umtx_regr
  real, dimension(:,:), allocatable :: umtx_rad
  logical :: dbg = .true.

CONTAINS

  subroutine statRetrInit( & !- not tested, pending for regression inputs
       genControl, &
       statRetrControl, &
       MWdefin, &
       IRdefin, &
       stateSpaceDefin &
       )
    TYPE(GenControl_t), intent(in) :: genControl
    TYPE(StatRetrControl_t), intent(in) :: statRetrControl
    TYPE(MWdefin_t), intent(in) :: MWdefin
    TYPE(IRdefin_t), intent(in) :: IRdefin
    TYPE(StateSpaceDefin_t), intent(inout) :: stateSpaceDefin
    
    !Local, although when statRetr is implemented, some of the
    !variables may need to be relocated to the module level
    character (len=mxLength) :: procName
    integer :: imol

    procName = ' [StatRetrModule::StatRetrInit]: '
    if (dbg) print *, trim(procName)//' starting ...'

    print *, trim(procName)//' returning before testing...'
    call getRegrParam( & !U/F_regr in RegrPrepModule::getRegrParam needs to be revised
         stateSpaceDefin%nParGatm, &
         nParRegr, &
         nChanRegr, &
         nEofRegr, &
         stateSpaceDefin%IG, &
         stateSpaceDefin%NG, &
         stateSpaceDefin%nLev, &
         prefRegr, &
         molIDRegr &
         )

    stateSpaceDefin%molID(1:nmolRegr) = molIDRegr(1:nmolRegr)
    if (.not. allocated(stateSpaceDefin%pRef)) allocate(stateSpaceDefin%pRef(stateSpaceDefin%nLev))
    stateSpaceDefin%pRef = prefRegr
    if (genControl%MWdataOn) then
       nChanMW = MWdefin%nChan
    else
       nChanMW = 0
    endif
    nChanIR = IRdefin%nChan
    nChan = nChanMW+nChanIR
    allocate(tbEquiv(nChan))
    iH2O=whereH2O(molIDregr)

    stateSpaceDefin%nChanMW = nChanMW
    stateSpaceDefin%NG%emMW = stateSpaceDefin%nChanMW
    stateSpaceDefin%IG = genIndices(stateSpaceDefin%NG)
    stateSpaceDefin%nParG = getVectorLength(stateSpaceDefin%NG)
    NGregr = stateSpaceDefin%NG
    IGregr = stateSpaceDefin%IG

    if (.not. allocated(emMW)) allocate(emMW(stateSpaceDefin%NG%emMW))
    if (.not. allocated(emRfIR)) allocate(emRfIR(stateSpaceDefin%NG%emRfIR))
    if (.not. allocated(xbackG_regr)) allocate(xbackG_regr(stateSpaceDefin%nParG))
    if (.not. allocated(rback_regr)) allocate(rback_regr(stateSpaceDefin%nParG))
    if (.not. allocated(xregr)) allocate(xregr(stateSpaceDefin%nParG))
    if (.not. allocated(wght_regr)) allocate(wght_regr(nChanRegr))
    if (.not. allocated(umtx_regr)) allocate(umtx_regr(stateSpaceDefin%nParG,NParRegr))
    if (.not. allocated(regr_coef)) allocate(regr_coef(nEofRegr,nParRegr)) 
    if (.not. allocated(umtx_rad)) allocate(umtx_rad(nChanRegr,nEofRegr)) 

    if (dbg) print *, trim(procName)//'ending ...'

  end subroutine StatRetrInit
  
  !
  ! Current statistical retrieval is one IR-FOR and one MW FOV based,
  ! which is evidenced by concatenating radiance of both MW and IR as
  ! one vector of tbEquiv, and the returned retrResult is on a single
  ! FOR. - Y. HE, 6/10/2016.
  ! 
  subroutine statRetr( & !- not tested, yet, pending for regression inputs
       genControl, &
       statRetrControl, &
       MWdefin, &
       MWob, &
       IRdefin, &
       IRFOR, &
       ancillaryFOR, &
       retrResult &
       )
    TYPE(GenControl_t), intent(in) :: genControl
    TYPE(StatRetrControl_t), intent(in) :: statRetrControl
    TYPE(MWdefin_t), intent(in) :: MWdefin
    TYPE(MWOb_t), intent(in) :: MWob
    TYPE(IRdefin_t), intent(in) :: IRdefin
    TYPE(IRFOR_t), intent(in) :: IRFOR
    TYPE(AncillaryFOR_t), intent(in) :: AncillaryFOR
    TYPE(RetrResult_t), intent(inout) :: retrResult
    !Local
    integer :: iClassAtmRegr
    integer :: i, ii
    !- these local variables ought to be from the input arguments
    integer :: iFOV = 1
    integer :: iFOR = 1
    integer :: igeo = 1
    character (len=mxLength) :: procName

    procName = ' [StatRetrModule::StatRetr]: '
    if (dbg) print *, trim(procName)//'starting ...'
    if (genControl%MWdataOn) then
       IF (MWdefin%radOrTb == 0) THEN
          tbEquiv(1:nChanMW) = MWob%rad(1:nChanMW)
       ELSE               ! Convert from radiance to radiance temperature
          tbEquiv(1:nChanMW) = &
               MWob%rad(1:nChanMW)*MWdefin%planckAlpha(1:nChanMW)+ &
               MWdefin%planckBeta(1:nChanMW)
       ENDIF
    end if
    tbEquiv(nChanMW+1:nChanMW+IRdefin%nChan) = &
         IRFOR%ob(iFOV)%rad(nChanMW+1:nChanMW+IRdefin%nChan)  !for regr

    call classatmRegr(tbEquiv,IRFOR%lat,IRFOR%lon,IRdefin%startTime, &
         IRFOR%EIA,iClassAtmRegr)
    !ancillaryFOR%pSfc or iFOR is not used within the callee
    call setFOVRegr(iClassAtmRegr,ancillaryFOR%pSfc,iFOR,&
         xbackG_regr,rback_regr,wght_regr,&
         umtx_regr,umtx_rad,regr_coef)
    call invrtRegr(tbEquiv(1:nChanRegr),&
         xbackG_regr(1:nParGRegr),&
         rback_regr(1:nParGRegr),&
         wght_regr(1:nChanRegr),&
         umtx_regr(1:nParGRegr,1:nParRegr),&
         umtx_rad(1:nChanRegr,1:nEofRegr),&
         regr_coef(1:nEofRegr,1:nParRegr),&
         nChanRegr,nParGRegr,&
         nEofRegr,nParRegr,xregr)
    xregr(IGRegr%pSfc) = ancillaryFOR%pSfc
    xregr(IGRegr%mol(ih2o):IGRegr%mol(ih2o)+NGRegr%mol(ih2o)-1) = &
         exp(xregr(IGRegr%mol(ih2o):IGRegr%mol(ih2o)+NGRegr%mol(ih2o)-1))

    retrResult%nParG = nParGRegr
    retrResult%stateV(1:nParGRegr) = xregr

    if (NGregr%emMW > 0) then
       emMW(1:NGregr%emMW) = xregr(IGregr%emMw:IGregr%emMW+NGregr%emMW-1)
       retrResult%nEmMW = NGregr%emMW
       retrResult%emMW(1:NGregr%emMW) = emMW
    end if

    if (NGregr%emRfIR > 0) then
       do i=1,NGregr%emRfIR
          ii=2*(i-1)+1
          emRfIR(ii)=xregr(IGregr%emRfIR+i-1)
          emRfIR(ii+1)=xregr(IGregr%emRfIR+NGregr%emRfIR+i-1)
       enddo
       retrResult%nEmRfIR = NGregr%emRfIR
       retrResult%emRfIR = emRfIR
    end if

    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine statRetr

  subroutine statRetrDestroy()
    !Local
    character (len=mxLength) :: procName
    procName = ' [StatRetrModule::StatRetrDestroy]: '
    if (dbg) print *, trim(procName)//'starting ...'
    if (allocated(xbackG_regr)) deallocate(xbackG_regr)
    if (allocated(rback_regr)) deallocate(rback_regr)
    if (allocated(xregr)) deallocate(xregr)
    if (allocated(wght_regr)) deallocate(wght_regr)
    if (allocated(umtx_regr)) deallocate(umtx_regr)
    if (allocated(regr_coef)) deallocate(regr_coef)
    if (allocated(umtx_rad)) deallocate(umtx_rad)
    if (allocated(emMW)) deallocate(emMW)
    if (allocated(emRfIR)) deallocate(emRfIR)
    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine StatRetrDestroy
       
END MODULE StatRetrModule
