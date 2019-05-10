MODULE AtmosClassModule
  !
  ! Module for atmospheric classification
  !
  ! Subroutines:
  !
  !     atmosClassInit
  !     loadAtmosClassData
  !     putAtmosClassData
  !     atmosClassify
  !     atmosClassDestroy
  !     classAtmNN
  !     rdthr
  !
  ! Derived data types:
  !
  !     AtmosClassData_t
  !
  ! USE:
  !
  !     ControlStructure
  !
  ! yhe@aer.com, 03/16/2016
  !
  USE constants, ONLY: &
       deg2rad, &
       MISSING_INT

  USE ToolboxModule, Only: &
       getUnit

  USE ControlStructure, Only: &
       GenControl_t

  USE MWobsStructure, Only: &
       MWdefin_t, &
       MWob_t

  USE IRobsStructure, Only: &
       IRdefin_t, &
       IRob_t

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: &
       atmosClassInit, &
       loadAtmosClassData, &
       putAtmosClassData, &
       atmosClassify, &
       atmosClassDestroy, &
       AtmosClassData_t, &
       MWmode, &
       IRmode, &
       MWIRmode

  integer, parameter :: dbl=SELECTED_REAL_KIND(12)
  real(kind=dbl), dimension(:),   allocatable :: profile
  real(kind=dbl), dimension(:),   allocatable :: a1
  integer :: iFstChan

  TYPE AtmosClassData_t
     logical :: disable
     logical :: channelMissingIsFatal=.false.
     integer :: nIn = 2       ! number of frequencies used in classification
     integer :: nThresh        ! number of thresholds used for classification 
     integer :: nAtmClasses  ! Number of atmosphere classes
     real, dimension(:),  allocatable :: threshVals
     integer, dimension(2) :: idIn
     real,    dimension(2) :: frqnn
     character(len=12), dimension(:),  allocatable :: idAtmClasses
  END TYPE AtmosClassData_t

  integer, parameter :: mxLength = 256
  integer, parameter :: MWmode=1
  integer, parameter :: IRmode=2
  integer, parameter :: MWIRmode=3
  integer, parameter :: missingID=-999
  integer :: nChan
  integer :: mode
  real,    dimension(:), allocatable :: rad
  real,    dimension(:), allocatable :: measErr
  integer, dimension(:), allocatable :: QC
  logical :: dbg = .false.

  TYPE(AtmosClassData_t), save :: classData
  logical, save :: dataReady = .FALSE.
  logical, save :: channelMissing

CONTAINS

  subroutine atmosClassInit(genControl,MWdefin,IRdefin,atmosClassMode)
    TYPE(GenControl_t), intent(in) :: genControl
    TYPE(MWdefin_t), intent(in) :: MWdefin
    TYPE(IRdefin_t), intent(in) :: IRdefin
    integer,         intent(in) :: atmosClassMode
    !Local
    real, parameter :: diffThreshold = 1.e-5
    character (len=mxLength) :: procName
    integer :: istatus
    integer :: i,j
    logical :: chancaExists
    real    :: reldif
    real,    dimension(:), allocatable :: frqall

    procName = '[AtmosClassModule::atmosClassInit]:'
    if (dbg) print *, trim(procName)//' starting ...'

    if (.not. dataReady) call loadAtmosClassData(genControl%atmosClassConfig)

    mode = atmosClassMode

    select case (mode)
    case (MWIRmode)
       nChan = MWdefin%nChan + IRdefin%nChan
    case (IRmode)
       nChan = IRdefin%nChan
    case (MWmode)
       nChan = MWdefin%nChan
    end select

    if (allocated(rad)) deallocate(rad)
    if (allocated(measErr)) deallocate(measErr)
    if (allocated(QC)) deallocate(QC)
    if (allocated(frqall)) deallocate(frqall)
    allocate(rad(nChan))
    allocate(measErr(nChan))
    allocate(QC(nChan))

    if (classData%disable) return

    allocate(frqall(nChan))

    select case (mode)
    case (MWIRmode)
       frqall(1:MWdefin%nChan) = MWdefin%frq
       frqall(MWdefin%nChan+1:nChan) = IRdefin%wvn
    case (IRmode)
       frqall(1:nChan) = IRdefin%wvn
    case (MWmode)
       frqall(1:nChan) = MWdefin%frq
    end select


  ! Check center frequencies for match
  ! with each channel used in classification

    iFstChan=1
    channelMissing=.false.
    do i=1,classData%nIn
      chancaExists=.FALSE.
      if (classData%idIn(i) == missingID) then ! Need to ID channel by frequency match
         do j=1,nChan
            reldif=abs(frqall(j)-classData%frqnn(i))/classData%frqnn(i)
            if (reldif < diffThreshold) then
               classData%idIn(i)=j
               chancaExists=.TRUE.
            endif
         enddo
      elseif (classData%idIn(i) > 0) then ! input is a channel identified by number
         reldif=abs(frqall(classData%idIn(i))-classData%frqnn(i))/ &
                classData%frqnn(i)
         if (reldif < diffThreshold) chancaExists=.TRUE.
      else
         iFstChan=iFstChan+1
      endif
      if ((i >= iFstChan) .and. .not. chancaExists) then
         channelMissing=.true.
         if (classData%channelMissingIsFatal) then
           print*,trim(procName)//'Error: ',&
                ' Found no match for Atm Class frequency',classData%frqnn(i)
           call exit(1)
         else
           print *
           print*,trim(procName)//'WARNING: ',&
                ' Found no match for Atm Class frequency',classData%frqnn(i)
           print *,'Cases will have class nAtmClasses+1'
           print *
         endif
      endif
    enddo

    deallocate(frqall)

    if (dbg) print *, trim(procName)//' ending ...'

  end subroutine atmosClassInit

  subroutine loadAtmosClassData(atmosClassConfigFile)
  ! This provides an option of a public method to load the static data
  ! at a level higher than atmosClassInit
    character(len=*), intent(in)    :: atmosClassConfigFile
    !Local
    character (len=mxLength) :: procName
    logical :: disable = .FALSE.
    logical :: channelMissingIsFatal = .FALSE.
    integer :: inu=12
    character(len=200) :: fileThresh
    TYPE(AtmosClassData_t) :: classDataIn
    namelist /classConfig/disable,channelMissingIsFatal,fileThresh

    procName = '[AtmosClassModule::loadAtmosClassData]:'
    if (dbg) print *, trim(procName)//' starting ...'

    open(inu,file=atmosClassConfigFile)
    read(inu,classConfig)
    close(inu)

    classDataIn%disable = disable
    classDataIn%channelMissingIsFatal = channelMissingIsFatal
    
    ! read data from file coefsClass and store it in classDataIn
    if (.not. disable) call rdthr(fileThresh,classDataIn)

    call putAtmosClassData(classDataIn)

    if (dbg) print *, trim(procName)//' ending ...'

  end subroutine loadAtmosClassData

  subroutine putAtmosClassData(classDataIn)
  ! This provides an option of a public method to insert the static data
    TYPE(AtmosClassData_t), intent(in) :: classDataIn
    !Local
    character (len=mxLength) :: procName

    procName = '[AtmosClassModule::putAtmosClassData]:'
    if (dbg) print *, trim(procName)//' starting ...'

    classData = classDataIn

    dataReady = .true.

    if (dbg) print *, trim(procName)//' ending ...'

  end subroutine putAtmosClassData

  subroutine atmosClassify( &
       MWdefin, &
       MWob, &
       IRdefin, &
       IRob, &
       surfaceClassMW, &
       surfaceClassIR, &
       atmosClass &   !integer, not structure
       )
    TYPE(MWdefin_t), intent(in) :: MWdefin
    TYPE(MWOb_t), intent(in) :: MWob
    TYPE(IRdefin_t), intent(in) :: IRdefin
    TYPE(IRob_t), intent(in) :: IRob
    integer, intent(in) :: surfaceClassMW
    integer, intent(in) :: surfaceClassIR
    integer, intent(out) :: atmosClass
    !Local
    character (len=mxLength) :: procName
    real, dimension(ClassData%nIn) :: bt_rad
    real*8, dimension(ClassData%nIn) :: freq
    real :: signal_est,snr,slope=12.54,offset=-0.991,bt_min=278.
    integer :: id
    
    

    procName = '[AtmosClassModule::atmosClassify]:'
    if (dbg) print *, trim(procName)//' starting ...'

!  set up for IR retrievals only
    select case (mode)
    case (MWIRmode)
       print *,'This code does not work with MW data'
       call errorHalt(1)
       rad(1:MWdefin%nChan) = MWob%rad
       QC(1:MWdefin%nChan) = MWob%QC
       rad(MWdefin%nChan+1:MWdefin%nChan+IRdefin%nChan) = IRob%rad
       QC(MWdefin%nChan+1:MWdefin%nChan+IRdefin%nChan) = IRob%QC
    case (IRmode) 
       rad(1:IRdefin%nChan) = IRob%rad
       QC(1:IRdefin%nChan) = IRob%QC
       freq = ClassData%frqnn
    case (MWmode)
       print *,'This code does not work with MW data'
       call errorHalt(1)
       rad(1:MWdefin%nChan) = MWob%rad
       QC(1:MWdefin%nChan) = MWob%QC
    end select

! Estimate NH3 BT signal 
    do id=1,ClassData%nIn
       bt_rad(id) = wnbrit(freq(id),rad(ClassData%idIn(id)))
    enddo
    signal_est = bt_rad(2)-bt_rad(1)
    snr = offset+slope*signal_est
    if (bt_rad(2).LT.bt_min) snr = 0.0    ! Cold surfaces or clouds make NH3 retrievals very unreliable

! Classify 
    if (abs(snr).GT.ClassData%threshVals(4)) atmosClass=1
    if (snr.GT.ClassData%threshVals(1).AND.snr.LT.ClassData%threshVals(2)) atmosClass=2
    if (snr.GT.ClassData%threshVals(3).AND.snr.LT.ClassData%threshVals(4)) atmosClass=2
    if (abs(snr).LT.ClassData%threshVals(3)) atmosClass=3


    if (dbg) print *, trim(procName)//' ending ...'

  end subroutine atmosClassify

  subroutine atmosClassDestroy()
    !Local
    character (len=mxLength) :: procName

    procName = '[AtmosClassModule::atmosClassDestroy]:'
    if (dbg) print *, trim(procName)//' starting ...'
    if (allocated(rad)) deallocate(rad)
    if (allocated(QC)) deallocate(QC)
    if (allocated(profile)) deallocate(profile)
    if (allocated(a1)) deallocate(a1)
    dataReady = .FALSE.
    if (dbg) print *, trim(procName)//' ending ...'
  end subroutine atmosClassDestroy

!----------------------------------------------------------------------------


  SUBROUTINE rdthr(dataFile,classDataIn)

    character(len=*),       intent(in)    :: dataFile
    TYPE(AtmosClassData_t), intent(inout) :: classDataIn
    INTEGER          :: iun,it
    INTEGER          :: nHeader
    INTEGER          :: nAtmClasses

    iun = getUnit()
    open(iun,file=dataFile)
    read(iun,*) classDataIn%frqnn
    classDataIn%idIn = missingID

    read(iun,*) classDataIn%nthresh
    if (allocated(classDataIn%threshVals)) deallocate(classDataIn%threshVals)
    allocate (classDataIn%threshVals(classDataIn%nthresh))
    do it=1,classDataIn%nthresh
       read(iun,*) classDataIn%threshVals(it)
    enddo
    read(iun,*) nAtmClasses
    if (allocated(classDataIn%idAtmClasses)) deallocate(classDataIn%idAtmClasses)
    allocate (classDataIn%idATmClasses(nAtmClasses))
    do it=1,nAtmClasses
       read(iun,*) classDataIn%idATmClasses(it)
    end  do
    close(iun)

    RETURN
  END SUBROUTINE rdthr
    
!----------------------------------------------------------------------------
     FUNCTION WNBRIT(VN,RAD)

!<f90Function>**********************************************************
!
! NAME:
!
!   WNBRIT
!
! PURPOSE:
!
!   Equivalent blackbody temperature of radiance (inverse Planck function).
!
! SYNTAX:
!
!   Results=WNBRIT(VN, RAD)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   VN      REAL*8  Wavenumber, units of cm^-1
!   RAD     REAL    Radiance
!
!   * OPTIONAL
!
! RETURN:
!
!     REAL    
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

! $ RADIANCE TO BRIGHTNESS TEMPERATURE (HMW)
! * 'NEW' PLANCK'S CONSTANT, VELOCITY OF LIGHT, BOLTZMANN'S CONSTANT
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL,   INTENT(IN)    :: RAD
      REAL*8, INTENT(IN)    :: VN
      REAL    :: WNBRIT
!
      PARAMETER (H = 6.626176D-27, C = 2.997925D+10, B = 1.380662D-16)
      PARAMETER (C1 = 2.D0*H*C*C)
      PARAMETER (C2 = H*C/B)
      TNB(X,Y,Z)=Y/DLOG(X/Z+1.D0)
!
      F1=C1*VN**3
      F2=C2*VN
      R=RAD
      TBB=TNB(F1,F2,R)
      WNBRIT=TBB
      RETURN
      END function wnbrit


END MODULE AtmosClassModule
