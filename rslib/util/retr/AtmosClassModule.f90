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
  !     rdwts
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
       classAtmNN, & ! This is public only for ClassAtmMod::classatm
       rdwts, & ! This is public only for ClassAtmMod::classatm
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
     integer :: nInNN        ! Actual number of inputs used in network
     integer :: nAtmClasses  ! Number of atmosphere classes
     integer :: ntrain       ! Actual number of profiles used in training
     integer, dimension(:),  allocatable :: idInNN
     real,    dimension(:),  allocatable :: frqnn
     real,    dimension(:),  allocatable :: scale
     ! inpwts: Array of input weights used in network layer 1
     ! bias1:  Array of biases used in network layer 2
     ! laywts: Array of layer weights used in network layer 2
     real(kind=dbl), dimension(:,:), allocatable :: inpwts
     real(kind=dbl), dimension(:),   allocatable :: bias1
     real(kind=dbl), dimension(:,:), allocatable :: laywts
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
    real, parameter :: diffThreshold = 1.e-4
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
    if (allocated(QC)) deallocate(QC)
    if (allocated(frqall)) deallocate(frqall)
    allocate(rad(nChan))
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
  ! 
  ! Allocate memory for workspace
  ! 
    if (allocated(profile)) deallocate(profile,a1)
    allocate(profile(classData%nInNN),a1(classData%ntrain), &
         stat=istatus)
    if (istatus /=0 ) then
       print*,trim(procName)//'Error: '
       print *,'Not enough memory for ', &
            'profile, a1; istatus:',istatus
       call exit(1)
    end if

  ! Check center frequencies for match
  ! with each channel used in classification

    iFstChan=1
    channelMissing=.false.
    do i=1,classData%nInNN
      chancaExists=.FALSE.
      if (classData%idInNN(i) == missingID) then ! Need to ID channel by frequency match
         do j=1,nChan
            reldif=abs(frqall(j)-classData%frqnn(i))/classData%frqnn(i)
            if (reldif < diffThreshold) then
               classData%idInNN(i)=j
               chancaExists=.TRUE.
            endif
         enddo
      elseif (classData%idInNN(i) > 0) then ! input is a channel identified by number
         reldif=abs(frqall(classData%idInNN(i))-classData%frqnn(i))/ &
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
    character(len=200) :: coefsClass
    integer :: inu=12
    TYPE(AtmosClassData_t) :: classDataIn
    namelist /classConfig/disable,channelMissingIsFatal,coefsClass

    procName = '[AtmosClassModule::loadAtmosClassData]:'
    if (dbg) print *, trim(procName)//' starting ...'

    open(inu,file=atmosClassConfigFile)
    read(inu,classConfig)
    close(inu)

    classDataIn%disable = disable
    classDataIn%channelMissingIsFatal = channelMissingIsFatal

    ! read data from file coefsClass and store it in classDataIn
    if (.not. disable) call rdwts(coefsClass,classDataIn)

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

    procName = '[AtmosClassModule::atmosClassify]:'
    if (dbg) print *, trim(procName)//' starting ...'

    select case (mode)
    case (MWIRmode)
       rad(1:MWdefin%nChan) = MWob%rad
       QC(1:MWdefin%nChan) = MWob%QC
       rad(MWdefin%nChan+1:MWdefin%nChan+IRdefin%nChan) = IRob%rad
       QC(MWdefin%nChan+1:MWdefin%nChan+IRdefin%nChan) = IRob%QC
       call classAtmNN(rad,QC,IRob%lat,IRob%lon,IRob%EIA,IRdefin%startTime, &
          surfaceClassIR,atmosClass)
    case (IRmode) ! with current dummy version of classatm, either MWIR or IR is OK
       rad(1:IRdefin%nChan) = IRob%rad
       QC(1:IRdefin%nChan) = IRob%QC
       call classAtmNN(rad,QC,IRob%lat,IRob%lon,IRob%EIA,IRdefin%startTime, &
          surfaceClassIR,atmosClass)
    case (MWmode)
       rad(1:MWdefin%nChan) = MWob%rad
       QC(1:MWdefin%nChan) = MWob%QC
       call classAtmNN(rad,QC,MWob%lat,MWob%lon,MWob%EIA,MWdefin%startTime, &
          surfaceClassMW,atmosClass)
    end select

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

  SUBROUTINE classAtmNN(rad,QC,xlat,xlon,eia,time,igeo,iclassatm)
    !------------------------------------------------------------      
    !   Subroutine to find the atmospheric class using a probablistic
    !   neural network.  An array
    !   (profile) of brightness temperatures is used as the input set,
    !   and the network classifies this array into a specific class.
    !   Upon the first entry to this routine,
    !   the results of that network training are read from a file.
    !   If any channels required for classification are missing, the
    !   subroutine returns iclassatm=nAtmClasses+1.
    !------------------------------------------------------------      
    !   Note that all local variables used in the network
    !   simulation need to be double precision.
    !------------------------------------------------------------      
    !
    !   Arguments (I=input, O=output):
    !
    !   rad:  (I) Array of brightness temperatures for current scene 
    !   QC:   (I) Array of quality control for each rad
    !   igeo: (I) Geographic land surface class (optional and not used; 
    !             provided only for back-compatibility; not needed operational)
    !
    !   iclassatm:   (O) Index of selected atmosphere class
    !
    !   U_Classatm: Logical unit number of network training file
    !   F_Classatm: Path of network training file
    !
    !   profile: Array of temperatures to be used in network
    !   simulation
    !
    !   The network has two layers.  Layer 1 has input weights given in
    !   array 'inpwts' and and a bias vector given in array 'bias1'.
    !   Layer 2 has layer weights given in array 'laywts'.  The output of
    !   layer 1 is given in array a1, and the output of layer 2 (i.e, the
    !   class of profile as determined by the network) is given in scalar
    !   iclassatm.
    !
    !---Arguments
    !---xlat,xlon,time,igeo are not used, and are only for interface 
    !---uniformity among classAtmNN versions
    REAL,    DIMENSION(:),           INTENT(IN)    :: rad
    INTEGER, DIMENSION(:),           INTENT(IN)    :: QC
    REAL,                  OPTIONAL, INTENT(IN)    :: xlat,xlon
    REAL,                  OPTIONAL, INTENT(IN)    :: eia
    INTEGER, DIMENSION(6), OPTIONAL, INTENT(IN)    :: time
    INTEGER,               OPTIONAL, INTENT(IN)    :: igeo
    INTEGER,                         INTENT(INOUT) :: iclassatm
    !---Local variables
    INTEGER                                   :: i
    REAL(kind=dbl)                            :: distvec,sumx,rmax
    REAL                                      :: secEIA
  !
  !   Load current scene brightness temperatures into classification
  !   input array
  !

    if (classData%disable) then
       iclassatm=1
       return
    endif

    if (channelMissing .or. &
       ANY(QC(classData%idInNN(iFstChan:classData%nInNN)) /= 0)) then
          iclassatm=classData%nAtmClasses+1
          RETURN
    endif

    if (classData%idInNN(1) == -1) then
       if (present(eia)) then
          secEIA=1./cos(eia*deg2rad)
          profile(1)=secEIA
       else
          print *,'err[AtmosClassModule::classAtmNN]: ',&
               ' Need EIA input but it is not present'
          call exit(1)
       endif
    endif

    profile(iFstChan:classData%nInNN)= &
        rad(classData%idInNN(iFstChan:classData%nInNN))
  !
  !   apply scale factors
  !
    profile=profile*classData%scale
  !
  !   layer 1
  !
    do i=1,classData%ntrain
       distvec=sum((classData%inpwts(i,:)-profile(:))**2)*classData%bias1(i)**2
       a1(i)=exp(distvec*(-1.))
    end do
  !
  !   layer 2
  !
    rmax=-1.
    iclassatm=MISSING_INT
    do i=1,classData%nAtmClasses
       sumx=dot_product(classData%laywts(i,:),a1(:))
       if (sumx > rmax) then
          iclassatm=i
          rmax=sumx
       end if
    end do
    RETURN
  END SUBROUTINE classAtmNN

!----------------------------------------------------------------------------
  
  SUBROUTINE rdwts(dataFile,classDataIn)
    !
    !   Subroutine to read the network training data.
    !   This is somewhat specific to the actual file in that it
    !   expects the structure of: 
    !     12 comment lines at the start of the file
    !     4 lines of header data
    !     1 line of layer 1 parameter information
    !     layer 1 data
    !     1 line of layer 2 parameter information
    !     layer 2 data.  
    !   As long as this structure is preserved, the routine should be able to
    !   successfully read files with different amounts of layer data.
    !
    !   Note that the data is read via list-directed I/O ("star" format)
    !   rather than through formatted I/O.  List directed I/O is not as
    !   robust as formatted, and can be a source of potential program
    !   failure.  
    !   
    !
    character(len=*),       intent(in)    :: dataFile
    TYPE(AtmosClassData_t), intent(inout) :: classDataIn
    ! Local constants:
    INTEGER, PARAMETER :: mxHeader=99
    CHARACTER(LEN=1), PARAMETER :: commentMarker='$'
    ! Local variables:
    INTEGER          :: iun
    INTEGER          :: nHeader
    INTEGER          :: nInNN,nAtmClasses
    INTEGER          :: istatus
    INTEGER          :: nlayers,ilayer,ntrain,ntrain2,i,j
    LOGICAL          :: oldFmt
    CHARACTER(LEN=6) :: tranfunc
    CHARACTER(LEN=8) :: lBuffer
    iun = getUnit()
    open(iun,file=dataFile)
    do i=1,mxHeader
       read(iun,'(a8)')lBuffer
       if (lBuffer(1:1) /= commentMarker) exit
    end do
    if (i > mxHeader) then
       print *,'err[AtmosClassModule::rdwts] '
       print *,'More header lines than expected; i,mxHeader:',i,mxHeader
       call exit(1)
    else if (i == 13) then ! This file has old format
       oldFmt=.true.
    else
       oldFmt=.false.
    end if
    read(lBuffer,*) nInNN,nAtmClasses
    if (allocated(classDataIn%frqnn)) deallocate(classDataIn%frqnn)
    if (allocated(classDataIn%idInNN)) deallocate(classDataIn%idInNN)
    if (allocated(classDataIn%scale)) deallocate(classDataIn%scale)
    if (allocated(classDataIn%idAtmClasses)) deallocate(classDataIn%idAtmClasses)
    allocate(classDataIn%idInNN(nInNN),classDataIn%scale(nInNN), &
       classDataIn%frqnn(nInNN), &
       classDataIn%idAtmClasses(nAtmClasses),STAT=istatus)
    if (istatus /=0 ) then
       print *,'err[AtmosClassModule::rdwts] '
       print *,'Not enough memory for idInNN&frqnn&idAtmClasses; istatus:',istatus
       call exit(1)
    end if
    if (oldFmt) then
       classDataIn%scale(1:nInNN)= 1. ! scale=1 leaves data unchanged
       classDataIn%idInNN(1:nInNN)= missingID
    else
       read(iun,*) classDataIn%scale(1:nInNN)
       read(iun,*) classDataIn%idInNN(1:nInNN)
    endif
    read(iun,*) classDataIn%frqnn(1:nInNN)
    read(iun,*) classDataIn%idAtmClasses(1:nAtmClasses)
    read(iun,*) nlayers
    if (nlayers /= 2) then
       print *,'err[AtmosClassModule::rdwts] '
       print *,' Training results file has',nlayers,' layers, not 2'
       call exit(1)
    end if
    read(iun,*) ilayer,nInNN,tranfunc,ntrain
    if (ilayer /= 1) then
       print *,'err[AtmosClassModule::rdwts] '
       print *,'Expected layer number 1 in training', &
              ' results file, got',ilayer,' instead'
       call exit(1)
    end if
    if (allocated(classDataIn%inpwts)) deallocate(classDataIn%inpwts)
    if (allocated(classDataIn%bias1)) deallocate(classDataIn%bias1)
    allocate(classDataIn%inpwts(ntrain,nInNN),classDataIn%bias1(ntrain), &
       STAT=istatus)
    if (istatus /=0 ) then
       print *,'err[AtmosClassModule::rdwts] '
       print *,'Not enough memory for inpwts, ', &
              'and bias1; istatus:',istatus
       call exit(1)
    end if
    read(iun,*) ((classDataIn%inpwts(i,j),i=1,ntrain),j=1,nInNN)
    read(iun,*) (classDataIn%bias1(i),i=1,ntrain)
    read(iun,*) ilayer,ntrain2,tranfunc,nAtmClasses
    if (ilayer /= 2) then
       print *,'err[AtmosClassModule::rdwts] '
       print *,'Expected layer number 2 in training', &
              ' results file, got',ilayer,' instead'
       call exit(1)
    end if
    if (ntrain2 /= ntrain) then
       print *,'err[AtmosClassModule::rdwts] '
       print *,'Found',ntrain2,' training channels used in', &
              ' training results file layer 2, not the',ntrain, &
              ' found in layer 1'
       call exit(1)
    end if
    if (allocated(classDataIn%laywts)) deallocate(classDataIn%laywts)
    allocate(classDataIn%laywts(nAtmClasses,ntrain),STAT=istatus)
    if (istatus /=0 ) then
       print *,'err[AtmosClassModule::rdwts] '
       print *,'Not enough memory for laywts; istatus:',istatus
       call exit(1)
    end if
    read(iun,*) ((classDataIn%laywts(i,j),i=1,nAtmClasses),j=1,ntrain)
    close(iun)

    classDataIn%nInNN = nInNN
    classDataIn%nAtmClasses = nAtmClasses
    classDataIn%ntrain = ntrain

    RETURN
  END SUBROUTINE rdwts

END MODULE AtmosClassModule
