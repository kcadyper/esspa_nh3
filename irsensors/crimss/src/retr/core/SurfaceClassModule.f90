MODULE SurfaceClassModule
  !
  ! Module for surface classification
  !
  ! Subroutines:
  !
  !     MWsurfaceClassInit
  !     MWsurfaceClassify
  !     loadMWsurfaceClassData
  !     putMWsurfaceClassData
  !     IRsurfaceClassInit
  !     IRsurfaceClassify
  !     loadIRsurfaceClassData
  !     putIRsurfaceClassData
  !
  ! Derived data types:
  !
  !     MWsurfaceClassData_t
  !     IRsurfaceClassData_t
  !
  ! USE:
  !
  !     None
  !
  ! yhe@aer.com, 03/16/2016
  !
  USE ControlStructure, Only: &
       GenControl_t

  USE MWobsStructure, Only: &
       MWdefin_t, &
       MWob_t

  USE IRobsStructure, Only: &
       IRdefin_t, &
       IRob_t

  USE SurfaceTypingModule, Only: &
       setSfc

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: &
       MWsurfaceClassInit, &
       MWsurfaceClassify, &
       loadMWsurfaceClassData, &
       putMWsurfaceClassData, &
       IRsurfaceClassInit, &
       IRsurfaceClassify, &
       loadIRsurfaceClassData, &
       putIRsurfaceClassData, &
       destroySurfaceClass, &
       MWsurfaceClassData_t, &
       IRsurfaceClassData_t

  TYPE MWsurfaceClassData_t
     real :: maxLandInOceanClass
  END TYPE MWsurfaceClassData_t

  TYPE IRsurfaceClassData_t
     real :: maxLandInOceanClass
  END TYPE IRsurfaceClassData_t

  integer, parameter :: mxLength=256
  logical :: dbg=.false.

  TYPE(MWsurfaceClassData_t), save :: MWclassData
  TYPE(IRsurfaceClassData_t), save :: IRclassData
  logical, save :: MWdataReady = .FALSE.
  logical, save :: IRdataReady = .FALSE.

CONTAINS

  subroutine MWsurfaceClassInit(genControl)
    TYPE(genControl_t), intent(in) :: genControl
    !Local
    character (len=mxLength) :: procName

    procName = '[SurfaceClassModule::MWsurfaceClassInit]:'
    if (dbg) print *, trim(procName)//'starting ...'

    if (.not. MWdataReady) &
       call loadMWsurfaceClassData(genControl%MWsurfaceClassConfig)

    if (dbg) print *, trim(procName)//' ending ...'
  end subroutine MWsurfaceClassInit

  subroutine loadMWsurfaceClassData(classConfigFile)
  ! This provides an option of a public method to load the static data
  ! at a level higher than MWsurfaceClassInit
    character(len=*), intent(in)    :: classConfigFile
    !Local
    character (len=mxLength) :: procName
    character(len=200) :: coefsClass
    real :: pLandMaxOc
    integer :: inu=12
    TYPE(MWsurfaceClassData_t) :: classDataIn
    namelist /classConfig/coefsClass,pLandMaxOc

    procName = '[SurfaceClassModule::loadMWsurfaceClassData]:'
    if (dbg) print *, trim(procName)//' starting ...'

    open(inu,file=classConfigFile)
    read(inu,classConfig)
    close(inu)

    ! read data from file coefsClass and store it in classDataIn
    classDataIn%maxLandInOceanClass = pLandMaxOc

    call putMWsurfaceClassData(classDataIn)

    if (dbg) print *, trim(procName)//' ending ...'

  end subroutine loadMWsurfaceClassData

  subroutine putMWsurfaceClassData(classDataIn)
  ! This provides an option of a public method to insert the static data
    TYPE(MWsurfaceClassData_t), intent(in) :: classDataIn
    !Local
    character (len=mxLength) :: procName

    procName = '[SurfaceClassModule::putMWsurfaceClassData]:'
    if (dbg) print *, trim(procName)//' starting ...'

    MWclassData = classDataIn

    MWdataReady = .TRUE.

    if (dbg) print *, trim(procName)//' ending ...'

  end subroutine putMWsurfaceClassData

  subroutine MWsurfaceClassify( &
       landFrac, &
       MWdefin, &
       MWob, &
       surfaceClassMW, &
       isLandMW &
       )
    real, intent(in) :: landFrac
    TYPE(MWdefin_t), intent(in) :: MWdefin
    TYPE(MWob_t), intent(in) :: MWob
    integer, intent(out) :: surfaceClassMW
    logical, intent(out) :: isLandMW
    !Local
    integer :: iexit
    character (len=mxLength) :: procName

    procName = '[SurfaceClassModule::MWsurfaceClassify]:'
    if (dbg) print *, trim(procName)//'starting ...'
    iexit=0
    call setSfc(surfaceClassMW,landFrac,MWclassData%maxLandInOceanClass)
    IF (landFrac > MWclassData%maxLandInOceanClass) THEN
       isLandMW=.TRUE.
    ELSE
       isLandMW=.FALSE.
    ENDIF
    if (dbg) print *, trim(procName)//' ending ...'
  end subroutine MWsurfaceClassify

  subroutine IRsurfaceClassInit(genControl)
    TYPE(genControl_t), intent(in) :: genControl
    !Local
    character (len=mxLength) :: procName

    procName = '[SurfaceClassModule::IRsurfaceClassInit]:'
    if (dbg) print *, trim(procName)//'starting ...'

    if (.not. IRdataReady) &
       call loadIRsurfaceClassData(genControl%IRsurfaceClassConfig)

    if (dbg) print *, trim(procName)//' ending ...'
  end subroutine IRsurfaceClassInit

  subroutine loadIRsurfaceClassData(classConfigFile)
  ! This provides an option of a public method to load the static data
  ! at a level higher than IRsurfaceClassInit
    character(len=*), intent(in)    :: classConfigFile
    !Local
    character (len=mxLength) :: procName
    character(len=200) :: coefsClass
    real :: pLandMaxOc
    integer :: inu=12
    TYPE(IRsurfaceClassData_t) :: classDataIn
    namelist /classConfig/coefsClass,pLandMaxOc

    procName = '[SurfaceClassModule::loadIRsurfaceClassData]:'
    if (dbg) print *, trim(procName)//' starting ...'

    open(inu,file=classConfigFile)
    read(inu,classConfig)
    close(inu)

    ! read data from file coefsClass and store it in classDataIn
    classDataIn%maxLandInOceanClass = pLandMaxOc

    call putIRsurfaceClassData(classDataIn)

    if (dbg) print *, trim(procName)//' ending ...'

  end subroutine loadIRsurfaceClassData

  subroutine putIRsurfaceClassData(classDataIn)
  ! This provides an option of a public method to insert the static data
    TYPE(IRsurfaceClassData_t), intent(in) :: classDataIn
    !Local
    character (len=mxLength) :: procName

    procName = '[SurfaceClassModule::putIRsurfaceClassData]:'
    if (dbg) print *, trim(procName)//' starting ...'

    IRclassData = classDataIn

    IRdataReady = .TRUE.

    if (dbg) print *, trim(procName)//' ending ...'

  end subroutine putIRsurfaceClassData

  subroutine IRsurfaceClassify( &
       landFrac, &
       IRdefin, &
       IRob, &
       surfaceClassIR, &
       isLandIR &
       )	
    real, intent(in) :: landFrac
    TYPE(IRdefin_t), intent(in) :: IRdefin
    TYPE(IRob_t), intent(in) :: IRob
    integer, intent(out) :: surfaceClassIR
    logical, intent(out) :: isLandIR
    !Local
    integer :: iexit
    character (len=mxLength) :: procName

    procName = '[SurfaceClassModule::IRsurfaceClassify]:'
    if (dbg) print *, trim(procName)//'starting ...'
    iexit=0
    call setSfc(surfaceClassIR,landFrac,IRclassData%maxLandInOceanClass)
    IF (landFrac > IRclassData%maxLandInOceanClass) THEN
       isLandIR=.TRUE.
    ELSE
       isLandIR=.FALSE.
    ENDIF
    if (dbg) print *, trim(procName)//' ending ...'
  end subroutine IRsurfaceClassify
  
  subroutine destroySurfaceClass()
    !Local
    character (len=mxLength) :: procName

    procName = '[SurfaceClassModule::destroySurfaceClass]'
    if (dbg) print *, trim(procName)//' starting ...'
    MWdataReady = .FALSE.
    IRdataReady = .FALSE.
    if (dbg) print *, trim(procName)//' ending ...'
  end subroutine destroySurfaceClass

END MODULE SurfaceClassModule
