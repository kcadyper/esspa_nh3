MODULE FiniteDiffModule
  !
  ! Module: Finite difference calculation for cloud scattering. All
  !         the subroutines are relocated from the main driver,
  !         retr_sig.f90 in SIGMA branch.
  !
  ! Subroutine:
  !
  !     initFiniteDiff
  !     loadFiniteDiff
  !     setProperties
  !     finiteDiff  (calls ossdrv_ir)
  !     destroyFiniteDiff
  ! 
  ! Derived data type:
  !
  !     None
  !
  ! USE:
  !
  !     oss_ir_module
  !     ToolboxModule
  !
  ! yhe@aer.com, 07/06/2016
  !
  USE oss_ir_module, Only: &
       ossdrv_ir

  USE ToolboxModule, Only: &
       getUnit

  USE StateIndexModule, Only: &
       StateIndex_t, &
       getNmol, &
       getAtmosVectorLength

  USE VertCoord, ONLY: mxCoordTyp,Pcoord


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: &
       initFiniteDiff, &
       loadFiniteDiff, &
       finiteDiff, &
       destroyFiniteDiff

  integer, parameter :: mxLength = 256
  integer :: nMol
  real :: delCldLiqTop,delCldLiqThk,delCldLiqAmt,delCldLiqDeff
  real :: delCldIceTop,delCldIceThk,delCldIceAmt,delCldIceDeff
  real :: delTskin
  real, dimension(:), allocatable :: xGFD,yFD
  real, dimension(:,:), allocatable :: xktGFD
  logical :: loadFiniteDiffDone=.false.
  logical :: dbg=.false.
  TYPE(StateIndex_t) :: IG
  TYPE(StateIndex_t) :: NG
  CHARACTER(LEN=mxCoordTyp) :: vCoordTyp
  integer :: nChan  !MW+IR
  integer :: nParGAtm

CONTAINS
  subroutine initFiniteDiff(IGin,NGin,vCoordTypin,nChanin)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   initFiniteDiff
!
! PURPOSE:
!
!   Initialize state vector, radiances, and Jacobians to zero for 
!   finite difference calculations
!
! SYNTAX:
!
!   CALL initFiniteDiff()
!
! ARGUMENTS:
!
!   
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
    TYPE(StateIndex_t), intent(in) :: IGin
    TYPE(StateIndex_t), intent(in) :: NGin
    CHARACTER(LEN=*), intent(in) :: vCoordTypin
    integer, intent(in) :: nChanin
    !Local
    character(len=mxLength) :: procName

    procName = '[FiniteDiffModule::initFiniteDiff]:'
    if (dbg) print *, trim(procName)//' starting ...'
    IG = IGin
    NG = NGin
    nMol = getNmol(NG)
    nChan = nChanIn
    nParGAtm = getAtmosVectorLength(NG)
    allocate(xktGFD(nParGAtm,nchan))
    xktGFD(:,:)=0.
    allocate(xGFD(nParGAtm))
    xGFD(:)=0.
    allocate(yFD(nchan))
    yFD=0.

    if (dbg) print *, trim(procName)//' ending ...'
    return

  end subroutine initFiniteDiff

  subroutine loadFiniteDiff(file_algcfg)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   loadFiniteDiff
!
! PURPOSE:
!
!   Load tunable controls on finite difference calculations
!
! SYNTAX:
!
!   CALL loadFiniteDiff(file_algcfg)
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
    real               :: delCldLiqTop,delCldLiqThk,delCldLiqAmt,delCldLiqDeff
    real               :: delCldIceTop,delCldIceThk,delCldIceAmt,delCldIceDeff
    real               :: delTskin

    !Local
    integer :: U_algconfig
    character(len=mxLength) :: procName
  
    namelist /finitediffdata/ &
         delCldLiqTop,delCldLiqThk,delCldLiqAmt,delCldLiqDeff,&
         delCldIceTop,delCldIceThk,delCldIceAmt,delCldIceDeff,&
         delTskin
  
    procName = '[FiniteDiffModule::loadFiniteDiff]:'
    if (dbg) print *, trim(procName)//' starting ...'
    if(loadFiniteDiffDone) return
    
    !------------------------------------------------------------------------
    ! Read namelist file items into local variables
    !------------------------------------------------------------------------
    U_algconfig = getUnit()
    open(U_algconfig,file=file_algcfg)
    read(U_algconfig,finitediffdata)
    close(U_algconfig)
    !------------------------------------------------------------------------
    ! Set global variables
    !------------------------------------------------------------------------
    call setPropertiesFiniteDiff(delCldLiqTop,delCldLiqThk,delCldLiqAmt,&
         delCldLiqDeff,delCldIceTop,delCldIceThk,delCldIceAmt,delCldIceDeff,&
         delTskin)
  
    if (dbg) print *, trim(procName)//' ending ...'
    return
  end subroutine loadFiniteDiff

  subroutine setPropertiesFiniteDiff(delCldLiqTop_in,delCldLiqThk_in, &
       delCldLiqAmt_in,delCldLiqDeff_in,delCldIceTop_in,delCldIceThk_in, &
       delCldIceAmt_in,delCldIceDeff_in,delTskin_in)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   setPropertiesFiniteDiff
!
! PURPOSE:
!
!   Set global variables from input arguments, for tunable controls on 
!   finite difference calculations
!
! SYNTAX:
!
!   CALL setPropertiesFiniteDiff(delCldLiqTop_in, delCldLiqThk_in, 
!      delCldLiqAmt_in, delCldLiqDeff_in, delCldIceTop_in, 
!      delCldIceThk_in, delCldIceAmt_in, delCldIceDeff_in, 
!      delTskin_in)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   delCldLiqTop_in   REAL  Finite diff step for liquid cloud top 
!                           in fraction of pres grid spacing 
!   delCldLiqThk_in   REAL  Finite diff step for liquid cloud 
!                           thickness in fraction of pres grid 
!                           spacing 
!   delCldLiqAmt_in   REAL  Finite diff step for liquid cloud 
!                           thickness (kg/m2) 
!   delCldLiqDeff_in  REAL  Finite diff step for liquid cloud 
!                           effective diameter, fractional 
!   delCldIceTop_in   REAL  Finite diff step for ice cloud top in 
!                           fraction of pres grid spacing 
!   delCldIceThk_in   REAL  Finite diff step for ice cloud 
!                           thickness in fraction of pres grid 
!                           spacing 
!   delCldIceAmt_in   REAL  Finite diff step for ice cloud 
!                           thickness (kg/m2) 
!   delCldIceDeff_in  REAL  Finite diff step for ice cloud 
!                           effective diameter, fractional 
!   delTskin_in       REAL  Finite diff step for surface skin 
!                           temperature (K)
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !--I/O variables
    real,    intent(in)        :: delCldLiqTop_in,delCldLiqThk_in
    real,    intent(in)        :: delCldLiqAmt_in,delCldLiqDeff_in
    real,    intent(in)        :: delCldIceTop_in,delCldIceThk_in
    real,    intent(in)        :: delCldIceAmt_in,delCldIceDeff_in
    real,    intent(in)        :: delTskin_in
    
    !Local
    character(len=mxLength) :: procName

    procName = '[FiniteDiffModule::setProperties]:'
    if (dbg) print *, trim(procName)//' starting ...'
    
    !------------------------------------------------------------------------
    ! Set global variables from input arguments
    !------------------------------------------------------------------------
    delCldLiqTop   = delCldLiqTop_in
    delCldLiqThk   = delCldLiqThk_in
    delCldLiqAmt   = delCldLiqAmt_in
    delCldLiqDeff  = delCldLiqDeff_in
    delCldIceTop   = delCldIceTop_in
    delCldIceThk   = delCldIceThk_in
    delCldIceAmt   = delCldIceAmt_in
    delCldIceDeff  = delCldIceDeff_in
    delTskin       = delTskin_in

    loadFiniteDiffDone=.true.
    if (dbg) print *, trim(procName)//' ending ...'
    return
  end subroutine setPropertiesFiniteDiff
  
  SUBROUTINE finiteDiff(xG,EmRfGrid,EmRf,obsang,sunang,azAngle,lat,obsLevel, &
                        yref,xkt,xkEmRf,pland,press,NR)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   finiteDiff
!
! PURPOSE:
!
!   Numerical computation of Jacobians.
!
! SYNTAX:
!
!   CALL finiteDiff(xG, EmRf, obsang, sunang, yref, xkt, xkEmRf, 
!   pland)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   EmRfGrid REAL IR emissivity/reflectivity  hinge points
!   EmRf     REAL IR emissivity/reflectivity at hinge points
!   obsang   REAL  viewing angle
!   sunang   REAL  solar angle
!   sunang   REAL  solar angle
!   azAngle  REAL  relative azimuth angle
!   lat      REAL  observation 
!   obsLevel INTEGER observer level
!   yref     REAL  Reference radiances
!   pland    REAL  Fraction of land (versus water) of the field of 
!                 view 
!   NR       StateIndex_t  Retrieval space state lengths
!   
!   INPUTS/OUTPUTS:
!   
!   xG      REAL  State vector in geophysical space
!   xkt     REAL  Jacobian matrix transpose 
!   xkEmRf  REAL  Jacobians for IR emissivity/reflectivity
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !---Input variables

    REAL, DIMENSION(:),     intent(inout) :: xG
    REAL, DIMENSION(:),     intent(in)    :: EmRfGrid
    REAL, DIMENSION(:,:),   intent(in)    :: EmRf
    REAL,                   intent(in)    :: pland,obsang,sunang
    REAL,                   intent(in)    :: azAngle, lat
    integer,                intent(inout) :: obsLevel
    REAL, DIMENSION(:),     intent(in)    :: yRef
    REAL, DIMENSION(:),     intent(in)    :: press
    TYPE(StateIndex_t),     intent(in)    :: NR
    !---Output variables
    REAL, DIMENSION(:,:),   intent(inout) :: xkt
    REAL, DIMENSION(:,:,:), intent(inout) :: XkEmRf
    !---Local variables
    INTEGER                          :: i
    INTEGER                          :: ilevLiq,ilevIce
    character(len=*), parameter      :: procName='[FiniteDiffModule::finiteDiff]:'
    if (dbg) print *, procName //' starting ...'

    xkt(:,:)=0.
    xGFD=0.
    if (NR%cldIce > 0) then
       i=1
       do while( press(i) < xG(IG%cldIce) )
          i=i+1
          if (i .GE. size(press)) exit 
       enddo
       ilevIce=i
       xGFD(IG%cldIce)=delCldIceTop*(press(ilevIce)-press(ilevIce-1))
       xGFD(IG%cldIce+2)=delCldIceAmt
       xGFD(IG%cldIce+3)=delCldIceDeff*xG(IG%cldIce+3)
    endif
    if (NR%cldLiq > 0) then
       i=1
       do while( press(i) < xG(IG%cldLiq) )
          i=i+1
          if (i .GE. size(press)) exit 
       enddo
       ilevLiq=i
       xGFD(IG%cldLiq)=delCldLiqTop*(press(ilevLiq)-press(ilevLiq-1))
       xGFD(IG%cldLiq+2)=delCldLiqAmt
       If (NG%cldLiq .GT. 3) xGFD(IG%cldLiq+3)=delCldLiqDeff*xG(IG%cldLiq+3)
    endif
    if (NR%Tskin > 0) then
       xGFD(IG%Tskin)=delTskin
    endif
       
    do i=1,size(xG)
       if(xGFD(i) == 0.) cycle
       xG(i)=xG(i)+xGFD(i)
       IF (TRIM(vCoordTyp) == Pcoord) THEN
          call ossdrv_ir(xG,EmRfGrid,EmRf,obsang,sunang,azAngle,obsLevel,yFD,xktGFD,xkEmRf,lat,&
            tempIndex=IG%temp,tSkinIndex=IG%tskin,pSurfIndex=IG%psfc,varMolIndex=IG%mol(1:nMol))
       ELSE
          call ossdrv_ir(xG,EmRfGrid,EmRf,obsang,sunang,azAngle,obsLevel,yFD,xktGFD,xkEmRf,lat,&
            tempIndex=IG%temp,tSkinIndex=IG%tskin,pSurfIndex=IG%psfc,varMolIndex=IG%mol(1:nMol),&
            puser=press)
       ENDIF
       xG(i)=xG(i)-xGFD(i)
       xkt(i,:)=(yFD(:)-yRef(:))/xGFD(i)
    enddo

    if (dbg) print *, trim(procName)//' ending ...'
    return

  END SUBROUTINE finiteDiff
 
  SUBROUTINE destroyFiniteDiff()

!<f90Subroutine>********************************************************
!
! NAME:
!
!   destroyFiniteDiff
!
! PURPOSE:
!
!   Deallocate arrays related to finite differencing.
!
! SYNTAX:
!
!   CALL destroyFiniteDiff()
!
! ARGUMENTS:
!
!   
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
    !Local
    character(len=mxLength) :: procName

    procName = '[FiniteDiffModule::destroyFiniteDiff]:'
    if (dbg) print *, trim(procName)//' starting ...'
    if (allocated(xktGFD))         deallocate(xktGFD)
    if (allocated(xGFD))           deallocate(xGFD)
    if (allocated(yFD))            deallocate(yFD)
    loadFiniteDiffDone=.false.
    if (dbg) print *, trim(procName)//' ending ...'
    return
  end subroutine  destroyFiniteDiff

END MODULE FiniteDiffModule
