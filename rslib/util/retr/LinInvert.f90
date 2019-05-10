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

MODULE LinInvert

! <f90Module>***********************************************************
!
! NAME:
!
!   LinInvert
!
! PURPOSE:
!
!   Tools for linear retrieval
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

! Tools for linear retrieval
  USE StateIndexModule, ONLY: StateIndex_t,getNMol,genIndices,getVectorLength
  USE LvlInterp
  USE InvertModule, ONLY: invrt1_GJ_factors
  USE constants, ONLY: MISSING_REAL
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: initBuildLin,buildLinearCoef,linRetr,getLinInvData,destroyLinInvert

! Global private data
  integer, parameter :: mxLength=256
  TYPE(StateIndex_t) :: NGR,IGR
  INTEGER, DIMENSION(:),   ALLOCATABLE :: mapGeo
  REAL,    DIMENSION(:),   ALLOCATABLE :: xBakCompr
  REAL,    DIMENSION(:,:), ALLOCATABLE :: uCompr

CONTAINS

  SUBROUTINE initBuildLin(IG,NG,nParG,NR_orig,IGRout,NGRout,nParGon)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   initBuildLin
!
! PURPOSE:
!
!   Define linear retrieval state vector content in geophysical space, and make
!   an index map from this content to the content of the full geo space
!
! SYNTAX:
!
!   CALL initBuildLin(IG, NG, nParG, NR_orig, IGRout, NGRout, nParGon)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   IG       TYPE(STATEINDEX_T)  Starting indices for sections of geophysical
!                                state vector
!   NG       TYPE(STATEINDEX_T)  Number of elements for sections of
!                                geophysical state vector
!   nParG    INTEGER             Total number of elements in geophysical
!                                state vector
!   NR_orig  TYPE(STATEINDEX_T)  Number of retrievable elements for sections
!                                of retrieval state vector
!
!   INPUTS/OUTPUTS:
!
!   IGRout   TYPE(STATEINDEX_T)  Starting indices for retrieved sections of
!                                geophysical state vector
!   NGRout   TYPE(STATEINDEX_T)  Number of elements for retrieved sections of
!                                geophysical state vector
!   nParGon  INTEGER             Total number of retrieved elements in
!                                geophysical state vector
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>


    ! Arguments:
    TYPE(StateIndex_t),      INTENT(IN)    :: IG      ! indices in full geo space
    TYPE(StateIndex_t),      INTENT(IN)    :: NG      ! sizes in full geo space
    INTEGER,                 INTENT(IN)    :: nParG   ! number of full geo-space
    TYPE(StateIndex_t),      INTENT(IN)    :: NR_orig ! which are retrievable
    TYPE(StateIndex_t),      INTENT(INOUT) :: IGRout  ! indices for linear retrieval
    TYPE(StateIndex_t),      INTENT(INOUT) :: NGRout  ! sizes for linear retrieval
    INTEGER,                 INTENT(INOUT) :: nParGon ! number of geo-space
                                                      !   parameters to retrieve

    ! Local variables:
    INTEGER :: i,j,nmol,ix

    ! Define linear retrieval state vector content in geophysical space,
    ! excluding anything not being retrieved
    NGR=NG
    IF (NR_orig%temp   <= 0) NGR%temp=0
    IF (NR_orig%Tskin  <= 0) NGR%Tskin=0
    IF (NR_orig%Psfc   <= 0) NGR%Psfc=0
    IF (NR_orig%cldLiq <= 0) NGR%cldLiq=0
    IF (NR_orig%cldIce <= 0) NGR%cldIce=0
    nmol=getNMol(NG)
    DO i=1,nmol
      IF (NR_orig%mol(i) <= 0) NGR%mol(i)=0
    ENDDO

    IGR = genIndices(NGR)

    NGRout=NGR
    IGRout=IGR

    nParGon=getVectorLength(NGR)
    ALLOCATE (mapGeo(nParGon),uCompr(nParGon,nParG),xBakCompr(nParGon))
    ! Make an index map from the content of the linear retrieval state vector
    ! to the content of the full geo space
    ! A coding implementation with one executable line with implicit do loops
    !  gave run-time errors
    ix=0
    DO i=IG%Temp,IG%Temp+NGR%Temp-1
      ix=ix+1
      mapGeo(ix)=i
    ENDDO
    DO i=IG%Tskin,IG%Tskin+NGR%Tskin-1
      ix=ix+1
      mapGeo(ix)=i
    ENDDO
    DO i=IG%Psfc,IG%Psfc+NGR%Psfc-1
      ix=ix+1
      mapGeo(ix)=i
    ENDDO
    DO j=1,nmol
      DO i=IG%mol(j),IG%mol(j)+NGR%mol(j)-1
        ix=ix+1
        mapGeo(ix)=i
      ENDDO
    ENDDO
    DO i=IG%cldLiq,IG%cldLiq+NGR%cldLiq-1
      ix=ix+1
      mapGeo(ix)=i
    ENDDO
    DO i=IG%cldIce,IG%cldIce+NGR%cldIce-1
      ix=ix+1
      mapGeo(ix)=i
    ENDDO
    DO i=IG%wind,IG%wind+NGR%wind-1
      ix=ix+1
      mapGeo(ix)=i
    ENDDO

    RETURN

  END SUBROUTINE initBuildLin


  SUBROUTINE buildLinearCoef(xkt,Sy,Sx,Y,kchan,xGes,vmtx,nPar,nChan,nParGon, &
       xBakG,xOffG,xGainG)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   buildLinearCoef
!
! PURPOSE:
!
!   Build the coefficients for a linear MAP retrieval
!
! SYNTAX:
!
!   CALL buildLinearCoef(xkt, Sy, Sx, Y, kchan, xGes, vmtx, nPar,
!      nChan, nParGon, xBakG, xOffG, xGainG)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   xkt      REAL     Jacobian matrix transpose
!   Sy       REAL     Measurement error variance
!   Sx       REAL     Background error covariance
!   Y        REAL     Radiometric data computed from state vector
!                     if idx(1:nch) that kchan(idx)>0 otherwise kchan=0
!                     then compression means that y(j) corresponds to Ym(idx(j))
!                     or kcnah(idx(j))>0
!                     the same for Jacobians
!
!   kchan    LOGICAL  Channel on/off mask
!   xGes     REAL     Retrieval state vector guess (current estimate)
!   vmtx     REAL     Transformation matrix between geophysical and retrieval
!                     spaces
!   nPar     INTEGER  Total number of elements in retrieval vector
!   nChan    INTEGER  Number of channels
!   nParGon  INTEGER  Total number of retrieved elements in geophysical state
!                     vector
!   xBakG    REAL     Background state vector in geophysical space
!
!   INPUTS/OUTPUTS:
!
!   xOffG    REAL     linear offset coefficients in geophysical space
!   xGainG   REAL     linear gain coefficients in geophysical space
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

  ! Build the coefficients for a linear MAP retrieval

    ! Arguments:
    REAL,    DIMENSION(:,:), INTENT(IN)    :: xkt     ! Jacobian matrix transposed
    REAL,    DIMENSION(:),   INTENT(IN)    :: Sy      ! measurement error variance
    REAL,    DIMENSION(:,:), INTENT(IN)    :: Sx      ! background error covariance
    REAL,    DIMENSION(:),   INTENT(IN)    :: Y       ! guess radiances
    LOGICAL, DIMENSION(:),   INTENT(IN)    :: kchan   ! flags for channels "on"
    REAL,    DIMENSION(:),   INTENT(IN)    :: xGes    ! guess state vector
                                                      !   in retrieval space
    REAL,    DIMENSION(:,:), INTENT(IN)    :: vmtx    ! transformation matrix
    INTEGER,                 INTENT(IN)    :: nPar    ! length of state vector
    INTEGER,                 INTENT(IN)    :: nChan   ! number of channels
    INTEGER,                 INTENT(IN)    :: nParGon ! number of geo-space
                                                      !   parameters to retrieve
    REAL,    DIMENSION(:),   INTENT(IN)    :: xBakG   ! background geo state vector
    REAL,    DIMENSION(:),   INTENT(INOUT) :: xOffG   ! linear offset coefficients
    REAL,    DIMENSION(:,:), INTENT(INOUT) :: xGainG  ! linear gain coefficients

    ! Local variables:
    INTEGER :: ich,j
    REAL,    DIMENSION(nPar)       :: xOff
    REAL,    DIMENSION(nPar,nchan) :: xGain
    REAL,    DIMENSION(nchan)      :: Yfilt

    ! Remove channels not turned on
    ! Y is compressed
    ich=0
    DO j=1,size(Y)
      IF (kchan(j)) THEN
        ich=ich+1
        Yfilt(ich)=Y(ich)
      ENDIF
    ENDDO
    IF (ich /= nChan) THEN
      print *,'err:[LinInvert::buildLinearCoef] ', &
         'Inconsistent number of channels: ',ich,nChan
      CALL exit(1)
    ENDIF

    ! Build coefficients in retrieval space
    CALL invrt1_GJ_factors(xkt,Sy,Sx,Yfilt,xGes,nPar,nChan,xOff,xGain)

    ! Eliminate un-retrieved geo-space state vector elements
    uCompr(:,1:nPar)=vmtx(mapGeo,1:nPar)
    xBakCompr=xBakG(mapGeo)

    ! Convert coefficients to geo space
    xOffG=MATMUL(uCompr(:,1:nPar),xOff(1:nPar))+xBakCompr
    xGainG=MATMUL(uCompr(:,1:nPar),xGain(1:nPar,1:nChan))
    RETURN

  END SUBROUTINE buildLinearCoef


  SUBROUTINE linRetr(xOffG,xGainG,QC,IGRin,NGRin,vCoordTyp,pref,pSfc,Tsfc,iH2O, &
      nMembers,radGroup,xGroup)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   linRetr
!
! PURPOSE:
!
!   Execute linear retrieval
!
! SYNTAX:
!
!   CALL linRetr(xOffG, xGainG, QC, IGRin, NGRin, pref, pSfc, Tsfc,
!      iH2O, nMembers, radGroup, xGroup)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   xOffG     REAL                linear offset coefficients in geophysical
!                                 space
!   xGainG    REAL                linear gain coefficients in geophysical
!                                 space
!   QC        INTEGER             quality control index
!   IGRin     TYPE(STATEINDEX_T)  Starting indices for retrieved sections of
!                                 geophysical state vector
!   NGRin     TYPE(STATEINDEX_T)  Number of elements for retrieved sections
!                                 of geophysical state vector
!   pref      REAL                Pressure on atmospheric grid levels
!   pSfc      REAL                Surface pressure
!   Tsfc      REAL                Surface-level air temperature
!   iH2O      INTEGER             Index for water vapor in molecules part of
!                                 state vector
!   nMembers  INTEGER             Number of members per group
!   radGroup  REAL                Radiometric data for the group
!
!   INPUTS/OUTPUTS:
!
!   xGroup    REAL                Retrieval state vector for the group
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

  ! Execute linear retrieval

    USE VertCoord, ONLY: mxCoordTyp,Pcoord

    ! Arguments:
    REAL,    DIMENSION(:),   INTENT(IN)    :: xOffG    ! linear offset coefficients
    REAL,    DIMENSION(:,:), INTENT(IN)    :: xGainG   ! linear gain coefficients
    INTEGER,                 INTENT(IN)    :: QC       ! are coefs valid?
    TYPE(StateIndex_t),      INTENT(IN)    :: IGRin    ! indices for linear retrieval
    TYPE(StateIndex_t),      INTENT(IN)    :: NGRin    ! sizes for linear retrieval
    CHARACTER(LEN=mxCoordTyp), INTENT(IN)  :: vCoordTyp ! vertical coordinate type
    REAL,    DIMENSION(:),   INTENT(IN)    :: pref     ! pressure grid; used only
                                                          ! if retrieving Tskin and T(p)
    REAL,                    INTENT(IN)    :: pSfc     ! surface pressure; used only
                                                          ! if retrieving Tskin and T(p)
    REAL,                    INTENT(IN)    :: Tsfc     ! sfc air temperature; used only
                                                          ! if retrieving Tskin, not T(p)
    INTEGER,                 INTENT(IN)    :: iH2O     ! where is H2O; used only
                                                          ! if retrieving H2O
    INTEGER,                 INTENT(IN)    :: nMembers ! # cases with common coef
    REAL,    DIMENSION(:,:), INTENT(IN)    :: radGroup ! radiances
    REAL,    DIMENSION(:,:), INTENT(INOUT) :: xGroup   ! retrievals

    ! Local variables:
    REAL, DIMENSION(SIZE(xOffG)) :: xTmp
    INTEGER :: n,nLev
    REAL    :: xSfc

    nLev=SIZE(pref)

    xGroup=MISSING_REAL  ! Overridden if rad data are valid
    IF (QC /= 0) RETURN

    DO n=1,nMembers

      xTmp=xOffG+MATMUL(xGainG,radGroup(:,n))

      ! Convert delta-Tskin to Tskin
      !   For retrieval in sigma coordinates, xOff and xGain can be reformulated
      !   to account for this without doing it at this stage
      IF ((TRIM(vCoordTyp) == Pcoord) .AND. (NGRin%Tskin > 0)) THEN
        IF (NGRin%temp > 0) THEN   ! must update Tsfc, if T(p) is retrieved
          CALL lvl_int(xTmp(IGRin%temp:IGRin%temp+NGRin%temp-1),pref,nLev, &
             pSfc,xSfc)
        ELSE
          xSfc=Tsfc
        ENDIF
        xTmp(IGRin%TSkin)=xTmp(IGRin%TSkin)+xSfc
      ENDIF

      ! Convert water vapor from log to linear
      IF (NGRin%mol(iH2O) > 0) &
        xTmp(IGRin%mol(iH2O):IGRin%mol(iH2O)+NGRin%mol(iH2O)-1)= &
             EXP(xTmp(IGRin%mol(iH2O):IGRin%mol(iH2O)+NGRin%mol(iH2O)-1))

      xGroup(:,n)=xTmp

    ENDDO

    RETURN

  END SUBROUTINE linRetr

  SUBROUTINE getLinInvData(IGout, NGout, nParGon)
!<f90Subroutine>********************************************************
!
! NAME:
!
!   getLinInvData
!
! PURPOSE:
!
!   Returns module global data
!
! SYNTAX:
!
!   CALL getLinInvData(NGout,IGout,nParGon)
!
! ARGUMENTS:
!
!   INPUTS/OUTPUTS:
!
!   IGRout   TYPE(STATEINDEX_T)  Starting indices for retrieved sections of
!                                geophysical state vector
!   NGRout   TYPE(STATEINDEX_T)  Number of elements for retrieved sections of
!                                geophysical state vector
!   nParGon  INTEGER             Total number of retrieved elements in
!                                geophysical state vector
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
    TYPE(StateIndex_t), INTENT(INOUT) :: IGout  ! indices for linear retrieval
    TYPE(StateIndex_t), INTENT(INOUT) :: NGout  ! sizes for linear retrieval
    INTEGER,            INTENT(INOUT) :: nParGon ! number of geo-space parameters to retrieve

    !Local
    character (len=mxLength) :: procName
    logical :: dbg = .true.

    procName = '[LinInver::getLinInvData]:'
    if (dbg) print *, trim(procName)//' starting...'
    IGout = IGR
    NGout = NGR
    nParGon=getVectorLength(NGR)
    if (dbg) print *, trim(procName)//' ending...'

  END SUBROUTINE getLinInvData

  SUBROUTINE destroyLinInvert()

!<f90Subroutine>********************************************************
!
! NAME:
!
!   destroyLinInvert
!
! PURPOSE:
!
!   Deallocate arrays
!
! SYNTAX:
!
!   CALL destroyLinInvert()
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

  ! Release memory

    DEALLOCATE (mapGeo,uCompr,xBakCompr)

    RETURN

  END SUBROUTINE destroyLinInvert

END MODULE LinInvert
