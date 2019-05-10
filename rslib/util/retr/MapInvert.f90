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
!   Copyright AER, Inc 2001-2009, All Right Reserved
!   See the file README-DATARIGHTS.txt included with this release
!   for additional details.
!
!*************************************************************</f90File>

!--------------------------------------------------------------
!
!  MODULE MapInvert: contains subroutines needed for the
!                retrieval mapping. AER Inc. 2004.!
!--------------------------------------------------------------
MODULE MapInvert

! <f90Module>***********************************************************
!
! NAME:
!
!   MapInvert
!
! PURPOSE:
!
!   Contains subroutines for transforming between retrieval vector space and
!   geophysical vector space.
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

  USE StateIndexModule
  USE LvlInterp
  USE VertCoord, ONLY: mxCoordTyp,Pcoord
  IMPLICIT NONE
  PRIVATE
  !------------------------------------------------------------------------
  !	Public items made available to calling program(s)
  !------------------------------------------------------------------------
  PUBLIC :: map_retr2geo,map_geo2retr,mapK_geo2retr,set_mw_invert
  PUBLIC :: MOL_TRAN_NONE,MOL_TRAN_LOG

  !------------------------------------------------------------------------
  !	Constants
  !------------------------------------------------------------------------
  INTEGER, PARAMETER :: MOL_TRAN_NONE = 0
  INTEGER, PARAMETER :: MOL_TRAN_LOG = 1

CONTAINS

  SUBROUTINE map_retr2geo(x,xBakG,xG,EmMw,xSfc,vmtx,NParG,NPar,&
       IG,NG,vCoordTyp,pref,nlev,nchan,IEmMwG,ih2o,nmol,molTran, &
       xOutRangeFlag,logSize)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   map_retr2geo
!
! PURPOSE:
!
!   Convert state vector from retrieval to geophysical space.
!
! SYNTAX:
!
!   CALL map_retr2geo(x, xBakG, xG, EmMw, xSfc, vmtx, NParG,
!      NPar, IG, NG, pref, nlev, nchan, IEmMwG, ih2o, nmol, molTran,
!      xOutRangeFlag, logSize)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   x        REAL                Retrieval state vector
!   xBakG    REAL                Background state vector in geophysical space
!   vmtx     REAL                Transformation matrix between geophysical
!                                and retrieval spaces
!   NParG    INTEGER             Total number of elements in geophysical
!                                state vector
!   NPar     INTEGER             Total number of elements in retrieval vector
!   IG       TYPE(STATEINDEX_T)  Starting indices for sections of geophysical
!                                state vector
!   NG       TYPE(STATEINDEX_T)  Number of elements for sections of
!                                geophysical state vector
!   pref     REAL                Pressure on atmospheric grid levels
!   nlev     INTEGER             Number of atmospheric levels
!   nchan    INTEGER             Number of channels
!   IEmMwG   INTEGER             Starting index for MW emissivity in
!                                geophysical state vector
!   ih2o     INTEGER             Index for water vapor in molecules part of
!                                state vector
!   nmol     INTEGER             Number of molecular species
!   logSize* LOGICAL             Flag for type of cloud size parameter
!                                retrieval: Logarithmic vs Linear
!
!   INPUTS/OUTPUTS:
!
!   xG       REAL                State vector in geophysical space
!   EmMw     REAL                MW emissivity
!   xSfc     REAL                Surface-level air temperature
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !---Input variables
    REAL, DIMENSION(:),   INTENT(IN)           :: x,xBakG,pref
    REAL, DIMENSION(:,:), INTENT(IN)           :: vmtx
    TYPE(StateIndex_t),   INTENT(IN)           :: IG,NG
    CHARACTER(LEN=mxCoordTyp),INTENT(IN)       :: vCoordTyp
    INTEGER,              INTENT(IN)           :: NparG,Npar,nlev
    INTEGER,              INTENT(IN)           :: nchan,IEmMwG,ih2o,nmol
    INTEGER, DIMENSION(:),INTENT(IN)           :: molTran
    LOGICAL,              INTENT(IN), OPTIONAL :: logSize
    !---Output variables
    REAL, DIMENSION(:),   INTENT(INOUT)        :: xG,EmMw
    REAL,                 INTENT(INOUT)        :: xSfc
    LOGICAL,              INTENT(INOUT)        :: xOutRangeFlag
    !---Local variables
    INTEGER              :: imol

    xG(1:NparG)=MATMUL(vmtx(1:NparG,1:Npar),x(1:Npar))
    !-- Remove the out-of-bound values from ill-conditioned inversion
    !-- to prevent from the exponential overfloating in moisture calculation.
    xOutRangeFlag=.FALSE.
    IF (MAXVAL(ABS(xG(1:NG%Temp))) > 100.) xOutRangeFlag=.TRUE.
    IF ((molTran(ih2o) == MOL_TRAN_LOG) .AND. &
         MAXVAL(ABS(xG(IG%mol(iH2o):IG%mol(iH2o)+NG%mol(iH2o)-1))) > 10.) &
         xOutRangeFlag=.TRUE.
    IF (xOutRangeFlag) THEN
       xG(1:NparG) = xBakG(1:NparG)
    ELSE
       xG(1:NparG) = xG(1:NparG) + xBakG(1:NparG)
    ENDIF
    IF (TRIM(vCoordTyp) == Pcoord) THEN
       CALL lvl_int(xG,pref,Nlev,xG(IG%Psfc),xSfc)
       xG(IG%TSkin)=xG(IG%TSkin)+xSfc
    ELSE
       xSfc=xG(IG%temp+NG%temp-1)
    ENDIF
    DO imol=1,nmol
       IF (molTran(imol) == MOL_TRAN_LOG) &
          xG(IG%mol(imol):IG%mol(imol)+NG%mol(imol)-1)=&
            EXP(xG(IG%mol(imol):IG%mol(imol)+NG%mol(imol)-1))
    ENDDO
    IF (present(logSize)) THEN
       IF (logSize) THEN
          IF (NG%cldLiq > 3) xG(IG%cldLiq+3)=EXP(xG(IG%cldLiq+3))
          IF (NG%cldIce > 3) xG(IG%cldIce+3)=EXP(xG(IG%cldIce+3))
       ENDIF
    ENDIF
    if (nchan > 0) EmMw(1:nchan)=xG(IEmMwG:IEmMwG+nchan-1)
    RETURN
  END SUBROUTINE map_retr2geo

  SUBROUTINE map_geo2retr(xG,xBakG,Tsfc,x,umtx,NparG,Npar,IG,NG,nmol,molTran, &
     vCoordTyp,logSize)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   map_geo2retr
!
! PURPOSE:
!
!   Convert state vector from geophysical to retrieval space.
!
! SYNTAX:
!
!   CALL map_geo2retr(xG, xBakG, Tsfc, x, umtx, NparG, Npar, IG, NG,
!      nmol, molTran, logSize)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   xG       REAL                State vector in geophysical space
!   xBakG    REAL                Background state vector in geophysical space
!   Tsfc     REAL                Surface-level air temperature
!   umtx     REAL                Transformation matrix between geophysical
!                                and retrieval spaces
!   NparG    INTEGER             Total number of elements in geophysical
!                                state vector
!   Npar     INTEGER             Total number of elements in retrieval vector
!   IG       TYPE(STATEINDEX_T)  Starting indices for sections of geophysical
!                                state vector
!   NG       TYPE(STATEINDEX_T)  Number of elements for sections of
!                                geophysical state vector
!   logSize* LOGICAL             Flag for type of cloud size parameter
!                                retrieval: Logarithmic vs Linear
!
!   INPUTS/OUTPUTS:
!
!   x        REAL                Retrieval state vector
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !---Input variables
    !---NparG: #parameters in geophysical space
    !---Npar : #parameters in retrieval space
    REAL, DIMENSION(:),   INTENT(IN)           :: xG,xBakG
    REAL, DIMENSION(:,:), INTENT(IN)           :: umtx
    REAL,                 INTENT(IN)           :: Tsfc
    INTEGER,              INTENT(IN)           :: NparG,Npar,nmol
    INTEGER, DIMENSION(:),INTENT(IN)           :: molTran
    CHARACTER(LEN=mxCoordTyp),INTENT(IN)       :: vCoordTyp
    LOGICAL,              INTENT(IN), OPTIONAL :: logSize
    TYPE(StateIndex_t),   INTENT(IN)           :: IG,NG
    !---Output variables
    REAL, DIMENSION(:),   INTENT(INOUT)        :: x
    !---Local variables
    REAL, DIMENSION(SIZE(xG)) :: xGTmp,xGdev
    INTEGER                   :: i
    INTEGER                   :: imol

    !----Transform variables
    xGTmp(1:NparG)=xG(1:NparG)
    IF (TRIM(vCoordTyp) == Pcoord) &
       xGTmp(IG%TSkin)=xG(IG%TSkin)-TSfc
    DO imol=1,nmol
       IF (molTran(imol) == MOL_TRAN_LOG) &
          xGTmp(IG%mol(imol):IG%mol(imol)+NG%mol(imol)-1)= &
          ALOG(xG(IG%mol(imol):IG%mol(imol)+NG%mol(imol)-1))
    ENDDO
    IF (present(logSize)) THEN
       IF (logSize) THEN
          IF (NG%cldLiq > 3) xGTmp(IG%cldLiq+3)=ALOG(xG(IG%cldLiq+3))
          IF (NG%cldIce > 3) xGTmp(IG%cldIce+3)=ALOG(xG(IG%cldIce+3))
       END IF
    END IF
    !----Compute deviation from background
    xGdev(1:NparG)=xGTmp(1:NparG)-xBakG(1:NparG)
    !----Apply basis function transformations
    x(1:Npar)=MATMUL(umtx(1:Npar,1:NparG),xGdev(1:NparG))
    RETURN
  END SUBROUTINE map_geo2retr


  !----------------------------------------------------------------------------
  ! This subroutine handles the conversion of dR/dw to dR/dlogw,
  ! and the computation of Tskin derivative without any compression,
  ! followed by conversion of xktTmp into EOF domain.
  ! This routine is intended to be used
  ! in MapInvert::set_MW_Invert()
  !----------------------------------------------------------------------------
  SUBROUTINE mapK_geo2retr(xktTmp,xkt,vmtx,IG,NG,xG,vCoordTyp,pref,Npar,NparG, &
             nch,nmol,molTran,logSize)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   mapK_geo2retr
!
! PURPOSE:
!
!   Convert Jacobians from geophysical to retrieval space.
!
! SYNTAX:
!
!   CALL mapK_geo2retr(xktTmp, xkt, vmtx, IG, NG, xG, pref, Npar,
!      NparG, nch, nmol, molTran, logSize)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   vmtx     REAL                Transformation matrix from retrieval space to
!                                geophysical space
!   IG       TYPE(STATEINDEX_T)  Starting indices for sections of geophysical
!                                state vector
!   NG       TYPE(STATEINDEX_T)  Number of elements for sections of
!                                geophysical state vector
!   xG       REAL                State vector in geophysical space
!   pref     REAL                Pressure on atmospheric grid levels
!   Npar     INTEGER             Total number of elements in retrieval vector
!   NparG    INTEGER             Total number of elements in geophysical
!                                state vector
!   nch      INTEGER             Number of channels turned on
!   logSize* LOGICAL             Flag for type of cloud size parameter
!                                retrieval: Logarithmic vs Linear
!
!   INPUTS/OUTPUTS:
!
!   xktTmp   REAL                Jacobian matrix transpose in geophysical
!                                space
!   xkt      REAL                Jacobian matrix transpose
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

     !---- input variables
     !--- nch represents the total number of compressed channels
     REAL, DIMENSION(:),    INTENT(IN)           :: xG
     REAL, DIMENSION(:),    INTENT(IN)           :: pref
     TYPE(StateIndex_t),    INTENT(IN)           :: IG, NG
     CHARACTER(LEN=mxCoordTyp),INTENT(IN)        :: vCoordTyp
     INTEGER,               INTENT(IN)           :: nch
     INTEGER,               INTENT(IN)           :: Npar, NparG
     INTEGER,               INTENT(IN)           :: nmol
     INTEGER, DIMENSION(:), INTENT(IN)           :: molTran
     LOGICAL,               INTENT(IN), OPTIONAL :: logSize
     REAL, DIMENSION(:, :), INTENT(IN)           :: vmtx
     !---- input/output variables
     REAL, DIMENSION(:, :), INTENT(INOUT)        :: xktTmp
     !--- output variables
     REAL, DIMENSION(:,:),  INTENT(INOUT)        :: xkt
     !--- local variables
     INTEGER               :: i, j
     INTEGER               :: imol
     INTEGER               :: Nsurf
     REAL                  :: xSfc,dTsdTupper,dTsdTlower

     xkt  = 0

     IF (TRIM(vCoordTyp) == Pcoord) THEN

        CALL lvl_int(xG,pref,NG%temp,xG(IG%Psfc),xSfc,ilvl=Nsurf, &
           dxlvldxu=dTsdTupper,dxlvldxl=dTsdTlower)

        DO j = 1, nch
           !---Treat Tskin derivative as Tskin-Tsfc
           xktTmp(Nsurf-1,j)=xktTmp(Nsurf-1,j)+xktTmp(IG%Tskin,j)*dTsdTupper
           xktTmp(Nsurf,  j)=xktTmp(Nsurf,  j)+xktTmp(IG%Tskin,j)*dTsdTlower
        END DO

     ENDIF

     DO j = 1, nch

        !---Convert dR/dw to dR/dlogw: No compression
        DO imol=1,nmol
           IF (molTran(imol) == MOL_TRAN_LOG) THEN
              DO i = IG%mol(imol), IG%mol(imol) + NG%mol(imol) - 1
                   xktTmp(i, j)  = xktTmp(i, j) * xG(i)
              END DO
           ENDIF
        END DO
     END DO

     !---For log-Deff's retrieval, transform derivatives
     IF (PRESENT(logSize)) THEN
        IF (logSize) THEN
           xktTmp(IG%cldLiq+3,1:nch) = xktTmp(IG%cldLiq+3,1:nch) * &
                xG(IG%cldLiq+3)
           xktTmp(IG%cldIce+3,1:nch) = xktTmp(IG%cldIce+3,1:nch) * &
                xG(IG%cldIce+3)
        END IF
     END IF

     ! Convert xktTmp into EOF domain
     xkt(1:Npar,1:nch) = MATMUL(TRANSPOSE(vmtx(1:NparG,1:Npar)), &
                        xktTmp(1:NparG, 1:nch))
     RETURN
  END SUBROUTINE mapK_geo2retr


  !---- This routine will be used in both MW and IR retrieval
  SUBROUTINE set_MW_Invert(xG,Ym,Y,xktG,xkEmMw,rerr,xkt,Sy,delY,nch,itermw,&
       kchan,vmtx,nchan,dradChan,alpha1slope,alpha1int,NparG,Npar,IG,NG,&
       iemmwG,nmol,molTran,vCoordTyp,pref,logSize)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   set_MW_Invert
!
! PURPOSE:
!
!   Prepare for inversion. Merge emissivity Jacobians into complete Jacobian
!   matrix. Convert radiance to the difference form. Tune measurement/model
!   error. Exclude unused channels.
!
! SYNTAX:
!
!   CALL set_MW_Invert(xG, Ym, Y, xktG, xkEmMw, rerr, xkt, Sy, delY,
!      nch, itermw, kchan, vmtx, nchan, dradChan, alpha1slope,
!      alpha1int, NparG, Npar, IG, NG, iemmwG, nmol, molTran, pref, logSize)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   xG           REAL                State vector in geophysical space
!   Ym           REAL                Radiometric measurements
!   Y            REAL                Radiometric data computed from state
!                                    vector
!   xkEmMw       REAL                Jacobians for MW emissivity
!   rerr         REAL                Instrument noise standard deviation
!   itermw       INTEGER             Iteration number
!   kchan        LOGICAL             Channel on/off mask
!   vmtx         REAL                Transformation matrix from retrieval space to
!                                    geophysical space
!   nchan        INTEGER             Number of channels
!   dradChan     REAL                Noise adjustment to stabilize
!                                    convergence
!   alpha1slope  REAL                Slope of factor to model linearization
!                                    error
!   alpha1int    REAL                Intercept of factor to model
!                                    linearization error
!   NparG        INTEGER             Total number of elements in geophysical
!                                    state vector
!   Npar         INTEGER             Total number of elements in retrieval
!                                    vector
!   IG           TYPE(STATEINDEX_T)  Starting indices for sections of
!                                    geophysical state vector
!   NG           TYPE(STATEINDEX_T)  Number of elements for sections of
!                                    geophysical state vector
!   iemmwG       INTEGER             Starting index for MW emissivity in
!                                    geophysical state vector
!   ih2o         INTEGER             Index for water vapor in molecules part
!                                    of state vector
!   pref         REAL                Pressure on atmospheric grid levels
!   logSize*     LOGICAL             Flag for type of cloud size parameter
!                                    retrieval: Logarithmic vs Linear
!
!   INPUTS/OUTPUTS:
!
!   xktG         REAL                Jacobian matrix transpose in geophysical
!                                    space
!   xkt          REAL                Jacobian matrix transpose
!   Sy           REAL                Measurement error variance
!   delY         REAL                Radiometric measurements minus computed
!   nch          INTEGER             Number of channels turned on
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !----Input/output variables
    REAL,    DIMENSION(:),   INTENT(IN)           :: xG,Ym,Y,rerr,xkEmMw
    REAL,    DIMENSION(:),   INTENT(INOUT)        :: Sy,delY
    REAL,    DIMENSION(:),   INTENT(IN)           :: dradChan
    REAL,    DIMENSION(:,:), INTENT(INOUT)        :: xktG
    REAL,    DIMENSION(:,:), INTENT(INOUT)        :: xkt
    REAL,    DIMENSION(:,:), INTENT(IN)           :: vmtx
    REAL,                    INTENT(IN)           :: alpha1slope,alpha1int
    INTEGER,                 INTENT(INOUT)        :: nch
    INTEGER,                 INTENT(IN)           :: itermw,NparG,Npar,nchan
    INTEGER,                 INTENT(IN)           :: iemmwG
    INTEGER,                 INTENT(IN)           :: nmol
    INTEGER, DIMENSION(:),   INTENT(IN)           :: molTran
    LOGICAL, DIMENSION(:),   INTENT(IN)           :: kchan
    TYPE(StateIndex_t),      INTENT(IN)           :: IG,NG
    CHARACTER(LEN=mxCoordTyp),INTENT(IN)          :: vCoordTyp
    REAL,    DIMENSION(:),   INTENT(IN)           :: pref
    LOGICAL,                 INTENT(IN), OPTIONAL :: logSize
    !---Local variables
    REAL,    DIMENSION(SIZE(xG),SIZE(Y)) :: xktTmp
    REAL                                 :: rerr2,alpha1
    INTEGER                              :: i,j,nchtmp

    nch=0
    xktTmp=0.

    DO j = 1,nchan
       IF(kchan(j))THEN
          nch=nch+1
          !---calculated Error covariance matix Sy:
          delY(nch)              = Ym(j)-Y(j)
          rerr2                  = rerr(j)**2+dradChan(j)**2
          alpha1                 = float(itermw)/alpha1slope+alpha1int
          Sy(nch)                = MAX(rerr2,(delY(nch)/alpha1)**2)
          !---setting up derivatives:
          xktTmp(1:NparG,nch)   = xktG(1:NparG,j)
          !---Map nchmw mw emis. into NEmMwG of Mw frq:
          xktTmp(IEmMwG+j-1,nch) = xkEmMw(j)
          xktG(IEmMwG+j-1,nch)   = xkEmMw(j)
       ENDIF
    END DO

    CALL mapK_geo2retr(xktTmp,xkt,vmtx,IG,NG,xG,vCoordTyp,pref,Npar,NparG,nch, &
            nmol,molTran,logSize)
    RETURN
  END SUBROUTINE set_MW_Invert

END MODULE MapInvert
