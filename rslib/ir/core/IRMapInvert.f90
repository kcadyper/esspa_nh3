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

!--------------------------------------------------------------
!
!  MODULE IRMapInvert: contains subroutines needed for the
!                algorithm retrieval mapping. AER Inc. 2004.!
!--------------------------------------------------------------
MODULE IRMapInvert

! <f90Module>***********************************************************
!
! NAME:
!
!   IRMapInvert
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
  USE MapInvert, ONLY : map_retr2geo, mapK_geo2retr, MOL_TRAN_LOG
  USE VertCoord, ONLY: mxCoordTyp,Pcoord
  IMPLICIT NONE
  PRIVATE
  !------------------------------------------------------------------------
  !	Public items made available to calling program(s)
  !------------------------------------------------------------------------
  PUBLIC :: IRmap_retr2geo,map_regr2retr_geo, &
       set_irmw_invert

CONTAINS

  SUBROUTINE IRmap_retr2geo(x,xBakG,xG,EmMw,xSfc,vmtx,NParG,NPar,&
       IG,NG,vCoordTyp,pref,nlev,nchan,IEmMwG,ih2o,nmol,molTran,IEmIRG,nemirG, &
       EmRf,xOutRangeFlag,logSize)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   IRmap_retr2geo
!
! PURPOSE:
!
!   Convert state vector from retrieval to geophysical space.
!
! SYNTAX:
!
!   CALL IRmap_retr2geo(x, xBakG, xG, EmMw, xSfc, vmtx, NParG, NPar,
!      IG, NG, pref, nlev, nchan, IEmMwG, ih2o, nmol, molTran, IEmIRG, nemirG, EmRf,
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
!   IEmIRG   INTEGER             Starting index for IR emissivity hinge
!                                points
!   nemirG   INTEGER             Number of IR emissivity hinge points
!   logSize  LOGICAL             Flag for type of cloud size parameter
!                                retrieval: Logarithmic vs Linear
!
!   INPUTS/OUTPUTS:
!
!   xG       REAL                State vector in geophysical space
!   EmMw     REAL                MW emissivity
!   xSfc     REAL                Surface-level air temperature
!   EmRf     REAL                IR emissivity/reflectivity at hinge points
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !---Input variables
    REAL, DIMENSION(:),   INTENT(IN)   :: x,xBakG
    REAL, DIMENSION(:),   INTENT(INOUT):: xG,EmMw
    REAL,                 INTENT(INOUT):: xSfc
    REAL, DIMENSION(:,:), INTENT(IN)   :: vmtx
    INTEGER,              INTENT(IN)   :: NparG,Npar
    TYPE(StateIndex_t),   INTENT(IN)   :: IG,NG
    CHARACTER(LEN=mxCoordTyp),INTENT(IN):: vCoordTyp
    REAL, DIMENSION(:),   INTENT(IN)   :: pref
    INTEGER,              INTENT(IN)   :: nlev,nchan,IEmMwG,ih2o,nmol
    INTEGER, DIMENSION(:),INTENT(IN)   :: molTran
    INTEGER,              INTENT(IN)   :: IEmIRG,nemirG
    REAL, DIMENSION(:),   INTENT(INOUT):: EmRf
    LOGICAL,              INTENT(INOUT):: xOutRangeFlag
    LOGICAL,              INTENT(IN)   :: logSize
    !---Local variables
    integer              :: i,ii
    INTEGER              :: imol

    call map_retr2geo(x,xBakG,xG,EmMw,xSfc,vmtx,NParG,NPar,&
       IG,NG,vCoordTyp,pref,nlev,nchan,IEmMwG,ih2o,nmol,molTran, &
       xOutRangeFlag,logSize)
    do i=1,nemirG
       ii=2*(i-1)+1
       emrf(ii)=xG(IEmIRG+i-1)
       emrf(ii+1)=xG(IEmIRG+nemirG+i-1)
    enddo
    RETURN
  END SUBROUTINE IRmap_retr2geo


  SUBROUTINE map_regr2retr_geo(xregr,npargregr,x,xBakG,xG,EmMw,&
       xSfc,umtx,NparG,Npar,&
       IG,NG,pref,nlev,nchan,IEmMwG,nmol,molTran,IEmIRG,nemirG,EmRf)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   map_regr2retr_geo
!
! PURPOSE:
!
!   For regression: convert state vector from retrieval to geophysical space.
!
! SYNTAX:
!
!   CALL map_regr2retr_geo(xregr, npargregr, x, xBakG, xG, EmMw,
!      xSfc, umtx, NparG, Npar, IG, NG, pref, nlev, nchan, IEmMwG,
!      nmol, molTran, IEmIRG, nemirG, EmRf)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   xregr      REAL                State vector for regression
!   npargregr  INTEGER             Total number of elements in retrieval
!                                  vector
!   xBakG      REAL                Background state vector in geophysical
!                                  space
!   umtx       REAL                Transformation matrix between geophysical
!                                  and retrieval spaces
!   NparG      INTEGER             Total number of elements in geophysical
!                                  state vector
!   Npar       INTEGER             Total number of elements in retrieval
!                                  vector
!   IG         TYPE(STATEINDEX_T)  Starting indices for sections of
!                                  geophysical state vector
!   NG         TYPE(STATEINDEX_T)  Number of elements for sections of
!                                  geophysical state vector
!   pref       REAL                Pressure on atmospheric grid levels
!   nlev       INTEGER             Number of atmospheric levels
!   nchan      INTEGER             Number of channels
!   IEmMwG     INTEGER             Starting index for MW emissivity in
!                                  geophysical state vector
!   IEmIRG     INTEGER             Starting index for IR emissivity hinge
!                                  points
!   nemirG     INTEGER             Number of IR emissivity hinge points
!
!   INPUTS/OUTPUTS:
!
!   x          REAL                Retrieval state vector
!   xG         REAL                State vector in geophysical space
!   EmMw       REAL                MW emissivity
!   xSfc       REAL                Surface-level air temperature
!   EmRf       REAL                IR emissivity/reflectivity at hinge points
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !---Input variables
    REAL, DIMENSION(:),   INTENT(IN)   :: xregr
    INTEGER,              INTENT(IN)   :: npargregr
    REAL, DIMENSION(:),   INTENT(INOUT):: x
    REAL, DIMENSION(:),   INTENT(IN)   :: xBakG
    REAL, DIMENSION(:),   INTENT(INOUT):: xG,EmMw
    REAL,                 INTENT(INOUT):: xSfc
    REAL, DIMENSION(:,:), INTENT(IN)   :: umtx
    INTEGER,              INTENT(IN)   :: NparG,Npar
    TYPE(StateIndex_t),   INTENT(IN)   :: IG,NG
    REAL, DIMENSION(:),   INTENT(IN)   :: pref
    INTEGER,              INTENT(IN)   :: nlev,nchan,IEmMwG,nmol
    INTEGER, DIMENSION(:),INTENT(IN)   :: molTran
    INTEGER,              INTENT(IN)   :: IEmIRG,nemirG
    REAL, DIMENSION(:),   INTENT(INOUT):: EmRf
    !---Local variables
    integer              :: i,ii
    INTEGER              :: imol
    REAL, DIMENSION(SIZE(xG)) :: xGTmp,xGdev
    REAL                 :: tsfc

    !----convert xG from geophy. into retrieval domain:
    CALL lvl_int(xregr,pref,Nlev,xregr(IG%Psfc),xSfc)
    TSfc=XSfc
    xGTmp=0.
    xGTmp(1:NparGRegr)=xregr(1:NparGRegr)
    xGTmp(IG%TSkin)=xregr(IG%TSkin)-TSfc
    DO imol=1,nmol
       IF (molTran(imol) == MOL_TRAN_LOG) &
          xGTmp(IG%mol(imol):IG%mol(imol)+NG%mol(imol)-1)= &
          ALOG(xregr(IG%mol(imol):IG%mol(imol)+NG%mol(imol)-1))
    ENDDO
    xGdev(1:NparGRegr)=xGTmp(1:NparGRegr)-xBakG(1:NparGRegr)
    x(1:Npar)=MATMUL(umtx(1:Npar,1:NparG),xGdev(1:NparG))

    xG(1:NparGRegr)=xregr(1:NparGRegr)
    xG((NparGRegr+1):NParG)=xBakG((NparGRegr+1):NParG)

    EmMw(1:nchan)=xG(IEmMwG:IEmMwG+nchan-1)
    do i=1,nemirG
       ii=2*(i-1)+1
       emrf(ii)=xG(IEmIRG+i-1)
       emrf(ii+1)=xG(IEmIRG+nemirG + i - 1)
    enddo
    RETURN
  END SUBROUTINE Map_regr2retr_geo


  subroutine set_IRMW_Invert(xG,Ym,Y,xktG,xkEmMw,rerr,AtmNoise,dradChan, &
       xkt,xkEmRf,Sy,delY,nch,iter,nmol,molTran,kchan,vmtx,nparG,nPar, &
       nchanmw,nchan,IG,NG,IEmMwG,IEmIrG,nemirG,vCoordTyp,pref,Airs_Err,logSize,&
       chanDataCompressed)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   set_IRMW_Invert
!
! PURPOSE:
!
!   Prepare for inversion. Merge emissivity Jacobians into complete Jacobian
!   matrix. Convert radiance to the difference form. Tune measurement/model
!   error. Exclude unused channels.
!
! SYNTAX:
!
!   CALL set_IRMW_Invert(xG, Ym, Y, xktG, xkEmMw, rerr,
!      AtmNoise, dradChan, xkt, xkEmRf, Sy, delY, nch, iter, nmol, molTran,
!      kchan, vmtx, nparG, nPar, nchanmw, nchan, IG, NG, IEmMwG,
!      IEmIrG, nemirG, pref, Airs_Err, logSize)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   xBakG     REAL                Background state vector in geophysical
!                                 space
!   Ym        REAL                Radiometric measurements
!   Y         REAL                Radiometric data computed from state vector
!                                 they are compressed
!                         if idx(1:nch) fpor all kchan(idx)==.true.
!                         then compression means that y(j) corresponds to Ym(idx(j))
!                         or kcnah(idx(j))=.true.
!                         the same for Jacobians
!
!   xkEmMw    REAL                Jacobians for MW emissivity
!   rerr      REAL                Instrument noise standard deviation
!   AtmNoise  REAL                Atmospheric model error standard deviation
!   dradChan  REAL                Noise adjustment to stabilize convergence
!   xkEmRf    REAL                Jacobians for IR emissivity/reflectivity
!   iter      INTEGER             Iteration number
!   kchan     logical             Channel on/off mask
!   vmtx      REAL                Transformation matrix from retrieval space to
!                                 geophysical space
!   nparG     INTEGER             Total number of elements in geophysical
!                                 state vector
!   nPar      INTEGER             Total number of elements in retrieval
!                                 vector
!   nchanmw   INTEGER             Number of MW channels
!   nchan     INTEGER             Number of channels
!   IG        TYPE(STATEINDEX_T)  Starting indices for sections of
!                                 geophysical state vector
!   NG        TYPE(STATEINDEX_T)  Number of elements for sections of
!                                 geophysical state vector
!   IEmMwG    INTEGER             Starting index for MW emissivity in
!                                 geophysical state vector
!   IEmIrG    INTEGER             Starting index for IR emissivity hinge
!                                 points
!   nemirG    INTEGER             Number of IR emissivity hinge points
!   pref      REAL                Pressure on atmospheric grid levels
!   logSize   LOGICAL             Flag for type of cloud size parameter
!                                 retrieval: Logarithmic vs Linear
!   chanDataCompressed  logical      .true. -  Y (computed radiometric) is compressed
!                                 .false.-  Y (computed radiometric) is uncompressed
!
!   INPUTS/OUTPUTS:
!
!   xG        REAL                State vector in geophysical space
!   xktG      REAL                Jacobian matrix transpose in geophysical
!                                 space
!   xkt       REAL                Jacobian matrix transpose
!   Sy        REAL                Measurement error variance
!   delY      REAL                Radiometric measurements minus computed
!   nch       INTEGER             Number of channels turned on
!   Airs_Err  REAL                Net error standard deviation, special case
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>


    real,    dimension(:),                      intent(inout):: xG
    real,    dimension(:),                      intent(in)   :: Ym,Y
    real,    dimension(:,:),                    intent(inout):: xktG
    real,    dimension(:),                      intent(in)   :: xkEmMw
    real,    dimension(:),                      intent(in)   :: rerr
    real,    dimension(:),                      intent(in)   :: AtmNoise,dradChan
    real,    dimension(:,:),                    intent(inout):: xkt
    real,    dimension(2,nemirG,nchan-nchanmw), intent(in)   :: xkEmRf
    real,    dimension(:),                      intent(inout):: Sy,delY
    integer,                                    intent(inout):: nch
    integer,                                    intent(in)   :: iter
    integer,                                    intent(in)   :: nmol
    integer, dimension(:),                      intent(in)   :: molTran
    logical, dimension(:),                      intent(in)   :: kchan
    real,    dimension(:,:),                    intent(in)   :: vmtx
    integer,                                    intent(in)   :: nparG,nPar
    integer,                                    intent(in)   :: nchanmw,nchan
    TYPE(StateIndex_t),                         intent(in)   :: IG,NG
    integer,                                    intent(in)   :: IEmMwG,IemIRG
    integer,                                    intent(in)   :: nemirG
    CHARACTER(LEN=mxCoordTyp),                  intent(in)   :: vCoordTyp
    real,    dimension(:),                      intent(in)   :: pref
    real,    dimension(:),                      intent(inout):: Airs_err
    logical,                                    intent(in)   :: logSize
    logical, optional,                          intent(in)   :: chanDataCompressed

    ! local variables:

    REAL,    DIMENSION(SIZE(xG),SIZE(Y)) :: xktTmp
    REAL,    DIMENSION(SIZE(Y)) :: rerrAtm
    integer :: i,j,k, jch
    real    :: alpha1,rerr2
    integer :: numCompChIR, numCompChMW, numCompCh
    logical :: chanDataCompressedLoc

    ! iatmnoise = 1 to add atmnoise
    integer :: iatmnoise=1

    ! - Combine emissivity and reflectivity Jacobians with
    !   state vector ones
    !---- Surface emissivity and reflectivity
    if (present(chanDataCompressed)) then
      chanDataCompressedLoc=chanDataCompressed
    else
      chanDataCompressedLoc=.false.
    end if

    numCompChMW=COUNT(kchan(1:nchanmw))
    numCompChIR=COUNT(kchan(nchanmw+1:nchan))
    numCompCh=numCompChIR+numCompChMW

    if (chanDataCompressedLoc) then
      ! add emissivity Jacobian into correct position(Jacobian is compressed)
      jch=0
      do j=numCompChMW+1,numCompCh
         jch=jch+1
         xktG(IEmIrG:IEmIrG+nemirG-1,j) = xkEmRf(1, 1:nemirG, jch)
         xktG(IEmIrG+nemirG:IEmIrG+2*nemirG-1,j)=xkEmRf(2,1:nemirG,jch)
      enddo
    else
      jch=0
      do j=nchanmw+1,nchan
         jch=jch+1
         xktG(IEmIrG:IEmIrG+nemirG-1,j) = xkEmRf(1, 1:nemirG, jch)
         xktG(IEmIrG+nemirG:IEmIrG+2*nemirG-1,j)=xkEmRf(2,1:nemirG,jch)
      end do
    end if

    !---- Add Atmospheric Noise to rerr if desired
    ! (note that it is already part of the MW portion of rerr)

    if (numCompChMW > 0) then
       ! Map nchmw mw emis. into NEmMwG of Mw frq:
       do j=1,numCompChMW
          xktG(IEmMwG+j-1,j) = xkEmMw(j)
       end do
       rerrAtm(1:nchanmw)=rerr(1:nchanmw)
    end if

    if (iatmnoise==1) then
       rerrAtm(nchanmw+1:nchan)=sqrt(rerr(nchanmw+1:nchan)**2+ &
            AtmNoise(nchanmw+1:nchan)**2)
    else
       rerrAtm(nchanmw+1:nchan)=rerr(nchanmw+1:nchan)
    end if

    nch=0
    alpha1=2.0

    do j = 1,nchan

       if(kchan(j))then

          nch=nch+1

          ! calculated Error covariance matix Sy:
          print *,'chanDataCompressedLoc ',chanDataCompressedLoc
          if (chanDataCompressedLoc) then
            delY(nch) = Ym(j)-Y(nch)
            xktTmp(1:NparG,nch) = xktG(1:NparG,nch)
            print *,ym(j),y(nch)
          else
            delY(nch) = Ym(j)-Y(j)
            xktTmp(1:NparG,nch) = xktG(1:NparG,j)
          end if

          rerr2=rerrAtm(j)**2

          if (j<=nchanmw) then

             Sy(nch)=max(rerr2+dradChan(j)**2,(delY(nch)/alpha1)**2)

          else

             Airs_Err(j)=rerr2-rerrAtm(j)**2+rerr(j)**2
             Airs_Err(j)=sqrt(Airs_Err(j))

             Sy(nch)=max(rerr2,(delY(nch)/alpha1)**2)
          end if

       endif
    end do

    CALL mapK_geo2retr(xktTmp,xkt,vmtx,IG,NG,xG,vCoordTyp,pref,Npar,NparG,nch, &
        nmol,molTran,logSize)

    return
  end subroutine set_IRMW_Invert

END MODULE IRMapInvert
