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
!   This software and data are covered by U.S. Patent No. 6,584,405
!   and are delivered with limited rights by AER, Inc.
!   See the file README-DATARIGHTS.txt included with this release
!   for additional details.
!
!*************************************************************</f90File>

!--------------------------------------------------------------
!
!  MODULE OSS_MW_MODULE: contains items needed for the
!                OSS Forward Model (MW). AER Inc. 2004.
!  How to use this module: Simply add the statement
!  "USE OSS_MW_MODULE" in the calling program/subroutine.
!  The "oss_init" call should occur only once, in the beginning
!  of the calling program. The "ossdrv_mw" call should occur
!  every time the simulation of the radiances is needed (usually
!  in the loop over the profiles).
!
!--------------------------------------------------------------
MODULE oss_mw_module
  USE OSSMWmoduleSubs, only: setpath, ImolG, cosmbk, my_ipol,  my_cFreq, &
          pRef, coef, coef1, coef2, nMols, nch, nChan,nfsmp, NParG, nMol,&
          mxLev, mxLay, mxMols, mxChan, mxParg, MxHmol, &
          ICldLiqG,ITempG,ITskinG,IPsfcG,ICldLiqG, &
          sphericalGeometryFlag, zIntegrationFlag, linInTauFlag, &
          kfix, dkfix, kh2o, dkh2o, kvar, dkvar, ichMap, ImolS, vfreq, &
          rEarth, myRadOrTb, numRefLev, fpathZ_mw, fpathP_mw, settabindx, &
          osstc_mw, OSStran, loadOSSselMW, loadOD, CosmicBackg, oss_destroy_base
  IMPLICIT NONE
  PRIVATE
  !------------------------------------------------------------------------
  ! Public items made available to calling program(s)
  !------------------------------------------------------------------------
  PUBLIC :: ossinit_mw,ossdrv_mw,oss_destroy

  !------------------------------------------------------------------------
  ! Overloaded subroutines
  !------------------------------------------------------------------------
  INTERFACE ossdrv_mw
     MODULE PROCEDURE ossdrv_mw_scal
     MODULE PROCEDURE ossdrv_mw_vect
  END INTERFACE

! Variable    Definition
! ---------   ------------------------------------------------------------------
! wdry        dry gas amount in layer (g/cm2)
! wvar        variable gas amount in layer (g/cm2)
! qvar        variable gas mass layer mixing ratio, except for H2O = specific humidity
! kvar        abs coef (cm2/g) intercept for variable gases (dry only) post-cull
! dkvar       abs coef (cm2/g) slope for variable gases (dry only) post-cull
! kvar_tmp    abs coef (cm2/g) intercept for variable gases (dry only) pre-cull
! dkvar_tmp   abs coef (cm2/g) slope for variable gases (dry only) pre-cull
! kfix        abs coef (cm2/g) intercept for fixed gases; then per unit dry gas
! dkfix       abs coef (cm2/g) slope for fixed gases; then per unit dry gas
! varDflt     default mixing ratio of variable gases (unitless)

CONTAINS

  !----------------------------------------------------------------------------
  !    PURPOSE: Performs the 1st interface  between the module and the
  !    calling program.
  !    It reads in the unit numbers, the files names, the geophysical vector
  !    pointers and returns back information that is useful to the user,
  !    pertaining to the frequency, the pressure grid, etc.
  !
  !----------------------------------------------------------------------------
  SUBROUTINE ossinit_mw(selFile,lutFile,nVarMol,varMolID,NparG_in,&
       nchan_out,chanFreq,nRef,pRef_out,ipol_out,iRadOrTb,&
       sphericalGeometry, zIntegration, linInTau)
    !---Input variables
    CHARACTER(LEN=*),      INTENT(IN)    :: selFile,lutFile
    INTEGER              , INTENT(IN)    :: nVarMol
    INTEGER, DIMENSION(:), INTENT(IN)    :: varMolID
    INTEGER,               INTENT(IN)    :: NparG_in

    ! optional input
    logical,                optional,intent(in) :: sphericalGeometry
    logical,                optional,intent(in) :: zIntegration
    logical,                optional,intent(in) :: linInTau

    !---Output variables
    INTEGER,                        INTENT(OUT) :: nchan_out,nRef
    REAL,            DIMENSION(:),INTENT(INOUT) :: chanFreq
    REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT) :: pRef_out
    INTEGER,         DIMENSION(:),INTENT(INOUT) :: ipol_out
    !---Input-Output variables
    INTEGER, INTENT(INOUT) :: iRadOrTb ! If input is 0 (TB) or 1 (rad), input
                                       ! iRadOrTb dictates units of computation
                                       ! and ossdrv_mw outputs y and Jacobians.
                                       ! Otherwise, units are same as in OSS
                                       ! training. On output, iRadOrTb reports
                                       ! which units were used.
    !---Local variables
    INTEGER                              :: ismp
    INTEGER                              :: ndim1
    REAL                                 :: alpha,beta

    !------------------------------------------------------------------------
    ! Set global variables from file variables
    !------------------------------------------------------------------------
    !---Read absorption coeffs LUTs and OSS coeffs File
    CALL loadOSSselMW(selFile)
    !---Sanity checks
    IF (nchan.GT.mxchan) THEN
       PRINT *,'err[oss_mw_module::oss_init_mw]:  Nchan > Mxchan :', &
          Nchan,Mxchan
       call exit(1)
    END IF

    CALL loadOD(lutFile,nVarMol,varMolID)

    !------------------------------------------------------------------------
    ! Set global variables from input arguments
    !------------------------------------------------------------------------
    if (.not. allocated(ImolG)) allocate(ImolG(nmol))
    NparG              = NparG_in

    !------------------------------------------------------------------------
    ! Copy to output arguments
    !------------------------------------------------------------------------
    !---Size conformity checks (safe interface with calling program/subroutine)
    IF (SIZE(chanFreq) < nchan) THEN
       PRINT *,'err[oss_mw_module::oss_init_mw]:  Insufficient Size for chanFreq'
       call exit(1)
    END IF
    if (.not. allocated(pRef_out)) allocate(pRef_out(numRefLev))

    IF (SIZE(ipol_out)  < nchan) THEN
       PRINT *,'err[oss_mw_module::oss_init_mw]:  Insufficient Size for IPOL_OUT'
       call exit(1)
    END IF
    nchan_out          = nchan
    nRef           = numRefLev
    chanFreq(1:nchan) = my_cFreq(1:nchan)
    pRef_out(1:numRefLev)   = Pref(1:numRefLev)
    ipol_out(1:nchan)  = my_ipol(1:nchan)
    if (iRadOrTb==0 .or. iRadOrTb==1) myRadOrTb=iRadOrTb
    iRadOrTb           = myRadOrTb

    !---Compute Cosmic Background
    if (.not. allocated(cosmbk))        allocate(cosmbk(nfsmp))
    ndim1=SIZE(coef,1)
    if (.not. allocated(coef1))         allocate(coef1(ndim1,nfsmp))
    if (.not. allocated(coef2))         allocate(coef2(ndim1,nfsmp))
    DO ISmp=1,NFSmp
       CALL CosmicBackg(vfreq(ismp),cosmbk(ismp),alpha,beta)
       IF (myRadOrTb == 0) THEN      !TB
          coef1(1:nch(ismp),ismp)=coef(1:nch(ismp),ismp)
          coef2(1:nch(ismp),ismp)=0.
       ELSE                          !Radiance
          coef1(1:nch(ismp),ismp)=coef(1:nch(ismp),ismp)/alpha
          coef2(1:nch(ismp),ismp)=coef(1:nch(ismp),ismp)*(-beta/alpha)
       ENDIF
    ENDDO

    !---Check adequacy of dimensions for ossdrv_mw
    IF (nVarMol > mxHmol) THEN
       PRINT *,'err[oss_mw_module::oss_init_mw]:  mxHmol not big enough'
       call exit(1)
    END IF
    IF (any(NmolS(:) > mxmolS)) THEN
       PRINT *,'err[oss_mw_module::oss_init_mw]:  mxmolS not big enough'
       call exit(1)
    END IF

    if (present(sphericalGeometry)) THEN
       sphericalGeometryFlag = sphericalGeometry
    else
       sphericalGeometryFlag = .false.
    end if

    if (present(zIntegration)) THEN
       zIntegrationFlag = zIntegration
    else
       zIntegrationFlag = .false.
    end if

    ! linInTauFlag is essential to pass to scatRT
    if (present(linInTau)) THEN
       linInTauFlag = linInTau
    else
       linInTauFlag = .false.
    end if
    RETURN
  END SUBROUTINE ossinit_mw

  !----------------------------------------------------------------------------
  ! PURPOSE: Driver for the microwave radiative transfer model.
  !          Computes both radiances and their jacobians wrt geophysical
  !          parameters.
  ! xG:      Profile vector of geophysical parameters (containing
  !          temperature and constituents profiles, surface
  !          pressure and skin temperature, and cloud parameters
  ! surfEm:    vector of MW surface emissivities.
  ! obsAngle:  Array of earth incidence angles (/channel).
  ! Y:       Vector of radiances or Tbs (depending on iRadOrTb).
  ! xkt:     Array of derivatives of the radiances wrt geophysical parameters.
  ! xkEmMw:  Vector of radiance derivatives wrt surface emissivities.
  ! lambertian:  Flag controling MW downweling emission reflection at the surface
  !              if not set, assumed to be .FALSE. - specular reflection
  !----------------------------------------------------------------------------
  SUBROUTINE ossdrv_mw_vect(xG,surfEm,obsAngle,obsLevel,y,xkt,xkEmMw,kchan,pland, &
       lat,tempIndex,tSkinIndex,pSurfIndex,varMolIndex,ICldLiqG_in,zSurf,zProf,pUser,&
       lambertian)
    !---Input variables
    REAL,    DIMENSION(:),          INTENT(IN) :: xG,surfEm,obsAngle
    LOGICAL, DIMENSION(:),          INTENT(IN) :: kchan
    REAL,                           INTENT(IN) :: pland
    INTEGER,                        INTENT(IN) :: tempIndex,tSkinIndex,pSurfIndex
    INTEGER, DIMENSION(:),          INTENT(IN) :: varMolIndex
    INTEGER,                        INTENT(IN) :: ICldLiqG_in
    INTEGER,                        INTENT(IN) :: obsLevel
    REAL,                           INTENT(IN) :: lat
    real              ,  optional, intent(in) :: zSurf
    real, dimension(:),  optional, intent(in) :: pUser
    real, dimension(:),  optional, intent(in) :: zProf
    logical,             optional, intent(IN) :: lambertian
    !---Output variables
    REAL,    DIMENSION(:),   INTENT(INOUT) :: y,xkEmMw
    REAL,    DIMENSION(:,:), INTENT(INOUT) :: xkt
    !---Local variables
    LOGICAL                          :: lookup
    REAL                             :: radTemp,bc,em,xkem
    INTEGER                          :: nn,ich,ich0,n,nSurf,n1,n2
    INTEGER, DIMENSION(MxLay)        :: indx_tmp_p1,indx_tmp_p2,indxp
    INTEGER, DIMENSION(MxLay)        :: indx_wvp_p1,indx_wvp_p2
    REAL,    DIMENSION(MxLay)        :: dtmp_p1,dtmp_p2
    REAL,    DIMENSION(MxLay)        :: at_p1,at_p2
    REAL,    DIMENSION(MxLay)        :: wdry,dtu,dtl,pavl,ap1,ap2
    REAL,    DIMENSION(MxHmol,MxLay) :: qvar,wvar
    REAL,    DIMENSION(MxLay,MxHmol) :: dwqu,dwql
    REAL,    DIMENSION(MxLay)        :: tavl,tautot,dtaudtmp
    REAL,    DIMENSION(MxmolS,MxLay) :: dtaudmol
    REAL,    DIMENSION(MxParG)       :: xkt_tmp
    REAL,    DIMENSION(mxchan)       :: umu
    REAL,    DIMENSION(MxLay,mxchan) :: abscld,dtaucdtop,dtaucdthick
    REAL,    DIMENSION(mxchan,Mxlay) :: taucld
    logical                          :: referenceGrid
    real,    DIMENSION(Mxlev)        :: pLoc
    real,    dimension(MxLay)        :: alt
    REAL,    DIMENSION(MxLay, mxchan):: umuLay
    real                             :: rRatioSq
    logical                          :: lambertianLoc
    logical                          :: interpSfc
    NparG = size(xG,1)
    !---Size conformity checks (safe interface with calling program/subroutine)
    IF (SIZE(xG)     < NparG) THEN
       PRINT *,'err[oss_mw_module::ossdrv_mw_vect]:  Insufficient Size for xG'
       STOP
    END IF
    IF (SIZE(surfEm)   < nchan) Then
       PRINT *,'err[oss_mw_module::ossdrv_mw_vect]:  Insufficient Size for surfEm'
       STOP
    END IF
    IF (SIZE(kchan)  < nchan) THEN
       PRINT *,'err[oss_mw_module::ossdrv_mw_vect]:  Insufficient Size for kchan'
       STOP
    END IF
    IF (SIZE(Y)      < nchan) THEN
       PRINT *,'err[oss_mw_module::ossdrv_mw_vect]:  Insufficient Size for Y'
       STOP
    END IF
    IF (SIZE(xkEmMw) < nchan) THEN
       PRINT *,'err[oss_mw_module::ossdrv_mw_vect]:  Insufficient Size for xkEmMw'
       STOP
    END IF
    IF (size(xkt,1)   < NparG) THEN
       PRINT *,'err[oss_mw_module::ossdrv_mw_vect]:  Insufficient 1st Dim for xkt'
       STOP
    END IF
    IF (size(xkt,2)   < nchan) THEN
       PRINT *,'err[oss_mw_module::ossdrv_mw_vect]:  Insufficient 2nd Dim for xkt'
       STOP
    END IF
    if (zIntegrationFlag) then
      if ( .not. (present(zSurf) .and. present(zProf))) then
         print*, 'Err[oss_mvr_module::ossdrv]: ',&
                 ' zIntegration requires zSurf and zProf to be provided'
         call exit(1)
      end if
    end if
    if (sphericalGeometryFlag) then
      if ( .not. (present(zSurf) .and. present(zProf))) then
         print*, 'Err[oss_mvr_module::ossdrv]: Spherical  geometry requires',&
                 ' zSurf and zProf to be provided'
         call exit(1)
      end if
    end if
    if (present(lambertian)) then
      lambertianLoc=lambertian
    else
      lambertianLoc=.false.
    end if
    !---Initialize radiance vector and k-matrix
   iTempG                 = tempIndex
   iTskinG                = tSkinIndex
   iPsfcG                 = pSurfIndex
   iMolG(1:nmol)          = varMolIndex(1:nmol)
   ICldLiqG               = ICldLiqG_in

    y         = 0.
    xkt       = 0.
    xkEmMw    = 0.
    lookup    = .false.
    if(present(pUser))then
       if (size(pUser).gt.mxlev) then
          print*,'Err[oss_mw_module::ossdrv]: Input pressure grid is too large.'
          call exit(1)
       endif
       !---Modify geophysical pointers based on size of input
       !     state vector:
       CALL setpath(xG(IPsfcG),pLoc,obsAngle(1:nchan),nSurf,obsLevel,lookup,n1,&
                    n2,umu(1:nchan), referenceGrid, interpSfc, pUser)
    else !initializes to zero to avoid NAN.
       !---Compute path variables
       call setpath(xG(IPsfcG),pLoc,obsAngle(1:nchan),nSurf,obsLevel,lookup,n1,&
                    n2,umu(1:nchan),referenceGrid, interpSfc)
    end if

    if ( present(zProf) .and. present(zSurf)) then
         alt(1:nSurf-1)=zProf(1:nSurf-1)
         alt(nSurf)=zSurf
    else
         alt(1:nSurf)=0.
    end if
    if (sphericalGeometryFlag) then
      do nn=1,nSurf-1
        rRatioSq = ((rEarth+alt(nSurf))/(rEarth+0.5*(alt(nn)+alt(nn+1))))**2
        umuLay(nn, 1:nchan) = sqrt(1.0 - rRatioSq*(1.0 - umu(1:nchan)**2))
      end do
    else
      do nn=1,nchan
        umuLay(1:nSurf-1, nn) = umu(nn)
      end do
    end if

    !---Compute average temperature and molecular amounts for the layers
    if (zIntegrationFlag) then
       CALL fpathZ_mw(xG,n2,tavl,pavl,wdry,qvar,wvar,dtu,dtl,dwqu,dwql,&
                  referenceGrid, interpSfc,pLoc, alt)
    else
      CALL fpathP_mw(xG,lat,n2,tavl,pavl,wdry,qvar,wvar,dtu,dtl,dwqu,dwql,&
                  referenceGrid, interpSfc,alt,pLoc)
    end if
    !---Compute coefficients for temperature interpolation of ODs
    CALL settabindx(tavl,qvar(1,1:N2),pavl,N2,referenceGrid,indx_tmp_p1,&
         indx_tmp_p2,indxp,dtmp_p1,dtmp_p2,at_p1,at_p2,ap1,ap2,indx_wvp_p1,&
         indx_wvp_p2)
    !---Compute cloud optical depth and derivatives
    CALL osstc_mw(nSurf,xG(IPsfcG),pLoc,xG,tavl,taucld,abscld,&
         dtaucdtop,dtaucdthick)
    !======================================================================
    !     Loop over spectral points
    !======================================================================
    DO nn=1,NFSmp
       !---Compute molecular optical depth for all atmospheric layers
       CALL OSStran(kfix(1,nn),dkfix(1,nn),kh2o(1,nn),dkh2o(1,nn),&
            kvar(1,nn),dkvar(1,nn),indx_tmp_p1,indx_tmp_p2,indxp, &
            dtmp_p1,dtmp_p2,at_p1,at_p2,ap1,ap2,indx_wvp_p1,indx_wvp_p2,&
            NmolS(nn),ImolS(1,nn),wdry,qvar,wvar,N2,&
            referenceGrid,tautot,dtaudtmp,dtaudmol)
       bc=cosmbk(nn)
       DO ich=1,nch(nn)
          ich0=ichMap(ich,nn)
          IF (kchan(ich0)) THEN
             em=surfEm(ich0)
             !---Perform RT calculations !clear sky model
             CALL OSSrad(tautot,dtaudtmp,dtaudmol,tavl,dtu,dtl,dwqu,dwql,&
                  NmolS(nn),ImolS(1,nn),xG,em,bc,N1,N2,umuLay(:,ich0),lookup,&
                  radTemp,xkt_tmp,xkem,taucld(ich0,n1:n2),abscld(n1:n2,ich0),&
                  dtaucdtop(n1:n2,ich0),dtaucdthick(n1:n2,ich0),lambertianLoc)
             y(ich0)=y(ich0)+radTemp*coef1(ich,nn)+coef2(ich,nn)
            DO n=1,NParG
                xkt(n,ich0)=xkt(n,ich0)+xkt_tmp(n)*coef1(ich,nn)
            END DO
             xkEmMw(ich0)=xkEmMw(ich0)+xkem*coef1(ich,nn)
          ENDIF
       END DO
    END DO
    RETURN
  END SUBROUTINE ossdrv_mw_vect

  !----------------------------------------------------------------------------
  ! PURPOSE: interface to the 'ossdrv_mw' subroutine for a scalar obsAngle.
  !----------------------------------------------------------------------------

  SUBROUTINE ossdrv_mw_scal(xG,surfEm,obsAngle,obsLevel,y,xkt,xkEmMw,kchan,pland, &
       lat,tempIndex,tSkinIndex,pSurfIndex,varMolIndex,ICldLiqG_in,zSurf,zProf,pUser,&
       lambertian)
    !---Input variables
    REAL,    DIMENSION(:),   INTENT(IN)    :: xG,surfEm
    LOGICAL, DIMENSION(:),   INTENT(IN)    :: kchan
    REAL,                    INTENT(IN)    :: pland,obsAngle
    INTEGER,                 INTENT(IN)    :: tempIndex,tSkinIndex,pSurfIndex
    INTEGER, DIMENSION(:), INTENT(IN) :: varMolIndex
    INTEGER, OPTIONAL,              INTENT(IN) :: ICldLiqG_in
    INTEGER,                 INTENT(IN)    :: obsLevel
    real, intent(in)                       :: lat

    real              ,  optional, intent(in) :: zSurf
    real, dimension(:),  optional, intent(in) :: pUser
    real, dimension(:),  optional, intent(in) :: zProf
    logical,             optional, intent(IN) :: lambertian
    !---Output variables
    REAL,    DIMENSION(:),   INTENT(INOUT) :: y,xkEmMw
    REAL,    DIMENSION(:,:), INTENT(INOUT) :: xkt
    !---Local variables
    REAL,    DIMENSION(SIZE(y))            :: viewang

    viewang=obsAngle
    CALL ossdrv_mw_vect(xG,surfEm,viewang,obsLevel,y,xkt,xkEmMw, kchan, pland, &
       lat,tempIndex,tSkinIndex,pSurfIndex,varMolIndex,ICldLiqG_in=ICldLiqG_in,&
       zSurf=zSurf,zProf=zProf,pUser=pUser,lambertian=lambertian)
    RETURN
  END SUBROUTINE ossdrv_mw_scal


  !----------------------------------------------------------------------------
  ! PURPOSE: Compute radiances and derivatives of radiances
  !          with respect to atmospheric and surface parameters
  !----------------------------------------------------------------------------
  SUBROUTINE ossrad(tautot,dtaudtmp,dtaudmol,tavl,dtu,dtl,dwqu,dwql,&
       NmolS,ImolS,xG,em,bc,N1,N2,umu,lookup,rad,xkt,xkem,taucld, &
       abscld,dtaucdtop,dtaucdthick,lambertian)
    !Parameters
    real,    parameter     :: LAMBERTIAN_REFL_SECANT=1.66
    !---Input variables
    LOGICAL,               INTENT(IN) :: lookup
    INTEGER,               INTENT(IN) :: N1,N2
    INTEGER,               INTENT(IN) :: NmolS
    INTEGER,               INTENT(IN) :: ImolS(NmolS)
    REAL,                  INTENT(IN) :: em,bc
    REAL,DIMENSION(MxLev), INTENT(IN) :: umu

    REAL,    DIMENSION(:), INTENT(IN) :: xG,tautot,dtaudtmp
    REAL,    DIMENSION(:), INTENT(IN) :: taucld,abscld,dtaucdtop,dtaucdthick
    REAL,    DIMENSION(:), INTENT(IN) :: tavl,dtu,dtl
    REAL,  DIMENSION(:,:), INTENT(IN) :: dtaudmol
    REAL,  DIMENSION(:,:), INTENT(IN) :: dwqu,dwql
    LOGICAL,               INTENT(IN) :: lambertian
    !---Output variables
    REAL,                  INTENT(INOUT) :: rad,xkem
    REAL,    DIMENSION(:), INTENT(INOUT) :: xkt

    !---Local variables
    INTEGER                :: l,i,ic2,ic1,lp1,kS,k,IXoff
    REAL                   :: sumtau_dwn,tausfc,draddrsfc,draddemis,ptop,aa
    REAL                   :: draddtskn,sumtau0_up,dtran,rsfc,ts
    REAL                   :: ttop,b1,b2,btop,txc,tautop,xx1,xx2
    REAL, DIMENSION(MxLev) :: txdn,txup
    REAL, DIMENSION(MxLay) :: drdw,draddtmp,draddtau
    real,DIMENSION(MxLev)  :: secRefl,sec


    !---Compute air temperature and level index of cloud top
    ptop=xG(ICldLiqG)
    DO i=1,N2+1
       IF(ptop.LT.pref(i))GOTO 10
    END DO
10  ic2              = i
    ic1              = ic2-1
    aa               = alog(ptop/pref(ic1))/alog(pref(ic2)/pref(ic1))
    ttop             = aa*(xG(ic2)-xG(ic1))+xG(ic1)
    b1               = 0.5*(xG(ic1)+ttop)
    b2               = 0.5*(xG(ic2)+ttop)
    btop             = ttop

    sec(1:N2)      = 1.0/umu(1:N2)
    IF (lambertian)  then
       secRefl(1:N2)=LAMBERTIAN_REFL_SECANT
    else
       secRefl(1:N2)=sec(1:N2)
    end if
    !---Compute transmittance profile along viewing path down to surface
    txdn(N1)         = 1.
    sumtau_dwn       = 0.
    DO l=N1,N2
       sumtau_dwn    = sumtau_dwn+(tautot(l)+taucld(l))*sec(l)
       txdn(l+1)     = EXP(-sumtau_dwn)
    END DO
    tausfc           = txdn(N2+1)
    !---Initialize radiance and derivative arrays
    rad              = 0.
    draddrsfc        = 0.
    draddtmp(1:n2)   = 0.
    draddtau(1:N2)   = 0.
    draddemis        = 0.
    draddtskn        = 0.
    !-----------------------------------------------------------------------
    !     Downwelling thermal radiance calculation:
    !-----------------------------------------------------------------------
    txup(N2+1)       = tausfc
    sumtau0_up       = 0.
    DO l=N2,1,-1
       sumtau0_up    = sumtau0_up+tautot(l)+taucld(l)
       txup(l)       = EXP(-(sumtau_dwn+sumtau0_up*secRefl(l)))
    END DO
    rad              = txup(1)*bc
    DO l=1,N2
       dtran         = txup(l+1)-txup(l)
       !---Derivative of downwelling emission wrt to temperature and constituents
       draddtmp(l)   = dtran
       draddtau(l)   = (txup(l)*tavl(l)-rad)*secRefl(l)
       rad           = rad+dtran*tavl(l)
    ENDDO
    !-----------------------------------------------------------------------
    !     Add Surface terms
    !-----------------------------------------------------------------------
    rsfc             = (1.-em)
    ts               = xG(ITskinG)
    !---Derivatives wrt emissivity and sfc skin temperature:
    draddemis        = tausfc*ts-rad
    draddtskn        = em*tausfc
    rad              = rad*rsfc+em*tausfc*ts
    !-----------------------------------------------------------------------
    !     Upwelling thermal radiance calculation
    !-----------------------------------------------------------------------
    DO l=N2,1,-1
       dtran         = txdn(l)-txdn(l+1)
       draddtau(l)   = draddtau(l)*rsfc+(txdn(l+1)*tavl(l)-rad)*sec(l)
       draddtmp(l)   = draddtmp(l)*rsfc+dtran+draddtau(l)*dtaudtmp(l)
       lp1=l+1
       IF(lp1.EQ.ic2)THEN
          txc        = EXP(-taucld(l)*sec(l))
          tautop     = aa*(txdn(ic2)-txdn(ic1))+txdn(ic1)
          IF(txc.LT..8)THEN
             xx1     = (txdn(l)-tautop)*b1
             xx2     = (tautop-txdn(lp1))*b2
             rad     = rad+(xx1+xx2)
          ELSE
             rad     = rad+dtran*tavl(l)
          END IF
       ELSE
          rad        = rad+dtran*tavl(l)
       END IF
    ENDDO
    !-----------------------------------------------------------------------
    !     Compute level derivatives and map to array XKT:
    !-----------------------------------------------------------------------
    xkt=0.
    !---Air Temperature
    xkt(ITempG+1:ITempG+N2-1) = draddtmp(2:N2)*dtu(2:N2)&
                               +draddtmp(1:N2-1)*dtl(1:N2-1)
    xkt(ITempG)       = draddtmp(1)*dtu(1)
    xkt(ITempG+N2)    = draddtmp(N2)*dtl(N2)
    !---Molecular concentrations
    DO kS=1,NmolS
       k=ImolS(kS)
       IXoff=ImolG(k)-1
       drdw(1:N2)            = draddtau(1:N2)*dtaudmol(kS,1:N2)
       xkt(IXoff+2:IXoff+N2) = drdw(2:N2)*dwqu(2:N2,k)+ &
                               drdw(1:N2-1)*dwql(1:N2-1,k)
       xkt(IXoff+1)            = drdw(1)*dwqu(1,k)
       xkt(IXoff+N2+1)         = drdw(N2)*dwql(N2,k)
    END DO
    !---Cloud parameters
    xkt(iCldLiqG)      = SUM(draddtau(1:N2)*dtaucdtop(1:N2))
    xkt(iCldLiqG+1)    = SUM(draddtau(1:N2)*dtaucdthick(1:N2))
    xkt(iCldLiqG+2)    = SUM(draddtau(1:N2)*abscld(1:N2))
    !---Surface terms
    xkt(ITskinG)      = draddtskn
    xkEm              = draddemis
    RETURN
  END SUBROUTINE ossrad

!--
  subroutine oss_destroy()
    call oss_destroy_base()
  end subroutine oss_destroy
END MODULE oss_mw_module


