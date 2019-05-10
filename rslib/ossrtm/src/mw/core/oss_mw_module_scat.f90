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
MODULE oss_mw_module_scat
  USE OSSMWmoduleSubs, only: setpath, ImolG, cosmbk, my_ipol,  my_cFreq, &
          pRef, coef, coef1, coef2, nMols, nch, nChan,nfsmp, NParG, nMol,&
          mxLev, mxLay, mxMols, mxChan, mxParg, MxHmol, &
          ICldLiqG,ITempG,ITskinG,IPsfcG,ICldLiqG, &
          sphericalGeometryFlag, zIntegrationFlag, linInTauFlag, &
          kfix, dkfix, kh2o, dkh2o, kvar, dkvar, ichMap, ImolS, vfreq, &
          rEarth, myRadOrTb, numRefLev, planck, GHzToHz, fpathZ_mw, fpathP_mw, &
          settabindx, OSStran, loadOSSselMW, loadOD, CosmicBackg, oss_destroy_base

  USE CloudDataStruct, only : CloudScene, CloudProfileType
  USE CloudModule, ONLY: initCloud, checkCloudModelType, destroyCloud, &
                setCloudOptSlab, setCloudOptLevel,setCloudOptLayer, &
                paramIndex, typeSlab, typeLevel, typeLayer, typeCloudDefault,&
                getCloudModelType

  USE oss_addbl, only : scattRT
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

  !------------------------------------------------------------------------
  ! Geophysical Vector Variables/Pointers
  !------------------------------------------------------------------------
  REAL,DIMENSION(:),POINTER :: spectPtsGrd
  INTEGER                   :: nSpcPts
  REAL                      :: pobs

  !------------------------------------------------------------------------
  ! Instrument parameters
  !------------------------------------------------------------------------
  REAL,    DIMENSION(:,:), ALLOCATABLE :: coef3
  REAL,    DIMENSION(:),   ALLOCATABLE :: vWvn

  !------------------------------------------------------------------------
  ! Variables saved for interpolation of cloud properties
  !------------------------------------------------------------------------
  INTEGER                               :: ip2
  REAL*8                                :: v2
  REAL,    DIMENSION(mxLay)             :: atabs,btabs,atscat,btscat,aasym,basym
  REAL,    DIMENSION(:,:), ALLOCATABLE  :: tabs2d,tscat2d,asym2d
  REAL,    DIMENSION(:,:,:,:), ALLOCATABLE ,save :: tabsJac,tscatJac,asymJac
  ! convert layer to level
  REAL,    DIMENSION(:,:,:),   ALLOCATABLE ,save :: jacUp,jacLw

  REAL,    DIMENSION(:,:), ALLOCATABLE ,save     :: tabsDat,tscatDat,asymDat
  CHARACTER(len=*), PARAMETER           :: modName='oss_mw_module_scat::'


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
  SUBROUTINE ossinit_mw(selFile,lutFile,ctab_file,nVarMol,varMolID,NparG_in,&
       cldTyp,cldScene,nHyd,nchan_out,chanFreq,nRef,pRef_out,ipol_out,iRadOrTb,&
       sphericalGeometry, zIntegration, linInTau)
    USE OSSPhysicalConstant, only : LtSp_cgs
    !---Input variables
    CHARACTER(LEN=*),                INTENT(IN) :: ctab_file
    CHARACTER(len=*),                INTENT(IN) :: selFile,lutFile
    INTEGER,                         INTENT(IN) :: nVarMol
    INTEGER,           DIMENSION(:), INTENT(IN) :: varMolID
    INTEGER,                         INTENT(IN) :: NparG_in
    type(CloudScene),  dimension(:), intent(in) :: cldScene
    INTEGER,                         INTENT(IN) :: nHyd
    INTEGER,                         INTENT(IN) :: cldTyp
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
    INTEGER                              :: ndim, nProp
    REAL                                 :: alpha,beta
    character(len=*), parameter          :: procName=modName//'ossinit_mw'

    !------------------------------------------------------------------------
    ! Set global variables from file variables
    !------------------------------------------------------------------------
    !---Read absorption coeffs LUTs and OSS coeffs File
    CALL loadOSSselMW(selFile)
    !---Sanity checks
    IF (nchan.GT.mxchan) THEN
       PRINT *,'err['//procName//']:  Nchan > Mxchan :', &
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
    ndim=SIZE(coef,1)
    if (.not. allocated(coef1))         allocate(coef1(ndim,nfsmp))
    if (.not. allocated(coef2))         allocate(coef2(ndim,nfsmp))
    if (.not. allocated(coef3))         allocate(coef3(ndim,nfsmp))
    if (.not. allocated(vWvn))          allocate(vWvn(nfsmp))
    DO ISmp=1,NFSmp
       CALL CosmicBackg(vfreq(ismp),cosmbk(ismp),alpha,beta)
       IF (myRadOrTb == 0) THEN      !TB
          coef1(1:nch(ismp),ismp)=coef(1:nch(ismp),ismp)
          coef2(1:nch(ismp),ismp)=0.
       ELSE                          !Radiance
          coef1(1:nch(ismp),ismp)=coef(1:nch(ismp),ismp)/alpha
          coef2(1:nch(ismp),ismp)=coef(1:nch(ismp),ismp)*(-beta/alpha)
       ENDIF
       ! Factors to convert radiance from wn units to freq units
       coef3(1:nch(ismp),ismp)=coef(1:nch(ismp),ismp)*GHzToHz/ltSp_cgs
       vWvn(ismp) = (vfreq(ismp)*GHzToHz)/ ltSp_cgs
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

    CALL initCloud(ctab_file,cldScene,spectPtsGrd,cldTyp,nHyd)

    nSpcPts=size(spectPtsGrd)

    allocate(tabs2d(mxLay,nSpcPts), tscat2d(mxLay,nSpcPts),&
         asym2d(mxLay,nSpcPts))

    allocate(tabsDat(mxLay,nSpcPts), &
               tscatDat(mxLay,nSpcPts), &
                asymDat(mxLay,nSpcPts))

    ndim=SIZE(cldScene)
    if (cldTyp == typeSlab) then
        allocate(tabsJac (mxLay,nSpcPts,0:4, ndim), &
                 tscatJac(mxLay,nSpcPts,0:4, ndim), &
                 asymJac (mxLay,nSpcPts,0:4, ndim))

    else if(cldTyp == typeLevel) then
        nProp=cldScene(1)%nProperty
        allocate(tabsJac (mxLay,nSpcPts,0:nProp, ndim), &
                 tscatJac(mxLay,nSpcPts,0:nProp, ndim), &
                 asymJac (mxLay,nSpcPts,0:nProp, ndim),&
                 jacUp(mxLay,0:nProp, ndim),jacLw(mxLay,0:nProp, ndim))

    else if(cldTyp == typeLayer) then
        nProp=cldScene(1)%nProperty
        allocate(tabsJac (mxLay,nSpcPts,0:nProp, ndim), &
                 tscatJac(mxLay,nSpcPts,0:nProp, ndim), &
                 asymJac (mxLay,nSpcPts,0:nProp, ndim))
    else
        PRINT *,'err['//procName//']: Unknown cloud model: ', cldTyp
        call exit(1)
    end if


    if (iRadOrTB == 0) then
       iRadOrTB=1
       print *,'warning:[oss_mw_module::oss_init_mw] '// &
            'OSS was trained in TB, but scattering is always in radiance'
       print *,'iRadOrTB is being reset to ', iRadOrTB
    end if

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


    ! linInTauFlag is passed to scatRT
    ! To compare with the previous version linInTauFlag=.TRUE.
    if (present(linInTau)) THEN
       linInTauFlag = linInTau
    else
       linInTauFlag = .TRUE.
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
  ! cldScene define the cloud property as a derived type (cloudScene)
  ! IcldLiq      INTEGER  Starting index for liquid cloud in state vector
  ! IcldIce      INTEGER  Starting index for ice cloud in state vector
  !----------------------------------------------------------------------------
   SUBROUTINE ossdrv_mw_vect(xG,surfEm,obsAngle,obsLevel,y,xkt,xkEmMw,kchan,pland,&
      lat,tempIndex,tSkinIndex,pSurfIndex,varMolIndex, zSurf,zProf,pUser,&
      cldScene, xktCloud,iCldLiq, iCldIce,lambertian)

     use CloudParameters, only: mxHydrometeor
     !---Input variables
     REAL,    DIMENSION(:),   INTENT(IN)       :: xG,surfEm,obsAngle
     LOGICAL, DIMENSION(:),   INTENT(IN)       :: kchan
     REAL,                    INTENT(IN)       :: pland
     REAL,                    INTENT(IN)       :: lat
     INTEGER,                 INTENT(IN)       :: obsLevel
     INTEGER                      , INTENT(IN) :: tempIndex,tSkinIndex,pSurfIndex
     integer,DIMENSION(:)      , intent(in) :: varMolIndex
     real              ,  optional, intent(in) :: zSurf
     real, dimension(:),  optional, intent(in) :: pUser
     real, dimension(:),  optional, intent(in) :: zProf
     logical,             optional, intent(IN) :: lambertian

     integer, intent(in),   optional           :: iCldLiq   ! location water cloud in xG
     integer, intent(in),   optional           :: iCldIce   ! location ice cloud in xG

     !---Output variables
     REAL,    DIMENSION(:),   INTENT(INOUT) :: y,xkEmMw
     REAL,    DIMENSION(:,:), INTENT(INOUT) :: xkt
     TYPE(CloudScene), DIMENSION(:), INTENT(INOUT), optional :: cldScene
     REAL,   DIMENSION(:,:,0:,:),  INTENT(INOUT), OPTIONAL :: xktCloud

     !---Local variables
     logical                          :: lambertianLoc
     LOGICAL                          :: lookup
     REAL                             :: rad,bc,em,solRefl,xkem,xkRefl
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

     REAL,    DIMENSION(mxLay)        :: tabs,tscat,asym
     REAL,    DIMENSION(mxLay)        :: draddtemp
     REAL,    DIMENSION(mxLay,0:5,mxHydrometeor)  :: jacAbs, jacScat, jacAsym

     REAL,    DIMENSION(mxLay)        :: flxnetMono
     REAL                             :: fbeam
     INTEGER, PARAMETER               :: nStr=4,iflux=0
     REAL,    PARAMETER               :: umu0=1.,delphi=0.
     LOGICAL                          :: sBeam=.false.
     logical                          :: isLayerMode
     integer                          :: k, ks
     integer                          :: IXoff
     REAL                             :: tskin_a,emis_a
     REAL, DIMENSION(MxLay)           :: temp_a
     REAL, DIMENSION(MxLay)           :: tautot_a,tabs_a,tscat_a,asym_a, drdw
     real                             :: cosmicBackground
     integer                          :: jj, jjHyd, jjJac
     logical                          :: referenceGrid
     real,    DIMENSION(Mxlev)        :: pLoc
     real,    dimension(MxLay)        :: alt
     integer                          :: cloudMode ! 0-'slab',1-'level',2-'layer'
     real, allocatable                :: xktCloudTemp(:,:,:)
     integer, dimension(2)            :: cldIndex
     type(CloudProfileType)           :: cloudProfile
     real                             :: dummy(MxLay)
     real                             :: ratio, amount
     integer                          :: idx
     integer                          :: iTempLoc
     logical                          :: interpSfc
     integer                          :: layerNum
     integer                          :: nHyd
     integer                          :: nProperty
     CHARACTER(len=*), PARAMETER      :: procName=modName//'ossdrv_mw_vect'
     CHARACTER(len=*), PARAMETER      :: errHeader='Err['//procName//']: '

     !---Size conformity checks (safe interface with calling program/subroutine)
     IF (SIZE(xG)       < NparG) THEN
        PRINT *,errHeader,' Insufficient Size for xG'
        call exit(1)
     END IF
     IF (SIZE(surfEm)   < nchan) Then
        PRINT *,errHeader,' Insufficient Size for surfEm'
        call exit(1)
     END IF
     IF (SIZE(kchan)  < nchan) THEN
        PRINT *,errHeader,' Insufficient Size for kchan'
        call exit(1)
     END IF
     IF (SIZE(Y)      < nchan) THEN
        PRINT *,errHeader,' Insufficient Size for Y'
        call exit(1)
     END IF
     IF (SIZE(xkEmMw) < nchan) THEN
        PRINT *,errHeader,' Insufficient Size for xkEmMw'
        call exit(1)
     END IF
     IF (size(xkt,1)   < NparG) THEN
        PRINT *,errHeader,' Insufficient 1st Dim for xkt'
        call exit(1)
     END IF
     IF (size(xkt,2)     < nchan) THEN
        PRINT *,errHeader,' Insufficient 2nd Dim for xkt'
        call exit(1)
     end if

    if (present(lambertian)) then
     lambertianLoc=lambertian
    else
       lambertianLoc=.false.
    end if

    cloudMode= getCloudModelType()
    select case(cloudMode)
      case(typeLevel)
        if (.not. present(cldScene)) then
           print*,errHeader,' For cloud model "level" variable cldScene must be set'
           call exit(1)
        endif
        nHyd = size(cldScene)
        nProperty = cldScene(1)%nProperty

      case(typeLayer)
        if (.not. present(cldScene)) then
           print*,errHeader,' For cloud model "layer" variable cldScene must be set'
           call exit(1)
        endif
        nHyd = size(cldScene)
        nProperty = cldScene(1)%nProperty

      case(typeSlab)
        if (.not. ( (present(iCldIce) ) .and. (present(iCldLiq))) )  then
           print*,errHeader,' For cloud model "slab" variables'&
                 ,' iCldIce and CldLiq must be set'
           call exit(1)
        endif
        nHyd = 2
        nProperty = 4
      case default
          print *,errHeader, ' Wrong cloudMode', cloudMode
          call exit(1)
    end select

    if (zIntegrationFlag) then
      if ( .not. (present(zSurf) .and. present(zProf))) then
         print*, errHeader,'zIntegration requires  zSurf and zProf to be provided'
         call exit(1)
      end if
    end if
    if (sphericalGeometryFlag) then
      if ( .not. (present(zSurf) .and. present(zProf))) then
         print*, errHeader,'sphericalGeometryFlag requires zSurf and zProf to be provided'
         call exit(1)
      end if
    end if
    !---Initialize radiance vector and k-matrix
    y         = 0.
    xkt       = 0.
    xkEmMw    = 0.
    lookup    = .false.

    ITempG=tempIndex
    ITskinG=tSkinIndex
    IPsfcG=pSurfIndex
    ImolG(1:nmol)=varMolIndex(1:nmol)

    if(present(pUser))then
       if (size(pUser)>mxlev) then
          print*,errHeader,'Input pressure grid is too large.'
          call exit(1)
       endif
       !---Modify geophysical pointers based on size of input
       !     state vector:
      CALL setpath(xG(IPsfcG),pLoc,obsAngle(1:nchan),nSurf,obsLevel,lookup,n1,n2,&
                  umu(1:nchan),referenceGrid, interpSfc, pUser)
    else !initializes to zero to avoid NAN.
      !---Compute path variables
      CALL setpath(xG(IPsfcG),pLoc,obsAngle(1:nchan),nSurf,obsLevel,lookup,n1,n2,&
                  umu(1:nchan),referenceGrid, interpSfc)
    endif

    if ( present(zProf) .and. present(zSurf)) then
         alt(1:nSurf-1)=zProf(1:nSurf-1)
         alt(nSurf)=zSurf
    else
         alt(1:nSurf)=0.
    end if
    !---Compute average temperature and molecular amounts for the layers
    if (zIntegrationFlag) then
       CALL fpathZ_mw(xG,n2,tavl,pavl,wdry,qvar,wvar,dtu,dtl,dwqu,dwql,&
             referenceGrid,interpSfc,pLoc, alt)
    else
      CALL fpathP_mw(xG,lat,n2,tavl,pavl,wdry,qvar,wvar,dtu,dtl,dwqu,dwql,&
             referenceGrid,interpSfc,alt,pLoc)
    end if
    !---Compute coefficients for temperature interpolation of ODs
    CALL settabindx(tavl,qvar(1,1:N2),pavl,N2,referenceGrid,indx_tmp_p1,indx_tmp_p2, &
         indxp,dtmp_p1,dtmp_p2,at_p1,at_p2,ap1,ap2,indx_wvp_p1,indx_wvp_p2)

    if (present(xktcloud)) then
       xktCloud  = 0.
       if (present(cldScene)) &
            allocate(xktCloudTemp(N2+1, 0:nProperty, nHyd))
    end if

    !---Compute cloud optical depth and derivatives
    layerNum=nSurf-1
    select case(cloudMode)
      case (typeSlab) !'slab'
        call setCloudOptSlab(xG,tavl,nSurf,pLoc,xG(pSurfIndex),iCldLiq,iCldIce,&
                            tabsDat,tscatDat,asymDat,cloudProfile,&
                            tabsJac,tscatJac,asymJac)
        layerNum=cloudProfile%levNum-1

      case (typeLevel) !'level'
        if (present(xktcloud)) then
          call setCloudOptLevel(xG,cldScene,nSurf, pLoc,xG(pSurfIndex),&
                     varMolIndex(1),tavl,tabsDat,tscatDat,asymDat, &
                     tabsJac,tscatJac,asymJac,jacUp,jacLw)
        else
          call setCloudOptLevel(xG,cldScene,nSurf,pLoc,xG(pSurfIndex),&
                       varMolIndex(1),tavl,tabsDat,tscatDat,asymDat)
        end if

      case (typeLayer) !'layer'
        if (present(xktcloud)) then
          call setCloudOptLayer(cldScene,nSurf, tavl, &
                                    tabsDat,tscatDat,asymDat, &
                                    tabsJac,tscatJac,asymJac)
        else
          call setCloudOptLayer(cldScene,nSurf, tavl, &
                                     tabsDat,tscatDat,asymDat)
        end if
    end select

    ip2=2
    v2=0.
    iTempLoc=paramIndex('Temp')
    !======================================================================
    !     Loop over spectral points
    !======================================================================
    DO nn=1,NFSmp

       !---Determine optical properties for nnth spectral point:
      CALL optIntLayer(vFreq(nn),layerNum,tabs,tscat,asym,     &
                        size(cldScene), jacAbs, jacScat, jacAsym)

       !---Compute molecular optical depth for all atmospheric layers
      CALL OSStran(kfix(1,nn),dkfix(1,nn),kh2o(1,nn),dkh2o(1,nn),&
            kvar(1,nn),dkvar(1,nn),indx_tmp_p1,indx_tmp_p2,indxp, &
            dtmp_p1,dtmp_p2,at_p1,at_p2,ap1,ap2,indx_wvp_p1,indx_wvp_p2,&
            NmolS(nn),ImolS(1,nn),wdry,qvar,wvar,N2,&
            referenceGrid,tautot,dtaudtmp,dtaudmol)

      if (cloudMode == typeSlab) then
        do ich= 1, layerNum
          cloudProfile%tau(ich) = cloudProfile%molAbs(ich)*tautot(cloudProfile%layerIndex(ich))
        end do
      end if
      bc=cosmbk(nn)
      DO ich=1,nch(nn)
          ich0=ichMap(ich,nn)
          IF (kchan(ich0)) THEN
             em=surfEm(ich0)
             solRefl=1.-em
             cosmicBackground = planck(vfreq(nn), bc)

            ! uncomment the line to ensure code compatibilities with previous version
            ! cosmicBackground=0.

            if (present(xktcloud)) then
               CALL scattRT(tautot,xG,xG(iTskinG),tabs,tscat,asym,em,solRefl, &
                    cosmicBackground, fBeam, vWvn(nn), nSurf-1, obsLevel, &
                    sBeam,nStr,umu(ich0),umu0, lookup,delphi,iflux,flxNetMono,&
                    rad,pland,linInTauFlag, lambertianLoc, tskin_a, emis_a, &
                    temp_a,tautot_a,tabs_a,tscat_a,asym_a)
            else
               CALL scattRT(tautot,xG,xG(iTskinG),tabs,tscat,asym,em,solRefl, &
                    cosmicBackground, fBeam, vWvn(nn), nSurf-1, obsLevel, &
                    sBeam,nStr,umu(ich0),umu0, lookup,delphi,iflux,flxNetMono,&
                    rad,pland,linInTauFlag, lambertianLoc)
            end if

            !------ Jacobians for cloudscene
            if (present(xktcloud)) then
               !--- Jacobians for emissivity
               xkem   = emis_a
               xkRefl = -emis_a

               if (cloudMode==0) then
                !++++++++++++++++++++++++++++++++++++++++++
                !'slab'
                !++++++++++++++++++++++++++++++++++++++++++

                xkt_tmp = 0.
                cldIndex =(/iCldLiq, iCldIce/)

                do jjHyd = 1, nHyd
                   if (cloudProfile%cloudPresent(jjHyd)) then
                      do jjJac = 0, 3
                        do jj=cloudProfile%cloudTop(jjHyd),cloudProfile%cloudBot(jjHyd)-1
                          idx=cldIndex(jjHyd) + jjJac
                          xkt_tmp(idx) = xkt_tmp(idx) + tabs_a(jj)*jacAbs(jj,jjJac,jjHyd) &
                                            + tscat_a(jj)*jacScat(jj,jjJac,jjHyd) + asym_a(jj)*jacAsym(jj,jjJac,jjHyd)
                        end do
                      end do
                  end if
                end do

                do jjHyd = 1, nHyd
                  if (cloudProfile%cloudPresent(jjHyd)) then
                    jj=cloudProfile%cloudBot(jjHyd)
                    ratio = cloudProfile%tau(jj-1)/(cloudProfile%presLevel(jj)-cloudProfile%presLevel(jj-1))
                    amount  = ratio*(tautot_a(jj-1)-tautot_a(jj))

                  ! contribution from temperature variation
                  ! cloud bottom
                    ratio = (cloudProfile%tempLevel(jj)-cloudProfile%tempLevel(jj-1)) &
                      /(cloudProfile%presLevel(jj)-cloudProfile%presLevel(jj-1))

                    amount  = amount + ratio*temp_a(jj)

                    idx=cldIndex(jjHyd) + 1
                    xkt_tmp(idx)=xkt_tmp(idx)+amount

                    jj=cloudProfile%cloudTop(jjHyd)
                    ratio = cloudProfile%tau(jj-1)/(cloudProfile%presLevel(jj)-cloudProfile%presLevel(jj-1))
                    amount  = amount +  ratio*(tautot_a(jj-1)-tautot_a(jj))

                   ! contribution from temperature variation
                   ! cloud top
                    ratio = (cloudProfile%tempLevel(jj)-cloudProfile%tempLevel(jj-1)) &
                           /(cloudProfile%presLevel(jj)-cloudProfile%presLevel(jj-1))

                    amount  = amount +ratio*temp_a(jj)

                    idx=cldIndex(jjHyd)
                    xkt_tmp(idx)=xkt_tmp(idx)+amount
                  end if
                end do

                dummy(1:layerNum)=0.
                do jj=1,N2
                  idx=cloudProfile%layerIndex(jj)
                  dummy(idx) =dummy(idx) + tautot_a(jj)*cloudProfile%molAbs(jj)
                end do
                draddtemp(1:N2)=dummy(1:N2)*dtaudtmp(1:N2)

                jjJac=4
                do jjHyd = 1, 2
                  do jj=1,N2
                    idx=cloudProfile%layerIndex(jj)
                    draddtemp(idx) = draddtemp(idx)&
                                   + tabs_a(jj)*jacAbs(jj,jjJac,jjHyd) &
                                   + tscat_a(jj)*jacScat(jj,jjJac,jjHyd)&
                                   + asym_a(jj)*jacAsym(jj,jjJac,jjHyd)
                  end do
                end do
                !
                xkt_tmp(ITempG+1:ITempG+N2-1)=draddtemp(2:N2)*dtu(2:N2)&
                                             +draddtemp(1:N2-1)*dtl(1:N2-1)
                xkt_tmp(ITempG)=draddtemp(1)*dtu(1)
                xkt_tmp(ITempG+N2)=draddtemp(N2)*dtl(N2)

                !------ Jacobians for molecular concentrations
                DO kS=1,NmolS(nn)
                  k=ImolS(kS,nn)
                  IXoff=ImolG(k)-1
                  drdw(1:N2)                = dummy(1:N2)*dtaudmol(kS,1:N2)
                  xkt_tmp(IXoff+2:IXoff+N2) = drdw(2:N2)*dwqu(2:N2,k)+ &
                                                   drdw(1:N2-1)*dwql(1:N2-1,k)
                  xkt_tmp(IXoff+1)          = drdw(1)*dwqu(1,k)
                  xkt_tmp(IXoff+N2+1)       = drdw(N2)*dwql(N2,k)
                END DO

                !   convert temp_a
                dummy=0.
                do jj=1, N2
                  idx=cloudProfile%layerIndex(jj)
                  dummy(idx)=dummy(idx) + temp_a(jj)*cloudProfile%ratio(jj)
                  dummy(idx+1)=dummy(idx+1) + temp_a(jj)*(1.0 -cloudProfile%ratio(jj))
                end do
                dummy(N2+1)= temp_a(N2 + 1)
                xkt_tmp(ITempG:ITempG+N2) = xkt_tmp(ITempG:ITempG+N2) + dummy(1:N2+1)
              else
                !--- Jacobians for temperature profiles
                draddtemp(1:N2)=tautot_a(1:N2)*dtaudtmp(1:N2)

                do jjHyd = 1, size(cldScene)
                  do jj =1, N2
                    drdw(jj) = tabs_a(jj) *jacAbs (jj,iTempLoc,jjHyd) + &
                            tscat_a(jj)*jacScat(jj,iTempLoc,jjHyd) + &
                            asym_a(jj) *jacAsym(jj,iTempLoc,jjHyd)
                  end do
                  draddtemp(1:N2)=draddtemp(1:N2)+drdw(1:N2)
                end do

                xkt_tmp=0.

                xkt_tmp(ITempG+1:ITempG+N2-1)=draddtemp(2:N2)*dtu(2:N2)&
                                             +draddtemp(1:N2-1)*dtl(1:N2-1)

                xkt_tmp(ITempG)=draddtemp(1)*dtu(1)
                xkt_tmp(ITempG+N2)=draddtemp(N2)*dtl(N2)

                xkt_tmp(ITempG:ITempG+N2) = xkt_tmp(ITempG:ITempG+N2)&
                                           +temp_a(1:N2+1)

                !------ Jacobians for molecular concentrations
                DO kS=1,NmolS(nn)
                  k=ImolS(kS,nn)
                  IXoff=ImolG(k)-1
                  drdw(1:N2)                = tautot_a(1:N2)*dtaudmol(kS,1:N2)
                  xkt_tmp(IXoff+2:IXoff+N2) = drdw(2:N2)*dwqu(2:N2,k)+ &
                                              drdw(1:N2-1)*dwql(1:N2-1,k)

                  xkt_tmp(IXoff+1)          = drdw(1)*dwqu(1,k)
                  xkt_tmp(IXoff+N2+1)       = drdw(N2)*dwql(N2,k)
                END DO

                if (cloudMode==typeLevel) then
                  !+++++++++++++++++++++++++++++++++++++++++++++++++++++
                  !!'level'
                  !++++++++++++++++++++++++++++++++++++
                  do jjHyd = 1, size(cldScene)
                    do jjJac = 0, cldScene(1)%nProperty
                      if (jjJac==iTempLoc) cycle
                        do jj =1, N2
                          drdw(jj) = tabs_a(jj) *jacAbs (jj,jjJac,jjHyd) + &
                                    tscat_a(jj)*jacScat(jj,jjJac,jjHyd) + &
                                     asym_a(jj) *jacAsym(jj,jjJac,jjHyd)
                        end do
                        xktCloudTemp(1,jjJac,jjHyd)    = drdw(1 )*jacUp(1,   jjJac,jjHyd)
                        xktCloudTemp(N2+1,jjJac,jjHyd) = drdw(N2)*jacLw(N2+1,jjJac,jjHyd)

                        do jj =2, N2
                          xktCloudTemp(jj,jjJac,jjHyd) = drdw(jj)  *jacUp(jj,jjJac,jjHyd)+&
                                                          drdw(jj-1)*jacLw(jj,jjJac,jjHyd)
                        end do
                        do jj =1, N2+1
                             xktCloud(jj,ich0,jjJac,jjHyd) =  xktCloud(jj,ich0,jjJac,jjHyd) + &
                                                    xktCloudTemp(jj,jjJac,jjHyd) * coef3(ich,nn)
                        end do
                      end do
                    end do

                elseif (cloudMode==typeLayer) then
               !+++++++++++++++++++++++++++++++++++++++++++++++++++++
               !!'layer'
               !++++++++++++++++++++++++++++++++++++
                  do jjHyd = 1, size(cldScene)
                    do jjJac = 0, cldScene(1)%nProperty
                      if (jjJac==iTempLoc) cycle
                      do jj =1, N2
                        drdw(jj) = tabs_a(jj) *jacAbs (jj,jjJac,jjHyd) + &
                                    tscat_a(jj)*jacScat(jj,jjJac,jjHyd)  + &
                                     asym_a(jj) *jacAsym(jj,jjJac,jjHyd)
                      end do

                      do jj =1, N2
                        xktCloud(jj,ich0,jjJac,jjHyd) =  xktCloud(jj,ich0,jjJac,jjHyd) + &
                                                 drdw(jj) * coef3(ich,nn)
                      end do
                    end do
                  end do
                end if

              end if

              xkt_tmp(iTskinG) = tskin_a

              DO n=1,NParG
                xkt(n,ich0)=xkt(n,ich0)+xkt_tmp(n)*coef3(ich,nn)
              END DO
              xkEmMw(ich0)=xkEmMw(ich0)+xkem*coef3(ich,nn)
            end if
            y(ich0)=y(ich0)+rad*coef3(ich,nn)
          ENDIF
       END DO
    END DO
    return
  END SUBROUTINE ossdrv_mw_vect

  !----------------------------------------------------------------------------
  ! PURPOSE: interface to the 'ossdrv_mw' subroutine for a scalar obsAngle.
  !----------------------------------------------------------------------------
  SUBROUTINE ossdrv_mw_scal(xG,surfEm,obsAngle,obsLevel,y,xkt,xkEmMw,kchan,pland,&
                    lat,tempIndex,tSkinIndex,pSurfIndex,varMolIndex,&
                    zSurf,zProf,pUser,&
                    cldScene, xktCloud,iCldLiq, iCldIce,lambertian)
    !---Input variables
    REAL,    DIMENSION(:),   INTENT(IN)    :: xG,surfEm
    LOGICAL, DIMENSION(:),   INTENT(IN)    :: kchan
    REAL,                    INTENT(IN)    :: lat
    INTEGER,                 INTENT(IN)    :: obsLevel
    integer                   , intent(in) :: tempIndex,tSkinIndex,pSurfIndex
    integer,DIMENSION(:)      , intent(in) :: varMolIndex
    real              ,  optional, intent(in) :: zSurf
    real, dimension(:),  optional, intent(in) :: pUser
    real, dimension(:),  optional, intent(in) :: zProf
    logical,             optional, intent(IN) :: lambertian
    integer, intent(in),   optional           :: iCldLiq   ! location water cloud in xG
    integer, intent(in),   optional           :: iCldIce   ! location ice cloud in xG

    TYPE(CloudScene), DIMENSION(:), INTENT(INOUT) :: cldScene
    REAL,    DIMENSION(:,:,0:,:), INTENT(INOUT), OPTIONAL :: xktCloud

    REAL,                    INTENT(IN)    :: pland,obsAngle
    !---Output variables
    REAL,    DIMENSION(:),   INTENT(INOUT) :: y,xkEmMw
    REAL,    DIMENSION(:,:), INTENT(INOUT) :: xkt
    !---Local variables
    REAL,    DIMENSION(SIZE(y))            :: viewang

    viewang=obsAngle
    CALL ossdrv_mw_vect(xG, surfEm, viewang, obsLevel,y, xkt, xkEmMw, kchan,pland,&
            lat,tempIndex,tSkinIndex,pSurfIndex,varMolIndex,&
            zSurf=zSurf,zProf=zProf,pUser=pUser,&
            cldScene=cldScene,xktCloud=xktcloud,lambertian=lambertian)
    RETURN
  END SUBROUTINE ossdrv_mw_scal


  !----------------------------------------------------------------------------
  ! PURPOSE: Compute radiances and derivatives of radiances
  !          with respect to atmospheric and surface parameters
  !----------------------------------------------------------------------------
  SUBROUTINE ossrad(tautot,dtaudtmp,dtaudmol,tavl,dtu,dtl,dwqu,dwql,&
       NmolS,ImolS,xG,em,bc,N1,N2,umu,lookup,rad,xkt,xkem,taucld, &
       abscld,dtaucdtop,dtaucdthick)
    !---Input variables
    LOGICAL,               INTENT(IN) :: lookup
    INTEGER,               INTENT(IN) :: N1,N2
    INTEGER,               INTENT(IN) :: NmolS
    INTEGER,               INTENT(IN) :: ImolS(NmolS)
    REAL,                  INTENT(IN) :: em,bc,umu
    REAL,    DIMENSION(:), INTENT(IN) :: xG,tautot,dtaudtmp
    REAL,    DIMENSION(:), INTENT(IN) :: taucld,abscld,dtaucdtop,dtaucdthick
    REAL,    DIMENSION(:), INTENT(IN) :: tavl,dtu,dtl
    REAL,  DIMENSION(:,:), INTENT(IN) :: dtaudmol
    REAL,  DIMENSION(:,:), INTENT(IN) :: dwqu,dwql
    !---Output variables
    REAL,                  INTENT(INOUT) :: rad,xkem
    REAL,    DIMENSION(:), INTENT(INOUT) :: xkt

    !---Local variables
    INTEGER                :: l,i,ic2,ic1,lp1,kS,k,IXoff
    REAL                   :: sumtau_dwn,tausfc,draddrsfc,draddemis,ptop,aa
    REAL                   :: draddtskn,sumtau0_up,dtran,rsfc,ts,sec,secdif
    REAL                   :: ttop,b1,b2,btop,txc,tautop,xx1,xx2
    REAL, DIMENSION(MxLev) :: txdn,txup
    REAL, DIMENSION(MxLay) :: drdw,draddtmp,draddtau
    !---Compute air temperature and level index of cloud top
    ptop=xG(ICldLiqG)
    DO i=1,N2+1
       IF(ptop.LT.pref(i))GOTO 10
      end do
10  ic2              = i
    ic1              = ic2-1
    aa               = alog(ptop/pref(ic1))/alog(pref(ic2)/pref(ic1))
    ttop             = aa*(xG(ic2)-xG(ic1))+xG(ic1)
    b1               = 0.5*(xG(ic1)+ttop)
    b2               = 0.5*(xG(ic2)+ttop)
    btop             = ttop

    sec              = 1./umu
    secdif           = sec
    !---Compute transmittance profile along viewing path down to surface
    txdn(N1)         = 1.
    sumtau_dwn       = 0.
    DO l=N1,N2
       sumtau_dwn    = sumtau_dwn+(tautot(l)+taucld(l))*sec
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
       txup(l)       = EXP(-(sumtau_dwn+sumtau0_up*secdif))
      END DO
    rad              = txup(1)*bc
    DO l=1,N2
       dtran         = txup(l+1)-txup(l)
       !---Derivative of downwelling emission wrt to temperature and constituents
       draddtmp(l)   = dtran
       draddtau(l)   = (txup(l)*tavl(l)-rad)*secdif
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
       draddtau(l)   = draddtau(l)*rsfc+(txdn(l+1)*tavl(l)-rad)*sec
       draddtmp(l)   = draddtmp(l)*rsfc+dtran+draddtau(l)*dtaudtmp(l)
       lp1=l+1
       IF(lp1==ic2)THEN
          txc        = EXP(-taucld(l)*sec)
          tautop     = aa*(txdn(ic2)-txdn(ic1))+txdn(ic1)
          IF(txc.LT..8)THEN
             xx1     = (txdn(l)-tautop)*b1
             xx2     = (tautop-txdn(lp1))*b2
             rad     = rad+(xx1+xx2)
    else
             rad     = rad+dtran*tavl(l)
    end if
       ELSE
          rad        = rad+dtran*tavl(l)
       END IF
    ENDDO
    !-----------------------------------------------------------------------
    !     Compute level derivatives and map to array XKT:
    !-----------------------------------------------------------------------
    xkt=0.
    !---Air Temperature
    xkt(ITempG+1:ITempG+N2-1) = draddtmp(2:N2)*dtu(2:N2)+draddtmp(1:N2-1)*dtl(1:N2-1)
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
    return
  END SUBROUTINE ossrad
!--
  !---------------------------------------------------------------------------------
  ! PURPOSE: Given a node frequency vn, perform spectral interpolation to determine
  !           cloud optical properties at vn based on stored values at bounding
  !           hinge points v1 and v2.
  !---------------------------------------------------------------------------------
  SUBROUTINE optInt(vn,nSurf,tabs,tscat,asym)
    !---Input variables
    REAL,               INTENT(IN)     :: vn
    INTEGER,            INTENT(IN)     :: nSurf
    !---Output variables
    REAL, DIMENSION(:), INTENT(INOUT)  :: tabs,tscat,asym
    !---Local variables
    INTEGER                :: nlay,ip,ip1
    REAL                   :: v1,dvint

    nlay = nSurf-1
    if(vn >= v2)then

       DO ip=ip2,nSpcPts-1
          IF(spectPtsGrd(ip).GT.vn) EXIT
    END DO
    ip1=ip-1
       ip2=ip
       v1=spectPtsGrd(ip1)
       v2=spectPtsGrd(ip2)
       dvint=v2-v1

       atabs(1:nlay)  = (tabs2d(1:nlay,ip2)-tabs2d(1:nlay,ip1)) / dvint
       btabs(1:nlay)  = ((v2*tabs2d(1:nlay,ip1))-(v1*tabs2d(1:nlay,ip2))) / dvint

       atscat(1:nlay)  = (tscat2d(1:nlay,ip2)-tscat2d(1:nlay,ip1)) / dvint
       btscat(1:nlay)  = ((v2*tscat2d(1:nlay,ip1))-(v1*tscat2d(1:nlay,ip2))) / dvint

       aasym(1:nlay)  = (asym2d(1:nlay,ip2)-asym2d(1:nlay,ip1)) / dvint
       basym(1:nlay)  = ((v2*asym2d(1:nlay,ip1))-(v1*asym2d(1:nlay,ip2))) / dvint
    endif

    tabs(1:nlay)     = vn * atabs(1:nlay)  + btabs(1:nlay)
    tscat(1:nlay)    = vn * atscat(1:nlay) + btscat(1:nlay)
    asym(1:nlay)     = vn * aasym(1:nlay)  + basym(1:nlay)

    RETURN

  END SUBROUTINE optInt

  !---------------------------------------------------------------------------------
  ! PURPOSE: Given a node frequency vn, perform spectral interpolation to determine
  !           cloud optical properties at vn based on stored values at bounding
  !           hinge points v1 and v2 in the layer mode.
  !---------------------------------------------------------------------------------
  !----
  SUBROUTINE optIntLayer(vn,nlay,tabs,tscat,asym,       &
                       nHydro, jacAbs, jacScat, jacAsym)

    !---Input variables
    REAL,               INTENT(IN)     :: vn
    INTEGER,            INTENT(IN)     :: nlay
    INTEGER,            INTENT(IN)     :: nHydro
    !---Output variables
    REAL, DIMENSION(:), INTENT(INOUT)  :: tabs,tscat,asym
    REAL, DIMENSION(:,0:,:), INTENT(INOUT), optional  :: jacAbs, jacScat, jacAsym
!    !---Local variables
    INTEGER                :: ip,ip1
    REAL*8                 :: v1, dvint
    integer                :: jj
    REAL                   :: slope, offset
    integer                :: jjJac, jjHyd
    logical                :: intepolateJacobian

    if (present(jacAbs) .and. &
        present(jacScat) .and. &
        present(jacAsym) ) then
        intepolateJacobian = .TRUE.
    else
        intepolateJacobian = .FALSE.
        end if

    DO ip=ip2, nSpcPts-1
      IF(spectPtsGrd(ip) > vn) EXIT
    END DO
    ip1=ip-1
    ip2=ip
    v1=spectPtsGrd(ip1)
    v2=spectPtsGrd(ip2)
    dvint=v2-v1
    jacAbs = 0.
    jacScat = 0.
    jacAsym = 0.
    do jj = 1, nlay
         slope   = (tabsDat(jj,ip2)-tabsDat(jj,ip1))
         offset  = (v2*tabsDat(jj,ip1)-v1*tabsDat(jj,ip2))
         tabs(jj)=  (vn * slope  + offset)/ dvint

         slope   = (tscatDat(jj,ip2)-tscatDat(jj,ip1))
         offset  = (v2*tscatDat(jj,ip1)-v1*tscatDat(jj,ip2))
         tscat(jj)=  (vn * slope  + offset)/ dvint

         slope   = (asymDat(jj,ip2)-asymDat(jj,ip1))
         offset  = (v2*asymDat(jj,ip1)-v1*asymDat(jj,ip2))
         asym(jj)=  (vn * slope  + offset)/ dvint
        if (intepolateJacobian) then
            do jjJac=0,size(tabsJac,3)-1
             do jjHyd=1,nHydro
                slope = (tabsJac(jj,ip2,jjJac,jjHyd)-tabsJac(jj,ip1,jjJac,jjHyd)) / dvint
                offset = (v2*tabsJac(jj,ip1,jjJac,jjHyd)-v1*tabsJac(jj,ip2,jjJac,jjHyd)) / dvint
                jacAbs(jj,jjJac,jjHyd)  = vn * slope  + offset

                slope = (tscatJac(jj,ip2,jjJac,jjHyd)-tscatJac(jj,ip1,jjJac,jjHyd)) / dvint
                offset = (v2*tscatJac(jj,ip1,jjJac,jjHyd)-v1*tscatJac(jj,ip2,jjJac,jjHyd)) / dvint
                jacScat(jj,jjJac,jjHyd)  = vn * slope  + offset

                slope = (asymJac(jj,ip2,jjJac,jjHyd)-asymJac(jj,ip1,jjJac,jjHyd)) / dvint
                offset = (v2*asymJac(jj,ip1,jjJac,jjHyd)-v1*asymJac(jj,ip2,jjJac,jjHyd)) / dvint
                jacAsym(jj,jjJac,jjHyd)  = vn * slope  + offset
          end do
       end do
        end if
    end do
    RETURN
  END SUBROUTINE optIntLayer

!--
  subroutine oss_destroy()
  call oss_destroy_base()
    if (allocated(coef3))         deallocate(coef3)
    if (allocated(vWvn))          deallocate(vWvn)
    if (associated(spectPtsGrd))  deallocate(spectPtsGrd)

    if (allocated(tabs2d))        deallocate(tabs2d)
    if (allocated(tscat2d))       deallocate(tscat2d)
    if (allocated(asym2d))        deallocate(asym2d)

    if (allocated(tabsDat))       deallocate(tabsDat)
    if (allocated(tscatDat))      deallocate(tscatDat)
    if (allocated(asymDat))       deallocate(asymDat)
    if (allocated(tabsJac))       deallocate(tabsJac)
    if (allocated(tscatJac))      deallocate(tscatJac)
    if (allocated(asymJac))       deallocate(asymJac)
    call destroyCloud()
  end subroutine oss_destroy
END MODULE oss_mw_module_scat


