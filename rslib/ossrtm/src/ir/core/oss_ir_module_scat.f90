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
!  MODULE OSS_IR_MODULE_SCAT: contains items needed for the
!                OSS Forward Model (IR). AER Inc. 2004.
!
!--------------------------------------------------------------
MODULE oss_ir_module_scat
  use params_module, ONLY: mxlev, mxlay, mxhmol, maxcmu, MxmolS, mxParG
  USE CloudDataStruct
  USE asolar_source_function
  USE oss_addbl, ONLY: scattRT
  IMPLICIT NONE
  PRIVATE
  !------------------------------------------------------------------------
  ! Public items made available to calling program(s)
  !------------------------------------------------------------------------
  PUBLIC :: ossinit_ir, ossdrv_ir, OSSdestroyIR

  !------------------------------------------------------------------------
  ! Variables saved for interpolation of cloud properties
  !------------------------------------------------------------------------
  INTEGER                                   :: nSpcPts
  REAL,    DIMENSION(:),       POINTER      :: spectPtsGrd
  REAL,    DIMENSION(:,:),     ALLOCATABLE  :: tabsDat
  REAL,    DIMENSION(:,:),     ALLOCATABLE  :: tscatDat
  REAL,    DIMENSION(:,:),     ALLOCATABLE  :: asymDat
  REAL,    DIMENSION(:,:,:,:), ALLOCATABLE  :: tabsJac
  REAL,    DIMENSION(:,:,:,:), ALLOCATABLE  :: tscatJac
  REAL,    DIMENSION(:,:,:,:), ALLOCATABLE  :: asymJac
  REAL,    DIMENSION(:,:,:),   ALLOCATABLE  :: jacUp  !layer to level
  REAL,    DIMENSION(:,:,:),   ALLOCATABLE  :: jacLw  !layer to level

  ! scattering regression
  INTEGER, dimension(:),       allocatable :: nchP_ir
  INTEGER, dimension(:),       allocatable :: iselP_ir
  REAL,    dimension(:,:),     allocatable :: coefP_ir
  INTEGER, dimension(:,:),     allocatable :: ichPmap_ir
  INTEGER                                  :: np_ir
  logical                                  :: getOSSscatRegrDone=.false.
  logical                                  :: ossScatRegr
  LOGICAL                                  :: scatterOn
  CHARACTER(len=*), PARAMETER              :: modName='oss_ir_module_scat::'

CONTAINS

  !----------------------------------------------------------------------------
  !    PURPOSE: Performs the 1st interface  between the module and the
  !    calling program.
  !    It reads in the unit numbers, the files names, the geophysical vector
  !    pointers and returns back information that is useful to the user,
  !    pertaining to the wavenumber, the pressure grid, etc.
  !
  !----------------------------------------------------------------------------
  SUBROUTINE ossinit_ir(selFile,lutFile,solFile,ctab_file,pc_file,cldTyp,nHyd,&
       cldScene, ossScatRegr_in, enableScatter,&
       defProfFile,nVarMol,varMolID,nChanOut, chanIndex, chanFreq, &
       nRef,pRefOut,sphericalGeometry, zIntegration,linInTau)

! ARGUMENTS:
!
!   INPUTS:
!
!   selFile         CHAR     OSS coefficients file path
!   lutFile         CHAR     OSS optical properties LUT file path
!   solFile         CHAR     solar data file path
!   ctab_file       CHAR     Cloud optical properties LUT file path
!   pc_file         CHAR
!   nVarMol         INTEGER  Number of molecular species
!   varMolID        INTEGER  List of IDs of relevant molecular species
!   ossScatRegr_in  LOGICAL
!   enableScatter   LOGICAL  Flag for including multiple scattering in
!                            calculations
!   cldTyp          INTEGER  Type of cloud model used
!   nHyd            INTEGER  Number of hydrometeor types
!   cldScene        type(CloudScene)   structure to store cloud properties
!
!   INPUTS OPTIONAL
!   sphericalGeometry LOGICAL
!   zIntegration      LOGICAL
!   linInTau          LOGICAL
!
!   INPUTS/OUTPUTS:
!
!   nChanOut        INTEGER  Number of channels
!   chanFreq        REAL     Center wavenumbers for instrument channels
!   nRef            INTEGER  Number of atmospheric levels
!   pRefOut         REAL     Pressure on atmospheric grid levels
!
!   * OPTIONAL


    use oss_ir_module, ossinit_irClear => ossinit_ir
    use CloudModule, only: initCloud, typeSlab, typeLevel, typeLayer, &
                           typeCloudDefault

      !---Input variables
    CHARACTER(len=*),                INTENT(IN)    :: selFile,lutFile
    CHARACTER(len=*),                INTENT(IN)    :: solFile,ctab_file
    CHARACTER(len=*),                INTENT(IN)    :: defProfFile
    LOGICAL,                         INTENT(IN)    :: enableScatter
    INTEGER,                         INTENT(IN)    :: cldTyp
    CHARACTER(LEN=*),                INTENT(IN)    :: pc_file
    INTEGER,                         INTENT(IN)    :: nHyd
    INTEGER,                         INTENT(IN)    :: nVarMol
    INTEGER,           DIMENSION(:), INTENT(IN)    :: varMolID
    LOGICAL,                         INTENT(IN)    :: ossScatRegr_in
    type(CloudScene),  dimension(:), intent(in)    :: cldScene

    !---Output variables
    INTEGER,                         INTENT(INOUT) :: nChanOut,nRef
    REAL,              DIMENSION(:), INTENT(INOUT) :: chanFreq
    integer,allocatable,dimension(:),intent(inout) :: chanIndex
    REAL,allocatable,  DIMENSION(:), INTENT(INOUT) :: pRefOut

    !---OPTIONAL Input variables
    LOGICAL,          OPTIONAL,INTENT(IN) :: linInTau
    logical,          optional,intent(in) :: sphericalGeometry
    logical,          optional,intent(in) :: zIntegration

    !local
    integer                               :: nProp
    character(len=*), parameter           :: procName=modName//'ossinit_ir'

    call ossinit_irClear(selFile,lutFile,solFile, &
       defProfFile,nVarMol,varMolID,nChanOut, chanIndex, chanFreq, &
       nRef,pRefOut,sphericalGeometry, zIntegration,linInTau)

    scatterOn          = enableScatter
    ossScatRegr        = ossScatRegr_in
    CALL getOSSScatRegr(pc_file)

    if (enableScatter) then
      call initCloud(ctab_file, cldScene, spectPtsGrd, cldTyp, nHyd)

      nSpcPts=size(spectPtsGrd)

       allocate(tabsDat(mxLay,nSpcPts), &
               tscatDat(mxLay,nSpcPts), &
                asymDat(mxLay,nSpcPts))

      nProp = cldScene(1)%nProperty
      select case(cldTyp)

        case(typeSlab)
          allocate(tabsJac (mxLay,nSpcPts,0:4, nHyd), &
                   tscatJac(mxLay,nSpcPts,0:4, nHyd), &
                   asymJac (mxLay,nSpcPts,0:4, nHyd))

        case(typeLevel)
          allocate(tabsJac (mxLay,nSpcPts,0:nProp, nHyd), &
                   tscatJac(mxLay,nSpcPts,0:nProp, nHyd), &
                   asymJac (mxLay,nSpcPts,0:nProp, nHyd),&
                   jacUp(mxLay,0:nProp, nHyd),&
                   jacLw(mxLay,0:nProp, nHyd))

        case(typeLayer)
         allocate(tabsJac (mxLay,nSpcPts,0:nProp, nHyd), &
                  tscatJac(mxLay,nSpcPts,0:nProp, nHyd), &
                  asymJac (mxLay,nSpcPts,0:nProp, nHyd))
        case default
         PRINT *,'err['//procName//']: unknown cloud model: ', cldTyp
         call exit(1)
      end select
    end if

    RETURN
  END SUBROUTINE ossinit_ir

  SUBROUTINE getOSSScatRegr(pfile)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   getOSSScatRegr
!
! PURPOSE:
!
!   Describe sub/fnc/main here
!
! SYNTAX:
!
!   CALL getOSSScatRegr(pfile)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   iu     INTEGER  file unit number for predictors
!   pfile  CHAR     predictor data file path
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
  use OSSIRmoduleSubs, only : nfsmp_ir,coef_ir,nch_ir,nchmax,ichmap_ir

    !---Input variables
    CHARACTER(len=*),INTENT(IN) :: pfile

    !---Local variables
    INTEGER                     :: ismp,nchPmax

    if (getOSSscatRegrDone) return
    !---open file
    if (ossScatRegr) then
       call loadOSSScatRegr(pfile)
    else !if predictors are not used,
         !re-load the oss-coefficients' info into variables and arrays
         !for proper use in ossdrv_ir
       np_ir=nfsmp_ir
       nchPmax=nchmax
       allocate (coefP_ir(nchPmax,np_ir),ichPmap_ir(nchPmax,np_ir))
       allocate (iselP_ir(np_ir),nchP_ir(np_ir))
       DO ismp=1,np_ir
          iselP_ir(ismp)=ismp
          nchP_ir(ismp)=nch_ir(ismp)
          coefP_ir(1:nchP_ir(ismp),ismp)=coef_ir(1:nchP_ir(ismp),ismp)
          ichPmap_ir(1:nchP_ir(ismp),ismp)=ichmap_ir(1:nchP_ir(ismp),ismp)
       ENDDO
    endif
    getOSSscatRegrDone = .true.
    RETURN
  END SUBROUTINE getOSSScatRegr

  SUBROUTINE loadOSSScatRegr(pfile)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   loadOSSScatRegr
!
! PURPOSE:
!
!   Describe sub/fnc/main here
!
! SYNTAX:
!
!   CALL loadOSSScatRegr(iu, pfile)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   iu     INTEGER  file unit number for predictors
!   pfile  CHAR     predictor data file path
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
    use OSSIRmoduleSubs, only : findFreeUnit

    CHARACTER(len=*),INTENT(IN) :: pfile
    !---Local variables
    INTEGER                     :: iu
    INTEGER                     :: ismp,nchPmax
    iu = findFreeUnit()
    OPEN(iu,file=pfile,form='unformatted',status='old', action='read')
    READ(iu)np_ir,nchPmax
    allocate (coefP_ir(nchPmax,np_ir),ichPmap_ir(nchPmax,np_ir))
    allocate (iselP_ir(np_ir),nchP_ir(np_ir))
    DO ismp=1,np_ir
       READ(iu)iselP_ir(ismp),nchP_ir(ismp)
       READ(iu)coefP_ir(1:nchP_ir(ismp),ismp),ichPmap_ir(1:nchP_ir(ismp),ismp)
    ENDDO
    CLOSE(iu)
    RETURN
  END SUBROUTINE loadOSSScatRegr

  !----------------------------------------------------------------------------
  ! PURPOSE: Driver for the IR radiative transfer model.
  !          Computes both radiances and their jacobians wrt geophysical
  !          parameters.
  !jacRequest set to true if Jacobian required to be computed
  !
  ! xG:      Profile vector of geophysical parmaters (containing
  !          temperature and constiutuents profiles, surface
  !          pressure and skin temperature, and cloud parameters
  ! surfEmRf:    vector of MW surface emissivities.
  ! Path:    Array of quantities characterizing the atmospehric
  !          path in plane-parallel atmosphere.
  ! Y:       Vector of radiances.
  ! xkt:     Array of derivatives of the radiances wrt geophysical parameters.
  ! xkEmRf:  Array of radiance derivatives wrt surface emissivities.
  !
  !----------------------------------------------------------------------------
  SUBROUTINE ossdrv_ir(xG,surfEmRfGrid,surfEmRf,obsAngle,solZenith,azAngle,&
                    obsLevel,y,xkt,xkEmRf,lat,tempIndex,tSkinIndex,pSurfIndex,&
                    varMolIndex,lambertian,zSurf,zProf,pUser,cldScene, &
                    xktCloud,iCldLiq, iCldIce)

    use OSSIRmoduleSubs, only : kself, kfix_ir, kvar_ir, kh2o_ir, dkh2o_ir,&
                    fpathZ, fpathP, layerAverage, OSStran, sunrad,&
                    myLinInTau, sphericalGeometryFlag,zIntegrationFlag,&
                    vWvn,ImolS, NmolS, iselS, coef_ir, nch_ir, &
                    f1_arr, f2_arr, ichMap_ir, setpath_ir, settabindx_ir,&
                    setIndex_selfCont, vinterp,GetOD,getCountChannel,&
                    getCountUsedNode

    use CloudModule, only: setCloudOptLevel, checkCloudModelType, &
                    setCloudOptSlab, setCloudOptLayer,typeCloudDefault,&
                    printCloudProfile,paramIndex,typeSlab,typeLevel,typeLayer,&
                    getCloudModelType

    use CloudParameters, only: mxHydrometeor, mxPhysicalProperty
    use oss_ir_module, only : OSSrad

    REAL,   DIMENSION(:)         , INTENT(INOUT) :: xG
    REAL,   DIMENSION(:)         , INTENT(IN) :: surfEmRfGrid
    REAL,   DIMENSION(:,:)       , INTENT(IN) :: surfEmRf
    REAL                         , INTENT(IN) :: obsAngle,solZenith,azAngle
    REAL                         , INTENT(IN) :: lat
    REAL,   DIMENSION(:),OPTIONAL, INTENT(IN) :: pUser
    REAL,   DIMENSION(:),OPTIONAL, INTENT(IN) :: zProf
    REAL                ,OPTIONAL, INTENT(IN) :: zSurf
    INTEGER                      , INTENT(IN) :: tempIndex,tSkinIndex,pSurfIndex
    type(CloudScene), dimension(:), intent(inout), optional :: cldScene
    ! location water cloud in xG
    integer, intent(in),   optional               :: iCldLiq
    ! location ice cloud in xG
    integer, intent(in),   optional               :: iCldIce

    INTEGER,DIMENSION(:)         , INTENT(IN) :: varMolIndex
    INTEGER, optional            , INTENT(IN) :: lambertian

    INTEGER                 ,     INTENT(IN)    :: obsLevel
    REAL,   DIMENSION(:)    ,     INTENT(INOUT) :: y
    REAL,   DIMENSION(:,:)  ,     INTENT(INOUT) :: xkt
    REAL,   DIMENSION(:, :,:) ,   INTENT(INOUT) :: xkEmRf
    REAL,   DIMENSION(:,:,0:,:),  INTENT(INOUT), OPTIONAL :: xktCloud

    !---Local variables
    INTEGER, PARAMETER               :: iflux=0,nStr=4 ! must be <=maxcmu
    LOGICAL                          :: lookup,sun,referenceGrid,interpSfc
    LOGICAL                          :: lambertianLoc
    INTEGER                          :: nsf,nparG
    REAL,    DIMENSION(Mxlev)        :: pLoc
    REAL,    DIMENSION(MxLay)        :: tavl
    REAL,    DIMENSION(MxLay)        :: pavl,wfix,dtu,dtl,ap1,ap2
    REAL,    DIMENSION(MxHmol,MxLay) :: q,w
    REAL,    DIMENSION(MxLay,0:MxHmol) :: dwqu,dwql
    real,    DIMENSION(MxLay,MxHmol) :: wvdwqu,wvdwql

    INTEGER, DIMENSION(MxLay)        :: indxt_p1,indxt_p2,indxp,indxt
    real,    DIMENSION(MxLay)        :: at1_p1,at2_p1,at1_p2,at2_p2,at1,at2
    REAL,    DIMENSION(MxLay)        :: adt1_p1,adt2_p1,adt1_p2,adt2_p2
    REAL,    DIMENSION(MxLay)        :: tautot,dtaudtmp
    REAL,    DIMENSION(0:MxmolS,MxLay) :: abso
    REAL,    DIMENSION(MxLay)        :: umuLay,umu0Lay
    REAL,    DIMENSION(mxLay)        :: tabs,tscat,asym
    REAL,    DIMENSION(mxLay,0:4,2)  :: jacAbs, jacScat, jacAsym

    REAL,    DIMENSION(mxLay)        :: flxnetMono
    REAL                             :: xkEmRf_tmp(2)
    REAL,    DIMENSION(mxParG)       :: xkt_tmp
    REAL                             :: umu,umu0,vn,fbeam,xx,rad
    REAL                             :: coefInt,tsfc,plogu,plogl
    INTEGER                          :: n1,n2,n,nn,ich,i,ich0,nSurf
    INTEGER                          :: k,ks,ixOff, jj
    real                             :: coefTmp
    INTEGER                          :: nNodes
    INTEGER                          :: nChSel, jjHyd, jjJac
    real,    DIMENSION(MxLay)        :: alt

    type(CloudProfileType)           :: cloudProfile

    REAL                             :: vEmRf(size(surfEmRf, dim=2))
    real                             :: reflectance, alpha
    real                             :: emissivity
    integer                          :: nEnd, cloudMode
    logical                          :: jacRequest = .true.
    real, parameter                  :: cosmicBackground=0.
    real, parameter                  :: pland=1.
    real, allocatable                :: xktCloudTemp(:,:,:)
    integer                          :: idx, cldIndex(2)
    REAL,    DIMENSION(mxLay)        :: dummy
    REAL,    DIMENSION(mxLay)        :: draddtemp
    REAL, DIMENSION(MxLay)           :: temp_a
    REAL, DIMENSION(MxLay)           :: tautot_a,tabs_a,tscat_a,asym_a, drdw
    real                             :: amount,emis_a, tskin_a,ratio

    integer                          :: iTempLoc
    integer                          :: layerNum
    integer                          :: obsLevelLoc
    INTEGER                          :: idxCloud, idxEmRf
    CHARACTER(len=*), PARAMETER      :: procName=modName//'ossdrv_ir'
    CHARACTER(len=*), PARAMETER      :: errHeader='Err['//procName//']:'

    nChSel=getCountChannel()
    obsLevelLoc=obsLevel
    nNodes = getCountUsedNode()

    nparG  = size(xG)

    nsf = size(surfEmRfGrid)
    if (size(surfEmRf,1) < nsf) then
       print*,errHeader,' Hinge point dimension of surfEmRf is too small', &
              size(surfEmRf,2), nsf
       call exit(1)
    endif
    if (size(xkEmRf,1) < nsf) then
       print*,errHeader,' Hinge point dimension of xkEmRf is too small', &
              size(surfEmRf,2), nsf
       call exit(1)
    endif
    if (size(xkEmRf,2) < nchSel) then
       print*,errHeader,' Channel dimension of xkEmRf is too small', &
              size(xkEmRf,2) , nchSel
       call exit(1)
    endif
    if (size(xkt,1) < nparG) then
       print*,errHeader,' Parameter dimension of xkt is too small', &
              size(xkt,1), size(xG)
       call exit(1)
    endif
    if (size(xkt,2) < nchSel) then
       print*,errHeader,' Channel dimension of xkt is too small', &
              size(xkt,2), nchSel
       call exit(1)
    endif

    if (present(lambertian)) then
      lambertianLoc=lambertian>0
    else
      lambertianLoc=.false.
    end if

    IF (SIZE(varMolIndex)> mxHmol) then
       print*,errHeader,' varMolIndex Vector too large', &
              size(varMolIndex,1), mxHmol
      call exit(1)
    end if

    cloudMode=getCloudModelType()
    select case(cloudMode)
      case(typeLevel)
        if (.not. present(cldScene)  ) then
           print*,errHeader,' For cloud model "level" ', &
                  'variable cldScene must be set'
           call exit(1)
        endif
      case(typeLayer)
        if (.not. present(cldScene)  ) then
           print*,errHeader,' For cloud model "layer" ',&
                  'variable cldScene must be set'
           call exit(1)
        endif
      case(typeSlab)
        if (.not. ( (present(iCldIce) ) .and. (present(iCldLiq))) )  then
           print*,errHeader,' For cloud model "slab" ', &
                  'variables iCldIce and CldLiq must be set'
           call exit(1)
        endif
    end select

    if (zIntegrationFlag) then
      if ( .not. (present(zSurf) .and. present(zProf))) then
         print*, errHeader,'zIntegration requires  zSurf and ',&
                'zProf to be provided'
         call exit(1)
      end if
    end if
    if (sphericalGeometryFlag) then
      if ( .not. (present(zSurf) .and. present(zProf))) then
         print*, errHeader,'sphericalGeometryFlag requires  zSurf and ', &
                'zProf to be provided'
         call exit(1)
      end if
    end if


    !======================================================================
    !     Initialize radiance vector and k-matrix
    !======================================================================
    y(1:NchSel)               = 0.

    DO ich=1,NchSel
        xkt(:,ich)     = 0.
        xkEmRf(1:Nsf,ich,1:2)  = 0.
    END DO

    !======================================================================
    !     Compute path variables
    !======================================================================

    if(present(pUser))then
       if (size(pUser)>mxlev) then
          print*,'Err[',procName,']: Input pressure grid is too large.'
          call exit(1)
       endif
       !---Modify geophysical pointers based on size of input
       !     state vector:
      CALL setpath_ir(xG(pSurfIndex),pLoc,obsLevelLoc,obsAngle,solZenith,nSurf,sun, &
                                     umu,umu0,lookup,referenceGrid,interpSfc,pUser)
    else !initializes to zero to avoid NAN.
      !---Compute path variables
      CALL setpath_ir(xG(pSurfIndex),pLoc,obsLevelLoc,obsAngle,solZenith,nSurf,sun, &
                                     umu,umu0,lookup,referenceGrid,interpSfc)
    endif

    n2=nSurf-1
    n1=obsLevelLoc
    if ( present(zProf) .and. present(zSurf)) then
         alt(1:N2)=zProf(1:N2)
         alt(N2+1)=zSurf
    else
         alt(1:N2+1)=0.
    end if

!---Compute average temperature, pressure and integrated molecular amounts for the layers
    call layerAverage(xG,N2,referenceGrid,interpSfc,&
             tempIndex,pSurfIndex,pLoc,pavl,&
             plogu,plogl, umu,umuLay,umu0,umu0Lay,alt, nEnd, alpha)

    if (zIntegrationFlag) then
      call fpathZ(xG,N2,nEnd,interpSfc,tempIndex, pSurfIndex,varMolIndex,pLoc,&
             tSfc,tavl,dtu,dtl,alpha,wfix,q,w,dwqu,dwql,wvdwqu,wvdwql,alt)
    else
      call fpathP(xG,lat,N2,nEnd,interpSfc,pSurfIndex,tempIndex,varMolIndex, &
             pLoc,tSfc,tavl,dtu,dtl,alpha, wfix,q,w,dwqu,dwql,wvdwqu,wvdwql,alt)
    end if

    if (present(xktcloud)) then
       xktCloud  = 0.
       if (present(cldScene)) allocate(xktCloudTemp(N2+1, 0:cldScene(1)%nProperty, size(cldScene)))
    end if

    tabs(:)=0.
    tscat(:)=0.
    asym(:)=0.
    !if (lookup .and. lambertianLoc /= 0) stop 'Incompatibility'
    !---Compute coefficients for temperature interpolation of ODs
    CALL settabindx_ir(tavl,pavl,N2,referenceGrid, &
         indxt_p1,indxt_p2,indxp,ap1,ap2,&
         at1_p1,at2_p1,at1_p2,at2_p2,adt1_p1,adt2_p1,adt1_p2,adt2_p2)

    CALL setIndex_selfCont(tavl,N2,indxt,at1,at2)

    !---Get profiles of optical properties for input to scattering module:
    !---Compute cloud optical depth and derivatives
    if (scatterOn) then
      layerNum=N2
      select case(cloudMode)
        case (typeSlab) !'slab'
          call setCloudOptSlab(xG,tavl,nSurf,pLoc, &
                             xG(pSurfIndex),iCldLiq,iCldIce,tabsDat,tscatDat,asymDat, &
                              cloudProfile, tabsJac,tscatJac,asymJac)

          layerNum=cloudProfile%levNum-1

        case (typeLevel) !'profile'
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
              tabsDat,tscatDat,asymDat,tabsJac,tscatJac,asymJac)
          else
            call setCloudOptLayer(cldScene,nSurf, tavl, &
               tabsDat,tscatDat,asymDat)
          end if
      end select
    end if
    !======================================================================
    !     Loop over spectral points
    !======================================================================
    idxEmRf = 2
    idxCloud = 2
    iTempLoc=paramIndex('Temp')

    NfsmpLoop: DO n=1, nNodes
      nn=iselS(n)
       !---Compute molecular optical depth for all atmospheric layers
      CALL OSStran(kfix_ir(1,nn),kh2o_ir(1,nn),dkh2o_ir(1,nn),kvar_ir(1,nn), &
            kself(:,nn),indxt_p1,indxt_p2,indxp,indxt, &
            ap1,ap2,at1_p1,at2_p1,at1_p2,at2_p2,at1,at2,  &
            adt1_p1,adt2_p1,adt1_p2,adt2_p2,N2,referenceGrid, &
            NmolS(nn),ImolS(1,nn),pavl,tavl,wfix,q,w,tautot,abso,dtaudtmp)

      vn=vWvn(nn)
      CALL vinterp(surfEmRf,surfEmRfGrid,vn,vEmrf,idxEmRf,coefInt)
      emissivity = vEmrf(1)
      reflectance= vEmrf(2)

      if (sun) then
        fbeam=sunrad(nn)
      else
        fbeam=0.
      end if
      if (scatterOn) then
        CALL optIntLayer(vn,layerNum+1,tabs,tscat,asym,jacAbs,jacScat,jacAsym, idxCloud)

        if (cloudMode==typeSlab) then
          do i= 1, layerNum
            cloudProfile%tau(i) = cloudProfile%molAbs(i)*tautot(cloudProfile%layerIndex(i))
          end do
          !---Perform RT calculations cloudy sky model

          CALL scattRT(cloudProfile%tau,cloudProfile%tempLevel,&
                xG(tSkinIndex),tabs,tscat,asym,emissivity,reflectance, &
                cosmicBackground, fBeam, vn, layerNum, obsLevelLoc, &
                sun,nStr,umu,umu0, lookup,azAngle,iflux,flxNetMono,&
                rad,pland,myLinInTau, lambertianLoc, tskin_a, emis_a, &
                temp_a,tautot_a,tabs_a,tscat_a,asym_a)

          xkt_tmp = 0.
          cldIndex =(/iCldLiq, iCldIce/)

          do jjHyd = 1, 2
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

          do jjHyd = 1, 2
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

             !          contribution from temperature variation
                        ! cloud top
              ratio = (cloudProfile%tempLevel(jj)-cloudProfile%tempLevel(jj-1)) &
                     /(cloudProfile%presLevel(jj)-cloudProfile%presLevel(jj-1))

              amount  = amount +ratio*temp_a(jj)

              idx=cldIndex(jjHyd)
              xkt_tmp(idx)=xkt_tmp(idx)+amount
            end if
          end do
          dummy(1:N2+1)=0.
          do jj=1,layerNum
            idx=cloudProfile%layerIndex(jj)
            dummy(idx) =dummy(idx) + tautot_a(jj)*cloudProfile%molAbs(jj)
          end do
          draddtemp(1:N2)=dummy(1:N2)*dtaudtmp(1:N2)

          jjJac=4
          do jjHyd = 1, 2
            do jj=1,N2
              idx=cloudProfile%layerIndex(jj)
              draddtemp(idx) = draddtemp(idx) + tabs_a(jj)*jacAbs(jj,jjJac,jjHyd) &
                                            + tscat_a(jj)*jacScat(jj,jjJac,jjHyd) + asym_a(jj)*jacAsym(jj,jjJac,jjHyd)
            end do
          end do
          !
          xkt_tmp(tempIndex+1:tempIndex+N2-1)=draddtemp(2:N2)*dtu(2:N2)+draddtemp(1:N2-1)*dtl(1:N2-1)
          xkt_tmp(tempIndex)=draddtemp(1)*dtu(1)
          xkt_tmp(tempIndex+N2)=draddtemp(N2)*dtl(N2)

          ! for WV
          ks = 1
          k=ImolS(kS,nn)  !variable gases indices ImolS(1,...,nmols)
          IXoff=varMolIndex(k)-1 !Pointer for output xkt_tmp
          ! contribution from fixed
          drdw(1:N2)=dummy(1:N2)*abso(0,1:N2)  !straight set in abso!
          xkt_tmp(IXoff+2:IXoff+N2)=drdw(2:N2)*dwqu(2:N2,0)+drdw(1:N2-1)*dwql(1:N2-1,0)
          !---level 1
          xkt_tmp(IXoff+1)=drdw(1)*dwqu(1,0)
          !---bottom level (N2+1)
          xkt_tmp(IXoff+N2+1)=drdw(N2)*dwql(N2,0)
          ! from WV
          drdw(1:N2)=dummy(1:N2)*abso(kS,1:N2)  !straight set in abso!
          xkt_tmp(IXoff+2:IXoff+N2)=xkt_tmp(IXoff+2:IXoff+N2)+drdw(2:N2)*dwqu(2:N2,k)+drdw(1:N2-1)*dwql(1:N2-1,k)
          !---level 1
          xkt_tmp(IXoff+1)=xkt_tmp(IXoff+1)+drdw(1)*dwqu(1,k)
          !---bottom level (N2+1)
          xkt_tmp(IXoff+N2+1)=xkt_tmp(IXoff+N2+1)+drdw(N2)*dwql(N2,k)

          DO kS=2,nmolS(nn)   !straight set of indices (1,...,nmols) for given wn
              k=ImolS(kS,nn)  !variable gases indices ImolS(1,...,nmols)
              drdw(1:N2)=-dummy(1:N2)*abso(k,1:N2)  !straight set in abso!
              xkt_tmp(IXoff+2:IXoff+N2)=xkt_tmp(IXoff+2:IXoff+N2) + &
                      drdw(2:N2)  *wvdwqu(2:N2,k) + &
                      drdw(1:N2-1)*wvdwql(1:N2-1,k)
              !---level 1
              xkt_tmp(IXoff+1)=xkt_tmp(IXoff+1)+drdw(1)*wvdwqu(1,k)
              !---bottom level (N2+1)
              xkt_tmp(IXoff+N2+1)=xkt_tmp(IXoff+N2+1)+drdw(N2)*wvdwql(N2,k)
          END DO

          DO kS=2,nmolS(nn)   !straight set of indices (1,...,nmols) for given wn
             k=ImolS(kS,nn)  !variable gases indices ImolS(1,...,nmols)
             IXoff=varMolIndex(k)-1 !Pointer for output xkt_tmp
             drdw(1:N2)=dummy(1:N2)*abso(kS,1:N2)  !straight set in abso!
             xkt_tmp(IXoff+2:IXoff+N2)=drdw(2:N2)*dwqu(2:N2,k)+drdw(1:N2-1)*dwql(1:N2-1,k)
             !---level 1
             xkt_tmp(IXoff+1)=drdw(1)*dwqu(1,k)
             !---bottom level (N2+1)
             xkt_tmp(IXoff+N2+1)=drdw(N2)*dwql(N2,k)
          END DO

          !   convert temp_a
          dummy=0.
          do jj=1, layerNum
            idx=cloudProfile%layerIndex(jj)
            dummy(idx)=dummy(idx) + temp_a(jj)*cloudProfile%ratio(jj)
            dummy(idx+1)=dummy(idx+1) + temp_a(jj)*(1.0 - cloudProfile%ratio(jj))
          end do
          dummy(N2+1)= temp_a(layerNum + 1)
          xkt_tmp(tempIndex:tempIndex+N2) = xkt_tmp(tempIndex:tempIndex+N2) +  dummy(1:N2+1)

        else
          if (present(xktcloud)) then
             CALL scattRT(tautot,xG,xG(tSkinIndex),tabs,tscat,asym,emissivity,reflectance, &
                  cosmicBackground, fBeam, vn, N2, obsLevelLoc, &
                  sun,nStr,umu,umu0, lookup,azAngle,iflux,flxNetMono,&
                  rad,pland,myLinInTau, lambertianLoc, tskin_a, emis_a, &
                  temp_a,tautot_a,tabs_a,tscat_a,asym_a)
          else
             CALL scattRT(tautot,xG,xG(tSkinIndex),tabs,tscat,asym,emissivity,reflectance, &
                  cosmicBackground, fBeam, vn, N2, obsLevelLoc, &
                  sun,nStr,umu,umu0, lookup,azAngle,iflux,flxNetMono,&
                  rad,pland,myLinInTau, lambertianLoc)
          end if

          !------ Jacobians for cloudscene
          if (present(xktcloud)) then
            !--- Jacobians for temperature profiles
            draddtemp(1:N2)=tautot_a(1:N2)*dtaudtmp(1:N2)

            do jjHyd = 1, size(cldScene)
              do jj =1, N2
                drdw(jj) = tabs_a(jj) *jacAbs (jj,iTempLoc,jjHyd) + &
                          tscat_a(jj) *jacScat(jj,iTempLoc,jjHyd) + &
                           asym_a(jj) *jacAsym(jj,iTempLoc,jjHyd)
              end do
              draddtemp(1:N2)=draddtemp(1:N2)+drdw(1:N2)
            end do

            xkt_tmp=0.

            xkt_tmp(tempIndex+1:tempIndex+N2-1)=draddtemp(2:N2)*dtu(2:N2)+draddtemp(1:N2-1)*dtl(1:N2-1)
            xkt_tmp(tempIndex)=draddtemp(1)*dtu(1)
            xkt_tmp(tempIndex+N2)=draddtemp(N2)*dtl(N2)
            xkt_tmp(tempIndex:tempIndex+N2) = xkt_tmp(tempIndex:tempIndex+N2) +  temp_a(1:N2+1)

                      ! for WV
            ks = 1
            k=ImolS(kS,nn)  !variable gases indices ImolS(1,...,nmols)
            IXoff=varMolIndex(k)-1 !Pointer for output xkt_tmp
            ! contribution from fixed
            drdw(1:N2)=tautot_a(1:N2)*abso(0,1:N2)  !straight set in abso!
            xkt_tmp(IXoff+2:IXoff+N2)=drdw(2:N2)*dwqu(2:N2,0)+drdw(1:N2-1)*dwql(1:N2-1,0)
            !---level 1
            xkt_tmp(IXoff+1)=drdw(1)*dwqu(1,0)
            !---bottom level (N2+1)
            xkt_tmp(IXoff+N2+1)=drdw(N2)*dwql(N2,0)
            ! from WV
            drdw(1:N2)=tautot_a(1:N2)*abso(kS,1:N2)  !straight set in abso!
            xkt_tmp(IXoff+2:IXoff+N2)=xkt_tmp(IXoff+2:IXoff+N2)+drdw(2:N2)*dwqu(2:N2,k)+drdw(1:N2-1)*dwql(1:N2-1,k)
            !---level 1
            xkt_tmp(IXoff+1)=xkt_tmp(IXoff+1)+drdw(1)*dwqu(1,k)
            !---bottom level (N2+1)
            xkt_tmp(IXoff+N2+1)=xkt_tmp(IXoff+N2+1)+drdw(N2)*dwql(N2,k)

            DO kS=2,nmolS(nn)   !straight set of indices (1,...,nmols) for given wn
                k=ImolS(kS,nn)  !variable gases indices ImolS(1,...,nmols)
                drdw(1:N2)=-tautot_a(1:N2)*abso(k,1:N2)  !straight set in abso!
                xkt_tmp(IXoff+2:IXoff+N2)=xkt_tmp(IXoff+2:IXoff+N2) + &
                        drdw(2:N2)  *wvdwqu(2:N2,k) + &
                        drdw(1:N2-1)*wvdwql(1:N2-1,k)
                !---level 1
                xkt_tmp(IXoff+1)=xkt_tmp(IXoff+1)+drdw(1)*wvdwqu(1,k)
                !---bottom level (N2+1)
                xkt_tmp(IXoff+N2+1)=xkt_tmp(IXoff+N2+1)+drdw(N2)*wvdwql(N2,k)
            END DO

            DO kS=2,nmolS(nn)   !straight set of indices (1,...,nmols) for given wn
               k=ImolS(kS,nn)  !variable gases indices ImolS(1,...,nmols)
               IXoff=varMolIndex(k)-1 !Pointer for output xkt_tmp
               drdw(1:N2)=tautot_a(1:N2)*abso(kS,1:N2)  !straight set in abso!
               xkt_tmp(IXoff+2:IXoff+N2)=drdw(2:N2)*dwqu(2:N2,k)+drdw(1:N2-1)*dwql(1:N2-1,k)
               !---level 1
               xkt_tmp(IXoff+1)=drdw(1)*dwqu(1,k)
               !---bottom level (N2+1)
               xkt_tmp(IXoff+N2+1)=drdw(N2)*dwql(N2,k)
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
                    end do
                  end do

            else if (cloudMode==typeLayer) then
            !+++++++++++++++++++++++++++++++++++++++++++++++++++++
            !!'layer'
            !++++++++++++++++++++++++++++++++++++
                do jjHyd = 1, size(cldScene)
                  do jjJac = 0, cldScene(1)%nProperty
                    if (jjJac==iTempLoc) cycle
                    do jj =1, N2
                      xktCloudTemp(jj,jjJac,jjHyd) = tabs_a(jj) *jacAbs (jj,jjJac,jjHyd) + &
                                  tscat_a(jj)*jacScat(jj,jjJac,jjHyd)  + &
                                   asym_a(jj) *jacAsym(jj,jjJac,jjHyd)
                    end do
                  end do
                end do
            else
                print *,'[ERR:', procName,']: wrong cloudModel type'
                call exit(1)
            end if
          end if
        end if

        xkt_tmp(tSkinIndex) = tskin_a
        nchLoop: DO ich=1, nch_ir(n)
          ich0 = ichMap_ir(ich,n)
          coefTmp = coef_ir(ich,n)
          if (jacRequest) then
            if (cloudMode==typeSlab) then
              do i=iCldLiq, iCldLiq+3
                xkt(i,ich0)=xkt(i,ich0) + xkt_tmp(i) * coefTmp
              end do
              do i=iCldIce, iCldIce+3
                xkt(i,ich0)=xkt(i,ich0) + xkt_tmp(i) * coefTmp
              end do
            else if (cloudMode==typeLevel) then
              do jjHyd = 1, size(cldScene)
                do jjJac = 0, cldScene(1)%nProperty
                  if (jjJac==iTempLoc) cycle
                  do jj =1, N2+1
                    xktCloud(jj,ich0,jjJac,jjHyd) = xktCloud(jj,ich0,jjJac,jjHyd) + &
                                               xktCloudTemp(jj,jjJac,jjHyd) * coefTmp
                  end do
                end do
              end do
            else if (cloudMode==typeLayer) then
              do jjHyd = 1, size(cldScene)
                do jjJac = 0, cldScene(1)%nProperty
                  if (jjJac==iTempLoc) cycle
                  do jj =1, N2
                    xktCloud(jj,ich0,jjJac,jjHyd) = xktCloud(jj,ich0,jjJac,jjHyd) + &
                                               xktCloudTemp(jj,jjJac,jjHyd) * coefTmp
                  end do
                end do
              end do
            end if
          end if
        END DO nchLoop

        !--- Jacobians for emissivity
        xkEmRf_tmp(1)= emis_a
        xkEmRf_tmp(2)=-emis_a

      else
        CALL OSSrad(jacRequest, tautot,abso,dtaudtmp,tavl,xG(tSkinIndex),dtu,dtl,&
                dwqu,dwql, wvdwqu,wvdwql,tempIndex,tSkinIndex,varMolIndex,plogu,&
                plogl,f1_arr(nn),f2_arr(nn),nmols(nn),imols(1,nn),xG,emissivity,&
                reflectance,fbeam,N1,N2,sun,umuLay,umu0Lay,umu0,lookup,&
                lambertianLoc,rad,xkt_tmp,xkEmRf_tmp)
      end if
      DO ich=1,nch_ir(n)
        ich0 = ichMap_ir(ich,n)
        coefTmp = coef_ir(ich,n)
        y(ich0)           = y(ich0)+rad*coefTmp
        if (.not. jacRequest) cycle
        xkt(tSkinIndex,ich0) = xkt(tSkinIndex,ich0)+xkt_tmp(tSkinIndex)*coef_ir(ich,n)
        xkt(tempIndex:tempIndex+N2,ich0) = xkt(tempIndex:tempIndex+N2,ich0) + &
                                           xkt_tmp(tempIndex:tempIndex+N2)*coefTmp
        DO k=1,NmolS(nn)
          ks=ImolS(k,nn)
          IXoff=varMolIndex(ks)
          xkt(IXoff:IXoff+N2,ich0)=xkt(IXoff:IXoff+N2,ich0) + xkt_tmp(IXoff:IXoff+N2) * coefTmp
        END DO
        DO i=1,2
         xx                    = xkEmRf_tmp(i)*coefTmp
         xkEmRf(idxEmRf-1,ich0, i) = xkEmRf(idxEmRf-1,ich0, i) + xx*(1.0-coefInt)
         xkEmRf(idxEmRf,  ich0, i) = xkEmRf(idxEmRf,  ich0, i) + xx*coefInt
        END DO
      END DO
    END DO NfsmpLoop
    RETURN
  END SUBROUTINE ossdrv_ir
!----------------------------------------------------------------------------
  ! PURPOSE: Interpolation in Log pressure
  !----------------------------------------------------------------------------
  SUBROUTINE lint_log(xInp,pGrid,N0,p0,x0)
    INTEGER           , INTENT(IN)   :: n0
    REAL              , INTENT(IN)   :: p0
    real, DIMENSION(:), INTENT(IN)   :: xInp,pGrid
    !---Output variable
    REAL              , INTENT(INOUT):: x0
    !---Local variable
    REAL              :: xx
    xx = log(xInp(N0)/xInp(N0-1))/ log(pGrid(N0)/pGrid(N0-1))
    x0 = xInp(N0-1)*(p0/pGrid(N0-1))**xx
    RETURN
  END SUBROUTINE lint_log

  !----
  SUBROUTINE optIntLayer(vn,nSurf,tabs,tscat,asym,       &
                       jacAbs, jacScat, jacAsym, idxCloud)

    !---Input variables
    REAL,               INTENT(IN)     :: vn
    INTEGER,            INTENT(IN)     :: nSurf
    INTEGER,            INTENT(INOUT)  :: idxCloud
    !---Output variables
    REAL, DIMENSION(:), INTENT(INOUT)  :: tabs,tscat,asym
    REAL, DIMENSION(:,0:,:), INTENT(INOUT)  :: jacAbs, jacScat, jacAsym
!    !---Local variables
    INTEGER            :: nlay,ip,ip1
    integer            :: jj
    real               :: ratio, invRatio
    integer            :: jjJac, jjHyd
    integer            :: nProperty, nHyd

    tabs = 0.
    asym = 0.
    tscat = 0.
    jacAbs = 0.
    jacScat = 0.
    jacAsym = 0.

    nProperty = size(tscatJac,dim=3)-1
    nHyd = size(tscatJac,dim=4)

    nlay = nSurf-1

    DO ip=idxCloud,nSpcPts-1
      IF(spectPtsGrd(ip) > vn) EXIT
    END DO

    idxCloud = ip
    ip1 = idxCloud-1
    ratio = (vn-spectPtsGrd(ip1))/(spectPtsGrd(idxCloud)-spectPtsGrd(ip1))

    if (ratio  < 0.) ratio = 0.
    if (ratio  > 1.) ratio = 1.

    invRatio = 1. - ratio

    do jj = 1, nlay
      tabs(jj) = invRatio*tabsDat(jj,ip1) + ratio*tabsDat(jj,idxCloud)
      tscat(jj) = invRatio*tscatDat(jj,ip1) + ratio*tscatDat(jj,idxCloud)
      if (tabs(jj) + tscat(jj) < epsilon(1.0) ) cycle

      asym(jj)  = invRatio*asymDat(jj,ip1)  + ratio*asymDat(jj,idxCloud)

      do jjHyd=1,nHyd
        do jjJac=0,nProperty
          jacAbs(jj,jjJac,jjHyd)  = invRatio*tabsJac(jj,ip1,jjJac,jjHyd)  + ratio*tabsJac(jj,idxCloud,jjJac,jjHyd)
          jacAsym(jj,jjJac,jjHyd) = invRatio*asymJac(jj,ip1,jjJac,jjHyd)  + ratio*asymJac(jj,idxCloud,jjJac,jjHyd)
          jacScat(jj,jjJac,jjHyd) = invRatio*tscatJac(jj,ip1,jjJac,jjHyd) + ratio*tscatJac(jj,idxCloud,jjJac,jjHyd)
        end do
      end do
    end do
    RETURN
  END SUBROUTINE optIntLayer

  SUBROUTINE OSSdestroyIR()
    use OSSIRmoduleSubs, only: OSSdestroyIRsubs
    call OSSdestroyIRsubs()
  END SUBROUTINE OSSdestroyIR

END MODULE oss_ir_module_scat
