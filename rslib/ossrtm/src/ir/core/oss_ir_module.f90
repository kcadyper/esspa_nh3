!--------------------------------------------------------------
!
!  MODULE OSS_IR_MODULE: contains items needed for the
!                OSS Forward Model (IR). AER Inc. 2004.
!
!--------------------------------------------------------------
MODULE oss_ir_module
  use params_module, ONLY: mxlev, mxlay, mxhmol, maxcmu, MxmolS, mxParG

  use OSSIRmoduleSubs, only : kself, kfix_ir, kvar_ir, kh2o_ir, dkh2o_ir,&
              fpathZ, fpathP, layerAverage, OSStran, sunrad, &
              myLinInTau, OSS_LIN_IN_TAU, sphericalGeometryFlag,&
              pi, nMol, UNDEFINED_INDEX,&
              zIntegrationFlag,nfsmp_ir,vWvn,&
              ImolS, NmolS, iselS, coef_ir, nch_ir, findFreeUnit, &
              f1_arr, f2_arr, ichMap_ir, rlin, setpath_ir, settabindx_ir,&
              setIndex_selfCont, vinterp, GetOD, planck

  IMPLICIT NONE
  PRIVATE
  !------------------------------------------------------------------------
  !	Public items made available to calling program(s)
  !------------------------------------------------------------------------
  PUBLIC :: ossinit_ir,ossdrv_ir,ossRad, OSSdestroyIR

CONTAINS
  !----------------------------------------------------------------------------
  !    PURPOSE: Performs the 1st interface  between the module and the
  !    calling program.
  !    It reads in the unit numbers, the files names, the geophysical vector
  !    pointers and returns back information that is useful to the user,
  !    pertaining to the wavenumber, the pressure grid, etc.
  !  selFile instrument chanel file
  !  lutFile  optical properties LUT file
  !  solFile solar spectrum file
  !  defProfFile standard profile file
  !  nVarMol number of variable molecules
  !  varMolID variable gas HITRAN indices
  !  chSelFile - chanel selection file name
  !  nChanOut - number of output channels
  !  chanIndex - (output) channel indices
  !  chanFreq  - output channel frequencies
  !  linInTau  - if to use linearInTau approximation
  !
  !----------------------------------------------------------------------------
  SUBROUTINE ossinit_ir(selFile,lutFile,solFile, &
       defProfFile,nVarMol,varMolID,nChanOut, chanIndex, chanFreq, &
       nRef,pRefOut,sphericalGeometry, zIntegration,linInTau)

    USE asolar_source_function, only : getSolar, InterpSolar
    use OSSIRmoduleSubs, only : getReferencePressure, getChannelFrequency,&
                                getCountChannel, getChannelIndex,&
                                getCurrentSelection

    !---Input variables
    CHARACTER(len=*), INTENT(IN) :: selFile,lutFile,solFile,defProfFile
    INTEGER,          INTENT(IN) :: nVarMol
    INTEGER, DIMENSION(:),          INTENT(IN) :: varMolID
    LOGICAL,          OPTIONAL,INTENT(IN) :: linInTau
    logical,          optional,intent(in) :: sphericalGeometry
    logical,          optional,intent(in) :: zIntegration

    !---Output variables
    INTEGER                            ,INTENT(INOUT) :: nChanOut
    real, DIMENSION(:)                 ,INTENT(INOUT) :: chanFreq
    INTEGER                            ,INTENT(INOUT) :: nRef
    integer, allocatable, dimension(:), intent(inout) :: chanIndex
    real,    allocatable, DIMENSION(:), INTENT(INOUT) :: pRefOut
    !local
    INTEGER                             :: iuSol
    integer                             :: selIndex

    !------------------------------------------------------------------------
    ! Set global variables from file variables
    !------------------------------------------------------------------------
    !---read solar file
    iuSol = findFreeUnit()
    call getSolar(iuSol,solFile)

    !---read OD and SEL data
     CALL GetOD(selFile,lutFile,defProfFile,nVarMol,varMolID)

    ! 0th channel selection index assumed
    selIndex = getCurrentSelection()
    nChanOut = getCountChannel(selIndex)

     !---Size conformity checks (safer interface with calling program/subroutine)
    IF (SIZE(chanFreq) < nChanOut) THEN
        print*, 'Err[oss_ir_module::oss_init_ir]: Insufficient Size for CWVN_OUT'
        call exit(1)
    END IF

    IF (SIZE(varMolID)  > mxHmol) THEN
       print*, 'Err[oss_ir_module::oss_init_ir]: MolID Vector too large'
       call exit(1)
    END IF

    IF (present(linInTau)) THEN
       myLinInTau = linInTau
    ENDIF

    if (present(zIntegration)) THEN
       zIntegrationFlag = zIntegration
    else
       zIntegrationFlag = .false.
    end if

    if (present(sphericalGeometry)) THEN
       sphericalGeometryFlag = sphericalGeometry
    else
       sphericalGeometryFlag = .false.
    end if

    !---interpolate solar
    if (allocated(sunrad)) deallocate (sunrad)
    allocate (sunrad(nfsmp_ir))
    call InterpSolar(NFSmp_ir,vwvn,Sunrad)

    !------------------------------------------------------------------------
    ! Copy to output arguments
    !------------------------------------------------------------------------
    if (.not.allocated(chanIndex)) allocate(chanIndex(nChanOut))

    call getChannelFrequency(chanFreq,selIndex)
    call getChannelIndex(chanIndex,selIndex)

    call getReferencePressure(pRefOut)
    nRef = size(pRefOut)

    RETURN
  END SUBROUTINE ossinit_ir

  !----------------------------------------------------------------------------
  ! PURPOSE: Driver for the IR radiative transfer model.
  !          Computes both radiances and their jacobians wrt geophysical
  !          parameters.
  ! xG:      Profile vector of geophysical parmaters (containing
  !          temperature and constiutuents profiles, surface
  !          pressure and skin temperature, and cloud parameters
  ! surfEmRf:    vector of MW surface emissivities.
  ! Path:    Array of quantities characterizing the atmospehric
  !          path in plane-parallel atmosphere.
  ! Y:       Vector of radiances.
  ! xkt:     Array of derivatives of the radiances wrt geophysical parameters.
  ! xkEmRf:  Array of radiance derivatives wrt surface emissivities.
  ! lat observer latitude (affects garvity acceleration calculation
  ! tempIndex,tSkinIndex,pSurfIndex,varMolIndex - indices of T, Tskin, Psurf, and Q in xG
  ! zProf  - altitude grid (optional)
  ! pUser    - pressure grid (optional)
  ! lambertian - if set to 0, specular reflection, otherwise, a diffusion approximation is used
  !              to calculate downwelling IR reflection from the surface
  !----------------------------------------------------------------------------
  SUBROUTINE ossdrv_ir(xG,surfEmRfGrid,surfEmRf,obsAngle,solZenith,azAngle,&
            obsLevel,y,xkt,xkEmRf,lat,tempIndex,tSkinIndex,pSurfIndex,&
            varMolIndex,lambertian, zSurf,zProf,pUser)
    USE OSSIRmoduleSubs, only : getCountChannel, getCountUsedNode

    !---Input variables
    real, DIMENSION(:)           , INTENT(IN) :: xG
    real, DIMENSION(:)           , INTENT(IN) :: surfEmRfGrid
    real, DIMENSION(:,:)         , INTENT(IN) :: surfEmRf
    REAL                         , INTENT(IN) :: obsAngle,solZenith,azAngle
    REAL                         , INTENT(IN) :: lat

    INTEGER, optional            , INTENT(IN) :: lambertian
    REAL              ,  OPTIONAL, INTENT(IN) :: zSurf
    real, DIMENSION(:),  OPTIONAL, INTENT(IN) :: pUser
    real, DIMENSION(:),  OPTIONAL, INTENT(IN) :: zProf

    INTEGER                      , INTENT(IN) :: tempIndex,tSkinIndex,pSurfIndex
    INTEGER,DIMENSION(:)         , INTENT(IN) :: varMolIndex
    !---Output variables
    INTEGER           ,            INTENT(IN) :: obsLevel
    real, DIMENSION(:)        , INTENT(INOUT) :: y
    real, DIMENSION(:,:)      , INTENT(INOUT) :: xkt
    real, DIMENSION(:,:,:)    , INTENT(INOUT) :: XkEmRf
    !---Local variables
    LOGICAL                          :: lookup,sun,referenceGrid,interpSfc
    INTEGER                          :: nsf,nparG
    real,    DIMENSION(Mxlev)        :: pLoc
    real,    DIMENSION(MxLay)        :: tavl,pavl,wfix,dtu,dtl,ap1,ap2
    real,    DIMENSION(MxHmol,MxLay) :: q,w
    real,    DIMENSION(MxLay,0:MxHmol) :: dwqu,dwql
    real,    DIMENSION(MxLay,MxHmol) :: wvdwqu,wvdwql

    INTEGER, DIMENSION(MxLay)        :: indxt_p1,indxt_p2,indxp,indxt
    real,    DIMENSION(MxLay)        :: at1_p1,at2_p1,at1_p2,at2_p2,at1,at2
    real,    DIMENSION(MxLay)        :: adt1_p1,adt2_p1,adt1_p2,adt2_p2
    real,    DIMENSION(MxLay)        :: tautot,dtaudtmp
    real,    DIMENSION(0:MxmolS,MxLay) :: abso
    real,    DIMENSION(MxLay)        :: umuLay,umu0Lay, alt
    REAL                             :: xkt_tmp(MxParG),xkEmRf_tmp(2)
    REAL                             :: umu,umu0,vn,fbeam,xx,rad
    REAL                             :: coefInt,Tsfc,plogu,plogl
    INTEGER                          :: n1,n2,ip0,n,nn,ich,i,ich0,nSurf
    INTEGER                          :: k,ks,ixOff
    INTEGER                          :: nChSel
    INTEGER                          :: nNodes
    REAL                             :: vEmRf(size(surfEmRf, dim=2))
    real                             :: reflectance
    real                             :: emissivity
    integer                          :: nEnd
    integer                          :: obsLevelLoc
    real                             :: alpha
    logical                          :: jacRequest = .true.
    logical                          :: lambertianLoc
    character(len=*), parameter      :: errHeader='Err[oss_ir_module::ossdrv]:'


    nChSel=getCountChannel()
    nNodes=getCountUsedNode()
    obsLevelLoc=obsLevel
    nsf = size(surfEmRfGrid)

    if (size(surfEmRf,1) < nsf) then
       print*,errHeader,' Hinge point dimension of surfEmRf is too small'
       call exit(1)
    endif
    if ( jacRequest) then
      if (size(XkEmRf,1) < nsf) then
         print*,errHeader,' Hinge point dimension of XkEmRf is too small'
         call exit(1)
      endif
      if (size(XkEmRf,2) < nchSel) then
         print*,errHeader,' Channel dimension of XkEmRf is too small'
         call exit(1)
      endif
      if (size(xkt,1)<size(xG)) then
         print*,errHeader,' Parameter dimension of xkt is too small'
         call exit(1)
      end if
      if (size(xkt,2)<nchSel) then
         print*,errHeader,' Channel dimension of xkt is too small'
         call exit(1)
      end if
      if (size(xkt,2)<nchSel) then
         print*,errHeader,' Channel dimension of xkt is too small'
         call exit(1)
      endif
    end if

    IF (SIZE(varMolIndex)> mxHmol) THEN
       print*, errHeader,' Imolind Vector too large'
       call exit(1)
    END IF
    if (zIntegrationFlag) then
      if ( .not. (present(zSurf) .and. present(zProf))) then
         print*, errHeader,' zIntegration requires zSurf and zProf ',&
                           'to be provided'
         call exit(1)
      end if
    end if
    if (sphericalGeometryFlag) then
      if ( .not. (present(zSurf) .and. present(zProf))) then
         print*, errHeader,' Spherical geometry requires zSurf and zProf ', &
                           'to be provided'
         call exit(1)
      end if
    end if
    if (present(lambertian)) THEN
      lambertianLoc=lambertian>0
      else
      lambertianLoc=.false.
    end if

    nparG = varMolIndex(nMol)-1 + tSkinIndex-tempIndex
    !======================================================================
    !     Initialize radiance vector and k-matrix
    !======================================================================
    y(1:NchSel)               = 0.

    if (jacRequest) then
      DO ich=1,NchSel
          xkt(:,ich)      = 0.
          xkEmRf(1:Nsf,ich, 1:2) = 0.
      END DO
    end if

    !======================================================================
    !     Compute path variables
    !======================================================================
    if(present(pUser))then
       if (size(pUser).gt.mxlev) then
          print*,errHeader,' Input pressure grid is too large.'
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

    if (interpSfc) then
        if (present(zProf) .AND. (.NOT. present(zSurf)) ) then
           print*,errHeader, ' zProf requires zSurf to be presented as well'
           call exit(1)
        end if
    end if

    n2=nSurf-1
    n1=obsLevelLoc

    if ( present(zProf) .and. present(zSurf)) then
         alt(1:N2)=zProf(1:N2)
         alt(N2+1)=zSurf
    else
         alt(1:N2+1)=0.
    end if

    !---Compute average temperature, pressure and integrated
    !   molecular amounts for the layers
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
    !---Compute coefficients for temperature interpolation of ODs
    CALL settabindx_ir(tavl,pavl,N2,referenceGrid, &
         indxt_p1,indxt_p2,indxp,ap1,ap2,&
         at1_p1,at2_p1,at1_p2,at2_p2,adt1_p1,adt2_p1,adt1_p2,adt2_p2)
    CALL setIndex_selfCont(tavl,N2,indxt,at1,at2)
    !======================================================================
    !     Loop over spectral points
    !======================================================================
    ip0 = 2
    NfsmpLoop: DO n = 1, nNodes
       nn=iselS(n)
       !---Compute molecular optical depth for all atmospheric layers
       CALL OSStran(kfix_ir(1,nn),kh2o_ir(1,nn),dkh2o_ir(1,nn),kvar_ir(1,nn), &
            kself(:,nn),indxt_p1,indxt_p2,indxp,indxt, &
            ap1,ap2,at1_p1,at2_p1,at1_p2,at2_p2,at1,at2,  &
            adt1_p1,adt2_p1,adt1_p2,adt2_p2,N2,referenceGrid, &
            NmolS(nn),ImolS(1,nn),pavl,tavl,wfix,q,w,tautot,abso,dtaudtmp)
       !---Interpolate input surface emissivity to node wavenumber
       vn=vWvn(nn)
       CALL vinterp(surfEmRf,surfEmRfGrid,vn,vEmrf,ip0,coefInt)
       ! lambertian
       emissivity = vEmrf(1)
       reflectance = vEmrf(2)

       IF (sun) then
         fbeam=sunrad(nn)
       else
          fbeam=0.
       end if
       !---Perform RT calculations !clear sky model

       CALL OSSrad(jacRequest, tautot,abso,dtaudtmp,tavl,Tsfc,dtu,dtl,dwqu,dwql,&
            wvdwqu, wvdwql,tempIndex,tSkinIndex,varMolIndex,plogu,plogl,&
            f1_arr(nn),f2_arr(nn),nmols(nn),imols(1,nn),xG,emissivity,&
            reflectance,fbeam,N1,N2,sun,umuLay,umu0Lay,umu0,lookup,lambertianLoc,&
            rad,xkt_tmp,xkEmRf_tmp)

       nchLoop: DO ich=1,nch_ir(n)
         ich0 = ichMap_ir(ich,n)
         y(ich0)           = y(ich0)+rad*coef_ir(ich,n)
         if (.not. jacRequest) cycle
         xkt(tempIndex:tempIndex+N2,ich0) = xkt(tempIndex:tempIndex+N2,ich0) + &
            xkt_tmp(tempIndex:tempIndex+N2)*coef_ir(ich,n)
         xkt(tSkinIndex,ich0) = xkt(tSkinIndex,ich0)+xkt_tmp(tSkinIndex)*coef_ir(ich,n)

         DO k=1,NmolS(nn)
            ks=ImolS(k,nn)
            IXoff=varMolIndex(ks)-1
            xkt(IXoff+1:IXoff+N2+1,ich0)=xkt(IXoff+1:IXoff+N2+1,ich0)+&
                                xkt_tmp(IXoff+1:IXoff+N2+1)*coef_ir(ich,n)
         END DO

          DO i=1,2
            xx                   = xkEmRf_tmp(i)*coef_ir(ich,n)
            xkEmRf(ip0-1,ich0, i) = xkEmRf(ip0-1,ich0, i) + xx*(1.0-coefInt)
            xkEmRf(ip0,  ich0, i) = xkEmRf(ip0,  ich0, i) + xx*coefInt
          END DO
      END DO nchLoop
    END DO NfsmpLoop
    RETURN
  END SUBROUTINE ossdrv_ir

  !----------------------------------------------------------------------------
  ! PURPOSE: Compute radiance (in mw/m2/str/cm-1) and derivatives of radiances
  !          with respect to atmospheric and surface parameters
  !------------------------------------------------scal----------------------------
  SUBROUTINE ossrad(jacRequest,tautot,abso,dtaudtmp,tavl,tsfc,dtu,dtl,dwqu,dwql, &
         wvdwqu,wvdwql,tempIndex,tSkinIndex,varMolIndex,plogu,plogl,f1,f2, &
         nmols,imols,xG,emissivity, reflectance,bs_sun,N1,N2,sun,&
         umuLay,umu0Lay,umu0, lookup,lambertian,rad,xkt,xkEmRf)
    !Parameters
    real,    PARAMETER     :: LAMBERTIAN_REFL_SECANT=1.66
    !---Input variables
    logical,                       intent(in) :: jacRequest
    LOGICAL             ,INTENT(IN) :: sun,lookup
    logical             ,INTENT(IN) :: lambertian
    INTEGER             ,INTENT(IN) :: N1,N2
    REAL                ,INTENT(IN) :: bs_sun,Tsfc
    INTEGER             ,INTENT(IN) :: NmolS,ImolS(NmolS)
    real, DIMENSION(:)  ,INTENT(IN) :: xG,tautot,dtaudtmp,tavl,dtu,dtl
    real, DIMENSION(0:,:),INTENT(IN):: abso
    real, DIMENSION(:,0:),INTENT(IN):: dwqu,dwql
    real, DIMENSION(:,:),INTENT(IN) :: wvdwqu,wvdwql
    real, DIMENSION(:)  ,INTENT(IN) :: umuLay,umu0Lay
    INTEGER             ,INTENT(IN) :: tempIndex,tSkinIndex
    INTEGER,DIMENSION(:),INTENT(IN) :: varMolIndex
    real                        :: emissivity, reflectance
    !---Output variables:
    REAL              ,INTENT(INOUT):: rad
    real, DIMENSION(:),INTENT(INOUT):: xkt,xkEmRf
    !---Local variables:
    INTEGER                :: N2prim
    INTEGER                :: l,ks,k,ixoff
    REAL                   :: f1,f2
    REAL                   :: bs,dbs,sumtau_dwn
    REAL                   :: tausfc,radsun,draddrsfc,draddemis
    REAL                   :: txsun,sumtau0_up,umu0
    REAL                   :: draddtskn,sumtau_up,dtran,rsfc
    REAL                   :: dbavgdb,dbavgdbdod
    REAL                   :: odsec,rlt,plogu,plogl,dTsfcdtu,dTsfcdtl
    real, DIMENSION(MxLev) :: txdn,txup,bbar,dbbar,draddtmp,draddtau,drdw
    real, DIMENSION(MxLev) :: blev,dblev,draddtmpdw,draddtmpuw,xGtmp
    real, DIMENSION(Mxlev) :: secRefl,sec,sec0,draddtau_sun
    real                   :: localRad
    rsfc = 0.0
    dTsfcdtl = 0.0
    dTsfcdtu = 0.0
    sumtau_dwn = 0.0
    tausfc = 1.0
    sumtau0_up = 0.0

    sec(1:N2)      = 1./umuLay(1:N2)
    sec0(1:N2)     = 1./umu0Lay(1:N2)
    IF (lambertian) then
       secRefl(1:N2)= LAMBERTIAN_REFL_SECANT
    ELSE
       secRefl(1:N2)= sec(1:N2)
    END IF
  !-----------------------------------------------------------------------
  !     Compute Planck function and its derivative wrt temperature
  !-----------------------------------------------------------------------

    if (myLinInTau) then
       xGtmp(1:n2)=xG(1:n2)
       xGtmp(n2+1)=tsfc
       DO l=1,n2+1
          CALL planck(f1,f2,xGtmp(l),blev(l),dblev(l))
       END DO
       call planck(f1,f2,tSfc,blev(n2+1),dblev(n2+1))
       dtsfcdtu=tsfc/xG(n2)*plogl/(plogl-plogu)
       dtsfcdtl=tsfc/xG(n2+1)*plogu/(plogu-plogl)
    endif
    DO l=1,n2
        CALL planck(f1,f2,tavl(l),bbar(l),dbbar(l))
    END DO
    CALL planck(f1,f2,xG(tSkinIndex),bs,dbs)

    !-----------------------------------------------------------------------
    !     Compute transmittance profile along viewing path down to surface
    !-----------------------------------------------------------------------
    if ( .not.lookup) then
       txdn(1:N1)        = 1.0
       sumtau_dwn        = 0.
       DO l=N1,N2
          sumtau_dwn     = sumtau_dwn+tautot(l)*sec(l)
          txdn(l+1)      = EXP(-sumtau_dwn)
       END DO
       tausfc            = txdn(N2+1)
    endif
    !---Initialize radiance and derivative arrays
    rad               = 0.
    radsun            = 0.
    draddrsfc         = 0.
    draddtau_sun(1:N2)= 0.
    draddtmp(1:N2)    = 0.
    draddtau(1:N2)    = 0.
    draddtmpdw(1:N2+1)= 0.
    draddtmpuw(1:N2+1)= 0.
    draddemis         = 0.
    draddtskn         = 0.
    !-----------------------------------------------------------------------
    !     1- Downwelling thermal radiance calculation:
    !-----------------------------------------------------------------------
    IF (lookup) then
        ! up-looking
       txup(N1:N2+1)  = 1.0
       sumtau_dwn     = 0.0
       sumtau_up      = 0.0
       sumtau0_up     = 0.0
       N2prim         = min(N1-1,N2)

       DO l=N2prim,1,-1
          sumtau_up  = sumtau_up+tautot(l)*sec(l)
          txup(l)     = EXP(-sumtau_up)
       END DO
       if (myLinInTau) then
          DO l=1,N2prim
             dtran       = txup(l+1)-txup(l)
             odsec=tautot(l)*sec(l)
             call rlin(odsec,dbavgdb,dbavgdbdod)
             if (OSS_LIN_IN_TAU) then
                rlt=blev(l+1)*(1-dbavgdb)+blev(l)*dbavgdb
                draddtmpuw(l)=draddtmpuw(l)+dtran*dblev(l)*dbavgdb
                draddtmpuw(l+1)=dtran*dblev(l+1)*(1.0-dbavgdb)
                draddtau(l)=(txup(l)*rlt - rad + &
                     dtran*dbavgdbdod*(blev(l)-blev(l+1)))*sec(l)
                draddtmp(l) = 0.0
             else
                rlt=2.0*bbar(l)*dbavgdb+blev(l+1)*(1.0-2.0*dbavgdb)
                draddtmpuw(l+1)=dtran*dblev(l+1)*(1.0-2.0*dbavgdb)
                draddtmp(l)=dtran*2.0*dbbar(l)*dbavgdb
                draddtau(l)=(txup(l)*rlt - rad + &
                     2.0*dtran*dbavgdbdod*(bbar(l)-blev(l+1)))*sec(l)
             end if
             rad=rad+dtran*rlt
          END DO
       else
          DO l=1,N2prim
             dtran       = txup(l+1)-txup(l)
          !---Derivative of upwelling emission wrt to temperature and optical thickness
             draddtau(l) = (txup(l)*bbar(l)-rad)*sec(l)
             draddtmp(l) = dtran*dbbar(l)
             rad         = dtran*bbar(l) + rad
          ENDDO
       endif
       DO l=1,N2prim
          draddtmp(l) = draddtmp(l) + draddtau(l)*dtaudtmp(l)
       end do
    ELSE
        ! downlooking
       IF (tausfc>1.e-06) THEN
          txup(N2+1)     = tausfc
          sumtau_up     = 0.
          DO l=N2,1,-1
             sumtau_up  = sumtau_up+tautot(l)*secRefl(l)
             sumtau0_up  = sumtau0_up+tautot(l)*sec0(l)
             txup(l)     = EXP(-(sumtau_dwn+sumtau_up))
          END DO

          if (myLinInTau) then
             DO l=1,N2
                dtran       = txup(l+1)-txup(l)
                odsec=tautot(l)*secRefl(l)
                call rlin(odsec,dbavgdb,dbavgdbdod)
                if (OSS_LIN_IN_TAU) then
                  rlt=blev(l+1)+(blev(l)-blev(l+1))*dbavgdb
                  draddtmpdw(l)=draddtmpdw(l)+dtran*dblev(l)*dbavgdb
                  draddtmpdw(l+1)=dtran*dblev(l+1)*(1.-dbavgdb)
                      draddtmp(l) =0.0
                  draddtau(l)=(txup(l)*rlt-rad+dtran*dbavgdbdod* &
                       (blev(l)-blev(l+1)))*secRefl(l)
                else
                    rlt=2.0*bbar(l)*dbavgdb+blev(l+1)*(1.0-2.0*dbavgdb)
                    draddtmpdw(l+1)=dtran*dblev(l+1)*(1.0-2.0*dbavgdb)
                    draddtmp(l)=dtran*2.0*dbbar(l)*dbavgdb
                    draddtau(l)=(txup(l)*rlt-rad+dtran*dbavgdbdod* &
                                        2.0*(bbar(l)-blev(l+1)))*secRefl(l)
                end if
                rad=rad+dtran*rlt

             ENDDO
             draddtmpdw(N2)=draddtmpdw(N2)+draddtmpdw(N2+1)*dtsfcdtu
             draddtmpdw(N2+1)=draddtmpdw(N2+1)*dtsfcdtl
          else
             DO l=1,N2
                dtran       = txup(l+1)-txup(l)
          !---Derivative of downwelling emission wrt to temperature and optical thickness
                draddtau(l) = (txup(l)*bbar(l)-rad)*secRefl(l)
                draddtmp(l) = dtran*dbbar(l)
                rad         = rad+dtran*bbar(l)
             ENDDO
          endif
       END IF
    !-----------------------------------------------------------------------
    !      Adjust on surface reflectivity
    !-----------------------------------------------------------------------
       rsfc              = (1.0-emissivity)
       DO l=1,N2
          draddtau(l) = draddtau(l)*rsfc
          draddtmp(l) = draddtmp(l)*rsfc
       END DO
    !-----------------------------------------------------------------------
    !     2- Add surface component
    !-----------------------------------------------------------------------
    !---Derivatives wrt emissivity and sfc skin temperature:
       draddemis         = tausfc*bs-rad
       draddtskn         = emissivity*tausfc*dbs
       localRad          = emissivity*tausfc*bs
       rad               = rad*rsfc +localRad

    !-----------------------------------------------------------------------
    !     3- Add solar component
    !-----------------------------------------------------------------------
        IF(sun)THEN
           txsun          = EXP(-(sumtau_dwn+sumtau0_up))
           !---Derivatives wrt sfc solar reflectance
           draddrsfc      = txsun*bs_sun*umu0/pi
       localRad         = reflectance*draddrsfc
       DO l=1,N2
          draddtau(l) = draddtau(l) - localRad*sec0(l)
       END DO
       rad            = rad + localRad
    endif

    !     4- Upwelling thermal radiance calculation
    !-----------------------------------------------------------------------
       if (myLinInTau) then
          DO l=N2,N1,-1
             dtran=txdn(l)-txdn(l+1)
             odsec=sec(l)*tautot(l)
             call rlin(odsec,dbavgdb,dbavgdbdod)

             if (OSS_LIN_IN_TAU) then
               rlt=blev(l)*(1.0-dbavgdb)+blev(l+1)*dbavgdb
               draddtmpuw(l)=dtran*dblev(l)*(1.0-dbavgdb)
               draddtmpuw(l+1)=draddtmpuw(l+1)+dtran*dblev(l+1)*dbavgdb
                 draddtau(l)=draddtau(l) + (txdn(l+1)*rlt-rad + &
                           dtran*dbavgdbdod*(blev(l+1)-blev(l)))*sec(l)
             else
                 rlt=2.0*bbar(l)*dbavgdb+blev(l)*(1.0-2.0*dbavgdb)
                 draddtmpuw(l)=dtran*dblev(l)*(1.0-2.0*dbavgdb)
                 draddtmp(l)=draddtmp(l) + dtran*2.0*dbbar(l)*dbavgdb
                 draddtau(l)=draddtau(l) + (txdn(l+1)*rlt-rad + &
                             dtran*dbavgdbdod*2.0*(bbar(l)-blev(l)))*sec(l)
             end if

             rad=rad+dtran*rlt
          ENDDO
          draddtmpuw(N2)=draddtmpuw(N2)+draddtmpuw(N2+1)*dtsfcdtu
          draddtmpuw(N2+1)=draddtmpuw(N2+1)*dtsfcdtl
       else
          DO l=N2,N1,-1
             dtran          = txdn(l)-txdn(l+1)
             draddtau(l)    = draddtau(l) + (txdn(l+1)*bbar(l)-rad)*sec(l)
             draddtmp(l)    = draddtmp(l) + dtran*dbbar(l)
             rad            = rad+dtran*bbar(l)
          ENDDO
       endif
       DO l=1,N2
           draddtmp(l) = draddtmp(l) + draddtau(l)*dtaudtmp(l)
       END DO
    endif

    xkt=0.
    if (jacRequest) then
      !-----------------------------------------------------------------------
      !     Compute level derivatives and and map to array XKT:
      !-----------------------------------------------------------------------
      !  Air Temperature
      xkt(tempIndex+1:tempIndex+N2-1)=draddtmp(2:N2)*dtu(2:N2)+draddtmp(1:N2-1)*dtl(1:N2-1)
      !---level 1
      xkt(tempIndex)=draddtmp(1)*dtu(1)
      !---bottom level (N2+1)
      xkt(tempIndex+N2)=draddtmp(N2)*dtl(N2)
      if (myLinInTau) then
         xkt(tempIndex+1:tempIndex+N2-1)=xkt(tempIndex+1:tempIndex+N2-1)+ &
             draddtmpuw(2:N2)+draddtmpdw(2:N2)*rsfc
         xkt(tempIndex)=xkt(tempIndex)+draddtmpuw(1)+draddtmpdw(1)*rsfc
         xkt(tempIndex+N2)=xkt(tempIndex+N2)+draddtmpuw(N2+1)+draddtmpdw(N2+1)*rsfc
      endif
      !---Molecular concentrations
      ks = 1
      k=ImolS(kS)  !variable gases indices ImolS(1,...,nmols)
      IXoff=varMolIndex(k)-1 !Pointer for output xkt
      ! contribution from fixed
      drdw(1:N2)=draddtau(1:N2)*abso(0,1:N2)  !straight set in abso!
      xkt(IXoff+2:IXoff+N2)=drdw(2:N2)*dwqu(2:N2,0)+drdw(1:N2-1)*dwql(1:N2-1,0)
      !---level 1
      xkt(IXoff+1)=drdw(1)*dwqu(1,0)
      !---bottom level (N2+1)
      xkt(IXoff+N2+1)=drdw(N2)*dwql(N2,0)

      drdw(1:N2)=draddtau(1:N2)*abso(kS,1:N2)  !straight set in abso!
      xkt(IXoff+2:IXoff+N2)=xkt(IXoff+2:IXoff+N2)+drdw(2:N2)*dwqu(2:N2,k)+drdw(1:N2-1)*dwql(1:N2-1,k)
      !---level 1
      xkt(IXoff+1)=xkt(IXoff+1)+drdw(1)*dwqu(1,k)
      !---bottom level (N2+1)
      xkt(IXoff+N2+1)=xkt(IXoff+N2+1)+drdw(N2)*dwql(N2,k)

      DO kS=2,nmolS   !straight set of indices (1,...,nmols) for given wn
          l=ImolS(kS)  !variable gases indices ImolS(1,...,nmols)
          drdw(1:N2)= -draddtau(1:N2)*abso(kS,1:N2)  !straight set in abso!
          xkt(IXoff+2:IXoff+N2)=xkt(IXoff+2:IXoff+N2) + &
                  drdw(2:N2)  *wvdwqu(2:N2,l) + &
                  drdw(1:N2-1)*wvdwql(1:N2-1,l)
          !---level 1
          xkt(IXoff+1)=xkt(IXoff+1)+drdw(1)*wvdwqu(1,l)
          !---bottom level (N2+1)
          xkt(IXoff+N2+1)=xkt(IXoff+N2+1)+drdw(N2)*wvdwql(N2,l)
      END DO

      DO kS=2,nmolS   !straight set of indices (1,...,nmols) for given wn
         k=ImolS(kS)  !variable gases indices ImolS(1,...,nmols)
         IXoff=varMolIndex(k)-1 !Pointer for output xkt
         drdw(1:N2)=draddtau(1:N2)*abso(kS,1:N2)  !straight set in abso!
         xkt(IXoff+2:IXoff+N2)=drdw(2:N2)*dwqu(2:N2,k)+drdw(1:N2-1)*dwql(1:N2-1,k)
         !---level 1
         xkt(IXoff+1)=drdw(1)*dwqu(1,k)
         !---bottom level (N2+1)
         xkt(IXoff+N2+1)=drdw(N2)*dwql(N2,k)
      END DO

      if (lookup) then
         xkt(tSkinIndex)=0.
         xkEmRf=0.
      else
         !---Surface terms
         xkt(tSkinIndex)=draddtskn
         xkEmRf(1)=draddemis
         xkEmRf(2)=draddrsfc
      endif
    end if
    RETURN
  END SUBROUTINE ossrad

  !-----------------------------------------------------------
  SUBROUTINE OSSdestroyIR()
    use OSSIRmoduleSubs, only: OSSdestroyIRsubs
    call OSSdestroyIRsubs()
  END SUBROUTINE OSSdestroyIR

END MODULE oss_ir_module
