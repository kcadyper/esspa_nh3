!--------------------------------------------------------------
!
!  MODULE BkgPrepModule: contains procedures related to processing
!          retrieval background data (a priori and error covariance)
!
!--------------------------------------------------------------
MODULE BkgPrepModule
  USE StateIndexModule
  USE ReadStdInputs
  USE VertCoord, ONLY: mxCoordTyp,Pcoord,Scoord,Hcoord, &
     putSigmDefin,putHybrDefin,getSigmPres,getHybrPres
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------
  !  List of Public subroutines (accessible from outside module) 
  !---------------------------------------------------------------
  PUBLIC :: loadBkgClim,setFOVbkg,combBkg,tuneBkg,transBkg,emisChgTest,closeBkg
  !---------------------------------------------------------------
  !  List of Public variables (accessible from outside module) 
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  !     Declarations of private data (private to the module)
  !---------------------------------------------------------------
  INCLUDE 'mxdims.incl'
  INTEGER, PARAMETER                      :: mxbakSfc=8,mxbakAtm=8
  INTEGER, DIMENSION(mxbakSfc)            :: SfcTypes
  INTEGER, DIMENSION(mxbakAtm)            :: AtmTypes
  INTEGER, DIMENSION(MxPar)               :: extrvec
  INTEGER                                 :: nSfcTypes,nAtmTypes,ih2o,nfor
  INTEGER                                 :: nParMw,nparGatm,nparmwG
  INTEGER                                 :: nemmw,iemmw,nemmwG,iemmwG
  INTEGER, DIMENSION(:), allocatable      :: myPres
  !---BackgTuning Namelist items
  REAL :: backcldtopoc,backcldtopld,backcldthkoc,backcldthkld
  REAL :: backcldamtoc,backcldamtld,varcldtopoc,varcldtopld
  REAL :: varcldthkoc,varcldthkld,varcldamtoc,varcldamtld
  REAL :: backicetopoc,backicetopld,backicethkoc,backicethkld
  REAL :: backiceamtoc,backiceamtld,varicetopoc,varicetopld
  REAL :: varicethkoc,varicethkld,variceamtoc,variceamtld
  !---Climo Background Matrix (atm and sfc)
  REAL, DIMENSION(MxParG,mxbakAtm)        :: back_clim_atm_all
  REAL, DIMENSION(MxParG,mxbakSfc)        :: back_clim_sfc_all
  !---Climo Covariance Matrix (atm and sfc)
  REAL, DIMENSION(MxParG,MxParG,mxbakAtm) :: ut_mw_cov_clim_atm_all
  REAL, DIMENSION(MxParG,MxParG,mxbakSfc) :: ut_mw_cov_clim_sfc_all
  !---Climo EOFs Matrices (atm and sfc)
  REAL, DIMENSION(MxParG,MxParG,mxbakAtm) :: u_mw_clim_atm_all
  REAL, DIMENSION(MxParG,MxParG,mxbakSfc) :: u_mw_clim_sfc_all
  !---EOF Matrix / Selected Climo - Cascade - External
  REAL, DIMENSION(MxParG,MxParG) :: ut_mw_cov_clim,ut_mw_cov_casc,  &
                                    ut_mw_cov_ext,ut_mw_cov_emissdb 
  !---Background Matrix / Selected Climo - Cascade - External
  REAL, DIMENSION(MxParG)        :: xbakg_clim,xbakg_casc,xbakg_ext,&
                                    xbakg_emissdb
  real,    dimension(:),   allocatable :: frq,emissDb,emChgMetric
  INTEGER, DIMENSION(:),   allocatable :: pol,qcEmDb
  real,    dimension(:,:), allocatable :: covEmDb,cvt,tmp
  integer, dimension(:),   allocatable :: chanMap
  real,    dimension(:),   allocatable :: tbwts,tbEmDb
  INTEGER :: IEMISEXT_FIRST=1,ncidEmisdb,nbkg,irec,nchEmDb
  INTEGER :: IBKGREXT_FIRST=1,ncidBkgrExt,nprf,ifov,nparGext
  INTEGER :: IEMCHANGE_FIRST=1
  !---Indices
  TYPE(StateIndex_t)             :: IG,NG,IR,NR

CONTAINS

  SUBROUTINE loadBkgClim(nparmwGout,nparmwout,IGout,NGout,nlevel,vCoordTyp,pref,&
       nchmw,freq,polarity,IR_in,NR_in,nemmw_in,iemmwGout,molID,ih2o_in,nfor_in)
    USE bkg_io_module
    INTEGER, DIMENSION(:)             :: molID
    REAL,    DIMENSION(mxparg,mxparg) :: cov, eof
    REAL,    DIMENSION(mxparg)        :: bak
    REAL,    DIMENSION(mxlev)         :: pref
    REAL,    DIMENSION(mxchan)        :: freq
    INTEGER, DIMENSION(mxchan)        :: polarity
    INTEGER                           :: nLevel,itype,nchmw,nparmwGout
    CHARACTER(LEN=mxCoordTyp)         :: vCoordTyp
    INTEGER                           :: nparmwout,nemmw_in,ih2o_in,nfor_in,iemmwGout
    REAL                              :: sigPtop
    REAL, DIMENSION(:),   ALLOCATABLE :: vCoordA,vCoordB
    !---Indices
    TYPE(StateIndex_t)                :: IGAtm,NGAtm,IGout,NGout,IR_in,NR_in,IGSfc,NGSfc
    NAMELIST /bkgtuning/                                                  &
         backcldtopoc,backcldtopld,backcldthkoc,backcldthkld,             &
         backcldamtoc,backcldamtld,varcldtopoc,varcldtopld,               &
         varcldthkoc,varcldthkld,varcldamtoc,varcldamtld,                 &
         backicetopoc,backicetopld,backicethkoc,backicethkld,             &
         backiceamtoc,backiceamtld,varicetopoc,varicetopld,               &
         varicethkoc,varicethkld,variceamtoc,variceamtld

    !---Load retr bkgtuning parameters
    OPEN(U_algconfig,file=F_algconfig)
    READ(U_algconfig,bkgtuning)
    CLOSE(U_algconfig)
    !------------------------------------------------------------------------
    !  Atmospheric covariance (assumed to be the 1st in covarFile(1:2))
    !------------------------------------------------------------------------
    CALL queryCov(F_covar(1),nLevel=nLevel,vCoordTyp=vCoordTyp)  
    allocate(vCoordA(nlevel),vCoordB(nlevel),myPres(nlevel))
    SELECT CASE (TRIM(vCoordTyp))
    CASE (Pcoord)
       CALL queryCov(F_covar(1),pressure=pref(1:nLevel))
       myPres=pref(1:nLevel)
    CASE (Scoord)
       CALL queryCov(F_covar(1),sigCoord=vCoordA,sigPtop=sigPtop)
       CALL putSigmDefin(nlevel,sigPtop,vCoordA)
    CASE (Hcoord)
       CALL queryCov(F_covar(1),hybCoordA=vCoordA,hybCoordB=vCoordB)
       CALL putHybrDefin(nlevel,vCoordA,vCoordB)
    CASE DEFAULT
       print*,'err[BkgPrepModule::loadBkgClimo]: ',&
       'Unrecognized vertical coordinate type in cov file: ', &
       TRIM(vCoordTyp)
       CALL errorHalt(1)
    END SELECT
    NGAtm = initLengths()
    CALL openCov(ncid=U_cov, file=F_covar(1),nbkg=nAtmTypes,molid=molID, &
         NParG=NParGatm,IG=IGAtm,NG=NGAtm)  
    !---Consistency checks 
    IF (nlevel > mxlev) THEN
       PRINT *,'err[BkgPrepModule::LoadBkgClim]:  Nlev in Bkg File too large.'
       PRINT *,'nlevel (', nlevel, ') >  mxlev (', mxlev,')'
       call errorHalt(1)
    END IF
    DO itype = 1, nAtmTypes
       CALL getCov(ncid=U_cov,irec=itype,dmean=bak(1:NParGatm),             &
            un_atm=eof(1:NParGatm,1:NParGatm),cov=cov(1:NParGatm,1:NParGatm),&
            TypeFlag=AtmTypes(itype))
       IF (AtmTypes(itype)<0) THEN
          PRINT *,'err[BkgPrepModule::LoadBkgClim]:  AtmTypes undefined in Atm Covariance'
          PRINT *,'itype: ', itype, ' ... AtmTypes(itype): ', AtmTypes(itype)
          call errorHalt(1)
       END IF
       ut_mw_cov_clim_atm_all(1:npargatm,1:npargatm,itype)=cov(1:npargatm,1:npargatm)
       back_clim_atm_all(1:npargatm,itype)                =bak(1:npargatm)
       u_mw_clim_atm_all(1:npargatm,1:npargatm,itype)     =eof(1:npargatm,1:npargatm)
    ENDDO
    CALL closeCov(ncid=U_cov)
    !------------------------------------------------------------------------
    !      Surface covariance (assumed to be the 2nd in covarFile(1:2))
    !------------------------------------------------------------------------
    call queryCov(file=F_covar(2),nchmw=nchmw)
    CALL openCov(ncid=U_cov, file=F_covar(2),nbkg=nSfcTypes, &
         freq=freq(1:nchmw),polarity=polarity(1:nchmw))
    allocate(frq(nchmw),pol(nchmw))
    frq=freq(1:nchmw)
    pol=polarity(1:nchmw)
    !---Consistency checks 
    IF (nchmw > mxchan) THEN
       PRINT *,'err[BkgPrepModule::LoadBkgClim]:  nchMw in Bkg File too large.'
       PRINT *,'nchmw (', nchmw, ') > mxchan (', mxchan, ')'
       call errorHalt(1)
    END IF
    DO itype = 1, nSfcTypes
       CALL getCov(ncid=U_cov,irec=itype,dmEmMw=bak(1:nchmw),&
            un_emmw=eof(1:nchmw,1:nchmw),&
            covEmMw=cov(1:nchmw,1:nchmw),TypeFlag=SfcTypes(itype))
       IF (SfcTypes(itype)<0) THEN
          PRINT *,'err[BkgPrepModule::LoadBkgClim]:  SfcTypes undefined in Sfc Covariance'
          PRINT *,'itype: ', itype, ' ... SfcTypes(itype): ', SfcTypes(itype)
          call errorHalt(1)
       END IF
       ut_mw_cov_clim_sfc_all(1:nchmw,1:nchmw,itype)=cov(1:nchmw,1:nchmw)
       back_clim_sfc_all(1:nchmw,itype)             =bak(1:nchmw)
       u_mw_clim_sfc_all(1:nchmw,1:nchmw,itype)     =eof(1:nchmw,1:nchmw)
    ENDDO
    CALL closeCov(ncid=U_cov)
    !-----------------------------------------------------------------
    ! Merge Atmospheric and Surface Indices
    !-----------------------------------------------------------------
    NG = initLengths()
    NR = initLengths()
    NGout = initLengths()
    IR            = IR_in
    NR            = NR_in
    ih2o          = ih2o_in
    nemmw         = nemmw_in
    iemmw         = IR%cldLiq+NR%cldLiq
    NparMw        = getVectorLength(NR)+nemmw
    NG%temp       = NGAtm%temp
    NG%tskin      = NGAtm%tskin
    NG%psfc       = NGAtm%psfc
    NG%mol(ih2o)  = NGAtm%mol(ih2o)
    NG%cldLiq     = NGAtm%cldLiq
    NG%cldIce     = NGAtm%cldIce
    NG%wind       = NGAtm%wind
    IG            = genIndices(NG)
    nparmwG       = getVectorLength(NG)+nchmw
    iemmwG        = IG%wind+NG%wind
    nemmwG        = nchmw
    IGout         = IG
    NGout         = NG
    nparmwGout    = nparmwG
    nparmwout     = nparmw
    iemmwGout     = iemmwG
    nfor=nfor_in
    RETURN
  END SUBROUTINE loadBkgClim


  SUBROUTINE setFOVbkg(icasc,iembkgflg,iextbkg,iclassatm,igeo, &
       openCascInErr,psfc,ifor,vCoordTyp,press,umtx)
    !---Input variables
    INTEGER,                 INTENT(IN)    :: icasc,iembkgflg,iextbkg,iclassatm
    INTEGER,                 INTENT(IN)    :: openCascInErr,igeo,ifor
    REAL,                    INTENT(IN)    :: psfc
    CHARACTER(LEN=mxCoordTyp),INTENT(IN)   :: vCoordTyp
    REAL,    DIMENSION(:),   INTENT(INOUT) :: press
    REAL,    DIMENSION(:,:), INTENT(INOUT) :: umtx
    !---Local variables
    INTEGER                                :: iatmtype,isfctype,iatmIndx,isfcIndx

    !---Check consistency between iclassatm & AtmTypes
    iatmIndx=-1
    DO iatmtype=1,nAtmTypes
       if (iclassatm == AtmTypes(iatmtype)) then
          iatmIndx=iatmtype
          exit 
       endif
    ENDDO
    IF(iatmIndx == -1) THEN
       PRINT *,'err[BkgPrepModule::setFOVbkg]:  Atmtype not specified'
       PRINT *,'iclassatm (', iclassatm, ') not in AtmTypes (', AtmTypes(:nAtmTypes), ')'
       call errorHalt(1)
    END IF
    !---Check consistency between igeo & SfcTypes
    isfcIndx=-1
    DO isfctype=1,nSfcTypes
       if (igeo == SfcTypes(isfctype)) then
          isfcIndx=isfctype
          exit 
       endif
    ENDDO
    IF(isfcIndx == -1) THEN
       PRINT *,'err[BkgPrepModule::setFOVbkg]:  Sfctype not specified'
       PRINT *,'igeo (', igeo, ') not in SfcTypes (', SfcTypes(:nSfcTypes), ')'
       call errorHalt(1)
    END IF
    !---Select climatology information first
    CALL getBkgClim(iatmIndx,isfcIndx,umtx)
    !---Get scene pressure profile
    SELECT CASE (TRIM(vCoordTyp))
    CASE (Pcoord)
       press=myPres
    CASE (Scoord)
       press=getSigmPres(psfc)
    CASE (Hcoord)
       press=getHybrPres(psfc)
    END SELECT
    !---Load cascade-generated bckgrd information 
    IF(icasc    .GT.1) CALL getBkgCasc(icasc,openCascInErr,ifor,vCoordTyp,umtx)
    !---Load externally-generated bckgrd information 
    IF(iembkgflg.GT.0) CALL getBkgEmis(iembkgflg,igeo,ifor)
    IF(iextbkg  .GT.0) CALL getBkgExt(iextbkg,iclassatm)
    RETURN
  END SUBROUTINE setFOVbkg


  SUBROUTINE getBkgClim(iclassatm,igeo,umtx)
    !---Input variables
    INTEGER,                 INTENT(IN)    :: iclassatm,igeo
    REAL,    DIMENSION(:,:), INTENT(INOUT) ::  umtx
    !---Local variables
    INTEGER :: i
    REAL    :: umtx_clim(mxparg,mxparg)

    !----Background selection
    xbakg_clim(1:npargatm)        =back_clim_atm_all(1:npargatm,iclassatm)
    xbakg_clim(npargatm+1:nparmwg)=back_clim_sfc_all(1:nemmwG,igeo)
    !----Covariance matrix selection
    ut_mw_cov_clim=0
    ut_mw_cov_clim(1:nparGatm,1:nparGatm)                = &
         ut_mw_cov_clim_atm_all(1:nparGatm,1:nparGatm,iclassatm)
    ut_mw_cov_clim(nparGatm+1:nparmwg,nparGatm+1:nparmwg)= &
         ut_mw_cov_clim_sfc_all(1:nemmwG,1:nemmwG,igeo)
    !----Transformation matrix selection
    umtx_clim=0
    umtx_clim(1:NParGatm,1:NParGAtm)                = &
         u_mw_clim_atm_all(1:npargatm,1:npargatm,iclassatm)
    umtx_clim(NParGatm+1:NparmwG,NParGAtm+1:NParmwg)= &
         u_mw_clim_sfc_all(1:nemmwG,1:nemmwG,igeo)
    !----Truncate (compress) the transformation matrix
    extrvec(1:nparMw)= (/(i,i=IG%Temp,IG%Temp+NR%Temp-1), &
         (i,i=IG%Tskin,IG%Tskin+NR%Tskin-1), (i, i=IG%Psfc,IG%Psfc+NR%Psfc-1), &
         (i,i=IG%mol(iH2o),IG%mol(iH2o)+NR%mol(iH2o)-1),   &
         (i,i=IG%cldLiq,IG%cldLiq+NR%cldLiq-1),              &
         (i,i=IG%cldIce,IG%cldIce+NR%cldIce-1),              &
         (i,i=IEmMwG,IEmMwG+NEmMw-1)/)
    umtx(:,1:NParMw)=umtx_clim(:,extrvec(1:NParMw))
    RETURN
  END SUBROUTINE getBkgClim

  SUBROUTINE getBkgCasc(icasc,openCascInErr,ifor,vCoordTyp,umtx)
    !---Input variables
    INTEGER,                 INTENT(IN)    :: icasc,openCascInErr,ifor
    CHARACTER(LEN=mxCoordTyp),INTENT(IN)   :: vCoordTyp
    REAL,    DIMENSION(:,:), INTENT(INOUT) ::  umtx
    !---Local variables
    INTEGER :: rdCascInErr,rdCascInErr1,rdCascInErr2,rdCascInErr3,rdCascInErr4
    INTEGER :: iread,i,j
    REAL    :: chisq_casc
    REAL, DIMENSION(Mxpar,Mxpar)     :: mtrx,SxMw_Dum
    REAL, DIMENSION(MxparG)          :: vector,xBakG_dum
    REAL, DIMENSION(NParMwG,NParMw)  :: tmplt

    xbakg_casc     = 0.
    ut_mw_cov_casc = 0.
    !---Read in Cascade data set from previous low-resolution retrieval
    IF (icasc.GT.1 .AND. openCascInErr .EQ. 0) THEN
       READ(U_cascIn,IOSTAT=rdCascInErr) chisq_casc
       IF (chisq_casc .LT. 0.2 .AND. rdCascInErr .EQ. 0) THEN
          iread=1
          READ(U_cascIn,IOSTAT=rdCascInErr1) vector(1:nParMwG)
          READ(U_cascIn,IOSTAT=rdCascInErr2) mtrx(1:nparmw,1:nparmw)
          IF (rdCascInErr1 .EQ. 0 .AND. rdCascInErr2 .EQ. 0) THEN
             xbakg_casc(1:nParMwG) = vector(1:nParMwG)
             !--Transform covariance matrix to geophysical space (subject to trades)
             tmplt=MATMUL(umtx(1:NParMwG,1:NParMw),mtrx(1:NParMw,1:NParMw))
             ut_mw_cov_casc(1:NParMwG,1:NParMwG)= &
                  MATMUL(tmplt,TRANSPOSE(umtx(1:NParMwG,1:NParMw)))
             DO i=IG%mol(iH2o),IG%mol(iH2o)+NG%mol(ih2o)-1
                xbakg_casc(i)=alog(xbakg_casc(i))
             END DO
             xbakg_casc(IG%cldLiq)=amin1(xbakg_casc(IG%cldLiq),900.0)
             !--set Tskin to delta Tskin scheme:
             IF (TRIM(vCoordTyp) == Pcoord) &
                xbakg_casc(IG%Tskin)= xbakg_casc(IG%Tskin)-xbakg_casc(NG%temp)
             xbakg_casc(IG%cldLiq+2)=0.0
             DO i=IG%mol(ih2o),IG%mol(ih2o)+NG%mol(ih2o)-1
                DO j=1,nparmwg
                   IF(i.NE.j)THEN
                      ut_mw_cov_casc(i,j)=0.0
                      ut_mw_cov_casc(j,i)=0.0
                   END IF
                END DO
                !-- Relax H2O cov:
                DO j=IG%mol(ih2o),IG%mol(ih2o)+NG%mol(ih2o)-1
                   ut_mw_cov_casc(j,i)=ut_mw_cov_casc(j,i)*1.7
                END DO
             END DO
             !---remove correlation between cloud amount and other parameters
             DO i=IG%cldLiq,IG%cldLiq+2 
                DO j=1,nparmwg
                   IF(i.NE.j)THEN
                      ut_mw_cov_casc(i,j)=0.0
                      ut_mw_cov_casc(j,i)=0.0
                   END IF
                END DO
             END DO
             ut_mw_cov_casc(IG%cldLiq+2,IG%cldLiq+2)=0.04
          ENDIF
       ELSE                ! just read over retrieved profile
          READ(U_cascIn,IOSTAT=rdCascInErr3) xBakG_dum(1:nParMwG)
          READ(U_cascIn,IOSTAT=rdCascInErr4) sxmw_dum(1:nparmw,1:nparmw)
          iread=0
       END IF
    END IF
    RETURN
  END SUBROUTINE getBkgCasc


  SUBROUTINE getBkgEmis(iembkgflg,igeo,ifor)
    USE bkg_io_module
    use ncdf_module
    use SpectralOperations
    !---Input variables
    INTEGER, INTENT(IN) :: iembkgflg,igeo,ifor
    !---Local variables
    integer                            :: nqc,nMtr,i
    real,    parameter                 :: dfrqLim=0.001 ! Freq match threshold
                                                        ! (fractional) 
    real,    dimension(:), allocatable :: frqdb
    integer, dimension(:), allocatable :: poldb
    
    if(IEMISEXT_FIRST == 1)then
       call queryCov(file=F_ldEmissExt,nchmw=nchEmDb)
       allocate(frqdb(nchEmDb),poldb(nchEmDb))
       call openCov(ncid=ncidemisdb,file=F_ldEmissExt,nqc=nqc,nbkg=nbkg, &
            freq=frqdb,polarity=poldb,status='old')
       if(nchEmDb /= nemmwG)then
          print *,'msg[BkgPrepModule::getBkgEmis]:  number of channels different ',&
               'in driver and EmisDb file - will convert to channel frequencies'
          print *,'nchEmDb (', nchEmDb,') /= nemmwG (',nemmwG,')'
          allocate(cvt(nemmwG,nchEmDb),tmp(nemmwG,nchEmDb)) !--Arrays for frq conversion
          call frqConvertGeneral(frqdb,poldb,frq,pol,cvt)
       else
          if (any((abs(frqdb-frq)/frqdb) > dfrqLim) .or. &
               any(poldb /= pol)) then
             print *,'err[BkgPrepModule::getBkgEmis]:  Channel set mismatch;'
             do i=1,nemmwG
                write(*,'(i4,2(f10.4,i4))')i,frqdb(i),frq(i),poldb(i),pol(i)
             enddo
             call errorHalt(1)    
          endif
       endif
       nMtr=readNcdfDim(ncidemisdb,'nMtr',silent=.true.)
       allocate(emissdb(nchEmDb),covEmDb(nchEmDb,nchEmDb),emChgMetric(nMtr),qcEmDb(nqc))
       IEMISEXT_FIRST=0
    endif
    
    if(nbkg /= nfor)then
       print *,'err[BkgPrepModule::getBkgEmis]:  numbers of FOVs ',&
            'in EmisDb and rad files do not match;'
       print *,'nbkg (',nbkg,') /= nfor (',nfor,')'
       call errorHalt(1)
    endif
    call getCov(ncid=ncidemisdb,irec=ifor,dmEmMw=emissDb,covEmMw=covEmDb,qc=qcEmDb)
    call readNcdfData(ncidemisdb,emChgMetric,varname='metric',record_no=ifor, &
         status='single_record')
    if(all(qcEmDb == 0))then
       if(nchEmDb /= nemmwG)then
          xbakg_emissdb(1:nemmwG)=matmul(cvt,emissdb)
          tmp = matmul(cvt,covEmDb)
          ut_mw_cov_emissdb(1:nemmwG,1:nemmwG) = matmul(tmp,transpose(cvt))
       else
          xbakg_emissdb(1:nemmwG)=emissdb
          ut_mw_cov_emissdb(1:nemmwG,1:nemmwG) = covEmDb
       endif
    endif
    
    if(ifor == nbkg)then
       call closeCov(ncid=ncidEmisdb)
       if(allocated(frqdb))deallocate(frqdb,poldb)
    endif
  
    RETURN
  END SUBROUTINE getBkgEmis


  SUBROUTINE getBkgExt(iextbkg,igeo)
    USE scene_io_module
    !---Input variables
    INTEGER, INTENT(IN) :: iextbkg,igeo
    !---Local variables
    INTEGER, dimension(maxMol) :: molIdExt
    !---Indices
    TYPE(StateIndex_t)         :: IGExt,NGExt
    
    NGExt = initLengths()
    if(IBKGREXT_FIRST == 1)then
       molIDExt=nullMolId
       call openScene(ncid=ncidBkgrext,file=F_bkgrdExtern,nprf=nprf, &
            molID=molIDExt,IG=IGExt,NG=NGExt,status='old')
       nparGExt=getVectorLength(NGExt)
       IBKGREXT_FIRST=0
       ifov=1
    endif
    !--Gross consistency check. Should be refined 
    !--to make sure that the structures are identical
    if(nparGExt /= getVectorLength(NG))then
       print *,'err[BkgPrepModule::getBkgExt]:  number of elements in ',&
            'NGs different in driver and BkgrdExtern file'
       print *,'nparGExt (',nparGExt,') /= getVectorLength(NG) (',getVectorLength(NG),')'
       call errorHalt(1)
    endif
    if(ifov > nprf)then
       print *,'err[BkgPrepModule::getBkgExt]:  number of profiles in ',&
            'BkgrExtern file exceeded'
       print *,"ifov (", ifov, ") > nprf (", nprf, ")"
       call errorHalt(1)
    endif
    call getScene(ncid=ncidBkgrext,irec=ifov,x=xbakg_ext(1:nparGExt))
    ifov=ifov+1
    if(ifov == nprf + 1)then
       call closeScene(ncid=ncidBkgrext)
    endif
    ut_mw_cov_ext =0.
    RETURN
  END SUBROUTINE getBkgExt

  SUBROUTINE emisChgTest(tbin,chanIDfull,useEmisDb)
    !---Input variables:
    real,     dimension(:),         intent(in)    :: tbin
    character(len=12),dimension(:), intent(in)    :: chanIDfull
    integer,                        intent(inout) :: useEmisDb
    
    !---Local variables:
    integer                            :: i,j
    real                               :: emChgThresh,tbmetric,metComp
    integer                            :: nchEmTest
    character(len=12)                  :: chanID_loc1,chanID_loc2
    character(len=12),dimension(:),allocatable :: chanIDwt
    
    if(any(qcEmDb /= 0)) then
       useEmisDb=0
       return
    endif
    
    if(IEMCHANGE_FIRST == 1)then
       !--- Open tuning file:
       open(U_emChange,file=F_emChange)
       read(U_emChange,'(i3)') nchEmTest
       allocate(tbwts(nchEmTest),chanIDwt(nchEmTest),tbEmDb(nchEmTest),chanMap(nchEmTest))
       !---Read change detection threshold, weighting coefficient data and
       !     channel subset definition:
       read(U_emChange,'(f7.4)') emChgThresh
       read(U_emChange,'(a12,1x,f11.8)')(chanIDwt(i),tbwts(i),i=1,nchEmTest)
       close(U_emChange)
       !---Compare channel subset with full channel set from driver.  The channel IDs
       !     would need to be present in radfile and need to be read and passed into
       !     this subroutine.  Have a check to make sure driver channel set and tuning
       !     file channel set are consistent.
       !---Determine indices of subset channels, and save into an integer vector.
       do i=1,nchEmTest
          do j=1,nemmwG
             chanID_loc1 = trim(adjustl(chanIDwt(i)))
             chanID_loc2 = trim(adjustl(chanIDfull(j)))
             if(chanID_loc1(1:12) == chanID_loc2(1:12))then
                chanMap(i)=j
                exit
             endif
          enddo
          if(j == (nemmwG+1))then
             print*,'err[BkgPrepModule::emisChgTest]:  channel in weighting file ',&
                  'not present in input radfile: '
             print*,chanIDwt(i)
             call errorHalt(1)
          endif
       enddo
       IEMCHANGE_FIRST=0
    endif
    !---Subset the input set of Tb:
    tbEmDb     = tbin(chanMap(1:nchEmTest))
    !---Compute emissivity change metric:
    tbmetric = dot_product(tbwts,tbEmDb)
    !---Compare metric to metric from dynamic database, and set flag accordingly:
    !    This assumes 'metric' is a real vector of length one.
    metComp  = tbmetric - emChgMetric(1)
    if(abs(metComp) <= emChgThresh)then
       useEmisDb=1
    else
       useEmisDb=0
    endif
    return
  END SUBROUTINE emisChgTest
  
  SUBROUTINE combBkg(icasc,iembkgflg,useEmisDb,iextbkg,dback0,ut_mw_cov_net,psfc)
    !---Input variables
    INTEGER :: icasc,iembkgflg,useEmisDb,iextbkg
    REAL    :: psfc
    !---Output variables
    REAL, DIMENSION(:),   INTENT(INOUT) :: dback0
    REAL, DIMENSION(:,:), INTENT(INOUT) :: ut_mw_cov_net
    !---Local variables
    INTEGER :: i,k,j,itype,ii,jj

    !---By default load the climatology background into the worskpace array
    dback0(1:NParmwG)=xbakg_clim(1:NParmwG)
    ut_mw_cov_net(1:NParMwG,1:NParMwG)=ut_mw_cov_clim(1:NParMwG,1:NParMwG)
    IF(iembkgflg > 0) THEN
       IF(useEmisDb > 0) THEN           !load local emissivity
          dback0(iemmwg:iemmwg+nemmwg-1)=xbakg_emissdb(1:nemmwG)
          ut_mw_cov_net(nParGatm+1:nparMwG,nParGatm+1:nparMwG) = &
               ut_mw_cov_emissdb(1:nemmwG,1:nemmwG)
       ENDIF
    END IF
    IF(icasc > 1) THEN
       PRINT*,'combine casc with climo'
    ELSE
       !print*,'load climo atm data'
    END IF
    IF(iextbkg > 0) THEN
       PRINT*,'combine climo with some external data in FOV workspace'
    END IF
    dback0(IG%Psfc)=psfc
    RETURN
  END SUBROUTINE combBkg


  SUBROUTINE tuneBkg(dback0,ut_mw_cov_net,iclassatm,igeo,clw_cov_mw,TightCLDcov,mwCld)
    !---Input variables
    REAL, DIMENSION(:,:), INTENT(INOUT) :: ut_mw_cov_net
    REAL, DIMENSION(:),   INTENT(INOUT) :: dback0
    REAL, DIMENSION(:,:), INTENT(OUT)   :: clw_cov_mw
    REAL                                :: psfc,TightCLDcov
    INTEGER,              INTENT(IN)    :: iclassatm,igeo,mwCld

    IF(igeo.LE.1)THEN                      !---Ocean
       dback0(IG%cldLiq)                       = backcldtopoc
       dback0(IG%cldLiq+1)                     = backcldthkoc
       dback0(IG%cldLiq+2)                     = backcldamtoc
       ut_mw_cov_net(IG%cldLiq,IG%cldLiq)      = varcldtopoc
       ut_mw_cov_net(IG%cldLiq+1,IG%cldLiq+1)  = varcldthkoc
       ut_mw_cov_net(IG%cldLiq+2,IG%cldLiq+2)  = varcldamtoc
    ELSE                                   !---Land
       dback0(IG%cldLiq)                       = backcldtopld
       dback0(IG%cldLiq+1)                     = backcldthkld
       dback0(IG%cldLiq+2)                     = backcldamtld
       ut_mw_cov_net(IG%cldLiq,IG%cldLiq)      = varcldtopld
       ut_mw_cov_net(IG%cldLiq+1,IG%cldLiq+1)  = varcldthkld
       ut_mw_cov_net(IG%cldLiq+2,IG%cldLiq+2)  = varcldamtld
    ENDIF

    !---In case the cloud is not retrieved, tighten the covariance
    IF(mwCld.EQ.0)THEN
       ut_mw_cov_net(IG%cldLiq:IG%cldLiq+NG%cldLiq-1,IG%cldLiq:IG%cldLiq+NG%cldLiq-1)=&
            TightCLDcov
    END IF
    !---Set the cloud covariance matrix
    clw_cov_mw(1:NR%cldLiq,1:NR%cldLiq) = &
         ut_mw_cov_net(IG%cldLiq:IG%cldLiq+NG%cldLiq-1,IG%cldLiq:IG%cldLiq+NG%cldLiq-1) 
    RETURN
  END SUBROUTINE tuneBkg


  SUBROUTINE transBkg(ut_mw_cov_net,ut_mw_cov,umtx)
    !---Input variables
    REAL, DIMENSION(:,:), INTENT(IN)    :: ut_mw_cov_net
    REAL, DIMENSION(:,:), INTENT(INOUT) :: umtx
    !---Output variables
    REAL, DIMENSION(:,:), INTENT(OUT)   :: ut_mw_cov
    !---Local variables
    INTEGER                             :: i
    REAL, DIMENSION(NParMw,NParMwG)     :: tmpl
    REAL, DIMENSION(NParMw,NParMw)      :: sxmwWnumerErr

    tmpl=MATMUL(TRANSPOSE(umtx(1:NParMwG,1:NParMw)), &
         ut_mw_cov_net(1:NParMwG,1:NParMwG))
    ut_mw_cov(1:NParMw,1:NParMw)=MATMUL(tmpl,umtx(1:NParMwG,1:NParMw))
    RETURN
  END SUBROUTINE transBkg

  SUBROUTINE closeBkg()
    if(allocated(frq))deallocate(frq,pol)
    if(allocated(cvt))deallocate(cvt,tmp)
    if(allocated(emissdb))deallocate(emissdb,covEmDb,emChgMetric,qcEmDb)
    if(allocated(tbwts))deallocate(tbwts,tbEmDb,chanMap)
    RETURN
  END SUBROUTINE closeBkg

END MODULE BkgPrepModule
