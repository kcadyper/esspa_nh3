!--------------------------------------------------------------
!
!  MODULE OSS_IR_MODULE: contains items needed for the
!                OSS Forward Model (IR). AER Inc. 2004.
!
! here data common for oss_ir_module and oss_ir_module_scat is stored
!--------------------------------------------------------------

MODULE OSSIRmoduleSubs
  USE asolar_source_function
  use params_module, ONLY: mxlev, mxlay, mxhmol, maxcmu, MxmolS, mxParG
  implicit none
  PRIVATE
  !------------------------------------------------------------------------
  !     Public items made available to calling program(s)
  !------------------------------------------------------------------------
  PUBLIC :: getOD
  PUBLIC :: OSSdestroyIRsubs
  PUBLIC :: addChanSelect, setChanSelect
  PUBLIC :: addChanSelectMask
  PUBLIC :: getMolecMass
  PUBLIC :: getChannelFrequency, getCountUsedNode
  PUBLIC :: getChannelIndex
  PUBLIC :: getCountChannel
  PUBLIC :: getUsedNodes, getLastSelectionIndex
  PUBLIC :: setPlaneParallelFlag
  PUBLIC :: getReferencePressure
  PUBLIC :: getSelectionIndex
  PUBLIC :: getConstant
  PUBLIC :: findFreeUnit
  PUBLIC :: planck
  PUBLIC :: rlin
  public :: fpathZ
  public :: fpathP
  public :: layerAverage
  public :: setpath_ir
  public :: settabindx_ir
  public :: osstran
  public :: vinterp
  public :: setIndex_selfCont
  public :: getCurrentSelection

  character(len=*), parameter :: modName = 'OSSmoduleSubsIR'

  !------------------------------------------------------------------------
  ! Parameters/Dimensions
  !------------------------------------------------------------------------
  real, public,   PARAMETER  :: pi=acos(-1.0), deg2rad = pi/180.0

  INTEGER, PARAMETER :: mxIndex=10
  real,    PARAMETER :: dbAvg2dtau=-1./12.0 !Weighting factor derivative wrt optical depth, for linear-in-tau, in low-OD limit
  real,    PARAMETER :: grav_const_req1= 9.80665 ! m / s^2
  real,    PARAMETER :: grav_const_req2=-0.02586 ! m / s^2
  real,    PARAMETER :: Rgas=8.3144621   ! Universal gas constant (J / mol / K)
  real,    PARAMETER :: drymwt=28.964 ! Dry air molecular mass (g / mol)
  real,    PARAMETER :: AVOGAD = 6.02214199E+23 ! Avogadro ( 1 / mol )
  real,   parameter  :: BOLTZ = 1.3806503E-16
  real,    parameter :: GASCON = 8.314472E+07

  INTEGER, PARAMETER :: MAXSPC=84
  real, DIMENSION(MAXSPC), PARAMETER  :: molWt=(/ &
                18.015, 44.010, 47.998,44.010, 28.011, 16.043, 31.999, &   ! 7
                30.010, 64.060, 46.010,17.030, 63.010, 17.000, 20.010, &   ! 14
                36.460, 80.920,127.910,51.450, 60.080, 30.030, 52.460, &   ! 21
                28.014, 27.030, 50.490,34.010, 26.030, 30.070, 34.000, &   ! 28
                66.010,146.050, 34.080,46.030,  0.   ,  0.   ,  0.   , &   ! 35
                 0.   ,  0.   , 28.053,32.042,  0.   ,  0.   ,  0.   , &   ! 42
                 0.   ,  0.   , 0.   ,0.   ,  0.   ,  0.   ,  0.   , &   ! 49
                 0.   ,153.820, 88.000,97.460,137.37,187.380,120.910 , &   ! 56
                 0.   , 86.470,108.010,  0.  , 68.12,121.05  ,  0.   , &   ! 63
                 0.   ,  0.   , 0.   ,0.   ,  0.   ,  0.   ,  0.   , &   ! 70
                 0.   ,  0.   , 0.   ,0.   ,  0.   ,  0.   ,  0.   , &   ! 77
                 0.   ,  0.   , 0.   ,19.02 ,  0.   ,  0.   ,  0.   /)   ! 84
  !  molecular weight normalized on dry air weight
   real, dimension(MAXSPC), parameter, public  :: normMolWt= molWt/drymwt

  real,    PARAMETER :: mperkm=1000. ! (m/km)
  real,    PARAMETER :: rairDflt=10.0/grav_const_req1 !10./g ->converts p(mb) in g/cm**2
  real,    PARAMETER :: rConst =  Rgas/(drymwt*1.E-3)/grav_const_req1/mperkm
  real,    PARAMETER :: rEarth = 6371.23
  real,    PARAMETER :: pMatchLim = 1.E-3

  !------------------------------------------------------------------------
  !     Lookup tables
  !------------------------------------------------------------------------
  INTEGER, PARAMETER :: expectedEndianCode = 123456789
  INTEGER, PARAMETER :: LUT_KIND = KIND(1)
  REAL(KIND = LUT_KIND), allocatable :: TmpSelf(:)

  INTEGER                                  :: nLayOD,nTmpOD,nWvpOD,nTself

  !------------------------------------------------------------------------
  !     OSS control and status flags
  !------------------------------------------------------------------------
  ! if compartibility with LBLRTM is required set to .false.
  LOGICAL, public, PARAMETER               :: OSS_LIN_IN_TAU = .FALSE.
  LOGICAL, save, public                    :: myLinInTau =.FALSE.
  logical, save, public                    :: sphericalGeometryFlag = .false.
  logical, save, public                    :: zIntegrationFlag      = .false.

  !------------------------------------------------------------------------
  ! loadDone implementation provides the option to call getOD() in advance
  LOGICAL, save                            :: getODdone  =.false.

  !------------------------------------------------------------------------
  ! Pressure Grid / Molecules / Geophysical Pointers
  !------------------------------------------------------------------------
  INTEGER                                               :: nNodes
  INTEGER                                               :: nChan
  INTEGER                                               :: numRefLev
  INTEGER                                               :: MolID(MxHmol)
  INTEGER,                                public        :: nMol
  INTEGER,                                public        :: nfsmp_ir
  INTEGER                                               :: nf_sel
  INTEGER,                                public        :: nChMax
  INTEGER, dimension (:),    allocatable, public        :: iselS
  INTEGER, dimension (:),    allocatable, public        :: NmolS
  INTEGER, dimension (:,:),  allocatable, public        :: ImolS
  real,    dimension (:,:),  allocatable, public        :: kfix_ir
  real,    dimension (:,:),  allocatable                :: dkFix_ir
  real,    dimension (:,:),  allocatable, public        :: kvar_ir
  real,    dimension (:,:),  allocatable, public        :: dkvar_ir
  real,    dimension (:,:),  allocatable, public        :: kh2o_ir
  real,    dimension (:,:),  allocatable, public        :: dkh2o_ir
  real,    dimension (:,:),  allocatable, public        :: kself
  real,    dimension (:,:),  allocatable, public        :: khdo_ir
  real,    dimension (:,:),  allocatable, public        :: dkhdo_ir
  real,    dimension (:),    allocatable                :: fixDMR
  real,    dimension (:,:),  allocatable                :: Tmptab
  real,    dimension (:,:),  allocatable                :: Wvptab
  logical                                               :: isHDOProfile

  !------------------------------------------------------------------------
  ! Planck function parameters
  real,    dimension (:),    allocatable, public        :: f1_arr
  real,    dimension (:),    allocatable, public        :: f2_arr

  !------------------------------------------------------------------------
  ! solar constant
  real,    dimension (:),    allocatable, public        :: sunrad

  !------------------------------------------------------------------------
  ! OSS selection file
  INTEGER*2, dimension(:), allocatable                  :: nch_Loc
  INTEGER,   dimension(:), allocatable, public          :: nch_ir

  !------------------------------------------------------------------------
  ! channel selection
  INTEGER, PARAMETER, public                            :: UNDEFINED_INDEX = -999
  integer                                               :: nchanAll
  INTEGER                                               :: iChSet=UNDEFINED_INDEX
  INTEGER                                               :: lastChSet ! Number of loaded chan sets
  integer, dimension(0:mxIndex)                         :: userIndex2chSetIndexMap

  INTEGER, dimension(:,:), allocatable, public          :: ichmap_ir
  real,    dimension(:,:), allocatable, public          :: coef_ir

  INTEGER, dimension(:),     allocatable                :: nChList_arr
  INTEGER, dimension(:,:),   allocatable                :: chanList_arr
  INTEGER, dimension(:),     allocatable                :: nNodes_arr
  INTEGER, dimension(:),     allocatable                :: isel_Loc
  INTEGER, dimension(:,:),   allocatable                :: iselS_arr
  INTEGER, dimension(:,:),   allocatable                :: nch_arr
  INTEGER, dimension(:,:),   allocatable                :: ichmap_Loc
  INTEGER, dimension(:,:,:), allocatable                :: ichMap_arr
  real,    dimension(:,:,:), allocatable                :: coef_arr


  !------------------------------------------------------------------------
  ! meteorology
  !------------------------------------------------------------------------
  real, dimension(:),           allocatable             :: pRef
  real, dimension(:),           allocatable             :: pavlref

  !------------------------------------------------------------------------
  ! Instrument parameters
  !------------------------------------------------------------------------
  REAL*8,  dimension(:),   allocatable, public          :: vWvn
  real,    dimension(:),   allocatable                  :: cWvn

  REAL(KIND=LUT_KIND),dimension(:,:), allocatable       :: coef_Loc


contains
  !----------------------------------------------------------------------------
  !    PURPOSE: Reads the pre-computed Look Up Tables needed later in the
  !    OSS radiative transfer model.
  !  selfile instrument chanel file
  !  odfile  optical properties LUT file
  !  defProfFile standard profile file
  !  nmol_in number of variable molecules
  !  molID_in variable gas HITRAN indices
  !----------------------------------------------------------------------------
  SUBROUTINE GetOD(selfile,odfile,defProfFile,nmol_in,molID_in)
    REAL , PARAMETER          :: SOfL =  29.9792458 !speed of light in km/sec, scaled by 1e-4

    !LBLRTM version
!   RADCN1 = 2.*PLANCK*CLIGHT*CLIGHT*1.E-07
!   RADCN2 = PLANCK*CLIGHT/BOLTZ
    real, parameter  :: RADCN1 = 1.191042722E-12
    real, parameter  :: RADCN2 = 1.4387752

    ! stands in numerator in dependence of radiance on TB,
    ! units: (mW/(m^2 ster cm^-1)) / [cm^-1)^3, scaled by 1e9
    REAL , PARAMETER          :: f1 = RADCN1 * 1e7
    !stands in exponential dependence of radiance on TB, like f2*wn*1e-2/TB
    REAL , PARAMETER          :: f2 = RADCN2

    !---Input variables
    CHARACTER(len=*),      INTENT(IN) :: selfile,odfile,defProfFile
    INTEGER,               INTENT(IN) :: nmol_in
    INTEGER, DIMENSION(:), INTENT(IN) :: molID_in
    !---Local variables
    INTEGER                   :: ius,iuo,iuF
    CHARACTER(len=100)        :: instr_info
    CHARACTER(len=24)         :: date_IDSel,date_IDLUT
    CHARACTER(len=20)         :: Sensor_ID,hdrOpt,Unit_Char_s,Unit_Char_l
    CHARACTER*12, allocatable :: chID(:)
    INTEGER                   :: fmtver,nHdrOpt
    INTEGER                   :: fmtLUT,osstran_opt
    INTEGER                   :: WMO_Satellite_ID,WMO_Sensor_ID,Sensor_Type
    INTEGER                   :: indxChSel
    INTEGER                   :: NmolS_tmp,ImolS_tmp(MxmolS),imol
    INTEGER                   :: NmolS_,ImolS_(MxmolS)
    INTEGER                   :: imols_indx(MxmolS)
    INTEGER                   :: endianCode,k,m,kk,ismp
    INTEGER                   :: i,nx,ks
    INTEGER                   :: nsize1,nsize2,nsize3,nsize4,nsize5
    INTEGER                   :: ihO,Spc_Units_s,Spc_Units_l
    INTEGER, allocatable      :: chIndx(:), iPol(:)
    REAL                      :: odfac = 0.
    REAL   , allocatable      :: kvar_tmp(:)
    INTEGER, allocatable      :: xid(:),xid_(:),xid_dpnd(:)
    REAL   , allocatable      :: vmol(:),qr(:)
    REAL   , allocatable      :: wvpTmp(:,:)
    INTEGER, DIMENSION(MxHmol):: map,mapS,iflag = 1
    INTEGER, allocatable      :: iSelFix(:),mFix(:),invMolID(:)
    INTEGER, allocatable      :: chanListTmp(:)
    INTEGER                   :: nwvpod_
    INTEGER                   :: nselFix,nFix,nLevF,icnt
    REAL(KIND=LUT_KIND)       :: V1m,V2m
    REAL(KIND=LUT_KIND),allocatable :: cWvnLoc(:),pRefLoc(:),fixLoc(:)
    REAL(KIND=LUT_KIND),allocatable :: tmpLoc(:,:),wvpLoc(:,:),dfltLoc(:,:)
    REAL(KIND=LUT_KIND),allocatable :: kfix_Loc(:),kh2o_Loc(:),kvar_Loc(:)
    REAL(KIND=LUT_KIND),allocatable :: dkfix_Loc(:),dkh2o_Loc(:)
    REAL(KIND=LUT_KIND),allocatable :: kself_Loc(:)

    real, dimension (:,:) , allocatable :: dkFix_ir
    real, dimension (:)   , allocatable :: fixDMR(:)
    real, dimension (:,:) , allocatable :: Wvptab
    INTEGER                             :: MolIDFix(MxHmol)
    INTEGER                             :: NmolFix

    real                                :: scale, dummyReal
    logical                             :: isHDO, isExist
    character(len=*), parameter           :: procName=modName//'GetOD'
    character(len=*), parameter :: errHeader='Err['//procName//']: '
    character(len=*), parameter :: wrnHeader='Wrn['//procName//']: '

    if(getODdone) return
    userIndex2chSetIndexMap = UNDEFINED_INDEX

    !---Check input molecular selection
    DO k=2,nmol_in
       IF (molID_in(k)<=molID_in(k-1)) THEN
          print*, errHeader,'Variable molecule IDs  in varMolID must be'&
                ,' provided in ascending order'
          call exit(1)
       END IF
    END DO
    IF (molID_in(1).NE.1) THEN
        print*, errHeader,' Water vapor must be included in list of'&
              ,' variable molecules'
        call exit(1)
    END IF
    !---open files
    ius = findFreeUnit()
    OPEN(ius,file=selfile,form='unformatted',status='old',convert='little_endian')
    !---Read selFile header
    READ(ius)endianCode  ! endian verification
    IF (endianCode /= expectedEndianCode) THEN
       close(ius)
       print *,'attampting to read SEL as little_endian'
       OPEN(ius,file=selfile,form='unformatted',status='old',convert='big_endian')
       READ(ius)endianCode  ! endian verification
       IF (endianCode /= expectedEndianCode) THEN
          print*, errHeader ,' Wrong Format of binary file: Sel', &
                endianCode
          call exit(1)
       END IF
    END IF
    READ(ius)date_IDsel,fmtver !date: YYYYMMDDHHMMSS.sss-0X000, X=4,or,5
    READ(ius)WMO_Satellite_ID,WMO_Sensor_ID !if specifically undefined: -1 -1
    READ(ius)Sensor_Type,Sensor_ID !Sensor_Type: 1 for IR
    READ(ius)instr_info ! any useful information
    READ(ius)nChan,nf_sel,nChMax
    allocate(chIndx(nChan),chID(nChan),iPol(nChan),cWvn(nChan),cWvnLoc(nChan))
    READ(ius)hdrOpt
    hdrOpt  = adjustl(hdrOpt)
    nHdrOpt = len_trim(hdrOpt)
    chIndx = (/(i,i=1,nChan)/)

    chID(1:nChan)='NA'
    Spc_Units_s= 0      !default due to IR
    Unit_Char_s= 'cm-1' !default due to IR
    iPol(1:nChan)=-1    !decide about default (?)
    DO ihO = 1, nHdrOpt
       SELECT CASE (HdrOpt(ihO:ihO))
          CASE ('K')
             READ(ius) chIndx(1:nChan)
          CASE ('I')
             READ(ius) chID(1:nChan)
          CASE ('U')   ! for spectral units
             READ(ius) Spc_Units_s, Unit_char_s
          CASE ('C')
             READ(ius)cWvnLoc(1:nChan)
             !From GHz to cm^-1
             IF (Spc_Units_s == 1) cWvnLoc(1:nChan)= cWvnLoc(1:nChan)/SOfL
             !From mu  to cm^-1
             IF (Spc_Units_s == 2) cWvnLoc(1:nChan)= 1E4/cWvnLoc(1:nChan)
             cWvn(:)=cWvnLoc(:)
          CASE ('P')
             READ(ius)iPol(1:nChan)
          CASE DEFAULT
             PRINT *,wrnHeader, 'Selection file header content flag not'&
                    ,' recognized:',HdrOpt(iHO:iHO)
       END SELECT
    END DO
    deallocate (cWvnLoc)
    !      LUT file
    iuo = findFreeUnit()
    OPEN(iuo,file=odfile,form='unformatted',status='old',convert='little_endian')
    READ(iuo)endianCode ! endian verification
    IF (endianCode /= expectedEndianCode) THEN
       close(iuo)
       print *,'attampting to read LUT as big_endian'
       OPEN(iuo,file=odfile,form='unformatted',status='old',convert='big_endian')
       READ(iuo)endianCode  ! endian verification
       IF (endianCode /= expectedEndianCode) THEN
          print*, 'Err[oss_ir_module::GetOD]: wrong Format of binary file: Sel', &
                endianCode
       call exit(1)
      END IF
    END IF

    READ(iuo)date_IDlut,fmtLUT !date: YYYYMMDDHHMMSS.sss-0X000 (X: 4, or, 5)
    IF (date_IDlut /= date_IDsel) THEN
       print*,errHeader,' IDs for sel and lut are inconsistent'
       call exit(1)
    END IF
    READ(iuo)HdrOpt
    hdrOpt  = adjustl(hdrOpt)
    nHdrOpt = len_trim(hdrOpt)
    !Defaults
    nx = 1 !For dK, number of dependent, e.g., dkFix, dkh2o, dkvar(:). Defualt
           !requires dkh2o only
    ALLOCATE (xid(nx), xid_dpnd(nx))
    xid(1)      = 1
    xid_dpnd(1) = 1
    ossTran_Opt = 1
    Spc_Units_l = 0
    Unit_Char_l = 'cm-1'
    DO iHO = 1, nHdrOpt
       SELECT CASE (HdrOpt(iHO:iHO))
          CASE ('F')
             READ(iuo)ossTran_Opt ! 1 is for IR
          CASE ('U')   ! for spectral units
             READ(iuo) Spc_Units_l, Unit_char_l
          CASE ('X')
             READ(iuo)nx
             IF ( ALLOCATED(xid) ) DEALLOCATE(xid)
             IF ( ALLOCATED(xid_dpnd) ) DEALLOCATE(xid_dpnd)
             IF (nx > 0) THEN
                ALLOCATE (xid(nx), xid_dpnd(nx))
                READ(iuo)xid(:), xid_dpnd(:)
                IF ( .not.ANY(xid(1:nx) == 1)) THEN
                   print*,errHeader,'If nx > 0, WV should be among dependent on'
                   call exit(1)
                END IF
                IF ( ANY(xid_dpnd(1:nx) /= 1)) THEN
                   print*, errHeader,'Currently, implementation is done'&
                         ,' for dependence on WV only'
                   call exit(1)
                END IF
             ELSE
                print*, errHeader&
                        ,'At least, one dependence should exist: WV on WV'
                call exit(1)
             END IF
          CASE ('N')
             READ(iuo)
          CASE DEFAULT
             PRINT *, wrnHeader,'Selection file header content flag'&
                    , ' not recognized:',HdrOpt(iHO:iHO)
       END SELECT
    END DO
    READ(iuo)V1m, V2m, nfsmp_ir

    READ(iuo)nfix,nMol !integer*4, conversion to int*2 is below

    READ(iuo) MolIDFix(1:nfix), molID(1:nMol)
    if (any(molID == 81)) then
       isHDO=.true.
    else
       isHDO=.false.
    end if

    nFix=nMol-nmol_in ! nFix is overwritten
    allocate (mFix(nFix)) !variable species to be transferred into Fix
    icnt=0
    !set mFix as those, which are not included into molID_in
    do k=1,nMol
       if ( .not. ANY(molID_in(1:nmol_in) == molID(k)) ) then
          icnt=icnt+1
          mFix(icnt)=molID(k)
       end if
    end do
    if (icnt .NE. nFix) THEN
       print*, errHeader,' Input Variable molecular list is inconsistent'&
             ,' with supplied LUT'
       call exit(1)
    end if
    READ(iuo)nLayOD,nTmpOD,nTself,nWvpOD
    nwvpod_=MAX(nWvpOD-1,1)
    !in this version we use nWvpOD=1: 2-points WV ODs plus linear dependence of Keff
    !(Keff=kh2o+dkh2o*q) has been used to determine kh2o and dkh2o
    IF (nwvpod_ > 1) THEN
       print*, errHeader,' Versions with NWVPOD_ > 1 is not supported'
       call exit(1)
    END IF
    numRefLev=nLayOD+1
    allocate (pRef(numRefLev),pRefLoc(numRefLev))
    allocate (tmpTab(nTmpOD,nLayOD),wvpTab(nwvpod_,nLayOD))
    allocate (tmpLoc(nTmpOD,numRefLev),wvpLoc(nwvpod_,nLayOD))
    allocate (fixDMR(nLayOD),TmpSelf(nTself))
    allocate (dfltLoc(nMol,nLayOD),fixLoc(nLayOD))
    allocate (pavlref(nLayOD))
    allocate (vwvn(nfsmp_ir),f1_arr(nfsmp_ir),f2_arr(nfsmp_ir))
    !wvpTab, mixing ratios for WV
    READ(iuo)pRefLoc(1:numRefLev),tmpLoc(1:nTmpOD,1:nLayOD),TmpSelf(1:nTself),wvpLoc(1:nwvpod_,1:nLayOD)
    pRef(1:numRefLev)=pRefLoc(1:numRefLev)
    tmptab(1:nTmpOD,1:nLayOD)=tmpLoc(1:nTmpOD,1:nLayOD)
    wvptab(1:nwvpod_,1:nLayOD)=wvpLoc(1:nwvpod_,1:nLayOD)

    pavlref(1:nLayOD)=(pRef(2:nLayOD+1)-pRef(1:nLayOD))/log(pRef(2:nLayOD+1)/pRef(1:nLayOD))
    ! Transform level reference temperatures onto a layer grid
    READ(iuo)fixLoc(1:nLayOD),dfltLoc(1:nMol,1:nLayOD)
    fixDMR(1:nLayOD)=fixLoc(1:nLayOD)
    allocate (invMolID(molID(nMol)))
    ! in order to count variable species interms of compressed indices, we set
    ! auxiliary array invMolID
    CALL invertMolID(nMol,molID(1:nMol),invMolID)
    ! File for re-setting defaults contains mixing ratios of variable species, which
    ! could be potentially used for making them fixed. Number of levels for this
    ! default species file and for LUT are to be identical
    ! go to specific quantities
    allocate (wvpTmp(nMol+2,nLayOD))
    wvpTmp = 0.0
    wvpTmp(1,1:nLayOD)=fixDMR(1:nLayOD)
    wvpTmp(2,1:nLayOD)=wvpTab(1,1:nLayOD)
    wvpTmp(3,1:nLayOD)=wvpTab(1,1:nLayOD)
    ! Define wvpTmp in order to use it without changes s/routines of GetOD
    wvpTmp(4:nMol+2,1:nLayOD)=dfltLoc(2:nMol,1:nLayOD)
    deallocate (pRefLoc,tmpLoc,wvpLoc,dfltLoc,fixLoc)
    ! process defProfFile if exists

    if (defProfFile /= 'NONE') then
      inquire (file=defProfFile, exist=isExist)
      if (.not. isExist) THEN
         print*, errHeader,' default profile file not found', trim(defProfFile)
         call exit(1)
      end if

      iuF = findFreeUnit()
      OPEN(iuF,file=defProfFile,status='old')
      READ(iuF,*)nlevF,nselFix
      IF (nlevF /= nLayOD+1) THEN
         print*, errHeader,'Layering in default profile file should be as in LUT'
         call exit(1)
      END IF
      allocate (iselFix(nselFix),vmol(nlevF),qr(nlevF))
      READ(iuF,*)iselFix(1:nselFix)
      DO m=1,nFix
         IF ( .not.ANY(iselFix(1:nselFix) == mFix(m))) THEN
            PRINT*,wrnHeader, ' defProfFile does not match Variable gas'&
                  ,' specified for being transferred to Fix'
            PRINT*, 'default concentration of LUT is unchanged'
         ELSE
            CALL amnt2change(iuF,mFix(m),vmol)
            DO k=1,nLayOD
                scale=drymwt/molWt(mFix(m))/(pRef(k+1) - pRef(k))
                call lpsum_log(pRef(k), pRef(k+1),vmol(k), vmol(k+1), scale, &
                         wvpTmp(2+invMolID(mFix(m)), k), dummyReal,dummyReal)
            end do
         end if
      end do
      close(iuF)
    end if

    !---Build mapping vector for user-selected variable molecules
    kk=0
    map(1:nMol)=0
    DO k=1,nMol
       !---map(k) contains new relative index for selected molecule and is zero for others
       IF ( ANY(molID_in(1:nmol_in) == molID(k))) THEN
          kk=kk+1
          map(k)=kk
       END IF
    END DO
    IF(kk.NE.nmol_in)THEN
       PRINT *,errHeader&
              ,'This database contains only the molecules with'
       PRINT *,'the following HITRAN ID as variable molecules: '
       PRINT *,molID(1:nMol)
       call exit(1)
    END IF
    !---Read oss parameters
    allocate (coef_ir(nChMax,nf_sel))
    allocate (ichmap_ir(nChMax,nf_sel))
    allocate (iselS(nf_sel),nch_ir(nf_sel))
    allocate (coef_Loc(nChMax,nf_sel))
    allocate (ichmap_Loc(nChMax,nf_sel))
    allocate (isel_Loc(nf_sel),nch_Loc(nf_sel))
    DO ismp=1,nf_sel
       READ(ius)isel_Loc(ismp),nch_Loc(ismp)
       READ(ius)coef_Loc(1:nch_Loc(ismp),ismp),ichmap_Loc(1:nch_Loc(ismp),ismp)
    END DO
    CLOSE(ius)
    allocate (nNodes_arr(0:mxIndex))
    allocate (iselS_arr(nf_sel, 0:mxIndex))
    allocate (nch_arr(nf_sel, 0:mxIndex))
    allocate (ichMap_arr(nChMax,nf_sel, 0:mxIndex))
    allocate (coef_arr(nChMax,nf_sel, 0:mxIndex))
    allocate (nChList_arr(0:mxIndex))
    allocate (chanList_arr(nChan, 0:mxIndex))

    lastChSet = UNDEFINED_INDEX
    nchanAll = nChan
    indxChSel=0
    CALL addChanSelect(chIndx, nChan, indxChSel)
    CALL setChanSelect(indxChSel)

    !---Read absorption coefficient tables
    nsize1=nLayOD*nTmpOD
    nsize2=nsize1*nWvpOD
    allocate (kfix_ir(nsize1,nfsmp_ir),kh2o_ir(nsize2,nfsmp_ir))
    allocate (kvar_ir(mxmols*nsize1,nfsmp_ir),kself(nTself,nfsmp_ir))
    allocate (nmols(nfsmp_ir),imols(mxmols,nfsmp_ir))
    ! xid includes 0 (if fixed gases are used as dependent) and HITRAN indices for
    ! other variable species; xid_ uses compressed indices instead
    ! The following cases are under consideration: 1. Fixed gases are in xid, WV is
    ! also in, other variable species could, or, could not be in. 2. Fixed gases are
    ! not in xid, however, WV is in, other variable species could, or, could not be in.
    allocate (xid_(nx))
    allocate (dkfix_ir(MAX(1,nsize1*(1-xid(1))),nfsmp_ir))
    allocate (dkh2o_ir(nsize1,nfsmp_ir))
    allocate (khdo_ir(nsize1,nfsmp_ir),dkhdo_ir(nsize1,nfsmp_ir))

    if (xid(1) == 0) then
       xid_(1) = 0
       xid_(2:nx)=invMolID(xid(2:nx))
       IF ( ANY(xid_(2:nx) == 0)) THEN
          print*, 'Err[oss_ir_module::GetOD]: xid contains indices out of molID'
          call exit(1)
       END IF
    end if
    if (xid(1) == 1) then
       xid_(1:nx) = invMolID(xid(1:nx))
       IF ( ANY(xid_(1:nx) == 0)) THEN
          print*, 'Err[oss_ir_module::GetOD]: xid contains indices out of molID'
          call exit(1)
       END IF
    end if
    nsize4=nsize1*(nx-2+xid(1))
    nsize5=nsize1*(1-xid(1))
    allocate (kfix_Loc(nsize1),dkh2o_Loc(nsize1),kh2o_Loc(nsize2),kself_Loc(nTself))
    allocate (dkfix_Loc(MAX(1,nsize5)))
    ! Transformation of xid from molID indices to compressed indices
    DO ISmp=1,Nfsmp_ir
       READ(iuo)vwvn(ismp),nmols_ !real*8, int*4
       IF (Spc_Units_l == 1) vWvn(ismp)= vWvn(ismp)/SOfL!From GHz 2 wn
       IF (Spc_Units_l == 2) vWvn(ismp)= 1e4/vWvn(ismp)!From mu 2 wn

       f1_arr(ismp)=f1*vwvn(ismp)**3
       f2_arr(ismp)=f2*vwvn(ismp)
       READ(iuo)Imols_(1:nmols_) !int*4
       ! transformation to int*2 is done for using nmols_Tmp and imols_tmp
       ! in odtresh, cum_fix and shrink_var
       nmols_Tmp = nmols_
       Imols_Tmp(1:nmols_)=Imols_(1:nmols_)
       nsize3=nsize1*(nmols_-1)
       allocate (kvar_tmp(max(1,nsize3)),kvar_Loc(max(1,nsize3)))
       READ(iuo)kfix_Loc(1:nsize1),kh2o_Loc(1:nsize2),kvar_Loc(1:nsize3),kself_Loc(1:nTself)

       kfix_ir(1:nsize1,ismp)=kfix_Loc(1:nsize1)
       kh2o_ir(1:nsize2,ismp)=kh2o_Loc(1:nsize2)
       kvar_tmp(1:nsize3)=kvar_Loc(1:nsize3)
       kself(1:nTself,ismp)=kself_Loc(1:nTself)
       ! For dkvar_tmp(:,:) and dkvar_Loc(:,:), the 2nd dimension is, in
       ! general smaller, than nsize4, because some of variable gases
       ! supposed to be within a group of nx-2-xid(1), may be off this group
       ! because of a compression due to imols. In this version dkvar is skipped
       READ(iuo)dkfix_Loc(1:nsize5),dkh2o_Loc(1:nsize1)

       !---code uses kfix_ir*rfix instead of kfix_ir
       dkfix_ir(1:nsize5,ismp)=dkfix_Loc(1:nsize5)
       dkh2o_ir(1:nsize1,ismp)=dkh2o_Loc(1:nsize1)
       if (isHDO) then
          read(iuo)kfix_Loc(1:nsize1), dkh2o_Loc(1:nsize1)
          khdo_ir(1:nsize1,ismp) =kfix_Loc(1:nsize1)
          dkhdo_ir(1:nsize1,ismp)=dkh2o_Loc(1:nsize1)
       end if
       IF (fmtLUT<2) then
         CALL init_kfix(kfix_ir(1,ismp),wvpTmp)
         IF (xid(1) == 0) THEN
            CALL init_kfix(dkfix_ir(1,ismp),wvpTmp)
         END IF
       end if
       !     Keep only user-selected molecules as variables molecules -
       !     merge remaining molecules with fixed gases
       !     First apply OD threshold to identify weakly absorbing constituents
       iflag(2:nMol) = 0
       ! to eleminate non-contributing species
       ! useless since odfac is set to 0
       CALL odthresh(kfix_ir(1,ismp),kvar_tmp,imols_tmp,nmols_tmp,wvpTmp,odfac,iflag)
       maps(1:nMol) = map(1:nMol) * iflag(1:nMol)
       kk = 0
       DO ks = 1,nmols_
          imol = imols_(ks)
          IF(mapS(imol) > 0)THEN
             kk = kk+1
             imols_indx(kk) = ks
          ELSE
              CALL cum_fix(kfix_ir(1,ismp),kvar_tmp,imol,nmols_tmp,ks,wvpTmp)
          END IF
       END DO
       nmols(ismp) = kk
       CALL shrink_var(kvar_tmp,imols_tmp,nmols_tmp,mapS,imols_indx,&
            kvar_ir(1,ISmp),imols(1,ISmp),nmols(ISmp))
       if (allocated(kvar_tmp)) deallocate(kvar_tmp)
       if (allocated(kvar_Loc)) deallocate(kvar_Loc)
    END DO
    nMol = nmol_in

    molID(1:nMol) = molID_in(1:nMol)
    CLOSE(iuo)
    deallocate (kfix_Loc,kh2o_Loc,dkfix_Loc,dkh2o_Loc)
    deallocate (dkfix_ir, fixDMR, Wvptab)
    if (allocated(chIndx))        deallocate(chIndx)
    if (allocated(chID))          deallocate(chID)
    if (allocated(iPol))          deallocate(iPol)
    getODdone=.true.
    RETURN
  END SUBROUTINE GetOD

!---------------------------------------------------------------------------------------------
! The subroutine odthresh finds species (water vapor excluded)
! which contribution to optical depth is less than odfac of total
! and set flag in array iflag

  SUBROUTINE odthresh(kfix0,kvar0,imols,nmols,w,odfac,iflag)
    !---Input variables
    INTEGER  ,                            INTENT(IN)   :: NmolS,ImolS(NmolS)
    real,DIMENSION(nLayOD,nTmpOD),        INTENT(IN)   :: kfix0
    real,DIMENSION(NmolS-1,nLayOD,nTmpOD),INTENT(IN)   :: kvar0
    real,                                 INTENT(IN)   :: odfac
    real, DIMENSION(:,:),                 INTENT(IN)   :: w
    !---Output variables
    INTEGER, DIMENSION(:), INTENT(INOUT)               :: iflag
    !---Local variables
    REAL                                   :: kbuf(50),ktot
    INTEGER                                :: nt,l,ks
    !     Compare optical depth of individual dry constituents to
    !     total optical depth and flag weakly absorbing molecules
    DO nt = 1,nTmpOD
       DO l = 1,nLayOD
          ktot = kfix0(l,nt)
          DO ks = 2,nmols
             kbuf(ks-1) = kvar0(ks-1,l,nt) * w(imols(ks)+2,l)
             ktot = ktot + kbuf(ks-1)
          END DO
          DO ks=2,nmols  ! could exit at first occurence of iflag = 1
             IF(kbuf(ks-1) > odfac * ktot) iflag(imols(ks)) = 1
          END DO
       END DO
    END DO
    RETURN
  END SUBROUTINE odthresh

!---------------------------------------------------------------------------------------------
! The subroutine cum_fix moves some variable species to the fixed part

  SUBROUTINE cum_fix(kfix0,kvar0,imol,nmols,ks,w)
    !---In/Out variables
    INTEGER  ,           INTENT(IN) :: ks
    INTEGER  ,           INTENT(IN) :: NmolS,Imol
    real, DIMENSION(:,:),INTENT(IN) :: w
    real,                INTENT(IN) :: kvar0(NmolS-1,nLayOD,nTmpOD)
    real,             INTENT(INOUT) :: kfix0(nLayOD,nTmpOD)
    !---Local variables
    INTEGER    :: nt
    DO nt=1,nTmpOD
       kfix0(1:nLayOD,nt) = kfix0(1:nLayOD,nt) +&
            kvar0(ks-1,1:nLayOD,nt) * w(imol+2,1:nLayOD)
    END DO
    RETURN
  END SUBROUTINE cum_fix

  SUBROUTINE init_kfix(kfix0,w)
    !---In/Out variables
    real, DIMENSION(:,:),INTENT(IN) :: w
    real,             INTENT(INOUT) :: kfix0(nLayOD,nTmpOD)
    !---Local variables
    INTEGER    :: nt
    DO nt=1,nTmpOD
       kfix0(1:nLayOD,nt) = kfix0(1:nLayOD,nt) * w(1,1:nLayOD)
    END DO
    RETURN
  END SUBROUTINE init_kfix

!---------------------------------------------------------------------------------------------
! The subroutine shrink_var reduce the size of kvar
! and makes a new map of molecular indices

  SUBROUTINE shrink_var(kvar0,imols0,nmols0,maps,imols_indx,kvar1,imols1,nmols1)
    !---In/Out variables
    INTEGER              , INTENT(IN) :: nmols0,imols0(nmols0),nmols1
    INTEGER, DIMENSION(:), INTENT(IN) :: maps,imols_indx
    REAL                 , INTENT(IN) :: kvar0(nmols0-1,nLayOD,nTmpOD)
    REAL              , INTENT(INOUT) :: kvar1(nmols1-1,nLayOD,nTmpOD)
    INTEGER           , INTENT(INOUT) :: imols1(nmols1)

    kvar1(1:nmols1-1,1:nLayOD,1:nTmpOD) = &
         kvar0(imols_indx(2:nmols1)-1,1:nLayOD,1:nTmpOD)
    imols1(1:nmols1) = maps(imols0(imols_indx(1:nmols1)))
    RETURN
  END SUBROUTINE shrink_var

  !----------------------------------------------------------------------------
  ! PURPOSE: Computes items related to the viewing geometry.
  !----------------------------------------------------------------------------

  SUBROUTINE setpath_ir(pSurf,pLoc,nObs,obsang,sunang,nSurf,sun,umu,umu0,lookup,&
             referenceGrid, interpSfc,pUser)
    real,                    PARAMETER     :: SUN_HORIZON = 85.0
    !---Input variables
    real,                       INTENT(IN) :: pSurf,obsang,sunang
    real, DIMENSION(:),OPTIONAL,INTENT(IN) :: pUser
    !---Output variables
    real, DIMENSION(:),      INTENT(INOUT) :: pLoc
    INTEGER,                 INTENT(INOUT) :: nObs
    INTEGER,                 INTENT(INOUT) :: nSurf
    LOGICAL,                 INTENT(INOUT) :: sun
    LOGICAL,                 INTENT(INOUT) :: lookup
    real,                    INTENT(INOUT) :: umu,umu0
    LOGICAL,                 INTENT(INOUT) :: referenceGrid,interpSfc
    !---Local variables
    REAL                 :: viewang
    INTEGER              :: npLev

    !---Set viewing geometry and pressure grid parameters
    if (present(pUser)) then
       referenceGrid=.false.
      nPLev=size(pUser)
       pLoc(1:npLev)=pUser(:)
      if (abs(pSurf-pUser(nPLev)) < (pUser(nPLev)-pUser(nPLev-1))* &
           pMatchLim) then
          interpSfc =.false.
      else
          interpSfc = .true.
      end if
    else
       referenceGrid = .true.
      nPLev=numRefLev
       pLoc(1:nplev)=pRef(:)
       interpSfc = .true.
    end if
    if (interpSfc) then
      DO nSurf=2,nPLev-1
         IF (pLoc(nSurf) >= pSurf) EXIT
      END DO
    else
      nSurf=nPLev
    end if

    if (nObs <= -1) then
       nObs = nSurf
    elseif (nObs == 0) then
       nObs = 1
    elseif (nObs > nSurf) then
       print*,'Err[oss_ir_module::setpath_ir]: nObs is beyond of its limit'
      call exit(1)
    end if

    viewang = obsang
    if (viewang > 90.0) then
       lookup=.TRUE.
       viewang=180.-viewang
    else
       lookup=.FALSE.
    end if
    if (lookup .and. nObs == 1) THEN
       print*, 'Err[oss_ir_module::setpath_ir]: ', &
               ' Lookup does not work for satellite'
       call exit(1)
    end if
    umu =COS(viewang*deg2rad)
    IF (sunang < SUN_HORIZON)THEN
       sun=.TRUE.
       umu0=COS(sunang*deg2rad)
    ELSE
       sun=.FALSE.
       umu0=1.
    END IF
    RETURN
  END SUBROUTINE setpath_ir

!-------------------------------------------------------------------------
! layerAverage  is responsible for calculus of
! average pressure, mu and mu_solar as a function of altitude
!
  SUBROUTINE layerAverage(xG,N2,referenceGrid,interpSfc,&
             tempIndex,pSurfIndex,pLoc,pavl, &
             plogu,plogl, umu,umuLay,umu0,umu0Lay,alt, nEnd, alpha)
    real, dimension(:)         , intent(in) :: xG
    real                       , intent(in) :: umu,umu0
    real, dimension(:)         , intent(in) :: pLoc
    integer                    , intent(in) :: N2
    logical                    , intent(in) :: referenceGrid,interpSfc
    integer                    , intent(in) :: tempIndex,pSurfIndex
    real, dimension(mxlay)     , intent(in) :: alt
    integer                    , intent(out) :: nEnd
    !---Output variables
    real, dimension(:)           ,intent(inout):: pavl
    real, dimension(:)           ,intent(inout):: umuLay,umu0Lay
    real                      ,intent(out) :: plogl
    real                      ,intent(out) :: plogu
    real                      ,intent(out) :: alpha
    !---Local variables
    integer                         :: l
    real                            :: rRatioSq

    if ( .not.referenceGrid) then
        pavl(1:n2) = &
           (pLoc(2:n2+1)-pLoc(1:n2))/log(pLoc(2:n2+1)/pLoc(1:n2))
    end if

    plogu = log(xG(pSurfIndex)/pLoc(n2))
    plogl = log(xG(pSurfIndex)/pLoc(n2+1))
    alpha = plogl/(plogl - plogu)

    if (interpSfc) THEN
        nEnd = N2-1
        pavl(n2)=(xG(pSurfIndex)-pLoc(n2))/plogu
    else
        nEnd = N2
    end if

    if (sphericalGeometryFlag) then
    ! Here set for "secants" is done in accordance with LBL
       do l=1,n2
         rRatioSq = ((rEarth+alt(n2+1))/(rEarth+0.5*(alt(l)+alt(l+1))))**2
         umuLay(l) = sqrt(1. - rRatioSq*(1. - umu**2))
         umu0Lay(l) = sqrt(1. - rRatioSq*(1. - umu0**2))
       end do
    else
        !plane parallel version
        umuLay(1:N2) = umu
        umu0Lay(1:N2) = umu0
    end if

  END SUBROUTINE layerAverage

!----------------------------------------------------------------------------
! PURPOSE:  This subroutine calculates the average temperature and molecular
!           amounts for all layers for given profiles of temperature and
!           molecular concentrations. It also calculates the derivatives of
!           tavl with respect to a change in the lower and upper boundary
!           temperatures and the derivatives of the molecular amounts with
!           respect to a change in the mixing ratios at the layer
!           boundaries. Molecular amounts are in molec./cm**2.
!           Integration assumes that T is linear in z (LnT linear in LnP)
!           and LnQ linear in LnP.
!           it implements integration over altitude
!----------------------------------------------------------------------------
  SUBROUTINE fpathZ(xG,N2,nEnd, interpSfc, tempIndex, pSurfIndex,varMolIndex,&
             pLoc, tSfc,tavl,dtu,dtl, alpha, wfix,q,w,dwqu,dwql,wvdwqu,wvdwql,alt)
    !---Input variables
    real, dimension(:)         , intent(in) :: xG
    real, dimension(:)         , intent(in) :: pLoc
    integer                    , intent(in) :: nEnd
    integer                    , intent(in) :: N2
    logical                    , intent(in) :: interpSfc
    integer                    , intent(in) :: tempIndex
    integer                    , intent(in) :: pSurfIndex
    integer, dimension(:)      , intent(in) :: varMolIndex
    real                       , intent(in) :: alpha
    !---Output variables
    real                         ,intent(out)  :: tSfc
    real, dimension(:)           ,intent(inout):: tavl,dtu,dtl
    real, dimension(:)           ,intent(inout):: wfix
    real, dimension(:,:)         ,intent(inout):: q,w
    real, dimension(:,0:)        ,intent(inout):: dwqu,dwql
    real, dimension(:,:)         ,intent(inout):: wvdwqu,wvdwql
    !---Local variables
    integer                         :: l,iH2o,k,iXoff,n
    real                            :: pSfc,qsfc,dzsfc
    real, dimension(mxlev)          :: qcor,qr
    real                            :: RHOTOT1,RHOTOT2
    real                            :: wtot
    real, dimension(mxlay)          :: alt

    !------------------------------------------------------------------------
    !     Compute Average Temperature for the layers
    !------------------------------------------------------------------------
    pSfc = xG(pSurfIndex)
    iXoff = tempIndex-1

    do l=1,nEnd
      call lzsum4T(xG(iXoff+l),xG(iXoff+l+1),tavl(l),dtu(l),dtl(l))
    end do

    if (interpSfc) then
      !---Surface layer
      tSfc=(xG(iXoff+N2)**alpha) * (xG(iXoff+N2+1)**(1.-alpha))
      call lzsum4T(xG(iXoff+N2),tSfc,tavl(N2),dtu(N2),dtl(N2))

      dtu(N2)=dtu(N2)+dtl(N2)*alpha*tSfc/xG(iXoff+N2)
      dtl(N2)=dtl(N2)*(1.-alpha)*tSfc/xG(iXoff+N2+1)
    else
      tSfc=xG(tempIndex + N2)
    end if

    iH2o=varMolIndex(1)-1

    !------------------------------------------------------------------------
    !     Calculate amounts for individual species and derivatives wrt mixing
    !     ratios for retrieved constituents.
    !------------------------------------------------------------------------
    qcor(1:N2+1)=1./(1. + xG(iH2o+1:iH2o+N2+1)/normMolWt(molID(1)) )
    do k=1,nMol
      iXoff=varMolIndex(k)-1
      !---Transform mix. ratios into volume mix fractions
      !  (relative to total air mol number)
      qr(1:N2+1)=qcor(1:N2+1)*xG(iXoff+1:iXoff+N2+1)/normMolWt(molID(k))

      call totalDensity(pLoc(1),xG(1),RHOTOT1)
      do l=1,nEnd
        call totalDensity(pLoc(l+1),xG(l+1),RHOTOT2)
        if (k == 1) THEN

             call lzsum4W(qr(l), qr(l + 1), alt(l) - alt(l + 1), &
                 RHOTOT1, RHOTOT2, w(k, l), dwqu(l, k), dwql(l, k), wtot)

             wfix(l) = wtot - w(k, l)
             q(1, l) = w(1, l)/wtot
          else
             call lzsum4W(qr(l), qr(l + 1), alt(l) - alt(l + 1), &
                 RHOTOT1, RHOTOT2, w(k, l), dwqu(l, k), dwql(l, k))
             q(k, l) = w(k, l)/wfix(l)
          end if
          RHOTOT1 = RHOTOT2
      end do
      if (interpSfc) then
        !---Surface Layer
        qsfc=(qr(N2)**alpha)*(qr(N2+1)**(1.-alpha))
        call totalDensity(pSfc,tSfc,RHOTOT2)
        if (k == 1) THEN
             dzsfc= alt(N2) - alt(N2+1)
             call lzsum4W(qr(N2), qsfc, dzsfc, &
                 RHOTOT1, RHOTOT2, w(k, N2), dwqu(N2, k), dwql(N2, k), wtot)
             wfix(N2) = wtot - w(1, N2)
             q(1, N2) = w(1, N2)/wtot
        else
             dzsfc= alt(N2) - alt(N2+1)
             call lzsum4W(qr(N2), qsfc, dzsfc, &
                 RHOTOT1, RHOTOT2, w(k, N2), dwqu(N2, k), dwql(N2, k))
             q(k, N2) = w(k, N2) / wfix(N2)
        end if
        dwqu(N2, k) = dwqu(N2, k) + dwql(N2, k) * alpha * qsfc/qr(N2)
        dwql(N2, k) = dwql(N2, k)*(1. - alpha) * qsfc/qr(N2 + 1)
      end if
    end do
    !---Derivatives of amount wrt dry mixing ratios
    do n=1,N2
      dwqu(n,1)=dwqu(n,1)*qcor(n)**2/normMolWt(molID(1))
      dwql(n,1)=dwql(n,1)*qcor(n+1)**2/normMolWt(molID(1))
      dwqu(n, 0) = - dwqu(n, 1)
      dwql(n, 0) = - dwql(n, 1)
      dwqu(n,2:nMol)=dwqu(n,2:nMol)*qcor(n)/normMolWt(molID(2:nMol))
      dwql(n,2:nMol)=dwql(n,2:nMol)*qcor(n+1)/normMolWt(molID(2:nMol))
    end do
    ! additional contribution to derivatives of
    !  amount wrt water vapor dry mixing ratio
    do k=2, nMol
       iXoff=varMolIndex(k)-1
       do n=1,N2
          wvdwqu(n,k)=dwqu(n,k)*qcor(n)  *xG(iXoff+n)/normMolWt(molID(1))
          wvdwql(n,k)=dwql(n,k)*qcor(n+1)*xG(iXoff+n+1)/normMolWt(molID(1))
       end do
    end do
    return
  END SUBROUTINE fpathZ

  !----------------------------------------------------------------------------
  ! PURPOSE:  This subroutine calculates the average temperature and molecular
  !           amounts for all layers for given profiles of temperature and
  !           molecular concentrations. It also calculates the derivatives of
  !           tavl with respect to a change in the lower and upper boundary
  !           temperatures and the derivatives of the molecular amounts with
  !           respect to a change in the mixing ratios at the layer
  !           boundaries. Molecular amounts are in molec./cm**2.
  !           Integration assumes that T is linear in z (LnT linear in LnP)
  !           and LnQ linear in LnP.
  !           it implements integration over P
  !----------------------------------------------------------------------------
  SUBROUTINE fpathP(xG,lat,N2,nEnd,interpSfc,&
             pSurfIndex, tempIndex, varMolIndex,pLoc, tSfc,tavl,dtu,dtl, alpha, &
             wfix,q,w,dwqu,dwql,wvdwqu,wvdwql,alt)
    !---Input variables
    real, dimension(:)         , intent(in) :: xG
    real                       , intent(in) :: lat
    real, dimension(:)         , intent(in) :: pLoc
    integer                    , intent(in) :: nEnd
    integer                    , intent(in) :: N2
    integer                    , intent(in) :: pSurfIndex
    integer                    , intent(in) :: tempIndex
    logical                    , intent(in) :: interpSfc
    integer, dimension(:)      , intent(in) :: varMolIndex
    real                       , intent(in) :: alpha
    !---Output variables
    real                         ,intent(out)  :: tSfc
    real, dimension(:)           ,intent(inout):: tavl,dtu,dtl
    real, dimension(:)           ,intent(inout):: wfix
    real, dimension(:,:)         ,intent(inout):: q,w
    real, dimension(:,0:)        ,intent(inout):: dwqu,dwql
    real, dimension(:,:)         ,intent(inout):: wvdwqu,wvdwql
    !---Local variables
    integer                         :: l,iH2o,k,iXoff,n
    real                            :: qsfc
    real, dimension(mxlev)          :: qcor,qr
    real, dimension(mxlay)          :: rairLay !converts p(mb) in g/cm**2
    real                            :: wtot
    real, dimension(mxlay)          :: alt
    real                            :: rair, scalBase
    real                            :: scal,  scall

    !------------------------------------------------------------------------
    !     Compute Average Temperature for the layers
    !------------------------------------------------------------------------
    iXoff = tempIndex - 1
    do l=1,nEnd
      scal = 1./(pLoc(l+1) - pLoc(l))
      call lpsum_log(pLoc(l), pLoc(l+1), xG(iXoff+l),xG(iXoff+l+1), &
             scal, tavl(l),dtu(l),dtl(l))
    end do

    if (interpSfc) then
      !---Surface layer
      tSfc=(xG(iXoff+N2)**alpha) * (xG(iXoff+N2+1)**(1.-alpha))
      scal = 1./(xG(pSurfIndex) - pLoc(N2))

      call lpsum_log(pLoc(N2), xG(pSurfIndex), xG(iXoff+N2),tSfc, &
             scal, tavl(N2),dtu(N2),dtl(N2))

      dtu(N2)=dtu(N2)+dtl(N2)*alpha*tSfc/xG(iXoff+N2)
      dtl(N2)=dtl(N2)*(1.-alpha)*tSfc/xG(iXoff+N2+1)
    else
      tSfc=xG(iXoff+N2+1)
    end if

    iH2o=varMolIndex(1)-1

    rair=10./(grav_const_req1 + grav_const_req2*COS(2.0*deg2rad*LAT))
    !---Derivatives of amount wrt dry mixing ratios
    qcor(1:N2+1)=1./(1.+xG(iH2o+1:iH2o+N2+1))
    rairLay(1:N2+1) = qcor(1:N2+1)*(1. + alt(1:N2+1)/rEarth)**2
    !------------------------------------------------------------------------
    !     Calculate amounts for individual species and derivatives wrt mixing
    !     ratios for retrieved constituents.
    !------------------------------------------------------------------------
    scalBase =  rair*AVOGAD/ drymwt
    do k=1,nMol
        scal = scalBase/normMolWt(molID(k))
        iXoff=varMolIndex(k)-1
        !---Transform mix. ratios into volume fractions (relative to total air molecular number)
        qr(1:N2+1)=rairLay(1:N2+1)*xG(iXoff+1:iXoff+N2+1)
        do l=1,nEnd
          call lpsum_log(pLoc(l), pLoc(l+1), qr(l), qr(l+1), scal, w(k, l), dwqu(l, k), dwql(l, k))
          if (k == 1) THEN
               call lpsum_log(pLoc(l), pLoc(l+1), rairLay(l), rairLay(l+1), scalBase, wtot, dwqu(l, 0), dwql(l, 0))
               wfix(l) = wtot
               q(1, l) = w(1, l)/(wtot+w(1, l))
            else
               q(k, l) = w(k, l)/wfix(l)
            end if
        end do
        if (interpSfc) then
          !---Surface Layer
          qsfc=(qr(N2)**alpha)*(qr(N2+1)**(1.-alpha))
          call lpsum_log(pLoc(N2), xG(pSurfIndex), qr(N2), qsfc, scal, w(k, N2), dwqu(N2, k), dwql(N2, k))
          if (k == 1) THEN
               scall=(rairLay(N2)**alpha)*(rairLay(N2+1)**(1.-alpha))
               call lpsum_log(pLoc(N2), xG(pSurfIndex), rairLay(N2), scall, scalBase, wtot, dwqu(N2, 0), dwql(N2, 0))
               wfix(N2) = wtot
               q(1, N2) = w(1, N2) / (wtot + w(1, N2))
            else
               q(k, N2) = w(k, N2) / wfix(N2)
            end if
            dwqu(N2, k) = dwqu(N2, k) + dwql(N2, k) * alpha * qsfc/qr(N2)
            dwql(N2, k) = dwql(N2, k)*(1. - alpha) * qsfc/qr(N2 + 1)
        end if
    end do

    do n=1,N2
       dwqu(n,0)= -dwqu(n,0)*rairLay(n)*qcor(n)
       dwql(n,0)= -dwql(n,0)*rairLay(n+1)*qcor(n+1)

       dwqu(n,1)=dwqu(n,1)*rairLay(n)*qcor(n)
       dwql(n,1)=dwql(n,1)*rairLay(n+1)*qcor(n+1)

       dwqu(n,2:nMol)=dwqu(n,2:nMol)*rairLay(n)
       dwql(n,2:nMol)=dwql(n,2:nMol)*rairLay(n+1)
    end do
    do k=2, nMol
        iXoff=varMolIndex(k)-1
        do n=1,N2
            wvdwqu(n,k)=dwqu(n,k)*qcor(n)  *xG(iXoff+n)
            wvdwql(n,k)=dwql(n,k)*qcor(n+1)*xG(iXoff+n+1)
        end do
    end do
    iXoff=varMolIndex(1)-1
    return
  end subroutine fpathP

  !----------------------------------------------------------------------------
  ! PURPOSE: Computes temperature/Water vapor indexes and interpolation
  !          coefficients.
  !----------------------------------------------------------------------------
  SUBROUTINE settabindx_ir(tavl,pavl,N2,referenceGrid, &
       indxt_p1,indxt_p2,indxp,ap1,ap2,&
       at1_p1,at2_p1,at1_p2,at2_p2,adt1_p1,adt2_p1,adt1_p2,adt2_p2)
    !---Input variables
    INTEGER                 ,INTENT(IN) :: n2
    real,    DIMENSION(:)   ,INTENT(IN) :: tavl,pavl
    LOGICAL                 ,INTENT(IN) :: referenceGrid
    !---Output variables
    INTEGER, DIMENSION(:),INTENT(INOUT) :: indxt_p1,indxt_p2,indxp
    real,    DIMENSION(:),INTENT(INOUT) :: at1_p1,at2_p1,at1_p2,at2_p2,ap1,ap2
    real,    DIMENSION(:),INTENT(INOUT) :: adt1_p1,adt2_p1,adt1_p2,adt2_p2
    !---Local variables
    INTEGER               :: l,i,j1,j2,lp1,lp2
    REAL                  :: dent1_p1,dent2_p1,dent1_p2,dent2_p2
    lp2 = 0
    !---Compute coefficients for temperature interpolation of ODs
    !!! INDXP: is related to the index of OD pressure layer located just above
    !!!        (by altitude) the current profile layer
    DO l=1,N2
       if (referenceGrid) then
          lp1    = l
          ap2(l) = 0.
          ap1(l) = 1.
       else
          lp1=1
          DO WHILE(pavlref(lp1)<pavl(l))
             lp1=lp1+1
          END DO
          lp2=lp1
          lp1=lp1-1
          if (lp1==0) then
             lp1=1
             lp2=2
          end if
          if (lp2>=numRefLev) then
             lp1=numRefLev-2
             lp2=numRefLev-1
          end if
          !!! if pavl(l) becomes located between lp1 and lp2:
          ap2(l)=(pavl(l)-pavlref(lp1))/(pavlref(lp2)-pavlref(lp1))
          ap1(l)=(pavlref(lp2)-pavl(l))/(pavlref(lp2)-pavlref(lp1))
       end if
       indxp(l)=lp1

       !!! INDXT_P1: is related to the midst temperature index within the
       !!!           upper (by altitude) OD pressure layer INDXP
       !!! INDXT_P2: is related to the midst temperature index within the
       !!!           lower (by altitude) OD pressure layer INDXP+1
       TempLoop1: DO i=3,nTmpOD-2
          IF(tmptab(i,lp1)>=tavl(l)) exit TempLoop1
       END DO TempLoop1
       IF (ABS(tavl(l)-tmptab(i-1,lp1))<=ABS(tavl(l)-tmptab(i,lp1))) THEN
          j1=i-1
       ELSE
          j1=i
       END IF
       indxt_p1(l)=j1

       dent1_p1=(tmptab(j1-1,lp1)-tmptab(j1,lp1))*(tmptab(j1-1,lp1)-tmptab(j1+1,lp1))
       dent2_p1=(tmptab(j1,lp1)-tmptab(j1-1,lp1))*(tmptab(j1,lp1)-tmptab(j1+1,lp1))
       at1_p1(l)=(tavl(l)-tmptab(j1,lp1))*(tavl(l)-tmptab(j1+1,lp1))/dent1_p1
       at2_p1(l)=(tavl(l)-tmptab(j1-1,lp1))*(tavl(l)-tmptab(j1+1,lp1))/dent2_p1
       adt1_p1(l)=(2*tavl(l)-tmptab(j1,lp1)-tmptab(j1+1,lp1))/dent1_p1
       adt2_p1(l)=(2*tavl(l)-tmptab(j1-1,lp1)-tmptab(j1+1,lp1))/dent2_p1
       if (.not.referenceGrid) then
          TempLoop2: DO i=3,nTmpOD-2
             IF(tmptab(i,lp2)>=tavl(l)) exit TempLoop2
          END DO TempLoop2
          IF (ABS(tavl(l)-tmptab(i-1,lp2))<=ABS(tavl(l)-tmptab(i,lp2))) THEN
             j2=i-1
          ELSE
             j2=i
          END IF
          indxt_p2(l)=j2
          !!! aty_px: interpolation coefficients, which are related to pressure
          !!!         layer (x=1,or,2, upper, or, lower) and temperature
          !!!         interpolation points (y=1,or,2, lower,or, midst)
          dent1_p2=(tmptab(j2-1,lp2)-tmptab(j2,lp2))*(tmptab(j2-1,lp2)-tmptab(j2+1,lp2))
          dent2_p2=(tmptab(j2,lp2)-tmptab(j2-1,lp2))*(tmptab(j2,lp2)-tmptab(j2+1,lp2))
          at1_p2(l)=(tavl(l)-tmptab(j2,lp2))*(tavl(l)-tmptab(j2+1,lp2))/dent1_p2
          at2_p2(l)=(tavl(l)-tmptab(j2-1,lp2))*(tavl(l)-tmptab(j2+1,lp2))/dent2_p2
          adt1_p2(l)=(2*tavl(l)-tmptab(j2,lp2)-tmptab(j2+1,lp2))/dent1_p2
          adt2_p2(l)=(2*tavl(l)-tmptab(j2-1,lp2)-tmptab(j2+1,lp2))/dent2_p2
       end if
    END DO
    RETURN
  END SUBROUTINE settabindx_ir

  !----------------------------------------------------------------------------
  ! PURPOSE:  Given the LUTs, the function computes the layer
  !           optical depths. Refer to ossrad_mw for the transmittance and
  !           brightness temperature calculation
  !----------------------------------------------------------------------------
  SUBROUTINE OSStran(kfix,kh2o,dkh2o,kvar,kself,indxt_p1,indxt_p2,indxp,indxt,&
    ap1,ap2,at1_p1,at2_p1,at1_p2,at2_p2,at1,at2,adt1_p1,adt2_p1,adt1_p2,adt2_p2,&
    N2,referenceGrid,NmolS,ImolS,pavl,tavl,wfix,q,w,tautot,abso,dtaudtmp)
    !---Input variables
    INTEGER                          ,INTENT(IN):: NmolS,ImolS(NmolS)
    INTEGER                          ,INTENT(IN):: n2
    INTEGER, DIMENSION(:)            ,INTENT(IN):: indxt_p1,indxt_p2,indxp
    INTEGER, DIMENSION(:)            ,INTENT(IN):: indxt
    real,    DIMENSION(:)            ,INTENT(IN):: at1_p1,at2_p1,at1_p2,at2_p2
    real,    DIMENSION(:)            ,INTENT(IN):: adt1_p1,adt2_p1
    real,    DIMENSION(:)            ,INTENT(IN):: adt1_p2,adt2_p2
    real,    DIMENSION(:)            ,INTENT(IN):: wfix,ap1,ap2
    real,    DIMENSION(:)            ,INTENT(IN):: at1,at2
    real,    DIMENSION(:)            ,INTENT(IN):: pavl,tavl
    real,    DIMENSION(:,:)          ,INTENT(IN):: q,w
    real,    DIMENSION(nLayOD,nTmpOD),INTENT(IN):: kfix,kh2o,dkh2o
    REAL                             ,INTENT(IN):: kvar(NmolS-1,nLayOD,nTmpOD)
    REAL                             ,INTENT(IN):: kself(nTself)
    LOGICAL                          ,INTENT(IN):: referenceGrid
    !---Output variables
    real,    DIMENSION(:)         ,INTENT(INOUT):: tautot,dtaudtmp
    real,    DIMENSION(0:,:)      ,INTENT(INOUT):: abso
    !---Local variables
    INTEGER                            :: indxt1_p1,indxt2_p1,indxt3_p1
    INTEGER                            :: indxt1_p2,indxt2_p2,indxt3_p2
    INTEGER                            :: l,ks,imol,indxp1,indxp2
    INTEGER                            :: indxt1,indxt2
    REAL                               :: dft1_p1,dft2_p1,dft1_p2,dft2_p2
    REAL                               :: abs0_p1,abs0_p2,dabs0_p1,dabs0_p2
    REAL                               :: at1_p1L,at2_p1L,adt1_p1L,adt2_p1L
    REAL                               :: at1_p2L,at2_p2L,adt1_p2L,adt2_p2L
    REAL                               :: ap1L,ap2L,at1L,at2L
    REAL                               :: k_wv_self,dtausdtmp
    REAL                               :: k_wv_p1, dkdt_wv_p1
    REAL                               :: k_wv_p2, dkdt_wv_p2
    REAL                               :: dk_wv_p1, ddkdt_wv_p1
    REAL                               :: dk_wv_p2, ddkdt_wv_p2
    REAL                               :: k_fix, k_wv, dk_wv
    REAL                               :: dkdt_fix, dkdt_wv, ddkdt_wv

    indxp2 = 0
    indxt1_p2 = 0
    indxt2_p2 = 0
    indxt3_p2 = 0

    DO l=1,N2

       indxt1_p1 = indxt_p1(l)-1
       indxt2_p1 = indxt1_p1+1
       indxt3_p1 = indxt1_p1+2
       indxp1   = indxp(l)
       indxt1   =indxt(l)
       indxt2   =indxt(l)+1
       at1_p1L  = at1_p1(l)
       at2_p1L  = at2_p1(l)
       at1_p2L  = at1_p2(l)
       at2_p2L  = at2_p2(l)
       adt1_p1L = adt1_p1(l)
       adt2_p1L = adt2_p1(l)
       adt1_p2L = adt1_p2(l)
       adt2_p2L = adt2_p2(l)
       ap1L     = ap1(l)
       ap2L     = ap2(l)
       at1L     = at1(l)
       at2L     = at2(l)
       !---Fixed gases

       call threePointInterpolation(at1_p1L, at2_p1L, adt1_p1L, adt2_p1L, &
                kfix(indxp1,indxt1_p1), kfix(indxp1,indxt2_p1),&
                kfix(indxp1,indxt3_p1), abs0_p1, dabs0_p1)

       call threePointInterpolation(at1_p1L, at2_p1L, adt1_p1L, adt2_p1L, &
                kh2o(indxp1,indxt1_p1), kh2o(indxp1,indxt2_p1),&
                kh2o(indxp1,indxt3_p1), k_wv_p1, dkdt_wv_p1)

       call threePointInterpolation(at1_p1L, at2_p1L, adt1_p1L, adt2_p1L, &
                dkh2o(indxp1,indxt1_p1), dkh2o(indxp1,indxt2_p1),&
                dkh2o(indxp1,indxt3_p1), dk_wv_p1, ddkdt_wv_p1)
       !---Water vapor

       if (referenceGrid) then
         k_fix = abs0_p1
         dkdt_fix = dabs0_p1

         k_wv  = k_wv_p1
         dk_wv = dk_wv_p1
         dkdt_wv = dkdt_wv_p1
         ddkdt_wv = ddkdt_wv_p1
         dabs0_p2=0.
         abs0_p2=0.
       else
          indxt1_p2 = indxt_p2(l)-1
          indxt2_p2 = indxt1_p2+1
          indxt3_p2 = indxt1_p2+2
          indxp2 = indxp1+1
          !---Fixed gases
          call threePointInterpolation(at1_p2L, at2_p2L, adt1_p2L, adt2_p2L, &
                kfix(indxp2,indxt1_p2), kfix(indxp2,indxt2_p2),&
                kfix(indxp2,indxt3_p2), abs0_p2, dabs0_p2)

          call threePointInterpolation(at1_p2L, at2_p2L, adt1_p2L, adt2_p2L, &
                kh2o(indxp2,indxt1_p2), kh2o(indxp2,indxt2_p2),&
                kh2o(indxp2,indxt3_p2), k_wv_p2, dkdt_wv_p2)

          call threePointInterpolation(at1_p2L, at2_p2L, adt1_p2L, adt2_p2L, &
                dkh2o(indxp2,indxt1_p2), dkh2o(indxp2,indxt2_p2),&
                dkh2o(indxp2,indxt3_p2), dk_wv_p2, ddkdt_wv_p2)

         k_fix = abs0_p1*ap1L + abs0_p2*ap2L
         dkdt_fix= dabs0_p1*ap1L + dabs0_p2*ap2L

         k_wv  = k_wv_p1*ap1L + k_wv_p2*ap2L
         dk_wv = dk_wv_p1*ap1L + dk_wv_p2*ap2L
         dkdt_wv = dkdt_wv_p1*ap1L + dkdt_wv_p2*ap2L
         ddkdt_wv = ddkdt_wv_p1*ap1L + ddkdt_wv_p2*ap2L
       end if

       ! kself
       k_wv_self   = (kself(indxt1)*at1L+kself(indxt2)*at2L)
       dtausdtmp = (kself(indxt2)-kself(indxt1))&
                    /(tmpself(indxt2)-tmpself(indxt1))

       if (referenceGrid) then
            k_wv_self = k_wv_self * pavlref(l)
            dtausdtmp = dtausdtmp * pavlref(l)
       else
            k_wv_self = k_wv_self * pavl(l)
            dtausdtmp = dtausdtmp * pavl(l)
       end if

       ! kfix
       abso(0,l) = k_fix-(dk_wv + k_wv_self)*q(1,l)**2

       tautot(l)   = k_fix * wfix(l)
       dtaudtmp(l) = dkdt_fix * wfix(l)

       tautot(l)  = tautot(l) + (k_wv + (dk_wv + k_wv_self)*q(1,l))*w(1,l)
       abso(1,l)  = (2.-q(1,l))*(dk_wv + k_wv_self)*q(1,l) + k_wv

       dtaudtmp(l) = dtaudtmp(l)+(dkdt_wv+(ddkdt_wv+dtausdtmp)*q(1,l))*w(1,l)

       !---Variable gases
       DO kS=2,nmolS
          Imol=ImolS(kS)
          dft1_p1 = kvar(kS-1,indxp1,indxt1_p1) - kvar(kS-1,indxp1,indxt3_p1)
          dft2_p1 = kvar(kS-1,indxp1,indxt2_p1) - kvar(kS-1,indxp1,indxt3_p1)
          abs0_p1=at1_p1L*dft1_p1+at2_p1L*dft2_p1+kvar(kS-1,indxp1,indxt3_p1)
          dabs0_p1= adt1_p1L*dft1_p1 + adt2_p1L*dft2_p1

          if (.not.referenceGrid) then
             dft1_p2 = kvar(kS-1,indxp2,indxt1_p2) - kvar(kS-1,indxp2,indxt3_p2)
             dft2_p2 = kvar(kS-1,indxp2,indxt2_p2) - kvar(kS-1,indxp2,indxt3_p2)
             abs0_p2=at1_p2L*dft1_p2+at2_p2L*dft2_p2+kvar(kS-1,indxp2,indxt3_p2)
             dabs0_p2= adt1_p2L*dft1_p2 + adt2_p2L*dft2_p2
          end if
          abso(kS,l)   = (abs0_p1*ap1L + abs0_p2*ap2L)
          tautot(l)   = tautot(l) + abso(kS,l) * w(Imol,l)
          dtaudtmp(l) = dtaudtmp(l)+(dabs0_p1*ap1L+dabs0_p2*ap2L)*w(Imol,l)
       END DO
    END DO
    RETURN
  END SUBROUTINE OSStran

!-------------------------------------------------------
! performes 3 point intepolation (cubic) using
! precomputed offsets d1 and d2
! as well as derivatives dd1, dd2
!
  subroutine threePointInterpolation(d1, d2, dd1, dd2, z1, z2, z3, y, dy)
     real, intent(in)   :: d1, d2, dd1, dd2
     real, intent(in)  :: z1,z2,z3
     real, intent(out)  :: y, dy

     ! local
     real        :: df1, df2

       df1 = z1 - z3
       df2 = z2 - z3
       y  = d1*df1 + d2*df2 + z3
       dy = dd1*df1 + dd2*df2
       return
  end subroutine threePointInterpolation

!-------------------------------------------------------
! find indices for self continum grid
!
  SUBROUTINE setIndex_selfCont(tavl,N2,indxt,at1,at2)
    !---Input variables
    INTEGER                 ,INTENT(IN) :: n2
    real,    DIMENSION(:)   ,INTENT(IN) :: tavl
    !---Output variables
    INTEGER, DIMENSION(:),INTENT(INOUT) :: indxt
    real,    DIMENSION(:),INTENT(INOUT) :: at1,at2
    !---Local variables
    INTEGER               :: l,lt,lt1,lt2

    DO l=1,N2
       lt=1
       DO WHILE(tmpself(lt) < tavl(l))
          lt=lt+1
          IF (lt > nTself) exit
       END DO

       if (lt == 1) then
          lt1=1
          lt2=2
       else if (lt > nTself) then
          lt1=nTself-1
          lt2=nTself
       else
          lt1=lt-1
          lt2=lt
       end if
       at1(l)=(tmpself(lt2)-tavl(l))/(tmpself(lt2)-tmpself(lt1))
       at2(l)=1.-at1(l)
       indxt(l)=lt1
    END DO
    RETURN
  END SUBROUTINE setIndex_selfCont

  !----------------------------------------------------------------------------
  ! PURPOSE: Interpolation in wavenumber
  !----------------------------------------------------------------------------
  SUBROUTINE vinterp(datin,gridin,vn,datout,ip0,coefInt)
    !---Input variables
    real, DIMENSION(:,:), INTENT(IN)  :: datin
    real, DIMENSION(:)  , INTENT(IN)  :: gridin
    REAL                , INTENT(IN)  :: vn
    !---Input/output variables
    INTEGER           ,INTENT(INOUT)  :: ip0
    !---Output variables
    real, DIMENSION(:),INTENT(INOUT)  :: datout
    REAL              ,INTENT(INOUT)  :: coefInt
    !---Local variables
    INTEGER :: nxdim,ngin,ip,i

    nxdim=SIZE(datin,2)
    ngin=SIZE(gridin)
    if (vn <= gridin(1)) then
      datout(1:nxdim)=datin(1,1:nxdim)
      return
    else if (vn > gridin(ngin)) then
      datout(1:nxdim)=datin(ngin,1:nxdim)
      return
    end if

    do ip=ip0,ngin-1
       if(gridin(ip)>vn) EXIT
    end do
    ip0=ip
    coefInt=(vn-gridin(ip0-1))/(gridin(ip0)-gridin(ip0-1))
    datout(1:nxdim)=coefInt*datin(ip0,1:nxdim) + (1.-coefInt)*datin(ip0-1,1:nxdim)
    RETURN
  END SUBROUTINE vinterp

!----------------------------------------------------------------------------

  SUBROUTINE totalDensity(p,t,rhotot, RHOW, x)
    real, PARAMETER     :: ALOSMT=2.6867775E+19!  Loschmidt number [1/cm^3]
    real, PARAMETER     :: PZERO=1013.25       !  standard pressure [hPa]
    real, PARAMETER     :: TZERO=273.15        !  standard temperature [K]
    real,           INTENT(IN)    :: p,t
    real, OPTIONAL, INTENT(IN)    :: x
    real, INTENT(OUT)             :: RHOTOT    ! total density in [molecules/cm^3]
    real, OPTIONAL, INTENT(OUT)   :: RHOW      ! density of WV in [molecules/cm^3]
    !Local
    real, parameter               :: ratio = drymwt/molWt(1)

    RHOTOT=ALOSMT*(p/PZERO)*(TZERO/t)
    IF (present(x)) THEN
       RHOW  =RHOTOT*x*ratio/(1.+x*ratio)
    END IF
    RETURN
  END SUBROUTINE totalDensity

!----------------------------------------------------------------------------
! computer average temperature assuming exponential decay with altitude
!
  SUBROUTINE lzsum4T(xu,xl,xint,dxu,dxl)
    real,   INTENT(IN) :: xu,xl
    real,   INTENT(OUT):: xint,dxu,dxl
    REAL               :: eps
    eps=xl/xu-1.
    if (abs(eps) > 1E-2) then
       xint    =(xl-xu)/log(xl/xu)
       dxu     =xint*(1./(xu-xl)+1./xu/log(xl/xu))
       dxl     =xint*(1./(xl-xu)+1./xl/log(xu/xl))
    else
       xint    =0.5*(xl+xu)-xu*eps*eps/12.0
       dxu     =0.5+eps/6.0
       dxl     =0.5-eps/6.0
    end if
    RETURN
  END SUBROUTINE lzsum4T

!----------------------------------------------------------------------------
! computer molecular amount assuming exponential decay with altitude
!
  SUBROUTINE lzsum4W(xu,xl,dz,RHOTOT1,RHOTOT2,xint,dxu,dxl,xint_tot)
    real,   INTENT(IN) :: xu,xl,dz
    real,   INTENT(IN) :: RHOTOT1     !total density (molec/cm^3) at upper level
    real,   INTENT(IN) :: RHOTOT2     !total density (molec/cm^3) at lower level
    real, INTENT(INOUT):: xint        !column amount (molec/cm^2)
    real, INTENT(INOUT):: dxu,dxl
    real, OPTIONAL, INTENT(INOUT) :: xint_tot !dry column amount (g/cm^2)

    REAL               :: RHOAIR_u ! density on the upper level (molec/cm^3)
    REAL               :: RHOAIR_l ! density on the lower level (molec/cm^3)
    REAL               :: HZ       ! inverse exponential
    RHOAIR_u=RHOTOT1*xu
    RHOAIR_l=RHOTOT2*xl

    HZ      =DZ/log(RHOAIR_u/RHOAIR_l)
    xint    =HZ*(RHOAIR_u-RHOAIR_l)*1E+5
    dxu     =xint*( RHOTOT1/(RHOAIR_u-RHOAIR_l)-1./log(RHOAIR_u/RHOAIR_l)/xu)
    dxl     =xint*(-RHOTOT2/(RHOAIR_u-RHOAIR_l)+1./log(RHOAIR_u/RHOAIR_l)/xl)

    IF (present(xint_tot)) THEN
       HZ      =DZ/log(RHOTOT1/RHOTOT2) !inverse exponential for total
       xint_tot=(HZ*(RHOTOT1-RHOTOT2)*1E+5)
    END IF
    RETURN
  END SUBROUTINE lzsum4W

  !----------------------------------------------------------------------------
  ! PURPOSE: Computes average layer quantities (or integrated amount)
  !          using a linear dependence on P.
  !----------------------------------------------------------------------------
  SUBROUTINE lpsum(pu,pl,xu,xl,scal,xint,dxu,dxl)
    !---Input variables
    real, INTENT(IN)  :: pu,pl,xu,xl,scal
    !---Output variables
    real, INTENT(OUT) :: xint,dxu,dxl

    xint = 0.5*(pl-pu)*scal
    dxu      = xint
    dxl      = xint
    xint     = xint*(xu+xl)
    RETURN
  END SUBROUTINE lpsum

  !----------------------------------------------------------------------------
  ! PURPOSE: Computes average layer quantities (or integrated amount)
  !          using a log-x dependence on log-p.
  !----------------------------------------------------------------------------
  SUBROUTINE lpsum_log(pu,pl,xu,xl,scal,xint,dxu,dxl)
    real, PARAMETER :: epsiln=1E-12
    !---Input variables
    real, INTENT(IN)  :: pu,pl,xu,xl,scal
    !---Output variables
    real, INTENT(OUT) :: xint,dxu,dxl
    !---Local variables
    REAL              :: hp,x0,zeta,alza,alpha
    hp       = log(pl/pu)
    x0       = pl*xl*hp
    zeta     = pu*xu/(pl*xl)
    IF(ABS(zeta-1.) > epsiln)THEN
       alza  = log(zeta)
       xint  = x0*(zeta-1.)/alza
       alpha = zeta/(zeta-1.)-1./alza
    ELSE
       xint  = x0*(1.+(zeta-1.)/2.)
       alpha = 0.5+(zeta-1.)/12.
    END IF
    xint     = xint*scal
    dxu      = xint*alpha/xu
    dxl      = xint*(1.-alpha)/xl
    RETURN
  END SUBROUTINE lpsum_log

  !----------------------------------------------------------------------------
  ! PURPOSE: Interpolation in linear pressure
  !----------------------------------------------------------------------------
  SUBROUTINE lint(xInp,pGrid,N0,p0,x0)
    INTEGER           , INTENT(IN)   :: n0
    REAL              , INTENT(IN)   :: p0
    real, DIMENSION(:), INTENT(IN)   :: xInp,pGrid
    !---Output variable
    REAL              , INTENT(INOUT):: x0
    !---Local variable
    REAL               :: xx
    xx = (xInp(N0)-xInp(N0-1))/(pGrid(N0)-pGrid(N0-1))
    x0 = xInp(N0-1)+(p0-pGrid(N0-1))*xx
    RETURN
  END SUBROUTINE lint

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


  !----------------------------------------------------------------------------
  ! PURPOSE: Calculates Planck function and its derivative  wrt temperature
  !----------------------------------------------------------------------------
  SUBROUTINE planck(f1,f2,t,rad,draddt)
    !---Input variables
    REAL     ,INTENT(IN)      :: t,f1,f2
    !---Output variables
    REAL  ,INTENT(INOUT)      :: rad,draddt
    !---Local variables
    REAL            :: f3
    f3=EXP(-f2/t)
    rad=f1/(1. - f3)
    draddt=(f2*f3/f1)*(rad/t)**2
    rad=rad*f3
    RETURN
  END SUBROUTINE planck

  !----------------------------------------------------------------------------
  ! PURPOSE: For linear-In-Tau approach, compute the Planck radiance weight
  !          (hlf) of the far boundary of the layer, and its derivative (dhlf)
  !----------------------------------------------------------------------------
  subroutine rlin(od,hlf,dhlf)
    !Input variables
    real,    INTENT(IN) :: od
    !Output variables
    real, INTENT(OUT) :: hlf,dhlf
    !Local variables
    real    tex
    if (od<pMatchLim) then
       hlf=.5+od*dbAvg2dtau
       dhlf=dbAvg2dtau
    else
       tex=exp(-od)
       hlf=1./od-tex/(1. - tex)
       dhlf=-1./od**2 + tex/(1.-tex)**2
    end if
    return
  end subroutine rlin

  !----------------------------------------------------------------------------
  ! PURPOSE: For replacing default mixing ratios of those molecules
  !          transferred from variable to fixed
  !----------------------------------------------------------------------------
  SUBROUTINE amnt2change(iuF,molID,vmol)
     !Input variables
     INTEGER                 ,INTENT(IN) :: iuF
     integer                 ,INTENT(IN) :: molID
     !Output variables
     real,    DIMENSION(:),INTENT(INOUT) :: vmol
     !Local variables
     INTEGER                :: k
     READ(iuF,*)k
     IF (molID == 1 .AND. k /= 1) THEN
        print*, 'Err[oss_ir_module::GetOD]: Include H2O mixing ratio into defProfFile'
        call exit(1)
     END IF
     READ(iuF,*)vmol(1:nLayOD+1)
     IF (k == molID) then
         RETURN
     end if

     DO WHILE (k /= molID)
        READ(iuF,*)k
        READ(iuF,*)vmol(1:nLayOD+1)
     end do
     RETURN
  END SUBROUTINE amnt2change

 !----------------------------------------------------------------------------
 ! PURPOSE: to add a channel subset for further use
 ! input
 !   iChList - subset indices (has to be the same format as in selection file)
 !   nChList - length of the subset
 !   optional input
 !   chSetIDin - user supplied ID to the subset
 !   optional output
 !   chSetIDout - system assigned ID to the subset
 !   in case chSetIDin is present then chSetIDout=chSetIDin
 !----------------------------------------------------------------------------
 subroutine addChanSelect(iChList, nChList, chSetIDin, chSetIDout)
   integer, dimension(:), intent(in)    :: iChList
   integer,               intent(in)    :: nChList
   integer,optional,      intent(in)    :: chSetIDin
   integer,optional,      intent(out)   :: chSetIDout

   integer                :: k,n1,ich,ismp,ismpSel, idx, lenBase, nCh
   integer                :: curChSet, chSetIDinLoc, chSetIDOutLoc

   character(len=*), parameter :: procName=modName//'::addChanSelect'
   character(len=*), parameter :: errHeader='Err['//procName//']: '
   character(len=*), parameter :: wrnHeader='Wrn['//procName//']: '

   if (lastChSet >= 0) then
      if (present(chSetIDin)) then
        if (chSetIDin <= 0) then
          print*, errHeader,' chSetIDin must be a positive integer.'
          call exit(1)
        end if
        if (present(chSetIDout)) chSetIDout=chSetIDin
        curChSet = userIndex2chSetIndex(chSetIDin)
        if (curChSet == UNDEFINED_INDEX) then
           if (lastChSet >= mxIndex) then
              print*, errHeader,' maximum number of channel subset has been reached.'
              call exit(1)
           end if
           lastChSet = lastChSet+1
           curChSet = lastChSet
        end if
      else
        lastChSet = lastChSet+1
        curChSet = lastChSet
        if (present(chSetIDout)) chSetIDout=curChSet
      end if

      lenBase = nChList_arr(0)
      do k=1,nChList
         if ( .not. ANY(iChList(k) == chanList_arr(1:lenBase,0))) THEN
            print*, errHeader,' In subset ',  lastChSet, &
               ' some channel indices do not present in the SEL file '
            call exit(1)
         end if
      end do
      do n1=1,nChList-1
         if ( ANY(iChList(n1+1:nChList) == iChList(n1))) THEN
            print*, errHeader,' There are duplicates in selected channels '
            call exit(1)
         end if
      end do
   else
      lastChSet = 0
      curChSet = 0
   end if

   if (present(chSetIDin)) then
     userIndex2chSetIndexMap(curChSet) = chSetIDin
   else
     userIndex2chSetIndexMap(curChSet) = curChSet
   end if

   ismpSel=0

   if (nChList == nchanAll) then
      nNodes_arr(curChSet)=nf_sel
      iselS_arr(1:nf_sel,curChSet)=isel_Loc(1:nf_sel)
      nch_arr(1:nf_sel,curChSet)=nch_Loc(1:nf_sel)
      do ismp=1,nf_sel
         nCh = nch_Loc(ismp)
         ichMap_arr(1:nCh,ismp,curChSet)= ichMap_Loc(1:nCh,ismp)
         coef_arr(1:nCh,ismp,curChSet)= coef_Loc(1:nCh,ismp)
      end do
   else
      do ismp=1,nf_sel
         ismpSel=ismpSel+1
         ich=0
         do k=1,nch_Loc(ismp)
            do idx = 1, nChList
               if (iChList(idx) == ichmap_Loc(k,ismp))  THEN
                   ich=ich+1
                   ichMap_arr(ich,ismpSel,curChSet)=idx
                   coef_arr(ich,ismpSel,curChSet)=coef_Loc(k,ismp)
                   EXIT
                end if
            end do
         end do
         if (ich == 0) THEN
            ismpSel=ismpSel-1
         else
            nch_arr(ismpSel,curChSet)=ich
            iselS_arr(ismpSel,curChSet)=isel_Loc(ismp)
         end if
         nNodes_arr(curChSet)=ismpSel
      end do
   end if

   nChList_arr(curChSet)=nChList
   chanList_arr(1:nChList,curChSet)=iChList(1:nChList)
   return
  end subroutine addChanSelect
 !----------------------------------------------------------------------------
 ! PURPOSE: to load channel subset for further use
 ! input
 !   iChList - mask (the same size as full set
 !   nChList - length of the mask
 !   optional input
 !   chSetIDin - user supplied ID to the subset
 !   optional output
 !   chSetIDout - system assigned ID to the subset
 !   in case chSetIDin is present then chSetIDout=chSetIDin
 !----------------------------------------------------------------------------
 subroutine addChanSelectMask(iChMask, nChMask, chSetIDin, chSetIDout)
    integer, dimension(:), intent(in)    :: iChMask
    integer,               intent(in)    :: nChMask
    integer,optional,      intent(in)    :: chSetIDin
    integer,optional,      intent(out)   :: chSetIDout

    integer                :: k,n1,ich,ismp,ismpSel, idx, lenBase
    integer                :: jj
    integer                :: curChSet
    integer                :: nChList
    integer, allocatable, dimension(:)   :: iChanIndex
    character(len=*), parameter :: procName=modName//'::addChanSelectMask'
    character(len=*), parameter :: errHeader='Err['//procName//']: '
    character(len=*), parameter :: wrnHeader='Wrn['//procName//']: '

    lenBase = nChList_arr(0)
    if ( nChMask /= lenBase) THEN
       print*, errHeader,' Supplied mask has to be of size all channel set', &
          nChMask, lenBase
       call exit(1)
    end if

    nChList=count(iChMask(1:nChMask) > 0)
    allocate(iChanIndex(nChList))
    iChanIndex = pack(chanList_arr(1:lenBase, 0),iChMask /= 0)
    call addChanSelect(iChanIndex, nChList, chSetIDin, chSetIDout)
    deallocate(iChanIndex)
    return
   end subroutine addChanSelectMask
  !------------------------------------------------------------------------
  SUBROUTINE setChanSelect(userChSet)
   !Input
    INTEGER              ,INTENT(IN) :: userChSet
   !local variable
    INTEGER                          :: nCh,k
    integer                          :: selectedChSet
    character(len=*), parameter :: procName=modName//'setChanSelect'
    character(len=*), parameter :: errHeader='Err['//procName//']: '
    character(len=*), parameter :: wrnHeader='Wrn['//procName//']: '

    selectedChSet=userIndex2chSetIndex(userChSet)

    IF (selectedChSet == iChSet) RETURN

    IF (lastChSet < 0) THEN
        print*, wrnHeader, 'Channel selections have not been loaded'&
              ,' channels selection is used '
        return
    END IF
    if (selectedChSet == UNDEFINED_INDEX ) then
        iChSet = UNDEFINED_INDEX
        return
    elseif (selectedChSet < 0 .OR. selectedChSet > lastChSet) THEN
        print*, wrnHeader, 'Invalid channel set selection.'&
              , ' Active selection index ', iChSet
        return
    END IF

    iChSet=selectedChSet

    nNodes=nNodes_arr(iChSet)
    iselS(1:nNodes)=iselS_arr(1:nNodes,iChSet)
    nch_ir(1:nNodes)=nch_arr(1:nNodes,iChSet)
    DO k=1,nNodes
      nCh=nch_ir(k)
      ichMap_ir(1:nCh,k)=ichMap_arr(1:nCh,k,iChSet)
      coef_ir(1:nCh,k)  =coef_arr(1:nCh,k,iChSet)
    END DO
    RETURN
  END SUBROUTINE setChanSelect

  !----------------------------------------------------------------------------
  ! PURPOSE: Create auxiliary array, mapping molID onto a position within
  !          array of absorption coefficients
  !----------------------------------------------------------------------------
  SUBROUTINE invertMolID(nmol,molID,invMolID)
    !Input variables
    INTEGER,    INTENT(IN) :: nmol
    INTEGER,   INTENT(IN) :: molID(nmol)
    !Output variables
    INTEGER,    INTENT(INOUT) :: invMolID(molID(nmol))
    !Local variables
    INTEGER k
    invMolID(:)=0

    DO k=1,nmol
       invMolID(molID(k))=k
    end do
    RETURN
  END SUBROUTINE invertMolID

  !----------------------------------------------------------------------------
  FUNCTION getMolecMass(mol)
    INTEGER, intent(in)   :: mol          ! ID of molecule
    REAL                  :: getMolecMass ! molecular mass (g/mol)
    getMolecMass=molWt(mol)
  END FUNCTION getMolecMass

  !--------------------------------------------------------------------------
  ! Helper function
  ! getLastSelectionIndex
  ! return the last channel selection index loaded from the selection file

  FUNCTION getLastSelectionIndex()
    integer                     ::  getLastSelectionIndex
    getLastSelectionIndex=lastChSet
    return
  END FUNCTION getLastSelectionIndex

  !--------------------------------------------------------------------------
  ! Helper function
  ! getLastSelectionIndex
  ! return the current channel selection index

  FUNCTION getSelectionIndex()
    integer                     ::  getSelectionIndex
    getSelectionIndex=iChSet
    return
  END FUNCTION getSelectionIndex

  !--------------------------------------------------------------------------
  ! Helper function
  ! getLastSelectionIndex
  ! convert user assigned selection ID to the local
  function userIndex2chSetIndex(selectionIndex)
    integer, intent(in) :: selectionIndex
    integer             :: userIndex2chSetIndex
    integer             :: kk

    userIndex2chSetIndex = UNDEFINED_INDEX

    if ( selectionIndex == UNDEFINED_INDEX) return

    if (selectionIndex == 0) then
       userIndex2chSetIndex = 0
       return
    else
       do kk=1, lastChSet
          if (userIndex2chSetIndexMap(kk) == selectionIndex) then
             userIndex2chSetIndex = kk
             return
          end if
       end do
    end if
    return
  end function userIndex2chSetIndex


  !--------------------------------------------------------------------------
  ! Helper function
  ! getCountChannel
  ! return channel count

  FUNCTION getCountChannel(selectionIndex)
    integer                       ::  getCountChannel
    integer, optional, intent(in) :: selectionIndex
    if (present(selectionIndex)) then
        if (checkSelectionIndex(selectionIndex) .and. &
            selectionIndex /= UNDEFINED_INDEX) then
            getCountChannel=nChList_arr(selectionIndex)
        else
            print*, 'Wrn[oss_ir_module:: getCountChannel]: ', &
                    'Invalid channel set selection.'
            getCountChannel=nChList_arr(iChSet)
        end if
    else
        getCountChannel=nChList_arr(iChSet)
    end if

    return
  END FUNCTION getCountChannel

  !--------------------------------------------------------------------------
  ! Helper function
  ! getCountUsedNode
  ! return used nodes count

  FUNCTION getCountUsedNode(selectionIndex)
    integer               ::  getCountUsedNode
    integer, optional, intent(in) :: selectionIndex
    ! local variables
    integer :: currentIndex
    currentIndex = iChSet

    if (present(selectionIndex)) then
        if (checkSelectionIndex(selectionIndex) .and. selectionIndex /= UNDEFINED_INDEX) then
            currentIndex = selectionIndex
        else
            print*, 'Wrn[oss_ir_module:: getCountUsedNode]: Invalid channel set selection: ', &
                 selectionIndex
            call exit(1)
        end if
    end if
    getCountUsedNode=nNodes
    return
  END FUNCTION getCountUsedNode

  !--------------------------------------------------------------------------
  ! Helper function
  ! getChannelFrequency
  ! return selected channel frequency

  SUBROUTINE getChannelFrequency(channelFrequency, selectionIndex)
    real, dimension(:), intent(out) :: channelFrequency
    integer, optional, intent(in) :: selectionIndex
    !local
    integer :: nchan_out
    integer :: currentIndex

    currentIndex = iChSet
    if (present(selectionIndex)) then
        if (checkSelectionIndex(selectionIndex) .and. selectionIndex .ne. UNDEFINED_INDEX) then
            currentIndex = selectionIndex
        else
            print*, 'Wrn[oss_ir_module:: getChannelFrequency]: Invalid channel set selection: ', &
                selectionIndex
               call exit(1)
            END IF
        end if
    nchan_out = nChList_arr(currentIndex)
        IF (SIZE(channelFrequency) < nchan_out) THEN
           print*, 'Err[oss_ir_module::getChannelFrequency]: channelFrequency vector too small'
           call exit(1)
        END IF
    channelFrequency(1:nchan_out) = cWvn(chanList_arr(1:nchan_out, currentIndex))
    return
  END SUBROUTINE getChannelFrequency
  !--------------------------------------------------------------------------
  ! Helper function
  ! getChannelIndex
  ! return selected channel indexes

  FUNCTION getCurrentSelection()
    integer getCurrentSelection
    getCurrentSelection = iChSet
  END FUNCTION getCurrentSelection
  !--------------------------------------------------------------------------
  ! Helper function
  ! getChannelIndex
  ! return selected channel indexes

  SUBROUTINE getChannelIndex(channelIndex, selectionIndex)
    integer, dimension(:), intent(out) :: channelIndex
    integer, optional, intent(in) :: selectionIndex

    !local
    integer :: nchan_out
    integer :: currentIndex
    currentIndex = iChSet

    if (present(selectionIndex)) then
        if (checkSelectionIndex(selectionIndex) &
            .and. selectionIndex /= UNDEFINED_INDEX) then
            nchan_out = nChList_arr(selectionIndex)
            currentIndex = selectionIndex
        else
            print*, 'Wrn[oss_ir_module:: getChannelIndex]: ', &
             'Invalid channel set selection: ',  selectionIndex
            call exit(1)
        end if
    end if

    nchan_out = nChList_arr(currentIndex)
    IF (SIZE(channelIndex) < nchan_out) THEN
       print*, 'Err[oss_ir_module::getChannelIndex]: channelIndex vector too small'
       call exit(1)
    END IF
    channelIndex(1:nchan_out) = chanList_arr(1:nchan_out, currentIndex)
    return
  END SUBROUTINE getChannelIndex

  !--------------------------------------------------------------------------
  ! Helper function
  ! checkSelectionIndex
  ! return  true if checkSelectionIndex is of valid value

  FUNCTION checkSelectionIndex(selectionIndex)
    logical                    ::checkSelectionIndex
    integer, intent(in)        :: selectionIndex
    if (selectionIndex == UNDEFINED_INDEX ) then
        checkSelectionIndex = .TRUE.
        return
    elseif ((selectionIndex >= 0) .AND. (selectionIndex <= lastChSet)) THEN
        checkSelectionIndex = .TRUE.
        return
    END IF
    checkSelectionIndex = .FALSE.
    return
  END FUNCTION checkSelectionIndex

  !--------------------------------------------------------------------------
  ! Helper function
  ! getUsedNodes
  ! return  used node frequency

  SUBROUTINE getUsedNodes(nodeFrequency)
    real, dimension(:), intent(out) :: nodeFrequency

    IF (SIZE(nodeFrequency) < nNodes) THEN
       print*, 'Err[oss_ir_module::getUsedNodes]: nodeFrequency vector too small'
       call exit(1)
    END IF

    nodeFrequency(1:nNodes) = vWvn(iselS(1:nNodes))
    return
  END SUBROUTINE getUsedNodes

  !--------------------------------------------------------------------------
  ! Helper function
  ! getReferenceGrid
  ! return reference pressure grid from LUT file

  SUBROUTINE getReferencePressure(pressure)
    real, allocatable, dimension(:), intent(inout) :: pressure
    if (allocated(pressure)) then
        if (size(pressure, dim = 1) == numRefLev) then
            pressure = pRef
            return
        end if
    end if
    allocate(pressure(numRefLev))
    pressure = pRef
    return
  END SUBROUTINE getReferencePressure

  !--------------------------------------------------------------------------
  ! Helper function
  ! FindFreeUnit
  ! return a free UNIT
  FUNCTION findFreeUnit()
    integer :: findFreeUnit
    ! local
    logical :: isOpened
    integer :: iostat
    DO findFreeUnit = 10, 200
            inquire(UNIT = FindFreeUnit, OPENED = isOpened, iostat = iostat)
            if (iostat .ne. 0 ) cycle
            if (.not. isOpened) return
    END DO
    STOP '[FindFreeUnit] No Free File Unit Found'
  END FUNCTION findFreeUnit

  !--------------------------------------------------------------------------
  ! getConstant
  ! return constants used on the module
  ! input name
  ! AVOGADRO Avogadro number
  ! RCONST gas constant
  ! RAIRDFLT - default rair
  !
  FUNCTION getConstant(name)
    real                         :: getConstant
    character(len=*), intent(in) :: name
    if (name(1:2) =='AVOGADRO') THEN
        getConstant = AVOGAD
    ELSE IF (name(1:2) =='RCONST') THEN
        getConstant = rConst
    ELSE IF (name(1:2) =='RGAS') THEN
        getConstant = Rgas
    ELSE IF (name(1:2) =='RAIRDFLT') THEN
        getConstant = rairDflt
    ELSE
        getConstant = - 999999.99
        print *, 'The name ', trim(name), ' is unknown'
    END IF
  END FUNCTION getConstant

  !--------------------------------------------------------------------------
  SUBROUTINE OSSdestroyIRsubs()
    if (allocated(nch_ir))        deallocate(nch_ir)
    if (allocated(vWvn))          deallocate(vWvn)
    if (allocated(f1_arr))        deallocate(f1_arr)
    if (allocated(f2_arr))        deallocate(f2_arr)
    if (allocated(coef_ir))       deallocate(coef_ir)
    if (allocated(cWvn))          deallocate(cWvn)
    if (allocated(pRef))          deallocate(pRef)
    if (allocated(pavlref))       deallocate(pavlref)
    if (allocated(sunrad))        deallocate(sunrad)
    if (allocated(ichmap_ir))     deallocate(ichmap_ir)
    if (allocated(NmolS))         deallocate(NmolS)
    if (allocated(ImolS))         deallocate(ImolS)
    if (allocated(Tmptab))        deallocate(Tmptab)

    if (allocated(kfix_ir))       deallocate(kfix_ir)
    if (allocated(kh2o_ir))       deallocate(kh2o_ir)
    if (allocated(dkh2o_ir))      deallocate(dkh2o_ir)
    if (allocated(kvar_ir))       deallocate(kvar_ir)
    if (allocated(kself))         deallocate(kself)

    if (allocated(nChList_arr))   deallocate(nChList_arr)
    if (allocated(chanList_arr))  deallocate(chanList_arr)
    if (allocated(iselS))         deallocate(iselS)
    if (allocated(nNodes_arr))    deallocate(nNodes_arr)
    if (allocated(iselS_arr))     deallocate(iselS_arr)
    if (allocated(isel_Loc))      deallocate(isel_Loc)
    if (allocated(nch_Loc))       deallocate(nch_Loc)
    if (allocated(nch_arr))       deallocate(nch_arr)
    if (allocated(coef_arr))      deallocate(coef_arr)
    if (allocated(coef_Loc))      deallocate(coef_Loc)
    if (allocated(ichMap_arr))    deallocate(ichMap_arr)
    if (allocated(ichmap_Loc))    deallocate(ichmap_Loc)

    if (allocated(TmpSelf))       deallocate(TmpSelf)

    getODdone=.false.
    myLinInTau=.FALSE.
    iChSet=UNDEFINED_INDEX
    sphericalGeometryFlag = .false.
    zIntegrationFlag = .false.
    call destroySolar()

  END SUBROUTINE OSSdestroyIRsubs

  !--------------------------------------------------------------------------
  !
  ! Set flag sphericalGeometryFlag to a desired state
  !
  ! depending on sphericalGeometryFlag all path calculations are performed assuming:
  !  plane parallel atmosphere (.FALSE.)
  !  spherical atmosphere (no refraction)  (.true.)

  SUBROUTINE setPlaneParallelFlag(flag)
      logical, intent(in) :: flag
      sphericalGeometryFlag = .not. flag
  END SUBROUTINE setPlaneParallelFlag

END MODULE OSSIRmoduleSubs
