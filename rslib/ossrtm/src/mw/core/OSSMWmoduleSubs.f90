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
MODULE OSSMWmoduleSubs
  IMPLICIT NONE
  PUBLIC
  !------------------------------------------------------------------------
  ! Public items made available to calling program(s)
  !------------------------------------------------------------------------
  PUBLIC :: osstc_mw,ABCOEFLIQ, &
            getOSSselMW,loadOSSselMW, &
            set2ndOrderPlanck,getPlanckLinCtr


  !----TEST: UNCOMMENT THIS IF YOU WANT TO PASS DATA BY USE ASSOCIATION
  !          INSTEAD OF ARGUMENTS TO OSSINIT
  !
  !PUBLIC :: NparG,ITempG,ITskinG,IPsfcG,ImolG,ICldLiqG
  !PUBLIC :: nchan,my_cFreq,nlev,Pref
  !--------

  !------------------------------------------------------------------------
  ! Parameters/Dimensions
  !------------------------------------------------------------------------
  INTEGER,PARAMETER,PUBLIC::OSS_LEN_ID=12       ! length of channel ID strings
  INTEGER,PARAMETER,PUBLIC::OSS_MAX_HDR_OPT=40  ! max optional header entries
  CHARACTER(len=OSS_LEN_ID),PARAMETER,PUBLIC::OSS_NULL_ID='NULL'
  INTEGER,PARAMETER::INACTIVE_UID=-1

  INTEGER,PARAMETER::mxlev=120,mxlay=mxlev-1,mxchan=80 ! set 142 for hymn test
  INTEGER,PARAMETER::mxHmol=2,mxmolS=2
  INTEGER,PARAMETER::mxparg=(mxHmol+1)*mxlev+mxchan+10
  INTEGER,PARAMETER::DOUBLE=SELECTED_REAL_KIND(14)
  !REAL,   PARAMETER::rair=1.020408163 !=10./g converts p(mb) to g/cm**2
  REAL(DOUBLE),PARAMETER::GHzToHz=1.e9     ! Units conversion
  !------------------------------------------------------------------------------
  !
  ! physical constant
  !
  ! REAL(DOUBLE),PARAMETER::Plnk=6.626176d-27,LtSp=2.997925d+10,Bltz=1.380662d-16
  !     constants were updated based on
  !     LBLRTM   $Revision: 27087 $ phys_consts.f90
  !     Constants from NIST May 2010
  REAL(DOUBLE),PARAMETER::Plnk=6.62606876D-27
  REAL(DOUBLE),PARAMETER::LtSp=2.99792458d+10
  REAL(DOUBLE),PARAMETER::Bltz=1.3806503d-16


  REAL, PARAMETER       :: deg2rad = acos(-1.0) / 180.0
  real,    parameter :: grav_const_req1= 9.80665 ! m / s^2
  real,    parameter :: grav_const_req2=-0.02586 ! m / s^2
  real,    parameter               :: rEarth = 6371.23 ! km
  real,    PARAMETER :: pMatchLim = 1.E-3 ! relative difference use to match
                                          ! pressure levels

  !------------------------------------------------------------------------
  ! Mode control flags
  !------------------------------------------------------------------------
  LOGICAL :: use2ndOrderPlanck=.FALSE. ! Default is 2-point linear fit
  !------------------------------------------------------------------------
  ! Geophysical Vector Variables/Pointers
  !------------------------------------------------------------------------
  INTEGER                   :: numRefLev
  INTEGER                   :: NparG,ITempG,ITskinG,IPsfcG,ICldLiqG
  INTEGER                   :: nmolfix,Nmol
  INTEGER, DIMENSION(:),   ALLOCATABLE :: molID,molIDfix
  INTEGER, DIMENSION(:),   ALLOCATABLE :: ImolG

  !------------------------------------------------------------------------
  ! Instrument parameters
  !------------------------------------------------------------------------
  REAL,    DIMENSION(:),   ALLOCATABLE :: pavlref
  INTEGER                              :: nchan,nfsmp,uid_sel
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: ichmap
  REAL,    DIMENSION(:,:), ALLOCATABLE :: coef,coef1,coef2
  INTEGER, DIMENSION(:),   ALLOCATABLE :: my_ipol,nch
  REAL,    DIMENSION(:),   ALLOCATABLE :: my_cFreq,vFreq,cosmbk
  REAL,    DIMENSION(:),   ALLOCATABLE :: my_side1,my_side2,my_hwhmx
  CHARACTER(len=OSS_LEN_ID), &
           DIMENSION(:),   ALLOCATABLE :: my_chanID
  INTEGER                              :: myRadOrTb=0 ! Mode switch: 0=Tb, 1=Rad
                                                      ! Tb assumed for old files

  !------------------------------------------------------------------------
  ! OD tables
  !------------------------------------------------------------------------
  INTEGER                               :: NLayOD,NTmpOD,NWvpOD
  INTEGER, DIMENSION(:),   ALLOCATABLE  :: NmolS
  INTEGER, DIMENSION(:,:), ALLOCATABLE  :: ImolS
  REAL,    DIMENSION(:),   ALLOCATABLE  :: fixtab,Pref
  REAL,    DIMENSION(:,:), ALLOCATABLE  :: varDflt
  REAL,    DIMENSION(:,:), ALLOCATABLE  :: Tmptab,Wvptab,kfix,dkfix,kh2o,dkh2o
  REAL,    DIMENSION(:,:), ALLOCATABLE  :: kvar,dkvar

  ! loadDone implementation provides the option to call load*() in advance
  logical                               :: loadOSSselMWdone=.false.
  logical                               :: loadODdone=.false.
  logical, save :: linInTauFlag          = .false.
  logical, save :: sphericalGeometryFlag = .false.
  logical, save :: zIntegrationFlag      = .false.

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

   real,    parameter :: drymwt=28.964 ! Dry air molecular mass (g / mol)
   real,    parameter :: AVOGAD = 6.02214199E+23 ! Avogadro ( 1 / mol )
   integer, parameter :: MAXSPC=84
   real, dimension(MAXSPC), parameter  :: molWt=(/ &
                18.015, 44.010, 47.998,44.010, 28.011, 16.043, 31.999, &   ! 7
                30.010, 64.060, 46.010,17.030, 63.010, 17.000, 20.010, &   ! 14
                36.460, 80.920,127.910,51.450, 60.080, 30.030, 52.460, &   ! 21
                28.014, 27.030, 50.490,34.010, 26.030, 30.070, 34.000, &   ! 28
                66.010,146.050, 34.080,46.030,  0.   ,  0.   ,  0.   , &   ! 35
                 0.   ,  0.   , 28.053,32.042,  0.   ,  0.   ,  0.   , &   ! 42
                 0.   ,  0.   , 0.0   ,0.0   ,  0.   ,  0.   ,  0.   , &   ! 49
                 0.   ,153.820, 88.000,97.460,137.37,187.380,120.910 , &   ! 56
                 0.   , 86.470,108.010,  0.  , 68.12,121.05  ,  0.   , &   ! 63
                 0.   ,  0.   , 0.0   ,0.0   ,  0.   ,  0.   ,  0.   , &   ! 70
                 0.   ,  0.   , 0.0   ,0.0   ,  0.   ,  0.   ,  0.   , &   ! 77
                 0.   ,  0.   , 0.0   ,19.02 ,  0.   ,  0.   ,  0.   /)    ! 84
   !  molecular weight normalized on dry air weight
   real, dimension(MAXSPC), parameter  :: normMolWt= molWt/drymwt

CONTAINS

 !----------------------------------------------------------------------------
  !    PURPOSE: Reads the file (if needed) and sets the channel definition data
  !    and OSS selected frequencies and weights
  !----------------------------------------------------------------------------
  SUBROUTINE getOSSselMW(selFile,nchanOSS,cFreq,ipol,side1,side2,hwhmx, &
    chanID,iRadOrTb)
    !---Arguments
    CHARACTER(len=*),      INTENT(IN)              :: selFile
    INTEGER,               INTENT(INOUT), OPTIONAL :: nchanOSS
    REAL,    DIMENSION(:), INTENT(INOUT), OPTIONAL :: cFreq
    INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: ipol
    REAL,    DIMENSION(:), INTENT(INOUT), OPTIONAL :: side1
    REAL,    DIMENSION(:), INTENT(INOUT), OPTIONAL :: side2
    REAL,    DIMENSION(:), INTENT(INOUT), OPTIONAL :: hwhmx
    CHARACTER(len=OSS_LEN_ID), &
             DIMENSION(:), INTENT(INOUT), OPTIONAL :: chanID
    INTEGER,               INTENT(INOUT), OPTIONAL :: iRadOrTb

    !------------------------------------------------------------------------
    ! Set global variables from file variables
    !------------------------------------------------------------------------
    call loadOSSselMW(selFile)

    !---Load optional output arguments
    IF (PRESENT(nchanOSS)) nchanOSS=nchan
    IF (PRESENT(cFreq)) THEN
       IF (SIZE(cFreq)  < nchan) THEN
          PRINT *,'err[oss_mw_module::getOSSselMW]: ', &
          'SIZE(cFreq)  < nchan: ',SIZE(cFreq),nchan
          STOP
       ELSE
          cFreq(1:nchan)=my_cFreq(1:nchan)
       END IF
    END IF
    IF (PRESENT(ipol)) THEN
       IF (SIZE(ipol)  < nchan) THEN
          PRINT *,'err[oss_mw_module::getOSSselMW]: ', &
          'SIZE(ipol)  < nchan: ',SIZE(ipol),nchan
          STOP
       ELSE
          ipol(1:nchan)=my_ipol(1:nchan)
       END IF
    END IF
    IF (PRESENT(side1)) THEN
       IF (SIZE(side1)  < nchan) THEN
          PRINT *,'err[oss_mw_module::getOSSselMW]: ', &
          'SIZE(side1)  < nchan: ',SIZE(side1),nchan
          STOP
       ELSE
          side1(1:nchan)=my_side1(1:nchan)
       END IF
    END IF
    IF (PRESENT(side2)) THEN
       IF (SIZE(side2)  < nchan) THEN
          PRINT *,'err[oss_mw_module::getOSSselMW]: ', &
          'SIZE(side2)  < nchan: ',SIZE(side2),nchan
          STOP
       ELSE
          side2(1:nchan)=my_side2(1:nchan)
       END IF
    END IF
    IF (PRESENT(hwhmx)) THEN
       IF (SIZE(hwhmx)  < nchan) THEN
          PRINT *,'err[oss_mw_module::getOSSselMW]: ', &
          'SIZE(hwhmx)  < nchan: ',SIZE(hwhmx),nchan
          STOP
       ELSE
          hwhmx(1:nchan)=my_hwhmx(1:nchan)
       END IF
    END IF
    IF (PRESENT(chanID)) THEN
       IF (SIZE(chanID)  < nchan) THEN
          PRINT *,'err[oss_mw_module::getOSSselMW]: ', &
          'SIZE(chanID)  < nchan: ',SIZE(chanID),nchan
          STOP
       ELSE
          chanID(1:nchan)=my_chanID(1:nchan)
       END IF
    END IF
    IF (PRESENT(iRadOrTb)) iRadOrTb=myRadOrTb

    RETURN
  END SUBROUTINE getOSSselMW

  !----------------------------------------------------------------------------
  !    PURPOSE: Reads the file that contains channel definition data and
  !    OSS selected frequencies and weights
  !----------------------------------------------------------------------------
  SUBROUTINE loadOSSselMW(selFile)
    !---Arguments
    CHARACTER(len=*),      INTENT(IN)              :: selFile
    !---Local variables
    CHARACTER*100     :: instr_info
    INTEGER           :: nf_sel,nchmax
    INTEGER           :: ismp,isels,j,istat
    INTEGER           :: nHdrOpt
    INTEGER           :: fmtver  ! file format version number
    CHARACTER(len=OSS_MAX_HDR_OPT)                       :: HdrOptFlag
    INTEGER           :: ius
    real*4,  dimension(:), allocatable :: realBuffer

    if(loadOSSselMWdone) return

    !---Open File
    ius = findFreeUnit()
    OPEN(ius,file=selFile,form='unformatted',status='old',err=800)
    !---Read the header of selection file
    READ(ius)uid_sel
    READ(ius,iostat=istat)fmtver,HdrOptFlag
    if (istat /= 0) then   ! for back compatibility with old format
       fmtver = 1
       HdrOptFlag = 'CPDEB'
       uid_sel=INACTIVE_UID
    endif
    if (fmtver /= 1) then
       PRINT *,'err[oss_mw_module::loadOSSselMW]:'
       PRINT *,'Format version of OSS sel file not recognized: ',fmtver
       STOP
    endif
    READ(ius)instr_info
    READ(ius)nchan,nf_sel,nchmax
    nfsmp=nf_sel
    if (.not. allocated(my_cFreq))      allocate(my_cFreq(nchan))
    if (.not. allocated(my_ipol))       allocate(my_ipol(nchan))
    if (.not. allocated(my_side1))      allocate(my_side1(nchan))
    if (.not. allocated(my_side2))      allocate(my_side2(nchan))
    if (.not. allocated(my_hwhmx))      allocate(my_hwhmx(nchan))
    if (.not. allocated(my_chanID))     allocate(my_chanID(nchan))
    if (.not. allocated(nch))           allocate(nch(nfsmp))
    if (.not. allocated(coef))          allocate(coef(nchmax,nfsmp))
    if (.not. allocated(ichmap))        allocate(ichmap(nchmax,nfsmp))
    !---Read optional header content
    HdrOptFlag = adjustl(HdrOptFlag)
    nHdrOpt = len_trim(HdrOptFlag)
    my_chanID = OSS_NULL_ID   ! Default, if absent from file
    allocate(realBuffer(nchan))
    do j=1,nHdrOpt
       select case (HdrOptFlag(j:j))
          case ('C')   ! Center Frequency (GHz)
             READ(ius)  realBuffer
             my_cFreq(1:nchan)=realBuffer

          case ('P')   ! Polarization (V=0, H=1, P=2, M=3, R=4, L=5)
             READ(ius)  my_ipol(1:nchan) ! integer

          case ('D')   ! First sideband freq (GHz)
             READ(ius)  realBuffer
             my_side1(1:nchan)=realBuffer
          case ('E')   ! Second sideband freq (GHz)
             READ(ius)  realBuffer
             my_side2(1:nchan)=realBuffer
          case ('B')   ! Half-Bandwidth (MHz)
             READ(ius)  realBuffer
             my_hwhmx(1:nchan)=realBuffer
          case ('I')   ! Channel ID
             READ(ius)  my_chanID(1:nchan) ! integer
          case ('T')   ! for Tb training
             myRadOrTb=0
          case ('R')   ! for Radiance training
             myRadOrTb=1
          case default
             PRINT *,'err[oss_mw_module::loadOSSselMW]:'
             PRINT *,'Selection file header content flag not recognized: ', &
                HdrOptFlag(j:j)
       STOP
       end select
    enddo
    deallocate(realBuffer)
    !---Read selection data
    DO ismp=1,nfsmp
       READ(ius)iselS,nch(ismp)   ! iselS not needed in present application
       allocate(realBuffer(nch(ismp)))
       READ(ius)realBuffer,ichmap(1:nch(ismp),ismp)
       coef(1:nch(ismp),ismp) = realBuffer
       deallocate(realBuffer)
    END DO
    CLOSE(ius)

    loadOSSselMWdone=.true.
    RETURN
800 PRINT *,'err[oss_mw_module::loadOSSselMW]:  Error opening file: ',TRIM(selFile)
    STOP
  END SUBROUTINE loadOSSselMW

  !----------------------------------------------------------------------------
  !    PURPOSE: Reads the pre-computed Look Up Tables needed later in the
  !    OSS radiative transfer model.
  !----------------------------------------------------------------------------
  SUBROUTINE loadOD(lutFile,nVarMol,varMolID)
    !---Input variables
    CHARACTER(len=*), INTENT(IN)  :: lutFile
    INTEGER,               INTENT(IN)  :: nVarMol
    INTEGER, DIMENSION(:), INTENT(IN)  :: varMolID
    !---Local variables
    INTEGER           :: iuo
    INTEGER           :: NmolS_tmp,imol
    INTEGER           :: uid_od,k,kk
    INTEGER           :: ismp,nsize1,nsize2,nsize3,nsize2_tmp,ks
    INTEGER           :: fmtLUT,istat
    REAL              :: odfac = 0.
    INTEGER, DIMENSION(:) , ALLOCATABLE :: ImolS_tmp
    INTEGER, DIMENSION(:) , ALLOCATABLE :: ImolS_indx
    REAL,    DIMENSION(:) , ALLOCATABLE :: kvar_tmp,dkvar_tmp
    INTEGER, DIMENSION(:) , ALLOCATABLE :: map,mapS,iflag
    ! buffers to read the files
    real*4,  dimension(:), allocatable   :: prefLoc, fixtabLoc
    real*4,  dimension(:,:), allocatable :: tmptabLoc, wvptabLoc, varDfltLoc
    real*4                               :: realVar
    real*4,  dimension(:), allocatable   :: kfixLoc, dkfixLoc, kh2oLoc, dkh2oLoc
    real*4,  dimension(:), allocatable   :: kvarLoc, dkvarLoc
    if(loadODdone) return

    !---Open File
    iuo = findFreeUnit()
    OPEN(iuo,file=lutFile,form='unformatted',status='old',err=800)
    !---Check input molecular selection
    DO k=2,nVarMol
       IF (varMolID(k).LE.varMolID(k-1)) THEN
          PRINT *,'err[oss_mw_module::loadOD]: '// &
            'HITRAN IDs must be in ascending order'
          STOP
       END IF
    END DO
    IF(varMolID(1).NE.1) THEN
          PRINT *,'err[oss_mw_module::loadOD]: '// &
            'Water vapor must be specified in molecular selection'
          STOP
       END IF
    !---Read header of Abs Coeff. file
    READ(iuo,iostat=istat) uid_od,fmtLUT
    if (istat /= 0) then   ! for back compatibility with old format
       fmtLUT = 0  ! format supports no variable gases but H2O
       REWIND iuo
       READ(iuo) uid_od
    endif
    if (fmtLUT < 0 .or. fmtLUT > 1) then
       PRINT *,'err[oss_mw_module::loadOD]:'
       PRINT *,'Format version of OSS LUT file not recognized: ',fmtLUT
       STOP
    endif
    READ(iuo)
    READ(iuo) nmolfix,nmol
    if (.not. allocated(molIDfix)) allocate(molIDfix(nmolfix))
    if (.not. allocated(molID))    allocate(molID(nmol))
    READ(iuo) molidfix(1:nmolfix),molid(1:nmol)
    READ(iuo) nlayod,ntmpod,nwvpod
    !---Sanity checks
    IF (nlayod.GT.mxlay) THEN
       PRINT *,'err[oss_mw_module::loadOD]:  NLay > Mxlay '
       STOP
    END IF
    IF (uid_sel.NE.INACTIVE_UID .AND. uid_sel.NE.uid_od) THEN
       PRINT *,'err[oss_mw_module::loadOD]:  Inconsistent Sel & OD files'
       STOP
    END IF

    nsize1=nlayod*ntmpod*nwvpod
    nsize3=nlayod*ntmpod
    nsize2=nsize3*(nmol-1)
    numRefLev=nlayod+1
    if (.not. allocated(kfix))   allocate(kfix(nsize3,nfsmp))
    if (.not. allocated(dkfix))  allocate(dkfix(nsize3,nfsmp))
    if (.not. allocated(kh2o))   allocate(kh2o(nsize1,nfsmp))
    if (.not. allocated(dkh2o))  allocate(dkh2o(nsize1,nfsmp))
    if (.not. allocated(wvptab)) allocate(wvptab(nwvpod,nlayod))
    if (.not. allocated(fixtab)) allocate(fixtab(nlayod))
    if (.not. allocated(varDflt))allocate(varDflt(nmol,nlayod))
    if (.not. allocated(tmptab)) allocate(tmptab(ntmpod,nlayod))
    if (.not. allocated(pref))   allocate(pref(numRefLev))
    if (.not. allocated(pavlref))   allocate(pavlref(nlayod))
    if (.not. allocated(vfreq))  allocate(vfreq(nfsmp))
    if (.not. allocated(NmolS))  allocate(NmolS(nfsmp))
    if (.not. allocated(ImolS))  allocate(ImolS(nmol,nfsmp))
    if (.not. allocated(ImolS_tmp)) allocate(ImolS_tmp(nmol))
    if (.not. allocated(ImolS_indx))allocate(ImolS_indx(nmol))
    if (.not. allocated(map))       allocate(map(nmol))
    if (.not. allocated(maps))      allocate(maps(nmol))
    if (.not. allocated(iflag))     allocate(iflag(nmol))
    if (.not. allocated(kvar))      allocate(kvar(max(1,nsize2),nfsmp))
    if (.not. allocated(dkvar))     allocate(dkvar(max(1,nsize2),nfsmp))
    if (.not. allocated(kvar_tmp))  allocate(kvar_tmp(max(1,nsize2)))
    if (.not. allocated(dkvar_tmp)) allocate(dkvar_tmp(max(1,nsize2)))

    IF (fmtLUT >= 1) THEN
       allocate(prefLoc(numRefLev),tmptabLoc(ntmpod,nlayod),&
            fixtabLoc(nlayod), wvptabLoc(nwvpod,nlayod),&
            varDfltLoc(nmol,nlayod))

       READ(iuo)prefLoc,tmptabLoc,fixtabLoc,wvptabLoc,varDfltLoc

       pref=prefLoc
       tmptab=tmptabLoc
       fixtab=fixtabLoc
       wvptab=wvptabLoc
       varDflt=varDfltLoc

       deallocate(prefLoc,tmptabLoc,&
            fixtabLoc, wvptabLoc,&
            varDfltLoc)

    ELSE  ! for back compatibility
       allocate(prefLoc(numRefLev),tmptabLoc(ntmpod,nlayod),&
            fixtabLoc(nlayod), wvptabLoc(nwvpod,nlayod))

       READ(iuo)prefLoc,tmptabLoc,fixtabLoc,wvptabLoc

       pref=prefLoc
       tmptab=tmptabLoc
       fixtab=fixtabLoc
       wvptab=wvptabLoc

       deallocate(prefLoc,tmptabLoc,&
            fixtabLoc, wvptabLoc)
    ENDIF
    pavlref(1:numRefLev-1)=0.5*(pref(1:numRefLev-1)+pref(2:numRefLev))

    !---Build mapping vector for user-selected variable molecules
    kk=0
    map(1:nmol)=0
    DO k=1,nmol
       !---map(k) contains new relative index for selected molecule
       ! and is zero for others
       IF(ANY (varMolID(1:nVarMol) == molid(k))) THEN
          kk=kk+1
          map(k)=kk
       END IF
    END DO
    IF(kk.NE.nVarMol)THEN
       PRINT *,'This database contains only the molecules with'
       PRINT *,'the following HITRAN ID as variable molecules: '
       PRINT *,molid(1:nmol)
       STOP 'loadOD: '
    END IF
    IF (fmtLUT >= 1) THEN
      allocate(kfixLoc(nsize3), dkfixLoc(nsize3),&
               kh2oLoc(nsize1), dkh2oLoc(nsize1),&
                 kvarLoc(nsize2), dkvarLoc(nsize2))
    else
      allocate(kfixLoc(nsize1), dkfixLoc(nsize1),&
               kh2oLoc(nsize1), dkh2oLoc(nsize1))
      kvar_tmp = 0.0
      dkvar_tmp = 0.0
      kfixLoc = 0.0
      dkfixLoc = 0.0
    end if
    DO ISmp=1,NFSmp
       IF (fmtLUT >= 1) THEN
          READ(iuo)realVar,NmolS_tmp
          vfreq(ismp) = realVar

          READ(iuo)ImolS_tmp(1:NmolS_tmp)
          nsize2_tmp=nsize3*(NmolS_tmp-1)

          READ(iuo)kfixLoc(1:nsize3),dkfixLoc(1:nsize3),&
                   kh2oLoc(1:nsize1),dkh2oLoc(1:nsize1),&
                   kvarLoc(1:nsize2_tmp),dkvarLoc(1:nsize2_tmp)
          kfix(1:nsize3,ismp)=kfixLoc(1:nsize3)
          dkfix(1:nsize3,ismp)=dkfixLoc(1:nsize3)
          kh2o(1:nsize1,ismp)=kh2oLoc(1:nsize1)
          dkh2o(1:nsize1,ismp)=dkh2oLoc(1:nsize1)
          kvar_tmp(1:nsize2_tmp)=kvarLoc(1:nsize2_tmp)
          dkvar_tmp(1:nsize2_tmp)=dkvarLoc(1:nsize2_tmp)
       ELSE  ! for back compatibility
          READ(iuo)realVar
          vfreq(ismp) = realVar

          READ(iuo)kfixLoc(1:nsize1),dkfixLoc(1:nsize1),&
                   kh2oLoc(1:nsize1),dkh2oLoc(1:nsize1)

          kh2o(1:nsize1,ismp)=kh2oLoc(1:nsize1)
          dkh2o(1:nsize1,ismp)=dkh2oLoc(1:nsize1)

          NmolS_tmp=1
          ImolS_tmp(1)=1 !  first per smp is first overall
          kfix(1:nsize3,ismp) =kfixLoc(1:nsize3)   ! Use vapor grid 1
          dkfix(1:nsize3,ismp)=dkfixLoc(1:nsize3)
       ENDIF
       CALL init_kfix(kfix(1,ismp),fixtab)
       CALL init_kfix(dkfix(1,ismp),fixtab)

       !     Keep only user-selected molecules as variables molecules -
       !     merge remaining molecules with fixed gases
       !     First apply OD threshold to identify weakly absorbing constituents
       iflag(1)=1
       iflag(2:nmol) = 0
       CALL odthresh(kfix(1,ismp),dkfix(1,ismp),kvar_tmp,dkvar_tmp, &
          ImolS_tmp,NmolS_tmp,wvptab,varDflt,odfac,iflag)
       maps(1:nmol) = map(1:nmol)*iflag(1:nmol)
       kk = 0
       DO ks = 1,NmolS_tmp
          imol = ImolS_tmp(ks)
          IF(mapS(imol) > 0)THEN
             kk = kk+1
             ImolS_indx(kk) = ks
          ELSE
             CALL cum_fix(kfix(1,ismp),kvar_tmp,imol,NmolS_tmp,ks,varDflt)
             CALL cum_fix(dkfix(1,ismp),dkvar_tmp,imol,NmolS_tmp,ks,varDflt)
          END IF
       END DO
       NmolS(ismp) = kk
       CALL shrink_var(kvar_tmp,ImolS_tmp,NmolS_tmp,mapS,ImolS_indx,&
            kvar(1,ISmp),ImolS(1,ISmp),NmolS(ISmp))
       CALL shrink_var(dkvar_tmp,ImolS_tmp,NmolS_tmp,mapS,ImolS_indx,&
            dkvar(1,ISmp),ImolS(1,ISmp),NmolS(ISmp))
    END DO
    deallocate(kfixLoc, dkfixLoc, kh2oLoc, dkh2oLoc)

    IF (fmtLUT >= 1) deallocate(kvarLoc, dkvarLoc)
    CLOSE(iuo)

    nmol = nVarMol
    molID(1:nmol) = varMolID(1:nmol)

    deallocate (kvar_tmp,dkvar_tmp,ImolS_tmp,ImolS_indx,map,maps,iflag)
    loadODdone=.true.
    RETURN
800 PRINT *,'err[oss_mw_module::loadOD]:  Error opening file: ',TRIM(lutFile)
    STOP
  END SUBROUTINE loadOD


  SUBROUTINE init_kfix(kfix0,rfix)
    !---In/Out variables
    REAL,             INTENT(INOUT) :: kfix0(NlayOD,NtmpOD)
    REAL, DIMENSION(:),  INTENT(IN) :: rfix
    !---Local variables
    INTEGER    :: nt

    DO nt=1,ntmpOD
       kfix0(1:NLayOD,nt) = kfix0(1:NLayOD,nt)*rfix(1:NlayOD)
    END DO
    RETURN
  END SUBROUTINE init_kfix


  SUBROUTINE odthresh(kfix0,dkfix0,kvar0,dkvar0,ImolS,NmolS,wvap,w,odfac,iflag)
    !---Input variables
    REAL, DIMENSION(NmolS-1,NlayOD,NtmpOD), &
                                    INTENT(IN) :: kvar0,dkvar0
    REAL, DIMENSION(NlayOD,NtmpOD), INTENT(IN) :: kfix0,dkfix0
    INTEGER,                        INTENT(IN) :: NmolS,ImolS(NmolS)
    REAL, DIMENSION(:,:),           INTENT(IN) :: wvap
    REAL, DIMENSION(:,:),           INTENT(IN) :: w
    REAL,                           INTENT(IN) :: odfac
    !---Output variables
    INTEGER, DIMENSION(:),       INTENT(INOUT) :: iflag
    !---Local variables
    REAL                                   :: wvref
    REAL                                   :: kbuf(50),ktot
    INTEGER                                :: nt,l,ks
    !     Compare optical depth of individual dry constituents to
    !     total optical depth and flag weakly absorbing molecules
    !     Uses abs coef*MR as proxy for OD.

    DO nt = 1,NtmpOD
       DO l = 1,NlayOD
          wvref=(wvap(1,l)+wvap(Nwvpod,l))/2.  ! use middle of h2o grid
          ktot = wvref*dkfix0(l,nt)+kfix0(l,nt)  ! already includes fixed MR
          DO ks = 2,NmolS
             kbuf(ks-1) = (wvref*dkvar0(ks-1,l,nt)+kvar0(ks-1,l,nt))*w(ImolS(ks),l)
             ktot = ktot + kbuf(ks-1)
          END DO
          DO ks=2,NmolS  ! could exit at first occurence of iflag = 1
             IF(kbuf(ks-1) > odfac*ktot) iflag(ImolS(ks)) = 1
          END DO
       END DO
    END DO
    RETURN
  END SUBROUTINE odthresh


  SUBROUTINE cum_fix(kfix0,kvar0,imol,NmolS,ks,w)
    !---In/Out variables
    INTEGER,               INTENT(IN)    :: ks
    INTEGER,               INTENT(IN)    :: NmolS,Imol
    REAL,                  INTENT(IN)    :: kvar0(NmolS-1,NlayOD,NtmpOD)
    REAL, DIMENSION(:,:),  INTENT(IN)    :: w
    REAL,                  INTENT(INOUT) :: kfix0(NlayOD,NtmpOD)
    !---Local variables
    INTEGER    :: nt

    DO nt=1,ntmpOD
       kfix0(1:NLayOD,nt) = kfix0(1:NLayOD,nt) +&
            kvar0(ks-1,1:NLayOD,nt)*w(imol,1:nlayod)
    END DO
    RETURN
  END SUBROUTINE cum_fix


  SUBROUTINE shrink_var(kvar0,ImolS0,NmolS0,maps,ImolS_indx,kvar1,ImolS1,NmolS1)
    !---In/Out variables
    INTEGER,               INTENT(IN)    :: NmolS0
    INTEGER,               INTENT(IN)    :: ImolS0(NmolS0)
    INTEGER,               INTENT(IN)    :: NmolS1
    INTEGER,               INTENT(INOUT) :: ImolS1(NmolS1)
    INTEGER, DIMENSION(:), INTENT(IN)    :: maps,ImolS_indx
    REAL,                  INTENT(IN)    :: kvar0(NmolS0-1,NlayOD,NtmpOD)
    REAL,                  INTENT(INOUT) :: kvar1(NmolS1-1,NlayOD,NtmpOD)

    kvar1(1:NmolS1-1,1:NLayOD,1:NtmpOD) = &
         kvar0(ImolS_indx(2:NmolS1)-1,1:NLayOD,1:NtmpOD)
    ImolS1(1:NmolS1) = maps(ImolS0(ImolS_indx(1:NmolS1)))
    RETURN
  END SUBROUTINE shrink_var

  !----------------------------------------------------------------------------
  ! PURPOSE: Computes the cosmic background, in a way consistent with
  !          linear function of radiance.
  !----------------------------------------------------------------------------
  SUBROUTINE CosmicBackg(f,cosmbk,alpha,beta)
    REAL, PARAMETER :: thetaC=2.73
    !---Input/Output variables
    REAL,    INTENT(IN) :: f
    REAL, INTENT(INOUT) :: cosmbk,alpha,beta
    !---Local variables
    REAL            :: radC
    radC        = Planck(f,thetaC)
    call compPlanckLin(f,alpha,beta)
    cosmbk = alpha*radC+beta
    RETURN
  END SUBROUTINE CosmicBackg

  FUNCTION Planck(f,t)
    !f: GHz
    !t: Kelvin
    !Planck: mW/m^2/str/GHz
    REAL(DOUBLE), PARAMETER :: c1=GHzToHz*2._DOUBLE*Plnk/(LtSp*LtSp), c2=Plnk/Bltz
    REAL(DOUBLE)            :: xnumer,xdenom,fHz
    REAL, INTENT(IN)        :: f,t
    REAL                    :: Planck
    fHz        = REAL(f,DOUBLE)*GHzToHz
    xnumer     = c1*(fHz**3)
    xdenom     = EXP((c2*fHz)/REAL(t,DOUBLE))-1._DOUBLE
    Planck = REAL(xnumer/xdenom)
    RETURN
  END FUNCTION Planck

  !----------------------------------------------------------------------------
  ! PURPOSE: Return the coefficients for linear function of radiance
  !          at the center frequency of each channel
  !----------------------------------------------------------------------------
  SUBROUTINE getPlanckLinCtr(alphaCtr,betaCtr)
    REAL, DIMENSION(:), INTENT(INOUT) :: alphaCtr,betaCtr
    INTEGER :: ich

    IF (.NOT. loadOSSselMWdone)  THEN
      PRINT *,'err:[oss_mw_module::getPlanckLinCtr] '// &
         'OSS data must be loaded first'
      STOP
    ENDIF
    IF (SIZE(alphaCtr) < nchan .or. SIZE(betaCtr) < nchan) THEN
      PRINT *,'err:[oss_mw_module::getPlanckLinCtr] Output arrays too small;'
      PRINT *,'SIZE(alphaCtr),SIZE(betaCtr),nchan :', &
         SIZE(alphaCtr),SIZE(betaCtr),nchan
      STOP
    ENDIF

    DO ich=1,nchan
      CALL compPlanckLin(my_cFreq(ich),alphaCtr(ich),betaCtr(ich))
    ENDDO
    RETURN

  END SUBROUTINE getPlanckLinCtr

  ! To use 2nd-order Taylor expansion of Planck, this must
  ! be called before oss_init_mw
  SUBROUTINE set2ndOrderPlanck()
    use2ndOrderPlanck=.TRUE.
    RETURN
  END SUBROUTINE set2ndOrderPlanck

  !----------------------------------------------------------------------------
  ! PURPOSE: Computes the coefficients for linear function of radiance.
  !----------------------------------------------------------------------------
  SUBROUTINE compPlanckLin(f,alpha,beta)
    ! These theta minimize difference from Planck over range 150-350 K:
    REAL,    INTENT(IN)    :: f
    REAL,    INTENT(INOUT) :: alpha,beta
    REAL, PARAMETER :: theta1=170.,theta2=310.
    REAL            :: rad1,rad2,fHz

    IF (.NOT. use2ndOrderPlanck) then  ! 2-point linear fit to Planck
      rad1        = Planck(f,theta1)
      rad2        = Planck(f,theta2)
      alpha       = (theta2-theta1)/(rad2-rad1)
      beta        = (theta1*rad2-theta2*rad1)/(rad2-rad1)
    ELSE                               ! 2nd-order Taylor expansion of Planck
      fHz         = f*REAL(GHzToHz)
      alpha       = REAL(LtSp**2)/(REAL(GHzToHz)*2.*(fHz)**2*REAL(Bltz))
      beta        = REAL(Plnk)*f*GHzToHz/(2.*REAL(Bltz))
    ENDIF
    RETURN
  END SUBROUTINE compPlanckLin

  !----------------------------------------------------------------------------
  ! PURPOSE: Computes moist air density
  !----------------------------------------------------------------------------
  subroutine lzsum4tot(p,t,rhotot)
    real, parameter     :: R=8.3144598           !  gas constant
    real, intent(in)    :: p,t
    real, intent(inout) :: RHOTOT                ! total density in g/cm^3 ???
    RHOTOT=p/t/R
    return
  end subroutine lzsum4tot
  !----------------------------------------------------------------------------
  ! PURPOSE: Computes layer average molecular amount
  !          using a zsum.
  !----------------------------------------------------------------------------
  subroutine lzsum4W(xu,xl,dz,RHOTOT1,RHOTOT2,xint,dxu,dxl)
    real,   intent(in) :: xu,xl,dz
    real,   intent(in) :: RHOTOT1     !total density (molec/cm^3) at upper level
    real,   intent(in) :: RHOTOT2     !total density (molec/cm^3) at lower level
    real, intent(inout):: xint        !column amount (molec/cm^2)
    real, intent(inout):: dxu,dxl

    real               :: RHOAIR_u ! density on the upper level (molec/cm^3)
    real               :: RHOAIR_l ! density on the lower level (molec/cm^3)
    real               :: HZ       ! inverse exponential
    RHOAIR_u=RHOTOT1*xu
    RHOAIR_l=RHOTOT2*xl

    HZ      =DZ/log(RHOAIR_u/RHOAIR_l)
    xint    =HZ*(RHOAIR_u-RHOAIR_l)*1.E+1
    dxu     =xint*(RHOTOT1/(RHOAIR_u-RHOAIR_l)-1.0/log(RHOAIR_u/RHOAIR_l)/xu)
    dxl     =xint*(-RHOTOT2/(RHOAIR_u-RHOAIR_l)+1.0/log(RHOAIR_u/RHOAIR_l)/xl)

    return
  end subroutine lzsum4W

  !----------------------------------------------------------------------------
  ! PURPOSE: Computes average layer temperature
  !          using a zsum.
  !----------------------------------------------------------------------------
  subroutine lzsum4T(xu,xl,xint,dxu,dxl)
    real,   intent(in) :: xu,xl
    real, intent(inout):: xint,dxu,dxl
    real               :: eps
    eps=xl/xu-1
    if (abs(eps) > 1.e-2) then
       xint    =(xl-xu)/log(xl/xu)
       dxu     =xint*(1.0/(xu-xl)+1.0/xu/log(xl/xu))
       dxl     =xint*(1.0/(xl-xu)+1.0/xl/log(xu/xl))
    else
       xint    =0.5*(xl+xu-xu*eps*eps/6.0)
       dxu     =0.5+eps/6.0
       dxl     =0.5-eps/6.0
    end if
    return
  end subroutine lzsum4T
  !----------------------------------------------------------------------------
  ! PURPOSE: Computes average layer quantities (or integrated amount)
  !          using a linear dependence on P.
  !----------------------------------------------------------------------------
  SUBROUTINE lpsum(pu,pl,xu,xl,scal,xint,dxu,dxl)
    !---Input variables
    REAL, INTENT(IN)  :: pu,pl,xu,xl,scal
    !---Output variables
    REAL, INTENT(OUT) :: xint,dxu,dxl

    xint = 0.5*(pl-pu)*scal
    dxu      = xint
    dxl      = xint
    xint     = xint*(xu+xl)
    RETURN
  END SUBROUTINE lpsum

  !----------------------------------------------------------------------------
  ! PURPOSE: Computes average layer quantities (or integrated amount)
  !          alza  = math.log(zeta)
  !          using a log-x dependence on log-p.
  !----------------------------------------------------------------------------
  SUBROUTINE lpsum_log(pu,pl,xu,xl,scal,xint,dxu,dxl)
    REAL, PARAMETER :: epsiln=1.e-12
    !---Input variables
    REAL, INTENT(IN)  :: pu,pl,xu,xl,scal
    !---Output variables
    REAL, INTENT(OUT) :: xint,dxu,dxl
    !---Local variables
    REAL              :: hp,x0,zeta,alza,alpha
    hp       = alog(pl/pu)
    x0       = pl*xl*hp
    zeta     = pu*xu/(pl*xl)
    IF(ABS(zeta-1.).GT.epsiln)THEN
       alza  = alog(zeta)
       xint  = x0*(zeta-1.)/alza
       alpha = zeta/(zeta-1.)-1./alza
    ELSE
       xint  = x0*2./(3.-zeta)
       alpha = zeta/(3.-zeta)
    END IF
    xint     = xint*scal
    dxu      = xint*alpha/xu
    dxl      = xint*(1.-alpha)/xl
    RETURN
  END SUBROUTINE lpsum_log

  SUBROUTINE lint(x,p,N0,p0,x0)
    INTEGER,            INTENT(IN) :: n0
    REAL,               INTENT(IN) :: p0
    REAL, DIMENSION(:), INTENT(IN) :: x,p
    !---Output variable
    REAL,            INTENT(INOUT) :: x0
    !---Local variable
    REAL                           :: xx
    xx = (x(N0)-x(N0-1))/(p(N0)-p(N0-1))
    x0 = x(N0-1)+(p0-p(N0-1))*xx
    RETURN
  END SUBROUTINE lint

  SUBROUTINE lint_log(x,p,N0,p0,x0)
    INTEGER,            INTENT(IN) :: n0
    REAL,               INTENT(IN) :: p0
    REAL, DIMENSION(:), INTENT(IN) :: x,p
    !---Output variable
    REAL,            INTENT(INOUT) :: x0
    !---Local variable
    REAL                           :: xx
    xx = alog(x(N0)/x(N0-1))/alog(p(N0)/p(N0-1))
    x0 = x(N0-1)*(p0/p(N0-1))**xx
    RETURN
  END SUBROUTINE lint_log

  !----------------------------------------------------------------------------
  ! PURPOSE: Computes items related to the viewing geometry.
  !----------------------------------------------------------------------------
  SUBROUTINE setpath(pSurf,pLoc,obsAngle,nSurf,nObs,lookup,n1,n2,umu,referenceGrid,&
                    interpSfc,pUser)
    !---Input variables
    REAL, INTENT(IN)                       :: pSurf
    REAL,    DIMENSION(:), INTENT(IN)      :: obsAngle
    LOGICAL,               INTENT(IN)      :: lookup
    LOGICAL,               INTENT(OUT)     :: referenceGrid
    LOGICAL,               INTENT(OUT)     :: interpSfc
    real, DIMENSION(:),OPTIONAL,INTENT(IN) :: pUser
    INTEGER,               INTENT(IN)      :: nObs

    !---Output variables
    INTEGER,               INTENT(INOUT) :: nSurf,n1,n2
    REAL,    DIMENSION(:), INTENT(INOUT) :: umu
    real, DIMENSION(:),    INTENT(INOUT) :: pLoc

    !---Local variables
    INTEGER                              :: npLev
    REAL, DIMENSION(SIZE(obsAngle))      :: viewang

    if (present(pUser)) then
      referenceGrid = .false.
      npLev = size(pUser)
      pLoc(1:npLev)  = pUser(1:npLev)
      if (abs(pSurf-pUser(npLev)) < (pUser(npLev)-pUser(npLev-1))* &
              pMatchLim) then
        interpSfc =.false.
      else
        interpSfc = .true.
      end if
    else
      interpSfc = .true.
      referenceGrid = .true.
      npLev=numRefLev
      pLoc(1:npLev)=pRef(1:npLev)
    end if

    DO nSurf=2,npLev-1
     IF(pLoc(nSurf)>=psurf) EXIT
    END DO
    if (nSurf>npLev) nSurf=npLev

    !---Set viewing geometry parameters
    viewang=obsAngle
    IF(lookup)THEN
       PRINT *,'err[oss_mw_module::setpath]:  Lookup=TRUE not supported yet.'
       STOP
       n1=nSurf
       n2=nObs-1
    ELSE
       IF(ANY(viewang > 90.)) THEN
          PRINT *,'err[oss_mw_module::setpath]:  EIA must be less than 90'
          STOP
       END IF
       n1=nObs
       n2=nSurf-1
    END IF
    umu =COS(viewang*deg2rad)
    RETURN
  END SUBROUTINE setpath

  !----------------------------------------------------------------------------
  ! PURPOSE: Based upon the features of the input cloud, i.e.,
  !          the cloud top, cloud thickness and amount, the
  !          function calculates MW cloud optical depths at
  !          each monochromatic channel for liquid cloud
  !----------------------------------------------------------------------------
  SUBROUTINE osstc_mw(NSurf,Psfc,press,xG,tavl,taucld,abscld, &
       dtaucdtop,dtaucdthick)
    REAL, PARAMETER :: gg=9.8,rdry=287.
    !---Input variables
    INTEGER,                 INTENT(IN)    :: NSurf
    REAL,                    INTENT(IN)    :: psfc
    REAL,    DIMENSION(:),   INTENT(IN)    :: press
    REAL,    DIMENSION(:),   INTENT(IN)    :: xG,tavl
    !---Output variables
    REAL,    DIMENSION(:,:), INTENT(INOUT) :: taucld,abscld,dtaucdtop,dtaucdthick
    !---Local variables
    INTEGER, PARAMETER :: nzinc=5
    REAL,    PARAMETER :: zinc=500.,tstd=280.,plinc=zinc*gg/(rdry*tstd)
    INTEGER, DIMENSION(nzinc)       :: j1,j2
    REAL,    DIMENSION(Mxlay)       :: absrt,sdtop,sdthk
    REAL,    DIMENSION(Mxlay,nzinc) :: delp,abs1
    INTEGER                         :: k,nlast,kl1,kl2,klo,kl,klok,j
    REAL                            :: frqtmp,ptop,pthick,clw
    REAL                            :: cfrac,pbot,ptopk,pbotk,abs0,swt
    !----Initializations
    taucld      = 0.
    abscld      = 0.
    dtaucdtop   = 0.
    dtaucdthick = 0.
    Nlast       = NSurf-1
    ptop        = xG(ICldLiqG)
    pthick      = xG(ICldLiqG+1)
    clw         = xG(ICldLiqG+2)
    cfrac       = 1.0
    pbot        = MIN(ptop+pthick,Psfc)

    !---pbot may cause index out of bound.
    pbot        = MIN(pbot, press(NSurf))
    pthick      = pbot-ptop

    !--Make comput. for nzinc cld positions, including actual one,
    !--to stabilize derivatives with respect to cld top and thick
    kl1=1
    kl2=nzinc
    klo=1+nzinc/2
    DO kl=kl1,kl2
       ptopk=EXP(alog(ptop)+float(kl-klo)*plinc)
       pbotk=ptopk+pthick
       IF(kl.LE.klo.OR.pbotk.LE.Psfc)THEN
          klOK=kl
          !---Find index of first level below cloud top
          DO j=1,NSurf
             IF(ptopk.LT.press(j)) go to 10
          END DO
10        j1(kl)=MIN(j,NSurf)
          !--Find index of first level at or below cloud base.
          DO j=j1(kl)-1,NSurf
             IF(pbotk.LE.press(j)) go to 20
          END DO
20        j2(kl)=MIN(j,NSurf)
          !---Calculate cloud thickness in each layer
          delp(j1(kl)-1,kl)=press(j1(kl))-ptopk ! Top layer
          DO j=j1(kl),j2(kl)-2
             delp(j,kl)=press(j+1)-press(j)
          ENDDO
          delp(j2(kl)-1,kl)=pbotk-amax1(press(j2(kl)-1),ptopk) !Bottom layer
       ENDIF
    ENDDO
    kl2=min0(kl2,klOK)
    !---Loop through channels
    DO k=1,nchan
       frqTmp=my_cfreq(k)
       !---Cloud layers
       DO j=j1(kl1)-1,j2(kl2)-1
          abs0=abcoefliq(frqTmp,tavl(j))
          absrt(j)=abs0*cfrac/pthick
       ENDDO
       sdtop(1:nlast)=0.
       sdthk(1:nlast)=0.
       swt=0.
       DO kl=kl1,kl2
          DO j=j1(kl)-1,j2(kl)-1
             abs1(j,kl)=absrt(j)*delp(j,kl)
             sdthk(j)=sdthk(j)-clw*abs1(j,kl)/pthick
          ENDDO
          sdtop(j1(kl)-1)=sdtop(j1(kl)-1)-clw*absrt(j1(kl)-1)
          sdtop(j2(kl)-1)=sdtop(j2(kl)-1)+clw*absrt(j2(kl)-1)
          sdthk(j2(kl)-1)=sdthk(j2(kl)-1)+clw*absrt(j2(kl)-1)
          swt=swt+1.
       ENDDO
       !---Derivatives are weighted avges for the ninc posits.
       DO j=1,NLast
          dtaucdtop(j,k)=sdtop(j)/swt
          dtaucdthick(j,k)=sdthk(j)/swt
       ENDDO
       !---Cloud layers for actual cloud position
       DO j=j1(klo)-1,j2(klo)-1
          abscld(j,k)=abs1(j,klo)
          taucld(k,j)=clw*abs1(j,klo)
       END DO
    ENDDO
    RETURN
  END SUBROUTINE osstc_mw
  !----------------------------------------------------------------------------
  ! PURPOSE:  This subroutine calculates the average temperature and molecular
  !           amounts for all layers for given profiles of temperature and
  !           molecular concentrations. It also calculates the derivatives of
  !           tavl with respect to a change in the lower and upper boundary
  !           temperatures and the derivatives of the molecular amounts with
  !           respect to a change in the mixing ratios at the layer
  !           boundaries.
  !           Integration assumes that T is linear in z (LnT linear in LnP)
  !           and LnQ linear in LnP.
  !----------------------------------------------------------------------------
  SUBROUTINE fpathP_mw(xG,lat,Nlast,tavl,pavl,wdry,qvar,wvar,dtu,dtl,dwqu,dwql,&
    referenceGrid,interpSfc,alt,pLoc)
    !---Input variables
    REAL, DIMENSION(:), INTENT(IN)    :: xG
    real              , intent(in)    :: lat
    INTEGER,            INTENT(IN)    :: nlast
    LOGICAL,            INTENT(IN)    :: referenceGrid
    LOGICAL,            INTENT(IN)    :: interpSfc
    REAL, DIMENSION(:), INTENT(IN)    :: alt
    REAL, DIMENSION(:), INTENT(IN)    :: pLoc
    !---Output variables
    REAL, DIMENSION(:), INTENT(INOUT)   :: tavl,pavl,wdry,dtu,dtl
    REAL, DIMENSION(:,:), INTENT(INOUT) :: qvar,wvar
    REAL, DIMENSION(:,:), INTENT(INOUT) :: dwqu,dwql
    !---Local variables
    INTEGER                           :: nSurf,l,k,IXoff,nEnd
    REAL                              :: psfc,tsfc,scal,scal2,qsfc,wtot
    REAL, DIMENSION(mxlev)            :: dp,qcor,qr
    REAL, DIMENSION(mxlay)            :: rairLay
    real                              :: rair

    nSurf=Nlast+1
    !------------------------------------------------------
    !     Compute Average Pressure for the layers
    !------------------------------------------------------
    Psfc=xG(IPsfcG)
    if (referenceGrid) then
       nEnd = Nlast-1
    else
       pavl(1:nlast) = 0.5*(pLoc(1:nlast)+pLoc(2:nSurf))
       nEnd = Nlast
    endif
    !------------------------------------------------------
    !     Compute Average Temperature for the layers
    !------------------------------------------------------
    DO l=1,nEnd
       dp(l)=(pLoc(l+1)-pLoc(l))
       scal=1./dp(l)
       CALL lpsum_log(pLoc(l),pLoc(l+1),xG(ITempG+l-1),xG(ITempG+l),&
            scal,tavl(l),dtu(l),dtl(l))
    END DO
    dp(Nlast)=(PSfc-pLoc(nSurf-1))
    if (interpSfc) then
    !---Surface layer
       scal=1./dp(Nlast)
       scal2=alog(PSfc/pLoc(nSurf))/alog(pLoc(nSurf-1)/pLoc(nSurf))
       CALL lint_log(xG(ITempG:ITempG+nSurf-1),pLoc,nSurf,Psfc,Tsfc)
       CALL lpsum_log(pLoc(nSurf-1),PSfc,xG(nSurf-1),TSfc,&
           scal,tavl(Nlast),dtu(Nlast),dtl(Nlast))
      dtu(Nlast)=dtu(Nlast)+dtl(Nlast)*scal2*tsfc/xG(nSurf-1)
      dtl(Nlast)=dtl(Nlast)*(1.0-scal2)*tsfc/xG(nSurf)
    end if
    !------------------------------------------------------------------------
    !     Calculate amounts for individual species and derivatives wrt level
    !     mixing ratios for retrieved constituents.
    !------------------------------------------------------------------------
    rair=10./(grav_const_req1 + grav_const_req2*COS(2.0*deg2rad*LAT))
    rairLay(1:Nlast) = rair*(1. + alt(1:Nlast)/rEarth)**2
    ! to run the tests uncomment two lines to be compitible with the old version
    ! rair=1.020408163
    ! rairLay(1:Nlast) = rair

    IXoff=ImolG(1)-1  ! water vapor is always molecule 1
    qcor(1:nSurf)=1.0/(1.+xG(IXoff+1:IXoff+nSurf))
    DO k=1,nmol
       IXoff=ImolG(k)-1
       !---Transform mix. ratios into mass fractions (relative to total air mass)
       qr(1:nSurf)=xG(IXoff+1:IXoff+nSurf)*qcor(1:nSurf)
       DO l=1,nEnd
        scal=rairLay(l)
        CALL lpsum_log(pLoc(l),pLoc(l+1),qr(l),qr(l+1), &
               scal,wvar(k,l),dwqu(l,k),dwql(l,k))
       END DO
       if (interpSfc) then
         !---Surface Layer
         CALL lint_log(qr,pLoc,nSurf,psfc,qsfc)
           scal=rairLay(Nlast)
         CALL lpsum_log(pLoc(nSurf-1),psfc,qr(nSurf-1),qsfc,&
            scal,wvar(k,Nlast),dwqu(Nlast,k),dwql(Nlast,k))
         dwqu(Nlast,k)=dwqu(Nlast,k)+ dwql(Nlast,k)*scal2*qsfc/qr(nSurf-1)
         dwql(Nlast,k)=dwql(Nlast,k)*(1.0-scal2)*qsfc/qr(nSurf)
       end if
    END DO
    !---Derivatives of amount wrt dry mixing ratios
    DO l=1,Nlast
       dwqu(l,1)=dwqu(l,1)*qcor(l)**2
       dwql(l,1)=dwql(l,1)*qcor(l+1)**2
       dwqu(l,2:nmol)=dwqu(l,2:nmol)*qcor(l)
       dwql(l,2:nmol)=dwql(l,2:nmol)*qcor(l+1)
    END DO
    !---Compute dry gas amount and average layer mixing ratios
    DO l=1,Nlast
       wtot=dp(l)*rair
       wdry(l)=wtot-wvar(1,l)
       qvar(1,l)=wvar(1,l)/wtot
       qvar(2:nmol,l)=wvar(2:nmol,l)/wdry(l)
    END DO
    RETURN
  END SUBROUTINE fpathP_mw

  SUBROUTINE fpathZ_mw(xG,Nlast,tavl,pavl,wdry,qvar,wvar,dtu,dtl,dwqu,dwql,&
    referenceGrid,interpSfc,pLoc,alt)
    !---Input variables
    REAL, DIMENSION(:), INTENT(IN)    :: xG
    INTEGER,            INTENT(IN)    :: nlast
    LOGICAL,            INTENT(IN)    :: referenceGrid
    LOGICAL,            INTENT(IN)    :: interpSfc
    REAL, DIMENSION(:), INTENT(IN)    :: pLoc

    real, dimension(mxlay)            :: alt

    !---Output variables
    REAL, DIMENSION(:), INTENT(INOUT)   :: tavl,pavl,wdry,dtu,dtl
    REAL, DIMENSION(:,:), INTENT(INOUT) :: qvar,wvar
    REAL, DIMENSION(:,:), INTENT(INOUT) :: dwqu,dwql
    !---Local variables
    INTEGER                           :: nSurf,l,k,IXoff,nEnd, iH2o, n
    REAL                              :: psfc,tsfc,scal,qsfc
    REAL, DIMENSION(mxlev)            :: qcor,qr
    REAL, DIMENSION(mxlay)            :: wtot
    real                              :: plogu, plogl, alpha, dzsfc
    real                              :: RHOTOT1, RHOTOT2, dummyReal

    nSurf=Nlast+1
    !------------------------------------------------------
    !     Compute Average Pressure for the layers
    !------------------------------------------------------
    Psfc=xG(IPsfcG)
    if (referenceGrid) then
       nEnd = Nlast-1
    else
       pavl(1:nlast) = 0.5*(pLoc(1:nlast)+pLoc(2:nSurf))
       nEnd = Nlast
    endif
    !------------------------------------------------------
    !     Compute Average Temperature for the layers
    !------------------------------------------------------
    iXoff = ITempG-1

    DO l=1,nEnd
       call lzsum4T(xG(iXoff+l),xG(iXoff+l+1),tavl(l),dtu(l),dtl(l))
    END DO
    if (interpSfc) then
    !---Surface layer
      plogu=log(Psfc/pLoc(nSurf-1))
      plogl=log(Psfc/pLoc(nSurf))
      alpha = plogl/(plogl-plogu)

      tSfc=(xG(iXoff+nSurf-1)**alpha)*(xG(iXoff+nSurf)**(1.0-alpha))
      call lzsum4T(xG(iXoff+nSurf-1),tSfc,tavl(Nlast),dtu(Nlast),dtl(Nlast))

      dtu(Nlast)=dtu(Nlast)+dtl(Nlast)*alpha*tsfc/xG(nSurf-1)
      dtl(Nlast)=dtl(Nlast)*(1.0-alpha)*tsfc/xG(nSurf)
    end if
    !------------------------------------------------------------------------
    !     Calculate amounts for individual species and derivatives wrt level
    !     mixing ratios for retrieved constituents.
    !------------------------------------------------------------------------
    iH2o=ImolG(1)-1

    ! correction factor from dry to moist
    ! and calculate the total air
    do l=1,nSurf
      qcor(l)=1.0/(1.0 + xG(iH2o+l)/normMolWt(molID(1)) )
      qr(l)=qcor(l)*(1.0 + xG(iH2o+l))
    end do
    call lzsum4tot(pLoc(1),xG(1),RHOTOT1)

    do l=1,nEnd
      call lzsum4tot(pLoc(l+1),xG(l+1),RHOTOT2)
      call lzsum4W(qr(l), qr(l + 1), alt(l) - alt(l + 1), &
               RHOTOT1, RHOTOT2, wtot(l), dummyReal, dummyReal)
      wtot(l) = wtot(l)* drymwt
      RHOTOT1 = RHOTOT2
    end do
    if (interpSfc) then
      qsfc=(qr(nSurf-1)**alpha)*(qr(nSurf)**(1.0-alpha))
      call lzsum4tot(pSfc,tSfc,RHOTOT2)

      dzsfc= alt(nSurf-1) - alt(nSurf)
      call lzsum4W(qr(nSurf-1), qsfc, dzsfc, &
                 RHOTOT1, RHOTOT2, wtot(nSurf), dummyReal, dummyReal)
      wtot(nSurf) = wtot(nSurf)* drymwt
    end if
    DO k=1,nmol
       IXoff=ImolG(k)-1
       !---Transform mix. ratios into mass fractions (relative to total air mass)
       do l=1,nSurf
          qr(l)=qcor(l)*xG(iXoff+l)
       end do
      call lzsum4tot(pLoc(1),xG(1),RHOTOT1)
      do l=1,nEnd
         call lzsum4tot(pLoc(l+1),xG(l+1),RHOTOT2)

         scal = (alt(l) - alt(l + 1)) *drymwt
         call lzsum4W(qr(l), qr(l + 1), scal, &
             RHOTOT1, RHOTOT2, wvar(k, l), dwqu(l, k), dwql(l, k))
         if (k == 1) THEN
            wdry(l) = wtot(l) - wvar(k, l)
            qvar(k, l) = wvar(k, l)/wtot(l)
         else
            qvar(k, l) = wvar(k, l)/wdry(l)
         end if
         RHOTOT1 = RHOTOT2
       end do
       if (interpSfc) then
          call lzsum4tot(pSfc,tSfc,RHOTOT2)
          qsfc=(qr(nSurf-1)**alpha)*(qr(nSurf)**(1.0-alpha))

         scal = dzsfc *drymwt
         call lzsum4W(qr(nSurf-1), qsfc, scal, &
                RHOTOT1, RHOTOT2, wvar(k, Nlast), dwqu(Nlast, k), dwql(Nlast, k))

         if (k == 1) THEN
            wdry(Nlast) = wtot(Nlast) - wvar(1, Nlast)
            qvar(1, Nlast) = wvar(1, Nlast)/wtot(Nlast)
         else
            qvar(k, Nlast) = wvar(k, Nlast) / wdry(Nlast)
         end if

         dwqu(Nlast, k) = dwqu(Nlast, k) + dwql(Nlast, k)*alpha*qsfc/qr(nSurf-1)
         dwql(Nlast, k) = dwql(Nlast, k)*(1.0 - alpha)*qsfc/qr(nSurf)

       end if
    END DO
    !---Derivatives of amount wrt dry mixing ratios
    DO n=1,Nlast
      dwqu(n,1)=dwqu(n,1)*qcor(n)**2
      dwql(n,1)=dwql(n,1)*qcor(n+1)**2
      dwqu(n,2:nMol)=dwqu(n,2:nMol)*qcor(n)
      dwql(n,2:nMol)=dwql(n,2:nMol)*qcor(n+1)
    END DO
    RETURN
  END SUBROUTINE fpathZ_mw

  !----------------------------------------------------------------------------
  ! PURPOSE: Computes temperature/Water vapor indexes and interpolation
  !          coefficients.
  !----------------------------------------------------------------------------
  SUBROUTINE settabindx(tavl,qh2o,pavl,N2,referenceGrid,indx_tmp_p1,indx_tmp_p2, &
       indxp,dtmp_p1,dtmp_p2,at_p1,at_p2,ap1,ap2,indx_wvp_p1,indx_wvp_p2)
    !---Input variables
    INTEGER,               INTENT(IN)    :: N2
    REAL,    DIMENSION(:), INTENT(IN)    :: tavl,qh2o,pavl
    LOGICAL,            INTENT(IN)    :: referenceGrid
    !---Output variables
    INTEGER, DIMENSION(:), INTENT(INOUT) :: indx_tmp_p1,indx_tmp_p2,indxp
    INTEGER, DIMENSION(:), INTENT(INOUT) :: indx_wvp_p1,indx_wvp_p2
    REAL,    DIMENSION(:), INTENT(INOUT) :: at_p1,dtmp_p1,ap1
    REAL,    DIMENSION(:), INTENT(INOUT) :: at_p2,dtmp_p2,ap2
    !---Local variables
    INTEGER                              :: i,l,lp1,lp2
    !---Compute temperature index and interpolation coefficient
    DO l=1,N2
       if (referenceGrid) then
          lp1=l
          ap2(l)=0.
          ap1(l)=1.
       else
          lp1=1
          DO WHILE(pavlref(lp1).LT.pavl(l))
             lp1=lp1+1
          ENDDO
          lp2=lp1
          lp1=lp1-1
          if (lp1 == 0) then
             lp1=1
             lp2=2
          endif
          if (lp2 .ge. numRefLev) then
             lp1=numRefLev-2
             lp2=numRefLev-1
          endif
          ap2(l)=(pavl(l)-pavlref(lp1))/(pavlref(lp2)-pavlref(lp1))
          ap1(l)=(pavlref(lp2)-pavl(l))/(pavlref(lp2)-pavlref(lp1))
       endif
       indxp(l)=lp1
       TempLoop1: DO i=2,ntmpod-1
          IF(tmptab(i,lp1).GE.tavl(l)) exit TempLoop1
       END DO TempLoop1
       indx_tmp_p1(l)=i-1
       dtmp_p1(l)=tmptab(i,lp1)-tmptab(i-1,lp1)
       at_p1(l)=(tavl(l)-tmptab(i-1,lp1))/dtmp_p1(l)
    !---Compute water vapor index
       WVLoop1: DO i=2,nwvpod-1
          IF(wvptab(i,lp1).GE.qh2o(l)) exit WVLoop1
       END DO WVLoop1
       indx_wvp_p1(l)=i-1
       if (.not.referenceGrid) then
          TempLoop2: DO i=2,ntmpod-1
             IF(tmptab(i,lp2).GE.tavl(l)) exit TempLoop2
          END DO TempLoop2
          indx_tmp_p2(l)=i-1
          dtmp_p2(l)=tmptab(i,lp2)-tmptab(i-1,lp2)
          at_p2(l)=(tavl(l)-tmptab(i-1,lp2))/dtmp_p2(l)
          WVLoop2: DO i=2,nwvpod-1
             IF(wvptab(i,lp2).GE.qh2o(l)) exit WVLoop2
          END DO WVLoop2
          indx_wvp_p2(l)=i-1
       endif
    ENDDO
    RETURN
  END SUBROUTINE settabindx

  !----------------------------------------------------------------------------
  ! PURPOSE:  Given the LUTs, the function computes its layer
  !           optical depths. Refer to ossrad_mw for the transmittance and
  !           brightness temperature calculation
  !----------------------------------------------------------------------------
  SUBROUTINE OSStran(kfix,dkfix,kh2o,dkh2o,kvar,dkvar,indx_tmp_p1, &
       indx_tmp_p2,indxp,dtmp_p1,dtmp_p2,at_p1,at_p2,ap1,ap2, &
       indx_wvp_p1,indx_wvp_p2,NmolS,ImolS,wdry,qvar,wvar,N2, &
       referenceGrid,tautot,dtaudtmp,dtaudmol)
    !---Input variables
    INTEGER,                                   INTENT(IN) :: NmolS
    INTEGER, DIMENSION(NmolS),                 INTENT(IN) :: ImolS
    REAL,    DIMENSION(NlayOD,NtmpOD),         INTENT(IN) :: kfix,dkfix
    REAL,    DIMENSION(NlayOD,NtmpOD,NwvpOD),  INTENT(IN) :: kh2o,dkh2o
    REAL,    DIMENSION(NmolS-1,NlayOD,NtmpOD), INTENT(IN) :: kvar,dkvar
    INTEGER, DIMENSION(:),    INTENT(IN)    :: indx_tmp_p1,indx_wvp_p1
    INTEGER, DIMENSION(:),    INTENT(IN)    :: indx_tmp_p2,indx_wvp_p2
    INTEGER, DIMENSION(:),    INTENT(IN)    :: indxp
    INTEGER,                  INTENT(IN)    :: n2
    REAL,    DIMENSION(:),    INTENT(IN)    :: dtmp_p1,dtmp_p2
    REAL,    DIMENSION(:),    INTENT(IN)    :: at_p1,at_p2,ap1,ap2,wdry
    REAL,    DIMENSION(:,:),  INTENT(IN)    :: qvar,wvar
    LOGICAL,            INTENT(IN)    :: referenceGrid

    !---Output variables
    REAL,    DIMENSION(:),    INTENT(INOUT) :: tautot,dtaudtmp
    REAL,    DIMENSION(:,:),  INTENT(INOUT) :: dtaudmol
    !---Local variables
    INTEGER:: l,indxw_p1,indxw_p2,indxt1_p1,indxt1_p2,indxt2_p1,indxt2_p2
    INTEGER:: indxp1,indxp2
    INTEGER:: kS,Imol
    REAL   :: qh2ol,wdryl,wh2ol
    REAL   :: atp1L,atp2L,ap1L,ap2L,dtmp1L,dtmp2L
    REAL   :: absfix1_p1,absfix1_p2,absfix2_p1,absfix2_p2,absfix_p1,absfix_p2
    REAL   :: absh2o1_p1,absh2o1_p2,absh2o2_p1,absh2o2_p2,absh2o_p1,absh2o_p2
    REAL   :: absvar1_p1,absvar1_p2,absvar2_p1,absvar2_p2,absvar_p2,absvar_p1
    REAL   :: dabsh2odq_p1,dabsh2odq_p2,dabsvardq_p1,dabsvardq_p2
    REAL   :: dabsfixdq_p1,dabsfixdq_p2,dtaudtmp_p1,dtaudtmp_p2
    DO l=1,N2
       !---Indexes and interpolation coeffs.
       qh2ol       = qvar(1,l)
       wdryl       = wdry(l)
       wh2ol       = wvar(1,l)
       indxt1_p1   = indx_tmp_p1(l)
       indxt2_p1   = indx_tmp_p1(l)+1
       indxw_p1    = indx_wvp_p1(l)
       indxp1      = indxp(l)
       atp1L       = at_p1(l)
       dtmp1L      = dtmp_p1(l)
       ap1L        = ap1(l)
       dtaudtmp_p1 = 0.
       !---Fixed gases
       absfix1_p1     = qh2ol*dkfix(indxp1,indxt1_p1) + kfix(indxp1,indxt1_p1)
       absfix2_p1     = qh2ol*dkfix(indxp1,indxt2_p1) + kfix(indxp1,indxt2_p1)
       absfix_p1      = atp1L*(absfix2_p1-absfix1_p1) + absfix1_p1
       dabsfixdq_p1   = atp1L*(dkfix(indxp1,indxt2_p1) - &
            dkfix(indxp1,indxt1_p1)) + dkfix(indxp1,indxt1_p1)
       !---Water vapor
       absh2o1_p1  = qh2ol*dkh2o(indxp1,indxt1_p1,indxw_p1) &
                    + kh2o(indxp1,indxt1_p1,indxw_p1)
       absh2o2_p1  = qh2ol*dkh2o(indxp1,indxt2_p1,indxw_p1) &
                    + kh2o(indxp1,indxt2_p1,indxw_p1)
       absh2o_p1   = atp1L*(absh2o2_p1-absh2o1_p1) + absh2o1_p1
       dabsh2odq_p1= atp1L*(dkh2o(indxp1,indxt2_p1,indxw_p1) - &
            dkh2o(indxp1,indxt1_p1,indxw_p1)) + dkh2o(indxp1,indxt1_p1,indxw_p1)
       dtaudtmp_p1 = ( (absfix2_p1 - absfix1_p1)*wdryl + &
            (absh2o2_p1 - absh2o1_p1)*wh2ol)*ap1L / dtmp1L

       if (referenceGrid) then
         absfix1_p2  =0.
         absfix2_p2  =0.
         absfix_p2   =0.
         absh2o1_p2  =0.
         absh2o2_p2  =0.
         absh2o_p2   =0.
         dabsfixdq_p2=0.
         dabsh2odq_p2=0.
         dtaudtmp_p2 =0.
         ap2L        =0.
       else
          indxt1_p2   = indx_tmp_p2(l)
          indxt2_p2   = indx_tmp_p2(l)+1
          indxw_p2    = indx_wvp_p2(l)
          indxp2      = indxp(l)+1
          atp2L       = at_p2(l)
          dtmp2L      = dtmp_p2(l)
          ap2L        = ap2(l)

          !---Fixed gases
          absfix1_p2  = qh2ol*dkfix(indxp2,indxt1_p2) + kfix(indxp2,indxt1_p2)
          absfix2_p2  = qh2ol*dkfix(indxp2,indxt2_p2) + kfix(indxp2,indxt2_p2)
          absfix_p2   = atp2L*(absfix2_p2-absfix1_p2) + absfix1_p2
          dabsfixdq_p2= atp2L*(dkfix(indxp2,indxt2_p2) - &
            dkfix(indxp2,indxt1_p2)) + dkfix(indxp2,indxt1_p2)
          !---Water vapor
          absh2o1_p2  = qh2ol*dkh2o(indxp2,indxt1_p2,indxw_p2) &
                      + kh2o(indxp2,indxt1_p2,indxw_p2)
          absh2o2_p2  = qh2ol*dkh2o(indxp2,indxt2_p2,indxw_p2) &
                      + kh2o(indxp2,indxt2_p2,indxw_p2)
          absh2o_p2   = atp2L*(absh2o2_p2-absh2o1_p2) + absh2o1_p2
          dabsh2odq_p2= atp2L*(dkh2o(indxp2,indxt2_p2,indxw_p2) - &
               dkh2o(indxp2,indxt1_p2,indxw_p2)) + &
               dkh2o(indxp2,indxt1_p2,indxw_p2)
          dtaudtmp_p2 = ( (absfix2_p2 - absfix1_p2)*wdryl + &
               (absh2o2_p2 - absh2o1_p2)*wh2ol)*ap2L / dtmp2L
       endif

       tautot(l)      = (absfix_p1*wdryl + absh2o_p1*wh2ol)*ap1L + &
            (absfix_p2*wdryl + absh2o_p2*wh2ol)*ap2L
       dtaudmol(1,l)  = (dabsfixdq_p1 - absfix_p1 + dabsh2odq_p1*qh2ol + &
            absh2o_p1)*ap1L + &
            (dabsfixdq_p2 - absfix_p2 + dabsh2odq_p2*qh2ol + &
            absh2o_p2)*ap2L
       dtaudtmp(l)    = ( (absfix2_p1 - absfix1_p1)*wdryl + &
            (absh2o2_p1 - absh2o1_p1)*wh2ol)*ap1L + &
            ( (absfix2_p2 - absfix1_p2)*wdryl + &
            (absh2o2_p2 - absh2o1_p2)*wh2ol)*ap2L
       !---Variable gases
       DO kS=2,nmolS
          Imol=ImolS(kS)
          absvar1_p1  = qh2ol*dkvar(kS-1,indxp1,indxt1_p1) + &
               kvar(kS-1,indxp1,indxt1_p1)
          absvar2_p1   = qh2ol*dkvar(kS-1,indxp1,indxt2_p1) + &
               kvar(kS-1,indxp1,indxt2_p1)
          absvar_p1   = atp1L*(absvar2_p1-absvar1_p1) + absvar1_p1
          dabsvardq_p1= atp1L*(dkvar(kS-1,indxp1,indxt2_p1) - &
               dkvar(kS-1,indxp1,indxt1_p1)) + dkvar(kS-1,indxp1,indxt1_p1)
          dtaudtmp_p1 = dtaudtmp_p1 + (absvar2_p1 - absvar1_p1)*&
               wvar(Imol,l)*ap1L / dtmp1L
          if (referenceGrid) then
            absvar1_p2  =0.
            absvar2_p2  =0.
            absvar_p2   =0.
            dabsvardq_p2=0.
          else
             absvar1_p2  = qh2ol*dkvar(kS-1,indxp2,indxt1_p2) + &
                  kvar(kS-1,indxp2,indxt1_p2)
             absvar2_p2  = qh2ol*dkvar(kS-1,indxp2,indxt2_p2) + &
                  kvar(kS-1,indxp2,indxt2_p2)
             absvar_p2   = atp2L*(absvar2_p2-absvar1_p2) + absvar1_p2
             dabsvardq_p2= atp2L*(dkvar(kS-1,indxp2,indxt2_p2) - &
                  dkvar(kS-1,indxp2,indxt1_p2)) + dkvar(kS-1,indxp2,indxt1_p2)
             dtaudtmp_p2 = dtaudtmp_p2 + (absvar2_p2 - absvar1_p2)*&
                  wvar(Imol,l) *ap2L / dtmp2L
          end if
          tautot(l)   = tautot(l) + (absvar_p1*ap1L + &
               absvar_p2*ap2L)*wvar(Imol,l)
          dtaudmol(1,l) = dtaudmol(1,l) + ((dabsvardq_p1  - absvar_p1)*ap1L + &
               (dabsvardq_p2  - absvar_p2)*ap2L)*qvar(Imol,l)
          dtaudmol(kS,l) = absvar_p1*ap1L + absvar_p2*ap2L
       END DO
       !---Scale derivative terms
       dtaudtmp(l) = dtaudtmp_p1 + dtaudtmp_p2
    END DO
    RETURN
  END SUBROUTINE OSStran


  REAL FUNCTION ABCOEFLIQ(FREQ,TEMP)
    !-------------------------------------------------------------------------
    !    COMPUTES ABSORPTION COEFFICIENT (m^2/kg) FOR SUSPENDED WATER DROPLETS
    !     FROM EQUATIONS OF LIEBE, HUFFORD AND MANABE
    !     (INT. J. IR & MM WAVES V.12(17) JULY 1991
    ! The code was written by Phil Rosenkranz 8/3/92.
    ! The alternative, exponential, form for FP was inserted by Alan Lipton
    ! 1/16/98, following K. F. Evans, and based on the above-referenced paper.
    ! The code was converted from a computation of absorption in nepers/km
    ! (WITH cloud water concentration as an input variable) to absorption
    ! coefficient by Alan Lipton 4/16/98.
    !     FREQ IN GHZ     (VALID FROM 0 TO 1000 GHZ)
    !     TEMP IN KELVIN
    !-------------------------------------------------------------------------
    !---Input variables
    REAL, INTENT(IN) :: FREQ,TEMP
    !---Local variables
    REAL             :: THETA1,EPS0,EPS1,EPS2,FP,FS
    COMPLEX          :: EPS,RE
    THETA1    = 1.-300./TEMP
    EPS0      = 77.66 - 103.3*THETA1
    EPS1      = .0671*EPS0
    EPS2      = 3.52 + 7.52*THETA1
    FP        = 20.1*EXP(7.88*THETA1)
    !FP       = (316.*THETA1 + 146.4)*THETA1 +20.20
    FS        = 39.8*FP
    EPS       = (EPS0-EPS1)/CMPLX(1.,FREQ/FP)+(EPS1-EPS2)/CMPLX(1.,FREQ/FS)+EPS2
    RE        = (EPS-1.)/(EPS+2.)
    ABCOEFLIQ = -.06286*AIMAG(RE)*FREQ
    RETURN
  END FUNCTION ABCOEFLIQ

!--
  FUNCTION Get_Lun() RESULT( Lun )
    INTEGER :: Lun
    Lun = 9
    Lun_Search: DO
      ! -- Increment logical unit number
      Lun = Lun + 1
      ! -- If file unit does not exist, set to -1 and exit
      IF ( .NOT. File_Unit_Exists( Lun ) ) THEN
        Lun = -1
        EXIT Lun_Search
        end if
      ! -- If the file is not open, we're done.
      IF ( .NOT. File_Open_by_Unit( Lun ) ) EXIT Lun_Search
    END DO Lun_Search
  END FUNCTION Get_Lun
! Helper function
! FindFreeUnit
! return a free UNIT
integer function findFreeUnit()
    logical isOpened
    integer iostat
    DO findFreeUnit = 10, 200
            inquire(UNIT = FindFreeUnit, OPENED = isOpened, iostat = iostat)
            if (iostat .ne. 0 ) cycle
            if (.not. isOpened) return
    END DO
    STOP '[FindFreeUnit] No Free File Unit Found'
end function findFreeUnit
!--
  FUNCTION File_Unit_Exists( FileID ) RESULT ( Existence )
    INTEGER, INTENT( IN ) :: FileID
    LOGICAL :: Existence
    INQUIRE( UNIT = FileID, EXIST = Existence )
  END FUNCTION File_Unit_Exists
!--
  FUNCTION File_Open_by_Unit( FileID ) RESULT ( Is_Open )
    INTEGER, INTENT( IN ) :: FileID
    LOGICAL :: Is_Open
    INQUIRE( UNIT = FileID, OPENED = Is_Open )
  END FUNCTION File_Open_by_Unit
!--
!--
  subroutine oss_destroy_base()
    if (allocated(molIDfix))      deallocate(molIDfix)
    if (allocated(molID))         deallocate(molID)
    if (allocated(ImolG))         deallocate(ImolG)
    if (allocated(ichmap))        deallocate(ichmap)
    if (allocated(nch))           deallocate(nch)
    if (allocated(coef))          deallocate(coef)
    if (allocated(coef1))         deallocate(coef1)
    if (allocated(coef2))         deallocate(coef2)
    if (allocated(my_cFreq))      deallocate(my_cFreq)
    if (allocated(my_ipol))       deallocate(my_ipol)
    if (allocated(vFreq))         deallocate(vFreq)
    if (allocated(cosmbk))        deallocate(cosmbk)
    if (allocated(fixtab))        deallocate(fixtab)
    if (allocated(varDflt))       deallocate(varDflt)
    if (allocated(pref))          deallocate(pref)
    if (allocated(NmolS))         deallocate(NmolS)
    if (allocated(ImolS))         deallocate(ImolS)
    if (allocated(Tmptab))        deallocate(Tmptab)
    if (allocated(Wvptab))        deallocate(Wvptab)
    if (allocated(kfix))          deallocate(kfix)
    if (allocated(dkfix))         deallocate(dkfix)
    if (allocated(kh2o))          deallocate(kh2o)
    if (allocated(dkh2o))         deallocate(dkh2o)
    if (allocated(kvar))          deallocate(kvar)
    if (allocated(dkvar))         deallocate(dkvar)
    if (allocated(my_side1))      deallocate(my_side1)
    if (allocated(my_side2))      deallocate(my_side2)
    if (allocated(my_hwhmx))      deallocate(my_hwhmx)
    if (allocated(my_chanID))     deallocate(my_chanID)
    use2ndOrderPlanck=.FALSE.
    myRadOrTb=0
    loadOSSselMWdone=.false.
    loadODdone=.false.
  end subroutine oss_destroy_base
END MODULE OSSMWmoduleSubs


