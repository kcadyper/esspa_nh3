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

!------------------------------------------------------------
!
!  MODULE OSS_ADDBL: contains items needed for the
!                OSS Scattering Adding-Doubling Module
!
!--------------------------------------------------------------
 MODULE oss_addbl

! <f90Module>***********************************************************
!
! NAME:
!
!   oss_addbl
!
! PURPOSE:
!
!   Contains items needed for the OSS scattering adding-doubling radiative
!   transfer.
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

  USE type_kinds, ONLY: Double, FP
  use params_module, ONLY: mxcang, maxcmu, mxang, mxang2, accur
  use subs_oss_addbl
  use oss_addbl_a
  implicit none

  PRIVATE

  public :: scattRT, ossscat

CONTAINS

  SUBROUTINE scattRT(tautotIn,tempIn,btempIn,tabsHydrIn,tscattIn,gasymIn,&
       emisIn,solReflIn, cosmicBackgroundIn,fBeamIn,vnIn,nLay,nObs,&
       sBeam,nStr,umuIn,umu0In,lookup,delphiIn,iflux,flxNetMonoIn,radOut,& 
       pland,linInTauFlag, lambertian, btempOut_a,emisOut_a,tempOut_a, &
       tautotOut_a,tabshydrOut_a,tscattOut_a,gasymOut_a)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   scattRT
!
! PURPOSE:
!
!   Wrapper to interface oss_module with ossScat.
!
! SYNTAX:
!
!   CALL scattRT(tautot,temp,btemp,tabsHydr,tscatt,gasym,emis,solRefl, &
!       cosmicBackground,fBeam,vn,nLay,nObs,&
!       sBeam,nStr,umu,umu0,lookup,delphi,iflux,flxNetMono,rad, &
!       linInTauFlag, btemp_a,emis_a,temp_a,tautot_a,tabshydr_a,tscatt_a,gasym_a,pland)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   tautot      REAL(FP)     Optical depth in layers: trace gases contribution
!   temp        REAL(FP)     Temperature
!   btemp       REAL(FP)     Surface skin temperature
!   tabsHydr    REAL(FP)     In-layer hydrometeor contribution to absorption
!                        optical depth
!   tscatt      REAL(FP)     In-layer cloud contribution to scattering optical
!                        depth
!   gasym       REAL(FP)     Asymmetry parameter
!   emis        REAL(FP)     Emissivity
!   solRefl     REAL(FP)     Surface reflectivity for solar radiation
!   vn          REAL(FP)     Wavenumber, units of cm^-1
!   nLay        INTEGER  Number of atmospheric layers
!   nObs        INTEGER  Index of atmospheric level of the observer
!   sBeam       LOGICAL  Flag for solar calculations
!   nStr        INTEGER  Number of streams
!   umu         REAL(FP)     Cosine of viewing path angle
!   umu0        REAL(FP)     Cosine of solar path angle
!   lookup      LOGICAL  Flag for up-looking instruments
!   delphi      REAL(FP)     Azimuthal angle offset
!   iflux       INTEGER  TBD
!   flxNetMono  REAL(FP)     TBD
!   cosmicBackground REAL(FP) cosmic background radiation
!   pland       REAL         Fraction of land
!   linInTauFlag  LOGICAL    IF LINER-IN-TAU SHOULBE BE USED
!
!   INPUTS/OUTPUTS:
!
!   fBeam       REAL(FP)     Solar beam
!   rad         REAL(FP)     Radiance
!   btemp_a, emis_a, temp_a, tautot_a, tabshydr_a, tscatt_a, gasym_a
!               REAL(FP)     derivative arrays
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    implicit none
! input
    INTEGER,               INTENT(IN)    :: nLay,nObs
    REAL,                  INTENT(IN)    :: tautotIn(nlay),tempIn(nlay+1)
    REAL,                  INTENT(IN)    :: btempIn
    REAL,                  INTENT(IN)    :: tabsHydrIn(nlay),tscattIn(nlay),gasymIn(nlay)
    REAL,                  INTENT(IN)    :: emisIn,solReflIn
    REAL,                  INTENT(IN)    :: fBeamIn, cosmicBackgroundIn
    REAL,                  INTENT(IN)    :: vnIn
    LOGICAL,               INTENT(IN)    :: sbeam
    INTEGER,               INTENT(IN)    :: nStr
    REAL,                  INTENT(IN)    :: umuIn,umu0In
    LOGICAL,               INTENT(IN)    :: lookup
    REAL,                  INTENT(IN)    :: delphiIn
    INTEGER,               INTENT(IN)    :: iflux
    REAL,                  INTENT(IN)    :: flxNetMonoIn(nlay)
    REAL,        OPTIONAL, INTENT(IN)    :: pland
    LOGICAL,     OPTIONAL, INTENT(IN)    :: linInTauFlag
    logical,     optional, intent(IN)    :: lambertian

    REAL,                  INTENT(OUT)   :: radOut
    real,    OPTIONAL,     INTENT(INOUT) :: btempOut_a
    real,    OPTIONAL,     INTENT(INOUT) :: tempOut_a(nlay+1)
    real,    OPTIONAL,     INTENT(INOUT) :: tautotOut_a(nlay), tabshydrOut_a(nlay), gasymOut_a(nlay), emisOut_a
    real,    OPTIONAL,     INTENT(INOUT) :: tscattOut_a(nlay)


    !---Wrapper to interface oss_module with ossScat.  Some calling arguments
    !      are placeholders.

    !--Input/output variables
    REAL(FP)    :: tautot(nlay),temp(nlay+1)
    REAL(FP)    :: btemp
    REAL(FP)    :: tabsHydr(nlay),tscatt(nlay),gasym(nlay)
    REAL(FP)    :: emis,solRefl
    REAL(FP)    :: fBeam
    REAL(Double):: vn
    REAL(FP)    :: umu,umu0
    REAL(FP)    :: delphi
    REAL(FP)    :: cosmicBackground
    REAL(FP)    :: flxNetMono(nlay)
    REAL(FP)    :: rad


    !---Internal variables
    INTEGER    :: lTyp(nlay)
    REAL(FP)                   :: tau(nlay), albedo
    REAL(FP)                   :: hl(0:maxcmu),pmom(0:maxcmu,nlay)
    REAL(FP),    PARAMETER     :: scatLB=1e-8, missingReal(Fp)=-999.
    INTEGER                :: l,k,ik,ikm1
    LOGICAL, parameter     :: deltam = .true.
    LOGICAL, save          :: lamber
    LOGICAL                :: plank
    LOGICAL     :: linInTauFlagLoc

!derivative arrays: d(rad)/dx
! for subroutine arguments:
    real(fp) :: rad_a !(drad/drad=1)
    real(fp) :: btemp_a
    real(fp) :: temp_a(nlay+1)
    real(fp) :: tautot_a(nlay), tabshydr_a(nlay), gasym_a(nlay), emis_a
    real(fp) :: tscatt_a(nlay)
! local variables (these are mapped to the above by pre_scattrt_a):
    real(fp) :: albedo_a, dummyReal
    real(fp) :: tau_a(nlay), pmom_a(0:maxcmu,nlay)

    !---Variable intialization
    tautot = tautotIn
    temp = tempIn
    btemp = btempIn
    tabsHydr = tabsHydrIn
    gasym = gasymIn
    emis = emisIn
    solRefl = solReflIn
    fBeam = fBeamIn
    cosmicBackground = cosmicBackgroundIn
    vn =vnIn
    umu = umuIn
    umu0 = umu0In
    delphi = delphiIn
    flxNetMono = flxNetMonoIn

    IF (present(lambertian)) THEN
       lamber=lambertian
    else
       lamber=.true.
    ENDIF

    plank=.false.
    if(vn.lt.5000.) plank=.true.

    if(.not.sbeam)fBeam=0.
    if (present(linInTauFlag)) then
    	linInTauFlagLoc = linInTauFlag
    else
    	linInTauFlagLoc = .false.
    end if

    if (nstr .gt. maxcmu) &
              call errmsg("nstr .gt. maxcmu in scattRT",.TRUE.)

    tscatt = tscattIn
    albedo = solRefl
    call pre_scattrt(tscatt,scatlb,tautot,tabsHydr,gasym,ltyp,tau,pmom,nlay,nstr)

    call ossscat( nlay, tau, tscatt, pmom, temp, lTyp, &
         vn, nstr, nobs, umu, fBeam, umu0, delphi, lamber, emis, albedo, hl, &
         btemp, deltam, plank, sbeam, lookup, rad, cosmicBackground, linInTauFlagLoc)
   radOut = rad

    IF (.not.(present(btempOut_a) .or. present(emisOut_a) .or.  present(tempOut_a) .or. &
              present(tautotOut_a) .or. present(tabshydrOut_a) .or. present(tscattOut_a) &
              .or.  present(gasymOut_a))) RETURN

    ! Initialize ADJ (Note: tautot_a, tabshydr_a, and gasym_a are initialized in pre_scattrt_a)
    tau_a(:)=0.
    tscatt_a(:)=0.
    pmom_a(:,:)=0.
    btemp_a=0.
    albedo_a=0.
    temp_a(:)=0.
    rad_a=1.0 !(drad/drad=1)
    tscatt = tscattIn
    albedo = solRefl
    call pre_scattrt(tscatt,scatlb,tautot,tabsHydr,gasym,ltyp,tau,pmom,nlay,nstr)
    call ossscat_a( nlay, tau,tau_a,  tscatt, tscatt_a, pmom, pmom_a, temp, temp_a, lTyp, &
         vn, nstr, nobs, umu, fBeam, umu0, delphi, lamber, emis, albedo, albedo_a, hl, &
         btemp, btemp_a, deltam, plank, sbeam, lookup, dummyReal, rad_a, linInTauFlagLoc)

    call pre_scattrt_a(tscatt, scatlb, tautot, tautot_a, tabshydr, tabshydr_a &
         , gasym, gasym_a, emis_a, ltyp, tau, tau_a, pmom, pmom_a &
         , albedo_a, nlay, nstr)

    btempOut_a = btemp_a
    emisOut_a = emis_a
    tempOut_a = temp_a
    tautotOut_a = tautot_a
    tabshydrOut_a = tabshydr_a
    tscattOut_a = tscatt_a
    gasymOut_a = gasym_a

    RETURN

  END SUBROUTINE SCATTRT

      SUBROUTINE OSSSCAT( NLAY, TAUABS, TAUSCAT, ALPHA, TEMP, LTYPIN, &
         VN, NSTR, LOBS, UMU, FBEAM, UMU0, DELPHI, LAMBER, emissivity, ALBEDO, HL, &
         TSFC, DELTAM, PLANKIN, SBEAMIN, LOOKUP, RADOUT, COSMICBACKGROUND, linInTauFlag)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   OSSSCAT
!
! PURPOSE:
!
!   Driver for scattering adding-doubling radiative transfer.
!
! SYNTAX:
!
!   CALL OSSSCAT(NLAY, TAUABS, TAUSCAT, ALPHA, TEMP, LTYPIN, VN,
!      NSTR, LOBS, UMU, FBEAM, UMU0, DELPHI, LAMBER, ALBEDO, HL,
!      TSFC, DELTAM, PLANKIN, SBEAMIN, LOOKUP, RADOUT)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   NLAY     INTEGER  Number of atmospheric layers
!   TAUABS   REAL(FP)     Array of cloud contribution to absorption optical
!                     depths
!   TEMP     REAL(FP)     Temperature
!   LTYPIN   INTEGER  Flag for pointing to cloud active layers
!   VN       REAL(Double)   Wavenumber, units of cm^-1
!   NSTR     INTEGER  Number of streams
!   LOBS     INTEGER  Observer level
!   UMU      REAL(FP)     Cosine of viewing path angle
!   FBEAM    REAL(FP)     Solar beam
!   UMU0     REAL(FP)     Cosine of solar path angle
!   DELPHI   REAL(FP)     Azimuthal angle offset
!   LAMBER   LOGICAL  Flag for Lambertian surface
!   ALBEDO   REAL(FP)     Albedo
!   HL       REAL(FP)     TBD
!   TSFC     REAL(FP)     Surface-level air temperature
!   DELTAM   LOGICAL  Flag for Delta-M scaling
!   PLANKIN  LOGICAL  Flag for applying Planck function
!   SBEAMIN  LOGICAL  Flag for solar calculations
!   LOOKUP   LOGICAL  Flag for up-looking instruments
!   COSMICBACKGROUND REAL(FP) cosmic background radiation
!   linInTauFlag  LOGICAL    IF LINER-IN-TAU SHOULBE BE USED
!
!   INPUTS/OUTPUTS:
!
!   TAUSCAT  REAL(FP)     Array of cloud contribution to scattering optical
!                     depths
!   ALPHA    REAL(FP)     Phase function
!   RADOUT   REAL(FP)     Radiance for output
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! **** ***********   VERSION 2.0 :    21 MARCH 1995   *****************
!
!                       ALGORITHM:    J.L. MONCET
!                                     S.A. CLOUGH
!
!                  IMPLEMENTATION:    J.L. MONCET
!                                     D.A. PORTMAN
!
!               ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!               840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139
!
!----------------------------------------------------------------------
!
!               WORK SUPPORTED BY:    THE ARM PROGRAM
!                                     OFFICE OF ENERGY RESEARCH
!                                     DEPARTMENT OF ENERGY
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
      INTEGER,                            INTENT(IN)    :: NLAY,NSTR
      REAL(Double),                             INTENT(IN)    :: VN
      REAL(FP),   DIMENSION(NLAY),           INTENT(IN)    :: TAUABS
      REAL(FP),   DIMENSION(NLAY),           INTENT(INOUT) :: TAUSCAT
      REAL(FP),   DIMENSION(0:MAXCMU,NLAY),INTENT(INOUT) :: ALPHA
      REAL(FP),   DIMENSION(0:NLAY),         INTENT(IN)    :: TEMP
!---- REAL(FP) UMU,UMU0,DELPHI,FBEAM,ACCUR
      REAL(FP),                               INTENT(IN)    :: UMU,UMU0
      REAL(FP),                               INTENT(IN)    :: DELPHI,FBEAM
      REAL(FP),                               INTENT(IN)    :: ALBEDO,HL(1)
      REAL(FP),                               INTENT(IN)    :: TSFC
      REAL(FP),                               INTENT(INOUT) :: RADOUT
      INTEGER,DIMENSION(NLAY),           INTENT(IN)    :: LTYPIN
      INTEGER,                            INTENT(IN)    :: LOBS
      LOGICAL,                            INTENT(IN)    :: DELTAM,PLANKIN
      LOGICAL,                            INTENT(IN)    :: SBEAMIN,LOOKUP
      LOGICAL,                            INTENT(IN)    :: LAMBER
      REAL(FP),                           INTENT(IN)    :: cosmicBackground
      LOGICAL                                           :: linInTauFlag
      REAL(FP),                           INTENT(IN)    :: emissivity
      integer                                           :: from
!
      REAL(FP) TAUEXT(NLAY),SSALB(NLAY),F0(0:NLAY)
      real(fp) alpha1(0:maxcmu),rad(mxang)
      LOGICAL REFLB,REFLT
      LOGICAL PLANK,SBEAM
      INTEGER LTYP(NLAY)
      REAL(FP) :: azerr,azterm, cl2inv,embot,emml,empl,emtop,f,gmu,gmu0,gmu0vec,gwt
      integer  :: i,init,isrfl,k,kconv,l,lnum,maz,n1,n2,n2ang,nang,nang2,naz,nbot
      integer  :: nbot0,ncang,ntop,ntop0
      REAL(FP) :: pi,pth,rbot,rlay,rtop,sumex,tlay,ylm,ylm0
!
      REAL(FP) :: TTOP(MXANG2),TBOT(MXANG2)
      dimension RLAY(MXANG2),TLAY(MXANG2), EMPL(MXANG),EMML(MXANG)
      dimension RTOP(MXANG2),EMTOP(MXANG)
      dimension RBOT(MXANG2),EMBOT(MXANG)

      COMMON /ANGL  / PTH(MXANG),GMU(MXANG),GWT(MXANG),GMU0,&
                     YLM(0:2*MXCANG-1,MXANG),YLM0(0:2*MXCANG-1),&
                     GMU0VEC(1)
      save /angl/
      COMMON /LPAR  / PLANK,SBEAM
      save /lpar/
      COMMON /PNLPAR/ NANG,N2ANG,NANG2
      save /pnlpar/
      COMMON /CONSTN/ PI,CL2INV
      save /constn/

!
      DATA INIT/1/
      SAVE INIT
!
      PLANK=PLANKIN
      SBEAM=SBEAMIN
      LTYP(1:NLAY)=LTYPIN(1:NLAY)
!
      IF(INIT.EQ.1)THEN
! ======================================================================
! === INITIALIZATION ===================================================
! ======================================================================
         PI=ACOS(-1.)
         CL2INV=1./LOG(2.)
!
         NCANG=NSTR/2
!
!-----Generate discrete angles and associated weights for
!-----Double-Gauss" quadrature (see Stamnes, 1981)
!
         CALL QGAUSN(NCANG,GMU,GWT)
!
         NANG=NCANG+1
         N2ANG=2*NANG
         NANG2=NANG*NANG
!
         DO 120 I=1,NCANG
  120    PTH(I)=1./GMU(I)
!
         INIT=0
      END IF
!
! *** Insert user angle.
!
      GMU(NANG)=UMU
      PTH(NANG)=1./UMU
!
! *** Set NAZ
!
      NAZ=0
      IF(SBEAM)THEN
         GMU0=UMU0
         IF((1.-GMU0).GT.1.E-05.AND. (1.-GMU(NANG)).GT.1.E-05)NAZ=N2ANG-3
      END IF
!
! *** Set control parameters for CHARTS run
!
      IF(LAMBER) THEN
         ISRFL=-1
         IF(ALBEDO.LE.1.E-6)ISRFL=0
      ELSE
         ISRFL=2
         IF(ALBEDO.LE.1.E-6)ISRFL=0
         !STOP 'OSSSCAT: This version only supports Lambertian surfaces'
      ENDIF
!
      CALL SETCTL(NLAY,LTYP,LOBS,SBEAM,PLANK,ISRFL,LOOKUP,NAZ,&
         NBOT0,NBOT,NTOP0,NTOP)
!
      N1=NTOP0
      N2=NBOT0
!
! *** Delta-M scaling
!
      F=0.
      DO I=1,NLAY
         IF(DELTAM)THEN
            F=ALPHA(N2ANG-2,I)
            TAUSCAT(I)=TAUSCAT(I)*(1.-F)
         END IF
!
         DO 50 L=0,N2ANG-3
            ALPHA(L,I)=(2*L+1)*(ALPHA(L,I)-F)/(1.-F)
 50      CONTINUE
      ENDDO

      DO 330 L=1,NLAY
         TAUEXT(L)=TAUABS(L)
         IF(LTYP(L).GT.0)THEN
            TAUEXT(L)=TAUEXT(L)+TAUSCAT(L)
            IF(TAUEXT(L).LE.0.)THEN
               SSALB(L)=0.
            ELSE
               SSALB(L)=TAUSCAT(L)/TAUEXT(L)
            END IF
         END IF
  330 CONTINUE
!
      IF(SBEAM)THEN
!
!-----Compute solar beam extinction (current version does not treat
!-----specularly reflected beam)
!
         SUMEX=0.
         F0(0)=FBEAM
         DO 350 L=1,NLAY
!-----Add scattering component to total extinction
            SUMEX=SUMEX+TAUEXT(L)
!-----Calculate solar irradiance at top of scattering layer
            F0(L)=FBEAM*EXP(-SUMEX/GMU0)
  350    CONTINUE
      END IF
!
! =====================================================================
! === BEGIN LOOP OVER AZIMUTHAL ANGLES ================================
! =====================================================================
      KCONV=0
      DO 700 MAZ=0,NAZ
       EMBOT=0.
       EMTOP=0.
       RTOP=0.
       RBOT=0.
       TTOP=0.
       TBOT=0.
       do k=1, nang2, nang+1
          TTOP(k) = 1.0
          TBOT(k) = 1.0
       end do
!
! *** Generate normalized Legendre polynomials for computational angles
! *** (GMU) and for incident solar beam angle cosine (GMU0).
!
      CALL LEPOLY_OSS(NANG,MAZ,2*MXCANG-1,N2ANG-3,GMU,YLM)
!     --Put GMU0 in vector form to prevent compiler error:
      GMU0VEC(1)=GMU0
      IF(SBEAM)CALL LEPOLY_OSS(1,MAZ,2*MXCANG-1,N2ANG-3,GMU0VEC,YLM0)
!      IF(SBEAM)CALL LEPOLY_OSS(1,MAZ,2*MXCANG-1,N2ANG-3,GMU0,YLM0)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     LOOP OVER LAYERS THE ABOVE OBSERVER LEVEL - MERGE DOWN CASE      +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REFLT=.FALSE.
!
      DO 500 LNUM=N1,LOBS-1
!
         IF(LTYP(LNUM).LE.0)THEN
!-----------------------------------------------------------------------
!     CLEAR LAYER CASE:
!-----------------------------------------------------------------------
! *** Accumulate partial results for successive clear layers "LAYDAT".
! *** Final merge with results contained in "ACMDWN" done at the end,
! *** when (LTYP=-99 OR -100).
            IF (LNUM .eq. N1) THEN
               CALL CLRLAY(TAUEXT(LNUM),TEMP(LNUM-1),&
                    TEMP(LNUM),REFLT,LTYP(LNUM),VN,tlay,empl,emml,cosmicBackground,linInTauFlag)
            ELSE
               CALL CLRLAY(TAUEXT(LNUM),TEMP(LNUM-1),&
                    TEMP(LNUM),REFLT,LTYP(LNUM),VN,tlay,empl,emml,0._fp,linInTauFlag)
           end if
            IF(LTYP(LNUM).LE.-99) then
               CALL ADDCLR(RTOP,TTOP,EMTOP,REFLT,PLANK, &
                    tlay, empl, emml)
            END IF
         ELSE
!-----------------------------------------------------------------------
!     CLOUDY LAYER CASE:
!-----------------------------------------------------------------------
! *** Calculate reflection/transmission matrices (in RLAY and TLAY) and
! *** internal sources (EMPL-upward- and EMML-downward) for current layer
!
            do k=0,n2ang-2
               alpha1(k)=alpha(k,lnum)
            end do
            CALL MSLAY(TAUEXT(LNUM),SSALB(LNUM),TEMP(LNUM-1),TEMP(LNUM),&
               F0(LNUM-1),RLAY,TLAY,EMPL,EMML,MAZ,VN,ALPHA1,0, linInTauFlag)
!
! *** Accumulate results in ACMWDN
!
            CALL ADDMS(RTOP,TTOP,EMTOP,REFLT,rlay,tlay,empl,emml)
         END IF
  500 CONTINUE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     LOOP OVER LAYERS BELOW THE OBSERVER LEVEL - MERGE UP CASE        +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! *** Initialize RBOT and EMBOT
!
      REFLB=.FALSE.
      CALL SURFCE(MAZ,ISRFL,TSFC,ALBEDO,F0(NLAY),RBOT,EMBOT,REFLB,VN, emissivity)
!
      if (maz .eq. 0) then
           from = nbot0
      else
           from = nbot
      end if

      DO 600 LNUM=from,LOBS,-1
         IF(LTYP(LNUM).LE.0)THEN
!-----------------------------------------------------------------------
!     CLEAR LAYER CASE:
!-----------------------------------------------------------------------
! *** Accumulate partial results for successive clear layers "LAYDAT".
! *** Final merge with results contained in "ACMUP" done at the end,
! *** when (LTYP=-99 OR -100).
            IF ( LNUM .eq. N1) THEN
                CALL CLRLAY(TAUEXT(LNUM),TEMP(LNUM),TEMP(LNUM-1), &
                     REFLB,LTYP(LNUM),VN,tlay, empl, emml, cosmicBackground, linInTauFlag)
            ELSE
                CALL CLRLAY(TAUEXT(LNUM),TEMP(LNUM),TEMP(LNUM-1), &
                     REFLB,LTYP(LNUM),VN,tlay, empl, emml, 0._fp, linInTauFlag)
            END IF
            IF(LTYP(LNUM).LE.-99) then
               CALL ADDCLR(RBOT,TBOT,EMBOT,REFLB,PLANK, &
                    tlay, empl, emml)
            END IF
         ELSE
!-----------------------------------------------------------------------
!     CLOUDY LAYER CASE:
!-----------------------------------------------------------------------
! *** Calculate reflection/transmission matrices (in RLAY and TLAY) and
! *** internal sources (EMPL-upward- and EMML-downward) for current layer
!
            do k=0,n2ang-2
               alpha1(k)=alpha(k,lnum)
            end do

            CALL MSLAY(TAUEXT(LNUM),SSALB(LNUM),TEMP(LNUM),TEMP(LNUM-1),&
               F0(LNUM-1),RLAY,TLAY,EMPL,EMML,MAZ,VN,ALPHA1,1, linInTauFlag)
!
! *** Accumulate results.
!
            CALL ADDMS(RBOT,TBOT,EMBOT,REFLB,rlay,tlay,empl,emml)
         END IF
  600 CONTINUE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     DO FINAL MERGE  - Output radiances in RAD                        +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (LOOKUP) THEN
         CALL MRGRAD(RTOP,EMTOP,RBOT,EMBOT,REFLT,REFLB,RAD)
      ELSE
         CALL MRGRAD(RBOT,EMBOT,RTOP,EMTOP,REFLB,REFLT,RAD)
      END IF
! =====================================================================
! === END OF PANEL LOOP ===============================================
! =====================================================================
      IF(MAZ.EQ.0)THEN
         RADOUT=RAD(NANG)
! *** (Reset control parameters for MAZ >0)
         PLANK=.FALSE.
         N1=NBOT
         N2=NTOP
      ELSE
! *** Compute radiance scaling for azimuthal order MAZ
         AZTERM=RAD(NANG)*COS(MAZ*DELPHI*PI/180.)
         RADOUT=RADOUT+AZTERM
!
         IF(RADOUT.NE.0.)THEN
            AZERR=ABS(AZTERM)/ABS(RADOUT)
            IF(AZERR.LE.ACCUR)KCONV=KCONV+1
         END IF
!
         IF (KCONV.GE.2) GOTO 800
!
      END IF
  700 CONTINUE
! =====================================================================
! === END OF LOOP OVER AZIMUTHAL ANGLES ===============================
! =====================================================================
  800 RETURN
!
      END subroutine ossscat
!
      subroutine pre_scattrt(tscatt,scatlb,tautot,tabsHydr,gasym,ltyp,tau,pmom,nlay,nstr)

        implicit none

        integer, intent(in) :: nlay, nstr
        REAL(FP),                  INTENT(IN)    :: tautot(nlay)
        REAL(FP),                  INTENT(IN)    :: tabsHydr(nlay),tscatt(nlay),gasym(nlay)
        REAL(FP),                  INTENT(IN)    :: scatlb

        integer, intent(out)                      :: ltyp(nlay)
        REAL(FP),                  INTENT(OUT)    :: tau(nlay)
        REAL(FP),intent(out) :: pmom(0:maxcmu,nlay)

        integer :: k,l

        pmom     = 0.
        lTyp     = 0

        !---Compute optical properties of default layers:
        do l=1,nLay
           !---Absorption optical depth (water vapor, fixed gases and hydrometeor):
           tau(l) =tautot(l) + tabsHydr(l)
           if (tscatt(l).gt.scatLB) lTyp(l)=1
           if (lTyp(l).gt.0.) then
              !---Phase function:
              pmom(0,l)=1.
              do k=1,nstr
                 pmom(k,l)=pmom(k-1,l)*gasym(l)
              end do
           end if
        end do
      end subroutine pre_scattrt

      SUBROUTINE PRE_SCATTRT_A(tscatt, scatlb, tautot, tautot_a, tabshydr, tabshydr_a &
           , gasym, gasym_a, emis_a, ltyp, tau, tau_a, pmom, pmom_a &
           , albedo_a, nlay, nstr)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: nlay, nstr
        REAL(fp) :: tautot(nlay)
        REAL(fp) :: tautot_a(nlay)
        REAL(fp) :: tabshydr(nlay), tscatt(nlay), gasym(nlay)
        REAL(fp) :: tabshydr_a(nlay), gasym_a(nlay)
        REAL(fp) :: scatlb
        REAL(fp) :: emis_a
        INTEGER, INTENT(OUT) :: ltyp(nlay)
        REAL(fp) :: tau(nlay)
        REAL(fp) :: tau_a(nlay)
        REAL(fp) :: pmom(0:maxcmu, nlay)
        REAL(fp) :: pmom_a(0:maxcmu, nlay)
        REAL(fp) :: albedo_a
        INTEGER :: k, l
        pmom = 0.
        ltyp = 0
        !---Compute optical properties of default layers:
        DO l=1,nlay
           !---Absorption optical depth (water vapor, fixed gases and hydrometeor):
           IF (tscatt(l) .GT. scatlb) ltyp(l) = 1
           IF (ltyp(l) .GT. 0.) THEN
              !---Phase function:
              pmom(0, l) = 1.
              DO k=1,nstr
                 pmom(k, l) = pmom(k-1, l)*gasym(l)
              END DO
           END IF
        END DO
        tautot_a(:) = 0.0
        tabshydr_a(:) = 0.0
        gasym_a(:) = 0.0
        DO l=nlay,1,-1
           IF (ltyp(l) .GT. 0.) THEN
              DO k=nstr,1,-1
                 pmom_a(k-1, l) = pmom_a(k-1, l) + gasym(l)*pmom_a(k, l)
                 gasym_a(l) = gasym_a(l) + pmom(k-1, l)*pmom_a(k, l)
                 pmom_a(k, l) = 0.0
              END DO
              pmom_a(0, l) = 0.0
           END IF
           tautot_a(l) = tautot_a(l) + tau_a(l)
           tabshydr_a(l) = tabshydr_a(l) + tau_a(l)
           tau_a(l) = 0.0
        END DO
        emis_a = -albedo_a
        albedo_a = 0.0
        pmom_a(0:maxcmu, :) = 0.0
      END SUBROUTINE PRE_SCATTRT_A

  end module oss_addbl
