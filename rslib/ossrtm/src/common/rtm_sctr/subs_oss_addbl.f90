module subs_oss_addbl

  use type_kinds, ONLY: Double, FP
  use params_module, ONLY: mxcang, maxcmu, mxang, mxang2, accur, one_dp
  implicit none
  public
    ! stands in numerator in dependence of radiance on TB,
    ! units: (mW/(m^2 ster cm^-1)) / [cm^-1)^3, scaled by 1e9
    real, parameter  :: RADCN1 = 1.191042722E-12* 1e7
    real, parameter  :: RADCN2 = 1.4387752

contains

  SUBROUTINE SETCTL(nlay, ltyp, lobs, sbeam, plank, isrfl, lookup, naz, &
&    nbot0, nbot, ntop0, ntop)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   SETCTL
!
! PURPOSE:
!
!   Set lowest and highest atmospheric layers for starting radiative transfer
!   calculations
!
! SYNTAX:
!
!   CALL SETCTL(NLAY, LTYP, LOBS, SBEAM, PLANK, ISRFL, LOOKUP, NAZ,
!      NBOT0, NBOT, NTOP0, NTOP)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   NLAY    INTEGER  Number of atmospheric layers
!   LOBS    INTEGER  Observer level
!   PLANK   LOGICAL  Flag for applying Planck function
!   ISRFL   INTEGER  Integer for setting reflectance flags
!   LOOKUP  LOGICAL  Flag for up-looking instruments
!
!   INPUTS/OUTPUTS:
!
!   LTYP    INTEGER  Identifier for clear/cloudy layers
!   SBEAM   LOGICAL  Flag for solar calculations
!   NAZ     INTEGER  Number of azimuthal directions in scattering
!                    calculations
!   NBOT0   INTEGER  Lowest atmospheric layer, at which RT calculation starts
!   NBOT    INTEGER  Same as NBOT0 for non-zero azimuthal angle
!   NTOP0   INTEGER  Highest atmospheric layer, at which RT calculation
!                    starts
!   NTOP    INTEGER  Same as NTOP0 for non-zero azimuthal angle
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     This subroutine adjusts LTYP for the clear layers and adjusts the
!     atmospheric boundaries to avoid in many cases treating systematically
!     all the input layers in the radiative transfer calculations (for
!     instance, in the solar case one can ignore the uppermost layers if
!     they are non-reflecting since they do not produce a diffuse radiation
!     field. Similarly, in presence of a Lambertian surface, there is no
!     contribution from the surface and the clear layers adjacent to the
!     surface for high order azimuthal components. CHARTS calculations
!     will in this case, be restricted to the portion of the atmosphere
!     between the top and bottom scattering layers). This subroutine also
!     resets ISBEAM, IPLANK and NAZ in order to save on the amount of
!     computation in some special cases (e.g. when the observer faces a
!     non-reflecting atmosphere).
!
!     Values of LTYP for clear layers are modified to indicate the
!     relative position of the each layer within a set of contiguous
!     clear layers above and below the observer level. LTYP values
!     are assigned values as follows:
!
!     LTYP = -1   - first clear layer
!          =  0   - middle layers
!          = -99  - last clear layer
!          = -100 - single clear layer
!
!     The input values of LTYP for cloudy layers indicate the actual cloud
!     model used for each layer. LTYP is modified in this routine to provide
!     a unique index value (>1) for each cloud layer.
!
!     Output variables:
!
!     NBOT0= lowest atmospheric layer below observer level at  which
!            one should start the radiative transfer calculations for
!            MAZ=0
!     NBOT = same as NBOT0 for MAZ >0
!     NTOP0= highest atmospheric layer above observer level at  which
!            one should start the radiative transfer calculations for
!            MAZ=0
!     NTOP = same as NTOP0 for MAZ >0
!
!     REVISED:  21 MARCH 1995
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    INTEGER, INTENT(IN) :: nlay, isrfl, lobs
    INTEGER, DIMENSION(nlay), INTENT(INOUT) :: ltyp
    LOGICAL, INTENT(IN) :: plank, lookup
    LOGICAL, INTENT(INOUT) :: sbeam
    INTEGER, INTENT(INOUT) :: nbot0, nbot, ntop0, ntop
    INTEGER, INTENT(INOUT) :: naz
    LOGICAL :: refls0, refls, reflt, reflb0, reflb
    INTEGER :: i
    DATA refls0, refls /.false., .false./
!
! *** Surface reflectance
!
    IF (isrfl .NE. 0) refls0 = .true.
    IF (isrfl .GT. 0) refls = .true.
!
! *** Set LTYP for clear layers above observer level (MERGE DOWN case)
!
    ntop = 1
    reflt = .false.
    DO i=1,lobs-1
      IF (ltyp(i) .EQ. 0) THEN
!
        IF (i .EQ. 1) THEN
          ltyp(i) = -1
        ELSE IF (ltyp(i-1) .GT. 0) THEN
          ltyp(i) = -1
        END IF
!
        IF (i .EQ. lobs - 1) THEN
          ltyp(i) = -99 + ltyp(i)
        ELSE IF (ltyp(i+1) .GT. 0) THEN
          ltyp(i) = -99 + ltyp(i)
        END IF
!
        IF (.NOT.reflt) ntop = ntop + 1
      ELSE
        reflt = .true.
      END IF
    END DO
!
! *** Set LTYP for clear layers below observer level (MERGE UP case)
!
    nbot = nlay
    reflb = .false.
    DO i=nlay,lobs,-1
      IF (ltyp(i) .EQ. 0) THEN
!
        IF (i .EQ. nlay) THEN
          ltyp(i) = -1
        ELSE IF (ltyp(i+1) .GT. 0) THEN
          ltyp(i) = -1
        END IF
!
        IF (i .EQ. lobs) THEN
          ltyp(i) = -99 + ltyp(i)
        ELSE IF (ltyp(i-1) .GT. 0) THEN
          ltyp(i) = -99 + ltyp(i)
        END IF
!
        IF (.NOT.reflb) nbot = nbot - 1
      ELSE
        reflb = .true.
      END IF
    END DO
!
    reflb0 = reflb .OR. refls0
    reflb = reflb .OR. refls
!
! *** Adjust boundaries for general solar case
!
    nbot0 = nbot
    ntop0 = ntop
    IF (refls0) nbot0 = nlay
    IF (refls) nbot = nlay
!
! *** Adjust boundaries for thermal case
!
    IF (plank) THEN
      nbot0 = nlay
      ntop0 = 1
    END IF
!
! *** Special case: No forward reflection - shrink atmospheric
! *** boundaries and force SBEAM = .FALSE. (assume that observer
! *** never looks directly into the sun.)
!
    IF (lookup .AND. (.NOT.reflt)) THEN
      sbeam = .false.
      naz = 0
!         ISRFL=-99
      nbot0 = lobs - 1
    END IF
!
    IF (.NOT.lookup) THEN
      IF (.NOT.reflb) naz = 0
      IF (.NOT.reflb0) THEN
        sbeam = .false.
        ntop0 = lobs
      END IF
    END IF
!
! *** Eliminate trivial case
!
!      IF ((LOOKUP.AND.LOBS.EQ.1).OR.&
!         (.NOT.LOOKUP.AND..NOT.REFLB0.AND.LOBS.EQ.NLAY+1))&
!         STOP 'CASE NOT TREATED'
    RETURN
  END SUBROUTINE SETCTL
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     OSS (Optimal Spectral Sampling)
!     $Name:  $
!     $Id: oss_addbl.f90,v 1.6 2009/04/03 20:21:19 crichard Exp $
!     Copyright AER, Inc., 2002, 2003. All rights Reserved.
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE SURFCE(maz, isrfl, tsfc, salb, f0, r, em, refls, vn, emissivity)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   SURFCE
!
! PURPOSE:
!
!   Subroutine initializes the reflection matrix and emission vector for the
!   lower atmosphere.
!
! SYNTAX:
!
!   CALL SURFCE(MAZ, ISRFL, TSFC, SALB, F0, R, EM, REFLS, VN)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   MAZ    INTEGER  Index for azimuthal direction in scattering calculations
!   ISRFL  INTEGER  Integer for setting reflectance flags
!   TSFC   REAL(FP)     Surface-level air temperature
!   SALB   REAL(FP)     Albedo
!   F0     REAL(FP)     solar irradiance at top of scattering layer
!   VN     REAL(Double)   Wavenumber, units of cm^-1
!
!   INPUTS/OUTPUTS:
!
!   R      REAL(FP)     Reflection matrix at NBOT0/NBOT
!   EM     REAL(FP)     Emission vector at NBOT0/NBOT
!   REFLS  LOGICAL  Flag for calculating surface reflectance
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     This subroutine initializes the reflection matrix (R) and emission
!     vector (EM) for the lower atmosphere by applying boundary conditions
!     at the surface. The current version treats Lambertian surfaces only.
!
!     MODIFIED:  21 MARCH 1995 to implement spectral interpolation
!                of surface properties within panel
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
    real(fp)                 :: emissivity
    INTEGER, INTENT(IN) :: maz, isrfl
    REAL(FP), INTENT(IN) :: tsfc, salb, f0
    INTEGER :: nang
    REAL(FP), DIMENSION(nang, nang), INTENT(INOUT) :: r
    REAL(FP), DIMENSION(nang), INTENT(INOUT) :: em
    REAL(Double), INTENT(IN) :: vn
    LOGICAL, INTENT(INOUT) :: refls
    LOGICAL :: plnck, sbeam
    REAL(FP) :: ylm0
    REAL(FP) :: pth
    REAL(FP) :: ylm
    REAL(FP) :: gmu
    REAL(FP) :: gwt
    REAL(FP) :: gmu0
    COMMON /angl/ pth(mxang), gmu(mxang), gwt(mxang), gmu0, ylm(0:2*&
&    mxcang-1, mxang), ylm0(0:2*mxcang-1)
    save /angl/
    REAL(FP) :: cl2inv
    REAL(FP) :: pi
    COMMON /constn/ pi, cl2inv
    save /constn/
    COMMON /lpar/ plnck, sbeam
    save /lpar/
    INTEGER :: nang2
    INTEGER :: n2ang
    COMMON /pnlpar/ nang, n2ang, nang2
    save /pnlpar/
    REAL(FP) :: rs
    INTEGER :: j
    REAL(FP) :: fac
    INTEGER :: i
    REAL(FP) :: b
    REAL(FP) :: wk
!
!
    refls = .false.
    em(1:nang)=0.
!
    IF (isrfl .EQ. -99) THEN
      RETURN
    ELSE IF (maz .GT. 0) THEN
!      IF(ISRFL.GT.0)GOTO 100
!
! *** LAMBERTIAN SURFACE
!
      RETURN
    ELSE
!
      rs = 0.
!
      IF (isrfl .EQ. -1) THEN
!
! *** Calculate surface reflectance across panel
!
        rs = salb
        refls = .true.
      ELSE IF (isrfl .EQ. 2) THEN
         rs = salb
         refls = .true.  !specular surface
        !STOP
      END IF
!
! *** Calculate reflection matrix
!
      IF (refls) THEN
         IF (isrfl .EQ. -1) THEN
            DO j=1,nang-1
               fac = 2._FP*gmu(j)*gwt(j)
               DO i=1,nang
                  r(j, i) = rs*fac
               END DO
            END DO
            DO i=1,nang
               r(nang, i) = 0.
            END DO
         ELSEIF (isrfl .EQ. 2) THEN
            r = 0
            DO j=1,nang
               r(j, j) = rs
            ENDDO
         ENDIF
      END IF

!
! *** Initialize EM by combining thermal emission and reflection of direct
! *** solar beam at the surface
!
      IF (plnck) THEN
        b = WNPLAN(vn, tsfc)
        wk = emissivity*b
!
        DO i=1,nang
          em(i) = wk
        END DO
      END IF
!
      IF (sbeam .AND. refls) THEN
        fac = gmu0/pi
        wk = rs*f0*fac
!
        DO i=1,nang
          em(i) = em(i) + wk
        END DO
      END IF
!
! *** General case (not implemented)
!
      RETURN
    END IF
  END SUBROUTINE SURFCE
!
  SUBROUTINE CLRLAY(tau, ta, tup, reflg, ltyp, vn, tlay, empl, emml, cosmicBackground, linInTauFlag)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   CLRLAY
!
! PURPOSE:
!
!   Calculate transmission and emission for a clear slab.
!
! SYNTAX:
!
!   CALL CLRLAY(TAU, TA, TUP, REFLG, LTYP, VN)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   TAU    REAL(FP)     Optical depth
!   TA     REAL(FP)     Lower level temperature
!   TUP    REAL(FP)     Upper level temperature
!   REFLG  LOGICAL  Flag for calculating surface reflectance
!   LTYP   INTEGER  Identifier for clear/cloudy layers
!   VN     REAL(Double)   Wavenumber, units of cm^-1
!   REAL(FP), INTENT(IN)    :: cosmicBackground - diffuse source at the top of atmosphere
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     This subroutine does the transmission and emission calculations
!     for a clear slab. The routine processes NLIM spectral points.
!     Transmission and forward and backward emission for each successive
!     clear layer are accumulated respectively in TLAY, EMPL and EMML.
!     If no thermal emission the routine accumulates the optical depths
!     in TLAY and calculates the transmission only when the last clear layer
!     of a series is reached. Linear-in-tau approximation is used for the
!     thermal source and it is assumed that the Planck function varies
!     linearly across the spectral interval spanned by the NLIM spectral
!     points.
!
!     REVISED:  11 OCTOBER 1993
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
    REAL(Double), INTENT(IN) :: vn
    INTEGER, INTENT(IN) :: ltyp
    LOGICAL, INTENT(IN) :: reflg
    REAL(FP), INTENT(IN) :: ta, tup, tau
    REAL(FP), INTENT(INOUT) :: tlay(:), empl(:), emml(:)
    REAL(FP),  INTENT(IN)    :: cosmicBackground
    LOGICAL                  :: linInTauFlag

!
    LOGICAL :: plnck, sbeam
    REAL(FP) :: ylm0
    REAL(FP) :: pth
    REAL(FP) :: ylm
    REAL(FP) :: gmu
    REAL(FP) :: gwt
    REAL(FP) :: gmu0
!
    COMMON /angl/ pth(mxang), gmu(mxang), gwt(mxang), gmu0, ylm(0:2*&
&    mxcang-1, mxang), ylm0(0:2*mxcang-1)
    save /angl/
    REAL(FP) :: cl2inv
    REAL(FP) :: pi
    COMMON /constn/ pi, cl2inv
    save /constn/
    COMMON /lpar/ plnck, sbeam
    save /lpar/
    INTEGER :: nang
    INTEGER :: nang2
    INTEGER :: n2ang
    COMMON /pnlpar/ nang, n2ang, nang2
    save /pnlpar/
    REAL(FP) :: tx
    REAL(FP) :: emm
    REAL(FP) :: emp
    REAL(FP) :: wk
    dimension tx(mxang), emp(mxang), emm(mxang), wk(mxang)
    INTEGER :: i
    REAL(FP) :: ba
    REAL(FP) :: bb
    REAL(FP) :: bav
    REAL(FP) :: delb
    INTEGER :: n
    REAL(FP) :: em0
    REAL(FP) :: cci
    REAL(FP) :: emd
    INTRINSIC EXP
!
!     Scale optical depths for all incident angles including observer
!     angle.
!
    DO i=1,nang
      wk(i) = tau*pth(i)
    END DO
!
! *** THERMAL CASE.
!
    IF (plnck) THEN
!
!     Calculate vertically averaged Planck function for the layer and
!     correction term for linear-in-Tau approximation at each spectral
!     point
      if (linInTauFlag) then
         ba = WNPLAN(vn, ta)
         bb = WNPLAN(vn, tup)
         bav = .5*(bb+ba)
         delb = .5*(bb-ba)
      else
! Or omit linear in tau
         bav = wnplan(vn,(ta+tup)/2._FP)
         delb = 0.
      end if
!
!     Compute transmittance and thermal emission
!
      DO n=1,nang
        tx(n) = EXP(-wk(n))
        em0 = 1. - tx(n)
        IF (wk(n) .GT. 1.e-03) THEN
          cci = 2._FP/wk(n)
          emd = cci*(tx(n)-1.) + (tx(n)+1.)
        ELSE
          emd = wk(n)*wk(n)*(1.-0.5*wk(n))/6.
        END IF
        emp(n) = em0*bav + emd*delb
        emm(n) = em0*bav - emd*delb
      END DO
!
!     Accumulate partial results for clear slab.
!
      IF (ltyp .EQ. -1 .OR. ltyp .EQ. -100) THEN
! *** First layer
        empl(1:nang)=emp(1:nang)
        IF (reflg) emml(1:nang)=emm(1:nang)
        tlay(1:nang)=tx(1:nang)
      ELSE
! *** Update forward and backward emission for subsequent layers
        wk(1:nang)=tx(1:nang)*empl(1:nang)
        empl(1:nang) = wk(1:nang) + emp(1:nang)
        IF (reflg) THEN
          wk(1:nang)=tlay(1:nang)*emm(1:nang)
          emml(1:nang) = wk(1:nang) + emml(1:nang)
        END IF
! *** Update transmission
        tlay(1:nang)=tlay(1:nang)*tx(1:nang)
        IF (reflg) THEN
          emml(1:nang)=emml(1:nang) + cosmicBackground*tlay(1:nang)
        END IF
      END IF
    ELSE
!
! *** SOLAR ONLY
!
      IF (ltyp .EQ. -1 .OR. ltyp .EQ. -100) THEN
        tlay(1:nang)=wk(1:nang)
      ELSE
        tlay(1:nang) = wk(1:nang) + tlay(1:nang)
      END IF
!
! *** Compute transmittance for combined consecutive clear layers
!
      IF (ltyp .LE. -99) THEN
        DO n=1,nang
          tlay(n) = EXP(-tlay(n))
        END DO
      END IF
    END IF
    RETURN
  END SUBROUTINE CLRLAY
!
  SUBROUTINE MSLAY(text, ssalb, ta, tup, f0, rlay, tlay, empl, emml, maz&
&    , vn, alpha, mrgup,linInTauFlag)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   MSLAY
!
! PURPOSE:
!
!   Apply doubling tables for a cloudy layer. Compute forward and backward
!   internal emission due to thermal and solar sources.
!
! SYNTAX:
!
!   CALL MSLAY(TEXT, SSALB, TA, TUP, F0, RLAY, TLAY, EMPL, EMML, MAZ,
!      VN, ALPHA, MRGUP)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   TEXT   REAL(FP)     Total extinction
!   SSALB  REAL(FP)     Albedo
!   TA     REAL(FP)     Lower level temperature
!   TUP    REAL(FP)     Upper level temperature
!   F0     REAL(FP)     solar irradiance at top of scattering layer
!   MAZ    INTEGER  Index for azimuthal direction in scattering calculations
!   VN     REAL(Double)   Wavenumber, units of cm^-1
!   ALPHA  REAL(FP)     Phase function
!   MRGUP  INTEGER  Flag for down/upward calculations
!
!   INPUTS/OUTPUTS:
!
!   RLAY   REAL(FP)     Layer refrection
!   TLAY   REAL(FP)     Layer transmission
!   EMPL   REAL(FP)     Forward emission
!   EMML   REAL(FP)     Backward emission
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     This subroutine does the interpolation of the doubling results
!     for a cloudy layer. The routine processes NLIM spectral points.
!     Reflection, transmission and forward and backward emission are output
!     respectively in RLAY, TLAY, EMPL and EMML. Linear-in-tau approximation
!     is used for the thermal source and it is assumed that the Planck
!     function varies linearly across the spectral interval spanned by the
!     NLIM spectral points.
!
!     LAST REVISED: MARCH 1995 to implement linear spectral interpolation
!                   of the doubling results.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
    REAL(FP), DIMENSION(0:2*mxcang), INTENT(IN) :: alpha
    INTEGER, INTENT(IN) :: mrgup, maz
    INTEGER :: nang2
    REAL(FP), DIMENSION(nang2), INTENT(INOUT) :: rlay, tlay
    INTEGER :: nang
    REAL(FP), DIMENSION(nang), INTENT(INOUT) :: empl, emml
    REAL(FP), INTENT(IN) :: text, ssalb
    REAL(FP), INTENT(IN) :: ta, tup, f0
    REAL(Double), INTENT(IN) :: vn
    LOGICAL                  :: linInTauFlag
!
    REAL(FP) :: em(2*mxang), sbm(2*mxang)
    LOGICAL :: plnck, sbeam
!
    COMMON /lpar/ plnck, sbeam
    save /lpar/
    INTEGER :: n2ang
    COMMON /pnlpar/ nang, n2ang, nang2
    save /pnlpar/
    REAL(FP) :: ba
    REAL(FP) :: bb
    REAL(FP) :: bav
    REAL(FP) :: delb
    INTEGER :: n
!
    empl(1:nang)=0.
    emml(1:nang)=0.
!
    IF (plnck) THEN
      if (linInTauFlag) then
      ba = WNPLAN(vn, ta)
      bb = WNPLAN(vn, tup)
      bav = .5*(bb+ba)
      delb = .5*(bb-ba)
    else
!       Or omit linear in tau
        bav = wnplan(vn,(ta+tup)/2_FP)
        delb = 0.
      end if
    END IF
!
!     Get doubling tables for current cloud layer
!
!CALL GENTAB(TABS,TSCAT,ALPHA,MAZ,RLAY,TLAY,EM,SBM)
    CALL GENTAB(text, ssalb, alpha, maz, rlay, tlay, em, sbm)
!
!     Compute forward and backward internal emission due to
!     thermal and solar sources
!
    IF (plnck) THEN
!
!     Calculate vertically averaged Planck function for the layer and
!     correction term for linear-in-Tau approximation at each spectral
!     point
!
      DO n=1,nang
        empl(n) = em(n)*bav + em(nang+n)*delb
        emml(n) = em(n)*bav - em(nang+n)*delb
      END DO
    END IF
!
    IF (sbeam) THEN
!
!     Add solar component
!
      IF (mrgup .EQ. 1) THEN
        DO n=1,nang
          empl(n) = empl(n) + sbm(n)*f0
          emml(n) = emml(n) + sbm(nang+n)*f0
        END DO
      ELSE
        DO n=1,nang
          empl(n) = empl(n) + sbm(nang+n)*f0
          emml(n) = emml(n) + sbm(n)*f0
        END DO
      END IF
    END IF
!
    RETURN
  END SUBROUTINE MSLAY
!
  SUBROUTINE GENTAB(text, ssalb, alpha, maz, r, trnsm, em, sbm)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   GENTAB
!
! PURPOSE:
!
!   Generate look-up tables of reflection, transmission matrices and source
!   vectors by doubling calculations.
!
! SYNTAX:
!
!   CALL GENTAB(TEXT, SSALB, ALPHA, MAZ, R, TRNSM, EM, SBM)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   TEXT   REAL(FP)     Total extinction
!   SSALB  REAL(FP)     Albedo
!   ALPHA  REAL(FP)     Phase function
!   MAZ    INTEGER  Index for azimuthal direction in scattering calculations
!
!   INPUTS/OUTPUTS:
!
!   R      REAL(FP)     Reflection matrix at NBOT0/NBOT
!   TRNSM  REAL(FP)     Transmission matrix at NBOT0/NBOT
!   EM     REAL(FP)     Emission vector at NBOT0/NBOT
!   SBM    REAL(FP)     upward and downward emission due to solar scattering
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Generate look up tables of reflexion, transmission matrices
!     and source vectors by doubling calculations for the range of
!     absorption optical depths specified by NMIN and NMAX. NMIN0
!     and NMAX0 indicate the range of the calculations covered in
!     previous calls to GENTAB so that calculations corresponding
!     to NMIN0<N<NMAX0 do not have to be redone. Table entries are
!     in the common block /TABX/.
!
!     NOTE:
!     Internal sources due to thermal emission and solar scattering
!     are stored in vectors EM and SBM respectively. The first NANG
!     elements of SBM correspond to the reflected solar radiation,
!     while the second NANG elements correspond to the transmitted
!     solar radiation. For the thermal sources, the first NANG elements
!     contain the emission for an isothermal layer which is the same
!     in the downward and upward directions. The second set of elements
!     are correction terms which account for the assumed linear-in-tau
!     dependence of the PLanck function within the layer. Upward and
!     downward thermal sources are derived in the subroutine MSLAY by
!     respectively adding and subtracting the correction terms from the
!     first set of elements.
!
!     REVISED:  21 MARCH 1995
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
    REAL(FP), DIMENSION(0:2*mxcang), INTENT(IN) :: alpha
    REAL(FP), INTENT(IN) :: text, ssalb
    INTEGER, INTENT(IN) :: maz
    INTEGER :: nang2
    REAL(FP), DIMENSION(nang2), INTENT(INOUT) :: r, trnsm
    INTEGER :: n2ang
    REAL(FP), DIMENSION(n2ang), INTENT(INOUT) :: em, sbm
!
    REAL(FP) :: pfp(mxang2), pfm(mxang2), pfp0(mxang), pfm0(mxang)
    LOGICAL :: plnck, sbeam
    REAL(FP) :: ylm0
    REAL(FP) :: pth
    REAL(FP) :: ylm
    REAL(FP) :: gmu
    REAL(FP) :: gwt
    REAL(FP) :: gmu0
!
    COMMON /angl/ pth(mxang), gmu(mxang), gwt(mxang), gmu0, ylm(0:2*&
&    mxcang-1, mxang), ylm0(0:2*mxcang-1)
    save /angl/
    COMMON /lpar/ plnck, sbeam
    save /lpar/
    INTEGER :: nang
    COMMON /pnlpar/ nang, n2ang, nang2
    save /pnlpar/
    REAL(FP) :: cl2inv
    REAL(FP) :: pi
    COMMON /constn/ pi, cl2inv
    save /constn/
    INTEGER :: i
    INTEGER :: j
    INTEGER :: ij
    REAL(FP) :: sgn
    INTEGER :: l
    REAL(FP) :: pij
    REAL(FP) :: scal
    REAL(FP) :: pi0
    REAL(FP) :: exp2
    INTEGER :: nit
    REAL(FP) :: dtau
    REAL(FP) :: att
    REAL(FP) :: pwr1
    INTRINSIC EXP
    INTRINSIC MAX
    INTRINSIC LOG
    INTRINSIC INT
    INTEGER :: y1
!
! *** Construct phase matrix.
! *** Pre-scale matrices for diamond initialization (DI) and apply
! *** angular weighting.
!
    DO i=1,nang
      DO j=1,nang-1
        ij = (i-1)*nang + j
        pfp(ij) = 0.
        pfm(ij) = 0.
        sgn = 1.
        DO l=maz,n2ang-3
          pij = alpha(l)*ylm(l, i)*ylm(l, j)
          pfp(ij) = pfp(ij) + pij
          pfm(ij) = pfm(ij) + sgn*pij
          sgn = -sgn
        END DO
        pfp(ij) = -(0.25*pfp(ij)*pth(i)*gwt(j))
        pfm(ij) = 0.25*pfm(ij)*pth(i)*gwt(j)
      END DO
! *** Observer angle
      pfp(ij+1) = 0.
      pfm(ij+1) = 0.
    END DO
!
! *** Construct phase matrix for solar source.
!
    IF (sbeam) THEN
      IF (maz .EQ. 0) THEN
        scal = .25*gmu0
      ELSE
        scal = .5*gmu0
      END IF
!
      DO i=1,nang
        pfp0(i) = 0.
        pfm0(i) = 0.
        sgn = 1.
        DO l=maz,n2ang-3
          pi0 = alpha(l)*ylm(l, i)*ylm0(l)
          pfp0(i) = pfp0(i) + pi0
          pfm0(i) = pfm0(i) + sgn*pi0
          sgn = -sgn
        END DO
        pfp0(i) = scal*pfp0(i)*pth(i)
        pfm0(i) = scal*pfm0(i)*pth(i)
      END DO
    END IF
!
! *** Use DI approximation when layer extinction optical depth
! *** is less than 2**-10.
!
!      TEXT=TABS+TSCAT
!      SSALB=TSCAT/TEXT
! *** Compute initial layer optical depth and number of doublings
! *** necessary to reach TEXT
    exp2 = LOG(text)*cl2inv
    IF (exp2 .LT. -30.) THEN
      exp2 = -30.
    ELSE
      exp2 = exp2
    END IF
    y1 = INT(exp2 + 4.)
    IF (0 .LT. y1) THEN
      nit = y1
    ELSE
      nit = 0
    END IF
    pwr1 = 2.**nit
    dtau = text/pwr1
    IF (sbeam) att = EXP(-(dtau/gmu0))
! *** Initialize doubling calculations
    CALL DI(dtau, ssalb, pfp, pfm, pfp0, pfm0, r, trnsm, em, sbm)
    CALL DOUBL(nit, att, r, trnsm, em, sbm)
!
    RETURN
  END SUBROUTINE GENTAB
!
  SUBROUTINE DI(deltau, ssalb, pfp, pfm, pfp0, pfm0, r, trnsm, em, sbm)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   DI
!
! PURPOSE:
!
!   Calculate reflection, transmission and sources for an infinitesimal layer by
!   using the Diamond Initialization scheme.
!
! SYNTAX:
!
!   CALL DI(DELTAU, SSALB, PFP, PFM, PFP0, PFM0, R, TRNSM, EM, SBM)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   DELTAU  REAL(FP)  TBD
!   SSALB   REAL(FP)  Albedo
!
!   INPUTS/OUTPUTS:
!
!   PFP     REAL(FP)  Phase matrix
!   PFM     REAL(FP)  Phase matrix
!   PFP0    REAL(FP)  Phase matrix for solar source
!   PFM0    REAL(FP)  Phase matrix for solar source
!   R       REAL(FP)  Reflection matrix at NBOT0/NBOT
!   TRNSM   REAL(FP)  Transmission matrix at NBOT0/NBOT
!   EM      REAL(FP)  Emission vector at NBOT0/NBOT
!   SBM     REAL(FP)  upward and downward emission due to solar scattering
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     This routine calculates reflection, transmission and sources for
!     an infinitesimal layer using the Diamond Initialization scheme.
!     (for the treatment of the thermal emission, see note above in
!     subroutine GENTAB)
!
!     REF: Wiscombe,1976, JQSRT, pp637-658.
!
!     REVISED:  04 OCTOBER 1993
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
    REAL(FP), INTENT(IN) :: deltau, ssalb
    INTEGER :: nang
    REAL(FP), DIMENSION(nang, nang), INTENT(INOUT) :: r, trnsm
    REAL(FP), DIMENSION(nang, 2), INTENT(INOUT) :: em, sbm
    REAL(FP), DIMENSION(nang, nang), INTENT(INOUT) :: pfp, pfm
    REAL(FP), DIMENSION(nang), INTENT(INOUT) :: pfp0, pfm0
!
    LOGICAL :: plnck, sbeam
    REAL(FP) :: ylm0
    REAL(FP) :: pth
    REAL(FP) :: ylm
    REAL(FP) :: gmu
    REAL(FP) :: gwt
    REAL(FP) :: gmu0
!
    COMMON /angl/ pth(mxang), gmu(mxang), gwt(mxang), gmu0, ylm(0:2*&
&    mxcang-1, mxang), ylm0(0:2*mxcang-1)
    save /angl/
    REAL(FP) :: cl2inv
    REAL(FP) :: pi
    COMMON /constn/ pi, cl2inv
    save /constn/
    COMMON /lpar/ plnck, sbeam
    save /lpar/
    INTEGER :: nang2
    INTEGER :: n2ang
    COMMON /pnlpar/ nang, n2ang, nang2
    save /pnlpar/
    REAL(FP) :: ti
    REAL(FP) :: gamma
    REAL(FP) :: rtinv
    REAL(FP) :: tinv
    REAL(FP) :: wkem
    REAL(FP) :: wk
    REAL(FP) :: wksbm
    dimension gamma(mxang2), ti(mxang2), tinv(mxang2), rtinv(mxang2&
&    ), wk(mxang2), wkem(mxang), wksbm(2*mxang)
    REAL(FP) :: c1
    REAL(FP) :: c2
    INTEGER :: i
    INTEGER :: ii
!!$    REAL(FP) :: det
    REAL(FP) :: cem
    INTEGER :: j
    INTEGER :: ij
    REAL(FP) :: csol
    INTRINSIC EXP
    REAL(FP) :: TWO=2.
    integer :: ierr
!
    c1 = deltau*ssalb
    c2 = 0.5*deltau
    ti(1:nang2) = reshape(pfp(:,:) * c1,(/nang2/))
    r(:,:) = pfm(:,:) * c1
    DO i=1,nang
      ii = (i-1)*nang + i
      ti(ii) = ti(ii) + pth(i)*c2 + 1.
    END DO
!
!     Build Gamma matrix - scale by a factor 2
!
    tinv(1:nang2)=ti(1:nang2)
!$    CALL sMINV(TINV,NANG,NANG,WK,DET,TOL,0,1)
    CALL SMINV(tinv, nang, ierr)
    CALL MXM(tinv, nang, r, nang, rtinv, nang)
    CALL MXM(r, nang, rtinv, nang, wk, nang)
    gamma(1:nang2)=ti(1:nang2)-wk(1:nang2)
!$    CALL sMINV(GAMMA,NANG,NANG,WK,DET,TOL,0,1)
    CALL SMINV(gamma, nang, ierr)
    gamma(1:nang2) = gamma(1:nang2) * two
!
!     Compute reflection matrix
!
    CALL MXM(rtinv, nang, gamma, nang, r, nang)
!
!     Compute transmission matrix
!
    trnsm(:,:)=reshape(gamma(1:nang2),(/nang,nang/))
    DO i=1,nang
      trnsm(i, i) = trnsm(i, i) - 1.
    END DO
!
!     Compute thermal emission terms
!
    IF (plnck) THEN
      cem = (deltau-c1)*0.5
      wkem(1:nang) = pth(1:nang) * cem
      wk(1:nang2) = gamma(1:nang2) + reshape(r(:,:),(/nang2/))
!        CALL MXV(WK,NANG,WKEM,NANG,EM)
      DO i=1,nang
        em(i, 1) = 0.
        DO j=1,nang
          ij = (i-1)*nang + j
          em(i, 1) = em(i, 1) + wk(ij)*wkem(j)
        END DO
      END DO
! *** No correction term for infinitesimal layer
      em(1:nang,2)=0.
    END IF
!
!     Compute upward and downward emission due to solar scattering
!
    IF (sbeam) THEN
      csol = .5*ssalb/pi*(1.-EXP(-(deltau/gmu0)))
      wksbm(1:nang) = pfm0(1:nang) * csol
      wksbm(nang+1:2*nang) = pfp0(1:nang) * csol
      DO i=1,nang
        sbm(i, 1) = 0.
        sbm(i, 2) = 0.
        DO j=1,nang
          ij = (i-1)*nang + j
          sbm(i, 1) = sbm(i, 1) + gamma(ij)*wksbm(j) + r(j, i)*wksbm(&
&            nang+j)
          sbm(i, 2) = sbm(i, 2) + gamma(ij)*wksbm(nang+j) + r(j, i)*&
&            wksbm(j)
        END DO
      END DO
    END IF
    RETURN
  END SUBROUTINE DI
!
  SUBROUTINE DOUBL(nit, att, r, trnsm, em, sbm)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   DOUBL
!
! PURPOSE:
!
!   Perform doubling calculations for a homogeneous scattering layer.
!
! SYNTAX:
!
!   CALL DOUBL(NIT, ATT, R, TRNSM, EM, SBM)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   NIT    INTEGER  number of doubling calculations
!
!   INPUTS/OUTPUTS:
!
!   ATT    REAL(FP)     TBD
!   R      REAL(FP)     Reflection matrix at NBOT0/NBOT
!   TRNSM  REAL(FP)     Transmission matrix at NBOT0/NBOT
!   EM     REAL(FP)     Emission vector at NBOT0/NBOT
!   SBM    REAL(FP)     upward and downward emission due to solar scattering
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     This subroutine does the doubling calculations for a homogeneous
!     scattering layer. The number of calculations is specified by NIT.
!     It is assumed that the reflection and transmission matrices and the
!     source vectors have already been initialized.
!
!     REVISED:  21 MARCH 1995
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
    INTEGER, INTENT(IN) :: nit
    REAL(FP), INTENT(INOUT) :: att
    INTEGER :: nang2
    REAL(FP), DIMENSION(nang2), INTENT(INOUT) :: r, trnsm
    INTEGER :: n2ang
    REAL(FP), DIMENSION(n2ang), INTENT(INOUT) :: em, sbm
!
    LOGICAL :: plnck, sbeam
!
    COMMON /lpar/ plnck, sbeam
    save /lpar/
    INTEGER :: nang
    COMMON /pnlpar/ nang, n2ang, nang2
    save /pnlpar/
    REAL(FP) :: wk1
    REAL(FP) :: wk2
    REAL(fp) :: wk3(mxang2)
    REAL(FP) :: rm
    REAL(FP) :: wkt
    dimension rm(mxang2), wk1(mxang2), wk2(mxang2), wkt(mxang2)
    INTEGER :: n
    INTEGER :: i
!!$    REAL(FP) :: det
    integer :: ierr
!
    DO n=1,nit
!
!     Compute multiple reflection term
!
      wk3(1:nang2)=r(1:nang2)
      CALL MXM(r, nang, wk3, nang, rm, nang)
      DO i=1,nang2
        rm(i) = -rm(i)
      END DO
      do i=1,nang2,nang+1
         rm(i)=rm(i)+1.
      end do

!$       CALL sMINV(RM,NANG,NANG,WK1,DET,TOL,0,1)
      CALL SMINV(rm, nang, ierr)
!
!     Compute new transmission matrix
!
      CALL MXM(rm, nang, trnsm, nang, wk1, nang)
      wkt(1:nang2)=trnsm(1:nang2)
      CALL MXM(wkt, nang, wk1, nang, trnsm, nang)
!
!     Double thermal and solar sources
!
      CALL MXM(r, nang, wk1, nang, wk2, nang)
      IF (plnck) CALL EMIS(em, wk1, wk2)
      IF (sbeam) CALL SUNB(sbm, wk1, wk2, att)
!
!     Compute new reflection matrix
!
      CALL MXM(wkt, nang, wk2, nang, wk1, nang)
      r(1:nang2) = r(1:nang2) + wk1(1:nang2)
    END DO
    RETURN
  END SUBROUTINE DOUBL
!
  SUBROUTINE EMIS(em, tm, tmr)
    IMPLICIT NONE
    INTEGER :: nang
!<f90Subroutine>********************************************************
!
! NAME:
!
!   EMIS
!
! PURPOSE:
!
!   Perform doubling of the thermal source
!
! SYNTAX:
!
!   CALL EMIS(EM, TM, TMR, WKEM)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   TM    REAL(FP)  TBD
!   TMR   REAL(FP)  TBD
!
!   INPUTS/OUTPUTS:
!
!   EM    REAL(FP)  Emission vector at NBOT0/NBOT
!   WKEM  REAL(FP)  TBD
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     This subroutine does the doubling of the thermal source
!     (see note in subroutine GENTAB)
!
!     REVISED:  04 OCTOBER 1993
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    REAL(FP), DIMENSION(nang, 2), INTENT(INOUT) :: em
    REAL(FP), DIMENSION(nang, nang), INTENT(IN) :: tm, tmr
    INTEGER :: nang2
    INTEGER :: n2ang
    REAL(FP), DIMENSION(nang, 2) :: wkem
!
    COMMON /pnlpar/ nang, n2ang, nang2
    save /pnlpar/
    INTEGER :: i
    INTEGER :: j
    REAL(FP) :: HALF=0.5
!
    wkem(:,:)=em(:,:)
    DO i=1,nang
      em(i, 2) = em(i, 2) + em(i, 1)
      DO j=1,nang
        em(i, 1) = em(i, 1) + (tm(j, i)+tmr(j, i))*wkem(j, 1)
        em(i, 2) = em(i, 2) + (tm(j, i)-tmr(j, i))*(wkem(j, 2)-wkem(j, 1&
&          ))
      END DO
    END DO
    em(1:nang, 2) = em(1:nang, 2)*half
    RETURN
  END SUBROUTINE EMIS
!
  SUBROUTINE SUNB(sbm, tm, tmr, att)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   SUNB
!
! PURPOSE:
!
!   Perform doubling of the solar source
!
! SYNTAX:
!
!   CALL SUNB(SBM, TM, TMR, WKEM, ATT)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   TM    REAL(FP)  TBD
!   TMR   REAL(FP)  TBD
!
!   INPUTS/OUTPUTS:
!
!   SBM   REAL(FP)  upward and downward emission due to solar scattering
!   WKEM  REAL(FP)  TBD
!   ATT   REAL(FP)  TBD
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     This subroutine does the doubling of the solar source
!     (see note in subroutine GENTAB)
!
!     REVISED:  04 OCTOBER 1993
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    REAL(FP), INTENT(INOUT) :: att
    INTEGER :: nang
    REAL(FP), DIMENSION(nang, 2), INTENT(INOUT) :: sbm
    REAL(FP), DIMENSION(nang, nang), INTENT(IN) :: tm, tmr
    INTEGER :: nang2
    INTEGER :: n2ang
    REAL(FP), DIMENSION(nang, 2) :: wkem
!
    COMMON /pnlpar/ nang, n2ang, nang2
    save /pnlpar/
    INTEGER :: i
    INTEGER :: j
!
    wkem(:,:)=sbm(:,:)
    sbm(1:nang, 2)=sbm(1:nang, 2)*att
    wkem(1:nang,1) = wkem(1:nang,1) * att
    DO i=1,nang
      DO j=1,nang
        sbm(i, 1) = sbm(i, 1) + tm(j, i)*wkem(j, 1) + tmr(j, i)*wkem(j, &
&          2)
        sbm(i, 2) = sbm(i, 2) + tm(j, i)*wkem(j, 2) + tmr(j, i)*wkem(j, &
&          1)
      END DO
    END DO
    att = att*att
    RETURN
  END SUBROUTINE SUNB
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     OSS (Optimal Spectral Sampling)
!     $Name:  $
!     $Id: oss_addbl.f90,v 1.6 2009/04/03 20:21:19 crichard Exp $
!     Copyright AER, Inc., 2002, 2003. All rights Reserved.
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ADDCLR(r, trnsm, em, reflg, plnck, tlay, empl, emml)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   ADDCLR
!
! PURPOSE:
!
!   Perform adding of a clear layer.
!
! SYNTAX:
!
!   CALL ADDCLR(R, TRNSM, EM, REFLG, PLNCK)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   REFLG  LOGICAL  Flag for calculating surface reflectance
!   PLNCK  LOGICAL  Flag for applying Planck function
!
!   INPUTS/OUTPUTS:
!
!   R      REAL(FP)     Reflection matrix at NBOT0/NBOT
!   TRNSM  REAL(FP)     Transmission matrix at NBOT0/NBOT
!   EM     REAL(FP)     Emission vector at NBOT0/NBOT
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     This subroutine does the adding of a clear layer - Transmission
!     and internal emission for the layer must be contained in common
!     block /LAYDAT/. Results are accumulated in arrays R, TRNSM and EM.
!     Transmittance is updated only if TXON = .TRUE., and reflection
!     matrix only if non-zero initially (REFLG=.FALSE.).
!
!     REVISED:  08 MARCH 1993
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
    REAL(FP), DIMENSION(:), INTENT(INOUT) :: r, trnsm, em
    LOGICAL, INTENT(IN) :: reflg, plnck
    REAL(FP), INTENT(IN) :: tlay(:), empl(:), emml(:)
    REAL(FP) :: cl2inv
    REAL(FP) :: pi
!
    COMMON /constn/ pi, cl2inv
    save /constn/
    INTEGER :: nang
    INTEGER :: nang2
    INTEGER :: n2ang
    COMMON /pnlpar/ nang, n2ang, nang2
    save /pnlpar/
    REAL(FP) :: tr
    REAL(FP) :: wkem
    dimension wkem(mxang), tr(mxang2)
!
!     Update emission
!
    wkem(1:nang)=tlay(1:nang)*em(1:nang)
!
   CALL DMXMV(tlay, trnsm, tr, nang)
   trnsm(1:nang2)=tr(1:nang2)

    IF (plnck) THEN
      em(1:nang) = wkem(1:nang) + empl(1:nang)
    ELSE
       em(1:nang)=wkem(1:nang)
    END IF
!
    IF (reflg) THEN
!
      CALL DMXMV(tlay, r, tr, nang)     ! tr(j,i) = tlay(i) * r(j,i)
!
      IF (plnck) THEN
        CALL MXVV(tr, emml, wkem, nang) ! wkem(i) = sum_j tr(i, j) * emml(j)
        em(1:nang) = wkem(1:nang) + em(1:nang)
      END IF
!
!     Update reflection matrix
!
      CALL MXDMV(tr, tlay, r, nang)     ! r(j, i) = tr(j, i)*tlay(j)
    END IF
    RETURN
  END SUBROUTINE ADDCLR
!
  SUBROUTINE ADDMS(r, trnsm, em, reflg, rlay, tlay, empl, emml)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   ADDMS
!
! PURPOSE:
!
!   General layer adding subroutine
!
! SYNTAX:
!
!   CALL ADDMS(R, TRNSM, EM, REFLG)
!
! ARGUMENTS:
!
!   INPUTS/OUTPUTS:
!
!   R      REAL(FP)     Reflection matrix at NBOT0/NBOT
!   TRNSM  REAL(FP)     Transmission matrix at NBOT0/NBOT
!   EM     REAL(FP)     Emission vector at NBOT0/NBOT
!   REFLG  LOGICAL  Flag for calculating surface reflectance
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     General layer adding routine - Reflection, transmission and
!     internal emission for the layer must be contained in common
!     block /LAYDAT/. Results are accumulated in arrays R, TRNSM and EM.
!     Transmittance is updated only if TXON = .TRUE., and reflection
!     matrix only if non-zero initially (REFLG=.FALSE.).
!
!     REVISED:  08 MARCH 1993
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
    REAL(FP), DIMENSION(:), INTENT(INOUT) :: r, trnsm, em
    REAL(FP), DIMENSION(:), INTENT(INOUT) :: rlay
    REAL(FP), DIMENSION(:), INTENT(IN) :: tlay, empl, emml
    LOGICAL, INTENT(INOUT) :: reflg
    INTEGER :: nang
    INTEGER :: nang2
    INTEGER :: n2ang
    COMMON /pnlpar/ nang, n2ang, nang2
    save /pnlpar/
    REAL(FP) :: wkem1
    REAL(FP) :: wk1
    REAL(FP) :: wkem2
    REAL(FP) :: wk2
    REAL(FP) :: rm
    dimension rm(mxang2), wk1(mxang2), wk2(mxang2), wkem1(mxang), &
&    wkem2(mxang)
!
!      EQUIVALENCE (WK0,RM)
!
    IF (reflg) THEN
!
!     Calculate multiple reflection term (result in RM)
!
      CALL MLREF(r, rlay, rm)
!
!     Update emission
!
      CALL MXMV(tlay, rm, wk1, nang)
      CALL MXMV(wk1, trnsm, wk2, nang)
      trnsm(1:nang2)=wk2(1:nang2)

      CALL MXMV(wk1, r, wk2, nang)
      CALL MXVV(wk1, em, wkem1, nang)
      CALL MXVV(wk2, emml, wkem2, nang)
      em(1:nang) = wkem1(1:nang) + wkem2(1:nang)
      em(1:nang) = em(1:nang) + empl(1:nang)
!
!     Update reflection matrix
!
!       -- old call:
!         CALL MXMV(WK2,TLAY,WK0,NANG)
!         CALL SVADD(WK0,RLAY,R,NANG2)
!       -- new call, removing equivalence:
      CALL MXMV(wk2, tlay, rm, nang)
      r(1:nang2) = rm(1:nang2) + rlay(1:nang2)
    ELSE
!
      CALL MXVV(tlay, em, wkem1, nang)
      em(1:nang) = empl(1:nang) + wkem1(1:nang)
!
      r(1:nang2)=rlay(1:nang2)
      trnsm(1:nang2)=tlay(1:nang2)
!
      reflg = .true.
!
    END IF
    RETURN
  END SUBROUTINE ADDMS
!
  SUBROUTINE MRGRAD(rf, emf, rb, emb, reflf, reflb, rad)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   MRGRAD
!
! PURPOSE:
!
!   Calculate the observed radiances by combining reflection and emission of the
!   atmosphere below and above the observer level.
!
! SYNTAX:
!
!   CALL MRGRAD(RF, EMF, RB, EMB, REFLF, REFLB, RAD)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   REFLF  LOGICAL  Flag to show whether in front of observer atmosphere
!                   reflecting
!   REFLB  LOGICAL  Flag to show whether behind observer atmosphere
!                   reflecting
!
!   INPUTS/OUTPUTS:
!
!   RF     REAL(FP)     Atmosphere Reflection in front of observer
!   EMF    REAL(FP)     Atmosphere Emission in front of observer
!   RB     REAL(FP)     Atmosphere Reflection behind observer
!   EMB    REAL(FP)     Atmosphere Emission behind observer
!   RAD    REAL(FP)     Radiance
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     This subroutine calculates the observed radiances by combining
!     reflection and emission of the atmosphere below and above the
!     observer level.
!
!     Flags REFLF/REFLB indicate whether the portion of the atmosphere
!     in front of/behind the observer is reflecting
!
!     REVISED:  11 OCTOBER 1993
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
    REAL(FP), DIMENSION(:), INTENT(INOUT) :: rf, emf, rb, emb, rad
    LOGICAL, INTENT(IN) :: reflf, reflb
    INTEGER :: nang
    INTEGER :: nang2
    INTEGER :: n2ang
!
    COMMON /pnlpar/ nang, n2ang, nang2
    save /pnlpar/
    REAL(FP) :: radtmp
    REAL(FP) :: rm
    dimension rm(mxang2), radtmp(mxang)
!
    IF (.NOT.reflf) THEN
       rad(1:nang)=emf(1:nang)
    ELSE IF (.NOT.reflb) THEN
      CALL MXVV(rf, emb, rad, nang)
      rad(1:nang) = emf(1:nang) + rad(1:nang)
    ELSE
      CALL MLREF(rf, rb, rm)
      CALL MXVV(rf, emb, radtmp, nang)
      radtmp(1:nang) = emf(1:nang) + radtmp(1:nang)
      CALL MXVV(rm, radtmp, rad, nang)
    END IF
    RETURN
  END SUBROUTINE MRGRAD
!
  SUBROUTINE MLREF(ra, rb, rm)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   MLREF
!
! PURPOSE:
!
!   Calculates multiple reflection term by truncated series.
!
! SYNTAX:
!
!   CALL MLREF(RA, RB)
!
! ARGUMENTS:
!
!   INPUTS/OUTPUTS:
!
!   RA  REAL(FP)  Atmosphere Reflection in front of observer
!   RB  REAL(FP)  Atmosphere Reflection behind observer
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Calculates multiple reflection term by truncated series.
!
!     REVISED:  11 OCTOBER 1993
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
    REAL(FP), DIMENSION(:), INTENT(IN) :: ra, rb
    REAL(FP), DIMENSION(:), INTENT(OUT) :: rm
    INTEGER :: nang
    INTEGER :: nang2
    INTEGER :: n2ang
!
    COMMON /pnlpar/ nang, n2ang, nang2
    save /pnlpar/
    REAL(FP) :: q1
    REAL(FP) :: q2
    REAL(FP) :: q3(mxang2)
    REAL(FP) :: rab
    dimension rab(mxang2), q1(mxang2), q2(mxang2)
    integer :: i
!
    CALL MXMV(ra, rb, rab, nang)
!
    q3(1:nang2)=rab(1:nang2)
    CALL MXMV(rab, q3, q1, nang)
    rm(1:nang2) = rab(1:nang2) + q1(1:nang2)
!
    CALL MXMV(rab, q1, q2, nang)
    rm(1:nang2) = rm(1:nang2) + q2(1:nang2)
!
    CALL MXMV(rab, q2, q1, nang)
    rm(1:nang2) = rm(1:nang2) + q1(1:nang2)
!
    CALL MXMV(rab, q1, q2, nang)
    CALL HOCOR(rm, q1, q2)
!
    ! to diagonal
    do i=1,nang2,nang+1
       rm(i)=rm(i)+1.
    end do
!
    RETURN
  END SUBROUTINE MLREF
!
  SUBROUTINE HOCOR(rm, q1, q2)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   HOCOR
!
! PURPOSE:
!
!   Correction for missing high order terms in MLREF
!
! SYNTAX:
!
!   CALL HOCOR(RM, Q1, Q2)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   Q1  REAL(FP)  TBD
!   Q2  REAL(FP)  TBD
!
!   INPUTS/OUTPUTS:
!
!   RM  REAL(FP)  TBD
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Correction for missing high order terms in MLREF
!
!     REVISED:  13 APRIL 1994
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    REAL(FP), DIMENSION(:), INTENT(INOUT) :: rm
    REAL(FP), DIMENSION(:), INTENT(IN) :: q1, q2
    INTEGER :: nang
    INTEGER :: nang2
    INTEGER :: n2ang
!
    COMMON /pnlpar/ nang, n2ang, nang2
    save /pnlpar/
    INTEGER :: l
    REAL(FP) :: epsil
    INTRINSIC MOD
    DATA epsil /1.e-16/
!
    DO l=1,nang2
      IF (MOD(l, nang) .NE. 0) rm(l) = rm(l) + q2(l)*q1(l)/(q1(l)-q2(l)+&
&          epsil)
    END DO
    RETURN
  END SUBROUTINE HOCOR
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     OSS (Optimal Spectral Sampling)
!     $Name:  $
!     $Id: oss_addbl.f90,v 1.6 2009/04/03 20:21:19 crichard Exp $
!     Copyright AER, Inc., 2002, 2003. All rights Reserved.
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     VECTOR OPERATIONS
!     REVISED:  04 OCTOBER 1993
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  SUBROUTINE MXM(a, nar, b, nac, c, nbc)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   MXM
!
! PURPOSE:
!
!   Product of two matrices
!
! SYNTAX:
!
!   CALL MXM(A, NAR, B, NAC, C, NBC)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   A    REAL(FP)     Generic array
!   NAR  INTEGER  Number of rows of matrix A
!   B    REAL(FP)     Generic array
!   NAC  INTEGER  Number of columns of matrix A
!   NBC  INTEGER  Number of columns of matrices B and C
!
!   INPUTS/OUTPUTS:
!
!   C    REAL(FP)     Generic array
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
! **********************************************************************
! *   CALCULATES PRODUCT OF TWO MATRICES--C=BA AND ASSUMES THE         *
! *   SKIP DISTANCE BETWEEN ELEMENTS TO BE 1.                          *
! *                                                                    *
! *   REVISED:  11-11-93 (D. PORTMAN)                                  *
! *                                                                    *
! *   ARGUMENTS:                                                       *
! *                                                                    *
! *   A       (input) Second matrix of product.                        *
! *   NAR     (input/output) Number of rows of matrices A and C.       *
! *   B       (input) First matrix of product.                         *
! *   NAC     (input) Number of columns of matrix A and the number of  *
! *           rows of matrix B.                                        *
! *   C       (output) Result matrix.                                  *
! *   NBC     (input/output) Number of columns of matrices B and C.    *
! *                                                                    *
! **********************************************************************
    INTEGER, INTENT(IN) :: nar, nac, nbc
    REAL(FP), DIMENSION(nar, nac), INTENT(IN) :: a
    REAL(FP), DIMENSION(nac, nbc), INTENT(IN) :: b
    REAL(FP), DIMENSION(nar, nbc), INTENT(OUT) :: c
    INTEGER :: i
    INTEGER :: j
    REAL(FP) :: sum
    INTEGER :: k
!
    DO i=1,nar
      DO j=1,nbc
        sum = 0.
        DO k=1,nac
          sum = sum + a(i, k)*b(k, j)
        END DO
        c(i, j) = sum
      END DO
    END DO
    RETURN
  END SUBROUTINE MXM
!
  SUBROUTINE SMINV(a, n, ierr)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   sMINV
!
! PURPOSE:
!
!   Invert a double precision matrix
!
! SYNTAX:
!
!   CALL sMINV(A, N, D, L, M)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   N  INTEGER  Dimension of array
!
!   INPUTS/OUTPUTS:
!
!   A  REAL(FP)     Generic array
!   D  REAL(FP)     Matrix determinant
!   L  REAL(FP)     Workspace vector
!   M  REAL(FP)     Workspace vector
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
! **********************************************************************
! ROUTINE          DMINV   SUBROUTINE  *****  [EYRE.RETCOF]
!
! PURPOSE          TO INVERT A DOUBLE PRECISION MATRIX
!
! VERSION          3.00,150385,J.R.EYRE
!
! DESCRIPTION      ROUTINE TO INVERT A DOUBLE PRECISION MATRIX.
!                  METHOD: THE STANDARD GAUSS-JORDAN METHOD IS USED.
!                  THE DETERMINANT IS ALSO CALCULATED. A DETERMINANT
!                  OF ZERO INDICATES THAT THE MATRIX IS SINGULAR.
!
! ARGUMENTS        A      R*8 INPUT MATRIX, DESTROYED IN COMPUTATION AND
!                             REPLACED BY RESULTANT INVERSE.
!                  N      I*4 ORDER OF MATRIX A
!                  D      R*8 RESULTANT DETERMINANT
!                  L      R*8 WORK VECTOR OF LENGTH N
!                  M      R*8 WORK VECTOR OF LENGTH N
!
!     ..................................................................
!
    INTEGER, INTENT(IN) :: n
    REAL(FP), DIMENSION(n*n), INTENT(INOUT) :: a
    integer, intent(out) :: ierr
    REAL(FP) :: d
    REAL(FP), DIMENSION(n) :: l, m
!
    REAL(FP) :: biga, hold
    INTEGER :: nk
    INTEGER :: k
    INTEGER :: kk
    INTEGER :: j
    INTEGER :: i
    INTEGER :: ij
    INTEGER :: ki
    INTEGER :: ji
    INTEGER :: jk
    INTEGER :: ik
    INTEGER :: kj
    INTEGER :: jq
    INTEGER :: jr
    INTRINSIC ABS
    REAL(FP) :: abs2
    REAL(FP) :: abs1
!
!     ..................................................................
!
!        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE
!        C$ IN COLUMNS 1-2 SHOULD BE REMOVED FROM THE DOUBLE PRECISION
!        STATEMENT WHICH FOLLOWS.
!
!$      REAL(Double) A,D,BIGA,HOLD
!        REAL(FP) A,D,BIGA,HOLD,L,M
!
!        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS
!        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS
!        ROUTINE.
!
!        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO
!        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  ABS IN STATEMENT
!        10 MUST BE CHANGED TO DABS.
!
! **********************************************************************
!
!        SEARCH FOR LARGEST ELEMENT
!
    ierr = 0
    d = 1.0
    DO k=1,n
      nk = (k-1)*n
      l(k) = k
      m(k) = k
      kk = nk + k
      biga = a(kk)
      DO j=k,n
        DO i=k,n
          ij = (j-1)*n + i
          IF (biga .GE. 0.) THEN
            abs1 = biga
          ELSE
            abs1 = -biga
          END IF
          IF (a(ij) .GE. 0.) THEN
            abs2 = a(ij)
          ELSE
            abs2 = -a(ij)
          END IF
!$ 10           IF(DABS(BIGA)-DABS(A(IJ))) 15,20,20
          IF (abs1 - abs2 .LT. 0) THEN
            biga = a(ij)
            l(k) = i
            m(k) = j
          END IF
        END DO
      END DO
      j = l(k)
!
!        INTERCHANGE ROWS
!
      if (j-k .ne. 0) then
!!$      IF (.NOT.j - k .LT. 0) THEN
!!$        IF (.NOT.j - k .EQ. 0) THEN
          ki = k - n
          DO i=1,n
            ki = k + (i-1)*n
            hold = -a(ki)
            ji = ki - k + j
            a(ki) = a(ji)
            a(ji) = hold
          END DO
!!$        END IF
      END IF
      i = m(k)
!
!        INTERCHANGE COLUMNS
!
      IF (i - k .ne. 0) THEN
         DO j=1,n
            jk = nk + j
            ji = n*(i-1) + j
            hold = -a(jk)
            a(jk) = a(ji)
            a(ji) = hold
         END DO
      END IF
!
!        DIVIDE COLUMN BY MINUS PIVOT (VALUE OF PIVOT ELEMENT IS
!        CONTAINED IN BIGA)
!
      IF (biga .EQ. 0) THEN
         d = 0.0
         ierr = 1
         RETURN
      END IF
      DO i=1,n
        if (i-k .ne. 0) then
!!$        IF (.NOT.i - k .LT. 0) THEN
!!$          IF (i - k .EQ. 0) GOTO 55
!!$        END IF
           ik = nk + i
           a(ik) = a(ik)/(-biga)
        endif
!!$ 55     CONTINUE
      END DO
!
!        REDUCE MATRIX
!
      DO i=1,n
        ik = nk + i
        hold = a(ik)
        DO j=1,n
          if (i-k .ne. 0 .and. j-k .ne. 0) then
!!$          IF (.NOT.i - k .LT. 0) THEN
!!$            IF (i - k .EQ. 0) GOTO 65
!!$          END IF
!!$          IF (.NOT.j - k .LT. 0) THEN
!!$            IF (j - k .EQ. 0) GOTO 65
!!$          END IF
             ij = (j-1)*n + i
             kj = ij - i + k
             a(ij) = hold*a(kj) + a(ij)
          endif
!!$ 65       CONTINUE
        END DO
      END DO
!
!        DIVIDE ROW BY PIVOT
!
      kj = k - n
      DO j=1,n
        kj = kj + n
        if (j-k .ne. 0) then
!!$        IF (.NOT.j - k .LT. 0) THEN
!!$          IF (j - k .EQ. 0) GOTO 75
!!$        END IF
           a(kj) = a(kj)/biga
        endif
!!$ 75     CONTINUE
      END DO
!
!        PRODUCT OF PIVOTS
!
      d = d*biga
!
!        REPLACE PIVOT BY RECIPROCAL
!
      a(kk) = 1.0/biga
    END DO
!
!        FINAL ROW AND COLUMN INTERCHANGE
!
    do k=n,1,-1
!!$    k = n
!!$ 100 k = k - 1
!!$    IF (.NOT.k .LT. 0) THEN
!!$      IF (.NOT.k .EQ. 0) THEN
       i = l(k)
       if (i-k .ne. 0) then
!!$       IF (.NOT.i - k .LT. 0) THEN
!!$          IF (.NOT.i - k .EQ. 0) THEN
          jq = n*(k-1)
          jr = n*(i-1)
          DO j=1,n
             jk = jq + j
             hold = a(jk)
             ji = jr + j
             a(jk) = -a(ji)
             a(ji) = hold
          END DO
!!$          END IF
       END IF
       j = m(k)
       if (j-k .ne. 0) then
!!$        IF (.NOT.j - k .LT. 0) THEN
!!$          IF (.NOT.j - k .EQ. 0) THEN
          ki = k - n
          DO i=1,n
             ki = ki + n
             hold = a(ki)
             ji = ki - k + j
             a(ki) = -a(ji)
             a(ji) = hold
          END DO
!!$          END IF
       END IF
    end do

!!$        GOTO 100
!!$      END IF
!!$    END IF
    RETURN
  END SUBROUTINE SMINV
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     VECTORIZED CHAIN MATRIX AND MATRIX/VECTOR OPERATIONS
!     REVISED:  04 OCTOBER 1993
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  SUBROUTINE MXMV(a, b, c, isiz)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   MXMV
!
! PURPOSE:
!
!   Product second matrix by first matrix, where both are square.
!
! SYNTAX:
!
!   CALL MXMV(A, B, C, ISIZ)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   A     REAL(FP)     Generic array
!   B     REAL(FP)     Generic array
!   ISIZ  INTEGER  Dimension of array
!
!   INPUTS/OUTPUTS:
!
!   C     REAL(FP)     Generic array
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
    INTEGER, INTENT(IN) :: isiz
    REAL(FP), DIMENSION(isiz, isiz), INTENT(IN) :: a, b
    REAL(FP), DIMENSION(isiz, isiz), INTENT(OUT) :: c
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
!
    DO i=1,isiz
      DO j=1,isiz
        c(j, i) = 0.
        DO k=1,isiz
          c(j, i) = c(j, i) + a(k, i)*b(j, k)
        END DO
      END DO
    END DO
    RETURN
  END SUBROUTINE MXMV
!
  SUBROUTINE DMXMV(a, b, c, isiz)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   DMXMV
!
! PURPOSE:
!
!   Product square matrix by diagonal matrix (stored as a vector). Result is
!   matrix.
!
! SYNTAX:
!
!   CALL DMXMV(A, B, C, ISIZ)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   A     REAL(FP)     Generic array
!   B     REAL(FP)     Generic array
!   ISIZ  INTEGER  Dimension of array
!
!   INPUTS/OUTPUTS:
!
!   C     REAL(FP)     Generic array
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
    INTEGER, INTENT(IN) :: isiz
    REAL(FP), DIMENSION(isiz), INTENT(IN) :: a
    REAL(FP), DIMENSION(isiz, isiz), INTENT(IN) :: b
    REAL(FP), DIMENSION(isiz, isiz), INTENT(INOUT) :: c
    INTEGER :: i
    INTEGER :: j
!
    DO i=1,isiz
      DO j=1,isiz
        c(j, i) = a(i)*b(j, i)
      END DO
    END DO
    RETURN
  END SUBROUTINE DMXMV
!
  SUBROUTINE MXDMV(a, b, c, isiz)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   MXDMV
!
! PURPOSE:
!
!   Product diagonal matrix (stored as a vector) by a square matrix. Result is
!   matrix.
!
! SYNTAX:
!
!   CALL MXDMV(A, B, C, ISIZ)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   A     REAL(FP)     Generic array
!   B     REAL(FP)     Generic array
!   ISIZ  INTEGER  Dimension of array
!
!   INPUTS/OUTPUTS:
!
!   C     REAL(FP)     Generic array
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
    INTEGER, INTENT(IN) :: isiz
    REAL(FP), DIMENSION(isiz), INTENT(IN) :: b
    REAL(FP), DIMENSION(isiz, isiz), INTENT(IN) :: a
    REAL(FP), DIMENSION(isiz, isiz), INTENT(INOUT) :: c
    INTEGER :: i
    INTEGER :: j
!
    DO i=1,isiz
      DO j=1,isiz
        c(j, i) = a(j, i)*b(j)
      END DO
    END DO
    RETURN
  END SUBROUTINE MXDMV
!
  SUBROUTINE MXVV(a, b, c, isiz)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   MXVV
!
! PURPOSE:
!
!   Product vector (transposed) by matrix. Result is vector.
!
! SYNTAX:
!
!   CALL MXVV(A, B, C, ISIZ)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   A     REAL(FP)     Generic array
!   B     REAL(FP)     Generic array
!   ISIZ  INTEGER  Dimension of array
!
!   INPUTS/OUTPUTS:
!
!   C     REAL(FP)     Generic array
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
    INTEGER, INTENT(IN) :: isiz
    REAL(FP), DIMENSION(isiz), INTENT(IN) :: b
    REAL(FP), DIMENSION(isiz, isiz), INTENT(IN) :: a
    REAL(FP), DIMENSION(isiz), INTENT(OUT) :: c
    INTEGER :: i
    INTEGER :: j
!
    DO i=1,isiz
      c(i) = 0.
      DO j=1,isiz
        c(i) = c(i) + a(j, i)*b(j)
      END DO
    END DO
    RETURN
  END SUBROUTINE MXVV
  SUBROUTINE LEPOLY_OSS(nmu, oylm, maxmu, twonm1, mu, ylm)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   LEPOLY_OSS
!
! PURPOSE:
!
!   Computes normalized associated Legendre polynomial, defined in terms of
!   associated Legendre polynomial
!
! SYNTAX:
!
!   CALL LEPOLY_OSS(NMU, oYLM, MAXMU, TWONM1, MU, YLM)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   NMU     INTEGER  Number of arguments of YLM
!   oYLM    INTEGER  Order of YLM
!   MAXMU   INTEGER  First dimension of YLM
!   TWONM1  INTEGER  Max Degree of YLM
!   MU      REAL(FP)     YLM arguments
!
!   INPUTS/OUTPUTS:
!
!   YLM     REAL(FP)     associated Legendre polynomials
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!       ORIGIN:     DISORT (Stamnes etal.) V1.1
!
!       COMPUTES THE NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL,
!       DEFINED IN TERMS OF THE ASSOCIATED LEGENDRE POLYNOMIAL
!       PLM = P-SUB-L-SUPER-oYLM AS
!
!             YLM(MU) = SQRT( (L-oYLM)!/(L+oYLM)! ) * PLM(MU)
!
!       FOR FIXED ORDER -oYLM- AND ALL DEGREES FROM L = oYLM TO TWONM1.
!       WHEN oYLM.GT.0, ASSUMES THAT Y-SUB(oYLM-1)-SUPER(oYLM-1) IS AVAILABLE
!       FROM A PRIOR CALL TO THE ROUTINE.
!
!       REFERENCE: Dave, J.V. and B.H. Armstrong, Computations of
!                  High-Order Associated Legendre Polynomials,
!                  J. Quant. Spectrosc. Radiat. Transfer 10,
!                  557-562, 1970.  (hereafter D/A)
!
!       METHOD: Varying degree recurrence relationship.
!
!       NOTE 1: The D/A formulas are transformed by
!               setting  oYLM = n-1; L = k-1.
!       NOTE 2: Assumes that routine is called first with  oYLM = 0,
!               then with  oYLM = 1, etc. up to  oYLM = TWONM1.
!       NOTE 3: Loops are written in such a way as to vectorize.
!
!  I N P U T     V A R I A B L E S:
!
!       NMU    :  NUMBER OF ARGUMENTS OF -YLM-
!       oYLM      :  ORDER OF -YLM-
!       MAXMU  :  FIRST DIMENSION OF -YLM-
!       TWONM1 :  MAX DEGREE OF -YLM-
!       MU(I)  :  I = 1 TO NMU, ARGUMENTS OF -YLM-
!       IF oYLM.GT.0, YLM(oYLM-1,I) FOR I = 1 TO NMU IS REQUIRED
!
!  O U T P U T     V A R I A B L E:
!
!       YLM(L,I) :  L = oYLM TO TWONM1, NORMALIZED ASSOCIATED LEGENDRE
!                   POLYNOMIALS EVALUATED AT ARGUMENT -MU(I)-
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    INTEGER, INTENT(IN) :: nmu, oylm, maxmu, twonm1
    REAL(FP), DIMENSION(:), INTENT(IN) :: mu
    REAL(FP), DIMENSION(0:maxmu, *), INTENT(INOUT) :: ylm
    INTEGER :: maxsqt
    PARAMETER (maxsqt=1000)
    REAL(FP) :: sqt(maxsqt)
    LOGICAL :: pass1
    SAVE sqt, pass1
    INTEGER :: ns
    INTEGER :: i
    INTEGER :: l
    REAL(FP) :: tmp1
    REAL(FP) :: tmp2
    REAL(Double) :: arg1
    REAL(FP) :: arg10
    REAL(FP) :: result1
    INTRINSIC FLOAT
    INTRINSIC SQRT
    DATA pass1 /.true./
!
    IF (pass1) THEN
      pass1 = .false.
      DO ns=1,maxsqt
        arg1 = FLOAT(ns)
        sqt(ns) = SQRT(arg1)
      END DO
    END IF
!
    IF (2*twonm1 .GT. maxsqt) THEN
      STOP
    ELSE
!
      IF (oylm .EQ. 0) THEN
!
!     UPWARD RECURRENCE FOR ORDINARY LEGENDRE POLYNOMIALS
!
        DO i=1,nmu
          ylm(0, i) = 1.
          ylm(1, i) = mu(i)
        END DO
!
        DO l=2,twonm1
          DO i=1,nmu
            ylm(l, i) = ((2*l-1)*mu(i)*ylm(l-1, i)-(l-1)*ylm(l-2, i))/l
          END DO
        END DO
      ELSE
        DO i=1,nmu
!
!     Y-SUB-M-SUPER-M; DERIVED FROM D/A EQS. (11,12)
!
          arg10 = 1. - mu(i)**2
          result1 = SQRT(arg10)
          ylm(oylm, i) = -(sqt(2*oylm-1)/sqt(2*oylm)*result1*ylm(oylm-1&
&            , i))
!
!     Y-SUB-(M+1)-SUPER-M; DERIVED FROM D/A EQS. (13,14) USING EQS. (11,12)
!
          ylm(oylm+1, i) = sqt(2*oylm+1)*mu(i)*ylm(oylm, i)
        END DO
!
!     UPWARD RECURRENCE; D/A EQ. (10)
!
        DO l=oylm+2,twonm1
          tmp1 = sqt(l-oylm)*sqt(l+oylm)
          tmp2 = sqt(l-oylm-1)*sqt(l+oylm-1)
          DO i=1,nmu
            ylm(l, i) = ((2*l-1)*mu(i)*ylm(l-1, i)-tmp2*ylm(l-2, i))/&
&              tmp1
          END DO
        END DO
      END IF
      RETURN
    END IF
  END SUBROUTINE LEPOLY_OSS
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     OSS (Optimal Spectral Sampling)
!     $Name:  $
!     $Id: oss_addbl.f90,v 1.6 2009/04/03 20:21:19 crichard Exp $
!     Copyright AER, Inc., 2002, 2003. All rights Reserved.
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  FUNCTION WDRADDT(vn, t)
    IMPLICIT NONE
!<f90Function>**********************************************************
!
! NAME:
!
!   wdraddt
!
! PURPOSE:
!
!   Derivative of radiance with respect to temperature.  Duplicate of DBDT.
!
! SYNTAX:
!
!   Results=wdraddt(vn, t)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   vn       REAL(Double)  Wavenumber, units of cm^-1
!   t        REAL(FP)    Temperature
!
!   * OPTIONAL
!
! RETURN:
!
!     REAL(FP)
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>
    REAL(Double), INTENT(IN) :: vn
    REAL(FP), INTENT(IN) :: t
    REAL(FP) :: wdraddt
    REAL(Double), parameter :: c1=RADCN1
    REAL(Double), parameter :: c2=RADCN2
    REAL(Double) :: t_1
    REAL(Double) :: ct1v1
    REAL(Double) :: ct2v1
    REAL(Double) :: ct3v1
    REAL(Double) :: bplan
    INTRINSIC EXP
!
    t_1 = 1/t
!
    ct1v1 = c1*vn**3
    ct2v1 = c2*vn*t_1
    ct3v1 = EXP(ct2v1)
!
    bplan = ct1v1/(ct3v1-1.)
    wdraddt = bplan/(ct3v1-1.)*ct2v1*ct3v1*t_1
!
    RETURN
  END FUNCTION WDRADDT
  FUNCTION WNBRIT(vn, rad)
    IMPLICIT NONE
!<f90Function>**********************************************************
!
! NAME:
!
!   WNBRIT
!
! PURPOSE:
!
!   Equivalent blackbody temperature of radiance (inverse Planck function).
!
! SYNTAX:
!
!   Results=WNBRIT(VN, RAD)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   VN      REAL(Double)  Wavenumber, units of cm^-1
!   RAD     REAL(FP)    Radiance
!
!   * OPTIONAL
!
! RETURN:
!
!     REAL(FP)
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>
! $ RADIANCE TO BRIGHTNESS TEMPERATURE (HMW)
! * 'NEW' PLANCK'S CONSTANT, VELOCITY OF LIGHT, BOLTZMANN'S CONSTANT
    REAL(FP), INTENT(IN) :: rad
    REAL(Double), INTENT(IN) :: vn
    REAL(FP) :: wnbrit
    REAL(Double), parameter :: c1=RADCN1
    REAL(Double), parameter :: c2=RADCN2
    REAL(Double) :: f1
    REAL(Double) :: f2
    REAL(Double) :: r
    REAL(Double) :: tbb
    INTRINSIC DLOG

!
    f1 = c1*vn**3
    f2 = c2*vn
    r = rad
    tbb = f2/DLOG(f1/r+1.d0)
    wnbrit = tbb
    RETURN
  END FUNCTION WNBRIT
  FUNCTION DBDT(vn, t)
    IMPLICIT NONE
!<f90Function>**********************************************************
!
! NAME:
!
!   DBDT
!
! PURPOSE:
!
!   Derivative of radiance with respect to temperature.  Duplicate of wdraddt
!
! SYNTAX:
!
!   Results=DBDT(VN, T)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   VN    REAL(Double)  Wavenumber, units of cm^-1
!   T     REAL(FP)    Temperature
!
!   * OPTIONAL
!
! RETURN:
!
!     REAL(FP)
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>
    REAL(FP), INTENT(IN) :: t
    REAL(Double), INTENT(IN) :: vn
    REAL(FP) :: dbdt
    REAL(Double), parameter :: c1=RADCN1
    REAL(Double), parameter :: c2=RADCN2
    REAL(Double) :: t_1
    REAL(Double) :: ct1v1
    REAL(Double) :: ct2v1
    REAL(Double) :: ct3v1
    REAL(Double) :: bplan
    INTRINSIC EXP
!
    t_1 = 1/t
!
    ct1v1 = c1*vn**3
    ct2v1 = c2*vn*t_1
    ct3v1 = EXP(ct2v1)
!
    bplan = ct1v1/(ct3v1-1.)
    dbdt = bplan/(ct3v1-1.)*ct2v1*ct3v1*t_1
!
    RETURN
  END FUNCTION DBDT
  FUNCTION WNMRAD(vn, tran, tem, ts, ls)
    IMPLICIT NONE
!<f90Function>**********************************************************
!
! NAME:
!
!   WNMRAD
!
! PURPOSE:
!
!   Radiance of absorbing/emitting atmosphere with a blackbody surface
!   background.
!
! SYNTAX:
!
!   Results=WNMRAD(VN, TRAN, TEM, TS, LS)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   VN      REAL(Double)   Wavenumber, units of cm^-1
!   TRAN    REAL(FP)     transmittance
!   TEM     REAL(FP)     level temperature
!   TS      REAL(FP)     Surface skin temperature
!   LS      INTEGER  index of the surface level
!
!   * OPTIONAL
!
! RETURN:
!
!     REAL(FP)
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>
! $ MONOCHROMATIC RADIANCE CALCULATION (HMW)
! Units are mw*cm/m2/sr
    REAL(Double), INTENT(IN) :: vn
    REAL(FP), DIMENSION(:), INTENT(IN) :: tran, tem
    REAL(FP), INTENT(IN) :: ts
    INTEGER, INTENT(IN) :: ls
    REAL(FP) :: wnmrad
    REAL(FP) :: tau1
    REAL(FP) :: t1
    REAL(FP) :: b1
    REAL(FP) :: rad
    INTEGER :: i
    REAL(FP) :: tau2
    REAL(FP) :: t2
    REAL(FP) :: b2
    REAL(FP) :: bsts
    REAL(FP) :: result1
!
    tau1 = tran(1)
    t1 = tem(1)
    b1 = WNPLAN(vn, t1)
    rad = 0.
    DO i=2,ls
      tau2 = tran(i)
      t2 = tem(i)
      b2 = WNPLAN(vn, t2)
      rad = rad + .5*(b1+b2)*(tau1-tau2)
      tau1 = tau2
      b1 = b2
    END DO
    result1 = WNPLAN(vn, ts)
    bsts = result1*tau1
    rad = rad + bsts
    wnmrad = rad
    RETURN
  END FUNCTION WNMRAD
  FUNCTION WNPLAN(vn, tem)
    IMPLICIT NONE
!<f90Function>**********************************************************
!
! NAME:
!
!   WNPLAN
!
! PURPOSE:
!
!   Planck function, with internal calculations in double precision.
!
! SYNTAX:
!
!   Results=WNPLAN(VN, TEM)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   VN      REAL(Double)  Wavenumber, units of cm^-1
!   TEM     REAL(FP)    level temperature
!
!   * OPTIONAL
!
! RETURN:
!
!     REAL(FP)
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>
! $ TEMPERATURE TO PLANCK RADIANCE (HMW)
! $ 'NEW' PLANCK'S CONSTANT, VELOCITY OF LIGHT, BOLTZMANN'S CONSTANT
! Units are mw*cm/m2/sr
    REAL(FP), INTENT(IN) :: tem
    REAL(Double), INTENT(IN) :: vn
    REAL(FP) :: wnplan
    REAL(Double), parameter :: c1=RADCN1
    REAL(Double), parameter :: c2=RADCN2
    REAL(Double) :: f1
    REAL(Double) :: f2
    REAL(Double) :: t
    REAL(Double) :: rad

!
!---To avoid dividing by zero
    IF (tem .EQ. 0.) THEN
      wnplan = 0.
      RETURN
    ELSE
      f1 = c1*vn**3
      f2 = c2*vn
      t = tem
      rad = f1/(EXP(f2/t)-one_dp)
      wnplan = rad
!---test
!WNPLAN=TEM
!------
      RETURN
    END IF
  END FUNCTION WNPLAN
  SUBROUTINE QGAUSN(oquad, gmu, gwt)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   QGAUSN
!
! PURPOSE:
!
!   Compute weights and abscissae for ordinary Gaussian quadrature on the
!   interval (0,1)
!
! SYNTAX:
!
!   CALL QGAUSN(oQuad, GMU, GWT)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   oQuad  INTEGER  order of quadrature rule
!
!   INPUTS/OUTPUTS:
!
!   GMU    REAL(FP)     array of abscissae
!   GWT    REAL(FP)     array of weights
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
!       Compute weights and abscissae for ordinary Gaussian quadrature
!       on the interval (0,1);  that is, such that
!           sum(i=1 to oQuad) ( GWT(i) f(GMU(i)) )
!       is a good approximation to
!           integral(0 to 1) ( f(x) dx )
!   INPUT :    oQuad       order of quadrature rule
!   OUTPUT :  GMU(I)   array of abscissae (I = 1 TO oQuad)
!             GWT(I)   array of weights (I = 1 TO oQuad)
!   REFERENCE:  Davis, P.J. and P. Rabinowitz, Methods of Numerical
!                   Integration, Academic Press, New York, pp. 87, 1975
!   METHOD:  Compute the abscissae as roots of the Legendre
!            polynomial P-sub-M using a cubically convergent
!            refinement of Newton's method.  Compute the
!            weights from EQ. 2.7.3.8 of Davis/Rabinowitz.  Note
!            that Newton's method can very easily diverge; only a
!            very good initial guess can guarantee convergence.
!            The initial guess used here has never led to divergence
!            even for oQuad up to 1000.
!   ACCURACY:  relative error no better than TOL or computer
!              precision (machine epsilon), whichever is larger
!   INTERNAL VARIABLES:
!    ITER      : number of Newton Method iterations
!    MAXIT     : maximum allowed iterations of Newton Method
!    PM2,PM1,P : 3 successive Legendre polynomials
!    PPR       : derivative of Legendre polynomial
!    P2PRI     : 2nd derivative of Legendre polynomial
!    TOL       : convergence criterion for Legendre poly root iteration
!    X,XI      : successive iterates in cubically-convergent version
!                of Newtons Method (seeking roots of Legendre poly.)
!   Called by- SETDIS, SURFAC
!   Calls- D1MACH, ERRMSG
! +-------------------------------------------------------------------+
!     .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: oquad
!     ..
!     .. Array Arguments ..
    REAL(FP), DIMENSION(oquad), INTENT(INOUT) :: gmu(oquad), gwt(oquad)
!     ..
!     .. Local Scalars ..
    INTEGER :: iter, k, lim, maxit, nn, np1
    REAL(FP) :: cona, pi, t
    REAL(Double) :: en, nnp1, one, p, p2pri, pm1, pm2, ppr, prod, &
&    tmp, tol, two, x, xi
!     ..
!     .. Intrinsic Functions ..
    INTRINSIC ABS, ASIN, COS, FLOAT, MOD, TAN
!     ..
    SAVE pi, tol
    REAL(Double) :: result1
    REAL(FP) :: arg1
    REAL(Double) :: abs1
    DATA pi /0.0/
    DATA maxit /1000/
    DATA one /1.d0/
    DATA two /2.d0/
    IF (pi .EQ. 0.0) THEN
      pi = 2.*ASIN(1.0)
      result1 = D1MACH(4)
      tol = 10.*result1
    END IF
    IF (oquad .LT. 1) CALL ERRMSG('QGAUSN--Bad value of oQuad', .true.)
    IF (oquad .EQ. 1) THEN
      gmu(1) = 0.5
      gwt(1) = 1.0
      RETURN
    ELSE
      en = oquad
      np1 = oquad + 1
      nnp1 = oquad*np1
      cona = FLOAT(oquad-1)/(8*oquad**3)
      lim = oquad/2
      DO k=1,lim
!                                        ** Initial guess for k-th root
!                                        ** of Legendre polynomial, from
!                                        ** Davis/Rabinowitz (2.7.3.3a)
        t = (4*k-1)*pi/(4*oquad+2)
        arg1 = t + cona/TAN(t)
        x = COS(arg1)
        iter = 0
!                                        ** Upward recurrence for
!                                        ** Legendre polynomials
 10     iter = iter + 1
        pm2 = one
        pm1 = x
        DO nn=2,oquad
          p = ((2*nn-1)*x*pm1-(nn-1)*pm2)/nn
          pm2 = pm1
          pm1 = p
        END DO
!                                              ** Newton Method
        tmp = one/(one-x**2)
        ppr = en*(pm2-x*p)*tmp
        p2pri = (two*x*ppr-nnp1*p)*tmp
        xi = x - p/ppr*(one+p/ppr*p2pri/(two*ppr))
        IF (xi - x .GE. 0.) THEN
          abs1 = xi - x
        ELSE
          abs1 = -(xi-x)
        END IF
!                                              ** Check for convergence
        IF (abs1 .GT. tol) THEN
          IF (iter .GT. maxit) CALL ERRMSG('QGAUSN--max iteration count'&
&                                     , .true.)
          x = xi
          GOTO 10
        ELSE
          gmu(k) = -x
!                             ** Iteration finished--calculate weights,
!                             ** abscissae for (-1,1)
          gwt(k) = two/(tmp*(en*pm2)**2)
          gmu(np1-k) = -gmu(k)
          gwt(np1-k) = gwt(k)
        END IF
      END DO
!                                    ** Set middle abscissa and weight
!                                    ** for rules of odd order
      IF (MOD(oquad, 2) .NE. 0) THEN
        gmu(lim+1) = 0.0
        prod = one
        DO k=3,oquad,2
          prod = prod*k/(k-1)
        END DO
        gwt(lim+1) = two/prod**2
      END IF
!                                        ** Convert from (-1,1) to (0,1)
      DO k=1,oquad
        gmu(k) = 0.5*gmu(k) + 0.5
        gwt(k) = 0.5*gwt(k)
      END DO
      RETURN
    END IF
  END SUBROUTINE QGAUSN
  FUNCTION D1MACH(i)
    IMPLICIT NONE
!<f90Function>**********************************************************
!
! NAME:
!
!   D1MACH
!
! PURPOSE:
!
!   Set Double-precision machine constants
!
! SYNTAX:
!
!   Results=D1MACH(I)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   I       INTEGER  Integer flag for setting machine constant
!
!   * OPTIONAL
!
! RETURN:
!
!     REAL(Double)
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>
!  Double-precision machine constants (see R1MACH for documentation).
!  By default, returns values appropriate for a computer with IEEE
!  arithmetic.  This is an abbreviated version of a routine widely
!  used for 20+ years by numerical analysts.  Most of the values in
!  the original version pertain to computers which went to computer
!  heaven years ago and are of little if any interest.
!
!  If the values herein do not work for any reason, just look in
!  your Fortran manual for the correct values (usually in the part
!  discussing representations of numbers) and insert them. The exact
!  values are not that important; they can be a factor of 2-3 off
!  without causing any harm.
!  Only I = 1,2,4 is actually used by DISORT.
!  This routine is superseded in Fortran-90 by the intrinsic numeric
!  inquiry functions HUGE(1.D0), TINY(1.D0), and EPSILON(1.D0).
!  The original version can be found on NetLib (search by name):
!      http://www.netlib.org/
! ====================================================================
    INTEGER, INTENT(IN) :: i
    REAL(Double) :: d1mach
!EXTERNAL  ERRMSG
    IF (i .EQ. 1) THEN
      d1mach = 2.3d-308
!        D1MACH = TINY(1.D0)
    ELSE IF (i .EQ. 2) THEN
      d1mach = 1.7d+308
!        D1MACH = HUGE(1.D0)
    ELSE IF (i .EQ. 4) THEN
      d1mach = 2.3d-16
!        D1MACH = EPSILON(1.D0)
    ELSE
       d1mach=0. !to suppress compiler warning
      CALL ERRMSG('D1MACH--argument incorrect', .true.)
    END IF
    RETURN
  END FUNCTION D1MACH
  SUBROUTINE ERRMSG(messag, fatal)
    IMPLICIT NONE
!<f90Subroutine>********************************************************
!
! NAME:
!
!   ErrMsg
!
! PURPOSE:
!
!   Print out a warning or error message
!
! SYNTAX:
!
!   CALL ErrMsg(MESSAG, FATAL)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   MESSAG  CHAR     text of message
!   FATAL   LOGICAL  Flag for fatal error
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
!        Print out a warning or error message;  abort if error
    LOGICAL, INTENT(IN) :: fatal
    CHARACTER(len=*), INTENT(IN) :: messag
!
    LOGICAL :: msglim
    INTEGER :: maxmsg, nummsg
    SAVE maxmsg, nummsg, msglim
    DATA nummsg /0/
    DATA maxmsg /100/
    DATA msglim /.false./
    IF (fatal) THEN
      WRITE(*, '(/,2A,/)') ' ******* ERROR >>>>>>  ', messag
      STOP
    ELSE
      nummsg = nummsg + 1
      IF (msglim) THEN
        RETURN
      ELSE
        IF (nummsg .LE. maxmsg) THEN
          WRITE(*, '(/,2A,/)') ' ******* WARNING >>>>>>  ', messag
        ELSE
          WRITE(*, 99)
          msglim = .true.
        END IF
        RETURN
      END IF
    END IF
 99 FORMAT( //,' >>>>>>  TOO MANY WARNING MESSAGES --  ', &
&           'They will no longer be printed  <<<<<<<', // )
  END SUBROUTINE ERRMSG

end module subs_oss_addbl
