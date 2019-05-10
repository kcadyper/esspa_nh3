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

MODULE CHKValid

! <f90Module>***********************************************************
!
! NAME:
!
!   CHKValid
!
! PURPOSE:
!
!   Subprograms related to checking and constraining the consistency of data
!   with physical limits.
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

  use ncdf_module, only : MAX_NAME_LENGTH
  use StateIndexModule
  use constants, only : MISSING_REAL,ratioMwt
  use MetFunctions, only : svp,mrFromVp
  implicit none
  PRIVATE
  ! loadDone implementation provides the option to call loadGuessLimits() in advance
  logical              :: loadGuessLimitsDone=.false.
  logical              :: loadBoundsDone=.false.
  !---GuessLimits namelist
  real                 :: maxcldbasebak,htol,tmax,tmin
  real                 :: ctoplow,icetoplow,deffLiqMin,deffIceMin
  real                 :: deffLiqMax,deffIceMax
  real                 :: thickmin,Pmin4SaturChk,liqCov,iceCov
  !---geoBounds namelist
  type physicalBounds
     character (len=MAX_NAME_LENGTH) :: ID
     real    :: lower, upper
  end type physicalBounds
  type(physicalBounds) :: &
       TempBnd = physicalBounds("temp",50.,400.), &
       TskinBnd = physicalBounds("tskin",50.,400.), &
       PsfcBnd = physicalBounds("psfc",100.,1200.), &
       H2OBnd  = physicalBounds("h2o",0.,0.050) ! g/g
  REAL, PARAMETER :: absoluteLimitSupersaturation=4.
  PUBLIC :: setPropertiesGuessLimits, loadGuessLimits, chkges, chkRad, &
       retrQC, getChiSq, checkProfile, setPropertiesBounds, loadBounds, &
       checkBounds, queryLowerBound, queryUpperBound, physicalBounds
  interface checkBounds
     module procedure checkBoundsVector
     module procedure checkBoundsScalar
  end interface
CONTAINS

  !--
  subroutine chkges(U_algconfig,F_algconfig, &
       prof,ifail,Sx,SxLiq,SxIce,nlev,press, &
       IG,NG,IR,NR,ih2o,debug, invCldCov_in)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   chkges
!
! PURPOSE:
!
!   Check and constraining the consistency of elements of a state vector with
!   physical limits.
!
! SYNTAX:
!
!   CALL chkges(U_algconfig, F_algconfig, prof, ifail, Sx, SxLiq,
!      SxIce, nlev, press, IG, NG, IR, NR, ih2o, debug, invCldCov_in)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   U_algconfig   INTEGER             Unit number for algorithm configuration
!                                     data file
!   F_algconfig   CHAR                File path for algorithm configuration
!                                     data
!   SxLiq         REAL                Background error covariance for liquid
!                                     cloud
!   SxIce         REAL                Background error covariance for ice
!                                     cloud
!   nlev          INTEGER             Number of atmospheric levels
!   press          REAL                Pressure on atmospheric grid levels
!   IG            TYPE(STATEINDEX_T)  Starting indices for sections of
!                                     geophysical state vector
!   NG            TYPE(STATEINDEX_T)  Number of elements for sections of
!                                     geophysical state vector
!   IR            TYPE(STATEINDEX_T)  Starting indices for sections of
!                                     retrieval state vector
!   NR            TYPE(STATEINDEX_T)  Number of retrieved elements for
!                                     sections of retrieval state vector
!   ih2o          INTEGER             Index for water vapor in molecules part
!                                     of state vector
!   debug         LOGICAL             Flag for debugging mode
!   invCldCov_in* LOGICAL             Flag for using inverse covariance
!                                     matrix in retrieval
!
!   OUTPUTS:
!
!   ifail         INTEGER             Flag that test was failed
!
!   INPUTS/OUTPUTS:
!
!   prof          REAL                Atmospheric profile
!   Sx            REAL                Background error covariance
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>


    !---In/Out variables
    INTEGER             ,INTENT(in)    :: U_algconfig
    CHARACTER(len=*)    ,INTENT(in)    :: F_algconfig
    INTEGER             ,INTENT(in)    :: nlev,ih2o
    INTEGER             ,INTENT(out)   :: ifail
    LOGICAL             ,INTENT(in)    :: debug
    REAL, DIMENSION(:)  ,INTENT(in)    :: press
    REAL, DIMENSION(:)  ,INTENT(inout) :: prof
    REAL, DIMENSION(:,:),INTENT(inout) :: Sx
    TYPE(StateIndex_t)  ,INTENT(in)    :: IG,NG,IR,NR
    REAL, DIMENSION(:,:),INTENT(in)    :: SxLiq,SxIce
    LOGICAL,OPTIONAL    ,INTENT(in)    :: invCldCov_in
    !---Local variables
    INTEGER              :: i
    REAL                 :: sv_pressure,cbaselim,cldtopcovMin=4.
    REAL                 :: qSat,qTol
    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: h2oBound
    INTEGER              :: ibit, qc
    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: self_temp, self_h2o
    real                 :: self_tskin, self_psfc
    logical              :: invCldCov
    real                 :: liqCovLoc, iceCovLoc

    !---Reading of the guesslimits parameters
    call loadGuessLimits(U_algconfig,F_algconfig)

    !---Set tightened covariance values according to whether or not the
    !    retrieval uses the inverse of the covariance matrix:
    invCldCov = .false.
    if (present(invCldCov_in)) then
       if ( invCldCov_in ) invCldCov = .true.
    endif
    liqCovLoc = merge( 1./liqCov, liqCov, invCldCov )
    iceCovLoc = merge( 1./iceCov, iceCov, invCldCov )

    if (.not. allocated(self_temp)) allocate(self_temp(NG%temp))
    if (.not. allocated(h2oBound)) allocate(h2oBound(NG%mol(iH2O)))
    if (.not. allocated(self_h2o)) allocate(self_h2o(NG%mol(iH2O)))

    ifail=0
    qc = 0

    ! check temperature: overwrite the default values from the namelist
    ibit = 0     !can be a different bit
    TempBnd%lower = tmin
    TempBnd%upper = tmax
    self_temp = prof(IG%temp:IG%temp+NG%temp-1)
    qc = checkBoundsVector(var=self_temp,bound=TempBnd,qc=qc,ibit=ibit,fill=.true.)

    ! check Tskin: overwrite the default values from the namelist
    TskinBnd%lower = tmin
    TskinBnd%upper = tmax
    self_tskin = prof(IG%tskin)
    qc = checkBoundsScalar(var=self_tskin,bound=TskinBnd,qc=qc,ibit=ibit,fill=.true.)

    ! check Psfc: use default values
    self_psfc = prof(IG%psfc)
    qc = checkBoundsScalar(var=self_psfc,bound=PsfcBnd,qc=qc,ibit=ibit,fill=.true.)

    if (btest(qc,ibit)) then
       ifail=1
       return
    endif

    ! calculate saturation vapor pressure
    H2oLoop: do i=0,NG%mol(iH2O)-1
       sv_pressure=svp(prof(IG%Temp+i))
       qSat=mrFromVp(sv_pressure,press(i+1))
       ! compute supersaturation tolerance; units of x(iH2O) are in g/g
       if (press(i+1) >= Pmin4SaturChk) then
          qTol=qSat*htol
          h2oBound(i+1)=minval((/qTol,H2Obnd%upper,prof(IG%mol(iH2o)+i)/))
       else
          qTol=qSat*absoluteLimitSupersaturation
          h2oBound(i+1) = min(qTol,H2Obnd%upper)
       endif
    enddo H2oLoop

    ! check mixing ratios for super-saturation
    self_h2o = prof(IG%mol(iH2O):IG%mol(iH2O)+NG%mol(iH2O)-1)
    qc = checkBoundsVector(var=self_h2o, &
         bound=H2OBnd,upperVector=h2oBound,qc=qc,ibit=ibit,cap=.true.)
    if (btest(qc,ibit)) & ! update h2o profile if super-saturated
         prof(IG%mol(iH2O):IG%mol(iH2O)+NG%mol(iH2O)-1) = self_h2o

    ! check cloud parameters if NR%cldLiq is nonzero
    IF (NR%Cldliq.gt.0 .AND. NR%Cldliq.lt.NLev-1) THEN
       if (prof(IG%Cldliq+1) .lt. thickmin) then
          prof(IG%Cldliq+1)=thickmin
       endif
       !---Check CLW:
       if (prof(IG%Cldliq+2) .lt. 0.) then
          prof(IG%Cldliq+2) = 0.
          Sx(IR%Cldliq+2,IR%Cldliq+2) = liqCovLoc
       else if (Sx(IR%Cldliq+2,IR%Cldliq+2) .eq. liqCovLoc) then
          Sx(IR%Cldliq+2,IR%Cldliq+2) = SxLiq(1+2,1+2)
       end if
       !---Check Deff:
       IF ( NR%Cldliq .GE. 4) THEN
          if (prof(IG%Cldliq+3) .lt. DeffLiqMin) then
             prof(IG%Cldliq+3) = DeffLiqMin
          end if
          if (prof(IG%Cldliq+3) .gt. DeffLiqMax) then
             prof(IG%Cldliq+3) = DeffLiqMax
          end if
       ENDIF
       !---Check cloud base:
       cbaselim=prof(IG%Psfc)
       if (prof(IG%Cldliq)+prof(IG%Cldliq+1) .gt. cbaselim) then
          if (debug) print *,'Liquid cloud base=', &
               prof(IG%Cldliq)+prof(IG%Cldliq+1),cbaselim
          if(Sx(IR%Cldliq+1,IR%Cldliq+1).le.cldtopcovMin)then   !cloud thickness not being retrieved
             prof(IG%CldLiq)=cbaselim-prof(IG%CldLiq+1)
          else
             if (prof(IG%CldLiq)+thickmin .gt. cbaselim) then
                if (debug) print *,'Mw cloud top low=', &
                     prof(IG%CldLiq)+thickmin,cbaselim
                prof(IG%CldLiq) = cbaselim-thickmin
             endif
             prof(IG%CldLiq+1) = cbaselim-prof(IG%CldLiq)
          endif
       end if
       !---Check CTP:
       if (prof(IG%CldLiq) .lt. ctopLow) then
          if (debug) print *,'Liquid cloud top high=',prof(IG%CldLiq),ctopLow
          prof(IG%CldLiq) = ctopLow
          Sx(IR%CldLiq,IR%CldLiq) = liqCovLoc
       else if (Sx(IR%CldLiq,IR%CldLiq) .eq. liqCovLoc) then
          Sx(IR%CldLiq,IR%CldLiq) = SxLiq(1,1)
       end if
    else
       !-------- placeholder for cloud profile retrieval checks
    endif
    ! check cloud parameters if NR%cldIce is nonzero
    IF (NR%CldIce.gt.0 .AND. NR%CldIce.lt.NLev-1) THEN
       if (prof(IG%CldIce+1) .lt. thickmin) then
          prof(IG%CldIce+1)=thickmin !if < thickmin, set thickness to thickmin
       endif
       !---Check IWP:
       if (prof(IG%CldIce+2) .lt. 0.) then
          prof(IG%CldIce+2) = 0. !if < 0, set IceAmt to 0.
          Sx(IR%CldIce+2,IR%CldIce+2) = iceCovLoc !and reset Cov
       else if (Sx(IR%CldIce+2,IR%CldIce+2) .eq. iceCovLoc) then
          Sx(IR%CldIce+2,IR%CldIce+2) = SxIce(1+2,1+2)
       end if
       !---Check Deff:
       if (prof(IG%CldIce+3) .lt. DeffIceMin) then
          prof(IG%CldIce+3) = DeffIceMin
       end if
       if (prof(IG%CldIce+3) .gt. DeffIceMax) then
          prof(IG%CldIce+3) = DeffIceMax
       end if
       !---Check cloud base:
       cbaselim=prof(IG%Psfc) !set cldBase at Surface
       if (prof(IG%CldIce)+prof(IG%CldIce+1) .gt. cbaselim) then
          !CldBottom below surface
          if (debug) print *,'Ice cloud base=', &
               prof(IG%CldIce)+prof(IG%CldIce+1),cbaselim
          if(Sx(IR%CldIce+1,IR%CldIce+1).le.cldtopcovMin)then   !cloud thickness not being retrieved
             prof(IG%CldIce)=cbaselim-prof(IG%CldIce+1) !reset cldTop
          else
             if (prof(IG%CldIce)+thickmin .gt. cbaselim) then
                if (debug) print *,'Mw cloud top low=', &
                     prof(IG%CldIce)+thickmin,cbaselim
                prof(IG%CldIce) = cbaselim-thickmin !reset cldTop with using thickmin
             endif
             prof(IG%CldIce+1) = cbaselim-prof(IG%CldIce)
          endif
       end if
       !---Check CTP:
       if (prof(IG%CldIce) .lt. icetopLow) then
          if (debug) print *,'Ice cloud top high=',prof(IG%CldIce),icetopLow
          prof(IG%CldIce) = icetopLow
          Sx(IR%CldIce,IR%CldIce) = iceCovLoc
       else if (Sx(IR%CldIce,IR%CldIce) .eq. iceCovLoc) then
          Sx(IR%CldIce,IR%CldIce) = SxIce(1,1)
       end if
    else
       !-------- placeholder for cloud profile retrieval checks
    endif

    return
  end subroutine chkges

  !--
  ! This subroutine is for re-writing covariance matrix. No clear vision
  ! exists yet about either it would be necessary to use
  !--
  subroutine resetCloudCovs(IR,NR,Sx,SxLiq,SxIce)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   resetCloudCovs
!
! PURPOSE:
!
!   Reset covariance matrix to a prior state.
!
! SYNTAX:
!
!   CALL resetCloudCovs(IR, NR, Sx, SxLiq, SxIce)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   IR     TYPE(STATEINDEX_T)  Starting indices for sections of retrieval
!                              state vector
!   NR     TYPE(STATEINDEX_T)  Number of retrieved elements for sections of
!                              retrieval state vector
!   SxLiq  REAL                Background error covariance for liquid cloud
!   SxIce  REAL                Background error covariance for ice cloud
!
!   INPUTS/OUTPUTS:
!
!   Sx     REAL                Background error covariance
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    TYPE(StateIndex_t),     INTENT(IN)    :: IR,NR
    REAL,   DIMENSION(:,:), INTENT(IN)    :: SxLiq,SxIce
    REAL,   DIMENSION(:,:), INTENT(INOUT) :: Sx

    Sx(IR%Cldliq:IR%Cldliq+NR%Cldliq-1,IR%Cldliq:IR%Cldliq+NR%Cldliq-1) = &
         SxLiq(1:NR%Cldliq,1:NR%Cldliq)
    Sx(IR%CldIce:IR%CldIce+NR%CldIce-1,IR%CldIce:IR%CldIce+NR%CldIce-1) = &
         SxIce(1:NR%CldIce,1:NR%CldIce)
    return
  end subroutine resetCloudCovs

  !--
  ! Set global variables from argument list
  !--
  subroutine setPropertiesGuessLimits(liqCov_in,iceCov_in, &
       maxcldbasebak_in,htol_in, &
       tmax_in,tmin_in,ctoplow_in,icetoplow_in, &
       deffLiqMin_in,deffIceMin_in,deffLiqMax_in,deffIceMax_in, &
       thickmin_in,Pmin4SaturChk_in)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   setPropertiesGuessLimits
!
! PURPOSE:
!
!   Set global variables from argument list, for physical limits.
!
! SYNTAX:
!
!   CALL setPropertiesGuessLimits(liqCov_in, iceCov_in,
!      maxcldbasebak_in, htol_in, tmax_in, tmin_in, ctoplow_in,
!      icetoplow_in, deffLiqMin_in, deffIceMin_in, deffLiqMax_in,
!      deffIceMax_in, thickmin_in, Pmin4SaturChk_in)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   liqCov_in         REAL  covariance to apply to clw when clw < 0
!   iceCov_in         REAL  covariance to apply to iwp when iwp < 0
!   maxcldbasebak_in  REAL  max pressure allowed for cloud base
!   htol_in           REAL  supersaturation tolerance
!   tmax_in           REAL  maximum allowed temperature
!   tmin_in           REAL  minimum allowed temperature
!   ctoplow_in        REAL  minimum allowed liquid cloud top pressure
!   icetoplow_in      REAL  minimum allowed ice cloud top pressure
!   deffLiqMin_in     REAL  minimum allowed liquid cloud particle size
!   deffIceMin_in     REAL  minimum allowed ice cloud particle size
!   deffLiqMax_in     REAL  maximum allowed liquid cloud particle size
!   deffIceMax_in     REAL  maximum allowed ice cloud particle size
!   thickmin_in       REAL  minimum allowed cloud thickness
!   Pmin4SaturChk_in  REAL  lowest pressure at which to check for
!                           supersaturation
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !--I/O variables
    real,  intent(in)           :: liqCov_in,iceCov_in
    real,  intent(in)           :: maxcldbasebak_in,htol_in,tmax_in,tmin_in
    real,  intent(in)           :: icetoplow_in,ctoplow_in
    real,  intent(in)           :: deffLiqMin_in,deffIceMin_in
    real,  intent(in)           :: deffLiqMax_in,deffIceMax_in
    real,  intent(in)           :: Pmin4SaturChk_in,thickmin_in

    !------------------------------------------------------------------------
    ! Set global variables from input arguments
    !------------------------------------------------------------------------
    liqCov        = liqCov_in
    iceCov        = iceCov_in
    maxcldbasebak = maxcldbasebak_in
    htol          = htol_in
    tmax          = tmax_in
    tmin          = tmin_in
    ctoplow       = ctoplow_in
    icetoplow     = icetoplow_in
    deffLiqMin    = deffLiqMin_in
    deffIceMin    = deffIceMin_in
    deffLiqMax    = deffLiqMax_in
    deffIceMax    = deffIceMax_in
    thickmin      = thickmin_in
    Pmin4SaturChk = Pmin4SaturChk_in

    loadGuessLimitsDone=.true.
    return
  end subroutine setPropertiesGuessLimits

  !--
  ! Read items from guessLimits namelist into local variables
  !--
  subroutine loadGuessLimits(unit_algcfg,file_algcfg)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   loadGuessLimits
!
! PURPOSE:
!
!   Load physical limits from a file.
!
! SYNTAX:
!
!   CALL loadGuessLimits(unit_algcfg, file_algcfg)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   unit_algcfg  INTEGER  file unit number for tuning data
!   file_algcfg  CHAR     Tuning data file
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !--I/O variables
    integer,                       intent(in)  :: unit_algcfg
    character(len=*),              intent(in)  :: file_algcfg

    !---Local variables
    real                 :: maxcldbasebak,htol,tmax,tmin
    real                 :: ctoplow,icetoplow,deffLiqMin,deffIceMin
    real                 :: deffLiqMax,deffIceMax
    real                 :: thickmin,Pmin4SaturChk,cldcov,icecov

    namelist /guesslimits/ cldcov,icecov,maxcldbasebak,htol,tmax,tmin, &
         ctoplow,icetoplow,deffLiqMin,deffIceMin,deffLiqMax,deffIceMax, &
         thickmin,Pmin4SaturChk

    if(loadGuessLimitsDone) return

    !------------------------------------------------------------------------
    ! Read namelist file items into local variables
    !------------------------------------------------------------------------
    open(unit_algcfg,file=file_algcfg)
    read(unit_algcfg,guesslimits)
    close(unit_algcfg)

    !------------------------------------------------------------------------
    ! Set global variables
    !------------------------------------------------------------------------
    call setPropertiesGuessLimits(cldcov,iceCov,maxcldbasebak,htol,tmax,tmin, &
         ctoplow,icetoplow,deffLiqMin,deffIceMin,deffLiqMax,deffIceMax, &
         thickmin,Pmin4SaturChk)

    return

  end subroutine loadGuessLimits

  !--
  SUBROUTINE chkRad(nchan,Ym,kchanIn,kchan,MinNChan,RadInvld,ValidRad,qc)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   chkRad
!
! PURPOSE:
!
!   Test validity of radiometric data, to set switches and QC flags.
!
! SYNTAX:
!
!   CALL chkRad(nchan, Ym, kchanIn, kchan, MinNChan, RadInvld, ValidRad, qc)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   nchan     INTEGER  Number of channels
!   Ym        REAL     Radiometric measurements
!   kchanIn   INTEGER  Channel on/off mask for input
!   MinNChan  INTEGER  lowest number of channels allowed
!   RadInvld  REAL     lowest radiometric measurement value allowed
!   qc        INTEGER  quality control index
!
!   INPUTS/OUTPUTS:
!
!   kchan     LOGICAL  Channel on/off mask
!   ValidRad  LOGICAL  Flag for valid/invalid radiometric measurements
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !---Input variables
    INTEGER,               INTENT(IN)    :: nchan,MinNChan
    INTEGER, DIMENSION(:), INTENT(IN)    :: qc
    LOGICAL, DIMENSION(:), INTENT(IN)    :: kchanIn
    REAL,    DIMENSION(:), INTENT(IN)    :: Ym
    REAL,                  INTENT(IN)    :: RadInvld
    !---Output variables
    LOGICAL, DIMENSION(:), INTENT(INOUT) :: kchan
    LOGICAL,               INTENT(INOUT) :: ValidRad
    !---Local variables
    INTEGER                              :: ichan,nTotChan,nQC

    !---By default kchan is initialized to kchanIn
    kchan(1:nchan)=kchanIn(1:nchan)
    !---Browse the radiance vector and set kchan to 0 if bad value
    DO ichan=1,nchan
       IF (Ym(ichan)<RadInvld) kchan(ichan) = .false.
    ENDDO
    !---Decide the validity of the radiances set
    ValidRad = .TRUE.
    nTotChan = COUNT(kchan(1:nchan))
    IF (nTotChan < MinNChan) ValidRad=.FALSE.
    nQC=size(qc)
    IF (nQC.NE.0) THEN
       IF (ANY(qc(1:nQC) .NE. 0)) ValidRad=.FALSE.
    ENDIF
    RETURN
  END SUBROUTINE chkRad

  !--
  ! Sets global variables from argument list--
  ! vectors that contain the ID, lower and upper bounds
  ! of physical parameters via dummy arguments.
  ! Each bounding variable should be of type "physicalBounds".
  ! Current parameters: temperature, skin temperature, pressure,
  ! and water vapor. This routine is public.
  !--
  subroutine setPropertiesBounds(bndT,bndTskin,bndP,bndH2O)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   setPropertiesBounds
!
! PURPOSE:
!
!   Sets global variables from argument list, for vectors that contain the ID,
!   lower and upper bounds of physical parameters.
!
! SYNTAX:
!
!   CALL setPropertiesBounds(bndT, bndTskin, bndP, bndH2O)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   bndT*     TYPE(PHYSICALBOUNDS)  structure with maximum and minimum bounds
!                                   for temperature
!   bndTskin* TYPE(PHYSICALBOUNDS)  structure with maximum and minimum bounds
!                                   for surface skin temperature
!   bndP*     TYPE(PHYSICALBOUNDS)  structure with maximum and minimum bounds
!                                   for surface pressure
!   bndH2O*   TYPE(PHYSICALBOUNDS)  structure with maximum and minimum bounds
!                                   for water vapor concentration
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !---Input variables
    type(physicalBounds), intent(in), optional :: bndT, bndTskin, bndP, bndH2O

    !------------------------------------------------------------------------
    ! Set global variables from input arguments
    !------------------------------------------------------------------------
    if (present(bndT) .and. trim(TempBnd%ID) == trim(bndT%ID)) then
       TempBnd = bndT
    endif
    if (present(bndTskin) .and. trim(TskinBnd%ID) == trim(bndTskin%ID)) then
       TskinBnd = bndTskin
    endif
    if (present(bndP) .and. trim(PsfcBnd%ID) == trim(bndP%ID)) then
       PsfcBnd = bndP
    endif
    if (present(bndH2O) .and. trim(H2OBnd%ID) == trim(bndH2O%ID)) then
       H2OBnd = bndH2O
    endif

    loadBoundsDone=.true.
    return
  end subroutine setPropertiesBounds

  !--
  ! Loads vectors that contain the ID, lower and upper bounds
  ! of physical parameters via namelist. Each field should
  ! contain a triplet, i.e., ID of type string, lower bound of type
  ! real, and upper bound of type real. This routine is public.
  !--
  subroutine loadBounds(U_bndconfig,F_bndconfig)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   loadBounds
!
! PURPOSE:
!
!   Loads vectors that contain the ID, lower and upper bounds of physical
!   parameters via namelist.
!
! SYNTAX:
!
!   CALL loadBounds(U_bndconfig, F_bndconfig)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   U_bndconfig  INTEGER  Unit number for physical bounds data
!   F_bndconfig  CHAR     File path for physical bounds data
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !---Input variables
    integer,              intent(in)           :: U_bndconfig
    character(len=*),     intent(in)           :: F_bndconfig

    !---Local variables
    type(physicalBounds) :: boundT, boundTskin, boundP, boundH2O

    namelist /geoBounds/boundT,boundTskin,boundP,boundH2O

    if(loadBoundsDone) return

    !------------------------------------------------------------------------
    ! Read namelist file items into local variables
    !------------------------------------------------------------------------
    open(U_bndconfig,file=F_bndconfig)
    READ(U_bndconfig,geoBounds)
    close(U_bndconfig)

    !------------------------------------------------------------------------
    ! Set global variables
    !------------------------------------------------------------------------
    call setPropertiesBounds(boundT,boundTskin,boundP,boundH2O)

    return
  end subroutine loadBounds

  !--
  ! Returns upper bound upon request with an input string.
  ! 'temp'  = temperature
  ! 'tskin' = skin temperature
  ! 'pre'   = pressure
  ! 'h2o'   = water vapor
  !--
  function queryUpperBound(varname)

!<f90Function>**********************************************************
!
! NAME:
!
!   queryUpperBound
!
! PURPOSE:
!
!   Returns upper bound upon request with an input string.
!
! SYNTAX:
!
!   Results=queryUpperBound(varname)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   varname          CHAR     Name of variable
!
!   * OPTIONAL
!
! RETURN:
!
!     INTEGER
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

    character (len=*), intent(in) :: varname
    integer                       :: queryUpperBound
    if (trim(varname) == trim(TempBnd%ID)) then
       queryUpperBound = TempBnd%upper
    else if (trim(varname) == trim(TskinBnd%ID)) then
       queryUpperBound = TskinBnd%upper
    else if (trim(varname) == trim(PsfcBnd%ID)) then
       queryUpperBound = PsfcBnd%upper
    else if (trim(varname) == trim(H2OBnd%ID)) then
       queryUpperBound = H2OBnd%upper
    else
       print*, 'err[ChkValid::queryUpperbound]: variable not in the list.'
    endif
  end function queryUpperBound

  !--
  ! Returns lower bound upon request with an input string.
  ! 'temp'  = temperature
  ! 'tskin' = skin temperature
  ! 'pre'   = pressure
  ! 'h2o'   = water vapor
  !--
  function queryLowerBound(varname)

!<f90Function>**********************************************************
!
! NAME:
!
!   queryLowerBound
!
! PURPOSE:
!
!   Returns lower bound upon request with an input string.
!
! SYNTAX:
!
!   Results=queryLowerBound(varname)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   varname          CHAR     Name of variable
!
!   * OPTIONAL
!
! RETURN:
!
!     INTEGER
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

    character (len=*), intent(in) :: varname
    integer                       :: queryLowerBound
    if (trim(varname) == trim(TempBnd%ID)) then
       queryLowerBound = TempBnd%lower
    else if (trim(varname) == trim(TskinBnd%ID)) then
       queryLowerBound = TskinBnd%lower
    else if (trim(varname) == trim(PsfcBnd%ID)) then
       queryLowerBound = PsfcBnd%lower
    else if (trim(varname) == trim(H2OBnd%ID)) then
       queryLowerBound = H2OBnd%lower
    else
       print*, 'err[ChkValid::queryLowerbound]: variable not in the list.'
    endif
  end function queryLowerBound

  !--
  ! Check all the physical bounds for an individual profile
  !--
  subroutine checkProfile(prof,IG,NG,molID,qc,ibit)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   checkProfile
!
! PURPOSE:
!
!   Check all the physical bounds for an individual profile
!
! SYNTAX:
!
!   CALL checkProfile(prof, IG, NG, molID, qc, ibit)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   IG     TYPE(STATEINDEX_T)  Starting indices for sections of geophysical
!                              state vector
!   NG     TYPE(STATEINDEX_T)  Number of elements for sections of geophysical
!                              state vector
!   molID  INTEGER             List of IDs of relevant molecular species
!   ibit*  INTEGER             index of bit on which to operate
!
!   INPUTS/OUTPUTS:
!
!   prof   REAL                Atmospheric profile
!   qc*    INTEGER             quality control index
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    real,    dimension(:), intent(inout)           :: prof
    type(StateIndex_t),    intent(in)              :: IG, NG
    integer, dimension(:), intent(in)              :: molID
    integer,               intent(inout), optional :: qc
    integer,               intent(in),    optional :: ibit
    real,    dimension(:), allocatable :: self_temp, self_h2o
    real                 :: self_tskin, self_psfc
    integer :: ibit_local=0, qc_local=0, i, iH2O

    if (present(qc)) qc_local = qc
    if (present(ibit)) ibit_local=ibit  !ibit_local may be different for each parameter

    if (.not. allocated(self_temp)) allocate(self_temp(NG%temp))
    self_temp = prof(IG%temp:IG%temp+NG%temp-1)
    qc_local = checkBoundsVector(self_temp,&
         bound=TempBnd,qc=qc_local,ibit=ibit_local,fill=.true.)

    self_tskin = prof(IG%tskin)
    qc_local = checkBoundsScalar(self_tskin,&
         bound=TskinBnd,qc=qc_local,ibit=ibit_local,fill=.true.)

    self_psfc = prof(IG%psfc)
    qc_local = checkBoundsScalar(self_psfc,&
         bound=PsfcBnd,qc=qc_local,ibit=ibit_local,fill=.true.)

    iH2O = whereH2O(molID)
    if (.not. allocated(self_h2o)) allocate(self_h2o(NG%mol(iH2O)))
    self_h2o = prof(IG%mol(iH2O):IG%mol(iH2O)+NG%mol(iH2O)-1)
    qc_local = checkBoundsVector(self_h2o, &
         bound=H2OBnd,qc=qc_local,ibit=ibit_local,fill=.true.)

    if (present(qc)) qc = qc_local

    if (allocated(self_temp)) deallocate(self_temp)
    if (allocated(self_h2o)) deallocate(self_h2o)

  end subroutine checkProfile

  !--
  function checkBoundsScalar(var,bound,qc,ibit,fill,cap,debug)

!<f90Function>**********************************************************
!
! NAME:
!
!   checkBoundsScalar
!
! PURPOSE:
!
!   Scalar interface for checkBoundsVector.
!
! SYNTAX:
!
!   Results=checkBoundsScalar(var, bound, qc, ibit, fill, cap, debug)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   var                REAL                  Generic variable
!   bound*             TYPE(PHYSICALBOUNDS)  structure with maximum and
!                                            minimum bounds
!   ibit               INTEGER               index of bit on which to operate
!   fill*              LOGICAL               Flag to fill out-of-bounds
!                                            values with missing data flag
!   cap*               LOGICAL               Flag reset values to the bounds
!                                            if they are out of bounds
!   debug*             LOGICAL               Flag for debugging mode
!
!   INPUTS/OUTPUTS:
!
!   qc*                INTEGER               quality control index
!
!   * OPTIONAL
!
! RETURN:
!
!     INTEGER
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

    real,                 intent(in)              :: var
    type(physicalBounds), intent(in),    optional :: bound
    integer,              intent(inout), optional :: qc
    integer,              intent(in)              :: ibit
    logical,              intent(in),    optional :: fill
    logical,              intent(in),    optional :: cap
    logical,              intent(in),    optional :: debug
    integer               :: checkBoundsScalar
    real, dimension(1)    :: var_local
    var_local = var
    checkBoundsScalar = checkBoundsVector(var_local,bound=bound,&
         qc=qc,ibit=ibit,fill=fill,cap=cap,debug=debug)
  end function checkBoundsScalar

  !--
  ! Check the bounds for an arbitrary physical parameter, and
  ! set the bit to 1 if out-of-bound is detected.
  !--
  function checkBoundsVector(var,bound,lowerVector,upperVector,qc, &
       ibit,fill,cap,debug)

!<f90Function>**********************************************************
!
! NAME:
!
!   checkBoundsVector
!
! PURPOSE:
!
!   Check the bounds for an arbitrary vector of physical parameters, and signal
!   if out-of-bound is detected.
!
! SYNTAX:
!
!   Results=checkBoundsVector(var, bound, lowerVector, upperVector,
!      qc, ibit, fill, cap, debug)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   bound              TYPE(PHYSICALBOUNDS)  structure with maximum and
!                                            minimum bounds
!   lowerVector*       REAL                  lower bounds vector
!   upperVector*       REAL                  upper bounds vector
!   qc*                INTEGER               quality control index
!   ibit*              INTEGER               index of bit on which to operate
!   fill*              LOGICAL               Flag to fill out-of-bounds
!                                            values with missing data flag
!   cap*               LOGICAL               Flag reset values to the bounds
!                                            if they are out of bounds
!   debug*             LOGICAL               Flag for debugging mode
!
!   INPUTS/OUTPUTS:
!
!   var                REAL                  Generic variable
!
!   * OPTIONAL
!
! RETURN:
!
!     INTEGER
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

    real,   dimension(:), intent(inout)         :: var
    type(physicalBounds), intent(in)            :: bound
    real,   dimension(:), intent(in),  optional :: lowerVector
    real,   dimension(:), intent(in),  optional :: upperVector
    integer,              intent(in),  optional :: qc
    integer,              intent(in),  optional :: ibit
    logical,              intent(in),  optional :: fill
    logical,              intent(in),  optional :: cap
    logical,              intent(in),  optional :: debug
    integer               :: checkBoundsVector
    integer :: ibit_local
    logical :: debug_local
    type(physicalBounds), dimension(size(var)) :: fill_local
    type(physicalBounds), dimension(size(var)) :: bnd_local
    integer :: i
    integer, dimension(size(var)) :: qc_local

    !initializes of QC
    if (present(qc)) then
       checkBoundsVector = qc
    else
       checkBoundsVector = 0
    endif

    !initializes bit position for evaluating QC
    if (present(ibit)) then
       ibit_local = ibit
    else
       ibit_local = 0
    endif

    !initializes debug flag
    if (present(debug)) then
       debug_local = debug
    else
       debug_local = .false.
    endif

    !initializes bound(s)
    bnd_local = bound
    if (present(lowerVector)) bnd_local%lower = lowerVector
    if (present(upperVector)) bnd_local%upper = upperVector

    !initializes fill values
    fill_local%lower=var
    fill_local%upper=var
    if (present(fill)) then
       if (fill) then
          fill_local%lower=MISSING_REAL
          fill_local%upper=MISSING_REAL
       endif
    endif
    if (present(cap)) then
       if (cap) then
          fill_local%lower=bnd_local%lower
          fill_local%upper=bnd_local%upper
       endif
    endif

    !checks bound(s) of all elements
    do i = 1, size(var)
       if (var(i) < bnd_local(i)%lower) then
          if (debug_local) &
               print*,"[warning - ChkValid::checkBoundsVector]: ", &
               bnd_local(i)%ID," of value ",var(i), " out of bound at level ", &
               i, " ...lower bound: ", bnd_local(i)%lower
          var(i) = fill_local(i)%lower
          qc_local(i) = 1
       else if (var(i) > bnd_local(i)%upper) then
          if (debug_local) &
               print*,"[warning - ChkValid::checkBoundsVector]: ", &
               bnd_local(i)%ID," of value ",var(i), " out of bound at level ", &
               i, " ...upper bound: ", bnd_local(i)%upper
          var(i) = fill_local(i)%upper
          qc_local(i) = 1
       else
          qc_local(i) = 0
       endif
    enddo
    if (any(qc_local /= 0)) then
       checkBoundsVector = ibset(checkBoundsVector,ibit_local)
    endif
  end function checkBoundsVector

  !--
  SUBROUTINE retrQC(validRad,qcRetr)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   retrQC
!
! PURPOSE:
!
!   Set an integer QC flag based on status of logical QC flag
!
! SYNTAX:
!
!   CALL retrQC(validRad, qcRetr)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   validRad  LOGICAL  Flag for valid/invalid radiometric measurements
!
!   INPUTS/OUTPUTS:
!
!   qcRetr    INTEGER  quality control index
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    LOGICAL,              INTENT(IN)    :: validRad
    INTEGER, DIMENSION(:),INTENT(INOUT) :: qcRetr
    qcRetr =0
    IF (.NOT.validRad) qcRetr(1)=1
  END SUBROUTINE retrQC

  !--
  SUBROUTINE getChiSq(Ym,y,rerr,chiSq,rms,kchan,nchanin,nch,chanDataCompressed)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   getChiSq
!
! PURPOSE:
!
!   Compute normalized chi-square and RMS difference values.
!
! SYNTAX:
!
!   CALL getChiSq(Ym, y, rerr, chiSq, rms, kchan, nchanin, nch)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   Ym       REAL     Radiometric measurements
!   y        REAL     Radiometric data computed from state vector
!   rerr     REAL     Instrument noise standard deviation
!   kchan    LOGICAL  Channel on/off mask
!            computed radiance is compressed
!            size(y) == sum(kchan>0)
!   nchanin  INTEGER  Number of channels for input
!   chanDataCompressed     logical      = .false. y is uncompressed
!                                          .true. y is compressed
!                         for every y(1:nch)>0,
!                         define idx(1:nch) that kchan(idx)>0
!                         then compression means that y(j) corresponds to Ym(idx(j))
!
!
!   INPUTS/OUTPUTS:
!
!   chiSq    REAL     chi-square metric
!   rms      REAL     RMS difference
!   nch      INTEGER  Number of channels turned on
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !----Input/Output variables
    REAL,    DIMENSION(:), INTENT(IN)    :: y,Ym,rerr
    REAL,                  INTENT(INOUT) :: chiSq,rms
    LOGICAL, DIMENSION(:), INTENT(IN)    :: kchan
    INTEGER,               INTENT(IN)    :: nchanin
    INTEGER,               INTENT(INOUT) :: nch
    logical, OPTIONAL,     INTENT(IN)    :: chanDataCompressed
    !----Local variables
    INTEGER               :: i
    REAL                  :: delY2
    logical               :: chanDataCompressedLoc

    chiSq = 0.0
    rms = 0.0
    nch = 0

    if (present(chanDataCompressed)) then
      chanDataCompressedLoc = chanDataCompressed
    else
      chanDataCompressedLoc = .false.
    end if

    if (chanDataCompressedLoc) then
      DO i=1,nchanIn
         IF(kchan(i)) THEN
            nch=nch+1
            delY2=(Ym(i)-Y(nch))**2
            chiSq=chiSq+delY2/rerr(i)**2
            rms=rms+dely2
         END IF
      END DO
    else
      DO i=1,nchanIn
         IF(kchan(i)) THEN
            nch=nch+1
            delY2=(Ym(i)-Y(i))**2
            chiSq=chiSq+delY2/rerr(i)**2
            rms=rms+dely2
         END IF
      END DO
    end if

    rms=SQRT(rms/float(nch))
    chiSq=(chiSq/float(nch))
    RETURN
  END SUBROUTINE getChiSq

END MODULE CHKVALID
