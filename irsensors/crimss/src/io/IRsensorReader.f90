!<f90File>**************************************************************
!
! CONTACT:
!
!   Atmospheric & Environmental Research, Inc
!   131 Hartwell Ave
!   Lexington ,MA 02421-3126 USA
!   Phone: 781.761.2288
!   E-mail: ipolonsk@aer.com
!
! COPYRIGHT NOTICE:
!
!   Copyright AER, Inc 2001-2009, All Rights Reserved
!   See the file README-DATARIGHTS.txt included with this release
!   for additional details.

MODULE IRsensorReader
  !
  ! Module that loads the CrIS L1B netCDF file
  !
  ! Subroutines
  !
  !     readIRsensorFile
  !
  ! USE:
  !     netcdf
  !
  ! ipolonsk@aer.com, 03/19/2018
  !
implicit none
public :: readIRsensorFile
private
  type :: crisGranule_t
    integer                              :: atrack
    integer                              :: xtrack
    integer                              :: nFOV
    integer                              :: nChan_lw
    integer                              :: nChan_mw
    integer                              :: nChan_sw
    real*4,  dimension(:,:,:,:), allocatable :: rad_lw           !(nChan_lw,nFOV,xtrack,atrack)
    real*4,  dimension(:,:,:,:), allocatable :: rad_mw           !(nChan_mw,nFOV,xtrack,atrack)
    real*4,  dimension(:,:,:,:), allocatable :: rad_sw           !(nChan_sw,nFOV,xtrack,atrack)
    real*4,  dimension(:,:),     allocatable :: nedn_lw          !(nChan_lw,nFOV)
    real*4,  dimension(:,:),     allocatable :: nedn_mw          !(nChan_mw,nFOV)
    real*4,  dimension(:,:),     allocatable :: nedn_sw          !(nChan_sw,nFOV)
    integer*8, dimension(:,:,:), allocatable :: l1b_qual         !(nFOV,xtrack,atrack)
    real*8,  dimension(:),       allocatable :: wnum_lw          !(nChan_lw)
    real*8,  dimension(:),       allocatable :: wnum_mw          !(nChan_mw)
    real*8,  dimension(:),       allocatable :: wnum_sw          !(nChan_sw)
    integer*8,dimension(:,:,:),  allocatable :: geo_qual         !(nFOV,xtrack,atrack)
    real*8,  dimension(:,:),     allocatable :: obs_time_tai     !(xtrack,atrack)

    real*4,  dimension(:,:,:),   allocatable :: lat              !(nFOV,xtrack,atrack)
    real*4,  dimension(:,:,:),   allocatable :: lon              !(nFOV,xtrack,atrack)
    real*4,  dimension(:,:,:),   allocatable :: sol_zen          !(nFOV,xtrack,atrack)
    real*4,  dimension(:,:,:),   allocatable :: sol_azi          !(nFOV,xtrack,atrack)
    real*4,  dimension(:,:,:),   allocatable :: sat_zen          !(nFOV,xtrack,atrack)
    real*4,  dimension(:,:,:),   allocatable :: sat_azi          !(nFOV,xtrack,atrack)
    real*4,  dimension(:,:,:),   allocatable :: sat_range        !(nFOV,xtrack,atrack)

    real*4,  dimension(:,:,:),   allocatable :: land_frac        !(nFOV,xtrack,atrack)
    real*4,  dimension(:,:,:),   allocatable :: surf_alt         !(nFOV,xtrack,atrack)
    real*4,  dimension(:,:,:),   allocatable :: surf_alt_sdev    !(nFOV,xtrack,atrack)

    real*4,  dimension(:),       allocatable :: sun_glint_lat    !(atrack)
    real*4,  dimension(:),       allocatable :: sun_glint_lon    !(atrack)
    integer, dimension(6)                    :: startTime
    integer, dimension(6)                    :: endTime
    integer, dimension(:,:,:),   allocatable :: obs_time_utc     !(8,xtrack,atrack)

  end type crisGranule_t

! l1b CRIS file include 2 extra spectral points at boundaries
logical,          parameter :: reduceGrid = .true.
integer,          parameter :: nPointRed = 2

character(len=*), parameter :: modName ='IRsensorReader::'
character(len=*), parameter :: modErr = 'ERR: ' // modName
type(crisGranule_t)         :: crisGran
logical         , parameter :: dbg = .false.
real, dimension(5), parameter  :: blackmanWeight = (/0.04, 0.25, 0.42, 0.25, 0.04/)
real, dimension(3), parameter  :: hammingWeight  = (/0.23, 0.54, 0.23/)

contains
  subroutine getCrisDims(nFOR, nFOV, nChan)
    integer, optional, intent(out) :: nFOR
    integer, optional, intent(out) :: nFOV
    integer, optional, intent(out) :: nChan

    if (present(nFOR)) nFOR = crisGran%atrack * abs(crisGran%xtrack)
    if (present(nFOV)) nFOV = crisGran%nFOV
    if (present(nChan))nChan= crisGran%nChan_lw+crisGran%nChan_mw +crisGran%nChan_sw
  end subroutine getCrisDims

  subroutine destroyCrisFile()
    call deAllocateCrisGran(crisGran)
  end subroutine destroyCrisFile

  function readIRsensorFile(fileName, IRgran, independentFOV,apodizationType)
    use IRobsStructure, only : IRgran_t
    USE ControlStructure, Only: GenControl_t

    character(len=*), intent(in)  :: fileName
    logical                       :: readIRsensorFile
    type(IRgran_t), intent(inout)  :: IRgran
    logical, optional, intent(in) :: independentFOV
    integer, optional,intent(in)   :: apodizationType  ! 0-sinc, 1-blackman, 2- hamming

    logical  :: independentFOVloc

    if (present(independentFOV)) then
      independentFOVloc=independentFOV
    else
      independentFOVloc=.false.
    end if

    readIRsensorFile =  readCrisFile(fileName,apodizationType)
    if (readIRsensorFile) &
         call convertCrisGranToIRgran(IRgran, independentFOVloc)

  end function readIRsensorFile

! verify that the file is SIPS netCDF file
  function readCrisFile(fname,apodizationType)
    use ncdfUtil, only: check, getDimLen, readVarNC, inqDim
    use netcdf, only: nf90_open, nf90_close, nf90_nowrite, NF90_GLOBAL, nf90_get_att
    character(len=*), intent(in)  :: fname
    integer, optional,intent(in)  :: apodizationType  ! 0-sinc, 1-blackman, 2- hamming
    logical                       :: readCrisFile
    integer     :: ncid
    integer     :: nAtrack, nXtrack, nChan_lw, nChan_mw, nChan_sw, nFOV
    integer     :: kk
    real, dimension(:), allocatable  :: bufReal
    character(len=100)               :: charBuf

    if (dbg) print *, 'readCrisFile: ', trim(fname)
    call check( nf90_open(fname, nf90_nowrite, ncid) )
    ! the CrIS l1b file has 'atrack' and 'xtrack' dimension which
    ! our simulation file does not so the presence of 'atrack' dimension was
    ! selected to differentiate  between the files
    if (.not. inqDim(ncid, "atrack",  dimLen=nAtrack)) then
      call check( nf90_close(ncid) )
      readCrisFile = .false.
      return
    end if
    if (dbg) print *, 'reading time_coverage_start'

    call check(nf90_get_att(ncid, NF90_GLOBAL, 'time_coverage_start', charBuf))
    read(charBuf,'(I4,5(1x,I2))') crisGran%startTime

    call check(nf90_get_att(ncid, NF90_GLOBAL, 'time_coverage_end', charBuf))
    read(charBuf,'(I4,5(1x,I2))') crisGran%endTime

    call getDimLen(ncid, "xtrack",  nXtrack)
    call getDimLen(ncid, "fov",     nFOV)
    call getDimLen(ncid, "chan_lw", nChan_lw)
    call getDimLen(ncid, "chan_mw", nChan_mw)
    call getDimLen(ncid, "chan_sw", nChan_sw)

    call allocateCrisGran(nAtrack, nXtrack, nChan_lw, nChan_mw, nChan_sw, nFOV, crisGran)

    call readVarNC(ncid,"rad_lw", crisGran%rad_lw)
    call readVarNC(ncid,"rad_mw", crisGran%rad_mw)
    call readVarNC(ncid,"rad_sw", crisGran%rad_sw)

    call readVarNC(ncid,"nedn_lw", crisGran%nedn_lw)
    call readVarNC(ncid,"nedn_mw", crisGran%nedn_mw)
    call readVarNC(ncid,"nedn_sw", crisGran%nedn_sw)

    call readVarNC(ncid,"wnum_lw", crisGran%wnum_lw)
    call readVarNC(ncid,"wnum_mw", crisGran%wnum_mw)
    call readVarNC(ncid,"wnum_sw", crisGran%wnum_sw)

    print *,'apod ',apodizationType
    if (apodizationType==1) then ! apply blackman Apodization
      call blackmanApodization(crisGran%rad_lw)
      call blackmanApodization(crisGran%rad_mw)
      call blackmanApodization(crisGran%rad_sw)
    else if(apodizationType==2) then ! apply hamming Apodization
      call hammingApodization(crisGran%rad_lw)
      call hammingApodization(crisGran%rad_mw)
      call hammingApodization(crisGran%rad_sw)
    end if

! In  the  V1.0  product,  the  Noise  Equivalent  Differential  Radiances  (NEdN;  nedn_lw,
! nedn_mw, and nedn_sw) are inadvertently provided in raw instrument counts units
! (inconsistent with the naming convention and attributes).
! To convert these values to their approximate radiance equivalents, multiply by a Planck radiance for a blackbody at
! a temperature of 280K. This issue will be corrected in the next releases
! correction
    allocate(bufReal(nChan_lw))
    call Planck(crisGran%wnum_lw, 280.0, bufReal)
    do kk=1, nFOV
      crisGran%nedn_lw(:,kk) = bufReal*crisGran%nedn_lw(:,kk)
    end do
    deallocate(bufReal)
    allocate(bufReal(nChan_mw))
    call Planck(crisGran%wnum_mw, 280.0, bufReal)
    do kk=1, nFOV
      crisGran%nedn_mw(:,kk) = bufReal*crisGran%nedn_mw(:,kk)
    end do
    deallocate(bufReal)

    allocate(bufReal(nChan_sw))
    call Planck(crisGran%wnum_sw, 280.0, bufReal)
    do kk=1, nFOV
      crisGran%nedn_sw(:,kk) = bufReal*crisGran%nedn_sw(:,kk)
    end do
    deallocate(bufReal)
!----------------------------------------------
    call readVarNC(ncid,"l1b_qual", crisGran%l1b_qual)
    call readVarNC(ncid,"geo_qual", crisGran%geo_qual)
    call readVarNC(ncid,"obs_time_tai", crisGran%obs_time_tai)

    call readVarNC(ncid,"lat", crisGran%lat)
    call readVarNC(ncid,"lon", crisGran%lon)
    call readVarNC(ncid,"land_frac", crisGran%land_frac)


    call readVarNC(ncid,"sol_zen", crisGran%sol_zen)
    call readVarNC(ncid,"sol_azi", crisGran%sol_azi)
    call readVarNC(ncid,"sat_zen", crisGran%sat_zen)
    call readVarNC(ncid,"sat_azi", crisGran%sat_azi)

    call readVarNC(ncid,"sat_range", crisGran%sat_range)
    call readVarNC(ncid,"surf_alt", crisGran%surf_alt)
    call readVarNC(ncid,"surf_alt_sdev", crisGran%surf_alt_sdev)
    call readVarNC(ncid,"obs_time_utc", crisGran%obs_time_utc)

    call check( nf90_close(ncid) )

    readCrisFile=.true.
    if (dbg) print *, 'readCrisFile: done'

    return
  end function readCrisFile

   ! allocation
   subroutine allocateCrisGran(nAtrack, nXtrack, nChan_lw, nChan_mw, nChan_sw, nFOV, this)
      integer,   intent(in) :: nAtrack, nXtrack, nChan_lw, nChan_mw, nChan_sw, nFOV
      type(crisGranule_t), intent(inout) :: this

      call deAllocateCrisGran(this)

      allocate(this%rad_lw(nChan_lw,nFOV,nXtrack,nAtrack))
      allocate(this%rad_mw(nChan_mw,nFOV,nXtrack,nAtrack))
      allocate(this%rad_sw(nChan_sw,nFOV,nXtrack,nAtrack))
      allocate(this%nedn_lw(nChan_lw,nFOV))
      allocate(this%nedn_mw(nChan_mw,nFOV))
      allocate(this%nedn_sw(nChan_sw,nFOV))

      allocate(this%wnum_lw(nChan_lw))
      allocate(this%wnum_mw(nChan_mw))
      allocate(this%wnum_sw(nChan_sw))

      allocate(this%l1b_qual(nFOV,nXtrack,nAtrack))
      allocate(this%geo_qual(nFOV,nXtrack,nAtrack))
      allocate(this%obs_time_tai(nXtrack,nAtrack))

      allocate(this%lat(nFOV,nXtrack,nAtrack))
      allocate(this%lon(nFOV,nXtrack,nAtrack))

      allocate(this%sol_zen(nFOV,nXtrack,nAtrack))
      allocate(this%sol_azi(nFOV,nXtrack,nAtrack))
      allocate(this%sat_zen(nFOV,nXtrack,nAtrack))
      allocate(this%sat_azi(nFOV,nXtrack,nAtrack))

      allocate(this%sat_range(nFOV,nXtrack,nAtrack))
      allocate(this%land_frac(nFOV,nXtrack,nAtrack))

      allocate(this%surf_alt(nFOV,nXtrack,nAtrack))
      allocate(this%surf_alt_sdev(nFOV,nXtrack,nAtrack))
      allocate(this%obs_time_utc(8,nXtrack,nAtrack))
      this%atrack = nAtrack
      this%xtrack = nXtrack
      this%nFOV   = nFOV

      this%nChan_lw = nChan_lw
      this%nChan_mw = nChan_mw
      this%nChan_sw = nChan_sw
   end subroutine allocateCrisGran

   subroutine deAllocateCrisGran(this)
      type(crisGranule_t), intent(inout) :: this
      if (allocated(this%rad_lw))         deallocate(this%rad_lw)
      if (allocated(this%rad_mw))         deallocate(this%rad_mw)
      if (allocated(this%rad_sw))         deallocate(this%rad_sw)

      if (allocated(this%nedn_lw))        deallocate(this%nedn_lw)
      if (allocated(this%nedn_mw))        deallocate(this%nedn_mw)
      if (allocated(this%nedn_sw))        deallocate(this%nedn_sw)
      if (allocated(this%wnum_lw))        deallocate(this%wnum_lw)
      if (allocated(this%wnum_mw))        deallocate(this%wnum_mw)
      if (allocated(this%wnum_sw))        deallocate(this%wnum_sw)

      if (allocated(this%sol_zen))        deallocate(this%sol_zen)
      if (allocated(this%sol_azi))        deallocate(this%sol_azi)
      if (allocated(this%sat_zen))        deallocate(this%sat_zen)
      if (allocated(this%sat_azi))        deallocate(this%sat_azi)

      if (allocated(this%l1b_qual))       deallocate(this%l1b_qual)
      if (allocated(this%obs_time_tai))   deallocate(this%obs_time_tai)
      if (allocated(this%lat))            deallocate(this%lat)
      if (allocated(this%lon))            deallocate(this%lon)
      if (allocated(this%sat_range))      deallocate(this%sat_range)
      if (allocated(this%geo_qual))       deallocate(this%geo_qual)
      if (allocated(this%land_frac))      deallocate(this%land_frac)

      if (allocated(this%surf_alt))       deallocate(this%surf_alt)
      if (allocated(this%surf_alt_sdev))  deallocate(this%surf_alt_sdev)
      if (allocated(this%obs_time_utc))   deallocate(this%obs_time_utc)

      this%atrack = -1
      this%xtrack = -1
      this%nFOV   = -1
      this%nChan_lw = -1
      this%nChan_mw = -1
      this%nChan_sw = -1
   end subroutine deAllocateCrisGran

   subroutine convertCrisGranToIRgran( IRgran, independentFOV)
    use constants, ONLY : MISSING_REAL
    use IRobsStructure, only : IRgran_t
    type(IRgran_t), intent(inout)  :: IRgran
    logical                        :: independentFOV
    !
    character(len=*), parameter :: procName = &
          modName // 'convertCrisGranToIRgran: '
    integer :: ifor, ifov, ich, ix, ia, jj, idx
    integer :: slw, smw, ssw
    integer :: elw, emw, esw
    !
    IRgran%nFOR      = crisGran%atrack * abs(crisGran%xtrack)
    if (independentFOV) then
      IRgran%nFOR      = IRgran%nFOR * crisGran%nFOV
      IRgran%def%nFOV  = 1
    else
      IRgran%def%nFOV  = crisGran%nFOV
    end if

    if (reduceGrid) then
      slw = 1
      elw = crisGran%nChan_lw - 2*nPointRed
      smw = elw + 1
      emw = elw + crisGran%nChan_mw - 2*nPointRed
      ssw = emw + 1
      esw = emw + crisGran%nChan_sw - 2*nPointRed
    else
      slw = 1
      elw = crisGran%nChan_lw
      smw = elw + 1
      emw = elw + crisGran%nChan_mw
      ssw = emw + 1
      esw = emw + crisGran%nChan_sw
    end if
    IRgran%def%nChan = esw

    if (.not. allocated(IRgran%def%wvn))     allocate(IRgran%def%wvn(IRgran%def%nChan))
    if (.not. allocated(IRgran%def%chanNum)) allocate(IRgran%def%chanNum(IRgran%def%nChan))
    if (.not. allocated(IRgran%def%chanID))  allocate(IRgran%def%chanID(IRgran%def%nChan))
    if (.not. allocated(IRgran%FOR))         allocate(IRgran%FOR(IRgran%nFOR))
    if (.not. allocated(IRgran%def%NEdN))    allocate(IRgran%def%NEdN(IRgran%def%nChan,IRgran%def%nFOV))


    do ifor = 1, IRgran%nFOR
       if (.not. allocated(IRgran%FOR(ifor)%ob)) &
            allocate(IRgran%FOR(ifor)%ob(IRgran%def%nFOV))
       do ifov = 1, IRgran%def%nFOV
          if (.not. allocated(IRgran%FOR(ifor)%ob(ifov)%rad)) &
               allocate(IRgran%FOR(ifor)%ob(ifov)%rad(IRgran%def%nChan))
          if (.not. allocated(IRgran%FOR(ifor)%ob(ifov)%QC)) &
               allocate(IRgran%FOR(ifor)%ob(ifov)%QC(IRgran%def%nChan))
       end do
    end do

    IRgran%def%startTime(:) = crisGran%startTime
    IRgran%def%endTime(:)   = crisGran%endTime

    forall (ich=1:IRgran%def%nChan) IRgran%def%chanNum(ich) = ich

    IRgran%def%wvn(slw:elw) = &
             crisGran%wnum_lw(1+nPointRed:crisGran%nChan_lw-nPointRed)
    IRgran%def%wvn(smw:emw) = &
             crisGran%wnum_mw(1+nPointRed:crisGran%nChan_mw-nPointRed)
    IRgran%def%wvn(ssw:esw) = &
             crisGran%wnum_sw(1+nPointRed:crisGran%nChan_sw-nPointRed)

    if (independentFOV) then
      IRgran%def%NEdN(:,1) = 0.
      do ifov = 1, crisGran%nFOV
        IRgran%def%NEdN(slw:elw,1)=IRgran%def%NEdN(slw:elw,1) + &
                                   crisGran%nedn_lw(1+nPointRed:crisGran%nChan_lw-nPointRed,ifov)

        IRgran%def%NEdN(smw:emw,1)=IRgran%def%NEdN(smw:emw,1) + &
                                   crisGran%nedn_mw(1+nPointRed:crisGran%nChan_mw-nPointRed,ifov)

        IRgran%def%NEdN(ssw:esw,1)=IRgran%def%NEdN(ssw:esw,1) + &
                                   crisGran%nedn_sw(1+nPointRed:crisGran%nChan_sw-nPointRed,ifov)
      end do
      IRgran%def%NEdN(:,1) = IRgran%def%NEdN(:,1)/real(crisGran%nFOV)
    else
      do ifov = 1, IRgran%def%nFOV
        IRgran%def%NEdN(slw:elw,ifov)=crisGran%nedn_lw(1+nPointRed:crisGran%nChan_lw-nPointRed,ifov)
        IRgran%def%NEdN(smw:emw,ifov)=crisGran%nedn_mw(1+nPointRed:crisGran%nChan_mw-nPointRed,ifov)
        IRgran%def%NEdN(ssw:esw,ifov)=crisGran%nedn_sw(1+nPointRed:crisGran%nChan_sw-nPointRed,ifov)
      end do
    end if
    do ia = 1, crisGran%atrack
      do ix = 1, crisGran%xtrack
        do idx = 1, crisGran%nFOV

          if (independentFOV) then
             ifor = idx + ((ix - 1) + (ia-1)*crisGran%xtrack)*crisGran%nFOV
             ifov = 1
          else
             ifor = ix  + (ia-1)*crisGran%xtrack
             ifov = idx
          end if

          IRgran%FOR(ifor)%ob(ifov)%rad(slw:elw) = &
                   crisGran%rad_lw(1+nPointRed:crisGran%nChan_lw-nPointRed,idx,ix,ia)
          IRgran%FOR(ifor)%ob(ifov)%rad(smw:emw) = &
                   crisGran%rad_mw(1+nPointRed:crisGran%nChan_mw-nPointRed,idx,ix,ia)
          IRgran%FOR(ifor)%ob(ifov)%rad(ssw:esw) = &
                   crisGran%rad_sw(1+nPointRed:crisGran%nChan_sw-nPointRed,idx,ix,ia)

          IRgran%FOR(ifor)%ob(ifov)%surf_alt   = crisGran%surf_alt(idx,ix,ia)
          IRgran%FOR(ifor)%ob(ifov)%land_frac  = crisGran%land_frac(idx,ix,ia)

          IRgran%FOR(ifor)%ob(ifov)%lat   = crisGran%lat(idx,ix,ia)
          IRgran%FOR(ifor)%ob(ifov)%lon   = crisGran%lon(idx,ix,ia)
          IRgran%FOR(ifor)%ob(ifov)%eia   = crisGran%sat_zen(idx,ix,ia)
          IRgran%FOR(ifor)%ob(ifov)%EAA   = crisGran%sat_azi(idx,ix,ia)
          IRgran%FOR(ifor)%ob(ifov)%SolIA = crisGran%sol_zen(idx,ix,ia)
          IRgran%FOR(ifor)%ob(ifov)%SolAA = crisGran%sol_azi(idx,ix,ia)
          IRgran%FOR(ifor)%ob(ifov)%QC    = crisGran%geo_qual(idx,ix,ia) + &
                                            crisGran%l1b_qual(idx,ix,ia)
          IRgran%FOR(ifor)%obs_time_utc   = crisGran%obs_time_utc(:,ix,ia)
        end do
        IRgran%FOR(ifor)%ob(:)%scanAng = MISSING_REAL
      end do
    end do
    if (independentFOV) then
       do  ifor = 1, IRgran%nFOR
         IRgran%FOR(ifor)%lat    = IRgran%FOR(ifor)%ob(1)%lat
         IRgran%FOR(ifor)%lon    = IRgran%FOR(ifor)%ob(1)%lon
         IRgran%FOR(ifor)%EIA    = IRgran%FOR(ifor)%ob(1)%EIA
         IRgran%FOR(ifor)%EAA    = IRgran%FOR(ifor)%ob(1)%EAA
         IRgran%FOR(ifor)%SolIA  = IRgran%FOR(ifor)%ob(1)%SolIA
         IRgran%FOR(ifor)%SolAA  = IRgran%FOR(ifor)%ob(1)%SolAA
         IRgran%FOR(ifor)%scanAng= IRgran%FOR(ifor)%ob(1)%scanAng
       enddo
    else
       do  ifor = 1, IRgran%nFOR
         IRgran%FOR(ifor)%lat    = IRgran%FOR(ifor)%ob(5)%lat
         IRgran%FOR(ifor)%lon    = IRgran%FOR(ifor)%ob(5)%lon
         IRgran%FOR(ifor)%EIA    = IRgran%FOR(ifor)%ob(5)%EIA
         IRgran%FOR(ifor)%EAA    = IRgran%FOR(ifor)%ob(5)%EAA
         IRgran%FOR(ifor)%SolIA  = IRgran%FOR(ifor)%ob(5)%SolIA
         IRgran%FOR(ifor)%SolAA  = IRgran%FOR(ifor)%ob(5)%SolAA
         IRgran%FOR(ifor)%scanAng= IRgran%FOR(ifor)%ob(5)%scanAng
       enddo
    end if
    if (dbg) print *, 'IRgran set'
    call deAllocateCrisGran(crisGran)

  end subroutine convertCrisGranToIRgran

  subroutine Planck(f, t, rad)
    use constants, ONLY : LtSp_cgs, Bltz_cgs, Plnk_cgs
      !  f: cm^{-1}
      !  t: Kelvin
      !rad: mW/m^2/str/cm^{-1}

      REAL*8,PARAMETER :: c1=2.*Plnk_cgs*LtSp_cgs**2
      REAL*8,PARAMETER :: c2=(Plnk_cgs*LtSp_cgs/Bltz_cgs)

      REAL*8,dimension(:), INTENT(IN)   :: f
      REAL               , INTENT(IN)   :: t
      REAL,  dimension(:), INTENT(INOUT):: rad

      rad = exp(-c2*f/t)
      rad = rad/(1. - rad)
      rad = c1*f**3*rad
      RETURN
  END subroutine Planck

  subroutine blackmanApodization(rad)
    real, dimension(:,:,:,:), intent(inout) :: rad
    real, dimension(size(rad,dim=1))        :: buf
    integer jj, sz, j2, j3, j4

    sz=size(rad, dim=1)-2

    do j4=1,size(rad, dim=4)
      do j3=1,size(rad, dim=3)
        do j2=1,size(rad, dim=2)
          buf=rad(:,j2,j3,j4)
          forall (jj=3:sz) rad(jj,j2,j3,j4)=sum(blackmanWeight*buf(jj-2:jj+2))
        enddo
      enddo
    enddo
  end subroutine blackmanApodization

  subroutine hammingApodization(rad)
    real, dimension(:,:,:,:), intent(inout) :: rad
    real, dimension(size(rad,dim=1))        :: buf
    integer jj, sz, j2, j3, j4

    sz=size(rad, dim=1)-1

    do j4=1,size(rad, dim=4)
      do j3=1,size(rad, dim=3)
        do j2=1,size(rad, dim=2)
          buf=rad(:,j2,j3,j4)
          forall (jj=2:sz) rad(jj,j2,j3,j4)=sum(hammingWeight*buf(jj-1:jj+1))
        enddo
      enddo
    enddo
  end subroutine hammingApodization

end MODULE IRsensorReader
