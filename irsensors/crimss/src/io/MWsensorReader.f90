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
!
!*************************************************************</f90File>
MODULE MWsensorReader
  !
  ! Module that loads the ATMS data and ancillary data
  ! Subroutines
  !
  !     readMWsensorFile
  !
  ! USE:
  !
  !     netcdf_module
  !     MWobsStructure
  implicit none
  public :: readMWsensorFile
  private
  type :: atmsGranule_t
    integer                              :: atrack
    integer                              :: xtrack
    integer                              :: nChan

    character(len=1),dimension(:),      allocatable :: polarization     !(channel)
    real*4,           dimension(:),     allocatable :: center_freq      !(channel)
    real*4,           dimension(:,:,:), allocatable :: antenna_temp     !(channel,xtrack,atrack)
    integer*4,        dimension(:,:),   allocatable :: calib_degraded   !(channel,atrack)
    real*8,           dimension(:,:),   allocatable :: obs_time_tai     !(channel,atrack)
    real*4,           dimension(:),     allocatable :: cold_nedt        !(channel)
    real*4,           dimension(:),     allocatable :: warm_nedt        !(channel)
    real*4,           dimension(:,:),   allocatable :: lat              !(xtrack,atrack)
    real*4,           dimension(:,:),   allocatable :: lon              !(xtrack,atrack)
    real*4,           dimension(:,:),   allocatable :: sol_zen          !(xtrack,atrack)
    real*4,           dimension(:,:),   allocatable :: sol_azi          !(xtrack,atrack)
    real*4,           dimension(:,:),   allocatable :: sat_zen          !(xtrack,atrack)
    real*4,           dimension(:,:),   allocatable :: sat_azi          !(xtrack,atrack)
    real*4,           dimension(:,:),   allocatable :: sat_range        !(xtrack,atrack)
    integer*2,        dimension(:,:),   allocatable :: geo_qual         !(xtrack,atrack)
    real*4,           dimension(:,:),   allocatable :: land_frac        !(xtrack,atrack)
    real*4,           dimension(:,:),   allocatable :: surf_alt         !(xtrack,atrack)
    real*4,           dimension(:,:),   allocatable :: surf_alt_sdev    !(xtrack,atrack)

    real*4,           dimension(:,:),   allocatable :: sun_glint_lat    !(xtrack,atrack)
    real*4,           dimension(:,:),   allocatable :: sun_glint_lon    !(xtrack,atrack)
    real*4,           dimension(:,:),   allocatable :: sun_glint_dist   !(xtrack,atrack)
    integer,          dimension(6)                  :: startTime
    integer,          dimension(6)                  :: endTime
    integer,          dimension(:,:,:), allocatable :: obs_time_utc     !(8,xtrack,atrack)
  end type atmsGranule_t

  character(len=*), parameter :: modName = 'atmsModule::'
  character(len=*), parameter :: modErr  = 'ERR: '//modName
  type(atmsGranule_t), save   :: atmsGran
  logical, parameter          :: dbg = .false.
contains

  function readMWsensorFile(fileName, mwGran)
    use MWobsStructure, only : MWgran_t
    character(len=*), intent(in)  :: fileName
    logical                       :: readMWsensorFile
    type(MWgran_t), intent(inout)  :: mwGran

    readMWsensorFile =  readAtmsFile(fileName)
    if (readMWsensorFile) &
         call convertAtmsGranToMWgran(mwGran)

  end function readMWsensorFile

  subroutine destroyAtmsFile()
    call deAllocateAtmsGran(atmsGran)
  end subroutine destroyAtmsFile

! verify that the file is SIPS netCDF file
  function readAtmsFile(fname)
    use ncdfUtil, only: check, getDimLen, readVarNC, inqDim, readStringAttr
    use netcdf, only: nf90_open, nf90_close, nf90_nowrite, NF90_GLOBAL, nf90_get_att

    character(len=*), intent(in)  :: fname
    logical                       :: readAtmsFile
    integer     :: ncid
    integer     :: nAtrack
    integer     :: nXtrack
    integer     :: nChan, kk
    character(len=100)               :: charBuf
    if (dbg) print *, 'reading: ', trim(fname)
    call check( nf90_open(fname, nf90_nowrite, ncid) )
    ! the ATMS l1b file has 'atrack' and 'xtrack' dimension which
    ! our simulation file does not so the presence of 'atrack' dimension was
    ! selected to differentiate  between the files
    if (.not. inqDim(ncid, "atrack",  dimLen=nAtrack)) then
      call check( nf90_close(ncid) )
      readAtmsFile = .false.
      return
    end if

    call readStringAttr(ncid, 'time_coverage_start', charBuf)
    read(charBuf,'(I4,5(1x,I2))') atmsGran%startTime

    call readStringAttr(ncid, 'time_coverage_end', charBuf)
    ! call check(nf90_get_att(ncid, NF90_GLOBAL, 'time_coverage_end', charBuf))
    read(charBuf,'(I4,5(1x,I2))') atmsGran%endTime

    call getDimLen(ncid, "xtrack",  nXtrack)
    call getDimLen(ncid, "channel", nChan)

    call allocateAtmsGran(nAtrack, nXtrack, nChan, atmsGran)

    call readVarNC(ncid,"center_freq", atmsGran%center_freq)
    call readVarNC(ncid,"antenna_temp", atmsGran%antenna_temp)
    call readVarNC(ncid,"cold_nedt", atmsGran%cold_nedt)
    call readVarNC(ncid,"warm_nedt", atmsGran%warm_nedt)
    call readVarNC(ncid,"lat", atmsGran%lat)
    call readVarNC(ncid,"lon", atmsGran%lon)

    call readVarNC(ncid,"land_frac", atmsGran%land_frac)

    call readVarNC(ncid,"calib_degraded", atmsGran%calib_degraded)
    call readVarNC(ncid,"obs_time_tai", atmsGran%obs_time_tai)
    call readVarNC(ncid,"geo_qual", atmsGran%geo_qual, fatal=.false.)

    call readVarNC(ncid,"sol_zen", atmsGran%sol_zen)
    call readVarNC(ncid,"sol_azi", atmsGran%sol_azi)
    call readVarNC(ncid,"sat_zen", atmsGran%sat_zen)
    call readVarNC(ncid,"sat_azi", atmsGran%sat_azi)

    call readVarNC(ncid,"sat_range", atmsGran%sat_range)
    call readVarNC(ncid,"surf_alt", atmsGran%surf_alt)
    call readVarNC(ncid,"surf_alt_sdev", atmsGran%surf_alt_sdev)

    call readVarNC(ncid,"polarization", atmsGran%polarization)
    call readVarNC(ncid,"obs_time_utc", atmsGran%obs_time_utc)

    call check( nf90_close(ncid) )
    readAtmsFile=.true.
    return
  end function readAtmsFile

   ! allocation
   subroutine allocateAtmsGran(nAtrack, nXtrack, nChan, this)
      integer,   intent(in) :: nAtrack, nXtrack, nChan
      type(atmsGranule_t), intent(inout) :: this

      call deAllocateAtmsGran(this)

      allocate(this%center_freq(nChan))
      allocate(this%polarization(nChan))

      allocate(this%antenna_temp(nChan,nXtrack,nAtrack))
      allocate(this%calib_degraded(nChan,nAtrack))
      allocate(this%obs_time_tai(nXtrack,nAtrack))

      allocate(this%cold_nedt(nChan))
      allocate(this%warm_nedt(nChan))

      allocate(this%lat(nXtrack,nAtrack))
      allocate(this%lon(nXtrack,nAtrack))

      allocate(this%sol_zen(nXtrack,nAtrack))
      allocate(this%sol_azi(nXtrack,nAtrack))
      allocate(this%sat_zen(nXtrack,nAtrack))
      allocate(this%sat_azi(nXtrack,nAtrack))

      allocate(this%sat_range(nXtrack,nAtrack))
      allocate(this%geo_qual(nXtrack,nAtrack))
      allocate(this%land_frac(nXtrack,nAtrack))

      allocate(this%surf_alt(nXtrack,nAtrack))
      allocate(this%surf_alt_sdev(nXtrack,nAtrack))
      allocate(this%sun_glint_lat(nXtrack,nAtrack))
      allocate(this%sun_glint_lon(nXtrack,nAtrack))
      allocate(this%sun_glint_dist(nXtrack,nAtrack))
      allocate(this%obs_time_utc(8,nXtrack,nAtrack))
      this%atrack = nAtrack
      this%xtrack = nXtrack
      this%nChan  = nChan
   end subroutine allocateAtmsGran

   subroutine deAllocateAtmsGran(this)
      type(atmsGranule_t), intent(inout) :: this

      if (allocated(this%center_freq))    deallocate(this%center_freq)
      if (allocated(this%polarization))   deallocate(this%polarization)

      if (allocated(this%antenna_temp))   deallocate(this%antenna_temp)
      if (allocated(this%calib_degraded)) deallocate(this%calib_degraded)
      if (allocated(this%obs_time_tai))   deallocate(this%obs_time_tai)

      if (allocated(this%cold_nedt))      deallocate(this%cold_nedt)
      if (allocated(this%warm_nedt))      deallocate(this%warm_nedt)

      if (allocated(this%lat))            deallocate(this%lat)
      if (allocated(this%lon))            deallocate(this%lon)

      if (allocated(this%sol_zen))        deallocate(this%sol_zen)
      if (allocated(this%sol_azi))        deallocate(this%sol_azi)
      if (allocated(this%sat_zen))        deallocate(this%sat_zen)
      if (allocated(this%sat_azi))        deallocate(this%sat_azi)

      if (allocated(this%sat_range))      deallocate(this%sat_range)
      if (allocated(this%geo_qual))       deallocate(this%geo_qual)
      if (allocated(this%land_frac))      deallocate(this%land_frac)

      if (allocated(this%surf_alt))       deallocate(this%surf_alt)
      if (allocated(this%surf_alt_sdev))  deallocate(this%surf_alt_sdev)
      if (allocated(this%sun_glint_lat))  deallocate(this%sun_glint_lat)
      if (allocated(this%sun_glint_lon))  deallocate(this%sun_glint_lon)
      if (allocated(this%sun_glint_dist)) deallocate(this%sun_glint_dist)
      if (allocated(this%obs_time_utc))   deallocate(this%obs_time_utc)
      this%atrack = -1
      this%xtrack = -1
      this%nChan  = -1
   end subroutine deAllocateAtmsGran

   subroutine convertAtmsGranToMWgran( mwGran)
    use constants, ONLY : MISSING_REAL
    use MWobsStructure, only : MWgran_t
    TYPE(MWgran_t), intent(inout) :: mwGran
    character(len=*), parameter :: procName = &
                      modName // 'convertAtmsGranToMWgran: '
    integer :: ifov, ich, ix, ia

    MWgran%def%nChan = atmsGran%nChan
    if (.not. allocated(MWgran%def%pol)) allocate(MWgran%def%pol(MWgran%def%nChan))
    if (.not. allocated(MWgran%def%chanID)) allocate(MWgran%def%chanID(MWgran%def%nChan))
    if (.not. allocated(MWgran%def%planckAlpha)) allocate(MWgran%def%planckAlpha(MWgran%def%nChan))
    if (.not. allocated(MWgran%def%planckBeta)) allocate(MWgran%def%planckBeta(MWgran%def%nChan))
    if (.not. allocated(MWgran%def%frq)) allocate(MWgran%def%frq(MWgran%def%nChan))
    ! convert from MHz to GHz
    MWgran%def%frq = atmsGran%center_freq/1e3
    where (atmsGran%polarization == 'V')
      MWgran%def%pol = 0
    elsewhere
      MWgran%def%pol = 1
    end where

    MWgran%def%startTime(:) = atmsGran%startTime
    MWgran%def%endTime(:)   = atmsGran%endTime

    MWgran%nFOV = atmsGran%atrack * atmsGran%xtrack

    if (dbg) print *, 'MWgran%nFOV:', MWgran%nFOV, atmsGran%atrack,  atmsGran%xtrack

    allocate(MWgran%ob(MWgran%nFOV))
    do ifov = 1, MWgran%nFOV
       if (.not. allocated(MWgran%ob(ifov)%rad)) allocate(MWgran%ob(ifov)%rad(MWgran%def%nChan))
       if (.not. allocated(MWgran%ob(ifov)%QC))  allocate(MWgran%ob(ifov)%QC(MWgran%def%nChan))
    end do

    MWgran%def%radOrTb = 0 !  radiance in BT
    do ia = 1, atmsGran%atrack
      do ix = 1, atmsGran%xtrack
        ifov = ix + (ia-1)*atmsGran%xtrack
        MWgran%ob(ifov)%lat = atmsGran%lat(ix,ia)
        MWgran%ob(ifov)%lon = atmsGran%lon(ix,ia)

        MWgran%ob(ifov)%surf_alt = atmsGran%surf_alt(ix,ia)
        MWgran%ob(ifov)%land_frac = atmsGran%land_frac(ix,ia)
        MWgran%ob(ifov)%EIA = atmsGran%sat_zen(ix,ia)
        MWgran%ob(ifov)%EAA = atmsGran%sat_azi(ix,ia)
        MWgran%ob(ifov)%rad = atmsGran%antenna_temp(:,ix,ia)
        MWgran%ob(ifov)%QC  = atmsGran%geo_qual(ix,ia) + atmsGran%calib_degraded(:,ia)
        MWgran%ob(ifov)%scanAng = 0.
        MWgran%ob(ifov)%scanPos = ix
        MWgran%ob(ifov)%obs_time_utc = atmsGran%obs_time_utc(:,ix,ia)
      end do
    end do
    call deAllocateAtmsGran(atmsGran)
    return
   end subroutine convertAtmsGranToMWgran
end MODULE MWsensorReader
