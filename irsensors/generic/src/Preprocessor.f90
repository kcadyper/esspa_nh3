program Preprocessor
!*********************************************************************
!     IRpreproc
!     1. calculate surface pressure at FOV centers. Uses NWP
!        surface pressure field and adjusts it for local topography.
!     2. create MWrad.nc from MW sensor (ATMS)
!     3. if irSensorFile == 'NONE' then AUX file based on MW sensor is
!        created
!
!     Syntax: preprocessor < run.in
!*********************************************************************

  use MWobsStructure, only : MWgran_t
  use MWsensorReader, only : readMWsensorFile
  use IRobsStructure, only : IRgran_t
  use IRsensorReader, only : readIRsensorFile
  use NWPprof       , only : nwpProfGet
  use ncdf_module   , only : nf_noerr, nf_open, nf_close
  use netcdf        , only : nf90_nowrite, nf90_open
  use ncdfUtil      , only : check, getDimLen, readVarNC

  use StateIndexModule,only: StateIndex_t, initLengths, getVectorLength,&
                             genIndices, nullMolId, idH2O, maxMol
  use rad_io_module , only : openRad, putRad
  use scene_io_module,only : openScene, queryScene, putScene, closeScene
  use date_module   , only : RelativeDay, Calendar

  use constants     , only : MISSING_REAL
  use MetFunctions  , only : adjSfcPresAlt
  implicit none
 !-- ModePsfc = 1 in the following argument list causes psfcFov to
 !-- be calculated using NWP data. The subroutine can also be run
 !-- with ModePsfc=2, in which case psfcFov is calculated starting
 !-- with 1013 mbar.
  integer, parameter                       :: ModePsfc = 1
  integer, parameter                       :: nAng = 1
  integer, parameter                       :: nGrids=4

  type(MWgran_t)                           :: mwGran
  type(IRgran_t)                           :: irGran
  type(StateIndex_t)                       :: IGgrid, NGgrid, IGfov, NGfov
  integer                                  :: nlev
  integer                                  :: nlpres
  integer                                  :: status,ncid
  integer                                  :: nfiles, nfil,ndayfov,iFor, jFovMV
  integer                                  :: NParGgrid, NParGfov
  real                                     :: hourfov,elevFov,plandFov,psfcFov,fracTime
  integer                                  :: ncidScene
  integer                                  :: nwpDt
  character(len=10)                        :: date
  character(len=3)                         :: day
  character (len=8)                        :: xid
  character(len=256), dimension(:), allocatable :: filenames
  integer                                  :: nday_start,nday_end,nd,n,nday_startValid
  integer                                  :: nhour_start,nhour_end
  real                                     :: latFov,lonFov
  integer                                  :: nXtrack, nChan
  integer                                  :: fovIdx
  real,       dimension(nAng)              :: eia, EAA
  integer,   dimension(6)                  :: time,time_start,time_end,time_act
  integer,   dimension(maxMol)             :: MolID,MolIDFov
  REAL,      dimension(nGrids)             :: latGrid4, lonGrid4, plandGrid4, elevGrid4
  real,      dimension(nGrids)             :: wt
  real,      dimension(:),     allocatable :: presGrid ! NWP grid from top to 1000 mbar
  real,      dimension(:),     allocatable :: presOrigGrid   ! original NWP grid (1000 mbar - top)
  real*4,    dimension(:),     allocatable :: radMW
  real,      dimension(:),     allocatable :: tim
  integer,   dimension(:),     allocatable :: ncidNWP
  REAL,      dimension(:,:),   allocatable :: profGrid4
  REAL,      dimension(:),     allocatable :: profGrid
  REAL,      dimension(:),     allocatable :: profFov
  real*4,    dimension(:,:),   allocatable :: biasCorr
  integer,   dimension(:),     allocatable :: qcNWP
  integer                                  :: npts
  integer                                  :: nQcNWP=0
  integer                                  :: qcScene
  integer                                  :: ncMWRad
  character (len=256)                      :: irSensorFile
  character (len=256)                      :: mwSensorFile
  character (len=256)                      :: mwBiasCorrFile
  character (len=256)                      :: nwpDir,  nwpPath
  character (len=256)                      :: AuxFile, MWradFile
  logical                                  :: errorFound, mwOnly
  logical                                  :: independentFOV = .TRUE.
  real                                     :: NWPres = 1.0
  logical                                  :: toRead

  namelist /preprocInput/ irSensorFile, mwSensorFile, nwpDir, AuxFile, &
                          MWradFile, mwBiasCorrFile, independentFOV, nwpDt, NWPres
  read(*,preprocInput)

  print *, 'irSensorFile:  ', trim(irSensorFile)
  print *, 'mwSensorFile:  ', trim(mwSensorFile)
  print *, 'nwpDir:        ', trim(nwpDir)
  print *, 'AuxFile:       ', trim(AuxFile)
  print *, 'MWradFile:     ', trim(MWradFile)
  print *, 'mwBiasCorrFile:', trim(mwBiasCorrFile)
  print *, 'independentFOV:', independentFOV
  print *, 'NWPres:        ', nwpRes
  print *, 'nwpDt:         ', nwpDt

  call check( nf90_open(mwBiasCorrFile, nf90_nowrite, ncid) )
  call getDimLen(ncid, "xtrack",  nXtrack)
  call getDimLen(ncid, "channel", nChan)
  allocate(biasCorr(nChan, nXtrack), radMW(nChan))
  call readVarNC(ncid,"biasCorr", biasCorr)
  call check(nf_close(ncid))

  molIDFov = nullMolId
  NGfov    = initLengths(Npsfc=1)
  IGfov    = genIndices(NGfov)
  NParGfov = getVectorLength(NGfov)

  allocate(profFov(NParGfov))
  profFov=0.
  nlev   =0
  if (trim(irSensorFile) == 'NONE') then
    mwOnly = .true.
  else
    mwOnly = .false.
  end if

  print *, 'reading MW file', trim(mwSensorFile)
  if (.not. readMWsensorFile(mwSensorFile, mwGran)) then
     print *, 'ERROR reading MW file: ', trim(mwSensorFile)
     call exit(1)
  end if

  if (.not. mwOnly) then
    print *, 'reading IR file', trim(irSensorFile)
    if (.not. readIRsensorFile(irSensorFile, irGran, independentFOV)) then
       print *, 'ERROR reading IR file: ', trim(irSensorFile)
       call exit(1)
    end if
    npts = size(irGran%FOR)
    ! Get the valid starting and ending date
    time_end   = IRgran%def%endTime
    time_start = IRgran%def%startTime
    fovIdx = 5
  else
    npts = size(mwGran%ob)
    time_end   = mwGran%def%endTime
    time_start = mwGran%def%startTime
    fovIdx = 1
  end if
  if (independentFOV) fovIdx = 1

  nday_start = RelativeDay(time_start(3),time_start(2),time_start(1))
  nday_end   = RelativeDay(time_end(3),time_end(2),time_end(1))
  ! start
  nhour_start = time_start(4)/nwpDt
  nhour_end   = ((nday_end - nday_start) * 24 +  time_end(4))/nwpDt + 1
  nfiles = nhour_end - nhour_start + 1
  allocate(filenames(nfiles),ncidNWP(nfiles),tim(nfiles))

  ! check files
  errorFound = .false.
  nd =  nday_start
  time_act(4) = nhour_start*nwpDt
  tim(1)=nd*24. + time_act(4)
  toRead = .true.
  do n = 1, nfiles
    call Calendar(nd,time_act(3),time_act(2),time_act(1),day)
    write(date,'(i4,2("/",i2.2))') time_act(1:3)
    nwpPath = trim(nwpDir) // '/'//trim(date) // '/gfsanl_4_'
    write(date,'(i4,3i2.2)') time_act(1:4)
    filenames(n)=trim(nwpPath)//date//'.nc'

    ! Examine the validity of the file

    status=nf_open(filenames(n),0,ncid)
    if (status /= nf_noerr) then
      print *, 'ERR:[Preprocessor] : netCDF file error - ', trim(filenames(n))
      call check(status, fatal=.FALSE.)
      errorFound = .true.
      cycle
    else
      call check(nf_close(ncid))
    endif
    print *, n, ' NWP file:', trim(filenames(n))

    if(toRead) then
      toRead = .false.
      nday_startValid=nd
      call queryScene(file=filenames(n),Nlevel=nlpres,MolID=MolID, &
           IG=IGgrid,NG=NGgrid,nqc=nQcNWP)
      if (nlpres /= NGgrid%temp) then
         print *, "ERR:[Preprocessor]:  inconsistent pressure levels in ",&
              trim(filenames(n))
         call errorHalt(1)
      endif
      if (allocated(presGrid)) deallocate(presGrid)
      allocate(presGrid(nlpres))
      if (allocated(presOrigGrid)) deallocate(presOrigGrid)
      allocate(presOrigGrid(nlpres))
      NParGgrid = getVectorLength(NGgrid)
      if (allocated(profGrid4)) deallocate(profGrid4)
      allocate(profGrid4(NParGgrid,nGrids))
      if (allocated(profGrid)) deallocate(profGrid)
      allocate(profGrid(NParGgrid))
      if (nQcNWP > 0) then
        if (allocated(qcNWP)) deallocate(qcNWP)
        allocate(qcNWP(nQcNWP))
      endif
      call queryScene(file=filenames(n),pressure=presGrid,&
                      MolID=MolID,IG=IGgrid,NG=NGgrid)
      ! Invert indexing for use in adjSfcPresAlt
      presOrigGrid(1:nlpres)=presGrid(nlpres:1:-1)
    endif
    call openScene(ncid=ncidNWP(n),file=filenames(n),&
         MolID=MolID,IG=IGgrid,NG=NGgrid)
    ! next file
    time_act(4) = time_act(4) + nwpDt
    if ( time_act(4) >= 24 ) then
      time_act(4) = time_act(4) - 24
      nd = nd + 1
    end if
    if (n<nfiles) tim(n+1)=tim(n) + nwpDt
  enddo

  molID=nullMolId
  molID(1)=idH2O

  if(errorFound) then
      print *,'ERR:[Preprocessor]: NWP files problem - quitting'
      call errorHalt(1)
  endif

  call openScene(ncidScene,Auxfile,Nlevel=nlev,&
                  MolID=MolIDFov,IG=IGfov,NG=NGfov,status='new')

  CALL openRad(ncid=ncMWRad,fname=MWradFile,pol=MWgran%def%pol(1:MWgran%def%nChan), &
             frq=MWgran%def%frq(1:MWgran%def%nChan),nchan=MWgran%def%nChan, status='new')

  PROFloop: do iFor = 1,npts
     if (mwOnly) then
       time = mwGran%ob(iFor)%obs_time_utc(1:6)
       time(6) = time(6)*1000 + mwGran%ob(iFor)%obs_time_utc(7) ! time(6): msec
     else
       time = irGran%FOR(iFor)%obs_time_utc(1:6)
       time(6) = time(6)*1000 + irGran%FOR(iFor)%obs_time_utc(7) ! time(6): msec
     end if

       qcScene = 0
       if (nQcNWP > 0) qcNWP = 0

       if (mwOnly) then
         latFov   = mwGran%ob(iFor)%lat
         lonFov   = mwGran%ob(iFor)%lon
         ElevFov  = mwGran%ob(iFor)%surf_alt
         PlandFov = mwGran%ob(iFor)%land_frac
       else
         latFov   = irGran%FOR(iFor)%ob(fovIdx)%lat
         lonFov   = irGran%FOR(iFor)%ob(fovIdx)%lon
         ElevFov  = irGran%FOR(iFor)%ob(fovIdx)%surf_alt
         PlandFov = irGran%FOR(iFor)%ob(fovIdx)%land_frac
       endif

       if (ElevFov  > 1e35) ElevFov  = MISSING_REAL
       if (PlandFov > 1e35) PlandFov = MISSING_REAL

       if (all(time(1:3) > 0)) then
         ndayfov=RelativeDay(time(3),time(2),time(1))
         hourfov=ndayfov*24. + time(4) + time(5)/60. + time(6)/3.6e6 ! time(6): msec

         do n=2, nfiles-1
          if (tim(n) > hourfov) exit
         end do
         nfil = n - 1
         fracTime = (hourfov-tim(nfil))/nwpDt
         call nwpProfGet(latFov,lonFov,ncidNWP,nfil,fracTime,profGrid4, &
             latGrid4,lonGrid4,plandGrid4,elevGrid4,wt,molID,NGgrid,IGgrid, &
             qcNWP,NWPres)

         ! if NWP interpolations are invalid
         if (btest(qcNWP(1),0)) then
            qcScene = ibset(qcScene,0)
         else
            qcScene = ibclr(qcScene,0)
         endif

         call adjSfcPresAlt(latFov,lonFov,elevFov,profGrid4,latGrid4,lonGrid4, &
              elevGrid4,presGrid,IGgrid,NGgrid,wt,molID,ModePsfc,psfcFov)

      else
         psfcFov = MISSING_REAL
      endif

      write(xid, '("IR",I6.6)') iFor

      profFov(IGfov%psfc)=psfcFov

      call putScene(ncidScene,xid,latFov,lonFov,pland=plandFov, &
                    time=time,x=profFov,surfalt=ElevFov)

  ! find MW closest sounding
      call getNearest(MWgran,latFov, lonFov, jFovMV)
      write(xid, '("MW",I6.6)') jFovMV
      eia=MWgran%ob(jFovMV)%EIA
      EAA=MWgran%ob(jFovMV)%EAA
  ! bias correction
      radMW = MWgran%ob(jFovMV)%rad - biasCorr(:, MWgran%ob(jFovMV)%scanPos)
      CALL putRad(ncid=ncMWRad, &
         xid=xId, &
         lat=MWgran%ob(jFovMV)%lat,&
         lon=MWgran%ob(jFovMV)%lon,&
         time=Time,&
         eia=eia,&
         EAA=EAA,&
         QC=MWgran%ob(jFovMV)%QC,&
         tb=radMW)

  enddo PROFloop
  print *, 'processed nFOR: ', npts

  do n=1,nfiles
   call closeScene(ncid=ncidNWP(n))
  enddo

  call closeScene(ncid=ncidScene)
  call closeScene(ncid=ncMWRad)

  deallocate(presOrigGrid,profGrid4,profGrid,filenames,&
             tim,ncidNWP,biasCorr,radMW)
  if (allocated(qcNWP))    deallocate(qcNWP)

contains
  subroutine getNearest(MWgranLoc, cLat, cLon, idx)
    type(MWgran_t)             :: MWgranLoc

    real,                    intent(in) :: cLat
    real,                    intent(in) :: cLon
    integer,                intent(out) :: idx

    real                                :: minR, r
    integer                             :: ja

    minR = 1e32
    do ja = 1, size(MWgranLoc%ob)
        r = abs(MWgranLoc%ob(ja)%lat-cLat) + abs(MWgranLoc%ob(ja)%lon-cLon)
        if (r < minR) then
          minR = r
          idx = ja
        end if
    end do
  end  subroutine getNearest

end program Preprocessor
