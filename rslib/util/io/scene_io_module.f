c $"$Id$"
      module scene_io_module
c**********************************************************************
c* Module name: scene_io_module
c* Purpose: Scene I/O
c* Description: The module consists of the following functions,
c*              openScene()
c*              closeScene()
c*              getScene()
c*              putScene()
c*
c* Var_Name  Type       Description
c* --------  ----       -----------
c* ncid      integer    file id
c* file      string     file name
c* nprf      integer    number of profiles
c* polarity  real       channel polarities
c* freq      real       channel frequencies (ghz)
c* pressure  real       pressure grid
c* casename  character  case identification string
c* nAngle    integer    number earth azimuth or incidence angles
c*                      (becoming obsolete, use nAng in the future)
c* nAng      integer    number earth azimuth or incidence angles
c* nLevel    integer    number of pressure levels
c* nchmw     integer    number of MW channels
c* nSfcGrid  integer    number of hinge points for emissivity or reflectivity
c* status    character  directs routine to open a new or old netCDF file
c* irec      integer    Record number of profile
c* xid       character  Profile ID. (e.g., IPO11111)
c* lat       real       latitude 
c* lon       real       longitude 
c* time      integer    year-month-day-hour-minute-msecond
c* x         real       Profile - T, Q, PSFC, TSKIN, clouds, sfc emissivity
c* landtype  integer    Land type 1=land 0=ocean
c* pland     real       percent land (0...1)
C* eaa       real       earth azimuth angle per channel
c* eia       real       earth incidence angle per channel
c* surfalt   real       Surface altitude
c* roUnits   integer    Units for radiom. output (RAD_UNITS, RADMW_UNITS, or TB_UNITS)
c* radSim    real       simulated or measured brightness temperatures or radiance
c* radRetr   real       retrieved brightness temperatures or radiance
c* radRes    real       retrieved brightness temperatures or radiance
c* iter      integer    current iteration number
c* rms       real       r.m.s array for all the iterations
c* chisq     real       Chi-Square values for all the iterations
c* noise     real       standard noise 
c* norbit    integer    Satellite orbit number
c* nline     integer    Scan line of orbit
c* nelement  integer    Element in scan line
c* msg       string     closing message 
c* IG        StateIndex_t  state vector indices in geophysical space
c* NG        StateIndex_t  state vector lengths in geophysical space
c* qc        integer    integer array of quality control flags
c* 
c* Externals:
c* File_Name       Description
c* ---------       -----------
c* ncdf_module     base module of performing netCDF I/Os
c* 
c* Developed by AER 2003
c* Copyright, AER, Inc., 2003
c**********************************************************************
      use StateIndexModule
      use ncdf_module
      implicit none
      private
      integer :: imol, fidIndex=0, index
      integer , dimension(MAXPAR,MAXFILE) :: mapFrom
      integer , dimension(MAXFILE) :: fid=-1, mapSize, xFileSize
      logical , dimension(MAXFILE) :: stateIndexEnabled
      integer, dimension(:), pointer :: mapBeamLoc ! to argument not present
      integer :: nBeamLoc, nAngLoc, eiaLenLoc
      real, dimension(:), allocatable :: eiaLoc
      type(StateIndex_t) :: IGF, NGF
      integer, parameter, public :: TB_UNITS    = 0
      integer, parameter, public :: RAD_UNITS   = 1
      integer, parameter, public :: RADMW_UNITS = 2
      
      public openScene, getScene, putScene, closeScene, queryScene,
     $     appendScene, replaceScene
      interface openScene
      module procedure openNcFile
      end interface

      interface queryScene
      module procedure queryNcScene
      end interface

      interface getscene
      module procedure getNcScene
      end interface

      interface putscene
      module procedure putNcScene
      end interface

      interface appendScene
      module procedure appendNcScene
      end interface

      interface replaceScene
      module procedure replaceNcScene
      end interface

      interface closeScene
      module procedure closeNcFile
      end interface
c----
      CONTAINS
c--
      subroutine openNcFile(ncid, file, creationDate, nprf, pressure,  
     $     polarity, freq, casename, nAng, nLevel, nchmw, nSfcGrid, 
     $     IG, NG, MolID, xidLen, nQc, nBeam, mapBeam, nInner, nIter,
     $     sfcGrid, edrName, edrDim,
     $     vCoordTyp,sigPtop,sigCoord,hybCoordA,hybCoordB,
     $     status)
      
      use MapBeamDefaults
      use constants, only: MISSING_INT

c-- Dummy variables
      integer, intent(inout) :: ncid
      character (len=*), intent (in) :: file
      character (len=*), intent (inout), optional :: creationDate
      integer, intent(inout), optional :: nprf, nAng, nLevel, nchmw
      integer, intent(inout), optional :: nSfcGrid, nInner
      integer, intent(inout), optional :: xidLen, nQc, nBeam, nIter
      integer, dimension(:), intent(inout), optional :: mapBeam
      real, dimension(:), intent (inout), optional :: pressure
      real, dimension(:), intent (inout), optional :: freq
      real, dimension(:), intent (inout), optional :: sfcGrid
      integer, dimension(:), intent (inout), optional :: polarity, MolID
      character (len=*), intent(inout), optional :: casename
      character (len=*), intent(in), optional :: edrName
      integer, intent(inout), optional :: edrDim
      character (len=*),   intent (inout), optional :: vCoordTyp
      real,                intent (inout), optional :: sigPtop
      real, dimension (:), intent (inout), optional :: sigCoord
      real, dimension (:), intent (inout), optional :: hybCoordA
      real, dimension (:), intent (inout), optional :: hybCoordB
      character (len=*), intent(in), optional :: status
      type(StateIndex_t), intent(inout), optional :: IG, NG

c-- Local variables
      character (len=MAX_NAME_LENGTH) :: lstatus
      integer, dimension(8) :: sysvalues
      character (len= 8) :: sysdate
      character (len=10) :: systime
      character (len= 5) :: syszone
      integer, dimension(2) :: myDim
      integer , dimension(maxMol) :: MolIDF,MolLenF
      integer :: nchmwloc, nAngle
      integer :: nSfcGridLoc
      logical :: vCfound

      if (.not. present(status)) then
         lstatus = 'old'
      else
         lstatus = trim(status)
      endif

      if (trim(lstatus) == 'new' .or. trim(lstatus) == 'NEW') then
         if ((present(freq) .or. present(polarity)) .and.
     $        .not. (present(nchmw))) then
            print *,'err[scene_io_module::openNcFile]: missing nchmw'
            call errorHalt(1)
         end if

         if ((present(sfcGrid)) .and. .not. (present(nSfcGrid))) then
            print *,'err[scene_io_module::openNcFile]: missing nSfcGrid'
            call errorHalt(1)
         end if

         call openNcdfFile(ncid, file, status='new', 
     $        unlimited_dim_name='nProfiles')

         call initFileID(ncid=ncid,fidIndex=fidIndex,
     $        stateIndexEnabled=stateIndexEnabled,IG=IG,NG=NG,fid=fid)

         if (present(casename))
     $        call writeNcdfAttr(ncid,attrName='case',
     $        attr=trim(casename))

         call date_and_time(sysdate, systime, syszone, sysvalues)

         call writeNcdfAttr(ncid,attrName='Creationtime',
     $        attr=trim(sysdate)//trim(systime))

         if (stateIndexEnabled(fidIndex)) then
            if (.not. present(MolID)) then
               print *,'err[scene_io_module::openNcFile]: ',
     $              ' missing MolID in output.'
               call errorHalt(1)
            end if

            NGF = expandMol(NG,MolID)
            IGF = genIndices(NGF)

            call writeNcdfAttr(ncid,attrName='Temperature',
     $           attr=(/IGF%temp,NGF%temp/))

            call writeNcdfAttr(ncid,attrName='Tskin',
     $           attr=(/IGF%tskin,NGF%tskin/))

            call writeNcdfAttr(ncid,attrName='SurfacePressure',
     $           attr=(/IGF%psfc,NGF%psfc/))
         
            call writeNcdfAttr(ncid,attrName='Mol',
     $           attr=IGF%mol(1:maxMol))
         
            call writeNcdfAttr(ncid,attrName='MolLen',
     $           attr=NGF%mol(1:maxMol))
         
            call writeNcdfAttr(ncid,attrName='LiqCloud',
     $           attr=(/IGF%cldLiq,NGF%cldLiq/))
         
            call writeNcdfAttr(ncid,attrName='IceCloud',
     $           attr=(/IGF%cldIce,NGF%cldIce/))
         
            call writeNcdfAttr(ncid,attrName='SurfaceWinds',
     $           attr=(/IGF%wind,NGF%wind/))
         endif
         
         if (present(nLevel)) 
     $        call writeNcdfAttr(ncid,attrName='nlev',attr=nLevel)

         if (present(nchmw)) then
            if (nchmw > 0) then
               nchmwLoc=nchmw
               call writeNcdfAttr(ncid,attrName='nchmw',attr=nchmwLoc)
               
               call mapBeamDefNew(nchmwLoc,nBeam,mapBeam,nBeamLoc,
     $              mapBeamLoc)

               if (present(nAng)) then
                  nAngLoc = nAng
               else
                  nAngLoc=nBeamLoc
               endif
               if (nAngLoc > 0)
     $              call writeNcdfDim(ncid,dimName='nAng',
     $              dimLen=nAngLoc)
               call writeNcdfAttr(ncid,attrName='nBeam',attr=nBeamLoc)
               call writeNcdfAttr(ncid,attrName='mapBeam', 
     $              attr=mapBeamLoc)
            endif
         endif

         if (present(nInner))
     $        call writeNcdfAttr(ncid,attrName='nInner',
     $        attr=nInner)

c         Number of hinge points for emissivity or reflectivity:            
         if (present(nSfcGrid)) 
     $           call writeNcdfAttr(ncid,attrName='nSfcGrid', 
     $           attr=nSfcGrid)

         if (present(vCoordTyp))
     &        call writeNcdfAttr(ncid,attrName='VertCoordType',
     &        attr=TRIM(vCoordTyp))

         if (present(pressure))
     $        call writeNcdfAttr(ncid,
     $        attrName='standardPressureGrid',
     $        attr=pressure)

         if (present(sigCoord))
     &        call writeNcdfAttr(ncid,
     &        attrName='SigmaCoord',
     &        attr=sigCoord)

         if (present(sigPtop))
     &        call writeNcdfAttr(ncid,
     &        attrName='SigmaTopPres',
     &        attr=sigPtop)

         if (present(hybCoordA))
     &        call writeNcdfAttr(ncid,
     &        attrName='HybridCoordA',
     &        attr=hybCoordA)

         if (present(hybCoordB))
     &        call writeNcdfAttr(ncid,
     &        attrName='HybridCoordB',
     &        attr=hybCoordB)

         if (present(sfcGrid)) then
            call writeNcdfAttr(ncid,attrName='sfcGrid',
     $           attr=sfcGrid)
         endif

         if (present(freq))
     $        call writeNcdfAttr(ncid,attrName='mwfrequencies',
     $        attr=freq)

         if (present(polarity))
     $        call writeNcdfAttr(ncid,attrName='mwpolarities',
     $        attr=polarity)

      else
         if (lstatus == 'replace' .or. lstatus == 'append') then
            call openNcdfFile(ncid, file, unlimited_dim_length=nprf,
     $           status='replace')
         else
            call openNcdfFile(ncid, file, unlimited_dim_length=nprf)
         endif

         call initFileID(ncid=ncid,fidIndex=fidIndex,
     $        stateIndexEnabled=stateIndexEnabled,IG=IG,NG=NG,fid=fid)

         if (present(creationDate))
     $        call readNcdfAttr(ncid,attrName='Creationtime',
     &        attr=creationDate)

         if (ncid <= 0) then 
            print *, 'warning[scene_io_module::openNcFile]: '//
     $           ' Invalid file ID'
            return
         endif

         if (present(xidLen)) xidLen=readNcdfDim(ncid,'nId',
     $        silent=.true.)

         if (present(nIter)) nIter=readNcdfDim(ncid,'TT_iteration',
     $        silent=.true.)

         if (present(nQc)) nQc=readNcdfDim(ncid,'nQC',silent=.true.)

         if (stateIndexEnabled(fidIndex)) then
            call readNcdfAttr(ncid,attr=myDim,attrName='Temperature')
            IGF%Temp = myDim(1)
            NGF%Temp = myDim(2)
         
            call readNcdfAttr(ncid,attrName='Tskin',attr=myDim)
            IGF%Tskin = myDim(1)
            NGF%Tskin = myDim(2)

            call readNcdfAttr(ncid,attrName='SurfacePressure',
     $           attr=myDim)
            IGF%Psfc = myDim(1)
            NGF%Psfc = myDim(2)

            call readNcdfAttr(ncid,attrName='Mol',attr=MolIDF)
            IGF%mol = MolIDF

            call readNcdfAttr(ncid,attrName='MolLen',attr=MolLenF)
            NGF%mol = MolLenF

            call readNcdfAttr(ncid,attrName='LiqCloud',attr=myDim)
            IGF%cldLiq = myDim(1)
            NGF%cldLiq = myDim(2)
         
            call readNcdfAttr(ncid,attrName='IceCloud',attr=myDim)
            IGF%cldIce = myDim(1)
            NGF%cldIce = myDim(2)

            call readNcdfAttr(ncid,attrName='SurfaceWinds',attr=myDim)
            IGF%Wind = myDim(1)
            NGF%Wind = myDim(2)
         endif

         call readNcdfAttr(ncid,attrName='nchmw',attr=nchmwLoc,
     $        silent=.true.)
         if (present(nchmw)) nchmw=nchmwLoc

         if (present(nInner))
     $        call readNcdfAttr(ncid,attrName='nInner',attr=nInner,
     $        silent=.true.)

         call readNcdfAttr(ncid,attrName='nBeam',attr=nBeamLoc, 
     $        silent=.true.)

         if (nchmwLoc /= 0) then
            allocate(mapBeamLoc(nchmwLoc))
            if (nBeamLoc == 0) then
               nBeamLoc = 1
               mapBeamLoc(1:nchmwLoc) = 1
            else
               call readNcdfAttr(ncid,attrName='mapBeam',
     $              attr=mapBeamLoc)
            endif
            if (present(mapBeam)) mapBeam(1:nchmwLoc) = mapBeamLoc
            deallocate(mapBeamLoc)
         else
            if (present(mapBeam)) mapBeam = MISSING_INT
         endif
         if (present(nBeam)) nBeam = nBeamLoc

         if (nBeamLoc == 0) then ! check the dependency of eia on beams
            eiaLenLoc=readNcdfDim(ncid,'eiaLen1',silent=.true.)
            if (eiaLenLoc > 0) then
               if (.not. allocated(eiaLoc))
     $              allocate(eiaLoc(eiaLenLoc))
               call readNcdfData(ncid,eiaLoc,varname='eia',
     $              record_no=1,status='single_record')
               if (any(eiaLoc(:) /= eiaLoc(1))) then
                  nAngLoc=size(eiaLoc)
                  if (allocated(eiaLoc)) deallocate(eiaLoc)
               else
                  nAngLoc = 1
               endif
            else
               nAngLoc = 0
            endif
         else
            nAngLoc = readNcdfDim(ncid,'nAng',silent=.TRUE.)
            if (nAngLoc == 0) nAngLoc=nBeamLoc
         endif
         
         if (present(nLevel))
     $        call readNcdfAttr(ncid,attrName='nlev',attr=nLevel)

         if (present(vCoordTyp)) then
            call readNcdfAttr(ncid,
     &        attrName='VertCoordType',
     &        attr=vCoordTyp,silent=.TRUE.,found=vCfound)
            if (.NOT. vCfound) vCoordTyp='pressureStd' ! for back compatibility
         endif
         
         if (present(pressure))
     $        call readNcdfAttr(ncid,
     $        attrName='standardPressureGrid',
     $        attr=pressure)

         if (present(sigCoord))
     &        call readNcdfAttr(ncid,
     &        attrName='SigmaCoord',
     &        attr=sigCoord)

         if (present(sigPtop))
     &        call readNcdfAttr(ncid,
     &        attrName='SigmaTopPres',
     &        attr=sigPtop)

         if (present(hybCoordA))
     &        call readNcdfAttr(ncid,
     &        attrName='HybridCoordA',
     &        attr=hybCoordA)

         if (present(hybCoordB))
     &        call readNcdfAttr(ncid,
     &        attrName='HybridCoordB',
     &        attr=hybCoordB)

         if (present(freq)) then
            call readNcdfAttr(ncid,attrName='mwfrequencies',attr=freq)
            if (size(freq) /= nchmwloc) then
               print *,'err[scene_io_module::openNcFile]: ',
     $              ' #freq and nchmw mismatched'
               call errorHalt(1)
            end if

         endif

         call readNcdfAttr(ncid,attrName='nSfcGrid',attr=nSfcGridLoc,
     $        silent=.true.)
         if (present(nSfcGrid)) nSfcGrid=nSfcGridLoc

         if (present(sfcGrid)) then
            call readNcdfAttr(ncid,attrName='sfcGrid',attr=sfcGrid)
            if (size(sfcGrid) /= nSfcGridloc) then
               print *,'err[scene_io_module::openNcFile]: ',
     $              ' #sfcGrid and nSfcGrid mismatched'
               call errorHalt(1)
            end if
         endif

         if (present(polarity)) then
            call readNcdfAttr(ncid,attrName='mwpolarities',
     $           attr=polarity)
            if (size(polarity) /= nchmwloc) then
               print *,'err[scene_io_module::openNcFile]: ',
     $              ' #mwpol and nchmw mismatched'
               call errorHalt(1)
            end if

         endif

         if (present(casename))
     $        call readNcdfAttr(ncid,attrName='case',
     $        attr=casename,silent=.true.)

         if (present(nAng)) then
            !The process of reading nAngle is for back-compatibility
            call readNcdfAttr(ncid,attrName='numberofangles',
     $           attr=nAngle,silent=.true.)
            if (nAngle /= 0) then
               if (nchmwLoc /= 0) then
                  if (nAngle /= nAngLoc) then
                     print*,'warning[scene_io_module::openNcFile]: '//
     $                    ' nAngle in the old file ignored'
                  endif
                  nAng = nAngLoc
               else
                  nAng = nAngle
               endif
            else
               if (nchmwLoc /= 0) then
                  nAng = nAngLoc
               else
                  nAng = readNcdfDim(ncid,'nAng',silent=.TRUE.)
               endif
            endif   
         endif
         
         if (present(edrName) .and. present(edrDim))
     $        edrDim=readNcdfNDim(fid=ncid,varName=edrName,
     $        silent=.TRUE.)

         if (stateIndexEnabled(fidIndex)) then
            call genMapFrom(IG=IG,NG=NG,IGF=IGF,NGF=NGF,MolID=MolID,
     $           map=mapFrom(:,fidIndex),mapSize=mapSize(fidIndex),
     $           xFileSize=xFileSize(fidIndex))
         endif
      endif

      end subroutine openNcFile
c--
      subroutine queryNcScene(file, creationDate, nprf, pressure,
     $     polarity, freq, casename, nLevel, nchmw, nSfcGrid, 
     $     IG, NG, Nmol, xidLen, nQc, MolID, nBeam, mapBeam, nang, 
     $     nInner, nIter, sfcGrid, edrName, edrDim,
     $     vCoordTyp,sigPtop,sigCoord,hybCoordA,hybCoordB)

c-- Dummy variables
      character (len=*), intent (in) :: file
      character (len=*), intent (inout), optional :: creationDate
      integer, intent(inout), optional :: nprf, nLevel, nchmw
      integer, intent(inout), optional :: nSfcGrid, nInner
      integer, intent(inout), optional :: xidLen,nQc,nBeam,nAng,nIter
      integer, dimension(:), intent(inout), optional :: mapBeam
      real, dimension(:), intent (inout), optional :: pressure
      real, dimension(:), intent (inout), optional :: freq
      real, dimension(:), intent (inout), optional :: sfcGrid
      integer, dimension(:), intent (inout), optional :: polarity, MolID
      character (len=*), intent(inout), optional :: casename
      type(StateIndex_t), intent(inout), optional :: IG, NG
      integer , intent(inout), optional :: NMol
      character (len=*), intent(in), optional :: edrName
      integer, intent(inout), optional :: edrDim
      character (len=*),   intent (inout), optional :: vCoordTyp
      real,                intent (inout), optional :: sigPtop
      real, dimension (:), intent (inout), optional :: sigCoord
      real, dimension (:), intent (inout), optional :: hybCoordA
      real, dimension (:), intent (inout), optional :: hybCoordB

c-- Local variables
      integer :: ncid
      
      if (present(MolID)) MolID=nullMolID
      if (present(NMol)) NMol=0
      call openNcFile(ncid=ncid, file=file, creationDate=creationDate, 
     $     nprf=nprf, pressure=pressure, polarity=polarity, freq=freq, 
     $     casename=casename, nAng=nAng, nLevel=nLevel, 
     $     nSfcGrid=nSfcGrid,  nchmw=nchmw, IG=IG, NG=NG, 
     $     xidLen=xidLen, nQc=nQc, MolID=MolID, nBeam=nBeam, 
     $     mapBeam=mapBeam, sfcGrid=sfcGrid, edrName=edrName, 
     $     edrDim=edrDim,  nInner=nInner, nIter=nIter, 
     $     vCoordTyp=vCoordTyp,sigPtop=sigPtop,
     $     sigCoord=sigCoord,hybCoordA=hybCoordA,hybCoordB=hybCoordB)
      call closeNcFile(ncid)
      end subroutine queryNcScene
c--
      subroutine getNcScene(ncid, irec, xid, lat, lon, time, 
     $     x, landtype, pland, eaa, eia, 
     $     pressure, radSim, radRetr, radRes,
     $     iter, rms, chisq, noise, norbit, nline, nelement,
     $     kchanF,surfalt,EmMw, EmIr, EmRf, qc, edr, edrName, edrV)

      integer, intent (in) :: ncid
      integer, intent (in) :: irec
      character (len=*), intent (inout), optional :: xid
      real, intent (inout), optional :: lat, lon, surfalt
      integer, dimension(:), intent (inout), optional :: time, kchanF
      real, dimension(:), intent (inout), optional :: x
      real, intent (inout), optional :: pland
      integer, intent (inout), optional :: landtype
      real, dimension (:), intent (inout), optional :: eaa, eia
      real, dimension (:), intent (inout), optional :: pressure
      real, dimension (:), intent (inout), optional :: 
     $     radSim, radRetr, radRes
      integer, intent (inout), optional :: iter
      real, dimension(:), intent (inout), optional :: rms, chisq
      real, intent (inout), optional :: noise
      integer, intent (inout), optional :: norbit, nline, nelement
      real , dimension (:), intent(inout), optional :: EmMW, EmIR, EmRf
      real , dimension (MAXPAR) :: xFile
      integer, dimension(:), intent(inout), optional :: qc
      real, intent(inout), optional :: edr
      character (len=*), intent(in), optional :: edrName
      real, intent(inout), dimension(:), optional :: edrV

      if (present(x)) then
         index = searchID(fid,ncid)
         if (.not. stateIndexEnabled(index)) then
            print *,'err[scene_io_module::getNcScene]: ',
     $           ' stateIndex not enabled.'
            call errorHalt(1)
         end if

         call readNcdfData(ncid,xFile(1:xFileSize(index)),
     $        varname='profiles',
     $        record_no=irec,status='single_record')
         if (mapSize(index) /= size(x)) then
            print *,'err[scene_io_module::getNcScene]: ',
     $           ' map size not comfortable'
            call errorHalt(1)
         end if

         x = xFile(mapFrom(1:mapSize(index),index))
      endif

      if (present(xid))
     $     call readNcdfData(ncid,xid,varname='id',
     $     record_no=irec,status='single_record')

      if (present(lat))
     $     call readNcdfData(ncid,lat,varname='lat',
     $     record_no=irec,status='single_record')

      if (present(lon))
     $     call readNcdfData(ncid,lon,varname='lon',
     $     record_no=irec,status='single_record')

      if (present(time))
     $     call readNcdfData(ncid,time,varname='date',
     $     record_no=irec,status='single_record')

      if (present(landtype))
     $     call readNcdfData(ncid,landtype,varname='Landtype',
     $     record_no=irec,status='single_record')

      if (present(pland))
     $     call readNcdfData(ncid,pland,varname='pland',
     $     record_no=irec,status='single_record')

      if (present(eia)) then
         eiaLenLoc=readNcdfDim(ncid,'eiaLen1',silent=.true.)
         if (size(eia) == 1 .and. eiaLenLoc > 1) then
            if (.not. allocated(eiaLoc)) allocate(eiaLoc(eiaLenLoc))
            call readNcdfData(ncid,eiaLoc,varname='eia',
     $           record_no=irec,status='single_record')
            eia = eiaLoc(1)
         else
            call readNcdfData(ncid,eia,varname='eia',
     $           record_no=irec,status='single_record')
         endif
      endif

      if (present(eaa))
     $     call readNcdfData(ncid,eaa,varname='eaa',
     $     record_no=irec,status='single_record')

      if (present(pressure))
     $     call readNcdfData(ncid,pressure,varname='pressure',
     $     record_no=irec,status='single_record')

      if (present(edr) .and. present(edrName)) then
         call readNcdfData(ncid,edr,varname=edrName,
     $        record_no=irec,status='single_record')
      endif

      if (present(edrV) .and. present(edrName)) then
         call readNcdfData(ncid,edrV,varname=edrName,
     $        record_no=irec,status='single_record')
      endif

      if (present(rms))
     $     call readNcdfData(ncid,rms,varname='rms',
     $     record_no=irec,status='single_record')

      if (present(chisq))
     $     call readNcdfData(ncid,chisq,varname='chisq',
     $     record_no=irec,status='single_record')

      if (present(iter))
     $     call readNcdfData(ncid,iter,varname='iteration',
     $     record_no=irec,status='single_record')

      if (present(noise))
     $     call readNcdfData(ncid,noise,varname='noise',
     $     record_no=irec,status='single_record')

      if (present(kchanF))
     $     call readNcdfData(ncid,kchanF,varname='kchanF',
     $     record_no=irec,status='single_record')

      if (present(radSim))
     $     call readNcdfData(ncid,radSim,varname='simrad',
     $     record_no=irec,status='single_record')

      if (present(radRetr))
     $     call readNcdfData(ncid,radRetr,varname='retrrad',
     $     record_no=irec,status='single_record')

      if (present(radRes))
     $     call readNcdfData(ncid,radRes,varname='resrad',
     $     record_no=irec,status='single_record')

      if (present(norbit))
     $     call readNcdfData(ncid,norbit,varname='norbit',
     $     record_no=irec,status='single_record')

      if (present(nline))
     $     call readNcdfData(ncid,nline,varname='nline',
     $     record_no=irec,status='single_record')

      if (present(nelement))
     $     call readNcdfData(ncid,nelement,varname='nelement',
     $     record_no=irec,status='single_record')

      if (present(surfalt)) 
     &     call readNcdfData(ncid,surfalt,varname='surfalt',
     &     record_no=irec,status='single_record')

      if (present(EmMw)) 
     &     call readNcdfData(ncid,EmMw,varname='EmMw',
     &     record_no=irec,status='single_record')

      if (present(EmIr)) 
     &     call readNcdfData(ncid,EmIr,varname='EmIr',
     &     record_no=irec,status='single_record')

      if (present(EmRf)) 
     &     call readNcdfData(ncid,EmRf,varname='EmRf',
     &     record_no=irec,status='single_record')

      if (present(qc))
     $     call readNcdfData(ncid,qc,varname='qc',
     $     record_no=irec,status='single_record')

      end subroutine getNcScene
c---- 
      subroutine putNcScene(ncid, xid, lat, lon, time, x,
     $     landtype, pland, eaa, eia, pressure,
     $     radSim, radRetr, radRes,
     $     roUnits, iter, rms, chisq, diagError,dofs,
     $     atmosclass, noise, 
     $     norbit, nline, nelement,
     $     kchanF,surfalt,EmMw,EmIr,EmRf,qc,edr,edrName,edrLongName,
     $     edrUnit,edrV)

      integer, intent (in) :: ncid
      character (len=*), intent (in), optional :: xid
      real, intent (in), optional :: lat, lon, surfalt
      integer, dimension(:), intent (in), optional :: time, kchanF
      real, dimension(:), intent (in), optional :: x
      integer, intent (in), optional :: landtype
      real, intent (in), optional :: pland
      real, dimension (:), intent (in), optional :: eaa, eia
      real, dimension(:), intent (in), optional :: pressure, diagError
      real, dimension(:), intent (in), optional :: 
     $     radSim, radRetr, radRes
      integer, intent (in), optional :: iter, atmosclass
      integer, intent (in), optional :: roUnits
      real, dimension(:), intent (in), optional :: rms, chisq
      real, intent (in), optional :: noise, dofs
      integer, intent (in), optional :: norbit, nline, nelement
      real , dimension (:), intent (in), optional :: EmMw, EmIr, EmRf
      integer, dimension(:), intent(in), optional :: qc
      real, intent (in), optional :: edr
      character (len=*), intent (in), optional :: 
     $     edrName, edrLongName, edrUnit
      real, intent (in), dimension(:), optional :: edrV
      character (len=MAX_NAME_LENGTH) :: localLongName, localUnit
      character (len=16) :: rUnits ! Units for any radiometric output variables

      rUnits='K'
      if (present(roUnits)) then
         select case (roUnits)
         case (RAD_UNITS)
            rUnits='mW/m2/str/cm-1'
         case (RADMW_UNITS)
            rUnits='mW/m2/str/GHz'
         case (TB_UNITS)
            rUnits='K'
         case default
            print *,'err[scene_io_module::openNcFile]:'//
     $           ' Invalid roUnits:',roUnits
            print *,'Choose from RAD_UNITS, RADMW_UNITS, or TB_UNITS'
            call errorHalt(1)
         end select
      endif

      if (present(x)) then
         call writeNcdfData(ncid,x,varname='profiles',
     $        varlenName=(/'nPar     ','nProfiles'/),
     $        varLongName='set of profile vectors',
     $        varUnit='Temp(Kelvin), H2O(g/g) in mol(1), psfc(mb)',
     $        status='append')
      endif

      if (present(xid))
     $     call writeNcdfData(ncid,xid,varname='id',
     $     strLenName='nId',
     $     varLongName='profile identification',
     $     status='append')
      
      if (present(lat))
     $     call writeNcdfData(ncid,lat,varname='lat',
     $     varLongName='latitude values',
     $     varUnit='degrees north',
     $     status='append')

      if (present(lon))
     $     call writeNcdfData(ncid,lon,varname='lon',
     $     varLongName='longitude',
     $     varUnit='degrees east',
     $     status='append')

      if (present(surfalt)) 
     &     call writeNcdfData(ncid,surfalt,varname='surfalt',
     &     varLongName='Surface Altitude',
     &     varUnit='m', status='append')
      
      if (present(landtype))
     $     call writeNcdfData(ncid,landtype,varname='Landtype',
     $     varLongName='land type classification (UNDEF=-1) ',
     $     status='append')

      if (present(pland))
     $     call writeNcdfData(ncid,pland,varname='pland',
     $     varLongName='% land (0..1) ',
     $     status='append')

      if (present(eaa))
     $     call writeNcdfData(ncid,eaa,varname='eaa',
     $     varlenName=(/'nAng     ','nProfiles'/),
     $     varLongName='Earth Azimuth Angle',
     $     varUnit='degrees',
     $     status='append')

      if (present(eia))
     $     call writeNcdfData(ncid,eia,varname='eia',
     $     varlenName=(/'nAng     ','nProfiles'/),
     $     varLongName='Earth Incidence Angle',
     $     varUnit='degrees',
     $     status='append')

      if (present(pressure))
     $     call writeNcdfData(ncid,pressure,varname='pressure',
     $     varLongName='Pressure',
     $     varUnit='mb',
     $     status='append')

      if (present(time))
     $     call writeNcdfData(ncid,time,varname='date',
     $     varlenName=(/'nTime    ','nProfiles'/),
     $     varLongName='date and time',
     $     varUnit='year month day hour min msec',
     $     status='append')
      
      if (present(EmMw)) 
     &     call writeNcdfData(ncid,EmMw,varname='EmMw',
     $     varlenName=(/'nchmw    ','nProfiles'/),
     &     varLongName='MW surface emissivity',
     &     varUnit='0..1', status='append')
      
      if (present(EmIr)) 
     &     call writeNcdfData(ncid,EmIr,varname='EmIr',
     $     varlenName=(/'nchir    ','nProfiles'/),
     &     varLongName='IR surface emissivity',
     &     varUnit='0..1', status='append')

      if (present(EmRf)) 
     &     call writeNcdfData(ncid,EmRf,varname='EmRf',
     $     varlenName=(/'nSfcGridx2','nProfiles '/),
     &     varLongName='IR surface emissivity and reflectivity',
     &     varUnit='0..1', status='append')

      if (present(edr) .and. present(edrName)) then
         if (present(edrLongName)) then
            localLongName = trim(edrLongName)
         else
            localLongName = trim(edrName)
         endif
         if (present(edrUnit)) then
            localUnit = trim(edrUnit)
         else
            localUnit = ""
         endif
         call writeNcdfData(ncid,edr,varname=edrName,
     $        varLongName=localLongName,
     $        varUnit=localUnit,status='append')
      endif

      if (present(edrV) .and. present(edrName))
     $     call writeNcdfData(ncid,edrV,varname=edrName,
     $     varLongName=edrLongName,
     $     varUnit=edrUnit,status='append')
      
      if (present(kchanF))
     $     call writeNcdfData(ncid,kchanF,varname='kchanF',
     $     varlenName=(/'nchmw    ','nProfiles'/),
     $     varLongName='Channel Selection Flag',
     $     status='append')

      if (present(radSim))
     $     call writeNcdfData(ncid,radSim,varname='simrad',
     $     varlenName=(/'nchmw    ','nProfiles'/),
     $     varLongName='Simulated Brightness Temperatures',
     $     varUnit=trim(rUnits),
     $     status='append')

      if (present(radRetr))
     $     call writeNcdfData(ncid,radRetr,varname='retrrad',
     $     varlenName=(/'nchmw    ','nProfiles'/),
     $     varLongName='Retrieved Brightness Temperatures',
     $     varUnit=trim(rUnits),
     $     status='append')

      if (present(radRes))
     $     call writeNcdfData(ncid,radRes,varname='resrad',
     $     varlenName=(/'nchmw    ','nProfiles'/),
     $     varLongName='Residual Brightness Temperatures',
     $     varUnit=trim(rUnits),
     $     status='append')

      if (present(norbit))
     $     call writeNcdfData(ncid,norbit,varname='norbit',
     $     varLongName='orbit number',
     $     status='append')
      
      if (present(nline))
     $     call writeNcdfData(ncid,nline,varname='nline',
     $     varLongName='scan line number of orbit',
     $     status='append')
      
      if (present(nelement))
     $     call writeNcdfData(ncid,nelement,varname='nelement',
     $     varLongName='element of scan line',
     $     status='append')
      
      if (present(rms))
     $     call writeNcdfData(ncid,rms,varname='rms',
     $     varlenName=(/'TT_iteration','nProfiles   '/),
     $     varLongName='root mean square residual ',
     $     varUnit=trim(rUnits),
     $     status='append')

      if (present(chisq))
     $     call writeNcdfData(ncid,chisq,varname='chisq',
     $     varlenName=(/'TT_iteration','nProfiles   '/),
     $     varLongName='Chi-Squared ',
     $     status='append')


      if (present(iter))
     $     call writeNcdfData(ncid,iter,varname='iteration',
     $     varLongName='number of iterations for convergence',
     $     status='append')

      if (present(atmosclass))
     $     call writeNcdfData(ncid,atmosclass,varname='atmosclass',
     $     varLongName='A prior profile classification',
     $     status='append')

      if (present(diagError))
     $     call writeNcdfData(ncid,diagError,varname='diagError',
     $     varLongName='Retrieved error covariance diagonal', 
     $     status='append')

      if (present(dofs))
     $     call writeNcdfData(ncid,dofs,varname='DOFS',
     $     varLongName='Degrees of freedom for signal', 
     $     status='append')

      if (present(noise))
     $     call writeNcdfData(ncid,noise,varname='noise',
     $     varLongName='total noise',
     $     varUnit='Kelvin',
     $     status='append')

      if (present(qc))
     $     call writeNcdfData(ncid,qc,varname='qc',
     $     varlenName=(/'nQC      ','nProfiles'/),
     $     varLongName='Quality control byte array',
     $     status='append')

      end subroutine putNcScene
c---- 
      subroutine appendNcScene(ncid, xid, lat, lon, time, x,
     $     landtype, pland, eaa, eia, pressure,
     $     radSim, radRetr, radRes,
     $     iter, rms, chisq, diagError, dofs,
     $     atmosclass, 
     $     noise, norbit, nline, nelement,
     $     kchanF,surfalt,EmMw,EmIr,qc,
     $     edr,edrName,edrLongName,edrUnit,edrV)

      integer, intent (in) :: ncid
      real, dimension(:), intent (in), optional :: x
      character (len=*), intent (in), optional :: xid
      real, intent (in), optional :: lat, lon, surfalt
      integer, dimension(:), intent (in), optional :: time, kchanF
      integer, intent (in), optional :: landtype
      real, intent (in), optional :: pland
      real, dimension (:), intent (in), optional :: eaa, eia
      real, dimension (:), intent (in), optional :: pressure,diagError
      real, dimension(:), intent (in), optional :: 
     $     radSim, radRetr, radRes
      integer, intent (in), optional :: iter, atmosclass
      real, dimension(:), intent (in), optional :: rms, chisq
      real, intent (in), optional :: noise, dofs
      integer, intent (in), optional :: norbit, nline, nelement
      real , dimension (:), intent (in), optional :: EmMw
      real , dimension (:), intent (in), optional :: EmIr
      integer, dimension(:), intent(in), optional :: qc
      real, intent (in), optional :: edr
      character (len=*), intent (in), optional :: 
     $     edrName, edrLongName, edrUnit
      real, intent (in), dimension(:), optional :: edrV
      integer :: rec

      rec = readNcdfDim(fid=ncid)+1 !get the new unlimited dimension

      call replaceNcScene(ncid=ncid, irec=rec, xid=xid, lat=lat, 
     $     lon=lon, time=time, x=x, landtype=landtype, pland=pland, 
     $     eaa=eaa, eia=eia, pressure=pressure,
     $     radSim=radSim,radRetr=radRetr,radRes=radRes,
     $     iter=iter, rms=rms, chisq=chisq, diagError=diagError,
     $     dofs=dofs, atmosclass=atmosclass,
     $     noise=noise, norbit=norbit, 
     $     nline=nline, nelement=nelement, kchanF=kchanF, 
     $     surfalt=surfalt, EmMw=EmMw, EmIr=EmIr, qc=qc, 
     $     edr=edr, edrName=edrName, 
     $     edrLongName=edrLongName, edrUnit=edrUnit, edrV=edrV)

      end subroutine appendNcScene
c---- 
      subroutine replaceNcScene(ncid, irec, xid, lat, lon, time, x,
     $     landtype, pland, eaa, eia, pressure,
     $     radSim, radRetr, radRes,
     $     iter, rms, chisq, diagError, dofs,
     $     atmosclass,
     $     noise, norbit, nline, nelement,
     $     kchanF,surfalt,EmMw,EmIr,qc,
     $     edr,edrName,edrLongName,edrUnit,edrV)

      integer, intent (in) :: ncid
      integer, intent (in), optional :: irec
      character (len=*), intent (in), optional :: xid
      real, intent (in), optional :: lat, lon, surfalt
      integer, dimension(:), intent (in), optional :: time, kchanF
      real, dimension(:), intent (in), optional :: x
      integer, intent (in), optional :: landtype
      real, intent (in), optional :: pland
      real, dimension (:), intent (in), optional :: eaa, eia
      real, dimension (:), intent (in), optional :: pressure,diagError
      real, dimension(:), intent (in), optional :: 
     $     radSim, radRetr, radRes
      integer, intent (in), optional :: iter, atmosclass
      real, dimension(:), intent (in), optional :: rms, chisq
      real, intent (in), optional :: noise, dofs
      integer, intent (in), optional :: norbit, nline, nelement
      real , dimension (:), intent (in), optional :: EmMw
      real , dimension (:), intent (in), optional :: EmIr
      integer, dimension(:), intent(in), optional :: qc
      real, intent (in), optional :: edr
      character (len=*), intent (in), optional :: 
     $     edrName, edrLongName, edrUnit
      real, intent (in), dimension(:), optional :: edrV

      if (present(x))
     $     call replaceNcdfData(fid=ncid,var=x,varname='profiles',
     $     rec=irec,status='append')

      if (present(xid))
     $     call replaceNcdfData(fid=ncid,var=xid,varname='id',
     $     rec=irec,status='append')
      
      if (present(lat))
     $     call replaceNcdfData(fid=ncid,var=lat,varname='lat',
     $     rec=irec,status='append')

      if (present(lon))
     $     call replaceNcdfData(fid=ncid,var=lon,varname='lon',
     $     rec=irec,status='append')

      if (present(surfalt)) 
     &     call replaceNcdfData(fid=ncid,var=surfalt,varname='surfalt',
     $     rec=irec,status='append')
      
      if (present(landtype))
     $     call replaceNcdfData(fid=ncid,var=landtype,
     $     varname='Landtype',rec=irec,status='append')

      if (present(pland))
     $     call replaceNcdfData(fid=ncid,var=pland,varname='pland',
     $     rec=irec,status='append')

      if (present(eaa))
     $     call replaceNcdfData(fid=ncid,var=eaa,varname='eaa',
     $     rec=irec,status='append')

      if (present(eia))
     $     call replaceNcdfData(fid=ncid,var=eia,varname='eia',
     $     rec=irec,status='append')

      if (present(pressure))
     $     call replaceNcdfData(fid=ncid,var=pressure,
     $     varname='pressure',
     $     rec=irec,status='append')

      if (present(time))
     $     call replaceNcdfData(fid=ncid,var=time,varname='date',
     $     rec=irec,status='append')
      
      if (present(EmMw)) 
     &     call replaceNcdfData(fid=ncid,var=EmMw,varname='EmMw',
     $     rec=irec,status='append')
      
      if (present(EmIr)) 
     &     call replaceNcdfData(fid=ncid,var=EmIr,varname='EmIr',
     $     rec=irec,status='append')

      if (present(edr) .and. present(edrName))
     $     call replaceNcdfData(fid=ncid,var=edr,varname=edrName,
     $     rec=irec,status='append')
      
      if (present(edrV) .and. present(edrName))
     $     call replaceNcdfData(fid=ncid,var=edrV,varname=edrName,
     $     rec=irec,status='append')
      
      if (present(kchanF))
     $     call replaceNcdfData(fid=ncid,var=kchanF,varname='kchanF',
     $     rec=irec,status='append')

      if (present(radSim))
     $     call replaceNcdfData(fid=ncid,var=radSim,varname='simrad',
     $     rec=irec,status='append')

      if (present(radRetr))
     $     call replaceNcdfData(fid=ncid,var=radRetr,varname='retrrad',
     $     rec=irec,status='append')

      if (present(radRes))
     $     call replaceNcdfData(fid=ncid,var=radRes,varname='resrad',
     $     rec=irec,status='append')

      if (present(norbit))
     $     call replaceNcdfData(fid=ncid,var=norbit,varname='norbit',
     $     rec=irec,status='append')
      
      if (present(nline))
     $     call replaceNcdfData(fid=ncid,var=nline,varname='nline',
     $     rec=irec,status='append')
      
      if (present(nelement))
     $     call replaceNcdfData(fid=ncid,var=nelement,
     $     varname='nelement',rec=irec,status='append')
      
      if (present(rms))
     $     call replaceNcdfData(fid=ncid,var=rms,varname='rms',
     $     rec=irec,status='append')

      if (present(chisq))
     $     call replaceNcdfData(fid=ncid,var=chisq,varname='chisq',
     $     rec=irec,status='append')

      if (present(iter))
     $     call replaceNcdfData(fid=ncid,var=iter,varname='iteration',
     $     rec=irec,status='append')

      if (present(diagError))
     $     call writeNcdfData(ncid,diagError,varname='diagError',
     $     varLongName='Retrieved error covariance diagonal', 
     $     status='append')

      if (present(dofs))
     $     call writeNcdfData(ncid,dofs,varname='DOFS',
     $     varLongName='Degrees of freedom for signal', 
     $     status='append')

      if (present(atmosclass))
     $     call replaceNcdfData(fid=ncid,var=atmosclass,
     $     varname='atmosclass',
     $     rec=irec,status='append')

      if (present(noise))
     $     call replaceNcdfData(fid=ncid,var=noise,varname='noise',
     $     rec=irec,status='append')

      if (present(qc))
     $     call replaceNcdfData(fid=ncid,var=qc,varname='qc',
     $     rec=irec,status='append')

      end subroutine replaceNcScene
c--
      subroutine closeNcFile(ncid,msg)
      character (len=*), intent (in), optional :: msg
      integer, intent (in) :: ncid

      call closeNcdfFile(ncid, msg)

      end subroutine closeNcFile
c----
      end module scene_io_module
