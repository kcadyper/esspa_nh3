module RDR_io_module
!**********************************************************************
!* Module: RDR_io_module
!* Purpose: Simulated Raw Data Record (RDR) I/O
!* Description: The module is used to open/close, read/write RDR data
!* to a netCDF file.
!* This module treats some geolocation data that would not be present in 
!* a real RDR file, but that are convenient to have in simulation
!*            
!* Inputs/Outputs:
!* Var_name      Type       Description
!* --------      ----       -----------
!* ncid          integer    NetCDF file ID (or descriptor)
!* fname         character  Filename for output.
!* norbit        integer    Satellite orbit number
!* nline         integer    Scan line of orbit
!* satalt        real       Satellite altitude
!* surfalt       real       Surface altitude
!* lat           real       FOV latitude
!* lon           real       FOV longitude                        
!* time          integer    Year-month-day-hour-minute-msecond
!* eaa           real       Earth azimuth angle per channel
!* eia           real       Earth incidence angle per channel
!* sza           real       Solar zenith angle
!* saa           real       Solar azimuth angle
!* pol           integer    channel polarities
!* frq           real       channel frequencies (GHz)
!* nchan         integer    number of sensor channels
!* nthmr         integer    number of thermistors
!* nthmTot       integer    totoal number of thermal readings per scan
!* nscid         integer    Spacecraft ID
!* instid        integer    Sensor ID
!* modver        integer    Sensor model configuration number
!* nfilever      integer    Sensor constants file version
!* qflg          integer    Data quality flag
!* npos          integer    Number of elements per scan line
!* speed         float      Satellite foward speed (nominal)
!* view          integer    Satellite viewing direction (forward/aft 1/0)
!* chanID        character  string that uniquely identifies each channel
!* thermID       character  string that uniquely identifies each thermistor

!* extended list for GeoRDRintfce
!* sel_sc        real       solar elevation angle w.r.t. S/C
!* saz_sc        real       solar azimuth angle w.r.t. S/C
!* mel_sc        real       lunar elevation angle w.r.t. S/C
!* maz_sc        real       lunar azimuth angle w.r.t. S/C
!* ang_mncs      real       lunar-viewing angle w.r.t. cold-sky bore sight
!*
!* Externals: ncdf_module
!* Copyright: 2003, AER, Inc.  
!* Developed by AER, Inc., 2003
!**********************************************************************

  implicit none
  private
  public :: openRDR,closeRDR,putRDR,getRDR,queryRDR
  
  integer , parameter            :: MAX_FILE=50
  integer , save                 :: irec=0
  logical , dimension (MAX_FILE) :: newForm=.false.
  real , parameter               :: defaultFloat=-99999.
  
  interface openRDR
    module procedure openNcRDR
  end interface

  interface queryRDR
    module procedure queryNcRDR
  end interface

  interface getRDR
    module procedure getNcRDR
  end interface

  interface putRDR
    module procedure putNcRDR
  end interface

  interface closeRDR
    module procedure closeNcRDR
  end interface

CONTAINS
      
  subroutine openNcRDR(ncid,fname,status,frq,pol,nrec,nchan, &
    nBeam,mapBeam,npos,nccals,nwcals,nthmrps,nthmr,nthmTot, &
    nqc,qflg,nscid,instid,modver,nfilever, &
    speed,satalt,view,createTime,chanID,thermID)

    use ncdf_module
    use MapBeamDefaults
  
    !---in/out variables
    integer,               intent(inout)           :: ncid
    character(len=*),      intent(in)              :: fname
    character(len=*),      intent(in),    optional :: status
    real,    dimension(:), intent(inout), optional :: frq
    integer, dimension(:), intent(inout), optional :: pol
    integer,               intent(inout), optional :: nrec
    integer,               intent(inout), optional :: nchan,nBeam
    integer, dimension(:), intent(inout), optional :: mapBeam,npos
    integer, dimension(:), intent(inout), optional :: nccals,nwcals
    integer,               intent(inout), optional :: nthmrps
    integer,               intent(inout), optional :: nthmr
    integer,               intent(inout), optional :: nthmTot
    integer,               intent(inout), optional :: nqc
    integer, dimension(:), intent(inout), optional :: qflg
    integer,               intent(inout), optional :: nscid
    integer,               intent(inout), optional :: instid,modver,nfilever
    real,                  intent(inout), optional :: speed,satalt
    integer,               intent(inout), optional :: view
    character(len=*),      intent(inout), optional :: createTime
    character(len=*),dimension(:),intent(inout),optional :: chanID
    character(len=*),dimension(:),intent(inout),optional :: thermID
    
    !---Local variables
    integer, dimension(8) :: sysvalues
    character(len=20)     :: sysdate, systime, syszone, myStatus
    integer               :: npts, nbmt, ich
    integer               :: nccas,nwcas
    integer               :: nchanloc, nBeamloc      ! Local since cannot write
    integer, dimension(:), pointer :: mapBeamloc     ! to argument not present
    integer               :: nthmTotloc
    integer :: lenID,lenTot,iid1,iid2
    integer, parameter    :: maxIDlen=4000
    character(len=maxIDlen) :: chanIDall,thermIDall
    integer :: i

    if (present(status)) then
       myStatus = status
    else
       myStatus='old'
    endif
    if (trim(myStatus) == 'new' .or. trim(myStatus) == 'NEW' .or. &
         trim(myStatus) == 'New') then

       if (.not. present(nchan)) then
          print *,'err[RDR_io_module::openNcRDR]: ', &
               ' missing nchan for new (write) mode.'
          call errorHalt(1)
       end if

       if (.not. present(nthmTot)) then
          print *,'err[RDR_io_module::openNcRDR]: ', &
               'missing nthmTot for new (write) mode.'
          call errorHalt(1)
       end if

       call mapBeamDefNew(nchan,nBeam,mapBeam,nBeamloc,mapBeamloc)
       
       if (.not. present(npos)) then
          print *,'err[RDR_io_module::openNcRDR]: ', &
               ' missing npos for new (write) mode.'
          call errorHalt(1)
       elseif (size(npos) < nBeamloc) then
          print *,'err[RDR_io_module::openNcRDR]: ', &
               ' size(npos) < nBeam :',size(npos),nBeamloc
          call errorHalt(1)
       end if
       if (.not. present(nccals)) then
          print *,'err[RDR_io_module::openNcRDR]: ', &
               ' missing nccals for new (write) mode.'
          call errorHalt(1)
       elseif (size(nccals) < nBeamloc) then
          print *,'err[RDR_io_module::openNcRDR]: ', &
               ' size(nccals) < nBeam :',size(nccals),nBeamloc
          call errorHalt(1)
       end if
       if (.not. present(nwcals)) then
          print *,'err[RDR_io_module::openNcRDR]: ', &
               ' missing nwcals for new (write) mode.'
          call errorHalt(1)
       elseif (size(nwcals) < nBeamloc) then
          print *,'err[RDR_io_module::openNcRDR]: ', &
               ' size(nwcals) < nBeam :',size(nwcals),nBeamloc
          call errorHalt(1)
       end if
       
       if (.not. present(nthmrps)) then
          print *,'err[RDR_io_module::openNcRDR]: ', &
               ' missing nthmrps for new (write) mode.'
          call errorHalt(1)
       end if
       
       call openNcdfFile(ncid, fname, status='new', &
            unlimited_dim_name='nRec')
       
       call date_and_time(sysdate, systime, syszone, sysvalues)
  
       call writeNcdfAttr(ncid,attrName='Creationtime', &
            attr=trim(sysdate//systime))
  
       if (present(nscid))  &
            call writeNcdfAttr(ncid,attrName='SpacecraftID', &
            attr=(/nscid/))
       
       if (present(instid))  &
            call writeNcdfAttr(ncid,attrName='SensorID',   &
            attr=(/instid/))
  
       if (present(modver))  &
            call writeNcdfAttr(ncid,attrName='SensorModelVersion', &
            attr=(/modver/))
  
       if (present(nfilever))  &
            call writeNcdfAttr(ncid,attrName='ConstantsFileVersion', &
            attr=(/nfilever/))

       if (present(qflg))  &
            call writeNcdfAttr(ncid,attrName='DataQualityFlags', &
            attr=int(qflg(1:nchan)))

       if (present(nthmr)) &
            call writeNcdfAttr(ncid,attrName='nthmr', &
            attr=(/nthmr/))

       call writeNcdfAttr(ncid,attrName='nchan',attr=(/nchan/))
       call writeNcdfAttr(ncid,attrName='nBeam',attr=(/nBeamloc/))
       call writeNcdfAttr(ncid,attrName='npos', attr=npos(1:nBeamloc))
       call writeNcdfAttr(ncid,attrName='nccals', attr=nccals(1:nBeamloc))
       call writeNcdfAttr(ncid,attrName='nwcals', attr=nwcals(1:nBeamloc))
       call writeNcdfAttr(ncid,attrName='nthmrps', attr=(/nthmrps/))
       call writeNcdfAttr(ncid,attrName='mapBeam', attr=mapBeamloc(1:nchan))

       ! explicitly determine summations to avoid potential
       !   compiler error:
       npts = 0
       do i=1,nchan
          npts=npts+npos(mapBeamloc(i))
       enddo
       nccas = 0
       do i=1,nchan
          nccas=nccas+nccals(mapBeamloc(i))
       enddo
       nwcas = 0
       do i=1,nchan
          nwcas=nwcas+nwcals(mapBeamloc(i))
       enddo

       nbmt = sum(npos(1:nBeamLoc))

       if (present(nQC)) &
            call writeNcdfDim(ncid,dimName='nQC',dimLen=nQC)
!!!       call writeNcdfDim(ncid,dimName='nRec',dimLen=nRec)
       call writeNcdfDim(ncid,dimName='nPts',dimLen=npts)
       call writeNcdfDim(ncid,dimName='nccas',dimLen=nccas)
       call writeNcdfDim(ncid,dimName='nwcas',dimLen=nwcas)
       call writeNcdfDim(ncid,dimName='nBmt',dimLen=nbmt)
       call writeNcdfDim(ncid,dimName='nthmTot',dimLen=nthmTot)
       
       if (present(frq))  &
            call writeNcdfAttr(ncid,attrName='mwfrequencies', &
            attr=frq(1:nchan))
  
       if (present(pol))  &
            call writeNcdfAttr(ncid,attrName='mwpolarities', &
            attr=pol(1:nchan))
  
       if (present(satalt))  &
            call writeNcdfAttr(ncid,attrName='SatelliteAltitude', &
            attr=(/satalt/))
  
       if (present(speed))  &
            call writeNcdfAttr(ncid, &
            attrName='SatelliteForwardSpeed', &
            attr=(/speed/))
  
       if (present(view))  &
            call writeNcdfAttr(ncid,attrName='ViewingDirection', &
            attr=(/view/))
  
       if (present(chanID)) then
          if (size(chanID) /= nchan) then
             print*, "err[RDR_io_module::openRDR]: ", &
                  "size(chanID) /= nchan:",size(chanID),nchan
             call errorHalt(1)
          endif
          lenID=len(chanID(1))
          lenTot=(lenID+1)*nchan
          if (lenTot > maxIDlen) then
             print *, "err[RDR_io_module::openRDR]: ", &
                  "Total length of chanID exceeds length of chanIDall"
             print*,'lenTot,maxIDlen:',lenTot,maxIDlen
             call errorHalt(1)
          endif
          call writeNcdfAttr(ncid,attrName='lenID',attr=lenID)
          chanIDall=repeat(' ',maxIDlen)
          do ich=1,nchan
             iid1=(ich-1)*(lenID+1)+1
             iid2=ich*(lenID+1)-1
             chanIDall(iid1:iid2)=chanID(ich)
          enddo
          call writeNcdfAttr(ncid,attrName='chanID', &
               attr=chanIDall(1:lenTot))
       endif

       if (present(thermID)) then
          lenID=len(thermID(1))
          lenTot=(lenID+1)*nthmTot/nthmrps
          if (lenTot > maxIDlen) then
             print *, "err[RDR_io_module::openRDR]: ", &
                  "Total length of thermID exceeds length of thermIDall"
             print*,'lenTot,maxIDlen:',lenTot,maxIDlen
             call errorHalt(1)
          endif
          call writeNcdfAttr(ncid,attrName='lenThmID',attr=lenID)
          thermIDall=repeat(' ',maxIDlen)
          do ich=1,nthmTot/nthmrps  
             iid1=(ich-1)*(lenID+1)+1
             iid2=ich*(lenID+1)-1
             thermIDall(iid1:iid2)=thermID(ich)
          enddo
          call writeNcdfAttr(ncid,attrName='thermID', &
               attr=thermIDall(1:lenTot))
       endif
  
       deallocate (mapBeamloc)

    else         !----------- read existing file -------------

       call openNcdfFile(ncid, trim(fname))

       call readNcdfAttr(ncid,attrName='nchan',attr=nchanloc)
       if (present(nchan)) nchan=nchanloc
       
       call readNcdfAttr(ncid,attrName='nBeam',attr=nBeamloc)
       if (present(nBeam)) nBeam=nBeamloc
  
       if (present(createTime)) &
            call readNcdfAttr(ncid,attrName='Creationtime', &
            attr=createTime)
  
       if (present(nrec)) &
         nrec=readNcdfDim(ncid,DimName='nRec')

       if (present(nqc)) nqc=readNcdfDim(ncid,'nQC')
       
       nthmTotLoc=readNcdfDim(ncid,'nthmTot')
       if (present(nthmTot)) nthmTot=nthmTotLoc
 
       if (present(nthmr)) &
            call readNcdfAttr(ncid,attrName='nthmr', &
            attr=nthmr)

       if (present(nscid))  &
            call readNcdfAttr(ncid,attrName='Spacecraft ID', &
            attr=nscid)

       if (present(instid))  &
            call readNcdfAttr(ncid,attrName='Sensor ID', &
            attr=instid)
  
       if (present(modver))  &
            call readNcdfAttr(ncid,attrName='Sensor Model Version', &
            attr=modver)
  
       if (present(nfilever))  &
            call readNcdfAttr(ncid,attrName='Constants File Version', &
            attr=nfilever)
  
       if (present(qflg)) then
          if (size(qflg) < nchanloc) then
             print *,'err[RDR_io_module::openNcRDR]: ', &
                  'qflg size too small',size(qflg),nchanloc
             call errorHalt(1)
          end if
          call readNcdfAttr(ncid,attrName='Data Quality Flags', &
          attr=qflg(1:nchanloc))
       endif
  
       if (present(npos)) then
          if (size(npos) < nBeamloc) then
             print *,'err[RDR_io_module::openNcRDR]: ', &
                  'npos size too small',size(npos),nBeamloc
             call errorHalt(1)
          end if
          call readNcdfAttr(ncid,attrName='npos',attr=npos(1:nBeamloc))
       endif
  
       if (present(nccals)) then
          if (size(nccals) < nBeamloc) then
             print *,'err[RDR_io_module::openNcRDR]: ', &
                  'nccals size too small',size(nccals),nBeamloc
             call errorHalt(1)
          end if
          call readNcdfAttr(ncid,attrName='nccals',attr=nccals(1:nBeamloc))
       endif
       if (present(nwcals)) then
          if (size(nwcals) < nBeamloc) then
             print *,'err[RDR_io_module::openNcRDR]: ', &
                  'nwcals size too small',size(nwcals),nBeamloc
             call errorHalt(1)
          end if
          call readNcdfAttr(ncid,attrName='nwcals',attr=nwcals(1:nBeamloc))
       endif

       if (present(nthmrps))  &
            call readNcdfAttr(ncid,attrName='nthmrps', &
            attr=nthmrps)
  
       if (present(mapBeam)) then
          if (size(mapBeam) < nchanloc) then
             print *,'err[RDR_io_module::openNcRDR]: ', &
                  'mapBeam size too small',size(mapBeam),nchanloc
             call errorHalt(1)
          end if
          call readNcdfAttr(ncid,attrName='mapBeam',attr=mapBeam(1:nchanloc))
       endif
  
       if (present(frq)) then
          if (size(frq) < nchanloc) then
             print *,'err[RDR_io_module::openNcRDR]: ', &
                  'frq size too small',size(frq),nchanloc
             call errorHalt(1)
          end if
          call readNcdfAttr(ncid,attrName='mwfrequencies', &
          attr=frq(1:nchanloc))
       endif
  
       if (present(pol)) then
          if (size(pol) < nchanloc) then
             print *,'err[RDR_io_module::openNcRDR]: ', &
                  'pol size too small',size(pol),nchanloc
             call errorHalt(1)
          end if
          call readNcdfAttr(ncid,attrName='mwpolarities', &
          attr=pol(1:nchanloc))
       endif
  
       if (present(satalt))  &
            call readNcdfAttr(ncid,attrName='SatelliteAltitude', &
            attr=satalt)
  
       if (present(speed))  &
            call readNcdfAttr(ncid, &
            attrName='Satellite Forward Speed', &
            attr=speed)
  
       if (present(view))  &
            call readNcdfAttr(ncid,attrName='Viewing Direction', &
            attr=view)
  
       if (present(chanID)) then
          if (size(chanID) /= nchanloc) then
             print*, "err[RDR_io_module::openRDR]: ", &
                  "size(chanID) /= nchan:",size(chanID),nchanloc
             call errorHalt(1)
          endif
          call readNcdfAttr(ncid,attrName='lenID',attr=lenID, &
               silent=.TRUE.)
          if (lenID == 0) then  ! no such attribute in the file
             chanID(:)=repeat(' ',len(chanID(1)))
          else
             if (len(chanID(1)) /= lenID) then
                print*, "err[RDR_io_module]: ", &
                     "Length of chanID does not agree with file content."
                print*,'len(chanID(1)), lenID: ',len(chanID(1)),lenID
                call errorHalt(1)
             endif
             lenTot=(lenID+1)*nchanloc
             chanIDall=repeat(' ',maxIDlen)
             call readNcdfAttr(ncid,attrName='chanID', &
                  attr=chanIDall(1:lenTot))
             do ich = 1, nchanloc
                iid1=(ich-1)*(lenID+1)+1
                iid2=ich*(lenID+1)-1
                chanID(ich)=chanIDall(iid1:iid2)
             enddo
          endif
       endif

       if (present(thermID) .and. present(nthmrps)) then
          call readNcdfAttr(ncid,attrName='lenID',attr=lenID, &
               silent=.TRUE.)
          if (lenID == 0) then  ! no such attribute in the file
             thermID(:)=repeat(' ',len(thermID(1)))
          else
             if (len(thermID(1)) /= lenID) then
                print*, "err[RDR_io_module]: ", &
                     "Length of thermID does not agree with file content."
                print*,'len(thermID(1)), lenID: ',len(thermID(1)),lenID
                call errorHalt(1)
             endif
             lenTot=(lenID+1)*nthmTotLoc/nthmrps 
             thermIDall=repeat(' ',maxIDlen)
             call readNcdfAttr(ncid,attrName='thermID', &
                  attr=thermIDall(1:lenTot))
             do ich = 1, nthmTotLoc/nthmrps     
                iid1=(ich-1)*(lenID+1)+1
                iid2=ich*(lenID+1)-1
                thermID(ich)=thermIDall(iid1:iid2)
             enddo
          endif
       endif

    endif

    return

  end subroutine openNcRDR

!-----------------------------------------------------------------------------

  subroutine queryNcRDR(fname,frq,pol,nrec,nchan, &
    nBeam,mapBeam,npos,nccals,nwcals,nthmrps,nthmr,nthmTot, &
    nqc,qflg,nscid,instid,modver,nfilever, &
    speed,satalt,view,createTime,chanID)
    
    !---in/out variables
    character(len=*),      intent(in),    optional :: fname
    real,    dimension(:), intent(inout), optional :: frq
    integer, dimension(:), intent(inout), optional :: pol
    integer,               intent(inout), optional :: nrec
    integer,               intent(inout), optional :: nchan,nBeam
    integer, dimension(:), intent(inout), optional :: mapBeam,npos
    integer, dimension(:), intent(inout), optional :: nccals,nwcals
    integer,               intent(inout), optional :: nthmrps
    integer,               intent(inout), optional :: nthmr
    integer,               intent(inout), optional :: nthmTot
    integer,               intent(inout), optional :: nqc
    integer, dimension(:), intent(inout), optional :: qflg
    integer,               intent(inout), optional :: nscid
    integer,               intent(inout), optional :: instid,modver,nfilever
    real,                  intent(inout), optional :: speed,satalt
    integer,               intent(inout), optional :: view
    character(len=*),      intent(inout), optional :: createTime
    character(len=*),dimension(:),intent(inout),optional :: chanID
  
    !---local variables
    integer :: ncid
    
    call openNcRDR(ncid=ncid,fname=fname,frq=frq,pol=pol,nrec=nrec, &
         nchan=nchan,nBeam=nBeam,mapBeam=mapBeam,npos=npos,nccals=nccals, &
         nwcals=nwcals,nthmrps=nthmrps,nthmr=nthmr,nthmTot=nthmTot, &
         nqc=nqc,qflg=qflg,nscid=nscid, &
         instid=instid,modver=modver,nfilever=nfilever,speed=speed, &
         satalt=satalt,view=view,createTime=createTime,chanID=chanID)
    
    call closeNcRDR(ncid=ncid)
  
  end subroutine queryNcRDR

  subroutine putNcRDR(ncid,vE,vC,vW,cThm, &
    norbit,nline,lat,lon,time,eaa,eia,sza,saa,surfalt,satalt,&
    sel_sc, saz_sc, mel_sc, maz_sc, ang_mncs, qc)

    use ncdf_module

    !---in/out variables
    integer,               intent(in)              :: ncid
    real,    dimension(:), intent(in),    optional :: vE,vC,vW
    integer, dimension(:), intent(in),    optional :: cThm
    integer,               intent(in),    optional :: norbit,nline
    real,    dimension(:), intent(in),    optional :: lat,lon
    real,    dimension(:), intent(in),    optional :: eaa,eia
    real,    dimension(:), intent(in),    optional :: sza,saa
    real,    dimension(:), intent(in),    optional :: surfalt
    real,                  intent(in),    optional :: satalt
    integer, dimension(:), intent(in),    optional :: qc
    integer, dimension(:), intent(in),    optional :: time
!!! extended for geoRDRintfce
    real,                  intent(in),    optional :: sel_sc,saz_sc
    real,    dimension(:), intent(in),    optional :: mel_sc,maz_sc
    real,    dimension(:), intent(in),    optional :: ang_mncs

    if (present(vE))  &
         call writeNcdfData(ncid,vE,varname='vE', &
         varLongName='Radiometric measurement of Earth', &
         varLenName=(/'nPts','nRec'/), &
         varUnit='counts', status='append')

    if (present(vC))  &
         call writeNcdfData(ncid,vC,varname='vC', &
         varLongName='Radiometric measurement of cold target', &
         varLenName=(/'nccas','nRec '/),&
         varUnit='counts', status='append')
    
    if (present(vW))  &
         call writeNcdfData(ncid,vW,varname='vW', &
         varLongName='Radiometric measurement of warm target', &
         varLenName=(/'nwcas','nRec '/),&
         varUnit='counts', status='append')
    
    if (present(cThm)) &
         call writeNcdfData(ncid,cThm,varname='cThm', &
         varLongName='Thermistor readings', &
         varLenName=(/'nthmTot','nRec   '/), &
         varUnit='counts',status='append')
    
    if (present(time)) &
         call writeNcdfData(ncid,time,varname='date', &
         varLongName='date and time', &
         varUnit='year month day hour min msec', &
         varLenName=(/'nTime','nrec '/), &
         status='append')

    if (present(norbit))  &
         call writeNcdfData(ncid,norbit,varname='norbit', &
         varLongName='orbit number', &
         status='append')
       
    if (present(nline))  &
         call writeNcdfData(ncid,nline,varname='nline', &
         varLongName='scan line number of orbit', &
         status='append')

    if (present(lat)) then
       if (any(lat(:) /= defaultFloat)) &
            call writeNcdfData(ncid,lat,varname='lat', &
            varLongName='latitude', &
            varUnit='degrees north', &
            varLenName=(/'nBmt','nRec'/), &
            status='append')
    endif
       
    if (present(lon)) then
       if (any(lon(:) /= defaultFloat)) &
            call writeNcdfData(ncid,lon,varname='lon', &
            varLongName='longitude', &
            varUnit='degrees east', &
            varLenName=(/'nBmt','nRec'/), &
            status='append')
    endif

    if (present(eia)) &
         call writeNcdfData(ncid,eia,varname='eia', &
         varLongName='Earth Incidence Angle', &
         varLenName=(/'nBmt','nRec'/), &
         varUnit='degrees', status='append')

    if (present(eaa))  &
         call writeNcdfData(ncid,eaa,varname='eaa', &
         varLongName='Earth Azimuth Angle', &
         varLenName=(/'nBmt','nRec'/), &
         varUnit='degrees', status='append')

    if (present(saa)) then
       if (any(saa(:) /= defaultFloat)) &
            call writeNcdfData(ncid,saa,varname='saa', &
            varLongName='Solar Azimuth Angle', &
            varLenName=(/'nBmt','nRec'/), &
            varUnit='degrees', status='append')
    endif

    if (present(sza)) then
       if (any(sza(:) /= defaultFloat)) &
            call writeNcdfData(ncid,sza,varname='sza', &
            varLongName='Solar Zenith Angle', &
            varLenName=(/'nBmt','nRec'/), &
            varUnit='degrees', status='append')
    endif

    if (present(surfalt)) then
       if (any(surfalt(:) /= defaultFloat)) &
            call writeNcdfData(ncid,surfalt,varname='surfalt', &
            varLongName='Surface Altitude', &
            varLenName=(/'nBmt','nRec'/), &
            varUnit='m', status='append')
    endif

    if (present(satalt)) then
       call writeNcdfData(ncid,satalt,varname='satalt', &
            varLongName='Satellite Altitude', &
            varUnit='km', status='append')
    endif

    if (present(qc)) &
         call writeNcdfData(ncid,qc,varname='qc', &
         varLongName='Quality control flag', &
         varLenName=(/'nQC ','nRec'/), &
         status='append')

!!! extended for GeoRDRintfce
    if (present(sel_sc))  &
         call writeNcdfData(ncid,sel_sc,varname='sel_sc', &
         varLongName='Solar zenith angle wrt S/C Coordinate', &
         varUnit='degrees', status='append')

    if (present(saz_sc))  &
         call writeNcdfData(ncid,saz_sc,varname='saz_sc', &
         varLongName='Solar azimuth angle wrt S/C Coordinate', &
         varUnit='degrees', status='append')

    if (present(mel_sc))  &
         call writeNcdfData(ncid,mel_sc,varname='mel_sc', &
         varLongName='Lunar zenith angle wrt S/C Coordinate', &
         varLenName=(/'nccas','nRec '/),&
         varUnit='degrees', status='append')
    
    if (present(maz_sc))  &
         call writeNcdfData(ncid,maz_sc,varname='maz_sc', &
         varLongName='Lunar azimuth angle wrt S/C Coordinate', &
         varLenName=(/'nccas','nRec '/),&
         varUnit='degrees', status='append')

    if (present(ang_mncs))  &
         call writeNcdfData(ncid,ang_mncs,varname='ang_mncs', &
         varLongName='Lunar viewing angle wrt cold-sky bore sight', &
         varLenName=(/'nccas','nRec '/),&
         varUnit='degrees', status='append')

    return

  end subroutine putNcRDR

!-----------------------------------------------------------------------------

  subroutine getNcRDR(ncid,irec,vE,vC,vW,cThm, & 
     norbit,nline,lat,lon,time,eaa,eia,sza,saa,surfalt,satalt, &
     sel_sc, saz_sc, mel_sc, maz_sc, ang_mncs, qc)

    use ncdf_module
                                                                                                      
    !---in/out variables
    integer,               intent(in)              :: ncid
    integer , intent(in), optional                 :: irec
    real,    dimension(:), intent(inout),    optional :: vE,vC,vW
    integer, dimension(:), intent(inout),    optional :: cThm
    integer,               intent(inout),    optional :: norbit,nline
    real,    dimension(:), intent(inout),    optional :: lat,lon
    integer, dimension(:), intent(inout),    optional :: time
    real,    dimension(:), intent(inout),    optional :: eaa,eia
    real,    dimension(:), intent(inout),    optional :: sza,saa
    real,    dimension(:), intent(inout),    optional :: surfalt
    real,                  intent(inout),    optional :: satalt
    integer, dimension(:), intent(inout),    optional :: qc
!!! extended for geoRDRintfce
    real,                  intent(inout),    optional :: sel_sc,saz_sc
    real,    dimension(:), intent(inout),    optional :: mel_sc,maz_sc
    real,    dimension(:), intent(inout),    optional :: ang_mncs

      if (present(vE))  &
          call readNcdfData(ncid,vE,varname='vE', &
          record_no=irec,status='single_record')

      if (present(vC))  &
          call readNcdfData(ncid,vC,varname='vC', &
          record_no=irec,status='single_record')

      if (present(vW))  &
          call readNcdfData(ncid,vW,varname='vW', &
          record_no=irec,status='single_record')

      if (present(cThm)) &
           call readNcdfData(ncid,cThm,varname='cThm', &
           record_no=irec,status='single_record')

      if (present(norbit))  &
          call readNcdfData(ncid,norbit,varname='norbit', &
          record_no=irec,status='single_record')
                                                                                              
      if (present(nline))  &
          call readNcdfData(ncid,nline,varname='nline', &
          record_no=irec,status='single_record')
                                                                                              
      if (present(lat))  &
          call readNcdfData(ncid,lat,varname='lat', &
          record_no=irec,status='single_record')
                                                                                              
      if (present(lon))  &
          call readNcdfData(ncid,lon,varname='lon', &
          record_no=irec,status='single_record')

      if (present(time))  &
          call readNcdfData(ncid,time,varname='date', &
          record_no=irec,status='single_record')
                                                                                              
      if (present(eaa))  &
          call readNcdfData(ncid,eaa,varname='eaa', &
          record_no=irec,status='single_record')
                                                                                              
      if (present(eia))  &
          call readNcdfData(ncid,eia,varname='eia', &
          record_no=irec,status='single_record')
                                                                                              
      if (present(sza))  &
          call readNcdfData(ncid,sza,varname='sza', &
          record_no=irec,status='single_record')

      if (present(saa))  &
          call readNcdfData(ncid,saa,varname='saa', &
          record_no=irec,status='single_record')
                                                                                              
      if (present(surfalt))  &
          call readNcdfData(ncid,surfalt,varname='surfalt', &
          record_no=irec,status='single_record')
                                                                                              
      if (present(satalt))  &
          call readNcdfData(ncid,satalt,varname='satalt', &
          record_no=irec,status='single_record')
                                                                                              
      if (present(qc))  &
          call readNcdfData(ncid,qc,varname='qc', &
          record_no=irec,status='single_record')

!!! extended for GeoRDRintfce
      if (present(sel_sc))  &
          call readNcdfData(ncid,sel_sc,varname='sel_sc', &
          record_no=irec,status='single_record')

      if (present(saz_sc))  &
          call readNcdfData(ncid,saz_sc,varname='saz_sc', &
          record_no=irec,status='single_record')

      if (present(mel_sc))  &
          call readNcdfData(ncid,mel_sc,varname='mel_sc', &
          record_no=irec,status='single_record')

      if (present(maz_sc))  &
          call readNcdfData(ncid,maz_sc,varname='maz_sc', &
          record_no=irec,status='single_record')

      if (present(ang_mncs))  &
          call readNcdfData(ncid,ang_mncs,varname='ang_mncs', &
          record_no=irec,status='single_record')
  
    return

  end subroutine getNcRDR

!------------------------------------------------------------------------------

  subroutine closeNcRDR(ncid,msg)

    use ncdf_module

    integer,           intent(in)           :: ncid
    character (len=*), intent(in), optional :: msg
    
    call closeNcdfFile(ncid,msg)

    return
      
  end subroutine closeNcRDR

end module RDR_io_module
