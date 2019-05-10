c "$Id$"
      module rad_io_module
c**********************************************************************
c* Module: rad_io_module
c* Purpose: Radiance/Brightness temperature I/O
c* Description: The module is used to open/close, read/write radiance
c*              or brightness temperature to the netCDF file.
c*            
c* Inputs/Outputs:
c* Var_name      Type       Description
c* --------      ----       -----------
c* ncid          integer    NetCDF file ID (or descriptor)
c* fname         character  Filename for output.
c* xid           character  Profile ID. (e.g., IPO11111)
c* norbit        integer    Satellite orbit number
c* nline         integer    Scan line of orbit
c* nelement      integer    Element in scan line
c* satalt        real       Satellite altitude
c* surfalt       real       Surface altitude
c* lat           real       FOV latitude
c* lon           real       FOV longitude                        
c* time          integer    Year-month-day-hour-minute-msecond
C* eaa           real       Earth azimuth angle per channel
C* eia           real       Earth incidence angle per channel
C* pra           real       polarization rotation angle due to attitude variation
C* sza           real       Solar Zenith Angle
c* tb (like tbsim)     real Brightness temperatures
c* pol           integer    channel polarities
c* frq           real       channel frequencies (ghz)
c* nchan       integer    number of microwave sensor channels
c* nscid         integer    Spacecraft ID
c* instid        integer    Sensor ID
c* modver        integer    Sensor model configuration number
c* nfilever      integer    Sensor constants file version
c* qflg          integer    Data quality flag
c* refaz         integer    Number of reference sample azimuths
c* eianom        float      Nominal Earth incidence angle (/channel)
c* nadnom        float      Nominal Nadir angle (npos x channel)
c* aznom         float      Nominal Azimuth angle (npos x channel)
c* npos          integer    Number of Pixels   
c* speed         float      Satellite foward speed (nominal)
c* view          integer    Satellite viewing direction (forward/aft 1/0)
c* FOVtype       integer    FOV organization flag (CFOV=0,EFOV=1,CFOVvarEIA=2)
c* chanID        character  string that uniquely identifies each channel
c*
c* Externals: ncdf_module
c* Copyright: 2003, AER, Inc.  
c* Developed by AER, Inc., 2003
c**********************************************************************
      use ncdf_module
      implicit none
      private
      public openRad,closeRad,getRad,putRad,getRad1Chan,putRad1Chan,
     $     queryRad,replaceRad1Chan,replaceRad,appendRad1Chan,appendRad

      integer , save :: irec=0
      logical , dimension (MAXFILE) :: newForm=.false.
      real , parameter :: defaultFloat=-99999.
      integer, dimension(MAXFILE) :: fid=-1, FovTypeLoc=-1
      integer :: fidIndex=0     !counter for # of opened files
      integer :: nBeamloc
      integer, dimension(:), pointer :: mapBeamloc ! to argument not present
      integer, dimension(:), allocatable :: nposloc
      interface openRad
      module procedure openNcRad
      end interface

      interface closeRad
      module procedure closeNcRad
      end interface

      interface queryRad
      module procedure queryNcRad
      end interface

      CONTAINS
      
      subroutine openNcRad(ncid,fname,satalt,pol,frq,nchan,nrec,
     $     nscid,instid,modver,nfilever,qflg,refaz,eianom,
     &     npos,nqc,speed,view,FOVtype,createTime,xidLen,
     $     nBeam,mapBeam,npts,nele,nang,nInner,status,
     $     nadnom,aznom,wvn,chanID)
      
      use MapBeamDefaults

      !---in/out variables
      integer , intent(inout) :: ncid
      character (len=*), intent (in), optional :: fname
      character (len=*), intent (in), optional :: status
      character (len=*), intent (inout), optional :: createTime
      integer , intent (inout), optional :: nchan, nrec
      integer , intent (inout), optional :: nqc,nBeam,nPts,nEle,nAng
      integer , intent (inout), optional :: nInner
      integer, dimension(:), intent(inout), optional :: npos, mapBeam
      integer , intent(inout), optional :: xidLen
      real , intent (inout), optional :: speed,satalt
      integer , intent (inout), optional :: view,FOVtype,nscid
      integer , intent (inout), optional :: instid,modver,nfilever
      integer , dimension (:), intent(inout), optional :: refaz,qflg
      real , dimension (:), intent (inout), optional :: frq,eianom
      real , dimension (:), intent (inout), optional :: wvn
      real , dimension (:,:), intent (inout), optional :: nadnom,aznom
      integer , dimension (:), intent (inout), optional :: pol
      character(len=*), dimension(:), intent (inout), optional :: chanID

      !---Local variables
      integer, dimension(8) :: sysvalues
      character (len= 8) :: sysdate
      character (len=10) :: systime
      character (len= 5) :: syszone
      character (len=20) :: myStatus
      integer :: nPtsLoc, nEleLoc, index, iBeam, nAngLoc, i
      integer :: nchanloc, ich, nrecLoc
      integer :: FovTypeTmp
      integer :: lenID,lenTot,iid1,iid2
      integer, parameter :: maxIDlen=4000
      character (len=maxIDlen) :: chanIDall

      if (.not.present(fname)) then
         print *,'err[rad_io_module::openNcRad]:  missing filename.'
         call errorHalt(1)
      end if

      if (present(status)) then
         myStatus = status
      else
         myStatus='old'
      endif
      if (trim(myStatus) == 'new' .or. trim(myStatus) == 'NEW' .or.
     $     trim(myStatus) == 'New') then

         if (.not. present(nchan)) then
            print*, "err[rad_io_module::openRad]: nchan not present."
            call errorHalt(1)
         endif
            
         call openNcdfFile(ncid, fname, status='new',
     &        unlimited_dim_name='nRec')

         call date_and_time(sysdate, systime, syszone, sysvalues)

         call writeNcdfAttr(ncid,attrName='Creationtime',
     &        attr=trim(sysdate//systime))

         if (present(nscid)) 
     &        call writeNcdfAttr(ncid,attrName='SpacecraftID',
     &        attr=(/nscid/))

         if (present(instid)) 
     &        call writeNcdfAttr(ncid,attrName='SensorID',  
     &        attr=(/instid/))

         if (present(modver)) 
     &        call writeNcdfAttr(ncid,attrName='SensorModelVersion',
     &        attr=(/modver/))

         if (present(nfilever)) 
     &        call writeNcdfAttr(ncid,attrName='ConstantsFileVersion',
     &        attr=(/nfilever/))

         if (present(qflg)) 
     &        call writeNcdfAttr(ncid,attrName='DataQualityFlags',
     &        attr=int(qflg(1:nchan)))

         if (.not. present(FOVtype)) then
            call initFOVtype(ncid,0)
         else
            call initFOVtype(ncid,FOVtype)
         endif

         index = searchID(ncid)
         call mapBeamDefNew(nchan,nBeam,mapBeam,nBeamloc,mapBeamloc)
         select case (FovTypeLoc(index))
         case (0)               !colocated with invariable eia with channels
            if (nBeamLoc /= 1) then
               print*,"err[rad_io_module::openNcRad]: ",
     $              "nBeamLoc should be 1 in the colocated mode"
               call errorHalt(1)
            endif

            if (.not. allocated(nposloc)) allocate(nposloc(nBeamloc))

            nposloc(:) = 1
            if (present(npos)) then 
               if (any(npos(:) /=1)) then
                  print*,"warning[rad_io_module]::openNcRad]: ",
     $                 "npos should be absent or a vector of 1's in ",
     $                 "the colocated mode"
               endif
            endif

            nPtsLoc=nchan
            nEleLoc=1
            nAngLoc=1
         case (1)               !un-colocated structure
            if (.not. allocated(nposloc)) allocate(nposloc(nBeamloc))

            if (.not. present(npos)) then
               print *,'err[rad_io_module::openNcRad]: ',
     &              ' missing npos for un-colocated mode.'
               call errorHalt(1)
            else
               nposloc = npos
            endif

            nPtsLoc = 0
            do ich = 1, nchan
               nPtsLoc=nPtsLoc+nposloc(mapBeamloc(ich))
            enddo
            nEleLoc=sum(nposloc(1:nBeamloc))
            nAngLoc=sum(nposloc(1:nBeamloc))
         case (2)               !colocated with varying eia's by channels
            if (present(npos)) then 
               if (any(npos(:) /=1)) then
                  print*,"warning[rad_io_module]::openNcRad]: ",
     $                 "npos should be absent or a vector of 1's in ",
     $                 "the colocated mode."
               endif
            endif

            if (.not. allocated(nposloc)) allocate(nposloc(nBeamloc))

            nposloc(:)=1
            nPtsLoc=nchan
            nEleLoc=1
            nAngLoc=nBeamloc
         case default
            print*, "err[rad_io_module::openNcRad]: ", "
     $           FOVtype exceeded the range (0-2)"
            call errorHalt(1)
         end select

         Call writeNcdfDim(ncid,dimName='nPts',dimLen=nPtsLoc)
         call writeNcdfDim(ncid,dimName='nEle',dimLen=nEleLoc)
         call writeNcdfDim(ncid,dimName='nAng',dimLen=nAngLoc)
         call writeNcdfAttr(ncid,attrName='FovOrganizationFlag',
     &        attr=(/FovTypeLoc(index)/))
         call writeNcdfAttr(ncid,attrName='nchan',attr=(/nchan/))

         if (present(nInner)) then
            if (FovTypeLoc(index) /= 1) then
               call writeNcdfAttr(ncid,attrName='nInner',
     $              attr=(/nInner/))
            else
               print*, 'warning[rad_io_module::openNcRad]: ',
     $              'nInner not allowed ',
     $              'for output when FOVtype 1. Use npos, instead'
            endif
         endif

         call writeNcdfAttr(ncid,attrName='nBeam',attr=nBeamloc)
         call writeNcdfAttr(ncid,attrName='mapBeam',attr=mapBeamloc)
         call writeNcdfAttr(ncid,attrName='npos',attr=nposloc)

         if (present(frq)) 
     &        call writeNcdfAttr(ncid,attrName='mwfrequencies',
     &        attr=frq(1:nchan))

         if (present(wvn)) 
     &        call writeNcdfAttr(ncid,attrName='wavenumbers',
     &        attr=wvn(1:nchan))

         if (present(pol)) 
     &        call writeNcdfAttr(ncid,attrName='mwpolarities',
     &        attr=pol(1:nchan))

         if (present(refaz)) 
     &        call writeNcdfAttr(ncid,attrName='RefChannelAzimuth',
     &        attr=int(refaz(1:nchan)))

         if (present(satalt)) 
     &        call writeNcdfAttr(ncid,attrName='SatelliteAltitude',
     &        attr=(/satalt/))

         if (present(eianom)) 
     &        call writeNcdfAttr(ncid,attrName='NomEarthIncAngl',
     &        attr=eianom(1:nchan))

         if (present(speed)) 
     &        call writeNcdfAttr(ncid,
     $        attrName='SatelliteForwardSpeed',
     &        attr=(/speed/))

         if (present(view)) 
     &        call writeNcdfAttr(ncid,attrName='ViewingDirection',
     &        attr=(/view/))

         if (present(nadnom))
     $        call writeNcdfData(ncid,nadnom,
     &        varname='nadnom',varUnit='degrees',
     &        varLongname='Nominal Nadir Angle')
            
         if (present(aznom))
     $        call writeNcdfData(ncid,aznom,
     &        varname='aznom',varUnit='degree',
     &        varLongname='Nominal Scan Azimuth Angle')
            
         if (present(chanID)) then
            if (size(chanID) /= nchan) then
               print*, "err[rad_io_module::openRad]: ",
     &         "size(chanID) /= nchan:",size(chanID),nchan
               call errorHalt(1)
            endif
            lenID=len(chanID(1))
            lenTot=(lenID+1)*nchan
            if (lenTot > maxIDlen) then
               print *, "err[rad_io_module::openRad]: ",
     &         "Total length of chanID exceeds length of chanIDall"
               print *,'lenTot,maxIDlen:',lenTot,maxIDlen
               call errorHalt(1)
            endif
            call writeNcdfAttr(ncid,attrName='lenID',attr=lenID)
            chanIDall=repeat(' ',maxIDlen)
            do ich = 1, nchan
               iid1=(ich-1)*(lenID+1)+1
               iid2=ich*(lenID+1)-1
               chanIDall(iid1:iid2)=chanID(ich)
            enddo
            call writeNcdfAttr(ncid,attrName='chanID',
     &        attr=chanIDall(1:lenTot))
         endif

      else                      !read existing file
         if (myStatus == 'replace' .or. myStatus == 'append') then
            call openNcdfFile(ncid, trim(fname), status='replace')
         else
            call openNcdfFile(ncid, trim(fname))
         endif
         nrecLoc=readNcdfDim(ncid,DimName='nRec',silent=.TRUE.)
         if (nrecLoc /= 0) then
            newForm(getNcidIndex(ncid))=.true.
         else                   !-- backward compatibility
            nrecLoc=readNcdfDim(ncid,DimName='nProfiles',silent=.TRUE.)
         endif
         
         if (present(nrec)) nrec=nrecLoc

         if (present(createTime))
     $        call readNcdfAttr(ncid,attrName='Creationtime',
     &        attr=createTime)
      
         if (present(xidLen)) xidLen=readNcdfDim(ncid,'nId')

         if (present(nqc)) nqc=readNcdfDim(ncid,'nQC',silent=.true.)

         call readNcdfAttr(ncid,attrName='nchan',attr=nchanloc)

         if (nchanloc == 0) then
            print*, "err[rad_io_module::openRad]: nchan equals zero."
            call errorHalt(1)
         endif

         if (present(nchan)) nchan=nchanloc

         if (present(nscid)) 
     &        call readNcdfAttr(ncid,attrName='Spacecraft ID',
     &        attr=nscid)
      
         if (present(instid)) 
     &        call readNcdfAttr(ncid,attrName='Sensor ID',
     &        attr=instid)
 
         if (present(modver)) 
     &        call readNcdfAttr(ncid,attrName='Sensor Model Version',
     &        attr=modver)

         if (present(nfilever)) 
     &        call readNcdfAttr(ncid,attrName='Constants File Version',
     &        attr=nfilever)

         if (present(qflg)) then 
            if (size(qflg) < nchanloc) then
               print*, "err[rad_io_module::openRad]: "//
     $              "length of vector qflg too small."
               call errorHalt(1)
            endif
            call readNcdfAttr(ncid,attrName='Data Quality Flags',
     &           attr=qflg(1:nchanloc))
         endif

         call readNcdfAttr(ncid,attrName='FovOrganizationFlag',
     &        attr=FovTypeTmp)
         call initFOVtype(ncid,FovTypeTmp)
         if (present(FOVtype)) FOVtype=FovTypeTmp

         nEleLoc = readNcdfDim(ncid,DimName='nEle')
         if (present(nEle)) nEle=nEleLoc

         call readNcdfAttr(ncid,attrName='nBeam',attr=nBeamLoc,
     $        silent=.true.)

         if (present(nInner)) then
            if (FovTypeTmp /= 1) then
               call readNcdfAttr(ncid,attrName='nInner',
     $              attr=nInner)
            else
               print*,'warning[openRad]: nInner not ',
     $              'allowed in input for FOVtype 1, ',
     $              'Use npos, instead.'
            endif
         endif

         nPtsloc = readNcdfDim(ncid,'nPts',silent=.TRUE.)
         if (present(nPts)) nPts=nPtsLoc

         allocate(mapBeamloc(nchanloc))

         if (nBeamloc == 0) then
            if (FovTypeTmp == 0) then
               nBeamloc = 1
               mapBeamloc(1:nchanloc) = 1
               if (.not. allocated(nposloc)) allocate(nposloc(nBeamloc))
               nposloc(:) = 1
            else
               nBeamloc = nchanloc
               mapBeamloc(1:nchanloc) = (/(i,i=1,nchanloc)/)
               if (.not. allocated(nposloc)) allocate(nposloc(nBeamloc))
               nposloc(:) = nEleLoc
            endif
            nAngLoc=sum(nposloc)
         else
            call readNcdfAttr(ncid,attrName='mapBeam',attr=mapBeamloc)
            if (.not. allocated(nposloc)) allocate(nposloc(nBeamloc))
            call readNcdfAttr(ncid,attrName='npos',attr=nposloc)
            nAngLoc = readNcdfDim(ncid,'nAng',silent=.TRUE.)
            if (nAngLoc == 0) nAngLoc=sum(nposloc)
         endif

         if (present(nBeam)) nBeam = nBeamloc
         if (present(mapBeam)) mapBeam = mapBeamloc
         if (present(npos)) npos = nposloc
         if (present(nAng)) nAng = nAngLoc

         if (present(frq)) then
            if (size(frq) < nchanloc) then
               print*, "err[rad_io_module::openRad]: "//
     $              "length of vector frq too small."
               call errorHalt(1)
            endif
            call readNcdfAttr(ncid,attrName='mwfrequencies',
     &           attr=frq(1:nchanloc))
         endif

         if (present(wvn)) then
            if (size(wvn) < nchanloc) then
               print*, "err[rad_io_module::openRad]: "//
     $              "length of vector wvn too small."
               call errorHalt(1)
            endif
            call readNcdfAttr(ncid,attrName='wavenumbers',
     &           attr=wvn(1:nchanloc))
         endif

         if (present(pol)) then
            if (size(pol) < nchanloc) then
               print*, "err[rad_io_module::openRad]: "//
     $              "length of vector pol too small."
               call errorHalt(1)
            endif
            call readNcdfAttr(ncid,attrName='mwpolarities',
     &           attr=pol(1:nchanloc))
         endif

         if (present(refaz)) then
            if (size(refaz) < nchanloc) then
               print*, "err[rad_io_module::openRad]: "//
     $              "length of vector refaz too small."
               call errorHalt(1)
            endif
            call readNcdfAttr(ncid,attrName='Ref Channel azimuth',
     &           attr=refaz(1:nchanloc))
         endif

         if (present(satalt)) 
     &        call readNcdfAttr(ncid,attrName='SatelliteAltitude',
     &        attr=satalt)

         if (present(eianom)) then
            if (size(eianom) < nchanloc) then
               print*, "err[rad_io_module::openRad]: "//
     $              "length of vector eianom too small."
               call errorHalt(1)
            endif
            call readNcdfAttr(ncid,attrName='Nom. Earth Inc Angl',
     &           attr=eianom(1:nchanloc))
         endif

         if (present(speed)) 
     &        call readNcdfAttr(ncid,
     $        attrName='Satellite Forward Speed',
     &        attr=speed)

         if (present(view)) 
     &        call readNcdfAttr(ncid,attrName='Viewing Direction',
     &        attr=view)

         if (present(nadnom))
     $        call readNcdfData(ncid,nadnom,varname='nadnom')
         
         if (present(aznom))
     $        call readNcdfData(ncid,aznom,varname='aznom')

         if (present(chanID)) then
            if (size(chanID) /= nchanloc) then
               print*, "err[rad_io_module::openRad]: ",
     &         "size(chanID) /= nchan:",size(chanID),nchanloc
               call errorHalt(1)
            endif
            call readNcdfAttr(ncid,attrName='lenID',attr=lenID,
     &         silent=.TRUE.)
            if (lenID == 0) then   ! No such attribute in file
               chanID(:)=repeat(' ',len(chanID(1)))
            else
               if (len(chanID(1)) /= lenID) then
                  print *, "err[rad_io_module::openRad]: ",
     &            "Length of chanID does not agree with file content."
                  print *,'len(chanID(1)),lenID:',len(chanID(1)),lenID
                  call errorHalt(1)
               endif
               lenTot=(lenID+1)*nchanloc
               chanIDall=repeat(' ',maxIDlen)
               call readNcdfAttr(ncid,attrName='chanID',
     &            attr=chanIDall(1:lenTot))
               do ich = 1, nchanloc
                  iid1=(ich-1)*(lenID+1)+1
                  iid2=ich*(lenID+1)-1
                  chanID(ich)=chanIDall(iid1:iid2)
               enddo
            endif
         endif

      endif

      deallocate(nposloc)
      deallocate(mapBeamLoc)

      return
      end subroutine openNcRad
c--        
      subroutine queryNcRad(fname,satalt,pol,frq,nchan,nrec,nqc,
     $     nscid,instid,modver,nfilever,qflg,refaz,eianom,
     &     npos,speed,view,FOVtype,createTime,xidLen,nBeam,
     $     mapBeam,npts,nele,nang,nInner,nadnom,aznom,wvn,chanID)

      !---in/out variables
      character (len=*), intent (in), optional :: fname
      character (len=*), intent (inout), optional :: createTime
      integer , intent (inout), optional :: nchan, nrec
      integer , intent (inout), optional :: nqc,nBeam,nPts,nEle,nAng
      integer , intent (inout), optional :: nInner
      integer , dimension(:), intent(inout), optional :: npos, mapBeam
      integer , intent(inout), optional :: xidLen
      real , intent (inout), optional :: speed,satalt
      integer , intent (inout), optional :: view,FOVtype,nscid
      integer , intent (inout), optional :: instid,modver,nfilever
      integer , dimension (:), intent(inout), optional :: refaz,qflg
      real , dimension (:), intent (inout), optional :: frq,eianom
      real , dimension (:), intent (inout), optional :: wvn
      integer , dimension (:), intent (inout), optional :: pol
      real , dimension (:,:), intent (inout), optional :: nadnom,aznom
      character(len=*), dimension(:), intent (inout), optional :: chanID

c-- Local variables
      integer :: ncid
      
      call openNcRad(ncid=ncid,fname=fname,satalt=satalt,pol=pol,
     $     frq=frq,nchan=nchan,nrec=nrec,nscid=nscid,nqc=nqc,
     $     instid=instid,modver=modver,nfilever=nfilever,qflg=qflg,
     $     refaz=refaz,eianom=eianom,
     &     npos=npos,speed=speed,view=view,FOVtype=FOVtype,
     $     createTime=createTime,xidLen=xidLen,nBeam=nBeam,
     $     mapBeam=mapBeam,npts=npts,nele=nele,nang=nang,
     $     nInner=nInner,nadnom=nadnom,aznom=aznom,wvn=wvn,
     $     chanID=chanID)
      call closeNcRad(ncid=ncid)
      end subroutine queryNcRad
c--
      subroutine putRad1Chan(ncid,xid,norbit,nline,nelement,
     $     lat,lon,time,eaa,eia,pra,tb,sza,saa,surfalt,satalt,ichan,qc,
     $     rad,radMW)

      !---in/out variables
      integer , intent(in) :: ncid
      character (len=*), intent (in), optional :: xid
      integer , intent (in), optional :: norbit,nline,nelement,ichan
      real , dimension (:), intent (in), optional :: lat,lon
      real , dimension (:), intent (in), optional :: surfalt
      real , intent (in), optional :: satalt
      real , dimension (:), intent (in), optional :: sza,saa
      real , dimension (:), intent (in), optional :: eaa,eia,tb
      real , dimension (:), intent (in), optional :: pra
      real , dimension (:), intent (in), optional :: rad
      real , dimension (:), intent (in), optional :: radMW
      integer , dimension (:), intent (in), optional :: time, qc
      integer :: index

      if (present(tb)) 
     &     call writeNcdfData(ncid,tb,varname='tb',
     &     varLongName='Brightness Temperatures',
     $     varLenName=(/'nPts','nRec'/),
     &     varUnit='Kelvin', status='append')
         
      if (present(rad)) 
     &     call writeNcdfData(ncid,rad,varname='rad',
     &     varLongName='Radiances',
     $     varLenName=(/'nPts','nRec'/),
     &     varUnit='mW/m2/str/cm-1', status='append')
         
      if (present(radMW)) 
     &     call writeNcdfData(ncid,radMW,varname='rad',
     &     varLongName='Radiances',
     &     varLenName=(/'nPts','nRec'/),
     &     varUnit='mW/m2/str/GHz', status='append')
         
      if (present(time))
     $     call writeNcdfData(ncid,time,varname='date',
     &     varLongName='date and time',
     &     varUnit='year month day hour min msec',
     $     varLenName=(/'nTime','nrec '/),
     &     status='append')

      if (present(xid)) 
     &     call writeNcdfData(ncid,xid,varname='id',
     &     varLongName='Sample identification string',
     $     strLenName='nId', status='append')

      if (present(norbit)) 
     &     call writeNcdfData(ncid,norbit,varname='norbit',
     &     varLongName='orbit number',
     &     status='append')
         
      if (present(nline)) 
     &     call writeNcdfData(ncid,nline,varname='nline',
     &     varLongName='scan line number of orbit',
     &     status='append')

      if (present(nelement)) 
     &     call writeNcdfData(ncid,nelement,varname='nelement',
     &     varLongName='element of scan line',
     &     status='append')

      if (present(ichan)) then
         index= searchID(ncid) 
         if (FovTypeLoc(index) == 0) then
            print *,'err[rad_io_module::putRad1Chan]: ',
     &           ' unmatched FOV type'
            call errorHalt(1)
         end if

         call writeNcdfData(ncid,ichan,varname='ichan',
     &        varLongName='channel numbers',
     &        status='append')
      endif

      if (present(lat)) then
         if (any(lat(:) /= defaultFloat))
     $        call writeNcdfData(ncid,lat,varname='lat',
     &        varLongName='latitude',
     &        varUnit='degrees north',
     $        varLenName=(/'nEle','nRec'/),
     &        status='append')
      endif
         
      if (present(lon)) then
         if (any(lon(:) /= defaultFloat))
     $        call writeNcdfData(ncid,lon,varname='lon',
     &        varLongName='longitude',
     &        varUnit='degrees east',
     $        varLenName=(/'nEle','nRec'/),
     &        status='append')
      endif

      if (present(eia))
     $     call writeNcdfData(ncid,eia,varname='eia',
     &     varLongName='Earth Incidence Angle',
     $     varLenName=(/'nAng','nRec'/),
     &     varUnit='degrees', status='append')

      if (present(eaa)) 
     &     call writeNcdfData(ncid,eaa,varname='eaa',
     &     varLongName='Earth Azimuth Angle',
     $     varLenName=(/'nAng','nRec'/),
     &     varUnit='degrees', status='append')
      
      if (present(pra))
     $     call writeNcdfData(ncid,pra,varname='pra',
     &     varLongName=
     $     'Polarization Rotation Angle due to Attitude variation',
     $     varLenName=(/'nAng','nRec'/),
     &     varUnit='degrees', status='append')

      if (present(saa)) then
         if (any(saa(:) /= defaultFloat))
     &        call writeNcdfData(ncid,saa,varname='saa',
     &        varLongName='Solar Azimuth Angle',
     $        varLenName=(/'nEle','nRec'/),
     &        varUnit='degrees', status='append')
      endif

      if (present(sza)) then
         if (any(sza(:) /= defaultFloat))
     &        call writeNcdfData(ncid,sza,varname='sza',
     &        varLongName='Solar Zenith Angle',
     $        varLenName=(/'nEle','nRec'/),
     &        varUnit='degrees', status='append')
      endif

      if (present(surfalt)) then
         if (any(surfalt(:) /= defaultFloat))
     &        call writeNcdfData(ncid,surfalt,varname='surfalt',
     &        varLongName='Surface Altitude',
     $        varLenName=(/'nEle','nRec'/),
     &        varUnit='m', status='append')
      endif

      if (present(satalt)) then
         call writeNcdfData(ncid,satalt,varname='satalt',
     &        varLongName='Satellite Altitude',
     &        varUnit='km', status='append')
      endif

      if (present(qc))
     $     call writeNcdfData(ncid,qc,varname='qc',
     &     varLongName='Quality control flag',
     $     varLenName=(/'nQC ','nRec'/),
     &     status='append')

      return
      end subroutine putRad1Chan
c--
      subroutine getRad1Chan(ncid,irec,xid,norbit,nline,nelement,
     &     lat,lon,time,eaa,eia,pra,tb,sza,saa,surfalt,satalt,ichan,qc,
     $     rad)

      !---in/out variables
      integer , intent(in) :: ncid
      integer , intent(in), optional :: irec
      character (len=*), intent(inout), optional :: xid
      real , dimension (:), intent(inout), optional :: lat,lon
      real , dimension (:), intent(inout), optional :: surfalt
      real , intent(inout), optional :: satalt
      real , dimension (:), intent(inout), optional :: sza,saa
      real , dimension (:), intent(inout), optional :: tb, eia, eaa
      real , dimension (:), intent(inout), optional :: pra
      real , dimension (:), intent(inout), optional :: rad
      integer , intent (inout), optional :: norbit,nline,nelement,ichan
      integer , dimension (:), intent (inout), optional :: time, qc
      integer :: index

      if (present(tb)) 
     &     call readNcdfData(ncid,tb,varname='tb',
     &     record_no=irec,status='single_record')

      if (present(rad)) 
     &     call readNcdfData(ncid,rad,varname='rad',
     &     record_no=irec,status='single_record')

      if (present(time)) 
     &     call readNcdfData(ncid,time,varname='date',
     &     record_no=irec,status='single_record')

      if (present(xid)) 
     &     call readNcdfData(ncid,xid,varname='id',
     &     record_no=irec,status='single_record')

      if (present(norbit)) 
     &     call readNcdfData(ncid,norbit,varname='norbit',
     &     record_no=irec,status='single_record')

      if (present(nline)) 
     &     call readNcdfData(ncid,nline,varname='nline',
     &     record_no=irec,status='single_record')

      if (present(nelement)) 
     &     call readNcdfData(ncid,nelement,varname='nelement',
     &     record_no=irec,status='single_record')

      if (present(ichan)) then
         index = searchID(ncid)
         if (FovTypeLoc(index) == 0) then
            print *,'err[rad_io_module::getRad1Chan]: ',
     &           ' unmatched FOV type'
            call errorHalt(1)
         end if

         call readNcdfData(ncid,nelement,varname='ichan',
     &        record_no=irec,status='single_record')
      endif

      if (present(lat)) then
         if (any(lat(:) /= defaultFloat))
     $        call readNcdfData(ncid,lat,varname='lat',
     &        record_no=irec,status='single_record')
      endif

      if (present(lon)) then
         if (any(lon(:) /= defaultFloat))
     $        call readNcdfData(ncid,lon,varname='lon',
     &        record_no=irec,status='single_record')
      endif

      if (present(eia)) then
         call readNcdfData(ncid,eia,varname='eia',
     &        record_no=irec,status='single_record')
      endif

      if (present(eaa)) then
         call readNcdfData(ncid,EAA,varname='eaa',
     &        record_no=irec,status='single_record')
      endif

      if (present(pra)) then
         call readNcdfData(ncid,pra,varname='pra',
     &        record_no=irec,status='single_record')
      endif

      if (present(saa)) then
         if (any(saa(:) /= defaultFloat))
     &        call readNcdfData(ncid,saa,varname='saa',
     &        record_no=irec,status='single_record')
      endif

      if (present(sza)) then
         if (any(sza(:) /= defaultFloat))
     &        call readNcdfData(ncid,sza,varname='sza',
     &        record_no=irec,status='single_record')
      endif
      
      if (present(surfalt)) then
         if (any(surfalt(:) /= defaultFloat))
     &        call readNcdfData(ncid,surfalt,varname='surfalt',
     &        record_no=irec,status='single_record')
      endif

      if (present(satalt)) then
         call readNcdfData(ncid,satalt,varname='satalt',
     &        record_no=irec,status='single_record')
      endif

      if (present(qc)) 
     &     call readNcdfData(ncid,qc,varname='qc',
     &     record_no=irec,status='single_record')
      
      return
      end subroutine getRad1Chan
c--
      subroutine closeNcRad(ncid,msg)
      integer, intent(in) :: ncid
      character (len=*), intent(in), optional :: msg

      call closeNcdfFile(ncid,msg)
      
      end subroutine closeNcRad
c--
      subroutine putRad(ncid,xid,norbit,nline,nelement,
     $     lat,lon,time,eaa,eia,pra,tb,sza,saa,surfalt,satalt,qc,
     $     rad,radMW)

      !---in/out variables
      integer , intent(in) :: ncid
      character (len=*), intent (in), optional :: xid
      integer , intent (in), optional :: norbit,nline,nelement
      real , intent (in), optional :: lat,lon
      real , intent (in), optional :: surfalt
      real , intent (in), optional :: satalt
      real , intent (in), optional :: sza,saa
      real , dimension(:), intent (in), optional :: tb
      real , dimension(:), intent (in), optional :: rad
      real , dimension(:), intent (in), optional :: radMW
      real , dimension(:), intent (in), optional :: eaa,eia
      real , dimension(:), intent (in), optional :: pra
      integer , dimension(:), intent (in), optional :: time,qc
      !-- local arguments
      real , dimension (1) :: myLat,myLon,mySurfalt
      real , dimension (1) :: mySza,mySaa
      
      if (present(lat)) then 
         mylat(1) = lat
      else
         mylat(1) = defaultFloat
      endif
      if (present(lon)) then
         mylon(1) = lon
      else
         mylon(1) = defaultFloat
      endif
      if (present(sza)) then
         mySza(1) = sza
      else
         mySza(1) = defaultFloat
      endif

      if (present(saa)) then
         mySaa(1) = saa
      else
         mySaa(1) = defaultFloat
      endif
      if (present(surfalt)) then 
         mySurfalt(1) = surfalt
      else
         mySurfalt(1) = defaultFloat
      endif

      call putRad1Chan(ncid=ncid,
     $     xid=xid,norbit=norbit,nline=nline,
     $     nelement=nelement,lat=myLat,lon=myLon,
     &     time=time,eaa=Eaa,eia=Eia,pra=pra,tb=tb,sza=mySza,
     $     saa=mySaa,surfalt=mySurfalt,satalt=satalt,qc=qc,
     $     rad=rad,radMW=radMW)

      return
      end subroutine putRad

c--
      subroutine getRad(ncid,irec,xid,norbit,nline,nelement,
     $     lat,lon,time,eaa,eia,pra,tb,sza,saa,surfalt,satalt,qc,
     $     rad)

      !---in/out variables
      integer , intent(in) :: ncid
      integer , intent (in), optional :: irec
      character (len=*), intent (inout), optional :: xid
      integer , intent (inout), optional :: norbit,nline,nelement
      real , intent (inout), optional :: lat,lon
      real , intent (inout), optional :: surfalt
      real , intent (inout), optional :: satalt
      real , dimension(:), intent (inout), optional :: sza,saa
      real , dimension(:), intent (inout), optional :: tb
      real , dimension(:), intent (inout), optional :: rad
      real , dimension(:), intent (inout), optional :: eaa,eia
      real , dimension(:), intent (inout), optional :: pra
      integer , dimension(:), intent (inout), optional :: time,qc
      !-- local arguments
      real , dimension (1) :: myLat,myLon,mySurfAlt
      real :: mySza,mySaa
      
      if (newForm(getNcidIndex(ncid))) then
         if (.not. present(lat)) then 
            myLat = defaultFloat
         else
            myLat = 0.
         endif
         if (.not. present(lon)) then
            myLon = defaultFloat
         else
            myLon = 0.
         endif
         if (.not. present(surfalt)) then 
            mySurfalt = defaultFloat
         else 
            mySurfalt = 0.
         endif
         call getRad1Chan(ncid=ncid,irec=irec,
     $        xid=xid,norbit=norbit,nline=nline,
     $        nelement=nelement,lat=myLat,lon=myLon,
     &        time=time,eaa=eaa,eia=eia,pra=pra,tb=tb,sza=sza,
     $        saa=saa,surfalt=mySurfAlt,satalt=satalt,qc=qc,
     $        rad=rad)
         if (present(lat)) lat = myLat(1)
         if (present(lon)) lon = myLon(1)
         if (present(surfalt)) surfalt = mySurfAlt(1)
      else
         call getNcOldRad(ncid=ncid,ifor=irec,
     $        xid=xid,norbit=norbit,nline=nline,
     $        nelement=nelement,lat=lat,lon=lon,
     &        time=time,eaa=eaa,eia=eia,tb=tb,
     $        sza=mysza,saa=mysaa,surfalt=surfalt)
         sza = mysza
         saa = mysaa
      endif

      return
      end subroutine getRad

c--
      subroutine getNcOldRad(ncid,ifor,xid,norbit,nline,nelement,
     &     lat,lon,time,eaa,eia,tb,sza,saa,surfalt)

      !---in/out variables
      integer , intent(in) :: ncid
      integer , intent(in), optional :: ifor
      character (len=*), intent(inout), optional :: xid
      real , intent(inout), optional :: lat,lon,surfalt
      real , intent(inout), optional :: sza,saa
      real , dimension (:), intent(inout), optional :: eaa,eia,tb
      integer , intent (inout), optional :: norbit,nline,nelement
      integer , dimension (:), intent (inout), optional :: time

      if (present(time)) 
     &     call readNcdfData(ncid,time,varname='date',
     &     record_no=ifor,status='single_record')

      if (present(xid)) 
     &     call readNcdfData(ncid,xid,varname='id',
     &     record_no=ifor,status='single_record')

      if (present(norbit)) 
     &     call readNcdfData(ncid,norbit,varname='norbit',
     &     record_no=ifor,status='single_record')

      if (present(nline)) 
     &     call readNcdfData(ncid,nline,varname='nline',
     &     record_no=ifor,status='single_record')

      if (present(nelement)) 
     &     call readNcdfData(ncid,nelement,varname='nelement',
     &     record_no=ifor,status='single_record')

      if (present(lat)) 
     &     call readNcdfData(ncid,lat,varname='lat',
     &     record_no=ifor,status='single_record')
      
      if (present(lon)) 
     &     call readNcdfData(ncid,lon,varname='lon',
     &     record_no=ifor,status='single_record')
      
      if (present(EIA)) 
     &     call readNcdfData(ncid,EIA,varname='eia',
     &     record_no=ifor,status='single_record')
      
      if (present(EAA)) 
     &     call readNcdfData(ncid,EAA,varname='eaa',
     &     record_no=ifor,status='single_record')
      
      if (present(saa)) 
     &     call readNcdfData(ncid,saa,varname='saa',
     &     record_no=ifor,status='single_record')
      
      if (present(sza)) 
     &     call readNcdfData(ncid,sza,varname='sza',
     &     record_no=ifor,status='single_record')
      
      if (present(surfalt)) 
     &     call readNcdfData(ncid,surfalt,varname='surfalt',
     &     record_no=ifor,status='single_record')

      if (present(tb)) 
     &     call readNcdfData(ncid,tb,varname='tb',
     &     record_no=ifor,status='single_record')

      return
      end subroutine getNcOldRad
c--
      subroutine replaceRad1Chan(ncid,irec,xid,norbit,nline,nelement,
     $     lat,lon,time,eaa,eia,pra,tb,sza,saa,surfalt,satalt,ichan,
     $     qc,rad)

      !---in/out variables
      integer, intent(in) :: ncid
      integer, intent(in) :: irec
      character (len=*), intent (in), optional :: xid
      integer, intent (in), optional :: norbit,nline,nelement,ichan
      real, dimension (:), intent (in), optional :: lat,lon
      real, dimension (:), intent (in), optional :: surfalt
      real, intent (in), optional :: satalt
      real, dimension (:), intent (in), optional :: sza,saa
      real, dimension (:), intent (in), optional :: eaa,eia,tb
      real, dimension (:), intent (in), optional :: rad
      real, dimension (:), intent (in), optional :: pra
      integer, dimension (:), intent (in), optional :: time, qc
      integer :: index

      if (present(tb))
     $     call replaceNcdfData(fid=ncid,var=tb,varname='tb',
     $     rec=irec,status='append')
         
      if (present(rad)) 
     &     call replaceNcdfData(fid=ncid,var=rad,varname='rad',
     $     rec=irec,status='append')
         
      if (present(time))
     $     call replaceNcdfData(fid=ncid,var=time,varname='date',
     &     rec=irec,status='append')

      if (present(xid)) 
     &     call replaceNcdfData(fid=ncid,var=xid,varname='id',
     $     rec=irec,status='append')

      if (present(norbit)) 
     &     call replaceNcdfData(fid=ncid,var=norbit,varname='norbit',
     &     rec=irec,status='append')
         
      if (present(nline)) 
     &     call replaceNcdfData(fid=ncid,var=nline,varname='nline',
     &     rec=irec,status='append')

      if (present(nelement)) 
     &     call replaceNcdfData(fid=ncid,var=nelement,
     &     varname='nelement',rec=irec,status='append')

      if (present(ichan)) then
         index= searchID(ncid) 
         if (FovTypeLoc(index) == 0) then
            print *,'err[rad_io_module::putRad1Chan]: ',
     &           ' unmatched FOV type'
            call errorHalt(1)
         end if

         call replaceNcdfData(fid=ncid,var=ichan,varname='ichan',
     &        rec=irec,status='append')
      endif
         
      if (present(lat)) then
         if (any(lat(:) /= defaultFloat))
     $        call replaceNcdfData(fid=ncid,var=lat,varname='lat',
     &        rec=irec,status='append')
      endif
         
      if (present(lon)) then
         if (any(lon(:) /= defaultFloat))
     $        call replaceNcdfData(fid=ncid,var=lon,varname='lon',
     &        rec=irec,status='append')
      endif

      if (present(eia))
     $     call replaceNcdfData(fid=ncid,var=eia,varname='eia',
     &     rec=irec,status='append')

      if (present(eaa)) 
     &     call replaceNcdfData(fid=ncid,var=eaa,varname='eaa',
     &     rec=irec,status='append')
      
      if (present(pra))
     $     call replaceNcdfData(fid=ncid,var=pra,varname='pra',
     &     rec=irec,status='append')

      if (present(saa)) then
         if (any(saa(:) /= defaultFloat))
     &        call replaceNcdfData(fid=ncid,var=saa,varname='saa',
     &        rec=irec,status='append')
      endif

      if (present(sza)) then
         if (any(sza(:) /= defaultFloat))
     &        call replaceNcdfData(fid=ncid,var=sza,varname='sza',
     &        rec=irec,status='append')
      endif

      if (present(surfalt)) then
         if (any(surfalt(:) /= defaultFloat))
     &        call replaceNcdfData(fid=ncid,var=surfalt,
     &        varname='surfalt',rec=irec,status='append')
      endif

      if (present(satalt))
     $     call replaceNcdfData(fid=ncid,var=satalt,varname='satalt',
     &     rec=irec,status='append')

      if (present(qc))
     $     call replaceNcdfData(fid=ncid,var=qc,varname='qc',
     &     rec=irec,status='append')

      return
      end subroutine replaceRad1Chan
c--
      subroutine replaceRad(ncid,irec,xid,norbit,nline,nelement,
     $     lat,lon,time,eaa,eia,pra,tb,sza,saa,surfalt,satalt,qc,rad)

      !---in/out variables
      integer, intent(in) :: ncid
      integer, intent(in) :: irec
      character (len=*), intent (in), optional :: xid
      integer, intent (in), optional :: norbit,nline,nelement
      real, dimension(:), intent (in), optional :: tb
      real, intent (in), optional :: lat,lon
      real, intent (in), optional :: surfalt
      real, intent (in), optional :: satalt
      real, intent (in), optional :: sza,saa
      real, dimension(:), intent (in), optional :: rad
      real, dimension(:), intent (in), optional :: eaa,eia
      real, dimension(:), intent (in), optional :: pra
      integer, dimension(:), intent (in), optional :: time,qc
      !-- local arguments
      real, dimension (1) :: myLat,myLon,mySurfalt
      real, dimension (1) :: mySza,mySaa
      
      if (present(lat)) then 
         mylat(1) = lat
      else
         mylat(1) = defaultFloat
      endif
      if (present(lon)) then
         mylon(1) = lon
      else
         mylon(1) = defaultFloat
      endif
      if (present(sza)) then
         mySza(1) = sza
      else
         mySza(1) = defaultFloat
      endif

      if (present(saa)) then
         mySaa(1) = saa
      else
         mySaa(1) = defaultFloat
      endif
      if (present(surfalt)) then 
         mySurfalt(1) = surfalt
      else
         mySurfalt(1) = defaultFloat
      endif

      call replaceRad1Chan(ncid=ncid,irec=irec,
     $     xid=xid,norbit=norbit,nline=nline,
     $     nelement=nelement,lat=myLat,lon=myLon,
     &     time=time,eaa=Eaa,eia=Eia,pra=pra,tb=tb,sza=mySza,
     $     saa=mySaa,surfalt=mySurfalt,satalt=satalt,qc=qc,rad=rad)

      return
      end subroutine replaceRad
c--
      subroutine appendRad1Chan(ncid,xid,norbit,nline,nelement,
     $     lat,lon,time,eaa,eia,pra,tb,sza,saa,surfalt,satalt,ichan,qc,
     $     rad)

      !---in/out variables
      integer, intent(in) :: ncid
      character (len=*), intent (in), optional :: xid
      integer, intent (in), optional :: norbit,nline,nelement,ichan
      real, dimension (:), intent (in), optional :: lat,lon
      real, dimension (:), intent (in), optional :: surfalt
      real, intent (in), optional :: satalt
      real, dimension (:), intent (in), optional :: sza,saa
      real, dimension (:), intent (in), optional :: eaa,eia,tb
      real, dimension (:), intent (in), optional :: pra
      real, dimension (:), intent (in), optional :: rad
      integer, dimension (:), intent (in), optional :: time, qc
      integer :: rec

      rec = readNcdfDim(fid=ncid)+1 !get the new unlimited dimension
         
      call replaceRad1Chan(ncid=ncid,irec=rec,
     $     xid=xid,norbit=norbit,nline=nline,
     $     nelement=nelement,lat=lat,lon=lon,
     &     time=time,eaa=eaa,eia=eia,pra=pra,tb=tb,sza=sza,
     $     saa=saa,surfalt=surfalt,satalt=satalt,qc=qc,rad=rad)

      return
      end subroutine appendRad1Chan
c--
      subroutine appendRad(ncid,xid,norbit,nline,nelement,
     $     lat,lon,time,eaa,eia,pra,tb,sza,saa,surfalt,satalt,qc,
     $     rad)

      !---in/out variables
      integer, intent(in) :: ncid
      character (len=*), intent (in), optional :: xid
      integer, intent (in), optional :: norbit,nline,nelement
      real, dimension(:), intent (in), optional :: tb
      real, intent (in), optional :: lat,lon
      real, intent (in), optional :: surfalt
      real, intent (in), optional :: satalt
      real, intent (in), optional :: sza,saa
      real, dimension(:), intent (in), optional :: rad
      real, dimension(:), intent (in), optional :: eaa,eia
      real, dimension(:), intent (in), optional :: pra
      integer, dimension(:), intent (in), optional :: time,qc
      !-- local arguments
      real, dimension (1) :: myLat,myLon,mySurfalt
      real, dimension (1) :: mySza,mySaa
      
      if (present(lat)) then 
         mylat(1) = lat
      else
         mylat(1) = defaultFloat
      endif
      if (present(lon)) then
         mylon(1) = lon
      else
         mylon(1) = defaultFloat
      endif
      if (present(sza)) then
         mySza(1) = sza
      else
         mySza(1) = defaultFloat
      endif

      if (present(saa)) then
         mySaa(1) = saa
      else
         mySaa(1) = defaultFloat
      endif
      if (present(surfalt)) then 
         mySurfalt(1) = surfalt
      else
         mySurfalt(1) = defaultFloat
      endif

      call appendRad1Chan(ncid=ncid,
     $     xid=xid,norbit=norbit,nline=nline,
     $     nelement=nelement,lat=myLat,lon=myLon,
     &     time=time,eaa=Eaa,eia=Eia,pra=pra,tb=tb,sza=mySza,
     $     saa=mySaa,surfalt=mySurfalt,satalt=satalt,qc=qc,
     $     rad=rad)

      return
      end subroutine appendRad
c----
      integer function searchID(ncid) result(index)
      integer, intent(in) :: ncid
      integer :: i
      index = -1
      do i=1,size(fid)
         if (fid(i) == ncid) then 
            index = i
            return
         endif
      enddo
      return
      end function searchID
c---- 
      subroutine initFOVtype(ncid,FOVtype)
      integer, intent (in) :: ncid, FOVtype
      integer :: index

      index = searchID(ncid)
      if (index < 0) then       !newly opened: create the map for new ncid
         fidIndex = fidIndex + 1
      else                      !previously opened: re-use the slot
         fidIndex = index
      endif
      fid(fidIndex) = ncid
      FovTypeLoc(fidIndex) = FOVtype
      end subroutine initFOVtype
c----
      end module rad_io_module
