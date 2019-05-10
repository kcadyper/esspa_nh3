c $"$Id$"
      module bkg_io_module
c**********************************************************************
c* Module name: bkg_io_module
c* Purpose: Handles covariance and background I/O
c* Functions:
c*           openCov
c*           getCov
c*           putCov
c*           closeCov
c*
c* Var_Name      Description
c* --------      -----------
c* ncid          file id
c* file          file name
c* creationDate  date of covariance is created
c* lcascd        cascade level
c* fovsz         fov cell size
c* nbkg          number of backgrounds
c* ifstruc       background file structure
c* iftran        geophys(0) or eigen(1) space flag
c* NParG         number of parameters in geophysical domain
c* nLevel        number of pressure levels
c* nchmw         number of MW channels
c* nemir         number of hinge-points for IR surface em & rlf
c* pressure      pressure grid
c* freq          channel frequencies (GHZ)
c* frqEmir       hinge-point frequencies for IR
c* polarity      channel polarities
c* status        directs routine to open a new or old netCDF file
c* irec          record number for directly accesing data record
c* dmean         background profile
c* un_mw         eigenvectors of parameter
c* cov           covariance matrix (dcov or ucov)
c* eiabkg        earth incidence angle per channel
c* TypeFlag      background type
c* bkgID         background identifier
c* ntran         number of parameters after transformation
c* date          date/time tag of background
c* msg           closing message 
c* IC        StateIndex_t  state vector indices in transformed space
c* NC        StateIndex_t  state vector lengths in transformed space
c* IG        StateIndex_t  state vector indices in geophysical space
c* NG        StateIndex_t  state vector lengths in geophysical space
c* qc            array of quality control words
c* 
c* Externals:
c* ncdf_module   base module for netCDF I/O
c* 
c* Developed by AER 2003
c* Copyright: AER, Inc., 2003
c**********************************************************************
      use ncdf_module
      use StateIndexModule

      implicit none
      integer :: imol, fidIndexG=0, fidIndexC=0, index
      integer , dimension(:,:), allocatable :: mapFromG, mapFromC
      integer , dimension(MAXFILE) :: fidG=-1, fidC=-1, 
     $     mapSizeG, mapSizeC, fileSizeG, fileSizeC
      logical , dimension(MAXFILE) :: stateIndexEnabled
      logical , dimension(MAXFILE) :: fileEOF
      integer, dimension(:), pointer :: mapBeamLoc
      integer :: nBeamLoc, nAngLoc
      type(StateIndex_t) :: IGF, NGF, ICF, NCF
      
      private
      integer, save :: record_no=0
      public openCov, closeCov, getCov, putCov, queryCov
c-----
      interface openCov
      module procedure openNcFile
      end interface

      interface queryCov
      module procedure queryNcCov
      end interface

      interface getCov
      module procedure getNcCov
      end interface

      interface putCov
      module procedure putNcCov
      end interface

      interface closeCov
      module procedure closeNcFile
      end interface
c---- 
      CONTAINS
c--   
      subroutine openNcFile(ncid, file, creationDate, lcascd, 
     $     fovsz, nbkg, ifstruc, iftran, invNotTrnsp, NParG, ntran,
     $     nLevel, nchmw, nemir, nTranMW, nTranIR, nAng, nBeam, mapBeam,
     $     pressure, freq, frqEmir, polarity,
     $     icell, eiabkg, bkgIdLen,IC,NC,IG,NG,molID,MolTran,
     $     nqc,vCoordTyp,sigPtop,sigCoord,hybCoordA,hybCoordB,status)

      use MapBeamDefaults
      use constants, only: MISSING_INT

c--   Dummy variables
      integer, intent(inout) :: ncid
      character (len=*), intent(in) :: file
      character (len=*), intent (inout), optional :: creationDate
      integer, intent(inout), optional :: lcascd, nbkg, nqc,
     $     ifstruc, iftran, invNotTrnsp, 
     $     NParG, nLevel, nchmw, nemir, icell
      integer, intent(inout), optional :: nBeam, nAng
      real, intent (inout), optional :: fovsz
      integer, dimension(:), intent(inout), optional :: mapBeam
      real, dimension(:), intent (inout), optional :: pressure
      real, dimension(:), intent (inout), optional :: freq,frqEmir
      integer, dimension(:), intent (inout), optional :: polarity,molID
      integer, dimension(:), intent (inout), optional :: MolTran
      real, dimension (:), intent (inout), optional :: eiabkg
      integer, intent (inout), optional :: bkgIdLen, ntran
      integer, intent (inout), optional :: nTranMW, nTranIR
      character (len=*),   intent (inout), optional :: vCoordTyp
      real,                intent (inout), optional :: sigPtop
      real, dimension (:), intent (inout), optional :: sigCoord
      real, dimension (:), intent (inout), optional :: hybCoordA
      real, dimension (:), intent (inout), optional :: hybCoordB
      character (len=*), intent(in), optional :: status
      type(StateIndex_t), intent (inout), optional :: IG, NG
      type(StateIndex_t), intent (inout), optional :: IC, NC

c--   Local variables
      character (len=MAX_NAME_LENGTH) :: lstatus
      integer, dimension(8) :: sysvalues
      character (len= 8) :: sysdate
      character (len=10) :: systime
      character (len= 5) :: syszone
      integer, dimension(2) :: myDim
      character(40) :: datestring
      integer , dimension (2,maxMol) :: mol
      integer , dimension(maxMol) :: MolIDF,MolLenF,TMolIDF,TMolLenF
      integer , dimension(maxMol) :: MolTranF
      integer :: nchmwloc,nemirloc, NMol
      logical :: molTranFound
      integer :: invNotTrnspLoc
      logical :: needIC
      logical :: vCfound
      integer :: ntmp
      type(StateIndex_t) :: ICtmp, NCtmp

      if (.not. present(status)) then
         lstatus = 'old'
      else
         lstatus = trim(status)
      endif

      needIC=.FALSE.

      if (trim(lstatus) == 'new' .or. trim(lstatus) == 'NEW') then

         if ((present(freq) .or. present(polarity)) .and.
     $        .not. (present(nchmw))) then
            print *,'err[bkg_io_module::openNcFile]: ',
     $           ' openCov: missing nchmw'
            call errorHalt(1)
         end if

         if ((present(frqEmir)).and.
     $        .not. (present(nemir))) then
            print *,'err[bkg_io_module::openNcFile}: ',
     $           ' openCov: missing nemir'
            call errorHalt(1)
         end if

         call openNcdfFile(ncid,file,status='new', 
     $        unlimited_dim_name='nbkg')

         call initFileID(ncid=ncid,fidIndex=fidIndexG,
     $        stateIndexEnabled=stateIndexEnabled,IG=IG,NG=NG,fid=fidG)

         call date_and_time(sysdate, systime, syszone, sysvalues)

         if (present(creationDate))
     $        call writeNcdfAttr(ncid,attrName='Creationtime',
     $        attr=trim(sysdate//systime))

         if (present(lcascd))
     $        call writeNcdfAttr(ncid,attrName='CascadeLevel',
     &        attr=(/lcascd/))

         if (present(fovsz))
     $        call writeNcdfAttr(ncid,attrName='FovCellSize',
     &        attr=(/fovsz/))

         if (present(ifstruc))
     $        call writeNcdfAttr(ncid,
     $        attrName='FileStructure',
     &        attr=(/ifstruc/))

         if (present(nLevel))
     $        call writeNcdfAttr(ncid,attrName='nlev',
     &        attr=(/nLevel/))

         if (present(nchmw)) then
            nchmwLoc=nchmw
            call writeNcdfAttr(ncid,attrName='nchmw',attr=nchmwLoc)
         
            call mapBeamDefNew(nchmwLoc,nBeam,mapBeam,nBeamLoc,
     $           mapBeamLoc)

            if (present(nAng)) then
               nAngLoc = nAng
            else
               nAngLoc=nBeamLoc
            endif
            if (nAngLoc > 0)
     $        call writeNcdfDim(ncid,dimName='nAng',dimLen=nAngLoc)
            call writeNcdfAttr(ncid,attrName='nBeam',attr=nBeamLoc)
            call writeNcdfAttr(ncid,attrName='mapBeam',attr=mapBeamLoc)
         endif

         if (present(nemir))
     $        call writeNcdfAttr(ncid,attr=nemir,attrName='nemir')

         if (present(icell))
     $        call writeNcdfAttr(ncid,attrName='icell',
     &        attr=(/icell/))

         if (present(invNotTrnsp)) then
            call writeNcdfAttr(ncid,attr=invNotTrnsp,
     &          attrName='invNotTrnspFlag')
            if (invNotTrnsp == 0) then
               fileEOF(fidIndexG)=.TRUE.
            else
               fileEOF(fidIndexG)=.FALSE.
               needIC=.TRUE.
            endif
         else
            fileEOF(fidIndexG)=.TRUE. ! assume trans matrix is EOF by default
         endif

         if (present(iftran)) then
            call writeNcdfAttr(ncid,attrName='TransFlag',
     $           attr=(/iftran/))
            if (iftran == 1) needIC=.TRUE.
         endif

         if (needIC) then
            if (present(IC) .and. present(NC)) then
               if (.not. stateIndexEnabled(fidIndexG)) then
                  print *,'err[bkg_io_module::openNcFile]: ',
     $                 ' openCov: stateIndex not enabled ',
     $                 'in needIC write mode'
                  call errorHalt(1)
               end if
               
               if (.not. present(MolID)) then
                  print *,'err[bkg_io_module::openNcFile]: ',
     $                 ' openCov: missing MolID in output.'
                  call errorHalt(1)
               end if

               NCF = expandMol(NC,MolID)
               ICF = genIndices(NCF)

               call writeNcdfAttr(ncid,attrName='TTemperature',
     &              attr=(/ICF%Temp,NCF%Temp/))
            
               call writeNcdfAttr(ncid,attrName='TTskin',
     &              attr=(/ICF%Tskin,NCF%Tskin/))

               call writeNcdfAttr(ncid,attrName='TSurfacePressure',
     &              attr=(/ICF%Psfc,NCF%Psfc/))
               
               call writeNcdfAttr(ncid,attrName='TMol',
     $              attr=ICF%mol(1:maxMol))
         
               call writeNcdfAttr(ncid,attrName='TMolLen',
     $              attr=NCF%mol(1:maxMol))
         
               call writeNcdfAttr(ncid,attrName='TLiqCloud',
     &              attr=(/ICF%cldLiq,NCF%cldLiq/))
               
               call writeNcdfAttr(ncid,attrName='TIceCloud',
     &              attr=(/ICF%cldIce,NCF%cldIce/))

               call writeNcdfAttr(ncid,attrName='TSurfaceWinds',
     &              attr=(/ICF%Wind,NCF%Wind/))
            else
               print*, 'warning[bkg_io_module::openNcFile]: ',
     $              'incomplete I/O args -- no transf index output.'
            endif
         endif

         if (present(bkgIdLen))
     $        call writeNcdfDim(ncid,dimLen=bkgIdLen,
     &        dimName='bkgIdLen')

         if (stateIndexEnabled(fidIndexG)) then
            if (.not. present(MolID)) then
               print *,'err[bkg_io_module::openNcFile]: ',
     $              ' openCov: missing MolID in output.'
               call errorHalt(1)
            end if

            NGF = expandMol(NG,MolID)
            IGF = genIndices(NGF)

            call writeNcdfAttr(ncid,attrName='Temperature',
     &           attr=(/IGF%Temp,NGF%Temp/))

            call writeNcdfAttr(ncid,attrName='Tskin',
     &           attr=(/IGF%Tskin,NGF%Tskin/))

            call writeNcdfAttr(ncid,attrName='SurfacePressure',
     &           attr=(/IGF%Psfc,NGF%Psfc/))

            call writeNcdfAttr(ncid,attrName='Mol',
     $           attr=IGF%mol(1:maxMol))
         
            call writeNcdfAttr(ncid,attrName='MolLen',
     $           attr=NGF%mol(1:maxMol))
         
            call writeNcdfAttr(ncid,attrName='LiqCloud',
     &           attr=(/IGF%cldLiq,NGF%cldLiq/))

            call writeNcdfAttr(ncid,attrName='IceCloud',
     &           attr=(/IGF%cldIce,NGF%cldIce/))

            call writeNcdfAttr(ncid,attrName='SurfaceWinds',
     &           attr=(/IGF%Wind,NGF%wind/))

            if (present(MolTran)) then
               MolTranF=0
               NMol = getNMol(NGF)
               MolTranF(MolID(1:NMol)) = MolTran(1:NMol)
               call writeNcdfAttr(ncid,attrName='MolTran',
     $              attr=MolTranF(1:maxMol))
            endif
         endif

         if (present(vCoordTyp))
     &        call writeNcdfAttr(ncid,attrName='VertCoordType',
     &        attr=TRIM(vCoordTyp))

         if (present(pressure))
     $        call writeNcdfAttr(ncid,
     $        attrName='StandardPressureGrid',
     &        attr=pressure)

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

         if (present(freq))
     $        call writeNcdfAttr(ncid,attrName='mwfrequencies',
     &        attr=freq)

         if (present(frqEmir))
     $        call writeNcdfAttr(ncid,attrName='wavenumbers',
     &        attr=frqEmir)

         if (present(polarity))
     $        call writeNcdfAttr(ncid,attrName='mwpolarities',
     &        attr=polarity)

         if (present(eiabkg)) 
     $        call writeNcdfAttr(ncid,attrName='eiabkg',attr=eiabkg)

      else
         call openNcdfFile(ncid,file,status='old')

         call initFileID(ncid=ncid,fidIndex=fidIndexG,
     $        stateIndexEnabled=stateIndexEnabled,IG=IG,NG=NG,fid=fidG)

         if (present(creationDate))
     $        call readNcdfAttr(ncid,attrName='Creationtime',
     &        attr=creationDate)

         if (ncid <= 0) then 
            print *, 'warning[bkg_io_module::openCov]:  Invalid file ID'
            return
         endif

         if (present(lcascd))
     $        call readNcdfAttr(ncid,attrName='CascadeLevel',
     &        attr=lcascd)

         if (present(fovsz))
     $        call readNcdfAttr(ncid,attrName='FovCellSize',
     &        attr=fovsz)

         if (present(nbkg)) nbkg=readNcdfDim(ncid,dimName='nbkg')

         if (present(nqc)) nqc=readNcdfDim(ncid,dimName='nQC')

         if (present(bkgIdLen)) 
     $        bkgIdLen=readNcdfDim(ncid,dimName='bkgIdLen')

         if (present(ifstruc))
     $        call readNcdfAttr(ncid,
     $        attrName='FileStructure',
     &        attr=ifstruc)

         call readNcdfAttr(ncid,attrName='invNotTrnspFlag',
     $        attr=invNotTrnspLoc,silent=.true.) ! assume inv is transpose (EOF)
         if (present(invNotTrnsp)) invNotTrnsp=invNotTrnspLoc
         if (invNotTrnspLoc == 0) then
            fileEOF(fidIndexG)=.TRUE.
         else
            fileEOF(fidIndexG)=.FALSE.
            needIC=.TRUE.
         endif

         if (present(iftran)) then
            call readNcdfAttr(ncid,attrName='TransFlag',
     &           attr=iftran)
            if (iftran == 1) needIC=.TRUE.
         endif

         ntmp=readNcdfDim(ncid,dimName='NParG',silent=.TRUE.)
         if (ntmp <= 0) needIC=.FALSE. ! No atmosphere data

         if (needIC) then
               
            call readNcdfAttr(ncid,attrName='TTemperature',
     $           attr=mydim)
            ICF%Temp = mydim(1)
            NCF%Temp = mydim(2)
            
            call readNcdfAttr(ncid,attrName='TTskin',attr=mydim)
            ICF%Tskin = mydim(1)
            NCF%Tskin = mydim(2)
            
            call readNcdfAttr(ncid,attrName='TSurfacePressure',
     $           attr=mydim)
            ICF%Psfc=mydim(1)
            NCF%Psfc=mydim(2)
            
            call readNcdfAttr(ncid,attrName='TMol',
     $           attr=TMolIDF)
            ICF%mol = TMolIDF
            
            call readNcdfAttr(ncid,attrName='TMolLen',
     $           attr=TMolLenF)
            NCF%mol = TMolLenF

            call readNcdfAttr(ncid,attrName='TLiqCloud',
     &           attr=mydim)
            ICF%cldLiq=mydim(1)
            NCF%cldLiq=mydim(2)
            
            call readNcdfAttr(ncid,attrName='TIceCloud',
     &           attr=mydim)
            ICF%cldIce=mydim(1)
            NCF%cldIce=mydim(2)

            call readNcdfAttr(ncid,attrName='TSurfaceWinds',
     $           attr=mydim)
            ICF%Wind=mydim(1)
            NCF%Wind=mydim(2)

         endif

         if (present(NParG))
     $        NParG=readNcdfDim(ncid,dimName='NParG')

         if (stateIndexEnabled(fidIndexG)) then
            call readNcdfAttr(ncid,attrName='Temperature',attr=mydim)
            IGF%Temp = mydim(1)
            NGF%Temp = mydim(2)
         
            call readNcdfAttr(ncid,attrName='Tskin',attr=mydim)
            IGF%Tskin = mydim(1)
            NGF%Tskin = mydim(2)
         
            call readNcdfAttr(ncid,attrName='SurfacePressure',
     &           attr=mydim)
            IGF%Psfc=mydim(1)
            NGF%Psfc=mydim(2)
         
            call readNcdfAttr(ncid,attrName='Mol',attr=MolIDF)
            IGF%mol = MolIDF

            call readNcdfAttr(ncid,attrName='MolLen',attr=MolLenF)
            NGF%mol = MolLenF

            call readNcdfAttr(ncid,attrName='LiqCloud',attr=mydim)
            IGF%cldLiq=mydim(1)
            NGF%cldLiq=mydim(2)
         
            call readNcdfAttr(ncid,attrName='IceCloud',attr=mydim)
            IGF%cldIce=mydim(1)
            NGF%cldIce=mydim(2)

            call readNcdfAttr(ncid,attrName='SurfaceWinds',attr=mydim)
            IGF%Wind=mydim(1)
            NGF%Wind=mydim(2)
         endif
         
         if (present(ntran))
     $        ntran=readNcdfDim(ncid,dimName='ntran')

         if (present(nTranMW))
     $        nTranMW=readNcdfDim(ncid,dimName='ntranEmMw')

         if (present(nTranIR))
     $        nTranIR=readNcdfDim(ncid,dimName='ntranEmIR')

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
     $        attrName='StandardPressureGrid',
     &        attr=pressure)

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
            call readNcdfAttr(ncid,attrName='mwfrequencies',
     &           attr=freq)
            call readNcdfAttr(ncid,attr=nchmwloc,attrName='nchmw')
            if (size(freq) /= nchmwloc) then
               print *,'err[bkg_io_module::openNcFile]: ',
     $              ' openCov: #freq and nchmw mismatched'
               call errorHalt(1)
            end if
         endif

         if (present(polarity)) then
            call readNcdfAttr(ncid,attrName='mwpolarities',
     &           attr=polarity)
            call readNcdfAttr(ncid,attr=nchmwloc,attrName='nchmw')
            if (size(freq) /= nchmwloc) then
               print *,'err[bkg_io_module::openNcFile]: ',
     $              ' openCov: #mwpol and nchmw mismatched'
               call errorHalt(1)
            end if
         endif

         call readNcdfAttr(ncid,attrName='nchmw',attr=nchmwLoc,
     $        silent=.true.)
         if (present(nchmw)) nchmw=nchmwLoc

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
            if (present(mapBeam)) mapBeam = mapBeamLoc
            deallocate(mapBeamLoc)
         else
            if (present(mapBeam)) mapBeam = MISSING_INT
         endif
         if (present(nBeam)) nBeam = nBeamLoc

         if (present(nAng)) then
            if (nBeamLoc == 0) then ! check the dependency of eia on beams
               nAng = 0
            else
               nAng = readNcdfDim(ncid,'nAng',silent=.TRUE.)
            endif
         endif

         if (present(frqEmir)) then
            call readNcdfAttr(ncid,attrName='wavenumbers',
     &           attr=frqEmir,silent=.true.)
            if (all(frqEmir == 0.))    ! for back compatibility
     &         call readNcdfAttr(ncid,attrName='irfrequencies',
     &              attr=frqEmir)
            call readNcdfAttr(ncid,attr=nemirloc,attrName='nemir')
            if (size(frqEmir) /= nemirloc) then
               print *,'err[bkg_io_module::openNcFile]: ',
     $              ' openCov: #freq and nemir mismatched'
               call errorHalt(1)
            end if
         endif

         if (present(nemir))
     $        call readNcdfAttr(ncid,attr=nemir,attrName='nemir')

         if (present(eiabkg)) 
     $        call readNcdfAttr(ncid,attrName='eiabkg',attr=eiabkg)

         if (needIC) then
            if (.not. allocated(mapFromC))
     $           allocate(mapFromC(MAXPAR,MAXFILE))
            index = searchID(fidC,ncid)
            if (index < 0) then 
               fidIndexC = fidIndexC + 1
            else
               fidIndexC = index
            endif
            fidC(fidIndexC) = ncid
            call genMapFrom(IG=ICtmp,NG=NCtmp,IGF=ICF,NGF=NCF,
     $           MolID=MolID,map=mapFromC(:,fidIndexC),
     $           mapSize=mapSizeC(fidIndexC),
     $           xFileSize=fileSizeC(fidIndexC))
         endif

         if (present(IC) .and. present(NC)) then
            if (.not. needIC) then
               print *,'err[bkg_io_module::openNcFile]: ',
     $              ' openCov: IC and NC not available in file '
               call errorHalt(1)
            end if
            IC = ICtmp
            NC = NCtmp
         endif

         if (stateIndexEnabled(fidIndexG)) then
            if (.not. allocated(mapFromG)) 
     $           allocate(mapFromG(MAXPAR,MAXFILE))
            call genMapFrom(IG=IG,NG=NG,IGF=IGF,NGF=NGF,
     $           MolID=MolID,map=mapFromG(:,fidIndexG),
     $           mapSize=mapSizeG(fidIndexG),
     $           xFileSize=fileSizeG(fidIndexG))
            if (.not. needIC) then
               if (.not. allocated(mapFromC))
     $              allocate(mapFromC(MAXPAR,MAXFILE))
               index = searchID(fidC,ncid)
               if (index < 0) then !newly opened
                  fidIndexC = fidIndexC + 1
               else             !previously opened
                  fidIndexC = index
               endif
               fidC(fidIndexC) = ncid
               mapFromC(:,fidIndexC) = mapFromG(:,fidIndexG)
               mapSizeC(fidIndexC) = mapSizeG(fidIndexG)
               fileSizeC(fidIndexC) = fileSizeG(fidIndexG)
            endif
         endif                  !stateIndexEnabled

         if (present(MolTran)) then
            if (.not. stateIndexEnabled(fidIndexG)) then
               print *,'err[bkg_io_module::openNcFile]: ',
     $              ' openCov: stateIndex not enabled ',
     $              'to read MolTran'
               call errorHalt(1)
            end if
            if (.not. present(MolID)) then
               print *,'err[bkg_io_module::openNcFile]: ',
     $              ' openCov: MolID is needed ',
     $              'to get MolTran'
               call errorHalt(1)
            end if
            call readNcdfAttr(ncid,attrName='MolTran',
     $           attr=MolTranF,silent=.true.,found=molTranFound)
            if (NGF%mol(1)>0 .and. (.not. molTranFound)) then
               MolTranF(1)=1 ! for back compatibility
               print *,'warning[bkg_io_module::openNcFile]: '//
     $            'molTran missing from file; assuming H2O is log'
            endif
            NMol = getNMol(NGF)
            MolTran(1:NMol) = MolTranF(MolID(1:NMol))
         endif
      endif                     !new file

      end subroutine openNcFile
c--
      subroutine queryNcCov(file, creationDate, lcascd, 
     $     fovsz, nbkg, ifstruc, iftran, invNotTrnsp, 
     $     NParG, ntran, nLevel,
     $     nchmw, nemir, nTranMW, nTranIR, nAng, nBeam, mapBeam, 
     $     pressure, freq, frqEmir, polarity, 
     $     icell, eiabkg, bkgIdLen, IC, NC, IG, NG, NMol, molID, 
     $     MolTran,nqc,vCoordTyp,sigPtop,sigCoord,hybCoordA,hybCoordB)

c--   Dummy variables
      character (len=*), intent(in) :: file
      character (len=*), intent (inout), optional :: creationDate
      integer, intent(inout), optional :: lcascd, nbkg, nqc,
     $     ifstruc, iftran, invNotTrnsp, 
     $     NParG, nLevel, nchmw, nemir, icell
      integer, intent(inout), optional :: nBeam, nAng
      real, intent (inout), optional :: fovsz
      integer, dimension(:), intent(inout), optional :: mapBeam
      real, dimension(:), intent (inout), optional :: pressure
      real, dimension(:), intent (inout), optional :: freq,frqEmir
      integer, dimension(:), intent (inout), optional :: polarity,molID
      integer, dimension(:), intent (inout), optional :: MolTran
      real, dimension (:), intent (inout), optional :: eiabkg
      integer, intent (inout), optional :: bkgIdLen, ntran
      integer, intent (inout), optional :: nTranMW, nTranIR
      type(StateIndex_t), intent (inout), optional :: IG, NG
      type(StateIndex_t), intent (inout), optional :: IC, NC
      integer , intent(inout), optional :: NMol
      character (len=*),   intent (inout), optional :: vCoordTyp
      real,                intent (inout), optional :: sigPtop
      real, dimension (:), intent (inout), optional :: sigCoord
      real, dimension (:), intent (inout), optional :: hybCoordA
      real, dimension (:), intent (inout), optional :: hybCoordB

c--   Local variables
      integer :: ncid

      if (present(molID)) molID=nullMolID
      if (present(NMol)) NMol=0
      call openNcFile(ncid=ncid, file=file, creationDate=creationDate, 
     $     lcascd=lcascd, fovsz=fovsz, nbkg=nbkg, ifstruc=ifstruc, 
     $     iftran=iftran, invNotTrnsp=invNotTrnsp, 
     $     NParG=NParG, ntran=ntran, nLevel=nLevel,
     $     nchmw=nchmw, nemir=nemir, nTranMW=nTranMW, nTranIR=nTranIR,
     $     nAng=nAng, nBeam=nBeam, mapBeam=mapBeam,
     $     pressure=pressure, freq=freq, 
     $     frqEmir=frqEmir, polarity=polarity, 
     $     icell=icell, eiabkg=eiabkg, 
     $     bkgIdLen=bkgIdLen, IC=IC, NC=NC, IG=IG, NG=NG, 
     $     molID=molID, MolTran=MolTran,nqc=nqc,vCoordTyp=vCoordTyp,
     $     sigPtop=sigPtop,
     $     sigCoord=sigCoord,hybCoordA=hybCoordA,hybCoordB=hybCoordB)
      call closeNcFile(ncid)
      end subroutine queryNcCov
c--         
      subroutine getNcCov(ncid, irec, dmean, un_atm, vn_atm, cov,
     $     dmEmMw, un_EmMw, vn_EmMw, covEmMw, 
     $     dmEmIR, un_EmIR, vn_EmIR, covEmIR,
     $     bkgID, TypeFlag, date, qc)

c-- Dummy arguments
      integer, intent (in) :: ncid
      integer, intent (in), optional :: irec
      real, dimension (:), intent (inout), optional :: dmean, dmEmMw,
     $     dmEmIR
      real, dimension (:,:), intent (inout), optional :: un_atm,un_EmMw,
     $     un_EmIR    
      real, dimension (:,:), intent (inout), optional :: vn_atm,vn_EmMw,
     $     vn_EmIR    
      real, dimension (:,:), intent (inout), optional :: cov,covEmMw,
     $     covEmIR
      integer, intent (inout), optional :: TypeFlag
      character (len=*), intent (inout), optional :: bkgID
      integer, dimension (:), intent (inout), optional :: date, qc
c-- Local arguments
      integer :: indexG, indexC
      real , dimension (:), allocatable :: dmeanF
      real , dimension (:,:), allocatable :: un_atmF,vn_atmF
      real , dimension (:,:), allocatable :: covF
      real , dimension (:,:), allocatable :: un_EmMwF,un_EmIRF

      if (.not. present(irec)) then
         record_no=record_no + 1
      else
         record_no = irec
      endif

      if (present(TypeFlag))
     $     call readNcdfData(ncid,TypeFlag,varname='TypeFlag',
     $     record_no=record_no,status='single_record')

      if (present(bkgID))
     $     call readNcdfData(ncid,bkgID,varname='backgrdID',
     $     record_no=record_no,status='single_record')

      if (present(date))
     $     call readNcdfData(ncid,date,varname='date',
     $     record_no=record_no,status='single_record')

      if (present(dmean)) then
         indexG = searchID(fidG,ncid)
         if (.not. stateIndexEnabled(indexG)) then
            print *,'err[bkg_io_module::getNcCov]: ',
     $           ' stateIndex not enabled for dmean.'
            call errorHalt(1)
         end if

         if (.not. allocated(dmeanF)) allocate(dmeanF(MAXPAR))
         call readNcdfData(ncid,dmeanF(1:fileSizeG(indexG)),
     $        varname='BackgroundVector',
     $        record_no=record_no,status='single_record')

         if (mapSizeG(indexG) /= size(dmean)) then
            print *,'err[bkg_io_module::getNcCov]: ',
     $           ' map size not comfortable with dmean'
            call errorHalt(1)
         end if
         dmean = dmeanF(mapFromG(1:mapSizeG(indexG),indexG))
         deallocate(dmeanF)
      endif

      if (present(un_atm)) then
         indexG = searchID(fidG,ncid)
         indexC = searchID(fidC,ncid)
         if (.not. stateIndexEnabled(indexG)) then
            print *,'err[bkg_io_module::getNcCov]  ',
     $           ' stateIndex not enabled for un_atm.'
            call errorHalt(1)
         end if

         if (.not. allocated(un_atmF)) 
     $        allocate(un_atmF(MAXPAR,MAXPAR))
         call readNcdfData(ncid,
     $        un_atmF(1:fileSizeG(indexG),1:fileSizeC(indexC)),
     $        varname='TransMatrix',
     $        record_no=record_no,status='single_record')

         if (mapSizeG(indexG) /= size(un_atm,dim=1) .or.
     $        mapSizeC(indexC) /= size(un_atm,dim=2)) then
            print *,'err[bkg_io_module::getNcCov]  ',
     $           ' map size not comfortable with un_atm'
            call errorHalt(1)
         end if

         un_atm = un_atmF(mapFromG(1:mapSizeG(indexG),indexG),
     $        mapFromC(1:mapSizeC(indexC),indexC))
         deallocate(un_atmF)
      endif

      if (present(vn_atm)) then
         indexG = searchID(fidG,ncid)
         indexC = searchID(fidC,ncid)
         if (mapSizeG(indexG) /= size(vn_atm,dim=2) .or.
     $        mapSizeC(indexC) /= size(vn_atm,dim=1)) then
            print *,'err[bkg_io_module::getNcCov]  ',
     $           ' map size not comfortable with vn_atm'
            call errorHalt(1)
         end if
         if (.not. stateIndexEnabled(indexG)) then
            print *,'err[bkg_io_module::getNcCov]  ',
     $           ' stateIndex not enabled for vn_atm.'
            call errorHalt(1)
         end if
         if (.not. allocated(vn_atmF)) allocate(vn_atmF(MAXPAR,MAXPAR))
         if (fileEOF(indexG)) then
            if (.not. allocated(un_atmF)) 
     $          allocate(un_atmF(MAXPAR,MAXPAR))
            call readNcdfData(ncid,
     $           un_atmF(1:fileSizeG(indexG),1:fileSizeC(indexC)),
     $           varname='TransMatrix',
     $           record_no=record_no,status='single_record')
            vn_atmF=TRANSPOSE(un_atmF)
            deallocate(un_atmF)
         else
            call readNcdfData(ncid,
     $           vn_atmF(1:fileSizeC(indexC),1:fileSizeG(indexG)),
     $           varname='InvTransMatrix',
     $           record_no=record_no,status='single_record')         
         endif
         vn_atm = vn_atmF(mapFromC(1:mapSizeC(indexC),indexC),
     $        mapFromG(1:mapSizeG(indexG),indexG))
         deallocate(vn_atmF)
      endif

      if (present(cov)) then
         indexG = searchID(fidG,ncid)
         if (.not. stateIndexEnabled(indexG)) then
            print *,'err[bkg_io_module::getNcCov]: ',
     $           ' stateIndex not enabled for cov'
            call errorHalt(1)
         end if

         if (.not. allocated(covF)) 
     $        allocate(covF(MAXPAR,MAXPAR))
         call readNcdfData(ncid,
     $        covF(1:fileSizeG(indexG),1:fileSizeG(indexG)),
     $        varname='CovMatrix',
     $        record_no=record_no,status='single_record')
         if (mapSizeG(indexG) /= size(cov,dim=1) .or.
     $        mapSizeG(indexG) /= size(cov,dim=2)) then
            print *,'err[bkg_io_module::getNcCov]: ',
     $           'map size not comfortable with cov'
            call errorHalt(1)
         end if

         cov = covF(mapFromG(1:mapSizeG(indexG),indexG),
     $        mapFromG(1:mapSizeG(indexG),indexG))
         deallocate(covF)
      endif

      if (present(dmEmMw))
     $     call readNcdfData(ncid,dmEmMw,varname='EmMwBackground',
     $     record_no=record_no,status='single_record')

      if (present(un_EmMw))
     $     call readNcdfData(ncid,un_EmMw,varname='EmMwTransMatrix',
     $     record_no=record_no,status='single_record')

      if (present(vn_EmMw)) then
         indexG = searchID(fidG,ncid)
         if (fileEOF(indexG)) then
            if (.not. allocated(un_EmMwF)) 
     $       allocate(un_EmMwF(size(vn_EmMw,dim=2),size(vn_EmMw,dim=1)))
            call readNcdfData(ncid,un_EmMw,varname='EmMwTransMatrix',
     $         record_no=record_no,status='single_record')
            vn_EmMw=TRANSPOSE(un_EmMwF)
            deallocate(un_EmMwF)
         else
            call readNcdfData(ncid,vn_EmMw,varname='EmMwInvTransMatrix',
     $      record_no=record_no,status='single_record')
         endif
      endif

      if (present(covEmMw))
     $     call readNcdfData(ncid,covEmMw,varname='EmMwCovMatrix',
     $     record_no=record_no,status='single_record')

      if (present(dmEmIR))
     $     call readNcdfData(ncid,dmEmIR,
     $     varname='EmIRBackground',
     $     record_no=record_no,status='single_record')

      if (present(un_EmIR))
     $     call readNcdfData(ncid,un_EmIR,
     $     varname='EmIRTransMatrix',
     $     record_no=record_no,status='single_record')

      if (present(vn_EmIR)) then
         indexG = searchID(fidG,ncid)
         if (fileEOF(indexG)) then
            if (.not. allocated(un_EmIRF)) 
     $       allocate(un_EmIRF(size(vn_EmIR,dim=2),size(vn_EmIR,dim=1)))
            call readNcdfData(ncid,un_EmIR,varname='EmIRTransMatrix',
     $         record_no=record_no,status='single_record')
            vn_EmIR=TRANSPOSE(un_EmIRF)
            deallocate(un_EmIRF)
         else
            call readNcdfData(ncid,vn_EmIR,varname='EmIRInvTransMatrix',
     $      record_no=record_no,status='single_record')
         endif
      endif

      if (present(covEmIR))
     $     call readNcdfData(ncid,covEmIR,
     $     varname='EmIRCovMatrix',
     $     record_no=record_no,status='single_record')

      if (present(qc))
     $     call readNcdfData(ncid,qc,varname='qc',
     $     record_no=record_no,status='single_record')

      return

      end subroutine getNcCov

      subroutine putNcCov(ncid, dmean, cov, un_atm, vn_atm,
     &     dmEmMw, covEmMw, un_EmMw, vn_EmMw, 
     &     dmEmIR, covEmIR, un_EmIR, vn_EmIR, 
     $     TypeFlag, bkgID, date, qc)
c--
      integer, intent(in) :: ncid
      real, dimension (:), intent (in), optional :: dmean, dmEmMw, 
     $     dmEmIR
      real, dimension (:,:), intent (in), optional :: cov,covEmMw,
     $     covEmIR
      real, dimension (:,:), intent (in), optional :: un_atm,un_EmMw,
     $     un_EmIR
      real, dimension (:,:), intent (in), optional :: vn_atm,vn_EmMw,
     $     vn_EmIR
      integer, intent (in), optional :: TypeFlag
      character (len=*), intent (in), optional :: bkgID
      integer, dimension (:), intent (in), optional :: date, qc

      integer :: indexG

      if (present(dmean))
     $     call writeNcdfData(ncid,dmean,
     $     varname='BackgroundVector',
     &     varlenName=(/'NParG','nbkg '/),
     &     varLongName='Mean Background Vector',
     $     varUnit='Temps(Kelvin), H2O(g/g) in mol(1), psfc(mb)',
     $     status='append')

      if (present(cov))
     $     call writeNcdfData(ncid,cov,
     $     varname='CovMatrix',
     $     varLenName=(/'NParG','NParG','nbkg '/),
     &     varLongName='Covariance Matrix',
     $     status='append')

      if (present(un_atm))
     $     call writeNcdfData(ncid,un_atm,
     $     varname='TransMatrix',
     &     varlenName=(/'NParG','ntran','nbkg '/),
     &     varLongName='Transformation Matrix',
     $     status='append')

      if (present(vn_atm)) then
         indexG = searchID(fidG,ncid)
         if (fileEOF(indexG)  .and. present(un_atm)) then
            print *,'warning[bkg_io_module::putNcCov]:  vn_atm not ',
     $           'written because un_atm is written ',
     $           'and is the transpose'
         else
           call writeNcdfData(ncid,vn_atm,
     $     varname='InvTransMatrix',
     &     varlenName=(/'ntran','NParG','nbkg '/),
     &     varLongName='Inverse Transformation Matrix',
     $     status='append')
         endif
      endif

      if (present(dmEmMw))
     $     call writeNcdfData(ncid,dmEmMw,
     $     varname='EmMwBackground',
     &     varlenName=(/'NGEmMw','nbkg  '/),
     &     varLongName='Mean MW emissivity Background Vector',
     $     status='append')

      if (present(covEmMw))
     $     call writeNcdfData(ncid,covEmMw,
     $     varname='EmMwCovMatrix',
     $     varLenName=(/'NGEmMw','NGEmMw','nbkg  '/),
     &     varLongName='MW Emissivity Covariance Matrix',
     $     status='append')

      if (present(un_EmMw))
     $     call writeNcdfData(ncid,un_EmMw,
     $     varname='EmMwTransMatrix',
     &     varlenName=(/'NGEmMw   ','ntranEmMw','nbkg     '/),
     &     varLongName='MW Emissivity Transformation Matrix',
     $     status='append')

      if (present(vn_EmMw)) then
         indexG = searchID(fidG,ncid)
         if (fileEOF(indexG)  .and. present(un_EmMw)) then
            print *,'warning[bkg_io_module::putNcCov]:  vn_EmMw not ',
     $           'written because un_EmMw is written ',
     $           'and is the transpose'
         else
           call writeNcdfData(ncid,vn_EmMw,
     $     varname='EmMwInvTransMatrix',
     &     varlenName=(/'ntranEmMw','NGEmMw   ','nbkg     '/),
     &     varLongName='MW Emissivity Inverse Transformation Matrix',
     $     status='append')
         endif
      endif

      if (present(dmEmIR))
     $     call writeNcdfData(ncid,dmEmIR,
     $     varname='EmIRBackground',
     &     varlenName=(/'NGEmIR','nbkg  '/),
     &     varLongName='Mean IR emissivity and reflectivity '
     $     //'Background Vector',
     $     status='append')

      if (present(covEmIR))
     $     call writeNcdfData(ncid,covEmIR,
     $     varname='EmIRCovMatrix',
     $     varLenName=(/'NGEmIR','NGEmIR','nbkg  '/),
     &     varLongName='IR Emissivity and Reflectivity '
     $     //'Covariance Matrix',
     $     status='append')

      if (present(un_EmIR))
     $     call writeNcdfData(ncid,un_EmIR,
     $     varname='EmIRTransMatrix',
     &     varlenName=(/'NGEmIR   ','ntranEmIR','nbkg     '/),
     &     varLongName='IR Emissivity and Reflectivity '
     $     //'Transformation Matrix',
     $     status='append')

      if (present(vn_EmIR)) then
         indexG = searchID(fidG,ncid)
         if (fileEOF(indexG)  .and. present(un_EmIR)) then
            print *,'warning[bkg_io_module::putNcCov]:  vn_EmIR not ',
     $           'written because un_EmIR is written ',
     $           'and is the transpose'
         else
           call writeNcdfData(ncid,vn_EmIR,
     $     varname='EmIRInvTransMatrix',
     &     varlenName=(/'ntranEmIR','NGEmIR   ','nbkg     '/),
     &     varLongName='IR Emissivity Inverse Transformation Matrix',
     $     status='append')
         endif
      endif

      if (present(TypeFlag))
     $     call writeNcdfData(ncid,TypeFlag,
     $     varname='TypeFlag',
     $     varLenName='nbkg',
     $     varLongName='Background Type',
     $     status='append')

      if (present(bkgID))
     $     call writeNcdfData(ncid,bkgID,
     $     varname='backgrdID',
     $     strLenName='bkgIdLen',
     $     varLenName='nbkg',
     &     varLongName='Background Identifier',
     $     status='append')

      if (present(date))
     $     call writeNcdfData(ncid,date,
     $     varname='date',
     &     varlenName=(/'nTime','nbkg '/),
     $     varLongname='Date/time tag of background',
     $     status='append')

      if (present(qc))
     $     call writeNcdfData(ncid,qc,
     $     varname='qc',
     &     varlenName=(/'nQC ','nbkg'/),
     $     varLongname='Quality Control Array',
     $     status='append')

      return
      end subroutine putNcCov

      subroutine closeNcFile(ncid,msg)
      character (len=*), intent(in), optional :: msg
      integer, intent (in) :: ncid

      call closeNcdfFile(ncid, msg)

      end subroutine closeNcFile

      end module bkg_io_module
