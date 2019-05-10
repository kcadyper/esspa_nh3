c $"$Id$"
      module regr_io_module
c**********************************************************************
c* Module name: regr_io_module
c* Purpose: Handles regression ceof I/O
c* Functions:
c*           openRegr
c*           getRegr
c*           closeRegr
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
c* IC        StateIndex_t  state vector indices in EOF space
c* NC        StateIndex_t  state vector lengths in EOF space
c* IG        StateIndex_t  state vector indices in geophysical space
c* NG        StateIndex_t  state vector lengths in geophysical space
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
      type(StateIndex_t) :: IGF, NGF, ICF, NCF
      
      private
      integer, save :: record_no=0
      public openRegr, closeRegr, getRegr, queryRegr
c-----
      interface openRegr
      module procedure openNcFile
      end interface

      interface queryRegr
      module procedure queryNcRegr
      end interface

      interface getRegr
      module procedure getNcRegr
      end interface

      interface closeRegr
      module procedure closeNcFile
      end interface
c---- 
      CONTAINS
c--   
      subroutine openNcFile(ncid, file, creationDate, lcascd, 
     $     fovsz, nbkg, ifstruc, iftran, NParG, ntran, nLevel,
     $     nchmw, nchir, neof, nemir, pressure, freq, 
     $     frqEmir, polarity, 
     $     icell, eiabkg, bkgIdLen,IC,NC,IG,NG,molID,MolTran,status)

c--   Dummy variables
      integer, intent(inout) :: ncid
      character (len=*), intent(in) :: file
      character (len=*), intent (inout), optional :: creationDate
      integer, intent(inout), optional :: lcascd, nbkg, neof,
     $     ifstruc, iftran, NParG, nLevel, nchmw, nemir, icell, nchir
      real, intent (inout), optional :: fovsz
      real, dimension(:), intent (inout), optional :: pressure
      real, dimension(:), intent (inout), optional :: freq,frqEmir
      integer, dimension(:), intent (inout), optional :: polarity,molID
      integer, dimension(:), intent (inout), optional :: MolTran
      real, dimension (:), intent (inout), optional :: eiabkg
      integer, intent (inout), optional :: bkgIdLen, ntran
      character (len=*), intent(in), optional :: status
      type(StateIndex_t), intent (inout), optional :: IG, NG
      type(StateIndex_t), intent (inout), optional :: IC, NC

c--   Local variables
      character (len=MAX_NAME_LENGTH) :: lstatus
      integer, dimension(2) :: myDim
      integer , dimension(maxMol) :: MolIDF,MolLenF,TMolIDF,TMolLenF
      integer , dimension(maxMol) :: MolTranF
      integer :: nchmwloc,nemirloc, NMol

      if (.not. present(status)) then
         lstatus = 'old'
      else
         lstatus = trim(status)
      endif

      if (trim(lstatus) == 'new' .or. trim(lstatus) == 'NEW') then

         print *,'err[regr_io_module::openNcFile]: '//
     $        ' Cannot create new regr file yet'
         call errorHalt(1)


      else
         call openNcdfFile(ncid,file,status='old')

         call initFileID(ncid=ncid,fidIndex=fidIndexG,
     $        stateIndexEnabled=stateIndexEnabled,IG=IG,NG=NG,fid=fidG)

         if (present(creationDate))
     $        call readNcdfAttr(ncid,attrName='Creationtime',
     &        attr=creationDate)

         if (ncid <= 0) then 
            print *,'warning[regr_io_module::openNcFile]: '
            print *,' Invalid file ID'
            return
         endif

         if (present(nbkg))
     $        nbkg=readNcdfDim(ncid,dimName='nbkg')
         
         if (present(bkgIdLen))
     $        bkgIdLen=readNcdfDim(ncid,dimName='bkgIdLen')

         if (present(ifstruc))
     $        call readNcdfAttr(ncid,
     $        attrName='FileStructure',
     &        attr=ifstruc)
         
         if (present(iftran)) then
            call readNcdfAttr(ncid,attrName='TransFlag',
     &           attr=iftran)
            if (iftran == 1) then
               if (present(IC) .and. present(NC)) then
                  if (.not. stateIndexEnabled(fidIndexG)) then
                     print *,'err[regr_io_module::openNcFile]: '//
     $                    ' stateIndex not enabled in iftran read mode.'
                     call errorHalt(1)
                  end if
                  
                  call readNcdfAttr(ncid,attrName='TTemperature',
     $                 attr=mydim)
                  ICF%Temp = mydim(1)
                  NCF%Temp = mydim(2)
               
                  call readNcdfAttr(ncid,attrName='TTskin',attr=mydim)
                  ICF%Tskin = mydim(1)
                  NCF%Tskin = mydim(2)
               
                  call readNcdfAttr(ncid,attrName='TSurfacePressure',
     $                 attr=mydim)
                  ICF%Psfc=mydim(1)
                  NCF%Psfc=mydim(2)
               
                  call readNcdfAttr(ncid,attrName='TMol',
     $                 attr=TMolIDF)
                  ICF%mol = TMolIDF
                  
                  call readNcdfAttr(ncid,attrName='TMolLen',
     $                 attr=TMolLenF)
                  NCF%mol = TMolLenF

                  call readNcdfAttr(ncid,attrName='TLiqCloud',
     &                 attr=mydim)
                  ICF%cldLiq=mydim(1)
                  NCF%cldLiq=mydim(2)
               
                  call readNcdfAttr(ncid,attrName='TIceCloud',
     &                 attr=mydim)
                  ICF%cldIce=mydim(1)
                  NCF%cldIce=mydim(2)

                  call readNcdfAttr(ncid,attrName='TSurfaceWinds',
     $                 attr=mydim)
                  ICF%Wind=mydim(1)
                  NCF%Wind=mydim(2)

                  if (.not. allocated(mapFromC))
     $                 allocate(mapFromC(MAXPAR,MAXFILE))
                  index = searchID(fidC,ncid)
                  if (index < 0) then 
                     fidIndexC = fidIndexC + 1
                  else
                     fidIndexC = index
                  endif
                  fidC(fidIndexC) = ncid
                  call genMapFrom(IG=IC,NG=NC,IGF=ICF,NGF=NCF,
     $                 MolID=MolIDF,map=mapFromC(:,fidIndexC),
     $                 mapSize=mapSizeC(fidIndexC),
     $                 xFileSize=fileSizeC(fidIndexC))
               else
                  print *,'err[regr_io_module::openNcFile]: '//
     $                 ' incomplete EOF I/O args - no index read'
                  call errorHalt(1)
               endif
            endif
         endif

         if (present(NParG))
     $        NParG=readNcdfDim(ncid,dimName='NParG')

         if (present(ntran))
     $        ntran=readNcdfDim(ncid,dimName='ntran')
         
         if (present(nchir))
     $        nchir=readNcdfDim(ncid,dimName='nChan')
         
         if (present(neof))
     $        neof=readNcdfDim(ncid,dimName='nEof')

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
         
         if (present(nLevel))
     $        call readNcdfAttr(ncid,attrName='nlev',attr=nLevel)

         if (present(pressure))
     $        call readNcdfAttr(ncid,
     $        attrName='StandardPressureGrid',
     &        attr=pressure)

         if (present(freq)) then
            call readNcdfAttr(ncid,attrName='mwfrequencies',
     &           attr=freq)
            call readNcdfAttr(ncid,attr=nchmwloc,attrName='nchmw')
            if (size(freq) /= nchmwloc) then
               print *,'err[regr_io_module::openNcFile]: '//
     $              ' #freq and nchmw mismatched'
               call errorHalt(1)
            end if
         endif

         if (present(polarity)) then
            call readNcdfAttr(ncid,attrName='mwpolarities',
     &           attr=polarity)
            call readNcdfAttr(ncid,attr=nchmwloc,attrName='nchmw')
            if (size(freq) /= nchmwloc) then
               print *,'err[regr_io_module::openNcFile]: '//
     $              ' #mwpol and nchmw mismatched'
               call errorHalt(1)
            end if
         endif

         if (present(nchmw))
     $        call readNcdfAttr(ncid,attr=nchmw,attrName='nchmw')

         if (present(frqEmir)) then
            call readNcdfAttr(ncid,attrName='irfrequencies',
     &           attr=frqEmir)
            call readNcdfAttr(ncid,attr=nemirloc,attrName='nemir')
            if (size(frqEmir) /= nemirloc) then
               print *,'err[regr_io_module::openNcFile]: '//
     $              ' #freq and nemir mismatched'
               call errorHalt(1)
            end if
         endif

         if (present(nemir))
     $        call readNcdfAttr(ncid,attr=nemir,attrName='nemir')

         if (present(eiabkg)) 
     $        call readNcdfAttr(ncid,attrName='eiabkg',attr=eiabkg)

         if (stateIndexEnabled(fidIndexG)) then

            if (.not. allocated(mapFromG)) 
     $           allocate(mapFromG(MAXPAR,MAXFILE))

            call genMapFrom(IG=IG,NG=NG,IGF=IGF,NGF=NGF,
     $           MolID=MolID,map=mapFromG(:,fidIndexG),
     $           mapSize=mapSizeG(fidIndexG),
     $           xFileSize=fileSizeG(fidIndexG))

            if (.not. present(NC)) then
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
            endif               !present(NC)

            if (present(MolTran)) then
               MolTranF=0
               call readNcdfAttr(ncid,attrName='MolTran',
     $              attr=MolTranF)
               NMol = getNMol(NGF)
               MolTran(1:NMol) = MolTranF(MolID(1:NMol))
            endif
         endif                  !stateIndexEnabled
      endif                     !new file

      end subroutine openNcFile
c--
      subroutine queryNcRegr(file, creationDate, lcascd, 
     $     fovsz, nbkg, ifstruc, iftran, NParG, ntran, nLevel,
     $     nchmw, nemir, pressure, freq, frqEmir, polarity, 
     $     icell, eiabkg, bkgIdLen, IC, NC, IG, NG, NMol, molID, 
     $     MolTran)

c--   Dummy variables
      character (len=*), intent(in) :: file
      character (len=*), intent (inout), optional :: creationDate
      integer, intent(inout), optional :: lcascd, nbkg, 
     $     ifstruc, iftran, NParG, nLevel, nchmw, nemir, icell
      real, intent (inout), optional :: fovsz
      real, dimension(:), intent (inout), optional :: pressure
      real, dimension(:), intent (inout), optional :: freq,frqEmir
      integer, dimension(:), intent (inout), optional :: polarity,molID
      integer, dimension(:), intent (inout), optional :: MolTran
      real, dimension (:), intent (inout), optional :: eiabkg
      integer, intent (inout), optional :: bkgIdLen, ntran
      type(StateIndex_t), intent (inout), optional :: IG, NG
      type(StateIndex_t), intent (inout), optional :: IC, NC
      integer , intent(inout), optional :: NMol

c--   Local variables
      integer :: ncid

      if (present(molID)) molID=nullMolID
      if (present(NMol)) NMol=0
      call openNcFile(ncid=ncid, file=file, creationDate=creationDate, 
     $     lcascd=lcascd, fovsz=fovsz, nbkg=nbkg, ifstruc=ifstruc, 
     $     iftran=iftran, NParG=NParG, ntran=ntran, nLevel=nLevel,
     $     nchmw=nchmw, nemir=nemir, pressure=pressure, freq=freq, 
     $     frqEmir=frqEmir, polarity=polarity, 
     $     icell=icell, eiabkg=eiabkg, 
     $     bkgIdLen=bkgIdLen, IC=IC, NC=NC, IG=IG, NG=NG, 
     $     molID=molID, MolTran=MolTran)
      call closeNcFile(ncid)
      end subroutine queryNcRegr
c--         
      subroutine getNcRegr(ncid, irec, dmean, un_atm, cov,
     $     dmEmMw, un_EmMw, covEmMw, dmEmIR, un_EmIR, covEmIR,
     $     un_rad, drmean, ntran, nparg, nchir, neof,
     $     wghtvec, regrcoef,
     $     bkgID, TypeFlag, date)

c-- Dummy arguments
      integer, intent (in) :: ncid
      integer, intent (in), optional :: irec
      real, dimension (:), intent (inout), optional :: drmean
      real, dimension (:), intent (inout), optional :: dmean, dmEmMw,
     $     dmEmIR, wghtvec
      real, dimension (:,:), intent (inout), optional :: un_atm,un_EmMw,
     $     un_EmIR, un_rad, regrcoef
      real, dimension (:,:), intent (inout), optional :: cov,covEmMw,
     $     covEmIR
      integer, intent (inout), optional :: TypeFlag
      integer, intent (inout), optional :: nparg, ntran, nchir, neof
      character (len=*), intent (inout), optional :: bkgID
      integer, dimension (:), intent (inout), optional :: date
c-- Local arguments
      integer :: indexG
      real , dimension (:), allocatable :: dmeanF
      real , dimension (:,:), allocatable :: un_atmF
      real , dimension (:,:), allocatable :: covF

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
            print *,'err[regr_io_module::getNcRegr]: '//
     $           ' stateIndex not enabled for dmean.'
            call errorHalt(1)
         end if
         
         if (.not. allocated(dmeanF)) allocate(dmeanF(MAXPAR))

         call readNcdfData(ncid,dmeanF(1:fileSizeG(indexG)),
     $        varname='BkgProfVector',
     $        record_no=record_no,status='single_record')
         if (mapSizeG(indexG) /= size(dmean)) then
            print *,'err[regr_io_module::getNcRegr]: '//
     $           ' map size not comfortable with dmean'
            call errorHalt(1)
         end if
         dmean = dmeanF(mapFromG(1:mapSizeG(indexG),indexG))
         deallocate(dmeanF)
      endif

      if (present(drmean)) then
         call readNcdfData(ncid,drmean(1:nchir),
     $        varname='BkgRadVector',
     $        record_no=record_no,status='single_record')
      endif

      if (present(wghtvec)) then
         call readNcdfData(ncid,wghtvec(1:nchir),
     $        varname='WeightVector',
     $        record_no=record_no,status='single_record')
      endif

      if (present(regrcoef)) then
         call readNcdfData(ncid,regrcoef(1:neof,1:ntran),
     $        varname='RegrCoef',
     $        record_no=record_no,status='single_record')
      endif

      if (present(un_atm)) then
         if (.not. allocated(un_atmF)) 
     $        allocate(un_atmF(MAXPAR,MAXPAR))
         call readNcdfData(ncid,
     $        un_atm(1:nparg,1:ntran),
     $        varname='EofProf',
     $        record_no=record_no,status='single_record')
         deallocate(un_atmF)
      endif

      if (present(un_rad)) then
         call readNcdfData(ncid,
     $        un_rad(1:nchir,1:neof),
     $        varname='EofRad',
     $        record_no=record_no,status='single_record')
      endif


      if (present(cov)) then
         indexG = searchID(fidG,ncid)
         if (.not. stateIndexEnabled(indexG)) then
            print *,'err[regr_io_module::getNcRegr]: '//
     $           ' stateIndex not enabled for cov.'
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
            print *,'err[regr_io_module::getNcRegr]: '//
     $           ' map size not comfortable with cov'
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

      if (present(covEmIR))
     $     call readNcdfData(ncid,covEmIR,
     $     varname='EmIRCovMatrix',
     $     record_no=record_no,status='single_record')

      return

      end subroutine getNcRegr

      subroutine closeNcFile(ncid,msg)
      character (len=*), intent(in), optional :: msg
      integer, intent (in) :: ncid

      call closeNcdfFile(ncid, msg)

      end subroutine closeNcFile

      end module Regr_io_module
