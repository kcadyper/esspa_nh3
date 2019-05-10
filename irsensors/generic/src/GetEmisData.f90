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

program GetEmisData

! <f90Main>*************************************************************
!
! NAME:
!
!   GetEmisData
!
! PURPOSE:
!
!   Program to access emissivity from a database
!
! INCLUDES:
!
!   None
!
!*************************************************************</f90Main>

!******************************************************
!     Program to access emissivity from a database
!******************************************************
  use scene_io_module
  use bkg_io_module
  use ncdf_module
  use AccessSfcGrid
  use constants, ONLY : MISSING_REAL
  use RadFile, ONLY: getDimRad,getRadData
  use ClusterFileIO, ONLY : getDimCluster,getClusterData

  implicit none

  real,    parameter :: dLim=0.001 ! Freq match threshold (fractional)
  real,    parameter :: extendFactor=2. ! for avg range
  real,    parameter :: latPolar=65. ! tolerate missing IR emis above this lat
  real,    parameter :: emisPolarStandard=0.985 ! see Dozier & Warren
                                                ! J. Water Res. 1982
  real,    dimension(:),   pointer      :: wvndb,wvndbCov
  integer                               :: nchdb,nchdbCov
  integer                               :: nchanSrc
  real,    dimension(:,:), allocatable  :: cov
  real,    dimension(:),   allocatable  :: emOut
  real,    dimension(:),   allocatable  :: emSum,emLoc
  real,    dimension(:,:), allocatable  :: emAll
  real,    dimension(:),   allocatable  :: emissdb,covdb
  integer, dimension(:),   allocatable  :: qcSrc,qcEmis
  integer                               :: nQcSrc=0
  integer, dimension(:,:), allocatable  :: iflg
  real,    dimension(:),   allocatable  :: latAll,lonAll
  real,    dimension(:),   allocatable  :: clatAll,clonAll
  integer, dimension(:),   allocatable  :: startIndices,nMembers
  integer, dimension(:),   allocatable  :: cloudyMembers,memberList

  integer :: i,irec,ncidSrc,ncidCls,nfov,ifov,ncidEmiss,time(6)
  integer :: nchanCls,nClust,nMembersTotal
  integer :: iMem,nSum
  real    :: lat,lon
  logical :: qcIn=.TRUE.
  integer :: nExtend=0,nFailed=0,nPolar=0
  real    :: acceptRangeExtended

  !---Control variables in namelist:
  integer :: nprofs              ! Number of locations to process (-1 = do all)
  real    :: acceptRange=1.2     ! Averaging range, as fraction of grid spacing
  logical :: sourceRad=.TRUE.    ! Chan set and locations from rad (vs scn) file?
  logical :: forceAverage=.FALSE.! Always spatial average, even if nearest OK?
  logical :: covIO=.FALSE.       ! Process gridded covariance data?
  logical :: clusterAvg=.FALSE.  ! Average emissivities over clusters before output
  character(len=256) :: emisDBaseBkgFile        ! Emis database
  character(len=256) :: emisDBaseCovFile='none' ! Covariance database
  character(len=256) :: srcFile                 ! Chan set & location source
  character(len=256) :: emisFile                ! Output emis file
  character(len=256) :: clusterFile             ! Cluster locations and map to members
  namelist /EmisInput/nprofs,acceptRange,emisDBaseBkgFile,covIO,clusterAvg, &
       emisDBaseCovFile,srcFile,sourceRad, &
       forceAverage,emisFile,clusterFile
  
  !---Read Namelist Inputs
  print *, 'Reading namelist...'
  read(*,emisinput)

  !---Average range were default average range fails
  acceptRangeExtended=acceptRange*extendFactor

  !---Read header of Source file
  if (sourceRad) then   !---source file is rad file
     call getDimRad(ncidSrc,srcFile,nchanSrc,nfov)
     nQcSrc=0
  else                  !---source file is scene file
     call openScene(ncid=ncidSrc,file=srcFile,nSfcGrid=nchanSrc,nprf=nfov, &
        nqc=nQcSrc)
  end if

  !---Read cluster data
  if (clusterAvg) then
     if (covIO) then
        print *,'err[GetEmisData]: cluster averaging not currently '// &
           'supported with covIO'
        call errorHalt(1)
     endif
     call getDimCluster(ncidCls,clusterFile,nchanCls,nClust,nMembersTotal)
     if (nMembersTotal /= nfov) then
        print *,'err[GetEmisData]: Mismatch nMembersTotal /= nfov:', &
             nMembersTotal,nfov
        call errorHalt(1)
     endif
     allocate (cloudyMembers(nfov),startIndices(nClust), &
                         nMembers(nClust),memberList(nfov))
     allocate (clatAll(nClust),clonAll(nClust))
     call getClusterData(ncidCls,clusterFile,cloudyMembers,startIndices, &
                         nMembers,memberList)
     call getRadData(ncidCls,clusterFile,latAll=clatAll,lonAll=clonAll)
  endif
    
  if(nfov < nprofs)then
     print *,'err[GetEmisData]:  # of profiles in the srcFile ',&
          'less than nprofs:',nfov,nprofs
     call errorHalt(1)
  elseif(nprofs < 0)then
     nprofs = nfov
  endif
  allocate(iflg(nprofs,4))
  !---Allocate memory for QC
  if (nQcSrc <= 0) then
     qcIn=.false.                  ! srcFile lacks QC
     allocate(qcSrc(1),qcEmis(1))
     qcSrc(1)=0                    ! assume all are OK and leave qcSrc fixed
  else
     allocate(qcSrc(nQcSrc),qcEmis(nQcSrc))
  endif
  !---Find the 1st valid time to feed loadSfcGrid
  if (sourceRad) then
     allocate (latAll(nfov),lonAll(nfov))
     call getRadData(ncidSrc,srcFile,latAll=latAll,lonAll=lonAll,time=time)
  else
     do irec = 1, nfov
        call getScene(ncid=ncidSrc,irec=irec,time=time)
        if (time(1) > 0 .and. time(2) > 0 .and. time(3) > 0) exit
     enddo
  end if
  !---Read emissivity Bkg database
  print *, 'Reading emissivity Bkg database...'
  call loadSfcGrid(ncfile=emisDBaseBkgFile,nchan=nchdb,wvn=wvndb, &
       time=time,acceptRange=acceptRange,forceAverage=forceAverage)

  !---allocate memory for emissiv vectors
  allocate(emissdb(nchdb),emOut(nchdb))
  if (clusterAvg) allocate (emAll(nchdb,nprofs),emSum(nchdb),emLoc(nchdb))
  !---Open Output File
  call openCov(ncidEmiss,emisFile,nemir=nchdb, &
          frqEmir=wvndb(1:nchdb),status='new')
  !---Loop over the number of profiles to process
  do ifov=1,nprofs
     if (sourceRad) then
        lat=latAll(ifov)
        lon=lonAll(ifov)
     else
        call getScene(ncid=ncidSrc,irec=ifov,lat=lat,lon=lon)
        if (qcIn) call getScene(ncid=ncidSrc,irec=ifov,qc=qcSrc)
     end if
     if (.not.(btest(qcSrc(1),0))) then
        iflg(ifov,1) = 1
        if(lon.gt.180.)lon=lon-360.
        call setSfcLoc(lat,lon,emissdb,iflg(ifov,2))
        emOut=emissdb
        if (iflg(ifov,2) == 0) then
        ! Try again with larger averaging radius
           call setSfcLoc(lat,lon,emissdb,iflg(ifov,2), &
              acceptRange=acceptRangeExtended)
           if (iflg(ifov,2) == 0) then
              if (abs(lat) >= latPolar) then
                emOut = emisPolarStandard
                nPolar=nPolar+1
              else
                emOut = MISSING_REAL
                nFailed=nFailed+1
              endif
           else
              nExtend=nExtend+1
           endif
        endif
     else
        iflg(ifov,1:2) = 0 ! QC/value depending on the Src QC
        emOut     = MISSING_REAL
     endif
     !---Set the emissiv 
     if (clusterAvg) then
        emAll(:,ifov)=emOut
     else
        call putCov(ncidEmiss,dmEmIR=emOut(1:nchdb))
     endif
  enddo
  !---Close db file
  call destroySfcGrid()

  if (nExtend > 0) print *,nExtend,'locations where extended average was used'
  if (nFailed > 0) print *,nFailed,'locations where extended average failed'
  if (nPolar  > 0) print *,nPolar,'locations where reverted to polar standard'

  if(covIO)then  !--Cov available
     !---Read emissivity Cov database
     print *, 'Reading emissivity Cov database...'
     call loadSfcGrid(ncfile=emisDBaseCovFile,nchan=nchdbCov,wvn=wvndbCov, &
          icov=1,time=time,acceptRange=acceptRange, &
          forceAverage=forceAverage)
     if (nchdbCov /= nchdb) then
        print *,'err[GetEmisData]:  Channel size mismatch; nchdbCov,nchdb:', &
             nchdbCov,nchdb
        call errorHalt(1)    
     elseif (any((abs(wvndbCov-wvndb)/wvndb) > dLim)) then
        print *,'err[GetEmisData]:  Channel set mismatch in Bkg and Cov files;'
        do i=1,nchdb
           write(*,'(i4,2(f10.4,i4))')i,wvndbCov(i),wvndb(i)
        enddo
        call errorHalt(1)    
     endif

     !---allocate memory four cov matrix
     allocate(covdb(nchdb*nchdb))
     !---Loop over the number of profiles to process
     do ifov=1,nprofs
        if (sourceRad) then
           lat=latAll(ifov)
           lon=lonAll(ifov)
        else
           call getScene(ncid=ncidSrc,irec=ifov,lat=lat,lon=lon)
           if (qcIn) call getScene(ncid=ncidSrc,irec=ifov,qc=qcSrc)
        end if
        if (.not.(btest(qcSrc(1),0))) then 
           if(lon.gt.180.)lon=lon-360.
           call setSfcLoc(lat,lon,covdb,iflg(ifov,3))
           if (iflg(ifov,3) == 0) then
           ! Try again with larger averaging radius
              call setSfcLoc(lat,lon,covdb,iflg(ifov,3), &
                 acceptRange=acceptRangeExtended)
              if (iflg(ifov,3) == 0) covdb = MISSING_REAL
           endif
        else
           iflg(ifov,3) = 0
           covdb = MISSING_REAL
        endif
        call putCov(ncidEmiss,covEmIR=reshape(covdb,(/size(wvndb),size(wvndb)/)))
     enddo
     !---Close files and deallocate memory
     call destroySfcGrid()
     deallocate(covdb)
  endif

  !--Prepare QC flags
  if (.not. clusterAvg) then
     do ifov=1,nprofs
        !---Set the emissiv QC/value depending on the Src QC and the iflg's
        qcEmis = 0                                                 ! By default, Emis QC is fine
        if(any(iflg(ifov,1:2) == 0))qcEmis(1) = ibset(qcEmis(1),0) ! bkg-related flags: summary
        if(iflg(ifov,1) == 0) qcEmis(1) = ibset(qcEmis(1),1)       !                    src
        if(iflg(ifov,2) == 0) qcEmis(1) = ibset(qcEmis(1),2)       !                    db
        if(covIO)then                                              ! cov-related flags
           if(iflg(ifov,3) == 0) then
              qcEmis(1) = ibset(qcEmis(1),0)                       ! summary
              qcEmis(1) = ibset(qcEmis(1),3)                       ! cov-specific bit
           endif
        else                                                       ! cov unavailable
           qcEmis(1) = ibset(qcEmis(1),5)
        endif
        qcEmis(1) = ibset(qcEmis(1),6)                             ! mtr unavailable
        call putCov(ncidEmiss,qc=qcEmis)
     enddo
  endif

  !--- Do cluster averaging of emissivities
  if (clusterAvg) then
     do i=1,nClust
        emSum=0.
        nSum=0
        do iMem=startIndices(i)+1,startIndices(i)+nMembers(i)
          emLoc=emAll(:,cloudyMembers(iMem)+1)
          if (any(emLoc /= MISSING_REAL)) then 
             emSum=emSum+emLoc
             nSum=nSum+1
          endif
        enddo
        qcEmis = 0 
        if (nSum > 0) then
           emOut = emSum/real(nSum)
        else
           emOut = MISSING_REAL
           qcEmis(1) = ibset(qcEmis(1),0)
        endif
        call putCov(ncidEmiss,dmEmIR=emOut(1:nchdb),qc=qcEmis)
     enddo
  endif

  !---Close files and deallocate memory
  call destroySfcGrid()
  if (.not. sourceRad) call closeScene(ncid=ncidSrc)
  call closeCov(ncidEmiss)
  deallocate(emissdb,emOut)
  deallocate(qcSrc,qcEmis,iflg)
  if (associated(wvndb)) deallocate(wvndb)
  if (associated(wvndbCov)) deallocate(wvndbCov)
  if (allocated(latAll)) deallocate(latAll,lonAll)
  if (clusterAvg) then
     deallocate (cloudyMembers,startIndices,nMembers,memberList)
     deallocate (clatAll,clonAll,emAll,emSum,emLoc)
  endif
  
end program GetEmisData
