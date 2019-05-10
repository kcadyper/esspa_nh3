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

program InterpNWPprof

! <f90Main>*************************************************************
!
! NAME:
!
!   InterpNWPprof
!
! PURPOSE:
!
!   Program to interpolate NWP vertical profiles at cluster centers.
!
! INCLUDES:
!
!   None
!
!*************************************************************</f90Main>

!*********************************************************************
!     Program to interpolate NWP vertical profiles at cluster centers. 
!     Syntax: InterpNWPprof < run.in
!*********************************************************************
  use StateIndexModule
  use gridProf
  use scene_io_module
  use ClusterFileIO, only : getDimCluster,getClusterData
  use interpNWPtools
  use constants,    only : MISSING_REAL
  implicit none

  character (len=256) :: NWPFile,clusterFile,sceneFile,auxFile,regrFile
  logical :: auxFlag=.true.

  integer :: nSamples,nlpres,nlwv   ! NWP file dimension
  logical :: invertPresGrid
  real, dimension(:),     allocatable :: presGrid  ! NWP pressure grid in mb
  real, dimension(:),     allocatable :: tskin,tair,psfc,pland,elev 
  real, dimension(:,:),   allocatable :: temp,h2o,rh

  integer :: nClusters,nPixels   ! cluster file dimension
  integer :: nChannels
  integer, dimension(:), allocatable :: nmembers,istart,memberList
  integer, dimension(:), allocatable :: fullMemberList

  real                        :: qCutOffPres=-9.  ! negative => gridProf uses default
  real, dimension(:), pointer :: pref   ! scene file pressure grid
  integer                     :: nlev   ! scene file pressure grid levels 
  type(StateIndex_t)          :: IGgrid, NGgrid, IGscene, NGscene, IGaux, NGaux
  integer                     :: NParGgrid, NParGscene,NParGaux
  integer, dimension(maxMol)  :: MolIDgrid,MolIDScene,MolIDAux

  real                                 :: plandClust,elevClust
  integer                              :: nQcScene = 1, nQcClust=1
  integer, dimension(:),   allocatable :: qcClust
  integer, dimension(:,:), allocatable :: qcScene
  real,    dimension(:),   allocatable :: profClust
  real,    dimension(:,:), allocatable :: profScene,profAux
  real,    dimension(:),   allocatable :: plandaux,elevaux
  integer                              :: ncidScene, ncidAux, ncidClust
  integer                              :: clust
  logical                              :: clustAvgFlag = .true.,debug
  integer                              :: nprofs
  namelist /enviofile/ NWPFile,clusterFile,regrFile,sceneFile,auxFile,auxFlag,clustAvgFlag,nprofs,debug

  read(*,enviofile)

  ! read NWP file dimension and cluster file dimension
  call getNWPdim(NWPFile,nSamples=nSamples,nlpres=nlpres,nlwv=nlwv)
  if (clustAvgFlag) then 
     call getDimCluster(ncidClust,clusterFile,nChannels,nClusters,nPixels)
  else
     nClusters = nSamples
     nPixels = nSamples
  endif
  if (nSamples /= nPixels) then
     print*,"Total number of members in NWP file and Cluster file don't agree:",nSamples, nPixels
     stop
  endif


  if(nprofs < 1) nprofs=nClusters
  nClusters=min(nClusters,nprofs)

  ! load cluster data     
  if (allocated(nmembers)) deallocate(nmembers)
  allocate(nmembers(nClusters))
  if (allocated(istart)) deallocate(istart)
  allocate(istart(nClusters))
  if (allocated(memberList)) deallocate(memberList)
  allocate(memberList(nPixels))
  if (allocated(fullMemberList)) deallocate(fullMemberList)
  allocate(fullMemberList(nPixels))
  if (clustAvgFlag) then
     call getClusterData(ncidClust,clusterFile,memberList,istart,nmembers,fullMemberList)
  else
     nmembers=1
     istart=(/(clust,clust=0,nClusters-1)/)
     memberList(1:nClusters)=istart(1:nClusters)
  endif

  ! load NWP data  
  if (allocated(presGrid)) deallocate(presGrid) 
  allocate(presGrid(nlpres))
  if (allocated(temp)) deallocate(temp) 
  allocate(temp(nlpres,nSamples))
  if (allocated(rh)) deallocate(rh) 
  allocate(rh(nlwv,nSamples))
  if (allocated(tskin)) deallocate(tskin) 
  allocate(tskin(nSamples))
  if (allocated(tair)) deallocate(tair) 
  allocate(tair(nSamples))
  if (allocated(psfc)) deallocate(psfc) 
  allocate(psfc(nSamples))
  call getNWPdata(NWPFile,presGrid=presGrid,tskin=tskin,tair=tair,psfc=psfc,temp=temp,rh=rh,nSamples=nSamples)
  if (auxFlag) then
     if (allocated(pland)) deallocate(pland)
     allocate(pland(nSamples))
     if (allocated(elev)) deallocate(elev)
     allocate(elev(nSamples))
     call getNWPdata(NWPFile,pland=pland,elev=elev,nSamples=nSamples)
     where (pland > 0.5)
        tskin=tair
     endwhere
  endif

  ! convert relative humidity to water vapor mixing ratio                                                                 
  if (allocated(h2o)) deallocate(h2o)
  allocate(h2o(nlwv,nSamples))
  call computeH2oMixr(temp,presGrid,rh,h2o)
  if (allocated(rh)) deallocate(rh)

  ! load vertical regrid coefficient
  if (associated(pref)) deallocate(pref)
  call setGridProf(regrFile,pref,qCutOffPres=qCutOffPres)
  nlev=size(pref)
 
  !Arrange NWP data array in vertical profile format
  invertPresGrid=(presGrid(1) > presGrid(nlpres))
  if (auxFlag) then
     call setStateIndex(nlpres,nlwv,nlev,IGgrid,NGgrid,molIDgrid,NParGgrid,&
          IGscene,NGscene,molIDScene,NParGscene,IGaux=IGaux,NGaux=NGaux,molIDaux=molIDaux,NParGaux=NParGaux)
  else
     call setStateIndex(nlpres,nlwv,nlev,IGgrid,NGgrid,molIDgrid,NParGgrid,&
          IGscene,NGscene,molIDScene,NParGscene)
  endif

  if (allocated(profClust)) deallocate(profClust)
  allocate(profClust(NParGgrid))
  if (allocated(qcClust)) deallocate(qcClust)
  allocate(qcClust(nQcClust))

  if (allocated(profScene)) deallocate(profScene)
  allocate(profScene(NParGscene,nClusters))    
  if (allocated(qcScene)) deallocate(qcScene)
  allocate(qcScene(nQcScene,nClusters))

  qcScene=0
  profScene=MISSING_REAL
  ! vertical grid interpolation
  CLUSTERloop: do clust=1,nClusters
     if (debug .and. (mod(clust,10000) ==0)) print*,'clust= ',clust
     call avgClusterProf(istart(clust),nmembers(clust),memberList,nSamples,IGgrid,NGgrid,molIDgrid, &
          temp=temp,h2o=h2o,tskin=tskin,psfc=psfc,prof=profClust,invertGrid=invertPresGrid,qc=qcClust)
     if (btest(qcClust(1),0)) then
        qcScene(1,clust)=ibset(qcScene(1,clust),0)
        profScene(:,clust)=MISSING_REAL
     else
        call regridProf(profClust,IGgrid,NGgrid,molIDgrid,profscene(:,clust),IGscene,NGscene)
        call checkCldliq(profscene(:,clust),IGscene,NGscene,qcScene(:,clust))
     endif
  enddo CLUSTERloop
  
  ! output sceneFile
  call openScene(ncid=ncidScene,file=sceneFile,Nlevel=nlev,pressure=pref,&
       MolID=MolIDscene,IG=IGscene,NG=NGscene,status='new')
  do clust=1,nClusters
     CALL putScene(ncid=ncidScene,x=profscene(:,clust),qc=qcScene(:,clust))      
  enddo
  call closeScene(ncidScene)
  
  !writes out the auxiliary data
  if (auxFlag) then
     if (allocated(profAux))  deallocate(profAux)
     allocate(profAux(NParGaux,nClusters))
     if (allocated(plandAux)) deallocate(plandAux)
     allocate(plandaux(nClusters))
     if (allocated(elevAux))  deallocate(elevAux)
     allocate(elevAux(nClusters))
     do clust=1,nClusters
     if (debug .and. (mod(clust,10000) ==0)) print*,'clust= ',clust
        call avgClusterProf(istart(clust),nmembers(clust),memberList,nSamples,IGgrid,NGgrid,molIDgrid, &
             pland=pland,elev=elev,plandAvg=plandClust,elevAvg=elevClust,qc=qcClust)
        plandAux(clust)=plandClust
        elevAux(clust)=elevClust
        if  (btest(qcClust(1),0)) then
           qcScene(1,clust)=ibset(qcScene(1,clust),0)
           profAux(:,clust)=MISSING_REAL
        else
           profAux(:,clust)=profscene(IGscene%psfc,clust)
        endif
     enddo
     ! output auxFile
     CALL openScene(ncid=ncidAux,file=auxfile,MolID=molIDaux,&
          IG=IGaux,NG=NGaux,status='new')
     do clust=1,nClusters
        CALL putScene(ncid=ncidAux,x=profAux(:,clust),qc=qcScene(:,clust),&
             pland=plandAux(clust),surfalt=elevAux(clust))
     enddo
     call closeScene(ncidAux)
     if (allocated(profAux))   deallocate(profAux)
     if (allocated(plandAux))  deallocate(plandAux)
     if (allocated(elevAux))   deallocate(elevAux)
  endif

  if (allocated(presGrid))  deallocate(presGrid) 
  if (allocated(temp))      deallocate(temp) 
  if (allocated(h2o))       deallocate(h2o) 
  if (allocated(tskin))     deallocate(tskin) 
  if (allocated(tair))      deallocate(tair) 
  if (allocated(psfc))      deallocate(psfc) 
  if (allocated(pland))     deallocate(pland) 
  if (allocated(elev))      deallocate(elev) 
  if (allocated(nmembers))  deallocate(nmembers) 
  if (allocated(istart))    deallocate(istart) 
  if (allocated(memberList))deallocate(memberList) 
  if (allocated(fullMemberList))deallocate(fullMemberList) 
  if (associated(pref))     deallocate(pref)
  if (allocated(profClust)) deallocate(profClust)
  if (allocated(qcClust))   deallocate(qcClust)
  if (allocated(profScene)) deallocate(profScene)
  if (allocated(qcScene))   deallocate(qcScene)

end program InterpNWPprof
  

