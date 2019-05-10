!********************************************************************
! $Id$
!********************************************************************
subroutine getNWPprof(latFov,lonFov,ncidmrf,nfil,fracTime,profGrid, &
     latGrid,lonGrid,plandGrid,elevGrid,wt,molID,NGgrid,IGgrid,qc,NWPres)
  !********************************************************************
  ! Subroutine to get 4 grid profiles surrounding a given lat/lon
  !
  ! Copyright, AER, Inc. 2004
  ! Developed at AER, Inc., 2004
  !********************************************************************

  use StateIndexModule
  use scene_io_module
  implicit none

  real, parameter :: dx_hfdeg=0.5        ! grid resolution for 0.5 deg
  real, parameter :: xming_hfdeg=0.0,xmaxg_hfdeg=359.5
  real, parameter :: yming_hfdeg=-90.0,ymaxg_hfdeg=90.0
  real, parameter :: dx_1deg=1.0         ! grid resolution for 1   deg
  real, parameter :: xming_1deg=0.0, xmaxg_1deg=359.0
  real, parameter :: yming_1deg=-90.0, ymaxg_1deg=90.0
  real            :: dx,xming,xmaxg,yming,ymaxg
  real, intent (in):: latFov,lonFov
  real :: lonFov_save,xpos,ypos
  integer :: ilonw,ilone,ilats,ilatn,lrecw,lrecs
  real :: xx,x1,y,y1,x1y1,xy1,xy,x1y
  integer :: nf1,nf2,np
  integer, intent(in) :: nfil
  real, intent(in) :: fracTime
  real, dimension(:,:), intent(inout) :: profGrid
  real, dimension(:), intent(inout) :: latGrid
  real, dimension(:), intent(inout) :: lonGrid
  real, dimension(:), intent(inout) :: plandGrid
  real, dimension(:), intent(inout) :: elevGrid
  real, dimension(:,:), allocatable :: profGrid_t1, profGrid_t2
  real, dimension(size(latGrid)) :: latGrid_t1, latGrid_t2
  real, dimension(size(lonGrid)) :: lonGrid_t1, lonGrid_t2
  real, dimension(size(plandGrid)) :: plandGrid_t1, plandGrid_t2
  real, dimension(size(elevGrid)) :: elevGrid_t1, elevGrid_t2
  real, dimension(:), intent(inout) :: wt
  integer, dimension(:), intent(in) :: molID
  integer, dimension(:), intent(in) :: ncidmrf
  integer, dimension(size(wt)) :: irecp
  type(StateIndex_t), intent(in) :: NGgrid, IGgrid
  integer :: iMol
  integer :: NParGgrid   !local state vector length
  integer :: nGrids
  integer, dimension(:), intent(inout), optional :: qc
  real, intent(in), optional :: NWPres
  
  if (present(NWPres)) then 
     dx=NWPres
  else
     dx=1.0
  endif

  if (dx == dx_hfdeg) then
     xming=xming_hfdeg
     xmaxg=xmaxg_hfdeg
     yming=yming_hfdeg
     ymaxg=ymaxg_hfdeg
  else
     xming=xming_1deg
     xmaxg=xmaxg_1deg
     yming=yming_1deg
     ymaxg=ymaxg_1deg
  endif

  nGrids = size(wt)
  lonFov_save=lonFov
  if(lonFov_save.lt.0.0) lonFov_save=lonFov_save+360
  !      
  !*Calc pos in file of recs east and west of point
  !        
  xpos=(lonFov_save-xming)/dx+1
  ilonw=xpos
  ilone=ilonw+1
  !
  !*Calc pos in file of recs north and south of point
  !
  ypos=(latFov-yming)/dx
  ilats=ypos
  ilatn=ilats+1

  lrecw=nint((xmaxg-xming)/dx+1.0)
  lrecs=nint((ymaxg-yming)/dx+1.0)

  if(ilone > lrecw)ilone=1

  xx=xpos-ilonw
  x1=1.0-xx
  y=ypos-ilats
  y1=1.0-y

  x1y1=x1*y1
  xy1=xx*y1
  xy=xx*y
  x1y=x1*y

  ! Store weights in an array to be passed out of the subroutine. The order 
  ! is important for the adjSfcPresAlt subroutine, make sure it is consistent.

  wt(1) = x1y
  wt(2) = xy
  wt(3) = xy1
  wt(4) = x1y1
  !
  !* Determine record numbers of the 4 profiles to be read in from the nwp file.
  !
  irecp(1)=ilatn*lrecw+ilonw  ! NW corner
  irecp(2)=ilatn*lrecw+ilone  ! NE corner
  irecp(3)=ilats*lrecw+ilone  ! SE corner
  irecp(4)=ilats*lrecw+ilonw  ! SW corner

  if(fracTime.lt.1.e-12)then
     nf1=nfil
     nf2=nfil
  else
     nf1=nfil
     nf2=nfil+1
  endif

  NParGgrid = getVectorLength(NGgrid)
  if (allocated(profGrid_t1)) deallocate(profGrid_t1)
  allocate(profGrid_t1(NParGgrid,nGrids))
  if (allocated(profGrid_t2)) deallocate(profGrid_t2)
  allocate(profGrid_t2(NParGgrid,nGrids))
  do np=1,nGrids
     call getScene(ncid=ncidmrf(nf1),irec=irecp(np), &
          lat=latGrid_t1(np),lon=lonGrid_t1(np), &
          surfalt=elevGrid_t1(np),pland=plandGrid_t1(np), &
          x=profGrid_t1(1:NParGgrid,np),qc=qc)
     call getScene(ncid=ncidmrf(nf2),irec=irecp(np), &
          lat=latGrid_t2(np),lon=lonGrid_t2(np), &
          surfalt=elevGrid_t2(np),pland=plandGrid_t2(np), &
          x=profGrid_t2(1:NParGgrid,np),qc=qc)
  enddo
  !-- Re-arrangement and time interpolation of NCEP fields
  !-- (indexing inverted for use in adjSfcPresAlt).

  if (NGgrid%temp > 0) &
       profGrid(IGgrid%temp:IGgrid%temp+NGgrid%temp-1,:) = &
       (1.0-fracTime)*profGrid_t1(IGgrid%temp:IGgrid%temp+NGgrid%temp-1,:)+ &
       fracTime*profGrid_t2(IGgrid%temp:IGgrid%temp+NGgrid%temp-1,:)
  
  do iMol=1,getNMol(molId)
     if (NGgrid%mol(iMol) > 0) &
          profGrid(IGgrid%mol(iMol):IGgrid%mol(iMol)+NGgrid%mol(iMol)-1,:) = &
          (1.0-fracTime)* &
          profGrid_t1(IGgrid%mol(iMol):IGgrid%mol(iMol)+NGgrid%mol(iMol)-1,:)+ &
          fracTime* &
          profGrid_t2(IGgrid%mol(iMol):IGgrid%mol(iMol)+NGgrid%mol(iMol)-1,:)
  enddo

  if (NGgrid%Tskin > 0) &
       profGrid(IGgrid%Tskin,:) = &
       (1.0-fracTime)*profGrid_t1(IGgrid%Tskin,:)+ &
       fracTime*profGrid_t2(IGgrid%Tskin,:)

  if (NGgrid%Psfc > 0) &
       profGrid(IGgrid%Psfc,:) = &
       (1.0-fracTime)*profGrid_t1(IGgrid%Psfc,:)+ &
       fracTime*profGrid_t2(IGgrid%Psfc,:)

  if (NGgrid%cldliq > 0) &
       profGrid(IGgrid%cldliq:IGgrid%cldliq+NGgrid%cldliq-1,:) = &
       (1.0-fracTime)*profGrid_t1(IGgrid%cldliq:IGgrid%cldliq+ &
       NGgrid%cldliq-1,:)+ &
       fracTime*profGrid_t2(IGgrid%cldliq:IGgrid%cldliq+NGgrid%cldliq-1,:)

  if (NGgrid%cldice > 0) &
       profGrid(IGgrid%cldice:IGgrid%cldice+NGgrid%cldice-1,:) = &
       (1.0-fracTime)*profGrid_t1(IGgrid%cldice:IGgrid%cldice+ &
       NGgrid%cldice-1,:)+ &
       fracTime*profGrid_t2(IGgrid%cldice:IGgrid%cldice+NGgrid%cldice-1,:)

  if (NGgrid%wind > 0) &
       profGrid(IGgrid%wind:IGgrid%wind+NGgrid%wind-1,:) = &
       (1.0-fracTime)*profGrid_t1(IGgrid%wind:IGgrid%wind+NGgrid%wind-1,:)+ &
       fracTime*profGrid_t2(IGgrid%wind:IGgrid%wind+NGgrid%wind-1,:)

  latGrid(:) = (1.0-fracTime)*latGrid_t1(:)+ fracTime*latGrid_t2(:)
  lonGrid(:) = (1.0-fracTime)*lonGrid_t1(:)+ fracTime*lonGrid_t2(:)
  plandGrid(:) = (1.0-fracTime)*plandGrid_t1(:)+ fracTime*plandGrid_t2(:)
  elevGrid(:) = (1.0-fracTime)*elevGrid_t1(:)+ fracTime*elevGrid_t2(:)
  return
end subroutine getNWPprof
