MODULE FOV2sfcgrd_module
  
  use constants, ONLY : ErthRad
  use mapConvert
  
  implicit none
  private
  public :: wrapLongitude,windowBounds
  
CONTAINS
  
SUBROUTINE wrapLongitude(latTile,lonTile,i,j,nrow,ncol,&
     latorigin,lonorigin,colGridOffset,rowGridOffset,mapRes, &
     iwrap)
  
  ! ####################################################################
  ! allow longitude wrap around if window goes off map (off hemisphere)
  ! assumptions:
  !   0 longitude is origin of map (wrap around occurs at +-180 degrees)
  !   ij2ll_ returns valid lat even off map as a given row's lat is constant
  !   ij2ll_ returns -180 -> +180 for valid or MISSING_REAL for off map
  ! ####################################################################
  
  real    :: latTile,lonTile ! tile center
  real    :: latorigin,lonorigin
  real    :: rowGridOffset,colGridOffset ! ctr of grid
  real    :: mapRes       ! grid resolution [km]
  real    :: col, row ! float col,row (0-based idx) of grid point
  
  real,parameter :: INVALID_LON=360.
        ! longitudes above this limit are declared invalid
  
  integer :: i,j
  integer :: nrow,ncol ! grid size (# of tiles)
  integer :: iwrap     ! actual col a window tile maps into (1-based idx)
  integer :: icol_pos180,icol_neg180 ! integer col (1-based idx) of dateline
  integer :: pos180valid,neg180valid ! valid(1) / invalid(1) flags
  integer :: wrapOffset ! number tiles for hemisphere wrap around
  
  ! off hemisphere, so remap tile on opposite hemisphere
  
  ! get icol_neg180:
  call ll2ij_sinusoid(lat=latTile,lon=-180.,col=col,row=row,&
       nx=ncol,ny=nrow,&
       colGridOffset=colGridOffset,rowGridOffset=rowGridOffset,&
       mapRes=mapRes,Re=ErthRad,lat0=latorigin,lon0=lonorigin)
  icol_neg180=int(col)+1

  ! opt check row
  if (int(row)+1 /= j) then
     print *,'err[fov2sfcgrd]: process row mismatch'
     call errorHalt(1)
  end if
  ! get icol_pos180
  call ll2ij_sinusoid(lat=latTile,lon=180.,col=col,row=row,&
       nx=ncol,ny=nrow,&
       colGridOffset=colGridOffset,rowGridOffset=rowGridOffset,&
       mapRes=mapRes,Re=ErthRad,lat0=latorigin,lon0=lonorigin)
  icol_pos180=int(col)+1
  
  ! opt check row
  if (int(row)+1 /= j) then
     print *,'err[fov2sfcgrd]:  process row mismatch'
     call errorHalt(1)
  end if

  ! check if the tile containing -+180 lon is
  !    usable as is (ctr has valid lon)
  ! if not, wrap one extra
  call ij2ll_sinusoid(col=real((icol_neg180-1)+.5),row=real((j-1)+.5),&
       lat=latTile,lon=lonTile,&
       nx=ncol,ny=nrow,&
       colGridOffset=colGridOffset,rowGridOffset=rowGridOffset,&
       mapRes=mapRes,Re=ErthRad,lat0=latorigin,lon0=lonorigin)
  if (abs(lonTile) > INVALID_LON) then 
     neg180valid=0
  else
     neg180valid=1
  end if
  
  call ij2ll_sinusoid(col=real((icol_pos180-1)+.5),row=real((j-1)+.5),&
       lat=latTile,lon=lonTile,&
       nx=ncol,ny=nrow,&
       colGridOffset=colGridOffset,rowGridOffset=rowGridOffset,&
       mapRes=mapRes,Re=ErthRad,lat0=latorigin,lon0=lonorigin)
  
  if (abs(lonTile) > INVALID_LON) then 
     pos180valid=0
  else
     pos180valid=1
  end if
  ! get wrapoffset and set new tile grid col
  ! western hemisphere wraps to eastern
  if (i <= (ncol/2)) then
     wrapOffset = icol_neg180-i
!     print *,'wrapLongitude: i,ncol,ncol/2,icol_neg180,wrapOffset(initial),neg180valid,=',&
!          i,ncol,ncol/2,icol_neg180,wrapOffset,neg180valid
     if (neg180valid == 0) wrapOffset = wrapOffset+1 ! invalid, add one
!     print *,'neg180valid, wrapLongitude: wrapOffset(final)=',neg180valid,wrapOffset
!     print *,'icol_pos180,wrapOffset(final),iwrap(initial)=',icol_pos180,wrapOffset,iwrap
     iwrap = icol_pos180-wrapOffset
     if (pos180valid == 1) iwrap = iwrap+1
!     print *,'pos180valid,iwrap(final)=',pos180valid,iwrap

     ! eastern hemisphere wraps to western
  else if (i > (ncol/2)) then
     wrapOffset = i-icol_pos180
     if (pos180valid == 0) wrapOffset = wrapOffset+1
     iwrap = icol_neg180+wrapOffset
     if (neg180valid == 1) iwrap = iwrap-1
  end if
  
  !get lon of new (hemisphere-wrapped) tile 
  call ij2ll_sinusoid(col=real((iwrap-1)+.5),row=real((j-1)+.5),&
       lat=latTile,lon=lonTile,&
       nx=ncol,ny=nrow,&
       colGridOffset=colGridOffset,rowGridOffset=rowGridOffset,&
       mapRes=mapRes,Re=ErthRad,lat0=latorigin,lon0=lonorigin)
  if (abs(lonTile) > INVALID_LON) then
!     print *,'warning[fov2sfcgrd]: unable to wrap process tile to opposite hemisphere'
!     print *,'probably due to window width wrapping a distance larger than earth circumferance at current lat'
!     print *,'err[fov2sfcgrd]: unable to wrap process tile to opposite hemisphere'
!     call errorHalt(1)
  end if
  
END SUBROUTINE wrapLongitude

! =================================================================

SUBROUTINE windowBounds(col,row,colsize,rowsize,nrow,ncol,&
     colmin,colmax,rowmin,rowmax)
  
  real    :: col, row  ! float col,row (0-based idx) of grid point
  integer :: nrow,ncol ! grid size (# of tiles)
  integer :: icol,irow ! integer col,row (1-based idx) of tile
  integer :: rowmin,rowmax,colmin,colmax
  integer :: colsize ! number of col tiles beyond ctr tile
  integer :: rowsize ! number of row tiles beyond ctr tile
  
  icol=int(col)+1
  irow=int(row)+1
  
  ! set window boundary (1-based idx tiles)
  colmin=icol-colsize
  colmax=icol+colsize
  rowmin=irow-rowsize
  rowmax=irow+rowsize
  
  ! keep window within grid bounds (rows only)
  ! (note: col off grid results in wrap hemisphere)
  if (rowmin < 1) rowmin=1
  if (rowmax > nrow) THEN
     rowmax=nrow
  endif
  
end subroutine windowBounds

END MODULE FOV2sfcgrd_module
