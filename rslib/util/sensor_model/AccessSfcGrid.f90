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

module AccessSfcGrid

! <f90Module>***********************************************************
!
! NAME:
!
!   AccessSfcGrid
!
! PURPOSE:
!
!   Tools for extracting data from a file of gridded surface grid.
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>


  use mapConvert
  use sfcgrd_io_module
  use date_module, only: julian
  use constants, only: rad2deg,deg2rad,ErthRad,MISSING_REAL,MISSING_INT
  use GeomSpheric, only: gcdist

  implicit none
  private
  public loadSfcGrid, setSfcLoc, destroySfcGrid, isLoadedSfcGrid

  real,    parameter :: latCartLimit=70.  ! latitude extent of Cartesian
                                          ! approximation for distances
  real,    parameter :: minValid=0.  ! Assumes any lower data are invalid
  integer, parameter :: mxpts=200    ! Max points to average
  integer :: ncid=-1  ! Negative indicates not open
  integer :: iTime, nTime, mynchan
  character (len=100) :: casename
  character (len=20)  :: proj
  character (len=5)   :: scaleunits,scaleunitsClean
  
  real    :: rowGridOffset,colGridOffset,latorigin,lonorigin,Rearth,scale
  integer :: itile, jtile, nitile, njtile
  integer, dimension(1) :: imask
  real, dimension(2) :: sto_sfc, sto_npts
  real    :: aRangeNormal
  real    :: scaleDeg
  logical :: myFrqPolGridded=.false.
  logical :: myForceAverage=.false.

  integer                                    :: ncol,nrow,nValsPerGrid
  integer,   dimension(:),       allocatable :: mapOrder
  integer*2, dimension(:),       allocatable :: year1d
  integer,   dimension(:),       allocatable :: ijulday1d
  integer*2, dimension(:),       allocatable :: datBuf_i2
  real,      dimension(:),       allocatable :: datBuf
  real,      dimension(:,:),     allocatable :: datBuf_avg

CONTAINS 

  subroutine loadSfcGrid(ncfile,nchan,frq,pol,wvn,icov,time,acceptRange, &
    forceAverage,frqPolGridded)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   loadSfcGrid
!
! PURPOSE:
!
!   Initialize access to the grid file by loading properties and creating 
!   conversion/interpolation maps.
!
! SYNTAX:
!
!   CALL loadSfcGrid(ncfile, nchan, frq, pol, wvn, icov, time, 
!      acceptRange, forceAverage, frqPolGridded)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   ncfile         CHAR     path of file
!   icov*          INTEGER  flag telling what covariance data are present
!   time*          INTEGER  Date/time
!   acceptRange*   REAL     range to accept data, when no valid data
!   forceAverage*  LOGICAL  Flag to force using average of valid surrounding 
!                           data 
!   
!   INPUTS/OUTPUTS:
!   
!   nchan          INTEGER  Number of channels
!   frq*           REAL     Channel frequencies
!   pol*           INTEGER  Polarization
!   wvn*           REAL     Wavenumbers
!   frqPolGridded* LOGICAL  Flag that data are on a regular 
!                           frequency/polarization grid
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    ! Load gridded data
    character (len=*),     intent(in)              :: ncfile
    integer,               intent(inout)           :: nchan
    real,    dimension(:), pointer,       optional :: frq
    integer, dimension(:), pointer,       optional :: pol
    real,    dimension(:), pointer,       optional :: wvn
    integer,               intent(in),    optional :: icov
    integer, dimension(6), intent(in),    optional :: time
    real,                  intent(in),    optional :: acceptRange
    logical,               intent(in),    optional :: forceAverage
    logical,               intent(inout), optional :: frqPolGridded

    ! acceptRange = range to accept data, when no valid data
    !               are at the nearest grid point, in terms of grid steps

    ! Accept data within aRangeDflt*scale when no valid data
    ! at nearest point

    ! forceAverage ==> Always use an average of valid surrounding points,
    !                  even if the nearest neighbor is valid
    real,    parameter :: aRangeDflt = 2.  
    integer :: i,myicov=0
    real    :: fPrevV,fPrevH
    integer :: nV,nH 
    real,    dimension(:), allocatable :: frqTmpV,frqTmpH
    real,    dimension(:), allocatable       :: julday1d

    if (present(icov)) myicov=icov
    if (present(forceAverage)) myForceAverage=forceAverage

    !---Free memory from any previous load
    call destroySfcGrid()

    !---Get database parameters
    call querySfcGrd(file=ncfile, &
         ncol=ncol, nrow=nrow, nValsPerGrid=nValsPerGrid, nTimeLevels=nTime, &
         casename=casename, nchan=mynchan, &
         proj=proj, latorigin=latorigin, lonorigin=lonorigin, &
         rowGridOffset=rowGridOffset, colGridOffset=colGridOffset, &
         scale=scale, scaleunits=scaleunits, re=Rearth, &
         itile=itile, jtile=jtile, nitile=nitile, njtile=njtile, &
         sto_x3d1_i2=sto_sfc)
    if(myicov == 0 .and. nValsPerGrid /= mynchan)then
       print *,'err[AccessSfcGrid::loadSfcGrid]:  error in 1D file: ',&
            'nValsPerGrid /= mynchan; ',nValsPerGrid, mynchan
       call errorHalt(1)
    elseif(myicov == 1 .and. nValsPerGrid /= mynchan*mynchan)then
       print *,'err[AccessSfcGrid::loadSfcGrid]:  error in 2D file: ',&
            'nValsPerGrid /= mynchan**2; ',nValsPerGrid, mynchan
       call errorHalt(1)
    elseif(myicov == 2 .and. nValsPerGrid /= 1)then
       print *,'err[AccessSfcGrid::loadSfcGrid]:  error in scalar file: ',&
            'nValsPerGrid /= 1; ',nValsPerGrid
       call errorHalt(1)
    endif
    if (present(acceptRange)) then
      aRangeNormal=acceptRange
    else
      aRangeNormal=aRangeDflt
    endif  
    ! A neighbor grid box center is never closer than 0.5 grid steps
    aRangeNormal=max(aRangeNormal,0.5)
    scaleunits=adjustl(scaleunits)
    ! netCDF pads with char(0)
    scaleunitsClean=scaleunits(1:verify(scaleunits,char(0),back=.true.))
    write(*,'(''setSfcLoc will accept data, as needed, within'',f7.2,1x,a5)') &
      aRangeNormal*scale,trim(scaleunitsClean)
    scaleDeg=scale*rad2deg/Rearth
    !---Now that nchan is known, open the file and get remaining parameters    
    nchan=mynchan
    if (allocated(julday1d)) deallocate(julday1d)
    if (present(frq) .and. present(pol)) then
       if (associated(frq)) deallocate(frq,pol)
       allocate(frq(mynchan),pol(mynchan))
       allocate(frqTmpV(mynchan),frqTmpH(mynchan))
       call openSfcGrd(ncid=ncid, file=ncfile,polarity=pol(1:mynchan), &
            freq=frq(1:mynchan),status='old')
    elseif (present(wvn)) then
       if (associated(wvn)) deallocate(wvn)
       allocate(wvn(mynchan))
       call openSfcGrd(ncid=ncid, file=ncfile, wvn=wvn(1:mynchan),status='old')
    endif
    allocate(mapOrder(mynchan))
    !---Allocate arrays and read in surface data
    allocate(year1d(nTime),julday1d(nTime),ijulday1d(nTime))
    call getSfcGrd(ncid=ncid, imask=imask, &
         year1d=year1d(1:nTime), &
         julday1d=julday1d(1:nTime))
    ijulday1d=nint(julday1d)
    !---Check for data came on frequency grid in standard order
    if(myicov >= 1)then
       deallocate(mapOrder)
       allocate(mapOrder(nValsPerGrid))
       mapOrder=(/(i,i=1,nValsPerGrid)/)
       myFrqPolGridded=.false.
    elseif (present(frq) .and. present(pol)) then
      if (all(frq(1:mynchan:2) == frq(2:mynchan:2)) .and. &
           all(pol(1:mynchan:2) == 0) .or. &
           all(pol(2:mynchan:2) == 1)) then 
        mapOrder=(/(i,i=1,mynchan)/)
        myFrqPolGridded=.true.
      !---Check whether data came on frequency grid in non-standard order
      else
        fPrevV=0.
        fPrevH=0.
        nV=0
        nH=0
        myFrqPolGridded=.false.         ! default
        do i=1,mynchan
          if (pol(i) == 0) then  ! V pol
            if (frq(i) <= fPrevV) exit
            nV=nV+1
            frqTmpV(nV)=frq(i)
            mapOrder(nV*2-1)=i
            fPrevV=frq(i)
          elseif (pol(i) == 1) then  ! H pol
            if (frq(i) <= fPrevH) exit
            nH=nH+1
            frqTmpH(nH)=frq(i)
            mapOrder(nH*2)=i
            fPrevH=frq(i)
          else
            exit
          endif
        enddo
        if (nV == nH) then
          if (all(frqTmpV(1:nV) == frqTmpH(1:nH))) then   ! frq/pol after mapping
            myFrqPolGridded=.true.
            frq(1:mynchan:2)=frqTmpV(1:mynchan/2)
            frq(2:mynchan:2)=frqTmpH(1:mynchan/2)
            pol(1:mynchan:2)=0
            pol(2:mynchan:2)=1
          endif
        endif
        if (.not. myFrqPolGridded) mapOrder=(/(i,i=1,mynchan)/)      ! default
      endif
    else
      mapOrder=(/(i,i=1,mynchan)/)
      myFrqPolGridded=.false.
    endif
    if (present(frqPolGridded)) frqPolGridded=myFrqPolGridded
    !--- Option to choose time level once and for all
    if (nTime == 1) then
      iTime=1
    elseif (present(time)) then
      iTime=iClosestDate(time,ijulday1d)
    else
      iTime=MISSING_INT
    endif
    deallocate(julday1d)
    if (allocated(frqTmpV)) deallocate(frqTmpV,frqTmpH)
    !--- Allocate memory for data, used in setSfcLoc
    allocate(datBuf_i2(nValsPerGrid),datBuf(nValsPerGrid), &
      datBuf_avg(mxpts,nValsPerGrid))

  end subroutine loadSfcGrid

!-----------------------------------------------------------------------------

  subroutine setSfcLoc(xlat,xlon,sfcData,qflag,time,acceptRange)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   setSfcLoc
!
! PURPOSE:
!
!   Obtain gridded data at or near a lat/lon location
!
! SYNTAX:
!
!   CALL setSfcLoc(xlat, xlon, sfcData, qflag, time, acceptRange)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   xlat         REAL     Latitude
!   xlon         REAL     Longitude
!   time*        INTEGER  Date/time
!   acceptRange* REAL     range to accept data, when no valid data
!   
!   INPUTS/OUTPUTS:
!   
!   sfcData      REAL     Surface data
!   qflag        INTEGER  Flag telling status of finding valid data
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

  ! Obtain gridded data at or near a lat/lon location
    real,                  intent(in)              :: xlat,xlon
    real,    dimension(:), intent(inout)           :: sfcData
    integer,               intent(inout)           :: qflag
    integer, dimension(6), intent(in),    optional :: time(6)
    real,                  intent(in),    optional :: acceptRange

    real,    parameter        :: lat2dist=ErthRad*deg2rad
    real,    parameter        :: exro=0.01 ! extra search area, for roundoff
    integer :: icol,irow,icol0m,icol0p,irow0m,irow0p
    integer :: i,jj,ii,npts,kk,js
    real    :: col,row,crow,dcol,wcol,ccol,dist,glat,glon,dlon
    real    :: aRange,cRange,coslat,latRange,lonRange,extram,extrap

    if (present(acceptRange)) then
      aRange=acceptRange
    else
      ! Retain range that was set in loadSfcGrid
      aRange=aRangeNormal
    endif
    ! A neighbor grid box center is never closer than 0.5 grid steps
    aRange=max(aRange,0.5)
    latRange=aRange*scale/lat2dist
    ! 0.5 excludes if grid box center not in range, exro is extra for round-off
    cRange=aRange-0.5+exro
    !----Find nearest grid point 
    call ll2ij_sinusoid(lat=xlat, lon=xlon, nx=ncol, ny=nrow, mapRes=scale, &
         colGridOffset=colGridOffset, rowGridOffset=rowGridOffset, &
         col=col, row=row)
    icol=floor(col)+1  ! Grid box is indexed by its corner
    irow=floor(row)+1
    if (nTime == 1) then
      iTime=1
    elseif (present(time)) then
      iTime=iClosestDate(time,ijulday1d)
    elseif (iTime == MISSING_INT) then
      print *,'err[AccessSfcGrid::setSfcLoc]: Need to select time'
      call errorHalt(1)
    endif
    !----By default, qflag=0 which means that nothing matches this (i,j) point
    qflag=0   
    !----Without interpolation: at least one point close by to use (qflag=1)
    call getSfcGrd(ncid=ncid, imask=imask, iloc=icol, jloc=irow, kloc=iTime, &
       x3d1_i2=datBuf_i2)
    datBuf=sto_sfc(1)+real(datBuf_i2(mapOrder))*sto_sfc(2)
    if(datBuf(1) .gt. minValid .and. .not. myForceAverage)then
       qflag=1
       sfcData=datBuf
       return
    endif
    npts=0    
    !----No point (icol,irow) found, so look for surrounding points
    if(qflag .eq. 0)then
       !----Find rows in range
       coslat=cos(xlat*deg2rad)
       if (coslat > 0.) then 
         lonRange=latRange/coslat
       else
         lonRange=360.
       endif
       ! When close to edge of grid, extra range accounts for grid not uniform 
       !   spacing across edge (roundoff when wrapping)
       extram=0.
       extrap=0.
       if (xlon+180. < lonRange) extram=1.
       if (180.-xlon < lonRange) extrap=1.
       irow0m=floor(row-cRange) 
       irow0p=floor(row+cRange)
       do js=irow0m,irow0p     ! zero-based
          crow=real(js)+0.5  ! Center of row
          jj=js+1
          if (jj < 1 .or. jj > nrow) cycle
          !---- Account for shifting alignment of points on sine grid
          call jlon2i_sinusoid(row=crow, lon=xlon, nx=ncol, ny=nrow, &
            mapRes=scale, colGridOffset=colGridOffset, &
            rowGridOffset=rowGridOffset, col=col)
       !----Find column in range
       !----In cartesian grid, column range narrows away from center row, but 
       !----here range is not narrowed because sine grid rows curve near poles.
       !----Exact distance check (below) eliminates points out of range.
          icol0m=floor(col-cRange-extram)
          icol0p=floor(col+cRange+extrap)
          do i=icol0m,icol0p     ! zero-based
             dcol=real(i)+0.5  ! center of column (real)
             !---- Account for periodicity of x in sine grid
             wcol=iWrap_sinusoid(row=crow, col=dcol, nx=ncol, ny=nrow, &
             mapRes=scale, colGridOffset=colGridOffset, &
             rowGridOffset=rowGridOffset)
             !--- After wrapping, reset to center of a grid box
             ccol=real(floor(wcol))+0.5
             !--- Check whether point is in range
             call ij2ll_sinusoid(col=ccol, row=crow, nx=ncol, ny=nrow, &
                mapRes=scale, colGridOffset=colGridOffset, &
                rowGridOffset=rowGridOffset, lat=glat, lon=glon)
             if (glat == MISSING_REAL .or. glon == MISSING_REAL) cycle
             if (abs(xlat) <= latCartLimit) then
               dlon=minval((/abs(glon-xlon),(360.+glon-xlon),(360.-glon+xlon)/))
               dist=sqrt((dlon*coslat)**2+abs(glat-xlat)**2)*lat2dist
             else
               dist=gcdist(glat,glon,xlat,xlon)
             endif
             if (dist > aRange*scale) cycle
             ii=floor(ccol)+1
             call getSfcGrd(ncid=ncid, imask=imask, &
                iloc=ii, jloc=jj, kloc=iTime, x3d1_i2=datBuf_i2)
             datBuf=sto_sfc(1)+real(datBuf_i2(mapOrder))*sto_sfc(2)
             if(datBuf(1) .gt. minValid)then
                npts=npts+1
                if (npts > mxpts) then
                  print *,'err[AccessSfcGrid::setSfcLoc]: Overran dimension;' &
                     ,' npts > mxpts:',npts,mxpts
                  call errorHalt(1)
                endif
                datBuf_avg(npts,:)=datBuf
             endif
          enddo
       enddo
       if(npts.ge.1) then !Found at least one surrounding point (qflag=2)
          qflag=2
          do kk=1,nValsPerGrid
             sfcData(kk)=sum(datBuf_avg(1:npts,kk))/float(npts)
          enddo
       endif
    endif
    return
  end subroutine setSfcLoc

!-----------------------------------------------------------------------------

  function iClosestDate(timeReq,jdays)

!<f90Function>**********************************************************
!
! NAME:
!
!   iClosestDate
!
! PURPOSE:
!
!   For a given date, find the index of the database array giving the closest 
!   temporal match.
!
! SYNTAX:
!
!   Results=iClosestDate(timeReq, jdays)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   timeReq       INTEGER  requested date/time
!   jdays         INTEGER  Julian days list
!
!   * OPTIONAL
!
! RETURN:
!
!     INTEGER  
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

     integer, dimension(:), intent(in)    :: timeReq
     integer, dimension(:), intent(in)    :: jdays
     integer  :: iClosestDate

     integer, parameter :: daysYr=365  ! days per year

     integer :: nDays
     integer :: jDiff,jdMin,jdayReq,k

     nDays=size(jdays)
     if (nDays == 1) then
       iClosestDate=1
       return
     endif

     jdayReq=julian(timeReq(3),timeReq(2),timeReq(1))
     jdMin=999
     ! Find nearest time, by season
     do k=1,nDays
       jDiff=min(abs(jdayReq-jdays(k)),(jdayReq+daysYr-jdays(k)))
       if (jDiff < jdMin) then
         jdMin=jDiff
         iClosestDate=k
       endif
     enddo
     if (jdMin >= 999) then
       print *,'err[AccessSfcGrid::iClosestDate]:  Time match failed;', &
           jdayReq,jdays(1:nDays)
       call errorHalt(1)
     endif

   end function iClosestDate

!-----------------------------------------------------------------------------

   subroutine destroySfcGrid()

!<f90Subroutine>********************************************************
!
! NAME:
!
!   destroySfcGrid
!
! PURPOSE:
!
!   Free dynamic memory of global variables and close netCDF file
!
! SYNTAX:
!
!   CALL destroySfcGrid()
!
! ARGUMENTS:
!
!   
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

   ! Free dynamic memory of global variables and close netCDF file
     
    !---Close database file
     if (ncid >= 0) then
        call closeSfcGrd(ncid=ncid)
        ncid=-1
     endif

     if (allocated(datBuf_i2)) deallocate(datBuf_i2)
     if (allocated(datBuf)) deallocate(datBuf)
     if (allocated(datBuf_avg)) deallocate(datBuf_avg)
     if (allocated(year1d)) deallocate(year1d,ijulday1d)
     if (allocated(mapOrder)) deallocate(mapOrder)
   end subroutine destroySfcGrid

!-----------------------------------------------------------------------------

   function isLoadedSfcGrid()

!<f90Function>**********************************************************
!
! NAME:
!
!   isLoadedSfcGrid
!
! PURPOSE:
!
!    Allow parent to know whether file is loaded, since can only do one at a 
!    time
!
! SYNTAX:
!
!   Results=isLoadedSfcGrid()
!
! ARGUMENTS:
!
!   
!
!   * OPTIONAL
!
! RETURN:
!
!     LOGICAL  
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

   ! Allow parent to know whether file is loaded, since can only
   ! do one at a time
   logical :: isLoadedSfcGrid

     isLoadedSfcGrid=.FALSE.
     if (ncid >= 0) isLoadedSfcGrid=.TRUE.
     return

   end function isLoadedSfcGrid

end module AccessSfcGrid
