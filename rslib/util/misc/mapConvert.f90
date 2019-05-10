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

module mapConvert

! <f90Module>***********************************************************
!
! NAME:
!
!   mapConvert
!
! PURPOSE:
!
!   Coordinate conversion utilities for sinusoidal grid and lat/lon
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>


  ! Coordinate conversion utilities for sinusoidal grid and lat/lon
  !
  !  Var_Name       Type       Description
  !  --------       ----       -----------
  !  col            real       grid column (grid location in x-direction)
  !  row            real       grid row (grid location in y-direction)
  !  lat            real       latitude (degrees, -90 to 90)
  !  lon            real       longitude (degrees, -180 to 180)
  !  nx             integer    number of global grid points in x-direction (i.e. West-East)
  !                            MUST BE AN EVEN NUMBER
  !  ny             integer    number of global grid points in y-direction (i.e. North-South)
  !                            MUST BE AN EVEN NUMBER
  !  mapRes         real       resolution of grid (km)
  !  Re             real       Earth radius (km)
  !  colGridOffset  real       column grid offset I.e. the offset
  !                            in columns from 0,0 in the relative tile
  !                            coordinate system to the origin of the grid 
  !                            specified by lon0
  !  rowGridOffset  real       row grid offset I.e. the offset
  !                            in rows from 0,0 in the relative tile
  !                            coordinate system to the origin of the grid 
  !                            specified by lat0
  !  lat0           real       latitude at origin of map projection  (degrees)
  !  lon0           real       longitude at origin of map projection (degrees) 
  !  colPointOffset real       column point offset I.e. the offset
  !                            in columns from the current lat/lon point
  !                            to the origin of the grid specified by lon0
  !  rowPointOffset real       row point offset I.e. the offset
  !                            in columns from the current lat/lon point
  !                            to the origin of the grid specified by lat0
  !
  ! colGridOffset and rowGridOffset determine the 
  ! absolute location of the grid/tile on the globe

  ! NOTE: Important to note that there are 2 ways of specifying the sine grid
  ! parameters: 1) mapRes(km), or 2) nx and ny (grid dimensions). If both are
  ! specified, then mapRes is used, and it is possible that only a portion of the
  ! globe will be covered by the nx by ny grid (e.g. for a single tile).
  ! Unless one is specifying a grid with an EXACT integer number of grid points
  ! between -180 deg lon and +180 deg lon at the equator, and between +90 and -90 lat
  ! then mapRes should be used, not nx and ny. 
  ! USING NX AND NY WILL ALWAYS FORCE AN EXACT INTEGER NUMBER OF GRIDS ONTO THE GLOBE
  ! AND ADJUST mapRes ACCORDINGLY.
  ! SPECIFYING mapRes WILL ADJUST NX AND NY ACCORDINGLY, ALLOWING FOR "PARTIAL GRIDS",
  ! i.e. THE STARTING AND ENDING ROWS AND COLUMNS MAY BE PARTIALLY OFF-HEMISPHERE
  !      AT THE EQUATOR.
  !
  ! LIMITATION: specified grid dimensions nx and ny must be even numbers to allow symmetry
  ! This is true even for smaller (i.e. tile) grids.

  use constants, only: ErthRad,pi,deg2rad,rad2deg,MISSING_REAL

  implicit none
  private
  public ij2ll_sinusoid, ll2ij_sinusoid, jlon2i_sinusoid, iWrap_sinusoid

  real, parameter :: llBoundTol=0.001 ! Tolerance for lat/lon out of bounds (deg)
  real, parameter :: ijBoundTol=1.e-5 ! Minimum tolerance for col/row out of bounds
  real, parameter :: ijRoundoffFactor=2. ! Maximum tolerance for col/row out of
                                         ! bounds is roundoff error times factor
  real :: phi, phi0, lambda, lambda0, x, y, rowabs, colabs

CONTAINS

  subroutine ij2ll_sinusoid(col, row, lat, lon, nx, ny, &
       colGridOffset, rowGridOffset, mapRes, Re, lat0, lon0)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   ij2ll_sinusoid
!
! PURPOSE:
!
!   Convert input column and row to lat/lon on sinusoidal grid. If result is off 
!   the Earth, returned values are MISSING_REAL
!
! SYNTAX:
!
!   CALL ij2ll_sinusoid(col, row, lat, lon, nx, ny, colGridOffset, 
!      rowGridOffset, mapRes, Re, lat0, lon0)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   col            REAL     grid column
!   row            REAL     grid row
!   nx*            INTEGER  number of global grid points in x-direction
!   ny*            INTEGER  number of global grid points in y-direction
!   colGridOffset* REAL     column grid offset
!   rowGridOffset* REAL     row grid offset
!   mapRes*        REAL     resolution of grid, for converting 
!                           spherical/sinusoidal grid 
!   Re*            REAL     Earth radius
!   lat0*          REAL     latitude at origin of map projection
!   lon0*          REAL     longitude at origin of map projection
!   
!   INPUTS/OUTPUTS:
!   
!   lat            REAL     Latitude
!   lon            REAL     Longitude
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>


    ! Convert input column and row to lat/lon on sinusoidal grid
    ! If result is off the Earth, returned values are MISSING_REAL

    real,    intent(in)           :: col, row
    integer, intent(in), optional :: nx, ny
    real,    intent(in), optional :: colGridOffset, rowGridOffset
    real,    intent(in), optional :: mapRes, Re, lat0, lon0
    real,    intent(inout)        :: lat, lon

    integer :: nx_loc, ny_loc
    real :: colGridOffset_loc, rowGridOffset_loc, mapRes_loc, re_loc, lat0_loc, lon0_loc
    real :: Recos,dlam

    call defineGrid(Re,mapRes,nx,ny,colGridOffset,rowGridOffset,lat0,lon0, &
    re_loc,mapRes_loc,nx_loc,ny_loc,colGridOffset_loc,rowGridOffset_loc, &
    lat0_loc,lon0_loc)

    ! convert origin lat/lon from degrees to radians

    phi0=lat0_loc*deg2rad
    lambda0=lon0_loc*deg2rad

    ! calculate latitude first
    
    y=(row-rowGridOffset_loc)*mapRes_loc
    phi=phi0-y/Re_loc
    lat=phi*rad2deg
    if (abs(lat) .gt. (90.+llBoundTol)) then
      lat=MISSING_REAL
    elseif (lat .gt. 90.) then
      lat=90.
    elseif (lat .lt. -90.) then
      lat=-90.
    endif

    ! calculate longitude, checking for off-hemisphere points

    x=(col-colGridOffset_loc)*mapRes_loc
    Recos=Re_loc*cos(phi)
    if (lat .eq. MISSING_REAL) then
       lon=MISSING_REAL
    elseif (abs(lat) .gt. (90.-llBoundTol)) then   ! longitude undefined at poles
       lon=0.
    else
       dlam=x/Recos
       if (abs(dlam) .gt. (pi+llBoundTol*deg2rad)) then
          lon=MISSING_REAL
       else
          if (dlam .gt. pi) then
             dlam=pi
          elseif (dlam .lt. -pi) then
             dlam=-pi
          endif
          lambda=lambda0+dlam
          lon=lambda*rad2deg
          ! 0 to 360 -> -180 to 180
          if(lon .gt. 180.)lon=lon-360.
       endif
    endif

    return
       
  end subroutine ij2ll_sinusoid

  subroutine ll2ij_sinusoid(lat, lon, col, row, nx, ny, colGridOffset, rowGridOffset, mapRes, Re, lat0, lon0, &
       colPointOffset, rowPointOffset)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   ll2ij_sinusoid
!
! PURPOSE:
!
!   Convert input lat/lon to column and row on sinusoidal grid.
!
! SYNTAX:
!
!   CALL ll2ij_sinusoid(lat, lon, col, row, nx, ny, colGridOffset, 
!      rowGridOffset, mapRes, Re, lat0, lon0, colPointOffset, 
!      rowPointOffset)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   lat             REAL     Latitude
!   lon             REAL     Longitude
!   nx*             INTEGER  number of global grid points in x-direction
!   ny*             INTEGER  number of global grid points in y-direction
!   colGridOffset*  REAL     column grid offset
!   rowGridOffset*  REAL     row grid offset
!   mapRes*         REAL     resolution of grid, for converting 
!                            spherical/sinusoidal grid 
!   Re*             REAL     Earth radius
!   lat0*           REAL     latitude at origin of map projection
!   lon0*           REAL     longitude at origin of map projection
!   
!   INPUTS/OUTPUTS:
!   
!   col             REAL     grid column
!   row             REAL     grid row
!   colPointOffset* REAL     column point offset
!   rowPointOffset* REAL     row point offset
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>


    ! Convert input lat/lon to column and row on sinusoidal grid.

    real,    intent(in)              :: lat, lon
    integer, intent(in),    optional :: nx, ny
    real,    intent(in),    optional :: colGridOffset, rowGridOffset
    real,    intent(in),    optional :: mapRes, Re, lat0, lon0
    real,    intent(inout)           :: col, row
    real,    intent(inout), optional :: colPointOffset, rowPointOffset

    integer :: nx_loc, ny_loc
    real :: colGridOffset_loc, rowGridOffset_loc, mapRes_loc, re_loc, lat0_loc, lon0_loc
    real :: lat_loc,lon180,rnxl,rnyl
    real :: cOffSmall,cBoundTol,rOffSmall,rBoundTol

    call defineGrid(Re,mapRes,nx,ny,colGridOffset,rowGridOffset,lat0,lon0, &
    re_loc,mapRes_loc,nx_loc,ny_loc,colGridOffset_loc,rowGridOffset_loc, &
    lat0_loc,lon0_loc)

    ! 0 to 360 -> -180 to 180
    lon180=lon
    if(lon .gt. 180.)lon180=lon180-360.

    ! check range of inputs

    lat_loc=lat
    if (abs(lat_loc) > (90.+llBoundTol)) then
       print *,'err[mapConvert::ll2ij_sinusoid]: ',&
            'Latitude out of bounds [-90:90]:',lat_loc
       call errorHalt(1)
    elseif (lat_loc > 90.) then
       lat_loc=90.
    elseif (lat_loc < -90.) then
       lat_loc=-90.
    endif
    if (abs(lon180) > (180.+llBoundTol)) then
       print *,'err[mapConvert::ll2ij_sinusoid]: ',&
            'Longitude out of bounds [-180:360]:',lon
       call errorHalt(1)
    elseif (lon180 > 180.) then
       lon180=180.
    elseif (lon180 < -180.) then
       lon180=-180.
    endif
    
    ! convert lat/lon from degrees to radians
    
    phi=lat_loc*deg2rad
    lambda=lon180*deg2rad
    phi0=lat0_loc*deg2rad
    lambda0=lon0_loc*deg2rad

    ! calculate physical distance x,y and then convert to col and row

    x=Re_loc*(lambda-lambda0)*cos(phi)
    y=-Re_loc*(phi-phi0)
    col=x/mapRes_loc + colGridOffset_loc
    row=y/mapRes_loc + rowGridOffset_loc

    ! prevent round-off error from giving result out of grid dimension

    rnxl=real(nx_loc)
    rnyl=real(ny_loc)
    ! estimate roundoff error in above equations for col and row
    cOffSmall=colGridOffset_loc*EPSILON(colGridOffset_loc)*ijRoundoffFactor
    cBoundTol=max(ijBoundTol,cOffSmall)
    rOffSmall=rowGridOffset_loc*EPSILON(rowGridOffset_loc)*ijRoundoffFactor
    rBoundTol=max(ijBoundTol,rOffSmall)
    if (col < -cBoundTol) then
       print *,'err[mapConvert::ll2ij_sinusoid]: ',&
            'Bound check failed to prevent negative col:',col
       call errorHalt(1)
    elseif (col < 0.) then
       col=0.
    elseif (col > rnxl+cBoundTol) then
       print *,'err[mapConvert::ll2ij_sinusoid]: ',&
            'Bound check failed to prevent out of range col:',col,nx
       call errorHalt(1)
    elseif (col >= rnxl) then
       col=rnxl-ijBoundTol
    endif
    if (row < -rBoundTol) then
       print *,'err[mapConvert::ll2ij_sinusoid]: ',&
            'Bound check failed to prevent negative row:',row
       call errorHalt(1)
    elseif (row < 0.) then
       row=0.
    elseif (row > rnyl+rBoundTol) then
       print *,'err[mapConvert::ll2ij_sinusoid]: ',&
            'Bound check failed to prevent out of range row:',row,ny
       call errorHalt(1)
    elseif (row >= rnyl) then
       row=rnyl-ijBoundTol
    endif

    ! special case at N and S Pole

    if(lat_loc .eq. 90. .or. lat_loc .eq. -90.)col=colGridOffset_loc

    ! Optionally compute colPointOffset and rowPointOffset
    ! These are the offset in columns and rows from the col and row of
    ! the input lat/lon to the origin of the grid specified by lon0 and lat0
    
    if(present(colPointOffset) .and. present(rowPointOffset))then
       colPointOffset=colGridOffset_loc-col
       rowPointOffset=rowGridOffset_loc-row
    endif

    return

  end subroutine ll2ij_sinusoid

  subroutine jlon2i_sinusoid(row, lon, col, nx, ny, &
       colGridOffset, rowGridOffset, mapRes, Re, lat0, lon0)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   jlon2i_sinusoid
!
! PURPOSE:
!
!   Given row (j) and longitiude on sinsusoidal grid, convert to corresponding 
!   column (i) in the row.
!
! SYNTAX:
!
!   CALL jlon2i_sinusoid(row, lon, col, nx, ny, colGridOffset, 
!      rowGridOffset, mapRes, Re, lat0, lon0)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   row            REAL     grid row
!   lon            REAL     Longitude
!   nx*            INTEGER  number of global grid points in x-direction
!   ny*            INTEGER  number of global grid points in y-direction
!   colGridOffset* REAL     column grid offset
!   rowGridOffset* REAL     row grid offset
!   mapRes*        REAL     resolution of grid, for converting 
!                           spherical/sinusoidal grid 
!   Re*            REAL     Earth radius
!   lat0*          REAL     latitude at origin of map projection
!   lon0*          REAL     longitude at origin of map projection
!   
!   INPUTS/OUTPUTS:
!   
!   col            REAL     grid column
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>


    ! Given row (j) and longitiude on sinsusoidal grid,
    ! convert to corresponding column (i) in the row.

    real,    intent(in)             :: row, lon
    integer, intent(in),   optional :: nx, ny
    real,    intent(in),   optional :: colGridOffset, rowGridOffset
    real,    intent(in),   optional :: mapRes, Re, lat0, lon0
    real,    intent(inout)          :: col

    integer :: nx_loc, ny_loc
    real :: colGridOffset_loc, rowGridOffset_loc, mapRes_loc, re_loc, lat0_loc, lon0_loc
    real :: lon180, lat_loc, rnxl
    real :: cOffSmall,cBoundTol

    call defineGrid(Re,mapRes,nx,ny,colGridOffset,rowGridOffset,lat0,lon0, &
    re_loc,mapRes_loc,nx_loc,ny_loc,colGridOffset_loc,rowGridOffset_loc, &
    lat0_loc,lon0_loc)

    ! convert origin lat/lon from degrees to radians

    phi0=lat0_loc*deg2rad
    lambda0=lon0_loc*deg2rad

    ! calculate latitude first
    
    y=(row-rowGridOffset_loc)*mapRes_loc
    phi=phi0-y/Re_loc
    lat_loc=phi*rad2deg

    ! Now that we have latitude of the input row,
    ! We can can convert the lat and lon to column within the row

    ! 0 to 360 -> -180 to 180
    lon180=lon
    if(lon .gt. 180.)lon180=lon180-360.

    ! check range of inputs

    if (abs(lat_loc) > (90.+llBoundTol)) then
       col=MISSING_REAL
       return
    elseif (lat_loc > 90.) then
       lat_loc=90.
    elseif (lat_loc < -90.) then
       lat_loc=-90.
    endif
    if (abs(lon180) > (180.+llBoundTol)) then
       print *,'err[mapConvert::jlon2i_sinusoid]: ',&
            'Longitude out of bounds [-180:360]:',lon
       call errorHalt(1)
    elseif (lon180 > 180.) then
       lon180=180.
    elseif (lon180 < -180.) then
       lon180=-180.
    endif      

    lambda=lon180*deg2rad

    ! calculate physical distance x and then convert to col

    x=Re_loc*(lambda-lambda0)*cos(phi)
    col=x/mapRes_loc + colGridOffset_loc

    ! prevent round-off error from giving result out of grid dimension

    rnxl=real(nx_loc)
    cOffSmall=colGridOffset_loc*EPSILON(colGridOffset_loc)*ijRoundoffFactor
    cBoundTol=max(ijBoundTol,cOffSmall)
    if (col < -cBoundTol) then
       print *,'err[mapConvert::jlon2i_sinusoid]: ',&
            'Bound check failed to prevent negative col:',col
       call errorHalt(1)
    elseif (col < 0.) then
       col=0.
    elseif (col > rnxl+cBoundTol) then
       print *,'err[mapConvert::jlon2i_sinusoid]: ',&
            'Bound check failed to prevent out of range col:',col,nx
       call errorHalt(1)
    elseif (col >= rnxl) then
       col=rnxl-ijBoundTol
    endif

    ! special case at N and S Pole

    if(lat_loc .eq. 90. .or. lat_loc .eq. -90.)col=colGridOffset_loc

    return
       
  end subroutine jlon2i_sinusoid

  function iWrap_sinusoid(row, col, nx, ny, &
       colGridOffset, rowGridOffset, mapRes, Re, lat0, lon0)

!<f90Function>**********************************************************
!
! NAME:
!
!   iWrap_sinusoid
!
! PURPOSE:
!
!   Given row (j) and column (i) on sinsusoidal grid, check whether i is out of 
!   range and, if so, find the equivalent i in range, taking account of the 
!   periodicity of the i dimension of a sinusoidal grid.
!
! SYNTAX:
!
!   Results=iWrap_sinusoid(row, col, nx, ny, colGridOffset, 
!      rowGridOffset, mapRes, Re, lat0, lon0)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   row             REAL     grid row
!   col             REAL     grid column
!   nx*             INTEGER  number of global grid points in x-direction
!   ny*             INTEGER  number of global grid points in y-direction
!   colGridOffset*  REAL     column grid offset
!   rowGridOffset*  REAL     row grid offset
!   mapRes*         REAL     resolution of grid, for converting 
!                            spherical/sinusoidal grid 
!   Re*             REAL     Earth radius
!   lat0*           REAL     latitude at origin of map projection
!   lon0*           REAL     longitude at origin of map projection
!
!   * OPTIONAL
!
! RETURN:
!
!     REAL     
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>


    ! Given row (j) and column (i) on sinsusoidal grid,
    ! check whether i is out of range and, if so, find the equivalent i in range,
    ! taking account of the periodicity of the i dimension of a sinusoidal grid

    real,    intent(in)           :: row, col
    integer, intent(in), optional :: nx, ny
    real,    intent(in), optional :: colGridOffset, rowGridOffset
    real,    intent(in), optional :: mapRes, Re, lat0, lon0
    real     :: iWrap_sinusoid

    integer :: nx_loc, ny_loc
    real :: colGridOffset_loc, rowGridOffset_loc, mapRes_loc, re_loc, lat0_loc, lon0_loc
    real :: lat_loc, xWrap,Recos,dlam

    call defineGrid(Re,mapRes,nx,ny,colGridOffset,rowGridOffset,lat0,lon0, &
    re_loc,mapRes_loc,nx_loc,ny_loc,colGridOffset_loc,rowGridOffset_loc, &
    lat0_loc,lon0_loc)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   defineGrid
!
! PURPOSE:
!
!   Define grid scale and size based on user inputs. Scale and size are not 
!   independent, so user has options about which parameters to provide. It is 
!   intended that all optional arguments are present in the calling statement, 
!   but may not be present in the parent of the calling subroutine.
!
! SYNTAX:
!
!   CALL defineGrid(Re, mapRes, nx, ny, colGridOffset, rowGridOffset, 
!      lat0, lon0, re_loc, mapRes_loc, nx_loc, ny_loc, 
!      colGridOffset_loc, rowGridOffset_loc, lat0_loc, lon0_loc)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   Re*                REAL     Earth radius
!   mapRes*            REAL     resolution of grid, for converting 
!                               spherical/sinusoidal grid 
!   nx*                INTEGER  number of global grid points in x-direction
!   ny*                INTEGER  number of global grid points in y-direction
!   colGridOffset*     REAL     column grid offset
!   rowGridOffset*     REAL     row grid offset
!   lat0*              REAL     latitude at origin of map projection
!   lon0*              REAL     longitude at origin of map projection
!   
!   INPUTS/OUTPUTS:
!   
!   re_loc             REAL     Earth radius
!   mapRes_loc         REAL     resolution of grid, for converting 
!                               spherical/sinusoidal grid 
!   nx_loc             INTEGER  number of global grid points in x-direction
!   ny_loc             INTEGER  number of global grid points in y-direction
!   colGridOffset_loc  REAL     column grid offset
!   rowGridOffset_loc  REAL     row grid offset
!   lat0_loc           REAL     latitude at origin of map projection
!   lon0_loc           REAL     longitude at origin of map projection
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>


    ! convert origin lat/lon from degrees to radians

    phi0=lat0_loc*deg2rad

    ! calculate latitude first
    
    y=(row-rowGridOffset_loc)*mapRes_loc
    phi=phi0-y/Re_loc
    lat_loc=phi*rad2deg

    ! check range of inputs

    if (abs(lat_loc) > (90.+llBoundTol)) then
       iWrap_sinusoid=MISSING_REAL
       return
    endif

    x=(col-colGridOffset_loc)*mapRes_loc
    Recos=Re_loc*cos(phi)
    if (abs(lat_loc) .gt. (90.-llBoundTol)) then   ! longitude undefined at poles
       iWrap_sinusoid=colGridOffset_loc
    else
       dlam=x/Recos
       if (dlam .gt. pi) then
          dlam=dlam-pi*2.
       elseif (dlam .lt. -pi) then
          dlam=dlam+pi*2.
       endif
    endif
    xWrap=dlam*Recos
    iWrap_sinusoid=xWrap/mapRes_loc+colGridOffset_loc

    return
       
  end function iWrap_sinusoid

  subroutine defineGrid(Re,mapRes,nx,ny,colGridOffset,rowGridOffset,lat0,lon0, &
    re_loc,mapRes_loc,nx_loc,ny_loc,colGridOffset_loc,rowGridOffset_loc, &
    lat0_loc,lon0_loc)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   defineGrid
!
! PURPOSE:
!
!   Define grid scale and size based on user inputs. Scale and size are not 
!   independent, so user has options about which parameters to provide. It is 
!   intended that all optional arguments are present in the calling statement, 
!   but may not be present in the parent of the calling subroutine.
!
! SYNTAX:
!
!   CALL defineGrid(Re, mapRes, nx, ny, colGridOffset, rowGridOffset, 
!      lat0, lon0, re_loc, mapRes_loc, nx_loc, ny_loc, 
!      colGridOffset_loc, rowGridOffset_loc, lat0_loc, lon0_loc)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   Re*                REAL     Earth radius
!   mapRes*            REAL     resolution of grid, for converting 
!                               spherical/sinusoidal grid 
!   nx*                INTEGER  number of global grid points in x-direction
!   ny*                INTEGER  number of global grid points in y-direction
!   colGridOffset*     REAL     column grid offset
!   rowGridOffset*     REAL     row grid offset
!   lat0*              REAL     latitude at origin of map projection
!   lon0*              REAL     longitude at origin of map projection
!   
!   INPUTS/OUTPUTS:
!   
!   re_loc             REAL     Earth radius
!   mapRes_loc         REAL     resolution of grid, for converting 
!                               spherical/sinusoidal grid 
!   nx_loc             INTEGER  number of global grid points in x-direction
!   ny_loc             INTEGER  number of global grid points in y-direction
!   colGridOffset_loc  REAL     column grid offset
!   rowGridOffset_loc  REAL     row grid offset
!   lat0_loc           REAL     latitude at origin of map projection
!   lon0_loc           REAL     longitude at origin of map projection
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>


    ! Define grid mapRes and size based on user inputs.  MapRes and size are
    ! not independent, so user has options about which parameters to provide.
    ! It is intended that all optional arguments are present in the calling statement,
    ! but may not be present in the parent of the calling subroutine.

    real,    intent(in),    optional :: Re
    real,    intent(in),    optional :: mapRes
    integer, intent(in),    optional :: nx,ny
    real,    intent(in),    optional :: colGridOffset,rowGridOffset
    real,    intent(in),    optional :: lat0,lon0
    real,    intent(inout)           :: re_loc
    real,    intent(inout)           :: mapRes_loc
    integer, intent(inout)           :: nx_loc,ny_loc
    real,    intent(inout)           :: colGridOffset_loc,rowGridOffset_loc
    real,    intent(inout)           :: lat0_loc,lon0_loc

    real :: rnx, rny

    if(present(Re))then
       re_loc=Re
    else
       re_loc=ErthRad
    endif

    ! If all three grid parameters present, use mapRes to define mapping
    if(present(nx) .and. present(ny) .and. present(mapRes))then
       nx_loc=nx
       ny_loc=ny
       mapRes_loc=mapRes
    elseif (present(mapRes)) then
    ! MapRes is present, calculate nx and ny which must be even numbers
       mapRes_loc=mapRes
       rnx=pi*re_loc/mapRes_loc
       if(mod(rnx,1.) .eq. 0.)then
          nx_loc=2*int(rnx)
       else
          nx_loc=2*(int(rnx)+1)
       endif
       rny=rnx/2.
       if(mod(rny,1.) .eq. 0.)then
          ny_loc=2*int(rny)
       else
          ny_loc=2*(int(rny)+1)
       endif
    elseif(present(nx) .and. present(ny))then
    ! MapRes is not present, calculate based on nx and ny
       nx_loc=nx
       ny_loc=ny
       mapRes_loc=2.*pi*re_loc/real(nx_loc)
    else
       print *,'err[mapConvert::defineGrid]: '
       print *,'Requires input mapRes or (nx,ny); presence of mapRes,nx,ny:', &
          present(mapRes),present(nx),present(ny) 
       call errorHalt(1)
    endif

    ! nx and ny must be even numbers
    ! so that grid is defined symmetrically
    ! such that lat0 and lon0 are mapped to col=nx/2 and row=ny/2

    if(mod(nx_loc,2) .ne. 0)then
       print *,'err[mapConvert::defineGrid]: ',&
            ' nx not an even number=',nx
       call errorHalt(1)
    endif
    if(mod(ny_loc,2) .ne. 0)then
       print *,'err[mapConvert::defineGrid]: ',&
            ' ny not an even number=',ny
       call errorHalt(1)
    endif

    ! if not present, set colGridOffset and rowGridOffset
    ! assuming grid 0,0 is at extreme top left (i.e. NW)
    ! corner of grid, defined by nx, ny, and mapRes
    ! This would generally be the case
    ! for a single tile, representing the entire globe

    if(.not. present(colGridOffset))then
       colGridOffset_loc=real(nx_loc)/2.
    else
       colGridOffset_loc=colGridOffset
    endif
    if(.not. present(rowGridOffset))then
       rowGridOffset_loc=real(ny_loc)/2.
    else
       rowGridOffset_loc=rowGridOffset
    endif

    if(.not. present(lat0))then
       lat0_loc=0.
    else
       lat0_loc=lat0
    endif
    if(.not. present(lon0))then
       lon0_loc=0.
    else
       lon0_loc=lon0
    endif

    return

  end subroutine defineGrid

end module mapConvert
