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

!-------------------------------------------------------------------------------
!
!  MODULE GeomSpheric: contains subroutines related to geometry on
!                     a sphere, including treatment of lat/lon coordinates
!
!     Copyright AER, Inc., 2004. All rights Reserved.
!-------------------------------------------------------------------------------
MODULE GeomSpheric

! <f90Module>***********************************************************
!
! NAME:
!
!   GeomSpheric
!
! PURPOSE:
!
!   Contains subroutines related to geometry on a sphere, including treatment of 
!   lat/lon coordinates, accurate at all latitudes.
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

  USE constants, only: PI,rad2deg,deg2rad,ErthRad
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------
  !  List of Public subroutines (accessible from outside module) 
  !---------------------------------------------------------------
  PUBLIC :: gcdist, mvltln, mvltlnd, angintr, angazim, azPolReverse, &
       lonReverse,setRotation,rot2xy

  !---------------------------------------------------------------
  !     Declarations of private data (private to the module)
  !---------------------------------------------------------------
  INTEGER, parameter :: DBL=SELECTED_REAL_KIND(15)
        
CONTAINS

  FUNCTION gcdist(rlat1,rlon1,rlat2,rlon2)

!<f90Function>**********************************************************
!
! NAME:
!
!   gcdist
!
! PURPOSE:
!
!   Distance between two lat-lon points on a great circle.
!
! SYNTAX:
!
!   Results=gcdist(rlat1, rlon1, rlat2, rlon2)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   rlat1   REAL  Latitude of point 1
!   rlon1   REAL  Longitude of point 1
!   rlat2   REAL  Latitude of point 2
!   rlon2   REAL  Longitude of point 2
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

  !-- Distance between two lat-lon points on a great circle.
  !-- It makes no difference which direction is positive/negative, as long
  !-- as the two pairs of inputs are consistent.

    real, intent(in)    :: rlat1  ! (degrees)
    real, intent(in)    :: rlon1  ! (degrees)
    real, intent(in)    :: rlat2  ! (degrees)
    real, intent(in)    :: rlon2  ! (degrees)

    real  :: gcdist ! (km)

    !-- Local variables
    real :: s1,s2,c1,c2,dlon,sth,cth

    !-- Check inputs
    if((abs(rlat1).gt.90.).or.(abs(rlat2).gt.90.).or. &
      (abs(rlon1).gt.360.).or.(abs(rlon2).gt.360.))then
      print *,'err[GeomSpheric::gcdist]:'
      print *,'Lat/lon out of range:',rlat1,rlon1,rlat2,rlon2
      call errorHalt(1)
    endif
    s1=sin(rlat1*deg2rad)
    s2=sin(rlat2*deg2rad)
    c1=cos(rlat1*deg2rad)
    c2=cos(rlat2*deg2rad)
    dlon=(rlon1-rlon2)*deg2rad
    sth=sin(dlon)
    cth=cos(dlon)
    gcdist=angintr(s1,c1,s2,c2,sth,cth)*ErthRad      
    return

  END FUNCTION gcdist

  !-----------------------------------------------------------------------------

  SUBROUTINE mvltlnd(rlat1,rlon1,dismv,dirmv,rlat2,rlon2)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   mvltlnd
!
! PURPOSE:
!
!   Latitude/longitude of the point arrived at by starting at rlat1, rlon1 and 
!   moving a distance dismv in the direction dirmov. This is an alternate 
!   interface for mvltln.
!
! SYNTAX:
!
!   CALL mvltlnd(rlat1, rlon1, dismv, dirmv, rlat2, rlon2)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   rlat1  REAL  Latitude of point 1
!   rlon1  REAL  Longitude of point 1
!   dismv  REAL  distance to move
!   dirmv  REAL  direction ro move, counted from South, degrees
!   
!   INPUTS/OUTPUTS:
!   
!   rlat2  REAL  Latitude of point 2
!   rlon2  REAL  Longitude of point 2
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

  !-- Latitude/longitude of the point arrived at by starting at rlat1, rlon1
  !-- and moving a distance dismv in the direction dirmov
  !-- This is an alternate interface for mvltln.

    real, intent(in)    :: rlat1  ! (degrees, positive North)
    real, intent(in)    :: rlon1  ! (degrees, positive West)
    real, intent(in)    :: dismv  ! (km)
    real, intent(in)    :: dirmv  ! (degrees, 0=South, positive counterclockwise)
    real, intent(inout) :: rlat2  ! (degrees, positive North)
    real, intent(inout) :: rlon2  ! (degrees, positive West)

    !-- Local variables
    real :: rlat,rlon,psimv,dirmov,xlat,xlon
    
    !-- Convert from degrees to radians
    rlat=rlat1*deg2rad
    rlon=rlon1*deg2rad
    dirmov=dirmv*deg2rad

    !-- Convert from distance to angular distance (radians)
    psimv=dismv/ErthRad

    call mvltln(rlat,rlon,psimv,dirmov,xlat,xlon)

    !-- Convert from radians to degrees
    rlat2=xlat*rad2deg
    rlon2=xlon*rad2deg

  END SUBROUTINE mvltlnd

  !-----------------------------------------------------------------------------

  SUBROUTINE mvltln(rlat,rlon,psimv,dirmov,xlat,xlon)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   mvltln
!
! PURPOSE:
!
!   Latitude/longitude of the point arrived at by starting at rlat, rlon and 
!   moving an angular distance psimv in the direction dirmov
!
! SYNTAX:
!
!   CALL mvltln(rlat, rlon, psimv, dirmov, xlat, xlon)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   rlat    REAL  Latitude
!   rlon    REAL  Longitude
!   psimv   REAL  angular distance to move, degrees
!   dirmov  REAL  direction ro move, counted from South, degrees
!   
!   INPUTS/OUTPUTS:
!   
!   xlat    REAL  Latitude
!   xlon    REAL  Longitude
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

  !-- Latitude/longitude of the point arrived at by starting at RLAT, RLON
  !-- and moving an angular distance PSIMV in the direction DIRMOV

    real, intent(in)    :: rlat   ! (radians, positive North)
    real, intent(in)    :: rlon   ! (radians, positive West)
    real, intent(in)    :: psimv  ! (radians)
    real, intent(in)    :: dirmov ! (radians, 0=South, positive counterclockwise)
    real, intent(inout) :: xlat   ! (radians, positive North)
    real, intent(inout) :: xlon   ! (radians, positive West)

    real, parameter :: pio2=PI/2.

    !-- Local variables
    real :: sinr1,cosr1,sinr2,cosr2,daz,sth,cth,pol,relaz

    sinr1=sin(rlat)
    cosr1=cos(rlat)
    !-- Note that sine of equatorial angle = cosine of polar angle, etc.
    sinr2=cos(psimv)
    cosr2=sin(psimv)
    !-- Note the North pole is at az=pi; (az1-az2)=pi-(2pi-aznod)=aznod-pi
    daz=dirmov-pi
    sth=sin(daz)
    cth=cos(daz)
    pol=angintr(sinr1,cosr1,sinr2,cosr2,sth,cth)
    xlat=pio2-pol
    relaz=angazim(sinr1,cosr1,sinr2,cosr2,sth,cth)
    xlon=rlon+pi-relaz
    return

  END SUBROUTINE mvltln

  !-----------------------------------------------------------------------------

  FUNCTION angintr(s1,c1,s2,c2,sth,cth)

!<f90Function>**********************************************************
!
! NAME:
!
!   angintr
!
! PURPOSE:
!
!   Compute the interior angle of two points given in spherical coordinates.
!
! SYNTAX:
!
!   Results=angintr(s1, c1, s2, c2, sth, cth)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   s1       REAL  Sine of latitude 1
!   c1       REAL  Cosine of latitude 1
!   s2       REAL  Sine of latitude 2
!   c2       REAL  Cosine of latitude 2
!   sth      REAL  Sine of longitudes difference
!   cth      REAL  Cosine of longitudes difference
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

  !-- Compute the interior angle of two points given in spherical coordinates.
  !-- The derivation is in Alan Lipton Research Log 3/31/93.
    !-- ANGINTR is given in radians.
    !-- Input parameters are:
    !--   S1 = SIN(LAT1)
    !--   S2 = SIN(LAT2)
    !--   C1 = COS(LAT1)
    !--   C2 = COS(LAT2)
    !--   STH = SIN(LON1-LON2)
    !--   CTH = COS(LON1-LON2)
    !--   longitudes are positive west

    real, intent(in)    :: s1
    real, intent(in)    :: c1
    real, intent(in)    :: s2
    real, intent(in)    :: c2
    real, intent(in)    :: sth
    real, intent(in)    :: cth

    real  :: angintr

    !-- Local variables
    real :: zp,rp

    zp=(cth*c2*c1+s2*s1)
    rp=sqrt((cth*c2*s1-s2*c1)**2+sth**2*c2**2)
    if(abs(zp).le.abs(rp))then
      angintr=acos(zp)
    else
      if(zp.ge.0.)then
        angintr=asin(rp)
      else
        angintr=pi-asin(rp)
      endif
    endif
    return

  END FUNCTION angintr

  !-----------------------------------------------------------------------------

  FUNCTION angazim(s1,c1,s2,c2,sth,cth)

!<f90Function>**********************************************************
!
! NAME:
!
!   angazim
!
! PURPOSE:
!
!   Compute the azimuth angle from one point to another point, in spherical 
!   coordinates.
!
! SYNTAX:
!
!   Results=angazim(s1, c1, s2, c2, sth, cth)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   s1       REAL  Sine of latitude 1
!   c1       REAL  Cosine of latitude 1
!   s2       REAL  Sine of latitude 2
!   c2       REAL  Cosine of latitude 2
!   sth      REAL  Sine of longitudes difference
!   cth      REAL  Cosine of longitudes difference
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

    !-- Compute the azimuth angle from one point to another point, in 
    !-- spherical coordinates.  Point 1 is the reference point.
    !-- The derivation is in Alan Lipton Research Log 3/31/93.
    !-- ANGAZIM is given in radians, south=0, increasing counterclockwise.
    !-- ANGAZIM is between 0. and +2*PI, except
    !-- if the points are coincident, ANGAZIM = 4*PI.
    !-- Input parameters are:
    !--   S1 = SIN(LAT1)
    !--   S2 = SIN(LAT2)
    !--   C1 = COS(LAT1)
    !--   C2 = COS(LAT2)
    !--   STH = SIN(LON1-LON2)
    !--   CTH = COS(LON1-LON2)
    !--   longitudes are positive west

    real, intent(in)    :: s1
    real, intent(in)    :: c1
    real, intent(in)    :: s2
    real, intent(in)    :: c2
    real, intent(in)    :: sth
    real, intent(in)    :: cth

    real  :: angazim

    real, parameter :: twopi=2.*PI, piovr2=PI/2.

    !-- Local variables
    real :: xp,yp

    xp=cth*c2*s1-s2*c1
    yp=sth*c2
    if(xp.eq.0..and.yp.eq.0.)then
      angazim=twopi*2.
    elseif(abs(xp).ge.abs(yp))then
      angazim=atan(yp/xp)
      if(xp.lt.0.)angazim=angazim+pi
    else
      angazim=piovr2-atan(xp/yp)
      if(yp.lt.0.)angazim=angazim-pi
    endif
    if(angazim.lt.0.)angazim=angazim+twopi
    return

  END FUNCTION angazim

  !-----------------------------------------------------------------------------

  FUNCTION azPolReverse(azIn)

!<f90Function>**********************************************************
!
! NAME:
!
!   azPolReverse
!
! PURPOSE:
!
!   Convert azimuth angle from South=0, positive counterclockwise to North=0, 
!   positive clockwise or vice versa
!
! SYNTAX:
!
!   Results=azPolReverse(azIn)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   azIn          REAL  Azimuthal angle at input
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


  !  Convert azimuth angle from South=0, positive counterclockwise to
  !  North=0, positive clockwise or vice versa
  !  Returned value is positive, as long as azSCCW < 3*PI
  !  Returned value is <= 2*PI, as long as azSCCW > -3*PI

    real, intent(in)    :: azIn      ! (radians)

    real :: azPolReverse ! (radians)

    real, parameter :: twopi=2.*PI

    azPolReverse=PI-azIn

    if (azPolReverse > twopi) azPolReverse=azPolReverse-twopi
    if (azPolReverse < 0.) azPolReverse=azPolReverse+twopi

  END FUNCTION azPolReverse
  
  FUNCTION lonReverse(lonIn,lonH)

!<f90Function>**********************************************************
!
! NAME:
!
!   lonReverse
!
! PURPOSE:
!
!   Convert longitude convention
!
! SYNTAX:
!
!   Results=lonReverse(lonIn, lonH)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   lonIn       REAL  Longitude at input
!   lonH        REAL  Longitude at which to cut the sphere
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

  
    !--Convert longitude convention
    !--LonH =   0: from [-180,180] to [0,360W]
    !--LonH = 180: from [0,360W] to [-180,180]
    
    real, intent(in)  :: lonIn,lonH     ! (deg)
    real              :: lonReverse     ! (deg)
    
    if(lonIn < lonH)then
       lonReverse=-lonIn
    else
       lonReverse=360.-lonIn
    endif
  
  END FUNCTION lonReverse

  subroutine setRotation(lon1,lon2,lat1,lat2,rotation,rotationInv)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   setRotation
!
! PURPOSE:
!
!   Construct rotation matrices.
!
! SYNTAX:
!
!   CALL setRotation(lon1, lon2, lat1, lat2, rotation, rotationInv)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   lon1         REAL  Longitude of point 1
!   lon2         REAL  Longitude of point 2
!   lat1         REAL  Latitude of point 1
!   lat2         REAL  Latitude of point 2
!   
!   INPUTS/OUTPUTS:
!   
!   rotation     REAL  Rotation matrix
!   rotationInv  REAL  Inverse rotation matrix
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    implicit none
    real,                 intent(in)    :: lon1,lon2,lat1,lat2
    real, dimension(3,3), intent(inout) :: rotation,rotationInv
    real                  :: lon1t,lon2t,s1,c1,s2,c2,sth,cth,rotang
    real, dimension(3,3)  :: rot1,rot2,rot3
    lon1t=lonReverse(lon1,0.)
    lon2t=lonReverse(lon2,0.)
    s1=sin(lat1*deg2rad)
    c1=cos(lat1*deg2rad)
    s2=sin(lat2*deg2rad)
    c2=cos(lat2*deg2rad)
    sth=sin((lon1t-lon2t)*deg2rad)
    cth=cos((lon1t-lon2t)*deg2rad)
    rotang=ANGAZIM(S1,C1,S2,C2,STH,CTH)*rad2deg
    
    !  Construct the three rotation matrices.
    !  Rot1 is about the x axis
    !  Rot2 is about the y axis
    !  Rot3 is about the z axis
    
    rot1(1,1:3) = (/ cos(lon1*deg2rad),   sin(lon1*deg2rad),                   0./)
    rot1(2,1:3) = (/-sin(lon1*deg2rad),   cos(lon1*deg2rad),                   0./)
    rot1(3,1:3) = (/                0.,                  0.,                   1./)
    
    rot2(1,1:3) = (/                c1,                  0.,                   s1/)
    rot2(2,1:3) = (/                0.,                  1.,                   0./)
    rot2(3,1:3) = (/               -s1,                  0.,                   c1/)
    
    rot3(1,1:3) = (/                1.,                  0.,                   0./)
    rot3(2,1:3) = (/                0., sin(rotang*deg2rad), -cos(rotang*deg2rad)/)
    rot3(3,1:3) = (/                0., cos(rotang*deg2rad),  sin(rotang*deg2rad)/)
    
    rotation=matmul(rot3,matmul(rot2,rot1))
    rotationInv=transpose(matmul(rot3,matmul(rot2,rot1)))
    return
  end subroutine setRotation
  
  subroutine rot2xy(lat,lon,npts,rotation)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   rot2xy
!
! PURPOSE:
!
!   Completes rotation of lat/lon coordinates to a new geographic origin.
!
! SYNTAX:
!
!   CALL rot2xy(lat, lon, npts, rotation)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   npts      INTEGER  number of points
!   
!   INPUTS/OUTPUTS:
!   
!   lat       REAL     Latitude
!   lon       REAL     Longitude
!   rotation  REAL     Rotation matrix
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    ! Completes rotation of lat/lon coordinates to a new geographic origin.
    ! All angles are in degrees but are converted to radians internally.
    ! Radians on input and output would be faster.
    
    ! lat      = original latitude array 
    ! lon      = original longitude array
    ! npts     = # of points
    ! rotation = rotation matrix defined by setRotation
    
    implicit none
    integer,               intent(in)    :: npts
    real, dimension(npts), intent(inout) :: lat,lon
    real, dimension(3,3),  intent(inout) :: rotation
    real, dimension(npts) :: x,y,z
    integer               :: i
    real                  :: epsilon
    
    ! Convert angles to radians:
    lat = lat*deg2rad
    lon = lon*deg2rad
    
    !  Compute the new x,y,z point in cartesian space
    
    x = rotation(1,1) * cos(lat)*cos(lon) + &
         rotation(1,2) * cos(lat)*sin(lon) + &
         rotation(1,3) * sin(lat)
    
    y = rotation(2,1) * cos(lat)*cos(lon) + &
         rotation(2,2) * cos(lat)*sin(lon) + &
         rotation(2,3) * sin(lat)
    
    z = rotation(3,1) * cos(lat)*cos(lon) + &
         rotation(3,2) * cos(lat)*sin(lon) + &
         rotation(3,3) * sin(lat)
    
    !  Points with essentially zero x and y will be treated as 0.
    
    epsilon = 1.0E-8
    do i=1,npts
       if (abs(x(i)) <= epsilon .and. abs(y(i)) <= epsilon)then
          x(i) = 0.0
          y(i) = 0.0
       end if
    end do
  
    !  Transform the cartesian point to spherical coordinates
  
    lat = atan2(z,sqrt(x**2+y**2))*rad2deg 
    lon = atan2(y,x)*rad2deg
    return
  end subroutine rot2xy

END MODULE GeomSpheric
