FUNCTION GRAV(Z, WINDE, WINDN, LAT, LON)
!***********************************************************************
!* Function name: grav
!* Purpose: Get acceleration due to gravity
!* Usage: call GRAV(Z, WINDE, WINDN, LAT, LON)
!* Description: Calculate the altitude adjusted accel. due to gravity
!* Inputs:
!* Var_Name      Type       Description
!* --------      ----       -----------
!* z             real       altitude
!* winde         real       wind velocity, east
!* windn         real       wind velocity, north
!* xlat          real       latitude in degrees
!* xlon          real       longitue in degrees
!*
!* Outputs:
!* Var_Name      Type       Description
!* --------      ----       -----------
!* grav          real       aceel due to gravity
!*
!* Copyright: n.a.
!* Modified by Atmospheric and Environmental Research, Inc.
!***********************************************************************
  use constants, ONLY : pi,deg2rad
  IMPLICIT NONE
  REAL, intent(in) :: Z, WINDE, WINDN, LAT, LON
  real :: G_SUR, R, COSLT, COSLT2, SINLT2, SIN2LT, COSLON, &
       LTRAD, C_SUR, C_Z, RTOT, GRAVZ, GRAV
  real, parameter :: B2=4.041031E+13, ABTERM=6.724285E-3
! Constants for normal gravity equation
! (see "American Institute of Physics Handbook", 1963, pg 2-102)
  real, parameter :: G0=9.780455,c1=5.30157E-3,c2=-5.85E-6,c3=6.40E-6
  real, parameter :: W=1.1574074E-5 ! Earth's rotational speed
  real, parameter :: pi2 = 2.*pi

! Calculate longitude term
! Add offset of 18 degrees, double it, convert to radians, and take the cosine

  COSLON=COS( deg2rad*2.0*( LON + 18.0 ) )
!
!     Calculate the latitude terms
  LTRAD=deg2rad*LAT ! Convert Latitude to radians
!     Calculate sine and cosine terms
  COSLT = COS(LTRAD)
  COSLT2 = COSLT**2
  SINLT2 = ( SIN(LTRAD ) )**2
  SIN2LT = ( SIN( 2.0*LTRAD ) )**2

  R = SQRT( B2/( 1.0 - COSLT2*ABTERM ) ) ! Earth's radius at this latitude
  RTOT = R + Z ! total distance from Earth's center
  G_SUR = G0*( 1.0 + C1*SINLT2 + C2*SIN2LT + C3*COSLT2*COSLON ) ! gravity at the Earth's surface
  C_SUR = COSLT2*R*(PI2*W)**2 ! centripetal term at the Earth's surface
  C_Z = ( ( PI2*RTOT*COSLT*W + WINDE )**2 + (WINDN)**2 )/RTOT ! centripetal term at altitude z (with wind)
  GRAVZ=(G_SUR + C_SUR)*(1.0 - R**2/RTOT**2) ! change in gravitation with altitude

  GRAV = G_SUR + (C_SUR - C_Z) - GRAVZ

  RETURN
END FUNCTION GRAV
