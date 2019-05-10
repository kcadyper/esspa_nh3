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

! &"$Id$"

module MetFunctions

! <f90Module>***********************************************************
!
! NAME:
!
!   MetFunctions
!
! PURPOSE:
!
!   Meteorological functions.
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

  use constants, only: ratioMwt,Rd,g0
  implicit none
  private
  ! estab contains svp from -90 to 65 deg.c at 1 deg intervals
  real, dimension(156) :: estab
  data estab/ &
       9.672e-5,1.160e-4,1.388e-4,1.658e-4,1.977e-4,2.353e-4, &
       2.796e-4,3.316e-4,3.925e-4,4.638e-4,5.472e-4,6.444e-4, &
       7.577e-4,8.894e-4,1.042e-3,1.220e-3,1.425e-3,1.662e-3, &
       1.936e-3,2.252e-3,2.615e-3,3.032e-3,3.511e-3,4.060e-3, &
       4.688e-3,5.406e-3,6.225e-3,7.159e-3,8.223e-3,9.432e-3, &
       1.080e-2,1.236e-2,1.413e-2,1.612e-2,1.838e-2,2.092e-2, &
       2.380e-2,2.703e-2,3.067e-2,3.476e-2,3.935e-2,4.449e-2, &
       5.026e-2,5.671e-2,6.393e-2,7.198e-2,8.097e-2,9.098e-2, &
       1.021e-1,1.145e-1,1.283e-1,1.436e-1,1.606e-1,1.794e-1, &
       2.002e-1,2.233e-1,2.488e-1,2.769e-1,3.079e-1,3.421e-1, &
       3.798e-1,4.213e-1,4.669e-1,5.170e-1,5.720e-1,6.323e-1, &
       6.985e-1,7.709e-1,8.502e-1,9.370e-1,   1.032,   1.135, &
       1.248,   1.371,   1.506,   1.652,      1.811,   1.984, &
       2.172,   2.376,   2.597,   2.839,      3.097,   3.522, &
       3.8619,  4.2148,  4.5451,  4.8981,    5.2753,   5.678, &
       6.1078,  6.5662,  7.0547,  7.5753,    8.1294,  8.7192, &
       9.3465,  10.013,  10.722,  11.474,    12.272,  13.119, &
       14.017,  14.969,  15.977,  17.044,    18.173,  19.367, &
       20.630,  21.964,  23.373,  24.861,    26.430,  28.086, &
       29.831,  31.671,  33.608,  35.649,    37.796,  40.055, &
       42.430,  44.927,  47.551,  50.307,    53.200,  56.236, &
       59.422,  62.762,  66.264,  69.934,    73.777,  77.803, &
       82.015,  86.423,  91.034,  95.855,    100.89,  106.16, &
       111.66,  117.40,  123.40,  129.65,    136.17,  142.98, &
       150.07,  157.46,  165.16,  173.18,    181.53,  190.22, &
       199.26,  208.67,  218.45,  228.61,    239.18,  250.16/
  public svp, mrFromVp, mrFromRHq, mrFromRHe, mrFromRH, &
       presToAlt, precipWater, adjSfcPresPoint, adjSfcPresAlt
  ! Define a version as the default
  interface mrFromRH
    module procedure mrFromRHe
  end interface
contains
  function svp(temp)

!<f90Function>**********************************************************
!
! NAME:
!
!   svp
!
! PURPOSE:
!
!   Calculate saturation vapor pressure.
!
! SYNTAX:
!
!   Results=svp(temp)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   temp  REAL  Temperature
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

    !***********************************************************************
    !* Function name: svp
    !* Purpose: Calculates saturation vapor pressure
    !* Usage: vp = svp(temp)
    !* Description: Calculates saturated vapour pressure in mb
    !*              at given temperature in K.
    !*              Value corresponds to svp over water for 
    !*              temp > 5 deg.C, to svp over ice for 
    !*              temp < 8 deg.C, and to transitional values between.
    !*              This is the formula used by met.o.11 and met.o.2b.
    !* Inputs:
    !* Var_name   Type    Description
    !* --------   ----    -----------
    !* temp       real    temperature in K
    !*
    !* Outputs:
    !* Var_name   Type    Description
    !* --------   ----    -----------
    !* svp        real    saturation vapor pressure (mb)
    !* 
    !* Includes: none
    !* Externals: none
    !*
    !* Copyright: Atmospheric and Environmental Research, Inc., 1997        
    !* Developed by Atmospheric and Environmental Research, Inc.            
    !***********************************************************************
    real, intent(in) :: temp
    real :: svp
    real :: tt, t0, e0, e1
    integer :: ind

    ! 183.15 = 273.15 - 90
    tt=temp-183.15	
    if (tt.le.0.) then
       svp=estab(1)
    else
       ind=int(tt)+1
       ind=min(ind,155)
       t0=ind-1
       e0=estab(ind)
       e1=estab(ind+1)
       svp=e0+(tt-t0)*(e1-e0)
    endif

    return
  end function svp

  !--
  function mrFromVp(vp,pres)

!<f90Function>**********************************************************
!
! NAME:
!
!   mrFromVp
!
! PURPOSE:
!
!   Calculates mixing ratio from vapor pressure and pressure.
!
! SYNTAX:
!
!   Results=mrFromVp(vp, pres)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   vp        REAL  vapor pressure
!   pres      REAL  Pressure
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

    !***********************************************************************
    !* Function name: mrFromVp
    !* Purpose: Calculates mixing ratio from vapor pressure and pressure
    !* Usage: mixingRatio = mrFromVp(vp,pres)
    !* Description: Applies a lower limit to the molecular fraction
    !*              of dry air.
    !* Inputs:
    !* Var_name   Type    Description
    !* --------   ----    -----------
    !* vp         real    vapor pressure (mb)
    !* pres       real    total pressure (mb)
    !*
    !* Outputs:
    !* Var_name   Type    Description
    !* --------   ----    -----------
    !* mrFromVp   real    mixing ratio (g/g)
    !* 
    !* Includes: none
    !* Externals: none
    !*
    !* Copyright: Atmospheric and Environmental Research, Inc., 2009        
    !* Developed by Atmospheric and Environmental Research, Inc.            
    !***********************************************************************
    ! Input arguments:
    real, intent(in) :: vp
    real, intent(in) :: pres
    real  :: mrFromVp
    ! Constants:
    real, parameter :: dryFracLimit=0.9  ! Minimum molecular fraction of dry air

    if ((pres-vp)/vp > dryFracLimit) then
       mrFromVp=ratioMwt*vp/(pres-vp)
    else
       mrFromVp=ratioMwt/dryFracLimit  ! Avoid zero or negative denominator
    endif

    return
  end function mrFromVp

  !--
  function mrFromRHq(rh,temp,pres)

!<f90Function>**********************************************************
!
! NAME:
!
!   mrFromRHq
!
! PURPOSE:
!
!   Calculates mixing ratio (q) from relative humidity (RH), 
!   based on the definition RH=q/qsat, where sat indicates saturation.
!
! SYNTAX:
!
!   Results=mrFromRHq(rh, temp, pres)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   rh         REAL  Relative humidity
!   temp       REAL  Temperature
!   pres       REAL  Pressure
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

    !***********************************************************************
    !* Function name: mrFromRHq
    !* Purpose: Calculates mixing ratio (q) from relative humidity (RH), based 
    !* on the definition RH=q/qsat, where sat indicates saturation.
    !* Usage: mixingRatio = mrFromRHq(rh,temp,pres)
    !* Inputs:
    !* Var_name   Type    Description
    !* --------   ----    -----------
    !* rh         real    relative humidity, as a fraction (not %)
    !* temp       real    temperature in K
    !* pres       real    pressure in mb
    !*
    !* Outputs:
    !* Var_name   Type    Description
    !* --------   ----    -----------
    !* mrFromRHq  real    mixing ratio (g/g)
    !* 
    !* Includes: none
    !* Externals: none
    !*
    !* Copyright: Atmospheric and Environmental Research, Inc., 1997        
    !* Developed by Atmospheric and Environmental Research, Inc.            
    !***********************************************************************
    ! Input arguments:
    real, intent(in) :: rh
    real, intent(in) :: temp
    real, intent(in) :: pres
    real  :: mrFromRHq 
    ! Local variables:
    real  :: vpSat,qSat

    vpSat=svp(temp)
    qSat=mrFromVp(vpSat,pres)
    mrFromRHq=rh*qSat

    return
  end function mrFromRHq

  !--
  function mrFromRHe(rh,temp,pres)

!<f90Function>**********************************************************
!
! NAME:
!
!   mrFromRHe
!
! PURPOSE:
!
!   Calculates mixing ratio (q) from relative humidity (RH), based on the 
!   definition RH=e/esat, where e is vapor pressure and sat indicates 
!   saturation.
!
! SYNTAX:
!
!   Results=mrFromRHe(rh, temp, pres)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   rh         REAL  Relative humidity
!   temp       REAL  Temperature
!   pres       REAL  Pressure
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

    !***********************************************************************
    !* Function name: mrFromRHe
    !* Purpose: Calculates mixing ratio (q) from relative humidity (RH), based 
    !* on the definition RH=e/esat, where e is vapor pressure and sat indicates
    !* saturation.
    !* Usage: mixingRatio = mrFromRHe(rh,temp,pres)
    !* Inputs:
    !* Var_name   Type    Description
    !* --------   ----    -----------
    !* rh         real    relative humidity, as a fraction (not %)
    !* temp       real    temperature in K
    !* pres       real    pressure in mb
    !*
    !* Outputs:
    !* Var_name   Type    Description
    !* --------   ----    -----------
    !* mrFromRHe  real    mixing ratio (g/g)
    !* 
    !* Includes: none
    !* Externals: none
    !*
    !* Copyright: Atmospheric and Environmental Research, Inc., 1997        
    !* Developed by Atmospheric and Environmental Research, Inc.            
    !***********************************************************************
    ! Input arguments:
    real, intent(in) :: rh
    real, intent(in) :: temp
    real, intent(in) :: pres
    real  :: mrFromRHe 
    ! Local variables:
    real :: vpSat,vp

    vpSat=svp(temp)
    vp=rh*vpSat
    mrFromRHe=mrFromVp(vp,pres)

    return
  end function mrFromRHe

  !--
  subroutine presToAlt(lat,lon,temp,h2o,pref,nref,surfalt,psfc,tsfc,hsfc,zref)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   presToAlt
!
! PURPOSE:
!
!   Calculate altitude profile using hydrostatic equation.
!
! SYNTAX:
!
!   CALL presToAlt(lat, lon, temp, h2o, pref, nref, surfalt, psfc, 
!      tsfc, hsfc, zref)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   lat      REAL     Latitude
!   lon      REAL     Longitude
!   temp     REAL     Temperature
!   h2o      REAL     water vapor amount
!   pref     REAL     Pressure on atmospheric grid levels
!   nref     INTEGER  index for the pressure level below surface
!   surfalt  REAL     surface altitude
!   psfc     REAL     Surface pressure
!   tsfc     REAL     Surface-level air temperature
!   hsfc     REAL     surface moisture
!   
!   INPUTS/OUTPUTS:
!   
!   zref     REAL     array of altitudes at pressure levels
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !*************************************************************************
    !* Subroutine name: presToAlt
    !* Purpose: Calculate altitude profile using hydrostatic equation
    !* Usage: call presToAlt(lat,lon,temp,h2o,pref,nlev,nref,surfalt,psfc,
    !*                       tsfc,hsfc,zref)
    !* Description: Use the hypsometric equation to compute altitude profile.
    !* Inputs:
    !* Var_name   Type     Description
    !* --------   ----     -----------
    !* lat        real     latitude  (to calculate gravity)
    !* lon        real     longitude (to calculate gravity)
    !* temp       real     temperature profile
    !* h2o        real     moisture profile [g/g]
    !* pref       real     pressure levels  [mbar]
    !* nref       integer  index for the pressure level below surface
    !* surfalt    real     surface altitude [km]
    !* psfc       real     surface pressure [mbar]
    !* tsfc       real     surface temperature
    !* hsfc       real     surface moisture [g/g]
    !*                     
    !* Outputs:
    !* Var_name   Type     Description
    !* --------   ----     -----------
    !* zref       real     vector of altitudes at each pressure level [km]
    !*
    !* Copyright: Atmospheric and Environmental Research, Inc., 2001-2005
    !* Developed by Atmospheric and Environmental Research, Inc.
    !*************************************************************************
    integer,            intent(in)     :: nref
    real,               intent(in)     :: lat,lon,surfalt,psfc,tsfc,hsfc
    real, dimension(:), intent(in)     :: temp, h2o, pref
    real, dimension(:), intent(inout)  :: zref

    !--Local variables
    integer                            :: i,nstart, nlev
    real                               :: grav,gor,tv1,tv2,tv,gort,dzp

    ! subroutine to convert pressure profile to altitude levels
    ! note that zref is in kilometers

    ! largest index (nref) is the 1st level below surface

    !--Calculation below the surface is an extrapolation based on tsfc and hsfc

    tv = tsfc * (1.+ratioMwt*hsfc)  !--Virtual temperature for all levels below surface
    gor=grav(surfalt,0.,0.,lat,lon)/Rd
    gort = gor/tv
    dzp = alog(pref(nref)/psfc) / gort
    zref(nref) = surfalt - dzp/1000. !--Level just below the surface
    nlev = size(pref)

    do i = nref+1,nlev
       gor=grav(zref(i-1),0.,0.,lat,lon)/Rd
       gort = gor/tv
       dzp = alog(pref(i)/pref(i-1)) / gort
       zref(i) = zref(i-1) - dzp/1000.
    end do

    !--Calculation above the surface is based on valid profile data.

    if(psfc < pref(nlev))then
       nstart=nref-1
    else
       nstart=nref
    endif

    tv1 = tsfc * (1.+ratioMwt*hsfc)
    tv2 = temp(nstart) * (1.+ratioMwt*h2o(nstart))
    tv = (tv1+tv2)*0.5
    gor=grav(surfalt,0.,0.,lat,lon)/Rd
    gort = gor/tv
    dzp = alog(psfc/pref(nstart)) / gort
    zref(nstart) = surfalt + dzp/1000. !--Level just above the surface

    do i = nstart,2,-1
       gor=grav(zref(i),0.,0.,lat,lon)/Rd
       tv1 = temp(i) * (1.+ratioMwt*h2o(i))
       tv2 = temp(i-1) * (1.+ratioMwt*h2o(i-1))

       tv = (tv1+tv2)*0.5
       gort = gor/tv

       dzp = alog(pref(i)/pref(i-1)) / gort
       zref(i-1) = zref(i) + dzp/1000.
    end do

    return
  end subroutine presToAlt

  !--
  subroutine precipWater(pref,psfc,h2o,pw)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   precipWater
!
! PURPOSE:
!
!   Compute precipitable water.
!
! SYNTAX:
!
!   CALL precipWater(pref, psfc, h2o, pw)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   pref  REAL  Pressure on atmospheric grid levels
!   psfc  REAL  Surface pressure
!   h2o   REAL  water vapor amount
!   
!   INPUTS/OUTPUTS:
!   
!   pw    REAL  precipitable water
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !**********************************************************************
    !* Function name: precipWater
    !* Purpose: Computes precipitable water
    !* Usage: call precipWater(pref,psfc,h2o,pw)
    !* Description: 
    !* Inputs:
    !*  Var_name  Type      Description
    !*  --------  ----      -----------
    !*  pref      real      pressure profile (mb)
    !*  psfc      real      surface pressure (mb)
    !*  h2o       real      mixing ratio (g/g)
    !*    
    !* Outputs:
    !*  Var_name  Type      Description
    !*  --------  ----      -----------
    !*  pw        real      precipitable water (kg/m2)
    !*
    !* History:
    !*     2000, First Version ......................................AER
    !*     May 17, 2005, Made independent of state vector............AER
    !*
    !* Copyright, AER, Inc., 2000-2005
    !**********************************************************************
    !dummy arguments
    real, dimension(:), intent(in)    :: h2o,pref
    real,               intent(in)    :: psfc
    real,               intent(inout) :: pw

    !local variables
    real :: dp,QINT,sfac,p0,p1,q0,q1,z,pu,pl,qu,ql,wu,wl,hp,x0,zeta,alza,xint
    integer :: i,Nlast,NSurf
    logical :: found_index
    data dp/0.01/

    if (size(pref) .ne. size(h2o)) then
       print*, "err[precipWater]: inconsistent levels between "// &
            "pressure profile and mixing ratio."
       call errorHalt(1)
    endif

    found_index = .false.
    do i=size(pref)-1,1,-1
       if (pref(i) .lt. psfc-dp) then
          found_index = .true.
          Nlast=i
          exit
       endif
    end do

    if (.not. found_index) then
       print *, "err[precipWater]: input pressure out of bounds"
       print *, "Pressure Grids: ", pref(1:size(pref))
       print *, "Psfc: ", psfc
       call errorHalt(1)
    else
       NSurf=Nlast+1
    endif

    QINT = 0.0
    P0=psfc
    Q0=h2o(NSurf)
    if (Q0 < 0.) Q0=0.000001
    if (pref(NSurf) /= psfc) then
       z=alog(h2o(Nlast)/Q0)/alog(pref(NLast)/pref(NSurf))
       q0=Q0*(psfc/pref(NSurf))**z
    endif

    do I=Nlast,1,-1
       P1=pref(i)
       Q1=h2o(i)
       qu=Q1
       ql=Q0
       wu=qu/(1.+qu)
       wl=ql/(1.+ql)
       hp=ALOG(P0/P1)
       x0=P0*wl*hp
       zeta=P1*wu/(P0*wl)
       if(abs(zeta-1.) > 0.00001)then
          alza=ALOG(zeta)
          xint=x0*(zeta-1.)/alza
       else
          xint=x0*2./(3.-zeta)
       endif
       QINT=QINT+xint*100./g0
       P0=P1
       Q0=Q1
    enddo

    pw = QINT

    return
  end subroutine precipWater

  !--
  subroutine adjSfcPresPoint(tvprof,press,Z0,Ps0,lat,lon,Elev,gElev,Ps)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   adjSfcPresPoint
!
! PURPOSE:
!
!   Generate altitude profile.
!
! SYNTAX:
!
!   CALL adjSfcPresPoint(tvprof, press, Z0, Ps0, lat, lon, Elev, gElev, Ps)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   tvprof  REAL  virtual temperature profile
!   press   REAL  Pressure on atmospheric levels
!   Z0      REAL  NWP surface elevation
!   Ps0     REAL  NWP surface pressure
!   lat     REAL  Latitude
!   lon     REAL  Longitude
!   Elev    REAL  Elevation of surface
!   gElev   REAL  local gravity
!   
!   INPUTS/OUTPUTS:
!   
!   Ps      REAL  height-corrected surface pressure
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !*************************************************************************
    !* Subroutine name: adjSfcPresPoint
    !* Purpose: Generates altitude profile 
    !* Usage: call adjSfcPresPoint(tvprof,press,z0,ps0,xlat,xlon,Elev,gElev,Ps)
    !* Description: Calculates height-corrected surface pressure at given lat, lon. 
    !* Inputs:
    !* Var_Name      Type       Description
    !* --------      ----       -----------
    !* tvprof        real       virtual temperature profile
    !* press         real       pressure profile
    !* Z0            real       NWP surface elevation
    !* Ps0           real       NWP surface pressure
    !* lat           real       latitude  in degrees
    !* lon           real       longitude in degrees
    !* Elev          real       local elevation
    !* gElev         real       local gravity
    !*
    !* Outputs:
    !* Var_Name      Type       Description
    !* --------      ----       -----------
    !* Ps            real       height-corrected surface pressure
    !*
    !* Externals:
    !* Name         Description
    !* ----         -----------
    !* grav         computes elevation corrected accel. due to grav.
    !*
    !* History:
    !*   May 26, 2005, imported from the independent subroutine, 
    !*                 adjSfcPressPoint, under rslib/util/sensor_model/
    !*                 .................................... Y. He, AER
    !*
    !* Copyright: 2004-2005, AER, Inc. All rights reserved.
    !* Developed by Atmospheric and Environmental Research, Inc.   
    !************************************************************************
    !dummy arguments
    real, dimension(:), intent(in)    :: tvprof,press
    real,               intent(in)    :: Z0, Ps0, lat, lon, Elev, gElev
    real,               intent(inout) :: Ps

    !-- local variables
    integer :: i,j,isurf
    real :: dz,Tv,tv1,tv2,hbar,p0,zref,diff,p,g,a,grav,tvbar
    real, dimension(0:size(press)) :: z

    !-- Stop if pressure profile does not match temperature profile
    if (size(press) /= size(tvprof)) then
       print*, "err[MetFunctions::adjSurfPressPoint]: mismatch #levels."
       print*, "pressure: ", size(press), "....temperature: ",size(tvprof)
       call errorHalt(1)
    endif

    !-- Determine first pressure level above NWP surface
    do j=1,size(press)
       if (Ps0.gt.press(j)) isurf=j
    enddo

    z(0) = Z0   !-- Initialize z array    

    !-- Calculate elevation for first two pressure levels above surface
    g    = grav(Z0,0.,0.,lat,lon)
    a    = 0.5*alog(Ps0*press(isurf))
    p    = (a-alog(press(isurf-1)))/alog(press(isurf)/press(isurf-1))
    Tv   = (1.-p)*(tvprof(isurf-1))+p*tvprof(isurf)
    dz   = ((Rd*Tv)/g)*(alog(Ps0/press(isurf))) 
    z(1) = z(0)+dz
    g    = grav(z(1),0.,0.,lat,lon)
    Tv   = (tvprof(isurf)+tvprof(isurf-1))*0.5
    dz   = ((Rd*Tv)/g)*(alog(press(isurf)/press(isurf-1))) 
    z(2) = z(1)+dz

    if(Elev.lt.z(2)) then ! Enough information to determine the Ps

       diff   = z(2)-z(1)
       tv1    = tvprof(isurf)
       tv2    = tvprof(isurf-1)
       if(Elev.lt.z(1)) then
          Hbar   = (Elev+Z0)*0.5
          p      = (Hbar-z(1))/diff
          P0     = Ps0
          Zref   = Z0
       else
          Hbar   = (Elev+z(1))*0.5
          p      = (Hbar-z(1))/diff
          P0     = press(isurf)
          Zref   = z(1)
       endif
       Tvbar  = (1.0-p)*tv1+p*tv2
       Ps     = P0*exp((gElev/(Rd*Tvbar))*(Zref-Elev))
       return

    else ! Need more info so integate up to the bounding levels

       do i=3,size(press)
          g    = grav(z(i-1),0.,0.,lat,lon)
          Tv   = (tvprof(isurf-(i-2))+tvprof(isurf-(i-1)))*0.5
          dz   = ((Rd*Tv)/g)*(alog(press(isurf-(i-2))/press(isurf-(i-1)))) 
          z(i) = z(i-1)+dz
          if(Elev.lt.z(i)) then       
             diff   = z(i)-z(i-1)
             tv1    = tvprof(isurf-(i-2))
             tv2    = tvprof(isurf-(i-1))
             Hbar   = (Elev+z(i-1))*0.5
             p      = (Hbar-z(i-1))/diff
             Tvbar  = (1.0-p)*tv1+p*tv2
             g=grav(Hbar,0.,0.,lat,lon) ! elevation corrected g
             Ps     = press(isurf-(i-2))*exp((g/(Rd*Tvbar))*(z(i-1)-Elev))
             return
          endif
       enddo
    endif
    return 
  end subroutine adjSfcPresPoint

  !--
  subroutine adjSfcPresAlt(latFov,lonFov,elevFov,profGrid4,latGrid4,lonGrid4, &
       elevGrid4,presGrid,IGgrid,NGgrid,wt,molID,itype,psfcFov)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   adjSfcPresAlt
!
! PURPOSE:
!
!   Adjust NWP-derived surface pressure for local topography.
!
! SYNTAX:
!
!   CALL adjSfcPresAlt(latFov, lonFov, elevFov, profGrid4, latGrid4, 
!      lonGrid4, elevGrid4, presGrid, IGgrid, NGgrid, wt, molID, 
!      itype, psfcFov)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   latFov     REAL                Latitude at FOV
!   lonFov     REAL                Longitude at FOV
!   elevFov    REAL                Elevation of surface at FOV
!   profGrid4  REAL                Atmospheric profile data from surrounding 
!                                  grid points 
!   latGrid4   REAL                Latitudes of surrounding grid points
!   lonGrid4   REAL                Longitudes of surrounding grid points
!   elevGrid4  REAL                Elevations of surfaceat the surrounding 
!                                  grid points 
!   presGrid   REAL                pressure grid
!   IGgrid     TYPE(STATEINDEX_T)  Starting indices for sections of 
!                                  geophysical state vector on grid 
!   NGgrid     TYPE(STATEINDEX_T)  Number of elements for sections of 
!                                  geophysical state vector on grid 
!   wt         REAL                bi-linear interpolation weights for four 
!                                  grid points surrounding FOV location 
!   molID      INTEGER             List of IDs of relevant molecular species
!   itype      INTEGER             integer flag to regulate profile input
!   
!   INPUTS/OUTPUTS:
!   
!   psfcFov    REAL                Surface pressure at FOV
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !***********************************************************************
    !* Subroutine name: 
    !*      adjSfcPresAlt
    !* Purpose: 
    !*      Adjust NWP-derived surface pressure for local topography
    !* Usage: 
    !*      call adjSfcPresAlt(latFov,lonFov,elevFov,profGrid4,latGrid4,
    !*                  lonGrid4,elevGrid4,presGrid,IGgrid,NGgrid,wt,molID,
    !*                  itype,psfcFov)
    !* Description: The subroutine maps surface pressure to xlat, xlon.
    !* Inputs:
    !* Var_Name      Type       Description
    !* --------      ----       -----------
    !* latFov        real       FOV latitude in degrees
    !* lonFov        real       FOV longitude in degrees
    !* elevFov       real       FOV elevation
    !* profGrid4     real       vector of nwp, climatology, background profiles
    !* latGrid4      real       latitudes at the horizontal grids
    !* lonGrid4      real       longitudes at the horizontal grids
    !* elevGrid4     real       elevations at the horizontal grids
    !* presGrid      real       pressure levels correponding to temperature 
    !*                          profile in profGrid
    !* IGgrid     StateIndex_t  index structure for the state vector at the 
    !*                          input gridded points
    !* NGgrid     StateIndex_t  length structure for the state vector at the 
    !*                          input gridded points
    !* wt            real       bi-linear interpolation weights for four grid 
    !*                          points surrounding FOV location
    !* molID         integer    array of molecular IDs
    !* itype         integer    1 - nwp or climatological profile input. 
    !*                          2 - no nwp or climatological input, start 
    !*                              with 1013 mb for surface pressure and 
    !*                              0 m for surface elevation
    !* Outputs:
    !* Var_Name      Type       Description
    !* --------      ----       -----------
    !* psfcFov       real       interpolated surface pressure at FOV location
    !*
    !* External Subroutines:
    !* Name             Description
    !* ----             -----------
    !* adjSfcPresPoint  computes Psfc
    !* grav             computes elevation corrected accel. due to grav.
    !*
    !* History:
    !*   May 26, 2005, imported from the independent subroutine, 
    !*                 adjSfcPressAlt, under rslib/util/sensor_model/ ...
    !*                 ....................................... Y. He, AER
    !*
    !* Copyright: 2002-2005, AER, Inc. All rights reserved.
    !* Developed by Atmospheric and Environmental Research, Inc.            
    !***********************************************************************
    use StateIndexModule
    !-- dummy arguments
    type(StateIndex_t),       intent(in)    :: IGgrid, NGgrid
    real,                     intent(in)    :: latFov,lonFov
    real,    dimension(:,:),  intent(in)    :: profGrid4
    real,    dimension(:),    intent(in)    :: latGrid4
    real,    dimension(:),    intent(in)    :: lonGrid4
    real,    dimension(:),    intent(in)    :: elevGrid4
    integer, dimension(:),    intent(in)    :: molID
    integer,                  intent(in)    :: itype
    real,    dimension(:),    intent(in)    :: presGrid
    real,                     intent(in)    :: elevFov
    real,    dimension(:),    intent(in)    :: wt
    real,                     intent(inout) :: psfcFov

    !-- local variables
    integer :: i,j
    real :: gelev,grav,latGrid,lonGrid,hnwp,psnwp,q,sig,tavg
    real, dimension(size(wt)) :: psi
    real, dimension(NGgrid%temp) :: Tv
    
    gElev=grav(elevFov,0.,0.,latFov,lonFov) !  elevFov corrected g 

    if(itype.le.1) then
       ! For each NWP grid point determine pressure at the FOV elevation
       do i=1,size(wt) 
          latGrid  = latGrid4(i)
          lonGrid  = lonGrid4(i)
          Hnwp  = elevGrid4(i)
          Psnwp  = profGrid4(IGgrid%Psfc,i)
          do j=1,NGgrid%temp
             if (j .ge. NGgrid%temp-NGgrid%mol(whereH2O(molID))+1) then
                q     = profGrid4(IGgrid%mol(whereH2O(molID))+j &
                     -(NGgrid%temp-NGgrid%mol(whereH2O(molID))+1),i)
                sig   = (1.+q/ratioMwt)/(1.+q)
                Tv(j) = sig*profGrid4(j,i)
             else 
                Tv(j) = profGrid4(j,i)
             endif
          enddo
          call adjSfcPresPoint(tv,presGrid,Hnwp,Psnwp, &
               latGrid,lonGrid,elevFov,gElev,Psi(i))
       enddo

       !-- Perform linear interpolation
       PsfcFov = 0
       do i = 1, size(wt)
          PsfcFov = PsfcFov + Psi(i)*wt(i)
       enddo

       return
    elseif(itype.eq.2)then

       ! Assume a standard profile with Tskin = 280 K and lapse rate of 6 K/km
       ! of 1000 mb corresponding to a surface elevation of 0

       psnwp=1013.
       tavg=280.-0.006*elevFov*0.5 ! mean temperature over the elevation
       psfcFov=1013.*exp(-elevFov*g0/rd/tavg)

    else

       ! Assume a surface pressure of 1000 mb 

       psfcFov=1000.

    end if

    RETURN
  END subroutine adjSfcPresAlt
end module MetFunctions
