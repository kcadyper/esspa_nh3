MODULE GeomSat
  
  implicit none
  public :: swathInfo, avgAlongScanSpacing,viewGeom,satSubpointHead
  public :: polarRotation, movePointing

CONTAINS

  subroutine swathInfo(satalt,az_max_in,nadir_ang_in,lat_c_in,sat_h_in,swathWidth,frel,s)
    use constants, ONLY : gravCnst,PI,deg2rad,ErthMass,ErthRad
    
    real, intent(in) :: satalt,az_max_in,nadir_ang_in,lat_c_in,sat_h_in
    real :: az_max,nadir_ang,lat_c,sat_h
    real :: radius_s,zenith,beta,swath_w,fres,lat_l,lat_r, &
         vel,orbit_p,vel_s,t_ce,eq_r,s_l,s_r
    real, intent(out) :: swathWidth,frel,s
    
    !--Convert to radians
    az_max=az_max_in*deg2rad
    nadir_ang=nadir_ang_in*deg2rad
    lat_c=lat_c_in*deg2rad
    sat_h=sat_h_in*deg2rad
    
    radius_s=satalt+ErthRad
    zenith=asin(sin(nadir_ang)*radius_s/ErthRad)
    beta=zenith-nadir_ang
    swath_w=acos(cos(2.*az_max)*(sin(beta)**2)+cos(beta)**2)
    swathWidth=swath_w*ErthRad
    fres=atan(cos(az_max)*tan(beta))
    frel=(beta-fres)*ErthRad
    lat_l=lat_c+sin(sat_h)*swath_w/2.
    lat_r=lat_c-sin(sat_h)*swath_w/2.
    vel=sqrt(gravCnst*ErthMass/radius_s/1000.)
    orbit_p=2.*PI*radius_s*1000./vel
    vel_s=vel*ErthRad/radius_s
    t_ce=frel*1000./vel_s
    eq_r=2.*PI*ErthRad*t_ce/24./3600.
    s_l=eq_r*cos(lat_l)*cos(sat_h)
    s_r=eq_r*cos(lat_r)*cos(sat_h)
    s=(s_l+s_r)*0.5
    return
  end subroutine swathInfo

  subroutine avgAlongScanSpacing(lat,lon,npos,a)
    use GeomSpheric
    use constants, ONLY : MISSING_REAL,deg2rad
    
    integer, intent(in)            :: npos
    real, dimension(:), intent(in) :: lat,lon
    real, intent(out)     :: a
    
    !--Local Variables
    integer :: nval,ipos
    real    :: lat1,lon1,lon1t,lon2t,lat2,lon2, &
         s1,c1,s2,c2,sth,cth
    
    !--Using the first valid scanline, 
    !--determine nominal CFOV spacing.
    
    nval=0
    a=0.
    do ipos=1,npos-1
       lat1=lat(ipos)
       lat2=lat(ipos+1)
       lon1=lon(ipos)
       lon2=lon(ipos+1)
       if(lat1 /= MISSING_REAL .and. lat2 /= MISSING_REAL)then
          lon1t=lonReverse(lon1,0.)
          lon2t=lonReverse(lon2,0.)
          s1=sin(lat1*deg2rad)
          c1=cos(lat1*deg2rad)
          s2=sin(lat2*deg2rad)
          c2=cos(lat2*deg2rad)
          sth=sin((lon1t-lon2t)*deg2rad)
          cth=cos((lon1t-lon2t)*deg2rad)
          a=a+gcdist(lat1,lon1t,lat2,lon2t)
          nval=nval+1
       endif
    enddo
    if(nval > 0)then
       a=a/float(nval)
    endif
    return
  end subroutine avgAlongScanSpacing
  
!==================================================================
!==================================================================
  
  SUBROUTINE viewGeom(ecross,otime,sscan,scone,slat,slon,zenang,azang,salt,calim,Incl,flat,flon)
    
    USE GeomSpheric, ONLY : mvltln,angazim
    USE constants, ONLY   : PI,ErthRad,deg2rad
    
! ************************************************************************
!   Input: 
!ecross : equator crossing of satellite on ascending orbit (ascending node)
!         (positive east) 
! otime : time since equator crossing (s)
! sscan : scan azimuth angle of observation (positive clockwise, zero fore)
! scone : cone angle of satellite scan (degrees)
!  slat : latitude of satellite
!  slon : longitude of satellite (positive east)
!  salt : satellite altitiude (km)
! calim : 
!  Incl : inclination of satellite orbit
!
!  Output:
!zenang : zenith angle of satellite from ground
! azang : azimuth angle of satellite from ground (positive clockwise,
!         zero at north)
!  flat : latitude of point viewed by satellite
!  flon : longitude of point viewed by satellite
!
!               all angles in degrees
! ************************************************************************
    
    REAL :: salt,calim
    REAL :: ecross,otime,sscan,scone,slat,slon,zenang,azang
    REAL :: flat,flon,sconrad,zenrad,cenang,sataz,sviewaz,rlat,rlon
    REAL :: xlat,xlon,s1,c1,s2,c2,sth,cth,azangr,Incl
    
    !Compute zenith angle using law of sines (noting symmetry of inv sine:
    !    inv gives zenith instead of its complement)
    sconrad=scone * deg2rad
    zenrad=asin(sin(sconrad)*(ErthRad+salt)/ErthRad)
    zenang=zenrad/deg2rad
    
    !Compute (earth centered) interior angle from satellite 
    !    subpoint to view spot
    cenang=zenrad-sconrad

    !Find lat/lon of satellite and local azimuth of satellite fore
    call satSubpointHead(ecross,otime,Incl,salt,slat,slon,sataz)

    !Find lat/lon of view spot using subroutine 'mvltln'
    sviewaz=(sataz+sscan)*deg2rad
    sviewaz=pi-sviewaz
    rlat=slat*deg2rad
    rlon=-slon*deg2rad
    call mvltln(rlat,rlon,cenang,sviewaz,xlat,xlon)
    flat=xlat/deg2rad
    flon=-xlon/deg2rad
    if(flon.le.-180.)flon=360.+flon
    if(flon.gt.180.)flon=flon-360.
    
    !Compute local azimith angle from ground to satellite
    s1=sin(xlat)
    c1=cos(xlat)
    s2=sin(rlat)
    c2=cos(rlat)
    sth=sin(xlon-rlon)
    cth=cos(xlon-rlon)
    if(cenang.ge.calim)then   ! preserve view azimuth when at nadir
       azangr=angazim(s1,c1,s2,c2,sth,cth)
    else
       azangr=sviewaz+pi
    endif
    azang=(pi-azangr)/deg2rad
    
    return
  end subroutine viewGeom


!==================================================================
!==================================================================

  SUBROUTINE satSubpointHead( Asc_Node,time_cnt,Inclination &
       ,Altitude,sub_lat,sub_lon,head_az ) 

    USE constants, ONLY : PI,ErthRad,deg2rad,ErthMass,gravCnst
    
! ************************************************************************
!      Written by Alan Lipton, AER, Inc. 8/28/98
!
!      Based on Subtrack, written by:
!
!	Michael K. Griffin   Robert P. d'Entremont   Gary B. Gustafson 
!
!	3 December 1993   SERCAA 
!
!	   This Program Computes the Subtrack point of a (sun-synchronous) 
!	Polar Orbit Given the Ascending Node of the Orbit. 
!
!     Input:
!      asc_node : ascending node; +ve East, (degrees)
!      time_cnt : time (s) since ascending equator crossing
!     Output:
!       sub_lat : latitude (degrees) of satellite subtrack point
!       sub_lon : longitude (degrees) of satellite subtrack point; +ve East
!       head_az : local azimuth angle (degrees) of satellite heading
!                 positive clockwise from North to the satellite fore 
! ************************************************************************
    
    Real Time_Cnt,Asc_Node,Tau,Tau_Prime,Lambda_Prime,Lambda_Temp

    REAL,PARAMETER :: pio2=PI/2.    
    REAL           :: Inclination,Altitude
    REAL    :: sub_lat,sub_lon,head_az
    REAL    :: period,period_1,period_3,pi_factor,earth_rotation_factor
    REAL    :: oinclin,sin_inc,cos_inc,cos_tau_prime,earth_spin
    REAL    :: phi_temp,oinc,bv,bh,azb,alpha
    REAL    :: tot_rad,periodsq           
    
!	Period of the Satellite Orbit 
    tot_rad=(ErthRad+altitude)*1.e3
    periodsq=4.*Pi**2*tot_rad**3/(gravcnst*ErthMass)
    
    Period=sqrt(periodsq)
    Period_1 = 0.25 * Period 
    Period_3 = 0.75 * Period 
    if(time_cnt.lt.0..or.time_cnt.gt.Period+1.)goto 910
 
    Pi_Factor = 2.00 * Pi / Period 

!	Earth's Rotation (Added to Longitude, + West), Radians/sec 
!
    Earth_Rotation_Factor = 2. * Pi / ( 24. * 3600. ) 
!
!	Inclination of the Satellite Orbit 
!
    oinclin=Inclination
    If( oinclin .Gt. 90.0 ) then
       oinclin = 180.0 - oinclin 
    else
       goto 920
    endif
    oinclin = oinclin * deg2rad 
    Sin_Inc = Sin( oinclin ) 
    Cos_Inc = Cos( oinclin ) 
!
!	Orbital Latitude in Radians (Based on Time from Ascending Node) 
!
    Tau = Pi_Factor * Time_Cnt 
!
!	First and Last Quarter of Orbit (Northbound)
!
    If( Time_Cnt .Le. Period_1 .or. Time_Cnt .ge. Period_3) Then
       
       Tau_Prime = Pi / 2. - Tau 
!
!	Middle half of Orbit (Southbound) 
!
    Else 
       
       Tau_Prime = Tau - Pi / 2. 
       
    End If
    
    Cos_Tau_Prime = Cos( Tau_Prime ) 
!
!	Westward Component of the Earth's Rotation 
!
    Earth_Spin = Earth_Rotation_Factor * Time_Cnt 
!
    Phi_Temp = Sin_Inc * Cos_Tau_Prime 
    Phi_Temp = Min( 1., Max( -1., Phi_Temp ) ) 
    sub_lat = Asin( Phi_Temp ) 
    
    Lambda_Temp = -( Cos_Inc * Cos_Tau_Prime ) / Cos( sub_lat ) 
    
    Lambda_Prime = Min( 1., Max( -1., Lambda_Temp ) ) 

    !	First Quarter of Orbit (Northbound) 

    If( Time_Cnt .Le. Period_1 .or. Time_Cnt .ge. Period_3) Then
 
       sub_lon = Acos( Lambda_Prime ) - Pi/2. + Earth_Spin 
       !
!	Second Quarter of Orbit (Southbound) 
!
    Else 
       
       sub_lon = 1.5 * Pi - Acos( Lambda_Prime ) + Earth_Spin 
       
    End If
!
!	Latitude in Degrees, +ve North 
!
    sub_lat = sub_lat / deg2rad 
    !
!	Convert From +ve West to +ve East 
!
    sub_lon = -sub_lon / deg2rad 
    sub_lon = sub_lon + asc_node
    If( sub_lon .Lt. -180.0 ) sub_lon = sub_lon + 360. 
    If( sub_lon .gt. 180.0 ) sub_lon = sub_lon - 360. 
!
!      Compute heading azimuth
!      Split up the computation to avoid numerical instability and
!      to ensure the correct sign.
    oinc=inclination*deg2rad
    bv=-sin(oinc)*cos(tau)
    bh=cos(oinc)
    if(abs(bv).le.abs(bh))then
       if(abs(bh).lt.1.e-7)then
          azb=0.
       else
          alpha=atan(bv/bh)
          if(bh.ge.0.)then
             azb=pio2+alpha
          else
             azb=pi+pio2+alpha
          endif
       endif
    else
       if(abs(bv).lt.1.e-7)then
          azb=0.
       else
          alpha=atan(bh/bv)
          if(bv.ge.0.)then
             azb=pio2+pio2-alpha
          else
             azb=-alpha
          endif
       endif
    endif
    if(azb.gt.pi)azb=azb-2.*pi
    head_az=azb/deg2rad
    
    Return 
    
910 write(*,915)time_cnt,Period
915 format(/1x,'orbit time out of range; time_cnt,Period:' &
         ,2f10.2)
    print *,'err[GeomSat::satSubpointHead]: STOP ERROR 910'
    call errorHalt(1)
920 write(*,925)Inclination
925 format(/1x,'Inclination must be greater than 90 for sun syc;' &
         ,'Inclination:',f10.1)
    print *,'err[GeomSat::satSubpointHead]: STOP ERROR 920'
    call errorHalt(1)
    
  End subroutine satSubpointHead

!==================================================================
!==================================================================

  SUBROUTINE polarRotation(nadir,azim,rot &
       ,roll_a2sens,pitch_a2sens &
       ,roll_sens2sc,pitch_sens2sc &
       ,roll_sc2orb,pitch_sc2orb &
     !  OUTPUT:
       ,rot_out)
    
    USE constants, ONLY : pi
    
    REAL :: nadir,azim,rot
    REAL :: nadircmp,azimrad,nadrad
    REAL :: cosbr,cosbp
    REAL :: pra_r_a2sens,pra_p_a2sens
    REAL :: pra_r_sens2sc,pra_p_sens2sc
    REAL :: pra_r_sc2orb,pra_p_sc2orb
    REAL :: roll_a2sens,pitch_a2sens,roll_sens2sc,pitch_sens2sc
    REAL :: roll_sc2orb,pitch_sc2orb
    REAL :: rot_out
    
    nadircmp = (pi/2) - (nadir * (pi/180))
    
    azimrad = azim * (pi / 180)
    nadrad  = nadircmp * (pi / 180)
    
    cosbr = COS(0-azimrad)*cos(-1*nadrad)
    cosbp = COS((-1*(pi/2))-azimrad)*COS(-1*nadrad)
    
    pra_r_a2sens  = roll_a2sens * cosbr
    pra_p_a2sens  = pitch_a2sens * cosbp
    
    pra_r_sens2sc = roll_sens2sc * cosbr
    pra_p_sens2sc = pitch_sens2sc * cosbp
    
    pra_r_sc2orb  = roll_sc2orb * cosbr
    pra_p_sc2orb  = pitch_sc2orb * cosbp
    
    rot_out = rot + pra_r_a2sens + pra_p_a2sens &
         + pra_r_sens2sc + pra_p_sens2sc &
         + pra_r_sc2orb + pra_p_sc2orb
    
  END SUBROUTINE polarRotation
  
!==================================================================
!==================================================================

  SUBROUTINE movePointing(nadirB,azimuthB,roll,pitch,yaw,nadirR,azimuthR)
    ! Moves (nadir,azimuth) pointing coordinates from body coordinate
    ! frame to reference frame when body is rotated thru roll-pitch-yaw
    ! angles relative to reference.
    !---In/Out variables
    REAL                 :: nadirB,azimuthB,roll,pitch,yaw
    REAL                 :: nadirR,azimuthR
    !---Local variables
    REAL,    DIMENSION(3) :: uBody,uRef,uSphere
    
    uSphere(1) = 1
    uSphere(2) = nadirB
    uSphere(3) = azimuthB
    call sphere2rect(uSphere,uBody)
    call moveRPY(uBody,uRef,roll,pitch,yaw)
    call rect2sphere(uRef,uSphere)
    nadirR = uSphere(2)
    azimuthR = uSphere(3)

    RETURN
  end SUBROUTINE movePointing

!==================================================================
!==================================================================

  SUBROUTINE moveRPY(uBody,uRef,roll,pitch,yaw)
    ! Transforms from body-centric coordinates (uBody=[x,y,z]) to 
    ! reference-system coordinates (uRef) when body has been rotated
    ! away from reference thru standard roll-pitch-yaw angles. 
    use constants, ONLY : deg2rad
    !---In/Out variables
    REAL,    DIMENSION(3) :: uBody,uRef
    REAL                  :: roll,pitch,yaw
    !---Local variables
    REAL,    DIMENSION(3,3) :: body2ref
    REAL                  :: rR,rP,rY

    rR = deg2rad*roll
    rP = deg2rad*pitch
    rY = deg2rad*yaw

    body2ref(1,1) = cos(rY)*cos(rP)-sin(rY)*sin(rR)*sin(rP)
    body2ref(1,2) = -sin(rY)*cos(rR)
    body2ref(1,3) = cos(rY)*sin(rP)+sin(rY)*sin(rR)*cos(rP)
    body2ref(2,1) = sin(rY)*cos(rP)+cos(rY)*sin(rR)*sin(rP)
    body2ref(2,2) = cos(rY)*cos(rR)
    body2ref(2,3) = sin(rY)*sin(rP)-cos(rY)*sin(rR)*cos(rP)
    body2ref(3,1) = -cos(rR)*sin(rP)
    body2ref(3,2) = sin(rR)
    body2ref(3,3) = cos(rR)*cos(rP)

    uRef = matmul(body2ref,uBody)
    
    RETURN
  END SUBROUTINE moveRPY

!==================================================================
!==================================================================

  SUBROUTINE sphere2rect(uSphere,uRect)
    ! Transforms from spherical to rectangular coordinates 
    use constants, ONLY : deg2rad
    !---In/Out variables
    REAL,    DIMENSION(3) :: uSphere,uRect
    !---Local variables
    REAL                  :: r, theta, phi

    r = uSphere(1)
    theta = deg2rad*uSphere(2)
    phi = deg2rad*uSphere(3)

    uRect(1) = r*sin(theta)*cos(phi)
    uRect(2) = r*sin(theta)*sin(phi)
    uRect(3) = r*cos(theta)
    
    RETURN
  END SUBROUTINE sphere2rect

!==================================================================
!==================================================================
  
  SUBROUTINE rect2sphere(uRect,uSphere)
    ! Transforms from rectangular to spherical coordinates
    use constants, ONLY : rad2deg
    !---In/Out variables
    REAL,    DIMENSION(3) :: uSphere,uRect
    !---Local variables
    ! REAL                  :: r, theta, phi

    ! r
    uSphere(1) = sqrt(uRect(1)**2 + uRect(2)**2 + uRect(3)**2)
    ! theta
    uSphere(2) = rad2deg*atan(sqrt(uRect(1)**2 + uRect(2)**2)/uRect(3))
    ! phi
    uSphere(3) = rad2deg*atan2(uRect(2),uRect(1))
    
    RETURN
  END SUBROUTINE rect2sphere
  
end MODULE GeomSat
