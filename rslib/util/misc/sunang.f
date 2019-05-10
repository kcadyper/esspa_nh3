c-------------------------------------------------------------------------------
c     $Name$ 
c     $Id$ 
c     Copyright AER, Inc., 2002, 2003. All rights Reserved.
c-------------------------------------------------------------------------------

      subroutine sunang(mode,time,alat,alon,sundec,sunzang,dosun)
c***********************************************************************
c* Function name: sunang
c* Purpose: Calculates solar zenith angle
c* Usage: call sunang(mode,time,alat,alon,sundec,sunzang,dosun)
c* Description: The function evaluates the solar zenith angle
c*              by the given time of the day, latitude and 
c*              longitude. 
c*              This routine is adapted from the pathfinder 
c*              routine sunang.f. This is approximate, but adequate.
c* Inputs:
c* Var_Name   Type          Description
c* --------   ----          -----------
c* mode       integer       0=solve for solar angle 
c* time       integer       (year,month,day,hour,minute,second)
c*
c* Outputs:
c* Var_Name      Type       Description
c* --------      ----       -----------
c* sundec        real       sub-solar latitude.
c* sunzang       real       solar zenith angle
c* dosun         logic      T=include the solar term, F=otherwise.
c*
c* Common blocks: none
c* Includes: none
c* Externals: none
c* 
c* Copyright: n.a.
c* Modified by Atmospheric and Environmental Research, Inc., 1998.
c***********************************************************************
      implicit none

      logical*4 dosun
      integer mode,yy,mm,dd,hh,mn,ss, time(*)
      real    alat, alon, sunzang, sundec

      real    declmax,daymaxd,dayspyr,pi2,rad2deg
      parameter (daymaxd=173.0)  ! day number of maximum declination
      parameter (declmax=23.45)  ! maximum solar declination(daymaxd)
      parameter (dayspyr=365.0)  ! solar declination cycle 
      parameter (pi2=6.2831853)  ! 2*pi
      parameter (rad2deg=57.2957795)

      integer i
      real    days(12)         ! julian day number of dd=0 for each month
      real    ztime, jday, alonx
      real    dec, sindec, cosdec, sinlat, coslat
      real    ha, cosha, cossun

      data  days/   0.0,  31.0,  59.0,  90.0, 120.0, 151.0,
     &            181.0, 212.0, 243.0, 273.0, 304.0, 334.0 /

      yy = time(1)
      mm = time(2)
      dd = time(3)
      hh = time(4)
      mn = time(5)
      ss = time(6)

      ztime = float(hh) + float(mn)/60.0 + float(ss)/3600.0  ! ut time

      jday  = days(mm) + float(dd)                  ! Julian day
      i = yy/4
      i = i*4                                  
      if(i.eq.yy.and.mm.eq.2) jday = jday + 1       ! leap year

      dec = declmax*cos(pi2*(jday-daymaxd)/dayspyr) ! sub-solar latitude
      sundec = dec
      cosdec = cos(dec/rad2deg)                     ! sine of sub-solar lat.
      sindec = sin(dec/rad2deg)                     ! cosine of sub-solar lat.

      sinlat = sin(alat/rad2deg)                    ! sine of latitude
      coslat = cos(alat/rad2deg)                    ! cosine of latitude

      ztime = ztime - 12.0                          ! 12:00 noon is 0

c--------------------------------------------
c     east longitude
c     e.g., -75 longitude = 5h w of greenwich
c          at 7 ut sun is at ha=0
c          at 8 ut sun is at ha=1
c--------------------------------------------

      if (mode.eq.0) then        ! solve for sunang
         alonx = alon
         if(alonx.gt.180.0) alonx = alonx - 360.0  ! 180 <= alonx <= 180
         
         ha = ztime + 24.0*alonx/360.0
         cosha = cos(pi2*ha/24.0)
         cossun = sinlat*sindec + coslat*cosdec*cosha
         sunzang = rad2deg*acos(cossun) ! solar zenith angle
         if(sunzang.le.85.0) then
            dosun = .true.
         else
            dosun = .false.
         endif
      else      !<-- solve for longitude
         cossun = cos(sunzang/rad2deg)    ! mu
         cosha = (cossun - sinlat*sindec)/(coslat*cosdec)
         ha = 24.0*acos(cosha)/pi2 ! this can be + or -
         alon = 360.0*(ha - ztime)/24.0

         alon = mod(alon,360.0)
         if(alon.lt.-180.0) alon = alon + 360.0
         if(alon.gt. 180.0) alon = alon - 360.0
c---------------------------------------------------
c-- Note: There are always 2 solutions at a given 
c--       ut time.  Without knowledge of the correct 
c--       hemisphere it is impossible to distinquish
c--       between these solutions.
c---------------------------------------------------
      endif

      return
      end





