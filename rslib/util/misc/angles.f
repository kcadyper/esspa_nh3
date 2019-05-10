c-------------------------------------------------------------------------------
c     $Name$ 
c     $Id$ 
c     Copyright AER, Inc., 2002, 2003. All rights Reserved.
c-------------------------------------------------------------------------------

      subroutine angles(ang,hite,zenang,path,smax,crh)
      use constants, ONLY : ErthRad,deg2rad
      !crh=ratio of satellite orbital radius to earth radius
      crh=(ErthRad+hite)/ErthRad
      !smax = nadir angle of view to earth edge
      smax=asin(1./crh)/deg2rad
      path=1.
      if (ang.eq.0.)then
         zenang=0.
         return
      endif
      !zenang = zenith angle of path at earth surface
      a=sin(ang*deg2rad)*crh
      if (abs(a).ge.1.) then
         print *,'err[angles::angles]:  STOP: 940'
         call errorHalt(1)
      end if
      a=asin(a)
      zenang=a/deg2rad
      !path = air mass factor
      path=path/cos(a)
      return
      end

      subroutine anglesInv(zenang,hite,ang)
      use constants, ONLY : ErthRad,deg2rad
      !crh=ratio of satellite orbital radius to earth radius
      crh=(ErthRad+hite)/ErthRad
      if (zenang.eq.0.)then
         ang=0.
         return
      endif
      !ang = nadir angle
      a=sin(zenang*deg2rad)/crh
      if (abs(a).ge.1.) then
         print *,'err[angles::angles]:  STOP: 940'
         call errorHalt(1)
      end if
      a=asin(a)
      ang=a/deg2rad
      return
      end

