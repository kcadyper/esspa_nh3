module AngleUtil

  use OSSPhysicalConstant, ONLY: ErthRad,deg2rad

  implicit none
  private
  public:: zenithFromNadir,nadirFromZenith

contains
   function zenithFromNadir(nadir,satAlt)
      real,    intent(in)  :: nadir  ! Satellite scan nadir angle (degrees)
      real,    intent(in)  :: satAlt ! Satellite altitude (km)
      real        :: zenithFromNadir ! Zenith angle at earth surface(degrees)

      real :: crh      ! ratio of satellite orbital radius to earth radius
      real :: sinZen
      real :: zenRad

      crh=(ErthRad+satAlt)/ErthRad
      if (nadir == 0.)then
         zenithFromNadir=0.
         return
      endif
      sinZen=sin(nadir*deg2rad)*crh
      if (abs(sinZen) >= 1.) then
         print *,'err[AngleUtil::zenithFromNadir]'
         call exit(1)
      end if
      zenRad=asin(sinZen)
      zenithFromNadir=zenRad/deg2rad
      return

   end function zenithFromNadir

   function nadirFromZenith(zenAng,satAlt)
      real,    intent(in)  :: zenAng  ! Zenith angle at earth surface(degrees)
      real,    intent(in)  :: satAlt ! Satellite altitude (km)
      real        :: nadirFromZenith ! Satellite scan nadir angle (degrees)

      real :: crh      ! ratio of satellite orbital radius to earth radius
      real :: sinNadir
      real :: nadirRad

      crh=(ErthRad+satAlt)/ErthRad
      if (zenAng == 0.)then
         nadirFromZenith=0.
         return
      endif
      sinNadir=sin(zenAng*deg2rad)/crh
      if (abs(sinNadir) >= 1.) then
         print *,'err[AngleUtil::nadirFromZenith]'
         call exit(1)
      end if
      nadirRad=asin(sinNadir)
      nadirFromZenith=nadirRad/deg2rad
      return

   end function nadirFromZenith

end module AngleUtil
