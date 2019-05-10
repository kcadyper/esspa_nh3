c &"$Id$"/

      subroutine emis_ocean_mw(scanang,zenang,frqc,ipol,nchan,temp,ws,
     &    em)
c**********************************************************************
c* Function name: emis_ocean_mw
c* Description: (see edmsp)
c* Usage: call emis_ocean_mw(scanang,zenang,frqc,ipol,nchan,temp,ws,em)
c* Input Args:
c*  var     i*n      decription
c*  ---     ---      ----------
c* Output Args:
c*  var     i*n      decription
c*  ---     ---      ----------
c*
c**********************************************************************
      integer ipol(*)
      real*4 frqc(*),ev,eh,em(*)
      data sal/34.0/
      save sal

      alpha=scanang/180.*acos(-1.)
      xxa=cos(alpha)*cos(alpha)
      xxb=sin(alpha)*sin(alpha)
      do i=1,nchan

         call rough_ks(frqc(i),zenang,temp,sal,ws,ev,eh)
c         call find_emiss_tot(frqc(i),zenang,temp,ws,sal,ev,eh)
c     
         select case (ipol(i))
            case (0)     ! quasi vertical
             em(i) =xxa*ev+xxb*eh
            case (1)     ! quasi horizontal
             em(i) =xxa*eh+xxb*ev
            case (2:5)   ! +45deg, -45deg, LCirc, RCirc
             em(i) =0.5*eh+0.5*ev
            case default
             print *,'ERROR in emis_ocean_mw: Unsupported polarization;'
     &               ,' i,ipol(i):',i,ipol(i)
         end select
         
      end do
      
      return
      end
      
