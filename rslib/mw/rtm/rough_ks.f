      subroutine rough_ks(f,angle,temp,sal,ws,vert,horiz)

c----------------------------------------------------------------------
c  rough calculates the surface reflectivities for a rough ocean
c  surface.  rough was written by tom wilheit, nasa gsfc but we have
c  substituted klein and swift's calculations for dielectric constant
c  of ocean water
c     input:
c     f      - (real) frequency in ghz.
c     angle  - (real) incidence angle in degrees.
c     temp   - (real) sea surface temperature iin degrees kelvin.
c     sal    - (real) salinity in parts per thousand.
c     ws     - (real) wind speed in m/s.
c
c     output:
c     vert   - (real) emissivity for vertical polarization.
c     horiz  - (real) emissivity for horizontal polarization.
c----------------------------------------------------------------------

c      implicit double precision (a-h,o-z)

      real f, angle, temp, sal, ws, vert, horiz

      var = 0.003 + 4.80*ws/1000.
      if (f .lt. 40.) var = var*(0.06852+0.03392*f)
      t = temp
      vert = 0.
      horiz = 0.
      ynorm = 0.
      sig = sqrt(var)
      step = sig/10.

      angr = angle/57.2957795
      sa = sin(angr)
      ca = cos(angr)

c index of refraction calculation    lane and saxton dielectric
c and international critical tables conductivity data

c set ppt to 35.0, it was never initialized (cel, 5/26/98)
c >>
      ppt=35.0
c <<

      s = ppt/58.4
      xlam = 51. - 4.835*s + 1.728*(t-273.) - 0.2037*s*(t-273.)
      if (s .gt. 1.) xlam = xlam - 0.35/s
      es = 190. - 81*s + 38*s*s - (03.75-2.*s+s*s)*t/10.
      t1 = 0.00199 * exp(2140/t)/t
      t2 = 0.00243 * exp(2060/t)/t - t1
      t3 = 0.00324 * exp(1968/t)/t - t1
      tau = t1 + (4.*t2 - t3)*s + (2.*t3 - 4.*t2)*s*s
c  tau in [nanoseconds], f in [ghz], t in deg [k]

c  compute dielectric constant of sea water.
      call nwat_ks(f,temp,sal,rpe,xipe)

      xmode = sqrt(rpe**2 + xipe**2)
      xmodn = sqrt(xmode)
      argn = 0.5*atan(xipe/rpe)
      rpn = xmodn*cos(argn)
      xipn = xmodn*sin(argn)

      do 100 j = 1,61
         do 100 kk = 1,31
            k = kk-1
            zx = step*(j-31)
            zy = step*k
            weight = exp(-((j-31)**2+k*k)/100.)
            if (kk .eq. 1) weight = weight/2.

c  facet hidden from view
            rdn = ca + zx*sa
            if (rdn .lt. 0.) go to 100

c123456789012345678901234567890123456789012345678901234567890123456789
            stp = (sa*sa + zx*zx + zy*zy - (zx*sa)**2 -
     +             2.*zx*sa*ca) / (1 + zx*zx+zy*zy)
C bug fix jdh 17 May 2001
c            if(stp.lt.0.)print*,'stp < 0',stp
            stp=max(stp,1.0e-4)
            si = sqrt(stp)
            if (si .lt. 0.001) si = 0.001
            em = zy * zy / (1. + zx*zx + zy*zy)
            em = em / (si*si)
            if (si .lt. 0.00101) em = .5
            ci = sqrt(1. - si*si)
            weight = weight*ci
            ynorm = ynorm + weight

c  projection of facet in view direction
            an2 = rpn**2 + xipn**2
            srr = rpn*si / an2
            sri = -xipn*si / an2
            xmod = ((1. - srr**2 + sri**2)**2 + 
     +              (2.*srr*sri)**2)**0.25
            arg = -0.5 * atan(2.*srr*sri / (1. - srr**2 + sri**2))
            crr = xmod * cos(arg)
            cri = xmod * sin(arg)

c.....calculate the reflected emissivity with fresnel equation
            rh = ((si*crr - ci*srr)**2 + 
     +            (si*cri - ci*sri)**2) / ((si*crr + ci*srr)**2 +
     +            (si*cri + sri*ci)**2)
            rv = rh * ((ci*crr - si*srr)**2 + 
     +                 (ci*cri - si*sri)**2) / ((ci*crr + si*srr)**2
     +                +(ci*cri + si*sri)**2)

            horiz = horiz + weight*((1.-em)*(1.-rh) + em*(1.-rv))
            vert = vert + weight*((1.-em)*(1.-rv) + em*(1.-rh))
  100 continue

      horiz = horiz/ynorm
      vert = vert/ynorm

c.....adds foam effect
c     no increase in foam effect above 37 ghz (ie use old foam model)
      if (f .lt. 38) then
         if (ws .le. 7.0 ) then
c23456789112345678921234567893123456789412345678951234567896123456789712
            foam  = 0.0012 * ws
            foamv = foam + foam*(-0.7783 + 0.1643*f - 0.005604*f*f
     &           + 0.0001066*f*f*f)
            foamh = foam + foam*(-0.36025 + 0.12506*f - 0.0051709*f*f
     &           + 0.000092046*f*f*f)

         else if ((ws .gt. 7.0 ) .and. (ws .lt. 17.0 )) then
            foam  = 0.000195*ws*ws - 0.00153*ws + 0.009555
            foamv = foam + foam*( -1.103 + 0.1866*f - 0.005541*f*f
     &           + 0.0001223*f*f*f +
     &           (0.03423 - 0.0006405*f - 0.000141519*f*f)*ws)
            foamh = foam + foam*(-0.1991 + 0.109*f - 0.004426*f*f
     &           + 0.00009429*f*f*f +
     &           (-0.02475 + 0.002654*f - 0.00012682*f*f)*ws)
         else
            foam  = 0.0051 * ws - 0.04587
            foamv = foam + foam*(-0.3146 + 0.1325*f - 0.005697*f*f
     &           + 0.00008418*f*f*f)
            foamh = foam + foam*(-0.58956 + 0.1478*f - 0.0062256*f*f
     &           + 0.00008869*f*f*f)
         endif

         horiz = foamh + (1. - foamh) * horiz
         vert  = foamv + (1. - foamv) * vert

      else
         if (ws .gt. 7.) then
            foam = .0060 * (ws-7.) * (1.-exp(-f/7.5))
            horiz = foam + (1.-foam)*horiz
            vert  = foam + (1.-foam)*vert
         else
            foam = 0
         endif
      endif
      return
      end

c****************************************************************************

      subroutine nwat_ks(freq,temp,sal,epsr,epsi)

c----------------------------------------------------------------------
c  nwat calculates the dielectric constant of water.
c     input:
c     freq   - frequency in ghz
c     t      - temperature in degrees kelvin.
c     sal    - salinity in ppt.
c
c     output:
c     epsr   - dielectric constant, real part
c     epsi   - dielectric constant, imaginary part
c----------------------------------------------------------------------

c      implicit double precision (a-h,o-z)

      real freq, temp, sal

      pi = acos(-1.)
      f = freq*1.e09
      t = temp - 273.15
      s = sal

c Calculate ionic conductivity (t in c and salinity in ppt)
      delta = 25. - t
      beta = 2.033e-02 + 1.266e-04*delta + 2.464e-06*delta**2
     -    -s*(1.849e-05 - 2.551e-07*delta + 2.551e-08*delta**2)
      sigma = s*(1.82521e-01 - 1.46192e-03*s + 2.09324e-5*s**2
     -     -1.28205e-07*s*s*s)
      sigma = sigma * exp(-delta*beta)

c Calculate static dielectric constant
      epss = 87.134 - 1.949e-01*t - 1.276e-02*t**2 + 2.491e-04*t*t*t
      ast = 1.0 + 1.613e-05*s*t - 3.656e-03*s + 3.210e-05*s**2
     -    - 4.232e-07*s*s*s
      epss = epss * ast

c Calculate relaxation time
      tau = 1.768e-11 - 6.086e-13*t + 1.104e-14*t**2 - 8.111e-17*t*t*t
      bst = 1.0 + 2.282e-05*s*t - 7.638e-04*s - 7.760e-06*s**2
     -    + 1.105e-08*s*s*s
      tau = tau*bst
      eps0 = 8.854e-12
      epsinf = 4.9
      alpha = 0.0
      omega = 2*pi*f
      xa = omega*tau
      xb = sigma / (omega*eps0)
      ccmod = 1. / (1. + xa**2)
      epsr = epsinf + (epss-epsinf)*ccmod
      epsi = -(epss-epsinf) * ccmod*xa - xb
      return
      end
