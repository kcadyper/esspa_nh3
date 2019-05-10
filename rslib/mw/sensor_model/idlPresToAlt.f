c $Id$

c     IDL interface to presToAlt.
c     2004/11/10:  John Galantowicz
c     Copyright AER, Inc. 2004

      subroutine idlPresToAlt(argc,argv)
c     IDL call must be lower case ('idlprestoalt_') even though upper is
c     used here

      Implicit None
      Integer*4 argc,argv(*),j
       	
      j = LOC(argc)

      call idlPresToAlt_core(%val(argv(1)),
     +     %val(argv(2)),
     +     %val(argv(3)),
     +     %val(argv(4)),
     +     %val(argv(5)),
     +     %val(argv(6)),
     +     %val(argv(7)),
     +     %val(argv(8)),
     +     %val(argv(9)))
      
      return
      end

c*****************************************************************************
      subroutine idlPresToAlt_core(lat,lon,profT,profH,pSurface,surfalt,
     $     nlev,pref,zref)
      use MetFunctions, only : presToAlt
      use LvlInterp
c     Calculates zref(pref) pressure to altitude
      
      use MapInvert

      implicit none
      
      integer nlev,nref
      real lat,lon,profT(nlev),profH(nlev),pSurface,surfalt,
     $     tsfc,hsfc,zref(nlev),pref(nlev),dxlvldp

c     nref = largest index (nref) is the 1st level below surface

      call lvl_int(profT,pref,nlev,pSurface,tsfc,dxlvldp,nref)
      call lvl_int(profH,pref,nlev,pSurface,hsfc,dxlvldp,nref)

c     Inputs and outputs are in km, water is g/g      
      call presToAlt(lat,lon,profT,profH,
     $     pref,nref,surfalt/1000.0,
     $     pSurface,tsfc,hsfc,zref)

      return
      end
