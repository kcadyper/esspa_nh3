c     IDL interface to putScene with pressure level conversions.
C     2004/11/10:  John Galantowicz
C     Copyright AER, Inc. 2004

      subroutine idl2scene(argc,argv)
c     Was slant2netcdf.f and slant2netcdf_core.f

      Implicit None
      Integer*4 argc,argv(*),j
       	
      j = LOC(argc)

      call idl2sceneCore(%val(argv(1)),
     +     %val(argv(2)),
     +     %val(argv(3)),
     +     %val(argv(4)),
     +     %val(argv(5)),
     +     %val(argv(6)),
     +     %val(argv(7)),
     +     %val(argv(8)),
     +     %val(argv(9)),
     +     %val(argv(10)),
     +     %val(argv(11)),
     +     %val(argv(12)),
     +     %val(argv(13)),
     +     %val(argv(14)),
     +     %val(argv(15)),
     +     %val(argv(16)),
     +     %val(argv(17)),
     +     %val(argv(18)),
     +     %val(argv(19)),
     +     %val(argv(20)),
     +     %val(argv(21)),
     +     %val(argv(22)),
     +     %val(argv(23)),
     +     %val(argv(24)),
     +     %val(argv(25)),
     +     %val(argv(26)),
     +     %val(argv(27)),
     +     %val(argv(28)),
     +     %val(argv(29)),
     +     %val(argv(30)),
     +     %val(argv(31)),
     +     %val(argv(32)),
     +     %val(argv(33)),
     +     %val(argv(34)),
     +     %val(argv(35)),
     +     %val(argv(36)),
     +     %val(argv(37)),
     +     %val(argv(38)))
      
      return
      end

c*****************************************************************************
      subroutine idl2sceneCore(icasename,ixid,lat,lon,time,pland,
     +     surfType,nLevelnwp,Tnwp,H2Onwp,hydrnwp,tskin,surfpres,uwind,
     +     vwind,ifname,eaaPath,eiaPath,nchan,frq,pol,
     +     first,last,pnwp,iBeam,ncasename,nxid,nfname,ncid,
     +     surfalt,t2m,q2m,nBeam,nomEIAxiBeam,nInner,nLevelH2O,
     +     nlevHy,nhydr)
      
      use ncdf_module
      use stateindexmodule
      use scene_io_module
      use constants

      implicit none
      
      include 'mxdims.incl'

      integer nLevelnwp,surfType,nchan,first,nparmwg,nLevelH2O
      integer nlevHy,nhydr
      integer i
      integer time(6),ncid,last,nBeam
      integer iBeam(nchan),ncasename,icasename(ncasename)
      integer nxid,ixid(nxid),nfname,ifname(nfname)
      integer MolID(maxMol),imH2O,nInner
      integer pol(nchan)
      
      character*20  casename
      character*12   xid
      character*200 fname
      
      real lat,lon,pland,surfalt,nomEIAxiBeam(nBeam)
      real t2m,q2m
      real Tnwp(nLevelnwp),H2Onwp(nLevelnwp),hydrnwp(nlevHy,nhydr)
      real eaaPath,eiaPath
      real eaaLocal(mxchan),eiaLocal(mxchan)
      real frq(nchan),surfpres,uwind
      real vwind,emmw(nchan),pnwp(nLevelnwp),ctop,intclw
      real cthk,itop,intice
      real tskin
      real xvars(MXPARG)

      type (stateindex_t) :: ig,ng

ccccc JG 12/2/03:  idl_2_fort string input method failed.  Replaced with 
ccccc integerized string arguments passed from IDL.

      if (ncasename.gt.20.or.nxid.gt.12.or.nfname.gt.200) then
         print *,'err[idl2scene::idl2sceneCore]: ',
     $        ' character array exceeds max dimension.'
         call errorHalt(1)
      endif
      do i=1,ncasename
         casename(i:i) = char(icasename(i))
      enddo
      do i=1,nxid
         xid(i:i) = char(ixid(i))
      enddo
      fname = "\0"
      do i=1,nfname
         fname(i:i) = char(ifname(i))
      enddo

      MolID(1) = 1
      imH2O = 1

c     Allocate xvars array passed to putscene.
c     JG: IDL user notes for getting dynamic memory indicate a problem
c     if IDL-written routines aren't used.  So we get rid of dynamic
c     allocation here and just declare a generously large array above.
c      allocate(xvars(nparmwg))

C     Find column integrated cloud liquid water in kg/m2

      if (nhydr >= 1) then
         call prof2CLW(pnwp,hydrnwp(1:nlevHy,1),nlevHy,
     $   surfpres,intclw,ctop,cthk)
      else
         intclw=0.
         ctop=0.
         cthk=0.
      endif

c     Fill the xvars array to pass to putscene.

      ng=initLengths(Ntemp=nLevelnwp,Ntskin=1,Npsfc=1,Ncldliq=3,Nwind=2)
      ng%mol(imH2O) = nLevelH2O

c     Automatic conversion of ng structure for number of elements per
c     parameter to starting index structure  ig.  Try "call printIndex(ig,ng)"
c     for details.
      nparmwg = getVectorLength(ng)
      ig = genindices(ng)

      xvars(IG%temp:IG%temp+NG%temp-1) = Tnwp ! temperature
      xvars(IG%mol(imH2O):IG%mol(imH2O)+NG%mol(imH2O)-1) = 
     $     H2Onwp(1:nLevelH2O)                ! water vapor
      xvars(IG%tskin) = tskin
      xvars(IG%psfc) = surfpres
      xvars(IG%cldliq) = ctop
      xvars(IG%cldliq+1) = cthk
      xvars(IG%cldliq+2) = intclw
c      xvars(IG%cldice) = 0.0 ! itop
c      xvars(IG%cldice+1) = 0.0 !intice
      xvars(IG%wind) = uwind
      xvars(IG%wind+1) = vwind

c     Emissivities are added in a subsequent process.
c     Fill in "missing" values needed as a placeholder in the file structure.

      emmw(:) = MISSING_REAL

c     Negative nomEIAxiBeam is flag to populate eia and eaa local with 
c     EIA and EAA of key channel obtained from orbsim.  These eia and eaa 
c     account for pointing error, via orbsim, but neglect alignment difference
c     from key channel.  Alternative use of nomEIA accounts for alignment 
c     difference from key channel, but neglects pointing error.
c     Could change to run orbsim for each channel for this purpose.
      if (any(nomEIAxiBeam(1:nBeam) < -9.)) then
         eiaLocal(1:nBeam) = eiaPath
      else
         eiaLocal(1:nBeam) = nomEIAxiBeam(1:nBeam)
      endif
      eaaLocal(:) = eaaPath

c     
c     Create netcdf scene file. See scene_io_module.f for 
c     openScene and putScene and retr.f and sim.f for example usage
c     
      if (first > 0) then 
         print*,'Before openScene: ',first,last
         print*,trim(fname)
         call openScene(ncid=ncid,file=fname(1:nfname),
     $        pressure=pnwp,polarity=pol(1:nchan),
     $        freq=frq(1:nchan),casename=casename,MolID=MolID,
     $        nLevel=nLevelnwp,nchmw=nchan,ig=ig,ng=ng,status='new',
     $        nInner=nInner,nBeam=nBeam,mapBeam=iBeam)
         first = 0
         print*,'After openScene: ',first,last
      endif
c     xid is meaningless but are needed by sim.f90 currently

      call putScene(ncid=ncid,xid=xid,lat=lat,lon=lon,
     $     x=xvars(1:nparmwg),pland=pland, 
     $     eia=eiaLocal(1:nBeam),emmw=emmw(1:nchan),landtype=surfType,
     $     time=time,surfalt=surfalt,eaa=eaaLocal(1:nBeam))
      if (t2m > 0.) then  ! 2m data not missing
         call writeNcdfData(ncid,t2m,varname='T2m',
     $     varlenName='nProfiles',varLongName='2-m Temperature',
     $     varUnit='K',status='append')
         call writeNcdfData(ncid,q2m,varname='q2m',
     $     varlenName='nProfiles',varLongName='2-m Mixing Ratio',
     $     varUnit='g/g',status='append')
      endif
      if (nhydr > 1) then
         call writeNcdfData(ncid,hydrnwp(1:nlevHy,1:nhydr),
     $     varname='hydr',
     $     varlenName=(/'nlevHydr ','nHydr    ','nProfiles'/),
     $     varLongName='cloud liquid, cloud ice, rain, snow, graupel',
     $     varUnit='mixing ratio (g/g)',status='append')
      endif

      if (last > 0) then 
         call closeScene(ncid)
         print*,'Closing: ',first,last
      endif

      return
      end
