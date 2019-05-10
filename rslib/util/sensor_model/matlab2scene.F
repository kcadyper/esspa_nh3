#include "fintrf.h"
C
#if 0
C     
C     matlab2scene.F
C     .F file needs to be preprocessed to generate .for equivalent
C     
#endif
C     
C     Matlab interface to putScene.
C     2004/11/10:  John Galantowicz
C     Copyright AER, Inc. 2004
      
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on the DEC Alpha
C     64-bit platform
C

      use stateindexmodule
      use scene_io_module
      
      implicit none

      include 'mxdims.incl'
c     MXLEV,MXCHAN,MXPARG

      integer plhs(*), prhs(*)
      mwpointer mxGetPr, mxCreateDoubleMatrix
      mwpointer pref_pr,freq_pr,nlevel_pr,pland_pr
      mwpointer nchan_pr,ngset_pr,xvars_pr,emis_pr,nInner_pr
      mwpointer eia_pr,surfalt_pr,first_pr,last_pr
      mwpointer lat_pr,lon_pr,pol_pr,eaa_pr
C-----------------------------------------------------------------------
C             

      integer nlhs, nrhs
      integer mxGetM, mxGetN,mxIsNumeric, mxIsString
      integer mxCreateString, mxGetString
      integer m, n, size, status,strlen
      integer ncid,nparmwg,imH2O

      integer*4 first,last,pol(MXCHAN),nlevel,nchan,landtype,nInner
      integer*4 ngset(7),MolID(maxMol),time(6)
      real*8    surfalt,lat,lon,eia(1),eaa(1)
      real*8    pland,pref(MXLEV),xvars(MXPARG),emis(MXCHAN)
      real*8    freq(MXCHAN)

      character*200  fname,casename
      character*5    xid

      type (stateindex_t) :: ig,ng

C     Check for proper number of arguments. 
      if(nrhs .ne. 19) then
         call mexErrMsgTxt('matlab2scene: 18 inputs required. ')
      elseif(nlhs .ne. 0) then
         call mexErrMsgTxt('N0 outputs required.')
      endif

C----------------------------------------------------
c     1 is str
      pref_pr = mxGetPr(prhs(2))
      pol_pr = mxGetPr(prhs(3))
      freq_pr = mxGetPr(prhs(4))
      nlevel_pr = mxGetPr(prhs(5))
      nchan_pr = mxGetPr(prhs(6))
c     7 is str
      ngset_pr = mxGetPr(prhs(8))
      xvars_pr = mxGetPr(prhs(9))
      pland_pr = mxGetPr(prhs(10))
      eia_pr = mxGetPr(prhs(11))
      surfalt_pr = mxGetPr(prhs(12))
      first_pr = mxGetPr(prhs(13))
      last_pr = mxGetPr(prhs(14))
      lat_pr = mxGetPr(prhs(15))
      lon_pr = mxGetPr(prhs(16))
      lon_pr = mxGetPr(prhs(16))
      eaa_pr = mxGetPr(prhs(17))
      emis_pr = mxGetPr(prhs(18))
      nInner_pr = mxGetPr(prhs(19))

      size=1
      call mxCopyPtrToInteger4(first_pr,first,size)
      call mxCopyPtrToInteger4(last_pr,last,size)
      call mxCopyPtrToInteger4(nlevel_pr,nlevel,size)
      call mxCopyPtrToInteger4(nchan_pr,nchan,size)
      call mxCopyPtrToInteger4(nInner_pr,nInner,size)

      call mxCopyPtrToReal8(surfalt_pr,surfalt,size)
      call mxCopyPtrToReal8(lat_pr,lat,size)
      call mxCopyPtrToReal8(lon_pr,lon,size)
      call mxCopyPtrToReal8(eia_pr,eia(1),size)
      call mxCopyPtrToReal8(pland_pr,pland,size)
      call mxCopyPtrToReal8(eaa_pr,eaa(1),size)

      status = mxGetString(prhs(1),fname,200)
      status = mxGetString(prhs(7),casename,200)

      m = mxGetM(prhs(2))
      n = mxGetN(prhs(2))
      size = m*n
      if (size.gt.MXLEV) then
         call mexErrMsgTxt('matlab2scene: pref too large.')
      endif
      call mxCopyPtrToReal8(pref_pr,pref,size)

      m = mxGetM(prhs(4))
      n = mxGetN(prhs(4))
      size = m*n
      if (size.gt.MXCHAN) then
         call mexErrMsgTxt('matlab2scene: freq too large.')
      endif
      call mxCopyPtrToReal8(freq_pr,freq,size)
      call mxCopyPtrToInteger4(pol_pr,pol,size)

      m = mxGetM(prhs(8))
      n = mxGetN(prhs(8))
      size = m*n
      if (size.gt.7) then
         call mexErrMsgTxt('matlab2scene: ngset too large.')
      endif
      call mxCopyPtrToInteger4(ngset_pr,ngset,size)

      m = mxGetM(prhs(9))
      n = mxGetN(prhs(9))
      size = m*n
      if (size.gt.MXPARG) then
         call mexErrMsgTxt('matlab2scene: xvars too large.')
      endif
      call mxCopyPtrToReal8(xvars_pr,xvars,size)
      nparmwg = size

      m = mxGetM(prhs(18))
      n = mxGetN(prhs(18))
      size = m*n
      if (size.gt.MXCHAN) then
         call mexErrMsgTxt('matlab2scene: emis too large.')
      endif
      call mxCopyPtrToReal8(emis_pr,emis,size)

c     Derive/fake parameters needed for netcdf scene file
      time(1) = 2009
      time(2) = 9
      time(3) = 9
      time(4) = 0
      time(5) = 0
      time(6) = 0
      
      landtype = 1
      xid = 'InSim'

c     Set state vector indexes
      ng%temp = ngset(1)
      ng%tskin = ngset(2)
      ng%psfc = ngset(3)
      imH2O = 1
      ng%mol(imH2O) = ngset(4)
      ng%cldliq = ngset(5)
      ng%cldice = ngset(6)
      ng%wind = ngset(7)

      ig = genindices(ng)
      MolID(1) = 1

c     Open (if first), write, and close (if last) profile
      
      if (first.eq.1) then 
         call openScene(ncid=ncid,file=fname,
     $        pressure=real(pref(1:nlevel)),polarity=pol(1:nchan),
     $        freq=real(freq(1:nchan)),casename=casename,
     $        nLevel=nlevel,nchmw=nchan,ig=ig,ng=ng,status='new',
     $        MolID=MolID,nInner=nInner)
      endif

      call putScene(ncid=ncid,lat=real(lat), 
     $     lon=real(lon),x=real(xvars(1:nparmwg)),pland=real(pland), 
     $     eia=real(eia(1:1)),time=time,eaa=real(eaa(1:1)),
     $     surfalt=real(surfalt),EmMw=real(emis(1:nchan)),
     $     landtype=landtype,xid=xid)
      
      if (last.eq.1) then 
         call closeScene(ncid)
      endif

      return
      end
