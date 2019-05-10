#include "fintrf.h"
C
#if 0
C     
C     tbScan2sdr.F
C     .F file needs to be preprocessed to generate .for equivalent
C     
#endif
C     
C     Matlab interface to putRad and putRad1chan.
C     2004/11/10:  John Galantowicz
C     Copyright AER, Inc. 2004
      
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on the DEC Alpha
C     64-bit platform
C
      implicit none

      include 'mxdims.incl'
c     MXLEV,MXCHAN,MXPARG

      integer plhs(*), prhs(*)
      mwpointer mxGetPr, mxCreateDoubleMatrix
      mwpointer nline_pr,nelement_pr,satalt_pr,surfalt_pr
      mwpointer lat_pr,lon_pr,eia_pr,eaa_pr,nchan_pr,nInner_pr
      mwpointer nscid_pr,instid_pr,sza_pr,saa_pr
      mwpointer modver_pr,nfilever_pr,FOVtype_pr,npos_pr
      mwpointer first_pr,last_pr,view_pr,speed_pr
      mwpointer tb_pr,surfpres_pr,pland_pr
      mwpointer frq_pr,pol_pr,aznom_pr,eianom_pr
      mwpointer nBeam_pr,mapBeam_pr
C-----------------------------------------------------------------------
C             

      integer nlhs, nrhs
      integer mxGetM, mxGetN,mxIsNumeric, mxIsString
      integer mxCreateString, mxGetString
      integer m, n, size, status,strlen, naznom,maznom
c     Max points = max(nchan)*max(npos)
      integer MXNPOS
      parameter(MXNPOS = 650)

      integer*4 nline,nelement(MXCHAN*MXNPOS),nchan,instid,ntb,nco
      integer*4 nscid,modver,nfilever,FOVtype,nInner
      integer*4 npos,view,first,last,pol(MXCHAN)
      integer*4 nBeam,mapBeam(MXCHAN),neia
      real*8    satalt,surfalt(MXCHAN*MXNPOS)
      real*8    lat(MXCHAN*MXNPOS),lon(MXCHAN*MXNPOS)
      real*8    eia(MXCHAN*MXNPOS),eaa(MXCHAN*MXNPOS)
      real*8    sza,saa,speed,tb(MXCHAN*MXNPOS),surfpres(MXCHAN*MXNPOS)
      real*8    pland(MXCHAN*MXNPOS)
      real*8    frq(MXCHAN),aznom(MXCHAN*MXNPOS),eianom(MXCHAN)

      character*8  xid
      character*200  fname

C     Check for proper number of arguments. 
      if(nrhs .ne. 33) then
         call mexErrMsgTxt('tbScan2sdr: Exact number inputs required.')
      elseif(nlhs .ne. 0) then
         call mexErrMsgTxt('tbScan2sdr: N0 outputs required.')
      endif

C----------------------------------------------------
      nline_pr = mxGetPr(prhs(1))
      nelement_pr = mxGetPr(prhs(2))
      satalt_pr = mxGetPr(prhs(3))
      surfalt_pr = mxGetPr(prhs(4))
      lat_pr = mxGetPr(prhs(5))
      lon_pr = mxGetPr(prhs(6))
      eia_pr = mxGetPr(prhs(7))
      eaa_pr = mxGetPr(prhs(8))
      saa_pr = mxGetPr(prhs(9))
      sza_pr = mxGetPr(prhs(10))
      nchan_pr = mxGetPr(prhs(11))
      nscid_pr = mxGetPr(prhs(12))
      instid_pr = mxGetPr(prhs(13))
      modver_pr = mxGetPr(prhs(14))
      nfilever_pr = mxGetPr(prhs(15))
      FOVtype_pr = mxGetPr(prhs(16))
      npos_pr = mxGetPr(prhs(17))
      view_pr = mxGetPr(prhs(18))
      first_pr = mxGetPr(prhs(19))
      last_pr = mxGetPr(prhs(20))
      speed_pr = mxGetPr(prhs(21))
      tb_pr = mxGetPr(prhs(22))
      pland_pr = mxGetPr(prhs(23))
      surfpres_pr = mxGetPr(prhs(24))
cc      prhs(25) = xid
cc      prhs(26) = fname
      frq_pr = mxGetPr(prhs(27))
      pol_pr = mxGetPr(prhs(28))
      aznom_pr = mxGetPr(prhs(29))
      eianom_pr = mxGetPr(prhs(30))
      nBeam_pr = mxGetPr(prhs(31))
      mapBeam_pr = mxGetPr(prhs(32))
      nInner_pr = mxGetPr(prhs(33))

c     SCALARS
      size=1
      call mxCopyPtrToInteger4(nline_pr,nline,size)
      call mxCopyPtrToInteger4(nchan_pr,nchan,size)
      call mxCopyPtrToInteger4(nscid_pr,nscid,size)
      call mxCopyPtrToInteger4(instid_pr,instid,size)
      call mxCopyPtrToInteger4(modver_pr,modver,size)
      call mxCopyPtrToInteger4(nfilever_pr,nfilever,size)
      call mxCopyPtrToInteger4(npos_pr,npos,size)
      call mxCopyPtrToInteger4(FOVtype_pr,FOVtype,size)
      call mxCopyPtrToInteger4(first_pr,first,size)
      call mxCopyPtrToInteger4(last_pr,last,size)
      call mxCopyPtrToInteger4(view_pr,view,size)
      call mxCopyPtrToInteger4(nBeam_pr,nBeam,size)
      call mxCopyPtrToInteger4(nInner_pr,nInner,size)
      
      call mxCopyPtrToReal8(satalt_pr,satalt,size)
      call mxCopyPtrToReal8(sza_pr,sza,size)
      call mxCopyPtrToReal8(saa_pr,saa,size)
      call mxCopyPtrToReal8(speed_pr,speed,size)

c     STRINGS
      status = mxGetString(prhs(25),xid,8)
      status = mxGetString(prhs(26),fname,200)

c     ARRAYS

c     Determine size of tb vector and convert
      m = mxGetM(prhs(22))
      n = mxGetN(prhs(22))
      ntb = m*n
      if (ntb.gt.MXCHAN*MXNPOS) then
         call mexErrMsgTxt('tbScan2sdr: ntb too large.')
      endif
c     Check that ntb is correct for FOVtype
      if (FOVtype.eq.1) then
         if (ntb.ne.npos) then
            call mexErrMsgTxt('tbScan2sdr: ntb ~= npos.')
         endif
      else
         if (ntb.ne.nchan) then
            call mexErrMsgTxt('tbScan2sdr: ntb ~= nchan.')
         endif            
      endif
      call mxCopyPtrToReal8(tb_pr,tb,ntb)

c     Determine size of data colocated with TBs
      m = mxGetM(prhs(2))
      n = mxGetN(prhs(2))
      nco = m*n
      if (nco.gt.MXCHAN*MXNPOS) then
         call mexErrMsgTxt('tbScan2sdr: nco too large.')
      endif
c     Check that nco is correct for FOVtype
      if (FOVtype.eq.1) then
         if (nco.ne.ntb) then
            call mexErrMsgTxt('tbScan2sdr: nco ~= ntb.')
         endif
      else
         if (nco.ne.1) then
            call mexErrMsgTxt('tbScan2sdr: nco ~= 1.')
         endif            
      endif
      call mxCopyPtrToInteger4(nelement_pr,nelement,nco)
      call mxCopyPtrToReal8(surfalt_pr,surfalt,nco)
      call mxCopyPtrToReal8(lat_pr,lat,nco)
      call mxCopyPtrToReal8(lon_pr,lon,nco)
      call mxCopyPtrToReal8(surfpres_pr,surfpres,nco)
      call mxCopyPtrToReal8(pland_pr,pland,nco)

c     Determin size of INPUT eia data (will be resized below)
      m = mxGetM(prhs(7))
      n = mxGetN(prhs(7))
      neia = m*n
      call mxCopyPtrToReal8(eia_pr,eia,neia)
      call mxCopyPtrToReal8(eaa_pr,eaa,neia)
      if (FOVtype.eq.1) then
         if (neia.ne.ntb) then
            call mexErrMsgTxt('tbScan2sdr: neia ~= ntb.')
         endif
      else
         if (neia.ne.nBeam) then
            call mexErrMsgTxt('tbScan2sdr: neia ~= nBeam.')
         endif            
      endif

c     Determine size of frq vector and convert all attribute vectors
      m = mxGetM(prhs(27))
      n = mxGetN(prhs(27))
      size = m*n
      if (size.gt.MXCHAN) then
         call mexErrMsgTxt('tbScan2sdr: frq too large.')
      endif
      if (size.ne.nchan) then
         call mexErrMsgTxt('tbScan2sdr: nchan mismatch.')
      endif
      call mxCopyPtrToReal8(frq_pr,frq,size)
      call mxCopyPtrToInteger4(pol_pr,pol,size)
      call mxCopyPtrToReal8(eianom_pr,eianom,size)
      call mxCopyPtrToInteger4(mapBeam_pr,mapBeam,size)

c     Determine size of aznom vector
      maznom = mxGetM(prhs(29))
      naznom = mxGetN(prhs(29))
      size = maznom*naznom
      if (size.gt.MXCHAN*MXNPOS) then
         call mexErrMsgTxt('tbScan2sdr: aznom too large.')
      endif
      if (naznom.ne.npos) then
         call mexErrMsgTxt('tbScan2sdr: naznom ~= npos.')
      endif
      call mxCopyPtrToReal8(aznom_pr,aznom,size)

C     Call the save data subroutine.

      call tbScan2sdr(xid,nline,nelement,satalt,surfalt,
     &     lat,lon,fname,eaa,eia,tb,sza,saa,nchan,frq,pol,
     &     nscid,instid,modver,nfilever,FOVtype,npos,view,
     &     speed,first,last,surfpres,pland,aznom,
     &     eianom,ntb,nco,nBeam,mapBeam,neia,naznom,maznom,nInner)

      return
      end

c*****************************************************************************
      subroutine tbScan2sdr(xid,nline,nelement,dsatalt,dsurfalt,
     $     dlat,dlon,fname,deaa,deia,dtb,dsza,dsaa,nchan,dfrq,pol,
     $     nscid,instid,modver,nfilever,FOVtype,npos,view,dspeed,
     $     first,last,dsurfpres,dpland,daznom,deianom,ntb,nco,nBeam,
     $     mapBeam,neia,naznom,maznom,nInner)

      use rad_io_module
      
      implicit none
      include 'mxdims.incl'
      
      integer mxnpos,mxpnts
      parameter (mxnpos = 1000, mxpnts = 1000)

      integer*4 norbit,nline,ntb,nco,nelement(nco)
      integer*4 nscid,instid,modver,nfilever
      integer*4 nchan,first,FOVtype,view,j,nInner
      integer*4 time(6),ncid,last,npos,i,pol(nchan)
      integer*4 length,qflg(mxchan),refaz(mxchan)
      integer*4 nposxBeam(mxchan),nBeam,mapBeam(nchan)
      integer*4 neia,naznom,maznom

      character*8   xid
      character*200 fname,auxfile
      
      real*8 dsatalt,dsurfalt(nco),dpland(nco),dspeed
      real*8 dlat(nco),dlon(nco),dsza,dsaa
      real*8 deaa(neia),deia(neia),dtb(ntb),dsurfpres(nco)
      real*8 dfrq(nchan),daznom(maznom*naznom),deianom(nchan)

      real*4 satalt,surfalt(mxpnts),pland(mxpnts),speed
      real*4 lat(mxpnts),lon(mxpnts),sza(mxpnts),saa(mxpnts)
      real*4 eaa(mxpnts),eia(mxpnts),tb(mxpnts),aznom(mxchan,mxpnts)
      real*4 frq(mxchan),eianom(mxchan),surfpres(mxpnts)
      
ccc   Map real*8 variables into real*4 to be ready to pass to putSDR         

      satalt = dsatalt
      surfalt(1:nco) = dsurfalt(1:nco)
      pland(1:nco) = dpland(1:nco)
      speed = dspeed
      lat(1:nco) = dlat(1:nco)
      lon(1:nco) = dlon(1:nco)
      sza(1:nco) = dsza
      saa(1:nco) = dsaa
      tb(1:ntb) = dtb(1:ntb)
      surfpres(1:nco) = dsurfpres(1:nco)
      frq(1:nchan) = dfrq(1:nchan)
      eianom(1:nchan) = deianom(1:nchan)
c     aznom: naznom=npos and maznom=1 for now
      do i=1,naznom
         do j=1,maznom
            aznom(i,j) = daznom((i-1)*maznom+j)
         enddo
      enddo
      
      norbit = 1
      time(1) = 2009
      time(2) = 9
      time(3) = 9
      time(4) = 0
      time(5) = 0
      time(6) = 0
         
c     create netcdf scene profiles for RT model to read
      
      if (first .eq. 1) then 
         print*,'tbScan2sdr: first time ',first,last
         if (FOVtype .eq. 1) then
c     NOTE: nchan = 1 for this case ALWAYS; nInner not allowed
            nposxBeam(1) = ntb
            call openRad(ncid=ncid,fname=fname,satalt=satalt,
     $           pol=pol,frq=frq(1:nchan),nchan=nchan,nscid=nscid,
     $           instid=instid,modver=modver,nfilever=nfilever,
     $           eianom=eianom(1:nchan),aznom=aznom(1:naznom,1:maznom),
     $           speed=speed,view=view,nBeam=nBeam,
     $           mapBeam=mapBeam(1:nchan),
     $           FOVtype=FOVtype,status='new',
     $           npos=nposxBeam(1:nBeam))
         else
c     Co-located; nInner allowed but npos must be absent or 1's
            nposxBeam(1:nBeam) = 1
            call openRad(ncid=ncid,fname=fname,satalt=satalt,
     $           pol=pol,frq=frq(1:nchan),nchan=nchan,nscid=nscid,
     $           instid=instid,modver=modver,nfilever=nfilever,
     $           eianom=eianom(1:nchan),aznom=aznom(1:naznom,1:maznom),
     $           speed=speed,view=view,nBeam=nBeam,
     $           mapBeam=mapBeam(1:nchan),
     $           FOVtype=FOVtype,status='new',
     $           npos=nposxBeam(1:nBeam),nInner=nInner)
         endif
      endif
      if (FOVtype .eq. 1) then  
c     Single channel EFOV data on scan
         eaa(1:ntb) = deaa(1:ntb)
         eia(1:ntb) = deia(1:ntb)
         call putRad1chan(ncid=ncid,xid=xid,norbit=norbit,
     $        nline=nline,lat=lat(1:ntb),lon=lon(1:ntb),
     $        time=time,eaa=eaa(1:ntb),eia=eia(1:ntb),
     $        tb=tb(1:ntb),surfalt=surfalt(1:ntb))
      else 
c     Collocated multi-channel CFOV data on point
c     sza=sza(1),saa=saa(1) removed until they can be simulated
         eaa(1:nBeam) = deaa(1:nBeam)
         eia(1:nBeam) = deia(1:nBeam)
         call putRad(ncid=ncid,xid=xid,norbit=norbit,nline=nline,
     $        nelement=nelement(1),lat=lat(1),lon=lon(1),time=time,
     $        eaa=eaa(1:nBeam),eia=eia(1:nBeam),tb=tb(1:ntb),
     $        surfalt=surfalt(1))
      endif
      if (last .eq. 1) then 
         print*,'tbScan2sdr: closing',first,last
         call closeRad(ncid=ncid)
      endif
      
      return
      end



c****************************************************************      
      FUNCTION LENGTH(STRING,MAX)
      INTEGER MAX
      CHARACTER*1 STRING(MAX)
      DO L=MAX,1,-1
         IF(STRING(L).NE.' ') GO TO 10
      END DO
      LENGTH=0
      RETURN
   10 LENGTH=L
      RETURN
      END










