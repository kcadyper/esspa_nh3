#include "fintrf.h"
C
#if 0
C     
C     getlatlon.F
C     .F file needs to be preprocessed to generate .for equivalent
C     
#endif
C     
C      
C     Matlab interface to viewGeom.
C     2004/11/10:  John Galantowicz
C     Copyright AER, Inc. 2004
      
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on the DEC Alpha
C     64-bit platform
C
      
      use GeomSat

      implicit none

      integer plhs(*), prhs(*)
      mwpointer mxGetPr, mxCreateDoubleMatrix
      mwpointer antCone_pr, antScan_pr, orbtime_pr, inlon2_pr
      mwpointer calim_pr, alt_pr, Incl_pr, inlon1_pr
      mwpointer ecross_pr, inlat1_pr, inlat2_pr
      mwpointer flat_pr,flon_pr,slat_pr,slon_pr
      mwpointer zenang_pr,azang_pr,a2sRPY_pr,s2scRPY_pr
C-----------------------------------------------------------------------
C             

      integer nlhs, nrhs
      integer mxGetM, mxGetN, mxIsNumeric, mxGetString
      integer m, n, size, status, i
      real*8  antCone,antScan,orbtime
      real*4  sensorCone,sensorScan,satCone,satScan
      real*8  alt,Incl,calim,inlat1,inlat2
      real*8  ecross,inlon1,inlon2
      real*8  a2sRPY(3),s2scRPY(3)
      real*4  flat,flon,slat,slon,zenang,azang

C     Check for proper number of arguments. 
      if(nrhs .ne. 13) then
         call mexErrMsgTxt('getlatlon.F: 13 inputs required.')
      elseif(nlhs .ne. 6) then
         call mexErrMsgTxt('getlatlon.F: 6 outputs required.')
      endif

C----------------------------------------------------

C     Check the size of the input arrays.
      do i=1,11
         m = mxGetM(prhs(i))
         n = mxGetN(prhs(i))
         size = m*n
         if (size.ne.1) then
            call mexErrMsgTxt('getlatlon.F: Inputs 1-11 size must = 1')
         endif
      enddo
      do i=12,13
         m = mxGetM(prhs(i))
         n = mxGetN(prhs(i))
         size = m*n
         if (size.ne.3) then
            call mexErrMsgTxt('getlatlon.F: Inputs 12-13 size must = 3')
         endif
      enddo
C----------------------------------------------------

C     Check to insure the input is a number.

      if(mxIsNumeric(prhs(1)) .eq. 0) then
         call mexErrMsgTxt('Input must be a number.')
      endif

      if(mxIsNumeric(prhs(2)) .eq. 0) then
         call mexErrMsgTxt('Input must be a number.')
      endif

      if(mxIsNumeric(prhs(3)) .eq. 0) then
         call mexErrMsgTxt('Input must be a number.')
      endif
C----------------------------------------------------

C     Create matrix for the return argument.
      do i=1,nlhs
         plhs(i) = mxCreateDoubleMatrix(1,1,0)
      enddo

c     Input pointers
      antCone_pr = mxGetPr(prhs(1))
      antScan_pr = mxGetPr(prhs(2))
      orbtime_pr = mxGetPr(prhs(3))
      alt_pr = mxGetPr(prhs(4))
      calim_pr = mxGetPr(prhs(5))
      Incl_pr = mxGetPr(prhs(6))
      ecross_pr = mxGetPr(prhs(7))
      inlat1_pr = mxGetPr(prhs(8))
      inlat2_pr = mxGetPr(prhs(9))
      inlon1_pr = mxGetPr(prhs(10))
      inlon2_pr = mxGetPr(prhs(11))
      a2sRPY_pr = mxGetPr(prhs(12))
      s2scRPY_pr = mxGetPr(prhs(13))

c     Output pointers
      flat_pr = mxGetPr(plhs(1))
      flon_pr = mxGetPr(plhs(2))
      slat_pr = mxGetPr(plhs(3))
      slon_pr = mxGetPr(plhs(4))
      zenang_pr = mxGetPr(plhs(5))
      azang_pr = mxGetPr(plhs(6))

c     Convert inputs
      size = 1
      call mxCopyPtrToReal8(antCone_pr,antCone,size)
      call mxCopyPtrToReal8(antScan_pr,antScan,size)
      call mxCopyPtrToReal8(orbtime_pr,orbtime,size)
      call mxCopyPtrToReal8(alt_pr,alt,size)
      call mxCopyPtrToReal8(calim_pr,calim,size)
      call mxCopyPtrToReal8(Incl_pr,Incl,size)
      call mxCopyPtrToReal8(ecross_pr,ecross,size)
      call mxCopyPtrToReal8(inlat1_pr,inlat1,size)
      call mxCopyPtrToReal8(inlat2_pr,inlat2,size)
      call mxCopyPtrToReal8(inlon1_pr,inlon1,size)
      call mxCopyPtrToReal8(inlon2_pr,inlon2,size)
      size = 3
      call mxCopyPtrToReal8(a2sRPY_pr,a2sRPY,size)
      call mxCopyPtrToReal8(s2scRPY_pr,s2scRPY,size)

c     Rotate nadir (antCone) and azimuth (antScan) from antenna to sensor
c     to spacecraft reference frame.

      call movePointing(real(antCone),real(antScan),real(a2sRPY(1)),
     $     real(a2sRPY(2)),real(a2sRPY(3)),sensorCone,sensorScan)
      call movePointing(sensorCone,sensorScan,real(s2scRPY(1)),
     $     real(s2scRPY(2)),real(s2scRPY(3)),satCone,satScan)
      
      call viewGeom(real(ecross),real(orbtime),satScan,satCone,slat, 
     $     slon,zenang,azang,real(alt),real(calim),real(Incl),flat,
     $     flon)

C     Load the data into y_pr, which is the output to MATLAB
      size = 1
      call mxCopyReal8ToPtr(dble(flat),flat_pr,size)     
      call mxCopyReal8ToPtr(dble(flon),flon_pr,size)     
      call mxCopyReal8ToPtr(dble(slat),slat_pr,size)     
      call mxCopyReal8ToPtr(dble(slon),slon_pr,size)     
      call mxCopyReal8ToPtr(dble(zenang),zenang_pr,size)     
      call mxCopyReal8ToPtr(dble(azang),azang_pr,size)     

      return
      end

