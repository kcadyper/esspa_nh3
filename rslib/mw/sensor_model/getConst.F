#include "fintrf.h"
C
#if 0
C     
C     getConst.F
C     .F file needs to be preprocessed to generate .for equivalent
C     
#endif
C     
C      
C     Matlab interface to sensor constants files and fortran constants.
C     2004/11/10:  John Galantowicz
C     Copyright AER, Inc. 2004
      
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on the DEC Alpha
C     64-bit platform
C
      
      use ReadStdInputs
      use SCFread
      use Constants

      implicit none

      integer plhs(*), prhs(*)
      mwpointer mxGetPr, mxCreateDoubleMatrix
      mwpointer out_pr(50),mxCreateString
C-----------------------------------------------------------------------
C             

      integer nlhs, nrhs, DESIGNLEN, IDLEN, iBeam, iChan
      parameter (DESIGNLEN = 20, IDLEN = 12)
      integer mxGetM, mxGetN, mxIsNumeric, mxGetString
      integer m, n, size, status, i
      integer U_ScanOrbit,U_Derived,U_Pointing,U_Cal
      parameter(U_ScanOrbit=20,U_Derived=30,U_Pointing=40,U_Cal=50)
c     Matlab is real*8 but subroutines are real*4
      real*4 orbitIncl,orbitAlt,earthRadius,muEarthKm
      real*4 rotSpeed,bore_az_t00,bore_az_t0
      real*4, dimension(:), pointer :: scan_start_time,integ_time
      real*4, dimension(:), pointer :: sample_time,sample_freq
      real*4 s2sc_roll,s2sc_pitch,s2sc_yaw,a2s_yaw, a2s_roll,a2s_pitch
      real*4, dimension(:), pointer :: beam_nadir0,bore2beam_az0
      real*4 beam_nadir,bore2beam_az
      real*4 bore2losAz,losNadir,fscaleinc,sc_yaw
      real*4 aview_spec,ascan_spec,nedt_spec,aview_nom,ascan_nom
      real*4 a2sRPY(3),s2scRPY(3)

      integer, dimension(:), pointer :: k_cal_cold,n_samples,k_cal_warm
      integer, dimension(:), pointer :: iBeams
      integer nBeam,rotDir,nChan
      
      character*200 F_ScanOrbit,F_Derived,F_Pointing,F_Cal
      character(IDLEN) chanIDIn,beamIDOut
      character(DESIGNLEN) design
      character(IDLEN), dimension(:), pointer :: beamIDs,chanIDs

C     Check for proper number of arguments. 
      if (.not.(
     $     (nrhs .eq. 0 .and. nlhs .eq. 3) 
     $     .or. 
     $     (nrhs .eq. 1 .and. nlhs .eq. 6) 
     $     .or. 
     $     (nrhs .eq. 5 .and. nlhs .eq. 28))) then
         call mexErrMsgTxt('getConst: 0 inputs -> 3 outputs or'//
     $        ' 1 input->6 outputs or '//
     $        ' 5 inputs->28 outputs required.')
      endif

C----------------------------------------------------

C     Get the size of the input array.
c      m = mxGetM(prhs(1))
c      n = mxGetN(prhs(1))
c      size = m*n
c      if (size.ne.1) then
c         call mexErrMsgTxt('getlatlon.F: Only supports array size=1')
c      endif
C----------------------------------------------------

C     Check to insure the input is a number.

c      if(mxIsNumeric(prhs(1)) .eq. 0) then
c         call mexErrMsgTxt('Input must be a number.')
c      endif

C----------------------------------------------------

C     Create matrix for the return arguments.
      do i=1,nlhs
         plhs(i) = mxCreateDoubleMatrix(1,1,0)
      enddo

c     Get input pointers

      if (nrhs .gt. 0) then
         status = mxGetString(prhs(1),F_ScanOrbit,200)
      endif
      if (nrhs .gt. 1) then
         status = mxGetString(prhs(2),F_Derived,200)
         status = mxGetString(prhs(3),F_Pointing,200)
         status = mxGetString(prhs(4),F_Cal,200)
         status = mxGetString(prhs(5),chanIDIn,IDLEN)
      endif
      
c     Get output pointers.  Note: plhs(8) and plhs(23) get redefined later
      do i=25,26
         plhs(i) = mxCreateDoubleMatrix(1,3,0)
      enddo
      do i=1,nlhs
         out_pr(i) = mxGetPr(plhs(i))
      enddo

c     Convert input pointers to values
c      call mxCopyPtrToReal8(scone_pr,scone,size)

      if (nrhs .gt. 0) then
         call getSCFscanOrbit(iu=U_ScanOrbit,fname=F_ScanOrbit,
     $        nbeam=nBeam)
         call getSCFscanOrbit(iu=U_ScanOrbit,fname=F_ScanOrbit,
     $        rotDir=rotDir,scanRate=rotSpeed,azRef=bore_az_t00,
     $        EarthRad=earthRadius,scAlt=orbitAlt,
     $        orbitIncl=orbitIncl,beamID=beamIDs,
     $        nsamp=n_samples,sTimeEv=scan_start_time,
     $        ncssamp=k_cal_cold,
     $        nwlsamp=k_cal_warm,tInteg=integ_time,
     $        sampInt=sample_time,sampFrq=sample_freq)
      else
         EarthRadius = ErthRad
      endif

c     Calculate muEarthKm in km from G and MEarth in Constants module
      muEarthKm = gravCnst*ErthMass/1000**3

C     Load the data into y_pr, which is the output to MATLAB
      size = 1
      if (nrhs .eq. 0) then
         call mxCopyReal8ToPtr(dble(earthRadius),out_pr(1),size)
         call mxCopyReal8ToPtr(dble(muEarthKm),out_pr(2),size)
         call mxCopyReal8ToPtr(dble(MISSING_REAL),out_pr(3),size)
      else
         call mxCopyReal8ToPtr(dble(orbitIncl),out_pr(1),size)
         call mxCopyReal8ToPtr(dble(orbitAlt),out_pr(2),size)
         call mxCopyReal8ToPtr(dble(earthRadius),out_pr(3),size)
         call mxCopyReal8ToPtr(dble(muEarthKm),out_pr(4),size)
         call mxCopyReal8ToPtr(dble(MISSING_REAL),out_pr(5),size)
         call mxCopyReal8ToPtr(dble(rotSpeed),out_pr(6),size)
      endif

      if (nrhs .gt. 1) then
         
c     Get iBeam and beamIDOut for chanIDIn

         call getSCFcalib(iu=U_Cal,fname=F_Cal,
     $        nchan=nChan,chanID=chanIDs,beamID=beamIDs,
     $        iBeam=iBeams)

         iChan = 0
         do i=1,nChan
            if (adjustl(chanIDs(i)) == adjustl(chanIDIn)) then
               ichan = i
               exit
            endif
         enddo
         if (iChan == 0) then
            call mexErrMsgTxt('getConst: Unmatched chanID.')
         endif
         iBeam = iBeams(iChan)
         beamIDOut = beamIDs(iBeam)

         call getSCFpointing(iu=U_Pointing,fname=F_Pointing,
     $        s2sc_roll=s2sc_roll,s2sc_pitch=s2sc_pitch,
     $        s2sc_yaw=s2sc_yaw,a2s_roll=a2s_roll,
     $        a2s_pitch=a2s_pitch,a2s_yaw=a2s_yaw,
     $        sc_yaw=sc_yaw,
     $        patt_Nd=beam_nadir0,patt_Az=bore2beam_az0)
         
cccc  Check if F_Derived was provided before loading
         if (LEN_TRIM(F_Derived) .gt. 3) then
            call getSCFAntDerived(iu=U_Derived,fname=F_Derived,
     $           chanIDIn=chanIDIn,design=design,
     $           bore2losAz=bore2losAz,losNadir=losNadir,
     $           fscaleinc=fscaleinc,aview_spec=aview_spec,
     $           ascan_spec=ascan_spec,nedt_spec=nedt_spec,
     $           aview_nom=aview_nom,ascan_nom=ascan_nom)
         else
            design = 'MISSING'
            bore2losAz = MISSING_REAL
            losNadir = MISSING_REAL
            fscaleinc = 1.0
            aview_spec = MISSING_REAL
            ascan_spec = MISSING_REAL
            nedt_spec = MISSING_REAL
            aview_nom = MISSING_REAL
            ascan_nom = MISSING_REAL
         endif

C     Skip RPY rotations -- Added separately in later applications
         bore2beam_az = bore2beam_az0(iBeam)
         beam_nadir = beam_nadir0(iBeam)
         bore_az_t0 = bore_az_t00-360.0

         a2sRPY(1) = a2s_roll
         a2sRPY(2) = a2s_pitch
         a2sRPY(3) = a2s_yaw
         s2scRPY(1) = s2sc_roll
         s2scRPY(2) = s2sc_pitch
         s2scRPY(3) = s2sc_yaw

C     Load the data into *_pr, which is the output to MATLAB
         size = 1
         call mxCopyReal8ToPtr(dble(bore_az_t0),out_pr(7),size)
         call mxCopyReal8ToPtr(dble(nchan),out_pr(8),size)
c         call mxCopyCharacterToPtr(design_array,out_pr(9),...)
         plhs(9) = mxCreateString(design)
         call mxCopyReal8ToPtr(dble(bore2beam_az),out_pr(10),size)
         call mxCopyReal8ToPtr(dble(scan_start_time(iBeam)),
     $        out_pr(11),size)
         call mxCopyReal8ToPtr(dble(integ_time(iBeam)),
     $        out_pr(12),size)
         call mxCopyReal8ToPtr(dble(sample_time(iBeam)),
     $        out_pr(13),size)
         call mxCopyReal8ToPtr(dble(n_samples(iBeam)),
     $        out_pr(14),size)
         call mxCopyReal8ToPtr(dble(beam_nadir),out_pr(15),size)
         call mxCopyReal8ToPtr(dble(aview_spec),out_pr(16),size)
         call mxCopyReal8ToPtr(dble(ascan_spec),out_pr(17),size)
         call mxCopyReal8ToPtr(dble(nedt_spec),out_pr(18),size)
         call mxCopyReal8ToPtr(dble(fscaleinc),
     $        out_pr(19),size)
         call mxCopyReal8ToPtr(dble(bore2losAz),out_pr(20),size)
         call mxCopyReal8ToPtr(dble(losNadir),out_pr(21),size)
         call mxCopyReal8ToPtr(dble(k_cal_cold(iBeam)),
     $        out_pr(22),size)
         plhs(23) = mxCreateString(beamIDOut)
         call mxCopyReal8ToPtr(dble(iBeam),out_pr(24),size)
         call mxCopyReal8ToPtr(dble(aview_nom),out_pr(27),size)
         call mxCopyReal8ToPtr(dble(ascan_nom),out_pr(28),size)
         size = 3
         call mxCopyReal8ToPtr(dble(a2sRPY),out_pr(25),size)
         call mxCopyReal8ToPtr(dble(s2scRPY),out_pr(26),size)
         
      endif

      return
      end

