      subroutine idlConstInit(argc,argv)
c     IDL call must be lower case ('idlconstinit_') even though upper is
c     used here

      Implicit None
      Integer*4 argc,argv(*),j
       	
      j = LOC(argc)
      
      call idlConstInit_core(%val(argv(1)),
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
     +     %val(argv(18)))
      return
      end

c*****************************************************************************
      subroutine idlConstInit_core(ifcoefsname,nfcoefsname,
     $     ifoptdpthname,nfoptdpthname,Frq,pol,
     $     ifScanOrbitName,nfScanOrbitName,
     $     ifPointingName,nfPointingName,
     $     ifCalibName,nfCalibName,
     $     earthRadius,orbitAlt,iBeamOSS,nBeamOSS,nchanIn,eiaBeamOSS)
c     Based on mwsensors01/generic/src/sim.f90
      
      use ReadStdInputs
      use oss_mw_module
      use SCFread
      use constants

      implicit none
      INCLUDE "mxdims.incl"
      
      integer nchanIn
      real Frq(nchanIn),earthRadius,orbitAlt,eiaSCF(nchanIn)
      real nadirrad,altplusR,eiaBeamOSS(mxchan)
      integer pol(nchanIn)
      integer i,nchanOSS,nfcoefsname,ifcoefsname(nfcoefsname)
      integer nfoptdpthname,ifoptdpthname(nfoptdpthname)
      integer nfScanOrbitName,ifScanOrbitName(nfScanOrbitName)
      integer nfPointingName,ifPointingName(nfPointingName)
      integer nfCalibName,ifCalibName(nfCalibName)
      integer nBeam,nBeamCal,iBeamOSS(nchanIn),iBeamOSS2iBeamSCF(mxchan)
      integer ossID2scfID(mxchan),nchanSCF,ioss,iscf,nBeamOSS
      integer ossID2iBeamSCF(mxchan)
      character(len=200) :: F_ScanOrbit,F_Pointing,F_Cal
      integer, parameter :: U_ScanOrbit=20,lID=OSS_LEN_ID
      integer, parameter :: U_Pointing=50,U_Cal=51
      character(lID), dimension(mxchan) :: ossID
      character(lID), dimension(:), pointer :: beamID,beamIDcal,scfID
      integer, dimension(:), pointer :: iBeamSCF
      real, dimension(:), pointer :: nadirSCF
      
      if (nfcoefsname.gt.200.or.nfoptdpthname.gt.200
     $     .or.nfScanOrbitName.gt.200.or.nfCalibName.gt.200
     $     .or.nfPointingName.gt.200) then
         print *,'err[idlConstInit::idlConstInit_core]: ',
     $        ' character array exceeds max dimension.'
         call errorHalt(1)
      endif
      do i=1,nfcoefsname
         F_osscoefs(i:i) = char(ifcoefsname(i))
      enddo
      do i=1,nfScanOrbitName
         F_ScanOrbit(i:i) = char(ifScanOrbitName(i))
      enddo
      do i=1,nfPointingName
         F_Pointing(i:i) = char(ifPointingName(i))
      enddo
      do i=1,nfCalibName
         F_Cal(i:i) = char(ifCalibName(i))
      enddo

c     Get CMIS channel IDs for nchan OSS channels
c     Later we will expect iBeamSCF to match OSS/inSim channel set
      call getOSSselMW(U_osscoefs,F_osscoefs,nchanOSS=nchanOSS,
     $     cFreq=Frq(1:nchanIn),ipol=pol(1:nchanIn),
     $     chanID=ossID(1:nchanIn))

c     Check that nchanOSS matches nchanIn
      if (nchanOSS /= nchanIn) then
         print *,'err[idlConsInit]:  nchanOSS ~= nchanIn '
         print *,nchanOSS,nchanIn
         call errorHalt(1)
      endif
      
c     Get scan/orbit sensor constants
      call getSCFscanOrbit(U_ScanOrbit,F_ScanOrbit,EarthRad=earthRadius,
     $     scAlt=orbitAlt,nbeam=nBeam,beamID=beamID)

c     Get pointing sensor constants per BEAM (e.g., nBeam)
      call getSCFpointing(U_Pointing,F_Pointing,patt_Nd=nadirSCF)

c     Get channel list, nchan, and iBeam from calibration SCF:
      call getSCFcalib(U_Cal,F_Cal,nchan=nchanSCF,chanID=scfID,
     $     iBeam=iBeamSCF,beamID=beamIDcal,nbeam=nbeamCal)

c1    Later we will expect iBeamSCF to match OSS/inSim channel set
c1      IF (nchanIn.ne.nchanSCF) THEN
c1         print*,'err[idlConstInit]: Inconsistency between # of chans '
c1         print*,' from inSim dictionary and calibration SCFs: '
c1         print*,'nchan in inSim dictionary: ',nchanIn
c1         print*,'nchan in Calibration file: ',nchanSCF
c1         call errorHalt(1)
c1      ENDIF
      IF (nBeam.ne.nbeamCal) THEN
         print*,'err[idlConstInit]: Inconsistency between # of beams '
         print*,' contained in Scan/Orbit and Calibration SCFs: '
         print*,'nbeams in Scan/Orbit file: ',nBeam
         print*,'nbeams in Calibration file: ',nbeamCal
         call errorHalt(1)
      ENDIF
      do i=1,nBeam
         if (beamID(i).ne.beamIDcal(i)) then
            print*,'err[idlConstInit]: Scan/Orbit beamID inconsistent '
            print*,'with Calib beam ID:'
            print*,'beam # | scan/orbit beamID | Calib beamID: '
            print*,i,beamID(i),beamIDcal(i)
            call errorHalt(1)
         endif
      enddo
      
c     Get CMIS SCF channel number matching each OSS channel (ossID2scfID)
c     Get CMIS SCF beam number matching each OSS channel
      ossID2scfID = 0
      ossID2iBeamSCF = 0
      do ioss=1,nchanOSS
         do iscf=1,nchanSCF
            if (ossID(ioss) == scfID(iscf)) then
               ossID2scfID(ioss) = iscf
               ossID2iBeamSCF(ioss) = iBeamSCF(iscf)
               exit
            endif
         enddo
         if (ossID2scfID(ioss) == 0) then
            print *,'err[idlConsInit]: ', 
     $           ' Could not find OSS channel ID in SCF list'
            print *,'ioss,ossID(ioss):  ',ioss,ossID(ioss)
            call errorHalt(1)
         endif
      enddo

c     Index each OSS channel to its own set of beam numbers
      iBeamOSS = 0
      iBeamOSS(1) = 1
c     and save mapping from each OSS beam to master SCF beam list
      iBeamOSS2iBeamSCF(1) = ossID2iBeamSCF(1)
      nBeamOSS = 1
      do ioss=2,nchanOSS
c     Check if channel's SCF-beam matches beam already in list.
c     If it is, then beam number is same as matching channel.
         do i=1,ioss-1
            if (ossID2iBeamSCF(ioss) == ossID2iBeamSCF(i)) then
               iBeamOSS(ioss) = iBeamOSS(i)
               exit
            endif
         enddo
c     If it's not, then add a new beam number.
         if (iBeamOSS(ioss) == 0) then
            nBeamOSS = nBeamOSS+1
            iBeamOSS(ioss) = nBeamOSS
            iBeamOSS2iBeamSCF(nBeamOSS) = ossID2iBeamSCF(ioss)
         endif
      enddo

c1    Later we will expect iBeamSCF to match OSS/inSim channel set
c1      nBeamOSS = nBeamCal
c1      iBeamOSS(1:nchanOSS) = iBeamSCF(1:nchanSCF)

c     Calculate the nominal eia for each beam. 
      altplusR = orbitAlt+earthRadius
      do i=1,nBeamOSS
         nadirrad = deg2rad*nadirSCF(iBeamOSS2iBeamSCF(i))
c1         nadirrad = deg2rad*nadirSCF(i)
         eiaBeamOSS(i) = asin(altplusR*sin(nadirrad)/earthRadius)
         eiaBeamOSS(i) = rad2deg*eiaBeamOSS(i)
      enddo
      
      return
      end
