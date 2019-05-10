      subroutine idlConstNoSCF(argc,argv)
c     IDL call must be lower case ('idlconstinit_') even though upper is
c     used here

      Implicit None
      Integer*4 argc,argv(*),j
       	
      j = LOC(argc)
      
      call idlConstNoSCF_core(%val(argv(1)),
     +     %val(argv(2)),
     +     %val(argv(3)),
     +     %val(argv(4)),
     +     %val(argv(5)),
     +     %val(argv(6)))
      return
      end

c*****************************************************************************
      subroutine idlConstNoSCF_core(ifcoefsname,nfcoefsname,
     $     Frq,pol,earthRadius,nchanIn)
c     Based on mwsensors01/generic/src/sim.f90
      
      use oss_mw_module
      use constants, only: ErthRad

      implicit none
      INCLUDE "mxdims.incl"
      
      integer nchanIn
      real Frq(nchanIn),earthRadius
      integer pol(nchanIn)
      integer i,nlev,nchanOSS,nfcoefsname,ifcoefsname(nfcoefsname)
      INTEGER,PARAMETER  :: U_osscoefs    = 11
      CHARACTER(LEN=200) :: F_osscoefs
      
      if (nfcoefsname.gt.200) then
         print *,'err[idlConstNoSCF::idlConstNoSCF_core]: ',
     $        ' character array exceeds max dimension.'
         call errorHalt(1)
      endif
      do i=1,nfcoefsname
         F_osscoefs(i:i) = char(ifcoefsname(i))
      enddo

c     Get CMIS channel IDs for nchan OSS channels
c     Later we will expect iBeamSCF to match OSS/inSim channel set
      call getOSSselMW(U_osscoefs,F_osscoefs,nchanOSS=nchanOSS,
     $     cFreq=Frq(1:nchanIn),ipol=pol(1:nchanIn))

c     Check that nchanOSS matches nchanIn
      if (nchanOSS /= nchanIn) then
         print *,'err[idlConsInit]:  nchanOSS ~= nchanIn '
         print *,nchanOSS,nchanIn
         call errorHalt(1)
      endif

      earthRadius=ErthRad
      
      return
      end
