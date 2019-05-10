MODULE CalibSim

  use PlanckTabMW

  implicit none
  private
  public:: initCalSim,calSimScan,getSensConst,saveOrigSC,errorSC
  public:: gainVarChan,NULLF

  integer,            parameter :: lID=12  ! length of ID strings
  character(len=200), parameter :: NULLF='NULL',BLANKF=''

!-- Global private data

  integer, dimension(:),        pointer :: iBeam
  integer                               :: iseedVC=9278473
  integer                               :: iseedVW=3479274
  integer                               :: iseedVE=3479274

! ScanOrbit sensor constants
  real                                  :: scanRate
  character(lID), dimension(:), pointer :: beamIDSO
  real,    dimension(:),        pointer :: sampInt
! Calib sensor constants
  integer                               :: nbeam
  character(lID), dimension(:), pointer :: beamID
  integer ,dimension(:),        pointer :: ncssamp,nwlsamp 
  real,    dimension(:,:),      pointer :: spocs       ,spocs0
  real,    dimension(:,:),      pointer :: spow        ,spow0
  real,    dimension(:),        pointer :: emis_refcs  ,emis_refcs0
  real,    dimension(:),        pointer :: emis_w      ,emis_w0
  integer                               :: nsccompc
  integer, dimension(:),        pointer :: sccompcID
  real,    dimension(:),        pointer :: emis_scc
  real,    dimension(:,:,:),    pointer :: vf_sccs     ,vf_sccs0
  real,    dimension(:,:,:),    pointer :: vf_scw      ,vf_scw0
  integer                               :: nsencompc
  integer, dimension(:),        pointer :: sencompcID
  real,    dimension(:),        pointer :: emis_senc
  real,    dimension(:,:,:),    pointer :: vf_sencs    ,vf_sencs0
  real,    dimension(:,:,:),    pointer :: vf_senw     ,vf_senw0
  real                                  :: tbcos
  real                                  :: tminw
  real                                  :: tmaxw
  real                                  :: tmaxwv
  integer                               :: nfft
  integer                               :: ifft
  integer                               :: nthmrps
  real,    dimension(:,:),      pointer :: dtthrm
! Calib Error sensor constants
  real,    dimension(:),        pointer :: spocsErr
  real,    dimension(:),        pointer :: spowErr
  real,    dimension(:),        pointer :: nonlin
  real,    dimension(:),        pointer :: gainVar
! Noise sensor constants
  integer                               :: nchan
  character(lID), dimension(:), pointer :: chanID
! Antenna View sensor constants
  integer, dimension(:),        pointer :: nsamp
  real,    dimension(:),        pointer :: spoev       ,spoev0
  real,    dimension(:),        pointer :: emis_ref    ,emis_ref0
  integer, dimension(:),        pointer :: refMap
  integer                               :: nsccomp
  integer, dimension(:),        pointer :: sccompID
  real,    dimension(:),        pointer :: emis_sc
  real,    dimension(:,:,:),    pointer :: vf_scev     ,vf_scev0
  integer                               :: nsencomp
  integer, dimension(:),        pointer :: sencompID
  real,    dimension(:),        pointer :: emis_sen
  real,    dimension(:,:,:),    pointer :: vf_senev    ,vf_senev0
! Antenna View Error sensor constants
  real,    dimension(:),        pointer :: spoevErr

! Actual errors by channel
  real,    dimension(:),    allocatable :: nonlinChan
  real,    dimension(:),    allocatable :: gainVarChan


! Variables that are static for repeated calls of calSimScan
  integer                              :: nchanParent
  integer, dimension(:),   allocatable :: chanMap ! parent list to full list
  integer, dimension(:),   allocatable :: nsampParent
  real,    dimension(:),   allocatable :: waven   ! wavenumber (cm-1) vector
  real,    dimension(:),   allocatable :: frots
! Arrays that are allocated at initialization in initCalSim
  real,    dimension(:),   allocatable :: rCosm
  real,    dimension(:,:), allocatable :: rSenc
  real,    dimension(:,:), allocatable :: rSen
  real,    dimension(:,:), allocatable :: rScftc
  real,    dimension(:,:), allocatable :: rScft

! public parameters
  integer, dimension(:), pointer,public :: ithmrw, ithmrc
  integer, dimension(:), pointer,public :: ithmrm, ithmrh
  integer, dimension(:), pointer,public :: ithmrsn,ithmrsnc
  integer, dimension(:), pointer,public :: ithmrsc,ithmrscc
  integer, dimension(:), pointer,public :: ithmrfft
  integer, public                       :: nthmrw, nthmrc
  integer, public                       :: nthmrm, nthmrh
  integer, public                       :: nthmrsc, nthmrscc
  integer, public                       :: nthmrsn, nthmrsnc
  integer, public                       :: nthmrfft
  integer, public                       :: nthmIn
  integer, public                       :: npnts, nthmr
  integer, dimension(:), pointer,public :: mapThm
  real,  dimension(:,:), pointer,public :: thermT ! nthmr x npnts

CONTAINS

  SUBROUTINE initCalSim(inChanID,freq,nsampPchan,nBeam,mapBeam,&
       ncssampP,nwlsampP)

    use constants, only: LtSp_cgs

!-- Calling arguments

    character(lID), dimension(:), intent(in) :: inChanID  ! parent's channels
    real, dimension(:),   intent(in) :: freq        ! parent's frequencies
    integer, dimension(:),intent(in) :: nsampPchan  ! #samp/chan for parent
    integer,              intent(in) :: nBeam
    integer, dimension(:),intent(in) :: mapBeam
    integer, dimension(:),pointer    :: ncssampP   ! #cold cal/beam for parent
    integer, dimension(:),pointer    :: nwlsampP   ! #warm cal/beam for parent

!-- Local variables

    integer            :: i,ich,j,iHold,iBeamPrev
    real               :: rotime

    real,           parameter :: GHzToHz=1.e9
    real,           parameter :: minTOms=60.e3
    character(lID), parameter :: NULL_PATTERN=repeat(' ',lID)
    if (allocated(waven)) deallocate(waven,frots,rCosm)
    allocate(waven(nchan),frots(nchan),rCosm(nchan))
    if (allocated(rSenc)) deallocate(rSenc,rSen,rScftc,rScft)
    allocate(rSenc(nchan,nsencompc), rSen(nchan,nsencomp), &
            rScftc(nchan,nsccompc), rScft(nchan,nsccomp) )

    nchanParent=size(inChanID)
    if (size(freq) /= nchanParent) then
      print *,'err[CalibSim::initCalSim]: Array size mismatch'
      print *,'size(freq),nchanParent:', &
               size(freq),nchanParent
      call errorHalt(1)
    endif
    if (allocated(chanMap)) deallocate(chanMap)
    allocate(chanMap(nchanParent))
    if (size(nsampPchan) /= nchanParent) then
      print *,'err[CalibSim::initCalSim]: Array size mismatch'
      print *,'size(nsampPchan),nchanParent:', &
               size(nsampPchan),nchanParent
      call errorHalt(1)
    endif
    if (allocated(nsampParent)) deallocate(nsampParent)
    allocate(nsampParent(nchanParent))
    nsampParent=nsampPchan
   
!-- Map from channels treated by parent program to full channel set

    chanMap=0
    do i=1,nchanParent
      do ich=1,nchan
        if (chanID(ich) == inChanID(i)) then
          chanMap(i)=ich
          exit
        endif
      enddo
      if (chanMap(i) == 0) then
        print *,'err[CalibSim::initCalSim]: Could not find channel ID in ', &
          'SCF list'
        print *,'i,inChanID(i):',i,inChanID(i)
        call errorHalt(1)
      endif
    enddo

!-- Find number of calibration samples for each beam treated in parent program 

    if (associated(ncssampP)) deallocate(ncssampP)
    allocate(ncssampP(nBeam))
    if (associated(nwlsampP)) deallocate(nwlsampP)
    allocate(nwlsampP(nBeam))
    do j=1,nBeam  
      iBeamPrev=0
      do i=1,nchanParent
        if (mapBeam(i) == j) then  ! map from beam to every channel using it
          iHold=iBeam(chanMap(i))  ! which is matching beam from full set
          if (iBeamPrev > 0 .and. iHold /= iBeamPrev) then
            print *,'err[CalibSim::initCalSim]: Parent has a pair of channels'
            print *,'on a common beam while SCF indicates separate beams'
            print *,'j,i,chanMap(i),iBeam(chanMap(i)),iBeamPrev:'
            print *, j,i,chanMap(i),iBeam(chanMap(i)),iBeamPrev
            call errorHalt(1)
          endif
          iBeamPrev=iHold
        endif
      enddo
      if (iBeamPrev == 0) then
        print *,'err[CalibSim::initCalSim]: Parent has a beam not used by ', &
          'any parent channels'
        print *,'j,mapBeam(1:nchanParent):',j,mapBeam(1:nchanParent)
        call errorHalt(1)
      endif
      ncssampP(j)=ncssamp(iBeamPrev)
      nwlsampP(j)=nwlsamp(iBeamPrev)
    enddo

!-- Convert frequency (GHz) to wavenumber (cm-1), for use in Planck formula

    do i=1,nchanParent
      waven(chanMap(i))=freq(i)*GHzToHz/LtSp_cgs
    enddo

!-- Fraction of rotation per scan interval (for gain variation)

    rotime=minTOms/scanRate      ! RPM to time for one rotation (ms)

    do ich=1,nchanParent
       frots(chanMap(ich))=sampInt(iBeam(chanMap(ich)))/rotime
    enddo
    do i=1,nchanParent
      if (nsampParent(i) > nsamp(ibeam(chanMap(i)))) then
        print *,'err[CalibSim::initCalSim]: Number of samples per scan in ', &
          'parent exceeds SCF number'
        print *,'i,chanMap(i),ibeam(chanMap(i)),nsampParent(i),', &
          'nsamp(ibeam(chanMap(i))):'
        print *, i,chanMap(i),ibeam(chanMap(i)),nsampParent(i),  &
           nsamp(ibeam(chanMap(i)))
        call errorHalt(1)
      endif

    enddo
        
    return

  END SUBROUTINE initCalSim

!----------------------------------------------------------------------------

  SUBROUTINE calSimScan(tEview,dt_fbw,tWtru,vCtru,vWtru, &
    t_refcs,t_reflf,t_refhf,tSenc,tSen,tScftc,tScft, &
    vEmsmt,vCmsmt,vWmsmt, &
    cPRTw,cPRTc,cPRTm,cPRTh,cPRTsenc,cPRTsen,cPRTscc,cPRTsc)

! Simulate the calibration of a scan line of sensor data.  The output is
! a set of simulated raw data record (RDR) parameters for a scan line.
!
! All temperatures are units of kelvin
! All radiances are units mW m-2 cm str-1
!
    use RTFquadratic
    use constants, only: CosmBckg,MISSING_REAL
    use SensorData, only: var_nedt

    real,    parameter :: lowTE=1.01,hiTE=400.! for screening flagged input data
    real,    parameter :: timeset=0.,bandwset=0.  ! For var_nedt
    integer, parameter :: icaset=0, icnvBW=1      ! For var_nedt

!-- Calling arguments

    real, dimension(:,:), intent(in)    :: tEview
    real,                 intent(in)    :: dt_fbw   ! front/back dtemp warm load
    real,                 intent(in)    :: tWtru    ! true temp of warm load
    real, dimension(:),   intent(in)    :: vCtru    ! true count of cold scene
    real, dimension(:),   intent(in)    :: vWtru    ! true count of warm scene
    real,                 intent(in)    :: t_refcs  ! true temp cold reflector
    real,                 intent(in)    :: t_reflf  ! true temp main reflector
    real,                 intent(in)    :: t_refhf  ! true temp hi-frq reflector
    real, dimension(:),   intent(in)    :: tSenc
    real, dimension(:),   intent(in)    :: tSen
    real, dimension(:),   intent(in)    :: tScftc
    real, dimension(:),   intent(in)    :: tScft
    real, dimension(:,:), intent(inout) :: vEmsmt
    real, dimension(:,:), intent(inout) :: vCmsmt
    real, dimension(:,:), intent(inout) :: vWmsmt
    integer, dimension(:,:), intent(inout) :: cPRTw
    integer, dimension(:,:), intent(inout) :: cPRTc
    integer, dimension(:,:), intent(inout) :: cPRTm
    integer, dimension(:,:), intent(inout) :: cPRTh
    integer, dimension(:,:), intent(inout) :: cPRTsenc
    integer, dimension(:,:), intent(inout) :: cPRTsen
    integer, dimension(:,:), intent(inout) :: cPRTscc
    integer, dimension(:,:), intent(inout) :: cPRTsc

!-- Local variables

    integer :: ich,ibm,icmp,isam,ith,its,mch,ical,iEd,iEr,limsam
    integer,dimension(1) :: mxval1,mxval2
    real    :: gauss     ! external library subprograms
    real    :: rotime
    real    :: rCamp,rWamp,rEamp,tCamp,tWamp,tEamp
    real    :: tWface,rWface,rCface,rELface,rEHface,rEface
    real    :: tCnoise,tWnoise,tEnoise,rCnoise,rWnoise,rEnoise
    real    :: rWall,rCall,rCref,vWend,tEall,rEall,rEref
    real    :: spoW_net,spoC_net,spoE_net,ref_net
    real    :: tmidlin,tmidqua,rmidlin,rmidqua,drmax
    real    :: rCmsmt,rWmsmt,rEmsmt,rEview
    real    :: RTFc0s,RTFc1s,RTFc2s,RTFc0e,RTFc1e,RTFc2e
    real    :: vEs,vEe
    character(lID) :: chID1

!-- Assume sensor component 1 is top of bucket and warm load reflection is
!-- from top of bucket
    integer, parameter :: isenW=1
!-- Assume sensor component 2 is reflector structure, seen directly from feed
    integer, parameter :: nsenEd=1
    integer, parameter, dimension(nsenEd) :: isenEd=(/2/)
!-- Assume sensor component 3 is edge of body, seen through reflecter
    integer, parameter :: nsenEr=1
    integer, parameter, dimension(nsenEr) :: isenEr=(/3/)

!-- Check inputs

    if (size(vCtru) /= nchan .or. size(vWtru) /= nchan) then
      print *,'err[CalibSim::calSimScan]: Array size mismatch'
      print *,'size(vCtru),size(vWtru),nchan:',size(vCtru),size(vWtru),nchan
      call errorHalt(1)
    endif
    if (size(tEview,1) /= nchanParent .or. size(vEmsmt,1) /= nchanParent .or. &
        size(vCmsmt,1) /= nchanParent .or. size(vWmsmt,1) /= nchanParent) then
      print *,'err[CalibSim::calSimScan]: Array size mismatch'
      print *,'size(tEview,1),size(vEmsmt,1),size(vCmsmt,1),size(vWmsmt,1),', &
        'nchanParent:', &
               size(tEview,1),size(vEmsmt,1),size(vCmsmt,1),size(vWmsmt,1), &
         nchanParent
      call errorHalt(1)
    endif

    mxval1 = maxval(nsampParent(:))
    if (size(tEview,2) /= mxval1(1) .or. &
         size(vEmsmt,2) < mxval1(1)) then
      print *,'err[CalibSim::calSimScan]: Array size mismatch'
      print *,'size(tEview,2),size(vEmsmt,2),maxval(nsampParent(:)):', &
        size(tEview,2),size(vEmsmt,2),mxval1(1)
      call errorHalt(1)
    endif
    
    mxval1 = maxval(ncssamp(:))
    mxval2 = maxval(nwlsamp(:))
    if ((size(vCmsmt,2) < mxval1(1) .or. &
         size(vWmsmt,2) < mxval2(1))) then
      print *,'err[CalibSim::calSimScan]: Array size too small'
      print *,'size(vCmsmt,2),maxval(ncssamp(:)):', &
        size(vCmsmt,2),mxval1(1)
      print *,'size(vWmsmt,2),maxval(nwlsamp(:)):', &
        size(vWmsmt,2),mxval2(1)
      call errorHalt(1)
    endif
    
    do mch=1,nchanParent
       ibm=ibeam(chanMap(mch))
       do ical=1,ncssamp(ibm)
          if (vf_sencs(ibm,ical,isenW) > 0.) then
             print *,'err[CalibSim::calSimScan]: The part of the ', &
                  'sensor seen in reflection off warm load', &
                  'should not be seen in direct cold sky spillover.'
             print *,'ibm,ical,isenW,vf_senw(ibm,ical,isenW),', &
                  'vf_sencs(ibm,ical,isenW):'
             print *, ibm,ical,isenW,vf_senw(ibm,ical,isenW), &
                  vf_sencs(ibm,ical,isenW)
             call errorHalt(1)
          endif
       enddo
       do ical=1,ncssamp(ibm)
          if (vf_senw(ibm,ical,isenW) > 0.) then
             print *,'err[CalibSim::calSimScan]: The part of the ', &
                  'sensor seen in reflection off warm load', &
                  'should not be seen in direct warm load spillover.'
             print *,'ibm,ical,isenW,vf_senw(ibm,ical,isenW),', &
                  'vf_sencs(ibm,ical,isenW):'
             print *, ibm,ical,isenW,vf_senw(ibm,ical,isenW), &
                  vf_sencs(ibm,ical,isenW)
             call errorHalt(1)
          endif
       enddo
    enddo
    if (any(isenEd(1:nsenEd) > nsencomp) .or. &
      any(isenEr(1:nsenEr) > nsencomp)) then
      print *,'err[CalibSim::calSimScan]: Sensor IDs out of range'
      print *,'nsencomp,isenEd(1:nsenEd),isenEr(1:nsenEr):', &
               nsencomp,isenEd(1:nsenEd),isenEr(1:nsenEr)
      call errorHalt(1)
    endif

!-- Thermometric measurement of warm calibration target and cold-sky reflector

    do its=1,nthmrps  
      do ith=1,nthmrw
        cPRTw(ith,its)=getCounts(tWtru,mapThm(ithmrw(ith)) )
      enddo
      do ith=1,nthmrc
        cPRTc(ith,its)=getCounts(t_refcs,mapThm(ithmrc(ith)) )
      enddo

!-- Thermometric measurement of main and high-frequency reflectors
      do ith=1,nthmrm
        cPRTm(ith,its)=getCounts(t_reflf,mapThm(ithmrm(ith)) )
      enddo
      do ith=1,nthmrh
        cPRTh(ith,its)=getCounts(t_refhf,mapThm(ithmrh(ith)) )
      enddo
    enddo

!-- Thermometric measurement of sensor and spacecraft components in view

    if (size(tSenc) < nsencompc .or. size(tSen) < nsencomp .or. &
        size(tScftc) < nsccompc .or. size(tScft) < nsccomp) then 
      print *,'err[CalibSim::calSimScan]: Component temperature array size', &
        ' mismatch'
      print *,'size(tSenc),nsencompc,size(tSen),nsencomp,', &
        'size(tScftc),nsccompc,size(tScft),nsccomp:'
      print *, size(tSenc),nsencompc,size(tSen),nsencomp,  &
         size(tScftc),nsccompc,size(tScft),nsccomp
      call errorHalt(1)
    endif

    do its=1,nthmrps  
      do icmp=1,nsencompc
        cPRTsenc(icmp,its)=getCounts(tSenc(icmp),mapThm(ithmrsnc(icmp)) )
      enddo
      do icmp=1,nsencomp
        cPRTsen(icmp,its)=getCounts(tSen(icmp), mapThm(ithmrsn(icmp)) ) 
      enddo
      do icmp=1,nsccompc
        cPRTscc(icmp,its)=getCounts(tScftc(icmp), mapThm(ithmrscc(icmp)) )
      enddo
      do icmp=1,nsccomp
        cPRTsc(icmp,its)=getCounts(tScft(icmp), mapThm(ithmrsc(icmp)) )
      enddo
    enddo

!-- Radiances of sources of extraneous radiation

    do mch=1,nchanParent
      ich=chanMap(mch)
      chID1 = chanID(ich)
      rCosm(ich) = fwdPTM(chID1=chID1, tin=CosmBckg)       ! cosmic

      do icmp=1,nsencompc
        rSenc(ich,icmp) = fwdPTM(chID1,tSenc(icmp))*emis_senc(icmp) + &
           rCosm(ich)*(1.-emis_senc(icmp))                     ! sensor
      enddo

      do icmp=1,nsencomp
        rSen(ich,icmp) = fwdPTM(chID1,tSen(icmp))*emis_sen(icmp) + &
          rCosm(ich)*(1.-emis_sen(icmp))                      ! sensor
      enddo

      do icmp=1,nsccompc
        rScftc(ich,icmp) = fwdPTM(chID1,tScftc(icmp))*emis_scc(icmp) + &
          rCosm(ich)*(1.-emis_scc(icmp))                      ! spaceraft
      enddo

      do icmp=1,nsccomp
        rScft(ich,icmp) = fwdPTM(chID1,tScft(icmp))*emis_sc(icmp) + &
          rCosm(ich)*(1.-emis_sc(icmp))                       ! spaceraft
      enddo
    enddo  ! EOL for nchanParent

    do mch=1,nchanParent
      ich=chanMap(mch)
      chID1  = chanID(ich)

!-- Face of warm calibration target

      tWface = tWtru +dt_fbw
      rWface = fwdPTM(chID1,tWface)

!-- Surface of cold calibration reflector
      rCface   = fwdPTM(chID1,t_refcs)

!-- Surface of Earth view reflectors
      rELface   = fwdPTM(chID1,t_reflf)
      rEHface   = fwdPTM(chID1,t_refhf)

!-- Noise amplitudes for calibration radiometric measurements
!--  - Assume target temperature ~= radiometric temperature of target

      tCamp=var_nedt(ich,CosmBckg,timeset,bandwset,icaset,icnvBW)

      tWamp=var_nedt(ich,tWface,timeset,bandwset,icaset,icnvBW)

      tWnoise=tWface+tWamp
      rWnoise=fwdPTM(chID1, tWnoise)
 
      rWamp=rWnoise-rWface
      rCamp = rWamp*tCamp/tWamp

      mxval1 = maxval((/ncssamp(ibeam(ich)),nwlsamp(ibeam(ich))/))

      !-- Loop over calibration samples 
      do ical=1,mxval1(1)
         
!-- Radiance from all components of warm target scene
         
         if (ical.le.nwlsamp(ibeam(ich))) then 
            
            ! load emission:
            rWall=rWface*(1.-spow(ibeam(ich),ical))*emis_w(ibeam(ich)) + &
            ! reflection off load from top of bucket:
                 rSenc(ich,isenW)* &
                 (1.-spow(ibeam(ich),ical))*(1.-emis_w(ibeam(ich)))
            
            spoW_net=spow(ibeam(ich),ical)
            
            ! spillover from sensor:
            do icmp=1,nsencompc
               rWall=rWall + rSenc(ich,icmp)*vf_senw(ibeam(ich),ical,icmp)
               spoW_net=spoW_net-vf_senw(ibeam(ich),ical,icmp)
            enddo
            
            ! spillover from spacecraft:
            do icmp=1,nsccompc
               rWall=rWall + rScftc(ich,icmp)*vf_scw(ibeam(ich),ical,icmp)
               spoW_net=spoW_net-vf_scw(ibeam(ich),ical,icmp)
            enddo
            ! spillover from space (direct):

            if (spoW_net < 0.) then
               print *,'err[CalibSim::calSimScan]: Warm load sensor', &
                    ' constants inconsistent'
               print*,'ich,ical,warm spillover record,spoW_net: ' , &
                    ich,ical,spow(ibeam(ich),ical),spoW_net
               call errorHalt(1)
            endif
            rWall=rWall+rCosm(ich)*spoW_net
            
         endif
!-- Radiance from all components of cold target scene
         if (ical.le.ncssamp(ibeam(ich))) then 
            
            ! spillover from sensor:
            spoC_net=spocs(ibeam(ich),ical)
            rCall=0.
            do icmp=1,nsencompc
               rCall=rCall + rSenc(ich,icmp)*vf_sencs(ibeam(ich),ical,icmp)
               spoC_net=spoC_net-vf_sencs(ibeam(ich),ical,icmp)
            enddo
            
            ! spillover from space (direct):
            if (spoC_net < 0.) then
               print *,'err[CalibSim::calSimScan]: Cold load sensor', &
                    ' constants inconsistent'
               print *,'ich,ical,spocs(ibeam(ich),ical),ical),spoC_net:', &
                    ich,ical,spocs(ibeam(ich),ical),spoC_net
               call errorHalt(1)
            endif
            rCall=rCall+rCosm(ich)*spoC_net
            
            ! radiance at cold sky reflector from spacecraft:
            ref_net=1.
            rCref=0.
            do icmp=1,nsccompc
               rCref=rCref + rScftc(ich,icmp)*vf_sccs(ibeam(ich),ical,icmp)
               ref_net=ref_net-vf_sccs(ibeam(ich),ical,icmp)
            enddo
            ! radiance at cold sky reflector from space:
            if (ref_net < 0.) then
               print *,'err[CalibSim::calSimScan]: Cold load sensor', &
                    ' constants inconsistent'
               print *,'ich,ical,ref_net:', &
                    ich,ical,ref_net
               call errorHalt(1)
            endif
            
            rCref=rCref+rCosm(ich)*ref_net
            rCall=rCall+rCref*(1.-spocs(ibeam(ich),ical))* &
                 (1.-emis_refcs(ibeam(ich)))
            
            ! reflector emission:
            rCall=rCall+rCface*(1.-spocs(ibeam(ich),ical))*emis_refcs(ibeam(ich))
            
         endif
         
!-- Define true radiometer transfer function on basis of first cal sample,
!-- once at beginning and once at end of scan period, accounting for gain
!-- variation
            
         if (ical == 1) then
            tmidlin=(CosmBckg+tWtru)/2.
            tmidqua=tmidlin+nonlinChan(ich)

            rmidlin = fwdPTM(chID1,tmidlin)
            rmidqua = fwdPTM(chID1,tmidqua)
            drmax = rmidqua - rmidlin
            call RTFquadraticFit(drmax,rCall,rWall,vCtru(ich),vWtru(ich), &
                 RTFc0s,RTFc1s,RTFc2s)
            vWend=vWtru(ich)+gainVarChan(ich)*(vWtru(ich)-vCtru(ich))
            call RTFquadraticFit(drmax,rCall,rWall,vCtru(ich),vWend, &
                 RTFc0e,RTFc1e,RTFc2e)
         endif
         
         !-- Radiometric measurement of calibration targets
         
         if (ical.le.ncssamp(ibeam(ich))) then 
            rCmsmt=rCall +gauss(iseedVC,rCamp,0.)
            vCmsmt(mch,ical)=RTFquadraticVal(rCmsmt,RTFc0s,RTFc1s,RTFc2s)
         endif
         if (ical.le.nwlsamp(ibeam(ich))) then 
            rWmsmt=rWall +gauss(iseedVW,rWamp,0.)
            vWmsmt(mch,ical)=RTFquadraticVal(rWmsmt,RTFc0s,RTFc1s,RTFc2s)
         endif

      enddo

!     Assume parent samples are the first positions of a full scan
      limsam=min(nsamp(ibeam(ich)),nsampParent(mch))
      do isam=1,limsam   !-- Loop over Earth view samples

!-- Radiance from all components of Earth scene
         
         ! check for unphysical (i.e. flagged as missing) Tb:
        if (tEview(mch,isam).lt.lowTE .OR. tEview(mch,isam).gt.hiTE)then 
           vEmsmt(mch,isam) = MISSING_REAL
           cycle
        endif
         
        ! spillover direct from sensor:
        spoE_net=spoev(ibeam(ich))
        rEall=0.
        do iEd=1,nsenEd
          rEall=rEall+rSen(ich,isenEd(iEd))* &
            vf_senev(ibeam(ich),isam,isenEd(iEd))
          spoE_net=spoE_net-vf_senev(ibeam(ich),isam,isenEd(iEd))
        enddo
        ! spillover from space (direct):
        if (spoE_net < 0.) then
          print *,'err[CalibSim::calSimScan]: Earth view sensor constants', &
            ' inconsistent'
          print *,'ich,isam,spoev(ibeam(ich)),spoE_net:', &
            ich,isam,spoev(ibeam(ich)),spoE_net
          call errorHalt(1)
        endif
        rEall=rEall+rCosm(ich)*spoE_net
        ! radiance at Earth view reflector from sensor:
        ref_net=1.
        rEref=0.
        do iEr=1,nsenEr
          rEref=rEref+rSen(ich,isenEr(iEr))* &
            vf_senev(ibeam(ich),isam,isenEr(iEr))
          ref_net=ref_net-vf_senev(ibeam(ich),isam,isenEr(iEr))
        enddo
        ! radiance at Earth view reflector from spacecraft:
        do icmp=1,nsccomp
          rEref=rEref + rScft(ich,icmp)*vf_scev(ibeam(ich),isam,icmp)
          ref_net=ref_net-vf_scev(ibeam(ich),isam,icmp)
        enddo
        ! radiance at Earth view reflector from Earth scene:
        if (ref_net < 0.) then
          print *,'err[CalibSim::calSimScan]: Earth view sensor constants', &
            ' inconsistent'
          print *,'ich,isam,ref_net:', &
            ich,isam,ref_net
          call errorHalt(1)
        endif

        rEview  = fwdPTM(chID1,tEview(mch,isam))

        rEref=rEref+rEview*ref_net
        rEall=rEall+rEref*(1.-spoev(ibeam(ich)))*(1.-emis_ref(ibeam(ich)))
        ! reflector emission:
        rEface=rELface
        if (refMap(ibeam(ich)) == 2) rEface=rEHface
        rEall=rEall+rEface*(1.-spoev(ibeam(ich)))*emis_ref(ibeam(ich))

!-- Noise for Earth view measurement

        tEall = invPTM(chID1, rEall)
        tEamp=var_nedt(ich,tEall,timeset,bandwset,icaset,icnvBW)
        tEnoise=tEall+tEamp

        rEnoise = fwdPTM(chID1, tEnoise)

        rEamp=rEnoise-rEall
        rEmsmt=rEall +gauss(iseedVE,rEamp,0.)
        vEs=RTFquadraticVal(rEmsmt,RTFc0s,RTFc1s,RTFc2s)
        vEe=RTFquadraticVal(rEmsmt,RTFc0e,RTFc1e,RTFc2e)
        vEmsmt(mch,isam)=vEs*(1.-frots(ich)*float(isam))+ &
          vEe*frots(ich)*float(isam)

      enddo
    enddo  ! End of loop over channels

    return

  END SUBROUTINE calSimScan

!----------------------------------------------------------------------------

  FUNCTION getCounts(Tin,ithm)

!  Inverse of Thermistor lookup table 

    integer                        :: getCounts ! Returned count
    real,               intent(in) :: Tin       ! Temps reported by PRT
    integer,            intent(in) :: ithm      ! row index in SCF/Thermistor

    integer :: i,nNodes

    nNodes=size(thermT,2)-1 

    if ((Tin < thermT(ithm,0)) .or. (Tin > thermT(ithm,nNodes))) then
      print *,'err[CalibSim::getCounts]: Temperature out of range'
      print *,'Tin,ithm,thermT(ithm,0),thermT(ithm,nNodes):',Tin,ithm, &
               thermT(ithm,0), thermT(ithm,nNodes)
      call errorHalt(1)
    endif

    do i=1,nNodes
      if (Tin < thermT(ithm,i)) then
        getCounts = i
        return
      endif
    enddo

  END FUNCTION getCounts

!----------------------------------------------------------------------------

  SUBROUTINE getSensConst(U_SO,U_Cal,U_CalErr,U_AntV,U_AntVErr,U_Bmap,U_Thm, &
    F_SO,F_Cal,F_CalErr,F_AntV,F_AntVErr,F_Bmap,F_Thm, &
    nsencompcD,nsencompD,nsccompcD,nsccompD,&
    nthmrpsD,thermID,thermT) 

! Read sensor constants files

    use SCFread

!-- Calling arguments

    integer,          intent(in)    :: U_SO
    integer,          intent(in)    :: U_Cal,U_CalErr
    integer,          intent(in)    :: U_AntV,U_AntVErr
    integer,          intent(in)    :: U_Bmap
    integer,          intent(in)    :: U_Thm 
    character(len=*), intent(in)    :: F_SO
    character(len=*), intent(in)    :: F_Cal,F_CalErr
    character(len=*), intent(in)    :: F_AntV,F_AntVErr
    character(len=*), intent(in)    :: F_Bmap
    character(len=*), intent(in)    :: F_Thm
    integer,          intent(inout) :: nthmrpsD
    integer,          intent(inout) :: nsencompcD
    integer,          intent(inout) :: nsencompD
    integer,          intent(inout) :: nsccompcD
    integer,          intent(inout) :: nsccompD
    character(lID), dimension(:),   pointer :: thermID
    real,           dimension(:,:), pointer :: thermT ! nthmr x npnts


!-- Local variables

    integer                                 :: i,j
    integer                                 :: nbeamCE,nchanCE
    integer                                 :: nbeamAV,nbeamVE
    real          , parameter               :: sigma3=3.
    character(lID), dimension(:), pointer   :: beamIDCE,chanIDCE
    character(lID), dimension(:), pointer   :: beamIDAV,beamIDVE
    
    call getSCFcalib(U_Cal,F_Cal, &
      nbeam=nbeam,beamID=beamID,chanID=chanID,nchan=nchan, &
      ncssamp=ncssamp,nwlsamp=nwlsamp,&  
      spocs=spocs,spow=spow,emis_refcs=emis_refcs,emis_w=emis_w, &
      nsccompc=nsccompc,sccompcID=sccompcID,emis_scc=emis_scc, &
      vf_sccs=vf_sccs,vf_scw=vf_scw, &
      nsencompc=nsencompc,sencompcID=sencompcID,emis_senc=emis_senc, &
      vf_sencs=vf_sencs,vf_senw=vf_senw, &
      nthmrps=nthmrps,iBeam=iBeam)
      
    call getSCFcalErr(U_CalErr,F_CalErr, &
      nbeamCE,beamIDCE,spocsErr,spowErr,nchanCE,chanIDCE, &
      nonlin,gainVar)
    
    ! check whether spillover data may likely become negative 
    !   during error application process:
    do i=1,nbeam
      do j=1,ncssamp(i)
        if (spocs(i,j).lt.(sigma3*spocsErr(i))) then
          print*,'three standard deviations about spillover value extend',&
               'to negative values for input spillover knowledge error'
          print*,'i, j, cold spillover, spillover knowledge error: '
          print*,i, j, spocs(i,j), spocsErr(i) 
          call errorHalt(1)
        endif
      enddo
    enddo
    do i=1,nbeam
      do j=1,nwlsamp(i)
        if (spow(i,j).lt.(sigma3*spowErr(i))) then
          print*,'three standard deviations about spillover value extend',&
               'to negative values for input spillover knowledge error'
          print*,'i, j, warm spillover, spillover knowledge error: '
          print*,i, j, spow(i,j), spowErr(i) 
          call errorHalt(1)
        endif
      enddo
    enddo
    
    if (nbeamCE /= nbeam) then
      print *,'err[CalibSim::getSensConst]: Cal Error data inconsistent', &
       ' with cal data'
      print *,'nbeamCE /= nbeam',nbeamCE,nbeam
      call errorHalt(1)
    endif
    if (any(beamIDCE /= beamID)) then
      print *,'err[CalibSim::getSensConst]: Cal Error beamID inconsistent', &
       ' with cal data'
      do i=1,nbeam
        print *,i,beamIDCE(i),beamID(i)
      enddo
      call errorHalt(1)
    endif
    if (nchanCE /= nchan) then
      print *,'err[CalibSim::getSensConst]: Cal Error data inconsistent', &
       ' with cal data'
      print *,'nchanCE /= nchan',nchanCE,nchan
      call errorHalt(1)
    endif    
    if (any(chanIDCE /= chanID)) then
      print *,'err[CalibSim::getSensConst]: Cal Error chanID inconsistent', &
       ' with cal data'
      do i=1,nchan
        print *,i,chanIDCE(i),chanID(i)
      enddo
      call errorHalt(1)
    endif
    
    call getSCFantView(U_AntV,F_AntV, &
      nbeamAV,beamIDAV,nsamp,spoev,emis_ref,refMap, &
      nsccomp,sccompID,emis_sc,vf_scev, &
      nsencomp,sencompID,emis_sen,vf_senev ) 

    call getSCFantViewErr(U_AntVErr,F_AntVErr, &
    nbeamVE,beamIDVE,spoevErr)
    
    do i=1,nBeam
      if (spoev(i).lt.(sigma3*spoevErr(i))) then
        print*,'three standard deviations about spillover value extend', & 
             'to negative values for input spillover knowledge error'
        print*,'i, j, earth view spillover, spillover knowledge error: '
        print*,i, j, spoev(i), spoevErr(i) 
        call errorHalt(1)
      endif
    enddo
    
    if (nbeamAV /= nbeam .or. nbeamVE /= nbeam) then
      print *,'err[CalibSim::getSensConst]: Antenna data inconsistent', &
       ' with cal data'
      print *,'nbeamAV,nbeam',nbeamAV,nbeam
      print *,'nbeamVE,nbeam',nbeamVE,nbeam
      call errorHalt(1)
    endif
    if (any(beamIDAV /= beamID) .or. any(beamIDVE /= beamID)) then
      print *,'err[CalibSim::getSensConst]: ', &
       'Antenna beamID inconsistent with cal data'
      do i=1,nbeam
        print *,i,beamIDAV(i),beamID(i),beamIDVE(i),beamID(i)
      enddo
      call errorHalt(1)
    endif

    call getSCFscanOrbit(U_SO,F_SO,scanRate=scanRate,beamID=beamIDSO, &
         sampInt=sampInt)
    if (any(beamIDSO /= beamID)) then
      print *,'err[CalibSim::getSensConst]: ', &
       'Scan/Orbit beamID inconsistent with cal data'
      do i=1,nbeam
        print *,i,beamIDSO(i),beamID(i)
      enddo
      call errorHalt(1)
    endif

    call getSCFthermistor(U_Thm,F_Thm,nthmr,npnts,thermID,thermT)

!-- Load dummy arguments for SC export
    nthmrpsD   = nthmrps
    nsencompcD = nsencompc
    nsencompD  = nsencomp
    nsccompcD  = nsccompc
    nsccompD   = nsccomp

    return

  END SUBROUTINE getSensConst

!----------------------------------------------------------------------------

  SUBROUTINE saveOrigSC()

! Save the originals of sensor constants that can get errors applied

    allocate(spocs0(size(spocs,1),size(spocs,2)))
    allocate(spow0(size(spow,1),size(spow,2)))
    allocate(emis_refcs0(size(emis_refcs)))
    allocate(emis_w0(size(emis_w)))
    allocate(vf_sccs0(size(vf_sccs,1),size(vf_sccs,2),size(vf_sccs,3)))
    allocate(vf_scw0(size(vf_scw,1),size(vf_scw,2),size(vf_scw,3)))
    allocate(vf_sencs0(size(vf_sencs,1),size(vf_sencs,2),size(vf_sencs,3)))
    allocate(vf_senw0(size(vf_senw,1),size(vf_senw,2),size(vf_senw,3)))
    allocate(spoev0(size(spoev)))
    allocate(emis_ref0(size(emis_ref)))
    allocate(vf_scev0(size(vf_scev,1),size(vf_scev,2),size(vf_scev,3)))
    allocate(vf_senev0(size(vf_senev,1),size(vf_senev,2),size(vf_senev,3)))
    
    spocs0      = spocs     
    spow0       = spow      
    emis_refcs0 = emis_refcs
    emis_w0     = emis_w    
    vf_sccs0    = vf_sccs   
    vf_scw0     = vf_scw    
    vf_sencs0   = vf_sencs  
    vf_senw0    = vf_senw   
    spoev0      = spoev     
    emis_ref0   = emis_ref  
    vf_scev0    = vf_scev   
    vf_senev0   = vf_senev  

    return

  END SUBROUTINE saveOrigSC

!----------------------------------------------------------------------------

  SUBROUTINE errorSC(RNnonlin,U_nonlin,F_nonlin, &
    RNgainVar,U_gainVar,F_gainVar,RNspoCalE,U_spoCalE,F_spoCalE, &
    RNspoevE,U_spoevE,F_spoevE)

! Add errors to sensor constants
! The errors are specified as scale factors.  Added error is scale factor 
! times SCF-defined error limit.  Factors may be negative.

!-- Arguments
    logical,            intent(in), optional :: RNnonlin
    integer,            intent(in), optional :: U_nonlin
    character(*),       intent(in), optional :: F_nonlin
    logical,            intent(in), optional :: RNgainVar
    integer,            intent(in), optional :: U_gainVar
    character(*),       intent(in), optional :: F_gainVar
    logical,            intent(in), optional :: RNspoCalE
    integer,            intent(in), optional :: U_spoCalE
    character(*),       intent(in), optional :: F_spoCalE
    logical,            intent(in), optional :: RNspoevE
    integer,            intent(in), optional :: U_spoevE
    character(*),       intent(in), optional :: F_spoevE

!-- Local variables
    real          :: gauss   ! external library function
    real          :: nonlinF,gainVarF
    real          :: spocsF,spowF,spoevF
    integer       :: i,j,nchanIn,ndatInCS,ndatInWL,ndatInEV
    integer, save :: ix=0

    if (allocated(nonlinChan)) deallocate(nonlinChan)
    if (allocated(gainVarChan)) deallocate(gainVarChan)
    allocate (nonlinChan(nchan),gainVarChan(nchan))
    nonlinChan=0.  ! initialize
    gainVarChan=0. ! initialize

    if (present(RNnonlin)) then
      if (RNnonlin) then
        do i=1,nchan
          call RANDOM_NUMBER(nonlinF)   ! Nonlinearity assumed positive always
          nonlinChan(i)=nonlin(iBeam(i))*nonlinF
        enddo
      endif
    endif

    if (present(F_nonlin)) then
      if (.not. present(U_nonlin)) then
        print *,'err[CalibSim::errorSC]: U_nonlin must come with F_nonlin'
        call errorHalt(1)
      endif
      if (F_nonlin /= NULLF .and. F_nonlin /= BLANKF) then
        if (any(nonlinChan(:) /= 0.)) &
          print *,'msg[CalibSim::errorSC]: F_nonlin overriding RNnonlin'
        open (U_nonlin,file=F_nonlin,status='old',action='read')
        read (U_nonlin,'(a4)')nchanIn
        if (nchanIn /= nchan) then
          print *,'err[CalibSim::errorSC]: nchan mismatch F_nonlin:', &
            nchan,nchanIn
          call errorHalt(1)
        endif
        do i=1,nchan
          read(U_nonlin,'(4x,f6.3)')nonlinF
          nonlinChan(i)=nonlin(iBeam(i))*nonlinF
        enddo
        close (U_nonlin)
      endif
    endif

    if (present(RNgainVar)) then
      if (RNgainVar) then
        do i=1,nchan
          call RANDOM_NUMBER(gainVarF)
          gainVarF=gainVarF*2.-1.  ! Expand range -1 to +1
          gainVarChan(i)=gainVar(iBeam(i))*gainVarF
        enddo
      endif
    endif

    if (present(F_gainVar)) then
      if (.not. present(U_gainVar)) then
        print *,'err[CalibSim::errorSC]: U_gainVar must come with F_gainVar'
        call errorHalt(1)
      endif
      if (F_gainVar /= NULLF .and. F_gainVar /= BLANKF) then
        if (any(gainVarChan(:) /= 0.)) &
          print *,'msg[CalibSim::errorSC]: F_gainVar overriding RNgainVar'
        open (U_gainVar,file=F_gainVar,status='old',action='read')
        read (U_gainVar,'(a4)')nchanIn
        if (nchanIn /= nchan) then
          print *,'err[CalibSim::errorSC]: nchan mismatch F_gainVar:',nchan,nchanIn
          call errorHalt(1)
        endif
        do i=1,nchan
          read(U_gainVar,'(4x,f6.3)')gainVarF
          gainVarChan(i)=gainVar(iBeam(i))*gainVarF
        enddo
        close (U_gainVar)
      endif
    endif

    if (present(RNspoCalE)) then
      if (RNspoCalE) then
        do i=1,nbeam
          do j=1,ncssamp(i)
             spocs(i,j)=spocs0(i,j)+spocsErr(i)*gauss(ix,1.,0.)
             if (spocs(i,j).lt.0.) spocs(i,j)=0.
          enddo
          do j=1,nwlsamp(i)
             spow(i,j)=spow0(i,j)+spowErr(i)*gauss(ix,1.,0.)
             if (spow(i,j).lt.0) spow(i,j)=0.
          enddo
        enddo
      endif
    endif

    if (present(F_spoCalE)) then
      if (.not. present(U_spoCalE)) then
        print *,'err[CalibSim::errorSC]: U_spoCalE must come with F_spoCalE'
        call errorHalt(1)
      endif
      if (F_spoCalE /= NULLF .and. F_spoCalE /= BLANKF) then
        if (any(spocs /= spocs0) .or. any(spow /= spow0)) &
          print *,'msg[CalibSim::errorSC]: F_spoCalE overriding RNspoCalE'
        open (U_spoCalE,file=F_spoCalE,status='old',action='read')
        read (U_spoCalE,'(2i4)') ndatInCS,ndatInWL
        
        if (ndatInCS /= sum(ncssamp(:))) then
          print *,'err[CalibSim::errorSC]: ncssamp mismatch F_spoCalE:', &
               '# of cold sky cal samples, # of records in cal bias file: ', &
               sum(ncssamp(:)),ndatInCS
          call errorHalt(1)
        endif
        if (ndatInWL /= sum(nwlsamp(:))) then
          print *,'err[CalibSim::errorSC]: nwlsamp mismatch F_spoCalE:', &
               '# of warm load cal samples, # of records in cal bias file: ',&
               sum(nwlsamp(:)),ndatInWL 
          call errorHalt(1)
        endif

        do i=1,nbeam
          do j=1,ncssamp(i)
            read(U_spoCalE,'(4x,f9.6)')spocsF
            spocs(i,j)=spocs0(i,j)+spocsErr(i)*spocsF
          enddo
        enddo
        do i=1,nbeam
          do j=1,nwlsamp(i)
            read(U_spoCalE,'(4x,f9.6)')spowF
            spow(i,j) =spow0(i,j)+spowErr(i)*spowF
          enddo
        enddo
        close (U_spoCalE)
      endif
    endif

    if (present(RNspoevE)) then
      if (RNspoevE) then
        do i=1,nbeam
          spoev(i)=spoev0(i)+spoevErr(i)*gauss(ix,1.,0.)
          if (spoev(i).lt.0.) spoev(i)=0.
        enddo
      endif
    endif

    if (present(F_spoevE)) then
      if (.not. present(U_spoevE)) then
        print *,'err[CalibSim::errorSC]: U_spoevE must come with F_spoevE'
        call errorHalt(1)
      endif
      if (F_spoevE /= NULLF .and. F_spoevE /= BLANKF) then
        if (any(spoev /= spoev0)) &
          print *,'msg[CalibSim::errorSC]: F_spoevE overriding RNspoevE'
        open (U_spoevE,file=F_spoevE,status='old',action='read')
        read (U_spoevE,'(i4)')ndatInEV
        
        if (ndatInEV /= nBeam) then
          print *,'err[CalibSim::errorSC]: nBeam mismatch F_spoevE:', &
               'nBeam, # of records in earth view bias file: ', &
               nBeam, ndatInEV
        endif
        do i=1,nbeam
          read(U_spoevE,'(4x,f9.6)')spoevF
          spoev(i)=spoev0(i)+spoevErr(i)*spoevF
        enddo
        close (U_spoevE)
      endif
    endif

    return

  END SUBROUTINE errorSC

!----------------------------------------------------------------------------

END MODULE CalibSim
