MODULE PlanckTabMW 

!**********************************************************************
!* Module : PlanckTabMW 
!* Purpose: Build a Planck Function table at set MW channel and temperature
!           so that accurate (better than 0.01 K) linear interpolations, 
!           for both forward and backward (inverse) calculations of the 
!           Planck function, can be achieved with reduced computation cost. 
!*
!* Inputs/Outputs:
!* 1. List of files containing channel-specific spectral response 
!     information, conforming with CMIS SCF/Spectral Response
!  2. Lists of channels with center frequency only
!  
!
!* nprim         integer    number of primary sensor channels
!* nctr          integer    number of channels with given center freq
!* nchanP        integer    total number of channels
!* chanIDp       character  channel ID
!* nNodes        integer    number of temperature nodes 
!* tPTM          real       nNodes temperature points
!                           (last one for CosmBckg)
!* rPTM          real       nchanP x nNodes radiances
!* tPTMmin       real       minimum PTM temperature
!* tPTMmax       real       maximum PTM temperature
!* rPTMlmt       real       radiance limits in in Planck Table Module
!
!* tslope        real       nchanP x (nNodes-2) temp gradients 
!* rslope        real       nchanP x (rNodes-2) radiance gradients 
!* F_list        character  list of SRF file names
!* F_ptable      character  Planck table file name
!*
!* Copyright: 2005, AER, Inc.
!* Developed by AER, Inc., 2005
!**********************************************************************

  use constants, only: CosmBckg

  implicit none
  private
  public  :: initPTM,buildPTM,fwdPTM,invPTM,freePTM

!-- Global private data
  real,           parameter :: GHzToHz=1.e9
  integer,        parameter :: MaxLine=500
  integer,        parameter :: lID    =12

  integer                                      :: nprim
  integer                                      :: nctr
  integer                                      :: nchanP  !nprim+nctr
  integer                                      :: nNodes
  character(lID), dimension(:),   pointer      :: chanIDp
  real,           dimension(:,:), allocatable  :: rPTM   !nchanP x nNodes  
  real,           dimension(:,:), allocatable  :: tslope !nchanP x (nNodes-2)
  real,           dimension(:,:), allocatable  :: rslope !nchanP x (nNodes-2)
  real,           dimension(:),   allocatable  :: cfreq  !effective center-f
  real                                         :: tin    !single temp in 
  real,           dimension(:),   allocatable  :: tPTM   !temp in K
  real                                         :: tgdPTM
  real,           public                         :: tPTMmin, tPTMmax
  character(lID), public, dimension(:),  pointer :: chanPTMid
  real,           public, dimension(:,:),pointer :: rPTMlmt

  integer                                      :: ich,i,j,k,n
  character(len=30)                            :: cfmt

  interface fwdPTM
     module procedure fwdPTMvec
     module procedure fwdPTMscl
  end interface

CONTAINS

  SUBROUTINE buildPTM(U_list,F_list,U_ctrF,F_ctrF,U_ptable,F_ptable, &
        tPTMmin,tPTMmax,nNodes,AddCtr)
  
  use scfread
  use Quadrature, only : intTrapezoid
  use constants,  only : MISSING_REAL, LtSp_cgs

!-- Calling arguments
  integer,            intent(in)     :: U_list
  integer,            intent(in)     :: U_ptable
  character(len=200), intent(in)     :: F_list
  character(len=200), intent(in)     :: F_ptable
  character(len=200), intent(in)     :: F_ctrF
  integer,            intent(in)     :: U_ctrF
  real,               intent(in)     :: tPTMmin
  real,               intent(in)     :: tPTMmax
  integer,            intent(in)     :: nNodes
  logical,            intent(in)     :: AddCtr

!-- Local variables
  character(len=30)                  :: cfmt
  character(len=200)                 :: F_srf
  character(len=300)                 :: dumCh
  integer                            :: U_srf=50
  integer                            :: nfreq, nband
  real                               :: offFreq
  real                               :: offFreq1, offFreq2
  real, dimension(:),   pointer      :: rFreq, rMag
  real, dimension(:),   allocatable  :: rMaglin,wavenR,radR ! Planck rad
  real, dimension(:),   allocatable  :: rFreqL, rFreqR, cfreqR
  real, dimension(:),   allocatable  :: crad, cradR
  real, dimension(:,:), allocatable  :: crad2d
  character(lID), dimension(:), allocatable:: cvec1
  real                               :: waven1, rad1
  real                               :: cfreq
  integer                            :: ipol
  integer                            :: ich,ifreq
  integer                            :: i,j,k
  integer, dimension(4)              :: sign1, sign2
  real                               :: int1, int2
  real                               :: planck   ! external subprogram

  if (allocated(tPTM)) deallocate(tPTM)
     allocate(tPTM(nNodes))
  if (allocated(cradR)) deallocate(cradR,crad,crad2d)
  allocate(cradR(nNodes),crad(nNodes),crad2d(4,nNodes))

  open(U_list,file=F_list,form='formatted',status='old')
     do ich =1, MaxLine 
       read(U_list, '(a200)', END=100)
     enddo
100  nprim = ich-1
     rewind(U_list)

  if (.not. AddCtr) then
    nctr = 0
  else
    open(U_ctrF,file=F_ctrF,form='formatted',status='old')
       do ich =1, MaxLine 
         read(U_ctrF, '(a200)', END=200)
       enddo
200    nctr = ich-1
       rewind(U_ctrF)
  endif
     nchanP=nprim+nctr

  if (allocated(cvec1)) deallocate(cvec1)
     allocate(cvec1(nchanP))

  offFreq2 = 0.0
  sign1 = (/1,-1, 1,-1/)   
  sign2 = (/1, 1,-1,-1/)

  tgdPTM = (tPTMmax -tPTMmin)/(nNodes -2)  
  do j=1,nNodes-1
    tPTM(j) = tPTMmin +(j-1) *tgdPTM
  enddo
    tPTM(nNodes) = CosmBckg 

  open(U_ptable, file=F_ptable, form='formatted',status='unknown')
  write(U_ptable,'(i4, 3(4x,i4))') nprim, nctr, nchanP, nNodes
  write(cfmt,'( ''('', i2, ''f9.2)'' )') nNodes 
  write(U_ptable,cfmt) tPTM

  write(U_ptable,*) 'No. chanIDp            ECF        Channel-integ radiance'

  do ich = 1,nprim
    read(U_list,'(a200)') F_srf
    call getSCFspecResponse(U_srf,F_srf,chanID=cvec1(ich),nresp=nfreq,&
         nband=nband,offFreq=offFreq1,respFreq=rFreq,respMagnitude=rMag)

    call parseID(cvec1(ich),ipol=ipol)

    ! variables for single passband 
    if (allocated(rMaglin)) deallocate(rMaglin,wavenR,radR,rFreqR)
       allocate(rMaglin(nfreq),wavenR(nfreq),radR(nfreq),rFreqR(nfreq))
    if (allocated(cfreqR)) deallocate(cfreqR)
       allocate(cfreqR(nband))

    do k=1, nband
       rFreqR = offFreq1+sign1(k)*rFreq + sign2(k)*offFreq2

       rMaglin = 10.0**(0.1*rMag)           ! Convert from dB to linear
       call intTrapezoid(rFreqR,rFreqR*rMaglin,int1)
       call intTrapezoid(rFreqR,rMaglin,int2)
       cfreqR(k)  = int1/int2 

    !-- Convert freq (GHz) to wavenumber (cm-1), for use in Planck function 
       wavenR=rFreqR*GHzToHz/LtSp_cgs
       do j=1,nNodes
         do i=1,nfreq
           radR(i) = planck(wavenR(i), tPTM(j))
         enddo
         call intTrapezoid(rFreqR,radR*rMaglin,cradR(j))
         crad2d(k,j)  = cradR(j)/int2
       enddo
     enddo

     cfreq = sum(cfreqR(1:nband))/nband
     do j=1, nNodes
        crad(j) = sum(crad2d(1:nband,j))/nband
     enddo
    write(cfmt,'( ''(i3,1x,a12,f14.9,'', i2, ''e14.6)'' )') nNodes
    write(U_ptable,cfmt) ich,cvec1(ich),cfreq,crad
  enddo

!-- channels with center frequency only

  if (nctr /= 0) then
     do ich=nprim+1,nchanP
       read(U_ctrF, '(a12,1x,f9.6)') cvec1(ich), cfreq 
       waven1=cfreq*GHzToHz/LtSp_cgs
       do j=1,nNodes
          crad(j) = planck(waven1, tPTM(j))
       enddo
       write(U_ptable,cfmt) ich,cvec1(ich),cfreq,crad
     enddo
  endif

  deallocate(rMaglin,wavenR,radR,rFreqR,cfreqR)
  deallocate(cradR,crad2d,crad,cvec1)
  
  close(U_ctrF)
  close(U_list)
  close(U_ptable)

  END SUBROUTINE buildPTM

!------------------------------------------------------

  SUBROUTINE initPTM(iu,fname) 

! Description: Initialize and read in PTM data

!-- Callling arguments

  integer,              intent(in)          :: iu 
  character(len=200),   intent(in)          :: fname 
  
!-- Local variables
  integer                                   :: ich, i, k
  character(lID)                            :: dumID
  character(len=300)                        :: dumCh


  open(iu,file=fname,form='formatted',status='old')

  read(iu,'(i4, 3(4x,i4) )') nprim, nctr, nchanP, nNodes
  
!-- Initialization of Planck Table variables 
 
  if (allocated(rPTM)) deallocate(rPTM)
     allocate(rPTM(nchanP,nNodes)) 
  if (allocated(tslope)) deallocate(tslope)
     allocate(tslope(nchanP,nNodes-2))
  if (allocated(rslope)) deallocate(rslope)
     allocate(rslope(nchanP,nNodes-2))
  if (allocated(cfreq)) deallocate(cfreq)
     allocate(cfreq(nchanP))
  if (allocated(tPTM)) deallocate(tPTM)
     allocate(tPTM(nNodes))
  if (associated(chanIDp)) deallocate(chanIDp)
     allocate(chanIDp(nchanP))

  if (associated(chanPTMid)) deallocate(chanPTMid)
     allocate(chanPTMid(nchanP))
  if (associated(rPTMlmt)) deallocate(rPTMlmt)
     allocate(rPTMlmt(nchanP,2))

  write(cfmt,'( ''('', i2, ''f9.2)'' )') nNodes
  read(iu,cfmt) tPTM
  tgdPTM = tPTM(2) - tPTM(1)
  read(iu,'(a300)') dumCh 

  write(cfmt,'( ''(i3,1x,a12,f14.9,'', i2, ''e14.6)'' )') nNodes
  do i = 1,nchanP
     read(iu,cfmt) ich,chanIDp(i),cfreq(i),rPTM(i,1:nNodes)
     do k=1, nNodes-2
       tslope(i,k) = tgdPTM/(rPTM(i,k+1) - rPTM(i,k))
       rslope(i,k) = 1/tslope(i,k) 
     enddo
     rPTMlmt(i,1) = rPTM(i,1)
     rPTMlmt(i,2) = rPTM(i,nNodes-1)
  enddo
  chanPTMid  = chanIDp

  close(iu)
  return

  END SUBROUTINE initPTM
                                                                                                                             
!------------------------------------------------------

  FUNCTION fwdPTMvec(ch1, ch2, tin)

!-- single-channel Planck calculation

  integer,        intent(in)    :: ch1       !beginning chan num
  integer,        intent(in)    :: ch2       !ending chan num
  real,           intent(in)    :: tin       !input temp in K
  real,   dimension(ch1:ch2)    :: fwdPTMvec !returned radiance vec
  
!-- Local variables
  integer                       :: ich, k

!-- input data validation and range-checking
  
  if (abs(tin - CosmBckg) <= 0.1*tgdPTM ) then
    do ich=ch1, ch2
       fwdPTMvec(ich) = rPTM(ich, nNodes)
    enddo
    return
  endif

  k = 1 + int((tin-tPTM(1))/tgdPTM)
  if (k == nNodes-1) k = nNodes-2
                                                                                  
  if ( (k < 1) .or. (tin > tPTM(nNodes-1)) ) then
    print *, 'err[PlanckTabMW::fwdPTM]: input temperature out of range'
    print *, 'tin, tPTM(1), tPTM(nNodes-1)=', tin, tPTM(1), tPTM(nNodes-1)
    call errorHalt(1)
    
  endif
                                                                                  
  do ich=ch1, ch2
   fwdPTMvec(ich) = rPTM(ich,k) + (tin - tPTM(k)) * rslope(ich,k)
  enddo
                                                                                  
  return

  END FUNCTION fwdPTMvec

!------------------------------------------------------

  FUNCTION fwdPTMscl(chID1, tin)

!-- Callling arguments

  real                            :: fwdPTMscl !returned radiance
  character(lID), intent(in)      :: chID1     !input single chanIDp
  real,              intent(in)   :: tin       !input temp in K

!-- Local variables
  logical                         :: flag1
  integer                         :: ich
  character(lID)                  :: chanID_loc1,chanID_loc2
  real, dimension(1)              :: rad

  chanID_loc1 = trim(adjustl(chID1))

  flag1 = .false.
  do ich=1,nchanP
    chanID_loc2 = trim(adjustl(chanIDp(ich)))
    if (chanID_loc1(1:12) == chanID_loc2(1:12)) then
          flag1 = .true.
          exit
    endif
  enddo

  if (.not. flag1 ) then
    print *, 'err[PlanckTabMW::fwdPTMscl]: input channel invalid'
    call errorHalt(1)
  endif

  rad       = fwdPTMvec(ich,ich,tin)
  fwdPTMscl = rad(1)

  return

  END FUNCTION fwdPTMscl

!------------------------------------------------------
                                                                                                                             
  FUNCTION invPTM(chID1, rin) 
                                                                                                                             
!-- Inverse Planck Function Table for MW 

  real                                     :: invPTM ! returned Tb
  character(lID), intent(in)               :: chID1  ! single chanIDp 
  real,              intent(in)            :: rin    ! radiance,

!-- Local variables
  logical                                  :: flag1 
  integer                                  :: ich, k
  real                                     :: rtmp1, rtmp2
  character(lID)                           :: chanID_loc1,chanID_loc2

!-- input data validation and range-checking

  flag1 = .false.
  chanID_loc1 = trim(adjustl(chID1))

  do ich=1,nchanP
    chanID_loc2 = trim(adjustl(chanIDp(ich)))
    if (chanID_loc1(1:12) == chanID_loc2(1:12)) then
       flag1 = .true.
       exit
    endif 
  enddo

  if (.not. flag1) then
    print *, 'err[PlanckTabMW::invPTM]: input channel invalid'
    call errorHalt(1)
  endif

  rtmp1 = 0.1* rPTM(ich,nNodes)
  if (abs(rin -rPTM(ich,nNodes)) < rtmp1 ) then
      invPTM = tPTM(nNodes)
      return
  endif

  rtmp1 = rPTM(ich,1)
  rtmp2 = rPTM(ich,nNodes-1)
  if ( (rin < rtmp1) .or. (rin > rtmp2) ) then   
    print *, 'err[PlanckTabMW::invPTM]: input value out of range'
    print *, 'rin,chan,rPTM(ich,1), rPTM(ich,nNodes-1)=', &
             rin, ich, rPTM(ich,1), rPTM(ich,nNodes-1)
    call errorHalt(1)
  endif

!-- find bins that bracket the input value and interpolate
  
  do k=1, nNodes-2
    if (rPTM(ich,k+1) >= rin) then
      invPTM = tPTM(k) + (rin - rPTM(ich,k)) * tslope(ich,k)
      return
    endif
  enddo

  return

  END FUNCTION invPTM

!------------------------------------------------------
                                                                                  
  subroutine freePTM()
                                                                                  
! Description: releases space allocated for global variables
 
  if (associated(chanIDp))   deallocate(chanIDp)
  if (associated(chanPTMid)) deallocate(chanPTMid)
  if (associated(rPTMlmt))   deallocate(rPTMlmt)
  if ( allocated(rPTM)  )    deallocate(rPTM)
  if ( allocated(tslope))    deallocate(tslope)
  if ( allocated(rslope))    deallocate(rslope)
  if ( allocated(cfreq) )    deallocate(cfreq)
  if ( allocated(tPTM)  )    deallocate(tPTM)
                                                                                  
  return
                                                                                  
  end subroutine freePTM

!------------------------------------------------------

end MODULE PlanckTabMW 
