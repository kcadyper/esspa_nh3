MODULE SCFread
  
  ! Read sensor constants files
  ! Note: Subroutines read optional arguments as dummy values to avoid memory
  !       violations when argument is not present.

  implicit none
  private
  public:: getSCFscanOrbit,getSCFcalib,getSCFcalErr,getSCFnoise
  public:: getSCFantView,getSCFantViewErr,getSCFAntDerived
  public:: getSCFpointing,getSCFpointingErr,getSCFpolarization,getSCFpolError
  public:: getSCFspecResponse,setBeamMap
  public:: getSCFpassbandErr,getSCFsolarSpillover
  public:: getSCFantPattern,getSCFthermistor
  public:: parseID

  ! Public parameters
                                                                                     
  integer, parameter, public :: LEN_ID_SCF=12  ! length of ID strings
                                                                                     
  ! Global private parameters
                                                                                     
  integer, parameter         :: lID=LEN_ID_SCF

  ! Global private variables: dummy variables for buffering in data before
  ! copying to optional calling argument variables

  real                                        :: rdum1,rdum2,rdum3,rdum4
  real                                        :: rdum5,rdum6
  real,           allocatable, dimension(:)   :: rvec1,rvec2,rvec3,rvec4
  real,           allocatable, dimension(:)   :: rvec5,rvec6,rvec7,rvec8
  real,           allocatable, dimension(:)   :: rvec9
  real,           allocatable, dimension(:,:) :: rmat1
  integer                                     :: idum1,idum2
  integer,        allocatable, dimension(:)   :: ivec1,ivec2,ivec3
  character(lID)                              :: cdum1
  character(lID), allocatable, dimension(:)   :: cvec1,cvec2
  character(len=20) :: cfmt,cfmtClear=repeat(' ',20)

  interface parseID
     module procedure parseIDscal
     module procedure parseIDvec
  end interface

CONTAINS

  SUBROUTINE getSCFscanOrbit(iu,fname, &
       rotDir,scanRate,azRef,nbeam,EarthRad,scAlt,orbitIncl, &
       beamID,nsamp,sTimeEv, &
       ncssamp,sTimeCs,nwlsamp,sTimeW, &
       tInteg,sampInt,sampFrq)

    integer,                      intent(in)              :: iu
    character(len=*),             intent(in)              :: fname
    integer,                      intent(inout), optional :: rotDir
    real,                         intent(inout), optional :: scanRate
    real,                         intent(inout), optional :: azRef
    integer,                      intent(inout), optional :: nbeam
    real,                         intent(inout), optional :: EarthRad
    real,                         intent(inout), optional :: scAlt
    real,                         intent(inout), optional :: orbitIncl
    character(lID), dimension(:), pointer,       optional :: beamID
    integer, dimension(:),        pointer,       optional :: nsamp
    real,    dimension(:),        pointer,       optional :: sTimeEv
    integer, dimension(:),        pointer,       optional :: ncssamp
    real,    dimension(:),        pointer,       optional :: sTimeCs
    integer, dimension(:),        pointer,       optional :: nwlsamp
    real,    dimension(:),        pointer,       optional :: sTimeW
    real,    dimension(:),        pointer,       optional :: tInteg
    real,    dimension(:),        pointer,       optional :: sampInt
    real,    dimension(:),        pointer,       optional :: sampFrq

    !-- Local variables

    integer :: ifmtver,i
    integer :: nbeaml

    !-- Read file

    open(iu,file=fname,status='old',action='read')
    read(iu,'(i4)')ifmtver
    if (ifmtver /=1) then
       print *,'err[SCFread::getSCFscanOrbit]: Format version not supported'
       print *,'ifmtver:',ifmtver
       call errorHalt(1)
    endif
    read(iu,*)
    read(iu,*)
    read(iu,*)

    read(iu,'(i2,f7.3,f7.2,i4,f9.3,f8.3,f7.3)') &
         idum1,rdum1,rdum2,nbeaml,rdum3,rdum4,rdum5
    if (present(rotDir)) rotDir=idum1                 ! ScanOrbit::rotDir
    if (present(scanRate)) scanRate=rdum1             ! ScanOrbit::scanRate
    if (present(azRef)) azRef=rdum2                   ! ScanOrbit::azRef
    if (present(nbeam)) nbeam=nbeaml                  ! ScanOrbit::nbeam
    if (present(EarthRad)) EarthRad=rdum3             ! ScanOrbit::EarthRad
    if (present(scAlt)) scAlt=rdum4                   ! ScanOrbit::scAlt
    if (present(orbitIncl)) orbitIncl=rdum5           ! ScanOrbit::orbitIncl

    allocate(cvec1(nbeaml), &
         ivec1(nbeaml),rvec1(nbeaml), &
         ivec2(nbeaml),rvec2(nbeaml), &
         ivec3(nbeaml),rvec3(nbeaml), &
         rvec4(nbeaml),rvec5(nbeaml),rvec6(nbeaml))
    do i=1,nbeaml
       read(iu,'(a12,1x,i4,f9.6,i3,f9.6,i3,3f9.6,f6.1)') &
            cvec1(i),ivec1(i),rvec1(i),ivec2(i),rvec2(i), &
            ivec3(i),rvec3(i),rvec4(i),rvec5(i),rvec6(i)
    enddo
    if (present(beamID)) then                         ! scanOrbit::beamID
       allocate(beamID(nbeaml))
       beamID=cvec1
    endif
    if (present(nsamp)) then                          ! scanOrbit::nsamp
       allocate(nsamp(nbeaml))
       nsamp=ivec1
    endif
    if (present(sTimeEv)) then                        ! scanOrbit::sTimeEv
       allocate(sTimeEv(nbeaml))
       sTimeEv=rvec1
    endif
    if (present(ncssamp)) then                        ! scanOrbit::ncssamp
       allocate(ncssamp(nbeaml))
       ncssamp=ivec2
    endif
    if (present(sTimeCs)) then                        ! scanOrbit::sTimeCs
       allocate(sTimeCs(nbeaml))
       sTimeCs=rvec2
    endif
    if (present(nwlsamp)) then                        ! scanOrbit::nwlsamp
       allocate(nwlsamp(nbeaml))
       nwlsamp=ivec3
    endif
    if (present(sTimeW)) then                         ! scanOrbit::sTimeW
       allocate(sTimeW(nbeaml))
       sTimeW=rvec3
    endif
    if (present(tInteg)) then                         ! scanOrbit::tInteg
       allocate(tInteg(nbeaml))
       tInteg=rvec4
    endif
    if (present(sampInt)) then                        ! scanOrbit::sampInt
       allocate(sampInt(nbeaml))
       sampInt=rvec5
    endif
    if (present(sampFrq)) then                        ! scanOrbit::sampFrq
       allocate(sampFrq(nbeaml))
       sampFrq=rvec6
    endif
    deallocate(cvec1,ivec1,rvec1,ivec2,rvec2,ivec3,rvec3,rvec4,rvec5,rvec6)

    close(iu)
    return

  END SUBROUTINE getSCFscanOrbit

  !----------------------------------------------------------------------------

  SUBROUTINE getSCFcalib(iu,fname, &
       nbeam,beamID,ncssamp,nwlsamp,&
       nchan,chanID,iBeam, &
       spocs,spow,emis_refcs,emis_w, &
       nsccompc,sccompcID,emis_scc,vf_sccs,vf_scw, &
       nsencompc,sencompcID,emis_senc,vf_sencs,vf_senw, &
       tbcos,tminw,tmaxw,tmaxwv, &
       nfft,ifft,nthmrps,dtthrm)

    integer,                      intent(in)              :: iu
    character(len=*),             intent(in)              :: fname
    integer,                      intent(inout), optional :: nbeam
    integer,                      intent(inout), optional :: nchan
    character(lID), dimension(:), pointer,       optional :: chanID
    character(lID), dimension(:), pointer,       optional :: beamID
    real,    dimension(:,:),      pointer,       optional :: spocs
    real,    dimension(:,:),      pointer,       optional :: spow
    real,    dimension(:),        pointer,       optional :: emis_refcs
    real,    dimension(:),        pointer,       optional :: emis_w
    integer,                      intent(inout), optional :: nsccompc
    integer, dimension(:),        pointer,       optional :: sccompcID
    real,    dimension(:),        pointer,       optional :: emis_scc
    real,    dimension(:,:,:),    pointer,       optional :: vf_sccs
    real,    dimension(:,:,:),    pointer,       optional :: vf_scw
    integer,                      intent(inout), optional :: nsencompc
    integer, dimension(:),        pointer,       optional :: sencompcID
    real,    dimension(:),        pointer,       optional :: emis_senc
    real,    dimension(:,:,:),    pointer,       optional :: vf_sencs
    real,    dimension(:,:,:),    pointer,       optional :: vf_senw
    real,                         intent(inout), optional :: tbcos
    real,                         intent(inout), optional :: tminw
    real,                         intent(inout), optional :: tmaxw
    real,                         intent(inout), optional :: tmaxwv
    integer,                      intent(inout), optional :: nfft
    integer,                      intent(inout), optional :: ifft
    integer,                      intent(inout), optional :: nthmrps
    real,    dimension(:,:),      pointer,       optional :: dtthrm
    integer, dimension(:),        pointer,       optional :: iBeam
    integer,dimension(:),         pointer,       optional :: ncssamp,nwlsamp

    !-- Local variables

    integer :: ifmtver,i,j,maxcscalsamp,maxwlcalsamp
    integer :: nbeaml,nsccompcl,nsencompcl
    integer :: nvfsccsl,nvfscwl,nvfsencsl,nvfsenwl
    integer :: nthmrpsl,nthmrwl,nchanl
    character(lID), allocatable, dimension(:) :: beamIDl
    character(lID), allocatable, dimension(:) :: chanIDMB,beamIDMB
    integer,        allocatable, dimension(:) :: ncssampl,nwlsampl
    integer                                   :: nwlspl,nchanMB,ncsspl

    !-- Read file

    open(iu,file=fname,status='old',action='read')
    read(iu,'(i4)')ifmtver
    if (ifmtver /=1) then
       print *,'err[SCFread::getSCFcalib]: Format version not supported'
       print *,'ifmtver:',ifmtver
       call errorHalt(1)
    endif
    read(iu,*)
    read(iu,*)
    read(iu,*)

    read(iu,'(i4)')nbeaml                          ! Calib::nbeam
    if (present(nbeam)) nbeam=nbeaml

    allocate(beamIDl(nbeaml),ncssampl(nbeaml),nwlsampl(nbeaml))
    do i=1,nbeaml
       read(iu,'(a12,1x,2i4)')beamIDl(i),ncssampl(i),nwlsampl(i)
    enddo
    maxcscalsamp = MAXVAL(ncssampl)
    maxwlcalsamp = MAXVAL(nwlsampl)
    if (present(beamID)) then                      ! Calib::beamID
       allocate(beamID(nbeaml))
       beamID=beamIDl
    endif
    if (present(ncssamp)) then                     ! Calib::ncssamp
       allocate(ncssamp(nbeaml))
       ncssamp=ncssampl
    endif
    if (present(nwlsamp)) then                     ! Calib::nwlsamp
       allocate(nwlsamp(nbeaml))
       nwlsamp=nwlsampl
    endif

    read(iu,'(i4)') nchanl
    allocate(cvec1(nchanl),cvec2(nchanl))
    do i=1,nchanl
       read(iu,'(a12,1x,a12)') cvec1(i),cvec2(i)
    enddo

    if (present(iBeam) .OR. present(nchan))  nchanMB =nchanl
    if (present(iBeam) .OR. present(chanID)) then
       allocate(chanIDMB(nchanl))
       chanIDMB=cvec1
    endif
    if (present(chanID)) then
       allocate(chanID(nchanl))
       chanID=chanIDMB
    endif
    
    if (present(nchan))  nchan =nchanMB
    if (present(iBeam)) then
       allocate(beamIDMB(nchanl))
       beamIDMB = cvec2
       call setBeamMap(nBeam=nBeaml,beamID=beamIDl,nchanMB=nchanMB, &
            chanIDMB=chanIDMB,beamIDMB=beamIDMB,iBeam=iBeam)

    endif
    deallocate(cvec1,cvec2)
    
    read(iu,'(i4)') ncsspl

    allocate(cvec1(ncsspl),ivec1(ncsspl), &
         rvec1(ncsspl))                            ! Calib::spocs
    do i=1,ncsspl
       read(iu,'(a12,1x,i4,f9.6)')cvec1(i),ivec1(i),rvec1(i)
    enddo
    
    if (present(spocs)) then
       allocate(spocs(nbeaml,maxcscalsamp))
       spocs=0.  ! Assume zero unless otherwise specified
       do j=1,ncsspl
          do i=1,nbeaml
             if (beamIDl(i) == cvec1(j)) then
                if (ivec1(j) < 1 .or. ivec1(j) > ncssampl(i)) then 
                   print *,'err[SCFread::getSCFcalib]: Invalid sample# of ', &
                        'spocs data:',cvec1(j),ivec1(j)
                   call errorHalt(1)
                endif
                spocs(i,ivec1(j))=rvec1(j)
                exit
             endif
          enddo
          if (i > nbeaml) then  ! no match was found
             print *,'err[SCFread::getSCFcalib]: Invalid ID of spocs data'
             print *,'beamID,sample#:',cvec1(j),ivec1(j)
             call errorHalt(1)
          endif
       enddo
    endif
    deallocate(cvec1,ivec1,rvec1)    
    
    read(iu,'(i4)') nwlspl
    allocate(cvec1(nwlspl),ivec1(nwlspl), &
         rvec1(nwlspl))                              
    do i=1,nwlspl
       read(iu,'(a12,1x,i4,f9.6)')cvec1(i),ivec1(i),rvec1(i)
    enddo
    
    if (present(spow)) then                        ! Calib::spow
       allocate(spow(nbeaml,maxwlcalsamp))
       spow=0.  ! Assume zero unless otherwise specified
       do j=1,nwlspl
          do i=1,nbeaml
             if (beamIDl(i) == cvec1(j)) then
                if (ivec1(j) < 1 .or. ivec1(j) > ncssampl(i)) then 
                   print *,'err[SCFread::getSCFcalib]: Invalid sample# of ', &
                        'spow data:',cvec1(j),ivec1(j)
                   call errorHalt(1)
                endif
                spow(i,ivec1(j))=rvec1(j)
                exit
             endif
          enddo
          if (i > nbeaml) then  ! no match was found
             print *,'err[SCFread::getSCFcalib]: Invalid ID of spow data'
             print *,'beamID,sample#:',cvec1(j),ivec1(j)
             call errorHalt(1)
          endif
       enddo
    endif
    deallocate(cvec1,ivec1,rvec1)    
    
    allocate(rvec1(nbeaml),rvec2(nbeaml))          ! Calib::emis_refcs
    do i=1,nbeaml
       read(iu,'(13x,2f8.6)')rvec1(i),rvec2(i)
    enddo
    if (present(emis_refcs)) then
       allocate(emis_refcs(nbeaml))
       emis_refcs=rvec1
    endif
    if (present(emis_w)) then                      ! Calib::emis_w
       allocate(emis_w(nbeaml))
       emis_w=rvec2
    endif
    deallocate(rvec1,rvec2)

    read(iu,'(i4)')nsccompcl                       ! Calib::nsccompc
    if (present(nsccompc)) nsccompc=nsccompcl

    allocate(ivec1(nsccompcl))                     ! Calib::sccompcID
    cfmt=cfmtClear
    write(cfmt,'(''('',i2,''i4)'')')nsccompcl
    read(iu,cfmt)ivec1(:)
    if (present(sccompcID)) then
       allocate(sccompcID(nsccompcl))
       sccompcID=ivec1
    endif
    deallocate(ivec1)

    allocate(rvec1(nsccompcl))                     ! Calib::emis_scc
    cfmt=cfmtClear
    write(cfmt,'(''('',i2,''f6.3)'')')nsccompcl
    read(iu,cfmt)rvec1(:)
    if (present(emis_scc)) then
       allocate(emis_scc(nsccompcl))
       emis_scc=rvec1
    endif
    deallocate(rvec1)

    read(iu,'(i4)')nvfsccsl
    allocate(cvec1(nvfsccsl),ivec1(nvfsccsl), &
         rmat1(nvfsccsl,nsccompcl))                ! Calib::vf_sccs
    cfmt=cfmtClear
    write(cfmt,'(''(a12,1x,i4,'',i2,''f9.5)'')')nsccompcl
    do i=1,nvfsccsl
       read(iu,cfmt)cvec1(i),ivec1(i),rmat1(i,:)
    enddo
    if (present(vf_sccs)) then
       allocate(vf_sccs(nbeaml,maxcscalsamp,nsccompcl))
       vf_sccs=0.  ! Assume zero unless otherwise specified
       do j=1,nvfsccsl
          do i=1,nbeaml
             if (beamIDl(i) == cvec1(j)) then
                if (ivec1(j) < 1 .or. ivec1(j) > ncssampl(i)) then
                   print *,'err[SCFread::getSCFcalib]: Invalid sample# of ', &
                        'vf_sccs data:',cvec1(j),ivec1(j)
                   call errorHalt(1)
                endif
                vf_sccs(i,ivec1(j),:)=rmat1(j,:)
                exit
             endif
          enddo
          if (i > nbeaml) then  ! no match was found
             print *,'err[SCFread::getSCFcalib]: Invalid ID of vf_sccs data'
             print *,'beamID,sample#:',cvec1(j),ivec1(j)
             call errorHalt(1)
          endif
       enddo
    endif
    deallocate(cvec1,ivec1,rmat1)

    read(iu,'(i4)')nvfscwl
    allocate(cvec1(nvfscwl),ivec1(nvfscwl), &
         rmat1(nvfscwl,nsccompcl))
    do i=1,nvfscwl
       read(iu,cfmt)cvec1(i),ivec1(i),rmat1(i,:)
    enddo
    if (present(vf_scw)) then                      ! Calib::vf_scw
       allocate(vf_scw(nbeaml,maxwlcalsamp,nsccompcl))
       vf_scw=0.  ! Assume zero unless otherwise specified
       do j=1,nvfscwl
          do i=1,nbeaml
             if (beamIDl(i) == cvec1(j)) then
                if (ivec1(j) < 1 .or. ivec1(j) > nwlsampl(i)) then
                   print *,'err[SCFread::getSCFcalib]: Invalid sample# of ', &
                        'vf_scw data:',cvec1(j),ivec1(j)
                   call errorHalt(1)
                endif
                vf_scw(i,ivec1(j),:)=rmat1(j,:)
                exit
             endif
          enddo
          if (i > nbeaml) then  ! no match was found
             print *,'err[SCFread::getSCFcalib]: Invalid ID of vf_scw data'
             print *,'beamID,sample#:',cvec1(j),ivec1(j)
             call errorHalt(1)
          endif
       enddo
    endif
    deallocate(cvec1,ivec1,rmat1)

    read(iu,'(i4)')nsencompcl                      ! Calib::nsencompc
    if (present(nsencompc)) nsencompc=nsencompcl

    allocate(ivec1(nsencompcl))                    ! Calib::sencompcID
    cfmt=cfmtClear
    write(cfmt,'(''('',i2,''i4)'')')nsencompcl
    read(iu,cfmt)ivec1(:)
    if (present(sencompcID)) then
       allocate(sencompcID(nsencompcl))
       sencompcID=ivec1
    endif
    deallocate(ivec1)

    allocate(rvec1(nsencompcl))                    ! Calib::emis_senc
    cfmt=cfmtClear
    write(cfmt,'(''('',i2,''f6.3)'')')nsencompcl
    read(iu,cfmt)rvec1(:)
    if (present(emis_senc)) then
       allocate(emis_senc(nsencompcl))
       emis_senc=rvec1
    endif
    deallocate(rvec1)

    read(iu,'(i4)')nvfsencsl
    allocate(cvec1(nvfsencsl),ivec1(nvfsencsl), &
         rmat1(nvfsencsl,nsencompcl))              ! Calib::vf_sencs
    cfmt=cfmtClear
    write(cfmt,'(''(a12,1x,i4,'',i2,''f9.5)'')')nsencompcl
    do i=1,nvfsencsl
       read(iu,cfmt)cvec1(i),ivec1(i),rmat1(i,:)
    enddo
    if (present(vf_sencs)) then
       allocate(vf_sencs(nbeaml,maxcscalsamp,nsencompcl))
       vf_sencs=0.  ! Assume zero unless otherwise specified
       do j=1,nvfsencsl
          do i=1,nbeaml
             if (beamIDl(i) == cvec1(j)) then
                if (ivec1(j) < 1 .or. ivec1(j) > ncssampl(i)) then 
                   print *,'err[SCFread::getSCFcalib]: Invalid sample# of ', &
                        'vf_sencs data:',cvec1(j),ivec1(j)
                   call errorHalt(1)
                endif
                vf_sencs(i,ivec1(j),:)=rmat1(j,:)
                exit
             endif
          enddo
          if (i > nbeaml) then  ! no match was found
             print *,'err[SCFread::getSCFcalib]: Invalid ID of vf_sencs data'
             print *,'beamID,sample#:',cvec1(j),ivec1(j)
             call errorHalt(1)
          endif
       enddo
    endif
    deallocate(cvec1,ivec1,rmat1)    

    read(iu,'(i4)')nvfsenwl
    allocate(cvec1(nvfsenwl),ivec1(nvfsenwl), &
         rmat1(nvfsenwl,nsencompcl))
    do i=1,nvfsenwl
       read(iu,cfmt)cvec1(i),ivec1(i),rmat1(i,:)
    enddo
    if (present(vf_senw)) then                     ! Calib::vf_senw
       allocate(vf_senw(nbeaml,maxwlcalsamp,nsencompcl))
       vf_senw=0.  ! Assume zero unless otherwise specified
       do j=1,nvfsenwl
          do i=1,nbeaml
             if (beamIDl(i) == cvec1(j)) then
                if (ivec1(j) < 1 .or. ivec1(j) > nwlsampl(i)) then
                   print *,'err[SCFread::getSCFcalib]: Invalid sample# of ', &
                        'vf_senw data:',cvec1(j),ivec1(j)
                   call errorHalt(1)
                endif
                vf_senw(i,ivec1(j),:)=rmat1(j,:)
                exit
             endif
          enddo
          if (i > nbeaml) then  ! no match was found
             print *,'err[SCFread::getSCFcalib]: Invalid ID of vf_senw data'
             print *,'beamID,sample#:',cvec1(j),ivec1(j)
             call errorHalt(1)
          endif
       enddo
    endif
    deallocate(cvec1,ivec1,rmat1)

    read(iu,'(f8.3)')rdum1                         ! Calib::tbcos
    if (present(tbcos)) tbcos=rdum1

    read(iu,'(2f8.3)')rdum1,rdum2 
    if (present(tminw)) tminw=rdum1                ! Calib::tminw
    if (present(tmaxw)) tmaxw=rdum2                ! Calib::tmaxw

    read(iu,'(f8.4)')rdum1                         ! Calib::tmaxwv
    if (present(tmaxwv)) tmaxwv=rdum1

    read(iu,'(i4)')idum1                           ! Calib::nfft
    if (present(nfft)) nfft=idum1

    read(iu,'(i4)')idum1                           ! Calib::ifft
    if (present(ifft)) ifft=idum1

    read(iu,'(i4)')nthmrpsl                        ! Calib::nthmrps
    if (present(nthmrps)) nthmrps=nthmrpsl

    if(allocated(rmat1)) deallocate(rmat1)
    allocate(rmat1(nthmrwl,nthmrpsl))              ! Calib::dtthrm
    cfmt=cfmtClear
    write(cfmt,'(''('',i2,''f8.3)'')') nthmrpsl
    
    do i=1,nthmrwl
       read(iu,cfmt) rmat1(i,:)
    enddo
    if (present(dtthrm)) then
       allocate(dtthrm(nthmrwl,nthmrpsl))
       dtthrm=rmat1
    endif
    deallocate(rmat1)

    close(iu)
    deallocate(beamIDl)
    return

  END SUBROUTINE getSCFcalib

  !----------------------------------------------------------------------------

  SUBROUTINE getSCFcalErr(iu,fname, &
       nbeam,beamID,spocsErr,spowErr,nchan,chanID,nonlin,gainVar)

    integer,                      intent(in)              :: iu
    character(len=*),             intent(in)              :: fname
    integer,                      intent(inout), optional :: nbeam
    integer,                      intent(inout), optional :: nchan
    character(lID), dimension(:), pointer,       optional :: chanID
    character(lID), dimension(:), pointer,       optional :: beamID
    real,    dimension(:),        pointer,       optional :: spocsErr
    real,    dimension(:),        pointer,       optional :: spowErr
    real,    dimension(:),        pointer,       optional :: nonlin
    real,    dimension(:),        pointer,       optional :: gainVar

    !-- Local variables

    integer :: ifmtver,i,nbeaml,nchanl

    !-- Read file

    open(iu,file=fname,status='old',action='read')
    read(iu,'(i4)')ifmtver
    if (ifmtver /=1) then
       print *,'err[SCFread::getSCFcalibErr]: Format version not supported'
       print *,'ifmtver:',ifmtver
       call errorHalt(1)
    endif
    read(iu,*)
    read(iu,*)
    read(iu,*)

    read(iu,'(i4)')nbeaml                          ! CalErr::nbeam
    if (present(nbeam)) nbeam=nbeaml

    allocate(cvec1(nbeaml))
    allocate(rvec1(nbeaml),rvec2(nbeaml))

    do i=1,nbeaml
       read(iu,'(a12,1x,2f9.6)') &
            cvec1(i),rvec1(i),rvec2(i)
    enddo
    if (present(beamID)) then                      ! CalErr::beamID
       allocate(beamID(nbeaml))
       beamID=cvec1
    endif
    if (present(spocsErr)) then                    ! CalErr::spocsErr
       allocate(spocsErr(nbeaml))
       spocsErr=rvec1
    endif
    if (present(spowErr)) then                     ! CalErr::spowErr
       allocate(spowErr(nbeaml))
       spowErr=rvec2
    endif
    deallocate(cvec1,rvec1,rvec2)

    read(iu,'(i4)')nchanl
    if (present(nchan))nchan=nchanl                ! CalErr::nchan
    
    allocate(cvec1(nchanl),rvec1(nchanl),rvec2(nchanl))
    do i=1,nchanl
       read(iu,'(a12,1x,2f6.3)') &
            cvec1(i),rvec1(i),rvec2(i)
    enddo
    if (present(chanID)) then 
       allocate(chanID(nchanl))
       chanID=cvec1                                ! CalErr::chanID
    endif
    if (present(nonlin)) then                      ! CalErr::nonlin
       allocate(nonlin(nchanl))
       nonlin=rvec1
    endif
    if (present(gainVar)) then                     ! CalErr::gainVar
       allocate(gainVar(nchanl))
       gainVar=rvec2
    endif
    deallocate(cvec1,rvec1,rvec2)

    close(iu)
    return

  END SUBROUTINE getSCFcalErr

  !----------------------------------------------------------------------------

  SUBROUTINE getSCFnoise(iu,fname, &
       nchan,ncount,chanID,noiseFig,rfbw,convbw,sampInterv, &
       integTime,dGG,videoDT,quantNoise,calAmp)

    integer,                      intent(in)              :: iu
    character(len=*),             intent(in)              :: fname
    integer,                      intent(inout), optional :: nchan
    integer,        dimension(:), pointer,       optional :: ncount
    character(lID), dimension(:), pointer,       optional :: chanID
    real,           dimension(:), pointer,       optional :: noiseFig
    real,           dimension(:), pointer,       optional :: rfbw
    real,           dimension(:), pointer,       optional :: convbw
    real,           dimension(:), pointer,       optional :: sampInterv
    real,           dimension(:), pointer,       optional :: integTime
    real,           dimension(:), pointer,       optional :: dGG
    real,           dimension(:), pointer,       optional :: videoDT
    real,           dimension(:), pointer,       optional :: quantNoise
    real,           dimension(:), pointer,       optional :: calAmp

    !-- Local variables

    integer :: ifmtver,i,nchanl

    !-- Read file

    open(iu,file=fname,status='old',action='read')
    read(iu,'(i4)')ifmtver
    if (ifmtver /=1) then
       print *,'err[SCFread::getSCFnoise]: Format version not supported'
       print *,'ifmtver:',ifmtver
       call errorHalt(1)
    endif
    read(iu,*)  
    read(iu,*)
    read(iu,*)

    read(iu,'(i4)') nchanl                            ! Noise::nchan
    read(iu,*)
    read(iu,*)
    if (present(nchan)) nchan=nchanl

    allocate(ivec1(nchanl),cvec1(nchanl),rvec1(nchanl),rvec2(nchanl))
    allocate(rvec3(nchanl),rvec4(nchanl),rvec5(nchanl),rvec6(nchanl))
    allocate(rvec7(nchanl),rvec8(nchanl),rvec9(nchanl))
    do i=1,nchanl
       read(iu,'(i4,1x,a12,1x,f6.2,2f8.2,2f6.2,f10.7,3f6.3)') &
            ivec1(i),cvec1(i),rvec1(i),rvec2(i),rvec3(i),rvec4(i), &
            rvec5(i),rvec6(i),rvec7(i),rvec8(i),rvec9(i)
    enddo
    if (present(ncount)) then                         ! Noise::ncount
       allocate(ncount(nchanl))
       ncount=ivec1
    endif
    if (present(chanID)) then                         ! Noise::chanID
       allocate(chanID(nchanl))
       chanID=cvec1
    endif
    if (present(noiseFig)) then                       ! Noise::noiseFig
       allocate(noiseFig(nchanl))
       noiseFig=rvec1
    endif
    if (present(rfbw)) then                           ! Noise::rfbw
       allocate(rfbw(nchanl))
       rfbw=rvec2
    endif
    if (present(convbw)) then                         ! Noise::convbw
       allocate(convbw(nchanl))
       convbw=rvec3
    endif
    if (present(sampInterv)) then                     ! Noise::sampInterv
       allocate(sampInterv(nchanl))
       sampInterv=rvec4
    endif
    if (present(integTime)) then                      ! Noise::integTime
       allocate(integTime(nchanl))
       integTime=rvec5
    endif
    if (present(dGG)) then                            ! Noise::dGG
       allocate(dGG(nchanl))
       dGG=rvec6
    endif
    if (present(videoDT)) then                        ! Noise::videoDT
       allocate(videoDT(nchanl))
       videoDT=rvec7
    endif
    if (present(quantNoise)) then                     ! Noise::quantNoise
       allocate(quantNoise(nchanl))
       quantNoise=rvec8
    endif
    if (present(calAmp)) then                         ! Noise::calAmp
       allocate(calAmp(nchanl))
       calAmp=rvec9
    endif
    deallocate(cvec1,ivec1,rvec1,rvec2,rvec3,rvec4)
    deallocate(rvec5,rvec6,rvec7,rvec8,rvec9)

    close(iu)
    return

  END SUBROUTINE getSCFnoise

  !---------------------------------------------------------------------------

  SUBROUTINE getSCFantView(iu,fname, &
       nbeam,beamID,nsamp,spoev,emis_ref,refMap, &
       nsccomp,sccompID,emis_sc,vf_scev, &
       nsencomp,sencompID,emis_sen,vf_senev)

    integer,                      intent(in)              :: iu
    character(len=*),             intent(in)              :: fname
    integer,                      intent(inout), optional :: nbeam
    character(lID), dimension(:), pointer,       optional :: beamID
    integer, dimension(:),        pointer,       optional :: nsamp
    real,    dimension(:),        pointer,       optional :: spoev
    real,    dimension(:),        pointer,       optional :: emis_ref
    integer, dimension(:),        pointer,       optional :: refMap
    integer,                      intent(inout), optional :: nsccomp
    integer, dimension(:),        pointer,       optional :: sccompID
    real,    dimension(:),        pointer,       optional :: emis_sc
    real,    dimension(:,:,:),    pointer,       optional :: vf_scev
    integer,                      intent(inout), optional :: nsencomp
    integer, dimension(:),        pointer,       optional :: sencompID
    real,    dimension(:),        pointer,       optional :: emis_sen
    real,    dimension(:,:,:),    pointer,       optional :: vf_senev

    !-- Local variables

    integer :: ifmtver,i,j,mxsamp,nsampall
    integer :: nbeaml,nsccompl,nsencompl,nvfscevl,nvfsenevl
    integer,        allocatable, dimension(:) :: nsampl
    character(lID), allocatable, dimension(:) :: beamIDl

    !-- Read file

    open(iu,file=fname,status='old',action='read')
    read(iu,'(i4)')ifmtver
    if (ifmtver /=1) then
       print *,'err[SCFread::getSCFantView]: Format version not supported'
       print *,'ifmtver:',ifmtver
       call errorHalt(1)
    endif
    read(iu,*)
    read(iu,*)
    read(iu,*)

    read(iu,'(i4)')nbeaml                             ! AntView::nbeam
    if (present(nbeam)) nbeam=nbeaml

    allocate(beamIDl(nbeaml),nsampl(nbeaml))
    allocate(rvec1(nbeaml),rvec2(nbeaml),ivec1(nbeaml))
    do i=1,nbeaml
       read(iu,'(a12,1x,i4,2f9.6,i4)') &
            beamIDl(i),nsampl(i),rvec1(i),rvec2(i),ivec1(i)
    enddo
    mxsamp=maxval(nsampl(:))
    nsampall=sum(nsampl(:))
    if (present(beamID)) then                         ! AntView::beamID
       allocate(beamID(nbeaml))
       beamID=beamIDl
    endif
    if (present(nsamp)) then                          ! AntView::nsamp
       allocate(nsamp(nbeaml))
       nsamp=nsampl
    endif
    if (present(spoev)) then                          ! AntView::spoev
       allocate(spoev(nbeaml))
       spoev=rvec1
    endif
    if (present(emis_ref)) then                       ! AntView::emis_ref
       allocate(emis_ref(nbeaml))
       emis_ref=rvec2
    endif
    if (present(refMap)) then                         ! AntView::refMap
       allocate(refMap(nbeaml))
       refMap=ivec1
    endif
    deallocate(rvec1,rvec2,ivec1)

    read(iu,'(i4)')nsccompl                           ! AntView::nsccomp
    if (present(nsccomp)) nsccomp=nsccompl

    allocate(ivec1(nsccompl))                         ! AntView::sccompID
    cfmt=cfmtClear
    write(cfmt,'(''('',i2,''i4)'')')nsccompl
    read(iu,cfmt)ivec1(:)
    if (present(sccompID)) then
       allocate(sccompID(nsccompl))
       sccompID=ivec1
    endif
    deallocate(ivec1)

    allocate(rvec1(nsccompl))                         ! AntView::emis_sc
    cfmt=cfmtClear
    write(cfmt,'(''('',i2,''f6.3)'')')nsccompl
    read(iu,cfmt)rvec1(:)
    if (present(emis_sc)) then
       allocate(emis_sc(nsccompl))
       emis_sc=rvec1
    endif
    deallocate(rvec1)

    read(iu,'(i4)')nvfscevl
    allocate(cvec1(nsampall),ivec1(nsampall), &
         rmat1(nsampall,nsccompl))                    ! AntView::vf_scev
    cfmt=cfmtClear

    write(cfmt,'(''(a12,1x,i4,'',i2,''f9.6)'')')nsccompl
    do i=1,nvfscevl
       read(iu,cfmt)cvec1(i),ivec1(i),rmat1(i,:)
    enddo
    if (present(vf_scev)) then
       allocate(vf_scev(nbeaml,mxsamp,nsccompl))
       vf_scev=0.  ! Assume zero unless otherwise specified
       do j=1,nvfscevl
          do i=1,nbeaml
             if (beamIDl(i) == cvec1(j)) then
                if (ivec1(j) < 1 .or. ivec1(j) > nsampl(i)) then
                   print *,'err[SCFread::getSCFantView]: Invalid sample# of ', &
                        'vf_scev data:',cvec1(j),ivec1(j)
                   call errorHalt(1)
                endif
                vf_scev(i,ivec1(j),:)=rmat1(j,:)
                exit
             endif
          enddo
          if (i > nbeaml) then  ! no match was found
             print *,'err[SCFread::getSCFantView]: Invalid ID of vf_scev data'
             print *,'beamID,sample#:',cvec1(j),ivec1(j)
             call errorHalt(1)
          endif
       enddo
    endif
    deallocate(cvec1,ivec1,rmat1)

    read(iu,'(i4)')nsencompl                          ! AntView::nsencomp
    if (present(nsencomp)) nsencomp=nsencompl

    allocate(ivec1(nsencompl))                        ! AntView::sencompID
    cfmt=cfmtClear
    write(cfmt,'(''('',i2,''i4)'')')nsencompl
    read(iu,cfmt)ivec1(:)
    if (present(sencompID)) then
       allocate(sencompID(nsencompl))
       sencompID=ivec1
    endif
    deallocate(ivec1)

    allocate(rvec1(nsencompl))                        ! AntView::emis_sen
    cfmt=cfmtClear
    write(cfmt,'(''('',i2,''f6.3)'')')nsencompl
    read(iu,cfmt)rvec1(:)
    if (present(emis_sen)) then
       allocate(emis_sen(nsencompl))
       emis_sen=rvec1
    endif
    deallocate(rvec1)

    read(iu,'(i4)')nvfsenevl
    allocate(cvec1(nsampall),ivec1(nsampall), &
         rmat1(nsampall,nsencompl))                   ! AntView::vf_senev
    cfmt=cfmtClear
    write(cfmt,'(''(a12,1x,i4,'',i2,''f9.6)'')')nsencompl
    do i=1,nvfsenevl
       read(iu,cfmt)cvec1(i),ivec1(i),rmat1(i,:)
    enddo
    if (present(vf_senev)) then
       allocate(vf_senev(nbeaml,mxsamp,nsencompl))
       vf_senev=0.  ! Assume zero unless otherwise specified
       do j=1,nvfsenevl
          do i=1,nbeaml
             if (beamIDl(i) == cvec1(j)) then
                if (ivec1(j) < 1 .or. ivec1(j) > nsampl(i)) then
                   print *,'err[SCFread::getSCFantView]: Invalid sample# of ', &
                        'vf_senev data:',cvec1(j),ivec1(j)
                   call errorHalt(1)
                endif
                vf_senev(i,ivec1(j),:)=rmat1(j,:)
                exit
             endif
          enddo
          if (i > nbeaml) then  ! no match was found
             print *,'err[SCFread::getSCFantView]: Invalid ID of vf_senev data'
             print *,'beamID,sample#:',cvec1(j),ivec1(j)
             call errorHalt(1)
          endif
       enddo
    endif
    deallocate(cvec1,ivec1,rmat1)

    close(iu)
    deallocate(nsampl,beamIDl)
    return

  END SUBROUTINE getSCFantView

  !----------------------------------------------------------------------------

  SUBROUTINE getSCFantViewErr(iu,fname, &
       nbeam,beamID,spoevErr)

    integer,                      intent(in)              :: iu
    character(len=*),             intent(in)              :: fname
    integer,                      intent(inout), optional :: nbeam
    character(lID), dimension(:), pointer,       optional :: beamID
    real,    dimension(:),        pointer,       optional :: spoevErr

    !-- Local variables

    integer :: ifmtver,i,nbeaml

    !-- Read file

    open(iu,file=fname,status='old',action='read')
    read(iu,'(i4)')ifmtver
    if (ifmtver /=1) then
       print *,'err[SCFread::getSCFantViewErr]: Format version not supported'
       print *,'ifmtver:',ifmtver
       call errorHalt(1)
    endif
    read(iu,*)
    read(iu,*)
    read(iu,*)

    read(iu,'(i4)')nbeaml                          ! AntViewErr::nbeam
    if (present(nbeam)) nbeam=nbeaml

    allocate(cvec1(nbeaml),rvec1(nbeaml))
    do i=1,nbeaml
       read(iu,'(a12,1x,f9.6)')cvec1(i),rvec1(i)
    enddo
    if (present(beamID)) then                      ! AntViewErr::beamID
       allocate(beamID(nbeaml))
       beamID=cvec1
    endif
    if (present(spoevErr)) then                    ! AntViewErr::spoevErr
       allocate(spoevErr(nbeaml))
       spoevErr=rvec1
    endif
    deallocate(cvec1,rvec1)

    close(iu)
    return

  END SUBROUTINE getSCFantViewErr

  !---------------------------------------------------------------------------

  SUBROUTINE getSCFpointing(iu,fname, &
       s2sc_roll,s2sc_pitch,s2sc_yaw, &
       a2s_roll,a2s_pitch,a2s_yaw, &
       sc_yaw, nbeam,beamID,patt_Nd,patt_Az)

    integer,                      intent(in)              :: iu
    character(len=*),             intent(in)              :: fname
    integer,                      intent(inout), optional :: nbeam
    character(lID), dimension(:), pointer,       optional :: beamID
    real,                         intent(inout), optional :: s2sc_roll
    real,                         intent(inout), optional :: s2sc_pitch
    real,                         intent(inout), optional :: s2sc_yaw
    real,                         intent(inout), optional :: a2s_roll
    real,                         intent(inout), optional :: a2s_pitch
    real,                         intent(inout), optional :: a2s_yaw
    real,                         intent(inout), optional :: sc_yaw
    real,           dimension(:), pointer,       optional :: patt_Nd
    real,           dimension(:), pointer,       optional :: patt_Az

    !-- Local variables

    integer :: ifmtver,nbeaml,i

    !-- Read file

    open(iu,file=fname,status='old',action='read')
    read(iu,'(i4)')ifmtver

    if (ifmtver /=1) then
       print *,'err[SCFread::getSCFpointing]: Format version not supported'
       print *,'ifmtver:',ifmtver
       call errorHalt(1)
    endif
    read(iu,*)
    read(iu,*)
    read(iu,*)

    read(iu,'(3f8.4)') rdum1,rdum2,rdum3
    if (present(s2sc_roll)) s2sc_roll=rdum1          ! Pointing :: s2sc_roll
    if (present(s2sc_pitch)) s2sc_pitch=rdum2        ! Pointing :: s2sc_pitch
    if (present(s2sc_yaw)) s2sc_yaw=rdum3            ! Pointing :: s2sc_yaw

    read(iu,'(3f8.4)') rdum1,rdum2,rdum3             
    if (present(a2s_roll)) a2s_roll=rdum1            ! Pointing :: a2s_roll
    if (present(a2s_pitch)) a2s_pitch=rdum2          ! Pointing :: a2s_pitch
    if (present(a2s_yaw)) a2s_yaw=rdum3              ! Pointing :: a2s_yaw

    read(iu,'(f8.3)') rdum1
    if (present(sc_yaw)) sc_yaw=rdum1                ! Pointing :: sc_yaw

    read(iu,'(i4)') nbeaml
    if (present(nbeam)) nbeam=nbeaml                 ! Pointing :: nbeam

    allocate(cvec1(nbeaml),rvec1(nbeaml),rvec2(nbeaml))
    do i=1,nbeaml
       read(iu,'(a12,1x,2f8.3)') cvec1(i),rvec1(i),rvec2(i) 
    enddo
    if (present(beamID)) then
       allocate(beamID(nbeaml))
       beamID=cvec1                                  ! Pointing :: beamID
    endif
    if (present(patt_Nd)) then
       allocate(patt_Nd(nbeaml))
       patt_Nd=rvec1                                 ! Pointing :: patt_Nd
    endif
    if (present(patt_Az)) then
       allocate(patt_Az(nbeaml))
       patt_Az=rvec2                                 ! Pointing :: patt_Az
    endif
    deallocate(cvec1,rvec1,rvec2)

    close(iu)
    return

  END SUBROUTINE getSCFpointing

  !--------------------------------------------------------------------------

  SUBROUTINE getSCFpointingErr(iu,fname, &
       boreAlongErr,boreCrossErr,nbeam,beamID, &
       beamRegAlongErr,beamRegCrossErr,scanTimeErr, &
       SCxposErr,SCyposErr,SCzposErr, &
       attErrQ1,attErrQ2,attErrQ3,attErrQ4)

    integer,                      intent(in)              :: iu
    character(len=*),             intent(in)              :: fname
    real,                         intent(inout), optional :: boreAlongErr
    real,                         intent(inout), optional :: boreCrossErr
    integer,                      intent(inout), optional :: nbeam
    character(lID), dimension(:), pointer,       optional :: beamID
    real,           dimension(:), pointer,       optional :: beamRegAlongErr
    real,           dimension(:), pointer,       optional :: beamRegCrossErr
    real,                         intent(inout), optional :: scanTimeErr
    real,                         intent(inout), optional :: SCxposErr
    real,                         intent(inout), optional :: SCyposErr
    real,                         intent(inout), optional :: SCzposErr
    real,                         intent(inout), optional :: attErrQ1
    real,                         intent(inout), optional :: attErrQ2
    real,                         intent(inout), optional :: attErrQ3
    real,                         intent(inout), optional :: attErrQ4

    !-- Local variables

    integer :: ifmtver,nbeaml,i

    !-- Read file

    open(iu,file=fname,status='old',action='read')
    read(iu,'(i4)')ifmtver
    if (ifmtver /=1) then
       print *,'err[SCFread::getSCFpointingErr]: Format version not supported'
       print *,'ifmtver:',ifmtver
       call errorHalt(1)
    endif
    read(iu,*)
    read(iu,*)
    read(iu,*)

    read(iu,'(f7.4)') rdum1
    if (present(boreAlongErr)) boreAlongErr=rdum1  
                                              ! PointingErr ::boreAlongErr

    read(iu,'(f7.4)') rdum1
    if (present(boreCrossErr)) boreCrossErr=rdum1  
                                              ! PointingErr ::boreCrossErr

    read(iu,'(i4)') nbeaml
    if (present(nbeam)) nbeam=nbeaml          ! PointingErr :: nbeam 

    allocate(cvec1(nbeaml),rvec1(nbeaml),rvec2(nbeaml))
    do i=1,nbeaml
       read(iu,'(a12,1x,2f7.4)') cvec1(i),rvec1(i),rvec2(i)
    enddo
    if (present(beamID)) then
       allocate(beamID(nbeaml))
       beamID=cvec1                           ! PointingErr :: beamID
    endif
    if (present(beamRegAlongErr)) then
       allocate(beamRegAlongErr(nbeaml))
       beamRegAlongErr=rvec1                  ! PointingErr :: beamRegAlongErr
    endif
    if (present(beamRegCrossErr)) then
       allocate(beamRegCrossErr(nbeaml))
       beamRegCrossErr=rvec2                 
                                              ! PointingErr :: beamRegCrossErr
    endif
    deallocate(cvec1,rvec1,rvec2)

    read(iu,'(f9.6)') rdum1
    if (present(scanTimeErr)) scanTimeErr=rdum1  
                                              ! PointingErr :: scanTimeErr

    read(iu,'(3f6.2)') rdum1,rdum2,rdum3
    if (present(SCxposErr)) SCxposErr=rdum1   ! PointingErr :: SCxposErr
    if (present(SCyposErr)) SCyposErr=rdum2   ! PointingErr :: SCyposErr
    if (present(SCzposErr)) SCzposErr=rdum3   ! PointingErr :: SCzposErr

    read(iu,'(4f6.2)') rdum1,rdum2,rdum3,rdum4
    if (present(attErrQ1)) attErrQ1=rdum1     ! PointingErr :: attErrQ1
    if (present(attErrQ2)) attErrQ2=rdum2     ! PointingErr :: attErrQ2
    if (present(attErrQ3)) attErrQ3=rdum3     ! PointingErr :: attErrQ3
    if (present(attErrQ4)) attErrQ4=rdum4     ! PointingErr :: attErrQ4

    close(iu)
    return

  END SUBROUTINE getSCFpointingErr

  !--------------------------------------------------------------------------
  SUBROUTINE getSCFpolarization(iu,fname, & 
       nchan,chanID,polMatrix,cpol,ipol,bandID)

    integer,                       intent(in)              :: iu
    character(len=*),              intent(in)              :: fname
    integer,                       intent(inout), optional :: nchan
    character(lID), dimension(:),  pointer,       optional :: chanID
    real,           dimension(:,:),pointer,       optional :: polMatrix
    character(len=*),dimension(:), pointer,       optional :: cpol
    integer         ,dimension(:), pointer,       optional :: ipol
    character(len=*),dimension(:), pointer,       optional :: bandID
    
    !-- Local variables

    integer :: ifmtver,nchanl,i
    character(lID),dimension(:),allocatable :: chanIDl
    
    !-- Read file

    open(iu,file=fname,status='old',action='read')
    read(iu,'(i4)')ifmtver
    if (ifmtver /=1) then
       print *,'err[SCFread::getSCFpolarization]: Format version not supported'
       print *,'ifmtver:',ifmtver
       call errorHalt(1)
    endif
    read(iu,*)
    read(iu,*)
    read(iu,*)

    read(iu,'(i4)') nchanl
    if (present(nchan)) nchan=nchanl              ! Polarization :: nchan
    
    allocate(chanIDl(nchanl))
    allocate(rmat1(nchanl,9))
    do i=1,nchanl
       read(iu,'(a12,1x,9f12.6)') chanIDl(i),rmat1(i,:)
    enddo
    if (present(chanID)) then
       allocate(chanID(nchanl))
       chanID=chanIDl                             ! Polarization :: chanID
    endif
    if (present(polMatrix)) then
       allocate(polMatrix(nchanl,9))
       polMatrix=rmat1                            ! Polarization :: polMatrix
    endif
    deallocate(rmat1)
    
    ! Parse channel ID strings if arguments are present:
    if (present(cpol)) then 
       allocate(cpol(nchanl),cvec1(nchanl))
       call parseID(chanIDstring=chanIDl,cpol=cvec1)
       cpol=cvec1
       deallocate(cvec1)
    endif
    
    if (present(ipol)) then 
       allocate(ipol(nchanl))
       call parseID(chanIDstring=chanIDl,ipol=ipol)
    endif
    
    if (present(bandID)) then 
       allocate(bandID(nchanl),cvec1(nchanl))
       call parseID(chanIDstring=chanIDl,bandID=cvec1)
       bandID=cvec1
       deallocate(cvec1)
    endif
    
    deallocate(chanIDl)
    
    close(iu)
    return
  END SUBROUTINE getSCFpolarization

  !--------------------------------------------------------------------------

  SUBROUTINE getSCFpolError(iu,fname, & 
       nchan,praErr,orthogErr,chanID,polMatrixErr)

    integer,                      intent(in)              :: iu
    character(len=*),             intent(in)              :: fname
    integer,                      intent(inout), optional :: nchan
    real,                         intent(inout), optional :: praErr
    real,                         intent(inout), optional :: orthogErr
    character(lID), dimension(:), pointer,       optional :: chanID
    real,           dimension(:,:),pointer,      optional :: polMatrixErr

    !-- Local variables

    integer :: ifmtver,nchanl,i

    !-- Read file

    open(iu,file=fname,status='old',action='read')
    read(iu,'(i4)')ifmtver
    if (ifmtver /=1) then
       print *,'err[SCFread::getSCFpolError]: Format version not supported'
       print *,'ifmtver:',ifmtver
       call errorHalt(1)
    endif
    read(iu,*)
    read(iu,*)
    read(iu,*)

    read(iu,'(i4)') nchanl
    if (present(nchan)) nchan=nchanl               ! PolError :: nchan

    read(iu,'(e14.6)') rdum1
    if (present(praErr)) praErr=rdum1              ! PolError :: praErr

    read(iu,'(f7.4)') rdum1
    if (present(orthogErr)) orthogErr=rdum1        ! PolError :: orthogErr

    allocate(cvec1(nchanl),rmat1(nchanl,9))
    do i=1,nchanl
       read(iu,'(a12,1x,9f12.6)') cvec1(i),rmat1(i,:)
    enddo
    if (present(chanID)) then
       allocate(chanID(nchanl)) 
       chanID=cvec1                                ! PolError :: chanID
    endif
    if (present(polMatrixErr)) then
       allocate(polMatrixErr(nchanl,9))
       polMatrixErr=rmat1                          ! PolError :: polMatrixErr
    endif
    deallocate(cvec1,rmat1)

    close(iu)
    return
  END SUBROUTINE getSCFpolError

  !--------------------------------------------------------------------------
  SUBROUTINE getSCFspecResponse(iu,fname, & 
       chanID,nband,offFreq,ntemp,itnom,nresp,nrespVec,temps,respFreq,respFreq2D, &
       respMagnitude,respMagnitude2D,cpol,ipol,bandID)
    
    USE constants, only : MISSING_REAL
    
    integer,                      intent(in)               :: iu
    character(len=*),             intent(in)               :: fname
    character(lID),               intent(inout), optional  :: chanID
    real,                         intent(inout), optional  :: offFreq
    integer,                      intent(inout), optional  :: nband
    integer,                      intent(inout), optional  :: ntemp
    integer,                      intent(inout), optional  :: itnom
    integer,       dimension(:),  pointer,       optional  :: nrespVec
    integer,                                     optional  :: nresp
    real,          dimension(:),  pointer,       optional  :: temps
    real,          dimension(:),  pointer,       optional  :: respFreq
    real,          dimension(:,:),pointer,       optional  :: respFreq2D
    real,          dimension(:),  pointer,       optional  :: respMagnitude
    real,          dimension(:,:),pointer,       optional  :: respMagnitude2D
    character(len=*),             intent(inout), optional  :: cpol
    integer         ,             intent(inout), optional  :: ipol
    character(len=*),             intent(inout), optional  :: bandID
    
    !-- Local variables
    
    integer                            :: ifmtver,i,j,nbandl,itnoml,ntempl
    integer                            :: nrows,pt1,pt2
    integer,dimension(:),  allocatable :: nrespl
    character(lID)                     :: chanIDl
    real                               :: offFreql
    
    !-- Read file
    open(iu,file=fname,status='old',action='read')
    read(iu,'(i4)')ifmtver
    if (ifmtver /=1) then
       print *,'err[SCFread::getSCFspecResponse]: Format version not supported'
       print *,'ifmtver:',ifmtver
       call errorHalt(1)
    endif
    read(iu,*)
    read(iu,*)
    read(iu,*)
    
    read(iu,'(a12,1x,i2,f14.9,i4,i4)') chanIDl,nbandl,offFreql,ntempl,itnoml
    if (present(chanID))  chanID  = chanIDl  ! specResponse :: chanID
    if (present(nband))   nband   = nbandl   ! specResponse :: nband
    if (present(offFreq)) offFreq = offFreql ! specResponse :: offFreq
    if (present(ntemp))   ntemp   = ntempl   ! specResponse :: ntemp
    if (present(itnom))   itnom   = itnoml   ! specResponse :: itnom
    
    allocate(nrespl(ntempl))
    read(iu,'(i5)')(nrespl(i),i=1,ntempl)

    if (present(nresp)) then
!       allocate(nresp(ntempl))
       nresp=nrespl(itnoml)                      ! specResponse :: nresp
    endif

    if (present(nrespVec)) then
       allocate(nrespVec(ntempl))                ! specResponse :: nresp
       nrespVec=nrespl
    endif
    
    allocate(rvec1(ntempl))
    read(iu,'(f7.2)')(rvec1(i),i=1,ntempl)
    if (present(temps)) then
       allocate(temps(ntempl))
       temps=rvec1                           ! specResponse :: temps
    endif
    
    nrows = sum(nrespl(1:ntempl))
    allocate(rvec2(nrows),rvec3(nrows))
    read(iu,'(f14.9,f8.2)')(rvec2(j),rvec3(j),j=1,nrows)
    
    ! Add offset frequency for channels that are not double sideband
    !  (for some such channels, data is given as IF):
    ! if (nbandl.eq.1) rvec2 = rvec2 + offFreql
    
    if (present(respFreq) .or. present(respMagnitude)) then
       if (itnoml.eq.1) then
          pt1 = 1
       else
          pt1 = sum(nrespl(1:(itnoml-1))) + 1
       endif
       pt2 = pt1 + nrespl(itnoml) - 1
    endif
    if (present(respFreq)) then
       allocate(respFreq(nrespl(itnoml)))
       respFreq = rvec2(pt1:pt2)             ! specResponse :: respFreq
    endif
    if (present(respMagnitude)) then
       allocate(respMagnitude(nrespl(itnoml)))
       respMagnitude = rvec3(pt1:pt2)        ! specResponse :: respMagnitude
    endif
    
    if (present(respFreq2D)) then 
       allocate(respFreq2D(MAXVAL(nrespl),ntempl))
       respFreq2D = MISSING_REAL
       pt1=1
       pt2=nrespl(1)
       do i=1,ntempl                         ! specResponse :: respFreq2D
          respFreq2D(1:nrespl(i),i) = rvec2(pt1:pt2)
          pt1=pt1+nrespl(i)
          if (i.lt.ntempl) pt2=pt2+nrespl(i+1)
       enddo
    endif
    
    if (present(respMagnitude2D)) then
       allocate(respMagnitude2D(MAXVAL(nrespl),ntempl))
       respMagnitude2D = MISSING_REAL
       pt1=1
       pt2=nrespl(1)
       do i=1,ntempl                         ! specResponse :: respMagnitude2D
          respMagnitude2D(1:nrespl(i),i) = rvec3(pt1:pt2)
          pt1=pt1+nrespl(i)
          if (i.lt.ntempl) pt2=pt2+nrespl(i+1)
       enddo
    endif
    
    ! Parse channel ID strings if arguments are present:
    if (present(cpol)) then 
       call parseID(chanIDstring=chanIDl,cpol=cdum1)
       cpol=cdum1                            ! specResponse :: cpol
    endif
    
    if (present(ipol)) then 
       call parseID(chanIDstring=chanIDl,ipol=ipol)
    endif                                    ! specResponse :: ipol
    
    if (present(bandID)) then 
       call parseID(chanIDstring=chanIDl,bandID=cdum1)
       bandID=cdum1                          ! specResponse :: bandID
    endif
    
    deallocate(rvec1,rvec2,rvec3,nrespl)
    
  END SUBROUTINE getSCFspecResponse
  
  !--------------------------------------------------------------------------

  SUBROUTINE getSCFpassbandErr(iu,fname,nchan,chanID,rippleErr,cFrqStabil, &
       bWidthStabil)

    integer,                      intent(in)              :: iu
    character(len=*),             intent(in)              :: fname
    integer,                      intent(inout), optional :: nchan
    character(lID),dimension(:),  pointer,       optional :: chanID
    real,          dimension(:),  pointer,      optional  :: rippleErr
    real,          dimension(:),  pointer,      optional  :: cFrqStabil
    real,          dimension(:),  pointer,      optional  :: bWidthStabil

    !-- Local variables

    integer :: ifmtver,nchanl,i

    !-- Read file

    open(iu,file=fname,status='old',action='read')
    read(iu,'(i4)')ifmtver
    if (ifmtver /=1) then
       print *,'err[SCFread::getSCFpassbandErr]: Format version not supported'
       print *,'ifmtver:',ifmtver
       call errorHalt(1)
    endif
    read(iu,*)
    read(iu,*)
    read(iu,*)

    read(iu,'(i4)') nchanl
    if (present(nchan)) nchan=nchanl         ! passbandErr :: nchan

    allocate(cvec1(nchanl),rvec1(nchanl),rvec2(nchanl),rvec3(nchanl))
!!! Note:  SCF definition file gives format code as
!!!         i4,f5.2,2f12.6,2f7.2
!!!        The i4/a12 difference is understood, but not the 
!!!         '2' appearing before the final two format statements.
    allocate(cvec1(nchanl),rvec1(nchanl),rvec2(nchanl),rvec3(nchanl))
    do i=1,nchanl
       read(iu,'(a12,1x,f5.2,f12.6,f7.2)') cvec1(i),rvec1(i),rvec2(i),rvec3(i)
    enddo

    if (present(chanID)) then
       allocate(chanID(nchanl))
       chanID=cvec1                          ! passbandErr :: chanID    
    endif
    if (present(rippleErr)) then
       allocate(rippleErr(nchanl))
       rippleErr=rvec1                       ! passbandErr :: rippleErr
    endif
    if (present(cfrqStabil)) then
       allocate(cfrqStabil(nchanl))
       cfrqStabil=rvec2                      ! passbandErr :: cfrqStabil
    endif
    if (present(bWidthStabil)) then
       allocate(bWidthStabil(nchanl))
       bWidthStabil=rvec3                    ! passbandErr :: bWidthStabil
    endif

    deallocate(cvec1,rvec1,rvec2,rvec3)

    close(iu)
    return

  END SUBROUTINE getSCFpassbandErr

  !--------------------------------------------------------------------------

  SUBROUTINE getSCFsolarSpillover(iu,fname,nchSSB,nZN,nAZ,idSSB,zenCoord, &
       azimCoord,solarSpilloverTb)

    integer,                      intent(in)              :: iu
    character(len=*),             intent(in)              :: fname
    integer,                      intent(inout),optional  :: nchSSB,nZN,nAZ
    character(lID),dimension(:),  pointer,      optional  :: idSSB
    real,          dimension(:),  pointer,      optional  :: zenCoord
    real,          dimension(:),  pointer,      optional  :: azimCoord
    real,          dimension(:),  pointer,      optional  :: solarSpilloverTb

    !-- Local variables

    integer :: ifmtver,nchSSBl,nZNl,nAZl,i,j,k,ct

    !-- Read file

    open(iu,file=fname,status='old',action='read')
    read(iu,'(i4)')ifmtver
    if (ifmtver /=1) then
       print *,'err[SCFread::getSCFsolarSpillover]: Format version not supported'
       print *,'ifmtver:',ifmtver
       call errorHalt(1)
    endif
    read(iu,*)
    read(iu,*)
    read(iu,*)

    read(iu,'(i4)') nchSSBl
    if (present(nchSSB)) nchSSB=nchSSBl   ! solarSpillover :: nchSSB

    read(iu,'(i4)') nZNl
    if (present(nZN)) nZN=nZNl            ! solarSpillover :: nZN

    read(iu,'(i4)') nAZl
    if (present(nAZ)) nAZ=nAZl            ! solarSpillover :: nAZ

    if (present(idSSB)) allocate(idSSB(nchSSBl))
    if (present(zenCoord)) allocate(zenCoord(nZNl))
    if (present(azimCoord)) allocate(azimCoord(nAZl))
    if (present(solarSpilloverTb)) allocate(solarSpilloverTb(nchSSBl*nZNl*nAZl))
    
    ct=1

    if (present(idSSB) .OR. present(zenCoord) .OR. present(azimCoord) &
         .OR. present(solarSpilloverTb)) then
       do i=1,nchSSBl
          do j=1,nZNl
             do k=1,nAZl
                read(iu,'(a12,3f8.3)') cdum1,rdum1,rdum2,rdum3
                if (j.EQ.1 .AND. k.EQ.1) then
                   if (present(idSSB)) idSSB(i)=cdum1 
                                          ! solarSpillover :: idSSB
                endif
                if (i.EQ.1 .AND. k.EQ.1) then
                   if (present(zenCoord)) zenCoord(j)=rdum1
                                          ! solarSpillover :: zenCoord
                endif
                if (i.EQ.1 .AND. j.EQ.1) then
                   if (present(azimCoord)) azimCoord(k)=rdum2
                                          ! solarSpillover :: azimCoord
                endif
                if (present(solarSpilloverTb)) then
                   solarSpilloverTb(ct)=rdum3      
                                          ! solarSpillover :: solarSpilloverTb
                   ct=ct+1
                endif
             enddo
          enddo
       enddo
    endif
    
    close(iu)
    return

  END SUBROUTINE getSCFsolarSpillover

  !--------------------------------------------------------------------------

  SUBROUTINE getSCFantPattern(iu,fname,pattFrq,nvv,nhh,powRad, &
       elev,azim,coPolGain,coPolPhase,xPolGain,xPolPhase)

    integer,                      intent(in)              :: iu
    character(len=*),             intent(in)              :: fname
    real,                         intent(inout), optional :: pattFrq
    integer,                      intent(inout), optional :: nvv,nhh
    real,                         intent(inout), optional :: powRad
    real,          dimension(:),  pointer,      optional  :: elev,azim
    real,          dimension(:),  pointer,      optional  :: coPolGain
    real,          dimension(:),  pointer,      optional  :: coPolPhase
    real,          dimension(:),  pointer,      optional  :: xPolGain
    real,          dimension(:),  pointer,      optional  :: xPolPhase

    !-- Local variables

    integer   :: ifmtver,i,j,nvl,nhl,ct
    character(LEN=2)  :: str1
    character(LEN=4)  :: str2,str3
    character(LEN=8)  :: str4
    character(LEN=15) :: str5
    character(LEN=32) :: str6

    !-- Read file

    open(iu,file=fname,status='old',action='read')
    read(iu,'(i4)')ifmtver
    if (ifmtver /=1) then
       print *,'err[SCFread::getSCFantPattern]: Format version not supported'
       print *,'ifmtver:',ifmtver
       call errorHalt(1)
    endif
    read(iu,*)
    read(iu,*)
    read(iu,*)

    read(iu,'(a2,f8.3,2(2x,a4,i4),a15,f9.5)') &
         str1,rdum1,str2,nvl,str3,nhl,str5,rdum2

    if (present(pattFrq)) pattFrq=rdum1        ! SCFantPattern :: pattFrq
    if (present(nvv)) nvv=nvl                  ! SCFantPattern :: nvv
    if (present(nhh)) nhh=nhl                  ! SCFantPattern :: nhh
    if (present(powRad)) powRad=rdum2          ! SCFantPattern :: powRad

    read(iu,'(a32)') str6
    read(iu,'(6a8)') str4,str4,str4,str4,str4,str4

    if (present(elev)) then
       allocate(elev(nvl))
       elev=0.
    endif
    if (present(azim)) then
       allocate(azim(nhl))
       azim=0.
    endif
    if (present(coPolGain)) then
       allocate(coPolGain(nvl*nhl))
       coPolGain=0.
    endif
    if (present(coPolPhase)) then
       allocate(coPolPhase(nvl*nhl))
       coPolPhase=0.
    endif
    if (present(xPolGain)) then
       allocate(xPolGain(nvl*nhl))
       xPolGain=0.
    endif
    if (present(xPolPhase)) then
       allocate(xPolPhase(nvl*nhl))
       xPolPhase=0.
    endif

    if (present(elev) .OR. present(azim) .OR. present(coPolGain) &
         .OR. present(coPolPhase) .OR. present(xPolGain) &
         .OR. present(xPolPhase)) then

       ct = 1
       do i=1,nvl
          do j=1,nhl
             read(iu,'(6f8.3)') rdum1,rdum2,rdum3,rdum4,rdum5,rdum6
             if (j.EQ.1) then
                if (present(elev)) elev(i)=rdum1 
                                             ! SCFantPattern::elev
             endif
             if (i.EQ.1) then
                if (present(azim)) azim(j)=rdum2 
                                             ! SCFantPattern::azim
             endif
             if (present(coPolGain)) coPolGain(ct)=rdum3
                                             ! SCFantPattern::coPolGain
             if (present(coPolPhase)) coPolPhase(ct)=rdum4
                                             ! SCFantPattern::coPolPhase
             if (present(xPolGain)) xPolGain(ct)=rdum5
                                             ! SCFantPattern::xPolGain
             if (present(xPolPhase)) xPolPhase(ct)=rdum6
                                             ! SCFantPattern::xPolPhase
             ct = ct + 1
          enddo
       enddo

    endif

    close(iu)
    return

  END SUBROUTINE getSCFantPattern
  !--------------------------------------------------------------------------
  SUBROUTINE getSCFthermistor(iu,fname,nthmr,npnts,thermID,thermT)
    
    integer,                      intent(in)      :: iu
    character(len=*),             intent(in)      :: fname
    integer,                      intent(inout)   :: nthmr
    integer,                      intent(inout)   :: npnts
    character(lID),dimension(:),  pointer         :: thermID
    real,          dimension(:,:),pointer         :: thermT
    
    !-- Local variables
    integer, parameter :: tPerLine=16
    integer   :: ifmtver,npntsl,nrect,pt1,pt2,i,j
    
    !-- Read file
    
    open(iu,file=fname,status='old',action='read')
    read(iu,'(i4)')ifmtver
    if (ifmtver /=1) then
       print *,'err[SCFread::getSCFthermistor]: Format version not supported'
       print *,'ifmtver:',ifmtver
       call errorHalt(1)
    endif
    read(iu,*)
    read(iu,*)
    read(iu,*)
    
    read(iu,'(i4)')nthmr                     ! SCFthermistor::nthmr
    
    read(iu,'(i4)')npnts                     ! SCFthermistor::npnts
    
    allocate(thermID(nthmr),thermT(nthmr,0:(npnts-1)))
    nrect = ceiling(npnts/real(tPerLine))
    
    cfmt=cfmtClear
    write(cfmt,'(''('',i2,''f7.2)'')') tPerLine
    
    do i=1,nthmr
       read(iu,'(a12)') thermID(i)           ! SCFthermistor::thermID
       pt1 = 0
       pt2 = tPerLine-1
       do j=1,nrect
          read(iu,cfmt) thermT(i,pt1:pt2)    ! SCFthermistor::thermT
          pt1 = pt1+tPerLine
          pt2 = pt2+tPerLine
       enddo
    enddo
    
    close(iu)
    return
    
  END SUBROUTINE getSCFthermistor

  !--------------------------------------------------------------------------

  SUBROUTINE setBeamMap(nbeam,beamID,nchanMB,chanIDMB,beamIDMB, &
       nchan,chanID,iBeam)

    ! Make a map that tells, for each channel, the index of the corresponding beam 

    !-- Calling arguments

    integer,                      intent(in)           :: nbeam
    character(lID), dimension(:), intent(in)           :: beamID
    integer,                      intent(in), optional :: nchan
    character(lID), dimension(:), intent(in), optional :: chanID
    integer,        dimension(:), pointer              :: iBeam
    integer,                      intent(in)           :: nchanMB
    character(lID), dimension(:), intent(in)           :: chanIDMB
    character(lID), dimension(:), intent(in)           :: beamIDMB

    !-- Local variables
    integer                                   :: i,j

    !-- Read the map, in terms of channel and beam IDs

    if (present(nchan)) then
       if (nchanMB /= nchan) then
          print *,'err[SCFread::setBeamMap]: Beam Map data inconsistent', &
               ' with prior data'
          print *,'nchanMB /= nchan',nchanMB,nchan
          call errorHalt(1)
       endif
    endif

    if (present(chanID)) then
       if (any(chanIDMB /= chanID)) then
          print *,'err[SCFread::setBeamMap]: Beam Map chanID inconsistent', &
               ' with prior data'
          do i=1,nchanMB
             print *,i,'>',chanIDMB(i),'<>',chanID(i),'<'
          enddo
          call errorHalt(1)
       endif
    endif

    !-- Make the integer map
    allocate(iBeam(nchanMB))
    iBeam=0
    do i=1,nchanMB
       do j=1,nbeam
          if (beamIDMB(i) == ADJUSTL(beamID(j))) then
             iBeam(i)=j
             exit
          endif
       enddo
       if (iBeam(i) == 0) then
          print *,'err[SCFread::setBeamMap]: Could not find beam ID in SCF data'
          print *,'                           that matches the map beam ID'
          print *,'i,chanIDMB(i),beamIDMB(i):',i,'  ',chanIDMB(i),'  ', &
               beamIDMB(i)
          call errorHalt(1)
       endif
    enddo

    return

  END SUBROUTINE setBeamMap

  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SUBROUTINE parseIDvec(chanIDstring,cpol,ipol,bandID,bandIDlong)
    
    ! Parse channel ID into character polarization ID, integer
    !    polarization ID and band ID tags.
    ! This subroutine is for vector input and output, and it loops
    !    over parseIDscal.
    
    IMPLICIT none
    
    INTEGER :: i
    INTEGER,PARAMETER  :: mxIDlen=30,npol=6
    CHARACTER(*),dimension(:),intent(IN)             :: chanIDstring
    CHARACTER(*),dimension(:),intent(inout),optional :: cpol    
    INTEGER, dimension(:),intent(inout),optional     :: ipol
    CHARACTER(*),dimension(:),intent(inout),optional :: bandID,bandIDlong
    character(lID) :: locCpol,locBandID,locBandIDlong
    INTEGER :: vecLen,vecLenTest,locIpol
    
    vecLen = SIZE(chanIDstring)
    
    IF(present(cpol)) THEN
       vecLenTest = size(cpol)
       IF (vecLenTest.lt.vecLen) THEN
          print*,'err[SCFread::parseIDvec]: '
          print*,'size mismatch of chanIDstring and cpol: ',vecLen,vecLenTest
          call errorHalt(1)
       ENDIF
    ENDIF
    IF (present(ipol)) THEN
       vecLenTest = size(ipol)
       IF (vecLenTest.lt.vecLen) THEN
          print*,'err[SCFread::parseIDvec]: '
          print*,'size mismatch of chanIDstring and ipol: ',vecLen,vecLenTest
          call errorHalt(1)
       ENDIF
    ENDIF       
    IF (present(bandID)) THEN
       vecLenTest = size(bandID)
       IF (vecLenTest.lt.vecLen) THEN
          print*,'err[SCFread::parseIDvec]: '
          print*,'size mismatch of chanIDstring and bandID: ',vecLen,vecLenTest
          call errorHalt(1)
       ENDIF
    ENDIF
    IF (present(bandIDlong)) THEN
       vecLenTest = size(bandIDlong)
       IF (vecLenTest.lt.vecLen) THEN
          print*,'err[SCFread::parseIDvec]: '
          print*,'size mismatch of chanIDstring and bandIDlong: ', &
               vecLen,vecLenTest
          call errorHalt(1)
       ENDIF
    ENDIF
    
    DO i=1,vecLen
       call parseIDscal(chanIDstring(i),locCpol,locIpol,locBandID,&
            locBandIDlong)
       IF (present(cpol))   cpol(i)   = locCpol
       IF (present(ipol))   ipol(i)   = locIpol
       IF (present(bandID)) bandID(i) = locBandID
       IF (present(bandIDlong)) bandIDlong(i) = locBandIDlong
    ENDDO
    
    return
    
  END SUBROUTINE parseIDvec
  
  SUBROUTINE parseIDscal(chanIDstring,cpol,ipol,bandID,bandIDlong)
    
    ! Parse channel ID into character polarization ID, integer
    !    polarization ID and band ID tags.
    ! This subroutine is for scalar input and output.
    
    IMPLICIT none
    
    CHARACTER(*),    INTENT(in)              :: chanIDstring
    CHARACTER(LEN=1),INTENT(inout), optional :: cpol
    INTEGER,         INTENT(inout), optional :: ipol
    CHARACTER(*),    INTENT(inout), optional :: bandID,bandIDlong
    
    INTEGER, PARAMETER  :: npol=8
    INTEGER             :: i,idxn,idxd,idxe
    INTEGER             :: IDlen,ipolLoc
    
    CHARACTER(LEN=1)    :: cpolLoc
    CHARACTER(npol)     :: cpolAll = 'VHPMRLUF'
    
    ! Find the end of the numeric part of channel ID:
    IDlen = LEN_TRIM(chanIDstring)
    DO i=1,IDlen
       IF (chanIDstring(i:i) < '0' .OR. chanIDstring(i:i) > '9') EXIT
    ENDDO
    idxn=i-1
    IF (idxn < 1) THEN
       print *,'err[SCFread::parseIDscal]: ', &
         'First character of chanID was not numeric.'
       print *,'chanIDstring: "',chanIDstring,'"'
       call errorHalt(1)
    ENDIF

    cpolLoc=chanIDstring(idxn+1:idxn+1)

    IF (present(cpol)) cpol = cpolLoc
             
    IF (present(ipol)) THEN
       ipolLoc=INDEX(cpolAll,cpolLoc)-1
       IF (ipolLoc < 0) THEN
          print *,'err[SCFread::parseIDscal]: ', &
            'Polarization character in channel ID not recognized'
          print *,'chanIDstring: "',chanIDstring,'", cpol: "',cpolLoc,'"'
          call errorHalt(1)
       ENDIF
       ipol = ipolLoc
    ENDIF

    IF (present(bandID)) THEN

    ! Find the end of the band part of channel ID:
      idxd=INDEX(chanIDstring,'-')
      IF (idxd > 0) THEN
        idxe=idxd-1
      ELSE
        idxe=IDlen
      ENDIF

      IF (LEN(bandID) < (idxe-1)) THEN
         print *,'err[SCFread::parseIDscal]: ', &
           'bandID calling argument was not long enough to hold result'
         print *,'LEN(bandID),actual:',LEN(bandID),idxe-1
         call errorHalt(1)
      ENDIF
      bandID=chanIDstring(1:idxn)//chanIDstring(idxn+2:idxe)

    ENDIF
    
    IF (present(bandIDlong)) THEN
       bandIDlong=chanIDstring(1:idxn)// &
            chanIDstring(idxn+2:LEN_TRIM(chanIDstring))
    ENDIF
    
    return
    
  END SUBROUTINE parseIDscal

  !--------------------------------------------------------------------------
  
  SUBROUTINE getSCFAntDerived(iu,fname,chanIDIn, &
       design,bore2losAz,losNadir,fscaleinc,aview_spec, &
       ascan_spec,nedt_spec,aview_nom,ascan_nom)

    integer,                      intent(in)              :: iu
    character(len=*),             intent(in)              :: fname
    character(lID),               intent(in)              :: chanIDIn
    character(len=20),            intent(inout), optional :: design
    real,                         intent(inout), optional :: bore2losAz
    real,                         intent(inout), optional :: losNadir
    real,                         intent(inout), optional :: fscaleinc
    real,                         intent(inout), optional :: aview_spec
    real,                         intent(inout), optional :: ascan_spec
    real,                         intent(inout), optional :: nedt_spec
    real,                         intent(inout), optional :: aview_nom
    real,                         intent(inout), optional :: ascan_nom

    !-- Local variables

    integer :: ifmtver,i
    integer :: nchanDerived,iChan
    character(lID), allocatable, dimension(:) :: chanIDderived
    character(20), dimension(1)   :: adum
    !character(20) :: adum

    !-- Read file

    open(iu,file=fname,status='old',action='read')
    read(iu,'(i4)')ifmtver
    if (ifmtver /=2) then
       print *,'err[SCFread::getSCFAntDerived]: Format version not supported'
       print *,'ifmtver:',ifmtver
       call errorHalt(1)
    endif
    read(iu,*)
    read(iu,*)
    read(iu,*)

    read(iu,'(a20,i3)') adum(1),nChanDerived

    if (present(design)) design=adum(1)                 !AntDerived::design
    
    read(iu,*)

    allocate(chanIDderived(nChanDerived), &
         rvec1(nChanDerived),rvec2(nChanDerived), &
         rvec3(nChanDerived),rvec4(nChanDerived), &
         rvec5(nChanDerived),rvec6(nChanDerived), &
         rvec7(nChanDerived),rvec8(nChanDerived))

    do i=1,nChanDerived
       read(iu,'(a12,1x,e12.5,7e12.5)') &
            chanIDderived(i),rvec1(i),rvec2(i),rvec3(i),rvec4(i), &
            rvec5(i),rvec6(i),rvec7(i),rvec8(i)
    enddo

    ! Find chanIDDerived matching chanIDIn
    iChan=0
    do i=1,nChanDerived
       if (chanIDDerived(i) == ADJUSTL(chanIDIn)) then
          iChan = i
          exit
       endif
    enddo
    
    if (present(bore2losAz)) bore2losAz = rvec1(iChan) !AntDerived::bore2losAz
    if (present(losNadir)) losNadir = rvec2(iChan) !AntDerived::losNadir
    if (present(fscaleinc)) fscaleinc = rvec3(iChan)    
    if (present(aview_spec)) aview_spec = rvec4(iChan) !AntDerived::aview_spec
    if (present(ascan_spec)) ascan_spec = rvec5(iChan) !AntDerived::ascan_spec
    if (present(nedt_spec)) nedt_spec = rvec6(iChan)   !AntDerived::nedt_spec
    if (present(aview_nom)) aview_nom = rvec7(iChan) !AntDerived::aview_nom
    if (present(ascan_nom)) ascan_nom = rvec8(iChan) !AntDerived::ascan_nom

    deallocate(chanIDderived,rvec1,rvec2,rvec3,rvec4,rvec5,rvec6,rvec7,rvec8)

    close(iu)
    return

  END SUBROUTINE getSCFAntDerived

  !--------------------------------------------------------------------------

END MODULE SCFread

