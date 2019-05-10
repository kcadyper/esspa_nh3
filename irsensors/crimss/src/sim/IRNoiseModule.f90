!-----------------------------------------------------
!
!  MODULE IRNoise: contains subroutines needed for the
!                IR noise generation. AER Inc. 2004.
!-----------------------------------------------------
MODULE IRNoiseModule
  IMPLICIT NONE
  PRIVATE
  !------------------------------------------------------------------------
  !	Public items made available to calling program(s)
  !------------------------------------------------------------------------
  PUBLIC :: initIRdeviceNoise,IRaddDeviceNoise,IRdeviceNoise,loadIRdeviceNoise, &
       destroyIRdeviceNoise

  !------------------------------------------------------------------------
  !     Declarations of private data (private to the module)
  !------------------------------------------------------------------------
  integer,            parameter       :: mxLength=256
  INTEGER,            parameter       :: mxchan=90000
  INTEGER                             :: nfrq,ntemp,myNch
  integer, dimension(:), allocatable  :: chanNumLocal
  REAL, DIMENSION(:), allocatable     :: bnfreq,AtmNoiseBT
  REAL, DIMENSION(:,:), allocatable   :: SceneTemp,IrDevNoise

  ! loadDone implementation provides the option to call loadIRdeviceNoise() in advance
  logical                             :: loadDone=.false.
  INTEGER                             :: initNoise=1
  logical                             :: dbg = .false.

CONTAINS

  !called by IRmeasErrModule::IRNoiseInit per ESDR module design
  SUBROUTINE initIRdeviceNoise(nChan,nednFile)
    integer, intent(in) :: nChan
    character (len=*), intent(in) :: nednFile
    !Local
    character (len=mxLength) :: procName

    procName = ' [IRNoiseModule::initIRdeviceNoise]: '
    if (dbg) print *, trim(procName)//'starting ...'

    call loadIRdeviceNoise(trim(nednFile))

    if(myNch /= nchan) then
       print*,procName//'err - Inconsistent number of channels: ', myNch,nChan
       call errorHalt(1)
    endif
    allocate(AtmNoiseBT(nChan))
    AtmNoiseBT=0.
    initNoise=0
    if (dbg) print *, trim(procName)//'ending ...'
  END SUBROUTINE initIRdeviceNoise

  SUBROUTINE IRaddDeviceNoise(frq,Y,DevNoise,sumIRDevNoise,kchan,nchan,chanIDreq)
    !---In/Out variables
    INTEGER               :: nchan
    LOGICAL, DIMENSION(:) :: kchan
    REAL,    DIMENSION(:) :: frq,Y,DevNoise
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN)    :: chanIDreq
    !---Local variables
    INTEGER, SAVE         :: iseed=1287365
    INTEGER               :: i
    REAL                  :: stdt,gauss
    REAL                  :: sumIRDevNoise

    if (.not. allocated(chanNumLocal)) allocate(chanNumLocal(nchan))
    call IRdeviceNoise(frq,Y,DevNoise,sumIRDevNoise=sumIRDevNoise, &
         kchan=kchan,nchan=nchan,chanIDreq=chanIDreq,chanNum=chanNumLocal)
    DO i=1,nchan
       stdt=devNoise(i)
       Y(i)=Y(i)+gauss(iseed,stdt,0.0)
    END DO
    RETURN
  END SUBROUTINE IRaddDeviceNoise

  SUBROUTINE IRdeviceNoise(frq,Y,DevNoise,AtmNoiseR,sumIRDevNoise,kchan,nchan, &
       chanIDreq,chanNum)
    use IRReadStdInputs, only: F_nedn
    !---In/Out variables
    INTEGER, intent(in), optional   :: nchan
    LOGICAL, DIMENSION(:), intent(in) :: kchan
    integer, dimension(:), intent(in) :: chanNum
    REAL,    DIMENSION(:), intent(in) :: frq,Y
    REAL,    DIMENSION(:), intent(out) :: DevNoise
    REAL,    DIMENSION(:), intent(inout), optional :: AtmNoiseR
    REAL,                  intent(out) :: sumIRDevNoise
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN), optional :: chanIDreq
    !---Local variables
    REAL                            :: v,Ytmp,tb,bt,draddt,stdi,stdt
    INTEGER                         :: i,j,nch,ib,iavg
    character (len=mxLength) :: procName

    procName = '[IRNoiseModule::IRdeviceNoise]:'
    if (dbg) print *, trim(procName)//' starting ...'
    !------------------------------------------------------------------------
    ! Set global variables from file variables
    !------------------------------------------------------------------------

    IF(initNoise == 1)THEN
       call initIRdeviceNoise(nChan,F_nedn)
    END IF

    !Comment out due to the introduction of initIRdeviceNoise per ESDR
    !modularization design
!!$    call loadIRdeviceNoise(F_nedn)
!!$
!!$    if(myNch /= nchan) then
!!$       print*,'err[IRNoiseModule::IRdeviceNoise]: '
!!$       print*,'Inconsistent number of channels: '
!!$       print*,myNch,nchan
!!$       call errorHalt(1)
!!$    endif
!!$
!!$    IF(init == 1)THEN
!!$       allocate(AtmNoiseBT(nchan))
!!$       AtmNoiseBT=0.
!!$
!!$       INIT=0
!!$    endif

    ! Channel IDs are currently unused for this sensor

    ! -----------------------------------------------------
    !  actual IR noise amplitude based on scene temperature
    !------------------------------------------------------

    sumIRDevNoise=0.
    nch=0
    YTmp=50.0
    iavg=0

!!$    do i=1,nchan
    do i=1,myNch
       v=frq(i)
       if(kchan(i))then
          YTmp=max(Y(i),1.e-3)
       else
          YTmp=YTmp
       end if
       tb=bt(v,YTmp)
       if(iavg == 0) then
          ! determine which instrument noise value to use
          ! based on scene temperature tb
          if(v < bnfreq(2))then
             ib=1
          else if(v >= bnfreq(nfrq)) then
             ib=nfrq
          else
             do j=2,nfrq-1
                if(v >= bnfreq(j) .and. v < bnfreq(j+1))then
                   ib=j
                endif
             enddo
          end if
          if(tb <= SceneTemp(ib,1)) then
             stdi=xlint(IrDevNoise(i,1),IrDevNoise(i,2), &
                  SceneTemp(ib,1),SceneTemp(ib,2),tb)
          else if(tb <= SceneTemp(ib,2)) then
             stdi=xlint(IrDevNoise(i,1),IrDevNoise(i,2), &
                  SceneTemp(ib,1),SceneTemp(ib,2),tb)
          else if(tb <= SceneTemp(ib,3)) then
             stdi=xlint(IrDevNoise(i,2),IrDevNoise(i,3), &
                  SceneTemp(ib,2),SceneTemp(ib,3),tb)
          else
             stdi=xlint(IrDevNoise(i,2),IrDevNoise(i,3), &
                  SceneTemp(ib,2),SceneTemp(ib,3),tb)
          end if
          DevNoise(i)=stdi
          if(present(AtmNoiseR))AtmNoiseR(i)=draddt(v,tb)*AtmNoiseBT(i)
       else                   ! use average value of all scene temperatures
          stdi=draddt(v,tb)*AtmNoiseBT(i)
       end if
       stdt=stdi
       if(kchan(i)) then
          nch=nch+1
          sumIRDevNoise=sumIRDevNoise+stdt**2
       end if
    end do
    sumIRDevNoise=sqrt(sumIRDevNoise/nch)

    if (dbg) print *, trim(procName)//' ending ...'
  END SUBROUTINE IRdeviceNoise

  subroutine loadIRdeviceNoise(file_nedn)
!!$    use IRReadStdInputs, only: U_nedn
    Use ToolboxModule, Only: getUnit

    !--I/O variables
    character(len=*),              intent(in)  :: file_nedn
    !--Local variables
    integer                                    :: myNfrq
    integer                                    :: i,j,ict
    integer                                    :: U_nedn
    real                                       :: dum
    character (len=mxLength) :: procName

    procName = '[IRNoiseModule::loadIRdeviceNoise]:'
    if (dbg) print *, trim(procName)//' starting ...'

    if(loadDone) return

    U_nedn = getUnit()
    !------------------------------------------------------------------------
    ! get number of channels in file (reading twice)
    !------------------------------------------------------------------------
    open(U_nedn,file=file_nedn,status='old')

    !count total and subtract out non-channel-related lines
    read(U_nedn,*)myNfrq,dum
    ict=1
    do i=1,mxchan
       read(U_nedn,*,end=100) dum
       ict=ict+1
    enddo

    close(U_nedn)
    print*,'err[IRNoiseModule::loadIRdeviceNoise]: '
    print*,'mxchan limit reached: '
    print*,ict,mxchan
    call errorHalt(1)

100 continue
    myNch=ict-(1+2*myNfrq)

    !------------------------------------------------------------------------
    ! Set global variables from file variables
    !------------------------------------------------------------------------
    rewind U_nedn

    !-- Read in end points of noise band calculations
    read(U_nedn,*)nfrq,ntemp

    allocate(bnfreq(nfrq))
    do i=1,nfrq
       read(U_nedn,*) bnfreq(i)
    end do

    !-- Read in scene temperatures for which NedNs were calculated
    allocate(SceneTemp(nfrq,ntemp))
    do i=1,nfrq
       read(U_nedn,*) (SceneTemp(i,j),j=1,ntemp)
    end do

    !-- Read in radiometric sensor noise values
    allocate(IrDevNoise(myNch,nfrq))
    do i=1,myNch
       read(U_nedn,*) (IrDevNoise(i,j),j=1,nfrq)
    end do
    close(U_nedn)

    loadDone=.true.
    if (dbg) print *, trim(procName)//' ending ...'
    return
  end subroutine loadIRdeviceNoise

  REAL function xlint(d1,d2,x1,x2,x)
    IMPLICIT NONE
    REAL :: d1,d2,x1,x2,x
    REAL :: w1,w2,xval
    if(x2.ne.x1) then
       w1=(x2-x)/(x2-x1)
       w2=1-w1
       xval=w1*d1+w2*d2
    else
       xval=d1
    endif
    xlint=xval

    return
  end function xlint

  !clean up and reset initial values
  subroutine destroyIRdeviceNoise()
    loadDone = .false.
    initNoise = 1
    if (allocated(AtmNoiseBT)) deallocate(AtmNoiseBT)
    if (allocated(chanNumLocal)) deallocate(chanNumLocal)
    if (allocated(SceneTemp)) deallocate(SceneTemp)
    if (allocated(IrDevNoise)) deallocate(IrDevNoise)
  end subroutine destroyIRdeviceNoise

END MODULE IRNoiseModule
