MODULE IRChannels

  implicit none
  private
  public :: initIRChannels, loadIRChannels

  !------------------------------------------------------------------------
  !     Declarations of private data (private to the module)
  !------------------------------------------------------------------------
  integer, parameter                 :: mxchan=90000
  integer, dimension(:), allocatable :: myKchan
  integer                            :: myNch
  ! loadDone implementation provides the option to call loadIRChannels() in advance
  logical                            :: loadDone=.false.

CONTAINS

  subroutine initIRChannels(nchanmw,nchan,frq,kchan)
    use IRReadStdInputs, only: F_chan

    !--I/O variables
    integer :: nchanmw,nchan
    real, dimension(:)    :: frq
    logical, dimension(:) :: kchan
    !--Local variables
    integer :: i

    !-----------------------------------
    !-- IR channel selection flags
    !-- They are read from a file
    !-- and can be modifed here as needed
    !-----------------------------------
    !------------------------------------------------------------------------
    ! Set global variables from file variables
    !------------------------------------------------------------------------
       call loadIRChannels(F_chan)

       if(myNch /= (nchan-nchanmw)) then
          print*,'err[IRChannels::initIRChannels]: '
          print*,'Inconsistent number of IRChannels: '
          print*,myNch,(nchan-nchanmw)
          call errorHalt(1)
       endif

    !------------------------------------------------------------------------
    ! Copy to output arguments
    !------------------------------------------------------------------------
    kchan(nchanmw+1:nchan)=myKchan(:) /= 0

  end subroutine initIRChannels

  subroutine loadIRChannels(file_chan)
    use IRReadStdInputs, only: U_chan

    !--I/O variables
    character(len=*),              intent(in)  :: file_chan
    !--Local variables
    real                                       :: flocal,dum
    integer                                    :: i,ict
    integer,                       parameter   :: mxchan=90000

    if(loadDone) return

    ! get number of channels in file

    open(U_chan, file=file_chan, status='old')
    ict=0
    do i=1,mxchan
       read(U_chan,*,end=100) flocal,dum,dum,dum,dum,dum,dum,dum
       ict=ict+1
    enddo

    close(U_chan)
    print*,'err[IRChannels::loadIRChannels]: '
    print*,'mxchan limit reached: '
    print*,ict,mxchan
    call errorHalt(1)

100 continue
    myNch=ict

    allocate(myKchan(myNch))

    rewind U_chan
    do i=1,myNch
       read(U_chan,*) flocal,dum,dum,dum,dum,dum,dum,myKchan(i)
    end do
    close(U_chan)

    loadDone=.true.

  end subroutine loadIRChannels

END MODULE IRChannels
