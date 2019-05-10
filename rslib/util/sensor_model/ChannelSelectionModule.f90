MODULE ChannelSelectionModule
  !
  ! Module for channel selection
  !
  ! Subroutines:
  !
  !      channelSelectionInit
  !      channelSetUnion
  !      channelSelectionDestroy
  !
  ! Derived data types:
  !
  !      ChannelSelections_t
  !      ChannelSet_t
  !
  ! USE
  !
  !      ControlStructure
  !      MWobsStructure
  !      IRobsStructure
  !
  ! yhe@aer.com, 03/16/2016
  !
  USE ControlStructure, Only: GenControl_t
  USE MWobsStructure, Only: MWdefin_t
  USE IRobsStructure, Only: IRdefin_t

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: &
       channelSelectionInit, &
       channelSetUnion, &
       channelSelectionDestroy, &
       ChannelSelections_t, &
       ChannelSet_t

  !-------------------------------------------------------------
  TYPE ChannelSet_t
     integer :: chanSetInd   ! channel subset index (use to set up oss ir module)
     logical :: defined = .False.
     integer :: nMW
     integer :: nIR
     integer, dimension(:), allocatable :: MW  !size=nMW
     integer, dimension(:), allocatable :: IR  !size=nIR
  END TYPE ChannelSet_t

  TYPE ChannelSelections_t
     TYPE(ChannelSet_t) :: MW
     TYPE(ChannelSet_t) :: CloudClearing
     TYPE(ChannelSet_t) :: PrimaryClear
     TYPE(ChannelSet_t) :: PrimaryCloudy
     TYPE(ChannelSet_t) :: Secondary
  END TYPE ChannelSelections_t

  logical, parameter          :: debug=.false.

CONTAINS

  subroutine channelSelectionInit(genControl,MWdefin,IRdefin,channelSelections)
    TYPE(GenControl_t), intent(in) :: genControl
    TYPE(MWdefin_t), intent(in) :: MWdefin
    TYPE(IRdefin_t), intent(in) :: IRdefin
    TYPE(ChannelSelections_t), intent(out) :: channelSelections
    ! local
    character (len=*), parameter :: procName = &
                           ' [ChannelSelectionModule::channelSelectionInit]: '

    if (debug) print *, procName, 'starting ...'

    ! currently file must be in netCDF fromat
    call channelSelectionNetCDFFileRead(genControl, &
                          MWdefin%nChan, IRdefin%nChan, &
                          channelSelections)

!  The adequacy of the channelSelection file is not checked here, it is deferred
!  to the appropriate place

    if (debug) call channelSelectionPrint(channelSelections)

    if (debug) print *, procName, 'ending ...'

  end subroutine channelSelectionInit

  function channelSetUnion(genControl,channelSelections)
! Compose a channel set that is the union of all the channel sets
    TYPE(GenControl_t),        intent(in) :: genControl
    TYPE(ChannelSelections_t), intent(in) :: channelSelections
    TYPE(ChannelSet_t) :: channelSetUnion
    integer :: nMW
    integer :: nIR
    integer :: nPrimMW
    integer :: nPrimIR
    character (len=*), parameter :: procName = &
                     ' [ChannelSelectionModule::channelSetUnion]: '

   ! Define vector sizes
    nMW = 0
    nIR = 0
    if (genControl%MWretrOn) nMW = channelSelections%MW%nMW
    if (genControl%primaryRetrOn) then
       if (genControl%MWdataOn) then
          nMW = max(nMW, channelSelections%PrimaryClear%nMW)
          nMW=max(nMW, channelSelections%PrimaryCloudy%nMW)
          nMW=max(nMW, channelSelections%CloudClearing%nMW)
          nMW=max(nMW,nPrimMW)
       end if

       nIR = max(nIR,channelSelections%PrimaryClear%nIR)
       nIR=max(nIR, channelSelections%PrimaryCloudy%nIR)
       nIR=max(nIR, channelSelections%CloudClearing%nIR)
     endif

     channelSetUnion%nMW=nMW
     channelSetUnion%nIR=nIR
     if (nMW > 0) then
        if (.not. allocated(channelSetUnion%MW)) &
            allocate(channelSetUnion%MW(nMW))
        channelSetUnion%MW=0
        if (genControl%MWretrOn) then
           where (channelSelections%MW%MW > 0)
              channelSetUnion%MW = 1
           end where
        endif
        if (genControl%primaryRetrOn .and. genControl%MWdataOn) then
           if (channelSelections%PrimaryClear%nMW == nMW) then
               where (channelSelections%PrimaryClear%MW > 0)
                  channelSetUnion%MW = 1
               end where
           endif
           if (channelSelections%CloudClearing%nMW == nMW) then
               where (channelSelections%CloudClearing%MW > 0)
                  channelSetUnion%MW = 1
               end where
           endif
           if (channelSelections%PrimaryCloudy%nMW == nMW) then
               where (channelSelections%PrimaryCloudy%MW > 0)
                  channelSetUnion%MW = 1
               end where
           endif
        endif
     endif
     if (nIR > 0) then
        if (.not. allocated(channelSetUnion%IR)) &
            allocate(channelSetUnion%IR(nIR))
        channelSetUnion%IR=0
        if (genControl%primaryRetrOn) then
           if (channelSelections%PrimaryClear%nIR == nIR) then
               where (channelSelections%PrimaryClear%IR > 0)
                  channelSetUnion%IR = 1
               end where
           end if
           if (channelSelections%CloudClearing%nIR == nIR) then
               where (channelSelections%CloudClearing%IR > 0)
                  channelSetUnion%IR = 1
               end where
           end if
           if (channelSelections%PrimaryCloudy%nIR == nIR) then
               where (channelSelections%PrimaryCloudy%IR > 0)
                  channelSetUnion%IR = 1
               end where
           end if
        end if
     end if
  end function channelSetUnion

  subroutine channelSelectionDestroy(channelSelections)
    TYPE(ChannelSelections_t), intent(inout) :: channelSelections
    character (len=*), parameter :: procName = &
                 ' [ChannelSelectionModule::channelSelectionDestroy]: '

    if (debug) print *, procName, 'starting ...'
    if (allocated(channelSelections%MW%MW)) &
         deallocate(channelSelections%MW%MW)
    if (allocated(channelSelections%PrimaryClear%MW)) &
         deallocate(channelSelections%PrimaryClear%MW)
    if (allocated(channelSelections%PrimaryClear%IR)) &
         deallocate(channelSelections%PrimaryClear%IR)
    if (allocated(channelSelections%PrimaryCloudy%MW)) &
         deallocate(channelSelections%PrimaryCloudy%MW)
    if (allocated(channelSelections%PrimaryCloudy%IR)) &
         deallocate(channelSelections%PrimaryCloudy%IR)
    if (allocated(channelSelections%CloudClearing%MW)) &
         deallocate(channelSelections%CloudClearing%MW)
    if (allocated(channelSelections%CloudClearing%IR)) &
         deallocate(channelSelections%CloudClearing%IR)
    if (allocated(channelSelections%Secondary%MW)) &
         deallocate(channelSelections%Secondary%MW)
    if (allocated(channelSelections%Secondary%IR)) &
         deallocate(channelSelections%Secondary%IR)
    if (debug) print *, procName, 'ending ...'
  end subroutine channelSelectionDestroy

  subroutine channelSelectionNetCDFFileRead(genControl, numChanMW, &
                                            numChanIR, channelSelections)
    use netcdf
    TYPE(GenControl_t),                 intent(in) :: genControl
    integer,                            intent(in)  :: numChanMW
    integer,                            intent(in)  :: numChanIR

    TYPE(ChannelSelections_t),          intent(out) :: channelSelections

    !local
    character (len=*), parameter :: procName = &
                 ' [ChannelSelectionModule::channelSelectionNetCDFFileRead]: '
    integer                               :: ncid
    integer                               :: listSize
    integer, dimension(nf90_max_var_dims) :: dimIDs

    integer, dimension(:), allocatable    :: list

    print *, procName, 'openning ...', &
            trim(genControl%channelSelectFiles(1))
    listSize=max(numChanMW,numChanIR)
    allocate(list(listSize))

    call handle_err(nf90_open(genControl%channelSelectFiles(1), NF90_NOWRITE, ncid))
    channelSelections%MW%defined=.false.
    channelSelections%MW%nIR = -1
    channelSelections%MW%chanSetInd = -1
    if (genControl%MWdataOn) call readFieldMW('MW',numChanMW,channelSelections%MW)

    channelSelections%CloudClearing%defined=.false.
    channelSelections%CloudClearing%chanSetInd = -1
    if (genControl%IRdataOn) &
         call readFieldIR('CloudClearingIR',numChanIR,channelSelections%CloudClearing)

    if (genControl%MWdataOn) &
         call readFieldMW('CloudClearingMW',numChanMW,channelSelections%CloudClearing)

    channelSelections%PrimaryClear%defined=.false.
    channelSelections%PrimaryClear%chanSetInd = -1
    if (genControl%IRdataOn)  &
      call readFieldIR('PrimaryClearIR',numChanIR,channelSelections%PrimaryClear)
    if (genControl%MWdataOn) &
      call readFieldMW('PrimaryClearMW',numChanMW,channelSelections%PrimaryClear)

    channelSelections%PrimaryCloudy%defined=.false.
    channelSelections%PrimaryCloudy%chanSetInd = -1
    if (genControl%IRdataOn) &
      call readFieldIR('PrimaryCloudyIR',numChanIR,channelSelections%PrimaryCloudy)
    if (genControl%MWdataOn) &
      call readFieldMW('PrimaryCloudyMW',numChanMW,channelSelections%PrimaryCloudy)

    channelSelections%Secondary%defined=.false.
    channelSelections%Secondary%chanSetInd = -1
    channelSelections%Secondary%nMW = -1
    if (genControl%IRdataOn) &
      call readFieldIR('SecondaryIR',numChanIR,channelSelections%Secondary)

    call handle_err(nf90_close(ncid))
    if (debug) print *, procName, 'ending ...'

  CONTAINS
    subroutine readFieldMW(baseName, numEl, chanSet)
      character(len=*),   intent(in)    :: baseName
      integer,            intent(in)    :: numEl
      TYPE(ChannelSet_t), intent(inout) :: chanSet
      integer                           :: varId, ncStatus, varSize, jj, allChan
      ncStatus = nf90_inq_varid(ncid, baseName, varid=varId)
      if (ncStatus==nf90_NoErr) then
        chanSet%defined=.true.
        call handle_err(nf90_get_att(ncid, varid, 'all_channels', allChan))
        chanSet%nMW=numEl
        allocate(chanSet%MW(numEl))
        if (allChan>0) then
           chanSet%MW = 1
           return
        else
           chanSet%MW = 0
        end if
        call handle_err(nf90_inquire_variable(ncid, varId, dimids = dimIDs))
        call handle_err(nf90_inquire_dimension(ncid, dimIDs(1), len = varSize))
        if (numEl<varSize) then
          print *, procName, 'MW grid size (', numEl,&
               ') less then channel selection list', varSize, ' for ', baseName
          call exit(1)
        end if
        call handle_err(nf90_get_var(ncid, varId, list(1:varSize)))
        do jj=1,varSize
          if (list(jj)>numEl) then
             print *, procName, 'Channel index (', &
                   list(jj),') exceeds the MW grid size (', numEl,&
                  ') for channel list', baseName
             call exit(1)
          end if
          chanSet%MW(list(jj)) = 1
        end do

      else
        chanSet%nMW=-1
      end if
    end subroutine readFieldMW

    subroutine readFieldIR(baseName, numEl, chanSet)
      character(len=*),   intent(in)    :: baseName
      integer,            intent(in)    :: numEl
      TYPE(ChannelSet_t), intent(inout) :: chanSet

      integer                           :: varId, ncStatus, varSize, jj, allChan
      ncStatus = nf90_inq_varid(ncid, baseName, varid=varId)
      if (ncStatus==nf90_NoErr) then
        chanSet%defined=.true.
        ! shortcut
        ! check if the attribute 'all_channel' set to 1
        call handle_err(nf90_get_att(ncid, varid, 'all_channels', allChan))
        chanSet%nIR=numEl
        allocate(chanSet%IR(numEl))
        if (allChan>0) then
           chanSet%IR = 1
           return
        else
           chanSet%IR = 0
        end if
        call handle_err(nf90_inquire_variable(ncid, varId, dimids = dimIDs))
        call handle_err(nf90_inquire_dimension(ncid, dimIDs(1), len = varSize))
        print *,'basename, varsize, numel ',basename,"  ",varsize,numel
        if (numEl<varSize) then
          print *, procName, 'IR grid size (', numEl,&
               ') less then channel selection list', varSize, ' for ', baseName
          call exit(1)
        end if
        call handle_err(nf90_get_var(ncid, varId, list(1:varSize)))
!!!!! Another Big fudge here KCP !!!!
        !do jj=1,varSize
        do jj=1,varSize-1
          if (list(jj)>numEl) then
             print *, procName, 'Channel index (', &
                   list(jj),') exceeds the IR grid size (', numEl,&
                  ') for channel list', baseName
             call exit(1)
          end if
!!!!! Big fudge here KCP !!!!
          !chanSet%IR(list(jj)) = 1
          chanSet%IR(list(jj)+2) = 1
        end do
      else
        chanSet%nIR=-1
      end if
    end subroutine readFieldIR

    !handle netcdf errors
    SUBROUTINE handle_err(status)
      integer, intent (in) :: status
      !local
      character(len=80) :: msg
      if (status==nf90_NoErr) return

      msg = nf90_strerror(status)
      print *, procName, trim(msg)
      if (index(msg,'warning') ==0) call exit(1)
    END SUBROUTINE handle_err

  end subroutine channelSelectionNetCDFFileRead

  subroutine channelSelectionPrint(chanSelect)
    TYPE(ChannelSelections_t),          intent(in) :: chanSelect
    print *,'ChannelSelections_t'

    if (chanSelect%MW%defined) then
      print *,'MW defined', chanSelect%MW%nMW
    else
      print *,'MW NOT defined'
    end if
    if (chanSelect%CloudClearing%defined) then
      print *,'CloudClearing defined', chanSelect%CloudClearing%nMW, &
              chanSelect%CloudClearing%nIR
    else
      print *,'CloudClearing NOT defined'
    end if
    if (chanSelect%PrimaryClear%defined) then
      print *,'PrimaryClear defined', chanSelect%PrimaryClear%nMW, &
              chanSelect%PrimaryClear%nIR
    else
      print *,'PrimaryClear NOT defined'
    end if
    if (chanSelect%PrimaryCloudy%defined) then
      print *,'PrimaryCloudy defined', chanSelect%PrimaryCloudy%nMW, &
              chanSelect%PrimaryCloudy%nIR
    else
      print *,'PrimaryCloudy NOT defined'
    end if
    if (chanSelect%Secondary%defined) then
      print *,'Secondary defined', chanSelect%Secondary%nMW, &
              chanSelect%Secondary%nIR
    else
      print *,'Secondary NOT defined'
    end if
  end subroutine channelSelectionPrint


END MODULE ChannelSelectionModule
