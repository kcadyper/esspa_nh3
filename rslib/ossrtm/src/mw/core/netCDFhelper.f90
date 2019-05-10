!<f90File>**************************************************************
!
! CONTACT:
!
!   Atmospheric & Environmental Research, Inc
!   131 Hartwell Ave
!   Lexington ,MA 02421-3126 USA
!   Phone: 781.761.2288
!   E-mail: guymin@aer.com
!
! COPYRIGHT NOTICE:
!
!   Copyright AER, Inc 2001-2009, All Rights Reserved
!   See the file README-DATARIGHTS.txt included with this release
!   for additional details.
!
!*************************************************************</f90File>

module netCDFhelper
  implicit none
  include 'netcdf.inc'
  !--
  public :: openNcdfFile, closeNcdfFile, readNcdfAttr

  integer, parameter :: &
       MAX_NAME_LENGTH = 200, &
       MAXFILE = 50,&
       MAXVAR = 100,&
       MAX_ATTR_LENGTH = 2000

  !--
  logical, dimension(MAXFILE), private, save :: &
       ncid_index_active=.false.
  integer, dimension(MAXFILE), private, save :: ncid=0
  integer, dimension(MAXFILE), private, save :: unlimited_dim_id=0
  character (len=MAX_NAME_LENGTH), private, save :: &
       myfile = 'myncdf.nc',fStatus
  character (len=MAX_NAME_LENGTH), &
       dimension(MAXFILE), private, save :: myUnlimited_dim_name=' '
  integer, dimension(MAXVAR,MAXFILE), private, save :: varid=0
  integer, private :: ncStatus
  integer, private, save :: ncid_cnt_w=0, ncid_index=0,&
       unlimited_length
  logical, dimension(MAXVAR,MAXFILE), private, save :: &
       existed=.false.
  integer, dimension(MAXVAR,MAXFILE), private, save :: &
       record_no=1
  integer, private :: io, iopen
  !--
  private :: openFile, &
       readFloat1Attribute, readFloatAttribute,&
       readInt1Attribute, readIntAttribute,&
       readText1Attribute,&
       closeFile, getNewNcidIndex,&
       handle_err
  !--
  interface openNcdfFile
     module procedure openFile
  end interface
  !--
  interface readNcdfAttr
     module procedure readFloat1Attribute
     module procedure readFloatAttribute
     module procedure readInt1Attribute
     module procedure readIntAttribute
     module procedure readText1Attribute
  end interface
  !--
  interface closeNcdfFile
     module procedure closeFile
  end interface
  !--
contains
  subroutine openFile(fid, file, status, unlimited_dim_name, &
       unlimited_dim_length)
    character (len=*), intent (in), optional :: file, status
    integer, intent (inout) :: fid
    character (len = MAX_NAME_LENGTH), save :: localFile
    character (len=*), intent (in), optional :: &
         unlimited_dim_name
    integer, intent (inout), optional :: unlimited_dim_length
    integer :: dimid, ncidtmp

    if (present(file)) myfile = file
    if (present(status)) then
       fStatus = status
    else
       fStatus = 'old'
    endif

    if (fStatus == 'new') then
       if (.not. present(file))&
            print *,'creating NCDF File: ', trim(myFile)
       ncStatus = nf_create(myfile,nf_write,ncidtmp)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,trim(myfile)//' in openFile')
       if (ncStatus < 0) then
          return
       endif
       call getNewNcidIndex(ncid_index)
       ncid(ncid_index) = ncidtmp
       fid = ncid(ncid_index)
       if (present(unlimited_dim_name)) then
          myUnlimited_dim_name(ncid_index) = unlimited_dim_name
       else
          myUnlimited_dim_name(ncid_index) = 'unlimitedDim'
       endif
       ncStatus = nf_def_dim(fid, &
            trim(myUnlimited_dim_name(ncid_index)), &
            nf_unlimited, unlimited_dim_id(ncid_index))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,myUnlimited_dim_name(ncid_index))
       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,myfile)
    else
       if (.not. present(file))&
            print *,'opening NCDF File: ', trim(myFile)

       if (fStatus == 'append' .or. fStatus == 'replace') then
          ncStatus = nf_open(trim(myFile),nf_write,ncidtmp)
       else
          ncStatus = nf_open(myFile, nf_nowrite, ncidtmp)
       endif
       if (ncStatus .ne. nf_noerr) then
          call handle_err(ncStatus,myFile)
          call exit(1)
       endif
       if (ncStatus < 0) then
          return
       endif
       call getNewNcidIndex(ncid_index)
       ncid(ncid_index) = ncidtmp
       fid = ncid(ncid_index)
       ncStatus = nf_inq_unlimdim(fid, unlimited_dim_id(ncid_index))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,myUnlimited_dim_name(ncid_index))
       if (unlimited_dim_id(ncid_index) > 0) then
          ncStatus = nf_inq_dim(fid, unlimited_dim_id(ncid_index), &
               myUnlimited_dim_name(ncid_index), unlimited_length)
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,&
               myUnlimited_dim_name(ncid_index))
          if (present(unlimited_dim_length)) then
             unlimited_dim_length = unlimited_length
          endif
       endif
    endif

  end subroutine openFile
  !--
  subroutine readFloat1Attribute(fid, attr, attrName, silent, found)
    integer, intent (in) :: fid
    real*4, intent (out) :: attr
    character (len=*), intent (in) :: attrName
    logical, intent (in), optional :: silent
    logical, intent (out), optional :: found

    if (present(found)) then
       found = .true.
    endif

    ncStatus = nf_get_att_real(fid, nf_global, &
         trim(attrName), attr)
    if (ncStatus .ne. nf_noerr) then
       if (present(silent)) then
          if (silent) then
             attr = 0.0
             if (present(found)) then
                found = .false.
             endif
             return
          endif
       endif
       call handle_err(ncStatus,trim(attrName))
       attr = 0.0
    endif
    return
  end subroutine readFloat1Attribute
  !--
  subroutine readFloatAttribute(fid, attr, attrName, silent, found)
    integer, intent (in) :: fid
    real*4, dimension(:), intent (out) :: attr
    character (len=*), intent (in) :: attrName
    logical, intent (in), optional :: silent
    logical, intent (out), optional :: found

    if (present(found)) then
       found = .true.
    endif

    ncStatus = nf_get_att_real(fid, nf_global, &
         trim(attrName), attr)

    if (ncStatus .ne. nf_noerr) then
       if (present(silent)) then
          if (silent) then
             attr = 0.0
             if (present(found)) then
                found = .false.
             endif
             return
          endif
       endif
       call handle_err(ncStatus,trim(attrName))
       attr = 0.0
    endif
    return
  end subroutine readFloatAttribute
 !--
  subroutine readInt1Attribute(fid, attr, attrName, silent, found)
    integer, intent (in) :: fid
    integer, intent (out) :: attr
    character (len=*), intent (in) :: attrName
    logical, intent (in), optional :: silent
    logical, intent (out), optional :: found

    if (present(found)) then
       found = .true.
    endif

    ncStatus = nf_get_att_int(fid, nf_global, &
         trim(attrName), attr)

    if (ncStatus .ne. nf_noerr) then
       if (present(silent)) then
          if (silent) then
             attr = 0
             if (present(found)) then
                found = .false.
             endif
             return
          endif
       endif
       call handle_err(ncStatus,trim(attrName))
       attr = 0
    endif
    return
  end subroutine readInt1Attribute
  !--
  subroutine readIntAttribute(fid, attr, attrName, silent, found)
    integer, intent (in) :: fid
    integer, dimension(:), intent (out) :: attr
    character (len=*), intent (in) :: attrName
    logical, intent (in), optional :: silent
    logical, intent (out), optional :: found

    if (present(found)) then
       found = .true.
    endif

    ncStatus = nf_get_att_int(fid, nf_global, &
         trim(attrName), attr)

    if (ncStatus .ne. nf_noerr) then
       if (present(silent)) then
          if (silent) then
             attr = 0
             if (present(found)) then
                found = .false.
             endif
             return
          endif
       endif
       call handle_err(ncStatus,trim(attrName))
       attr = 0
    endif
    return
  end subroutine readIntAttribute
  !--
  subroutine readText1Attribute(fid, attr, attrName, silent, found)
    integer, intent (in) :: fid
    character (len=*), intent (out) :: attr
    character (len=*), intent (in) :: attrName
    logical, intent (in), optional :: silent
    logical, intent (out), optional :: found
    integer :: attrLen, copyLen
    character (len=MAX_ATTR_LENGTH) :: attrLocal

    if (present(found)) then
       found = .true.
    endif

    ncStatus = nf_inq_attlen(fid, nf_global, &
         trim(attrName), attrLen)

    if (ncStatus .ne. nf_noerr) then
       if (present(silent)) then
          if (silent) then
             attr = '\0'
             if (present(found)) then
                found = .false.
             endif
             return
          endif
       endif
       call handle_err(ncStatus,trim(attrName))
       attr = '\0'
       return
    endif

    if (attrLen > MAX_ATTR_LENGTH) then
       print *,'err[netCDFhelper::readText1Attribute]: ', &
          'Attribute length exceeds max;'
       print *,'attrName,attrLen,MAX_ATTR_LENGTH: ', &
          trim(attrName),attrLen,MAX_ATTR_LENGTH
       call exit(1)
    endif

    ncStatus = nf_get_att_text(fid, nf_global, &
         trim(attrName), attrLocal)

    if (ncStatus .ne. nf_noerr) then
       call handle_err(ncStatus,trim(attrName))
       attr = '\0'
    endif

    attr=REPEAT(' ',LEN(attr))   ! Pad with blanks
    copyLen=min(attrLen,len(attr))
    attr(1:copyLen)=attrLocal(1:copyLen)

    return
  end subroutine readText1Attribute
  !--
  subroutine closeFile(fid, message)
    character (len=*), intent (in), optional :: message
    integer, intent (in) :: fid

    if (present(message)) then
       print *, '------------------------'
       print *, 'closing unit: ', fid
       print *, message(1:len_trim(message))
       print *, '------------------------'
    endif
    ncStatus = nf_close(fid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,'closing')

    do io = 1, ncid_cnt_w
       if (ncid(io) == fid) then
          iopen = io
          exit
       endif
    enddo

    !--   Clear the index slot so it can be reused
    ncid_index_active(iopen) = .false.
    ncid(iopen) = 0
    unlimited_dim_id(iopen) = 0
    myUnlimited_dim_name(iopen)=' '
    varid(:,iopen) = 0
    existed(:,iopen) = .false.
    record_no(:,iopen) = 1
  end subroutine closeFile
  !--
  subroutine handle_err(ncStatus,msg)
    integer :: ncStatus
    character (len=*), intent (in), optional :: msg
    character (len=256) :: lmsg

    if (present(msg)) then
       if (len(msg) > 256) print*,'message is truncated'
       lmsg=trim(msg)
    else
       lmsg=':::'
    endif

    SELECT CASE (ncStatus)
    CASE (2)
       print *, 'err[netCDFhelper::handle_err]: ',&
            ' ('//trim(msg)//') '//trim(nf_strerror(ncStatus))
       call exit(1)
    CASE DEFAULT
       print *, 'warning[netCDFhelper::handle_err]: ',&
            ' ('//trim(msg)//') ',trim(nf_strerror(ncStatus)),ncStatus
    END SELECT
  end subroutine handle_err
  !--
  subroutine getNewNcidIndex(ncid_index)

    implicit none

    integer, intent (out) :: ncid_index

    !--   Local variables
    integer :: io

    do io = 1, ncid_cnt_w
       if (.not. ncid_index_active(io)) then
          ncid_index = io
          ncid_index_active(io) = .true.
          return
       endif
    enddo

    ncid_cnt_w = ncid_cnt_w + 1
    if (ncid_cnt_w .gt. MAXFILE) then
       print *,'err[netCDFhelper::getNewNcidIndex]: ',&
            ' ncid_cnt_w overran MAXFILE'
       call exit(1)
    end if
    ncid_index = ncid_cnt_w
    ncid_index_active(ncid_index) = .true.
    return
  end subroutine getNewNcidIndex

end module netCDFhelper


