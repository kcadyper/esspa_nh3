!<f90File>**************************************************************
!
! CONTACT:
!
!   Atmospheric & Environmental Research, Inc
!   131 Hartwell Ave
!   Lexington ,MA 02421-3126 USA
!   Phone: 781.761.2288
!   E-mail: yhe@aer.com
!
! COPYRIGHT NOTICE:
!
!   Copyright AER, Inc 2016-, All Rights Reserved
!   See the file README-DATARIGHTS.txt included with this release
!   for additional details.
!
!*************************************************************</f90File>

MODULE CloudNameDictionary
! <f90Module>***********************************************************
!
! NAME:
!
!   CloudNameDictionary
!
! PURPOSE:
!
!   This is a utility module that translates or maps the optical field
!   names in a cloud optical table file for the current implementation.
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>
  USE CloudParameters
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: matchName

  character (len=mxMsgLength) :: msg

CONTAINS

  LOGICAL FUNCTION matchName(srcName,key,prefix)
    character (len=*), intent(in) :: srcName, key
    character (len=*), intent(in), optional :: prefix
    logical :: found
    integer :: ix
    
    found = .false.
    if (present(prefix)) then
       ix = index(trim(srcName),trim(prefix)//trim(lookup(trim(key))))
    else
       ix = index(trim(srcName),trim(lookup(trim(key))))
    end if
    if (ix /= 0) then 
       found = .true.
    end if
    matchName = found

  END FUNCTION matchName

  character (len=mxNameLength) function lookup(key)
    character (len=*), intent(in) :: key
    
    select case (to_upper(trim(key)))
      case('SIZP')
         lookup = 'SizP'
      case('SZP')
         lookup = 'SizP'
      case('SIZS')
         lookup = 'SizS'
      case('SZS')
         lookup = 'SizS'
      case('T')
         lookup = 'Temp'
      case('TEMP')
         lookup = trim(key)
      case('DENS')
         lookup = trim(key)
      case('MELTF')
         lookup = trim(key)
      case default
         print *, '[error in CloudNameDictionary::lookup] ... '//trim(key)//' not recognized'
         call exit(101)
    end select

  end function lookup

  !===============================
  !Converts a string to upper case
  !===============================
  function to_upper (strIn) result (strOut)

    character(*), intent(in) :: strIn
    character(LEN(strIn))    :: strOut

    integer :: ic, i

    character(26), parameter :: CAP = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

    strOut = strIn
    do i = 1, len_trim(strIn)
        ic = index(low, strIn(i:i))
        if (ic > 0) strOut(i:i) = CAP(ic:ic)
    end do

  end function to_upper
  
  !===============================
  !Converts a string to lower case
  !===============================
  function to_lower (strIn) result (strOut)

    character(*), intent(in) :: strIn
    character(len(strIn))    :: strOut

    Integer :: ic, i

    character(26), parameter :: CAP = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

    strOut = strIn
    do i = 1, len_trim(strIn)
        ic = index(CAP, strIn(i:i))
        if (ic > 0) strOut(i:i) = low(ic:ic)
    end do

  end function to_lower

END MODULE CloudNameDictionary
