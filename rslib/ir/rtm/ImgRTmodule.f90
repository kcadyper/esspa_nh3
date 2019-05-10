MODULE ImgRTmodule
  !
  ! Module for executing L1b image retrieval process
  !
  ! Subroutines:
  !
  !     ImgRTinit
  !     ImgRTexec
  !     ImgRTdestroy
  ! 
  ! Derived data type:
  !
  !     ImgRTcontrol_t
  !
  ! USE:
  !
  !     None
  !
  ! yhe@aer.com, 03/16/2016
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: &
       ImgRTinit, &
       ImgRTexec, &
       ImgRTdestroy, &
       ImgRTcontrol_t

  integer, parameter :: mxLength=256
  TYPE ImgRTcontrol_t
     !TBD
  END TYPE ImgRTcontrol_t

CONTAINS

  subroutine ImgRTinit(imgRTcontrol)
    TYPE(ImgRTcontrol_t), intent(in) :: imgRTcontrol
    !Local
    character (len=mxLength) :: procName

    procName = ' [ImgRTModule::ImgRTInit]: '
  end subroutine ImgRTinit

  subroutine ImgRTexec()
    !Local
    character (len=mxLength) :: procName

    procName = ' [ImgRTModule::ImgRTexec]: '
    print *, trim(procName)//'starting ...'
    print *, trim(procName)//'ending ...'

  end subroutine ImgRTexec

  subroutine ImgRTdestroy()
    !Local
    character (len=mxLength) :: procName

    procName = ' [ImgRTModule::ImgRTdestroy]: '
    print *, trim(procName)//'starting ...'
    print *, trim(procName)//'ending ...'

  end subroutine ImgRTdestroy
       
END MODULE ImgRTmodule
