MODULE CloudMitigateModule
  !
  ! Module for imposing limits over physical properties
  !
  ! Subroutines:
  !
  !     CloudClassifyInitial
  !
  ! Derived data types:
  !
  !     None
  !
  ! USE:
  !
  !     ControlStructure
  !     ImgObsStructure
  ! 
  ! yhe@aer.com, 03/16/2016
  !
  USE ControlStructure, Only: &
       GenControl_t

  USE ImgObsStructure, Only: &
       ImgOb_t

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: &
       CloudClassifyInitial

  integer, parameter :: mxLength = 256

CONTAINS

  subroutine CloudClassifyInitial(&
       genControl, &
       imgOb, &
       cloudClass &
       )
    TYPE(GenControl_t), intent(in) :: genControl
    TYPE(ImgOb_t), dimension(:), intent(in) :: imgOb
    integer, dimension(:), intent(out) :: cloudClass
    !Local
    character (len=mxLength) :: procName
    logical :: debug=.true.

    procName = '[CloudMitigateModule::CloudClassifyInitial]:'
    !if (debug) print *, trim(procName)//' starting ...'
    if (genControl%imgDataOn) cloudClass = imgOb(:)%cloudType !temporary
    !if (debug) print *, trim(procName)//' ending ...'
  end subroutine CloudClassifyInitial

END MODULE CloudMitigateModule
