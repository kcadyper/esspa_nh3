MODULE ImgObsStructure
  !
  ! Module that defines all the imagery data structures
  !
  ! yhe@aer.com, 06/14/2016
  !
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: &
       ImgDef_t, &
       ImgOb_t, &
       ImgGran_t

  TYPE ImgDef_t
     integer :: nFOV
     integer, dimension(6) :: startTime ! year,month,day,hour,minute,second
     integer, dimension(6) :: endTime
  END TYPE ImgDef_t

  TYPE ImgOb_t
     integer :: cloudType
  END TYPE ImgOb_t

  TYPE ImgGran_t
     integer :: nFOR
     TYPE(ImgDef_t) :: def
     TYPE(ImgOb_t), dimension(:,:), allocatable :: ob  !size=(nFOV,nFOR)
  END TYPE ImgGran_t

END MODULE ImgObsStructure
