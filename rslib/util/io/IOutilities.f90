!----------------------------------------------------------------------------
!
! MODULE IOutilities
!
! This module contains utilities to support other modules that do IO  
!
!----------------------------------------------------------------------------
MODULE IOutilities

  USE ncdf_module, ONLY: openNcdfFile,closeNcdfFile,readNcdfDim
  
  IMPLICIT none
  
  PRIVATE
  
  PUBLIC :: protectNonStdRadfiles
  
  CONTAINS

    ! -------------------------------------------------------------------------
    ! Prevent overwriting existing files, unless they are in the standard rad
    ! file format used by the AER remote sensing group for simulated data
    ! -------------------------------------------------------------------------
    SUBROUTINE protectNonStdRadfiles(fileName)
    ! Arguments
    character(len=*), intent(in) :: fileName
    ! Local variables
    logical :: fileExists
    integer :: ncid
    integer :: nPts

    inquire(FILE=trim(fileName),EXIST=fileExists)
    if (fileExists) then
      call openNcdfFile(ncid, trim(fileName))
      nPts=readNcdfDim(ncid,'nPts',silent=.true.)
      call closeNcdfFile(ncid)
      if (nPts == 0) then
         print *,'err[IOutilities::protectNonStdRadfiles]:'
         print *,'  Rad output file ',trim(fileName)
         print *,'  exists and is not in AER simulation file format.'
         print *,'  Stopping without overwriting'
         call exit(1)
      endif
    endif

    END SUBROUTINE protectNonStdRadfiles

END MODULE IOutilities
