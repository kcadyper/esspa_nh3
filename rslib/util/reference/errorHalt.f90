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

! Halt execution when an error occurs, while setting the process status
! to an integer value.  
! This procedure provides a means for a script that is executing a 
! FORTRAN program to distinguish between its handling of fatal error 
! conditions and normal program termination.

! The statement "call exit(i)" is not ISO standard.  It is used because 
! there is no ISO standard method for FORTAN programs to set a termination
! status.  If this code is applied on a platform that does not accept this
! statement, it can be replaced by a platform-specific version or by the 
! standard "stop", which does not generally set a termination status other
! than 0.

subroutine errorHalt(iStatus)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   errorHalt
!
! PURPOSE:
!
!   Halt execution when an error occurs, while setting the process status to an 
!   integer value.
!
! SYNTAX:
!
!   CALL errorHalt(iStatus)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   iStatus  INTEGER  Termination status code
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>


  implicit none

  integer, intent(in) :: iStatus    ! Termination status code

  call exit(iStatus)

end subroutine errorHalt
