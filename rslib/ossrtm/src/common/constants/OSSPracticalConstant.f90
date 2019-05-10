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

!------------------------------------------------------------------------------
!
!  MODULE CONSTANTS: It contains OSS parameters needed acrosss a number
!                    of applications. AER Inc. 2004.
!                    Using this module will enforce consistency.
!                    To use it, type:
!                    USE OSSPracticalConstant, ONLY : list_of_constants_needed
!                    in the program/subroutine that needs the constants.
!
!  NOTE: - All constants are declared public in this module.
!        - The use of the 'ONLY' keyword with the USE statement, will
!          help avoid having naming conflicts.
!
!------------------------------------------------------------------------------
MODULE OSSPracticalConstant

! <f90Module>***********************************************************
!
! NAME:
!
!   constants
!
! PURPOSE:
!
!   Contains physical constants needed acrosss a number of applications
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

  IMPLICIT NONE

  !----------------------------------------------------------------------------
  !                       MAGIC NUMBERS CONSTANTS
  !----------------------------------------------------------------------------
  !TYPE                 | Name       | Value | Description            |  Units
  !----------------------------------------------------------------------------
  CHARACTER, PARAMETER::MISSING_CHAR = '*'   !For missing characters  |   -
  INTEGER,   PARAMETER::MISSING_INT  = -999  !For missing integers    |   -
  REAL,      PARAMETER::MISSING_REAL = -999. !For missing reals       |   -

END MODULE OSSPracticalConstant
