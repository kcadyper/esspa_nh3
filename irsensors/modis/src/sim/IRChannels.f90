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

MODULE IRChannels

! <f90Module>***********************************************************
!
! NAME:
!
!   IRChannels
!
! PURPOSE:
!
!   Manipulate settings for channels.
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>


  implicit none
  private
  public :: initIRChannels

CONTAINS

  subroutine initIRChannels(nchanmw,nchan,wvn,kchan)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   initIRChannels
!
! PURPOSE:
!
!   Turn on/off channels.
!
! SYNTAX:
!
!   CALL initIRChannels(nchanmw, nchan, wvn, kchan)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   nchanmw  INTEGER  Number of MW channels
!   nchan    INTEGER  Number of channels
!   wvn      REAL     Wavenumbers
!
!   INPUTS/OUTPUTS:
!
!   kchan    LOGICAL  Channel on/off mask
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    USE IRReadStdInputs
    implicit NONE
    !--I/O variables
    integer,               intent(in)    :: nchanmw,nchan
    real,    dimension(:), intent(in)    :: wvn
    LOGICAL, dimension(:), intent(inout) :: kchan
    !--Local variables
    real    :: flocal,dum
    integer :: i

    ! This subroutine is used for turning on or off channels
    !   falling within various spectral ranges.
    ! For small channel sets, the user may simply edit the
    !   kchan vector in the input file, rather than using this
    !   subroutine.  It that case, this becomes a dummy subroutine
    !   that remains present to maintain compatibility.

    return

  end subroutine initIRChannels

END MODULE IRChannels
