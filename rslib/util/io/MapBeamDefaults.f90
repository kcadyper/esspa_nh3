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

MODULE MapBeamDefaults

  implicit none
  private
  public :: mapBeamDefNew

CONTAINS
      
  subroutine mapBeamDefNew(nchan,nBeam,mapBeam,nBeamloc,mapBeamloc)

! Check inputs and set defaults for nBeam and mapBeam for opening a new file
! nBeam   = number of beams, at any given scan position
! mapBeam = tells which beam each channel is on
! Local versions of nBeam and mapBeam are used because cannot write to an
! optional argument, which has no assigned memory if not present
    
    !---in/out variables
    integer,               intent(in)              :: nchan
    integer,               intent(in),    optional :: nBeam
    integer, dimension(:), intent(in),    optional :: mapBeam
    integer,               intent(inout)           :: nBeamloc
    integer, dimension(:), pointer                 :: mapBeamloc

    !---local variables
    integer                                        :: i

    if (present(nBeam)) then
       nBeamloc=nBeam
    else
       nBeamloc=1                        ! assume all channels on same beam*
    end if
    if (present(mapBeam)) then
       if (.not. present(nBeam)) then
          print *,'err[MapBeamDefaults::mapBeamDefNew]: ', &
               ' User provided mapBeam, but not nBeam'
          call errorHalt(1)
       end if
       if (any(mapBeam(1:nchan) < 1) .or. any(mapBeam(1:nchan) > nBeamloc)) &
          then
          print *,'err[MapBeamDefaults::mapBeamDefNew]: ', &
               ' mapBeam out of range; nBeam,min(mapBeam),max(mapBeam):', &
               nBeamloc,minval(mapBeam(1:nchan)),maxval(mapBeam(1:nchan))
          call errorHalt(1)
       end if
       allocate (mapBeamloc(size(mapBeam)))
       mapBeamloc=mapBeam
    else   ! mapBeam not present
       if (.not. present(nBeam)) then    ! assume all channels on same beam*
          allocate (mapBeamloc(nchan))
          mapBeamloc(1:nchan)=1
       elseif (nBeam == 1) then          ! assume all channels on same beam
          allocate (mapBeamloc(nchan))
          mapBeamloc(1:nchan)=1
       elseif (nBeam == nchan) then      ! assume each channel has own beam
          allocate (mapBeamloc(nchan))
          mapBeamloc(1:nchan)=(/(i,i=1,nchan)/)
       else
          print *,'err[MapBeamDefaults::mapBeamDefNew]: ', &
           ' mapBeam not given and no default consistent with nBeam',nBeam
          call errorHalt(1)
       end if
    end if

    return

  end subroutine mapBeamDefNew

END MODULE MapBeamDefaults
