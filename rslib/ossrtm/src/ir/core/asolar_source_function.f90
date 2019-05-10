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

!--------------------------------------------------------------
!
!  MODULE ASOLAR_SOURCE_FUNCTION: contains items needed to include
!     the sun in OSS Forward Model calculations. AER Inc. 2004.
!
!--------------------------------------------------------------
MODULE asolar_source_function

! <f90Module>***********************************************************
!
! NAME:
!
!   asolar_source_function
!
! PURPOSE:
!
!   Contains items needed to include the sun in OSS Forward Model calculations.
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

  IMPLICIT NONE
  PRIVATE
  !------------------------------------------------------------------------
  !	Public items made available to calling program(s)
  !------------------------------------------------------------------------
  PUBLIC :: GetSolar, InterpSolar, destroySolar
  !------------------------------------------------------------------------
  REAL,    DIMENSION(:), ALLOCATABLE  :: sunrad,sunwvn
  INTEGER                             :: nptsSun

  ! loadDone implementation provides the option to call GetSolar() in advance
  logical, save                       :: loadDone=.false.

CONTAINS

SUBROUTINE GetSolar(Lsol,solFile)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   GetSolar
!
! PURPOSE:
!
!   Read from file for obtaining solar source function.
!
! SYNTAX:
!
!   CALL GetSolar(Lsol, solFile)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   Lsol     INTEGER  file unit number
!   solFile  CHAR     solar file path
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

!c***********************************************************************
!c* Function name: getsolar
!c* Purpose: obtain solar source function
!c* Usage: call getsolar
!c*
!c* Description: read from file
!c*
!c* Inputs:
!c* Var_name     Type       Description
!c* --------     ----       -----------
!c* Lsol         integer    file unit number
!c* solFile      character  solar file
!c*
!c* Outputs: (kept in module)
!c* Var_name     Type       Description
!c* --------     ----       -----------
!c* sunrad       real       solar radiance spectrum
!c* sunwvn       real       wavenumber grid for sunrad
!c* nptsSun      integer    number of points for solar spectrum
!c*
!c* Common blocks: defined the includes
!c* Includes:
!c* Name        Description
!c* ----        -----------
!c* Externals: none
!c*
!c*<copyright_OSS>=====================================================
!c*
!c* Developed by Atmospheric and Environmental Research, Inc.
!c*
!c* Copyright: Atmospheric and Environmental Inc., 1997-2007
!c*                 All Rights Reserved
!c*
!c====================================================</copyright_OSS>
!c***********************************************************************
      INTEGER,             INTENT(IN)     :: Lsol
      CHARACTER(len=*),    INTENT(IN)     :: solFile
      !---Local variables
      CHARACTER*100        :: filename
      CHARACTER*10         :: header  ! change format 10 if this changes
      INTEGER              :: ipt

      if(loadDone) return

      nptsSun=0

! open solar file
      filename=TRIM(solFile)
      open(Lsol,file=filename,status='old',err=200)

! read file if it exists
      read(Lsol,*) nptsSun
      read(Lsol,10) header
      read(Lsol,10) header
10    format(a10)

      ALLOCATE(sunwvn(nptsSun),sunrad(nptsSun))

      do 100 ipt=1,nptsSun
         read(Lsol,*) sunwvn(ipt),sunrad(ipt)
 100  continue

      sunrad=sunrad*1.0e+7  ! convert from (W cm-2 / cm-1) to (mW m-2 / cm-1)

      close(Lsol)

      loadDone=.true.
      return

! jump here if error on open
 200  continue
      write(*,*) 'solar function not used'
      return

END SUBROUTINE getSolar

!***********************************************************************
SUBROUTINE InterpSolar(nptout,wvnout,sunout)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   InterpSolar
!
! PURPOSE:
!
!   Interpolate solar source data from input grid to output grid.
!
! SYNTAX:
!
!   CALL InterpSolar(nptout, wvnout, sunout)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   nptout  INTEGER  Number of output grid
!   wvnout  REAL*8   Wavenumbers on output grid
!
!   INPUTS/OUTPUTS:
!
!   sunout  REAL     Interpolated solar contribution on output grid
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

!c***********************************************************************
!c* Function name: InterpSolar
!c* Purpose: Interpolation routine
!c* Usage: call InterpSolar(nptout,wvnout,sunout)
!c*
!c* Description: Interpolate sunin from input grid to output grid
!c*              Set sunout to zero if nptsSun equal zero
!c*
!c* Inputs:
!c* Var_name     Type       Description
!c* --------     ----       -----------
!c* sunin        real       input data (already in module)
!c* nptin        integer    dimension of input (already in module)
!c* nptout       integer    dimension of input
!c*
!c* Outputs:
!c* Var_name     Type       Description
!c* --------     ----       -----------
!c* sunout       real       input data
!c*
!c* Common blocks: defined the includes
!c* Includes:
!c* Name        Description
!c* ----        -----------
!c* Externals: none
!c*
!c*<copyright_OSS>=====================================================
!c*
!c* Developed by Atmospheric and Environmental Research, Inc.
!c*
!c* Copyright: Atmospheric and Environmental Inc., 1997-2007
!c*                 All Rights Reserved
!c*
!c====================================================</copyright_OSS>
!c***********************************************************************
  INTEGER,                 INTENT(IN)     :: nptout
  REAL*8,  DIMENSION(:),   INTENT(IN)     :: wvnout
  REAL,    DIMENSION(:),   INTENT(INOUT)  :: sunout
  !---Local variables
  INTEGER                  :: i,k

  if (nptsSun.eq.0) then
     write(*,*) 'interpsolar: sunout=0.'
     sunout=0.
     return
  endif

  do i=1,nptout
     if (wvnout(i).le.sunwvn(1))then
        sunout(i)=sunrad(1)+(sunrad(2)-sunrad(1))*(wvnout(i)-sunwvn(1))/ &
             (sunwvn(2)-sunwvn(1))
     else if (wvnout(i).ge.sunwvn(nptsSun)) then
        sunout(i)=sunrad(nptsSun-1)+(sunrad(nptsSun)-sunrad(nptsSun-1))* &
             (wvnout(i)-sunwvn(nptsSun-1)) &
             /(sunwvn(nptsSun)-sunwvn(nptsSun-1))
     else
        do k=1,nptsSun-1
           if (wvnout(i).ge.sunwvn(k).and. &
                wvnout(i).lt.sunwvn(k+1))then
              sunout(i)=sunrad(k)+ &
                   (sunrad(k+1)-sunrad(k))*(wvnout(i)-sunwvn(k))/ &
                   (sunwvn(k+1)-sunwvn(k))
           endif
        end do
     end if
  end do
  return
END SUBROUTINE InterpSolar
!c***********************************************************************

subroutine destroySolar()
    if (allocated(sunrad))       deallocate(sunrad)
    if (allocated(sunwvn))       deallocate(sunwvn)
    loadDone=.false.
end subroutine destroySolar

END MODULE asolar_source_function

