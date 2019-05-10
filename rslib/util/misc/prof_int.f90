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

!                    Input    Input Input  Inputs Output   Inputs
!                      |        |     |    |    |  |   |   |     |
SUBROUTINE Prof_Int( nLev_In, p_In, x_In, nLev, p, x, ip, xMin, xMax )
!
! Inputs:
! Variable   Type   Contents
! --------   ----   --------
! nLev_In    INT    # of Levels in Grid "p_In"
! p_In       REAL   Input Pressure Grid on Which "x_In" is Defined
! x_In       REAL   Input Atmospheric Profile to Use in Interpolating
!
! nLev       INT    # of Gridpoints in the Output Grid "p"
! p          REAL   Valid Pressure at the Output Gridpoints
!
! xMin       REAL   Minimum x Value for Levels Above p_In(1)
! xMax       REAL   Maximum x Value for Levels Above p_In(1)
!
! Output:
! Variable   Type   Contents
! --------   ----   --------
! x          REAL   Interpolated Profile, Valid on Grid "p"
! ip         INT    # Levels Interpolated TO
!
! Usage Notes: 1. The Pressure Grids "p" and "p_In" Must be 
!                 Monotonically Increasing.
!              2. This Routine Extrapolates Above the Highest
!                 "p_In" Levels, If Necessary, But Does NOT
!                 Extrapolate Below pSfc.
!
  IMPLICIT None

  INTEGER :: nLev_In, nLev, ip, im

  REAL(Kind=4) :: p_In(nLev_In), x_In(nLev_In)
  REAL(Kind=4) :: p(nLev), x(nLev)
  REAL(kind=4) :: xMin, xMax

  x(1:nLev) = -999.  ! To Flag Points Below pSfc
!
! Interpolate the atmospheric profile "x_In," valid on
! the input pressure grid "p_In," and return the
! interpolated results of "x" at points in "p"
!
  ip = 1  ! Index for "p" (the Output Pressure Grid)
  im = 2  ! Index for the Level Just Below "ip"
!
! Starting with the Lowest Pressure (at Level ip=1), 
! Continue on As Long As the the Output Pressure Gridpoint
! at Level "ip" is Less Than the Input Surface Pressure
!
  DO WHILE( p(ip) < p_In(nLev_In) .AND. ip <= nLev )
!
! Find the Lower-Bounding Level "im" for p(ip)
!
     DO WHILE( p_In(im) < p(ip) .AND. im < nLev_In )
        im = im+1
     END DO
!
! Boundaries for Pressure "p(ip)" Have Been Found.
!
     CALL Lint( p_In(im-1), p_In(im), x_In(im-1), x_In(im), p(ip), x(ip) )
     IF( im == 2 ) x(ip) = Max( Min(x(ip),xMax), xMin )

     ip = ip+1

  END DO

  IF( ip <= nLev ) THEN
     CALL Lint( p_In(nLev_In), p(ip-1), x_In(nLev_In), x(ip-1), p(ip), x(ip) )
  ELSE
     ip = ip-1
  END IF

  RETURN

END SUBROUTINE Prof_Int

SUBROUTINE Lint( p_In1, p_In2, x_In1, x_In2, p, x )  ! Check your Navel
!
! Inputs
! ------
! p_In1 - Pressure Grid Points to Interpolate Between
! p_In2
!
! x_In1 - Profile Data to Interpolate With
! x_In2
!
! p - Pressure Level to Interpolate AT
!
! Output: x - Interpolated Profile Value at Level "p"
! ------
!
  REAL(Kind=4) :: p_In1, p_In2, x_In1, x_In2, p, x, aa

! Add special handling when x_In1 or x_In2 eq 0 (e.g. cloud liq water)
  if(x_In1 .eq. 0.)then
     if(x_In2 .eq. 0.)then
        x=0.
        return
     else
        x_In1=1.E-20
     endif
  elseif (x_In2 .eq. 0.)then
        x_In2=1.E-20
  endif
  aa = LOG( x_In2/x_In1 ) / LOG( p_In2/p_In1 )
  x = x_In1 * (p/p_In1)**aa

  RETURN
END SUBROUTINE Lint

SUBROUTINE Lint0( p_In1, p_In2, x_In1, x_In2, p, x, iatm )  ! Check your Navel
!
! Inputs
! ------
! p_In1 - Pressure Grid Points to Interpolate Between
! p_In2
!
! x_In1 - Profile Data to Interpolate With
! x_In2
!
! p - Pressure Level to Interpolate AT
!
! Output: x - Interpolated Profile Value at Level "p"
! ------
!
  REAL(Kind=4) :: p_In1, p_In2, x_In1, x_In2, p, x, aa
  integer :: iatm

! Add special handling when x_In1 or x_In2 eq 0 (e.g. cloud liq water)
  if(x_In1 .eq. 0.)then
     if(x_In2 .eq. 0.)then
        x=0.
        return
     else
        x_In1=1.E-20
     endif
  elseif (x_In2 .eq. 0.)then
        x_In2=1.E-20
  endif
  aa = LOG( x_In2/x_In1 ) / LOG( p_In2/p_In1 )
  if(iatm .ge. 10078 .and. iatm .le. 10080)then
     print *,'LOG( x_In2/x_In1 ), LOG( p_In2/p_In1 ),aa,(p/p_In1)**aa,x_In1 * (p/p_In1)**aa,=',&
          LOG( x_In2/x_In1 ), LOG( p_In2/p_In1 ),aa,(p/p_In1)**aa,x_In1 * (p/p_In1)**aa
  endif
  x = x_In1 * (p/p_In1)**aa

  RETURN
END SUBROUTINE Lint0

