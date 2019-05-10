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

!-------------------------------------------------------------------------------
!
!  MODULE LvlInterp: contains subroutines related to interpolation over
!                    grid levels
!
!     Copyright AER, Inc., 2008. All rights Reserved.
!-------------------------------------------------------------------------------
MODULE LvlInterp

! <f90Module>***********************************************************
!
! NAME:
!
!   LvlInterp
!
! PURPOSE:
!
!   Contains subroutines related to interpolation over grid levels
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>


  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------
  !  List of Public subroutines (accessible from outside module) 
  !---------------------------------------------------------------
  PUBLIC :: lvl_int

  !---------------------------------------------------------------
  !     Declarations of private data (private to the module)
  !---------------------------------------------------------------
        
CONTAINS

  SUBROUTINE lvl_int(xInp,pGrid,nxdim,plvl,xlvl,dxlvldp,ilvl,dxlvldxu,dxlvldxl)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   lvl_int
!
! PURPOSE:
!
!   performs logarithmic interpolation of physical quantities within pressure 
!   layers
!
! SYNTAX:
!
!   CALL lvl_int(xInp, pGrid, nxdim, plvl, xlvl, dxlvldp, ilvl, 
!      dxlvldxu, dxlvldxl)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   xInp      REAL     Input data
!   pGrid     REAL     Pressure grid
!   nxdim     INTEGER  Array dimension
!   plvl      REAL     Pressure at which interpolation to be done
!   
!   INPUTS/OUTPUTS:
!   
!   xlvl      REAL     interpolated value
!   dxlvldp*  REAL     derivative of xlvl wrt plvl
!   ilvl*     INTEGER  index of layer used in interpolation
!   dxlvldxu* REAL     derivative of xlvl wrt upper level x
!   dxlvldxl* REAL     derivative of xlvl wrt upper level x
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !---Input variables
    REAL,    DIMENSION(*), INTENT(IN)              :: xInp,pGrid
    INTEGER,               INTENT(IN)              :: nxdim
    REAL,                  INTENT(IN)              :: plvl
    !---Output variables
    REAL,                  INTENT(INOUT)           :: xlvl
    REAL,                  INTENT(INOUT), OPTIONAL :: dxlvldp
    INTEGER,               INTENT(INOUT), OPTIONAL :: ilvl
    REAL,                  INTENT(INOUT), OPTIONAL :: dxlvldxu,dxlvldxl
    !---Local constants
    REAL,    PARAMETER :: dp=.001
    INTEGER, PARAMETER :: iMissing=-1
    !---Local variables
    REAL               :: xx
    REAL               :: prlog,psexp,rlog
    INTEGER            :: iu,il,i
    
    iu=iMissing
    LevLoop: DO i=nxdim-1,1,-1
       IF (pGrid(i).LT.plvl-dp) THEN
          iu=i
          EXIT
       ENDIF
    END DO LevLoop
    IF (iu == iMissing) THEN
       PRINT *,'err[LvlInterp::lvl_int]: input pressure out of bounds'
       PRINT *, pGrid(1:nxdim-1), "...",plvl
       call errorHalt(1)
    ENDIF
    il=iu+1
    if (present(ilvl)) ilvl=il
    prlog=log(pGrid(iu)/pGrid(il))
    xx=log(xInp(iu)/xInp(il))/prlog
    psexp=(plvl/pGrid(il))**xx
    xlvl=xInp(il)*psexp
    if (present(dxlvldp)) dxlvldp=xx*xlvl/plvl
  
    !Derivatives with respect to xInp(iu) and xInp(il)
    if (present(dxlvldxu) .and. present(dxlvldxl)) then
       rlog=log(plvl/pGrid(il))/prlog
       dxlvldxu=(xInp(il)/xInp(iu)) * psexp * rlog
       dxlvldxl=psexp*(real(1)-rlog)
    endif
  
    RETURN
  END SUBROUTINE lvl_int
  
END MODULE LvlInterp
