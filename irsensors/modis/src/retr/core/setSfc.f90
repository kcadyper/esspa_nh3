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

MODULE SurfaceTypingModule

! <f90Module>***********************************************************
!
! NAME:
!
!   SurfaceTypingModule
!
! PURPOSE:
!
!   For surface classification
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

  USE StateIndexModule
  IMPLICIT NONE
  PRIVATE
  !------------------------------------------------------
  !	Public items made available to calling program(s)   
  !------------------------------------------------------
  PUBLIC :: setSfc
  !---------------------------------------------------------
  !     Declarations of private data (private to the module)
  !---------------------------------------------------------
  
CONTAINS

  subroutine setSfc(iStrat,iexit,x,igeo,plandAvg,IG,NG,IemMwG,NEmMwG, &
      plandMaxOc)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   setSfc
!
! PURPOSE:
!
!   Set surface classification based on inputs.
!
! SYNTAX:
!
!   CALL setSfc(iStrat, iexit, x, igeo, plandAvg, IG, NG, IemMwG, 
!      NEmMwG, plandMaxOc)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   iStrat      INTEGER             Stratification level flag
!   iexit       INTEGER             Inicator of state of surface 
!                                   classification 
!   x           REAL                Retrieval state vector
!   plandAvg    REAL                pland representing average
!   IG          TYPE(STATEINDEX_T)  Starting indices for sections 
!                                   of geophysical state vector 
!   NG          TYPE(STATEINDEX_T)  Number of elements for sections 
!                                   of geophysical state vector 
!   IemMwG      INTEGER             Starting index for MW 
!                                   emissivity in geophysical state 
!                                   vector 
!   NEmMwG      INTEGER             Number of Mw emissivities
!   plandMaxOc  REAL                Upper limit on pland for 
!                                   classifiying as water surface 
!   
!   INPUTS/OUTPUTS:
!   
!   igeo        INTEGER             Surface (geography) class index
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !--I/O variables
    integer,            intent(in)    :: iStrat,iexit
    real, dimension(:), intent(in)    :: x
    real,               intent(in)    :: plandAvg
    real,               intent(in)    :: plandMaxOc
    TYPE(StateIndex_t), intent(in)    :: IG, NG
    integer,            intent(in)    :: IemMwG,NEmMwG
    integer,            intent(inout) :: igeo    
    !--Local variables
    real                              :: Ts,demMw

    if(plandAvg < plandMaxOc)then
       igeo=1
    else
       igeo=2
    end if

    return
  end subroutine setSfc
  
end MODULE SurfaceTypingModule

