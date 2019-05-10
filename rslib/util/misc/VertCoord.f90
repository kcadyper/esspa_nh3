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
MODULE VertCoord

  IMPLICIT NONE

  PRIVATE

  !---------------------------------------------------------------
  !     Vertical grid types
  !---------------------------------------------------------------
  INTEGER,                   PARAMETER, PUBLIC :: mxCoordTyp=12
  CHARACTER(LEN=mxCoordTyp), PARAMETER, PUBLIC :: Pcoord='pressureStd'
  CHARACTER(LEN=mxCoordTyp), PARAMETER, PUBLIC :: Scoord='sigma'
  CHARACTER(LEN=mxCoordTyp), PARAMETER, PUBLIC :: Hcoord='hybrid-sigma'

  !---------------------------------------------------------------
  !  Public subroutines (accessible from outside module) 
  !---------------------------------------------------------------
  PUBLIC :: putSigmDefin,putHybrDefin,getSigmPres,getHybrPres,destroyVertCoord
  PUBLIC :: getPresByType

  !---------------------------------------------------------------
  !     Declarations of private data (private to the module)
  !---------------------------------------------------------------
  INTEGER                         :: myNlevel
  REAL                            :: mySigPtop
  REAL, DIMENSION(:), ALLOCATABLE :: myVcoordA,myVcoordB

CONTAINS

  SUBROUTINE putSigmDefin(nlevel,sigPtop,vCoordA)
    INTEGER,            INTENT(IN) :: nlevel
    REAL,               INTENT(IN) :: sigPtop
    REAL, DIMENSION(:), INTENT(IN) :: vCoordA

    if (.not. allocated(myVcoordA)) ALLOCATE (myVcoordA(nlevel))
    myNlevel=nlevel
    mySigPtop=sigPtop
    myVcoordA=vCoordA

  END SUBROUTINE putSigmDefin

  SUBROUTINE putHybrDefin(nlevel,vCoordA,vCoordB)
    INTEGER,            INTENT(IN) :: nlevel
    REAL, DIMENSION(:), INTENT(IN) :: vCoordA
    REAL, DIMENSION(:), INTENT(IN) :: vCoordB

    if (.not. allocated(myVcoordA)) ALLOCATE (myVcoordA(nlevel))
    if (.not. allocated(myVcoordB)) ALLOCATE (myVcoordB(nlevel))
    myNlevel=nlevel
    myVcoordA=vCoordA
    myVcoordB=vCoordB

  END SUBROUTINE putHybrDefin

  FUNCTION getSigmPres(psfc)
    REAL,  INTENT(IN)         :: psfc
    REAL, DIMENSION(myNlevel) :: getSigmPres

    getSigmPres=mySigPtop+myVcoordA*(psfc-mySigPtop)

  END FUNCTION getSigmPres

  FUNCTION getHybrPres(psfc)
    REAL,  INTENT(IN)         :: psfc
    REAL, DIMENSION(myNlevel) :: getHybrPres

    getHybrPres=myVcoordA+myVcoordB*psfc

  END FUNCTION getHybrPres

  SUBROUTINE destroyVertCoord()

    IF (ALLOCATED(myVcoordA)) DEALLOCATE(myVcoordA)
    IF (ALLOCATED(myVcoordB)) DEALLOCATE(myVcoordB)

  END SUBROUTINE destroyVertCoord

  FUNCTION getPresByType(vCoordTyp,nLev,pSfc, &
                         press,pTop,coordA,coordB)

    character(len=mxCoordTyp),  intent(in) :: vCoordTyp
    integer,                    intent(in) :: nLev
    real,                       intent(in) :: pSfc
    real, dimension(:),         intent(in) :: press
    real,                       intent(in) :: pTop
    real, dimension(:),         intent(in) :: coordA
    real, dimension(:),         intent(in) :: coordB
    real, dimension(nLev)                  :: getPresByType

    SELECT CASE (TRIM(vCoordTyp))
    CASE (Pcoord)
      getPresByType(1:nLev)=press(1:nLev)
    CASE (Scoord)
      getPresByType(1:nLev)=pTop+coordA(1:nLev)*(pSfc-pTop)
    CASE (Hcoord)
      getPresByType(1:nLev)=coordA(1:nLev)+coordB(1:nLev)*pSfc
    CASE DEFAULT
      PRINT *,'err[VertCoord::getPresByType]: Unrecognized coordinate type: ', &
         TRIM(vCoordTyp)
      CALL EXIT(1)
    END SELECT

  END FUNCTION getPresByType

END MODULE VertCoord
