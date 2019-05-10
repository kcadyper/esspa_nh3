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
!  MODULE SensorData: This module provides generic interface to 
!         sensor-specific functionality.
!         This version does nothing but provide a dummy 
!         interface, and can be used with sensors for which we 
!         are not using dynamic noise.
!         For sensors that use dynamic noise, a sensor-specific 
!         version can be stored in the sensor's src directory 
!         and made locally to take precedence over this library 
!         version.  In this case, the library versions of any 
!         modules that use SensorData must also be made locally 
!         after SensorData is compiled.
!         AER Inc. 2007.
!--------------------------------------------------------------
MODULE SensorData
  IMPLICIT NONE
  PRIVATE
  !------------------------------------------------------------------------
  !	Public items made available to calling program(s)   
  !------------------------------------------------------------------------
  PUBLIC :: setDynNoise, queryDynNoise, dynNoise, queryBeamMap,&
       querySpillover, destroySensorData,var_nedt
  INTEGER, PARAMETER, PUBLIC :: LEN_ID_SensorData=1  ! length of ID strings

  !------------------------------------------------------------------------
  !     Declarations of private data (private to the module)
  !------------------------------------------------------------------------
  INTEGER, PARAMETER                    :: lID=LEN_ID_SensorData  
CONTAINS
  
  SUBROUTINE setDynNoise(iun,nedt_file)
    INTEGER,          INTENT(IN) :: iun
    CHARACTER(LEN=*), INTENT(IN) :: nedt_file
    RETURN
  END SUBROUTINE setDynNoise

  !----------------------------------------------------------------------------

  SUBROUTINE queryDynNoise(nchan,chanID)
    !---In/Out variables
    INTEGER                                :: nchan
    CHARACTER(lID), DIMENSION(:), POINTER  :: chanID
   RETURN
  END SUBROUTINE queryDynNoise

  !----------------------------------------------------------------------------

  FUNCTION dynNoise(ispec,tAntn,icaset)
    REAL       :: dynNoise
    !---In/Out variables
    INTEGER, INTENT(IN)    :: ispec,icaset
    REAL,    INTENT(IN)    :: tAntn
    RETURN
  END FUNCTION dynNoise

  !----------------------------------------------------------------------------
  
  SUBROUTINE queryBeamMap(iu,fname,chanID,beamID,iBeam)
    !---In/Out variables
    INTEGER,                      INTENT(IN)  :: iu
    CHARACTER(LEN=*),             INTENT(IN)  :: fname
    CHARACTER(lID), DIMENSION(:), POINTER     :: chanID
    CHARACTER(lID), DIMENSION(:), POINTER     :: beamID  
    INTEGER, DIMENSION(:),        POINTER     :: iBeam
    RETURN
  END SUBROUTINE queryBeamMap

  !----------------------------------------------------------------------------

  SUBROUTINE querySpillover(iu,fname,beamID,spoev)
    !---In/Out variables
    INTEGER,                      INTENT(IN) :: iu
    CHARACTER(LEN=*),             INTENT(IN) :: fname
    CHARACTER(lID), DIMENSION(:), POINTER    :: beamID  
    REAL,           DIMENSION(:), POINTER    :: spoev
    RETURN
  END SUBROUTINE querySpillover

  !----------------------------------------------------------------------------
  FUNCTION var_nedt(ispec,tAntn,time_set,bw_set,ica_set,icnvBW)
    REAL :: var_nedt
    !---In/Out variables
    INTEGER                       :: ispec,ica_set,icnvBW
    REAL                          :: tAntn,time_set,bw_set
    RETURN
  END FUNCTION var_nedt

  !----------------------------------------------------------------------------

  SUBROUTINE destroySensorData()  
    RETURN
  END SUBROUTINE destroySensorData

  !----------------------------------------------------------------------------

END MODULE SensorData
