MODULE CloudParameters
!
! Parameter module that contains all the constants used throughtout
! the cloud process
!
! yhe@aer.com, 12/17/2015
!
  IMPLICIT NONE

  INTEGER, PARAMETER :: mxNameLength = 128
  INTEGER, PARAMETER :: mxMsgLength = 256
  INTEGER, PARAMETER :: mxHydrometeor = 6      !e.g., liquid, ice, rain, snow, graupel, melt
  INTEGER, PARAMETER :: mxPhysicalProperty = 5 !e.g., primary size, secondary size, temperature, density, melt
  INTEGER, PARAMETER :: mxOpticalProperty = 3  !e.g., absoprtion, scattering, asymmetry
  CHARACTER(mxNameLength),DIMENSION(mxOpticalProperty),PARAMETER :: &
       knownOptProp=(/'kabs ','kscat','gasym'/)
  CHARACTER(mxNameLength),DIMENSION(mxOpticalProperty),PARAMETER :: & 
       knownOptUnits=(/'micrometers','NULL       ','K          '/)
  INTEGER, PARAMETER :: nHydNonPcp=2  !fallback number of hydrometeors
  INTEGER, PARAMETER :: mxIntplType=2 !number of interploation types
  

END MODULE CloudParameters
