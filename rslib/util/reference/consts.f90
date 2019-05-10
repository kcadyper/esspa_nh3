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
!  MODULE CONSTANTS: It contains physical constants needed acrosss a number
!                    of applications. AER Inc. 2004.
!                    Using this module will enforce consistency.
!                    To use it, type:
!                    USE constants, ONLY : list_of_constants_needed
!                    in the program/subroutine that needs the constants.
!
!  NOTE: - All constants are declared public in this module.
!        - The use of the 'ONLY' keyword with the USE statement, will
!          help avoid having naming conflicts.
!
!------------------------------------------------------------------------------
MODULE constants

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
  !                   FUNDAMENTAL PHYSICAL CONSTANTS
  !----------------------------------------------------------------------------
  !TYPE         | Name     | Value        | Description             |    Units
  !----------------------------------------------------------------------------
  !     constants were updated based on
  !     LBLRTM   $Revision: 27087 $ phys_consts.f90
  !     Constants from NIST May 2010
  ! OLD version
  ! REAL,PARAMETER::Plnk_cgs = 6.626176e-27 !Planck constant          | [g.cm^2/s]
  ! REAL,PARAMETER::LtSp_cgs = 2.997925e+10 !Speed of Light           | [cm/s]
  ! REAL,PARAMETER::Bltz_cgs = 1.380662e-16 !Boltzman Constant        | [g.cm^2/(s^2.K)]

  REAL,PARAMETER::Plnk_cgs = 6.62606876E-27!Planck constant         | [g.cm^2/s]
  REAL,PARAMETER::LtSp_cgs = 2.99792458E+10!Speed of Light          | [cm/s]
  REAL,PARAMETER::Bltz_cgs = 1.3806503e-16!Boltzman Constant        | [g.cm^2/(s^2.K)]
  REAL,PARAMETER::xNa      = 6.0238E+23   !Avogadro's Number        | [1/(mol.g.mole)]
  REAL,PARAMETER::gravCnst = 6.672e-11    !Gravity constant         | [N.m^2/kg^2]
  REAL,PARAMETER::Rd       = 287.         !Dry Air Gas const.       | [J/(kg.K)]
  REAL,PARAMETER::cp       = 1004.        !Spec.heat of d.air @P=ct | [J/(kg.K)]
  REAL,PARAMETER::cv       = 717.         !Spec.heat of d.air @V=ct | [J/(kg.K)]
  REAL,PARAMETER::T0       = 273.15       !Standard Temperature     | [K]
  REAL,PARAMETER::P0       = 1013.        !Standard Pressure        | [mb]
  !----------------------------------------------------------------------------
  !                       MATHEMATICAL CONSTANTS
  !----------------------------------------------------------------------------
  !TYPE         | Name     | Value        | Description             |   Units
  !----------------------------------------------------------------------------
  REAL,PARAMETER::PI       = 3.141592654  !PI number                | [rad]
  REAL,PARAMETER::rad2deg  = 180./PI      !Radian 2 Degree factor   | [deg/rad]
  REAL,PARAMETER::deg2rad  = PI/180.      !Degree 2 Radian factor   | [rad/deg]
  !----------------------------------------------------------------------------
  !                       GEOPHYSICAL CONSTANTS
  !----------------------------------------------------------------------------
  !TYPE         | Name     | Value        | Description             |   Units
  !----------------------------------------------------------------------------
  REAL,PARAMETER::g0       = 9.81         !Gravity @ sea level      | [m/s^2]
  REAL,PARAMETER::wvmwt    = 18.016       !Water vapor molec weight | [g]
  REAL,PARAMETER::drymwt   = 28.964       !Dry Air molecular weight | [g]
  REAL,PARAMETER::ratioMwt = wvmwt/drymwt !used as 0.622,0.621      | [-]
  REAL,PARAMETER::ErthMass = 5.9763e24    !Mass of the earth        | [kg]
  REAL,PARAMETER::ErthRad  = 6371.2       !Earth radius             | [km]
  REAL,PARAMETER::CosmBckg = 2.73         !Cosmic Background Tb     | [K]
  !----------------------------------------------------------------------------
  !                       MAGIC NUMBERS CONSTANTS
  !----------------------------------------------------------------------------
  !TYPE                 | Name       | Value | Description            |  Units
  !----------------------------------------------------------------------------
  CHARACTER, PARAMETER::MISSING_CHAR = '*'   !For missing characters  |   -
  INTEGER,   PARAMETER::MISSING_INT  = -999  !For missing integers    |   -
  REAL,      PARAMETER::MISSING_REAL = -999. !For missing reals       |   -
  !----------------------------------------------------------------------------
END MODULE constants
