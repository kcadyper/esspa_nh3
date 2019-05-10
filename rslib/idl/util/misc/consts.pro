;------------------------------------------------------------------------------
;
; CONSTANTS FILE: This file is equivalent to the F90 version consts.f90.
;                 It contains constants needed acrosss a number of IDL 
;                 applications. AER Inc. 2004.
;                 Using this file will enforce consistency.
;                 To use it, type:
;                 @consts.pro
;                 in the program/subroutine that needs the constants.
;
;------------------------------------------------------------------------------

;----------------------------------------------------------------------------
;                   FUNDAMENTAL PHYSICAL CONSTANTS
;----------------------------------------------------------------------------
;Name    | Value        | Description             |    Units
;----------------------------------------------------------------------------
Plnk_cgs = 6.626176D-27 ;Planck constant          | [g.cm^2/s]
LtSp_cgs = 2.997925D+10 ;Speed of Light           | [cm/s]
Bltz_cgs = 1.380662D-16 ;Boltzman Constant        | [g.cm^2/(s^2.K)]
xNa      = 6.0238D+23   ;Avogadro's Number        | [1/(mol.g.mole)]
gravCnst = 6.672D-11    ;Gravity constant         | [N.m^2/kg^2]
Rd       = 287.         ;Dry Air Gas const.       | [J/(kg.K)]
cp       = 1004.        ;Spec.heat of d.air @P=ct | [J/(kg.K)]
cv       = 717.         ;Spec.heat of d.air @V=ct | [J/(kg.K)]
T0       = 273.15       ;Standard Temperature     | [K]
P0       = 1013.        ;Standard Pressure        | [mb]
;----------------------------------------------------------------------------
;                       MATHEMATICAL CONSTANTS
;----------------------------------------------------------------------------
;Name    | Value        | Description             |   Units
;----------------------------------------------------------------------------
PI       = 3.141592654D0;PI number                | [rad]
rad2deg  = 180.D0/PI    ;Radian 2 Degree factor   | [deg/rad]
deg2rad  = PI/180.D0    ;Degree 2 Radian factor   | [rad/deg]
;----------------------------------------------------------------------------
;                       GEOPHYSICAL CONSTANTS
;----------------------------------------------------------------------------
;Name    | Value        | Description             |   Units
;----------------------------------------------------------------------------
g0       = 9.81         ;Gravity @ sea level      | [m/s^2]
wvmwt    = 18.016       ;Water vapor molec weight | [g]
drymwt   = 28.964       ;Dry Air molecular weight | [g]
ratioMwt = wvmwt/drymwt ;used as 0.622,0.621      | [-]
ErthMass = 5.9763D24    ;Mass of the earth        | [kg]
ErthRad  = 6371.2       ;Earth radius             | [km]
CosmBckg = 2.73         ;Cosmic Background Tb     | [K]
;----------------------------------------------------------------------------
;                       MAGIC NUMBERS CONSTANTS
;----------------------------------------------------------------------------
;Name         | Value   | Description             |   Units
;----------------------------------------------------------------------------
MISSING       = -999    ;For missing values       |   -
MISSING_CHAR  = '*'     ;For missing characters   |   -
;----------------------------------------------------------------------------
