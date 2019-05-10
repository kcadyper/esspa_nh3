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

c-------------------------------------------------------------------------------
c     $Name$ 
c     $Id$ 
c     Copyright AER, Inc., 2002, 2003. All rights Reserved.
c-------------------------------------------------------------------------------

      subroutine classatm(Ym0,xlat,xlon,time,igeo,iclassatm)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   classatm
!
! PURPOSE:
!
!   Set atmosphere class.
!
! SYNTAX:
!
!   CALL classatm(Ym0, xlat, xlon, time, igeo, iclassatm)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   Ym0        REAL     Radiance
!   xlat       REAL     Latitude
!   xlon       REAL     Longitude
!   time       INTEGER  Date/time
!   igeo       INTEGER  Surface (geography) class index
!   
!   INPUTS/OUTPUTS:
!   
!   iclassatm  INTEGER  Atmospheric class index
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

      integer,               intent(in)    :: time(6),igeo
      real,                  intent(in)    :: xlat,xlon
      real,   dimension(*),  intent(in)    :: Ym0
      integer,               intent(inout) :: iclassatm
      !---dummy value
      iclassatm=1
      !iclassatm=igeo

      return
      end


      subroutine classatmRegr(Ym0,xlat,xlon,
     $     time,igeo,scanAng,iclassatmRegr)
      integer,               intent(in)    :: time(6),igeo

!<f90Subroutine>********************************************************
!
! NAME:
!
!   classatmRegr
!
! PURPOSE:
!
!   Set atmosphere class, when using regression.
!
! SYNTAX:
!
!   CALL classatmRegr(Ym0, xlat, xlon, time, igeo, scanAng, 
!   iclassatmRegr)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   Ym0            REAL     Radiance
!   xlat           REAL     Latitude
!   xlon           REAL     Longitude
!   time           INTEGER  Date/time
!   igeo           INTEGER  Surface (geography) class index
!   scanAng        REAL     scan Angle
!   
!   INPUTS/OUTPUTS:
!   
!   iclassatmRegr  INTEGER  Atmospheric class index for regression
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

      real,                  intent(in)    :: xlat,xlon,scanAng
      real,   dimension(*),  intent(in)    :: Ym0
      integer,               intent(inout) :: iclassatmRegr

      !---dummy value
      iclassatmRegr=1


      return
      end


