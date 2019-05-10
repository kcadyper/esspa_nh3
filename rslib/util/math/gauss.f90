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

function gauss(iseed,sdv,mean)

!<f90Function>**********************************************************
!
! NAME:
!
!   gauss
!
! PURPOSE:
!
!   Compute Gaussian distribution.
!
! SYNTAX:
!
!   Results=gauss(iseed, sdv, mean)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   sdv    REAL     Desired standard deviation of the data distribution
!   mean   REAL     Desired mean of the data distribution
!   
!   INPUTS/OUTPUTS:
!   
!   iseed  INTEGER  seed for random number generator
!
!   * OPTIONAL
!
! RETURN:
!
!     REAL     
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>


!*********************************************************************
! Function Name: gauss
! Purpose: Computes Gaussian distribution
! Usage: call gauss(ix,sdv,mean) 
! Description: Computes a normally distributed random number with a 
!              given mean and standard deviation.
!              This routine uses the Box-Muller method for generating
!              random deviates with a normal (Gaussian) distribution
!              from a uniform distribution.  Uniform distribution generated
!              by RandomMT.f90 module functions.
! Inputs:
! Var_Name      Type      Description
! --------      ----      -----------
! iseed         integer   The seed. Initializes the random number
!                         generator on first call to gauss. If zero 
!                         then date/time function is used to set seed
!                         to a variable number. It is reset to zero
!                         inside gauss after the random number 
!                         generator has been initialized.
!  sdv          real      The desired standard deviation of the normal 
!                         distribution.
!  mean         real      The desired mean of the normal distribution.
!
! Outputs:
! Var_name    Type      Description
! --------    ----      -----------
! gauss       real      random number with Gaussian distribution
!
! Copyright: Atmospheric and Environmental Research, Inc., 1997, 2006
! Developed by Atmospheric and Environmental Research, Inc.
! Updated from f77 to f90 March, 2006
!*********************************************************************

  use RandomMT

  implicit none

!--- Calling Arguments      
  integer, intent(inout)    :: iseed
  real,    intent(in)       :: sdv
  real,    intent(in)       :: mean
  real                      :: gauss
 
!--- Local variables
  integer, save             :: ifirst=1 
  real                      :: dist,xx(2),zz
  integer                   :: DATE_TIME(8)
      
!--- On first call, set the random number generator seed
!--- use date/time function to assign seed if iseed=0
  if(ifirst .eq. 1) then
     if(iseed.ne.0) then 
        call initRandMT(iseed)
        iseed = 0
     else
        call date_and_time(values=DATE_TIME)
        iseed = DATE_TIME(8)*10000
        call date_and_time(values=DATE_TIME)
        iseed = iseed + DATE_TIME(8)
        call initRandMT(iseed)
        iseed = 0
     endif
  endif
  ifirst = 0

!--- generate the random numbers
  do
     xx(1) = randMT()
     xx(2) = randMT()
     xx = 2.0*xx - 1.0
     dist = xx(1)**2 + xx(2)**2
     if (dist < 1.0) exit
  end do
  zz = xx(1)*sqrt(-2.0*log(dist)/dist)
      
  gauss = zz*sdv + mean
  return

end function gauss

