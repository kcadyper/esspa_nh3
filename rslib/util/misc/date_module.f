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

      MODULE date_module

! <f90Module>***********************************************************
!
! NAME:
!
!   date_module
!
! PURPOSE:
!
!   Tools for converting date information
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

      implicit none
      private
      public Julian, RelativeDay, Calendar, Julian2Calendar
      interface Julian
      module procedure JDAY
      end interface
      interface RelativeDay
      module procedure NDAY
      end interface
      interface Calendar
      module procedure GDATE
      end interface
      interface Julian2Calendar
      module procedure julday2calday
      end interface

      contains
******************************************************************************
***** Calendar subprograms                                              ******
******************************************************************************
      FUNCTION JDAY(IDOM,IM,IY)

!<f90Function>**********************************************************
!
! NAME:
!
!   JDAY
!
! PURPOSE:
!
!   Convert calendar date to Julian date.
!
! SYNTAX:
!
!   Results=JDAY(IDOM, IM, IY)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   IDOM  INTEGER  Day of the month
!   IM    INTEGER  Index of month
!   IY    INTEGER  Index of year
!
!   * OPTIONAL
!
! RETURN:
!
!     INTEGER  
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

***** Input is calendar date; day, month, year (year is four digits).
***** Returned value is Julian date.
      integer, intent(in) :: IDOM, IM, IY
      integer :: jday
      integer, dimension(12) :: NDS,ND
      integer :: i, j, ny
      DATA NDS/31,28,31,30,31,30,31,31,30,31,30,31/
***** Alter month calendar for leap years
      DO 130 J=1,12
         ND(J)=NDS(J)
 130  CONTINUE
      IF(MOD(IY,4).EQ.0)ND(2)=29
***** Check range of inputs
      IF (IM.LT.1 .OR. IM.GT.12) THEN
        PRINT *,'err[date_module::JDAY]: ',
     &       ' Month out of range; IDOM,IM,IY:',IDOM,IM,IY
        CALL errorHalt(1)
      ENDIF
      IF (IDOM.LT.1 .OR. IDOM.GT.ND(IM)) THEN
        PRINT *,'err[date_module::JDAY]: ',
     &       ' Day out of range; IDOM,IM,IY:',IDOM,IM,IY
        CALL errorHalt(1)
      ENDIF
      NY=0
***** Months
      DO 240 I=1,IM-1
         NY=NY+ND(I)
 240  CONTINUE
***** Days
      NY=NY+IDOM
      JDAY=NY
      RETURN
      END FUNCTION JDAY
      FUNCTION NDAY(IDOM,IM,IY)

!<f90Function>**********************************************************
!
! NAME:
!
!   NDAY
!
! PURPOSE:
!
!   Convert calendar date to relative day count.
!
! SYNTAX:
!
!   Results=NDAY(IDOM, IM, IY)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   IDOM  INTEGER  Day of the month
!   IM    INTEGER  Index of month
!   IY    INTEGER  Index of year
!
!   * OPTIONAL
!
! RETURN:
!
!     INTEGER  
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

***** Input is calendar date (year is four digits).
***** Returned value is relative day count. Day 1 is 1 Jan 1989.
      integer, intent(in) :: idom, im, iy
      integer :: nday
      integer, dimension(12) :: NDS,ND
      integer :: i, j, ny
      DATA NDS/31,28,31,30,31,30,31,31,30,31,30,31/
***** Alter month calendar for leap years
      DO 130 J=1,12
         ND(J)=NDS(J)
 130  CONTINUE
      IF (MOD(IY,4).EQ.0)ND(2)=29
      NY=0
***** Years
      DO 220 I=1989,IY-1
         NY=NY+365
         IF(MOD(I,4).EQ.0)NY=NY+1
 220  CONTINUE
***** Months
      DO 240 I=1,IM-1
         NY=NY+ND(I)
 240  CONTINUE
***** Days
      NY=NY+IDOM
      NDAY=NY
      RETURN
      END FUNCTION NDAY
      SUBROUTINE GDATE(IDAY,IDOM,IMM,IY,DAY)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   GDATE
!
! PURPOSE:
!
!   Convert relative day count to calendar date.
!
! SYNTAX:
!
!   CALL GDATE(IDAY, IDOM, IMM, IY, DAY)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   IDAY  INTEGER  Relative day count
!   
!   INPUTS/OUTPUTS:
!   
!   IDOM  INTEGER  Day of the month
!   IMM   INTEGER  Index of month
!   IY    INTEGER  Index of year
!   DAY   CHAR     Character form for day of the week
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

***** Input is relative day count.  Day 1 is Sun, 1 Jan 1989.
***** Output is date information
      integer, dimension(12) :: NDS,ND
      integer, intent(in)             :: iday
      integer, intent(inout)          :: idom, imm, iy
      character(len=*), intent(inout) :: DAY
      integer :: left, julian, npy, j, im, npm
      CHARACTER MON*3,DAYL(0:6)*3,MONL(12)*3
      DATA NDS/31,28,31,30,31,30,31,31,30,31,30,31/
      DATA DAYL/'Sun','Mon','Tue','Wed','Thu','Fri','Sat'/
      DATA MONL/'Jan','Feb','Mar','Apr','May','Jun','Jul',
     &     'Aug','Sep','Oct','Nov','Dec'/
***** Day of the week string
      DAY=DAYL(MOD(IDAY+6,7))
      LEFT=IDAY
***** Subtract off the whole years
      DO 100 IY=1989,1989+30
         NPY=365
         IF(MOD(IY,4).EQ.0)NPY=366
         LEFT=LEFT-NPY
         IF(LEFT.LE.0)GOTO 200
 100  CONTINUE
      PRINT *,'err[date_module::GDATE]: ',
     &     ' OVERRAN LOOP 100 IN SUBROUTINE DATE'
      CALL errorHalt(1)
 200  JULIAN=LEFT+NPY
***** Alter month calendar for leap years
      DO 230 J=1,12
         ND(J)=NDS(J)
 230  CONTINUE
      IF(MOD(IY,4).EQ.0)ND(2)=29
***** Identify the month
      IDOM=JULIAN
      DO 300 IM=1,12
         NPM=ND(IM)
         IDOM=IDOM-NPM
         IF(IDOM.LE.0)GOTO 310
 300  CONTINUE
      PRINT *,'err[date_module::GDATE]: ',
     &     ' OVERRAN LOOP 300 IN SUBROUTINE DATE'
      CALL errorHalt(1)
 310  IDOM=IDOM+NPM
c      MON=MONL(IM)
      IMM=IM
      RETURN
      END SUBROUTINE GDATE

      subroutine julday2calday(yr,jday,mth,idom)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   julday2calday
!
! PURPOSE:
!
!   Converts Julian day and year to calendar date.
!
! SYNTAX:
!
!   CALL julday2calday(yr, jday, mth, idom)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   yr    INTEGER  year
!   jday  INTEGER  day of a year
!   
!   INPUTS/OUTPUTS:
!   
!   mth   INTEGER  month
!   idom  INTEGER  Day of the month
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

c**************************************************************
c* Converts Julian day to calendar by given year and Julian day
c* 
c* History:
c*     11/18/2003 - Imported from write_amsusdr by Y. He.
c**************************************************************
      integer , intent(in)    :: yr, jday
      integer , intent(inout) :: mth,idom
      integer , dimension(12) ::
     $     mdays=(/31,28,31,30,31,30,31,31,30,31,30,31/),
     $     mdaysLeap=(/31,29,31,30,31,30,31,31,30,31,30,31/)
      integer :: im, sumday, i

      sumday=0
      im=0

      if (mod(yr,4).ne.0) then  ! not leap year
         do i=1,12 
            im=im+1
            sumday=sumday+mdays(i)
            if(sumday.ge.jday) exit
         end do
         sumday=sumday - mdays(im)
      else                      ! leap year
         do i=1,12 
            im=im+1
            sumday=sumday+mdaysLeap(i)
            if(sumday.ge.jday) exit
         end do
         sumday=sumday - mdaysLeap(im)
      end if

      mth = im
      idom = jday - sumday

      return
      end subroutine julday2calday

      end module date_module
