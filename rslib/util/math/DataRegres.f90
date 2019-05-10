MODULE DataRegres 

! <f90Module>***********************************************************
!
! NAME:
!
!   DATAREGRES
!
! PURPOSE:
!
!   Collection of 1-D regression subroutines 
!
! INCLUDES:
!
!   None
!
! COMMENTS: Extrapolation is allowed 
!
!   None
!
! MODIFICATION_HISTORY:
!
!   Written by: Haijun Hu, Fri Jan 13, 2006
!
!***********************************************************</f90Module>


  implicit none

  private
  public:: curveFit
  public:: linearFit !quadraticFit

  interface curveFit
     module procedure linearFit
  end interface

CONTAINS

!************************************************************************
  subroutine linearFit(xi,yi,xins,youts,flg1d) 

!<f90Subroutine>********************************************************
!
! NAME:
!
!   LINEARFIT
!
! PURPOSE:
!
!  Fit pairs of input data (xi, yi) into a straight line
!       y= a0 + a1* x
!  using least-squares method.
!
! SYNTAX:
!
!   CALL LINEARFIT(XI, YI, XINS, YOUTS, FLG1D)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   XI     REAL     x values for curve fitting
!   YI     REAL     y values for curve fitting
!   XINS   REAL     x values for interpolation
!   FLG1D  INTEGER  vector flags for data filtering
!   
!   OUTPUTS:
!   
!   YOUTS  REAL     intepolated y values 
!
! INCLUDES:
!
!   None
!
! COMMENTS:
!
!   None
!
! MODIFICATION_HISTORY:
!
!   Written by: Haijun Hu, Fri Jan 13, 2006
!
!*******************************************************</f90Subroutine>

!-- Calling Arguments
    real,    dimension(:), intent(in)        :: xi
    real,    dimension(:), intent(in)        :: yi 
    real,    dimension(:), intent(in)        :: xins 
    real,    dimension(:), intent(out)       :: youts
    integer, dimension(:), intent(in)        :: flg1d


! Local variables
    integer                         :: i,N,M
    real                            :: a0,a1
    real                            :: xsum,ysum,anorm,dnorm
    real                            :: Nval

    N = size(xi)
    M = size(xins)
    if ( (size(yi) /= N) .or. (size(flg1d) /= N) &
         .or. (size(youts) /= M) ) then
       print*, 'err[DataRegres]: X, Y sizes mismatch'
       print*, 'size(xi),size(yi),size(flg1d)=',size(xi),size(yi),size(flg1d)
       print*, 'size(xins), size(youts)=', size(xins), size(youts)
       call errorHalt(1)
    endif


    !if ( (xins(M) > xi(N)) .or. (xins(1) < xi(1)) ) then
    !   print*, 'warning[DataRegres]: Xins out of range'
    !   print*, 'N,xi(1)...xi(N)=', N,xi(1), xi(N) 
    !   print*, 'M,xins(1)...xins(M)=',M,xins(1),xins(M)
    !endif

    select case (N)
    case (1)
       a1 = 0.0
       a0 = yi(1)
       if (flg1d(1) /= 0) then 
          print*, 'err[DataRegres]: Only Input Bad'
          call errorHalt(1)
       endif
    case (2)
       a1 = (yi(2)-yi(1))/(xi(2)-xi(1))
       a0 = yi(1)
       if ((flg1d(1) /= 0) .or. (flg1d(2) /= 0)) then
          print*, 'err[DataRegres]: 1 of 2 Inputs Bad'
          call errorHalt(1)
       endif
    case (3:)
       xsum  = sum(xi, dim=1, mask=(flg1d==0) )
       ysum  = sum(yi, dim=1, mask=(flg1d==0) )
       Nval  = float(N - sum(flg1d, dim=1, mask=(flg1d /=0 )))
       anorm = sum(xi*yi, dim=1, mask=(flg1d==0)) - xsum*ysum/Nval
       dnorm = sum(xi*xi, dim=1, mask=(flg1d==0)) - xsum*xsum/Nval
        if (dnorm == 0.0) then
          print*, 'err[DataRegres]: divided by zero'
          call errorHalt(1)
        endif
       a1 = anorm/dnorm
       a0 = (ysum - a1*xsum)/Nval
    end select

    youts = a0 + a1*xins
     
    return

  end subroutine linearFit

!-----------------------------------------------------------------

END MODULE DataRegres 
