MODULE Quadrature 

!**********************************************************************
!* Module : Quadrature 
!*
!* Description: Collection of quadrature subroutines
!*
!*
!* Copyright: 2005, AER, Inc.
!* Developed by AER, Inc., 2005
!**********************************************************************

  implicit none

  private
  public  :: intTrapezoid

CONTAINS

!------------------------------------------------------

  subroutine intTrapezoid(x,y,intg)

! Description:  Integration of y = f(x) using trapezoid rule 

    real, dimension(:), intent(in)   :: x,y
    real              , intent(out)  :: intg

!-- Local variables

    real, allocatable, dimension(:)  :: xx,yy
    integer                          :: N

    N = size(x)
    allocate(xx(N-1),yy(N-1))

    yy = y(2:N)+y(1:N-1);
    xx = abs(x(2:N)-x(1:N-1));

    intg = sum(yy*xx/2.0,dim=1,mask=(xx /= 0.0))

    deallocate(xx,yy)
    return

  end subroutine intTrapezoid

end MODULE Quadrature 
