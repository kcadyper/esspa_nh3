MODULE RTFquadratic

  implicit none
  private
  public:: RTFquadraticFit, RTFquadraticVal

! Global private variables
  real :: c0    ! zero order coefficient
  real :: c1    ! 1st order coefficient
  real :: c2    ! 2nd order coefficient
  
CONTAINS

!************************************************************************
  subroutine RTFquadraticFit(rmax,r1,r2,v1,v2,c0r,c1r,c2r)
!  Compute the coefficients of a quadratic radiometer transfer function.
!  Nonlinearity is modeled in terms of maximum departure of count from
!  the linear result, deriving the maximum count departure from a given
!  maximum radiance departure.
!  The fit is made as a linear part with a quadratic departure.
! 
!  rmax  = maximum possible departure of radiance from linearity
!  r1    = radiance at first node of RTF
!  r2    = radiance at second node of RTF
!  v1    = count at first node of RTF
!  v2    = count at second node of RTF
!  

! Arguments
    real, intent(in)              :: rmax
    real, intent(in)              :: r1,r2
    real, intent(in)              :: v1,v2
    real, intent(inout), optional :: c0r,c1r,c2r

! Local variables
    real :: l0,l1,rdif,rsum,vmax
    real :: q0,q1,q2

! Compute the linear fit

    rdif=r2-r1
    l0=(v1*r2-v2*r1)/rdif
    l1=(v2-v1)/rdif

! Convert maximum non-linearity from radiance to count, making a 1st-order
! approximation

    vmax=rmax*l1

! Compute the quadratic departure

    rsum=r1+r2
    q0=vmax*(1-(rsum/rdif)**2)
    q1=vmax*4*rsum/(rdif**2)
    q2=-vmax*4/(rdif**2)

! Total coefficients

    c0=l0+q0
    c1=l1+q1
    c2=q2

    if (present(c0r)) c0r=c0
    if (present(c1r)) c1r=c1
    if (present(c2r)) c2r=c2

    return
  end subroutine RTFquadraticFit

!************************************************************************
  function RTFquadraticVal(rad,c0r,c1r,c2r)
!  Compute the value of a quadratic radiometer transfer function for a 
!  given radiance.
! !  
! Return value
  real :: RTFquadraticVal  ! Count equivalent to the input radiance

! Arguments
    real                          :: rad  ! Input radiance
    real, intent(inout), optional :: c0r,c1r,c2r

    if (present(c0r)) c0=c0r
    if (present(c1r)) c1=c1r
    if (present(c2r)) c2=c2r

    RTFquadraticVal=c0+c1*rad+c2*rad**2

    return
  end function RTFquadraticVal

END MODULE RTFquadratic
