subroutine setEmisQuasiPol(frq,ipol,nchan,nadirAng,ev,eh,em)
  use constants, ONLY : deg2rad
  implicit none
  real :: nadirAng,alpha,xxa,xxb
  integer :: nchan,n
  real,   dimension(nchan)  :: ev,eh,em,frq
  integer, dimension(nchan) :: ipol
  
  alpha=nadirAng*deg2rad
  xxa=cos(alpha)*cos(alpha)
  xxb=sin(alpha)*sin(alpha)
!!  where-elsewhere-end where can't be handled by 
!!  version 7.30 of fortran compiler on SGI  
!!$  where(ipol == 0)  ! quasi-vertical
!!$     em =xxa*ev+xxb*eh
!!$  elsewhere                   ! quasi-horizontal
!!$     em =xxa*eh+xxb*ev
!!$  end where
!!$ 
  do n = 1, nchan
     if (ipol(n) == 0) then
        em(n) =xxa*ev(n)+xxb*eh(n)
     else
        em(n) =xxa*eh(n)+xxb*ev(n)
     endif
  enddo
  return
end subroutine setEmisQuasiPol
