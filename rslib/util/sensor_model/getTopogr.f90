subroutine getTopogr(topoFile,dlat,dlon,lat,lon,hght,pland)
!*************************************************************************
!* Function name: getTopogr
!* Purpose: Get surface elevation 
!* Usage: call getTopogr(topoFile,dlat,dlon,lat,lon,hght,pland)
!* Description: Use digital elevation map to determine elevation and pland
!*              at lat, lon.
!* Inputs:
!* Var_Name      Type       Description
!* --------      ----       -----------
!* topoFile      character  topography file
!* dlat          real       latitude range
!* dlon          real       longitude range
!* lat           real       latitude in degrees
!* lon           real       longitue in degrees
!*
!* Outputs:
!* Var_Name      Type       Description
!* --------      ----       -----------
!* hght          real       elevation
!* pland         real       land type index
!*
!* Externals:
!* Name         Description
!* ----         -----------
!* Gtopo        reads elevation from direct access file
!*
!* Copyright: Atmospheric and Environmental Research, Inc., 1999-2004        
!* Developed by Atmospheric and Environmental Research, Inc.            
!*************************************************************************
  implicit none
  integer, parameter :: nrla=20,nrlo=20
  character(len=*), intent(inout) :: topoFile
  real, intent(in) :: dlat,dlon,lat,lon
  real, intent(out) :: hght,pland
  integer :: nlat,nlon,k,l
  real :: resolat,resolon,lat2,lon2,elevG,plandG
  
  resolat=1./float(nrla)
  resolon=1./float(nrlo)
  nlat=max(1,int(dlat/resolat)+1)
  nlon=max(1,int(dlon/resolon)+1)
  hght=0.
  pland=0.
  do k=1,nlon
     do l=1,nlat
        lat2=lat-(dlat/2.)+(l-1)*resolat
        lon2=lon-(dlon/2.)+(k-1)*resolon
        call Gtopo(topoFile,lat2,lon2,ElevG,PlandG)
        hght=hght+ElevG
        pland=pland+PlandG
     enddo
  enddo
  hght=hght/float(nlon*nlat)
  pland=pland/float(nlon*nlat)
  return
end subroutine getTopogr

