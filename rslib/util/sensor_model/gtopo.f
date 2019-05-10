      subroutine gTopo(topoFile,Lat,Lon,Elev,Pland)
!**********************************************************************
!* Function Name: Gtopo
!* Purpose: Determine surface elevation and land-fraction/type 
!* Usage: call  Gtopo(lat,lon,elev,pland)
!* Description: For a given lat/lon, the surface elevation(elev) and 
!*              land-fractio/type(pland) are determined. The variable pland
!*              will vary between 0-1 or have value of 2. 
!*              pland    type
!*              -----    ----
!*              0        ocean
!*              1        land
!*              0->1     coast
!*              2        lake
!*              for ocean the elevation is 0. The elevation data are from
!*              USGS GTOPO30 and the surface type is from unknown source
!*
!* Inputs:
!* Var_name    Type   Description
!* --------    ----   -----------
!* Lat         real   latitude of point
!* Lon         real   longditude of point
!* 
!* Outputs:
!* Var_name    Type   Description
!* --------    ----   -----------
!* Elev        real   surface elevation at point
!* Pland       real   land-fractiontype at point
!* 
!* Copyright: AER, Inc., 1999-2004
!* Developed by AER, Inc., 1999-2004
!**********************************************************************
      implicit none
      integer, parameter :: lrec=7200
      integer*4, parameter :: Cols=7200,Rows=3600
      character(len=*), intent(in) :: topoFile
      real, intent(in) :: Lat,Lon
      real, intent(out) :: elev,pland
      real :: Slope=20.
      integer*2, dimension(lrec*2) :: Ige_Lat
      integer :: Ifirst, Lun,Ilat,Ilon
      data Ifirst/1/
      data Lun /21/
      save Lun
      
      if(Ifirst.eq.1) then
         Ifirst = 0
         open( Lun,file=topoFile,status='old',
     $        access='Direct',recl=lrec*4 )
      endif
      
      Ilat = Min( Max( Int(-Slope*( Lat- 90)+1), 1 ), Rows )
      Ilon = Min( Max( Int( Slope*(Lon+180)+1), 1 ), Cols )
      
      Read(Lun, Rec=ilat ) Ige_Lat
      
      Elev  = Ige_Lat(Ilon*2-1)
      Pland = Ige_Lat(Ilon*2)/float(1000)
      
      return
      end
      
