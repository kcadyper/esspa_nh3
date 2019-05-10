!----------------------------------------------------------------------------
!
! MODULE RadFile
!
! This module contains routines for loading radfile dimensions, attributes
!   and variables, thereby compartmentalizing the data input.  After calling 
!   'getDimRad' and allocating necessary arrays, the user should call 
!   'getRadData' to load all data prior to looping over profiles.  
!
! This version is designed for the standard-format radfiles, for which
!   the I/O is handled using rad_io_module.f.  This version may be overridden
!   by a sensor-specific version.
!
!----------------------------------------------------------------------------
MODULE RadFile
  
  USE rad_io_module
  
  IMPLICIT none
  
  PRIVATE
  
  PUBLIC :: getDimRad, getRadData, getCloudData
  
  CONTAINS
    
    SUBROUTINE getDimRad(ncid,fName,nChan,nRec,nAng)
      
      ! Read in and return the rad file dimensions.
      
      !---I/O variables:
      integer         ,intent(inout) :: ncid
      character(len=*),intent( in  ) :: fname
      integer         ,intent(inout) :: nChan,nRec,nAng

      call queryRad(fName=fName,nChan=nChan,nRec=nRec,nAng=nAng)
      
    END SUBROUTINE getDimRad

    !--------------------------------------------------------------------------

    SUBROUTINE getRadData(ncid,fName,radAll,eiaAll,latAll,lonAll,wvn,chanID)
      
      ! Read in rad file variables and attributes and load into arrays.
      
      !---I/O variables:
      integer            ,intent(inout) :: ncid
      character (len=*)  ,intent(in)    :: fname
      real,dimension(:,:),intent(inout) :: radAll
      real,dimension(:  ),intent(inout) :: eiaAll,latAll,lonAll
      real,dimension(:  ),intent(inout) :: wvn
      character (len=*),dimension(:),intent(inout),optional :: chanID

      !---Local variables:
      integer                           :: ifor,nfor,nchan
      real                              :: lat,lon
      real,dimension(1)                 :: eia
      real,dimension(:),allocatable     :: rVec
      
      call openRad(ncid=ncid,fname=fName,wvn=wvn,nrec=nfor,nchan=nchan,status='old')
      
      if (allocated(rVec)) deallocate(rVec)
      allocate(rVec(nchan))
      
      do ifor=1,nfor
         call getRad(ncid=ncid,irec=ifor,lat=lat,lon=lon, &
              eia=eia,rad=rVec)
         radAll(:,ifor) = rVec
         latAll(ifor)   = lat
         lonAll(ifor)   = lon
         eiaAll(ifor)   = eia(1)
      enddo

      ! Channel IDs are currently not handled
      if (present(chanID)) chanID(:)=repeat(' ',len(chanID(1)))

      if (allocated(rVec)) deallocate(rVec)
      call closeRad(ncid)
      
    END SUBROUTINE getRadData

    !--------------------------------------------------------------------------

    SUBROUTINE getCloudData(ncid,fName,cldTypAll)
      
      ! Read in cloud-related variables from a rad file and load into arrays.
      ! This is currently a dummy; just to provide an interface that is
      ! consistent with other versions.
      
      !---I/O variables:
      integer                                    ,intent(inout) :: ncid
      character (len=*)                          ,intent(in)    :: fname
      integer(selected_int_kind(1)),dimension(:) ,intent(inout) :: cldTypAll

      return
      
    END SUBROUTINE getCloudData
    
  END MODULE RadFile
