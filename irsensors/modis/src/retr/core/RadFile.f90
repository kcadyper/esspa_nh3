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

!----------------------------------------------------------------------------
!
! MODULE RadFile
!
! This module contains routines for loading radfile dimensions, attributes
!   and variables, thereby compartmentalizing the data input.  After calling 
!   'getDimRad' and allocating necessary arrays, the user should call 
!   'getRadData' to load all data prior to looping over profiles.  
!
!----------------------------------------------------------------------------
MODULE RadFile

! <f90Module>***********************************************************
!
! NAME:
!
!   RadFile
!
! PURPOSE:
!
!   Contains routines for loading radfile dimensions, attributes and 
!   variables, thereby compartmentalizing the data input.
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

  
  USE ncdf_module
  
  IMPLICIT none
  
  PRIVATE
  
  PUBLIC :: getDimRad,getRadData,getCloudData, &
            getGeogData
  
  CONTAINS
    
    SUBROUTINE getDimRad(ncid,fName,nChan,nRec,nAng)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   getDimRad
!
! PURPOSE:
!
!   Read in and return the rad file dimensions.
!
! SYNTAX:
!
!   CALL getDimRad(ncid, fName, nChan, nRec, nAng)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   fName  CHAR     File name
!   
!   INPUTS/OUTPUTS:
!   
!   ncid   INTEGER  ID for ncdf file 
!   nChan  INTEGER  Number of channels
!   nRec   INTEGER  Number of records
!   nAng*  INTEGER  number of angles
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

      
      ! Read in and return the rad file dimensions.
      
      !---I/O variables:
      integer         ,intent(inout)           :: ncid
      character(len=*),intent(in)              :: fname
      integer         ,intent(inout)           :: nChan,nRec
      integer         ,intent(inout), optional :: nAng

      !---Local variables:
      integer,parameter              :: nAngPar=1
      
      if (present(nAng)) nAng = nAngPar

      call openNcdfFile(ncid, fname, status='old')
      nChan = readNcdfDim(ncid,'nChannels')
      nRec  = readNcdfDim(ncid,'nRec')
      call closeNcdfFile(ncid)
      
    END SUBROUTINE getDimRad

    !--------------------------------------------------------------------------

    SUBROUTINE getRadData(ncid,fName,radAll,eiaAll,latAll,lonAll,wvn,time,chanID)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   getRadData
!
! PURPOSE:
!
!   Read in rad file variables and attributes and load into arrays.
!
! SYNTAX:
!
!   CALL getRadData(ncid, fName, radAll, eiaAll, latAll, lonAll, wvn, 
!      time, chanID)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   fName   CHAR     File name
!   
!   INPUTS/OUTPUTS:
!   
!   ncid    INTEGER  ID for ncdf file 
!   radAll* REAL     Radiances for all locations and channels
!   eiaAll* REAL     Satellite Earth incidence (zenith) angles
!   latAll* REAL     Latitudes
!   lonAll* REAL     Longitudes
!   wvn*    REAL     Wavenumbers
!   time*   INTEGER  Date/time
!   chanID* CHAR     Channel identifier string
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

      
      ! Read in rad file variables and attributes and load into arrays.
      
      !---I/O variables:
      integer             ,intent(inout)          :: ncid
      character (len=*)   ,intent(in)             :: fname
      real,dimension(:,:) ,intent(inout),optional :: radAll
      real,dimension(:  ) ,intent(inout),optional :: eiaAll,latAll,lonAll
      real,dimension(:  ) ,intent(inout),optional :: wvn
      integer,dimension(6),intent(inout),optional :: time
      character (len=*),dimension(:),intent(inout),optional :: chanID
      
      !---Local variables:
      integer(selected_int_kind(8)),dimension(:),allocatable   :: iVec
      integer, parameter :: maxIDlen=4000
      character (len=maxIDlen) :: chanIDall
      character(len=20) :: radVarName='rad',eiaVarName='SATZEN'
      character(len=20) :: latVarName='lat',lonVarName='lon'
      character(len=19) :: sotime
      real              :: scale,offset
      integer           :: ncStatus,tmpVarID
      integer           :: nPrf
      integer           :: lenID
      integer           :: nChanLoc
      integer           :: lenTot
      integer           :: ich,iid1,iid2
      
      call openNcdfFile(ncid, fname, status='old')

      nPrf     = readNcdfDim(ncid,'nRec')
      
      if (allocated(iVec)) deallocate(iVec)
      allocate(iVec(nPrf))

      !---Read in radiance data:
      if (present(radAll)) then
         ncStatus = nf_inq_varid(ncid,trim(radVarName),tmpvarid)
         ncStatus = nf_get_var_real(ncid,tmpvarid,radAll)
         if (ncStatus /= nf_noErr) call handle_nc_err(ncStatus,'getRadData', &
              'Rad file does not contain variable '//trim(radVarName))
         ncStatus = nf_get_att_real(ncid, tmpvarid, 'scale' ,scale)
         ncStatus = nf_get_att_real(ncid, tmpvarid, 'offset',offset)
         radAll   = radAll*scale + offset
      endif

      !---Read in EIA data:
      if (present(eiaAll)) then
         ncStatus = nf_inq_varid(ncid,trim(eiaVarName),tmpvarid)
         ncStatus = nf_get_var_int(ncid,tmpvarid,iVec)
         if (ncStatus /= nf_noErr) call handle_nc_err(ncStatus,'getRadData', &
              'Rad file does not contain variable '//trim(eiaVarName))
         ncStatus = nf_get_att_real(ncid, tmpvarid, 'scale' ,scale)
         ncStatus = nf_get_att_real(ncid, tmpvarid, 'offset',offset)
         eiaAll   = iVec*scale + offset
      endif

      !---Read in latitude data:
      if (present(latAll)) then
         ncStatus = nf_inq_varid(ncid,trim(latVarName),tmpvarid)
         ncStatus = nf_get_var_int(ncid,tmpvarid,iVec)
         if (ncStatus /= nf_noErr) call handle_nc_err(ncStatus,'getRadData', &
              'Rad file does not contain variable '//trim(latVarName))
         ncStatus = nf_get_att_real(ncid, tmpvarid, 'scale' ,scale)
         ncStatus = nf_get_att_real(ncid, tmpvarid, 'offset',offset)
         latAll   = iVec*scale + offset
      endif

      !---Read in longitude data:
      if (present(lonAll)) then
         ncStatus = nf_inq_varid(ncid,trim(lonVarName),tmpvarid)
         ncStatus = nf_get_var_int(ncid,tmpvarid,iVec)
         if (ncStatus /= nf_noErr) call handle_nc_err(ncStatus,'getRadData', &
              'Rad file does not contain variable '//trim(lonVarName))
         ncStatus = nf_get_att_real(ncid, tmpvarid, 'scale' ,scale)
         ncStatus = nf_get_att_real(ncid, tmpvarid, 'offset',offset)
         lonAll   = iVec*scale + offset
      endif

      !---Read in wavenumber grid:
      if (present(wvn)) call readNcdfAttr(ncid,attrName='wavenumbers',attr=wvn)

      !---Read in time string and parse into integer elements:
      if (present(time)) then
         call readNcdfAttr(ncid,attrName='SensorObservationTime',attr=sotime)
         read(sotime,'(i4,5(1x,i2))')time
         time(6)=time(6)*1000 ! seconds to milliseconds
      endif

      !---Read channel ID strings
      if (present(chanID)) then
         call readNcdfAttr(ncid,attrName='lenID',attr=lenID,silent=.TRUE.)
         if (lenID == 0) then   ! No such attribute in file
            chanID(:)=repeat(' ',len(chanID(1)))
         else
            nChanLoc = readNcdfDim(ncid,'nChannels')
            if (len(chanID(1)) /= lenID) then
               print *, "err[RadFile::getRadData]: ", &
               "Length of chanID does not agree with file content."
               print *,'len(chanID(1)),lenID:',len(chanID(1)),lenID
               call errorHalt(1)
            endif
            lenTot=(lenID+1)*nChanLoc-1
            if (lenTot > maxIDlen) then
               print *, "err[RadFile::getRadData]: ", &
               "Total length of chanID exceeds length of chanIDall"
               print *,'lenTot,maxIDlen:',lenTot,maxIDlen
               call errorHalt(1)
            endif
            chanIDall=repeat(' ',maxIDlen)
            call readNcdfAttr(ncid,attrName='chanID',attr=chanIDall(1:lenTot))
            do ich = 1, nChanLoc
               iid1=(ich-1)*(lenID+1)+1
               iid2=ich*(lenID+1)-1
               chanID(ich)=adjustl(chanIDall(iid1:iid2))
            enddo
         endif
      endif
      
      call closeNcdfFile(ncid)
      
      if (allocated(iVec)) deallocate(iVec)
      
    END SUBROUTINE getRadData

    !--------------------------------------------------------------------------
    
    SUBROUTINE getCloudData(ncid,fname,cldTypAll)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   getCloudData
!
! PURPOSE:
!
!   Load cloud-related variables.
!
! SYNTAX:
!
!   CALL getCloudData(ncid, fname, cldTypAll)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   fname      CHAR     File name
!   
!   INPUTS/OUTPUTS:
!   
!   ncid       INTEGER  ID for ncdf file 
!   cldTypAll  INTEGER  IDs of cloud type class
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>


      ! Load cloud-related variables.

      !---I/O variables:
      integer,                                   intent(inout) :: ncid
      character (len=*),                         intent(in)    :: fname
      integer(selected_int_kind(1)),dimension(:),intent(inout) :: cldTypAll

      !---Local variables:
      character(len=20) :: typVarName='PHASE'
      integer           :: ncStatus,tmpVarID

      call openNcdfFile(ncid, fname, status='old')

      !---Read in cloud type data:
      ncStatus = nf_inq_varid(ncid,trim(typVarName),tmpvarid)
      ncStatus = nf_get_var_int1(ncid,tmpvarid,cldTypAll)
      if (ncStatus /= nf_noErr) call handle_nc_err(ncStatus,'getRadData', &
           'Rad file does not contain variable '//trim(typVarName))

      call closeNcdfFile(ncid)

    END SUBROUTINE getCloudData
    
    !--------------------------------------------------------------------------
    
    SUBROUTINE getGeogData(ncid,fname,geography)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   getGeogData
!
! PURPOSE:
!
!   Read geographysical data stored as short
!
! SYNTAX:
!
!   CALL getGeogData(ncid, fname, geography)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   fname      CHAR     File name
!   
!   INPUTS/OUTPUTS:
!   
!   ncid       INTEGER  ID for ncdf file 
!   geography  INTEGER  Geography IDs
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    ! Read geography data stored as short
    ! The I/O parameters
      INTEGER,                                     INTENT(INOUT) :: ncid
      CHARACTER(LEN=*),                            INTENT(IN)    :: fname
      INTEGER(SELECTED_INT_KIND(4)), DIMENSION(:), INTENT(INOUT) :: geography

    ! Local variables
      INTEGER            :: ncStatus,varID
      CHARACTER(LEN=15)  :: varName = 'GEOGRAPHY'

      CALL openNcdfFile(ncid,fname,status='old')
      ncStatus = nf_inq_varid(ncid,varName,varID)
      IF (ncStatus .NE. nf_noErr) CALL handle_nc_err(ncStatus,'getGeogData', &
         "Rad file does not contain variable " // varName)
      ncStatus = nf_get_var_int(ncid,varID,geography)

      IF (ncStatus .NE. nf_noErr) CALL handle_nc_err(ncStatus,'getGeogData', &
          'Failed to read data ' // varName)

      CALL closeNcdfFile(ncid)
    END SUBROUTINE getGeogData

    ! --------------------------------------------------------------------------

    SUBROUTINE handle_nc_err(ncStatus,subproName,msg)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   handle_nc_err
!
! PURPOSE:
!
!   Report error condition.
!
! SYNTAX:
!
!   CALL handle_nc_err(ncStatus, subproName, msg)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   ncStatus    INTEGER  Error status
!   subproName  CHAR     Name of subroutine that called
!   msg*        CHAR     Text of message
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>


      integer,           intent(in)           :: ncStatus
      character (len=*), intent(in)           :: subproName
      character (len=*), intent(in), optional :: msg
      character (len=256) :: lmsg

      if (present(msg)) then
         if (len(msg) > 256) print*,'message is truncated'
         lmsg=trim(msg)
      else
         lmsg=':::'
      endif
      print *, 'err[RadFile::handle_nc_err(',trim(subproName),')]: ',&
           ' ('//trim(msg)//') '//trim(nf_strerror(ncStatus))
      call errorHalt(1)

    END SUBROUTINE handle_nc_err

END MODULE RadFile
