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

MODULE CloudTab

! <f90Module>***********************************************************
!
! NAME:
!
!   CloudTab
!
! PURPOSE:
!
!   The module produces cloud properties from tables of gridded data.
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>


  ! This module produces cloud properties from tables of gridded data.
  ! The user provides a netCDF file that contains the data. The netCDF
  ! file must contain these anciliary data: - the nodes of the grid
  ! for each dimension of the data - strings that define the type of
  ! interpolation (e.g. linear) to use for each dimension (it is
  ! assumed that the table was built to achieve required accuracy with
  ! a certain type of interpolation) - a flag for each grid that
  ! indicates whether the grid is uniformly spaced, which speeds
  ! computation. The data may represent scattering and absorption
  ! coefficients, asymmetry factor, etc. There may be more than one
  ! type of data in the file. Each type of data has a distinct netCDF
  ! variable name. The data may be on a grid with a number of
  ! dimensions that the current code is designed to handle, where the
  ! dimensions may cover: primary size parameter (may be effective
  ! diameter, mass-diameter, etc.), secondary size parameter (may be
  ! variance or width of size ditribution), temperature (K), density,
  ! melt fraction, etc., and spectral point (may be wavenumber,
  ! frequency, etc.). With the exception of spectral point, the
  ! dimension may be 1, whereby the property is assumed not to vary in
  ! that dimension. The subroutine getCloudTab provides cloud property
  ! data for all spectral points interpolated to the specified values
  ! of physical properties such as size parameters mentioned above.

  USE NETCDF
  USE CloudParameters
  USE CloudDataStruct
  USE CloudNameDictionary, Only: matchName

  implicit none
  private
  public:: initCloudTab, getCloudTab, destroyCloudTab

  ! Global private data
  character(len=mxNameLength), dimension(mxIntplType)  :: validInterpTypes=(/'none  ','linear'/)
  CHARACTER(mxNameLength),DIMENSION(mxOpticalProperty) :: loadedOptProp
  logical, dimension(mxHydrometeor)  :: initialized=.false.

  ! loadDone implementation provides the option to call loadCloudTab() in advance
  logical                               :: loadDone=.false.
  logical                               :: first_init=.true.
  logical                               :: outBndsExit=.false.
  ! To allow cloud data tabulation to use different grid for each
  ! data type (see EmisTab):
  !  * move myspecPtsGrd into structure Variables
  !  * add function to interpolate from spectral point grid to input
  !    spectral position, such that parent program does not need to
  !    know spectral point grid, the grid is
  !  * remove nSpcPts and spectPtsGrd from cloudTabInit calling arguments
  !-
  !- Comments:
  !-
  !-  May need to come back to address or revise the requests above
  !-
  real, dimension(:), allocatable, target  :: myspectPtsGrd
  integer                                  :: nSPsave=0
  integer                                  :: mynSpcPts
  integer :: npProperty  !obtained from the cloud table file
  integer :: noProperty  !obtained from the cloud table file
  integer :: nHydr
  integer, dimension(mxHydrometeor,mxPhysicalProperty) :: indMap

  TYPE(CloudTable), DIMENSION(:), allocatable :: cldTbl
  TYPE(CloudOpticalSet), DIMENSION(:), allocatable :: cldOptSet

CONTAINS

  SUBROUTINE initCloudTab(fname,dataName_in,cldScene,outBndsExit_in,nSpcPts,spectPtsGrd)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   cloudTabInit
!
! PURPOSE:
!
!   Initialize access to tabulated cloud data.
!
! SYNTAX:
!
!   CALL cloudTabInit(fname, dataName_in, outBndsExit_in, nSpcPts,
!      spectPtsGrd)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   fname           CHAR     File name
!   dataName_in     CHAR     String combining hydrometeor and optical
!                            property ID
!   outBndsExit_in* LOGICAL  Flag to exit if value is out of bounds
!
!   INPUTS/OUTPUTS:
!
!   nSpcPts*        INTEGER  Number of spectral grid points
!   spectPtsGrd*    REAL     Spectral grid
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>


    !  Initialize access to tabulated cloud data. Data table is assumed
    !  to represent a spectral point set and the user may verify that
    !  the set is as expected by returning spectPtsGrd and nSpcPts.
    use CloudNameDictionary, only : matchName

    !---Input/Output variables:
    character(len=*),               intent(in)              :: fname
    character(len=*), dimension(:), intent(in)              :: dataName_in
    TYPE(CloudScene), DIMENSION(:), INTENT(IN)              :: cldScene
    logical,                        intent(in),    optional :: outBndsExit_in
    integer,                        intent(inout), optional :: nSpcPts
    real,             dimension(:), pointer,       optional :: spectPtsGrd
    integer                                                 :: iop,iopp,iHydr,jHydr,ipp,ispp
    integer                                                 :: currentHydr
    logical                                                 :: hydrFound
    logical                                                 :: optFound
    character (len=*), parameter                            :: procName=" [CloudTab::initCloudTab]: "
    !-- Read file

    !  The job of loadCloudTab() can optionally be called from the top-level
    !  program and executed at initialization by first making calls to:
    !    Retr::getNumHydr (to return nHydr)
    !    allocate memory for dataType vector
    !    Retr::getPropertiesHydrType (to return dataType)
    !    CloudModule::buildCldDataNames
    !    allocate memory for dataName vector
    !    CloudModule::getCldDataNames
    !  to construct the dataName_in input argument for loadCloudTab()

    call loadCloudTab(fname,dataName_in)

    if (noProperty > mxOpticalProperty) then
       print *, procName, 'err - loaded number of optical properties greater than allowed'
       call exit(55)
    end if

    do iop = 1, noProperty
       optFound = .false.
       do iopp = 1, mxOpticalProperty
          if (trim(knownOptProp(iopp)) /= trim(loadedOptProp(iop))) cycle
          optFound = .true.
       end do
       if (.not. optFound) then
          print *, procName, 'err - unknown optical table - '//trim(loadedOptProp(iop))
          call exit(56)
       end if
    end do

    if (present(nSpcPts)) nSpcPts=mynSpcPts

    if (present(spectPtsGrd)) then
       if (associated(spectPtsGrd)) deallocate(spectPtsGrd)
       allocate (spectPtsGrd(mynSpcPts))
       spectPtsGrd=mySpectPtsGrd
    endif

    if(present(outBndsExit_in))outBndsExit=outBndsExit_in

    !
    !Check if the caller-requested cloud physical properties are
    !provided and in the right order. Otherwise, stop. (YHE: may
    !consider more flexibility of using the optical table in the
    !future, if the requested physical properties are out of order
    !from the ones provided by the optical table.
    !
    if (size(cldScene) > nHydr) then
       print *, procName, 'err - requested number of hydrometeors exceeded'
       call exit(71)
    end if

    do iHydr = 1, size(cldScene)
       hydrFound = .false.
       do jHydr = 1, size(cldTbl)
          if (trim(cldScene(iHydr)%name) == trim(cldTbl(jHydr)%name)) then
             currentHydr = jHydr
             hydrFound = .true.
             exit
          end if
       end do
       if (.not. hydrFound) then
          print *, procName, 'err - unrecognized hydrometeor name'//&
               trim(cldScene(iHydr)%name)//','//trim(cldTbl(jHydr)%name)
          call exit(72)
       end if
       do ipp = 1, cldTbl(currentHydr)%npProperty
          do ispp = 1, cldScene(iHydr)%nProperty
             if ( matchName(cldScene(iHydr)%physDescr(ispp)%name,&
                            cldTbl(currentHydr)%pProperty(ipp)%name)) then
                indMap(currentHydr,ipp) = ispp
                exit
             end if
          end do
       end do

       if (cldScene(iHydr)%nProperty > cldTbl(currentHydr)%npProperty) then
          print *, procName, 'err - number of physical properties exceeded'
          call exit(73)
       end if
    end do

    return

  END SUBROUTINE initCloudTab

  SUBROUTINE loadCloudTab(fname,dataName_in)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   loadCloudTab
!
! PURPOSE:
!
!   Load tabulated cloud data.
!
! SYNTAX:
!
!   CALL loadCloudTab(fname, dataName_in)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   fname        CHAR  File name
!   dataName_in  CHAR  String combining hydrometeor and optical property ID
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !  Initialize access to tabulated cloud data. Data table is assumed
    !  to represent a spectral point set.

    !---Input/Output variables:
    character(len=*),               intent(in) :: fname
    character(len=*), dimension(:), intent(in) :: dataName_in

    !-- Local variables
    integer, dimension(:), allocatable    :: mydimid
    integer                               :: tmpvarid
    integer                               :: ncid,ncStatus,istatus
    character(len=mxNameLength)           :: htName

    INTEGER, DIMENSION(:), ALLOCATABLE :: pdimArr
    CHARACTER (LEN=mxNameLength) :: pGridPreFix = 'grid'
    CHARACTER (LEN=mxNameLength) :: intplPrefix = 'interp'
    CHARACTER (LEN=mxNameLength) :: uGridPreFix = 'unifGrid'
    CHARACTER (LEN=mxNameLength) :: unitPreFix = 'unit'
    CHARACTER (len=mxNameLength), DIMENSION(:), ALLOCATABLE :: tblNames
    CHARACTER (len=mxNameLength), DIMENSION(:), ALLOCATABLE :: gAttrNames
    CHARACTER (len=mxNameLength), DIMENSION(mxHydrometeor)  :: hydrNames
    CHARACTER (len=mxNameLength) :: dimName
    CHARACTER (len=mxNameLength) :: spcName
    CHARACTER (len=mxNameLength) :: optName

    INTEGER :: iHydr, ih, itbl, iattr, ipp, iop, iopp, ix, iy, iz
    INTEGER :: nDims, nGlobalDims, nTbls, nAttrs
    INTEGER :: unlimDimId
    INTEGER :: dimLen
    LOGICAL :: attrFound

    real*4, allocatable :: buffer1D(:)
    real*4, allocatable :: buffer2D(:,:)
    real*4, allocatable :: buffer3D(:,:,:)
    real*4, allocatable :: buffer4D(:,:,:,:)
    real*4, allocatable :: buffer5D(:,:,:,:,:)
    real*4, allocatable :: buffer6D(:,:,:,:,:,:)

    if(loadDone) return

    !Start loading data from the cloud table
    ncStatus = nf90_open(fname, NF90_NOWRITE, ncid)
    if (ncStatus /= nf90_NoErr) call handle_err(ncStatus)
    ncStatus = nf90_Inquire(ncid, nDimensions=nGlobalDims, nVariables=nTbls, &
         nAttributes=nAttrs, unlimitedDimId=unlimDimId)
    if(ncStatus /= nf90_NoErr) call handle_err(ncStatus)

    if (.not. allocated(gAttrNames)) allocate(gAttrNames(nAttrs))
    do iattr = 1, nAttrs
       ncStatus = nf90_inq_attname(ncid, NF90_GLOBAL, iattr, gAttrNames(iattr))
       if (ncStatus /= nf90_NoErr) call handle_err(ncStatus)
    end do

    if (.not. allocated(tblNames)) allocate(tblNames(nTbls))
    do itbl = 1, nTbls
       ncStatus = nf90_inquire_variable(ncid, itbl, name=tblNames(itbl))
       if(ncStatus /= nf90_NoErr) call handle_err(ncStatus)
    end do

    ! Collect the hydrometeor names in the file
    iHydr = 0
    do itbl = 1, nTbls
       ix = index(trim(tblNames(itbl)), '_')
       htName = tblNames(itbl)(1:ix-1)
       if (iHydr == 0) then
          iHydr = iHydr + 1
          hydrNames(iHydr) = trim(htName)
       else
          do ih = 1, iHydr
             ix = index(trim(tblNames(itbl)),trim(hydrNames(ih)))
             if (ix /= 0) exit
          end do
          if (ix == 0) then
             iHydr = iHydr + 1
             hydrNames(iHydr) = trim(htName)
          end if
       end if
    end do
    nHydr = iHydr

    ! Collect the optical table names in the file
    iop = 0
    do itbl = 1, nTbls
       ix = index(trim(tblNames(itbl)), '_')
       optName = tblNames(itbl)(ix+1:)
       if (iop == 0) then
          iop = iop + 1
          loadedOptProp(iop) = trim(optName)
       else
          do iopp = 1, iop
             ix = index(trim(tblNames(itbl)),trim(loadedOptProp(iopp)))
             if (ix /= 0) exit
          end do
          if (ix == 0) then
             iop = iop + 1
             loadedOptProp(iop) = trim(optName)
          end if
       end if
    end do

    noProperty = iop
    if (.not. allocated(cldTbl)) allocate(cldTbl(nHydr))
    if (.not. allocated(cldOptSet)) allocate(cldOptSet(nHydr))

    iHydr = 0
    do itbl = 1, nTbls
       !Any one of the optical names can be used to identify the
       !hydrometeor names
       ix =index(tblNames(itbl),trim(loadedOptProp(1)))
       if (ix  /= 0) then
          iHydr = iHydr + 1
       end if
       ncStatus = nf90_inquire_variable(ncid, itbl, ndims=nDims)
       npProperty = nDims - 1  !exclude the dimension of spectral points
       cldTbl(iHydr)%name = trim(hydrNames(iHydr))
       cldTbl(iHydr)%npProperty = npProperty
       cldOptSet(iHydr)%name = trim(hydrNames(iHydr))
       cldOptSet(iHydr)%nProperty = noProperty
       if (.not. allocated(mydimid)) allocate(mydimid(nDims))
       if (.not. allocated(cldTbl(iHydr)%pProperty)) &
            allocate(cldTbl(iHydr)%pProperty(npProperty))
       if (.not. allocated(cldTbl(iHydr)%oProperty)) &
            allocate(cldTbl(iHydr)%oProperty(noProperty))
       if (.not. allocated(cldOptSet(iHydr)%property)) &
            allocate(cldOptSet(iHydr)%property(noProperty))
    end do

    do itbl = 1, nTbls
       do iHydr = 1, nHydr
          htName = trim(hydrNames(iHydr))
          ix =index(tblNames(itbl),trim(hydrNames(iHydr)))
          if (ix  /= 0) exit
       end do

       ncStatus = nf90_inquire_variable(ncid, itbl, ndims=nDims)
       if(ncStatus /= nf90_NoErr) call handle_err(ncStatus)

       ncStatus = nf90_inq_varid(ncid,trim(tblNames(itbl)),tmpvarid)
       if(ncStatus /= nf90_NoErr) call handle_err(ncStatus)

       ncStatus = nf90_inquire_variable(ncid,tmpvarid,dimids=mydimid)
       if(ncStatus /= nf90_NoErr) call handle_err(ncStatus)

       npProperty = cldTbl(iHydr)%npProperty
       if (.not. allocated(pdimArr)) allocate(pdimArr(npProperty))
       pdimArr = 0
       do ipp = 1, npProperty
          ncStatus = nf90_inquire_dimension(ncid,mydimid(ipp),name=dimName,len=dimLen)
          if(ncStatus /= nf90_NoErr) call handle_err(ncStatus)
          if (.not. allocated(cldTbl(iHydr)%pProperty(ipp)%grid)) &
               allocate(cldTbl(iHydr)%pProperty(ipp)%grid(dimLen))
          ix =index(dimName,trim(htName))

          if (ix  /= 0) then
             cldTbl(iHydr)%pProperty(ipp)%name = trim(dimName(ix+len(trim(htName)):))
          end if
          cldTbl(iHydr)%pProperty(ipp)%size = dimLen
          pdimArr(ipp) = dimLen

          do iattr = 1, nAttrs
             attrFound = matchName(gAttrNames(iattr),trim(cldTbl(iHydr)%pProperty(ipp)%name), &
                  trim(intplPrefix)//trim(htName))
             if (attrFound) then
                ncStatus = nf90_get_att(ncid, nf90_global, trim(gAttrNames(iattr)), &
                     cldTbl(iHydr)%pProperty(ipp)%intplMethod)
                call cleanStr(cldTbl(iHydr)%pProperty(ipp)%intplMethod)
                if(ncStatus /= nf90_NoErr) call handle_err(ncStatus)
                cycle
             end if
             attrFound = matchName(gAttrNames(iattr),trim(cldTbl(iHydr)%pProperty(ipp)%name), &
                  trim(uGridPrefix)//trim(htName))
             if (attrFound) then
                ncStatus = nf90_get_att(ncid, nf90_global, trim(gAttrNames(iattr)), &
                     cldTbl(iHydr)%pProperty(ipp)%unifGrid)
                if(ncStatus /= nf90_NoErr) call handle_err(ncStatus)
                cycle
             end if

             attrFound = matchName(gAttrNames(iattr),trim(cldTbl(iHydr)%pProperty(ipp)%name), &
                  trim(pGridPrefix)//trim(htName))
             if (attrFound) then
                allocate( buffer1D( size( cldTbl(iHydr)%pProperty(ipp)%grid) ) )

                ncStatus = nf90_get_att(ncid, nf90_global, trim(gAttrNames(iattr)), &
                              buffer1D)

                if(ncStatus /= nf90_NoErr) call handle_err(ncStatus)
                cldTbl(iHydr)%pProperty(ipp)%grid = buffer1D
                deallocate(buffer1D)
                cycle
             end if

             ! units contain the same substrings for physical
             ! properties as in the dimension names, which are
             ! different substrings from the rest of the global
             ! attributes
             iy = index(gAttrNames(iattr),trim(unitPrefix)//trim(cldTbl(iHydr)%pProperty(ipp)%name))
             if (iy /=0) then
                ncStatus = nf90_get_att(ncid, nf90_global, trim(gAttrNames(iattr)), &
                     cldTbl(iHydr)%pProperty(ipp)%units)
                if(ncStatus /= nf90_NoErr) call handle_err(ncStatus, 'warning: units not found')
                cycle
             end if
          end do
       end do
       ncStatus = nf90_inquire_dimension(ncid,mydimid(nDims),name=spcName,len=mynSpcPts)
       if(ncStatus /= nf90_NoErr) call handle_err(ncStatus)

       if (.not. allocated(mySpectPtsGrd)) then
          allocate(mySpectPtsGrd(mynSpcPts))
          allocate( buffer1D(mynSpcPts) )

          ncStatus = nf90_get_att(ncid, nf90_global, 'spectralPoints', buffer1D)
          if(ncStatus /= nf90_NoErr) call handle_err(ncStatus)
          mySpectPtsGrd = buffer1D
          nSPsave=mynSpcPts
          deallocate( buffer1D)
       else
          if (mynSpcPts /= nSPsave) then
             print *,'err[CloudTab::loadCloudTab]: '// &
                  'Inconsistent number of spectral points'
             print *,'mynSpcPts /= nSPsave:',mynSpcPts,nSPsave
             call exit(1)
          endif
       endif

       cldTbl(iHydr)%noProperty = noProperty
       do iop = 1, noProperty
          iz = index(tblNames(itbl),trim(htName)//'_'//trim(loadedOptProp(iop)))
          if (iz  /= 0) exit
       end do
       cldTbl(iHydr)%oProperty(iop)%name = trim(tblNames(itbl))
       cldOptSet(iHydr)%property(iop)%name = trim(tblNames(itbl))

       ncStatus = nf90_get_att(ncid, tmpvarid, "units", cldTbl(iHydr)%oProperty(iop)%units)
       if (ncStatus /= nf90_NoErr) then
          call handle_err(ncStatus,'warning: units not found')
       else
          cldOptSet(iHydr)%property(iop)%units = trim(cldTbl(iHydr)%oProperty(iop)%units)
       end if
       if (.not. allocated(cldOptSet(iHydr)%property(iop)%values)) &
            allocate(cldOptSet(iHydr)%property(iop)%values(mynSpcPts))

       if (npProperty >= 1) then
          if (.not. allocated(cldOptSet(iHydr)%property(iop)%deriv)) &
               allocate(cldOptSet(iHydr)%property(iop)%deriv(npProperty,mynSpcPts), stat=istatus)
          if (istatus /=0 ) then
             print *,'err[CloudTab::loadCloudTab]: ', &
                  'Not enough memory for CloudTab::cldOptSet%property%deriv; istatus:',istatus
             call exit(1)
          end if
       end if

       if (npProperty >= 1) then
          if (.not. allocated(cldTbl(iHydr)%oProperty(iop)%tbl2)) &
               allocate(cldTbl(iHydr)%oProperty(iop)%tbl2(pdimArr(npProperty),mynSpcPts),stat=istatus)
          if (istatus /=0 ) then
             print *,'err[CloudTab::loadCloudTab]: ', &
                  'Not enough memory for CloudTab::cldTbl%oProperty; istatus:',istatus
             call exit(1)
          end if
       end if
       if (npProperty >= 2) then
          if (.not. allocated(cldTbl(iHydr)%oProperty(iop)%tbl3)) &
               allocate(cldTbl(iHydr)%oProperty(iop)%tbl3(pdimArr(npProperty-1),pdimArr(npProperty), &
               mynSpcPts), stat=istatus)
          if (istatus /=0 ) then
             print *,'err[CloudTab::loadCloudTab]: ', &
                  'Not enough memory for CloudTab::cldTbl%oProperty; istatus:',istatus
             call exit(1)
          endif
       end if
       if (npProperty >= 3) then
          if (.not. allocated(cldTbl(iHydr)%oProperty(iop)%tbl4)) &
               allocate(cldTbl(iHydr)%oProperty(iop)%tbl4(pdimArr(npProperty-2),pdimArr(npProperty-1), &
               pdimArr(npProperty),mynSpcPts), stat=istatus)
          if (istatus /=0 ) then
             print *,'err[CloudTab::loadCloudTab]: ', &
                  'Not enough memory for CloudTab::cldTbl%oProperty; istatus:',istatus
             call exit(1)
          endif
       end if
       if (npProperty >= 4) then
          if (.not. allocated(cldTbl(iHydr)%oProperty(iop)%tbl5)) &
               allocate(cldTbl(iHydr)%oProperty(iop)%tbl5(pdimArr(npProperty-3), &
               pdimArr(npProperty-2),pdimArr(npProperty-1),pdimArr(npProperty),mynSpcPts),stat=istatus)
          if (istatus /=0 ) then
             print *,'err[CloudTab::loadCloudTab]: ', &
                  'Not enough memory for CloudTab::cldTbl%oProperty; istatus:',istatus
             call exit(1)
          endif
       end if
       if (npProperty >= 5) then
          if (.not. allocated(cldTbl(iHydr)%oProperty(iop)%tbl6)) &
               allocate(cldTbl(iHydr)%oProperty(iop)%tbl6(pdimArr(npProperty-4),pdimArr(npProperty-3), &
               pdimArr(npProperty-2),pdimArr(npProperty-1),pdimArr(npProperty),mynSpcPts), stat=istatus)
          if (istatus /=0 ) then
             print *,'err[CloudTab::loadCloudTab]: ', &
                  'Not enough memory for CloudTab::cldTbl%oProperty; istatus:',istatus
             call exit(1)
          endif
       end if

       ! load table according to the number of indepdenent physical
       ! properties available
       optTable: select case (npProperty)
       case(1)
          allocate(buffer2D( pdimArr(1), mynSpcPts) ,stat=istatus)
          ncStatus = nf90_get_var(ncid, tmpvarid, buffer2D)
          if(ncStatus /= nf90_NoErr) call handle_err(ncStatus,'[error in CloudTab::loadCloudTab]: tbl2')
          cldTbl(iHydr)%oProperty(iop)%tbl2 = buffer2D
          deallocate(buffer2D)
       case(2)
          allocate(buffer3D( pdimArr(1), pdimArr(2),mynSpcPts) ,stat=istatus)
          ncStatus = nf90_get_var(ncid, tmpvarid, buffer3D)
          if(ncStatus /= nf90_NoErr) call handle_err(ncStatus,'[error in CloudTab::loadCloudTab]: tbl3')
          cldTbl(iHydr)%oProperty(iop)%tbl3 = buffer3D
          deallocate(buffer3D)
       case(3)
          allocate(buffer4D( pdimArr(1), pdimArr(2), pdimArr(3), mynSpcPts) ,stat=istatus)
          ncStatus = nf90_get_var(ncid, tmpvarid, buffer4D)
          if(ncStatus /= nf90_NoErr) call handle_err(ncStatus,'[error in CloudTab::loadCloudTab]: tbl4')
          cldTbl(iHydr)%oProperty(iop)%tbl4 = buffer4D
          deallocate(buffer4D)
       case(4)
          allocate(buffer5D( pdimArr(1), pdimArr(2), pdimArr(3), pdimArr(4), mynSpcPts) ,stat=istatus)
          ncStatus = nf90_get_var(ncid, tmpvarid, buffer5D)
          if(ncStatus /= nf90_NoErr) call handle_err(ncStatus,'[error in CloudTab::loadCloudTab]: tbl5')
          cldTbl(iHydr)%oProperty(iop)%tbl5 = buffer5D
          deallocate(buffer5D)
       case(5)
          allocate(buffer6D( pdimArr(1), pdimArr(2), pdimArr(3), pdimArr(4), pdimArr(5), mynSpcPts) ,stat=istatus)
          ncStatus = nf90_get_var(ncid, tmpvarid, buffer6D)
          if(ncStatus /= nf90_NoErr) call handle_err(ncStatus,'[error in CloudTab::loadCloudTab]: tbl6')
          cldTbl(iHydr)%oProperty(iop)%tbl6 = buffer6D
          deallocate(buffer6D)
       case default
          call handle_err(99,"[error]: reached max. number of physical properties in cloud optical table")
       end select optTable

       ! Check the interpolation method for each physical property
       do ipp = 1, npProperty
          call interpTypeCheck(trim(htName)//trim(cldTbl(iHydr)%pProperty(ipp)%name), &
               cldTbl(iHydr)%pProperty(ipp)%size, &
               cldTbl(iHydr)%pProperty(ipp)%intplMethod)
       end do

       ! Compute the interpolation step for each physical property
       do ipp = 1, npProperty
          cldTbl(iHydr)%pProperty(ipp)%step = &
               stepUniform(trim(htName)//trim(cldTbl(iHydr)%pProperty(ipp)%name), &
               cldTbl(iHydr)%pProperty(ipp)%unifGrid, &
               cldTbl(iHydr)%pProperty(ipp)%grid)
       end do

       initialized(iHydr)=.true.

    end do

    ncStatus = nf90_close(ncid)
    if (ncStatus /= nf90_NoErr) call handle_err(ncStatus)

    if (allocated(tblNames)) deallocate(tblNames)
    if (allocated(gAttrNames)) deallocate(gAttrNames)
    if (allocated(pdimArr)) deallocate(pdimArr)
    if (allocated(mydimid)) deallocate(mydimid)
    loadDone=.true.
    return

  END SUBROUTINE loadCloudTab

  !--------------------------------------------------------------------------

  SUBROUTINE interpTypeCheck(name,nelem,type)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   interpTypeCheck
!
! PURPOSE:
!
!   Check that the software can handle the interpolation type for which the data
!   were stored.
!
! SYNTAX:
!
!   CALL interpTypeCheck(name, nelem, type)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   name   CHAR     Generic name
!   nelem  INTEGER  Number of elements
!
!   INPUTS/OUTPUTS:
!
!   type   CHAR     Type identifier
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>


    character(len=*), intent(in)    :: name
    integer,          intent(in)    :: nelem
    character(len=*), intent(inout) :: type
    integer                         :: itype

    character(len=mxNameLength)     :: tmptype
    integer             :: i

    itype=0
    if (nelem == 1) then
       type='none'      ! Override interpolation type when only one node
       print *,'msg[CloudTab::interpTypeCheck]: '
       print *,trim(name),' grid has only 1 node; will be assumed constant'
    endif

    ! Get integer index of interpolation type
    do i=1,size(validInterpTypes)
       tmptype=validInterpTypes(i)
       if (trim(adjustl(type)) == trim(adjustl(tmptype))) then
          itype=i
          exit
       endif
    enddo
    if (itype == 0) then
       print *,'err[CloudTab::interpTypeCheck]: '
       print *,'For '//trim(name)//', interp type "'//trim(adjustl(type))// &
            '" not in valid types:'
       do i=1,size(validInterpTypes)
          print *,trim(validInterpTypes(i))
       end do
       call exit(1)
    endif

    return

  END SUBROUTINE interpTypeCheck

  !--------------------------------------------------------------------------

  FUNCTION stepUniform(name,unig,grid)

!<f90Function>**********************************************************
!
! NAME:
!
!   stepUniform
!
! PURPOSE:
!
!   Compute the grid step size, expecting that spacing is uniform.
!
! SYNTAX:
!
!   Results=stepUniform(name, unig, grid)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   name         CHAR     Generic name
!   unig         INTEGER  Integer flag for uniform grid
!   grid         REAL     Grid nodes
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


    character(len=*),   intent(in)  :: name
    integer,            intent(in)  :: unig
    real, dimension(:), intent(in)  :: grid
    real                            :: stepUniform

    real, parameter :: deplim=0.001
    integer         :: i,ngrid
    real            :: step,dep

    stepUniform=0.
    if (unig /= 1) return
    ngrid=size(grid)
    if (ngrid == 1) return

    step=(grid(ngrid)-grid(1))/float(ngrid-1)

    do i=2,size(grid)
       dep=abs(grid(i)-grid(1)-step*float(i-1))/step
       if (dep > deplim) then
          print*,grid
          print *,'err[CloudTab::stepUniform]: ', &
               name,' grid claimed to be uniform but was not; i,dep:',i,dep
          call exit(1)
       endif
    enddo

    stepUniform=step
    return

  END FUNCTION stepUniform

  !--------------------------------------------------------------------------

  SUBROUTINE getCloudTab(hydrName, cldSet, cldOpt)

    !<f90Subroutine>**************************d******************************
    !
    ! NAME:
    !
    !   getCloudTab
    !
    ! PURPOSE:
    !
    !   Provide cloud property data for all spectral points interpolated to the
    !   specified values of temperature and primary and secondary parameters of the
    !   particle size distribution model.
    !
    ! SYNTAX:
    !
    !   CALL getCloudTab(hydrName, cldSet, cldOpt)
    !
    ! ARGUMENTS:
    !
    !   INPUTS:
    !
    !   hydrName  character               The name of the hydrometeor to be processed
    !   cldSet    TYPE(CloudPhysicalSet)  All cloud physical properties at a point in space
    !
    !   OUTPUTS:
    !
    !   cldOpt  TYPE(CloudOpticalProperty)
    !                                     Cloud optical properties at all spectral points
    ! INCLUDES:
    !
    !   None
    !
    !*******************************************************</f90Subroutine>

    CHARACTER (LEN=*), INTENT(IN)                         :: hydrName
    TYPE(CloudPhysicalSet), INTENT(IN)                    :: cldSet
    TYPE (CloudOpticalProperty), dimension(:), allocatable, intent(INOUT) :: cldOpt

    integer, dimension(mxPhysicalProperty) :: lowerB, upperB
    integer :: iHydr, ipp, iop
    integer :: currentHydr
    integer :: lowerBound, upperBound
    integer :: ix
    character (len=*), parameter :: procName = " [CloudTab::getCloudTab]: "

    do iHydr = 1, size(cldTbl)
       if (.not. initialized(iHydr)) then
          print *, procName, "err - uninitialized output hydrometeor table "&
                 //trim(cldTbl(iHydr)%name)
          call exit(12)
       end if
    end do

    ! Compute the index for the matching hydrometeor type
    currentHydr = -1
    do iHydr = 1, size(cldTbl)
       ix = index(cldTbl(iHydr)%name, trim(hydrName))
       if (ix /= 0) then
          currentHydr = iHydr
          exit
       end if
    end do

    if (currentHydr == -1) then
         print *, procName,"[error]: "//trim(hydrName)//" not loaded in the table"
         call exit(16)
    end if

    if (cldSet%nProperty /= cldTbl(currentHydr)%npProperty) then
       if (cldSet%nProperty > cldTbl(currentHydr)%npProperty) then
          print *, procName, 'err - requested number of physical properties ',&
                'greater than provided (', cldSet%nProperty, &
                ' vs. ', cldTbl(currentHydr)%npProperty, ')'
          call exit(99)
       end if
       print *, procName, 'warning - number of physical properties ', &
            'differ b/w scene and table (', cldSet%nProperty, &
            ' vs. ', cldTbl(currentHydr)%npProperty, ')'
       call exit(88)
    end if

    do iop = 1, cldOptSet(currentHydr)%nProperty
       if (size(cldOptSet(currentHydr)%property(iop)%values) /= mynSpcPts) then
          print *, procName, 'err - mismatch number of spectral points for ', &
               trim(cldOptSet(currentHydr)%property(iop)%name)
          call exit(14)
       endif
    end do

    ! Find bounding grid indices
    do ipp = 1, cldSet%nProperty
       call boundGrid(cldTbl(currentHydr)%pProperty(ipp)%name, &
            cldTbl(currentHydr)%pProperty(ipp)%grid, &
            cldTbl(currentHydr)%pProperty(ipp)%unifGrid, &
            cldTbl(currentHydr)%pProperty(ipp)%step, &
            cldSet%value(indMap(currentHydr,ipp)),lowerBound,upperBound)
       lowerB(ipp) = lowerBound
       upperB(ipp) = upperBound
    end do

    do iop = 1, cldTbl(currentHydr)%noProperty
       calcOptTable: select case (cldSet%nProperty)
       case (1)
          call calcCloudTab2D(lowerB,upperB,iop,currentHydr,cldSet)
       case (2)
          call calcCloudTab3D(lowerB,upperB,iop,currentHydr,cldSet)
       case (3)
          call calcCloudTab4D(lowerB,upperB,iop,currentHydr,cldSet)
       case (4)
          call calcCloudTab5D(lowerB,upperB,iop,currentHydr,cldSet)
       case (5)
          call calcCloudTab6D(lowerB,upperB,iop,currentHydr,cldSet)
       case default
          call handle_err(99, trim(procName)//"error - calcCloudTab not defined when greater than 6-D")
       end select calcOptTable
       cldOpt(iop)%values = cldOptSet(currentHydr)%property(iop)%values
       cldOpt(iop)%units = trim(cldOptSet(currentHydr)%property(iop)%units)
       cldOpt(iop)%deriv = cldOptSet(currentHydr)%property(iop)%deriv
    end do

  END SUBROUTINE getCloudTab

  SUBROUTINE calcCloudTab2D(lowerB,upperB,iop,iHydr,cldSet)
    integer, dimension(:), intent(IN)  :: lowerB, upperB
    integer, intent(in)                :: iop, iHydr
    integer                            :: ix
    TYPE(CloudPhysicalSet), INTENT(IN) :: cldSet

    ix = cldSet%nProperty
    call interpolate(trim(cldTbl(iHydr)%pProperty(ix)%intplMethod), &
         cldTbl(iHydr)%pProperty(ix)%grid, &
         lowerB(ix),upperB(ix), &
         cldSet%value(indMap(iHydr,ix)), &  !replaced cldTbl(iHydr)%pProperty(ix)%grid
         cldTbl(iHydr)%oProperty(iop)%tbl2(:,:), &
         cldOptSet(iHydr)%property(iop)%values(:), &
         cldOptSet(iHydr)%property(iop)%deriv(ix,:))

  END SUBROUTINE calcCloudTab2D

  SUBROUTINE calcCloudTab3D(lowerB,upperB,iop,iHydr,cldSet)
    integer, dimension(:), intent(IN) :: lowerB, upperB
    integer, intent(in) :: iop, iHydr
    integer :: ib, ix
    TYPE(CloudPhysicalSet), INTENT(IN) :: cldSet
    character (len=*), parameter       :: procName = " [CloudTab::calcCloudTab3D]: "
    ix = cldSet%nProperty
    do ib = lowerB(ix), upperB(ix)
       call interpolate(trim(cldTbl(iHydr)%pProperty(ix-1)%intplMethod), &
            cldTbl(iHydr)%pProperty(ix-1)%grid, &
            lowerB(ix-1),upperB(ix-1), &
            cldSet%value(indMap(iHydr,ix-1)), &
            cldTbl(iHydr)%oProperty(iop)%tbl3(:,ib,:), &
            cldTbl(iHydr)%oProperty(iop)%tbl2(ib,:), &
            cldOptSet(iHydr)%property(iop)%deriv(ix-1,:))
    end do
    call interpolate(trim(cldTbl(iHydr)%pProperty(ix)%intplMethod), &
         cldTbl(iHydr)%pProperty(ix)%grid, &
         lowerB(ix),upperB(ix), &
         cldSet%value(indMap(iHydr,ix)), &
         cldTbl(iHydr)%oProperty(iop)%tbl2(:,:), &
         cldOptSet(iHydr)%property(iop)%values(:), &
         cldOptSet(iHydr)%property(iop)%deriv(ix,:))

  END SUBROUTINE calcCloudTab3D

  SUBROUTINE calcCloudTab4D(lowerB,upperB,iop,iHydr,cldSet)
    integer, dimension(:), intent(IN) :: lowerB, upperB
    integer, intent(IN) :: iop, iHydr
    integer :: ib1, ib2, ix
    TYPE(CloudPhysicalSet), INTENT(IN) :: cldSet
    character (len=*), parameter       :: procName= " [CloudTab::calcCloudTab4D]: "

    ix = cldSet%nProperty
    do ib2 = lowerB(ix),upperB(ix)
       do ib1 = lowerB(ix-1), upperB(ix-1)
          call interpolate(trim(cldTbl(iHydr)%pProperty(ix-2)%intplMethod), &
               cldTbl(iHydr)%pProperty(ix-2)%grid, &
               lowerB(ix-2),upperB(ix-2), &
               cldSet%value(indMap(iHydr,ix-2)), &
               cldTbl(iHydr)%oProperty(iop)%tbl4(:,ib1,ib2,:), &
               cldTbl(iHydr)%oProperty(iop)%tbl3(ib1,ib2,:), &
               cldOptSet(iHydr)%property(iop)%deriv(ix-2,:))
       end do
       call interpolate(trim(cldTbl(iHydr)%pProperty(ix-1)%intplMethod), &
            cldTbl(iHydr)%pProperty(ix-1)%grid, &
            lowerB(ix-1),upperB(ix-1), &
            cldSet%value(indMap(iHydr,ix-1)), &
            cldTbl(iHydr)%oProperty(iop)%tbl3(:,ib2,:), &
            cldTbl(iHydr)%oProperty(iop)%tbl2(ib2,:), &
            cldOptSet(iHydr)%property(iop)%deriv(ix-1,:))
    end do
    call interpolate(trim(cldTbl(iHydr)%pProperty(ix)%intplMethod), &
         cldTbl(iHydr)%pProperty(ix)%grid, &
         lowerB(ix),upperB(ix), &
         cldSet%value(indMap(iHydr,ix)), &
         cldTbl(iHydr)%oProperty(iop)%tbl2(:,:), &
         cldOptSet(iHydr)%property(iop)%values(:), &
         cldOptSet(iHydr)%property(iop)%deriv(ix,:))

  END SUBROUTINE calcCloudTab4D

  SUBROUTINE calcCloudTab5D(lowerB,upperB,iop,iHydr,cldSet)
    integer, dimension(:), intent(IN) :: lowerB, upperB
    integer, intent(IN) :: iop, iHydr
    integer :: ib1, ib2, ib3, ix
    TYPE(CloudPhysicalSet), INTENT(IN) :: cldSet

    ix = cldSet%nProperty
    do ib3 = lowerB(ix), upperB(ix)
       do ib2 = lowerB(ix-1), upperB(ix-1)
          do ib1 = lowerB(ix-2), upperB(ix-2)
             call interpolate(trim(cldTbl(iHydr)%pProperty(ix-3)%intplMethod), &
                  cldTbl(iHydr)%pProperty(ix-3)%grid, &
                  lowerB(ix-3),upperB(ix-3), &
                  cldSet%value(indMap(iHydr,ix-3)), &
                  cldTbl(iHydr)%oProperty(iop)%tbl5(:,ib1,ib2,ib3,:), &
                  cldTbl(iHydr)%oProperty(iop)%tbl4(ib1,ib2,ib3,:), &
                  cldOptSet(iHydr)%property(iop)%deriv(ix-3,:))
          end do
          call interpolate(trim(cldTbl(iHydr)%pProperty(ix-2)%intplMethod), &
               cldTbl(iHydr)%pProperty(ix-2)%grid, &
               lowerB(ix-2),upperB(ix-2), &
               cldSet%value(indMap(iHydr,ix-2)), &
               cldTbl(iHydr)%oProperty(iop)%tbl4(:,ib2,ib3,:), &
               cldTbl(iHydr)%oProperty(iop)%tbl3(ib2,ib3,:), &
               cldOptSet(iHydr)%property(iop)%deriv(ix-2,:))
       end do
       call interpolate(trim(cldTbl(iHydr)%pProperty(ix-1)%intplMethod), &
            cldTbl(iHydr)%pProperty(ix-1)%grid, &
            lowerB(ix-1),upperB(ix-1), &
            cldSet%value(indMap(iHydr,ix-1)), &
            cldTbl(iHydr)%oProperty(iop)%tbl3(:,ib3,:), &
            cldTbl(iHydr)%oProperty(iop)%tbl2(ib3,:), &
            cldOptSet(iHydr)%property(iop)%deriv(ix-1,:))
    enddo
    call interpolate(trim(cldTbl(iHydr)%pProperty(ix)%intplMethod), &
         cldTbl(iHydr)%pProperty(ix)%grid, &
         lowerB(ix),upperB(ix), &
         cldSet%value(indMap(iHydr,ix)), &
         cldTbl(iHydr)%oProperty(iop)%tbl2(:,:), &
         cldOptSet(iHydr)%property(iop)%values(:), &
         cldOptSet(iHydr)%property(iop)%deriv(ix,:))

  END SUBROUTINE calcCloudTab5D

  SUBROUTINE calcCloudTab6D(lowerB,upperB,iop,iHydr,cldSet)
    integer, dimension(:), intent(IN) :: lowerB, upperB
    integer, intent(IN) :: iop, iHydr
    integer :: ib1, ib2, ib3, ib4, ix
    TYPE(CloudPhysicalSet), INTENT(IN) :: cldSet

    ix = cldSet%nProperty
    do ib4 = lowerB(ix), upperB(ix)
       do ib3 = lowerB(ix-1), upperB(ix-1)
          do ib2 = lowerB(ix-2), upperB(ix-2)
             do ib1 = lowerB(ix-3), upperB(ix-3)
                call interpolate(trim(cldTbl(iHydr)%pProperty(ix-4)%intplMethod), &
                     cldTbl(iHydr)%pProperty(ix-4)%grid, &
                     lowerB(ix-4),upperB(ix-4), &
                     cldSet%value(indMap(iHydr,ix-4)), &
                     cldTbl(iHydr)%oProperty(iop)%tbl6(:,ib1,ib2,ib3,ib4,:), &
                     cldTbl(iHydr)%oProperty(iop)%tbl5(ib1,ib2,ib3,ib4,:), &
                     cldOptSet(iHydr)%property(iop)%deriv(ix-4,:))
             end do
             call interpolate(trim(cldTbl(iHydr)%pProperty(ix-3)%intplMethod), &
                  cldTbl(iHydr)%pProperty(ix-3)%grid, &
                  lowerB(ix-3),upperB(ix-3), &
                  cldSet%value(indMap(iHydr,ix-3)), &
                  cldTbl(iHydr)%oProperty(iop)%tbl5(:,ib2,ib3,ib4,:), &
                  cldTbl(iHydr)%oProperty(iop)%tbl4(ib2,ib3,ib4,:), &
                  cldOptSet(iHydr)%property(iop)%deriv(ix-3,:))
          end do
          call interpolate(trim(cldTbl(iHydr)%pProperty(ix-2)%intplMethod), &
               cldTbl(iHydr)%pProperty(ix-2)%grid, &
               lowerB(ix-2),upperB(ix-2), &
               cldSet%value(indMap(iHydr,ix-2)), &
               cldTbl(iHydr)%oProperty(iop)%tbl4(:,ib3,ib4,:), &
               cldTbl(iHydr)%oProperty(iop)%tbl3(ib3,ib4,:), &
               cldOptSet(iHydr)%property(iop)%deriv(ix-2,:))
       enddo
       call interpolate(trim(cldTbl(iHydr)%pProperty(ix-1)%intplMethod), &
            cldTbl(iHydr)%pProperty(ix-1)%grid, &
            lowerB(ix-1),upperB(ix-1), &
            cldSet%value(indMap(iHydr,ix-1)), &
            cldTbl(iHydr)%oProperty(iop)%tbl3(:,ib4,:), &
            cldTbl(iHydr)%oProperty(iop)%tbl2(ib4,:), &
            cldOptSet(iHydr)%property(iop)%deriv(ix-1,:))
    enddo
    call interpolate(trim(cldTbl(iHydr)%pProperty(ix)%intplMethod), &
         cldTbl(iHydr)%pProperty(ix)%grid, &
         lowerB(ix),upperB(ix), &
         cldSet%value(indMap(iHydr,ix)), &
         cldTbl(iHydr)%oProperty(iop)%tbl2(:,:), &
         cldOptSet(iHydr)%property(iop)%values(:), &
         cldOptSet(iHydr)%property(iop)%deriv(ix,:))

  END SUBROUTINE calcCloudTab6D

  SUBROUTINE boundGrid(name,grid,unig,step,val,i1,i2)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   boundGrid
!
! PURPOSE:
!
!   Find the indices of a physical quantity grid that bound its given value.
!
! SYNTAX:
!
!   CALL boundGrid(name, grid, unig, step, val, i1, i2)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   name  CHAR     Generic name
!   grid  REAL     Grid nodes
!   unig  INTEGER  Integer flag for uniform grid
!   step  REAL     Uniform grid spacing
!   val   REAL     Data value
!
!   INPUTS/OUTPUTS:
!
!   i1    INTEGER  Grid index
!   i2    INTEGER  Grid index
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>


  ! Find the indices of the grid that bound the given value

    character(len=*),   intent(in)    :: name
    real, dimension(:), intent(in)    :: grid
    integer,            intent(in)    :: unig
    real,               intent(in)    :: step
    real,               intent(in)    :: val
    integer,            intent(inout) :: i1
    integer,            intent(inout) :: i2

    integer :: i,ngrid

    ngrid=size(grid)

    ! Grid has size 1
    if (ngrid == 1) then
      i1=1
      i2=1
    ! Grid is uniformly spaced; saves time without loop over if block
    elseif (unig == 1) then
      i1=floor((val-grid(1))/step)+1
      i2=i1+1
      if (i1 < 1) then
        if (outBndsExit) then
          print *,'err[CloudTab::boundGrid]: ', &
          name,' value',val,' is out of uniform range:',grid(1),grid(ngrid)
          call exit(1)
        elseif (val < 0.) then
          !print *,'err[CloudTab::boundGrid]: ', &
          !name,' value',val,' is unphysical negative; min grid=',grid(1)
          !call exit(1)
          i1=1
          i2=2
        else
          i1=1
          i2=2
        endif
      endif
      if (i2 > ngrid) then
        if (val <= grid(ngrid)) then  ! Accommodate roundoff
          i2=ngrid
          i1=i2-1
        elseif (outBndsExit) then
          print *,'err[CloudTab::boundGrid]: ', &
          name,' value',val,' is out of uniform range:',grid(1),grid(ngrid)
          call exit(1)
        else
          i2=ngrid
          i1=i2-1
        endif
      endif
    else
    ! Grid is not uniformly spaced
      if ((val < grid(1) .or. val > grid(ngrid)) .and. outBndsExit) then
        print *,'err[CloudTab::boundGrid]: ', &
        name,' value',val,' is out of range:',grid(1),grid(ngrid)
        call exit(1)
      endif

      do i=2,ngrid
         if (grid(i) >= val) exit
      enddo
      i2=min(i,ngrid)
      i1=i2-1

   end if

   return

  END SUBROUTINE boundGrid

  !--------------------------------------------------------------------------

  SUBROUTINE interpolate(type,grid,i1,i2,xin,arr,arrIntrp,arrDer)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   interpolate
!
! PURPOSE:
!
!   Interpolate from grid, according the specified method.
!
! SYNTAX:
!
!   CALL interpolate(type, grid, i1, i2, xin, arr, arrIntrp, arrDer)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   type      CHAR     Type identifier
!   grid      REAL     Grid nodes
!   i1        INTEGER  Grid index
!   i2        INTEGER  Grid index
!   xin       REAL     Interpolation point
!   arr       REAL     Array to be processed
!
!   INPUTS/OUTPUTS:
!
!   arrIntrp  REAL     Interpolated array
!   arrDer*   REAL     Array derivatives within interpolation range
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>


    character(len=*),     intent(in)              :: type
    real, dimension(:),   intent(in)              :: grid
    integer,              intent(in)              :: i1
    integer,              intent(in)              :: i2
    real,                 intent(in)              :: xin
    real, dimension(:,:), intent(in)              :: arr
    real, dimension(:),   intent(inout)           :: arrIntrp
    real, dimension(:),   intent(inout), optional :: arrDer
    real, dimension(size(arr,2))                  :: slope

    select case(type)
    case ('none')
       arrIntrp=arr(i1,:)
       if(present(arrDer))arrDer=0.
    case ('linear')
       if (xin>grid(i2)) then
          slope = 0.
          arrIntrp=arr(i2,:)
       else if (xin<grid(i1)) then
          slope = 0.
          arrIntrp=arr(i1,:)
       else
         slope = (arr(i2,:)-arr(i1,:))/(grid(i2)-grid(i1))
         arrIntrp=arr(i1,:)+(xin-grid(i1))*slope
       endif
!       uncomment the lines if extrapolation is desired
!       slope = (arr(i2,:)-arr(i1,:))/(grid(i2)-grid(i1))
!       arrIntrp=arr(i1,:)+(xin-grid(i1))*slope
       if(present(arrDer))arrDer = slope
    case default
       print *,'err[CloudTab::interpolate]: ', &
            'Unsupported interpolation method: ',type
       call exit(1)
    end select

    return

  END SUBROUTINE interpolate

  !--------------------------------------------------------------------------

  SUBROUTINE destroyCloudTab()

    !<f90Subroutine>********************************************************
    !
    ! NAME:
    !
    !   destroyCloudTab
    !
    ! PURPOSE:
    !
    !   Deallocate memory and reset the default values
    !
    ! SYNTAX:
    !
    !   CALL destroyCloudTab()
    !
    ! ARGUMENTS:
    !
    !   INPUTS:
    !
    !     None
    !
    ! INCLUDES:
    !
    !   None
    !
    !*******************************************************</f90Subroutine>

    integer  :: iHydr
    integer  :: ipp, iop

    do iHydr = 1, nHydr
       do ipp = 1, cldTbl(iHydr)%npProperty
          if (allocated(cldTbl(iHydr)%pProperty(ipp)%grid)) deallocate(cldTbl(iHydr)%pProperty(ipp)%grid)
       end do
       if (allocated(cldTbl(iHydr)%pProperty)) deallocate(cldTbl(iHydr)%pProperty)
       do iop = 1, cldTbl(iHydr)%noProperty
          if (cldTbl(iHydr)%npProperty >= 1) then
             if (allocated(cldTbl(iHydr)%oProperty(iop)%tbl2)) deallocate(cldTbl(iHydr)%oProperty(iop)%tbl2)
          end if
          if (cldTbl(iHydr)%npProperty >= 2) then
             if (allocated(cldTbl(iHydr)%oProperty(iop)%tbl3)) deallocate(cldTbl(iHydr)%oProperty(iop)%tbl3)
          end if
          if (cldTbl(iHydr)%npProperty >= 3) then
             if (allocated(cldTbl(iHydr)%oProperty(iop)%tbl4)) deallocate(cldTbl(iHydr)%oProperty(iop)%tbl4)
          end if
          if (cldTbl(iHydr)%npProperty >= 4) then
             if (allocated(cldTbl(iHydr)%oProperty(iop)%tbl5)) deallocate(cldTbl(iHydr)%oProperty(iop)%tbl5)
          end if
          if (cldTbl(iHydr)%npProperty >= 5) then
             if (allocated(cldTbl(iHydr)%oProperty(iop)%tbl6)) deallocate(cldTbl(iHydr)%oProperty(iop)%tbl6)
          end if
       end do
       do iop = 1, cldOptSet(iHydr)%nProperty
          if (allocated(cldOptSet(iHydr)%property(iop)%values)) deallocate(cldOptSet(iHydr)%property(iop)%values)
          if (allocated(cldOptSet(iHydr)%property(iop)%deriv)) deallocate(cldOptSet(iHydr)%property(iop)%deriv)
       end do
       if (allocated(cldTbl(iHydr)%oProperty)) deallocate(cldTbl(iHydr)%oProperty)
       if (allocated(cldOptSet(iHydr)%property)) deallocate(cldOptSet(iHydr)%property)
    end do

    if (allocated(cldTbl)) deallocate(cldTbl)
    if (allocated(cldOptSet)) deallocate(cldOptSet)
    if (allocated(mySpectPtsGrd)) deallocate(mySpectPtsGrd)

    ! Reset to the initial values
    initialized = .false.
    loadDone = .false.
    first_init = .true.
    outBndsExit = .false.
    nSPsave = 0
    npProperty = mxPhysicalProperty
    noProperty = mxOpticalProperty

    return

  END SUBROUTINE destroyCloudTab

  !handle netcdf errors
  SUBROUTINE handle_err(status, msg)
    integer, intent (in) :: status
    character (len=*), intent(in), optional :: msg
    integer :: iwarn

    if(status /= nf90_noerr) then
       if (present(msg)) then
          iwarn = index(trim(nf90_strerror(status))//trim(msg), 'warning')
       else
          iwarn = index(trim(nf90_strerror(status)),'warning')
       end if
       print *, trim(nf90_strerror(status))
       if (iwarn /= 0) then
          print *, trim(msg)
       else
          call exit(1)
       end if
    else
       if (present(msg)) then
          iwarn = index(msg,'warning')
          if (iwarn /= 0) then
             print *, "[warning - ", status, "]: ", trim(msg)
          else
             print *, "[error - ", status, "]: ", trim(msg)
             call exit(1)
          end if
       else
          print *, "[info]: ", status
       end if
    end if

  END SUBROUTINE handle_err

  !Remove the null character from the string
  subroutine cleanStr(str)

    CHARACTER (len=*), intent(inout) :: str
    integer :: istr
    integer :: ascCode

    do istr = 1, len(str)
       ascCode = iachar(str(istr:istr))
       if (ascCode == 0) ascCode = 32
       str(istr:istr) = achar(ascCode)
    end do

  END subroutine cleanStr
END MODULE CloudTab
