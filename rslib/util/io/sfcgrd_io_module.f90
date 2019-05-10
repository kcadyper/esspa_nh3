! "$Id$"
module sfcgrd_io_module
  ! ---------------------------------------------------------------------
  !  Module name: sfcgrd_io_module
  !  Purpose: Surface (Earth) Gridded Data I/O
  !  Description: The module consists of the following functions,
  !               openSfcGrd()
  !               closeSfcGrd()
  !               getSfcGrd()
  !               putSfcGrd()
  !               querySfcGrd()
  !  NOTE: This version is a substantial revision of original. It supports 
  !        a different and more flexible data format and structure. This is
  !        NOT backward compatible with the original sfcgrd_io_module.
  ! 
  !  Var_Name            Type       Description
  !  --------            ----       -----------
  !  ncid                integer    file id
  !  file                string     file name
  !  nCol                integer    number of columns in grid (this tile only, not globally unless global tile)
  !  nRow                integer    number of rows in grid (this tile only, not globally unless global tile)
  !  nUnlim              integer    length of unlimited dimension (obtained from openNcdfFile)
  !  dimUnlim            integer    length of unlimited dimension (obtained from readNcdfDim)
  !  dimFixed            integer    length of  fixed dimension (obtained from readNcdfDim)
  !  ndimUnlim           integer    number of individual dims folded into file unlimited dimension (global attribute)
  !  dimUnlimDims        integer    vector of individual dim lengths folded into file unlimited dimension (global attribute)
  !  dimNamesUnlim       character  vector of dimension names for dimUnlimDims (global attribute)
  !  ndimFixed           integer    number of individual fixed dims (global attribute)
  !  dimFixedDims        integer    vector of individual fixed dim lengths (global attribute)
  !  dimNamesFixed       character  vector of dimension names for dimFixedDims (global attribute)
  !  dimUnlimName        character  name of unlimited dimension
  !  nTimeLevels         integer    number of time levels
  !  polarity            real       channel polarizations
  !  freq                real       channel frequencies (ghz)
  !  wvn                 real       channel wavenumbers (cm-1)
  !  eia                 real       earth incidence angles for entire data set (degrees)
  !  eaa                 real       earth azimuth angles for entire data set (degrees)
  !  casename            character  case identification string
  !  nchmw               integer    number of MW channels (for back compatibility)
  !  nchan               integer    number of channels
  !  status              character  directs routine to open a new or old netCDF file
  !  proj                character  map projection type
  !  reflat              real       map projection reference latitude (reserved for future use)
  !  reflon              real       map projection reference longitude (reserved for future use)
  !  latorigin           real       global map origin latitude
  !  lonorigin           real       global map origin longitude
  !  rowGridOffset       real       row offset (from 0,0 of local tile grid) to global map origin
  !  colGridOffset       real       column offset (from 0,0 of local tile grid) to global map origin
  !  scale               real       map scale (i.e. grid resolution)
  !  scaleunits          character  units of variable scale 
  !  re                  real       earth radius (units same as scaleunits, normally km)
  !  itile               integer    tile i (column) coordinate of this dataset (upper left corner column=0)
  !  jtile               integer    tile j (row) coordinate of this dataset (upper left corner row=0)
  !  nitile              integer    number of global tiles in i-direction 
  !  njtile              integer    number of global tiles in j-direction
  !  xorigin             real       starting x coordinate (e.g. for regular lat/lon grids)
  !  delx                real       x coordinate grid spacing (e.g. for regular lat/lon grids)
  !  xunits              character  units of variables xorigin and delx
  !  yorigin             real       starting y coordinate (e.g. for regular lat/lon grids)
  !  dely                real       y coordinate grid spacing (e.g. for regular lat/lon grids)
  !  xunits              character  units of variables yorigin and dely
  !  timeLevelIncrement  real       delta-time represented by a single time level
  !  timeLevelUnits      character  units of timeLevelIncrement
  !  sto_x3d1_i1         real       storage offset and scale for data array x3d1_i1 (vector length of 2)
  !  sto_x3d1_i2         real       storage offset and scale for data array x3d1_i2 (vector length of 2)
  !  sto_x3d1_r4         real       storage offset and scale for data array x3d1_r4 (vector length of 2)
  !  sto_x2d1_i1         real       storage offset and scale for data array x2d1_i1 (vector length of 2)
  !  sto_x2d1_i2         real       storage offset and scale for data array x2d1_i2 (vector length of 2)
  !  sto_x2d1_r4         real       storage offset and scale for data array x2d1_r4 (vector length of 2)
  !  NOTE: offset and scale conventions are such that to convert between 
  !        stored values and physical units: 
  !        x_physical = x_stored*scale+offset = x_stored*sto(2)+sto(1)
  !        x_stored = (x_physical-offset)/scale = (x_physical-sto(1))/sto(2)                          
  !
  !  iloc                integer    "i-coordinate" (1st array index) location of grid (normally grid column)
  !                                 to be read/written
  !  jloc                integer    "j-coordinate" (2nd array index) location of grid (normally grid row)
  !                                 to be read/written
  !  kloc                integer    "k-coordinate" (3rd array index) location of grid (normally time level)
  !                                 to be read/written
  !  xid                 character  gridded data ID. (e.g., PRIjul97)
  !  xidlen              integer    string length of gridded data ID. (e.g., PRIjul97)
  !  x3d1_i1             byte       gridded 3-d (normally ncol*nrow*ntimelevels,nvalspergrid) data (byte integer) 
  !  x3d1_i2             integer*2  gridded 3-d data (2-byte integer) 
  !  x3d1_r4             real       gridded 3-d data (4-byte float) 
  !  x2d1_i1             byte       gridded 2-d data (normally ncol*nrow*ntimelevels) (byte integer) 
  !  x2d1_i2             integer*2  gridded 2-d data (2-byte integer) 
  !  x2d1_r4             real       gridded 2-d data (4-byte float) 
  !  x3d1_i1_longname    character  long name of data stored in variable x3d1_i1
  !  x3d1_i2_longname    character  long name of data stored in variable x3d1_i2
  !  x3d1_r4_longname    character  long name of data stored in variable x3d1_r4
  !  x2d1_i1_longname    character  long name of data stored in variable x2d1_i1
  !  x2d1_i2_longname    character  long name of data stored in variable x2d1_i2
  !  x2d1_r4_longname    character  long name of data stored in variable x2d1_r4
  !  x3d1_i1_units       character  units of data stored in variable x3d1_i1
  !  x3d1_i2_units       character  units of data stored in variable x3d1_i2
  !  x3d1_r4_units       character  units of data stored in variable x3d1_r4
  !  x2d1_i1_units       character  units of data stored in variable x2d1_i1
  !  x2d1_i2_units       character  units of data stored in variable x2d1_i2
  !  x2d1_r4_units       character  units of data stored in variable x2d1_r4
  !  year                integer    valid year (2-dimensional, i.e. ncol,nrow)
  !  time                real       valid julian day, including fractional amount (2-dimensional)
  !  msg                 character  closing message 
  !  lgetsto             logical    logical mask defining which storage offset and scale
  !                                 attributes to read from file
  !  lputsto             logical    logical mask defining which storage and scale
  !                                 attributes to write to file
  !  lgetdata            logical    logical mask defining which of the 6 possible 
  !                                 datasets to read from file
  !  lputdata            logical    logical mask defining which of the 6 possible 
  !                                 datasets to write to file
  !  lgetyear            logical    logical for reading (or not) year dataset from file
  !  lputyear            logical    logical for writing (or not) year dataset to file
  !  lgetjulday          logical    logical for reading (or not) julday dataset from file
  !  lputjulday          logical    logical for writing (or not) julday dataset to file
  !  lgetyear1d          logical    logical for reading (or not) 1-d (one per time level) year dataset from file
  !  lputyear1d          logical    logical for writing (or not) 1-d (one per time level) year dataset to file
  !  lgetjulday1d        logical    logical for reading (or not) 1-d (one per time level) julday dataset from file
  !  lputjulday1d        logical    logical for writing (or not) 1-d (one per time level) julday dataset to file
  !  
  !  Externals:
  !  File_Name       Description
  !  ---------       -----------
  !  ncdf_module     base module for performing netCDF I/Os
  !  
  !  Developed by AER 2005
  !  Copyright, AER, Inc., 2005
  ! ------------------------------------------------------------
  use ncdf_module
  implicit none
  private

  public openSfcGrd, getSfcGrd, putSfcGrd, closeSfcGrd, querySfcGrd
  interface openSfcGrd
     module procedure openNcSfcGrd
  end interface

  interface querySfcGrd
     module procedure queryNcSfcGrd
  end interface

  interface getSfcGrd
     module procedure get4DSfcGrd
     module procedure get1DSfcGrd
  end interface

  interface putSfcGrd
     module procedure put4DSfcGrd
     module procedure put1DSfcGrd
  end interface

  interface closeSfcGRd
     module procedure closeNcGFile
  end interface
  ! ---
CONTAINS
  ! -
  subroutine openNcSfcGrd(ncid, file, creationDate, &
       nCol, nRow, nValsPerGrid, nTimeLevels, &
       nUnlim, &
       dimUnlim, dimFixed, &
       ndimUnlim, dimUnlimDims, dimNamesUnlim, &
       ndimFixed, dimFixedDims, dimNamesFixed, &
       dimUnlimName, &
       polarity, freq, wvn, eia, eaa, casename, nchmw, nchan, proj, reflat, reflon, &
       latorigin, lonorigin, rowGridOffset, colGridOffset, scale, scaleunits, re, &
       itile, jtile, nitile, njtile, xorigin, delx, xunits, yorigin, dely, yunits, &
       timeLevelIncrement, timeLevelUnits, &
       lgetsto, lputsto, &
       sto_x3d1_i1, sto_x3d1_i2, sto_x3d1_r4, &
       sto_x2d1_i1, sto_x2d1_i2, sto_x2d1_r4, &
       xidLen, status)

    ! - Dummy variables
    integer, intent(inout) :: ncid 
    character (len=*), intent (in) :: file
    character (len=*), intent (inout), optional :: creationDate
    integer, intent(inout), optional :: xidlen
    integer, intent(inout), optional :: nrow, ncol, nValsPerGrid, nTimeLevels, nUnlim
    real, intent(inout), optional :: reflat, reflon, latorigin, lonorigin
    real, intent(inout), optional :: xorigin, delx, yorigin, dely
    real, intent(inout), optional :: rowGridOffset, colGridOffset
    integer, intent(inout), optional :: itile, jtile, nitile, njtile
    integer, intent(inout), optional :: nchmw
    integer, intent(inout), optional :: nchan
    real, intent (inout), optional :: scale, re
    real, dimension(:), intent (inout), optional :: freq, eia, eaa
    real, dimension(:), intent (inout), optional :: wvn
    integer, dimension(:), intent (inout), optional :: polarity
    logical, dimension(:), intent(in), optional :: lputsto, lgetsto
    real, dimension(2), intent (inout), optional :: sto_x3d1_i1, sto_x3d1_i2, sto_x3d1_r4
    real, dimension(2), intent (inout), optional :: sto_x2d1_i1, sto_x2d1_i2, sto_x2d1_r4
    integer, intent(inout), optional :: dimUnlim, dimFixed
    integer, intent(inout), optional :: ndimUnlim, ndimFixed
    integer, dimension(:), intent(inout), optional :: dimUnlimDims, dimFixedDims

    real, intent(inout), optional :: timeLevelIncrement
    character (len=*), intent(inout), optional :: timeLevelUnits
    character (len=*), intent(inout), optional :: casename
    character (len=*), intent(inout), optional :: proj
    character (len=*), intent(inout), optional :: xunits, yunits, scaleunits
    character (len=*), intent(in), optional :: status
    character (len=*), dimension(:), intent(inout), optional :: dimNamesUnlim
    character (len=*), dimension(:), intent(inout), optional :: dimNamesFixed
    character (len=*), intent(inout), optional :: dimUnlimName


    ! NOTE: add other dummy arguments here as needed
    ! to support other global attributes


    ! - Local variables
    integer :: i
    character (len=MAX_NAME_LENGTH) :: lstatus
    integer, dimension(8) :: sysvalues
    character (len= 8) :: sysdate
    character (len=10) :: systime
    character (len= 5) :: syszone
    integer :: nchanloc
    logical, dimension(6) :: lputsto1, lgetsto1
    integer :: ndimUnlim1, ndimFixed1, dimUnlim1, dimFixed1, nUnlim1
    integer, allocatable, dimension(:) :: dimUnlimDims1, dimFixedDims1
    character(len=20), allocatable, dimension(:) :: dimNamesUnlim1, dimNamesFixed1
    character(len=50) :: dimUnlimName1
    character(len=1) :: sAttrNum


!    print *,'Begin openSfcGrd: status=',trim(status)
    if (.not. present(status)) then
       lstatus = 'old'
    else
       lstatus = trim(status)
    endif

! Open new file and write global attributes

    if (trim(lstatus) == 'new' .or. trim(lstatus) == 'NEW') then 

! Abort if required arguments not present for new file creation

       if(.not. present(nCol)) &
            stop 'err[openSfcGrd]: missing nCol'
       if(.not. present(nRow)) &
            stop 'err[openSfcGrd]: missing nRow'
       if(.not. present(nTimeLevels)) &
            stop 'err[openSfcGrd]: missing nTimeLevels'
       if(.not. present(nValsPerGrid)) &
            stop 'err[openSfcGrd]: missing nValsPerGrid'

! Define defaults if arguments not present
! Note: For now, keep it simple:
! force global dimension attributes to hard-coded values based on "ncol,nrow,ntimelevels" model
! Eventually, these could be user-defined with checks for consistency with actual variable dimensions, 
! readjusting dimension ordering if necessary, etc.

!       print *,'Define defaults if arguments not present'

       nDimUnlim1=3
       print *,'nDimUnlim1=',nDimUnlim1
       allocate(dimUnlimDims1(ndimUnlim1),dimNamesUnlim1(ndimUnlim1))
       nDimFixed1=1
       print *,'nDimFixed1=',nDimFixed1
       allocate(dimFixedDims1(ndimFixed1),dimNamesFixed1(ndimFixed1))
       dimUnlimDims1(1)=nCol
       dimUnlimDims1(2)=nRow
       dimUnlimDims1(3)=nTimeLevels
       dimNamesUnlim1(1)='nCol'
       dimNamesUnlim1(2)='nRow'
       dimNamesUnlim1(3)='nTimeLevels'
       dimFixedDims1(1)=nValsPerGrid
       dimNamesFixed1(1)='nValsPerGrid'

       if(.not. present(dimUnlimName))then
! construct output dimname and dims for unlimited dim
          dimUnlimName1=''
          do i=1,ndimUnlim1
             if(i .ne. ndimUnlim1)then 
                dimUnlimName1=trim(dimUnlimName1)//trim(dimNamesUnlim1(i))//'_'
             else
                dimUnlimName1=trim(dimUnlimName1)//trim(dimNamesUnlim1(i))
             endif
          enddo
       else
          dimUnlimName1=trim(dimUnlimName)
       endif
       print *,'dimUnlimName1=',dimUnlimName1

       lputsto1(:)=.TRUE.
       if(present(lputsto))lputsto1=lputsto
       if ((present(wvn) .or. present(freq) .or. &
            present(polarity) .or. present(eia) .or. present(eaa)) .and. &
            .not. (present(nchmw) .or. present(nchan))) &
            stop 'err[openSfcGrd]: missing nchan'

       call openNcdfFile(ncid, file, status='new',  &
            unlimited_dim_name=trim(dimUnlimName1))

       if (present(casename)) &
            call writeNcdfAttr(ncid,attrName='case', &
            attr=trim(casename))
       call date_and_time(sysdate, systime, syszone, sysvalues)
       call writeNcdfAttr(ncid,attrName='CreationTime', &
            attr=trim(sysdate)//trim(systime))

       print *,'Begin writing dimension attributes'
       call writeNcdfAttr(ncid,attrName='nDimUnlim',attr=nDimUnlim1)
       call writeNcdfAttr(ncid,attrName='dimUnlimDims',attr=dimUnlimDims1(1:nDimUnlim1))
       print *,'finished nDimUnlim and DimUnlimDims'
       deallocate(dimUnlimDims1)
       do i=1,ndimUnlim1
          write(sAttrNum,'(i1)')i
          call writeNcdfAttr(ncid,attrName='dimNamesUnlim'//sAttrNum,&
               attr=trim(dimNamesUnlim1(i)))
       enddo
       deallocate(dimNamesUnlim1)
       print *,'finished dimNamesUnlim1'

       call writeNcdfAttr(ncid,attrName='nDimFixed',attr=nDimFixed1)
       call writeNcdfAttr(ncid,attrName='dimFixedDims',attr=dimFixedDims1(1:nDimFixed1))
       deallocate(dimFixedDims1)
       print *,'finished nDimFixed and dimFixedDims'
       do i=1,ndimFixed1
          write(sAttrNum,'(i1)')i
          call writeNcdfAttr(ncid,attrName='dimNamesFixed'//sAttrNum,&
               attr=trim(dimNamesFixed1(i)))
       enddo
       deallocate(dimNamesFixed1)
       print *,'finished dimNamesFixed1'
       call writeNcdfAttr(ncid,attrName='dimUnlimName',&
            attr=trim(dimUnlimName1))
       print *,'finished dimUnlimName1'
       print *,'End writing dimension attributes'

       if (present(nchmw)) &
            call writeNcdfAttr(ncid,attrName='nchmw',attr=nchmw)

       if (present(nchan)) &
            call writeNcdfAttr(ncid,attrName='nchan',attr=nchan)

       if (present(freq)) &
            call writeNcdfAttr(ncid,attrName='mwfrequencies', &
            attr=freq)

       if (present(wvn)) &
            call writeNcdfAttr(ncid,attrName='wavenumbers', &
            attr=wvn)

       if (present(eia)) &
            call writeNcdfAttr(ncid,attrName='eia', &
            attr=eia)

!       if (present(eaa)) &
!            call writeNcdfAttr(ncid,attrName='eaa', &
!            attr=eaa)

       if (present(polarity)) &
            call writeNcdfAttr(ncid,attrName='mwpolarities',&
            attr=polarity)

       if (present(sto_x3d1_i1) .and. lputsto1(1)) &
            call writeNcdfAttr(ncid,attrName='sto_x3d1_i1',&
            attr=sto_x3d1_i1)

       if (present(sto_x3d1_i2) .and. lputsto1(2)) &
            call writeNcdfAttr(ncid,attrName='sto_x3d1_i2',&
            attr=sto_x3d1_i2)

       if (present(sto_x3d1_r4) .and. lputsto1(3)) &
            call writeNcdfAttr(ncid,attrName='sto_x3d1_r4',&
            attr=sto_x3d1_r4)

       if (present(sto_x2d1_i1) .and. lputsto1(4)) &
            call writeNcdfAttr(ncid,attrName='sto_x2d1_i1',&
            attr=sto_x2d1_i1)

       if (present(sto_x2d1_i2) .and. lputsto1(5)) &
            call writeNcdfAttr(ncid,attrName='sto_x2d1_i2',&
            attr=sto_x2d1_i2)

       if (present(sto_x2d1_r4) .and. lputsto1(6)) &
            call writeNcdfAttr(ncid,attrName='sto_x2d1_r4',&
            attr=sto_x2d1_r4)

       if (present(proj)) &
            call writeNcdfAttr(ncid,attrName='map_projection_type',attr=proj)

       if (present(reflat)) &
            call writeNcdfAttr(ncid,attrName='map_reference_latitude',attr=reflat)

       if (present(reflon)) &
            call writeNcdfAttr(ncid,attrName='map_reference_longitude',attr=reflon)

       if (present(latorigin)) &
            call writeNcdfAttr(ncid,attrName='map_origin_latitude',attr=latorigin)

       if (present(lonorigin)) &
            call writeNcdfAttr(ncid,attrName='map_origin_longitude',attr=lonorigin)

       if (present(rowGridOffset)) &
            call writeNcdfAttr(ncid,attrName='grid_origin_offset_row',attr=rowGridOffset)

       if (present(colGridOffset)) &
            call writeNcdfAttr(ncid,attrName='grid_origin_offset_col',attr=colGridOffset)

       if (present(scale)) &
            call writeNcdfAttr(ncid,attrName='map_scale',attr=scale)

       if (present(scaleunits)) &
            call writeNcdfAttr(ncid,attrName='map_scale_units',attr=scaleunits)

       if (present(re)) &
            call writeNcdfAttr(ncid,attrName='earth_radius',attr=re)

       if (present(itile)) &
            call writeNcdfAttr(ncid,attrName='tile_column_index',attr=itile)

       if (present(jtile)) &
            call writeNcdfAttr(ncid,attrName='tile_row_index',attr=jtile)

       if (present(nitile)) &
            call writeNcdfAttr(ncid,attrName='ncol_globaltiles',attr=nitile)

       if (present(njtile)) &
            call writeNcdfAttr(ncid,attrName='nrow_globaltiles',attr=njtile)

       if (present(xorigin)) &
            call writeNcdfAttr(ncid,attrName='xorigin',attr=xorigin)

       if (present(delx)) &
            call writeNcdfAttr(ncid,attrName='delta_x',attr=delx)

       if (present(xunits)) &
            call writeNcdfAttr(ncid,attrName='xunits', &
            attr=trim(xunits))

       if (present(yorigin)) &
            call writeNcdfAttr(ncid,attrName='yorigin',attr=yorigin)

       if (present(dely)) &
            call writeNcdfAttr(ncid,attrName='delta_y',attr=dely)

       if (present(yunits)) &
            call writeNcdfAttr(ncid,attrName='yunits', &
            attr=trim(yunits))

       if (present(timeLevelIncrement)) &
            call writeNcdfAttr(ncid,attrName='timeLevelIncrement',attr=timeLevelIncrement)

       if (present(timeLevelUnits)) &
            call writeNcdfAttr(ncid,attrName='timeLevelUnits',attr=timeLevelUnits)

    else
       lgetsto1(:)=.TRUE.
       if(present(lgetsto))lgetsto1=lgetsto

! open existing file in 'replace' mode, and return

       if(status == 'replace' .or. status == 'REPLACE')then
          call openNcdfFile(ncid, file, status='replace')
          return

! open existing file, and read global attributes

       else
!          print *,'open existing file for read-only access'
          call openNcdfFile(ncid, file, unlimited_dim_length=nUnlim1)

          if (ncid <= 0) then 
             print *, 'err[openSfcGrd]: Invalid file ID'
             return
          endif

!          print *,'nUnlim1=',nUnlim1
          if(present(nUnlim))then
             nUnlim=nUnlim1
!             print *,'nUnlim=',nUnlim
          endif
       endif

! check for consistency between dimension attributes and actual dimension lengths

! unlimited dim
       
!       print *,'Checking unlim dim consistency'
       call readNcdfAttr(ncid,attrName='nDimUnlim',attr=nDimUnlim1)
!       print *,'nDimUnlim1=',nDimUnlim1
       allocate(dimUnlimDims1(ndimUnlim1))
       call readNcdfAttr(ncid,attrName='dimUnlimDims',attr=dimUnlimDims1)
!       print *,'dimUnlimDims1=',dimUnlimDims1
       if(product(dimUnlimDims1) .ne. nUnlim1)then
          print *,'err[openSfcGrd]: unlimited_dim_length inconsistent with product(dimUnlimDims1)'
          print *,'nUnlim1,product(dimUnlimDims1)=',nUnlim1,product(dimUnlimDims1)
          stop
       endif
       allocate(dimNamesUnlim1(ndimUnlim1))
       do i=1,nDimUnlim1
          dimNamesUnlim1(i)=''
          write(sAttrNum,'(i1)')i
          call readNcdfAttr(ncid,attrName='dimNamesUnlim'//sAttrNum, &
               attr=dimNamesUnlim1(i))
       enddo

! fixed dims (currently limited to 1)

!       print *,'Checking fixed dim consistency'
       call readNcdfAttr(ncid,attrName='nDimFixed',attr=nDimFixed1)
!       print *,'nDimFixed1=',nDimFixed1
       if(nDimFixed1 .gt. 1)then
          stop 'err[openSfcGrd]: nDimFixed .gt. 1 not currently supported'
       endif
       allocate(dimFixedDims1(nDimFixed1),dimNamesFixed1(nDimFixed1))
       if(nDimFixed1 .gt. 0)call readNcdfAttr(ncid,attrName='dimFixedDims',attr=dimFixedDims1)
!       print *,'dimFixedDims1=',dimFixedDims1
       if(any(dimFixedDims1 .gt. 0))then
       do i=1,nDimFixed1
          dimNamesFixed1(i)=''
          write(sAttrNum,'(i1)')i
          call readNcdfAttr(ncid,attrName='dimNamesFixed'//sAttrNum, &
               attr=dimNamesFixed1(i))
!          print *,'i,dimNamesFixed1(i)=',i,trim(dimNamesFixed1(i))
          dimFixed1=readNcdfDim(ncid,trim(dimNamesFixed1(i)))
          if(dimFixed1 .ne. dimFixedDims1(i))then
             print *,'err[openSfcGrd]: actual file fixed dimension inconsistent with dimension attribute'
             print *,'Fixed dim index, dimNamesFixed, dimFixed1, dimFixedDims(i)=',&
                  i,trim(dimNamesFixed1(i)), dimFixed1, dimFixedDims1(i)
             stop
          endif
       enddo
       endif
       if (present(dimFixed)) dimFixed=dimFixed1

! unlimited dimname 

!       print *,'Checking unlim dimname and dim consistency'
       dimUnlimName1=''
       call readNcdfAttr(ncid,attrName='dimUnlimName',attr=dimUnlimName1)
!       print *,'dimUnlimName1=',trim(dimUnlimName1)
       dimUnlim1=readNcdfDim(ncid,trim(dimUnlimName1))
       if(dimUnlim1 .ne. nUnlim1)then
          print *,'err[openSfcGrd]: unlimited dimension inconsistency'
          print *,'dimUnlimName, dimUnlim1, nUnlim1=',trim(dimUnlimName1),dimUnlim1, nUnlim1
          stop
       endif
       if(present(dimUnlim)) dimUnlim=dimUnlim1

       if(allocated(dimUnlimDims1))deallocate(dimUnlimDims1)
       if(allocated(dimNamesUnlim1))deallocate(dimNamesUnlim1)
       if(allocated(dimFixedDims1))deallocate(dimFixedDims1)
       if(allocated(dimNamesFixed1))deallocate(dimNamesFixed1)

       if (present(creationDate)) &
            call readNcdfAttr(ncid,attrName='CreationTime', &
            attr=creationDate)

       if (present(nrow)) then
          nrow=0
          call readNcdfAttr(ncid,attrName='nDimUnlim',attr=nDimUnlim1)
          allocate(dimUnlimDims1(ndimUnlim1),dimNamesUnlim1(ndimUnlim1))
          call readNcdfAttr(ncid,attrName='dimUnlimDims',attr=dimUnlimDims1)
          do i=1,nDimUnlim1
             dimNamesUnlim1(i)=''
             write(sAttrNum,'(i1)')i
             call readNcdfAttr(ncid,attrName='dimNamesUnlim'//sAttrNum, &
                  attr=dimNamesUnlim1(i))
             if(dimNamesUnlim1(i)(1:4) .eq. 'nRow')nrow=dimUnlimDims1(i)
          enddo
          if(nrow .eq. 0)then
             print *,'err[openSfcGrd]: nrow attribute not found'
             stop
          endif
         deallocate(dimUnlimDims1,dimNamesUnlim1)
       endif

       if (present(ncol)) then
          ncol=0
          call readNcdfAttr(ncid,attrName='nDimUnlim',attr=nDimUnlim1)
          allocate(dimUnlimDims1(ndimUnlim1),dimNamesUnlim1(ndimUnlim1))
          call readNcdfAttr(ncid,attrName='dimUnlimDims',attr=dimUnlimDims1)
          do i=1,nDimUnlim1
             dimNamesUnlim1(i)=''
             write(sAttrNum,'(i1)')i
             call readNcdfAttr(ncid,attrName='dimNamesUnlim'//sAttrNum, &
                  attr=dimNamesUnlim1(i))
             if(dimNamesUnlim1(i)(1:4) .eq. 'nCol')ncol=dimUnlimDims1(i)
          enddo
          if(ncol .eq. 0)then
             print *,'err[openSfcGrd]: ncol attribute not found'
             stop
          endif
         deallocate(dimUnlimDims1,dimNamesUnlim1)
       endif

       if (present(nTimeLevels)) then
          nTimeLevels=0
          call readNcdfAttr(ncid,attrName='nDimUnlim',attr=nDimUnlim1)
          allocate(dimUnlimDims1(ndimUnlim1),dimNamesUnlim1(ndimUnlim1))
          call readNcdfAttr(ncid,attrName='dimUnlimDims',attr=dimUnlimDims1)
          do i=1,nDimUnlim1
             dimNamesUnlim1(i)=''
             write(sAttrNum,'(i1)')i
             call readNcdfAttr(ncid,attrName='dimNamesUnlim'//sAttrNum, &
                  attr=dimNamesUnlim1(i))
             if(dimNamesUnlim1(i)(1:11) .eq. 'nTimeLevels') &
                  nTimeLevels=dimUnlimDims1(i)
          enddo
          if(nTimeLevels .eq. 0)then
             print *,'err[openSfcGrd]: nTimeLevels attribute not found'
             stop
          endif
         deallocate(dimUnlimDims1,dimNamesUnlim1)
       endif

       if (present(xidlen)) xidlen=readNcdfDim(ncid,'nxId')
       if (present(nvalspergrid)) nvalspergrid=readNcdfDim(ncid,'nValsPerGrid')

!       if (present(dimUnlim)) then
!          dimUnlimName1=''
!          call readNcdfAttr(ncid,attrName='dimUnlimName',attr=dimUnlimName1)
!          print *,'dimUnlimName1=',trim(dimUnlimName1)
!          dimUnlim=readNcdfDim(ncid,trim(dimUnlimName1))
!       endif

       if (present(dimUnlimDims))then
          call readNcdfAttr(ncid,attrName='dimUnlimDims',attr=dimUnlimDims)
          call readNcdfAttr(ncid,attrName='nDimUnlim',attr=nDimUnlim)
          if (size(dimUnlimDims) /= nDimUnlim) &
               stop 'err[openSfcGrd]: # dimUnlimDims and nDimUnlim mismatched'
       endif
       if (present(nDimUnlim)) &
            call readNcdfAttr(ncid,attrName='nDimUnlim',attr=nDimUnlim)

       if (present(nDimFixed)) &
            call readNcdfAttr(ncid,attrName='nDimFixed',attr=nDimFixed)

       if (present(dimFixedDims))then
          call readNcdfAttr(ncid,attrName='nDimFixed',attr=nDimFixed)
          if(nDimFixed .gt. 0) then
             call readNcdfAttr(ncid,attrName='dimFixedDims',attr=dimFixedDims)
             if (size(dimFixedDims) /= nDimFixed) &
                  stop 'err[openSfcGrd]: # dimFixedDims and nDimFixed mismatched'
          else
             print *,'warning[openSfcGrd]: dimFixedDims present, but nDimFixed .eq. 0'
          endif
       endif

       if (present(dimNamesUnlim)) then
          call readNcdfAttr(ncid,attrName='nDimUnlim',attr=nDimUnlim1)
          do i=1,nDimUnlim1
             dimNamesUnlim(i)=''
             write(sAttrNum,'(i1)')i
             call readNcdfAttr(ncid,attrName='dimNamesUnlim'//sAttrNum, &
                  attr=dimNamesUnlim(i))
          enddo
       endif

       if (present(dimNamesFixed)) then
          call readNcdfAttr(ncid,attrName='nDimFixed',attr=nDimFixed1)
          do i=1,nDimFixed1
             dimNamesFixed(i)=''
             write(sAttrNum,'(i1)')i
             call readNcdfAttr(ncid,attrName='dimNamesFixed'//sAttrNum, &
                  attr=dimNamesFixed(i))
          enddo
       endif

       if (present(dimUnlimName)) then
          dimUnlimName1=''
          call readNcdfAttr(ncid,attrName='dimUnlimName',attr=dimUnlimName)
       endif


       if (present(freq)) then
          call readNcdfAttr(ncid,attrName='mwfrequencies',attr=freq)
          call readNcdfAttr(ncid,attrName='nchmw',attr=nchanloc,silent=.true.)
          if (nchanloc == 0) &
             call readNcdfAttr(ncid,attrName='nchan',attr=nchanloc)
          if (size(freq) /= nchanloc) &
               stop 'err[openSfcGrd]: #freq and nchan mismatched'
       endif

       if (present(wvn)) then
          call readNcdfAttr(ncid,attrName='wavenumbers',attr=wvn)
          call readNcdfAttr(ncid,attrName='nchan',attr=nchanloc)
          if (size(wvn) /= nchanloc) &
               stop 'err[openSfcGrd]: #wvn and nchan mismatched'
       endif

       if (present(eia)) then
          call readNcdfAttr(ncid,attrName='eia',attr=eia)
          call readNcdfAttr(ncid,attrName='nchmw',attr=nchanloc,silent=.true.)
          if (nchanloc == 0) &
             call readNcdfAttr(ncid,attrName='nchan',attr=nchanloc)
          if (size(eia) /= nchanloc) &
               stop 'err[openSfcGrd]: #eia and nchan mismatched'
       endif

       if (present(eaa)) then
          call readNcdfAttr(ncid,attrName='eaa',attr=eaa)
          call readNcdfAttr(ncid,attrName='nchmw',attr=nchanloc,silent=.true.)
          if (nchanloc == 0) &
             call readNcdfAttr(ncid,attrName='nchan',attr=nchanloc)
          if (size(eaa) /= nchanloc) &
               stop 'err[openSfcGrd]: #eaa and nchan mismatched'
       endif

       if (present(polarity)) then
          call readNcdfAttr(ncid,attrName='mwpolarities', &
               attr=polarity)
          call readNcdfAttr(ncid,attrName='nchmw',attr=nchanloc,silent=.true.)
          if (nchanloc == 0) &
             call readNcdfAttr(ncid,attrName='nchan',attr=nchanloc)
          if (size(polarity) /= nchanloc) &
               stop 'err[openSfcGrd]: #mwpol and nchan mismatched'
       endif

       if (present(nchmw)) &
               call readNcdfAttr(ncid,attrName='nchmw',attr=nchmw)

       if (present(nchan)) then
           call readNcdfAttr(ncid,attrName='nchan',attr=nchan,silent=.true.)
           if (nchan == 0) &  ! for back compatibility
             call readNcdfAttr(ncid,attrName='nchmw',attr=nchan,silent=.true.)
       endif

       if (present(casename)) &
            call readNcdfAttr(ncid,attrName='case', &
            attr=casename)

       if (present(proj)) &
            call readNcdfAttr(ncid,attrName='map_projection_type',attr=proj)

       if (present(reflat)) &
            call readNcdfAttr(ncid,attrName='map_reference_latitude',attr=reflat)

       if (present(reflon)) &
            call readNcdfAttr(ncid,attrName='map_reference_longitude',attr=reflon)

       if (present(latorigin)) &
            call readNcdfAttr(ncid,attrName='map_origin_latitude',attr=latorigin)

       if (present(lonorigin)) &
            call readNcdfAttr(ncid,attrName='map_origin_longitude',attr=lonorigin)

       if (present(rowGridOffset)) &
            call readNcdfAttr(ncid,attrName='grid_origin_offset_row',attr=rowGridOffset)

       if (present(colGridOffset)) &
            call readNcdfAttr(ncid,attrName='grid_origin_offset_col',attr=colGridOffset)

       if (present(scale)) &
            call readNcdfAttr(ncid,attrName='map_scale',attr=scale)

       if (present(scaleunits)) &
            call readNcdfAttr(ncid,attrName='map_scale_units', &
            attr=scaleunits)

       if (present(re)) &
            call readNcdfAttr(ncid,attrName='earth_radius',attr=re)

       if (present(itile)) &
            call readNcdfAttr(ncid,attrName='tile_column_index',attr=itile)

       if (present(jtile)) &
            call readNcdfAttr(ncid,attrName='tile_row_index',attr=jtile)

       if (present(nitile)) &
            call readNcdfAttr(ncid,attrName='ncol_globaltiles',attr=nitile)

       if (present(njtile)) &
            call readNcdfAttr(ncid,attrName='nrow_globaltiles',attr=njtile)

       if (present(xorigin)) &
            call readNcdfAttr(ncid,attrName='xorigin',attr=xorigin)

       if (present(delx)) &
            call readNcdfAttr(ncid,attrName='delta_x',attr=delx)

       if (present(xunits)) &
            call readNcdfAttr(ncid,attrName='xunits', &
            attr=xunits)

       if (present(yorigin)) &
            call readNcdfAttr(ncid,attrName='yorigin',attr=yorigin)

       if (present(dely)) &
            call readNcdfAttr(ncid,attrName='delta_y',attr=dely)

       if (present(yunits)) &
            call readNcdfAttr(ncid,attrName='yunits', &
            attr=yunits)

       if (present(timeLevelIncrement)) &
            call readNcdfAttr(ncid,attrName='timeLevelIncrement',attr=timeLevelIncrement)

       if (present(timeLevelUnits)) &
            call readNcdfAttr(ncid,attrName='timeLevelUnits', &
            attr=timeLevelUnits)

       if (present(sto_x3d1_i1) .and. lgetsto1(1)) &
            call readNcdfAttr(ncid,attrName='sto_x3d1_i1',attr=sto_x3d1_i1)

       if (present(sto_x3d1_i2) .and. lgetsto1(2)) &
            call readNcdfAttr(ncid,attrName='sto_x3d1_i2',attr=sto_x3d1_i2)

       if (present(sto_x3d1_r4) .and. lgetsto1(3)) &
            call readNcdfAttr(ncid,attrName='sto_x3d1_r4',attr=sto_x3d1_r4)

       if (present(sto_x2d1_i1) .and. lgetsto1(4)) &
            call readNcdfAttr(ncid,attrName='sto_x2d1_i1',attr=sto_x2d1_i1)

       if (present(sto_x2d1_i2) .and. lgetsto1(5)) &
            call readNcdfAttr(ncid,attrName='sto_x2d1_i2',attr=sto_x2d1_i2)

       if (present(sto_x2d1_r4) .and. lgetsto1(6)) &
            call readNcdfAttr(ncid,attrName='sto_x2d1_r4',attr=sto_x2d1_r4)

    endif

  end subroutine openNcSfcGrd

  ! -
  subroutine queryNcSfcGrd(file, creationDate, ncol, nrow, nValsPerGrid, nTimeLevels, &
       nUnlim, dimUnlimName, &
       ndimUnlim, dimUnlimDims, dimNamesUnlim, &
       ndimFixed, dimFixedDims, dimNamesFixed, &       
       dimUnlim, dimFixed, &
       polarity, freq, wvn, eia, eaa, casename, nchmw, nchan, reflat, reflon, &
       proj, latorigin, lonorigin, rowGridOffset, colGridOffset, scale, scaleunits, re, &
       itile, jtile, nitile, njtile, xorigin, delx, xunits, yorigin, dely, yunits, &
       lgetsto, &
       sto_x3d1_i1, sto_x3d1_i2, sto_x3d1_r4, &
       sto_x2d1_i1, sto_x2d1_i2, sto_x2d1_r4, &
       timeLevelIncrement, timeLevelUnits, &
       xidlen)

    ! - Dummy variables
    character (len=*), intent (in) :: file
    character (len=*), intent (inout), optional :: creationDate
    integer, intent(inout), optional :: xidlen
    integer, intent(inout), optional :: nrow, ncol, nValsPerGrid, nTimeLevels, nUnlim
    real, intent(inout), optional :: reflat, reflon, latorigin, lonorigin
    real, intent(inout), optional :: xorigin, delx, yorigin, dely
    real, intent(inout), optional :: rowGridOffset, colGridOffset
    integer, intent(inout), optional :: itile, jtile, nitile, njtile
    integer, intent(inout), optional :: nchmw
    integer, intent(inout), optional :: nchan
    real, intent (inout), optional :: scale, re
    real, dimension(:), intent (inout), optional :: freq, eia, eaa
    real, dimension(:), intent (inout), optional :: wvn
    integer, dimension(:), intent (inout), optional :: polarity
    logical, dimension(:), intent(in), optional :: lgetsto
    real, dimension(2), intent (inout), optional :: sto_x3d1_i1, sto_x3d1_i2, sto_x3d1_r4
    real, dimension(2), intent (inout), optional :: sto_x2d1_i1, sto_x2d1_i2, sto_x2d1_r4
    integer, intent(inout), optional :: ndimUnlim, ndimFixed, dimUnlim, dimFixed
    integer, dimension(:), intent(inout), optional :: dimUnlimDims, dimFixedDims
    real, intent(inout), optional :: timeLevelIncrement
    character (len=*), intent(inout), optional :: casename
    character (len=*), intent(inout), optional :: xunits, yunits, scaleunits
    character (len=*), intent(inout), optional :: proj
    character (len=*), dimension(:), intent(inout), optional :: dimNamesUnlim
    character (len=*), dimension(:), intent(inout), optional :: dimNamesFixed
    character (len=*), intent(inout), optional :: dimUnlimName
    character (len=*), intent(inout), optional :: timeLevelUnits
    ! NOTE: add other dummy arguments here as needed
    ! to support other global attributes

! - Local variables
    integer :: ncid
    logical, dimension(6) :: lgetsto1
    lgetsto1=.TRUE.
    if(present(lgetsto))lgetsto1=lgetsto
    print *,'inside querySfcGrd: call openNcSfcGrd'
      
    call openNcSfcGrd(ncid=ncid, file=file, creationDate=creationDate, &
         ncol=ncol, nrow=nrow, nValsPerGrid=nValsPerGrid, nTimeLevels=nTimeLevels, &
         nUnlim=nUnlim, dimUnlimName=dimUnlimName, &
         ndimUnlim=ndimUnlim, dimUnlimDims=dimUnlimDims, dimNamesUnlim=dimNamesUnlim, &
         ndimFixed=ndimFixed, dimFixedDims=dimFixedDims, dimNamesFixed=dimNamesFixed, &       
         dimUnlim=dimUnlim, dimFixed=dimFixed, &
         polarity=polarity, freq=freq, wvn=wvn, eia=eia, eaa=eaa, casename=casename, &
         nchmw=nchmw, nchan=nchan, proj=proj, &
         reflat=reflat, reflon=reflon, latorigin=latorigin, lonorigin=lonorigin, &
         rowGridOffset=rowGridOffset, colGridOffset=colGridOffset, scale=scale, scaleunits=scaleunits, re=re, &
         itile=itile, jtile=jtile, nitile=nitile, njtile=njtile, &
         xorigin=xorigin, delx=delx, xunits=xunits, yorigin=yorigin, dely=dely, yunits=yunits, &
         lgetsto=lgetsto1, &
         sto_x3d1_i1=sto_x3d1_i1, sto_x3d1_i2=sto_x3d1_i2, sto_x3d1_r4=sto_x3d1_r4, &
         sto_x2d1_i1=sto_x2d1_i1, sto_x2d1_i2=sto_x2d1_i2, sto_x2d1_r4=sto_x2d1_r4, &
         xidlen=xidlen, status='old')
    call closeNcGFile(ncid)
  end subroutine queryNcSfcGrd

! -

! Cover routine for getNcSfcGrd
! Used for arrays passed in with "original" rank 
! (i.e. 4 for (nValsPerGrid,nCol,nRow,nTime) arrays, 
!       3 for (nCol,nRow,nTime) arrays)

  subroutine get4DSfcGrd(ncid, imask, iloc, jloc, kloc, &
       xid, &
       x3d1_i1, x3d1_i2, x3d1_r4, &
       x2d1_i1, x2d1_i2, x2d1_r4, &
       year, &
       julday, &
       year1d, &
       julday1d)
    integer, intent(in) :: ncid
    integer, dimension(:,:,:,:), intent(in) :: imask
    integer, intent(in), optional :: iloc, jloc, kloc
    character (len=*), dimension(:), intent(inout), optional :: xid
    integer*1, dimension(:,:,:,:), intent (inout), optional :: x3d1_i1
    integer*2, dimension(:,:,:,:), intent (inout), optional :: x3d1_i2
    real, dimension(:,:,:,:), intent (inout), optional :: x3d1_r4
    integer*1, dimension(:,:,:), intent (inout), optional :: x2d1_i1
    integer*2, dimension(:,:,:), intent (inout), optional :: x2d1_i2
    real, dimension(:,:,:), intent (inout), optional :: x2d1_r4
    integer*2, dimension(:,:,:), intent (inout), optional :: year
    real, dimension(:,:,:), intent (inout), optional :: julday
    integer*2, dimension(:), intent (inout), optional :: year1d
    real, dimension(:), intent (inout), optional :: julday1d

    ! NOTE: add other dummy arguments here as needed
    ! to support other gridded fields (e.g. Tskin, soil moisture, etc.)

! Local variables

! Data arrays used for interface with getNcSfcGrd
    
    integer*1, allocatable, dimension(:) :: x3d1vec_i1
    integer*1, allocatable, dimension(:) :: x2d1vec_i1
    integer*2, allocatable, dimension(:) :: x3d1vec_i2
    integer*2, allocatable, dimension(:) :: x2d1vec_i2
    real, allocatable, dimension(:) :: x3d1vec_r4
    real, allocatable, dimension(:) :: x2d1vec_r4
    integer*2, allocatable, dimension(:) :: yearvec
    real, allocatable, dimension(:) :: juldayvec
    integer*2, allocatable, dimension(:) :: year1dvec
    real, allocatable, dimension(:) :: julday1dvec
    character (len=8), allocatable, dimension(:) :: xidvec


!    integer*1, allocatable, dimension(:,:,:,:) :: x3d1_i1_1
!    integer*2, allocatable, dimension(:,:,:,:) :: x3d1_i2_1
!    real, allocatable, dimension(:,:,:,:) :: x3d1_r4_1
!    integer*1, allocatable, dimension(:,:,:) :: x2d1_i1_1
!    integer*2, allocatable, dimension(:,:,:) :: x2d1_i2_1
!    real, allocatable, dimension(:,:,:) :: x2d1_r4_1
!    integer*2, allocatable, dimension(:,:,:) :: year_1
!    real, allocatable, dimension(:,:,:) :: julday_1

!
    logical :: lgetid1, lgetyear1, lgetjulday1, lgetyear1d1, lgetjulday1d1
    logical, dimension(6) :: lgetdata1

    integer :: iloc1, jloc1, kloc1
    integer :: nElem, ni, nj, nk, nFixed

!
    print *,'get4DSfcGrd:'
!
    nElem=0
    ni=0
    nj=0
    nk=0
    nFixed=0
    
! Allocate 1d arrays and initialize to missing 
! (needed as placeholders in call to getNcSfcGrd if not actually present)

    allocate(x3d1vec_i1(1))
    allocate(x3d1vec_i2(1))
    allocate(x3d1vec_r4(1))
    allocate(x2d1vec_i1(1))
    allocate(x2d1vec_i2(1))
    allocate(x2d1vec_r4(1))
    allocate(yearvec(1))
    allocate(juldayvec(1))
    allocate(year1dvec(1))
    allocate(julday1dvec(1))
    allocate(xidvec(1))

! Initialize local variables iloc1,jloc1,kloc1 and reset if corresponding optional argument
! is present in argument list

    iloc1=0
    jloc1=0
    kloc1=0
    if(present(iloc))iloc1=iloc
    if(present(jloc))jloc1=jloc
    if(present(kloc))kloc1=kloc

! Set defaults for logical variables to FALSE 
! and then reset to TRUE if corresponding optional argument is present

    lgetid1=.FALSE.
    lgetyear1=.FALSE.
    lgetjulday1=.FALSE.
    lgetdata1=.FALSE.
    lgetyear1d1=.FALSE.
    lgetjulday1d1=.FALSE.

!    print *,'getSfcGrd: lgetdata1=',lgetdata1


    if (present(xid)) then
       lgetid1=.TRUE.
       if(allocated(xidvec))deallocate(xidvec)
       nElem=size(xid)
       allocate(xidvec(nElem))
    endif
    if (present(year1d)) then
       lgetyear1d1=.TRUE.
       if(allocated(year1dvec))deallocate(year1dvec)
       nElem=size(year1d)
       allocate(year1dvec(nElem))
    endif
    if (present(julday1d)) then
       lgetjulday1d1=.TRUE.
       if(allocated(julday1dvec))deallocate(julday1dvec)
       nElem=size(julday1d)
       allocate(julday1dvec(nElem))
    endif

    if (present(year)) then
       lgetyear1=.TRUE.
       if(allocated(yearvec))deallocate(yearvec)
       nElem=size(year)
       allocate(yearvec(nElem))
    endif
    if (present(julday)) then
       lgetjulday1=.TRUE.
       if(allocated(juldayvec))deallocate(juldayvec)
       nElem=size(julday)
       allocate(juldayvec(nElem))
    endif


    if (present(x3d1_i1)) then
       lgetdata1(1)=.TRUE.
       if(allocated(x3d1vec_i1))deallocate(x3d1vec_i1)
       nElem=size(x3d1_i1)
       allocate(x3d1vec_i1(nElem))
    endif

    if (present(x3d1_i2)) then
       lgetdata1(2)=.TRUE.
       if(allocated(x3d1vec_i2))deallocate(x3d1vec_i2)
       nElem=size(x3d1_i2)
       allocate(x3d1vec_i2(nElem))
    endif

    if (present(x3d1_r4)) then
       lgetdata1(3)=.TRUE.
       if(allocated(x3d1vec_r4))deallocate(x3d1vec_r4)
       nElem=size(x3d1_r4)
       allocate(x3d1vec_r4(nElem))
    endif

    if (present(x2d1_i1)) then
       lgetdata1(4)=.TRUE.
       if(allocated(x2d1vec_i1))deallocate(x2d1vec_i1)
       nElem=size(x2d1_i1)
       allocate(x2d1vec_i1(nElem))
    endif

    if (present(x2d1_i2)) then
       lgetdata1(5)=.TRUE.
       if(allocated(x2d1vec_i2))deallocate(x2d1vec_i2)
       nElem=size(x2d1_i2)
       allocate(x2d1vec_i2(nElem))
    endif

    if (present(x2d1_r4)) then
       lgetdata1(6)=.TRUE.
       if(allocated(x2d1vec_r4))deallocate(x2d1vec_r4)
       nElem=size(x2d1_r4)
       allocate(x2d1vec_r4(nElem))
    endif

    call getNcSfcGrd(ncid=ncid, iloc=iloc1, jloc=jloc1, kloc=kloc1, &
         lgetid=lgetid1, xid=xidvec, &
         lgetdata=lgetdata1, &
         x3d1_i1=x3d1vec_i1, &
         x3d1_i2=x3d1vec_i2, &
         x3d1_r4=x3d1vec_r4, &
         x2d1_i1=x2d1vec_i1, &
         x2d1_i2=x2d1vec_i2, &
         x2d1_r4=x2d1vec_r4, &
         lgetyear=lgetyear1, year=yearvec, &
         lgetjulday=lgetjulday1, julday=juldayvec, &
         lgetyear1d=lgetyear1d1, year1d=year1dvec, &
         lgetjulday1d=lgetjulday1d1, julday1d=julday1dvec)
    
    if(present(x3d1_i1))then
       nFixed=size(x3d1_i1,1)
       ni=size(x3d1_i1,2)
       nj=size(x3d1_i1,3)
       nk=size(x3d1_i1,4)
       x3d1_i1=reshape(x3d1vec_i1,(/nFixed,ni,nj,nk/))
    endif

    if(present(x3d1_i2))then
       nFixed=size(x3d1_i2,1)
       ni=size(x3d1_i2,2)
       nj=size(x3d1_i2,3)
       nk=size(x3d1_i2,4)
       print *,'x3d1_i2 present: nFixed,ni,nj,nk=',nFixed,ni,nj,nk
       print *,'size(x3d1_i2)=',size(x3d1_i2)
       print *,'size(x3d1vec_i2)=',size(x3d1vec_i2)
       x3d1_i2=reshape(x3d1vec_i2,(/nFixed,ni,nj,nk/))
    endif

    if(present(x3d1_r4))then
       nFixed=size(x3d1_r4,1)
       ni=size(x3d1_r4,2)
       nj=size(x3d1_r4,3)
       nk=size(x3d1_r4,4)
       print *,'x3d1_r4 present: nFixed,ni,nj,nk=',nFixed,ni,nj,nk
       x3d1_r4=reshape(x3d1vec_r4,(/nFixed,ni,nj,nk/))
    endif

    if(present(x2d1_i1))then
       ni=size(x2d1_i1,1)
       nj=size(x2d1_i1,2)
       nk=size(x2d1_i1,3)
       x2d1_i1=reshape(x2d1vec_i1,(/ni,nj,nk/))
    endif

    if(present(x2d1_i2))then
       ni=size(x2d1_i2,1)
       nj=size(x2d1_i2,2)
       nk=size(x2d1_i2,3)
       print *,'x2d1_i2 present: ni,nj,nk=',ni,nj,nk
       print *,'size(x2d1_i2)=',size(x2d1_i2)
       print *,'size(x2d1vec_i2)=',size(x2d1vec_i2)
       x2d1_i2=reshape(x2d1vec_i2,(/ni,nj,nk/))
    endif

    if(present(x2d1_r4))then
       ni=size(x2d1_r4,1)
       nj=size(x2d1_r4,2)
       nk=size(x2d1_r4,3)
       x2d1_r4=reshape(x2d1vec_r4,(/ni,nj,nk/))
    endif

    if(present(year))then
       print *,'year present'
       ni=size(year,1)
       nj=size(year,2)
       nk=size(year,3)
       year=reshape(yearvec,(/ni,nj,nk/))
    endif

    if(present(julday))then
       print *,'julday present'
       ni=size(julday,1)
       nj=size(julday,2)
       nk=size(julday,3)
       julday=reshape(juldayvec,(/ni,nj,nk/))
    endif

    if(present(year1d))then
       print *,'year1d present'
       year1d=year1dvec
    endif

    if(present(julday1d))then
       print *,'julday1d present'
       julday1d=julday1dvec
    endif

    if(present(xid))then
       print *,'xid present'
       xid=xidvec
    endif

    deallocate(x3d1vec_i1, x3d1vec_i2, x3d1vec_r4)
    deallocate(x2d1vec_i1, x2d1vec_i2, x2d1vec_r4)
    deallocate(yearvec, juldayvec)
    deallocate(year1dvec, julday1dvec)
    deallocate(xidvec)

    return


  end subroutine get4DSfcGrd

! --- 

! --- 
! Cover routine for getNcSfcGrd
! Used for arrays passed in with reduced rather than "original" rank
! using reshape intrinsic function
! (i.e. 1 for (nValsPerGrid,nCol,nRow,nTime) arrays, 
!       1 for (nCol,nRow,nTime) arrays)


  subroutine get1DSfcGrd(ncid, imask, iloc, jloc, kloc, &
       xid, &
       x3d1_i1, x3d1_i2, x3d1_r4, &
       x2d1_i1, x2d1_i2, x2d1_r4, &
       year, &
       julday, &
       year1d, &
       julday1d)
    integer, intent(in) :: ncid
    integer, dimension(:), intent(in) :: imask
    integer, intent(in), optional :: iloc, jloc, kloc
    character (len=*), dimension(:), intent(inout), optional :: xid
    integer*1, dimension(:), intent (inout), optional :: x3d1_i1
    integer*2, dimension(:), intent (inout), optional :: x3d1_i2
    real, dimension(:), intent (inout), optional :: x3d1_r4
    integer*1, dimension(:), intent (inout), optional :: x2d1_i1
    integer*2, dimension(:), intent (inout), optional :: x2d1_i2
    real, dimension(:), intent (inout), optional :: x2d1_r4
    integer*2, dimension(:), intent (inout), optional :: year
    real, dimension(:), intent (inout), optional :: julday
    integer*2, dimension(:), intent (inout), optional :: year1d
    real, dimension(:), intent (inout), optional :: julday1d

    ! NOTE: add other dummy arguments here as needed
    ! to support other gridded fields (e.g. Tskin, soil moisture, etc.)

! Local variables

! Data arrays used for interface with getNcSfcGrd
    
    integer*1, allocatable, dimension(:) :: x3d1vec_i1
    integer*1, allocatable, dimension(:) :: x2d1vec_i1
    integer*2, allocatable, dimension(:) :: x3d1vec_i2
    integer*2, allocatable, dimension(:) :: x2d1vec_i2
    real, allocatable, dimension(:) :: x3d1vec_r4
    real, allocatable, dimension(:) :: x2d1vec_r4
    integer*2, allocatable, dimension(:) :: yearvec
    real, allocatable, dimension(:) :: juldayvec
    integer*2, allocatable, dimension(:) :: year1dvec
    real, allocatable, dimension(:) :: julday1dvec
    character (len=8), allocatable, dimension(:) :: xidvec


!    integer*1, allocatable, dimension(:,:,:,:) :: x3d1_i1_1
!    integer*2, allocatable, dimension(:,:,:,:) :: x3d1_i2_1
!    real, allocatable, dimension(:,:,:,:) :: x3d1_r4_1
!    integer*1, allocatable, dimension(:,:,:) :: x2d1_i1_1
!    integer*2, allocatable, dimension(:,:,:) :: x2d1_i2_1
!    real, allocatable, dimension(:,:,:) :: x2d1_r4_1
!    integer*2, allocatable, dimension(:,:,:) :: year_1
!    real, allocatable, dimension(:,:,:) :: julday_1

!
    logical :: lgetid1, lgetyear1, lgetjulday1, lgetyear1d1, lgetjulday1d1
    logical, dimension(6) :: lgetdata1

    integer :: iloc1, jloc1, kloc1
    integer :: nElem, ni, nj, nk, nFixed

!
    nElem=0
    ni=0
    nj=0
    nk=0
    nFixed=0
    
! Allocate 1d arrays and initialize to missing 
! (needed as placeholders in call to getNcSfcGrd if not actually present)

    allocate(x3d1vec_i1(1))
    allocate(x3d1vec_i2(1))
    allocate(x3d1vec_r4(1))
    allocate(x2d1vec_i1(1))
    allocate(x2d1vec_i2(1))
    allocate(x2d1vec_r4(1))
    allocate(yearvec(1))
    allocate(juldayvec(1))
    allocate(year1dvec(1))
    allocate(julday1dvec(1))
    allocate(xidvec(1))

! Initialize local variables iloc1,jloc1,kloc1 and reset if corresponding optional argument
! is present in argument list

    iloc1=0
    jloc1=0
    kloc1=0
    if(present(iloc))iloc1=iloc
    if(present(jloc))jloc1=jloc
    if(present(kloc))kloc1=kloc

! Set defaults for logical variables to FALSE 
! and then reset to TRUE if corresponding optional argument is present

    lgetid1=.FALSE.
    lgetyear1=.FALSE.
    lgetjulday1=.FALSE.
    lgetdata1=.FALSE.
    lgetyear1d1=.FALSE.
    lgetjulday1d1=.FALSE.

!    print *,'getSfcGrd: lgetdata1=',lgetdata1


    if (present(xid)) then
       lgetid1=.TRUE.
       if(allocated(xidvec))deallocate(xidvec)
       nElem=size(xid)
       allocate(xidvec(nElem))
    endif
    if (present(year1d)) then
       lgetyear1d1=.TRUE.
       if(allocated(year1dvec))deallocate(year1dvec)
       nElem=size(year1d)
       allocate(year1dvec(nElem))
    endif
    if (present(julday1d)) then
       lgetjulday1d1=.TRUE.
       if(allocated(julday1dvec))deallocate(julday1dvec)
       nElem=size(julday1d)
       allocate(julday1dvec(nElem))
    endif

    if (present(year)) then
       lgetyear1=.TRUE.
       if(allocated(yearvec))deallocate(yearvec)
       nElem=size(year)
       allocate(yearvec(nElem))
    endif
    if (present(julday)) then
       lgetjulday1=.TRUE.
       if(allocated(juldayvec))deallocate(juldayvec)
       nElem=size(julday)
       allocate(juldayvec(nElem))
    endif

    if (present(x3d1_i1)) then
       lgetdata1(1)=.TRUE.
       if(allocated(x3d1vec_i1))deallocate(x3d1vec_i1)
       nElem=size(x3d1_i1)
       allocate(x3d1vec_i1(nElem))
    endif

    if (present(x3d1_i2)) then
       lgetdata1(2)=.TRUE.
       if(allocated(x3d1vec_i2))deallocate(x3d1vec_i2)
       nElem=size(x3d1_i2)
       allocate(x3d1vec_i2(nElem))
    endif

    if (present(x3d1_r4)) then
       lgetdata1(3)=.TRUE.
       if(allocated(x3d1vec_r4))deallocate(x3d1vec_r4)
       nElem=size(x3d1_r4)
       allocate(x3d1vec_r4(nElem))
    endif

    if (present(x2d1_i1)) then
       lgetdata1(4)=.TRUE.
       if(allocated(x2d1vec_i1))deallocate(x2d1vec_i1)
       nElem=size(x2d1_i1)
       allocate(x2d1vec_i1(nElem))
    endif

    if (present(x2d1_i2)) then
       lgetdata1(5)=.TRUE.
       if(allocated(x2d1vec_i2))deallocate(x2d1vec_i2)
       nElem=size(x2d1_i2)
       allocate(x2d1vec_i2(nElem))
    endif

    if (present(x2d1_r4)) then
       lgetdata1(6)=.TRUE.
       if(allocated(x2d1vec_r4))deallocate(x2d1vec_r4)
       nElem=size(x2d1_r4)
       allocate(x2d1vec_r4(nElem))
    endif

    call getNcSfcGrd(ncid=ncid, iloc=iloc1, jloc=jloc1, kloc=kloc1, &
         lgetid=lgetid1, xid=xidvec, &
         lgetdata=lgetdata1, &
         x3d1_i1=x3d1vec_i1, &
         x3d1_i2=x3d1vec_i2, &
         x3d1_r4=x3d1vec_r4, &
         x2d1_i1=x2d1vec_i1, &
         x2d1_i2=x2d1vec_i2, &
         x2d1_r4=x2d1vec_r4, &
         lgetyear=lgetyear1, year=yearvec, &
         lgetjulday=lgetjulday1, julday=juldayvec, &
         lgetyear1d=lgetyear1d1, year1d=year1dvec, &
         lgetjulday1d=lgetjulday1d1, julday1d=julday1dvec)
    
    if(present(x3d1_i1))then
       x3d1_i1=x3d1vec_i1
    endif

    if(present(x3d1_i2))then
       x3d1_i2=x3d1vec_i2
    endif

    if(present(x3d1_r4))then
       x3d1_r4=x3d1vec_r4
    endif

    if(present(x2d1_i1))then
       x2d1_i1=x2d1vec_i1
    endif

    if(present(x2d1_i2))then
       x2d1_i2=x2d1vec_i2
    endif

    if(present(x2d1_r4))then
       x2d1_r4=x2d1vec_r4
    endif

    if(present(year))then
       year=yearvec
    endif

    if(present(julday))then
       julday=juldayvec
    endif

    if(present(year1d))then
       print *,'year1d present'
       year1d=year1dvec
    endif

    if(present(julday1d))then
       print *,'julday1d present'
       julday1d=julday1dvec
    endif

    if(present(xid))then
       print *,'xid present'
       xid=xidvec
    endif

    deallocate(x3d1vec_i1, x3d1vec_i2, x3d1vec_r4)
    deallocate(x2d1vec_i1, x2d1vec_i2, x2d1vec_r4)
    deallocate(yearvec, juldayvec)
    deallocate(year1dvec, julday1dvec)
    deallocate(xidvec)

    return


  end subroutine get1DSfcGrd



! --- 

! -
  subroutine getNcSfcGrd(ncid, iloc, jloc, kloc, &
       lgetid, xid, &
       lgetdata, &
       x3d1_i1, x3d1_i2, x3d1_r4, &
       x2d1_i1, x2d1_i2, x2d1_r4, &
       lgetyear, year, &
       lgetjulday, julday, &
       lgetyear1d, year1d, &
       lgetjulday1d, julday1d)
    integer, intent(in) :: ncid
    integer, intent(in), optional :: iloc, jloc, kloc
    logical, intent(in), optional :: lgetid
    logical, intent(in), optional :: lgetyear
    logical, intent(in), optional :: lgetjulday
    logical, intent(in), optional :: lgetyear1d
    logical, intent(in), optional :: lgetjulday1d
    logical, dimension(:), intent(in), optional :: lgetdata
    character (len=*), dimension(:), intent(inout), optional :: xid
    integer*1, dimension(:), intent (inout), optional :: x3d1_i1
    integer*2, dimension(:), intent (inout), optional :: x3d1_i2
    real, dimension(:), intent (inout), optional :: x3d1_r4
    integer*1, dimension(:), intent (inout), optional :: x2d1_i1
    integer*2, dimension(:), intent (inout), optional :: x2d1_i2
    real, dimension(:), intent (inout), optional :: x2d1_r4
    integer*2, dimension(:), intent (inout), optional :: year
    real, dimension(:), intent (inout), optional :: julday
    integer*2, dimension(:), intent (inout), optional :: year1d
    real, dimension(:), intent (inout), optional :: julday1d

    ! NOTE: add other dummy arguments here as needed
    ! to support other gridded fields (e.g. Tskin, soil moisture, etc.)

    integer :: ni, nj, nk, nelemFixed, nelemUnlim, i, irec
    integer :: ndimUnlim1, ndimFixed1
    integer, allocatable, dimension(:) :: dimUnlimDims1, dimFixedDims1

! Local data arrays used for input from nCDF file
    
    integer*1, allocatable, dimension(:,:) :: x3d1vec_i1
    integer*1, allocatable, dimension(:) :: x2d1vec_i1
    integer*2, allocatable, dimension(:,:) :: x3d1vec_i2
    integer*2, allocatable, dimension(:) :: x2d1vec_i2
    real, allocatable, dimension(:,:) :: x3d1vec_r4
    real, allocatable, dimension(:) :: x2d1vec_r4
    integer*2, allocatable, dimension(:) :: yearvec
    real, allocatable, dimension(:) :: juldayvec

    integer*2 x2d1scalar_i2

!
    logical :: lgetid1, lgetyear1, lgetjulday1, lgetyear1d1, lgetjulday1d1
    logical, dimension(6) :: lgetdata1

    lgetid1=.TRUE.
    lgetyear1=.TRUE.
    lgetjulday1=.TRUE.
    lgetdata1=.TRUE.
    lgetyear1d1=.TRUE.
    lgetjulday1d1=.TRUE.

    if(present(lgetid))lgetid1=lgetid
    if(present(lgetyear))lgetyear1=lgetyear
    if(present(lgetjulday))lgetjulday1=lgetjulday
    if(present(lgetdata))lgetdata1=lgetdata
    if(present(lgetyear))lgetyear1d1=lgetyear1d
    if(present(lgetjulday))lgetjulday1d1=lgetjulday1d
!    print *,'getSfcGrd: lgetdata1=',lgetdata1

! Before reading any data from the file,
! we need to determine the dimension attributes of the data in the file.
! Then we can compare these with the dimensions of the input arguments
! to check for consistency and to determine if we are reading a single record
! or the entire data array stored in the file.

!    print *,'getSfcGrd: Reading dimension attributes'
    call readNcdfAttr(ncid,attrName='nDimUnlim',attr=nDimUnlim1)
!    print *,'nDimUnlim1=',nDimUnlim1
    allocate(dimUnlimDims1(nDimUnlim1))
    call readNcdfAttr(ncid,attrName='dimUnlimDims',attr=dimUnlimDims1)

    call readNcdfAttr(ncid,attrName='nDimFixed',attr=nDimFixed1)
!    print *,'nDimFixed1=',nDimFixed1
    allocate(dimFixedDims1(nDimFixed1))
    if(nDimFixed1 .gt. 0)call readNcdfAttr(ncid,attrName='dimFixedDims',attr=dimFixedDims1)
!    print *,'dimFixedDims1=',dimFixedDims1
!    print *,'Finished reading dimension attributes'

    nelemFixed=product(dimFixedDims1)

! Set array size parameters needed for converting from:
! rank-4 arrays passed in to corresponding rank-2 arrays read from file
! rank-3 arrays passed in to corresponding rank-1 arrays read from file
! Internal representation of gridded data folds unlimited dims (e.g. row,col,time)
! into a single dimension to allow single-record read and replace

    ni=1
    nj=1
    nk=1
    ni=dimUnlimDims1(1)
    if(ndimUnlim1 .gt. 1)nj=dimUnlimDims1(2)
    if(ndimUnlim1 .gt. 2)nk=dimUnlimDims1(3)
    if(ndimUnlim1 .gt. 3)then
       stop 'err[putSfcGrd]: ndimUnlim1 .gt. 3 not supported'
    endif
    nelemUnlim=ni*nj*nk
!    print *,'ni,nj,nk,nelemUnlim,nelemFixed=',ni,nj,nk,nelemUnlim,nelemFixed

    if (present(xid) .and. lgetid1) &
         call readNcdfData(ncid,xid,varname='xid')

    if (present(year1d) .and. lgetyear1d1) &
         call readNcdfData(ncid,year1d,varname='year1d')

    if (present(julday1d) .and. lgetjulday1d1) &
         call readNcdfData(ncid,julday1d,varname='julday1d')
    
    if (present(x3d1_i1) .and. lgetdata1(1)) then
       if(size(x3d1_i1) .ne. nelemFixed .and. &
            size(x3d1_i1) .ne. nelemUnlim*nelemFixed)then
          print *,'err[getSfcGrd]: size(x3d1_i1) inconsistent with nelemFixed and nelemUnlim'
          print *,'size(x3d1_i1),nelemUnlim,nelemFixed=',size(x3d1_i1),nelemUnlim,nelemFixed
          stop
       endif
! based on array size, determine if we are reading entire array
! or simply reading a single record

! single record read
       if(size(x3d1_i1) .eq. nelemFixed)then
          if(.not. present(iloc) .or. .not. present(jloc) .or. .not. present(kloc))then
             print *,'err[getSfcGrd]: missing iloc,jloc, or kloc in single record read mode'
             stop
          endif
          irec=(kloc-1)*ni*nj+(jloc-1)*ni+iloc
          call readNcdfData(ncid,var=x3d1_i1,varname='x3d1_i1',&
               record_no=irec, status='single_record')
! read entire array
       else
          if(allocated(x3d1vec_i1))deallocate(x3d1vec_i1)
          allocate(x3d1vec_i1(nelemFixed,nelemUnlim))
          x3d1vec_i1=reshape(x3d1_i1,(/nelemFixed,nelemUnlim/))
          call readNcdfData(ncid,x3d1vec_i1,varname='x3d1_i1')
          x3d1_i1=reshape(x3d1vec_i1,(/nelemFixed*nelemUnlim/))
          deallocate(x3d1vec_i1)
       endif
    endif

    if (present(x3d1_i2) .and. lgetdata1(2)) then
       if(size(x3d1_i2) .ne. nelemFixed .and. &
            size(x3d1_i2) .ne. nelemUnlim*nelemFixed)then
          print *,'err[getSfcGrd]: size(x3d1_i2) inconsistent with nelemFixed and nelemUnlim'
          print *,'size(x3d1_i2),nelemUnlim,nelemFixed=',size(x3d1_i2),nelemUnlim,nelemFixed
          stop
       endif
! based on array size, determine if we are reading entire array
! or simply reading a single record
!       print *,'size(x3d1_i2)=',size(x3d1_i2)

! single record read
       if(size(x3d1_i2) .eq. nelemFixed)then
          if(.not. present(iloc) .or. .not. present(jloc) .or. .not. present(kloc))then
             print *,'err[getSfcGrd]: missing iloc,jloc, or kloc in single record read mode'
             stop
          endif
          irec=(kloc-1)*ni*nj+(jloc-1)*ni+iloc
!          print *,'iloc,jloc,kloc,irec=',iloc,jloc,kloc,irec
          call readNcdfData(ncid,var=x3d1_i2,varname='x3d1_i2',&
               record_no=irec, status='single_record')
! read entire array
       else
          if(allocated(x3d1vec_i2))deallocate(x3d1vec_i2)
          allocate(x3d1vec_i2(nelemFixed,nelemUnlim))
          x3d1vec_i2=reshape(x3d1_i2,(/nelemFixed,nelemUnlim/))
          call readNcdfData(ncid,x3d1vec_i2,varname='x3d1_i2')
          x3d1_i2=reshape(x3d1vec_i2,(/nelemFixed*nelemUnlim/))
          deallocate(x3d1vec_i2)
       endif
    endif

    if (present(x3d1_r4) .and. lgetdata1(3)) then
       if(size(x3d1_r4) .ne. nelemFixed .and. &
            size(x3d1_r4) .ne. nelemUnlim*nelemFixed)then
          print *,'err[getSfcGrd]: size(x3d1_r4) inconsistent with nelemFixed and nelemUnlim'
          print *,'size(x3d1_r4),nelemUnlim,nelemFixed=',size(x3d1_r4),nelemUnlim,nelemFixed
          stop
       endif
! based on array size, determine if we are reading entire array
! or simply reading a single record

! single record read
       if(size(x3d1_r4) .eq. nelemFixed)then
          if(.not. present(iloc) .or. .not. present(jloc) .or. .not. present(kloc))then
             print *,'err[getSfcGrd]: missing iloc,jloc, or kloc in single record read mode'
             stop
          endif
          irec=(kloc-1)*ni*nj+(jloc-1)*ni+iloc
          call readNcdfData(ncid,var=x3d1_r4,varname='x3d1_r4',&
               record_no=irec, status='single_record')
! read entire array
       else
          if(allocated(x3d1vec_r4))deallocate(x3d1vec_r4)
          allocate(x3d1vec_r4(nelemFixed,nelemUnlim))
          x3d1vec_r4=reshape(x3d1_r4,(/nelemFixed,nelemUnlim/))
          call readNcdfData(ncid,x3d1vec_r4,varname='x3d1_r4')
          x3d1_r4=reshape(x3d1vec_r4,(/nelemFixed*nelemUnlim/))
          deallocate(x3d1vec_r4)
       endif
    endif

    if (present(x2d1_i1) .and. lgetdata1(4)) then
       if(size(x2d1_i1) .ne. 1 .and. &
            size(x2d1_i1) .ne. nelemUnlim)then
          print *,'err[getSfcGrd]: size(x2d1_i1) inconsistent with nelemUnlim'
          print *,'size(x2d1_i1),nelemUnlim=',size(x2d1_i1),nelemUnlim
          stop
       endif
! based on array size, determine if we are writing entire array
! or simply replacing a single record

! single record read
       if(size(x2d1_i1) .eq. 1)then
          if(.not. present(iloc) .or. .not. present(jloc) .or. .not. present(kloc))then
             print *,'err[getSfcGrd]: missing iloc,jloc, or kloc in single record read mode'
             stop
          endif
          irec=(kloc-1)*ni*nj+(jloc-1)*ni+iloc
          call readNcdfData(ncid,var=x2d1_i1(1),varname='x2d1_i1',&
               record_no=irec,status='single_record')
! write entire array
       else
          if(allocated(x2d1vec_i1))deallocate(x2d1vec_i1)
          allocate(x2d1vec_i1(nelemUnlim))
          x2d1vec_i1=reshape(x2d1_i1,(/nelemUnlim/))
          call readNcdfData(ncid,x2d1vec_i1,varname='x2d1_i1')
          x2d1_i1=reshape(x2d1vec_i1,(/nelemUnlim/))
          deallocate(x2d1vec_i1)
       endif
    endif

    if (present(x2d1_i2) .and. lgetdata1(5)) then
       if(size(x2d1_i2) .ne. 1 .and. &
            size(x2d1_i2) .ne. nelemUnlim)then
          print *,'err[getSfcGrd]: size(x2d1_i2) inconsistent with nelemUnlim'
          print *,'size(x2d1_i2),nelemUnlim=',size(x2d1_i2),nelemUnlim
          stop
       endif
! based on array size, determine if we are writing entire array
! or simply replacing a single record

! single record read
       if(size(x2d1_i2) .eq. 1)then
          if(.not. present(iloc) .or. .not. present(jloc) .or. .not. present(kloc))then
             print *,'err[getSfcGrd]: missing iloc,jloc, or kloc in single record read mode'
             stop
          endif
          irec=(kloc-1)*ni*nj+(jloc-1)*ni+iloc
          call readNcdfData(ncid,var=x2d1_i2(1),varname='x2d1_i2',&
               record_no=irec,status='single_record')
! write entire array
       else
          if(allocated(x2d1vec_i2))deallocate(x2d1vec_i2)
          allocate(x2d1vec_i2(nelemUnlim))
          x2d1vec_i2=reshape(x2d1_i2,(/nelemUnlim/))
          call readNcdfData(ncid,x2d1vec_i2,varname='x2d1_i2')
          x2d1_i2=reshape(x2d1vec_i2,(/nelemUnlim/))
          deallocate(x2d1vec_i2)
       endif
    endif

    if (present(x2d1_r4) .and. lgetdata1(6)) then
       if(size(x2d1_r4) .ne. 1 .and. &
            size(x2d1_r4) .ne. nelemUnlim)then
          print *,'err[getSfcGrd]: size(x2d1_r4) inconsistent with nelemUnlim'
          print *,'size(x2d1_r4),nelemUnlim=',size(x2d1_r4),nelemUnlim
          stop
       endif
! based on array size, determine if we are writing entire array
! or simply replacing a single record

! single record read
       if(size(x2d1_r4) .eq. 1)then
          if(.not. present(iloc) .or. .not. present(jloc) .or. .not. present(kloc))then
             print *,'err[getSfcGrd]: missing iloc,jloc, or kloc in single record read mode'
             stop
          endif
          irec=(kloc-1)*ni*nj+(jloc-1)*ni+iloc
          call readNcdfData(ncid,var=x2d1_r4(1),varname='x2d1_r4',&
               record_no=irec,status='single_record')
! write entire array
       else
          if(allocated(x2d1vec_r4))deallocate(x2d1vec_r4)
          allocate(x2d1vec_r4(nelemUnlim))
          x2d1vec_r4=reshape(x2d1_r4,(/nelemUnlim/))
          call readNcdfData(ncid,x2d1vec_r4,varname='x2d1_r4')
          x2d1_r4=reshape(x2d1vec_r4,(/nelemUnlim/))
          deallocate(x2d1vec_r4)
       endif
    endif

    if (present(year) .and. lgetyear1) then 
       if(size(year) .ne. 1 .and. &
            size(year) .ne. nelemUnlim)then
          print *,'err[getSfcGrd]: size(year) inconsistent with nelemUnlim'
          print *,'size(year),nelemUnlim=',size(year),nelemUnlim
          stop
       endif
! based on array size, determine if we are writing entire array
! or simply replacing a single record

! single record read
       if(size(year) .eq. 1)then
          if(.not. present(iloc) .or. .not. present(jloc) .or. .not. present(kloc))then
             print *,'err[getSfcGrd]: missing iloc,jloc, or kloc in single record read mode'
             stop
          endif
          irec=(kloc-1)*ni*nj+(jloc-1)*ni+iloc
          call readNcdfData(ncid,var=year(1),varname='year',&
               record_no=irec,status='single_record')
! write entire array
       else
          if(allocated(yearvec))deallocate(yearvec)
          allocate(yearvec(nelemUnlim))
          yearvec=reshape(year,(/nelemUnlim/))
          call readNcdfData(ncid,yearvec,varname='year')
          year=reshape(yearvec,(/nelemUnlim/))
          deallocate(yearvec)
       endif
    endif

    if (present(julday) .and. lgetjulday1) then 
       if(size(julday) .ne. 1 .and. &
            size(julday) .ne. nelemUnlim)then
          print *,'err[getSfcGrd]: size(julday) inconsistent with nelemUnlim'
          print *,'size(julday),nelemUnlim=',size(julday),nelemUnlim
          stop
       endif
! based on array size, determine if we are writing entire array
! or simply replacing a single record

! single record read
       if(size(julday) .eq. 1)then
          if(.not. present(iloc) .or. .not. present(jloc) .or. .not. present(kloc))then
             print *,'err[getSfcGrd]: missing iloc,jloc, or kloc in single record read mode'
             stop
          endif
          irec=(kloc-1)*ni*nj+(jloc-1)*ni+iloc
          call readNcdfData(ncid,var=julday(1),varname='julday',&
               record_no=irec,status='single_record')
! write entire array
       else
          if(allocated(juldayvec))deallocate(juldayvec)
          allocate(juldayvec(nelemUnlim))
          juldayvec=reshape(julday,(/nelemUnlim/))
          call readNcdfData(ncid,juldayvec,varname='julday')
          julday=reshape(juldayvec,(/nelemUnlim/))
          deallocate(juldayvec)
       endif
    endif

  end subroutine getNcSfcGrd

! --- 

! Cover routine for putNcSfcGrd
! Used for arrays passed in with "original" rank 
! (i.e. 4 for (nValsPerGrid,nCol,nRow,nTime) arrays, 
!       3 for (nCol,nRow,nTime) arrays)

  subroutine put4DSfcGrd(ncid, imask, &
       ndimFixed, dimFixedDims, ndimUnlim, dimUnlimDims, &
       dimNamesFixed, dimNamesUnlim, dimUnlimName, dimNameId, &
       iloc, jloc, kloc, &
       xid, &
       x3d1_i1, x3d1_i1_longname, x3d1_i1_units, &
       x3d1_i2, x3d1_i2_longname, x3d1_i2_units, &
       x3d1_r4, x3d1_r4_longname, x3d1_r4_units, &
       x2d1_i1, x2d1_i1_longname, x2d1_i1_units, &
       x2d1_i2, x2d1_i2_longname, x2d1_i2_units, &
       x2d1_r4, x2d1_r4_longname, x2d1_r4_units, &
       year, &
       julday, &
       year1d, &
       julday1d)
    
    integer, intent (in) :: ncid
    integer, dimension(:,:,:,:), intent(in) :: imask
    integer, intent(in), optional :: ndimFixed, ndimUnlim
    integer, dimension(:), intent(in), optional :: dimFixedDims, dimUnlimDims
    character (len=*), dimension(:), intent(in), optional :: dimNamesFixed,dimNamesUnlim
    character (len=*), intent(in), optional :: dimUnlimName, dimNameId
    integer, intent(inout), optional  :: iloc, jloc, kloc

    character (len=*), intent (in), optional :: xid(:)
    character (len=*), intent (in), optional :: x3d1_i1_longname, x3d1_i1_units
    character (len=*), intent (in), optional :: x3d1_i2_longname, x3d1_i2_units
    character (len=*), intent (in), optional :: x3d1_r4_longname, x3d1_r4_units
    character (len=*), intent (in), optional :: x2d1_i1_longname, x2d1_i1_units
    character (len=*), intent (in), optional :: x2d1_i2_longname, x2d1_i2_units
    character (len=*), intent (in), optional :: x2d1_r4_longname, x2d1_r4_units
    integer*1, dimension(:,:,:,:), intent (in), optional :: x3d1_i1
    integer*1, dimension(:,:,:), intent (in), optional :: x2d1_i1
    integer*2, dimension(:,:,:,:), intent (in), optional :: x3d1_i2
    integer*2, dimension(:,:,:), intent (in), optional :: x2d1_i2
    real, dimension(:,:,:,:), intent (in), optional :: x3d1_r4
    real, dimension(:,:,:), intent (in), optional :: x2d1_r4
    integer*2, dimension(:,:,:), intent (in), optional :: year
    integer*2, dimension(:), intent (in), optional :: year1d
    real, dimension(:,:,:), intent (in), optional :: julday
    real, dimension(:), intent (in), optional :: julday1d

! Local variables

    integer*1, allocatable, dimension(:) :: x3d1vec_i1
    integer*1, allocatable, dimension(:) :: x2d1vec_i1
    integer*2, allocatable, dimension(:) :: x3d1vec_i2
    integer*2, allocatable, dimension(:) :: x2d1vec_i2
    real, allocatable, dimension(:) :: x3d1vec_r4
    real, allocatable, dimension(:) :: x2d1vec_r4
    integer*2, allocatable, dimension(:) :: yearvec
    real, allocatable, dimension(:) :: juldayvec
    integer*2, allocatable, dimension(:) :: year1dvec
    real, allocatable, dimension(:) :: julday1dvec
    character (len=8), allocatable, dimension(:) :: xidvec

    logical, dimension(6) :: lputdata1
    logical :: lputid1, lputyear1, lputjulday1, lputyear1d1, lputjulday1d1

! Dimension attributes reserved for future use
    integer :: ndimFixed1, ndimUnlim1
    integer, allocatable, dimension(:) :: dimFixedDims1, dimUnlimDims1
    character (len=200), allocatable, dimension(:) :: dimNamesFixed1,dimNamesUnlim1
    character (len=200) :: dimUnlimName1, dimNameId1
!

    integer :: iloc1, jloc1, kloc1

    character (len=200) :: x3d1_i1_longname1, x3d1_i1_units1
    character (len=200) :: x3d1_i2_longname1, x3d1_i2_units1
    character (len=200) :: x3d1_r4_longname1, x3d1_r4_units1
    character (len=200) :: x2d1_i1_longname1, x2d1_i1_units1
    character (len=200) :: x2d1_i2_longname1, x2d1_i2_units1
    character (len=200) :: x2d1_r4_longname1, x2d1_r4_units1

    integer :: nElem

! Allocate 1d arrays and initialize to missing (needed as placeholders if not actually present)

    allocate(x3d1vec_i1(1))
    allocate(x3d1vec_i2(1))
    allocate(x3d1vec_r4(1))
    allocate(x2d1vec_i1(1))
    allocate(x2d1vec_i2(1))
    allocate(x2d1vec_r4(1))
    allocate(yearvec(1))
    allocate(juldayvec(1))
    allocate(year1dvec(1))
    allocate(julday1dvec(1))
    allocate(xidvec(1))

! Initialize local variables iloc1,jloc1,kloc1 and reset if corresponding optional argument
! is present in argument list

    iloc1=0
    jloc1=0
    kloc1=0
    if(present(iloc))iloc1=iloc
    if(present(jloc))jloc1=jloc
    if(present(kloc))kloc1=kloc

! Set defaults for logical variables to FALSE 
! and then reset to TRUE if corresponding optional argument is present

    lputdata1(:)=.FALSE.
    lputid1=.FALSE.
    lputyear1=.FALSE.
    lputjulday1=.FALSE.
    lputyear1d1=.FALSE.
    lputjulday1d1=.FALSE.

    x3d1_i1_longname1='Missing' 
    x3d1_i2_longname1='Missing' 
    x3d1_r4_longname1='Missing' 
    x3d1_i1_units1='Missing'
    x3d1_i2_units1='Missing'
    x3d1_r4_units1='Missing'

    if(present(xid))then
       lputid1=.TRUE.
       if(allocated(xidvec))deallocate(xidvec)
       nElem=size(xid)
       allocate(xidvec(nElem))
       xidvec=xid
    endif
    if(present(year1d))then
       lputyear1d1=.TRUE.
       if(allocated(year1dvec))deallocate(year1dvec)
       nElem=size(year1d)
       allocate(year1dvec(nElem))
       year1dvec=year1d
    endif
    if(present(julday1d))then
       lputjulday1d1=.TRUE.
       if(allocated(julday1dvec))deallocate(julday1dvec)
       nElem=size(julday1d)
       allocate(julday1dvec(nElem))
       julday1dvec=julday1d
    endif

    if(present(year))then
       lputyear1=.TRUE.
       if(allocated(yearvec))deallocate(yearvec)
       nElem=size(year)
       allocate(yearvec(nElem))
       yearvec=reshape(year,(/nElem/))
    endif

    if(present(julday))then
       lputjulday1=.TRUE.
       if(allocated(juldayvec))deallocate(juldayvec)
       nElem=size(julday)
       allocate(juldayvec(nElem))
       juldayvec=reshape(julday,(/nElem/))
    endif

    if(present(x3d1_i1))then
       lputdata1(1)=.TRUE.
       if(present(x3d1_i1_longname))x3d1_i1_longname1=x3d1_i1_longname
       if(present(x3d1_i1_units))x3d1_i1_units1=x3d1_i1_units
       if(allocated(x3d1vec_i1))deallocate(x3d1vec_i1)
       nElem=size(x3d1_i1)
       allocate(x3d1vec_i1(nElem))
       x3d1vec_i1=reshape(x3d1_i1,(/nElem/))
    endif

    if(present(x3d1_i2))then
       lputdata1(2)=.TRUE.
       if(present(x3d1_i2_longname))x3d1_i2_longname1=x3d1_i2_longname
       if(present(x3d1_i2_units))x3d1_i2_units1=x3d1_i2_units
       if(allocated(x3d1vec_i2))deallocate(x3d1vec_i2)
       nElem=size(x3d1_i2)
       allocate(x3d1vec_i2(nElem))
       x3d1vec_i2=reshape(x3d1_i2,(/nElem/))
    endif

    if(present(x3d1_r4))then
       lputdata1(3)=.TRUE.
       if(present(x3d1_r4_longname))x3d1_r4_longname1=x3d1_r4_longname
       if(present(x3d1_r4_units))x3d1_r4_units1=x3d1_r4_units
       if(allocated(x3d1vec_r4))deallocate(x3d1vec_r4)
       nElem=size(x3d1_r4)
       allocate(x3d1vec_r4(nElem))
       x3d1vec_r4=reshape(x3d1_r4,(/nElem/))
    endif

    if(present(x2d1_i1))then
       lputdata1(4)=.TRUE.
       if(present(x2d1_i1_longname))x2d1_i1_longname1=x2d1_i1_longname
       if(present(x2d1_i1_units))x2d1_i1_units1=x2d1_i1_units
       if(allocated(x2d1vec_i1))deallocate(x2d1vec_i1)
       nElem=size(x2d1_i1)
       allocate(x2d1vec_i1(nElem))
       x2d1vec_i1=reshape(x2d1_i1,(/nElem/))
    endif

    if(present(x2d1_i2))then
       lputdata1(5)=.TRUE.
       if(present(x2d1_i2_longname))x2d1_i2_longname1=x2d1_i2_longname
       if(present(x2d1_i2_units))x2d1_i2_units1=x2d1_i2_units
       if(allocated(x2d1vec_i2))deallocate(x2d1vec_i2)
       nElem=size(x2d1_i2)
       allocate(x2d1vec_i2(nElem))
       x2d1vec_i2=reshape(x2d1_i2,(/nElem/))
    endif

    if(present(x2d1_r4))then
       lputdata1(6)=.TRUE.
       if(present(x2d1_r4_longname))x2d1_r4_longname1=x2d1_r4_longname
       if(present(x2d1_r4_units))x2d1_r4_units1=x2d1_r4_units
       if(allocated(x2d1vec_r4))deallocate(x2d1vec_r4)
       nElem=size(x2d1_r4)
       allocate(x2d1vec_r4(nElem))
       x2d1vec_r4=reshape(x2d1_r4,(/nElem/))
    endif

    call putNcSfcGrd(ncid=ncid, &
       iloc=iloc1, jloc=jloc1, kloc=kloc1, &
       lputid=lputid1, xid=xidvec, &
       lputdata=lputdata1, &
       x3d1_i1=x3d1vec_i1, x3d1_i1_longname=x3d1_i1_longname1, x3d1_i1_units=x3d1_i1_units1, &
       x3d1_i2=x3d1vec_i2, x3d1_i2_longname=x3d1_i2_longname1, x3d1_i2_units=x3d1_i2_units1, &
       x3d1_r4=x3d1vec_r4, x3d1_r4_longname=x3d1_r4_longname1, x3d1_r4_units=x3d1_r4_units1, &
       x2d1_i1=x2d1vec_i1, x2d1_i1_longname=x2d1_i1_longname1, x2d1_i1_units=x2d1_i1_units1, &
       x2d1_i2=x2d1vec_i2, x2d1_i2_longname=x2d1_i2_longname1, x2d1_i2_units=x2d1_i2_units1, &
       x2d1_r4=x2d1vec_r4, x2d1_r4_longname=x2d1_r4_longname1, x2d1_r4_units=x2d1_r4_units1, &
       lputyear=lputyear1, year=yearvec, &
       lputjulday=lputjulday1, julday=juldayvec, &
       lputyear1d=lputyear1d1, year1d=year1dvec, &
       lputjulday1d=lputjulday1d1, julday1d=julday1dvec)

    deallocate(x3d1vec_i1, x3d1vec_i2, x3d1vec_r4)
    deallocate(x2d1vec_i1, x2d1vec_i2, x2d1vec_r4)
    deallocate(yearvec, juldayvec)
    deallocate(year1dvec, julday1dvec)
    deallocate(xidvec)
    
    return

  end subroutine put4DSfcGrd

! ---

! --- 

! Cover routine for putNcSfcGrd
! Used for arrays passed in with reduced rather than "original" rank
! using reshape intrinsic function
! (i.e. 1 for (nValsPerGrid,nCol,nRow,nTime) arrays, 
!       1 for (nCol,nRow,nTime) arrays)

  subroutine put1DSfcGrd(ncid, imask, &
       ndimFixed, dimFixedDims, ndimUnlim, dimUnlimDims, &
       dimNamesFixed, dimNamesUnlim, dimUnlimName, dimNameId, &
       iloc, jloc, kloc, &
       xid, &
       x3d1_i1, x3d1_i1_longname, x3d1_i1_units, &
       x3d1_i2, x3d1_i2_longname, x3d1_i2_units, &
       x3d1_r4, x3d1_r4_longname, x3d1_r4_units, &
       x2d1_i1, x2d1_i1_longname, x2d1_i1_units, &
       x2d1_i2, x2d1_i2_longname, x2d1_i2_units, &
       x2d1_r4, x2d1_r4_longname, x2d1_r4_units, &
       year, &
       julday, &
       year1d, &
       julday1d)
    
    integer, intent (in) :: ncid
    integer, dimension(:), intent(in) :: imask
    integer, intent(in), optional :: ndimFixed, ndimUnlim
    integer, dimension(:), intent(in), optional :: dimFixedDims, dimUnlimDims
    character (len=*), dimension(:), intent(in), optional :: dimNamesFixed,dimNamesUnlim
    character (len=*), intent(in), optional :: dimUnlimName, dimNameId
    integer, intent(inout), optional  :: iloc, jloc, kloc

    character (len=*), intent (in), optional :: xid(:)
    character (len=*), intent (in), optional :: x3d1_i1_longname, x3d1_i1_units
    character (len=*), intent (in), optional :: x3d1_i2_longname, x3d1_i2_units
    character (len=*), intent (in), optional :: x3d1_r4_longname, x3d1_r4_units
    character (len=*), intent (in), optional :: x2d1_i1_longname, x2d1_i1_units
    character (len=*), intent (in), optional :: x2d1_i2_longname, x2d1_i2_units
    character (len=*), intent (in), optional :: x2d1_r4_longname, x2d1_r4_units
    integer*1, dimension(:), intent (in), optional :: x3d1_i1
    integer*1, dimension(:), intent (in), optional :: x2d1_i1
    integer*2, dimension(:), intent (in), optional :: x3d1_i2
    integer*2, dimension(:), intent (in), optional :: x2d1_i2
    real, dimension(:), intent (in), optional :: x3d1_r4
    real, dimension(:), intent (in), optional :: x2d1_r4
    integer*2, dimension(:), intent (in), optional :: year
    integer*2, dimension(:), intent (in), optional :: year1d
    real, dimension(:), intent (in), optional :: julday
    real, dimension(:), intent (in), optional :: julday1d

! Local variables
    integer*1, allocatable, dimension(:) :: x3d1vec_i1
    integer*1, allocatable, dimension(:) :: x2d1vec_i1
    integer*2, allocatable, dimension(:) :: x3d1vec_i2
    integer*2, allocatable, dimension(:) :: x2d1vec_i2
    real, allocatable, dimension(:) :: x3d1vec_r4
    real, allocatable, dimension(:) :: x2d1vec_r4
    integer*2, allocatable, dimension(:) :: yearvec
    real, allocatable, dimension(:) :: juldayvec
    integer*2, allocatable, dimension(:) :: year1dvec
    real, allocatable, dimension(:) :: julday1dvec
    character (len=8), allocatable, dimension(:) :: xidvec

    logical, dimension(6) :: lputdata1
    logical :: lputid1, lputyear1, lputjulday1, lputyear1d1, lputjulday1d1

! Dimension attributes reserved for future use
    integer :: ndimFixed1, ndimUnlim1
    integer, allocatable, dimension(:) :: dimFixedDims1, dimUnlimDims1
    character (len=200), allocatable, dimension(:) :: dimNamesFixed1,dimNamesUnlim1
    character (len=200) :: dimUnlimName1, dimNameId1
!
    integer :: iloc1, jloc1, kloc1

    character (len=200) :: x3d1_i1_longname1, x3d1_i1_units1
    character (len=200) :: x3d1_i2_longname1, x3d1_i2_units1
    character (len=200) :: x3d1_r4_longname1, x3d1_r4_units1
    character (len=200) :: x2d1_i1_longname1, x2d1_i1_units1
    character (len=200) :: x2d1_i2_longname1, x2d1_i2_units1
    character (len=200) :: x2d1_r4_longname1, x2d1_r4_units1

    integer :: nElem

! Allocate 1d arrays and initialize to missing (needed as placeholders if not actually present)

    allocate(x3d1vec_i1(1))
    allocate(x3d1vec_i2(1))
    allocate(x3d1vec_r4(1))
    allocate(x2d1vec_i1(1))
    allocate(x2d1vec_i2(1))
    allocate(x2d1vec_r4(1))
    allocate(yearvec(1))
    allocate(juldayvec(1))
    allocate(year1dvec(1))
    allocate(julday1dvec(1))
    allocate(xidvec(1))

! Initialize local variables iloc1,jloc1,kloc1 and reset if corresponding optional argument
! is present in argument list

    iloc1=0
    jloc1=0
    kloc1=0
    if(present(iloc))iloc1=iloc
    if(present(jloc))jloc1=jloc
    if(present(kloc))kloc1=kloc

! Set defaults for logical variables to FALSE 
! and then reset to TRUE if corresponding optional argument is present

    lputdata1(:)=.FALSE.
    lputid1=.FALSE.
    lputyear1=.FALSE.
    lputjulday1=.FALSE.
    lputyear1d1=.FALSE.
    lputjulday1d1=.FALSE.

    x3d1_i1_longname1='Missing' 
    x3d1_i2_longname1='Missing' 
    x3d1_r4_longname1='Missing' 
    x3d1_i1_units1='Missing'
    x3d1_i2_units1='Missing'
    x3d1_r4_units1='Missing'

    if(present(xid))then
       lputid1=.TRUE.
       if(allocated(xidvec))deallocate(xidvec)
       nElem=size(xid)
       allocate(xidvec(nElem))
       xidvec=xid
    endif
    if(present(year1d))then
       lputyear1d1=.TRUE.
       if(allocated(year1dvec))deallocate(year1dvec)
       nElem=size(year1d)
       allocate(year1dvec(nElem))
       year1dvec=year1d
    endif
    if(present(julday1d))then
       lputjulday1d1=.TRUE.
       if(allocated(julday1dvec))deallocate(julday1dvec)
       nElem=size(julday1d)
       allocate(julday1dvec(nElem))
       julday1dvec=julday1d
    endif

    if(present(year))then
       lputyear1=.TRUE.
       if(allocated(yearvec))deallocate(yearvec)
       nElem=size(year)
       allocate(yearvec(nElem))
       yearvec=reshape(year,(/nElem/))
    endif

    if(present(julday))then
       lputjulday1=.TRUE.
       if(allocated(juldayvec))deallocate(juldayvec)
       nElem=size(julday)
       allocate(juldayvec(nElem))
       juldayvec=reshape(julday,(/nElem/))
    endif

    if(present(x3d1_i1))then
       lputdata1(1)=.TRUE.
       if(present(x3d1_i1_longname))x3d1_i1_longname1=x3d1_i1_longname
       if(present(x3d1_i1_units))x3d1_i1_units1=x3d1_i1_units
       if(allocated(x3d1vec_i1))deallocate(x3d1vec_i1)
       nElem=size(x3d1_i1)
       allocate(x3d1vec_i1(nElem))
       x3d1vec_i1=reshape(x3d1_i1,(/nElem/))
    endif

    if(present(x3d1_i2))then
       lputdata1(2)=.TRUE.
       if(present(x3d1_i2_longname))x3d1_i2_longname1=x3d1_i2_longname
       if(present(x3d1_i2_units))x3d1_i2_units1=x3d1_i2_units
       if(allocated(x3d1vec_i2))deallocate(x3d1vec_i2)
       nElem=size(x3d1_i2)
       allocate(x3d1vec_i2(nElem))
       x3d1vec_i2=reshape(x3d1_i2,(/nElem/))
    endif

    if(present(x3d1_r4))then
       lputdata1(3)=.TRUE.
       if(present(x3d1_r4_longname))x3d1_r4_longname1=x3d1_r4_longname
       if(present(x3d1_r4_units))x3d1_r4_units1=x3d1_r4_units
       if(allocated(x3d1vec_r4))deallocate(x3d1vec_r4)
       nElem=size(x3d1_r4)
       allocate(x3d1vec_r4(nElem))
       x3d1vec_r4=reshape(x3d1_r4,(/nElem/))
    endif

    if(present(x2d1_i1))then
       lputdata1(4)=.TRUE.
       if(present(x2d1_i1_longname))x2d1_i1_longname1=x2d1_i1_longname
       if(present(x2d1_i1_units))x2d1_i1_units1=x2d1_i1_units
       if(allocated(x2d1vec_i1))deallocate(x2d1vec_i1)
       nElem=size(x2d1_i1)
       allocate(x2d1vec_i1(nElem))
       x2d1vec_i1=reshape(x2d1_i1,(/nElem/))
    endif

    if(present(x2d1_i2))then
       lputdata1(5)=.TRUE.
       if(present(x2d1_i2_longname))x2d1_i2_longname1=x2d1_i2_longname
       if(present(x2d1_i2_units))x2d1_i2_units1=x2d1_i2_units
       if(allocated(x2d1vec_i2))deallocate(x2d1vec_i2)
       nElem=size(x2d1_i2)
       allocate(x2d1vec_i2(nElem))
       x2d1vec_i2=reshape(x2d1_i2,(/nElem/))
    endif

    if(present(x2d1_r4))then
       lputdata1(6)=.TRUE.
       if(present(x2d1_r4_longname))x2d1_r4_longname1=x2d1_r4_longname
       if(present(x2d1_r4_units))x2d1_r4_units1=x2d1_r4_units
       if(allocated(x2d1vec_r4))deallocate(x2d1vec_r4)
       nElem=size(x2d1_r4)
       allocate(x2d1vec_r4(nElem))
       x2d1vec_r4=reshape(x2d1_r4,(/nElem/))
    endif

    print *,'put1DSfcGrd: call putNcSfcGrd'
    print *,'iloc1,jloc1,kloc1=',iloc1,jloc1,kloc1
    print *,'year,julday=',year,julday
    call putNcSfcGrd(ncid=ncid, &
       iloc=iloc1, jloc=jloc1, kloc=kloc1, &
       lputid=lputid1, xid=xidvec, &
       lputdata=lputdata1, &
       x3d1_i1=x3d1vec_i1, x3d1_i1_longname=x3d1_i1_longname1, x3d1_i1_units=x3d1_i1_units1, &
       x3d1_i2=x3d1vec_i2, x3d1_i2_longname=x3d1_i2_longname1, x3d1_i2_units=x3d1_i2_units1, &
       x3d1_r4=x3d1vec_r4, x3d1_r4_longname=x3d1_r4_longname1, x3d1_r4_units=x3d1_r4_units1, &
       x2d1_i1=x2d1vec_i1, x2d1_i1_longname=x2d1_i1_longname1, x2d1_i1_units=x2d1_i1_units1, &
       x2d1_i2=x2d1vec_i2, x2d1_i2_longname=x2d1_i2_longname1, x2d1_i2_units=x2d1_i2_units1, &
       x2d1_r4=x2d1vec_r4, x2d1_r4_longname=x2d1_r4_longname1, x2d1_r4_units=x2d1_r4_units1, &
       lputyear=lputyear1, year=yearvec, &
       lputjulday=lputjulday1, julday=juldayvec, &
       lputyear1d=lputyear1d1, year1d=year1dvec, &
       lputjulday1d=lputjulday1d1, julday1d=julday1dvec)

    deallocate(x3d1vec_i1, x3d1vec_i2, x3d1vec_r4)
    deallocate(x2d1vec_i1, x2d1vec_i2, x2d1vec_r4)
    deallocate(yearvec, juldayvec)
    deallocate(year1dvec, julday1dvec)
    deallocate(xidvec)

    return

  end subroutine put1DSfcGrd

! ---

  subroutine putNcSfcGrd(ncid, &
       ndimFixed, dimFixedDims, ndimUnlim, dimUnlimDims, &
       dimNamesFixed, dimNamesUnlim, dimUnlimName, dimNameId, &
       iloc, jloc, kloc, &
       lputid, xid, &
       lputdata, &
       x3d1_i1, x3d1_i1_longname, x3d1_i1_units, &
       x3d1_i2, x3d1_i2_longname, x3d1_i2_units, &
       x3d1_r4, x3d1_r4_longname, x3d1_r4_units, &
       x2d1_i1, x2d1_i1_longname, x2d1_i1_units, &
       x2d1_i2, x2d1_i2_longname, x2d1_i2_units, &
       x2d1_r4, x2d1_r4_longname, x2d1_r4_units, &
       lputyear, year, &
       lputjulday, julday, &
       lputyear1d, year1d, &
       lputjulday1d, julday1d)
    
    integer, intent (in) :: ncid
    integer, intent(in), optional :: ndimFixed, ndimUnlim
    integer, dimension(:), intent(in), optional :: dimFixedDims, dimUnlimDims
    character (len=*), dimension(:), intent(in), optional :: dimNamesFixed,dimNamesUnlim
    character (len=*), intent(in), optional :: dimUnlimName, dimNameId
    integer, intent(inout), optional  :: iloc, jloc, kloc

    logical, dimension(:), optional, intent (in) :: lputdata
    logical, optional, intent(in) :: lputid, lputyear, lputjulday, lputyear1d, lputjulday1d
    character (len=*), intent (in), optional :: xid(:)
    character (len=*), intent (in), optional :: x3d1_i1_longname, x3d1_i1_units
    character (len=*), intent (in), optional :: x3d1_i2_longname, x3d1_i2_units
    character (len=*), intent (in), optional :: x3d1_r4_longname, x3d1_r4_units
    character (len=*), intent (in), optional :: x2d1_i1_longname, x2d1_i1_units
    character (len=*), intent (in), optional :: x2d1_i2_longname, x2d1_i2_units
    character (len=*), intent (in), optional :: x2d1_r4_longname, x2d1_r4_units
    integer*1, dimension(:), intent (in), optional :: x3d1_i1
    integer*1, dimension(:), intent (in), optional :: x2d1_i1
    integer*2, dimension(:), intent (in), optional :: x3d1_i2
    integer*2, dimension(:), intent (in), optional :: x2d1_i2
    real, dimension(:), intent (in), optional :: x3d1_r4
    real, dimension(:), intent (in), optional :: x2d1_r4
    integer*2, dimension(:), intent (in), optional :: year, year1d
    real, dimension(:), intent (in), optional :: julday, julday1d
    ! NOTE: add other dummy arguments here as needed
    ! to support other gridded fields (e.g. Tskin, soil moisture, etc.)

    character (len=200) :: x3d_longname1, x2d_longname1
    character (len=200) :: x3d_units1, x2d_units1
    logical, dimension(6) :: lputdata1
    logical :: lputid1, lputyear1, lputjulday1, lputyear1d1, lputjulday1d1
    integer :: ni, nj, nk, nelemFixed, nelemUnlim, i, irec
    integer :: ndimUnlim1, ndimFixed1
    integer, allocatable, dimension(:) :: dimUnlimDims1, dimFixedDims1
    character(len=20), allocatable, dimension(:) :: dimNamesFixed1, dimNamesUnlim1
    character(len=50), allocatable, dimension(:) :: varLenName
    character(len=50) :: dimUnlimName1
    character(len=50) :: dimNameId1
    character(len=1) :: sAttrNum

! Local data arrays used for output to nCDF file
    
    integer*1, allocatable, dimension(:,:) :: x3d1vec_i1
    integer*1, allocatable, dimension(:) :: x2d1vec_i1
    integer*2, allocatable, dimension(:,:) :: x3d1vec_i2
    integer*2, allocatable, dimension(:) :: x2d1vec_i2
    real, allocatable, dimension(:,:) :: x3d1vec_r4
    real, allocatable, dimension(:) :: x2d1vec_r4
    integer*2, allocatable, dimension(:) :: yearvec
    real, allocatable, dimension(:) :: juldayvec


    x3d_longname1='Missing'
    x3d_units1='Missing'
    x2d_longname1='Missing'
    x2d_units1='Missing'

    lputdata1(:)=.TRUE.
    lputid1=.TRUE.
    lputyear1=.TRUE.
    lputjulday1=.TRUE.
    lputyear1d1=.TRUE.
    lputjulday1d1=.TRUE.
    if(present(lputdata))lputdata1=lputdata
    if(present(lputid))lputid1=lputid
    if(present(lputyear))lputyear1=lputyear
    if(present(lputjulday))lputjulday1=lputjulday
    if(present(lputyear1d))lputyear1d1=lputyear1d
    if(present(lputjulday1d))lputjulday1d1=lputjulday1d
    
! Set default values for dimension names if missing

    print *,'putSfcGrd: Reading dimension attributes'
    call readNcdfAttr(ncid,attrName='nDimUnlim',attr=nDimUnlim1)
    print *,'nDimUnlim1=',nDimUnlim1
    allocate(dimUnlimDims1(nDimUnlim1),dimNamesUnlim1(nDimUnlim1))
    call readNcdfAttr(ncid,attrName='dimUnlimDims',attr=dimUnlimDims1)
    print *,'dimUnlimDims1=',dimUnlimDims1
    do i=1,nDimUnlim1
       dimNamesUnlim1(i)=''
       write(sAttrNum,'(i1)')i
       call readNcdfAttr(ncid,attrName='dimNamesUnlim'//sAttrNum, &
            attr=dimNamesUnlim1(i))
       print *,'i,dimNamesUnlim1(i)=',i,trim(dimNamesUnlim1(i))
    enddo

    call readNcdfAttr(ncid,attrName='nDimFixed',attr=nDimFixed1)
    print *,'nDimFixed1=',nDimFixed1
    allocate(dimFixedDims1(nDimFixed1),dimNamesFixed1(nDimFixed1))
    call readNcdfAttr(ncid,attrName='dimFixedDims',attr=dimFixedDims1)
    print *,'dimFixedDims1=',dimFixedDims1
    do i=1,nDimFixed1
       dimNamesFixed1(i)=''
       write(sAttrNum,'(i1)')i
       call readNcdfAttr(ncid,attrName='dimNamesFixed'//sAttrNum, &
            attr=dimNamesFixed1(i))
       print *,'i,dimNamesUnlim1(i)=',i,trim(dimNamesFixed1(i))
    enddo
    dimUnlimName1=''
    call readNcdfAttr(ncid,attrName='dimUnlimName',attr=dimUnlimName1)
    print *,'dimUnlimName1=',trim(dimUnlimName1)
    print *,'Finished reading dimension attributes'
    nelemFixed=product(dimFixedDims1)

! Future extension: Insert setting dimension attributes based on input arg list here
! Currently we read directly from existing global attributes in the file
!


! Set array size parameters needed for converting from:
! rank-4 arrays passed in to corresponding rank-2 arrays written to file
! rank-3 arrays passed in to corresponding rank-1 arrays written to file
! internal representation of gridded data folds unlimited dims (e.g. row,col,time)
! into a single dimension to allow single-record read and replace

    ni=1
    nj=1
    nk=1
    ni=dimUnlimDims1(1)
    if(ndimUnlim1 .gt. 1)nj=dimUnlimDims1(2)
    if(ndimUnlim1 .gt. 2)nk=dimUnlimDims1(3)
    if(ndimUnlim1 .gt. 3)then
       stop 'err[putSfcGrd]: ndimUnlim1 .gt. 3 not supported'
    endif
    nelemUnlim=ni*nj*nk
    print *,'ni,nj,nk,nelemUnlim,nelemFixed=',ni,nj,nk,nelemUnlim,nelemFixed

    print *,'putNcSfcGrd: iloc,jloc,kloc=',iloc,jloc,kloc

    if (present(xid) .and. lputid1) then
       print *,'sfcgrd_io_module:xid present, writing xid'
       
       if(.not. present(dimNameId))then
          dimNameId1='nTimeLevels'
       else
          dimNameId1=trim(dimNameId)
       endif

       call writeNcdfData(ncid,xid,varname='xid', &
            strLenName='nxId', &
            varLenName=dimNameId1,&
            varLongName='gridded field identifier', &
            varUnit='None')
    endif

    if (present(year1d) .and. lputyear1d1) then
       print *,'sfcgrd_io_module:year1d present, writing year1d'
       
       if(.not. present(dimNameId))then
          dimNameId1='nTimeLevels'
       else
          dimNameId1=trim(dimNameId)
       endif

       call writeNcdfData(ncid,year1d,varname='year1d', &
            varLenName=dimNameId1,&
            varLongName='1-d year, one per time level', &
            varUnit='years')
    endif

    if (present(julday1d) .and. lputjulday1d1) then
       print *,'sfcgrd_io_module:julday1d present, writing julday1d'
       
       if(.not. present(dimNameId))then
          dimNameId1='nTimeLevels'
       else
          dimNameId1=trim(dimNameId)
       endif

       call writeNcdfData(ncid,julday1d,varname='julday1d', &
            varLenName=dimNameId1,&
            varLongName='1-d julian day, one per time level', &
            varUnit='julian days')
    endif

    
    if (present(x3d1_i1) .and. lputdata1(1)) then
       print *,'sfcgrd_io_module: x3d1_i1 present, writing data'

! based on array size, determine if we are writing entire array
! or simply replacing a single record

       if(size(x3d1_i1) .ne. nelemFixed .and. &
            size(x3d1_i1) .ne. nelemUnlim*nelemFixed)then
          print *,'err[putSfcGrd]: size(x3d1_i1) inconsistent with nelemFixed and nelemUnlim'
          print *,'size(x3d1_i1),nelemUnlim,nelemFixed=',size(x3d1_i1),nelemUnlim,nelemFixed
          stop
       endif
! single record replacement
       if(size(x3d1_i1) .eq. nelemFixed)then
          if(.not. present(iloc) .or. .not. present(jloc) .or. .not. present(kloc))then
             print *,'err[putSfcGrd]: missing iloc,jloc, or kloc in single record replacement mode'
             stop
          endif
          if(allocated(x3d1vec_i1))deallocate(x3d1vec_i1)
          allocate(x3d1vec_i1(nelemFixed,1))
          x3d1vec_i1=reshape(x3d1_i1,(/nelemFixed,1/))
          irec=(kloc-1)*ni*nj+(jloc-1)*ni+iloc
          call replaceNcdfData(ncid,var=x3d1vec_i1,varname='x3d1_i1',rec=irec)
          deallocate(x3d1vec_i1)
! write entire array
       else
          if(allocated(x3d1vec_i1))deallocate(x3d1vec_i1)
          allocate(x3d1vec_i1(nelemFixed,nelemUnlim))
          x3d1vec_i1=reshape(x3d1_i1,(/nelemFixed,nelemUnlim/))

          if(present(x3d1_i1_longname))x3d_longname1=trim(x3d1_i1_longname)
          if(present(x3d1_i1_units))x3d_units1=trim(x3d1_i1_units)
          print *,'x3d_units1=',trim(x3d_units1)
          print *,'x3d_longname1=',trim(x3d_longname1)
! construct string array for varLenName 
          if(allocated(varLenName))deallocate(varLenName)
          allocate(varLenName(ndimFixed1+1))
          do i=1,ndimFixed1
             varLenName(i)=trim(dimNamesFixed1(i))
          enddo
          varLenName(ndimFixed1+1)=trim(dimUnlimName1)
          call writeNcdfData(ncid,x3d1vec_i1,varname='x3d1_i1', &
               varLenName=varLenName,&
               varLongName=trim(x3d_longname1), &
               varUnit=trim(x3d_units1))
          deallocate(x3d1vec_i1)
          deallocate(varLenName)
       endif
    endif

    if (present(x3d1_i2) .and. lputdata1(2)) then
       print *,'sfcgrd_io_module: x3d1_i2 present, writing data'

! based on array size, determine if we are writing entire array
! or simply replacing a single record

       if(size(x3d1_i2) .ne. nelemFixed .and. &
            size(x3d1_i2) .ne. nelemUnlim*nelemFixed)then
          print *,'err[putSfcGrd]: size(x3d1_i2) inconsistent with nelemFixed and nelemUnlim'
          print *,'size(x3d1_i2),nelemUnlim,nelemFixed=',size(x3d1_i2),nelemUnlim,nelemFixed
          stop
       endif
! single record replacement
       if(size(x3d1_i2) .eq. nelemFixed)then
          print *,'putNcSfcGrd:single record replacement'
          print *,'size(x3d1_i2)=',size(x3d1_i2)
          if(.not. present(iloc) .or. .not. present(jloc) .or. .not. present(kloc))then
             print *,'err[putSfcGrd]: missing iloc,jloc, or kloc in single record replacement mode'
             stop
          endif
          if(allocated(x3d1vec_i2))deallocate(x3d1vec_i2)
          allocate(x3d1vec_i2(nelemFixed,1))
          x3d1vec_i2=reshape(x3d1_i2,(/nelemFixed,1/))
          irec=(kloc-1)*ni*nj+(jloc-1)*ni+iloc
          call replaceNcdfData(ncid,var=x3d1vec_i2,varname='x3d1_i2',rec=irec)
          deallocate(x3d1vec_i2)
! write entire array
       else
          if(allocated(x3d1vec_i2))deallocate(x3d1vec_i2)
          allocate(x3d1vec_i2(nelemFixed,nelemUnlim))
          x3d1vec_i2=reshape(x3d1_i2,(/nelemFixed,nelemUnlim/))

          if(present(x3d1_i2_longname))x3d_longname1=trim(x3d1_i2_longname)
          if(present(x3d1_i2_units))x3d_units1=trim(x3d1_i2_units)
          print *,'x3d_units1=',trim(x3d_units1)
          print *,'x3d_longname1=',trim(x3d_longname1)
! construct string array for varLenName 
          if(allocated(varLenName))deallocate(varLenName)
          allocate(varLenName(ndimFixed1+1))
          do i=1,ndimFixed1
             varLenName(i)=trim(dimNamesFixed1(i))
          enddo
          varLenName(ndimFixed1+1)=trim(dimUnlimName1)
          call writeNcdfData(ncid,x3d1vec_i2,varname='x3d1_i2', &
               varLenName=varLenName,&
               varLongName=trim(x3d_longname1), &
               varUnit=trim(x3d_units1))
          deallocate(x3d1vec_i2)
          deallocate(varLenName)
       endif
    endif

    if (present(x3d1_r4) .and. lputdata1(3)) then
       print *,'sfcgrd_io_module: x3d1_r4 present, writing data'

! based on array size, determine if we are writing entire array
! or simply replacing a single record

       if(size(x3d1_r4) .ne. nelemFixed .and. &
            size(x3d1_r4) .ne. nelemUnlim*nelemFixed)then
          print *,'err[putSfcGrd]: size(x3d1_r4) inconsistent with nelemFixed and nelemUnlim'
          print *,'size(x3d1_r4),nelemUnlim,nelemFixed=',size(x3d1_r4),nelemUnlim,nelemFixed
          stop
       endif
! single record replacement
       if(size(x3d1_r4) .eq. nelemFixed)then
          print *,'putNcSfcGrd:single record replacement'
          print *,'size(x3d1_r4)=',size(x3d1_r4)
          if(.not. present(iloc) .or. .not. present(jloc) .or. .not. present(kloc))then
             print *,'err[putSfcGrd]: missing iloc,jloc, or kloc in single record replacement mode'
             stop
          endif
          if(allocated(x3d1vec_r4))deallocate(x3d1vec_r4)
          allocate(x3d1vec_r4(nelemFixed,1))
          x3d1vec_r4=reshape(x3d1_r4,(/nelemFixed,1/))
          irec=(kloc-1)*ni*nj+(jloc-1)*ni+iloc
          call replaceNcdfData(ncid,var=x3d1vec_r4,varname='x3d1_r4',rec=irec)
          deallocate(x3d1vec_r4)
! write entire array
       else
          print *,'writing entire array'
          if(allocated(x3d1vec_r4))deallocate(x3d1vec_r4)
          allocate(x3d1vec_r4(nelemFixed,nelemUnlim))
          x3d1vec_r4=reshape(x3d1_r4,(/nelemFixed,nelemUnlim/))

          if(present(x3d1_r4_longname))x3d_longname1=trim(x3d1_r4_longname)
          if(present(x3d1_r4_units))x3d_units1=trim(x3d1_r4_units)
          print *,'x3d_units1=',trim(x3d_units1)
          print *,'x3d_longname1=',trim(x3d_longname1)
! construct string array for varLenName 
          print *,'construct string array for varLenName'
          if(allocated(varLenName))deallocate(varLenName)
          allocate(varLenName(ndimFixed1+1))
          do i=1,ndimFixed1
             varLenName(i)=trim(dimNamesFixed1(i))
          enddo
          varLenName(ndimFixed1+1)=trim(dimUnlimName1)
          print *,'varLenName(1)=',trim(varLenName(1))
          print *,'varLenName(2)=',trim(varLenName(2))
          call writeNcdfData(ncid,x3d1vec_r4,varname='x3d1_r4', &
               varLenName=varLenName,&
               varLongName=trim(x3d_longname1), &
               varUnit=trim(x3d_units1))
          deallocate(x3d1vec_r4)
          deallocate(varLenName)
       endif
    endif


    if (present(x2d1_i1) .and. lputdata1(4)) then
       print *,'sfcgrd_io_module: x2d1_i1 present, writing data'

! based on array size, determine if we are writing entire array
! or simply replacing a single record

       if(size(x2d1_i1) .ne. 1 .and. &
            size(x2d1_i1) .ne. nelemUnlim)then
          print *,'err[putSfcGrd]: size(x2d1_i1) inconsistent with nelemUnlim'
          print *,'size(x2d1_i1),nelemUnlim=',size(x2d1_i1),nelemUnlim
          stop
       endif
! single record replacement
       if(size(x2d1_i1) .eq. 1)then
          if(.not. present(iloc) .or. .not. present(jloc) .or. .not. present(kloc))then
             print *,'err[putSfcGrd]: missing iloc,jloc, or kloc in single record replacement mode'
             stop
          endif
          if(allocated(x2d1vec_i1))deallocate(x2d1vec_i1)
          allocate(x2d1vec_i1(1))
          x2d1vec_i1=reshape(x2d1_i1,(/1/))
          irec=(kloc-1)*ni*nj+(jloc-1)*ni+iloc
          call replaceNcdfData(ncid,var=x2d1vec_i1,varname='x2d1_i1',rec=irec)
          deallocate(x2d1vec_i1)
! write entire array
       else
          if(allocated(x2d1vec_i1))deallocate(x2d1vec_i1)
          allocate(x2d1vec_i1(nelemUnlim))
          x2d1vec_i1=reshape(x2d1_i1,(/nelemUnlim/))

          if(present(x2d1_i1_longname))x2d_longname1=trim(x2d1_i1_longname)
          if(present(x2d1_i1_units))x2d_units1=trim(x2d1_i1_units)
          print *,'x2d_units1=',trim(x2d_units1)
          print *,'x2d_longname1=',trim(x2d_longname1)
! construct string array for varLenName 
          if(allocated(varLenName))deallocate(varLenName)
          allocate(varLenName(1))
          varLenName(1)=trim(dimUnlimName1)
          call writeNcdfData(ncid,x2d1vec_i1,varname='x2d1_i1', &
               varLenName=varLenName(1),&
               varLongName=trim(x2d_longname1), &
               varUnit=trim(x2d_units1))
          deallocate(x2d1vec_i1)
          deallocate(varLenName)
       endif
    endif

    if (present(x2d1_i2) .and. lputdata1(5)) then
       print *,'sfcgrd_io_module: x2d1_i2 present, writing data'

! based on array size, determine if we are writing entire array
! or simply replacing a single record

       if(size(x2d1_i2) .ne. 1 .and. &
            size(x2d1_i2) .ne. nelemUnlim)then
          print *,'err[putSfcGrd]: size(x2d1_i2) inconsistent with nelemUnlim'
          print *,'size(x2d1_i2),nelemUnlim=',size(x2d1_i2),nelemUnlim
          stop
       endif
! single record replacement
       if(size(x2d1_i2) .eq. 1)then
          print *,'putNcSfcGrd:single record replacement'
          print *,'size(x2d1_i2)=',size(x2d1_i2)
          if(.not. present(iloc) .or. .not. present(jloc) .or. .not. present(kloc))then
             print *,'err[putSfcGrd]: missing iloc,jloc, or kloc in single record replacement mode'
             stop
          endif
          if(allocated(x2d1vec_i2))deallocate(x2d1vec_i2)
          allocate(x2d1vec_i2(1))
          x2d1vec_i2=reshape(x2d1_i2,(/1/))
          irec=(kloc-1)*ni*nj+(jloc-1)*ni+iloc
          call replaceNcdfData(ncid,var=x2d1vec_i2,varname='x2d1_i2',rec=irec)
          deallocate(x2d1vec_i2)
! write entire array
       else
          if(allocated(x2d1vec_i2))deallocate(x2d1vec_i2)
          allocate(x2d1vec_i2(nelemUnlim))
          x2d1vec_i2=reshape(x2d1_i2,(/nelemUnlim/))

          if(present(x2d1_i2_longname))x2d_longname1=trim(x2d1_i2_longname)
          if(present(x2d1_i2_units))x2d_units1=trim(x2d1_i2_units)
          print *,'x2d_units1=',trim(x2d_units1)
          print *,'x2d_longname1=',trim(x2d_longname1)
! construct string array for varLenName 
          if(allocated(varLenName))deallocate(varLenName)
          allocate(varLenName(1))
          varLenName(1)=trim(dimUnlimName1)
          call writeNcdfData(ncid,x2d1vec_i2,varname='x2d1_i2', &
               varLenName=varLenName(1),&
               varLongName=trim(x2d_longname1), &
               varUnit=trim(x2d_units1))
          deallocate(x2d1vec_i2)
          deallocate(varLenName)
       endif
    endif

    if (present(x2d1_r4) .and. lputdata1(6)) then
       print *,'sfcgrd_io_module: x2d1_r4 present, writing data'

! based on array size, determine if we are writing entire array
! or simply replacing a single record

       if(size(x2d1_r4) .ne. 1 .and. &
            size(x2d1_r4) .ne. nelemUnlim)then
          print *,'err[putSfcGrd]: size(x2d1_r4) inconsistent with nelemUnlim'
          print *,'size(x2d1_r4),nelemUnlim=',size(x2d1_r4),nelemUnlim
          stop
       endif
! single record replacement
       if(size(x2d1_r4) .eq. 1)then
          if(.not. present(iloc) .or. .not. present(jloc) .or. .not. present(kloc))then
             print *,'err[putSfcGrd]: missing iloc,jloc, or kloc in single record replacement mode'
             stop
          endif
          if(allocated(x2d1vec_r4))deallocate(x2d1vec_r4)
          allocate(x2d1vec_r4(1))
          x2d1vec_r4=reshape(x2d1_r4,(/1/))
          irec=(kloc-1)*ni*nj+(jloc-1)*ni+iloc
          call replaceNcdfData(ncid,var=x2d1vec_r4,varname='x2d1_r4',rec=irec)
          deallocate(x2d1vec_r4)
! write entire array
       else
          if(allocated(x2d1vec_r4))deallocate(x2d1vec_r4)
          allocate(x2d1vec_r4(nelemUnlim))
          x2d1vec_r4=reshape(x2d1_r4,(/nelemUnlim/))

          if(present(x2d1_r4_longname))x2d_longname1=trim(x2d1_r4_longname)
          if(present(x2d1_r4_units))x2d_units1=trim(x2d1_r4_units)
          print *,'x2d_units1=',trim(x2d_units1)
          print *,'x2d_longname1=',trim(x2d_longname1)
! construct string array for varLenName 
           if(allocated(varLenName))deallocate(varLenName)
          allocate(varLenName(1))
          varLenName(1)=trim(dimUnlimName1)
          call writeNcdfData(ncid,x2d1vec_r4,varname='x2d1_r4', &
               varLenName=varLenName(1),&
               varLongName=trim(x2d_longname1), &
               varUnit=trim(x2d_units1))
          deallocate(x2d1vec_r4)
          deallocate(varLenName)
       endif
    endif

    if (present(year) .and. lputyear1) then
       print *,'sfcgrd_io_module:year present, writing year'
 
! based on array size, determine if we are writing entire array
! or simply replacing a single record

       if(size(year) .ne. 1 .and. &
            size(year) .ne. nelemUnlim)then
          print *,'err[putSfcGrd]: size(year) inconsistent with nelemUnlim'
          print *,'size(year),nelemUnlim=',size(year),nelemUnlim
          stop
       endif
! single record replacement
       if(size(year) .eq. 1)then
          if(.not. present(iloc) .or. .not. present(jloc) .or. .not. present(kloc))then
             print *,'err[putSfcGrd]: missing iloc,jloc, or kloc in single record replacement mode'
             stop
          endif
          if(allocated(yearvec))deallocate(yearvec)
          allocate(yearvec(1))
          yearvec=reshape(year,(/1/))
          irec=(kloc-1)*ni*nj+(jloc-1)*ni+iloc
          call replaceNcdfData(ncid,var=yearvec,varname='year',rec=irec)
          deallocate(yearvec)
! write entire array
       else
          if(allocated(yearvec))deallocate(yearvec)
          allocate(yearvec(nelemUnlim))
          yearvec=reshape(year,(/nelemUnlim/))

! construct string array for varLenName 
          if(allocated(varLenName))deallocate(varLenName)
          allocate(varLenName(1))
          varLenName(1)=trim(dimUnlimName1)
          call writeNcdfData(ncid,yearvec,varname='year', &
               varLenName=varLenName(1),&
               varLongName='valid time of gridded data (year)', &
               varUnit='years')
          deallocate(yearvec)
          deallocate(varLenName)
       endif
    endif

    if (present(julday) .and. lputjulday1) then
       print *,'sfcgrd_io_module:julday present, writing julday'
 
! based on array size, determine if we are writing entire array
! or simply replacing a single record

       if(size(julday) .ne. 1 .and. &
            size(julday) .ne. nelemUnlim)then
          print *,'err[putSfcGrd]: size(julday) inconsistent with nelemUnlim'
          print *,'size(julday),nelemUnlim=',size(julday),nelemUnlim
          stop
       endif
! single record replacement
       if(size(julday) .eq. 1)then
          if(.not. present(iloc) .or. .not. present(jloc) .or. .not. present(kloc))then
             print *,'err[putSfcGrd]: missing iloc,jloc, or kloc in single record replacement mode'
             stop
          endif
          if(allocated(juldayvec))deallocate(juldayvec)
          allocate(juldayvec(1))
          juldayvec=reshape(julday,(/1/))
          irec=(kloc-1)*ni*nj+(jloc-1)*ni+iloc
          call replaceNcdfData(ncid,var=juldayvec,varname='julday',rec=irec)
          deallocate(juldayvec)
! write entire array
       else
          if(allocated(juldayvec))deallocate(juldayvec)
          allocate(juldayvec(nelemUnlim))
          juldayvec=reshape(julday,(/nelemUnlim/))

! construct string array for varLenName 
          if(allocated(varLenName))deallocate(varLenName)
          allocate(varLenName(1))
          varLenName(1)=trim(dimUnlimName1)
          call writeNcdfData(ncid,juldayvec,varname='julday', &
               varLenName=varLenName(1),&
               varLongName='valid time of gridded data (julday)', &
               varUnit='julian day, including fraction')
          deallocate(juldayvec)
          deallocate(varLenName)
       endif
    endif


  end subroutine putNcSfcGrd
! -
  subroutine closeNcGFile(ncid,msg)
    character (len=*), intent (in), optional :: msg
    integer, intent (in) :: ncid
    
    call closeNcdfFile(ncid, msg)
    
  end subroutine closeNcGFile
! ---
end module sfcgrd_io_module

