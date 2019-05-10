 
! "$Id$"

! Module for managing interface with surface grid database files

module sfcgrd_DB_module

  use sfcgrd_io_module
  use sfcgrdDBtype_module
  use bilintSintype_module
  use mapConvert
  use FOV2sfcgrd_module
  use constants, ONLY : MISSING_REAL, MISSING_INT, MISSING_CHAR

  implicit none
  private
  public getSfcDataBase, getSfcDBparams
  public getSfcDBlocs
  public getSfcLandType, getSfcWaterType, getSfcElev, getSfcSeaIce
  public getBilintParams
  public getTileLatRange
  public getSfcDBWindowLocs
  public destroySfcgrdDB

! params for "equivalent global grid" shared within module here? (NOT passed as params)

  ! ---
CONTAINS
  ! ---

  subroutine getSfcDataBase(fnameDBFileList,lgetDBparams,lgetDBlocs,lverbose,&
       sfcgrdDB_current,sfcgrdDB_path)
    character(len=256), intent(in) :: fnameDBFileList
    character(len=256), optional, intent(in) :: sfcgrdDB_path  !path to sfc type data files
    logical, intent(in) :: lgetDBparams, lgetDBlocs
    logical, optional, intent(in) :: lverbose

    type(sfcgrdDB), optional, intent(inout) :: sfcgrdDB_current

! local variables

    logical :: lverbose0

    lverbose0=.FALSE.
    if(present(lverbose))lverbose0=lverbose

    if(lgetDBparams) call getSfcDBparams(fnameDBFileList,sfcgrdDB_current=sfcgrdDB_current, &
         sfcgrdDB_path=sfcgrdDB_path,lverbose=lverbose0)

!    if(lgetDBlocs) call getSfcDBlocs()


  end subroutine getSfcDataBase


  subroutine getSfcDBparams(fnameDBFileList,&
       DBTileListCol,DBTileListRow,lverbose, &
       sfcgrdDB_current,sfcgrdDB_path)

    character(len=256), optional, intent(in) :: fnameDBFileList
    integer, dimension(:), optional, intent(in) :: DBTileListCol
    integer, dimension(:), optional, intent(in) :: DBTileListRow
    character(len=256), optional, intent(in) :: sfcgrdDB_path  !path to sfc type data files
    logical, optional, intent(in) :: lverbose

    type(sfcgrdDB), intent(out) :: sfcgrdDB_current

! local variables
  
    character(len=512) :: fnameDBFileList0
    character(len=256) :: sfcgrdDB_path0
    character(len=256) :: filename
    integer :: i
    integer :: iuFileList=11

! use existing DB tilefiles
    if(present(fnameDBFileList))then

       if (present(sfcgrdDB_path)) then
          sfcgrdDB_path0 = sfcgrdDB_path
       else
          sfcgrdDB_path0 = "."
       endif

       fnameDBFileList0=trim(sfcgrdDB_path0)//"/"//trim(fnameDBFileList)

! first read list of DB files from ascii file

       print *,'getSfcDBparams: Reading file lists'
! Read in list of input tile files  
       print *,'getSfcDBparams: Tile file list located in: ',trim(fnameDBFileList0)
       open(iuFileList,file=trim(fnameDBFileList0),form='formatted',status='old')
       sfcgrdDB_current%listfile_tilefiles=trim(fnameDBFileList0)
       read(iuFileList,*) sfcgrdDB_current%ntilefiles
       print *, 'getSfcDBparams: sfcgrdDB_current%ntilefiles=',sfcgrdDB_current%ntilefiles

       if(allocated(sfcgrdDB_current%tilefile))deallocate(sfcgrdDB_current%tilefile)
       allocate(sfcgrdDB_current%tilefile(sfcgrdDB_current%ntilefiles))
       do i=1,sfcgrdDB_current%ntilefiles
          read(iuFileList,'(a)') filename
          sfcgrdDB_current%tilefile(i) = trim(sfcgrdDB_path0)//"/"//trim(filename)
          if(present(lverbose))then
             if(lverbose)print *,'i,sfcgrdDB_current%tilefile(i)=',i,trim(sfcgrdDB_current%tilefile(i))
          endif
       enddo
       close(iuFileList)

! allocate memory for output data (derived type)

       if(allocated(sfcgrdDB_current%tile_colglobal))deallocate(sfcgrdDB_current%tile_colglobal)
       allocate(sfcgrdDB_current%tile_colglobal(sfcgrdDB_current%ntilefiles))
       if(allocated(sfcgrdDB_current%tile_rowglobal))deallocate(sfcgrdDB_current%tile_rowglobal)
       allocate(sfcgrdDB_current%tile_rowglobal(sfcgrdDB_current%ntilefiles))
       if(allocated(sfcgrdDB_current%tile_ncolglobal))deallocate(sfcgrdDB_current%tile_ncolglobal)
       allocate(sfcgrdDB_current%tile_ncolglobal(sfcgrdDB_current%ntilefiles))
       if(allocated(sfcgrdDB_current%tile_nrowglobal))deallocate(sfcgrdDB_current%tile_nrowglobal)
       allocate(sfcgrdDB_current%tile_nrowglobal(sfcgrdDB_current%ntilefiles))
       if(allocated(sfcgrdDB_current%tile_ncol))deallocate(sfcgrdDB_current%tile_ncol)
       allocate(sfcgrdDB_current%tile_ncol(sfcgrdDB_current%ntilefiles))
       if(allocated(sfcgrdDB_current%tile_nrow))deallocate(sfcgrdDB_current%tile_nrow)
       allocate(sfcgrdDB_current%tile_nrow(sfcgrdDB_current%ntilefiles))
       if(allocated(sfcgrdDB_current%tile_colGridOffset))deallocate(sfcgrdDB_current%tile_colGridOffset)
       allocate(sfcgrdDB_current%tile_colGridOffset(sfcgrdDB_current%ntilefiles))
       if(allocated(sfcgrdDB_current%tile_rowGridOffset))deallocate(sfcgrdDB_current%tile_rowGridOffset)
       allocate(sfcgrdDB_current%tile_rowGridOffset(sfcgrdDB_current%ntilefiles))
       if(allocated(sfcgrdDB_current%tile_Re))deallocate(sfcgrdDB_current%tile_Re)
       allocate(sfcgrdDB_current%tile_Re(sfcgrdDB_current%ntilefiles))
       if(allocated(sfcgrdDB_current%tile_scale))deallocate(sfcgrdDB_current%tile_scale)
       allocate(sfcgrdDB_current%tile_scale(sfcgrdDB_current%ntilefiles))

! then get grid params for each tile file and return in output variables

       call getTileParams(tilefile=sfcgrdDB_current%tilefile, &
            tile_col=sfcgrdDB_current%tile_colglobal, &
            tile_row=sfcgrdDB_current%tile_rowglobal, &
            nx_globaltiles=sfcgrdDB_current%tile_ncolglobal, &
            ny_globaltiles=sfcgrdDB_current%tile_nrowglobal, &
            nx_tile=sfcgrdDB_current%tile_ncol, &
            ny_tile=sfcgrdDB_current%tile_nrow, &
            colGridOffset=sfcgrdDB_current%tile_colGridOffset, &
            rowGridOffset=sfcgrdDB_current%tile_rowGridOffset, &
            Re=sfcgrdDB_current%tile_Re, &
            scale=sfcgrdDB_current%tile_scale)


       if(present(lverbose))then
          if(lverbose)then
             print *,'getSfcDBparams: sfcgrdDB_current%listfile_tilefiles=',&
                  trim(sfcgrdDB_current%listfile_tilefiles)
             print *,'getSfcDBparams: sfcgrdDB_current%ntilefiles=',sfcgrdDB_current%ntilefiles
             do i=1,sfcgrdDB_current%ntilefiles
                print *,'i,sfcgrdDB_current%tilefile(i)=',i,trim(sfcgrdDB_current%tilefile(i))
             enddo
             print *,'sfcgrdDB_current%tile_ncolglobal=',sfcgrdDB_current%tile_ncolglobal
             print *,'sfcgrdDB_current%tile_nrowglobal=',sfcgrdDB_current%tile_nrowglobal
             print *,'sfcgrdDB_current%tile_ncol=',sfcgrdDB_current%tile_ncol
             print *,'sfcgrdDB_current%tile_nrow=',sfcgrdDB_current%tile_nrow
             print *,'sfcgrdDB_current%tile_colglobal=',sfcgrdDB_current%tile_colglobal
             print *,'sfcgrdDB_current%tile_rowglobal=',sfcgrdDB_current%tile_rowglobal
             print *,'sfcgrdDB_current%tile_colGridOffset=',sfcgrdDB_current%tile_colGridOffset
             print *,'sfcgrdDB_current%tile_rowGridOffset=',sfcgrdDB_current%tile_rowGridOffset
             print *,'sfcgrdDB_current%tile_Re=',sfcgrdDB_current%tile_Re
             print *,'sfcgrdDB_current%tile_scale=',sfcgrdDB_current%tile_scale
          endif
       endif
       return

! else create new DB tilefiles
    elseif(present(DBTileListCol) .and. present(DBTileListCol))then
       print *,'err[getSfcDBparams]: DBTileListCol and DBTileListRow argument option not implemented'
       call errorHalt(1)       
!       call getTileList()

! then set grid params for each tile file
!       call setTileParams()

    else
       print *,'err[getSfcDBparams]: both fnameDBFileList and fnameDBTileList missing'
       call errorHalt(1)       
    endif



  end subroutine getSfcDBparams


  subroutine getSfcDBlocs(sfcgrdDB_current, lat, lon, bilintSinusoid, &
       ifileDB, ifileDB_bilint, &
       ncolGlobalGrid, nrowGlobalGrid, &
       colGridOffsetGlobal, rowGridOffsetGlobal, &
       igrdGlobal, jgrdGlobal, &
       igrdLocal, jgrdLocal, &
       bilintParams, &
       colTile, rowTile)

    type(sfcgrdDB), intent(in) :: sfcgrdDB_current
    real, dimension(:), intent(in) :: lat, lon
    logical, intent(in) :: bilintSinusoid

    integer, optional, intent(out) :: ncolGlobalGrid, nrowGlobalGrid
    real, optional, intent(out) :: colGridOffsetGlobal, rowGridOffsetGlobal

    integer, dimension(:), pointer, optional :: ifileDB
    integer, dimension(:,:), pointer, optional :: ifileDB_bilint
    real, dimension(:), pointer, optional :: igrdGlobal, jgrdGlobal
    real, dimension(:), pointer, optional :: igrdLocal, jgrdLocal
    integer, dimension(:), pointer, optional :: colTile, rowTile

! needed for bilinear interpolation option
    type(bilintSin), optional, intent(out) :: bilintParams
! two sets of arrays needed due to: 
! (1) shift of grid values w/r to array locations (i.e. value stored at array(1,1) is valid at grid(0.5,0.5))
! (2) possibility of wraparound at dateline. Required array locations will not be consistent with grid values
! arrays for storing *grid* locations (non-integer values, needed for finite differences)
!    real, dimension(:), pointer, optional :: i0floorG, i0ceilG, jfloorG, jceilG
!    real, dimension(:), pointer, optional :: i_jfloorG, i_jceilG
! arrays for storing *array* index locations in local tiles (from tile files)
!    integer, dimension(:), pointer, optional :: i0floorALocal, i1floorALocal
!    integer, dimension(:), pointer, optional :: i0ceilALocal, i1ceilALocal
!    integer, dimension(:), pointer, optional :: jfloorALocal, jceilALocal

    
! Local variables

    integer :: i, j, nlatlon
    integer :: ncolGlobalGrid0, nrowGlobalGrid0
    real :: colGridOffsetGlobal0, rowGridOffsetGlobal0
    real :: scale0, Re0
    real, dimension(:), allocatable :: igrdGlobal0, jgrdglobal0
!    real, dimension(:), allocatable :: igrdLocal0, jgrdLocal0
!    integer, dimension(:), allocatable :: colTile0, rowTile0

    integer, parameter :: nBilint=4
    

!    print *,'Begin getSfcDBLocs:'
!    print *,'sfcgrdDB_current%ntilefiles=',sfcgrdDB_current%ntilefiles

! Make sure database has at least one tile

    if(sfcgrdDB_current%ntilefiles .lt. 1)then
       print *,'err[getSfcDBLocs]: ntilefiles must be ge 1'
       call errorHalt(1)       
    endif


! Calculate the grid parameters for an input grid with same scale as 
! an individual input tile, but with global coverage

    ncolGlobalGrid0=sfcgrdDB_current%tile_ncolglobal(1)*sfcgrdDB_current%tile_ncol(1)
    nrowGlobalGrid0=sfcgrdDB_current%tile_nrowglobal(1)*sfcgrdDB_current%tile_nrow(1)
    colGridOffsetGlobal0=(real(ncolGlobalGrid0)/2.)
    rowGridOffsetGlobal0=(real(nrowGlobalGrid0)/2.)
    scale0=sfcgrdDB_current%tile_scale(1)
    Re0=sfcgrdDB_current%tile_Re(1)

! Total number of lat/lon locations

    nlatlon=size(lat)

! Allocate local arrays

    if(allocated(igrdGlobal0))deallocate(igrdGlobal0)
    allocate(igrdGlobal0(nlatlon))
    if(allocated(jgrdGlobal0))deallocate(jgrdGlobal0)
    allocate(jgrdGlobal0(nlatlon))
!    if(allocated(igrdLocal0))deallocate(igrdLocal0)
!    allocate(igrdLocal0(nlatlon))
!    if(allocated(jgrdLocal0))deallocate(jgrdLocal0)
!    allocate(jgrdLocal0(nlatlon))
!    if(allocated(colTile0))deallocate(colTile0)
!    allocate(colTile0(nlatlon))
!    if(allocated(rowTile0))deallocate(rowTile0)
!    allocate(rowTile0(nlatlon))

! Allocate/Initialize pointer arrays if present

! global grid, not optional. no need to check if present
    if(present(igrdGlobal))then
       if(associated(igrdGlobal))deallocate(igrdGlobal)
       allocate(igrdGlobal(nlatlon))
       igrdGlobal=MISSING_REAL
    endif
    if(present(jgrdGlobal))then
       if(associated(jgrdGlobal))deallocate(jgrdGlobal)
       allocate(jgrdGlobal(nlatlon))
       jgrdGlobal=MISSING_REAL
    endif
    if(present(igrdLocal))then
       if(associated(igrdLocal))deallocate(igrdLocal)
       allocate(igrdLocal(nlatlon))
       igrdLocal=MISSING_REAL
    endif
    if(present(jgrdLocal))then
       if(associated(jgrdLocal))deallocate(jgrdLocal)
       allocate(jgrdLocal(nlatlon))
       jgrdLocal=MISSING_REAL
    endif
    if(present(colTile))then
       if(associated(colTile))deallocate(colTile)
       allocate(colTile(nlatlon))
       colTile=MISSING_INT
    endif
    if(present(rowTile))then
       if(associated(rowTile))deallocate(rowTile)
       allocate(rowTile(nlatlon))
       rowTile=MISSING_INT
    endif

! Loop over lat/lon locations calculating corresponding grid locations
! and tile file indices


! Make sure optional pointer arrays present if bilinear interpolation
! is set.
! If so, allocate memory

    if(bilintSinusoid .and. (.not. present(bilintParams) .or. .not. present(ifileDB_bilint)))then
       print *,'err[getSfcDBLocs]: ',&
            'bilintSinusoid=.True. with missing optional args bilintParams or ifileDB_bilint'
       call errorHalt(1)       
    endif

!    print *,'Begin loop over lat/lon locations'
    do i=1,nlatlon

! Get global col,row for current lat/lon       

! Skip if lat and/or lon missing

       if(lat(i) .eq. MISSING_REAL .or. lon(i) .eq. MISSING_REAL)cycle

!       print *,'call ll2ij_sinusoid'
!       print *,'i,lat(i),lon(i),igrdGlobal(i),jgrdGlobal(i),',&
!            'ncolGlobalGrid0,nrowGlobalGrid0,colGridOffsetGlobal0,rowGridOffsetGlobal0,',&
!            'sfcgrdDB_current%tile_scale(i),sfcgrdDB_current%tile_Re(i)=',&
!            i,lat(i),lon(i),igrdGlobal(i),jgrdGlobal(i),&
!            ncolGlobalGrid0,nrowGlobalGrid0,colGridOffsetGlobal0,rowGridOffsetGlobal0,&
!            sfcgrdDB_current%tile_scale(i),sfcgrdDB_current%tile_Re(i)

       call ll2ij_sinusoid(lat=lat(i), lon=lon(i), &
            col=igrdGlobal0(i), row=jgrdGlobal0(i), &
            nx=ncolGlobalGrid0, ny=nrowGlobalGrid0, &
            colGridOffset=colGridOffsetGlobal0, &
            rowGridOffset=rowGridOffsetGlobal0, &
            mapRes=scale0, &
            Re=Re0)
       if(present(igrdLocal))then
          igrdLocal(i)=mod(igrdGlobal0(i),real(sfcgrdDB_current%tile_ncol(1)))
       endif
       if(present(jgrdLocal))then
          jgrdLocal(i)=mod(jgrdGlobal0(i),real(sfcgrdDB_current%tile_nrow(1)))
       endif
       if(present(colTile))then
          colTile(i)=int(igrdGlobal0(i)/sfcgrdDB_current%tile_ncol(1))
       endif
       if(present(rowTile))then
          rowTile(i)=int(jgrdGlobal0(i)/sfcgrdDB_current%tile_nrow(1))
       endif

    enddo
!    print *,'End loop over lat/lon locations'
!    print *,'getSfcDBlocs:igrdGlobal0,jgrdGlobal0,igrdLocal,jgrdLocal,colTile,rowTile=',&
!         igrdGlobal0,jgrdGlobal0,igrdLocal,jgrdLocal,colTile,rowTile

! Find corresponding tile file indices, if present. First initialize to missing
    if(present(ifileDB))then
       if(associated(ifileDB))deallocate(ifileDB)
       allocate(ifileDB(nlatlon))
       ifileDB=MISSING_INT
       loop_nlatlon: do i=1,nlatlon
!          print *,'getSfcDBlocs: current lat, lon, computed tile col,row=',&
!               lat(i),lon(i),&
!               int(igrdGlobal0(i)/sfcgrdDB_current%tile_ncol(1)),&
!               int(jgrdGlobal0(i)/sfcgrdDB_current%tile_nrow(1))
          loop_ntilefiles: do j=1,sfcgrdDB_current%ntilefiles
             if(int(igrdGlobal0(i)/sfcgrdDB_current%tile_ncol(1)) .eq. sfcgrdDB_current%tile_colglobal(j) &
                  .and. &
                int(jgrdGlobal0(i)/sfcgrdDB_current%tile_nrow(1)) .eq. sfcgrdDB_current%tile_rowglobal(j))then
                ifileDB(i)=j
!                print *,'getSfcDBlocs: found matching tile file in DB: ifileDB(i),sfcgrdDB_current%tile_colglobal(j),sfcgrdDB_current%tile_rowglobal(j)=',&
!                     ifileDB(i),sfcgrdDB_current%tile_colglobal(j),sfcgrdDB_current%tile_rowglobal(j)
                cycle loop_nlatlon
             endif
          enddo loop_ntilefiles
!          print *,'warning[getSfcDBlocs]: no tile file index match found; index set to missing'
!          print *,'i,lat,lon,igrdGlobal0,jgrdGlobal0=',&
!               i,lat(i),lon(i),igrdGlobal0(i),jgrdGlobal0(i)
       enddo loop_nlatlon
    endif

! If bilinear interpolation option is on, get the corresponding parameters

    if(bilintSinusoid)then
       call getBilintParams(igrdGlobal=igrdGlobal0, &
            jgrdGlobal=jgrdGlobal0, &
            lonFOV=lon, &
            ncolGlobalGrid=ncolGlobalGrid0, &
            nrowGlobalGrid=nrowGlobalGrid0, &
            colGridOffsetGlobal=colGridOffsetGlobal0, &
            rowGridOffsetGlobal=rowGridOffsetGlobal0, &
            tile_ncol=sfcgrdDB_current%tile_ncol(1), &
            tile_nrow=sfcgrdDB_current%tile_nrow(1), &
            scale=scale0, Re=Re0, &
            bilintParams=bilintParams)

!       print *,'begin call getBilintFileInd'
       call getBilintFileInd(sfcgrdDB_current=sfcgrdDB_current,&
            bilintParams=bilintParams,&
            lonFOV=lon, &
            ifileDB_bilint=ifileDB_bilint)
!       print *,'end call getBilintFileInd'
!       print *,'size(ifileDB_bilint)=',size(ifileDB_bilint)
!       print *,'ifileDB_bilint(:,1)=',ifileDB_bilint(:,1)
!       print *,'ifileDB_bilint(:,2)=',ifileDB_bilint(:,2)
!       print *,'ifileDB_bilint(:,3)=',ifileDB_bilint(:,3)
!       print *,'ifileDB_bilint(:,4)=',ifileDB_bilint(:,4)
    endif

    if(present(ncolGlobalGrid))ncolGlobalGrid=ncolGlobalGrid0
    if(present(nrowGlobalGrid))nrowGlobalGrid=nrowGlobalGrid0
    if(present(colGridOffsetGlobal))colGridOffsetGlobal=colGridOffsetGlobal0
    if(present(rowGridOffsetGlobal))rowGridOffsetGlobal=rowGridOffsetGlobal0
    if(present(igrdGlobal))igrdGlobal=igrdGlobal0
    if(present(jgrdGlobal))jgrdGlobal=jgrdGlobal0

    if(allocated(igrdGlobal0))deallocate(igrdGlobal0)
    if(allocated(jgrdGlobal0))deallocate(jgrdGlobal0)
!    if(allocated(igrdLocal0))deallocate(igrdLocal0)
!    if(allocated(jgrdLocal0))deallocate(jgrdLocal0)
!    if(allocated(colTile0))deallocate(colTile0)
!    if(allocated(rowTile0))deallocate(rowTile0)

!    print *,'End getSfcDBLocs:'

  end subroutine getSfcDBlocs


!  subroutine getFileList(fnameDBFileList,FileList)
!    character (len=256) :: GaussWtFile




!  end subroutine getFileList







!================================================================================
!
! Read in input tile files, returning vectors of the tile grid parameters 
!
!================================================================================

  subroutine getTileParams(tilefile, &
       tile_col, tile_row, &
       nx_globaltiles, ny_globaltiles, &
       nx_tile, ny_tile, &
       colGridOffset, rowGridOffset, &
    latorigin, lonorigin, &
       Re, scale)
    
    character(len=*), dimension(:) :: tilefile
    integer, dimension(:), optional :: tile_col, tile_row
    integer, dimension(:), optional :: nx_tile, ny_tile
    integer, dimension(:), optional :: nx_globaltiles, ny_globaltiles
    real, dimension(:), optional :: colGridOffset, rowGridOffset
    real, dimension(:), optional :: latorigin, lonorigin
    real, dimension(:), optional :: Re, scale

! Local variables

    integer :: i, ntiles, ncid
    integer :: itile, jtile, ncol, nrow, nitile, njtile
    real :: ioffset, joffset, lat0, lon0, re0, scale0
    
    ntiles=size(tilefile)
    
! Loop over all files in list, reading tile and grid parameters
    
    do i=1,ntiles
       call openSfcGrd(ncid=ncid, file=tilefile(i), &
            ncol=ncol, nrow=nrow, &
            latorigin=lat0, lonorigin=lon0, &
            colGridOffset=ioffset, rowGridOffset=joffset, &
            scale=scale0, re=re0, &
            itile=itile, jtile=jtile, &
            nitile=nitile, njtile=njtile, &
            status='old')
       if(present(tile_col))tile_col(i)=itile
       if(present(tile_row))tile_row(i)=jtile
       if(present(nx_tile))nx_tile(i)=ncol
       if(present(ny_tile))ny_tile(i)=nrow
       if(present(nx_globaltiles))nx_globaltiles(i)=nitile
       if(present(ny_globaltiles))ny_globaltiles(i)=njtile
       if(present(colGridOffset))colGridOffset(i)=ioffset
       if(present(rowGridOffset))rowGridOffset(i)=joffset
       if(present(latorigin))latorigin(i)=lat0
       if(present(lonorigin))lonorigin(i)=lon0
       if(present(Re))Re(i)=re0
       if(present(scale))scale(i)=scale0
       call closeSfcGrd(ncid)
    enddo
    
    return
    
  end subroutine getTileParams

  subroutine getSfcLandType(sfcgrdDB_Land, lat, lon, resetMISSING, &
!  subroutine getSfcLandType(sfcgrdDB_Land, lat, lon, &
       itypeLand,dompcntLand)

    type(sfcgrdDB), intent(in) :: sfcgrdDB_Land
    real, dimension(:), intent(in) :: lat, lon
    logical, optional, intent(in) :: resetMISSING

    integer, dimension(:,:), pointer, optional :: itypeLand
    integer, dimension(:,:), pointer, optional :: dompcntLand
    
! local variables

    integer :: nlatlon, nValsPerGrid
    integer*1, allocatable, dimension(:) :: itypeLandV
    integer*2, allocatable, dimension(:) :: dompcntLandV
    logical :: resetMISSING0

    integer, dimension(:), pointer :: ifileDB_Land
    real, dimension(:), pointer :: igrdLocal_Land
    real, dimension(:), pointer :: jgrdLocal_Land

    integer :: ilocLand, jlocLand, i, ncidLand
    integer, dimension(1) :: imask

!

!    print *,'Start getSfcLandType'

    if(present(resetMISSING))then
       resetMISSING0=resetMISSING
    else
       resetMISSING0=.TRUE.
    endif

    nlatlon=size(lat)
    nValsPerGrid=2
!    print *,'nlatlon=',nlatlon

    if(allocated(itypeLandV))deallocate(itypeLandV)
    allocate(itypeLandV(nValsPerGrid))
    if(allocated(dompcntLandV))deallocate(dompcntLandV)
    allocate(dompcntLandV(nValsPerGrid))

    if(associated(ifileDB_Land))deallocate(ifileDB_Land)
    allocate(ifileDB_Land(nlatlon))
    if(associated(igrdLocal_Land))deallocate(igrdLocal_Land)
    allocate(igrdLocal_Land(nlatlon))
    if(associated(jgrdLocal_Land))deallocate(jgrdLocal_Land)
    allocate(jgrdLocal_Land(nlatlon))

    if(present(itypeLand))then
       if(associated(itypeLand))deallocate(itypeLand)
       allocate(itypeLand(nlatlon,nValsPerGrid))
       itypeLand=MISSING_INT
! Initialize to 0=water (so if tile file missing, assume water point)
       if(resetMISSING0)itypeLand(:,1)=0
    endif
    if(present(dompcntLand))then
       if(associated(dompcntLand))deallocate(dompcntLand)
       allocate(dompcntLand(nlatlon,nValsPerGrid))
       dompcntLand=MISSING_INT
! Initialize to 100 percent water (so if tile file missing, assume water point)
       if(resetMISSING0)dompcntLand(:,1)=10000
    endif

!    print *,'call getSfcDBlocs'
    call getSfcDBlocs(sfcgrdDB_current=sfcgrdDB_Land, &
         lat=lat(1:nlatlon), &
         lon=lon(1:nlatlon),&
         bilintSinusoid=.False., &
         ifileDB=ifileDB_Land, &
         igrdLocal=igrdLocal_Land, &
         jgrdLocal=jgrdLocal_Land)

!    print *,'getsfcLandType: lat=',lat
!    print *,'getsfcLandType: lon=',lon
!    print *,'getsfcLandType: ifileDB_Land=',ifileDB_Land
!    print *,'getsfcLandType: igrdLocal_Land=',igrdLocal_Land
!    print *,'getsfcLandType: jgrdLocal_Land=',jgrdLocal_Land
    
    do i=1,nlatlon

!       print *,'i=',i
       if(ifileDB_Land(i) .gt. 0) then
          
!          print *,'call openSfcGrd'
          call openSfcGrd(ncid=ncidLand, &
               file=sfcgrdDB_Land%tilefile(ifileDB_Land(i)), &
               status='old')
          ilocLand=min(int(igrdLocal_Land(i))+1,sfcgrdDB_Land%tile_ncol(1))
          jlocLand=min(int(jgrdLocal_Land(i))+1,sfcgrdDB_Land%tile_nrow(1))
!          print *,'ilocLand,jlocLand=',ilocLand,jlocLand

          call getSfcGrd(ncid=ncidLand, imask=imask, &
               iloc=ilocLand, jloc=jlocLand, kloc=1, &
               x3d1_i1=itypeLandV,x3d1_i2=dompcntLandV)

!          print *,'getsfcLandType: itypeLandV=',itypeLandV
!          print *,'getsfcLandType: dompcntLandV=',dompcntLandV

          if(any(itypeLandV .lt. 0 .or. itypeLandV .gt. 16))then
             if(resetMISSING0)then
                where(itypeLandV .lt. 0)
                   itypeLandV=0
                   dompcntLandV=10000
                endwhere
!                itypeLandV(1)=0
!                dompcntLandV(1)=10000
             endif
          endif
          
          if(present(itypeLand))itypeLand(i,:)=itypeLandV
          if(present(dompcntLand))dompcntLand(i,:)=dompcntLandV

          call closeSfcGrd(ncidLand)

       endif

    enddo

    deallocate(itypeLandV)
    deallocate(dompcntLandV)
    deallocate(ifileDB_Land)
    deallocate(igrdLocal_Land)
    deallocate(jgrdLocal_Land)


  end subroutine getSfcLandType


  subroutine getSfcWaterType(sfcgrdDB_Water, lat, lon, resetMISSING, bilintSinusoid, &
!  subroutine getSfcWaterType(sfcgrdDB_Water, lat, lon, &
       pcntWater,fracLand)

    type(sfcgrdDB), intent(in) :: sfcgrdDB_Water
    real, dimension(:), intent(in) :: lat, lon
    logical, intent(in) ::  bilintSinusoid
    logical, optional, intent(in) :: resetMISSING

    integer, dimension(:,:), pointer, optional :: pcntWater
    real, dimension(:), pointer, optional :: fracLand
    
! local variables

    integer :: nlatlon, nValsPerGrid
    integer*2, allocatable, dimension(:) :: pcntWaterV
    logical :: resetMISSING0
    logical :: useBilint

    integer, dimension(:), pointer :: ifileDB_Water
    real, dimension(:), pointer :: igrdLocal_Water
    real, dimension(:), pointer :: jgrdLocal_Water
    real, dimension(:), pointer :: jgrdGlobal_Water

    integer, dimension(:,:), pointer :: ifileDB_bilint
    real, allocatable, dimension(:) :: fracLand_i0jf, fracLand_i1jf
    real, allocatable, dimension(:) :: fracLand_i0jc, fracLand_i1jc
    type(bilintSin) :: bilintParams

    integer :: ilocWater, jlocWater, i, ncidWater
    real, allocatable, dimension(:) :: jgrdGlob, fracLand_bilint
    real, allocatable, dimension(:) :: i0floorG,i_jfloorG,i0ceilG,i_jceilG,jfloorG

    integer, dimension(1) :: imask

!
    nlatlon=size(lat)
    nValsPerGrid=4

    if(present(resetMISSING))then
       resetMISSING0=resetMISSING
    else
       resetMISSING0=.TRUE.
    endif

    useBilint=.FALSE.
    if(bilintSinusoid)then
       useBilint=.TRUE.
       if(allocated(fracLand_i0jf))deallocate(fracLand_i0jf)
       allocate(fracLand_i0jf(nlatlon))
       if(allocated(fracLand_i1jf))deallocate(fracLand_i1jf)
       allocate(fracLand_i1jf(nlatlon))
       if(allocated(fracLand_i0jc))deallocate(fracLand_i0jc)
       allocate(fracLand_i0jc(nlatlon))
       if(allocated(fracLand_i1jc))deallocate(fracLand_i1jc)
       allocate(fracLand_i1jc(nlatlon))
    endif


    if(allocated(pcntWaterV))deallocate(pcntWaterV)
    allocate(pcntWaterV(nValsPerGrid))

    if(associated(ifileDB_Water))deallocate(ifileDB_Water)
    allocate(ifileDB_Water(nlatlon))
    if(associated(igrdLocal_Water))deallocate(igrdLocal_Water)
    allocate(igrdLocal_Water(nlatlon))
    if(associated(jgrdLocal_Water))deallocate(jgrdLocal_Water)
    allocate(jgrdLocal_Water(nlatlon))

    if(present(pcntWater))then
       if(associated(pcntWater))deallocate(pcntWater)
       allocate(pcntWater(nlatlon,nValsPerGrid))
       pcntWater=MISSING_INT
    endif
    if(present(fracLand))then
       if(associated(fracLand))deallocate(fracLand)
       allocate(fracLand(nlatlon))
! Eventually Initialize to 0. (ocean point) instead?
       fracLand=MISSING_REAL
       if(resetMISSING0)fracLand=0.
    endif

! No bilinear interpolation requested
    if(.not. useBilint)then
       call getSfcDBlocs(sfcgrdDB_current=sfcgrdDB_Water, &
            lat=lat(1:nlatlon), &
            lon=lon(1:nlatlon),&
            bilintSinusoid=.False., &
            ifileDB=ifileDB_Water, &
            igrdLocal=igrdLocal_Water, &
            jgrdLocal=jgrdLocal_Water)

!    print *,'getsfcWaterType: lat=',lat
!    print *,'getsfcWaterType: lon=',lon
!    print *,'getsfcWaterType: ifileDB_Water=',ifileDB_Water
!    print *,'getsfcWaterType: igrdLocal_Water=',igrdLocal_Water
!    print *,'getsfcWaterType: jgrdLocal_Water=',jgrdLocal_Water
    
       do i=1,nlatlon
          
          if(ifileDB_Water(i) .gt. 0) then
          
             call openSfcGrd(ncid=ncidWater, &
                  file=sfcgrdDB_Water%tilefile(ifileDB_Water(i)), &
                  status='old')
             ilocWater=min(int(igrdLocal_Water(i))+1,sfcgrdDB_Water%tile_ncol(1))
             jlocWater=min(int(jgrdLocal_Water(i))+1,sfcgrdDB_Water%tile_nrow(1))

             call getSfcGrd(ncid=ncidWater, imask=imask, &
                  iloc=ilocWater, jloc=jlocWater, kloc=1, &
                  x3d1_i2=pcntWaterV)

!          print *,'getsfcWaterType: pcntWaterV=',pcntWaterV

             if(present(pcntWater))pcntWater(i,:)=pcntWaterV
             if(present(fracLand))then
!             print *,'any(pcntWaterV .lt. 0)=',any(pcntWaterV .lt. 0)
                if(any(pcntWaterV .lt. 0))then
!                print *,'getsfcWaterType: WARNING, element(s) of pcntWaterV .lt. 0'
                   fracLand(i)= MISSING_REAL
                   if(resetMISSING0)fracLand(i)=0.
                else
                   fracLand(i)= 1.-real(sum(pcntWaterV))/10000.
                endif
             endif
!          print *,'getsfcWaterType: fracLand=',fracLand

             call closeSfcGrd(ncidWater)
             
          endif

       enddo

! Use bilinear interpolation rather than nearest neighbor for land fraction
    else
!       print *,'getSfcWaterType: bilint option'
!       print *,'call getSfcDBlocs'
       call getSfcDBlocs(sfcgrdDB_current=sfcgrdDB_Water, &
            bilintSinusoid=.True.,&
            lat=lat(1:nlatlon), &
            lon=lon(1:nlatlon),&
            ifileDB=ifileDB_Water, &
            ifileDB_bilint=ifileDB_bilint, &
            igrdLocal=igrdLocal_Water, &
            jgrdLocal=jgrdLocal_Water, &
            jgrdGlobal=jgrdGlobal_Water, &
            bilintParams=bilintParams)

!       print *,'finished call getSfcDBlocs'
!       print *,'allocate mirror arrays'
       if(allocated(jgrdGlob))deallocate(jgrdGlob)
       allocate(jgrdGlob(nlatlon))
       if(allocated(i0floorG))deallocate(i0floorG)
       allocate(i0floorG(nlatlon))
       if(allocated(i_jfloorG))deallocate(i_jfloorG)
       allocate(i_jfloorG(nlatlon))
       if(allocated(i0ceilG))deallocate(i0ceilG)
       allocate(i0ceilG(nlatlon))
       if(allocated(i_jceilG))deallocate(i_jceilG)
       allocate(i_jceilG(nlatlon))
       if(allocated(jfloorG))deallocate(jfloorG)
       allocate(jfloorG(nlatlon))
       if(allocated(fracLand_bilint))deallocate(fracLand_bilint)
       allocate(fracLand_bilint(nlatlon))

!       print *,'copy into mirror arrays'
       jgrdGlob=jgrdGlobal_Water
       i0floorG=bilintParams%i0floorG
       i_jfloorG=bilintParams%i_jfloorG
       i0ceilG=bilintParams%i0ceilG
       i_jceilG=bilintParams%i_jceilG
       jfloorG=bilintParams%jfloorG

       do i=1,nlatlon

!          print *,'i=',i
! Even if bilinear interp is requested for land fraction, still use nearest neighbor for pcntWater
          if(present(pcntWater) .and. ifileDB_Water(i) .gt. 0)then
!             print *,'pcntWater present: using NN interp'
             call openSfcGrd(ncid=ncidWater, &
                  file=sfcgrdDB_Water%tilefile(ifileDB_Water(i)), &
                  status='old')
             ilocWater=min(int(igrdLocal_Water(i))+1,sfcgrdDB_Water%tile_ncol(1))
             jlocWater=min(int(jgrdLocal_Water(i))+1,sfcgrdDB_Water%tile_nrow(1))

             call getSfcGrd(ncid=ncidWater, imask=imask, &
                  iloc=ilocWater, jloc=jlocWater, kloc=1, &
                  x3d1_i2=pcntWaterV)
             pcntWater(i,:)=pcntWaterV
             call closeSfcGrd(ncidWater)
          endif
          if(any(ifileDB_bilint(i,:) .le. 0))cycle
! Special case: one or more params are missing. 
! OR
! Special case: within one-half grid unit from top or bottom of grid. 
! Use nearest neighbor for land fraction instead
          if(BilintParams%jfloorALocal(i) .eq. MISSING_INT .or. &
             BilintParams%jceilALocal(i) .eq. MISSING_INT .or. &
             BilintParams%i0floorALocal(i) .eq. MISSING_INT .or. &
             BilintParams%i1floorALocal(i) .eq. MISSING_INT .or. &
             BilintParams%i0ceilALocal(i) .eq. MISSING_INT .or. &
             BilintParams%i1ceilALocal(i) .eq. MISSING_INT .or. &
             (BilintParams%jfloorALocal(i) .eq. BilintParams%jceilALocal(i)))then
!             print *,'special case: use NN interp'
! If pcntWater has already been read in
             if(ifileDB_Water(i) .gt. 0)then
                if(present(pcntWater))then
                   if(any(pcntWaterV .lt. 0.))cycle
                   fracLand(i)= 1.-real(sum(pcntWaterV))/10000.
                else
                   call openSfcGrd(ncid=ncidWater, &
                        file=sfcgrdDB_Water%tilefile(ifileDB_Water(i)), &
                        status='old')
                   ilocWater=min(int(igrdLocal_Water(i))+1,sfcgrdDB_Water%tile_ncol(1))
                   jlocWater=min(int(jgrdLocal_Water(i))+1,sfcgrdDB_Water%tile_nrow(1))

                   call getSfcGrd(ncid=ncidWater, imask=imask, &
                        iloc=ilocWater, jloc=jlocWater, kloc=1, &
                        x3d1_i2=pcntWaterV)
                   call closeSfcGrd(ncidWater)
                   if(any(pcntWaterV .lt. 0.))cycle
                   fracLand(i)= 1.-real(sum(pcntWaterV))/10000.
                endif

             endif

! Normal case: all params valid. do bilinear interpolation
          else
!             print *,'Normal case: bilinear interp'

! Get the data from the grid at the 4 specified locations

!  1: (i0floorALOcal,jfloorAlocal)
!             print *,'1: (i0floorALOcal,jfloorAlocal)'
             call openSfcGrd(ncid=ncidWater, &
                  file=sfcgrdDB_Water%tilefile(ifileDB_bilint(i,1)), &
                  status='old')
             ilocWater=BilintParams%i0floorALocal(i)
             jlocWater=BilintParams%jfloorALocal(i)

             call getSfcGrd(ncid=ncidWater, imask=imask, &
                  iloc=ilocWater, jloc=jlocWater, kloc=1, &
                  x3d1_i2=pcntWaterV)
             if(any(pcntWaterV .lt. 0))then
!                print *,'getsfcWaterType: WARNING, element(s) of pcntWaterV .lt. 0'
                fracLand_i0jf(i)= MISSING_REAL
                if(resetMISSING0)fracLand_i0jf(i)=0.
             else
                fracLand_i0jf(i)=1.-real(sum(pcntWaterV))/10000.
             endif
             call closeSfcGrd(ncidWater)

!  2: (i1floorALOcal,jfloorAlocal)
!             print *,'2: (i1floorALOcal,jfloorAlocal)'
             call openSfcGrd(ncid=ncidWater, &
                  file=sfcgrdDB_Water%tilefile(ifileDB_bilint(i,2)), &
                  status='old')
             ilocWater=BilintParams%i1floorALocal(i)
             jlocWater=BilintParams%jfloorALocal(i)

             call getSfcGrd(ncid=ncidWater, imask=imask, &
                  iloc=ilocWater, jloc=jlocWater, kloc=1, &
                  x3d1_i2=pcntWaterV)
             if(any(pcntWaterV .lt. 0))then
!                print *,'getsfcWaterType: WARNING, element(s) of pcntWaterV .lt. 0'
                fracLand_i1jf(i)= MISSING_REAL
                if(resetMISSING0)fracLand_i1jf(i)=0.
             else
                fracLand_i1jf(i)=1.-real(sum(pcntWaterV))/10000.
             endif
             call closeSfcGrd(ncidWater)

!  3: (i0ceilALOcal,jceilAlocal)
!             print *,'3: (i0ceilALOcal,jceilAlocal)'
             call openSfcGrd(ncid=ncidWater, &
                  file=sfcgrdDB_Water%tilefile(ifileDB_bilint(i,3)), &
                  status='old')
             ilocWater=BilintParams%i0ceilALocal(i)
             jlocWater=BilintParams%jceilALocal(i)

             call getSfcGrd(ncid=ncidWater, imask=imask, &
                  iloc=ilocWater, jloc=jlocWater, kloc=1, &
                  x3d1_i2=pcntWaterV)
             if(any(pcntWaterV .lt. 0))then
!                print *,'getsfcWaterType: WARNING, element(s) of pcntWaterV .lt. 0'
                fracLand_i0jc(i)= MISSING_REAL
                if(resetMISSING0)fracLand_i0jc(i)=0.
             else
                fracLand_i0jc(i)=1.-real(sum(pcntWaterV))/10000.
             endif
             call closeSfcGrd(ncidWater)

!  4: (i1ceilALOcal,jceilAlocal)
!             print *,'4: (i1ceilALOcal,jceilAlocal)'
             call openSfcGrd(ncid=ncidWater, &
                  file=sfcgrdDB_Water%tilefile(ifileDB_bilint(i,4)), &
                  status='old')
             ilocWater=BilintParams%i1ceilALocal(i)
             jlocWater=BilintParams%jceilALocal(i)

             call getSfcGrd(ncid=ncidWater, imask=imask, &
                  iloc=ilocWater, jloc=jlocWater, kloc=1, &
                  x3d1_i2=pcntWaterV)
             if(any(pcntWaterV .lt. 0))then
!                print *,'getsfcWaterType: WARNING, element(s) of pcntWaterV .lt. 0'
                fracLand_i1jc(i)= MISSING_REAL
                if(resetMISSING0)fracLand_i1jc(i)=0.
             else
                fracLand_i1jc(i)=1.-real(sum(pcntWaterV))/10000.
             endif
             call closeSfcGrd(ncidWater)

! Bilinearly interpolate on the sinusoid grid

!             print *,'call calcBilint'
             call calcBilint(data2d_i0jf=fracLand_i0jf(i:i),&
                  data2d_i1jf=fracLand_i1jf(i:i),&
                  data2d_i0jc=fracLand_i0jc(i:i),&
                  data2d_i1jc=fracLand_i1jc(i:i),&
                  i0floorG=i0floorG(i:i),&
                  i_jfloorG=i_jfloorG(i:i),&
                  i0ceilG=i0ceilG(i:i),&
                  i_jceilG=i_jceilG(i:i),&
                  jfloorG=jfloorG(i:i),&
                  jgrd=jgrdGlob(i:i),&
                  data2d_bilint=fracLand_bilint(i:i))
!             print *,'finished call calcBilint'
             fracLand(i)=fracLand_bilint(i)
!             print *,'finished calc for this i=',i
          endif
          
       enddo


    endif

    deallocate(pcntWaterV)
    deallocate(igrdLocal_Water)
    deallocate(jgrdLocal_Water)
    if(associated(ifileDB_Water))deallocate(ifileDB_Water)
    if(associated(ifileDB_bilint))deallocate(ifileDB_bilint)

    if(allocated(jgrdGlob))deallocate(jgrdGlob)
    if(allocated(i0floorG))deallocate(i0floorG)
    if(allocated(i_jfloorG))deallocate(i_jfloorG)
    if(allocated(i0ceilG))deallocate(i0ceilG)
    if(allocated(i_jceilG))deallocate(i_jceilG)
    if(allocated(jfloorG))deallocate(jfloorG)
    if(allocated(fracLand_bilint))deallocate(fracLand_bilint)

  end subroutine getSfcWaterType


  subroutine getSfcElev(sfcgrdDB_Elev, lat, lon, bilintSinusoid, &
       sfcElev)

    type(sfcgrdDB), intent(in) :: sfcgrdDB_Elev
    real, dimension(:), intent(in) :: lat, lon
    logical, intent(in) ::  bilintSinusoid

    integer, dimension(:), pointer :: sfcElev
    
! local variables

    integer :: nlatlon, nValsPerGrid
    integer*2, dimension(1) :: sfcElevV
    logical :: useBilint

    integer, dimension(:), pointer :: ifileDB_Elev
    real, dimension(:), pointer :: igrdLocal_Elev
    real, dimension(:), pointer :: jgrdLocal_Elev

    integer, dimension(:,:), pointer :: ifileDB_bilint
    real, allocatable, dimension(:) :: sfcElev_i0jf, sfcElev_i1jf
    real, allocatable, dimension(:) :: sfcElev_i0jc, sfcElev_i1jc
    type(bilintSin) :: bilintParams
    real, allocatable, dimension(:) :: jgrdGlob, sfcElev_bilint
    real, allocatable, dimension(:) :: i0floorG,i_jfloorG,i0ceilG,i_jceilG,jfloorG


    integer :: ilocElev, jlocElev, i, ncidElev
    integer, dimension(1) :: imask
    integer :: ilocMod

!

!    print *,'Start getSfcElev'

    nlatlon=size(lat)
!    print *,'nlatlon=',nlatlon

    useBilint=.FALSE.
    if(bilintSinusoid)then
       useBilint=.TRUE.
       if(allocated(sfcElev_i0jf))deallocate(sfcElev_i0jf)
       allocate(sfcElev_i0jf(nlatlon))
       if(allocated(sfcElev_i1jf))deallocate(sfcElev_i1jf)
       allocate(sfcElev_i1jf(nlatlon))
       if(allocated(sfcElev_i0jc))deallocate(sfcElev_i0jc)
       allocate(sfcElev_i0jc(nlatlon))
       if(allocated(sfcElev_i1jc))deallocate(sfcElev_i1jc)
       allocate(sfcElev_i1jc(nlatlon))
    endif

    if(associated(sfcElev))deallocate(sfcElev)
    allocate(sfcElev(nlatlon))
    if(associated(ifileDB_Elev))deallocate(ifileDB_Elev)
    allocate(ifileDB_Elev(nlatlon))
    if(associated(igrdLocal_Elev))deallocate(igrdLocal_Elev)
    allocate(igrdLocal_Elev(nlatlon))
    if(associated(jgrdLocal_Elev))deallocate(jgrdLocal_Elev)
    allocate(jgrdLocal_Elev(nlatlon))

! Initialize to missing
    sfcElev=MISSING_REAL

!    print *,'call getSfcDBlocs'
    if(.not. useBilint)then
    call getSfcDBlocs(sfcgrdDB_current=sfcgrdDB_Elev, &
         lat=lat(1:nlatlon), &
         lon=lon(1:nlatlon),&
         bilintSinusoid=.False., &
         ifileDB=ifileDB_Elev, &
         igrdLocal=igrdLocal_Elev, &
         jgrdLocal=jgrdLocal_Elev)

!    print *,'getsfcElev: lat=',lat
!    print *,'getsfcElev: lon=',lon
!    print *,'getsfcElev: ifileDB_Elev=',ifileDB_Elev
!    print *,'getsfcElev: igrdLocal_Elev=',igrdLocal_Elev
!    print *,'getsfcElev: jgrdLocal_Elev=',jgrdLocal_Elev
    
    do i=1,nlatlon

!       print *,'i=',i
       if(ifileDB_Elev(i) .gt. 0) then
          
!          print *,'call openSfcGrd'
          call openSfcGrd(ncid=ncidElev, &
               file=sfcgrdDB_Elev%tilefile(ifileDB_Elev(i)), &
               status='old')
          ilocElev=min(int(igrdLocal_Elev(i))+1,sfcgrdDB_Elev%tile_ncol(1))
          jlocElev=min(int(jgrdLocal_Elev(i))+1,sfcgrdDB_Elev%tile_nrow(1))
!          print *,'ilocElev,jlocElev=',ilocElev,jlocElev

          call getSfcGrd(ncid=ncidElev, imask=imask, &
               iloc=ilocElev, jloc=jlocElev, kloc=1, &
               x2d1_i2=sfcElevV)

! Check for off-hemisphere missing, due to roundoff
! Shift one-grid onto hemisphere, and retry
          if(sfcElevV(1) .lt. -900.)then
             if(lon(i) .lt. 0.)then
                ilocMod=min(ilocElev+1,sfcgrdDB_Elev%tile_ncol(1))
                call getSfcGrd(ncid=ncidElev, imask=imask, &
                     iloc=ilocMod, jloc=jlocElev, kloc=1, &
                     x2d1_i2=sfcElevV)
             elseif(lon(i) .gt. 0.)then
                ilocMod=max(ilocElev-1,1)
                call getSfcGrd(ncid=ncidElev, imask=imask, &
                     iloc=ilocMod, jloc=jlocElev, kloc=1, &
                     x2d1_i2=sfcElevV)
             endif
          endif
!          print *,'getsfcElev: sfcElevV=',sfcElevV

          sfcElev(i)=sfcElevV(1)
!          print *,'finished: sfcElev(i)=sfcElevV(1)'
          call closeSfcGrd(ncidElev)
!          print *,'finished:closeSfcGrd(ncidElev)'

       endif

    enddo

! Use bilinear interpolation rather than nearest neighbor for surface elevation
    else
!       print *,'getSfcElev: bilint option'
!       print *,'call getSfcDBlocs'
       print *,'err[getSfcElev]: bilint option not implemented'
       call errorHalt(1)       
    endif


    deallocate(igrdLocal_Elev)
    deallocate(jgrdLocal_Elev)
    if(associated(ifileDB_Elev))deallocate(ifileDB_Elev)
    if(associated(ifileDB_bilint))deallocate(ifileDB_bilint)


  end subroutine getSfcElev

  subroutine getSfcSeaIce(sfcgrdDB_SeaIce, lat, lon, &
       sfcSeaIce)

    type(sfcgrdDB), intent(in) :: sfcgrdDB_SeaIce
    real, dimension(:), intent(in) :: lat, lon

    integer, dimension(:), pointer :: sfcSeaIce
    
! local variables

    integer :: nlatlon, nValsPerGrid
    integer*1, dimension(1) :: sfcSeaIceV

    integer, dimension(:), pointer :: ifileDB_SeaIce
    real, dimension(:), pointer :: igrdLocal_SeaIce
    real, dimension(:), pointer :: jgrdLocal_SeaIce

    integer :: ilocSeaIce, jlocSeaIce, i, ncidSeaIce
    integer, dimension(1) :: imask
    integer :: ilocMod

!

!    print *,'Start getSfcSeaIce'

    nlatlon=size(lat)
!    print *,'nlatlon=',nlatlon

    if(associated(sfcSeaIce))deallocate(sfcSeaIce)
    allocate(sfcSeaIce(nlatlon))
    if(associated(ifileDB_SeaIce))deallocate(ifileDB_SeaIce)
    allocate(ifileDB_SeaIce(nlatlon))
    if(associated(igrdLocal_SeaIce))deallocate(igrdLocal_SeaIce)
    allocate(igrdLocal_SeaIce(nlatlon))
    if(associated(jgrdLocal_SeaIce))deallocate(jgrdLocal_SeaIce)
    allocate(jgrdLocal_SeaIce(nlatlon))

! Initialize to missing
    sfcSeaIce=MISSING_INT

!    print *,'call getSfcDBlocs'
    call getSfcDBlocs(sfcgrdDB_current=sfcgrdDB_SeaIce, &
         lat=lat(1:nlatlon), &
         lon=lon(1:nlatlon),&
         bilintSinusoid=.False., &
         ifileDB=ifileDB_SeaIce, &
         igrdLocal=igrdLocal_SeaIce, &
         jgrdLocal=jgrdLocal_SeaIce)

!    print *,'getsfcSeaIce: lat=',lat
!    print *,'getsfcSeaIce: lon=',lon
!    print *,'getsfcSeaIce: ifileDB_SeaIce=',ifileDB_SeaIce
!    print *,'getsfcSeaIce: igrdLocal_SeaIce=',igrdLocal_SeaIce
!    print *,'getsfcSeaIce: jgrdLocal_SeaIce=',jgrdLocal_SeaIce
    
    do i=1,nlatlon

!       print *,'i=',i
       if(ifileDB_SeaIce(i) .gt. 0) then
          
!          print *,'call openSfcGrd'
          call openSfcGrd(ncid=ncidSeaIce, &
               file=sfcgrdDB_SeaIce%tilefile(ifileDB_SeaIce(i)), &
               status='old')
          ilocSeaIce=min(int(igrdLocal_SeaIce(i))+1,sfcgrdDB_SeaIce%tile_ncol(1))
          jlocSeaIce=min(int(jgrdLocal_SeaIce(i))+1,sfcgrdDB_SeaIce%tile_nrow(1))
!          print *,'ilocSeaIce,jlocSeaIce=',ilocSeaIce,jlocSeaIce

          call getSfcGrd(ncid=ncidSeaIce, imask=imask, &
               iloc=ilocSeaIce, jloc=jlocSeaIce, kloc=1, &
               x2d1_i1=sfcSeaIceV)

! Check for off-hemisphere missing, due to roundoff
! Shift one-grid onto hemisphere, and retry
!          if(sfcSeaIceV(1) .gt. 100.)then
!             if(lon(i) .lt. 0.)then
!                ilocMod=min(ilocSeaIce+1,sfcgrdDB_SeaIce%tile_ncol(1))
!                call getSfcGrd(ncid=ncidSeaIce, imask=imask, &
!                     iloc=ilocMod, jloc=jlocSeaIce, kloc=1, &
!                     x2d1_i1=sfcSeaIceV)
!             elseif(lon(i) .gt. 0.)then
!                ilocMod=max(ilocSeaIce-1,1)
!                call getSfcGrd(ncid=ncidSeaIce, imask=imask, &
!                     iloc=ilocMod, jloc=jlocSeaIce, kloc=1, &
!                     x2d1_i1=sfcSeaIceV)
!             endif
!          endif
!          print *,'getsfcSeaIce: sfcSeaIceV=',sfcSeaIceV

          if(sfcSeaIceV(1) .ge. 0 .and. sfcSeaIceV(1) .le. 100)sfcSeaIce(i)=sfcSeaIceV(1)
!          print *,'finished: sfcSeaIce(i)=sfcSeaIceV(1)'
          call closeSfcGrd(ncidSeaIce)
!          print *,'finished:closeSfcGrd(ncidSeaIce)'

       endif

    enddo

    deallocate(igrdLocal_SeaIce)
    deallocate(jgrdLocal_SeaIce)
    if(associated(ifileDB_SeaIce))deallocate(ifileDB_SeaIce)


  end subroutine getSfcSeaIce


  subroutine getBilintParams(igrdGlobal, jgrdGlobal, lonFOV, &
       ncolGlobalGrid, nrowGlobalGrid, &
       colGridOffsetGlobal, rowGridOffsetGlobal, &
       tile_ncol, tile_nrow, &
       scale, Re, &
       bilintParams)
    
    real, dimension(:), intent(in) :: igrdGlobal, jgrdGlobal
    real, dimension(:), intent(in) :: lonFOV
    integer, intent(in) :: ncolGlobalGrid, nrowGlobalGrid
    real, intent(in) :: colGridOffsetGlobal, rowGridOffsetGlobal
    integer, intent(in) :: tile_ncol, tile_nrow
    real, intent(in) :: scale, Re

    type(bilintSin), intent(out) :: bilintParams

! needed for bilinear interpolation option
! two sets of arrays needed due to: 
! (1) shift of grid values w/r to array locations (i.e. value stored at array(1,1) is valid at grid(0.5,0.5))
! (2) possibility of wraparound at dateline. Required array locations will not be consistent with grid values
! vars for storing *grid* locations (non-integer values, needed for finite differences)
! vars for storing *array* index locations in local tiles (from tile files)

! local variables

    integer :: nfov,i
    real :: i180E, i180W, di180, i1floorG, i1ceilG, igrdWrap, j180
    real :: latGrid, lonGrid
    integer :: i0floorA, i1floorA
    integer :: i0ceilA, i1ceilA
    integer :: jfloorA, jceilA

! allocate pointer arrays

    nfov=size(igrdGlobal)

    if(associated(bilintParams%i0floorG))deallocate(bilintParams%i0floorG)
    allocate(bilintParams%i0floorG(nfov))
    bilintParams%i0floorG=MISSING_REAL
    if(associated(bilintParams%i0ceilG))deallocate(bilintParams%i0ceilG)
    allocate(bilintParams%i0ceilG(nfov))
    bilintParams%i0ceilG=MISSING_REAL
    if(associated(bilintParams%jfloorG))deallocate(bilintParams%jfloorG)
    allocate(bilintParams%jfloorG(nfov))
    bilintParams%jfloorG=MISSING_REAL
    if(associated(bilintParams%jceilG))deallocate(bilintParams%jceilG)
    allocate(bilintParams%jceilG(nfov))
    bilintParams%jceilG=MISSING_REAL
    if(associated(bilintParams%i_jfloorG))deallocate(bilintParams%i_jfloorG)
    allocate(bilintParams%i_jfloorG(nfov))
    bilintParams%i_jfloorG=MISSING_REAL
    if(associated(bilintParams%i_jceilG))deallocate(bilintParams%i_jceilG)
    allocate(bilintParams%i_jceilG(nfov))
    bilintParams%i_jceilG=MISSING_REAL
    if(associated(bilintParams%i0floorALocal))deallocate(bilintParams%i0floorALocal)
    allocate(bilintParams%i0floorALocal(nfov))
    bilintParams%i0floorALocal=MISSING_INT
    if(associated(bilintParams%i1floorALocal))deallocate(bilintParams%i1floorALocal)
    allocate(bilintParams%i1floorALocal(nfov))
    bilintParams%i1floorALocal=MISSING_INT
    if(associated(bilintParams%i0ceilALocal))deallocate(bilintParams%i0ceilALocal)
    allocate(bilintParams%i0ceilALocal(nfov))
    bilintParams%i0ceilALocal=MISSING_INT
    if(associated(bilintParams%i1ceilALocal))deallocate(bilintParams%i1ceilALocal)
    allocate(bilintParams%i1ceilALocal(nfov))
    bilintParams%i1ceilALocal=MISSING_INT
    if(associated(bilintParams%jfloorALocal))deallocate(bilintParams%jfloorALocal)
    allocate(bilintParams%jfloorALocal(nfov))
    bilintParams%jfloorALocal=MISSING_INT
    if(associated(bilintParams%jceilALocal))deallocate(bilintParams%jceilALocal)
    allocate(bilintParams%jceilALocal(nfov))
    bilintParams%jceilALocal=MISSING_INT
    if(associated(bilintParams%i0floor_tilecol))deallocate(bilintParams%i0floor_tilecol)
    allocate(bilintParams%i0floor_tilecol(nfov))
    bilintParams%i0floor_tilecol=MISSING_INT
    if(associated(bilintParams%i1floor_tilecol))deallocate(bilintParams%i1floor_tilecol)
    allocate(bilintParams%i1floor_tilecol(nfov))
    bilintParams%i1floor_tilecol=MISSING_INT
    if(associated(bilintParams%i0ceil_tilecol))deallocate(bilintParams%i0ceil_tilecol)
    allocate(bilintParams%i0ceil_tilecol(nfov))
    bilintParams%i0ceil_tilecol=MISSING_INT
    if(associated(bilintParams%i1ceil_tilecol))deallocate(bilintParams%i1ceil_tilecol)
    allocate(bilintParams%i1ceil_tilecol(nfov))
    bilintParams%i1ceil_tilecol=MISSING_INT
    if(associated(bilintParams%jfloor_tilerow))deallocate(bilintParams%jfloor_tilerow)
    allocate(bilintParams%jfloor_tilerow(nfov))
    bilintParams%jfloor_tilerow=MISSING_INT
    if(associated(bilintParams%jceil_tilerow))deallocate(bilintParams%jceil_tilerow)
    allocate(bilintParams%jceil_tilerow(nfov))
    bilintParams%jceil_tilerow=MISSING_INT

! Loop over all input grids calculating interpolation parameters one element at a time

    do i=1,nfov

       if(lonFOV(i) .eq. MISSING_REAL)cycle
! Grid coordinates of the rows directly neighboring the current point
! Factors of 0.5 are needed because array(1,1) represents data centered at grid(0.5,0.5)
! max,min used to make sure we stay on the grid
! Later on we'll need to check for this special case when jfloorG and jceilG are the same
       bilintParams%jfloorG(i)=max(real(int(jgrdGlobal(i)+0.5))-0.5,0.5)
       bilintParams%jceilG(i)=min(real(int(jgrdGlobal(i)+0.5))+0.5,nrowGlobalGrid-0.5)

! Corresponding global array row indices 
       jfloorA=int(bilintParams%jfloorG(i))+1
       jceilA=int(bilintParams%jceilG(i))+1
! Converted to local tile indices for output
       bilintParams%jfloorALocal(i)=mod((jfloorA-1),tile_nrow)+1
       bilintParams%jceilALocal(i)=mod((jceilA-1),tile_nrow)+1
! Corresponding tile row indices
       bilintParams%jfloor_tilerow(i)=int((jfloorA-1)/tile_nrow)
       bilintParams%jceil_tilerow(i)=int((jceilA-1)/tile_nrow)

! Column grid coordinate for lonFOV in upper neighboring row
       call jlon2i_sinusoid(row=bilintParams%jfloorG(i), &
            lon=lonFOV(i), col=bilintParams%i_jfloorG(i), &
            nx=ncolGlobalGrid, ny=nrowGlobalGrid, &
            colGridOffset=colGridOffsetGlobal, &
            rowGridOffset=rowGridOffsetGlobal, &
            mapRes=scale, Re=Re)

! Column grid coordinate for lonFOV in lower neighboring row
       call jlon2i_sinusoid(row=bilintParams%jceilG(i), &
            lon=lonFOV(i), col=bilintParams%i_jceilG(i), &
            nx=ncolGlobalGrid, ny=nrowGlobalGrid, &
            colGridOffset=colGridOffsetGlobal, &
            rowGridOffset=rowGridOffsetGlobal, &
            mapRes=scale, Re=Re)

! Neighboring column coordinates surrounding i_jfloorG and i_jceilG
       bilintParams%i0floorG(i)=real(int(bilintParams%i_jfloorG(i)+0.5))-0.5
       i1floorG=real(int(bilintParams%i_jfloorG(i)+0.5))+0.5
       bilintParams%i0ceilG(i)=real(int(bilintParams%i_jceilG(i)+0.5))-0.5
       i1ceilG=real(int(bilintParams%i_jceilG(i)+0.5))+0.5
!       print *,'bilintParams%i0floorG(i),i1floorG=',bilintParams%i0floorG(i),i1floorG
!       print *,'bilintParams%i0ceilG(i),i1ceilG=',bilintParams%i0ceilG(i),i1ceilG

! Check for any column grid coordinates that are off-hemisphere, and 
! try to wrap around to other side of grid/dateline. Otherwise,
! convert directly to corresponding array indices

! Checking i0floorG:
       call ij2ll_sinusoid(col=bilintParams%i0floorG(i), row=bilintParams%jfloorG(i),&
            lat=latGrid, lon=lonGrid, &
            nx=ncolGlobalGrid, ny=nrowGlobalGrid, &
            colGridOffset=colGridOffsetGlobal, &
            rowGridOffset=rowGridOffsetGlobal, &
            mapRes=scale, Re=Re)
       
       if(lonGrid .le. -900.)then
          if(lonFOV(i) .gt. 0. .and. lonFOV(i) .le. 180.)then
             print *,'err[getBilintParams]: i0floorG is off-hemisphere on right side of grid'
             call errorHalt(1)       
          endif
          call ll2ij_sinusoid(lat=latGrid, lon=-180., col=i180W, row=j180, &
               nx=ncolGlobalGrid, ny=nrowGlobalGrid, &
               colGridOffset=colGridOffsetGlobal, &
               rowGridOffset=rowGridOffsetGlobal, &
               mapRes=scale, Re=Re)
          di180=bilintParams%i0floorG(i)-i180W
          call ll2ij_sinusoid(lat=latGrid, lon=180., col=i180E, row=j180, &
               nx=ncolGlobalGrid, ny=nrowGlobalGrid, &
               colGridOffset=colGridOffsetGlobal, &
               rowGridOffset=rowGridOffsetGlobal, &
               mapRes=scale, Re=Re)
          igrdWrap=i180E+di180
! Check if wrapped grid is on-hemisphere
          call ij2ll_sinusoid(col=igrdWrap, row=bilintParams%jfloorG(i),&
               lat=latGrid, lon=lonGrid, &
               nx=ncolGlobalGrid, ny=nrowGlobalGrid, &
               colGridOffset=colGridOffsetGlobal, &
               rowGridOffset=rowGridOffsetGlobal, &
               mapRes=scale, Re=Re)
          if(lonGrid .le. -900.)then
             i0floorA=MISSING_INT
             bilintParams%i0floorALocal(i)=MISSING_INT
          else
             i0floorA=int(igrdWrap)+1
             bilintParams%i0floorALocal(i)=mod((i0floorA-1),tile_ncol)+1
          endif

       else
          i0floorA=int(bilintParams%i0floorG(i))+1
          bilintParams%i0floorALocal(i)=mod((i0floorA-1),tile_ncol)+1
       endif
! Corresponding tile column index
       if(i0floorA .ne. MISSING_INT)then
          bilintParams%i0floor_tilecol(i)=int((i0floorA-1)/tile_ncol)
       endif

! Checking i0ceilG:
       call ij2ll_sinusoid(col=bilintParams%i0ceilG(i), row=bilintParams%jceilG(i),&
            lat=latGrid, lon=lonGrid, &
            nx=ncolGlobalGrid, ny=nrowGlobalGrid, &
            colGridOffset=colGridOffsetGlobal, &
            rowGridOffset=rowGridOffsetGlobal, &
            mapRes=scale, Re=Re)

       if(lonGrid .le. -900.)then
          if(lonFOV(i) .gt. 0. .and. lonFOV(i) .le. 180.)then
             print *,'err[getBilintParams]: i0ceilG is off-hemisphere on right side of grid'
             call errorHalt(1)       
          endif
          call ll2ij_sinusoid(lat=latGrid, lon=-180., col=i180W, row=j180, &
               nx=ncolGlobalGrid, ny=nrowGlobalGrid, &
               colGridOffset=colGridOffsetGlobal, &
               rowGridOffset=rowGridOffsetGlobal, &
               mapRes=scale, Re=Re)
          di180=bilintParams%i0ceilG(i)-i180W
          call ll2ij_sinusoid(lat=latGrid, lon=180., col=i180E, row=j180, &
               nx=ncolGlobalGrid, ny=nrowGlobalGrid, &
               colGridOffset=colGridOffsetGlobal, &
               rowGridOffset=rowGridOffsetGlobal, &
               mapRes=scale, Re=Re)
          igrdWrap=i180E+di180
! Check if wrapped grid is on-hemisphere
          call ij2ll_sinusoid(col=igrdWrap, row=bilintParams%jceilG(i),&
               lat=latGrid, lon=lonGrid, &
               nx=ncolGlobalGrid, ny=nrowGlobalGrid, &
               colGridOffset=colGridOffsetGlobal, &
               rowGridOffset=rowGridOffsetGlobal, &
               mapRes=scale, Re=Re)
          if(lonGrid .le. -900.)then
             i0ceilA=MISSING_INT
             bilintParams%i0ceilALocal(i)=MISSING_INT
          else
             i0ceilA=int(igrdWrap)+1
             bilintParams%i0ceilALocal(i)=mod((i0ceilA-1),tile_ncol)+1
          endif

       else
          i0ceilA=int(bilintParams%i0ceilG(i))+1
          bilintParams%i0ceilALocal(i)=mod((i0ceilA-1),tile_ncol)+1
       endif
! Corresponding tile column index
       if(i0ceilA .ne. MISSING_INT)then
          bilintParams%i0ceil_tilecol(i)=int((i0ceilA-1)/tile_ncol)
       endif
       
! Checking i1floorG:
       call ij2ll_sinusoid(col=i1floorG, row=bilintParams%jfloorG(i),&
            lat=latGrid, lon=lonGrid, &
            nx=ncolGlobalGrid, ny=nrowGlobalGrid, &
            colGridOffset=colGridOffsetGlobal, &
            rowGridOffset=rowGridOffsetGlobal, &
            mapRes=scale, Re=Re)

       if(lonGrid .le. -900.)then
          if(lonFOV(i) .lt. 0. .and. lonFOV(i) .ge. -180.)then
             print *,'err[getBilintParams]: i1floorG is off-hemisphere on left side of grid'
             call errorHalt(1)       
          endif
!          print *,'NOTE: i1floorG is off-hemisphere on right side of grid; attempting wraparound'
          call ll2ij_sinusoid(lat=latGrid, lon=180., col=i180E, row=j180, &
               nx=ncolGlobalGrid, ny=nrowGlobalGrid, &
               colGridOffset=colGridOffsetGlobal, &
               rowGridOffset=rowGridOffsetGlobal, &
               mapRes=scale, Re=Re)
          di180=i1floorG-i180E
          call ll2ij_sinusoid(lat=latGrid, lon=-180., col=i180W, row=j180, &
               nx=ncolGlobalGrid, ny=nrowGlobalGrid, &
               colGridOffset=colGridOffsetGlobal, &
               rowGridOffset=rowGridOffsetGlobal, &
               mapRes=scale, Re=Re)
          igrdWrap=i180W+di180
!          print *,'i1floorG,i180E,i180W,di180,igrdWrap=',i1floorG,i180E,i180W,di180,igrdWrap
! Check if wrapped grid is on-hemisphere
          call ij2ll_sinusoid(col=igrdWrap, row=bilintParams%jfloorG(i),&
               lat=latGrid, lon=lonGrid, &
               nx=ncolGlobalGrid, ny=nrowGlobalGrid, &
               colGridOffset=colGridOffsetGlobal, &
               rowGridOffset=rowGridOffsetGlobal, &
               mapRes=scale, Re=Re)
          if(lonGrid .le. -900.)then
!             print *,'i1floorG still off-hemisphere: igrdWrap=',igrdWrap
             i1floorA=MISSING_INT
             bilintParams%i1floorALocal(i)=MISSING_INT
          else
             i1floorA=int(igrdWrap)+1
             bilintParams%i1floorALocal(i)=mod((i1floorA-1),tile_ncol)+1
          endif
          
       else
          i1floorA=int(i1floorG)+1
          bilintParams%i1floorALocal(i)=mod((i1floorA-1),tile_ncol)+1
       endif

! Corresponding tile column index
       if(i1floorA .ne. MISSING_INT)then
          bilintParams%i1floor_tilecol(i)=int((i1floorA-1)/tile_ncol)
       endif

! Checking i1ceilG:
       call ij2ll_sinusoid(col=i1ceilG, row=bilintParams%jceilG(i),&
            lat=latGrid, lon=lonGrid, &
            nx=ncolGlobalGrid, ny=nrowGlobalGrid, &
            colGridOffset=colGridOffsetGlobal, &
            rowGridOffset=rowGridOffsetGlobal, &
            mapRes=scale, Re=Re)

       if(lonGrid .le. -900.)then
          if(lonFOV(i) .lt. 0. .and. lonFOV(i) .ge. -180.)then
             print *,'err[getBilintParams]: i1ceilG is off-hemisphere on left side of grid'
             call errorHalt(1)       
          endif
!          print *,'NOTE: i1ceilG is off-hemisphere on right side of grid; attempting wraparound'
          call ll2ij_sinusoid(lat=latGrid, lon=180., col=i180E, row=j180, &
               nx=ncolGlobalGrid, ny=nrowGlobalGrid, &
               colGridOffset=colGridOffsetGlobal, &
               rowGridOffset=rowGridOffsetGlobal, &
               mapRes=scale, Re=Re)
          di180=i1ceilG-i180E
          call ll2ij_sinusoid(lat=latGrid, lon=-180., col=i180W, row=j180, &
               nx=ncolGlobalGrid, ny=nrowGlobalGrid, &
               colGridOffset=colGridOffsetGlobal, &
               rowGridOffset=rowGridOffsetGlobal, &
               mapRes=scale, Re=Re)
          igrdWrap=i180W+di180
!          print *,'i1ceilG,i180E,i180W,di180,igrdWrap=',i1ceilG,i180E,i180W,di180,igrdWrap
! Check if wrapped grid is on-hemisphere
          call ij2ll_sinusoid(col=igrdWrap, row=bilintParams%jceilG(i),&
               lat=latGrid, lon=lonGrid, &
               nx=ncolGlobalGrid, ny=nrowGlobalGrid, &
               colGridOffset=colGridOffsetGlobal, &
               rowGridOffset=rowGridOffsetGlobal, &
               mapRes=scale, Re=Re)
          if(lonGrid .le. -900.)then
!             print *,'i1ceilG still off-hemisphere: igrdWrap=',igrdWrap
             i1ceilA=MISSING_INT
             bilintParams%i1ceilALocal(i)=MISSING_INT
          else
             i1ceilA=int(igrdWrap)+1
             bilintParams%i1ceilALocal(i)=mod((i1ceilA-1),tile_ncol)+1
          endif

       else
          i1ceilA=int(i1ceilG)+1
          bilintParams%i1ceilALocal(i)=mod((i1ceilA-1),tile_ncol)+1
       endif

! Corresponding tile column index
       if(i1ceilA .ne. MISSING_INT)then
          bilintParams%i1ceil_tilecol(i)=int((i1ceilA-1)/tile_ncol)
       endif
       
    enddo

!    print *,'getBilintParams: i0floorG, i0ceilG, jfloorG, jceilG=',i0floorG, i0ceilG, jfloorG, jceilG
!    print *,'getBilintParams: i_jfloorG, i_jceilG=',i_jfloorG, i_jceilG
!    print *,'getBilintParams: i0floorALocal, i1floorALocal',i0floorALocal, i1floorALocal
!    print *,'getBilintParams: i0ceilALocal, i1ceilALocal=',i0ceilALocal, i1ceilALocal
!    print *,'getBilintParams: jfloorALocal, jceilALocal=',jfloorALocal, jceilALocal

  end subroutine getBilintParams

  subroutine getBilintFileInd(sfcgrdDB_current,bilintParams,lonFOV,ifileDB_bilint)

    type(sfcgrdDB), intent(in) :: sfcgrdDB_current
    type(bilintSin), intent(in) :: bilintParams
    real, dimension(:), intent(in) :: lonFOV
    integer, dimension(:,:), pointer :: ifileDB_bilint

! local variables

    integer :: nlatlon, i, j, nBilint

! Find corresponding tile file indices, if present. First initialize to missing

    nlatlon=size(lonFOV)
    if(associated(ifileDB_bilint))then
       nBilint=size(ifileDB_bilint(1,:))
       if(nBilint .ne. 4)then
          print *,'err[getBilintFileInd]: nBilint .ne. 4'
          call errorHalt(1)       
       endif
    else
       nBilint=4
!       print *,'allocating ifileDB_bilint: nlatlon,nBilint=',nlatlon,nBilint
       allocate(ifileDB_bilint(nlatlon,nBilint))
!       print *,'size(ifileDB_bilint)=',size(ifileDB_bilint)
    endif
    ifileDB_bilint=MISSING_INT

! Need to do 4 separate nested loops; each pair of nested loops
! corresponding to the following interpolation coordinates:
! 1:(i0floor,jfloor)
! 2:(i1floor,jfloor)
! 3:(i0ceil,jceil)
! 4:(i1ceil,jceil)

! 1:(i0floor,jfloor)
    loop_nlatlon1: do i=1,nlatlon
       loop_ntilefiles1: do j=1,sfcgrdDB_current%ntilefiles
          if(bilintParams%i0floor_tilecol(i) .eq. sfcgrdDB_current%tile_colglobal(j) .and. &
             bilintParams%jfloor_tilerow(i) .eq. sfcgrdDB_current%tile_rowglobal(j))then
             ifileDB_bilint(i,1)=j
             cycle loop_nlatlon1
          endif
       enddo loop_ntilefiles1
    enddo loop_nlatlon1

! 2:(i1floor,jfloor)
    loop_nlatlon2: do i=1,nlatlon
!       print *,'Searching for matching tile file at:bilintParams%i1floor_tilecol(i),bilintParams%jfloor_tilerow(i)=',&
!            bilintParams%i1floor_tilecol(i),bilintParams%jfloor_tilerow(i)
       loop_ntilefiles2: do j=1,sfcgrdDB_current%ntilefiles
!          print *,'Checking tile: sfcgrdDB_current%tile_colglobal(j),sfcgrdDB_current%tile_rowglobal(j)=',&
!               sfcgrdDB_current%tile_colglobal(j),sfcgrdDB_current%tile_rowglobal(j)
          if(bilintParams%i1floor_tilecol(i) .eq. sfcgrdDB_current%tile_colglobal(j) .and. &
             bilintParams%jfloor_tilerow(i) .eq. sfcgrdDB_current%tile_rowglobal(j))then
             ifileDB_bilint(i,2)=j
!             print *,'getBilintFileInd: found matching tile file: '
!             print *,'bilintParams%i1floor_tilecol(i),bilintParams%jfloor_tilerow(i)=',&
!                  bilintParams%i1floor_tilecol(i),bilintParams%jfloor_tilerow(i)
!             print *,'j,sfcgrdDB_current%tile_colglobal(j),sfcgrdDB_current%tile_rowglobal(j)=',&
!                  j,sfcgrdDB_current%tile_colglobal(j),sfcgrdDB_current%tile_rowglobal(j)
             cycle loop_nlatlon2
          endif
       enddo loop_ntilefiles2
    enddo loop_nlatlon2

! 3:(i0ceil,jceil)
    loop_nlatlon3: do i=1,nlatlon
       loop_ntilefiles3: do j=1,sfcgrdDB_current%ntilefiles
          if(bilintParams%i0ceil_tilecol(i) .eq. sfcgrdDB_current%tile_colglobal(j) .and. &
             bilintParams%jceil_tilerow(i) .eq. sfcgrdDB_current%tile_rowglobal(j))then
             ifileDB_bilint(i,3)=j
             cycle loop_nlatlon3
          endif
       enddo loop_ntilefiles3
    enddo loop_nlatlon3

! 4:(i1ceil,jceil)
    loop_nlatlon4: do i=1,nlatlon
       loop_ntilefiles4: do j=1,sfcgrdDB_current%ntilefiles
          if(bilintParams%i1ceil_tilecol(i) .eq. sfcgrdDB_current%tile_colglobal(j) .and. &
             bilintParams%jceil_tilerow(i) .eq. sfcgrdDB_current%tile_rowglobal(j))then
             ifileDB_bilint(i,4)=j
             cycle loop_nlatlon4
          endif
       enddo loop_ntilefiles4
    enddo loop_nlatlon4

!    print *,'size(ifileDB_bilint)=',size(ifileDB_bilint)
!    print *,'ifileDB_bilint=',ifileDB_bilint


  end subroutine getBilintFileInd

  subroutine calcBilint(data2d_i0jf,data2d_i1jf,data2d_i0jc,data2d_i1jc,&
       i0floorG,&
       i_jfloorG,&
       i0ceilG,&
       i_jceilG,&
       jfloorG,&
       jgrd,&
       data2d_bilint)

    real, dimension(:), intent(in) :: data2d_i0jf, data2d_i1jf, data2d_i0jc, data2d_i1jc
    real, dimension(:), intent(in) :: i0floorG, i_jfloorG, i0ceilG, i_jceilG, jfloorG
    real, dimension(:) :: jgrd
    real, dimension(:) :: data2d_bilint

! local variables

    integer :: i, nptsOut, nptsParams
    real :: data2d_jceil, data2d_jfloor

    nptsOut=size(data2d_bilint)
    nptsParams=size(i_jfloorG)
    if(nptsOut .ne. nptsParams)then
          print *,'err[calcBilint]: nptsOut .ne. nptsParams=', nptsOut, nptsParams
          call errorHalt(1)       
    endif
!    print *,'calcBilint: nptsOut,nptsParams=',nptsOut,nptsParams
    data2d_bilint=MISSING_REAL
    do i=1,nptsOut
!       print *,'calcBilint: i=',i
!       print *,'jfloorG(i)=',jfloorG(i)
       if(data2d_i0jf(i) .eq. MISSING_REAL .or. &
          data2d_i1jf(i) .eq. MISSING_REAL .or. &
          data2d_i0jc(i) .eq. MISSING_REAL .or. &
          data2d_i1jc(i) .eq. MISSING_REAL)cycle

!       print *,'calc data2d_jfloor'
       data2d_jfloor=data2d_i0jf(i)+ &
            (i_jfloorG(i)-i0floorG(i))*(data2d_i1jf(i)-data2d_i0jf(i))

!       print *,'calc data2d_jceil'
       data2d_jceil=data2d_i0jc(i)+ &
            (i_jceilG(i)-i0ceilG(i))*(data2d_i1jc(i)-data2d_i0jc(i))

!       print *,'calc data2d_bilint'
!       print *,'data2d_jfloor=',data2d_jfloor
!       print *,'data2d_jceil=',data2d_jceil
!       print *,'jgrd(i)=',jgrd(i)
!       print *,'jfloorG(i)=',jfloorG(i)

       data2d_bilint(i)=data2d_jfloor+(jgrd(i)-jfloorG(i))*(data2d_jceil-data2d_jfloor)
    enddo

!    print *,'calcBilint: end'
  end subroutine calcBilint

  subroutine getTileLatRange(sfcgrdDB_current, minLat, maxLat, minAbsLat, maxAbsLat)

    type(sfcgrdDB), intent(in) :: sfcgrdDB_current
    real, optional, intent(out) :: minLat, maxLat, minAbsLat, maxAbsLat

    real :: minLat0=1.E+03, maxLat0=-1.E+03, minAbsLat0=1.E+03, maxAbsLat0=-1.E+03
    integer :: itile, i, j
    real :: latGrd, lonGrd

    if(present(minLat))minLat=MISSING_REAL
    if(present(maxLat))maxLat=MISSING_REAL
    if(present(minAbsLat))minAbsLat=MISSING_REAL
    if(present(maxAbsLat))maxAbsLat=MISSING_REAL

    
    print *,'sfcgrdDB_current%ntilefiles=',sfcgrdDB_current%ntilefiles
    do itile=1,sfcgrdDB_current%ntilefiles
       i=sfcgrdDB_current%tile_ncol(itile)/2
!       print *,'',
       do j=0,sfcgrdDB_current%tile_nrow(itile)-1
          call ij2ll_sinusoid(col=real(i)+0.5, row=real(j)+0.5, &
               lat=latGrd, lon=lonGrd, &
               nx=sfcgrdDB_current%tile_ncol(itile), ny=sfcgrdDB_current%tile_nrow(itile), &
               colGridOffset=sfcgrdDB_current%tile_colGridOffset(itile), &
               rowGridOffset=sfcgrdDB_current%tile_rowGridOffset(itile), &
               mapRes=sfcgrdDB_current%tile_scale(itile), Re=sfcgrdDB_current%tile_Re(itile))
          if(latGrd .lt. minLat0)minLat0=latGrd
          if(latGrd .gt. maxLat0)maxLat0=latGrd
          if(abs(latGrd) .lt. minAbsLat0)minAbsLat0=abs(latGrd)
          if(abs(latGrd) .gt. maxAbsLat0)maxAbsLat0=abs(latGrd)

       enddo
    enddo

    if(present(minLat) .and. abs(minLat0) .le. 90.)minLat=minLat0
    if(present(maxLat) .and. abs(maxLat0) .le. 90.)maxLat=maxLat0
    if(present(minAbsLat) .and. abs(minAbsLat0) .le. 90.)minAbsLat=minAbsLat0
    if(present(maxAbsLat) .and. abs(maxAbsLat0) .le. 90.)maxAbsLat=maxAbsLat0



  end subroutine getTileLatRange


  subroutine getSfcDBWindowLocs(sfcgrdDB_current, lat, lon, colsizeWin, rowsizeWin, &
       latsWin, lonsWin, &
       ifileDB, &

       ncolGlobalGrid, nrowGlobalGrid, &
       colGridOffsetGlobal, rowGridOffsetGlobal, &
       igrdGlobal, jgrdGlobal, &
       igrdLocal, jgrdLocal, &
       colTile, rowTile)

    type(sfcgrdDB), intent(in) :: sfcgrdDB_current
    real, intent(in) :: lat, lon
    integer, intent(in) :: colsizeWin, rowsizeWin

    real, dimension(:), pointer :: latsWin, lonsWin
    integer, dimension(:), pointer, optional :: ifileDB
    integer, optional, intent(out) :: ncolGlobalGrid, nrowGlobalGrid
    real, optional, intent(out) :: colGridOffsetGlobal, rowGridOffsetGlobal
    real, dimension(:), pointer, optional :: igrdGlobal, jgrdGlobal
    real, dimension(:), pointer, optional :: igrdLocal, jgrdLocal
    integer, dimension(:), pointer, optional :: colTile, rowTile

    
! Local variables

    integer :: i, j, nlatlon, ii, colmin, colmax, rowmin, rowmax, iwrap
    integer :: ncolGlobalGrid0, nrowGlobalGrid0
    real :: col, row
    real :: colGridOffsetGlobal0, rowGridOffsetGlobal0
    real :: scale0, Re0
    real, allocatable, dimension(:) :: igrdGlobal0, jgrdglobal0

!    print *,'Begin getSfcDBWindowLocs:'
!    print *,'sfcgrdDB_current%ntilefiles=',sfcgrdDB_current%ntilefiles

! Make sure database has at least one tile

    if(sfcgrdDB_current%ntilefiles .lt. 1)then
       print *,'err[getSfcDBWindowLocs]: ntilefiles must be ge 1'
       call errorHalt(1)       
    endif


! Calculate the grid parameters for an input grid with same scale as 
! an individual input tile, but with global coverage

    ncolGlobalGrid0=sfcgrdDB_current%tile_ncolglobal(1)*sfcgrdDB_current%tile_ncol(1)
    nrowGlobalGrid0=sfcgrdDB_current%tile_nrowglobal(1)*sfcgrdDB_current%tile_nrow(1)
    colGridOffsetGlobal0=(real(ncolGlobalGrid0)/2.)
    rowGridOffsetGlobal0=(real(nrowGlobalGrid0)/2.)
    scale0=sfcgrdDB_current%tile_scale(1)
    Re0=sfcgrdDB_current%tile_Re(1)

! Total number of lat/lon locations (i.e. total number of grid locations within the search window)

! Determine actual search window boundaries and therefore, the total number of lats and lons to be returned
! This can be different from simply: nlatlon=colsizeWin*rowsizeWin
! because search window may extend off edges of grid 

! First, determine grid location of FOV on global grid

    call ll2ij_sinusoid(lat=lat, lon=lon, &
         col=col, row=row, &
         nx=ncolGlobalGrid0, ny=nrowGlobalGrid0, &
         colGridOffset=colGridOffsetGlobal0, &
         rowGridOffset=rowGridOffsetGlobal0, &
         mapRes=scale0, &
         Re=Re0)

    call windowBounds(col,row,colsizeWin,rowsizeWin,&
         nrowGlobalGrid0,ncolGlobalGrid0,&
         colmin,colmax,rowmin,rowmax)
    
!  This behavior may occur if search window exceeds map boundaries;

    if (rowmin.gt.rowmax)then
       print *,'err[getSfcDBWindowLocs]: rowmin.gt.rowmax'
       call errorHalt(1)       
    endif

    nlatlon=(colmax-colmin+1)*(rowmax-rowmin+1)
!    print *,'getSfcDBWindowLocs: colmin,colmax,rowmin,rowmax=',colmin,colmax,rowmin,rowmax
!    print *,'getSfcDBWindowLocs: nlatlon=',nlatlon


! Allocate local arrays

    if(allocated(igrdGlobal0))deallocate(igrdGlobal0)
    allocate(igrdGlobal0(nlatlon))
    if(allocated(jgrdGlobal0))deallocate(jgrdGlobal0)
    allocate(jgrdGlobal0(nlatlon))

! Allocate/Initialize required output lats/lons

    if(associated(latsWin))deallocate(latsWin)
    allocate(latsWin(nlatlon))
    latsWin=MISSING_REAL
    if(associated(lonsWin))deallocate(lonsWin)
    allocate(lonsWin(nlatlon))
    lonsWin=MISSING_REAL

! Allocate/Initialize pointer arrays if present

    if(present(igrdGlobal))then
       if(associated(igrdGlobal))deallocate(igrdGlobal)
       allocate(igrdGlobal(nlatlon))
       igrdGlobal=MISSING_REAL
    endif
    if(present(jgrdGlobal))then
       if(associated(jgrdGlobal))deallocate(jgrdGlobal)
       allocate(jgrdGlobal(nlatlon))
       jgrdGlobal=MISSING_REAL
    endif
    if(present(igrdLocal))then
       if(associated(igrdLocal))deallocate(igrdLocal)
       allocate(igrdLocal(nlatlon))
       igrdLocal=MISSING_REAL
    endif
    if(present(jgrdLocal))then
       if(associated(jgrdLocal))deallocate(jgrdLocal)
       allocate(jgrdLocal(nlatlon))
       jgrdLocal=MISSING_REAL
    endif
    if(present(colTile))then
       if(associated(colTile))deallocate(colTile)
       allocate(colTile(nlatlon))
       colTile=MISSING_INT
    endif
    if(present(rowTile))then
       if(associated(rowTile))deallocate(rowTile)
       allocate(rowTile(nlatlon))
       rowTile=MISSING_INT
    endif



! Loop over col/row locations calculating corresponding lats/lons, grid locations
! and tile file indices

!    print *,'Begin loop over col/row locations'

    ii=0
!    do i=1,nlatlon
    do j=rowmin,rowmax
       do i=colmin,colmax
          ii=ii+1


! Get global lat/lon for current for current global col/row

!          print *,'call ij2ll_sinusoid'
!          print *,'ii,i,j,real((i-1)+.5),real((j-1)+.5),',&
!               'ncolGlobalGrid0,nrowGlobalGrid0,colGridOffsetGlobal0,rowGridOffsetGlobal0,',&
!               'scale0,Re0=',&
!               ii,i,j,real((i-1)+.5),real((j-1)+.5),&
!               ncolGlobalGrid0,nrowGlobalGrid0,colGridOffsetGlobal0,rowGridOffsetGlobal0,&
!               scale0,Re0
          call ij2ll_sinusoid(col=real((i-1)+.5),row=real((j-1)+.5),&
               lat=latsWin(ii),lon=lonsWin(ii),&
               nx=ncolGlobalGrid0,ny=nrowGlobalGrid0,&
               colGridOffset=colGridOffsetGlobal0,&
               rowGridOffset=rowGridOffsetGlobal0,&
               mapRes=scale0,Re=Re0)
!          print *,'latsWin(ii),lonsWin(ii)=',latsWin(ii),lonsWin(ii)

! If lon is missing (i.e. off-hemisphere) attempt wrap-around

          if(lonsWin(ii) .eq. MISSING_REAL)then
             call wrapLongitude(latsWin(ii),lonsWin(ii),i,j,&
                  nrowGlobalGrid0,ncolGlobalGrid0,&
                  0.,0.,colGridOffsetGlobal0,rowGridOffsetGlobal0,&
                  scale0,iwrap)
          endif
! Still possible for wrapped lon to result in missing (large window relative to local earth circumference). 
! If true, skip to next grid locatoin
          if(lonsWin(ii) .eq. MISSING_REAL)cycle
          call ll2ij_sinusoid(lat=latsWin(ii), lon=lonsWin(ii), &
               col=igrdGlobal0(ii), row=jgrdGlobal0(ii), &
               nx=ncolGlobalGrid0, ny=nrowGlobalGrid0, &
               colGridOffset=colGridOffsetGlobal0, &
               rowGridOffset=rowGridOffsetGlobal0, &
               mapRes=scale0, &
               Re=Re0)
          if(present(igrdLocal))then
             igrdLocal(ii)=mod(igrdGlobal0(ii),real(sfcgrdDB_current%tile_ncol(1)))
          endif
          if(present(jgrdLocal))then
             jgrdLocal(ii)=mod(jgrdGlobal0(ii),real(sfcgrdDB_current%tile_nrow(1)))
          endif
          if(present(colTile))then
             colTile(ii)=int(igrdGlobal0(ii)/sfcgrdDB_current%tile_ncol(1))
          endif
          if(present(rowTile))then
             rowTile(ii)=int(jgrdGlobal0(ii)/sfcgrdDB_current%tile_nrow(1))
          endif
       enddo
       
    enddo
!    print *,'End loop over col/row locations'
!    print *,'getSfcDBlocs:igrdGlobal0,jgrdGlobal0,igrdLocal,jgrdLocal,colTile,rowTile=',&
!         igrdGlobal0,jgrdGlobal0,igrdLocal,jgrdLocal,colTile,rowTile

! Find corresponding tile file indices, if present. First initialize to missing
    if(present(ifileDB))then
       if(associated(ifileDB))deallocate(ifileDB)
       allocate(ifileDB(nlatlon))
       ifileDB=MISSING_INT
       loop_nlatlon: do i=1,nlatlon
!          print *,'getSfcDBlocs: current lat, lon, computed tile col,row=',&
!               lat(i),lon(i),&
!               int(igrdGlobal0(i)/sfcgrdDB_current%tile_ncol(1)),&
!               int(jgrdGlobal0(i)/sfcgrdDB_current%tile_nrow(1))
          loop_ntilefiles: do j=1,sfcgrdDB_current%ntilefiles
             if(int(igrdGlobal0(i)/sfcgrdDB_current%tile_ncol(1)) .eq. sfcgrdDB_current%tile_colglobal(j) &
                  .and. &
                int(jgrdGlobal0(i)/sfcgrdDB_current%tile_nrow(1)) .eq. sfcgrdDB_current%tile_rowglobal(j))then
                ifileDB(i)=j
!                print *,'getSfcDBlocs: found matching tile file in DB: ifileDB(i),sfcgrdDB_current%tile_colglobal(j),sfcgrdDB_current%tile_rowglobal(j)=',&
!                     ifileDB(i),sfcgrdDB_current%tile_colglobal(j),sfcgrdDB_current%tile_rowglobal(j)
                cycle loop_nlatlon
             endif
          enddo loop_ntilefiles
!          print *,'warning[getSfcDBlocs]: no tile file index match found; index set to missing'
!          print *,'i,lat,lon,igrdGlobal0,jgrdGlobal0=',&
!               i,lat(i),lon(i),igrdGlobal0(i),jgrdGlobal0(i)
       enddo loop_nlatlon
    endif

    if(present(ncolGlobalGrid))ncolGlobalGrid=ncolGlobalGrid0
    if(present(nrowGlobalGrid))nrowGlobalGrid=nrowGlobalGrid0
    if(present(colGridOffsetGlobal))colGridOffsetGlobal=colGridOffsetGlobal0
    if(present(rowGridOffsetGlobal))rowGridOffsetGlobal=rowGridOffsetGlobal0
    if(present(igrdGlobal))igrdGlobal=igrdGlobal0
    if(present(jgrdGlobal))jgrdGlobal=jgrdGlobal0

    if(allocated(igrdGlobal0))deallocate(igrdGlobal0)
    if(allocated(jgrdGlobal0))deallocate(jgrdGlobal0)
!    if(allocated(igrdLocal0))deallocate(igrdLocal0)
!    if(allocated(jgrdLocal0))deallocate(jgrdLocal0)
!    if(allocated(colTile0))deallocate(colTile0)
!    if(allocated(rowTile0))deallocate(rowTile0)

!    print *,'End getSfcDBWindowLocs:'

  end subroutine getSfcDBWindowLocs


  subroutine destroySfcgrdDB(sfcgrdDB_current)

!    Deallocate memory when done with database

     type(sfcgrdDB), intent(inout) :: sfcgrdDB_current

     if (allocated(sfcgrdDB_current%tilefile)) &
         deallocate(sfcgrdDB_current%tilefile)
     if (allocated(sfcgrdDB_current%tile_colglobal)) &
         deallocate(sfcgrdDB_current%tile_colglobal)
     if (allocated(sfcgrdDB_current%tile_rowglobal)) &
         deallocate(sfcgrdDB_current%tile_rowglobal)
     if (allocated(sfcgrdDB_current%tile_ncolglobal)) &
         deallocate(sfcgrdDB_current%tile_ncolglobal)
     if (allocated(sfcgrdDB_current%tile_nrowglobal)) &
         deallocate(sfcgrdDB_current%tile_nrowglobal)
     if (allocated(sfcgrdDB_current%tile_ncol)) &
         deallocate(sfcgrdDB_current%tile_ncol)
     if (allocated(sfcgrdDB_current%tile_nrow)) &
         deallocate(sfcgrdDB_current%tile_nrow)
     if (allocated(sfcgrdDB_current%tile_colGridOffset)) &
         deallocate(sfcgrdDB_current%tile_colGridOffset)
     if (allocated(sfcgrdDB_current%tile_rowGridOffset)) &
         deallocate(sfcgrdDB_current%tile_rowGridOffset)
     if (allocated(sfcgrdDB_current%tile_Re)) &
         deallocate(sfcgrdDB_current%tile_Re)
     if (allocated(sfcgrdDB_current%tile_scale)) &
         deallocate(sfcgrdDB_current%tile_scale)

  end subroutine destroySfcgrdDB

end module sfcgrd_DB_module
