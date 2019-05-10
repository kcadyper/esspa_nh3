module typedEmis_module

  implicit none
  private
  public initTypedEmis, getRandTypedEmis
  
  integer, public :: LANDTYPE_UNKNOWN=-1

  integer :: nSfcTyp, nTypedEmis
  integer :: nchmwIn
  real, allocatable, dimension(:,:,:) :: typedEmisTable
  integer, allocatable, dimension(:) :: SfcTypTable,nobsEmisTable
     
CONTAINS 

  subroutine initTypedEmis(ncfile, frqIn, polIn, frqPolGridded)

    use ncdf_module
    
    character (len=*), intent (in)       :: ncfile
    logical, intent(inout), optional     :: frqPolGridded
    
    integer :: maxSfctyp=50, maxTypedEmis=5000, maxnchmw=300
    real, allocatable, dimension(:,:,:) :: typedEmisTable1
    integer, allocatable, dimension(:) :: SfcTypTable1,nobsEmisTable1
    logical :: typeExists
    logical :: myFrqPolGridded=.false.
    
    integer                                :: ncidEmis,nobsEmis
    integer                                :: landtypeEmis
    real, allocatable, dimension(:)        :: emIn
    real, dimension(:), pointer            :: frqIn
    integer, dimension(:), pointer         :: polIn
    real, dimension(:), allocatable        :: freqIn
    integer                                :: nfreqIn
    integer                                :: varlen1,varlen3
    integer, dimension(2)                  :: varlen2
    integer                                :: i,k

    if(allocated(typedEmisTable1))deallocate(typedEmisTable1)
    allocate(typedEmisTable1(maxnchmw,maxTypedEmis,maxSfctyp))
    if(allocated(SfcTypTable1))deallocate(SfcTypTable1)
    allocate(SfcTypTable1(maxSfctyp))
    if(allocated(nobsEmisTable1))deallocate(nobsEmisTable1)
    allocate(nobsEmisTable1(maxSfctyp))
    
!  initialize arrays

    typedEmisTable1=-1.
    SfcTypTable1=LANDTYPE_UNKNOWN
    nobsEmisTable1=0
    nSfcTyp=0

    typeExists=.FALSE.

!  Begin reading typed emissivity file

    call openNcdfFile(ncidEmis,trim(ncfile),status='old',unlimited_dim_length=nobsEmis)
    nchmwIn=readNcdfDim(fid=ncidEmis,dimName='nEmissivities')

    if(associated(frqIn))deallocate(frqIn)
    if(associated(polIn))deallocate(polIn)
    if(allocated(emIn))deallocate(emIn)
    allocate(frqIn(nchmwIn))
    allocate(polIn(nchmwIn))
    allocate(emIn(nchmwIn))

    call readNcdfAttr(fid=ncidEmis,attrName='mwSurfaceFrequencies', &
                      attr=frqIn)
    call readNcdfAttr(fid=ncidEmis,attrName='mwSurfaceFrequencies', &
                      attr=frqIn)
    call readNcdfAttr(fid=ncidEmis,attrName='mwSurfacePolarities', &
                      attr=polIn)

    !---Check for standard frequencies
    nfreqIn=nchmwIn/2
    if(allocated(freqIn))deallocate(freqIn)
    allocate(freqIn(nfreqIn))
    freqIn(1:nfreqIn)=frqIn(1:nchmwIn:2)
    if(any(freqIn(1:nfreqIn) /= frqIn(2:nchmwIn:2)) .or. &
         any(polIn(1:nchmwIn:2) /= 0) .or. &
         any(polIn(2:nchmwIn:2) /= 1)) then 
       print*,'err[typeEmis_module]: Expecting '// &
            'V&H data for each of standard frequencies'
       call errorHalt(1)
    endif
    myFrqPolGridded=.true.
    if (present(frqPolGridded)) frqPolGridded=myFrqPolGridded

    do i=1,nobsEmis

       call readNcdfData(fid=ncidEmis,var=emIn(1:nchmwIn),varname='Emissivities',&
            varlen=varlen2,record_no=i,status='single_record')
       emIn=emIn/10000.
       call readNcdfData(fid=ncidEmis,var=landtypeEmis,varname='landSurfaceType',&
            varlen=varlen1,record_no=i,status='single_record')
       
       if(i .eq. 1)then
          nSfcTyp=nSfcTyp+1
          nobsEmisTable1(nSfcTyp)=nobsEmisTable1(nSfcTyp)+1
          typedEmisTable1(1:nchmwIn,nobsEmisTable1(nSfcTyp),nSfcTyp)=emIn
          SfcTypTable1(nSfcTyp)=landtypeEmis
       else
          if(any(SfcTypTable1(1:nsfcTyp) .eq. landtypeEmis))then
             do k=1,nsfcTyp
                if(landtypeEmis .eq. SfcTypTable1(k))then
                   typeExists=.TRUE.
                   nobsEmisTable1(k)=nobsEmisTable1(k)+1
                   typedEmisTable1(1:nchmwIn,nobsEmisTable1(k),k)=emIn
                   exit
                endif
             enddo
             if(.not. typeExists)stop 'Error in initTypedEmis: typeExists unexpected .FALSE.'
             typeExists=.FALSE.
          else
             nSfcTyp=nSfcTyp+1
             nobsEmisTable1(nSfcTyp)=nobsEmisTable1(nSfcTyp)+1
             typedEmisTable1(1:nchmwIn,nobsEmisTable1(nSfcTyp),nSfcTyp)=emIn
             SfcTypTable1(nsfcTyp)=landtypeEmis
          endif
       endif

    enddo

    nTypedEmis=maxval(nobsEmisTable1(1:nSfcTyp))
    

    if(allocated(typedEmisTable))deallocate(typedEmisTable)
    allocate(typedEmisTable(nchmwIn,nTypedEmis,nSfcTyp))
    if(allocated(SfcTypTable))deallocate(SfcTypTable)
    allocate(SfcTypTable(nSfcTyp))
    if(allocated(nobsEmisTable))deallocate(nobsEmisTable)
    allocate(nobsEmisTable(nSfcTyp))

    typedEmisTable=typedEmisTable1(1:nchmwIn,1:nTypedEmis,1:nSfcTyp)
    SfcTypTable=SfcTypTable1(1:nSfcTyp)
    nobsEmisTable=nobsEmisTable1(1:nSfcTyp)

    deallocate(typedEmisTable1)
    deallocate(SfcTypTable1)
    deallocate(nobsEmisTable1)
    deallocate(emIn)
    deallocate(freqIn)
    
    call closeNcdfFile(fid=ncidEmis,message='Closed netCDF file with typed emissivities')

    print *,'SfcTypTable=',SfcTypTable(1:nSfcTyp)
    print *,'nobsEmisTable=',nobsEmisTable(1:nSfcTyp)

  end subroutine initTypedEmis



  subroutine getRandTypedEmis(landtype,typedEmis,randnum)

    integer, intent(inout)       :: landtype
    real, dimension(:), pointer  :: typedEmis
    real, intent(in)             :: randnum
    
    integer                      :: k, irandobs
    logical                      :: foundType

    foundType=.FALSE.
    if (associated(typedEmis)) deallocate(typedEmis)
    allocate(typedEmis(nchmwIn))
    

! Patch:
! Reset surface type=13 (Urban) to type=16 (Barren)
! Low res smoothed database does not contain urban scenes
    if(landtype .eq. 13)landtype=16
    
    do k=1,nSfcTyp
       if(landtype .eq. SfcTypTable(k))then
          foundType=.TRUE.
! select profile randomly from all profiles of this sfc type        
          irandobs=nint( randnum*(nobsEmisTable(k)-1.) + 1.)
          typedEmis=typedEmisTable(1:nchmwIn,irandobs,k)
      endif
    enddo

    if(.not. foundType)then
       print *,'WARNING in getRandTypedEmis: landtype=',landtype,&
            ' not found in look up table. Setting typedEmis=-1'
       typedEmis=-1.
    else
       foundType=.FALSE.
    endif
  
  end subroutine getRandTypedEmis

end module typedEmis_module

