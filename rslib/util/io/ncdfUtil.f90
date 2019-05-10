MODULE ncdfUtil
  !
  ! Module that provides FORTRAN 90 style netCDF reading support
  use netcdf
  implicit none

  public :: check
  public :: getDimLen
  public :: readVarNC
  public :: inqDim
  public :: readStringAttr
  interface readVarNC
     module procedure readReal
     module procedure readReal1D
     module procedure readReal2D
     module procedure readReal3D
     module procedure readReal4D
     module procedure readDouble
     module procedure readDouble1D
     module procedure readDouble2D
     module procedure readDouble3D
     module procedure readDouble4D
     module procedure readInt
     module procedure readInt1D
     module procedure readInt2D
     module procedure readInt3D
     module procedure readInt4D
     module procedure readLongInt
     module procedure readLongInt1D
     module procedure readLongInt2D
     module procedure readLongInt3D
     module procedure readLongInt4D
     module procedure readShortInt
     module procedure readShortInt1D
     module procedure readShortInt2D
     module procedure readShortInt3D
     module procedure readShortInt4D

     module procedure readChar1D
  end interface

  private
  character(len=*), parameter :: moduleErr = 'ERR: ncdfUtil::'
  character(len=*), parameter :: moduleWrn = 'WRN: ncdfUtil::'
  logical, parameter :: dbg = .false.


  contains
    ! aux function
   subroutine check(status, varName, fatal)
     integer,                     intent(in) :: status
     character(len=*),  optional, intent(in) :: varName
     logical,           optional, intent(in) :: fatal

     logical             :: fatalLoc
     if (present(fatal)) then
        fatalLoc = fatal
      else
        fatalLoc = .true.
      end if

     if(status /= nf90_noerr) then
       if (present(varName)) print *,'Processing: ', varName
       if (fatalLoc) then
         print *, moduleErr, 'netCDF error: ', status, ' : ', trim(nf90_strerror(status))
        call exit(1)
      else
         print *, moduleWrn, 'netCDF WARNING: ', status, ' : ', trim(nf90_strerror(status))
         print *, 'status', status
    end if
     end if
   end subroutine check

   function inqDim(id, dimName, dimLen)
      integer,                intent(in)   :: id
      character(len=*),       intent(in)   :: dimName
      integer, optional,      intent(inout):: dimLen
      integer                              :: dimId
      logical                              :: inqDim

      if (dbg) print*, ' ncdfUtil::inqDim '
      inqDim = (nf90_noerr == nf90_inq_dimid(id, dimName, dimId))
      if (present(dimLen)) then
        if (inqDim) then
          call check( nf90_inquire_dimension(id, dimId,len=dimLen) )
        else
          dimLen = -1
        end if
      end if
   end function inqDim

   subroutine getDimLen(id, dimName, dimLen)
      integer,                intent(in)   :: id
      character(len=*),       intent(in)   :: dimName
      integer,                intent(inout):: dimLen
      integer                              :: dimId

      if (dbg) print*, ' ncdfUtil::getDimLen '
      call check( nf90_inq_dimid(id, dimName, dimId), varName=dimName )
      call check( nf90_inquire_dimension(id, dimId,len=dimLen) )
   end subroutine getDimLen

! read real*4
   subroutine readReal(id, varName, val, fatal)
      integer,                intent(in)   :: id
      character(len=*),       intent(in)   :: varName
      real*4,                 intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
     if (dbg) print*, ' ncdfUtil::readReal '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readReal
   subroutine readReal1D(id, varName, val, fatal)
      integer,                intent(in)   :: id
      character(len=*),       intent(in)   :: varName
      real*4, dimension(:),   intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readReal1D '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readReal1D
   subroutine readReal2D(id, varName, val, fatal)
      integer,                intent(in)   :: id
      character(len=*),       intent(in)   :: varName
      real*4, dimension(:,:), intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readReal2D '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readReal2D
   subroutine readReal3D(id, varName, val, fatal)
      integer,                intent(in)   :: id
      character(len=*),       intent(in)   :: varName
      real*4,dimension(:,:,:),intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readReal3D '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readReal3D
   subroutine readReal4D(id, varName, val, fatal)
      integer,                intent(in)   :: id
      character(len=*),       intent(in)   :: varName
      real*4,dimension(:,:,:,:),intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readReal4D '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readReal4D

! --    integer*4
   subroutine readInt(id, varName, val, fatal)
      integer,                   intent(in)   :: id
      character(len=*),          intent(in)   :: varName
      integer*4,                 intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readInt '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readInt
   subroutine readInt1D(id, varName, val, fatal)
      integer,                   intent(in)   :: id
      character(len=*),          intent(in)   :: varName
      integer*4, dimension(:),   intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readInt1D '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readInt1D
   subroutine readInt2D(id, varName, val, fatal)
      integer,                   intent(in)   :: id
      character(len=*),          intent(in)   :: varName
      integer*4, dimension(:,:), intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readInt2D '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readInt2D
   subroutine readInt3D(id, varName, val, fatal)
      integer,                   intent(in)   :: id
      character(len=*),          intent(in)   :: varName
      integer*4,dimension(:,:,:),intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readInt3D '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readInt3D
   subroutine readInt4D(id, varName, val, fatal)
      integer,                   intent(in)   :: id
      character(len=*),          intent(in)   :: varName
      integer*4,dimension(:,:,:,:),intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readInt4D '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readInt4D

! --    real*8
   subroutine readDouble(id, varName, val, fatal)
      integer,                intent(in)   :: id
      character(len=*),       intent(in)   :: varName
      real*8,                 intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readDouble '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readDouble
   subroutine readDouble1D(id, varName, val, fatal)
      integer,                intent(in)   :: id
      character(len=*),       intent(in)   :: varName
      real*8, dimension(:),   intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readDouble1D '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readDouble1D
   subroutine readDouble2D(id, varName, val, fatal)
      integer,                intent(in)   :: id
      character(len=*),       intent(in)   :: varName
      real*8, dimension(:,:), intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readDouble2D '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readDouble2D
   subroutine readDouble3D(id, varName, val, fatal)
      integer,                intent(in)   :: id
      character(len=*),       intent(in)   :: varName
      real*8,dimension(:,:,:),intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readDouble3D '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readDouble3D
   subroutine readDouble4D(id, varName, val, fatal)
      integer,                intent(in)   :: id
      character(len=*),       intent(in)   :: varName
      real*8,dimension(:,:,:,:),intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readDouble4D '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readDouble4D

! --    integer*8
   subroutine readLongInt(id, varName, val, fatal)
      integer,                   intent(in)   :: id
      character(len=*),          intent(in)   :: varName
      integer*8,                 intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readLongInt '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readLongInt
   subroutine readLongInt1D(id, varName, val, fatal)
      integer,                   intent(in)   :: id
      character(len=*),          intent(in)   :: varName
      integer*8, dimension(:),   intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readLongInt1D '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readLongInt1D
   subroutine readLongInt2D(id, varName, val, fatal)
      integer,                   intent(in)   :: id
      character(len=*),          intent(in)   :: varName
      integer*8, dimension(:,:), intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readLongInt2D
   subroutine readLongInt3D(id, varName, val, fatal)
      integer,                   intent(in)   :: id
      character(len=*),          intent(in)   :: varName
      integer*8,dimension(:,:,:),intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readLongInt3D '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readLongInt3D
   subroutine readLongInt4D(id, varName, val, fatal)
      integer,                   intent(in)   :: id
      character(len=*),          intent(in)   :: varName
      integer*8,dimension(:,:,:,:),intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readLongInt4D '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readLongInt4D

! --    integer*2
   subroutine readShortInt(id, varName, val, fatal)
      integer,                   intent(in)   :: id
      character(len=*),          intent(in)   :: varName
      integer*2,                 intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readShortInt
   subroutine readShortInt1D(id, varName, val, fatal)
      integer,                   intent(in)   :: id
      character(len=*),          intent(in)   :: varName
      integer*2, dimension(:),   intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readShortInt '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readShortInt1D
   subroutine readShortInt2D(id, varName, val, fatal)
      integer,                   intent(in)   :: id
      character(len=*),          intent(in)   :: varName
      integer*2, dimension(:,:), intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readShortInt2D '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readShortInt2D
   subroutine readShortInt3D(id, varName, val, fatal)
      integer,                   intent(in)   :: id
      character(len=*),          intent(in)   :: varName
      integer*2,dimension(:,:,:),intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readShortInt3D '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readShortInt3D
   subroutine readShortInt4D(id, varName, val, fatal)
      integer,                   intent(in)   :: id
      character(len=*),          intent(in)   :: varName
      integer*2,dimension(:,:,:,:),intent(inout):: val
      integer                              :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readShortInt4D '
      call check(nf90_inq_varid(id, varName, varid=varId), varName, fatal)
      call check(nf90_get_var(id, varId, val), varName, fatal)
   end subroutine readShortInt4D

   subroutine readChar1D(id, varName, val)
      use iso_c_binding, only: c_char, C_NULL_CHAR
      interface
        subroutine readCharData(ncid, varName, buf) bind(C, name="readCharData")
          use iso_c_binding, only: c_char
          integer, value, intent(in)        :: ncid
          character(kind=c_char)            :: varName(*)
          character(kind=c_char)            :: buf(*)
        end subroutine readCharData
      end interface

      integer,                   intent(in)   :: id
      character(len=*),          intent(in)   :: varName
      character(len=1), dimension(:),   intent(inout):: val

      if (dbg) print*, ' ncdfUtil::readChar1D '
      call readCharData(id, C_CHAR_""//varName//C_NULL_CHAR, val)
   end subroutine readChar1D

   subroutine readStringAttr(id, attrName, val)
      use iso_c_binding, only: c_char, C_NULL_CHAR
      interface
        subroutine readCharAttribute(ncid, attrName, buf) bind(C, name="readCharAttribute")
          use iso_c_binding, only: c_char
          integer, value, intent(in)        :: ncid
          character(kind=c_char)            :: attrName(*)
          character(kind=c_char)            :: buf(*)
        end subroutine readCharAttribute
      end interface

      integer,                   intent(in)   :: id
      character(len=*),          intent(in)   :: attrName
      character(len=*),          intent(inout):: val

      if (dbg) print*, ' ncdfUtil::readCharAttribute '
      call readCharAttribute(id, C_CHAR_""//attrName//C_NULL_CHAR, val)
   end subroutine readStringAttr

end module ncdfUtil