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

module ncdf_module
  implicit none
  include 'netcdf.inc'
  !--
  public :: openNcdfFile, closeNcdfFile, &
       writeNcdfAttr, writeNcdfDim, writeNcdfData,&
       readNcdfAttr, readNcdfDim, readNcdfData, readNcdfNDim, &
       replaceNcdfData, appendNcdfData, getNcidIndex, varExists
  integer, parameter, public :: &
       MAX_NAME_LENGTH = 200, &
       MAXFILE = 50,&
       MAXVAR = 100,&
       MAX_DIM_NAME_LEN = 10,&
       MAX_NO_GLOBAL_DIM = 20,&
       MAX_ATTR_LENGTH = 2000

  !--
  integer, parameter, private :: dimA=1, dimB=2, dimC=3, dimD=4
  logical, dimension(MAXFILE), private, save :: &
       ncid_index_active=.false.
  integer, dimension(MAXFILE), private, save :: ncid=0
  integer, dimension(MAXFILE), private, save :: unlimited_dim_id=0
  character (len=MAX_NAME_LENGTH), private, save :: &
       myfile = 'myncdf.nc', myVarName, myVarLenName,&
       myVarLenName1, myVarLenName2, myVarLenName3, myVarLenName4,&
       fStatus, myStrLenName
  character (len=MAX_NAME_LENGTH), &
       dimension(MAXFILE), private, save :: myUnlimited_dim_name=' '
  character (len=MAX_DIM_NAME_LEN), &
       dimension(MAX_NO_GLOBAL_DIM), private, save ::&
       globalVarLenName
  character (len=MAX_NAME_LENGTH), &
       dimension(MAXVAR), private, save :: tmpVarName='\0'
  integer, dimension(MAX_NO_GLOBAL_DIM), private, save ::&
       globalVarLenId
  integer, dimension(MAXVAR,MAXFILE), private, save :: varid=0
  integer, private, save :: iname=0
  integer, private :: ncStatus
  integer, private, save :: ncid_cnt_w=0, ncid_index=0,&
       varLen1, varLen2, varLen3, varLen4, strlen, unlimited_length
  integer, dimension(dimA), private :: varDim
  integer, dimension(dimB), private :: var2Dim, szB
  integer, dimension(dimC), private :: var3Dim, szC
  integer, dimension(dimD), private :: var4Dim, szD
  logical, dimension(MAXVAR,MAXFILE), private, save :: &
       existed=.false.
  integer, dimension(MAXVAR,MAXFILE), private, save :: &
       record_no=1
  integer, private :: irec, io, iopen, i, j, k, m
  !--
  private :: openFile, &
       writeFloat1Attribute, writeFloatAttribute, &
       writeInt1Attribute, writeIntAttribute, &
       writeText1Attribute, writeGlobalDim, &
       readFloat1Attribute, readFloatAttribute,&
       readInt1Attribute, readIntAttribute,readVarDim,&
       readText1Attribute, readGlobalDim,&
       writeInt1DParser, writeInt2DParser,&
       writeInt3DParser, writeInt4DParser,&
       writeInt1Byte1DParser, &
       writeInt1Byte2DParser, writeInt1Byte3DParser,&
       writeInt1Byte4DParser,&
       writeInt2Byte1DParser, writeInt2Byte2DParser, &
       writeInt2Byte3DParser, writeInt2Byte4DParser,&
       writeFlt1DParser,&
       writeFlt2DParser, writeFlt3DParser, writeFlt4DParser,&
       writeDFlt1DParser, writeDFlt2DParser,&
       writeDFlt3DParser, writeStr1DParser,&
       writeStr2DParser, writeStr3DParser,&
       writeInteger1D, writeInteger2D, writeInteger3D,&
       writeInteger4D,&
       writeInteger1Byte1D,writeInteger1Byte2D,&
       writeInteger1Byte3D, writeInteger1Byte4D,&
       writeInteger2Byte1D,writeInteger2Byte2D,&
       writeInteger2Byte3D, writeInteger2Byte4D,&
       writeFloat1D, writeFloat2D, writeFloat3D, writeFloat4D,&
       writeDFloat1D, writeDFloat2D, writeDFloat3D,writeDFloat4D,&
       writeString1D, writeString2D, writeString3D,&
       readInt1DParser, readInt2DParser, readInt3DParser,&
       readInt4DParser,&
       readFlt1DParser, readFlt2DParser, readFlt3DParser,&
       readFlt4DParser,&
       readDFlt1DParser, readDFlt2DParser, readDFlt3DParser,&
       readStr1DParser, readStr2DParser, readStr3DParser,&
       readInteger1D, readInteger2D, readInteger3D, readInteger4D,&
       readInteger1Byte1D, readInteger1Byte2D,&
       readInteger1Byte3D, readInteger1Byte4D,&
       readInteger2Byte1D, readInteger2Byte2D,&
       readInteger2Byte3D, readInteger2Byte4D,&
       readFloat1D, readFloat2D, readFloat3D,&
       readFloat4D, readDFloat1D,&
       readDFloat2D, readDFloat3D, readString1D, &
       readString2D, readString3D, &
       replaceInteger1D, replaceInteger2D, replaceInteger3D, replaceInteger4D, &
       replaceInteger1Byte1D, replaceInteger1Byte2D, &
       replaceInteger1Byte3D, replaceInteger1Byte4D, &
       replaceInteger2Byte1D, replaceInteger2Byte2D, &
       replaceInteger2Byte3D, replaceInteger2Byte4D, &
       replaceFloat1D, replaceFloat2D, replaceFloat3D, replaceFloat4D, &
       replaceDFloat1D, replaceDFloat2D, replaceDFloat3D, replaceDFloat4D, &
       replaceString1D, replaceString2D, replaceString3D, &
       appendInteger1D, appendInteger2D, appendInteger3D, appendInteger4D, &
       appendInteger1Byte1D, appendInteger1Byte2D, &
       appendInteger1Byte3D, appendInteger1Byte4D, &
       appendInteger2Byte1D, appendInteger2Byte2D, &
       appendInteger2Byte3D, appendInteger2Byte4D, &
       appendFloat1D, appendFloat2D, appendFloat3D, appendFloat4D, &
       appendDFloat1D, appendDFloat2D, appendDFloat3D, appendDFloat4D, &
       appendString1D, appendString2D, appendString3D, &
       replaceInt1DParser, replaceInt2DParser, replaceInt3DParser, replaceInt4DParser, &
       replaceInt1Byte1DParser, replaceInt1Byte2DParser, &
       replaceInt1Byte3DParser, replaceInt1Byte4DParser, &
       replaceInt2Byte1DParser, replaceInt2Byte2DParser, &
       replaceInt2Byte3DParser, replaceInt2Byte4DParser, &
       replaceFlt1DParser, replaceFlt2DParser, replaceFlt3DParser, replaceFlt4DParser, &
       replaceDFlt1DParser, replaceDFlt2DParser, &
       replaceDFlt3DParser, replaceDFlt4DParser, &
       replaceStr1DParser, replaceStr2DParser, replaceStr3DParser, &
       appendInt1DParser, appendInt2DParser, appendInt3DParser, appendInt4DParser, &
       appendInt1Byte1DParser, appendInt1Byte2DParser, &
       appendInt1Byte3DParser, appendInt1Byte4DParser, &
       appendInt2Byte1DParser, appendInt2Byte2DParser, &
       appendInt2Byte3DParser, appendInt2Byte4DParser, &
       appendFlt1DParser, appendFlt2DParser, appendFlt3DParser, appendFlt4DParser, &
       appendDFlt1DParser, appendDFlt2DParser, &
       appendDFlt3DParser, appendDFlt4DParser, &
       appendStr1DParser, appendStr2DParser, appendStr3DParser, &
       closeFile, FindString, configureName, getNewNcidIndex,&
       handle_err,&
       writeDouble1Attribute, writeDoubleAttribute, &
       readDouble1Attribute, readDoubleAttribute
  !--
  interface openNcdfFile
     module procedure openFile
  end interface
  !--
  interface writeNcdfAttr
     module procedure writeDouble1Attribute
     module procedure writeDoubleAttribute
     module procedure writeFloat1Attribute
     module procedure writeFloatAttribute
     module procedure writeInt1Attribute
     module procedure writeIntAttribute
     module procedure writeText1Attribute
  end interface
  !--
  interface writeNcdfDim
     module procedure writeGlobalDim
  end interface
  !--
  interface readNcdfAttr
     module procedure readDouble1Attribute
     module procedure readDoubleAttribute
     module procedure readFloat1Attribute
     module procedure readFloatAttribute
     module procedure readInt1Attribute
     module procedure readIntAttribute
     module procedure readText1Attribute
  end interface
  !--
  interface readNcdfDim
     module procedure readGlobalDim
     module procedure readUnlimitDim
  end interface
  !--
  interface readNcdfNDim
     module procedure readVarDim
  end interface
  !--
  interface writeNcdfData
     module procedure writeInt1DParser
     module procedure writeInt2DParser
     module procedure writeInt3DParser
     module procedure writeInt4DParser
     module procedure writeInt1Byte1DParser
     module procedure writeInt1Byte2DParser
     module procedure writeInt1Byte3DParser
     module procedure writeInt1Byte4DParser
     module procedure writeInt2Byte1DParser
     module procedure writeInt2Byte2DParser
     module procedure writeInt2Byte3DParser
     module procedure writeInt2Byte4DParser

     module procedure writeFlt1DParser
     module procedure writeFlt2DParser
     module procedure writeFlt3DParser
     module procedure writeFlt4DParser

     module procedure writeDFlt1DParser
     module procedure writeDFlt2DParser
     module procedure writeDFlt3DParser
     module procedure writeDFlt4DParser

     module procedure writeStr1DParser
     module procedure writeStr2DParser
     module procedure writeStr3DParser
     module procedure writeInteger1D
     module procedure writeInteger2D
     module procedure writeInteger3D
     module procedure writeInteger4D
     module procedure writeInteger1Byte1D
     module procedure writeInteger1Byte2D
     module procedure writeInteger1Byte3D
     module procedure writeInteger1Byte4D
     module procedure writeInteger2Byte1D
     module procedure writeInteger2Byte2D
     module procedure writeInteger2Byte3D
     module procedure writeInteger2Byte4D

     module procedure writeFloat1D
     module procedure writeFloat2D
     module procedure writeFloat3D
     module procedure writeFloat4D

     module procedure writeDFloat1D
     module procedure writeDFloat2D
     module procedure writeDFloat3D
     module procedure writeDFloat4D

     module procedure writeString1D
     module procedure writeString2D
     module procedure writeString3D
  end interface
  !--
  interface replaceNcdfData
     module procedure replaceInteger1D
     module procedure replaceInteger2D
     module procedure replaceInteger3D
     module procedure replaceInteger4D
     module procedure replaceInteger1Byte1D
     module procedure replaceInteger1Byte2D
     module procedure replaceInteger1Byte3D
     module procedure replaceInteger1Byte4D
     module procedure replaceInteger2Byte1D
     module procedure replaceInteger2Byte2D
     module procedure replaceInteger2Byte3D
     module procedure replaceInteger2Byte4D
     module procedure replaceFloat1D
     module procedure replaceFloat2D
     module procedure replaceFloat3D
     module procedure replaceFloat4D
     module procedure replaceDFloat1D
     module procedure replaceDFloat2D
     module procedure replaceDFloat3D
     module procedure replaceDFloat4D
     module procedure replaceString1D
     module procedure replaceString2D
     module procedure replaceString3D
     module procedure replaceInt1DParser
     module procedure replaceInt2DParser
     module procedure replaceInt3DParser
     module procedure replaceInt4DParser
     module procedure replaceInt1Byte1DParser
     module procedure replaceInt1Byte2DParser
     module procedure replaceInt1Byte3DParser
     module procedure replaceInt1Byte4DParser
     module procedure replaceInt2Byte1DParser
     module procedure replaceInt2Byte2DParser
     module procedure replaceInt2Byte3DParser
     module procedure replaceInt2Byte4DParser
     module procedure replaceFlt1DParser
     module procedure replaceFlt2DParser
     module procedure replaceFlt3DParser
     module procedure replaceFlt4DParser
     module procedure replaceDFlt1DParser
     module procedure replaceDFlt2DParser
     module procedure replaceDFlt3DParser
     module procedure replaceDFlt4DParser
     module procedure replaceStr1DParser
     module procedure replaceStr2DParser
     module procedure replaceStr3DParser
  end interface
  !--
  interface appendNcdfData
     module procedure appendInteger1D
     module procedure appendInteger2D
     module procedure appendInteger3D
     module procedure appendInteger4D
     module procedure appendInteger1Byte1D
     module procedure appendInteger1Byte2D
     module procedure appendInteger1Byte3D
     module procedure appendInteger1Byte4D
     module procedure appendInteger2Byte1D
     module procedure appendInteger2Byte2D
     module procedure appendInteger2Byte3D
     module procedure appendInteger2Byte4D
     module procedure appendFloat1D
     module procedure appendFloat2D
     module procedure appendFloat3D
     module procedure appendFloat4D
     module procedure appendDFloat1D
     module procedure appendDFloat2D
     module procedure appendDFloat3D
     module procedure appendDFloat4D
     module procedure appendString1D
     module procedure appendString2D
     module procedure appendString3D
     module procedure appendInt1DParser
     module procedure appendInt2DParser
     module procedure appendInt3DParser
     module procedure appendInt4DParser
     module procedure appendInt1Byte1DParser
     module procedure appendInt1Byte2DParser
     module procedure appendInt1Byte3DParser
     module procedure appendInt1Byte4DParser
     module procedure appendInt2Byte1DParser
     module procedure appendInt2Byte2DParser
     module procedure appendInt2Byte3DParser
     module procedure appendInt2Byte4DParser
     module procedure appendFlt1DParser
     module procedure appendFlt2DParser
     module procedure appendFlt3DParser
     module procedure appendFlt4DParser
     module procedure appendDFlt1DParser
     module procedure appendDFlt2DParser
     module procedure appendDFlt3DParser
     module procedure appendDFlt4DParser
     module procedure appendStr1DParser
     module procedure appendStr2DParser
     module procedure appendStr3DParser
  end interface
  !--
  interface readNcdfData
     module procedure readInt1DParser
     module procedure readInt2DParser
     module procedure readInt3DParser
     module procedure readInt4DParser
     module procedure readInt1Byte1DParser
     module procedure readInt1Byte2DParser
     module procedure readInt1Byte3DParser
     module procedure readInt1Byte4DParser
     module procedure readInt2Byte1DParser
     module procedure readInt2Byte2DParser
     module procedure readInt2Byte3DParser
     module procedure readInt2Byte4DParser
     module procedure readFlt1DParser
     module procedure readFlt2DParser
     module procedure readFlt3DParser
     module procedure readFlt4DParser
     module procedure readDFlt1DParser
     module procedure readDFlt2DParser
     module procedure readDFlt3DParser
     module procedure readStr1DParser
     module procedure readStr2DParser
     module procedure readStr3DParser
     module procedure readInteger1D
     module procedure readInteger2D
     module procedure readInteger3D
     module procedure readInteger1Byte1D
     module procedure readInteger1Byte2D
     module procedure readInteger1Byte3D
     module procedure readInteger1Byte4D
     module procedure readInteger2Byte1D
     module procedure readInteger2Byte2D
     module procedure readInteger2Byte3D
     module procedure readInteger2Byte4D
     module procedure readFloat1D
     module procedure readFloat2D
     module procedure readFloat3D
     module procedure readFloat4D
     module procedure readDFloat1D
     module procedure readDFloat2D
     module procedure readDFloat3D
     module procedure readString1D
     module procedure readString2D
     module procedure readString3D
  end interface
  !--
  interface closeNcdfFile
     module procedure closeFile
  end interface
  !--
contains
  subroutine openFile(fid, file, status, unlimited_dim_name, &
       unlimited_dim_length)
    character (len=*), intent (in), optional :: file, status
    integer, intent (inout) :: fid
    character (len = MAX_NAME_LENGTH), save :: localFile
    character (len=*), intent (in), optional :: &
         unlimited_dim_name
    integer, intent (inout), optional :: unlimited_dim_length
    integer :: dimid, ncidtmp

    if (present(file)) myfile = file
    if (present(status)) then
       fStatus = status
    else
       fStatus = 'old'
    endif

    if (fStatus == 'new') then
       if (.not. present(file))&
            print *,'creating NCDF File: ', trim(myFile)
       ncStatus = nf_create(myfile,nf_write,ncidtmp)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,trim(myfile)//' in openFile')
       if (ncStatus < 0) then
          return
       endif
       call getNewNcidIndex(ncid_index)
       ncid(ncid_index) = ncidtmp
       fid = ncid(ncid_index)
       if (present(unlimited_dim_name)) then
          myUnlimited_dim_name(ncid_index) = unlimited_dim_name
       else
          myUnlimited_dim_name(ncid_index) = 'unlimitedDim'
       endif
       ncStatus = nf_def_dim(fid, &
            trim(myUnlimited_dim_name(ncid_index)), &
            nf_unlimited, unlimited_dim_id(ncid_index))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,myUnlimited_dim_name(ncid_index))
       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,myfile)
    else
       if (.not. present(file))&
            print *,'opening NCDF File: ', trim(myFile)

       if (fStatus == 'append' .or. fStatus == 'replace') then
          ncStatus = nf_open(trim(myFile),nf_write,ncidtmp)
       else
          ncStatus = nf_open(myFile, nf_nowrite, ncidtmp)
       endif
       if (ncStatus .ne. nf_noerr) then
          call handle_err(ncStatus,myFile)
          call exit(1)
       endif
       if (ncStatus < 0) then
          return
       endif
       call getNewNcidIndex(ncid_index)
       ncid(ncid_index) = ncidtmp
       fid = ncid(ncid_index)
       ncStatus = nf_inq_unlimdim(fid, unlimited_dim_id(ncid_index))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,myUnlimited_dim_name(ncid_index))
       if (unlimited_dim_id(ncid_index) > 0) then
          ncStatus = nf_inq_dim(fid, unlimited_dim_id(ncid_index), &
               myUnlimited_dim_name(ncid_index), unlimited_length)
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,&
               myUnlimited_dim_name(ncid_index))
          if (present(unlimited_dim_length)) then
             unlimited_dim_length = unlimited_length
          endif
       endif
    endif

  end subroutine openFile
  !--
  subroutine writeFloat1Attribute(fid, attr, attrName)
    integer, intent (in) :: fid
    real*4, intent (in) :: attr
    character (len=*), intent (in) :: attrName

    ncStatus = nf_redef(fid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,trim(attrName))

    ncStatus = nf_put_att_real(fid, nf_global, &
         trim(attrName), nf_float, 1, attr)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(attrName))

    ncStatus = nf_enddef(fid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,attrName)

  end subroutine writeFloat1Attribute
  !--
  subroutine writeFloatAttribute(fid, attr, attrName)
    integer, intent (in) :: fid
    real*4, dimension(:), intent (in) :: attr
    character (len=*), intent (in) :: attrName

    ncStatus = nf_redef(fid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,trim(attrName))

    ncStatus = nf_put_att_real(fid, nf_global, &
         trim(attrName), nf_float, size(attr), attr)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(attrName))

    ncStatus = nf_enddef(fid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,attrName)

  end subroutine writeFloatAttribute
  !--
  subroutine writeDouble1Attribute(fid, attr, attrName)
    integer, intent (in) :: fid
    real*8, intent (in) :: attr
    character (len=*), intent (in) :: attrName

    real*4              :: attrLoc
    attrLoc = attr
    ncStatus = nf_redef(fid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,trim(attrName))

    ncStatus = nf_put_att_real(fid, nf_global, &
         trim(attrName), nf_float, 1, attrLoc)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(attrName))

    ncStatus = nf_enddef(fid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,attrName)

  end subroutine writeDouble1Attribute
  !--

  subroutine writeDoubleAttribute(fid, attr, attrName)
    integer, intent (in) :: fid
    real*8, dimension(:), intent (in) :: attr
    character (len=*), intent (in) :: attrName

    real*4, dimension(size(attr)) :: attrLoc

    attrLoc = attr
    ncStatus = nf_redef(fid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,trim(attrName))

    ncStatus = nf_put_att_real(fid, nf_global, &
         trim(attrName), nf_float, size(attr), attrLoc)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(attrName))

    ncStatus = nf_enddef(fid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,attrName)

  end subroutine writeDoubleAttribute
  !--
  subroutine writeInt1Attribute(fid, attr, attrName)
    integer, intent (in) :: fid
    integer, intent (in) :: attr
    character (len=*), intent (in) :: attrName

    ncStatus = nf_redef(fid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,trim(attrName))

    ncStatus = nf_put_att_int(fid, nf_global, &
         trim(attrName), nf_int, 1, attr)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(attrName))

    ncStatus = nf_enddef(fid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,attrName)

  end subroutine writeInt1Attribute
  !--
  subroutine writeIntAttribute(fid, attr, attrName)
    integer, intent (in) :: fid
    integer, dimension(:), intent (in) :: attr
    character (len=*), intent (in) :: attrName

    ncStatus = nf_redef(fid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,trim(attrName))

    ncStatus = nf_put_att_int(fid, nf_global, &
         trim(attrName), nf_int, size(attr), attr)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(attrName))

    ncStatus = nf_enddef(fid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,attrName)

  end subroutine writeIntAttribute
  !--
  subroutine writeText1Attribute(fid, attr, attrName)
    integer, intent (in) :: fid
    character (len=*), intent (in) :: attr
    character (len=*), intent (in) :: attrName

    ncStatus = nf_redef(fid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,trim(attrName))

    ncStatus = nf_put_att_text(fid, nf_global, &
         trim(attrName), len_trim(attr), trim(attr))
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(attrName))

    ncStatus = nf_enddef(fid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,attrName)

  end subroutine writeText1Attribute
  !--
  subroutine writeGlobalDim(fid, dimlen, dimName)
    integer, intent (in) :: fid
    integer, intent (in) :: dimlen
    character (len=*), intent (in) :: dimName
    integer :: tmpid
    character (len=256) :: fidstr

    ncStatus = nf_redef(fid)
    write(fidstr,*) fid
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,trim(fidstr))

    ncStatus = nf_def_dim(fid,trim(dimName),dimlen,tmpid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(dimName))

    ncStatus = nf_enddef(fid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,dimName)

  end subroutine writeGlobalDim
  !--
  subroutine writeInteger1D(fid, var, varName, varLenName, &
       varUnit, varLongName, append, unlimited_dim)
    integer, dimension(:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varLenName, varUnit
    logical, intent (in), optional :: append
    integer, intent (in) :: fid
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer :: dimid, iopen
    integer, dimension(1) :: myVarDim

    call configureName(fid,iopen,define,varName=varName)

    if (present(varLenName)) then
       myvarLenName = trim(varLenName)
    else
       myvarLenName = trim(tmpVarName(iname))//'Len'
    endif

    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(unlimited_dim)) then
          myvarDim(1) = unlimited_dim_id(iopen)
       else
          if (trim(myVarLenName) == &
               trim(myUnlimited_dim_name(iopen))) then
             myvarDim(1) = unlimited_dim_id(iopen)
          else
             ncStatus = nf_def_dim(fid,myVarLenName,&
                  size(var),dimid)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName),dimid)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,tmpVarName(iname))
             endif
             myvarDim(1) = dimid
          endif
       endif

       ncStatus = nf_def_var(fid, trim(tmpVarName(iname)), nf_int, &
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,tmpVarName(iname))
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),'units', &
               len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,tmpVarName(iname))
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       ncStatus = nf_put_vara_int(fid, &
            varid(iname,iopen),1,size(var),var)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    else
       ncStatus = nf_put_var1_int(fid, &
            varid(iname,iopen), record_no(iname,iopen),var(1))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeInteger1D
  !--
  subroutine replaceInteger1D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    integer, dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then   !writes all records
       ncStatus = nf_put_vara_int(fid,vid,start_rec,size(var),var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
    else                              !writes one record a time
       ncStatus = nf_put_var1_int(fid,vid,start_rec,var(1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
    endif
  end subroutine replaceInteger1D
  !--
  subroutine appendInteger1D(fid, var, varName, append)
    integer, intent(in) :: fid
    integer, dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceInteger1D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendInteger1D
  !--
  subroutine writeInteger2D(fid, var, varName, varLenName, &
       varUnit, varLongName, append, unlimited_dim)
    integer, dimension(:,:), intent (in) :: var
    character (len=*), intent (in), optional ::&
         varName, varLongName, varUnit
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    logical, intent (in), optional :: append
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer:: dimid1, dimid2, iopen
    integer, dimension(2) :: myVarDim

    call configureName(fid,iopen,define,varName=varName)

    if (present(varLenName)) then
       if (size(varLenName) /= 2) then
          print *,'err[ncdf_module::writeInteger2D]: ',&
               ' 2-element array required.'
          call exit(1)
       endif
       myvarLenName1 = varLenName(1)
       myvarLenName2 = varLenName(2)
    else
       myvarLenName1 = trim(tmpVarName(iname))//'Len1'
       myvarLenName2 = trim(tmpVarName(iname))//'Len2'
    endif

    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(unlimited_dim)) then
          select case (unlimited_dim)
          case (1)
             myvarDim(1) = unlimited_dim_id(iopen)
             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var(1,:)),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,tmpVarName(iname))
             endif
             myvarDim(2) = dimid2
          case (2)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var(:,1)),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,tmpVarName(iname))
             endif
             myvarDim(1) = dimid1
             myvarDim(2) = unlimited_dim_id(iopen)
          case default
             print *,'err[ncdf_module::writeInteger2D]: ',&
                  ' unlimited_dim for 2D integer: 1 or 2'
             call exit(1)
          end select
       else
          ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
               size(var(:,1)),dimid1)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName1),dimid1)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName1)
          endif
          myvarDim(1) = dimid1

          ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
               size(var(1,:)),dimid2)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName2),dimid2)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myvarLenName2)
          endif
          myvarDim(2) = dimid2
       endif
       ncStatus = nf_def_var(fid, trim(tmpVarName(iname)),nf_int, &
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varLongName)
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),'units', &
               len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varUnit)
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       do irec = 1, size(var(1,:))
          ncStatus = nf_put_vara_int(fid,varid(iname,iopen),&
               (/1,irec/),(/size(var(:,1)),1/),var(:,irec))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,tmpVarName(iname))
       enddo
    else
       ncStatus = nf_put_vara_int(fid,varid(iname,iopen),&
            (/1,record_no(iname,iopen)/),(/size(var,1),1/),var(:,1))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,trim(tmpVarName(iname)))
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeInteger2D
  !--
  subroutine replaceInteger2D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    integer, dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then   !writes all records
       do irec = 1, size(var,2)
          ncStatus = nf_put_vara_int(fid,vid,(/1,start_rec+irec-1/), &
               (/size(var,1),1/),var(:,irec))
          if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
       enddo
    else                              !writes one record a time
       ncStatus = nf_put_vara_int(fid,vid,(/1,start_rec/), &
            (/size(var,1),1/),var(:,1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
    endif
  end subroutine replaceInteger2D
  !--
  subroutine appendInteger2D(fid, var, varName, append)
    integer, intent(in) :: fid
    integer, dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceInteger2D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendInteger2D
  !--
  subroutine writeInteger3D(fid, var, varName, varLenName, &
       varUnit, varLongName, append, unlimited_dim)
    integer, dimension(:,:,:), intent (in) :: var
    character (len=*), intent (in), optional ::&
         varName, varLongName, varUnit
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    logical, intent (in), optional :: append
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer :: dimid1, dimid2, dimid3, iopen
    integer, dimension(3) :: myVarDim

    call configureName(fid,iopen,&
         define,varName=varName)

    if (present(varLenName)) then
       myvarLenName1 = varLenName(1)
       myvarLenName2 = varLenName(2)
       myvarLenName3 = varLenName(3)
    else
       myvarLenName1 = trim(tmpVarName(iname))//'Len1'
       myvarLenName2 = trim(tmpVarName(iname))//'Len2'
       myvarLenName3 = trim(tmpVarName(iname))//'Len3'
    endif

    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(unlimited_dim)) then
          select case (unlimited_dim)
          case (1)
             myvarDim(1) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName3)
             endif
             myvarDim(3) = dimid3
          case (2)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             myvarDim(2) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName3)
             endif
             myvarDim(3) = dimid3
          case (3)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             myvarDim(3) = unlimited_dim_id(iopen)
          case default
             print *,'err[ncdf_module::writeInteger3D]: ',&
                  ' unlimited_dim for 2D integer: 1 or 2'
             call exit(1)
          end select
       else
          ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
               size(var,1),dimid1)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName1),dimid1)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName1)
          endif
          myvarDim(1) = dimid1

          ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
               size(var,2),dimid2)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName2),dimid2)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName2)
          endif
          myvarDim(2) = dimid2

          ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
               size(var,3),dimid3)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName3),dimid3)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName3)
          endif
          myvarDim(3) = dimid3
       endif

       ncStatus = nf_def_var(fid, trim(tmpVarName(iname)), nf_int, &
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varLongName)
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'units', len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varUnit)
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       do j = 1, size(var,3)
          do i = 1, size(var,2)
             ncStatus = nf_put_vara_int(fid,varid(iname,iopen),&
                  (/1,i,j/),(/size(var(:,i,j)),1,1/),var(:,i,j))
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,tmpVarName(iname))
          enddo
       enddo
    else
       do i = 1, size(var,2)
          ncStatus = nf_put_vara_int(fid,varid(iname,iopen),&
               (/1,i,record_no(iname,iopen)/),&
               (/size(var(:,i,1)),1,1/),&
               var(:,i,1))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,tmpVarName(iname))
       enddo
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeInteger3D
  !--
  subroutine replaceInteger3D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    integer, dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then   !writes all records
       do j = 1, size(var,3)
          do i = 1, size(var,2)
             ncStatus = nf_put_vara_int(fid,vid,(/1,i,start_rec+j-1/), &
                  (/size(var,1),1,1/),var(:,i,j))
             if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
          enddo
       enddo
    else                              !writes one record a time
       do i = 1, size(var,2)
          ncStatus = nf_put_vara_int(fid,vid, &
               (/1,i,start_rec/), &
               (/size(var,1),1,1/),var(:,i,1))
          if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
       enddo
    endif
  end subroutine replaceInteger3D
  !--
  subroutine appendInteger3D(fid, var, varName, append)
    integer, intent(in) :: fid
    integer, dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceInteger3D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendInteger3D
  !--
  subroutine writeInteger4D(fid, var, varName, varLenName, &
       varUnit, varLongName, append, unlimited_dim)
    integer, dimension(:,:,:,:), intent (in) :: var
    character (len=*), intent (in), optional ::&
         varName, varLongName, varUnit
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    logical, intent (in), optional :: append
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer :: dimid1, dimid2, dimid3, dimid4, iopen
    integer, dimension(4) :: myVarDim

    call configureName(fid,iopen,&
         define,varName=varName)

    if (present(varLenName)) then
       myvarLenName1 = varLenName(1)
       myvarLenName2 = varLenName(2)
       myvarLenName3 = varLenName(3)
       myvarLenName4 = varLenName(4)
    else
       myvarLenName1 = trim(tmpVarName(iname))//'Len1'
       myvarLenName2 = trim(tmpVarName(iname))//'Len2'
       myvarLenName3 = trim(tmpVarName(iname))//'Len3'
       myvarLenName4 = trim(tmpVarName(iname))//'Len4'
    endif

    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(unlimited_dim)) then
          select case (unlimited_dim)
          case (1)
             myvarDim(1) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName3)
             endif
             myvarDim(3) = dimid3

             ncStatus = nf_def_dim(fid, trim(myvarLenName4),&
                  size(var,4),dimid4)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName4),dimid4)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName4)
             endif
             myvarDim(4) = dimid4
          case (2)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             myvarDim(2) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName3)
             endif
             myvarDim(3) = dimid3

             ncStatus = nf_def_dim(fid, trim(myvarLenName4),&
                  size(var,4),dimid4)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName4),dimid4)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName4)
             endif
             myvarDim(4) = dimid4
          case (3)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             myvarDim(3) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName4),&
                  size(var,4),dimid4)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName4),dimid4)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName4)
             endif
             myvarDim(4) = dimid4
          case (4)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName3)
             endif
             myvarDim(3) = dimid3

             myvarDim(4) = unlimited_dim_id(iopen)
          case default
             print *,'err[ncdf_module::writeInteger4D]: ',&
                  ' unlimited_dim for 2D integer: 1 or 2'
             call exit(1)
          end select
       else
          ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
               size(var,1),dimid1)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName1),dimid1)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName1)
          endif
          myvarDim(1) = dimid1

          ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
               size(var,2),dimid2)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName2),dimid2)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName2)
          endif
          myvarDim(2) = dimid2

          ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
               size(var,3),dimid3)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName3),dimid3)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName3)
          endif
          myvarDim(3) = dimid3

          ncStatus = nf_def_dim(fid, trim(myvarLenName4),&
               size(var,4),dimid4)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName4),dimid4)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName4)
          endif
          myvarDim(4) = dimid4
       endif

       ncStatus = nf_def_var(fid, trim(tmpVarName(iname)), nf_int, &
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varLongName)
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'units', len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varUnit)
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       do k = 1, size(var,4)
          do j = 1, size(var,3)
             do i = 1, size(var,2)
                ncStatus = nf_put_vara_int(fid,varid(iname,iopen),&
                     (/1,i,j,k/),(/size(var(:,i,j,k)),1,1,1/),&
                     var(:,i,j,k))
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,tmpVarName(iname))
             enddo
          enddo
       enddo
    else
       do j = 1, size(var,3)
          do i = 1, size(var,2)
             ncStatus = nf_put_vara_int(fid,varid(iname,iopen),&
                  (/1,i,j,record_no(iname,iopen)/),&
                  (/size(var(:,i,j,1)),1,1,1/),&
                  var(:,i,j,1))
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,tmpVarName(iname))
          enddo
       enddo
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeInteger4D
  !--
  subroutine replaceInteger4D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    integer, dimension(:,:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then   !writes all records
       do k = 1, size(var,4)
          do j = 1, size(var,3)
             do i = 1, size(var,2)
                ncStatus = nf_put_vara_int(fid,vid,(/1,i,j,start_rec+k-1/), &
                     (/size(var,1),1,1,1/),var(:,i,j,k))
                if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
             enddo
          enddo
       enddo
    else                              !writes one record a time
       do j = 1, size(var,3)
          do i = 1, size(var,2)
             ncStatus = nf_put_vara_int(fid,vid, &
                  (/1,i,j,start_rec/), &
                  (/size(var,1),1,1,1/),var(:,i,j,1))
             if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
          enddo
       enddo
    endif
  end subroutine replaceInteger4D
  !--
  subroutine appendInteger4D(fid, var, varName, append)
    integer, intent(in) :: fid
    integer, dimension(:,:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceInteger4D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendInteger4D
  !--
  subroutine writeInteger1Byte1D(fid, var, varName, varLenName, &
       varUnit, varLongName, append, unlimited_dim)
    integer (kind=1), dimension(:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varLenName, varUnit
    logical, intent (in), optional :: append
    integer, intent (in) :: fid
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer :: dimid, iopen
    integer, dimension(1) :: myVarDim

    call configureName(fid,iopen,define,varName=varName)

    if (present(varLenName)) then
       myvarLenName = trim(varLenName)
    else
       myvarLenName = trim(tmpVarName(iname))//'Len'
    endif

    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(unlimited_dim)) then
          myvarDim(1) = unlimited_dim_id(iopen)
       else
          if (trim(myVarLenName) == &
               trim(myUnlimited_dim_name(iopen))) then
             myvarDim(1) = unlimited_dim_id(iopen)
          else
             ncStatus = nf_def_dim(fid,myVarLenName,&
                  size(var),dimid)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName),dimid)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,tmpVarName(iname))
             endif
             myvarDim(1) = dimid
          endif
       endif

       ncStatus = nf_def_var(fid, trim(tmpVarName(iname)), nf_byte, &
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,tmpVarName(iname))
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),'units', &
               len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,tmpVarName(iname))
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       ncStatus = nf_put_vara_int1(fid, &
            varid(iname,iopen),1,size(var),var)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    else
       ncStatus = nf_put_var1_int1(fid, &
            varid(iname,iopen), record_no(iname,iopen),var(1))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeInteger1Byte1D
  !--
  subroutine replaceInteger1Byte1D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    integer (kind=1), dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then   !writes all records
       ncStatus = nf_put_vara_int1(fid,vid,start_rec,size(var),var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
    else                              !writes one record a time
       ncStatus = nf_put_var1_int1(fid,vid,start_rec,var(1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
    endif
  end subroutine replaceInteger1Byte1D
  !--
  subroutine appendInteger1Byte1D(fid, var, varName, append)
    integer, intent(in) :: fid
    integer (kind=1), dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceInteger1Byte1D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendInteger1Byte1D
  !--
  subroutine writeInteger1Byte2D(fid, var, varName, varLenName, &
       varUnit, varLongName, append, unlimited_dim)
    integer (kind=1), dimension(:,:), intent (in) :: var
    character (len=*), intent (in), optional ::&
         varName, varLongName, varUnit
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    logical, intent (in), optional :: append
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer:: dimid1, dimid2, iopen
    integer, dimension(2) :: myVarDim

    call configureName(fid,iopen,define,varName=varName)

    if (present(varLenName)) then
       if (size(varLenName) /= 2) then
          print *,'err[ncdf_module::writeInteger2D]: ',&
               ' 2-element array required.'
          call exit(1)
       endif
       myvarLenName1 = varLenName(1)
       myvarLenName2 = varLenName(2)
    else
       myvarLenName1 = trim(tmpVarName(iname))//'Len1'
       myvarLenName2 = trim(tmpVarName(iname))//'Len2'
    endif

    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(unlimited_dim)) then
          select case (unlimited_dim)
          case (1)
             myvarDim(1) = unlimited_dim_id(iopen)
             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var(1,:)),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,tmpVarName(iname))
             endif
             myvarDim(2) = dimid2
          case (2)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var(:,1)),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,tmpVarName(iname))
             endif
             myvarDim(1) = dimid1
             myvarDim(2) = unlimited_dim_id(iopen)
          case default
             print *,'err[ncdf_module::writeInteger2D]: ',&
                  ' unlimited_dim for 2D integer: 1 or 2'
             call exit(1)
          end select
       else
          ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
               size(var(:,1)),dimid1)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName1),dimid1)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName1)
          endif
          myvarDim(1) = dimid1

          ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
               size(var(1,:)),dimid2)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName2),dimid2)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myvarLenName2)
          endif
          myvarDim(2) = dimid2
       endif
       ncStatus = nf_def_var(fid, trim(tmpVarName(iname)),nf_byte, &
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varLongName)
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),'units', &
               len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varUnit)
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       do irec = 1, size(var(1,:))
          ncStatus = nf_put_vara_int1(fid,varid(iname,iopen),&
               (/1,irec/),(/size(var(:,irec)),1/),var(:,irec))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,tmpVarName(iname))
       enddo
    else
       ncStatus = nf_put_vara_int1(fid,varid(iname,iopen),&
            (/1,record_no(iname,iopen)/),(/size(var,1),1/),var(:,1))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,trim(tmpVarName(iname)))
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeInteger1Byte2D
  !--
  subroutine replaceInteger1Byte2D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    integer (kind=1), dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then   !writes all records
       do irec = 1, size(var,2)
          ncStatus = nf_put_vara_int1(fid,vid,(/1,start_rec+irec-1/), &
               (/size(var,1),1/),var(:,irec))
          if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
       enddo
    else                              !writes one record a time
       ncStatus = nf_put_vara_int1(fid,vid,(/1,start_rec/), &
            (/size(var,1),1/),var(:,1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
    endif
  end subroutine replaceInteger1Byte2D
  !--
  subroutine appendInteger1Byte2D(fid, var, varName, append)
    integer, intent(in) :: fid
    integer (kind=1), dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceInteger1Byte2D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendInteger1Byte2D
  !--
  subroutine writeInteger1Byte3D(fid, var, varName, varLenName, &
       varUnit, varLongName, append, unlimited_dim)
    integer (kind=1), dimension(:,:,:), intent (in) :: var
    character (len=*), intent (in), optional ::&
         varName, varLongName, varUnit
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    logical, intent (in), optional :: append
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer :: dimid1, dimid2, dimid3, iopen
    integer, dimension(3) :: myVarDim

    call configureName(fid,iopen,&
         define,varName=varName)

    if (present(varLenName)) then
       myvarLenName1 = varLenName(1)
       myvarLenName2 = varLenName(2)
       myvarLenName3 = varLenName(3)
    else
       myvarLenName1 = trim(tmpVarName(iname))//'Len1'
       myvarLenName2 = trim(tmpVarName(iname))//'Len2'
       myvarLenName3 = trim(tmpVarName(iname))//'Len3'
    endif

    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(unlimited_dim)) then
          select case (unlimited_dim)
          case (1)
             myvarDim(1) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName3)
             endif
             myvarDim(3) = dimid3
          case (2)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             myvarDim(2) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName3)
             endif
             myvarDim(3) = dimid3
          case (3)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             myvarDim(3) = unlimited_dim_id(iopen)
          case default
             print *,'err[ncdf_module::writeInteger3D]: ',&
                  ' unlimited_dim for 2D integer: 1 or 2'
             call exit(1)
          end select
       else
          ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
               size(var,1),dimid1)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName1),dimid1)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName1)
          endif
          myvarDim(1) = dimid1

          ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
               size(var,2),dimid2)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName2),dimid2)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName2)
          endif
          myvarDim(2) = dimid2

          ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
               size(var,3),dimid3)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName3),dimid3)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName3)
          endif
          myvarDim(3) = dimid3
       endif

       ncStatus = nf_def_var(fid, trim(tmpVarName(iname)), nf_byte, &
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varLongName)
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'units', len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varUnit)
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       do j = 1, size(var,3)
          do i = 1, size(var,2)
             ncStatus = nf_put_vara_int1(fid,varid(iname,iopen),&
                  (/1,i,j/),(/size(var(:,i,j)),1,1/),var(:,i,j))
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,tmpVarName(iname))
          enddo
       enddo
    else
       do i = 1, size(var,2)
          ncStatus = nf_put_vara_int1(fid,varid(iname,iopen),&
               (/1,i,record_no(iname,iopen)/),&
               (/size(var(:,i,1)),1,1/),&
               var(:,i,1))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,tmpVarName(iname))
       enddo
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeInteger1Byte3D
  !--
  subroutine replaceInteger1Byte3D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    integer (kind=1), dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then   !writes all records
       do j = 1, size(var,3)
          do i = 1, size(var,2)
             ncStatus = nf_put_vara_int1(fid,vid,(/1,i,start_rec+j-1/), &
                  (/size(var,1),1,1/),var(:,i,j))
             if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
          enddo
       enddo
    else                              !writes one record a time
       do i = 1, size(var,2)
          ncStatus = nf_put_vara_int1(fid,vid, &
               (/1,i,start_rec/), &
               (/size(var,1),1,1/),var(:,i,1))
          if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
       enddo
    endif
  end subroutine replaceInteger1Byte3D
  !--
  subroutine appendInteger1Byte3D(fid, var, varName, append)
    integer, intent(in) :: fid
    integer (kind=1), dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceInteger1Byte3D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendInteger1Byte3D
  !--
  subroutine writeInteger1Byte4D(fid, var, varName, varLenName, &
       varUnit, varLongName, append, unlimited_dim)
    integer*1, dimension(:,:,:,:), intent (in) :: var
    character (len=*), intent (in), optional ::&
         varName, varLongName, varUnit
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    logical, intent (in), optional :: append
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer :: dimid1, dimid2, dimid3, dimid4, iopen
    integer, dimension(4) :: myVarDim

    call configureName(fid,iopen,&
         define,varName=varName)

    if (present(varLenName)) then
       myvarLenName1 = varLenName(1)
       myvarLenName2 = varLenName(2)
       myvarLenName3 = varLenName(3)
       myvarLenName4 = varLenName(4)
    else
       myvarLenName1 = trim(tmpVarName(iname))//'Len1'
       myvarLenName2 = trim(tmpVarName(iname))//'Len2'
       myvarLenName3 = trim(tmpVarName(iname))//'Len3'
       myvarLenName4 = trim(tmpVarName(iname))//'Len4'
    endif

    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(unlimited_dim)) then
          select case (unlimited_dim)
          case (1)
             myvarDim(1) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName3)
             endif
             myvarDim(3) = dimid3

             ncStatus = nf_def_dim(fid, trim(myvarLenName4),&
                  size(var,4),dimid4)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName4),dimid4)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName4)
             endif
             myvarDim(4) = dimid4
          case (2)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             myvarDim(2) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName3)
             endif
             myvarDim(3) = dimid3

             ncStatus = nf_def_dim(fid, trim(myvarLenName4),&
                  size(var,4),dimid4)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName4),dimid4)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName4)
             endif
             myvarDim(4) = dimid4
          case (3)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             myvarDim(3) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName4),&
                  size(var,4),dimid4)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName4),dimid4)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName4)
             endif
             myvarDim(4) = dimid4
          case (4)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName3)
             endif
             myvarDim(3) = dimid3

             myvarDim(4) = unlimited_dim_id(iopen)
          case default
             print *,'err[ncdf_module::writeInteger1Byte4D]: ',&
                  ' unlimited_dim for 2D integer: 1 or 2'
             call exit(1)
          end select
       else
          ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
               size(var,1),dimid1)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName1),dimid1)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName1)
          endif
          myvarDim(1) = dimid1

          ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
               size(var,2),dimid2)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName2),dimid2)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName2)
          endif
          myvarDim(2) = dimid2

          ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
               size(var,3),dimid3)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName3),dimid3)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName3)
          endif
          myvarDim(3) = dimid3

          ncStatus = nf_def_dim(fid, trim(myvarLenName4),&
               size(var,4),dimid4)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName4),dimid4)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName4)
          endif
          myvarDim(4) = dimid4
       endif

       ncStatus = nf_def_var(fid, trim(tmpVarName(iname)), nf_byte, &
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varLongName)
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'units', len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varUnit)
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       do k = 1, size(var,4)
          do j = 1, size(var,3)
             do i = 1, size(var,2)
                ncStatus = nf_put_vara_int1(fid,varid(iname,iopen),&
                     (/1,i,j,k/),(/size(var(:,i,j,k)),1,1,1/),&
                     var(:,i,j,k))
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,tmpVarName(iname))
             enddo
          enddo
       enddo
    else
       do j = 1, size(var,3)
          do i = 1, size(var,2)
             ncStatus = nf_put_vara_int1(fid,varid(iname,iopen),&
                  (/1,i,j,record_no(iname,iopen)/),&
                  (/size(var(:,i,j,1)),1,1,1/),&
                  var(:,i,j,1))
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,tmpVarName(iname))
          enddo
       enddo
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeInteger1Byte4D
  !--
  subroutine replaceInteger1Byte4D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    integer (kind=1), dimension(:,:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then   !writes all records
       do k = 1, size(var,4)
          do j = 1, size(var,3)
             do i = 1, size(var,2)
                ncStatus = nf_put_vara_int1(fid,vid,(/1,i,j,start_rec+k-1/), &
                     (/size(var,1),1,1,1/),var(:,i,j,k))
                if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
             enddo
          enddo
       enddo
    else                              !writes one record a time
       do j = 1, size(var,3)
          do i = 1, size(var,2)
             ncStatus = nf_put_vara_int1(fid,vid, &
                  (/1,i,j,start_rec/), &
                  (/size(var,1),1,1,1/),var(:,i,j,1))
             if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
          enddo
       enddo
    endif
  end subroutine replaceInteger1Byte4D
  !--
  subroutine appendInteger1Byte4D(fid, var, varName, append)
    integer, intent(in) :: fid
    integer (kind=1), dimension(:,:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceInteger1Byte4D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendInteger1Byte4D
  !--
  subroutine writeInteger2Byte1D(fid, var, varName, varLenName, &
       varUnit, varLongName, append, unlimited_dim)
    integer (kind=2), dimension(:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varLenName, varUnit
    logical, intent (in), optional :: append
    integer, intent (in) :: fid
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer :: dimid, iopen
    integer, dimension(1) :: myVarDim

    call configureName(fid,iopen,define,varName=varName)

    if (present(varLenName)) then
       myvarLenName = trim(varLenName)
    else
       myvarLenName = trim(tmpVarName(iname))//'Len'
    endif

    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(unlimited_dim)) then
          myvarDim(1) = unlimited_dim_id(iopen)
       else
          if (trim(myVarLenName) == &
               trim(myUnlimited_dim_name(iopen))) then
             myvarDim(1) = unlimited_dim_id(iopen)
          else
             ncStatus = nf_def_dim(fid,myVarLenName,&
                  size(var),dimid)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName),dimid)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,tmpVarName(iname))
             endif
             myvarDim(1) = dimid
          endif
       endif

       ncStatus = nf_def_var(fid, trim(tmpVarName(iname)), nf_short, &
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,tmpVarName(iname))
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),'units', &
               len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,tmpVarName(iname))
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       ncStatus = nf_put_vara_int2(fid, &
            varid(iname,iopen),1,size(var),var)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    else
       ncStatus = nf_put_var1_int2(fid, &
            varid(iname,iopen), record_no(iname,iopen),var(1))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeInteger2Byte1D
  !--
  subroutine replaceInteger2Byte1D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    integer (kind=2), dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then   !writes all records
       ncStatus = nf_put_vara_int2(fid,vid,start_rec,size(var),var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
    else                              !writes one record a time
       ncStatus = nf_put_var1_int2(fid,vid,start_rec,var(1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
    endif
  end subroutine replaceInteger2Byte1D
  !--
  subroutine appendInteger2Byte1D(fid, var, varName, append)
    integer, intent(in) :: fid
    integer (kind=2), dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceInteger2Byte1D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendInteger2Byte1D
  !--
  subroutine writeInteger2Byte2D(fid, var, varName, varLenName, &
       varUnit, varLongName, append, unlimited_dim)
    integer (kind=2), dimension(:,:), intent (in) :: var
    character (len=*), intent (in), optional ::&
         varName, varLongName, varUnit
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    logical, intent (in), optional :: append
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer:: dimid1, dimid2, iopen
    integer, dimension(2) :: myVarDim

    call configureName(fid,iopen,define,varName=varName)

    if (present(varLenName)) then
       if (size(varLenName) /= 2) then
          print *,'err[ncdf_module::writeInteger2D]: ',&
               ' 2-element array required.'
          call exit(1)
       endif
       myvarLenName1 = varLenName(1)
       myvarLenName2 = varLenName(2)
    else
       myvarLenName1 = trim(tmpVarName(iname))//'Len1'
       myvarLenName2 = trim(tmpVarName(iname))//'Len2'
    endif

    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(unlimited_dim)) then
          select case (unlimited_dim)
          case (1)
             myvarDim(1) = unlimited_dim_id(iopen)
             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var(1,:)),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,tmpVarName(iname))
             endif
             myvarDim(2) = dimid2
          case (2)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var(:,1)),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,tmpVarName(iname))
             endif
             myvarDim(1) = dimid1
             myvarDim(2) = unlimited_dim_id(iopen)
          case default
             print *,'err[ncdf_module::writeInteger2D]: ',&
                  ' unlimited_dim for 2D integer: 1 or 2'
             call exit(1)
          end select
       else
          ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
               size(var(:,1)),dimid1)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName1),dimid1)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName1)
          endif
          myvarDim(1) = dimid1

          ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
               size(var(1,:)),dimid2)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName2),dimid2)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myvarLenName2)
          endif
          myvarDim(2) = dimid2
       endif
       ncStatus = nf_def_var(fid, trim(tmpVarName(iname)),nf_short, &
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varLongName)
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),'units', &
               len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varUnit)
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       do irec = 1, size(var(1,:))
          ncStatus = nf_put_vara_int2(fid,varid(iname,iopen),&
               (/1,irec/),(/size(var(:,irec)),1/),var(:,irec))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,tmpVarName(iname))
       enddo
    else
       ncStatus = nf_put_vara_int2(fid,varid(iname,iopen),&
            (/1,record_no(iname,iopen)/),(/size(var,1),1/),var(:,1))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,trim(tmpVarName(iname)))
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeInteger2Byte2D
  !--
  subroutine replaceInteger2Byte2D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    integer (kind=2), dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then   !writes all records
       do irec = 1, size(var,2)
          ncStatus = nf_put_vara_int2(fid,vid,(/1,start_rec+irec-1/), &
               (/size(var,1),1/),var(:,irec))
          if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
       enddo
    else                              !writes one record a time
       ncStatus = nf_put_vara_int2(fid,vid,(/1,start_rec/), &
            (/size(var,1),1/),var(:,1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
    endif
  end subroutine replaceInteger2Byte2D
  !--
  subroutine appendInteger2Byte2D(fid, var, varName, append)
    integer, intent(in) :: fid
    integer (kind=2), dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceInteger2Byte2D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendInteger2Byte2D
  !--
  subroutine writeInteger2Byte3D(fid, var, varName, varLenName, &
       varUnit, varLongName, append, unlimited_dim)
    integer (kind=2), dimension(:,:,:), intent (in) :: var
    character (len=*), intent (in), optional ::&
         varName, varLongName, varUnit
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    logical, intent (in), optional :: append
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer :: dimid1, dimid2, dimid3, iopen
    integer, dimension(3) :: myVarDim

    call configureName(fid,iopen,&
         define,varName=varName)

    if (present(varLenName)) then
       myvarLenName1 = varLenName(1)
       myvarLenName2 = varLenName(2)
       myvarLenName3 = varLenName(3)
    else
       myvarLenName1 = trim(tmpVarName(iname))//'Len1'
       myvarLenName2 = trim(tmpVarName(iname))//'Len2'
       myvarLenName3 = trim(tmpVarName(iname))//'Len3'
    endif

    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(unlimited_dim)) then
          select case (unlimited_dim)
          case (1)
             myvarDim(1) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName3)
             endif
             myvarDim(3) = dimid3
          case (2)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             myvarDim(2) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName3)
             endif
             myvarDim(3) = dimid3
          case (3)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             myvarDim(3) = unlimited_dim_id(iopen)
          case default
             print *,'err[ncdf_module::writeInteger3D]: ',&
                  ' unlimited_dim for 2D integer: 1 or 2'
             call exit(1)
          end select
       else
          ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
               size(var,1),dimid1)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName1),dimid1)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName1)
          endif
          myvarDim(1) = dimid1

          ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
               size(var,2),dimid2)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName2),dimid2)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName2)
          endif
          myvarDim(2) = dimid2

          ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
               size(var,3),dimid3)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName3),dimid3)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName3)
          endif
          myvarDim(3) = dimid3
       endif

       ncStatus = nf_def_var(fid, trim(tmpVarName(iname)), nf_short, &
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varLongName)
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'units', len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varUnit)
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       do j = 1, size(var,3)
          do i = 1, size(var,2)
             ncStatus = nf_put_vara_int2(fid,varid(iname,iopen),&
                  (/1,i,j/),(/size(var(:,i,j)),1,1/),var(:,i,j))
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,tmpVarName(iname))
          enddo
       enddo
    else
       do i = 1, size(var,2)
          ncStatus = nf_put_vara_int2(fid,varid(iname,iopen),&
               (/1,i,record_no(iname,iopen)/),&
               (/size(var(:,i,1)),1,1/),&
               var(:,i,1))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,tmpVarName(iname))
       enddo
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeInteger2Byte3D
  !--
  subroutine replaceInteger2Byte3D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    integer (kind=2), dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then   !writes all records
       do j = 1, size(var,3)
          do i = 1, size(var,2)
             ncStatus = nf_put_vara_int2(fid,vid,(/1,i,start_rec+j-1/), &
                  (/size(var,1),1,1/),var(:,i,j))
             if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
          enddo
       enddo
    else                              !writes one record a time
       do i = 1, size(var,2)
          ncStatus = nf_put_vara_int2(fid,vid, &
               (/1,i,start_rec/), &
               (/size(var,1),1,1/),var(:,i,1))
          if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
       enddo
    endif
  end subroutine replaceInteger2Byte3D
  !--
  subroutine appendInteger2Byte3D(fid, var, varName, append)
    integer, intent(in) :: fid
    integer (kind=2), dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceInteger2Byte3D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendInteger2Byte3D
  !--
  subroutine writeInteger2Byte4D(fid, var, varName, varLenName, &
       varUnit, varLongName, append, unlimited_dim)
    integer*2, dimension(:,:,:,:), intent (in) :: var
    character (len=*), intent (in), optional ::&
         varName, varLongName, varUnit
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    logical, intent (in), optional :: append
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer :: dimid1, dimid2, dimid3, dimid4, iopen
    integer, dimension(4) :: myVarDim

    call configureName(fid,iopen,&
         define,varName=varName)

    if (present(varLenName)) then
       myvarLenName1 = varLenName(1)
       myvarLenName2 = varLenName(2)
       myvarLenName3 = varLenName(3)
       myvarLenName4 = varLenName(4)
    else
       myvarLenName1 = trim(tmpVarName(iname))//'Len1'
       myvarLenName2 = trim(tmpVarName(iname))//'Len2'
       myvarLenName3 = trim(tmpVarName(iname))//'Len3'
       myvarLenName4 = trim(tmpVarName(iname))//'Len4'
    endif

    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(unlimited_dim)) then
          select case (unlimited_dim)
          case (1)
             myvarDim(1) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName3)
             endif
             myvarDim(3) = dimid3

             ncStatus = nf_def_dim(fid, trim(myvarLenName4),&
                  size(var,4),dimid4)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName4),dimid4)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName4)
             endif
             myvarDim(4) = dimid4
          case (2)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             myvarDim(2) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName3)
             endif
             myvarDim(3) = dimid3

             ncStatus = nf_def_dim(fid, trim(myvarLenName4),&
                  size(var,4),dimid4)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName4),dimid4)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName4)
             endif
             myvarDim(4) = dimid4
          case (3)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             myvarDim(3) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName4),&
                  size(var,4),dimid4)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName4),dimid4)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName4)
             endif
             myvarDim(4) = dimid4
          case (4)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName3)
             endif
             myvarDim(3) = dimid3

             myvarDim(4) = unlimited_dim_id(iopen)
          case default
             print *,'err[ncdf_module::writeInteger2Byte4D]: ',&
                  ' unlimited_dim for 2D integer: 1 or 2'
             call exit(1)
          end select
       else
          ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
               size(var,1),dimid1)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName1),dimid1)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName1)
          endif
          myvarDim(1) = dimid1

          ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
               size(var,2),dimid2)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName2),dimid2)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName2)
          endif
          myvarDim(2) = dimid2

          ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
               size(var,3),dimid3)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName3),dimid3)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName3)
          endif
          myvarDim(3) = dimid3

          ncStatus = nf_def_dim(fid, trim(myvarLenName4),&
               size(var,4),dimid4)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName4),dimid4)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName4)
          endif
          myvarDim(4) = dimid4
       endif

       ncStatus = nf_def_var(fid, trim(tmpVarName(iname)), nf_short, &
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varLongName)
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'units', len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varUnit)
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       do k = 1, size(var,4)
          do j = 1, size(var,3)
             do i = 1, size(var,2)
                ncStatus = nf_put_vara_int2(fid,varid(iname,iopen),&
                     (/1,i,j,k/),(/size(var(:,i,j,k)),1,1,1/),&
                     var(:,i,j,k))
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,tmpVarName(iname))
             enddo
          enddo
       enddo
    else
       do j = 1, size(var,3)
          do i = 1, size(var,2)
             ncStatus = nf_put_vara_int2(fid,varid(iname,iopen),&
                  (/1,i,j,record_no(iname,iopen)/),&
                  (/size(var(:,i,j,1)),1,1,1/),&
                  var(:,i,j,1))
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,tmpVarName(iname))
          enddo
       enddo
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeInteger2Byte4D
  !--
  subroutine replaceInteger2Byte4D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    integer (kind=2), dimension(:,:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then   !writes all records
       do k = 1, size(var,4)
          do j = 1, size(var,3)
             do i = 1, size(var,2)
                ncStatus = nf_put_vara_int2(fid,vid,(/1,i,j,start_rec+k-1/), &
                     (/size(var,1),1,1,1/),var(:,i,j,k))
                if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
             enddo
          enddo
       enddo
    else                              !writes one record a time
       do j = 1, size(var,3)
          do i = 1, size(var,2)
             ncStatus = nf_put_vara_int2(fid,vid, &
                  (/1,i,j,start_rec/), &
                  (/size(var,1),1,1,1/),var(:,i,j,1))
             if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
          enddo
       enddo
    endif
  end subroutine replaceInteger2Byte4D
  !--
  subroutine appendInteger2Byte4D(fid, var, varName, append)
    integer, intent(in) :: fid
    integer (kind=2), dimension(:,:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceInteger2Byte4D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendInteger2Byte4D
  !--
  subroutine writeFloat1D(fid, var, varName, varLenName, &
       varUnit, varLongName, append, unlimited_dim)
    real*4, dimension(:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varLenName, varUnit
    integer, intent (in) :: fid
    logical, intent (in), optional :: append
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer :: dimid, iopen
    integer, dimension(1) :: myVarDim

    call configureName(fid,iopen,&
         define,varName=varName)

    if (present(varLenName)) then
       myvarLenName = trim(varLenName)
    else
       myvarLenName = trim(tmpVarName(iname))//'Len'
    endif

    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(unlimited_dim)) then
          myvarDim(1) = unlimited_dim_id(iopen)
       else
          if (trim(myVarLenName) == &
               trim(myUnlimited_dim_name(iopen))) then
             myvarDim(1) = unlimited_dim_id(iopen)
          else
             ncStatus = nf_def_dim(fid,trim(myVarLenName),&
                  size(var),dimid)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName),dimid)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName)
             endif
             myvarDim(1) = dimid
          endif
       endif

       ncStatus = nf_def_var(fid,trim(tmpVarName(iname)),nf_float,&
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,trim(tmpVarName(iname)))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,trim(varLongName))
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'units', len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,trim(varUnit))
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       ncStatus = nf_put_vara_real(fid,varid(iname,iopen),&
            1,size(var),var)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    else
       ncStatus = nf_put_var1_real(fid,varid(iname,iopen),&
            record_no(iname,iopen),var(1))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeFloat1D
  !--
  subroutine replaceFloat1D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    real*4, dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then   !writes all records
       ncStatus = nf_put_vara_real(fid,vid,start_rec,size(var),var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
    else                              !writes one record a time
       ncStatus = nf_put_var1_real(fid,vid,start_rec,var(1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
    endif
  end subroutine replaceFloat1D
  !--
  subroutine appendFloat1D(fid, var, varName, append)
    integer, intent(in) :: fid
    real*4, dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceFloat1D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendFloat1D
  !--
  subroutine writeFloat2D(fid, var, varName, varLenName, &
       varUnit, varLongName, append, unlimited_dim)
    real*4, dimension(:,:), intent (in) :: var
    character (len=*), intent (in), optional ::&
         varName, varLongName, varUnit
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    logical, intent (in), optional :: append
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer :: dimid1, dimid2, iopen
    integer, dimension(2) :: myVarDim

    call configureName(fid,iopen,&
         define,varName=varName)

    if (present(varLenName)) then
       if (size(varLenName) /= 2) then
          print *,'err[ncdf_module::writeFloat2D]: ',&
               ' 2-element array required.'
          call exit(1)
       endif
       myvarLenName1 = varLenName(1)
       myvarLenName2 = varLenName(2)
    else
       myvarLenName1 = trim(tmpVarName(iname))//'Len1'
       myvarLenName2 = trim(tmpVarName(iname))//'Len2'
    endif

    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(unlimited_dim)) then
          select case (unlimited_dim)
          case (1)
             myvarDim(1) = unlimited_dim_id(iopen)
             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var(1,:)),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myvarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,tmpVarName(iname))
             endif
             myvarDim(2) = dimid2
          case (2)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var(:,1)),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myvarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,tmpVarName(iname))
             endif
             myvarDim(1) = dimid1
             myvarDim(2) = unlimited_dim_id(iopen)
          case default
             print *,'err[ncdf_module::writeFloat2D]: ',&
                  ' unlimited_dim for 2D float: 1 or 2'
             call exit(1)
          end select
       else
          ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
               size(var(:,1)),dimid1)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName1),dimid1)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,tmpVarName(iname))
          endif
          myvarDim(1) = dimid1

          ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
               size(var(1,:)),dimid2)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName2),dimid2)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,tmpVarName(iname))
          endif
          myvarDim(2) = dimid2
       endif
       ncStatus = nf_def_var(fid,trim(tmpVarName(iname)),nf_float,&
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,trim(tmpVarName(iname)))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,trim(varLongName))
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'units', len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,trim(varUnit))
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       do irec = 1, size(var(1,:))
          ncStatus = nf_put_vara_real(fid,varid(iname,iopen),&
               (/1,irec/),(/size(var(:,irec)),1/),var(:,irec))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,tmpVarName(iname))
       enddo
    else
       ncStatus = nf_put_vara_real(fid,varid(iname,iopen),&
            (/1,record_no(iname,iopen)/),(/size(var,1),1/),var(:,1))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeFloat2D
  !--
  subroutine replaceFloat2D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    real*4, dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then   !writes all records
       do irec = 1, size(var,2)
          ncStatus = nf_put_vara_real(fid,vid,(/1,start_rec+irec-1/), &
               (/size(var,1),1/),var(:,irec))
          if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
       enddo
    else                              !writes one record a time
       ncStatus = nf_put_vara_real(fid,vid,(/1,start_rec/), &
            (/size(var,1),1/),var(:,1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
    endif
  end subroutine replaceFloat2D
  !--
  subroutine appendFloat2D(fid, var, varName, append)
    integer, intent(in) :: fid
    real*4, dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceFloat2D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendFloat2D
  !--
  subroutine writeFloat3D(fid, var, varName, varLenName, &
       varUnit, varLongName, append, unlimited_dim)
    real*4, dimension(:,:,:), intent (in) :: var
    character (len=*), intent (in), optional ::&
         varName, varLongName, varUnit
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    logical, intent (in), optional :: append
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer :: dimid1, dimid2, dimid3, iopen
    integer, dimension(3) :: myVarDim

    call configureName(fid,iopen,&
         define,varName=varName)

    if (present(varLenName)) then
       myvarLenName1 = varLenName(1)
       myvarLenName2 = varLenName(2)
       myvarLenName3 = varLenName(3)
    else
       myvarLenName1 = trim(tmpVarName(iname))//'Len1'
       myvarLenName2 = trim(tmpVarName(iname))//'Len2'
       myvarLenName3 = trim(tmpVarName(iname))//'Len3'
    endif

    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(unlimited_dim)) then
          select case (unlimited_dim)
          case (1)
             myvarDim(1) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName3)
             endif
             myvarDim(3) = dimid3
          case (2)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             myvarDim(2) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName3)
             endif
             myvarDim(3) = dimid3
          case (3)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             myvarDim(3) = unlimited_dim_id(iopen)
          case default
             print *,'err[ndcf_module::writeFloat3D]: ',&
                  ' unlimited_dim for 2D integer: 1 or 2'
             call exit(1)
          end select
       else
          ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
               size(var,1),dimid1)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName1),dimid1)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName1)
          endif
          myvarDim(1) = dimid1

          ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
               size(var,2),dimid2)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName2),dimid2)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName2)
          endif
          myvarDim(2) = dimid2

          ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
               size(var,3),dimid3)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName3),dimid3)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName3)
          endif
          myvarDim(3) = dimid3
       endif

       ncStatus = nf_def_var(fid,trim(tmpVarName(iname)),nf_float, &
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varLongName)
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'units', len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varUnit)
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       do j = 1, size(var,3)
          do i = 1, size(var,2)
             ncStatus = nf_put_vara_real(fid,varid(iname,iopen),&
                  (/1,i,j/),(/size(var(:,i,j)),1,1/),var(:,i,j))
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,tmpVarName(iname))
          enddo
       enddo
    else
       do i = 1, size(var,2)
          ncStatus = nf_put_vara_real(fid,varid(iname,iopen),&
               (/1,i,record_no(iname,iopen)/),&
               (/size(var(:,i,1)),1,1/),&
               var(:,i,1))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,tmpVarName(iname))
       enddo
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeFloat3D
  !--
  subroutine replaceFloat3D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    real*4, dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then   !writes all records
       do j = 1, size(var,3)
          do i = 1, size(var,2)
             ncStatus = nf_put_vara_real(fid,vid,(/1,i,start_rec+j-1/), &
                  (/size(var,1),1,1/),var(:,i,j))
             if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
          enddo
       enddo
    else                              !writes one record a time
       do i = 1, size(var,2)
          ncStatus = nf_put_vara_real(fid,vid, &
               (/1,i,start_rec/), &
               (/size(var,1),1,1/),var(:,i,1))
          if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
       enddo
    endif
  end subroutine replaceFloat3D
  !--
  subroutine appendFloat3D(fid, var, varName, append)
    integer, intent(in) :: fid
    real*4, dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceFloat3D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendFloat3D
  !--
  subroutine writeFloat4D(fid, var, varName, varLenName, &
       varUnit, varLongName, append, unlimited_dim)
    real*4, dimension(:,:,:,:), intent (in) :: var
    character (len=*), intent (in), optional ::&
         varName, varLongName, varUnit
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    logical, intent (in), optional :: append
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer :: dimid1, dimid2, dimid3, dimid4, iopen
    integer, dimension(4) :: myVarDim

    call configureName(fid,iopen,&
         define,varName=varName)

    if (present(varLenName)) then
       myvarLenName1 = varLenName(1)
       myvarLenName2 = varLenName(2)
       myvarLenName3 = varLenName(3)
       myvarLenName4 = varLenName(4)
    else
       myvarLenName1 = trim(tmpVarName(iname))//'Len1'
       myvarLenName2 = trim(tmpVarName(iname))//'Len2'
       myvarLenName3 = trim(tmpVarName(iname))//'Len3'
       myvarLenName4 = trim(tmpVarName(iname))//'Len4'
    endif

    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(unlimited_dim)) then
          select case (unlimited_dim)
          case (1)
             myvarDim(1) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName3)
             endif
             myvarDim(3) = dimid3

             ncStatus = nf_def_dim(fid, trim(myvarLenName4),&
                  size(var,4),dimid4)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName4),dimid4)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName4)
             endif
             myvarDim(4) = dimid4
          case (2)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             myvarDim(2) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName3)
             endif
             myvarDim(3) = dimid3

             ncStatus = nf_def_dim(fid, trim(myvarLenName4),&
                  size(var,4),dimid4)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName4),dimid4)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName4)
             endif
             myvarDim(4) = dimid4
          case (3)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             myvarDim(3) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName4),&
                  size(var,4),dimid4)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName4),dimid4)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName4)
             endif
             myvarDim(4) = dimid4
          case (4)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName3)
             endif
             myvarDim(3) = dimid3

             myvarDim(4) = unlimited_dim_id(iopen)
          case default
             print *,'err[ncdf_module::writeFloat4D]: ',&
                  ' unlimited_dim for 2D integer: 1 or 2'
             call exit(1)
          end select
       else
          ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
               size(var,1),dimid1)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName1),dimid1)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName1)
          endif
          myvarDim(1) = dimid1

          ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
               size(var,2),dimid2)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName2),dimid2)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName2)
          endif
          myvarDim(2) = dimid2

          ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
               size(var,3),dimid3)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName3),dimid3)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName3)
          endif
          myvarDim(3) = dimid3

          ncStatus = nf_def_dim(fid, trim(myvarLenName4),&
               size(var,4),dimid4)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName4),dimid4)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName4)
          endif
          myvarDim(4) = dimid4
       endif

       ncStatus = nf_def_var(fid, trim(tmpVarName(iname)), nf_float, &
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varLongName)
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'units', len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varUnit)
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       do k = 1, size(var,4)
          do j = 1, size(var,3)
             do i = 1, size(var,2)
                ncStatus = nf_put_vara_real(fid,varid(iname,iopen),&
                     (/1,i,j,k/),(/size(var(:,i,j,k)),1,1,1/),&
                     var(:,i,j,k))
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,tmpVarName(iname))
             enddo
          enddo
       enddo
    else
       do j = 1, size(var,3)
          do i = 1, size(var,2)
             ncStatus = nf_put_vara_real(fid,varid(iname,iopen),&
                  (/1,i,j,record_no(iname,iopen)/),&
                  (/size(var(:,i,j,1)),1,1,1/),&
                  var(:,i,j,1))
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,tmpVarName(iname))
          enddo
       enddo
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeFloat4D
  !--
  subroutine replaceFloat4D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    real*4, dimension(:,:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then   !writes all records
       do k = 1, size(var,4)
          do j = 1, size(var,3)
             do i = 1, size(var,2)
                ncStatus = nf_put_vara_real(fid,vid,(/1,i,j,start_rec+k-1/), &
                     (/size(var,1),1,1,1/),var(:,i,j,k))
                if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
             enddo
          enddo
       enddo
    else                              !writes one record a time
       do j = 1, size(var,3)
          do i = 1, size(var,2)
             ncStatus = nf_put_vara_real(fid,vid, &
                  (/1,i,j,start_rec/), &
                  (/size(var,1),1,1,1/),var(:,i,j,1))
             if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
          enddo
       enddo
    endif
  end subroutine replaceFloat4D
  !--
  subroutine appendFloat4D(fid, var, varName, append)
    integer, intent(in) :: fid
    real*4, dimension(:,:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceFloat4D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendFloat4D
  !--
  subroutine writeDFloat1D(fid, var, varName, varLenName, &
       varUnit, varLongName, append, unlimited_dim)
    double precision, dimension(:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varLenName, varUnit
    integer, intent (in) :: fid
    logical, intent (in), optional :: append
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer :: dimid, iopen
    integer, dimension(1) :: myVarDim

    call configureName(fid,iopen,&
         define,varName=varName)

    if (present(varLenName)) then
       myvarLenName = trim(varLenName)
    else
       myvarLenName = trim(tmpVarName(iname))//'Len'
    endif

    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(unlimited_dim)) then
          myvarDim(1) = unlimited_dim_id(iopen)
       else
          if (trim(myVarLenName) == &
               trim(myUnlimited_dim_name(iopen))) then
             myvarDim(1) = unlimited_dim_id(iopen)
          else
             ncStatus = nf_def_dim(fid,trim(myVarLenName),&
                  size(var),dimid)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName),dimid)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName)
             endif
             myvarDim(1) = dimid
          endif
       endif

       ncStatus = nf_def_var(fid,trim(tmpVarName(iname)),nf_double, &
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varLongName)
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'units', len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varUnit)
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       ncStatus = nf_put_vara_double(fid,varid(iname,iopen),&
            1,size(var),var)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    else
       ncStatus = nf_put_var1_double(fid,varid(iname,iopen),&
            record_no(iname,iopen),var(1))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeDFloat1D
  !--
  subroutine replaceDFloat1D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    double precision, dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then   !writes all records
       ncStatus = nf_put_vara_double(fid,vid,start_rec,size(var),var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
    else                              !writes one record a time
       ncStatus = nf_put_var1_double(fid,vid,start_rec,var(1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
    endif
  end subroutine replaceDFloat1D
  !--
  subroutine appendDFloat1D(fid, var, varName, append)
    integer, intent(in) :: fid
    double precision, dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceDFloat1D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendDFloat1D
  !--
  subroutine writeDFloat2D(fid, var, varName, varLenName, &
       varUnit, varLongName, append, unlimited_dim)
    double precision, dimension(:,:), intent (in) :: var
    character (len=*), intent (in), optional ::&
         varName, varLongName, varUnit
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    logical, intent (in), optional :: append
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer :: dimid1, dimid2, iopen
    integer, dimension(2) :: myVarDim

    call configureName(fid,iopen,&
         define,varName=varName)

    if (present(varLenName)) then
       if (size(varLenName) /= 2) then
          print *,'err[ncdf_module::writeDFloat2D]: ',&
               ' 2-element array required.'
          call exit(1)
       endif
       myvarLenName1 = varLenName(1)
       myvarLenName2 = varLenName(2)
    else
       myvarLenName1 = trim(tmpVarName(iname))//'Len1'
       myvarLenName2 = trim(tmpVarName(iname))//'Len2'
    endif

    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(unlimited_dim)) then
          select case (unlimited_dim)
          case (1)
             myvarDim(1) = unlimited_dim_id(iopen)
             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var(1,:)),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2
          case (2)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var(:,1)),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1
             myvarDim(2) = unlimited_dim_id(iopen)
          case default
             print *,'err[ncdf_module::writeDFloat2D]: ',&
                  ' unlimited_dim for 2D float: 1 or 2'
             call exit(1)
          end select
       else
          ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
               size(var(:,1)),dimid1)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName1),dimid1)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName1)
          endif
          myvarDim(1) = dimid1

          ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
               size(var(1,:)),dimid2)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName2),dimid2)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName2)
          endif
          myvarDim(2) = dimid2
       endif
       ncStatus = nf_def_var(fid,trim(tmpVarName(iname)),nf_double, &
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,trim(tmpVarName(iname)))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varLongName)
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'units', len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varUnit)
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       do irec = 1, size(var(1,:))
          ncStatus = nf_put_vara_double(fid,varid(iname,iopen),&
               (/1,irec/),(/size(var(:,irec)),1/),var(:,irec))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,tmpVarName(iname))
       enddo
    else
       ncStatus = nf_put_vara_double(fid,varid(iname,iopen),&
            (/1,record_no(iname,iopen)/),(/size(var,1),1/),var(:,1))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeDFloat2D
  !--
  subroutine replaceDFloat2D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    double precision, dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then   !writes all records
       do irec = 1, size(var,2)
          ncStatus = nf_put_vara_double(fid,vid,(/1,start_rec+irec-1/), &
               (/size(var,1),1/),var(:,irec))
          if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
       enddo
    else                              !writes one record a time
       ncStatus = nf_put_vara_double(fid,vid,(/1,start_rec/), &
            (/size(var,1),1/),var(:,1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
    endif
  end subroutine replaceDFloat2D
  !--
  subroutine appendDFloat2D(fid, var, varName, append)
    integer, intent(in) :: fid
    double precision, dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceDFloat2D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendDFloat2D
  !--
  subroutine writeDFloat3D(fid, var, varName, varLenName, &
       varUnit, varLongName, append, unlimited_dim)
    double precision, dimension(:,:,:), intent (in) :: var
    character (len=*), intent (in), optional ::&
         varName, varLongName, varUnit
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    logical, intent (in), optional :: append
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer :: dimid1, dimid2, dimid3, iopen
    integer, dimension(3) :: myVarDim

    call configureName(fid,iopen,&
         define,varName=varName)

    if (present(varLenName)) then
       myvarLenName1 = varLenName(1)
       myvarLenName2 = varLenName(2)
       myvarLenName3 = varLenName(3)
    else
       myvarLenName1 = trim(tmpVarName(iname))//'Len1'
       myvarLenName2 = trim(tmpVarName(iname))//'Len2'
       myvarLenName3 = trim(tmpVarName(iname))//'Len3'
    endif

    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(unlimited_dim)) then
          select case (unlimited_dim)
          case (1)
             myvarDim(1) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName3)
             endif
             myvarDim(3) = dimid3
          case (2)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             myvarDim(2) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName3)
             endif
             myvarDim(3) = dimid3
          case (3)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             myvarDim(3) = unlimited_dim_id(iopen)
          case default
             print *,'err[ncdf_module::writeDFloat3D]: ',&
                  ' unlimited_dim for 2D integer: 1 or 2'
             call exit(1)
          end select
       else
          ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
               size(var,1),dimid1)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName1),dimid1)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName1)
          endif
          myvarDim(1) = dimid1

          ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
               size(var,2),dimid2)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName2),dimid2)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName2)
          endif
          myvarDim(2) = dimid2

          ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
               size(var,3),dimid3)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName3),dimid3)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName3)
          endif
          myvarDim(3) = dimid3
       endif

       ncStatus = nf_def_var(fid,trim(tmpVarName(iname)),nf_double, &
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varLongName)
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'units', len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varUnit)
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       do j = 1, size(var,3)
          do i = 1, size(var,2)
             ncStatus = nf_put_vara_double(fid,varid(iname,iopen),&
                  (/1,i,j/),(/size(var(:,i,j)),1,1/),var(:,i,j))
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,tmpVarName(iname))
          enddo
       enddo
    else
       do i = 1, size(var,2)
          ncStatus = nf_put_vara_double(fid,varid(iname,iopen),&
               (/1,i,record_no(iname,iopen)/),&
               (/size(var(:,i,1)),1,1/),&
               var(:,i,1))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,tmpVarName(iname))
       enddo
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeDFloat3D
  !--
  subroutine replaceDFloat3D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    double precision, dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then   !writes all records
       do j = 1, size(var,3)
          do i = 1, size(var,2)
             ncStatus = nf_put_vara_double(fid,vid,(/1,i,start_rec+j-1/), &
                  (/size(var,1),1,1/),var(:,i,j))
             if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
          enddo
       enddo
    else                              !writes one record a time
       do i = 1, size(var,2)
          ncStatus = nf_put_vara_double(fid,vid, &
               (/1,i,start_rec/), &
               (/size(var,1),1,1/),var(:,i,1))
          if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
       enddo
    endif
  end subroutine replaceDFloat3D
  !--
  subroutine appendDFloat3D(fid, var, varName, append)
    integer, intent(in) :: fid
    double precision, dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceDFloat3D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendDFloat3D
  !--
  subroutine writeDFloat4D(fid, var, varName, varLenName, &
       varUnit, varLongName, append, unlimited_dim)
    double precision, dimension(:,:,:,:), intent (in) :: var
    character (len=*), intent (in), optional ::&
         varName, varLongName, varUnit
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    logical, intent (in), optional :: append
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer :: dimid1, dimid2, dimid3, dimid4, iopen
    integer, dimension(4) :: myVarDim

    call configureName(fid,iopen,&
         define,varName=varName)

    if (present(varLenName)) then
       myvarLenName1 = varLenName(1)
       myvarLenName2 = varLenName(2)
       myvarLenName3 = varLenName(3)
       myvarLenName4 = varLenName(4)
    else
       myvarLenName1 = trim(tmpVarName(iname))//'Len1'
       myvarLenName2 = trim(tmpVarName(iname))//'Len2'
       myvarLenName3 = trim(tmpVarName(iname))//'Len3'
       myvarLenName4 = trim(tmpVarName(iname))//'Len4'
    endif

    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(unlimited_dim)) then
          select case (unlimited_dim)
          case (1)
             myvarDim(1) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName3)
             endif
             myvarDim(3) = dimid3

             ncStatus = nf_def_dim(fid, trim(myvarLenName4),&
                  size(var,4),dimid4)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName4),dimid4)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName4)
             endif
             myvarDim(4) = dimid4
          case (2)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             myvarDim(2) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName3)
             endif
             myvarDim(3) = dimid3

             ncStatus = nf_def_dim(fid, trim(myvarLenName4),&
                  size(var,4),dimid4)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName4),dimid4)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName4)
             endif
             myvarDim(4) = dimid4
          case (3)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             myvarDim(3) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myvarLenName4),&
                  size(var,4),dimid4)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName4),dimid4)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myvarLenName4)
             endif
             myvarDim(4) = dimid4
          case (4)
             ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
                  size(var,1),dimid1)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid1)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(1) = dimid1

             ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
                  size(var,2),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(2) = dimid2

             ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
                  size(var,3),dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName3)
             endif
             myvarDim(3) = dimid3

             myvarDim(4) = unlimited_dim_id(iopen)
          case default
             print *,'err[ncdf_module::writeDFloat4D]: ',&
                  ' unlimited_dim for 2D integer: 1 or 2'
             call exit(1)
          end select
       else
          ncStatus = nf_def_dim(fid, trim(myvarLenName1),&
               size(var,1),dimid1)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName1),dimid1)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName1)
          endif
          myvarDim(1) = dimid1

          ncStatus = nf_def_dim(fid, trim(myvarLenName2),&
               size(var,2),dimid2)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName2),dimid2)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName2)
          endif
          myvarDim(2) = dimid2

          ncStatus = nf_def_dim(fid, trim(myvarLenName3),&
               size(var,3),dimid3)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName3),dimid3)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName3)
          endif
          myvarDim(3) = dimid3

          ncStatus = nf_def_dim(fid, trim(myvarLenName4),&
               size(var,4),dimid4)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName4),dimid4)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName4)
          endif
          myvarDim(4) = dimid4
       endif

       ncStatus = nf_def_var(fid, trim(tmpVarName(iname)), nf_double, &
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varLongName)
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'units', len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varUnit)
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       do k = 1, size(var,4)
          do j = 1, size(var,3)
             do i = 1, size(var,2)
                ncStatus = nf_put_vara_double(fid,varid(iname,iopen),&
                     (/1,i,j,k/),(/size(var(:,i,j,k)),1,1,1/),&
                     var(:,i,j,k))
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,tmpVarName(iname))
             enddo
          enddo
       enddo
    else
       do j = 1, size(var,3)
          do i = 1, size(var,2)
             ncStatus = nf_put_vara_double(fid,varid(iname,iopen),&
                  (/1,i,j,record_no(iname,iopen)/),&
                  (/size(var(:,i,j,1)),1,1,1/),&
                  var(:,i,j,1))
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,tmpVarName(iname))
          enddo
       enddo
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeDFloat4D
  !--
  subroutine replaceDFloat4D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    double precision, dimension(:,:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then   !writes all records
       do k = 1, size(var,4)
          do j = 1, size(var,3)
             do i = 1, size(var,2)
                ncStatus = nf_put_vara_double(fid,vid,(/1,i,j,start_rec+k-1/), &
                     (/size(var,1),1,1,1/),var(:,i,j,k))
                if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
             enddo
          enddo
       enddo
    else                              !writes one record a time
       do j = 1, size(var,3)
          do i = 1, size(var,2)
             ncStatus = nf_put_vara_double(fid,vid, &
                  (/1,i,j,start_rec/), &
                  (/size(var,1),1,1,1/),var(:,i,j,1))
             if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
          enddo
       enddo
    endif
  end subroutine replaceDFloat4D
  !--
  subroutine appendDFloat4D(fid, var, varName, append)
    integer, intent(in) :: fid
    double precision, dimension(:,:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceDFloat4D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendDFloat4D
  !--
  subroutine writeString1D(fid, var, varName, strLenName,&
       varLenName, varUnit, varLongName, append, unlimited_dim)
    character (len=*), dimension(:), intent (in) :: var
    character (len=*), intent (in), optional ::&
         varName, varLongName, varUnit
    character (len=*), intent (in), optional :: varLenName
    character (len=*), intent (in), optional :: strLenName
    integer, intent (in) :: fid
    logical, intent (in), optional :: append
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer :: dimid1, dimid2, iopen
    integer, dimension(2) :: myVarDim

    call configureName(fid,iopen,&
         define,varName=varName)

    if (present(varLenName)) then
       myvarLenName = trim(varLenName)
    else
       myvarLenName = trim(tmpVarName(iname))//'Len'
    endif

    if (present(strLenName)) then
       myStrLenName = trim(strLenName)
    else
       myStrLenName = trim(tmpVarName(iname))//'StrLen'
    endif

    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       ncStatus = nf_def_dim(fid,&
            myStrLenName, len_trim(var(1)), dimid1)
       if (ncStatus == -42) then
          ncStatus = nf_inq_dimid(fid,&
               trim(myStrLenName),dimid1)
       else
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,myStrLenName)
       endif
       myvarDim(1) = dimid1

       if (present(unlimited_dim)) then
          myvarDim(2) = unlimited_dim_id(iopen)
       else
          if (trim(myVarLenName) == &
               trim(myUnlimited_dim_name(iopen))) then
             myvarDim(2) = unlimited_dim_id(iopen)
          else
             ncStatus = nf_def_dim(fid, trim(myVarLenName),&
                  size(var),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName)
             endif
             myvarDim(2) = dimid2
          endif
       endif

       ncStatus = nf_def_var(fid,trim(tmpVarName(iname)),nf_char, &
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varLongName)
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'units', len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varUnit)
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       do irec = 1, size(var)
          ncStatus = nf_put_vara_text(fid,varid(iname,iopen),&
               (/1,irec/),(/len_trim(var(irec)),1/), var(irec))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,tmpVarName(iname))
       enddo
    else
       ncStatus = nf_put_vara_text(fid, varid(iname,iopen),&
            (/1,record_no(iname,iopen)/),&
            (/len_trim(var(1)),1/), &
            trim(var(1)))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeString1D
  !--
  subroutine replaceString1D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    character (len=*), dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then   !writes all records
       do irec = 1, size(var)
          ncStatus = nf_put_vara_text(fid,vid,(/1,start_rec+irec-1/),&
               (/len_trim(var(irec)),1/),var(irec))
          if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
       enddo
    else                              !writes one record a time
       ncStatus = nf_put_vara_text(fid,vid,&
            (/1,start_rec/),&
            (/len_trim(var(1)),1/),trim(var(1)))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
    endif
  end subroutine replaceString1D
  !--
  subroutine appendString1D(fid, var, varName, append)
    integer, intent(in) :: fid
    character (len=*), dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceString1D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendString1D
  !--
  subroutine writeString2D(fid, var, varName, strLenName, &
       varLenName, varUnit, varLongName, append, unlimited_dim)
    character (len=*), dimension(:,:), intent (in) :: var
    character (len=*), intent (in), optional ::&
         varName, varLongName, varUnit
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    character (len=*), intent (in), optional :: strLenName
    integer, intent (in) :: fid
    logical, intent (in), optional :: append
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer :: dimid1, dimid2, dimid3, iopen
    integer, dimension(3) :: myVarDim

    call configureName(fid,iopen,&
         define,varName=varName)

    if (present(varLenName)) then
       myvarLenName1 = varLenName(1)
       myvarLenName2 = varLenName(2)
    else
       myvarLenName1 = trim(tmpVarName(iname))//'Len1'
       myvarLenName2 = trim(tmpVarName(iname))//'Len2'
    endif

    if (present(strLenName)) then
       myStrLenName = trim(strLenName)
    else
       myStrLenName = trim(tmpVarName(iname))//'StrLen'
    endif

    szB(1) = size(var,1)
    szB(2) = size(var,2)
    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
       ncStatus = nf_def_dim(fid,&
            myStrLenName,len_trim(var(1,1)),dimid1)
       if (ncStatus == -42) then
          ncStatus = nf_inq_dimid(fid,&
               trim(myStrLenName),dimid1)
       else
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,myStrLenName)
       endif
       myvarDim(1) = dimid1

       if (present(unlimited_dim)) then
          select case (unlimited_dim)
          case (1)
             myvarDim(2) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myVarLenName2),&
                  szB(2), dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myvarDim(3) = dimid3
          case (2)
             ncStatus = nf_def_dim(fid, trim(myVarLenName1),&
                  szB(1),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myvarDim(2) = dimid2
             myvarDim(3) = unlimited_dim_id(iopen)
          case default
             print *,'err[ncdf_module::writeString]: ',&
                  ' unlimited_dim for 2D string: 1 or 2'
             call exit(1)
          end select
       else
          ncStatus = nf_def_dim(fid, trim(myVarLenName1),&
               szB(1),dimid2)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName1),dimid2)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName1)
          endif
          myvarDim(2) = dimid2

          ncStatus = nf_def_dim(fid, trim(myVarLenName2),&
               szB(2),dimid3)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName2),dimid3)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName2)
          endif
          myvarDim(3) = dimid3
       endif
       ncStatus = nf_def_var(fid,trim(tmpVarName(iname)),nf_char, &
            size(myvarDim), myvarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varLongName)
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'units', len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varUnit)
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       do j = 1, szB(2)
          do i = 1, szB(1)
             ncStatus = nf_put_vara_text(fid,varid(iname,iopen),&
                  (/1,i,j/),(/len_trim(var(i,j)),1,1/),var(i,j))
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,tmpVarName(iname))
          enddo
       enddo
    else
       do i = 1, szB(1)
          ncStatus = nf_put_vara_text(fid,varid(iname,iopen),&
               (/1,i,record_no(iname,iopen)/),&
               (/len_trim(var(i,1)),1,1/),&
               var(i,1))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,tmpVarName(iname))
       enddo
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeString2D
  !--
  subroutine replaceString2D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    character (len=*), dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then
       do j = 1, size(var,2)
          do i = 1, size(var,1)
             ncStatus = nf_put_vara_text(fid,vid,(/1,i,start_rec+j-1/),&
                  (/len_trim(var(i,j)),1,1/),var(i,j))
             if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
          enddo
       enddo
    else
       do i = 1, size(var,1)
          ncStatus = nf_put_vara_text(fid,vid,&
               (/1,i,start_rec/),&
               (/len_trim(var(i,1)),1,1/),var(i,1))
          if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
       enddo
    endif
  end subroutine replaceString2D
  !--
  subroutine appendString2D(fid, var, varName, append)
    integer, intent(in) :: fid
    character (len=*), dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceString2D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendString2D
  !--
  subroutine writeString3D(fid, var, varName, strLenName,&
       varLenName, varUnit, varLongName, append, unlimited_dim)
    character (len=*), dimension(:,:,:), intent (in) :: var
    character (len=*), intent (in), optional ::&
         varName, varLongName, varUnit
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    character (len=*), intent (in), optional :: strLenName
    integer, intent (in) :: fid
    logical, intent (in), optional :: append
    logical :: define
    integer, intent (in), optional :: unlimited_dim
    integer :: dimid1, dimid2, dimid3, dimid4, iopen
    integer, dimension(4) :: myVarDim

    call configureName(fid,iopen,&
         define,varName=varName)

    if (present(varLenName)) then
       myvarLenName1 = varLenName(1)
       myvarLenName2 = varLenName(2)
       myvarLenName3 = varLenName(3)
    else
       myvarLenName1 = trim(tmpVarName(iname))//'Len1'
       myvarLenName2 = trim(tmpVarName(iname))//'Len2'
       myvarLenName3 = trim(tmpVarName(iname))//'Len3'
    endif

    if (present(strLenName)) then
       myStrLenName = trim(strLenName)
    else
       myStrLenName = trim(tmpVarName(iname))//'StrLen'
    endif

    szC(1) = size(var,1)
    szC(2) = size(var,2)
    szC(3) = size(var,3)
    if (define) then
       ncStatus = nf_redef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       ncStatus = nf_def_dim(fid,&
            myStrLenName,len(var(1,1,1)),dimid1)
       if (ncStatus == -42) then
          ncStatus = nf_inq_dimid(fid,&
               trim(myStrLenName),dimid1)
       else
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,myStrLenName)
       endif
       myVarDim(1) = dimid1

       if (present(unlimited_dim)) then
          select case (unlimited_dim)
          case (1)
             myVarDim(2) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myVarLenName2),&
                  szC(2), dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myVarDim(3) = dimid3

             ncStatus = nf_def_dim(fid, trim(myVarLenName3),&
                  szC(3), dimid4)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid4)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName3)
             endif
             myVarDim(4) = dimid4
          case (2)
             ncStatus = nf_def_dim(fid, trim(myVarLenName1),&
                  szC(1), dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myVarDim(2) = dimid2

             myVarDim(3) = unlimited_dim_id(iopen)

             ncStatus = nf_def_dim(fid, trim(myVarLenName3),&
                  szC(3), dimid4)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName3),dimid4)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName3)
             endif
             myVarDim(4) = dimid4
          case (3)
             ncStatus = nf_def_dim(fid, trim(myVarLenName1),&
                  szC(1),dimid2)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName1),dimid2)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName1)
             endif
             myVarDim(2) = dimid2

             ncStatus = nf_def_dim(fid, trim(myVarLenName2),&
                  szC(2), dimid3)
             if (ncStatus == -42) then
                ncStatus = nf_inq_dimid(fid,&
                     trim(myVarLenName2),dimid3)
             else
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,myVarLenName2)
             endif
             myVarDim(3) = dimid3

             myVarDim(4) = unlimited_dim_id(iopen)
          case default
             print *,'err[ncdf_module::writeString3D]: ',&
                  ' unlimited_dim for 2D string: 1 or 2'
             call exit(1)
          end select
       else
          ncStatus = nf_def_dim(fid, trim(myVarLenName1),&
               szC(1),dimid2)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName1),dimid2)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName1)
          endif
          myVarDim(2) = dimid2

          ncStatus = nf_def_dim(fid, trim(myVarLenName2),&
               szC(2), dimid3)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName2),dimid3)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName2)
          endif
          myVarDim(3) = dimid3

          ncStatus = nf_def_dim(fid, trim(myVarLenName3),&
               szC(3), dimid4)
          if (ncStatus == -42) then
             ncStatus = nf_inq_dimid(fid,&
                  trim(myVarLenName3),dimid4)
          else
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,myVarLenName3)
          endif
          myVarDim(4) = dimid4
       endif

       ncStatus = nf_def_var(fid,trim(tmpVarName(iname)),nf_char, &
            size(myVarDim), myVarDim, varid(iname,iopen))
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))

       if (present(varLongName)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'long_name',len_trim(varLongName),trim(varLongName))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varLongName)
       endif

       if (present(varUnit)) then
          ncStatus  = nf_put_att_text(fid,varid(iname,iopen),&
               'units', len_trim(varUnit),trim(varUnit))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varUnit)
       endif

       ncStatus = nf_enddef(fid)
       if (ncStatus .ne. nf_noerr) &
            call handle_err(ncStatus,tmpVarName(iname))
    endif

    if (.not. present(append)) then
       do k = 1, szC(3)
          do j = 1, szC(2)
             do i = 1, szC(1)
                ncStatus = nf_put_vara_text(fid,varid(iname,iopen),&
                     (/1,i,j,k/),(/len(var(i,j,k)),1,1,1/),&
                     var(i,j,k))
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,tmpVarName(iname))
             enddo
          enddo
       enddo
    else
       do j = 1, szC(2)
          do i = 1, szC(1)
             ncStatus = nf_put_vara_text(fid,&
                  varid(iname,iopen),(/1,i,j,&
                  record_no(iname,iopen)/),&
                  (/len(var(1,1,1)),1,1,1/),var(i,j,1))
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,tmpVarName(iname))
          enddo
       enddo
       record_no(iname,iopen) = record_no(iname,iopen) + 1
    endif

  end subroutine writeString3D
  !--
  subroutine replaceString3D(fid,var,varName,rec,append)
    integer, intent(in) :: fid
    character (len=*), dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    integer, intent(in), optional :: rec
    logical, intent (in), optional :: append
    integer :: vid, unlimid, unlimlen, start_rec

    if (present(rec)) then     !replaces
       start_rec = rec
    else                       !appends
       ncStatus = nf_inq_unlimdim(fid,unlimid)
       ncStatus = nf_inq_dimlen(fid,unlimid,unlimlen)
       start_rec = unlimlen+1
    endif
    ncStatus = nf_inq_varid(fid,trim(varname),vid)
    if (.not. present(append)) then
       do k = 1, size(var,3)
          do j = 1, size(var,2)
             do i = 1, size(var,1)
                ncStatus = nf_put_vara_text(fid,vid,(/1,i,j,start_rec+k-1/),&
                     (/len_trim(var(i,j,k)),1,1,1/),var(i,j,k))
                if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
             enddo
          enddo
       enddo
    else
       do j = 1, size(var,2)
          do i = 1, size(var,1)
             ncStatus = nf_put_vara_text(fid,vid,&
                  (/1,i,j,start_rec/),&
                  (/len_trim(var(1,1,1)),1,1,1/),var(i,j,1))
             if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varname)
          enddo
       enddo
    endif
  end subroutine replaceString3D
  !--
  subroutine appendString3D(fid, var, varName, append)
    integer, intent(in) :: fid
    character (len=*), dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    logical, intent (in), optional :: append

    call replaceString3D(fid=fid,var=var,varName=varName,append=append)

  end subroutine appendString3D
  !--
  subroutine readFloat1Attribute(fid, attr, attrName, silent, found)
    integer, intent (in) :: fid
    real*4, intent (out) :: attr
    character (len=*), intent (in) :: attrName
    logical, intent (in), optional :: silent
    logical, intent (out), optional :: found

    if (present(found)) then
       found = .true.
    endif

    ncStatus = nf_get_att_real(fid, nf_global, &
         trim(attrName), attr)
    if (ncStatus .ne. nf_noerr) then
       if (present(silent)) then
          if (silent) then
             attr = 0.0
             if (present(found)) then
                found = .false.
             endif
             return
          endif
       endif
       call handle_err(ncStatus,trim(attrName))
       attr = 0.0
    endif
    return
  end subroutine readFloat1Attribute
  !--
  subroutine readFloatAttribute(fid, attr, attrName, silent, found)
    integer, intent (in) :: fid
    real*4, dimension(:), intent (out) :: attr
    character (len=*), intent (in) :: attrName
    logical, intent (in), optional :: silent
    logical, intent (out), optional :: found

    if (present(found)) then
       found = .true.
    endif

    ncStatus = nf_get_att_real(fid, nf_global, &
         trim(attrName), attr)

    if (ncStatus .ne. nf_noerr) then
       if (present(silent)) then
          if (silent) then
             attr = 0.0
             if (present(found)) then
                found = .false.
             endif
             return
          endif
       endif
       call handle_err(ncStatus,trim(attrName))
       attr = 0.0
    endif
    return
  end subroutine readFloatAttribute
  !--
  subroutine readDouble1Attribute(fid, attr, attrName, silent, found)
    integer, intent (in) :: fid
    real*8, intent (out) :: attr
    character (len=*), intent (in) :: attrName
    logical, intent (in), optional :: silent
    logical, intent (out), optional :: found

    real*4              :: attrLoc

    if (present(found)) then
       found = .true.
    endif

    ncStatus = nf_get_att_real(fid, nf_global, &
         trim(attrName), attrLoc)

    if (ncStatus .ne. nf_noerr) then
       if (present(silent)) then
          if (silent) then
             attr = 0.0
             if (present(found)) then
                found = .false.
             endif
             return
          endif
       endif
       call handle_err(ncStatus,trim(attrName))
       attr = 0.0
    else
       attr = attrLoc
    endif
    return
  end subroutine readDouble1Attribute
  !--
  subroutine readDoubleAttribute(fid, attr, attrName, silent, found)
    integer, intent (in) :: fid
    real*8, dimension(:), intent (out) :: attr
    character (len=*), intent (in) :: attrName
    logical, intent (in), optional :: silent
    logical, intent (out), optional :: found

    real*4, dimension(size(attr)) :: attrLoc

    if (present(found)) then
       found = .true.
    endif

    ncStatus = nf_get_att_real(fid, nf_global, &
         trim(attrName), attrLoc)

    if (ncStatus .ne. nf_noerr) then
       if (present(silent)) then
          if (silent) then
             attr = 0.0
             if (present(found)) then
                found = .false.
             endif
             return
          endif
       endif
       call handle_err(ncStatus,trim(attrName))
       attr = 0.0
    else
       attr = attrLoc
    endif
    return
  end subroutine readDoubleAttribute
 !--
  subroutine readInt1Attribute(fid, attr, attrName, silent, found)
    integer, intent (in) :: fid
    integer, intent (out) :: attr
    character (len=*), intent (in) :: attrName
    logical, intent (in), optional :: silent
    logical, intent (out), optional :: found

    if (present(found)) then
       found = .true.
    endif

    ncStatus = nf_get_att_int(fid, nf_global, &
         trim(attrName), attr)

    if (ncStatus .ne. nf_noerr) then
       if (present(silent)) then
          if (silent) then
             attr = 0
             if (present(found)) then
                found = .false.
             endif
             return
          endif
       endif
       call handle_err(ncStatus,trim(attrName))
       attr = 0
    endif
    return
  end subroutine readInt1Attribute
  !--
  subroutine readIntAttribute(fid, attr, attrName, silent, found)
    integer, intent (in) :: fid
    integer, dimension(:), intent (out) :: attr
    character (len=*), intent (in) :: attrName
    logical, intent (in), optional :: silent
    logical, intent (out), optional :: found

    if (present(found)) then
       found = .true.
    endif

    ncStatus = nf_get_att_int(fid, nf_global, &
         trim(attrName), attr)

    if (ncStatus .ne. nf_noerr) then
       if (present(silent)) then
          if (silent) then
             attr = 0
             if (present(found)) then
                found = .false.
             endif
             return
          endif
       endif
       call handle_err(ncStatus,trim(attrName))
       attr = 0
    endif
    return
  end subroutine readIntAttribute
  !--
  subroutine readText1Attribute(fid, attr, attrName, silent, found)
    integer, intent (in) :: fid
    character (len=*), intent (out) :: attr
    character (len=*), intent (in) :: attrName
    logical, intent (in), optional :: silent
    logical, intent (out), optional :: found
    integer :: attrLen, copyLen
    character (len=MAX_ATTR_LENGTH) :: attrLocal

    if (present(found)) then
       found = .true.
    endif

    ncStatus = nf_inq_attlen(fid, nf_global, &
         trim(attrName), attrLen)

    if (ncStatus .ne. nf_noerr) then
       if (present(silent)) then
          if (silent) then
             attr = '\0'
             if (present(found)) then
                found = .false.
             endif
             return
          endif
       endif
       call handle_err(ncStatus,trim(attrName))
       attr = '\0'
       return
    endif

    if (attrLen > MAX_ATTR_LENGTH) then
       print *,'err[ncdf_module::readText1Attribute]: ', &
          'Attribute length exceeds max;'
       print *,'attrName,attrLen,MAX_ATTR_LENGTH: ', &
          trim(attrName),attrLen,MAX_ATTR_LENGTH
       call exit(1)
    endif

    ncStatus = nf_get_att_text(fid, nf_global, &
         trim(attrName), attrLocal)

    if (ncStatus .ne. nf_noerr) then
       call handle_err(ncStatus,trim(attrName))
       attr = '\0'
    endif

    attr=REPEAT(' ',LEN(attr))   ! Pad with blanks
    copyLen=min(attrLen,len(attr))
    attr(1:copyLen)=attrLocal(1:copyLen)

    return
  end subroutine readText1Attribute
  !--
  integer function readUnlimitDim(fid) result(dimLen)
    integer, intent (in) :: fid
    integer :: unlimid

    ncStatus = nf_inq_unlimdim(fid,unlimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,'unlimited dim')

    ncStatus = nf_inq_dimlen(fid,unlimid,dimLen)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,'unlimited dim')

  end function readUnlimitDim
  !--
  integer function readGlobalDim(fid,dimName,silent) result(dimLen)
    integer, intent (in) :: fid
    character (len=*), intent (in) :: dimName
    logical, intent (in), optional :: silent
    integer :: tmpid

    ncStatus = nf_inq_dimid(fid,trim(dimName),tmpid)

    if (ncStatus .ne. nf_noerr) then
       dimLen=0
       if (present(silent)) then
          if (silent) then
             return
          endif
       endif
       call handle_err(ncStatus,trim(dimName))
       return
    endif

    ncStatus = nf_inq_dimlen(fid,tmpid,dimlen)
    if (ncStatus .ne. nf_noerr) then
       dimLen=0
       call handle_err(ncStatus,dimName)
    endif

  end function readGlobalDim
  !--
  integer function readVarDim(fid, varName, silent) result(nDim)
    integer, intent (in) :: fid
    character (len=*), intent (in) :: varName
    logical, intent (in), optional :: silent
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),tmpvarid)
    if (ncStatus .ne. nf_noerr) then
       nDim=0
       if (present(silent)) then
          if (silent) then
             return
          endif
       endif
       call handle_err(ncStatus,'<readVarDim>'//trim(varName))
       return
    endif

    ncStatus = nf_inq_varndims(fid,tmpvarid,nDim)
    if (ncStatus .ne. nf_noerr) then
       nDim=0
       call handle_err(ncStatus,'<readVarDim>'//varName)
    endif

  end function readVarDim
  !--
  subroutine readInteger1D(fid, var, varName, varLen, &
       varUnit, varLongName, record_no)
    integer, dimension(:), intent (out) :: var
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional :: &
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, intent (out), optional :: varLen
    integer, intent (in), optional :: record_no
    integer :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,varName)
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid,varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen = varLen1

    if (present(record_no)) then
       ncStatus = nf_get_var1_int(fid,tmpvarid,&
            record_no,var(1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    else
       if (varLen1 /= size(var)) then
          print *,'err[ncdf_module::readInteger1D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_int(fid,tmpvarid,1,&
            varLen1,var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    endif

  end subroutine readInteger1D
  !--
  subroutine readInteger2D(fid, var, varName, varLen,&
       varUnit, varLongName, record_no)
    integer, dimension(:,:), intent (out) :: var
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional :: &
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(2), intent (out), optional :: varLen
    integer, dimension(2) :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),&
         tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(varName))
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(1),varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(2),varLen2)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen = (/varLen1, varLen2/)

    if (present(record_no)) then
       if (varLen1 /= size(var,1)) then
          print *,'err[ncdf_module::readInteger2D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_int(fid,tmpvarid,&
            (/1,record_no/),(/varlen1,1/),var(:,1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    else
       if (varLen1 /= size(var,1) .or. &
            varLen2 /= size(var,2)) then
          print *,'err[ncdf_module::readInteger2D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_int(fid,tmpvarid,&
            (/1,1/),(/varlen1,varlen2/),var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    endif

  end subroutine readInteger2D
  !--
  subroutine readInteger3D(fid, var, varName, varLen,&
       varUnit, varLongName, record_no)
    integer, dimension(:,:,:), intent (out) :: var
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional :: &
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, intent (out), dimension(3), optional :: varLen
    integer, dimension(3) :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),&
         tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(varName))
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(1),varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(2),varLen2)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(3),varLen3)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen = (/varLen1, varLen2, varLen3/)

    if (present(record_no)) then
       if (varLen1 /= size(var,1) .or. &
            varLen2 /= size(var,2)) then
          print *,'err[ncdf_module::readInteger3D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_int(fid,tmpvarid,&
            (/1,1,record_no/),(/varlen1,varLen2,1/),&
            var(:,:,1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    else
       if (varLen1 /= size(var,1) .or. &
            varLen2 /= size(var,2) .or.&
            varLen3 /= size(var,3)) then
          print *,'err[ncdf_module::readInteger3D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_int(fid,tmpvarid,&
            (/1,1,1/),(/varlen1,varlen2,varlen3/),var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    endif

  end subroutine readInteger3D
  !--
  subroutine readInteger4D(fid, var, varName, varLen,&
       varUnit, varLongName, record_no)
    integer, dimension(:,:,:,:), intent (out) :: var
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional :: &
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, intent (out), dimension(4), optional :: varLen
    integer, dimension(4) :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),&
         tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(varName))
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(1),varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(2),varLen2)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(3),varLen3)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(4),varLen4)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen =&
         (/varLen1, varLen2, varLen3, varLen4/)

    if (present(record_no)) then
       if (varLen1 /= size(var,1) .or. &
            varLen2 /= size(var,2) .or. &
            varLen3 /= size(var,3)) then
          print *,'err[ncdf_module::readInteger4D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_int(fid,tmpvarid,&
            (/1,1,1,record_no/),(/varlen1,varLen2,varLen3,1/),&
            var(:,:,:,1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    else
       if (varLen1 /= size(var,1) .or. &
            varLen2 /= size(var,2) .or.&
            varLen3 /= size(var,3) .or.&
            varLen4 /= size(var,4)) then
          print *,'err[ncdf_module::readInteger4D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_int(fid,tmpvarid,&
            (/1,1,1,1/),(/varlen1,varlen2,varlen3,varlen4/),var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    endif

  end subroutine readInteger4D
  !--
  subroutine readInteger1Byte1D(fid, var, varName, varLen, &
       varUnit, varLongName, record_no)
    integer (kind=1), dimension(:), intent (out) :: var
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional :: &
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, intent (out), optional :: varLen
    integer, intent (in), optional :: record_no
    integer :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,varName)
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid,varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen = varLen1

    if (present(record_no)) then
       ncStatus = nf_get_var1_int1(fid,tmpvarid,&
            record_no,var(1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    else
       if (varLen1 /= size(var)) then
          print *,'err[ncdf_module::readInteger1Byte1D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_int1(fid,tmpvarid,1,&
            varLen1,var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    endif

  end subroutine readInteger1Byte1D
  !--
  subroutine readInteger1Byte2D(fid, var, varName, varLen,&
       varUnit, varLongName, record_no)
    integer (kind=1), dimension(:,:), intent (out) :: var
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional :: &
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(2), intent (out), optional :: varLen
    integer, dimension(2) :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),&
         tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(varName))
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(1),varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(2),varLen2)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen = (/varLen1, varLen2/)

    if (present(record_no)) then
       if (varLen1 /= size(var,1)) then
          print *,'err[ncdf_module::readInteger1Byte2D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_int1(fid,tmpvarid,&
            (/1,record_no/),(/varlen1,1/),var(:,1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    else
       if (varLen1 /= size(var,1) .or. &
            varLen2 /= size(var,2)) then
          print *,'err[ncdf_module::readInteger1Byte2D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_int1(fid,tmpvarid,&
            (/1,1/),(/varlen1,varlen2/),var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    endif

  end subroutine readInteger1Byte2D
  !--
  subroutine readInteger1Byte3D(fid, var, varName, varLen,&
       varUnit, varLongName, record_no)
    integer (kind=1), dimension(:,:,:), intent (out) :: var
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional :: &
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, intent (out), dimension(3), optional :: varLen
    integer, dimension(3) :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),&
         tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(varName))
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(1),varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(2),varLen2)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(3),varLen3)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen = (/varLen1, varLen2, varLen3/)

    if (present(record_no)) then
       if (varLen1 /= size(var,1) .or. &
            varLen2 /= size(var,2)) then
          print *,'err[ncdf_module::readInteger1Byte3D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_int1(fid,tmpvarid,&
            (/1,1,record_no/),(/varlen1,varLen2,1/),&
            var(:,:,1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    else
       if (varLen1 /= size(var,1) .or. &
            varLen2 /= size(var,2) .or.&
            varLen3 /= size(var,3)) then
          print *,'err[ncdf_module::readInteger1Byte3D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_int1(fid,tmpvarid,&
            (/1,1,1/),(/varlen1,varlen2,varlen3/),var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    endif

  end subroutine readInteger1Byte3D
  !--
  subroutine readInteger1Byte4D(fid, var, varName, varLen,&
       varUnit, varLongName, record_no)
    integer*1, dimension(:,:,:,:), intent (out) :: var
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional :: &
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, intent (out), dimension(4), optional :: varLen
    integer, dimension(4) :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),&
         tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(varName))
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(1),varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(2),varLen2)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(3),varLen3)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(4),varLen4)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen =&
         (/varLen1, varLen2, varLen3, varLen4/)

    if (present(record_no)) then
       if (varLen1 /= size(var,1) .or. &
            varLen2 /= size(var,2) .or. &
            varLen3 /= size(var,3)) then
          print *,'err[ncdf_module::readInteger1Byte4D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_int1(fid,tmpvarid,&
            (/1,1,1,record_no/),(/varlen1,varLen2,varLen3,1/),&
            var(:,:,:,1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    else
       if (varLen1 /= size(var,1) .or. &
            varLen2 /= size(var,2) .or.&
            varLen3 /= size(var,3) .or.&
            varLen4 /= size(var,4)) then
          print *,'err[ncdf_module::readInteger1Byte4D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_int1(fid,tmpvarid,&
            (/1,1,1,1/),(/varlen1,varlen2,varlen3,varlen4/),var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    endif

  end subroutine readInteger1Byte4D
  !--
  subroutine readInteger2Byte1D(fid, var, varName, varLen, &
       varUnit, varLongName, record_no)
    integer (kind=2), dimension(:), intent (out) :: var
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional :: &
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, intent (out), optional :: varLen
    integer, intent (in), optional :: record_no
    integer :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,varName)
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid,varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen = varLen1

    if (present(record_no)) then
       ncStatus = nf_get_var1_int2(fid,tmpvarid,&
            record_no,var(1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    else
       if (varLen1 /= size(var)) then
          print *,'err[ncdf_module::readInteger2Byte1D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_int2(fid,tmpvarid,1,&
            varLen1,var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    endif

  end subroutine readInteger2Byte1D
  !--
  subroutine readInteger2Byte2D(fid, var, varName, varLen,&
       varUnit, varLongName, record_no)
    integer (kind=2), dimension(:,:), intent (out) :: var
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional :: &
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(2), intent (out), optional :: varLen
    integer, dimension(2) :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),&
         tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(varName))
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(1),varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(2),varLen2)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen = (/varLen1, varLen2/)

    if (present(record_no)) then
       if (varLen1 /= size(var,1)) then
          print *,'err[ncdf_module::readInteger2Byte2D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_int2(fid,tmpvarid,&
            (/1,record_no/),(/varlen1,1/),var(:,1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    else
       if (varLen1 /= size(var,1) .or. &
            varLen2 /= size(var,2)) then
          print *,'err[ncdf_module::readInteger2Byte2D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_int2(fid,tmpvarid,&
            (/1,1/),(/varlen1,varlen2/),var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    endif

  end subroutine readInteger2Byte2D
  !--
  subroutine readInteger2Byte3D(fid, var, varName, varLen,&
       varUnit, varLongName, record_no)
    integer (kind=2), dimension(:,:,:), intent (out) :: var
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional :: &
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, intent (out), dimension(3), optional :: varLen
    integer, dimension(3) :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),&
         tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(varName))
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(1),varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(2),varLen2)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(3),varLen3)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen = (/varLen1, varLen2, varLen3/)

    if (present(record_no)) then
       if (varLen1 /= size(var,1) .or. &
            varLen2 /= size(var,2)) then
          print *,'err[ncdf_module::readInteger2Byte3D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_int2(fid,tmpvarid,&
            (/1,1,record_no/),(/varlen1,varLen2,1/),&
            var(:,:,1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    else
       if (varLen1 /= size(var,1) .or. &
            varLen2 /= size(var,2) .or.&
            varLen3 /= size(var,3)) then
          print *,'err[ncdf_module::readInteger2Byte3D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_int2(fid,tmpvarid,&
            (/1,1,1/),(/varlen1,varlen2,varlen3/),var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    endif

  end subroutine readInteger2Byte3D
  !--
  subroutine readInteger2Byte4D(fid, var, varName, varLen,&
       varUnit, varLongName, record_no)
    integer*2, dimension(:,:,:,:), intent (out) :: var
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional :: &
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, intent (out), dimension(4), optional :: varLen
    integer, dimension(4) :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),&
         tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(varName))
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(1),varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(2),varLen2)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(3),varLen3)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(4),varLen4)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen =&
         (/varLen1, varLen2, varLen3, varLen4/)

    if (present(record_no)) then
       if (varLen1 /= size(var,1) .or. &
            varLen2 /= size(var,2) .or. &
            varLen3 /= size(var,3)) then
          print *,'err[ncdf_module::readInteger2Byte4D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_int2(fid,tmpvarid,&
            (/1,1,1,record_no/),(/varlen1,varLen2,varLen3,1/),&
            var(:,:,:,1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    else
       if (varLen1 /= size(var,1) .or. &
            varLen2 /= size(var,2) .or.&
            varLen3 /= size(var,3) .or.&
            varLen4 /= size(var,4)) then
          print *,'err[ncdf_module::readInteger2Byte4D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_int2(fid,tmpvarid,&
            (/1,1,1,1/),(/varlen1,varlen2,varlen3,varlen4/),var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    endif

  end subroutine readInteger2Byte4D
  !--
  subroutine readFloat1D(fid, var, varName, varLen, &
       varUnit, varLongName, record_no)
    real*4, dimension(:), intent (out) :: var
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional :: &
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, intent (out), optional :: varLen
    integer, intent (in), optional :: record_no
    integer :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),&
         tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(varName))
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid,varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen = varLen1

    if (present(record_no)) then
       ncStatus = nf_get_var1_real(fid,tmpvarid,&
            record_no,var(1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    else
       if (varLen1 /= size(var)) then
          print *,'err[ncdf_module::readFloat1D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_real(fid,tmpvarid,1,&
            varLen1,var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    endif

  end subroutine readFloat1D
  !--
  subroutine readFloat2D(fid, var, varName, varLen,&
       varUnit, varLongName, record_no)
    real*4, dimension(:,:), intent (out) :: var
    integer, dimension(2), intent (out), optional :: varLen
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional :: &
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(2) :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),&
         tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(varName))
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(1),varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(2),varLen2)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen = (/varLen1, varLen2/)

    if (present(record_no)) then
       if (varLen1 /= size(var,1)) then
          print *,'err[ncdf_module::readFloat2D]: ',&
               ' ('//trim(varName)//') variable size mismatched', varLen1,  size(var,1)
          call exit(1)
       end if
       ncStatus = nf_get_vara_real(fid,tmpvarid,&
            (/1,record_no/),(/varlen1,1/),var(:,1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    else
       if (varLen1 /= size(var,1) .or.&
            varLen2 /= size(var,2)) then
          print *,'err[ncdf_module::readFloat2D]: ',&
               ' ('//trim(varName)//') variable size mismatched', varLen2,  size(var,2)
          call exit(1)
       end if
       ncStatus = nf_get_vara_real(fid,tmpvarid,&
            (/1,1/),(/varlen1,varlen2/),var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    endif

  end subroutine readFloat2D
  !--
  subroutine readFloat3D(fid, var, varName, varLen,&
       varUnit, varLongName, record_no)
    real*4, dimension(:,:,:), intent (out) :: var
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional :: &
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(3), intent (out), optional :: varLen
    integer, dimension(3) :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),&
         tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(varName))
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(1),varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(2),varLen2)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(3),varLen3)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen = (/varLen1, varLen2, varLen3/)

    if (present(record_no)) then
       if (varLen1 /= size(var,1) .or.&
            varLen2 /= size(var,2)) then
          print *,'err[ncdf_module::readFloat3D]: ',&
               ' ('//trim(varName)//') variable size mismatched',&
               varLen1,  size(var,1), varLen2,  size(var,2)
          call exit(1)
       end if
       ncStatus = nf_get_vara_real(fid,tmpvarid,&
            (/1,1,record_no/),(/varlen1,varlen2,1/),var(:,:,1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    else
       if (varLen1 /= size(var,1) .or.&
            varLen2 /= size(var,2) .or.&
            varLen3 /= size(var,3)) then
          print *,'err[ncdf_module::readFloat3D]: ',&
               ' ('//trim(varName)//') variable size mismatched', &
               varLen1,  size(var,1), varLen2,  size(var,2), varLen3,  size(var,3)
          call exit(1)
       end if
       ncStatus = nf_get_vara_real(fid,tmpvarid,&
            (/1,1,1/),(/varlen1,varlen2,varlen3/),var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    endif

  end subroutine readFloat3D
  !--
  subroutine readFloat4D(fid, var, varName, varLen,&
       varUnit, varLongName, record_no)
    real*4, dimension(:,:,:,:), intent (out) :: var
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional :: &
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, intent (out), dimension(4), optional :: varLen
    integer, dimension(4) :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),&
         tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(varName))
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(1),varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(2),varLen2)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(3),varLen3)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(4),varLen4)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen =&
         (/varLen1, varLen2, varLen3, varLen4/)

    if (present(record_no)) then
       if (varLen1 /= size(var,1) .or. &
            varLen2 /= size(var,2) .or. &
            varLen3 /= size(var,3)) then
          print *,'err[ncdf_module::readFloat4D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_real(fid,tmpvarid,&
            (/1,1,1,record_no/),(/varlen1,varLen2,varLen3,1/),&
            var(:,:,:,1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    else
       if (varLen1 /= size(var,1) .or. &
            varLen2 /= size(var,2) .or.&
            varLen3 /= size(var,3) .or.&
            varLen4 /= size(var,4)) then
          print *,'err[ncdf_module::readFloat4D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_real(fid,tmpvarid,&
            (/1,1,1,1/),(/varlen1,varlen2,varlen3,varlen4/),var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    endif

  end subroutine readFloat4D
  !--
  subroutine readDFloat1D(fid, var, varName, varLen, &
       varUnit, varLongName, record_no)
    double precision, dimension(:), intent (out) :: var
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional :: &
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, intent (out), optional :: varLen
    integer :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),&
         tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(varName))
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid,varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen = varLen1

    if (present(record_no)) then
       ncStatus = nf_get_var1_double(fid,&
            tmpvarid,record_no,var(1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    else
       if (varLen1 /= size(var)) then
          print *,'err[ncdf_module::readDFloat1D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_double(fid,&
            tmpvarid,1,varLen1,var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    endif

  end subroutine readDFloat1D
  !--
  subroutine readDFloat2D(fid, var, varName, varLen,&
       varUnit, varLongName, record_no)
    double precision, dimension(:,:), intent (out) :: var
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional ::&
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(2), intent (out), optional :: varLen
    integer, dimension(2) :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),&
         tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(varName))
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(1),varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(2),varLen2)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen = (/varLen1, varLen2/)

    if (present(record_no)) then
       if (varLen1 /= size(var,1)) then
          print *,'err[ncdf_module::readDFloat2D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_double(fid,tmpvarid,&
            (/1,record_no/),(/varlen1,1/),var(:,1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    else
       if (varLen1 /= size(var,1) .or.&
            varLen2 /= size(var,2)) then
          print *,'err[ncdf_module::readDFloat2D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_double(fid,tmpvarid,&
            (/1,1/),(/varlen1,varlen2/),var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    endif

  end subroutine readDFloat2D
  !--
  subroutine readDFloat3D(fid, var, varName, varLen,&
       varUnit, varLongName, record_no)
    double precision, dimension(:,:,:), intent (out) :: var
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional :: &
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, dimension(3), intent (out), optional :: varLen
    integer, intent (in), optional :: record_no
    integer, dimension(3) :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),&
         tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(varName))
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(1),varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(2),varLen2)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(3),varLen3)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen = (/varLen1, varLen2, varLen3/)

    if (present(record_no)) then
       if (varLen1 /= size(var,1) .or.&
            varLen2 /= size(var,2)) then
          print *,'err[ncdf_module::readDFloat3D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_double(fid,tmpvarid,&
            (/1,1,record_no/),(/varlen1,varlen2,1/),var(:,:,1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    else
       if (varLen1 /= size(var,1) .or.&
            varLen2 /= size(var,2) .or.&
            varLen3 /= size(var,3)) then
          print *,'err[ncdf_module::readDFloat3D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_double(fid,tmpvarid,&
            (/1,1,1/),(/varlen1,varlen2,varlen3/),var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    endif

  end subroutine readDFloat3D
  !--
  subroutine readString1D(fid, var, varName, varLen,&
       varUnit, varLongName, record_no)
    character (len=*), dimension(:), intent (out) :: var
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional :: &
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, intent (out), optional :: varLen
    integer, dimension(2) :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),&
         tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(varName))
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(1),strLen)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (len(var(1)) < strLen) then
       print *,'err[ncdf_module::readString1D]:  string length'
       call exit(1)
    end if

    ncStatus = nf_inq_dimlen(fid,mydimid(2),varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen = varLen1

    if (present(record_no)) then
       ncStatus = nf_get_vara_text(fid,tmpvarid,&
            (/1,record_no/),(/strlen,1/),var(1))
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    else
       if (varLen1 /= size(var)) then
          print *,'err[ncdf_module::readString1D: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       ncStatus = nf_get_vara_text(fid,tmpvarid,&
            (/1,1/),(/strlen,varlen1/),var)
       if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)
    endif

  end subroutine readString1D
  !--
  subroutine readString2D(fid, var, varName, varLen,&
       varUnit, varLongName, record_no)
    character (len=*), dimension(:,:), intent (out) :: var
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional :: &
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(2), intent (out), optional :: varLen
    integer, dimension(3) :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),&
         tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(varName))
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(1),strLen)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (len(var(1,1)) < strLen) then
       print *,'err[ncdf_module::readString2D]:  string length'
       call exit(1)
    end if

    ncStatus = nf_inq_dimlen(fid,mydimid(2),varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(3),varLen2)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen = (/varLen1, varLen2/)

    if (present(record_no)) then
       if (varLen1 /= size(var,1)) then
          print *,'err[readString2D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       do i = 1, varlen1
          ncStatus = nf_get_vara_text(fid,tmpvarid,&
               (/1,i,record_no/),(/strlen,1,1/),var(i,1))
          if (ncStatus .ne. nf_noerr) &
               call handle_err(ncStatus,varName)
       enddo
    else
       if (varLen1 /= size(var,1) .or.&
            varLen2 /= size(var,2)) then
          print *,'err[ncdf_module::readString2D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       do j = 1, varlen2
          do i = 1, varlen1
             ncStatus = nf_get_vara_text(fid,tmpvarid,&
                  (/1,i,j/),(/strlen,1,1/),var(i,j))
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,varName)
          enddo
       enddo
    endif

  end subroutine readString2D
  !--
  subroutine readString3D(fid, var, varName, varLen,&
       varUnit, varLongName, record_no)
    character (len=*), dimension(:,:,:), intent (out) :: var
    character (len=*), intent (in) :: varName
    character (len=*), intent (out), optional :: &
         varUnit, varLongName
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(3), intent (out), optional :: varLen
    integer, dimension(4) :: mydimid
    integer :: tmpvarid

    ncStatus = nf_inq_varid(fid,trim(varName),&
         tmpvarid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,trim(varName))
    if (ncStatus < 0) return

    ncStatus = nf_inq_vardimid(fid,tmpvarid,mydimid)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(1),strLen)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (len(var(1,1,1)) < strLen) then
       print *,'err[ncdf_module::readString3D]:  string length'
       call exit(1)
    end if

    ncStatus = nf_inq_dimlen(fid,mydimid(2),varLen1)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(3),varLen2)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    ncStatus = nf_inq_dimlen(fid,mydimid(4),varLen3)
    if (ncStatus .ne. nf_noerr) call handle_err(ncStatus,varName)

    if (present(varLen)) varLen = (/varLen1, varLen2, varLen3/)

    if (present(record_no)) then
       if (varLen1 /= size(var,1) .or.&
            varLen2 /= size(var,2)) then
          print *,'err[ncdf_module::readString3D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       do j = 1, varlen2
          do i = 1, varlen1
             ncStatus = nf_get_vara_text(fid,tmpvarid,&
                  (/1,i,j,record_no/),(/strlen,1,1,1/),var(i,j,1))
             if (ncStatus .ne. nf_noerr) &
                  call handle_err(ncStatus,varName)
          enddo
       enddo
    else
       if (varLen1 /= size(var,1) .or.&
            varLen2 /= size(var,2) .or.&
            varLen3 /= size(var,3)) then
          print *,'err[ncdf_module::readString3D]: ',&
               ' ('//trim(varName)//') variable size mismatched'
          call exit(1)
       end if
       do k = 1, varlen3
          do j = 1, varlen2
             do i = 1, varlen1
                ncStatus = nf_get_vara_text(fid,tmpvarid,&
                     (/1,i,j,k/),(/strlen,1,1,1/),var(i,j,k))
                if (ncStatus .ne. nf_noerr) &
                     call handle_err(ncStatus,varName)
             enddo
          enddo
       enddo
    endif

  end subroutine readString3D
  !--
  subroutine writeInt1DParser(fid, var, varName, varLenName, &
       varUnit, varLongName, status, unlimited_dim)
    integer, intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varLenName, varUnit
    character (len=*), intent (in) :: status
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    integer, dimension (:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=1
    else
       unlimited=unlimited_dim
    endif

    if (.not. allocated(buf)) allocate(buf(1))
    buf(1) = var
    call writeInteger1D(fid, buf, &
         varName=varName,&
         varLenName=varLenname, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine writeInt1DParser
  !--
  subroutine replaceInt1DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    integer, intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    integer, dimension (:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(1))
    buf(1) = var
    call replaceInteger1D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)
  end subroutine replaceInt1DParser
  !--
  subroutine appendInt1DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    integer, intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    integer, dimension (:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(1))
    buf(1) = var

    call appendInteger1D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendInt1DParser
  !--
  subroutine writeInt2DParser(fid, var, varName, varLenName, &
       varUnit, varLongName, status, unlimited_dim)
    integer, dimension(:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varUnit
    character (len=*), intent (in) :: status
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    integer, dimension (:,:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=2
    else
       unlimited=unlimited_dim
    endif

    if (.not. allocated(buf)) allocate(buf(size(var),1))
    buf(:,1) = var
    call writeInteger2D(fid, buf, &
         varName=varName,&
         varLenName=varLenname, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine writeInt2DParser
  !--
  subroutine replaceInt2DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    integer, dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    integer, dimension (:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var),1))
    buf(:,1) = var
    call replaceInteger2D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)
  end subroutine replaceInt2DParser
  !--
  subroutine appendInt2DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    integer, dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    integer, dimension (:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var),1))
    buf(:,1) = var

    call appendInteger2D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendInt2DParser
  !--
  subroutine writeInt3DParser(fid, var, varName, varLenName, &
       varUnit, varLongName, status, unlimited_dim)
    integer, dimension(:,:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varUnit
    character (len=*), intent (in) :: status
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    integer, dimension (:,:,:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=3
    else
       unlimited=unlimited_dim
    endif

    if (.not. allocated(buf)) &
         allocate(buf(size(var,1),size(var,2),1))
    buf(:,:,1) = var
    call writeInteger3D(fid, buf, &
         varName=varName,&
         varLenName=varLenname, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine writeInt3DParser
  !--
  subroutine replaceInt3DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    integer, dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    integer, dimension (:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var,1),size(var,2),1))
    buf(:,:,1) = var
    call replaceInteger3D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine replaceInt3DParser
  !--
  subroutine appendInt3DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    integer, dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    integer, dimension (:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var,1),size(var,2),1))
    buf(:,:,1) = var

    call appendInteger3D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendInt3DParser
  !--
  subroutine writeInt4DParser(fid, var, varName, varLenName, &
       varUnit, varLongName, status, unlimited_dim)
    integer, dimension(:,:,:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varUnit
    character (len=*), intent (in) :: status
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    integer, dimension (:,:,:,:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=4
    else
      unlimited=unlimited_dim
    endif

    if (.not. allocated(buf)) &
         allocate(buf(size(var,1),size(var,2),size(var,3),1))
    buf(:,:,:,1) = var
    call writeInteger4D(fid, buf, &
         varName=varName,&
         varLenName=varLenname, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine writeInt4DParser
  !--
  subroutine replaceInt4DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    integer, dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    integer, dimension (:,:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) &
         allocate(buf(size(var,1),size(var,2),size(var,3),1))
    buf(:,:,:,1) = var
    call replaceInteger4D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine replaceInt4DParser
  !--
  subroutine appendInt4DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    integer, dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    integer, dimension (:,:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) &
         allocate(buf(size(var,1),size(var,2),size(var,3),1))
    buf(:,:,:,1) = var

    call appendInteger4D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendInt4DParser
  !--
  subroutine writeInt1Byte1DParser(fid, var, varName, varLenName, &
       varUnit, varLongName, status, unlimited_dim)
    integer (kind=1), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varLenName, varUnit
    character (len=*), intent (in) :: status
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    integer (kind=1), dimension (:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=1
    else
       unlimited=unlimited_dim
    endif

    if (.not. allocated(buf)) allocate(buf(1))
    buf(1) = var
    call writeInteger1Byte1D(fid, buf, &
         varName=varName,&
         varLenName=varLenname, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine writeInt1Byte1DParser
  !--
  subroutine replaceInt1Byte1DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    integer (kind=1), intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    integer (kind=1), dimension (:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(1))
    buf(1) = var
    call replaceInteger1Byte1D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)
  end subroutine replaceInt1Byte1DParser
  !--
  subroutine appendInt1Byte1DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    integer (kind=1), intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    integer (kind=1), dimension (:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(1))
    buf(1) = var

    call appendInteger1Byte1D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendInt1Byte1DParser
  !--
  subroutine writeInt1Byte2DParser(fid, var, varName, varLenName, &
       varUnit, varLongName, status, unlimited_dim)
    integer (kind=1), dimension(:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varUnit
    character (len=*), intent (in) :: status
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    integer (kind=1), dimension (:,:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=2
    else
       unlimited=unlimited_dim
    endif

    if (.not. allocated(buf)) allocate(buf(size(var),1))
    buf(:,1) = var
    call writeInteger1Byte2D(fid, buf, &
         varName=varName,&
         varLenName=varLenname, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine writeInt1Byte2DParser
  !--
  subroutine replaceInt1Byte2DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    integer (kind=1), dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    integer (kind=1), dimension (:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var),1))
    buf(:,1) = var
    call replaceInteger1Byte2D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)
  end subroutine replaceInt1Byte2DParser
  !--
  subroutine appendInt1Byte2DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    integer (kind=1), dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    integer (kind=1), dimension (:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var),1))
    buf(:,1) = var

    call appendInteger1Byte2D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendInt1Byte2DParser
  !--
  subroutine writeInt1Byte3DParser(fid, var, varName, varLenName, &
       varUnit, varLongName, status, unlimited_dim)
    integer (kind=1), dimension(:,:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varUnit
    character (len=*), intent (in) :: status
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    integer (kind=1), dimension (:,:,:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=3
    else
       unlimited=unlimited_dim
    endif

    if (.not. allocated(buf)) &
         allocate(buf(size(var,1),size(var,2),1))
    buf(:,:,1) = var
    call writeInteger1Byte3D(fid, buf, &
         varName=varName,&
         varLenName=varLenname, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine writeInt1Byte3DParser
  !--
  subroutine replaceInt1Byte3DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    integer (kind=1), dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    integer (kind=1), dimension (:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var,1),size(var,2),1))
    buf(:,:,1) = var
    call replaceInteger1Byte3D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine replaceInt1Byte3DParser
  !--
  subroutine appendInt1Byte3DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    integer (kind=1), dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    integer (kind=1), dimension (:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var,1),size(var,2),1))
    buf(:,:,1) = var

    call appendInteger1Byte3D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendInt1Byte3DParser
  !--
  subroutine writeInt1Byte4DParser(fid, var, varName, varLenName, &
       varUnit, varLongName, status, unlimited_dim)
    integer*1, dimension(:,:,:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varUnit
    character (len=*), intent (in) :: status
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    integer*1, dimension (:,:,:,:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=4
    else
       unlimited=unlimited_dim
    endif

    if (.not. allocated(buf)) &
         allocate(buf(size(var,1),size(var,2),size(var,3),1))
    buf(:,:,:,1) = var
    call writeInteger1Byte4D(fid, buf, &
         varName=varName,&
         varLenName=varLenname, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine writeInt1Byte4DParser
  !--
  subroutine replaceInt1Byte4DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    integer (kind=1), dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    integer (kind=1), dimension (:,:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) &
         allocate(buf(size(var,1),size(var,2),size(var,3),1))
    buf(:,:,:,1) = var
    call replaceInteger1Byte4D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine replaceInt1Byte4DParser
  !--
  subroutine appendInt1Byte4DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    integer (kind=1), dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    integer (kind=1), dimension (:,:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) &
         allocate(buf(size(var,1),size(var,2),size(var,3),1))
    buf(:,:,:,1) = var

    call appendInteger1Byte4D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendInt1Byte4DParser

  !--
  subroutine writeInt2Byte1DParser(fid, var, varName, varLenName, &
       varUnit, varLongName, status, unlimited_dim)
    integer (kind=2), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varLenName, varUnit
    character (len=*), intent (in) :: status
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    integer (kind=2), dimension (:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=1
    else
       unlimited=unlimited_dim
    endif

    if (.not. allocated(buf)) allocate(buf(1))
    buf(1) = var
    call writeInteger2Byte1D(fid, buf, &
         varName=varName,&
         varLenName=varLenname, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine writeInt2Byte1DParser
  !--
  subroutine replaceInt2Byte1DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    integer (kind=2), intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    integer (kind=2), dimension (:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(1))
    buf(1) = var
    call replaceInteger2Byte1D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)
  end subroutine replaceInt2Byte1DParser
  !--
  subroutine appendInt2Byte1DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    integer (kind=2), intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    integer (kind=2), dimension (:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(1))
    buf(1) = var

    call appendInteger2Byte1D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendInt2Byte1DParser
  !--
  subroutine writeInt2Byte2DParser(fid, var, varName, varLenName, &
       varUnit, varLongName, status, unlimited_dim)
    integer (kind=2), dimension(:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varUnit
    character (len=*), intent (in) :: status
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    integer (kind=2), dimension (:,:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=2
    else
       unlimited=unlimited_dim
    endif

    if (.not. allocated(buf)) allocate(buf(size(var),1))
    buf(:,1) = var
    call writeInteger2Byte2D(fid, buf, &
         varName=varName,&
         varLenName=varLenname, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine writeInt2Byte2DParser
  !--
  subroutine replaceInt2Byte2DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    integer (kind=2), dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    integer (kind=2), dimension (:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var),1))
    buf(:,1) = var
    call replaceInteger2Byte2D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)
  end subroutine replaceInt2Byte2DParser
  !--
  subroutine appendInt2Byte2DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    integer (kind=2), dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    integer (kind=2), dimension (:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var),1))
    buf(:,1) = var

    call appendInteger2Byte2D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendInt2Byte2DParser
  !--
  subroutine writeInt2Byte3DParser(fid, var, varName, varLenName, &
       varUnit, varLongName, status, unlimited_dim)
    integer (kind=2), dimension(:,:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varUnit
    character (len=*), intent (in) :: status
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    integer (kind=2), dimension (:,:,:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=3
    else
       unlimited=unlimited_dim
    endif

    if (.not. allocated(buf)) &
         allocate(buf(size(var,1),size(var,2),1))
    buf(:,:,1) = var
    call writeInteger2Byte3D(fid, buf, &
         varName=varName,&
         varLenName=varLenname, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine writeInt2Byte3DParser
  !--
  subroutine replaceInt2Byte3DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    integer (kind=2), dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    integer (kind=2), dimension (:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var,1),size(var,2),1))
    buf(:,:,1) = var
    call replaceInteger2Byte3D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine replaceInt2Byte3DParser
  !--
  subroutine appendInt2Byte3DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    integer (kind=2), dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    integer (kind=2), dimension (:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var,1),size(var,2),1))
    buf(:,:,1) = var

    call appendInteger2Byte3D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendInt2Byte3DParser
  !--
  subroutine writeInt2Byte4DParser(fid, var, varName, varLenName, &
       varUnit, varLongName, status, unlimited_dim)
    integer*2, dimension(:,:,:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varUnit
    character (len=*), intent (in) :: status
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    integer*2, dimension (:,:,:,:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=4
    else
       unlimited=unlimited_dim
    endif

    if (.not. allocated(buf)) &
         allocate(buf(size(var,1),size(var,2),size(var,3),1))
    buf(:,:,:,1) = var
    call writeInteger2Byte4D(fid, buf, &
         varName=varName,&
         varLenName=varLenname, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine writeInt2Byte4DParser
  !--
  subroutine replaceInt2Byte4DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    integer (kind=2), dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    integer (kind=2), dimension (:,:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) &
         allocate(buf(size(var,1),size(var,2),size(var,3),1))
    buf(:,:,:,1) = var
    call replaceInteger2Byte4D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine replaceInt2Byte4DParser
  !--
  subroutine appendInt2Byte4DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    integer (kind=2), dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    integer (kind=2), dimension (:,:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) &
         allocate(buf(size(var,1),size(var,2),size(var,3),1))
    buf(:,:,:,1) = var

    call appendInteger2Byte4D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendInt2Byte4DParser
  !--
  subroutine writeFlt1DParser(fid, var, varName, varLenName, &
       varUnit, varLongName, status, unlimited_dim)
    real*4, intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varLenName, varUnit
    character (len=*), intent (in) :: status
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    real*4, dimension (:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=1
    else
       unlimited=unlimited_dim
    endif

    if (.not. allocated(buf)) allocate(buf(1))
    buf(1) = var
    call writeFloat1D(fid, buf, &
         varName=varName,&
         varLenName=varLenname, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine writeFlt1DParser
  !--
  subroutine replaceFlt1DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    real*4, intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    real*4, dimension (:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(1))
    buf(1) = var
    call replaceFloat1D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)
  end subroutine replaceFlt1DParser
  !--
  subroutine appendFlt1DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    real*4, intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    real*4, dimension (:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(1))
    buf(1) = var

    call appendFloat1D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendFlt1DParser
  !--
  subroutine writeFlt2DParser(fid, var, varName, varLenName, &
       varUnit, varLongName, status, unlimited_dim)
    real*4, dimension(:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varUnit
    character (len=*), intent (in) :: status
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    real*4, dimension (:,:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=2
    else
       unlimited=unlimited_dim
    endif
    if (.not. allocated(buf)) allocate(buf(size(var),1))

    buf(:,1) = var
    call writeFloat2D(fid, buf,&
         varName=varName,&
         varLenName=varLenName, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)

    if (allocated(buf)) deallocate(buf)

  end subroutine writeFlt2DParser
  !--
  subroutine replaceFlt2DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    real*4, dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    real*4, dimension (:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var),1))
    buf(:,1) = var
    call replaceFloat2D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)
  end subroutine replaceFlt2DParser
  !--
  subroutine appendFlt2DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    real*4, dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    real*4, dimension (:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var),1))
    buf(:,1) = var

    call appendFloat2D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendFlt2DParser
  !--
  subroutine writeFlt3DParser(fid, var, varName, varLenName, &
       varUnit, varLongName, status, unlimited_dim)
    real*4, dimension(:,:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varUnit
    character (len=*), intent (in) :: status
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    real*4, dimension (:,:,:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=3
    else
       unlimited=unlimited_dim
    endif
    if (.not. allocated(buf)) &
         allocate(buf(size(var,1),size(var,2),1))
    buf(:,:,1) = var
    call writeFloat3D(fid, buf, &
         varName=varName,&
         varLenName=varLenname, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine writeFlt3DParser
  !--
  subroutine replaceFlt3DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    real*4, dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    real*4, dimension (:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var,1),size(var,2),1))
    buf(:,:,1) = var
    call replaceFloat3D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine replaceFlt3DParser
  !--
  subroutine appendFlt3DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    real*4, dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    real*4, dimension (:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var,1),size(var,2),1))
    buf(:,:,1) = var

    call appendFloat3D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendFlt3DParser
  !--
  subroutine writeFlt4DParser(fid, var, varName, varLenName, &
       varUnit, varLongName, status, unlimited_dim)
    real*4, dimension(:,:,:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varUnit
    character (len=*), intent (in) :: status
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    real*4, dimension (:,:,:,:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=4
    else
       unlimited=unlimited_dim
    endif
    if (.not. allocated(buf)) &
         allocate(buf(size(var,1),size(var,2),size(var,3),1))
    buf(:,:,:,1) = var

    call writeFloat4D(fid, buf, &
         varName=varName,&
         varLenName=varLenname, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine writeFlt4DParser
  !--
  subroutine replaceFlt4DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    real*4, dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    real*4, dimension (:,:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) &
         allocate(buf(size(var,1),size(var,2),size(var,3),1))
    buf(:,:,:,1) = var
    call replaceFloat4D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine replaceFlt4DParser
  !--
  subroutine appendFlt4DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    real*4, dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    real*4, dimension (:,:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) &
         allocate(buf(size(var,1),size(var,2),size(var,3),1))
    buf(:,:,:,1) = var

    call appendFloat4D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendFlt4DParser
  !--
  subroutine writeDFlt1DParser(fid, var, varName, varLenName, &
       varUnit, varLongName, status, unlimited_dim)
    double precision, intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varLenName, varUnit
    character (len=*), intent (in) :: status
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    double precision, dimension (:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=1
    else
       unlimited=unlimited_dim
    endif
    if (.not. allocated(buf)) allocate(buf(1))
    buf(1) = var
    call writeDFloat1D(fid, buf, &
         varName=varName,&
         varLenName=varLenname, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine writeDFlt1DParser
  !--
  subroutine replaceDFlt1DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    double precision, intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    double precision, dimension (:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(1))
    buf(1) = var
    call replaceDFloat1D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)
  end subroutine replaceDFlt1DParser
  !--
  subroutine appendDFlt1DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    double precision, intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    double precision, dimension (:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(1))
    buf(1) = var

    call appendDFloat1D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendDFlt1DParser
  !--
  subroutine writeDFlt2DParser(fid, var, varName, varLenName, &
       varUnit, varLongName, status, unlimited_dim)
    double precision, dimension(:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varUnit
    character (len=*), intent (in) :: status
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    double precision, dimension (:,:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=2
    else
       unlimited=unlimited_dim
    endif

    if (.not. allocated(buf)) allocate(buf(size(var),1))
    buf(:,1) = var
    call writeDFloat2D(fid, buf,&
         varName=varName,&
         varLenName=varLenName, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited,&
         append=append)

    if (allocated(buf)) deallocate(buf)

  end subroutine writeDFlt2DParser
  !--
  subroutine replaceDFlt2DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    double precision, dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    double precision, dimension (:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var),1))
    buf(:,1) = var
    call replaceDFloat2D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)
  end subroutine replaceDFlt2DParser
  !--
  subroutine appendDFlt2DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    double precision, dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    double precision, dimension (:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var),1))
    buf(:,1) = var

    call appendDFloat2D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendDFlt2DParser
  !--
  subroutine writeDFlt3DParser(fid, var, varName, varLenName, &
       varUnit, varLongName, status, unlimited_dim)
    double precision, dimension(:,:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varUnit
    character (len=*), intent (in) :: status
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    double precision, dimension (:,:,:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=3
    else
       unlimited=unlimited_dim
    endif
    if (.not. allocated(buf)) &
         allocate(buf(size(var,1),size(var,2),1))
    buf(:,:,1) = var
    call writeDFloat3D(fid, buf, &
         varName=varName,&
         varLenName=varLenname, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine writeDFlt3DParser
  !--
  subroutine replaceDFlt3DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    double precision, dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    double precision, dimension (:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var,1),size(var,2),1))
    buf(:,:,1) = var
    call replaceDFloat3D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine replaceDFlt3DParser
  !--
  subroutine appendDFlt3DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    double precision, dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    double precision, dimension (:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var,1),size(var,2),1))
    buf(:,:,1) = var

    call appendDFloat3D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendDFlt3DParser
  !--
  subroutine writeDFlt4DParser(fid, var, varName, varLenName, &
       varUnit, varLongName, status, unlimited_dim)
    double precision, dimension(:,:,:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varUnit
    character (len=*), intent (in) :: status
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    double precision, dimension (:,:,:,:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=4
    else
       unlimited=unlimited_dim
    endif
    if (.not. allocated(buf)) &
         allocate(buf(size(var,1),size(var,2),size(var,3),1))
    buf(:,:,:,1) = var

    call writeDFloat4D(fid, buf, &
         varName=varName,&
         varLenName=varLenname, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine writeDFlt4DParser
  !--
  subroutine replaceDFlt4DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    double precision, dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    double precision, dimension (:,:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) &
         allocate(buf(size(var,1),size(var,2),size(var,3),1))
    buf(:,:,:,1) = var
    call replaceDFloat4D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine replaceDFlt4DParser
  !--
  subroutine appendDFlt4DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    double precision, dimension(:,:,:), intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    double precision, dimension (:,:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) &
         allocate(buf(size(var,1),size(var,2),size(var,3),1))
    buf(:,:,:,1) = var

    call appendDFloat4D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendDFlt4DParser
  !--
  subroutine writeStr1DParser(fid, var, varName, strLenName,&
       varLenName, varUnit, varLongName, status, unlimited_dim)
    character (len=*), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varLenName, varUnit
    character (len=*), intent (in), optional :: strLenName
    character (len=*), intent (in) :: status
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    character (len=len(var)), dimension (:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=1
    else
       unlimited=unlimited_dim
    endif
    if (.not. allocated(buf)) allocate(buf(1))
    buf(1) = var
    call writeString1D(fid, buf, &
         varName=varName,&
         varLenName=varLenName, &
         strLenName=strLenName, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine writeStr1DParser
  !--
  subroutine replaceStr1DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    character (len=*), intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    character (len=len(var)), dimension (:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(1))
    buf(1) = var
    call replaceString1D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)
  end subroutine replaceStr1DParser
  !--
  subroutine appendStr1DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    character (len=*), intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    character (len=len(var)), dimension (:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(1))
    buf(1) = var

    call appendString1D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendStr1DParser
  !--
  subroutine writeStr2DParser(fid, var, varName, strLenName,&
       varLenName, varUnit, varLongName, status, unlimited_dim)
    character (len=*), dimension(:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varUnit
    character (len=*), intent (in) :: status
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    character (len=*), intent (in), optional :: strLenName
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    character (len=len(var(1))), &
         dimension (:,:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=2
    else
       unlimited=unlimited_dim
    endif
    if (.not. allocated(buf)) allocate(buf(size(var),1))
    buf(:,1) = var
    call writeString2D(fid, buf, &
         varName=varName,&
         varLenName=varLenname, &
         strLenName=strLenname, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine writeStr2DParser
  !--
  subroutine replaceStr2DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    character (len=*), dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    character (len=len(var(1))), dimension (:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var),1))
    buf(:,1) = var
    call replaceString2D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)
  end subroutine replaceStr2DParser
  !--
  subroutine appendStr2DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    character (len=*), dimension(:), intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    character (len=len(var(1))), dimension (:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var),1))
    buf(:,1) = var

    call appendString2D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendStr2DParser
  !--
  subroutine writeStr3DParser(fid, var, varName, strLenName, &
       varLenName, varUnit, varLongName, status, unlimited_dim)
    character (len=*), dimension(:,:), intent (in) :: var
    character (len=*), intent (in), optional :: &
         varName, varLongName, varUnit
    character (len=*), intent (in), optional :: strLenName
    character (len=*), intent (in) :: status
    character (len=*), dimension(:), intent (in), optional :: &
         varLenName
    integer, intent (in) :: fid
    integer, intent (in), optional :: unlimited_dim
    integer :: unlimited
    character (len=len(var(1,1))), &
         dimension (:,:,:), allocatable :: buf
    logical :: append

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. present(unlimited_dim)) then
       unlimited=3
    else
       unlimited=unlimited_dim
    endif
    if (.not. allocated(buf)) &
         allocate(buf(size(var,1),size(var,2),1))

    buf(:,:,1) = var
    call writeString3D(fid, buf, &
         varName=varName,&
         varLenName=varLenname, &
         strLenName=strLenname, &
         varUnit=varUnit, &
         varLongName=varLongName,&
         unlimited_dim=unlimited, &
         append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine writeStr3DParser
  !--
  subroutine replaceStr3DParser(fid,var,varName,rec,status)
    integer, intent(in) :: fid
    character (len=*), dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname
    character (len=*), intent(in) :: status
    integer, intent(in), optional :: rec
    logical :: append
    character (len=len(var(1,1))), dimension (:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var,1),size(var,2),1))
    buf(:,:,1) = var
    call replaceString3D(fid=fid,var=buf,varName=varName,rec=rec,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine replaceStr3DParser
  !--
  subroutine appendStr3DParser(fid, var, varName, status)
    integer, intent(in) :: fid
    character (len=*), dimension(:,:), intent(in) :: var
    character (len=*), intent(in) :: varname, status
    logical :: append
    character (len=len(var(1,1))), dimension (:,:,:), allocatable :: buf

    if (status == 'append') then
       append = .true.
    else
       append = .false.
    endif
    if (.not. allocated(buf)) allocate(buf(size(var,1),size(var,2),1))
    buf(:,:,1) = var

    call appendString3D(fid=fid,var=buf,varName=varName,append=append)
    if (allocated(buf)) deallocate(buf)

  end subroutine appendStr3DParser
  !--
  subroutine readInt1DParser(fid, var, varName, varLen, &
       varUnit, record_no, status)
    integer, intent (in) :: fid
    integer, intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (out), optional :: varLen
    integer, intent (in) :: record_no
    character (len=*), intent (in) :: status
    integer, dimension(:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(1))

    call readInteger1D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var = myvar(1)
    deallocate(myvar)

  end subroutine readInt1Dparser
  !--
  subroutine readInt2DParser(fid, var, varName, varLen,&
       varUnit, record_no, status)
    integer, dimension(:), intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(2), intent (out), optional :: varLen
    character (len=*), intent (in) :: status
    integer, dimension(:,:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(size(var,1),1))

    call readInteger2D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var(:) = myvar(:,1)
    deallocate(myvar)

  end subroutine readInt2DParser
  !--
  subroutine readInt3DParser(fid, var, varName, varLen,&
       varUnit, record_no, status)
    integer, dimension(:,:), intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(3), intent (out), optional :: varLen
    character (len=*), intent (in) :: status
    integer, dimension(:,:,:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(size(var,1),size(var,2),1))

    call readInteger3D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var = myvar(:,:,1)
    deallocate(myvar)

  end subroutine readInt3DParser
  !--
  subroutine readInt4DParser(fid, var, varName, varLen,&
       varUnit, record_no, status)
    integer, dimension(:,:,:), intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(4), intent (out), optional :: varLen
    character (len=*), intent (in) :: status
    integer, dimension(:,:,:,:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(size(var,1),size(var,2),size(var,3),1))

    call readInteger4D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var = myvar(:,:,:,1)
    deallocate(myvar)

  end subroutine readInt4DParser
  !--
  subroutine readInt1Byte1DParser(fid, var, varName, varLen, &
       varUnit, record_no, status)
    integer, intent (in) :: fid
    integer (kind=1), intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (out), optional :: varLen
    integer, intent (in) :: record_no
    character (len=*), intent (in) :: status
    integer (kind=1), dimension(:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(1))

    call readInteger1Byte1D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var = myvar(1)
    deallocate(myvar)

  end subroutine readInt1Byte1Dparser
  !--
  subroutine readInt1Byte2DParser(fid, var, varName, varLen,&
       varUnit, record_no, status)
    integer (kind=1), dimension(:), intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(2), intent (out), optional :: varLen
    character (len=*), intent (in) :: status
    integer (kind=1), dimension(:,:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(size(var,1),1))

    call readInteger1Byte2D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var(:) = myvar(:,1)
    deallocate(myvar)

  end subroutine readInt1Byte2DParser
  !--
  subroutine readInt1Byte3DParser(fid, var, varName, varLen,&
       varUnit, record_no, status)
    integer (kind=1), dimension(:,:), intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(3), intent (out), optional :: varLen
    character (len=*), intent (in) :: status
    integer (kind=1), dimension(:,:,:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(size(var,1),size(var,2),1))

    call readInteger1Byte3D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var = myvar(:,:,1)
    deallocate(myvar)

  end subroutine readInt1Byte3DParser
  !--
  subroutine readInt1Byte4DParser(fid, var, varName, varLen,&
       varUnit, record_no, status)
    integer*1, dimension(:,:,:), intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(4), intent (out), optional :: varLen
    character (len=*), intent (in) :: status
    integer*1, dimension(:,:,:,:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(size(var,1),size(var,2),size(var,3),1))

    call readInteger1Byte4D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var = myvar(:,:,:,1)
    deallocate(myvar)

  end subroutine readInt1Byte4DParser
  !--
  subroutine readInt2Byte1DParser(fid, var, varName, varLen, &
       varUnit, record_no, status)
    integer, intent (in) :: fid
    integer (kind=2), intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (out), optional :: varLen
    integer, intent (in) :: record_no
    character (len=*), intent (in) :: status
    integer (kind=2), dimension(:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(1))

    call readInteger2Byte1D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var = myvar(1)
    deallocate(myvar)

  end subroutine readInt2Byte1Dparser
  !--
  subroutine readInt2Byte2DParser(fid, var, varName, varLen,&
       varUnit, record_no, status)
    integer (kind=2), dimension(:), intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(2), intent (out), optional :: varLen
    character (len=*), intent (in) :: status
    integer (kind=2), dimension(:,:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(size(var,1),1))

    call readInteger2Byte2D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var(:) = myvar(:,1)
    deallocate(myvar)

  end subroutine readInt2Byte2DParser
  !--
  subroutine readInt2Byte3DParser(fid, var, varName, varLen,&
       varUnit, record_no, status)
    integer (kind=2), dimension(:,:), intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(3), intent (out), optional :: varLen
    character (len=*), intent (in) :: status
    integer (kind=2), dimension(:,:,:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(size(var,1),size(var,2),1))

    call readInteger2Byte3D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var = myvar(:,:,1)
    deallocate(myvar)

  end subroutine readInt2Byte3DParser
  !--
  subroutine readInt2Byte4DParser(fid, var, varName, varLen,&
       varUnit, record_no, status)
    integer*2, dimension(:,:,:), intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(4), intent (out), optional :: varLen
    character (len=*), intent (in) :: status
    integer*2, dimension(:,:,:,:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(size(var,1),size(var,2),size(var,3),1))

    call readInteger2Byte4D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var = myvar(:,:,:,1)
    deallocate(myvar)

  end subroutine readInt2Byte4DParser
  !--
  subroutine readFlt1DParser(fid, var, varName, varLen, &
       varUnit, record_no, status)
    integer, intent (in) :: fid
    real*4, intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (out), optional :: varLen
    integer, intent (in), optional :: record_no
    character (len=*), intent (in) :: status
    real*4, dimension(:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(1))

    call readFloat1D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var = myvar(1)
    deallocate(myvar)

  end subroutine readFlt1Dparser
  !--
  subroutine readFlt2DParser(fid, var, varName, varLen,&
       varUnit, record_no, status)
    real*4, dimension(:), intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(2), intent (out), optional :: varLen
    character (len=*), intent (in) :: status
    real*4, dimension(:,:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(size(var),1))

    call readFloat2D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var = myvar(:,1)
    deallocate(myvar)

  end subroutine readFlt2DParser
  !--
  subroutine readFlt3DParser(fid, var, varName, varLen,&
       varUnit, record_no, status)
    real*4, dimension(:,:), intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(3), intent (out), optional :: varLen
    character (len=*), intent (in) :: status
    real*4, dimension(:,:,:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(size(var,1),size(var,2),1))

    call readFloat3D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var = myvar(:,:,1)
    deallocate(myvar)

  end subroutine readFlt3DParser
  !--
  subroutine readFlt4DParser(fid, var, varName, varLen,&
       varUnit, record_no, status)
    real*4, dimension(:,:,:), intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(4), intent (out), optional :: varLen
    character (len=*), intent (in) :: status
    real*4, dimension(:,:,:,:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(size(var,1),size(var,2),size(var,3),1))

    call readFloat4D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var = myvar(:,:,:,1)
    deallocate(myvar)

  end subroutine readFlt4DParser
  !--
  subroutine readDFlt1DParser(fid, var, varName, varLen, &
       varUnit, record_no, status)
    integer, intent (in) :: fid
    double precision, intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (out), optional :: varLen
    integer, intent (in) :: record_no
    character (len=*), intent (in) :: status
    double precision, dimension(:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(1))

    call readDFloat1D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var = myvar(1)
    deallocate(myvar)

  end subroutine readDFlt1Dparser
  !--
  subroutine readDFlt2DParser(fid, var, varName, varLen,&
       varUnit, record_no, status)
    double precision, dimension(:), intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(2), intent (out), optional :: varLen
    character (len=*), intent (in) :: status
    double precision, dimension(:,:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(size(var,1),1))

    call readDFloat2D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var = myvar(:,1)
    deallocate(myvar)

  end subroutine readDFlt2DParser
  !--
  subroutine readDFlt3DParser(fid, var, varName, varLen,&
       varUnit, record_no, status)
    double precision, dimension(:,:), intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(3), intent (out), optional :: varLen
    character (len=*), intent (in) :: status
    double precision, dimension(:,:,:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(size(var,1),size(var,2),1))

    call readDFloat3D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var = myvar(:,:,1)
    deallocate(myvar)

  end subroutine readDFlt3DParser
  !--
  subroutine readStr1DParser(fid, var, varName, varLen, &
       varUnit, record_no, status)
    integer, intent (in) :: fid
    character (len=*), intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (in) :: record_no
    integer, intent (out), optional :: varLen
    character (len=*), intent (in) :: status
    character (len=MAX_NAME_LENGTH), &
         dimension(:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(1))

    call readString1D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var = myvar(1)
    deallocate(myvar)

  end subroutine readStr1Dparser
  !--
  subroutine readStr2DParser(fid, var, varName, varLen,&
       varUnit, record_no, status)
    character (len=*), dimension(:), intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(2), intent (out), optional :: varLen
    character (len=*), intent (in) :: status
    character (len=MAX_NAME_LENGTH), &
         dimension(:,:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(size(var,1),1))

    call readString2D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var = myvar(:,1)
    deallocate(myvar)

  end subroutine readStr2DParser
  !--
  subroutine readStr3DParser(fid, var, varName, varLen,&
       varUnit, record_no, status)
    character (len=*), dimension(:,:), intent (out) :: var
    character (len=*), intent (in), optional :: varName
    character (len=*), intent (out), optional :: varUnit
    integer, intent (in) :: fid
    integer, intent (in), optional :: record_no
    integer, dimension(3), intent (out), optional :: varLen
    character (len=*), intent (in) :: status
    character (len=MAX_NAME_LENGTH), &
         dimension(:,:,:), allocatable :: myvar

    if (allocated(myvar)) deallocate(myvar)
    allocate(myvar(size(var,1),size(var,2),1))

    call readString3D(fid, myvar, &
         varName=varName,&
         varLen=varLen,&
         varUnit=varUnit, &
         record_no=record_no)

    var = myvar(:,:,1)
    deallocate(myvar)

  end subroutine readStr3DParser
  !--
  subroutine closeFile(fid, message)
    character (len=*), intent (in), optional :: message
    integer, intent (in) :: fid

    if (present(message)) then
       print *, '------------------------'
       print *, 'closing unit: ', fid
       print *, message(1:len_trim(message))
       print *, '------------------------'
    endif
    ncStatus = nf_close(fid)
    if (ncStatus .ne. nf_noerr) &
         call handle_err(ncStatus,'closing')

    do io = 1, ncid_cnt_w
       if (ncid(io) == fid) then
          iopen = io
          exit
       endif
    enddo

    !--   Clear the index slot so it can be reused
    ncid_index_active(iopen) = .false.
    ncid(iopen) = 0
    unlimited_dim_id(iopen) = 0
    myUnlimited_dim_name(iopen)=' '
    varid(:,iopen) = 0
    existed(:,iopen) = .false.
    record_no(:,iopen) = 1
  end subroutine closeFile
  !--
  subroutine handle_err(ncStatus,msg)
    integer :: ncStatus
    character (len=*), intent (in), optional :: msg
    character (len=256) :: lmsg

    if (present(msg)) then
       if (len(msg) > 256) print*,'message is truncated'
       lmsg=trim(msg)
    else
       lmsg=':::'
    endif

    SELECT CASE (ncStatus)
    CASE (2)
       print *, 'err[ncdf_module::handle_err]: ',&
            ' ('//trim(msg)//') '//trim(nf_strerror(ncStatus))
       call exit(1)
    CASE DEFAULT
       print *, 'warning[ncdf_module::handle_err]: ',&
            ' ('//trim(msg)//') ',trim(nf_strerror(ncStatus)),ncStatus
    END SELECT
  end subroutine handle_err
  !--
  integer function FindString(name, nameArray)
    character (len=*), intent (in) :: name
    character (len=*), dimension (:), intent (in) :: nameArray
    integer :: i

    FindString = 0
    do i = 1, size(nameArray)
       if (trim(name) == trim(nameArray(i))) then
          FindString = i
          exit
       endif
    enddo
  end function FindString
  !--
  subroutine configureName(fid,iopen,define,varName)
    integer, intent (in) :: fid
    character (len=*), intent (in), optional :: varName
    logical, intent (out) :: define
    integer :: index
    integer, intent (out) :: iopen
    integer, save :: total_name=0

    do io = 1, ncid_cnt_w
       if (ncid(io) == fid) then
          iopen = io
          exit
       endif
    enddo

    if (present(varName)) then
       index = FindString(varName, tmpVarName)
       if (index == 0) then
          total_name = total_name + 1
          iname = total_name
          if (iname .gt. MAXVAR) then
             print *,'err[ncdf_module::configureName]: ',&
                  ' present(varName): iname overran MAXVAR'
             call exit(1)
          end if
          tmpVarName(iname) = varName
          define = .true.
          existed(iname,iopen) = .true.
       else
          iname = index
          if (existed(iname,iopen)) then
             define = .false.
          else
             define = .true.
             existed(iname,iopen) = .true.
          endif
       endif
    else
       total_name = total_name + 1
       iname = total_name
       if (iname .gt. MAXVAR) then
          print *,'err[ncdf_module::configureName]: ',&
               ' iname overran MAXVAR'
          call exit(1)
       end if
       tmpVarName(iname) = 'newVar'
       define = .true.
       existed(iname,iopen) = .true.
    endif
  end subroutine configureName
  !--
  subroutine getNewNcidIndex(ncid_index)

    implicit none

    integer, intent (out) :: ncid_index

    !--   Local variables
    integer :: io

    do io = 1, ncid_cnt_w
       if (.not. ncid_index_active(io)) then
          ncid_index = io
          ncid_index_active(io) = .true.
          return
       endif
    enddo

    ncid_cnt_w = ncid_cnt_w + 1
    if (ncid_cnt_w .gt. MAXFILE) then
       print *,'err[ncdf_module::getNewNcidIndex]: ',&
            ' ncid_cnt_w overran MAXFILE'
       call exit(1)
    end if
    ncid_index = ncid_cnt_w
    ncid_index_active(ncid_index) = .true.
    return
  end subroutine getNewNcidIndex
  !--
  integer function getNcidIndex(fid)

    implicit none

    integer, intent (in) :: fid
    ! integer :: getNcIdIndex

    !--   Local variables
    integer :: io
    logical :: found
    character (len=MAX_NAME_LENGTH) :: procName

    procName = "[ncdf_module::getNcidIndex]:"
    found = .FALSE.
    do io = 1, ncid_cnt_w
       if (fid == ncid(io)) then
          getNcIdIndex = io
          found = .TRUE.
          return
       endif
    enddo

    if (.not. found) then
       print *, trim(procName)//' index not found... returning 999'
       getNcIdIndex = 999
       return
    end if
  end function getNcidIndex
  !--
  logical function varExists(fid,varName)

    implicit none

    integer, intent (in) :: fid
    character (len=*), intent (in) :: varName

    !--   Local variables
    character (len=MAX_NAME_LENGTH) :: procName
    integer :: tmpvarid

    procName = "[ncdf_module::varExists]:"

    ncStatus = nf_inq_varid(fid,trim(varName),tmpvarid)
    if (ncStatus == nf_noerr) then
       varExists=.TRUE.
    else
       varExists=.FALSE.
    endif

  end function varExists
end module ncdf_module


