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

MODULE OriginalFile

! <f90Module>***********************************************************
!
! NAME:
!
!   OriginalFile
!
! PURPOSE:
!
!   Get various information from the original ncdf file.
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>


  USE ncdf_module

  IMPLICIT NONE

  PUBLIC:: getDataScaleAndOffset,getDimensions

  PRIVATE:: getByteDataScaleAndOffset, &
            getIntegerDataScaleAndOffset, &
            getCharacterVariableAttribute, &
            getRealVariableAttribute


  INTERFACE getDataScaleAndOffset 
    MODULE PROCEDURE getByteDataScaleAndOffset
    MODULE PROCEDURE getIntegerDataScaleAndOffset
  END INTERFACE 

  INTERFACE getVariableAttribute
    MODULE PROCEDURE getCharacterVariableAttribute
    MODULE PROCEDURE getRealVariableAttribute
  END INTERFACE

  CONTAINS

    ! --------------------------------------------------------------------
    ! This subroutine reads in byte type of data array from the
    ! original file specified in fname. It returns the data and optionally
    ! returns scale, offset, unit, valid_min, valid_max and missing_value.
    ! --------------------------------------------------------------------
    SUBROUTINE getByteDataScaleAndOffset(ncid,fname,byteArray,varName, &
               theOffset,theScale,theUnit,theMin,theMax,theMissing)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   getByteDataScaleAndOffset
!
! PURPOSE:
!
!   Read byte-type data array and return it. Optionally return the 
!   data attributes.
!
! SYNTAX:
!
!   CALL getByteDataScaleAndOffset(ncid, fname, byteArray, varName, 
!      theOffset, theScale, theUnit, theMin, theMax, theMissing)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   fname       CHAR     File name
!   varName     CHAR     Name of variable
!   
!   INPUTS/OUTPUTS:
!   
!   ncid        INTEGER  ID for ncdf file 
!   byteArray   INTEGER  Generic byte array
!   theOffset   REAL     Offset attribute
!   theScale    REAL     Scale factor attribute
!   theUnit*    CHAR     Unit attribute
!   theMin*     REAL     Miniumum value attribute
!   theMax*     REAL     Maximum value attribute
!   theMissing* REAL     Missing value attribute
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    ! The I/O parameters
      INTEGER,                                     INTENT(INOUT)          :: ncid
      CHARACTER(LEN=*),                            INTENT(IN)             :: fname,varName
      INTEGER(SELECTED_INT_KIND(1)),DIMENSION(:,:),INTENT(INOUT)          :: byteArray
      REAL,                                        INTENT(INOUT)          :: theOffset,theScale 
      CHARACTER(LEN=*),                            INTENT(INOUT),OPTIONAL :: theUnit
      REAL,                                        INTENT(INOUT),OPTIONAL :: theMin,theMax,theMissing
    ! Local variables
      INTEGER     :: varID,ncStatus

      CALL openNcdfFile(ncid,fname,status='old')
      ncStatus = nf_inq_varid(ncid,trim(varName),varID)
      IF (ncStatus .NE. nf_noErr) THEN
         CALL handle_nc_err(ncStatus,"getDataScaleAndOffset::" &
           // varName // " not found")
      ENDIF
      ncStatus = nf_get_var_int1(ncid,varID,byteArray)
      IF (ncStatus .NE. nf_noErr) &
         CALL handle_nc_err(ncStatus,"getDataScaleAndOffset::" &
             // "Failed to read data.")
      !CALL readNcdfData(ncid,byteArray,varName)
      CALL getVariableAttribute(ncid,varName,'offset',theOffset)
      CALL getVariableAttribute(ncid,varName,'scale',theScale)
      IF (PRESENT(theUnit)) &
         CALL getVariableAttribute(ncid,varName,'units',theUnit)
      IF (PRESENT(theMin)) &
         CALL getVariableAttribute(ncid,varName,'valid_min',theMin)
      IF (PRESENT(theMax)) &
         CALL getVariableAttribute(ncid,varName,'valid_max',theMax)
      IF (PRESENT(theMissing)) &
         CALL getVariableAttribute(ncid,varName,'missing_value',theMissing)
      CALL closeNcdfFile(ncid)
    END SUBROUTINE getByteDataScaleAndOffset


    ! ------------------------------------------------------------------
    ! Gets the attributes of the variable specified in varName that
    ! are real values.
    ! ------------------------------------------------------------------
    SUBROUTINE getRealVariableAttribute(ncid,varName,attName,attrib)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   getRealVariableAttribute
!
! PURPOSE:
!
!   Get a real-type attribute of the variable specified in varName.
!
! SYNTAX:
!
!   CALL getRealVariableAttribute(ncid, varName, attName, attrib)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   varName  CHAR     Name of variable
!   attName  CHAR     attribute name
!   
!   INPUTS/OUTPUTS:
!   
!   ncid     INTEGER  ID for ncdf file 
!   attrib   REAL     value of an attribue
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    ! The I/O variables
      INTEGER,          INTENT(INOUT) :: ncid
      CHARACTER(LEN=*), INTENT(IN)    :: varName,attName
      REAL,             INTENT(INOUT) :: attrib

    ! Local variables
      INTEGER   :: varID,ncStatus,status1

      ncStatus = nf_inq_varid(ncid,trim(varName),varID)

      IF (ncStatus .NE. nf_noErr) THEN
         PRINT *, "Could not locate ", varName
      ELSE
         status1 = nf_get_att_real(ncid,varID,trim(attName),attrib)
         IF (status1 .NE. nf_noErr) THEN
            PRINT *, "Error reading ", attName, " for ",varName
            CALL errorHalt(1)
         ENDIF
      ENDIF
    END SUBROUTINE getRealVariableAttribute


    ! ------------------------------------------------------------------
    ! Gets the attributes of the variable specified in varName that
    ! are string values.
    ! ------------------------------------------------------------------

    SUBROUTINE getCharacterVariableAttribute(ncid,varName,attName,attrib)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   getCharacterVariableAttribute
!
! PURPOSE:
!
!   Get a string-type attribute of the variable specified in varName.
!
! SYNTAX:
!
!   CALL getCharacterVariableAttribute(ncid, varName, attName, attrib)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   varName  CHAR     Name of variable
!   attName  CHAR     attribute name
!   
!   INPUTS/OUTPUTS:
!   
!   ncid     INTEGER  ID for ncdf file 
!   attrib   CHAR     value of an attribue
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    ! The I/O variables
      INTEGER,          INTENT(INOUT) :: ncid
      CHARACTER(LEN=*), INTENT(IN)    :: varName,attName
      CHARACTER(LEN=*), INTENT(INOUT) :: attrib

    ! Local variables
      INTEGER   :: varID,ncStatus,status1

      ncStatus = nf_inq_varid(ncid,trim(varName),varID)

      IF (ncStatus .NE. nf_noErr) THEN
         PRINT *, "Could not locate ", varName
      ELSE
         status1 = nf_get_att_text(ncid,varID,attName,attrib)
         IF (status1 .NE. nf_noErr) THEN
            PRINT *, "Error reading ", attName, " for ",varName
            CALL errorHalt(1)
         ENDIF
      ENDIF
    END SUBROUTINE getCharacterVariableAttribute

    ! --------------------------------------------------------------------
    ! This subroutine reads in integer type of data array from the
    ! original file specified in fname. It returns the data and optionally
    ! returns scale, offset, unit, valid_min, valid_max and missing_value.
    ! --------------------------------------------------------------------
    SUBROUTINE getIntegerDataScaleAndOffset(ncid,fname,intArray,attrName, &
               theOffset,theScale,theUnit,theMin,theMax,theMissing)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   getIntegerDataScaleAndOffset
!
! PURPOSE:
!
!   Read integer-type data array and return it. Optionally return the 
!   data attributes.
!
! SYNTAX:
!
!   CALL getIntegerDataScaleAndOffset(ncid, fname, intArray, 
!      attrName, theOffset, theScale, theUnit, theMin, theMax, 
!      theMissing)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   fname       CHAR     File name
!   attrName    CHAR     attribute name
!   
!   INPUTS/OUTPUTS:
!   
!   ncid        INTEGER  ID for ncdf file 
!   intArray    INTEGER  Generic integer array
!   theOffset   REAL     Offset attribute
!   theScale    REAL     Scale factor attribute
!   theUnit*    CHAR     Unit attribute
!   theMin*     REAL     Miniumum value attribute
!   theMax*     REAL     Maximum value attribute
!   theMissing* REAL     Missing value attribute
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    ! The I/O parameters
      INTEGER,                 INTENT(INOUT)           :: ncid
      CHARACTER(LEN=*),        INTENT(IN)              :: fname,attrName
      INTEGER, DIMENSION(:,:), INTENT(INOUT)           :: intArray
      REAL,                    INTENT(INOUT)           :: theOffset,theScale 
      CHARACTER(LEN=*),        INTENT(INOUT), OPTIONAL :: theUnit
      REAL,                    INTENT(INOUT), OPTIONAL :: theMin,theMax,theMissing
    ! Local variables
      INTEGER      :: ncStatus,varID

      CALL openNcdfFile(ncid,fname,status='old')
      ncStatus = nf_inq_varid(ncid,trim(attrName),varID)

      IF (ncStatus .EQ. nf_noErr) THEN
          ncStatus = nf_get_var_int(ncid,varID,intArray)
          IF (ncStatus .NE. nf_noErr) CALL handle_nc_err(ncStatus, &
            "OriginalFile::getDataScaleAndOffset: " // attrName // " not found!")

          ncStatus = nf_get_att_real(ncid,varID,'offset',theOffset)
          ncStatus = nf_get_att_real(ncid,varID,'scale',theScale)

          IF (PRESENT(theUnit)) ncStatus = nf_get_att_text(ncid,varID, &
               'units',theUnit)
          IF (PRESENT(theMin)) ncStatus = nf_get_att_real(ncid,varID, &
               'valid_min',theMin)
          IF (PRESENT(theMax)) ncStatus = nf_get_att_real(ncid,varID, &
               'valid_max',theMax)
          IF (PRESENT(theMissing)) ncStatus = nf_get_att_real(ncid,varID, &
               'missing_value',theMissing)
      ENDIF
     
      CALL closeNcdfFile(ncid)
    END SUBROUTINE getIntegerDataScaleAndOffset


    ! -------------------------------------------------------------------
    ! This subroutine returns the lines and samples dimensions of
    ! the data arrays in original data file specified in fname.
    ! -------------------------------------------------------------------
    SUBROUTINE getDimensions(ncid,fname,lines,samples)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   getDimensions
!
! PURPOSE:
!
!   Return the lines and samples dimensions of the data arrays in 
!   original data file specified in fname.
!
! SYNTAX:
!
!   CALL getDimensions(ncid, fname, lines, samples)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   fname    CHAR     File name
!   
!   INPUTS/OUTPUTS:
!   
!   ncid     INTEGER  ID for ncdf file 
!   lines    INTEGER  number of lines
!   samples  INTEGER  number of samples
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    ! THE I/O parameters
      INTEGER,          INTENT(INOUT) :: ncid
      CHARACTER(LEN=*), INTENT(IN)    :: fname
      INTEGER,          INTENT(INOUT) :: lines,samples

      CALL openNcdfFile(ncid,fname,status='old')
      lines = readNcdfDim(ncid,'lines')
      samples = readNcdfDim(ncid,'samples')
      CALL closeNcdfFile(ncid)
    END SUBROUTINE getDimensions

    SUBROUTINE handle_nc_err(ncStatus,msg)

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
!   CALL handle_nc_err(ncStatus, msg)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   ncStatus  INTEGER  Error status
!   msg       CHAR     Text of message
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    ! The I/O parameters
      INTEGER,          INTENT(IN) :: ncStatus
      CHARACTER(LEN=*), INTENT(IN) :: msg

      PRINT *, msg
      CALL errorHalt(1)
    END SUBROUTINE handle_nc_err

END MODULE OriginalFile
