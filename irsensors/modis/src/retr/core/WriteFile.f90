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

MODULE WriteFile

! <f90Module>***********************************************************
!
! NAME:
!
!   WriteFile
!
! PURPOSE:
!
!   Provide services to write the results of linear retrieval to a 
!   specified netCDF file.
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

    USE ncdf_module
    USE StateIndexModule, ONLY: maxMol,StateIndex_t

    IMPLICIT NONE

    PUBLIC:: writeRetrievalResults,addMoreAttributes

    PRIVATE:: writeFloatResults,writeIntResults,writeWarning


    INTERFACE writeRetrievalResults
       MODULE PROCEDURE writeFloatResults
       MODULE PROCEDURE writeIntResults
    END INTERFACE


   CONTAINS

     ! ----------------------------------------------------------------
     ! The WriteFile module provides the services to write the results
     ! of stage 2 retrieval to a specified netCDF file.
     ! ----------------------------------------------------------------



      SUBROUTINE writeFloatResults(fileID,xGroupGrid,latGrid,lonGrid, &
              clusteridGrid,satZenGrid,geogGrid,latUnit,lonUnit, &
              satZenUnit,geogUnit,nElements,nLines,nPar,IGR,NGR) 

!<f90Subroutine>********************************************************
!
! NAME:
!
!   writeFloatResults
!
! PURPOSE:
!
!   Write data as unscaled float variables
!
! SYNTAX:
!
!   CALL writeFloatResults(fileID, xGroupGrid, latGrid, lonGrid, 
!      clusteridGrid, satZenGrid, geogGrid, latUnit, lonUnit, 
!      satZenUnit, geogUnit, nElements, nLines, nPar, IGR, NGR)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   fileID         INTEGER             Index of file
!   xGroupGrid     REAL                Grid of retrieved parameters
!   latGrid        REAL                Grid of latitudes
!   lonGrid        REAL                Grid of longitudes
!   clusteridGrid  INTEGER             Grid of cluster indices
!   satZenGrid     REAL                Grid of satellite-zenith 
!                                      angles 
!   geogGrid       INTEGER             Grid of geography IDs
!   latUnit        CHAR                units for latitude
!   lonUnit        CHAR                units for longitude
!   satZenUnit     CHAR                units for satellite zenith 
!                                      angle 
!   geogUnit       CHAR                units for geography ID
!   nElements      INTEGER             Number of elements
!   nLines         INTEGER             number of lines
!   nPar           INTEGER             Total number of elements in 
!                                      retrieval vector 
!   IGR            TYPE(STATEINDEX_T)  Starting indices for 
!                                      retrieved sections of 
!                                      geophysical state vector
!   NGR            TYPE(STATEINDEX_T)  Number of elements for  
!                                      retrieved sections of 
!                                      geophysical state vector
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

      ! The parameters
         INTEGER, INTENT(IN)                 :: fileID
         REAL,                         DIMENSION(:,:),  INTENT(IN) :: latGrid
         REAL,                         DIMENSION(:,:),  INTENT(IN) :: lonGrid
         REAL,                         DIMENSION(:,:),  INTENT(IN) :: satZenGrid
         INTEGER(SELECTED_INT_KIND(4)),DIMENSION(:,:),  INTENT(IN) :: geogGrid
         INTEGER,                      DIMENSION(:,:),  INTENT(IN) :: clusteridGrid
         REAL,                         DIMENSION(:,:,:),INTENT(IN) :: xGroupGrid
         CHARACTER(LEN=10),                             INTENT(IN) :: latUnit
         CHARACTER(LEN=10),                             INTENT(IN) :: lonUnit
         CHARACTER(LEN=10),                             INTENT(IN) :: satZenUnit
         CHARACTER(LEN=10),                             INTENT(IN) :: geogUnit
         INTEGER,                                       INTENT(IN) :: nPar,nLines,nElements
         TYPE(StateIndex_t),                            INTENT(IN) :: IGR,NGR

         CALL writeNcdfDim(fileID,dimName='nLines',dimLen=nLines)
         CALL writeNcdfDim(fileID,dimName='nElements',dimLen=nElements)
         CALL writeNcdfDim(fileID,dimName='nPar',dimLen=nPar)
         CALL writeNcdfData(fileID,latGrid,varName='lat', &
                varLenName=(/'nElements', 'nLines   '/), &
                varLongName='Latitude data',varUnit=latUnit)
         CALL writeNcdfData(fileID,lonGrid,varName='lon', &
                varLenName=(/'nElements', 'nLines   '/), &
                varLongName='Longitude data',varUnit=lonUnit)
         CALL writeNcdfData(fileID,satZenGrid,varName='SATZEN', &
                varLenName=(/'nElements', 'nLines   '/), &
                varLongName='Satellite zenith angle data',varUnit=satZenUnit)
         CALL writeNcdfData(fileID,geogGrid,varName='GEOGRAPHY', &
                varLenName=(/'nElements', 'nLines   '/), &
                varLongName='The geography data',varUnit=geogUnit)
         CALL writeNcdfData(fileID,clusteridGrid,varName='ClusterID', &
                varLenName=(/'nElements', 'nLines   '/), &
                varLongName='ClusterID of the cluster the image point belonged to',&
                varUnit='unitless')
         CALL writeXGroupResults(fileID,xGroupGrid,nPar,IGR,NGR)
         
      END SUBROUTINE writeFloatResults

      
      SUBROUTINE writeIntResults(fileID,xGroupGrid,latGrid,lonGrid, &
              clusteridGrid,satZenGrid,geogGrid,latUnit,lonUnit, &
              satZenUnit,geogUnit,samples,lines,nPar,IGR,NGR) 

!<f90Subroutine>********************************************************
!
! NAME:
!
!   writeIntResults
!
! PURPOSE:
!
!   Write data as scaled integer variables
!
! SYNTAX:
!
!   CALL writeIntResults(fileID, xGroupGrid, latGrid, lonGrid, 
!      clusteridGrid, satZenGrid, geogGrid, latUnit, lonUnit, 
!      satZenUnit, geogUnit, samples, lines, nPar, IGR, NGR)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   fileID         INTEGER             Index of file
!   xGroupGrid     REAL                Grid of retrieved parameters
!   latGrid        INTEGER             Grid of latitudes
!   lonGrid        INTEGER             Grid of longitudes
!   clusteridGrid  INTEGER             Grid of cluster indices
!   satZenGrid     INTEGER             Grid of satellite-zenith 
!                                      angles 
!   geogGrid       INTEGER             Grid of geography IDs
!   latUnit        CHAR                units for latitude
!   lonUnit        CHAR                units for longitude
!   satZenUnit     CHAR                units for satellite zenith 
!                                      angle 
!   geogUnit       CHAR                units for geography ID
!   samples        INTEGER             number of samples
!   lines          INTEGER             number of lines
!   nPar           INTEGER             Total number of elements in 
!                                      retrieval vector 
!   IGR            TYPE(STATEINDEX_T)  Starting indices for 
!                                      retrieved sections of 
!                                      geophysical state vector
!   NGR            TYPE(STATEINDEX_T)  Number of elements for  
!                                      retrieved sections of 
!                                      geophysical state vector
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

      ! The parameters
         INTEGER,                                       INTENT(IN) :: fileID
         INTEGER, DIMENSION(:,:),                       INTENT(IN) :: latGrid
         INTEGER, DIMENSION(:,:),                       INTENT(IN) :: lonGrid
         INTEGER, DIMENSION(:,:),                       INTENT(IN) :: satZenGrid
         INTEGER(SELECTED_INT_KIND(1)),DIMENSION(:,:),  INTENT(IN) :: geogGrid
         INTEGER,                      DIMENSION(:,:),  INTENT(IN) :: clusteridGrid
         REAL,                         DIMENSION(:,:,:),INTENT(IN) :: xGroupGrid
         CHARACTER(LEN=10),                             INTENT(IN) :: latUnit
         CHARACTER(LEN=10),                             INTENT(IN) :: lonUnit
         CHARACTER(LEN=10),                             INTENT(IN) :: satZenUnit
         CHARACTER(LEN=10),                             INTENT(IN) :: geogUnit
         INTEGER,                                       INTENT(IN) :: nPar,lines,samples
         TYPE(StateIndex_t),                            INTENT(IN) :: IGR, NGR

         CALL writeNcdfDim(fileID,dimName='lines',dimLen=lines)
         CALL writeNcdfDim(fileID,dimName='samples',dimLen=samples)
         CALL writeNcdfDim(fileID,dimName='nPar',dimLen=nPar)
         CALL writeNcdfData(fileID,latGrid,varName='lat', &
                varLenName=(/'samples', 'lines  '/), &
                varLongName='Latitude data',varUnit=latUnit)
         CALL writeNcdfData(fileID,lonGrid,varName='lon', &
                varLenName=(/'samples', 'lines  '/), &
                varLongName='Longitude data',varUnit=lonUnit)
         CALL writeNcdfData(fileID,satZenGrid,varName='SATZEN', &
                varLenName=(/'samples', 'lines  '/), &
                varLongName='Satellite zenith angle data',varUnit=satZenUnit)
         CALL writeNcdfData(fileID,geogGrid,varName='GEOGRAPHY', &
                varLenName=(/'samples', 'lines  '/), &
                varLongName='The geography data',varUnit=geogUnit)
         CALL writeNcdfData(fileID,clusteridGrid,varName='ClusterID', &
                varLenName=(/'samples', 'lines  '/), &
                varLongName='ClusterID of the cluster the image point belonged to',&
                varUnit='unitless')
         CALL writeXGroupResults(fileID,xGroupGrid,nPar,IGR,NGR)
      END SUBROUTINE writeIntResults



      ! ---------------------------------------------------------------
      !  Adds additional float attributes to an already defined netCDF 
      ! variable.
      ! attNames: The new attribute names to be added
      ! attribs: Corresponding real values for the attributes
      ! variableName : The name of the variable for which the attributes 
      ! are added.
      ! ---------------------------------------------------------------
      SUBROUTINE addMoreAttributes(fileID,variableName,attNames,attribs)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   addMoreAttributes
!
! PURPOSE:
!
!   Adds additional float attributes to an already defined netCDF 
!   variable.
!
! SYNTAX:
!
!   CALL addMoreAttributes(fileID, variableName, attNames, attribs)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   fileID        INTEGER  Index of file
!   variableName  CHAR     variable name
!   attNames      CHAR     attribute names
!   attribs       REAL     attributes
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

      ! The parameters
         INTEGER,                        INTENT(IN) :: fileID
         CHARACTER(LEN=*),               INTENT(IN) :: variableName
         CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: attNames
         REAL, DIMENSION(:),             INTENT(IN) :: attribs
      ! Local variables
         INTEGER    :: i,variableID,ncStatus,status1,status2
   
         ncStatus = nf_redef(fileID)
         IF (ncStatus .EQ. nf_noErr) THEN
            status1 = nf_inq_varid(fileID,variableName,variableID)
         IF ((status1 .EQ. nf_noErr) .AND. (SIZE(attNames) .EQ. SIZE(attribs))) THEN
             DO i = 1,SIZE(attribs)
                status2 = nf_put_att_real(fileID,variableID,attNames(i),NF_FLOAT,1,&
                          attribs(i))
                IF (status2 .NE. nf_noErr) THEN
                    PRINT *, "Error writing ", attNames(i), " to ", variableName
                    CALL errorHalt(1)
                ENDIF
             ENDDO
          ENDIF
          ELSE
             PRINT *, "Failed to redefine ", fileID
          ENDIF
          ncStatus = nf_enddef(fileID)
       END SUBROUTINE addMoreAttributes

       ! -----------------------------------------------------------------
       ! Writes the retrieved xGroup parameters cldLiq and cldIce
       ! if both of these are available, or just the one that is available. 
       ! -----------------------------------------------------------------
       SUBROUTINE writeXGroupResults(fileID,xGroupGrid,nPar,IGR,NGR)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   writeXGroupResults
!
! PURPOSE:
!
!   Writes the linearly retrieved parameters to a file, translating 
!   from state vector index storage to individual variable names.
!
! SYNTAX:
!
!   CALL writeXGroupResults(fileID, xGroupGrid, nPar, IGR, NGR)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   fileID      INTEGER             Index of file
!   xGroupGrid  REAL                Grid of retrieved parameters
!   nPar        INTEGER             Total number of elements in 
!                                   retrieval vector 
!   IGR         TYPE(STATEINDEX_T)  Starting indices for retrieved 
!                                   sections of geophysical state 
!                                   vector 
!   NGR         TYPE(STATEINDEX_T)  Number of elements for 
!                                   retrieved sections of 
!                                   geophysical state vector
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

       ! The parameters
         INTEGER,                INTENT(IN) :: fileID,nPar
         REAL, DIMENSION(:,:,:), INTENT(IN) :: xGroupGrid
         TYPE(StateIndex_t),     INTENT(IN) :: IGR, NGR
       ! Local variables
         INTEGER              :: i
         LOGICAL              :: unknownParam
         INTEGER, PARAMETER   :: lenN=20
         INTEGER, PARAMETER   :: lenU=12
         CHARACTER(LEN=lenU)  :: nullName
         CHARACTER(LEN=lenN), DIMENSION(:), ALLOCATABLE :: parNames
         CHARACTER(LEN=lenU), DIMENSION(:), ALLOCATABLE :: parUnits

         CHARACTER(LEN=lenN),               PARAMETER :: namesTskin = &
          "SKIN_TEMPERATURE    "
         CHARACTER(LEN=lenN), DIMENSION(4), PARAMETER :: namesCldLiq = (/ &
          "LIQ_CLOUD_TOP_PRES  ","LIQ_CLOUD_THICK_PRES", &
          "LIQ_WATER_PATH      ","LIQ_CLOUD_EFFEC_DIAM"/)
         CHARACTER(LEN=lenN), DIMENSION(4), PARAMETER :: namesCldIce = (/ &
          "ICE_CLOUD_TOP_PRES  ","ICE_CLOUD_THICK_PRES", &
          "ICE_WATER_PATH      ","ICE_CLOUD_EFFEC_DIAM"/)
         CHARACTER(LEN=lenU),               PARAMETER :: unitsTskin = &
          "K           "
         CHARACTER(LEN=lenU), DIMENSION(4), PARAMETER :: unitsCldLiq = (/ &
          "mb          ","mb          ","kg/m^2      ","micrometers "/)
         CHARACTER(LEN=lenU), DIMENSION(4), PARAMETER :: unitsCldIce = (/ &
          "mb          ","mb          ","kg/m^2      ","micrometers "/)

         ALLOCATE(parNames(nPar),parUnits(nPar))
         nullName=REPEAT('*',lenN)
         parNames=nullName
       ! write the skin temperature if present
         IF (NGR%tskin == 1) THEN
            parNames(IGR%tskin) = namesTskin
            parUnits(IGR%tskin) = unitsTskin
         ELSEIF (NGR%tskin > 0) THEN
            CALL writeWarning("Tskin")
         ENDIF
      ! Write the cldLiq parameters if present
         IF ((NGR%cldLiq >= 1) .AND. (NGR%cldLiq <= 4)) THEN
            DO i = 1, NGR%cldLiq
               parNames(IGR%cldLiq+i-1) = namesCldLiq(i)
               parUnits(IGR%cldLiq+i-1) = unitsCldLiq(i)
            ENDDO
         ELSEIF (NGR%cldLiq > 0) THEN
            CALL writeWarning("cldLiq")
         ENDIF
      ! Write the cldIce parameters if present
         IF ((NGR%cldIce >= 1) .AND. (NGR%cldIce <= 4)) THEN
            DO i = 1, NGR%cldIce
               parNames(IGR%cldIce+i-1) = namesCldIce(i)
               parUnits(IGR%cldIce+i-1) = unitsCldIce(i)
            ENDDO
         ELSEIF (NGR%cldIce > 0) THEN
            CALL writeWarning("cldIce")
         ENDIF
       
         IF (NGR%temp .GT. 0) CALL writeWarning("temp")
         IF (NGR%psfc .GT. 0) CALL writeWarning("psfc")
         IF (NGR%wind .GT. 0) CALL writeWarning("wind")
         IF (SUM(NGR%mol(1:maxMol)) .GT. 0) CALL writeWarning("mol")

         DO i = 1, nPar
           IF (parNames(i) == nullName) CYCLE
           CALL WriteNcdfData(fileID,xGroupGrid(:,:,i), &
               varName=trim(parNames(i)), &
               varLenName=(/'samples', 'lines  '/), &
               varLongName="Retrieved "// trim(parNames(i)), &
               varUnit=trim(parUnits(i)))
         ENDDO 

         IF (ALLOCATED(parNames)) DEALLOCATE(parNames)
         IF (ALLOCATED(parUnits)) DEALLOCATE(parUnits)

       END SUBROUTINE writeXGroupResults

       SUBROUTINE writeWarning(procName)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   writeWarning
!
! PURPOSE:
!
!   Write a warning message.
!
! SYNTAX:
!
!   CALL writeWarning(procName)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   procName  CHAR  Name of variable, subject of warning message. 
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

       ! The parameters
          CHARACTER(LEN=*), INTENT(IN) :: procName

           PRINT *, "::Warning[", procName, &
                    "Do not know how to write the parameter!]"
       END SUBROUTINE writeWarning

 END MODULE writeFile
