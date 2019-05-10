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

PROGRAM LinRetrGroup

! <f90Main>*************************************************************
!
! NAME:
!
!   LinRetrGroup
!
! PURPOSE:
!
!   Main program for the linear second stage in two-stage retrieval
!
! INCLUDES:
!
!   None
!
!*************************************************************</f90Main>

    USE ncdf_module
    USE StateIndexModule, ONLY: maxMol,StateIndex_t,whereH2O
    USE LinInvert, ONLY: linRetr
    USE constants, ONLY: MISSING_REAL,MISSING_INT
    USE VertCoord, ONLY: mxCoordTyp
    USE RadFile, ONLY: getDimRad,getRadData, &
                       getGeogData
    USE ClusterFileIO, ONLY : getDimCluster,getClusterData
    USE CoefFile, ONLY: getCoefDimensions,getLinearCoefftsData
    USE OriginalFile
    USE WriteFile

    IMPLICIT NONE

    ! ------------------------------------------------------------
    ! This program reads in the radiance data, cluster data and
    ! linear coefficients data. Then it performs stage 2 retrieval
    ! and writes the retrieved results to a netCDF file.
    ! Some intermediate results are also written out to another
    ! netCDF file purely for diagnostic purpose.
    ! ----------------------------------------------------------

    ! define the buffers to read into
    CHARACTER(LEN=200) :: radFileName,clusFileName,linCoefFileName
    CHARACTER(LEN=200) :: originalFileName,tempFileName,stage2OutFile
    INTEGER            :: radFileID,clusFileID,linCoefFileID
    INTEGER            :: originalFileID,tempFileID,stage2OutID
    INTEGER            :: nChannels, nMembersTotal, nRec, nchan, &
                          clusChannels,clusMembersTotal !just error check
    INTEGER            :: nPar,nGroup,nChOn
    INTEGER            :: i,j,ih2o,nmol
    REAL               :: chiSqThresh
    REAL               :: QCfactor,waveThreshold
    CHARACTER(LEN=*), PARAMETER :: msg1 = "Wrote the read data and interim results to dummy file."
    CHARACTER(LEN=*), PARAMETER :: msg2 = "Wrote retrieval results."

    ! variables to be used in extracting only cloudy pixels' data
    ! using cluster info and rad file
    REAL, DIMENSION(:,:), ALLOCATABLE  :: cldRadiance
    REAL, DIMENSION(:), ALLOCATABLE    :: radWaveNumbers
    INTEGER, DIMENSION(:), ALLOCATABLE :: startIndices
    INTEGER, DIMENSION(:), ALLOCATABLE :: cldyMemberList
    INTEGER, DIMENSION(:), ALLOCATABLE :: memberList
    INTEGER, DIMENSION(:), ALLOCATABLE :: numMembers
    INTEGER                            :: nLines,nElements,maxClusSize
    ! variables to be used for reading lat, lon etc from original file
    ! instead of extracting cloudy pixels alone from radiance file.
    INTEGER                     :: lines,samples
    INTEGER, DIMENSION(:,:),ALLOCATABLE :: origLats, &
                                           origLons,origSatzen
    INTEGER(SELECTED_INT_KIND(1)), DIMENSION(:,:),ALLOCATABLE :: origGeog
    ! Additional attributes of netCDF variables
    REAL     :: latScale,lonScale,satZenScale,geogScale
    REAL     :: latOffset,lonOffset,satZenOffset,geogOffset
    REAL     :: latMin,lonMin,satZenMin,geogMin
    REAL     :: latMax,lonMax,satZenMax,geogMax
    REAL     :: latMissing,lonMissing,satZenMissing,geogMissing
    CHARACTER(LEN=20) :: latUnit,lonUnit,satZenUnit,geogUnit
    ! linear coefficients
    REAL, DIMENSION(:), ALLOCATABLE     :: linWaveNumbers
    REAL, DIMENSION(:,:), ALLOCATABLE   :: xOffset
    REAL, DIMENSION(:,:,:), ALLOCATABLE :: xGain
    REAL, DIMENSION(:), ALLOCATABLE     :: pSfc, tSfc
    REAL, DIMENSION(:), ALLOCATABLE     :: chiSqFinal
    REAL, DIMENSION(:), ALLOCATABLE     :: pressure
    REAL, DIMENSION(:,:),ALLOCATABLE    :: xGroup,clusRadGroup
    LOGICAL, DIMENSION(:), ALLOCATABLE  :: kchan
    INTEGER, DIMENSION(:), ALLOCATABLE  :: subscriptMap
    TYPE(StateIndex_t)                  :: IGR, NGR
    INTEGER                             :: nlev
    CHARACTER(LEN=mxCoordTyp)           :: vCoordTyp

    ! Computed results
    ! Grid in Satellite projection space for final retrieved results
    ! Both rad file and original data file represents lat,lon,satzen
    ! as shorts and GEOGRAPHY as Byte. Get offset and scale parameters.
    REAL, DIMENSION(:,:,:), ALLOCATABLE  :: xGroupGrid
    INTEGER, DIMENSION(:,:), ALLOCATABLE  :: clusteridGrid

    ! some dummys to be used in RadFile::getDimRad()
    ! and RadFile::getRadData()
    INTEGER    :: dupAng
    REAL, DIMENSION(:), ALLOCATABLE    :: cldyLats,cldyLons,cldySatZens
    ! some buffers for writing results
    CHARACTER(LEN=15), DIMENSION(5) :: nameArray = &
            (/'scale          ','offset         ','valid_min      ', &
            'valid_max      ','missing_value  '/)
    REAL, DIMENSION(5)                :: valueArray
    ! Default values for scale, offset,units and such if not read in
    latScale    = 1.0
    lonScale    = 1.0
    satZenScale = 1.0
    geogScale   = 1.0
    latOffset   = 0.0
    lonOffset   = 0.0
    satZenOffset = 0.0
    geogOffset   = 0.0
    latMin = 0.0
    lonMin = 0.0
    satZenMin = 0.0
    geogMin   = 0
    latMax = 90.0
    lonMax = 180.
    satZenMax = 180.0
    geogMax = 4
    latMissing = MISSING_REAL
    lonMissing = MISSING_REAL
    satZenMissing = MISSING_REAL
    geogMissing = 5
    latUnit = "degrees"
    lonUnit = "degrees"
    satZenUnit = "degrees"
    geogUnit   = "Unknown"

    CALL readConfigParameters(radFileName,clusFileName,tempFileName, &
         linCoefFileName,originalFileName,stage2OutFile, &
         QCfactor,waveThreshold)
    CALL openNcdfFile(tempFileID, tempFileName,status='new', &
         unlimited_dim_name='nProfiles')
    CALL openNcdfFile(stage2OutID, stage2OutFile,status='new',&
         unlimited_dim_name='nProfiles')
    CALL getDimRad(radFileID,radFileName,nChannels,nMembersTotal,dupAng)
    ! get necessary dimensions from cluster map file
    CALL getDimCluster(clusFileID,clusFileName,clusChannels,nRec, &
         clusMembersTotal,nLines,nElements)
    ! get necessary dimensions from linear coeffts file
    CALL getCoefDimensions(linCoefFileID,linCoefFileName,nPar,nGroup,nChOn, &
                 nchan,nlev)
    ! Get the dimensions from original data file grid sizes of latitude,
    ! longitude etc instead of cloudy only values.
    CALL getDimensions(originalFileID,originalFileName,lines,samples)

    IF (nMembersTotal .NE. clusMembersTotal) THEN
       PRINT *, "Dimensions did not match for radiance and cluster map files."
       CALL errorHalt(1)
    ENDIF
    IF ((samples .NE. nElements) .OR. (lines .NE. nLines)) THEN
       PRINT *, "Grid dimensions of original file did not match that ", &
               "of rad file."
       CALL errorHalt(1)
    ENDIF
    IF (nGroup .NE. nRec) THEN
       PRINT *, "nGroup in coeffts file did not match nRec in cluster data."
       CALL errorHalt(1)
    ENDIF
    ! Memory allocation for cluster and radiance input data
    ALLOCATE(cldRadiance(nChannels, 0:nMembersTotal-1) )
    ALLOCATE(radWaveNumbers(nChannels))
    ALLOCATE(startIndices(nRec))
    ALLOCATE(numMembers(nRec))
    ALLOCATE(cldyMemberList(nMembersTotal))
    ALLOCATE(memberList(nMembersTotal))
    ALLOCATE(xOffset(nPar,nGroup))
    ALLOCATE(xGain(nPar,nChOn,nGroup))
    ALLOCATE(pSfc(nGroup))
    ALLOCATE(tSfc(nGroup))
    ALLOCATE(kchan(nchan))
    ALLOCATE(linWaveNumbers(nchan))
    ALLOCATE(chiSqFinal(nGroup))
    ALLOCATE(pressure(nlev))
    ALLOCATE(subscriptMap(nchOn))
    ALLOCATE(origLats(samples,lines))
    ALLOCATE(origLons(samples,lines))
    ALLOCATE(origSatzen(samples,lines))
    ALLOCATE(origGeog(samples,lines))
    ALLOCATE(cldyLats(nMembersTotal))
    ALLOCATE(cldyLons(nMembersTotal))
    ALLOCATE(cldySatZens(nMembersTotal))

    ! Allocate the final output Grid arrays
    ALLOCATE(clusteridGrid(nElements,nLines))
    ALLOCATE(xGroupGrid(nElements,nLines,nPar))

    xGroupGrid = 0.
    clusteridGrid = MISSING_INT
    print *, "nPar: ", nPar, " nGroup: ", nGroup, " nChOn: ", nChOn
    CALL getRadData(radFileID,radFileName,radAll=cldRadiance, &
                    eiaAll=cldySatZens,latAll=cldyLats,lonAll=cldyLons, &
                    wvn=radWaveNumbers)
    CALL getClusterData(clusFileID,clusFileName,cldyMemberList,startIndices, &
         numMembers,memberList)

    maxClusSize = maxVal(numMembers)
    print *, "Read cluster map data."
    ALLOCATE(clusRadGroup(nChOn,maxClusSize))
    ALLOCATE(xGroup(nPar,maxClusSize))

    CALL getLinearCoefftsData(linCoefFileID,linCoefFileName,xOffset, &
            xGain,pSfc,tSfc,chiSqFinal,kchan,chiSqThresh, &
            vCoordTyp,pressure,linWaveNumbers,IGR,NGR)
    print *, "Read  coefft data"
    ih2o=1    ! assume water vapor is in the first slot, because MolID is not available

    CALL createSubscriptMap(kchan,nchan,nChannels,radWaveNumbers, &
                 linWaveNumbers,waveThreshold,subscriptMap)
    print *, "Created subscript map."
    print *, "QCfactor (const): ", QCfactor

    CALL drvRetr(nRec,nPar,startIndices,subscriptMap,numMembers, &
           cldyMemberList,memberList,cldRadiance,chiSqFinal,chiSqThresh, &
           QCfactor,xOffset(1:nPar,1:nGroup),xGain(1:nPar,1:nChOn,1:nGroup), &
           pSfc,tSfc,vCoordTyp,pressure,ih2o,IGR,NGR,clusRadGroup, &
           xGroup(1:nPar,:),xGroupGrid,clusteridGrid)

    ! Read latitude,longitude,satzen, and geography from original file
    CALL getDataScaleAndOffset(originalFileID,originalFileName,origLats, &
          'LATITUDE',latOffset,latScale,theUnit=latUnit,theMin=latMin, &
            theMax=latMax,theMissing=latMissing)
    CALL getDataScaleAndOffset(originalFileID,originalFileName,origLons, &
           'LONGITUDE',lonOffset,lonScale,theUnit=lonUnit,theMin=lonMin, &
             theMax=lonMax,theMissing=lonMissing)
    CALL getDataScaleAndOffset(originalFileID,originalFileName,origSatZen, &
           'SATZEN',satZenOffset,satZenScale,theUnit=satZenUnit, &
             theMin=satZenMin,theMax=satZenMax,theMissing=satZenMissing)
    CALL getDataScaleAndOffset(originalFileID,originalFileName,origGeog, &
           'GEOGRAPHY',geogOffset,geogScale,theUnit=geogUnit, &
             theMin=geogMin,theMax=geogMax,theMissing=geogMissing)

    CALL writeRetrievalResults(stage2OutID,xGroupGrid,origLats, &
          origLons,clusteridGrid,origSatZen,origGeog, &
          latUnit,lonUnit,satZenUnit,geogUnit,samples,lines,nPar, &
          IGR,NGR)

    valueArray = (/latScale,latOffset,latMin,latMax,latMissing/)
    CALL addMoreAttributes(stage2OutID,'lat',nameArray,valueArray)
    valueArray = (/lonScale,lonOffset,lonMin,lonMax,lonMissing/)
    CALL addMoreAttributes(stage2OutID,'lon',nameArray,valueArray)
    valueArray = (/satZenScale,satZenOffset,satZenMin,satZenMax,satZenMissing/)
    CALL addMoreAttributes(stage2OutID,'SATZEN',nameArray,valueArray)
    valueArray = (/geogScale,geogOffset,geogMin,geogMax,geogMissing/)
    CALL addMoreAttributes(stage2OutID,'GEOGRAPHY',nameArray,valueArray)
    CALL writeToFile(tempFileID, nMembersTotal,nChannels,nRec,cldyMemberList, &
         cldRadiance,radWaveNumbers,startIndices, &
         numMembers,linWaveNumbers,xOffset,xGain, &
         pSfc,tSfc,nchan,kchan,nPar,nGroup,nChOn, &
         chiSqThresh,chiSqFinal,subscriptMap)
    CALL closeNcdfFile(tempFileID,msg1)
    CALL closeNcdfFile(stage2OutID,msg2)
    ! Release all allocated memory
    IF (ALLOCATED(cldRadiance)) DEALLOCATE(cldRadiance)
    IF (ALLOCATED(radWaveNumbers)) DEALLOCATE(radWaveNumbers)
    IF (ALLOCATED(startIndices)) DEALLOCATE(startIndices)
    IF (ALLOCATED(numMembers)) DEALLOCATE(numMembers)
    IF (ALLOCATED(cldyMemberList)) DEALLOCATE(cldyMemberList)
    IF (ALLOCATED(memberList)) DEALLOCATE(memberList)
    IF (ALLOCATED(xOffset)) DEALLOCATE(xOffset)
    IF (ALLOCATED(xGain)) DEALLOCATE(xGain)
    IF (ALLOCATED(pSfc)) DEALLOCATE(pSfc)
    IF (ALLOCATED(tSfc)) DEALLOCATE(tSfc)
    IF (ALLOCATED(kchan)) DEALLOCATE(kchan)
    IF (ALLOCATED(linWaveNumbers)) DEALLOCATE(linWaveNumbers)
    IF (ALLOCATED(chiSqFinal)) DEALLOCATE(chiSqFinal)
    IF (ALLOCATED(pressure)) DEALLOCATE(pressure)
    IF (ALLOCATED(subscriptMap)) DEALLOCATE(subscriptMap)
    IF (ALLOCATED(clusteridGrid)) DEALLOCATE(clusteridGrid)
    IF (ALLOCATED(xGroupGrid)) DEALLOCATE(xGroupGrid)
    IF (ALLOCATED(xGroup)) DEALLOCATE(xGroup)
    IF (ALLOCATED(clusRadGroup)) DEALLOCATE(clusRadGroup)
    IF (ALLOCATED(cldyLats)) DEALLOCATE(cldyLats)
    IF (ALLOCATED(cldyLons)) DEALLOCATE(cldyLons)
    IF (ALLOCATED(cldySatZens)) DEALLOCATE(cldySatZens)
    IF (ALLOCATED(origLats)) DEALLOCATE(origLats)
    IF (ALLOCATED(origLons)) DEALLOCATE(origLons)
    IF (ALLOCATED(origSatZen)) DEALLOCATE(origSatZen)
    IF (ALLOCATED(origGeog)) DEALLOCATE(origGeog)


CONTAINS


    ! -----------------------------------------------------------------
    ! Use the cldyMemberList values to index into the actual lat, lon,
    ! and SATZEN value in the SDR data in rad file. Store the values to
    ! the corresponding grid arrays. Use this subroutine if the objective
    ! is to map only the lats, lons and satzens of cloudy pixels.
    ! -----------------------------------------------------------------
    SUBROUTINE mapLatLon(nRec,nElements,lats,lons,satZens,geography, &
               memberList,cldyMemberList,startIndices,numMembers, &
               latGrid,lonGrid,satZenGrid,geogGrid)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   mapLatLon
!
! PURPOSE:
!
!   Transfer geographic information from cluster indexed storage to
!   line/element grid storage.
!
! SYNTAX:
!
!   CALL mapLatLon(nRec, nElements, lats, lons, satZens, geography,
!      memberList, cldyMemberList, startIndices, numMembers, latGrid,
!      lonGrid, satZenGrid, geogGrid)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   nRec            INTEGER  Number of records
!   nElements       INTEGER  Number of elements
!   lats            REAL     Array of latitudes
!   lons            REAL     Array of longitudes
!   satZens         REAL     Array of satellite zenith angles
!   geography       INTEGER  Geography IDs
!   memberList      INTEGER  List of cluster members
!   cldyMemberList  INTEGER  Pointer (0-based) to cloudy-pixel file
!                            location of each pixel in cluster
!   startIndices    INTEGER  Index into member lists for first
!                            member of each cluster (0-based)
!   numMembers      INTEGER  Numbers of members per group
!
!   INPUTS/OUTPUTS:
!
!   latGrid         REAL     Grid of latitudes
!   lonGrid         REAL     Grid of longitudes
!   satZenGrid      REAL     Grid of satellite-zenith angles
!   geogGrid        INTEGER  Grid of geography IDs
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    ! The parameters
        INTEGER,                                     INTENT(IN)    :: nRec,nElements
        REAL,                         DIMENSION(:),  INTENT(IN)    :: lats,lons
        REAL,                         DIMENSION(:),  INTENT(IN)    :: satZens
        INTEGER(SELECTED_INT_KIND(4)),DIMENSION(:),  INTENT(IN)    :: geography
        INTEGER,                      DIMENSION(:),  INTENT(IN)    :: memberList
        INTEGER,                      DIMENSION(:),  INTENT(IN)    :: cldyMemberList
        INTEGER,                      DIMENSION(:),  INTENT(IN)    :: startIndices
        INTEGER,                      DIMENSION(:),  INTENT(IN)    :: numMembers
        REAL,                         DIMENSION(:,:),INTENT(INOUT) :: latGrid
        REAL,                         DIMENSION(:,:),INTENT(INOUT) :: lonGrid
        REAL,                         DIMENSION(:,:),INTENT(INOUT) :: satZenGrid
        INTEGER(SELECTED_INT_KIND(4)),DIMENSION(:,:),INTENT(INOUT) :: geogGrid
    ! Local variables
        INTEGER     :: i,j,indx,startIndex,line,element
        REAL        :: theLat, theLon

        DO i = 1, nRec
           startIndex = startIndices(i)
           DO j = 1,numMembers(i)
              indx = memberList(startIndex + j)
              line = indx / nElements
              element = indx - line * nElements
              theLat = lats(indx)
              theLon = lons(indx)
              latGrid(line,element) = theLat
              lonGrid(line,element) = theLon
              IF (theLat .GT. 99.0) THEN
                 latGrid(line,element) = MISSING_REAL
              ENDIF
              IF (theLon .GT. 180.0)  THEN
                 lonGrid(line,element) = MISSING_REAL
              ENDIF
              satZenGrid(line,element) = satZens(indx)
              geogGrid(line,element) = geography(indx)
           ENDDO
        ENDDO
      END SUBROUTINE mapLatLon

    ! -----------------------------------------------------------------------
    ! Returns a subscriptMap of size nChON if the wavenumbers of linear
    ! coeffts data and radiance wavenumbers match within tolerance waveThresh.
    ! If an ON channel in kchan does not have a matching channel in radiance
    ! data then we have a fatal error.
    ! -----------------------------------------------------------------------
    SUBROUTINE createSubscriptMap(kchan,nchan,radChannels,radWaveNumbers, &
            linWaveNumbers,waveThresh,subscriptMap)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   createSubscriptMap
!
! PURPOSE:
!
!   Create an index list to map from radiance data file channels to
!   channels used with the linear coefficients.
!
! SYNTAX:
!
!   CALL createSubscriptMap(kchan, nchan, radChannels,
!      radWaveNumbers, linWaveNumbers, waveThresh, subscriptMap)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   kchan           LOGICAL  Channel on/off mask
!   nchan           INTEGER  Number of channels
!   radChannels     INTEGER  Number of channels from rad file
!   radWaveNumbers  REAL     Center wavenumbers of channels from
!                            rad file
!   linWaveNumbers  REAL     Center wavenumbers of channels from
!                            linear coefficients file
!   waveThresh      REAL     Threshold channel wavenumber
!                            difference to accept as a match
!
!   INPUTS/OUTPUTS:
!
!   subscriptMap    INTEGER  Index map to rad file channel for each
!                            linear coeffient file channel
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    ! The parameters
        LOGICAL, DIMENSION(:), INTENT(IN)    :: kchan
        INTEGER,               INTENT(IN)    :: nchan, radChannels
        REAL,    DIMENSION(:), INTENT(IN)    :: radWaveNumbers ! from rad data
        REAL,    DIMENSION(:), INTENT(IN)    :: linWaveNumbers ! from coeffts file
        REAL,                  INTENT(IN)    :: waveThresh ! for consistency check
        INTEGER, DIMENSION(:), INTENT(INOUT) :: subscriptMap
    ! Local variables
        INTEGER    :: i, ii, j
        LOGICAL    :: found

        ii = 0
        DO i = 1, nchan
           IF (kchan(i)) THEN
              found = .FALSE.
              DO j = 1, radChannels
                 IF (found) EXIT
                 IF (ABS(linWaveNumbers(i) - radWaveNumbers(j)) .LE. waveThresh) THEN
                    ii = ii + 1
                    subScriptMap(ii) = j
                    found = .TRUE.
                 ENDIF
              ENDDO
              IF (.NOT.(found)) THEN
                 PRINT *, "Err[createSubScriptMap::Matching rad channel ", &
                       "not found for ", linWaveNumbers(i),"]"
                 CALL errorHalt(1)
              ENDIF
           ENDIF
        ENDDO
     END SUBROUTINE createSubScriptMap


    ! ------------------------------------------------------------------
    ! For each cluster, pack the radiance data of all cloudy
    ! cluster members' ON channels. Compute QC vector and set up the
    ! necessary parameters to invoke linRetr(). Finally do the
    ! bounds testing
    ! ------------------------------------------------------------------
    SUBROUTINE drvRetr(nRec,nPar,startIndices,subscriptMap,numMembers, &
          cldyMemberList,memberList,cldRadiance,chiSqFinal,chiSqThresh, &
          qcFactor,xOff,xGain,pSfc,tSfc,vCoordTyp,pressure,ih2o,IGR,NGR,&
          clusRadGroup,xGroup,xGroupGrid,clusteridGrid)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   drvRetr
!
! PURPOSE:
!
!   Driver for linear retrieval, by group. Do the bounds testing.
!
! SYNTAX:
!
!   CALL drvRetr(nRec, nPar, startIndices, subscriptMap, numMembers,
!      cldyMemberList, memberList, cldRadiance, chiSqFinal,
!      chiSqThresh, qcFactor, xOff, xGain, pSfc, tSfc, pressure,
!      ih2o, IGR, NGR, clusRadGroup, xGroup, xGroupGrid,
!      clusteridGrid)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   nRec            INTEGER             Number of records
!   nPar            INTEGER             Total number of elements in
!                                       retrieval vector
!   startIndices    INTEGER             Index into member lists for
!                                       first member of each
!                                       cluster (0-based)
!   subscriptMap    INTEGER             Index map to rad file
!                                       channel for each linear
!                                       coeffient file channel
!   numMembers      INTEGER             Numbers of members per
!                                       group
!   cldyMemberList  INTEGER             Pointer (0-based) to
!                                       cloudy-pixel file location
!                                       of each pixel in cluster
!   memberList      INTEGER             List of cluster members
!   cldRadiance     REAL                Radiance data at cloudy
!                                       pixels
!   chiSqFinal      REAL                Final chi-squared metric
!   chiSqThresh     REAL                Threshold value of
!                                       chi-squared needed for
!                                       convergence
!   qcFactor        REAL                Multiplier of chi-squared
!                                       threshold to accept as
!                                       passing QC
!   xOff            REAL                linear offset coefficients
!   xGain           REAL                linear gain coefficients
!   pSfc            REAL                Surface pressure
!   tSfc            REAL                Surface-level air
!                                       temperature
!   pressure        REAL                Pressure on atmospheric
!                                       levels
!   ih2o            INTEGER
!   IGR             TYPE(STATEINDEX_T)  Starting indices for
!                                       retrieved sections of
!                                       geophysical state vector
!   NGR             TYPE(STATEINDEX_T)  Number of elements for
!                                       retrieved sections of
!                                       geophysical state vector
!
!   INPUTS/OUTPUTS:
!
!   clusRadGroup    REAL                Radiances for the group
!   xGroup          REAL                Retrieval state vector for
!                                       the group
!   xGroupGrid      REAL                Grid of retrieved
!                                       parameters
!   clusteridGrid   INTEGER             Grid of cluster indices
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    ! The parameters
        INTEGER,                   INTENT(IN)    :: nRec,nPar,ih2o
        INTEGER, DIMENSION(:),     INTENT(IN)    :: startIndices
        INTEGER, DIMENSION(:),     INTENT(IN)    :: subscriptMap
        INTEGER, DIMENSION(:),     INTENT(IN)    :: numMembers
        INTEGER, DIMENSION(:),     INTENT(IN)    :: cldyMemberList
        INTEGER, DIMENSION(:),     INTENT(IN)    :: memberList
        REAL,    DIMENSION(:,0:),  INTENT(IN)    :: cldRadiance
        REAL,    DIMENSION(:),     INTENT(IN)    :: chiSqFinal
        REAL,                      INTENT(IN)    :: chiSqThresh
        REAL,                      INTENT(IN)    :: qcFactor
        REAL,    DIMENSION(:,:),   INTENT(IN)    :: xOff
        REAL,    DIMENSION(:,:,:), INTENT(IN)    :: xGain
        REAL,    DIMENSION(:),     INTENT(IN)    :: pSfc,tSfc
        CHARACTER(LEN=mxCoordTyp), INTENT(IN)    :: vCoordTyp
        REAL,    DIMENSION(:),     INTENT(IN)    :: pressure
        TYPE(StateIndex_t),        INTENT(IN)    :: IGR,NGR
        REAL,    DIMENSION(:,:),   INTENT(INOUT) :: clusRadGroup,xGroup
        REAL,    DIMENSION(:,:,:), INTENT(INOUT) :: xGroupGrid
        INTEGER, DIMENSION(:,:),   INTENT(INOUT) :: clusteridGrid
    ! Local variables
        INTEGER     :: i,j,k

        INTEGER, DIMENSION(:), ALLOCATABLE :: qcVector
        ALLOCATE(qcVector(nRec))

        DO i = 1, nRec
           clusRadGroup = 0.0
           CALL packClusterRadData(subscriptMap,numMembers(i), &
                   startIndices(i),cldyMemberList,cldRadiance, &
                   clusRadGroup(:,1:numMembers(i)))
           qcVector(i) = computeQCvalue(chiSqFinal(i),chiSqThresh,qcFactor)
           xGroup = 0.0
           CALL linRetr(xOff(:,i),xGain(:,:,i),qcVector(i),IGR,NGR, &
                vCoordTyp,pressure,pSfc(i),tSfc(i),ih2o,numMembers(i), &
                clusRadGroup(:,1:numMembers(i)),xGroup(1:nPar,1:numMembers(i)))

           ! Do bounds testing for xGroup parameters
           DO j = 1,numMembers(i)
              ! Ensure that the value of particle equivalent diameter
              ! of ice or liquid is atleast 0.0 and not negative.
              IF (NGR%cldLiq >= 4) THEN
                 IF (xGroup(IGR%cldLiq+3, j) .LT. 0.0) &
                     xGroup(IGR%cldLiq+3, j) = 0.0
              ENDIF
              IF (NGR%cldIce >= 4) THEN
                 IF (xGroup(IGR%cldIce+3, j) .LT. 0.0) &
                    xGroup(IGR%cldIce+3, j) = 0.0
              ENDIF
              ! Ensure that the value
              ! of total water (ice or liquid) is atleast 0.0 and not negative.
              IF (NGR%cldLiq >= 3) THEN
                 IF (xGroup(IGR%cldLiq+2, j) .LT. 0.0) &
                    xGroup(IGR%cldLiq+2, j) = 0.0
              ENDIF
              IF (NGR%cldIce >= 3) THEN
                 IF (xGroup(IGR%cldIce+2, j) .LT. 0.0) &
                    xGroup(IGR%cldIce+2, j) = 0.0
              ENDIF
           ENDDO
           CALL saveParameters(xGroup(1:nPar,1:numMembers(i)),memberList, &
                nLines,nElements,startIndices(i),numMembers(i),i, &
                xGroupGrid,clusteridGrid)
        ENDDO

        IF (ALLOCATED(qcVector)) DEALLOCATE(qcVector)

    END SUBROUTINE drvRetr


    ! ------------------------------------------------------------
    ! Store all xGroup values to the nLinesxnElements grid
    ! and the clusterID info to the nLinesxnElements clusteridGrid.
    ! ------------------------------------------------------------

    SUBROUTINE saveParameters(xGroup,memberList,nLines,nElements, &
             startIndex,nMembers,clusterID, &
             xGroupGrid,clusteridGrid)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   saveParameters
!
! PURPOSE:
!
!   Transfer retrieved results from cluster indexed storage to
!   line/element grid storage.
!
! SYNTAX:
!
!   CALL saveParameters(xGroup, memberList, nLines, nElements,
!      startIndex, nMembers, clusterID, xGroupGrid, clusteridGrid)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   xGroup         REAL     Retrieval state vector for the group
!   memberList     INTEGER  List of cluster members
!   nLines         INTEGER  number of lines
!   nElements      INTEGER  Number of elements
!   startIndex     INTEGER  Starting index for cluster members
!   nMembers       INTEGER  Number of members per group
!   clusterID      INTEGER  cluster index
!
!   INPUTS/OUTPUTS:
!
!   xGroupGrid     REAL     Grid of retrieved parameters
!   clusteridGrid  INTEGER  Grid of cluster indices
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    ! The parameters
        REAL,    DIMENSION(:,:),   INTENT(IN)    :: xGroup
        INTEGER, DIMENSION(:),     INTENT(IN)    :: memberList
        INTEGER,                   INTENT(IN)    :: nLines,nElements,nMembers
        REAL,    DIMENSION(:,:,:), INTENT(INOUT) :: xGroupGrid
        INTEGER, DIMENSION(:,:),   INTENT(INOUT) :: clusteridGrid
        INTEGER,                   INTENT(IN)    :: clusterID,startIndex
     ! Local variables
       INTEGER     :: i, j
       INTEGER     :: index,line,element

     ! Adust line and element to 1-based index.
       DO i = 1, nMembers
          index = memberList(startIndex + i)
          line = index / nElements
          element = index - (line * nElements) + 1
          line = line + 1
          clusteridGrid(element,line) = clusterID
          DO j = 1, nPar
             xGroupGrid(element,line,j) = xGroup(j,i)
          ENDDO
       ENDDO
    END SUBROUTINE saveParameters

    ! -------------------------------------------------------------
    ! Packs rad data from all ON channels for all members of
    ! cluster into clusterRadGroup.
    ! -------------------------------------------------------------
    SUBROUTINE packClusterRadData(subscriptMap,nMembers, &
             startIndex,cldyMemberList,cldRadiance, &
             clusterRadGroup)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   packClusterRadData
!
! PURPOSE:
!
!   Pack rad data from all ON channels for all members of cluster into
!   clusterRadGroup.
!
! SYNTAX:
!
!   CALL packClusterRadData(subscriptMap, nMembers, startIndex,
!      cldyMemberList, cldRadiance, clusterRadGroup)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   subscriptMap     INTEGER  Index map to rad file channel for
!                             each linear coeffient file channel
!   nMembers         INTEGER  Number of members per group
!   startIndex       INTEGER  Starting index for cluster members
!   cldyMemberList   INTEGER  Pointer (0-based) to cloudy-pixel
!                             file location of each pixel in
!                             cluster
!   cldRadiance      REAL     Radiance data at cloudy pixels
!
!   INPUTS/OUTPUTS:
!
!   clusterRadGroup  REAL     Radiances for the group
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    ! The parameters
        INTEGER, DIMENSION(:),   INTENT(IN)    :: subscriptMap
        INTEGER,                 INTENT(IN)    :: nMembers,startIndex
        REAL,    DIMENSION(:,0:), INTENT(IN)   :: cldRadiance
        INTEGER, DIMENSION(:),   INTENT(IN)    :: cldyMemberList
        REAL,    DIMENSION(:,:), INTENT(INOUT) :: clusterRadGroup

       ! Local variables
        INTEGER         :: i, j
        INTEGER         :: radIndex
       ! CLDY_MEMBER_LIST is 0-based and assumes
       ! rad array to be 0-based as well. So adjust it to 1-based
        DO i = 1, nMembers
           radIndex = cldyMemberList(startIndex + i)
           DO j = 1, SIZE(subscriptMap)
              clusterRadGroup(j, i) = cldRadiance(subscriptMap(j),radIndex)
           ENDDO
        ENDDO

    END SUBROUTINE packClusterRadData

    ! -------------------------------------------------------------
    ! Compute QC as good (0) if chiSqFinal < chiSqThresh * qcFactor
    ! -------------------------------------------------------------
    FUNCTION computeQCvalue(chiSqFinal,chiSqThresh,qcFactor)

!<f90Function>**********************************************************
!
! NAME:
!
!   computeQCvalue
!
! PURPOSE:
!
!   Compute QC
!
! SYNTAX:
!
!   Results=computeQCvalue(chiSqFinal, chiSqThresh, qcFactor)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   chiSqFinal      REAL     Final chi-squared metric
!   chiSqThresh     REAL     Threshold value of chi-squared needed
!                            for convergence
!   qcFactor        REAL     Multiplier of chi-squared threshold to
!                            accept as passing QC
!
!   * OPTIONAL
!
! RETURN:
!
!     INTEGER
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

       INTEGER  :: computeQCvalue
      ! parameters
       REAL, INTENT(IN)    :: chiSqFinal,chiSqThresh
       REAL, INTENT(IN)    :: qcFactor
      ! local variables

       IF (chiSqFinal .LT. (chiSqThresh * qcFactor)) THEN
         computeQCvalue = 0
       ELSE
             computeQCvalue = 1
       ENDIF
    END FUNCTION computeQCvalue


   ! For Intermediate data verification only
   SUBROUTINE writeToFile(fileID,nMembersTotal,nChannels,nRec,cldyMemberList, &
                radianceData,radWaveNumbers,startIndices, &
                numMembers,linWaveNumbers, &
                xOff,xGain,pSfc,tSfc,nchan,kchan,nPar,nGroup,nChOn, &
                chiSqThres,chiSqFinal,subscriptMap)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   writeToFile
!
! PURPOSE:
!
!   For Intermediate data verification only
!
! SYNTAX:
!
!   CALL writeToFile(fileID, nMembersTotal, nChannels, nRec,
!      cldyMemberList, radianceData, radWaveNumbers, startIndices,
!      numMembers, linWaveNumbers, xOff, xGain, pSfc, tSfc, nchan,
!      kchan, nPar, nGroup, nChOn, chiSqThres, chiSqFinal,
!      subscriptMap)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   fileID          INTEGER  Index of file
!   nMembersTotal   INTEGER  Number of members of all clusters
!                            combined
!   nChannels       INTEGER  Number of channels
!   nRec            INTEGER  Number of records
!   cldyMemberList  INTEGER  Pointer (0-based) to cloudy-pixel file
!                            location of each pixel in cluster
!   radianceData    REAL     Radiance data
!   radWaveNumbers  REAL     Center wavenumbers of channels from
!                            rad file
!   startIndices    INTEGER  Index into member lists for first
!                            member of each cluster (0-based)
!   numMembers      INTEGER  Numbers of members per group
!   linWaveNumbers  REAL     Center wavenumbers of channels from
!                            linear coefficients file
!   xOff            REAL     linear offset coefficients
!   xGain           REAL     linear gain coefficients
!   pSfc            REAL     Surface pressure
!   tSfc            REAL     Surface-level air temperature
!   nchan           INTEGER  Number of channels
!   kchan           LOGICAL  Channel on/off mask
!   nPar            INTEGER  Total number of elements in retrieval
!                            vector
!   nGroup          INTEGER  TBD
!   nChOn           INTEGER  Number of channels turned on
!   chiSqThres      REAL
!   chiSqFinal      REAL     Final chi-squared metric
!   subscriptMap    INTEGER  Index map to rad file channel for each
!                            linear coeffient file channel
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

   ! The parameters
     INTEGER,                   INTENT(IN) :: fileID
     INTEGER,                   INTENT(IN) :: nMembersTotal, nChannels, nRec
     INTEGER,                   INTENT(IN) :: nchan
     INTEGER, DIMENSION(:),     INTENT(IN) :: cldyMemberList
     INTEGER, DIMENSION(:),     INTENT(IN) :: startIndices, numMembers
     LOGICAL, DIMENSION(:),     INTENT(IN) :: kchan
     REAL,    DIMENSION(:,:),   INTENT(IN) :: xOff, radianceData
     REAL,    DIMENSION(:,:,:), INTENT(IN) :: xGain
     REAL,    DIMENSION(:),     INTENT(IN) :: pSfc,tSfc,radWaveNumbers
     REAL,    DIMENSION(:),     INTENT(IN) :: linWaveNumbers
     INTEGER,                   INTENT(IN) :: nPar,nGroup,nChOn
     REAL,                      INTENT(IN) :: chiSqThres
     REAL,    DIMENSION(:),     INTENT(IN) :: chiSqFinal
     INTEGER, DIMENSION(:),     INTENT(IN) :: subscriptMap

     ! Local variables
     CALL writeNcdfDim(fileID, dimName='nMembersTotal',dimLen=nMembersTotal)
     CALL writeNcdfDim(fileID, dimName='nChannels',dimLen=nChannels)
     CALL writeNcdfDim(fileID, dimName='nRec', dimLen=nRec)
     CALL writeNcdfDim(fileID, dimName='nPar', dimLen=nPar)
     CALL writeNcdfDim(fileID, dimName='nChOn', dimLen=nChOn)
     CALL writeNcdfDim(fileID, dimName='nGroup', dimLen=nGroup)
     CALL writeNcdfData(fileID, cldyMemberList, varName='CLDY_MEMBER_LIST', &
             varLenName = ('nMembersTotal'), &
             varLongName = 'Cluster membership data', &
             varUnit = 'unitless')
     CALL writeNcdfData(fileID, radianceData, varName='rad', &
             varLenName=(/'nChannels    ', 'nMembersTotal'/), &
             varLongName='Radiance data', &
             varUnit = "mW / (m**2 cm**-1 ster)" )
     CALL writeNcdfData(fileID, startIndices, varName='START_INDEX', &
             varLenName = ('nRec'), &
             varLongName ='Cluster starting index', &
             varUnit='unitless')
     CALL writeNcdfData(fileID, numMembers, varName='N_MEMBERS', &
             varLenName = ('nRec'), &
             varLongName ='#Cluster members', varUnit='unitless')
     CALL writeNcdfData(fileID, xOff, varName='xOff', &
             varLenName = (/'nPar  ', 'nGroup'/), &
             varLongName = 'Linear inversion offset', &
             varUnit = 'unitless')
     CALL writeNcdfData(fileID, xGain, varName='xGain', &
             varLenName = (/'nPar  ', 'nChOn ', 'nGroup'/), &
             varLongName = 'Linear inversion gain',  &
             varUnit = 'unitless')
     CALL writeNcdfData(fileID, pSfc, varName='pSfc', &
             varLenName = ('nGroup'), &
             varLongName='Surface pressure', varUnit='mb')
     CALL writeNcdfData(fileID, tSfc, varName='tSfc', &
             varLenName = ('nGroup'), &
             varLongName='Surface temperature', varUnit='K')
     CALL writeNcdfAttr(fileID, attrName='radWavenumbers',attr=radWaveNumbers)
     CALL writeNcdfAttr(fileID, attrName='nchan', attr=nchan)
     CALL writeNcdfAttr(fileID, attrName='kchan', attr=kchan)
     CALL writeNcdfAttr(fileID, attrName='chiSqThres', attr=chiSqThres)
     CALL writeNcdfAttr(fileID, attrName='linWaveNumbers',attr=linWaveNumbers)
     CALL writeNcdfData(fileID, chiSqFinal, varName='chiSqFinal', &
             varLenName=('nGroup'), varUnit='unitless')
     CALL writeNcdfData(fileID, subscriptMap, varName='SubscriptMap', &
             varLenName = ('nChOn'), &
             varLongName='Subscripts of ON channels', varUnit='unitless')
  END SUBROUTINE writeToFile

  ! -------------------------------------------------------------------------
  ! Reads the control parameters to execute the program.
  ! -------------------------------------------------------------------------
  SUBROUTINE readConfigParameters(radFileName,clusFileName,tempFileName, &
                linCoefFileName,originalFileName,stage2OutFile, &
                QCfactor,waveThreshold)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   readConfigParameters
!
! PURPOSE:
!
!   Reads the control parameters to execute the program.
!
! SYNTAX:
!
!   CALL readConfigParameters(radFileName, clusFileName,
!      tempFileName, linCoefFileName, originalFileName,
!      stage2OutFile, QCfactor, waveThreshold)
!
! ARGUMENTS:
!
!   INPUTS/OUTPUTS:
!
!   radFileName       CHAR  File path for radiance file
!   clusFileName      CHAR  File path for cluster data file
!   tempFileName      CHAR  File path for debug file
!   linCoefFileName   CHAR  File path for linear coefficients data
!   originalFileName  CHAR  File path for radiance file in original
!                           gridded form
!   stage2OutFile     CHAR  File path for output of linear
!                           retrieval
!   QCfactor          REAL  Multiplier of chi-squared threshold to
!                           accept as passing QC
!   waveThreshold     REAL  Threshold channel wavenumber difference
!                           to accept as a match
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

   ! The parameters
      CHARACTER(LEN=200), INTENT(INOUT) :: radFileName,clusFileName
      CHARACTER(LEN=200), INTENT(INOUT) :: linCoefFileName,originalFileName
      CHARACTER(LEN=200), INTENT(INOUT) :: tempFileName,stage2OutFile
      REAL,               INTENT(INOUT) :: QCfactor,waveThreshold

      NAMELIST /Stage2Control/ QCfactor,waveThreshold
      NAMELIST /RadFiles/ radFileName
      NAMELIST /ClusterFiles/ clusFileName
      NAMELIST /CoefFiles/ linCoefFileName
      NAMELIST /OrigFiles/ originalFileName
      NAMELIST /Stage2OutFiles/ stage2OutFile
      NAMELIST /TempFiles/ tempFileName
      READ (*,Stage2Control)
      READ (*,RadFiles)
      READ (*,ClusterFiles)
      READ (*,CoefFiles)
      READ (*,OrigFiles)
      READ (*,Stage2OutFiles)
      READ (*,TempFiles)
  END SUBROUTINE readConfigParameters

END PROGRAM LinRetrGroup

