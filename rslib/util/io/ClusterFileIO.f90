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
! MODULE ClusterFileIO
!
! This module contains routines for loading cluster file dimensions, attributes
!   and variables, thereby compartmentalizing the data input.  After calling 
!   'getDimCluster' and allocating necessary arrays, the user should call 
!   'getClusterData' to load all data prior to looping over profiles.  
!
!----------------------------------------------------------------------------
MODULE ClusterFileIO

! <f90Module>***********************************************************
!
! NAME:
!
!   ClusterFileIO
!
! PURPOSE:
!
!   Routines for loading cluster file dimensions, attributes and variables, 
!   thereby compartmentalizing the data input.
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

  
  USE ncdf_module
  
  IMPLICIT none
  
  PRIVATE
  
  PUBLIC :: getDimCluster,getClusterData
  
  CONTAINS

    ! -------------------------------------------------------------------------
    ! Reads and returns the dimensions and attributes that serve as actual
    ! dimensions for the variable data in cluster file.
    ! -------------------------------------------------------------------------
    SUBROUTINE getDimCluster(ncid,fname,nChannels,nClust,nMembersTotal,nLines, &
                      nElements)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   getDimCluster
!
! PURPOSE:
!
!   Reads and returns the dimensions and attributes for the variable data in 
!   cluster file
!
! SYNTAX:
!
!   CALL getDimCluster(ncid, fname, nChannels, nClust, nMembersTotal, 
!      nLines, nElements)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   fname          CHAR     File name
!   
!   INPUTS/OUTPUTS:
!   
!   ncid           INTEGER  ID for ncdf file 
!   nChannels      INTEGER  Number of channels
!   nClust         INTEGER  Number of clusters
!   nMembersTotal  INTEGER  Number of members of all clusters combined
!   nLines*        INTEGER  number of lines
!   nElements*     INTEGER  Number of elements
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !-- I/O variables:
      INTEGER,          INTENT(INOUT)           :: ncid
      CHARACTER(LEN=*), INTENT(IN)              :: fname
      INTEGER,          INTENT(INOUT)           :: nChannels ! #channels
      INTEGER,          INTENT(INOUT)           :: nClust      ! #clusters
      INTEGER,          INTENT(INOUT)           :: nMembersTotal !#cloudy pixels
      INTEGER,          INTENT(INOUT), OPTIONAL :: nLines, nElements
      
      CALL openNcdfFile(ncid, fname, status='old')
      nChannels = readNcdfDim(ncid,'nChannels')
      nClust  = readNcdfDim(ncid,'nRec')
      nMembersTotal = readNcdfDim(ncid,'nMembersTotal')
      IF (PRESENT(nLines)) &
         CALL readNcdfAttr(ncid,attrName='nLines',attr=nLines)
      IF (PRESENT(nElements)) &
         CALL readNcdfAttr(ncid,attrName='nElements',attr=nElements)
      CALL closeNcdfFile(ncid)

    END SUBROUTINE getDimCluster
    
    ! -------------------------------------------------------------------------
    ! Reads and returns all cluster information necessary for stage 2 retrieval
    ! from the cluster output file specified in fname.
    ! -------------------------------------------------------------------------
    SUBROUTINE getClusterData(ncid,fname,cloudyMembers,startIndices,nMembers, &
                         memberList)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   getClusterData
!
! PURPOSE:
!
!   Reads and returns all cluster information necessary for stage 2 retrieval.
!
! SYNTAX:
!
!   CALL getClusterData(ncid, fname, cloudyMembers, startIndices, 
!      nMembers, memberList)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   fname          CHAR     File name
!   
!   INPUTS/OUTPUTS:
!   
!   ncid           INTEGER  ID for ncdf file 
!   cloudyMembers  INTEGER  Pointer (0-based) to cloudy-pixel file location 
!                           of each pixel in cluster 
!   startIndices   INTEGER  Index into member lists for first member of each 
!                           cluster (0-based) 
!   nMembers       INTEGER  Number of members per group
!   memberList     INTEGER  List of cluster members
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !-- I/O variables
      INTEGER,               INTENT(INOUT) :: ncid
      CHARACTER(LEN=*),      INTENT(IN)    :: fname
      INTEGER, DIMENSION(:), INTENT(INOUT) :: cloudyMembers ! cloudy member list
      INTEGER, DIMENSION(:), INTENT(INOUT) :: startIndices 
      INTEGER, DIMENSION(:), INTENT(INOUT) :: nMembers  ! #members per cluster
      INTEGER, DIMENSION(:), INTENT(INOUT) :: memberList
  
      CALL openNcdfFile(ncid,fname,status='old')
      CALL readNcdfData(ncid,cloudyMembers,'CLDY_MEMBER_LIST')
      CALL readNcdfData(ncid,startIndices,'START_INDEX')
      CALL readNcdfData(ncid,nMembers,'N_MEMBERS')
      CALL readNcdfData(ncid,memberList,'MEMBER_LIST')
      CALL closeNcdfFile(ncid)
    END SUBROUTINE getClusterData

END MODULE ClusterFileIO
