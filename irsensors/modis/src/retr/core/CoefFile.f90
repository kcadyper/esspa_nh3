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

MODULE CoefFile

! <f90Module>***********************************************************
!
! NAME:
!
!   CoefFile
!
! PURPOSE:
!
!   Manipulate the linear coefficient file.
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

   USE ncdf_module, ONLY: openNcdfFile,closeNcdfFile,readNcdfDim, &
         readNcdfData,readNcdfAttr
   USE StateIndexModule, ONLY: maxMol,StateIndex_t

   IMPLICIT NONE

   PRIVATE

   PUBLIC:: getCoefDimensions,getLinearCoefftsData

   CONTAINS

      ! -----------------------------------------------------------------------
      ! This subroutine reads the dimensions and attributes that serve as
      ! actual dimensions for the various data in linear coefficients file.
      ! -----------------------------------------------------------------------
      SUBROUTINE getCoefDimensions(ncid,fileName,nPar,nGroup,nChOn,nchan,nlev)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   getCoefDimensions
!
! PURPOSE:
!
!   Read dimensions and attributes that serve as actual dimensions for
!   the various data in linear coefficients file.
!
! SYNTAX:
!
!   CALL getCoefDimensions(ncid, fileName, nPar, nGroup, nChOn,
!      nchan, nlev)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   fileName  CHAR     File name
!
!   INPUTS/OUTPUTS:
!
!   ncid      INTEGER  ID for ncdf file
!   nPar      INTEGER  Total number of elements in retrieval vector
!   nGroup    INTEGER  TBD
!   nChOn     INTEGER  Number of channels turned on
!   nchan     INTEGER  Number of channels
!   nlev      INTEGER  Number of atmospheric levels
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

      ! The parameters
         INTEGER,          INTENT(INOUT) :: ncid
         INTEGER,          INTENT(INOUT) :: nPar,nGroup,nChOn
         INTEGER,          INTENT(INOUT) :: nchan,nlev ! attributes
         CHARACTER (LEN=*),INTENT(IN)    :: fileName

         CALL openNcdfFile(ncid,fileName,status='old')
         nPar = readNcdfDim(ncid, 'nPar')
         nGroup = readNcdfDim(ncid, 'nGroup')
         nChOn = readNcdfDim(ncid, 'nChOn')
         CALL readNcdfAttr(ncid, attrName='nchan', attr=nchan)
         CALL readNcdfAttr(ncid, attrName='nlev', attr=nlev)
         CALL closeNcdfFile(ncid)
      END SUBROUTINE getCoefDimensions

      ! ---------------------------------------------------------------------
      ! This subroutine loads the various linear coefficients data from the
      ! coefficients file, fileName and builds the IGR and NGR structures
      ! ---------------------------------------------------------------------
      SUBROUTINE getLinearCoefftsData(ncid,fileName,xOff,xGain,pSfc,tSfc, &
                     chiSqFinal,kchan,chiSqThres,pressure,linWaveNumbers, &
                     IGR,NGR)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   getLinearCoefftsData
!
! PURPOSE:
!
!   Load various linear coefficients data from the coefficients file
!   and build the IGR and NGR structures
!
! SYNTAX:
!
!   CALL getLinearCoefftsData(ncid, fileName, xOff, xGain, pSfc,
!      tSfc, chiSqFinal, kchan, chiSqThres, pressure, linWaveNumbers,
!      IGR, NGR)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   fileName        CHAR                File name
!
!   INPUTS/OUTPUTS:
!
!   ncid            INTEGER             ID for ncdf file
!   xOff            REAL                linear offset coefficients
!   xGain           REAL                linear gain coefficients
!   pSfc            REAL                Surface pressure
!   tSfc            REAL                Surface-level air
!                                       temperature
!   chiSqFinal      REAL                Final chi-squared metric
!   kchan           logical             Channel on/off mask
!   chiSqThres      REAL
!   pressure        REAL                Pressure on atmospheric
!                                       levels
!   linWaveNumbers  REAL                Center wavenumbers of
!                                       channels from linear
!                                       coefficients file
!   IGR             TYPE(STATEINDEX_T)  Starting indices for
!                                       retrieved sections of
!                                       geophysical state vector
!   NGR             TYPE(STATEINDEX_T)  Number of elements for
!                                       retrieved sections of
!                                       geophysical state vector
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

      ! The parameters
         INTEGER,                   INTENT(INOUT) :: ncid
         CHARACTER(LEN=*),          INTENT(IN)    :: fileName
         REAL,    DIMENSION(:,:),   INTENT(INOUT) :: xOff
         REAL,    DIMENSION(:,:,:), INTENT(INOUT) :: xGain
         REAL,    DIMENSION(:),     INTENT(INOUT) :: pSfc,tSfc,chiSqFinal
         LOGICAL, DIMENSION(:),     INTENT(INOUT) :: kchan
         REAL,                      INTENT(INOUT) :: chiSqThres
         REAL,    DIMENSION(:),     INTENT(INOUT) :: linWaveNumbers
         TYPE(StateIndex_t),        INTENT(INOUT) :: IGR,NGR
         REAL,    DIMENSION(:),     INTENT(INOUT) :: pressure
      ! Local Variables
         INTEGER, DIMENSION(2)  :: temperature,tskin,sfcPress,liqCld, &
                                        iceCld,sfcWnd
         INTEGER, DIMENSION(:), ALLOCATABLE :: mol, molLen

         CALL openNcdfFile(ncid, fileName, status='old')
         ALLOCATE(mol(maxMol))
         ALLOCATE(molLen(maxMol))
         CALL readNcdfAttr(ncid, attrName='Temperature', attr=temperature)
         CALL readNcdfAttr(ncid, attrName='Tskin', attr=tskin)
         CALL readNcdfAttr(ncid, attrName='SurfacePressure', attr=sfcPress)
         CALL readNcdfAttr(ncid, attrName='LiqCloud', attr=liqCld)
         CALL readNcdfAttr(ncid, attrName='IceCloud', attr=iceCld)
         CALL readNcdfAttr(ncid, attrName='SurfaceWinds', attr=sfcWnd)
         CALL readNcdfAttr(ncid, attrName='kchan', attr=kchan)
         CALL readNcdfAttr(ncid, attrName='standardPressureGrid', attr=pressure)
         CALL readNcdfAttr(ncid, attrName='Mol', attr=mol)
         CALL readNcdfAttr(ncid, attrName='MolLen', attr=molLen)
         CALL readNcdfAttr(ncid, attrName='chiSqThres', attr=chiSqThres)
         CALL readNcdfAttr(ncid, attrName='wavenumbers', attr=linWaveNumbers)
         ! Build IGR and NGR structures
         IGR%temp    = temperature(1)
         NGR%temp    = temperature(2)
         IGR%tskin   = tskin(1)
         NGR%tskin   = tskin(2)
         IGR%psfc    = sfcPress(1)
         NGR%psfc    = sfcPress(2)
         IGR%cldliq  = liqCld(1)
         NGR%cldliq  = liqCld(2)
         IGR%cldice  = iceCld(1)
         NGR%cldice  = iceCld(2)
         IGR%wind    = sfcWnd(1)
         NGR%wind    = sfcWnd(2)
         IGR%mol     = mol
         NGR%mol     = molLen

         ! Read all data
         ! Read Tsfc only if temp data is present.
         CALL readNcdfData(ncid, xOff, 'xOff')
         CALL readNcdfData(ncid, xGain, 'xGain')
         CALL readNcdfData(ncid, pSfc, 'pSfc')
         CALL readNcdfData(ncid, chiSqFinal, 'chiSqFinal')

         IF (NGR%tskin > 0)  CALL readNcdfData(ncid, tSfc, 'Tsfc')

         IF (ALLOCATED(mol)) DEALLOCATE(mol)
         IF (ALLOCATED(molLen)) DEALLOCATE(molLen)
         CALL closeNcdfFile(ncid)
     END SUBROUTINE getLinearCoefftsData

END MODULE CoefFile
