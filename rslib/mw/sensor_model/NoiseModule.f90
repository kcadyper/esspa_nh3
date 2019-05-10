!*************************************************************
!   Copyright: AER, Inc., All Right Reserved

!--------------------------------------------------------------
!
!  MODULE Noise: contains subroutines needed for the
!                noise generation. AER Inc. 2004.
!--------------------------------------------------------------
MODULE NoiseModule
  USE constants, ONLY : MISSING_REAL
  USE OSSMWmoduleSubs, ONLY : getOSSselMW,OSS_LEN_ID,OSS_NULL_ID
  USE SensorData, ONLY : setDynNoise,queryDynNoise,dynNoise,queryBeamMap,&
       querySpillover, destroySensorData, LEN_ID_SensorData
  IMPLICIT NONE
  PRIVATE
  !------------------------------------------------------------------------
  !	Public items made available to calling program(s)
  !------------------------------------------------------------------------
  PUBLIC :: addDeviceNoise,deviceNoise,loadNedtStatic,loadNAF,destroyNoise
  !------------------------------------------------------------------------
  !     Declarations of private data (private to the module)
  !------------------------------------------------------------------------
  INTEGER, PARAMETER :: mxscans=365,mxcell=7,mxchan=80
  integer, parameter :: mxLength=256
  INTEGER                                              :: ncellNAF,nscanNAF
  INTEGER                                              :: nchanNedt,nchanNAF
  REAL,                  DIMENSION(mxchan)             :: NedtStatic
  REAL,                  DIMENSION(:,:,:), ALLOCATABLE :: nrfBuf
  CHARACTER(OSS_LEN_ID), DIMENSION(:),     ALLOCATABLE :: NAFchanID

  ! load*Done implementation provides the option to call load*() in advance
  logical                            :: loadNedtStaticDone=.false.
  logical                            :: loadNAFdone=.false.

CONTAINS


  SUBROUTINE addDeviceNoise(Y,devNoise,nchan)
    !---In/Out variables
    INTEGER               :: nchan
    REAL,    DIMENSION(:) :: Y,devNoise
    !---Local variables
    INTEGER, SAVE         :: iseed=1287365
    INTEGER               :: i
    REAL                  :: stdt,gauss

    ! -- Loop over all nchan regardless of kchan flags so that the sequence
    ! -- of gauss calls remains the same no matter the active channel set.
    DO i=1,nchan
       stdt=devNoise(i)
       Y(i)=Y(i)+gauss(iseed,stdt,0.0)
    END DO
    RETURN
  END SUBROUTINE addDeviceNoise


  SUBROUTINE deviceNoise(icellres,tb,rerr,sum_rerr,kchan,nchan,&
       DynamicNoise,spilloverNoise,iscanPosIn,calibAmp,Fnedt,Fnrf)

    USE ReadStdInputs, ONLY: U_osscoefs,F_osscoefs,U_nedt,F_nedt, &
         U_beam,F_beam,U_spov,F_spov
    USE constants, ONLY : CosmBckg

    !---In/Out variables
    LOGICAL               :: DynamicNoise
    LOGICAL, OPTIONAL     :: spilloverNoise
    LOGICAL, OPTIONAL     :: calibAmp
    INTEGER               :: icellres,nchan
    INTEGER, OPTIONAL     :: iscanPosIn
    LOGICAL, DIMENSION(:) :: kchan
    REAL,    DIMENSION(:) :: tb,rerr
    REAL                  :: sum_rerr
    character(len=*), intent(in), optional :: Fnedt
    character(len=*), intent(in), optional :: Fnrf
    !---Local variables
    INTEGER            :: icaset     ! 1 => use the tabulated calibration amplification factors
    REAL               :: tAntn,avefac
    INTEGER            :: nfr,k,i,j,ispec,nch,nChanAll,nchanOSS,iscanPos=1
    INTEGER            :: myicellres
    INTEGER, SAVE      :: init=1,nscan,ncell
    INTEGER, SAVE, DIMENSION(mxchan)                 :: mapnedt
    REAL,    SAVE, DIMENSION(mxchan,mxcell,mxscans)  :: nrf
    CHARACTER(OSS_LEN_ID), DIMENSION(mxchan)         :: ossChanID
    CHARACTER(LEN_ID_SensorData), DIMENSION(:), POINTER     :: tChanID
    CHARACTER(LEN_ID_SensorData), DIMENSION(:), POINTER     :: chanIDc
    CHARACTER(LEN_ID_SensorData), DIMENSION(:), POINTER     :: beamIDc,beamIDav
    LOGICAL                                          :: mySpilloverNoise
    LOGICAL                                          :: myCalibAmp
    INTEGER                                          :: nnoise
    INTEGER, DIMENSION(:),               POINTER     :: iBeam
    REAL,    DIMENSION(:),               POINTER     :: spoev
    REAL,    SAVE, DIMENSION(mxchan)                 :: spoch
    character (len=mxLength)                         :: nedtFile
    character (len=mxLength)                         :: nrfFile

    myicellres = 1
    IF (icellres > 0) myicellres = icellres
    mySpilloverNoise=.FALSE.
    IF (present(spilloverNoise)) mySpilloverNoise=spilloverNoise
    myCalibAmp=.FALSE.
    IF (present(calibAmp)) myCalibAmp=calibAmp
    icaset = 0
    IF (myCalibAmp) icaset=1
    iscanPos=1
    IF (present(iscanPosIn))  iscanPos=iscanPosIn
    !---Consistency checks
    IF (nchan      .GT. mxchan) THEN
       PRINT *,'err[NoiseModule::deviceNoise]:  nchan > mxchan'
       call errorHalt(1)
    END IF
    IF (SIZE(tb)   .LT. nchan) THEN
       PRINT *,'err[NoiseModule::deviceNoise]:  size(tb) and nchan inconsistent'
       call errorHalt(1)
    END IF
    IF (SIZE(kchan).LT. nchan) THEN
       PRINT *,'err[NoiseModule::deviceNoise]:  size(kchan) and nchan inconsistent'
       call errorHalt(1)
    END IF
    IF (SIZE(rerr) .LT. nchan) THEN
       PRINT *,'err[NoiseModule::deviceNoise]:  size(rerr) and nchan inconsistent'
       call errorHalt(1)
    END IF
    !---One-time reading of the noise parameters
    IF(init .EQ. 1)THEN
       if (present(Fnedt)) then
          nedtFile = trim(adjustl(Fnedt))
       else
          nedtFile = trim(adjustl(F_nedt))
       end if
       call getOSSselMW(F_osscoefs,nchanOSS=nchanOSS, &
          chanID=ossChanID(1:nchan))
       IF (nchanOSS /= nchan) THEN
          PRINT *,'err[NoiseModule::deviceNoise]:  '
          PRINT *,'nchanOSS and nchan inconsistent: ',nchanOSS,nchan
          call errorHalt(1)
       END IF
       IF (DynamicNoise) THEN
          !--- obtain identifiers of channels
          IF (OSS_LEN_ID /= LEN_ID_SensorData) THEN
             PRINT *,'err[NoiseModule::deviceNoise]:  '
             PRINT *,'OSS_LEN_ID and LEN_ID_SensorData inconsistent: ', &
                  OSS_LEN_ID,LEN_ID_SensorData
             call errorHalt(1)
          END IF

          !--- obtain indicies that map channels to NEDT parameter records
          CALL setDynNoise(U_nedt,nedtFile)
          CALL queryDynNoise(nnoise,tChanID)
          CALL getMapNedt(nnoise,tChanID,ossChanID,nchan,mapnedt)
          IF (mySpilloverNoise) THEN
             !--- obtain map from channel to beam
             call queryBeamMap(U_beam,F_beam,chanIDc,beamIDc,iBeam)
             !--- obtain spillover per channel, mapping to full channel set
             call querySpillover(U_spov,F_spov,beamIDav,spoev)
             IF (any(chanIDc /= tChanID)) THEN  ! ensure mapnedt<=size(iBeam)
                PRINT *,'err[NoiseModule::deviceNoise]:  ', &
                  'Channel ID inconsistency'
                call errorHalt(1)
             END IF
             IF (any(beamIDc /= beamIDav)) THEN
                PRINT *,'err[NoiseModule::deviceNoise]:  ', &
                  'Beam ID inconsistency'
                call errorHalt(1)
             END IF
             spoch(1:nchan)=spoev(iBeam(mapnedt(1:nchan)))
          ELSE
             spoch=0.
          ENDIF
          IF (associated(tchanID))  deallocate(tchanID)
          IF (associated(chanIDc))  deallocate(chanIDc)
          IF (associated(beamIDc))  deallocate(beamIDc)
          IF (associated(beamIDav)) deallocate(beamIDav)
          IF (associated(ibeam))    deallocate(ibeam)
          IF (associated(spoev))    deallocate(spoev)
       ENDIF
       IF (.not.(DynamicNoise)) THEN
          !--- get the static nedts when they do not depend on antenna temp
          call loadNedtStatic(nedtFile)
          IF (nchanNedt .NE. nchan) THEN
             PRINT *,'err[NoiseModule::deviceNoise]: ',&
                  'Inconsistency in the number of Nedt channels:', nchanNedt,nchan
             call errorHalt(1)
          ENDIF
       ENDIF
       IF (icellres > 0) THEN
          !--- obtain noise averaging factor DATA
          CALL getNAF(nchan,ossChanID(1:nchan),ncell,nscan,nrf,Fnrf=Fnrf)
       ELSE
          ncell = 0
          nscan = 0
          nrf = 1.0
       ENDIF
       init=0
    END IF
    !---Computation of nedt (individual and total)
    IF (icellres > 0) THEN
       IF (iscanPos > nscan) THEN
          PRINT *,'err[NoiseModule::deviceNoise]:  iscan > nscan'
          call errorHalt(1)
       END IF
       IF (myiCellRes > ncell) THEN
          PRINT *,'err[NoiseModule::deviceNoise]:  iCell > ncell'
          call errorHalt(1)
       END IF
    ENDIF
    nch      = 0
    sum_rerr = 0.
    DO i=1,nchan
       ispec   = mapnedt(i)
       tAntn  = tb(i)*(1.-spoch(i))+CosmBckg*spoch(i)
       IF (DynamicNoise) THEN
          rerr(i) = dynNoise(ispec,tAntn,icaset)
          rerr(i) = rerr(i)/(1.-spoch(i))
       ENDIF
       IF (.not.(DynamicNoise)) rerr(i) = NedtStatic(i)
       avefac  = nrf(i,myiCellRes,iscanPos)
       rerr(i) = rerr(i)*avefac
       IF (kchan(i)) THEN
          nch      = nch + 1
          sum_rerr = sum_rerr+rerr(i)**2
       ENDIF
    ENDDO
    if (nch == 0) then
       sum_rerr = MISSING_REAL
    else
       sum_rerr=SQRT(sum_rerr/float(nch))
    endif
    RETURN
  END SUBROUTINE deviceNoise

  subroutine loadNedtStatic(file_nedt)
!!$    use ReadStdInputs, only: U_nedt
    USE ToolboxModule, Only: getUnit
    !---Input/Output variables
    character(len=*),               intent(in)    :: file_nedt
    !---Local variables
    integer                 :: j
    integer                 :: U_nedt

    if(loadNedtStaticDone) return

    U_nedt = getUnit()
    !------------------------------------------------------------------------
    ! Set global variables from file variables
    !------------------------------------------------------------------------
    open(U_nedt,file = file_nedt,status = 'old')
    read(U_nedt,*)
    read(U_nedt,*) nchanNedt

    read(U_nedt,*) (NedtStatic(j), j=1,nchanNedt)
    close(U_nedt)

    loadNedtStaticDone=.true.
    return
  end subroutine loadNedtStatic

    SUBROUTINE getNAF(nchan,chanID,ncell,nscan,nrf,Fnrf)
      use ReadStdInputs, only: F_nrf
      !---Input/Output variables
      INTEGER,                        INTENT(IN)    :: nchan
      CHARACTER(LEN=*), DIMENSION(:), INTENT(IN)    :: ChanID
      INTEGER,                        INTENT(INOUT) :: ncell,nscan
      REAL, DIMENSION(:,:,:),         INTENT(INOUT) :: nrf
      character(len=*), intent(in), optional :: Fnrf
      !---Local variables
      INTEGER                 :: i,k
      LOGICAL                 :: assumeOrder   ! for back compatibility
      character(len=mxLength) :: nrfFile

      !------------------------------------------------------------------------
      ! Set global variables from file variables
      !------------------------------------------------------------------------
      if (present(Fnrf)) then
         nrfFile = trim(adjustl(Fnrf))
      else
         nrfFile = trim(adjustl(F_nrf))
      end if
      call loadNAF(nrfFile)

      assumeOrder=.FALSE.
      IF (ALL(chanID(1:nchan) == OSS_NULL_ID)) assumeOrder=.TRUE.

      !------------------------------------------------------------------------
      ! Consistency checks
      !------------------------------------------------------------------------
      IF (assumeOrder .and. nchanNAF /= nchan) THEN
         PRINT *,'err[NoiseModule::getNAF]:  nchanNAF /= nchan: ',nchanNAF,nchan
         call errorHalt(1)
      END IF
      IF (nscanNAF > mxscans) THEN
         PRINT *,'err[NoiseModule::getNAF]:  mxscans not big enough'
         call errorHalt(1)
      END IF
      IF (ncellNAF >  mxcell) THEN
         PRINT *,'err[NoiseModule::getNAF]:  mxcell not big enough'
         call errorHalt(1)
      END IF

      !------------------------------------------------------------------------
      ! Copy to output arguments
      !------------------------------------------------------------------------
      ncell=ncellNAF
      nscan=nscanNAF

      !------------------------------------------------------------------------
      ! Find NAF data for each channel, by matching channel IDs
      !------------------------------------------------------------------------
      IF (assumeOrder) THEN
         nrf(1:nchanNAF,1:ncellNAF,1:nscanNAF)=nrfBuf
      ELSE
         nrf=MISSING_REAL
         DO k=1,nchan
            DO i=1,nchanNAF
               IF (NAFchanID(i) == chanID(k)) THEN
                  nrf(k,1:ncellNAF,1:nscanNAF)=nrfBuf(i,:,:)
                  EXIT
               END IF
            END DO
            IF (nrf(k,1,1) == MISSING_REAL) THEN
               PRINT *,'err[NoiseModule::getNAF]: ', &
                    'No NAF data found for channel ',chanID(k)
               call errorHalt(1)
            END IF
         END DO
      END IF
      DEALLOCATE(NAFchanID,nrfBuf)

      RETURN
    END SUBROUTINE getNAF

    subroutine loadNAF(file_nrf)
!!$      use ReadStdInputs, only: U_nrf
      USE ToolboxModule, Only: getUnit
      !---Input/Output variables
      character(len=*), intent(in)    :: file_nrf
      !---Local variables
      integer                         :: iscan,iscan0,i,j,ict
      real                            :: dum
      character(len=80)               :: cl
      character(len=7)                :: vid
      integer :: U_nrf
      if(loadNAFdone) return

      U_nrf = getUnit()
      ! Setting global variables from file variables
      OPEN(U_nrf,file=file_nrf,status='old')
      READ(U_nrf,'(a7,i4)')vid,nchanNAF
      IF (vid /= 'nchans=') THEN
         PRINT *,'err[NoiseModule::loadNAF]: ', &
              ' nchanNAF not found--incompatible file format'
         call errorHalt(1)
      END IF

      READ(U_nrf,'(7x,i4)')nscanNAF
      READ(U_nrf,'(7x,i4)')ncellNAF
      IF (ALLOCATED(NAFchanID)) DEALLOCATE(NAFchanID)
      IF (ALLOCATED(nrfBuf)) DEALLOCATE(nrfBuf)
      ALLOCATE(NAFchanID(nchanNAF),nrfBuf(nchanNAF,ncellNAF,nscanNAF))
      DO iscan=1,nscanNAF
         READ(U_nrf,'(a80)')cl
         READ(U_nrf,'(6x,i4)')iscan0
         IF (iscan .ne. iscan0) THEN
            PRINT *,'err[NoiseModule::getNAF]:  scans must be in order'
            call errorHalt(1)
         ENDIF
         READ(U_nrf,'(a80)')cl
         READ(U_nrf,'(a80)')cl
         DO i=1,nchanNAF
            READ(U_nrf,'(a12,4x,7(f5.3,2x))') NAFchanID(i), &
                 (nrfBuf(i,j,iscan),j=1,ncellNAF)
         END DO
      ENDDO
      CLOSE(U_nrf)

      loadNAFdone=.true.
      return
    end subroutine loadNAF

    SUBROUTINE getMapNedt(nnoise,tChanID,ChanID,nchan,mapnedt)
      !---Input/Output variables
      INTEGER,                             INTENT(IN)    :: nnoise
      CHARACTER(LEN_ID_SensorData), DIMENSION(:), POINTER       :: tChanID
      CHARACTER(LEN=*), DIMENSION(:),      INTENT(IN)    :: ChanID
      INTEGER,                             INTENT(IN)    :: nchan
      INTEGER,          DIMENSION(:),      INTENT(INOUT) :: mapnedt
      !---Local variables
      INTEGER      :: i,nChanAll,iformat

      mapnedt(:)=0
      DO i=1,nnoise
         WHERE (ChanID == tChanID(i))
            mapnedt = i
         END WHERE
      ENDDO
      DO i=1,nchan
         if (mapnedt(i) == 0) then
            print *,'err[NoiseModule::getMapNedt]: ', &
                 ' Noise data were not found for channel',i,chanID(i)
            call errorHalt(1)
         endif
      ENDDO
    END SUBROUTINE getMapNedt

    SUBROUTINE destroyNoise()
      CALL destroySensorData()
    END SUBROUTINE destroyNoise

  END MODULE NoiseModule
