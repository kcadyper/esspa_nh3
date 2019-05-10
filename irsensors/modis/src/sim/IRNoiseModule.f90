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

!-----------------------------------------------------
!
!  MODULE IRNoise: Contains subroutines needed for the
!                    IR noise generation.
!                  This version is specific to MODIS.
!                  AER Inc. 2008.
!-----------------------------------------------------
MODULE IRNoiseModule

! <f90Module>***********************************************************
!
! NAME:
!
!   IRNoiseModule
!
! PURPOSE:
!
!   Computations related to noise in IR channels.
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

  IMPLICIT NONE
  PRIVATE
  !------------------------------------------------------------------------
  !	Public items made available to calling program(s)
  !------------------------------------------------------------------------
  PUBLIC :: IRaddDeviceNoise,IRdeviceNoise,loadIRdeviceNoise

  !------------------------------------------------------------------------
  !     Declarations of private data (private to the module)
  !------------------------------------------------------------------------
  INTEGER,            parameter       :: mxchan=50
  INTEGER                             :: myNch
  REAL, DIMENSION(:), allocatable     :: AtmNoiseBT
  REAL, DIMENSION(:), allocatable     :: irDevNoise
  CHARACTER(len=12),DIMENSION(:),allocatable :: chanID

  ! loadDone implementation provides the option to call loadIRdeviceNoise() in advance
  logical                             :: loadDone=.false.

CONTAINS

  SUBROUTINE IRaddDeviceNoise(wvn,Y,DevNoise,sumIRDevNoise,kchan,nchan,chanIDreq)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   IRaddDeviceNoise
!
! PURPOSE:
!
!   Add Gaussian random noise to radiances.
!
! SYNTAX:
!
!   CALL IRaddDeviceNoise(wvn, Y, DevNoise, sumIRDevNoise, kchan,
!      nchan, chanIDreq)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   wvn            REAL     Wavenumbers
!   kchan          LOGICAL  Channel on/off mask
!   nchan          INTEGER  Number of channels
!   chanIDreq      CHAR     Channel identifier string for required
!                           channels
!
!   INPUTS/OUTPUTS:
!
!   Y              REAL     Radiometric data computed from state
!                           vector
!   DevNoise       REAL     Instrument noise standard deviation
!   sumIRDevNoise  REAL     Total instrument noise for channels
!                           turned on
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !---In/Out variables
    INTEGER,                        INTENT(IN)    :: nchan
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN)    :: chanIDreq
    LOGICAL,          DIMENSION(:), INTENT(IN)    :: kchan
    REAL,             DIMENSION(:), INTENT(IN)    :: wvn
    REAL,             DIMENSION(:), INTENT(INOUT) :: Y,DevNoise
    REAL,                           INTENT(INOUT) :: sumIRDevNoise
    !---Local variables
    INTEGER, SAVE         :: iseed=1287365
    INTEGER               :: i
    REAL                  :: stdt,gauss

    call IRdeviceNoise(wvn,Y,DevNoise,sumIRDevNoise=sumIRDevNoise, &
         kchan=kchan,nchan=nchan,chanIDreq=chanIDreq)
    DO i=1,nchan
       stdt=devNoise(i)
       Y(i)=Y(i)+gauss(iseed,stdt,0.0)
    END DO
    RETURN
  END SUBROUTINE IRaddDeviceNoise

  ! ======================================================================

  SUBROUTINE IRdeviceNoise(wvn,Yin,DevNoise,AtmNoiseR,sumIRDevNoise,kchan,nchan, &
        chanIDreq)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   IRdeviceNoise
!
! PURPOSE:
!
!   Determine noise standard deviation for each channel.
!
! SYNTAX:
!
!   CALL IRdeviceNoise(wvn, Yin, DevNoise, AtmNoiseR, sumIRDevNoise,
!      kchan, nchan, chanIDreq)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   wvn            REAL     Wavenumbers
!   Yin            REAL     Radiometric data for input
!   kchan          LOGICAL  Channel on/off mask
!   nchan          INTEGER  Number of channels
!   chanIDreq      CHAR     Channel identifier string for required
!                           channels
!
!   INPUTS/OUTPUTS:
!
!   DevNoise       REAL     Instrument noise standard deviation
!   AtmNoiseR*     REAL
!   sumIRDevNoise  REAL     Total instrument noise for channels
!                           turned on
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    use IRReadStdInputs, only: F_nedn
    !---In/Out variables
    INTEGER,                        INTENT(IN)              :: nchan
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN)              :: chanIDreq
    LOGICAL, DIMENSION(:),          INTENT(IN)              :: kchan
    REAL,    DIMENSION(:),          INTENT(IN)              :: wvn,Yin
    REAL,    DIMENSION(:),          INTENT(INOUT)           :: DevNoise
    REAL,    DIMENSION(:),          INTENT(INOUT), OPTIONAL :: AtmNoiseR
    REAL,                           INTENT(INOUT)           :: sumIRDevNoise
    !---Local variables
    REAL                            :: v,Ytmp,tb,bt,draddt
    INTEGER                         :: i,j,nch
    INTEGER, SAVE                   :: init=1
    INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: mapChan

    !------------------------------------------------------------------------
    ! Set global variables from file variables
    !------------------------------------------------------------------------

    call loadIRdeviceNoise(F_nedn)

    if(myNch /= nchan) then
       print*,'err[IRNoiseModule::IRdeviceNoise]: '
       print*,'Inconsistent number of channels: '
       print*,myNch,nchan
       call errorHalt(1)
    endif

    IF(init == 1)THEN
       allocate(AtmNoiseBT(nchan))
       AtmNoiseBT = 0.
       ! Map to the noise data for each "on" channel
       allocate (mapChan(nChan))
       if (all((.not.kchan) .or. (len_trim(chanIDreq) == 0))) then
          print *,'msg:[IRNoiseModule::IRdeviceNoise]: '// &
             'No channel IDs provided'
          print *,'  Assuming that noise '// &
              'data come in the same order as the OSS channels'
          mapChan=(/(i,i=1,nChan)/)
       else
          mapChan=0
          do i=1,nChan
             if (.not.kchan(i)) cycle
             do j=1,myNch
                if (trim(chanID(j)) == trim(chanIDreq(i))) then
                  mapChan(i)=j
                  exit
                endif
             enddo
             if (mapChan(i) == 0) then
                print*,'err[IRNoiseModule::IRdeviceNoise]: '
                print*,'No noise data found for channel:',i,chanIDreq(i)
                call errorHalt(1)
             endif
          enddo
       endif
       init       = 0
    endif

    !--------------------------
    !  Get IR noise amplitude:
    !--------------------------

    sumIRDevNoise=0.
    nch=0

    do i=1,nChan
       v=wvn(i)
       if(kchan(i))then
          YTmp=max(Yin(i),1.e-3)
          tb=bt(v,YTmp)
          devNoise(i) = irDevNoise(mapChan(i))
          nch=nch+1
          sumIRDevNoise=sumIRDevNoise+devNoise(i)**2
          if(present(AtmNoiseR))AtmNoiseR(i)=draddt(v,tb)*AtmNoiseBT(i)
       else
          devNoise(i) = 0.
          if(present(AtmNoiseR))AtmNoiseR(i)=0.
       end if
    enddo
    sumIRDevNoise=sqrt(sumIRDevNoise/nch)

    return

  END SUBROUTINE IRdeviceNoise

  ! ======================================================================

  SUBROUTINE loadIRdeviceNoise(F_nedn)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   loadIRdeviceNoise
!
! PURPOSE:
!
!   Load noise data into global arrays.
!
! SYNTAX:
!
!   CALL loadIRdeviceNoise(F_nedn)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   F_nedn  CHAR  File path for channel noise data
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>


    USE IRReadStdInputs, ONLY: U_nedn

    ! Load noise data into global arrays.

    CHARACTER(len=*), INTENT(IN) :: F_nedn

    INTEGER           :: iChan

    if(loadDone) return

    open(U_nedn,file=F_nedn,status='old')
    read(U_nedn,*) myNch
    allocate(irDevNoise(myNch))
    allocate(chanID(myNch))

    DO iChan=1,myNch
       read(U_nedn,*) chanID(iChan),irDevNoise(iChan)
    ENDDO

    close(U_nedn)
    loadDone=.true.

    return

  END SUBROUTINE loadIRdeviceNoise

END MODULE IRNoiseModule
