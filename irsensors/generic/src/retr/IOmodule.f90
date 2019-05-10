MODULE IOmodule
  !
  ! Module that loads the sensor data and ancillary data, and write out
  ! retrieval results
  !
  ! Subroutines
  !
  !     getMWsensorData
  !     getIRsensorData
  !     getImgSensorData
  !     getAncillaryData
  !     destroyMWsensorData
  !     destroyIRsensorData
  !     destroyImgSensorData
  !     destroyancillaryData
  !     putRetr
  !
  ! USE:
  !
  !     netcdf_module
  !     RetrModule
  !     MWobsStructure
  !     rad_io_module
  !     RadFile
  !
  ! yhe@aer.com, 03/16/2016
  !
  USE ncdf_module, Only: &
       nf_open, &
       nf_inq_varid, &
       nf_nowrite, &
       nf_noerr, &
       writeNcdfData, &
       writeNcdfAttr, &
       closencdffile

  USE scene_io_module, Only: &
       openScene, &
       getScene, &
       putScene, &
       closeScene, &
       TB_UNITS, &
       RADMW_UNITS

  USE rad_io_module, Only: &
       openRad, &
       closeRad, &
       getRad, &
       queryRad

  USE bkg_io_module, Only: &
       queryCov, &
       openCov,  &
       getCov,   &
       closeCov

  USE RadFile, Only: &
       getDimRad, &
       getRadData

  USE ControlStructure, Only: &
       GenControl_t, &
       StateSpaceDefin_t

  USE OutputStructure, Only: &
       OutputGran_t, &
       LinearGran_t

  USE PrimaryRetrModule, Only: &
       PrimaryRetrControl_t

  USE MWobsStructure, Only: &
       MWgran_t

  USE IRobsStructure, Only: &
       IRgran_t

  USE ImgObsStructure, Only: &
       ImgDef_t, &
       ImgOb_t, &
       ImgGran_t

  USE AncillaryStructure, Only: &
       AncillaryGran_t

  USE StateIndexModule, Only: &
       initLengths, &
       getVectorLength, &
       StateIndex_t, &
       maxMol, &
       nullMolID, &
       compressMol, &
       genIndices

  USE constants, Only: &
       MISSING_INT, &
       MISSING_REAL, &
       MISSING_CHAR

  USE ToolboxModule, Only: &
       getUnit

  USE VertCoord, ONLY: &
       mxCoordTyp, &
       Pcoord

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: &
       getMWsensorData, &
       getIRsensorData, &
       getImgSensorData, &
       getAncillaryData, &
       putRetr, &
       putLinear, &
       destroyMWsensorData, &
       destroyIRsensorData, &
       destroyImgSensorData, &
       destroyancillaryData

  integer, parameter :: mxLength = 256
  logical :: loadAuxDone=.false.
  logical :: loadExtDone=.false.
  logical :: loadEmisDone=.false.
  logical :: dbg = .false.
  logical :: prn = .false.

CONTAINS

  ! Would the getMWsensorData, etc, be made flexible to load any
  ! number of FORs/FOVs by the control parameters?
  subroutine getMWsensorData(genControl, MWgran)

    TYPE(genControl_t), intent(in) :: genControl
    TYPE(MWgran_t), intent(out) :: MWgran
    integer :: U_Rad
    integer :: tmpvarid
    integer :: ncStatus
    integer :: ifov
    integer :: ich
    integer :: nAng
    real, dimension(:), allocatable :: eia
    character (len=mxLength) :: procName

    procName = ' [IOmodule::getMWsensorData]: '
    if (dbg) print *, trim(procName)//'starting ...'

    if (.not. genControl%MWdataOn) then
       MWgran%nFOV = 0
       MWgran%def%nChan = 0
       allocate(MWgran%ob(1))
       return
    endif

    call queryRad(fname=trim(genControl%MWobsFile),&
         nchan=MWgran%def%nChan,nang=nAng)
    if (.not. allocated(MWgran%def%pol)) allocate(MWgran%def%pol(MWgran%def%nChan))
    if (.not. allocated(MWgran%def%frq)) allocate(MWgran%def%frq(MWgran%def%nChan))
    if (.not. allocated(MWgran%def%chanID)) allocate(MWgran%def%chanID(MWgran%def%nChan))
    if (.not. allocated(MWgran%def%planckAlpha)) allocate(MWgran%def%planckAlpha(MWgran%def%nChan))
    if (.not. allocated(MWgran%def%planckBeta)) allocate(MWgran%def%planckBeta(MWgran%def%nChan))
    !temporary
    if (nAng > 1) then
       print *, trim(procName)//'error - nAng greater than 1'
       stop
    end if
    allocate(eia(nAng))

    if (dbg .and. prn) then
       print *, trim(procName)//' nchan = ',MWgran%def%nChan
       print *, trim(procName)//' nAng = ',nAng
    end if

    call openRad(ncid=U_Rad,fname=trim(genControl%MWobsFile), &
         pol=MWgran%def%pol, &
         frq=MWgran%def%frq, &
         nrec=MWgran%nFOV,status='old')

    if (dbg .and. prn) then
       print *, trim(procName)//' nFOV = ',MWgran%nFOV
       print *, trim(procName)//' ich, pol, frq = '
       do ich = 1, MWgran%def%nChan
          print *, ich, MWgran%def%pol(ich),MWgran%def%frq(ich)
       end do
    end if

    !place-holders for Jacobians and emissivity for later
    !processing
    allocate(MWgran%ob(MWgran%nFOV))
    do ifov = 1, MWgran%nFOV
       if (.not. allocated(MWgran%ob(ifov)%rad)) allocate(MWgran%ob(ifov)%rad(MWgran%def%nChan))
       if (.not. allocated(MWgran%ob(ifov)%QC)) allocate(MWgran%ob(ifov)%QC(MWgran%def%nChan))
    end do

    !Decide if it is brightness or radiance stored in the file by
    !inquiry of the variable name in the file - temporary, esp. it is
    !hardwired name.
    ncStatus = nf_inq_varid(U_Rad,"tb",tmpvarid)
    if (ncStatus == nf_noerr) then
       MWgran%def%radOrTb=0
    else
       ncStatus = nf_inq_varid(U_Rad,"rad",tmpvarid)
       if (ncStatus == nf_noerr) then
          MWgran%def%radOrTb=1
       else
          print *, trim(procName)//'error - no Rad or TB found'
          stop
       endif
    endif

    if (dbg .and. prn) print *, trim(procName)//'RadOrTb = ',MWgran%def%radOrTb

    call getRad(ncid=U_Rad,irec=1,time=MWgran%def%startTime)
    call getRad(ncid=U_Rad,irec=MWgran%nFOV,time=MWgran%def%endTime)

    if (dbg .and. prn) then
       print *, trim(procName)//'startTime = ',MWgran%def%startTime
       print *, trim(procName)//'endTime = ',MWgran%def%endTime
    end if

    do ifov = 1, MWgran%nFOV
       call getRad(ncid=U_Rad,irec=ifov,&
            lat=MWgran%ob(ifov)%lat, &
            lon=MWgran%ob(ifov)%lon, &
            eia=EIA)
       MWgran%ob(ifov)%EIA = eia(1)
       IF (MWgran%def%radOrTb == 0) THEN
          CALL getRad(ncid=U_Rad,irec=ifov,tb=MWgran%ob(ifov)%rad)
       ELSE
          CALL getRad(ncid=U_Rad,irec=ifov,rad=MWgran%ob(ifov)%rad)
       ENDIF
       MWgran%ob(ifov)%QC = 0
       MWgran%ob(ifov)%EAA = 0.
       MWgran%ob(ifov)%scanAng = 0.
       MWgran%ob(ifov)%scanPos = 0

       if (dbg .and. prn) then
          print *, trim(procName), &
               MWgran%ob(ifov)%lat, &
               MWgran%ob(ifov)%lon, &
               MWgran%ob(ifov)%EIA
          print *, trim(procName)//'Rad/TB = '
          do ich = 1, MWgran%def%nChan
             print *, ich, MWgran%def%pol(ich),MWgran%def%frq(ich), MWgran%ob(ifov)%rad(ich)
          end do
       end if

    end do

    call closeRad(U_rad)

    deallocate(eia)

    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine getMWsensorData

  subroutine getIRsensorData(genControl,IRgran)
    use IRsensorReader, only: readIRsensorFile
    TYPE(genControl_t), intent(in) :: genControl
    TYPE(IRgran_t), intent(out) :: IRgran
    integer :: U_Rad
    integer :: nAng
    real, dimension(:), allocatable :: eiaAll
    real, dimension(:), allocatable :: latAll
    real, dimension(:), allocatable :: lonAll
    real, dimension(:,:), allocatable :: radTmp
    integer :: ifov, ifor, ich
    character (len=mxLength) :: procName

    procName = ' [IOmodule::getIRsensorData]: '

    if (dbg) print *, trim(procName)//'starting ...'

    if (.not. genControl%IRdataOn) then
       IRgran%nFOR = 0
       IRgran%def%nChan = 0
       if (.not. allocated(IRgran%FOR)) allocate(IRgran%FOR(1))
       if (.not. allocated(IRgran%FOR(1)%ob)) allocate(IRgran%FOR(1)%ob(1))
       return
    endif
    if (readIRsensorFile(genControl%IRobsFile, IRgran, genControl%independentFOV,genControl%apodizationType)) then
      print *,'IR sensor read'
      print *,'NFOR',  size(IRgran%FOR)
      print *,'NFOV',  size(IRgran%FOR(1)%ob(:))
      return
    end if  

    !---Get radiance file dimensions and allocate arrays:
    call getDimRad(U_Rad,trim(genControl%IRobsFile), &
         IRgran%def%nChan, &
         IRgran%nFOR,nAng) !Is nAng necessary? It is 1 in the input file.

    if (dbg .and. prn) then
       print *, trim(procName)//"nFOR  = " , IRgran%nFOR
       print *, trim(procName)//"nChan = " , IRgran%def%nChan
    end if

    IRgran%def%nFOV = 1  !temporary until available in the input file
    if (.not. allocated(IRgran%def%wvn)) allocate(IRgran%def%wvn(IRgran%def%nChan))
    if (.not. allocated(IRgran%def%chanNum)) allocate(IRgran%def%chanNum(IRgran%def%nChan))
    if (.not. allocated(IRgran%def%chanID)) allocate(IRgran%def%chanID(IRgran%def%nChan))
    if (.not. allocated(IRgran%FOR)) allocate(IRgran%FOR(IRgran%nFOR))
    do ifor = 1, IRgran%nFOR
       if (.not. allocated(IRgran%FOR(ifor)%ob)) &
            allocate(IRgran%FOR(ifor)%ob(IRgran%def%nFOV))
       do ifov = 1, IRgran%def%nFOV
          if (.not. allocated(IRgran%FOR(ifor)%ob(ifov)%rad)) &
               allocate(IRgran%FOR(ifor)%ob(ifov)%rad(IRgran%def%nChan))
          if (.not. allocated(IRgran%FOR(ifor)%ob(ifov)%QC)) &
               allocate(IRgran%FOR(ifor)%ob(ifov)%QC(IRgran%def%nChan))
       end do
    end do

    ! this part may not be necessary, if decided that genControl does not
    ! overwrite the specifications in the file, same as in genMWsensorData

!!$    nFOR = genControl%nFORmax
!!$    if (nFOR < 1) then
!!$       nFOR = IRgran%nFOR
!!$       print *, trim(procName)//"reset nFOR in genControl to comply with IR"
!!$    end if
!!$
!!$    if (IRgran%nFOR > nFOR) then
!!$       IRgran%nFOR = nFOR
!!$       print *, trim(procName)//"reset IR nFOR to comply with genControl"
!!$    end if

    ! fill with dummy values while getRadData does not support reading times
    IRgran%def%startTime(:)=9999
    IRgran%def%endTime(:)=9999

    !temporary until the operational data is available
    allocate(radTmp(IRgran%def%nChan,IRgran%nFOR))
    allocate(eiaAll(IRgran%nFOR),latAll(IRgran%nFOR),lonAll(IRgran%nFOR))

    !---Load radiance data: missing variables in test = EAA, SolIA, SolAA, scanAng, QC
    call getRadData(U_Rad,trim(genControl%IRobsFile), &
         radTmp,eiaAll,latAll,lonAll,IRgran%def%wvn, &
         chanID=IRgran%def%chanID)

    IRgran%def%chanNum = (/(ich, ich=1,IRgran%def%nChan)/)

    do ifor = 1, IRgran%nFOR
       do ifov = 1, IRgran%def%nFOV
          IRgran%FOR(ifor)%ob(ifov)%rad = radTmp(:,ifor)
          IRgran%FOR(ifor)%ob(ifov)%lat = latAll(ifor)
          IRgran%FOR(ifor)%ob(ifov)%lon = lonAll(ifor)
          IRgran%FOR(ifor)%ob(ifov)%eia = eiaAll(ifor)
          !---missing variables in test = EAA, SolIA, SolAA, scanAng, QC
          IRgran%FOR(ifor)%ob(ifov)%QC = 0
          IRgran%FOR(ifor)%ob(ifov)%EAA = 0.
          IRgran%FOR(ifor)%ob(ifov)%SolIA = 90.
          IRgran%FOR(ifor)%ob(ifov)%SolAA = 0.
          IRgran%FOR(ifor)%ob(ifov)%scanAng = 0.
       end do
    end do

    if (dbg .and. prn) then
       print *, trim(procName)//'ich, chanID, wvn = '
       do ich = 1, IRgran%def%nChan
          print *, ich, trim(IRgran%def%chanID(ich)), IRgran%def%wvn(ich)
       end do

       print *, trim(procName)//'ifor, ifov, lat, lon, eia rad = '
       do ifor = 1, 2
          do ifov = 1, IRgran%def%nFOV
             print *, ifor, ifov, &
                  IRgran%FOR(ifor)%ob(ifov)%lat, &
                  IRgran%FOR(ifor)%ob(ifov)%lon, &
                  IRgran%FOR(ifor)%ob(ifov)%eia
             do ich = 1, IRgran%def%nChan
                print *, ich, IRgran%FOR(ifor)%ob(ifov)%rad(ich)
             end do
          end do
       end do
    end if

    if (dbg) print *, trim(procName)//'ending ...'
    deallocate(radTmp)
    deallocate(eiaAll)
    deallocate(latAll)
    deallocate(lonAll)
  end subroutine getIRsensorData

  subroutine getImgSensorData(genControl, ImgGran)
    TYPE(genControl_t), intent(in) :: genControl
    TYPE(ImgGran_t), intent(out) :: ImgGran
    !Local
    integer :: U_Rad
    integer, dimension(:,:), allocatable :: cldTypes
    character (len=mxLength) :: procName

    procName = ' [IOmodule::getImgSensorData]: '
    if (dbg) print *, trim(procName)//'starting ...'

    if (.not. genControl%imgDataOn) then
       imgGran%nFOR=0
       imgGran%def%nFOV=0
       if (.not. allocated(imgGran%ob)) allocate(imgGran%ob(1,1))
       return
    endif

    !temporary until cloud specific reader is available and FOV is
    !available
    print *, trim(procName)//'openRad ...'
    call openRad(ncid=U_Rad,fname=trim(genControl%IRobsFile), &
         nRec=imgGran%nFOR)
    imgGran%def%nFOV = 1  !temporary, could be made the same as IRgran%def%nFOV, if needed
    if (.not. allocated(imgGran%ob)) &
         allocate(imgGran%ob(imgGran%def%nFOV,imgGran%nFOR))
    if (.not. allocated(cldTypes)) allocate(cldTypes(imgGran%def%nFOV,imgGran%nFOR))
    !currently there is no subroutine to load 2-D cloud image
    !data. The subroutine getCloudData in the RadFile handles 1-D
    !cloud type vector with length of nFOR
    print *, '[IOmodule::getImgSensorData]: calling getCloudData ... to be implemented'
!!$    call getCloudData(U_Rad,trim(genControl%IRobsFile),cldTypes)
!!$    do iFOR = 1, imgGran%nFOR
!!$       do iFOV = 1, imgGran%def%nFOV
!!$          imgGran%ob(iFOV,iFOR)%cloudType = cldTypes(iFOV,iFOR)
!!$       end do
!!$    end do

    if (allocated(cldTypes)) deallocate(cldTypes)
    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine getImgSensorData

  subroutine getAncillaryData(genControl, ancillaryGran)
    TYPE(genControl_t), intent(in) :: genControl
    TYPE(AncillaryGran_t), intent(out) :: ancillaryGran
    !Local
    TYPE(StateIndex_t) :: IGAux,NGAux
    logical :: is_ncdf
    logical :: file_exists = .FALSE.
    integer :: ncidtmp,openStatus
    integer :: nFORaux
    integer :: U_Aux
    integer :: U_bkgrdExtern
    integer :: U_IRemissExt
    integer :: nforExt
    integer :: nforEmis,nemirExt=0
    integer :: nqcEmis
    integer :: ifor, ifov
    integer, dimension (:), allocatable  :: qcEmis
    real :: pland
    real, dimension(:), allocatable :: psfcTmp
    real, dimension(:), allocatable :: xExtTmp
    real, dimension(:), allocatable :: xEmisTmp
    character(len=8) :: dummy
    character (len=mxLength) :: procName

    procName = ' [IOmodule::getAncillaryData]: '
    if (dbg) print *, trim(procName)//'starting ...'

    ancillaryGran%nFOR = MISSING_INT !initialized
    ancillaryGran%def%nFOV = 1 !temporary until provided by the file

    ancillaryGran%def%MolID = nullMolID
    if (.not. loadExtDone) then
       !-- loadExt: atmosphere profiles
       if (dbg .and. prn) print *, trim(procName)//' loading ... '//trim(genControl%ancillaryAtmosFile)

       inquire(FILE=trim(genControl%ancillaryAtmosFile),EXIST=file_exists)
       if (.not. file_exists) then
          print *, trim(procName)//' warning - genControl%ancillaryAtmosFile not exist : ',&
                   trim(genControl%ancillaryAtmosFile)
          ancillaryGran%def%NG=initLengths()
       end if

       if (ANY((/genControl%externDataFlags%temp, &
            genControl%externDataFlags%mol, &
            genControl%externDataFlags%pSfc, &
            genControl%externDataFlags%wind, &
            genControl%externDataFlags%cloudLiq, &
            genControl%externDataFlags%cloudIce/)) .and. file_exists) then
          call openScene(ncid=U_bkgrdExtern, &
               file=trim(genControl%ancillaryAtmosFile), &
               nprf=nforExt, &
               molID=ancillaryGran%def%MolID, &
               status='old', &
               IG=ancillaryGran%def%IG, &
               NG=ancillaryGran%def%NG)

          if (nforExt < genControl%nFORmax) then
             print *,trim(procName)//' err - Not enough external data;  nforExt, genControl%nFORmax', &
                  nforExt, genControl%nFORmax
             call errorHalt(1)
          endif

          ancillaryGran%nFOR=nforExt
          ancillaryGran%def%nParG = getVectorLength(ancillaryGran%def%NG)
          if (.not. allocated(ancillaryGran%FOR)) &
               allocate(ancillaryGran%FOR(ancillaryGran%nFOR))
          do ifor = 1, ancillaryGran%nFOR
             if (.not. allocated(ancillaryGran%FOR(ifor)%data)) &
                  allocate(ancillaryGran%FOR(ifor)%data(ancillaryGran%def%nFOV))
             do ifov = 1, ancillaryGran%def%nFOV
                if (.not. allocated(ancillaryGran%FOR(ifor)%data(ifov)%stateV)) &
                     allocate(ancillaryGran%FOR(ifor)%data(ifov)%stateV(ancillaryGran%def%nParG))
             end do
          end do
          allocate (xExtTmp(ancillaryGran%def%nParG))
          IGaux = ancillaryGran%def%IG
          do ifor=1,ancillaryGran%nFOR
             CALL getScene(ncid=U_bkgrdExtern,irec=ifor,x=xExtTmp,pland=pland)
             do ifov = 1, ancillaryGran%def%nFOV
                ancillaryGran%FOR(ifor)%data(ifov)%stateV(:)=xExtTmp(:)
                ancillaryGran%FOR(ifor)%data(ifov)%psfc = xExtTmp(IGAux%psfc)
                ancillaryGran%FOR(ifor)%data(ifov)%landFrac = pland
                ancillaryGran%FOR(ifor)%data(ifov)%elev = MISSING_REAL
             end do
          enddo
          call closeScene(U_bkgrdExtern)
          deallocate(xExtTmp)
          loadExtDone = .true.
          loadAuxDone = .true.
       endif
    end if

    if (.not. loadAuxDone) then
       !-- loadAux: surface data
       inquire(FILE=trim(genControl%ancillarySfcFile),EXIST=file_exists)
       if (.not. file_exists) then
          print *, trim(procName)//' err - genControl%ancillarySfcFile not exist'
          call errorHalt(1)
       end if
       if (file_exists) then
          if (dbg .and. prn) print *, trim(procName)//' loading ... '// trim(genControl%ancillarySfcFile)
          openStatus = nf_open(trim(genControl%ancillarySfcFile),nf_nowrite,ncidtmp)
          is_ncdf = (openStatus == nf_noerr)
          if ( is_ncdf ) call closeNcdfFile(ncidtmp)
       end if

       if ( is_ncdf ) then
          !---Reading of netCDF auxiliary file:
          ! temporary solution until a more suitable I/O interface for
          ! auxiliary data is in place
          call openScene(ncid=U_Aux,file=trim(genControl%ancillarySfcFile), &
               nprf=nFORaux,IG=IGAux,NG=NGAux,status='old')
          ancillaryGran%nFOR=nFORaux
          if (.not. allocated(ancillaryGran%FOR)) &
               allocate(ancillaryGran%FOR(ancillaryGran%nFOR))
          do ifor = 1, ancillaryGran%nFOR
             if (.not. allocated(ancillaryGran%FOR(ifor)%data)) &
                  allocate(ancillaryGran%FOR(ifor)%data(ancillaryGran%def%nFOV))
          end do
          if (.not. allocated(psfcTmp)) allocate(psfcTmp(getVectorLength(NGAux)))
          do ifor = 1, ancillaryGran%nFOR
             call getScene(ncid=U_Aux,irec=ifor,x=psfcTmp,pland=ancillaryGran%FOR(ifor)%landFrac)
             ancillaryGran%FOR(ifor)%psfc = psfcTmp(IGAux%pSfc)
             do ifov = 1, ancillaryGran%def%nFOV
                ancillaryGran%FOR(ifor)%data(ifov)%psfc = ancillaryGran%FOR(ifor)%psfc
                ancillaryGran%FOR(ifor)%data(ifov)%landFrac = ancillaryGran%FOR(ifor)%landFrac
                ancillaryGran%FOR(ifor)%data(ifov)%elev = MISSING_REAL
             end do
          enddo
          deallocate(pSfcTmp)
          call closeScene(U_Aux)
       else
          !---Reading of ASCII auxiliary file, for backward compatibility:
          U_Aux = getUnit()
          open(U_Aux, file=trim(genControl%ancillarySfcFile), status='old')
          read(U_Aux,'(2i6)') nFORaux
          ancillaryGran%nFOR=nFORaux
          if (.not. allocated(ancillaryGran%FOR)) &
               allocate(ancillaryGran%FOR(ancillaryGran%nFOR))
          do ifor = 1, ancillaryGran%nFOR
             if (.not. allocated(ancillaryGran%FOR(ifor)%data)) &
                  allocate(ancillaryGran%FOR(ifor)%data(ancillaryGran%def%nFOV))
          end do
          do ifor = 1, ancillaryGran%nFOR
             read(U_Aux,'(a8,2f10.3)') &
                  dummy, &
                  ancillaryGran%FOR(ifor)%psfc, &
                  ancillaryGran%FOR(ifor)%landFrac
             !-temporary until values at each FOV is available
             do ifov = 1, ancillaryGran%def%nFOV
                ancillaryGran%FOR(ifor)%data(ifov)%psfc = ancillaryGran%FOR(ifor)%psfc
                ancillaryGran%FOR(ifor)%data(ifov)%elev = MISSING_REAL
                ancillaryGran%FOR(ifor)%data(ifov)%landFrac = ancillaryGran%FOR(ifor)%landFrac
             end do
          enddo
          close(U_Aux)
       endif

       loadAuxDone = .true.
    end if
    if (.not. loadEmisDone) then
       !-- omit emMW
       ancillaryGran%def%nEmMW = 0
       if (ancillaryGran%def%nEmMW /= 0) then
          if (.not. allocated(ancillaryGran%def%frq)) &
               allocate(ancillaryGran%def%frq(ancillaryGran%def%nEmMW))
          if (.not. allocated(ancillaryGran%def%pol)) &
               allocate(ancillaryGran%def%pol(ancillaryGran%def%nEmMW))
       end if

       !-- loadEmis: IR emissivity
       if (dbg .and. prn) print *, trim(procName)//' loading ... '//trim(genControl%ancillaryIRemisFile)
       inquire(FILE=trim(genControl%ancillaryIRemisFile),EXIST=file_exists)
       if (.not. file_exists) then
          print *, trim(procName)//' warning - genControl%ancillaryIRemisFile not exist : ', &
                   trim(genControl%ancillaryIRemisFile)
       end if

       if (genControl%externDataFlags%emRfIR .and. file_exists) then
          call queryCov(file=trim(genControl%ancillaryIRemisFile), &
               nbkg=nforEmis,nemir=nemirExt,nqc=nqcEmis)

          if (nforEmis < genControl%nFORmax) then
             print *,trim(procName)//' err - Not enough external emissivity data; '//&
                  'nforEmis, genControl%nFORmax',nforExt, genControl%nFORmax
             call errorHalt(1)
          endif

          if (nforEmis /= ancillaryGran%nFOR) then
             print *, trim(procName)//' err - inconsistent nFOR b/w emissivity and '//&
                  'surface data', nforEmis, ancillaryGran%nFOR
             call errorHalt(1)
          end if

          ancillaryGran%def%nEmIR = nemirExt
          !Purposely expand the length to twice as long for storing
          !reflectivity later on in BackgroundModule
          ancillaryGran%def%nEmRfIR = 2*nemirExt !per instruction 6/22/2016
          if (.not. allocated(ancillaryGran%def%wvn)) &
               allocate(ancillaryGran%def%wvn(nemirExt))
          call openCov(ncid=U_IRemissExt, &
               file=trim(genControl%ancillaryIRemisFile), &
               frqEmir=ancillaryGran%def%wvn,status='old')

          !-- Not a replicate of profiles or surface data, because it
          !-- is an independent process.
          if (.not. allocated(ancillaryGran%FOR)) &
               allocate(ancillaryGran%FOR(ancillaryGran%nFOR))
          do ifor = 1, ancillaryGran%nFOR
             if (.not. allocated(ancillaryGran%FOR(ifor)%data)) &
                  allocate(ancillaryGran%FOR(ifor)%data(ancillaryGran%def%nFOV))
             if (ancillaryGran%def%nEmMW /= 0) then
                if (.not. allocated(ancillaryGran%FOR(ifor)%emMW)) &
                     allocate(ancillaryGran%FOR(ifor)%emMW(ancillaryGran%def%nEmMW))
             end if
             do ifov = 1, ancillaryGran%def%nFOV
                if (.not. allocated(ancillaryGran%FOR(ifor)%data(ifov)%emRfIR)) &
                     allocate(ancillaryGran%FOR(ifor)%data(ifov)%emRfIR(ancillaryGran%def%nEmRfIR))
             end do
          end do
          !genControl%ancillaryIRemisFile contains emissivity with
          !length of nemirExt, which is the same as ancillaryGran%def%wvn
          allocate (xEmisTmp(2*nemirExt),qcEmis(nqcEmis))
          do ifor = 1, ancillaryGran%nFOR
             CALL getCov(ncid=U_IRemissExt,irec=ifor,dmEmIR=xEmisTmp,qc=qcEmis)
             do ifov = 1, ancillaryGran%def%nFOV
                ancillaryGran%FOR(ifor)%data(ifov)%emRfIR(1:2*nemirExt) = xEmisTmp(:)
                if (btest(qcEmis(1),0)) &
                     ancillaryGran%FOR(ifor)%data(ifov)%emRfIR(1:2*nemirExt) = MISSING_REAL
             end do
          enddo
          call closeCov(U_IRemissExt)
          deallocate(xEmisTmp,qcEmis)
          loadEmisDone = .true.
       endif
    end if

    if (dbg) print *, trim(procName)//' ending ...'

  end subroutine getAncillaryData

  subroutine putRetr(&
       inFile, &
       stateSpaceDefin, &
       retrOutput &
       )
    character (len=*), intent(in) :: inFile
    TYPE(StateSpaceDefin_t), intent(in) :: stateSpaceDefin
    TYPE(OutputGran_t), intent(in) :: retrOutput
    !Local
    TYPE(StateIndex_t) :: NGatm
    TYPE(StateIndex_t) :: IGatm
    integer :: nParGatm
    integer :: nChanMW
    integer :: nChanIR
    integer :: nChan
    integer :: nEmIR
    integer :: iFOR
    integer :: iFOV
    integer :: ncid
    integer :: nLev
    integer, dimension(6) :: time
    integer, dimension(:), allocatable :: MolID
    integer, dimension(:), allocatable :: polMW
!!$    real :: sigPtop
    real, dimension(1) :: EIA
    real, dimension(1) :: rms
    real, dimension(1) :: chiSq
    real, dimension(:), allocatable :: frqMW
    real, dimension(:), allocatable :: wvnIR
    character (len=mxCoordTyp) :: vCoordTyp
!!$    character (len=mxLength) :: vCoordA
!!$    character (len=mxLength) :: vCoordB
    character (len=mxLength) :: ncFile
    character (len=mxLength) :: procName

    procName = '[IOmodule::putRetr]:'
    if (dbg) print *, trim(procName)//' starting ...'

    if (.not. allocated(MolID)) allocate(MolID(size(stateSpaceDefin%MolID)))
    ncFile = trim(infile)
    MolID = stateSpaceDefin%MolID
    NGatm = initLengths(Natmos=retrOutput%def%NG)
    IGatm = genIndices(NGatm)
    nParGAtm = getVectorLength(NGatm)
    vCoordTyp = trim(stateSpaceDefin%vCoordTyp)
    nChanMW = stateSpaceDefin%nChanMW
    nChanIR = size(stateSpaceDefin%wvnIR)
    nChan = nChanMW + nChanIR
    nEmIR = stateSpaceDefin%nEmIR
    !-need to be integrated into the header definition
!!$    vCoordA=stateSpaceDefin%vCoordA
!!$    vCoordB=stateSpaceDefin%vCoordB
!!$    sigPtop=stateSpaceDefin%sigPtop

    !Define a local variable to comply with the required (inout)
    !attribute in openScene call
    nLev = stateSpaceDefin%nLev
    if (allocated(stateSpaceDefin%polMW)) then
       if (.not. allocated(polMW)) allocate(polMW(size(stateSpaceDefin%polMW)))
       polMW = stateSpaceDefin%polMW
    end if
    if (allocated(stateSpaceDefin%frqMW)) then
       if (.not. allocated(frqMW)) allocate(frqMW(size(stateSpaceDefin%frqMW)))
       frqMW = stateSpaceDefin%frqMW
    end if
    if (allocated(stateSpaceDefin%wvnIR)) then
       if (.not. allocated(wvnIR)) allocate(wvnIR(size(stateSpaceDefin%wvnIR)))
       wvnIR = stateSpaceDefin%wvnIR
    end if

    print *,'ncfile =',ncfile
    call openScene(ncid=ncid,file=ncFile, &
         MolId=MolID, &
         polarity=polMW, &
         freq=frqMW, &
         casename=ncFile, &
         nLevel=nLev, &
         nchmw=nChanMW, &
         vCoordTyp=vCoordTyp, &
         IG=IGatm, &
         NG=NGatm, &
         sfcGrid=wvnIR, &
         nSfcGrid=nEmIR, &
         status='new')

    if (trim(vCoordTyp) == Pcoord) &
         call writeNcdfAttr(ncid,attrName='standardPressureGrid',attr=stateSpaceDefin%pRef)

    !- temporary, until all the sources are resolved
    time = (/2016,7,28,16,21,00/) !will be in the header for the beginning and ending times of a granuale
    if (allocated(retrOutput%FOR)) then
       do iFOR=1, retrOutput%nFOR
          EIA(1) = retrOutput%FOR(iFOR)%EIA
          rms(1) = retrOutput%FOR(iFOR)%rms
          chiSq(1) = retrOutput%FOR(iFOR)%chiSq
          call putScene(ncid=ncid, &
               lat=retrOutput%FOR(iFOR)%lat, &
               lon=retrOutput%FOR(iFOR)%lon, &
               x=retrOutput%FOR(iFOR)%state(1:nParGatm), &
               eia=EIA, &
               pland=retrOutput%FOR(iFOR)%landFrac, &
               landType=retrOutput%FOR(iFOR)%landType &
               )
          if (retrOutput%def%nEmMW > 0) &
             call putScene(ncid=ncid, &
                  emMW=retrOutput%FOR(iFOR)%emMW &
                  )
          call putScene(ncid=ncid, &
               rms=rms, &
               chisq=chiSq, &
               iter=retrOutput%FOR(iFOR)%nIter, &
               roUnits=retrOutput%def%roUnits &
               )
          if (trim(vCoordTyp) /= Pcoord) call putScene(ncid=ncid,pressure=retrOutput%FOR(iFOR)%press)
       end do
    end if

    !- corresponding writing subprogram
    if (allocated(retrOutput%FOR2D)) then
       do iFOR=1, retrOutput%nFOR
          do iFOV=1,retrOutput%def%nFOV
             EIA(1) = retrOutput%FOR2D(iFOV,iFOR)%EIA
             rms(1) = retrOutput%FOR2D(iFOV,iFOR)%rms
             chiSq(1) = retrOutput%FOR2D(iFOV,iFOR)%chiSq
             call putScene(ncid=ncid, &
                  lat=retrOutput%FOR2D(iFOV,iFOR)%lat, &
                  lon=retrOutput%FOR2D(iFOV,iFOR)%lon, &
                  x=retrOutput%FOR2D(iFOV,iFOR)%state(1:nParGatm), &
                  eia=EIA, &
                  pland=retrOutput%FOR2D(iFOV,iFOR)%landFrac, &
                  landType=retrOutput%FOR2D(iFOV,iFOR)%landType &
                  )
             if (retrOutput%def%nEmMW > 0) &
                call putScene(ncid=ncid, &
                     emMW=retrOutput%FOR2D(iFOV,iFOR)%emMW &
                     )
             call putScene(ncid=ncid, &
                  emIR=retrOutput%FOR2D(iFOV,iFOR)%emisIR, &
                  rms=rms, &
                  chisq=chiSq, &
                  iter=retrOutput%FOR2D(iFOV,iFOR)%nIter, &
                  diagError =retrOutput%FOR2D(iFOV,iFOR)%diagError(1:nParGatm), &
                  dofs =retrOutput%FOR2D(iFOV,iFOR)%dofs, &
                  atmosclass = retrOutput%FOR2D(iFOV,iFOR)%atmosclass &
                  )
          if (trim(vCoordTyp) /= Pcoord) call putScene(ncid=ncid,pressure=retrOutput%FOR2D(iFOV,iFOR)%press)
          end do
       end do
    end if

    call closeScene(ncid)
    if (allocated(MolID)) deallocate(MolID)
    if (allocated(polMW)) deallocate(polMW)
    if (allocated(frqMW)) deallocate(frqMW)
    if (allocated(wvnIR)) deallocate(wvnIR)
    if (dbg) print *, trim(procName)//' ending ...'

  end subroutine putRetr

  subroutine putLinear(&
       ncFile, &
       stateSpaceDefin, &
       primaryRetrControl, &
       linearOutput &
       )
    !<f90Subroutine>********************************************************
    !
    ! NAME:
    !
    !   putLinear
    !
    ! PURPOSE:
    !
    !   Write data for linear inversion to a file
    !
    ! SYNTAX:
    !
    !   CALL putLinear(fileLinInv)
    !
    ! ARGUMENTS:
    !
    !   INPUTS:
    !
    !   ncFile  CHAR  Path of the output file for linear inversion coefficients
    !
    !   * OPTIONAL
    !
    ! INCLUDES:
    !
    !   None
    !
    !*******************************************************</f90Subroutine>
    character(len=*), intent(in) :: ncFile
    TYPE(StateSpaceDefin_t), intent(in) :: stateSpaceDefin
    TYPE(PrimaryRetrControl_t), intent(in) :: primaryRetrControl
    TYPE(LinearGran_t), intent(in) :: linearOutput

    !Local
    integer :: ncid
    integer :: iFOR
    integer :: iParG
    integer :: iChan
    integer :: nLevel
    integer, dimension(:), allocatable :: MolID
    real, dimension(:), allocatable :: Tsfc
    real, dimension(:), allocatable :: chiSq
    real, dimension(:,:), allocatable :: xOffG
    real, dimension(:,:,:), allocatable :: xGainG
    TYPE(StateIndex_t) :: IGR
    TYPE(StateIndex_t) :: NGR
    character (len=mxLength) :: procName

    procName = '[IOmodule::putLinear]:'
    if (dbg) print *, trim(procName)//' starting ...'

    !  Write data for linear inversion to a file

    nLevel = stateSpaceDefin%nLev
    if (.not. allocated(MolID)) allocate(MolID(size(stateSpaceDefin%MolID)))
    MolID = stateSpaceDefin%MolID
    NGR = linearOutput%def%NG
    IGR = genIndices(NGR)
    call openScene(ncid,ncFile,MolId=MolID,nLevel=nLevel,status='new',IG=IGR,NG=NGR)
    CALL writeNcdfAttr(ncid, &
         attrName='VertCoordType', &
         attr=TRIM(stateSpaceDefin%vCoordTyp))
    IF (TRIM(stateSpaceDefin%vCoordTyp) == Pcoord) &
         CALL writeNcdfAttr(ncid,attrName='standardPressureGrid',attr=stateSpaceDefin%pRef)

    call writeNcdfAttr(ncid,attrName='nchan',attr=linearOutput%def%nChan)
!!$    call writeNcdfAttr(ncid,attrName='kchan',attr=kchan(1:nchan))
    call writeNcdfAttr(ncid,attrName='chiSqThres',attr=primaryRetrControl%chiSqConv)
    call writeNcdfAttr(ncid,attrName='wavenumbers', &
         attr=stateSpaceDefin%wvnIR(1:stateSpaceDefin%nEmIR))

    if (.not. allocated(ChiSq)) allocate(ChiSq(linearOutput%nFOR))
    if (.not. allocated(Tsfc)) allocate(Tsfc(linearOutput%nFOR))
    if (.not. allocated(xOffG)) allocate(xOffG(linearOutput%def%nParG,linearOutput%nFOR))
    if (.not. allocated(xGainG)) &
         allocate(xGainG(linearOutput%def%nParG,linearOutput%def%nChan,linearOutput%nFOR))

    forall (iFOR=1:linearOutput%nFOR)
       ChiSq(iFOR) = linearOutput%FOR(iFOR)%chiSq
       Tsfc(iFOR) = linearOutput%FOR(iFOR)%Tsfc
    end forall

    forall (iFOR=1:linearOutput%nFOR,iParG=1:linearOutput%def%nParG)
       xOffG(iParG,iFOR) = linearOutput%FOR(iFOR)%xOffG(iParG)
    end forall

    forall (iFOR=1:linearOutput%nFOR,iChan=1:linearOutput%def%nChan,iParG=1:linearOutput%def%nParG)
       xGainG(iParG,iChan,iFOR) = linearOutput%FOR(iFOR)%xGainG(iChan,iParG)
    end forall

    call writeNcdfData(ncid, &
         xOffG, &
         varname='xOff', &
         varlenName=(/'nPar  ','nGroup'/), &
         varLongName='Linear inversion offset state vector')

    call writeNcdfData(ncid, &
         xGainG, &
         varname='xGain', &
         varlenName=(/'nPar  ','nChOn ','nGroup'/), &
         varLongName='Linear inversion gain matrix')

    call writeNcdfData(ncid, &
         chiSq, &
         varname='chiSqFinal', &
         varlenName=('nGroup'), &
         varLongName='Chi-Squared at last iteration', &
         varUnit='unitless')

!!$    call writeNcdfData(ncid,psfc, &
!!$         varname='pSfc', &
!!$         varlenName=('nGroup'), &
!!$         varLongName='Surface pressure', &
!!$         varUnit='mb')

    call writeNcdfData(ncid, &
         Tsfc, &
         varname='Tsfc', &
         varlenName=('nGroup'), &
         varLongName='Air temperature at the surface', &
         varUnit='K')

    call closeScene(ncid)

    if (allocated(MolID)) deallocate(MolID)
    if (allocated(ChiSq)) deallocate(ChiSq)
    if (allocated(Tsfc)) deallocate(Tsfc)
    if (allocated(xOffG)) deallocate(xOffG)
    if (allocated(xGainG)) deallocate(xGainG)

    if (dbg) print *, trim(procName)//' ending ...'

  end subroutine putLinear

  subroutine destroyMWsensorData(MWgran)
    TYPE(MWgran_t), intent(inout) :: MWgran
    !Local
    integer :: ifov
    character (len=mxLength) :: procName

    procName = ' [IOmodule::destroyMWsensorData]: '
    if (dbg) print *, trim(procName)//'starting ...'
    if (allocated(MWgran%def%pol)) deallocate(MWgran%def%pol)
    if (allocated(MWgran%def%frq)) deallocate(MWgran%def%frq)
    if (allocated(MWgran%def%chanID)) deallocate(MWgran%def%chanID)
    if (allocated(MWgran%def%planckAlpha)) deallocate(MWgran%def%planckAlpha)
    if (allocated(MWgran%def%planckBeta)) deallocate(MWgran%def%planckBeta)
    do ifov = 1, MWgran%nFOV
       if (allocated(MWgran%ob(ifov)%rad)) deallocate(MWgran%ob(ifov)%rad)
       if (allocated(MWgran%ob(ifov)%QC)) deallocate(MWgran%ob(ifov)%QC)
    end do
    if (allocated(MWgran%ob)) deallocate(MWgran%ob)
    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine destroyMWsensorData

  subroutine destroyIRsensorData(IRgran)
    TYPE(IRgran_t), intent(inout) :: IRgran
    !Local
    integer :: ifor, ifov
    character (len=mxLength) :: procName

    procName = ' [IOmodule::destroyIRsensorData]: '
    if (dbg) print *, trim(procName)//'starting ...'
    do ifor = 1, IRgran%nFOR
       do ifov = 1, IRgran%def%nFOV
          if (allocated(IRgran%FOR(ifor)%ob(ifov)%rad)) &
               deallocate(IRgran%FOR(ifor)%ob(ifov)%rad)
          if (allocated(IRgran%FOR(ifor)%ob(ifov)%QC)) &
               deallocate(IRgran%FOR(ifor)%ob(ifov)%QC)
       end do
       if (allocated(IRgran%FOR(ifor)%ob)) deallocate(IRgran%FOR(ifor)%ob)
    end do
    if (allocated(IRgran%FOR)) deallocate(IRgran%FOR)
    if (allocated(IRgran%def%wvn)) deallocate(IRgran%def%wvn)
    if (allocated(IRgran%def%chanID)) deallocate(IRgran%def%chanID)
    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine destroyIRsensorData

  subroutine destroyImgSensorData(imgGran)
    TYPE(ImgGran_t), intent(inout) :: imgGran
    !Local
    character (len=mxLength) :: procName

    procName = ' [IOmodule::destroyImgSensorData]: '
    if (dbg) print *, trim(procName)//'starting ...'
    if (allocated(imgGran%ob)) deallocate(imgGran%ob)
    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine destroyImgSensorData

  subroutine destroyancillaryData(ancillaryGran)
    TYPE(AncillaryGran_t), intent(inout) :: ancillaryGran
    integer :: ifor, ifov
    !Local
    character (len=mxLength) :: procName

    procName = ' [IOmodule::destroyAncillaryData]: '
    if (dbg) print *, trim(procName)//'starting ...'
    do ifor = 1, ancillaryGran%nFOR
       do ifov = 1, ancillaryGran%def%nFOV
          if (allocated(ancillaryGran%FOR(ifor)%data(ifov)%stateV)) &
               deallocate(ancillaryGran%FOR(ifor)%data(ifov)%stateV)
       end do
       if (allocated(ancillaryGran%FOR(ifor)%data)) &
            deallocate(ancillaryGran%FOR(ifor)%data)
       if (allocated(ancillaryGran%FOR(ifor)%emMW)) &
            deallocate(ancillaryGran%FOR(ifor)%emMW)
    end do
    if (allocated(ancillaryGran%FOR)) deallocate(ancillaryGran%FOR)
    loadAuxDone = .FALSE.
    loadExtDone = .FALSE.
    loadEmisDone = .FALSE.
    if (dbg) print *, trim(procName)//'ending ...'
  end subroutine destroyancillaryData

END MODULE IOmodule
