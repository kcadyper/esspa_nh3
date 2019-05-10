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

MODULE CloudModule

! <f90Module>***********************************************************
!
! NAME:
!
!   CloudModule
!
! PURPOSE:
!
!   Convert between cloud physical properties and optical properties.
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>
  use params_module, ONLY: mxlay

  USE CloudParameters
  USE CloudDataStruct
  USE CloudTab

  IMPLICIT none
  PRIVATE

  PUBLIC :: initCloud, setCloudOptLevel, setCloudOptSlab, destroyCloud
  PUBLIC :: checkCloudModelType
  PUBLIC :: getCloudModelType
  public :: setCloudOptLayer
  public :: printCloudProfile
  public :: paramIndex

  real, public, parameter            :: minCloudAmount = epsilon(1.0) ! 1d-11
  INTEGER,PARAMETER,PUBLIC           :: iCLW=1,iIce=2,iRain=3,iSnow=4,iGraupel=5

  INTEGER,PARAMETER                  :: nHydDflt=5
  INTEGER,PARAMETER                  :: iCTPice=0,iThkIce=1,iAmtIce=2,iDeffIce=3
  INTEGER,PARAMETER                  :: iCTPliq=0,iThkLiq=1,iAmtLiq=2,iDeffLiq=3
  INTEGER,PARAMETER,PUBLIC           :: typeSlab=0,typeLevel=1,typeLayer=2
  INTEGER,PARAMETER,PUBLIC           :: typeCloudDefault=-1

  INTEGER,PARAMETER                  :: topIndx(2) = (/iCTPliq,iCTPice/)
  INTEGER,PARAMETER                  :: thkIndx(2) = (/iThkLiq,iThkIce/)
  INTEGER,PARAMETER                  :: amtIndx(2) = (/iAmtLiq,iAmtIce/)
  INTEGER,PARAMETER                  :: sizIndx(2) = (/iDeffLiq,iDeffIce/)

  INTEGER,PARAMETER                  :: mxLev=120
  INTEGER,PARAMETER                  :: nOptical=mxOpticalProperty
  REAL,DIMENSION(:),ALLOCATABLE      :: kabs,kscat,asymi,itscat
  INTEGER,DIMENSION(:),ALLOCATABLE   :: hydIndx

  INTEGER,DIMENSION(:,:),ALLOCATABLE :: indxArr
  INTEGER                            :: nSpcPts,nHydLoc
  INTEGER                            :: npProperty
  REAL,DIMENSION(:),ALLOCATABLE      :: intAmt
  CHARACTER(LEN=20), DIMENSION(:), ALLOCATABLE :: mydataName
  character(len=20), dimension(:), allocatable :: hydName
  TYPE(CloudTable), DIMENSION(:), allocatable :: cldTbl
  TYPE(CloudOpticalProperty), dimension(:), allocatable :: cldOpt
  TYPE(CloudPhysicalSet)             :: cldSet
  INTEGER                            :: iSizP, iSizS, iTemp, iDens, iMeltF
  INTEGER                            :: ikabs, ikscat, igasym
  INTEGER                            :: moduleCloudModelType ! 'level', 'layer' or slab
  !------------------------------------------------
  !  old convention
  !  CHARACTER(LEN=7),PARAMETER, public :: moduleCloudModelList(3) = (/'profile',
  !                                                         'proflay', 'layer  '/)
  !
  !------------------------------------------------
  integer,parameter, public         :: moduleCloudModelList(3) = (/typeSlab, &
                                                                   typeLevel, &
                                                                   typeLayer/)

  character(len=*), parameter        :: modName='CloudModule::'

CONTAINS
  function paramIndex(paramName)
    character(len=*) paramName
    integer paramIndex
    ! local
    character(len=len(paramName)) str
    integer j, k
    character(len=*), parameter        :: progName=modName//'paramIndex'

    ! convert to upper case
    do j=1,len(paramName)
      k=iachar(paramName(j:j))
      if (k>=iachar('a') .and. k<=iachar('z')) then
        str(j:j)=achar(k-32)
      else
        str(j:j)=paramName(j:j)
      end if
    end do

    select case(str)
    case('SIZP')
      paramIndex=iSizP
    case('SIZS')
      paramIndex=iSizS
    case('TEMP')
      paramIndex=iTemp
    case('IDENS')
      paramIndex=iDens
    case('IMELTF')
      paramIndex=iMeltF
    case default
      print *, progName, 'error: unknown parameter name: ', paramName, str
      call exit(1)
    end select
  end function paramIndex

  SUBROUTINE initCloud(ctab_file,cldScene,spectPtsGrd,cldTyp,nHyd)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   initCloud
!
! PURPOSE:
!
!   Initialize conversion between cloud physical properties and optical
!   properties, specifying which hydrometeors to treat and the type of profile
!   model.
!
! SYNTAX:
!
!   CALL initCloud(ctab_file, dataType, spectPtsGrd, iCldLiq,
!      iCldIce, cldTyp, nHyd)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   ctab_file    CHAR     Cloud optical properties LUT file path
!   iCldLiq*     INTEGER  Starting index for liquid cloud in state vector
!   iCldIce*     INTEGER  Starting index for ice cloud in state vector
!   cldTyp*      INTEGER  Type of cloud model used
!   nHyd*        INTEGER  Number of hydrometeor types
!
!   INPUTS/OUTPUTS:
!
!   spectPtsGrd  REAL     Spectral grid
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>


    USE OSSPracticalConstant, ONLY: MISSING_CHAR

    !---Input/Output variables:
    CHARACTER(LEN=*),              INTENT(IN)   :: ctab_file
    TYPE(CloudScene), DIMENSION(:),INTENT(IN)   :: cldScene
    REAL, DIMENSION(:), POINTER,   INTENT(INOUT):: spectPtsGrd
    INTEGER,                       INTENT(IN)   :: cldTyp
    INTEGER,            OPTIONAL,  INTENT(IN)   :: nHyd

    !---Local variables:
    CHARACTER(LEN=20),DIMENSION(:),ALLOCATABLE  :: dataName
    INTEGER                                     :: dSize
    INTEGER                                     :: i,dim1,ict
    INTEGER,DIMENSION(:),ALLOCATABLE            :: indxArr1d
    integer                                     :: nNames
    integer                                     :: iHydr
    integer                                     :: iop
    integer                                     :: ipp
    character(len=*), parameter                 :: procName = " [CloudModule::initCloud]: "
    character(len=*), parameter                 :: errHeader='Err['//procName//']: '

    moduleCloudModelType = cldTyp

    if(present(nHyd))then
       if(nHyd > nHydDflt)then
          print*,errHeader,'Input nHyd exceeds the default (maximum) nHyd value: '
          print*,nHyd,nHydDflt
          call exit(1)
       endif
       nHydLoc=nHyd
    else
       nHydLoc=nHydDflt
    endif
! supported models
! profile - define on levels
! proflay - paraneters defined on layers
! layer   - 4 parameter slab only 2 types of clouds can be used

    if( .not. any(cldTyp == moduleCloudModelList) ) then
        print*, errHeader, 'Unknown mode: ', cldTyp
        print *, 'supported modes:'
        do i=1,size(moduleCloudModelList)
            print *, moduleCloudModelList(i)
        end do
        call exit(1)
    end if

    !TODO user has no control. It is always preset #
    if (cldTyp == typeSlab) then
      nHydLoc=nHydNonPcp
      allocate(hydIndx(nHydLoc))
    end if

    ict=0
    do i=1,size(cldScene)
       if (cldScene(i)%name(1:1) /= MISSING_CHAR) ict=ict+1
    enddo
    if(ict < nHydLoc)then
       print*, errHeader,'nHydLoc parameter greater than number of requested '
       print*,' hydrometeors in dataType character vector.'
       print*,'nHydLoc = ',nHydLoc
       print*,'# of requested hydrometeors: ',ict
       call exit(1)
    endif

    dSize = nHydLoc*nOptical
    allocate(indxArr1d(dSize))
    indxArr1d = (/(i,i=1,dSize)/)
    dim1 = dsize/nOptical
    allocate(indxArr(dim1,nOptical))
    indxArr = reshape(indxArr1d,(/dim1,nOptical/),order=(/2,1/))

    allocate(hydName(nHydLoc))
    do iHydr=1,nHydLoc
       hydName(iHydr) = trim(cldScene(iHydr)%name)
    end do


    ! indxArr logic assumes all three optical properties are requested
    !   for each hydrometeor type; thus, 'nOptical' evenly divides 'dSize'

    !--Prepare dataName by combining hydrometeor and optical property ID strings:
    call buildCldDataNames(cldScene(1:nHydLoc)%name,nNames)
    allocate(dataName(nNames))

    call getCldDataNames(dataName)

    CALL initCloudTab(ctab_file,dataName,cldScene,outBndsExit_in=.false.,&
                      nSpcPts=nSpcPts, spectPtsGrd=spectPtsGrd)

    ! Prepare the optical properties. Only the properties interpolated
    ! from all physical properties to the spectral points for each
    ! hydrometeor is needed.
    ! TODO possible complication that this datatype allows that different
    ! cloud types may have a different number of physical properties

    npProperty = cldScene(1)%nProperty
    if (allocated(cldOpt)) deallocate(cldOpt)
    allocate(cldOpt(nOptical))
    do iop = 1, nOptical
       cldOpt(iop)%name = trim(knownOptProp(iop))
       cldOpt(iop)%units = trim(knownOptUnits(iop))
       select case (trim(cldOpt(iop)%name))
         case ('kabs')
            ikabs = iop
         case ('kscat')
            ikscat = iop
         case ('gasym')
            igasym = iop
         case default
            print *, errHeader,' Unrecognized optical table:'&
                   ,trim(cldOpt(iop)%name)
            call exit(99)
       end select
       allocate(cldOpt(iop)%values(nSpcPts))
       allocate(cldOpt(iop)%deriv(npProperty,nSpcPts))
    end do

    do ipp = 1, npProperty
      select case(trim(cldScene(1)%physDescr(ipp)%name))
        case ('SizP')
          iSizP = ipp
        case ('SizS')
          iSizS = ipp
        case ('Temp')
          iTemp = ipp
        case ('Dens')
          iDens = ipp
        case ('MeltF')
          iMeltF = ipp
        case default
          print *,  'Warning: ', procName
          print *, trim(cldScene(1)%physDescr(ipp)%name), ' property NOT found'
      end select
    end do

    allocate(kabs(nSpcPts),kscat(nSpcPts),asymi(nSpcPts),itscat(nSpcPts))
    allocate(intAmt(nHydLoc))
    cldSet%nProperty = npProperty
    allocate(cldSet%value(npProperty))
    deallocate(indxArr1d,dataName)

    return

  END SUBROUTINE initCloud

  !===================================

  SUBROUTINE buildCldDataNames(dataType,nNames)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   buildCldDataNames
!
! PURPOSE:
!
!   Prepare names of data needed from optical property tables by combining
!   hydrometeor and optical property ID strings.
!
! SYNTAX:
!
!   CALL buildCldDataNames(dataType, nNames)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   dataType  CHAR     List of hydrometeor types
!
!   INPUTS/OUTPUTS:
!
!   nNames    INTEGER  Number of hydrometeor and optical property types
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !--This is to give access to information that allows top-level program to call
    !--loadCloudTab.  See comment in CloudTab::cloudTabInit.

    !---Input/Output variables:
    character(len=*), dimension(:), intent(in)    :: dataType
    integer,                        intent(inout) :: nNames

    !---Local variables:
    integer                                       :: i,j,ict,mynHyd
    character(len=1),               parameter     :: spc='_'

    mynHyd=size(dataType)
    nNames=mynHyd*nOptical
    if (allocated(mydataName)) deallocate(mydataName)
    allocate (mydataName(nNames))

    !--Prepare dataName by combining hydrometeor and optical property ID strings:
    ict=1
    do i=1,mynHyd
       do j=1,nOptical
          mydataName(ict) = trim(dataType(i)) // spc // trim(knownOptProp(j))
          ict=ict+1
       enddo
    enddo

    return

  END SUBROUTINE buildCldDataNames

  ! ===================================

  SUBROUTINE getCldDataNames(dataName)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   getCldDataNames
!
! PURPOSE:
!
!   Get names of data needed from optical property tables, from global memory.
!   This is to give access to information that allows top-level program to call
!   loadCloudTab.
!
! SYNTAX:
!
!   CALL getCldDataNames(dataName)
!
! ARGUMENTS:
!
!   INPUTS/OUTPUTS:
!
!   dataName  CHAR  String combining hydrometeor and optical property ID
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !--This is to give access to information that allows top-level program to call
    !--loadCloudTab.  See comment in CloudTab::cloudTabInit.

    !---Input/Output variables:
    character(len=*), dimension(:), intent(inout) :: dataName

    if (size(dataName) /= size(mydataName)) then
       print *,'err[CloudModule::getCldDataNames]: ', &
            'vector length mismatch; size(dataName) /= size(mydataName):', &
            size(dataName) /= size(mydataName)
       call exit(1)
    endif
    dataName = mydataName

    return

  END SUBROUTINE getCldDataNames

  ! ===================================
  SUBROUTINE setCloudOptLayer(cldScene,nSurf,tavl,tabs,tscat,asym,&
                            tabsJac,tscatJac,asymJac,jacUp,jacLw)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   setCloudOptLayer
!
! PURPOSE:
!
!   Convert from cloud physical properties to optical properties, with input as
!   a profile defined on layers.
!
! SYNTAX:
!
!   setCloudOptLayer(xG,cldScene,nSurf,pref,psfc,ih2o,tabs,tscat,asym, &
!               tabsJacUp,tscatJacUp,asymJacUp,tabsJacLw,tscatJacLw,asymJacLw)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   xG     REAL     State vector in geophysical space
!   cldScene        derived type (CloudScene)
!   nSurf  INTEGER  Index of first atmospheric level at or below surface
!   pref   REAL     Pressure on atmospheric grid levels
!   psfc   REAL     Surface pressure
!   ih2o   INTEGER  Index for water vapor in molecules part of state vector
!   tavl   real     layer average temperature
!
!   INPUTS/OUTPUTS:
!
!   tabs   REAL     In-layer cloud contribution to absorption optical depth
!   tscat  REAL     In-layer cloud contribution to scattering optical depth
!   asym   REAL     Asymmetry parameter
!
!   optional
!   tabsJacXX   REAL     cloud layer absorption optical depth Jacobian
!   tscatJacXX  REAL     cloud layer scattering optical depth Jacobian
!   asymJacXX   REAL     cloud layer asymmetry parameter Jacobian
!   where XX is 'Up', 'Lw' denoting effect of upper or lower level
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    USE OSSPhysicalConstant, ONLY : rd,g0,ratioMwt
    !---Input variables:
    INTEGER,              INTENT(IN)    :: nSurf
    REAL, DIMENSION(:),   INTENT(IN)    :: tavl
    TYPE(CloudScene), DIMENSION(:), INTENT(INOUT) :: cldScene
    !---Output variables:
    REAL, DIMENSION(:,:), INTENT(INOUT) :: tabs,tscat,asym

    REAL, DIMENSION(:,:,0:,:), INTENT(INOUT), optional :: tabsJac,tscatJac,asymJac   ! layer Jacobians 0 - amount, 1, 2, ... properties
    REAL, DIMENSION(:,:,:), INTENT(INOUT), optional :: jacUp,jacLw                   ! convert layer to level
    !---Local variables:
    REAL, DIMENSION(mxLev,nHydLoc)     :: uAmt
    INTEGER                            :: iLyr,iHydr,iSpc,n2,ipp
    REAL,PARAMETER                     :: amtScl=1e3  ! from kg m-2 to g m-2
    INTEGER                            :: nLevel, nHydr

    logical                            :: JacobianCalculation
    real                               :: deltaG
    integer                            :: cloudPropNum

    if ( present(tabsJac) .and. present(tscatJac) .and. present(asymJac) ) then
        JacobianCalculation = .TRUE.
        tabsJac=0.
        tscatJac=0.
        asymJac=0.
    else
        JacobianCalculation = .FALSE.
    end if

    uAmt = 0
    n2 = nSurf-1
    cloudPropNum = cldScene(1)%nProperty


    nHydr = size(cldScene)
    if (nHydr > nHydLoc) then
      print *, "[error]: max number of hydrometeors exceeded"
      call exit(67)
    end if
    do iHydr = 1, nHydr
       nLevel = cldScene(iHydr)%nLevel
       if (nLevel > mxLev) then
          print *, "[error]: max number of hydrometeors exceeded"
          call exit(66)
       end if
    end do
    !---Compute physical property from level to layer
    do iHydr = 1, nHydLoc
        do iLyr=1,cldScene(iHydr)%nLayer
            cldScene(iHydr)%physProp(iLyr)%value(iTemp) = tavl(iLyr)
        end do
    end do

    itscat = 0.
    tabs=0.0
    tscat=0.0
    asym=0.0
    do iLyr=1,n2
       do iHydr=1,nHydLoc
           uAmt(iLyr,iHydr) = amtScl * cldScene(iHydr)%hydrAmt(iLyr)
       end do
    end do
    do iLyr=1,n2
       do iHydr=1,nHydLoc
          !---Get tabulated optical properties, convert to optical depths, and sum
          !    over hydrometeor classes:
         if(cldScene(iHydr)%hydrAmt(iLyr) < minCloudAmount) cycle
         call getCloudTab(cldScene(iHydr)%name,cldScene(iHydr)%physProp(iLyr),cldOpt)

         ! sum(a_n w_n)
         tabs(iLyr,:)  = tabs(iLyr,:)  + cldOpt(ikabs)%values * uAmt(iLyr,iHydr)
         itscat = cldOpt(ikscat)%values * uAmt(iLyr,iHydr)
         tscat(iLyr,:) = tscat(iLyr,:) + itscat(:)
         asym(iLyr,:)  = asym(iLyr,:)  + itscat(:) * cldOpt(igasym)%values

          if (JacobianCalculation) then
              ! a_n
              tabsJac(iLyr,:,0,iHydr) = cldOpt(ikabs)%values*amtScl
              ! b_n
              tscatJac(iLyr,:,0,iHydr) = cldOpt(ikscat)%values*amtScl
              ! g_n
              asymJac(iLyr,:,0,iHydr)= cldOpt(igasym)%values
              do ipp=1, cldScene(iHydr)%nProperty                       ! w_n da_n/dp
                 ! w_n da_n/dp
                 tabsJac(iLyr,:,ipp,iHydr)=cldOpt(ikabs)%deriv(ipp,:) * uAmt(iLyr,iHydr)
                 ! w_n db_n/dp
                 tscatJac(iLyr,:,ipp,iHydr)=cldOpt(ikscat)%deriv(ipp,:)*uAmt(iLyr,iHydr)
                 ! w_n*b_n*dg_n/dp
                 asymJac(iLyr,:,ipp,iHydr) = itscat(:) * cldOpt(igasym)%deriv(ipp,:)
              end do
          end if
      end do
!
!       !---The following loop could be replaced with a WHERE construct, although the
!       !     current version of the SGI compiler does not have the proper libraries
!       !     to support WHERE in DEBUG mode.
       do iSpc=1,nSpcPts
          if(tscat(iLyr,iSpc) > minCloudAmount) then
             ! g = sum(s_n g_n)/ sum(s_n)
             asym(iLyr,iSpc) = asym(iLyr,iSpc) / tscat(iLyr,iSpc)
             if (JacobianCalculation) then
                ! dg/dw_n = (g_n - g)*b_n/ sum(s_n)
                ! dg/dp = (w_n (g_n - g) *db_n/dp + s_n dg_n/dp)/ sum(s_n)

                 do iHydr=1,nHydLoc
                    deltaG = asymJac(iLyr,iSpc,0,iHydr) -  asym(iLyr,iSpc)
                    asymJac(iLyr,iSpc,0,iHydr)=  deltaG *tscatJac(iLyr,iSpc,0,iHydr)

                    do ipp=1, cloudPropNum
                       asymJac(iLyr,iSpc,ipp,iHydr) = &
                          tscatJac(iLyr,iSpc,ipp,iHydr) * deltaG + asymJac(iLyr,iSpc,ipp,iHydr)
                    end do

                 end do
                 asymJac(iLyr,iSpc,:,:)=asymJac(iLyr,iSpc,:,:) / tscat(iLyr,iSpc)
             end if
          else
             asym(iLyr,iSpc) = 0.
             if (JacobianCalculation) then
                 do iHydr=1,nHydLoc
                    asymJac(iLyr,iSpc,:,iHydr)=  0.0
                 end do
             end if
          endif
       end do
    end do
    RETURN
  END SUBROUTINE setCloudOptLayer
 ! ===================================

  SUBROUTINE setCloudOptLevel(xG,cldScene,nSurf,pref,psfc,ih2o,tavl,tabs,tscat,asym,&
                            tabsJac,tscatJac,asymJac,JacUp,JacLw)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   setCloudOptLevel
!
! PURPOSE:
!
!   Convert from cloud physical properties to optical properties, with input as
!   a profile.
!
! SYNTAX:
!
!   setCloudOptLevel(xG,cldScene,nSurf,pref,psfc,ih2o,tabs,tscat,asym, &
!               tabsJacUp,tscatJacUp,asymJacUp,tabsJacLw,tscatJacLw,asymJacLw)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   xG     REAL     State vector in geophysical space
!   cldScene        derived type (CloudScene)
!   nSurf  INTEGER  Index of first atmospheric level at or below surface
!   pref   REAL     Pressure on atmospheric grid levels
!   psfc   REAL     Surface pressure
!   ih2o   INTEGER  Index for water vapor in molecules part of state vector
!   tavl   real     layer average temperature
!
!   INPUTS/OUTPUTS:
!
!   tabs   REAL     In-layer cloud contribution to absorption optical depth
!   tscat  REAL     In-layer cloud contribution to scattering optical depth
!   asym   REAL     Asymmetry parameter
!
!   optional
!   tabsJacXX   REAL     cloud layer absorption optical depth Jacobian
!   tscatJacXX  REAL     cloud layer scattering optical depth Jacobian
!   asymJacXX   REAL     cloud layer asymmetry parameter Jacobian
!   where XX is 'Up', 'Lw' denoting effect of upper or lower level
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    USE OSSPhysicalConstant, ONLY : rd,g0,ratioMwt
    !---Input variables:
    INTEGER,              INTENT(IN)    :: nSurf,ih2o
    REAL, DIMENSION(:),   INTENT(IN)    :: xG,pref
    REAL, DIMENSION(:),   INTENT(IN)    :: tavl
    TYPE(CloudScene), DIMENSION(:), INTENT(INOUT) :: cldScene
    REAL,                 INTENT(IN)    :: psfc
    !---Output variables:
    REAL, DIMENSION(:,:), INTENT(INOUT) :: tabs,tscat,asym

    ! layer Jacobians 0 - amount, 1, 2, ... properties
    REAL, DIMENSION(:,:,0:,:), INTENT(INOUT), optional :: tabsJac,tscatJac,asymJac

    ! convert layer to level
    REAL, DIMENSION(:,0:,:), INTENT(INOUT), optional :: jacUp,jacLw
    !---Local variables:
    REAL, DIMENSION(mxLev,nHydLoc)     :: uAmt,mcon
    real                               :: dAmtU,dAmtL
    REAL, DIMENSION(nHydLoc)           :: mconSfc
    INTEGER                            :: iLyr,iLev,iHydr,iSpc,nLayer, ipp
    REAL                               :: scal
    REAL,PARAMETER                     :: amtScl=1e5
    INTEGER                            :: nLevel, nHydr, cloudPropNum

    logical                            :: JacobianCalculation
    real                               :: deltaG
    real                               :: denum
    character(len=*), parameter        :: procName = " [CloudModule::setCloudOptLevel]: "

    if ( present(tabsJac) .and. present(tabsJac) .and. present(tabsJac)) then
        JacobianCalculation = .TRUE.
    else
        JacobianCalculation = .FALSE.
    end if

    uAmt = 0
    nLayer = nSurf-1
    cloudPropNum = cldScene(1)%nProperty

    nHydr = size(cldScene)
    if (nHydr > nHydLoc) then
        print *, procName, " max number of hydrometeors exceeded"
        call exit(67)
    end if
    do iHydr = 1, nHydr
       nLevel = cldScene(iHydr)%nLevel
       if (nLevel > mxLev) then
         print *, procName, " max number of levels for hydrometeor ", &
                  iHydr, cldScene(iHydr)%name, " exceeded"
         call exit(66)
       end if
    end do

    !---Compute mass concentration of each hydrometeor class on the fixed
    !     pressure levels:
    do iLev=1,nSurf
       denum = 1.0 + xG(ih2o+iLev-1)
       do iHydr = 1, nHydLoc
          denum = denum + cldScene(iHydr)%hydrAmt(iLev)
       end do
       do iHydr=1,nHydLoc
         mcon(iLev,iHydr) = cldScene(iHydr)%hydrAmt(iLev) / denum
         if (JacobianCalculation) then
             jacUp(iLev,0,iHydr) = (1.0 - mcon(iLev,iHydr)) / denum
             jacLw(iLev,0,iHydr) = jacUp(iLev,0,iHydr)
         end if
       enddo
    enddo
    !---Compute hydrometeor layer amount:
    ! -- near surface layer will be treated separately

    scal = amtScl/g0 ! convert p(mb) to g/m^2 ; (1mb = 100 kg/(m s^2))
    if (JacobianCalculation) jacLw(1,0,:) = 0.0
    do iLyr=1,nLayer-1
       do iHydr = 1,nHydLoc
          if (any(mcon(iLyr:iLyr+1,iHydr) > minCloudAmount)) then
             call lpsum(pref(iLyr),pref(iLyr+1),mcon(iLyr,iHydr), &
                  mcon(iLyr+1,iHydr),scal, &
                  uAmt(iLyr,iHydr),dAmtU,dAmtL)
             if (JacobianCalculation) then
                  jacLw(iLyr+1,0,iHydr) = jacLw(iLyr+1,0,iHydr)*dAmtL
                  jacUp(iLyr,0,iHydr)   = jacUp(iLyr,0,iHydr)*dAmtU
             end if
          else
             uAmt(iLyr,iHydr)=0.
             if (JacobianCalculation) then
               jacLw(iLyr+1,0,iHydr) = 0.0
               jacUp(iLyr,0,iHydr) = 0.0
             end if
          endif
       enddo
    enddo
    !---Surface layer:
    do iHydr=1,nHydLoc
       if(max(mcon(nSurf,iHydr),mcon(nSurf-1,iHydr))>epsilon(1.))then
          call lint(mcon(:,iHydr),pref,nSurf,psfc,mconSfc(iHydr))
          call lpsum(pref(nLayer),psfc,mcon(nLayer,iHydr),mconSfc(iHydr),scal,uAmt(nLayer,iHydr), &
               dAmtU,dAmtL)

         if (JacobianCalculation) then
           jacLw(nSurf ,0,iHydr) = jacLw(nSurf ,0,iHydr)*dAmtL
           jacUp(nLayer,0,iHydr) = jacUp(nLayer,0,iHydr)*dAmtU
         end if
       else
          uAmt(nLayer,iHydr)=0.
          if (JacobianCalculation) then
            jacLw(nSurf, 0,iHydr) = 0.0
            jacUp(nLayer,0,iHydr) = 0.0
          end if
       end if
    end do

    !---Interpolate the physical properties from level to layer
    !---Since the layer cloud amount is non zero only if any of
    !---level cloud amounts is
    do iHydr = 1, nHydLoc
      ! first clean up
      if (allocated(cldScene(iHydr)%physPropLayer) ) then
          do iLev=1,cldScene(iHydr)%nLayer
              if ( allocated(cldScene(iHydr)%physPropLayer(iLev)%value) ) &
                         deallocate(cldScene(iHydr)%physPropLayer(iLev)%value)
          end do
          deallocate(cldScene(iHydr)%physPropLayer)
      end if
      allocate(cldScene(iHydr)%physPropLayer(cldScene(iHydr)%nLayer))
      do iLev=1,nLayer
        cldScene(iHydr)%physPropLayer(iLev)%nProperty = cldScene(iHydr)%nProperty

        allocate(cldScene(iHydr)%physPropLayer(iLev)%value(cldScene(iHydr)%nProperty))

        do ipp=1,cloudPropNum
           if (ipp == iTemp) cycle
            cldScene(iHydr)%physPropLayer(iLev)%value(ipp) = &
               (cldScene(iHydr)%physProp(iLev)%value(ipp) + &
                cldScene(iHydr)%physProp(iLev+1)%value(ipp)) / 2.0
        end do

        ! special treatment for layer temperature
        cldScene(iHydr)%physProp(iLev)%value(iTemp) = tavl(iLev)
        cldScene(iHydr)%physPropLayer(iLev)%value(iTemp) = tavl(iLev)

        if (JacobianCalculation) then
          jacUp(iLev,  1:cloudPropNum,iHydr) = 0.5
          jacLw(iLev+1,1:cloudPropNum,iHydr) = 0.5
        end if
      end do
      if (JacobianCalculation) then
        jacLw(1,    1:cloudPropNum,iHydr) = 0.
        jacUp(nSurf,1:cloudPropNum,iHydr) = 0.
      end if
    end do

    itscat = 0.
    tabs=0.
    tscat=0.
    asym=0.
    if (JacobianCalculation) then
      tabsJac=0.
      asymJac=0.
      tscatJac=0.
    end if

    do iLyr=1,nLayer
      do iHydr=1,nHydLoc
        !---Get tabulated optical properties, convert to optical depths, and sum
        !    over hydrometeor classes:
        if(uAmt(iLyr,iHydr) < minCloudAmount) cycle

        call getCloudTab(cldScene(iHydr)%name,cldScene(iHydr)%physPropLayer(iLyr),cldOpt)

        ! sum(a_n w_n)
        tabs(iLyr,:)  = tabs(iLyr,:)  + cldOpt(ikabs)%values * uAmt(iLyr,iHydr)
        itscat = cldOpt(ikscat)%values * uAmt(iLyr,iHydr)        ! s_n=b_n w_n
        tscat(iLyr,:) = tscat(iLyr,:) + itscat(:)                ! sum(s_n)
         asym(iLyr,:)  = asym(iLyr,:)  + itscat(:) * cldOpt(igasym)%values    ! sum(s_n g_n)

        if (JacobianCalculation) then
          tabsJac(iLyr,:,0,iHydr) =cldOpt(ikabs) %values            ! a_n
          tscatJac(iLyr,:,0,iHydr)=cldOpt(ikscat)%values            ! b_n
          ! g_n
          asymJac(iLyr,:,0,iHydr)= cldOpt(igasym)%values ! * cldOpt(ikscat)%values
          do ipp=1, cldScene(iHydr)%nProperty
            ! w_n*da_n/dp
            tabsJac(iLyr,:,ipp,iHydr) =cldOpt(ikabs) %deriv(ipp,:)*uAmt(iLyr,iHydr)
            ! w_n*db_n/dp
            tscatJac(iLyr,:,ipp,iHydr)=cldOpt(ikscat)%deriv(ipp,:)*uAmt(iLyr,iHydr)
            ! w_n*b_n*dg_n/dp
            asymJac(iLyr,:,ipp,iHydr) = itscat(:) * cldOpt(igasym)%deriv(ipp,:)
          end do
        end if
      end do

      do iSpc=1,nSpcPts
        if(tscat(iLyr,iSpc) > minCloudAmount) then
          asym(iLyr,iSpc) = asym(iLyr,iSpc) / tscat(iLyr,iSpc)   ! g = sum(s_n g_n)/ sum(s_n)
          if (JacobianCalculation) then
            ! dg/dw_n = (g_n - g)*b_n/ sum(s_n)
            do iHydr=1,nHydLoc
              deltaG = asymJac(iLyr,iSpc,0,iHydr) -  asym(iLyr,iSpc)
              asymJac(iLyr,iSpc,0,iHydr)=  deltaG *tscatJac(iLyr,iSpc,0,iHydr)

              do ipp=1, cldScene(iHydr)%nProperty
                asymJac(iLyr,iSpc,ipp,iHydr) = &
                    tscatJac(iLyr,iSpc,ipp,iHydr) * deltaG + asymJac(iLyr,iSpc,ipp,iHydr)
              end do
            end do
            asymJac(iLyr,iSpc,:,:)=asymJac(iLyr,iSpc,:,:) / tscat(iLyr,iSpc)
          end if
        else
          asym(iLyr,iSpc) = 0.
          if (JacobianCalculation) asymJac(iLyr,iSpc,:,:)=0.
        end if
      end do
    end do
    return
  END SUBROUTINE setCloudOptLevel

!==================================================================================================
  SUBROUTINE setCloudOptSlab(xG,tavl,nSurf,pref,psfc,iCldLiq,iCldIce,tabs,tscat,&
                             asym,cloudProfile,tabsJac,tscatJac,asymJac)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   setCloudOptSlab
!
! PURPOSE:
!
!   Convert from cloud physical properties to optical properties, with input
!   from a layer model. If necessary it adds up to 2 additional layer
!
! SYNTAX:
!
!   CALL setCloudOptSlab(xG, cldScene, nSurf, pref, psfc, iCldLiq,
!      iCldIce, tabs, tscat, asym)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   xG       REAL     State vector in geophysical space
!   tavl     REAL     layer-averaged temperature profile
!   nSurf    INTEGER  Index of first atmospheric level at or below surface
!   pref     REAL     Pressure on atmospheric grid levels
!   psfc     REAL     Surface pressure
!   iCldLiq  INTEGER  Starting index for liquid cloud in state vector
!   iCldIce  INTEGER  Starting index for ice cloud in state vector
!
!   INPUTS/OUTPUTS:
!
!   tabs     REAL     In-layer cloud contribution to absorption optical depth
!   tscat    REAL     In-layer cloud contribution to scattering optical depth
!   asym     REAL     Asymmetry parameter
!
!   * OPTIONAL
!   tabsJac   REAL     cloud layer absorption optical depth Jacobian
!   tscatJac  REAL     cloud layer scattering optical depth Jacobian
!   asymJac   REAL     cloud layer asymmetry parameter Jacobian
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !---Input variables:
    INTEGER,              INTENT(IN)    :: nSurf,iCldLiq,iCldIce
    REAL, DIMENSION(:),   INTENT(IN)    :: xG,pref
    REAL, DIMENSION(:),   INTENT(IN)    :: tavl
    REAL,                 INTENT(IN)    :: psfc
    !
    type(CloudProfileType), intent(inout)    :: cloudProfile
    !---Output variables:
    REAL, DIMENSION(:,:), INTENT(INOUT) :: tabs,tscat,asym
    REAL, DIMENSION(:,:,0:,:), INTENT(INOUT), optional :: tabsJac,tscatJac,asymJac

    !---Local variables:
    REAL, DIMENSION(mxLev,nHydLoc)     :: cFrac
    INTEGER,DIMENSION(mxLev)           :: ltyp
    INTEGER                            :: iLyr,iHydr,iSpc,nLayer
    INTEGER                            :: topID,thkID,amtID
    INTEGER                            :: j,j1,j2, k
    integer                            :: indexLayer(2, nHydLoc) ! contain layer index for top and bottom layers containing the cloud
    real                               :: cldThickness(nHydLoc)
    REAL                               :: cldTop,cldBot,cldThk
    REAL,PARAMETER,DIMENSION(nHydNonPcp) :: sclPath=(/(1000.,j=1,nHydNonPcp)/)

    real                               :: deltaG, ratio
    logical                            :: JacobianCalculation
    integer                            :: curLevNum, idx

    if (nHydLoc /= 2) then
         print *,'Cloud slab model must contain only 2 component: liquid and ice'
         print *, 'number of cloud models', nHydLoc
         call exit(1)
    end if

    if ( present(tabsJac) .and. present(tscatJac) .and. present(asymJac) ) then
        JacobianCalculation = .TRUE.
        tabsJac=0.
        tscatJac=0.
        asymJac=0.
    else
        JacobianCalculation = .FALSE.
    end if

    tabs  = 0.
    tscat = 0.
    asym  = 0.
    nLayer = nSurf-1

    ltyp(:)    = 0
    cFrac(:,:) = 0.
    intAmt(:)  = 0.
    hydIndx = (/iCldLiq,iCldIce/)
    !-- Input: several-parameter cloud model
    !  - Determine fraction of cloud falling within each layer (for each hydro type).

    ! initialize
    cloudProfile%presLevel(1:nSurf) = pref(1:nSurf)
    cloudProfile%layerIndex(1:nLayer)   = (/ (k, k=1,nLayer) /)
    cloudProfile%toUse = .false.
    cloudProfile%levNum = nSurf

    cloudProfile%cloudTop = 0
    cloudProfile%cloudBot = 0
    cloudProfile%tau=0.0
    cloudProfile%ratio=1.0
    cloudProfile%cloudPresent=.false.

    curLevNum = nSurf
    indexLayer=0
    do iHydr=1,nHydLoc
       amtID = hydIndx(iHydr)+amtIndx(iHydr)
       if (xG(amtID) > minCloudAmount) then
          topID = hydIndx(iHydr)+topIndx(iHydr)
          thkID = hydIndx(iHydr)+thkIndx(iHydr)
          intAmt(iHydr) = xG(amtID) * sclPath(iHydr)

          cloudProfile%cloudPresent(iHydr)=.true.

          cldTop = xG(topID)
          cldThk = xG(thkID)
          if (cldThk <= epsilon(1.)) then
             print *, 'cloud module does not handle clouds of pressure thickness less then ', epsilon(1.)
             print *, 'cloud pressure thickness for model ', iHydr, ' is ', cldThk
             call exit(1)
          end if

          cldBot = cldThk + cldTop

          if(levelTheSame(cldTop, cldBot)) then
             print *, 'cloud module considers the pressure levels for cloud top and cloud bottom the same'
             print *, 'cloud pressure top    ', cldTop
             print *, 'cloud pressure bottom ', cldBot
             call exit(1)
          end if

          if (cldBot > psfc  + 1.0) then
             print *,'Cloud Model cannot treat clouds with bottom boundary below surface'
             print *, 'surface pressure', psfc
             print *, 'cloud bottom pressure', cldBot
             call exit(1)
          end if

          cldThickness(iHydr) = cldThk
          call insertLevel(cldTop, cloudProfile)
          call insertLevel(cldBot, cloudProfile)
       end if
    end do
    do j=1, cloudProfile%levNum-1
       idx = cloudProfile%layerIndex(j)
       ! set layer quantities
       cloudProfile%tempLayer(j) =  tavl(idx)

       ! set level quantities
       ratio = (cloudProfile%presLevel(j) -pref(idx)) / (pref(idx+1)- pref(idx))
       cloudProfile%tempLevel(j) =  ratio*(xG(idx+1)-xG(idx)) + xG(idx)
       cloudProfile%ratio(j) = 1.0 - ratio

       cloudProfile%molAbs(j) = (cloudProfile%presLevel(j+1) -cloudProfile%presLevel(j))/(pref(idx+1)  -pref(idx))
    enddo

    cloudProfile%tempLevel(cloudProfile%levNum) = xG(nSurf)

    ! identify the cloud slab location

    do iHydr=1,nHydLoc
       if (cloudProfile%cloudPresent(iHydr))then
          amtID = hydIndx(iHydr)+amtIndx(iHydr)
          topID = hydIndx(iHydr)+topIndx(iHydr)
          thkID = hydIndx(iHydr)+thkIndx(iHydr)

          cldTop = xG(topID)
          cldThk = xG(thkID)
          cldBot = cldThk + cldTop
          cldThickness(iHydr) = cldThk

          do j1=1, cloudProfile%levNum
             if(levelTheSame(cldTop, cloudProfile%presLevel(j1))) exit
          enddo
          cloudProfile%cloudTop(iHydr) = j1

          do j2=j1, cloudProfile%levNum
             if(levelTheSame(cldBot, cloudProfile%presLevel(j2))) then
                exit
             end if
          enddo
          cloudProfile%cloudBot(iHydr) = j2

          !-----Calculate cloud amt in each layer
          do j=j1, j2-1
            cFrac(j, iHydr)=(cloudProfile%presLevel(j + 1)-cloudProfile%presLevel(j))/cldThk
            ltyp(j) = iHydr
          end do
          indexLayer(1, iHydr) = j1
          indexLayer(2, iHydr) = j2-1
       endif
    enddo
    !   - Given constant Deff, and T from profile, get opt. depth profiles
   do iHydr=1,nHydLoc
      if (cloudProfile%cloudPresent(iHydr)) then
         do iLyr=indexLayer(1, iHydr), indexLayer(2, iHydr)

            cldSet%value(iSizP) = xG(hydIndx(iHydr)+sizIndx(iHydr))
            cldSet%value(iSizS) = 0.0
            cldSet%value(iTemp) = cloudProfile%tempLayer(iLyr)

            call getCloudTab(hydName(iHydr),cldSet,cldOpt)

            tabs(iLyr,:)  = tabs(iLyr,:) + cldOpt(ikabs)%values * cFrac(iLyr,iHydr) * intAmt(iHydr)

            itscat = cldOpt(ikscat)%values * cFrac(iLyr,iHydr) * intAmt(iHydr)
            tscat(iLyr,:) = tscat(iLyr,:) + itscat(:)

            asym(iLyr,:)  = asym(iLyr,:)  + itscat(:) * cldOpt(igasym)%values

            if (JacobianCalculation) then
               ! we have 4 Jacobians with respect to amount, cloud top, thickness, size
               ! a_n
               tabsJac(iLyr,:,amtIndx(iHydr),iHydr) = &
                      cldOpt(ikabs)%values*cFrac(iLyr,iHydr)*sclPath(iHydr)
               ! w_n da_n/dp
               tabsJac(iLyr,:,sizIndx(iHydr),iHydr) = &
                      cldOpt(ikabs)%deriv(iSizP,:)*cFrac(iLyr,iHydr)*intAmt(iHydr)
               tabsJac(iLyr,:,4,iHydr)=&
                      cldOpt(ikabs)%deriv(iTemp,:)*cFrac(iLyr,iHydr)*intAmt(iHydr)

               ! b_n
               tscatJac(iLyr,:,amtIndx(iHydr),iHydr) = &
                      cldOpt(ikscat)%values*cFrac(iLyr,iHydr)*sclPath(iHydr)
               ! w_n db_n/dp
               tscatJac(iLyr,:,sizIndx(iHydr),iHydr) = &
                      cldOpt(ikscat)%deriv(iSizP,:)*cFrac(iLyr,iHydr)*intAmt(iHydr)
               tscatJac(iLyr,:,4             ,iHydr) = &
                      cldOpt(ikscat)%deriv(iTemp,:)*cFrac(iLyr,iHydr)*intAmt(iHydr)

               ! g_n
               asymJac(iLyr,:,amtIndx(iHydr),iHydr)= &
                      cldOpt(igasym)%values ! * cldOpt(ikscat)%values
               ! w_n*b_n*dg_n/dp
               asymJac(iLyr,:,sizIndx(iHydr),iHydr) = &
                      itscat(:) * cldOpt(igasym)%deriv(iSizP,:)
               asymJac(iLyr,:,4             ,iHydr) = &
                      itscat(:) * cldOpt(igasym)%deriv(iTemp,:)

               if (indexLayer(1, iHydr) < indexLayer(2, iHydr)) then
                  if (iLyr == indexLayer(2, iHydr)) then
                      cldTop = xG(hydIndx(iHydr)+topIndx(iHydr))
                      tabsJac(iLyr,:,thkIndx(iHydr),iHydr) =  &
                            cldOpt(ikabs)%values * (cloudProfile%presLevel(iLyr) - cldTop) &
                                                  * intAmt(iHydr)/cldThickness(iHydr)**2
                      tscatJac(iLyr,:,thkIndx(iHydr),iHydr) = &
                            cldOpt(ikscat)%values *(cloudProfile%presLevel(iLyr) - cldTop) &
                                                  * intAmt(iHydr)/cldThickness(iHydr)**2
                  else
                      tabsJac(iLyr,:,thkIndx(iHydr),iHydr) = - cldOpt(ikabs)%values  &
                            * cFrac(iLyr,iHydr) * intAmt(iHydr)/cldThickness(iHydr)
                      tscatJac(iLyr,:,thkIndx(iHydr),iHydr) = - cldOpt(ikscat)%values  &
                            * cFrac(iLyr,iHydr) * intAmt(iHydr)/cldThickness(iHydr)
                  end if
               end if
            end if
         end do
        if (JacobianCalculation) then
            if (indexLayer(1, iHydr) < indexLayer(2, iHydr)) then
               iLyr = indexLayer(1, iHydr)
               tabsJac(iLyr,:,topIndx(iHydr),iHydr)  = &
                        -cldOpt(ikabs)%values/cldThickness(iHydr) *intAmt(iHydr)
               tscatJac(iLyr,:,topIndx(iHydr),iHydr) = &
                        -cldOpt(ikscat)%values/cldThickness(iHydr)*intAmt(iHydr)

               iLyr = indexLayer(2, iHydr)
               tabsJac(iLyr,:,topIndx(iHydr),iHydr)  =  &
                         cldOpt(ikabs)%values/cldThickness(iHydr) *intAmt(iHydr)
               tscatJac(iLyr,:,topIndx(iHydr),iHydr) =  &
                         cldOpt(ikscat)%values/cldThickness(iHydr)*intAmt(iHydr)
            end if
        end if
      end if
   enddo

       !---The following loop could be replaced with a WHERE construct, although the
       !     current version of the SGI compiler does not have the proper libraries
       !     to support WHERE in DEBUG mode.
    do iLyr=1,cloudProfile%levNum-1
       do iSpc=1,nSpcPts
          if(tscat(iLyr,iSpc) > 0.0) then
             ! g = sum(s_n g_n)/ sum(s_n)
             asym(iLyr,iSpc) = asym(iLyr,iSpc) / tscat(iLyr,iSpc)
             if (JacobianCalculation) then
                ! dg/dw_n = (g_n - g)*b_n/ sum(s_n)
                ! dg/dp = (w_n (g_n - g) *db_n/dp + s_n dg_n/dp)/ sum(s_n)

                 do iHydr=1,nHydLoc
                    deltaG = asymJac(iLyr,iSpc,amtIndx(iHydr),iHydr) - asym(iLyr,iSpc)
                    asymJac(iLyr,iSpc,amtIndx(iHydr),iHydr)=  &
                               deltaG *tscatJac(iLyr,iSpc,amtIndx(iHydr),iHydr)

                    if (indexLayer(1, iHydr) .ne. indexLayer(2, iHydr)) then
                        asymJac(iLyr,iSpc,topIndx(iHydr),iHydr)=  &
                              deltaG *tscatJac(iLyr,iSpc,topIndx(iHydr),iHydr)

                        asymJac(iLyr,iSpc,thkIndx(iHydr),iHydr)=  &
                              deltaG *tscatJac(iLyr,iSpc,thkIndx(iHydr),iHydr)
                    end if

                    asymJac(iLyr,iSpc,4,iHydr) = &
                         tscatJac(iLyr,iSpc,4,iHydr)*deltaG + asymJac(iLyr,iSpc,4,iHydr)

                    asymJac(iLyr,iSpc,sizIndx(iHydr),iHydr) = &
                          tscatJac(iLyr,iSpc,sizIndx(iHydr),iHydr)*deltaG  &
                          + asymJac(iLyr,iSpc,sizIndx(iHydr),iHydr)
                 end do
                 asymJac(iLyr,iSpc,:,:)=asymJac(iLyr,iSpc,:,:)/tscat(iLyr,iSpc)
             end if
          else
             asym(iLyr,iSpc) = 0.
             if (JacobianCalculation) asymJac(iLyr,iSpc,:,:)=  0.0
          endif
       enddo
    enddo
    ! ltyp may be passed to ossdrv and input to scattRT,
    !  to avoid repeatedly determining ltyp there.
    return
  END SUBROUTINE setCloudOptSlab
  ! ===================================

   subroutine insertLevel(pIn,cloudProfile)
!<f90Function>********************************************************
!
! NAME:
!
!   insertLevel
!
! PURPOSE:
!
!   insert a level to CloudProfileType
!
!*******************************************************</f90Function>
      real, intent(in)   :: pIn
      type(CloudProfileType), intent(inout)  :: cloudProfile

      integer j, k, idx

      if (pIn < cloudProfile%presLevel(1)) then
         print *,'pressure level to be inserted must be within pUser pressure grid range '
         print *,'pressure level to be inserted: ', pIn
         print *,'grid top pressure: ', cloudProfile%presLevel(1)
         call exit(1)
      end if
      if (pIn > cloudProfile%presLevel(cloudProfile%levNum)) then
         print *,'pressure level to be inserted must be within pUser pressure grid range '
         print *,'pressure level to be inserted: ', pIn
         print *,'grid bottom pressure: ', cloudProfile%presLevel(cloudProfile%levNum)
         call exit(1)
      end if

      do j=2, cloudProfile%levNum-1
         if(pIn < cloudProfile%presLevel(j)) exit
      enddo

      if (levelTheSame(pIn, cloudProfile%presLevel(j-1))) return

      cloudProfile%toUse = .true.
      do k=cloudProfile%levNum, j, -1
         idx = k+1
         cloudProfile%layerIndex(idx) = cloudProfile%layerIndex(k)
         cloudProfile%presLevel(idx)  = cloudProfile%presLevel(k)
      end do
      cloudProfile%presLevel(j)  = pIn
      cloudProfile%layerIndex(j) = cloudProfile%layerIndex(j-1)
      cloudProfile%levNum = cloudProfile%levNum + 1
      return
  end subroutine insertLevel

!----------------------------------------------------------------------------
  logical function levelTheSame(pIn, pTest)
!<f90Function>********************************************************
!
! NAME:
!
!   levelTheSame
!
! PURPOSE:
!
!   check if the pressure levels is the same with tolerance.
!
!*******************************************************</f90Function>
     real, intent(in)     :: pIn, pTest
     levelTheSame = abs(pTest - pIn) < 1d-5*pTest
  end function levelTheSame


!<f90Subroutine>********************************************************
!
! NAME:
!
!   destroyCloud
!
! PURPOSE:
!
!   Free-up memory.
!
! SYNTAX:
!
!   CALL destroyCloud()
!
! ARGUMENTS:
!
!   INPUTS:
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
  SUBROUTINE destroyCloud()

    integer :: iop

    if (allocated(kabs))   deallocate(kabs)
    if (allocated(kscat))  deallocate(kscat)
    if (allocated(asymi))  deallocate(asymi)
    if (allocated(itscat)) deallocate(itscat)
    if (allocated(indxArr))deallocate(indxArr)
    if (allocated(intAmt)) deallocate(intAmt)
    if (allocated(hydIndx))deallocate(hydIndx)
    if (allocated(cldOpt)) then
      do iop = 1, nOptical
         if (allocated(cldOpt(iop)%values)) deallocate(cldOpt(iop)%values)
         if (allocated(cldOpt(iop)%deriv)) deallocate(cldOpt(iop)%deriv)
      end do
      deallocate(cldOpt)
    end if
    if (allocated(cldSet%value)) deallocate(cldSet%value)
    if (allocated(hydName)) deallocate(hydName)
    call destroyCloudTab()
    return

  END SUBROUTINE destroyCloud

!=========================================================================
  function checkCloudModelType(testCloudModelType)
!<f90Subroutine>********************************************************
!
! NAME:
!
!   checkCloudModelType
!
! PURPOSE:
!
!   Check that the cloud model is the same as passed to initCloud
!
! SYNTAX:
!
!   checkCloudModelType(testCloudModelType)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   testCloudModelType    CHARACTER(len=*)  cloud model name to test
!
!   INPUTS/OUTPUTS:
!
!   .true. if modle supported
!   .false. if it is not saupported
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
    integer, intent(in)    :: testCloudModelType !
    logical                :: checkCloudModelType

    checkCloudModelType = testCloudModelType == moduleCloudModelType
  end function checkCloudModelType

!=========================================================================
  function getCloudModelType()
!<f90Subroutine>********************************************************
!
! NAME:
!
!   getCloudModelType
!
! PURPOSE:
!
!   return the cloud model ID that passed to initCloud
!
! SYNTAX:
!
!   getCloudModelType(testCloudModelType)
!
! ARGUMENTS:
!
!   INPUTS:
!
!
!   INPUTS/OUTPUTS:
!
!   cloud model ID
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
    integer           :: getCloudModelType
    getCloudModelType = moduleCloudModelType
  end function getCloudModelType


!=========================================================================
  subroutine printCloudProfile (cp)
!<f90Subroutine>********************************************************
!
! NAME:
!
!   printCloudProfile
!
! PURPOSE:
!
!   print function for CloudProfileType
!
! SYNTAX:
!
!   printCloudProfile (cp)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   cp    type(CloudProfileType)  CloudProfileType
!
!   INPUTS/OUTPUTS:
!
!   None
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>
     type(CloudProfileType), intent(in) :: cp

     integer jj

     print *,'toUse  :', cp%toUse
     print *,'levNum :', cp%levNum

     print '("present   :", 2L6)', cp%cloudPresent(1:2)
     print '("cloudTop  :", 2I6)', cp%cloudTop(1:2)
     print '("cloudBot  :", 2I6)', cp%cloudBot(1:2)

     do jj=1, cp%levNum-1
        if (abs(cp%molAbs(jj)-1.) < epsilon(1.0)) cycle
        print '("lev: ", I3, "   ", 2F8.2, 1p100E12.4)',  jj, cp%presLevel(jj), &
                              cp%tempLevel(jj), cp%ratio(jj)
        print '("lay: ", "   ", I3, 2F8.2, 1p100E12.4 )', cp%layerIndex(jj), &
                              cp%presLevel(jj+1) - cp%presLevel(jj), &
                              cp%tempLayer(jj),   cp%molAbs(jj), cp%tau(jj)
     end do
     jj=cp%levNum
     print '("lev: ", I3, "   ", 2F8.2, 1p100E12.4)',  jj, cp%presLevel(jj), &
                              cp%tempLevel(jj), cp%ratio(jj)
  end subroutine printCloudProfile
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! The following interpolation and layer averaging routines were copied here from
  !   oss_mw as a short-term fix.  They should either become public subroutines of
  !   oss_mw, or be placed in a more general module.

  !----------------------------------------------------------------------------
  ! PURPOSE: Computes average layer quantities (or integrated amount)
  !          using a linear dependence on P.
  !----------------------------------------------------------------------------

  !YHE: This is a general process that should be part of the
  !mathmatics module with "public" attribute to be available for all
  !algorithms that need it.
  SUBROUTINE lpsum(pu,pl,xu,xl,scal,xint,dxu,dxl)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   lpsum
!
! PURPOSE:
!
!   Computes average layer quantities (or integrated amount) using a linear
!   dependence on P.
!
! SYNTAX:
!
!   CALL lpsum(pu, pl, xu, xl, scal, xint, dxu, dxl)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   pu    REAL  Upper level pressure
!   pl    REAL  Lower level pressure
!   xu    REAL  Quantity to be integrated on upper level
!   xl    REAL  Quantity to be integrated on lower level
!   scal  REAL  Scale factor
!
!   INPUTS/OUTPUTS:
!
!   xint  REAL  Integrated quantity
!   dxu   REAL  Layer-to-upper-level derivative
!   dxl   REAL  Layer-to-lower-level derivative
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !---Input variables
    REAL, INTENT(IN)    :: pu,pl
    REAL, INTENT(IN)    :: xu,xl,scal

    !---Output variables
    REAL, INTENT(OUT)   :: xint,dxu,dxl

    xint = 0.5*(pl-pu)*scal
    dxu      = xint
    dxl      = xint
    xint     = xint*(xu+xl)
    return
  END SUBROUTINE lpsum

  !----------------------------------------------------------------------------
  ! PURPOSE: Computes average layer quantities (or integrated amount)
  !          using a log-x dependence on log-p.
  !----------------------------------------------------------------------------

  !YHE: same argument as in lpsum, i.e., should be in a more general module.

  SUBROUTINE lpsum_log(pu,pl,xu,xl,scal,xint,dxu,dxl)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   lpsum_log
!
! PURPOSE:
!
!   Computes average layer quantities (or integrated amount) using a log-x
!   dependence on log-p.
!
! SYNTAX:
!
!   CALL lpsum_log(pu, pl, xu, xl, scal, xint, dxu, dxl)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   pu    REAL  Upper level pressure
!   pl    REAL  Lower level pressure
!   xu    REAL  Quantity to be integrated on upper level
!   xl    REAL  Quantity to be integrated on lower level
!   scal  REAL  Scale factor
!
!   INPUTS/OUTPUTS:
!
!   xint  REAL  Integrated quantity
!   dxu   REAL  Layer-to-upper-level derivative
!   dxl   REAL  Layer-to-lower-level derivative
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !---Input variables
    REAL, INTENT(IN)    :: pu,pl,xu,xl,scal
    !---Output variables
    REAL, INTENT(INOUT) :: xint,dxu,dxl
    !---Local variables
    REAL, PARAMETER   :: epsiln=1.e-05
    REAL              :: hp,x0,zeta,alza,alpha
    hp       = alog(pl/pu)
    x0       = pl*xl*hp
    zeta     = pu*xu/(pl*xl)
    IF(ABS(zeta-1.).GT.epsiln)THEN
       alza  = alog(zeta)
       xint  = x0*(zeta-1.)/alza
       alpha = zeta/(zeta-1.)-1./alza
    ELSE
       xint  = x0*2./(3.-zeta)
       alpha = zeta/(3.-zeta)
    END IF
    xint     = xint*scal
    dxu      = xint*alpha/xu
    dxl      = xint*(1.-alpha)/xl
    return
  END SUBROUTINE lpsum_log

  !YHE: same argument as in lpsum.

  SUBROUTINE lint(xInp,pGrid,N0,p0,x0)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   lint
!
! PURPOSE:
!
!   Interpolation in linear pressure
!
! SYNTAX:
!
!   CALL lint(xInp, pGrid, N0, p0, x0)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   xInp   REAL     Input data
!   pGrid  REAL     Pressure grid
!   N0     INTEGER  Level index pointing to interpolation range
!   p0     REAL     Pressure at which interpolation to be done
!
!   INPUTS/OUTPUTS:
!
!   x0     REAL     interpolated value
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    INTEGER,            INTENT(IN)    :: n0
    REAL,               INTENT(IN)    :: p0
    REAL, DIMENSION(:), INTENT(IN)    :: xInp
    REAL, DIMENSION(:), INTENT(IN)    :: pGrid
    !---Output variable
    REAL,               INTENT(INOUT) :: x0
    !---Local variable
    REAL                           :: xx
    xx = (xInp(N0)-xInp(N0-1))/(pGrid(N0)-pGrid(N0-1))
    x0 = xInp(N0-1)+(p0-pGrid(N0-1))*xx
    return
  END SUBROUTINE lint

  !YHE: same argument as for lpsum

  SUBROUTINE lint_log(xInp,pGrid,N0,p0,x0)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   lint_log
!
! PURPOSE:
!
!   Interpolation in log pressure
!
! SYNTAX:
!
!   CALL lint_log(xInp, pGrid, N0, p0, x0)
!
! ARGUMENTS:
!
!   INPUTS:
!
!   xInp   REAL     Input data
!   pGrid  REAL     Pressure grid
!   N0     INTEGER  Level index pointing to interpolation range
!   p0     REAL     Pressure at which interpolation to be done
!
!   INPUTS/OUTPUTS:
!
!   x0     REAL     interpolated value
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    INTEGER,            INTENT(IN)    :: n0
    REAL,               INTENT(IN)    :: p0
    REAL, DIMENSION(:), INTENT(IN)    :: xInp,pGrid
    !---Output variable
    REAL,               INTENT(INOUT) :: x0
    !---Local variable
    REAL                           :: xx
    xx = log(xInp(N0)/xInp(N0-1))/log(pGrid(N0)/pGrid(N0-1))
    x0 = xInp(N0-1)*(p0/pGrid(N0-1))**xx
    return
  END SUBROUTINE lint_log
END MODULE CloudModule
