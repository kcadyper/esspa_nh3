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

!--------------------------------------------------------------
!
!  MODULE RegrPrepModule: contains procedures related to processing
!          retrieval regressionback data (background and regression coefs)
!
!--------------------------------------------------------------
MODULE RegrPrepModule
  USE StateIndexModule
  USE IRReadStdInputs, only: U_regr, F_regr
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------
  !  List of Public subroutines (accessible from outside module) 
  !---------------------------------------------------------------
  PUBLIC :: getRegrParam,loadRegr,setFOVRegr
  !---------------------------------------------------------------
  !  List of Public variables (accessible from outside module) 
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  !     Declarations of private data (private to the module)
  !---------------------------------------------------------------
  INTEGER, DIMENSION(:), allocatable      :: AtmTypes
  INTEGER                                 :: nAtmTypes
  INTEGER                                 :: nparGatm
  INTEGER                                 :: nchirregr,ntranregr
  INTEGER                                 :: neofregr
  INTEGER                                 :: mynLevel
  REAL,    DIMENSION(:), pointer          :: mypressure

  !---Regression Mean Matrix (atm and rad)
  REAL, DIMENSION(:,:), allocatable       :: back_regr_atm_all,back_regr_rad
  !---Regression Coef Matrix (atm)
  REAL, DIMENSION(:,:,:), allocatable     :: regr_coef_atm_all
  !---EOF Matrices (atm and rad)
  REAL, DIMENSION(:,:,:), allocatable     :: u_regr_atm_all,u_rad_all
  !---Radiance Weight Vector
  REAL, DIMENSION(:,:), allocatable       :: wght_regr_rad

  !---Indices
  TYPE(StateIndex_t)                      :: IG,NG

  ! loadDone implementation provides the option to call loadRegr() in advance
  logical                            :: loadDone=.false.

CONTAINS

  SUBROUTINE getRegrParam(nparGout,nparout,nchanout,neofout,&
       IGout,NGout,nlevel,pressure,molID)

    !--I/O variables
    integer,                  intent(inout) :: nparGout,nparout
    integer,                  intent(inout) :: nchanout,neofout
    TYPE(StateIndex_t),       intent(inout) :: IGout,NGout
    integer,                  intent(inout) :: nLevel
    real,    dimension(:),    pointer       :: pressure
    integer, dimension(:),    intent(in)    :: molid


    !------------------------------------------------------------------------
    ! Set global variables from input arguments
    !------------------------------------------------------------------------
    call loadRegr(F_regr,molID)

    !------------------------------------------------------------------------
    ! Copy to output arguments
    !------------------------------------------------------------------------
    allocate(pressure(mynlevel))
    IGout         = IG
    NGout         = NG
    nparGout      = nparGAtm
    nparout       = ntranregr
    nchanout      = nchirregr
    neofout       = neofregr
    nLevel        = mynLevel
    pressure      = mypressure

    RETURN
  END SUBROUTINE getRegrParam

  subroutine loadRegr(file_regr,molID)
    use Regr_io_module

    !--I/O variables
    character(len=*),              intent(in)  :: file_regr
    integer, dimension(:)                      :: molID
    !---Local variables
    integer                                    :: itype

    if(loadDone) return

    call queryRegr(file_regr,nLevel=mynLevel)
    allocate(mypressure(mynlevel))
    CALL openRegr(ncid=U_regr, file=file_regr,nbkg=nAtmTypes, molid=molID,&
         NParG=NParGatm,nchir=nchirregr,ntran=ntranregr,neof=neofregr, &
         pressure=mypressure,IG=IG,NG=NG)
    allocate(back_regr_atm_all(NParGAtm,nAtmTypes))
    allocate(u_regr_atm_all(NParGAtm,NParGAtm,nAtmTypes))
    allocate(u_rad_all(nchirregr,neofregr,nAtmTypes))
    allocate(back_regr_rad(nchirregr,nAtmTypes),wght_regr_rad(nchirregr,nAtmTypes))
    allocate(regr_coef_atm_all(neofregr,ntranregr,nAtmTypes))
    allocate(AtmTypes(nAtmTypes))
    DO itype = 1, nAtmTypes
       CALL getRegr(ncid=U_regr,irec=itype, nparg=npargatm, &
            ntran=ntranregr, nchir=nchirregr, neof=neofregr,&
            dmean=back_regr_atm_all(1:npargatm,itype),             &
            un_atm=u_regr_atm_all(1:npargatm,1:ntranregr,itype), &
            un_rad=u_rad_all(1:nchirregr,1:neofregr,itype), &
            drmean=back_regr_rad(1:nchirregr,itype), &
            wghtvec=wght_regr_rad(1:nchirregr,itype), &
            regrcoef=regr_coef_atm_all(1:neofregr,1:ntranregr,itype), &
            TypeFlag=AtmTypes(itype))
       IF (AtmTypes(itype) < 0) then
          print*,'err[RegrPrepModule::loadRegr]: AtmTypes undefined in Regr File'
          call errorHalt(1)
       endif
    ENDDO
    CALL closeRegr(ncid=U_regr)

    loadDone=.true.
    RETURN
  end subroutine loadRegr

  SUBROUTINE setFOVRegr(iclassatm,psfc,ifor,&
          xbackG,rback,wght,umtx,umtxrad,regrcoef)
    !---Input variables
    INTEGER,                 INTENT(IN)    :: iclassatm
    INTEGER,                 INTENT(IN)    :: ifor
    REAL,                    INTENT(IN)    :: psfc
    REAL,    DIMENSION(:),   INTENT(INOUT) :: xbackG,rback,wght
    REAL,    DIMENSION(:,:), INTENT(INOUT) :: umtx,umtxrad
    REAL,    DIMENSION(:,:), INTENT(INOUT) :: regrcoef

    wght(1:nchirregr)=wght_regr_rad(1:nchirregr,iclassatm)
    xbackG(1:npargatm)=back_regr_atm_all(1:npargatm,iclassatm)
    rback(1:nchirregr)=back_regr_rad(1:nchirregr,iclassatm)
    umtx(1:npargatm,1:ntranregr)=&
         u_regr_atm_all(1:npargatm,1:ntranregr,iclassatm)
    umtxrad(1:nchirregr,1:neofregr)=&
         u_rad_all(1:nchirregr,1:neofregr,iclassatm)
    regrcoef(1:neofregr,1:ntranregr)=&
         regr_coef_atm_all(1:neofregr,1:ntranregr,iclassatm)

    RETURN
  END SUBROUTINE setFOVRegr
  
END MODULE RegrPrepModule
