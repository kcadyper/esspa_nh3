MODULE LimitsModule
  !
  ! Module for imposing limits over physical properties
  !
  ! Subroutines:
  !
  !     limitsInit
  !     imposeLimits (call chkges)
  !     limitsDestroy
  !
  ! Derived data types:
  !
  !     None
  !
  ! USE:
  !
  !     ControlStructure
  ! 
  ! yhe@aer.com, 03/16/2016
  !
  USE ToolboxModule, Only: &
       getUnit

  USE ControlStructure, Only: &
       GenControl_t, &
       StateSpaceDefin_t
  
  USE ChkValid, Only: &
       loadGuessLimits

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: &
       limitsInit, &
       imposeLimits, &
       limitsDestroy

  integer, parameter :: mxLength = 256
  real, dimension(:,:), allocatable :: clw_cov
  real, dimension(:,:), allocatable :: ciw_cov
  logical :: dbg = .false.

CONTAINS

  subroutine limitsInit( &
       genControl &
!       genControl, &
!       stateSpaceDefin &
       )

    TYPE(GenControl_t), intent(in) :: genControl
!    TYPE(StateSpaceDefin_t), intent(in) :: stateSpaceDefin
    !Local
    character (len=mxLength) :: procName
    integer :: nCloudLiq
    integer :: nCloudIce
    integer :: U_limCfg

    procName = "[LimitsModule::limitsInit]:"
    if (dbg) print *, trim(procName)//"starting ..."
! Regarding the lines commented out, see irsensors issue #2 in GitLab
! By commenting the lines out, the code calling limitsInit does not
! need to provide a stateSpaceDefin while this code is not modular enough
! to make appropriate use of a stateSpaceDefin.
!    nCloudLiq = stateSpaceDefin%NG%cldLiq
!    nCloudIce = stateSpaceDefin%NG%cldIce
    U_limCfg = getUnit() 
    call loadGuessLimits(U_limCfg,genControl%limitsConfig)
!    if (.not. allocated(clw_cov)) allocate(clw_cov(nCloudLiq, nCloudLiq))
!    if (.not. allocated(ciw_cov)) allocate(ciw_cov(nCloudIce, nCloudIce))

    if (dbg) print *, trim(procName)//"ending ..."

  end subroutine limitsInit

  subroutine imposeLimits()
    !Local
    character (len=mxLength) :: procName

    procName = "[LimitsModule::imposeLimits]:"
    if (dbg) print *, trim(procName)//"starting ..."
    if (dbg) print *, trim(procName)//"ending ..."

  end subroutine ImposeLimits

  subroutine limitsDestroy()
    !Local
    character (len=mxLength) :: procName

    procName = "[LimitsModule::limitsDestroy]:"
    if (dbg) print *, trim(procName)//"starting ..."
    if (allocated(clw_cov)) deallocate(clw_cov)
    if (allocated(ciw_cov)) deallocate(ciw_cov)
    if (dbg) print *, trim(procName)//"ending ..."

  end subroutine limitsDestroy
       
END MODULE LimitsModule
