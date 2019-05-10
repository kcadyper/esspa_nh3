!--------------------------------------------------------------
!
!  MODULE ClassSfcMod: contains procedures for surface classification
!>>This crimss version is temporary. It allows radflt in scncov
!  to build. Later, SurfaceClassModule should be revised to handle
!  surface classification, and that should be done in a way that
!  supports radflt, similarly to the way util/retr/AtmosClassModule.f90
!  does atmosphere classification.
!
!--------------------------------------------------------------
MODULE ClassSfcMod
  IMPLICIT NONE
  PRIVATE

  !---------------------------------------------------------------
  !  List of Public subroutines (accessible from outside module)
  !---------------------------------------------------------------
  PUBLIC :: initHighEmisTest,highEmisTest,sfcIceTest

  !------------------------------------------------------------------------
  ! Overloaded subroutines
  !------------------------------------------------------------------------
  ! Integer kchan is used by scncov, which is built with mwsensors 
  ! and irsensors
  INTERFACE highEmisTest
     MODULE PROCEDURE highEmisTest_integer
     MODULE PROCEDURE highEmisTest_logical
  END INTERFACE
  INTERFACE sfcIceTest
     MODULE PROCEDURE sfcIceTest_integer
     MODULE PROCEDURE sfcIceTest_logical
  END INTERFACE

  !---------------------------------------------------------------
  !     Declarations of private data (private to the module)
  !---------------------------------------------------------------

CONTAINS

  subroutine initHighEmisTest(frq,pol)
    !------------------------------------------------------------
    !   Initialize highEmisTest
    !   Either called explicitly or invoked automatically on first call
    !   to highEmisTest.
    !------------------------------------------------------------

    real,              dimension(:), optional, intent(in)    :: frq
    integer,           dimension(:), optional, intent(in)    :: pol

    return
  end subroutine initHighEmisTest

!----------------------------------------------------------------------------

  subroutine highEmisTest_logical(Y,kchan,iHiEmFlag)

    real,    dimension(*), intent(in)    :: Y
    logical, dimension(*), intent(in)    :: kchan
    integer,               intent(inout) :: iHiEmFlag

    iHiEmFlag=0

    return
  end subroutine highEmisTest_logical

!----------------------------------------------------------------------------

!----------------------------------------------------------------------------

  subroutine highEmisTest_integer(Y,kchan,iHiEmFlag)

    real,    dimension(:), intent(in)    :: Y
    integer, dimension(:), intent(in)    :: kchan
    integer,               intent(inout) :: iHiEmFlag

    logical, dimension(size(kchan)) :: kchanLogical

    kchanLogical(:)=kchan(:)/=0
    call highEmisTest_logical(Y,kchanLogical,iHiEmFlag)

    return
  end subroutine highEmisTest_integer

!----------------------------------------------------------------------------

  subroutine sfcIceTest_logical(Y, iSfcIceFlag,kchan,frq,pol)
!************************************************************
!* Function Name: sfcIceTest
!* Purpose: Test brightness temperatures for signature of
!*          sea ice or fresh water ice
!* Usage: call sfcIceTest(Y, kchan, iSfcIceFlag)
!* Description: Test brightness temperatures for signature of
!*              sea ice or fresh water ice
!*
!* Inputs:
!*  Var_Name     Type       Description
!*  --------     ----       -----------
!*  Y            real       radiance array.
!*  kchan        logical    channel on/off flag array.
!*
!* Outputs:
!*  Var_Name     Type       Description
!*  --------     ----       -----------
!*  iSfcIceFlag  integer    flag is 1 if ice surface
!*                          is detected and is 0 otherwise
!* Includes:
!* File_Name       Description
!* ---------       -----------
!* crims.incl      defines various dimensions and parameters for the
!*                 retrieval
!* ossdrv_mw.com   maximum #'s of OSS points for MW, common blocks for
!*                 gas optical depths and reference profiles
!*
!* Copyright: AER, Inc., 2000
!* Developed by Atmospheric and Environmental Research, Inc.
!************************************************************


    real, dimension(*):: Y
    integer   :: iSfcIceFlag
    logical   :: kchan(*)
    real      :: frq(*)
    integer   :: pol(*)

    iSfcIceFlag=0

    return
  end subroutine sfcIceTest_logical

!----------------------------------------------------------------------------

  subroutine sfcIceTest_integer(Y, iSfcIceFlag,kchan,frq,pol)

    real, dimension(*):: Y
    integer   :: iSfcIceFlag
    integer   :: kchan(:)
    real      :: frq(:)
    integer   :: pol(:)

    logical, dimension(size(kchan)) :: kchanLogical

    kchanLogical(:)=kchan(:)/=0
    call sfcIceTest_logical(Y, iSfcIceFlag,kchanLogical,frq,pol)

    return
  end subroutine sfcIceTest_integer

END MODULE ClassSfcMod
