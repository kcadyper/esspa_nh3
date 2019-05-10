MODULE PolSim
  
  USE constants, only: deg2rad
  implicit none
  private
  public :: rotationAlign,rotationAlignUncert,attitudeUncert
  
CONTAINS
  
  SUBROUTINE rotationAlignUncert(rotalign_bias,rotalign_sdev,alignKnowl &
     ,    rot_uncert_out)
  
    REAL :: rot_uncert,rot_uncert_out,gauss, rotalign_bias
    REAL ::  rotalign_sdev,alignKnowl
    REAL,PARAMETER :: mean=0.0
    INTEGER,PARAMETER :: seed=110037428
    LOGICAL :: rand_rot
    
    IF (rotalign_bias.EQ.0) THEN
       ! If applying randomized rotational uncertainty, assign an 
       !   uncertainty determined using gaussian noise approach:
       rot_uncert_out = gauss(seed,rotalign_sdev,mean) * alignKnowl 
    ELSE
       ! If not using a random uncertainty, just set the output to 
       !  the value specified in the namelist:
       rot_uncert_out = rotalign_bias * alignKnowl  
    ENDIF
    
  END SUBROUTINE rotationAlignUncert
  
! ============================================================  
  
  SUBROUTINE attitudeUncert(att_bias,att_sdev,dt_rot,    &
       ant_rpy_err,sens_rpy_err,sc_rpy_err)
    
    REAL                 :: att_bias,att_sdev,gauss
    REAL,DIMENSION(3)    :: ant_rpy_err,sens_rpy_err
    REAL,DIMENSION(3)    :: sc_rpy_err
    REAL,PARAMETER :: mean=0.0
    INTEGER,PARAMETER :: seed=110037428
    REAL :: dt_rot
    INTEGER :: i
    
    IF (att_bias.EQ.0) THEN
       
       DO i=1,3
          ! Using a scale factor that is small compared to the
          !   typical rotations (typcially factor is ~0.01 deg)
          ant_rpy_err(i) = gauss(seed,att_sdev,mean) * dt_rot
          sens_rpy_err(i) = gauss(seed,att_sdev,mean) * dt_rot
          sc_rpy_err(i) = gauss(seed,att_sdev,mean) * dt_rot
       ENDDO
       
    ELSE
       ant_rpy_err = att_bias  
       sens_rpy_err = att_bias 
       sc_rpy_err = att_bias   
    ENDIF
    
  END SUBROUTINE attitudeUncert

! ============================================================  
  
  SUBROUTINE rotationAlign(tbstokes,pra,     TAP)
    
    REAL(4)       :: TBSTOKES(4),TA(6)
    REAL(4)       :: PRA,COS2X,SIN2X
    REAL(4)       :: TASTOKESP(4),TAP(9)
    REAL(4)       :: QUVEC(2),QUVECP(2),PMAT(2,2)
    
    !  INTEGER       :: i,j,k
    
    QUVEC = (/TBSTOKES(1)-TBSTOKES(2),TBSTOKES(3)/) 
    
    SIN2X = SIN(2.0*PRA*deg2rad)
    COS2X = COS(2.0*PRA*deg2rad)
    
    PMAT(1,1) = COS2X
    PMAT(1,2) =+SIN2X
    PMAT(2,1) =-SIN2X
    PMAT(2,2) = COS2X
    
    QUVECP(1:2) = MATMUL(PMAT(1:2,1:2),QUVEC(1:2)) 
    
    TASTOKESP(1) =  ((TBSTOKES(1)+TBSTOKES(2)) + QUVECP(1))/2.0
    TASTOKESP(2) =  ((TBSTOKES(1)+TBSTOKES(2)) - QUVECP(1))/2.0
    TASTOKESP(3) =  QUVECP(2)   ! rotated P-M
    TASTOKESP(4) =  TBSTOKES(4) ! L-R 
    
    TAP(1:2) = TASTOKESP(1:2)
    TAP(3)   = 0.5*(TASTOKESP(1)+TASTOKESP(2)+TASTOKESP(3))
    TAP(4)   = 0.5*(TASTOKESP(1)+TASTOKESP(2)-TASTOKESP(3))
    TAP(5)   = 0.5*(TASTOKESP(1)+TASTOKESP(2)+TASTOKESP(4))
    TAP(6)   = 0.5*(TASTOKESP(1)+TASTOKESP(2)-TASTOKESP(4))
    TAP(7)   = TASTOKESP(1) - TASTOKESP(2)  ! TQ
    TAP(8)   = TASTOKESP(3)                 ! TU
    TAP(9)   = TASTOKESP(4)                 ! T4
    
  END SUBROUTINE rotationAlign
  
END MODULE PolSim
