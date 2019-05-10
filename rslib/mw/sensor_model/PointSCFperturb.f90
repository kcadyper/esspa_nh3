MODULE PointSCFperturb
  
  implicit none
  private
  public :: getPointSCF, perturbPointing
  
  integer,parameter :: lID=12
  integer,parameter :: cmtLine=200
  character(len=200),parameter :: BLANKF=''
  
!-- Global private data  
  
  INTEGER :: iu
  CHARACTER(LEN=50) :: fname
  
  !-- Pointing SCF:
  integer  :: ifmtver
  character(cmtLine) :: cmt1,cmt2,cmt3
  integer            :: nbeam
  real               :: s2sc_roll,s2sc_pitch,s2sc_yaw
  real               :: a2s_roll,a2s_pitch,a2s_yaw
  real               :: sc_yaw
  character(lID),dimension(:),pointer :: beamID
  real,          dimension(:),pointer :: patt_Nd
  real,          dimension(:),pointer :: patt_Az
  
  !-- Pointing Errors SCF:
  real     :: boreAlongErr,boreCrossErr
  real,dimension(:),pointer :: beamRegAlongErr,beamRegCrossErr
  
  !-- Scan Orbit SCF:
  real     :: azRef,scanRate,EarthRad,scAlt,orbitIncl
  integer  :: rotDir,nbeamSOrb
  real,dimension(:),pointer :: sTimeEv,sTimeCs,sTimeW,tInteg,sampInt,sampFrq
  integer,dimension(:),pointer :: nsamp,ncssamp,nwlsamp
  character(lID),dimension(:),pointer :: beamIDsOrb
  
CONTAINS
  
  SUBROUTINE getPointSCF(U_Point,U_PointErr,U_ScanOrbit, &
       F_Point,F_PointErr,F_ScanOrbit)
    
    USE SCFread
    
    integer, intent(in)  :: U_Point,U_PointErr,U_ScanOrbit
    character(len=*),intent(in) :: F_Point,F_PointErr,F_ScanOrbit
    
    CALL getSCFpointing(iu=U_Point,fname=TRIM(F_Point), &
         s2sc_roll=s2sc_roll,s2sc_pitch=s2sc_pitch,s2sc_yaw=s2sc_yaw, &
         a2s_roll=a2s_roll,a2s_pitch=a2s_pitch,a2s_yaw=a2s_yaw, &
         sc_yaw=sc_yaw,nbeam=nbeam,beamID=beamID, &
         patt_Nd=patt_Nd,patt_Az=patt_Az)
    
    CALL getSCFpointingErr(iu=U_PointErr,fname=TRIM(F_PointErr), &
         boreAlongErr=boreAlongErr,boreCrossErr=boreCrossErr, &
         beamRegAlongErr=beamRegAlongErr,beamRegCrossErr=beamRegCrossErr)
    
    CALL getSCFscanOrbit(iu=U_ScanOrbit,fname=TRIM(F_ScanOrbit), &
         rotDir=rotDir,scanRate=scanRate,azRef=azRef, &
         nbeam=nbeamSOrb,EarthRad=EarthRad,scAlt=scAlt, &
         orbitIncl=orbitIncl,beamID=beamIDsOrb,nsamp=nsamp, &
         sTimeEv=sTimeEv,ncssamp=ncssamp,sTimeCs=sTimeCs, &
         nwlsamp=nwlsamp,sTimeW=sTimeW,tInteg=tInteg, &
         sampInt=sampInt,sampFrq=sampFrq)    
    
  END SUBROUTINE getPointSCF
  
  SUBROUTINE perturbPointing(dt_rpy,RN_rpyS2SC,RN_rpyA2S,RN_scYaw, &
       RN_boreAlong,RN_boreCross,RN_sTimeEv,RN_azRef, &
       U_rpyS2SC,U_rpyA2S,U_scYaw, &
       U_boreAlong,U_boreCross,U_sTimeEv,U_azRef,F_rpyS2SC, &
       F_rpyA2S,F_SCyaw, &
       F_boreAlong,F_boreCross,F_sTimeEv,F_azRef, &
       U_point,F_point,U_ScanOrbit,F_ScanOrbit,U_outPoint,F_outPoint, &
       U_outScOrb,F_outScOrb)
    
    real,        intent(in),optional :: dt_rpy
    logical,     intent(in),optional :: RN_rpyS2SC
    logical,     intent(in),optional :: RN_rpyA2S
    logical,     intent(in),optional :: RN_scYaw
    logical,     intent(in),optional :: RN_boreAlong
    logical,     intent(in),optional :: RN_boreCross
    logical,     intent(in),optional :: RN_sTimeEv
    logical,     intent(in),optional :: RN_azRef
    integer,     intent(in),optional :: U_rpyS2SC
    integer,     intent(in),optional :: U_rpyA2S
    integer,     intent(in),optional :: U_scYaw
    integer,     intent(in),optional :: U_boreAlong
    integer,     intent(in),optional :: U_boreCross
    integer,     intent(in),optional :: U_sTimeEv
    integer,     intent(in),optional :: U_azRef
    integer,     intent(in),optional :: U_point
    integer,     intent(in),optional :: U_ScanOrbit
    integer,     intent(in),optional :: U_outPoint
    integer,     intent(in),optional :: U_outScOrb
    character(*),intent(in),optional :: F_rpyS2SC
    character(*),intent(in),optional :: F_rpyA2S
    character(*),intent(in),optional :: F_SCyaw
    character(*),intent(in),optional :: F_boreAlong
    character(*),intent(in),optional :: F_boreCross
    character(*),intent(in),optional :: F_sTimeEv
    character(*),intent(in),optional :: F_azRef
    character(*),intent(in),optional :: F_point
    character(*),intent(in),optional :: F_ScanOrbit
    character(*),intent(in),optional :: F_outPoint
    character(*),intent(in),optional :: F_outScOrb
    
!-- Local variables
    real      :: gauss
    real      :: s2sc_rollF,s2sc_pitchF,s2sc_yawF
    real      :: a2s_rollF,a2s_pitchF,a2s_yawF
    real      :: sc_yawF,sTimeEvF,azRefF
    real      :: beamRegAlongF,beamRegCrossF
    real,dimension(:),allocatable :: sTimeEv0,patt_Nd0,patt_Az0
    real      :: s2sc_roll0,s2sc_pitch0,s2sc_yaw0
    real      :: a2s_roll0,a2s_pitch0,a2s_yaw0,sc_yaw0
    real      :: azRef0
    integer,parameter :: ix=1133   ! seed for gauss function
    integer   :: i
    
    s2sc_roll0  = s2sc_roll
    s2sc_pitch0 = s2sc_pitch
    s2sc_yaw0   = s2sc_yaw
    a2s_roll0  = a2s_roll
    a2s_pitch0 = a2s_pitch
    a2s_yaw0   = a2s_yaw
    sc_yaw0 = sc_yaw    
    allocate(patt_Nd0(nbeam),patt_Az0(nbeam))
    patt_Nd0 = patt_Nd
    patt_Az0 = patt_Az
    allocate(sTimeEv0(nbeam))
    sTimeEv0 = sTimeEv
    azRef0   = azRef
    
    !  Apply errors to [sensor to spacecraft] roll, pitch, yaw:
    if (present(RN_rpyS2SC)) then
       if (RN_rpyS2SC) then
          s2sc_roll   = s2sc_roll+dt_rpy*gauss(ix,1.,0.)
          s2sc_pitch  = s2sc_pitch+dt_rpy*gauss(ix,1.,0.)
          s2sc_yaw    = s2sc_yaw+dt_rpy*gauss(ix,1.,0.)
       endif
    endif
    
    !  Apply errors to [antenna to sensor] roll, pitch, yaw:
    if (present(RN_rpyA2S)) then
       if (RN_rpyA2S) then
          a2s_roll   = a2s_roll+dt_rpy*gauss(ix,1.,0.)
          a2s_pitch  = a2s_pitch+dt_rpy*gauss(ix,1.,0.)
          a2s_yaw    = a2s_yaw+dt_rpy*gauss(ix,1.,0.)
       endif
    endif
    
    ! Apply error to nominal spacecraft yaw:
    if (present(RN_scYaw)) then
       if (RN_scYaw) then
          sc_yaw  =  sc_yaw+dt_rpy*gauss(ix,1.,0.)
       endif
    endif
    
    !----------------------------------------------------------
    ! Apply errors using fixed scaling factor (i.e. bias), if the
    !  option is selected:    
    
    if (present(F_rpyS2SC)) then
       if (.not. present(U_rpyS2SC)) then
          print *,'err[PointSFCPerturb::perturbPointing]: ',&
               ' U_rpyS2SC must come with F_rpyS2SC'
          call errorHalt(1)
       endif
       
       if (F_rpyS2SC /= BLANKF) then
          if (present(RN_rpyS2SC)) then
             if (RN_rpyS2SC) then
                print*,'message from perturbPointing: '
                print*,'overriding random SCF perturbation'
             endif
          endif
          
          open(U_rpyS2SC,file=F_rpyS2SC,status='old',action='read')
          read(U_rpyS2SC,'(3f9.6)') s2sc_rollF,s2sc_pitchF,s2sc_yawF
          
          s2sc_roll  = s2sc_roll0+s2sc_rollF
          s2sc_pitch = s2sc_pitch0+s2sc_pitchF
          s2sc_yaw   = s2sc_yaw0+s2sc_yawF
          close(U_rpyS2SC)
       endif
    endif

    if (present(F_rpyA2S)) then
       if (.not. present(U_rpyA2S)) then
          print *,'err[PointSFCPerturb::perturbPointing]: ',&
               ' U_rpyA2S must come with F_rpyA2S'
          call errorHalt(1)
       endif
       
       if (F_rpyA2S /= BLANKF) then
          open(U_rpyA2S,file=F_rpyA2S,status='old',action='read')
          read(U_rpyA2S,'(3f9.6)') a2s_rollF,a2s_pitchF,a2s_yawF
          
          a2s_roll  = a2s_roll0+a2s_rollF
          a2s_pitch = a2s_pitch0+a2s_pitchF
          a2s_yaw   = a2s_yaw0+a2s_yawF
          close(U_rpyA2S)
       endif
    endif    
    
    if (present(F_SCyaw)) then
       if (.not. present(U_scYaw)) then
          print *,'err[PointSFCPerturb::perturbPointing]: ',&
               ' U_scYaw must come with F_scYaw'
          call errorHalt(1)
       endif

       if (F_SCyaw /= BLANKF) then
          open(U_scYaw,file=F_SCyaw,status='old',action='read')
          read(U_scYaw,'(f9.6)') sc_yawF
          
          sc_yaw  = sc_yaw0+sc_yawF
       endif
    endif

    ! Apply errors to boresight along and cross-scan pointing:    
    if (present(RN_boreAlong) .OR. present(RN_boreCross)) then
       if (RN_boreAlong .OR. RN_boreCross) then
          do i=1,nbeam
             if (RN_boreCross) &
                  patt_Nd(i) = patt_Nd(i) + &
                  boreCrossErr*gauss(ix,1.,0.) + &
                  beamRegCrossErr(i)*gauss(ix,1.,0.)
             
             if (RN_boreAlong) &
                  patt_Az(i) = patt_Az(i) + beamRegAlongErr(i)*gauss(ix,1.,0.)
          enddo
       endif
    endif
    
    if (present(F_boreAlong)) then
       if (.not. present(U_boreAlong)) then
          print *,'err[PointSFCPerturb::perturbPointing]: ',&
               ' U_scYaw must come with F_scYaw'
          call errorHalt(1)
       endif
       if (F_boreAlong /= BLANKF) then
          open(U_boreAlong,file=F_boreAlong,status='old',action='read')
          read(U_boreAlong,'(f7.4)') beamRegAlongF
          do i=1,nbeam
             patt_Az(i) = patt_Az0(i) + beamRegAlongErr(i)*beamRegAlongF
          enddo
          close (U_boreAlong)
       endif
    endif
    
    if (present(F_boreCross)) then
       if (.not. present(U_boreCross)) then
          print *,'err[PointSFCPerturb::perturbPointing]: ',&
               ' U_scYaw must come with F_scYaw'
          call errorHalt(1)
       endif
       if (F_boreCross /= BLANKF) then
          open(U_boreCross,file=F_boreCross,status='old',action='read')
          read(U_boreCross,'(f7.4)') beamRegCrossF
          do i=1,nbeam
             patt_Nd(i) = patt_Nd0(i) + (boreCrossErr*beamRegCrossF) + &
                  (beamRegCrossErr(i)*beamRegCrossF)
          enddo
          close (U_boreCross)
       endif
    endif    
    deallocate(patt_Nd0,patt_Az0)
    
    !-------------
    !  Scan Orbit SCF:
    if (present(RN_sTimeEv)) then
       if (RN_sTimeEv) then
          do i=1,nbeam
             sTimeEv(i) = sTimeEv(i) + dt_rpy*gauss(ix,1.,0.)
          enddo
       endif
    endif
    if (present(RN_azRef)) then 
       if (RN_azRef) azRef = azRef + dt_rpy*gauss(ix,1.,0.)
    endif
    
    if (present(F_sTimeEv)) then
       if (.not. present(U_sTimeEv)) then
          print *,'err[PointSFCPerturb::perturbPointing]: ',&
               ' U_sTimeEv must come with F_sTimeEv'
          call errorHalt(1)
       endif

       if (F_sTimeEv /= BLANKF) then
          if (present(RN_sTimeEv)) then
             if (RN_sTimeEv) then
                print*,'message from perturbPointing: '
                print*,'overriding random SCF perturbation of sTimeEv'
             endif
          endif
          
          open(U_sTimeEv,file=F_sTimeEv,status='old',action='read')
          read(U_sTimeEv,'(f9.6)') sTimeEvF
          do i=1,nbeam
             sTimeEv(i) = sTimeEv0(i) + sTimeEvF
          enddo
          close (U_sTimeEv)
       endif
       
    endif
    deallocate(sTimeEv0)
    
    if (present(F_azRef)) then
       if (.not. present(U_azRef)) then
          print *,'err[PointSFCPerturb::perturbPointing]: ',&
               ' U_azRef must come with F_azRef'
          call errorHalt(1)
       endif
       
       if (F_azRef /= BLANKF) then
          if (present(RN_azRef)) then
             if (RN_azRef) then
                print*,'message from perturbPointing: '
                print*,'overriding random SCF perturbation of azRef'
             endif
          endif
          
          open(U_azRef,file=F_azRef,status='old',action='read')
          read(U_azRef,'(f7.2)') azRefF
          azRef = azRef0 + azRefF
          close (U_azRef)
       endif
    endif
    
    ! ------------
    !  Write perturbed Pointing SCF:

    open(U_point,file=F_point,status='old',action='read')
    read(U_point,'(i4)') ifmtver
    read(U_point,*) cmt1
    read(U_point,*) cmt2
    read(U_point,*) cmt3
    close(U_point)
    
    open(U_outPoint,file=F_outPoint,status='replace',action='write')
    write(U_outPoint,'(i4)') ifmtver
    write(U_outPoint,*) TRIM(cmt1)
    write(U_outPoint,*) TRIM(cmt2)
    write(U_outPoint,*) TRIM(cmt3)
    
    write(U_outPoint,'(3f8.4)') s2sc_roll,s2sc_pitch,s2sc_yaw
    write(U_outPoint,'(3f8.4)') a2s_roll,a2s_pitch,a2s_yaw
    write(U_outPoint,'(3f8.3)') sc_yaw
    write(U_outPoint,'(i4)') nbeam
    do i=1,nbeam
       write(U_outPoint,'(a12,2f8.3)') beamID(i),patt_Nd(i),patt_Az(i)
    enddo
    
    close(U_outPoint)  
    
    ! ------------
    !  Write perturbed Scan/Orbit SCF:
    
    open(U_ScanOrbit,file=F_ScanOrbit,status='old',action='read')
    read(U_ScanOrbit,'(i4)') ifmtver
    read(U_ScanOrbit,'(a150)') cmt1
    read(U_ScanOrbit,'(a150)') cmt2
    read(U_ScanOrbit,'(a150)') cmt3
    close(U_ScanOrbit)
    
    open(U_outScOrb,file=F_outScOrb,status='replace',action='write')
    write(U_outScOrb,'(i4)') ifmtver
    write(U_outScOrb,*) TRIM(cmt1)
    write(U_outScOrb,*) TRIM(cmt2)
    write(U_outScOrb,*) TRIM(cmt3)
    
    write(U_outScOrb,'(i2,f7.3,f7.2,i4,f9.3,f8.3,f7.3)') &
         rotDir,scanRate,azRef,nbeamSOrb,EarthRad,scAlt,orbitIncl
    
    do i=1,nbeamSOrb
       write(U_outScOrb,'(a12,1x,i4,f9.6,i3,f9.6,i3,3f9.6,f6.1)') &
            beamIDsOrb(i),nsamp(i),sTimeEv(i),ncssamp(i), &
            sTimeCs(i),nwlsamp(i),sTimeW(i),tInteg(i),sampInt(i), &
            sampFrq(i)
    enddo
    
    close(U_outScOrb)
    
  end SUBROUTINE perturbPointing
  
end module PointSCFperturb
