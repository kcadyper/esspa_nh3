!--------------------------------------------------------------
!
!  MODULE ClassAtmMod: contains procedures for atmosphere classification
!                      based on a neural network
!                      This is for compatibility with mwsensors
!
!--------------------------------------------------------------
MODULE ClassAtmMod

  USE AtmosClassModule, ONLY: &
     atmosClassInit, &
     putAtmosClassData, &
     classAtmNN, &
     rdwts, &
     AtmosClassData_t, &
     MWmode

  USE ControlStructure, Only: &
       GenControl_t

  USE MWobsStructure, Only: &
       MWdefin_t

  USE IRobsStructure, Only: &
       IRdefin_t

  IMPLICIT NONE
  PRIVATE

  !---------------------------------------------------------------
  !  List of Public subroutines (accessible from outside module)
  !---------------------------------------------------------------
  PUBLIC :: classatm

  !---------------------------------------------------------------
  !     Declarations of private data (private to the module)
  !---------------------------------------------------------------
  INTEGER, DIMENSION(:), ALLOCATABLE :: QC
  LOGICAL, SAVE :: first=.true.

CONTAINS

!----------------------------------------------------------------------------

  SUBROUTINE classatm(Ym,kchan,xlat,xlon,eia,time,igeo,iclassatm)
    !------------------------------------------------------------
    !   Subroutine to find the atmospheric class using a probablistic
    !   neural network.  An array
    !   (profile) of brightness temperatures is used as the input set,
    !   and the network classifies this array into a specific class.
    !   Upon the first entry to this routine,
    !   the results of that network training are read from a file.
    !   If any channels required for classification are missing, the
    !   subroutine returns iclassatm=natmclasses+1.
    !------------------------------------------------------------
    !   Note that all local variables used in the network
    !   simulation need to be double precision.
    !------------------------------------------------------------
    !
    !   Arguments (I=input, O=output):
    !
    !   Ym:   (I) Array of brightness temperatures for current scene
    !   kchan:(I) valid channel mask for current scene
    !   igeo: (I) Geographic land surface class (optional and not used;
    !             provided only for back-compatibility; not needed operational)
    !
    !   iclassatm:   (O) Index of selected atmosphere class
    !
    !   U_Classatm: Logical unit number of network training file
    !   F_Classatm: Path of network training file
    !
    !   profile: Array of temperatures to be used in network
    !   simulation
    !
    !   The network has two layers.  Layer 1 has input weights given in
    !   array 'inpwts' and and a bias vector given in array 'bias1'.
    !   Layer 2 has layer weights given in array 'laywts'.  The output of
    !   layer 1 is given in array a1, and the output of layer 2 (i.e, the
    !   class of profile as determined by the network) is given in scalar
    !   iclassatm.
    !
    USE ReadStdInputs
    USE OSSMWmoduleSubs, ONLY : getOSSselMW

    !---Arguments
    !---xlat,xlon,time,igeo are not used, and are only for interface
    !---uniformity among classatm versions
    REAL,    DIMENSION(:),           INTENT(IN)    :: Ym
    LOGICAL, DIMENSION(:),           INTENT(IN)    :: kchan
    REAL,                  OPTIONAL, INTENT(IN)    :: xlat,xlon
    REAL,                  OPTIONAL, INTENT(IN)    :: eia
    INTEGER, DIMENSION(6), OPTIONAL, INTENT(IN)    :: time
    INTEGER,               OPTIONAL, INTENT(IN)    :: igeo
    INTEGER,                         INTENT(INOUT) :: iclassatm
    !---Local variables
    TYPE(GenControl_t) :: genControl
    TYPE(MWdefin_t), SAVE :: MWdefin
    TYPE(IRdefin_t) :: IRdefin
    TYPE(AtmosClassData_t) :: classDataIn
    INTEGER :: atmosClassMode
    INTEGER :: i
  !
    if (first) then
       classDataIn%disable = .false.
       call rdwts(F_Classatm,classDataIn)
       call putAtmosClassData(classDataIn)
       genControl%atmosClassConfig=' '
       CALL getOSSselMW(F_osscoefs,nchanOSS=MWdefin%nChan)
       allocate(MWdefin%frq(MWdefin%nchan))
       CALL getOSSselMW(F_osscoefs,cFreq=MWdefin%frq(1:MWdefin%nChan))
       atmosClassMode=MWmode
       call atmosClassInit(genControl,MWdefin,IRdefin,atmosClassMode)
       allocate(QC(MWdefin%nChan))
       first=.false.
    endif

    where (kchan(1:MWdefin%nChan))
      QC=0
    elsewhere
      QC=1
    end where
    call classAtmNN(Ym,QC,xlat,xlon,eia,time,igeo,iclassatm)

  END SUBROUTINE classatm

END MODULE ClassAtmMod
