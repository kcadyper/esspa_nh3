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

module StateIndexModule

! <f90Module>***********************************************************
!
! NAME:
!
!   StateIndexModule
!
! PURPOSE:
!
!   Module to define the state vector index structure and associated operations
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

  !--------------------------------------------------------------------------
  !"$Id$"
  !--------------------------------------------------------------------------
  ! Module to define the state vector index structure and associated
  ! operations.
  !
  ! Data type name: StateIndex_t
  !      temp: temperature
  !     tskin: skin temperature
  !      psfc: surface pressure
  !       mol: molecular species
  !    cldliq: Liquid cloud parameters 
  !            (cloud top, cloud thickness, cloud fraction)
  !    cldice: Ice cloud parameters
  !      wind: surface wind
  !
  ! Operation name:
  !  initLengths: subroutine to initialize all the lengths of each 
  !               individual parameter in the state vector to be zero.
  !   genIndices: 
  !             function to generate indices with the given vector length 
  !             for each variable.
  !   getVectorLength: 
  !             function to return the total vector length given the
  !             collective variable lengths.
  !   getAtmosVectorLength: 
  !             function to return the total vector length given the
  !             collective variable lengths, accounting for atmosphere only
  !   getNMol:
  !             The overloading interface that is consistes of two 
  !             functions. One is to get the number of molecules by
  !             from a state vector (or profile). The other is to get
  !             the number of molecules by a given molecule ID vector.
  !   whereH2O:
  !             returns the H2O ID, stops if not found
  !   searchID: 
  !             searches the given ID in an ID array and returns the 
  !             index of the ID in the ID array if found; Otherwise, 
  !             returns -1.
  !   compressMol: 
  !             returns the compressed structure (length or index) 
  !             for state vector according to the given compressed MolID
  !   printIndex: 
  !             subroutine to print out tabulated indices and lengths
  !             of all the variables.
  !   molecStringToIndex:
  !             Convert a string molecule ID to the index where the 
  !             molecule is stored in the standard convention
  !   charToMolID:
  !             Convert list of string molecule IDs to an integer vector of 
  !             indices where molecules are stored in the standard convention
  !
  !    Variable       Type             Description
  !    --------       ----             -----------  
  !    NV             StateIndex_t     parameter lengths
  !    IX             StateIndex_t     parameter indices
  !
  ! Copyright, AER, Inc., 2004
  ! Developed at AER, Inc., January, 2004
  !--------------------------------------------------------------------------
  use constants, only: MISSING_CHAR
  implicit none
  integer, parameter :: maxMol=35       !maximum number of molecule species
  integer, parameter :: nullMolID=-1    !default value for init molID
  integer, parameter :: MAXLEV=101      !maximum number of vertical grids
  integer, parameter :: MAXPAR=(maxMol+3)*MAXLEV+4 
  !maximum number of paramters in 
  !the state vector. MAXLEV is 
  !not a hard limit on the number of 
  !levels. It is used only as a basis 
  !for static allocation of memory.
  integer, parameter :: &
       idH2O  =1, idCO2   =2, idO3  =3, idN2O    =4, idCO    =5, &
       idCH4  =6, idO2    =7, idNO  =8, idSO2    =9, idNO2  =10, &
       idNH3 =11, idHNO3 =12, idOH =13, idHF    =14, idHCL  =15, &
       idHBR =16, idHI   =17, idCLO=18, idOCS   =19, idH2CO =20, &
       idHOCL=21, idN2   =22, idHCN=23, idCH3CL =24, idH2O2 =25, &
       idC2H2=26, idC2H6 =27, idPH3=28, idCOF2  =29, idSF6  =30, &
       idH2S =31, idHCOOH=32, idHO2=33, idATMOXY=34, idCLNO3=35
  integer, private :: i                 !loop index
  logical, private :: dbg = .true.
  integer, parameter, private :: mxLength = 256
  private getNMolfromProf, getNMolfromMolID, getMolIDfromExpandedMolN

  TYPE StateFlags_t
     logical :: temp
     logical :: Tskin
     logical, dimension(maxMol) :: mol
     logical :: pSfc
     logical :: wind
     logical :: cloudLiq
     logical :: cloudIce
     logical :: emMW
     logical :: emRfIR
  END TYPE StateFlags_t

!!$  type StateIndex_t
!!$     integer :: temp, tskin, psfc, mol(maxMol), cldliq, cldice, wind
!!$  end type StateIndex_t

  TYPE StateIndex_t
     integer :: temp
     integer :: Tskin
     integer, dimension(maxMol) :: mol
     integer :: pSfc
     integer :: wind
!!$     integer :: cloudLiq
!!$     integer :: cloudIce
     integer :: cldLiq
     integer :: cldIce
     integer :: emMW
     integer :: emRfIR
  END TYPE StateIndex_t

  interface getNMol
     module procedure getNMolfromProf
     module procedure getNMolfromMolID
  end interface

  interface initFileID
     module procedure initFileID
  end interface

contains
  function initLengths(Ntemp,Ntskin,Npsfc,NmolLev,&
       Ncldliq,Ncldice,Nwind,NemMW,NemRfIR,Natmos)

!<f90Function>**********************************************************
!
! NAME:
!
!   initLengths
!
! PURPOSE:
!
!   Initialize all the lengths of each individual parameter in the state vector 
!   to zero.
!
! SYNTAX:
!
!   Results=initLengths(Ntemp, Ntskin, Npsfc, NmolLev, Ncldliq, 
!      Ncldice, Nwind, NemMW, NemRfIR, Natmos)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   Ntemp*       INTEGER             Number of elements for temperature
!   Ntskin*      INTEGER             Number of elements for Tskin
!   Npsfc*       INTEGER             Number of elements for surface pressure
!   NmolLev*     INTEGER             Number of elements for molecular 
!                                    concentrations 
!   Ncldliq*     INTEGER             Number of elements for liquid cloud
!   Ncldice*     INTEGER             Number of elements for ice cloud
!   Nwind*       INTEGER             Number of elements for wind
!   NemMW*       INTEGER             Number of elements for MW emissivity
!   NemRfIR*     INTEGER             Number of elements for IR emissivity/reflectivity
!   Natmos*      StateIndex_t        Numbers for all atmosphere elements
!
!   * OPTIONAL
!
! RETURN:
!
!     TYPE(STATEINDEX_T)  
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

    integer,               intent(in), optional :: Ntemp,Ntskin,Npsfc
    integer,               intent(in), optional :: Ncldliq,Ncldice,Nwind
    integer, dimension(:), intent(in), optional :: NmolLev
    integer,               intent(in), optional :: NemMW,NemRfIR
    type(StateIndex_t),    intent(in), optional :: Natmos
    type(StateIndex_t)     :: initLengths
    initLengths%temp = 0
    initLengths%tskin = 0
    initLengths%psfc = 0
    initLengths%mol = 0
    initLengths%cldLiq = 0
    initLengths%cldIce = 0
    initLengths%wind = 0
    initLengths%emMW = 0
    initLengths%emRfIR = 0
    if (present(Ntemp))   initLengths%temp   = Ntemp
    if (present(Ntskin))  initLengths%tskin  = Ntskin
    if (present(Npsfc))   initLengths%psfc   = Npsfc
    if (present(NmolLev)) initLengths%mol(1:size(NmolLev)) = NmolLev
    if (present(Ncldliq)) initLengths%cldLiq  = Ncldliq
    if (present(Ncldice)) initLengths%cldIce = Ncldice
    if (present(Nwind))   initLengths%wind   = Nwind
    if (present(NemMW))   initLengths%emMW   = NemMW
    if (present(NemRfIR)) initLengths%emRfIR = NemRfIR
    if (present(Natmos)) then
       initLengths%temp   = Natmos%temp
       initLengths%tskin  = Natmos%tskin
       initLengths%psfc   = Natmos%psfc
       initLengths%mol    = Natmos%mol
       initLengths%cldLiq = Natmos%cldLiq
       initLengths%cldIce = Natmos%cldIce
       initLengths%wind   = Natmos%wind
    endif
  end function initLengths

  function genIndices(NV)

!<f90Function>**********************************************************
!
! NAME:
!
!   genIndices
!
! PURPOSE:
!
!   Generate indices with given vector length for each variable.
!
! SYNTAX:
!
!   Results=genIndices(NV)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   NV          TYPE(STATEINDEX_T)  Number of elements for sections of state 
!                                   vector
!
!   * OPTIONAL
!
! RETURN:
!
!     TYPE(STATEINDEX_T)  
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

    type(StateIndex_t), intent(in) :: NV
    type(StateIndex_t)  :: genIndices
    genIndices%temp = 1
    genIndices%tskin = genIndices%temp + NV%temp
    genIndices%psfc = genIndices%tskin + NV%tskin
    genIndices%mol(1) = genIndices%psfc + NV%psfc
    do i = 2, maxMol
       genIndices%mol(i) = genIndices%mol(i-1) + NV%mol(i-1)
    enddo
    genIndices%cldLiq = genIndices%mol(maxMol) + NV%mol(maxMol)
    genIndices%cldIce = genIndices%cldLiq + NV%cldLiq
    genIndices%wind = genIndices%cldIce + NV%cldIce
    genIndices%emMW = genIndices%wind + NV%wind
    genIndices%emRfIR = genIndices%emMW + NV%emMW
  end function genIndices

  function getVectorLength(NV)

!<f90Function>**********************************************************
!
! NAME:
!
!   getVectorLength
!
! PURPOSE:
!
!   Return the length of state vector
!
! SYNTAX:
!
!   Results=getVectorLength(NV)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   NV               TYPE(STATEINDEX_T)  Number of elements for sections of 
!                                        state vector
!
!   * OPTIONAL
!
! RETURN:
!
!     INTEGER             
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

    type(StateIndex_t), intent(in) :: NV
    integer  :: getVectorLength

    getVectorLength = NV%temp + NV%tskin + NV%psfc + sum(NV%mol) + &
         NV%cldLiq + NV%cldIce + NV%wind + NV%emMW + NV%emRfIR
  end function getVectorLength

  function getAtmosVectorLength(NV)

!<f90Function>**********************************************************
!
! NAME:
!
!   getAtmosVectorLength
!
! PURPOSE:
!
!   Return the length of state vector
!
! SYNTAX:
!
!   Results=getAtmosVectorLength(NV)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   NV               TYPE(STATEINDEX_T)  Number of elements for sections of 
!                                        state vector
!
!   * OPTIONAL
!
! RETURN:
!
!     INTEGER             
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

    type(StateIndex_t), intent(in) :: NV
    integer  :: getAtmosVectorLength

    getAtmosVectorLength = NV%temp + NV%tskin + NV%psfc + sum(NV%mol) + &
         NV%cldLiq + NV%cldIce + NV%wind
  end function getAtmosVectorLength

  function getNMolfromProf(NV)

!<f90Function>**********************************************************
!
! NAME:
!
!   getNMolfromProf
!
! PURPOSE:
!
!   Get the number of molecules from a state vector (or profile) structure.
!
! SYNTAX:
!
!   Results=getNMolfromProf(NV)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   NV               TYPE(STATEINDEX_T)  Number of elements for sections of 
!                                        state vector
!
!   * OPTIONAL
!
! RETURN:
!
!     INTEGER             
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

    type(StateIndex_t), intent(in) :: NV
    integer :: getNMolfromProf
    integer :: imol
    getNMolfromProf = 0
    do imol = 1, size(NV%mol)
       if (NV%mol(imol) > 0) then
          getNMolfromProf = getNMolfromProf + 1
       endif
    enddo
  end function getNMolfromProf

  function getNMolfromMolID(MolID)

!<f90Function>**********************************************************
!
! NAME:
!
!   getNMolfromMolID
!
! PURPOSE:
!
!   Get the number of molecules from a given molecule ID vector
!
! SYNTAX:
!
!   Results=getNMolfromMolID(MolID)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   MolID             INTEGER  List of IDs of relevant molecular species
!
!   * OPTIONAL
!
! RETURN:
!
!     INTEGER  
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

    integer, dimension(:), intent(in) :: MolID
    integer :: getNMolfromMolID
    integer :: i
    getNMolfromMolID = 0
    do i = 1, size(MolID)
       if (MolID(i) <= 0 .or. MolID(i) > maxMol) exit
       getNMolfromMolID = getNMolfromMolID + 1
    enddo
  end function getNMolfromMolID

  subroutine getMolIDfromExpandedMolN(NV,molID,nMol)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   getMolIDfromExpandedMolN
!
! PURPOSE:
!
!   Create a molecule ID vector from a state vector structure.
!   This works properly only if the  molecule sizes in NV are
!   in expanded form
!
! SYNTAX:
!
!   CALL getMolIDfromExpandedMolN(NV, molID, nMol)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   NV     TYPE(STATEINDEX_T)  Number of elements for sections of state 
!                              vector 
!   
!   INPUTS/OUTPUTS:
!   
!   molID  INTEGER             List of IDs of relevant molecular species
!   nMol*  INTEGER             Number of molecular species
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    type(StateIndex_t),    intent(in)              :: NV
    integer, dimension(:), intent(inout)           :: molID
    integer,               intent(inout), optional :: nMol
    type(StateIndex_t) :: IX
    integer :: imol, cntMol
    IX = genIndices(NV)
    cntMol = 0
    do imol = 1, size(NV%mol)
       if (NV%mol(imol) > 0) then
          cntMol = cntMol + 1
          MolID(cntMol) = imol
       endif
    enddo
    if (present(nMol)) nMol = cntMol
    return
  end subroutine getMolIDfromExpandedMolN

  
  ! -------------------------------------------------
  ! Returns the index of H2O if found, otherwise -1
  ! The caller should check for the returned result.
  ! -------------------------------------------------
  function whereH2O(molID)

!<f90Function>**********************************************************
!
! NAME:
!
!   whereH2O
!
! PURPOSE:
!
!   Returns the index of H2O if found, otherwise -1.
!
! SYNTAX:
!
!   Results=whereH2O(molID)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   molID     INTEGER  List of IDs of relevant molecular species
!
!   * OPTIONAL
!
! RETURN:
!
!     INTEGER  
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

    integer, dimension(:), intent(in) :: molID
    integer  :: whereH2O
    do i=1,size(molID)
       if (molID(i) == idH2O) then
          whereH2O=i
          return
       endif
    enddo
    whereH2O = -1  ! H2O not found
    print *,'Warning[StateIndexModule::whereH2O]:  H2O not found in molID'
  end function whereH2O

  function searchID(IDarr,ID)

!<f90Function>**********************************************************
!
! NAME:
!
!   searchID
!
! PURPOSE:
!
!   Searches for a certain ID in an array of IDs and returns the index of the ID 
!   in the ID array if found; Otherwise, returns -1.
!
! SYNTAX:
!
!   Results=searchID(IDarr, ID)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   IDarr     INTEGER  Array of identifier indices
!   ID        INTEGER  Identifier index
!
!   * OPTIONAL
!
! RETURN:
!
!     INTEGER  
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

    integer, dimension(:), intent(in) :: IDarr
    integer,               intent(in) :: ID
    integer  :: searchID
    searchID = -1
    do i=1,size(IDarr)
       if (IDarr(i) == ID) then 
          searchID = i
          return
       endif
    enddo
    return
  end function searchID

  function compressMol(Nin,MolID)

!<f90Function>**********************************************************
!
! NAME:
!
!   compressMol
!
! PURPOSE:
!
!   Converts the molecule-related portion of a state vector size structure so 
!   that the active molecules are moved from the order specified in MolID to the 
!   front.
!
! SYNTAX:
!
!   Results=compressMol(Nin, MolID)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   Nin          TYPE(STATEINDEX_T)  Number of elements for sections of state 
!                                    vector at input 
!   MolID        INTEGER             List of IDs of relevant molecular 
!                                    species
!
!   * OPTIONAL
!
! RETURN:
!
!     TYPE(STATEINDEX_T)  
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

    type(StateIndex_t),    intent(in) :: Nin
    integer, dimension(:), intent(in) :: MolID
    type(StateIndex_t)     :: compressMol
    integer :: NMol
    NMol = getNMol(MolID)
    compressMol = Nin
    compressMol%mol(1:maxMol) = 0
    compressMol%mol(1:NMol) = Nin%mol(MolID(1:NMol))
  end function compressMol

  function expandMol(Nin,MolID)

!<f90Function>**********************************************************
!
! NAME:
!
!   expandMol
!
! PURPOSE:
!
!   Converts the molecule-related portion of a state vector size structure so 
!   that the active molecules are moved from the front, in the order specified 
!   in MolID.
!
! SYNTAX:
!
!   Results=expandMol(Nin, MolID)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   Nin        TYPE(STATEINDEX_T)  Number of elements for sections of state 
!                                  vector at input 
!   MolID      INTEGER             List of IDs of relevant molecular species
!
!   * OPTIONAL
!
! RETURN:
!
!     TYPE(STATEINDEX_T)  
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

    type(StateIndex_t),    intent(in) :: Nin
    integer, dimension(:), intent(in) :: MolID
    type(StateIndex_t)     :: expandMol
    integer :: NMol
    NMol = getNMol(Nin)
    if (size(MolID) < NMol) then
       print *,"err[StateModuleIndex::expandMol]:  #Mol exceedes size(MolID)."
       call errorHalt(1)
    end if
    expandMol = Nin
    expandMol%mol(1:maxMol) = 0
    expandMol%mol(MolID(1:NMol)) = Nin%mol(1:NMol)
  end function expandMol

  subroutine genMapFrom(IG,NG,IGF,NGF,MolID,map,mapSize,xFileSize)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   genMapFrom
!
! PURPOSE:
!
!   Generate state vector size and index structures listing molecules according 
!   to the content of a source file (default), optionally reduced according to 
!   user selections (if MolID is present). Generate an index map for remapping a 
!   state vector, according to the size structures of the source and the output.
!
! SYNTAX:
!
!   CALL genMapFrom(IG, NG, IGF, NGF, MolID, map, mapSize, xFileSize)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   IGF        TYPE(STATEINDEX_T)  Starting indices for sections of 
!                                  geophysical state vector in file 
!   NGF        TYPE(STATEINDEX_T)  Number of elements for sections of state 
!                                  vector in file 
!   
!   INPUTS/OUTPUTS:
!   
!   IG         TYPE(STATEINDEX_T)  Starting indices for sections of 
!                                  geophysical state vector 
!   NG         TYPE(STATEINDEX_T)  Number of elements for sections of 
!                                  geophysical state vector 
!   MolID*     INTEGER             List of IDs of relevant molecular species
!   map        INTEGER             List of indices to use in remapping data
!   mapSize    INTEGER             Map dimension
!   xFileSize  INTEGER             number of constituents in NGF
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    type(StateIndex_t),    intent(in)               :: IGF,NGF
    type(StateIndex_t),    intent(inout)            :: IG,NG
    integer, dimension(:), intent(inout), optional  :: MolID
    integer, dimension(:), intent(inout)            :: map
    integer,               intent(inout)            :: mapSize
    integer,               intent(inout)            :: xFileSize
    integer, dimension(:), allocatable :: MolIDF
    integer, dimension(maxMol) :: MolIDloc
    integer :: NMol, NMolF, i, imol, index, NParG

    xFileSize = getAtmosVectorLength(NGF)
    NMolF = getNMol(NGF)
    if (allocated(MolIDF)) deallocate(MolIDF)
    allocate(MolIDF(NMolF))
    call getMolIDfromExpandedMolN(NGF,MolIDF)
    MolIDloc=nullMolID

    if (present(MolID)) then
       if (any(MolID /= nullMolID)) then
          !detect if user's molID's are in the file
          do imol = 1, size(MolID)
             if (MolID(imol) <= 0 .or. MolID(imol) > maxMol) exit
             index = searchID(MolIDF,MolID(imol))
             if (index < 0) then
                print *,'err[StateIndexModule::genMapFrom]: ',&
                     ' MolID(',MolID(imol),') not found'
                call errorHalt(1)
             endif
          enddo
       else
          MolID(1:NMolF) = MolIDF(1:NMolF)
       endif
       MolIDloc(1:NMolF)=MolID(1:NMolF)
    else
       MolIDloc(1:NMolF) = MolIDF(1:NMolF)
    endif
    ! build mapFrom from user requested MolID
    NG = compressMol(NGF,MolIDloc)
    IG = genIndices(NG)
    deallocate(MolIDF)
    NParG = getAtmosVectorLength(NG)
    if (NParG > size(map)) then
       print *,'err[StateIndexModule::genMapFrom]:  map size exceeded.'
       call errorHalt(1)
    end if
    mapSize = NParG
    NMol = getNMol(NG)

    map(IG%temp:IG%temp+NG%temp-1) = &
         (/(i,i=IGF%temp,IGF%temp+NGF%temp-1)/)
    map(IG%tskin:IG%tskin+NG%tskin-1) = &
         (/(i,i=IGF%tskin,IGF%tskin+NGF%tskin-1)/)
    map(IG%psfc:IG%psfc+NG%psfc-1) = &
         (/(i,i=IGF%psfc,IGF%psfc+NGF%psfc-1)/)
    do imol=1,NMol
       map(IG%mol(imol):IG%mol(imol)+NG%mol(imol)-1) = &
            (/(i, i=IGF%mol(MolIDloc(imol)), &
            IGF%mol(MolIDloc(imol))+NGF%mol(MolIDloc(imol))-1)/) 
    enddo
    map(IG%cldLiq:IG%cldLiq+NG%cldLiq-1) = &
         (/(i,i=IGF%cldLiq,IGF%cldLiq+NGF%cldLiq-1)/)
    map(IG%cldIce:IG%cldIce+NG%cldIce-1) = &
         (/(i,i=IGF%cldIce,IGF%cldIce+NGF%cldIce-1)/)
    map(IG%wind:IG%wind+NG%wind-1) = &
         (/(i,i=IGF%wind,IGF%wind+NGF%wind-1)/)
  end subroutine genMapFrom

  subroutine initFileID(ncid,fidIndex,stateIndexEnabled,IG,NG,fid)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   initFileID
!
! PURPOSE:
!
!   Keep list of IDs of all open netCDF files, for keeping track of the 
!   properties of each file.
!
! SYNTAX:
!
!   CALL initFileID(ncid, fidIndex, stateIndexEnabled, IG, NG, fid)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   ncid               INTEGER             ID for ncdf file 
!   IG*                TYPE(STATEINDEX_T)  Starting indices for sections of 
!                                          geophysical state vector 
!   NG*                TYPE(STATEINDEX_T)  Number of elements for sections of 
!                                          geophysical state vector 
!   
!   INPUTS/OUTPUTS:
!   
!   fidIndex           INTEGER             File identifier index
!   stateIndexEnabled  LOGICAL             Flag to enable/disable stateIndex
!   fid                INTEGER             File identifier index list
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    integer,               intent(in)           :: ncid
    integer,               intent(inout)        :: fidIndex
    logical, dimension(:), intent(inout)        :: stateIndexEnabled
    type(StateIndex_t),    intent(in), optional :: IG,NG
    integer, dimension(:), intent(inout)        :: fid
    integer :: index

    index = searchID(fid,ncid)
    if (index < 0) then       !newly opened: create the map for new ncid
       fidIndex = fidIndex + 1
    else                      !previously opened: re-use the slot
       fidIndex = index
    endif
    fid(fidIndex) = ncid

    if (present(IG) .and. present(NG)) then
       stateIndexEnabled(fidIndex) = .TRUE.
    else
       stateIndexEnabled(fidIndex) = .FALSE.
    endif
  end subroutine initFileID

  subroutine printIndex(NV)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   printIndex
!
! PURPOSE:
!
!   Print state vector size and index structures
!
! SYNTAX:
!
!   CALL printIndex(NV)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   NV  TYPE(STATEINDEX_T)  Number of elements for sections of state vector
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    type(StateIndex_t), intent(in) :: NV
    type(StateIndex_t) :: IX
    IX = genIndices(NV)
    print '(a8,2a7)','--------','-------','-------'
    print '(a8,2a7)','Var_Name','  Index',' Length'
    print '(a8,2a7)','--------','-------','-------'
    print '(a8,2i7)','Temp    ', IX%temp,   NV%temp
    print '(a8,2i7)','Tskin   ', IX%tskin,  NV%tskin
    print '(a8,2i7)','Psfc    ', IX%psfc,   NV%psfc
    do i = 1, maxMol
       print '(a4,i2,a2,2i7)','mol(', i,') ',IX%mol(i), NV%mol(i)
    enddo
    print '(a8,2i7)','CldLiq  ', IX%cldLiq,  NV%cldLiq
    print '(a8,2i7)','CldIce  ', IX%cldIce, NV%cldIce
    print '(a8,2i7)','sfcWind ', IX%wind,   NV%wind
    print '(a8,2i7)','MW emis ', IX%emMW,   NV%emMW
    print '(a8,2i7)','IR EmRf ', IX%emRfIR,   NV%emRfIR
    print '(a8,2a7)','--------','-------','-------'
    print '(a8,2a7)','--------','-------','-------'
    print '(a17,i5)', 'Total #Elements: ', getVectorLength(NV)
  end subroutine printIndex

  function molecStringToIndex(molecString)

!<f90Function>**********************************************************
!
! NAME:
!
!   molecStringToIndex
!
! PURPOSE:
!
!   Convert a character string that identifies a molecule to the index where the 
!   molecule is stored in the standard convention
!
! SYNTAX:
!
!   Results=molecStringToIndex(molecString)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   molecString         CHAR     Charcter string that identifies a molecular 
!                                species
!
!   * OPTIONAL
!
! RETURN:
!
!     INTEGER  
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

  ! Convert a character string that identifies a molecule to the 
  ! index where the molecule is stored in the standard convention
  ! 0 is returned if no match is found

    character(len=*),    intent(in)  :: molecString
    integer                          :: molecStringToIndex

    integer, parameter               :: maxLenMolecString=12
    integer                          :: ii,icoff,luse
    character(len=maxLenMolecString) :: slocal

    ! Convert to all upper case for comparisons
    icoff=ichar('A')-ichar('a')
    slocal=repeat(' ',maxLenMolecString)
    lUse=min(len(molecString),maxLenMolecString)
    slocal(1:lUse)=molecString(1:lUse)
    do i=1,lUse
      if (slocal(i:i) >= 'a' .and. slocal(i:i) <= 'z') &
        slocal(i:i) = char(ichar(slocal(i:i))+icoff)
    enddo

    ! Find the right index
    select case(trim(slocal))
    case ('H2O')
      molecStringToIndex=1
    case ('CO2')
      molecStringToIndex=2
    case ('O3')
      molecStringToIndex=3
    case ('N2O')
      molecStringToIndex=4
    case ('CO')
      molecStringToIndex=5
    case ('CH4')
      molecStringToIndex=6
    case ('O2')
      molecStringToIndex=7
    case ('NO')
      molecStringToIndex=8
    case ('SO2')
      molecStringToIndex=9
    case ('NO2')
      molecStringToIndex=10
    case ('NH3')
      molecStringToIndex=11
    case ('HNO3')
      molecStringToIndex=12
    case ('OH')
      molecStringToIndex=13
    case ('HF')
      molecStringToIndex=14
    case ('HCL')
      molecStringToIndex=15
    case ('HBR')
      molecStringToIndex=16
    case ('HI')
      molecStringToIndex=17
    case ('CLO')
      molecStringToIndex=18
    case ('OCS')
      molecStringToIndex=19
    case ('H2CO')
      molecStringToIndex=20
    case ('HOCL')
      molecStringToIndex=21
    case ('N2')
      molecStringToIndex=22
    case ('HCN')
      molecStringToIndex=23
    case ('CH3CL')
      molecStringToIndex=24
    case ('H2O2')
      molecStringToIndex=25
    case ('C2H2')
      molecStringToIndex=26
    case ('C2H6')
      molecStringToIndex=27
    case ('PH3')
      molecStringToIndex=28
    case ('COF2')
      molecStringToIndex=29
    case ('SF6')
      molecStringToIndex=30
    case ('H2S')
      molecStringToIndex=31
    case ('HCOOH')
      molecStringToIndex=32
    case ('HO2')
      molecStringToIndex=33
    case ('ATMOXY')
      molecStringToIndex=34
    case ('CLNO3')
      molecStringToIndex=35
    case default
      molecStringToIndex=0
    end select

  end function molecStringToIndex

  subroutine charToMolID(molOn,molID)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   charToMolID
!
! PURPOSE:
!
!   Convert list of string molecule IDs to an integer vector of indices where 
!   molecules are stored in the standard convention
!
! SYNTAX:
!
!   CALL charToMolID(molOn, molID)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   molOn  CHAR     List of string molecular IDs
!   
!   INPUTS/OUTPUTS:
!   
!   molID  INTEGER  List of IDs of relevant molecular species
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

  ! Convert list of string molecule IDs to an integer vector of
  ! indices where molecules are stored in the standard convention
    character(len=*), dimension(:), intent(in)    :: molOn
    integer,          dimension(:), intent(inout) :: molID
    integer :: imol, cntMol

    molID=nullMolID
    cntMol=0
    do imol=1,size(molOn)
      if (molOn(imol) /= MISSING_CHAR) then
        cntMol=cntMol+1
        if (cntMol > size(molID)) then
           print *,'err:[StateIndexModule::charToMolID] cntMol > size(molID):',&
              cntMol,size(molID)
           call errorHalt(1)
        endif
        molID(cntMol)=molecStringToIndex(molOn(imol))
        if (molID(cntMol) <= 0) then
           print *,'err:[StateIndexModule::charToMolID] '// &
             'Invalid molecule selected; '// &
             'imol,molOn(imol):',imol,molOn(imol)
           call errorHalt(1)
        endif
      endif
    enddo
    return
  end subroutine charToMolID

  logical function stateIndicesMatch(stateIndexA,stateIndexB,atmos,sfcMW,sfcIR)
    TYPE(StateIndex_t), intent(in) :: stateIndexA
    TYPE(StateIndex_t), intent(in) :: stateIndexB
    logical, intent(in), optional :: atmos
    logical, intent(in), optional :: sfcMW
    logical, intent(in), optional :: sfcIR
    !local
    character (len=mxLength) :: procName

    procName = ' [StateIndexModule::stateIndicesMatch]: '
    if (dbg) print *, trim(procName)//'starting ...'
    stateIndicesMatch = .TRUE.
    if (present(atmos)) then
       if (atmos) then
          if (stateIndexA%temp /= stateIndexB%temp) then
             print*,trim(procName)//'error - temp inconsistent', &
                stateIndexA%temp,stateIndexB%temp
             stateIndicesMatch = .FAlSE.
          end if
          if (stateIndexA%Tskin /= stateIndexB%Tskin) then
             print*,trim(procName)//'error - Tskin inconsistent', &
                stateIndexA%Tskin,stateIndexB%Tskin
             stateIndicesMatch = .FAlSE.
          end if
          if (stateIndexA%pSfc /= stateIndexB%pSfc) then
             print*,trim(procName)//'error - pSfc inconsistent', &
                stateIndexA%pSfc,stateIndexB%pSfc
             stateIndicesMatch = .FAlSE.
          end if
          if (any((stateIndexA%mol - stateIndexB%mol) /= 0)) then
             print*,trim(procName)//'error - mol inconsistent', &
                stateIndexA%mol,stateIndexB%mol
             stateIndicesMatch = .FAlSE.
          end if
          if (stateIndexA%cldLiq /= stateIndexB%cldLiq) then
             print*,trim(procName)//'error - cldLiq inconsistent', &
                stateIndexA%cldLiq,stateIndexB%cldLiq
             stateIndicesMatch = .FAlSE.
          end if
          if (stateIndexA%cldIce /= stateIndexB%cldIce) then
             print*,trim(procName)//'error - cldIce inconsistent', &
                stateIndexA%cldIce,stateIndexB%cldIce
             stateIndicesMatch = .FAlSE.
          end if
       end if
    end if
    if (present(sfcMW)) then
       if (sfcMW) then
          if (stateIndexA%emMW /= stateIndexB%emMW) then
             print*,trim(procName)//'error - emMW inconsistent', &
                stateIndexA%emMW,stateIndexB%emMW
             stateIndicesMatch = .FAlSE.
          end if
       end if
    end if
    if (present(sfcIR)) then
       if (sfcIR) then
          if (stateIndexA%emRfIR /= stateIndexB%emRfIR) then
             print*,trim(procName)//'error - emRfIR inconsistent', &
                stateIndexA%emRfIR,stateIndexB%emRfIR
             stateIndicesMatch = .FAlSE.
          end if
       end if
    end if

    if (dbg) print *, trim(procName)//'ending ...'
  end function stateIndicesMatch

end module StateIndexModule
