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

subroutine SetXPtrR(nTskin,NTemp,Ngas,NCldLiq,NCldIce,molID,nmol,NG,IR,NR)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   SetXPtrR
!
! PURPOSE:
!
!   Set the retrieval vector structure based on input sizes of the components, 
!   with consistency checks.
!
! SYNTAX:
!
!   CALL SetXPtrR(nTskin, NTemp, Ngas, NCldLiq, NCldIce, molID, nmol, 
!      NG, IR, NR)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   nTskin   INTEGER             Number of elements for Tskin
!   NTemp    INTEGER             Number of elements for temperature
!   Ngas     INTEGER             Number of elements for molecular 
!                                concentrations 
!   NCldLiq  INTEGER             Number of elements for liquid cloud
!   NCldIce  INTEGER             Number of elements for ice cloud
!   molID    INTEGER             List of IDs of relevant molecular species
!   nmol     INTEGER             Number of molecular species
!   NG       TYPE(STATEINDEX_T)  Number of elements for sections of 
!                                geophysical state vector 
!   
!   INPUTS/OUTPUTS:
!   
!   IR       TYPE(STATEINDEX_T)  Starting indices for sections of retrieval 
!                                state vector 
!   NR       TYPE(STATEINDEX_T)  Number of retrieved elements for sections of 
!                                retrieval state vector
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

  USE StateIndexModule
  implicit none
  !---In/Out variables
  INTEGER,                   INTENT(IN)    :: nTskin
  INTEGER,                   INTENT(IN)    :: NTemp
  INTEGER, DIMENSION(maxMol),INTENT(IN)    :: ngas
  INTEGER,                   INTENT(IN)    :: NCldLiq,NCldIce
  INTEGER, DIMENSION(nmol),  INTENT(IN)    :: molID
  INTEGER,                   INTENT(IN)    :: nmol
  TYPE(StateIndex_t),        INTENT(IN)    :: NG
  TYPE(StateIndex_t),        INTENT(INOUT) :: IR,NR
  !--Local variables
  INTEGER            :: i
  
  NR = initLengths()
  NR%temp         = NTemp
  if (NR%temp .GT. NG%temp) then
     print*,"err[SetXPtrR]: Retr versus Geo: NR%Temp exceedes NG%temp"
     call errorHalt(1)
  endif
  NR%Tskin = nTskin
  if (NR%temp .GT. 0) NR%Tskin = 1  ! Always retrieve Tskin if retrieving T(p)
  if (NR%Tskin .GT. NG%Tskin) then
     print*,"err[SetXPtrR]: Retr versus Geo: NR%Tskin exceedes NG%Tskin"
     call errorHalt(1)
  endif
  NR%Psfc         = 0 

  do i=1,nmol
     NR%mol(i)=ngas(molID(i))
     if (NR%mol(i) .GT. NG%mol(i)) then
        print*, 'MolID=',molID(i)
        print*,"err[SetXPtrR]: Retr versus Geo: NR%mol exceedes NG%mol"
        call errorHalt(1)
     endif
  enddo
  NR%Cldliq       = NCldLiq
  if (NR%Cldliq .NE. NG%Cldliq  .AND.  NR%Cldliq .NE. 0) then
     ! NR%CldLiq must equal NG%CldLiq if cloud liquid is to be retrieved.
     print*,"err[SetXPtrR]: Retr versus Geo: CldLiq is retrieved, "
     print*,"but NR%CldLiq .ne. NG%CldLiq"
     print*,'NR%CldLiq,NG%CldLiq:',NR%CldLiq,NG%CldLiq
     call errorHalt(1)
  endif  
  NR%CldIce       = NCldIce
  if (NR%CldIce .NE. NG%CldIce  .AND.  NR%CldIce .NE. 0) then
     ! NR%CldIce must equal NG%CldIce if cloud ice is to be retrieved.
     print*,"err[SetXPtrR]: Retr versus Geo: CldIce is retrieved, "
     print*,"but NR%CldIce .ne. NG%CldIce"
     print*,'NR%CldIce,NG%CldIce:',NR%CldIce,NG%CldIce
     call errorHalt(1)
  endif  
  IR              = genIndices(NR)
  
  return
end subroutine SetXPtrR
