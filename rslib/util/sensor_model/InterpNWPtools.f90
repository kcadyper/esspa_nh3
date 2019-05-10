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

module interpNWPtools

! <f90Module>***********************************************************
!
! NAME:
!
!   interpNWPtools
!
! PURPOSE:
!
!   Module to read NWP and cluster data, average NWP to cluster centers and 
!   produce vectical profiles.
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

  !********************************************************************
  ! Module to read NWP and cluster data, 
  ! average NWP to cluster centers and produce vectical profiles
  !********************************************************************
  use StateIndexModule
  use ncdf_module
  use constants, ONLY : MISSING_REAL
  use MetFunctions, only : mrFromRH
  implicit none
  private
  public getNWPdim,getNWPdata,setStateIndex,&
       avgClusterProf,checkCldliq,computeH2oMixr
contains

!-------------------------------------------------------------------------------------------
subroutine getNWPdim(file,nSamples,nlpres,nlwv)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   getNWPdim
!
! PURPOSE:
!
!   Get NWP related dimensions from the input ncdf file
!
! SYNTAX:
!
!   CALL getNWPdim(file, nSamples, nlpres, nlwv)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   file      CHAR     file path
!   
!   INPUTS/OUTPUTS:
!   
!   nSamples* INTEGER  number of samples
!   nlpres*   INTEGER  pressure grid: Number of levels
!   nlwv*     INTEGER  Number of WV levels
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

  character (len=*), intent(in)              :: file
  integer,           intent(inout), optional :: nSamples,nlpres,nlwv
  integer :: ncid
  
  call openNcdfFile(ncid, file, status='old')
  if (present(nSamples)) nSamples=readNcdfDim(ncid,'nMembersTotal')
  if (present(nlpres))   nlpres=readNcdfDim(ncid,'vertical_pressure_grid')
  if (present(nlwv))     nlwv=readNcdfDim(ncid,'h2o_nLevels')
  call closeNcdfFile(ncid)
  return
end subroutine getNWPdim

!-------------------------------------------------------------------------------------------
subroutine getNWPdata(file,nSamples,nlpres,nlwv,presGrid,tair,tskin,psfc,temp,rh,pland,elev)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   getNWPdata
!
! PURPOSE:
!
!   Get NWP data from the input ncdf file
!
! SYNTAX:
!
!   CALL getNWPdata(file, nSamples, nlpres, nlwv, presGrid, tair, 
!      tskin, psfc, temp, rh, pland, elev)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   file      CHAR     file path
!   nSamples* INTEGER  number of samples
!   nlpres*   INTEGER  pressure grid: Number of levels
!   nlwv*     INTEGER  Number of WV levels
!   
!   INPUTS/OUTPUTS:
!   
!   presGrid* REAL     pressure grid
!   tair*     REAL     shelter-height temperature
!   tskin*    REAL     Surface skin temperature
!   psfc*     REAL     Surface pressure
!   temp*     REAL     Temperature
!   rh*       REAL     Relative humidity
!   pland*    REAL     Fraction of land (versus water) of the field of view
!   elev*     REAL     Elevation of surface
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

  character (len=*),    intent(in)              :: file
  integer,              intent(in),    optional :: nSamples,nlpres,nlwv
  real, dimension(:),   intent(inout), optional :: presGrid
  real, dimension(:),   intent(inout), optional :: tair,tskin,psfc,pland,elev
  real, dimension(:,:), intent(inout), optional :: temp,rh

  integer :: ncid
  integer(kind=2),dimension(:), allocatable     :: buf1D
  integer(kind=2),dimension(:,:), allocatable   :: buf2D

  integer :: nLines1,nSamples1,nlpres1,nlwv1
  real    :: missing_value,valid_min,valid_max,scale,offset
  
  call openNcdfFile(ncid, file, status='old')
  nSamples1=readNcdfDim(ncid,'nMembersTotal')
  if (present(nSamples)) then
     if (nSamples1 /= nSamples) then
        print*,'err[InterpNWPtools::GetNWPdata]: number of samples =',nSamples1
        return
     endif
  endif

  nlpres1=readNcdfDim(ncid,'vertical_pressure_grid')
  if (present(nlpres)) then
     if (nlpres1 /= nlpres) then
        print*,'err[InterpNWPtools::GetNWPdata]: number of press grids =',nlpres1
        return
     endif
  endif

  nlwv1=readNcdfDim(ncid,'h2o_nLevels')
  if (present(nlwv)) then
     if (nlwv1 /= nlwv) then
        print*,'err[InterpNWPtools::GetNWPdata]: number of rh grids =',nlwv1
        return
     endif
  endif

!----
  if (present(presGrid)) then
     if (allocated(buf1D)) deallocate(buf1D) 
     allocate(buf1D(nlpres1))
     call readNcdfData(ncid,buf1D,'PRESSURE_GRID')
     call readDataAtt(ncid,'PRESSURE_GRID',missing_value,valid_min,valid_max,scale,offset)
     presGrid(1:nlpres1)=buf1D*scale+offset
     where(presGrid == missing_value .or. presGrid < valid_min-scale .or. &
          presGrid > valid_max+scale)
     presGrid = MISSING_REAL
     endwhere
  endif
!----

  if (allocated(buf1D)) deallocate(buf1D) 
  allocate(buf1D(nSamples1))
  if (present(tskin)) then  
     call readNcdfData(ncid,buf1D,'TSKIN')
     call readDataAtt(ncid,'TSKIN',missing_value,valid_min,valid_max,scale,offset)
     tskin(1:nSamples1) = buf1D*scale+offset
     where(tskin == missing_value .or. tskin < valid_min-scale .or. &
          tskin > valid_max+scale)
        tskin = MISSING_REAL
     endwhere
  endif
!----

  if (present(tair)) then  
     call readNcdfData(ncid,buf1D,'SHELTER_TEMP')
     call readDataAtt(ncid,'SHELTER_TEMP',missing_value,valid_min,valid_max,scale,offset)
     tair(1:nSamples1) = buf1D*scale+offset
     where(tair == missing_value .or. tair < valid_min-scale .or. &
          tair > valid_max+scale)
        tair = MISSING_REAL
     endwhere
  endif
!----

  if (present(psfc)) then  
     call readNcdfData(ncid,buf1D,'PSFC')
     call readDataAtt(ncid,'PSFC',missing_value,valid_min,valid_max,scale,offset)
     psfc(1:nSamples1) = buf1D*scale+offset
     where(psfc == missing_value .or. psfc < valid_min-scale .or. &
          psfc > valid_max+scale)
        psfc = MISSING_REAL
     endwhere
  endif
!----

  if (present(pland)) then  
     call readNcdfData(ncid,buf1D,'pLand')
     call readDataAtt(ncid,'pLand',missing_value,valid_min,valid_max,scale,offset)
     pland(1:nSamples1) = buf1D*scale+offset
     where(pland == missing_value .or. pland < valid_min-scale .or. &
          pland > valid_max+scale)
        pland = MISSING_REAL
     endwhere
  endif
!----

  if (present(elev)) then  
     call readNcdfData(ncid,buf1D,'TERRAIN_HEIGHT')
     call readDataAtt(ncid,'TERRAIN_HEIGHT',missing_value,valid_min,valid_max,scale,offset)
     elev(1:nSamples1) = buf1D*scale+offset
     where(elev == missing_value .or. elev < valid_min-scale .or. &
          elev > valid_max+scale)
        elev = MISSING_REAL
     endwhere
  endif
!----
  
  if (present(temp)) then  
     if (allocated(buf2D)) deallocate(buf2D) 
     allocate(buf2D(nlpres1,nSamples1))
     call readNcdfData(ncid,buf2D,'TEMP')
     call readDataAtt(ncid,'TEMP',missing_value,valid_min,valid_max,scale,offset)
     temp(1:nlpres1,1:nSamples1) = buf2D*scale+offset
     where(temp == missing_value .or. temp < valid_min-scale .or. &
          temp > valid_max+scale)
        temp = MISSING_REAL
     endwhere
  endif
!----

  if (present(rh)) then  
     if (allocated(buf2D)) deallocate(buf2D) 
     allocate(buf2D(nlwv1,nSamples1))
     call readNcdfData(ncid,buf2D,'H2O')
     call readDataAtt(ncid,'H2O',missing_value,valid_min,valid_max,scale,offset)
     rh(1:nlwv1,1:nSamples1) = buf2D*scale+offset
     where(rh == missing_value .or. rh < valid_min-scale .or. &
          rh > valid_max+scale)
        rh = MISSING_REAL
     endwhere
  endif

  call closeNcdfFile(ncid)
  if (allocated(buf1D)) deallocate(buf1D) 
  if (allocated(buf2D)) deallocate(buf2D) 
  return
end subroutine getNWPdata
!-------------------------------------------------------------------------------------------

subroutine setStateIndex(nlpres,nlwv,nlev,IGsrc,NGsrc,molIDsrc,NParGsrc,&
       IGout,NGout,molIDout,NParGout,IGaux,NGaux,molIDaux,NParGaux)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   setStateIndex
!
! PURPOSE:
!
!   Set state vector lengths and indices
!
! SYNTAX:
!
!   CALL setStateIndex(nlpres, nlwv, nlev, IGsrc, NGsrc, molIDsrc, 
!      NParGsrc, IGout, NGout, molIDout, NParGout, IGaux, NGaux, 
!      molIDaux, NParGaux)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   nlpres    INTEGER             pressure grid: Number of levels
!   nlwv      INTEGER             Number of WV levels
!   nlev      INTEGER             Number of atmospheric levels
!   
!   INPUTS/OUTPUTS:
!   
!   IGsrc     TYPE(STATEINDEX_T)  Starting indices for sections of 
!                                 geophysical state vector of source data 
!   NGsrc     TYPE(STATEINDEX_T)  Number of elements for sections of 
!                                 geophysical state vector of source data 
!   molIDsrc  INTEGER             List of IDs of relevant molecular species 
!                                 of source data 
!   NParGsrc  INTEGER             Total number of elements in geophysical 
!                                 state vector of source data 
!   IGout     TYPE(STATEINDEX_T)  Starting indices for sections of 
!                                 geophysical state vector for output 
!   NGout     TYPE(STATEINDEX_T)  Number of elements for sections of 
!                                 geophysical state vector for output 
!   molIDout  INTEGER             List of IDs of relevant molecular species 
!                                 for output 
!   NParGout  INTEGER             Total number of elements in geophysical 
!                                 state vector for output 
!   IGaux*    TYPE(STATEINDEX_T)  Starting indices for sections of 
!                                 geophysical state vector for aux file 
!   NGaux*    TYPE(STATEINDEX_T)  Number of elements for sections of 
!                                 geophysical state vector for aux file 
!   molIDaux* INTEGER             List of IDs of relevant molecular species 
!                                 for aux file 
!   NParGaux* INTEGER             Total number of elements in geophysical 
!                                 state vector for aux file
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

  integer,               intent(in)              :: nlpres,nlwv,nlev
  type(StateIndex_t),    intent(inout)           :: IGsrc, NGsrc, IGout, NGout
  integer, dimension(:), intent(inout)           :: MolIDsrc,molIDout
  integer,               intent(inout)           :: NParGsrc,NParGout
  type(StateIndex_t),    intent(inout), optional :: IGaux, NGaux
  integer, dimension(:), intent(inout), optional :: molIDaux
  integer,               intent(inout), optional :: NParGaux
  integer                                         :: Ncldliq

!----
  molIDsrc(1)=idH2O
  NGsrc=initLengths(Ntemp=nlpres,Ntskin=1,Npsfc=1,NmolLev=(/nlwv/))
  IGsrc = genIndices(NGsrc)
  NParGsrc=getVectorLength(NGsrc)
!----
  if (NGsrc%cldliq>0) then
     Ncldliq=nlev
  else
     Ncldliq=3
  endif
  molIDout=molIDsrc
  NGout=initLengths(Ntemp=nlev,Ntskin=NGsrc%tskin,NmolLev=(/nlev/), &
       Npsfc=NGsrc%psfc,Ncldliq=Ncldliq,Ncldice=NGsrc%cldice,Nwind=NGsrc%wind)
  IGout = genIndices(NGout)
  NParGout=getVectorLength(NGout)
!----  
  if (present(molIDaux)) molIDaux=nullMolId
  if (present(NGaux))    NGaux = initLengths(Npsfc=1)
  if (present(IGaux)) then
     if (present(NGaux)) IGaux = genIndices(NGaux)
  endif
  if (present(NParGaux)) then
     if (present(NGaux)) NParGaux=getVectorLength(NGaux)
  endif
  return
end subroutine setStateIndex

!-------------------------------------------------------------------------------------------   
subroutine computeH2oMixr(temp,pres,rh,h2o)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   computeH2oMixr
!
! PURPOSE:
!
!   Loop over samples for conversion from relative humidity to mixing ratio.
!
! SYNTAX:
!
!   CALL computeH2oMixr(temp, pres, rh, h2o)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   temp  REAL  Temperature
!   pres  REAL  Pressure
!   rh    REAL  Relative humidity
!   
!   INPUTS/OUTPUTS:
!   
!   h2o   REAL  water vapor amount
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

  real, dimension(:,:), intent(in)    :: temp,rh
  real, dimension(:),   intent(in)    :: pres
  real, dimension(:,:), intent(inout) :: h2o
  integer                      :: nlpres,nlwv,i,j
  logical                      :: incrPresGrid

  nlpres=size(pres,1)
  nlwv=size(rh,1)
  incrPresGrid = (pres(1) < pres(nlpres))
  h2o=0
  do i=1,size(rh,2)
     if (incrPresGrid) then
        h2o(1:nlwv,i)= &
             h2oMixr(temp(nlpres-nlwv+1:nlpres,i),pres(nlpres-nlwv+1:nlpres),rh(1:nlwv,i),nlwv)
     else
        h2o(1:nlwv,i)= &
             h2oMixr(temp(1:nlwv,i),pres(1:nlwv),rh(1:nlwv,i),nlwv)
     endif
  enddo
  return
end subroutine computeH2oMixr

!-------------------------------------------------------------------------------------------
subroutine avgClusterProf(clustStart,clustNmem,memberList,nSamples,IG,NG,molID,&
     temp,h2o,tskin,psfc,pland,elev,prof,plandAvg,elevAvg,invertGrid,qc)
  integer,                 intent(in)              :: clustStart,clustNmem  !0-based

!<f90Subroutine>********************************************************
!
! NAME:
!
!   avgClusterProf
!
! PURPOSE:
!
!   Compute average over the constituents of a cluster.
!
! SYNTAX:
!
!   CALL avgClusterProf(clustStart, clustNmem, memberList, nSamples, 
!      IG, NG, molID, temp, h2o, tskin, psfc, pland, elev, prof, 
!      plandAvg, elevAvg, invertGrid, qc)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   clustStart  INTEGER             Starting pixel index within Cluster
!   clustNmem   INTEGER             Number of pixels within Cluster
!   memberList  INTEGER             List of cluster members
!   nSamples    INTEGER             number of samples
!   IG          TYPE(STATEINDEX_T)  Starting indices for sections of 
!                                   geophysical state vector 
!   NG          TYPE(STATEINDEX_T)  Number of elements for sections of 
!                                   geophysical state vector 
!   molID       INTEGER             List of IDs of relevant molecular species
!   temp*       REAL                Temperature
!   h2o*        REAL                water vapor amount
!   tskin*      REAL                Surface skin temperature
!   psfc*       REAL                Surface pressure
!   pland*      REAL                Fraction of land (versus water) of the 
!                                   field of view 
!   elev*       REAL                Elevation of surface
!   invertGrid* LOGICAL             Flag for reversing the order of the grid
!   
!   INPUTS/OUTPUTS:
!   
!   prof*       REAL                Atmospheric profile
!   plandAvg*   REAL                pland representing average
!   elevAvg*    REAL                elevation averaged
!   qc          INTEGER             quality control index
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

  integer, dimension(:),   intent(in)              :: memberList !0-based
  integer,                 intent(in)              :: nSamples
  type(StateIndex_t),      intent(in)              :: IG, NG 
  integer, dimension(:),   intent(in)              :: MolID
  real,    dimension(:,:), intent(in),    optional :: temp,h2o
  real,    dimension(:),   intent(in),    optional :: tskin,psfc,pland,elev
  logical,                 intent(in),    optional :: invertGrid  
  integer, dimension(:),   intent(inout)           :: qc
  real,    dimension(:),   intent(inout), optional :: prof
  real,                    intent(inout), optional :: plandAvg,elevAvg

  integer                                           :: i,iSamp
  integer                                           :: H2oID
  logical                                           :: invert

  if (present(invertGrid)) then
     invert=invertGrid
  else
     invert=.FALSE.
  endif

  qc(1)=0
  H2oID = whereH2O(MolID)
  if (present(prof)) prof=0
  if (present(plandAvg)) plandAvg=0
  if (present(elevAvg))  elevAvg=0
  
  do i=clustStart+1,clustStart+clustNmem
     iSamp=memberList(i)+1

     if (present(temp) .and. present(prof))  then
        if (any(temp(:,iSamp) == MISSING_REAL)) then
           qc(1)=ibset(qc(1),0)
           prof(IG%temp:IG%temp+NG%temp-1)=MISSING_REAL
        else
           prof(IG%temp:IG%temp+NG%temp-1)= &
                prof(IG%temp:IG%temp+NG%temp-1)+temp(:,iSamp)
        endif
     endif

     if (present(h2o) .and. present(prof)) then
        if (any(h2o(:,iSamp)== MISSING_REAL)) then
           qc(1)=ibset(qc(1),0)
           prof(IG%mol(H2oID):IG%mol(H2oID)+NG%mol(H2oID)-1)=MISSING_REAL
        else
           prof(IG%mol(H2oID):IG%mol(H2oID)+NG%mol(H2oID)-1)= &
                prof(IG%mol(H2oID):IG%mol(H2oID)+NG%mol(H2oID)-1)+ h2o(:,iSamp)
        endif
     endif

     if (present(tskin) .and. present(prof)) then
        if (tskin(iSamp) == MISSING_REAL) then
           qc(1)=ibset(qc(1),0)
           prof(IG%tskin)=MISSING_REAL
        else
           prof(IG%tskin)=prof(IG%tskin)+ tskin(iSamp)
        endif
     endif

     if (present(psfc) .and. present(prof)) then
        if (psfc(iSamp) == MISSING_REAL) then
           qc(1)=ibset(qc(1),0)
           prof(IG%psfc)=MISSING_REAL
        else
           prof(IG%psfc)=prof(IG%psfc)+ psfc(iSamp)
        endif
     endif

     if (present(pland) .and. present(plandAvg)) then
        if (pland(iSamp) == MISSING_REAL) then
           qc(1)=ibset(qc(1),0)
           plandAvg=MISSING_REAL
        else
           plandAvg=plandAvg+ pland(iSamp)
        endif
     endif

     if (present(elev) .and. present(elevAvg)) then
        if (elev(iSamp) == MISSING_REAL) then
           qc(1)=ibset(qc(1),0)
           elevAvg=MISSING_REAL
        else
           elevAvg=elevAvg+ elev(iSamp)
        endif
     endif
     
     if (btest(qc(1),0)) then
        if (present(prof)) prof=MISSING_REAL
        if (present(plandAvg)) plandAvg=MISSING_REAL
        if (present(elevAvg))  elevAvg=MISSING_REAL
        return
     endif
  enddo

  if (present(prof)) then
     prof=prof/clustNmem
     if (invert) then
        prof(IG%temp:IG%temp+NG%temp-1)=prof(IG%temp+NG%temp-1:IG%temp:-1)
        prof(IG%mol(H2oID):IG%mol(H2oID)+NG%mol(H2oID)-1)=prof(IG%mol(H2oID)+NG%mol(H2oID)-1:IG%mol(H2oID):-1)
     endif
     ! Convert to log(mixing ratio)
     prof(IG%mol(H2oID):IG%mol(H2oID)+NG%mol(H2oID)-1) = &
          alog(max(prof(IG%mol(H2oID):IG%mol(H2oID)+NG%mol(H2oID)-1),1.e-7))
  endif
  
  if (present(plandAvg)) plandAvg=plandAvg/clustNmem
  if (present(elevAvg))  elevAvg=elevAvg/clustNmem
  return
end subroutine avgClusterProf
!-------------------------------------------------------------------------------------------
subroutine checkCldliq(prof,IG,NG,qc)  

!<f90Subroutine>********************************************************
!
! NAME:
!
!   checkCldliq
!
! PURPOSE:
!
!   QC for cloud liquid top and thickness.
!
! SYNTAX:
!
!   CALL checkCldliq(prof, IG, NG, qc)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   IG    TYPE(STATEINDEX_T)  Starting indices for sections of geophysical 
!                             state vector 
!   NG    TYPE(STATEINDEX_T)  Number of elements for sections of geophysical 
!                             state vector 
!   
!   INPUTS/OUTPUTS:
!   
!   prof  REAL                Atmospheric profile
!   qc    INTEGER             quality control index
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

  type(StateIndex_t),    intent(in)    :: IG, NG 
  real,    dimension(:), intent(inout) :: prof
  integer, dimension(:), intent(inout) :: qc
  real,    parameter            :: Ptop_def = 500.0
  real,    parameter            :: Pthick_def = 50.0
  real,    parameter            :: Ptop_lower = 0.0
  real,    parameter            :: Ptop_upper = 1200.
  real,    parameter            :: Pthick_lower = 0.0
  real,    parameter            :: Pthick_upper = 1200.0
  
  if (prof(IG%cldliq+2) <= 0.0) then 
     prof(IG%cldliq) = Ptop_def
     prof(IG%cldliq+1) = Pthick_def
     prof(IG%cldliq+2) = 0.0
  else 
     if ( prof(IG%cldliq)   <= Ptop_lower   .or. &
          prof(IG%cldliq)   >= Ptop_upper   .or. &
          prof(IG%cldliq+1) <= Pthick_lower .or. &
          prof(IG%cldliq+1) >= Pthick_upper ) then
        print*, prof(IG%cldliq:IG%cldliq+NG%cldliq-1)
        print *,'err[InterpNWPtools::checkCldliq]:  invalid cloud top and thickness.'
        qc(1)=ibset(qc(1),0)
     endif
  endif
  return
end subroutine checkCldliq
!-------------------------------------------------------------------------------------------
function h2oMixr(T,pGrid,rh,nlev)

!<f90Function>**********************************************************
!
! NAME:
!
!   h2oMixr
!
! PURPOSE:
!
!   Manipulate relative humidity and loop over levels for conversion from 
!   relative humidity to mixing ratio.
!
! SYNTAX:
!
!   Results=h2oMixr(T, pGrid, rh, nlev)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   T        REAL     Temperature
!   pGrid    REAL     Pressure grid
!   rh       REAL     Relative humidity
!   nlev     INTEGER  Number of atmospheric levels
!
!   * OPTIONAL
!
! RETURN:
!
!     REAL     
!
! INCLUDES:
!
!   None
!
!*********************************************************</f90Function>

  real,    dimension(nlev)           :: h2oMixr ! q in g/g
  real,    dimension(:), intent (in) :: T,pGrid,rh  ! T in K, pGrid in mb, rh in %
  integer,               intent (in) :: nlev
  integer                            :: i
  real                               :: rhAdj
  real, parameter :: pctToFrac=0.01
  real, parameter :: minRHpct=1.
  do i=1,nlev
     rhAdj=pctToFrac*max(rh(i),minRHpct)
     h2oMixr(i)=mrFromRH(rhAdj,T(i),pGrid(i))
  enddo
  return
end function h2oMixr
!-------------------------------------------------------------------------------------------
subroutine readDataAtt(findex, varName, missing_value,valid_min,valid_max,scale,offset)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   readDataAtt
!
! PURPOSE:
!
!   Read attributes from the input ncdf file.
!
! SYNTAX:
!
!   CALL readDataAtt(findex, varName, missing_value, valid_min, 
!      valid_max, scale, offset)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   findex         INTEGER  file ID index
!   varName        CHAR     Name of variable
!   
!   INPUTS/OUTPUTS:
!   
!   missing_value  REAL     Value representing missing data
!   valid_min      REAL     Minimum valid value
!   valid_max      REAL     Maximum valid value
!   scale          REAL     Data compression scale factor
!   offset         REAL     Data compression offset
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

  integer,           intent(in)    :: findex
  character (len=*), intent(in)    :: varName
  real,              intent(inout) :: missing_value,valid_min,valid_max
  real,              intent(inout) :: scale,offset
  integer                        :: status,vid

  status = nf_inq_varid(findex,trim(varName),vid)
  if (status .ne. nf_noerr) print*,nf_strerror(status)
  status = nf_get_att_real(findex, vid, 'missing_value',missing_value)
  if (status .ne. nf_noerr) print*,nf_strerror(status)
  status = nf_get_att_real(findex, vid, 'valid_min',valid_min)
  if (status .ne. nf_noerr) print*,nf_strerror(status)
  status = nf_get_att_real(findex, vid, 'valid_max',valid_max)
  if (status .ne. nf_noerr) print*,nf_strerror(status)
  status = nf_get_att_real(findex, vid, 'scale',scale)
  if (status .ne. nf_noerr) print*,nf_strerror(status)
  status = nf_get_att_real(findex, vid, 'offset',offset)
  if (status .ne. nf_noerr) print*,nf_strerror(status)
  return
end subroutine readDataAtt
end module interpNWPtools
