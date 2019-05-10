MODULE EmisTab

  ! This module produces emissivity data from tables of gridded data.
  ! The user provides a netCDF file that contains the data.
  ! The data may represent the total emissivity, the wind-induced component
  ! of emissivity, etc.  There may be more than one type of data in the file. 
  ! Each type of data has a distinct netCDF variable name, and the 
  ! implemented options are identified in the "select case(dataType)" block.
  ! The data must be on a grid of 6 dimensions, where the dimensions 
  ! respectively cover: salinity(ppt), temperature(K), wind speed(m/s),
  ! relative wind direction(degrees), earth incidence angle(degrees),
  ! and channel.  With the exception of channel, the dimension may be 1,
  ! whereby the emissivity is assumed not to vary in that dimension.
  ! Channels are described by frequency-polarization pairs.  Polarization is
  ! identified by standard indices, defined in character cpolAll.
  ! The data may be provided for a particular channel set, or on a standard
  ! frequency grid.  If the data are on standard frequencies, there should be
  ! 2 (V,H) or 4 (V,H,U,F) Stokes parameters per frequency, where the 3rd 
  ! and 4th are assumed 0 when it contains 2 Stokes.  The data must be ordered
  ! to cycle through all Stokes parameters before incrementing frequency.
  ! If the user wants data for a certain channel set, they may provide the 
  ! channel definitions at initialization, whereby the table will be re-gridded
  ! to the given channels and stored, so no further channel re-gridding
  ! is needed when each emissivity is requested.
  ! The subroutine that provides emissivity data is overloaded with these
  ! options: 1) user supplies a list of channel indices and obtains a vector
  ! of emissivities, 2) user supplies a single (scalar) channel index and 
  ! obtains a scalar emissivity, 3) user supplies a vectors of frequency and
  ! polarization and obtains a vector of emissivities, and 4) user supplies 
  ! a frequency-polarization pair and obtains a scalar emissivity.  
  ! If the user is doing many emissivity computations for a static 
  ! channel set, it is more efficient not to use options 3 or 4, 
  ! but to supply the static channel definitions as the 
  ! frqOut and polOut arguments to emisTabInit and then use options 1 or 2
  ! to obtain the interpolated emissivities.
  ! The netCDF file must contain these anciliary data: 
  ! - the nodes of the grid for each dimension of the data
  ! - strings that define the type of interpolation (e.g. linear) to use
  !   for each dimension (it is assumed that the table was built to achieve
  !   required accuracy with a certain type of interpolation)
  ! - a flag for each grid that indicates whether the grid is uniformly
  !   spaced, which speeds computation
  ! - a flag that indicates whether the data are provided on a standard
  !   frequency grid, and are thus suitable for interpolation to frequencies
  !   of a channel set

  use netCDFhelper, only: openNcdfFile, closeNcdfFile, readNcdfAttr, &
     nf_noerr, nf_strerror, nf_inq_varname, nf_inq_varid, nf_inq_vardimid, &
     nf_inq_dimlen, nf_get_var_real
  implicit none
  private
  public:: emisTabInit,getEmisTab,emisTabDestroy
  public:: rotMapPol,mapPolDefault

  ! Global private data
  
  integer, parameter    :: mxiv=12  ! length of variable types
  integer, parameter    :: mxty=5   ! max # number of data types
  character(len=mxiv), dimension(mxty) :: dataType

  integer, parameter                    :: mxim=12  ! length of interp types
  integer, parameter                    :: mxit=2   ! number of interp types
  integer, parameter                    :: mxfd=5   ! number of fields
  character(len=mxim), dimension(mxfd)  :: tySal,tyT,tySp,tyDr,tyAn
  integer                               :: imSal,imT,imSp,imDr,imAn
  character(len=mxim), dimension(mxit)  :: validInterpTypes=(/'none  ','linear'/)
  integer,             dimension(mxit)  :: i1add=0
  integer,             dimension(mxit)  :: i2add=0
  
  integer,             dimension(mxfd)  :: unigSal,unigT,unigSp,unigDr,unigAn
  integer,             dimension(mxfd)  :: nGauss_in
  integer                               :: frqPolGridded
  real,                dimension(mxfd)  :: stepSal,stepT,stepSp,stepDr,stepAn
  integer,             dimension(mxfd)  :: nStokes
  logical,             dimension(mxfd)  :: initialized=.false.
  
  logical :: first_init=.true.
  logical :: outBndsExit=.true.
  
  integer,parameter :: nPols=8
  character(nPols)  :: cpolAll = 'VHPMRLUF'  ! vertical, horizontal, +45deg,
!                                 -45deg, right circular, left circular,
!                                 3rd Stokes, 4th Stokes
  real,    dimension(4,nPols) :: mapPol
  real,    dimension(4,nPols) :: mapPolDefault=reshape((/1. ,0. , 0. , 0. , &
                                                         0. ,1. , 0. , 0. , &
                                                         0.5,0.5, 0.5, 0. , &
                                                         0.5,0.5,-0.5, 0. , &
                                                         0.5,0.5, 0. , 0.5, &
                                                         0.5,0.5, 0. ,-0.5, &
                                                         0. ,0. , 1. , 0. , &
                                                         0., 0. , 0. , 1.   /)&
                                                         ,(/4,nPols/))
  
  real, dimension(:,:,:,:,:,:),    pointer     :: datArr,datArr2
  real, dimension(:,:,:,:,:),      pointer     :: tmp5,tmp5ds
  real, dimension(:,:,:,:),        pointer     :: tmp4,tmp4ds,tmp4dt
  real, dimension(:,:,:),          pointer     :: tmp3,tmp3ds,tmp3dt,tmp3dw
  real, dimension(:,:),            pointer     :: tmp2,tmp2ds,tmp2dt,tmp2dw,tmp2dd
  real, dimension(:),              pointer     :: grSal,grT,grSp,grDr,grAn
  
  type Variables
     real, dimension(:,:,:,:,:,:), pointer     :: datArr
     real, dimension(:,:,:,:,:),   pointer     :: tmp5,tmp5ds
     real, dimension(:,:,:,:),     pointer     :: tmp4,tmp4ds,tmp4dt
     real, dimension(:,:,:),       pointer     :: tmp3,tmp3ds,tmp3dt,tmp3dw
     real, dimension(:,:),         pointer     :: tmp2,tmp2ds,tmp2dt,tmp2dw,tmp2dd
     real, dimension(:),           pointer     :: grSal,grT,grSp,grDr,grAn
     real,    dimension(:),        pointer     :: emisOn,emisGrd
     integer, dimension(:),        pointer     :: mapToGr
     real,    dimension(:),        pointer     :: myfrqGrd
     integer, dimension(:),        pointer     :: mypolGrd

  end type Variables
  type(Variables), dimension(mxfd)             :: allVar
  
  interface getEmisTab
     module procedure getEmisTabChanVec
     module procedure getEmisTabChanSngl
     module procedure getEmisLiveChanVec
     module procedure getEmisLiveChanSngl
  end interface
  
CONTAINS
  
  SUBROUTINE emisTabInit(indx,fname,dataType_in, &
       frqOut,polOut,outBndsExit_in,polRotates,satAlt,frqGrd,polGrd, &
       nChan,errMax,nGaussOut,gaussAngOut)
    
    ! Initialize access to tabulated emissivity data
    !  If data table is on a standard frequency/polarization grid, the user
    !  can specify a channel set on which output will be provided 
    !  (frqOut, polOut), and can further specify whether polarization should
    !  rotate with nadir angle (as in cross-track scanners).  If rotated
    !  polarization is selected, the user must supply the satellite altitude
    !  (satAlt).
    !  If data table represents a particular channel set, the user may verify 
    !  that the channel set is as expected by returning frqGrd, polGrd, and
    !  nChan.  If it is expected that the tabulated data have polarization
    !  rotating with nadir angle, that property can be verified and the
    !  associated satellite altitude can be checked by using polRotates and
    !  satAlt as output arguments.

    use OSSPracticalConstant, only: MISSING_REAL

    integer,                      intent(in)              :: indx
    character(len=*),             intent(in)              :: fname
    character(len=*),             intent(in)              :: dataType_in
    real,    dimension(:),        intent(in),    optional :: frqOut
    integer, dimension(:),        intent(in),    optional :: polOut
    logical,                      intent(in),    optional :: outBndsExit_in
    logical,                      intent(inout), optional :: polRotates
    real,                         intent(inout), optional :: satAlt
    real,    dimension(:),        pointer,       optional :: frqGrd
    integer, dimension(:),        pointer,       optional :: polGrd
    integer,                      intent(inout), optional :: nChan
    real,                         intent(inout), optional :: errMax
    integer,                      intent(inout), optional :: nGaussOut
    real,  dimension(:), pointer,                optional :: gaussAngOut

    !-- Local variables
    
    integer, dimension(6) :: mydimid
    integer               :: tmpvarid
    real                  :: w1,w2
    integer               :: ncid,ncStatus,istatus,i,j1,j2,k,nChanOut
    integer               :: neSal,neT,neSp,neDr,neAn
    logical               :: myPolRotates  ! Polarization rotates with 
                                           !   cross-track nadir angle
    integer               :: iAn,tabPolRotates
    real                  :: angle,tabSatAlt
    integer               :: mynChan
    
    
    if (indx > mxfd) then
       print *,'err[EmisTab::emisTabInit]: ', &
            'Trying to access field beyond the rank of the global table'
       call errorHalt(1)
    endif
    
    call emisTabDestroy(indx)  ! Start fresh with memory allocation
    
    !-- Read file
    
    call openNcdfFile(ncid, fname)
    
    !--The "1" in the statement below indicates that 
    !--the field is the 1st (and only) variable in the file
    
    ncStatus = nf_inq_varname(ncid,1,dataType(indx))
    
    ncStatus = nf_inq_varid(ncid,trim(dataType(indx)),tmpvarid)
    if (trim(adjustl(dataType_in)) /=  trim(adjustl(dataType(indx)))) then
       print *,'err[EmisTab::emisTabInit]: File does not contain dataType '// &
            trim(dataType_in)//' ('//trim(dataType(indx))//') '
       call errorHalt(1)
    endif
    ! Dimensions depend ultimately on dataType; file could contain > 1 datatype
    ncStatus = nf_inq_vardimid(ncid,tmpvarid,mydimid)
    if (ncStatus /= nf_noerr) call handle_nc_err(ncStatus,'emisTabInit', &
      'Failed to get dimension IDs')
    ncStatus = nf_inq_dimlen(ncid,mydimid(1),neSal)
    if (ncStatus /= nf_noerr) call handle_nc_err(ncStatus,'emisTabInit', &
      'Failed to get dimension of Sal')
    ncStatus = nf_inq_dimlen(ncid,mydimid(2),neT)
    if (ncStatus /= nf_noerr) call handle_nc_err(ncStatus,'emisTabInit', &
      'Failed to get dimension of T')
    ncStatus = nf_inq_dimlen(ncid,mydimid(3),neSp)
    if (ncStatus /= nf_noerr) call handle_nc_err(ncStatus,'emisTabInit', &
      'Failed to get dimension of Sp')
    ncStatus = nf_inq_dimlen(ncid,mydimid(4),neDr)
    if (ncStatus /= nf_noerr) call handle_nc_err(ncStatus,'emisTabInit', &
      'Failed to get dimension of Dr')
    ncStatus = nf_inq_dimlen(ncid,mydimid(5),neAn)
    if (ncStatus /= nf_noerr) call handle_nc_err(ncStatus,'emisTabInit', &
      'Failed to get dimension of An')
    ncStatus = nf_inq_dimlen(ncid,mydimid(6),mynChan)
    if (ncStatus /= nf_noerr) call handle_nc_err(ncStatus,'emisTabInit', &
      'Failed to get nChan')
    if (present(nChan)) nChan=mynChan
    
    allocate(allVar(indx)%grSal(neSal),allVar(indx)%grT(neT), &
         allVar(indx)%grSp(neSp),allVar(indx)%grDr(neDr),allVar(indx)%grAn(neAn)) 
    allocate(allVar(indx)%datArr(neSal,neT,neSp,neDr,neAn,mynChan))
    datArr => allVar(indx)%datArr 
    
    allocate(allVar(indx)%tmp2(neAn,mynChan),allVar(indx)%tmp2ds(neAn,mynChan), &
         allVar(indx)%tmp2dt(neAn,mynChan), &
         allVar(indx)%tmp2dw(neAn,mynChan),allVar(indx)%tmp2dd(neAn,mynChan), &
         allVar(indx)%tmp3(neDr,neAn,mynChan),allVar(indx)%tmp3ds(neDr,neAn,mynChan), &
         allVar(indx)%tmp3dt(neDr,neAn,mynChan),allVar(indx)%tmp3dw(neDr,neAn,mynChan), &
         allVar(indx)%tmp4(neSp,neDr,neAn,mynChan),allVar(indx)%tmp4ds(neSp,neDr,neAn,mynChan), &
         allVar(indx)%tmp4dt(neSp,neDr,neAn,mynChan), & 
         allVar(indx)%tmp5(neT,neSp,neDr,neAn,mynChan), &
         allVar(indx)%tmp5ds(neT,neSp,neDr,neAn,mynChan), &
         stat=istatus)
    
    if (istatus /=0 ) then
       print *,'err[EmisTab::emisTabInit]: ', &
            'Not enough memory for tmp2,tmp3,tmp4,tmp5; istatus:',istatus
       call errorHalt(1)
    endif
    
    ncStatus = nf_get_var_real(ncid,tmpvarid,allVar(indx)%datArr)
    if (ncStatus /= nf_noerr) then
       print *,'err[EmisTab::emisTabInit]: (get_var_real:'//trim(dataType(indx))// &
            ') '//trim(nf_strerror(ncStatus))
       call errorHalt(1)
    endif
    
    call readNcdfAttr(ncid,attrName='nGauss',attr=nGauss_in(indx),silent=.true.)
    if(nGauss_in(indx) == 0)nGauss_in(indx)=1 ! --When nGauss_in(indx) is not present
    if (present(frqOut) .and. present(polOut) .and. nGauss_in(indx) > 1)then
       print *,'err[EmisTab::emisTabInit]: ', &
            'Cannot have nGauss > 1 when frqOut and polOut present;',nGauss_in(indx)
       call errorHalt(1)
    endif

    if(present(nGaussOut) .and. present(gaussAngOut))then
       nGaussOut=nGauss_in(indx)
       allocate(gaussAngOut(nGaussOut))
       call readNcdfAttr(ncid,attrName='GaussAngles',attr=gaussAngOut)
    endif
    
    call readNcdfAttr(ncid,attrName='gridSalin',attr=allVar(indx)%grSal)
    call readNcdfAttr(ncid,attrName='gridTemp',attr=allVar(indx)%grT)
    if(trim(dataType(indx)) == 'windSpeed')then
       call readNcdfAttr(ncid,attrName='gridEmiss',attr=allVar(indx)%grSp)
    else
       call readNcdfAttr(ncid,attrName='gridSpeed',attr=allVar(indx)%grSp)
    endif
    call readNcdfAttr(ncid,attrName='gridDirec',attr=allVar(indx)%grDr)
    call readNcdfAttr(ncid,attrName='gridAngle',attr=allVar(indx)%grAn)
    
    allocate(allVar(indx)%myfrqGrd(mynChan/nGauss_in(indx)), &
             allVar(indx)%mypolGrd(mynChan/nGauss_in(indx)))
    
    call readNcdfAttr(ncid,attrName='frequencies',attr=allVar(indx)%myfrqGrd)
    call readNcdfAttr(ncid,attrName='polarizations',attr=allVar(indx)%mypolGrd)
    
    call readNcdfAttr(ncid,attrName='interpSalin',attr=tySal(indx))
    call readNcdfAttr(ncid,attrName='interpTemp',attr=tyT(indx))
    if(trim(dataType(indx)) == 'windSpeed')then
       call readNcdfAttr(ncid,attrName='interpEmiss',attr=tySp(indx))
    else
       call readNcdfAttr(ncid,attrName='interpSpeed',attr=tySp(indx))
    endif
    call readNcdfAttr(ncid,attrName='interpDirec',attr=tyDr(indx))
    call readNcdfAttr(ncid,attrName='interpAngle',attr=tyAn(indx))
    call interpTypeCheck('Salin',neSal,tySal(indx),imSal)
    call interpTypeCheck('Temp', neT,  tyT(indx),  imT  )
    if(trim(dataType(indx)) == 'windSpeed')then
       call interpTypeCheck('Emiss',neSp, tySp(indx), imSp )
    else
       call interpTypeCheck('Speed',neSp, tySp(indx), imSp )
    endif
    call interpTypeCheck('Direc',neDr, tyDr(indx), imDr )
    call interpTypeCheck('Angle',neAn, tyAn(indx), imAn )
    
    call readNcdfAttr(ncid,attrName='unifGridSalin',attr=unigSal(indx))
    call readNcdfAttr(ncid,attrName='unifGridTemp',attr=unigT(indx))
    if(trim(dataType(indx)) == 'windSpeed')then
       call readNcdfAttr(ncid,attrName='unifGridEmiss',attr=unigSp(indx))
    else
       call readNcdfAttr(ncid,attrName='unifGridSpeed',attr=unigSp(indx))
    endif
    call readNcdfAttr(ncid,attrName='unifGridDirec',attr=unigDr(indx))
    call readNcdfAttr(ncid,attrName='unifGridAngle',attr=unigAn(indx))
    stepSal(indx)=stepUniform('Salin',unigSal(indx),allVar(indx)%grSal)
    stepT(indx)  =stepUniform('Temp', unigT(indx),allVar(indx)%grT)
    if(trim(dataType(indx)) == 'windSpeed')then
       stepSp(indx) =stepUniform('Emiss',unigSp(indx),allVar(indx)%grSp)
    else
       stepSp(indx) =stepUniform('Speed',unigSp(indx),allVar(indx)%grSp)
    endif
    stepDr(indx) =stepUniform('Direc',unigDr(indx),allVar(indx)%grDr)
    stepAn(indx) =stepUniform('Angle',unigAn(indx),allVar(indx)%grAn)
    
    call readNcdfAttr(ncid,attrName='frqPolGridded',attr=frqPolGridded)
    call readNcdfAttr(ncid,attrName='polRotatesWithAngle',attr=tabPolRotates)
    if (tabPolRotates == 1) then
      call readNcdfAttr(ncid,attrName='satelliteAltitude',attr=tabSatAlt)
    else
      tabSatAlt=MISSING_REAL
    endif
    
    if (present(frqGrd)) then
      if (associated(frqGrd)) deallocate(frqGrd)
      allocate (frqGrd(mynChan))
      frqGrd=allVar(indx)%myfrqGrd
    endif
    if (present(polGrd)) then
      if (associated(polGrd)) deallocate(polGrd)
      allocate (polGrd(mynChan))
      polGrd=allVar(indx)%mypolGrd
    endif
    if (present(errMax)) then
      call readNcdfAttr(ncid,attrName='errorLimit',attr=errMax)
    endif

    call closeNcdfFile(ncid)

    ! Handle pre-process for data that come on a frequency grid

    if (frqPolGridded /= 0) then
      ! Check input grid polarization sequence
      call checkPolGrid(indx,allVar(indx)%mypolGrd,nStokes(indx))
      ! Weights for mapping from tabulated data to each polarization
      mapPol=mapPolDefault

    ! Convert from gridded frequencies to channel frequencies
      if (present(frqOut) .and. present(polOut)) then
        nChanOut=size(frqOut)
        if (size(polOut) /= nChanOut) then
          print *,'err[EmisTab::emisTabInit]: ', &
          'Mismatch sizes of frqOut and polOut;',size(frqOut),size(polOut)
          call errorHalt(1)
        endif
        if (tabPolRotates /= 0) then
          print *,'err[EmisTab::emisTabInit]: ', &
          'Does not support regridding from data tabulated in rotated form'
          call errorHalt(1)
        endif
        if (present(polRotates)) then
          myPolRotates=polRotates
          if (myPolRotates .and. (.not. present(satAlt))) then
            print *,'err[EmisTab::emisTabInit]: ', &
            'Cannot rotate polarization without satellite altitude'
            call errorHalt(1)
          endif
        else
          myPolRotates=.false.
        endif
        allocate(datArr2(neSal,neT,neSp,neDr,neAn,nChanOut),stat=istatus)
        if (istatus /=0 ) then
          print *,'err[EmisTab::emisTabInit]: ', &
          'Not enough memory for datArr2; istatus:',istatus
          call errorHalt(1)
        endif

        ! Channel position in frequency grid
        if (minval(frqOut) < allVar(indx)%myfrqGrd(1) .or. &
            maxval(frqOut) > allVar(indx)%myfrqGrd(mynChan)) then
          print *,'err[EmisTab::emisTabInit]: ', &
            'Requested freq exceeds range of tabulated data'
          print *,'Requested:',minval(frqOut),maxval(frqOut), &
            '  Tabulated:', &
            allVar(indx)%myfrqGrd(1),allVar(indx)%myfrqGrd(mynChan)
          call errorHalt(1)
        endif
        datArr2=0.
        do i=1,nChanOut
          call getFrqGridDat(frqOut(i),allVar(indx)%myfrqGrd,nStokes(indx), &
            j1,j2,w1,w2)
        ! Execute frequency interpolation and polarization mapping
          do iAn=1,neAn
            angle=allVar(indx)%grAn(iAn)
            if (myPolRotates) mapPol=rotMapPol(angle,satAlt)
            do k=1,nStokes(indx)
              datArr2(:,:,:,:,iAn,i)=datArr2(:,:,:,:,iAn,i)+ &
                datArr(:,:,:,:,iAn,j1+k-1)*mapPol(k,polOut(i)+1)*w1+ &
                datArr(:,:,:,:,iAn,j2+k-1)*mapPol(k,polOut(i)+1)*w2
            enddo
          enddo
        enddo
        mynChan=nChanOut
        deallocate(allVar(indx)%datArr)
        allVar(indx)%datArr=>datArr2
        nullify(datArr2)
        frqPolGridded=0
      else  ! Returning data on tabulated frequency grid
        if (present(polRotates) .or. present(satAlt)) then
          print *,'warning[EmisTab::emisTabInit]: '// &
          'Optional arguments polRotates and satAlt '
          print *,'should not be present when table is on standard grid '// &
          'and no channel set '
          print *,'conversion is specified by frqOut and polOut'
        endif
      endif
    else   ! Tabulated data were not on a frequency grid
      if (present(frqOut) .or. present(polOut)) then
        if (any(frqOut /= allVar(indx)%myfrqGrd) .or. &
          any(polOut /= allVar(indx)%mypolGrd)) then
          print *,'err[EmisTab::emisTabInit]: ', &
          'Requested channel regrid but source is not gridded'
          call errorHalt(1)
        endif
      endif
      if (present(polRotates)) then
        if (tabPolRotates == 1) then
          polRotates=.true.
        else
          polRotates=.false.
        endif
      endif
      if (present(satAlt)) then
        satAlt=tabSatAlt
      endif
    endif

    initialized(indx)=.true.

    if(present(outBndsExit_in))outBndsExit=outBndsExit_in
    
    return
    
  END SUBROUTINE emisTabInit
  
  !--------------------------------------------------------------------------

  SUBROUTINE handle_nc_err(ncStatus,subproName,msg)
    integer,           intent(in)           :: ncStatus
    character (len=*), intent(in)           :: subproName
    character (len=*), intent(in), optional :: msg
    character (len=256) :: lmsg

    if (present(msg)) then
       if (len(msg) > 256) print*,'message is truncated'
       lmsg=trim(msg)
    else
       lmsg=':::'
    endif
    print *, 'err[ncdf_module::handle_nc_err(',trim(subproName),')]: ',&
      ' ('//trim(msg)//') '//trim(nf_strerror(ncStatus))
    call errorHalt(1)

  END SUBROUTINE handle_nc_err

  !--------------------------------------------------------------------------

  SUBROUTINE interpTypeCheck(name,nelem,type,itype)

    character(len=*), intent(in)    :: name
    integer,          intent(in)    :: nelem
    character(len=*), intent(inout) :: type
    integer,          intent(inout) :: itype

    character(len=mxim) :: tmptype
    integer             :: i

    itype=0
    if (nelem == 1) then
      type='none'      ! Override interpolation type when only one node
      print *,'warning[EmisTab::interpTypeCheck]: '
      print *,name,' grid has only 1 node; will be assumed constant'
    endif

    ! Get integer index of interpolation type
    do i=1,size(validInterpTypes)
      tmptype=validInterpTypes(i)
      if (trim(adjustl(type)) == trim(adjustl(tmptype))) then
        itype=i
        exit
      endif
    enddo
    if (itype == 0) then
      print *,'err[EmisTab::interpTypeCheck]: '
      print *,'For '//trim(name)//', interp type "'//trim(adjustl(type))// &
        '" not in valid types:'
      print *,validInterpTypes
      call errorHalt(1)
    endif
    
    return

  END SUBROUTINE interpTypeCheck

  !--------------------------------------------------------------------------

  FUNCTION stepUniform(name,unig,grid)

    character(len=*),   intent(in)  :: name
    integer,            intent(in)  :: unig
    real, dimension(:), intent(in)  :: grid
    real                            :: stepUniform
                        
    real, parameter :: deplim=0.001
    integer         :: i,ngrid
    real            :: step,dep

    step=0.
    if (unig /= 1) return
    ngrid=size(grid)
    if (ngrid == 1) return

    step=(grid(ngrid)-grid(1))/float(ngrid-1)

    do i=2,size(grid)
      dep=abs(grid(i)-grid(1)-step*float(i-1))/step
      if (dep > deplim) then
        print *,'err[EmisTab::stepUniform]: ', &
        name,' grid claimed to be uniform but was not; i,dep:',i,dep
        print *,grid
        call errorHalt(1)
      endif
    enddo

    stepUniform=step
    return

  END FUNCTION stepUniform

  !--------------------------------------------------------------------------

  SUBROUTINE checkPolGrid(indx,mypolGrd,nStokes)

  ! Check the sequence of polarizations in the frequency/polarization
  ! grid of the input tabulated data, to ensure it is valid

    integer,               intent(in)    :: indx
    integer, dimension(:), intent(in)    :: mypolGrd
    integer,               intent(inout) :: nStokes

    integer :: i,nChanGr

    ! Find repeat cycle of polarizations
    nChanGr=size(mypolGrd)
    do i=2,nChanGr
      if (mypolGrd(i) == mypolGrd(1)) exit
    enddo
    nStokes=i-1

    ! Verify that repeat cycle is consistent for entire grid
    do i=nStokes+1,nChanGr,nStokes
      if (any(myPolGrd(i:i-1+nStokes) /= myPolGrd(1:nStokes))) then
        print *,'err[EmisTab::checkPolGrid]: ', &
        'Polarization grid of input is not consistent; '
        print *,'First set:',myPolGrd(1:nStokes)
        print *,i,i-1+nStokes,'range set:',myPolGrd(i:i-1+nStokes)
        call errorHalt(1)
      endif
    enddo

    ! Verify that the polarization sequence matches expected for mapPol
    if (nStokes == 2) then
      if (any(myPolGrd(1:nStokes) /= &
        (/index(cpolAll,'V')-1,index(cpolAll,'H')-1/))) then
        print *,'err[EmisTab::checkPolGrid]: ', &
        'Database polarization indexes not as expected; ',myPolGrd(1:nStokes)
        print *,'Should be for V,H'
        call errorHalt(1)
      endif
      print *,'Database contains only V and H; 3rd and 4th Stokes assumed zero'
    elseif (nStokes == 4) then
      if (any(myPolGrd(1:nStokes) /= &
        (/index(cpolAll,'V')-1,index(cpolAll,'H')-1, &
        index(cpolAll,'U')-1,index(cpolAll,'F')-1/))) then
        print *,'err[EmisTab::checkPolGrid]: ', &
        'Database polarization indices not as expected; ',myPolGrd(1:nStokes)
        print *,'Should be for V, H, and 3rd and 4th Stokes (U,F)'
        call errorHalt(1)
      endif
    else
      print *,'err[EmisTab::checkPolGrid]: '
      print *,'Expected database to have V,H or V,H,U,F polarizations;', &
        myPolGrd(1:nStokes)
      call errorHalt(1)
    endif

    allocate (allVar(indx)%mapToGr(nChanGr), &
              allVar(indx)%emisOn(nChanGr), &
              allVar(indx)%emisGrd(nChanGr))

    return

  END SUBROUTINE checkPolGrid

  !--------------------------------------------------------------------------

  FUNCTION rotMapPol(angle,satAlt)

  ! Produce map (weights) for converting from V H polarization to the
  ! quasi-V and quasi-H that occurs for some cross-track scanners, where
  ! polarization rotates with nadir angle.
  ! The map does not neglect the 3rd Stokes (U), but it does not produce
  ! the map elements for any rotated channels other than V and H.

  use OSSPhysicalConstant, ONLY: deg2rad
  use AngleUtil, ONLY: nadirFromZenith

  real,                        intent(in) :: angle
  real,                        intent(in) :: satAlt
  real,    dimension(4,nPols)             :: rotMapPol

  real :: nadir,theta,cth,sth,c2th,s2th,sc

  ! Convert from Earth incidence angle to nadir angle
  nadir=nadirFromZenith(angle,satAlt)
  theta=nadir*deg2rad

  ! Map the rotation
  cth=cos(theta)
  sth=sin(theta)
  c2th=cth**2
  s2th=sth**2
  sc=sth*cth
  
  rotMapPol=0.
  rotMapPol(:,1)=(/c2th,s2th, 2.*sc,0./)   ! Quasi-V from V,H,U,F
  rotMapPol(:,2)=(/s2th,c2th,-2.*sc,0./)   ! Quasi-H from V,H,U,F

  return

  END FUNCTION rotMapPol

  !--------------------------------------------------------------------------

  SUBROUTINE getFrqGridDat(frq,myfrqGrd,nStokes,j1,j2,w1,w2)

    real,                  intent(in)    :: frq
    real,    dimension(:), intent(in)    :: myfrqGrd
    integer,               intent(in)    :: nStokes
    integer,               intent(inout) :: j1
    integer,               intent(inout) :: j2
    real,                  intent(inout) :: w1
    real,                  intent(inout) :: w2

    integer :: nChanGr

    ! Find indices of bounding frequencies
    ! Assume already have done check for out of bounds frq
    nChanGr=size(myfrqGrd)
    do j2=1+nStokes,nChanGr,nStokes
      if (myfrqGrd(j2) >= frq) exit
    enddo
    j1=j2-nStokes
  
 ! Weights for frequency interpolation
    w2=(frq-myfrqGrd(j1))/(myfrqGrd(j2)-myfrqGrd(j1))
    w1=1.-w2

    return

  END SUBROUTINE getFrqGridDat

  !--------------------------------------------------------------------------

  SUBROUTINE getEmisTabChanVec(indx,intrpScheme,listChans,salin,temp,speed,direc,angle, &
       emis,demdsal,demdtemp,demdspeed,demddir)

  ! Base subroutine to get data for a set of particular channels in
  ! the tabulated channel grid. The input variable intrpScheme determines
  ! the interpolation scheme used: 
  ! intrpScheme = 0: The most general form. It can be slower than 1 and 2 if
  !                  several dimensions are degenerate. Used by the table generator.
  ! intrpScheme = 1: only T dependency is included (e.g., specular)
  ! intrpScheme = 2: only wind dependency is included 
  !                  (e.g., wind-induced emissivity and reflectivity in the AER model)
  ! Options 1 and 2 are used by the surface module.

    integer,               intent(in)              :: indx
    integer,               intent(in)              :: intrpScheme
    real,                  intent(in)              :: angle
    integer, dimension(:), intent(in)              :: listChans
    real,                  intent(in),    optional :: salin
    real,                  intent(in),    optional :: temp
    real,                  intent(in),    optional :: speed
    real,                  intent(in),    optional :: direc
    real,    dimension(:), intent(inout)           :: emis
    real,    dimension(:), intent(inout), optional :: demdsal
    real,    dimension(:), intent(inout), optional :: demdtemp
    real,    dimension(:), intent(inout), optional :: demdspeed
    real,    dimension(:), intent(inout), optional :: demddir
    
    integer :: iSal1,iSal2,iT1,iT2,iSp1,iSp2,iDr1,iDr2,iAn1,iAn2,iEm1,iEm2
    integer :: iSal,iT,iSp,iDr,iAn,ich,lch,nChanGr
    
    if (indx > mxfd) then
       print *,'err[EmisTab::getEmisTabChanVec]: ', &
            'Trying to access field beyond the rank of the global table'
       call errorHalt(1)
    endif
    
    if (.not. initialized(indx)) then
       print *,'err[EmisTab::getEmisTabChanVec]: ', &
            'Must call emisTabInit before first call to getEmisTab'
       call errorHalt(1)
    endif
    
    ! Check validity of inputs
    nChanGr=size(allVar(indx)%datArr,dim=6)  ! 6th dimension is # channels
    if (any(listChans < 1) .or. any(listChans > nChanGr)) then
       print *,'err[EmisTab::getEmisTabChanVec]: ', &
            'Invalid channel number in list; nChan=',nChanGr
       print *,'Requested channels:',listChans
       call errorHalt(1)
    endif
    if (size(listChans) /= size(emis)) then
       print *,'err[EmisTab::getEmisTabChanVec]: ', &
            'Mismatch sizes of listChans and emis;',size(listChans),size(emis)
       call errorHalt(1)
    endif

    datArr => allVar(indx)%datArr 
    tmp5   => allVar(indx)%tmp5
    tmp5ds => allVar(indx)%tmp5ds
    tmp4   => allVar(indx)%tmp4
    tmp4ds => allVar(indx)%tmp4ds
    tmp4dt => allVar(indx)%tmp4dt
    tmp3   => allVar(indx)%tmp3
    tmp3ds => allVar(indx)%tmp3ds
    tmp3dt => allVar(indx)%tmp3dt
    tmp2   => allVar(indx)%tmp2
    tmp2ds => allVar(indx)%tmp2ds
    tmp2dt => allVar(indx)%tmp2dt
    tmp2dd => allVar(indx)%tmp2dd
    grSal  => allVar(indx)%grSal
    grT    => allVar(indx)%grT
    grSp   => allVar(indx)%grSp
    grDr   => allVar(indx)%grDr
    grAn   => allVar(indx)%grAn
    
    ! Find bounding grid indices
    call boundGrid('Salin',grSal,unigSal(indx),stepSal(indx),salin,iSal1,iSal2)
    call boundGrid('Temp', grT,  unigT(indx),  stepT(indx)  ,temp, iT1,  iT2)
    call boundGrid('Direc',grDr, unigDr(indx), stepDr(indx) ,direc,iDr1, iDr2)
    call boundGrid('Angle',grAn, unigAn(indx), stepAn(indx) ,angle,iAn1, iAn2)
       
    if(intrpScheme == 0)then 

       tmp3dw => allVar(indx)%tmp3dw
       tmp2dw => allVar(indx)%tmp2dw
       
       ! Find bounding grid indices for speed
       call boundGrid('Speed',grSp, unigSp(indx),  stepSp(indx) ,speed,iSp1, iSp2)
       
       do lch=1,size(listChans)
          ich=listChans(lch)
          do iAn=iAn1-i1add(imAn),iAn2+i2add(imAn)
             do iDr=iDr1-i1add(imDr),iDr2+i2add(imDr)
                do iSp=iSp1-i1add(imSp),iSp2+i2add(imSp)
                   do iT=iT1-i1add(imT),iT2+i2add(imT)
                      tmp5(iT,iSp,iDr,iAn,ich)= interpolate(tySal(indx),grSal,iSal1,iSal2, &
                           i1add(imSal),i2add(imSal),salin,datArr(:,iT,iSp,iDr,iAn,ich))
                      tmp5ds(iT,iSp,iDr,iAn,ich)= differentiate(tySal(indx),grSal,iSal1,iSal2, &
                           i1add(imSal),i2add(imSal),salin,datArr(:,iT,iSp,iDr,iAn,ich))
                   enddo
                   tmp4(iSp,iDr,iAn,ich)= interpolate(tyT(indx),grT,iT1,iT2, &
                        i1add(imT),i2add(imT),temp,tmp5(:,iSp,iDr,iAn,ich))
                   tmp4ds(iSp,iDr,iAn,ich)= interpolate(tyT(indx),grT,iT1,iT2, &
                        i1add(imT),i2add(imT),temp,tmp5ds(:,iSp,iDr,iAn,ich))
                   tmp4dt(iSp,iDr,iAn,ich)= differentiate(tyT(indx),grT,iT1,iT2, &
                        i1add(imT),i2add(imT),temp,tmp5(:,iSp,iDr,iAn,ich))
                enddo
                tmp3(iDr,iAn,ich)= interpolate(tySp(indx),grSp,iSp1,iSp2, &
                     i1add(imSp),i2add(imSp),speed,tmp4(:,iDr,iAn,ich))
                tmp3ds(iDr,iAn,ich)= interpolate(tySp(indx),grSp,iSp1,iSp2, &
                     i1add(imSp),i2add(imSp),speed,tmp4ds(:,iDr,iAn,ich))
                tmp3dt(iDr,iAn,ich)= interpolate(tySp(indx),grSp,iSp1,iSp2, &
                     i1add(imSp),i2add(imSp),speed,tmp4dt(:,iDr,iAn,ich))
                tmp3dw(iDr,iAn,ich)= differentiate(tySp(indx),grSp,iSp1,iSp2, &
                     i1add(imSp),i2add(imSp),speed,tmp4(:,iDr,iAn,ich))
             enddo
             tmp2(iAn,ich)=interpolate(tyDr(indx),grDr,iDr1,iDr2, &
                  i1add(imDr),i2add(imDr),direc,tmp3(:,iAn,ich))
             tmp2ds(iAn,ich)=interpolate(tyDr(indx),grDr,iDr1,iDr2, &
                  i1add(imDr),i2add(imDr),direc,tmp3ds(:,iAn,ich))
             tmp2dt(iAn,ich)=interpolate(tyDr(indx),grDr,iDr1,iDr2, &
                  i1add(imDr),i2add(imDr),direc,tmp3dt(:,iAn,ich))
             tmp2dw(iAn,ich)=interpolate(tyDr(indx),grDr,iDr1,iDr2, &
                  i1add(imDr),i2add(imDr),direc,tmp3dw(:,iAn,ich))
             tmp2dd(iAn,ich)=differentiate(tyDr(indx),grDr,iDr1,iDr2, &
                  i1add(imDr),i2add(imDr),direc,tmp3(:,iAn,ich))
          enddo
          emis(lch)=interpolate(tyAn(indx),grAn,iAn1,iAn2, &
               i1add(imAn),i2add(imAn),angle,tmp2(:,ich))
          if(present(demdsal))then
             demdsal(lch)=interpolate(tyAn(indx),grAn,iAn1,iAn2, &
                  i1add(imAn),i2add(imAn),angle,tmp2ds(:,ich))
          elseif(present(demdtemp))then
             demdtemp(lch)=interpolate(tyAn(indx),grAn,iAn1,iAn2, &
                  i1add(imAn),i2add(imAn),angle,tmp2dt(:,ich))
          elseif(present(demdspeed))then
             demdspeed(lch)=interpolate(tyAn(indx),grAn,iAn1,iAn2, &
                  i1add(imAn),i2add(imAn),angle,tmp2dw(:,ich))
          elseif(present(demddir))then
             demddir(lch)=interpolate(tyAn(indx),grAn,iAn1,iAn2, &
                  i1add(imAn),i2add(imAn),angle,tmp2dd(:,ich))
          endif
       enddo
       
       nullify(tmp3dw,tmp2dw)
       
    elseif(intrpScheme == 1)then 
       
       do lch=1,size(listChans)
          ich=listChans(lch)
          do iAn=iAn1-i1add(imAn),iAn2+i2add(imAn)
             tmp2(iAn,ich)=interpolate(tyT(indx),grT,iT1,iT2, &
                  i1add(imT),i2add(imT),temp,datArr(1,:,1,1,iAn,ich))
             tmp2dt(iAn,ich)=differentiate(tyT(indx),grT,iT1,iT2, &
                  i1add(imT),i2add(imT),temp,datArr(1,:,1,1,iAn,ich))
          enddo
          emis(lch)=interpolate(tyAn(indx),grAn,iAn1,iAn2, &
               i1add(imAn),i2add(imAn),angle,tmp2(:,ich))
          if(present(demdtemp))then
             demdtemp(lch)=interpolate(tyAn(indx),grAn,iAn1,iAn2, &
                  i1add(imAn),i2add(imAn),angle,tmp2dt(:,ich))
          endif
       enddo
       
    elseif(intrpScheme == 2)then

       tmp3dw => allVar(indx)%tmp3dw
       tmp2dw => allVar(indx)%tmp2dw
       
       ! Find bounding grid indices for emiss/speed
       if(trim(dataType(indx)) == 'windSpeed')then  ! This is for emis retrieval
          call boundGrid('Emiss',grSp, unigSp(indx),  stepSp(indx) ,speed,iSp1, iSp2)
       else                                         ! This is for wspd retrieval
          call boundGrid('Speed',grSp, unigSp(indx),  stepSp(indx) ,speed,iSp1, iSp2)
       endif
       
       do lch=1,size(listChans)
          ich=listChans(lch)
          do iAn=iAn1-i1add(imAn),iAn2+i2add(imAn)
             tmp2(iAn,ich)=interpolate(tySp(indx),grSp,iSp1,iSp2, &
                  i1add(imSp),i2add(imSp),speed,datArr(1,1,:,1,iAn,ich))
             tmp2dw(iAn,ich)=differentiate(tySp(indx),grSp,iSp1,iSp2, &
                  i1add(imSp),i2add(imSp),speed,datArr(1,1,:,1,iAn,ich))
          enddo
          emis(lch)=interpolate(tyAn(indx),grAn,iAn1,iAn2, &
               i1add(imAn),i2add(imAn),angle,tmp2(:,ich))
          if(present(demdspeed))then
             demdspeed(lch)=interpolate(tyAn(indx),grAn,iAn1,iAn2, &
                  i1add(imAn),i2add(imAn),angle,tmp2dw(:,ich))
          endif
       enddo

       nullify(tmp2dw)
       
    endif
    
    nullify(datArr)
    nullify(tmp5,tmp5ds,tmp4,tmp4ds,tmp4dt,tmp3,tmp3ds,tmp3dt,tmp2,tmp2ds,tmp2dt,tmp2dd)
    nullify(grSal,grT,grSp,grDr,grAn)
    
    return

  END SUBROUTINE getEmisTabChanVec

  !--------------------------------------------------------------------------
  
  SUBROUTINE getEmisTabChanSngl(indx,intrpScheme,iChan,salin,temp,speed,direc,angle, &
       emis,demdsal,demdtemp,demdspeed,demddir)
    
    ! Subroutine to get data for a particular single channel in
    ! the tabulated channel grid

    integer, intent(in)              :: indx
    integer, intent(in)              :: intrpScheme
    integer, intent(in)              :: iChan
    real,    intent(in)              :: angle
    real,    intent(in),    optional :: salin
    real,    intent(in),    optional :: temp
    real,    intent(in),    optional :: speed
    real,    intent(in),    optional :: direc
    real,    intent(inout)           :: emis
    real,    intent(inout), optional :: demdsal
    real,    intent(inout), optional :: demdtemp
    real,    intent(inout), optional :: demdspeed
    real,    intent(inout), optional :: demddir

    integer, dimension(1)            :: listChans
    real,    dimension(1)            :: myemis,mydemdsal,mydemdtemp, &
                                        mydemdspeed,mydemddir

    listChans=iChan
    if(present(speed))then
       call getEmisTabChanVec(indx,intrpScheme,listChans,salin,temp, &
            speed=speed,direc=direc,angle=angle,emis=myemis,demdsal=mydemdsal, &
            demdtemp=mydemdtemp,demdspeed=mydemdspeed,demddir=mydemddir)
       emis=myemis(1)
       if(present(demdsal))demdsal=mydemdsal(1)
       if(present(demdtemp))demdsal=mydemdtemp(1)
       if(present(demdspeed))demdspeed=mydemdspeed(1)
       if(present(demddir))demdsal=mydemddir(1)
    endif
    
    return    

  END SUBROUTINE getEmisTabChanSngl

  !--------------------------------------------------------------------------
  
  SUBROUTINE getEmisLiveChanVec(indx,intrpScheme,frq,pol,salin,temp,speed, &
       direc,angle,emis,polRotates,satAlt)

  ! Subroutine to get data for a set of specified frequencies and polarizations,
  ! where channelization is computed live (on the fly) from tabulated data.
  ! Optional argument satAlt should be present if and only if polRotates is 
  ! present and is true 

    integer,               intent(in)              :: indx
    integer,               intent(in)              :: intrpScheme
    real,    dimension(:), intent(in)              :: frq
    integer, dimension(:), intent(in)              :: pol
    real,                  intent(in)              :: salin
    real,                  intent(in)              :: temp
    real,                  intent(in)              :: speed
    real,                  intent(in)              :: direc
    real,                  intent(in)              :: angle
    real,    dimension(:), intent(inout)           :: emis
    logical,               intent(in),   optional  :: polRotates
    real,                  intent(in),   optional  :: satAlt

    integer :: i,k,l,nChanLoc,nOn,nChanGr
    integer, parameter :: mxcho=200  ! Maximum number of channels to output
    integer, dimension(mxcho) :: iflg,j1,j2
    real,    dimension(mxcho) :: w1,w2
    real                      :: absAngle

    if (nGauss_in(indx) > 1) then
      print *,'err[EmisTab::getEmisLiveChanVec]: ',&
        'nGauss > 1 not supported in live mode;',nGauss_in(indx)
      call errorHalt(1)
    endif

    if (frqPolGridded == 0) then
      print *,'err[EmisTab::getEmisLiveChanVec]: ',&
        'Source data are not gridded, so cannot '
      print *,'use option to select by frequency and polarization index'
      call errorHalt(1)
    endif

    nChanLoc=size(frq)
    if (nChanLoc /= size(pol) .or. nChanLoc /= size(emis)) then
      print *,'err[EmisTab::getEmisLiveChanVec]: ',&
        'Mismatch sizes of frq/pol/emis;',size(frq),size(pol),size(emis)
      call errorHalt(1)
    endif
    if (nChanLoc > mxcho) then  ! Dynamic memory could reallocate many times
      print *,'err[EmisTab::getEmisLiveChanVec]: ',&
        'Exceeded max number of channels to out data; nChanLoc,mxcho', &
        nChanLoc,mxcho
      call errorHalt(1)
    endif

  ! Allow user to input negative angle, interpreted as positive EIA and
  ! negative nadir angle in cross-track coordinates, to support rotating
  ! polarization with non-zero 3rd Stokes
    absAngle=abs(angle)

  ! Weights for mapping from tabulated data to each polarization
    mapPol=mapPolDefault
    if (present(polRotates)) then
      if (polRotates) then
        if (present(satAlt)) then
          mapPol=rotMapPol(angle,satAlt)
        else
          print *,'err[EmisTab::getEmisLiveChanVec]: ',&
            'Cannot rotate polarization without satellite altitude'
          call errorHalt(1)
        endif
      endif
    endif

  ! Check that frequency is in range
    nChanGr=size(allVar(indx)%myfrqGrd)
    if (minval(frq) < allVar(indx)%myfrqGrd(1) .or. &
        maxval(frq) > allVar(indx)%myfrqGrd(nChanGr)) then
      print *,'err[EmisTab::getEmisLiveChanVec]: ', &
        'Requested freq exceeds range of tabulated data'
      print *,'Requested:',minval(frq),maxval(frq), &
        '  Tabulated:',allVar(indx)%myfrqGrd(1),allVar(indx)%myfrqGrd(nChanGr)
      call errorHalt(1)
    endif

  ! Flag relevant channels (relevant polarizations of bracketing frequencies)
    iflg=0
    do l=1,nChanLoc
      call getFrqGridDat(frq(l),allVar(indx)%myfrqGrd,nStokes(indx), &
        j1(l),j2(l),w1(l),w2(l))
      do k=1,nStokes(indx)
        if (mapPol(k,pol(l)+1) /= 0.) then
          iflg(j1(l)+k-1)=1
          iflg(j2(l)+k-1)=1
        endif
      enddo
    enddo

  ! Convert flags into map from relevant channel list to channel grid
    nOn=1
    allVar(indx)%mapToGr=0
    do i=1,nChanGr
      allVar(indx)%mapToGr(nOn)=i    ! Overwrite slot until nOn increments
      nOn=nOn+iflg(i)   ! Increment to next slot when flagged channel
    enddo
    nOn=nOn-1

  ! Get emissivities for all relevant channels at once, for efficiency
    call getEmisTabChanVec(indx,intrpScheme,allVar(indx)%mapToGr(1:nOn),salin,temp,speed=speed, &
         direc=direc,angle=absAngle,emis=allVar(indx)%emisOn(1:nOn))

  ! Load emis data back to channel grid
    allVar(indx)%emisGrd=0.
    allVar(indx)%emisGrd(allVar(indx)%mapToGr(1:nOn))=allVar(indx)%emisOn(1:nOn)

  ! Execute frequency interpolation and polarization mapping
    emis=0.
    do l=1,nChanLoc
      ! Any irrelevant polarizations have emis and map zero
      do k=1,nStokes(indx)
        emis(l)=emis(l)+ &
          allVar(indx)%emisGrd(j1(l)+k-1)*mapPol(k,pol(l)+1)*w1(l)+ &
          allVar(indx)%emisGrd(j2(l)+k-1)*mapPol(k,pol(l)+1)*w2(l)
      enddo
    enddo

    return

  END SUBROUTINE getEmisLiveChanVec

  !--------------------------------------------------------------------------
  
  SUBROUTINE getEmisLiveChanSngl(indx,intrpScheme,frq,pol,salin,temp,speed, &
       direc,angle,emis,polRotates,satAlt)

  ! Subroutine to get data for a single specified frequency and polarization,
  ! where channelization is computed live (on the fly) from tabulated data.
  ! Optional argument satAlt should be present if and only if polRotates is 
  ! present and is true 

    integer, intent(in)             :: indx
    integer, intent(in)             :: intrpScheme
    real,    intent(in)             :: frq
    integer, intent(in)             :: pol
    real,    intent(in)             :: salin
    real,    intent(in)             :: temp
    real,    intent(in)             :: speed
    real,    intent(in)             :: direc
    real,    intent(in)             :: angle
    real,    intent(inout)          :: emis
    logical, intent(in),   optional :: polRotates
    real,    intent(in),   optional :: satAlt

    integer, dimension(1) :: myPol
    real,    dimension(1) :: myFrq,myemis

    if (nGauss_in(indx) > 1) then
      print *,'err[EmisTab::getEmisLiveChanSngl]: ',&
        'nGauss > 1 not supported in live mode;',nGauss_in(indx)
      call errorHalt(1)
    endif

    myFrq=frq
    myPol=pol
    call getEmisLiveChanVec(indx,intrpScheme,myFrq,myPol,salin,temp,speed, &
         direc=direc,angle=angle,emis=myemis,polRotates=polRotates,satAlt=satAlt)
    emis=myemis(1)

    return

  END SUBROUTINE getEmisLiveChanSngl

  !--------------------------------------------------------------------------

  SUBROUTINE boundGrid(name,grid,unig,step,val,i1,i2)
    
  ! Find the indices of the grid that bound the given value

    character(len=*),   intent(in)    :: name
    real, dimension(:), intent(in)    :: grid
    integer,            intent(in)    :: unig
    real,               intent(in)    :: step
    real,               intent(in)    :: val
    integer,            intent(inout) :: i1
    integer,            intent(inout) :: i2

    integer :: i,ngrid

    ngrid=size(grid)

    ! Grid has size 1
    if (ngrid == 1) then
      i1=1
      i2=1
    ! Grid is uniformly spaced; saves time without loop over if block
    elseif (unig == 1) then
      i1=floor((val-grid(1))/step)+1
      i2=i1+1
      if (i1 < 1) then
        if (outBndsExit) then
          print *,'err[EmisTab::boundGrid]: ', &
          name,' value',val,' is out of uniform range:',grid(1),grid(ngrid)
          call errorHalt(1)
        elseif (val < 0.) then
          print *,'err[EmisTab::boundGrid]: ', &
          name,' value',val,' is unphysical negative; min grid=',grid(1)
          call errorHalt(1)
        else
          i1=1
          i2=2
        endif
      endif
      if (i2 > ngrid) then
        if (val <= grid(ngrid)) then  ! Accommodate roundoff 
          i2=ngrid
          i1=i2-1
        elseif (outBndsExit) then
          print *,'err[EmisTab::boundGrid]: ', &
          name,' value',val,' is out of uniform range:',grid(1),grid(ngrid)
          call errorHalt(1)
        else
          i2=ngrid
          i1=i2-1
        endif
      endif
    else
    ! Grid is not uniformly spaced
      if ((val < grid(1) .or. val > grid(ngrid)) .and. outBndsExit) then
        print *,'err[EmisTab::boundGrid]: ', &
        name,' value',val,' is out of range:',grid(1),grid(ngrid)
        call errorHalt(1)
      elseif (val < 0.) then
        print *,'err[EmisTab::boundGrid]: ', &
        name,' value',val,' is unphysical negative; min grid=',grid(1)
        call errorHalt(1)
      endif
      
      do i=2,ngrid
         if (grid(i) >= val) exit
      enddo
      i2=min(i,ngrid)
      i1=i2-1
      
   endif
   
   return
   
  END SUBROUTINE boundGrid

  !--------------------------------------------------------------------------

  FUNCTION interpolate(type,grid,i1,i2,ia1,ia2,xin,vec)

    character(len=*),   intent(in)  :: type
    real, dimension(:), intent(in)  :: grid
    integer,            intent(in)  :: i1
    integer,            intent(in)  :: i2
    integer,            intent(in)  :: ia1
    integer,            intent(in)  :: ia2
    real,               intent(in)  :: xin
    real, dimension(:), intent(in)  :: vec
    real                            :: interpolate

    select case(type)
    case ('none')
      interpolate=vec(i1)
    case ('linear')
      interpolate=(vec(i1)*(grid(i2)-xin)+vec(i2)*(xin-grid(i1)))/ &
        (grid(i2)-grid(i1))
    case default
      print *,'err[EmisTab::interpolate]: ', &
      'Unsupported interpolation method: ',type
      call errorHalt(1)
    end select

    return

  END FUNCTION interpolate

  !--------------------------------------------------------------------------

  FUNCTION differentiate(type,grid,i1,i2,ia1,ia2,xin,vec)

    character(len=*),   intent(in)  :: type
    real, dimension(:), intent(in)  :: grid
    integer,            intent(in)  :: i1
    integer,            intent(in)  :: i2
    integer,            intent(in)  :: ia1
    integer,            intent(in)  :: ia2
    real,               intent(in)  :: xin
    real, dimension(:), intent(in)  :: vec
    real                            :: differentiate

    select case(type)
    case ('none')
      differentiate = 0.
    case ('linear')
      differentiate =(vec(i2)-vec(i1))/(grid(i2)-grid(i1))
    case default
      print *,'err[EmisTab::differentiate]: ', &
      'Unsupported differentiation method: ',type
      call errorHalt(1)
    end select

    return
    
  END FUNCTION differentiate
  
  !--------------------------------------------------------------------------
  
  SUBROUTINE emisTabDestroy(indx)
    
    integer, intent(in)  :: indx
    integer              :: ix
    
    ! Deallocate memory that can be allocated by emisTabInit, except
    ! for optional arguments of emisTabInit that are handled there
    
    if(first_init)then
       
       do ix=1,mxfd
          
          nullify(allVar(ix)%grSal,allVar(ix)%grT, &
               allVar(ix)%grSp,allVar(ix)%grDr,allVar(ix)%grAn)
          nullify(allVar(ix)%datArr)
          nullify(allVar(ix)%tmp2,allVar(ix)%tmp2ds, &
               allVar(ix)%tmp2dt,allVar(ix)%tmp2dd, &
               allVar(ix)%tmp3,allVar(ix)%tmp3ds,allVar(ix)%tmp3dt, &
               allVar(ix)%tmp4,allVar(ix)%tmp4ds,allVar(ix)%tmp4dt, & 
               allVar(ix)%tmp5,allVar(ix)%tmp5ds)
          nullify(allVar(ix)%tmp2dw,allVar(ix)%tmp3dw)
          nullify(allVar(ix)%mapToGr,allVar(ix)%emisOn,allVar(ix)%emisGrd)
          nullify(allVar(ix)%myfrqGrd,allVar(ix)%mypolGrd)
          
       enddo
       
       first_init=.false.
       
    endif
    
    if(indx == 0)then ! destroy all arrays
       
       do ix=1,mxfd
          
          if(associated(allVar(ix)%grSal))deallocate(allVar(ix)%grSal)
          if(associated(allVar(ix)%grT))deallocate(allVar(ix)%grT)
          if(associated(allVar(ix)%grSp))deallocate(allVar(ix)%grSp)
          if(associated(allVar(ix)%grDr))deallocate(allVar(ix)%grDr)
          if(associated(allVar(ix)%grAn))deallocate(allVar(ix)%grAn)
          if(associated(allVar(ix)%datArr))deallocate(allVar(ix)%datArr)
          if(associated(allVar(ix)%tmp2))deallocate(allVar(ix)%tmp2)
          if(associated(allVar(ix)%tmp2ds))deallocate(allVar(ix)%tmp2ds)
          if(associated(allVar(ix)%tmp2dt))deallocate(allVar(ix)%tmp2dt)
          if(associated(allVar(ix)%tmp2dw))deallocate(allVar(ix)%tmp2dw)
          if(associated(allVar(ix)%tmp2dd))deallocate(allVar(ix)%tmp2dd)
          if(associated(allVar(ix)%tmp3))deallocate(allVar(ix)%tmp3)
          if(associated(allVar(ix)%tmp3ds))deallocate(allVar(ix)%tmp3ds)
          if(associated(allVar(ix)%tmp3dt))deallocate(allVar(ix)%tmp3dt)
          if(associated(allVar(ix)%tmp3dw))deallocate(allVar(ix)%tmp3dw)
          if(associated(allVar(ix)%tmp4))deallocate(allVar(ix)%tmp4)
          if(associated(allVar(ix)%tmp4ds))deallocate(allVar(ix)%tmp4ds)
          if(associated(allVar(ix)%tmp4dt))deallocate(allVar(ix)%tmp4dt)
          if(associated(allVar(ix)%tmp5))deallocate(allVar(ix)%tmp5)
          if(associated(allVar(ix)%tmp5ds))deallocate(allVar(ix)%tmp5ds)
          if(associated(allVar(ix)%mapToGr))deallocate(allVar(ix)%mapToGr)
          if(associated(allVar(ix)%emisOn))deallocate(allVar(ix)%emisOn)
          if(associated(allVar(ix)%emisGrd))deallocate(allVar(ix)%emisGrd)
          if(associated(allVar(ix)%myfrqGrd))deallocate(allVar(ix)%myfrqGrd)
          if(associated(allVar(ix)%mypolGrd))deallocate(allVar(ix)%mypolGrd)
                    
       enddo
       
    else  ! destroy arrays pertaining to a particular file
       
       if(associated(allVar(indx)%grSal))deallocate(allVar(indx)%grSal,allVar(indx)%grT, &
            allVar(indx)%grSp,allVar(indx)%grDr,allVar(indx)%grAn)
       if(associated(allVar(indx)%datArr))deallocate(allVar(indx)%datArr)
       if(associated(allVar(indx)%tmp2))deallocate(allVar(indx)%tmp2,allVar(indx)%tmp2ds, &
            allVar(indx)%tmp2dt,allVar(indx)%tmp2dd, &
            allVar(indx)%tmp3,allVar(indx)%tmp3ds,allVar(indx)%tmp3dt, &
            allVar(indx)%tmp4,allVar(indx)%tmp4ds,allVar(indx)%tmp4dt, & 
            allVar(indx)%tmp5,allVar(indx)%tmp5ds)
       if(associated(allVar(indx)%tmp2dw))deallocate(allVar(indx)%tmp2dw,allVar(indx)%tmp3dw)
       if(associated(allVar(indx)%mapToGr))deallocate(allVar(indx)%mapToGr, &
            allVar(indx)%emisOn,allVar(indx)%emisGrd)
       if(associated(allVar(indx)%myfrqGrd))deallocate(allVar(indx)%myfrqGrd, &
            allVar(indx)%mypolGrd)
       
    endif
    
    nullify(datArr)
    
    return
    
  END SUBROUTINE emisTabDestroy

END MODULE EmisTab
