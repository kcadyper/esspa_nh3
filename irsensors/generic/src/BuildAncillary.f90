PROGRAM BuildAncillary

   use ncdf_module
   use StateIndexModule
   use scene_io_module
   use GranuleReadModule
   USE Dimensions, ONLY: MXLEV,MXCHAN,MXPARG
   USE VertCoord, ONLY: mxCoordTyp,Pcoord,Scoord,Hcoord, &
     putSigmDefin,putHybrDefin,getSigmPres,getHybrPres
   
   implicit none
   
!-------------------------------------------------------------------------------
   
!  Constants

   integer, parameter                  :: mxStrLen=12
   character(len=mxStrLen)             :: NULL='NULL'
   integer, parameter                  :: lentype=12
   integer                             :: iuC=11
   real, parameter                     :: max_freq=2600.

!  Scene variables

   character(len=8)                    :: xid
   character(len=20)                   :: casename='Regridded'
   real, dimension(MXPARG)             :: x
   real, dimension(MXCHAN)             :: emmw
   real, dimension(MXCHAN*2)           :: EmRf
   real, dimension(:),allocatable      :: eaa, eia
   real, dimension(MXLEV)              :: pref
   real, dimension(MXLEV)              :: prefIn
   real, dimension(MXCHAN)             :: wvn
   integer                             :: nprofiles
   integer                             :: NParG
   integer                             :: nlev
   integer                             :: nlevIn
   integer                             :: nchmw
   integer                             :: nchir
   integer                             :: neang
   integer, dimension(maxMol)          :: MolID
   integer, dimension(6)               :: time
   type(StateIndex_t)                  :: IG, NG
   real                                :: pTop
   character(len=mxCoordTyp)           :: vCoordTyp
   character(len=mxCoordTyp)           :: vCoordTypIn
   real, dimension(MXLEV)              :: sigCoord
   real, dimension(MXLEV)              :: hybCoordA
   real, dimension(MXLEV)              :: hybCoordB
   real, dimension(MXLEV)              :: sigCoordIn
   real, dimension(:),   allocatable   :: plog,dlog
   real, dimension(:),   allocatable   :: presOut,plogOut
   real, dimension(:),   allocatable   :: xlog
   real, dimension(:,:), allocatable   :: trns

! Input scene variables not passed directly to output

      real*4, dimension(:),allocatable :: lat_arr, lon_arr,pland_arr      
      real*4, dimension(:),allocatable :: sfcp_arr,sfct_arr
      integer, dimension(:),allocatable :: sfcp_index
      real*4, dimension(:,:),allocatable :: temp_arr, wv_arr,emrf_arr,emrf_qc
      real*4, dimension(:,:),allocatable :: emrf_out,emrf_qc_out
   
!  Other variables
   
   character(len=200)                  :: coordFile
   character(len=200)                  :: fIn
   character(len=200)                  :: fOut,fOutEmRf
   integer                             :: ncidIn
   integer                             :: ncidOut
   integer                             :: i,ii,k,n,nn,m,idum,last,il
   integer                             :: nMol
   logical                             :: exists
   integer                             :: ncStatus
   integer                             :: nCoordPar

   integer                             :: nout

!  Namelist

   namelist /runGetL2/ vCoordTyp,coordFile,fIn,fOut,fOutEmRf,nout

!-------------------------------------------------------------------------------

!  Read namelist
   
   read(*,runGetL2)

!  Read file with new coordinates
   open(iuC,file=coordFile,action='READ')
   read(iuC,*)nlev,nCoordPar
   select case (TRIM(vCoordTyp))
   case (Pcoord)
      do k=1,nlev
         read(iuC,*)idum,pref(k)
      enddo
   case (Scoord)
      do k=1,nlev
         read(iuC,*)idum,sigCoord(k)
      enddo
   case (Hcoord)
      do k=1,nlev
         read(iuC,*)idum,hybCoordA(k),hybCoordB(k)
      enddo
   case default
      print *,'ERROR in BuildAncillary: '// &
      'Unrecognized vertical coordinate type requested: ', &
      TRIM(vCoordTyp)
      call exit(1)
   end select
   close(iuC)

!  Get contents of input files

   ncidIn=0
   call GranuleRead(fname=fIn,ncid=ncidIn,MolID=MolID,nAir_pres=nlevIn,&
      nchmw=nchmw,nSurf_wnum_ir=nchir,vCoordTyp=vCoordTypIn,nprf=nprofiles, &
      air_pres=prefIn(1:nLevIn),sfcgrid=wvn(1:2*nchir),&
      lat_arr=lat_arr,lon_arr=lon_arr,land_frac_arr=pland_arr, &
      temp_arr=temp_arr,wv_arr=wv_arr,emrf_arr=emrf_arr,emrf_qc=emrf_qc, &
      sfcp_arr=sfcp_arr, sfct_arr=sfct_arr,sfcp_index=sfcp_index)

   allocate(lat_arr(nprofiles),lon_arr(nprofiles),pland_arr(nprofiles))
   allocate(temp_arr(nLevIn,nprofiles),wv_arr(nLevIn,nprofiles))
   allocate(emrf_arr(2*nchir,nprofiles),emrf_qc(nchir,nprofiles))
   allocate(sfct_arr(nprofiles),sfcp_arr(nprofiles),sfcp_index(nprofiles))

   call GranuleRead(fname=fIn,ncid=ncidIn,MolID=MolID,nAir_pres=nlevIn, &
      nchmw=nchmw,nSurf_wnum_ir=nchir,vCoordTyp=vCoordTypIn,nprf=nprofiles, &
      air_pres=prefIn(1:nLevIn),sfcgrid=wvn(1:2*nchir),&
      lat_arr=lat_arr,lon_arr=lon_arr,land_frac_arr=pland_arr, &
      temp_arr=temp_arr,wv_arr=wv_arr,emrf_arr=emrf_arr,emrf_qc=emrf_qc, &
      sfcp_arr=sfcp_arr, sfct_arr=sfct_arr,sfcp_index=sfcp_index)

! May want to change this later: need to set these variables for simulated runs,
! but L2 does not provide them.  
   neang = 1
   allocate(eia(neang),eaa(neang))
   eia(1:neang) = 0.
   eaa(1:neang) = 0.

   if(nlevIn > MXLEV)then
     print *,'ERROR in BuildAncillary: nlevIn exceeds MXLEV;',nlevIn,MXLEV
     stop
   elseif(nchmw > MXCHAN)then
     print *,'ERROR in BuildAncillary: nchmw exceeds MXCHAN;',nchmw,MXCHAN
     stop
   elseif(nchir > MXCHAN*2)then
     print *,'ERROR in BuildAncillary: nchir exceeds MXCHAN*2;',nchir,MXCHAN*2
     stop
   endif

   if (nlevIn > 0) then
      select case (TRIM(vCoordTypIn))
      case (Pcoord)
         ALLOCATE(plog(nlevIn),dlog(nlevIn-1))
         plog=log(prefIn(1:nlevIn))
         dlog=plog(2:nlevIn)-plog(1:nlevIn-1)
         pTop=prefIn(1) ! Treat the top of the standard grid as sigma=0 level
      case default
         print *,'ERROR in BuildAncillary: '// &
         'Unrecognized vertical coordinate type in file: ', &
         TRIM(vCoordTypIn)
         call exit(1)
      end select
   else
      print *,'ERROR in BuildAncillary: Expecting non-zero nlev in file: ',nlevIn
      call exit(1)
   endif

   nMol=getNMol(MolID)
   NG%temp=nlev
   NG%tskin=1
   NG%psfc=1
   do m=1,nMol
      NG%mol(m)=nlev
   enddo
   NG%mol(nmol+1:maxmol)=0
   NG%wind = 0
   NG%cldliq = 0
   NG%cldice = 0
   NG%emMW = 0
   NG%emRfIR = 0
   IG=genIndices(NG)
   NParG   = getVectorLength(NG)
   if (NParG > MXPARG) then
      print *,'ERROR in BuildAncillary: NParG exceeds dimension: ',NParG,MXPARG
      call exit(1)
   endif

!  Check for existing file

   INQUIRE(FILE=fOut,EXIST=exists)
   if (exists) then
       print *,'WARNING from BuildAncillary: Overwriting file'
       print *,fOut
   endif

   ALLOCATE(presOut(nlev),plogOut(nlev))
   ALLOCATE(xlog(nlev),trns(nlev,nlevIn))

!  Open netCDF scene file and write attributes

   call openScene(ncidOut,fOut,casename=casename,nAng=neang, &
      vCoordTyp=vCoordTyp, nLevel=nlev,nSfcGrid=nchir, &
      IG=IG,NG=NG,MolID=MolID,status='new')

   if (nchir > 0) then
      last = 0
      do ii=1,nchir
         if (wvn(ii).LT.1.0e6) then
            last = last+1
         else
            exit
        endif
      end do
      call writeNcdfAttr(ncidOut,attrName='sfcGrid',  attr=wvn(1:last))
   endif

   select case (TRIM(vCoordTyp))
   case (Pcoord)
      call writeNcdfAttr(ncidOut,attrName='standardPressureGrid', &
        attr=pref(1:nlev))
   case (Scoord)
      call writeNcdfAttr(ncidOut,attrName='SigmaCoord', &
        attr=sigCoord(1:nlev))
      call writeNcdfAttr(ncidOut,attrName='SigmaTopPres', &
        attr=pTop)
      call putSigmDefin(nlev,ptop,sigCoord(1:nlev))
   case (Hcoord)
      call writeNcdfAttr(ncidOut,attrName='HybridCoordA', &
        attr=hybCoordA(1:nlev))
      call writeNcdfAttr(ncidOut,attrName='HybridCoordB', &
        attr=hybCoordB(1:nlev))
      call putHybrDefin(nlev,hybCoordA(1:nlev),hybCoordB(1:nlev))
   end select
          
!  Copy data
 
   print *,'Regridding',nout,' profiles'

   do ii = 1,nout  ! Loop over scenes


     !---make transformation matrix
     if (TRIM(vCoordTypIn) == Pcoord .and. TRIM(vCoordTyp) == Scoord) then
        presOut(1:nlev)=getSigmPres(sfcp_arr(ii))
        plogOut(1:nlev)=log(presOut(1:nlev))
     elseif (TRIM(vCoordTypIn) == Pcoord .and. TRIM(vCoordTyp) == Hcoord) then
        presOut(1:nlev)=getHybrPres(sfcp_arr(ii))
        plogOut(1:nlev)=log(presOut(1:nlev))
     else
        print *,'ERROR in BuildAncillary: Unsupported regrid: ', &
           TRIM(vCoordTypIn),' to ',TRIM(vCoordTyp)
        call exit(1)
     endif
     if (TRIM(vCoordTypIn) == Pcoord .and. &
         (TRIM(vCoordTyp) == Scoord .or. TRIM(vCoordTyp) == Hcoord)) then
        trns=0.
        do i=1,nlev
           do n=2,sfcp_index(ii)
              if ((plog(n) >= plogOut(i)) .OR. (n == sfcp_index(ii))) then
                 trns(i,n)=(plogOut(i)-plog(n-1))/dlog(n-1)
                 trns(i,n-1)=(plog(n)-plogOut(i))/dlog(n-1)
                 EXIT
              endif
           enddo
        enddo
        last=0
        do n=1,nlevIn
           if (SUM(trns(:,n)) > 0.)last=n
        enddo
     endif

!    Map data into state vector
     if (NG%temp > 0) then
        xlog=MATMUL(trns(:,1:last),log(temp_arr(1:last,ii)))
        x(IG%temp:IG%temp+NG%temp-1)=exp(xlog)
     endif
     do m=1,nMol  ! restricts to molecules with NG%mol(mm)>0
!    Flip negative values to positive (need better approach)
        forall (il=1:last, wv_arr(il,ii) < 0.)
              wv_arr(il,ii) = abs(wv_arr(il,ii))
        endforall
        xlog=MATMUL(trns(:,1:last),log(wv_arr(1:last,ii)))
        x(IG%mol(m):IG%mol(m)+NG%mol(m)-1)=exp(xlog)
     enddo
     if (NG%tskin > 0) x(IG%tskin:IG%tskin+NG%tskin-1)= sfct_arr(ii)
     if (NG%psfc > 0) x(IG%psfc:IG%psfc+NG%psfc-1)=  sfcp_arr(ii)

!    Write data to NetCDF scene file
     call putScene(ncidOut,xid=xid, &
        lat=lat_arr(ii),lon=lon_arr(ii),pland=pland_arr(ii))
     !if (neang > 0) call putScene(ncidOut,eaa=eaa(1:neang),eia=eia(1:neang))
     if (NParG > 0) call putScene(ncidOut,x=x(1:NParG))
     if (nchmw > 0) call putScene(ncidOut,emmw=emmw(1:nchmw))
     !if (nchir > 0) call putScene(ncidOut,EmRf=EmRf_arr(1:nchir*2,ii))
     if (TRIM(vCoordTyp) == Scoord .or. TRIM(vCoordTyp) == Hcoord) &
        call putScene(ncidOut,pressure=presOut(1:nlev))

   enddo  ! End loop over scenes

   call closeScene(ncidOut)

   deallocate(plog,dlog,presOut,plogOut)
   deallocate(xlog,trns)

   print *,'Wrote',nout,'profiles'

! Write emissivity in a separate file
   INQUIRE(FILE=fOutEmRf,EXIST=exists)
   if (exists) then
       print *,'WARNING from BuildAncillary: Overwriting file'
       print *,fOutEmRf
   endif

! Only store useful frequencies for now (LE 2600.)
   if (nchir > 0) then
       last = 0
       do ii=1,nchir
          if (wvn(ii).LE.max_freq) then
             last = last+1
          else
             exit
         endif
       end do
   endif


   allocate(emrf_out(2*last,nout),emrf_qc_out(last,nout))
   emrf_out(1:last,:) = emrf_arr(1:2*last-1:2,1:nout)
   emrf_out(last+1:2*last,:) = emrf_arr(2:2*last:2,1:nout)
   emrf_qc_out = emrf_qc(1:last,1:nout)

! Create and fill emissivity file
   call openNcdfFile(ncidOut, fOutEmRf, status='new',unlimited_dim_name='nbkg')
   call writeNcdfDim(ncidOut,dimName='NGEmIR',dimLen=2*last)
   call writeNcdfDim(ncidOut,dimName='nQC',dimLen=last)

   call writeNcdfAttr(ncidOut,attrName='wavenumbers', &
     attr=wvn(1:last))

   call writeNcdfAttr(ncidOut,attrName='FileStructure',  attr=1)
   call writeNcdfAttr(ncidOut,attrName='nemir',attr=last)
   call writeNcdfAttr(ncidOut,attrName='TransFlag',  attr=0)



     call writeNcdfData(ncidOut,emrf_out,varname='EmIRBackground', &
           varLenName=(/ 'NGEmIR','nbkg  '/), varLongName='IR Emissivity and Reflectivity Vector ')

     call writeNcdfData(ncidOut,emrf_qc_out,varname='qc', &
           varLenName=(/ 'nQC ','nbkg'/) , varLongName='EmRf QC')

   call closeNcdfFile(ncidOut)


END PROGRAM BuildAncillary
