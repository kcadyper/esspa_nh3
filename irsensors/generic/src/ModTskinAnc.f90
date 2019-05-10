PROGRAM ModTskinAnc

   use ncdf_module
   use StateIndexModule
   use scene_io_module
   USE Dimensions, ONLY: MXLEV,MXCHAN,MXPARG
   use ncdfUtil, only:getDimLen, readVarNC, check
   use netcdf



!  Scene variables

   type(StateIndex_t)                  :: IG, NG
   real, dimension(:),allocatable      :: xstate
   real, dimension(:,:),allocatable    :: profiles
   integer                             :: nprf,nParRetr,nParG

!  Other variables

   character(len=200)                  :: fInOut,fRetv
   integer                             :: ncidIn
   integer                             :: ncidRetr
   real, dimension(:),allocatable      :: new_var
   integer                             :: ip
   integer                             :: kount
   real                                :: delta

!  Namelist

   namelist /runModL2/ fInOut, fRetv

!  Read namelist
   
   read(*,runModL2)

! Read state vector from existing file,edit and output

! Open Scene files
     print *,'fscene ',fInOut
     call openScene(ncidIn,fInOut,nprf=nprf,IG=IG,NG=NG,status='replace')
     !call openScene(ncidOut,fOut,status='replace')
     NParG   = getVectorLength(NG)
     print *,'ncidin ',ncidin,' nprf =',nprf
     print *,'IG%tskin ',IG%tskin,' NG%tskin =',ng%tskin, 'NParG ',nParG

! Read in retrieved SFCT and store
     call check( nf90_open(fRetv, nf90_nowrite, ncidRetr) )
      call getDimLen(ncidRetr, "nPar",  nParRetr)
     allocate(profiles(nParRetr,nprf))
     print *,'nParRetr ',nParRetr
     print *,'will read'
     call readVarNC(ncidRetr,"profiles",profiles)
     print *,'have read'

     call check( nf90_close(ncidRetr) )
     allocate(new_var(nprf))
     new_var(:) = profiles(102,:)

! Create new state vector
     allocate(xstate(nParG))
     print *,'NParG ',NParG

! For each profile; read in existing state vector and modify SFCT
     do ip=1,nprf
        call getScene(ncidIn,ip,x=xstate)
        xstate(IG%tskin) = new_var(ip)   !new_var comes from a different file
        !print *, ip,xstate(IG%tskin) 
        call replaceScene(ncidIn,ip,x=xstate)
     enddo
     call closeScene(ncidIn)

     deallocate (new_var)

     end
     
