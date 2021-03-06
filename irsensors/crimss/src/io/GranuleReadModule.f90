Module  GranuleReadModule

    use ncdfUtil, only:getDimLen, readVarNC, check
    use netcdf
    use StateIndexModule
    use constants, only: g0,Rd

    implicit none
    public :: GranuleRead


  contains

   subroutine GranuleRead (fname,ncid,MolID,nAir_pres, nchmw,nSurf_wnum_ir,vCoordTyp,nprf, &
              air_pres,sfcgrid, lat_arr,lon_arr,land_frac_arr, &
               temp_arr,wv_arr,emrf_arr,emrf_qc,sfcp_arr,sfct_arr,sfcp_index)

    character(len=*),   intent(in)      :: fname 
    integer,            intent(inout)   :: ncid
    integer,dimension(:),intent(inout)  :: MolID
    integer,            intent(out)     :: nAir_pres,nchmw,nsurf_wnum_ir,nprf
    character(len=*),   intent(out)     :: vCoordTyp
    real*4,dimension(:),intent(inout)   :: air_pres
    real*4,dimension(:),intent(inout)   :: sfcgrid
    real*4,dimension(:), intent(inout)  :: lat_arr,lon_arr,land_frac_arr,sfcp_arr,sfct_arr
    integer,dimension(:), intent(inout) :: sfcp_index
    real*4,dimension(:,:),intent(inout) :: temp_arr,wv_arr,emrf_arr,emrf_qc


! Local variables


    integer                                  :: nAtrack,nXtrack, nFov,ncid_orig
    integer                                  :: nCld_top_pres,nAir_pres_wv,wv_index

    real*4 , allocatable, dimension(:,:,:)   :: fov_lon
    real*4 , allocatable, dimension(:,:,:)   :: fov_lat

    real*4 , allocatable, dimension(:,:,:)   :: fov_land_frac
    real*4 , allocatable, dimension(:,:,:)   :: fov_surf_alt
    real*4 , allocatable, dimension(:,:,:,:) :: cld_frac

    integer, allocatable, dimension(:,:)     :: air_pres_nsurf
    real*4                                   :: for_sfcp
    real*4, allocatable, dimension(:,:,:)    :: fov_sfcp
    integer                                   ::for_sfcp_index
    integer, allocatable, dimension(:,:,:)    ::fov_sfcp_index

    real*4 , allocatable, dimension(:,:,:)   :: air_temp
    real*4 , allocatable, dimension(:,:,:,:) :: fov_air_temp
    real*4 , allocatable, dimension(:,:)     :: surf_air_temp
 
    real*4 , allocatable, dimension(:,:)     :: surf_temp
    real*4 , allocatable, dimension(:,:,:)   :: fov_surf_temp

    real*4 , allocatable, dimension(:,:,:)   :: spec_hum
    real*4 , allocatable, dimension(:,:,:,:) :: fov_spec_hum
    real*4 , allocatable, dimension(:,:)     :: surf_spec_hum

    real*4 , allocatable, dimension(:,:,:)   :: surf_ir_wnum
    real*4 , allocatable, dimension(:,:,:)   :: surf_ir_emis,surf_ir_emis_qc
    real*4 , allocatable, dimension(:,:,:,:) :: fov_surf_ir_emrf
    real*4 , allocatable, dimension(:,:,:,:) :: fov_surf_ir_emrf_qc

    integer                                  :: ix,iy,iv,isurf

    ncid_orig=999 
    print *, 'check on ', trim(fname)
    if (ncid.eq.0) ncid_orig=0
    call check( nf90_open(fname, nf90_nowrite, ncid) )
    

    MolId(1) = 1
    MolId(2:) = 0
    vCoordtyp = 'pressureStd'
    nchmw=0

    ! read dimensions
    call getDimLen(ncid, "atrack",  nAtrack)
    call getDimLen(ncid, "xtrack",  nXtrack)
    call getDimLen(ncid, "fov",  nFov)
    call getDimLen(ncid, "air_pres",  nAir_pres)
    call getDimLen(ncid, "air_pres_h2o",  nAir_pres_wv)
    call getDimLen(ncid, "cld_top_pres",  nCld_top_pres)
    call getDimLen(ncid, "surf_wnum_ir",  nSurf_wnum_ir)


! Return if initial call
    nprf = nAtrack*nXtrack*nFov
    if (ncid_orig.eq.0) then
       call check( nf90_close(ncid) )
       return
    endif

! Read all desired variables
    allocate(fov_lon(nFov,nXtrack,nAtrack),fov_lat(nFov,nXtrack,nAtrack))
    call readVarNC(ncid,"fov_lon", fov_lon)
    call readVarNC(ncid,"fov_lat", fov_lat)

    allocate(fov_land_frac(nFov,nXtrack,nAtrack),fov_surf_alt(nFov,nXtrack,nAtrack),cld_frac(nCld_top_pres,nFov,nXtrack,nAtrack))
    call readVarNC(ncid,"fov_land_frac", fov_land_frac)
    call readVarNC(ncid,"fov_surf_alt", fov_surf_alt)
    call readVarNC(ncid,"cld_frac", cld_frac)

    call readVarNC(ncid,"air_pres", air_pres)
! Convert to mbar
    air_pres = air_pres/100.
! Do not need to read in h2o levels, as they correspond to the last 66 levels of
! the regular levels


    allocate(air_pres_nsurf(nXtrack,nAtrack),fov_sfcp(nFov,nXtrack,nAtrack),fov_sfcp_index(nFov,nXtrack,nAtrack))
    call readVarNC(ncid,"air_pres_nsurf", air_pres_nsurf)

    allocate(air_temp(nAir_pres,nXtrack,nAtrack),surf_air_temp(nXtrack,nAtrack))
    allocate (fov_air_temp (nAir_pres,nFov,nXtrack,nAtrack))
    call readVarNC(ncid,"air_temp", air_temp)
    call readVarNC(ncid,"surf_air_temp", surf_air_temp)

    allocate(surf_temp(nXtrack,nAtrack),fov_surf_temp(nFov,nXtrack,nAtrack))
    call readVarNC(ncid,"surf_temp", surf_temp)

    allocate(spec_hum(nAir_pres_wv,nXtrack,nAtrack),surf_spec_hum(nXtrack,nAtrack))
    allocate(fov_spec_hum(nAir_pres,nFov,nXtrack,nAtrack))
    call readVarNC(ncid,"spec_hum", spec_hum)
    call readVarNC(ncid,"surf_spec_hum", surf_spec_hum)

    allocate(surf_ir_wnum(nSurf_wnum_ir,nXtrack,nAtrack),surf_ir_emis(nSurf_wnum_ir,nXtrack,nAtrack))
    allocate(surf_ir_emis_qc(nSurf_wnum_ir,nXtrack,nAtrack))
    allocate(fov_surf_ir_emrf(2*nSurf_wnum_ir,nFov,nXtrack,nAtrack))
    allocate(fov_surf_ir_emrf_qc(nSurf_wnum_ir,nFov,nXtrack,nAtrack))

    call readVarNC(ncid,"surf_ir_wnum", surf_ir_wnum)
! Surface emissivity grid is the same for all pixels
    sfcgrid = surf_ir_wnum(:,1,1)

    call readVarNC(ncid,"surf_ir_emis", surf_ir_emis)
    call readVarNC(ncid,"surf_ir_emis_qc", surf_ir_emis_qc)

    call check( nf90_close(ncid) )

! Build lat,lon,pland output arrays
    lat_arr = reshape (fov_lat,(/ nprf /))
    lon_arr = reshape (fov_lon,(/ nprf /))
    land_frac_arr = reshape (fov_land_frac,(/ nprf /))

!Build surface temperature output arrays: assume SFCT going into retrieval is uniform over FOR
! thus the need for SFCT retrieval

    do ix=1,nXtrack
       do iy = 1,nAtrack
          fov_surf_temp(1:nFov,ix,iy) = surf_temp(ix,iy)
       end do
    end do
    sfct_arr = reshape(fov_surf_temp,(/ nprf /))
 

! Build surface pressure arrays
    do ix=1,nXtrack
       do iy=1,nAtrack
          for_sfcp = air_pres(air_pres_nsurf(ix,iy))
          for_sfcp_index = air_pres_nsurf(ix,iy)
          do iv=1,nFov
             fov_sfcp(iv,ix,iy) = for_sfcp*exp(-g0*(fov_surf_alt(iv,ix,iy)-fov_surf_alt(5,ix,iy))/(Rd*surf_air_temp(ix,iy)))
             fov_sfcp_index(iv,ix,iy) = for_sfcp_index
          end do
       end do
    end do
    sfcp_arr = reshape(fov_sfcp,(/ nprf /))
    sfcp_index = reshape(fov_sfcp_index, (/ nprf /))

! Build temperature arrays
    do ix=1,nXtrack
       do iy=1,nAtrack
          isurf= air_pres_nsurf(ix,iy)
          do iv=1,nFov
             fov_air_temp(:,iv,ix,iy) = air_temp(:,ix,iy)
             fov_air_temp(isurf,iv,ix,iy) = surf_air_temp(ix,iy)
          end do
       end do
    end do
    temp_arr = reshape (fov_air_temp,(/ nAir_pres,nprf /))
    
! Build wv  arrays
    wv_index = nAir_pres-nAir_pres_wv+1
    do ix=1,nXtrack
       do iy=1,nAtrack
          isurf= air_pres_nsurf(ix,iy)
          do iv=1,nFov
             fov_spec_hum(wv_index:nAir_pres,iv,ix,iy) = spec_hum(:,ix,iy)
             fov_spec_hum(isurf,iv,ix,iy) = surf_spec_hum(ix,iy)
             fov_spec_hum(1:wv_index-1,iv,ix,iy) = spec_hum(1,ix,iy)
          end do
       end do
    end do
    wv_arr = reshape (fov_spec_hum,(/ nAir_pres,nprf /))

! Build emiss/refl  arrays
    do ix=1,nXtrack
       do iy=1,nAtrack
          do iv=1,nFov
             fov_surf_ir_emrf(1:2*nSurf_wnum_ir-1:2,iv,ix,iy) = surf_ir_emis(:,ix,iy)
             fov_surf_ir_emrf(2:2*nSurf_wnum_ir:2,iv,ix,iy) = 1.0-surf_ir_emis(:,ix,iy)
             fov_surf_ir_emrf_qc(1:nSurf_wnum_ir,iv,ix,iy) = surf_ir_emis_qc(:,ix,iy)
          end do
       end do
    end do
    emrf_arr = reshape (fov_surf_ir_emrf,(/ 2*nSurf_wnum_ir,nprf /))
    emrf_qc = reshape (fov_surf_ir_emrf_qc,(/ nSurf_wnum_ir,nprf /))


    deallocate (fov_lon,fov_lat,fov_land_frac,fov_surf_alt,cld_frac)
    deallocate (air_pres_nsurf,air_temp,surf_air_temp,surf_temp,spec_hum,surf_spec_hum)
    deallocate (surf_ir_wnum,surf_ir_emis,surf_ir_emis_qc)
    deallocate (fov_surf_temp,fov_sfcp,fov_air_temp,fov_spec_hum,fov_surf_ir_emrf,fov_surf_ir_emrf_qc)

   end subroutine GranuleRead
end module GranuleReadModule
