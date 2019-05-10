; Finds distribution of NH3 a priori in cloudy pixles
; Also finds distribution of DOFS in cloudy pixels
; Finally finds percentage of unp pixels that are clear and cloudy


; Read VIIRS  data
read_ncdf_jh,datav,filename='/project/p1913/esspa4_l2/l2_data/crisimg.snpp.20150422T0800.nc'
id_cloud = 2
id_band = 3  ; viirs reflectance band M04
viirs_subsets = ['All','Clear','Cloudy']
cris_lat_arr = datav.lat
cris_lon_arr = datav.lon
ncdf_cris = n_elements(cris_lat_arr)
viirs_refl = datav.viirs_refl
viirs_refl_sdev = datav.viirs_refl_sdev

dim_names = datav.dimensions.names
id_viirs_sub = where(dim_names eq 'viirs_subset')
nsubs = datav.dimensions.sizes(id_viirs_sub)
id_viirs_rband = where(dim_names eq 'viirs_refl_band')
nrbands = datav.dimensions.sizes(id_viirs_rband)

viirs_refl = reform(viirs_refl,nsubs,nrbands,ncdf_cris)
viirs_refl_sdev = reform(viirs_refl_sdev,nsubs,nrbands,ncdf_cris)

 ;Read CrIS data
read_ncdf_jh,data,filename='irsensors/crimss/run/aer_nc/retr.nh3.ir.nc'
cris_lat_arr = data.lat
cris_lon_arr = data.lon
ncdf = n_elements(cris_lat_arr)
nh3_id = data.attributes.global.mol[10]
nh3_len = data.attributes.global.mollen[10]

airmwt = 28.9
nh3mwt = 17.03
mrat = airmwt/nh3mwt

cris_sfc_arr = fltarr(ncdf)
; Convert from g/g to ppbv
cris_sfc_arr = mrat*data.profiles[nh3_id-1+nh3_len-1,*]*1.0e+9
xa_type = data.atmosclass
dofs = data.dofs

viirs_refl = viirs_refl[id_cloud,id_band,*]
viirs_refl_sdev = viirs_refl_sdev[id_cloud,id_band,*]

cndex = where(viirs_refl gt 0.3,complement=lndex,ncld)
undex = where(xa_type[cndex] eq 3,ncld_unp)
mndex = where(xa_type[cndex] eq 2,ncld_mod)
pndex = where(xa_type[cndex] eq 1,ncld_pol)

ncld = float(ncld)
print,'Fraction  of cloudy pixels for that belong  to each profile type (unp,mod,pol) '
print, float(ncld_unp)/ncld,float(ncld_mod)/ncld,float(ncld_pol)/ncld

undex = where(xa_type eq 3,nunp)
lndex = where(viirs_refl[undex] lt 0.3,nclear)
print, 'Fraction of unpolluted pixels that are clear'
print, float(nclear)/float(nunp)

end
