; Plots CrIS NH3, VIIRS derived  cloudy reflectance and  VIIRS cirrus clouds
; This version reads in coordinates of pixel edges from L2 file

; Set up plot
plot_type = 'PS'
if (plot_type eq 'PS') then begin
   set_plot,'PS'
   device,/color
   device,/Helvetica,/Bold
   device,/landscape
   !p.font = 0
endif


loadct, 13

!x.thick = 1
!y.thick = 1
!p.thick = 1
!p.charsize = 1.0

!p.multi = [0,1,3]

airmwt = 28.9
nh3mwt = 17.03
mrat = airmwt/nh3mwt


; Read CrIS data
read_ncdf_jh,data,filename='irsensors/crimss/run/aer_nc/retr.nh3.ir.nc'
cris_lat_arr = data.lat
cris_lon_arr = data.lon
ncdf = n_elements(cris_lat_arr)
nh3_id = data.attributes.global.mol[10]
nh3_len = data.attributes.global.mollen[10]

cris_sfc_arr = fltarr(ncdf)
; Convert from g/g to ppbv
cris_sfc_arr = mrat*data.profiles[nh3_id-1+nh3_len-1,*]*1.0e+9
xa_type = data.atmosclass
chisq = data.chisq


a = findgen(17)*(!pi*2/16.)
clevels_n = [0.,1.0,2.0,findgen(8)*3.0+5.,40.]
cores_n = indgen(12)*20+35


; Read VIIRS  data
read_ncdf_jh,datav,filename='/project/p1913/esspa4_l2/l2_data/crisimg.snpp.20150422T0800.nc'
viirs_subsets = ['All','Clear','Cloudy']
id_band_refl = 3  ; viirs reflectance band M04
id_band_m15 = 3
id_band_m16 = 4
id_viirs_sel=0  ; all
cris_lat_arr = datav.lat
cris_lon_arr = datav.lon
ncdf_cris = n_elements(cris_lat_arr)
viirs_refl = datav.viirs_refl
viirs_refl_sdev = datav.viirs_refl_sdev
viirs_bt = datav.viirs_bt

dim_names = datav.dimensions.names
id_viirs_sub = where(dim_names eq 'viirs_subset')
nsubs = datav.dimensions.sizes(id_viirs_sub)
id_viirs_rband = where(dim_names eq 'viirs_refl_band')
nrbands = datav.dimensions.sizes(id_viirs_rband)
id_viirs_bband = where(dim_names eq 'viirs_emis_band')
nbbands = datav.dimensions.sizes(id_viirs_bband)

viirs_refl = reform(viirs_refl,nrbands,nsubs,ncdf_cris)
viirs_refl_sdev = reform(viirs_refl_sdev,nrbands,nsubs,ncdf_cris)
viirs_bt = reform(viirs_bt,nbbands,nsub,ncdf_cris)

a = findgen(17)*(!pi*2/16.)
clevels_r = findgen(11)*0.1
cores_r = indgen(11)*22+35

; Read in pixel bounds
read_ncdf_jh,l2,filename='/project/p1913/esspa4_l2/l2_data/SNDR.SNPP.CRIMSS.20150422T0800.m06.g081.L2_RET_CLIMCAPS_NSR.std.v01_23_01.J.180327125004.nc'
lat_bnds = l2.fov_lat_bnds
lon_bnds = l2.fov_lon_bnds

dim_names = l2.dimensions.names
id_poly = where(dim_names eq 'fov_poly')
npoly = l2.dimensions.sizes[id_poly]

lat_bnds = reform(lat_bnds,npoly,ncdf_cris)
lon_bnds = reform(lon_bnds,npoly,ncdf_cris)


; Plot CrIS coords
!p.position = [0.20,0.15,0.83,0.90]

min_lat = 20.0
max_lat = 50.0
min_lon = 55.0
max_lon = 95.0

nlab_lon = 5
lon_lab = indgen(nlab_lon)*10.0+min_lon
nlab_lat = 7
lat_lab = indgen(nlab_lat)*5.0+min_lat


!x.style=8
!y.style=8

p_arr = fltarr(3,4)
p_arr[0,*] = [0.05,0.4,0.32,0.80]
p_arr[1,*] = [0.37,0.4,0.64,0.80]
p_arr[2,*] = [0.69,0.4,0.96,0.80]

; Calculate BT diff (M15-M16) and use in cirrus cloud test

; Map all NH3, VIIRS 0.55 refl,  NH3 clear based on VIIRS refl
is = id_viirs_sel
titulo = ['NH3','VIIRS 0.55 Refl' ,'NH3 ']
device,filename='plot_cris_viirs_comb.ps'
for ip=0,2 do begin

   !p.position = p_arr[ip,*]
   map_set,0.,90.,limit=[min_lat,min_lon,max_lat,max_lon],$
       /cylin,/continent,mlinethick=3,/noborder,title = titulo[ip],/advance

   axis,xaxis=0,xticks=nlab_lon-1,xtickv=[findgen(nlab_lon)/(nlab_lon-1)],xtickname=string(lon_lab,format='(f6.2)'),xtitle='Longitude'
   axis,xaxis=1,xticks=1,xtickname=[' ',' '],xtickv=[min_lon,max_lon]

   axis,yaxis=0,yticks=nlab_lat-1,ytickv=[findgen(nlab_lat)/(nlab_lat-1)],ytickname=string(lat_lab,format='(f4.1)'),ytitle='Latitude'
   axis,yaxis=1,yticks=1,ytickname=[' ',' '],ytickv=[min_lat,max_lat]

   case ip of
   0: begin
      neta = npoly[0]
      for ic=0,ncdf_cris-1 do begin
	 lat_eta = fltarr(neta)
	 lon_eta = fltarr(neta)
	 for ie=0,neta-1 do begin
	    lat_eta[ie] = lat_bnds[ie,ic]
	    lon_eta[ie] = lon_bnds[ie,ic]
	endfor
        cndex = where(clevels_n gt cris_sfc_arr[ic],ng)
	if (ng gt 0) then id = cndex[0]-1 else id = 11 
	if (cris_sfc_arr[ic] gt 0.0) then begin
	   oplot,lon_eta,lat_eta
	   polyfill,lon_eta,lat_eta,color=cores_n[id]
	endif
      endfor
      draw_scale_bar,'ppbv',[.05,.32,.32,.34],cores_n[0:11],clevels_n[0:11],['(f3.0)'],[''],[0.0,-0.02],charsize=0.7
   end
   1:begin
      neta = npoly[0]
      for ic=0,ncdf_cris-1 do begin
	 lat_eta = fltarr(neta)
	 lon_eta = fltarr(neta)
	 for ie=0,neta-1 do begin
	    lat_eta[ie] = lat_bnds[ie,ic]
	    lon_eta[ie] = lon_bnds[ie,ic]
	endfor
	cndex = where(clevels_r gt viirs_refl[id_band_refl,id_viirs_sel,ic],ng)
	if (ng gt 0) then id = cndex[0]-1 else id = 10 
	if (viirs_refl[id_band_refl,id_viirs_sel,ic] gt 0.0) then begin
	   oplot,lon_eta,lat_eta
	   polyfill,lon_eta,lat_eta,color=cores_n[id]
	endif
      endfor
      draw_scale_bar,'refl',[.37,.32,.64,.34],cores_r[0:10],clevels_r[0:10],['(f3.1)'],[''],[0.0,-0.02],charsize=0.7
   end
   2: begin
      neta = npoly[0]
      for ic=0,ncdf_cris-1 do begin
	 lat_eta = fltarr(neta)
	 lon_eta = fltarr(neta)
	 for ie=0,neta-1 do begin
	    lat_eta[ie] = lat_bnds[ie,ic]
	    lon_eta[ie] = lon_bnds[ie,ic]
	endfor
        cndex = where(clevels_n gt cris_sfc_arr[ic],ng)
	if (ng gt 0) then id = cndex[0]-1 else id = 11 
	if (cris_sfc_arr[ic] gt 0.0 AND viirs_refl[id_band_refl,id_viirs_sel,ic] lt 0.3) then begin
	   oplot,lon_eta,lat_eta
	   polyfill,lon_eta,lat_eta,color=cores_n[id]
	endif
      endfor
      draw_scale_bar,'ppbv',[.69,.32,.96,.34],cores_n[0:11],clevels_n[0:11],['(f3.0)'],[''],[0.0,-0.02],charsize=0.7
   end
   endcase
endfor
device,/close_file



end
