; plot CrIS  coordinates and values for all observations in  a file
; This version reads in coordinates of pixel edges from L2 file

; Set up plot
plot_type = 'PS'
if (plot_type eq 'PS') then begin
   set_plot,'PS'
   device,/color
   device,/Helvetica,/Bold
   device,filename='plot_cris_map.ps'
   !p.font = 0
endif


loadct, 13

!x.thick = 1
!y.thick = 1
!p.thick = 1
!p.charsize = 1.0

a = findgen(17)*(!pi*2/16.)
clevels = [0.,1.0,2.0,findgen(8)*3.0+5.,40.]
cores = indgen(12)*20+35

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
ncdf_cris = ncdf
;sndex = sort(cris_lat_arr)
;cris_lat_arr = cris_lat_arr[sndex]
;cris_lon_arr = cris_lon_arr[sndex]
;cris_sfc_arr = cris_sfc_arr[sndex]

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

map_set,0.,90.,limit=[min_lat,min_lon,max_lat,max_lon],$
    /cylin,/continent,mlinethick=3,/noborder,title = 'South Asia - April 22, 2015'

axis,xaxis=0,xticks=nlab_lon-1,xtickv=[findgen(nlab_lon)/(nlab_lon-1)],xtickname=string(lon_lab,format='(f6.2)'),xtitle='Longitude'
axis,xaxis=1,xticks=1,xtickname=[' ',' '],xtickv=[min_lon,max_lon]

axis,yaxis=0,yticks=nlab_lat-1,ytickv=[findgen(nlab_lat)/(nlab_lat-1)],ytickname=string(lat_lab,format='(f4.1)'),ytitle='Latitude'
axis,yaxis=1,yticks=1,ytickname=[' ',' '],ytickv=[min_lat,max_lat]

neta = npoly[0]
for ic=0,ncdf_cris-1 do begin
   lat_eta = fltarr(neta)
   lon_eta = fltarr(neta)
   for ie=0,neta-1 do begin
      lat_eta[ie] = lat_bnds[ie,ic]
      lon_eta[ie] = lon_bnds[ie,ic]
  endfor
  cndex = where(clevels gt cris_sfc_arr[ic],ng)
  if (ng gt 0) then id = cndex[0]-1 else id = 11 
  ;if (cris_sfc_arr[ic] gt 15.) then stop
  ;if (cris_sfc_arr[ic] gt 1.0) then begin
  ;if (xa_type[ic] le 2 ) then begin
  ;if (chisq[ic] lt 1000. AND cris_sfc_arr[ic] lt 200.) then begin
  if (cris_sfc_arr[ic] gt 0.0) then begin
     oplot,lon_eta,lat_eta
     polyfill,lon_eta,lat_eta,color=cores[id]
  endif
endfor

draw_scale_bar,'ppbv',[.87,.15,.92,.9],cores[0:11],clevels[0:11],['(f3.0)'],[''],[0.01,0.],charsize=0.7
device,/close_file


end
