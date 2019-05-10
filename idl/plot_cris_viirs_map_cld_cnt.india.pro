; Plot VIIRS derived fields 
; This version reads in coordinates of pixel edges from L2 file

; Set up plot
plot_type = 'PS'
if (plot_type eq 'PS') then begin
   set_plot,'PS'
   device,/color
   device,/Helvetica,/Bold
   !p.font = 0
endif


loadct, 13

!x.thick = 1
!y.thick = 1
!p.thick = 1
!p.charsize = 1.0

a = findgen(17)*(!pi*2/16.)
clevels = findgen(11)*100
cores = indgen(11)*22+35


; Read VIIRS  data
read_ncdf_jh,datav,filename='/project/p1913/esspa4_l2/l2_data/crisimg.snpp.20150422T0800.nc'
id_band = 3  ; viirs reflectance band M04
viirs_subsets = ['All','Clear','Cloudy']
cris_lat_arr = datav.lat
cris_lon_arr = datav.lon
ncdf_cris = n_elements(cris_lat_arr)
viirs_count = datav.viirs_count

dim_names = datav.dimensions.names
id_viirs_sub = where(dim_names eq 'viirs_subset')
nsubs = datav.dimensions.sizes(id_viirs_sub)

viirs_count = reform(viirs_count,nsubs,ncdf_cris)

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

; Loop over pixel sets (all, clear, cloudy)
for is=0,2 do begin
   device,filename='plot_cris_viirs_map_cld_count.'+viirs_subsets[is]+'.ps'
   map_set,0.,90.,limit=[min_lat,min_lon,max_lat,max_lon],$
       /cylin,/continent,mlinethick=3,/noborder,title = 'South Asia - April 22, 2015 '+viirs_subsets[is]

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
     cndex = where(clevels gt viirs_count[is,ic],ng)
     if (ng gt 0) then id = cndex[0]-1 else id = 11 
     if (viirs_count[is,ic] ge 0.0) then begin
	oplot,lon_eta,lat_eta
	polyfill,lon_eta,lat_eta,color=cores[id]
     endif
   endfor
   draw_scale_bar,'count',[.87,.15,.92,.9],cores[0:10],clevels[0:10],['(i4)'],[''],[0.01,0.],charsize=0.7
   device,/close_file
endfor



end
