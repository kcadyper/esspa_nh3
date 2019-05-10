
; Set up plot
plot_type = 'PS'
if (plot_type eq 'PS') then begin
   set_plot,'PS'
   device,/color
   device,/Helvetica,/Bold
   device,filename='plot_cris_gridded.ps'
   !p.font = 0
endif


loadct, 13

!x.thick = 3
!y.thick = 3
!p.thick = 2
!p.charsize = 1.0

!p.multi = [0,1,3]


a = findgen(17)*(!pi*2/16.)
clevels = [findgen(11),15.]
cores = indgen(12)*20+35

!x.style=8
!y.style=8

cmon = 'April 22, 2015'
viirs_dir = '/project/p1913/esspa4_l2/l2_data/'
retr_dir = '/nas/project/p1913/esspa11/viirs_cloud_test/'

min_lat = 10.
max_lat = 65.
min_lon = 255.0 
max_lon = 300.0 

nlab_lat = 12
lat_lab = indgen(nlab_lat)*5.+min_lat
nlab_lon = 10
lon_lab = indgen(nlab_lon)*5.+min_lon

; zoom in
min_lat = 10.
max_lat = 50.
min_lon = 265.0 
max_lon = 295.0 

nlab_lat = 9
lat_lab = indgen(nlab_lat)*5.+min_lat
nlab_lon = 7
lon_lab = indgen(nlab_lon)*5.+min_lon


airmwt = 28.9
nh3mwt = 17.03
mrat = airmwt/nh3mwt

;goto, skip_read


kount = 0
; Read CrIS data
spawn,'ls '+retr_dir+'retr*nc',cdf_files
spawn,'ls '+viirs_dir+'crisimg.snpp.*T18*nc',viirs_files

ncdf = n_elements(cdf_files)
for ic=0,ncdf-1 do begin
   read_ncdf_jh,data,filename=cdf_files[ic]
   if (ic eq 0) then begin
      nh3_id = data.attributes.global.mol[10]
      nh3_len = data.attributes.global.mollen[10]

      cris_lat_arr = data.lat
      cris_lon_arr = data.lon+360.
      cris_sfc_arr = reform(mrat*data.profiles[nh3_id-1+nh3_len-1,*]*1.0e+9)
      xa_type = data.atmosclass
      chisq = data.chisq
      npts_cris = n_elements(cris_lat_arr)
   endif else begin
      cris_lat_arr = [cris_lat_arr,data.lat]
      cris_lon_arr = [cris_lon_arr,data.lon+360.]
      cris_sfc_arr =[cris_sfc_arr, reform(mrat*data.profiles[nh3_id-1+nh3_len-1,*]*1.0e+9)]
      xa_type = [xa_type,data.atmosclass]
      chisq = [chisq,data.chisq]
   endelse

; Read VIIRS  data
   read_ncdf_jh,datav,filename=viirs_files[ic]
   id_band = 3  ; viirs reflectance band M04
   id_viirs_sel=0  ; all pixels
   viirs_subsets = ['All','Clear','Cloudy']
   viirs_refl = datav.viirs_refl

   dim_names = datav.dimensions.names
   id_viirs_sub = where(dim_names eq 'viirs_subset')
   nsubs = datav.dimensions.sizes(id_viirs_sub)
   id_viirs_rband = where(dim_names eq 'viirs_refl_band')
   nrbands = datav.dimensions.sizes(id_viirs_rband)

   if (ic eq 0)  then $
      viirs_refl_arr = reform(viirs_refl[id_band,id_viirs_sel,*,*,*],npts_cris) $ 
   else viirs_refl_arr = [viirs_refl_arr,reform(viirs_refl[id_band,id_viirs_sel,*,*,*],npts_cris) ]

endfor      

      
; Plot maps
p_arr = fltarr(3,4)
p_arr[0,*] = [0.05,0.4,0.32,0.80]
p_arr[1,*] = [0.37,0.4,0.64,0.80]
p_arr[2,*] = [0.69,0.4,0.96,0.80]

; Map NH3, clouds, NH3 without clouds
titulo = ['NH3:pol+mod','NH3:all','NH3:all-clouds']
for ip=0,2 do begin
      
   !p.position = p_arr[ip,*]
   map_set,0.,260.,limit=[min_lat,min_lon,max_lat,max_lon],/advance,$
       /cylin,/continent,/usa,mlinethick=3,/noborder,title = cmon

   axis,xaxis=0,xticks=nlab_lon-1,xtickv=[findgen(nlab_lon)/(nlab_lon-1)],xtickname=string(lon_lab,format='(i3)'),xtitle='Longitude' 
   axis,xaxis=1,xticks=1,xtickname=[' ',' '],xtickv=[min_lon,max_lon]

   axis,yaxis=0,yticks=nlab_lat-1,ytickv=[findgen(nlab_lat)/(nlab_lat-1)],ytickname=string(lat_lab,format='(i2)'),ytitle='Latitude' 
   axis,yaxis=1,yticks=1,ytickname=[' ',' '],ytickv=[min_lat,max_lat]

   ; Average over deltaxdelta deg boxes
   delta = 0.2
   nlat = (max_lat-min_lat)/delta+1
   nlon = (max_lon-min_lon)/delta+1
   lat_arr = findgen(nlat)*delta+min_lat
   lon_arr = findgen(nlon)*delta+min_lon

   for ilat=0,nlat-2 do begin
      andex = where(cris_lat_arr ge lat_arr[ilat] AND cris_lat_arr lt lat_arr[ilat+1],na)
      for ilon=0,nlon-2 do begin 
	 if (na gt 0) then begin
	    ondex = where(cris_lon_arr[andex] ge lon_arr[ilon] AND cris_lon_arr[andex] lt lon_arr[ilon+1],no)
	    if (no gt 0) then begin
	       case ip of
		  0: vndex = where(viirs_refl_arr[andex[ondex]]  gt 0.,nval)
                  1: vndex = where((viirs_refl_arr[andex[ondex]]  gt 0.) $
                             AND (xa_type[andex[ondex]] le 2),nval)
		  2: vndex = where(viirs_refl_arr[andex[ondex]]  lt 0.3,nval)
               endcase
               ;print, ip, ilat,ilon,no,nval
	       if (nval gt 0) then begin
		  mean_nh3 = mean(cris_sfc_arr[andex[ondex[vndex]]])
		  cndex = where(clevels gt mean_nh3,ng)
		  if (ng gt 0) then icolor = cndex[0]-1 else icolor = 11 
		  lat_box = [lat_arr[ilat],lat_arr[ilat],lat_arr[ilat+1],lat_arr[ilat+1]]
		  lon_box = [lon_arr[ilon],lon_arr[ilon+1],lon_arr[ilon+1],lon_arr[ilon]]
		  ;oplot,lon_box,lat_box,color=1
		  polyfill,lon_box,lat_box,color=cores[icolor]
               endif
	    endif
	 endif
      endfor
   endfor

   map_continents
endfor
draw_scale_bar,'ppbv',[.37,.32,.64,.34],cores[0:11],clevels[0:11],['(i2)'],[''],[0.0,-0.02],charsize=0.7
device,/close_file


end
