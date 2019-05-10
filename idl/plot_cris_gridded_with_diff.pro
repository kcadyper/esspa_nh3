; !!!!! not complete !!!!
; Set up plot
plot_type = 'PS'
if (plot_type eq 'PS') then begin
   set_plot,'PS'
   device,/color
   device,/Helvetica,/Bold
   device,filename='plot_cris_gridded_diff.ps'
   !p.font = 0
endif


loadct, 13

!x.thick = 3
!y.thick = 3
!p.thick = 2
!p.charsize = 0.7

!p.multi = [0,2,2]
p_arr = fltarr(4,4)
p_arr[0,*] = [0.05,0.55,0.45,0.95]
p_arr[1,*] = [0.55,0.55,0.95,0.95]
p_arr[2,*] = [0.05,0.10,0.45,0.5]
p_arr[3,*] = [0.55,0.10,0.95,0.5]

b_arr = fltarr(4,4)
b_arr[0,*] = [0.45,0.55,0.47,0.95]
b_arr[1,*] = [0.95,0.55,0.97,0.95]
b_arr[2,*] = [0.45,0.10,0.47,0.5]
b_arr[3,*] = [0.95,0.10,0.97,0.5]


a = findgen(17)*(!pi*2/16.)
clevels = [findgen(11),15.]
cores = indgen(12)*20+35

dlevels = findgen(12)-1.0

!x.style=8
!y.style=8

cmon = 'April 22, 2015'
viirs_dir = '/project/p1913/esspa4_l2/l2_data/'
retr_dir = '/nas/project/p1913/esspa11/viirs_cloud_test/'

min_lat = 0.0 
max_lat = 50.0
min_lon = 55.0
max_lon = 95.0

nlab_lat = 6 
lat_lab = indgen(nlab_lat)*10.0+min_lat
nlab_lon = 5 
lon_lab = indgen(nlab_lon)*10.0+min_lon


airmwt = 28.9
nh3mwt = 17.03
mrat = airmwt/nh3mwt

;goto, skip_read


kount = 0
; Read CrIS data
spawn,'ls '+retr_dir+'retr*india*nc',cdf_files
spawn,'ls '+viirs_dir+'crisimg.snpp.*T0*nc',viirs_files

ncdf = n_elements(cdf_files)
for ic=0,ncdf-1 do begin
   read_ncdf_jh,data,filename=cdf_files[ic]
   if (ic eq 0) then begin
      nh3_id = data.attributes.global.mol[10]
      nh3_len = data.attributes.global.mollen[10]

      cris_lat_arr = data.lat
      cris_lon_arr = data.lon
      cris_sfc_arr = reform(mrat*data.profiles[nh3_id-1+nh3_len-1,*]*1.0e+9)
      xa_type = data.atmosclass
      chisq = data.chisq
      npts_cris = n_elements(cris_lat_arr)
   endif else begin
      cris_lat_arr = [cris_lat_arr,data.lat]
      cris_lon_arr = [cris_lon_arr,data.lon]
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

; Map NH3, Nh3 no unp, NH3 all without clouds
delta = 1.0
nlat = (max_lat-min_lat)/delta+1
nlon = (max_lon-min_lon)/delta+1
lat_arr = findgen(nlat)*delta+min_lat
lon_arr = findgen(nlon)*delta+min_lon
mean_nh3_arr_cris = fltarr(3,nlat,nlon)
mean_nh3_arr_cris[*,*,*] = -999.

titulo = ['NH3:pol+mod','NH3:all','NH3:all-clouds']
for ip=0,3 do begin
  !p.position = p_arr[ip,*]
      
   map_set,0.,75.,limit=[min_lat,min_lon,max_lat,max_lon],/advance,$
       /cylin,/continent,mlinethick=3,/noborder

   axis,xaxis=0,xticks=nlab_lon-1,xtickv=[findgen(nlab_lon)/(nlab_lon-1)],xtickname=[string(lon_lab[0:nlab_lon-2],format='(i3)'),' ']
   axis,xaxis=1,xticks=1,xtickname=[' ',' '],xtickv=[min_lon,max_lon]

   axis,yaxis=0,yticks=nlab_lat-1,ytickv=[findgen(nlab_lat)/(nlab_lat-1)],ytickname=string(lat_lab,format='(i2)')
   axis,yaxis=1,yticks=1,ytickname=[' ',' '],ytickv=[min_lat,max_lat]

   if (ip le 2) then begin
   ; Average over deltaxdelta deg boxes

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
		     mean_nh3_arr_cris[ip,ilat,ilon] = mean_nh3
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
      draw_scale_bar,'ppb',b_arr[ip,*],cores[0:11],clevels[0:11],['(i2)'],[''],[0.01,-0.0],charsize=0.5
   endif else begin
      for ilat=0,nlat-2 do begin
	 for ilon=0,nlon-2 do begin 
            if ((mean_nh3_arr_cris[1,ilat,ilon] gt 0.) AND (mean_nh3_arr_cris[2,ilat,ilon] gt 0.)) then begin
               diff = mean_nh3_arr_cris[1,ilat,ilon]-mean_nh3_arr_cris[2,ilat,ilon] 
               print, ilat,ilon,mean_nh3_arr_cris[1,ilat,ilon],mean_nh3_arr_cris[2,ilat,ilon],diff
               cndex = where(dlevels gt diff)
	       if (ng gt 0) then icolor = cndex[0]-1 else icolor = 11 
	       lat_box = [lat_arr[ilat],lat_arr[ilat],lat_arr[ilat+1],lat_arr[ilat+1]]
	       lon_box = [lon_arr[ilon],lon_arr[ilon+1],lon_arr[ilon+1],lon_arr[ilon]]
	       ;oplot,lon_box,lat_box,color=1
	       polyfill,lon_box,lat_box,color=cores[icolor]
	    endif
         endfor
      endfor
      map_continents
      draw_scale_bar,'ppb',b_arr[ip,*],cores[0:11],dlevels[0:11],['(i2)'],[''],[0.01,-0.0],charsize=0.5
   endelse

endfor
device,/close_file


end
