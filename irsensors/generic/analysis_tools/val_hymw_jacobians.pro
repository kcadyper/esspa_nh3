;; base = '/project/p1855/alipton/mwsensors/hymw/run/wrf_4cases_oc_Jacobians_v4_REDO.nc'
;; base = '/nas/home/yhe/projects/OSScloudInterface/mwsensors_cloud/hymw/run/wrf_4cases_oc_Jacobians_v4_REDO.nc'
;; base = '/nas/home/yhe/projects/OSScloudInterface/hymw_test/hymw/run/wrf_4cases_oc_Jacobians_v4_REDO.baseline
;; .nc'
;; test = '/nas/home/yhe/projects/OSScloudInterface/mwsensors_cloud/hymw/run/wrf_4cases_oc_Jacobians_v4_REDO.nc'
base = 'run/wrf_4cases_oc_Jacobians_v4_REDO.TEST.nc'
test = 'run/wrf_4cases_oc_Jacobians_v4_REDO.nc'

bData = rd_Jacobians(ncfile=base)
tData = rd_Jacobians(ncfile=test)

ntagBase = n_tags(bData)
tagsBase = tag_names(bData)

ntagTest = n_tags(tData)
tagsTest = tag_names(tData)

if ntagBase ne ntagTest then begin
   print, "[err]: mismatch number of tags", ntagBase, ' != ', ntagTest
   stop
endif

for itg = 0, ntagTest-1 do begin
   if tagsBase[itg] ne tagsTest[itg] then begin
      print, "[err]: mismatch tag - ", tagsBase[itg], ' != ', tagsTest[itg]
      stop
   endif
endfor

for itg = 0, ntagTest-1 do begin
   delta = tData.(itg) - bData.(itg)
   print, itg, ',', tagsTest[itg], ',  min/max(test-base), ', min(delta), ',', max(delta)
endfor

;; for itg = 0, ntagTest-1 do begin
;;    delta = histo(tData.(itg) - bData.(itg))
;;    print, itg, ' ', tagsTest[itg]
;;    plot, delta
;;    pause
;; endfor

;; for ich = 0, (size(tdata.kicmr,/dim))[1]-1 do begin
;;    for ilev = 0, (size(tdata.kicmr,/dim))[0]-1 do begin 
;;       delta = tdata.kicmr[ilev,ich] - bdata.kicmr[ilev,ich]
;;       if (delta ne 0.0) then begin
;;          print, ich, ilev, delta
;;       endif
;;    endfor
;;    pause
;; endfor
end
