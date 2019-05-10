baselineFile = '../data/baseline/rad_sim_rain_melt.nc'
;; testFile = '../data/output/rad_sim_rain_melt.20160129.11.nc'
;; baselineFile = '../data/baseline/rad_sim_all5.nc'
;; testFile = '../data/output/rad_sim_all5.20160129.12.nc'
;; testFile = '../mwsensors/amsu/run/rad_sim_all5.nc'
testFile = '../mwsensors/amsu/run/rad_sim_rain+melt.nc'

;; baselineFile = '../oss_cloud_struct_benchmark/rad_sim_all5.nc'
;; testFile = 'rad_sim_rain_melt.20160129.01.nc'

tb = rd_tb_sim(ncfile=baselineFile)
tb_r = rd_tb_sim(ncfile=testFile)

cnt_diff = 0
cnt_nan = 0
for i = 0, n_elements(tb)-1 do begin
   if (abs(tb[i]-tb_r[i]) gt 1.e-5) then begin
      cnt_diff = cnt_diff + 1
   endif

   if (finite(tb_r[i],/nan)) then begin
      cnt_nan = cnt_nan + 1
   endif
endfor
print, "#diff pts = ", cnt_diff, 100*float(cnt_diff)/float(n_elements(tb)),"%"
print, "#NaN = ", cnt_nan, 100*float(cnt_nan)/float(n_elements(tb)),"%"

end
