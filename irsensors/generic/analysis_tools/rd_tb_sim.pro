;*************************************************************
; IDL script for reading NetCDF file:
;*************************************************************
function rd_tb_sim, ncfile=ncfile
;; print, 'Warning:  If you have moved rad_sim_rain_melt.20160128.01.nc from ../data/output/'+string(10B)+$
;; 	string(9B)+'idl will not be able to open the file unless you modify'+$
;; 	string(10B)+string(9B)+'the NCDF_OPEN line in this script to reflect the new path.'

if not Keyword_Set(ncfile) then ncfile = '../data/output/rad_sim_rain_melt.20160128.01.nc'

ncid = NCDF_OPEN(ncfile)            ; Open The NetCDF file

NCDF_VARGET, ncid,  0, lat      ; Read in variable 'lat'

NCDF_VARGET, ncid,  1, lon      ; Read in variable 'lon'

NCDF_VARGET, ncid,  2, eia      ; Read in variable 'eia'

NCDF_VARGET, ncid,  3, qc      ; Read in variable 'qc'

NCDF_VARGET, ncid,  4, tb      ; Read in variable 'tb'

NCDF_CLOSE, ncid      ; Close the NetCDF file

return, tb

end
