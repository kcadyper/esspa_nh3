function rd_Jacobians, ncfile=ncfile
;+
;*************************************************************
; Reader for loading Jacobians generated by mwsensors/hymw
;
; Returns a data the structure with all the fields
; 
; yhe@aer.com, 02/02/2016
;*************************************************************
;-
if not Keyword_Set(ncfile) then ncfile = 'run/wrf_4cases_oc_Jacobians_v4_REDO.nc'

ncid = NCDF_OPEN(ncfile)            ; Open The NetCDF file

NCDF_VARGET, ncid,  0, date      ; Read in variable 'date'

NCDF_VARGET, ncid,  1, id      ; Read in variable 'id'

NCDF_VARGET, ncid,  2, lat      ; Read in variable 'lat'

NCDF_VARGET, ncid,  3, lon      ; Read in variable 'lon'

NCDF_VARGET, ncid,  4, eia      ; Read in variable 'eia'

NCDF_VARGET, ncid,  5, qc      ; Read in variable 'qc'

NCDF_VARGET, ncid,  6, rad      ; Read in variable 'rad'

NCDF_VARGET, ncid,  7, Ktemp      ; Read in variable 'Ktemp'

NCDF_VARGET, ncid,  8, Kh2o      ; Read in variable 'Kh2o'

NCDF_VARGET, ncid,  9, Ktskin      ; Read in variable 'Ktskin'

NCDF_VARGET, ncid,  10, Kemis      ; Read in variable 'Kemis'

NCDF_VARGET, ncid,  11, KlqMR      ; Read in variable 'KlqMR'

NCDF_VARGET, ncid,  12, KicMR      ; Read in variable 'KicMR'

NCDF_VARGET, ncid,  13, KrnMR      ; Read in variable 'KrnMR'

NCDF_VARGET, ncid,  14, KsnMR      ; Read in variable 'KsnMR'

NCDF_VARGET, ncid,  15, KgpMR      ; Read in variable 'KgpMR'

NCDF_VARGET, ncid,  16, KlqReff      ; Read in variable 'KlqReff'

NCDF_VARGET, ncid,  17, KicReff      ; Read in variable 'KicReff'

NCDF_VARGET, ncid,  18, KrnReff      ; Read in variable 'KrnReff'

NCDF_VARGET, ncid,  19, KsnReff      ; Read in variable 'KsnReff'

NCDF_VARGET, ncid,  20, KgpReff      ; Read in variable 'KgpReff'

NCDF_CLOSE, ncid      ; Close the NetCDF file

data = {Jacobians, $
        date: date, $
        id: id, $
        lat: lat, $
        lon: lon, $
        eia: eia, $
        qc: qc, $
        rad: rad, $
        Ktemp: Ktemp, $
        Kh2o: Kh2o, $
        Ktskin: Ktskin, $
        Kemis: Kemis, $
        KlqMR: KlqMR, $
        KicMR: KicMR, $
        KrnMR: KrnMR, $
        KsnMR: KsnMR, $
        KgpMR: KgpMR, $
        KlqReff: KlqReff, $
        KicReff: KicReff, $
        KrnReff: KrnReff, $
        KsnReff: KsnReff, $
        KgpReff: KgpReff $
       }

return, data

end
