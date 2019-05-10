pro wr_stats, file, data
;*************************************************************
; IDL script for writing NetCDF file:
;*************************************************************
sz = size(data.forindex)

ncid = NCDF_CREATE(file,/clob)            ; Create the NetCDF file
NCDF_CONTROL, ncid, /FILL   

idnfor = NCDF_DIMDEF(ncid, 'nfor',sz[1])

idfor = NCDF_VARDEF(ncid, "nprofs", idnfor, /SHORT)
NCDF_ATTPUT, ncid, idfor, "long_name", "profile index"
NCDF_ATTPUT, ncid, idfor, "units", "none"
NCDF_ATTPUT, ncid, idfor, "scale_factor", 0.0
NCDF_ATTPUT, ncid, idfor, "add_offset", 0.0
NCDF_ATTPUT, ncid, idfor, "_fillValue", -999
NCDF_ATTPUT, ncid, idfor, "validRange", [0,500]

idtemp = NCDF_VARDEF(ncid,'temp',idnfor,/DOUBLE)
NCDF_ATTPUT, ncid, idtemp, "long_name", "absolute value of relative temperature difference"
NCDF_ATTPUT, ncid, idtemp, "units", "K"
NCDF_ATTPUT, ncid, idtemp, "_fillValue", -999.0
NCDF_ATTPUT, ncid, idtemp, "validRange", [0.0,500.0]

idtskin = NCDF_VARDEF(ncid,'tskin',idnfor,/DOUBLE)
NCDF_ATTPUT, ncid, idtskin, "long_name", "absolute value of relative skin temperature difference"
NCDF_ATTPUT, ncid, idtskin, "units", "K"
NCDF_ATTPUT, ncid, idtskin, "_fillValue", -999.0
NCDF_ATTPUT, ncid, idtskin, "validRange", [0.0,500.0]

idh2o = NCDF_VARDEF(ncid,'h2o',idnfor,/DOUBLE)
NCDF_ATTPUT, ncid, idh2o, "long_name", "absolute value of relative water vapor difference"
NCDF_ATTPUT, ncid, idh2o, "units", "g g-1"
NCDF_ATTPUT, ncid, idh2o, "_fillValue", -999.0
NCDF_ATTPUT, ncid, idh2o, "validRange", [0.0,10.0]

ido3 = NCDF_VARDEF(ncid,'o3',idnfor,/DOUBLE)
NCDF_ATTPUT, ncid, ido3, "long_name", "absolute value of relative ozone difference"
NCDF_ATTPUT, ncid, ido3, "units", "ppt"
NCDF_ATTPUT, ncid, ido3, "_fillValue", -999.0
NCDF_ATTPUT, ncid, ido3, "validRange", [0.0,10.0]

idclw = NCDF_VARDEF(ncid,'clw',idnfor,/DOUBLE)
NCDF_ATTPUT, ncid, idclw, "long_name", "absolute value of relative cloud liquid water difference"
NCDF_ATTPUT, ncid, idclw, "units", "g g-1"
NCDF_ATTPUT, ncid, idclw, "_fillValue", -999.0
NCDF_ATTPUT, ncid, idclw, "validRange", [0.0,10.0]

idrms = NCDF_VARDEF(ncid,'rms',idnfor,/DOUBLE)
NCDF_ATTPUT, ncid, idrms, "long_name", "absolute value of RMS difference"
NCDF_ATTPUT, ncid, idrms, "units", "none"
NCDF_ATTPUT, ncid, idrms, "_fillValue", -999.0
NCDF_ATTPUT, ncid, idrms, "validRange", [0.0,10.0]

idchisq = NCDF_VARDEF(ncid,'chisq',idnfor,/DOUBLE)
NCDF_ATTPUT, ncid, idchisq, "long_name", "absolute value of ChiSq difference"
NCDF_ATTPUT, ncid, idchisq, "units", "none"
NCDF_ATTPUT, ncid, idchisq, "_fillValue", -999.0
NCDF_ATTPUT, ncid, idchisq, "validRange", [0.0,10.0]

iditer = NCDF_VARDEF(ncid,'iter',idnfor,/SHORT)
NCDF_ATTPUT, ncid, iditer, "long_name", "absolute value of maximum iteration difference"
NCDF_ATTPUT, ncid, iditer, "units", "none"
NCDF_ATTPUT, ncid, iditer, "_fillValue", -999.0
NCDF_ATTPUT, ncid, iditer, "validRange", [0.0,10.0]

NCDF_ATTPUT, ncid, /GLOBAL, "CreationDate", systime(/UTC)

NCDF_CONTROL, ncid, /ENDEF ;

NCDF_VARPUT, ncid, idfor, data.forindex
NCDF_VARPUT, ncid, idtemp, data.temp
NCDF_VARPUT, ncid, idtskin, data.tskin
NCDF_VARPUT, ncid, idh2o, data.h2o
NCDF_VARPUT, ncid, ido3, data.o3
NCDF_VARPUT, ncid, idclw, data.clw
NCDF_VARPUT, ncid, idrms, data.rms
NCDF_VARPUT, ncid, idchisq, data.chisq
NCDF_VARPUT, ncid, iditer, data.iter

NCDF_CLOSE, ncid      ; Close the NetCDF file
print, "<<< wr_stats completed >>>"

end
