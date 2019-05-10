function defNcVar, cdfID, name, dim, dtype
  flt=1
  case dtype of
     1: byt=1
     2: shr=1
     3: lng=1
     4: flt=1
     5: dbl=1
     7: str=1
     12: ush=1
     13: ulg=1
     15: u64=1
     else: print, "[msg]... using default data type float..."
  endcase
  vid = NCDF_VARDEF(cdfID,name,dim, $
                    BYTE=byt, $
                    CHAR=chr, $
                    DOUBLE=dbl, $
                    FLOAT=flt, $
                    LONG=lng, $
                    SHORT=shr, $
                    STRING=str, $
                    UBYTE=uby, $
                    UINT64=u64, $
                    ULONG=ulg, $
                    USHORT=ush)
  return, vid
end

function defFill, dtype
  _fill = -999.
  case dtype of
     1: _fill=255
     2: shr=-999
     3: lng=-9999
     4: flt=-999.
     5: dbl=-999.D0
     7: str='X'
     12: ush=65536
     13: ulg=999999L
     15: u64=999999L
     else: print, "[msg]... using default fillValue for float..."
  endcase
  return, _fill
end

pro wr_xy, file, x, y
;*************************************************************
; create a netCDF file for a quick plot
;*************************************************************
if n_elements(x) ne n_elements(y) then stop, "[err]: ...x and y in different size..."
nx = n_elements(x)
ny = n_elements(y)

ncid = NCDF_CREATE(file,/clob)            ; Create the NetCDF file
NCDF_CONTROL, ncid, /FILL   

idnx = NCDF_DIMDEF(ncid, 'nx', nx)
idny = NCDF_DIMDEF(ncid, 'ny', ny)

;; idx = NCDF_VARDEF(ncid,'x',idnx,/DOUBLE)
idx = defNcVar(ncid,'x', idnx, size(x,/type))
NCDF_ATTPUT, ncid, idx, "long_name", "x"
NCDF_ATTPUT, ncid, idx, "units", "none"
NCDF_ATTPUT, ncid, idx, "_fillValue", defFill(size(x,/type))
NCDF_ATTPUT, ncid, idx, "validRange", [min(x),max(x)]

;; idy = NCDF_VARDEF(ncid,'y',idny,/DOUBLE)
idy = defNcVar(ncid,'y', idny, size(y,/type))
NCDF_ATTPUT, ncid, idy, "long_name", "y"
NCDF_ATTPUT, ncid, idy, "units", "none"
NCDF_ATTPUT, ncid, idy, "_fillValue", defFill(size(y,/type))
NCDF_ATTPUT, ncid, idy, "validRange", [min(x),max(x)]

NCDF_ATTPUT, ncid, /GLOBAL, "CreationDate", systime(/UTC)

NCDF_CONTROL, ncid, /ENDEF ;

NCDF_VARPUT, ncid, idx, x
NCDF_VARPUT, ncid, idy, y

NCDF_CLOSE, ncid      ; Close the NetCDF file

print, "<<< wr_xy completed >>>"

end
