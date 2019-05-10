function rd_retr, ncfile=ncfile
;+
;*************************************************************
; Reader for loading retrieval results irsensors/modis
;
; Returns a data the structure with all the fields
; 
; yhe@aer.com, 03/02/2016
;*************************************************************
;-
  if not Keyword_Set(ncfile) then ncfile = 'run/retr.ir.nc'

  ncid = NCDF_OPEN(ncfile)      ; Open The NetCDF file
  profid = NCDF_VARID(ncid,'profiles')
  idid = NCDF_VARID(ncid,'id')
  latid = NCDF_VARID(ncid,'lat')
  lonid = NCDF_VARID(ncid,'lon')
  topoid = NCDF_VARID(ncid,'Landtype')
  plandid = NCDF_VARID(ncid,'pland')
  angleid = NCDF_VARID(ncid,'eia')
  timeid = NCDF_VARID(ncid,'date')
  emrfid = NCDF_VARID(ncid,'EmRf')
  rmsid = NCDF_VARID(ncid,'rms')
  chisqid = NCDF_VARID(ncid,'chisq')
  iterid = NCDF_VARID(ncid,'iteration')
  noiseid = NCDF_VARID(ncid,'noise')
  emmwid = NCDF_VARID(ncid,'EmMw')
  emirid = NCDF_VARID(ncid,'EmIr')

  if profid ge 0 then NCDF_VARGET, ncid, profid, profiles ; Read in variable 'profiles'
  if idid ge 0 then  NCDF_VARGET, ncid, idid, id          ; Read in variable 'id'
  if latid ge 0 then NCDF_VARGET, ncid, latid, lat        ; Read in variable 'lat'
  if lonid ge 0 then NCDF_VARGET, ncid, lonid, lon        ; Read in variable 'lon'
  if topoid ge 0 then NCDF_VARGET, ncid, topoid, Landtype ; Read in variable 'Landtype'
  if plandid ge 0 then NCDF_VARGET, ncid, plandid, pland  ; Read in variable 'pland'
  if angleid ge 0 then NCDF_VARGET, ncid, angleid, eia    ; Read in variable 'eia'
  if timeid ge 0 then NCDF_VARGET, ncid, timeid, date     ; Read in variable 'date'
  if emrfid ge 0 then NCDF_VARGET, ncid, emrfid, EmRf     ; Read in variable 'EmRf'
  if rmsid ge 0 then NCDF_VARGET, ncid, rmsid, rms        ; Read in variable 'rms'
  if chisqid ge 0 then NCDF_VARGET, ncid, chisqid, chisq  ; Read in variable 'chisq'
  if iterid ge 0 then NCDF_VARGET, ncid, iterid, iteration ; Read in variable 'iteration'
  if noiseid ge 0 then NCDF_VARGET, ncid, noiseid, noise   ; Read in variable 'noise'
  if emmwid ge 0 then NCDF_VARGET, ncid, emmwid, EmMw       ; Read in variable 'EmMw'
  if emirid ge 0 then NCDF_VARGET, ncid, emirid, EmIr       ; Read in variable 'EmIr'

  ;; NCDF_ATTGET, ncid, 'SurfacePressure', psfcdim, /GLOBAL
  NCDF_ATTGET,ncid,/GLOBAL,'Temperature',          Itemp
  NCDF_ATTGET,ncid,/GLOBAL,'Tskin',                Itskin
  NCDF_ATTGET,ncid,/GLOBAL,'SurfacePressure',      Ipsfc
  NCDF_ATTGET,ncid,/GLOBAL,'Mol',                  Imol
  NCDF_ATTGET,ncid,/GLOBAL,'MolLen',               Lmol
  NCDF_ATTGET,ncid,/GLOBAL,'LiqCloud',             Iliq
  NCDF_ATTGET,ncid,/GLOBAL,'IceCloud',             Iice
  NCDF_ATTGET,ncid,/GLOBAL,'SurfaceWinds',         Iwind
  shapep=size(profiles,/dimensions)
  if (size(profiles))[0] eq 1 then nprof=1 else nprof=shapep[1]
  vTyp='pressureStd'
  ncinfo=ncdf_inquire(ncid)
  for i=0,ncinfo.ngatts-1 do begin
    attname=ncdf_attname(ncid,i,/GLOBAL)
    if (attname eq 'VertCoordType') then begin
      NCDF_ATTGET,ncid,/GLOBAL,'VertCoordType',      vTypx
      vTyp=string(vTypx)
    endif
  endfor
  if strtrim(vTyp,2) eq 'pressureStd' then begin
     NCDF_ATTGET,ncid,/GLOBAL,'standardPressureGrid', pref
     pres=fltarr(n_elements(pref),nprof)
     for i=0,nprof-1 do pres[*,i]=pref
  endif else begin
    ncdf_varget,ncid, ncdf_varid( ncid,'pressure'),    pres
  endelse

  NCDF_CLOSE, ncid              ; Close the NetCDF file

  if idid lt 0 then id = 'id'
  if timeid lt 0 then date = '0000-00-00'
  if noiseid lt 0 then noise = [0.0]
  if rmsid lt 0 then rms = [0.0]
  if chisqid lt 0 then chisq = [0.0]
  if iterid lt 0 then iteration = [0]
  if angleid lt 0 then eia = [0.0]

  emis = [0.0]
  if emmwid ge 0 then emis = EmMw
  if emirid ge 0 then emis = EmIr
  if emrfid ge 0 then emis = EmRf

  index = {temp: {start: Itemp[0]-1,   len: Itemp[1]}, $
           tskin:{start: Itskin[0]-1,  len: Itskin[1]}, $
           psfc: {start: Ipsfc[0]-1,   len: Ipsfc[1]}, $
           H2O:  {start: Imol[1-1]-1,  len: Lmol[1-1]}, $
           O3:   {start: Imol[3-1]-1,  len: Lmol[3-1]}, $
           clw:  {start: Iliq[0]-1,    len: Iliq[1]}, $
           ciw:  {start: Iice[0]-1,    len: Iice[1]}, $
           wsfc: {start: Iwind[0]-1,   len: Iwind[1]} $
          }

  data = {profs: profiles, $
          index: index, $
          lat: lat, $
          lon: lon, $
          pre: pres, $
          landtype: landType, $
          pland: pland, $
          eia: eia, $
          rms: rms, $
          chisq: chisq, $
          iter: iteration, $
          emis: emis $
         }

  return, data

end
