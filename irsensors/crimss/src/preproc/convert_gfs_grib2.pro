pro  convert_gfs_grib2,file,outpath=outpath

if not keyword_set(file) then file='/home/ipolonsk/sandbox/python/cloud_clearing/gfsanl_4_20150701_0600_000.grb2'
if not keyword_set(outpath) then outpath='./'

print, '========================= RUN ==========================================='
if keyword_set(file) then print, 'file:', file
if keyword_set(outpath) then print, 'outpath:', outpath

excluded_keys = ['gribSection0', 'template not found', 'Parameter information', 'grib 2 Section 5 DATA REPRESENTATION SECTION', $
                 'grib 2 Section 6 BIT-MAP SECTION', 'grib 2 Section 7 data', 'Derived forecast',  'g2grid', 'latLonValues' ,'latitudes', $
                 'longitudes','distinctLatitudes','distinctLongitudes']

n_records = grib_count(file)
i_record = 0
fid = grib_open(file)
parameternames = strarr(n_records)
levels = intarr(n_records)
leveltypes = strarr(n_records)
nT = 0
while (i_record lt n_records) do begin
   h = grib_new_from_file(fid)
   iter = grib_keys_iterator_new(h, /all)
;  Loop over keys in record, looking for the parameter key.
   while grib_keys_iterator_next(iter) do begin
      dum=''
      key = grib_keys_iterator_get_name(iter)
      if total(excluded_keys eq key) ge 1 then continue
      
      if strlowcase(key) eq 'parametername' then begin
        pname = strtrim(grib_get(h, key))
        if pname eq 'Temperature' then print, i_record, ' : ', key, ' : ', pname
      endif
      if strlowcase(key) eq 'level' then begin
        ;print, i_record, ' : ', key, ' : ', grib_get(h, key)
;        if pname eq 'Temperature' then print, i_record, ' : ', key, ' : ', grib_get(h, key)
      endif
      if strlowcase(key) eq 'typeoflevel' then begin
        tname = strtrim(grib_get(h, key))
        if pname eq 'Temperature' and tname eq 'isobaricInhPa' then begin
          nT = nT +1
          print, i_record, ' : ', key, ' : ', grib_get(h, key), nT
        endif
      endif
;      if strlowcase(key) eq 'typeoflevel' then print, i_record, ' : ', key, ' : ', grib_get(h, key)
      
      if strlowcase(key) eq 'parametername' then parameternames[i_record] = strtrim(grib_get(h, key),2)
      if strlowcase(key) eq 'level' then levels[i_record] = grib_get(h, key)
      if strlowcase(key) eq 'typeoflevel' then leveltypes[i_record] = strtrim(grib_get(h, key),2)
   endwhile
   grib_keys_iterator_delete, iter
   grib_release, h
   i_record++
endwhile
grib_close, fid
print,'n_records: ',n_records

;;;  v4, 31 levels 0.5 degree grids
nT0=31L
nRH0=31L
nO30=17L
nx=720L
ny=361L
nprof = nx*ny
Temprecord = where(parameternames eq "Temperature" and leveltypes eq 'isobaricInhPa',nT)
;for i=0,NT-1 do print,parameternames[Temprecord[i]],' ',leveltypes[Temprecord[i]],levels[Temprecord[i]]
RHrecord = where(parameternames eq "Relative humidity" and leveltypes eq 'isobaricInhPa',NRH)
;for i=0,NRH-1 do print,parameternames[RHrecord[i]],' ',leveltypes[RHrecord[i]],levels[RHrecord[i]]
O3record = where(parameternames eq '192' and leveltypes eq 'isobaricInhPa',NO3)

nT0=nT
NRH0=NRH
nO30=NO3

if nT0 ne NT then stop
if nRH0 ne NRH then stop
;for i=0,NO3-1 do print,parameternames[O3record[i]],' ',leveltypes[O3record[i]],levels[O3record[i]]
if nO3 ne NO3 then stop
isoT = double(levels[Temprecord])
isoRH = double(levels[RHrecord])
ind=where(isoT ne isoRH)
if ind[0] ne -1 then stop

nlev = nT0 + nRH0 + nO3  + 2 + 2 + 2 + 1 + 1  
;      Temp + wv + O3 + tsfc + psfc + winds + sfc hgt + pland + tair + rhair
datarecord = make_array(nlev,/long,value=-1)
;; Temp profile
IGTemp0 = 0
for j = 0, nT0-1 do begin
   datarecord[IGtemp0+j]=Temprecord[j]
endfor

;; Tsfc and Psfc
IGTsfc0 = IGTemp0 + nT0
Tsfcrecord = where(parameternames eq "Temperature" and leveltypes eq 'surface' and levels eq 0,cnt)
if cnt ne 1 then stop
datarecord[IGTsfc0] = Tsfcrecord[0]
IGPsfc0 = IGTsfc0 + 1
psfcrecord = where(parameternames eq "Pressure" and leveltypes eq 'surface' and levels eq 0,cnt)
if cnt ne 1 then stop
datarecord[IGPsfc0] = psfcrecord[0]


;; WV profile
IGH2o0 = IGPsfc0 + 1
for j = 0, nRH0-1 do begin
   datarecord[IGH2o0+j]=RHrecord[j]                
endfor

;; O3
IGO30 = IGH2o0 + nRH0
for j = 0, nO30-1 do begin
   datarecord[IGO30+j]=O3record[j]
endfor

;; Wind
IGWind0 = IGO30 + nO30
uwindrecord =where(parameternames eq "u-component of wind" and leveltypes eq 'heightAboveGround' and levels eq 10,cnt)
if cnt ne 1 then stop
datarecord[IGWind0] = uwindrecord[0]
vwindrecord =where(parameternames eq "v-component of wind" and leveltypes eq 'heightAboveGround' and levels eq 10,cnt)
if cnt ne 1 then stop
datarecord[IGWind0+1] = vwindrecord[0]

;; surface hgt
shgtrecord = where(parameternames eq "Geopotential height" and leveltypes eq 'surface',cnt)
if cnt ne 1 then stop
datarecord[IGWind0+2] = shgtrecord[0]

; pland
plandrecord = where(parameternames eq "Land cover" and leveltypes eq 'surface',cnt)
if cnt ne 1 then stop
datarecord[IGWind0+3] = plandrecord[0]

;Tair
tairrecord = where(parameternames eq "Temperature" and leveltypes eq 'heightAboveGround' and levels eq 2,cnt)
if cnt ne 1 then stop
datarecord[IGWind0+4] = tairrecord[0]
; RHair
rhrecord = where(parameternames eq "Relative humidity" and leveltypes eq 'heightAboveGround' and levels eq 2,cnt)
if cnt ne 1 then stop
datarecord[IGWind0+5] = rhrecord[0]

;;  read data file
lbl=fltarr(15)
databuf= dblarr(nlev,nx,ny)
fid = grib_open(file)
h_record = lon64arr(n_records)
for i=0, max(datarecord) do h_record[i] = grib_new_from_file(fid)
for l=0,nlev-1 do begin
    id = datarecord[l]
    h = h_record[id]
    iter = grib_keys_iterator_new(h, /all)
    ;print,l,id,h,' ',iter,'   ',parameternames[id],' ',leveltypes[id],levels[id]
    while grib_keys_iterator_next(iter) do begin
       key = grib_keys_iterator_get_name(iter)
       ;print,key
       if total(excluded_keys eq key) ge 1 then continue
       if strcmp(key, 'values') then begin
          val = grib_get_values(h) ; preserves dimensionality
          databuf[l, *,*]=val
       endif else if l eq 0 then begin
          if strcmp(key, 'Ni') then begin
             lbl[1] = grib_get(h, key)
          endif else if strcmp(key, 'Nj') then begin
             lbl[2] = grib_get(h, key)
          endif else if strcmp(key, 'latitudeOfFirstGridPoint') then begin
             lbl[3] = grib_get(h, key)
          endif else if strcmp(key, 'longitudeOfFirstGridPoint') then begin
             lbl[4] = grib_get(h, key)
          endif else if strcmp(key, 'latitudeOfLastGridPoint') then begin
             lbl[6] = grib_get(h, key)
          endif else if strcmp(key, 'longitudeOfLastGridPoint') then begin
             lbl[7] = grib_get(h, key)
          endif else if strcmp(key, 'iDirectionIncrement') then begin
             lbl[8] = grib_get(h, key)
          endif else if strcmp(key, 'jDirectionIncrement') then begin
             lbl[9] = grib_get(h, key)
          endif else if strcmp(key, 'julianDay') then begin
             lbl[13] = grib_get(h, key)
          endif else if strcmp(key, 'forecastTime') then begin
             lbl[14] = grib_get(h, key)
          endif else if strcmp(key, 'grib2divider') then begin
             grib2divider = grib_get(h, key)
          endif else if strcmp(key, 'gridDefinition') then begin
             gridDefinition = grib_get(h, key)
          endif else if strcmp(key, 'gridDefinitionTemplateNumber') then begin
             lbl[0] = grib_get(h, key)
          endif
      endif
   endwhile
endfor
; Release all the handles and close the file.
for id=0, max(datarecord) do begin
   grib_release, h_record[id]
endfor
grib_close, fid


;  check consistency
if lbl[1] ne nx then stop
if lbl[2] ne ny then stop

if lbl[3] ne 9e7 then stop  ; North 90 deg
if lbl[7] ne 3.595e8  then stop  ; East 395.5 deg 
if lbl[8] ne 500000 then stop  ; 0.5 deg increment
if lbl[9] ne 500000  then stop  ; 0.5 deg increment

;  convert RH to mixing ratio
for i=0,nx-1 do for j=0,ny-1 do begin
   for l=IGh2o0,IGh2o0+nRH0-1 do begin
      RH =  databuf[l,i,j]
      RHT = databuf[l-nT0-2,i,j]
      a = 0.01 * RH *svp(RHT)
      
      press=isoRH[l-nT0-2]
      Y = 0.622*a/(press-a)     ;in g/g
      Y =max([Y,2e-6])
      databuf[l,i,j]=y
   endfor

   ; RH Air
   l = nlev-1
   RH =  databuf[l,i,j]
   RHT = databuf[l-1,i,j]
   a = 0.01 * RH *svp(RHT)
   
   press=databuf[nT0+1,i,j]*0.01
   Y = 0.622*a/(press-a)        ;in g/g
   Y =max([Y,2e-6])
   databuf[l,i,j]=y
endfor
print,max(databuf[nT0+1,*,*],min=pp),pp
databuf[nT0+1,*,*]=databuf[nT0+1,*,*]*0.01  ; convert to mbar

vCoordTyp=string([112B,114B,101B,115B,115B,117B,114B,101B,83B,116B,100B,32B])
vCoordTyp='pressureStd'


profread = fltarr(nlev,nx,ny)
;  reverse lat grid, start from -90 deg
for l=0,nlev-1 do for i=0,nx-1 do profread[l,i,*]=reverse(databuf[l,i,*],3)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;     finish reading ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 ; add the bottom layer at 1100mb
nTout = nT0+1
nRHout = nRH0+1
nO3out = nTout
nlevout = nTout+nRHout+nO3out+4  ; Temp, Tsfc, Psfc, WV, O3, wind

;; fill in the output arrays
profout = fltarr(nlevout,nx,ny)
profout[0:nT0-1,*,*] = profread[0:nT0-1,*,*]   ; Temp
profout[nTout:nTout+1,*,*] = profread[nT0:nT0+1,*,*]  ; Tsfc + Psfc
profout[nTout+2:nTout+nRH0+1,*,*] = profread[nT0+2:nT0+nRH0+1,*,*]  ; WV
profout[nTout+NRHout+2:nTout+nRHout+NO30+1,*,*] = profread[nT0+nRH0+2:nT0+nRH0+NO30+1,*,*] ; O3
profout[nTout+NRHout+nO3out+2:nTout+nRHout+nO3out+3,*,*] = profread[nT0+nRH0+nO30+2:nT0+nRH0+nO30+3,*,*]  ; wind

for i=0,nx-1 do for j=0,ny-1 do begin
   profout[nTout+nRHout+NO30+2:nTout+nRHout+NO3out+1,i,j] = profout[nTout+nRHout+NO30+1,i,j]/2
endfor
pbot=1100.
deltaP=1
p0 = double(pbot)
for i=0,nx-1 do for j=0,ny-1 do begin
   ;print,i,j

   ind=where(isoT lt profread[nT0+1,i,j],cnt)   ;psfc

   if cnt eq nT0 then begin      
      p1 = double(isoT[nT0-1])
      p2 = double(profread[nT0+1,i,j]) ; psfc
      x1 = double(profread[nT0-1,i,j])
      x2 = double(profread[nlev-2,i,j])      ; Tair
   endif else begin
      p1 = double(isoT[nT0-2])
      p2 = double(isoT[nT0-1])
      x1 = double(profread[nT0-2,i,j])
      x2 = double(profread[nT0-1,i,j])
   endelse
   if p2-p1 ge deltaP then x0=lint_log(x1,x2,p1,p2,p0) else x0=x2
   ;print,x0,x2,reform(profread[cnt-1:nT0-1,i,j])
   ;print,p0,p2,reform(isoT[cnt-1:nT0-1])
   profout[nTout-1,i,j]=x0

   if cnt eq nT0 then begin
      p1 = double(isoT[nT0-1])
      p2 = double(profread[nT0+1,i,j])    ; psfc
      x1 = double(profread[nT0+nRH0+1,i,j])
      x2 = double(profread[nlev-1,i,j])           ; RHair
   endif else begin
      p1 = double(isoT[nT0-2])
      p2 = double(isoT[nT0-1])
      x1 = double(profread[nT0+nRH0,i,j])
      x2 = double(profread[nT0+nRH0+1,i,j])
   endelse      
   if (p2-p1) ge deltaP then x0=lint_log(x1,x2,p1,p2,p0) else x0=x2
   ;print,x0,x1,x2,reform(profread[cnt-1+nRH0+2:nT0+nRH0+1,i,j])
   ;print,p0,p1,p2,reform(isoT[cnt-1:nT0-1])
   profout[nTout+nRHout+1,i,j]=x0
  

   continue   
   if cnt eq nT0 then begin
      p1 = double(isoT[nT0-1])
      p2 = double(profread[nT0+1,i,j])    ; psfc
      x1 = double(profread[nT0+nRH0+1,i,j])
      x2 = double(profread[nlev-1,i,j])           ; RHair
   endif else begin
      p1 = double(isoT[nT0-2])
      p2 = double(isoT[nT0-1])
      x1 = double(profread[nT0+nRH0,i,j])
      x2 = double(profread[nT0+nRH0+1,i,j])
   endelse
   x1 = x1/(1+x1)                   ;; mixting ratio to mass fraction
   x2 = x2/(1+x2)      
   if (p2-p1) ge deltaP then x0=lint_log(x1,x2,p1,p2,p0) else x0=x2
   print,x0,x1,x2,reform(profread[cnt-1+nRH0+2:nT0+nRH0+1,i,j]),x0/(1-x0)
   print,p0,p1,p2,reform(isoT[cnt-1:nT0-1])
   x0 = x0/(1-x0)
   profout[nT0+nRH0+3,i,j]=x0  
endfor

qc=intarr(2,nprof)
ind=where((profout[nTout+nRHout+1,*,*] le 0 or profout[nTout+nRHout+1,*,*] gt 0.1) and profread[nTout+1,*,*] le 1000,cnt1 )
if ind[0] ne -1 then begin
   tp = reform(profout[nTout+nRHout+1,*,*])   
   print,'WV1',cnt1,max(tp[ind]),min(tp[ind])
   qc[0,ind]=1
endif
ind=where((profout[nTout+nRHout+1,*,*] le 0 or profout[nTout+nRHout+1,*,*] gt 0.1) and profread[nTout+1,*,*] ge 1000,cnt1 )
if ind[0] ne -1 then begin
   tp = reform(profout[nTout+nRHout+1,*,*])   
   print,'WV2',cnt1,max(tp[ind]),min(tp[ind])
   qc[0,ind]=1
endif
ind=where((profout[nTout-1,*,*] lt 200 or profout[nTout-1,*,*] gt 400)  and profread[nTout+1,*,*] le 1000,cnt1)
if ind[0] ne -1 then begin
   tp = reform(profout[nTout-1,*,*])   
   print,'T1',cnt1,max(tp[ind]),min(tp[ind])
   qc[0,ind]=1
endif
ind=where((profout[nTout-1,*,*] lt 200 or profout[nTout-1,*,*] gt 400)  and profread[nTout+1,*,*] gt 1000,cnt1)
if ind[0] ne -1 then begin
   tp = reform(profout[nTout-1,*,*])   
   print,'T2',cnt1,max(tp[ind]),min(tp[ind])
   qc[0,ind]=1
endif
ind=where(qc[0,*] eq 1,cnt1)
print,'qc',cnt1

;qc[*]=0
;  output grids
res=0.5
lat = indgen(ny)*res-90
lon = indgen(nx)*res
latarr = fltarr(nx,ny)
lonarr = fltarr(nx,ny)
for i=0,ny-1 do lonarr[*,i]=lon
for i=0,nx-1 do latarr[i,*]=lat
caldat, lbl[13],month,day,year,hh
    
jday = julday(month,day,year,0,0,0) - julday(1,1,year,0,0,0)+1
datestr = string(year,'(I4.4)')+string(month,'(I2.2)')+string(day,'(I2.2)')+string(hh,'(I2.2)')
;datestr = string(year,'(I4.4)')+string(jday,'(I3.3)')  ; jday
outputfile = outpath+'/gfsanl_4_'+datestr[0]+'.nc'

ncid=ncdf_create(outputfile,/clobber)
recId = NCDF_DIMDEF(ncid, "nProfiles", /unlimited)
parID= NCDF_DIMDEF(ncid, "nPar", nlevout)
nQcid = NCDF_DIMDEF(ncid, "nQC",  2)

profilesId = NCDF_VARDEF(ncid, 'profiles', [parId, recId], /FLOAT)
NCDF_ATTPUT,ncid,profilesId,"long_name","set of profile vectors"
NCDF_ATTPUT,ncid,profilesId,"units","Temp(Kelvin), H2O(g/g) in mol(1), psfc(mb)" 

latId = NCDF_VARDEF(ncid, "lat", [recId], /FLOAT)
NCDF_ATTPUT,ncid,latId,"long_name","latitude"
NCDF_ATTPUT,ncid,latId,"Units","degrees north"

lonId = NCDF_VARDEF(ncid, "lon", [recId], /FLOAT)
NCDF_ATTPUT,ncid,lonId,"long_name","longitude"
NCDF_ATTPUT,ncid,lonId,"Units","degrees east"

surfaltID=NCDF_VARDEF(ncid, "surfalt", [recId], /FLOAT)
NCDF_ATTPUT,ncid,surfaltID,"long_name","Surface Altitude"
NCDF_ATTPUT,ncid,surfaltID,"Units","m"

plandId = NCDF_VARDEF(ncid, "pland", [recId], /FLOAT)
NCDF_ATTPUT,ncid,plandId,"long_name","% land mask"

qcId = NCDF_VARDEF(ncid, "qc", [nQcId,recId], /long)
NCDF_ATTPUT,ncid,qcId,"long_name","Quality control byte array"

ncdf_attput,ncid,/GLOBAL, "Creationtime", systime()
ncdf_attput,ncid,/GLOBAL, "Temperature", [long(1),long(nTout)]
ncdf_attput,ncid,/GLOBAL, "Tskin", [long(nTout+1),1]
ncdf_attput,ncid,/GLOBAL, "SurfacePressure",[long(nTout+2),1]
ncdf_attput,ncid,/GLOBAL, "Mol", [long(nTout+3),long(nTout+nRHout+3),long(nTout+nRHout+3),lonarr(32)+nTout+nRHout+nO3out+3]
ncdf_attput,ncid,/GLOBAL, "MolLen",[long(nRHout),0,long(nO3out),lonarr(32)]
ncdf_attput,ncid,/GLOBAL, "LiqCloud",[long(nTout+nRHout+nO3out+3),long(0)]
ncdf_attput,ncid,/GLOBAL, "IceCloud",[long(nTout+nRHout+nO3out+3),long(0)]
ncdf_attput,ncid,/GLOBAL, "SurfaceWinds",[long(nTout+nRHout+nO3out+3),long(2)]
ncdf_attput,ncid,/GLOBAL, "nlev",long(nTout)
;ncdf_attput,ncid,/GLOBAL, "VertCoordType",vCoordTyp
ncdf_attput,ncid,/GLOBAL, "standardPressureGrid",[float(isoT),pbot]
ncdf_attput,ncid,/GLOBAL, "Source",file_basename(file)
ncdf_control,ncid,/endef


ncdf_varput,ncid,ncdf_varid(ncid,'profiles'),reform(profout,nlevout,nprof)
ncdf_varput,ncid,ncdf_varid(ncid,'lat'),reform(latarr,nprof)
ncdf_varput,ncid,ncdf_varid(ncid,'lon'),reform(lonarr,nprof)
ncdf_varput,ncid,ncdf_varid(ncid,'surfalt'),reform(profread[IGWind0+2,*,*],nprof)
ncdf_varput,ncid,ncdf_varid(ncid,'pland'),reform(profread[IGWind0+3,*,*],nprof)
ncdf_varput,ncid,ncdf_varid(ncid,'qc'),qc*0
ncdf_close,ncid
print, '========================= DONE ==========================================='

end
