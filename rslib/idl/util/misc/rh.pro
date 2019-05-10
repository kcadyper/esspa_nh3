;-------------------------------------------------------------------
; $Id$
;-------------------------------------------------------------------

;+
; Function: relative_humidity
; Purpose: Calculates relative humidity
; Usage: s = relative_humidity(T,p,q)
; 
; Input: 
;   T = temperature (K)
;   p = pressure (mb)
;   q = specific humidity (g/g)
;
; History:
; 
;   May 16, 2005: Extracted from fdaat_MeteorFct.pro, which was
;                 written in 2004 by S.A. Boukabara to be compatible
;                 with svp, and modified to test the comfortablility
;                 of input variables .................... Y. He, AER
;
; Copyright, 2003-2005, AER. Inc., All rights reserved.
;-

function relative_humidity,T,p,q
a  = svp(t)
qs = 0.622*a/(p-a)
return, 100.*q/qs 
end

;+
; Function: rh
; Purpose: Calculates relative humidity
; Usage: s = rh(T,p,q)
; 
; Input: 
;   T = temperature or temperature profile (K)
;   p = pressure or pressure levels (mb)
;   q = specific humidity or specific humidity array (g/g)
;
; History:
; 
;   May 16, 2005: Extracted from fdaat_MeteorFct.pro, which was
;                 written in 2004 by S.A. Boukabara to be compatible
;                 with svp, and modified to test the comfortablility
;                 of input variables .................... Y. He, AER
;
; Copyright, 2003-2005, AER. Inc., All rights reserved.
;-

FUNCTION rh, T, p, q
if n_elements(T) ne n_elements(p) or n_elements(T) ne n_elements(q) or $
  n_elements(p) ne n_elements(q) then $
  message, "err[rh]: variable sizes not comfortable."

if n_elements(T) gt 1 then begin
    n=n_elements(T)
    RH=fltarr(n)
    for i=0L,n-1 do begin
        RH[i] = relative_humidity(T(i),p(i),q(i))
    endfor
endif else begin
    RH = relative_humidity(T,p,q)
endelse
Return, RH
END
