;-------------------------------------------------------------------
; $Id$
;-------------------------------------------------------------------

;+
; Function: precipwater
; Purpose: Calculates precipitable water (kg/m^2)
; Usage: precipWater,p,psfc,h2o,qint
; 
; Inputs: 
;    p = pressure profile (mb)
; psfc = surface pressure (mb)
;  h2o = water vapor profile (g/g)
;
; Output:
; qint = integrated precipitable water (kg/m^2)
;
; History:
;   May 16, 2005: Extracted from fdaat_MeteorFct.pro ... Y. He, AER
;
; Copyright, 2003-2005, AER. Inc., All rights reserved.
;-

pro precipWater,p,psfc,h2o,qint

ind_=where(p ge psfc,icnt)
if (icnt ge 1) then nsfc=min(ind_) else nsfc=n_elements(p)-1
QINT=0.
;--case where P(NSFC)=PSFC
P0=PSFC
Q0=h2o(nsfc)
;--case where P(NSFC)<>PSFC
if(psfc ne p(nsfc))then begin
    ;----(log q / log P)= cte
    x=alog(h2o(nsfc-1)/h2o(nsfc))/alog(p(nsfc-1)/p(nsfc))
    q0=h2o(nsfc)*(psfc/p(nsfc))^x
                                ;----( q /  P)= cte
                                ;pint=p(nsfc-1)-p(nsfc)
                                ;a=(p(nsfc-1)-psfc)/pint
                                ;b=(psfc-p(nsfc))/pint
                                ;q0=(1-a)*h2o(nsfc-1)+(1-b)*h2o(nsfc)
endif   

for I=nsfc-1,0,-1 do begin
    P1=P(I)
    Q1=H2O(I)
    qu=Q1
    ql=Q0
    wu=qu/(1.+qu)
    wl=ql/(1.+ql)
    hp=ALOG(P0/P1)
    x0=P0*wl*hp
    zeta=P1*wu/(P0*wl)
    if(abs(zeta-1.) gt 0.00001)then begin
        alza=ALOG(zeta)
        xint=x0*(zeta-1.)/alza
    endif else begin
        xint=x0*2./(3.-zeta)
    endelse
    QINT=QINT+xint*100./9.81
    P0=P1
    Q0=Q1
endfor
return
end
