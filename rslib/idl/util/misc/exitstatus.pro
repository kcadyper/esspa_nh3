; $Id$

;+ 
; procedure to return exit code. It is designed to work with
; upper-level script.
;
; History: 
;   May 20, 2005, First Version.................S. Zaccheo, AER
;
; Copyright, 2005, AER, Inc. All rights reserved.
;-
pro exitStatus,CODE=code
if keyword_set(code) then print,"<PYIEXE ERROR:",code,">" $
else print,"<PYIEXE:OK>" 
end
