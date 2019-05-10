;+
; Relative error
; yhe@aer, 2016
;-
function rerror, test, base
  if not array_equal(size(test,/dim),size(base,/dim)) then begin
     stop, "[rerror]: test and base in different shapes"
  endif
  
  delta = float(test) - float(base)
  error = fltarr(size(test,/dim))
  idx0 = where (test eq 0 and base eq 0, cnt0)
  if cnt0 gt 0 then error[idx0] = 0
  idx1 = where (test ne 0 and base eq 0, cnt1)
  if cnt1 gt 0 then error[idx1] = abs(delta[idx1])/test[idx1]
  idx2 = where (test ne 0 and base ne 0, cnt2)
  if cnt2 gt 0 then error[idx2] = abs(delta[idx2])/base[idx2]
  return, error
end
