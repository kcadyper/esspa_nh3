pro val_chisq, tData, bData, index, chisq, errIdx, nfor, iuo

;- tag: ChiSQ
print, 'Chisq - '
printf, iuo, 'Chisq - '
tChisq = fltarr(nfor)
bChisq = fltarr(nfor)
;;;Use this block if both baseline and test files are "sig" version
;; for ifor = 0, (size(tData.chisq,/dim))[1]-1 do begin
;;    tChisq[errIdx[ifor]] = tData.chisq[tData.iter[errIdx[ifor]]-1,errIdx[ifor]]
;;    bChisq[errIdx[ifor]] = bData.chisq[bData.iter[errIdx[ifor]]-1,errIdx[ifor]]
;; endfor
;; print, "ifor, test-Chisq, base-Chisq"
;; print, "ifor, abs(sig-esdr), sig_Chisq, esdr_Chisq"
print, "ifor, base, test, (test-base), abs(test-base)/base"
printf, iuo, "ifor, base, test, (test-base), abs(test-base)/base"
for ifor = 0, n_elements(errIdx)-1 do begin
   if (size(bData.chisq))[1] gt 1 then $
      bChisq[errIdx[ifor]] = bData.chisq[bData.iter[errIdx[ifor]]-1,errIdx[ifor]] $ ;"sig" version
   else $
      bChisq[errIdx[ifor]] = bData.chisq[0,errIdx[ifor]] ;"esdr" version
   if (size(tData.chisq))[1] gt 1 then $
      tChisq[errIdx[ifor]] = tData.chisq[tData.iter[errIdx[ifor]]-1,errIdx[ifor]] $ ;"sig" version
   else $
      tChisq[errIdx[ifor]] = tData.chisq[0,errIdx[ifor]] ;"esdr" version
   chisq[errIdx[ifor]] = tChisq[errIdx[ifor]]-bChisq[errIdx[ifor]]
   ;; print, format='(i4,f10.6)', errIdx[ifor], chisq[errIdx[ifor]]
   ;; if chisq[errIdx[ifor]] le threshold then continue
   print, format='(i4,a1,f16.8,a1,f16.8,a1,f16.8,a1,f16.8,a1,f16.8)', $
          errIdx[ifor]+1, ',', bChisq[errIdx[ifor]], ',', tChisq[errIdx[ifor]], ',', $
          chisq[errIdx[ifor]], ',', abs(chisq[errIdx[ifor]])/bChisq[errIdx[ifor]]
   printf, iuo, format='(i4,a1,f16.8,a1,f16.8,a1,f16.8,a1,f16.8,a1,f16.8)', $
          errIdx[ifor]+1, ',', bChisq[errIdx[ifor]], ',', tChisq[errIdx[ifor]], ',', $
          chisq[errIdx[ifor]], ',', abs(chisq[errIdx[ifor]])/bChisq[errIdx[ifor]]
endfor

; as specified in standard_crimss_config.in
chisqconvMW = 0.9   ; value of chisq needed for convergence in MW
chisqconv   = 0.7   ; value of chisq needed for convergence in IR
threshold = chisqconvMW ;chisqconv
sidx = where (bChisq ge threshold, scnt)
eidx = where (tChisq ge threshold, ecnt)
print, "=== Total #ChiSq > ", threshold, " ===="
print, scnt, ecnt
print, "=== diverging indices ===="
print, iuo, "=== Total #ChiSq > ", threshold, " ===="
print, iuo, scnt, ecnt
print, iuo, "=== diverging indices ===="
;increment of 1 to start off from 1-based
if scnt ne ecnt then stop, "[err]: base and test mismatched in ChiSq outliers"
for ix = 0, n_elements(sidx)-1 do begin
   print, format='(i4,a1,i4,a1,f12.6,a1,f12.6)',sidx[ix]+1,',',eidx[ix]+1,',',bChisq[sidx[ix]],',',tChisq[eidx[ix]] 
   printf, iuo, format='(i4,a1,i4,a1,f12.6,a1,f12.6)',sidx[ix]+1,',',eidx[ix]+1,',',bChisq[sidx[ix]],',',tChisq[eidx[ix]] 
endfor

;; errIdx = [40,84,97,192,204,209,270,344,435]
;; print, "=== remaining profiles with H2O error greater than 1E-5 ===="
;; for ix = 0, n_elements(errIdx)-1 do begin
;;    print, format='(i4,a1,i4,a1,i4,a1,f12.6,a1,f12.6)', $
;;           errIdx[ix]+1,',',sData.iter[h2oidx[ix]],',',eData.iter[h2oidx[ix]], $
;;           ',', bchisq[h2oidx[ix]], ',', tchisq[h2oidx[ix]] 
;; endfor
;; plot, bchisq[sidx], psym=2
;; oplot, tchisq[eidx], psym=5

print, "<< val_chisq done >>"

end
