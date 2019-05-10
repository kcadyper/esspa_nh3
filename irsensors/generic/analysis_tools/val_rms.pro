pro val_rms, tData, bData, index, rms, errIdx, nfor, iuo

;; baseline run by SIG version contains RMS and ChiSq as arrays of
;; dimension of (nIter, nFOR), while the modularized ESDR code
;; generates RMS and ChiSq as arrays of dimension of (1, nFOR)
;- tag: RMS
print, 'RMS - '
printf, iuo, 'RMS - '
tRMS = fltarr(nfor)
bRMS = fltarr(nfor)
;;;Use this block if both baseline and test files are "sig" version
;; for ifor = 0, n_elements(errIdx)-1 do begin
;;    tRMS[errIdx[ifor]] = tData.rms[tData.iter[errIdx[ifor]]-1,errIdx[ifor]]
;;    bRMS[errIdx[ifor]] = bData.rms[bData.iter[errIdx[ifor]]-1,errIdx[ifor]]
;; endfor
;; print, "errIdx[ifor]+1, test-RMS,  base-RMS[errIdx[ifor]]"
print, "ifor, base, esdr, (esdr-base), abs(esdr-base)/base"
printf, iuo, "ifor, base, esdr, (esdr-base), abs(esdr-base)/base"
for ifor = 0, n_elements(errIdx)-1 do begin
   if (size(bData.rms))[1] gt 1 then $
      bRMS[errIdx[ifor]] = bData.rms[bData.iter[errIdx[ifor]]-1,errIdx[ifor]] $;"sig" version
   else $
      bRMS[errIdx[ifor]] = bData.rms[0,errIdx[ifor]] ;"esdr" version
   if (size(tData.rms))[1] gt 1 then $
      tRMS[errIdx[ifor]] = tData.rms[tData.iter[errIdx[ifor]]-1,errIdx[ifor]] $;"sig" version
   else $
      tRMS[errIdx[ifor]] = tData.rms[0,errIdx[ifor]] ;"esdr" version
   rms[errIdx[ifor]] = abs(tRMS[errIdx[ifor]]-bRMS[errIdx[ifor]])
   ;; print, format='(i4,f10.6)', errIdx[ifor]+1, rms[errIdx[ifor]]
   print, errIdx[ifor]+1, ',', $
          bRMS[errIdx[ifor]], ',', tRMS[errIdx[ifor]], ',', $
          (bRMS[errIdx[ifor]]-tRMS[errIdx[ifor]]), ',', $
          abs(bRMS[errIdx[ifor]]-tRMS[errIdx[ifor]])/bRMS[errIdx[ifor]]
   printf, iuo, errIdx[ifor]+1, ',', $
          bRMS[errIdx[ifor]], ',', tRMS[errIdx[ifor]], ',', $
          (bRMS[errIdx[ifor]]-tRMS[errIdx[ifor]]), ',', $
          abs(bRMS[errIdx[ifor]]-tRMS[errIdx[ifor]])/bRMS[errIdx[ifor]]
endfor

print, "<< val_rms done >>"

end
