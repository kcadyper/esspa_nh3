pro val_iter, tData, bData, index, iter, errIdx, nfor, iuo

;- tag: Iter (10)
print, 'Iter - '
printf, iuo, 'Iter - '
;; print, "ifor, tData.iter[ifor], bData.iter[ifor]"
for ifor = 0, n_elements(errIdx)-1 do begin
   iter[errIdx[ifor]] = abs(tData.iter[errIdx[ifor]]-bData.iter[errIdx[ifor]])
   print, errIdx[ifor]+1, ',', iter[errIdx[ifor]]
   printf, iuo, errIdx[ifor]+1, ',', iter[errIdx[ifor]]
   ;; print, format='(3i6)', errIdx[ifor], tData.iter[errIdx[ifor]], bData.iter[errIdx[ifor]]
endfor

;- tag: noise (10)
;; print, 'Noise - '
;; print, max(rerror(tData.noise,bData.noise))

print, "<< val_iter done >>"

end
