pro val_temp, tData, bData, index, temp, errIdx, nfor, iuo

print, 'Temperature profiles - '
printf, iuo, 'Temperature profiles - '
sidx = index.temp.start
eidx = index.temp.start + index.temp.len - 1
for ifor = 0, n_elements(errIdx)-1 do begin
   temp[errIdx[ifor]] = max(rerror(tData.profs[sidx:eidx,errIdx[ifor]], $
                                   bData.profs[sidx:eidx,errIdx[ifor]]))
   print, errIdx[ifor]+1, ',', temp[errIdx[ifor]]
   printf, iuo, errIdx[ifor]+1, ',', temp[errIdx[ifor]]
endfor

if n_elements(errIdx) le 10 then begin
   for ifor = 0, n_elements(errIdx)-1 do begin
      print, "ifor = ", errIdx[ifor]+1
      print, "level, base, test, (test-base), abs(test-base)/base"
      for idx = sidx, eidx do begin
         print, idx-sidx+1,',', (bData.pre)[idx-sidx,errIdx[ifor]], ',', $
                bData.profs[idx,errIdx[ifor]], ',', tData.profs[idx,errIdx[ifor]], ',', $
                (bData.profs[idx,errIdx[ifor]]-tData.profs[idx,errIdx[ifor]]), ',', $
                abs(bData.profs[idx,errIdx[ifor]]-tData.profs[idx,errIdx[ifor]])/bData.profs[idx,errIdx[ifor]]
      endfor
   endfor
endif
;; for ifor = 0, n_elements(errIdx)-1 do begin
;;    print, "ifor = ", errIdx[ifor]+1
;;    print, "level, pressure, T_sig, T_esdr, abs(sig-esdr)/sig"
;;    for idx = sidx, eidx do begin
;;       ;; if abs(tData.profs[idx,errIdx[ifor]]-bData.profs[idx,errIdx[ifor]])/bData.profs[idx,errIdx[ifor]] le threshold then continue
;;       ;; print, idx-sidx+1, (bData.pre)[idx-sidx,errIdx[ifor]], $
;;       ;;        abs(tData.profs[idx,errIdx[ifor]]-bData.profs[idx,errIdx[ifor]])/bData.profs[idx,errIdx[ifor]], $
;;       ;;        abs(tData.profs[idx,errIdx[ifor]]-bData.profs[idx,errIdx[ifor]])
;;       print, idx-sidx+1, ',', (bData.pre)[idx-sidx,errIdx[ifor]], ',',$
;;              bData.profs[idx,errIdx[ifor]],',', tData.profs[idx,errIdx[ifor]], ',', $
;;              abs(tData.profs[idx,errIdx[ifor]]-bData.profs[idx,errIdx[ifor]])/bData.profs[idx,errIdx[ifor]]
;;    endfor
;; endfor

print, "<< val_temp done >>"

end
