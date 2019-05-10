pro val_ozone, tData, bData, index, o3, errIdx, nfor, iuo, psfc=psfc, gData=gData

print, 'Ozone profiles - '
printf, iuo, 'Ozone profiles - '
sidx = index.O3.start
eidx = index.O3.start + index.O3.len - 1

for ifor = 0, n_elements(errIdx)-1 do begin
   o3[errIdx[ifor]] = max(rerror(tData.profs[sidx:eidx,errIdx[ifor]],bData.profs[sidx:eidx,errIdx[ifor]]))
   print, errIdx[ifor]+1, ',', o3[errIdx[ifor]]
   printf, iuo, errIdx[ifor]+1, ',', o3[errIdx[ifor]]
endfor

if n_elements(errIdx) le 10 then begin
   for ifor = 0, n_elements(errIdx)-1 do begin
      print, "ifor = ", errIdx[ifor]+1
      ;; print, "testData.profs, baseData.profs"
      print, "level, pressure, base, test, abs(base-test)/base"
      for idx = sidx, eidx do begin
         ;; if abs(tData.profs[idx,errIdx[ifor]]-bData.profs[idx,errIdx[ifor]])/bData.profs[idx,errIdx[ifor]] le threshold then continue
         ;; print, idx-sidx+1,(bData.pre)[idx-sidx,errIdx[ifor]],abs(tData.profs[idx,errIdx[ifor]]-bData.profs[idx,errIdx[ifor]])
         print, idx-sidx+1, ',', (bData.pre)[idx-sidx,errIdx[ifor]], ',',$
                bData.profs[idx,errIdx[ifor]],',', tData.profs[idx,errIdx[ifor]], ',', $
                bData.profs[idx,errIdx[ifor]]-tData.profs[idx,errIdx[ifor]], ',', $
                abs(bData.profs[idx,errIdx[ifor]]-tData.profs[idx,errIdx[ifor]])/bData.profs[idx,errIdx[ifor]]
      endfor
   endfor
endif

print, "<< val_ozone done >>"

end
