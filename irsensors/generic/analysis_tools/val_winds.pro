pro val_winds, tData, bData, index, wsfc, errIdx, nfor, iuo

print, 'Surface Winds - '
printf, iuo, 'Surface Winds - '
sidx = index.wsfc.start
eidx = index.wsfc.start + index.wsfc.len - 1
for ifor = 0, n_elements(errIdx)-1 do begin
   print, errIdx[ifor]+1, ',', max(rerror(tData.profs[sidx:eidx,errIdx[ifor]], $
                                     bData.profs[sidx:eidx,errIdx[ifor]]))
   printf, iuo, errIdx[ifor]+1, ',', max(rerror(tData.profs[sidx:eidx,errIdx[ifor]], $
                                     bData.profs[sidx:eidx,errIdx[ifor]]))
   ;; print, "ifor = ", ifor
   ;; for idx = sidx, eidx do begin
   ;;    print, tData.profs[idx,ifor],bData.profs[idx,ifor]
   ;; endfor
endfor

print, "<< val_winds done >>"

end
