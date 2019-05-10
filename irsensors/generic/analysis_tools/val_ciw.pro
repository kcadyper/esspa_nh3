pro val_ciw, tData, bData, index, ciw, errIdx, nfor, iuo

print, 'Ice Cloud - '
printf, iuo, 'Ice Cloud - '
sidx = index.ciw.start
eidx = index.ciw.start + index.ciw.len - 1
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

print, "<< val_ciw done >>"

end
