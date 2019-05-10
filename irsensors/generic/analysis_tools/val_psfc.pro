pro val_psfc, tData, bData, index, psfc, errIdx, nfor, iuo

print, 'Psfc - '
printf, iuo, 'Psfc - '
sidx = index.psfc.start
eidx = index.psfc.start + index.psfc.len - 1
Psfc = bData.profs[sidx:eidx,*]
for ifor = 0, n_elements(errIdx)-1 do begin
   print, errIdx[ifor]+1, ',', $
          abs(tData.profs[sidx:eidx,ifor]-bData.profs[sidx:eidx,ifor])/bData.profs[sidx:eidx,ifor]
   printf, iuo, errIdx[ifor]+1, ',', $
          abs(tData.profs[sidx:eidx,ifor]-bData.profs[sidx:eidx,ifor])/bData.profs[sidx:eidx,ifor]
   ;; print, format='(i5,a1,f10.4,a1,f10.4)', errIdx[ifor]+1, ',', $
   ;;        bData.profs[sidx:eidx,errIdx[ifor]], ',', tData.profs[sidx:eidx,errIdx[ifor]]
endfor
print, "<< val_psfc done >>"

end
