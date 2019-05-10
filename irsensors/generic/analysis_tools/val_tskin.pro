pro val_tskin, tData, bData, index, tskin, errIdx, nfor, iuo

print, 'Tskin - '
printf, iuo, 'Tskin - '
sidx = index.tskin.start
eidx = index.tskin.start + index.tskin.len - 1
for ifor = 0, n_elements(errIdx)-1 do begin
   tskin[errIdx[ifor]] = abs(tData.profs[sidx:eidx,errIdx[ifor]]- $
                             bData.profs[sidx:eidx,errIdx[ifor]])/ $
                         bData.profs[sidx:eidx,errIdx[ifor]]
   print, errIdx[ifor]+1, ',', tskin[errIdx[ifor]]
   printf, iuo, errIdx[ifor]+1, ',', tskin[errIdx[ifor]]
endfor

print, "<< val_tskin done >>"

end
