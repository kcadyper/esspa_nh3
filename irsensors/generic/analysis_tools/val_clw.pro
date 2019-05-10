pro val_clw, tData, bData, index, clw, errIdx, nfor, iuo

print, 'Liquid Cloud - '
printf, iuo, 'Liquid Cloud - '
;; errIdx = [17,40,97,223,298,315,344,347,378,388,398]
sidx = index.clw.start
eidx = index.clw.start + index.clw.len - 1
;; for ifor = 0, nfor-1 do begin
;;    clw[ifor] = max(rerror(tData.profs[sidx:eidx,ifor],bData.profs[sidx:eidx,ifor]))
;;    print, ifor, clw[ifor]
   ;; print, "ifor = ", ifor 
   ;; for idx = sidx, eidx do begin
   ;;    print, tData.profs[idx,ifor],bData.profs[idx,ifor]
   ;; endfor
;; endfor

print, "ifor, max(abs(test-base)), max(abs(test-base)/base)"
printf, iuo, "ifor, max(abs(test-base)), max(abs(test-base)/base)"
for ifor = 0, n_elements(errIdx)-1 do begin
   rerr = rerror(tData.profs[sidx:eidx,errIdx[ifor]],bData.profs[sidx:eidx,errIdx[ifor]])
   clw[errIdx[ifor]] = max(rerr)
   ix = where (rerr eq max(rerr))
   print, errIdx[ifor]+1, ',', clw[errIdx[ifor]], ',', rerr[ix]*(bData.profs[sidx:eidx,errIdx[ifor]])[ix]
   printf, iuo, errIdx[ifor]+1, ',', clw[errIdx[ifor]], ',', rerr[ix]*(bData.profs[sidx:eidx,errIdx[ifor]])[ix]
endfor

print, "<< val_clw done >>"

end
