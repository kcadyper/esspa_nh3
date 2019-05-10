pro val_h2o, tData, bData, index, h2o, errIdx, nfor, iuo, psfc=psfc, gData=gData, uBound=uBound

  print, 'Water vapor profiles - '
  printf, iuo, 'Water vapor profiles - '
  sidx = index.H2O.start
  eidx = index.H2O.start + index.H2O.len - 1

  if     Keyword_Set(uBound) then pTop=uBound
  if not Keyword_Set(uBound) then pTop=0.

  if not Keyword_Set(gData) then begin
     for ifor = 0, n_elements(errIdx)-1 do begin
        eif=errIdx[ifor]
        press=reform(tData.pre[*,eif])
        levIn=where((press ge pTop) and (press le tData.profs[index.psfc.start,eif]))
        indH2O=index.H2O.start+levIn
        h2o[eif] = max(rerror(tData.profs[indH2O,eif], $
                                       bData.profs[indH2O,eif]))
        idmx=where(rerror(tData.profs[indH2O,eif],bData.profs[indH2O,eif]) eq h2o[eif])
        print, eif+1, ',', h2o[eif], ',', press[levIn[idmx[0]]]
        printf, iuo, eif+1, ',', h2o[eif], ',', press[levIn[idmx[0]]]
     endfor

     if n_elements(errIdx) le 10 then begin
        for ifor = 0, n_elements(errIdx)-1 do begin
           print, "ifor = ", errIdx[ifor]+1
           print, "level, pressure, base, test, (test-base), abs(test-base)/base"
           for idx = sidx, eidx do begin
              ;; if abs(tData.profs[idx,errIdx[ifor]]-bData.profs[idx,errIdx[ifor]])/bData.profs[idx,errIdx[ifor]] le threshold then continue
              ;; print, idx-sidx+1, (bData.pre)[idx-sidx,errIdx[ifor]], $
              ;;        abs(tData.profs[idx,errIdx[ifor]]-bData.profs[idx,errIdx[ifor]])/bData.profs[idx,errIdx[ifor]], $
              ;;        abs(tData.profs[idx,errIdx[ifor]]-bData.profs[idx,errIdx[ifor]])
              print, idx-sidx+1, ',', (bData.pre)[idx-sidx,errIdx[ifor]], ',',$
                     bData.profs[idx,errIdx[ifor]],',', tData.profs[idx,errIdx[ifor]], ',', $
                     (bData.profs[idx,errIdx[ifor]]-tData.profs[idx,errIdx[ifor]]), ',', $
                     abs(tData.profs[idx,errIdx[ifor]]-bData.profs[idx,errIdx[ifor]])/bData.profs[idx,errIdx[ifor]]
           endfor
        endfor
     endif
  endif 

  if Keyword_Set(gData) then begin
     h2o_sig_bot = fltarr(n_elements(errIdx))
     h2o_sig_bot_bkg_aerr = fltarr(n_elements(errIdx))
     h2o_sig_bot_rerr = fltarr(n_elements(errIdx))
     h2o_sig_bot_aerr = fltarr(n_elements(errIdx))
     print, "ifor, bottom_index, Psfc, pressure<=Psfc, sig, (sig-bkg), sig_aerr, sig_rerr"
     for ifor = 0, n_elements(errIdx)-1 do begin
        pres = bData.pre[*,errIdx[ifor]]
        bidx = where (pres le Psfc[ifor])
        h2o_sig_bot[ifor] = bData.profs[sidx+bidx[-1],errIdx[ifor]]
        h2o_sig_bot_bkg_aerr[ifor] = bData.profs[sidx+bidx[-1],errIdx[ifor]]-gData.profs[sidx+bidx[-1],errIdx[ifor]]
        h2o_sig_bot_aerr[ifor] = abs(tData.profs[sidx+bidx[-1],errIdx[ifor]]-bData.profs[sidx+bidx[-1],errIdx[ifor]])
        h2o_sig_bot_rerr[ifor] = h2o_sig_bot_aerr[ifor]/bData.profs[sidx+bidx[-1],errIdx[ifor]]
        ;; if h2o_sig_bot_aerr[ifor] le 4.E-5 then continue
        print, format='(i4,a1,i4,a1,f9.3,a1,f9.3,a1,f14.9,a1,f14.9,a1,e14.6,a1,e14.6)', $
               errIdx[ifor]+1,",",bidx[-1],",",Psfc[ifor],",",(pres)[bidx[-1]], ",", $
               h2o_sig_bot[ifor], ",", h2o_sig_bot_bkg_aerr[ifor], ",", $
               h2o_sig_bot_aerr[ifor], ",", h2o_sig_bot_rerr[ifor]
     endfor
     set_plot, 'ps'
     device, filename='rerr.vs.h2o_sig_bot.ps'
     plot, h2o_sig_bot, h2o_sig_bot_rerr, psym=1, xtitle='H2O_sig(bottom)',ytitle='Relative Error'
     device,/close
     device, filename='rerr.vs.h2o_sig_bot-h2o_bkg_bot.ps'
     plot, h2o_sig_bot_bkg_aerr, h2o_sig_bot_rerr, psym=1, xtitle='H2O_sig(bottom)-H2O_bkg(bottom)',ytitle='Relative Error'
     device,/close
     device, filename='aerr.vs.h2o_sig_bot.ps'
     plot, h2o_sig_bot, h2o_sig_bot_aerr, psym=1, xtitle='H2O_sig(bottom)',ytitle='Absolute Error'
     device,/close
     set_plot,'x'
  endif
  print, "<< val_h2o done >>"

end
