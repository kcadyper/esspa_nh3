;+
;
; IDL procedure to test the IR/MW retrievals against the
; baseline. Strictly speaking, baseline is what has been created by a
; standard version of the IR/MW retrieval code. But this script could
; be used in a loosely defined terms to compare any pairs of retrieval
; results. 
;
; Inputs:
;
;   base = baseline file
;   test = test file
;   bkg  = background or guess file (optional)
; 
; Outputs:
;    
;   None, except that when outfile is specified the statistics of
;   analysis is written out to the file in netCDF format.
;
; History:
; 
;   01/13/2017, initial, yhe@aer.com
;
;-
pro val_retr, base, test, bkg=bkg, outfile=outfile
;
;-- Baseline under p1855 works with the interpolation algorithm before
;   optimization, while local baseline works for the optimized one.
;   Baseline results were generated using igeo=2 (land) and "sig"
;   version of the retrieval code, of which the results are the same
;   as those in retr_sig.20161129.02.igeo=2.nc under the testPath.
;--
  outasc='./val_retr.out'
  print, 'writing ASCCI to ',outasc
  openw,iuo,outasc,/get_lun
  print, 'baseline file = ', base
  print, '    test file = ', test
  printf, iuo, 'baseline file = '+strtrim(base,2)
  printf, iuo, '    test file = '+strtrim(test,2)
  if Keyword_Set(bkg) then $
     print, '     bkg file = ', bkg
;--
  threshold = 1.E-5
;--
  bData = rd_Retr(ncfile=base)
  tData = rd_Retr(ncfile=test)
  if Keyword_Set(bkg) then gData = rd_Retr(ncfile=bkg)

;;may be read from the files
  index = bData.index

  ntagBase = n_tags(bData)
  tagsBase = tag_names(bData)

  ntagTest = n_tags(tData)
  tagsTest = tag_names(tData)

  if Keyword_Set(gData) then begin
     ntagBkg = n_tags(gData)
     tagsBkg = tag_names(gData)
  endif

  if ntagBase ne ntagTest then begin
     print, "[err]: mismatch number of tags", ntagBase, ' != ', ntagTest
     stop
  endif

  for itg = 0, ntagTest-1 do begin
     if tagsBase[itg] ne tagsTest[itg] then begin
        print, "[err]: mismatch tag - ", tagsBase[itg], ' != ', tagsTest[itg]
        stop
     endif
  endfor


  if Keyword_Set(gData) then begin
     if ntagBase ne ntagBkg then begin
        print, "[err]: mismatch number of tags", ntagBase, ' != ', ntagBkg
        stop
     endif
  endif

  if Keyword_Set(gData) then begin
     for itg = 0, ntagTest-1 do begin
        if tagsTest[itg] ne tagsBkg[itg] then begin
           print, "[err]: mismatch tag - ", tagsTest[itg], ' != ', tagsBkg[itg]
           stop
        endif
     endfor
  endif

  print, "calculating relative errors ... "

; Define vectors for storing differences
  if (size(tData.profs,/dim))[1] ne (size(bData.profs,/dim))[1] then $
     print, "[warning]: #profiles different b/w test and baseline ... min() selected", $
            (size(tData.profs,/dim))[1],(size(bData.profs,/dim))[1]

  nfor = min([(size(tData.profs,/dim))[1],(size(bData.profs,/dim))[1]])
  forIndex = indgen(min([(size(tData.profs,/dim))[1],(size(bData.profs,/dim))[1]]))
  errIdx = forIndex
;;
;;;- specific tests
;; errIdx = [17,40,47,62,84,97,118,180,192,204,209,223,227,270,298,315,344,367,378,388,435,452,457]
;; errIdx = [227]
  temp = fltarr(nfor)
  tskin = fltarr(nfor)
  psfc = fltarr(nfor)
  h2o = fltarr(nfor)
  o3 = fltarr(nfor)
  clw = fltarr(nfor)
  wsfc = fltarr(nfor)
  rms = fltarr(nfor)
  chisq = fltarr(nfor)
  iter = intarr(nfor)

  val_temp, tData, bData, index, temp, errIdx, nfor, iuo
  val_tskin, tData, bData, index, tskin, errIdx, nfor, iuo
  val_psfc, tData, bData, index, psfc, errIdx, nfor, iuo
  val_h2o, tData, bData, index, h2o, errIdx, nfor, iuo, psfc=psfc, gData=gData, uBound=300.
;  if index.O3.len gt 0 then val_ozone, tData, bData, index, o3, errIdx, nfor, iuo
;  val_clw, tData, bData, index, clw, errIdx, nfor, iuo
;  val_ciw, tData, bData, index, ciw, errIdx, nfor, iuo
;  val_winds, tData, bData, index, wsfc, errIdx, nfor, iuo
;  val_rms, tData, bData, index, rms, errIdx, nfor, iuo
;  val_chisq, tData, bData, index, chisq, errIdx, nfor, iuo
;  val_iter, tData, bData, index, iter, errIdx, nfor, iuo

;- tag: noise (10)
;; print, 'Noise - '
;; print, max(rerror(tData.noise,bData.noise))

  if Keyword_Set(outfile) then begin ;write out to a file
     stats = {forindex: forindex, $
              temp: temp, $
              tskin: tskin, $
              h2o: h2o, $
              o3: o3, $
              clw: clw, $
              rms: rms, $
              chisq: chisq, $
              iter: iter $
             }

     wr_stats, outfile, stats
  endif
  close,iuo
  print, "<< val_retr done >>"
end
