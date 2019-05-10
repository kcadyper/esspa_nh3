; Code compares SFCT retrieval with Gennady tables (supposedly unapodized)  and no LM here in esspa10 and SFCT retrieval
;with Hamming weighted tables and data  and LM (run_sfct_ham_lm_nh3_ham_lm in esspa9)

; Code compares also compares SFCT retrieval with Gennady tables (supposedly unapodized)  and no LM here in esspa10 and SFCT retrieval
; with Gennady tables  tables and unapodized data  and LM (run in esspa9)
; This last test is to determine if the differences I see in the first test are due to the LM or the Hamming
; weighting or maybe both.



fn_no_apo = 'irsensors/crimss/run/retr.sfct.ir.nc'

fn_apo = '/project/p1913/esspa9/irsensors/crimss/run_sfct_ham_lm_nh3_ham_lm/retr.sfct.ir.nc'
;fn_no_apo_lm = '/project/p1913/esspa9/irsensors/crimss/run/retr.sfct.ir.nc'

; From here down apo means no_apo_lm for the second test
;fn_apo = fn_no_apo_lm

plot_type = 'PS'
if (plot_type eq 'PS') then begin
   set_plot,'PS'
   device,/color,xsize=20.,ysize=22.,yoffset=1,xoffset=1
   device,/Helvetica,/Bold
   !p.font = 0
endif


bhohls

!x.thick = 3
!y.thick = 3
!p.thick = 5
!p.charsize = 1.5

read_ncdf_jh,data_no_apo,filename=fn_no_apo
read_ncdf_jh,data_apo,filename=fn_apo

sfct_no_apo = data_no_apo.profiles[101,*]
sfct_apo = data_apo.profiles[101,*]
chisq_no_apo = data_no_apo.chisq
chisq_apo = data_apo.chisq


; plot SFCT diff 
; SFCT first
!x.title = 'Obs '
!y.title = 'Apo-No Apo SFCT'

!x.range = 0.
!y.range = 0.
plot,sfct_apo-sfct_no_apo,psym=1

; Limit analysis to cases with "good" SFCT retrieval (chisq lt 100.)
index = where (chisq_no_apo lt 100.,ngood)
print , 'no apo  ',ngood
plot,sfct_apo[index]-sfct_no_apo[index],psym=1

index = where (chisq_apo lt 100.,ngood)
print , ' apo  ',ngood
plot,sfct_apo[index]-sfct_no_apo[index],psym=1


; Now chisq scatter plots
!x.title = 'no apo chisq'
!y.title = 'apo chisq'

!x.range = 0.
!y.range = 0.
plot,chisq_no_apo[index],chisq_apo[index],psym=1

goto, skip
chisq_sfct = chisq_sfct[index]
chisq_nh3 = chisq_nh3[index]
atmos_class = data_nh3.atmosclass
atmos_class = atmos_class[index]
tcon = tcon[index]

!x.range = [0.,100.]
!y.range = [0.,300.]

; Create separate scatter plots based on a priori selection
iclass = [3,2,1]
class_name = ['unp','mod','pol']

for ic=0,2 do begin
   andex = where(atmos_class eq iclass[ic],nm)
   print, class_name[ic],nm
   plot,chisq_sfct[andex],chisq_nh3[andex],psym=2,title=class_name[ic]
   print,correlate(chisq_sfct[andex],chisq_nh3[andex])
   if (ic eq 2) then begin
      sndex = where(tcon[andex] gt 10 AND chisq_nh3[andex] lt 5.,ngood) 
      for ig=0,ngood-1 do print, index[andex[sndex[ig]]],tcon[andex[sndex[ig]]],chisq_nh3[andex[sndex[ig]]]
   endif
endfor

skip:

if (plot_type eq 'PS') then device,/close_file
end
