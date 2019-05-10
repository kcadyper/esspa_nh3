; Code plots histograms for each Nh3 apriori of SFCT and surface NH3 from two algorithms:
; 1. SFCT run with LM and Hamming and SFTC with new sinc
; 2. SFCT run with no LM and Hamming and SFTC with new sinc


fn_lm = 'irsensors/crimss/run_sfct_hamming_nh3_new_sinc/retr.ir.exp.nc'
;fn_no_lm = 'irsensors/crimss/run_no_lm_sfct_hamming_nh3_new_sinc/retr.ir.exp.nc'
fn_no_lm = 'irsensors/crimss/run/retr.ir.exp.nc'

plot_type = 'PS'
if (plot_type eq 'PS') then begin
   set_plot,'PS'
   device,/color,xsize=20.,ysize=21.,yoffset=1,xoffset=1
   device,/Helvetica,/Bold
   !p.font = 0
endif


bhohls

!x.thick = 3
!y.thick = 3
!p.thick = 5
!p.charsize = 1.

!p.multi = [0,2,3]


; Read data
read_ncdf_jh,data_lm,filename=fn_lm
read_ncdf_jh,data_no_lm,filename=fn_no_lm

sfct_lm = data_lm.profiles[101,*]
nh3_lm = data_lm.profiles[304,*]*(28.9/17.0)*1.0e9
nh3_class = data_lm.atmosclass

sfct_no_lm = data_no_lm.profiles[101,*]
nh3_no_lm = data_no_lm.profiles[304,*]*(28.9/17.0)*1.0e9

; Create separate plots based on a priori selection
iclass = [3,2,1]
class_name = ['unp','mod','pol']
nh3_max = [5.,20.,200.]
nh3_bin = [0.5,2.,10.]

for ic=0,2 do begin
   andex = where(nh3_class eq iclass[ic],nm)

   histo_lm = histogram(sfct_lm[andex],binsize=5.,min=265.,max=335,loc=xarr)
   histo_no_lm = histogram(sfct_no_lm[andex],binsize=5.,min=265.,max=335,loc=xarr)

   plot,xarr,histo_lm,title='SFCT',/nodata,psym=10
   oplot,xarr,histo_lm,color=2,psym=10
   oplot,xarr,histo_no_lm,color=10,psym=10
   xyouts,0.4,0.9-0.35*ic,string(median(sfct_lm[andex]),format='(f5.1)')+'K',/normal,color=2
   xyouts,0.4,0.85-0.35*ic,string(median(sfct_lm[andex]),format='(f5.1)')+'K',/normal,color=10
   xyouts,0.5,0.95-0.35*ic,class_name[ic],/normal

   histo_lm = histogram(nh3_lm[andex],binsize=nh3_bin[ic],min=0.,max=nh3_max[ic],loc=xarr)
   histo_no_lm = histogram(nh3_no_lm[andex],binsize=nh3_bin[ic],min=0.,max=nh3_max[ic],loc=xarr)
   plot,xarr,histo_lm,title='NH3',/nodata,psym=10
   oplot,xarr,histo_lm,color=2,psym=10
   oplot,xarr,histo_no_lm,color=10,psym=10
   xyouts,0.8,0.9-0.35*ic,string(median(nh3_lm[andex]),format='(f4.1)')+'ppbv',/normal,color=2
   xyouts,0.8,0.85-0.35*ic,string(median(nh3_no_lm[andex]),format='(f4.1)')+'ppbv',/normal,color=10
endfor


if (plot_type eq 'PS') then device,/close_file
end
