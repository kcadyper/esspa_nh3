; Code compares chisq from the SFCT retrieval with chisq from the Nh3 retrievals; hoping they are correlated


fn_nh3 = 'irsensors/crimss/run/retr.ir.exp.nc'
fn_sfct = '/project/p1913/esspa9/irsensors/crimss/run/retr.sfct.ir.nc'

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

read_ncdf_jh,data_nh3,filename=fn_nh3
read_ncdf_jh,data_sfct,filename=fn_sfct

chisq_sfct = data_sfct.chisq
chisq_nh3 = data_nh3.chisq
tcon = data_sfct.profiles[101,*]-data_sfct.profiles[100,*]


; plot scatter plot
!x.title = ' SFCT CHISQ '
!y.title = 'NH3 CHISQ'

!x.range = 0.
!y.range = 0.
plot,chisq_sfct,chisq_nh3,psym=2

; Limit analysis to cases with "good" SFCT retrieval (chisq lt 100.)
index = where (chisq_sfct lt 100.)

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


if (plot_type eq 'PS') then device,/close_file
end
