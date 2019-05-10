; Code plots the chisq of the SFCT retrieval vs the difference between the L2 SFCT from CLIMCAPS
; and the retrieved SFCT
; Assuming tha large CHISQ are associated with large changes, indicating the presence of clouds.
; Actually it seems that large chisq indicate no convergent solution, and therefore possibly clouds,
; since the temperature would have to change by more than I am allowing it to to get a correct solution.

; While playing around with the covariance calculation for the SFCT I found I had used the wrong file  
; as input to BuildCovariance.in.

fn_l2 = 'irsensors/crimss/run/l2.safe.nc'
fn_retv = 'irsensors/crimss/run/retr.sfct.ir.nc'

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

read_ncdf_jh,data_l2,filename=fn_l2
read_ncdf_jh,data_retv,filename=fn_retv

sfct_l2 = data_l2.profiles[101,*]
sfct_retv = data_retv.profiles[101,*]
diff = sfct_retv-sfct_l2

chisq = data_retv.chisq

; plot histogram of chisq
histo = histogram(chisq,loc=xarr,min=0.,max=100.,binsize=5.)
!x.title = 'CHISQ'
!y.title = 'Count'
plot,xarr,histo,psym=10
print, total(histo)

histo = histogram(chisq,loc=xarr,min=0.,max=1000.,binsize=100.)
!x.title = 'CHISQ'
!y.title = 'Count'
plot,xarr,histo,psym=10
print, total(histo)


; plot scatter plot
!x.title = ' Temperature Change [K]'
!y.title = 'CHISQ'

plot,diff,chisq,psym=2


if (plot_type eq 'PS') then device,/close_file
end
