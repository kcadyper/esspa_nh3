function qmin = findscale1(xAndyint,scalechanstr,npnto2,zenrad_peak, ...
    beta_peak,alt,Rearth,xpntref,ypntref,incref,refhsr,antpat_format, ...
    antdir);
% Used to find the scaling factors needed to scale a given
% antenna pattern to the desired FOV size

if (strcmp(antpat_format,'elaz'))
  refpat = cmis_refpat_elaz(antdir,scalechanstr,npnto2, ...
      incref,xAndyint(1),xAndyint(2),zenrad_peak, ...
      beta_peak,alt,Rearth);
else
  error('findscale1: Only elaz antenna format supported')
  [refpat tmp] = cmis_effant(scalechanstr,npnto2, ...
      incscale,tysamp,vxinttime_in,scan2view,zenrad_peak, ...
      beta_peak,aviewscalechan,ascanscalechan,alt,Rearth,1);
end

[aview,ascan] = gethsr(refpat,xpntref,ypntref,incref,npnto2);

qmin = sum((refhsr - [aview ascan]).^2);
%disp(['findscale1:' num2str([aview ascan qmin  xAndyint])])

%FMINSEARCH uses
%    these options: Display, TolX, TolFun, MaxFunEvals, and MaxIter.
% Display - Level of display [ off | iter | notify | final ]  
% TolX - Termination tolerance on X [ positive scalar ]
% TolFun - Termination tolerance on the function value [ positive scalar ]
% MaxFunEvals - Maximum number of function evaluations allowed 
%                      [ positive integer ]
% MaxIter - Maximum number of iterations allowed [ positive integer ]
