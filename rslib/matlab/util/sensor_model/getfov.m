function [aview,ascan]=getfov(pat,xpnt,ypnt);
% Total vertical and horizontal extent of 3db contour of footprint.
% Size of rectangular box circumscribing contour.  See also gethsr.m

c = contourc(xpnt,ypnt,pat/max(pat(:)),[0.5 0.5]);

% Should be only ONE contour at 0.5; Fail if more
if ((size(c,2) > c(2,1)+1) | (c(1,1) ~= 0.5))
  error('getfov: Half-max contour level failed.')
end

xc = c(1,2:end);
yc = c(2,2:end);

aview = max(yc)-min(yc); 
ascan = max(xc)-min(xc);
