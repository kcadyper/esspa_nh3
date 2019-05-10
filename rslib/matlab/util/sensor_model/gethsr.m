function [aview,ascan]=gethsr(pat,xpnt,ypnt,inc,npnto2);
% Maximum extent of 3dB contour of footprint along the vertical and
% horizontal lines that cross at the footprint origin (where
% xpnt and ypnt = 0).  See also getfov.m.  Runs faster than 
% getfov.m for quicker intercomparison of footprints.

peak = max(pat(:));

cut = pat(:,npnto2+1)'/peak;
[gsort,igsort] = sort(abs(cut(1:npnto2)-0.5));
ig = sort(igsort(1:2));
i1 = ig(1); i2 = ig(2);
a1 = ypnt(i1) + (0.5-cut(i1))*(ypnt(i2)-ypnt(i1))/(cut(i2)-cut(i1));

[gsort,igsort] = sort(abs(cut(npnto2:end)-0.5));
ig = sort(igsort(1:2));
i1 = npnto2+ig(1)-1; i2 = npnto2+ig(2)-1;
a2 = ypnt(i1) + (0.5-cut(i1))*(ypnt(i2)-ypnt(i1))/(cut(i2)-cut(i1));

aview = abs(a1)+abs(a2);

cut = pat(npnto2+1,:)'/peak;
[gsort,igsort] = sort(abs(cut(1:npnto2)-0.5));
ig = sort(igsort(1:2));
i1 = ig(1); i2 = ig(2);
a1 = xpnt(i1) + (0.5-cut(i1))*(xpnt(i2)-xpnt(i1))/(cut(i2)-cut(i1));

[gsort,igsort] = sort(abs(cut(npnto2:end)-0.5));
ig = sort(igsort(1:2));
i1 = npnto2+ig(1)-1; i2 = npnto2+ig(2)-1;
a2 = xpnt(i1) + (0.5-cut(i1))*(xpnt(i2)-xpnt(i1))/(cut(i2)-cut(i1));

ascan = abs(a1)+abs(a2);

