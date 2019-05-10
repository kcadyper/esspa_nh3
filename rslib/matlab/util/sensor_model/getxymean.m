function [x,y] = getxymean(pat,xmat,ymat);

% Center of gravity within 3dB contour.  Finds LOS projected location.

block = pat/max(pat(:)) >= 0.5;
tmppat = pat.*block;
tmpx = xmat.*tmppat; tmpy = ymat.*tmppat;
x = sum(tmpx(:))/sum(tmppat(:));
y = sum(tmpy(:))/sum(tmppat(:));
