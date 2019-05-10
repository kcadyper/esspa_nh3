function [cout,aout] = xysphere(alpha,beta,gam,a,b,c,Rearth);
% Gives aout = x, bout = y where is in gam=0 direction and x is 
% perpendicular to y. Alpha, a, and c are unused. abc in rad.
% alpha, beta, gam in deg.

% See notes in "CMIS Footprint Matching" 11/19/2000

betarad = deg2rad(beta); %90 degrees
gamrad = deg2rad(gam);

cout = asin(sin(b).*sin(gamrad)./sin(betarad));

aout = acos( ...
    (cos(b).*cos(cout) - sin(b).*sin(cout).*cos(betarad).*cos(gamrad)) ...
      ./(1 - sin(b).*sin(cout).*sin(betarad).*sin(gamrad)) );

cout = cout*Rearth;
aout = aout*Rearth;

return