function [cfov_anglexpos,cfov_npos,cfov_xxpos,cfov_yxpos,nomcsr,swxpos, ...
      cfov_nadirrad,cfov_angleinc,cfov_sample1_time, ...
      cfov_sample_time,scanStep,a2sRPY,s2scRPY,cfov_tysampxpos,cfov_chans] ...
    = getcfovCs(scanOrbitPath,pointingPath,antDerivedPath, ...
    calPath,chan2antPath,antDir,antType,csrID,nominalChannel,nomSatAlt)

% nomcsr is a 2-element vector.
lowalt = nomSatAlt-11;
hialt = nomSatAlt+11;
%swathwidth0 = 1700+25;%Extra distance added for optional extra CFOV positions
swathwidth0 = 1700; % Now trying to get minimum swath only

% Define parameter lookup tables for each nomcsr
csrIDList = {'15' '20' '25' '40' '50' 'WV' 'SM' 'ST' 'UA' '86'};
ncsr = 10;
nomcsr1set = [15 20 25 40 50 56 70 86 200 86];
nomcsr2set = [15 20 25 40 50 35 40 52 200 86];
% exagerated to get all channels by frequency
%special_freq_range = {[] [] [] [] [] [9 40] [5 90] [5 40]}; 
% now needed because there are too many special cases
special_channels = ...
    {{'36V' '36H' '60V' '89V' '89H' '166V' '183V'} ...
    {'18V' '18H' '18P' '18M' '18R' '18L' '23V' '23H' ...
    '36V' '36H' '36P' '36M' '60V' '89V' '89H' '166V' '183V'} ...
    {'18V' '18H' '23V' '23H' ...
    '36V' '36H' '60V' '89V' '89H' '166V' '183V'} ...
    {'10V' '10H' '18V' '18H' '23V' '23H' ...
    '36V' '36H' '60V' '89V' '89H' '166V' '183V'} ...
    {'6V' '6H' '10V' '10H' '18V' '18H' '23V' '23H' ...
    '36V' '36H' '60V' '89V' '89H' '166V' '183V'} ...
    {'10V' '10H' '10R' '10L' '18V' '18H' '18P' '18M' '18R' '18L' ...
    '23V' '23H' '36V' '36H' '36P' '36M' '89V' '89H'} ...
    {'6V' '6H' '10V' '10H' '18V' '18H' '36V' '36H' '89V' '89H'} ...
    {'6V' '6H' '10V' '10H' '10R' '10L' '18V' '18H' '18P' '18M' '18R' '18L' ...
    '23V' '23H' '36V' '36H' '36P' '36M'} ...
    {'60L'} ...
    {'6V' '6H' '10V' '10H'}};
products = {'CM' 'CM' 'CM' 'CM' 'CM' 'WindVect' 'SoilMoist' 'SST' 'UAVTP' ...
    'SST'};

% At 25 km, 50 km NOMCSR is barely Nyquist sampled.  Is this OK?  
% Is 24/12/6 km sampling OK?  Make nominal for now.  Must be integer
% multiples of minimum spacing.
sampset_nomSatAlt = [6 6 12 12 24 12 12 24 12 24];
% All CFOVs are sampled in every scan except 50 km
scanStepSet = [1 1 1 1 2 1 1 2 2 2];
% edrAltSet = [17 20 27 27 27 0 0 0 80 0]; % Derived set
edrAltSet = [15 15 15 15 15 0 0 0 80 0]; % More realistic (less stressing)
% Cell size is used for minimum swath width definition. 
% Wording in spec is "outermost edge of edr reporting cells"
cellWidthSet = min(nomcsr1set,nomcsr2set);

cellWidthSet(10) = 52;

icsr = find(strcmp(csrIDList,csrID));
if (isempty(icsr))
  error(['getcfovCs: invalid csrID: ' csrID])
end

nomcsr = [nomcsr1set(icsr) nomcsr2set(icsr)];
scanStep = scanStepSet(icsr);  
%cfov_freqs = special_freq_range{icsr};
cfov_chans = special_channels{icsr};
edrAlt = edrAltSet(icsr);
cellWidth = cellWidthSet(icsr);

% CFOV altitude selection:

% -Setting sampling using 833 alt is OK.  Sampling and along-scan footprint
% size scale with altitude in the same way:
% Alt:          816    833    850
% Ascan/vxsamp: 2.7293 2.7286 2.7279
% Ascan         66.2   68.0   69.9  -> +/-3%
% -Must set sampset_nomSatAlt to slightly less than Nyquist to insure good
% sampling despite variances in CFOV size.
% -CFOV size will also scale with alt for a fixed coef set.  If coefs
% are tuned such that nominal CFOV is met at 833 km, then true size
% will be +/-3% depending of true altitude.  Make 833 nominal for now.
% -1700 km swath width still has to be met at 816 km, so use 816 to
% regulate number of CFOV positions to meet swath requirement.

npnto2 = 100; nopatCk = 1; patsize = 3.5;
[beam_scananglexpos,cfov_nadirrad, ...
      jpatv,incj,xpntj,xmatj,ypntj,ymatj, ...
      tysamp,xzxpos,yzxpos,xlosxpos,ylosxpos,los_scananglexpos, ...
      aview_spec,ascan_spec,vxsamp,vxint,npos,design,aview,ascan, ...
      earthRadius,cfov_beta_nomSatAlt,nedt,k_cal,scan_start_time, ...
      sample_time,rps,integ_time,los_eia,beamID,iBeam, ...
      a2sRPY,s2scRPY,tysampxpos,dyintegxpos, ...
      antfiles,nantfiles,antfilefreq,freq2use] ...
    = getSensorCs_SCF(scanOrbitPath,pointingPath,antDerivedPath, ...
    calPath,chan2antPath,antDir,antType,nomSatAlt,nominalChannel,npnto2, ...
    nopatCk,patsize);

% Define sampling for CFOVs. CFOVs must be spaced <1/2 nomcsr, lowalt & hialt
% What scan angle gives 1700 km swathwidth? Lowest sat alt is limiting case

% lowest satellite distance to Earth center
lowaltplusR = lowalt + earthRadius; 
% highest EDR distance to Earth center
edrAltplusR = edrAlt + earthRadius; 
% Effective "EIA" for these conditions (i.e., zenith angle at edrAlt)
edrAltEIA = asin(lowaltplusR*sin(cfov_nadirrad)/edrAltplusR);
edrAltBeta = edrAltEIA - cfov_nadirrad;

cfov_beta_lowalt = ...
    asin(lowaltplusR*sin(cfov_nadirrad)/earthRadius) - cfov_nadirrad;
% Before 4/2/06 was:
% (equivalently since sin(cfov_beta_lowalt)=cfov_perp2arc/earthRadius):
%sample_sw = swathwidth0-sampset_nomSatAlt(icsr);
%swathangle ...
%    = rad2deg(asin( sin(sample_sw/(2*earthRadius)) /sin(cfov_beta_lowalt)));

% As of 4/2/06: Changed name to field-of-regard angle (twice swathangle!)
sample_sw = swathwidth0 - cellWidth;
forangle = 2*rad2deg(asin( sin(sample_sw/(2*earthRadius)) /sin(edrAltBeta)));

% To maintain nyquist along-scan, high alt is limiting case
% BUT:  Sampling at nomSatAlt is already better than Nyquist so Nyquist
% will be maintained at all altitudes. USE NOMSATALT for sampling.
% 4/2/06: Make angleinc round down to 100ths
cfov_perp2arc = earthRadius*sin(cfov_beta_nomSatAlt);
cfov_scancircle = cfov_perp2arc*2*pi;
cfov_angleinc = 0.01*floor(100*360*sampset_nomSatAlt(icsr)/cfov_scancircle);
cfov_npos = 1+2*ceil(0.5*forangle./cfov_angleinc);
cfov_anglexpos = ([0:cfov_npos-1]-(cfov_npos-1)/2)*cfov_angleinc;

if (1 == 0)
  % This is to re-calculate swathwidth on Earth for vertical prof EDR at alt
  edrsw = asin(sin(deg2rad(cfov_anglexpos))*sin(edrAltBeta))*2*earthRadius;
  edrsw = abs(edrsw)+cellWidth;
end

% Convert angle data into cfov_sample1_time and cfov_sample_time
% See Testbed notes 2/26/04:  Time/angle conversion uses nominalChannel
% data as time/angle reference point.
azdot = 360*rps; % degress/sec
cfov_sample_time = cfov_angleinc/azdot;
cfov_sample1_time = (cfov_anglexpos(1) - los_scananglexpos(1))/azdot ...
    + (scan_start_time + integ_time/2);

% Other parameters are at for nominal altitude.
[cfov_xxpos,cfov_yxpos] = ...
     xysphere([],90,cfov_anglexpos,[],cfov_beta_nomSatAlt,[],earthRadius);
% add satelite motion to yz
%cfov_yxpos = ...
%    cfov_yxpos+tysamp*[0:cfov_npos-1]*cfov_angleinc/360;
% Revised 12/23/05 to match getSensorCs_SCF
mu_earth = 3.986E5;
nomSatAltplusR = nomSatAlt+earthRadius;
ground_speed = earthRadius*sqrt(mu_earth/nomSatAltplusR^3);
cfov_yxpos = cfov_yxpos + ground_speed*(scan_start_time ...
    + cfov_sample_time/2 + [0:cfov_npos-1]*cfov_sample_time);

% align with nominalChannel LOS along-track at center of scan
% 2006/05/02:  Now only (yMidnom-yMidcfov)=-0.08-0.05 km so skip correction
%yMidnom = interp1(los_scananglexpos,ylosxpos,0.0);
%yMidcfov = interp1(cfov_anglexpos,cfov_yxpos,0.0);
%cfov_yxpos = cfov_yxpos + (yMidnom-yMidcfov);

ground_speed_xpos = ground_speed*cos(abs(cfov_xxpos)/earthRadius);
cfov_tysampxpos = ground_speed_xpos/rps;

% Swathwidth for low altitude orbit
% nomcsr(2)/2 is added to either side of swath to include up to edge of cell
[x,y] = ...
    xysphere([],90,cfov_anglexpos,[],cfov_beta_lowalt,[],earthRadius);
swxpos = 2*abs(x)+cellWidth;
