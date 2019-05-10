function [beam_scananglexpos,beam_nadirrad, ...
      jpatv,incj,xpntj,xmatj,ypntj,ymatj, ...
      tysamp,xzxpos,yzxpos,xlosxpos,ylosxpos,los_scananglexpos, ...
      aview_spec,ascan_spec,vxsamp,vxint,npos,design,aview,ascan, ...
      earthRadius,beam_beta,nedt_spec,k_cal,scan_start_time, ...
      sample_time,rps,integ_time,los_eia,beamID,iBeam,a2sRPY,s2scRPY, ...
      tysampxpos,dyintegxpos,antfiles,nantfiles,antfilefreq,freq2use] ...
    = getSensorCs_SCF(scanOrbitPath,pointingPath,antDerivedPath,calPath, ...
    chan2antPath,antDir,antType,alt,chanID,npnto2,nopatCk,patsize)
% Notes: Returns aview_spec, ascan_spec because these are *not* alt dependent.

% Get sensor constants for chanID from files

[orbitIncl,orbitAlt,earthRadius,muEarthKm,missingVal,rpm,bore_az_t0, ...
      nchan,design,bore2beam_scanangle,scan_start_time,integ_time, ...
      sample_time,n_samples,beam_nadir,aview_spec,ascan_spec, ...
      nedt_spec,fscaleincj,bore2los_scanangle,los_nadir, ...
      k_cal,beamID,iBeam,a2sRPY,s2scRPY,aview_nom,ascan_nom] ...
    = getConst(scanOrbitPath,antDerivedPath,pointingPath,calPath,chanID);

% Convert ms to s
integ_time = 1.0E-3*integ_time;
sample_time = 1.0E-3*sample_time;

% Derive constants

rps = rpm/60;
azdot = 360*rps;
altplusR = earthRadius+alt;
mu_earth = 3.986E5;
ground_speed = earthRadius*sqrt(mu_earth/altplusR^3);
tysamp = ground_speed/rps;
beam_nadirrad = deg2rad(beam_nadir);
beam_eiarad = asin(altplusR*sin(beam_nadirrad)/earthRadius);
% beta (angle at earth center connecting satellite to earth point)
beam_beta =  beam_eiarad - beam_nadirrad;
% perpendicular distance from earth-center-to-satellite line to arc
perp2arc = earthRadius*sin(beam_beta);
% circumference of entire scan arc (i.e., satellite standing still)
scancircle = perp2arc*2*pi;
% integration distance (not time) along scan
vxint = scancircle*integ_time*rps;

% Get channel peak scan angles for positioning footprints

% Calculate angle over which sample is taken
integ_scanangle = integ_time*azdot;
% Find scan angle of first sample beam center
beam_sample1_az = scan_start_time*azdot + bore_az_t0 ...
  + bore2beam_scanangle + integ_scanangle/2;

% Calculate all scan angles 
sample_scanangle = sample_time*azdot;
beam_scananglexpos = beam_sample1_az + sample_scanangle*[0:n_samples-1];

% X and Y spherical earth surface coords at each position along arc.
% Y is measured from sub-satellite point along satellite ground track.
% X is measured perp. to Y axis to position.
% Analogy: Y = distance along meridian, X = dist along parallel. 
%
% Remember: pattern is positioned based on pattern peak. 
[xzxpos,yzxpos] = ...
   xysphere([],90,beam_scananglexpos,[],beam_beta,[],earthRadius);
%% add satelite motion to yz
%%yzxpos = yzxpos+tysamp*[0:n_samples-1]*sample_time*rps;
%% add reflector phase to yz
%%yzxpos = yzxpos+tysamp*(reflector-1)/2;
% Revised:  Motion is relative to time0
yzxpos = yzxpos + ground_speed*(scan_start_time ...
  + sample_time/2 + [0:n_samples-1]*sample_time);

% dyinteg and tysampxpos (new var.) now include fact
% that earth is spherical, not cylindrical.  Allows for better 
% LOCAL footprint positioning.  See notes 3/18/05
ground_speed_xpos = ground_speed*cos(abs(xzxpos)/earthRadius);
dyintegxpos = ground_speed_xpos*integ_time;
tysampxpos = ground_speed_xpos/rps;

% Derive positions of LOS, which are offset from beam peak positions
los_scananglexpos = beam_scananglexpos ...
   + bore2los_scanangle - bore2beam_scanangle;
los_nadirrad = deg2rad(los_nadir);
los_eiarad = asin(altplusR*sin(los_nadirrad)/earthRadius);
los_beta = los_eiarad - los_nadirrad;
[xlosxpos,ylosxpos] = ...
   xysphere([],90,los_scananglexpos,[],los_beta,[],earthRadius);
ylosxpos = ylosxpos + ground_speed*(scan_start_time ...
  + sample_time/2 + [0:n_samples-1]*sample_time);

% Read mapping from channels to frequencies etc. needed for antenna pattern
fid = fopen(chan2antPath,'r');
tmp = fgets(fid); tmp = fgets(fid); tmp = fgets(fid); tmp = fgets(fid);
nchanLoc = fscanf(fid,'%d',1);
clear chanIDsLoc chanIDs2use freq aviewSpecxChan ascanSpecxChan
for ichan=1:nchanLoc
  chanIDsLoc{ichan} = fscanf(fid,'%s',1);
  chanIDs2use{ichan} = fscanf(fid,'%s',1);
  freq(ichan) = fscanf(fid,'%f',1);
  aviewSpecxChan(ichan) = fscanf(fid,'%f',1);
  ascanSpecxChan(ichan) = fscanf(fid,'%f',1);
end
fclose(fid);

jchan = strmatch(chanID,chanIDsLoc);
% aview/ascan_spec info seems redundant but may not be provided by
% getConst (when antDerivedPath is undefined).
aview_spec = aviewSpecxChan(jchan);
ascan_spec = ascanSpecxChan(jchan);

chanID2use = chanIDs2use{jchan};
jchan2use = strmatch(chanID2use,chanIDsLoc);
if (~strmatch(chanID2use,chanID))
  disp(['warn[getSensorCs_SCF]: Using ' chanID2use ' antenna pattern for ' ...
	chanID]);
end
freq2use = freq(jchan2use);

% Confirm that chanID2use is at least on the same beam as chanID

[dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,dum10,dum11,dum12, ...
      dum13,dum14,dum15,dum16,dum17,dum18,dum19,dum20,dum21, ...
      dum22,beamID2use,dum23,dum24,dum25,dum26,dum27] ...
    = getConst(scanOrbitPath,antDerivedPath,pointingPath,calPath,chanID2use);

if (~strcmp(beamID,beamID2use))
  error(['err[getSensorCs_SCF]: beamID ~= beamID2use:' ...
	beamID ', ' beamID2use]);
end
  
% Get list of files in antDir
fileRoot = [antDir '/AntPatt_' antType];
fileType = [fileRoot '*.dat'];
eval('dirStruct = dir(fileType);','noFiles = 1;')
nfiles = length(dirStruct);
if (nfiles == 0)
  error(['err[getSensorCs_SCF]: No antenna files. antdir=' antDir ...
	', antType = ' antType]);
end
for i=1:nfiles
  fileNames{i} = [antDir '/' dirStruct(i).name];
end

% Find filenames with beamID4ant
beamID4ant = beamID;
if (~isempty(strfind(beamID,'red')))
  beamID4ant = beamID(1:end-4);
end
i0set = strfind(fileNames,['_' beamID4ant '_']);
iantset = [];
for i=1:nfiles
  if (~isempty(i0set{i}))
    iantset = [iantset i];
  end
end
nantset = length(iantset);

if (nantset > 1) 
  % Get frequencies of each file
  for j=1:nantset
    i = iantset(j);
    istart = i0set{i}+2+length(beamID4ant);
    iend = istart-2+strfind(fileNames{i}(istart:end),'_');
    antfreq(j) = str2num(fileNames{i}(istart:iend));
  end
  [mn,imn] = min(abs(antfreq-freq2use));
  if (mn <= 0.01)
    % file is ~exact match to channel frequency
    iantfiles = iantset(imn);
    nantfiles = 1;
    antfilefreq = antfreq(imn);
  else
    % find bracketing frequencies
    nantfiles = 0;
    [srt,isrt] = sort(antfreq);
    for j=1:nantset-1
      if (freq2use <= antfreq(isrt(j+1)))
	if (freq2use > antfreq(isrt(j)))
	  iantfiles = iantset(isrt([j j+1]));
	  nantfiles = 2;
	  antfilefreq = antfreq(isrt([j j+1]));
	end
	break
      end
    end
    if (nantfiles == 0)
      error(['err[getSensorCs_SCF]: Antenna files with bracketing ' ...
	    'frequencies not found. antfreq = ' num2str(antfreq) ...
	    ', freq2use = ' num2str(freq2use)]);
    end
  end
  for i=1:nantfiles
    antfiles{i} = fileNames{iantfiles(i)};
  end
else
  antfilefreq = 0;
  nantfiles = 1;
  antfiles{1} = fileNames{iantset};
end

incj = []; xpntj = []; ypntj = []; incscale = []; jpatv = [];
xmatj = []; ymatj = []; 
if (nopatCk == 0) 
  % Get projected antenna pattern OPTIONALLY SCALED TO SPEC SIZE
  incj = max([[aview_spec*patsize 40]/(2*npnto2)]);
  xpntj = incj*[-npnto2:npnto2]; xmatj = repmat(xpntj,length(xpntj),1);
  ypntj = xpntj; ymatj = repmat(ypntj',1,length(xpntj));

  % no-scale option is now built into antDerived SCF file choice
  incscale = fscaleincj*incj*[1 1];

  [jpatv tmp] =  cmis_effant_elaz(antfiles,nantfiles, ...
      antfilefreq,freq2use,npnto2,incscale,vxint,0, ...
      beam_nadirrad,beam_beta,aview_spec,ascan_spec,alt,earthRadius);
  u = incj*trapz(incj*trapz(jpatv));
  jpatv = jpatv/u;
  [aview,ascan] = getfov(jpatv,xpntj,ypntj);
else
  aview = aview_nom;
  ascan = ascan_nom;
end

% sample distance along scan
perp2arc = earthRadius*sin(los_beta);
scancircle = perp2arc*2*pi;
vxsamp = scancircle*sample_time*rps;
npos = n_samples;
los_eia = rad2deg(los_eiarad);