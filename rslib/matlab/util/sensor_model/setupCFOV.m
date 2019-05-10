csrID = csrIDs{icsr};
[cfov_anglexpos,cfov_npos,cfov_xxpos,cfov_yxpos,csr,swxpos, ...
      cfov_nadirrad,sample_angle,sample1_time,sample_time, ...
      scanStep,a2sRPY,s2scRPY,cfov_tysampxpos,cfov_chans] ...
    = getcfovCs(scanOrbitPath,pointingPath,antDerivedPath, ...
    calPath,chan2antPath,antDir,antType,csrID,nominalChannel,nomSatAlt);

cfov_nom_size = csr(1);
csrname = ['csr' csrID];

if (cfovCk == 1) 
  if (rfovSizeMethod == 1)
    xid4CFOV = [csrname 'std'];
  elseif (rfovSizeMethod == 2)
    xid4CFOV = [csrname 'frcnom'];
  else
    error(spatAvg:setupCFOV:: Unexpected rfovSizeMethod.')
  end
else % squareCk == 1
  xid4CFOV = ['hcs' num2str(cfov_nom_size)];
end
if (rfovSize2Force ~= 0)
  xid4CFOV = [xid4CFOV '_size' num2str(rfovSize2Force)];
end
% Define output file names
avgGeoFile = [avgGeoRoot '_' xid4CFOV '.nc'];
avgradfile = [avgRadRoot '_' xid4CFOV '.nc'];

if (~(vertAvgFlag == 1 & vertAvgLocFlag == 1))
  % Default: Load ascii cfov orbital data file from orbsim.
  filename = ['orbitdata_' csrname '_' speedfactorstr '.txt'];
  eval(['load ' workDir '/' filename]);
  datastr = filename(1:end-4);
  eval(['orbdata = ' datastr ';']); 
else
  % Special case with geolocation coming from input grid.
  avgGeoFile = [avgGeoRoot '_grid_' xid4CFOV '.nc'];

  % Read ascii vertical profile location file
  filename = [vertAvgLocRoot csrID '.dat'];
  fid = fopen(filename,'r');
  sz = fscanf(fid,'%d',2);
  orbdata = zeros(sz(1)*sz(2),12); ipnt = 0;
  for iscan=1:sz(2)
    for ipos=1:sz(1)
      ipnt = ipnt+1;
      orbdata(ipnt,1) = iscan;
      orbdata(ipnt,7:8) = fscanf(fid,'%f',2)';
    end
  end
  fclose(fid);
  maxScans2Avg = Inf;
  npos = sz(1);
end

% Initialize tb array
nscans = max(orbdata(:,1));
nscans2avg = min(nscans,maxScans2Avg);

numprofs = npos*nscans2avg/npos_skip;
avgtbs4icsr = missingVal*ones(numprofs,nchans4icsr);
centerij_array = missingVal*ones(numprofs,2);

% Redo mapBeam for chanlist4icsr
mapBeamtmp = mapBeam(chanlist4icsr);
% Shift mapping for beams removed when set was reduced
oldset = [];
for i=1:nBeam
  if (any(mapBeamtmp == i))
    oldset = [oldset i];
  end
end
nBeam4icsr = length(oldset);
mapBeam4icsr = zeros(1,nchans4icsr);
for i=1:nBeam4icsr
  iset = find(mapBeamtmp == oldset(i));
  mapBeam4icsr(iset) = i;
end