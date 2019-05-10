function [refpat,incref,xmatref,ymatref] ...
    = getrefpatSM(refpatFileRoot,rfovScaleChanDef,nomcsr, ...
    scanOrbitPath,pointingPath,antDerivedPath,calPath,chan2antPath, ...
    antDir,antType,nomSatAlt,chanID,npnto2,patsize,satAlt, ...
    rfovSize2Force,rfovSizeMethod,squareCk)
% For calling from Testbed spatavg.m etc.
% Limited functionality compared to prior versions.  See genrefpat.,m
% for off-line calculations of xAndyint and incScale parameters. Also has
% simple square (tophat) pattern option.
% rev. 20060301: Now rfovSizeMethod=1 only uses EFOV size not
% spec to generate RFOV, although a spec-based version may be
% run by changing the defaults to input the spec SCFs.

% Use satellite altitute for output reference pattern.
[beam_scananglexpos,beam_nadirradj, ...
      jpatv,incj,xpntj,xmatj,ypntj,ymatj, ...
      tysamp,xzxpos,yzxpos,xlosxpos,ylosxpos,los_scananglexpos, ...
      aviewj_spec,ascanj_spec,vxsampj,vxintj,nposj,design,aviewj,ascanj,...
      earthRadius,beam_betaj,nedtj,k_calj,scan_start_timej, ...
      sample_timej,rps,integ_timej,los_eiaj,beamIDj,iBeamj, ...
      a2sRPY,s2scRPY,tysampxpos,dyintegxpos, ...
      antfiles,nantfiles,antfilefreq,freq2use] ...
    = getSensorCs_SCF(scanOrbitPath,pointingPath,antDerivedPath, ...
    calPath,chan2antPath,antDir,antType,nomSatAlt,chanID,npnto2, ...
    0,patsize);

hsrj = [aviewj ascanj];
hsrj_spec = [aviewj_spec ascanj_spec];
% 12/23/05: Was refhsr = max(1.01*hsrj_spec,nomcsr)
% Now rounds to 1 decimal place and uses max of spec and nom hsr
% since some nom hsr may be larger than spec!
hsrjmax = 0.1*ceil(1.01*max(hsrj_spec,hsrj)*10);
% 3/1/06: Always use hsrj as noted at top regardless of spec.
hsrjmax = 0.1*ceil(1.01*hsrj*10);

% Default rfovScaleChan
rfovScaleChan = rfovScaleChanDef;

% Determine refhsr and perhaps change rfovScaleChan
if (rfovSize2Force > 0)
  refhsr = rfovSize2Force*[1 1];
elseif (rfovSizeMethod == 2 | squareCk == 1)
  % Force nominal CFOV size
  refhsr = nomcsr;
elseif (rfovSizeMethod == 1)
  % Refhsr is always at least a little bigger than sensor channel hsr.
  refhsr = max(hsrjmax,nomcsr);
  if (any(refhsr > nomcsr) | (nomcsr == 15))
    rfovScaleChan = chanID;
  end
else
  error('getrefpatSM: rfovSizeMethod not supported.')
end

incref = max([[refhsr*patsize 40]/(2*npnto2)]);
xpntref = incref*[-npnto2:npnto2];  
xmatref = repmat(xpntref,length(xpntref),1);
ypntref = xpntref; 
ymatref = repmat(ypntref',1,length(xpntref));

if (squareCk == 0)

  % Use scalechan (SC suffix) as the basis for fabricating refpat. 
  % Only need to get nadir and beta angles and antenna file name
  [beam_scananglexposSC,beam_nadirradSC, ...
	jpatvSC,incSC,xpntSC,xmatSC,ypntSC,ymatSC,tysamp, ...
	xzxposSC,yzxposSC,xlosxposSC,ylosxposSC,los_scananglexposSC, ...
	aviewSC_spec,ascanSC_spec,vxsampSC,vxintSC,nposSC,design, ...
	aviewSC,ascanSC,...
	earthRadius,beam_betaSC,nedtSC,k_calSC,scan_start_timeSC, ...
	sample_timeSC,rps,integ_timeSC,los_eiaSC,beamIDSC,iBeamSC, ...
	a2sRPY,s2scRPY,tysampxposSC,dyintegxposSC, ...
	antfilesSC,nantfilesSC,antfilefreqSC,freq2useSC] ...
      = getSensorCs_SCF(scanOrbitPath,pointingPath,antDerivedPath, ...
      calPath,chan2antPath,antDir,antType,nomSatAlt,rfovScaleChan, ...
      npnto2,0,patsize);

  % Find case in input file
  fid = fopen([refpatFileRoot design '_' num2str(rfovSizeMethod) ...
	'.txt'],'r');
  irow = [];
  if (fid ~= -1)
    count = 1; refhsrIn = []; scalePrmIn = []; scaleChanIn = {}; irow = 0;
    while (count > 0)
      [rdum1,count] = fscanf(fid,'%e',2);
      [rdum2,count] = fscanf(fid,'%e',3);
      [adum,count] = fscanf(fid,'%s',1);
      if (count > 0)
	irow = irow+1;
	refhsrIn = [refhsrIn; rdum1'];
	scalePrmIn = [scalePrmIn; rdum2'];
	scaleChanIn{irow} = adum;
      end
    end
    fclose(fid);
    scaleChanCk = zeros(irow,1);
    imatch = strmatch(rfovScaleChan,scaleChanIn);
    scaleChanCk(imatch) = 1;
    irow = find(abs(refhsrIn(:,1) - refhsr(1)) <= 0.2 ...
	& abs(refhsrIn(:,2) - refhsr(2)) <= 0.2 & scaleChanCk);
    if (length(irow) > 1)
      %error('getrefpatSM: multiple matches for refhsr/scalechan')
      [mx,imx] = max(refhsrIn(irow,1) - refhsr(1));
      irow = irow(imx);
    end
    if (isempty(irow))
      error(['getrefpatSM: no match for refhsr/scalechan: ' ...
	    num2str(refhsr) ' ' rfovScaleChan])
    end
  else
    error(['getrefpatSM: ' [refpatFileRoot design '_' ...
	    num2str(rfovSizeMethod) '.txt'] ' missing']);
  end
  
  % Using archived scaling data
  xAndyint = scalePrmIn(irow,1:2);
  incScale = scalePrmIn(irow,3);

  % Check reference footprint pattern at nomSatAlt
  
  refpat = cmis_refpat_elaz(antfilesSC,nantfilesSC, ...
      antfilefreqSC,freq2useSC,npnto2,incScale*incref, ...
      xAndyint(1),xAndyint(2),beam_nadirradSC,beam_betaSC,nomSatAlt, ...
      earthRadius);

  [av,as] = getfov(refpat,xpntref,ypntref);
  % less precise if npnto2 is smaller or incref is larger than defaults
  incref0 = max([[refhsr*2.0 40]/(2*200)]);
  error0 = 0.3; if (max(refhsr) >= 50) error0 = 0.4; end
  dFOVgoal = error0*max(200/npnto2,1)*refhsr(1)/40.0 + max(incref-incref0,0); 
  if ((abs(av-refhsr(1)) > dFOVgoal | abs(as-refhsr(2)) > dFOVgoal))
    error(['getrefpatSM: Refpat HSR failed to match reference goal. ' ...
	  'dFOVgoal, refhsr, av, as: ' num2str([dFOVgoal refhsr av as])])
  end

  if (satAlt ~= nomSatAlt)
    refpat = cmis_refpat_elaz(antfilesSC,nantfilesSC, ...
	antfilefreqSC,freq2useSC,npnto2,incScale*incref, ...
	xAndyint(1),xAndyint(2),beam_nadirradSC,beam_betaSC,satAlt, ...
	earthRadius);
  end
    
  % Center refpat
  [x0,y0] = getxymean(refpat,xmatref,ymatref);
  refpat = interp2(xmatref-x0,ymatref-y0,refpat,xmatref,ymatref, ...
      'bilinear',0.0);
  refpat(find(isnan(refpat))) = 0;
    
else
  % Get square pattern. Doesn't depend on satellite altitude
  refpat = (abs(xmatref) < refhsr(2)/2) & (abs(ymatref) < refhsr(1)/2);
end
