function [eff,instant,omegaE] = cmis_effant_elaz(antfiles,nantfiles, ...
    antfilefreq,freq2use,npnto2,inc,vxinttime_in,circCk, ...
    nadirrad,beta,aview,ascan,alt,earthRadius)

  npnt = 2*npnto2+1;
  altplusR = earthRadius+alt;

  % Note that we assume the az,el = 0 location is the geometric boresight
  % of the antenna.  The direction of maximum gain for a particular
  % channel (i.e., 19H) may not be directly along this boresight 
  % as measured in the files.  

  % EL over AZ: See notes, 10/20/2000.  Az moves then el to define coordinate.
  % EL can remain constant while az changes to define each image row.
  
  % Input data
  for ifile=1:nantfiles
    fid = fopen(antfiles{ifile},'r');
    if (fid <= 0)
      error(['cmis_effant_elaz: Can''t open ' antfiles{ifile} ' file'])
    end
    tmp = fgets(fid);     tmp = fgets(fid);     tmp = fgets(fid); 
    freqLoc = fscanf(fid,'%*1s%8f',1);
    nv = fscanf(fid,'%*6s%4d',1);
    nh = fscanf(fid,'%*6s%4d',1);
    tmp = fgets(fid); tmp = fgets(fid); tmp = fgets(fid);
    col123 = [];
    for i=1:(nv*nh)
      col123(i,:) = fscanf(fid,'%8f',3);
      tmp = fgets(fid);
    end
    fclose(fid);

    % Check that freqLoc matches antfilefreq(ifile)
    if (abs(freqLoc-antfilefreq(ifile)) > 0.01) 
      error(['err[cmis_effant_elaz]: Antenna pattern file frequencies do ' ...
	    ' not match.  freqLoc, antfilefreq, ifile: ' ...
	    num2str([freqLoc antfilefreq(ifile) ifile])]);
    end
    
    gainant = 10.^(reshape(col123(:,3),nv,nh)'/10); 
    gainant = gainant/max(gainant(:));    
    if (ifile == 1)
      el2ant = deg2rad(reshape(col123(:,1),nv,nh))'; 
      az2ant = deg2rad(reshape(col123(:,2),nv,nh))';
      gainant1 = gainant;
    else
      gainant2 = gainant;
    end
  end
    
  if (nantfiles == 1)
    gainant = gainant1;
  else
    % interpolate gain in frequency
    gainant = gainant1 + (freq2use-antfilefreq(1))*(gainant2-gainant1) ...
	/(antfilefreq(2)-antfilefreq(1));
  end

  % OmegaE is the integral of normalized gain over earth solid angle
  % Differential solid angle in ELoverAz coordinates is cosEldEldAz
  
  dEldAz = (el2ant(2,1)-el2ant(1,1))*(az2ant(1,2)-az2ant(1,1));
  omegaE = dEldAz*trapz(trapz(gainant.*cos(el2ant)));
  
  % Set up
  xpnt = inc(1)*[-npnto2:npnto2]; ypnt = inc(2)*[-npnto2:npnto2];
  xmat = repmat(xpnt,npnt,1);ymat = repmat(ypnt',1,npnt);

  if (1 == 1)
    % Get az2, el2 of footprint points
    el2lim = deg2rad(10.0); el2ref = [-el2lim:0.0001:el2lim];
    nadirref = nadirrad+el2ref;
    eiaref = asin(altplusR*sin(nadirref)/earthRadius);
    betaref = eiaref - nadirref;
    ypref = (betaref - beta)*earthRadius;
    sat2scanref = earthRadius*sin(betaref)./sin(nadirref);

    % Elevation is constant along rows and with y.
    % x,y earth grid is defined as move y first then x.
    el2 = repmat(interp1(ypref,el2ref,ypnt)',1,npnt);

    xtan = earthRadius*tan(xmat/earthRadius);
    sat2scan_ymat = repmat(interp1(ypref,sat2scanref,ypnt)',1,npnt);
    az2 = asin(xtan./sat2scan_ymat);

  else %another approximation, not as accurate?
    % Get az,el of footprint points
    ymat_beta = beta+ymat/earthRadius;
    ymat_nadir2scan = earthRadius*sin(ymat_beta);
    ymat_nadirrad = atan(ymat_nadir2scan./(a-earthRadius*cos(ymat_beta)));
    ymat_sat2scan = ymat_nadir2scan./sin(ymat_nadirrad);
    ymat_eiarad = ymat_beta + ymat_nadirrad;
    eiarad = beta+nadirrad
    el = rad2deg(atan(sqrt((ymat.*cos(ymat_eiarad)).^2+xmat.^2) ...
		     ./ymat_sat2scan));  
	% -> curved earth approximation
    az = rad2deg(-atan2(xmat,ymat.*cos(ymat_eiarad)));
    az(find(az < 0)) = 360.0+az(find(az < 0));
  end

  % Interpolate gain for instantaneous pattern
  instant = interp2(az2ant,el2ant,gainant,az2,el2,'linear',0.0);
  instant(find(isnan(instant))) = 0;

  % Convert to gain per dA (was dOmega or solid angle)
  coseia = cos(asin(altplusR*sin(el2+nadirrad)/earthRadius));
  instant = instant.*coseia./(sat2scan_ymat.^2);

  % Convolve instant pattern for EFOV

  if (circCk < 0)
    % Just return instant
    eff = NaN*size(instant);
  else
    %Define the smoothing function
    if (circCk == 0)
      ivxinttime = round(vxinttime_in/inc(1));
    elseif (circCk == 1)
      ivxinttime = round((aview-(aview-ascan)/10)/inc(1));
    else
      % need different settings depending on scalechan pattern
      %    ivxinttime =round((aview*circCk-(aview-ascan) ...
      %                /max([(22.1*circCk-13.1) 1]))/inc(1));
      ivxinttime =round((aview*circCk-(aview-ascan) ...
	  /max([(22.1*circCk-10) 1]))/inc(1));
    end
    rem = npnt-ivxinttime;
    intbox = [zeros(1,round(rem/2)) ones(1,ivxinttime) ...
	  zeros(1,npnt-ivxinttime-round(rem/2))];
    eff = conv2(1,intbox,instant,'same');
  end

  return