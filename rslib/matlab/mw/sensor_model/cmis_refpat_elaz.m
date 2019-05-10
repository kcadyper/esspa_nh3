function eff=cmis_refpat_elaz(antfiles,nantfiles, ...
    antfilefreq,freq2use,npnto2,inc, ...
    vxint,vyint,nadirrad,beta,alt,earthRadius)

  % Get just instant pattern for scalechan
  [dum,instant] = cmis_effant_elaz(antfiles,nantfiles, ...
      antfilefreq,freq2use,npnto2,inc*[1 1], ...
      0,-1,nadirrad,beta,0,0,alt,earthRadius);

  % Set up
  npnt = 2*npnto2+1;
  xpnt = inc*[-npnto2:npnto2]; ypnt = inc*[-npnto2:npnto2];
  xmat = repmat(xpnt,npnt,1);ymat = repmat(ypnt',1,npnt);

  % Convolves along-scan (x) and cross-scan (y) to match 
  % reference pattern to desired size

  nincsample = round(vxint/inc);
  rem = npnt-nincsample;
  intboxrow = [zeros(1,round(rem/2)) ones(1,nincsample) ...
	zeros(1,npnt-nincsample-round(rem/2))];
  if (vyint > 0)
    nincsample = round(vyint/inc);
    rem = npnt-nincsample;
    intboxcol = [zeros(1,round(rem/2)) ones(1,nincsample) ...
	  zeros(1,npnt-nincsample-round(rem/2))];
  else
    intboxcol = 1;
  end
  % convolve instant in row (x) direction and column (y) direction (if needed)
  eff = conv2(intboxcol,intboxrow,instant,'same');
  
  % Conv2 changes peak of eff so renormalize to peak=1
  eff = eff/max(eff(:));

return