% Core of orbsim shared by EFOV and CFOV calls

% Set up orbit time loop
otime0 = eq2start_time;
if (ascCk == 1)
  flat = -10000.0;
  slat = flat;
else
  flat = 10000.0;
  slat = flat;
end
iscans = 1;

ofid = fopen(filename,'w');

% Loop thru scans
while ((flat < maxlat & ascCk == 1) | (slat > minlat & ascCk == 0))   
  % Loop thru samples in each scan 
  for isamp=1:npos
    % Calculate scan angle and time for each sample 
    sscan = sample1_angle + (isamp-1)*sample_angle;
    otime= sample1_time + otime0 + (isamp-1)*sample_time;
    
    %  pass scone, sscan and otime to fortran to get lat/lon on ground
    
    [flat,flon,slat,slon,zenang,azang] = ...
	getLatLon(scone,sscan,otime,nomSatAlt,centerAngleLim, ...
	orbitIncl,ecross,minlat,maxlat,minlon,maxlon,a2sRPY,s2scRPY);
    fprintf(ofid,['%4i %4i %6.2f %4i %6.2f %8.4f %8.4f %9.4f ' ...
	  '%6.2f %6.2f %6.2f %6.2f %6.2f\n'], ...
	iscans,ichan,scone,isamp,sscan,otime,flat,flon, ...
	slat,slon,zenang,azang,otime0);
    
  end     %  end sample loop isamp
  otime0 = otime0 + scanStep/rps;
  iscans = iscans + scanStep;
end         %  end time loop flat < maxlat (ascending) or flat > minlat
fclose(ofid);
