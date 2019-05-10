function avg=getfovavg2(subcube,dataimg,i1img,i2img,j1img,j2img, ...
    defavg,dataimgx,delx,dataimgy,dely,rotang,xmatj,ymatj,antpat)

%%% Translate the data to the footprint grid ORIGIN.
% datatrans is on a grid with has same coordinate system origin 
% as the footprint pattern but rotated relative to footprint
% pattern grid orientation.
% dataimgX+delX = orginal data grid coordinates relative to
% footprint center.  dataimgX = Features (0,0) value at grid
% center point for easy rotation.

%%%  Rotate antpat on its grid
antpatrot = imrotate(antpat,-rotang,'bilinear','crop');
antpatrot(find(antpatrot == 0)) = NaN;

%%% Translate antpat to center on data grid
%%% AND interp antpat to data grid locations
%antpattrans = interp2(xmatj-delx,ymatj-dely,antpatrot,xmatj,ymatj);
%antatdata = interp2(xmatj,ymatj,antpatrot,dataimgx,dataimgy);

antatdata = interp2(xmatj-delx,ymatj-dely,antpatrot,dataimgx,dataimgy);
% Renomalize antenna pattern to sum to 1
antCk = isfinite(antatdata);
antatdata = antatdata/sum(antatdata(find(antCk)));

% Check for bad values (NaNs) and fill with 0s.  This should be ok, 
% since our data subset is larger than the footprint.  NaNs are coming
% from the interpolation of the rotated data subset.

nvar = size(subcube,3);
avg = defavg*ones(1,nvar);
for ivar=1:nvar
  % Fill dataimg with data near footprint
  dataimg(i1img:i2img,j1img:j2img)= subcube(:,:,ivar);

  igood = find(isfinite(dataimg) & antCk);
  % sum weights; sum=0 if igood is empty
  total_good = sum(antatdata(igood));
  if (total_good >= .9)
    %%% Sum weighted avg to get 1 value for the footprint
    weights = antatdata(igood)/total_good;
    avg(ivar) = sum(dataimg(igood).*weights);
  end   %  end test on percent good
end
