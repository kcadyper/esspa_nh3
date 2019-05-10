function [avg,antpatrot]=getfovavg(subimg,dataimg,i1img,i2img,j1img,j2img, ...
    defVal,dataimgx,delx,dataimgy,dely,rotang,xmatj,ymatj,antpat,antpatrot)

% Fill dataimg with data near footprint
dataimg(i1img:i2img,j1img:j2img)= subimg;

%%% Translate the data to the footprint grid ORIGIN.
% datatrans is on a grid with has same coordinate system origin 
% as the footprint pattern but rotated relative to footprint
% pattern grid orientation.
% dataimgX+delX = orginal data grid coordinates relative to
% footprint center.  dataimgX = Features (0,0) value at grid
% center point for easy rotation.
datatrans = interp2(dataimgx+delx,dataimgy+dely, ...
    dataimg,dataimgx,dataimgy);

%%% Rotate data image.  Imrotate sets invalid values on the periphery
%%% to 0, so add delta=1.0E5 to move data away from 0 values.  
%%% Surfalt, pland, etc. may have 0 or <0 values.
dataimgrot = imrotate(datatrans+1.0E5,rotang,'bilinear','crop');
dataimgrot(find(dataimgrot == 0)) = NaN;
dataimgrot = dataimgrot-1.0E5;

%%%  Interp data image to footprint pattern grid locations
dataatant = interp2(dataimgx,dataimgy,dataimgrot,xmatj,ymatj);

% Check for bad values (NaNs) and avoid.  This should be ok, 
% since our data subset is larger than the footprint.  NaNs are coming
% from the interpolation of the rotated data subset.

igood = find(isfinite(dataatant));
% sum weights; sum=0 if igood is empty []
total_good = sum(antpat(igood));
if (total_good >= .9)
  %%% Sum weighted avg to get 1 value for the footprint
  weights = antpat(igood)/total_good;
  avg = sum(dataatant(igood).*weights);
else
  avg = defVal;
end   %  end test on percent good


    
    