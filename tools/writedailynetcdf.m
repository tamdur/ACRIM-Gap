%Write netcdf file containing output from ACRIM Gap experiment 
% Ted Amdur
% 2023/03/20
load('ar2_23_03_17.mat')
load('obs_23_02_01.mat','offsets')
tsi=xAll;
fillvalue = 1.0e+36;
nccreate('btsi_19840615_19961016.nc','TSI_ensemble','Dimensions', {'Date',149,'realization',10000},...
         'FillValue',fillvalue);
ncwrite('btsi_19840615_19961016.nc','TSI_ensemble',tsi);
nccreate('btsi_19840615_19961016.nc','TSI','Dimensions', {'Date',149},...
         'FillValue',fillvalue);
ncwrite('btsi_19840615_19961016.nc','TSI',prctile(tsi',50));
JD=juliandate(outDat.opts.dateM);
nccreate('btsi_19840615_19961016.nc','JD','Dimensions', {'Julian_Day',149},...
         'FillValue',fillvalue);
ncwrite('btsi_19840615_19961016.nc','JD',zeros(149,1)');
nccreate('btsi_19840615_19961016.nc','Assumed_Solar_Constant','Dimensions', {'value',1},...
         'FillValue',fillvalue);
ncwrite('btsi_19840615_19961016.nc','Assumed_Solar_Constant',offsets(~outDat.oindex));

%Write global attributes
ncwriteatt('btsi_19840615_19961016.nc','/','title','Monthly TSI for the ACRIM Gap as calculated by BTSI');
ncwriteatt('btsi_19840615_19961016.nc','/','source','BTSI_v1');
ncwriteatt('btsi_19840615_19961016.nc','/','institution','Harvard University Department of Earth and Planetary Sciences');
ncwriteatt('btsi_19840615_19961016.nc','/','created_by','Ted Amdur (amdur@g.harvard.edu)');
ncwriteatt('btsi_19840615_19961016.nc','/','creation_date',datestr(now));
ncwriteatt('btsi_19840615_19961016.nc','/','DOI','');

ncwriteatt('btsi_19840615_19961016.nc','/','comments',strcat("Published here are the full ensemble of 10,000 realizations of TSI,", ...
    " as well as the maximim likelihood estimate.",...
    " TSI values are published as anomalies from the ACRIM2 mean over this time interval."));

%TSI attributes
ncwriteatt('btsi_19840615_19961016.nc','/TSI_ensemble','standard_name','TSI_anomaly');
ncwriteatt('btsi_19840615_19961016.nc','/TSI_ensemble','long_name','Total_solar_irradiance_anomaly');
ncwriteatt('btsi_19840615_19961016.nc','/TSI_ensemble','units','Watts_per_meter_squared');
ncwriteatt('btsi_19840615_19961016.nc','/TSI_ensemble','missing_value','1.0e+36');

ncwriteatt('btsi_19840615_19961016.nc','/TSI','standard_name','TSI_anomaly');
ncwriteatt('btsi_19840615_19961016.nc','/TSI','long_name','Total_solar_irradiance_anomaly');
ncwriteatt('btsi_19840615_19961016.nc','/TSI','units','Watts_per_meter_squared');
ncwriteatt('btsi_19840615_19961016.nc','/TSI','missing_value','1.0e+36');

%JD attributes
ncwriteatt('btsi_19840615_19961016.nc','/JD','standard_name','Julian_day');
ncwriteatt('btsi_19840615_19961016.nc','/JD','units','int');

%Assumed_Solar_Constant attributes
ncwriteatt('btsi_19840615_19961016.nc','/Assumed_Solar_Constant','standard_name','Solar Constant');
ncwriteatt('btsi_19840615_19961016.nc','/Assumed_Solar_Constant','units','Watts_per_meter_squared');




