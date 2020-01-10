clear
clc

ncdisp('file_precip.nc');        % tells you what fields are in each .nc file
ncdisp('file_sat.nc');

% precip has fields LONT, LATT, and PRECIP
% load each into separate arrays
% - you can give them the same name or not
% - I usually give them the same name, but both files have LONT and LATT
% matlab array = ncread('filename.nc','field name')

lonp = ncread('file_precip.nc','LONT');
latp = ncread('file_precip.nc','LATT');
precip = ncread('file_precip.nc','PRECIP');

% sat has fields LONT, LATT, and SAT

lons = ncread('file_sat.nc','LONT');
lats = ncread('file_sat.nc','LATT');
sat = ncread('file_sat.nc','SAT');

% the end


