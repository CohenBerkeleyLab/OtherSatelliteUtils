function [ band_3_alb ] = get_mcd43_albedo( lon, lat, date_in )
%GET_MCD43_ALBEDO Get an albedo for a given location and date
%   Searches for the MODIS MCD43C3 albedo and returns the band 3 value.
%   Meant for single locations, to be used interactively.

E = JLLErrors;

satshare = '/Volumes/share-sat';
modis_path = fullfile(satshare, 'SAT','MODIS','MCD43C3');

[modis_lons, modis_lats] = mcd43_lon_lats;


% Figure out which file to load - remembering that these are produced every
% 8 days.

date_n = datenum(date_in);

year_in = year(date_n);
jul_day = modis_date_to_day(date_in);

% Go to the appropriate folder and make a list of the day-of-year of each
% file. We'll use that to calculate the closest one.

folder = fullfile(modis_path,num2str(year_in));
modis_pattern = sprintf('MCD43C3.A%04d*.hdf',year_in);
F = dir(fullfile(folder, modis_pattern));
files_doy = zeros(size(F));
for a=1:numel(F)
    files_doy(a) = str2double(F(a).name(14:16));
end

delta = abs(files_doy - jul_day);
% If two days are equally close, this will pick the earlier one.
[~,closest_day] = min(delta(:));

% Load the band 3 albedo (wavelengths that matter for NO2)
mcd43_info = hdfinfo(fullfile(folder, F(closest_day).name));
band3 = hdfread(hdf_dsetID(mcd43_info,1,1,'Albedo_BSA_Band3'));
band3 = double(band3);
band3 = flipud(band3);
band3(band3==32767)=NaN; % 32767 is the fill value for this data set
band3 = band3 * 1e-3; % Albedo matrix needs to have the scale factor applied

% Double check that band3 has the same shape as the lon/lat matrices
if any(size(band3) ~= size(modis_lons)) || any(size(band3) ~= size(modis_lats))
    E.sizeMismatch('modis_lons', 'modis_lats', 'band3');
end

% Interpolate the albedo to the given lat/lon. Check that it stays in the
% expected range of 0 to 1.

band_3_alb = interp2(modis_lons, modis_lats, band3, lon, lat);
if band_3_alb < 0 || band_3_alb > 1
    E.badvar('band_3_alb','must be between 0 and 1');
end

end

