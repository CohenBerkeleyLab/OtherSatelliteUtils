function [  ] = modis_truecolor_image(  )
%modis_truecolor Creates a truecolor image at 500 m resolution using MOD09
%   This function will loop through all the MOD09 files in a folder for the
%   date given and look for those that fall in the specified lon/lat
%   boundaries. It will interpolate the data to a 500 m square grid in
%   order to combine multiple granules and in order to allow the image
%   command to display lat/lon values

%%%%%%%%%%%%%%%%%%%%%%
%%%%% USER INPUT %%%%%
%%%%%%%%%%%%%%%%%%%%%%

lonlim = [-125, -114];
latlim = [32, 42];

mod_date = '06/22/2008';

mod09_dir = '/Volumes/share-sat/SAT/MODIS/MYD09/MODIS_CA_Fires_Jun2008/';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% STATIC VARIABLES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mod09_pattern = 'MYD09.A%s%03d*.hdf';
L = 0;

%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN LOOP %%%%%
%%%%%%%%%%%%%%%%%%%%%

the_date = datestr(mod_date,29);
mod_year = the_date(1:4);
mod_day = modis_date_to_day(the_date);
mod09_file = sprintf(mod09_pattern,mod_year,mod_day);

mod09_files = dir(fullfile(mod09_dir, mod09_file));
for a=1:numel(mod09_files)
    fname = mod09_files(a).name;
    modf = fullfile(mod09_dir, fname);
    fprintf('Now examining %s\n',modf);
    lon_temp = hdfread(modf,'Geolocation Fields/Longitude');
    lat_temp = hdfread(modf,'Geolocation Fields/Latitude');
    
    % If the file has no points within our boundaries of interest, skip it
    if all(lon_temp(:) < min(lonlim)) || all(lon_temp(:) > max(lonlim))
        continue
    elseif all(lat_temp(:) < min(latlim)) || all(lat_temp(:) > max(latlim))
        continue
    end
    
    L = L+1;
    % Interpolate the longitude and latitude to 500 m increments and append
    % them to the running matrix of lons/lats
    [lon,lat] = interp_mod09_latlon(lon_temp, lat_temp);
    if ~exist('lons','var')
        lons = lon(:); lats = lat(:);
    else
        lons = [lons; lon(:)]; lats = [lats; lat(:)];
    end
    
    % Read in the 500m surface reflectance bands 1, 3, and 4 (red, blue,
    % and green respectively)
    r = hdfread(modf, 'Data Fields/500m Surface Reflectance Band 1');
    g = hdfread(modf, 'Data Fields/500m Surface Reflectance Band 4');
    b = hdfread(modf, 'Data Fields/500m Surface Reflectance Band 3');
    
    % The reflectances need scaled properly to make the images the
    % appropriate brightness.  
    [r,g,b] = scale_modis(r,g,b);
    
    % Add these to the running list of color values
    if ~exist('reds','var')
        reds = r(:);
        greens = g(:);
        blues = b(:);
    else
        reds = [reds; r(:)];
        greens = [greens; g(:)];
        blues = [blues; b(:)];
    end
end

% Create grid lat/lon vectors, assuming 500 m = 0.005 degrees
grid_lonvec = min(lonlim):0.005:max(lonlim);
grid_latvec = min(latlim):0.005:max(latlim);

[latmat, lonmat] = meshgrid(grid_latvec, grid_lonvec);

% Reshape the color matrices into vectors and interpolate to the grid
% points defined by the lat/lon matrices above.

latvec = double(lats);
lonvec = double(lons);


xx = latvec > min(latlim) & latvec < max(latlim) & lonvec > min(lonlim) & lonvec < max(lonlim);


redmat = griddata(latvec(xx),lonvec(xx),reds(xx),latmat,lonmat);
greenmat = griddata(latvec(xx),lonvec(xx),greens(xx),latmat,lonmat);
bluemat = griddata(latvec(xx),lonvec(xx),blues(xx),latmat,lonmat);

% The image function displays these upside down
redmat = redmat';
greenmat = greenmat';
bluemat = bluemat';

rgb = cat(3,redmat,greenmat,bluemat);

if min(rgb(:))<0;
    rgb = rgb + min(rgb(:));
end
if max(rgb(:))>1;
    rgb = rgb ./ max(rgb(:));
end

figure; image(lonlim,latlim,rgb)
set(gca,'ydir','normal');
putvar(rgb,lonlim,latlim);

end

