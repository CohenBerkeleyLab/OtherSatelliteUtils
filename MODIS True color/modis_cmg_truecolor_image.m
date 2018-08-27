function [  ] = modis_cmg_truecolor_image( modis_file, lonlim, latlim )
%MODIS_CMG_TRUECOLOR_IMAGE Create a truecolor image from a CMG MODIS file
%   MODIS data is (or at least was) available in two formats: 1) what I'll
%   call "granule level" that provides about 5 minutes of MODIS data in
%   each file, and 2) data gridded to a fixed (usually equirectangular)
%   grid that contains data from a whole day or multiple days. This
%   latter format is what I'll call "cmg" for "climate modeling grid."
%
%   MODIS_CMG_TRUECOLOR_IMAGE( MODIS_FILE ) takes a path to a MODIS M?D09
%   CMG file and plots a truecolor image. If you want to plot granule level
%   data, see MODIS_TRUECOLOR_IMAGE(). This function does not accept
%   multiple files, since averaging an image in time is not super useful.

E = JLLErrors;

if ~exist(modis_file, 'file')
    E.badinput('%s does not exist', modis_file);
end

if ~exist('lonlim', 'var')
    lonlim = [-180 180];
end
if ~exist('latlim', 'var')
    latlim = [-90 90];
end

modi = hdfinfo(modis_file);

% Read in the 500m surface reflectance bands 1, 3, and 4 (red, blue,
% and green respectively). Currently, I think there's a bug in the product,
% because the file specification says that the parameter should =
% (file_data - add_offset) * scale_factor, but the scale factor in the
% specification is 0.0001, not 10000 like it is in the files (file space:
%   https://ladsweb.modaps.eosdis.nasa.gov/api/v1/filespec/collection=6&product=MYD09CMG)
% Here, it doesn't actually matter because we need to rescale to [0, 255]
% before applying the brightness scale
fill_val = -28672; % when removing fill values, (below) I usually give a little bit of room for error, hence the <=*0.9)
r = double(hdfread(modi.Filename, '/MOD09CMG/Data Fields/Coarse Resolution Surface Reflectance Band 1'));
r(r<=fill_val*0.9)=nan;
%r = r / 10000;
g = double(hdfread(modi.Filename, '/MOD09CMG/Data Fields/Coarse Resolution Surface Reflectance Band 4'));
g(g<=fill_val*0.9)=nan;
%g = g / 10000;
b = double(hdfread(modi.Filename, '/MOD09CMG/Data Fields/Coarse Resolution Surface Reflectance Band 3'));
b(b<=fill_val*0.9)=nan;
%b = b / 10000;

% Convert these to the 0 to 255 range
color_range = [0 255];
r = range_squeeze(r, color_range);
g = range_squeeze(g, color_range);
b = range_squeeze(b, color_range);

% The reflectances need scaled properly to make the images the
% appropriate brightness. This assumes that the RGB values are given on the
% range [0, 255].
[r,g,b] = scale_modis(r,g,b);

rgb = cat(3, r, g, b);

if min(rgb(:))<0;
    rgb = rgb + min(rgb(:));
end
if max(rgb(:))>1;
    rgb = rgb ./ max(rgb(:));
end

[lon, lat, xx, yy] = modis_cmg_latlon(0.05, lonlim, latlim);
% The image data will be upside down by default. We can subset it at the
% same time.
%rgb = flipud(rgb(yy,xx,:));
rgb = rgb(yy,xx,:);

figure;
image(lon, lat, rgb);
set(gca,'ydir','normal');
end

