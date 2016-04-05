function [ lon2, lat2 ] = interp_mod09_latlon( lon, lat )
%interp_mod09_latlon Interpolate MOD09 lat/lon to 500 m resolution
%   Geolocation data for MOD09 data sets is given only for the 1000 m
%   resolution.  To use it for 500 m resolution, the lat/lon data needs to
%   be interpolated as seen at
%   https://earthdata.nasa.gov/sites/default/files/field/document/MODIS_True_Color.pdf
%   on page 11.

narginchk(2,2);
E = JLLErrors;

% Verify that lon and lat have the same dimensions
s = size(lon);
if ~all(s==size(lat))
    E.dimMismatch('lon','lat');
end

% Assuming the 1 km resolution data falls on the grid at 1,2,...,n then in
% the y-direction (first dimension), the points need to be interpolated to
% 0.75, 1.25, 1.75,...,(n+0.25) and to 1, 1.5, 2,...,(n+0.5) along the
% x-direction (second dimension).


xqv = 1:0.5:s(2)+0.5;
yqv = 0.75:0.5:s(1)+0.25;

[xq,yq] = meshgrid(xqv,yqv);

lon2 = interp2(lon,xq,yq);
lat2 = interp2(lat,xq,yq);

% The right edge and top and bottom of the matrix requires extrapolation,
% which interp2 won't do. Do a very simple guesstimate by adding the
% closest difference to the closest point.

lon2(2:end-1,end) = (lon2(2:end-1,end-1)-lon2(2:end-1,end-2)) + lon2(2:end-1,end-1);
lat2(2:end-1,end) = (lat2(2:end-1,end-1)-lat2(2:end-1,end-2)) + lat2(2:end-1,end-1);

lon2(1,:) = (lon2(2,:) - lon2(3,:)) + lon2(2,:);
lat2(1,:) = (lat2(2,:) - lat2(3,:)) + lat2(2,:);

lon2(end,:) = (lon2(end-1,:) - lon2(end-2,:)) + lon2(end-1,:);
lat2(end,:) = (lat2(end-1,:) - lat2(end-2,:)) + lat2(end-1,:);

end

