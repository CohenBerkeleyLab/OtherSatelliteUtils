function [  ] = fix_gcavg_latlon( path_in, path_out )
%FIX_GCAVG_LATLON Fix GC_avg variables having NaNs for geographic coords
%   I derped. One time fix script.

F = dir(fullfile(path_in,'*.mat'));

D = load(fullfile(path_in, F(1).name));
sz = size(D.GC_avg.Longitude);
xres = 360/sz(1);
yres = 180/sz(2);
lonvec = -180+xres/2:xres:180-xres/2;
latvec = -90+yres/2:yres:90-yres/2;
[LON, LAT] = meshgrid(lonvec, latvec);
LON = LON';
LAT = LAT';

parfor f=1:numel(F)
    E = load(fullfile(path_in, F(f).name));
    GC_avg = E.GC_avg;
    GC_avg.Longitude = LON;
    GC_avg.Latitude = LAT;
    saveData(fullfile(path_out, F(f).name), GC_avg)
end


end

function saveData(filepath, GC_avg)
save(filepath, 'GC_avg');
end
