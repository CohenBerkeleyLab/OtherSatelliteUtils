function [  ] = load_and_grid_OMNO2( start_date, end_date )
%LOAD_AND_GRID_OMNO2 Read in and grid level 2 OMNO2 data.
%   Current implementation will simply assume that you wish to load
%   worldwide data.  This will simplify the coding slightly, but will need
%   to be modified later if you want to only load e.g. US data (which would
%   be important if you want to use a finer gridding resolution).

E = JLLErrors;
DEBUG_LEVEL = 2;

global onCluster
if isempty(onCluster)
    onCluster = false;
end

shareroot = getenv('SYNOMNT');

if onCluster
    omno2_path = '/global/home/users/laughner/myscratch/SAT/OMI/OMNO2';
    save_path = '/global/home/users/laughner/myscratch/MATLAB/Data/OMI/OMNO2/2.5x2.0-avg-newweight-allvcd';
elseif ~isempty(shareroot)
    omno2_path = fullfile(shareroot,'share-sat','SAT','OMI','OMNO2');
    save_path = fullfile(shareroot,'share2','USERS','LaughnerJ','DOMINO-OMNO2_comparison','OMNO2','2.5x2.0-avg-test');
else
    omno2_path = '/Volumes/share-sat/SAT/OMI/OMNO2';
    save_path = '/Volumes/share2/USERS/LaughnerJ/DOMINO-OMNO2_comparision/OMNO2/0.25x0.25';
end

% Existing files must be more recent than this to be considered complete
min_date = today;

% These should be nx2 cell arrays; the first column lists the dataset name
% in the HDFv5 file, the second what you want it to be called in the Data
% structure. Make the second column an empty string if you want it to be
% the same as the first. data_fields will be read from the 'Data Fields'
% group, and geo_fields from the "Geolocation Fields" group.
data_fields = { 'AmfTrop','';...
                'ScatteringWeight','';...
                'CloudFraction','';...
                'CloudPressure','';...
                'SlantColumnAmountNO2','';...
                'SlantColumnAmountNO2Destriped','';...
                'SlantColumnAmountNO2Std','';...
                'TerrainReflectivity','';...
                'TerrainHeight','';...
                'VcdQualityFlags','';...
                'XTrackQualityFlags','';...
                'ColumnAmountNO2Trop','';...
                'ColumnAmountNO2TropStd','';...
                'ColumnAmountNO2Strat','';...
                'ColumnAmountNO2StratStd',''};
geo_fields = {  'Longitude','';...
                'FoV75CornerLongitude','Loncorn';...
                'Latitude','';...
                'FoV75CornerLatitude','Latcorn';...
                'GroundPixelQualityFlags',''};
            
for a=1:size(data_fields,1)
    if isempty(data_fields{a,2})
        data_fields{a,2} = data_fields{a,1};
    end
end

for a=1:size(geo_fields)
    if isempty(geo_fields{a,2})
        geo_fields{a,2} = geo_fields{a,1};
    end
end

all_struct_fields = cat(1,geo_fields(:,2),data_fields(:,2));

xres = 2.5; yres = 2.0;
%[loncorn, latcorn] = grid_corners(xres, yres);
[loncorn, latcorn] = geos_chem_corners;
loncorn = loncorn';
lon = loncorn(1:end-1,1:end-1) + xres/2;
loncorn(end,:) = 181;
latcorn = latcorn';
lat = latcorn(1:end-1,1:end-1) + yres/2;
latcorn(:,end) = 91;

parfor d=datenum(start_date):datenum(end_date)
    t = getCurrentTask();
    if isempty(t)
        t.ID = 0;
    end
    if DEBUG_LEVEL > 0; fprintf('W%d: Loading/gridding data for %s\n',t.ID,datestr(d)); end

    % Skip this day if the file has already been created and is new enough
    save_name = sprintf('OMI_OMNO2_%04d%02d%02d.mat',year(d),month(d),day(d));
    if exist(fullfile(save_path, save_name),'file')
        F = dir(fullfile(save_path, save_name));
        if F.datenum > datenum(min_date)
            fprintf('W%d: Skipping %s, already complete\n',t.ID,save_name);
            continue
        end
    end

    if DEBUG_LEVEL > 2; fprintf('W%d: Findings files for %s\n',t.ID,datestr(d)); end
    % OMNO2 data should be organized in subfolders by year and month
    full_path = fullfile(omno2_path, sprintf('%04d',year(d)), sprintf('%02d',month(d)));
    file_pat = sprintf('OMI-Aura_L2-OMNO2_%04dm%02d%02d*.he5',year(d),month(d),day(d));
    F = dir(fullfile(full_path, file_pat));
    Data = make_empty_struct_from_cell(all_struct_fields);
    Data = repmat(Data,1,numel(F));
    
    if DEBUG_LEVEL > 2; fprintf('W%d: Initializing GC structure\n', t.ID); end
    omi_mat = nan(size(loncorn)-1);
    omi_cell = cell(size(loncorn)-1);
    GC = struct('AmfTrop', omi_mat, 'CloudFraction', omi_mat, 'CloudPressure', omi_mat, 'SlantColumnAmountNO2', omi_mat,...
        'SlantColumnAmountNO2Destriped', omi_mat, 'SlantColumnAmountNO2Std', omi_mat, 'TerrainReflectivity', omi_mat,...
        'TerrainHeight', omi_mat, 'ColumnAmountNO2Trop', omi_mat, 'ColumnAmountNO2TropStd', omi_mat,...
        'ColumnAmountNO2Strat', omi_mat, 'ColumnAmountNO2StratStd', omi_mat, 'Areaweight', omi_mat,...
         'VcdQualityFlags', {omi_cell}, 'XTrackQualityFlags', {omi_cell}, 'GroundPixelQualityFlags', {omi_cell});
    GC = repmat(GC,1,numel(F));
    
    %lonmin = -180;  lonmax = 180;
    %latmin = -90;   latmax = 90;
    %reslat = 2; reslon = 2.5;
    DB = struct('a', [], 'swath', [], 'lon',[],'lat',[],'aw',[],'Q',[],'W',[]); % debugging only
    for s=1:numel(F)
        if DEBUG_LEVEL > 1; fprintf('\t W%d: Handling file %d of %d\n',t.ID,s,numel(F)); end
        hi = h5info(fullfile(full_path, F(s).name));
        Data(s).Filename = F(s).name;
        Data(s) = read_fields(Data(s), hi, 2, geo_fields);
        Data(s) = read_fields(Data(s), hi, 1, data_fields);
        %OMI(s) = add2grid_DOMINO(Data(s),OMI(s),reslat,reslon,[lonmin, lonmax],[latmin, latmax]);
        [GC(s), DB] = grid_to_gc(Data(s), GC(s),loncorn,latcorn,DEBUG_LEVEL,DB,s);
    end

    % Average over the days' files to get a single day's data
    if DEBUG_LEVEL > 2; fprintf('W%d: Averaging GC structure\n',t.ID); end
    GC_avg = struct(GC(1));
    GC_avg = rmfield(GC_avg,'Areaweight');
    fns = fieldnames(GC_avg);
    aw = cat(3,GC.Areaweight);
    if DEBUG_LEVEL > 3; fprintf('W%d: Percent nans areaweight = %.4f\n', t.ID, sum(isnan(aw(:)))/numel(aw)*100); end
    total_aw = nansum2(aw,3);
    if DEBUG_LEVEL > 3; fprintf('W%d: Percent nans total areaweight = %.4f\n', t.ID, sum(isnan(total_aw(:)))/numel(total_aw)*100); end
    for f=1:numel(fns)
        if DEBUG_LEVEL > 3; fprintf('\tW%d: Avg. field %s\n',t.ID,fns{f}); end
        if ~iscell(GC_avg.(fns{f}))
            GC_avg.(fns{f}) = nansum2(cat(3, GC.(fns{f})) .* aw, 3) ./ total_aw;
        else
            GC_avg = rmfield(GC_avg, fns{f});
        end
    end

    % Add these now so that they are not changed into NaNs by multiplying with areaweight,
    % or are attempted to be gridded.
    GC_avg.Longitude = lon;
    GC_avg.Latitude = lat;
    GC_avg.Weight = total_aw;

    save_name = sprintf('OMI_OMNO2_%04d%02d%02d.mat',year(d),month(d),day(d));
    if DEBUG_LEVEL > 0; fprintf('W%d: Saving file %s\n', t.ID, fullfile(save_path, save_name)); end
    saveData(fullfile(save_path,save_name),GC_avg);

end


end

function saveData(filename,GC_avg)
    save(filename,'GC_avg','-v7.3')
end

function [loncorn, latcorn] = grid_corners(lonres, latres)
    loncorn = -180:lonres:180;
    latcorn = -90:latres:90;
    % Add a little give on the last corner since we test >= the lower index
    % corner but < the higher index one; this will prevent issues with pixels
    % at exactly 180 E or 90 N
    loncorn(end)=180.1;
    latcorn(end)=90.1;
    [loncorn, latcorn] = meshgrid(loncorn, latcorn);
end

function [GC, DB] = grid_to_gc(Data, GC, gloncorn, glatcorn, DEBUG_LEVEL, DB, swathnum)
%[gloncorn, glatcorn] = geos_chem_corners();

sz = size(gloncorn)-1;

for p=1:numel(Data.Longitude)
    [Data.Loncorn(:,p), Data.Latcorn(:,p)] = uncross_corners(Data.Loncorn(:,p), Data.Latcorn(:,p));
end

Data = calc_pixel_area(Data);
TotalAreaweight = zeros(size(GC.ColumnAmountNO2Trop));

% Filter for clouds, row anomaly and albedo (DOMINO recommends only using
% albedo < 0.3 as snow/ice can interfere in cloud retrieval)
row_anom = Data.XTrackQualityFlags > 0;
vcd_qual = Data.VcdQualityFlags > 0;
clds = Data.CloudFraction > 0.3;
alb = Data.TerrainReflectivity > 0.3;
%negvcds = Data.ColumnAmountNO2Trop < 0;
%extreme_vcd = Data.ColumnAmountNO2Trop > 1e17;

fns = fieldnames(GC);
for c=1:numel(fns)
    if DEBUG_LEVEL > 2; fprintf('File %s: Data.%s is %.4f %% nans before fill removal\n',Data.Filename(19:32), fns{c}, sum(isnan(Data.(fns{c})(:)))/numel(Data.(fns{c}))); end
    Data.(fns{c})(row_anom | vcd_qual | clds | alb) = nan;
    if DEBUG_LEVEL > 2; fprintf('File %s: Data.%s is %.4f %% nans after fill removal\n',Data.Filename(19:32), fns{c}, sum(isnan(Data.(fns{c})(:)))/numel(Data.(fns{c}))); end
end

glon_vec = gloncorn(:,1);
glat_vec = glatcorn(1,:);

earth_ellip = referenceEllipsoid('wgs84','kilometer');
gc_area = nan(numel(glon_vec)-1, numel(glat_vec)-1);
for i=1:(numel(glon_vec)-1)
    for j=1:(numel(glat_vec)-1)
        xall = [glon_vec(i), glon_vec(i), glon_vec(i+1), glon_vec(i+1), glon_vec(i)];
        yall = [glat_vec(j), glat_vec(j+1), glat_vec(j+1), glat_vec(j), glat_vec(j)];
        [xall, yall] = poly2cw(xall, yall);
        gc_area(i,j) = areaint(yall, xall, earth_ellip);
    end
end

fprintf('Starting loop\n');
tval = tic;
for a=1:numel(Data.Longitude)
    %if mod(a,100) == 1
     %   fprintf('Pixel %d: time elapsed = %f\n',a,toc(tval));
    %end
    [ax, ay] = ind2sub(size(Data.Longitude), a);
    if Data.Longitude(a) < -200 || Data.Latitude(a) < -200 || isnan(Data.Longitude(a)) || isnan(Data.Latitude(a))
        % Fill values are ~1x10^-30, this will skip any pixels with
        % fill values for coordinates
        if DEBUG_LEVEL > 3; fprintf('Skipping pixel in %s b/c lat/lon is fill\n',Data.Filename); end
        continue
    end
    xx = Data.Longitude(a) >= glon_vec(1:end-1) & Data.Longitude(a) < glon_vec(2:end);
    yy = Data.Latitude(a) >= glat_vec(1:end-1) & Data.Latitude(a) < glat_vec(2:end);
    if sum(xx) > 1 || sum(yy) > 1
        error('load_and_grid_OMNO2:pixel_assignment_error', 'Pixel %d in %s fell into more than 1 grid cell', a, Data.Filename);
    elseif sum(xx) < 1 || sum(yy) < 1
        error('load_and_grid_OMNO2:pixel_assignment_error', 'Pixel %d in %s could not be assigned to a grid cell', a, Data.Filename);
    end
    xxf = find(xx); yyf = find(yy);
    Q = calc_pix_grid_overlap(glon_vec, glat_vec, xxf, yyf, gc_area(xx,yy), squeeze(Data.Loncorn(:,ax,ay)), squeeze(Data.Latcorn(:,ax,ay)), earth_ellip);
    Weight = Data.Areaweight(a) * Q;
    TotalAreaweight(xx,yy) = nansum2([TotalAreaweight(xx,yy), Weight]);    
    for c=1:numel(fns)
        if iscell(GC.(fns{c}))
           GC.(fns{c}){xx,yy} = cat(1,GC.(fns{c}){xx,yy},Data.(fns{c})(a));
        else
            % The areaweight will be divided out at the end of the loop
            GC.(fns{c})(xx,yy) = nansum2([GC.(fns{c})(xx,yy), Data.(fns{c})(a) .* Weight]);
        end 
    end
    % Debugging only
    if xxf == 5 && yyf == 17
        DB.a = [DB.a, a];
        DB.swath = [DB.swath, swathnum];
        DB.lon = [DB.lon, Data.Longitude(a)];
        DB.lat = [DB.lat, Data.Latitude(a)];
        DB.aw = [DB.aw, Data.Areaweight(a)];
        DB.Q = [DB.Q, Q];
        DB.W = [DB.W, Weight];
    end
end

for c=1:numel(fns)
    if ~iscell(GC.(fns{c}))
        GC.(fns{c}) = GC.(fns{c}) ./ TotalAreaweight;
    end
end
GC.Areaweight = TotalAreaweight;
fprintf('\t Time for one file = %f\n',toc(tval));

%for a=1:sz(1)
%    if mod(a,25) == 1
%        fprintf('Now on %d of %d\n',a,sz(1));
%    end
%    tval_a = tic;
%    for b=1:sz(2)
%        tval=tic;
%        xx = Data.Longitude >= gloncorn(a,b) & Data.Longitude < gloncorn(a+1,b+1);
%        yy = Data.Latitude >= glatcorn(a,b) & Data.Latitude < glatcorn(a+1,b+1);
%        fprintf('\t Timer: xx yy = %f\n',toc(tval));
%        for c=1:numel(fns)
%            if sum(xx(:)&yy(:)) > 0
%                if iscell(GC.(fns{c}))
%                    tval=tic;
%                    GC.(fns{c}){a,b} = Data.(fns{c})(xx&yy);
%                    fprintf('\t Timer: binning cell arrays = %f\n',toc(tval));
%                else
%                    % Calculate an area-weighted mean
%                    tval=tic;
%                    GC.(fns{c})(a,b) = nansum2(Data.(fns{c})(xx&yy) .* Data.Areaweight(xx&yy)) ./ nansum2(Data.Areaweight(xx&yy));
%                    fprintf('\t Timer: averaging = %f\n',toc(tval));
%                end
%            end
%        end
%    end
%    fprintf('  Timer: one loop over b = %f\n',toc(tval_a));
%end

end

function Q = calc_pix_grid_overlap(glon_vec, glat_vec, xx, yy, gc_area, pixloncorn, pixlatcorn, earth_ellip)
% get the grid cell corners
t = getCurrentTask;
if isempty(t)
    t.ID = 0;
end
if any(isnan(pixloncorn)) || any(isnan(pixlatcorn)) || any(pixloncorn < -180) || any(pixloncorn > 180) || any(pixlatcorn < -90) || any(pixlatcorn > 90)
    Q = 0;
    return;
end
lonc = glon_vec(xx:xx+1);
latc = glat_vec(yy:yy+1);
gc_xall = [lonc(1), lonc(1), lonc(2), lonc(2), lonc(1)];
gc_yall = [latc(1), latc(2), latc(2), latc(1), latc(1)];
% ensure both are clockwise
[gc_xall, gc_yall] = poly2cw(gc_xall, gc_yall);
% now handled at the beginning of grid_to_gc
%[pixloncorn, pixlatcorn] = uncross_corners(pixloncorn, pixlatcorn);
[pixloncorn, pixlatcorn] = poly2cw(pixloncorn, pixlatcorn);
% create a polygon that represents the area of overlap and calculate its area
% in km.
[xt, yt] = polybool('intersection',gc_xall,gc_yall,pixloncorn,pixlatcorn);
if isempty(xt) || isempty(yt) || (any(sign(pixloncorn)~=sign(pixloncorn(1))) && (any(abs(pixloncorn)>90) || any(abs(pixlatcorn)>80)))
    % the last test handles issues where a pixel straddles the international
    % date line. there's better ways to handle it (wrap the pixel corner around
    % to be the same sign as the GC corners) but I don't feel like doing that atm.
    Q = 0;
    return;
elseif any(isnan(xt))
    error('load_and_grid_domino:calc_pix_grid_overlap','W%d: The pixel corners are wrong (%s, %s) at [%d, %d]',t.ID,mat2str(pixloncorn),mat2str(pixlatcorn), xx, yy);
end
overlap_area = areaint(yt,xt,earth_ellip);
Q = overlap_area/gc_area;
end

function [x,y]= uncross_corners(x,y)
    m1 = (y(3) - y(2))/(x(3) - x(2));
    b1 = y(2) - m1*x(2);
    m2 = (y(4) - y(1))/(x(4) - x(1));
    b2 = y(1) - m2*x(1);
    flip_bool = false;
    if ~isinf(m1) && ~isinf(m2)
        % As long as neither slope is infinite, solve for the x-coordinate
        % of the intercept and see if it falls inside the polygon - if so,
        % the corners need flipped.
        inpt = (b2-b1)/(m1-m2);
        if inpt > min(x(2:3)) && inpt < max(x(2:3))
            flip_bool = true;
        end
    elseif isinf(m1) && ~isinf(m2)
        % If one is infinite, fine the y-coord where the other one is at
        % it's x-coordinate and do the same test
        inpt = m2*x(2)+b2;
        if inpt > min(y(2:3)) && inpt < max(y(2:3))
            flip_bool = true;
        end
    elseif isinf(m2) && ~isinf(m1)
        inpt = m1*x(1) + b1;
        if inpt > min(y([1,4])) && inpt < max(y([1,4]))
            flip_bool = true;
        end
        % If both are infinite, they are parallel and the corners do not
        % need flipped.
    end
    if flip_bool
        tmp = x(4);
        x(4) = x(3);
        x(3) = tmp;
        tmp = y(4);
        y(4) = y(3);
        y(3) = tmp;
    end
end




function Data = calc_pixel_area(Data)
Lon1 = squeeze(Data.Loncorn(1,:,:));
Lon2 = squeeze(Data.Loncorn(2,:,:));
Lon3 = squeeze(Data.Loncorn(3,:,:));
Lon4 = squeeze(Data.Loncorn(4,:,:));

Lat1 = squeeze(Data.Latcorn(1,:,:));
Lat2 = squeeze(Data.Latcorn(2,:,:));
Lat3 = squeeze(Data.Latcorn(3,:,:));
Lat4 = squeeze(Data.Latcorn(4,:,:));

Amin = 312; % 13 km x 24 km
Amax = 7500; % approximate max seen from BEHR
Data.Areaweight = nan(size(Data.Longitude));

for x=1:size(Lon1,1)
    for y=1:size(Lon1,2)
        pixelarea = (m_lldist([Lon1(x,y)-180 Lon2(x,y)-180],[Lat1(x,y) Lat2(x,y)]))*(m_lldist([Lon1(x,y)-180, Lon4(x,y)-180],[Lat1(x,y), Lat4(x,y)]));
        %Data.Areaweight(x,y) = 1/pixelarea;
        Data.Areaweight(x,y) = 1 - (pixelarea - Amin)/Amax;
    end
end
end

function Data = read_fields(Data, hi, gnum, field_names)
% Reads in specified fields from the .he5 file given. hi should be the
% result of h5info.
for f=1:size(field_names,1)
    fillval = double(h5readatt(hi.Filename, h5dsetname(hi,1,2,1,gnum ,field_names{f,1}), 'MissingValue'));
    offset = double(h5readatt(hi.Filename, h5dsetname(hi,1,2,1,gnum ,field_names{f,1}), 'Offset'));
    scalefactor = double(h5readatt(hi.Filename, h5dsetname(hi,1,2,1,gnum ,field_names{f,1}), 'ScaleFactor'));
    
    
    Data.(field_names{f,2}) = double(Data.(field_names{f,2}));
    Data.(field_names{f,2}) = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,gnum ,field_names{f,1})));
    Data.(field_names{f,2})(Data.(field_names{f,2}) == fillval) = nan;
    Data.(field_names{f,2}) = (Data.(field_names{f,2}) - offset) * scalefactor;
end
end

