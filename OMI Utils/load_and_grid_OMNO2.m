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

if onCluster
    omno2_path = '/global/home/users/laughner/myscratch/SAT/OMI/OMNO2';
    save_path = '/global/home/users/laughner/myscratch/MATLAB/Data/OMI/OMNO2/0.25x0.25';
else
    omno2_path = '/Volumes/share-sat/SAT/OMI/OMNO2';
    save_path = '/Volumes/share-sat/SAT/OMI/OMNO2/MatFiles';
end

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

xres = 0.25; yres = 0.25;
[loncorn, latcorn] = grid_corners(xres, yres);
loncorn = loncorn';
lon = loncorn(1:end-1,1:end-1) + xres/2;
latcorn = latcorn';
lat = latcorn(1:end-1,1:end-1) + yres/2;

for d=datenum(start_date):datenum(end_date)
    if DEBUG_LEVEL > 0; fprintf('Loading/gridding data for %s\n',datestr(d)); end

    % OMNO2 data should be organized in subfolders by year and month
    full_path = fullfile(omno2_path, sprintf('%04d',year(d)), sprintf('%02d',month(d)));
    file_pat = sprintf('OMI-Aura_L2-OMNO2_%04dm%02d%02d*.he5',year(d),month(d),day(d));
    F = dir(fullfile(full_path, file_pat));
    Data = make_empty_struct_from_cell(all_struct_fields);
    Data = repmat(Data,1,numel(F));
    
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
    
    for s=1:numel(F)
        if DEBUG_LEVEL > 1; fprintf('\tHandling file %d of %d\n',s,numel(F)); end
        hi = h5info(fullfile(full_path, F(s).name));
        Data(s).Filename = F(s).name;
        Data(s) = read_fields(Data(s), hi, 2, geo_fields);
        Data(s) = read_fields(Data(s), hi, 1, data_fields);
        %OMI(s) = add2grid_DOMINO(Data(s),OMI(s),reslat,reslon,[lonmin, lonmax],[latmin, latmax]);
        GC(s) = grid_to_gc(Data(s), GC(s),loncorn,latcorn);
    end

    for s=1:numel(GC)
        % Add these now so that grid_to_gc doesn't try to grid them
        GC(s).Longitude = lon;
        GC(s).Latitude = lat;
    end

    save_name = sprintf('OMI_OMNO2_%04d%02d%02d.mat',year(d),month(d),day(d));
    saveData(fullfile(save_path,save_name),Data,GC);

end


end

function saveData(filename,Data, GC)
    save(filename,'GC','Data','-v7.3')
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

function GC = grid_to_gc(Data, GC, gloncorn, glatcorn)
%[gloncorn, glatcorn] = geos_chem_corners();

sz = size(gloncorn)-1;

Data = calc_pixel_area(Data);
TotalAreaweight = zeros(size(GC.ColumnAmountNO2Trop));

% Filter for clouds, row anomaly and albedo (DOMINO recommends only using
% albedo < 0.3 as snow/ice can interfere in cloud retrieval)
row_anom = Data.XTrackQualityFlags > 0;
vcd_qual = Data.VcdQualityFlags > 0;
clds = Data.CloudFraction > 0.3;
alb = Data.TerrainReflectivity > 0.3;
negvcds = Data.ColumnAmountNO2Trop < 0;
extreme_vcd = Data.ColumnAmountNO2Trop > 1e17;

fns = fieldnames(GC);
for c=1:numel(fns)
    Data.(fns{c})(row_anom | vcd_qual | clds | alb | negvcds | extreme_vcd) = nan;
end

glon_vec = gloncorn(:,1);
glat_vec = glatcorn(1,:);

fprintf('Starting loop\n');
tval = tic;
for a=1:numel(Data.Longitude)
    %if mod(a,100) == 1
     %   fprintf('Pixel %d: time elapsed = %f\n',a,toc(tval));
    %end
    if Data.Longitude(a) < -200 || Data.Latitude(a) < -200 || isnan(Data.Longitude(a)) || isnan(Data.Latitude(a))
        % Fill values are ~1x10^-30, this will skip any pixels with
        % fill values for coordinates
        continue
    end
    xx = Data.Longitude(a) >= glon_vec(1:end-1) & Data.Longitude(a) < glon_vec(2:end);
    yy = Data.Latitude(a) >= glat_vec(1:end-1) & Data.Latitude(a) < glat_vec(2:end);
    if sum(xx) > 1 || sum(yy) > 1
        error('load_and_grid_OMNO2:pixel_assignment_error', 'Pixel %d in %s fell into more than 1 grid cell', a, Data.Filename);
    elseif sum(xx) < 1 || sum(yy) < 1
        error('load_and_grid_OMNO2:pixel_assignment_error', 'Pixel %d in %s could not be assigned to a grid cell', a, Data.Filename);
    end
    TotalAreaweight(xx,yy) = nansum2([TotalAreaweight(xx,yy), Data.Areaweight(a)]);
    for c=1:numel(fns)
        if iscell(GC.(fns{c}))
           GC.(fns{c}){xx,yy} = cat(1,GC.(fns{c}){xx,yy},Data.(fns{c})(a));
        else
            % The areaweight will be divided out at the end of the loop
            GC.(fns{c})(xx,yy) = nansum2([GC.(fns{c})(xx,yy), Data.(fns{c})(a) .* Data.Areaweight(a)]);
        end 
    end
end

for c=1:numel(fns)
    if ~iscell(GC.(fns{c}))
        GC.(fns{c}) = GC.(fns{c}) ./ TotalAreaweight;
    end
end
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

function Data = calc_pixel_area(Data)
Lon1 = Data.Loncorn(:,:,1);
Lon2 = Data.Loncorn(:,:,2);
Lon3 = Data.Loncorn(:,:,3);
Lon4 = Data.Loncorn(:,:,4);

Lat1 = Data.Latcorn(:,:,1);
Lat2 = Data.Latcorn(:,:,2);
Lat3 = Data.Latcorn(:,:,3);
Lat4 = Data.Latcorn(:,:,4);

Data.Areaweight = nan(size(Data.Longitude));

for x=1:size(Lon1,1)
    for y=1:size(Lon1,2)
        pixelarea = (m_lldist([Lon1(x,y)-180 Lon2(x,y)-180],[Lat1(x,y) Lat2(x,y)]))*(m_lldist([Lon1(x,y)-180, Lon4(x,y)-180],[Lat1(x,y), Lat4(x,y)]));
        Data.Areaweight(x,y) = 1/pixelarea;
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
    Data.(field_names{f,2}) = (double(h5read(hi.Filename, h5dsetname(hi,1,2,1,gnum ,field_names{f,1}))) - offset) * scalefactor;
    Data.(field_names{f,2})(Data.(field_names{f,2}) == fillval) = nan;
end
end

