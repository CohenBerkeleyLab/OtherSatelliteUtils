function [  ] = load_and_grid_DOMINO( start_date, end_date )
%LOAD_AND_GRID_DOMINO Read in and grid DOMINO data.
%   Current implementation will simply assume that you wish to load
%   worldwide data.  This will simplify the coding slightly, but will need
%   to be modified later if you want to only load e.g. US data (which would
%   be important if you want to use a finer gridding resolution).

E = JLLErrors;

global onCluster
if isempty(onCluster)
    onCluster = false;
end

if onCluster
    domino_path = '/global/home/users/laughner/myscratch/SAT/OMI/DOMINOv2';
    save_path = '/global/home/users/laughner/myscratch/MATLAB/Data/OMI/DOMINO';
else
    domino_path = '/Volumes/share-sat/SAT/OMI/DOMINOv2.0';
    save_path = '/Volumes/share-sat/SAT/OMI/DOMINOv2.0/MatFiles';
end

% These should be nx2 cell arrays; the first column lists the dataset name
% in the HDFv5 file, the second what you want it to be called in the Data
% structure. Make the second column an empty string if you want it to be
% the same as the first. data_fields will be read from the 'Data Fields'
% group, and geo_fields from the "Geolocation Fields" group.
data_fields = { 'AirMassFactorTropospheric','';...
                'AveragingKernel','';...
                'CloudFraction','';...
                'CloudPressure','';...
                'GhostColumn','';...
                'SlantColumnAmountNO2','';...
                'SlantColumnAmountNO2Std','';...
                'SurfaceAlbedo','';...
                'TerrainHeight','';...
                'TroposphericVerticalColumn','';...
                'TroposphericVerticalColumnError',''};   
geo_fields = {  'Longitude','';...
                'LongitudeCornerpoints','Loncorn';...
                'Latitude','';...
                'LatitudeCornerpoints','Latcorn';...
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

for d=datenum(start_date):datenum(end_date)
    % DOMINO data should be organized in subfolders by year and month
    full_path = fullfile(domino_path, sprintf('%04d',year(d)), sprintf('%02d',month(d)));
    file_pat = sprintf('OMI-Aura_L2-OMDOMINO_%04dm%02d%02d*.he5',year(d),month(d),day(d));
    F = dir(fullfile(full_path, file_pat));
    Data = make_empty_struct_from_cell(all_struct_fields);
    Data = repmat(Data,1,numel(F));
    
    omi_mat = nan(144,91);
    omi_cell = cell(144,91);
    GC = struct('AirMassFactorTropospheric', omi_mat, 'AveragingKernel', omi_mat, 'CloudFraction', omi_mat, 'CloudPressure', omi_mat,...
        'GhostColumn', omi_mat, 'SlantColumnAmountNO2', omi_mat, 'SlantColumnAmountNO2Std', omi_mat, 'SurfaceAlbedo', omi_mat,...
        'TerrainHeight', omi_mat, 'TroposphericVerticalColumn', omi_mat, 'TroposphericVerticalColumnError', omi_mat,...
        'Area', omi_mat,'Areaweight', omi_mat, 'TroposphericColumnFlag', omi_cell, 'GroundPixelQualityFlags', omi_cell);
    GC = repmat(GC,1,numel(F));
    
    %lonmin = -180;  lonmax = 180;
    %latmin = -90;   latmax = 90;
    %reslat = 2; reslon = 2.5;
    
    for s=1:numel(F)
        hi = h5info(fullfile(full_path, F(s).name));
        Data(s) = read_fields(Data(s), hi, 2, geo_fields);
        Data(s) = read_fields(Data(s), hi, 1, data_fields);
        %OMI(s) = add2grid_DOMINO(Data(s),OMI(s),reslat,reslon,[lonmin, lonmax],[latmin, latmax]);
        GC(s) = grid_to_gc(Data(s), GC(s));
    end
    
    save_name = sprintf('OMI_DOMINO_%04d%02d%02d.mat',year(d),month(d),day(d));
    saveData(fullfile(save_path,save_name),Data,GC);

end


end

function saveData(filename,Data, GC)
    save(filename,'GC','Data')
end

function GC = grid_to_gc(Data, GC)
[gloncorn, glatcorn] = geos_chem_corners();
gloncorn = gloncorn';
glatcorn = glatcorn';

sz = size(gloncorn)-1;

Data = calc_pixel_area(Data);

% Filter for clouds, row anomaly and albedo (DOMINO recommends only using
% albedo < 0.3 as snow/ice can interfere in cloud retrieval)
row_anom = Data.TroposphericColumnFlag < 0;
clds = Data.CloudFraction > 0.3;
alb = Data.SurfaceAlbedo > 0.3;

for a=1:sz(1)
    for b=1:sz(2)
        xx = Data.Longitude >= gloncorn(a,b) & Data.Longitude < gloncorn(a+1,b+1);
        yy = Data.Latitude >= glatcorn(a,b) & Data.Latitude < glatcorn(a+1,b+1);
        fns = fieldnames(GC);
        for c=1:numel(fns)
            Data.(fns{c})(row_anom | clds | alb) = nan;
            if iscell(GC.(fns{c}))
                GC.(fns{c}){a,b} = Data.(fns{c})(xx&yy);
            else
                % Calculate an area-weighted mean
                GC.(fns{c})(a,b) = nansum2(Data.(fns{c}).(xx&yy) .* Data.Areaweight(xx&yy)) ./ nansum2(Data.Areaweight(xx&yy));
            end
        end
    end
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

Data.Areaweight = nan(size(Data.Longitude));

for x=1:size(Lon1,1)
    for y=size(Lon1,2)
        pixelarea = (m_lldist([Lon1(x,y)-180 Lon2(x,y)-180],[Lat1(x,y) Lat2(x,y)]))*(m_lldist([Lon1(x,y)-180, Lon4(x,y)-180],[Lat1(x,y), Lat4(x,y)]));
        Data.Areaweight(x,y) = 1/pixelarea;
    end
end
end

function Data = read_fields(Data, hi, gnum, field_names, keep_int)
% Reads in specified fields from the .he5 file given. hi should be the
% result of h5info.
for f=1:size(field_names,1)
    fillval = h5readatt(hi.Filename, h5dsetname(hi,1,2,1,gnum ,field_names{f,1}), 'MissingValue');
    offset = h5readatt(hi.Filename, h5dsetname(hi,1,2,1,gnum ,field_names{f,1}), 'Offset');
    scalefactor = h5readatt(hi.Filename, h5dsetname(hi,1,2,1,gnum ,field_names{f,1}), 'ScaleFactor');
    
    Data.(field_names{f,2}) = (h5read(hi.Filename, h5dsetname(hi,1,2,1,gnum ,field_names{f,1})) - offset) * scalefactor;
    Data.(field_names{f,2})(Data.(field_names{f,2}) == fillval) = nan;
    
    if  ~keep_int || ~isinteger(Data.(field_names{f,2}))
        Data.(field_names{f,2}) = double(Data.(field_names{f,2}));
    end
end
end

