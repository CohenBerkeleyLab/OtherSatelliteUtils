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
    OMI = struct('AirMassFactorTropospheric', [], 'AveragingKernel', [], 'CloudFraction', [], 'CloudPressure', [],...
        'GhostColumn', [], 'SlantColumnAmountNO2', [], 'SlantColumnAmountNO2Std', [], 'SurfaceAlbedo', [],...
        'TerrainHeight', [], 'TroposphericVerticalColumn', [], 'TroposphericVerticalColumnError', [],...
        'TroposphericColumnFlag', {{}}, 'GroundPixelQualityFlags', {{}});
    OMI = repmat(OMI,1,numel(F));
    
    lonmin = -180;  lonmax = 180;
    latmin = -90;   latmax = 90;
    reslat = 2; reslon = 2.5;
    
    for s=1:numel(F)
        hi = h5info(fullfile(full_path, F(s).name));
        Data(s) = read_fields(Data(s), hi, 2, geo_fields);
        Data(s) = read_fields(Data(s), hi, 1, data_fields);
        OMI(s) = add2grid_DOMINO(Data(s),OMI(s),reslat,reslon,[lonmin, lonmax],[latmin, latmax]);
    end
    
    save_name = sprintf('OMI_DOMINO_%04d%02d%02d.mat',year(d),month(d),day(d));
    saveData(fullfile(save_path,save_name),Data,OMI);

end


end

function saveData(filename,Data,OMI)
    save(filename,'OMI','Data')
end

function read_fields(Data, hi, gnum, field_names)
% Reads in specified fields from the .he5 file given. hi should be the
% result of h5info.
for f=1:size(field_names,1)
    Data.(field_names{f,2}) = h5read(hi.Filename, h5dsetname(hi,1,2,1,gnum ,field_names{f,1}));
    if ~isinteger(Data.(field_names{f,2}))
        Data.(field_names{f,2}) = double(Data.(field_names{f,2}));
    end
end
end

