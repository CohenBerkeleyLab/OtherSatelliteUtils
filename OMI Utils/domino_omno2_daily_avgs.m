function [  ] = domino_omno2_daily_avgs( start_date, end_date )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   Note that this function has now been wrapped into the load_and_grid functions
%   making it obsolete.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% USER DEFINED VALUES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Begin function domino_omno2_daily_avgs\n');
warning('The functionality of this script should be in the various load_and_grid functions. Are you sure you should be using this one?');


% Define paths to the gridded files. This should be the path that has two
% subfolders, DOMINO and OMNO2, which each have under them data organised
% by resolution and year. The data should be the output of the
% "load_and_grid" functions, and so should have the variables "Data" and
% "GC" in them.
global onCluster
if isempty(onCluster)
    onCluster = false;
end
if onCluster
    data_path = '/global/home/users/laughner/myscratch/MATLAB/Data/OMI/';
    save_path = '/global/home/users/laughner/myscratch/MATLAB/Data/OMI/DOMINO-OMNO2-Comparison';
else
    data_path = '/Volumes/share2/USERS/LaughnerJ/DOMINO-OMNO2_comparision';
    save_path = '/Users/Josh/Documents/MATLAB/Non BEHR Satellite/Workspaces/DOMINO-OMNO2-Comparison';
end

% Define the resolution. Used for the path and the initialization of the
% matrices.
lonres = 0.25;
latres = 0.25;

% These should be the names of the fields in the GC structure that you want
% to compare across. The first column should be the DOMINO field name, the
% second the OMNO2 field name.

fields_to_plot = {  'TroposphericVerticalColumn', 'ColumnAmountNO2Trop';...
                    'AirMassFactorTropospheric', 'AmfTrop';...
                    'AssimilatedStratosphericVerticalColumn', 'ColumnAmountNO2Strat';...
                    'SlantColumnAmountNO2', 'SlantColumnAmountNO2Destriped'};

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Main function beginning\n');

% We need to load all the days within the date range given, average the
% swaths down to a single matrix for that day, and put that average into
% the larger matrix.

%dates = datenum(sprintf('%04d-01-01',yr)):datenum(sprintf('%04d-12-31',yr));
dates = datenum(start_date):datenum(end_date);

ndays = numel(dates);
nfields = size(fields_to_plot,1);
nlons = 360/lonres;
nlats = 180/latres;

fprintf('Initializing large matrices\n');

daily_domino = nan(nlons,nlats,ndays,nfields);
daily_omno2 = nan(nlons,nlats,ndays,nfields);

% Loading the variables will take a long time so run it in parallel. (This
% may not actually save time, it depends on if the files system is
% optimized for parallel reads).

fprintf('Starting parfor loop\n');
%parfor d=1:ndays
for d=1:ndays
    fprintf('Loading files for %s\n',datestr(dates(d),'yyyy-mm-dd'));
    tmp_mats_domino = load_and_avg('DOMINO', dates(d), lonres, nlons, latres, nlats, fields_to_plot(:,1), data_path); %#ok<PFBNS>
    tmp_mats_omno2 = load_and_avg('OMNO2', dates(d), lonres, nlons, latres, nlats, fields_to_plot(:,2), data_path);
    for a = 1:nfields
        daily_domino(:,:,d,a) = tmp_mats_domino(:,:,a);
        daily_omno2(:,:,d,a) = tmp_mats_omno2(:,:,a);
    end
end

for a=1:nfields
    DOMINO.(fields_to_plot{a,1}) = daily_domino(:,:,:,a);
    OMNO2.(fields_to_plot{a,2}) = daily_omno2(:,:,:,a);
end
savename = sprintf('domino-omno2_avgs_%s-%s.mat',datestr(start_date,'yyyymmdd'),datestr(end_date,'yyyymmdd'));
fprintf('Saving as %s\n',fullfile(save_path,savename));
save(fullfile(save_path, savename), 'DOMINO', 'OMNO2','-v7.3');

end

function avg_mats = load_and_avg(product, date_i, lonres, nlons, latres, nlats, fns, data_path)
fname = sprintf('OMI_%s_%04d%02d%02d.mat',product,year(date_i),month(date_i),day(date_i));
fullpath = fullfile(data_path,product,sprintf('%sx%s',num2str(lonres),num2str(latres)),fname); % num2str intelligently trims trailing zeros

avg_mats = nan(nlons, nlats, numel(fns));

if exist(fullpath,'file')
    D = load(fullpath,'GC');
    aw = cat(3, D.GC.Areaweight);
    
    for a=1:numel(fns)
        data = cat(3, D.GC.(fns{a}));
        avg_mats(:,:,a) = nansum2(data .* aw, 3) ./ nansum2(aw,3);
    end
else
    fprintf('File %s not found\n', fullpath)
end
end
