function [  ] = domino_omno2_comparison_plots( )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% USER DEFINED VALUES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Set to true if you want to save the matrices of daily averaged values to
% work with separately later.
save_daily_avg_mats = true;

% Define the year to compare here.
yr = 2012;

% Define the resolution. Used for the path and the initialization of the
% matrices.
lonres = 0.25;
latres = 0.25;

% These should be the names of the fields in the GC structure that you want
% to compare across. The first column should be the DOMINO field name, the
% second the OMNO2 field name.

fields_to_plot = {  'TroposphericVerticalColumn', 'ColumnAmountNO2Trop'};%;...
    %'AirMassFactorTropospheric', 'AmfTrop';...
    %'AssimilatedStratosphericVerticalColumn', 'ColumnAmountNO2Strat';...
    %'SlantColumnAmountNO2', 'SlantColumnAmountNO2Destriped'};

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% The first step will be to load all the days of the year, average the
% swaths down to a single matrix for that day, and put that average into
% the larger matrix.

%dates = datenum(sprintf('%04d-01-01',yr)):datenum(sprintf('%04d-12-31',yr));
dates = datenum(sprintf('%04d-01-01',yr)):datenum(sprintf('%04d-01-03',yr));

ndays = numel(dates);
nfields = size(fields_to_plot,1);
nlons = 360/lonres;
nlats = 180/latres;

daily_domino = nan(nlons,nlats,ndays,nfields);
daily_omno2 = nan(nlons,nlats,ndays,nfields);

% Loading the variables will take a long time so run it in parallel. (This
% may not actually save time, it depends on if the files system is
% optimized for parallel reads).

%parfor d=1:ndays
for d=1:ndays
    tmp_mats_domino = load_and_avg('DOMINO', dates(d), lonres, nlons, latres, nlats, fields_to_plot(:,1), data_path); %#ok<PFBNS>
    tmp_mats_omno2 = load_and_avg('OMNO2', dates(d), lonres, nlons, latres, nlats, fields_to_plot(:,2), data_path);
    for a = 1:nfields
        daily_domino(:,:,d,a) = tmp_mats_domino(:,:,a);
        daily_omno2(:,:,d,a) = tmp_mats_omno2(:,:,a);
    end
end

% Now actually go through and for each lat/lon box and each variable,
% calculate the slope, intercept, R2, and std. deviations of the OMNO2
% variables vs. the DOMINO variables

slopes = nan(nlons, nlats, nfields);
intercepts = nan(nlons, nlats, nfields);
r2s = nan(nlons, nlats, nfields);
stddevms = nan(nlons, nlats, nfields);
stddevbs = nan(nlons, nlats, nfields);

w=warning('off','all');
%parfor a=1:nlons
for a=1:nlons
    for b=1:nlats
        for f=1:nfields
            domino_data = squeeze(daily_domino(a,b,:,f));
            omno2_data = squeeze(daily_omno2(a,b,:,f));
            [~,~,~,L] = calc_fit_line(domino_data, omno2_data, 'regression', 'rma');
            slopes(a,b,f) = L.P(1);
            intercepts(a,b,f) = L.P(2);
            r2s(a,b,f) = L.R2;
            stddevms(a,b,f) = L.StdDevM;
            stddevbs(a,b,f) = L.StdDevB;
        end
    end
end
warning(w);

% Finally, save the results. Do not overwrite existing files.
savename_spec = 'domino-omno2_results-%d.mat';
i=0;
savename = sprintf(savename_spec, i);
while exist(fullfile(save_path, savename),'file')
    i=i+1;
    savename = sprintf(savename_spec, i);
end
save(fullfile(save_path, savename), 'slopes', 'intercepts', 'r2s', 'stddevms', 'stddevbs', 'fields_to_plot');

% If the user wants to save the averaged matrices, include those as a
% separate file. Do some variable rearraging to make it a little more
% obvious what everything is.
if save_daily_avg_mats
    for a=1:nfields
        DOMINO.(fields_to_plot{a,1}) = daily_domino(:,:,:,a);
        OMNO2.(fields_to_plot{a,2}) = daily_omno2(:,:,:,a);
    end
    savename2_spec = 'domino-omno2_daily_mats-%d.mat';
    i=0;
    savename2 = sprintf(savename2_spec, i);
    while exist(fullfile(save_path, savename2), 'file')
        i=i+1;
        savename2 = sprintf(savename2_spec, i);
    end
    save(fullfile(save_path, savename2), 'DOMINO', 'OMNO2','-v7.3');
end


end

function avg_mats = load_and_avg(product, date_i, lonres, nlons, latres, nlats, fns, data_path)
fname = sprintf('OMI_%s_%04d%02d%02d.mat',product,year(date_i),month(date_i),day(date_i));
fullpath = fullfile(data_path,product,sprintf('%sx%s',num2str(lonres),num2str(latres)),fname); % num2str intelligently trims trailing zeros

avg_mats = nan(nlons, nlats, numel(fns));

if exist(fullpath,'file')
    D = load(fullpath,'GC');
    aw = cat(3, D.GC.Areaweight);
    
    % DEBUG ONLY %
    aw = aw(1:end-1,1:end-1,:);
    %            %
    
    for a=1:numel(fns)
        data = cat(3, D.GC.(fns{a}));
        % DEBUG ONLY %
        data = data(1:end-1,1:end-1,:);
        %            %
        avg_mats(:,:,a) = nansum2(data .* aw, 3) ./ nansum2(aw,3);
    end
else
    fprintf('File %s not found\n', fullpath)
end
end