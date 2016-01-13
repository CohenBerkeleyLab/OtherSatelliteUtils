function [  ] = domino_omno2_comparison_plots( domino_dir, omno2_dir, fields_to_plot )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ischar(domino_dir) || ~ischar(omno2_dir) || ~exist(domino_dir, 'dir') || ~exist(omno2_dir, 'dir')
    error('domino_omno2_comparison_plots:bad_input','The inputs must be a directories given as strings');
end
if ~iscellstr(fields_to_plot) || size(fields_to_plot,2) ~= 2
    error('domino_omno2_comparison_plots:bad_input','fields_to_plot must be an n-by-2 cell array of strings');
end

% to be implemented later if necessary: use to restrict data to given time
% frame.
% if ~exist('start_date','var')
%     start_date = datestr(0,'yyyymmdd');
% else
%     start_date = datestr(start_date,'yyyymmdd');
% end
% 
% if ~exist('end_date','var')
%     end_date = datestr(0,'yyyymmdd');
% else
%     end_date = datestr(end_date,'yyyymmdd');
% end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DATA LOADING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

% Find all .mat files in the directories and load each one.
domino_files = dir(fullfile(domino_dir,'*.mat'));
omno2_files = dir(fullfile(omno2_dir,'*.mat'));

nfields = size(fields_to_plot,1);
fns_domino = fields_to_plot(:,1);
fns_omno2 = fields_to_plot(:,2);

% Each file will contain the variable GC_avg, which in turn has 

DOMINO = nan(1440, 720, numel(domino_files), nfields);
OMNO2 = nan(1440, 720, numel(omno2_files), nfields);

parfor a = 1:numel(domino_files)
    F_D = load(fullfile(domino_dir, domino_files(a).name));
    for f = 1:nfields
        DOMINO(:,:,a,f) = F_D.GC_avg.(fns_domino{f});
    end
    F_O = load(fullfile(omno2_dir, omno2_files(a).name));
    for f = 1:nfields
        OMNO2(:,:,a,f) = F_O.GC_avg.(fns_omno2{f});
    end
end



if numel(fns_domino) ~= numel(fns_omno2)
    error('domino_omno2_comparison_plots:bad_data','The DOMINO and OMNO2 variables have a different number of fields');
end

nfields = numel(fns_domino);

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
            domino_data = squeeze(DOMINO(a,b,:,f));
            omno2_data = squeeze(OMNO2(a,b,:,f));
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
    
end


end

