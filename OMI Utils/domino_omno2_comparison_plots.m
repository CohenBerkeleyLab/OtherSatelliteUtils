function [  ] = domino_omno2_comparison_plots( files )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscellstr(files)
    error('domino_omno2_comparison_plots:bad_input','The input "files" must be a cell array of filenames as strings');
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DATA LOADING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

% The files should all contain two variables: DOMINO and OMNO2, both of
% which will be structures with the respective variables to compare as the
% fields. These will need to be concatenated along the third dimension, as
% the function that put them together can be split up in time to make the
% computation manageable.

D(numel(files)).DOMINO = struct;
D(numel(files)).OMNO2 = struct;

for a = 1:numel(files)
    D(a) = load(files{a});
end

fns_domino = fieldnames(D(1).DOMINO);
fns_omno2 = fieldnames(D(1).OMNO2);

if numel(fns_domino) ~= numel(fns_omno2)
    error('domino_omno2_comparison_plots:bad_data','The DOMINO and OMNO2 variables have a different number of fields');
end

nfields = numel(fns_domino);

for a = 1:numel(fns_domino)
    DOMINO.(fns_domino{a}) = cat(3, D.DOMINO.(fns_domino{a}));
    OMNO2.(fns_omno2{a}) = cat(3, D.OMNO2.(fns_domino{a}));
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
            domino_data = squeeze(DOMINO.(fns_domino{f})(a,b,:));
            omno2_data = squeeze(OMNO2.(fns_omno2{f})(a,b,:));
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

