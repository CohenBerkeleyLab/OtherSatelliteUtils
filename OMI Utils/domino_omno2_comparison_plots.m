function [  ] = domino_omno2_comparison_plots( product1, product2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

product1 = lower(product1);
product2 = lower(product2);

allowed_products = {'domino','omno2','behr'};

if ~ismember(product1, allowed_products) || ~ismember(product2, allowed_products)
    error('domino_omno2_comparison_plots:bad_input','The allowed products are %s.',strjoin(allowed_products, ', '));
end


% Always make BEHR the second product - this'll make it easier to match the smaller lat/lon array
% down the road
if strcmpi(product1,'behr')
    product1 = product2;
    product2 = 'behr';
end

[product1_dir, product1_fields] = set_product(product1);
[product2_dir, product2_fields] = set_product(product2);

if numel(product2_fields) ~= numel(product1_fields)
    error('domino_omno2_comparison_plots:different_nfields','The two products have different numbers of fields, remember to fill the cell arrays with empty strings (see comments in subfunction set_product');
end

xx = ~iscellcontents(product1_fields, 'isempty') & ~iscellcontents(product2_fields, 'isempty');
product1_fields = product1_fields(xx);
product2_fields = product2_fields(xx);

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

save_path = '/global/home/users/laughner/myscratch/MATLAB/Data/OMI/DOMINO-OMNO2-Comparison';

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DATA LOADING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

% Find all .mat files in the directories and load each one.
p1_files = dir(fullfile(product1_dir,'*.mat'));
p2_files = dir(fullfile(product2_dir,'*.mat'));

nfields = numel(product1_fields);


% Each file will contain the variable GC_avg, which in turn has all the averaged
% fields in it.

% We'll need to load the first files to determine how to line up the 
% matrices if using the BEHR product.

F1 = load(fullfile(product1_dir, p1_files(1).name));
F2 = load(fullfile(product2_dir, p2_files(1).name));

[xx,yy] = match_geo(F1,F2);

clear('F1','F2');

nlons = sum(xx);
nlats = sum(yy);

P1 = nan(nlons, nlats, numel(p1_files), nfields);
P2 = nan(nlons, nlats, numel(p2_files), nfields);

fprintf('Starting loop over files\n');

for a = 1:numel(p1_files)
    fprintf('Loading files %s and %s\n', p1_files(a).name, p2_files(a).name);
    fprintf('File 1\n')
    F1 = load(fullfile(product1_dir, p1_files(a).name));
    Ptmp = nan(sum(xx), sum(yy), nfields);
    for f = 1:nfields
        if numel(xx) ~= size(F1.GC_avg.Longitude,1) || numel(yy) ~= size(F1.GC_avg.Latitude,2)
            error('domino_omno2_comparison_plots:inconsistent_coordinates','The dimensions of the matrices in product 1 are different in different files');
        end
        Ptmp(:,:,f) = F1.GC_avg.(product1_fields{f})(xx,yy);
    end
    P1(:,:,a,:) = Ptmp;
    
    fprintf('File 2\n')
    F2 = load(fullfile(product2_dir, p2_files(a).name));
    Ptmp = nan(sum(xx), sum(yy), nfields);
    for f = 1:nfields
        if sum(xx) ~= size(F2.GC_avg.Longitude,1) || sum(yy) ~= size(F2.GC_avg.Latitude,2)
            error('domino_omno2_comparison_plots:inconsistent_coordinates','The dimensions of the matrices in product 2 are different in different files');
        end
        Ptmp(:,:,f) = F2.GC_avg.(product2_fields{f});
    end
    P2(:,:,a,:) = Ptmp;
end


% Now actually go through and for each lat/lon box and each variable,
% calculate the slope, intercept, R2, and std. deviations of the P2
% variables vs. the P1 variables

slopes = nan(nlons, nlats, nfields);
intercepts = nan(nlons, nlats, nfields);
r2s = nan(nlons, nlats, nfields);
stddevms = nan(nlons, nlats, nfields);
stddevbs = nan(nlons, nlats, nfields);

w=warning('off','all');
parfor a=1:nlons
%for a=1:nlons
    warning('off','all')
    P1_slice = P1(a,:,:,:);
    P2_slice = P2(a,:,:,:);
    slopes_a = nan(1,nlats,nfields);
    int_a = nan(1,nlats,nfields);
    r2_a = nan(1,nlats,nfields);
    sdm_a = nan(1,nlats,nfields);
    sdb_a = nan(1,nlats,nfields);
    for b=1:nlats
        for f=1:nfields
            p1_data = squeeze(P1_slice(1,b,:,f));
            p2_data = squeeze(P2_slice(1,b,:,f));
            [~,~,~,L] = calc_fit_line(p1_data, p2_data, 'regression', 'rma');
            slopes_a(1,b,f) = L.P(1);
            int_a(1,b,f) = L.P(2);
            r2_a(1,b,f) = L.R2;
            sdm_a(1,b,f) = L.StdDevM;
            sdb_a(1,b,f) = L.StdDevB;
        end
    end
    slopes(a,:,:) = slopes_a;
    intercepts(a,:,:) = int_a;
    r2s(a,:,:) = r2_a;
    stddevms(a,:,:) = sdm_a;
    stddevbs(a,:,:) = sdb_a;
end
warning(w);

% Finally, save the results. Do not overwrite existing files.
savename_spec = '%s-%s_results-%d.mat';
i=0;
savename = sprintf(savename_spec, product1, product2, i);
while exist(fullfile(save_path, savename),'file')
    i=i+1;
    savename = sprintf(savename_spec, product1, product2, i);
end
save(fullfile(save_path, savename), 'slopes', 'intercepts', 'r2s', 'stddevms', 'stddevbs', 'product1_fields', 'product2_fields');

end

function [p_dir, p_fields] = set_product(product)
% Sets both the directory and the fields for each product. Each product
% must have the same number of fields; if one product does not have an
% equivalent field (e.g. BEHR does not have stratospheric columns) put an
% empty string there for that product.
global onCluster;
if isempty(onCluster); onCluster = false; end 
switch product
    case 'domino'
        if onCluster
            p_dir = '/global/scratch/laughner/MATLAB/Data/OMI/DOMINO/0.25x0.25-avg';
        else
            error('domino_omno2_comparison_plots:not_implemented','Off cluster paths are not defined');
        end
        p_fields = {'TroposphericVerticalColumn'; 'AirMassFactorTropospheric'; 'AssimilatedStratosphericVerticalColumn'; };
    case 'omno2'
        if onCluster
            p_dir = '/global/scratch/laughner/MATLAB/Data/OMI/OMNO2/0.25x0.25-avg';
        else
            error('domino_omno2_comparison_plots:not_implemented','Off cluster paths are not defined');
        end
        p_fields = {'ColumnAmountNO2Trop'; 'AmfTrop'; 'ColumnAmountNO2Strat'};
    case 'behr'
        if onCluster
            p_dir = '/global/scratch/laughner/MATLAB/Data/OMI/BEHR/0.25x0.25-avg';
        else
            error('domino_omno2_comparison_plots:not_implemented','Off cluster paths are not defined');
        end
        p_fields = {'BEHRColumnAmountNO2Trop'; 'BEHRAMFTrop'; ''};
end
end

function [xx,yy] = match_geo(F1, F2)
% F2 must contain the smaller matrices
lon1 = F1.GC_avg.Longitude(:,1);
lat1 = F1.GC_avg.Latitude(1,:);
lon2 = F2.GC_avg.Longitude(:,1);
lat2 = F2.GC_avg.Latitude(1,:);

xx = ismember(lon1, lon2);
yy = ismember(lat1, lat2);

if sum(xx) ~= numel(lon2) || sum(yy) ~= numel(lat2)
    error('domino_omno2_comparison_plots:grid_mismatch','Lat/lon grid mismatch. It is assumed that the lat/lon coordinates in product2 are a subset of those in product1, and this is not true.')
end
end
