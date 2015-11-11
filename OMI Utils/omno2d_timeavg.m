function [ no2 ] = omno2d_timeavg( start_date, end_date, varargin )
%OMNO2D_TIMEAVG Averages cloud screened OMNO2 data over a given time period 
%   OMNO2d is the gridded, Level 3 product for OMI NO2. Since it is already
%   gridded, it is easy to average over time. Inputs are:
%     (optional)
%       start_date, end_date - first and last days to average. Any MATLAB
%       date format is acceptible.
%
%     (parameters)
%       dayofweek - can be any of 'all', 'weekday', 'weekend',
%       'r_weekday', 'r_weekend' which determine which days of week to
%       include.
%           'all' - all days
%           'weekday', - M through F
%           'weekend', - Sat & Sun
%           'r_weekday', - Tu through F
%           'r_weekend', - Sun only
%       Useful for weekend effect stuff.
%
%       cloud_screened - true/false whether to use cloud screened data or
%       not. Defaults to true. OMNO2d screens for cloud fraction < 30% (I
%       believe this is geometric, not radiance, but the documentation is
%       unclear).  Must be a scalar boolean or numeric value.
%
%       ignore_missing - true/false whether or not to error if a day's file
%       is missing. Defaults to false, i.e. will error if a file is
%       missing. Must be a scalar boolean or numeric value.
%
%   This will load each file and its tropospheric NO2 column, and do a
%   weighted average using the "weight" field.
%
%   The path to the OMNO2d folder is defined below; this may need changed
%   on your system to point to the appropriate file server location.
%
%   Josh Laughner <joshlaugh5@gmail.com> 10 Nov 2015
%
%   Dependencies: 
%       BEHR_MatlabClasses_GitRepo.git/JLLErrors.m
%       BEHR_GitRepo.git/Utils/h5dsetname.m

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;
DEBUG_LEVEL = 1;

% start and end dates must be an understood date format
try
    s_datenum = datenum(start_date);
    e_datenum = datenum(end_date);
catch err
    if strcmp(err.identifier,'MATLAB:datenum:ConvertDateString')
        E.badinput('start_date and end_date must be in MATLAB-understood date formats')
    else
        rethrow(err)
    end
end

% parse the parameter arguments
p = inputParser;
p.addParameter('dayofweek','all',@(x) ismember(lower(x), {'all', 'weekday', 'weekend', 'r_weekday', 'r_weekend'}));
p.addParameter('cloud_screened',true,@(x) (isscalar(x) && (isnumeric(x) || islogical(x))));
p.addParameter('ignore_missing',false,@(x) (isscalar(x) && (isnumeric(x) || islogical(x))));
p.parse(varargin{:});
pout = p.Results;

dayofweek = lower(pout.dayofweek);
cloud_screened = pout.cloud_screened;
ignore_missing = pout.ignore_missing;

% define the path to the OMNO2d data. This should be in volume 1 of the
% Synology file server (IP = 128.32.208.13 as of 10 Nov 2015). This should
% be a folder containing folders by year
omno2d_path = '/Volumes/share-sat/SAT/OMI/OMNO2d';
if ~exist(omno2d_path,'dir')
    E.dir_dne('omno2d_path')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% setup which days of week to include. weekday in Matlab 2014b defines 1 =
% Sun, 7 = Sat
switch dayofweek
    case 'all'
        days = 1:7;
    case 'weekday'
        days = 2:6;
    case 'weekend'
        days = [1,7];
    case 'r_weekday'
        days = 3:6;
    case 'r_weekend'
        days = 1;
end

dates = s_datenum:e_datenum;

% set the dataset name based on whether or not to use cloud screen
if cloud_screened
    dataset_name = 'ColumnAmountNO2TropCloudScreened';
else
    dataset_name = 'ColumnAmountNO2Trop';
end

% OMNO2d is gridded on a 1440x720 grid of 0.25 degree cells. Set up the
% matrices to receive the data from each day.

weighted_no2 = zeros(1440,720);
weights = zeros(1440,720);

for d=1:numel(dates)
    if ~ismember(weekday(dates(d)),days)
        % skip based on day-of-week
        continue
    end
    
    if DEBUG_LEVEL > 0
        fprintf('Adding %s\n',datestr(dates(d)));
    end
    
    % define the OMNO2d filename. Check that exactly one file exists,
    % unless ignore_missing is on, then then check that no more than 1
    % exists.
    fname = sprintf('OMI-Aura_L3-OMNO2d_%04dm%02d%02d_v003*.he5', year(dates(d)), month(dates(d)), day(dates(d)));
    F = dir(fullfile(omno2d_path, sprintf('%04d',year(dates(d))),fname));
    if isempty(F) && ~ignore_missing
        E.filenotfound(fname);
    elseif numel(F) > 1
        E.toomanyfiles(fname);
    end
    
    hinfo = h5info(fullfile(omno2d_path, sprintf('%04d',year(dates(d))), F(1).name));
    todays_no2 = h5read(hinfo.Filename, h5dsetname(hinfo,1,2,1,1,dataset_name));
    todays_weight = h5read(hinfo.Filename, h5dsetname(hinfo,1,2,1,1,'Weight'));
    
    % Ignore negative values - they are either fill values or unphysical
    todays_weight(todays_no2<0) = 0;
    
    weighted_no2 = weighted_no2 + todays_no2 .* todays_weight;
    weights = weights + todays_weight;
    
end

% Finish the weighted averaging.
no2 = weighted_no2 ./ weights;

end

