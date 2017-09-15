function [ residuals, lon, lat ] = calc_sw_residuals( filepath, startdate, enddate, varargin )
%RESIDUALS = CALC_SW_RESIDUALS(FILEPATH, STARTDATE, ENDDATE)
%   This function is a one-off script to calculate residual differences
%   between BEHR and OMNO2 scattering weights. FILEPATH specifies the path
%   to the OMI_BEHR and OMI_SP files with the scattering weights you want
%   to compare. STARTDATE and ENDDATE specify which days to compare. The
%   output, RESIDUALS will be the sum of squares of the differences in
%   scattering weights after interpolation to the OMNO2 pressure levels. It
%   will output in a matrix of swaths concatenated along the along track
%   dimension.
%
%   Additional arguments are 'land' or 'ocean' which will restrict the
%   comparison to over land or ocean only. In this case, residuals in
%   pixels that do not meet the criteria are set to nans.

E = JLLErrors;

if ~ischar(filepath) || ~exist(filepath,'dir')
    E.badinput('FILEPATH must be a directory')
end

% First we find the BEHR files in the given directory between the dates
% specified.

behrnames = files_in_dates(filepath, 'OMI_BEHR*.mat', startdate, enddate);

% Load the first BEHR file, see if it has the fields ScatteringWeight and
% ScatteringWtPressure, if not, try finding all the OMI_SP files and see if
% they're in there.

B = load(fullfile(filepath, behrnames{1}), 'Data');
get_sw_in_sp = false;
if any(~isfield(B.Data,{'ScatteringWeight','ScatteringWtPressure'}))
    fprintf('OMNO2 scattering weights not in BEHR files, checking for SP files in same directory..\n')
    spnames = files_in_dates(filepath, 'OMI_SP*.mat', startdate, enddate);
    if isempty(spnames)
        E.callError('no_sw','OMNO2 scattering weights not present in BEHR files and no SP files included in same directory');
    elseif numel(spnames) ~= numel(behrnames)
        E.callError('no_sw','Different number of BEHR and SP files; must have same number to match up BEHR and OMNO2 scattering weights between the two.')
    end
    S = load(fullfile(filepath, spnames{1}));
    if any(~isfield(S.Data, {'ScatteringWeight','ScatteringWtPressure'}))
        E.callError('no_sw','OMNO2 scattering weights not found in either BEHR or SP files');
    end
    fprintf('Weights found in SP files.\n')
    get_sw_in_sp = true;
    omno2_sw_pres = S.Data(1).ScatteringWtPressure;
else
    omno2_sw_pres = B.Data(1).ScatteringWtPressure;
end
trop_sw = omno2_sw_pres >= 200;
omno2_sw_pres(~trop_sw) = [];

% Now load each file. Get the BEHR scattering weights and interpolate to
% the OMNO2 pressures. Calculate the sum of squared residuals for each
% pixel and report this. Concatenate across days, and also concatenate
% pixel lat/lon to allow mapping if desired.

residuals = [];
lon = [];
lat = [];

for a=1:numel(behrnames)
    B = load(fullfile(filepath, behrnames{a}), 'Data');
    if get_sw_in_sp
        S = load(fullfile(filepath, spnames{a}));
        if numel(B.Data) ~= numel(S.Data)
            fprintf('%s has different number of swaths than %s, skipping\n', behrnames{a}, spnames{a});
            continue
        end
    end
    for b=1:numel(B.Data)
        lon = cat(1,lon, B.Data(b).Longitude);
        lat = cat(1,lat, B.Data(b).Latitude);
        this_behrsw = B.Data(b).BEHRScatteringWeights;
        this_behrswpres = B.Data(b).BEHRPressureLevels;
        this_resid = nan(size(B.Data(b).Longitude));
        if ~get_sw_in_sp
            this_omno2sw = B.Data(b).ScatteringWeight;
        else
            this_omno2sw = S.Data(b).ScatteringWeight;
        end
        for c=1:numel(this_resid)
            behr_sw_vec = this_behrsw(:,c);
            behr_swp_vec = this_behrswpres(:,c);
            nans = isnan(behr_sw_vec);
            behr_sw_vec(nans) = [];
            behr_swp_vec(nans) = [];
            if numel(behr_sw_vec) > 1
                interp_behrsw = interp1(behr_swp_vec, behr_sw_vec, omno2_sw_pres);
                this_resid(c) = nansum2((interp_behrsw - this_omno2sw(trop_sw,c)).^2);
            end
        end
        residuals = cat(1, residuals, this_resid);
    end
end

% Finally NaN out any points that aren't over ocean or land, but only if
% one has been requested.
if any(ismember({'ocean','land'},varargin))
    xx = island(lon,lat);
    if ismember('land',varargin)
        xx = ~xx;
    end
    residuals(xx) = nan;
end

end

function files = files_in_dates(filepath, pattern, startdate, enddate)
F = dir(fullfile(filepath, pattern));
files = {F.name}';
xx = false(size(files));
for a=1:numel(files)
    [s,e] = regexp(files{a},'\d\d\d\d\d\d\d\d');
    dstr = files{a}(s:e);
    if datenum(dstr,'yyyymmdd') >= datenum(startdate) && datenum(dstr,'yyyymmdd') <= datenum(enddate)
        xx(a) = true;
    end
end
files(~xx) = [];
end
