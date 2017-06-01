function [ Data ] = read_mls_data( filename )
%READ_MLS_DATA Returns a structure with filtered data from MLS files
%   The Microwave Limb Sounder (MLS) onboard the Aura satellite uses
%   naturally emitted microwaves to probe the concentration of various
%   species in the upper troposphere and stratosphere. There are a number
%   of screening parameters that are to be applied in order to filter off
%   poor quality data. These vary with the compound being retrieved.
%
%   This function will takes as its only argument a filename (with path if
%   necessary) to the MLS file you wish to load.  It will apply the proper
%   filtering for the given compound (detected automatically from the HDF
%   file structure) and return a Matlab structure with the matrix of
%   concentrations, the pressure level vector, and the lat/lon coordinates
%   of each measurement.
%
%   Josh Laughner <joshlaugh5@gmail.com> 10 Nov 2015

E = JLLErrors;
DEBUG_LEVEL = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(filename,'file')
    E.badinput('The file %s does not exist',filename)
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the HDF structure and determine which compound it is. Aura HDF-EOS
% files are structured as /HDFEOS/Swaths/<product> and in MLS it seems that
% the first product is always the main one, so we sneakily use fileparts to
% get just the last part of the swath name (after the last /)
hi = h5info(filename);
[~, cmpd] = fileparts(hi.Groups(1).Groups(2).Groups(1).Name);

% Load the spatial information
Data.Longitude = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'Longitude')));
Data.Latitude = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'Latitude')));
Data.Pressure = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'Pressure')));
Data.SolarZenithAngle = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'SolarZenithAngle')));
Data.OrbitGeodeticAngle = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'OrbitGeodeticAngle')));
Data.LineOfSightAngle = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'LineOfSightAngle')));
Data.LocalSolarTime = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'LocalSolarTime')));

% Returns the data filtered such that poor quality data is replaced with
% NaNs.
Data.(cmpd) = double(filter_meas(hi, cmpd, Data.Pressure, DEBUG_LEVEL));



end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SUBFUNCTIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

function data = filter_meas(hi, cmpd_name, pressure, DEBUG_LEVEL)
[data, precision, quality, status, convergence] = load_data(hi);
E = JLLErrors;

% Set the criteria based on which compound we're loading. Generally, there
% are at least 4 criteria to check: 
%   1) is precision negative
%   2) is quality above the threshold value
%   3) is status even (or is status 0, depending)
%   4) is convergence sufficiently close to unity
%   5) some have additional criteria. Use the "addn_matrix" variable to
%   describe this.
% Specifics can be found in the MLS Level 2 data quality and description
% document. Current version of MLS data is 4.2 as of 10 Nov 2015, so see
% https://mls.jpl.nasa.gov/data/v4-2_data_quality_document.pdf for
% specifics.
%
% Create matrices indicating which data points should be set to NaN.
switch cmpd_name
    case 'HNO3'
        precision_mat = precision < 0;
        quality_mat = quality <= 1.0; %1.0 was used in Miyazaki 2012 %0.8 is what's recommended in the data document;
        status_mat = status > 0;
        convergence_mat = convergence >= 1.03;
        fill_val_mat = data < -900; % Fill value is -999.99 (approximately)
        addn_matrix = hno3_addn_fxn(data,pressure);
    case 'O3'
        precision_mat = precision < 0;
        quality_mat = quality <= 1.0;
        status_mat = mod(status,2) > 0;
        convergence_mat = convergence >= 1.03;
        fill_val_mat = data < -900; % Fill value is -999.99 (approximately)
        addn_matrix = false(size(data));
    otherwise
        E.notimplemented('Filtering criteria for %s has not been set up',cmpd_name);
end

data(precision_mat | addn_matrix | fill_val_mat) = nan;
data(:, quality_mat | status_mat | convergence_mat) = nan;

end


function set_to_nan = hno3_addn_fxn(data, pressure)
% An additional filtering function for HNO3 looks at the volume mixing
% ratios at 316 - 215 hPa. The data quality document indicates that
% filtering off overly negative values at these pressure levels will
% remove profiles that result from poor retrieval fits resulting in
% oscillatory behavior. Because this affects the profile as a whole, these
% pressure levels are used to diagnose the whole profile.

p316 = abs(pressure - 316) < 1; % allow for floating point error
p215to68 = pressure <= 215 & pressure >= 68;

set_to_nan = false(size(data));

% Profiles load as (pressure level) x (ground coordinate), so iterate over
% each ground coordinate and check the pixel.
for a=1:size(data,2)
    if data(p316,a) < -2e-9 || any(data(p215to68,a)) < -1.6e-9
        set_to_nan(:,a) = true;
    end
end
end

function [data, precision, quality, status, convergence] = load_data(hi)
% Will load the datasets from the file. The structure should be the same
% regardless of the compound in question.
data = h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'L2gpValue'));
precision = h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'L2gpPrecision'));
quality = h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'Quality'));
status = h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'Status'));
convergence = h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'Convergence'));
end

