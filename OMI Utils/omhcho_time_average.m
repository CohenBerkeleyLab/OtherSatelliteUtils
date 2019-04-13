function [ values, lon_grid, lat_grid, weights, stddev ] = omhcho_time_average(start_dates, end_dates, varargin)
%OMHCHO_TIME_AVERAGE Grid and average OMI HCHO data in time
%   [ VALUES, LON_GRID, LAT_GRID ] = OMHCHO_TIME_AVERAGE( START_DATES,
%   END_DATES ) Loads OMI HCHO data between START_DATES and END_DATES.
%   These can be date strings, date numbers, or cell arrays of either. If
%   given as cell arrays, both cell arrays must have the same number of
%   elements. Passing as cell arrays allows you to average noncontiguous
%   time periods. For example:
%
%       OMHCHO_TIME_AVERAGE( {'2012-04-01', '2013-04-01'}, {'2012-09-30',
%       '2013-09-30'} )
%
%   would create one average for all data between April 1st and Sept 30th
%   from both 2012 and 2013.
%
%   Parameters:
%
%       'domain' - specify the domain to grid the data to as a GlobeGrid
%       instance. Default is GlobeGrid(0.05, 'domain', 'us').
E = JLLErrors;
p = advInputParser;
p.addParameter('dayofweek', 'UMTWRFS');
p.addParameter('holidays', false);
p.addParameter('grid', GlobeGrid(0.05, 'domain', 'us'));
p.addParameter('uncertainty', 'weight');
p.addParameter('DEBUG_LEVEL',1);

p.parse(varargin{:});
pout = p.Results;

start_datenums = validate_date(start_dates);
end_datenums = validate_date(end_dates);

if numel(start_datenums) ~= numel(end_datenums)
    E.badinput('Number of start and end dates must be equal')
end

the_grid = pout.grid;
days_of_week = pout.dayofweek;

if (~islogical(pout.holidays) && ~isnumeric(pout.holidays)) || ~isscalar(pout.holidays)
    E.badinput('"holidays" must be a scalar logical or numeric value');
elseif pout.holidays
    holidays = [];
else
    holidays = 0;
end

uncertainty_mode = pout.uncertainty;

DEBUG_LEVEL = pout.DEBUG_LEVEL;

% Convert 1 letter day abbreviations to a vector that isbusday can
% understand
day_abbrevs = {'U', 'M', 'T', 'W', 'R', 'F', 'S'};
weekend = true(size(day_abbrevs));
for a=1:numel(day_abbrevs)
    weekend(a) = isempty(strfind(upper(days_of_week), day_abbrevs{a}));
end

values_avg = RunningAverage();
stddev_run = RunningAverage();
omhcho_files = [];
last_path = '';

for i_range = 1:numel(start_datenums)
    for dnum = start_datenums(i_range):end_datenums(i_range)
        if ~isbusday(dnum, holidays, weekend)
            if DEBUG_LEVEL > 0
                fprintf('Skipping %s due to day of week specification (%s)\n', datestr(dnum), days_of_week);
            end
            continue
        end
        fprintf('Now gridding %s\n', datestr(dnum));
        curr_files = get_file_for_date(dnum);
        for i_file = 1:numel(curr_files)
            this_file = fullfile(last_path, curr_files(i_file).name);
            xx_rows = in_domain(this_file);
            if ~any(xx_rows)
                continue
            end
            
            Data = load_omhcho_h5(this_file, xx_rows);
            Data.Weight = reject_bad_pix(Data);
            
            loncorn = convert_pix_corners(Data.PixelCornerLongitudes);
            latcorn = convert_pix_corners(Data.PixelCornerLatitudes);
            [gridded_val, gridded_wt] = cvm_generic_wrapper(loncorn, latcorn, Data.ReferenceSectorCorrectedVerticalColumn, the_grid, 'weights', Data.Weight);
            gridded_uncert = cvm_generic_wrapper(loncorn, latcorn, Data.ColumnUncertainty, the_grid, 'weights', Data.Weight);
            values_avg.addData(gridded_val, gridded_wt);
            stddev_run.addData(gridded_uncert.^2, gridded_wt);
        end
    end
end

values = values_avg.getWeightedAverage();
% I added the column uncertainties in quadrature, so we take
% the square root here to get back to the final value, divided
% by sqrt(n) to calculate the standard error of the mean
stddev = sqrt(stddev_run.values) ./ sqrt(stddev_run.weights); 
weights = values_avg.weights;
lon_grid = the_grid.GridLon;
lat_grid = the_grid.GridLat;

    function F = get_file_for_date(dnum)
        year_str = datestr(dnum, 'yyyy');
        month_str = datestr(dnum, 'mm');
        curr_path = fullfile(behr_paths.omno2_dir, '..', '..', 'OMHCHO', year_str, month_str);
        if ~strcmp(curr_path, last_path)
            % Use cached list of files if we can to speed things up
            omhcho_files = dir(fullfile(curr_path, '*.he5'));
        end
        last_path = curr_path;
        
        xx = regcmp({omhcho_files.name}, sprintf('^OMI-Aura_L2-OMHCHO_%sm%s%s', year_str, month_str, datestr(dnum, 'dd')));
        F = omhcho_files(xx);
    end
    
    function xx = in_domain(h5filename)
        hi = h5info(h5filename);
        lon = h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'Longitude'));
        lat = h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'Latitude'));
        is_in = lon >= the_grid.DomainLon(1) & lon <= the_grid.DomainLon(2) & lat >= the_grid.DomainLat(1) & lat <= the_grid.DomainLat(2);
        xx = any(is_in, 1);
    end

    function Data = load_omhcho_h5(h5filename, xx_rows)
        hi = h5info(h5filename);
        group1_fields = {'AMFCloudFraction','ReferenceSectorCorrectedVerticalColumn','ColumnUncertainty','MainDataQualityFlag','FittingRMS','PixelCornerLatitudes','PixelCornerLongitudes'};
        group2_fields = {'XtrackQualityFlags'};
        Data = make_empty_struct_from_cell([group1_fields, group2_fields]);
        
        xx_corner_rows = false(1,length(xx_rows)+1);
        xx_corner_rows(xx_rows) = true;
        xx_corner_rows(find(xx_rows,1,'last')+1) = true;
        
        for i = 1:numel(group1_fields)
            fname = group1_fields{i};
            Data.(fname) = h5readomi(hi.Filename, h5dsetname(hi,1,2,1,1, fname));
            if size(Data.(fname),2) == length(xx_rows)
                Data.(fname)(:,~xx_rows) = [];
            elseif size(Data.(fname),2) == length(xx_corner_rows)
                Data.(fname)(:,~xx_corner_rows) = [];
            else
                E.notimplemented('size(dataset,2) = %d', size(Data.(fname),2));
            end
        end
        
        for i = 1:numel(group2_fields)
            fname = group2_fields{i};
            Data.(fname) = h5readomi(hi.Filename, h5dsetname(hi,1,2,1,2, fname), 'keep_type', true);
            Data.(fname)(:,~xx_rows) = [];
        end
    end

    function weights = reject_bad_pix(Data)
        weights = zeros(size(Data.ReferenceSectorCorrectedVerticalColumn));
        xx = true(size(weights));
        xx = xx & Data.AMFCloudFraction <= 0.2;
        xx = xx & Data.MainDataQualityFlag == 0;
        xx = xx & bitand(Data.XtrackQualityFlags, int8(3)) == 0;
        if strcmpi(uncertainty_mode, 'cutoff')
            xx = xx & Data.FittingRMS < 0.02;
        end
        weights(xx) = 1;
        if strcmpi(uncertainty_mode, 'weight')
            weights = weights * 1 ./ Data.FittingRMS;
        end
    end
end



function corners = convert_pix_corners(tiled_corners)
corners = nan([4, size(tiled_corners)-1]);
for i = 1:size(tiled_corners,1)-1
    for j=1:size(tiled_corners,2)-1
        pix_corners = tiled_corners(i:i+1,j:j+1);
        pix_corners = pix_corners([1 2 4 3]); % this should put them in (counter)clockwise order
        corners(:,i,j) = pix_corners;
    end
end
end
