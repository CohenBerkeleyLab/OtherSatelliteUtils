function [ varargout ] = plot_myd14( start_date, end_date, lonlim, latlim )
%PLOT_MYD14(start_date, end_date) Plot thermal anomalies between two days
%   Function will identify all thermal anomalies between the two days
%   specified and plot their locations and color them by FRP in MW.
%
%   If only one argument is given or the second argument is an empty
%   string, then only the day given in the first argument will be plotted.
%   Any date string format is acceptable.
%
%   You may also specify a manual longitude and latitude limit as inputs 3
%   and 4; these must be of class double and be 1-by-2 vectors
%
%   Josh Laughner <joshlaugh5@gmail.com> 30 Jul 2015

%%%%% INPUT CHECKING %%%%%
E=JLLErrors;

narginchk(1,4)
if nargin < 2 || isempty(end_date)
    end_date = start_date;
end

if nargin > 2 && (~isrow(lonlim) || any(size(lonlim) ~= [1 2]) || ~isa(lonlim,'double'))
    E.badinput('lonlim, if given, must be a 1-by-2 vector of doubles')
end
if nargin > 3 && (~isrow(latlim) || any(size(latlim) ~= [1 2]) || ~isa(latlim, 'double'))
    E.badinput('latlim, if given, must be a 1-by-2 vector of doubles')
end

%%%%% DEFINED FILE PATHS %%%%%
myd14_root = fullfile('/Volumes','share-sat','SAT','MODIS','MYD14');

%%%%% MAIN FUNCTION %%%%%
dates = datenum(start_date):datenum(end_date);
i=1;
for d=1:numel(dates)
    curr_year = sprintf('%04d',year(dates(d)));
    curr_month = sprintf('%02d',month(dates(d)));
    
    curr_path = fullfile(myd14_root, curr_year, curr_month);
    filepat = sprintf('MYD14.A%s%03d.*.hdf',curr_year,modis_date_to_day(dates(d)));
    
    F = dir(fullfile(curr_path, filepat));
    for a=2:numel(F)
        [~, frp{i}, lon{i}, lat{i}] = read_myd14(fullfile(curr_path,F(a).name));
    end
    i=i+1;
end

frp = [frp{:}];
lon = [lon{:}];
lat = [lat{:}];

if ~exist('lonlim','var')
    lonlim = double([floor(min(lon)), ceil(max(lon))]);
end
if ~exist('latlim','var')
    latlim = double([floor(min(lat)), ceil(max(lat))]);
end

fig = figure;
worldmap(latlim, lonlim);

scatterm(lat,lon,18,frp);
cb=colorbar;

set(get(cb,'Label'),'String','FRP (MW)');
state_outlines(fig,'all');

if nargout > 0
    varargout{1} = frp;
    varargout{2} = lon;
    varargout{3} = lat;
end

end

