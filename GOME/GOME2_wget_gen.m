function [ varargout ] = GOME2_wget_gen( filepath, startdate, enddate )
% GOME2_WGET_GEN(FILEPATH, STARTDATE, ENDDATE) 
%   Generates a list of http links for wget to retrieve GOME 2 files. Requires
%   three inputs: the file name to output (can be a full path), and the start 
%   and end dates.
%
%   STATUS = GOME2_wget_gen(___) will return the closing status of the file.

% Open a new file to write the links to
fid = fopen(filepath,'w');

% This is the format of the URL where the GOME2 NO2 files can be found.
% %1$s will be replaced by year, %2$s by month, and %3$s by day
gome_string = 'http://www.temis.nl/airpollution/no2col/data/gome2_v2/%1$s/%2$s/no2track%1$s%2$s%3$s.hdf\n';

dates = datenum(startdate):datenum(enddate);
for d=1:numel(dates)
    curr_date = datestr(dates(d),29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    
    % Output the link for each day to the file
    fprintf(fid,gome_string,year,month,day);
end

% Close the file
status = fclose(fid);
if nargout > 0
    varargout{1} = status;
end

end

