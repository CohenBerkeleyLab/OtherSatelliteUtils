function read_GOME2

% read_GOME2_v1
%
%   Script that will read GOME-2 NO2 data into a Matlab file containing a
%   structure "Data" where the top-level indicies represent different
%   swaths for a that day

date_start = '2/05/2012';
date_end = '12/31/2012';

latbdy = [-90 90];
lonbdy = [-180 180];

GOME_dir = '/Volumes/share-sat/SAT/GOME-2/HDF Files/';
mat_dir = '/Volumes/share-sat/SAT/GOME-2/Matlab Files/Pixels/Global/';

GOME_prefix = 'no2track';
save_prefix = 'GOME2_';

DEBUG_LEVEL = 2;

datenums = datenum(date_start):datenum(date_end);
for d=1:numel(datenums);
    curr_date = datestr(datenums(d),29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    
    D=0; % Keeps track of what index of Data we're entering
    
    fname = sprintf('%s%s%s%s.hdf',GOME_prefix,year,month,day);
    filepath = fullfile(GOME_dir,year,fname);
    if ~exist(filepath,'file');
        if DEBUG_LEVEL > 0; fprintf('No data available for %s\n',datestr(curr_date)); end
        continue
    else
        if DEBUG_LEVEL > 0; fprintf('Now on %s\n',curr_date); end
        
        % If the file exists, get its structure
        hdfi = hdfinfo(filepath);
        
        Data = struct('Date',datestr(curr_date,29),'Time',[],'Latitude',[],'Longitude',[],'Latcorn',[],'Loncorn',[],'SolarZenithAngle',[],'ViewingZenithAngle',[],...
            'RelativeAzimuthAngle',[],'NO2ColumnAmount',[],'NO2ColumnAbsError',[],'NO2ColumnAbsError_NoProfError',[],'NO2ColumnAmountTrop',[],...
            'NO2ColumnTropAbsError',[],'NO2ColumnTropAbsError_NoProfError',[],'NO2ColumnAmountStrat',[],'NO2ColumStratAbsError',[],'TropFlag',[],...
            'SurfacePressure',[],'AvgKernel',[],'GhostColumn',[],'NO2SlantColumn',[],'CloudFraction',[],'CloudPressure',[],'CloudRadianceFraction',[],...
            'AMF',[],'AMFTrop',[],'AMFGeometric',[],'Tropopause',[],'Albedo',[]);
        
        % Loop through the field names, looking for ones that have "NO2" in
        % the name. Since the first field is pressure levels, we'll skip
        % that one right off.
        n = length(hdfi.Vdata);
        for a=2:n
            fieldname = hdfi.Vdata(a).Name;
            if isempty(regexp(fieldname,'NO2','ONCE'));
                continue
            else
                
                
                no2cell = hdfread(filepath,hdfi.Vdata(a).Name);
                Latitude = no2cell{4}; 
                Longitude = no2cell{3}; 
                Longitude(Longitude > 180) = Longitude(Longitude > 180) - 360;
                
                xx = Latitude >= min(latbdy) & Latitude <= max(latbdy) & Longitude >= min(lonbdy) & Longitude <= max(lonbdy);
                
                % If no pixels fall within our lat/lon boundary, skip this
                % field
                if sum(xx)==0
                    continue
                end
                
                % As long as pixels have been found, load the other two
                % fields that refer to this swath.                
                geocell = hdfread(filepath,hdfi.Vdata(a+1).Name);
                anccell = hdfread(filepath, hdfi.Vdata(a+2).Name);
                
                % Save all the fields we care about in the data structure,
                % removing points outside the lat/lon boundary.  Also,
                % check fields that shouldn't be negative for negative
                % values, these indicate likely fills
                
                D=D+1;
                Data(D).Date = datestr(curr_date,29);
                Data(D).Time = double(no2cell{2}(xx));
                Data(D).Latitude = Latitude(xx);
                Data(D).Longitude = Longitude(xx);
                Data(D).Latcorn = double(geocell{6}(:,xx));
                Data(D).Loncorn = double(geocell{5}(:,xx))-360;
                Data(D).SolarZenithAngle = double(geocell{1}(xx));
                Data(D).ViewingZenithAngle = double(geocell{2}(xx));
                Data(D).RelativeAzimuthAngle = double(geocell{3}(xx));
                Data(D).NO2ColumnAmount = double(no2cell{5}(xx))*1e15;       Data(D).NO2ColumnAmount(Data(D).NO2ColumnAmount<0) = NaN;
                Data(D).NO2ColumnAbsError = double(no2cell{6}(xx))*1e15;
                Data(D).NO2ColumnAbsError_NoProfError = double(no2cell{13}(xx))*1e15;
                Data(D).NO2ColumnAmountTrop = double(no2cell{7}(xx))*1e15;   Data(D).NO2ColumnAmountTrop(Data(D).NO2ColumnAmountTrop<0) = NaN;
                Data(D).NO2ColumnTropAbsError = double(no2cell{8}(xx))*1e15;
                Data(D).NO2ColumnTropAbsError_NoProfError = double(no2cell{14}(xx))*1e15;
                Data(D).NO2ColumnAmountStrat = double(no2cell{9}(xx))*1e15;  Data(D).NO2ColumnAmountStrat(Data(D).NO2ColumnAmountStrat<0) = NaN;
                Data(D).NO2ColumStratAbsError = double(no2cell{10}(xx))*1e15;
                Data(D).TropFlag = no2cell{11}(xx);
                Data(D).SurfacePressure = no2cell{12}(xx);
                Data(D).AvgKernel = (double(no2cell{15}(:,xx)));
                Data(D).GhostColumn = double(no2cell{16}(xx))*1e15;          Data(D).GhostColumn(Data(D).GhostColumn<0) = NaN;
                Data(D).NO2SlantColumn = double(anccell{1}(xx))*1e15;        Data(D).NO2SlantColumn(Data(D).NO2SlantColumn<0) = NaN;
                Data(D).CloudFraction = double(anccell{6}(xx));         Data(D).CloudFraction(Data(D).CloudFraction<0) = NaN;
                Data(D).CloudPressure = double(anccell{7}(xx));         Data(D).CloudPressure(Data(D).CloudPressure<0) = NaN;
                % Cloud radiance fraction seems to be given as a percent,
                % e.g. 100 = full cloud fraction.  Geometric cloud fraction
                % seems to be given as a decimal.
                Data(D).CloudRadianceFraction = double(anccell{9}(xx))/100; Data(D).CloudRadianceFraction(Data(D).CloudRadianceFraction<0) = NaN;
                Data(D).AMF = double(anccell{2}(xx));                   Data(D).AMF(Data(D).AMF<0) = NaN;
                Data(D).AMFTrop = double(anccell{3}(xx));               Data(D).AMFTrop(Data(D).AMFTrop<0) = NaN;
                Data(D).AMFGeometric = double(anccell{4}(xx));          Data(D).AMFGeometric(Data(D).AMFGeometric<0) = NaN;
                Data(D).Tropopause = double(anccell{10}(xx));           Data(D).Tropopause(Data(D).Tropopause<0) = NaN;
                Data(D).Albedo = double(anccell{9}(xx));                Data(D).Albedo(Data(D).Albedo<0) = NaN;
                
            end
        end
    end
    
    % Save the data structure
    savename = sprintf('%s%s%s%s.mat',save_prefix,year,month,day);
    savepath = fullfile(mat_dir,year,savename);
    save(savepath,'Data');
    clear Data
end
end