function OMI = hdf_quadrangle_GOME(Data, maxx, minx, maxy, miny, lCoordLon, lCoordLat, Lon1, Lon2, Lon4, Lat1, Lat2, Lat4)

% hdf_quadrangle_general: Version of hdf_quadrangle_5km_new by Ashley
% Russell from 11/19/2009 that serves as the template for oversampled
% grids.
%
%  ### This function is used by the Python script "hdf_quad_gen.py" to
%  generate final functions with all the fields filled in.  The python
%  script requires this template, plus a list a fields (each field on its
%  own line.  Any lines with the first non-whitespace characters %$f will
%  be replicated, replacing $field with every listed field. Those starting
%  with %$c will use the cell fields file, replacing $cellfield. Lines
%  beginning with %# will only have $keyfield replaced with the first
%  listed field and will not be replicated.  Be sure to rename the function
%  in the generated file.
%
%  Requires functions exchange_coord, calcline, clip.
%
%   This function should not be called directly; it is intended to be
%   called only from add2grid_general, which prepares geographic data for
%   this function.  It oversamples satellite data, mapping pixel
%   information to a fixed grid of smaller pixels.  These grids are
%   returned as the structure "OMI" (so named because it was first written
%   for NO2 data from the OMI instrument onboard the Aura satellite).  The
%   structure can be used in the no2_column_map_2014 function to produce
%   maps.
%
%   Josh Laughner <joshlaugh5@gmail.com> 18 Jul 2014


% Prepare empty matrices to receive the information that will make up the
% oversampled grid.  Those expected to receive quality flags will be cell
% arrays instead.

Time=zeros(maxx,maxy);
SolarZenithAngle=zeros(maxx,maxy);
ViewingZenithAngle=zeros(maxx,maxy);
RelativeAzimuthAngle=zeros(maxx,maxy);
NO2ColumnAmount=zeros(maxx,maxy);
NO2ColumnAbsError=zeros(maxx,maxy);
NO2ColumnAbsError_NoProfError=zeros(maxx,maxy);
NO2ColumnAmountTrop=zeros(maxx,maxy);
NO2ColumnTropAbsError=zeros(maxx,maxy);
NO2ColumnTropAbsError_NoProfError=zeros(maxx,maxy);
NO2ColumnAmountStrat=zeros(maxx,maxy);
NO2ColumStratAbsError=zeros(maxx,maxy);
SurfacePressure=zeros(maxx,maxy);
AvgKernel=zeros(maxx,maxy);
GhostColumn=zeros(maxx,maxy);
NO2SlantColumn=zeros(maxx,maxy);
CloudFraction=zeros(maxx,maxy);
CloudPressure=zeros(maxx,maxy);
CloudRadianceFraction=zeros(maxx,maxy);
AMF=zeros(maxx,maxy);
AMFTrop=zeros(maxx,maxy);
AMFGeometric=zeros(maxx,maxy);
Tropopause=zeros(maxx,maxy);
Albedo=zeros(maxx,maxy);
Count = zeros(maxx, maxy);
Area = nan(maxx, maxy);
Areaweight = nan(maxx, maxy);

TropFlag=cell(maxx,maxy);


%JLL 2-14-2014: Loads all the relevant fields from Data (the file
%loaded from reading the OMI_SP file)
% Time_i=Data(d).Time;
% ViewingZenithAngle_i=Data(d).ViewingZenithAngle;
% SolarZenithAngle_i=Data(d).SolarZenithAngle;
% ViewingAzimuthAngle_i=Data(d).ViewingAzimuthAngle;
% SolarAzimuthAngle_i=Data(d).SolarAzimuthAngle;
% CloudFraction_i=Data(d).CloudFraction;
% CloudRadianceFraction_i=Data(d).CloudRadianceFraction;
% ColumnAmountNO2_i=Data(d).ColumnAmountNO2;
% SlantColumnAmountNO2_i=Data(d).SlantColumnAmountNO2;
% ColumnAmountNO2Trop_i=Data(d).ColumnAmountNO2Trop;
% TerrainHeight_i=Data(d).TerrainHeight;
% TerrainPressure_i=Data(d).TerrainPressure;
% TerrainReflectivity_i=Data(d).TerrainReflectivity;
% CloudPressure_i=Data(d).CloudPressure;
% RelativeAzimuthAngle_i=Data(d).RelativeAzimuthAngle;
% Latitude_i=Data(d).Latitude;
% Longitude_i=Data(d).Longitude;
% GLOBETerpres_i=Data(d).GLOBETerpres;
% MODISAlbedo_i=Data(d).MODISAlbedo;
% BEHRAMFTrop_i=Data(d).BEHRAMFTrop;
% BEHRColumnAmountNO2Trop_i=Data(d).BEHRColumnAmountNO2Trop;
% MODISCloud_i=Data(d).MODISCloud;
% Row_i=Data(d).Row;
% Swath_i=Data(d).Swath;
% AMFTrop_i=Data(d).AMFTrop;
% AMFStrat_i=Data(d).AMFStrat;
% vcdQualityFlags_i=Data(d).vcdQualityFlags;
% XTrackQualityFlags_i=Data(d).XTrackQualityFlags;
% TropopausePressure_i=Data(d).TropopausePressure;
%
% Data.Count = zeros(size(Data.Latitude));
% Data.Area = nan(size(Data.Latitude));
% Data.Areaweight = nan(size(Data.Latitude));


Dimensions = size(Data.Time);
for x=1:1:Dimensions(1)*Dimensions(2); %JLL 18 Mar 2014: Loop over each NO2 column in Data(d)
    y=1;
    
    pixelarea=(m_lldist([Lon1(x,y)-180 Lon2(x,y)-180],[Lat1(x,y) Lat2(x,y)]))*(m_lldist([Lon1(x,y)-180, Lon4(x,y)-180],[Lat1(x,y), Lat4(x,y)])); %JLL 20 Mar 2014: This calculates the area of the 50% pixel response area in km. (Removed /1000 b/c this function should return in km already)
    
    %JLL 18 Mar 2014: lCoordLat/Lon are defined in add2grid_5km_new, they
    %are the lat/lon values (here only the corners are used) as multiples
    %of the resolution away from the lat/lon minimum.
    x1=round(lCoordLat(x,y,1)); x1=clip(x1,minx,maxx);
    y1=round(lCoordLon(x,y,1)); y1=clip(y1,miny,maxy);
    x2=round(lCoordLat(x,y,2)); x2=clip(x2,minx,maxx);
    y2=round(lCoordLon(x,y,2)); y2=clip(y2,miny,maxy);
    x3=round(lCoordLat(x,y,3)); x3=clip(x3,minx,maxx);
    y3=round(lCoordLon(x,y,3)); y3=clip(y3,miny,maxy);
    x4=round(lCoordLat(x,y,4)); x4=clip(x4,minx,maxx);
    y4=round(lCoordLon(x,y,4)); y4=clip(y4,miny,maxy);
    
    
    if y2<y1; [x1,y1,x2,y2]=exchange_coord(x1,y1,x2,y2); end
    if y3<y1; [x1,y1,x3,y3]=exchange_coord(x1,y1,x3,y3); end
    if y4<y1; [x1,y1,x4,y4]=exchange_coord(x1,y1,x4,y4); end
    if y2>y3; [x2,y2,x3,y3]=exchange_coord(x2,y2,x3,y3); end
    if y4>y3; [x4,y4,x3,y3]=exchange_coord(x4,y4,x3,y3); end
    if x4>x2; [x4,y4,x2,y2]=exchange_coord(x4,y4,x2,y2); end
    %JLL 2-14-2014: x1 to x4 and y1 to y4 are now integers, between their
    %respective min and max values (derived from the lat/lon bound
    %specified in the main BEHR file) and arranged so that the points are in order
    %going around the outside (i.e., pt. 2 will not be caddycorner to pt. 1). Further,
    %the points usually end up counterclockwise, with 1 as the bottom
    %point.
    
    %JLL 2-14-2014: Load, in turn, each value of the fields loaded from the
    %OMI standard product
    
    Time_val = Data.Time(x);
    SolarZenithAngle_val = Data.SolarZenithAngle(x);
    ViewingZenithAngle_val = Data.ViewingZenithAngle(x);
    RelativeAzimuthAngle_val = Data.RelativeAzimuthAngle(x);
    NO2ColumnAmount_val = Data.NO2ColumnAmount(x);
    NO2ColumnAbsError_val = Data.NO2ColumnAbsError(x);
    NO2ColumnAbsError_NoProfError_val = Data.NO2ColumnAbsError_NoProfError(x);
    NO2ColumnAmountTrop_val = Data.NO2ColumnAmountTrop(x);
    NO2ColumnTropAbsError_val = Data.NO2ColumnTropAbsError(x);
    NO2ColumnTropAbsError_NoProfError_val = Data.NO2ColumnTropAbsError_NoProfError(x);
    NO2ColumnAmountStrat_val = Data.NO2ColumnAmountStrat(x);
    NO2ColumStratAbsError_val = Data.NO2ColumStratAbsError(x);
    SurfacePressure_val = Data.SurfacePressure(x);
    AvgKernel_val = Data.AvgKernel(x);
    GhostColumn_val = Data.GhostColumn(x);
    NO2SlantColumn_val = Data.NO2SlantColumn(x);
    CloudFraction_val = Data.CloudFraction(x);
    CloudPressure_val = Data.CloudPressure(x);
    CloudRadianceFraction_val = Data.CloudRadianceFraction(x);
    AMF_val = Data.AMF(x);
    AMFTrop_val = Data.AMFTrop(x);
    AMFGeometric_val = Data.AMFGeometric(x);
    Tropopause_val = Data.Tropopause(x);
    Albedo_val = Data.Albedo(x);
    TropFlag_val = Data.TropFlag(x);
    
    
    %dim=[maxx maxy];
    bottom=y1+1; %JLL 18 Mar 2014: Having the bottom advance by one ensures that data points right on the bottom/top border don't get double counted (at least, I think that's the point here)
    top=y3;
    if (bottom<maxy) && (top>=1); %JLL 2-14-2014: Why are these not separate conditions, i.e. "if bottom < maxy; bottom = clip(...); end; if top >= 1..."
        bottom=clip(bottom,1,maxy);
        top=clip(top,1,maxy);
    end
    
    for y_quad=bottom:1:top; %JLL 19 Mar 2014:
        if (x2>=calcline(y2,x1,y1,x3,y3)) && (x4<=calcline(y4,x1,y1,x3,y3)); %JLL 18 Mar 2014: Tests if the points are arranged counterclockwise
            if y_quad<y4; %JLL 19 Mar 2014: y1 and y3 will be the bottom and top point of the quadrangle, so y2 and y4 are vertices on the sides
                left=calcline(y_quad,x1,y1,x4,y4); %JLL 19 Mar 2014: Use the line between y1 and y4 to calc the left side for the given row...
            else
                left=calcline(y_quad,x4,y4,x3,y3); %JLL 19 Mar 2014: ...unless the row is above y4, then use the y4-y3 line
            end
            if y_quad<y2;
                right=calcline(y_quad,x1,y1,x2,y2); %JLL 19 Mar 2014: Same thought process for the right
            else
                right=calcline(y_quad,x2,y2,x3,y3);
            end
        else %JLL 19 Mar 2014: This section *should* handle any cases in which the corners are not arranged counterclockwise
            left=calcline(y_quad,x1,y1,x3,y3);
            if y2>y4;
                [x4,y4,x2,y2]=exchange_coord(x4,y4,x2,y2);
            end
            if y_quad<y2;
                right=calcline(y_quad,x1,y1,x2,y2);
            elseif y_quad<y4;
                right=calcline(y_quad,x2,y2,x4,y4);
            else
                right=calcline(y_quad,x4,y4,x3,y3);
            end
            if (x2<=calcline(y2,x1,y1,x3,y3)) && (x4<=calcline(y4,x1,y1,x3,y3));
                placeholder=left;
                left=right;
                right=placeholder;
                clear placeholder
            end
        end
        right=right-1; %JLL 19 Mar 2014: Like with the bottom, decrement this by 1 to avoid double counting points
        if left<=0;
        elseif (maxx>=right) && (right>=1) && (1<=left) && (left<maxx); %JLL 19 Mar 2014: Make sure the left and right bounds are inside the permissible limits
            clip(left,1,maxx); %JLL 19 Mar 2014: Kind of redundant...
            clip(right,1,maxx);
            for x_quad=left+1:right;
                % JLL 18 Jul 2014: If there is a valid trace gas column,
                % add it to the grid.  If there is already a measurement
                % there, average the new and existing measurements,
                % weighted by the number of measurments that went into the
                % existing average
                if NO2ColumnAmountTrop(x_quad, y_quad) ~= 0 && ~isnan(Data.NO2ColumnAmountTrop(x));
                    % Count, area, and areaweight require special handling
                    Count(x_quad,y_quad)=Count(x_quad,y_quad)+1;
                    Area(x_quad,y_quad)=nansum([Area(x_quad,y_quad)*(Count(x_quad,y_quad)-1), pixelarea])/(Count(x_quad, y_quad));
                    Areaweight(x_quad,y_quad)=Count(x_quad, y_quad)/Area(x_quad,y_quad);
                    
                    % Regular fields will be a running average
                    Time(x_quad, y_quad) = sum([Time(x_quad, y_quad)*(Count(x_quad, y_quad)-1), Time_val])/Count(x_quad,y_quad);
                    SolarZenithAngle(x_quad, y_quad) = sum([SolarZenithAngle(x_quad, y_quad)*(Count(x_quad, y_quad)-1), SolarZenithAngle_val])/Count(x_quad,y_quad);
                    ViewingZenithAngle(x_quad, y_quad) = sum([ViewingZenithAngle(x_quad, y_quad)*(Count(x_quad, y_quad)-1), ViewingZenithAngle_val])/Count(x_quad,y_quad);
                    RelativeAzimuthAngle(x_quad, y_quad) = sum([RelativeAzimuthAngle(x_quad, y_quad)*(Count(x_quad, y_quad)-1), RelativeAzimuthAngle_val])/Count(x_quad,y_quad);
                    NO2ColumnAmount(x_quad, y_quad) = sum([NO2ColumnAmount(x_quad, y_quad)*(Count(x_quad, y_quad)-1), NO2ColumnAmount_val])/Count(x_quad,y_quad);
                    NO2ColumnAbsError(x_quad, y_quad) = sum([NO2ColumnAbsError(x_quad, y_quad)*(Count(x_quad, y_quad)-1), NO2ColumnAbsError_val])/Count(x_quad,y_quad);
                    NO2ColumnAbsError_NoProfError(x_quad, y_quad) = sum([NO2ColumnAbsError_NoProfError(x_quad, y_quad)*(Count(x_quad, y_quad)-1), NO2ColumnAbsError_NoProfError_val])/Count(x_quad,y_quad);
                    NO2ColumnAmountTrop(x_quad, y_quad) = sum([NO2ColumnAmountTrop(x_quad, y_quad)*(Count(x_quad, y_quad)-1), NO2ColumnAmountTrop_val])/Count(x_quad,y_quad);
                    NO2ColumnTropAbsError(x_quad, y_quad) = sum([NO2ColumnTropAbsError(x_quad, y_quad)*(Count(x_quad, y_quad)-1), NO2ColumnTropAbsError_val])/Count(x_quad,y_quad);
                    NO2ColumnTropAbsError_NoProfError(x_quad, y_quad) = sum([NO2ColumnTropAbsError_NoProfError(x_quad, y_quad)*(Count(x_quad, y_quad)-1), NO2ColumnTropAbsError_NoProfError_val])/Count(x_quad,y_quad);
                    NO2ColumnAmountStrat(x_quad, y_quad) = sum([NO2ColumnAmountStrat(x_quad, y_quad)*(Count(x_quad, y_quad)-1), NO2ColumnAmountStrat_val])/Count(x_quad,y_quad);
                    NO2ColumStratAbsError(x_quad, y_quad) = sum([NO2ColumStratAbsError(x_quad, y_quad)*(Count(x_quad, y_quad)-1), NO2ColumStratAbsError_val])/Count(x_quad,y_quad);
                    SurfacePressure(x_quad, y_quad) = sum([SurfacePressure(x_quad, y_quad)*(Count(x_quad, y_quad)-1), SurfacePressure_val])/Count(x_quad,y_quad);
                    AvgKernel(x_quad, y_quad) = sum([AvgKernel(x_quad, y_quad)*(Count(x_quad, y_quad)-1), AvgKernel_val])/Count(x_quad,y_quad);
                    GhostColumn(x_quad, y_quad) = sum([GhostColumn(x_quad, y_quad)*(Count(x_quad, y_quad)-1), GhostColumn_val])/Count(x_quad,y_quad);
                    NO2SlantColumn(x_quad, y_quad) = sum([NO2SlantColumn(x_quad, y_quad)*(Count(x_quad, y_quad)-1), NO2SlantColumn_val])/Count(x_quad,y_quad);
                    CloudFraction(x_quad, y_quad) = sum([CloudFraction(x_quad, y_quad)*(Count(x_quad, y_quad)-1), CloudFraction_val])/Count(x_quad,y_quad);
                    CloudPressure(x_quad, y_quad) = sum([CloudPressure(x_quad, y_quad)*(Count(x_quad, y_quad)-1), CloudPressure_val])/Count(x_quad,y_quad);
                    CloudRadianceFraction(x_quad, y_quad) = sum([CloudRadianceFraction(x_quad, y_quad)*(Count(x_quad, y_quad)-1), CloudRadianceFraction_val])/Count(x_quad,y_quad);
                    AMF(x_quad, y_quad) = sum([AMF(x_quad, y_quad)*(Count(x_quad, y_quad)-1), AMF_val])/Count(x_quad,y_quad);
                    AMFTrop(x_quad, y_quad) = sum([AMFTrop(x_quad, y_quad)*(Count(x_quad, y_quad)-1), AMFTrop_val])/Count(x_quad,y_quad);
                    AMFGeometric(x_quad, y_quad) = sum([AMFGeometric(x_quad, y_quad)*(Count(x_quad, y_quad)-1), AMFGeometric_val])/Count(x_quad,y_quad);
                    Tropopause(x_quad, y_quad) = sum([Tropopause(x_quad, y_quad)*(Count(x_quad, y_quad)-1), Tropopause_val])/Count(x_quad,y_quad);
                    Albedo(x_quad, y_quad) = sum([Albedo(x_quad, y_quad)*(Count(x_quad, y_quad)-1), Albedo_val])/Count(x_quad,y_quad);
                    
                    % Flag fields will append the flag value to a matrix in
                    % a cell corresponding to this grid cell
                    
                    TropFlag(x_quad, y_quad) = {[TropFlag{x_quad, y_quad}, TropFlag_val]};
                    
                    % If there is no existing field
                elseif ~isnan(Data.NO2ColumnAmountTrop(x)) %JLL 19 Mar 2014: I added the logical test here, before this was just an 'else' statement, but it would make sense not to add a value if there was no valid NO2 column.
                    % Count, area, and areaweight require special handling
                    Count(x_quad,y_quad)=Count(x_quad,y_quad)+1;
                    Area(x_quad,y_quad)=pixelarea;
                    Areaweight(x_quad,y_quad)=1/pixelarea;
                    
                    % Regular fields will be a running average
                    Time(x_quad, y_quad) = Time_val;
                    SolarZenithAngle(x_quad, y_quad) = SolarZenithAngle_val;
                    ViewingZenithAngle(x_quad, y_quad) = ViewingZenithAngle_val;
                    RelativeAzimuthAngle(x_quad, y_quad) = RelativeAzimuthAngle_val;
                    NO2ColumnAmount(x_quad, y_quad) = NO2ColumnAmount_val;
                    NO2ColumnAbsError(x_quad, y_quad) = NO2ColumnAbsError_val;
                    NO2ColumnAbsError_NoProfError(x_quad, y_quad) = NO2ColumnAbsError_NoProfError_val;
                    NO2ColumnAmountTrop(x_quad, y_quad) = NO2ColumnAmountTrop_val;
                    NO2ColumnTropAbsError(x_quad, y_quad) = NO2ColumnTropAbsError_val;
                    NO2ColumnTropAbsError_NoProfError(x_quad, y_quad) = NO2ColumnTropAbsError_NoProfError_val;
                    NO2ColumnAmountStrat(x_quad, y_quad) = NO2ColumnAmountStrat_val;
                    NO2ColumStratAbsError(x_quad, y_quad) = NO2ColumStratAbsError_val;
                    SurfacePressure(x_quad, y_quad) = SurfacePressure_val;
                    AvgKernel(x_quad, y_quad) = AvgKernel_val;
                    GhostColumn(x_quad, y_quad) = GhostColumn_val;
                    NO2SlantColumn(x_quad, y_quad) = NO2SlantColumn_val;
                    CloudFraction(x_quad, y_quad) = CloudFraction_val;
                    CloudPressure(x_quad, y_quad) = CloudPressure_val;
                    CloudRadianceFraction(x_quad, y_quad) = CloudRadianceFraction_val;
                    AMF(x_quad, y_quad) = AMF_val;
                    AMFTrop(x_quad, y_quad) = AMFTrop_val;
                    AMFGeometric(x_quad, y_quad) = AMFGeometric_val;
                    Tropopause(x_quad, y_quad) = Tropopause_val;
                    Albedo(x_quad, y_quad) = Albedo_val;
                    
                    % Flag fields will append the flag value to a matrix in
                    % a cell corresponding to this grid cell
                    TropFlag(x_quad, y_quad) = {TropFlag_val};
                end
            end
        end
    end
end

% Create the OMI structure for output
OMI.Time = Time;
OMI.SolarZenithAngle = SolarZenithAngle;
OMI.ViewingZenithAngle = ViewingZenithAngle;
OMI.RelativeAzimuthAngle = RelativeAzimuthAngle;
OMI.NO2ColumnAmount = NO2ColumnAmount;
OMI.NO2ColumnAbsError = NO2ColumnAbsError;
OMI.NO2ColumnAbsError_NoProfError = NO2ColumnAbsError_NoProfError;
OMI.NO2ColumnAmountTrop = NO2ColumnAmountTrop;
OMI.NO2ColumnTropAbsError = NO2ColumnTropAbsError;
OMI.NO2ColumnTropAbsError_NoProfError = NO2ColumnTropAbsError_NoProfError;
OMI.NO2ColumnAmountStrat = NO2ColumnAmountStrat;
OMI.NO2ColumStratAbsError = NO2ColumStratAbsError;
OMI.SurfacePressure = SurfacePressure;
OMI.AvgKernel = AvgKernel;
OMI.GhostColumn = GhostColumn;
OMI.NO2SlantColumn = NO2SlantColumn;
OMI.CloudFraction = CloudFraction;
OMI.CloudPressure = CloudPressure;
OMI.CloudRadianceFraction = CloudRadianceFraction;
OMI.AMF = AMF;
OMI.AMFTrop = AMFTrop;
OMI.AMFGeometric = AMFGeometric;
OMI.Tropopause = Tropopause;
OMI.Albedo = Albedo;
OMI.Count = Count;
OMI.Area = Area;
OMI.Areaweight = Areaweight;
OMI.TropFlag = TropFlag;

end
