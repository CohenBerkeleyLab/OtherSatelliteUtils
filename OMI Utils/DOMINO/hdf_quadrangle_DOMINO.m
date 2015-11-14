function OMI = hdf_quadrangle_DOMINO(Data, OMI, maxx, minx, maxy, miny, lCoordLon, lCoordLat, Lon1, Lon2, Lon4, Lat1, Lat2, Lat4)

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

fill_val = -9e9;
AirMassFactorTropospheric=fill_val * ones(maxx,maxy);
AveragingKernel=fill_val * ones(maxx,maxy);
CloudFraction=fill_val * ones(maxx,maxy);
CloudPressure=fill_val * ones(maxx,maxy);
GhostColumn=fill_val * ones(maxx,maxy);
SlantColumnAmountNO2=fill_val * ones(maxx,maxy);
SlantColumnAmountNO2Std=fill_val * ones(maxx,maxy);
SurfaceAlbedo=fill_val * ones(maxx,maxy);
TerrainHeight=fill_val * ones(maxx,maxy);
TroposphericVerticalColumn=fill_val * ones(maxx,maxy);
TroposphericVerticalColumnError=fill_val * ones(maxx,maxy);
Count = zeros(maxx, maxy);
Area = nan(maxx, maxy);
Areaweight = nan(maxx, maxy);

AirMassFactorTropospheric=cell(maxx,maxy);
AveragingKernel=cell(maxx,maxy);
CloudFraction=cell(maxx,maxy);
CloudPressure=cell(maxx,maxy);
GhostColumn=cell(maxx,maxy);
SlantColumnAmountNO2=cell(maxx,maxy);
SlantColumnAmountNO2Std=cell(maxx,maxy);
SurfaceAlbedo=cell(maxx,maxy);
TerrainHeight=cell(maxx,maxy);
TroposphericVerticalColumn=cell(maxx,maxy);
TroposphericVerticalColumnError=cell(maxx,maxy);


 Dimensions = size(Data.AirMassFactorTropospheric);
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
    
    AirMassFactorTropospheric_val = Data.AirMassFactorTropospheric(x);
    AveragingKernel_val = Data.AveragingKernel(x);
    CloudFraction_val = Data.CloudFraction(x);
    CloudPressure_val = Data.CloudPressure(x);
    GhostColumn_val = Data.GhostColumn(x);
    SlantColumnAmountNO2_val = Data.SlantColumnAmountNO2(x);
    SlantColumnAmountNO2Std_val = Data.SlantColumnAmountNO2Std(x);
    SurfaceAlbedo_val = Data.SurfaceAlbedo(x);
    TerrainHeight_val = Data.TerrainHeight(x);
    TroposphericVerticalColumn_val = Data.TroposphericVerticalColumn(x);
    TroposphericVerticalColumnError_val = Data.TroposphericVerticalColumnError(x);
    AirMassFactorTropospheric_val = Data.AirMassFactorTropospheric(x);
    AveragingKernel_val = Data.AveragingKernel(x);
    CloudFraction_val = Data.CloudFraction(x);
    CloudPressure_val = Data.CloudPressure(x);
    GhostColumn_val = Data.GhostColumn(x);
    SlantColumnAmountNO2_val = Data.SlantColumnAmountNO2(x);
    SlantColumnAmountNO2Std_val = Data.SlantColumnAmountNO2Std(x);
    SurfaceAlbedo_val = Data.SurfaceAlbedo(x);
    TerrainHeight_val = Data.TerrainHeight(x);
    TroposphericVerticalColumn_val = Data.TroposphericVerticalColumn(x);
    TroposphericVerticalColumnError_val = Data.TroposphericVerticalColumnError(x);
    
    
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
                 if AirMassFactorTropospheric(x_quad, y_quad) ~= fill_val && ~isnan(AirMassFactorTropospheric_val);
                % Count, area, and areaweight require special handling
                Count(x_quad,y_quad)=Count(x_quad,y_quad)+1;
                Area(x_quad,y_quad)=nansum([Area(x_quad,y_quad)*(Count(x_quad,y_quad)-1), pixelarea])/(Count(x_quad, y_quad));
                Areaweight(x_quad,y_quad)=Count(x_quad, y_quad)/Area(x_quad,y_quad);
                
                % Regular fields will be a running average
                AirMassFactorTropospheric(x_quad, y_quad) = sum([AirMassFactorTropospheric(x_quad, y_quad)*(Count(x_quad, y_quad)-1), AirMassFactorTropospheric_val])/Count(x_quad,y_quad);
                AveragingKernel(x_quad, y_quad) = sum([AveragingKernel(x_quad, y_quad)*(Count(x_quad, y_quad)-1), AveragingKernel_val])/Count(x_quad,y_quad);
                CloudFraction(x_quad, y_quad) = sum([CloudFraction(x_quad, y_quad)*(Count(x_quad, y_quad)-1), CloudFraction_val])/Count(x_quad,y_quad);
                CloudPressure(x_quad, y_quad) = sum([CloudPressure(x_quad, y_quad)*(Count(x_quad, y_quad)-1), CloudPressure_val])/Count(x_quad,y_quad);
                GhostColumn(x_quad, y_quad) = sum([GhostColumn(x_quad, y_quad)*(Count(x_quad, y_quad)-1), GhostColumn_val])/Count(x_quad,y_quad);
                SlantColumnAmountNO2(x_quad, y_quad) = sum([SlantColumnAmountNO2(x_quad, y_quad)*(Count(x_quad, y_quad)-1), SlantColumnAmountNO2_val])/Count(x_quad,y_quad);
                SlantColumnAmountNO2Std(x_quad, y_quad) = sum([SlantColumnAmountNO2Std(x_quad, y_quad)*(Count(x_quad, y_quad)-1), SlantColumnAmountNO2Std_val])/Count(x_quad,y_quad);
                SurfaceAlbedo(x_quad, y_quad) = sum([SurfaceAlbedo(x_quad, y_quad)*(Count(x_quad, y_quad)-1), SurfaceAlbedo_val])/Count(x_quad,y_quad);
                TerrainHeight(x_quad, y_quad) = sum([TerrainHeight(x_quad, y_quad)*(Count(x_quad, y_quad)-1), TerrainHeight_val])/Count(x_quad,y_quad);
                TroposphericVerticalColumn(x_quad, y_quad) = sum([TroposphericVerticalColumn(x_quad, y_quad)*(Count(x_quad, y_quad)-1), TroposphericVerticalColumn_val])/Count(x_quad,y_quad);
                TroposphericVerticalColumnError(x_quad, y_quad) = sum([TroposphericVerticalColumnError(x_quad, y_quad)*(Count(x_quad, y_quad)-1), TroposphericVerticalColumnError_val])/Count(x_quad,y_quad);
                
                % Flag fields will append the flag value to a matrix in
                % a cell corresponding to this grid cell
                
                AirMassFactorTropospheric(x_quad, y_quad) = {[AirMassFactorTropospheric{x_quad, y_quad}, AirMassFactorTropospheric_val]};
                AveragingKernel(x_quad, y_quad) = {[AveragingKernel{x_quad, y_quad}, AveragingKernel_val]};
                CloudFraction(x_quad, y_quad) = {[CloudFraction{x_quad, y_quad}, CloudFraction_val]};
                CloudPressure(x_quad, y_quad) = {[CloudPressure{x_quad, y_quad}, CloudPressure_val]};
                GhostColumn(x_quad, y_quad) = {[GhostColumn{x_quad, y_quad}, GhostColumn_val]};
                SlantColumnAmountNO2(x_quad, y_quad) = {[SlantColumnAmountNO2{x_quad, y_quad}, SlantColumnAmountNO2_val]};
                SlantColumnAmountNO2Std(x_quad, y_quad) = {[SlantColumnAmountNO2Std{x_quad, y_quad}, SlantColumnAmountNO2Std_val]};
                SurfaceAlbedo(x_quad, y_quad) = {[SurfaceAlbedo{x_quad, y_quad}, SurfaceAlbedo_val]};
                TerrainHeight(x_quad, y_quad) = {[TerrainHeight{x_quad, y_quad}, TerrainHeight_val]};
                TroposphericVerticalColumn(x_quad, y_quad) = {[TroposphericVerticalColumn{x_quad, y_quad}, TroposphericVerticalColumn_val]};
                TroposphericVerticalColumnError(x_quad, y_quad) = {[TroposphericVerticalColumnError{x_quad, y_quad}, TroposphericVerticalColumnError_val]};
                
                % If there is no existing field
                 elseif ~isnan(AirMassFactorTropospheric_val) %JLL 19 Mar 2014: I added the logical test here, before this was just an 'else' statement, but it would make sense not to add a value if there was no valid NO2 column.
                    % Count, area, and areaweight require special handling
                    Count(x_quad,y_quad)=Count(x_quad,y_quad)+1;
                    Area(x_quad,y_quad)=pixelarea;
                    Areaweight(x_quad,y_quad)=1/pixelarea;
                    
                    % Regular fields will be a running average
                    AirMassFactorTropospheric(x_quad, y_quad) = AirMassFactorTropospheric_val;
                    AveragingKernel(x_quad, y_quad) = AveragingKernel_val;
                    CloudFraction(x_quad, y_quad) = CloudFraction_val;
                    CloudPressure(x_quad, y_quad) = CloudPressure_val;
                    GhostColumn(x_quad, y_quad) = GhostColumn_val;
                    SlantColumnAmountNO2(x_quad, y_quad) = SlantColumnAmountNO2_val;
                    SlantColumnAmountNO2Std(x_quad, y_quad) = SlantColumnAmountNO2Std_val;
                    SurfaceAlbedo(x_quad, y_quad) = SurfaceAlbedo_val;
                    TerrainHeight(x_quad, y_quad) = TerrainHeight_val;
                    TroposphericVerticalColumn(x_quad, y_quad) = TroposphericVerticalColumn_val;
                    TroposphericVerticalColumnError(x_quad, y_quad) = TroposphericVerticalColumnError_val;
                    
                    % Flag fields will append the flag value to a matrix in
                    % a cell corresponding to this grid cell
                    AirMassFactorTropospheric(x_quad, y_quad) = {AirMassFactorTropospheric_val};
                    AveragingKernel(x_quad, y_quad) = {AveragingKernel_val};
                    CloudFraction(x_quad, y_quad) = {CloudFraction_val};
                    CloudPressure(x_quad, y_quad) = {CloudPressure_val};
                    GhostColumn(x_quad, y_quad) = {GhostColumn_val};
                    SlantColumnAmountNO2(x_quad, y_quad) = {SlantColumnAmountNO2_val};
                    SlantColumnAmountNO2Std(x_quad, y_quad) = {SlantColumnAmountNO2Std_val};
                    SurfaceAlbedo(x_quad, y_quad) = {SurfaceAlbedo_val};
                    TerrainHeight(x_quad, y_quad) = {TerrainHeight_val};
                    TroposphericVerticalColumn(x_quad, y_quad) = {TroposphericVerticalColumn_val};
                    TroposphericVerticalColumnError(x_quad, y_quad) = {TroposphericVerticalColumnError_val};
                 end
            end
    	end
    end
end

% Create the OMI structure for output
OMI.AirMassFactorTropospheric = AirMassFactorTropospheric;
OMI.AveragingKernel = AveragingKernel;
OMI.CloudFraction = CloudFraction;
OMI.CloudPressure = CloudPressure;
OMI.GhostColumn = GhostColumn;
OMI.SlantColumnAmountNO2 = SlantColumnAmountNO2;
OMI.SlantColumnAmountNO2Std = SlantColumnAmountNO2Std;
OMI.SurfaceAlbedo = SurfaceAlbedo;
OMI.TerrainHeight = TerrainHeight;
OMI.TroposphericVerticalColumn = TroposphericVerticalColumn;
OMI.TroposphericVerticalColumnError = TroposphericVerticalColumnError;
OMI.Count = Count;
OMI.Area = Area;
OMI.Areaweight = Areaweight;
OMI.AirMassFactorTropospheric = AirMassFactorTropospheric;
OMI.AveragingKernel = AveragingKernel;
OMI.CloudFraction = CloudFraction;
OMI.CloudPressure = CloudPressure;
OMI.GhostColumn = GhostColumn;
OMI.SlantColumnAmountNO2 = SlantColumnAmountNO2;
OMI.SlantColumnAmountNO2Std = SlantColumnAmountNO2Std;
OMI.SurfaceAlbedo = SurfaceAlbedo;
OMI.TerrainHeight = TerrainHeight;
OMI.TroposphericVerticalColumn = TroposphericVerticalColumn;
OMI.TroposphericVerticalColumnError = TroposphericVerticalColumnError;

% Replace fill values with NaNs. Of course, we can only do this for numeric
% fields, not cells or structs
fns = fieldnames(OMI);
for a=1:numel(fns)
    if isnumeric(OMI.(fns{a}))
        OMI.(fns{a})(OMI.(fns{a}) == fill_val) = NaN;
    end
end

end
