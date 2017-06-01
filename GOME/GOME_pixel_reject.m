function [ omi, RejectFlags ] = GOME_pixel_reject( omi_in, cloud_type, cloud_frac )
%[ OMI, REJECTFLAGS ] = GOME_pixel_reject( OMI_IN, CLOUD_TYPE, CLOUD_FRAC )
% Set areaweight to 0 for any pixels that will adversely
% affect the accuracy of an NO2 map.
%   There are a number of criteria that need to be evaluated for an OMI
%   pixel before it can be reliably used as an NO2 measurement.  This
%   function will set the areaweight value to 0 for any pixel which fails
%   these criteria.
%
%   Inputs:
%       omi_in: An OMI structure (gridded data) which can be obtained from
%       hdf_quadrangle_GOME.
%       cloud_type: 'geo' or 'rad'. Sets whether to use geometric or
%       radiative cloud fraction.
%       cloud_frac: The maximum allowed cloud fraction in a pixel.  If
%          using geometric cloud product ~0.1-0.2 is recommended.  For
%          radiance cloud fraction, 0.3-0.5 is a good range.
%
%   The rejection criteria are:
%       Column amount < 0: Likely a fill value or unphysical result.
%       VCD Quality Flag is not a even number: the VCD quality flag is a
%          bit array, with the least significant bit as a summary.  If this
%          bit is 1, then some error occured during the NASA retrieval and
%          we should ignore this pixel.
%       Cloud fraction too great: Cloudy pixels will have much of the
%          tropospheric column obscured, and should not contribute to the
%          average.
%       Column amount > 1E17: This magnitude of tropospheric column is
%          known to be affected by the row anomaly.
%       Column amount is NaN: NaN indicates some failure, either averaging
%          0 values, or another mathematical mistake.



omi = omi_in;

RejectFlags = uint8(zeros(size(omi.Latitude)));

xx = omi.NO2ColumnAmountTrop<=0;
omi.Areaweight(xx) = 0; %Do not average in negative tropospheric column densities
RejectFlags(xx) = bitset(RejectFlags(xx),1*ones(size(RejectFlags(xx))));

if iscell(omi.TropFlag) % The flags may be a cell array or not, depending on whether this is for Data or OMI (gridded) structure
    for a=1:numel(omi.TropFlag);
        if any([omi.TropFlag{a}]==-1)
            omi.Areaweight(a) = 0; %If any of the troposphere flags is -1, the tropospheric retrieval is not meaningful
            RejectFlags(a) = bitset(RejectFlags(a),2);
        end
    end
else
    xx = omi.TropFlag==-1;
    omi.Areaweight(xx) = 0;
    RejectFlags(xx) = bitset(RejectFlags(xx),2*ones(size(RejectFlags(xx))));
end

if strcmpi(cloud_type,'geo'); 
    % Reject pixels with NaN for cloud fraction or fill values for cloud fraction
    xx = omi.CloudFraction > cloud_frac | isnan(omi.CloudFraction) | omi.CloudFraction < 0;
elseif strcmpi(cloud_type,'rad'); 
    % Reject pixels with NaN for cloud fraction or with fill values for cloud fraction
    xx = omi.CloudRadianceFraction > cloud_frac | isnan(omi.CloudRadianceFraction) | omi.CloudRadianceFraction < 0;
end %Do not include the element if the cloud fraction is greater than the allowable criteria
omi.Areaweight(xx) = 0;
RejectFlags(xx) = bitset(RejectFlags(xx),3*ones(size(RejectFlags(xx))));

xx=isnan(omi.NO2ColumnAmountTrop); omi.NO2ColumnAmountTrop(xx)=0; omi.Areaweight(xx)=0; %Set any column NaNs to 0 and do not include them in the average
RejectFlags(xx) = bitset(RejectFlags(xx),4*ones(size(RejectFlags(xx))));

end

