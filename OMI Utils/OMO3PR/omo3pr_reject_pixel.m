function [ rejects ] = omo3pr_reject_pixel( Data )
%OMO3PR_REJECT_PIXEL Returns logical array TRUE for pixels to ignore
%   REJECTS = OMO3PR_REJECT_PIXEL( DATA ) Given the structure DATA, return
%   REJECTS, a logical array the same size as DATA.Longitude, that is TRUE
%   for any pixel that should not be used.

cldfraccrit = 0.2;
aotcrit = 1e-5;

% Reject any pixels with the cloud fraction too great in either UV channel
rejects = false(size(Data.Longitude));
rejects(Data.EffectiveCloudFractionUV1 > cldfraccrit | Data.EffectiveCloudFractionUV2 > cldfraccrit) = true;

% Elevated aerosol layers cause too much O3 to be put in the troposphere,
% so reject any pixel with basically detectable aerosol
rejects(Data.AerosolOpticalThickness > aotcrit) = true;

% Also reject any pixel affected by the row anomaly (see
% https://disc.gsfc.nasa.gov/Aura/data-holdings/OMI/documents/v003/OMO3PRO_README.shtml)
rejects(Data.ReflectanceCostFunction > 30) = true;

% No mention of filtering on the MeasurementQualityFlags was made in the
% readme

end

