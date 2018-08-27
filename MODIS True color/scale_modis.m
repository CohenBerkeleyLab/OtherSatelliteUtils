function [ varargout ] = scale_modis( varargin )
%SCALE_MODIS Scales MODIS truecolor data to display correctly as an image
%   See http://www.idlcoyote.com/ip_tips/brightmodis.html
%   and http://glasslab.engr.ccny.cuny.edu/u/fhs/gclust/doc/_modules/glasslab_cluster/io/modis.html#rgb

E = JLLErrors;

input_vec = [0, 30, 60, 120, 190, 255];
output_vec = [0, 110, 160, 210, 240, 255];

% Allow the user to input a single rgb 3D matrix or separate r, g, and b
% matrices.
if nargin == 1;
    rgb = varargin{1};
elseif nargin == 3;
    rgb = cat(3,varargin{:});
else
    E.callError('wrong_number_args','This function takes either 1 m x n x 3 RGB matrix or 3 m x n matrices for R, G, and B components');
end


% Remove any fill values.
rgb(rgb==-28672) = NaN;

% Scale the values to [0, 1.1] then divide by the greatest component (r, g,
% or b) per pixel to remove the offset due to broadband absorbance and
% convert reflectance into color values
rgb = double(rgb);
rgb = rgb./10000;
rgb(rgb>1.1) = 1.1; rgb(rgb<0)=0;

for a = 1:3
    rgb(:,:,a) = rgb(:,:,a) ./ max(max(rgb(:,:,a)));
end
rgb = round(rgb .* 255);

% Then rescale the new matrix using the two vectors defined above.
s = size(rgb);
rgb = rgb(:);
% rgb_scaled = zeros(size(rgb));

rgb_scaled = interp1(input_vec,output_vec,rgb);

% for a = 1:numel(input_vec)-1
%     xx = rgb >= input_vec(a) & rgb < input_vec(a+1);
%     conv = rgb(xx);
%     
%     lb_in = input_vec(a); ub_in = input_vec(a+1)-lb_in;
%     lb_out = output_vec(a); ub_out = output_vec(a+1)-lb_out;
%     
%     conv = (conv - lb_in) ./ ub_in;
%     conv = (conv .* ub_out) + lb_out;
%     
%     rgb_scaled(xx) = conv;
% end

rgb_scaled = reshape(rgb_scaled,s)./255;

if nargin == 1
    varargout{1} = rgb_scaled;
elseif nargin == 3
    varargout{1} = rgb_scaled(:,:,1);
    varargout{2} = rgb_scaled(:,:,2);
    varargout{3} = rgb_scaled(:,:,3);
end

end

