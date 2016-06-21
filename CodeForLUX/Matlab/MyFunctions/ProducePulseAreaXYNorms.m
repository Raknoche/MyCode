%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This function produces XY normalization maps from XY pulse area mean maps
%
%  Inputs:
%       - xbins              - 1xN vector defining x binning of XY means map
%       - ybins              - 1xM vector defining y binning of XY means map
%       - xymap              - MxN map of pulse area means
%       - xymap_sigma        - MxN map of 1-sigma error on pulse area means
%       - x_center           - X coordinate to normalize the XY map to 
%       - y_center           - Y coordinate to normalize the XY map to 
%
%  Outputs:
%       - center             - Value of the XY mean map at (x_center,y_center)
%       - sigma_center       - One sigma error on the value at the XY mean map at (x_center, y_center)
%       - norm               - MxN matrix of normalization factors
%       - sigma_norm         - MxN matrix of errors on the normalization factors
% 
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [center, sigma_center, norm, sigma_norm] = ProducePulseAreaXYNorms(xbins, ybins, xymap, xymap_sigma, x_center, y_center)

center                          = interp2(xbins,ybins,xymap,x_center,y_center,'cubic');%Normalize to the center
norm                            = center./xymap; 
norm(isinf(norm))               = 1;
norm(isnan(norm))               = 1;

sigma_center                    = interp2(xbins,ybins,xymap_sigma,x_center,y_center,'cubic');%Normalize to the center
sigma_norm                      = sqrt((sigma_center./xymap).^2+(xymap_sigma.*center./xymap.^2).^2);
sigma_norm(isinf(sigma_norm))   = 1;
sigma_norm(isnan(sigma_norm))   = 1;

     
end

