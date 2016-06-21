%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This function produces XY dependence maps of the S2 Pulse width
%
%  Inputs:
%       - s2x                  - X position of events after golden selection
%       - s2y                  - Y position of events after golden selection
%       - s2_width             - S2 pulse width after golden selection
%       - xbin_min             - minimum value of the xbins
%       - xbin_max             - maximum value of the xbins
%       - xbin_min             - minimum value of the ybins
%       - xbin_max             - maximum value of the ybins
%       - xy_bin_size          - bin size of the xbins and ybins
%       - cut                  - Optional cut to clean up the pulse areas (defaults to all ones if not provided)
%
%  Outputs:
%       - s2_xbins             - 1xN vector defining x binning in XY map
%       - s2_ybins             - 1xM vector defining y binning in XY map
%       - means_s2_xy_width    - MxN matrix of pulse width means
%       - sigma_s2_xy_width    - MxN matrix of error on the pulse width means
% 
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s2_xbins, s2_ybins, means_s2_xy_width, sigma_s2_xy_width ] = MeasureS2PulseWidthXYDep(s2x, s2y, s2_width, xbin_min, xbin_max, ybin_min, ybin_max, xy_bin_size, cut)

    %Assign default cut if it is not provided by user
    if (nargin < 9)
        cut = logical(zeros(1,length(s2_width)) + 1).';
    end

    %Set up bins
    s2_xbins = (xbin_max-xbin_min)./xy_bin_size;
    s2_ybins = (ybin_max-ybin_min)./xy_bin_size;
    
    %Preallocate memory
    means_s2_xy_width = zeros(floor(s2_ybins),floor(s2_xbins));
    sigma_s2_xy_width = zeros(floor(s2_ybins),floor(s2_xbins));

    %Measure the XY map
    for x_bin=xbin_min:xy_bin_size:(xbin_max-xy_bin_size);
        for y_bin=ybin_min:xy_bin_size:(ybin_max-xy_bin_size);          
            x_min   = x_bin; 
            x_max   = x_bin+xy_bin_size;
            y_min   = y_bin;
            y_max   = y_bin+xy_bin_size;      

            %Matrix indices
            x_count = int32(1+x_bin/xy_bin_size+xbin_max/xy_bin_size);
            y_count = int32(1+y_bin/xy_bin_size+ybin_max/xy_bin_size); 
            
            %Select events in the XY bin
            bin_cut = inrange(s2x,[x_min,x_max]) & inrange(s2y,[y_min,y_max]) & cut; 

            %Fill in the map.  Shape is non Gaussian, so just use mean. 
            %X and Y indices swapped intentionally due to how interp2D works
            means_s2_xy_width(y_count,x_count) = mean(s2_width(bin_cut));
            sigma_s2_xy_width(y_count,x_count) = std(s2_width(bin_cut));
        end
    end

    %set = 0 if NaN is returned for mean.
    means_s2_xy_width(isnan(means_s2_xy_width))=0; 
    sigma_s2_xy_width(isnan(sigma_s2_xy_width))=0;
     
end

