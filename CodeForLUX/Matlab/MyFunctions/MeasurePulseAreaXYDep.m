%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This function produces XY dependence maps of the given pulse area
%
%  Inputs:
%       - bin_min              - Minimum range of the Gaussian fit histograms
%       - bin_max              - Maximum range of the Gaussian fit histograms
%       - bin_size             - Bin size of the Gaussian fit histograms
%       - xbin_min             - minimum value of the xbins 
%       - xbin_max             - maximum value of the xbins
%       - xbin_min             - minimum value of the ybins
%       - xbin_max             - maximum value of the ybins
%       - xy_bin_size          - bin size of the xbins and ybins
%       - pulse_area           - Pulse area after golden event selection
%       - s2x                  - X positions after golden event selection
%       - s2y                  - Y position after golden event selection
%       - cut                  - Optional cut to clean up the pulse areas (defaults to all ones if not provided)
%
%  Outputs:
%       - mean_both_xy         - MxN matrix of pulse area means
%       - Sigma_both           - MxN matrix of one-sigma error on the means (NOT gaussian width)
%       - Count_both           - Number of events in each bin
% 
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mean_both_xy, Sigma_both, Count_both] = ...
    MeasurePulseAreaXYDep(bin_min, bin_max, bin_size, xbin_min, xbin_max, ybin_min, ybin_max, xy_bin_size, pulse_area, s2x, s2y, cut)

    %Assign default cut if it is not provided by user
    if (nargin < 12)
        cut = logical(zeros(1,length(pulse_area)) + 1).';
    end
    
    %set up energy histgoram bins
    x2           = bin_min:bin_size:bin_max;
    x2           = x2';

    %Preallocating Memory
    mean_both_xy = zeros(floor((ybin_max-ybin_min)./xy_bin_size),floor((xbin_max-xbin_min)./xy_bin_size));
    Sigma_both   = zeros(floor((ybin_max-ybin_min)./xy_bin_size),floor((xbin_max-xbin_min)./xy_bin_size));
    Count_both   = zeros(floor((ybin_max-ybin_min)./xy_bin_size),floor((xbin_max-xbin_min)./xy_bin_size));

    %Make the XY Map
    for x_bin=xbin_min:xy_bin_size:(xbin_max-xy_bin_size);
        for y_bin=ybin_min:xy_bin_size:(ybin_max-xy_bin_size);          
            x_min = x_bin; x_max = x_bin+xy_bin_size;
            y_min = y_bin; y_max = y_bin+xy_bin_size;      

            %Indicies for the matrix
            x_count=int32(1+x_bin/xy_bin_size+xbin_max/xy_bin_size);
            y_count=int32(1+y_bin/xy_bin_size+ybin_max/xy_bin_size); 

            %Select events in the XY bin
            bin_cut = inrange(s2x,[x_min,x_max]) & inrange(s2y,[y_min,y_max]) & cut;            

            hist_bin = hist(pulse_area(bin_cut) ,x2)'/bin_size; 
            Count_both(y_count,x_count)=sum(hist_bin)*bin_size; % X and Y flip ... because MATLAB

            %If there are enough events, make the map entry. Otherwise fill the map with zero
            %X and Y indicies swapped intentionally, due to how interp2D works
            if Count_both(y_count,x_count)>350;                         
                
                %Fit a Gaussian to the data
                Fit_both_xy= fit(x2(2:end-1),hist_bin(2:end-1),'gauss1'); 
                mean_both_xy(y_count,x_count)=Fit_both_xy.b1; 
                Sig_b=Fit_both_xy.c1/sqrt(2)/sqrt(Count_both(y_count,x_count));
                
                %Fit error checking
                if(isnan(Sig_b))
                    mean_both_xy(y_count,x_count)=0;
                    Sig_b=0;
                end

                %uncertainty in mean
                Sigma_both(y_count,x_count)=Sig_b;

            end             
        end     
    end

    
    mean_both_xy(isnan(mean_both_xy)) = 0;

     
end

