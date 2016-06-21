%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This function produces Z dependence maps of the S2 pulse width
%
%  Inputs:
%       - s2_width                 - S2 Width RQ after golden selection
%       - drift_time               - Drift time after golden selection
%       - det_edge                 - Edge of the drift time histogram
%       - cut                      - cut to clean up events before fitting. Defaults to all events passing
%
%  Outputs:
%       - mean_index_s2_z_width    - 1xN matrix of z-slice bin centers
%       - means_s2_z_width         - 1xN matrix of S2_width means 
%       - means_error_s2_z_width   - 1xN matrix of S2_width mean errors
%       - P_s2_width               - Polynomial fit to S2 Width z dependence
%       - sigma_P_s2_width         - Errors of Polynomial fit to S2 Width z dependence
%       - dT_fit                   - Drift time bins used in the polynomial fit
% 
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mean_index_s2_z_width, means_s2_z_width, means_error_s2_z_width, P_s2_width, sigma_P_s2_width, dT_fit] = ...
    MeasureS2PulseWidthZDep(s2_width, drift_time, det_edge, cut)
  
    %Assign default cut if it is not provided by user
    if (nargin < 4)
        cut = logical(zeros(1,length(s2_width)) + 1).';
    end
    

    %Sort events so they match expected input
    A                                = sort(size(s2_width(cut))); 
    
    %Set bin size
    min_evts_per_bin                 = 300;
    max_num_bins                     = 65;
    bins                             = ceil(A(2)./min_evts_per_bin); 
    if bins>max_num_bins
       bins = max_num_bins; 
    end
    bin_size                         = (floor(det_edge/bins*4)/4);
   
    %Preallocate Memory
    means_s2_z_width                 = zeros(bins,1);
    means_error_s2_z_width           = zeros(bins,1);
    mean_index_s2_z_width            = zeros(bins,1);
        
    %Calculate Z Dependence map
    for bin=1:bins
        binStart                     = ((bin-1).*bin_size);  binEnd = (bin.*bin_size);   
        time_cut                     = inrange(drift_time,[binStart,binEnd]) & cut==1;

        %distribution vs. depth in non Gaussian, so just use the mean
        means_s2_z_width(bin)        = mean(s2_width(time_cut));
        means_error_s2_z_width(bin)  = std(s2_width(time_cut))/sqrt(length(s2_width(time_cut))-1); % sigma/sqrt(N-1)
        mean_index_s2_z_width(bin)   = (binStart+binEnd)/2;         
    end
    
    %Fit polynomial to Z dependence
    dT_range                         = inrange( mean_index_s2_z_width,[0, det_edge] );
    dT_fit                           = 0:1:det_edge;    
    [P_s2_width S]                   = polyfit(mean_index_s2_z_width(dT_range),means_s2_z_width(dT_range),3);
    sigma_P_s2_width                 = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df)'; 
    
end

