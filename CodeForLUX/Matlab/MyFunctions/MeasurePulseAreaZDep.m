%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This function measures the Z dependence of a pulse area RQ
%
%  Inputs:
%       - bin_min              - Minimum of the histograms used during Gaussian fitting
%       - bin_max              - Maximum of the histograms used during Gaussian fitting
%       - bin_size             - Bin width of the histograms used during Gaussian fitting
%       - pulse_area           - Pulse areas to be fit
%       - drift_time           - Drift time of the pulse areas
%       - cut                  - Optional cut to clean up the data (defaults to every event passing if not specified)
%
%  Outputs:
%       - means                - 1xM vector of mean pulse areas in each drift time
%       - means_error          - 1xM vector of error on the mean pulse areas in each drift time
%       - mean_index           - 1xM vector of containing center of each drift time bin
%       - Fit_both             - 1xM cell containing Gaussian fits to each of the drift time bins
%       - bins                 - total number of drift time bins
%       - x                    - pulse area binning of the histogram used for each Gaussian fit
%       - hist_both            - histogram data for each of the drift time bins.  
%                                  (First index is the count in the pulse area bins, second index is the corresponding drift time bin)
% 
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [means,means_error,mean_index,Fit_both, bins, x, hist_both] = MeasurePulseAreaZDep(bin_min, bin_max, bin_size, pulse_area, drift_time, det_edge, cut)
    
    %Assign default cut if it is not provided by user
    if (nargin < 7)
        cut = logical(zeros(1,length(pulse_area)) + 1).';
    end

    %set up energy bins for the gaussian fit
    x            = bin_min:bin_size:bin_max;
    x            = x';
    
    %sort, incase 1xn or nx1 vector is defined.
    A            = sort(size(pulse_area(cut))); 
    
    %get the number of points, adjust binning if stats are low.
    evts_per_bin = 300;
    max_num_bins = 65;
    bins         = ceil(A(2)./evts_per_bin);
    if bins>max_num_bins
       bins      = max_num_bins;
    end
    bin_size     = (floor(det_edge/bins*4)/4);%span up to the det_edge

    %Preallocate memory
    hist_both    = zeros(length(x),bins);
    Fit_both     = cell(1,bins);
    means        = zeros(bins,1);
    means_error  = zeros(bins,1);
    mean_index   = zeros(bins,1);
    
    %Calculate the z-dependence of the data
    for bin=1:bins
        
        %Select data in the bin
        binStart           = ((bin-1).*bin_size);  
        binEnd             = (bin.*bin_size);   
        time_cut           = inrange(drift_time,[binStart,binEnd]) & cut==1;
       
        %Do a Gaussian fit to the bin
        hist_both(:,bin)   = hist(pulse_area(time_cut),x)'/bin_size;
        temp_hist          = hist_both(:,bin);   
        Fit_both{bin}      = fit(x(2:end-1),temp_hist(2:end-1),'gauss1');%remove bin edges

        %Fill in the map
        means(bin)         = Fit_both{bin}.b1;
        means_error(bin)   = Fit_both{bin}.c1/sqrt(2)/sqrt(sum(hist_both(:,bin))*bin_size); % == sigma/sqrt(N);
        mean_index(bin)    = (binStart+binEnd)/2;     
    end

end

