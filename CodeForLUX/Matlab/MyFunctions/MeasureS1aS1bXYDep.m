%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This function produces XY dependence maps of the S1a, S1b, and S1a/S1b means
%
%  Inputs:
%       - s1ab_x               - S1a/S1b X coordinates after golden selection
%       - s1ab_y               - S1a/S1b Y coordinates after golden selection
%       - s1ab_z               - S1a/S1b drift time after golden selection
%       - s1a_phe_both_z       - Z corrected S1a pulse area after golden selection
%       - s1b_phe_both_z       - Z corrected S1b pulse area after golden selection
%       - s1ab_z_cut           - S1a/S1b drift time cut from GetS1aS1bFiducialCuts
%
%  Outputs:
%       - s1ab_xbins           - 1xN vector defining x binning in XY map
%       - s1ab_ybins           - 1xM vector defining y binning in XY map
%       - s1a_xy_mean          - MxN matrix of S1a means
%       - s1a_xy_mean_err      - MxN matrix of S1a mean errors
%       - s1b_xy_mean          - MxN matrix of S1b means
%       - s1b_xy_mean_err      - MxN matrix of S1b mean errors
%       - s1ab_xy_mean         - MxN matrix of S1a/S1ab means
%       - s1ab_xy_mean_err     - MxN matrix of S1a/S1b mean errors
% 
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s1ab_xbins, s1ab_ybins, s1a_xy_mean, s1a_xy_mean_err, s1b_xy_mean, s1b_xy_mean_err, s1ab_xy_mean, s1ab_xy_mean_err] = MeasureS1aS1bXYDep(s1ab_x, s1ab_y, s1ab_z, s1a_phe_both_z, s1b_phe_both_z, s1ab_z_cut)

%Define bounds of map
s1ab_xbin_min  = -25;
s1ab_ybin_min  = -25;
s1ab_xbin_max  = 25;
s1ab_ybin_max  = 25;
s1ab_xybinsize = sqrt((50*50)/(length(s1ab_z)/300));
s1a_hist_bins  = (0:2:400);
s1b_hist_bins  = (0:2:300);

%Set up XY bins for map
s1ab_xbins     = s1ab_xbin_min+s1ab_xybinsize/2:s1ab_xybinsize:s1ab_xbin_max;
s1ab_ybins     = s1ab_ybin_min+s1ab_xybinsize/2:s1ab_xybinsize:s1ab_ybin_max; 

%Preallocate memory
s1a_xy_mean     = zeros(length(s1ab_xbins),length(s1ab_ybins));
s1b_xy_mean     = zeros(length(s1ab_xbins),length(s1ab_ybins));
s1a_xy_mean_err = zeros(length(s1ab_xbins),length(s1ab_ybins));
s1b_xy_mean_err = zeros(length(s1ab_xbins),length(s1ab_ybins));

%Measure S1a/S1b XY Dependence
for x_bin=s1ab_xbin_min:s1ab_xybinsize:(s1ab_xbin_max-s1ab_xybinsize/2);
	for y_bin=s1ab_ybin_min:s1ab_xybinsize:(s1ab_ybin_max-s1ab_xybinsize/2);          
        x_min = x_bin; x_max = x_bin+s1ab_xybinsize;
        y_min = y_bin; y_max = y_bin+s1ab_xybinsize;      
        
        %Indicies for the map
        x_count=int32(1+x_bin/s1ab_xybinsize+s1ab_xbin_max/s1ab_xybinsize);
        y_count=int32(1+y_bin/s1ab_xybinsize+s1ab_ybin_max/s1ab_xybinsize); 
        
        %Cut to select events in the XY bin
        bin_cut = s1ab_z_cut & ...
                  inrange(s1ab_z,[10 330]) & ...
                  s1a_phe_both_z>0 & ...
                  s1b_phe_both_z>0 & ...
                  inrange(s1ab_x,[x_min,x_max]) & ...
                  inrange(s1ab_y,[y_min,y_max]); 
        
        %If there are enough events, calculate the S1a and S1b mean for the bin     
        if length(s1a_phe_both_z(bin_cut))>100;
            s1a_hist_cut                     = bin_cut & inrange(s1a_phe_both_z,[0 400]);
            s1a_z_fit                        = fit(s1a_hist_bins.',hist(s1a_phe_both_z(s1a_hist_cut),s1a_hist_bins).','gauss1');
            s1a_xy_mean(y_count,x_count)     = s1a_z_fit.b1;
            s1a_xy_mean_err(y_count,x_count) = s1a_z_fit.c1/sqrt(2)/length(s1a_phe_both_z(s1a_hist_cut));

            s1b_hist_cut                     = bin_cut & inrange(s1b_phe_both_z,[0 300]);           
            s1b_z_fit                        = fit(s1b_hist_bins.',hist(s1b_phe_both_z(s1b_hist_cut),s1b_hist_bins).','gauss1');
            s1b_xy_mean(y_count,x_count)     = s1b_z_fit.b1;
            s1b_xy_mean_err(y_count,x_count) = s1b_z_fit.c1/sqrt(2)/length(s1b_phe_both_z(s1b_hist_cut));
        else
            s1a_xy_mean(y_count,x_count)     = 0;
            s1a_xy_mean_err(y_count,x_count) = 0;
            s1b_xy_mean(y_count,x_count)     = 0;
            s1b_xy_mean_err(y_count,x_count) = 0;   
        end
                
	end
end
 
     %Calculate S1a/S1b maps from S1a and S1b maps
     s1ab_xy_mean     = s1a_xy_mean./s1b_xy_mean;
     s1ab_xy_mean_err = sqrt( (s1b_xy_mean_err.*s1a_xy_mean./(s1b_xy_mean.^2)).^2 + (s1a_xy_mean_err./s1b_xy_mean).^2);
     
end

