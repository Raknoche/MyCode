%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This function produces two dimensional R^2 v Z maps of the S1a, S1b, and S1a/S1b  means
%
%  Inputs:
%       - s1ab_radius          - S1a/S1b radius after golden selection
%       - s1ab_z               - S1a/S1b drift time after golden selection
%       - s1a_phe_both         - S1a pulse areas after golden selection
%       - s1b_phe_both         - S1b pulse areas after golden selection
%
%  Outputs:
%       - s1ab_r2bins          - 1xN vector defining r^2 binning of the R^2 v Z map
%       - s1ab_zbins           - 1xN vector defining z binning of the R^2 v Z map
%       - s1a_r2z_mean         - MxN matrix of the S1a means
%       - s1a_r2z_mean_err     - MxN matrix of the S1a mean errors
%       - s1b_r2z_mean         - MxN matrix of the S1b means
%       - s1b_r2z_mean_err     - MxN matrix of the S1a mean errors
%       - s1ab_r2z_mean        - MxN matrix of the S1a/S1b means
%       - s1ab_r2z_mean_err    - MxN matrix of the S1a/S1b mean errors
% 
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s1ab_r2bins, s1ab_zbins, s1a_r2z_mean, s1a_r2z_mean_err, s1b_r2z_mean, s1b_r2z_mean_err, s1ab_r2z_mean, s1ab_r2z_mean_err] = MeasureS1aS1bRvZDep(s1ab_radius, s1ab_z, s1a_phe_both, s1b_phe_both)

%Define bin edges
s1ab_r2bin_min    = 0;
s1ab_r2bin_max    = 25*25;
s1ab_zbin_min     = 10;
s1ab_zbin_max     = 330;

%Define bin sizes (Aiming for ~300 events per bin)
s1ab_r2z_binnum   = (length(s1ab_z))/300;
s1ab_r2_binnum    = floor(sqrt(s1ab_r2z_binnum))-1;
s1ab_z_binnum     = floor(sqrt(s1ab_r2z_binnum));
s1ab_r2_binsize   = (s1ab_r2bin_max-s1ab_r2bin_min)/s1ab_r2_binnum;
s1ab_z_binsize    = (s1ab_zbin_max-s1ab_zbin_min)/s1ab_z_binnum;

%Define bins
s1ab_r2bins       = s1ab_r2bin_min+s1ab_r2_binsize/2:s1ab_r2_binsize:s1ab_r2bin_max;
s1ab_zbins        = s1ab_zbin_min+s1ab_z_binsize/2:s1ab_z_binsize:s1ab_zbin_max; 

%Preallocate memory
s1a_r2z_mean     = zeros(length(s1ab_zbins),length(s1ab_r2bins));
s1b_r2z_mean     = zeros(length(s1ab_zbins),length(s1ab_r2bins));
s1a_r2z_mean_err = zeros(length(s1ab_zbins),length(s1ab_r2bins));
s1b_r2z_mean_err = zeros(length(s1ab_zbins),length(s1ab_r2bins));

%Calculate the R^2vZ map
r2_counter=1;
 for r2_bin_max=s1ab_r2bin_min+s1ab_r2_binsize:s1ab_r2_binsize:s1ab_r2bin_max;
     z_counter=1;
        for z_bin_max=s1ab_zbin_min+s1ab_z_binsize:s1ab_z_binsize:s1ab_zbin_max;          
           
            %Cut to select events within the bin
            bin_cut = s1a_phe_both>0 & s1b_phe_both>0 & inrange(s1ab_radius.^2,[r2_bin_max-s1ab_r2_binsize,r2_bin_max]).' & inrange(s1ab_z,[z_bin_max-s1ab_z_binsize,z_bin_max]).'; 
            
            %If there are enough events, calculate the S1a and S1b means
            if length(s1a_phe_both(bin_cut))>100;
                clear s1a_confint  s1b_confint
                s1a_hist_max                           = 400;
                s1a_r2z_hist_bins                      = (0:2:s1a_hist_max);
                s1a_r2z_hist_cut                       = bin_cut & inrange(s1a_phe_both,[0 s1a_hist_max]); 
                s1a_r2z_fit                            = fit(s1a_r2z_hist_bins.',hist(s1a_phe_both(s1a_r2z_hist_cut),s1a_r2z_hist_bins).','gauss1');
                s1a_r2z_mean(z_counter,r2_counter)     = s1a_r2z_fit.b1;
                s1a_confint                            = confint(s1a_r2z_fit,0.68);
                s1a_r2z_mean_err(z_counter,r2_counter) = abs(s1a_confint(1,2)-s1a_r2z_fit.b1);
                
                s1b_hist_max                           = 300;
                s1b_r2z_hist_bins                      = (0:2:s1b_hist_max);
                s1b_r2z_hist_cut                       = bin_cut & inrange(s1b_phe_both,[0 s1b_hist_max]);
                s1b_r2z_fit                            = fit(s1b_r2z_hist_bins.',hist(s1b_phe_both(s1b_r2z_hist_cut),s1b_r2z_hist_bins).','gauss1');
                s1b_r2z_mean(z_counter,r2_counter)     = s1b_r2z_fit.b1;
                s1b_confint                            = confint(s1b_r2z_fit,0.68);
                s1b_r2z_mean_err(z_counter,r2_counter) = abs(s1b_confint(1,2)-s1b_r2z_fit.b1);
            else
                s1a_r2z_mean(z_counter,r2_counter)     = 0;
                s1a_r2z_mean_err(z_counter,r2_counter) = 0;
                s1b_r2z_mean(z_counter,r2_counter)     = 0;
                s1b_r2z_mean_err(z_counter,r2_counter) = 0;               
            end
            
           z_counter=z_counter+1;    
        end
        
        r2_counter=r2_counter+1;
 end
     s1ab_r2z_mean     = s1a_r2z_mean./s1b_r2z_mean;
     s1ab_r2z_mean_err = sqrt( (s1b_r2z_mean_err.*s1a_r2z_mean./(s1b_r2z_mean.^2)).^2 + (s1a_r2z_mean_err./s1b_r2z_mean).^2);
     
end
