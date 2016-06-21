%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This function produces fiducial volume cuts for Kr S1a and S1b pulses 
%
%  Inputs:
%       - s1ab_z             - S1a/S1b z coordinates after golden cut is applied
%       - s2radius           - S1a/S1b radius coordinates after golden cut is applied
%
%  Outputs:
%       - s1ab_z_cut           - Fidiucal cut for making XY dependence maps with S1a/S1b
%       - s1ab_r_cut           - Fidiucal cut for making Z dependence maps with S1a/S1b
%       - s1ab_z_cut_upper     - Upper limit on S1a/S1b drift time histogram 
%
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ s1ab_z_cut, s1ab_r_cut, s1ab_z_cut_upper ] = GetS1aS1bFiducialCuts(s1ab_z, s1ab_radius)

%Additional Z cut for s1ab_xy measurement
s1ab_z_hist_binning                = (0:1:400); %Drift time range
s1ab_z_hist_cut                    = inrange(s1ab_z,[0,400]);
[s1ab_z_hist s1ab_z_hist_bins]     = hist(s1ab_z(s1ab_z_hist_cut),s1ab_z_hist_binning);
frac_hist                          = cumsum(s1ab_z_hist)./sum(s1ab_z_hist);
s1ab_z_cut_lower                   = min(s1ab_z_hist_bins(frac_hist>.10));
s1ab_z_cut_upper                   = max(s1ab_z_hist_bins(frac_hist<.90));
s1ab_z_cut                         = inrange(s1ab_z,[s1ab_z_cut_lower,s1ab_z_cut_upper]);

%Additional radius cut for s1ab_Z measurement
s1ab_r_hist_binning                = (0:1:26);
s1ab_r_upperhist_cut               = (s1ab_z>s1ab_z_cut_upper & s1ab_radius<=26);
s1ab_r_lowerhist_cut               = (s1ab_z<s1ab_z_cut_lower & s1ab_radius<=26);
[upperz_r_hist upperz_r_hist_bins] = hist(s1ab_radius(s1ab_r_upperhist_cut),s1ab_r_hist_binning);
[lowerz_r_hist lowerz_r_hist_bins] = hist(s1ab_radius(s1ab_r_lowerhist_cut),s1ab_r_hist_binning);
frac_lowerz_r_hist                 = cumsum(lowerz_r_hist)./sum(lowerz_r_hist); %Find 80% radial edge of bottom 10% drift time
lowerz_rbound                      = max(lowerz_r_hist_bins(frac_lowerz_r_hist<0.80));
frac_upperz_r_hist                 = cumsum(upperz_r_hist)./sum(upperz_r_hist); %Find 90% radial edge of top 90% drift time
upperz_rbound                      = max(upperz_r_hist_bins(frac_upperz_r_hist<0.90));   
rslope                             = (s1ab_z_cut_lower-s1ab_z_cut_upper)/(lowerz_rbound-upperz_rbound); %find slope and intercept of radial cut line
rintercept                         = s1ab_z_cut_lower-rslope.*lowerz_rbound;
s1ab_r_cut                         = (s1ab_z-rslope.*s1ab_radius)<rintercept;    

end

