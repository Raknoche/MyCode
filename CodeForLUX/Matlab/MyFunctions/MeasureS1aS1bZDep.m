%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This function produces Z dependence maps and polynomial fits for the S1a, S1b, and S1a/S1b means
%
%  Inputs:
%       - det_edge             - Edge of the drift time histogram (typically 330 uSec in Run04)
%       - s1ab_z               - S1a/S1b drift time after golden selection
%       - s1a_phe_both         - S1a pulse area after golden selection
%       - s1b_phe_both         - S1b pulse area after golden selection
%       - s1ab_r_cut           - S1a/S1b radial cut from GetS1aS1bFiducialCuts
%       - s1ab_z_cut_upper     - Upper limit on S1a/S1b drift time histogram from GetS1aS1bFiducialCuts
%
%  Outputs:
%       - s1ab_bincenters      - 1xN matrix of z-slice bin centers
%       - s1a_z_means          - 1xN matrix of S1a means 
%       - s1a_z_means_err      - 1xN matrix of S1a mean errors
%       - s1b_z_means          - 1xN matrix of S1b means 
%       - s1b_z_means_err      - 1xN matrix of S1b mean errors
%       - s1ab_z_means         - 1xN matrix of S1a/S1b means 
%       - s1ab_z_means_err     - 1xN matrix of S1a/S1b mean errors
%       - s1a_P                - Polynomial fit to S1a z dependence
%       - s1a_S                - Errors of Polynomial fit to S1a z dependence
%       - s1b_P                - Polynomial fit to S1b z dependence
%       - s1b_S                - Errors of Polynomial fit to S1b z dependence
%       - s1ab_P               - Polynomial fit to S1a/S1b z dependence
%       - s1ab_S               - Errors of Polynomial fit to S1a/S1b z dependence
% 
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s1ab_bincenters, s1a_z_means, s1a_z_means_err, s1b_z_means, s1b_z_means_err, s1ab_z_means, s1ab_z_means_err, s1a_P, s1a_S, s1b_P, s1b_S, s1ab_P, s1ab_S] = MeasureS1aS1bZDep(det_edge, s1ab_z, s1a_phe_both, s1b_phe_both, s1ab_r_cut, s1ab_z_cut_upper)

%Set up bin size for Z dependence
dT_step         = det_edge/(length(s1ab_z)/300); %300 events per z bin    

%Set up variables for fitting
s1a_fit_bins    = (0:2:400);
s1b_fit_bins    = (0:2:300);

%Preallocate memory
map_size        = length(10+dT_step:dT_step:s1ab_z_cut_upper);
s1a_z_means     = zeros(1,map_size);
s1a_z_means_err = zeros(1,map_size);
s1b_z_means     = zeros(1,map_size);
s1b_z_means_err = zeros(1,map_size);
s1ab_bincenters = zeros(1,map_size);

%Measure Z Dependence of S1a/S1b
i=1;
for z_max=10+dT_step:dT_step:s1ab_z_cut_upper;
    s1a_fit_cut        = (s1ab_r_cut.' & s1a_phe_both>0 & s1a_phe_both<400 & inrange(s1ab_z,[z_max-dT_step,z_max]).');
    s1b_fit_cut        = (s1ab_r_cut.' & s1b_phe_both>0 & s1b_phe_both<300 & inrange(s1ab_z,[z_max-dT_step,z_max]).');
    s1a_fit            = fit(s1a_fit_bins.',hist(s1a_phe_both(s1a_fit_cut),s1a_fit_bins).','gauss1');
    s1b_fit            = fit(s1b_fit_bins.',hist(s1b_phe_both(s1b_fit_cut),s1b_fit_bins).','gauss1');
    s1a_z_means(i)     = s1a_fit.b1;
    s1a_z_means_err(i) = s1a_fit.c1/sqrt(2)/sqrt(length(s1a_phe_both(s1a_fit_cut)));
    s1b_z_means(i)     = s1b_fit.b1;
    s1b_z_means_err(i) = s1b_fit.c1/sqrt(2)/sqrt(length(s1b_phe_both(s1b_fit_cut)));
    s1ab_bincenters(i) = z_max-dT_step/2;
    i=i+1;
end

s1ab_z_means     = s1a_z_means./s1b_z_means;
s1ab_z_means_err = sqrt( (s1b_z_means_err.*s1a_z_means./(s1b_z_means.^2)).^2 + (s1a_z_means_err./s1b_z_means).^2);
    
%fit polynomial to s1az and s1bz, to remove z dependence later in xy fit
[s1a_P, s1a_S]   = polyfit(s1ab_bincenters,s1a_z_means,3);
[s1b_P, s1b_S]   = polyfit(s1ab_bincenters,s1b_z_means,3);
[s1ab_P, s1ab_S] = polyfit(s1ab_bincenters,s1ab_z_means,2);   

end

