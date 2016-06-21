%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This function produces field removal (Z dependence only) matrices for Kr83m S1 and S2 data.
%       Note this will ONLY work for Kr83m, and not for other sources
%
%  Inputs:
%       - s1ab_P               - Polynomial fit to S1a/S1b z dependence
%       - drift_time           - Drift time after golden event selection
%
%  Outputs:
%       - s1_fieldremoval      - 1xN matrix of S1 field removal (apply by S1./s1_field_removal)
%       - s2_fieldremoval      - 1xN matrix of S2 field removal (apply by S2./s2_field_removal) 
% 
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s1_fieldremoval, s2_fieldremoval] = RemoveFieldDependence_1D(s1ab_P, drift_time)

%Hardcoded numbers from Energy Chi2 Fit           
s1ab_s2_xyz_fit_map.p1 = -0.4988; 
s1ab_s2_xyz_fit_map.p2 = 1.484;
s1ab_s2_xyz_fit_map.p3 = 1-s1ab_s2_xyz_fit_map.p1.*2.7303.^2-2.7303.*s1ab_s2_xyz_fit_map.p2;

s1ab_s1_xyz_fit_map.p1 = 0.065;
s1ab_s1_xyz_fit_map.p2 = 0.02;
s1ab_s1_xyz_fit_map.p3 = 1-s1ab_s1_xyz_fit_map.p1.*2.7303.^2-2.7303.*s1ab_s1_xyz_fit_map.p2;
 
%Getting field removal using Z dependence polynomials
s1as1b_zvalues        = polyval(s1ab_P,drift_time);
s2_fieldremoval       = (s1ab_s2_xyz_fit_map.p1.*s1as1b_zvalues.^2+s1ab_s2_xyz_fit_map.p2.*s1as1b_zvalues+s1ab_s2_xyz_fit_map.p3);
s1_fieldremoval       = (s1ab_s1_xyz_fit_map.p1.*s1as1b_zvalues.^2+s1ab_s1_xyz_fit_map.p2.*s1as1b_zvalues+s1ab_s1_xyz_fit_map.p3); 

%Filling in NaN values with (Z dependence should never be NaN anyway)
s2_fieldremoval(isnan(s2_fieldremoval))=1;
s1_fieldremoval(isnan(s1_fieldremoval))=1;


end

