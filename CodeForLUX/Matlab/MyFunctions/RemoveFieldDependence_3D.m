%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This function produces field removal (3D dependence) matrices for Kr83m S1 and S2 data.
%       Note this will ONLY work for Kr83m, and not for other sources
%
%  Inputs:
%       - s1ab_xyz_xbins       - 1xN matrix for X binning for 3D matrix of S1a/S1b means
%       - s1ab_xyz_ybins       - 1xN matrix for Y binning for 3D matrix of S1a/S1b means
%       - s1ab_xyz_ybins       - 1xN matrix for Z binning for 3D matrix of S1a/S1b means
%       - s1ab_xyz_mean        - 3D matrix of S1a/S1b means
%       - s1ab_P               - Polynomial fit to S1a/S1b z dependence
%       - s2x                  - X position after golden event selection
%       - s2y                  - Y position after golden event selection
%       - drift_time           - Drift time after golden event selection
%
%  Outputs:
%       - s1_fieldremoval      - 1xN matrix of S1 field removal (apply by S1./s1_field_removal)
%       - s2_fieldremoval      - 1xN matrix of S2 field removal (apply by S2./s2_field_removal) 
%       - s1as1b_at_center     - value of S1a/S1b at the center of the detector
% 
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s1_fieldremoval, s2_fieldremoval, s1as1b_at_center] = RemoveFieldDependence_3D(s1ab_xyz_xbins, s1ab_xyz_ybins, s1ab_xyz_zbins, s1ab_xyz_mean, s1ab_P, s2x, s2y, drift_time)

%Hardcoded numbers from Energy Chi2 Fit           
s1ab_s2_xyz_fit_map.p1 = -0.4988; 
s1ab_s2_xyz_fit_map.p2 = 1.484;
s1ab_s2_xyz_fit_map.p3 = 1-s1ab_s2_xyz_fit_map.p1.*2.7303.^2-2.7303.*s1ab_s2_xyz_fit_map.p2;

s1ab_s1_xyz_fit_map.p1 = 0.065;
s1ab_s1_xyz_fit_map.p2 = 0.02;
s1ab_s1_xyz_fit_map.p3 = 1-s1ab_s1_xyz_fit_map.p1.*2.7303.^2-2.7303.*s1ab_s1_xyz_fit_map.p2;
 
%Removing zeros from map with nearest neighbor method - needed for interpolation purposes
s1ab_xyz_mean_old = s1ab_xyz_mean;
temp_map          = s1ab_xyz_mean;
q                 = ~isnan(temp_map);
r                 = nearestpoint(find(~q),find(q));
s1ab_xyz_mean     = temp_map;
temp_mapQ         = temp_map(q); 
s1ab_xyz_mean(~q) = temp_mapQ(r);
s1ab_xyz_mean     = reshape(s1ab_xyz_mean,size(s1ab_xyz_mean_old,1),size(s1ab_xyz_mean_old,2),size(s1ab_xyz_mean_old,3));

%interpolating with cubic method so that we get NaN when extrapolating
s1as1b_values     = interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1ab_xyz_mean,s2x,s2y,drift_time,'cubic');
s1as1b_at_center  = interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1ab_xyz_mean,x_center,y_center,z_center,'cubic');
s2_fieldremoval   = (s1ab_s2_xyz_fit_map.p1.*s1as1b_values.^2+s1ab_s2_xyz_fit_map.p2.*s1as1b_values+s1ab_s2_xyz_fit_map.p3); 
s1_fieldremoval   = (s1ab_s1_xyz_fit_map.p1.*s1as1b_values.^2+s1ab_s1_xyz_fit_map.p2.*s1as1b_values+s1ab_s1_xyz_fit_map.p3); 

%Filling in NaN values with the s1a/s1b Z map polynomial 
s1as1b_zvalues   = polyval(s1ab_P,drift_time);
s2_zfieldremoval = (s1ab_s2_xyz_fit_map.p1.*s1as1b_zvalues.^2+s1ab_s2_xyz_fit_map.p2.*s1as1b_zvalues+s1ab_s2_xyz_fit_map.p3);
s1_zfieldremoval = (s1ab_s1_xyz_fit_map.p1.*s1as1b_zvalues.^2+s1ab_s1_xyz_fit_map.p2.*s1as1b_zvalues+s1ab_s1_xyz_fit_map.p3); 
s2_fieldremoval(isnan(s2_fieldremoval)) = s2_zfieldremoval(isnan(s2_fieldremoval));
s1_fieldremoval(isnan(s1_fieldremoval)) = s1_zfieldremoval(isnan(s1_fieldremoval));

end

