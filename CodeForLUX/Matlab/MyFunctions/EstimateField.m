%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%   This function takes s2radius and s2drift inputs in the form of 1D vectors.  It
%   uses Scott Hertel's field simulation to determine the field at the location of each data point.
%   The output of this function is a 1D vector of the same length of the original data which 
%   has corresponding field values for each data point in the units of V/cm.
% 
%   We use a cubic interpolation of Scott's map within the bounds of his simulated data
%   and use a nearest neighbor interpolation outside of the bounds.  The bounds are defined by 
%   fitting a spline to a rvZ scatterplot in the simulated data, with the additional constraint
%   that any simulated data with drift_time<20 uSec is thrown out since the simualtion seems 
%   inaccurate in that region.  Additionally, any real data that has radius < 1 cm is 
%   considered to be outside of the simulation bounds, since a cubic interpolation fails in that region
%   despite the data being within the simulation bounds.
%
%  Inputs:
%       - s2radius             - 1xN Event radii after golden selection
%       - s2drift              - 1xN Event drift times after golden selection
%       - year                 - A string defining which of Scott's simulations to use.  ('2014' or '2015)
%
%  Outputs:
%       - FieldValues          - 1xN Estimated field values in the vicinity of each event
% 
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ FieldValues ] = EstimateField( s2radius, s2drift, year )

%% Load simualted field data from Scott

switch(year)
    case '2014'
        load('C:\Program Files\MATLAB\R2012a\bin\Run04_ERBand\Kr2p22\sep2014_field_map.mat');
    case '2015'
        load('C:\Program Files\MATLAB\R2012a\bin\Run04_ERBand\Kr2p22\sep2015_field_map.mat')
end
        

%% Clean up simulated field data

%reshape the data
sim_radius  = reshape(sep2014_field_map.r_end_cm,1,size(sep2014_field_map.r_end_cm,1)*size(sep2014_field_map.r_end_cm,2));
sim_drift   = reshape(sep2014_field_map.t_end_us,1,size(sep2014_field_map.t_end_us,1)*size(sep2014_field_map.t_end_us,2));
sim_field   = reshape(sep2014_field_map.E_start_Vcm,1,size(sep2014_field_map.E_start_Vcm,1)*size(sep2014_field_map.E_start_Vcm,2));

% Remove NaN from simulated field data
nan_cut     = (~isnan(sim_radius) & ~isnan(sim_drift) & ~isnan(sim_field));
sim_radius  = sim_radius(nan_cut);
sim_drift   = sim_drift(nan_cut);
sim_field   = sim_field(nan_cut);

%% Make a 2D map of the field, with radius_bin and drift_bin vectors

% Can use griddata to interpolate the Field Map, but it doesn't work well for points outside the boundary
% Instead, we first find the boundary, use nearest neighbor interp outside of it, and cubic interp inside of it

%Finds the boundary of the radius, drift_time scatter plot
clear r_center max_z
i       = 1;
r_step  = 0.5;

for r_max=[0.45:r_step:25];
    r_center(i) = mean(sim_radius(inrange(sim_radius,[r_max-r_step,r_max]))); %#ok<AGROW>
    max_z(i)    = max(sim_drift(inrange(sim_radius,[r_max-r_step,r_max]))); %#ok<AGROW>
    i=i+1;
end
 
%Construct the boundary from a spline fit to the maximum values found above
boundary_spline                  = csapi(r_center,max_z);

%Evaluates if data is inside or outside the boundary
data_boundary                    = fnval(boundary_spline, s2radius); 

%A cut that is one in inside, zero if outside.  Has extra specifications to avoid bad simulation data at low drift and radius
data_inside_bound                = (s2drift<=data_boundary & s2drift>=20 & s2radius>=1); 

%To get field value for the data, first split into inside/outside bound data vectors
data_inbound_radius              = s2radius(data_inside_bound);
data_inbound_drift               = s2drift(data_inside_bound);
data_outbound_radius             = s2radius(~data_inside_bound);
data_outbound_drift              = s2drift(~data_inside_bound);

%Interpolate: outside bound use nearest neighbor, inside bound use cubic
field_inbound_value              = griddata(sim_drift,sim_radius,sim_field,data_inbound_drift,data_inbound_radius,'cubic');

%only uses nearest for sim_drift>=20 uSec to avoid problem at low drift in sim
field_outbound_value             = griddata(sim_drift(sim_drift>=20),sim_radius(sim_drift>=20),sim_field(sim_drift>=20),data_outbound_drift,data_outbound_radius,'nearest'); 

%Patch inbound and outbound field values together
FieldValues                      = zeros(length(field_inbound_value)+length(field_outbound_value),1);
FieldValues(data_inside_bound)   = field_inbound_value;
FieldValues(~data_inside_bound)  = field_outbound_value;

end
