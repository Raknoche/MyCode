%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ____    ____       _____    _____      _____      _____         _____          ____    ____        
%|    |  |    |  ___|\    \  |\    \    /    /| ___|\    \    ___|\    \    ____|\   \  |    |       
%|    |  |    | |    |\    \ | \    \  /    / ||    |\    \  /    /\    \  /    /\    \ |    |       
%|    | /    // |    | |    ||  \____\/    /  /|    | |    ||    |  |    ||    |  |    ||    |       
%|    |/ _ _//  |    |/____/  \ |    /    /  / |    |/____/||    |  |____||    |__|    ||    |  ____ 
%|    |\    \'  |    |\    \   \|___/    /  /  |    ||    |||    |   ____ |    .--.    ||    | |    |
%|    | \    \  |    | |    |      /    /  /   |    ||____|/|    |  |    ||    |  |    ||    | |    |
%|____|  \____\ |____| |____|     /____/  /    |____|       |\ ___\/    /||____|  |____||____|/____/|
%|    |   |    ||    | |    |    |`    | /     |    |       | |   /____/ ||    |  |    ||    |     ||
%|____|   |____||____| |____|    |_____|/      |____|        \|___|    | /|____|  |____||____|_____|/
%  \(       )/    \(     )/         )/           \(            \( |____|/   \(      )/    \(    )/   
%
%
%
%  function [ lug_val ] = LUXkrypCal(d,user_name,dp_version,algorithm_name, submit,inpaint_on)
%
%  This function will process Kr83m data, extract light response maps (in XYZ) and IQs, and output them to the LUG.  
%
%
%  Inputs
%       - d             - the data structure to be analysed.
%       - user_name     - string identify user on the LUG
%       - dp_version    - sting corresponding to the DP version. Ex. '2.22'
%       - submit        - 1 or 0 to submit to LUG or not. By default a 0
%       - inpaint_on    - Option to inpaint empty matrix elements in correction maps.
%                           (default is 0)
%
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [lug_val] = KrypCal(d,user_name,dp_version,algorithm_name, submit)

%% Checking the input arguments, and filling in defaults

if (nargin <= 4)
    submit     = 0;
    
elseif(nargin <= 3)
    algorithm_name = 'test_KrypCal';
    submit         = 0;
    
elseif (nargin <= 2)
    dp_version     = '9.99'; %9.99 is dp number assigned to tests
    algorithm_name = 'test_KrypCal';
    submit         = 0;

elseif (nargin <= 1)
    user_name      = 'RichardKnoche';
    dp_version     = '9.99'; %9.99 is dp number assigned to tests
    algorithm_name = 'test_KrypCal';
    submit         = 0;
    
end

%Check that the data struct isn't empty
if ~exist('d','var') || isempty(d)
    error('No valid data set provided');
end

%Calibrations was originally a user input to specify which calibrations to do
%The feature was never used, so now it is hard coded to do all calibrations
calibrations=ones(6,1); 


%% Defining constant variables to be used later in the code

%Extract relevant file name and cp information
file       = char(d.files_loaded(1));
file_id    = file(1:19);
cp_start   = strfind(file,'cp')+2;       %find location of 'cp' in string
cp         = file(cp_start:cp_start+4);  %cp number has four digits
file_id_cp = [file_id '_cp' cp];
        
%Query the LUG to get gs number for the dataset
query_str  = ['select * from lug_complete_process_record where cp="' cp '" ;' ];
data2      = MySQLQuery_UMD(query_str,'read'); 
if str2double(cp) == data2.cp 
    gs = data2.gs;
else
    error('CP number not found on the LUG');
end

%Make a folder to save the plots
mkdir(strcat('LUX_corrections/',file_id_cp)) 

%Define correction map limits
%Hard coded since these don't change
grid_size        = 3;     %cm - This is a default, and isn't used if there are enough events
rcut_min         = 0;     %cm
rcut_max         = 25;    %cm
s1area_bound_min = 50;    %phe
s1area_bound_max = 650;   %phe
s2area_bound_min = 1000;  %phe
s2area_bound_max = 70000; %phe


%Edges of the S1, S2, and SE maps
S1_xbin_min = -25;
S1_xbin_max = 25;
S1_ybin_min = -25;
S1_ybin_max = 25;

S2_xbin_min = -25;
S2_xbin_max = 25;
S2_ybin_min = -25;
S2_ybin_max = 25;

%Edges of the drift time maps
det_edge    = 330;           %Edge of the drift time histogram (found by eye)
zcut_min    = 10;            %Extraction field starts at 4 uSec
zcut_max    = 0.95*det_edge; %Do not want to calculate maps with data near the Cathode


%% Event selection with golden cuts

[s1_single_cut, s2_single_cut] = GetGoldenCuts(d, s1area_bound_min, s1area_bound_max, s2area_bound_min, s2area_bound_max);


%% Extracting pulse information

%Cleaning up the data struct
d.z_drift_samples(isnan(d.z_drift_samples)) = 0.0;     %get rid of NaN in RQ       
d.phe_bottom   = d.pulse_area_phe./(1+d.top_bottom_ratio); %bottom PMT pulse area

%Pulse area quantities
s1_phe_both    = d.pulse_area_phe(s1_single_cut);
s1_phe_bottom  = d.phe_bottom(s1_single_cut);
s2_phe_both    = d.pulse_area_phe(s2_single_cut);
s2_phe_bottom  = d.phe_bottom(s2_single_cut);

%Position Reconstruction
drift_time     = d.z_drift_samples(s2_single_cut)/100;  % Converts samples (in units of 10 nS) to uSeconds
s2x            = d.x_cm(s2_single_cut);
s2y            = d.y_cm(s2_single_cut);   
s2radius       = (s2x.^2+s2y.^2).^(0.5);

%Event timing
d.livetime_sec = sum(d.livetime_end_samples-d.livetime_latch_samples)/1e8;
    
%Pulse width information
s2_width       = (d.aft_t2_samples(s2_single_cut)-d.aft_t0_samples(s2_single_cut)); %cut above 800 samples
s1_width       = d.pulse_end_samples(s1_single_cut)-d.pulse_start_samples(s1_single_cut);
  
%Number of events
kr_energy_cut  = inrange(s1_phe_both,[s1area_bound_min s1area_bound_max]) & inrange(s2_phe_both,[s2area_bound_min s2area_bound_max]);%for counting Kr events
Kr_events      = length(drift_time(kr_energy_cut));% Are there enough events for the XY or XYZ map?

%Center of the detector in uncorrected coordinates
x_center       = mean(s2x(drift_time>4)); %exclude extraction field region at uSec<4
y_center       = mean(s2y(drift_time>4));
z_center       = mean(drift_time(drift_time>4));


%% Calculate XY binning for correction maps

%If there are enough events we use variable bin sizes
minimum_events_for_binning = 30000;
if Kr_events < minimum_events_for_binning 
    S1xybinsize = grid_size;
    S2xybinsize = grid_size;           
else
    S1xybinsize = sqrt((50*50)/(Kr_events/600));%Note: Due to the radial field, bins near the edge won't actually have 600 events in them
    S2xybinsize = sqrt((50*50)/(Kr_events/600));%cm
end

%Set up XY binning for pulse area maps
s2xbins = (S2_xbin_min+S2xybinsize/2):S2xybinsize:(S2_xbin_max-S2xybinsize/2);
s2ybins = (S2_ybin_min+S2xybinsize/2):S2xybinsize:(S2_ybin_max-S2xybinsize/2);
s1xbins = (S1_xbin_min+S1xybinsize/2):S1xybinsize:(S1_xbin_max-S1xybinsize/2);
s1ybins = (S1_ybin_min+S1xybinsize/2):S1xybinsize:(S1_ybin_max-S1xybinsize/2);


%% Calculate Raw S1 and S2 means - used to track variation over time on LUG

% Raw S1 at center
s1_fit_bins                           = (0:10:2000);
s1_fit_cut                            = inrange(drift_time,[z_center-5, z_center+5]) & inrange(s1_phe_both,[0 2000]);
s1_phe_both_raw_gauss_fit             = fit(s1_fit_bins.',hist(s1_phe_both(s1_fit_cut),s1_fit_bins).','gauss1');
s1_phe_both_raw_gauss_fit_mean        = s1_phe_both_raw_gauss_fit.b1;
s1_phe_both_raw_gauss_fit_sigma       = s1_phe_both_raw_gauss_fit.c1/sqrt(2);
s1_phe_both_raw_gauss_fit_confint     = confint(s1_phe_both_raw_gauss_fit,0.68);
s1_phe_both_raw_gauss_fit_mean_error  = abs(s1_phe_both_raw_gauss_fit_confint(1,2)-s1_phe_both_raw_gauss_fit_mean);
s1_phe_both_raw_gauss_fit_sigma_error = abs(s1_phe_both_raw_gauss_fit_confint(1,3)/sqrt(2)-s1_phe_both_raw_gauss_fit_sigma);
 
% Raw S2 at top
s2_fit_bins                           = (0:100:70000);
s2_fit_cut                            = inrange(drift_time,[5, 15]) & inrange(s2_phe_both,[0 70000]);
s2_phe_both_raw_gauss_fit             = fit(s2_fit_bins.',hist(s2_phe_both(s2_fit_cut),s2_fit_bins).','gauss1');
s2_phe_both_raw_gauss_fit_mean        = s2_phe_both_raw_gauss_fit.b1;
s2_phe_both_raw_gauss_fit_sigma       = s2_phe_both_raw_gauss_fit.c1/sqrt(2);
s2_phe_both_raw_gauss_fit_confint     = confint(s2_phe_both_raw_gauss_fit,0.68);
s2_phe_both_raw_gauss_fit_mean_error  = abs(s2_phe_both_raw_gauss_fit_confint(1,2)-s2_phe_both_raw_gauss_fit_mean);
s2_phe_both_raw_gauss_fit_sigma_error = abs(s2_phe_both_raw_gauss_fit_confint(1,3)/sqrt(2)-s2_phe_both_raw_gauss_fit_sigma); 


%% Calculate S1a/S1b Z map

%Set up cuts
[ s1ab_cut ] = GetS1aS1bCuts(d, s1_single_cut, s2_single_cut);

%Get S1a and S1b quantities
s1ab_x                             = d.x_cm(s2_single_cut & repmat(s1ab_cut,10,1));
s1ab_y                             = d.y_cm(s2_single_cut & repmat(s1ab_cut,10,1));
s1ab_z                             = d.z_drift_samples(s2_single_cut & repmat(s1ab_cut,10,1))/100;
s1ab_radius                        = sqrt(s1ab_x.^2+s1ab_y.^2);
s1a_phe_both                       = d.Kr83fit_s1a_area_phe(s1ab_cut);
s1b_phe_both                       = d.Kr83fit_s1b_area_phe(s1ab_cut);

%Get additional fiducial cuts for S1a/S1b maps
[ s1ab_z_cut, s1ab_r_cut, s1ab_z_cut_upper ]         = GetS1aS1bFiducialCuts(s1ab_z, s1ab_radius);

%Measure S1a/S1b Z dependence
[s1ab_bincenters, s1a_z_means, s1a_z_means_err, s1b_z_means, s1b_z_means_err, s1ab_z_means, ...
    s1ab_z_means_err, s1a_P, s1a_S, s1b_P, s1b_S, s1ab_P, s1ab_S] = MeasureS1aS1bZDep(det_edge, s1ab_z, s1a_phe_both, s1b_phe_both, s1ab_r_cut, s1ab_z_cut_upper);

%Make and save S1a/S1b Z dependence plot    
s1az_plot_bins = (10:1:320);

s1az_plot = figure;
hold on;
errorbar(s1ab_bincenters,s1a_z_means,s1a_z_means_err,'.k')   
xlabel('Drift Time (uSec)');
ylabel('S1a Mean (phe)'); 
title('S1a Z Dependence'); 
myfigview(16);
plot(s1az_plot_bins,polyval(s1a_P,s1az_plot_bins),'-r','LineWidth',2)
saveas(s1az_plot,strcat('LUX_corrections/',file_id_cp,'/s1a_mean_z',file_id_cp),'fig'); 
saveas(s1az_plot,strcat('LUX_corrections/',file_id_cp,'/s1a_mean_z',file_id_cp),'jpg'); 
  
s1bz_plot = figure;
hold on;
errorbar(s1ab_bincenters,s1b_z_means,s1b_z_means_err,'.k')   
xlabel('Drift Time (uSec)');
ylabel('S1b Mean (phe)'); 
title('S1b Z Dependence'); 
myfigview(16);
plot(s1az_plot_bins,polyval(s1b_P,s1az_plot_bins),'-r','LineWidth',2)
saveas(s1bz_plot,strcat('LUX_corrections/',file_id_cp,'/s1b_mean_z',file_id_cp),'fig'); 
saveas(s1bz_plot,strcat('LUX_corrections/',file_id_cp,'/s1b_mean_z',file_id_cp),'jpg');    
    
s1a_over_bz_plot=figure;
hold on;
errorbar(s1ab_bincenters,s1ab_z_means,s1ab_z_means_err,'.k')   
xlabel('Drift Time (uSec)');
ylabel('S1a/b'); 
title('S1 a/b Z Dependence'); 
myfigview(16);
plot(s1az_plot_bins,polyval(s1ab_P,s1az_plot_bins),'-r','LineWidth',2)
saveas(s1a_over_bz_plot,strcat('LUX_corrections/',file_id_cp,'/s1ab_ratio_z',file_id_cp),'fig'); 
saveas(s1a_over_bz_plot,strcat('LUX_corrections/',file_id_cp,'/s1ab_ratio_z',file_id_cp),'jpg');    


%% Calculate S1a/S1b XY map (after z correction)    

%Correcting z dependence to get XY dependence
s1a_phe_both_z=s1a_phe_both.'.*polyval(s1a_P,z_center)./polyval(s1a_P,s1ab_z);
s1b_phe_both_z=s1b_phe_both.'.*polyval(s1b_P,z_center)./polyval(s1b_P,s1ab_z);

%Measure S1a/S1b XY map
[s1ab_xbins, s1ab_ybins, s1a_xy_mean, s1a_xy_mean_err, s1b_xy_mean, s1b_xy_mean_err, s1ab_xy_mean, s1ab_xy_mean_err] = MeasureS1aS1bXYDep(s1ab_x, s1ab_y, s1ab_z, s1a_phe_both_z, s1b_phe_both_z, s1ab_z_cut);
   
%Make and Save Plots of S1a XY means
s1a_xy_fig = figure;
color_range_max=max(max(s1a_xy_mean));
color_range_min=min(min(s1a_xy_mean(s1a_xy_mean>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s1ab_xbins,s1ab_ybins,s1a_xy_mean,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. S1a Mean vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'S1a Mean (phe)','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);
saveas(s1a_xy_fig,strcat('LUX_corrections/',file_id_cp,'/s1a_xy',file_id_cp),'fig');
saveas(s1a_xy_fig,strcat('LUX_corrections/',file_id_cp,'/s1a_xy',file_id_cp),'jpg');
    
%Plot s1a XY means
s1b_xy_fig = figure;
color_range_max=max(max(s1b_xy_mean));
color_range_min=min(min(s1b_xy_mean(s1b_xy_mean>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s1ab_xbins,s1ab_ybins,s1b_xy_mean,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. S1b Mean vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'S1b Mean (phe)','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);
saveas(s1b_xy_fig,strcat('LUX_corrections/',file_id_cp,'/s1b_xy',file_id_cp),'fig');
saveas(s1b_xy_fig,strcat('LUX_corrections/',file_id_cp,'/s1b_xy',file_id_cp),'jpg');
    
%s1a/b XY Ratio
s1ab_xy_fig = figure;
color_range_max=max(max(s1ab_xy_mean));
color_range_min=min(min(s1ab_xy_mean(s1ab_xy_mean>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s1ab_xbins,s1ab_ybins,s1ab_xy_mean,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. S1a/b Ratio vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'S1a/b Ratio','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);
saveas(s1ab_xy_fig,strcat('LUX_corrections/',file_id_cp,'/s1ab_ratio_xy',file_id_cp),'fig'); 
saveas(s1ab_xy_fig,strcat('LUX_corrections/',file_id_cp,'/s1ab_ratio_xy',file_id_cp),'jpg');    
    

%% Calculate S1a/S1b 3D map


if length(s1ab_z)>100000;
    [s1ab_xyz_xbins, s1ab_xyz_ybins, s1ab_xyz_zbins, s1a_xyz_mean, s1a_xyz_mean_err, ...
        s1b_xyz_mean, s1b_xyz_mean_err, s1ab_xyz_mean, s1ab_xyz_mean_err] = ...
        MeasureS1aS1b3DDep(s1ab_x, s1ab_y, s1ab_z, det_edge, s1a_phe_both, s1b_phe_both);
end


%% Calculate S1a/S1b 3D map R^2 v Z map 

[s1ab_r2bins, s1ab_zbins, s1a_r2z_mean, s1a_r2z_mean_err, s1b_r2z_mean, s1b_r2z_mean_err, ...
    s1ab_r2z_mean, s1ab_r2z_mean_err] = MeasureS1aS1bRvZDep(s1ab_radius, s1ab_z, s1a_phe_both, s1b_phe_both);

%%Make and Save Plots of S1a R^2 v Z maps
s1a_r2z_fig = figure;
color_range_max=max(max(s1a_r2z_mean));
color_range_min=min(min(s1a_r2z_mean(s1a_r2z_mean>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s1ab_r2bins,s1ab_zbins,s1a_r2z_mean,vc,'LineColor','none');
xlabel('R^2 (cm^2)','fontsize', 18);
ylabel('Z (uSec)','fontsize', 18);
title(strcat(file_id_cp, '. S1a Mean vs. R^2 v Z.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'S1a Mean (phe)','FontSize',16)
set(g,'FontSize',18)
hold on
axis([0 625 0 330])
myfigview(16);
set(gca,'Ydir','reverse')
saveas(s1a_r2z_fig,strcat('LUX_corrections/',file_id_cp,'/s1a_r2z_fig',file_id_cp),'fig'); 
saveas(s1a_r2z_fig,strcat('LUX_corrections/',file_id_cp,'/s1a_r2z_fig',file_id_cp),'jpg');    
   
%Plot s1a R2Z means   
s1b_r2z_fig = figure;
color_range_max=max(max(s1b_r2z_mean));
color_range_min=min(min(s1b_r2z_mean(s1b_r2z_mean>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s1ab_r2bins,s1ab_zbins,s1b_r2z_mean,vc,'LineColor','none');
xlabel('R^2 (cm^2)','fontsize', 18);
ylabel('Z (uSec)','fontsize', 18);
title(strcat(file_id_cp, '. S1b Mean vs. R^2 v Z.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'S1b Mean (phe)','FontSize',16)
set(g,'FontSize',18)
hold on
axis([0 625 0 330])
myfigview(16);
set(gca,'Ydir','reverse')
saveas(s1b_r2z_fig,strcat('LUX_corrections/',file_id_cp,'/s1b_r2z_fig',file_id_cp),'fig'); 
saveas(s1b_r2z_fig,strcat('LUX_corrections/',file_id_cp,'/s1b_r2z_fig',file_id_cp),'jpg');    
    
%s1a/b R2Z Ratio
s1ab_r2z_fig = figure;
color_range_max=max(max(s1ab_r2z_mean));
color_range_min=min(min(s1ab_r2z_mean(s1ab_r2z_mean>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s1ab_r2bins,s1ab_zbins,s1ab_r2z_mean,vc,'LineColor','none');
xlabel('R^2 (cm^2)','fontsize', 18);
ylabel('Z (uSec)','fontsize', 18);
title(strcat(file_id_cp, '. S1a/S1b Mean vs. R^2 v Z.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'S1a/b Ratio','FontSize',16)
set(g,'FontSize',18)
hold on
axis([0 625 0 330])
myfigview(16);
set(gca,'Ydir','reverse')
saveas(s1ab_r2z_fig,strcat('LUX_corrections/',file_id_cp,'/s1ab_r2z_fig',file_id_cp),'fig'); 
saveas(s1ab_r2z_fig,strcat('LUX_corrections/',file_id_cp,'/s1ab_r2z_fig',file_id_cp),'jpg');    
 

%% Convert S1a/S1b R^2 v Z to field measurement based on Scott's Chi by Eye

load('2015_Field_to_S1aS1b.mat');
[field_map_2015, field_map_2015_staterr]=polyval(S1aS1b_fit_2015,s1ab_r2z_mean,S1aS1b_fit_2015_error);

load('2014_Field_to_S1aS1b.mat');
[field_map_2014, field_map_2014_staterr]=polyval(S1aS1b_fit_2014,s1ab_r2z_mean,S1aS1b_fit_2014_error);

field_map_syserr=((field_map_2015-(field_map_2014+field_map_2015)/2).^2 + (field_map_2014-(field_map_2014+field_map_2015)/2).^2).^(1/2);
mean_field_map=(field_map_2014+field_map_2015)/2;
mean_field_map_totalerr=sqrt(field_map_syserr.^2+((field_map_2014_staterr+field_map_2015_staterr)/2).^2);


%Mean R^2vZ Field Map
field_r2z_fig = figure;
color_range_max=max(max(mean_field_map));
color_range_min=min(min(mean_field_map(mean_field_map>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s1ab_r2bins,s1ab_zbins,mean_field_map,vc,'LineColor','none');
xlabel('R^2 (cm^2)','fontsize', 18);
ylabel('Z (uSec)','fontsize', 18);
title(strcat(file_id_cp, '. Electric Field vs. R^2 v Z.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'Electric Field (V/cm)','FontSize',16)
set(g,'FontSize',18)
hold on
axis([0 625 0 330])
myfigview(16);
set(gca,'Ydir','reverse')
saveas(field_r2z_fig,strcat('LUX_corrections/',file_id_cp,'/field_r2z_fig',file_id_cp),'fig'); 
saveas(field_r2z_fig,strcat('LUX_corrections/',file_id_cp,'/field_r2z_fig',file_id_cp),'jpg'); 

%Mean R^2vZ Field Map Error
field_error_r2z_fig = figure;
color_range_max=max(max(mean_field_map_totalerr));
color_range_min=min(min(mean_field_map_totalerr(mean_field_map_totalerr>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s1ab_r2bins,s1ab_zbins,mean_field_map_totalerr,vc,'LineColor','none');
xlabel('R^2 (cm^2)','fontsize', 18);
ylabel('Z (uSec)','fontsize', 18);
title(strcat(file_id_cp, '. Electric Field Error vs. R^2 v Z.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'Electric Field Error (V/cm)','FontSize',16)
set(g,'FontSize',18)
hold on
axis([0 625 0 330])
myfigview(16);
set(gca,'Ydir','reverse')
saveas(field_error_r2z_fig,strcat('LUX_corrections/',file_id_cp,'/field_error_r2z_fig',file_id_cp),'fig'); 
saveas(field_error_r2z_fig,strcat('LUX_corrections/',file_id_cp,'/field_error_r2z_fig',file_id_cp),'jpg'); 

%Mean R^2vZ Field Map Fractional Error
field_frac_error_r2z_fig = figure;
color_range_max=max(max(mean_field_map_totalerr./mean_field_map));
color_range_min=min(min(mean_field_map_totalerr(mean_field_map_totalerr./mean_field_map>0)./mean_field_map(mean_field_map_totalerr./mean_field_map>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s1ab_r2bins,s1ab_zbins,mean_field_map_totalerr./mean_field_map,vc,'LineColor','none');
xlabel('R^2 (cm^2)','fontsize', 18);
ylabel('Z (uSec)','fontsize', 18);
title(strcat(file_id_cp, '. Electric Field Fractional Error vs. R^2 v Z.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'Electric Field Fractional Error','FontSize',16)
set(g,'FontSize',18)
hold on
axis([0 625 0 330])
myfigview(16);
set(gca,'Ydir','reverse')
saveas(field_frac_error_r2z_fig,strcat('LUX_corrections/',file_id_cp,'/field_frac_error_r2z_fig',file_id_cp),'fig'); 
saveas(field_frac_error_r2z_fig,strcat('LUX_corrections/',file_id_cp,'/field_frac_error_r2z_fig',file_id_cp),'jpg'); 

%2015 R^2vZ Field Map
field_r2z_2015_fig = figure;
color_range_max=max(max(field_map_2015));
color_range_min=min(min(field_map_2015(field_map_2015>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s1ab_r2bins,s1ab_zbins,field_map_2015,vc,'LineColor','none');
xlabel('R^2 (cm^2)','fontsize', 18);
ylabel('Z (uSec)','fontsize', 18);
title(strcat(file_id_cp, '. Electric Field (Based on 2015) vs. R^2 v Z.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'2015 Based Field Value (V/cm)','FontSize',16)
set(g,'FontSize',18)
hold on
axis([0 625 0 330])
myfigview(16);
set(gca,'Ydir','reverse')
saveas(field_r2z_2015_fig,strcat('LUX_corrections/',file_id_cp,'/field_r2z_2015_fig',file_id_cp),'fig'); 
saveas(field_r2z_2015_fig,strcat('LUX_corrections/',file_id_cp,'/field_r2z_2015_fig',file_id_cp),'jpg'); 

%2014 R^2vZ Field Map
field_r2z_2014_fig = figure;
color_range_max=max(max(field_map_2014));
color_range_min=min(min(field_map_2014(field_map_2014>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s1ab_r2bins,s1ab_zbins,mean_field_map,vc,'LineColor','none');
xlabel('R^2 (cm^2)','fontsize', 18);
ylabel('Z (uSec)','fontsize', 18);
title(strcat(file_id_cp, '. Electric Field (Based on 2014) vs. R^2 v Z.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'2014 Based Field Value (V/cm)','FontSize',16)
set(g,'FontSize',18)
hold on
axis([0 625 0 330])
myfigview(16);
set(gca,'Ydir','reverse')
saveas(field_r2z_2014_fig,strcat('LUX_corrections/',file_id_cp,'/field_r2z_2014_fig',file_id_cp),'fig'); 
saveas(field_r2z_2014_fig,strcat('LUX_corrections/',file_id_cp,'/field_r2z_2014_fig',file_id_cp),'jpg'); 


%% Removing field dependence from S1 and S2 data

%Using the value of s1a/s1b for each event, apply the field map to remove recombination fluctuations from the drift field
%Note: s1as1b field map is normalized such that field normalization equals one at the center of the detector AT THE TIME THE MAP WAS MADE!  

%Remove 3D dependence of field if there are enough events
%Fall back to removing Z dependence only if there aren't enough events
if length(s1ab_z)>100000;
    [s1_fieldremoval, s2_fieldremoval, s1as1b_at_center] = RemoveFieldDependence_3D(s1ab_xyz_xbins, s1ab_xyz_ybins, s1ab_xyz_zbins, s1ab_xyz_mean, s1ab_P, s2x, s2y, drift_time);
    
else 
    [s1_fieldremoval, s2_fieldremoval] = RemoveFieldDependence_1D(s1ab_P, drift_time);

end 

s2_phe_bottom=s2_phe_bottom./(s2_fieldremoval);
s2_phe_both=s2_phe_both./(s2_fieldremoval);
s1_phe_bottom=s1_phe_bottom./(s1_fieldremoval);
s1_phe_both=s1_phe_both./(s1_fieldremoval); 


%% Select SE events

%Produce cut for single electron selection
[ SE_good_cut ]      = GetSECuts( d, s1_single_cut, s2_single_cut);

%Define pulse quantities
SE_phe_both          = d.pulse_area_phe(SE_good_cut);
SE_phe_bottom        = d.phe_bottom(SE_good_cut);
SE_phe_top           = SE_phe_both-SE_phe_bottom;
SE_x                 = d.x_cm(SE_good_cut);
SE_y                 = d.y_cm(SE_good_cut);
SE_radius            = (SE_x.^2+SE_y.^2).^(1/2);
SE_max_radius        = 17; %#ok<*NASGU> %Used for the maps in the next blocks of code
num_SE_events        = length(SE_phe_both);


%% Calculate average single electron size from Kr events, using both PMT arrays

%Fit SE Size using both PMT arrays
SE_both_hist_binning = (0:0.2:60);
[ SE_mu_both, SE_sig_mu_both, SE_sig_both, SE_sig_sig_both, skew_both, SE_cut, xfit, yfit_both ] = MeasureSESize( SE_both_hist_binning, SE_phe_both, SE_radius);

%Make plot of the result
both_SE_fig=figure;
hold on;
step(SE_cut,SE_both_hist_binning,'k');  
h_fit=plot(xfit,yfit_both,'-r','linewidth',2);
box on;
myfigview(18);
xlabel('Pulse Area Phe (Phe) '); ylabel('Count/Phe ');
legend([h_fit],strcat('\mu = ',num2str(SE_mu_both,'%2.2f'),'\newline \sigma= ',num2str(SE_sig_both,'%2.2f'),'\newline skew=', num2str(skew_both,'%2.2f') ), 'location','northeast'); %#ok<NBRAK>
saveas(both_SE_fig,strcat('LUX_corrections/',file_id_cp,'/both_SE_',file_id_cp),'jpg');
saveas(both_SE_fig,strcat('LUX_corrections/',file_id_cp,'/both_SE_',file_id_cp),'fig'); 
   

%% Calculate XY map of single electron size from Kr events, using both PMT arrays

[se_xbins, se_ybins, mean_SE_both_xy, skew_SE_both_xy, sigma_SE_both_xy, ...
    kurt_SE_both_xy, mean_err_SE_both_xy,sigma_err_SE_both_xy] = MeasureSEXYDep(SE_phe_both, SE_x, SE_y, SE_both_hist_binning);


%Calculate corrections matrix and 1 sigma corrections matrix for SE_both
SE_both_mu_center                                   = interp2(se_xbins,se_ybins,mean_SE_both_xy,x_center,y_center,'cubic');%Normalize to the center 
norm_SE_both_mu                                     = SE_both_mu_center./mean_SE_both_xy; 
norm_SE_both_mu(isinf(norm_SE_both_mu))             = 1;
norm_SE_both_mu(isnan(norm_SE_both_mu))             = 1;

%Calculate error on SE center
error_SE_both_mu_center                             = interp2(se_xbins,se_ybins,mean_err_SE_both_xy,x_center,y_center,'cubic');%Normalize to the center 
error_SE_both_mu_norm                               = sqrt((error_SE_both_mu_center./mean_SE_both_xy).^2+(mean_err_SE_both_xy.*SE_both_mu_center./mean_SE_both_xy.^2).^2);
error_SE_both_mu_norm(isinf(error_SE_both_mu_norm)) = 1;
error_SE_both_mu_norm(isnan(error_SE_both_mu_norm)) = 1;

%Calculate kurt,sigma, and skew at center
SE_both_kurt_center                                 = interp2(se_xbins,se_ybins,kurt_SE_both_xy,x_center,y_center,'cubic');
SE_both_skew_center                                 = interp2(se_xbins,se_ybins,skew_SE_both_xy,x_center,y_center,'cubic');
SE_both_sigma_center                                = interp2(se_xbins,se_ybins,sigma_SE_both_xy,x_center,y_center,'cubic');

%Plot the mean SE_XY both
SE_xy_both_fig = figure;
color_range_max=max(max(mean_SE_both_xy));
color_range_min=min(min(mean_SE_both_xy(mean_SE_both_xy>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(se_xbins,se_ybins,mean_SE_both_xy,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. SE Mean (Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);
saveas(SE_xy_both_fig,strcat('LUX_corrections/',file_id_cp,'/SE_XY_both_',file_id_cp),'jpg');
saveas(SE_xy_both_fig,strcat('LUX_corrections/',file_id_cp,'/SE_XY_both_',file_id_cp),'fig');

%Plot the 1 sigma of the mean
sigma_SE_xy_both_fig = figure;
color_range_max=max(max(mean_err_SE_both_xy));
color_range_min=min(min(mean_err_SE_both_xy(mean_err_SE_both_xy>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(se_xbins,se_ybins,mean_err_SE_both_xy,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. 1-Sigma SE Mean(Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);
saveas(sigma_SE_xy_both_fig,strcat('LUX_corrections/',file_id_cp,'/Sigma_SE_XY_both_',file_id_cp),'jpg');
saveas(sigma_SE_xy_both_fig,strcat('LUX_corrections/',file_id_cp,'/Sigma_SE_XY_both_',file_id_cp),'fig');

    
%% Calculate average single electrons size from Kr events, using bottom PMT array

%Fit SE Size using bottom PMT arrays
SE_bot_hist_binning = (0:0.2:30);
[ SE_mu_bot, SE_sig_mu_bot, SE_sig_bot, SE_sig_sig_bot, skew_bot, SE_cut_bot, xfit, yfit_bot ] = MeasureSESize( SE_bot_hist_binning, SE_phe_bottom, SE_radius);
    
bottom_SE_fig=figure;
hold on;
step(SE_cut_bot,SE_bot_hist_binning,'k');  
h_fit=plot(xfit,yfit_bot,'-r','linewidth',2);
box on;
myfigview(18);
xlabel('Pulse Area Phe (Phe) '); ylabel('Count/Phe ');
legend([h_fit],strcat('\mu = ',num2str(SE_mu_bot,'%2.2f'),'\newline \sigma= ',num2str(SE_sig_bot,'%2.2f'),'\newline skew=', num2str(skew_bot,'%2.2f') ), 'location','northeast');  %#ok<NBRAK>
saveas(bottom_SE_fig,strcat('LUX_corrections/',file_id_cp,'/bottom_SE_',file_id_cp),'jpg');
saveas(bottom_SE_fig,strcat('LUX_corrections/',file_id_cp,'/bottom_SE_',file_id_cp),'fig');   
    
    
%% Calculate XY map of single electron size from Kr events, using bottom PMT arrays

[se_xbins, se_ybins, mean_SE_bot_xy, skew_SE_bot_xy, sigma_SE_bot_xy, ...
    kurt_SE_bot_xy, mean_err_SE_bot_xy,sigma_err_SE_bot_xy] = MeasureSEXYDep(SE_phe_bottom, SE_x, SE_y, SE_bot_hist_binning);


%Calculate corrections matrix and 1 sigma corrections matrix for SE_bot
SE_bot_mu_center                                   = interp2(se_xbins,se_ybins,mean_SE_bot_xy,x_center,y_center,'cubic');%Normalize to the center 
norm_SE_bot_mu                                     = SE_bot_mu_center./mean_SE_bot_xy; 
norm_SE_bot_mu(isinf(norm_SE_bot_mu))              = 1;
norm_SE_bot_mu(isnan(norm_SE_bot_mu))              = 1;

%Calculate error on SE center
error_SE_bot_mu_center                             = interp2(se_xbins,se_ybins,mean_err_SE_bot_xy,x_center,y_center,'cubic');%Normalize to the center 
error_SE_bot_mu_norm                               = sqrt((error_SE_bot_mu_center./mean_SE_bot_xy).^2+(mean_err_SE_bot_xy.*SE_bot_mu_center./mean_SE_bot_xy.^2).^2);
error_SE_bot_mu_norm(isinf(error_SE_bot_mu_norm))  = 1;
error_SE_bot_mu_norm(isnan(error_SE_bot_mu_norm))  = 1;

%Calculate kurt,sigma, and skew at center
SE_bot_kurt_center                                 = interp2(se_xbins,se_ybins,kurt_SE_bot_xy,x_center,y_center,'cubic');
SE_bot_skew_center                                 = interp2(se_xbins,se_ybins,skew_SE_bot_xy,x_center,y_center,'cubic');
SE_bot_sigma_center                                = interp2(se_xbins,se_ybins,sigma_SE_bot_xy,x_center,y_center,'cubic');

%%Plot the mean SE_XY bot %%%%%%%%%%%%%%%%
SE_xy_bot_fig = figure;
color_range_max=max(max(mean_SE_bot_xy));
color_range_min=min(min(mean_SE_bot_xy(mean_SE_bot_xy>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(se_xbins,se_ybins,mean_SE_bot_xy,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. SE Mean (Bot PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);
saveas(SE_xy_bot_fig,strcat('LUX_corrections/',file_id_cp,'/SE_XY_bot_',file_id_cp),'jpg');
saveas(SE_xy_bot_fig,strcat('LUX_corrections/',file_id_cp,'/SE_XY_bot_',file_id_cp),'fig');

%%Plot the 1 sigma of the mean
sigma_SE_xy_bot_fig = figure;
color_range_max=max(max(mean_err_SE_bot_xy));
color_range_min=min(min(mean_err_SE_bot_xy(mean_err_SE_bot_xy>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(se_xbins,se_ybins,mean_err_SE_bot_xy,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. 1-Sigma SE Mean(Bot PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);
saveas(sigma_SE_xy_bot_fig,strcat('LUX_corrections/',file_id_cp,'/Sigma_SE_XY_bot_',file_id_cp),'jpg');
saveas(sigma_SE_xy_bot_fig,strcat('LUX_corrections/',file_id_cp,'/Sigma_SE_XY_bot_',file_id_cp),'fig');

   
%% Calculate average single electrons size from Kr events, using top PMT array

%Fit SE Size using top PMT arrays
SE_top_hist_binning = (0:0.2:30);
[ SE_mu_top, SE_sig_mu_top, SE_sig_top, SE_sig_sig_top, skew_top, SE_cut_top, xfit, yfit_top ] = MeasureSESize( SE_top_hist_binning, SE_phe_top, SE_radius);
        
top_SE_fig=figure;
hold on;
step(SE_cut_top,SE_top_hist_binning,'k');
h_fit=plot(xfit,yfit_top,'-r','linewidth',2);
box on;
myfigview(18);
xlabel('Pulse Area Phe (Phe) '); ylabel('Count/Phe ');
legend([h_fit],strcat('\mu = ',num2str(SE_mu_top,'%2.2f'),'\newline \sigma= ',num2str(SE_sig_top,'%2.2f'),'\newline skew=', num2str(skew_top,'%2.2f') ), 'location','northeast'); %#ok<NBRAK>
saveas(top_SE_fig,strcat('LUX_corrections/',file_id_cp,'/top_SE_',file_id_cp),'jpg');
saveas(top_SE_fig,strcat('LUX_corrections/',file_id_cp,'/top_SE_',file_id_cp),'fig');
    

%% Calculate XY map of single electron size from KR events, using top PMT array

[se_xbins, se_ybins, mean_SE_top_xy, skew_SE_top_xy, sigma_SE_top_xy, ...
    kurt_SE_top_xy, mean_err_SE_top_xy,sigma_err_SE_top_xy] = MeasureSEXYDep(SE_phe_top, SE_x, SE_y, SE_top_hist_binning);


%Calculate corrections matrix and 1 sigma corrections matrix for SE_top
SE_top_mu_center                                   = interp2(se_xbins,se_ybins,mean_SE_top_xy,x_center,y_center,'cubic');%Normalize to the center 
norm_SE_top_mu                                     = SE_top_mu_center./mean_SE_top_xy; 
norm_SE_top_mu(isinf(norm_SE_top_mu))              = 1;
norm_SE_top_mu(isnan(norm_SE_top_mu))              = 1;

%Calculate error on SE center
error_SE_top_mu_center                             = interp2(se_xbins,se_ybins,mean_err_SE_top_xy,x_center,y_center,'cubic');%Normalize to the center 
error_SE_top_mu_norm                               = sqrt((error_SE_top_mu_center./mean_SE_top_xy).^2+(mean_err_SE_top_xy.*SE_top_mu_center./mean_SE_top_xy.^2).^2);
error_SE_top_mu_norm(isinf(error_SE_top_mu_norm))  = 1;
error_SE_top_mu_norm(isnan(error_SE_top_mu_norm))  = 1;

%Calculate kurt,sigma, and skew at center
SE_top_kurt_center                                 = interp2(se_xbins,se_ybins,kurt_SE_top_xy,x_center,y_center,'cubic');
SE_top_skew_center                                 = interp2(se_xbins,se_ybins,skew_SE_top_xy,x_center,y_center,'cubic');
SE_top_sigma_center                                = interp2(se_xbins,se_ybins,sigma_SE_top_xy,x_center,y_center,'cubic');

%%Plot the mean SE_XY top %%%%%%%%%%%%%%%%
SE_xy_top_fig = figure;
color_range_max=max(max(mean_SE_top_xy));
color_range_min=min(min(mean_SE_top_xy(mean_SE_top_xy>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(se_xbins,se_ybins,mean_SE_top_xy,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. SE Mean (top PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);
saveas(SE_xy_top_fig,strcat('LUX_corrections/',file_id_cp,'/SE_XY_top_',file_id_cp),'jpg');
saveas(SE_xy_top_fig,strcat('LUX_corrections/',file_id_cp,'/SE_XY_top_',file_id_cp),'fig');

%%Plot the 1 sigma of the mean
sigma_SE_xy_top_fig = figure;
color_range_max=max(max(mean_err_SE_top_xy));
color_range_min=min(min(mean_err_SE_top_xy(mean_err_SE_top_xy>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(se_xbins,se_ybins,mean_err_SE_top_xy,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. 1-Sigma SE Mean(top PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);
saveas(sigma_SE_xy_top_fig,strcat('LUX_corrections/',file_id_cp,'/Sigma_SE_XY_top_',file_id_cp),'jpg');
saveas(sigma_SE_xy_top_fig,strcat('LUX_corrections/',file_id_cp,'/Sigma_SE_XY_top_',file_id_cp),'fig');


%% Make Density plot of R^2 v Z distribution of events
 
%densityplot is computationally expensive, cut sample size down to 40k
    
s2radius_d=s2radius(s2radius<25);
drift_time_d=drift_time(s2radius<25);

if length(s2radius_d) > 40000;
    index_perm= randperm(length(s2radius_d)); % randomly redistrubute the event numbers
    s2radius_d=s2radius_d(index_perm(1:40000),:); %cut out 60k events
    drift_time_d=drift_time_d(index_perm(1:40000),:); %cut out 60k events
end

MAP=colormap;
MAP(1,1:3)=1; % turn 0 from dark blue to white.

radial_field_fig=figure;
densityplot(s2radius_d.^2,drift_time_d,'rectangle',[15,15],'lin');
set(gca,'YDir','reverse');
colormap(MAP); box on;
c_h=colorbar;
xlabel('Radius^2 [cm^2] '); ylabel('Drift Time [\mus] ');
title(strcat(file_id_cp,'.83mKr. dT vs R^2'), 'FontSize', 16,'Interpreter','none');
myfigview(18);
ylabel(c_h,'Count/(\mus\times cm^2)','FontSize',16); 
xlim([0 650]);
saveas(radial_field_fig,strcat('LUX_corrections/',file_id_cp,'/RadialField_',file_id_cp),'jpg');
saveas(radial_field_fig,strcat('LUX_corrections/',file_id_cp,'/RadialField_',file_id_cp),'fig');


%% Calculating the S1 Z-dependence, using both PMT arrays

%edges of histogram binning for fit
s1_bin_min          = 100;
s1_bin_max          = 500;
s1_bin_size         = 8;

%cut to clean up the S1 
cut                 = inrange(s2radius,[rcut_min,rcut_max]) & inrange(s1_phe_both,[s1_bin_min,s1_bin_max]) & inrange(s2_width,[100,1000]) & inrange(s1_width,[30,400]);   

%Do the fit
[means,means_error,mean_index,Fit_s1_both, bins, x, hist_s1_both] = MeasurePulseAreaZDep(s1_bin_min, s1_bin_max, s1_bin_size, s1_phe_both, drift_time, det_edge, cut);

%Fit to the Z depedence
dT_range            = inrange( mean_index,[zcut_min, zcut_max] );
dT_fit              = 0:1:det_edge;    
[P_s1_both S]       = polyfit(mean_index(dT_range),means(dT_range),2);
sigma_P_s1_both     = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df)';    
s1_both_means       = means;
s1_both_means_sigma = means_error;
s1_both_means_bin   = mean_index;

%Correct Z dependenece of S1 pulse areas
s1_z_correction    = polyval(P_s1_both,z_center)./polyval(P_s1_both,drift_time);
s1_phe_both_z      = s1_phe_both.*s1_z_correction;


%Fit a Gaussian to the z-corrected data for tracking over time in the LUG
s1_phe_both_z_fit_binning           = (0:10:1000);
s1_phe_both_z_fit_cut               = inrange(s1_phe_both_z,[min(s1_phe_both_z_fit_binning) max(s1_phe_both_z_fit_binning)]);
s1_phe_both_z_gauss_fit             = fit(s1_phe_both_z_fit_binning.',hist(s1_phe_both_z(s1_phe_both_z_fit_cut),s1_phe_both_z_fit_binning).','gauss1');
s1_phe_both_z_gauss_fit_mean        = s1_phe_both_z_gauss_fit.b1;
s1_phe_both_z_gauss_fit_sigma       = s1_phe_both_z_gauss_fit.c1/sqrt(2);
s1_phe_both_z_gauss_fit_confint     = confint(s1_phe_both_z_gauss_fit,0.68);
s1_phe_both_z_gauss_fit_mean_error  = abs(s1_phe_both_z_gauss_fit_confint(1,2)-s1_phe_both_z_gauss_fit_mean);
s1_phe_both_z_gauss_fit_sigma_error = abs(s1_phe_both_z_gauss_fit_confint(1,3)/sqrt(2) - s1_phe_both_z_gauss_fit_sigma);

%Plot the S1 z dependence
light_correction_fig = figure;
errorbar(mean_index,means,means_error,'.k');
hold on     
[yplot_s1_both] = polyval(P_s1_both,dT_fit);
[center_s1_both, center_s1_both_err] = polyval(P_s1_both,z_center, S); %160[us]*1.51[\mus/us]
plot(dT_fit,yplot_s1_both,'-r','LineWidth',1);    
xlabel('Depth (\mus)', 'FontSize',16); ylabel('Mean (Pulse Area Phe)', 'FontSize', 16);
title(strcat(file_id_cp,'.Kr83m. S1 Mean vs. Z'), 'FontSize', 16,'Interpreter','none');
line([0 det_edge],[center_s1_both center_s1_both],'linestyle','--');
line([z_center; z_center;],[min(means) max(means)],'linestyle','--');
legend('S1\_both Mean', strcat('y=', num2str(P_s1_both(1),'%10.2e'), '*Z^2 + ', num2str(P_s1_both(2),'%10.2e')...
, '*Z + ', num2str(P_s1_both(3),4)),strcat('Det center = ',num2str(center_s1_both,4), ' [Phe]'),'location','northwest')
box on;
myfigview(16);
xlim([0 det_edge]);
ylim([0.8*min(means) 1.15*max(means)]);
saveas(light_correction_fig,strcat('LUX_corrections/',file_id_cp,'/S1_both_LC_',file_id_cp),'jpg');
saveas(light_correction_fig,strcat('LUX_corrections/',file_id_cp,'/S1_both_LC_',file_id_cp),'fig');

%Plot 7 of the Gaussian histograms     
hist_fig = figure;
hold on
plot_step = floor ((bins-1)/7);    
j=1;
t = {num2str(zeros(7,1))};
for i=2:plot_step:bins  %first bin contains anode... skip it. start at i=2.
    if j<= 7
        xfit=100:2:500;
        yfit=Fit_s1_both{i}.a1.*exp(-((xfit-Fit_s1_both{i}.b1)./Fit_s1_both{i}.c1).^2);

        if j==1
            mark = '+'; color = 'k'; t{j}=num2str(mean_index(i),4);
        elseif j==2
            mark = '.'; color = 'b'; t{j}=num2str(mean_index(i),4);    
        elseif j==3
            mark = '*'; color = 'r'; t{j}=num2str(mean_index(i),4);
        elseif j==4
            mark = 'x'; color = 'g'; t{j}=num2str(mean_index(i),4);
        elseif j==5
            mark = '+'; color = 'm'; t{j}=num2str(mean_index(i),4);
        elseif j==6
            mark = '.'; color = 'c'; t{j}=num2str(mean_index(i),4);                
        elseif j==7
            mark = '*'; color = 'y'; t{j}=num2str(mean_index(i),4);                
        end                  

        plot(x,hist_s1_both(:,i),strcat(mark,color));
        plot(xfit,yfit,color);
    end

j=j+1;
end
title(strcat(file_id_cp,'.Kr83m. S1_both Hist Data'), 'FontSize', 16,'Interpreter','none');
legend(strcat(t{1},' \mus'),strcat(t{1},' \mus Fit'),strcat(t{2},' \mus'),strcat(t{2},' \mus Fit')...
    ,strcat(t{3},' \mus'),strcat(t{3},' \mus Fit'),strcat(t{4},' \mus'),strcat(t{4},' \mus Fit')...
    ,strcat(t{5},' \mus'),strcat(t{5},' \mus Fit'),strcat(t{6},' \mus'),strcat(t{6},' \mus Fit')...
    ,strcat(t{7},' \mus'),strcat(t{7},' \mus Fit') ); 
xlabel('Pulse Area Phe','FontSize',16); ylabel('Counts/Pulse-Area-Phe','FontSize',16);
myfigview(16);
xlim([s1_bin_min-100 s1_bin_max+100]);
box on;
saveas(hist_fig,strcat('LUX_corrections/',file_id_cp,'/S1_both_hist_',file_id_cp),'jpg');
saveas(hist_fig,strcat('LUX_corrections/',file_id_cp,'/S1_both_hist_',file_id_cp),'fig');

    
%% Calculating the S1 Z-dependence, using bottom PMT array

%edges of histogram binning for fit
s1_bin_min             = 50;
s1_bin_max             = 400;
s1_bin_size            = 8;

%cut to clean up the S1 
cut                    = inrange(s2radius,[rcut_min,rcut_max]) & inrange(s1_phe_bottom,[s1_bin_min,s1_bin_max]) & inrange(s2_width,[100,1000]) & inrange(s1_width,[30,400]);   

%Do the fit
[means, means_error, mean_index, Fit_s1_bottom, bins, x, hist_s1_bottom] = MeasurePulseAreaZDep(s1_bin_min, s1_bin_max, s1_bin_size, s1_phe_bottom, drift_time, det_edge, cut);

%Fit to the Z depedence
dT_range               = inrange( mean_index,[zcut_min, zcut_max] );
dT_fit                 = 0:1:det_edge;    
[P_s1_bottom S]        = polyfit(mean_index(dT_range),means(dT_range),2);
sigma_P_s1_bottom      = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df)';   
s1_bottom_means        = means;
s1_bottom_means_sigma  = means_error;
s1_bottom_means_bin    = mean_index;

%Correct Z dependenece of S1 pulse areas
s1_z_correction        = polyval(P_s1_bottom,z_center)./polyval(P_s1_bottom,drift_time);
s1_phe_bottom_z        = s1_phe_bottom.*s1_z_correction;


%Plots the S1 z dependence
light_correction_fig_bottom = figure;
errorbar(mean_index,means,means_error,'.k');
hold on     
[yplot_s1_bottom] = polyval(P_s1_bottom,dT_fit);
[center_s1_bottom, center_s1_bottom_err] = polyval(P_s1_bottom,z_center,S); %160[us]*1.51[\mus/us]
plot(dT_fit,yplot_s1_bottom,'-r','LineWidth',1);    
xlabel('Depth (\mus)', 'FontSize',16); ylabel('Mean (Pulse Area Phe)', 'FontSize', 16);
title(strcat(file_id_cp,'.Kr83m. S1_bottom Mean vs. Z'), 'FontSize', 16,'Interpreter','none');
line([0 det_edge],[center_s1_bottom center_s1_bottom],'linestyle','--');
line([z_center; z_center;],[min(means) max(means)],'linestyle','--');
legend('S1\_bottom Mean', strcat('y=', num2str(P_s1_bottom(1),'%10.2e'), '*Z^2 + ', num2str(P_s1_bottom(2),'%10.2e')...
    , '*Z + ', num2str(P_s1_bottom(3),4)),strcat('Det center = ',num2str(center_s1_bottom,4), ' [Phe]'),'location','northwest')
box on;
myfigview(16);
xlim([0 det_edge]);
ylim([0.8*min(means) 1.15*max(means)]);
saveas(light_correction_fig_bottom,strcat('LUX_corrections/',file_id_cp,'/S1_bottom_LC_',file_id_cp),'jpg');
saveas(light_correction_fig_bottom,strcat('LUX_corrections/',file_id_cp,'/S1_bottom_LC_',file_id_cp),'fig');
    
%Plot 7 of the histograms     
hist_s1_bottom_fig = figure;
hold on
plot_step = floor((bins-1)/7);    
j=1;
t = {num2str(zeros(7,1))};
for i=2:plot_step:bins  %first bin contains anode... skip it. start at i=2.
    if j<=7
        xfit=50:2:500;
        yfit=Fit_s1_bottom{i}.a1.*exp(-((xfit-Fit_s1_bottom{i}.b1)./Fit_s1_bottom{i}.c1).^2);
        if j==1
          mark = '+'; color = 'k'; t{j}=num2str(mean_index(i),4);
        elseif j==2
            mark = '.'; color = 'b'; t{j}=num2str(mean_index(i),4);    
        elseif j==3
            mark = '*'; color = 'r'; t{j}=num2str(mean_index(i),4);
        elseif j==4
            mark = 'x'; color = 'g'; t{j}=num2str(mean_index(i),4);
        elseif j==5
            mark = '+'; color = 'm'; t{j}=num2str(mean_index(i),4);
        elseif j==6
            mark = '.'; color = 'c'; t{j}=num2str(mean_index(i),4);                
        elseif j==7
           mark = '*'; color = 'y'; t{j}=num2str(mean_index(i),4);                
        end                  
        plot(x,hist_s1_bottom(:,i),strcat(mark,color));
        plot(xfit,yfit,color);
    end
    j=j+1;
end
title(strcat(file_id_cp,'.Kr83m. S1_bottom Hist Data'), 'FontSize', 16,'Interpreter','none');
legend(strcat(t{1},' \mus'),strcat(t{1},' \mus Fit'),strcat(t{2},' \mus'),strcat(t{2},' \mus Fit')...
    ,strcat(t{3},' \mus'),strcat(t{3},' \mus Fit'),strcat(t{4},' \mus'),strcat(t{4},' \mus Fit')...
    ,strcat(t{5},' \mus'),strcat(t{5},' \mus Fit'),strcat(t{6},' \mus'),strcat(t{6},' \mus Fit')...
    ,strcat(t{7},' \mus'),strcat(t{7},' \mus Fit') ); 
xlabel('Pulse Area Phe','FontSize',16); ylabel('Counts/Pulse-Area-Phe','FontSize',16);
myfigview(16);
xlim([s1_bin_min-50 s1_bin_max+50]);
box on;
saveas(hist_s1_bottom_fig,strcat('LUX_corrections/',file_id_cp,'/S1_bottom_hist_',file_id_cp),'jpg');
saveas(hist_s1_bottom_fig,strcat('LUX_corrections/',file_id_cp,'/S1_bottom_hist_',file_id_cp),'fig');


%% Calculating the S1 Z-dependence, using top PMT array

dT_range                           = inrange( mean_index,[zcut_min, zcut_max] );
[P_s1_top S]                       = polyfit(s1_both_means_bin(dT_range),s1_both_means(dT_range)-s1_bottom_means(dT_range),2);
sigma_P_s1_top                     = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df)';
[center_s1_top, center_s1_top_err] = polyval(P_s1_top,z_center,S); %160[us]*1.51[\mus/us]


%% Calculating the S2 Z-dependence, using both PMT arrays

%edges of histogram binning for fit
s2_bin_min                         = 1000;
s2_bin_max                         = 40000;
s2_bin_size                        = 500;

%cut to clean up the S2 pulses
cut                                = inrange(s2radius,[rcut_min,rcut_max]) & inrange(s2_phe_both,[s2_bin_min,s2_bin_max]) & inrange(s2_width,[100,1000]) & inrange(s1_width,[30,400]);

%Do the fit
[means,means_error,mean_index,Fit_s2_both, bins, x2, hist_s2_both] = MeasurePulseAreaZDep(s2_bin_min, s2_bin_max, s2_bin_size, s2_phe_both, drift_time, det_edge, cut);

%Fit to the Z depedence
dT_range                           = inrange( mean_index,[5+mean_index(end)-mean_index(end-1),500] );%from 50 to 300us for the fit
s2_both_norm_z_index               = mean_index(dT_range);
s2_both_norm_z_means               = means(dT_range);
s2_at_top_both                     = RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,4); %normalize to right below the gate
s2_both_norm_z                     = RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,4)./s2_both_norm_z_means;
s2_both_means                      = means;
s2_both_means_sigma                = means_error;
s2_both_means_bin                  = mean_index;        

%%Tracker for pseudo-lifetime
s2_at_bottom_both                  = RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,320);
pseudo_lifetime_both               = -320/log(s2_at_bottom_both/s2_at_top_both);
dT_fit                             = 0:1:det_edge;    
Fit_EL                             = fit(mean_index(dT_range),means(dT_range),'exp1');
yfitEL                             = Fit_EL.a.*exp(Fit_EL.b.*dT_fit);
e_lifetime_both                    = -1/Fit_EL.b;
s2_z0_both                         = Fit_EL.a;
b                                  = confint(Fit_EL, 0.683);%1 sigma matrix
sigma_e_lifetime_both              =(1/b(1,2)-1/b(2,2))/2;
sigma_s2_z0_both                   =(b(2,1)-b(1,1))/2;
 
%Make plot of Z dependence
lifetime_fig_both = figure;
errorbar(mean_index(2:end),means(2:end),means_error(2:end),'.k');  %first bin is anode, skip it
hold on     
plot(dT_fit,yfitEL,'-r','LineWidth',1);   
xlabel('Time (\mus)', 'FontSize',16); ylabel('Mean S2 (Pulse Area Phe)', 'FontSize', 16);
title(strcat(file_id_cp,'.Kr83m. S2_both Electron Lifetime'), 'FontSize', 16,'Interpreter','none');
legend('S2\_both Mean', strcat( '\lambda= ', num2str(e_lifetime_both,4), ' \pm ', num2str(sigma_e_lifetime_both,2), ' \mus' )...
    ,'location','northeast');
myfigview(16);
xlim([0 det_edge]);
saveas(lifetime_fig_both,strcat('LUX_corrections/',file_id_cp,'/S2_both_lifetime_',file_id_cp),'jpg');
saveas(lifetime_fig_both,strcat('LUX_corrections/',file_id_cp,'/S2_both_lifetime_',file_id_cp),'fig');

%Plot 7 of the Gaussian fit histograms     
hist_s2_both_fig = figure;
hold on
plot_step = floor((bins-1)/7);    
j=1;
t = {num2str(zeros(7,1))};
for i=2:plot_step:bins  %first bin contains anode... skip it. start at i=2.
    if j<=7
    xfit=s2_bin_min:s2_bin_size:s2_bin_max;
    yfit=Fit_s2_both{i}.a1.*exp(-((xfit-Fit_s2_both{i}.b1)./Fit_s2_both{i}.c1).^2);

        if j==1
            mark = '+'; color = 'k'; t{j}=num2str(mean_index(i),4);
        elseif j==2
            mark = '.'; color = 'b'; t{j}=num2str(mean_index(i),4);    
        elseif j==3
            mark = '*'; color = 'r'; t{j}=num2str(mean_index(i),4);
        elseif j==4
            mark = 'x'; color = 'g'; t{j}=num2str(mean_index(i),4);
        elseif j==5
            mark = '+'; color = 'm'; t{j}=num2str(mean_index(i),4);
        elseif j==6
            mark = '.'; color = 'c'; t{j}=num2str(mean_index(i),4);                
        elseif j==7
            mark = '*'; color = 'y'; t{j}=num2str(mean_index(i),4);                
        end                  

    plot(x2,hist_s2_both(:,i)./d.livetime_sec,strcat(mark,color));
    plot(xfit,yfit./d.livetime_sec,color);
    end
    j=j+1;
end
    title(strcat(file_id_cp,'.Kr83m. S2_both'), 'FontSize', 16,'Interpreter','none');
    legend(strcat(t{1},' \mus'),strcat(t{1},' \mus Fit'),strcat(t{2},' \mus'),strcat(t{2},' \mus Fit')...
        ,strcat(t{3},' \mus'),strcat(t{3},' \mus Fit'),strcat(t{4},' \mus'),strcat(t{4},' \mus Fit')...
        ,strcat(t{5},' \mus'),strcat(t{5},' \mus Fit'),strcat(t{6},' \mus'),strcat(t{6},' \mus Fit')...
        ,strcat(t{7},' \mus'),strcat(t{7},' \mus Fit') ); 
    xlabel('Pulse Area Phe','FontSize',16); ylabel('Counts/Phe/sec','FontSize',16);
myfigview(16);
xlim([s2_bin_min-50 s2_bin_max+500]);
saveas(hist_s2_both_fig,strcat('LUX_corrections/',file_id_cp,'/S2_both_hist_EL_',file_id_cp),'jpg');
saveas(hist_s2_both_fig,strcat('LUX_corrections/',file_id_cp,'/S2_both_hist_EL_',file_id_cp),'fig');
    
%Correct Z dependenece of S1 pulse areas
s2_phe_both_z                       = s2_phe_both.*RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z, drift_time);

%Fit a Gaussian to the z-corrected data for tracking over time in the LUG 
Fit_s2_both_z                       = fit(x2,hist(s2_phe_both_z(inrange(s2_phe_both,[s2_bin_min,s2_bin_max])),x2)'/s2_bin_size,'gauss1');%remove bin edges
s2_both_center_z                    = Fit_s2_both_z.b1;
temp_conf                           = confint(Fit_s2_both_z,0.68);
s2_both_center_error_z              = abs(s2_both_center_z-temp_conf(1,2));

s2_phe_both_z_gauss_fit_mean        = Fit_s2_both_z.b1;
s2_phe_both_z_gauss_fit_sigma       = Fit_s2_both_z.c1/sqrt(2);
s2_phe_both_z_gauss_fit_mean_error  = abs(temp_conf(1,2)-s2_phe_both_z_gauss_fit_mean);
s2_phe_both_z_gauss_fit_sigma_error = abs(temp_conf(1,3)/sqrt(2)-s2_phe_both_z_gauss_fit_sigma);


%% Calculating the S2 Z-dependence, using bottom PMT arrays

%edges of histogram binning for fit
s2_bin_min                         = 50;
s2_bin_max                         = 10000;
s2_bin_size                        = 300;

%cut to clean up the S2 pulses
cut                                = inrange(s2radius,[rcut_min,rcut_max]) & inrange(s2_phe_bottom,[s2_bin_min,s2_bin_max]) & inrange(s2_width,[100,1000]) & inrange(s1_width,[30,400]) ;

%Do the fit
[means,means_error,mean_index,Fit_s2_bottom, bins, x2, hist_s2_bottom] = MeasurePulseAreaZDep(s2_bin_min, s2_bin_max, s2_bin_size, s2_phe_bottom, drift_time, det_edge, cut);

%Fit to the Z depedence
dT_range                           = inrange( mean_index,[5+mean_index(end)-mean_index(end-1),500] );%from 50 to 300us for the fit
s2_bottom_norm_z_index             = mean_index(dT_range);
s2_bottom_norm_z_means             = means(dT_range);
s2_at_top_bottom                   = RK_1DCubicInterp_LinearExtrap(s2_bottom_norm_z_index,s2_bottom_norm_z_means,4); %normalize to right below the gate
s2_bottom_norm_z                   = RK_1DCubicInterp_LinearExtrap(s2_bottom_norm_z_index,s2_bottom_norm_z_means,4)./s2_bottom_norm_z_means;
s2_bottom_means                    = means;
s2_bottom_means_sigma              = means_error;
s2_bottom_means_bin                = mean_index;   

%Tracker for pseudo-lifetime IQ
s2_at_bottom_bottom                = RK_1DCubicInterp_LinearExtrap(s2_bottom_norm_z_index,s2_bottom_norm_z_means,320);
pseudo_lifetime_bottom             = -320/log(s2_at_bottom_bottom/s2_at_top_bottom);
dT_fit                             = 0:1:det_edge;    
Fit_EL                             = fit(mean_index(dT_range),means(dT_range),'exp1');
yfitEL                             = Fit_EL.a.*exp(Fit_EL.b.*dT_fit);
e_lifetime                         = -1/Fit_EL.b;
s2_z0_bottom                       = Fit_EL.a;
b                                  = confint(Fit_EL, 0.683);%1 sigma matrix
sigma_e_lifetime                   = (1/b(1,2)-1/b(2,2))/2;
sigma_s2_z0_bottom                 = (b(2,1)-b(1,1))/2;


%Make plot of Z dependence
lifetime_fig_bottom = figure;
errorbar(mean_index(2:end),means(2:end),means_error(2:end),'.k');  %first bin is anode, skip it
hold on     
plot(dT_fit,yfitEL,'-r','LineWidth',1);   
xlabel('Depth (\mus)', 'FontSize',16); ylabel('Mean S2 (Pulse Area Phe)', 'FontSize', 16);
title(strcat(file_id_cp,'.Kr83m. S2_bottom Electron Lifetime'), 'FontSize', 16,'Interpreter','none');
legend('S2\_bottom Mean', strcat( '\lambda= ', num2str(e_lifetime,4), ' \pm ', num2str(sigma_e_lifetime,2), ' \mus' )...
    ,'location','northeast');
myfigview(16);
xlim([0 det_edge]);
saveas(lifetime_fig_bottom,strcat('LUX_corrections/',file_id_cp,'/S2_bottom_lifetime_',file_id_cp),'jpg');
saveas(lifetime_fig_bottom,strcat('LUX_corrections/',file_id_cp,'/S2_bottom_lifetime_',file_id_cp),'fig');

%Plot 7 of the Gaussian fit histograms     
hist_s2_bottom_fig = figure;
hold on
plot_step = floor((bins-1)/7);    
j=1;
t = {num2str(zeros(7,1))};
for i=2:plot_step:bins  %first bin contains anode... skip it. start at i=2.
    if j<=7
    xfit=100:50:10000;
    yfit=Fit_s2_bottom{i}.a1.*exp(-((xfit-Fit_s2_bottom{i}.b1)./Fit_s2_bottom{i}.c1).^2);

        if j==1
            mark = '+'; color = 'k'; t{j}=num2str(mean_index(i),4);
        elseif j==2
            mark = '.'; color = 'b'; t{j}=num2str(mean_index(i),4);    
        elseif j==3
            mark = '*'; color = 'r'; t{j}=num2str(mean_index(i),4);
        elseif j==4
            mark = 'x'; color = 'g'; t{j}=num2str(mean_index(i),4);
        elseif j==5
            mark = '+'; color = 'm'; t{j}=num2str(mean_index(i),4);
        elseif j==6
            mark = '.'; color = 'c'; t{j}=num2str(mean_index(i),4);                
        elseif j==7
            mark = '*'; color = 'y'; t{j}=num2str(mean_index(i),4);                
        end                  

    plot(x2,hist_s2_bottom(:,i)./d.livetime_sec,strcat(mark,color));
    plot(xfit,yfit/d.livetime_sec,color);
    end
    j=j+1;
end
title(strcat(file_id_cp,'.Kr83m.S2_Bottom'), 'FontSize', 16,'Interpreter','none');
legend(strcat(t{1},' \mus'),strcat(t{1},' \mus Fit'),strcat(t{2},' \mus'),strcat(t{2},' \mus Fit')...
    ,strcat(t{3},' \mus'),strcat(t{3},' \mus Fit'),strcat(t{4},' \mus'),strcat(t{4},' \mus Fit')...
    ,strcat(t{5},' \mus'),strcat(t{5},' \mus Fit'),strcat(t{6},' \mus'),strcat(t{6},' \mus Fit')...
    ,strcat(t{7},' \mus'),strcat(t{7},' \mus Fit')); 
xlabel('Pulse Area Phe','FontSize',16); ylabel('Counts/Phe/sec','FontSize',16);
myfigview(16);
xlim([s2_bin_min-50 s2_bin_max+500]);
saveas(hist_s2_bottom_fig,strcat('LUX_corrections/',file_id_cp,'/S2_bottom_hist_EL_',file_id_cp),'jpg');
saveas(hist_s2_bottom_fig,strcat('LUX_corrections/',file_id_cp,'/S2_bottom_hist_EL_',file_id_cp),'fig');

%Corrected Z dependence in the data
s2_phe_bottom_z                   = s2_phe_bottom.*RK_1DCubicInterp_LinearExtrap(s2_bottom_norm_z_index,s2_bottom_norm_z, drift_time);


%Fit a Gaussian for tracking in the LUG
Fit_s2_bot_z                     = fit(x2,hist(s2_phe_bottom_z,x2)'/s2_bin_size,'gauss1');%remove bin edges
s2_bot_center_z                  = Fit_s2_bot_z.b1;
temp_conf                        = confint(Fit_s2_bot_z);
s2_bot_center_error_z            = abs(s2_bot_center_z-temp_conf(1,2));


%% Calculating the S2 XY-dependence, using both PMT arrays

%Define variables for making the XY map
s2_bin_min                                             = 1000;
s2_bin_max                                             = 40000;
s2_bin_size                                            = 500;

cleanup_cut                                            = inrange(drift_time,[zcut_min,zcut_max]) & inrange(s1_phe_both_z,[s1area_bound_min,s1area_bound_max]) & ...
                                                         inrange(s2_phe_both_z,[s2_bin_min,s2_bin_max]) & inrange(s2radius,[rcut_min,rcut_max]) & ...
                                                         inrange(s2_width,[100,1000]) & inrange(s1_width,[30,400]);
%Make the XY map
[mean_S2_both_xy, Sigma_S2_both]                       = MeasurePulseAreaXYDep(s2_bin_min, s2_bin_max, s2_bin_size, S2_xbin_min,...
                                                            S2_xbin_max, S2_ybin_min, S2_ybin_max, S2xybinsize, s2_phe_both_z, s2x, s2y, cleanup_cut);

%Calculate corrections matrix
[center, sigma_center, norm_S2_all, sigma_norm_S2_all] = ProducePulseAreaXYNorms(s2xbins, s2ybins, mean_S2_both_xy, Sigma_S2_both, x_center, y_center);

%Define variables for tracking mean at center over time in the LUG
s2_both_center_xyz                                     = center;
s2_both_center_xyz_error                               = sigma_center;

%Correct the Z corrected data in XY
s2xy_correction                                        = interp2(s2xbins,s2ybins,norm_S2_all,s2x,s2y,'cubic');
s2xy_correction(find(s2xy_correction==0))              = 1.0; %#ok<FNDSB>
s2xy_correction(find(isnan(s2xy_correction)))          = 1.0; %#ok<FNDSB>
s2_phe_both_xyz                                        = s2_phe_both_z.*s2xy_correction;

%Fit a Gaussian tothe corrected data for tracking in the LUG
s2_phe_both_xyz_fit_bins                               = (0:100:70000);
s2_phe_both_xyz_fit_cut                                = inrange(s2_phe_both_xyz,[min(s2_phe_both_xyz_fit_bins) max(s2_phe_both_xyz_fit_bins)]);
s2_phe_both_xyz_gauss_fit                              = fit(s2_phe_both_xyz_fit_bins.',hist(s2_phe_both_xyz(s2_phe_both_xyz_fit_cut),s2_phe_both_xyz_fit_bins).','gauss1');
s2_phe_both_xyz_gauss_fit_mean                         = s2_phe_both_xyz_gauss_fit.b1;
s2_phe_both_xyz_gauss_fit_sigma                        = s2_phe_both_xyz_gauss_fit.c1/sqrt(2);
s2_phe_both_xyz_gauss_fit_confint                      = confint(s2_phe_both_xyz_gauss_fit,0.68);
s2_phe_both_xyz_gauss_fit_mean_error                   = abs(s2_phe_both_xyz_gauss_fit_confint(1,2)-s2_phe_both_xyz_gauss_fit_mean);
s2_phe_both_xyz_gauss_fit_sigma_error                  = abs(s2_phe_both_xyz_gauss_fit_confint(1,3)/sqrt(2)-s2_phe_both_xyz_gauss_fit_sigma);

%Plot the mean S2_XY both
s2_xy_both_fig = figure;
color_range_max=max(max(mean_S2_both_xy));
color_range_min=min(min(mean_S2_both_xy(mean_S2_both_xy>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s2xbins,s2ybins,mean_S2_both_xy,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. S2 Mean (Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);
saveas(s2_xy_both_fig,strcat('LUX_corrections/',file_id_cp,'/S2_XY_both_',file_id_cp),'jpg');
saveas(s2_xy_both_fig,strcat('LUX_corrections/',file_id_cp,'/S2_XY_both_',file_id_cp),'fig');

%Plot the 1 sigma error of the mean
sigma_s2_xy_both_fig = figure;
color_range_max=max(max(Sigma_S2_both));
color_range_min=min(min(Sigma_S2_both(Sigma_S2_both>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s2xbins,s2ybins,Sigma_S2_both,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. 1-Sigma S2 Mean(Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);
saveas(sigma_s2_xy_both_fig,strcat('LUX_corrections/',file_id_cp,'/Sigma_S2_XY_both_',file_id_cp),'jpg');
saveas(sigma_s2_xy_both_fig,strcat('LUX_corrections/',file_id_cp,'/Sigma_S2_XY_both_',file_id_cp),'fig');


%% Calculating the S2 XY-dependence, using bottom PMT array

%Define variables for making the XY map
s2_bin_min                                             = 50;
s2_bin_max                                             = 10000;
s2_bin_size                                            = 300;

cleanup_cut                                            = inrange(drift_time,[zcut_min,zcut_max]) & inrange(s1_phe_both_z,[100,s1area_bound_max]) & ...
                                                         inrange(s2_phe_bottom_z,[s2_bin_min,s2_bin_max]) & inrange(s2radius,[rcut_min,rcut_max]) & ...
                                                         inrange(s2_width,[100,1000]) & inrange(s1_width,[30,400]);

%Make the XY map
[mean_S2_bottom_xy, Sigma_S2_bottom, Count_S2_bottom]  = MeasurePulseAreaXYDep(s2_bin_min, s2_bin_max, s2_bin_size, S2_xbin_min,...
                                                            S2_xbin_max, S2_ybin_min, S2_ybin_max, S2xybinsize, s2_phe_bottom_z, s2x, s2y, cleanup_cut);

%Calculate corrections matrix
[center, sigma_center, norm_S2_bot, sigma_norm_S2_bot] = ProducePulseAreaXYNorms(s2xbins, s2ybins, mean_S2_bottom_xy, Sigma_S2_bottom, x_center, y_center);

%Define variables for tracking mean at center over time in the LUG
s2_bot_center_xyz                                      = center;
s2_bot_center_xyz_error                                = sigma_center;


%Plot the mean S2_XY bottom
s2_xy_bottom_fig = figure;
color_range_max=  max(max(mean_S2_bottom_xy)) ;
color_range_min= min(min(mean_S2_bottom_xy(mean_S2_bottom_xy>0))) ;
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s2xbins,s2ybins,mean_S2_bottom_xy,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp,'. S2 Mean (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);
saveas(s2_xy_bottom_fig,strcat('LUX_corrections/',file_id_cp,'/S2_XY_bottom_',file_id_cp),'jpg');
saveas(s2_xy_bottom_fig,strcat('LUX_corrections/',file_id_cp,'/S2_XY_bottom_',file_id_cp),'fig');


%Plot the 1 sigma error of the mean
sigma_s2_xy_bottom_fig = figure;
color_range_max=max(max(Sigma_S2_bottom));
color_range_min=min(min(Sigma_S2_bottom(Sigma_S2_bottom>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s2xbins,s2ybins,Sigma_S2_bottom,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. 1-Sigma S2 Mean (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);
saveas(sigma_s2_xy_bottom_fig,strcat('LUX_corrections/',file_id_cp,'/Sigma_S2_XY_bottom_',file_id_cp),'jpg');
saveas(sigma_s2_xy_bottom_fig,strcat('LUX_corrections/',file_id_cp,'/Sigma_S2_XY_bottom_',file_id_cp),'fig');


%% Calculating the S1 XY-dependence, using both PMT arrays

%Define variables for making the XY map
s1_bin_min                                             = 100;
s1_bin_max                                             = 500;
s1_bin_size                                            = 8;

cleanup_cut                                            = inrange(drift_time,[zcut_min,zcut_max]) & inrange(s1_phe_both_z,[s1_bin_min,s1_bin_max]) & ...
                                                         inrange(s2_phe_both_z,[1000,s2area_bound_max]) & inrange(s2radius,[rcut_min,rcut_max]) & ...
                                                         inrange(s2_width,[100 1000]) & inrange(s1_width,[30 400]);
                                                     
%Make the XY map
[mean_S1_both_xy, Sigma_S1_both, Count_S1_both]        = MeasurePulseAreaXYDep(s1_bin_min, s1_bin_max, s1_bin_size, S1_xbin_min,...
                                                            S1_xbin_max, S1_ybin_min, S1_ybin_max, S1xybinsize, s1_phe_both_z, s2x, s2y, cleanup_cut);

%Calculate corrections matrix
[center, sigma_center, norm_S1_all, sigma_norm_S1_all] = ProducePulseAreaXYNorms(s1xbins, s1ybins, mean_S1_both_xy, Sigma_S1_both, x_center, y_center);

%Define variables for tracking mean at center over time in the LUG
s1_both_center_xyz                                     = center;
s1_both_center_xyz_error                               = sigma_center;

%Correct the Z corrected data in XY
s1xy_correction                                        = interp2(s1xbins,s1ybins,norm_S1_all,s2x,s2y,'cubic');
s1xy_correction(find(s1xy_correction==0))              = 1.0;   %#ok<FNDSB>
s1xy_correction(find(isnan(s1xy_correction)))          = 1.0;   %#ok<FNDSB>
s1_phe_both_xyz                                        = s1_phe_both_z.*s1xy_correction;

%Fit a Gaussian to the corrected data for tracking in the LUG
s1_phe_both_xyz_fit_binning                            = (0:10:1000);
s1_phe_both_xyz_fit_cut                                = inrange(s1_phe_both_xyz,[min(s1_phe_both_xyz_fit_binning) max(s1_phe_both_xyz_fit_binning)]);
s1_phe_both_xyz_gauss_fit                              = fit(s1_phe_both_xyz_fit_binning.',hist(s1_phe_both_xyz(s1_phe_both_xyz_fit_cut),s1_phe_both_xyz_fit_binning).','gauss1');
s1_phe_both_xyz_gauss_fit_mean                         = s1_phe_both_xyz_gauss_fit.b1;
s1_phe_both_xyz_gauss_fit_sigma                        = s1_phe_both_xyz_gauss_fit.c1/sqrt(2);
s1_phe_both_xyz_gauss_fit_confint                      = confint(s1_phe_both_xyz_gauss_fit,0.68);
s1_phe_both_xyz_gauss_fit_mean_error                   = abs(s1_phe_both_xyz_gauss_fit_confint(1,2)-s1_phe_both_xyz_gauss_fit_mean);
s1_phe_both_xyz_gauss_fit_sigma_error                  = abs(s1_phe_both_xyz_gauss_fit_confint(1,3)/sqrt(2)-s1_phe_both_xyz_gauss_fit_sigma);

%Plot the mean S1_XY both
s1_xy_both_fig = figure;
color_range_max=max(max(mean_S1_both_xy));
color_range_min=min(min(mean_S1_both_xy(mean_S1_both_xy>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s1xbins,s1ybins,mean_S1_both_xy,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp,'. S1 Mean (Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);
saveas(s1_xy_both_fig,strcat('LUX_corrections/',file_id_cp,'/S1_XY_both_',file_id_cp),'jpg');
saveas(s1_xy_both_fig,strcat('LUX_corrections/',file_id_cp,'/S1_XY_both_',file_id_cp),'fig');

%Plot the 1 sigma error of the mean
sigma_s1_xy_both_fig = figure;
color_range_max=max(max(Sigma_S1_both));
color_range_min=min(min(Sigma_S1_both(Sigma_S1_both>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s1xbins,s1ybins,Sigma_S1_both,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. 1-Sigma S1 Mean(Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);
saveas(sigma_s1_xy_both_fig,strcat('LUX_corrections/',file_id_cp,'/Sigma_S1_XY_both_',file_id_cp),'jpg');
saveas(sigma_s1_xy_both_fig,strcat('LUX_corrections/',file_id_cp,'/Sigma_S1_XY_both_',file_id_cp),'fig');

%Plot Number of events vs. XY
s1_xy_both_count_fig = figure;
color_range_max=max(max(Count_S1_both));
color_range_min=min(min(Count_S1_both(Count_S1_both>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s1xbins,s1ybins,Count_S1_both,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp,'. Count vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'Count','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16)
saveas(s1_xy_both_count_fig,strcat('LUX_corrections/',file_id_cp,'/S1_Count_XY_both_',file_id_cp),'jpg'); 
saveas(s1_xy_both_count_fig,strcat('LUX_corrections/',file_id_cp,'/S1_Count_XY_both_',file_id_cp),'fig');


%% Calculating the S1 XY-dependence, using bottom PMT arrays

%Define variables for making the XY map
s1_bin_min                                             = 50;
s1_bin_max                                             = 400;
s1_bin_size                                            = 8;

cleanup_cut                                            = inrange(drift_time,[zcut_min,zcut_max]) & inrange(s1_phe_both_z,[s1_bin_min,s1_bin_max]) & ...
                                                         inrange(s2_phe_both_z,[1000,s2area_bound_max]) & inrange(s2radius,[rcut_min,rcut_max]) & ...
                                                         inrange(s2_width,[100 1000]) & inrange(s1_width,[30 400]);
                                                     
%Make the XY map
[mean_S1_bottom_xy, Sigma_S1_bottom]                   = MeasurePulseAreaXYDep(s1_bin_min, s1_bin_max, s1_bin_size, S1_xbin_min,...
                                                            S1_xbin_max, S1_ybin_min, S1_ybin_max, S1xybinsize, s1_phe_bottom_z, s2x, s2y, cleanup_cut);

%Calculate corrections matrix
[center, sigma_center, norm_S1_bot, sigma_norm_S1_bot] = ProducePulseAreaXYNorms(s1xbins, s1ybins, mean_S1_bottom_xy, Sigma_S1_bottom, x_center, y_center);

%Define variables for tracking mean at center over time in the LUG
s1_bot_center_xyz                                      = center;
s1_bot_center_xyz_error                                = sigma_center;


%Plot the mean S1_XY bottom
s1_xy_bottom_fig = figure;
color_range_max=max(max(mean_S1_bottom_xy));
color_range_min=min(min(mean_S1_bottom_xy(mean_S1_bottom_xy>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s1xbins,s1ybins,mean_S1_bottom_xy,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp,'. S1 Mean (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);
saveas(s1_xy_bottom_fig,strcat('LUX_corrections/',file_id_cp,'/S1_XY_bottom_',file_id_cp),'jpg');
saveas(s1_xy_bottom_fig,strcat('LUX_corrections/',file_id_cp,'/S1_XY_bottom_',file_id_cp),'fig');

%Plot the 1 sigma error of the mean
sigma_s1_xy_bottom_fig = figure;
color_range_max=max(max(Sigma_S1_bottom));
color_range_min=min(min(Sigma_S1_bottom(Sigma_S1_bottom>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s1xbins,s1ybins,Sigma_S1_bottom,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. 1-Sigma S1 Mean(Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);
saveas(sigma_s1_xy_bottom_fig,strcat('LUX_corrections/',file_id_cp,'/Sigma_S1_XY_bottom_',file_id_cp),'jpg');
saveas(sigma_s1_xy_bottom_fig,strcat('LUX_corrections/',file_id_cp,'/Sigma_S1_XY_bottom_',file_id_cp),'fig');


%% Calculating the S1 XY-dependence, using top PMT array  

%Set up top array quantities
mean_S1_top_xy                                         = mean_S1_both_xy-mean_S1_bottom_xy;
Sigma_S1_top                                           = sqrt(Sigma_S1_both.^2+Sigma_S1_bottom.^2);

%Make the correction maps
[center, sigma_center, norm_S1_top, sigma_norm_S1_top] = ProducePulseAreaXYNorms(s1xbins, s1ybins, mean_S1_top_xy, Sigma_S1_top, x_center, y_center);

%Create variables for tracking over time in LUG
s1_top_center_xyz                                      = center;
s1_top_center_xyz_error                                = sigma_center;


%% Mapping S1 in XYZ (Only if there are enough statistics)

if Kr_events>100000
    KrypCal_S1_3D(s1_phe_both,s1_phe_bottom,s2_phe_both,s2x,s2y,drift_time,x_center,y_center,z_center,det_edge,user_name,file_id_cp,dp_version, algorithm_name, submit);
end    
        
    
%% Calculate S2 Pulse Width Z Dependence (Pulse_End_Samples -Pulse_Start_Samples)

%Redefine s2_width as pulse_end_samples - pulse_start_samples 
s2_width       = d.pulse_end_samples(s2_single_cut)-d.pulse_start_samples(s2_single_cut);

%Make cut to clean up data and fit for the Z dependence
cleanup_cut    = inrange(s1_phe_both_z,[100,s1area_bound_max]) & inrange(s2_phe_both_z,[1000,s2area_bound_max]) & inrange(s2radius,[rcut_min,rcut_max]) ...
                 & inrange(s2_width,[450 1200]) & inrange(s1_width,[30 400]);  

[mean_index_s2_z_width_full, means_s2_z_width_full, means_error_s2_z_width_full, P_s2_width_full, sigma_P_s2_width_full, dT_fit] = ...
                 MeasureS2PulseWidthZDep(s2_width, drift_time, det_edge, cleanup_cut); %#ok<ASGLU>

%apply Z dep S2 width correction
norm_s2_width  = RK_1DCubicInterp_LinearExtrap(mean_index_s2_z_width_full,means_s2_z_width_full,4); %normalize to 0 us, top
s2_width_z     = s2_width.*norm_s2_width./RK_1DCubicInterp_LinearExtrap(mean_index_s2_z_width_full,means_s2_z_width_full,drift_time);

%Make plot of Z dependence
s2_z_width_full_fig = figure;
errorbar( mean_index_s2_z_width_full(2:end),means_s2_z_width_full(2:end),means_error_s2_z_width_full(2:end),'.k');  %first bin is anode, skip it
hold on    
[yplot_s2_width] = polyval(P_s2_width_full,dT_fit);
plot(dT_fit,yplot_s2_width,'-r','LineWidth',1);   
xlabel('Depth (\mus)', 'FontSize',16); ylabel('S2 Width (End-Start Samples)', 'FontSize', 16);
title(strcat(file_id_cp,'.Kr83m. S2 width vs. Z'), 'FontSize', 16,'Interpreter','none');
legend('S2 Width Mean', strcat('y=', num2str(P_s2_width_full(1),'%10.2e'), '*Z^3 + ', num2str(P_s2_width_full(2),'%10.2e')...
    , '*Z^2 + ', num2str(P_s2_width_full(3),'%10.2e'), '*Z + ',num2str(P_s2_width_full(4),4)),'location','northwest')
myfigview(16);
xlim([0 det_edge]);
saveas(s2_z_width_full_fig,strcat('LUX_corrections/',file_id_cp,'/S2_z_width_end_',file_id_cp),'jpg');
saveas(s2_z_width_full_fig,strcat('LUX_corrections/',file_id_cp,'/S2_z_width_end_',file_id_cp),'fig');


%% Calculate XY dependence of S2 Pulse Width (Pulse_End_Samples -Pulse_Start_Samples)

%Make cut to clean up data and fit for the Z dependence   
cleanup_cut = inrange(drift_time,[zcut_min,zcut_max])  & inrange(s2_width_z,[200 1400])& inrange(s2radius,[rcut_min,rcut_max]) & ...
              inrange(s2_phe_both_z,[1000,s2area_bound_max]) & inrange(s2radius,[0 24]) ;

[s2_xbins, s2_ybins, means_s2_xy_width_full, sigma_s2_xy_width_full ] = MeasureS2PulseWidthXYDep(s2x, s2y, s2_width_z, S2_xbin_min, S2_xbin_max, S2_ybin_min, S2_ybin_max, S2xybinsize, cleanup_cut); %#ok<ASGLU>

%Make plot of XY dependence
s2_width_xy_end_fig = figure;
color_range_max=max(max(means_s2_xy_width_full));
color_range_min=min(min(means_s2_xy_width_full(means_s2_xy_width_full>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s2xbins,s2ybins,means_s2_xy_width_full,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. S2 Width vs. XY (end-start).'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'S2 width Samples. Normalized to z=0 \mus ','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);
saveas(s2_width_xy_end_fig,strcat('LUX_corrections/',file_id_cp,'/S2_Width_xy_end',file_id_cp),'fig');
saveas(s2_width_xy_end_fig,strcat('LUX_corrections/',file_id_cp,'/S2_Width_xy_end',file_id_cp),'jpg');
    

%% Calculate S2 Pulse Width Z Dependence (AFT_t1 - AFT_t0)x2 
%Redefine s2_width as aft_t1 - aft_t0 x2
s2_width       = (d.aft_t1_samples(s2_single_cut)-d.aft_t0_samples(s2_single_cut))*2;

%Make cut to clean up data and fit for the Z dependence
cleanup_cut    = inrange(s1_phe_both_z,[100,s1area_bound_max]) & inrange(s2_phe_both_z,[1000,s2area_bound_max]) & inrange(s2radius,[rcut_min,rcut_max]) ...
                 & inrange(s2_width,[100 1000]) & inrange(s1_width,[30 400]);

[mean_index_s2_z_width_aft1, means_s2_z_width_aft1, means_error_s2_z_width_aft1, P_s2_width_aft1, sigma_P_s2_width_aft1, dT_fit] = ...
                 MeasureS2PulseWidthZDep(s2_width, drift_time, det_edge, cleanup_cut); %#ok<ASGLU>
    
%apply Z dep S2 width correction
norm_s2_width = RK_1DCubicInterp_LinearExtrap(mean_index_s2_z_width_aft1,means_s2_z_width_aft1,4); %normalize to 0 us (the top)
s2_width_z    = s2_width.*norm_s2_width./RK_1DCubicInterp_LinearExtrap(mean_index_s2_z_width_aft1,means_s2_z_width_aft1,drift_time);
 
%Make plot of Z dependence
 
s2_z_width_aft1_fig = figure;
errorbar( mean_index_s2_z_width_aft1(2:end),means_s2_z_width_aft1(2:end),means_error_s2_z_width_aft1(2:end),'.k');  %first bin is anode, skip it
hold on    
[yplot_s2_width] = polyval(P_s2_width_aft1,dT_fit);
plot(dT_fit,yplot_s2_width,'-r','LineWidth',1);   
xlabel('Depth (\mus)', 'FontSize',16); ylabel('S2 Width (AFT1-AFT0)x2 (Samples)', 'FontSize', 16);
title(strcat(file_id_cp,'.Kr83m. S2 width vs. Z'), 'FontSize', 16,'Interpreter','none');
legend('S2 Width Mean', strcat('y=', num2str(P_s2_width_aft1(1),'%10.2e'), '*Z^3 + ', num2str(P_s2_width_aft1(2),'%10.2e')...
    , '*Z^2 + ', num2str(P_s2_width_aft1(3),'%10.2e'), '*Z + ',num2str(P_s2_width_aft1(4),4)),'location','northwest')
myfigview(16);
xlim([0 det_edge]);
saveas(s2_z_width_aft1_fig,strcat('LUX_corrections/',file_id_cp,'/S2_z_width_aft1_',file_id_cp),'jpg');
saveas(s2_z_width_aft1_fig,strcat('LUX_corrections/',file_id_cp,'/S2_z_width_aft1_',file_id_cp),'fig');
    

%% Calculate XY dependence of S2 Pulse Width (AFT_t1 - AFT_t0)x2 

%Make cut to clean up data and fit for the Z dependence   
cleanup_cut = inrange(drift_time,[zcut_min,zcut_max])  & inrange(s2_width_z,[50 700])& inrange(s2radius,[rcut_min,rcut_max]) & ...
              inrange(s2_phe_both_z,[1000,s2area_bound_max]) & inrange(s2radius,[0 24]) ;

[s2_xbins, s2_ybins, means_s2_xy_width_aft1, sigma_s2_xy_width_aft1 ] = MeasureS2PulseWidthXYDep(s2x, s2y, s2_width_z, S2_xbin_min, S2_xbin_max, S2_ybin_min, S2_ybin_max, S2xybinsize, cleanup_cut); %#ok<ASGLU>

%Make plot of XY dependence
s2_width_xy_aft_fig = figure;
color_range_max=max(max(means_s2_xy_width_aft1));
color_range_min=min(min(means_s2_xy_width_aft1(means_s2_xy_width_aft1>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
contourf(s2xbins,s2ybins,means_s2_xy_width_aft1,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. S2 Width vs. XY (AFT1-AFT0)x2.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'S2 width Samples. Normalized to z=0 \mus ','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);
saveas(s2_width_xy_aft_fig,strcat('LUX_corrections/',file_id_cp,'/S2_Width_xy_aft1',file_id_cp),'fig');
saveas(s2_width_xy_aft_fig,strcat('LUX_corrections/',file_id_cp,'/S2_Width_xy_aft1',file_id_cp),'jpg');


%% Save the data locally
   
save(strcat('LUX_corrections/',file_id_cp,'/',file_id_cp,'_IQs'),'e_lifetime','sigma_e_lifetime','Kr_events','P_s1_both','P_s1_bottom','sigma_P_s1_both','sigma_P_s1_bottom'...
    , 's1_both_means','s1_both_means_sigma','s1_both_means_bin','s1_bottom_means','s1_bottom_means_sigma','s1_bottom_means_bin'...
    , 's2_bottom_means','s2_bottom_means_sigma','s2_bottom_means_bin', 's2xbins','s2ybins','norm_S2_all','norm_S2_bot'...
    , 'mean_S2_both_xy','Sigma_S2_both','mean_S2_bottom_xy','Sigma_S2_bottom','Count_S2_bottom', 's1xbins','s1ybins'...
    ,'norm_S1_all','norm_S1_bot', 'mean_S1_both_xy','Sigma_S1_both','mean_S1_bottom_xy','Sigma_S1_bottom','Count_S1_both'...
    ,'sigma_norm_S1_all','sigma_norm_S1_bot','sigma_norm_S2_all','sigma_norm_S2_bot','file_id','file_id_cp','cp','s2_z0_both'...
    ,'sigma_s2_z0_both','e_lifetime_both','sigma_e_lifetime_both','s2_z0_bottom','sigma_s2_z0_bottom','s2_both_means','s2_both_means_bin'...
    ,'s2_both_means_sigma','mean_index_s2_z_width_full','means_s2_z_width_full','means_error_s2_z_width_full','mean_index_s2_z_width_aft1','means_s2_z_width_aft1'...
    ,'means_error_s2_z_width_aft1','means_s2_xy_width_full','sigma_s2_xy_width_full','means_s2_xy_width_aft1','sigma_s2_xy_width_aft1'...
    ,'P_s1_top','sigma_P_s1_top','norm_S1_top', 'sigma_norm_S1_top','P_s2_width_full','sigma_P_s2_width_full','P_s2_width_aft1','sigma_P_s2_width_aft1','det_edge'...
    ,'SE_mu_both','SE_sig_mu_both','SE_sig_both','SE_sig_sig_both','skew_both','SE_mu_bot','SE_sig_mu_bot','SE_sig_bot','SE_sig_sig_bot','skew_bot','SE_mu_top'...
    ,'SE_sig_mu_top','SE_sig_top','SE_sig_sig_top','skew_top','SE_both_mu_center','error_SE_both_mu_center','SE_both_kurt_center','SE_both_skew_center','SE_both_sigma_center'...
    ,'SE_bot_mu_center','error_SE_bot_mu_center','SE_bot_kurt_center','SE_bot_skew_center','SE_bot_sigma_center','SE_top_mu_center','error_SE_top_mu_center'...
    ,'SE_top_kurt_center','SE_top_skew_center','SE_top_sigma_center','se_xbins','se_ybins'...
    ,'mean_SE_both_xy','skew_SE_both_xy','sigma_SE_both_xy','kurt_SE_both_xy','mean_err_SE_both_xy','sigma_err_SE_both_xy','norm_SE_both_mu','error_SE_both_mu_norm'...
    ,'mean_SE_bot_xy','skew_SE_bot_xy','sigma_SE_bot_xy','kurt_SE_bot_xy','mean_err_SE_bot_xy','sigma_err_SE_bot_xy','norm_SE_bot_mu','error_SE_bot_mu_norm'...
    ,'mean_SE_top_xy','skew_SE_top_xy','sigma_SE_top_xy','kurt_SE_top_xy','mean_err_SE_top_xy','sigma_err_SE_top_xy','norm_SE_top_mu','error_SE_top_mu_norm'...
    ,'s1a_z_means','s1a_z_means_err','s1b_z_means','s1b_z_means_err','s1ab_bincenters','s1a_P','s1a_S','s1b_P','s1b_S','s1ab_P','s1ab_S'...
    ,'s1a_xy_mean','s1a_xy_mean_err','s1b_xy_mean','s1b_xy_mean_err','s1ab_xy_mean','s1ab_xy_mean_err','s1ab_xbins','s1ab_ybins'...
    ,'x_center','y_center','z_center'...
    ,'s2_bottom_norm_z_index','s2_bottom_norm_z_means','s2_at_top_bottom','s2_at_bottom_bottom','s2_bottom_norm_z','pseudo_lifetime_bottom'...
    ,'s2_both_norm_z_index','s2_both_norm_z_means','s2_at_top_both','s2_at_bottom_both','s2_both_norm_z','pseudo_lifetime_both'...
    ,'s1_phe_both_z_gauss_fit_mean','s1_phe_both_z_gauss_fit_sigma','s1_phe_both_z_gauss_fit_mean_error','s1_phe_both_z_gauss_fit_sigma_error'...
    ,'s2_phe_both_z_gauss_fit_mean','s2_phe_both_z_gauss_fit_sigma','s2_phe_both_z_gauss_fit_mean_error','s2_phe_both_z_gauss_fit_sigma_error'...
    ,'s1_phe_both_xyz_gauss_fit_mean','s1_phe_both_xyz_gauss_fit_sigma','s1_phe_both_xyz_gauss_fit_mean_error','s1_phe_both_xyz_gauss_fit_sigma_error'...
    ,'s2_phe_both_xyz_gauss_fit_mean','s2_phe_both_xyz_gauss_fit_sigma','s2_phe_both_xyz_gauss_fit_mean_error','s2_phe_both_xyz_gauss_fit_sigma_error'...
    ,'s1ab_r2bins','s1ab_zbins','s1a_r2z_mean','s1a_r2z_mean_err','s1b_r2z_mean','s1b_r2z_mean_err','s1ab_r2z_mean','s1ab_r2z_mean_err'...
    ,'field_map_2015','field_map_2015_staterr','field_map_2014','field_map_2014_staterr','mean_field_map','field_map_syserr','mean_field_map_totalerr');
    
 if length(s1ab_z)>100000;
    save(strcat('LUX_corrections/',file_id_cp,'/',file_id_cp,'_S1ab3DIQs'),'s1a_xyz_mean','s1a_xyz_mean_err','s1b_xyz_mean','s1b_xyz_mean_err','s1ab_xyz_mean','s1ab_xyz_mean_err'...
        ,'s1ab_xyz_xbins','s1ab_xyz_ybins','s1ab_xyz_zbins','s1as1b_at_center');
 else
    save(strcat('LUX_corrections/',file_id_cp,'/',file_id_cp,'_S1ab1D2DIQs'),'s1as1b_z_at_center','s1as1b_xy_at_center');     
 end

    
%% Submitting the results to the LUG

user    = user_name;
source  = 'Kr-83';
version = dp_version;
clear lug_val Image_paths xmlinput

if submit == 1;
%% S2 Z Dependence
    if calibrations(1)==1;
        lug_val(1)  = {file_id};
        lug_val(2)  = {cp}; %Proccessing version number
        lug_val(3)  = {gs}; %Proccessing gs number   
        lug_val(4)  = {user};
        lug_val(5)  = {algorithm_name};
        lug_val(6)  = {version};    
        lug_val(7)  = {e_lifetime};         % Electron lifetime (in us)
        lug_val(8)  = {sigma_e_lifetime};             % Electron lifetime - 1 sigma plus (in us)        
        lug_val(9)  = {s2_z0_bottom};        %S2_bottom intercept
        lug_val(10) = {sigma_s2_z0_bottom}; %S2_bottom intercept 1 sigma
        lug_val(11) = {e_lifetime_both};   %electron lifetime from Both PMT
        lug_val(12) = {sigma_e_lifetime_both};       %1 sigma electron lifetime from Both PMT
        lug_val(13) = {s2_z0_both};         %S2_both intercept
        lug_val(14) = {sigma_s2_z0_both};   %S2_both intercept 1 sigma  
        lug_val(19) = {Kr_events};          %Number of Kr events
        lug_val(20) = {s2_bottom_means_bin};% the bin center
        lug_val(21) = {s2_bottom_means};    % the mean of the gaussian fit
        lug_val(22) = {s2_bottom_means_sigma};% one sigma of the mean 
        lug_val(23) = {s2_both_means_bin};  % the bin center
        lug_val(24) = {s2_both_means};      % the mean of the gaussian fit
        lug_val(25) = {s2_both_means_sigma};% one sigma of the mean   
        lug_val(29) = {det_edge};
        lug_val(30) = {s2_both_center_z}; %Not actually center... it's the value at the top
        lug_val(31) = {s2_both_center_error_z};
        lug_val(32) = {s2_bot_center_z};
        lug_val(33) = {s2_bot_center_error_z};
        lug_val(34) = {s2_bottom_norm_z_index};
        lug_val(35) = {s2_bottom_norm_z_means};
        lug_val(36) = {s2_bottom_norm_z};
        lug_val(37) = {s2_at_top_bottom};
        lug_val(38) = {s2_at_bottom_bottom};
        lug_val(39) = {pseudo_lifetime_bottom};
        lug_val(40) = {s2_both_norm_z_index};
        lug_val(41) = {s2_both_norm_z_means};
        lug_val(42) = {s2_both_norm_z};
        lug_val(43) = {s2_at_top_both};
        lug_val(44) = {s2_at_bottom_both};
        lug_val(45) = {pseudo_lifetime_both};  
        lug_val(46) = {s2_phe_both_z_gauss_fit_mean};
        lug_val(47) = {s2_phe_both_z_gauss_fit_sigma};
        lug_val(48) = {s2_phe_both_z_gauss_fit_mean_error};
        lug_val(49) = {s2_phe_both_z_gauss_fit_sigma_error};
        lug_val(50) = {s2_phe_both_raw_gauss_fit_mean};
        lug_val(51) = {s2_phe_both_raw_gauss_fit_sigma};
        lug_val(52) = {s2_phe_both_raw_gauss_fit_mean_error};
        lug_val(53) = {s2_phe_both_raw_gauss_fit_sigma_error};

        Image_paths = {strcat('LUX_corrections/',file_id_cp,'/S2_bottom_lifetime_',file_id_cp,'.jpg'),strcat('LUX_corrections/',file_id_cp,'/S2_bottom_hist_EL_',file_id_cp,'.jpg'),...
                       strcat('LUX_corrections/',file_id_cp,'/S2_both_lifetime_',file_id_cp,'.jpg'),strcat('LUX_corrections/',file_id_cp,'/S2_both_hist_EL_',file_id_cp,'.jpg')}; 

        xmlinput.iq.global.filename_prefix                               = lug_val(1);
        xmlinput.iq.global.cp_number                                     = lug_val(2);
        xmlinput.iq.global.gs_number                                     = lug_val(3);
        xmlinput.iq.global.computed_by                                   = lug_val(4);
        xmlinput.iq.global.computed_date                                 = datestr(now,'yyyymmddTHHMM');
        xmlinput.iq.global.source                                        = source;
        xmlinput.iq.global.source_id                                     = 'unknown';
        xmlinput.iq.global.notes                                         = ' ';
        xmlinput.iq.global.algorithm_name                                = lug_val(5);
        xmlinput.iq.global.algorithm_version                             = lug_val(6);   
        xmlinput.iq.correction.fit.electron_attenuation_us               = lug_val(7); 
        xmlinput.iq.correction.fit.sigma_electron_attenuation_us         = lug_val(8);        
        xmlinput.iq.correction.fit.s2_bottom_z0_phe                      = lug_val(9);
        xmlinput.iq.correction.fit.sigma_s2_bottom_z0_phe                = lug_val(10);
        xmlinput.iq.correction.fit.electron_attenuation_us_both          = lug_val(11); 
        xmlinput.iq.correction.fit.sigma_electron_attenuation_us_both    = lug_val(12);        
        xmlinput.iq.correction.fit.s2_both_z0_phe                        = lug_val(13);
        xmlinput.iq.correction.fit.sigma_s2_both_z0_phe                  = lug_val(14);       
        xmlinput.iq.correction.fit.number_of_kr_events                   = lug_val(19);
        xmlinput.iq.correction.fit.bin_means_s2_bottom                   = lug_val(20);
        xmlinput.iq.correction.fit.means_s2_bottom                       = lug_val(21);
        xmlinput.iq.correction.fit.sigma_means_s2_bottom                 = lug_val(22); 
        xmlinput.iq.correction.fit.bin_means_s2_both                     = lug_val(23);
        xmlinput.iq.correction.fit.means_s2_both                         = lug_val(24);
        xmlinput.iq.correction.fit.sigma_means_s2_both                   = lug_val(25);
        xmlinput.iq.correction.fit.detector_edge                         = lug_val(29);
        xmlinput.iq.correction.fit.s2_both_mean_at_top_z                 = lug_val(30);
        xmlinput.iq.correction.fit.s2_both_mean_at_top_err_z             = lug_val(31);
        xmlinput.iq.correction.fit.s2_bot_mean_at_top_z                  = lug_val(32);
        xmlinput.iq.correction.fit.s2_bot_mean_at_top_err_z              = lug_val(33);
        xmlinput.iq.correction.fit.s2_bottom_norm_z_index                = lug_val(34);
        xmlinput.iq.correction.fit.s2_bottom_norm_z_means                = lug_val(35);
        xmlinput.iq.correction.fit.s2_bottom_norm_z                      = lug_val(36);
        xmlinput.iq.correction.fit.s2_at_top_bottom                      = lug_val(37);
        xmlinput.iq.correction.fit.s2_at_bottom_bottom                   = lug_val(38);
        xmlinput.iq.correction.fit.pseudo_lifetime_bottom                = lug_val(39);
        xmlinput.iq.correction.fit.s2_both_norm_z_index                  = lug_val(40);
        xmlinput.iq.correction.fit.s2_both_norm_z_means                  = lug_val(41);
        xmlinput.iq.correction.fit.s2_both_norm_z                        = lug_val(42);
        xmlinput.iq.correction.fit.s2_at_top_both                        = lug_val(43);
        xmlinput.iq.correction.fit.s2_at_bottom_both                     = lug_val(44);
        xmlinput.iq.correction.fit.pseudo_lifetime_both                  = lug_val(45);
        xmlinput.iq.correction.fit.s2_phe_both_z_gauss_fit_mean          = lug_val(46);
        xmlinput.iq.correction.fit.s2_phe_both_z_gauss_fit_sigma         = lug_val(47);
        xmlinput.iq.correction.fit.s2_phe_both_z_gauss_fit_mean_error    = lug_val(48);
        xmlinput.iq.correction.fit.s2_phe_both_z_gauss_fit_sigma_error   = lug_val(49);    
        xmlinput.iq.correction.fit.s2_phe_both_raw_gauss_fit_mean        = lug_val(50);
        xmlinput.iq.correction.fit.s2_phe_both_raw_gauss_fit_sigma       = lug_val(51);
        xmlinput.iq.correction.fit.s2_phe_both_raw_gauss_fit_mean_error  = lug_val(52);
        xmlinput.iq.correction.fit.s2_phe_both_raw_gauss_fit_sigma_error = lug_val(53);    

        LUXSubmitIQ(user,'electron_lifetime',algorithm_name,version,file_id,cp,xmlinput,'Kr83 Calibration','',Image_paths);

    end

    clear lug_val Image_paths xmlinput
   
%% S1 Z Dependence
    if calibrations(2)==1;

        lug_val(1)  = {file_id};
        lug_val(2)  = {cp}; %Proccessing version number
        lug_val(3)  = {gs};
        lug_val(4)  = {user};
        lug_val(5)  = {algorithm_name};
        lug_val(6)  = {version};
        lug_val(7)  = {P_s1_both};         % Poly paramaters for S1 vs. dT
        lug_val(8)  = {sigma_P_s1_both};   % Sigma Poly paramaters for S1 vs. dT
        lug_val(9)  = {P_s1_bottom};       % Poly paramaters for S1 vs. dT
        lug_val(10) = {sigma_P_s1_bottom}; % Sigma Poly paramaters for S1 vs. dT
        lug_val(11) = {P_s1_top};       % Poly paramaters for S1 vs. dT
        lug_val(12) = {sigma_P_s1_top}; % Sigma Poly paramaters for S1 vs. dT
        lug_val(13) = {Kr_events}; 
        lug_val(14) = {s1_both_means_bin};% the bin center
        lug_val(15) = {s1_both_means};  % the mean of the gaussian fit
        lug_val(16) = {s1_both_means_sigma};% one sigma
        lug_val(17) = {s1_bottom_means_bin};% the bin center
        lug_val(18) = {s1_bottom_means};  % the mean of the gaussian fit
        lug_val(19) = {s1_bottom_means_sigma};% one sigma    
        lug_val(20) = {center_s1_both};
        lug_val(21) = {center_s1_both_err};
        lug_val(22) = {center_s1_bottom};
        lug_val(23) = {center_s1_bottom_err};    
        lug_val(24) = {center_s1_top};
        lug_val(25) = {center_s1_top_err};
        lug_val(26) = {s1_phe_both_z_gauss_fit_mean};
        lug_val(27) = {s1_phe_both_z_gauss_fit_sigma};
        lug_val(28) = {s1_phe_both_z_gauss_fit_mean_error};
        lug_val(29) = {s1_phe_both_z_gauss_fit_sigma_error};
        lug_val(30) = {s1_phe_both_raw_gauss_fit_mean};
        lug_val(31) = {s1_phe_both_raw_gauss_fit_sigma};
        lug_val(32) = {s1_phe_both_raw_gauss_fit_mean_error};
        lug_val(33) = {s1_phe_both_raw_gauss_fit_sigma_error};    

        Image_paths = {strcat('LUX_corrections/',file_id_cp,'/S1_both_LC_',file_id_cp,'.jpg'),strcat('LUX_corrections/',file_id_cp,'/S1_both_hist_',file_id_cp,'.jpg'),...
                       strcat('LUX_corrections/',file_id_cp,'/S1_bottom_LC_',file_id_cp,'.jpg'),strcat('LUX_corrections/',file_id_cp,'/S1_bottom_hist_',file_id_cp,'.jpg') };  

        xmlinput.iq.global.filename_prefix                               = lug_val(1);
        xmlinput.iq.global.cp_number                                     = lug_val(2);
        xmlinput.iq.global.gs_number                                     = lug_val(3);
        xmlinput.iq.global.computed_by                                   = lug_val(4);
        xmlinput.iq.global.computed_date                                 = datestr(now,'yyyymmddTHHMM');
        xmlinput.iq.global.source                                        = source;
        xmlinput.iq.global.source_id                                     = 'unknown';
        xmlinput.iq.global.notes                                         = ' ';
        xmlinput.iq.global.algorithm_name                                = lug_val(5);
        xmlinput.iq.global.algorithm_version                             = lug_val(6);   
        xmlinput.iq.correction.fit.s1_both_zdep_quad_fit                 = lug_val(7);
        xmlinput.iq.correction.fit.s1_both_sigma_zdep_quad_fit           = lug_val(8);
        xmlinput.iq.correction.fit.s1_bottom_zdep_quad_fit               = lug_val(9);
        xmlinput.iq.correction.fit.s1_bottom_sigma_zdep_quad_fit         = lug_val(10);
        xmlinput.iq.correction.fit.s1_top_zdep_quad_fit                  = lug_val(11);
        xmlinput.iq.correction.fit.s1_top_sigma_zdep_quad_fit            = lug_val(12);
        xmlinput.iq.correction.fit.number_of_kr_events                   = lug_val(13);
        xmlinput.iq.correction.fit.bin_means_s1_both                     = lug_val(14);
        xmlinput.iq.correction.fit.means_s1_both                         = lug_val(15);
        xmlinput.iq.correction.fit.sigma_means_s1_both                   = lug_val(16);
        xmlinput.iq.correction.fit.bin_means_s1_bottom                   = lug_val(17);
        xmlinput.iq.correction.fit.means_s1_bottom                       = lug_val(18);
        xmlinput.iq.correction.fit.sigma_means_s1_bottom                 = lug_val(19);            
        xmlinput.iq.correction.fit.center_s1_both_z                      = lug_val(20);            
        xmlinput.iq.correction.fit.center_s1_both_err_z                  = lug_val(21);            
        xmlinput.iq.correction.fit.center_s1_bot_z                       = lug_val(22);            
        xmlinput.iq.correction.fit.center_s1_bot_err_z                   = lug_val(23);   
        xmlinput.iq.correction.fit.center_s1_top_z                       = lug_val(24);            
        xmlinput.iq.correction.fit.center_s1_top_err_z                   = lug_val(25);   
        xmlinput.iq.correction.fit.s1_phe_both_z_gauss_fit_mean          = lug_val(26);            
        xmlinput.iq.correction.fit.s1_phe_both_z_gauss_fit_sigma         = lug_val(27);   
        xmlinput.iq.correction.fit.s1_phe_both_z_gauss_fit_mean_error    = lug_val(28);            
        xmlinput.iq.correction.fit.s1_phe_both_z_gauss_fit_sigma_error   = lug_val(29);   
        xmlinput.iq.correction.fit.s1_phe_both_raw_gauss_fit_mean        = lug_val(30);            
        xmlinput.iq.correction.fit.s1_phe_both_raw_gauss_fit_sigma       = lug_val(31);   
        xmlinput.iq.correction.fit.s1_phe_both_raw_gauss_fit_mean_error  = lug_val(32);            
        xmlinput.iq.correction.fit.s1_phe_both_raw_gauss_fit_sigma_error = lug_val(33);      


        LUXSubmitIQ(user,'z_dep_s1_correction',algorithm_name,version,file_id,cp,xmlinput,'Kr83 Calibration','',Image_paths);

    end

    clear lug_val Image_paths xmlinput
  
%% S2 XY Correction
    if calibrations(3)==1;

        lug_val(1)  = {file_id};
        lug_val(2)  = {user};
        lug_val(3)  = {algorithm_name};
        lug_val(4)  = {version};    
        lug_val(5)  = {s2xbins};              
        lug_val(6)  = {s2ybins};           
        lug_val(7)  = {norm_S2_all};
        lug_val(8)  = {sigma_norm_S2_all};
        lug_val(9)  = {norm_S2_bot};
        lug_val(10) = {sigma_norm_S2_bot};
        lug_val(11) = {mean_S2_both_xy};
        lug_val(12) = {Sigma_S2_both};
        lug_val(13) = {mean_S2_bottom_xy};
        lug_val(14) = {Sigma_S2_bottom};
        lug_val(15) = {Count_S2_bottom};
        lug_val(16) = {cp}; %Proccessing version number
        lug_val(17) = {gs}; %Proccessing gs number
        lug_val(18) = {Kr_events};
        lug_val(20) = {s2_both_center_xyz}; %Not actually "center"... it's the value at the top
        lug_val(21) = {s2_both_center_xyz_error};
        lug_val(22) = {s2_bot_center_xyz};
        lug_val(23) = {s2_bot_center_xyz_error};   
        lug_val(24) = {s2_phe_both_xyz_gauss_fit_mean};
        lug_val(25) = {s2_phe_both_xyz_gauss_fit_sigma};
        lug_val(26) = {s2_phe_both_xyz_gauss_fit_mean_error};
        lug_val(27) = {s2_phe_both_xyz_gauss_fit_sigma_error};

        Image_paths = {strcat('LUX_corrections/',file_id_cp,'/S2_XY_both_',file_id_cp,'.jpg'), strcat('LUX_corrections/',file_id_cp,'/Sigma_S2_XY_both_',file_id_cp,'.jpg'),...
                       strcat('LUX_corrections/',file_id_cp,'/S2_XY_bottom_',file_id_cp,'.jpg'), strcat('LUX_corrections/',file_id_cp,'/Sigma_S2_XY_bottom_',file_id_cp,'.jpg') }; 

        xmlinput.iq.global.filename_prefix                               = lug_val(1);
        xmlinput.iq.global.cp_number                                     = lug_val(16);
        xmlinput.iq.global.gs_number                                     = lug_val(17);
        xmlinput.iq.global.computed_by                                   = lug_val(2);
        xmlinput.iq.global.computed_date                                 = datestr(now,'yyyymmddTHHMM');
        xmlinput.iq.global.source                                        = source;
        xmlinput.iq.global.source_id                                     = 'unknown';
        xmlinput.iq.global.notes                                         = ' ';
        xmlinput.iq.global.algorithm_name                                = lug_val(3);
        xmlinput.iq.global.algorithm_version                             = lug_val(4);
        xmlinput.iq.correction.fit.number_of_kr_events                   = lug_val(18);
        xmlinput.iq.correction.fit.x_bin_center                          = lug_val(5); 
        xmlinput.iq.correction.fit.y_bin_center                          = lug_val(6);
        xmlinput.iq.correction.fit.norm_s2_both                          = lug_val(7);
        xmlinput.iq.correction.fit.sigma_norm_s2_both                    = lug_val(8);
        xmlinput.iq.correction.fit.norm_s2_bottom                        = lug_val(9);
        xmlinput.iq.correction.fit.sigma_norm_s2_bottom                  = lug_val(10);
        xmlinput.iq.correction.fit.mean_s2_both                          = lug_val(11);
        xmlinput.iq.correction.fit.sigma_s2_both                         = lug_val(12);
        xmlinput.iq.correction.fit.mean_s2_bottom                        = lug_val(13);
        xmlinput.iq.correction.fit.sigma_s2_bottom                       = lug_val(14);
        xmlinput.iq.correction.fit.count_s2                              = lug_val(15);
        xmlinput.iq.correction.fit.s2_both_mean_at_top_xyz               = lug_val(20);
        xmlinput.iq.correction.fit.s2_both_mean_at_top_err_xyz           = lug_val(21);
        xmlinput.iq.correction.fit.s2_bot_mean_at_top_xyz                = lug_val(22);
        xmlinput.iq.correction.fit.s2_bot_mean_at_top_err_xyz            = lug_val(23);
        xmlinput.iq.correction.fit.s2_phe_both_xyz_gauss_fit_mean        = lug_val(24);
        xmlinput.iq.correction.fit.s2_phe_both_xyz_gauss_fit_sigma       = lug_val(25);
        xmlinput.iq.correction.fit.s2_phe_both_xyz_gauss_fit_mean_error  = lug_val(26);
        xmlinput.iq.correction.fit.s2_phe_both_xyz_gauss_fit_sigma_error = lug_val(27);

        LUXSubmitIQ(user,'s2_xy_correction',algorithm_name,version,file_id,cp,xmlinput,'Kr83 Calibration','',Image_paths);
    end

    clear lug_val Image_paths xmlinput
    
%% S1 XY Correction
    if calibrations(4)==1;
        lug_val(1)  = {file_id};
        lug_val(2)  = {cp}; %Proccessing version number
        lug_val(3)  = {gs};
        lug_val(4)  = {user};
        lug_val(5)  = {algorithm_name};
        lug_val(6)  = {version}; 
        lug_val(7)  = {Kr_events};
        lug_val(8)  = {s1xbins};                
        lug_val(9)  = {s1ybins};           
        lug_val(10) = {norm_S1_all};
        lug_val(11) = {sigma_norm_S1_all};
        lug_val(12) = {norm_S1_bot};
        lug_val(13) = {sigma_norm_S1_bot};
        lug_val(14) = {norm_S1_top};
        lug_val(15) = {sigma_norm_S1_top};
        lug_val(16) = {mean_S1_both_xy};
        lug_val(17) = {Sigma_S1_both};
        lug_val(18) = {mean_S1_bottom_xy};
        lug_val(19) = {Sigma_S1_bottom};
        lug_val(20) = {Count_S1_both};
        lug_val(21) = {s1_both_center_xyz};
        lug_val(22) = {s1_both_center_xyz_error};
        lug_val(23) = {s1_bot_center_xyz};
        lug_val(24) = {s1_bot_center_xyz_error};
        lug_val(25) = {s1_top_center_xyz};
        lug_val(26) = {s1_top_center_xyz_error};
        lug_val(27) = {s1_phe_both_xyz_gauss_fit_mean};
        lug_val(28) = {s1_phe_both_xyz_gauss_fit_sigma};
        lug_val(29) = {s1_phe_both_xyz_gauss_fit_mean_error};
        lug_val(30) = {s1_phe_both_xyz_gauss_fit_sigma_error};

        Image_paths = {strcat('LUX_corrections/',file_id_cp,'/S1_XY_both_',file_id_cp,'.jpg'), strcat('LUX_corrections/',file_id_cp,'/Sigma_S1_XY_both_',file_id_cp,'.jpg')...
                      ,strcat('LUX_corrections/',file_id_cp,'/S1_XY_bottom_',file_id_cp,'.jpg'), strcat('LUX_corrections/',file_id_cp,'/Sigma_S1_XY_bottom_',file_id_cp,'.jpg')...
                      ,strcat('LUX_corrections/',file_id_cp,'/S1_Count_XY_both_',file_id_cp,'.jpg') }; 


        xmlinput.iq.global.filename_prefix                               = lug_val(1);
        xmlinput.iq.global.cp_number                                     = lug_val(2);
        xmlinput.iq.global.gs_number                                     = lug_val(3);
        xmlinput.iq.global.computed_by                                   = lug_val(4);
        xmlinput.iq.global.computed_date                                 = datestr(now,'yyyymmddTHHMM');
        xmlinput.iq.global.source                                        = source;
        xmlinput.iq.global.source_id                                     = 'unknown';
        xmlinput.iq.global.notes                                         = ' ';
        xmlinput.iq.global.algorithm_name                                = lug_val(5);
        xmlinput.iq.global.algorithm_version                             = lug_val(6);
        xmlinput.iq.correction.fit.number_of_kr_events                   = lug_val(7);
        xmlinput.iq.correction.fit.x_bin_center                          = lug_val(8); 
        xmlinput.iq.correction.fit.y_bin_center                          = lug_val(9);
        xmlinput.iq.correction.fit.norm_s1_both                          = lug_val(10);
        xmlinput.iq.correction.fit.sigma_norm_s1_both                    = lug_val(11);
        xmlinput.iq.correction.fit.norm_s1_bottom                        = lug_val(12);
        xmlinput.iq.correction.fit.sigma_norm_s1_bottom                  = lug_val(13);
        xmlinput.iq.correction.fit.norm_s1_top                           = lug_val(14);
        xmlinput.iq.correction.fit.sigma_norm_s1_top                     = lug_val(15);
        xmlinput.iq.correction.fit.mean_s1_both                          = lug_val(16);
        xmlinput.iq.correction.fit.sigma_s1_both                         = lug_val(17);
        xmlinput.iq.correction.fit.mean_s1_bottom                        = lug_val(18);
        xmlinput.iq.correction.fit.sigma_s1_bottom                       = lug_val(19);
        xmlinput.iq.correction.fit.count_s1                              = lug_val(20);
        xmlinput.iq.correction.fit.center_s1_both_z                      = lug_val(21);
        xmlinput.iq.correction.fit.center_s1_both_err_z                  = lug_val(22);
        xmlinput.iq.correction.fit.center_s1_bot_z                       = lug_val(23);
        xmlinput.iq.correction.fit.center_s1_bot_err_z                   = lug_val(24);
        xmlinput.iq.correction.fit.center_s1_top_z                       = lug_val(25);
        xmlinput.iq.correction.fit.center_s1_top_err_z                   = lug_val(26);
        xmlinput.iq.correction.fit.s1_phe_both_xyz_gauss_fit_mean        = lug_val(27);
        xmlinput.iq.correction.fit.s1_phe_both_xyz_gauss_fit_sigma       = lug_val(28);
        xmlinput.iq.correction.fit.s1_phe_both_xyz_gauss_fit_mean_error  = lug_val(29);
        xmlinput.iq.correction.fit.s1_phe_both_xyz_gauss_fit_sigma_error = lug_val(30);

        LUXSubmitIQ(user,'s1_xy_correction',algorithm_name,version,file_id,cp,xmlinput,'Kr83 Calibration','',Image_paths);
    end   

    clear lug_val Image_paths xmlinput
        
%% S2 Width
    if calibrations(6)==1;

        lug_val(1)  = {file_id};
        lug_val(2)  = {user};
        lug_val(3)  = {algorithm_name};
        lug_val(4)  = {version};   
        lug_val(5)  = {cp}; %Proccessing version number
        lug_val(6)  = {gs}; %Proccessing gs number
        lug_val(7)  = {Kr_events};   
        lug_val(8)  = {mean_index_s2_z_width_full};
        lug_val(9)  = {means_s2_z_width_full};
        lug_val(10) = {means_error_s2_z_width_full};
        lug_val(11) = {mean_index_s2_z_width_aft1};
        lug_val(12) = {means_s2_z_width_aft1};
        lug_val(13) = {means_error_s2_z_width_aft1};
        lug_val(14) = {s2xbins};              
        lug_val(15) = {s2ybins};           
        lug_val(16) = {means_s2_xy_width_full};
        lug_val(17) = {sigma_s2_xy_width_full};
        lug_val(18) = {means_s2_xy_width_aft1};
        lug_val(19) = {sigma_s2_xy_width_aft1};

        Image_paths = {strcat('LUX_corrections/',file_id_cp,'/S2_z_width_end_',file_id_cp,'.jpg'), strcat('LUX_corrections/',file_id_cp,'/S2_Width_xy_end',file_id_cp,'.jpg')...
                      ,strcat('LUX_corrections/',file_id_cp,'/S2_z_width_aft1_',file_id_cp,'.jpg'), strcat('LUX_corrections/',file_id_cp,'/S2_Width_xy_aft1',file_id_cp,'.jpg')... 
                      ,strcat('LUX_corrections/',file_id_cp,'/RadialField_',file_id_cp,'.jpg') }; 

        xmlinput.iq.global.filename_prefix                     = lug_val(1);          
        xmlinput.iq.global.computed_by                         = lug_val(2);
        xmlinput.iq.global.computed_date                       = datestr(now,'yyyymmddTHHMM');
        xmlinput.iq.global.source                              = source;
        xmlinput.iq.global.source_id                           = 'unknown';
        xmlinput.iq.global.notes                               = ' ';
        xmlinput.iq.global.algorithm_name                      = lug_val(3);
        xmlinput.iq.global.algorithm_version                   = lug_val(4);
        xmlinput.iq.global.cp_number                           = lug_val(5);
        xmlinput.iq.global.gs_number                           = lug_val(6);
        xmlinput.iq.correction.fit.number_of_kr_events         = lug_val(7);
        xmlinput.iq.correction.fit.bin_means_s2_z_width_full   = lug_val(8);
        xmlinput.iq.correction.fit.means_s2_z_width_full       = lug_val(9);
        xmlinput.iq.correction.fit.sigma_means_s2_z_width_full = lug_val(10);
        xmlinput.iq.correction.fit.bin_means_s2_z_width_aft1   = lug_val(11);
        xmlinput.iq.correction.fit.means_s2_z_width_aft1       = lug_val(12);
        xmlinput.iq.correction.fit.sigma_means_s2_z_width_aft1 = lug_val(13);                        
        xmlinput.iq.correction.fit.x_bin_center                = lug_val(14); 
        xmlinput.iq.correction.fit.y_bin_center                = lug_val(15);            
        xmlinput.iq.correction.fit.means_s2_xy_width_full      = lug_val(16);
        xmlinput.iq.correction.fit.sigma_s2_xy_width_full      = lug_val(17);
        xmlinput.iq.correction.fit.means_s2_xy_width_aft1      = lug_val(18);
        xmlinput.iq.correction.fit.sigma_s2_xy_width_aft1      = lug_val(19);

        LUXSubmitIQ(user,'s2_width_kr',algorithm_name,version,file_id,cp,xmlinput,'Kr83 Calibration','',Image_paths);

    end
    clear lug_val Image_paths xmlinput    

%% Single electron size 
    if calibrations(6)==1;

        lug_val(1)  = {file_id};
        lug_val(2)  = {user};
        lug_val(3)  = {algorithm_name};
        lug_val(4)  = {version};   
        lug_val(5)  = {cp}; %Proccessing version number
        lug_val(6)  = {gs}; %Proccessing gs number
        lug_val(7)  = {Kr_events};   
        lug_val(8)  = {SE_mu_both};
        lug_val(9)  = {SE_sig_mu_both};
        lug_val(10) = {SE_sig_both};
        lug_val(11) = {SE_sig_sig_both};
        lug_val(12) = {skew_both};
        lug_val(13) = {SE_mu_bot};
        lug_val(14) = {SE_sig_mu_bot};
        lug_val(15) = {SE_sig_bot};
        lug_val(16) = {SE_sig_sig_bot};
        lug_val(17) = {skew_bot};
        lug_val(18) = {SE_mu_top};
        lug_val(19) = {SE_sig_mu_top};
        lug_val(20) = {SE_sig_top};
        lug_val(21) = {SE_sig_sig_top};
        lug_val(22) = {skew_top};
        lug_val(23) = {num_SE_events};

        Image_paths = {strcat('LUX_corrections/',file_id_cp,'/both_SE_',file_id_cp,'.jpg'), strcat('LUX_corrections/',file_id_cp,'/bottom_SE_',file_id_cp,'.jpg')...
                      ,strcat('LUX_corrections/',file_id_cp,'/top_SE_',file_id_cp,'.jpg') };

        xmlinput.iq.global.filename_prefix             = lug_val(1);          
        xmlinput.iq.global.computed_by                 = lug_val(2);
        xmlinput.iq.global.computed_date               = datestr(now,'yyyymmddTHHMM');
        xmlinput.iq.global.source                      = source;
        xmlinput.iq.global.source_id                   = 'unknown';
        xmlinput.iq.global.notes                       = ' ';
        xmlinput.iq.global.algorithm_name              = lug_val(3);
        xmlinput.iq.global.algorithm_version           = lug_val(4);
        xmlinput.iq.global.cp_number                   = lug_val(5);
        xmlinput.iq.global.gs_number                   = lug_val(6);
        xmlinput.iq.correction.fit.number_of_kr_events = lug_val(7);
        xmlinput.iq.correction.fit.e_mean_both         = lug_val(8);
        xmlinput.iq.correction.fit.e_mean_err_both     = lug_val(9);
        xmlinput.iq.correction.fit.e_sig_both          = lug_val(10);
        xmlinput.iq.correction.fit.e_sig_err_both      = lug_val(11);
        xmlinput.iq.correction.fit.e_fit_skew_both     = lug_val(12);
        xmlinput.iq.correction.fit.e_mean_bottom       = lug_val(13);
        xmlinput.iq.correction.fit.e_mean_err_bottom   = lug_val(14);
        xmlinput.iq.correction.fit.e_sig_bottom        = lug_val(15);
        xmlinput.iq.correction.fit.e_sig_err_bottom    = lug_val(16);
        xmlinput.iq.correction.fit.e_fit_skew_bottom   = lug_val(17);
        xmlinput.iq.correction.fit.e_mean_top          = lug_val(18);
        xmlinput.iq.correction.fit.e_mean_err_top      = lug_val(19);
        xmlinput.iq.correction.fit.e_sig_top           = lug_val(20);
        xmlinput.iq.correction.fit.e_sig_err_top       = lug_val(21);
        xmlinput.iq.correction.fit.e_fit_skew_top      = lug_val(22);
        xmlinput.iq.correction.fit.number_SE_events    = lug_val(23); 

        LUXSubmitIQ(user,'single_e_kr',algorithm_name,version,file_id,cp,xmlinput,'Kr83 Calibration','',Image_paths);
        clear lug_val Image_paths xmlinput    

        %% SE XY Maps
        lug_val(1)  = {file_id};
        lug_val(2)  = {user};
        lug_val(3)  = {algorithm_name};
        lug_val(4)  = {version};   
        lug_val(5)  = {cp}; %Proccessing version number
        lug_val(6)  = {gs}; %Proccessing gs number
        lug_val(7)  = {Kr_events}; 
        lug_val(8)  = {SE_both_mu_center};
        lug_val(9)  = {error_SE_both_mu_center};
        lug_val(10) = {SE_both_kurt_center};
        lug_val(11) = {SE_both_skew_center};
        lug_val(12) = {SE_both_sigma_center};
        lug_val(13) = {SE_bot_mu_center};
        lug_val(14) = {error_SE_bot_mu_center};
        lug_val(15) = {SE_bot_kurt_center};
        lug_val(16) = {SE_bot_skew_center};
        lug_val(17) = {SE_bot_sigma_center};
        lug_val(18) = {SE_top_mu_center};
        lug_val(19) = {error_SE_top_mu_center};
        lug_val(20) = {SE_top_kurt_center};
        lug_val(21) = {SE_top_skew_center};
        lug_val(22) = {SE_top_sigma_center};
        lug_val(23) = {se_xbins};
        lug_val(24) = {se_ybins};
        lug_val(25) = {mean_SE_both_xy};
        lug_val(26) = {skew_SE_both_xy};
        lug_val(27) = {sigma_SE_both_xy};
        lug_val(28) = {kurt_SE_both_xy};
        lug_val(29) = {mean_err_SE_both_xy};
        lug_val(30) = {sigma_err_SE_both_xy};
        lug_val(31) = {norm_SE_both_mu};
        lug_val(32) = {error_SE_both_mu_norm};
        lug_val(33) = {mean_SE_bot_xy};
        lug_val(34) = {skew_SE_bot_xy};
        lug_val(35) = {sigma_SE_bot_xy};
        lug_val(36) = {kurt_SE_bot_xy};
        lug_val(37) = {mean_err_SE_bot_xy};
        lug_val(38) = {sigma_err_SE_bot_xy};
        lug_val(39) = {norm_SE_bot_mu};
        lug_val(40) = {error_SE_bot_mu_norm};
        lug_val(41) = {mean_SE_top_xy};
        lug_val(42) = {skew_SE_top_xy};
        lug_val(43) = {sigma_SE_top_xy};
        lug_val(44) = {kurt_SE_top_xy};
        lug_val(45) = {mean_err_SE_top_xy};
        lug_val(46) = {sigma_err_SE_top_xy};
        lug_val(47) = {norm_SE_top_mu};
        lug_val(48) = {error_SE_top_mu_norm};

        Image_paths = {strcat('LUX_corrections/',file_id_cp,'/SE_XY_both_',file_id_cp,'.jpg')...
                      ,strcat('LUX_corrections/',file_id_cp,'/SE_XY_bot_',file_id_cp,'.jpg')...
                      ,strcat('LUX_corrections/',file_id_cp,'/SE_XY_top_',file_id_cp,'.jpg')};

        xmlinput.iq.global.filename_prefix                  = lug_val(1);          
        xmlinput.iq.global.computed_by                      = lug_val(2);
        xmlinput.iq.global.computed_date                    = datestr(now,'yyyymmddTHHMM');
        xmlinput.iq.global.source                           = source;
        xmlinput.iq.global.source_id                        = 'unknown';
        xmlinput.iq.global.notes                            = ' ';
        xmlinput.iq.global.algorithm_name                   = lug_val(3);
        xmlinput.iq.global.algorithm_version                = lug_val(4);
        xmlinput.iq.global.cp_number                        = lug_val(5);
        xmlinput.iq.global.gs_number                        = lug_val(6);
        xmlinput.iq.correction.fit.number_of_kr_events      = lug_val(7);
        xmlinput.iq.correction.fit.SE_both_mu_center        = lug_val(8);
        xmlinput.iq.correction.fit.SE_both_mu_center_error  = lug_val(9);
        xmlinput.iq.correction.fit.SE_both_kurt_center      = lug_val(10);
        xmlinput.iq.correction.fit.SE_both_skew_center      = lug_val(11);
        xmlinput.iq.correction.fit.SE_both_sigma_center     = lug_val(12);
        xmlinput.iq.correction.fit.SE_bot_mu_center         = lug_val(13);
        xmlinput.iq.correction.fit.SE_bot_mu_center_error   = lug_val(14);
        xmlinput.iq.correction.fit.SE_bot_kurt_center       = lug_val(15);
        xmlinput.iq.correction.fit.SE_bot_skew_center       = lug_val(16);
        xmlinput.iq.correction.fit.SE_bot_sigma_center      = lug_val(17);
        xmlinput.iq.correction.fit.SE_top_mu_center         = lug_val(18);
        xmlinput.iq.correction.fit.SE_top_mu_center_error   = lug_val(19);
        xmlinput.iq.correction.fit.SE_top_kurt_center       = lug_val(20);
        xmlinput.iq.correction.fit.SE_top_skew_center       = lug_val(21);
        xmlinput.iq.correction.fit.SE_top_sigma_center      = lug_val(22);
        xmlinput.iq.correction.fit.SE_xbin_centers          = lug_val(23);
        xmlinput.iq.correction.fit.SE_ybin_centers          = lug_val(24);
        xmlinput.iq.correction.fit.SE_both_mean_xy          = lug_val(25);
        xmlinput.iq.correction.fit.SE_both_skew_xy          = lug_val(26);
        xmlinput.iq.correction.fit.SE_both_sigma_xy         = lug_val(27);
        xmlinput.iq.correction.fit.SE_both_kurt_xy          = lug_val(28);
        xmlinput.iq.correction.fit.SE_both_mean_err_xy      = lug_val(29);
        xmlinput.iq.correction.fit.SE_both_sigma_err_xy     = lug_val(30);
        xmlinput.iq.correction.fit.SE_both_mean_xy_norm     = lug_val(31);
        xmlinput.iq.correction.fit.SE_both_mean_err_xy_norm = lug_val(32);
        xmlinput.iq.correction.fit.SE_bot_mean_xy           = lug_val(33);
        xmlinput.iq.correction.fit.SE_bot_skew_xy           = lug_val(34);
        xmlinput.iq.correction.fit.SE_bot_sigma_xy          = lug_val(35);
        xmlinput.iq.correction.fit.SE_bot_kurt_xy           = lug_val(36);
        xmlinput.iq.correction.fit.SE_bot_mean_err_xy       = lug_val(37);
        xmlinput.iq.correction.fit.SE_bot_sigma_err_xy      = lug_val(38);
        xmlinput.iq.correction.fit.SE_bot_mean_xy_norm      = lug_val(39);
        xmlinput.iq.correction.fit.SE_bot_mean_err_xy_norm  = lug_val(40);
        xmlinput.iq.correction.fit.SE_top_mean_xy           = lug_val(41);
        xmlinput.iq.correction.fit.SE_top_skew_xy           = lug_val(42);
        xmlinput.iq.correction.fit.SE_top_sigma_xy          = lug_val(43);
        xmlinput.iq.correction.fit.SE_top_kurt_xy           = lug_val(44);
        xmlinput.iq.correction.fit.SE_top_mean_err_xy       = lug_val(45);
        xmlinput.iq.correction.fit.SE_top_sigma_err_xy      = lug_val(46);
        xmlinput.iq.correction.fit.SE_top_mean_xy_norm      = lug_val(47);
        xmlinput.iq.correction.fit.SE_top_mean_err_xy_norm  = lug_val(48);                                    

        LUXSubmitIQ(user,'single_e_xy_kr',algorithm_name,version,file_id,cp,xmlinput,'Kr83 Calibration','',Image_paths);

    end
    clear lug_val Image_paths xmlinput   
       
%% Kr s1a/s1b Z maps

    %Added after submit == 1 was hard coded, so this doesn't check for that condition
    lug_val(1)  = {file_id};
    lug_val(2)  = {user};
    lug_val(3)  = {algorithm_name};
    lug_val(4)  = {version};   
    lug_val(5)  = {cp}; %Proccessing version number
    lug_val(6)  = {gs}; %Proccessing gs number
    lug_val(7)  = {Kr_events};   
    lug_val(8)  = {s1a_z_means};
    lug_val(9)  = {s1a_z_means_err};
    lug_val(10) = {s1b_z_means};
    lug_val(11) = {s1b_z_means_err};
    lug_val(12) = {s1ab_z_means};
    lug_val(13) = {s1ab_z_means_err};
    lug_val(14) = {s1ab_bincenters};
    lug_val(15) = {s1a_P};
    lug_val(16) = {s1a_S.R};    
    lug_val(17) = {s1a_S.df};             
    lug_val(18) = {s1a_S.normr};              
    lug_val(19) = {s1b_P};       
    lug_val(20) = {s1b_S.R};
    lug_val(21) = {s1b_S.df};
    lug_val(22) = {s1b_S.normr};   
    lug_val(23) = {s1ab_P};
    lug_val(24) = {s1ab_S.R};
    lug_val(25) = {s1ab_S.df};
    lug_val(26) = {s1ab_S.normr};

    Image_paths = {strcat('LUX_corrections/',file_id_cp,'/s1a_mean_z',file_id_cp,'.jpg'), strcat('LUX_corrections/',file_id_cp,'/s1b_mean_z',file_id_cp,'.jpg')...
                  ,strcat('LUX_corrections/',file_id_cp,'/s1ab_ratio_z',file_id_cp,'.jpg')}; 

    xmlinput.iq.global.filename_prefix             = lug_val(1);          
    xmlinput.iq.global.computed_by                 = lug_val(2);
    xmlinput.iq.global.computed_date               = datestr(now,'yyyymmddTHHMM');
    xmlinput.iq.global.source                      = source;
    xmlinput.iq.global.source_id                   = 'unknown';
    xmlinput.iq.global.notes                       = ' ';
    xmlinput.iq.global.algorithm_name              = lug_val(3);
    xmlinput.iq.global.algorithm_version           = lug_val(4);
    xmlinput.iq.global.cp_number                   = lug_val(5);
    xmlinput.iq.global.gs_number                   = lug_val(6);
    xmlinput.iq.correction.fit.number_of_kr_events = lug_val(7);
    xmlinput.iq.correction.fit.s1a_zmeans          = lug_val(8);
    xmlinput.iq.correction.fit.s1a_zmeans_err      = lug_val(9);
    xmlinput.iq.correction.fit.s1b_zmeans          = lug_val(10);
    xmlinput.iq.correction.fit.s1b_zmeans_err      = lug_val(11);
    xmlinput.iq.correction.fit.s1ab_zratios        = lug_val(12);
    xmlinput.iq.correction.fit.s1ab_zratios_err    = lug_val(13);
    xmlinput.iq.correction.fit.s1ab_zbincenters    = lug_val(14);
    xmlinput.iq.correction.fit.s1a_P               = lug_val(15);                        
    xmlinput.iq.correction.fit.s1a_S.R             = lug_val(16);
    xmlinput.iq.correction.fit.s1a_S.df            = lug_val(17); 
    xmlinput.iq.correction.fit.s1a_S.normr         = lug_val(18); 
    xmlinput.iq.correction.fit.s1b_P               = lug_val(19);            
    xmlinput.iq.correction.fit.s1b_S.R             = lug_val(20);
    xmlinput.iq.correction.fit.s1b_S.df            = lug_val(21); 
    xmlinput.iq.correction.fit.s1b_S.normr         = lug_val(22); 
    xmlinput.iq.correction.fit.s1ab_P              = lug_val(23);
    xmlinput.iq.correction.fit.s1ab_S.R            = lug_val(24);
    xmlinput.iq.correction.fit.s1ab_S.df           = lug_val(25); 
    xmlinput.iq.correction.fit.s1ab_S.normr        = lug_val(26); 

    LUXSubmitIQ(user,'s1ab_z',algorithm_name,version,file_id,cp,xmlinput,'Kr83 Calibration','',Image_paths);
    clear lug_val Image_paths xmlinput      
    
%% Kr s1a/s1b XY maps

    %Added after submit == 1 was hard coded, so this doesn't check for that condition
    lug_val(1)  = {file_id};
    lug_val(2)  = {user};
    lug_val(3)  = {algorithm_name};
    lug_val(4)  = {version};   
    lug_val(5)  = {cp}; %Proccessing version number
    lug_val(6)  = {gs}; %Proccessing gs number
    lug_val(7)  = {Kr_events};   
    lug_val(8)  = {s1a_xy_mean};
    lug_val(9)  = {s1a_xy_mean_err};
    lug_val(10) = {s1b_xy_mean};
    lug_val(11) = {s1b_xy_mean_err};
    lug_val(12) = {s1ab_xy_mean};
    lug_val(13) = {s1ab_xy_mean_err};
    lug_val(14) = {s1ab_xbins};
    lug_val(15) = {s1ab_ybins};
    lug_val(16) = {x_center};
    lug_val(17) = {y_center};
    lug_val(18) = {z_center};

    Image_paths = {strcat('LUX_corrections/',file_id_cp,'/s1a_xy',file_id_cp,'.jpg'), strcat('LUX_corrections/',file_id_cp,'/s1b_xy',file_id_cp,'.jpg')...
                  ,strcat('LUX_corrections/',file_id_cp,'/s1ab_ratio_xy',file_id_cp,'.jpg')}; 

    xmlinput.iq.global.filename_prefix             = lug_val(1);          
    xmlinput.iq.global.computed_by                 = lug_val(2);
    xmlinput.iq.global.computed_date               = datestr(now,'yyyymmddTHHMM');
    xmlinput.iq.global.source                      = source;
    xmlinput.iq.global.source_id                   = 'unknown';
    xmlinput.iq.global.notes                       = ' ';
    xmlinput.iq.global.algorithm_name              = lug_val(3);
    xmlinput.iq.global.algorithm_version           = lug_val(4);
    xmlinput.iq.global.cp_number                   = lug_val(5);
    xmlinput.iq.global.gs_number                   = lug_val(6);
    xmlinput.iq.correction.fit.number_of_kr_events = lug_val(7);
    xmlinput.iq.correction.fit.s1a_xymeans         = lug_val(8);
    xmlinput.iq.correction.fit.s1a_xymeans_err     = lug_val(9);
    xmlinput.iq.correction.fit.s1b_xymeans         = lug_val(10);
    xmlinput.iq.correction.fit.s1b_xymeans_err     = lug_val(11);
    xmlinput.iq.correction.fit.s1ab_xyratios       = lug_val(12);
    xmlinput.iq.correction.fit.s1ab_xyratios_err   = lug_val(13);
    xmlinput.iq.correction.fit.s1ab_xbincenters    = lug_val(14);
    xmlinput.iq.correction.fit.s1ab_ybincenters    = lug_val(15);
    xmlinput.iq.correction.fit.x_center            = lug_val(16);
    xmlinput.iq.correction.fit.y_center            = lug_val(17);
    xmlinput.iq.correction.fit.z_center            = lug_val(18);

    LUXSubmitIQ(user,'s1ab_xy',algorithm_name,version,file_id,cp,xmlinput,'Kr83 Calibration','',Image_paths);

    clear lug_val Image_paths xmlinput    
    
 %% Kr s1a/s1b XYZ maps

    %Added after submit == 1 was hard coded, so this doesn't check for that condition
    if length(s1ab_z)>100000; %require 100,000 events to do 3D binning
        lug_val(1)  = {file_id};
        lug_val(2)  = {user};
        lug_val(3)  = {algorithm_name};
        lug_val(4)  = {version};   
        lug_val(5)  = {cp}; %Proccessing version number
        lug_val(6)  = {gs}; %Proccessing gs number
        lug_val(7)  = {Kr_events};   
        lug_val(8)  = {s1a_xyz_mean};
        lug_val(9)  = {s1a_xyz_mean_err};
        lug_val(10) = {s1b_xyz_mean};
        lug_val(11) = {s1b_xyz_mean_err};
        lug_val(12) = {s1ab_xyz_mean};
        lug_val(13) = {s1ab_xyz_mean_err};
        lug_val(14) = {s1ab_xyz_xbins};
        lug_val(15) = {s1ab_xyz_ybins};
        lug_val(16) = {s1ab_xyz_zbins};
        lug_val(17) = {x_center};
        lug_val(18) = {y_center};
        lug_val(19) = {z_center};
        lug_val(20) = {s1as1b_at_center};

        xmlinput.iq.global.filename_prefix             = lug_val(1);          
        xmlinput.iq.global.computed_by                 = lug_val(2);
        xmlinput.iq.global.computed_date               = datestr(now,'yyyymmddTHHMM');
        xmlinput.iq.global.source                      = source;
        xmlinput.iq.global.source_id                   = 'unknown';
        xmlinput.iq.global.notes                       = ' ';
        xmlinput.iq.global.algorithm_name              = lug_val(3);
        xmlinput.iq.global.algorithm_version           = lug_val(4);
        xmlinput.iq.global.cp_number                   = lug_val(5);
        xmlinput.iq.global.gs_number                   = lug_val(6);
        xmlinput.iq.correction.fit.number_of_kr_events = lug_val(7);
        xmlinput.iq.correction.fit.s1a_xyzmeans        = lug_val(8);
        xmlinput.iq.correction.fit.s1a_xyzmeans_err    = lug_val(9);
        xmlinput.iq.correction.fit.s1b_xyzmeans        = lug_val(10);
        xmlinput.iq.correction.fit.s1b_xyzmeans_err    = lug_val(11);
        xmlinput.iq.correction.fit.s1ab_xyzratios      = lug_val(12);
        xmlinput.iq.correction.fit.s1ab_xyzratios_err  = lug_val(13);
        xmlinput.iq.correction.fit.s1ab_xbincenters    = lug_val(14);
        xmlinput.iq.correction.fit.s1ab_ybincenters    = lug_val(15);
        xmlinput.iq.correction.fit.s1ab_zbincenters    = lug_val(16);
        xmlinput.iq.correction.fit.x_center            = lug_val(17);
        xmlinput.iq.correction.fit.y_center            = lug_val(18);
        xmlinput.iq.correction.fit.z_center            = lug_val(19);
        xmlinput.iq.correction.fit.s1as1b_at_center    = lug_val(20);

        LUXSubmitIQ(user,'s1ab_3D',algorithm_name,version,file_id,cp,xmlinput,'Kr83 Calibration','','');

    end
    clear lug_val Image_paths xmlinput     

%% S1a/S1b R^2vZ and Field Maps

    %Added after submit == 1 was hard coded, so this doesn't check for that condition
    lug_val(1)  = {file_id};
    lug_val(2)  = {user};
    lug_val(3)  = {algorithm_name};
    lug_val(4)  = {version};   
    lug_val(5)  = {cp}; %Proccessing version number
    lug_val(6)  = {gs}; %Proccessing gs number
    lug_val(7)  = {Kr_events};   
    lug_val(8)  = {s1ab_r2bins};
    lug_val(9)  = {s1ab_zbins};
    lug_val(10) = {s1a_r2z_mean};
    lug_val(11) = {s1a_r2z_mean_err};
    lug_val(12) = {s1b_r2z_mean};
    lug_val(13) = {s1b_r2z_mean_err};
    lug_val(14) = {s1ab_r2z_mean};
    lug_val(15) = {s1ab_r2z_mean_err};
    lug_val(16) = {field_map_2015};
    lug_val(17) = {field_map_2015_staterr};
    lug_val(18) = {field_map_2014};
    lug_val(19) = {field_map_2014_staterr};
    lug_val(20) = {mean_field_map};
    lug_val(21) = {field_map_syserr};
    lug_val(22) = {mean_field_map_totalerr};

    Image_paths = {strcat('LUX_corrections/',file_id_cp,'/s1a_r2z_fig',file_id_cp,'.jpg'), strcat('LUX_corrections/',file_id_cp,'/s1b_r2z_fig',file_id_cp,'.jpg')...
                  ,strcat('LUX_corrections/',file_id_cp,'/s1ab_r2z_fig',file_id_cp,'.jpg'),strcat('LUX_corrections/',file_id_cp,'/field_r2z_fig',file_id_cp,'.jpg')...
                  ,strcat('LUX_corrections/',file_id_cp,'/field_frac_error_r2z_fig',file_id_cp,'.jpg')}; 

    xmlinput.iq.global.filename_prefix                 = lug_val(1);          
    xmlinput.iq.global.computed_by                     = lug_val(2);
    xmlinput.iq.global.computed_date                   = datestr(now,'yyyymmddTHHMM');
    xmlinput.iq.global.source                          = source;
    xmlinput.iq.global.source_id                       = 'unknown';
    xmlinput.iq.global.notes                           = ' ';
    xmlinput.iq.global.algorithm_name                  = lug_val(3);
    xmlinput.iq.global.algorithm_version               = lug_val(4);
    xmlinput.iq.global.cp_number                       = lug_val(5);
    xmlinput.iq.global.gs_number                       = lug_val(6);
    xmlinput.iq.correction.fit.number_of_kr_events     = lug_val(7);
    xmlinput.iq.correction.fit.s1ab_r2bins             = lug_val(8);
    xmlinput.iq.correction.fit.s1ab_zbins              = lug_val(9);
    xmlinput.iq.correction.fit.s1a_r2z_mean            = lug_val(10);
    xmlinput.iq.correction.fit.s1a_r2z_mean_err        = lug_val(11);
    xmlinput.iq.correction.fit.s1b_r2z_mean            = lug_val(12);
    xmlinput.iq.correction.fit.s1b_r2z_mean_err        = lug_val(13);
    xmlinput.iq.correction.fit.s1ab_r2z_mean           = lug_val(14);
    xmlinput.iq.correction.fit.s1ab_r2z_mean_err       = lug_val(15);                        
    xmlinput.iq.correction.fit.field_map_2015          = lug_val(16);
    xmlinput.iq.correction.fit.field_map_2015_staterr  = lug_val(17); 
    xmlinput.iq.correction.fit.field_map_2014          = lug_val(18); 
    xmlinput.iq.correction.fit.field_map_2014_staterr  = lug_val(19);            
    xmlinput.iq.correction.fit.mean_field_map          = lug_val(20);
    xmlinput.iq.correction.fit.field_map_syserr        = lug_val(21); 
    xmlinput.iq.correction.fit.mean_field_map_totalerr = lug_val(22); 


    LUXSubmitIQ(user,'S1ab_Field',algorithm_name,version,file_id,cp,xmlinput,'Kr83 Calibration','',Image_paths);

    clear lug_val Image_paths xmlinput      
    
end

clear;

fprintf('KrypCal Finished \n');

end

