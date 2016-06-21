%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This script compares Kr S1a/S1b data to Scott's field simulations in Sep2015 so that
%       we can use S1a/S1b as a measure of field in Kr calibrations. Sorry for the mess...
%       this is mostly pieced together from other code I've already written.  I doubt anybody else
%       will look at this script anyway.
% 
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loading Data

path='C:\Users\Richard\Desktop\LUX Analysis\lux10_20150929T1905_cp17548';
   
rqs_to_load = {'pulse_area_phe','event_timestamp_samples'...
   ,'pulse_classification' ...
   ,'z_drift_samples' , 's1s2_pairing','golden'...
   ,'selected_s1_s2','top_bottom_ratio','x_cm','y_cm'...
   ,'full_evt_area_phe',...
   'event_number','chi2','prompt_fraction','aft_t1_samples','pulse_start_samples',...
   'pulse_end_samples','top_bottom_asymmetry','aft_t0_samples','aft_t2_samples',...
   'full_evt_area_phe','admin','Kr83fit_s1a_area_phe','Kr83fit_s1b_area_phe',...
   'hft_t50r_samples','hft_t50l_samples','Kr83fit_dt_samples'};

d = LUXLoadMultipleRQMs_framework(path,rqs_to_load);  

distinguish_string='test';

%Hardcoding inputs for ease of use during testing
user_name='RichardKnoche';
dp_version='2.15_test';
algorithm_name='test_KrypCal';
submit=0;
delay_min=0;
inpaint_on=1;

file=char(d.files_loaded(1));
file_id=file(1:19);
cp_start=strfind(file,'cp')+2;%find location of 'cp' in string
cp=file(cp_start:cp_start+4);%cp number has four digits
file_id_cp = [file_id '_cp' cp '_' distinguish_string];


gs=0;
   
grid_size=3; %cm XY plane
rcut_min = 0;
rcut_max = 25;%cm
s1area_bound_min = 50;%100
s1area_bound_max = 650;%600
s2area_bound_min = 1000;%200
s2area_bound_max = 32000;%30000
min_evts_per_bin = 300;
max_num_bins = 65;

S1xybinsize = grid_size;%cm
S1_xbin_min = -25;
S1_xbin_max = 25;
S1_ybin_min = -25;
S1_ybin_max = 25;


SE_xbin_min = -25;
SE_xbin_max = 25;
SE_ybin_min = -25;
SE_ybin_max = 25;
% SExybinsize = grid_size; Defined later, based on number of SE events

S2xybinsize = grid_size;%cm
S2_xbin_min = -25;
S2_xbin_max = 25;
S2_ybin_min = -25;
S2_ybin_max = 25;

 
d.z_drift_samples(isnan(d.z_drift_samples)) = 0.0; % get rid of NaN        
events=sort(size(d.pulse_area_phe)); %#ok<*NASGU> %The first element is the number of pulses. sometimes it's 5, sometimes it's 10

s1_area_cut= inrange(d.pulse_area_phe,[s1area_bound_min,s1area_bound_max]);
s2_area_cut= inrange(d.pulse_area_phe,[s2area_bound_min,s2area_bound_max]);

s1_class=(d.pulse_classification==1 )& s1_area_cut ; %area cut for Kr events
s2_class=(d.pulse_classification==2) & s2_area_cut ;   
s4_class=(d.pulse_classification==4) ;

events=sort(size(d.pulse_area_phe)); %The first element is the number of pulses. sometimes it's 5, sometimes it's 10
cut_pulse_s1 = d.pulse_classification == 1 | d.pulse_classification == 9;
cut_pulse_s2 = d.pulse_classification == 2;
cut_s2_with_threshold = d.pulse_area_phe.*cut_pulse_s2 > 100; % subset of cut_pulse_s2
cut_legit_s2_in_legit_event = d.s1s2_pairing.*cut_s2_with_threshold; % this should be used as s2 classification cuts
cut_golden_event = sum(cut_legit_s2_in_legit_event) == 1; %defines golden events to be events which have one and only one paired S2 above the threshold of 100 phe - there can be multiple S1s still
cut_s2_in_golden_events = logical(repmat(cut_golden_event,[10,1]).*cut_legit_s2_in_legit_event); %Selects S2 that is in a golden event
cut_s1_in_golden_events = logical(repmat(cut_golden_event,[10,1]).*cut_pulse_s1.*d.s1s2_pairing); %Selects first S1 that is in a golden event
% select Kr83 events with cut on S2
cut_s2_area = inrange(d.pulse_area_phe, [s2area_bound_min, s2area_bound_max]);
cut_s1_area = inrange(d.pulse_area_phe, [s1area_bound_min, s1area_bound_max]);
cut_s2_for = cut_s2_in_golden_events.*cut_s2_area; %Selects S2 that is in a golden event and in Kr area bounds
cut_s1_for = cut_s1_in_golden_events.*cut_s1_area; %Selects first S1 that is in a golden event and in Kr area bounds
cut_selected_events = sum(cut_s2_for) == 1 & sum(cut_s1_for) == 1 & sum(d.pulse_classification==1)==1; %Requires that "good" golden events have only one S1, that the S1 be within area bounds, and the S2 be within area bounds
%Note sum(cut_s1_for) == 1 above only requires that the first of the S1 in an event be within area bounds, since the S1S2pairing part of cut_s1_in_golden_events is 0 for all subsequent S1s in the events
s1_single_cut = logical(repmat(cut_selected_events,[10,1]).*cut_s1_in_golden_events);
s2_single_cut = logical(repmat(cut_selected_events,[10,1]).*cut_s2_in_golden_events);
 
%% Set up cuts and variables
s1ab_cut=logical(sum(s1_single_cut(:,:),1));
s1a_phe_both=d.Kr83fit_s1a_area_phe(s1ab_cut);
s1b_phe_both=d.Kr83fit_s1b_area_phe(s1ab_cut);
s1ab_timing=d.Kr83fit_dt_samples(s1ab_cut);
s1ab_x=d.x_cm(s2_single_cut);
s1ab_y=d.y_cm(s2_single_cut);
s1ab_z=d.z_drift_samples(s2_single_cut)/100;
s1ab_radius=sqrt(s1ab_x.^2+s1ab_y.^2);
s1ab_timing_cut=s1ab_timing>13 & s1b_phe_both>30 & s1ab_radius.'<25;

s1ab_x=s1ab_x(s1ab_timing_cut);
s1ab_y=s1ab_y(s1ab_timing_cut);
s1ab_z=s1ab_z(s1ab_timing_cut);
s1a_phe_both=s1a_phe_both(s1ab_timing_cut);
s1b_phe_both=s1b_phe_both(s1ab_timing_cut);
s1ab_radius=sqrt(s1ab_x.^2+s1ab_y.^2);

%% S1a/b R^2 v Z map 
s1ab_r2bin_min=0;
s1ab_r2bin_max=25*25;
s1ab_zbin_min=10;
s1ab_zbin_max=330;

s1ab_r2z_binnum=(length(s1ab_z))/300;
s1ab_r2_binnum=floor(sqrt(s1ab_r2z_binnum))-1;
s1ab_z_binnum=floor(sqrt(s1ab_r2z_binnum));

s1ab_r2_binsize=(s1ab_r2bin_max-s1ab_r2bin_min)/s1ab_r2_binnum;
s1ab_z_binsize=(s1ab_zbin_max-s1ab_zbin_min)/s1ab_z_binnum;

s1ab_r2bins=s1ab_r2bin_min+s1ab_r2_binsize/2:s1ab_r2_binsize:s1ab_r2bin_max;
s1ab_zbins=s1ab_zbin_min+s1ab_z_binsize/2:s1ab_z_binsize:s1ab_zbin_max; 

s1a_r2z_mean=zeros(length(s1ab_zbins),length(s1ab_r2bins));
s1b_r2z_mean=zeros(length(s1ab_zbins),length(s1ab_r2bins));
s1a_r2z_mean_err=zeros(length(s1ab_zbins),length(s1ab_r2bins));
s1b_r2z_mean_err=zeros(length(s1ab_zbins),length(s1ab_r2bins));
field_map_2014=zeros(length(s1ab_zbins),length(s1ab_r2bins));
field_map_2015=zeros(length(s1ab_zbins),length(s1ab_r2bins));
field_map_error=zeros(length(s1ab_zbins),length(s1ab_r2bins));

r2_counter=1;
 for r2_bin_max=s1ab_r2bin_min+s1ab_r2_binsize:s1ab_r2_binsize:s1ab_r2bin_max;
     z_counter=1;
        for z_bin_max=s1ab_zbin_min+s1ab_z_binsize:s1ab_z_binsize:s1ab_zbin_max;          
           
          bin_cut = s1a_phe_both>0 & s1b_phe_both>0 & inrange(s1ab_radius.^2,[r2_bin_max-s1ab_r2_binsize,r2_bin_max]).' & inrange(s1ab_z,[z_bin_max-s1ab_z_binsize,z_bin_max]).'; 
            if length(s1a_phe_both(bin_cut))>100;
                clear s1a_confint  s1b_confint
                s1a_r2z_fit=fit([0:2:400].',hist(s1a_phe_both(bin_cut & inrange(s1a_phe_both,[0 400])),[0:2:400]).','gauss1'); %#ok<*NBRAK>
                s1a_r2z_mean(z_counter,r2_counter)=s1a_r2z_fit.b1;
                s1a_confint=confint(s1a_r2z_fit,0.68);
                s1a_r2z_mean_err(z_counter,r2_counter)=abs(s1a_confint(1,2)-s1a_r2z_fit.b1);
                field_map_2014(z_counter,r2_counter)=EstimateField_Sep2014((r2_bin_max-s1ab_r2_binsize/2)^(1/2),z_bin_max-s1ab_z_binsize/2);
                field_map_2015(z_counter,r2_counter)=EstimateField_Sep2015((r2_bin_max-s1ab_r2_binsize/2)^(1/2),z_bin_max-s1ab_z_binsize/2);
                field_map_error(z_counter,r2_counter)=std([EstimateField_Sep2015((r2_bin_max-s1ab_r2_binsize/2)^(1/2),z_bin_max-s1ab_z_binsize/2) EstimateField_Sep2014((r2_bin_max-s1ab_r2_binsize/2)^(1/2),z_bin_max-s1ab_z_binsize/2)]);
                
                s1b_r2z_fit=fit([0:2:300].',hist(s1b_phe_both(bin_cut & inrange(s1b_phe_both,[0 300])),[0:2:300]).','gauss1');
                s1b_r2z_mean(z_counter,r2_counter)=s1b_r2z_fit.b1;
                s1b_confint=confint(s1b_r2z_fit,0.68);
                s1b_r2z_mean_err(z_counter,r2_counter)=abs(s1b_confint(1,2)-s1b_r2z_fit.b1);
            else
                s1a_r2z_mean(z_counter,r2_counter)=0;
                s1a_r2z_mean_err(z_counter,r2_counter)=0;
                s1b_r2z_mean(z_counter,r2_counter)=0;
                s1b_r2z_mean_err(z_counter,r2_counter)=0; 
                field_map_2014(z_counter,r2_counter)=0;
                field_map_2015(z_counter,r2_counter)=0;
                
            end
           z_counter=z_counter+1;     
        end
        r2_counter=r2_counter+1;
 end
s1ab_r2z_mean=s1a_r2z_mean./s1b_r2z_mean;
s1ab_r2z_mean_err=sqrt( (s1b_r2z_mean_err.*s1a_r2z_mean./(s1b_r2z_mean.^2)).^2 + (s1a_r2z_mean_err./s1b_r2z_mean).^2);
   
%Plot s1a R2Z means
color_range_max=max(max(s1a_r2z_mean));
color_range_min=min(min(s1a_r2z_mean(s1a_r2z_mean>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;

s1a_r2z_fig = figure;
contourf(s1ab_r2bins,s1ab_zbins,s1a_r2z_mean,vc,'LineColor','none');
xlabel('R^2 (cm^2)','fontsize', 18);
ylabel('Z (uSec)','fontsize', 18);
title(strcat('2015 S1a Mean vs. R^2 v Z.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'S1a Mean (phe)','FontSize',16)
set(g,'FontSize',18)
hold on
axis([0 625 0 330])
myfigview(16);
set(gca,'Ydir','reverse')

%Plot s1a R2Z means
color_range_max=max(max(s1b_r2z_mean));
color_range_min=min(min(s1b_r2z_mean(s1b_r2z_mean>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;

s1b_r2z_fig = figure;
contourf(s1ab_r2bins,s1ab_zbins,s1b_r2z_mean,vc,'LineColor','none');
xlabel('R^2 (cm^2)','fontsize', 18);
ylabel('Z (uSec)','fontsize', 18);
title(strcat('2015 S1b Mean vs. R^2 v Z.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'S1b Mean (phe)','FontSize',16)
set(g,'FontSize',18)
hold on
axis([0 625 0 330])
myfigview(16);
set(gca,'Ydir','reverse')

%s1a/b R2Z Ratio
color_range_max=max(max(s1ab_r2z_mean));
color_range_min=min(min(s1ab_r2z_mean(s1ab_r2z_mean>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;

s1ab_r2z_fig = figure;
contourf(s1ab_r2bins,s1ab_zbins,s1ab_r2z_mean,vc,'LineColor','none');
xlabel('R^2 (cm^2)','fontsize', 18);
ylabel('Z (uSec)','fontsize', 18);
title(strcat('2015 S1a/S1b Mean vs. R^2 v Z.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'S1a/b Ratio','FontSize',16)
set(g,'FontSize',18)
hold on
axis([0 625 0 330])
myfigview(16);
set(gca,'Ydir','reverse')

%Field R2Z Map
color_range_max=max(max(field_map_2015));
color_range_min=min(min(field_map_2015(field_map_2015>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;

field_map_fig = figure;
contourf(s1ab_r2bins,s1ab_zbins,field_map_2015,vc,'LineColor','none');
xlabel('R^2 (cm^2)','fontsize', 18);
ylabel('Z (uSec)','fontsize', 18);
title(strcat('2015 Field Value vs. R^2 v Z.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'Field (V/cm)','FontSize',16)
set(g,'FontSize',18)
hold on
axis([0 625 0 330])
myfigview(16);
set(gca,'Ydir','reverse')    
    
%% Comparing S1a/S1b RvZ to Scott's field map
%     figure
%     errorbar(field_map_2015(:),s1ab_r2z_mean(:),s1ab_r2z_mean_err(:),'.k');
%     xlabel('Field (V/cm)'); ylabel('S1a/S1b Ratio'); myfigview(16);
%     
%     figure
%     rkploterr(s1ab_r2z_mean(:),field_map_2015(:),s1ab_r2z_mean_err(:),field_map_error(:),[1 0 0],[],[],2);
%     ylabel('Field (V/cm)'); xlabel('S1a/S1b Ratio'); myfigview(16);
    
%% Fit a polynomial
% Call fminsearch with a random starting point.
cut= ~isnan(s1ab_r2z_mean) & ~isnan(field_map_2015) & s1ab_r2z_mean>0 & field_map_2015>0;
data.xdata=s1ab_r2z_mean(cut);
data.ydata=field_map_2015(cut);
data.xdataerr=s1ab_r2z_mean_err(cut);
data.ydataerr=field_map_error(cut);
[S1aS1b_fit_2015, S1aS1b_fit_2015_error]=polyfit(data.xdata,data.ydata,2); %given X
%example usage [a b]=polyval(S1aS1b_fit_2015,[2.2:0.01:3],S1aS1b_fit_2015_error);

%This is used for fitting and taking x and y errors into account... but since errors are underestimated on some points I'm not going to factor them in
% 
% inv_sigma2=@(theta,data) 1./(data.ydataerr.^2);
% % inv_sigma2 = @(theta,data) 1./( ( (2.*theta(1).*data.xdata+theta(2)).*data.xdataerr).^2 + data.ydataerr.^2);
% 
% modelfun = @(x,theta) theta(1).*x.^2+theta(2).*x+theta(3);
% ssfun = @(theta,data) 0.5.*( sum( (data.ydata - modelfun(data.xdata,theta)).^2.*inv_sigma2(theta,data) - log(inv_sigma2(theta,data)) ));
% polyval(data.xdata,data.ydata,3)
% 
% init=polyfit(data.xdata,data.ydata,2);
% [tmin,ssmin]=fminsearch(ssfun,[init(1);init(2);init(3)],[],data); %minimizing sum of squares

%% Comparing S1a/S1b RvZ to Scott's field map
%     figure
%     errorbar(field_map_2014(:),s1ab_r2z_mean(:),s1ab_r2z_mean_err(:),'.k');
%     xlabel('Field (V/cm)'); ylabel('S1a/S1b Ratio'); myfigview(16);
%     
    figure
    rkploterr(s1ab_r2z_mean(:),field_map_2015(:),s1ab_r2z_mean_err(:),field_map_error(:),[ ],[],[],2);
    ylabel('Field (V/cm)'); xlabel('S1a/S1b Ratio'); myfigview(16);
    x=[0:0.1:3];
    y=polyval(S1aS1b_fit_2015,x);
    plot(x,y,'Color',[0.6 0.6 0.6])
    
%% saving

save('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal\KrypCal_2p20_Code\2015_Field_to_S1aS1b','S1aS1b_fit_2015','S1aS1b_fit_2015_error');
