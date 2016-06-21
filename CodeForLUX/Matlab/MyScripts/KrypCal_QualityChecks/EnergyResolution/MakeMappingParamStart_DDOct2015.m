path='C:\Users\Richard\Desktop\LUX Analysis\lux10_20151014T0010_cp17711';

   
rqs_to_load = {'pulse_area_phe','event_timestamp_samples'...
   ,'pulse_classification' ...
   ,'z_drift_samples' , 's1s2_pairing'...
   ,'top_bottom_ratio','x_cm','y_cm'...
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
    events=sort(size(d.pulse_area_phe)); %The first element is the number of pulses. sometimes it's 5, sometimes it's 10
          
    s1_area_cut= inrange(d.pulse_area_phe,[s1area_bound_min,s1area_bound_max]);
    s2_area_cut= inrange(d.pulse_area_phe,[s2area_bound_min,s2area_bound_max]);
    
    s1_class=(d.pulse_classification==1 )& s1_area_cut ; %area cut for Kr events
    s2_class=(d.pulse_classification==2) & s2_area_cut ;   
    s4_class=(d.pulse_classification==4) ;
   

%     s1_single_cut =logical( (s1_class & d.golden & d.selected_s1_s2).*repmat(sum(s2_class & d.golden & d.selected_s1_s2)==1,events(1),1) ); % was s1_before_s2_cut before using golden
%     s2_single_cut =logical( (s2_class & d.golden & d.selected_s1_s2).*repmat(sum(s1_class & d.golden & d.selected_s1_s2)==1,events(1),1) );


 
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

    drift_time = d.z_drift_samples(s2_single_cut)/100;  % us
        
    d.phe_bottom=d.pulse_area_phe./(1+d.top_bottom_ratio); %bottom PMT pulse area

    s1_phe_both = d.pulse_area_phe(s1_single_cut);
    s1_phe_bottom = d.phe_bottom(s1_single_cut);
    
    s2_phe_both = d.pulse_area_phe(s2_single_cut);
    s2_phe_bottom = d.phe_bottom(s2_single_cut);
    
    s2x = d.x_cm(s2_single_cut);
    s2y = d.y_cm(s2_single_cut);   

    d.livetime_sec=sum(d.livetime_end_samples-d.livetime_latch_samples)/1e8;
    evt_cut=logical(sum(s2_single_cut));%Cut for all the events passing the single S1 & S2 cut
    event_timestamp_samples=d.event_timestamp_samples(evt_cut);
    
    time_wait_cut=event_timestamp_samples/1e8/60 > delay_min; %allow x min for Kr mixing
    
    s2_width=(d.aft_t2_samples(s2_single_cut)-d.aft_t0_samples(s2_single_cut)); %cut above 800 samples
    s1_width=d.pulse_end_samples(s1_single_cut)-d.pulse_start_samples(s1_single_cut);
    
    s2radius = (s2x.^2+s2y.^2).^(0.5);
    clean_cut = (log10(s1_phe_both)+0.6*log10(s2_phe_both)>4.5) & s2radius<25;

         s1_phe_bottom=s1_phe_bottom(clean_cut);
         s1_phe_both=s1_phe_both(clean_cut);
         s2_phe_bottom=s2_phe_bottom(clean_cut);
         s2_phe_both=s2_phe_both(clean_cut);
         drift_time=drift_time(clean_cut);
         s2x=s2x(clean_cut);
         s2y=s2y(clean_cut);
         s2radius = (s2x.^2+s2y.^2).^(0.5);
         event_timestamp_samples=event_timestamp_samples(clean_cut);
         s2_width=s2_width(clean_cut);
         s1_width=s1_width(clean_cut);

kr_energy_cut = inrange(s1_phe_both,[s1area_bound_min s1area_bound_max]) & inrange(s2_phe_both,[s2area_bound_min s2area_bound_max]);%for counting Kr events
Kr_events=length(drift_time(kr_energy_cut));% Are there enough events for the XY or XYZ map?

%%%% Detector Center Detection %%%%%%%
x_center=mean(s2x(drift_time>4)); %exclude extraction field region at uSec<4
y_center=mean(s2y(drift_time>4));
z_center=mean(drift_time(drift_time>4));
det_edge=330;

%%
    zcut_min = 10; %field changes at 4 uSec
    zcut_max = 0.95*det_edge;


    if Kr_events < 30000 %then create a 10x10 grid. Otherwise 25x25. Want about 30 evts per bin
           
            S1xybinsize = sqrt((50*50)/(Kr_events/150));
            S1_xbin_min = -25;
            S1_xbin_max = 25;
            S1_ybin_min = -25;
            S1_ybin_max = 25;
            SE_ybin_min = -25;
            SE_ybin_max = 25;

            S2xybinsize = sqrt((50*50)/(Kr_events/150));%cm
%             SExybinsize = 5; DEFINED LATER, based on number of SE events
            S2_xbin_min = -25;
            S2_xbin_max = 25;
            S2_ybin_min = -25;
            S2_ybin_max = 25;
            SE_xbin_min = -25;
            SE_xbin_max = 25;
            
                    
    end
    
    s1_xbins = (S1_xbin_max-S1_xbin_min)./S1xybinsize;
    s1_ybins = (S1_ybin_max-S1_ybin_min)./S1xybinsize;

    s2_xbins = (S2_xbin_max-S2_xbin_min)./S2xybinsize;
    s2_ybins = (S2_ybin_max-S2_ybin_min)./S2xybinsize;
        

 if Kr_events < 4000; % skip krpCal if less than 4000 events
     calibrations=zeros(6,1);
 end
 
%%
%Set up cuts and variables
    s1ab_cut=logical(sum(s1_single_cut(:,:),1));
    s1a_phe_both=d.Kr83fit_s1a_area_phe(s1ab_cut);
    s1b_phe_both=d.Kr83fit_s1b_area_phe(s1ab_cut);
    s1ab_timing=d.Kr83fit_dt_samples(s1ab_cut);
    s1ab_x=d.x_cm(s2_single_cut);
    s1ab_y=d.y_cm(s2_single_cut);
    s1ab_z=d.z_drift_samples(s2_single_cut)/100;
    s1ab_radius=sqrt(s1ab_x.^2+s1ab_y.^2);

    s1ab_timing_cut=s1ab_timing>13 & s1b_phe_both>30 & s1ab_radius.'<=25;
       %plotting x=s1a v y=s1a/s1b cuts: y<0.035x-1, y>-0.035.*x+5; y>0.01.*x; %y<4.2; x<255
%     s1ab_areacut= (s1a_phe_both./s1b_phe_both-0.035.*s1a_phe_both<-1) & (s1a_phe_both./s1b_phe_both+0.035.*s1a_phe_both>5)...
%         & (s1a_phe_both./s1b_phe_both > 0.01.* s1a_phe_both) & (s1a_phe_both./s1b_phe_both<4.2) & (s1a_phe_both<255);
    
    s1ab_x=s1ab_x(s1ab_timing_cut);
    s1ab_y=s1ab_y(s1ab_timing_cut);
    s1ab_z=s1ab_z(s1ab_timing_cut);
    s1a_phe_both=s1a_phe_both(s1ab_timing_cut);
    s1b_phe_both=s1b_phe_both(s1ab_timing_cut);
    
    dT_step=det_edge/(length(s1ab_z)/1000); %2000 events per z bin
    
    
    clear s1a_z_means s1a_z_means_err s1b_z_means s1b_z_means_err s1ab_bincenters s1ab_z_means s1ab_z_means_err
%Z Dependence of S1a/S1b
i=1;
    for z_max=10+dT_step:dT_step:det_edge;
        s1a_fit=fit([0:2:400].',hist(s1a_phe_both(s1a_phe_both>0 & s1a_phe_both<400 & inrange(s1ab_z,[z_max-dT_step,z_max]).'),[0:2:400]).','gauss1');
        s1b_fit=fit([0:2:300].',hist(s1b_phe_both(s1b_phe_both>0 & s1b_phe_both<300 & inrange(s1ab_z,[z_max-dT_step,z_max]).'),[0:2:300]).','gauss1');
        s1a_z_means(i)=s1a_fit.b1;
        s1a_z_means_err(i)=s1a_fit.c1/sqrt(2)/sqrt(length(s1a_phe_both(s1a_phe_both>0 & s1a_phe_both<400 & inrange(s1ab_z,[z_max-dT_step,z_max]).')));
        s1b_z_means(i)=s1b_fit.b1;
        s1b_z_means_err(i)=s1b_fit.c1/sqrt(2)/sqrt(length(s1b_phe_both(s1b_phe_both>0 & s1b_phe_both<300 & inrange(s1ab_z,[z_max-dT_step,z_max]).')));
        s1ab_bincenters(i)=z_max-dT_step/2;
        i=i+1;
    end

s1ab_z_means=s1a_z_means./s1b_z_means;
s1ab_z_means_err=sqrt( (s1b_z_means_err.*s1a_z_means./(s1b_z_means.^2)).^2 + (s1a_z_means_err./s1b_z_means).^2);
    
%fit polynomial to s1az and s1bz, to remove z dependence later in xy fit
[s1a_P, s1a_S]=polyfit(s1ab_bincenters,s1a_z_means,3);
[s1b_P, s1b_S]=polyfit(s1ab_bincenters,s1b_z_means,3);
[s1ab_P, s1ab_S]=polyfit(s1ab_bincenters,s1ab_z_means,2);

%S1a z dependence plot    
s1az_plot=figure;
errorbar(s1ab_bincenters,s1a_z_means,s1a_z_means_err,'.k')   
xlabel('Drift Time (uSec)');ylabel('S1a Mean (phe)'); title('S1a Z Dependence'); myfigview(16);
hold on;
plot([10:1:320],polyval(s1a_P,[10:1:320]),'-r','LineWidth',2)


%S1b z dependence plot    
s1bz_plot=figure;
errorbar(s1ab_bincenters,s1b_z_means,s1b_z_means_err,'.k')   
xlabel('Drift Time (uSec)');ylabel('S1b Mean (phe)'); title('S1b Z Dependence'); myfigview(16);
hold on;
plot([10:1:320],polyval(s1b_P,[10:1:320]),'-r','LineWidth',2)
  
%S1a/b z dependence plot
s1a_over_bz_plot=figure;
errorbar(s1ab_bincenters,s1ab_z_means,s1ab_z_means_err,'.k')   
xlabel('Drift Time (uSec)');ylabel('S1a/b'); title('S1 a/b Z Dependence'); myfigview(16);
hold on;
plot([10:1:320],polyval(s1ab_P,[10:1:320]),'-r','LineWidth',2)
  
%Correcting z dependence to get XY dependence in same manner that Kr/CH3T
%map does
    s1a_phe_both_z=s1a_phe_both.'.*polyval(s1a_P,z_center)./polyval(s1a_P,s1ab_z);
    s1b_phe_both_z=s1b_phe_both.'.*polyval(s1b_P,z_center)./polyval(s1b_P,s1ab_z);

%% S1a/b XY map (after z correction)    
s1ab_xbin_min=-25;
s1ab_ybin_min=-25;
s1ab_xbin_max=25;
s1ab_ybin_max=25;
s1ab_xybinsize=sqrt((50*50)/(length(s1ab_z)/300));
clear s1a_xy_mean s1a_xy_mean_err s1b_xy_mean s1b_xy_mean_err

s1ab_xbins=s1ab_xbin_min+s1ab_xybinsize/2:s1ab_xybinsize:s1ab_xbin_max;
s1ab_ybins=s1ab_ybin_min+s1ab_xybinsize/2:s1ab_xybinsize:s1ab_ybin_max; 

s1a_xy_mean=zeros(length(s1ab_ybins),length(s1ab_xbins));
s1b_xy_mean=zeros(length(s1ab_ybins),length(s1ab_xbins));

 for x_bin=s1ab_xbin_min:s1ab_xybinsize:(s1ab_xbin_max-s1ab_xybinsize/2);
        for y_bin=s1ab_ybin_min:s1ab_xybinsize:(s1ab_ybin_max-s1ab_xybinsize/2);          
          x_min = x_bin; x_max = x_bin+s1ab_xybinsize;
          y_min = y_bin; y_max = y_bin+s1ab_xybinsize;      

            x_count=int32(1+x_bin/s1ab_xybinsize+s1ab_xbin_max/s1ab_xybinsize);
            y_count=int32(1+y_bin/s1ab_xybinsize+s1ab_ybin_max/s1ab_xybinsize); 

          bin_cut = s1a_phe_both>0 & s1b_phe_both>0 & inrange(s1ab_x,[x_min,x_max]).' & inrange(s1ab_y,[y_min,y_max]).'; 
            if length(s1a_phe_both_z(bin_cut))>100;
                s1a_z_fit=fit([0:2:400].',hist(s1a_phe_both_z(bin_cut.' & inrange(s1a_phe_both_z,[0 400])),[0:2:400]).','gauss1');
                s1a_xy_mean(y_count,x_count)=s1a_z_fit.b1;
                s1a_xy_mean_err(y_count,x_count)=s1a_z_fit.c1/sqrt(2)/length(s1a_phe_both_z(bin_cut.' & inrange(s1a_phe_both_z,[0 400])));

                
                s1b_z_fit=fit([0:2:300].',hist(s1b_phe_both_z(bin_cut.' & inrange(s1b_phe_both_z,[0 300])),[0:2:300]).','gauss1');
                s1b_xy_mean(y_count,x_count)=s1b_z_fit.b1;
                s1b_xy_mean_err(y_count,x_count)=s1b_z_fit.c1/sqrt(2)/length(s1b_phe_both_z(bin_cut.' & inrange(s1b_phe_both_z,[0 300])));
            else
                s1a_xy_mean(y_count,x_count)=0;
                s1a_xy_mean_err(y_count,x_count)=0;
                s1b_xy_mean(y_count,x_count)=0;
                s1b_xy_mean_err(y_count,x_count)=0;   
            end
                
        end
 end
     s1ab_xy_mean=s1a_xy_mean./s1b_xy_mean;
     s1ab_xy_mean_err=sqrt( (s1b_xy_mean_err.*s1a_xy_mean./(s1b_xy_mean.^2)).^2 + (s1a_xy_mean_err./s1b_xy_mean).^2);
   

    %Plot s1a XY means
    color_range_max=max(max(s1a_xy_mean));
    color_range_min=min(min(s1a_xy_mean(s1a_xy_mean>0)));
    vc_step=(color_range_max-color_range_min)/50;
    vc=color_range_min:vc_step:color_range_max;
    
    s1a_xy_fig = figure;
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

    %Plot s1a XY means
    color_range_max=max(max(s1b_xy_mean));
    color_range_min=min(min(s1b_xy_mean(s1b_xy_mean>0)));
    vc_step=(color_range_max-color_range_min)/50;
    vc=color_range_min:vc_step:color_range_max;
    
    s1b_xy_fig = figure;
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
 
    %s1a/b XY Ratio
    color_range_max=max(max(s1ab_xy_mean));
    color_range_min=min(min(s1ab_xy_mean(s1ab_xy_mean>0)));
    vc_step=(color_range_max-color_range_min)/50;
    vc=color_range_min:vc_step:color_range_max;
    
    s1ab_xy_fig = figure;
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

  
%% 3D Kr s1a/s1b map

if length(s1ab_z)>90000; %require 100,000 events to do 3D binning
s1ab_xyz_numbins=floor((length(s1ab_z)/250)^(1/3)); %number of bins in one direction
s1ab_xyz_zstep=det_edge/s1ab_xyz_numbins;
s1ab_xyz_xstep=50/s1ab_xyz_numbins;
s1ab_xyz_ystep=50/s1ab_xyz_numbins;
s1ab_xyz_xmax=25;
r_max=25;

s1ab_xyz_zbins=10+s1ab_xyz_zstep/2:s1ab_xyz_zstep:10+s1ab_xyz_zstep/2+s1ab_xyz_zstep*s1ab_xyz_numbins; %20 us
s1ab_xyz_xbins=(-r_max+s1ab_xyz_xstep/2):s1ab_xyz_xstep:r_max;
s1ab_xyz_ybins=(-r_max+s1ab_xyz_ystep/2):s1ab_xyz_ystep:r_max;

s1a_xyz_mean=zeros(length(s1ab_xyz_ybins),length(s1ab_xyz_xbins),length(s1ab_xyz_zbins));
s1a_xyz_mean_err=zeros(length(s1ab_xyz_ybins),length(s1ab_xyz_xbins),length(s1ab_xyz_zbins));
s1b_xyz_mean=zeros(length(s1ab_xyz_ybins),length(s1ab_xyz_xbins),length(s1ab_xyz_zbins));
s1b_xyz_mean_err=zeros(length(s1ab_xyz_ybins),length(s1ab_xyz_xbins),length(s1ab_xyz_zbins));

for k = s1ab_xyz_zbins; %s1zbins are the center of the bin    
    for i = (-r_max+s1ab_xyz_ystep):s1ab_xyz_ystep:r_max %map to rows going down (y_cm) -- s1x bins are the max edge of the bin
        for j = (-r_max+s1ab_xyz_xstep):s1ab_xyz_xstep:r_max %map columns across (x_cm) -- s1y bins are the max edge of the bin
                       
            l=int8(i/s1ab_xyz_ystep+s1ab_xyz_xmax/s1ab_xyz_ystep);
            m=int8(j/s1ab_xyz_xstep+s1ab_xyz_xmax/s1ab_xyz_xstep);
            n=int8((k-(10+s1ab_xyz_zstep/2))/s1ab_xyz_zstep + 1);
            
           %sort, and make the cut. using the variable q
           q = s1ab_x<j & s1ab_x>(j-s1ab_xyz_xstep) & s1ab_y<i & s1ab_y>(i-s1ab_xyz_ystep) & inrange(s1ab_z,[(k-s1ab_xyz_zstep/2) , (k+s1ab_xyz_zstep/2)]); %no 0th element!
                    %swap x,y due to matrix x-row, y-column definition

            %Count the number of events per bin
            Count_S1ab_3D(l,m,n)=length(s1a_phe_both(q));
            
            if (Count_S1ab_3D(l,m,n) >= 100) % at least 100 counts before fitting. 
                s1a_z_fit=fit([0:2:400].',hist(s1a_phe_both(q.' & inrange(s1a_phe_both,[0 400])),[0:2:400]).','gauss1');
                s1a_xyz_mean(l,m,n)=s1a_z_fit.b1;
                s1a_xyz_mean_err(l,m,n)=s1a_z_fit.c1/sqrt(2)/length(s1a_phe_both(q.' & inrange(s1a_phe_both,[0 400])));

                
                s1b_z_fit=fit([0:2:300].',hist(s1b_phe_both(q.' & inrange(s1b_phe_both,[0 300]) ),[0:2:300]).','gauss1');
                s1b_xyz_mean(l,m,n)=s1b_z_fit.b1;
                s1b_xyz_mean_err(l,m,n)=s1b_z_fit.c1/sqrt(2)/length(s1b_phe_both(q.' & inrange(s1b_phe_both,[0 300])));
             else %not enough stats to do the fit
                s1a_xyz_mean(l,m,n)=0;
                s1a_xyz_mean_err(l,m,n)=0;
                s1b_xyz_mean(l,m,n)=0;
                s1b_xyz_mean_err(l,m,n)=0;
            end
                                      
        end
    end
k
end

     s1ab_xyz_mean=s1a_xyz_mean./s1b_xyz_mean;
     s1ab_xyz_mean_err=sqrt( (s1b_xyz_mean_err.*s1a_xyz_mean./(s1b_xyz_mean.^2)).^2 + (s1a_xyz_mean_err./s1b_xyz_mean).^2);

end
 
save('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p21\KrypCal2p21_QualityChecks\EnergyResolution\MappingParametersFit_Start_Oct2015DD_CutDown','s1ab_xyz_xbins','s1ab_xyz_ybins','s1ab_xyz_zbins','s1ab_xyz_mean','Kr_events','S1_xbin_max','S1_xbin_min','S1_ybin_max','S1_ybin_min','S1xybinsize','S2_xbin_max','S2_xbin_min','S2_ybin_max','S2_ybin_min','S2xybinsize','det_edge','drift_time','g','grid_size','i','max_num_bins','min_evts_per_bin','rcut_max','rcut_min','s1_phe_both','s1_phe_bottom','s1_width','s1_xbins','s1_ybins','s1a_P','s1a_S','s1a_fit','s1a_over_bz_plot','s1a_xy_fig','s1a_xy_mean','s1a_xy_mean_err','s1a_z_fit','s1a_z_means','s1a_z_means_err','s1ab_P','s1ab_S','s1ab_bincenters','s1ab_x','s1ab_xbin_max','s1ab_xbin_min','s1ab_xbins','s1ab_xy_fig','s1ab_xy_mean','s1ab_xy_mean_err','s1ab_xybinsize','s1ab_y','s1ab_ybin_max','s1ab_ybin_min','s1ab_ybins','s1ab_z','s1ab_z_means','s1ab_z_means_err','s1area_bound_max','s1area_bound_min','s1az_plot','s1b_P','s1b_S','s1b_fit','s1b_xy_fig','s1b_xy_mean','s1b_xy_mean_err','s1b_z_fit','s1b_z_means','s1b_z_means_err','s1bz_plot','s2_phe_both','s2_phe_bottom','s2_width','s2_xbins','s2_ybins','s2area_bound_max','s2area_bound_min','s2radius','s2x','s2y','x_bin','x_center','x_count','x_max','x_min','y_bin','y_center','y_count','y_max','y_min','z_center','z_max','zcut_max','zcut_min');