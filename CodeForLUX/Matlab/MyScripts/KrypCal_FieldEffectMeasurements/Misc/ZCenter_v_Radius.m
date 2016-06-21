path='C:\Users\Richard\Desktop\LUX Analysis\lux10_20150929T1905_cp18197';

   
rqs_to_load = {'pulse_area_phe','event_timestamp_samples'...
   ,'pulse_classification' ...
   ,'z_drift_samples' , 's1s2_pairing'...
   ,'top_bottom_ratio','x_cm','y_cm'...
   ,'full_evt_area_phe',...
   'event_number','chi2','prompt_fraction','aft_t1_samples','pulse_start_samples',...
   'pulse_end_samples','top_bottom_asymmetry','aft_t0_samples','aft_t2_samples',...
   'full_evt_area_phe','admin','Kr83fit_s1a_area_phe','Kr83fit_s1b_area_phe',...
   'hft_t50r_samples','hft_t50l_samples','Kr83fit_dt_samples','x_cm_del','y_cm_del'};

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
    s2x_del=d.x_cm_del(s2_single_cut);
    s2y_del=d.y_cm_del(s2_single_cut);

    d.livetime_sec=sum(d.livetime_end_samples-d.livetime_latch_samples)/1e8;
    evt_cut=logical(sum(s2_single_cut));%Cut for all the events passing the single S1 & S2 cut
    event_timestamp_samples=d.event_timestamp_samples(evt_cut);
    
    time_wait_cut=event_timestamp_samples/1e8/60 > delay_min; %allow x min for Kr mixing
    
    s2_width=(d.aft_t2_samples(s2_single_cut)-d.aft_t0_samples(s2_single_cut)); %cut above 800 samples
    s1_width=d.pulse_end_samples(s1_single_cut)-d.pulse_start_samples(s1_single_cut);
    
    s2radius = (s2x.^2+s2y.^2).^(0.5);
    s2radius_del=(s2x_del.^2+s2y_del.^2).^(0.5);
    clean_cut = (log10(s1_phe_both)+0.6*log10(s2_phe_both)>4.5) & s2radius<25;

         s1_phe_bottom=s1_phe_bottom(clean_cut);
         s1_phe_both=s1_phe_both(clean_cut);
         s2_phe_bottom=s2_phe_bottom(clean_cut);
         s2_phe_both=s2_phe_both(clean_cut);
         drift_time=drift_time(clean_cut);
         s2x=s2x(clean_cut);
         s2y=s2y(clean_cut);
         s2x_del=s2x_del(clean_cut);
         s2y_del=s2y_del(clean_cut);
         s2radius_del=s2radius_del(clean_cut);
         s2radius = (s2x.^2+s2y.^2).^(0.5);
         event_timestamp_samples=event_timestamp_samples(clean_cut);
         s2_width=s2_width(clean_cut);
         s1_width=s1_width(clean_cut);

kr_energy_cut = inrange(s1_phe_both,[s1area_bound_min s1area_bound_max]) & inrange(s2_phe_both,[s2area_bound_min s2area_bound_max]);%for counting Kr events
Kr_events=length(drift_time(kr_energy_cut));% Are there enough events for the XY or XYZ map?

%%%% Detector Center Detection %%%%%%%
x_center=mean(s2x(drift_time>4)); %exclude extraction field region at uSec<4
y_center=mean(s2y(drift_time>4));
% z_center=mean(drift_time(drift_time>4));

for i=1:1:250;
    z_center(i)=mean(drift_time(drift_time>4 & s2radius<=(i/10)));
end

figure
plot([0.1:0.1:25],z_center,'.k')

for i=1:1:250;
    z_center_del(i)=mean(drift_time(drift_time>4 & s2radius_del<=(i/10)));
end

figure
plot([0.1:0.1:25],z_center_del,'.k')
title('Delensed Radius')

for i=1:1:250;
    z_count(i)=length(drift_time(drift_time>4 & s2radius_del.^2<=(i/10).^2));
end


figure
plot([0.1:0.1:25],z_count,'.k')
title('Delensed Radius - Z count')

figure %changing radial cut here shows the efffect of the field effect on radial cut result
step(drift_time(s2radius<20),[0:2:360])



%Nonuniform drift time even in s2radius_del... leads to artificially high drift time z_center
figure
step(drift_time(s2radius_del<15),[0:2:360])

%% What Lucie's map tells us
load('C:\Program Files\MATLAB\R2012a\bin\Run04_ERBand\Kr2p20\LucieSep2015_Map');
true_x= LucieSep2015Map(:,1);
true_y= LucieSep2015Map(:,2);
true_z= LucieSep2015Map(:,3);
x= LucieSep2015Map(:,5); %uncorrected x
y= LucieSep2015Map(:,6); %uncorrected y
field = LucieSep2015Map(:,4)/100; %Field in V/cm instead of V/m
z=LucieSep2015Map(:,7); %drift time in uSec
r=sqrt(x.^2+y.^2);
event_class=LucieSep2015Map(:,8);

%True z center of the detector (in lucie's coordinates) is ~-10cm
 mean(true_z)
 
%The drift time corresponding to this, near xy=0 is 143
mean(z(~isnan(z) & true_radius<2 & inrange(true_z,[-10.5,-9.5])))

%However, if we look at the drift time histogram from Lucie's map, at low radius
%it is clear that the "edge" of the detector is around 330 uSec
figure
step(z(~isnan(z) & true_radius<7),[0:6:360])

%I think the low "average" drift time value comes from the holes in the dT histo at
%higher dT, and the spike around 30 uSec.  This is a problem with Lucie's map, and is not physical.
%going off of the "edge" of the histogram around 330 uSec, the drift time corresponding to the detector
%center is about 160 uSec... which is consistent with doing no radial cut in the data