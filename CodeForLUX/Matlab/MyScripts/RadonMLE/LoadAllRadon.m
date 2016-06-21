rqs_to_load = {'pulse_area_phe','event_timestamp_samples'...
   ,'pulse_classification' ...
   ,'z_drift_samples' , 's1s2_pairing'...
   ,'top_bottom_ratio','x_cm','y_cm','x_cm_tmplt','y_cm_tmplt'...
   ,'full_evt_area_phe',...
   'event_number','chi2','prompt_fraction','aft_t1_samples','pulse_start_samples',...
   'pulse_end_samples','top_bottom_asymmetry','aft_t0_samples','aft_t2_samples',...
   'full_evt_area_phe','admin',...
   'hft_t50r_samples','hft_t50l_samples','z_corrected_pulse_area_all_phe','xyz_corrected_pulse_area_all_phe',...
   'z_corrected_pulse_area_bot_phe','xyz_corrected_pulse_area_bot_phe','x_corrected','y_corrected'};


% defining cuts
zcut_min = 30;%us
zcut_max = 300;%us
rcut_min = 0;
rcut_max = 25;%cm
s1area_bound_min = 0;%2 PMT hits
s1area_bound_max = 10^5;%include Kr 83 data and gammas

s2area_bound_min = 100;%100 ~ 0.5 keVee
s2area_bound_max = 2*10^6;% %both PMT Arrays and Kr83


folders_to_load=dir('C:\Users\Richard\Desktop\Run03Radon\');
folders_to_load(1:2)=[];
   
    
ii=1;
path=strcat('C:\Users\Richard\Desktop\Run03Radon\',folders_to_load(ii).name);
d = LUXLoadMultipleRQMs_framework(path,rqs_to_load);


d.z_drift_samples(isnan(d.z_drift_samples)) = 0.0; % get rid of NaN        

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

drift_time = d.z_drift_samples(s2_single_cut)/100; %units of us


d.phe_bottom=d.pulse_area_phe./(1+d.top_bottom_ratio); %bottom PMT pulse area

s1_phe_both = d.pulse_area_phe(s1_single_cut);
s1_phe_bottom = d.phe_bottom(s1_single_cut);
s1_phe_both_z=d.z_corrected_pulse_area_all_phe(s1_single_cut);
s1_phe_both_xyz=d.xyz_corrected_pulse_area_all_phe(s1_single_cut);

s2_phe_both = d.pulse_area_phe(s2_single_cut);
s2_phe_bottom = d.phe_bottom(s2_single_cut);
s2_phe_bottom_z=d.z_corrected_pulse_area_bot_phe(s2_single_cut); 
s2_phe_bottom_xyz=d.xyz_corrected_pulse_area_bot_phe(s2_single_cut); 

s2x = d.x_cm(s2_single_cut);
s2y = d.y_cm(s2_single_cut);
s2radius = (s2x.^2+s2y.^2).^(0.5);
s2x_c = d.x_corrected(s2_single_cut);
s2y_c = d.y_corrected(s2_single_cut);
s2radius_c = (s2x_c.^2+s2y_c.^2).^(0.5);


timestamp_vec_1 = sort([d.livetime_latch_samples d.livetime_end_samples]);
livetime_sec = sum(timestamp_vec_1(2:2:end) - timestamp_vec_1(1:2:end)) / 1e8;

evt_cut=logical(sum(s2_single_cut));%Cut for all the events passing the single S1 & S2 cut
event_number=d.event_number(evt_cut)';
event_timestamp_samples=d.event_timestamp_samples(evt_cut)';



clean_cut=-(s1_phe_both+s2_phe_both)' + d.full_evt_area_phe(logical(squeeze(sum(s2_single_cut,1)))) < 100 ...
& drift_time' > 0;%more than 1/2 the area is S1+S2. And S1 before S2
clean_cut=clean_cut';


s1_phe_both_2=s1_phe_both;
s1_phe_both_z_2=s1_phe_both_z ;
s1_phe_both_xyz_2=s1_phe_both_xyz ;        
s2_phe_bottom_2=s2_phe_bottom ;
s2_phe_bottom_xyz_2= s2_phe_bottom_xyz ;
s2_phe_bottom_z_2= s2_phe_bottom_z ;
s2_phe_both_2=s2_phe_both ;
drift_time_2=drift_time ;
event_timestamp_samples_2=event_timestamp_samples ;
event_number_2=event_number;
s2x_2=s2x ;
s2y_2=s2y ;
s2x_c_2 = s2x_c ;
s2y_c_2 = s2y_c ;
clean_cut_2=clean_cut;

        
time_diff_samples_1=0; %no time difference from first file
file_start_samples(ii)=time_diff_samples_1;
file_end_samples(ii)=time_diff_samples_1+d.event_timestamp_samples(end);

for  ii=2:length(folders_to_load) % 
    tic; 
    path=strcat('C:\Users\Richard\Desktop\Run03Radon\',folders_to_load(ii).name);
    clear d1;
    d1=d;
    clear d;

    d = LUXLoadMultipleRQMs_framework(path,rqs_to_load);


    d.z_drift_samples(isnan(d.z_drift_samples)) = 0.0; % get rid of NaN        

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

    drift_time = d.z_drift_samples(s2_single_cut)/100; %units of us

    d.phe_bottom=d.pulse_area_phe./(1+d.top_bottom_ratio); %bottom PMT pulse area

    s1_phe_both = d.pulse_area_phe(s1_single_cut);
    s1_phe_bottom = d.phe_bottom(s1_single_cut);
    s1_phe_both_z=d.z_corrected_pulse_area_all_phe(s1_single_cut);
    s1_phe_both_xyz=d.xyz_corrected_pulse_area_all_phe(s1_single_cut);

    s2_phe_both = d.pulse_area_phe(s2_single_cut);
    s2_phe_bottom = d.phe_bottom(s2_single_cut);
    s2_phe_bottom_z=d.z_corrected_pulse_area_bot_phe(s2_single_cut); 
    s2_phe_bottom_xyz=d.xyz_corrected_pulse_area_bot_phe(s2_single_cut); 

    s2x = d.x_cm(s2_single_cut);
    s2y = d.y_cm(s2_single_cut);
    s2radius = (s2x.^2+s2y.^2).^(0.5);
    s2x_c = d.x_corrected(s2_single_cut);
    s2y_c = d.y_corrected(s2_single_cut);
    s2radius_c = (s2x_c.^2+s2y_c.^2).^(0.5);

    timestamp_vec_1 = sort([d.livetime_latch_samples d.livetime_end_samples]);
    livetime_sec = livetime_sec + sum(timestamp_vec_1(2:2:end) - timestamp_vec_1(1:2:end)) / 1e8;

    evt_cut=logical(sum(s2_single_cut));%Cut for all the events passing the single S1 & S2 cut
    event_number=d.event_number(evt_cut)';
    event_timestamp_samples=d.event_timestamp_samples(evt_cut)';

    clean_cut=-(s1_phe_both+s2_phe_both)' + d.full_evt_area_phe(logical(squeeze(sum(s2_single_cut,1)))) < 100 ...
    & drift_time' > 0;%more than 1/2 the area is S1+S2. And S1 before S2
    clean_cut=clean_cut';


    s1_phe_both=[s1_phe_both ; s1_phe_both_2]; %#ok<*AGROW>
    s1_phe_both_z=[s1_phe_both_z ; s1_phe_both_z_2];
    s1_phe_both_xyz=[s1_phe_both_xyz ; s1_phe_both_xyz_2] ;        
    s2_phe_bottom=[s2_phe_bottom ; s2_phe_bottom_2] ;
    s2_phe_bottom_xyz= [s2_phe_bottom_xyz ; s2_phe_bottom_xyz_2];
    s2_phe_bottom_z= [s2_phe_bottom_z ; s2_phe_bottom_z_2];
    s2_phe_both=[s2_phe_both ; s2_phe_both_2];
    drift_time=[drift_time ; drift_time_2];
    event_timestamp_samples=[event_timestamp_samples_2 ; event_timestamp_samples+time_diff_samples_1;];
    event_number=[event_number; event_number_2];
    s2x=[s2x ; s2x_2];
    s2y=[s2y ; s2y_2];
    s2radius = (s2x.^2+s2y.^2).^(0.5);
    s2x_c = [s2x_c ; s2x_c_2];
    s2y_c = [s2y_c ; s2y_c_2];
    s2radius_c = (s2x_c.^2+s2y_c.^2).^(0.5);
    clean_cut=[clean_cut ; clean_cut_2 ];

    s1_phe_both_2=s1_phe_both;
    s1_phe_both_z_2=s1_phe_both_z ;
    s1_phe_both_xyz_2=s1_phe_both_xyz ;        
    s2_phe_bottom_2=s2_phe_bottom ;
    s2_phe_bottom_xyz_2= s2_phe_bottom_xyz ;
    s2_phe_bottom_z_2= s2_phe_bottom_z ;
    s2_phe_both_2=s2_phe_both ;
    drift_time_2=drift_time ;
    event_timestamp_samples_2=event_timestamp_samples ;
    event_number_2=event_number;
    s2x_2=s2x ;
    s2y_2=s2y ;
    s2x_c_2 = s2x_c ;
    s2y_c_2 = s2y_c ;
    clean_cut_2=clean_cut;

    save('ALL_alpha_data','drift_time','s2radius','s1_phe_both_xyz','s1_phe_both','s2_phe_bottom_xyz','event_timestamp_samples'...
    ,'s2radius_c','s1_phe_both','s2_phe_bottom','s2_phe_both','s2_phe_bottom','event_number','s2x','s2y','s2x_c','s2y_c','livetime_sec','clean_cut'...
    , 'file_start_samples','file_end_samples');              

    fprintf(strcat(num2str(ii),' out of ', num2str(length(folders_to_load)),' data sets loaded') );

    toc;

end
    

fprintf('Finished');