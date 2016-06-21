path='C:\Users\Richard\Desktop\LUX Analysis\lux10_20150929T1905_cp18197';

   
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
s2area_bound_max = 70000;%30000
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
    s1ab_radius=s1ab_radius(s1ab_timing_cut);
    
    dT_step=det_edge/(length(s1ab_z)/1000); %2000 events per z bin
    
%% load lucie's sep 2015 field map

load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p21\TritiumMethod\LucieSep2015_Map_v1p0');
field_s2x=LucieSep2015_Map_v1p0(:,5);
field_s2y=LucieSep2015_Map_v1p0(:,6);
field_drift_time=LucieSep2015_Map_v1p0(:,7);
field_values=0.01.*LucieSep2015_Map_v1p0(:,4);

field_cut=~isnan(field_s2x) & ~isnan(field_s2y) & ~isnan(field_drift_time) & ~isnan(field_values);

field_s2x=field_s2x(field_cut);
field_s2y=field_s2y(field_cut);
field_drift_time=field_drift_time(field_cut);
field_values=field_values(field_cut);


%% 3D Kr s1a, s1b, and field map

s1ab_xyz_numbins=floor((length(s1ab_z)/250)^(1/3)); %number of bins in one direction
s1ab_xyz_zstep=det_edge/s1ab_xyz_numbins;
s1ab_xyz_xstep=50/s1ab_xyz_numbins;
s1ab_xyz_ystep=50/s1ab_xyz_numbins;
s1ab_xyz_xmax=25;
r_max=25;

s1ab_xyz_zbins=10+s1ab_xyz_zstep/2:s1ab_xyz_zstep:10+s1ab_xyz_zstep/2+s1ab_xyz_zstep*s1ab_xyz_numbins; %20 us
s1ab_xyz_xbins=(-r_max+s1ab_xyz_xstep/2):s1ab_xyz_xstep:r_max;
s1ab_xyz_ybins=(-r_max+s1ab_xyz_ystep/2):s1ab_xyz_ystep:r_max;

s1a_xyz_mean=zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));
s1a_xyz_mean_err=zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));
s1b_xyz_mean=zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));
s1b_xyz_mean_err=zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));

field_xyz_mean=zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));
field_xyz_mean_error=zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));


for k = s1ab_xyz_zbins; %s1zbins are the center of the bin    
    for i = (-r_max+s1ab_xyz_xstep):s1ab_xyz_xstep:r_max %map to rows going down (y_cm) -- s1x bins are the max edge of the bin
        for j = (-r_max+s1ab_xyz_ystep):s1ab_xyz_ystep:r_max %map columns across (x_cm) -- s1y bins are the max edge of the bin
                       
            l=int8(i/s1ab_xyz_xstep+s1ab_xyz_xmax/s1ab_xyz_xstep);
            m=int8(j/s1ab_xyz_ystep+s1ab_xyz_xmax/s1ab_xyz_ystep);
            n=int8((k-(10+s1ab_xyz_zstep/2))/s1ab_xyz_zstep + 1);
            
           %sort, and make the cut. using the variable q
           q = s1ab_x<j & s1ab_x>(j-s1ab_xyz_xstep) & s1ab_y<i & s1ab_y>(i-s1ab_xyz_ystep) & inrange(s1ab_z,[(k-s1ab_xyz_zstep/2) , (k+s1ab_xyz_zstep/2)]); %no 0th element!
                    %swap x,y due to matrix x-row, y-column definition
           q_field = field_s2x<j & field_s2x>(j-s1ab_xyz_xstep) & field_s2y<i & field_s2y>(i-s1ab_xyz_ystep) & inrange(field_drift_time,[(k-s1ab_xyz_zstep/2) , (k+s1ab_xyz_zstep/2)]); %no 0th element!

            %Count the number of events per bin
            Count_S1ab_3D(l,m,n)=length(s1a_phe_both(q));
            
            if length(field_values(q_field)) >=10
                field_xyz_mean(l,m,n)=mean(field_values(q_field));
                field_xyz_mean_error(l,m,n)=std(field_values(q_field));
            else
                field_xyz_mean(l,m,n)=nan;
                field_xyz_mean_error(l,m,n)=nan;                
            end
            
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


% %Removing zeros from the maps with nearest neighbor.
s1ab_xyz_mean_old=s1ab_xyz_mean;
temp_map=s1ab_xyz_mean;
q= ~isnan(temp_map);
r=nearestpoint(find(~q),find(q));
s1ab_xyz_mean=temp_map; temp_mapQ=temp_map(q); s1ab_xyz_mean(~q)=temp_mapQ(r);s1ab_xyz_mean=reshape(s1ab_xyz_mean,size(s1ab_xyz_mean_old,1),size(s1ab_xyz_mean_old,2),size(s1ab_xyz_mean_old,3));

 
s1a_xyz_mean_old=s1a_xyz_mean;
temp_map=s1a_xyz_mean;
q= ~isnan(temp_map);
r=nearestpoint(find(~q),find(q));
s1a_xyz_mean=temp_map; temp_mapQ=temp_map(q); s1a_xyz_mean(~q)=temp_mapQ(r);s1a_xyz_mean=reshape(s1a_xyz_mean,size(s1a_xyz_mean_old,1),size(s1a_xyz_mean_old,2),size(s1a_xyz_mean_old,3));

s1b_xyz_mean_old=s1b_xyz_mean;
temp_map=s1b_xyz_mean;
q= ~isnan(temp_map);
r=nearestpoint(find(~q),find(q));
s1b_xyz_mean=temp_map; temp_mapQ=temp_map(q); s1b_xyz_mean(~q)=temp_mapQ(r);s1b_xyz_mean=reshape(s1b_xyz_mean,size(s1b_xyz_mean_old,1),size(s1b_xyz_mean_old,2),size(s1b_xyz_mean_old,3));

field_xyz_mean_old=field_xyz_mean;
temp_map=field_xyz_mean;
q= ~isnan(temp_map);
r=nearestpoint(find(~q),find(q));
field_xyz_mean=temp_map; temp_mapQ=temp_map(q); field_xyz_mean(~q)=temp_mapQ(r);field_xyz_mean=reshape(field_xyz_mean,size(field_xyz_mean_old,1),size(field_xyz_mean_old,2),size(field_xyz_mean_old,3));

field_xyz_mean_error_old=field_xyz_mean_error;
temp_map=field_xyz_mean_error;
q= ~isnan(temp_map);
r=nearestpoint(find(~q),find(q));
field_xyz_mean_error=temp_map; temp_mapQ=temp_map(q); field_xyz_mean_error(~q)=temp_mapQ(r);field_xyz_mean_error=reshape(field_xyz_mean_error,size(field_xyz_mean_error_old,1),size(field_xyz_mean_error_old,2),size(field_xyz_mean_error_old,3));


%% Estimate the field throughout the detector using S1a/S1b r^2vZ maps

load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal\KrypCal_2p20_Code\2015_Field_to_S1aS1b.mat');
[field_map_2015, field_map_2015_staterr]=polyval(S1aS1b_fit_2015,s1ab_xyz_mean,S1aS1b_fit_2015_error);

load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal\KrypCal_2p20_Code\2014_Field_to_S1aS1b.mat');
[field_map_2014, field_map_2014_staterr]=polyval(S1aS1b_fit_2014,s1ab_xyz_mean,S1aS1b_fit_2014_error);

field_map_lucie=field_xyz_mean;
field_map_lucie_err=field_xyz_mean_error;
field_map_syserr=((field_map_2015-(field_map_2014+field_map_2015+field_map_lucie)/3).^2 + (field_map_2014-(field_map_2014+field_map_2015+field_map_lucie)/3).^2 + (field_map_lucie-(field_map_2014+field_map_2015+field_map_lucie)/3).^2).^(1/2);

%field_map_syserr is ~16.5% 
mean_field_map=field_map_lucie; %(field_map_2014+field_map_2015+field_map_lucie)/3;
mean_field_map_totalerr=sqrt(field_map_syserr.^2+((field_map_2014_staterr.^2+field_map_2015_staterr.^2+field_map_lucie_err.^2)/(3^2)));


%% Use the field map to determine number of S1 photons in XYZ 
%photons per keV = 118.90 * [electric field in V/cm] ^ -0.15626

% num_phot=118.9.*mean_field_map.^-0.15626; %num_phot/keV
num_phot=(-0.0004895.*mean_field_map.*log(1 + 1/((8.9e-4).*mean_field_map)) + 1).*55.2;

%errors 118.9 +/- 4.14, -0.15626 +/- 0.006
% num_phot_stat_err=sqrt( (4.14.*mean_field_map.^-0.15626).^2 + (118.9.*(mean_field_map.^-0.15626).*log(mean_field_map).*0.006).^2 + (118.9.*(-0.15626).*(mean_field_map.^(-0.15626-1)).*mean_field_map_totalerr).^2);
% num_phot_err=sqrt( num_phot_stat_err.^2 + (0.02.*num_phot).^2); %adding 2% systematic error from Matthew

a1=-0.55; del_a1=0.03;
a2=(8.9e-4); del_a2=(1.6e-4);
c=55.2; del_c=1.4904;

dphot_dc=a1.*a2.*mean_field_map.*log(1+1./(a2.*mean_field_map))+1;
dphot_da1=a2.*c.*mean_field_map.*log(1+1./(a2.*mean_field_map));
dphot_da2=a1.*c.*mean_field_map.*( ((mean_field_map.*a2+1).*log(1./(mean_field_map.*a2)+1)-1)./(mean_field_map.*a2+1));
dphot_df=a1.*a2.*c.*( ((mean_field_map.*a2+1).*log(1./(mean_field_map.*a2)+1)-1)./(mean_field_map.*a2+1));

num_phot_err=sqrt( (dphot_dc.*del_c).^2 + (dphot_da1.*del_a1).^2+(dphot_da2.*del_a2).^2+(dphot_df.*mean_field_map_totalerr).^2);

%% We have an estimate of the expected number of photons v XYZ based on field strength alone. FOR 3D Map
% Therefore, num_phot(xyz)/num_phot(center) tells us the strength of the field effect in S1a data
num_phot_center=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,num_phot,x_center,y_center,z_center,'cubic');
num_phot_center_err=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,num_phot_err,x_center,y_center,z_center,'cubic');


S1fieldeffect=num_phot./num_phot_center;
S1fieldeffect_err=sqrt( (num_phot_err./num_phot_center).^2 + ( (num_phot.*num_phot_center_err)/(num_phot_center.^2)).^2 );

s1as1b_at_center=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1ab_xyz_mean,x_center,y_center,z_center,'cubic');

% % Plot S1fieldeffect versus S1a/S1b... hopefully it's similar to what KrypCal measures
% figure
% rkploterr(s1ab_xyz_mean,S1fieldeffect,s1ab_xyz_mean_err,S1fieldeffect_err,[0 0 0],'.',100,1);
% xlabel('S1a/S1b');
% ylabel('Number of S1a photons (xyz) / Number of S1a photons (center)');
% myfigview(16);
% 
% %Overlay Kr2.20 result, renormalized to the new S1aS1b center value
% s1ab_s1_xyz_fit_map.p1=0.34;
% s1ab_s1_xyz_fit_map.p2=-1.4;
% s1ab_s1_xyz_fit_map.p3=2.2874;
% x=[2.2:0.1:3.1]; 
% y=  (s1ab_s1_xyz_fit_map.p1.*x.^2+s1ab_s1_xyz_fit_map.p2.*x+s1ab_s1_xyz_fit_map.p3)./(s1ab_s1_xyz_fit_map.p1.*s1as1b_at_center.^2+s1ab_s1_xyz_fit_map.p2.*s1as1b_at_center+s1ab_s1_xyz_fit_map.p3); 
% plot(x,y,'-b','LineWidth',3);
% 
% 
% %Overlay Kr2.18 Result
% s1ab_s1_xyz_fit_map.p1=0.46;
% b2=1-2.7347.*s1ab_s1_xyz_fit_map.p1;
% s1ab_s1_xyz_fit_map.p2=b2;
% x=[2.2:0.1:3.1]; 
% y= (s1ab_s1_xyz_fit_map.p1.*x+s1ab_s1_xyz_fit_map.p2)./(s1ab_s1_xyz_fit_map.p1.*s1as1b_at_center+s1ab_s1_xyz_fit_map.p2); 
% plot(x,y,'-r','LineWidth',3);
% ylim([0.5 2]);
% xlim([2.2 3.1]);
% 
% % Fit a polynomial to the s1 field effect result
% S1fieldpoly=polyfit(s1ab_xyz_mean(S1fieldeffect~=0 & s1ab_xyz_mean~=0),S1fieldeffect(S1fieldeffect~=0 & s1ab_xyz_mean~=0),3);
% plot(x,polyval(S1fieldpoly,x),'Color',[0.7 0.7 0.7],'LineWidth',3);
% title('These shouldnt necessarily be the same, Grey is S1a(XYZ)/S1(C), and others are (S1(XYZ)/S1(center))')
% %SHOULD SWITCH TO A SPLINE FIT!!


%% estimate g1 from S1a at center and field at center.  Since S1 is normalized to the center, detector inefficiency correction=1 here

s1a_center=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1a_xyz_mean,x_center,y_center,z_center,'cubic');
s1a_center_err=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1a_xyz_mean_err,x_center,y_center,z_center,'cubic');

g1_interp_measurement=s1a_center./(32.1.*num_phot_center);
g1_interp_measurement_err=sqrt( (s1a_center_err./(32.1.*num_phot_center)).^2 + ( (s1a_center.*num_phot_center_err)/(32.1.*num_phot_center.^2)).^2);


%% Now remove the field effect from raw S1a data, and measure the residual detector inefficency variation -- THIS SHOULD BE THE SAME AS THE ONE MEASURED FROM S1_E
% 
%     s1ab_x=s1ab_x(s1ab_timing_cut);
%     s1ab_y=s1ab_y(s1ab_timing_cut);
%     s1ab_z=s1ab_z(s1ab_timing_cut);
%     s1a_phe_both=s1a_phe_both(s1ab_timing_cut);
%     s1b_phe_both=s1b_phe_both(s1ab_timing_cut);
%     s1ab_radius=s1ab_radius(s1ab_timing_cut);
%     



%Normalization of S1a field effect to center: S1fieldeffect(xyz/c) is equivalent to 1/norm_s1_both_xyz in KrypCal
S1fieldeffect_norm=1./S1fieldeffect; %This is really s1a field effect
S1fieldeffect_norm(isinf(S1fieldeffect_norm))=1;%remove infinity. no correction outside 25cm
S1fieldeffect_norm(isnan(S1fieldeffect_norm))=1;%remove nan. no correction outside 25cm
S1fieldeffect_norm_err=S1fieldeffect_err./(S1fieldeffect.^2);
%Removing field effect from raw S1a data

s1a_phe_both_f=s1a_phe_both.*interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,S1fieldeffect_norm,s1ab_x,s1ab_y,s1ab_z,'spline').';
s1a_phe_both_f_err=s1a_phe_both.*interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,S1fieldeffect_norm_err,s1ab_x,s1ab_y,s1ab_z,'cubic').';

%% S1a_f Z dep
   
    clear y_s1_both y_s1_bottom Fit_s1_both Fit_s1_bottom mean_s1_both mean_s1_both
    clear sigma_s1_both sigma_s1_both mean_index means mean_index means_error x temp_hist
    clear hist_s1_both x xfit S
    
    %set up energy bins for the gaussian fit
    bin_s1=8;%5
    s1_bin_min=50;
    s1_bin_max=500;
    x=s1_bin_min:bin_s1:s1_bin_max; %150:bin_s1:500;
    x=x';
        
    cut = inrange(s1ab_radius,[rcut_min,rcut_max]) & inrange(s1a_phe_both_f,[s1_bin_min,s1_bin_max]).';   
    
    A = sort(size(s1_phe_both(cut))); %%sort, incase 1xn or nx1 vector is defined.
    bins = ceil(A(2)./min_evts_per_bin); %get the number of points, adjust binning if stats are low.
    if bins>max_num_bins
       bins = max_num_bins; % usually we will have 65 bin. 0:5:325
    end
%     bin_size = (488/bins);%span the entire detector. ~49.0 cm. Bin centers every 4\mus
     bin_size = (floor(det_edge/bins*4)/4);%span the entire detector. 325us. typically = 5us. Round to nearest 0.25.

   
    hist_s1_both = zeros(length(x),bins);
    Fit_s1_both = cell(1,bins);
    means = zeros(bins,1);
    means_error = zeros(bins,1);
    mean_index = zeros(bins,1);
    sys_error_z=zeros(bins,1);
    
      
    for bin=1:bins

     binStart = ((bin-1).*bin_size);  binEnd = (bin.*bin_size);   

%      time_cut = inrange(drift_distance_\mus,[binStart,binEnd]) & cut==1;
       time_cut = inrange(s1ab_z,[binStart,binEnd]) & cut==1;
     
        hist_s1_both(:,bin)= hist(s1a_phe_both_f(time_cut),x)'/bin_s1;
        temp_hist=hist_s1_both(:,bin);   
            amp_start=max(temp_hist);
            mean_start=sum(temp_hist.*x)/sum(temp_hist);
            sigma_start=std(temp_hist.*x);
        Fit_s1_both{bin}=fit(x(2:end-1),temp_hist(2:end-1),'gauss1');%remove bin edges

     means(bin)=Fit_s1_both{bin}.b1;
     sys_error_z(bin)=mean(s1a_phe_both_f_err(time_cut)); %average systematic error for the bin
     means_error(bin) = Fit_s1_both{bin}.c1/sqrt(2)/sqrt(sum(hist_s1_both(:,bin))*bin_s1); % == sigma/sqrt(N);
     mean_index(bin) = (binStart+binEnd)/2;     

    end
    
 % //////////////////////////////Make Plots/////////////////////////////////
 
    dT_range= inrange( mean_index,[0, 350] );
    dT_fit=0:1:det_edge;    
    [P_s1_both S] = polyfit(mean_index(dT_range),means(dT_range),2);
    sigma_P_s1_both = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df)';    
    
    light_correction_fig = figure;
    errorbar(mean_index,means,means_error,'.k');
    hold on     
    [yplot_s1_both] = polyval(P_s1_both,dT_fit);
    [center_s1_both, center_s1_both_err] = polyval(P_s1_both,z_center, S); %160[us]*1.51[\mus/us]
    plot(dT_fit,yplot_s1_both,'-r','LineWidth',1);    
    xlabel('Depth (\mus)', 'FontSize',16); ylabel('Mean (Pulse Area Phe)', 'FontSize', 16);
    line([0 det_edge],[center_s1_both center_s1_both],'linestyle','--');
    line([z_center; z_center;],[min(means) max(means)],'linestyle','--');
    legend('S1\_both Mean', strcat('y=', num2str(P_s1_both(1),'%10.2e'), '*Z^2 + ', num2str(P_s1_both(2),'%10.2e')...
        , '*Z + ', num2str(P_s1_both(3),4)),strcat('Det center = ',num2str(center_s1_both,4), ' [Phe]'),'location','northwest')
    
    box on;
    myfigview(16);
    xlim([0 det_edge]);
    ylim([0.8*min(means) 1.15*max(means)]);
    

%% Measuring the residual S1a_F XYZ variation to get efficiency corrections - use same binning as S1a/S1b 3D map

s1a_f_xyz_mean=zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));
s1a_f_xyz_mean_err=zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));
s1a_f_xyz_sys_err=zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));

for k = s1ab_xyz_zbins; %s1zbins are the center of the bin    
    for i = (-r_max+s1ab_xyz_xstep):s1ab_xyz_xstep:r_max %map to rows going down (y_cm) -- s1x bins are the max edge of the bin
        for j = (-r_max+s1ab_xyz_ystep):s1ab_xyz_ystep:r_max %map columns across (x_cm) -- s1y bins are the max edge of the bin
                       
            l=int8(i/s1ab_xyz_xstep+s1ab_xyz_xmax/s1ab_xyz_xstep);
            m=int8(j/s1ab_xyz_ystep+s1ab_xyz_xmax/s1ab_xyz_ystep);
            n=int8((k-(10+s1ab_xyz_zstep/2))/s1ab_xyz_zstep + 1);
            
           %sort, and make the cut. using the variable q
           q = s1ab_x<j & s1ab_x>(j-s1ab_xyz_xstep) & s1ab_y<i & s1ab_y>(i-s1ab_xyz_ystep) & inrange(s1ab_z,[(k-s1ab_xyz_zstep/2) , (k+s1ab_xyz_zstep/2)]); %no 0th element!
                    %swap x,y due to matrix x-row, y-column definition

            %Count the number of events per bin
            Count_S1a_f_3D(l,m,n)=length(s1a_phe_both_f(q));

            
            if (Count_S1a_f_3D(l,m,n) >= 100) % at least 100 counts before fitting. 
                s1a_f_z_fit=fit([0:2:400].',hist(s1a_phe_both_f(q.' & inrange(s1a_phe_both_f,[0 400])),[0:2:400]).','gauss1');
                s1a_f_xyz_mean(l,m,n)=s1a_f_z_fit.b1;
                s1a_f_xyz_mean_err(l,m,n)=s1a_f_z_fit.c1/sqrt(2)/length(s1a_phe_both_f(q.' & inrange(s1a_phe_both_f,[0 400])));                
                s1a_f_xyz_sys_err(l,m,n)=mean(s1a_phe_both_f_err(q.' & inrange(s1a_phe_both_f,[0 400]))); %average systematic error for the bin
             else %not enough stats to do the fit
                s1a_f_xyz_mean(l,m,n)=0;
                s1a_f_xyz_mean_err(l,m,n)=0;
                s1a_f_xyz_sys_err(l,m,n)=0;
            end
                                      
        end
    end
k
end

total_xyz_err=sqrt(s1a_f_xyz_sys_err.^2 + s1a_f_xyz_mean_err.^2);

%Producing S1 detector inefficiency correction fro S1a_f map
center_s1a_f=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1a_f_xyz_mean,x_center,y_center,z_center,'spline');%Normalize to the center (x=y=0. dT=160us)
center_s1a_f_staterr=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1a_f_xyz_mean_err,x_center,y_center,z_center,'spline');
center_s1a_f_syserr=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1a_f_xyz_sys_err,x_center,y_center,z_center,'spline');
center_s1a_f_totalerr=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,total_xyz_err,x_center,y_center,z_center,'spline');

norm_s1_both_xyz=center_s1a_f./s1a_f_xyz_mean; %Normalize to center (x=y=0. dT=160)
norm_s1_both_xyz(isinf(norm_s1_both_xyz))=1;%remove infinity. no correction outside 25cm
norm_s1_both_xyz(isnan(norm_s1_both_xyz))=1;%remove nan. no correction outside 25cm
norm_s1_both_xyz_staterr=sqrt( (center_s1a_f_staterr./s1a_f_xyz_mean).^2 + ((center_s1a_f.*s1a_f_xyz_mean_err)./(s1a_f_xyz_mean.^2)).^2);
norm_s1_both_xyz_totalerr=sqrt( (center_s1a_f_totalerr./s1a_f_xyz_mean).^2 + ((center_s1a_f.*total_xyz_err)./(s1a_f_xyz_mean.^2)).^2);

%Fill in NaN errors with nearest neighbors
norm_s1_both_xyz_staterr_old=norm_s1_both_xyz_staterr;
temp_map=norm_s1_both_xyz_staterr;
q= ~isnan(norm_s1_both_xyz_staterr);
r=nearestpoint(find(~q),find(q));
norm_s1_both_xyz_staterr=temp_map; temp_mapQ=temp_map(q); norm_s1_both_xyz_staterr(~q)=temp_mapQ(r);norm_s1_both_xyz_staterr=reshape(norm_s1_both_xyz_staterr,size(norm_s1_both_xyz_staterr_old,1),size(norm_s1_both_xyz_staterr_old,2),size(norm_s1_both_xyz_staterr_old,3));

norm_s1_both_xyz_totalerr_old=norm_s1_both_xyz_totalerr;
temp_map=norm_s1_both_xyz_totalerr;
q= ~isnan(norm_s1_both_xyz_totalerr);
r=nearestpoint(find(~q),find(q));
norm_s1_both_xyz_totalerr=temp_map; temp_mapQ=temp_map(q); norm_s1_both_xyz_totalerr(~q)=temp_mapQ(r);norm_s1_both_xyz_totalerr=reshape(norm_s1_both_xyz_totalerr,size(norm_s1_both_xyz_totalerr_old,1),size(norm_s1_both_xyz_totalerr_old,2),size(norm_s1_both_xyz_totalerr_old,3));
 


%% Now remove the detector inefficiency variation from the raw S1a data and determine g1 using total g1=S1a_E/(num_phot_center*keV).  CROSS CHECK OF KRYPCAL G1
s1a_phe_both_e=s1a_phe_both.*interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,norm_s1_both_xyz,s1ab_x,s1ab_y,s1ab_z,'spline').'; %this is s1a_E
s1a_phe_both_e_totalerr=s1a_phe_both.*interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,norm_s1_both_xyz_totalerr,s1ab_x,s1ab_y,s1ab_z,'cubic').'; %this is s1a_E
s1a_phe_both_e_staterr=s1a_phe_both.*interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,norm_s1_both_xyz_staterr,s1ab_x,s1ab_y,s1ab_z,'cubic').'; %this is s1a_E

s1a_e_xyz_mean=zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));
s1a_e_xyz_mean_err=zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));

for k = s1ab_xyz_zbins; %s1zbins are the center of the bin    
    for i = (-r_max+s1ab_xyz_xstep):s1ab_xyz_xstep:r_max %map to rows going down (y_cm) -- s1x bins are the max edge of the bin
        for j = (-r_max+s1ab_xyz_ystep):s1ab_xyz_ystep:r_max %map columns across (x_cm) -- s1y bins are the max edge of the bin
                       
            l=int8(i/s1ab_xyz_xstep+s1ab_xyz_xmax/s1ab_xyz_xstep);
            m=int8(j/s1ab_xyz_ystep+s1ab_xyz_xmax/s1ab_xyz_ystep);
            n=int8((k-(10+s1ab_xyz_zstep/2))/s1ab_xyz_zstep + 1);
            
           %sort, and make the cut. using the variable q
           q = s1ab_x<j & s1ab_x>(j-s1ab_xyz_xstep) & s1ab_y<i & s1ab_y>(i-s1ab_xyz_ystep) & inrange(s1ab_z,[(k-s1ab_xyz_zstep/2) , (k+s1ab_xyz_zstep/2)]); %no 0th element!
                    %swap x,y due to matrix x-row, y-column definition

            %Count the number of events per bin
            Count_S1a_e_3D(l,m,n)=length(s1a_phe_both_e(q));

            
            if (Count_S1a_e_3D(l,m,n) >= 100) % at least 100 counts before fitting. 
                s1a_e_z_fit=fit([0:2:400].',hist(s1a_phe_both_e(q.' & inrange(s1a_phe_both_e,[0 400])),[0:2:400]).','gauss1');
                s1a_e_xyz_mean(l,m,n)=s1a_e_z_fit.b1;
                s1a_e_xyz_mean_err(l,m,n)=s1a_e_z_fit.c1/sqrt(2)/length(s1a_phe_both_e(q.' & inrange(s1a_phe_both_e,[0 400])));                
                s1a_e_xyz_stat_err(l,m,n)=mean(s1a_phe_both_e_staterr(q.'));
                s1a_e_xyz_total_err(l,m,n)=mean(s1a_phe_both_e_totalerr(q.'));
             else %not enough stats to do the fit
                s1a_e_xyz_mean(l,m,n)=0;
                s1a_e_xyz_mean_err(l,m,n)=0;
                s1a_e_xyz_stat_err(l,m,n)=0;
                s1a_e_xyz_total_err(l,m,n)=0;
            end
                                      
        end
    end
k
end

%Producing S1 detector inefficiency correction fro S1a_f map
center_s1a_e=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1a_e_xyz_mean,x_center,y_center,z_center,'spline');%Normalize to the center (x=y=0. dT=160us)
center_s1a_e_err=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1a_e_xyz_mean_err,x_center,y_center,z_center,'spline');%Normalize to the center (x=y=0. dT=160us)
center_s1a_e_statsyserr=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1a_e_xyz_stat_err,x_center,y_center,z_center,'spline');
center_s1a_e_totalstaterr=sqrt(center_s1a_e_statsyserr.^2 + center_s1a_e_err.^2);
center_s1a_e_totalerr=sqrt(center_s1a_e_err.^2  + interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1a_e_xyz_total_err,x_center,y_center,z_center,'spline').^2);


%total signal fits, to see if it is consistent with center_s1a_e
% s1a_E_fit=fit([0:1:400].',hist(s1a_phe_both_e,[0:1:400]).','gauss1');
% s1a_E_fit_err=confint(s1a_E_fit,0.68);
% s1a_E_mean=s1a_E_fit.b1;
% s1a_E_mean_err=abs(s1a_E_fit.b1-s1a_E_fit_err(1,2));

% g1_measurement=s1a_E_mean./(32.1.*num_phot_center)
% g1_measurement_err=sqrt( (s1a_E_mean_err./(32.1.*num_phot_center)).^2 + ( (s1a_E_mean.*num_phot_center_err)/(32.1.*num_phot_center.^2)).^2);

g1_measurement=center_s1a_e./(32.1.*num_phot_center)
g1_measurement_err=sqrt( (center_s1a_e_totalerr./(32.1.*num_phot_center)).^2 + ( (center_s1a_e.*num_phot_center_err)/(32.1.*num_phot_center.^2)).^2); %was just _err before, not _totalerr

%% Now remove the detector inefficency variation from the raw S1 (not S1a) data, and determine the field effect in S1 data.  Compare to KrypCal results. CROSS CHECK OF KRYPCAL S1 RELATIONSHIP

s1_phe_both_xyz=s1_phe_both.*interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,norm_s1_both_xyz,s2x,s2y,drift_time,'spline');
s1_phe_both_xyz_staterr=s1_phe_both.*interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,norm_s1_both_xyz_staterr,s2x,s2y,drift_time,'spline',mean(norm_s1_both_xyz_staterr(:)));
s1_phe_both_xyz_totalerr=s1_phe_both.*interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,norm_s1_both_xyz_totalerr,s2x,s2y,drift_time,'spline',mean(norm_s1_both_xyz_totalerr(:)));

s1_phe_both_xyz_mean=zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));
s1_phe_both_xyz_mean_err=zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));

for k = s1ab_xyz_zbins; %s1zbins are the center of the bin    
    for i = (-r_max+s1ab_xyz_xstep):s1ab_xyz_xstep:r_max %map to rows going down (y_cm) -- s1x bins are the max edge of the bin
        for j = (-r_max+s1ab_xyz_ystep):s1ab_xyz_ystep:r_max %map columns across (x_cm) -- s1y bins are the max edge of the bin
                       
            l=int8(i/s1ab_xyz_xstep+s1ab_xyz_xmax/s1ab_xyz_xstep);
            m=int8(j/s1ab_xyz_ystep+s1ab_xyz_xmax/s1ab_xyz_ystep);
            n=int8((k-(10+s1ab_xyz_zstep/2))/s1ab_xyz_zstep + 1);
            
           %sort, and make the cut. using the variable q
           q = s2x<j & s2x>(j-s1ab_xyz_xstep) & s2y<i & s2y>(i-s1ab_xyz_ystep) & inrange(drift_time,[(k-s1ab_xyz_zstep/2) , (k+s1ab_xyz_zstep/2)]); %no 0th element!
                    %swap x,y due to matrix x-row, y-column definition

            %Count the number of events per bin
            Count_3D(l,m,n)=length(s1_phe_both_xyz(q));

            
            if (Count_3D(l,m,n) >= 100) % at least 100 counts before fitting. 
                s1_phe_both_xyz_fit=fit([0:2:400].',hist(s1_phe_both_xyz(q & inrange(s1_phe_both_xyz,[0 400])),[0:2:400]).','gauss1');
                s1_phe_both_xyz_mean(l,m,n)=s1_phe_both_xyz_fit.b1;
                s1_phe_both_xyz_mean_err(l,m,n)=s1_phe_both_xyz_fit.c1/sqrt(2)/length(s1_phe_both_xyz(q & inrange(s1_phe_both_xyz,[0 400])));                
                s1_phe_both_xyz_stat_err(l,m,n)=mean(s1_phe_both_xyz_staterr(q));
                s1_phe_both_xyz_total_err(l,m,n)=mean(s1_phe_both_xyz_totalerr(q));
                
             else %not enough stats to do the fit
                s1_phe_both_xyz_mean(l,m,n)=0;
                s1_phe_both_xyz_mean_err(l,m,n)=0;
                s1_phe_both_xyz_stat_err(l,m,n)=0;
                s1_phe_both_xyz_total_err(l,m,n)=0;
                
            end
                                      
        end
    end
k
end

s1xyz_total_stat_err=sqrt(s1_phe_both_xyz_mean_err.^2 + s1_phe_both_xyz_stat_err.^2);
s1xyz_total_err=sqrt( s1_phe_both_xyz_mean_err.^2 + s1_phe_both_xyz_total_err.^2);

%CHANGE THIS TO BE BASED OFF TOTAL FROM ABOVE!?
%Producing S1 detector inefficiency correction fro S1a_f map - ERROR DOES NOT TAKE INTO ACCOUNT THE UNCERTAINITY ON THE S1 EFFICIENCY CORRECTION!
center_s1E_phe_both=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1_phe_both_xyz_mean,x_center,y_center,z_center,'spline');
center_s1E_phe_both_err=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1_phe_both_xyz_mean_err,x_center,y_center,z_center,'spline');
center_s1E_phe_both_statsyserr=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1_phe_both_xyz_stat_err,x_center,y_center,z_center,'spline');
center_s1E_phe_both_totalstaterr=sqrt(center_s1E_phe_both_err.^2 + center_s1E_phe_both_statsyserr.^2);
center_s1E_phe_both_totalerr=sqrt(center_s1E_phe_both_err.^2 + interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1_phe_both_xyz_total_err,x_center,y_center,z_center,'spline').^2);

CombinedS1FieldEffect=s1_phe_both_xyz_mean./center_s1E_phe_both; 
CombinedS1FieldEffect_err=sqrt ( (s1_phe_both_xyz_mean_err./center_s1E_phe_both).^2 + ( (s1_phe_both_xyz_mean.*center_s1E_phe_both_err)./(center_s1E_phe_both.^2)).^2);
CombinedS1FieldEffect_totalstaterr=sqrt ( (s1xyz_total_stat_err./center_s1E_phe_both).^2 + ( (s1_phe_both_xyz_mean.*center_s1E_phe_both_totalstaterr)./(center_s1E_phe_both.^2)).^2);
CombinedS1FieldEffect_totalerr=sqrt ( (s1xyz_total_err./center_s1E_phe_both).^2 + ( (s1_phe_both_xyz_mean.*center_s1E_phe_both_totalerr)./(center_s1E_phe_both.^2)).^2);

%Make S1a/S1b v S1 field effect plot
figure
rkploterr(s1ab_xyz_mean,CombinedS1FieldEffect,s1ab_xyz_mean_err,CombinedS1FieldEffect_totalerr,[0 0 0],'.',100,1);
xlabel('S1a/S1b');
ylabel('S1_E(XYZ)/S1_E(C)');
myfigview(16);

%Overlay Kr2.20 result, renormalized to the new S1aS1b center value
s1ab_s1_xyz_fit_map.p1=0.34;
s1ab_s1_xyz_fit_map.p2=-1.4;
s1ab_s1_xyz_fit_map.p3=2.2874;
x=[2.2:0.1:3.1]; 
y=  (s1ab_s1_xyz_fit_map.p1.*x.^2+s1ab_s1_xyz_fit_map.p2.*x+s1ab_s1_xyz_fit_map.p3)./(s1ab_s1_xyz_fit_map.p1.*s1as1b_at_center.^2+s1ab_s1_xyz_fit_map.p2.*s1as1b_at_center+s1ab_s1_xyz_fit_map.p3); 
plot(x,y,'-b','LineWidth',3);


%Overlay Kr2.18 Result
s1ab_s1_xyz_fit_map.p1=0.46;
b2=1-2.7347.*s1ab_s1_xyz_fit_map.p1;
s1ab_s1_xyz_fit_map.p2=b2;
x=[2.2:0.1:3.1]; 
y= (s1ab_s1_xyz_fit_map.p1.*x+s1ab_s1_xyz_fit_map.p2)./(s1ab_s1_xyz_fit_map.p1.*s1as1b_at_center+s1ab_s1_xyz_fit_map.p2); 
plot(x,y,'-r','LineWidth',3);
ylim([0.5 1.6]);
xlim([2.2 3.1]);

% Fit a polynomial to the s1 field effect result
S1fieldpoly=polyfit(s1ab_xyz_mean(CombinedS1FieldEffect~=0 & s1ab_xyz_mean~=0),CombinedS1FieldEffect(CombinedS1FieldEffect~=0 & s1ab_xyz_mean~=0),2);
plot(x,polyval(S1fieldpoly,x),'Color',[0.7 0.7 0.7],'LineWidth',3);
%SHOULD SWITCH TO A SPLINE FIT!!


%% Use recomb physics to determine S2 relation from this direct measurement
alpha=0.11;
Rkr_center=0.7708; %from S2 measurement code
s2_field_effect=( 1 - ( alpha + Rkr_center)*polyval(S1fieldpoly,[2:0.01:3]) + alpha)/(1 - Rkr_center);

S2fieldpoly=polyfit([2:0.01:3],s2_field_effect,2);

figure
x=[2.2:0.01:3];
y=polyval(S1fieldpoly,x);
plot(x,y,'-r','LineWidth',2);
hold on;
y=polyval(S2fieldpoly,x);
plot(x,y,'-k','LineWidth',2);
xlabel('S1a/S1b'); 
ylabel('S1_E(XYZ)/S1_E(C) or S2_E(XYZ)/S2_E(C)');
legend('S1','S2');
myfigview(16);

save('S1FieldPolynomial_LucieMaps','S1fieldpoly','S2fieldpoly');
