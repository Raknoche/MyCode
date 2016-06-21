%originally this used 2.13 cp17*** to define SE size.  Very small difference, but not significant.

path='C:\Users\Richard\Desktop\LUX Analysis\lux10_20150929T1905_cp18197'; %this is 2.18 corrected

   
rqs_to_load = {'pulse_area_phe','event_timestamp_samples'...
   ,'pulse_classification' ...
   ,'z_drift_samples' , 's1s2_pairing'...
   ,'top_bottom_ratio','x_cm','y_cm'...
   ,'full_evt_area_phe',...
   'event_number','chi2','prompt_fraction','aft_t1_samples','pulse_start_samples',...
   'pulse_end_samples','top_bottom_asymmetry','aft_t0_samples','aft_t2_samples',...
   'full_evt_area_phe','admin','Kr83fit_s1a_area_phe','Kr83fit_s1b_area_phe',...
   'hft_t50r_samples','hft_t50l_samples','Kr83fit_dt_samples','xyz_corrected_pulse_area_all_phe',...
   'correction_s1_z_dependence','correction_s1_xyz_dependence','correction_s2_xy_dependence','correction_s1_z_dependence_bot',...
   'correction_s1_xy_dependence','correction_s1_xy_dependence_bot','correction_s2_xy_dependence_bot','z_corrected_pulse_area_all_phe'};

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
s2_c=d.xyz_corrected_pulse_area_all_phe(s2_single_cut);
s1_c=d.xyz_corrected_pulse_area_all_phe(s1_single_cut);

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
         s2_c=s2_c(clean_cut);
         s1_c=s1_c(clean_cut);

kr_energy_cut = inrange(s1_phe_both,[s1area_bound_min s1area_bound_max]) & inrange(s2_phe_both,[s2area_bound_min s2area_bound_max]);%for counting Kr events
Kr_events=length(drift_time(kr_energy_cut));% Are there enough events for the XY or XYZ map?

%%%% Detector Center Detection %%%%%%%
x_center=mean(s2x(drift_time>4)); %exclude extraction field region at uSec<4
y_center=mean(s2y(drift_time>4));
z_center=mean(drift_time(drift_time>4));
det_edge=330;

z_center
%%

 %   Calculate Single Electron Size 
    s4_class=(d.pulse_classification==4) ; %Use or class=2 to fix 33 phe cut off... have to edit se_wrt_first_s1_func tho
    d = se_wrt_first_s1_func(d);
    cut_se_s1_z = d.t_btw_se_first_s1 > 100*10;

    SE_good_cut =logical( cut_se_s1_z & s4_class & ( cumsum(s1_single_cut)-cumsum(s2_single_cut) ) ); % only events between good S1 and S2s    

    SE_phe_both=d.pulse_area_phe(SE_good_cut);
    SE_phe_bot=d.phe_bottom(SE_good_cut);
    SE_phe_top=SE_phe_both-SE_phe_bot;
    
    SE_x=d.x_cm(SE_good_cut);
    SE_y=d.y_cm(SE_good_cut);
    SE_radius=(SE_x.^2+SE_y.^2).^(1/2);
    
    %Both PMT first fit
    SE_fit=fit([0:1:60]',hist(SE_phe_both(SE_radius<17 & inrange(SE_phe_both,[0 60])),[0:1:60])','gauss1');
    SE_Size=SE_fit.b1;
    SE_sig=SE_fit.c1/sqrt(2);
    SE_both_conf=confint(SE_fit,0.683);
    SE_sig_mu_both=abs(SE_Size-SE_both_conf(1,2)); % one sigma
    SE_sig_sig_both=abs(SE_sig-SE_both_conf(1,3)/sqrt(2)); % one sigma
    
    %Both PMT Second fit
    n_sig_fit=2;
    
    fit_start=SE_Size-n_sig_fit*SE_sig;
    fit_end=SE_Size+n_sig_fit*SE_sig;
    xfit=fit_start:0.1:fit_end;
    
    SE_cut=SE_phe_both(SE_radius<17);
    [xxo, yyo, xo, no] = step(SE_cut(inrange(SE_cut,[fit_start fit_end])),[0:0.2:60]);
    s = fitoptions('Method','NonlinearLeastSquares','Robust','On','Algorithm','Trust-Region','Lower',[-5,0,0,0],'Upper',[5,Inf,Inf,Inf],...
        'Startpoint',[1,max(no),10,5]);
    skewGauss1 = fittype('a1/c1*exp(-((x-b1)/c1)^2)*(1+erf(a*(x-b1)/c1))','options',s);

    [f_both,g_both] = fit(xo',no',skewGauss1);
    delta_both = f_both.a/sqrt(1+f_both.a^2);
    SE_mu_both = f_both.b1 + f_both.c1/sqrt(2)*delta_both*sqrt(2/pi);
    SE_sig_both = sqrt((f_both.c1/sqrt(2))^2 * (1 - 2 * (delta_both)^2/pi));
    skew_both = (4-pi)/2 * (delta_both * sqrt(2/pi))^3 / (1 - 2 * delta_both^2/pi)^(3/2);
    kurt_both = 2*(pi-3) * (delta_both * sqrt(2/pi))^3 / (1 - 2 * delta_both^2/pi)^2;
    ci_both = confint(f_both,0.683);
    y_fit_both = f_both(xfit);   
   
    both_SE_fig=figure;
    hold on;
    step(SE_cut,[0:0.2:60],'k');  
    h_fit=plot(xfit,y_fit_both,'-r','linewidth',2);
    box on;
    myfigview(18);
    xlabel('Pulse Area Phe (Phe) '); ylabel('Count/Phe ');
    legend([h_fit],strcat('\mu = ',num2str(SE_mu_both,'%2.2f'),'\newline \sigma= ',num2str(SE_sig_both,'%2.2f'),'\newline skew=', num2str(skew_both,'%2.2f') ), 'location','northeast');

    
    %% SE Both XY Map
    SExybinsize=sqrt((50*50)/(length(SE_phe_both)/100));
    num_SE_events=length(SE_phe_both);
   se_xbins=SE_xbin_min+SExybinsize/2:SExybinsize:SE_xbin_max;
   se_ybins=SE_ybin_min+SExybinsize/2:SExybinsize:SE_ybin_max;
   
   mean_SE_both_xy = zeros(length(se_xbins),length(se_ybins));
   mean_err_SE_both_xy = zeros(length(se_xbins),length(se_ybins));
   skew_SE_both_xy = zeros(length(se_xbins),length(se_ybins));
   sigma_SE_both_xy = zeros(length(se_xbins),length(se_ybins));
   kurt_SE_both_xy = zeros(length(se_xbins),length(se_ybins));
   sigma_err_SE_both_xy = zeros(length(se_xbins),length(se_ybins));
   
  for x_bin=SE_xbin_min:SExybinsize:(SE_xbin_max-SExybinsize);
    for y_bin=SE_ybin_min:SExybinsize:(SE_ybin_max-SExybinsize);          
      x_min = x_bin; x_max = x_bin+SExybinsize;
      y_min = y_bin; y_max = y_bin+SExybinsize;      
      
        x_count=int32(1+x_bin/SExybinsize+SE_xbin_max/SExybinsize);
        y_count=int32(1+y_bin/SExybinsize+SE_ybin_max/SExybinsize); 
      
      
      bin_cut = (SE_phe_both>0) & inrange(SE_x,[x_min,x_max]) & inrange(SE_y,[y_min,y_max]);            

    if length(SE_phe_both(bin_cut)) > 60; 
        %Both PMT first fit
        clear temp_SE
        temp_SE=SE_phe_both(bin_cut);
        SE_fit_xy=fit([0:1:60]',hist(temp_SE(inrange(temp_SE,[0 60])),[0:1:60])','gauss1');
        SE_Size_xy=SE_fit_xy.b1;
        SE_sig_xy=SE_fit_xy.c1/sqrt(2);
        SE_both_conf_xy=confint(SE_fit_xy,0.683);
        SE_sig_mu_both_xy=abs(SE_Size_xy-SE_both_conf_xy(1,2)); % one sigma
        SE_sig_sig_both_xy=abs(SE_sig_xy-SE_both_conf_xy(1,3)/sqrt(2)); % one sigma
    
        %Both PMT Second fit
        n_sig_fit=2;
    
        fit_start=SE_Size_xy-n_sig_fit*SE_sig_xy;
        fit_end=SE_Size_xy+n_sig_fit*SE_sig_xy;
        xfit=fit_start:0.1:fit_end;
    
        [xxo, yyo, xo, no] = step(temp_SE(inrange(temp_SE,[fit_start fit_end])),[0:0.2:60]);
        s = fitoptions('Method','NonlinearLeastSquares','Robust','On','Algorithm','Trust-Region','Lower',[-5,0,0,0],'Upper',[5,Inf,Inf,Inf],...
        'Startpoint',[1,max(no),10,5]);
        skewGauss1 = fittype('a1/c1*exp(-((x-b1)/c1)^2)*(1+erf(a*(x-b1)/c1))','options',s);

    [f_both_xy,g_both_xy] = fit(xo',no',skewGauss1);
    delta_both_xy = f_both_xy.a/sqrt(1+f_both_xy.a^2);
    SE_mu_both_xy = f_both_xy.b1 + f_both_xy.c1/sqrt(2)*delta_both_xy*sqrt(2/pi);
    SE_sig_both_xy = sqrt((f_both_xy.c1/sqrt(2))^2 * (1 - 2 * (delta_both_xy)^2/pi));
    skew_both_xy = (4-pi)/2 * (delta_both_xy * sqrt(2/pi))^3 / (1 - 2 * delta_both_xy^2/pi)^(3/2);
    kurt_both_xy = 2*(pi-3) * (delta_both_xy * sqrt(2/pi))^3 / (1 - 2 * delta_both_xy^2/pi)^2;
    ci_both_xy = confint(f_both_xy,0.683);
    y_fit_both_xy = f_both_xy(xfit);   
      
    mean_SE_both_xy(y_count,x_count)=SE_mu_both_xy;  
    skew_SE_both_xy(y_count,x_count)=skew_both_xy; 
    sigma_SE_both_xy(y_count,x_count)=SE_sig_both_xy;
    kurt_SE_both_xy(y_count,x_count)=kurt_both_xy;
    mean_err_SE_both_xy(y_count,x_count)=SE_sig_mu_both_xy;
    sigma_err_SE_both_xy(y_count,x_count)=SE_sig_sig_both_xy;    
    

      else
          
        mean_SE_both_xy(y_count,x_count)=0;  
        skew_SE_both_xy(y_count,x_count)=0; 
        sigma_SE_both_xy(y_count,x_count)=0;
        kurt_SE_both_xy(y_count,x_count)=0;
        mean_err_SE_both_xy(y_count,x_count)=0;
        sigma_err_SE_both_xy(y_count,x_count)=0;    
    
      end           
    end    
  end

%%Plot the mean SE_XY both %%%%%%%%%%%%%%%%
color_range_max=max(max(mean_SE_both_xy));
color_range_min=min(min(mean_SE_both_xy(mean_SE_both_xy>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
 
SE_xy_both_fig = figure;
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


%%Plot the 1 sigma of the mean

color_range_max=max(max(mean_err_SE_both_xy));
color_range_min=min(min(mean_err_SE_both_xy(mean_err_SE_both_xy>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
 
sigma_SE_xy_both_fig = figure;
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

%Calculate corrections matrix and 1 sigma corrections matrix for SE_both
% if inpaint_on==1;  mean_SE_both_xy(mean_SE_both_xy==0)=nan; mean_SE_both_xy=inpaint_nans(mean_SE_both_xy,3); end;
SE_both_mu_center=interp2(se_xbins,se_ybins,mean_SE_both_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
norm_SE_both_mu=SE_both_mu_center./mean_SE_both_xy; 
norm_SE_both_mu(isinf(norm_SE_both_mu))=1;%no infinity! Do not correct events outside r=25cm
norm_SE_both_mu(isnan(norm_SE_both_mu))=1;

%Calculate error on SE center
% if inpaint_on==1;  mean_err_SE_both_xy(mean_err_SE_both_xy==0)=nan; mean_err_SE_both_xy=inpaint_nans(mean_err_SE_both_xy,3); end;
error_SE_both_mu_center=interp2(se_xbins,se_ybins,mean_err_SE_both_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
error_SE_both_mu_norm=sqrt((error_SE_both_mu_center./mean_SE_both_xy).^2+(mean_err_SE_both_xy.*SE_both_mu_center./mean_SE_both_xy.^2).^2);
error_SE_both_mu_norm(isinf(error_SE_both_mu_norm))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
error_SE_both_mu_norm(isnan(error_SE_both_mu_norm))=1;%no nan. Set Sigma=1, which is 100% uncertainty

%Calculate kurt,sigma, and skew at center
% if inpaint_on==1;  kurt_SE_both_xy(kurt_SE_both_xy==0)=nan; kurt_SE_both_xy=inpaint_nans(kurt_SE_both_xy,3); end;
% if inpaint_on==1;  skew_SE_both_xy(skew_SE_both_xy==0)=nan; skew_SE_both_xy=inpaint_nans(skew_SE_both_xy,3); end;
% if inpaint_on==1;  sigma_SE_both_xy(sigma_SE_both_xy==0)=nan; sigma_SE_both_xy=inpaint_nans(sigma_SE_both_xy,3); end;
SE_both_kurt_center=interp2(se_xbins,se_ybins,kurt_SE_both_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
SE_both_skew_center=interp2(se_xbins,se_ybins,skew_SE_both_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
SE_both_sigma_center=interp2(se_xbins,se_ybins,sigma_SE_both_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)


%% Assuming working corrections, find the optimal g2 with g1=0.105 from Kr2.18

s2_c=d.xyz_corrected_pulse_area_all_phe(s2_single_cut);
s1_c=d.xyz_corrected_pulse_area_all_phe(s1_single_cut);
% s2_c=d.z_corrected_pulse_area_all_phe(s2_single_cut);
% s1_c=d.pulse_area_phe(s1_single_cut).*d.correction_s1_z_dependence(s1_single_cut).*d.correction_s1_xy_dependence(s1_single_cut);

% cut=inrange(log10(s2_c),[4.1,5.1]) & inrange(log10(s1_c),[2 2.6]);
    
% for g2=15:0.01:22;
    