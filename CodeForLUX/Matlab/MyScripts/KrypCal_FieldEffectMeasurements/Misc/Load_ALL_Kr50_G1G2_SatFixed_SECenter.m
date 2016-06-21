dir_path='/scratch4/LUXData/G1G2_Data/Kr_50Vcm/';

potential_folders_to_load={' '};
folders_to_load={' '};

% Check that we are only loading folders with lux10_
% and that they contain at least 10 MB of data
list=dir(dir_path);

    k=1;
for ii=1:length(list)
    
    if length(list(ii).name) > 6
        
        if  strcmp(list(ii).name(1:6),'lux10_') 
            
            folder_contents=dir(strcat(dir_path,'/',list(ii).name));
            dir_size=sum([folder_contents.bytes]);
            
            if dir_size > 10^7 % 10 MB
             
                potential_folders_to_load{k}=list(ii).name;
                k=k+1;
            
            end
                
                
        end
    end
    
end


% Once we determined good folders to load, check their multiplicity
% Keep higher CP unless duplicated folder with a higher CP number contains no data
folder_names_char=char(potential_folders_to_load);
folder_name_no_cp=cellstr([folder_names_char(:,1:19)]); %convert the char back to a cell array of chars

k=1;
good_index=1;
while good_index<=length(potential_folders_to_load)
        
        
       folder_multiplicity=sum(strcmp(potential_folders_to_load{good_index}(1:19),folder_name_no_cp ));
       good_index=good_index+folder_multiplicity-1;
            
       folders_to_load{k}=potential_folders_to_load{good_index};
       
       k=k+1;     
       good_index=good_index+1;
    
end

%%

rqs_to_load = {'pulse_area_phe','event_timestamp_samples'...
   ,'pulse_classification' ...
   ,'z_drift_samples' , 's1s2_pairing','golden'...
   ,'top_bottom_ratio','z_corrected_pulse_area_all_phe','event_timestamp_samples'...
   ,'x_corrected','y_corrected','selected_s1_s2','top_bottom_ratio',...
   'correction_electron_lifetime','correction_s2_xy_dependence','correction_s1_xyz_dependence',...
   'correction_s2_xy_dependence_bot','correction_s1_xyz_dependence_bot','aft_t0_samples',...
   'correction_s2_xy_dependence','spike_count','xyz_corrected_pulse_area_all_phe','xyz_corrected_pulse_area_bot_phe','s2_rec','x_cm','y_cm'};
    
ii=1;

    path=strcat(dir_path,'/',folders_to_load{ii});
    
    d = LUXLoadMultipleRQMs_framework(path,rqs_to_load);
    d = se_wrt_first_s1_func(d);

%% VUV and Spike Corrections
    % These correction factors are put in here by hand; they are the ratio of
% the two gain versions, per channel.
gainCorr = [ 1.07852885,  1.02147302,  1.14704873,  1.06599796,  1.,1.06565116,  1.17260696,  1.23915269,  1.05237586,  1.04581427,1.19078187,  1.01920447,  1.10806574,  1.04368312,  1.09359576,        1.20713643,  1.01269677,  0.97496447,  1.15328135,  1.04632389,1.09657976,  1.0540953 ,  1.06966811,  1.18564767,  1.03347816,1.12664979,  1.21384238,  1.04435383,  1.096459  ,  1.07002231,1.13197335,  1.        ,  0.99399306,  1.01098228,  1.2029547 ,1.07295668,  1.0696273 ,  1.00686858,  1.08148498,  1.16649075,1.04850371,  1.05249835,  1.14049429,  1.00785073,  1.0943977 ,1.02962349,  1.08931321,  1.18235792,  1.03934669,  1.00709524,1.14178315,  1.04650374,  1.07221255,  1.03081791,  1.1338022 ,1.23879521,  0.9849912 ,  1.02989339,  1.20126104,  1.0419401 ,1.08283583,  1.05111627,  1.16620644,  1.19375442,  1.06260808,1.01510865,  1.13129887,  1.09737434,  1.05748364,  1.08007699,1.07439593,  1.2135306 ,  1.05225111,  1.06547247,  1.13095721,1.07715678,  1.0901655 ,  1.06751444,  1.12705449,  1.17057193,1.09581789,  1.09417162,  1.23750879,  1.0706705 ,  1.18163888,1.10125419,  1.14468798,  1.19880273,  1.12981092,  1.06592446,1.19612131,  1.14547149,  1.        ,  1.06916959,  1.1169941 ,1.15057597,  1.08348327,  1.08020233,  1.13099728,  1.08311348,1.11291194,  1.04754566,  1.10311997,  1.18206575,  1.0658937 ,1.03623549,  1.17802513,  1.07419023,  1.04725227,  1.08861405,1.10805692,  1.19684977,  1.07895852,  1.05339989,  1.15828588,1.08408482,  1.02837954,  1.05048981,  1.09229397,  1.14161773,1.01265111,  1.06664276];

    
n_sig=2;

rcut_min = 0;
rcut_max = 25;%cm
        s1area_bound_min = 0;
        s1area_bound_max = 10^5;

        s2area_bound_min = 50;
        s2area_bound_max = 2*10^6;
   
    time_wait=d.event_timestamp_samples/1e8/60 > 60; %cut out first hour of high field data
    time_wait_cut=repmat(time_wait,10,1);
   
    d.z_drift_samples(isnan(d.z_drift_samples)) = 0.0; % get rid of NaN
    d.phe_bottom=d.pulse_area_phe./(1+d.top_bottom_ratio); %bottom PMT pulse area

    events=sort(size(d.pulse_area_phe)); %The first element is the number of pulses. sometimes it's 5, sometimes it's 10

    s1_area_cut= inrange(d.pulse_area_phe,[s1area_bound_min,s1area_bound_max]);
    s2_area_cut= inrange(d.pulse_area_phe,[s2area_bound_min,s2area_bound_max]);

    s1_class=time_wait_cut & (d.pulse_classification==1 )& s1_area_cut ; %area cut for Kr events
    s2_class=time_wait_cut & (d.pulse_classification==2) & s2_area_cut ;

    s1_single_cut = logical( (s1_class & d.golden & d.selected_s1_s2).*repmat(sum(s2_class & d.golden & d.selected_s1_s2)==1,events(1),1) ); % was s1_before_s2_cut before using golden
    s2_single_cut = logical( (s2_class & d.golden & d.selected_s1_s2).*repmat(sum(s1_class & d.golden & d.selected_s1_s2)==1,events(1),1) );
   
    s2_rec_c=d.s2_rec;

    s1=d.pulse_area_phe(s1_single_cut);
    s2=d.pulse_area_phe(s2_single_cut);
    s2_bot=d.phe_bottom(s2_single_cut);
    s2_c=s2_rec_c(s2_single_cut);
                  
    s1_phe_both_xyz_VUV=mean(gainCorr).*s1;
    s2_phe_both_xyz_VUV=mean(gainCorr).*s2;
    s2_phe_bot_xyz_VUV=mean(gainCorr).*s2_bot;
    s2_phe_both_xyz_VUV_SatFixed=mean(gainCorr).*s2_c;

    s2x_c=d.x_corrected(s2_single_cut);
    s2y_c=d.y_corrected(s2_single_cut);
    drift_time = d.z_drift_samples(s2_single_cut)/100; %units of us
             
    %% No tracking of XeAct means, sigma, errors, and number of events since we don't do the cuts in this code

%%

 %   Calculate Single Electron Size 
    s4_class=time_wait_cut & (d.pulse_classification==4) ;
    cut_se_s1_z = d.t_btw_se_first_s1 > 100*10;

    SE_good_cut =logical( cut_se_s1_z & s4_class & ( cumsum(s1_single_cut)-cumsum(s2_single_cut) ) ); % only events between good S1 and S2s    
    SE_phe_both_VUV=mean(gainCorr).*d.xyz_corrected_pulse_area_all_phe(SE_good_cut);

    SE_x=d.x_cm(SE_good_cut);
    SE_y=d.y_cm(SE_good_cut);
    SE_radius=(SE_x.^2+SE_y.^2).^(1/2);
    
    num_SE(ii)=length(SE_phe_both_VUV(SE_radius<12));

    SE_fit_VUV=fit([0:0.5:45]',hist(SE_phe_both_VUV(SE_radius<12),[0:0.5:45])','gauss1');
    SE_Size_VUV(ii)=SE_fit_VUV.b1;
    SE_sig_VUV(ii)=SE_fit_VUV.c1/sqrt(2);
    SE_Size_err_VUV(ii)=SE_sig_VUV(ii)/sqrt(length(SE_phe_both_VUV(SE_radius<12)));
    
    %Both PMT Second fit
    n_sig_fit=2;
    
    fit_start=SE_Size_VUV(ii)-n_sig_fit*SE_sig_VUV(ii);
    fit_end=SE_Size_VUV(ii)+n_sig_fit*SE_sig_VUV(ii);
    xfit=fit_start:0.1:fit_end;
    
    SE_cut=SE_phe_both_VUV(SE_radius<12);
    [xxo, yyo, xo, no] = step(SE_cut(inrange(SE_cut,[fit_start fit_end])));
    s = fitoptions('Method','NonlinearLeastSquares','Robust','On','Algorithm','Trust-Region','Lower',[-5,0,0,0],'Upper',[5,Inf,Inf,Inf],...
        'Startpoint',[1,max(no),10,5]);
    skewGauss1 = fittype('a1/c1*exp(-((x-b1)/c1)^2)*(1+erf(a*(x-b1)/c1))','options',s);

    [f_bot,g_bot] = fit(xo',no',skewGauss1);
    delta_bot = f_bot.a/sqrt(1+f_bot.a^2);
    SE_mu_both = f_bot.b1 + f_bot.c1/sqrt(2)*delta_bot*sqrt(2/pi);
    SE_sig_both = sqrt((f_bot.c1/sqrt(2))^2 * (1 - 2 * (delta_bot)^2/pi));
    skew_both = (4-pi)/2 * (delta_bot * sqrt(2/pi))^3 / (1 - 2 * delta_bot^2/pi)^(3/2);
    kurt_both = 2*(pi-3) * (delta_bot * sqrt(2/pi))^3 / (1 - 2 * delta_bot^2/pi)^2;
    ci_both = confint(f_bot);
    y_fit_both = f_bot(xfit);   
    
    SE_Size_VUV(ii)=SE_mu_both;
    SE_sig_VUV(ii)=SE_sig_both;
    
    S2_numelectrons_VUV=s2_phe_both_xyz_VUV./SE_Size_VUV(ii);
    S2_numelectrons_VUV_SatFixed=s2_phe_both_xyz_VUV_SatFixed./SE_Size_VUV(ii);
    
    %%
    
    %Bottom SE Size
     SE_phe_bot_VUV=mean(gainCorr).*d.xyz_corrected_pulse_area_bot_phe(SE_good_cut);

    num_SE_bot(ii)=length(SE_phe_bot_VUV(SE_radius<12));

    SE_fit_VUV=fit([0:0.2:30]',hist(SE_phe_bot_VUV(SE_radius<12),[0:0.2:30])','gauss1');
    SE_Size_bot_VUV(ii)=SE_fit_VUV.b1;
    SE_sig_bot_VUV(ii)=SE_fit_VUV.c1/sqrt(2);
    SE_Size_bot_err_VUV(ii)=SE_sig_VUV(ii)/sqrt(length(SE_phe_both_VUV(SE_radius<12)));
    
        %Bottom PMT Second fit
    n_sig_fit=2;
    
    fit_start=SE_Size_bot_VUV(ii)-n_sig_fit*SE_sig_VUV(ii);
    fit_end=SE_Size_bot_VUV(ii)+n_sig_fit*SE_sig_VUV(ii);
    xfit=fit_start:0.1:fit_end;
        
    SE_cut_bot=SE_phe_bot_VUV(SE_radius<12);

    [xxo, yyo, xo, no] = step(SE_cut_bot(inrange(SE_cut_bot,[fit_start fit_end])));
    s = fitoptions('Method','NonlinearLeastSquares','Robust','On','Algorithm','Trust-Region','Lower',[-5,0,0,0],'Upper',[5,Inf,Inf,Inf],...
        'Startpoint',[1,max(no),10,5]);
    skewGauss1 = fittype('a1/c1*exp(-((x-b1)/c1)^2)*(1+erf(a*(x-b1)/c1))','options',s);

    [f_bot,g_bot] = fit(xo',no',skewGauss1);
    delta_bot = f_bot.a/sqrt(1+f_bot.a^2);
    SE_mu_bot = f_bot.b1 + f_bot.c1/sqrt(2)*delta_bot*sqrt(2/pi);
    SE_sig_bot = sqrt((f_bot.c1/sqrt(2))^2 * (1 - 2 * (delta_bot)^2/pi));
    skew_bot = (4-pi)/2 * (delta_bot * sqrt(2/pi))^3 / (1 - 2 * delta_bot^2/pi)^(3/2);
    kurt_bot = 2*(pi-3) * (delta_bot * sqrt(2/pi))^3 / (1 - 2 * delta_bot^2/pi)^2;
    ci_bot = confint(f_bot);
    y_fit_bot = f_bot(xfit);  
    
    SE_Size_bot_VUV(ii)=SE_mu_both;
    SE_sig_bot_VUV(ii)=SE_sig_both;
    
    S2_numelectrons_bot_VUV=s2_phe_bot_xyz_VUV./SE_Size_bot_VUV(ii);

      save('/scratch4/LUXData/G1G2_Data/Kr_50Vcm/Run03Kr_50Vcm_SatCorr_SECenter_TimeCut','num_SE','SE_Size_VUV','SE_Size_err_VUV','SE_sig_VUV',...
      's1_phe_both_xyz_VUV',...
      's2_phe_both_xyz_VUV','S2_numelectrons_bot_VUV','SE_Size_bot_VUV','SE_sig_bot_VUV','SE_Size_bot_err_VUV','S2_numelectrons_VUV_SatFixed','s2_phe_both_xyz_VUV_SatFixed','s2x_c','s2y_c','drift_time','S2_numelectrons_VUV');        
