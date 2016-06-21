load('Run03_Dec_Golden_CH3T.mat')

radius=(x_corr.^2 + y_corr.^2).^(1/2);

s1_min=0;
s1_max=50;
    
%     liquid_cut= inrange(drift_time,[30 300]) & inrange(s2radius_c,[0,18]) & inrange(s1area,[s1_min s1_max]) ...
%     & (log10(s1area)+0.5*log10(s2area)<3.8) & s2_phe_both>100 ; %used to be drift time 100 to 250

   
    band_sigma=4.5; %4.5 sigma = 1 in 150,000 chance of being outside this
    
             
            s2_s1_bin= 1:.05:5;
           ER_bin_size=2; %3 works
           index_i=1;
           ER_bins=ER_bin_size/2+s1_min:ER_bin_size:s1_max;
               lower_bound=ones(size(ER_bins));
               mean_ER_band=ones(size(ER_bins));
               upper_bound=ones(size(ER_bins));
                counts_ER_hist=ones(size(ER_bins));
               ER_sigma=ones(size(ER_bins));
               ER_sigma_sigma=ones(size(ER_bins));
        
%             temporary, for calculating correct E-lifetime   
        for ER_bin = (ER_bin_size/2)+s1_min:ER_bin_size:(s1_max); %HIST uses bin centers
   
        
            ER_fit_cut=inrange(s1area,[ER_bin-ER_bin_size/2,ER_bin+ER_bin_size/2]);
             
            ER_band_hist=hist(log10(s2area(ER_fit_cut)./s1area(ER_fit_cut)),s2_s1_bin);
            counts_ER_hist(index_i)=sum(ER_band_hist);
           
            Fit_ER=fit(s2_s1_bin',ER_band_hist','gauss1');
             
            rms_fit=fit_ML_normal(log10(s2area(ER_fit_cut)./s1area(ER_fit_cut)),s2_s1_bin);
            
            %band=confint(Fit_ER,.1);
            %ER_sigma(index_i)= Fit_ER.c1/sqrt(2);
            ER_sigma(index_i)= sqrt(rms_fit.sig2);
            ER_sigma_sigma(index_i)=ER_sigma(index_i)/sqrt(2*counts_ER_hist(index_i));
            
            lower_bound(index_i)=Fit_ER.b1-ER_sigma(index_i)*band_sigma;
            mean_ER_band(index_i)=Fit_ER.b1;
            upper_bound(index_i)=Fit_ER.b1+ER_sigma(index_i)*band_sigma;
            
            %mean_ER_band_sim(index_i)=Fit_ER_sim.b1;
            
            index_i=index_i+1;
         end
        
         
         % Fit power law to the mean
         g=fittype('a*x^b');
         ER_mean_power_fit=fit(ER_bins(ER_bins>2)',mean_ER_band(ER_bins>2)',g,'startpoint',[ 2.5 -.1]);
         ER_mean_fit=ER_mean_power_fit.a.*ER_bins.^(ER_mean_power_fit.b);
            b=confint(ER_mean_power_fit,.68);
            sigma_a_ER_mean=(b(2,1)-b(1,1))/2;
            sigma_b_ER_mean=(b(2,2)-b(1,2))/2;
            
         ER_lower_power_fit=fit(ER_bins(ER_bins>2)',lower_bound(ER_bins>2)',g,'startpoint',[ 2.5 -.1]);
         ER_lower_fit=ER_lower_power_fit.a.*ER_bins.^(ER_lower_power_fit.b);
            b=confint(ER_lower_power_fit,.68);
            sigma_a_ER_lower=(b(2,1)-b(1,1))/2;
            sigma_b_ER_lower=(b(2,2)-b(1,2))/2;
            
         ER_upper_power_fit=fit(ER_bins(ER_bins>2)',upper_bound(ER_bins>2)',g,'startpoint',[ 2.5 -.1]);
         ER_upper_fit=ER_upper_power_fit.a.*ER_bins.^(ER_upper_power_fit.b);
            b=confint(ER_upper_power_fit,.68);
            sigma_a_ER_upper=(b(2,1)-b(1,1))/2;
            sigma_b_ER_upper=(b(2,2)-b(1,2))/2;

            
  %Make Plot    

 no_fid_new=figure;
 hold on

 scatter(s1area,log10(s2area./s1area),'.k');
 ylabel('Corrected log10(S2/S1)'); xlabel('S1 corrected (Pulse Area Phe)'); 
 title('Corrected CH3T','fontsize',16,'Interpreter','none');
 myfigview(16);
 ylim([0 3]);xlim([0 s1_max]);
%  plot(ER_bins,lower_bound,'ok','markersize',6,'linewidth',2)
 plot(ER_bins,ER_lower_fit,'--r','markersize',14,'linewidth',2)
%  plot(ER_bins,mean_ER_band,'ok','markersize',6,'linewidth',2)%plot mean of ER band
 plot(ER_bins,ER_mean_fit,'-r','markersize',20,'linewidth',2)
%  plot(ER_bins,upper_bound,'ok','markersize',6,'linewidth',2)
 plot(ER_bins,ER_upper_fit,'--r','markersize',14,'linewidth',2)
 bin_step=2;

            total_count=length(log10(s2area));
 
    
 box;
  text(10,2.8,strcat('Total Count= ',num2str(total_count)),'fontsize',16);
  %text(5,.7,strcat('NEST\_SE_{bot} Fit=',num2str(ER_mean_NEST_fit.SE_ER_fit*73/64*27/24.5,3),'\pm ', num2str(sigma_SE_ER_fit*73/64*27/24.5,2),'[Phe]'),'fontsize',16,'color','m','interpreter','tex');
  %text(5,.8,'NEST Platinum V4-b', 'fontsize',16,'color','m','interpreter','tex' );
  text(5,.55,strcat('+',num2str(band_sigma),'\sigma = (',num2str(ER_upper_power_fit.a,3),'\pm', num2str(sigma_a_ER_upper,1) ,')\times S1_c ^{', num2str(ER_upper_power_fit.b,3),'\pm',num2str(sigma_b_ER_upper,1),'}'),'fontsize',16,'interpreter','tex');
  text(5,.3,strcat('Mean=(',num2str(ER_mean_power_fit.a,3),'\pm', num2str(sigma_a_ER_mean,1),')\times S1_c ^{', num2str(ER_mean_power_fit.b,3),'\pm', num2str(sigma_b_ER_mean,1),'}'),'fontsize',16,'interpreter','tex');
  text(5,.05,strcat('-',num2str(band_sigma),'\sigma = (',num2str(ER_lower_power_fit.a,3),'\pm', num2str(sigma_a_ER_lower,1),')\times S1_c ^{', num2str(ER_lower_power_fit.b,3),'\pm', num2str(sigma_b_ER_lower,1),'}'),'fontsize',16,'interpreter','tex');
  

  
%% Cutting on populations above and below the band width

below_band_cut=log10(s2area./s1area)< (ER_lower_power_fit.a.*(s1area.^ER_lower_power_fit.b));
above_band_cut=log10(s2area./s1area)> (ER_upper_power_fit.a.*(s1area.^ER_upper_power_fit.b));
hold on;
scatter(s1area(below_band_cut),log10(s2area(below_band_cut)./s1area(below_band_cut)),'.m');
scatter(s1area(above_band_cut),log10(s2area(above_band_cut)./s1area(above_band_cut)),'.b');

%% Plotting XY dist of below band population

figure
scatter(x_corr(below_band_cut),y_corr(below_band_cut),'.k')
hold on;
circle(0,0,25);
xlabel('Corrected x (cm)');ylabel('Corrected y (cm)');
myfigview(16);

%% Checking corrections are okay by comparing raw and corrected pulse area
figure
scatter(s1area_raw,s1area,'.k');
hold on;
scatter(s1area_raw(below_band_cut),s1area(below_band_cut),'.m');
scatter(s1area_raw(above_band_cut),s1area(above_band_cut),'.b');
xlabel('Uncorrected S1');ylabel('Corrected S1'); myfigview(16);
legend('All Events','Below Band','Above Band');


figure
scatter(s2area_raw,s2area,'.k');
hold on;
scatter(s2area_raw(below_band_cut),s2area(below_band_cut),'.m');
scatter(s2area_raw(above_band_cut),s2area(above_band_cut),'.b');
xlabel('Uncorrected S2');ylabel('Corrected S2'); myfigview(16);
legend('All Events','Below Band','Above Band');

%% Checking that corrections are okay by making the raw ER band and looking at events

            s2_s1_bin= 1:.05:5;
           ER_bin_size=2; %3 works
           index_i=1;
           ER_bins=ER_bin_size/2+s1_min:ER_bin_size:s1_max;
               lower_bound=ones(size(ER_bins));
               mean_ER_band=ones(size(ER_bins));
               upper_bound=ones(size(ER_bins));
                counts_ER_hist=ones(size(ER_bins));
               ER_sigma=ones(size(ER_bins));
               ER_sigma_sigma=ones(size(ER_bins));
        
%             temporary, for calculating correct E-lifetime   
        for ER_bin = (ER_bin_size/2)+s1_min:ER_bin_size:(s1_max); %HIST uses bin centers
   
        
            ER_fit_cut=inrange(s1area_raw,[ER_bin-ER_bin_size/2,ER_bin+ER_bin_size/2]);
             
            ER_band_hist=hist(log10(s2area_raw(ER_fit_cut)./s1area_raw(ER_fit_cut)),s2_s1_bin);
            counts_ER_hist(index_i)=sum(ER_band_hist);
           
            Fit_ER=fit(s2_s1_bin',ER_band_hist','gauss1');
             
            rms_fit=fit_ML_normal(log10(s2area_raw(ER_fit_cut)./s1area_raw(ER_fit_cut)),s2_s1_bin);
            
            %band=confint(Fit_ER,.1);
            %ER_sigma(index_i)= Fit_ER.c1/sqrt(2);
            ER_sigma(index_i)= sqrt(rms_fit.sig2);
            ER_sigma_sigma(index_i)=ER_sigma(index_i)/sqrt(2*counts_ER_hist(index_i));
            
            lower_bound(index_i)=Fit_ER.b1-ER_sigma(index_i)*band_sigma;
            mean_ER_band(index_i)=Fit_ER.b1;
            upper_bound(index_i)=Fit_ER.b1+ER_sigma(index_i)*band_sigma;
            
            %mean_ER_band_sim(index_i)=Fit_ER_sim.b1;
            
            index_i=index_i+1;
         end
        
         
         % Fit power law to the mean
         g=fittype('a*x^b');
         ER_mean_power_fit=fit(ER_bins(ER_bins>2)',mean_ER_band(ER_bins>2)',g,'startpoint',[ 2.5 -.1]);
         ER_mean_fit=ER_mean_power_fit.a.*ER_bins.^(ER_mean_power_fit.b);
            b=confint(ER_mean_power_fit,.68);
            sigma_a_ER_mean=(b(2,1)-b(1,1))/2;
            sigma_b_ER_mean=(b(2,2)-b(1,2))/2;
            
         ER_lower_power_fit=fit(ER_bins(ER_bins>2)',lower_bound(ER_bins>2)',g,'startpoint',[ 2.5 -.1]);
         ER_lower_fit=ER_lower_power_fit.a.*ER_bins.^(ER_lower_power_fit.b);
            b=confint(ER_lower_power_fit,.68);
            sigma_a_ER_lower=(b(2,1)-b(1,1))/2;
            sigma_b_ER_lower=(b(2,2)-b(1,2))/2;
            
         ER_upper_power_fit=fit(ER_bins(ER_bins>2)',upper_bound(ER_bins>2)',g,'startpoint',[ 2.5 -.1]);
         ER_upper_fit=ER_upper_power_fit.a.*ER_bins.^(ER_upper_power_fit.b);
            b=confint(ER_upper_power_fit,.68);
            sigma_a_ER_upper=(b(2,1)-b(1,1))/2;
            sigma_b_ER_upper=(b(2,2)-b(1,2))/2;

            
  %Make Plot    

 no_fid_new=figure;
 hold on

 scatter(s1area_raw,log10(s2area_raw./s1area_raw),'.k');
 ylabel('Uncorrected log10(S2/S1)'); xlabel('S1 Uncorrected (Pulse Area Phe)'); 
 title('Uncorrected CH3T','fontsize',16,'Interpreter','none');
 myfigview(16);
 ylim([0 3]);xlim([0 s1_max]);
%  plot(ER_bins,lower_bound,'ok','markersize',6,'linewidth',2)
 plot(ER_bins,ER_lower_fit,'--r','markersize',14,'linewidth',2)
%  plot(ER_bins,mean_ER_band,'ok','markersize',6,'linewidth',2)%plot mean of ER band
 plot(ER_bins,ER_mean_fit,'-r','markersize',20,'linewidth',2)
%  plot(ER_bins,upper_bound,'ok','markersize',6,'linewidth',2)
 plot(ER_bins,ER_upper_fit,'--r','markersize',14,'linewidth',2)
 bin_step=2;
scatter(s1area_raw(below_band_cut),log10(s2area_raw(below_band_cut)./s1area_raw(below_band_cut)),'.m');
scatter(s1area_raw(above_band_cut),log10(s2area_raw(above_band_cut)./s1area_raw(above_band_cut)),'.b');

%% Printing Luxstamps
belowband_events=luxstamp(below_band_cut);
aboveband_events=luxstamp(above_band_cut);

