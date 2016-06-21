% Load in the S1a/S1b -> ER Band map

load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p22\KrypCal2p22_QualityChecks\ERBand_Prediction\S1aS1b_To_ERBand_RZVoxels.mat')


%% Load in Sep2015 ER S1a/S1b

load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p22\KrypCal2p22_QualityChecks\ERBand_Prediction\PredictVoxelizedBands\Feb2016_S1aS1b_Vox.mat')


SE_size_Sep2015=25.64;
SE_size_Sep2014=27.21;
SE_size_Nov2014=26.13;
SE_size_Feb2015=26.23;
SE_size_Feb2016=25.05;


voxel_counter=1;
for i=1:size(s1ab_z_means,2)
    MC_events=[];
    MC_events_x=[];

            if ~isnan(s1ab_z_means(i))
                [powerlaw_a powerlaw_a_err]=polyval(S1aS1b_to_mean_a_fit,s1ab_z_means(i),S1aS1b_to_mean_a_fit_err);
                [powerlaw_b powerlaw_b_err]=polyval(S1aS1b_to_mean_b_fit,s1ab_z_means(i),S1aS1b_to_mean_b_fit_err);
                [powerlaw_a_upper powerlaw_a_upper_err]=polyval(S1aS1b_to_upper_a_fit,s1ab_z_means(i),S1aS1b_to_upper_a_fit_err);
                [powerlaw_b_upper powerlaw_b_upper_err]=polyval(S1aS1b_to_upper_b_fit,s1ab_z_means(i),S1aS1b_to_upper_b_fit_err);
                [powerlaw_a_lower powerlaw_a_lower_err]=polyval(S1aS1b_to_lower_a_fit,s1ab_z_means(i),S1aS1b_to_lower_a_fit_err);
                [powerlaw_b_lower powerlaw_b_lower_err]=polyval(S1aS1b_to_lower_b_fit,s1ab_z_means(i),S1aS1b_to_lower_b_fit_err);
            
        %Make MC of predicted band for the voxel
                for counter=1:10000;
                     xx=50*rand(1); %pick random S1 value
                     yy= log10((10.^(powerlaw_a.*xx.^powerlaw_b)).*(SE_size_Feb2016/SE_size_Sep2015)); %Get associate log10(S1/S2) mean from power law fit
                     yy_sigma=abs(yy-log10((10.^(powerlaw_a_upper.*xx.^powerlaw_b_upper)).*(SE_size_Feb2016/SE_size_Sep2015)) )/1.28; %Get associate log10(S1/S2) gaussian width from power law fit
                     MC_events=[MC_events; normrnd(yy,yy_sigma,1)];
                     MC_events_x=[MC_events_x; xx];
               end
            end
            
        %Fit to the Zslice MC Band
        s1_min=0;
        s1_max=50;
        band_sigma=1.28; %90% confidence bounds
            
           s2_s1_bin= 1:.05:5;
           ER_bin_size=2; %3 works
           index_i=1;
           ER_bins=ER_bin_size/2+s1_min:ER_bin_size:s1_max;
               lower_bound=ones(size(ER_bins));
               mean_ER_band=ones(size(ER_bins));
               upper_bound=ones(size(ER_bins));
               MC_counts_ER_hist=ones(size(ER_bins));
               ER_sigma=ones(size(ER_bins));
               ER_sigma_sigma=ones(size(ER_bins));
        
%             temporary, for calculating correct E-lifetime   
        for ER_bin = (ER_bin_size/2)+s1_min:ER_bin_size:(s1_max); %HIST uses bin centers
   
        
            ER_fit_cut= inrange(MC_events_x,[ER_bin-ER_bin_size/2,ER_bin+ER_bin_size/2]);
      
            ER_band_hist=hist(MC_events(ER_fit_cut),s2_s1_bin);
            MC_counts_ER_hist(index_i)=sum(ER_band_hist);
           
            Fit_ER=fit(s2_s1_bin',ER_band_hist','gauss1');
             
            rms_fit=fit_ML_normal(MC_events(ER_fit_cut),s2_s1_bin);
            
            %band=confint(Fit_ER,.1);
            %ER_sigma(index_i)= Fit_ER.c1/sqrt(2);
            ER_sigma(index_i)= sqrt(rms_fit.sig2);
            ER_sigma_sigma(index_i)=ER_sigma(index_i)/sqrt(2*MC_counts_ER_hist(index_i));
            
            lower_bound(index_i)=Fit_ER.b1-ER_sigma(index_i)*band_sigma;
            mean_ER_band(index_i)=Fit_ER.b1;
            upper_bound(index_i)=Fit_ER.b1+ER_sigma(index_i)*band_sigma;
            
            %mean_ER_band_sim(index_i)=Fit_ER_sim.b1;
            
            index_i=index_i+1;
         end
        
         
         % Fit power law to the mean
         g=fittype('a*x^b');
         MC_ER_mean_power_fit=fit(ER_bins(ER_bins>2)',mean_ER_band(ER_bins>2)',g,'startpoint',[ 2.5 -.1]);
         MC_ER_mean_fit=MC_ER_mean_power_fit.a.*ER_bins.^(MC_ER_mean_power_fit.b);
            b=confint(MC_ER_mean_power_fit,.68);
            MC_sigma_a_ER_mean=(b(2,1)-b(1,1))/2;
            MC_sigma_b_ER_mean=(b(2,2)-b(1,2))/2;
            
         MC_ER_lower_power_fit=fit(ER_bins(ER_bins>2)',lower_bound(ER_bins>2)',g,'startpoint',[ 2.5 -.1]);
         MC_ER_lower_fit=MC_ER_lower_power_fit.a.*ER_bins.^(MC_ER_lower_power_fit.b);
            b=confint(MC_ER_lower_power_fit,.68);
            MC_sigma_a_ER_lower=(b(2,1)-b(1,1))/2;
            MC_sigma_b_ER_lower=(b(2,2)-b(1,2))/2;
            
         MC_ER_upper_power_fit=fit(ER_bins(ER_bins>2)',upper_bound(ER_bins>2)',g,'startpoint',[ 2.5 -.1]);
         MC_ER_upper_fit=MC_ER_upper_power_fit.a.*ER_bins.^(MC_ER_upper_power_fit.b);
            b=confint(MC_ER_upper_power_fit,.68);
            MC_sigma_a_ER_upper=(b(2,1)-b(1,1))/2;
            MC_sigma_b_ER_upper=(b(2,2)-b(1,2))/2;
            
            
     MC_ER_lower_power_fit_vox{voxel_counter}=MC_ER_lower_power_fit;
     MC_ER_upper_power_fit_vox{voxel_counter}=MC_ER_upper_power_fit;
     MC_ER_mean_power_fit_vox{voxel_counter}=MC_ER_mean_power_fit;
     MC_ydata{voxel_counter}=MC_events;     
     MC_xdata{voxel_counter}=MC_events_x;
     
     voxel_counter=voxel_counter+1;
end

%% Plot Each Z Slice band and prediction

load('Feb2016_ER_Vox.mat');

x=[1:0.1:50];
for i =1:length(MC_ER_upper_power_fit_vox)
    
    ER_mean_power_fit=ER_mean_power_fit_vox{i};
    ER_upper_power_fit=ER_upper_power_fit_vox{i};
    ER_lower_power_fit=ER_lower_power_fit_vox{i};
    
    MC_ER_mean_power_fit=MC_ER_mean_power_fit_vox{i};
    MC_ER_upper_power_fit=MC_ER_upper_power_fit_vox{i};
    MC_ER_lower_power_fit=MC_ER_lower_power_fit_vox{i};
    MC_y=MC_ydata{i};     
    MC_x=MC_xdata{i};
    
%Real data band
band_mean_data=ER_mean_power_fit.a*x.^ER_mean_power_fit.b;
band_upper_data=ER_upper_power_fit.a.*x.^ER_upper_power_fit.b;
band_lower_data=ER_lower_power_fit.a.*x.^ER_lower_power_fit.b;
band_width_data=(band_mean_data-band_lower_data)./(1.28);



%KrypCal Band - uses slightly different ER band than the one in the figures
band_mean_mc=MC_ER_mean_power_fit.a.*x.^MC_ER_mean_power_fit.b;
band_upper_mc=MC_ER_upper_power_fit.a*x.^MC_ER_upper_power_fit.b;
band_lower_mc=MC_ER_lower_power_fit.a*x.^MC_ER_lower_power_fit.b;
band_width_mc=(band_mean_mc-band_lower_mc)./(1.28);


%Plot the bands
figure
hold on;
plot(MC_x,MC_y,'.k');
plot(x,band_mean_mc,'-b','LineWidth',2);
plot(x,band_mean_data,'-r','LineWidth',2);
legend('MC Data','MC Fit','Data Fit');
xlabel('S1 (phe)'); ylabel('log10(S2/S1)'); myfigview(16);
plot(x,band_upper_mc,'--b','LineWidth',2);
plot(x,band_lower_mc,'--b','LineWidth',2);
plot(x,band_upper_data,'--r','LineWidth',2);
plot(x,band_lower_data,'--r','LineWidth',2);
ylim([1.5 3.5])


%Improvement in bands: Negative is worse resolution
width_difference=mean((band_width_data-band_width_mc)./(band_width_data))

%Shift in bands: Positive is upward shift
mean_shift=-mean((band_mean_data-band_mean_mc)./(band_mean_data))

%Chi2 for band mean
temp=confint(ER_mean_power_fit);
ER_mean_power_fit_a_err=abs(ER_mean_power_fit.a-temp(1,1));
ER_mean_power_fit_b_err=abs(ER_mean_power_fit.b-temp(1,2));
band_mean_data_err=sqrt( ((x.^ER_mean_power_fit.b).*(ER_mean_power_fit_a_err)).^2 + (ER_mean_power_fit_b_err.*ER_mean_power_fit.a.*(x.^ER_mean_power_fit.b).*log(x)).^2);
mean_chi2=sum( ((band_mean_mc-band_mean_data)./band_mean_data_err).^2)
mean_p_value=1-chi2cdf(mean_chi2,length(band_mean_data)-1)

%Chi2 for upper 90%
temp=confint(ER_upper_power_fit);
ER_upper_power_fit_a_err=abs(ER_upper_power_fit.a-temp(1,1));
ER_upper_power_fit_b_err=abs(ER_upper_power_fit.b-temp(1,2));
band_upper_data_err=sqrt( ((x.^ER_upper_power_fit.b).*(ER_upper_power_fit_a_err)).^2 + (ER_upper_power_fit_b_err.*ER_upper_power_fit.a.*(x.^ER_upper_power_fit.b).*log(x)).^2);
upper_chi2=sum( ((band_upper_mc-band_upper_data)./band_upper_data_err).^2)
upper_p_value=1-chi2cdf(upper_chi2,length(band_mean_data)-1)

%Chi2 for lower 90%
temp=confint(ER_lower_power_fit);
ER_lower_power_fit_a_err=abs(ER_lower_power_fit.a-temp(1,1));
ER_lower_power_fit_b_err=abs(ER_lower_power_fit.b-temp(1,2));
band_lower_data_err=sqrt( ((x.^ER_lower_power_fit.b).*(ER_lower_power_fit_a_err)).^2 + (ER_lower_power_fit_b_err.*ER_lower_power_fit.a.*(x.^ER_lower_power_fit.b).*log(x)).^2);
lower_chi2=sum( ((band_lower_mc-band_lower_data)./band_lower_data_err).^2)
lower_p_value=1-chi2cdf(lower_chi2,length(band_mean_data)-1)


%Chi2 for total width
expected_width=(band_upper_data-band_lower_data);
expected_width_err=sqrt(band_upper_data_err.^2 + band_lower_data_err.^2);
observed_width=(band_upper_mc-band_lower_mc);
width_chi2=sum( ((observed_width-expected_width)./(expected_width_err)).^2)
width_p_value=1-chi2cdf(width_chi2,length(band_mean_data)-1)



end

