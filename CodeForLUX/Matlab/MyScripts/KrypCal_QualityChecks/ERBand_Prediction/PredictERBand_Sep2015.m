% Load in the S1a/S1b -> ER Band map

load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p22\KrypCal2p22_QualityChecks\ERBand_Prediction\S1aS1b_To_ERBand_RZVoxels.mat')


%% Load in Sep2015 ER S1a/S1b

load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p22\KrypCal2p22_QualityChecks\ERBand_Prediction\Sep2015_S1aS1b.mat')


SE_size_Sep2015=25.64;
SE_size_Sep2014=27.21;
SE_size_Nov2014=26.13;
SE_size_Feb2015=26.23;


MC_events=[];
MC_events_x=[];

% 
% for i=1:size(s1ab_xyz_mean,1)
%     for j=1:size(s1ab_xyz_mean,2)
%         for k=1:size(s1ab_xyz_mean,3)
%             if ~isnan(s1ab_xyz_mean(i,j,k))
%                 [powerlaw_a powerlaw_a_err]=polyval(S1aS1b_to_mean_a_fit,s1ab_xyz_mean(i,j,k),S1aS1b_to_mean_a_fit_err);
%                 [powerlaw_b powerlaw_b_err]=polyval(S1aS1b_to_mean_b_fit,s1ab_xyz_mean(i,j,k),S1aS1b_to_mean_b_fit_err);
%                 [powerlaw_a_upper powerlaw_a_upper_err]=polyval(S1aS1b_to_upper_a_fit,s1ab_xyz_mean(i,j,k),S1aS1b_to_upper_a_fit_err);
%                 [powerlaw_b_upper powerlaw_b_upper_err]=polyval(S1aS1b_to_upper_b_fit,s1ab_xyz_mean(i,j,k),S1aS1b_to_upper_b_fit_err);
%                 [powerlaw_a_lower powerlaw_a_lower_err]=polyval(S1aS1b_to_lower_a_fit,s1ab_xyz_mean(i,j,k),S1aS1b_to_lower_a_fit_err);
%                 [powerlaw_b_lower powerlaw_b_lower_err]=polyval(S1aS1b_to_lower_b_fit,s1ab_xyz_mean(i,j,k),S1aS1b_to_lower_b_fit_err);
%             
%         %Make MC of predicted band for the voxel
%                 for counter=1:Kr_count_3D(i,j,k)/10;
%                      xx=50*rand(1); %pick random S1 value
%                      yy= powerlaw_a.*xx.^powerlaw_b; %Get associate log10(S1/S2) mean from power law fit
%                      yy_sigma=abs(yy-(powerlaw_a_upper.*xx.^powerlaw_b_upper))/1.28; %Get associate log10(S1/S2) gaussian width from power law fit
%                      MC_events=[MC_events; normrnd(yy,yy_sigma,1)];
%                      MC_events_x=[MC_events_x; xx];
%                end
%             end
%         end
%     end
% end



% for i=1:size(s1ab_r2z_mean,1)
%     for j=1:size(s1ab_r2z_mean,2)
%             if ~isnan(s1ab_r2z_mean(i,j))
%                 [powerlaw_a powerlaw_a_err]=polyval(S1aS1b_to_mean_a_fit,s1ab_r2z_mean(i,j),S1aS1b_to_mean_a_fit_err);
%                 [powerlaw_b powerlaw_b_err]=polyval(S1aS1b_to_mean_b_fit,s1ab_r2z_mean(i,j),S1aS1b_to_mean_b_fit_err);
%                 [powerlaw_a_upper powerlaw_a_upper_err]=polyval(S1aS1b_to_upper_a_fit,s1ab_r2z_mean(i,j),S1aS1b_to_upper_a_fit_err);
%                 [powerlaw_b_upper powerlaw_b_upper_err]=polyval(S1aS1b_to_upper_b_fit,s1ab_r2z_mean(i,j),S1aS1b_to_upper_b_fit_err);
%                 [powerlaw_a_lower powerlaw_a_lower_err]=polyval(S1aS1b_to_lower_a_fit,s1ab_r2z_mean(i,j),S1aS1b_to_lower_a_fit_err);
%                 [powerlaw_b_lower powerlaw_b_lower_err]=polyval(S1aS1b_to_lower_b_fit,s1ab_r2z_mean(i,j),S1aS1b_to_lower_b_fit_err);
%             
%         %Make MC of predicted band for the voxel
%                 for counter=1:Kr_count_2D(i,j)/10;
%                      xx=50*rand(1); %pick random S1 value
%                      yy= powerlaw_a.*xx.^powerlaw_b; %Get associate log10(S1/S2) mean from power law fit
%                      yy_sigma=abs(yy-(powerlaw_a_upper.*xx.^powerlaw_b_upper))/1.28; %Get associate log10(S1/S2) gaussian width from power law fit
%                      MC_events=[MC_events; normrnd(yy,yy_sigma,1)];
%                      MC_events_x=[MC_events_x; xx];
%                end
%             end
%      end
% end



for i=1:size(s1ab_z_means,2)
            if ~isnan(s1ab_z_means(i))
                [powerlaw_a powerlaw_a_err]=polyval(S1aS1b_to_mean_a_fit,s1ab_z_means(i),S1aS1b_to_mean_a_fit_err);
                [powerlaw_b powerlaw_b_err]=polyval(S1aS1b_to_mean_b_fit,s1ab_z_means(i),S1aS1b_to_mean_b_fit_err);
                [powerlaw_a_upper powerlaw_a_upper_err]=polyval(S1aS1b_to_upper_a_fit,s1ab_z_means(i),S1aS1b_to_upper_a_fit_err);
                [powerlaw_b_upper powerlaw_b_upper_err]=polyval(S1aS1b_to_upper_b_fit,s1ab_z_means(i),S1aS1b_to_upper_b_fit_err);
                [powerlaw_a_lower powerlaw_a_lower_err]=polyval(S1aS1b_to_lower_a_fit,s1ab_z_means(i),S1aS1b_to_lower_a_fit_err);
                [powerlaw_b_lower powerlaw_b_lower_err]=polyval(S1aS1b_to_lower_b_fit,s1ab_z_means(i),S1aS1b_to_lower_b_fit_err);
            
        %Make MC of predicted band for the voxel
                for counter=1:Kr_count_1D(i);
                     xx=50*rand(1); %pick random S1 value
                     yy= powerlaw_a.*xx.^powerlaw_b; %Get associate log10(S1/S2) mean from power law fit
                     yy_sigma=abs(yy-(powerlaw_a_upper.*xx.^powerlaw_b_upper))/1.28; %Get associate log10(S1/S2) gaussian width from power law fit
                     MC_events=[MC_events; normrnd(yy,yy_sigma,1)];
                     MC_events_x=[MC_events_x; xx];
               end
            end
end

%% Fit total MC band
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

        MC_powerlaw_a=MC_ER_mean_power_fit.a;
        MC_powerlaw_b=MC_ER_mean_power_fit.b;
        MC_powerlaw_a_upper=MC_ER_upper_power_fit.a;
        MC_powerlaw_b_upper=MC_ER_upper_power_fit.b;
        MC_powerlaw_a_lower=MC_ER_lower_power_fit.a;
        MC_powerlaw_b_lower=MC_ER_lower_power_fit.b;
        MC_powerlaw_a_err=MC_sigma_a_ER_mean;
        MC_powerlaw_b_err=MC_sigma_b_ER_mean;
        MC_powerlaw_a_upper_err=MC_sigma_a_ER_upper;
        MC_powerlaw_b_upper_err=MC_sigma_b_ER_upper;
        MC_powerlaw_a_lower_err=MC_sigma_a_ER_lower;
        MC_powerlaw_b_lower_err=MC_sigma_b_ER_lower;
        
  %% Compare to actual band

load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p22\KrypCal2p22_QualityChecks\ERBand_Prediction\Sep2015_ERFit.mat')

  %Make Plot    

 no_fid_new=figure;
 hold on
 plot(MC_events_x,MC_events,'.k','MarkerSize',1)
%  cc=colorbar;
%     ylabel(cc,'Count','FontSize',16);
%  colormap;
 ylabel('Corrected log10(S2/S1)'); xlabel('S1 corrected (Pulse Area Phe)'); 
 title(sprintf('Run 04 CH3T, MonteCarlo Band'),'fontsize',16,'Interpreter','none');
 myfigview(16);
 ylim([0 3.5]);xlim([0 s1_max]);
%  plot(ER_bins,lower_bound,'ok','markersize',6,'linewidth',2)
plot(ER_bins,MC_ER_mean_fit,'-b','markersize',20,'linewidth',2)
plot(ER_bins,ER_mean_fit,'-r','markersize',14,'linewidth',2)
%  plot(ER_bins,mean_ER_band,'ok','markersize',6,'linewidth',2)%plot mean of ER band
 plot(ER_bins,ER_lower_fit,'--r','markersize',20,'linewidth',2)
%  plot(ER_bins,upper_bound,'ok','markersize',6,'linewidth',2)
 plot(ER_bins,ER_upper_fit,'--r','markersize',14,'linewidth',2)
plot(ER_bins,MC_ER_lower_fit,'--b','markersize',20,'linewidth',2)
plot(ER_bins,MC_ER_upper_fit,'--b','markersize',20,'linewidth',2)
  legend('Data','MC Band','Measured ER Band');

 
     bin_step=2;

            total_count=length(log10(MC_events_x));
 
%     
%  box;
%   text(10,2.8,strcat('Total Count= ',num2str(total_count)),'fontsize',16);
%   text(5,.55,strcat('+',num2str(band_sigma),'\sigma = (',num2str(ER_upper_power_fit.a,3),'\pm', num2str(sigma_a_ER_upper,1) ,')\times S1_c ^{', num2str(ER_upper_power_fit.b,3),'\pm',num2str(sigma_b_ER_upper,1),'}'),'fontsize',16,'interpreter','tex');
%   text(5,.3,strcat('Mean=(',num2str(ER_mean_power_fit.a,3),'\pm', num2str(sigma_a_ER_mean,1),')\times S1_c ^{', num2str(ER_mean_power_fit.b,3),'\pm', num2str(sigma_b_ER_mean,1),'}'),'fontsize',16,'interpreter','tex');
%   text(5,.05,strcat('-',num2str(band_sigma),'\sigma = (',num2str(ER_lower_power_fit.a,3),'\pm', num2str(sigma_a_ER_lower,1),')\times S1_c ^{', num2str(ER_lower_power_fit.b,3),'\pm', num2str(sigma_b_ER_lower,1),'}'),'fontsize',16,'interpreter','tex');
%   title('Sep2015 CH3T')
  
   box;

  text(5,.55,strcat('+',num2str(band_sigma),'\sigma = (',num2str(MC_ER_upper_power_fit.a,3),'\pm', num2str(MC_sigma_a_ER_upper,1) ,')\times S1_c ^{', num2str(MC_ER_upper_power_fit.b,3),'\pm',num2str(MC_sigma_b_ER_upper,1),'}'),'fontsize',16,'interpreter','tex');
  text(5,.3,strcat('Mean=(',num2str(MC_ER_mean_power_fit.a,3),'\pm', num2str(MC_sigma_a_ER_mean,1),')\times S1_c ^{', num2str(MC_ER_mean_power_fit.b,3),'\pm', num2str(MC_sigma_b_ER_mean,1),'}'),'fontsize',16,'interpreter','tex');
  text(5,.05,strcat('-',num2str(band_sigma),'\sigma = (',num2str(MC_ER_lower_power_fit.a,3),'\pm', num2str(MC_sigma_a_ER_lower,1),')\times S1_c ^{', num2str(MC_ER_lower_power_fit.b,3),'\pm', num2str(MC_sigma_b_ER_lower,1),'}'),'fontsize',16,'interpreter','tex');
  title('Sep2015 CH3T')
  

   
   
%% Stats about the difference


%Real Data Band 
x=[1:0.1:50];
band_mean_data=ER_mean_power_fit.a*x.^ER_mean_power_fit.b;
band_upper_data=ER_upper_power_fit.a.*x.^ER_upper_power_fit.b;
band_lower_data=ER_lower_power_fit.a.*x.^ER_lower_power_fit.b;
band_width_data=(band_mean_data-band_lower_data)./(1.28);

%KrypCal Band - uses slightly different ER band than the one in the figures
band_mean_mc=MC_ER_mean_power_fit.a.*x.^MC_ER_mean_power_fit.b;
band_upper_mc=MC_ER_upper_power_fit.a*x.^MC_ER_upper_power_fit.b;
band_lower_mc=MC_ER_lower_power_fit.a*x.^MC_ER_lower_power_fit.b;
band_width_mc=(band_mean_mc-band_lower_mc)./(1.28);

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
