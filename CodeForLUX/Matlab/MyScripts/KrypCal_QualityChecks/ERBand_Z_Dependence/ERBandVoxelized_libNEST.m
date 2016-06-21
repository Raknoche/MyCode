load('LibNEST_ERData.mat')

s2radius_del=s2x_del+s2y_del;
s1_min=0;
s1_max=50;

max_r=25;

%% Loop over the 3x3x3 sections of the detetctor

z_step=floor((300-40)/3);
x_step=floor(50/3);
y_step=floor(50/3);

voxel_counter=1;

  ER_lower_power_fit_vox={};
  ER_upper_power_fit_vox={};
  ER_mean_power_fit_vox={};

for x_max=-25+x_step:x_step:25;
    for y_max=-25+y_step:y_step:25;
        for z_max=40+z_step:z_step:300;
            %note: voxels index in z, then y, then x
            
fiducial_cut=inrange(s2x_del,[x_max-x_step x_max]) & inrange(s2y_del,[y_max-y_step y_max]) & inrange(drift_time,[z_max-z_step z_max]);
   
tritium_cut= fiducial_cut;




%% Make Band Plot
   band_sigma=1.28; %90% confidence bounds
    
             
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
   
        
            ER_fit_cut=tritium_cut & inrange(s1_phe_both_xyz,[ER_bin-ER_bin_size/2,ER_bin+ER_bin_size/2]);
             
            ER_band_hist=hist(log10(s2_phe_both_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            counts_ER_hist(index_i)=sum(ER_band_hist);
           
            Fit_ER=fit(s2_s1_bin',ER_band_hist','gauss1');
             
            rms_fit=fit_ML_normal(log10(s2_phe_both_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            
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

            
  ER_lower_power_fit_vox{voxel_counter}=ER_lower_power_fit;
  ER_upper_power_fit_vox{voxel_counter}=ER_upper_power_fit;
  ER_mean_power_fit_vox{voxel_counter}=ER_mean_power_fit;
  x_center_vox{voxel_counter}=x_max-x_step/2;
  y_center_vox{voxel_counter}=y_max-y_step/2;
  z_center_vox{voxel_counter}=z_max-z_step/2;
   
  voxel_counter=voxel_counter+1;
        end
    end
end

%% Plot the ER bands
save('Sep2015_2p22_LibNEST_ERFit_Voxels','ER_lower_power_fit_vox','ER_upper_power_fit_vox','ER_mean_power_fit_vox','x_center_vox','y_center_vox','z_center_vox');

% load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p22\KrypCal2p22_QualityChecks\ERBand_Prediction\Sep2015_ERFit.mat');

figure
hold on;
x=[0:0.1:50];
% y_total=ER_mean_power_fit.a.*x.^ER_mean_power_fit.b;
% y_total_upper=ER_upper_power_fit.a.*x.^ER_upper_power_fit.b;
% y_total_lower=ER_lower_power_fit.a.*x.^ER_lower_power_fit.b;
plot(x,y_total,'-k','LineWidth',2,'Color',[0 0 0])

for i=1:(length(ER_mean_power_fit_vox));
    clear y
    y=ER_mean_power_fit_vox{i}.a.*x.^ER_mean_power_fit_vox{i}.b;
    plot(x,y,'-','LineWidth',1,'Color',[0.8 0.8 0.8])
    
    y=ER_lower_power_fit_vox{i}.a.*x.^ER_lower_power_fit_vox{i}.b;
    plot(x,y,'--','LineWidth',1,'Color',[0.8 0.8 0.8])

    y=ER_upper_power_fit_vox{i}.a.*x.^ER_upper_power_fit_vox{i}.b;
    plot(x,y,'--','LineWidth',1,'Color',[0.8 0.8 0.8])    
end
xlabel('Corrected S1 (phe)'); ylabel('Corrected log10(S2/S1)'); myfigview(16);
ylim([1.8 3])    
plot(x,y_total,'-k','LineWidth',2,'Color',[0 0 0])
plot(x,y_total_upper,'--k','LineWidth',2,'Color',[0 0 0])
plot(x,y_total_lower,'--k','LineWidth',2,'Color',[0 0 0])
legend('Total Band','Voxelized Band Means');

 
