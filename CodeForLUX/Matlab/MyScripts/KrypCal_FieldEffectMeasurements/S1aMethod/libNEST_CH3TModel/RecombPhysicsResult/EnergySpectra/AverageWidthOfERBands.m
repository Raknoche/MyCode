%% Sep 2014 Bands

%Field Neglected Band 
x=[1:0.1:50];
band_mean_fieldneg=3.17.*x.^-0.101;
band_upper_fieldneg=3.37.*x.^-0.103;
band_width_fieldneg=(band_upper_fieldneg-band_mean_fieldneg)./(1.28);
% average_band_width_fieldneg=mean(band_width_fieldneg);

%KrypCal Band - uses slightly different ER band than the one in the figures
band_mean_kr=2.99.*x.^-0.11;
band_upper_kr=3.19.*x.^-0.11;
% band_width_kr=(10.^band_upper_kr-10.^band_mean_kr)./(10.^band_mean_kr);
% band_width_kr=(10.^band_mean_kr-10.^band_lower_kr)./(10.^band_mean_kr);
% band_width_kr=(band_mean_kr-band_lower_kr)./(band_mean_kr);
band_width_kr=(band_upper_kr-band_mean_kr)./(1.28);

% average_band_width_kr=mean(band_width_kr);

%Improvement in bands: Negative is worse resolution
% kr_width_improvement=(average_band_width_fieldneg-average_band_width_kr)/(average_band_width_fieldneg)
kr_width_improvement=(band_width_fieldneg-band_width_kr)/(band_width_fieldneg);

mean(kr_width_improvement)

%Shift in bands: Negative is upward shift
kr_shift=(band_mean_fieldneg-band_mean_kr)/(band_mean_fieldneg)

%% Sep 2015 Bands
clear all
x=[1:0.1:50];

%Field Neglected Band - Mean: 3.23*s1^-0.101, +1.28sig: 3.43*s1^-0.101

band_mean_fieldneg=3.23.*x.^-0.101;
band_upper_fieldneg=3.43.*x.^-0.101;
band_width_fieldneg=(band_upper_fieldneg-band_mean_fieldneg)./(1.28);

%kr Band - Mean: 

band_mean_kr=2.99.*x.^-0.113;
band_upper_kr=3.16.*x.^-0.108;
% band_width_kr=(10.^band_upper_kr-10.^band_mean_kr)./(10.^band_mean_kr);
% band_width_kr=(10.^band_mean_kr-10.^band_lower_kr)./(10.^band_mean_kr);
% band_width_kr=(band_mean_kr-band_lower_kr)./(band_mean_kr);
band_width_kr=(band_upper_kr-band_mean_kr)./(1.28);

average_band_width_kr=mean(band_width_kr);

%Improvement in bands: Negative is worse resolution
kr_width_improvement=(band_width_fieldneg-band_width_kr)/(band_width_fieldneg);

mean(kr_width_improvement)

%Shift in bands: Negative is upward shift
kr_shift=(band_mean_fieldneg-band_mean_kr)/(band_mean_fieldneg)
