% SatFixed_AllSys_SEDQCenter_Cov_RPlot_Kr50TimeCut_NewCorr

%NOTE: WHEN PRESENTED TO AWG, THIS HAD A BUG IN THE 50 V/cm S2 XY
%CORRECTIONS -- IT WAS USING S1 XY CORRECTIONS INSTEAD.  THIS MADE G1 =
%0.12 and EE =0.432 instead of G1=0.117 and EE =0.445 (this is a 1 sigma
%diff EE)

%Fixes the s2 error = sigma/sqrt(num s1) bug
%Fixes the Energy to 203 instead of 208 bug 
%Fixes the hour of high field in 50 V/cm data bug
%Changes 50 V/cm det_middle from (463-4.23)/2= 229.4 to 184.5 based on counting events
%Changes 100 V/cm det_middle from (351-4.23)/2= 173.4 to 171.5 based on counting events
%Fixes 50 V/cm and 100 V/cm corrections to be based on uncorrected xy

%Does NOT fix the fact that r<18 cm cut is different for 50 and 100 V/cm data at different depths

%% Find g1 and g2 by fiting to multiple mono energetic peaks

     
     %% Kr83m at 42.6 keV and 170 V/cm

     clear radius s2x_c s2y_c drift_time s1_phe_both_xyz s2_phe_bottom_xyz S2_numelectrons_VUV
load('C:\Program Files\MATLAB\R2012a\bin\DokeWithDQSE\Run03Kr_170Vcm_DQ_SE.mat'); %Load Patched Kr Data

n_sig=3;
     radius=(s2x_c.^2+s2y_c.^2).^(1/2);
     
     
     %REDO ALL OF THIS WITH SAT FIXED AND S2 BOTH
     s1_phe_both_xyz_kr170=s1_phe_both_xyz_VUV(inrange(drift_time,[30 300]) & inrange(radius,[0 18]) & inrange(s2_phe_both_xyz_VUV_SatFixed,[1585,20000]) & inrange(s1_phe_both_xyz_VUV,[150,400]));
     s2_phe_bottom_xyz_kr170=S2_numelectrons_VUV_SatFixed(inrange(drift_time,[30 300]) & inrange(radius,[0 18]) & inrange(s2_phe_both_xyz_VUV_SatFixed,[1585,20000]) & inrange(s1_phe_both_xyz_VUV,[150,400]));
    s2_phe_bottom_xyz_kr170_phe=s2_phe_both_xyz_VUV_SatFixed(inrange(drift_time,[30 300]) & inrange(radius,[0 18]) & inrange(s2_phe_both_xyz_VUV_SatFixed,[1585,20000]) & inrange(s1_phe_both_xyz_VUV,[150,400]));
     s1_phe_both_xyz_kr170_phe=s1_phe_both_xyz_VUV(inrange(drift_time,[30 300]) & inrange(radius,[0 18]) & inrange(s2_phe_both_xyz_VUV_SatFixed,[1585,20000]) & inrange(s1_phe_both_xyz_VUV,[150,400]));

     
          clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2   
     Fit_kr_s1=fit([0:1:600]',hist(s1_phe_both_xyz_kr170,[0:1:600])','gauss1');
     Fit_kr_s2=fit([0:2:1000]',hist(s2_phe_bottom_xyz_kr170,[0:2:1000])','gauss1');    
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_xyz_kr170)); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_xyz_kr170)); % == sigma/sqrt(N); 
     %cut and fit around 3 sigma of the mean, after the first guess 
     cut_kr_sig_s1= inrange(s1_phe_both_xyz_kr170,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([0:1:600]',hist(s1_phe_both_xyz_kr170(cut_kr_sig_s1),[0:1:600])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_bottom_xyz_kr170,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([0:2:1000]',hist(s2_phe_bottom_xyz_kr170(cut_kr_sig_s2),[0:2:1000])','gauss1');     
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Kr_170_S1Mu=Fit_kr_s1.b1;
     Kr_170_S1Sig=Fit_kr_s1_sig;
     Kr_170_S2Mu=Fit_kr_s2.b1;
     Kr_170_S2Sig=Fit_kr_s2_sig;
     Kr_170_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_xyz_kr170(cut_kr_sig_s1))); % == sigma/sqrt(N);
     Kr_170_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_xyz_kr170(cut_kr_sig_s2))); % == sigma/sqrt(N);
     
     Kr_170_SE=mean(SE_Size_VUV);

       %% Kr83m at 42.6 keV and 100 V/cm  -- NEEDS ITS OWN XYZ Corrections (LOAD DATA GIVES NO XYZ CORRECTION)
         clear radius s2x_c s2y_c drift_time s1_phe_both_xyz s2_phe_bottom_xyz S2_numelectrons_VUV
load('C:\Program Files\MATLAB\R2012a\bin\DokeWithDQSE\Run03Kr_100Vcm_DQ_SE.mat'); %Load Patched Kr Data

n_sig=3;
     radius_c=(s2x_c.^2+s2y_c.^2).^(1/2);
     radius=(s2x.^2+s2y.^2).^(1/2);
%      Getting s1 z correction
    
dT_step=15;
dT_max=30+dT_step;
i=1;
% det_edge=351;

for dT_max=30+dT_step:dT_step:310-dT_step;
s1z_mean_fit=fit([0:3:500]',hist(s1_phe_both_xyz_VUV(inrange(drift_time,[dT_max-dT_step dT_max])),[0:3:500])','gauss1'); 
s1z_mean(i)=s1z_mean_fit.b1;
s1z_dT_loc(i)=dT_max-dT_step/2;
i=i+1;
end
mean_fit_s1z=polyfit(s1z_dT_loc,s1z_mean,2);

s1_phe_both_xyz_z=s1_phe_both_xyz_VUV.*polyval(mean_fit_s1z,171.5)./polyval(mean_fit_s1z,drift_time); %171.5 found from counting events


%Getting s1 xy correction

clear x2 hist_bin

xx_step=2;
yy_step=2;
s1_both_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(s1_phe_both_xyz_z(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step]))) > 100;
            s1_both_xymeans_fit=fit([0:3:500]',hist(s1_phe_both_xyz_z(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step])),[0:3:500])','gauss1'); 
            s1_both_xymeans(j,i)=s1_both_xymeans_fit.b1;
           else 
          
        s1_both_xymeans(j,i)=0;     
        
      end
      j=j+1;
    end
    i=i+1;
   end
   
s1_both_xymeans(isnan(s1_both_xymeans)) = 0;

s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;

    %Calculate corrections matrix and 1 sigma corrections matrix
center=interp2(s2x_bins,s2y_bins,s1_both_xymeans,0,0,'spline',1);%Normalize to the center (x=y=0)
norm_S1_both=center./s1_both_xymeans; 
norm_S1_both(isinf(norm_S1_both))=1;%no infinity! Do not correct events outside r=25cm
norm_S1_both(isnan(norm_S1_both))=1;

s1_phe_both_xyz_VUV_old=s1_phe_both_xyz_VUV;
clear s1_phe_both_xyz_VUV
s1_phe_both_xyz_VUV=s1_phe_both_xyz_z.*interp2(s2x_bins,s2y_bins,norm_S1_both,s2x,s2y,'spline');


%Getting s2 z correction
    
dT_step=15;
dT_max=30+dT_step;
i=1;
% det_edge=351;

for dT_max=30+dT_step:dT_step:310-dT_step;
s2z_mean_fit=fit([0:100:15000]',hist(s2_phe_both_xyz_VUV_SatFixed(inrange(drift_time,[dT_max-dT_step dT_max])),[0:100:15000])','gauss1'); 
s2z_mean(i)=s2z_mean_fit.b1;
s2z_dT_loc(i)=dT_max-dT_step/2;
i=i+1;
end

s2z_mean_fit=fit(s2z_dT_loc.',s2z_mean.','exp1');
electron_lifetime=-1/s2z_mean_fit.b;
s2_phe_both_z=s2_phe_both_xyz_VUV_SatFixed.*exp(drift_time./electron_lifetime);

%Getting s2 xy correction

clear x2 hist_bin

xx_step=2;
yy_step=2;
s2_bot_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(s2_phe_both_z(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step]))) > 300;
            s2_bot_xymeans_fit=fit([0:200:20000]',hist(s2_phe_both_z(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step])),[0:200:20000])','gauss1'); 
            s2_bot_xymeans(j,i)=s2_bot_xymeans_fit.b1;
           else 
          
        s2_bot_xymeans(j,i)=0;     
        
      end
      j=j+1;
    end
    i=i+1;
   end
   
s2_bot_xymeans(isnan(s2_bot_xymeans)) = 0;

s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;

    %Calculate corrections matrix and 1 sigma corrections matrix
center=interp2(s2x_bins,s2y_bins,s2_bot_xymeans,0,0,'spline',1);%Normalize to the center (x=y=0)
norm_S2_bot=center./s2_bot_xymeans; 
norm_S2_bot(isinf(norm_S2_bot))=1;%no infinity! Do not correct events outside r=25cm
norm_S2_bot(isnan(norm_S2_bot))=1;

s2_phe_both_xyz_VUV_SatFixed_old=s2_phe_both_xyz_VUV_SatFixed;
clear s2_phe_both_xyz_VUV_SatFixed
s2_phe_both_xyz_VUV_SatFixed=s2_phe_both_z.*interp2(s2x_bins,s2y_bins,norm_S2_bot,s2x,s2y,'spline');


%%%%%%
     S2_numelectrons_VUV=s2_phe_both_xyz_VUV_SatFixed./SE_Size_VUV;
     s2_phe_both_xyz_kr100_phe=s2_phe_both_xyz_VUV_SatFixed(inrange(drift_time,[30 300]) & inrange(radius_c,[0 18]) & inrange(s2_phe_both_xyz_VUV_SatFixed,[1585,20000]) & inrange(s1_phe_both_xyz_VUV,[150,400]));
     s1_phe_both_xyz_kr100_phe=s1_phe_both_xyz_VUV(inrange(drift_time,[30 300]) & inrange(radius_c,[0 18]) & inrange(s2_phe_both_xyz_VUV_SatFixed,[1585,20000]) & inrange(s1_phe_both_xyz_VUV,[150,400]));
     
     s1_phe_both_xyz_kr100=s1_phe_both_xyz_VUV(inrange(drift_time,[30 300]) & inrange(radius_c,[0 18]) & s2_phe_both_xyz_VUV_SatFixed > 1585);
     s2_phe_both_xyz_kr100=S2_numelectrons_VUV(inrange(drift_time,[30 300]) & inrange(radius_c,[0 18]) & s2_phe_both_xyz_VUV_SatFixed > 1585);
    
     
     
          clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2   
     Fit_kr_s1=fit([0:1:600]',hist(s1_phe_both_xyz_kr100(s1_phe_both_xyz_kr100<600),[0:1:600])','gauss1');
     Fit_kr_s2=fit([0:2:1000]',hist(s2_phe_both_xyz_kr100(s2_phe_both_xyz_kr100<1000),[0:2:1000])','gauss1');  
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_xyz_kr100)); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_both_xyz_kr100)); % == sigma/sqrt(N);
     %cut and fit around 3 sigma of the mean, after the first guess 
     cut_kr_sig_s1= inrange(s1_phe_both_xyz_kr100,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([0:1:600]',hist(s1_phe_both_xyz_kr100(cut_kr_sig_s1),[0:1:600])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_both_xyz_kr100,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([0:2:1000]',hist(s2_phe_both_xyz_kr100(cut_kr_sig_s2),[0:2:1000])','gauss1');        
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Kr_100_S1Mu=Fit_kr_s1.b1;
     Kr_100_S1Sig=Fit_kr_s1_sig;
     Kr_100_S2Mu=Fit_kr_s2.b1;
     Kr_100_S2Sig=Fit_kr_s2_sig;
     Kr_100_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_xyz_kr100(cut_kr_sig_s1))); % == sigma/sqrt(N);
     Kr_100_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_both_xyz_kr100(cut_kr_sig_s2))); % == sigma/sqrt(N);
         
     Kr_100_SE=mean(SE_Size_VUV);
%% Load Kr 50 V/cm
         clear radius radius_c s2x_c s2y_c drift_time s1_phe_both_xyz s2_phe_bottom_xyz S2_numelectrons_VUV
load('C:\Program Files\MATLAB\R2012a\bin\DokeWithDQSE\Run03Kr_50Vcm_DQ_SE.mat'); %Load Patched Kr Data

n_sig=3;
     radius=(s2x_c.^2+s2y_c.^2).^(1/2);
     
%      Getting s1 z correction
    
dT_step=15;
dT_max=30+dT_step;
i=1;
% det_edge=463;

for dT_max=30+dT_step:dT_step:310-dT_step;
s1z_mean_fit=fit([0:3:500]',hist(s1_phe_both_xyz_VUV(inrange(drift_time,[dT_max-dT_step dT_max])),[0:3:500])','gauss1'); 
s1z_mean(i)=s1z_mean_fit.b1;
s1z_dT_loc(i)=dT_max-dT_step/2;
i=i+1;
end
mean_fit_s1z=polyfit(s1z_dT_loc,s1z_mean,2);

s1_phe_both_xyz_z=s1_phe_both_xyz_VUV.*polyval(mean_fit_s1z,184.5)./polyval(mean_fit_s1z,drift_time);


%Getting s1 xy correction

clear x2 hist_bin

xx_step=2;
yy_step=2;
s1_both_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(s1_phe_both_xyz_z(s1_phe_both_xyz_z<500 & inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step]))) > 100;
            s1_both_xymeans_fit=fit([0:3:500]',hist(s1_phe_both_xyz_z(s1_phe_both_xyz_z<500 & inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step])),[0:3:500])','gauss1'); 
            s1_both_xymeans(j,i)=s1_both_xymeans_fit.b1;
           else 
          
        s1_both_xymeans(j,i)=0;     
        
      end
      j=j+1;
    end
    i=i+1;
   end
   
s1_both_xymeans(isnan(s1_both_xymeans)) = 0;

s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;

    %Calculate corrections matrix and 1 sigma corrections matrix
center=interp2(s2x_bins,s2y_bins,s1_both_xymeans,0,0,'spline',1);%Normalize to the center (x=y=0)
norm_S1_both=center./s1_both_xymeans; 
norm_S1_both(isinf(norm_S1_both))=1;%no infinity! Do not correct events outside r=25cm
norm_S1_both(isnan(norm_S1_both))=1;

s1_phe_both_xyz_VUV_old=s1_phe_both_xyz_VUV;
clear s1_phe_both_xyz_VUV
s1_phe_both_xyz_VUV=s1_phe_both_xyz_z.*interp2(s2x_bins,s2y_bins,norm_S1_both,s2x,s2y,'spline');


%Getting s2 z correction
    
dT_step=15;
dT_max=30+dT_step;
i=1;
% det_edge=463;

for dT_max=30+dT_step:dT_step:310-dT_step;
s2z_mean_fit=fit([0:200:20000]',hist(s2_phe_both_xyz_VUV_SatFixed(inrange(drift_time,[dT_max-dT_step dT_max])),[0:200:20000])','gauss1'); 
s2z_mean(i)=s2z_mean_fit.b1;
s2z_dT_loc(i)=dT_max-dT_step/2;
i=i+1;
end

s2z_mean_fit=fit(s2z_dT_loc.',s2z_mean.','exp1');
electron_lifetime=-1/s2z_mean_fit.b;
s2_phe_both_z=s2_phe_both_xyz_VUV_SatFixed.*exp(drift_time./electron_lifetime);

%Getting s2 xy correction

clear x2 hist_bin

xx_step=2;
yy_step=2;
s2_bot_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(s2_phe_both_z(s2_phe_both_z<20000 & inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step]))) > 100;
            s2_bot_xymeans_fit=fit([0:200:20000]',hist(s2_phe_both_z(s2_phe_both_z<20000 & inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step])),[0:200:20000])','gauss1'); 
            s2_bot_xymeans(j,i)=s2_bot_xymeans_fit.b1;
           else 
          
        s2_bot_xymeans(j,i)=0;     
        
      end
      j=j+1;
    end
    i=i+1;
   end
   
s2_bot_xymeans(isnan(s2_bot_xymeans)) = 0;

s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;

    %Calculate corrections matrix and 1 sigma corrections matrix
center=interp2(s2x_bins,s2y_bins,s2_bot_xymeans,0,0,'spline',1);%Normalize to the center (x=y=0)
norm_S2_bot=center./s2_bot_xymeans; 
norm_S2_bot(isinf(norm_S2_bot))=1;%no infinity! Do not correct events outside r=25cm
norm_S2_bot(isnan(norm_S2_bot))=1;

s2_phe_both_xyz_VUV_SatFixed_old=s2_phe_both_xyz_VUV_SatFixed;
clear s2_phe_both_xyz_VUV_SatFixed
s2_phe_both_xyz_VUV_SatFixed=s2_phe_both_z.*interp2(s2x_bins,s2y_bins,norm_S2_bot,s2x,s2y,'spline');


%%%%%%
     S2_numelectrons_VUV=s2_phe_both_xyz_VUV_SatFixed./SE_Size_VUV;
     
     
     s1_phe_both_xyz_kr50=s1_phe_both_xyz_VUV(inrange(drift_time,[30 400]) & inrange(radius,[0 18]) & s2_phe_both_xyz_VUV_SatFixed > 1585);
     s2_phe_bottom_xyz_kr50=S2_numelectrons_VUV(inrange(drift_time,[30 400]) & inrange(radius,[0 18]) & s2_phe_both_xyz_VUV_SatFixed > 1585); %really this is S2 both, just never renamed it
    
     
     
          clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2   
     Fit_kr_s1=fit([0:1:600]',hist(s1_phe_both_xyz_kr50(s1_phe_both_xyz_kr50<600),[0:1:600])','gauss1');
     Fit_kr_s2=fit([0:2:1000]',hist(s2_phe_bottom_xyz_kr50(s2_phe_bottom_xyz_kr50<1000),[0:2:1000])','gauss1');  
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_xyz_kr50)); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_xyz_kr50)); % == sigma/sqrt(N);
     %cut and fit around 3 sigma of the mean, after the first guess 
     cut_kr_sig_s1= inrange(s1_phe_both_xyz_kr50,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([0:1:600]',hist(s1_phe_both_xyz_kr50(cut_kr_sig_s1),[0:1:600])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_bottom_xyz_kr50,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([0:2:1000]',hist(s2_phe_bottom_xyz_kr50(cut_kr_sig_s2),[0:2:1000])','gauss1');        
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Kr_50_S1Mu=Fit_kr_s1.b1;
     Kr_50_S1Sig=Fit_kr_s1_sig;
     Kr_50_S2Mu=Fit_kr_s2.b1;
     Kr_50_S2Sig=Fit_kr_s2_sig;
     Kr_50_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_xyz_kr50(cut_kr_sig_s1))); % == sigma/sqrt(N);
     Kr_50_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_xyz_kr50(cut_kr_sig_s2))); % == sigma/sqrt(N);
         
     Kr_50_SE=mean(SE_Size_VUV);

     
     %% Cesium
     
     clear radius s2x_c s2y_c drift_time s1_phe_both_xyz s2_phe_bottom_xyz S2_numelectrons_VUV   
load('C:\Program Files\MATLAB\R2012a\bin\DokeWithDQSE\Run03_Cs_DQ_SE.mat');

n_sig=2;
     radius=(s2x_c.^2+s2y_c.^2).^(1/2);
          S2_numelectrons_VUV_SatFixed=(25.13/24.78).*S2_numelectrons_VUV_SatFixed; %converting S2/SE to be SE from Kr instead of Cs
     
     s1_phe_both_xyz_cs=s1_phe_both_xyz_VUV(inrange(drift_time,[30 300]) & inrange(radius,[0 18]));
     s2_phe_both_xyz_cs=s2_phe_both_xyz_VUV_SatFixed(inrange(drift_time,[30 300]) & inrange(radius,[0 18]));
     s2_phe_both_xyz_cs_e=S2_numelectrons_VUV_SatFixed(inrange(drift_time,[30 300]) & inrange(radius,[0 18]));        
     s2_phe_both_xyz_cs_e_bot=S2_numelectrons_bot_VUV(inrange(drift_time,[30 300]) & inrange(radius,[0 18]));        


     s1_phe_both_cut_cs=s1_phe_both_xyz_cs( ((log10(s1_phe_both_xyz_cs)+1.*log10(s2_phe_both_xyz_cs)) < 9.1) & ((log10(s1_phe_both_xyz_cs)+1.*log10(s2_phe_both_xyz_cs)) > 8.83));
     s2_phe_both_cut_cs=s2_phe_both_xyz_cs_e( ((log10(s1_phe_both_xyz_cs)+1.*log10(s2_phe_both_xyz_cs)) < 9.1) & ((log10(s1_phe_both_xyz_cs)+1.*log10(s2_phe_both_xyz_cs)) > 8.83));

     s2_phe_both_cut_cs_bot=s2_phe_both_xyz_cs_e_bot( ((log10(s1_phe_both_xyz_cs)+1.*log10(s2_phe_both_xyz_cs)) < 9.1) & ((log10(s1_phe_both_xyz_cs)+1.*log10(s2_phe_both_xyz_cs)) > 8.83));

     
     clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2   
     Fit_kr_s1=fit([0:50:6000]',hist(s1_phe_both_cut_cs,[0:50:6000])','gauss1');
     Fit_kr_s2=fit([0:300:30000]',hist(s2_phe_both_cut_cs,[0:300:30000])','gauss1');    
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_cut_cs)); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_both_cut_cs)); % == sigma/sqrt(N);    
     %cut and fit around 3 sigma of the mean, after the first guess 
     cut_kr_sig_s1= inrange(s1_phe_both_cut_cs,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([0:50:6000]',hist(s1_phe_both_cut_cs(cut_kr_sig_s1),[0:50:6000])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_both_cut_cs,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([0:300:30000]',hist(s2_phe_both_cut_cs(cut_kr_sig_s2),[0:300:30000])','gauss1');     
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Cs_S1Mu=Fit_kr_s1.b1;
     Cs_S1Sig=Fit_kr_s1_sig;
     Cs_S2Mu=Fit_kr_s2.b1;
     Cs_S2Sig=Fit_kr_s2_sig;
     Cs_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_cut_cs(cut_kr_sig_s1))); % == sigma/sqrt(N);
     Cs_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_both_cut_cs(cut_kr_sig_s2))); % == sigma/sqrt(N);
     
     Cs_SE=mean(SE_Size_VUV);
     
     %Get Old Cesium point v New cesium point CHANGE THE HISTO RANGE
     %HERE!!!!
     clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2   
     Fit_kr_s1=fit([0:50:6000]',hist(s1_phe_both_cut_cs,[0:50:6000])','gauss1');
     Fit_kr_s2=fit([0:300:30000]',hist(s2_phe_both_cut_cs_bot,[0:300:30000])','gauss1');    
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_cut_cs)); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_both_cut_cs_bot)); % == sigma/sqrt(N);    
     %cut and fit around 3 sigma of the mean, after the first guess 
     cut_kr_sig_s1= inrange(s1_phe_both_cut_cs,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([0:50:6000]',hist(s1_phe_both_cut_cs(cut_kr_sig_s1),[0:50:6000])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_both_cut_cs_bot,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([0:300:30000]',hist(s2_phe_both_cut_cs_bot(cut_kr_sig_s2),[0:300:30000])','gauss1');     
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Cs_S1Mu_bot=Fit_kr_s1.b1;
     Cs_S1Sig_bot=Fit_kr_s1_sig;
     Cs_S2Mu_bot=Fit_kr_s2.b1;
     Cs_S2Sig_bot=Fit_kr_s2_sig;
     Cs_s1_meanerr_bot = Fit_kr_s1_sig/sqrt(length(s1_phe_both_cut_cs(cut_kr_sig_s1))); % == sigma/sqrt(N);
     Cs_s2_meanerr_bot = Fit_kr_s2_sig/sqrt(length(s2_phe_both_cut_cs(cut_kr_sig_s2))); % == sigma/sqrt(N);
     
     %Extra error from total - bot difference
%      Cs_extra_err=abs(Cs_S2Mu-Cs_S2Mu_bot);

     
            %% Xe Activation Lines
         clear radius s2x_c s2y_c drift_time s1_phe_both_xyz s2_phe_bottom_xyz S2_numelectrons_VUV
load('C:\Program Files\MATLAB\R2012a\bin\DokeWithDQSE\Run03_XeAct_DQ_SE.mat'); %Load Patched Kr Data

%Converting to SE from Kr
S2_numelectrons_VUV_SatFixed=(25.81/24.66).*S2_numelectrons_VUV_SatFixed;

n_sig=3;
     radius=(s2x_c.^2+s2y_c.^2).^(1/2);
     
     s1_phe_both_xyz_act=s1_phe_both_xyz_VUV(inrange(drift_time,[30 300]) & inrange(radius,[0 18]));
     s2_phe_bottom_xyz_act=s2_phe_both_xyz_VUV_SatFixed(inrange(drift_time,[30 300]) & inrange(radius,[0 18]));
     s2_phe_bottom_xyz_act_e=S2_numelectrons_VUV_SatFixed(inrange(drift_time,[30 300]) & inrange(radius,[0 18])); %this is really s2 both, just wasnt renamed

     
     clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2 
     %Xe131m at 163.9 keV DONE
%   
     s1_phe_both_131=s1_phe_both_xyz_act( ((log10(s1_phe_both_xyz_act)+0.6*log10(s2_phe_bottom_xyz_act)) < 5.85) & ((log10(s1_phe_both_xyz_act)+1.6.*log10(s2_phe_bottom_xyz_act)) > 5.65) & s2_phe_bottom_xyz_act > 10000);
     s2_phe_bottom_131=s2_phe_bottom_xyz_act_e( ((log10(s1_phe_both_xyz_act)+0.6*log10(s2_phe_bottom_xyz_act)) < 5.85) & ((log10(s1_phe_both_xyz_act)+1.6.*log10(s2_phe_bottom_xyz_act)) > 5.65) & s2_phe_bottom_xyz_act > 10000);   
     Fit_kr_s1=fit([500:5:1500]',hist(s1_phe_both_131(inrange(s1_phe_both_131,[500 1500])),[500:5:1500])','gauss1');
     Fit_kr_s2=fit([500:30:4000]',hist(s2_phe_bottom_131(inrange(s2_phe_bottom_131,[500 4000])),[500:30:4000])','gauss1');    
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_131(inrange(s1_phe_both_131,[500 1500])))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_131(inrange(s2_phe_bottom_131,[500 4000])))); % == sigma/sqrt(N);   
     %cut and fit around 3 sigma of the mean, after the first guess 
     cut_kr_sig_s1= inrange(s1_phe_both_131,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([500:5:1500]',hist(s1_phe_both_131(cut_kr_sig_s1),[500:5:1500])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_bottom_131,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([500:30:4000]',hist(s2_phe_bottom_131(cut_kr_sig_s2),[500:30:4000])','gauss1');      
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Xe131_S1Mu=Fit_kr_s1.b1;
     Xe131_S1Sig=Fit_kr_s1_sig;
     Xe131_S2Mu=Fit_kr_s2.b1;
     Xe131_S2Sig=Fit_kr_s2_sig;
     Xe131_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_131(cut_kr_sig_s1))); % == sigma/sqrt(N);
     Xe131_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_131(cut_kr_sig_s2))); % == sigma/sqrt(N);
    
          clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2 
     % Xe127 at 207.3 keV Very shitty stats DONE
     s1_phe_both_127=s1_phe_both_xyz_act( ((log10(s1_phe_both_xyz_act)+0.6.*log10(s2_phe_bottom_xyz_act)) < 5.96) & ((log10(s1_phe_both_xyz_act)+0.6.*log10(s2_phe_bottom_xyz_act)) > 5.85)  & s2_phe_bottom_xyz_act > 10000);
     s2_phe_bottom_127=s2_phe_bottom_xyz_act_e( ((log10(s1_phe_both_xyz_act)+0.6.*log10(s2_phe_bottom_xyz_act)) < 5.96) & ((log10(s1_phe_both_xyz_act)+0.6.*log10(s2_phe_bottom_xyz_act)) > 5.85)  & s2_phe_bottom_xyz_act > 10000);        
     Fit_kr_s1=fit([600:10:1800]',hist(s1_phe_both_127(inrange(s1_phe_both_127,[600 1800])),[600:10:1800])','gauss1');
     Fit_kr_s2=fit([500:40:5500]',hist(s2_phe_bottom_127(inrange(s2_phe_bottom_127,[500 5500])),[500:40:5500])','gauss1');     
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_127(inrange(s1_phe_both_127,[600 1800])))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_127(inrange(s2_phe_bottom_127,[500 5500])))); % == sigma/sqrt(N);    
     %cut and fit around 3 sigma of the mean, after the first guess 
     cut_kr_sig_s1= inrange(s1_phe_both_127,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([600:10:1800]',hist(s1_phe_both_127(cut_kr_sig_s1),[600:10:1800])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_bottom_127,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([500:40:5500]',hist(s2_phe_bottom_127(cut_kr_sig_s2),[500:40:5500])','gauss1');          
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Xe127_S1Mu=Fit_kr_s1.b1;
     Xe127_S1Sig=Fit_kr_s1_sig;
     Xe127_S2Mu=Fit_kr_s2.b1;
     Xe127_S2Sig=Fit_kr_s2_sig;
     Xe127_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_127(cut_kr_sig_s1))); % == sigma/sqrt(N);
     Xe127_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_127(cut_kr_sig_s2))); % == sigma/sqrt(N);
     
          clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2   
     % Xe129m at 236.1 keV DONE
     s1_phe_both_129=s1_phe_both_xyz_act( ((log10(s1_phe_both_xyz_act)+0.6.*log10(s2_phe_bottom_xyz_act)) < 6.12) & ((log10(s1_phe_both_xyz_act)+0.6.*log10(s2_phe_bottom_xyz_act)) > 5.96)  & s2_phe_bottom_xyz_act > 10000 );
     s2_phe_bottom_129=s2_phe_bottom_xyz_act_e( ((log10(s1_phe_both_xyz_act)+0.6.*log10(s2_phe_bottom_xyz_act)) < 6.12) & ((log10(s1_phe_both_xyz_act)+0.6.*log10(s2_phe_bottom_xyz_act)) > 5.96) & s2_phe_bottom_xyz_act > 10000 );
     Fit_kr_s1=fit([600:10:2000]',hist(s1_phe_both_129(inrange(s1_phe_both_129,[600 2000])),[600:10:2000])','gauss1');
     Fit_kr_s2=fit([500:30:6000]',hist(s2_phe_bottom_129(inrange(s2_phe_bottom_129,[500 6000])),[500:30:6000])','gauss1');     
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_129(inrange(s1_phe_both_129,[600 2000])))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_129(inrange(s2_phe_bottom_129,[500 6000])))); % == sigma/sqrt(N);     
     %cut and fit around 3 sigma of the mean, after the first guess 
     cut_kr_sig_s1= inrange(s1_phe_both_129,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([600:10:2000]',hist(s1_phe_both_129(cut_kr_sig_s1),[600:10:2000])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_bottom_129,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([500:30:6000]',hist(s2_phe_bottom_129(cut_kr_sig_s2),[500:30:6000])','gauss1');         
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Xe129_S1Mu=Fit_kr_s1.b1;
     Xe129_S1Sig=Fit_kr_s1_sig;
     Xe129_S2Mu=Fit_kr_s2.b1;
     Xe129_S2Sig=Fit_kr_s2_sig;
     Xe129_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_129(cut_kr_sig_s1))); % == sigma/sqrt(N);
     Xe129_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_129(cut_kr_sig_s2))); % == sigma/sqrt(N);
     
     
     
          clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2 
     % Xe127 at 410 keV DONE
     s1_phe_both_127_2=s1_phe_both_xyz_act( ((log10(s1_phe_both_xyz_act)+0.6.*log10(s2_phe_bottom_xyz_act)) < 6.5) & ((log10(s1_phe_both_xyz_act)+0.6.*log10(s2_phe_bottom_xyz_act)) > 6.3) & s2_phe_bottom_xyz_act > 10000 );
     s2_phe_bottom_127_2=s2_phe_bottom_xyz_act_e( ((log10(s1_phe_both_xyz_act)+0.6.*log10(s2_phe_bottom_xyz_act)) < 6.5) & ((log10(s1_phe_both_xyz_act)+0.6.*log10(s2_phe_bottom_xyz_act)) > 6.3) & s2_phe_bottom_xyz_act > 10000 ); 
     Fit_kr_s1=fit([600:20:3500]',hist(s1_phe_both_127_2(inrange(s1_phe_both_127_2,[600 3500])),[600:20:3500])','gauss1');
     Fit_kr_s2=fit([1000:50:15000]',hist(s2_phe_bottom_127_2(inrange(s2_phe_bottom_127_2,[1000 15000])),[1000:50:15000])','gauss1');    
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_127_2(inrange(s1_phe_both_127_2,[600 3500])))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_127_2(inrange(s2_phe_bottom_127_2,[1000 15000])))); % == sigma/sqrt(N);
     %cut and fit around 3 sigma of the mean, after the first guess 
     cut_kr_sig_s1= inrange(s1_phe_both_127_2,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([600:20:3500]',hist(s1_phe_both_127_2(cut_kr_sig_s1),[600:20:3500])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_bottom_127_2,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([1000:50:15000]',hist(s2_phe_bottom_127_2(cut_kr_sig_s2),[1000:50:15000])','gauss1');          
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Xe127_2_S1Mu=Fit_kr_s1.b1;
     Xe127_2_S1Sig=Fit_kr_s1_sig;
     Xe127_2_S2Mu=Fit_kr_s2.b1;
     Xe127_2_S2Sig=Fit_kr_s2_sig;   
     Xe127_2_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_127_2(cut_kr_sig_s1))); % == sigma/sqrt(N);
     Xe127_2_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_127_2(cut_kr_sig_s2))); % == sigma/sqrt(N);

%      Bi-214 at 607 keV (shitty stats) DONE
     s1_phe_both_214=s1_phe_both_xyz_act( ((log10(s1_phe_both_xyz_act)+0.8.*log10(s2_phe_bottom_xyz_act)) < 7.85) & ((log10(s1_phe_both_xyz_act)+0.8.*log10(s2_phe_bottom_xyz_act)) > 7.65) & s2_phe_bottom_xyz_act > 10000 );
     s2_phe_bottom_214=s2_phe_bottom_xyz_act_e( ((log10(s1_phe_both_xyz_act)+0.8.*log10(s2_phe_bottom_xyz_act)) < 7.85) & ((log10(s1_phe_both_xyz_act)+0.8.*log10(s2_phe_bottom_xyz_act)) > 7.65) & s2_phe_bottom_xyz_act > 10000 );
     Fit_kr_s1=fit([0:50:5000]',hist(s1_phe_both_214(inrange(s1_phe_both_214,[0 5000])),[0:50:5000])','gauss1');
     Fit_kr_s2=fit([0:200:20000]',hist(s2_phe_bottom_214(inrange(s2_phe_bottom_214,[0 20000])),[0:200:20000])','gauss1');    
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_214(inrange(s1_phe_both_214,[0 5000])))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_214(inrange(s2_phe_bottom_214,[0 20000])))); % == sigma/sqrt(N);
     %cut and fit around 3 sigma of the mean, after the first guess 
     cut_kr_sig_s1= inrange(s1_phe_both_214,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([0:50:5000]',hist(s1_phe_both_214(cut_kr_sig_s1),[0:50:5000])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_bottom_214,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([0:200:20000]',hist(s2_phe_bottom_214(cut_kr_sig_s2),[0:200:20000])','gauss1');          
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Bi214_S1Mu=Fit_kr_s1.b1;
     Bi214_S1Sig=Fit_kr_s1_sig;
     Bi214_S2Mu=Fit_kr_s2.b1;
     Bi214_S2Sig=Fit_kr_s2_sig;   
     Bi214_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_214(cut_kr_sig_s1))); % == sigma/sqrt(N);
     Bi214_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_214(cut_kr_sig_s2))); % == sigma/sqrt(N);
     
  
      
Tot_bot_diffs = [16.2959   24.5270    0.8642  629.6990   54.5295  106.4575 138.0665  124.6813  552.9506]; %S2 Total v S2 bot differences

%Inlcuded errors:
%S1 variation over time, based on Kr170: 1.9799
%S2/SE variation over time, based on Kr170: 9.3597
%Stat error on Mean fit (_s1_meanerr and s2_meanerr)
%Error in using top v bot S2 (ToT_bot_diffs)
%Saturation error, such as (Cs_S2Mu*0.055).^2
     
     for count=1:5; %attila does 5 iterations
     %% Make Doke plot     
     S1_OverE_keV=[Kr_170_S1Mu/41.55 Kr_100_S1Mu/41.55 Kr_50_S1Mu/41.55 Cs_S1Mu/661.657 Xe131_S1Mu/163.9 Xe129_S1Mu/236.1 Xe127_S1Mu/208.3 Xe127_2_S1Mu/408.791 Bi214_S1Mu/609.312]; %Bi214_S1Mu/607 Xe127_S1Mu/207.3
     S2_OverE_keV=[Kr_170_S2Mu/41.55 Kr_100_S2Mu/41.55 Kr_50_S2Mu/41.55 Cs_S2Mu/661.657 Xe131_S2Mu/163.9 Xe129_S2Mu/236.1 Xe127_S2Mu/208.3 Xe127_2_S2Mu/408.791 Bi214_S2Mu/609.312]; %Bi214_S2Mu/607 Xe127_S2Mu/207.3
       
     S1_mean_sig=[(Kr_170_s1_meanerr.^2 + 1.9799.^2).^(1/2) (Kr_100_s1_meanerr.^2 + (1.9799*Kr_100_S1Mu/Kr_170_S1Mu).^2).^(1/2) (Kr_50_s1_meanerr.^2 + (1.9799*Kr_50_S1Mu/Kr_170_S1Mu).^2).^(1/2) (Cs_s1_meanerr.^2 + (1.9799*Cs_S1Mu/Kr_170_S1Mu).^2).^(1/2)...
         (Xe131_s1_meanerr.^2 + (1.9799*Xe131_S1Mu/Kr_170_S1Mu).^2).^(1/2) (Xe129_s1_meanerr.^2 + (1.9799*Xe129_S1Mu/Kr_170_S1Mu).^2).^(1/2) (Xe127_s1_meanerr.^2 + (1.9799*Xe127_S1Mu/Kr_170_S1Mu).^2).^(1/2)...
          (Xe127_2_s1_meanerr + (1.9799*Xe127_2_S1Mu/Kr_170_S1Mu).^2).^(1/2) (Bi214_s1_meanerr + (1.9799*Bi214_S1Mu/Kr_170_S1Mu).^2).^(1/2)]; % 
     S2_mean_sig=[((Kr_170_S2Mu*0.01).^2 + Kr_170_s2_meanerr.^2 + Tot_bot_diffs(1).^2 + 9.3597.^2).^(1/2) ((Kr_100_S2Mu*0.01).^2 + Kr_100_s2_meanerr.^2 + Tot_bot_diffs(2).^2 + (9.3597*Kr_100_S2Mu/Kr_170_S2Mu).^2).^(1/2) ((Kr_50_S2Mu*0.01).^2 + Kr_50_s2_meanerr.^2 + Tot_bot_diffs(3).^2 + (9.3597*Kr_50_S2Mu/Kr_170_S2Mu).^2).^(1/2) ((Cs_S2Mu*0.01).^2 + (Cs_S2Mu*0.055).^2 + Cs_s2_meanerr.^2 + Tot_bot_diffs(4).^2 + (9.3597*Cs_S2Mu/Kr_170_S2Mu).^2).^(1/2)...
         ((Xe131_S2Mu*0.01).^2 + Xe131_s2_meanerr.^2 + Tot_bot_diffs(5).^2 + (9.3597*Xe131_S2Mu/Kr_170_S2Mu).^2).^(1/2) ((Xe129_S2Mu*0.01).^2 + Xe129_s2_meanerr.^2 + Tot_bot_diffs(6).^2 + (9.3597*Xe129_S2Mu/Kr_170_S2Mu).^2).^(1/2) ((Xe127_S2Mu*0.01).^2 + Xe127_s2_meanerr.^2 + Tot_bot_diffs(7).^2 + (9.3597*Xe127_S2Mu/Kr_170_S2Mu).^2).^(1/2)...
          ((Xe127_2_S2Mu*0.01).^2 + Xe127_2_s2_meanerr.^2 + Tot_bot_diffs(8).^2 + (9.3597*Xe127_2_S2Mu/Kr_170_S2Mu).^2).^(1/2) ((Bi214_S2Mu*0.01).^2 + Bi214_s2_meanerr.^2 + Tot_bot_diffs(9).^2+ (9.3597*Bi214_S2Mu/Kr_170_S2Mu).^2).^(1/2)]; %
     
     Energies=[41.55;41.55;41.55;661.657;163.9;236.1;208.3;408.791;609.312]; %in keV %607 is Bi ;207.3 is Xe127

      %Doing MCMC
     clear data model options
     data.xdata=S2_OverE_keV'; data.ydata=S1_OverE_keV'; %setting up data
     data.xdataerr=(S2_mean_sig./(Energies'))'; data.ydataerr=(S1_mean_sig./(Energies'))';
     
modelfun = @(x,theta) theta(1).*x+theta(2);
inv_sigma2 = @(theta,data) 1./( (theta(1).*data.xdataerr).^2 + data.ydataerr.^2);
ssfun = @(theta,data) 0.5.*( sum( (data.ydata - modelfun(data.xdata,theta)).^2.*inv_sigma2(theta,data) - log(inv_sigma2(theta,data)) ));

[tmin,ssmin]=fminsearch(ssfun,[-0.1;7.3],[],data); %minimizing sum of squares
     n = length(data.xdata); %number of data points
     p = 2; %number of parameters in model (slope and intercept)
     mse = ssmin/(n-p); % estimate for the error variance

     %setting up variables for MCMC code
     params = {{'theta1', tmin(1), -1} {'theta2', tmin(2), 0}};
     model.ssfun  = ssfun;
     options.nsimu = 50000;
    [res,chain,s2chain] = mcmcrun(model,data,params,options);
      
    

     g1=mean(chain(:,2))/73 % W=1/73 [keV/Phe]
     EE=-mean(chain(:,2))/mean(chain(:,1))/73

%Previous way of getting errors -- INCORRECT METHOD     
%      sig_g1=std(chain(:,2))/73
%      sig_EE=1/73 * sqrt( (std(chain(:,2))/mean(chain(:,1)))^2 + (std(chain(:,1))*mean(chain(:,2))/mean(chain(:,1))^2)^2 ) 

%Getting Errors 
    clear g1_chain EE_chain
    g1_chain=(chain(:,2))./73; % W=1/73 [keV/Phe]
    EE_chain=-chain(:,2)./chain(:,1)./73;
    g2_chain=24.66.*EE_chain;

%cov_g1_EE and cov_g1_g2 are consistent when no SE systematic is added in
cov_g1_EE=cov(g1_chain,EE_chain);
cov_g1_1g2=cov(g1_chain,1./g2_chain);
    cov_g1_1g2(2,2)=cov_g1_1g2(2,2)+(0.258/(25.2*25.2))^2; %Adding in error on SE
cov_g1_g2=cov(g1_chain,g2_chain);
    cov_g1_g2(2,2)=cov_g1_g2(2,2)+0.258^2;%Adding in error on SE
cov_m_b=cov(chain(:,1),chain(:,2));
sig_g1=sqrt(cov_g1_EE(1,1));
sig_g2=sqrt(cov_g1_g2(2,2));
sig_EE=sqrt(cov_g1_EE(2,2));

%      %Doke Plot
%      Doke_fig=figure;
%      x = linspace(0,250)';
%     modelfun = @(x,theta) theta(1)*x + theta(2); %setting up model
%     out = mcmcpred(res,chain,[],x,modelfun);
%     mcmcpredplot(out); %plots 50,90,95, and 99% bands
%     hold on;
%     %Making Doke plot
%     scatter(S2_OverE_keV, S1_OverE_keV,150,Energies,'o','fill','linewidth',1,'markeredgecolor','k');
%     xlabel('S2/E = n_e/(n_e+n_\gamma) \times g2/(W*SE) [electrons/keV] ');
%     ylabel('S1/E = n_\gamma/(n_e+n_\gamma) \times g1/W [Phe/keV] ');
% %     set(gca,'ytick',[3.8:.2:5.8]);set(gca,'xtick',[80:20:200]);
%     h_c=colorbar;caxis([0 700]); box on;
%     myfigview(16); ylabel(h_c,'Energy [keV] ','fontsize',20);    
%     for i=1:length(Energies)
%     plot([S2_OverE_keV(i)-S2_mean_sig(i)./(Energies(i)) S2_OverE_keV(i)+S2_mean_sig(i)./(Energies(i))],[S1_OverE_keV(i) S1_OverE_keV(i)],'-k','LineWidth',2)
%     plot([S2_OverE_keV(i) S2_OverE_keV(i)],[S1_OverE_keV(i)-S1_mean_sig(i)/(Energies(i)) S1_OverE_keV(i)+S1_mean_sig(i)/(Energies(i))],'-k','LineWidth',2)
%     end     
%     xlim([5 30]);ylim([3.8 7.5]);
%     set(gca,'ytick',[3.5:.25:8]);set(gca,'xtick',[0:2:30]);
%     text(16,6.5,strcat('Best fit:  \newline g1=',num2str(g1,3),'\pm',num2str(sig_g1,1),'\newline','EE=',num2str(EE,3),'\pm',num2str(sig_EE,'%1.2f')),'FontSize',16)
%  
    clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2
    
 %Redo Kr170
 n_sig_stat=4;
 n_sig=3;
 E_com_kr170=1/73.*(s1_phe_both_xyz_kr170/g1 + s2_phe_bottom_xyz_kr170/(EE));
 %Cuts on combined energy
 Fit_kr_E=fit([30:0.2:60]',hist(E_com_kr170,[30:0.2:60])','gauss1'); 
 cut_kr_sig_E= inrange(E_com_kr170,[Fit_kr_E.b1-n_sig*Fit_kr_E.c1/sqrt(2) Fit_kr_E.b1+n_sig*Fit_kr_E.c1/sqrt(2)]);
%Fits to energy cut
     Fit_kr_s1=fit([0:1:600]',hist(s1_phe_both_xyz_kr170(cut_kr_sig_E ),[0:1:600])','gauss1');
     Fit_kr_s2=fit([0:2:1000]',hist(s2_phe_bottom_xyz_kr170(cut_kr_sig_E ),[0:2:1000])','gauss1');
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_xyz_kr170(cut_kr_sig_E ))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_xyz_kr170(cut_kr_sig_E ))); % == sigma/sqrt(N);
%Fits to energy cut and sigma cut    
     cut_kr_sig_s1= inrange(s1_phe_both_xyz_kr170,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([0:1:600]',hist(s1_phe_both_xyz_kr170(cut_kr_sig_E & cut_kr_sig_s1),[0:1:600])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_bottom_xyz_kr170,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([0:2:1000]',hist(s2_phe_bottom_xyz_kr170(cut_kr_sig_E & cut_kr_sig_s2),[0:2:1000])','gauss1');              
     Kr_170_S1Mu=Fit_kr_s1.b1;
     Kr_170_S1Sig=Fit_kr_s1_sig;
     Kr_170_S2Mu=Fit_kr_s2.b1;
     Kr_170_S2Sig=Fit_kr_s2_sig;
     Kr_170_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_xyz_kr170(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Kr_170_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_xyz_kr170(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Kr170_fit_s1=Fit_kr_s1;
     Kr170_fit_s2=Fit_kr_s2;   
     Kr170_fit_E=Fit_kr_E;
     
     clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2    
 %Redo Kr100
 E_com_kr100=1/73.*(s1_phe_both_xyz_kr100/g1 + s2_phe_both_xyz_kr100/(EE));
 n_sig=3;
 %Cuts on combined energy
 Fit_kr_E=fit([30:0.2:60]',hist(E_com_kr100,[30:0.2:60])','gauss1'); 
 cut_kr_sig_E= inrange(E_com_kr100,[Fit_kr_E.b1-n_sig*Fit_kr_E.c1/sqrt(2) Fit_kr_E.b1+n_sig*Fit_kr_E.c1/sqrt(2)]);
%Fits to energy cut
     Fit_kr_s1=fit([0:1:600]',hist(s1_phe_both_xyz_kr100(cut_kr_sig_E ),[0:1:600])','gauss1');
     Fit_kr_s2=fit([0:2:1000]',hist(s2_phe_both_xyz_kr100(cut_kr_sig_E ),[0:2:1000])','gauss1');
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_xyz_kr100(cut_kr_sig_E ))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_both_xyz_kr100(cut_kr_sig_E ))); % == sigma/sqrt(N);
%Fits to energy cut and sigma cut    
     cut_kr_sig_s1= inrange(s1_phe_both_xyz_kr100,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([0:1:600]',hist(s1_phe_both_xyz_kr100(cut_kr_sig_E & cut_kr_sig_s1),[0:1:600])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_both_xyz_kr100,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([0:2:1000]',hist(s2_phe_both_xyz_kr100(cut_kr_sig_E & cut_kr_sig_s2),[0:2:1000])','gauss1');              
     Kr_100_S1Mu=Fit_kr_s1.b1;
     Kr_100_S1Sig=Fit_kr_s1_sig;
     Kr_100_S2Mu=Fit_kr_s2.b1;
     Kr_100_S2Sig=Fit_kr_s2_sig;
     Kr_100_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_xyz_kr100(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Kr_100_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_both_xyz_kr100(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Kr100_fit_s1=Fit_kr_s1;
     Kr100_fit_s2=Fit_kr_s2;
     Kr100_fit_E=Fit_kr_E;
     
         clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2    
 %Redo Kr50
 E_com_kr50=1/73.*(s1_phe_both_xyz_kr50/g1 + s2_phe_bottom_xyz_kr50/(EE));
 n_sig=3;
 %Cuts on combined energy
 Fit_kr_E=fit([30:0.2:60]',hist(E_com_kr50,[30:0.2:60])','gauss1'); 
 cut_kr_sig_E= inrange(E_com_kr50,[Fit_kr_E.b1-n_sig*Fit_kr_E.c1/sqrt(2) Fit_kr_E.b1+n_sig*Fit_kr_E.c1/sqrt(2)]);
%Fits to energy cut
     Fit_kr_s1=fit([0:1:600]',hist(s1_phe_both_xyz_kr50(cut_kr_sig_E ),[0:1:600])','gauss1');
     Fit_kr_s2=fit([0:2:1000]',hist(s2_phe_bottom_xyz_kr50(cut_kr_sig_E ),[0:2:1000])','gauss1');
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_xyz_kr50(cut_kr_sig_E ))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_xyz_kr50(cut_kr_sig_E ))); % == sigma/sqrt(N);
%Fits to energy cut and sigma cut    
     cut_kr_sig_s1= inrange(s1_phe_both_xyz_kr50,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([0:1:600]',hist(s1_phe_both_xyz_kr50(cut_kr_sig_E & cut_kr_sig_s1),[0:1:600])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_bottom_xyz_kr50,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([0:2:1000]',hist(s2_phe_bottom_xyz_kr50(cut_kr_sig_E & cut_kr_sig_s2),[0:2:1000])','gauss1');              
     Kr_50_S1Mu=Fit_kr_s1.b1;
     Kr_50_S1Sig=Fit_kr_s1_sig;
     Kr_50_S2Mu=Fit_kr_s2.b1;
     Kr_50_S2Sig=Fit_kr_s2_sig;
     Kr_50_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_xyz_kr50(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Kr_50_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_xyz_kr50(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Kr50_fit_s1=Fit_kr_s1;
     Kr50_fit_s2=Fit_kr_s2;  
     Kr50_fit_E=Fit_kr_E;
     clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2 
 
     %Redo Cs
 n_sig=2;
 E_com_cs=1/73.*(s1_phe_both_cut_cs/g1 + s2_phe_both_cut_cs/(EE));
 %Cuts on combined energy
 Fit_kr_E=fit([600:2:800]',hist(E_com_cs,[600:2:800])','gauss1'); 
 cut_kr_sig_E= inrange(E_com_cs,[Fit_kr_E.b1-n_sig*Fit_kr_E.c1/sqrt(2) Fit_kr_E.b1+n_sig*Fit_kr_E.c1/sqrt(2)]);
%Fits to energy cut
     Fit_kr_s1=fit([0:50:6000]',hist(s1_phe_both_cut_cs(cut_kr_sig_E),[0:50:6000])','gauss1');
     Fit_kr_s2=fit([0:300:30000]',hist(s2_phe_both_cut_cs(cut_kr_sig_E ),[0:300:30000])','gauss1');
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_cut_cs(cut_kr_sig_E ))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_both_cut_cs(cut_kr_sig_E ))); % == sigma/sqrt(N);
%Fits to energy cut and sigma cut    
     cut_kr_sig_s1= inrange(s1_phe_both_cut_cs,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([0:50:6000]',hist(s1_phe_both_cut_cs(cut_kr_sig_E & cut_kr_sig_s1),[0:50:6000])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_both_cut_cs,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([0:300:30000]',hist(s2_phe_both_cut_cs(cut_kr_sig_E & cut_kr_sig_s2),[0:300:30000])','gauss1');              
     Cs_S1Mu=Fit_kr_s1.b1;
     Cs_S1Sig=Fit_kr_s1_sig;
     Cs_S2Mu=Fit_kr_s2.b1;
     Cs_S2Sig=Fit_kr_s2_sig;
     Cs_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_cut_cs(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Cs_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_both_cut_cs(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Cs_fit_s1=Fit_kr_s1;
     Cs_fit_s2=Fit_kr_s2;     
     Cs_fit_E=Fit_kr_E;
     clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2      
  
%Redo Xe131m at 163.9 keV (shitty stats)
n_sig=2;
 E_com_131=1/73.*(s1_phe_both_131/g1 + s2_phe_bottom_131/(EE));
 %Cuts on combined energy
 Fit_kr_E=fit([130:1:200]',hist(E_com_131,[130:1:200])','gauss1'); 
 cut_kr_sig_E= inrange(E_com_131,[Fit_kr_E.b1-n_sig*Fit_kr_E.c1/sqrt(2) Fit_kr_E.b1+n_sig*Fit_kr_E.c1/sqrt(2)]);
%Fits to energy cut
     Fit_kr_s1=fit([500:5:1500]',hist(s1_phe_both_131(cut_kr_sig_E),[500:5:1500])','gauss1');
     Fit_kr_s2=fit([500:30:4000]',hist(s2_phe_bottom_131(cut_kr_sig_E ),[500:30:4000])','gauss1');
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_131(cut_kr_sig_E ))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_131(cut_kr_sig_E ))); % == sigma/sqrt(N);
%Fits to energy cut and sigma cut    
     cut_kr_sig_s1= inrange(s1_phe_both_131,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([500:5:1500]',hist(s1_phe_both_131(cut_kr_sig_E & cut_kr_sig_s1),[500:5:1500])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_bottom_131,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([500:30:4000]',hist(s2_phe_bottom_131(cut_kr_sig_E & cut_kr_sig_s2),[500:30:4000])','gauss1');              
     Xe131_S1Mu=Fit_kr_s1.b1;
     Xe131_S1Sig=Fit_kr_s1_sig;
     Xe131_S2Mu=Fit_kr_s2.b1;
     Xe131_S2Sig=Fit_kr_s2_sig;
     Xe131_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_131(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Xe131_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_131(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Xe131_fit_s1=Fit_kr_s1;
     Xe131_fit_s2=Fit_kr_s2;
     Xe131_fit_E=Fit_kr_E;     
     clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2    

%Redo Xe127m at 207.3
n_sig=1.4;
 E_com_127=1/73.*(s1_phe_both_127/g1 + s2_phe_bottom_127/(EE));
 %Cuts on combined energy
 Fit_kr_E=fit([160:1:250]',hist(E_com_127,[160:1:250])','gauss1'); 
 cut_kr_sig_E= inrange(E_com_127,[Fit_kr_E.b1-n_sig*Fit_kr_E.c1/sqrt(2) Fit_kr_E.b1+n_sig*Fit_kr_E.c1/sqrt(2)]);
%Fits to energy cut
     Fit_kr_s1=fit([600:10:1800]',hist(s1_phe_both_127(cut_kr_sig_E),[600:10:1800])','gauss1');
     Fit_kr_s2=fit([500:30:5500]',hist(s2_phe_bottom_127(cut_kr_sig_E ),[500:30:5500])','gauss1');
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_127(cut_kr_sig_E ))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_127(cut_kr_sig_E ))); % == sigma/sqrt(N);
%Fits to energy cut and sigma cut    
     cut_kr_sig_s1= inrange(s1_phe_both_127,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([600:10:1800]',hist(s1_phe_both_127(cut_kr_sig_E & cut_kr_sig_s1),[600:10:1800])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_bottom_127,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([500:40:5500]',hist(s2_phe_bottom_127(cut_kr_sig_E & cut_kr_sig_s2),[500:40:5500])','gauss1');              
     Xe127_S1Mu=Fit_kr_s1.b1;
     Xe127_S1Sig=Fit_kr_s1_sig;
     Xe127_S2Mu=Fit_kr_s2.b1;
     Xe127_S2Sig=Fit_kr_s2_sig;
     Xe127_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_127(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Xe127_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_127(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Xe127_fit_s1=Fit_kr_s1;
     Xe127_fit_s2=Fit_kr_s2;
     Xe127_fit_E=Fit_kr_E;     
     clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2    

% Redo Xe129 at 236.1
n_sig=1.5; 
 E_com_129=1/73.*(s1_phe_both_129/g1 + s2_phe_bottom_129/(EE));
 %Cuts on combined energy
 Fit_kr_E=fit([200:1:275]',hist(E_com_129,[200:1:275])','gauss1'); 
 cut_kr_sig_E= inrange(E_com_129,[Fit_kr_E.b1-n_sig*Fit_kr_E.c1/sqrt(2) Fit_kr_E.b1+n_sig*Fit_kr_E.c1/sqrt(2)]);
%Fits to energy cut
     Fit_kr_s1=fit([600:10:2000]',hist(s1_phe_both_129(cut_kr_sig_E),[600:10:2000])','gauss1');
     Fit_kr_s2=fit([500:30:6000]',hist(s2_phe_bottom_129(cut_kr_sig_E ),[500:30:6000])','gauss1');
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_129(cut_kr_sig_E ))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_129(cut_kr_sig_E ))); % == sigma/sqrt(N);
%Fits to energy cut and sigma cut    
     cut_kr_sig_s1= inrange(s1_phe_both_129,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([600:10:2000]',hist(s1_phe_both_129(cut_kr_sig_E & cut_kr_sig_s1),[600:10:2000])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_bottom_129,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([500:30:6000]',hist(s2_phe_bottom_129(cut_kr_sig_E & cut_kr_sig_s2),[500:30:6000])','gauss1');              
     Xe129_S1Mu=Fit_kr_s1.b1;
     Xe129_S1Sig=Fit_kr_s1_sig;
     Xe129_S2Mu=Fit_kr_s2.b1;
     Xe129_S2Sig=Fit_kr_s2_sig;
     Xe129_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_129(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Xe129_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_129(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Xe129_fit_s1=Fit_kr_s1;
     Xe129_fit_s2=Fit_kr_s2;
     Xe129_fit_E=Fit_kr_E;
     clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2    
     
         
% Redo 127 at 410
  n_sig=0.7;
 E_com_127_2=1/73.*(s1_phe_both_127_2/g1 + s2_phe_bottom_127_2/(EE));
 %Cuts on combined energy
 Fit_kr_E=fit([320:1:500]',hist(E_com_127_2,[320:1:500])','gauss1'); 
 cut_kr_sig_E= inrange(E_com_127_2,[Fit_kr_E.b1-n_sig*Fit_kr_E.c1/sqrt(2) Fit_kr_E.b1+n_sig*Fit_kr_E.c1/sqrt(2)]);
%Fits to energy cut
     Fit_kr_s1=fit([600:20:3500]',hist(s1_phe_both_127_2(cut_kr_sig_E),[600:20:3500])','gauss1');
     Fit_kr_s2=fit([1000:50:15000]',hist(s2_phe_bottom_127_2(cut_kr_sig_E ),[1000:50:15000])','gauss1');
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_127_2(cut_kr_sig_E ))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_127_2(cut_kr_sig_E ))); % == sigma/sqrt(N);
%Fits to energy cut and sigma cut    
     cut_kr_sig_s1= inrange(s1_phe_both_127_2,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([600:20:3500]',hist(s1_phe_both_127_2(cut_kr_sig_E & cut_kr_sig_s1),[600:20:3500])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_bottom_127_2,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([1000:50:15000]',hist(s2_phe_bottom_127_2(cut_kr_sig_E & cut_kr_sig_s2),[1000:50:15000])','gauss1');              
     Xe127_2_S1Mu=Fit_kr_s1.b1;
     Xe127_2_S1Sig=Fit_kr_s1_sig;
     Xe127_2_S2Mu=Fit_kr_s2.b1;
     Xe127_2_S2Sig=Fit_kr_s2_sig;
     Xe127_2_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_127_2(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Xe127_2_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_127_2(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);           
     Xe127_2_fit_s1=Fit_kr_s1;
     Xe127_2_fit_s2=Fit_kr_s2;
     Xe127_2_fit_E=Fit_kr_E;
     clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2    

     % Redo Bi214 at 601     
  n_sig=1.5;
 E_com_214=1/73.*(s1_phe_both_214/g1 + s2_phe_bottom_214/EE);
  %Cuts on combined energy
 Fit_kr_E=fit([500:5:1000]',hist(E_com_214,[500:5:1000])','gauss1'); 
 cut_kr_sig_E= inrange(E_com_214,[Fit_kr_E.b1-n_sig*Fit_kr_E.c1/sqrt(2) Fit_kr_E.b1+n_sig*Fit_kr_E.c1/sqrt(2)]);
%Fits to energy cut
     Fit_kr_s1=fit([0:50:5000]',hist(s1_phe_both_214(cut_kr_sig_E),[0:50:5000])','gauss1');
     Fit_kr_s2=fit([0:200:20000]',hist(s2_phe_bottom_214(cut_kr_sig_E ),[0:200:20000])','gauss1');
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_214(cut_kr_sig_E ))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_214(cut_kr_sig_E ))); % == sigma/sqrt(N);
%Fits to energy cut and sigma cut    
     cut_kr_sig_s1= inrange(s1_phe_both_214,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([0:50:5000]',hist(s1_phe_both_214(cut_kr_sig_E & cut_kr_sig_s1),[0:50:5000])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_bottom_214,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([0:200:20000]',hist(s2_phe_bottom_214(cut_kr_sig_E & cut_kr_sig_s2),[0:200:20000])','gauss1');              
     Bi214_S1Mu=Fit_kr_s1.b1;
     Bi214_S1Sig=Fit_kr_s1_sig;
     Bi214_S2Mu=Fit_kr_s2.b1;
     Bi214_S2Sig=Fit_kr_s2_sig;
     Bi214_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_214(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Bi214_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_214(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);           
     Bi214_fit_s1=Fit_kr_s1;
     Bi214_fit_s2=Fit_kr_s2;
     Bi214_fit_E=Fit_kr_E;
     clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2    
    
     end
    
     %% Making plots
     
     %Doke Plot
     Doke_fig=figure;
     x = linspace(0,250)';
     dimc=[0.9 0.9 0.9];
     xpoints=0:2:30;
    modelfun = @(x,theta) theta(1)*x + theta(2); %setting up model
%     out = mcmcpred(res,chain,[],x,modelfun);
%     mcmcpredplot(out); %plots 50,90,95, and 99% bands
    hold on;
    %Making Doke plot
onesig_upper=(mean(chain(:,1))+std(chain(:,1))).*xpoints+(mean(chain(:,2))-std(chain(:,2)));
onesig_upper_2=(mean(chain(:,1))-std(chain(:,1))).*xpoints+(mean(chain(:,2)));
onesig_upper_3=(mean(chain(:,1))).*xpoints+(mean(chain(:,2))-std(chain(:,2)));
onesig_lower_3=(mean(chain(:,1))).*xpoints+(mean(chain(:,2))+std(chain(:,2)));
onesig_lower_2=(mean(chain(:,1))+std(chain(:,1))).*xpoints+(mean(chain(:,2)));
onesig_lower=(mean(chain(:,1))-std(chain(:,1))).*xpoints+(mean(chain(:,2))+std(chain(:,2)));
twosig_upper=(mean(chain(:,1))+2*std(chain(:,1))).*xpoints+(mean(chain(:,2))-2*std(chain(:,2)));
twosig_lower=(mean(chain(:,1))-2*std(chain(:,1))).*xpoints+(mean(chain(:,2))+2*std(chain(:,2)));
twosig_upper_2=(mean(chain(:,1))+2*std(chain(:,1))).*xpoints+(mean(chain(:,2)));
twosig_lower_2=(mean(chain(:,1))-2*std(chain(:,1))).*xpoints+(mean(chain(:,2)));
twosig_upper_3=(mean(chain(:,1))).*xpoints+(mean(chain(:,2))-2*std(chain(:,2)));
twosig_lower_3=(mean(chain(:,1))).*xpoints+(mean(chain(:,2))+2*std(chain(:,2)));

fillyy(xpoints,twosig_lower_3,twosig_upper_3,dimc)
fillyy(xpoints,twosig_lower_2,twosig_upper_2,dimc)
fillyy(xpoints,onesig_lower_2,onesig_upper_2,0.5*dimc)
fillyy(xpoints,onesig_lower_3,onesig_upper_3,0.5*dimc)

% plot(xpoints,onesig_upper,'--b','LineWidth',2)
% plot(xpoints,onesig_lower,'--b','LineWidth',2)
% plot(xpoints,twosig_upper,'--r','LineWidth',2)
% plot(xpoints,twosig_lower,'--r','LineWidth',2)
% plot(xpoints,onesig_upper_2,'--b','LineWidth',2)
% plot(xpoints,onesig_lower_2,'--b','LineWidth',2)
% plot(xpoints,twosig_upper_2,'--r','LineWidth',2)
% plot(xpoints,twosig_lower_2,'--r','LineWidth',2)
% plot(xpoints,onesig_upper_3,'--g','LineWidth',2)
% plot(xpoints,onesig_lower_3,'--g','LineWidth',2)
% plot(xpoints,twosig_upper_3,'--m','LineWidth',2)
% plot(xpoints,twosig_lower_3,'--m','LineWidth',2)

    scatter(S2_OverE_keV, S1_OverE_keV,150,Energies,'o','fill','linewidth',1,'markeredgecolor','k');
    xlabel('S2/E = n_e/(n_e+n_\gamma) \times g2/(W*SE) [electrons/keV] ');
    ylabel('S1/E = n_\gamma/(n_e+n_\gamma) \times g1/W [Phe/keV] ');
%     set(gca,'ytick',[3.8:.2:5.8]);set(gca,'xtick',[80:20:200]);
    h_c=colorbar;caxis([0 700]); box on;
    myfigview(16); ylabel(h_c,'Energy [keV] ','fontsize',20);    
    for i=1:length(Energies)
    plot([S2_OverE_keV(i)-S2_mean_sig(i)./(Energies(i)) S2_OverE_keV(i)+S2_mean_sig(i)./(Energies(i))],[S1_OverE_keV(i) S1_OverE_keV(i)],'-k','LineWidth',2)
    plot([S2_OverE_keV(i) S2_OverE_keV(i)],[S1_OverE_keV(i)-S1_mean_sig(i)/(Energies(i)) S1_OverE_keV(i)+S1_mean_sig(i)/(Energies(i))],'-k','LineWidth',2)
    end     
xlim([5 30]);ylim([3 7.5]);
xpoints=0:2:30;
ypoints=mean(chain(:,1)).*xpoints+mean(chain(:,2));
plot(xpoints,ypoints,'-k','LineWidth',2)
set(gca,'ytick',[3.25:.25:8]);set(gca,'xtick',[0:2:30]);
text(16,6.5,strcat('Best fit:  \newline g1=',num2str(g1,3),'\pm',num2str(sig_g1,1),'\newline','EE=',num2str(EE,3),'\pm',num2str(sig_EE,'%1.2f')),'FontSize',16)
 
    
% Plot theta1 and theta2 for each simulation
     figure; clf
    mcmcplot(chain,[],res,'chainpanel');myfigview(16);


%Plot theta1 v theta2 scatter to show covariance
figure; clf
densityplot(chain(:,1),chain(:,2),'rectangle',[0.001 0.01])
xlabel('Slope (S1/S2)'),ylabel('Intercept (S1/E)');myfigview(16);
c=colorbar;
ylabel(c,'Density of points');

% Plot g1 v g2
 figure; clf
densityplot(g1_chain,g2_chain,'rectangle',[0.001 0.01])
xlabel('G1'),ylabel('G2');myfigview(16);
c=colorbar;
ylabel(c,'Density of points');

%Plot g1 v 1/g2
figure; clf
densityplot(g1_chain,1./g2_chain,'rectangle',[0.001 0.001])
xlabel('G1'),ylabel('1/G2');myfigview(16);
c=colorbar;
ylabel(c,'Density of points');



%Save correlation of g1 and g2
save('G1G2Cov_DQ_KrSE','cov_m_b','cov_g1_EE','cov_g1_g2','cov_g1_1g2');
     

% figure; clf
% mcmcplot(sqrt(s2chain),[],[],'hist')
% title('Error std posterior');
% myfigview(16);

%All Energy Peaks
figure
hold on;
step(E_com_kr170,[20:0.1:70],'k')
xlabel('Energy (keV)');
ylabel('Count');
title('Kr83 170V/cm Energy Reconstruction');
line([41.55 41.55],[0 length(E_com_kr170)],'LineWidth',2)
myfigview(16);

figure
hold on;
step(E_com_kr100,[20:0.1:70],'k')
xlabel('Energy (keV)');
ylabel('Count');
title('Kr83 100V/cm Energy Reconstruction');
line([41.55 41.55],[0 length(E_com_kr100)],'LineWidth',2)
myfigview(16);

figure
hold on;
step(E_com_kr50,[20:0.1:70],'k')
xlabel('Energy (keV)');
ylabel('Count');
title('Kr83 50V/cm Energy Reconstruction');
line([41.55 41.55],[0 length(E_com_kr50)],'LineWidth',2)
myfigview(16);

figure
hold on;
step(E_com_127,[150:2:250],'k')
xlabel('Energy (keV)');
ylabel('Count');
title('Xe127 Energy Reconstruction');
line([208.3 208.3],[0 length(E_com_127)],'LineWidth',2)
myfigview(16);

figure
hold on;
step(E_com_127_2,[350:3:500],'k')
xlabel('Energy (keV)');
ylabel('Count');
title('Xe127 Energy Reconstruction');
line([408.791 408.791],[0 length(E_com_127_2)],'LineWidth',2)
myfigview(16);


figure
hold on;
step(E_com_129,[200:1:300],'k')
xlabel('Energy (keV)');
ylabel('Count');
title('Xe129 Energy Reconstruction');
line([236.1 236.1],[0 length(E_com_127_2)],'LineWidth',2)
myfigview(16);

figure
hold on;
step(E_com_131,[100:1:200],'k')
xlabel('Energy (keV)');
ylabel('Count');
title('Xe131 Energy Reconstruction');
line([163.9 163.9],[0 length(E_com_131)],'LineWidth',2)
myfigview(16);

figure
hold on;
step(E_com_214,[500:5:800],'k');
xlabel('Energy (keV)');
ylabel('Count');
title('Bi214 Energy Reconstruction');
line([609.312 609.312],[0 length(E_com_214)],'LineWidth',2)
myfigview(16);

figure
hold on;
step(E_com_cs,[600:2:800],'k');
xlabel('Energy (keV)');
ylabel('Count');
title('Cs137 Energy Reconstruction');
line([661.657 661.657],[0 length(E_com_cs)],'LineWidth',2)
myfigview(16);


%Event Selection
all_s1=[s1_phe_both_xyz_act;s1_phe_both_xyz_cs;s1_phe_both_xyz_kr100_phe];
all_s2=[s2_phe_bottom_xyz_act;s2_phe_both_xyz_cs;s2_phe_both_xyz_kr100_phe];
% 
% % all_s1_e=[s1_phe_both_127;s1_phe_both_127_2;s1_phe_both_129;s1_phe_both_131;s1_phe_both_214;s1_phe_both_cut_cs; s1_phe_both_xyz_kr100];
% % all_s2_e=[s2_phe_bottom_127;s2_phe_bottom_127_2; s2_phe_bottom_129; s2_phe_bottom_131; s2_phe_bottom_214; s2_phe_both_cut_cs; s2_phe_both_xyz_kr100];
% figure
% density(log10(all_s2),log10(all_s1),[3:0.015:6.5],[2:0.015:4.5],'binscale','log10');
% xlabel('log10(S2 phe)');
% ylabel('log10(S1 phe)')
% myfigview(16);
% 
% figure 
% hold on;
% density(log10(s2_phe_bottom_xyz_act),log10(s1_phe_both_xyz_act),[4:0.015:6.5],[2.5:0.015:4.5],'binscale','log10')
% colormapeditor
% 
% xdata=[4.1:0.1:5];
% ydata=-0.6*xdata+5.75;
% plot(xdata,ydata,'-k','Linewidth',1)
% ydata=-0.6*xdata+5.66;
% plot(xdata,ydata,'-k','Linewidth',1)
% ydata=-0.6*xdata+5.475;
% plot(xdata,ydata,'-k','Linewidth',1)
% ydata=-0.6*xdata+5.61;
% plot(xdata,ydata,'-k','Linewidth',1)
% ydata=-0.6*xdata+5.87;
% plot(xdata,ydata,'-k','Linewidth',1)
% ydata=-0.65*xdata+6.52;
% plot(xdata,ydata,'-k','Linewidth',1)
% 
% xdata=[4.5:0.1:5.3];
% ydata=-0.65*xdata+6.52;
% plot(xdata,ydata,'-k','Linewidth',1)
% ydata=-0.65*xdata+6.34;
% plot(xdata,ydata,'-k','Linewidth',1)
% 
% xdata=[4.6:0.1:5.4];
% ydata=-0.8*xdata+7.55;
% plot(xdata,ydata,'-k','Linewidth',1)
% ydata=-0.8*xdata+7.35;
% plot(xdata,ydata,'-k','Linewidth',1)
% xlabel('log10(S2 phe)');
% ylabel('log10(S1 phe)')
% myfigview(16);
% 



%% Redo fits for Rplot

 %Redo Kr170
 n_sig_stat=4;
 n_sig=3;
 E_com_kr170=1/73.*(s1_phe_both_xyz_kr170/g1 + s2_phe_bottom_xyz_kr170/(EE));
 %Cuts on combined energy
 Fit_kr_E=fit([30:0.2:60]',hist(E_com_kr170,[30:0.2:60])','gauss1'); 
 cut_kr_sig_E= inrange(E_com_kr170,[Fit_kr_E.b1-n_sig*Fit_kr_E.c1/sqrt(2) Fit_kr_E.b1+n_sig*Fit_kr_E.c1/sqrt(2)]);
%Fits to energy cut
     Fit_kr_s1=fit([0:1:600]',hist(s1_phe_both_xyz_kr170(cut_kr_sig_E ),[0:1:600])','gauss1');
     Fit_kr_s2=fit([0:2:1000]',hist(s2_phe_bottom_xyz_kr170(cut_kr_sig_E ),[0:2:1000])','gauss1');
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_xyz_kr170(cut_kr_sig_E ))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_xyz_kr170(cut_kr_sig_E ))); % == sigma/sqrt(N);
%Fits to energy cut and sigma cut    
     cut_kr_sig_s1= inrange(s1_phe_both_xyz_kr170,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([0:1:600]',hist(s1_phe_both_xyz_kr170(cut_kr_sig_E & cut_kr_sig_s1),[0:1:600])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_bottom_xyz_kr170,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([0:2:1000]',hist(s2_phe_bottom_xyz_kr170(cut_kr_sig_E & cut_kr_sig_s2),[0:2:1000])','gauss1');              
     Kr_170_S1Mu=Fit_kr_s1.b1;
     Kr_170_S1Sig=Fit_kr_s1_sig;
     Kr_170_S2Mu=Fit_kr_s2.b1;
     Kr_170_S2Sig=Fit_kr_s2_sig;
     Kr_170_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_xyz_kr170(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Kr_170_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_xyz_kr170(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Kr170_fit_s1=Fit_kr_s1;
     Kr170_fit_s2=Fit_kr_s2;   
     Kr170_fit_E=Fit_kr_E;
     
     clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2    
 %Redo Kr100
 E_com_kr100=1/73.*(s1_phe_both_xyz_kr100/g1 + s2_phe_both_xyz_kr100/(EE));
 n_sig=3;
 %Cuts on combined energy
 Fit_kr_E=fit([30:0.2:60]',hist(E_com_kr100,[30:0.2:60])','gauss1'); 
 cut_kr_sig_E= inrange(E_com_kr100,[Fit_kr_E.b1-n_sig*Fit_kr_E.c1/sqrt(2) Fit_kr_E.b1+n_sig*Fit_kr_E.c1/sqrt(2)]);
%Fits to energy cut
     Fit_kr_s1=fit([0:1:600]',hist(s1_phe_both_xyz_kr100(cut_kr_sig_E ),[0:1:600])','gauss1');
     Fit_kr_s2=fit([0:2:1000]',hist(s2_phe_both_xyz_kr100(cut_kr_sig_E ),[0:2:1000])','gauss1');
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_xyz_kr100(cut_kr_sig_E ))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_both_xyz_kr100(cut_kr_sig_E ))); % == sigma/sqrt(N);
%Fits to energy cut and sigma cut    
     cut_kr_sig_s1= inrange(s1_phe_both_xyz_kr100,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([0:1:600]',hist(s1_phe_both_xyz_kr100(cut_kr_sig_E & cut_kr_sig_s1),[0:1:600])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_both_xyz_kr100,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([0:2:1000]',hist(s2_phe_both_xyz_kr100(cut_kr_sig_E & cut_kr_sig_s2),[0:2:1000])','gauss1');              
     Kr_100_S1Mu=Fit_kr_s1.b1;
     Kr_100_S1Sig=Fit_kr_s1_sig;
     Kr_100_S2Mu=Fit_kr_s2.b1;
     Kr_100_S2Sig=Fit_kr_s2_sig;
     Kr_100_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_xyz_kr100(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Kr_100_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_both_xyz_kr100(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Kr100_fit_s1=Fit_kr_s1;
     Kr100_fit_s2=Fit_kr_s2;
     Kr100_fit_E=Fit_kr_E;
     
         clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2    
 %Redo Kr50
 E_com_kr50=1/73.*(s1_phe_both_xyz_kr50/g1 + s2_phe_bottom_xyz_kr50/(EE));
 n_sig=3;
 %Cuts on combined energy
 Fit_kr_E=fit([30:0.2:60]',hist(E_com_kr50,[30:0.2:60])','gauss1'); 
 cut_kr_sig_E= inrange(E_com_kr50,[Fit_kr_E.b1-n_sig*Fit_kr_E.c1/sqrt(2) Fit_kr_E.b1+n_sig*Fit_kr_E.c1/sqrt(2)]);
%Fits to energy cut
     Fit_kr_s1=fit([0:1:600]',hist(s1_phe_both_xyz_kr50(cut_kr_sig_E ),[0:1:600])','gauss1');
     Fit_kr_s2=fit([0:2:1000]',hist(s2_phe_bottom_xyz_kr50(cut_kr_sig_E ),[0:2:1000])','gauss1');
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_xyz_kr50(cut_kr_sig_E ))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_xyz_kr50(cut_kr_sig_E ))); % == sigma/sqrt(N);
%Fits to energy cut and sigma cut    
     cut_kr_sig_s1= inrange(s1_phe_both_xyz_kr50,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([0:1:600]',hist(s1_phe_both_xyz_kr50(cut_kr_sig_E & cut_kr_sig_s1),[0:1:600])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_bottom_xyz_kr50,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([0:2:1000]',hist(s2_phe_bottom_xyz_kr50(cut_kr_sig_E & cut_kr_sig_s2),[0:2:1000])','gauss1');              
     Kr_50_S1Mu=Fit_kr_s1.b1;
     Kr_50_S1Sig=Fit_kr_s1_sig;
     Kr_50_S2Mu=Fit_kr_s2.b1;
     Kr_50_S2Sig=Fit_kr_s2_sig;
     Kr_50_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_xyz_kr50(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Kr_50_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_xyz_kr50(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Kr50_fit_s1=Fit_kr_s1;
     Kr50_fit_s2=Fit_kr_s2;  
     Kr50_fit_E=Fit_kr_E;
     clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2 
 
     %Redo Cs
 n_sig=1.5;
 E_com_cs=1/73.*(s1_phe_both_cut_cs/g1 + s2_phe_both_cut_cs/(EE));
 %Cuts on combined energy
 Fit_kr_E=fit([600:2:800]',hist(E_com_cs,[600:2:800])','gauss1'); 
 cut_kr_sig_E= inrange(E_com_cs,[Fit_kr_E.b1-n_sig*Fit_kr_E.c1/sqrt(2) Fit_kr_E.b1+n_sig*Fit_kr_E.c1/sqrt(2)]);
%Fits to energy cut
     Fit_kr_s1=fit([0:50:6000]',hist(s1_phe_both_cut_cs(cut_kr_sig_E),[0:50:6000])','gauss1');
     Fit_kr_s2=fit([0:300:30000]',hist(s2_phe_both_cut_cs(cut_kr_sig_E ),[0:300:30000])','gauss1');
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_cut_cs(cut_kr_sig_E ))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_both_cut_cs(cut_kr_sig_E ))); % == sigma/sqrt(N);
%Fits to energy cut and sigma cut    
     cut_kr_sig_s1= inrange(s1_phe_both_cut_cs,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([0:50:6000]',hist(s1_phe_both_cut_cs(cut_kr_sig_E & cut_kr_sig_s1),[0:50:6000])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_both_cut_cs,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([0:300:30000]',hist(s2_phe_both_cut_cs(cut_kr_sig_E & cut_kr_sig_s2),[0:300:30000])','gauss1');              
     Cs_S1Mu=Fit_kr_s1.b1;
     Cs_S1Sig=Fit_kr_s1_sig;
     Cs_S2Mu=Fit_kr_s2.b1;
     Cs_S2Sig=Fit_kr_s2_sig;
     Cs_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_cut_cs(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Cs_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_both_cut_cs(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Cs_fit_s1=Fit_kr_s1;
     Cs_fit_s2=Fit_kr_s2;     
     Cs_fit_E=Fit_kr_E;
     clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2      
  
%Redo Xe131m at 163.9 keV (shitty stats)
n_sig=1;
 E_com_131=1/73.*(s1_phe_both_131/g1 + s2_phe_bottom_131/(EE));
 %Cuts on combined energy
 Fit_kr_E=fit([130:1:200]',hist(E_com_131,[130:1:200])','gauss1'); 
 cut_kr_sig_E= inrange(E_com_131,[Fit_kr_E.b1-n_sig*Fit_kr_E.c1/sqrt(2) Fit_kr_E.b1+n_sig*Fit_kr_E.c1/sqrt(2)]);
%Fits to energy cut
     Fit_kr_s1=fit([500:5:1500]',hist(s1_phe_both_131(cut_kr_sig_E),[500:5:1500])','gauss1');
     Fit_kr_s2=fit([500:30:4000]',hist(s2_phe_bottom_131(cut_kr_sig_E ),[500:30:4000])','gauss1');
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_131(cut_kr_sig_E ))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_131(cut_kr_sig_E ))); % == sigma/sqrt(N);
%Fits to energy cut and sigma cut    
     cut_kr_sig_s1= inrange(s1_phe_both_131,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([500:5:1500]',hist(s1_phe_both_131(cut_kr_sig_E & cut_kr_sig_s1),[500:5:1500])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_bottom_131,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([500:30:4000]',hist(s2_phe_bottom_131(cut_kr_sig_E & cut_kr_sig_s2),[500:30:4000])','gauss1');              
     Xe131_S1Mu=Fit_kr_s1.b1;
     Xe131_S1Sig=Fit_kr_s1_sig;
     Xe131_S2Mu=Fit_kr_s2.b1;
     Xe131_S2Sig=Fit_kr_s2_sig;
     Xe131_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_131(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Xe131_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_131(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Xe131_fit_s1=Fit_kr_s1;
     Xe131_fit_s2=Fit_kr_s2;
     Xe131_fit_E=Fit_kr_E;     
     clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2    

%Redo Xe127m at 207.3
n_sig=1;
 E_com_127=1/73.*(s1_phe_both_127/g1 + s2_phe_bottom_127/(EE));
 %Cuts on combined energy
 Fit_kr_E=fit([160:1:250]',hist(E_com_127,[160:1:250])','gauss1'); 
 cut_kr_sig_E= inrange(E_com_127,[Fit_kr_E.b1-n_sig*Fit_kr_E.c1/sqrt(2) Fit_kr_E.b1+n_sig*Fit_kr_E.c1/sqrt(2)]);
%Fits to energy cut
     Fit_kr_s1=fit([600:10:1800]',hist(s1_phe_both_127(cut_kr_sig_E),[600:10:1800])','gauss1');
     Fit_kr_s2=fit([500:30:5500]',hist(s2_phe_bottom_127(cut_kr_sig_E ),[500:30:5500])','gauss1');
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_127(cut_kr_sig_E ))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_127(cut_kr_sig_E ))); % == sigma/sqrt(N);
%Fits to energy cut and sigma cut    
     cut_kr_sig_s1= inrange(s1_phe_both_127,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([600:10:1800]',hist(s1_phe_both_127(cut_kr_sig_E & cut_kr_sig_s1),[600:10:1800])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_bottom_127,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([500:40:5500]',hist(s2_phe_bottom_127(cut_kr_sig_E & cut_kr_sig_s2),[500:40:5500])','gauss1');              
     Xe127_S1Mu=Fit_kr_s1.b1;
     Xe127_S1Sig=Fit_kr_s1_sig;
     Xe127_S2Mu=Fit_kr_s2.b1;
     Xe127_S2Sig=Fit_kr_s2_sig;
     Xe127_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_127(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Xe127_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_127(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Xe127_fit_s1=Fit_kr_s1;
     Xe127_fit_s2=Fit_kr_s2;
     Xe127_fit_E=Fit_kr_E;     
     clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2    

% Redo Xe129 at 236.1
n_sig=1; 
 E_com_129=1/73.*(s1_phe_both_129/g1 + s2_phe_bottom_129/(EE));
 %Cuts on combined energy
 Fit_kr_E=fit([200:1:275]',hist(E_com_129,[200:1:275])','gauss1'); 
 cut_kr_sig_E= inrange(E_com_129,[Fit_kr_E.b1-n_sig*Fit_kr_E.c1/sqrt(2) Fit_kr_E.b1+n_sig*Fit_kr_E.c1/sqrt(2)]);
%Fits to energy cut
     Fit_kr_s1=fit([600:10:2000]',hist(s1_phe_both_129(cut_kr_sig_E),[600:10:2000])','gauss1');
     Fit_kr_s2=fit([500:30:6000]',hist(s2_phe_bottom_129(cut_kr_sig_E ),[500:30:6000])','gauss1');
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_129(cut_kr_sig_E ))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_129(cut_kr_sig_E ))); % == sigma/sqrt(N);
%Fits to energy cut and sigma cut    
     cut_kr_sig_s1= inrange(s1_phe_both_129,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([600:10:2000]',hist(s1_phe_both_129(cut_kr_sig_E & cut_kr_sig_s1),[600:10:2000])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_bottom_129,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([500:30:6000]',hist(s2_phe_bottom_129(cut_kr_sig_E & cut_kr_sig_s2),[500:30:6000])','gauss1');              
     Xe129_S1Mu=Fit_kr_s1.b1;
     Xe129_S1Sig=Fit_kr_s1_sig;
     Xe129_S2Mu=Fit_kr_s2.b1;
     Xe129_S2Sig=Fit_kr_s2_sig;
     Xe129_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_129(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Xe129_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_129(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Xe129_fit_s1=Fit_kr_s1;
     Xe129_fit_s2=Fit_kr_s2;
     Xe129_fit_E=Fit_kr_E;
     clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2    
     
         
% Redo 127 at 410
  n_sig=1.5;
 E_com_127_2=1/73.*(s1_phe_both_127_2/g1 + s2_phe_bottom_127_2/(EE));
 %Cuts on combined energy
 Fit_kr_E=fit([320:1:500]',hist(E_com_127_2,[320:1:500])','gauss1'); 
 cut_kr_sig_E= inrange(E_com_127_2,[Fit_kr_E.b1-n_sig*Fit_kr_E.c1/sqrt(2) Fit_kr_E.b1+n_sig*Fit_kr_E.c1/sqrt(2)]);
%Fits to energy cut
     Fit_kr_s1=fit([600:20:3500]',hist(s1_phe_both_127_2(cut_kr_sig_E),[600:20:3500])','gauss1');
     Fit_kr_s2=fit([1000:50:15000]',hist(s2_phe_bottom_127_2(cut_kr_sig_E ),[1000:50:15000])','gauss1');
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_127_2(cut_kr_sig_E ))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_127_2(cut_kr_sig_E ))); % == sigma/sqrt(N);
%Fits to energy cut and sigma cut    
     cut_kr_sig_s1= inrange(s1_phe_both_127_2,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([600:20:3500]',hist(s1_phe_both_127_2(cut_kr_sig_E & cut_kr_sig_s1),[600:20:3500])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_bottom_127_2,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([1000:50:15000]',hist(s2_phe_bottom_127_2(cut_kr_sig_E & cut_kr_sig_s2),[1000:50:15000])','gauss1');              
     Xe127_2_S1Mu=Fit_kr_s1.b1;
     Xe127_2_S1Sig=Fit_kr_s1_sig;
     Xe127_2_S2Mu=Fit_kr_s2.b1;
     Xe127_2_S2Sig=Fit_kr_s2_sig;
     Xe127_2_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_127_2(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Xe127_2_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_127_2(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);           
     Xe127_2_fit_s1=Fit_kr_s1;
     Xe127_2_fit_s2=Fit_kr_s2;
     Xe127_2_fit_E=Fit_kr_E;
     clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2    

     % Redo Bi214 at 601     
  n_sig=1.5;
 E_com_214=1/73.*(s1_phe_both_214/g1 + s2_phe_bottom_214/EE);
  %Cuts on combined energy
 Fit_kr_E=fit([500:5:1000]',hist(E_com_214,[500:5:1000])','gauss1'); 
 cut_kr_sig_E= inrange(E_com_214,[Fit_kr_E.b1-n_sig*Fit_kr_E.c1/sqrt(2) Fit_kr_E.b1+n_sig*Fit_kr_E.c1/sqrt(2)]);
%Fits to energy cut
     Fit_kr_s1=fit([0:50:5000]',hist(s1_phe_both_214(cut_kr_sig_E),[0:50:5000])','gauss1');
     Fit_kr_s2=fit([0:200:20000]',hist(s2_phe_bottom_214(cut_kr_sig_E ),[0:200:20000])','gauss1');
     Fit_kr_s1_sig = Fit_kr_s1.c1/sqrt(2);
     Fit_kr_s2_sig = Fit_kr_s2.c1/sqrt(2);
     Fit_kr_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_214(cut_kr_sig_E ))); % == sigma/sqrt(N);
     Fit_kr_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_214(cut_kr_sig_E ))); % == sigma/sqrt(N);
%Fits to energy cut and sigma cut    
     cut_kr_sig_s1= inrange(s1_phe_both_214,[Fit_kr_s1.b1-n_sig*Fit_kr_s1_sig Fit_kr_s1.b1+n_sig*Fit_kr_s1_sig]);
     Fit_kr_s1=fit([0:50:5000]',hist(s1_phe_both_214(cut_kr_sig_E & cut_kr_sig_s1),[0:50:5000])','gauss1');
     cut_kr_sig_s2= inrange(s2_phe_bottom_214,[Fit_kr_s2.b1-n_sig*Fit_kr_s2_sig Fit_kr_s2.b1+n_sig*Fit_kr_s2_sig]);
     Fit_kr_s2=fit([0:200:20000]',hist(s2_phe_bottom_214(cut_kr_sig_E & cut_kr_sig_s2),[0:200:20000])','gauss1');              
     Bi214_S1Mu=Fit_kr_s1.b1;
     Bi214_S1Sig=Fit_kr_s1_sig;
     Bi214_S2Mu=Fit_kr_s2.b1;
     Bi214_S2Sig=Fit_kr_s2_sig;
     Bi214_s1_meanerr = Fit_kr_s1_sig/sqrt(length(s1_phe_both_214(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);
     Bi214_s2_meanerr = Fit_kr_s2_sig/sqrt(length(s2_phe_bottom_214(cut_kr_sig_E & cut_kr_sig_s2))); % == sigma/sqrt(N);           
     Bi214_fit_s1=Fit_kr_s1;
     Bi214_fit_s2=Fit_kr_s2;
     Bi214_fit_E=Fit_kr_E;
     clear Fit_kr_E cut_kr_sig_E Fit_kr_E Fit_kr_s1 Fit_kr_s2 cut_kr_sig_s1 cut_kr_sig_s2 
     
     %% Doke plot LY measurements
     S1_OverE_keV=[Kr_170_S1Mu/41.55 Kr_100_S1Mu/41.55 Kr_50_S1Mu/41.55 Cs_S1Mu/661.657 Xe131_S1Mu/163.9 Xe129_S1Mu/236.1 Xe127_S1Mu/208.3 Xe127_2_S1Mu/408.791 Bi214_S1Mu/609.312]; %Bi214_S1Mu/607 Xe127_S1Mu/207.3
     S2_OverE_keV=[Kr_170_S2Mu/41.55 Kr_100_S2Mu/41.55 Kr_50_S2Mu/41.55 Cs_S2Mu/661.657 Xe131_S2Mu/163.9 Xe129_S2Mu/236.1 Xe127_S2Mu/208.3 Xe127_2_S2Mu/408.791 Bi214_S2Mu/609.312]; %Bi214_S2Mu/607 Xe127_S2Mu/207.3
       
    S1_mean_sig=[(Kr_170_s1_meanerr.^2 + 1.9799.^2).^(1/2) (Kr_100_s1_meanerr.^2 + (1.9799*Kr_100_S1Mu/Kr_170_S1Mu).^2).^(1/2) (Kr_50_s1_meanerr.^2 + (1.9799*Kr_50_S1Mu/Kr_170_S1Mu).^2).^(1/2) (Cs_s1_meanerr.^2 + (1.9799*Cs_S1Mu/Kr_170_S1Mu).^2).^(1/2)...
         (Xe131_s1_meanerr.^2 + (1.9799*Xe131_S1Mu/Kr_170_S1Mu).^2).^(1/2) (Xe129_s1_meanerr.^2 + (1.9799*Xe129_S1Mu/Kr_170_S1Mu).^2).^(1/2) (Xe127_s1_meanerr.^2 + (1.9799*Xe127_S1Mu/Kr_170_S1Mu).^2).^(1/2)...
          (Xe127_2_s1_meanerr + (1.9799*Xe127_2_S1Mu/Kr_170_S1Mu).^2).^(1/2) (Bi214_s1_meanerr + (1.9799*Bi214_S1Mu/Kr_170_S1Mu).^2).^(1/2)]; % 
    S2_mean_sig=[((Kr_170_S2Mu*0.01).^2 + Kr_170_s2_meanerr.^2 + Tot_bot_diffs(1).^2 + 9.3597.^2).^(1/2) ((Kr_100_S2Mu*0.01).^2 + Kr_100_s2_meanerr.^2 + Tot_bot_diffs(2).^2 + (9.3597*Kr_100_S2Mu/Kr_170_S2Mu).^2).^(1/2) ((Kr_50_S2Mu*0.01).^2 + Kr_50_s2_meanerr.^2 + Tot_bot_diffs(3).^2 + (9.3597*Kr_50_S2Mu/Kr_170_S2Mu).^2).^(1/2) ((Cs_S2Mu*0.01).^2 + (Cs_S2Mu*0.055).^2 + Cs_s2_meanerr.^2 + Tot_bot_diffs(4).^2 + (9.3597*Cs_S2Mu/Kr_170_S2Mu).^2).^(1/2)...
         ((Xe131_S2Mu*0.01).^2 + Xe131_s2_meanerr.^2 + Tot_bot_diffs(5).^2 + (9.3597*Xe131_S2Mu/Kr_170_S2Mu).^2).^(1/2) ((Xe129_S2Mu*0.01).^2 + Xe129_s2_meanerr.^2 + Tot_bot_diffs(6).^2 + (9.3597*Xe129_S2Mu/Kr_170_S2Mu).^2).^(1/2) ((Xe127_S2Mu*0.01).^2 + Xe127_s2_meanerr.^2 + Tot_bot_diffs(7).^2 + (9.3597*Xe127_S2Mu/Kr_170_S2Mu).^2).^(1/2)...
          ((Xe127_2_S2Mu*0.01).^2 + Xe127_2_s2_meanerr.^2 + Tot_bot_diffs(8).^2 + (9.3597*Xe127_2_S2Mu/Kr_170_S2Mu).^2).^(1/2) ((Bi214_S2Mu*0.01).^2 + Bi214_s2_meanerr.^2 + Tot_bot_diffs(9).^2+ (9.3597*Bi214_S2Mu/Kr_170_S2Mu).^2).^(1/2)]; %
     
     photonsratio_OverE_keV=(1/73).*S1_OverE_keV./0.123;
     photonsratio_OverE_keV_err=(1/73)*(((S1_mean_sig./(Energies.'))./0.123).^2+(0.003*S1_OverE_keV./(0.123^2)).^2).^(1/2); %hardcoded with the official g1 numbers and errors
  
     electronsratio_OverE_keV=(1/73).*S2_OverE_keV./(0.423);
     electronsratio_OverE_keV_err=(1/73).*(((S2_mean_sig./(Energies.'))./0.423).^2 + (0.019.*S2_OverE_keV./(0.423^2)).^2).^(1/2);
    
     %Remake doke plot in units of quanta over total quanta
     Doke_quanta_fig=figure;
     hold on;
    scatter(electronsratio_OverE_keV,photonsratio_OverE_keV,150,Energies,'o','fill','linewidth',1,'markeredgecolor','k');
    xlabel('n_e/(n_e+n_\gamma)');
    ylabel('n_\gamma/(n_e+n_\gamma)');
    h_c=colorbar;caxis([0 700]); box on;
    myfigview(16); ylabel(h_c,'Energy [keV] ','fontsize',20);    
    for i=1:length(Energies)
    plot([electronsratio_OverE_keV(i)-electronsratio_OverE_keV_err(i) electronsratio_OverE_keV(i)+electronsratio_OverE_keV_err(i)],[photonsratio_OverE_keV(i) photonsratio_OverE_keV(i)],'-k','LineWidth',2)
    plot([electronsratio_OverE_keV(i) electronsratio_OverE_keV(i)],[photonsratio_OverE_keV(i)-photonsratio_OverE_keV_err(i) photonsratio_OverE_keV(i)+photonsratio_OverE_keV_err(i)],'-k','LineWidth',2)
    end     
xlim([0 1]);ylim([0 1]);
xpoints=0:0.05:1;
yfit=polyfit(electronsratio_OverE_keV,photonsratio_OverE_keV,1);
ypoints=polyval(yfit,xpoints);
plot(xpoints,ypoints,'-k','LineWidth',2)
set(gca,'ytick',[0:0.2:1]);set(gca,'xtick',[0:0.2:1]);
text(16,6.5,strcat('Best fit:  \newline g1=',num2str(g1,3),'\pm',num2str(sig_g1,1),'\newline','EE=',num2str(EE,3),'\pm',num2str(sig_EE,'%1.2f')),'FontSize',16)
 
   %Charge yield calculation      
     photons_OverE_keV=S1_OverE_keV./0.123;
     photons_OverE_keV_err=(((S1_mean_sig./(Energies.'))./0.123).^2+(0.003*S1_OverE_keV./(0.123^2)).^2).^(1/2); %hardcoded with the official g1 numbers and errors
     LY_plot=figure;
     errorbar(Energies(4:end),photons_OverE_keV(4:end),photons_OverE_keV_err(4:end),'.k')
 
 %%  Calculate sigmaR sigmaNg and sigmaNe in units of quanta

     W=1/73;    
     S1_mean=S1_OverE_keV.*Energies.';
     S2_mean=S2_OverE_keV.*Energies.';
    %% Make Doke plot          
    sigR(1)=sqrt(1/2)*sqrt( (Kr170_fit_s1.c1/sqrt(2)/g1)^2 + (Kr170_fit_s2.c1/sqrt(2)/EE)^2 - (Kr170_fit_E.c1/sqrt(2)/W)^2 );
    sigR(2)=sqrt(1/2)*sqrt( (Kr100_fit_s1.c1/sqrt(2)/g1)^2 + (Kr100_fit_s2.c1/sqrt(2)/EE)^2 - (Kr100_fit_E.c1/sqrt(2)/W)^2 );
    sigR(3)=sqrt(1/2)*sqrt( (Kr50_fit_s1.c1/sqrt(2)/g1)^2 + (Kr50_fit_s2.c1/sqrt(2)/EE)^2 - (Kr50_fit_E.c1/sqrt(2)/W)^2 );
    sigR(4)=sqrt(1/2)*sqrt( (Cs_fit_s1.c1/sqrt(2)/g1)^2 + (Cs_fit_s2.c1/sqrt(2)/EE)^2 - (Cs_fit_E.c1/sqrt(2)/W)^2 );
    sigR(5)=sqrt(1/2)*sqrt( (Xe131_fit_s1.c1/sqrt(2)/g1)^2 + (Xe131_fit_s2.c1/sqrt(2)/EE)^2 - (Xe131_fit_E.c1/sqrt(2)/W)^2 );
    sigR(6)=sqrt(1/2)*sqrt( (Xe129_fit_s1.c1/sqrt(2)/g1)^2 + (Xe129_fit_s2.c1/sqrt(2)/EE)^2 - (Xe129_fit_E.c1/sqrt(2)/W)^2 );
    sigR(7)=sqrt(1/2)*sqrt( (Xe127_fit_s1.c1/sqrt(2)/g1)^2 + (Xe127_fit_s2.c1/sqrt(2)/EE)^2 - (Xe127_fit_E.c1/sqrt(2)/W)^2 );
    sigR(8)=sqrt(1/2)*sqrt( (Xe127_2_fit_s1.c1/sqrt(2)/g1)^2 + (Xe127_2_fit_s2.c1/sqrt(2)/EE)^2 - (Xe127_2_fit_E.c1/sqrt(2)/W)^2 );
    sigR(9)=sqrt(1/2)*sqrt( (Bi214_fit_s1.c1/sqrt(2)/g1)^2 + (Bi214_fit_s2.c1/sqrt(2)/EE)^2 - (Bi214_fit_E.c1/sqrt(2)/W)^2 );
    
    %figure out sig_sig for s1,s2, and E fits
    clear temp; temp=confint(Kr170_fit_s1,.683);
    Kr170_fit_s1_sig_sig=abs((Kr170_fit_s1.c1-temp(1,3)))/sqrt(2);
    clear temp; temp=confint(Kr170_fit_s2,.683);
    Kr170_fit_s2_sig_sig=abs((Kr170_fit_s2.c1-temp(1,3)))/sqrt(2);
    clear temp; temp=confint(Kr170_fit_E,.683);
    Kr170_fit_E_sig_sig=abs((Kr170_fit_E.c1-temp(1,3)))/sqrt(2);

    clear temp; temp=confint(Kr100_fit_s1,.683);
    Kr100_fit_s1_sig_sig=abs((Kr100_fit_s1.c1-temp(1,3)))/sqrt(2);
    clear temp; temp=confint(Kr100_fit_s2,.683);
    Kr100_fit_s2_sig_sig=abs((Kr100_fit_s2.c1-temp(1,3)))/sqrt(2);
    clear temp; temp=confint(Kr100_fit_E,.683);
    Kr100_fit_E_sig_sig=abs((Kr100_fit_E.c1-temp(1,3)))/sqrt(2);

    clear temp; temp=confint(Kr50_fit_s1,.683);
    Kr50_fit_s1_sig_sig=abs((Kr50_fit_s1.c1-temp(1,3)))/sqrt(2);
    clear temp; temp=confint(Kr50_fit_s2,.683);
    Kr50_fit_s2_sig_sig=abs((Kr50_fit_s2.c1-temp(1,3)))/sqrt(2);
    clear temp; temp=confint(Kr50_fit_E,.683);
    Kr50_fit_E_sig_sig=abs((Kr50_fit_E.c1-temp(1,3)))/sqrt(2);    
    
    clear temp; temp=confint(Cs_fit_s1,.683);
    Cs_fit_s1_sig_sig=abs((Cs_fit_s1.c1-temp(1,3)))/sqrt(2);
    clear temp; temp=confint(Cs_fit_s2,.683);
    Cs_fit_s2_sig_sig=abs((Cs_fit_s2.c1-temp(1,3)))/sqrt(2);
    clear temp; temp=confint(Cs_fit_E,.683);
    Cs_fit_E_sig_sig=abs((Cs_fit_E.c1-temp(1,3)))/sqrt(2);    

    clear temp; temp=confint(Xe131_fit_s1,.683);
    Xe131_fit_s1_sig_sig=abs((Xe131_fit_s1.c1-temp(1,3)))/sqrt(2);
    clear temp; temp=confint(Xe131_fit_s2,.683);
    Xe131_fit_s2_sig_sig=abs((Xe131_fit_s2.c1-temp(1,3)))/sqrt(2);
    clear temp; temp=confint(Xe131_fit_E,.683);
    Xe131_fit_E_sig_sig=abs((Xe131_fit_E.c1-temp(1,3)))/sqrt(2);    

    clear temp; temp=confint(Xe129_fit_s1,.683);
    Xe129_fit_s1_sig_sig=abs((Xe129_fit_s1.c1-temp(1,3)))/sqrt(2);
    clear temp; temp=confint(Xe129_fit_s2,.683);
    Xe129_fit_s2_sig_sig=abs((Xe129_fit_s2.c1-temp(1,3)))/sqrt(2);
    clear temp; temp=confint(Xe129_fit_E,.683);
    Xe129_fit_E_sig_sig=abs((Xe129_fit_E.c1-temp(1,3)))/sqrt(2);    
 
    clear temp; temp=confint(Xe127_fit_s1,.683);
    Xe127_fit_s1_sig_sig=abs((Xe127_fit_s1.c1-temp(1,3)))/sqrt(2);
    clear temp; temp=confint(Xe127_fit_s2,.683);
    Xe127_fit_s2_sig_sig=abs((Xe127_fit_s2.c1-temp(1,3)))/sqrt(2);
    clear temp; temp=confint(Xe127_fit_E,.683);
    Xe127_fit_E_sig_sig=abs((Xe127_fit_E.c1-temp(1,3)))/sqrt(2);    

    clear temp; temp=confint(Xe127_2_fit_s1,.683);
    Xe127_2_fit_s1_sig_sig=abs((Xe127_2_fit_s1.c1-temp(1,3)))/sqrt(2);
    clear temp; temp=confint(Xe127_2_fit_s2,.683);
    Xe127_2_fit_s2_sig_sig=abs((Xe127_2_fit_s2.c1-temp(1,3)))/sqrt(2);
    clear temp; temp=confint(Xe127_2_fit_E,.683);
    Xe127_2_fit_E_sig_sig=abs((Xe127_2_fit_E.c1-temp(1,3)))/sqrt(2);  
    
    clear temp; temp=confint(Bi214_fit_s1,.683);
    Bi214_fit_s1_sig_sig=abs((Bi214_fit_s1.c1-temp(1,3)))/sqrt(2);
    clear temp; temp=confint(Bi214_fit_s2,.683);
    Bi214_fit_s2_sig_sig=abs((Bi214_fit_s2.c1-temp(1,3)))/sqrt(2);
    clear temp; temp=confint(Bi214_fit_E,.683);
    Bi214_fit_E_sig_sig=abs((Bi214_fit_E.c1-temp(1,3)))/sqrt(2);        
    
    sig_sigR(1)=sqrt(1/2)*sqrt( ((Kr170_fit_s1.c1/sqrt(2)+Kr170_fit_s1_sig_sig)/g1)^2 + ((Kr170_fit_s2.c1/sqrt(2)+Kr170_fit_s2_sig_sig)/EE)^2 - ((Kr170_fit_E.c1/sqrt(2)-Kr170_fit_E_sig_sig)/W)^2 ) - sigR(1);
    sig_sigR(2)=sqrt(1/2)*sqrt( ((Kr100_fit_s1.c1/sqrt(2)+Kr100_fit_s1_sig_sig)/g1)^2 + ((Kr100_fit_s2.c1/sqrt(2)+Kr100_fit_s2_sig_sig)/EE)^2 - ((Kr100_fit_E.c1/sqrt(2)-Kr100_fit_E_sig_sig)/W)^2 ) - sigR(2);
    sig_sigR(3)=sqrt(1/2)*sqrt( ((Kr50_fit_s1.c1/sqrt(2)+Kr50_fit_s1_sig_sig)/g1)^2 + ((Kr50_fit_s2.c1/sqrt(2)+Kr50_fit_s2_sig_sig)/EE)^2 - ((Kr50_fit_E.c1/sqrt(2)-Kr50_fit_E_sig_sig)/W)^2 ) - sigR(3);
    sig_sigR(4)=sqrt(1/2)*sqrt( ((Cs_fit_s1.c1/sqrt(2)+Cs_fit_s1_sig_sig)/g1)^2 + ((Cs_fit_s2.c1/sqrt(2)+Cs_fit_s2_sig_sig)/EE)^2 - ((Cs_fit_E.c1/sqrt(2)-Cs_fit_E_sig_sig)/W)^2 ) - sigR(4);
    sig_sigR(5)=sqrt(1/2)*sqrt( ((Xe131_fit_s1.c1/sqrt(2)+Xe131_fit_s1_sig_sig)/g1)^2 + ((Xe131_fit_s2.c1/sqrt(2)+Xe131_fit_s2_sig_sig)/EE)^2 - ((Xe131_fit_E.c1/sqrt(2)-Xe131_fit_E_sig_sig)/W)^2 ) - sigR(5);
    sig_sigR(6)=sqrt(1/2)*sqrt( ((Xe129_fit_s1.c1/sqrt(2)+Xe129_fit_s1_sig_sig)/g1)^2 + ((Xe129_fit_s2.c1/sqrt(2)+Xe129_fit_s2_sig_sig)/EE)^2 - ((Xe129_fit_E.c1/sqrt(2)-Xe129_fit_E_sig_sig)/W)^2 ) - sigR(6);
    sig_sigR(7)=sqrt(1/2)*sqrt( ((Xe127_fit_s1.c1/sqrt(2)+Xe127_fit_s1_sig_sig)/g1)^2 + ((Xe127_fit_s2.c1/sqrt(2)+Xe127_fit_s2_sig_sig)/EE)^2 - ((Xe127_fit_E.c1/sqrt(2)-Xe127_fit_E_sig_sig)/W)^2 ) - sigR(7);
    sig_sigR(8)=sqrt(1/2)*sqrt( ((Xe127_2_fit_s1.c1/sqrt(2)+Xe127_2_fit_s1_sig_sig)/g1)^2 + ((Xe127_2_fit_s2.c1/sqrt(2)+Xe127_2_fit_s2_sig_sig)/EE)^2 - ((Xe127_2_fit_E.c1/sqrt(2)-Xe127_2_fit_E_sig_sig)/W)^2 ) - sigR(8);
    sig_sigR(9)=sqrt(1/2)*sqrt( ((Bi214_fit_s1.c1/sqrt(2)+Bi214_fit_s1_sig_sig)/g1)^2 + ((Bi214_fit_s2.c1/sqrt(2)+Bi214_fit_s2_sig_sig)/EE)^2 - ((Bi214_fit_E.c1/sqrt(2)-Bi214_fit_E_sig_sig)/W)^2 ) - sigR(9);
  
    s1_ng_stat(1)=sqrt( (Kr170_fit_s1.c1/sqrt(2)/g1)^2-sigR(1)^2);
    s1_ng_stat(2)=sqrt( (Kr100_fit_s1.c1/sqrt(2)/g1)^2-sigR(2)^2);
    s1_ng_stat(3)=sqrt( (Kr50_fit_s1.c1/sqrt(2)/g1)^2-sigR(3)^2);
    s1_ng_stat(4)=sqrt( (Cs_fit_s1.c1/sqrt(2)/g1)^2-sigR(4)^2);
    s1_ng_stat(5)=sqrt( (Xe131_fit_s1.c1/sqrt(2)/g1)^2-sigR(5)^2);
    s1_ng_stat(6)=sqrt( (Xe129_fit_s1.c1/sqrt(2)/g1)^2-sigR(6)^2);
    s1_ng_stat(7)=sqrt( (Xe127_fit_s1.c1/sqrt(2)/g1)^2-sigR(7)^2);
    s1_ng_stat(8)=sqrt( (Xe127_2_fit_s1.c1/sqrt(2)/g1)^2-sigR(8)^2);
    s1_ng_stat(9)=sqrt( (Bi214_fit_s1.c1/sqrt(2)/g1)^2-sigR(9)^2);

    
    sig_s1_ng_stat(1)=sqrt( ((Kr170_fit_s1.c1/sqrt(2)+Kr170_fit_s1_sig_sig)/g1)^2-sigR(1)^2) - s1_ng_stat(1);
    sig_s1_ng_stat(2)=sqrt( ((Kr100_fit_s1.c1/sqrt(2)+Kr100_fit_s1_sig_sig)/g1)^2-sigR(2)^2) - s1_ng_stat(2);
    sig_s1_ng_stat(3)=sqrt( ((Kr50_fit_s1.c1/sqrt(2)+Kr50_fit_s1_sig_sig)/g1)^2-sigR(3)^2) - s1_ng_stat(3);
    sig_s1_ng_stat(4)=sqrt( ((Cs_fit_s1.c1/sqrt(2)+Cs_fit_s1_sig_sig)/g1)^2-sigR(4)^2) - s1_ng_stat(4);
    sig_s1_ng_stat(5)=sqrt( ((Xe131_fit_s1.c1/sqrt(2)+Xe131_fit_s1_sig_sig)/g1)^2-sigR(5)^2) - s1_ng_stat(5);
    sig_s1_ng_stat(6)=sqrt( ((Xe129_fit_s1.c1/sqrt(2)+Xe129_fit_s1_sig_sig)/g1)^2-sigR(6)^2) - s1_ng_stat(6);
    sig_s1_ng_stat(7)=sqrt( ((Xe127_fit_s1.c1/sqrt(2)+Xe127_fit_s1_sig_sig)/g1)^2-sigR(7)^2) - s1_ng_stat(7);
    sig_s1_ng_stat(8)=sqrt( ((Xe127_2_fit_s1.c1/sqrt(2)+Xe127_2_fit_s1_sig_sig)/g1)^2-sigR(8)^2) - s1_ng_stat(8);
    sig_s1_ng_stat(9)=sqrt( ((Bi214_fit_s1.c1/sqrt(2)+Bi214_fit_s1_sig_sig)/g1)^2-sigR(9)^2) - s1_ng_stat(9);

    
    
    s2_ne_stat(1)=sqrt( (Kr170_fit_s2.c1/sqrt(2)/EE)^2-sigR(1)^2);
    s2_ne_stat(2)=sqrt( (Kr100_fit_s2.c1/sqrt(2)/EE)^2-sigR(2)^2);
    s2_ne_stat(3)=sqrt( (Kr50_fit_s2.c1/sqrt(2)/EE)^2-sigR(3)^2);
    s2_ne_stat(4)=sqrt( (Cs_fit_s2.c1/sqrt(2)/EE)^2-sigR(4)^2);
    s2_ne_stat(5)=sqrt( (Xe131_fit_s2.c1/sqrt(2)/EE)^2-sigR(5)^2);
    s2_ne_stat(6)=sqrt( (Xe129_fit_s2.c1/sqrt(2)/EE)^2-sigR(6)^2);
    s2_ne_stat(7)=sqrt( (Xe127_fit_s2.c1/sqrt(2)/EE)^2-sigR(7)^2);
    s2_ne_stat(8)=sqrt( (Xe127_2_fit_s2.c1/sqrt(2)/EE)^2-sigR(8)^2);
    s2_ne_stat(9)=sqrt( (Bi214_fit_s2.c1/sqrt(2)/EE)^2-sigR(9)^2);


    sig_s2_ne_stat(1)=sqrt( ((Kr170_fit_s2.c1/sqrt(2)+Kr170_fit_s2_sig_sig)/EE)^2-sigR(1)^2) - s2_ne_stat(1);
    sig_s2_ne_stat(2)=sqrt( ((Kr100_fit_s2.c1/sqrt(2)+Kr100_fit_s2_sig_sig)/EE)^2-sigR(2)^2) - s2_ne_stat(2);
    sig_s2_ne_stat(3)=sqrt( ((Kr50_fit_s2.c1/sqrt(2)+Kr50_fit_s2_sig_sig)/EE)^2-sigR(3)^2) - s2_ne_stat(3);
    sig_s2_ne_stat(4)=sqrt( ((Cs_fit_s2.c1/sqrt(2)+Cs_fit_s2_sig_sig)/EE)^2-sigR(4)^2) - s2_ne_stat(4);
    sig_s2_ne_stat(5)=sqrt( ((Xe131_fit_s2.c1/sqrt(2)+Xe131_fit_s2_sig_sig)/EE)^2-sigR(5)^2) - s2_ne_stat(5);
    sig_s2_ne_stat(6)=sqrt( ((Xe129_fit_s2.c1/sqrt(2)+Xe129_fit_s2_sig_sig)/EE)^2-sigR(6)^2) - s2_ne_stat(6);
    sig_s2_ne_stat(7)=sqrt( ((Xe127_fit_s2.c1/sqrt(2)+Xe127_fit_s2_sig_sig)/EE)^2-sigR(7)^2) - s2_ne_stat(7);
    sig_s2_ne_stat(8)=sqrt( ((Xe127_2_fit_s2.c1/sqrt(2)+Xe127_2_fit_s2_sig_sig)/EE)^2-sigR(8)^2) - s2_ne_stat(8);
    sig_s2_ne_stat(9)=sqrt( ((Bi214_fit_s2.c1/sqrt(2)+Bi214_fit_s2_sig_sig)/EE)^2-sigR(9)^2) - s2_ne_stat(9);
 %% Fit to S1_stat vs. ng, and S2_stat vs. ng. Do not fit to the Xray. Kr is ok
 %
 % Then plot. S1_stat, S2_stat, sigR vs. Energy and Quanta
 %
 
   g = fittype('a*x');
     
   Fit_ng=fit(S1_mean([1:8])'/g1,s1_ng_stat([1:8])',g,'startpoint',[0.005]);  % not fitting to the X-ray
   Fit_ne=fit(S2_mean([1:8])'/EE,s2_ne_stat([1:8])',g,'startpoint',[0.005]);
   Fit_R=fit(Energies([1:8])/W,sigR([1:8])',g,'startpoint',[0.01]);
   
   fluc_E_fig= figure;
    hold on;
     her=errorbar(Energies([1:8]),sigR([1:8]),sig_sigR([1:8]),'ok','markerfacecolor','k','markersize',8,'markeredgecolor','k');
        removeErrorBarEnds(her);
     her2=errorbar(Energies([1:8]),s1_ng_stat([1:8]),sig_s1_ng_stat([1:8]),'sb','markerfacecolor','b','markersize',8,'markeredgecolor','k');
        removeErrorBarEnds(her2);
     her3=errorbar(Energies([1:8]),s2_ne_stat([1:8]),sig_s2_ne_stat([1:8]),'^r','markerfacecolor','r','markersize',8,'markeredgecolor','k');
        removeErrorBarEnds(her3);
     box on;
     ylabel('stdev (#quanta) ');
     xlabel('Energy [keV_{ee}] ');
     myfigview(18);
     ylim([0 3500]);
     legend([her her2 her3],'R Flucs','S1(n_\gamma) Flucs','S2(n_e) Flucs','location','northwest');
           
  %////////////Quanta////////////////////////////////////////////////         
        sig_W=0.001;

    fluc_Q_fig=figure;
    hold on;
    
         x_fit=0:500:50000;
     hR=plot(x_fit,x_fit.*Fit_R.a,'--k','linewidth',2);
     h_S1=plot(x_fit,x_fit.*Fit_ng.a,'--b','linewidth',2);
     h_S2=plot(x_fit,x_fit.*Fit_ne.a,'--r','linewidth',2);
     
     
     HH=ploterr((Energies([1:8])/W),(sigR([1:8])),(sig_W*Energies([1:8])/W^2), (sig_sigR([1:8])));
          set(HH(1),'Marker','o','Color','k','LineStyle','none','markersize',8,'linewidth',1,'markerfacecolor','k','markeredgecolor','k');set(HH(2),'Color','k','linewidth',1);set(HH(3),'Color','k','linewidth',1);
     HH2=ploterr(S1_mean([1:8])/g1,s1_ng_stat([1:8]),S1_mean_sig([1:8])/g1,sig_s1_ng_stat([1:8]));
          set(HH2(1),'Marker','s','Color','b','LineStyle','none','markersize',8,'linewidth',1,'markerfacecolor','b','markeredgecolor','k');set(HH2(2),'Color','b','linewidth',1);set(HH2(3),'Color','b','linewidth',1);
     HH3=ploterr(S2_mean([1:8])/EE,s2_ne_stat([1:8]),S2_mean_sig([1:8])/EE,sig_s2_ne_stat([1:8]));  
          set(HH3(1),'Marker','^','Color','r','LineStyle','none','markersize',8,'linewidth',1,'markerfacecolor','r','markeredgecolor','k');set(HH3(2),'Color','r','linewidth',1);set(HH3(3),'Color','r','linewidth',1);

     
     ylim([0 3500]); 
     box on;
     ylabel('stdev (#quanta) ');
     xlabel('Quanta ');
     
     myfigview(18);
     
     legend([HH2(1) HH3(1) HH(1)],'S1(n_\gamma) Flucs','S2(n_e) Flucs','R Flucs','location','northwest');


%  %% Plot in Discrimination space (ne/ng) showing all calibrations, and calculate ovals of R,S1,S2 fluctuations.
%  
%       Ng_mean=S1_mean/g1;
%       Ne_mean=S2_mean./EE;
%       E_keV_fits=[Kr170_fit_E.b1 Kr100_fit_E.b1 Kr50_fit_E.b1 Cs_fit_E.b1 Xe131_fit_E.b1 Xe129_fit_E.b1 Xe127_fit_E.b1 Xe127_2_fit_E.b1 Bi214_fit_E.b1];
%         
%       sigR_logNeNg_u=log((Ne_mean+sigR)./(Ng_mean-sigR))-log(Ne_mean./Ng_mean);
%       sigR_logNeNg_l=log((Ne_mean-sigR)./(Ng_mean+sigR))-log(Ne_mean./Ng_mean);
%       
%       sigNg_logNeNg_u=log((Ne_mean)./(Ng_mean+s1_ng_stat));
%       sigNg_logNeNg_uE=E_keV_fits+s1_ng_stat*W;
%       sigNg_logNeNg_l=log((Ne_mean)./(Ng_mean-s1_ng_stat));
%       sigNg_logNeNg_lE=E_keV_fits-s1_ng_stat*W;
%       
%       sigNe_logNeNg_u=log((Ne_mean+s2_ne_stat)./(Ng_mean));
%       sigNe_logNeNg_uE=E_keV_fits+s2_ne_stat*W;
%       sigNe_logNeNg_l=log((Ne_mean-s2_ne_stat)./(Ng_mean));
%       sigNe_logNeNg_lE=E_keV_fits-s2_ne_stat*W;
%       
%      %  
%        All_disc_figure = figure;
%        hold on;
%         density(E_com(cut & E_com<800 ), logNeNg(cut & E_com<800 ), [0:5:800],[-2:.02:.5],'binscale','log10');
%         colormap(MAP);
%         myfigview(20); box on; set(gca,'linewidth',2);
%         hc=colorbar;
%         ylabel(hc,'log_{10}[Count/(n_e/n_\gamma)] ','fontsize',20);
%          xlabel('Energy [keV] '); ylabel('Count/5keV ');
%          ylabel('log(n_e/n_\gamma) ');
%          set(gca,'xtick',[0:100:800]);
%          caxis([0 3]);
%            
%     % add all fits and save plot     
%      H_Rec=errorbar( E_keV_fits,log(Ne_mean./Ng_mean),sigR_logNeNg_l, sigR_logNeNg_u ,'xk','markersize',8,'linewidth',2,'linestyle','none');
%         removeErrorBarEnds(H_Rec);
%      %h_Rec=line([E_keV_fits;E_keV_fits],[sigR_logNeNg_l;sigR_logNeNg_u]+[log(Ne_mean./Ng_mean);log(Ne_mean./Ng_mean)],'color','k','linewidth',2);
%      h_S2_stat=line([sigNe_logNeNg_lE;sigNe_logNeNg_uE],[sigNe_logNeNg_l;sigNe_logNeNg_u],'color','r','linewidth',2);
%      h_S1_stat=line([sigNg_logNeNg_lE;sigNg_logNeNg_uE],[sigNg_logNeNg_l;sigNg_logNeNg_u],'color','w','linewidth',2);
%     
%      % Plot ellipses
%      rot_s2s1=[0 0 1/100 1/100 1/100 1/200 1/300 1/500];
%      for ii=1:length(E_keV_fits);
%          
%       center_logS2S2=log(Ne_mean(ii)/Ng_mean(ii));
%       sig_E_ellipse=(sigNg_logNeNg_uE(ii)-sigNg_logNeNg_lE(ii))/2;
%      plotellipse([E_keV_fits(ii), center_logS2S2, sig_E_ellipse, sigR_logNeNg_u(ii), rot_s2s1(ii)],'-k')
%      
%      end
%          

      %% Calculate quanta and recombination fraciton

 

 Ng_mono= S1_mean/g1;
 Ne_mono= S2_mean/EE;

data_alpha=0.20;
N_ions_mono=(Ne_mono+Ng_mono)./(1+data_alpha);
N_quanta_mono=(Ne_mono+Ng_mono);
r_mono=(Ng_mono./Ne_mono-data_alpha)./(Ng_mono./Ne_mono+1);
R_mono=r_mono.*N_ions_mono;
bino_V_mono=(1-r_mono).*r_mono.*N_ions_mono;
y_mono=Ne_mono./(Ne_mono+Ng_mono);
s_mono=Ng_mono./(Ne_mono+Ng_mono);


% save('sigR_ng_ne','sigR','sig_sigR','s1_ng_stat','sig_s1_ng_stat',...
%     's2_ne_stat','sig_s2_ne_stat','Ng_mono','Ne_mono','N_ions_mono',...
%     'N_quanta_mono','r_mono','R_mono','bino_V_mono','y_mono','s_mono');

% save('G1G2_DQ_KrSE_Session');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% using: 
%   29.7 keV Xe-xray (only near the skin)
%   42.6 keV combinded Kr83m decay (be careful for time dependance)
%   662 keV Cs137
%   236 keV Xe-129 (late March, early April 2013)
%   163 keV Xe-131 (late March, early April 2013)
%   236.1, 207.3 Xe-127 (late March, early April 2013)
%   410 Xe-127
%   607 Bi-214
%
%   First make a density plot s1 vs s2 then apply a simple box cut
%       - find mean s2 and s1 for a given energy, then make Doke plot
%       -plot S1/E vs. S2/E, use intercepts to find g1,g2
%
%       - Using an R_c<20 cm cut to increase the stats for Xe_activation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


