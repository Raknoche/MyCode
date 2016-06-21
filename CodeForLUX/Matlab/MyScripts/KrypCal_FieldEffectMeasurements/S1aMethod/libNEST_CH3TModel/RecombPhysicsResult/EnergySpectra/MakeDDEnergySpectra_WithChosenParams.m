% NEED TO UPDATE THESE POLYNOMIALS AND G1 AND G2

s1ab_s2_xyz_fit_map.p1=-0.0097;
s1ab_s2_xyz_fit_map.p2=-1.7;
s1ab_s2_xyz_fit_map.p3=5.7225;


s1ab_s1_xyz_fit_map.p1=0.0025;
s1ab_s1_xyz_fit_map.p2=0.4424;
s1ab_s1_xyz_fit_map.p3=-0.2289;

inpaint_on=1;

%HAVE TO MAKE THIS .MAT FILE STILL FROM THE DD KR DATASET
load('C:\Program Files\MATLAB\R2012a\bin\Run04_ERBand\Kr2p14_FixedZHisto\MappingParametersFit_Start_Oct2015DDKr.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate 2015 DD Corrections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %If we can use 3D map, do so.  If there aren't enough events do 1D then 2D
 if length(s1ab_z)>100000; %require 100,000 events to do 3D binning
%Removing zeros from map with nearest neighbor method
s1ab_xyz_mean_old=s1ab_xyz_mean;
temp_map=s1ab_xyz_mean;
q= ~isnan(temp_map);
r=nearestpoint(find(~q),find(q));
s1ab_xyz_mean=temp_map; temp_mapQ=temp_map(q); s1ab_xyz_mean(~q)=temp_mapQ(r);s1ab_xyz_mean=reshape(s1ab_xyz_mean,size(s1ab_xyz_mean_old,1),size(s1ab_xyz_mean_old,2),size(s1ab_xyz_mean_old,3));

%interpolating with cubic so that we get NaN when extrapolating
 s1as1b_values=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1ab_xyz_mean,s2x,s2y,drift_time,'cubic');
 s1as1b_at_center=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1ab_xyz_mean,x_center,y_center,z_center,'cubic');
 s2_3Dfieldremoval= (s1ab_s2_xyz_fit_map.p1.*s1as1b_values.^2+s1ab_s2_xyz_fit_map.p2.*s1as1b_values+s1ab_s2_xyz_fit_map.p3); 
 s1_3Dfieldremoval= (s1ab_s1_xyz_fit_map.p1.*s1as1b_values.^2+s1ab_s1_xyz_fit_map.p2.*s1as1b_values+s1ab_s1_xyz_fit_map.p3); 
 
%Filling in NaN values with the s1a/s1b Z map polynomial 
  s1as1b_zvalues= polyval(s1ab_P,drift_time);
  s1as1b_z_at_center=s1as1b_at_center;   
  s2_zfieldremoval= (s1ab_s2_xyz_fit_map.p1.*s1as1b_zvalues.^2+s1ab_s2_xyz_fit_map.p2.*s1as1b_zvalues+s1ab_s2_xyz_fit_map.p3);
  s1_zfieldremoval= (s1ab_s1_xyz_fit_map.p1.*s1as1b_zvalues.^2+s1ab_s1_xyz_fit_map.p2.*s1as1b_zvalues+s1ab_s1_xyz_fit_map.p3); 
  s2_3Dfieldremoval(isnan(s2_3Dfieldremoval))=s2_zfieldremoval(isnan(s2_3Dfieldremoval));
  s1_3Dfieldremoval(isnan(s1_3Dfieldremoval))=s1_zfieldremoval(isnan(s1_3Dfieldremoval));
 
s2_phe_bottom=s2_phe_bottom./s2_3Dfieldremoval;
s2_phe_both=s2_phe_both./s2_3Dfieldremoval;
s1_phe_bottom=s1_phe_bottom./s1_3Dfieldremoval;
s1_phe_both=s1_phe_both./s1_3Dfieldremoval;


 else 
%Removing zeros from map with nearest neighbor method
s1ab_xy_mean_old=s1ab_xy_mean;
temp_map=s1ab_xy_mean;
q= ~isnan(temp_map);
r=nearestpoint(find(~q),find(q));
s1ab_xy_mean=temp_map; temp_mapQ=temp_map(q); s1ab_xy_mean(~q)=temp_mapQ(r);s1ab_xy_mean=reshape(s1ab_xy_mean,size(s1ab_xy_mean_old,1),size(s1ab_xy_mean_old,2));

%interpolating with cubic so that we get NaN when extrapolating
 s1as1b_zvalues= polyval(s1ab_P,drift_time);
 s1as1b_zcorr_xyvalues=interp2(s1ab_xbins,s1ab_ybins,s1ab_xy_mean,s2x,s2y,'cubic');
 s1as1b_z_at_center=polyval(s1ab_P,z_center); 
 s1as1b_xy_at_center=interp2(s1ab_xbins,s1ab_ybins,s1ab_xy_mean,x_center,y_center,'cubic');
 
 s2_zfieldremoval= (s1ab_s2_xyz_fit_map.p1.*s1as1b_zvalues.^2+s1ab_s2_xyz_fit_map.p2.*s1as1b_zvalues+s1ab_s2_xyz_fit_map.p3);
 s2_xyfieldremoval= (s1ab_s2_xyz_fit_map.p1.*s1as1b_zcorr_xyvalues.^2+s1ab_s2_xyz_fit_map.p2.*s1as1b_zcorr_xyvalues+s1ab_s2_xyz_fit_map.p3); 
 s1_zfieldremoval= (s1ab_s1_xyz_fit_map.p1.*s1as1b_zvalues.^2+s1ab_s1_xyz_fit_map.p2.*s1as1b_zvalues+s1ab_s1_xyz_fit_map.p3); 
 s1_xyfieldremoval= (s1ab_s1_xyz_fit_map.p1.*s1as1b_zcorr_xyvalues.^2+s1ab_s1_xyz_fit_map.p2.*s1as1b_zcorr_xyvalues+s1ab_s1_xyz_fit_map.p3); 

 %Filling in NaN values with no XY correction.  (Z dependence should never be NaN)
 s2_xyfieldremoval(isnan(s2_xyfieldremoval))=1;
 s2_zfieldremoval(isnan(s2_zfieldremoval))=1;
 s1_xyfieldremoval(isnan(s1_xyfieldremoval))=1;
 s1_zfieldremoval(isnan(s1_zfieldremoval))=1;

 
s2_phe_bottom=s2_phe_bottom./(s2_zfieldremoval.*s2_xyfieldremoval);
s2_phe_both=s2_phe_both./(s2_zfieldremoval.*s2_xyfieldremoval);
s1_phe_bottom=s1_phe_bottom./(s1_zfieldremoval.*s1_xyfieldremoval);
s1_phe_both=s1_phe_both./(s1_zfieldremoval.*s1_xyfieldremoval);



 end 

    clear y_s1_both y_s1_bottom Fit_s1_both Fit_s1_bottom mean_s1_both mean_s1_both
    clear sigma_s1_both sigma_s1_both mean_index means mean_index means_error x temp_hist
    clear hist_s1_both x xfit S
    
    %set up energy bins for the gaussian fit
    bin_s1=8;%5
    s1_bin_min=100;
    s1_bin_max=1000;
    x=s1_bin_min:bin_s1:s1_bin_max; %150:bin_s1:500;
    x=x';
        
    cut = inrange(s2radius,[rcut_min,rcut_max]) & inrange(s1_phe_both,[s1_bin_min,s1_bin_max]) & inrange(s2_width,[100,1500]) & inrange(s1_width,[30,400]);   
    
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
    
      
    for bin=1:bins

     binStart = ((bin-1).*bin_size);  binEnd = (bin.*bin_size);   

%      time_cut = inrange(drift_distance_\mus,[binStart,binEnd]) & cut==1;
       time_cut = inrange(drift_time,[binStart,binEnd]) & cut==1;
     
        hist_s1_both(:,bin)= hist(s1_phe_both(time_cut),x)'/bin_s1;
        temp_hist=hist_s1_both(:,bin);   
            amp_start=max(temp_hist);
            mean_start=sum(temp_hist.*x)/sum(temp_hist);
            sigma_start=std(temp_hist.*x);
        Fit_s1_both{bin}=fit(x(2:end-1),temp_hist(2:end-1),'gauss1');%remove bin edges

     means(bin)=Fit_s1_both{bin}.b1;
     means_error(bin) = Fit_s1_both{bin}.c1/sqrt(2)/sqrt(sum(hist_s1_both(:,bin))*bin_s1); % == sigma/sqrt(N);
     mean_index(bin) = (binStart+binEnd)/2;     

    end
    
 % //////////////////////////////Make Plots/////////////////////////////////
 
    dT_range= inrange( mean_index,[zcut_min, zcut_max] );
    dT_fit=0:1:det_edge;    
    [P_s1_both S] = polyfit(mean_index(dT_range),means(dT_range),2);
    sigma_P_s1_both = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df)';    

    s1_both_means=means;
    s1_both_means_sigma=means_error;
    s1_both_means_bin=mean_index;
    
    
%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   Calculating the S2 Z-dependence (i.e. electron lifetime)
%   Using Both PMT arrays
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     s1_z_correction=polyval(P_s1_both,241.6)./polyval(P_s1_both,drift_distance_\mus);%normalize s1 to 160 ns.
    s1_z_correction=polyval(P_s1_both,z_center)./polyval(P_s1_both,drift_time);
    
    %s1_both_z_correction=interp1(s1_both_means_bin,s1_both_means,160)./interp1(s1_both_means_bin,s1_both_means,drift_time);
    s1_phe_both_z=s1_phe_both.*s1_z_correction;

    clear y_s1_both y_s1_both Fit_s1_both Fit_EL xEL yfitEL xfit
    clear sigma_s1_both sigma_s1_both mean_index means mean_index means_error x temp_hist
    clear hist_s1_both dT_range dT_fit S x2
    
    %set up energy bins for the gaussian fit. Using both S2 only
   bin_s2=500;
   s2_bin_min=1000;
   s2_bin_max=40000;
   x2=s2_bin_min:bin_s2:s2_bin_max;
   x2=x2';
        
    cut = inrange(s2radius,[rcut_min,rcut_max]) & inrange(s2_phe_both,[s2_bin_min,s2_bin_max]) & inrange(s2_width,[100,1000]) & inrange(s1_width,[30,400]) ;
    A = sort(size(s2_phe_both(cut))); %%sort, incase 1xn or nx1 vector is defined.
    bins = ceil(A(2)./min_evts_per_bin); %get the number of points, adjust binning if stats are low.
    if bins>max_num_bins
       bins = max_num_bins; % usually we will have 65 bin. 0:5:325
    end
%     bin_size = (488/bins);%span the entire detector. 325us. typically = 5us
    bin_size = (floor(det_edge/bins*4)/4);%span the entire detector. 325us. typically = 5us. Round to nearest 0.25.
   
        hist_s2_both = zeros(length(x2),bins);
        Fit_s2_both = cell(1,bins);
        means=zeros(bins,1);
        means_error = zeros(bins,1);
        mean_index= zeros(bins,1);
        
    for bin=1:bins

     binStart = ((bin-1).*bin_size);  binEnd = (bin.*bin_size);   

%      time_cut = inrange(drift_distance_\mus,[binStart,binEnd]) & cut==1;
       time_cut = inrange(drift_time,[binStart,binEnd]) & cut==1; 
       
        hist_s2_both(:,bin)= hist(s2_phe_both(time_cut),x2)'/bin_s2;
        temp_hist=hist_s2_both(:,bin);
            amp_start=max(temp_hist);
            mean_start=sum(temp_hist.*x2)/sum(temp_hist);
            sigma_start=std(temp_hist.*x2);
        Fit_s2_both{bin}=fit(x2(2:end-1),temp_hist(2:end-1),'gauss1');%remove bin edges
          

         means(bin)=Fit_s2_both{bin}.b1;
         means_error(bin) = Fit_s2_both{bin}.c1/sqrt(2)/sqrt(sum(hist_s2_both(:,bin))*bin_s2); % == sigma/sqrt(N);
         mean_index(bin) = (binStart+binEnd)/2;     
     
    end
    
 % //////////////////////////////Make Plots/////////////////////////////////
 
  dT_range= inrange( mean_index,[5+bin_size, 500] );%from 50 to 300us for the fit, based on when s1as1b values get strange

    s2_both_norm_z_index=mean_index(dT_range);
    s2_both_norm_z_means=means(dT_range);
    s2_at_top=RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,4); %normalize to right below the gate
    s2_both_norm_z=RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,4)./s2_both_norm_z_means;

    s2_both_means=means;
    s2_both_means_sigma=means_error;
    s2_both_means_bin=mean_index;      

    %%Tracker for pseudo-lifetime
    s2_at_bottom=RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,320);
    pseudo_lifetime_2015=-320/log(s2_at_bottom/s2_at_top);
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Mapping S2 in XY (after correcting for Z-dep). Both PMT
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s2_phe_both_z=s2_phe_both.*RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z, drift_time);

    clear x2 hist_bin

    %set up energy bins for the gaussian fit. Using both S2 only
   bin_s2=500;
   s2_bin_min=1000;
   s2_bin_max=40000;
   x2=s2_bin_min:bin_s2:s2_bin_max;
   x2=x2';

%set up matricies to be filled
mean_S2_both_xy = zeros(floor(s2_xbins),floor(s2_ybins));
Sigma_S2_both = zeros(floor(s2_xbins),floor(s2_ybins));
Count_S2_both = zeros(floor(s2_xbins),floor(s2_ybins));

time_energy_cut = inrange(drift_time,[zcut_min,zcut_max]) & inrange(s1_phe_both_z,[100,s1area_bound_max]) & inrange(s2_phe_both_z,[s2_bin_min,s2_bin_max]) & inrange(s2radius,[rcut_min,rcut_max]) ...
    & inrange(s2_width,[100,1000]) & inrange(s1_width,[30,400]);

%%%%variables for uber turbo boost!
s2_phe_both_z_2=s2_phe_both_z(time_energy_cut);
s2x_2=s2x(time_energy_cut);
s2y_2=s2y(time_energy_cut);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   for x_bin=S2_xbin_min:S2xybinsize:(S2_xbin_max-S2xybinsize);
      %fprintf('X = %d \n',x_bin)
    for y_bin=S2_ybin_min:S2xybinsize:(S2_ybin_max-S2xybinsize);          
      x_min = x_bin; x_max = x_bin+S2xybinsize;
      y_min = y_bin; y_max = y_bin+S2xybinsize;      
      
        x_count=int32(1+x_bin/S2xybinsize+S2_xbin_max/S2xybinsize);
        y_count=int32(1+y_bin/S2xybinsize+S2_ybin_max/S2xybinsize); 
      
      
        bin_cut = inrange(s2x_2,[x_min,x_max]) & inrange(s2y_2,[y_min,y_max]);            
      
       hist_bin = hist(s2_phe_both_z_2(bin_cut) ,x2)'/bin_s2; 
       Count_S2_both(y_count,x_count)=sum(hist_bin)*bin_s2; % X and Y flip ... because MATLAB
       
       %%%%%%%%Truncate variables after binning is complete
        s2_phe_both_z_2=s2_phe_both_z_2(~bin_cut);
        s2x_2=s2x_2(~bin_cut);
        s2y_2=s2y_2(~bin_cut);
        clear bin_cut;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
      if Count_S2_both(y_count,x_count)>30;                         
                
                   amp_start=max(hist_bin);
                   mean_start=sum(hist_bin.*x2)/sum(hist_bin);
                   sigma_start=std(hist_bin.*x2);    
                Fit_S2_both_xy= fit(x2(2:end-1),hist_bin(2:end-1),'gauss1'); %remove edges
                
                %Peak location (mean)
                mean_S2_both_xy(y_count,x_count)=Fit_S2_both_xy.b1; % X and Y flip ... because MATLAB
                
                %1-sigma of peak position
                Sig_b=Fit_S2_both_xy.c1/sqrt(2)/sqrt(Count_S2_both(y_count,x_count)); %sigma/sqrt(N)
   

                    %Fit error checking
                        if(isnan(Sig_b))
                         mean_S2_both_xy(y_count,x_count)=0;
                         Sig_b=0;
                        end

                    %uncertainty in mean
                    Sigma_S2_both(y_count,x_count)=Sig_b;% X and Y flip ... because MATLAB
       
      else
          
        mean_S2_both_xy(y_count,x_count)=0;  %don't so gaussin fit if less than 30 events
       
      end             
    end     
   end
   
   mean_S2_both_xy(isnan(mean_S2_both_xy)) = 0;
       
   
%%Plot the mean S2_XY both %%%%%%%%%%%%%%%%

s2xbins=(S2_xbin_min+S2xybinsize/2):S2xybinsize:(S2_xbin_max-S2xybinsize/2);
s2ybins=(S2_ybin_min+S2xybinsize/2):S2xybinsize:(S2_ybin_max-S2xybinsize/2);


%Calculate corrections matrix and 1 sigma corrections matrix
% if inpaint_on==1;  mean_S2_both_xy(mean_S2_both_xy==0)=nan; mean_S2_both_xy=inpaint_nans(mean_S2_both_xy,3); end;
center=interp2(s2xbins,s2ybins,mean_S2_both_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
norm_S2_all=center./mean_S2_both_xy; 
norm_S2_all(isinf(norm_S2_all))=1;%no infinity! Do not correct events outside r=25cm
norm_S2_all(isnan(norm_S2_all))=1;

% if inpaint_on==1;  Sigma_S2_both(Sigma_S2_both==0)=nan; Sigma_S2_both=inpaint_nans(Sigma_S2_both,3); end;
sigma_center=interp2(s2xbins,s2ybins,Sigma_S2_both,x_center,y_center,'spline');%Normalize to the center (x=y=0)
sigma_norm_S2_all=sqrt((sigma_center./mean_S2_both_xy).^2+(Sigma_S2_both.*center./mean_S2_both_xy.^2).^2);
sigma_norm_S2_all(isinf(sigma_norm_S2_all))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
sigma_norm_S2_all(isnan(sigma_norm_S2_all))=1;%no nan. Set Sigma=1, which is 100% uncertainty


s2_both_center_xyz=center;
s2_both_center_xyz_error=sigma_center;



        
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Mapping S1 in XY (after correction for S1 Z-dependence). Both PMT
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear hist_bin x

   
    %set up energy bins for the gaussian fit
    bin_s1=8;%5
    s1_bin_min=100;
    s1_bin_max=500;
    x=s1_bin_min:bin_s1:s1_bin_max; %150:bin_s1:500;
    x=x';

%set up matricies to be filled
mean_S1_both_xy = zeros(floor(s1_xbins),floor(s1_ybins));
Sigma_S1_both = zeros(floor(s1_xbins),floor(s1_ybins));
Count_S1_both = zeros(floor(s1_xbins),floor(s1_ybins));

% s1_z_correction=polyval(P_s1_both,241.6)./polyval(P_s1_both,drift_distance_\mus);%normalize s1 to 160 ns.
s1_z_correction=polyval(P_s1_both,z_center)./polyval(P_s1_both,drift_time);%normalize s1 to 160 ns.

s1_phe_both_z=s1_phe_both.*s1_z_correction;
% s2_phe_both_z=s2_phe_both.*exp(drift_distance_\mus/e_lifetime); %Correct S2 for lifeimte  %s2_phe_both_z[5000,18000]
% s2_phe_both_z=s2_phe_both.*exp(drift_time/e_lifetime_both); %Correct S2 for lifeimte  %s2_phe_both_z[5000,18000]



time_energy_cut = inrange(drift_time,[zcut_min,zcut_max]) & inrange(s1_phe_both_z,[s1_bin_min,s1_bin_max]) & inrange(s2_phe_both_z,[1000,s2area_bound_max]) & inrange(s2radius,[rcut_min,rcut_max]) ...
    &inrange(s2_width,[100 1000]) & inrange(s1_width,[30 400]);


%%%%variables for uber turbo boost!
s1_phe_both_z_2=s1_phe_both_z(time_energy_cut);
s2x_2=s2x(time_energy_cut);
s2y_2=s2y(time_energy_cut);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   for x_bin=S2_xbin_min:S2xybinsize:(S2_xbin_max-S2xybinsize);
      %fprintf('X = %d \n',x_bin)
    for y_bin=S2_ybin_min:S2xybinsize:(S2_ybin_max-S2xybinsize);          
      x_min = x_bin; x_max = x_bin+S2xybinsize;
      y_min = y_bin; y_max = y_bin+S2xybinsize;      
      
        x_count=int32(1+x_bin/S2xybinsize+S2_xbin_max/S2xybinsize);
        y_count=int32(1+y_bin/S2xybinsize+S2_ybin_max/S2xybinsize); 
      
      
      bin_cut = inrange(s2x_2,[x_min,x_max]) & inrange(s2y_2,[y_min,y_max]);            
      
       hist_bin = hist(s1_phe_both_z_2(bin_cut) ,x)'/bin_s1; 
       Count_S1_both(y_count,x_count)=sum(hist_bin)*bin_s1; % X and Y flip ... because MATLAB
       
       %%%%%%%%Truncate variables after binning is complete
        s1_phe_both_z_2=s1_phe_both_z_2(~bin_cut);
        s2x_2=s2x_2(~bin_cut);
        s2y_2=s2y_2(~bin_cut);
        clear bin_cut;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
      if Count_S1_both(y_count,x_count)>30; 
                        
                    amp_start=max(hist_bin);
                    mean_start=sum(hist_bin.*x)/sum(hist_bin);
                    sigma_start=std(hist_bin.*x);                     
                Fit_S1_both_xy= fit(x(2:end-1),hist_bin(2:end-1),'gauss1'); %remove edges
                
                %Peak location (mean)
                mean_S1_both_xy(y_count,x_count)=Fit_S1_both_xy.b1; % X and Y flip ... because MATLAB
                
                %1-sigma of peak position
                Sig_b=Fit_S1_both_xy.c1/sqrt(2)/sqrt(Count_S1_both(y_count,x_count)); %sigma/sqrt(N)
   
                    %Fit error checking
                        if(isnan(Sig_b))
                          mean_S1_both_xy(y_count,x_count)=0;
                         Sig_b=0;
                        end
                                
                    %uncertainty in mean
                    Sigma_S1_both(y_count,x_count)=Sig_b;% X and Y flip ... because MATLAB   
     
      else
          
        mean_S1_both_xy(y_count,x_count)=0;  
       
      end            
    end    
   end
   
   mean_S1_both_xy(isnan(mean_S1_both_xy)) = 0;
   
   
%%Plot the correction S1_XY both %%%%%%%%%%%%

s1xbins=(S1_xbin_min+S1xybinsize/2):S1xybinsize:(S1_xbin_max-S1xybinsize/2);
s1ybins=(S1_ybin_min+S1xybinsize/2):S1xybinsize:(S1_ybin_max-S1xybinsize/2);

%Calculate corrections matrix and 1 sigma corrections matrix
% if inpaint_on==1;  mean_S1_both_xy(mean_S1_both_xy==0)=nan; mean_S1_both_xy=inpaint_nans(mean_S1_both_xy,3); end;
center=interp2(s1xbins,s1ybins,mean_S1_both_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
norm_S1_all=center./mean_S1_both_xy; 
norm_S1_all(isinf(norm_S1_all))=1;%no infinity! Do not correct events outside r=25cm
norm_S1_all(isnan(norm_S1_all))=1;

% if inpaint_on==1;  Sigma_S1_both(Sigma_S1_both==0)=nan; Sigma_S1_both=inpaint_nans(Sigma_S1_both,3); end;
sigma_center=interp2(s1xbins,s1ybins,Sigma_S1_both,x_center,y_center,'spline');%Normalize to the center (x=y=0)
sigma_norm_S1_all=sqrt((sigma_center./mean_S1_both_xy).^2+(Sigma_S1_both.*center./mean_S1_both_xy.^2).^2);
sigma_norm_S1_all(isinf(sigma_norm_S1_all))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
sigma_norm_S1_all(isnan(sigma_norm_S1_all))=1;%no nan. Set Sigma=1, which is 100% uncertainty

s1_both_center_xyz=center;
s1_both_center_xyz_error=sigma_center;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Mapping S1 in XYZ. 3D Map, caluculate and Submit IQ
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% if Kr_events>150000
%    
% %     LUXkrypCal_S1_3D(s1_phe_both,s1_phe_bottom,s2_phe_both,s2x,s2y,drift_distance_\mus,user_name,file_id_cp,dp_version,submit);
%     LUXkrypCal_S1_3D_Run04_local_2p14FixedZHisto(s1_phe_both,s1_phe_bottom,s2_phe_both,s2x,s2y,drift_time,x_center,y_center,z_center,det_edge,inpaint_on,user_name,file_id_cp,dp_version, algorithm_name, submit);
% 
% 
% end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Apply Corrections to 2015 DD Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% NEED TO LOAD THE DD DATA .MAT FILE HERE
clear s1_phe_both s2_phe_both
load('C:\Program Files\MATLAB\R2012a\bin\Run04_ERBand\Kr2p14_FixedZHisto\Oct2015_DD.mat');

% Apply the corrections
%-----------------------------------
% Finding the electron lifetime
%-----------------------------------

%  electron_lifetime = e_lifetime_both;

%-----------------------------------
% Finding the S1 Z-dep values
%-----------------------------------

z_dep_both_values = P_s1_both;
z_dep_par_all = z_dep_both_values;      

%------------------------------------------
% Finding the S2 xy correction map values
%------------------------------------------

   s2_x_bins = s2xbins;  
   s2_y_bins = s2ybins; 
   s2_map_all = norm_S2_all; 
   
%------------------------------------------
% Finding the S1 xy correction map values
%------------------------------------------

      
   s1_x_bins = s1xbins; 
   s1_y_bins = s1ybins; 
   s1_map_all = norm_S1_all; 

   
% %------------------------------------------
% % Finding the S1 xyz correction map values
% %------------------------------------------
   
% 
% if nargin == 3      
%    s1_xyz_x_bins = s1xbins_3D;
%    s1_xyz_y_bins = s1ybins_3D;
%    s1_xyz_z_bins = s1zbins_3D;
%    
% %    xx = size(s1_xyz_x_bins); yy = size(s1_xyz_y_bins); zz = size(s1_xyz_z_bins);
%    
%    s1_xyz_map_all = norm_s1_both_xyz; %lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.norm_s1_both_xyz; s1_xyz_map_all = reshape(s1_xyz_map_all,xx(2),yy(2),zz(2));
%    s1_xyz_map_bottom = norm_s1_bottom_xyz; %lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.norm_s1_bottom_xyz;   s1_xyz_map_bottom = reshape(s1_xyz_map_bottom,xx(2),yy(2),zz(2));
% end



%% Calculating corrections

%% Finding the drift time, x and y associated with S1 - this isn't how they are actually applied in DP, since DP estimates the S1 and here we know the exact S1.  
% It's the best we can do for ease here and it shouldn't matter much


s1_drift_ns = drift_time.*1000;
s1_ref_z_ns=z_center*1000;
s1_x_cm = s2x;
s1_y_cm = s2y;

%--------------------------------------------------------------------------
% Calculate Z corrections
%--------------------------------------------------------------------------

% Calculate S1 Z-correction
s1_z_correction = polyval(z_dep_par_all,s1_ref_z_ns./1000)./polyval(z_dep_par_all,s1_drift_ns./1000);
s1_z_correction(find(isnan(s1_z_correction))) = 1.0;

s2_z_correction = RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z, drift_time);

% % Calculate electron lifetime correction (S2 Z-correction)
% electron_lifetime_correction = exp(drift_time./electron_lifetime);
% electron_lifetime_correction(find(isnan(electron_lifetime_correction))) = 1.0;
 
%--------------------------------------------------------------------------
% Calculate XY corrections
%--------------------------------------------------------------------------

% Calculate S1 XY-only corrections
s1xy_correction = interp2(s1_x_bins,s1_y_bins,s1_map_all,s1_x_cm,s1_y_cm,'spline',1);
s1xy_correction(find(s1xy_correction==0))=1.0;  s1xy_correction(find(isnan(s1xy_correction)))=1.0;

% Calculate S2 XY corrections
s2xy_correction = interp2(s2_x_bins,s2_y_bins,s2_map_all,s2x,s2y,'spline',1);
s2xy_correction(find(s2xy_correction==0))=1.0;  s2xy_correction(find(isnan(s2xy_correction)))=1.0;

%--------------------------------------------------------------------------
% Calculate XYZ corrections
%--------------------------------------------------------------------------

s1xyz_correction=s1xy_correction.*s1_z_correction;

%% add RQs of correction factors

% d.correction_electron_lifetime = single(electron_lifetime_correction);
d.correction_s1_z_dependence = single(s1_z_correction);
d.correction_s1_xyz_dependence = single(s1xyz_correction); %this is either 3D or 1D+2D
d.correction_s1_xy_dependence = single(s1xy_correction);
d.correction_s2_xy_dependence = single(s2xy_correction);

% d.admin.corrections.electron_lifetime = electron_lifetime;
d.admin.corrections.s1_z_dependence.all = z_dep_par_all;

%--------------------------------------------------------------------------
% Apply Z corrections
%--------------------------------------------------------------------------

s1_phe_both_z=s1_phe_both.*s1_z_correction;
s2_phe_both_z=s2_phe_both.*s2_z_correction;
  

s1_phe_both_xyz=s1_phe_both.*s1xyz_correction;
s2_phe_both_xyz=s2_phe_both.*s2_z_correction.*s2xy_correction;
  %% Redo with r<10 cm

fid_cut=(drift_time< -s2radius.^2+280) & drift_time>50 & drift_time<250;

DD_s1_2015=s1_phe_both_xyz(fid_cut);
DD_s2_2015=s2_phe_both_xyz(fid_cut);

%% NEED TO COME UP WITH CUTS TO SELECT NON DD DATA... FOR INSTANCE, NR CUT?  Also a slopped radial cut to cut out wall
 
%% Make the plots


g1=0.10; ee=0.69;
%           g2_sep2014=27.21*ee; %converting to g2 with SE size
          g2=25.28*ee;
%Kr 2014 data

DD_E_2015=(1/73).*(DD_s1_2015./g1+DD_s2_2015./g2);

        
figure
step(DD_E_2015,[0:3:500],'k')
hold on;
line([163.9 163.9],[0 length(DD_E_2015)],'Color','r','LineWidth',2) %Xe131
% line([236.1 236.1],[0 length(DD_E_2015)],'Color','r','LineWidth',2) %Xe129
% % line([280 280],[0 length(DD_E_2015)],'Color','r','LineWidth',2) %Xe125?
% line([40 40],[0 length(DD_E_2015)],'Color','r','LineWidth',2) 
% line([80 80],[0 length(DD_E_2015)],'Color','r','LineWidth',2)
% % line([33.2 33.2],[0 length(DD_E_2015)],'Color','r','LineWidth',2) %Xray (127?)
% % line([5.3 5.3],[0 length(DD_E_2015)],'Color','r','LineWidth',2) %xray (127?)
% line([275 275],[0 length(DD_E_2015)],'Color','r','LineWidth',2)
ylim([0 100]);
xlim([20 400]);
xlabel('Energy (keV)'); ylabel('Counts'); myfigview(16); 
legend('Data','Expected');

ht=text(0,0,'Xe125 (275 keV)');
set(ht,'FontSize',18)
set(ht,'Rotation',90)
set(ht,'Color',[1 0 0])

ht=text(0,0,'Xe129 (236 keV)');
set(ht,'FontSize',18)
set(ht,'Rotation',90)
set(ht,'Color',[1 0 0])

ht=text(0,0,'Xe 131 (163.9 keV)');
set(ht,'FontSize',18)
set(ht,'Rotation',90)
set(ht,'Color',[1 0 0])

ht=text(0,0,'Xe 129 (40 keV)');
set(ht,'FontSize',18)
set(ht,'Rotation',90)
set(ht,'Color',[1 0 0])

ht=text(0,0,'Xe 131 (80 keV)');
set(ht,'FontSize',18)
set(ht,'Rotation',90)
set(ht,'Color',[1 0 0])

%% Fit to the mose distinguishable DD peaks

cut_40kev=inrange(DD_E_2015,[30 60]);
cut_80kev=inrange(DD_E_2015,[70 95]);
cut_164kev=inrange(DD_E_2015,[140 180]);
cut_236kev=inrange(DD_E_2015,[210,250]);
cut_275kev=inrange(DD_E_2015,[265,285]);

fit_40kev=fit([30:3:60].',hist(DD_E_2015(cut_40kev),[30:3:60]).','gauss1');
fit_80kev=fit([70:3:95].',hist(DD_E_2015(cut_80kev),[70:3:95]).','gauss1');
fit_164kev=fit([140:3:180].',hist(DD_E_2015(cut_164kev),[140:3:180]).','gauss1');
fit_236kev=fit([210:3:250].',hist(DD_E_2015(cut_236kev),[210:3:250]).','gauss1');
fit_275kev=fit([265:3:285].',hist(DD_E_2015(cut_275kev),[265:3:285]).','gauss1');

resolution_40keV=(fit_40kev.c1/sqrt(2))/fit_40kev.b1
resolution_80keV=(fit_80kev.c1/sqrt(2))/fit_80kev.b1
resolution_164keV=(fit_164kev.c1/sqrt(2))/fit_164kev.b1
resolution_236keV=(fit_236kev.c1/sqrt(2))/fit_236kev.b1
resolution_275keV=(fit_275kev.c1/sqrt(2))/fit_275kev.b1

%show the fits
% figure
% step(DD_E_2015,[0:3:500],'-k');
hold on;
x=[25:1:60]; y=fit_40kev.a1.*exp(-((x-fit_40kev.b1)./fit_40kev.c1).^2);
plot(x,y,'-b','LineWidth',2);
x=[60:3:120]; y=fit_80kev.a1.*exp(-((x-fit_80kev.b1)./fit_80kev.c1).^2);
plot(x,y,'-b','LineWidth',2);
x=[130:3:190]; y=fit_164kev.a1.*exp(-((x-fit_164kev.b1)./fit_164kev.c1).^2);
plot(x,y,'-b','LineWidth',2);
x=[200:3:260]; y=fit_236kev.a1.*exp(-((x-fit_236kev.b1)./fit_236kev.c1).^2);
plot(x,y,'-b','LineWidth',2);
x=[260:3:295]; y=fit_275kev.a1.*exp(-((x-fit_275kev.b1)./fit_275kev.c1).^2);
plot(x,y,'-b','LineWidth',2);
line([236.1 236.1],[0 length(DD_E_2015)],'Color','r','LineWidth',2) %Xe129
line([40 40],[0 length(DD_E_2015)],'Color','r','LineWidth',2) 
line([80 80],[0 length(DD_E_2015)],'Color','r','LineWidth',2)
line([275 275],[0 length(DD_E_2015)],'Color','r','LineWidth',2)
ylim([0 100]);
xlim([20 400]);
legend('Data','Expected','Fits');

fit_40kev.b1
fit_40kev.c1/sqrt(2)
temp=confint(fit_40kev,0.68);
abs(fit_40kev.b1-temp(1,2))
abs(fit_40kev.c1-temp(1,3))/sqrt(2)

fit_80kev.b1
fit_80kev.c1/sqrt(2)
temp=confint(fit_80kev,0.68);
abs(fit_80kev.b1-temp(1,2))
abs(fit_80kev.c1-temp(1,3))/sqrt(2)

fit_164kev.b1
fit_164kev.c1/sqrt(2)
temp=confint(fit_164kev,0.68);
abs(fit_164kev.b1-temp(1,2))
abs(fit_164kev.c1-temp(1,3))/sqrt(2)


fit_236kev.b1
fit_236kev.c1/sqrt(2)
temp=confint(fit_236kev,0.68);
abs(fit_236kev.b1-temp(1,2))
abs(fit_236kev.c1-temp(1,3))/sqrt(2)

fit_275kev.b1
fit_275kev.c1/sqrt(2)
temp=confint(fit_275kev,0.68);
abs(fit_275kev.b1-temp(1,2))
abs(fit_275kev.c1-temp(1,3))/sqrt(2)