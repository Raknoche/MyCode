%% First Get the DD Band in Z Slices

load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p21\KrypCal2p21_QualityChecks\Discrimination\MappingParametersFit_Start_DD');
         %%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 % Removing field dependence from S1 and S2
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  %Hardcoded numbers from Energy Chi2 Fit
           
             
          
s1ab_s2_xyz_fit_map.p1=-0.4988; 
s1ab_s2_xyz_fit_map.p2=1.484;
c1=1-s1ab_s2_xyz_fit_map.p1.*2.7303.^2-2.7303.*s1ab_s2_xyz_fit_map.p2;
s1ab_s2_xyz_fit_map.p3=c1;

s1ab_s1_xyz_fit_map.p1=0.065;
s1ab_s1_xyz_fit_map.p2=0.02;
c2=1-s1ab_s1_xyz_fit_map.p1.*2.7303.^2-2.7303.*s1ab_s1_xyz_fit_map.p2;
s1ab_s1_xyz_fit_map.p3=c2;

 
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


 else %Only do 1D in this case
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

 
s2_phe_bottom=s2_phe_bottom./(s2_zfieldremoval);
s2_phe_both=s2_phe_both./(s2_zfieldremoval);
s1_phe_bottom=s1_phe_bottom./(s1_zfieldremoval);
s1_phe_both=s1_phe_both./(s1_zfieldremoval);



 end 


    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Calculating the S1 Z-dependence (Both PMT arrays)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    clear y_s1_both y_s1_bottom Fit_s1_both Fit_s1_bottom mean_s1_both mean_s1_both
    clear sigma_s1_both sigma_s1_both mean_index means mean_index means_error x temp_hist
    clear hist_s1_both x xfit S
    
    %set up energy bins for the gaussian fit
    bin_s1=8;%5
    s1_bin_min=100;
    s1_bin_max=500;
    x=s1_bin_min:bin_s1:s1_bin_max; %150:bin_s1:500;
    x=x';
        
    cut = inrange(s2radius,[rcut_min,rcut_max]) & inrange(s1_phe_both,[s1_bin_min,s1_bin_max]) & inrange(s2_width,[100,1000]) & inrange(s1_width,[30,400]);   
    
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
    
%     
%     light_correction_fig = figure;
%     errorbar(mean_index,means,means_error,'.k');
%     hold on     
%     [yplot_s1_both] = polyval(P_s1_both,dT_fit);
%     [center_s1_both, center_s1_both_err] = polyval(P_s1_both,z_center, S); %160[us]*1.51[\mus/us]
%     plot(dT_fit,yplot_s1_both,'-r','LineWidth',1);    
%     xlabel('Depth (\mus)', 'FontSize',16); ylabel('Mean (Pulse Area Phe)', 'FontSize', 16);
%     title(strcat('.Kr83m. S1 Mean vs. Z'), 'FontSize', 16,'Interpreter','none');
%     line([0 det_edge],[center_s1_both center_s1_both],'linestyle','--');
%     line([z_center; z_center;],[min(means) max(means)],'linestyle','--');
%     legend('S1\_both Mean', strcat('y=', num2str(P_s1_both(1),'%10.2e'), '*Z^2 + ', num2str(P_s1_both(2),'%10.2e')...
%         , '*Z + ', num2str(P_s1_both(3),4)),strcat('Det center = ',num2str(center_s1_both,4), ' [Phe]'),'location','northwest')
%     
%     box on;
%     myfigview(16);
%     xlim([0 det_edge]);
%     ylim([0.8*min(means) 1.15*max(means)]);

      
    s1_both_means=means;
    s1_both_means_sigma=means_error;
    s1_both_means_bin=mean_index;
    

    
    
       
    
%     s1_z_correction=polyval(P_s1_both,241.6)./polyval(P_s1_both,drift_distance_\mus);%normalize s1 to 160 ns.
    s1_z_correction=polyval(P_s1_both,z_center)./polyval(P_s1_both,drift_time);
    
    %s1_both_z_correction=interp1(s1_both_means_bin,s1_both_means,160)./interp1(s1_both_means_bin,s1_both_means,drift_time);
    s1_phe_both_z=s1_phe_both.*s1_z_correction;


%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   Calculating the S2 Z-dependence (i.e. electron lifetime)
%   Using Both PMT arrays
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
     
     
     %%calculate S2/S1 also, just to double check the lifetime from S2 only
     % hist_s2_both_s1(:,bin)= hist(s2_phe_both(time_cut)./s1_phe_both_z(time_cut),x3)';
     % temp_hist=hist_s2_both_s1(:,bin);   
      %Fit_s2_both_s1{bin}=fit(x3(2:end-1),temp_hist(2:end-1),'gauss1');%remove bin edges        
       % means_s2_s1(bin)=Fit_s2_both_s1{bin}.b1;
       % means_error_s2_s1(bin) = Fit_s2_both_s1{bin}.c1/2/sqrt(sum(hist_s2_both_s1(:,bin))); % == sigma/sqrt(N);
       % mean_index_s2_s1(bin) = (binStart+binEnd)/2;     

    end
    
 % //////////////////////////////Make Plots/////////////////////////////////
 
  
     dT_range= inrange( mean_index,[5+bin_size,500] );%from 50 to 300us for the fit


    s2_both_norm_z_index=mean_index(dT_range);
    s2_both_norm_z_means=means(dT_range);
    s2_at_top_both=RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,4); %normalize to right below the gate
    s2_both_norm_z=RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,4)./s2_both_norm_z_means;
  
    %%Tracker for pseudo-lifetime
    s2_at_bottom_both=RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,320);
    pseudo_lifetime_both=-320/log(s2_at_bottom_both/s2_at_top_both);
 
    dT_fit=0:1:det_edge;    
    
    Fit_EL= fit(mean_index(dT_range),means(dT_range),'exp1');
    yfitEL = Fit_EL.a.*exp(Fit_EL.b.*dT_fit);
    
    e_lifetime_both=-1/Fit_EL.b;
    s2_z0_both=Fit_EL.a;
       
    b=confint(Fit_EL, 0.683);%1 sigma matrix
    sigma_e_lifetime_both=(1/b(1,2)-1/b(2,2))/2;
    sigma_s2_z0_both=(b(2,1)-b(1,1))/2;
%  
%     lifetime_fig_both = figure;
%     errorbar(mean_index(2:end),means(2:end),means_error(2:end),'.k');  %first bin is anode, skip it
%     hold on     
%     plot(dT_fit,yfitEL,'-r','LineWidth',1);   
%     xlabel('Time (\mus)', 'FontSize',16); ylabel('Mean S2 (Pulse Area Phe)', 'FontSize', 16);
%     title(strcat('.Kr83m. S2_both Electron Lifetime'), 'FontSize', 16,'Interpreter','none');
%     legend('S2\_both Mean', strcat( '\lambda= ', num2str(e_lifetime_both,4), ' \pm ', num2str(sigma_e_lifetime_both,2), ' \mus' )...
%         ,'location','northeast');
%     myfigview(16);
%     xlim([0 det_edge]);   
% 
%     
        
    s2_both_means=means;
    s2_both_means_sigma=means_error;
    s2_both_means_bin=mean_index;    
    
    s2_phe_both_z=s2_phe_both.*RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z, drift_time);

    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Mapping S2 in XY (after correcting for Z-dep). Both PMT
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    clear x2 hist_bin

    %set up energy bins for the gaussian fit. Using both S2 only
   bin_s2=500;
   s2_bin_min=1000;
   s2_bin_max=40000;
   x2=s2_bin_min:bin_s2:s2_bin_max;
   x2=x2';

%set up matricies to be filled
mean_S2_both_xy = zeros(floor(s2_ybins),floor(s2_xbins));
Sigma_S2_both = zeros(floor(s2_ybins),floor(s2_xbins));
Count_S2_both = zeros(floor(s2_ybins),floor(s2_xbins));

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
       
      if Count_S2_both(y_count,x_count)>350;                         
                
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

color_range_max=max(max(mean_S2_both_xy));
color_range_min=min(min(mean_S2_both_xy(mean_S2_both_xy>0)));
vc_step=(color_range_max-color_range_min)/50;

vc=color_range_min:vc_step:color_range_max;
%  
% s2_xy_both_fig = figure;
% contourf(s2xbins,s2ybins,mean_S2_both_xy,vc,'LineColor','none');
% xlabel('x (cm)','fontsize', 18);
% ylabel('y (cm)','fontsize', 18);
% title(strcat('. S2 Mean (Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
% caxis([color_range_min color_range_max])
% g=colorbar;
% ylabel(g,'phe','FontSize',16)
% set(g,'FontSize',18)
% hold on
% LUXPlotTopPMTs
% axis([-25 25 -25 25])
% myfigview(16);


%Calculate corrections matrix and 1 sigma corrections matrix
% if inpaint_on==1;  mean_S2_both_xy(mean_S2_both_xy==0)=nan; mean_S2_both_xy=inpaint_nans(mean_S2_both_xy,3); end;
center=interp2(s2xbins,s2ybins,mean_S2_both_xy,x_center,y_center,'cubic');%Normalize to the center (x=y=0)
norm_S2_all=center./mean_S2_both_xy; 
norm_S2_all(isinf(norm_S2_all))=1;%no infinity! Do not correct events outside r=25cm
norm_S2_all(isnan(norm_S2_all))=1;

% if inpaint_on==1;  Sigma_S2_both(Sigma_S2_both==0)=nan; Sigma_S2_both=inpaint_nans(Sigma_S2_both,3); end;
sigma_center=interp2(s2xbins,s2ybins,Sigma_S2_both,x_center,y_center,'cubic');%Normalize to the center (x=y=0)
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
mean_S1_both_xy = zeros(floor(s1_ybins),floor(s1_xbins));
Sigma_S1_both = zeros(floor(s1_ybins),floor(s1_xbins));
Count_S1_both = zeros(floor(s1_ybins),floor(s1_xbins));

% s1_z_correction=polyval(P_s1_both,241.6)./polyval(P_s1_both,drift_distance_\mus);%normalize s1 to 160 ns.
s1_z_correction=polyval(P_s1_both,z_center)./polyval(P_s1_both,drift_time);%normalize s1 to 160 ns.

s1_phe_both_z=s1_phe_both.*s1_z_correction;
% s2_phe_both_z=s2_phe_both.*exp(drift_distance_\mus/e_lifetime); %Correct S2 for lifeimte  %s2_phe_both_z[5000,18000]
% s2_phe_both_z=s2_phe_both.*exp(drift_time/e_lifetime_both); %Correct S2 for lifeimte  %s2_phe_both_z[5000,18000]
s2_phe_both_z=s2_phe_both.*RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z, drift_time);



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
       
      if Count_S1_both(y_count,x_count)>350; 
                        
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

color_range_max=max(max(mean_S1_both_xy));
color_range_min=min(min(mean_S1_both_xy(mean_S1_both_xy>0)));
vc_step=(color_range_max-color_range_min)/50;

vc=color_range_min:vc_step:color_range_max;
 
%  
% s1_xy_both_fig = figure;
% contourf(s1xbins,s1ybins,mean_S1_both_xy,vc,'LineColor','none');
% xlabel('x (cm)','fontsize', 18);
% ylabel('y (cm)','fontsize', 18);
% title(strcat('. S1 Mean (Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
% caxis([color_range_min color_range_max])
% g=colorbar;
% ylabel(g,'phe','FontSize',16)
% set(g,'FontSize',18)
% hold on
% LUXPlotTopPMTs
% axis([-25 25 -25 25])
% myfigview(16);


%Calculate corrections matrix and 1 sigma corrections matrix
% if inpaint_on==1;  mean_S1_both_xy(mean_S1_both_xy==0)=nan; mean_S1_both_xy=inpaint_nans(mean_S1_both_xy,3); end;
center=interp2(s1xbins,s1ybins,mean_S1_both_xy,x_center,y_center,'cubic');%Normalize to the center (x=y=0)
norm_S1_all=center./mean_S1_both_xy; 
norm_S1_all(isinf(norm_S1_all))=1;%no infinity! Do not correct events outside r=25cm
norm_S1_all(isnan(norm_S1_all))=1;

% if inpaint_on==1;  Sigma_S1_both(Sigma_S1_both==0)=nan; Sigma_S1_both=inpaint_nans(Sigma_S1_both,3); end;
sigma_center=interp2(s1xbins,s1ybins,Sigma_S1_both,x_center,y_center,'cubic');%Normalize to the center (x=y=0)
sigma_norm_S1_all=sqrt((sigma_center./mean_S1_both_xy).^2+(Sigma_S1_both.*center./mean_S1_both_xy.^2).^2);
sigma_norm_S1_all(isinf(sigma_norm_S1_all))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
sigma_norm_S1_all(isnan(sigma_norm_S1_all))=1;%no nan. Set Sigma=1, which is 100% uncertainty

s1_both_center_xyz=center;
s1_both_center_xyz_error=sigma_center;

%% Apply corrections to 7p5 DD data and calculate the NRband

load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p21\KrypCal2p21_QualityChecks\NRBand_Z_Dependence\MutliZ_DD_7p5.mat');
clear s1_phe_both_xyz s2_phe_both_xyz s1_phe_both_z s2_phe_both_z

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


%------------------------------------------
% Calculating corrections
%------------------------------------------

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
 
%--------------------------------------------------------------------------
% Calculate XY corrections
%--------------------------------------------------------------------------

% Calculate S1 XY-only corrections
s1xy_correction = interp2(s1_x_bins,s1_y_bins,s1_map_all,s1_x_cm,s1_y_cm,'cubic');
s1xy_correction(find(s1xy_correction==0))=1.0;  s1xy_correction(find(isnan(s1xy_correction)))=1.0;

% Calculate S2 XY corrections
s2xy_correction = interp2(s2_x_bins,s2_y_bins,s2_map_all,s2x,s2y,'cubic');
s2xy_correction(find(s2xy_correction==0))=1.0;  s2xy_correction(find(isnan(s2xy_correction)))=1.0;

%--------------------------------------------------------------------------
% Calculate XYZ corrections
%--------------------------------------------------------------------------

s1xyz_correction=s1xy_correction.*s1_z_correction;

%--------------------------------------------------------------------------
% Apply corrections
%--------------------------------------------------------------------------

s1_phe_both_z=s1_phe_both.*s1_z_correction;
s2_phe_both_z=s2_phe_both.*s2_z_correction;

% Calculate Energy spectra
s1_phe_both_xyz=s1_phe_both.*s1xyz_correction;
s2_phe_both_xyz=s2_phe_both.*s2_z_correction.*s2xy_correction;


%% Make NR Band
max_r=18;
s1_min=0;
s1_max=50;


pulse_area_cut=s2_phe_both_xyz>150 & s2_phe_both<10000 & s1_phe_both<300;
NR_beam_cut=inrange(drift_time,[35 75]) & inrange(s2radius_del,[0,max_r]) & ...
    (s2y < -8.*s2x+70) & (s2y>-20.*s2x+10);


bin_counter=1;
  
DD_cut= clean_cut & pulse_area_cut  & NR_beam_cut;

DD_s1_all=s1_phe_both_xyz(DD_cut);
DD_s2_all=s2_phe_both_xyz(DD_cut);

%%
load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p21\KrypCal2p21_QualityChecks\NRBand_Z_Dependence\MutliZ_DD_15.mat');
clear s1_phe_both_xyz s2_phe_both_xyz s1_phe_both_z s2_phe_both_z

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


%------------------------------------------
% Calculating corrections
%------------------------------------------

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
 
%--------------------------------------------------------------------------
% Calculate XY corrections
%--------------------------------------------------------------------------

% Calculate S1 XY-only corrections
s1xy_correction = interp2(s1_x_bins,s1_y_bins,s1_map_all,s1_x_cm,s1_y_cm,'cubic');
s1xy_correction(find(s1xy_correction==0))=1.0;  s1xy_correction(find(isnan(s1xy_correction)))=1.0;

% Calculate S2 XY corrections
s2xy_correction = interp2(s2_x_bins,s2_y_bins,s2_map_all,s2x,s2y,'cubic');
s2xy_correction(find(s2xy_correction==0))=1.0;  s2xy_correction(find(isnan(s2xy_correction)))=1.0;

%--------------------------------------------------------------------------
% Calculate XYZ corrections
%--------------------------------------------------------------------------

s1xyz_correction=s1xy_correction.*s1_z_correction;

%--------------------------------------------------------------------------
% Apply corrections
%--------------------------------------------------------------------------

s1_phe_both_z=s1_phe_both.*s1_z_correction;
s2_phe_both_z=s2_phe_both.*s2_z_correction;

% Calculate Energy spectra
s1_phe_both_xyz=s1_phe_both.*s1xyz_correction;
s2_phe_both_xyz=s2_phe_both.*s2_z_correction.*s2xy_correction;


%% Make NR Band
max_r=18;
s1_min=0;
s1_max=50;


pulse_area_cut=s2_phe_both_xyz>150 & s2_phe_both<10000 & s1_phe_both<300;
NR_beam_cut=inrange(drift_time,[75 125]) & inrange(s2radius_del,[0,max_r]) & ...
    (s2y < -8.*s2x+70) & (s2y>-20.*s2x+10);



  
DD_cut= clean_cut & pulse_area_cut & NR_beam_cut;

DD_s1_all=[DD_s1_all; s1_phe_both_xyz(DD_cut)];
DD_s2_all=[DD_s2_all; s2_phe_both_xyz(DD_cut)];


%%
load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p21\KrypCal2p21_QualityChecks\NRBand_Z_Dependence\MutliZ_DD_22p5.mat');
clear s1_phe_both_xyz s2_phe_both_xyz s1_phe_both_z s2_phe_both_z

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


%------------------------------------------
% Calculating corrections
%------------------------------------------

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
 
%--------------------------------------------------------------------------
% Calculate XY corrections
%--------------------------------------------------------------------------

% Calculate S1 XY-only corrections
s1xy_correction = interp2(s1_x_bins,s1_y_bins,s1_map_all,s1_x_cm,s1_y_cm,'cubic');
s1xy_correction(find(s1xy_correction==0))=1.0;  s1xy_correction(find(isnan(s1xy_correction)))=1.0;

% Calculate S2 XY corrections
s2xy_correction = interp2(s2_x_bins,s2_y_bins,s2_map_all,s2x,s2y,'cubic');
s2xy_correction(find(s2xy_correction==0))=1.0;  s2xy_correction(find(isnan(s2xy_correction)))=1.0;

%--------------------------------------------------------------------------
% Calculate XYZ corrections
%--------------------------------------------------------------------------

s1xyz_correction=s1xy_correction.*s1_z_correction;

%--------------------------------------------------------------------------
% Apply corrections
%--------------------------------------------------------------------------

s1_phe_both_z=s1_phe_both.*s1_z_correction;
s2_phe_both_z=s2_phe_both.*s2_z_correction;

% Calculate Energy spectra
s1_phe_both_xyz=s1_phe_both.*s1xyz_correction;
s2_phe_both_xyz=s2_phe_both.*s2_z_correction.*s2xy_correction;


%% Make NR Band
max_r=21;
s1_min=0;
s1_max=50;


pulse_area_cut=s2_phe_both_xyz>150 & s2_phe_both<10000 & s1_phe_both<300;
NR_beam_cut=inrange(drift_time,[125 175]) & inrange(s2radius_del,[0,max_r]) & ...
    (s2y < -8.*s2x+50) & (s2y>-20.*s2x+10);



  
DD_cut= clean_cut & pulse_area_cut & NR_beam_cut;

DD_s1_all=[DD_s1_all; s1_phe_both_xyz(DD_cut)];
DD_s2_all=[DD_s2_all; s2_phe_both_xyz(DD_cut)];

%%
load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p21\KrypCal2p21_QualityChecks\NRBand_Z_Dependence\MutliZ_DD_30.mat');
clear s1_phe_both_xyz s2_phe_both_xyz s1_phe_both_z s2_phe_both_z

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


%------------------------------------------
% Calculating corrections
%------------------------------------------

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
 
%--------------------------------------------------------------------------
% Calculate XY corrections
%--------------------------------------------------------------------------

% Calculate S1 XY-only corrections
s1xy_correction = interp2(s1_x_bins,s1_y_bins,s1_map_all,s1_x_cm,s1_y_cm,'cubic');
s1xy_correction(find(s1xy_correction==0))=1.0;  s1xy_correction(find(isnan(s1xy_correction)))=1.0;

% Calculate S2 XY corrections
s2xy_correction = interp2(s2_x_bins,s2_y_bins,s2_map_all,s2x,s2y,'cubic');
s2xy_correction(find(s2xy_correction==0))=1.0;  s2xy_correction(find(isnan(s2xy_correction)))=1.0;

%--------------------------------------------------------------------------
% Calculate XYZ corrections
%--------------------------------------------------------------------------

s1xyz_correction=s1xy_correction.*s1_z_correction;

%--------------------------------------------------------------------------
% Apply corrections
%--------------------------------------------------------------------------

s1_phe_both_z=s1_phe_both.*s1_z_correction;
s2_phe_both_z=s2_phe_both.*s2_z_correction;

% Calculate Energy spectra
s1_phe_both_xyz=s1_phe_both.*s1xyz_correction;
s2_phe_both_xyz=s2_phe_both.*s2_z_correction.*s2xy_correction;


%% Make NR Band
max_r=18;
s1_min=0;
s1_max=50;


pulse_area_cut=s2_phe_both_xyz>150 & s2_phe_both<10000 & s1_phe_both<300;
NR_beam_cut=inrange(drift_time,[180 225]) & inrange(s2radius_del,[0,max_r]) & ...
    (s2y < -8.*s2x+36) & (s2y>-20.*s2x+5);



  
DD_cut= clean_cut & pulse_area_cut & NR_beam_cut;

DD_s1_all=[DD_s1_all; s1_phe_both_xyz(DD_cut)];
DD_s2_all=[DD_s2_all; s2_phe_both_xyz(DD_cut)];

%%
load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p21\KrypCal2p21_QualityChecks\NRBand_Z_Dependence\MutliZ_DD_37p5.mat');
clear s1_phe_both_xyz s2_phe_both_xyz s1_phe_both_z s2_phe_both_z

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


%------------------------------------------
% Calculating corrections
%------------------------------------------

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
 
%--------------------------------------------------------------------------
% Calculate XY corrections
%--------------------------------------------------------------------------

% Calculate S1 XY-only corrections
s1xy_correction = interp2(s1_x_bins,s1_y_bins,s1_map_all,s1_x_cm,s1_y_cm,'cubic');
s1xy_correction(find(s1xy_correction==0))=1.0;  s1xy_correction(find(isnan(s1xy_correction)))=1.0;

% Calculate S2 XY corrections
s2xy_correction = interp2(s2_x_bins,s2_y_bins,s2_map_all,s2x,s2y,'cubic');
s2xy_correction(find(s2xy_correction==0))=1.0;  s2xy_correction(find(isnan(s2xy_correction)))=1.0;

%--------------------------------------------------------------------------
% Calculate XYZ corrections
%--------------------------------------------------------------------------

s1xyz_correction=s1xy_correction.*s1_z_correction;

%--------------------------------------------------------------------------
% Apply corrections
%--------------------------------------------------------------------------

s1_phe_both_z=s1_phe_both.*s1_z_correction;
s2_phe_both_z=s2_phe_both.*s2_z_correction;

% Calculate Energy spectra
s1_phe_both_xyz=s1_phe_both.*s1xyz_correction;
s2_phe_both_xyz=s2_phe_both.*s2_z_correction.*s2xy_correction;


%% Make NR Band
max_r=18;
s1_min=0;
s1_max=50;


pulse_area_cut=s2_phe_both_xyz>150 & s2_phe_both<10000 & s1_phe_both<300;
NR_beam_cut=inrange(drift_time,[235 280]) & inrange(s2radius_del,[0,max_r]) & ...
    (s2y < -8.*s2x+36) & (s2y>-15.*s2x-5);



  
DD_cut= clean_cut & pulse_area_cut & NR_beam_cut;

DD_s1_all=[DD_s1_all; s1_phe_both_xyz(DD_cut)];
DD_s2_all=[DD_s2_all; s2_phe_both_xyz(DD_cut)];

%% Make Band Plot
   band_sigma=1.28; %90% confidence bounds
    
             
            s2_s1_bin= 1:0.05:5;
           NR_bin_size=3; %3 works
           index_i=1;
           NR_bins=NR_bin_size/2+s1_min:NR_bin_size:s1_max;
               lower_bound=ones(size(NR_bins));
               mean_NR_band=ones(size(NR_bins));
               upper_bound=ones(size(NR_bins));
                counts_NR_hist=ones(size(NR_bins));
               NR_sigma=ones(size(NR_bins));
               NR_sigma_sigma=ones(size(NR_bins));
            clear mean_NR_band_err lower_bound_err upper_bound_err lower_bound mean_NR_band upper_bound
        
%             temporary, for calculating correct E-lifetime   
        for NR_bin = (NR_bin_size/2)+s1_min:NR_bin_size:(s1_max); %HIST uses bin centers
   
        
            NR_fit_cut=inrange(DD_s1_all,[NR_bin-NR_bin_size/2,NR_bin+NR_bin_size/2]);
             
            NR_band_hist=hist(log10(DD_s2_all(NR_fit_cut)./DD_s1_all(NR_fit_cut)),s2_s1_bin);
            counts_NR_hist(index_i)=sum(NR_band_hist);
           
            Fit_NR=fit(s2_s1_bin',NR_band_hist','gauss1');
             
            rms_fit=fit_ML_normal(log10(DD_s2_all(NR_fit_cut)./DD_s1_all(NR_fit_cut)),s2_s1_bin);
            
            %band=confint(Fit_NR,.1);
%             NR_sigma(index_i)= Fit_NR.c1/sqrt(2);
            NR_sigma(index_i)= sqrt(rms_fit.sig2);
            NR_sigma_sigma(index_i)=NR_sigma(index_i)/sqrt(2*counts_NR_hist(index_i));
            
            lower_bound(index_i)=Fit_NR.b1-NR_sigma(index_i)*band_sigma;
            mean_NR_band(index_i)=Fit_NR.b1;
            upper_bound(index_i)=Fit_NR.b1+NR_sigma(index_i)*band_sigma;

            temp_err=confint(Fit_NR,0.68);            
            mean_NR_band_err(index_i)=abs(temp_err(1,2)-Fit_NR.b1);
            lower_bound_err(index_i)=sqrt(mean_NR_band_err(index_i).^2 + (band_sigma.*abs(temp_err(1,3)-Fit_NR.c1)/sqrt(2)).^2);
            upper_bound_err(index_i)=sqrt(mean_NR_band_err(index_i).^2 + (band_sigma.*abs(temp_err(1,3)-Fit_NR.c1)/sqrt(2)).^2);

            
            %mean_NR_band_sim(index_i)=Fit_NR_sim.b1;
            
            index_i=index_i+1;
         end
        
         
         % Fit power law to the mean
         g=fittype('a*x^b');
         NR_mean_power_fit=fit(NR_bins(NR_bins>1)',mean_NR_band(NR_bins>1)',g,'startpoint',[ 2.5 -.1]);
         NR_mean_fit=NR_mean_power_fit.a.*NR_bins.^(NR_mean_power_fit.b);
            
         NR_lower_power_fit=fit(NR_bins(NR_bins>1)',lower_bound(NR_bins>1)',g,'startpoint',[ 2.5 -.1]);
         NR_lower_fit=NR_lower_power_fit.a.*NR_bins.^(NR_lower_power_fit.b);

         NR_upper_power_fit=fit(NR_bins(NR_bins>1)',upper_bound(NR_bins>1)',g,'startpoint',[ 2.5 -.1]);
         NR_upper_fit=NR_upper_power_fit.a.*NR_bins.^(NR_upper_power_fit.b);
         
         NR_mean_fit_points=mean_NR_band;
         NR_lower_fit_points=lower_bound;
         NR_upper_fit_points=upper_bound;
                 
         NR_mean_fit_points_err=mean_NR_band_err;
         NR_lower_fit_points_err=lower_bound_err;
         NR_upper_fit_points_err=upper_bound_err;

         NR_bins_tracked=NR_bins;

  %%          
  %Make Plot    

 figure;
 hold on
 ylabel('Corrected log10(S2/S1)'); xlabel('S1 corrected (Pulse Area Phe)'); 
 title(sprintf('Run 04 NR Band, R<%d cm',max_r),'fontsize',16,'Interpreter','none');
 myfigview(16);
 ylim([0 3.5]);xlim([0 s1_max]);
 
 plot(DD_s1_all,log10(DD_s2_all./DD_s1_all),'.k');
 plot(NR_bins,NR_lower_fit,'--r','markersize',14,'linewidth',2)
 plot(NR_bins,NR_mean_fit,'-r','markersize',20,'linewidth',2)
 plot(NR_bins,NR_upper_fit,'--r','markersize',14,'linewidth',2)

%% Next, get the ER band in the same z-slices from the closest CH3T data set (Sep 2014)

clearvars -except NR_bins_tracked NR_lower_fit NR_mean_fit NR_upper_fit ...
         NR_mean_fit NR_lower_fit NR_upper_fit NR_mean_fit_points NR_lower_fit_points NR_upper_fit_points ...
         NR_mean_fit_points_err NR_lower_fit_points_err NR_upper_fit_points_err


load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p21\KrypCal2p21_QualityChecks\Discrimination\MappingParametersFit_Start_CH3T');

         %%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 % Removing field dependence from S1 and S2
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  %Hardcoded numbers from Energy Chi2 Fit
           
             
          
s1ab_s2_xyz_fit_map.p1=-0.4988; 
s1ab_s2_xyz_fit_map.p2=1.484;
c1=1-s1ab_s2_xyz_fit_map.p1.*2.7303.^2-2.7303.*s1ab_s2_xyz_fit_map.p2;
s1ab_s2_xyz_fit_map.p3=c1;

s1ab_s1_xyz_fit_map.p1=0.065;
s1ab_s1_xyz_fit_map.p2=0.02;
c2=1-s1ab_s1_xyz_fit_map.p1.*2.7303.^2-2.7303.*s1ab_s1_xyz_fit_map.p2;
s1ab_s1_xyz_fit_map.p3=c2;

 
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


 else %Only do 1D in this case
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

 
s2_phe_bottom=s2_phe_bottom./(s2_zfieldremoval);
s2_phe_both=s2_phe_both./(s2_zfieldremoval);
s1_phe_bottom=s1_phe_bottom./(s1_zfieldremoval);
s1_phe_both=s1_phe_both./(s1_zfieldremoval);



 end 


    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Calculating the S1 Z-dependence (Both PMT arrays)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    clear y_s1_both y_s1_bottom Fit_s1_both Fit_s1_bottom mean_s1_both mean_s1_both
    clear sigma_s1_both sigma_s1_both mean_index means mean_index means_error x temp_hist
    clear hist_s1_both x xfit S
    
    %set up energy bins for the gaussian fit
    bin_s1=8;%5
    s1_bin_min=100;
    s1_bin_max=500;
    x=s1_bin_min:bin_s1:s1_bin_max; %150:bin_s1:500;
    x=x';
        
    cut = inrange(s2radius,[rcut_min,rcut_max]) & inrange(s1_phe_both,[s1_bin_min,s1_bin_max]) & inrange(s2_width,[100,1000]) & inrange(s1_width,[30,400]);   
    
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
    
%     
%     light_correction_fig = figure;
%     errorbar(mean_index,means,means_error,'.k');
%     hold on     
%     [yplot_s1_both] = polyval(P_s1_both,dT_fit);
%     [center_s1_both, center_s1_both_err] = polyval(P_s1_both,z_center, S); %160[us]*1.51[\mus/us]
%     plot(dT_fit,yplot_s1_both,'-r','LineWidth',1);    
%     xlabel('Depth (\mus)', 'FontSize',16); ylabel('Mean (Pulse Area Phe)', 'FontSize', 16);
%     title(strcat('.Kr83m. S1 Mean vs. Z'), 'FontSize', 16,'Interpreter','none');
%     line([0 det_edge],[center_s1_both center_s1_both],'linestyle','--');
%     line([z_center; z_center;],[min(means) max(means)],'linestyle','--');
%     legend('S1\_both Mean', strcat('y=', num2str(P_s1_both(1),'%10.2e'), '*Z^2 + ', num2str(P_s1_both(2),'%10.2e')...
%         , '*Z + ', num2str(P_s1_both(3),4)),strcat('Det center = ',num2str(center_s1_both,4), ' [Phe]'),'location','northwest')
%     
%     box on;
%     myfigview(16);
%     xlim([0 det_edge]);
%     ylim([0.8*min(means) 1.15*max(means)]);

      
    s1_both_means=means;
    s1_both_means_sigma=means_error;
    s1_both_means_bin=mean_index;
    

    
    
       
    
%     s1_z_correction=polyval(P_s1_both,241.6)./polyval(P_s1_both,drift_distance_\mus);%normalize s1 to 160 ns.
    s1_z_correction=polyval(P_s1_both,z_center)./polyval(P_s1_both,drift_time);
    
    %s1_both_z_correction=interp1(s1_both_means_bin,s1_both_means,160)./interp1(s1_both_means_bin,s1_both_means,drift_time);
    s1_phe_both_z=s1_phe_both.*s1_z_correction;


%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   Calculating the S2 Z-dependence (i.e. electron lifetime)
%   Using Both PMT arrays
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
     
     
     %%calculate S2/S1 also, just to double check the lifetime from S2 only
     % hist_s2_both_s1(:,bin)= hist(s2_phe_both(time_cut)./s1_phe_both_z(time_cut),x3)';
     % temp_hist=hist_s2_both_s1(:,bin);   
      %Fit_s2_both_s1{bin}=fit(x3(2:end-1),temp_hist(2:end-1),'gauss1');%remove bin edges        
       % means_s2_s1(bin)=Fit_s2_both_s1{bin}.b1;
       % means_error_s2_s1(bin) = Fit_s2_both_s1{bin}.c1/2/sqrt(sum(hist_s2_both_s1(:,bin))); % == sigma/sqrt(N);
       % mean_index_s2_s1(bin) = (binStart+binEnd)/2;     

    end
    
 % //////////////////////////////Make Plots/////////////////////////////////
 
  
     dT_range= inrange( mean_index,[5+bin_size,500] );%from 50 to 300us for the fit


    s2_both_norm_z_index=mean_index(dT_range);
    s2_both_norm_z_means=means(dT_range);
    s2_at_top_both=RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,4); %normalize to right below the gate
    s2_both_norm_z=RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,4)./s2_both_norm_z_means;
  
    %%Tracker for pseudo-lifetime
    s2_at_bottom_both=RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,320);
    pseudo_lifetime_both=-320/log(s2_at_bottom_both/s2_at_top_both);
 
    dT_fit=0:1:det_edge;    
    
    Fit_EL= fit(mean_index(dT_range),means(dT_range),'exp1');
    yfitEL = Fit_EL.a.*exp(Fit_EL.b.*dT_fit);
    
    e_lifetime_both=-1/Fit_EL.b;
    s2_z0_both=Fit_EL.a;
       
    b=confint(Fit_EL, 0.683);%1 sigma matrix
    sigma_e_lifetime_both=(1/b(1,2)-1/b(2,2))/2;
    sigma_s2_z0_both=(b(2,1)-b(1,1))/2;
%  
%     lifetime_fig_both = figure;
%     errorbar(mean_index(2:end),means(2:end),means_error(2:end),'.k');  %first bin is anode, skip it
%     hold on     
%     plot(dT_fit,yfitEL,'-r','LineWidth',1);   
%     xlabel('Time (\mus)', 'FontSize',16); ylabel('Mean S2 (Pulse Area Phe)', 'FontSize', 16);
%     title(strcat('.Kr83m. S2_both Electron Lifetime'), 'FontSize', 16,'Interpreter','none');
%     legend('S2\_both Mean', strcat( '\lambda= ', num2str(e_lifetime_both,4), ' \pm ', num2str(sigma_e_lifetime_both,2), ' \mus' )...
%         ,'location','northeast');
%     myfigview(16);
%     xlim([0 det_edge]);   

    
        
    s2_both_means=means;
    s2_both_means_sigma=means_error;
    s2_both_means_bin=mean_index;    
    
    s2_phe_both_z=s2_phe_both.*RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z, drift_time);

    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Mapping S2 in XY (after correcting for Z-dep). Both PMT
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    clear x2 hist_bin

    %set up energy bins for the gaussian fit. Using both S2 only
   bin_s2=500;
   s2_bin_min=1000;
   s2_bin_max=40000;
   x2=s2_bin_min:bin_s2:s2_bin_max;
   x2=x2';

%set up matricies to be filled
mean_S2_both_xy = zeros(floor(s2_ybins),floor(s2_xbins));
Sigma_S2_both = zeros(floor(s2_ybins),floor(s2_xbins));
Count_S2_both = zeros(floor(s2_ybins),floor(s2_xbins));

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
       
      if Count_S2_both(y_count,x_count)>350;                         
                
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

color_range_max=max(max(mean_S2_both_xy));
color_range_min=min(min(mean_S2_both_xy(mean_S2_both_xy>0)));
vc_step=(color_range_max-color_range_min)/50;

vc=color_range_min:vc_step:color_range_max;
 
% s2_xy_both_fig = figure;
% contourf(s2xbins,s2ybins,mean_S2_both_xy,vc,'LineColor','none');
% xlabel('x (cm)','fontsize', 18);
% ylabel('y (cm)','fontsize', 18);
% title(strcat('. S2 Mean (Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
% caxis([color_range_min color_range_max])
% g=colorbar;
% ylabel(g,'phe','FontSize',16)
% set(g,'FontSize',18)
% hold on
% LUXPlotTopPMTs
% axis([-25 25 -25 25])
% myfigview(16);


%Calculate corrections matrix and 1 sigma corrections matrix
% if inpaint_on==1;  mean_S2_both_xy(mean_S2_both_xy==0)=nan; mean_S2_both_xy=inpaint_nans(mean_S2_both_xy,3); end;
center=interp2(s2xbins,s2ybins,mean_S2_both_xy,x_center,y_center,'cubic');%Normalize to the center (x=y=0)
norm_S2_all=center./mean_S2_both_xy; 
norm_S2_all(isinf(norm_S2_all))=1;%no infinity! Do not correct events outside r=25cm
norm_S2_all(isnan(norm_S2_all))=1;

% if inpaint_on==1;  Sigma_S2_both(Sigma_S2_both==0)=nan; Sigma_S2_both=inpaint_nans(Sigma_S2_both,3); end;
sigma_center=interp2(s2xbins,s2ybins,Sigma_S2_both,x_center,y_center,'cubic');%Normalize to the center (x=y=0)
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
mean_S1_both_xy = zeros(floor(s1_ybins),floor(s1_xbins));
Sigma_S1_both = zeros(floor(s1_ybins),floor(s1_xbins));
Count_S1_both = zeros(floor(s1_ybins),floor(s1_xbins));

% s1_z_correction=polyval(P_s1_both,241.6)./polyval(P_s1_both,drift_distance_\mus);%normalize s1 to 160 ns.
s1_z_correction=polyval(P_s1_both,z_center)./polyval(P_s1_both,drift_time);%normalize s1 to 160 ns.

s1_phe_both_z=s1_phe_both.*s1_z_correction;
% s2_phe_both_z=s2_phe_both.*exp(drift_distance_\mus/e_lifetime); %Correct S2 for lifeimte  %s2_phe_both_z[5000,18000]
% s2_phe_both_z=s2_phe_both.*exp(drift_time/e_lifetime_both); %Correct S2 for lifeimte  %s2_phe_both_z[5000,18000]
s2_phe_both_z=s2_phe_both.*RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z, drift_time);



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
       
      if Count_S1_both(y_count,x_count)>350; 
                        
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

color_range_max=max(max(mean_S1_both_xy));
color_range_min=min(min(mean_S1_both_xy(mean_S1_both_xy>0)));
vc_step=(color_range_max-color_range_min)/50;

vc=color_range_min:vc_step:color_range_max;
 
%  
% s1_xy_both_fig = figure;
% contourf(s1xbins,s1ybins,mean_S1_both_xy,vc,'LineColor','none');
% xlabel('x (cm)','fontsize', 18);
% ylabel('y (cm)','fontsize', 18);
% title(strcat('. S1 Mean (Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
% caxis([color_range_min color_range_max])
% g=colorbar;
% ylabel(g,'phe','FontSize',16)
% set(g,'FontSize',18)
% hold on
% LUXPlotTopPMTs
% axis([-25 25 -25 25])
% myfigview(16);


%Calculate corrections matrix and 1 sigma corrections matrix
% if inpaint_on==1;  mean_S1_both_xy(mean_S1_both_xy==0)=nan; mean_S1_both_xy=inpaint_nans(mean_S1_both_xy,3); end;
center=interp2(s1xbins,s1ybins,mean_S1_both_xy,x_center,y_center,'cubic');%Normalize to the center (x=y=0)
norm_S1_all=center./mean_S1_both_xy; 
norm_S1_all(isinf(norm_S1_all))=1;%no infinity! Do not correct events outside r=25cm
norm_S1_all(isnan(norm_S1_all))=1;

% if inpaint_on==1;  Sigma_S1_both(Sigma_S1_both==0)=nan; Sigma_S1_both=inpaint_nans(Sigma_S1_both,3); end;
sigma_center=interp2(s1xbins,s1ybins,Sigma_S1_both,x_center,y_center,'cubic');%Normalize to the center (x=y=0)
sigma_norm_S1_all=sqrt((sigma_center./mean_S1_both_xy).^2+(Sigma_S1_both.*center./mean_S1_both_xy.^2).^2);
sigma_norm_S1_all(isinf(sigma_norm_S1_all))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
sigma_norm_S1_all(isnan(sigma_norm_S1_all))=1;%no nan. Set Sigma=1, which is 100% uncertainty

s1_both_center_xyz=center;
s1_both_center_xyz_error=sigma_center;

%% Apply corrections to ER data and calculate the ERband

load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\Run04_CH3T\KrypCal2p18_IQs\MakePlots\CH3T_Sep2014_Kr2p18.mat');
clear s1_phe_both_xyz s2_phe_both_xyz s1_phe_both_z s2_phe_both_z

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


%------------------------------------------
% Calculating corrections
%------------------------------------------

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
 
%--------------------------------------------------------------------------
% Calculate XY corrections
%--------------------------------------------------------------------------

% Calculate S1 XY-only corrections
s1xy_correction = interp2(s1_x_bins,s1_y_bins,s1_map_all,s1_x_cm,s1_y_cm,'cubic');
s1xy_correction(find(s1xy_correction==0))=1.0;  s1xy_correction(find(isnan(s1xy_correction)))=1.0;

% Calculate S2 XY corrections
s2xy_correction = interp2(s2_x_bins,s2_y_bins,s2_map_all,s2x,s2y,'cubic');
s2xy_correction(find(s2xy_correction==0))=1.0;  s2xy_correction(find(isnan(s2xy_correction)))=1.0;

%--------------------------------------------------------------------------
% Calculate XYZ corrections
%--------------------------------------------------------------------------

s1xyz_correction=s1xy_correction.*s1_z_correction;

%--------------------------------------------------------------------------
% Apply corrections
%--------------------------------------------------------------------------

s1_phe_both_z=s1_phe_both.*s1_z_correction;
s2_phe_both_z=s2_phe_both.*s2_z_correction;

% Calculate Energy spectra
s1_phe_both_xyz=s1_phe_both.*s1xyz_correction;
s2_phe_both_xyz=s2_phe_both.*s2_z_correction.*s2xy_correction;


%% Make ER Band
max_r=18;
s1_min=0;
s1_max=50;


pulse_area_cut= (log10(s1_phe_both)+0.7*log10(s2_phe_both)<4.85) & (log10(s1_phe_both)+0.7*log10(s2_phe_both)>2) & ...
    (log10(s1_phe_both)-2*log10(s2_phe_both_xyz)<-4.7) & (log10(s1_phe_both)-2*log10(s2_phe_both)>-6.7);
fiducial_cut=inrange(drift_time,[35 75]) & inrange(s2radius_del,[0,max_r]);

tritium_cut= clean_cut & pulse_area_cut  & fiducial_cut;

bin_counter=1;
  


%% Make Band Plot
   band_sigma=1.28; %90% confidence bounds
    
             
            s2_s1_bin= 1:0.05:5;
           ER_bin_size=3; %3 works
           index_i=1;
           ER_bins=ER_bin_size/2+s1_min:ER_bin_size:s1_max;
               lower_bound=ones(size(ER_bins));
               mean_ER_band=ones(size(ER_bins));
               upper_bound=ones(size(ER_bins));
                counts_ER_hist=ones(size(ER_bins));
               ER_sigma=ones(size(ER_bins));
               ER_sigma_sigma=ones(size(ER_bins));
        
%             temporary, for calculating correct E-lifetime   
            clear mean_ER_band_err lower_bound_err upper_bound_err lower_bound mean_ER_band upper_bound  ...
                gauss_sig_list_err gauss_amp_list gauss_amp_list_err
            
                     %Refit the NR band - ER binning is same as NR binning
             NR_mean_data=NR_mean_fit_points;
             NR_mean_data_err=NR_mean_fit_points_err;
                  
            
        for ER_bin = (ER_bin_size/2)+s1_min:ER_bin_size:(s1_max); %HIST uses bin centers
   
        
            ER_fit_cut=tritium_cut & inrange(s1_phe_both_xyz,[ER_bin-ER_bin_size/2,ER_bin+ER_bin_size/2]);
             
            ER_band_hist=hist(log10(s2_phe_both_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            counts_ER_hist(index_i)=sum(ER_band_hist);
           
            Fit_ER=fit(s2_s1_bin',ER_band_hist','gauss1');
             
            rms_fit=fit_ML_normal(log10(s2_phe_both_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            
            %band=confint(Fit_ER,.1);
%             ER_sigma(index_i)= Fit_ER.c1/sqrt(2);
            ER_sigma(index_i)= sqrt(rms_fit.sig2);
            ER_sigma_sigma(index_i)=ER_sigma(index_i)/sqrt(2*counts_ER_hist(index_i));
            
            lower_bound(index_i)=Fit_ER.b1-ER_sigma(index_i)*band_sigma;
            mean_ER_band(index_i)=Fit_ER.b1;
            upper_bound(index_i)=Fit_ER.b1+ER_sigma(index_i)*band_sigma;


            ER_bin_count(index_i)=length(s1_phe_both_xyz(ER_fit_cut));
            ER_bin_belowNR_count(index_i)=length(s1_phe_both_xyz(ER_fit_cut & log10(s2_phe_both_xyz./s1_phe_both_xyz)<NR_mean_data(index_i)));
            
            temp_err=confint(Fit_ER,0.68);            
            mean_ER_band_err(index_i)=abs(temp_err(1,2)-Fit_ER.b1);
            lower_bound_err(index_i)=sqrt(mean_ER_band_err(index_i).^2 + (band_sigma.*abs(temp_err(1,3)-Fit_ER.c1)/sqrt(2)).^2);
            upper_bound_err(index_i)=sqrt(mean_ER_band_err(index_i).^2 + (band_sigma.*abs(temp_err(1,3)-Fit_ER.c1)/sqrt(2)).^2);

            amp(index_i)=Fit_ER.a1;
            amp_err(index_i)=abs(temp_err(1,1)-Fit_ER.a1);
            
            gauss_sig_list_err(index_i)=abs(temp_err(1,3) - Fit_ER.c1)./sqrt(2);
            gauss_amp_list(index_i)=Fit_ER.a1;
            gauss_amp_list_err(index_i)=abs(temp_err(1,1)-Fit_ER.a1);
            %mean_ER_band_sim(index_i)=Fit_ER_sim.b1;
            
            index_i=index_i+1;
         end
        
         
         % Fit power law to the mean
         g=fittype('a*x^b');
         ER_mean_power_fit=fit(ER_bins(ER_bins>1)',mean_ER_band(ER_bins>1)',g,'startpoint',[ 2.5 -.1]);
         ER_mean_fit{bin_counter}=ER_mean_power_fit.a.*ER_bins.^(ER_mean_power_fit.b);
            
         ER_lower_power_fit=fit(ER_bins(ER_bins>1)',lower_bound(ER_bins>1)',g,'startpoint',[ 2.5 -.1]);
         ER_lower_fit{bin_counter}=ER_lower_power_fit.a.*ER_bins.^(ER_lower_power_fit.b);

         ER_upper_power_fit=fit(ER_bins(ER_bins>1)',upper_bound(ER_bins>1)',g,'startpoint',[ 2.5 -.1]);
         ER_upper_fit{bin_counter}=ER_upper_power_fit.a.*ER_bins.^(ER_upper_power_fit.b);
         
         ER_mean_fit_points{bin_counter}=mean_ER_band;
         ER_lower_fit_points{bin_counter}=lower_bound;
         ER_upper_fit_points{bin_counter}=upper_bound;
                 
         ER_mean_fit_points_err{bin_counter}=mean_ER_band_err;
         ER_lower_fit_points_err{bin_counter}=lower_bound_err;
         ER_upper_fit_points_err{bin_counter}=upper_bound_err;

         ER_bins_tracked{bin_counter}=ER_bins;
         
         
         

             NR_bin_data=NR_bins_tracked;
             NR_mean_refit=fit(NR_bin_data(NR_bin_data>1)',NR_mean_data(NR_bin_data>1)',g,'startpoint',[ 2.5 -.1]);
         
             
         %Gaussian Assumption            
         clear gauss_mean gauss_sig gauss_amp gauss_amp_err gauss_mean_err gauss_sig_err leakage_frac_err leakage_frac NR_mean
         for i=1:length(mean_ER_band)
             
             NR_mean=NR_mean_data(i);            
             gauss_mean= mean_ER_band(i);
             gauss_sig= ER_sigma(i);
             gauss_amp=amp(i); %estimate gaussian amplitude from the individual bins
             gauss_amp_err=amp_err(i); %estimate gaussian amplitude error from the individual bins
             gauss_mean_err=mean_ER_band_err(i); %estimate mean err from the individual bins
             gauss_sig_err=ER_sigma_sigma(i); %estimate sigma error from the individual bins

         %Using solution of definite integral so that I can propogate errors with analytic approach
         total_count=sqrt(pi./2).*gauss_amp.*gauss_sig.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) + 1 );
         leakage_count=sqrt(pi./2).*gauss_amp.*gauss_sig.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) - erf((gauss_mean-NR_mean)./(sqrt(2).*gauss_sig)));
         
         %derivative w.r.t amplitude
         total_count_damp = sqrt(pi/2).*gauss_sig.*(erf(gauss_mean/(sqrt(2)*gauss_sig)) +1);
         leakage_count_damp= sqrt(pi/2).*gauss_sig.*(erf(gauss_mean/(sqrt(2)*gauss_sig)) - erf((gauss_mean - NR_mean)/(sqrt(2)*gauss_sig)));
         
         %derivative w.r.t mean
         total_count_dmean = gauss_amp.*exp(-(gauss_mean.^2)./(2.*gauss_sig.^2));
         leakage_count_dmean = sqrt(pi/2).*gauss_amp.*gauss_sig.*( ( sqrt(2/pi).*exp(-(gauss_mean.^2)./(2.*gauss_sig.^2)))./gauss_sig - sqrt(2/pi).*exp(-((gauss_mean-NR_mean).^2)./(2.*gauss_sig.^2))./gauss_sig);
         
         %derivative w.r.t sigma
         total_count_dsig=sqrt(pi/2).*gauss_amp.* (erf ( gauss_mean./(sqrt(2).*gauss_sig)) + 1) - (gauss_amp .* gauss_mean .* exp( - (gauss_mean.^2)./(2.*gauss_sig.^2)))./gauss_sig;
         leakage_count_dsig=sqrt(pi/2).*gauss_amp.*gauss_sig.*( (sqrt(2/pi).*(gauss_mean-NR_mean).*exp( (-(gauss_mean-NR_mean).^2)./(2.*gauss_sig.^2)))./(gauss_sig.^2) ...
             - (sqrt(2/pi).*(gauss_mean).*exp( (-(gauss_mean).^2)./(2.*gauss_sig.^2)))./(gauss_sig.^2)) + ...
             sqrt(pi/2).*gauss_amp.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) - erf((gauss_mean-NR_mean)./(sqrt(2).*gauss_sig)));
         
         total_count_err=sqrt( (total_count_damp.*gauss_amp_err).^2 + (total_count_dmean.*gauss_mean_err).^2 + (total_count_dsig.*gauss_sig_err).^2);
         leakage_count_err=sqrt( (leakage_count_damp.*gauss_amp_err).^2 + (leakage_count_dmean.*gauss_mean_err).^2 + (leakage_count_dsig.*gauss_sig_err).^2);
         
         leakage_frac(i)=leakage_count./total_count;
         leakage_frac_err(i)=sqrt( (leakage_count_err./total_count).^2 + (leakage_count.*total_count_err./(total_count.^2)).^2);

         end


         gaussian_leakage_frac_tracker{bin_counter}=leakage_frac;
         gaussian_leakage_frac_err_tracker{bin_counter}=leakage_frac_err;         
              

        %Counting Data Approach          
         clear gauss_mean gauss_sig gauss_amp gauss_amp_err gauss_mean_err gauss_sig_err leakage_frac_err leakage_frac NR_mean
         for i=1:length(mean_ER_band)

         %Number of events in S1 bin
         total_count=ER_bin_count(i);
         leakage_count=ER_bin_belowNR_count(i);
  
         total_count_err=sqrt(total_count);
         leakage_count_err=sqrt(leakage_count);
         
         leakage_frac(i)=leakage_count./total_count;
         leakage_frac_err(i)=sqrt( (leakage_count_err./total_count).^2 + (leakage_count.*total_count_err./(total_count.^2)).^2);

         end

        
         counting_leakage_frac_tracker{bin_counter}=leakage_frac;
         counting_leakage_frac_err_tracker{bin_counter}=leakage_frac_err;            

         %Based on ER Fit Approach
         i=1;
         clear gauss_mean gauss_sig gauss_amp gauss_amp_err gauss_mean_err gauss_sig_err leakage_frac_err leakage_frac NR_mean
         for S1_val=1:1:49
             
                NR_mean=NR_mean_refit.a.*S1_val.^(NR_mean_refit.b);
             
             gauss_mean= ER_mean_power_fit.a.*S1_val.^(ER_mean_power_fit.b);
             gauss_sig= abs((ER_upper_power_fit.a.*S1_val.^(ER_upper_power_fit.b)) - gauss_mean)/1.28;
             gauss_amp=interp1(ER_bins,gauss_amp_list,S1_val); %estimate gaussian amplitude from the individual bins
             gauss_amp_err=interp1(ER_bins,gauss_amp_list_err,S1_val); %estimate gaussian amplitude error from the individual bins
             gauss_mean_err=interp1(ER_bins,mean_ER_band_err,S1_val); %estimate mean err from the individual bins
             gauss_sig_err=interp1(ER_bins,gauss_sig_list_err,S1_val); %estimate sigma error from the individual bins

         %Using solution of definite integral so that I can propogate errors with analytic approach
         total_count=sqrt(pi./2).*gauss_amp.*gauss_sig.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) + 1 );
         leakage_count=sqrt(pi./2).*gauss_amp.*gauss_sig.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) - erf((gauss_mean-NR_mean)./(sqrt(2).*gauss_sig)));
         
         %derivative w.r.t amplitude
         total_count_damp = sqrt(pi/2).*gauss_sig.*(erf(gauss_mean/(sqrt(2)*gauss_sig)) +1);
         leakage_count_damp= sqrt(pi/2).*gauss_sig.*(erf(gauss_mean/(sqrt(2)*gauss_sig)) - erf((gauss_mean - NR_mean)/(sqrt(2)*gauss_sig)));
         
         %derivative w.r.t mean
         total_count_dmean = gauss_amp.*exp(-(gauss_mean.^2)./(2.*gauss_sig.^2));
         leakage_count_dmean = sqrt(pi/2).*gauss_amp.*gauss_sig.*( ( sqrt(2/pi).*exp(-(gauss_mean.^2)./(2.*gauss_sig.^2)))./gauss_sig - sqrt(2/pi).*exp(-((gauss_mean-NR_mean).^2)./(2.*gauss_sig.^2))./gauss_sig);
         
         %derivative w.r.t sigma
         total_count_dsig=sqrt(pi/2).*gauss_amp.* (erf ( gauss_mean./(sqrt(2).*gauss_sig)) + 1) - (gauss_amp .* gauss_mean .* exp( - (gauss_mean.^2)./(2.*gauss_sig.^2)))./gauss_sig;
         leakage_count_dsig=sqrt(pi/2).*gauss_amp.*gauss_sig.*( (sqrt(2/pi).*(gauss_mean-NR_mean).*exp( (-(gauss_mean-NR_mean).^2)./(2.*gauss_sig.^2)))./(gauss_sig.^2) ...
             - (sqrt(2/pi).*(gauss_mean).*exp( (-(gauss_mean).^2)./(2.*gauss_sig.^2)))./(gauss_sig.^2)) + ...
             sqrt(pi/2).*gauss_amp.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) - erf((gauss_mean-NR_mean)./(sqrt(2).*gauss_sig)));
         
         total_count_err=sqrt( (total_count_damp.*gauss_amp_err).^2 + (total_count_dmean.*gauss_mean_err).^2 + (total_count_dsig.*gauss_sig_err).^2);
         leakage_count_err=sqrt( (leakage_count_damp.*gauss_amp_err).^2 + (leakage_count_dmean.*gauss_mean_err).^2 + (leakage_count_dsig.*gauss_sig_err).^2);
         
         leakage_frac(i)=leakage_count./total_count;
         leakage_frac_err(i)=sqrt( (leakage_count_err./total_count).^2 + (leakage_count.*total_count_err./(total_count.^2)).^2);
         i=i+1;
        
         end
         
         Fit_leakage_frac_tracker{bin_counter}=leakage_frac;
         Fit_leakage_frac_err_tracker{bin_counter}=leakage_frac_err; 
         
         bin_counter=bin_counter+1;
         
         %% Move onto next Z slice
max_r=18;
s1_min=0;
s1_max=50;


pulse_area_cut= (log10(s1_phe_both)+0.7*log10(s2_phe_both)<4.85) & (log10(s1_phe_both)+0.7*log10(s2_phe_both)>2) & ...
    (log10(s1_phe_both)-2*log10(s2_phe_both_xyz)<-4.7) & (log10(s1_phe_both)-2*log10(s2_phe_both)>-6.7);
fiducial_cut=inrange(drift_time,[75 125]) & inrange(s2radius_del,[0,max_r]);

tritium_cut= clean_cut & pulse_area_cut  & fiducial_cut;

  


%% Make Band Plot
   band_sigma=1.28; %90% confidence bounds
    
             
            s2_s1_bin= 1:0.05:5;
           ER_bin_size=3; %3 works
           index_i=1;
           ER_bins=ER_bin_size/2+s1_min:ER_bin_size:s1_max;
               lower_bound=ones(size(ER_bins));
               mean_ER_band=ones(size(ER_bins));
               upper_bound=ones(size(ER_bins));
                counts_ER_hist=ones(size(ER_bins));
               ER_sigma=ones(size(ER_bins));
               ER_sigma_sigma=ones(size(ER_bins));
        
%             temporary, for calculating correct E-lifetime   
            clear mean_ER_band_err lower_bound_err upper_bound_err lower_bound mean_ER_band upper_bound  ...
                gauss_sig_list_err gauss_amp_list gauss_amp_list_err
            
                     %Refit the NR band - ER binning is same as NR binning
             NR_mean_data=NR_mean_fit_points;
             NR_mean_data_err=NR_mean_fit_points_err;
                  
            
        for ER_bin = (ER_bin_size/2)+s1_min:ER_bin_size:(s1_max); %HIST uses bin centers
   
        
            ER_fit_cut=tritium_cut & inrange(s1_phe_both_xyz,[ER_bin-ER_bin_size/2,ER_bin+ER_bin_size/2]);
             
            ER_band_hist=hist(log10(s2_phe_both_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            counts_ER_hist(index_i)=sum(ER_band_hist);
           
            Fit_ER=fit(s2_s1_bin',ER_band_hist','gauss1');
             
            rms_fit=fit_ML_normal(log10(s2_phe_both_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            
            %band=confint(Fit_ER,.1);
%             ER_sigma(index_i)= Fit_ER.c1/sqrt(2);
            ER_sigma(index_i)= sqrt(rms_fit.sig2);
            ER_sigma_sigma(index_i)=ER_sigma(index_i)/sqrt(2*counts_ER_hist(index_i));
            
            lower_bound(index_i)=Fit_ER.b1-ER_sigma(index_i)*band_sigma;
            mean_ER_band(index_i)=Fit_ER.b1;
            upper_bound(index_i)=Fit_ER.b1+ER_sigma(index_i)*band_sigma;


            ER_bin_count(index_i)=length(s1_phe_both_xyz(ER_fit_cut));
            ER_bin_belowNR_count(index_i)=length(s1_phe_both_xyz(ER_fit_cut & log10(s2_phe_both_xyz./s1_phe_both_xyz)<NR_mean_data(index_i)));
            
            temp_err=confint(Fit_ER,0.68);            
            mean_ER_band_err(index_i)=abs(temp_err(1,2)-Fit_ER.b1);
            lower_bound_err(index_i)=sqrt(mean_ER_band_err(index_i).^2 + (band_sigma.*abs(temp_err(1,3)-Fit_ER.c1)/sqrt(2)).^2);
            upper_bound_err(index_i)=sqrt(mean_ER_band_err(index_i).^2 + (band_sigma.*abs(temp_err(1,3)-Fit_ER.c1)/sqrt(2)).^2);

            amp(index_i)=Fit_ER.a1;
            amp_err(index_i)=abs(temp_err(1,1)-Fit_ER.a1);
            
            gauss_sig_list_err(index_i)=abs(temp_err(1,3) - Fit_ER.c1)./sqrt(2);
            gauss_amp_list(index_i)=Fit_ER.a1;
            gauss_amp_list_err(index_i)=abs(temp_err(1,1)-Fit_ER.a1);
            %mean_ER_band_sim(index_i)=Fit_ER_sim.b1;
            
            index_i=index_i+1;
         end
        
         
         % Fit power law to the mean
         g=fittype('a*x^b');
         ER_mean_power_fit=fit(ER_bins(ER_bins>1)',mean_ER_band(ER_bins>1)',g,'startpoint',[ 2.5 -.1]);
         ER_mean_fit{bin_counter}=ER_mean_power_fit.a.*ER_bins.^(ER_mean_power_fit.b);
            
         ER_lower_power_fit=fit(ER_bins(ER_bins>1)',lower_bound(ER_bins>1)',g,'startpoint',[ 2.5 -.1]);
         ER_lower_fit{bin_counter}=ER_lower_power_fit.a.*ER_bins.^(ER_lower_power_fit.b);

         ER_upper_power_fit=fit(ER_bins(ER_bins>1)',upper_bound(ER_bins>1)',g,'startpoint',[ 2.5 -.1]);
         ER_upper_fit{bin_counter}=ER_upper_power_fit.a.*ER_bins.^(ER_upper_power_fit.b);
         
         ER_mean_fit_points{bin_counter}=mean_ER_band;
         ER_lower_fit_points{bin_counter}=lower_bound;
         ER_upper_fit_points{bin_counter}=upper_bound;
                 
         ER_mean_fit_points_err{bin_counter}=mean_ER_band_err;
         ER_lower_fit_points_err{bin_counter}=lower_bound_err;
         ER_upper_fit_points_err{bin_counter}=upper_bound_err;

         ER_bins_tracked{bin_counter}=ER_bins;
         
         
         

             NR_bin_data=NR_bins_tracked;
             NR_mean_refit=fit(NR_bin_data(NR_bin_data>1)',NR_mean_data(NR_bin_data>1)',g,'startpoint',[ 2.5 -.1]);
         
             
         %Gaussian Assumption            
         clear gauss_mean gauss_sig gauss_amp gauss_amp_err gauss_mean_err gauss_sig_err leakage_frac_err leakage_frac NR_mean
         for i=1:length(mean_ER_band)
             
             NR_mean=NR_mean_data(i);            
             gauss_mean= mean_ER_band(i);
             gauss_sig= ER_sigma(i);
             gauss_amp=amp(i); %estimate gaussian amplitude from the individual bins
             gauss_amp_err=amp_err(i); %estimate gaussian amplitude error from the individual bins
             gauss_mean_err=mean_ER_band_err(i); %estimate mean err from the individual bins
             gauss_sig_err=ER_sigma_sigma(i); %estimate sigma error from the individual bins

         %Using solution of definite integral so that I can propogate errors with analytic approach
         total_count=sqrt(pi./2).*gauss_amp.*gauss_sig.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) + 1 );
         leakage_count=sqrt(pi./2).*gauss_amp.*gauss_sig.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) - erf((gauss_mean-NR_mean)./(sqrt(2).*gauss_sig)));
         
         %derivative w.r.t amplitude
         total_count_damp = sqrt(pi/2).*gauss_sig.*(erf(gauss_mean/(sqrt(2)*gauss_sig)) +1);
         leakage_count_damp= sqrt(pi/2).*gauss_sig.*(erf(gauss_mean/(sqrt(2)*gauss_sig)) - erf((gauss_mean - NR_mean)/(sqrt(2)*gauss_sig)));
         
         %derivative w.r.t mean
         total_count_dmean = gauss_amp.*exp(-(gauss_mean.^2)./(2.*gauss_sig.^2));
         leakage_count_dmean = sqrt(pi/2).*gauss_amp.*gauss_sig.*( ( sqrt(2/pi).*exp(-(gauss_mean.^2)./(2.*gauss_sig.^2)))./gauss_sig - sqrt(2/pi).*exp(-((gauss_mean-NR_mean).^2)./(2.*gauss_sig.^2))./gauss_sig);
         
         %derivative w.r.t sigma
         total_count_dsig=sqrt(pi/2).*gauss_amp.* (erf ( gauss_mean./(sqrt(2).*gauss_sig)) + 1) - (gauss_amp .* gauss_mean .* exp( - (gauss_mean.^2)./(2.*gauss_sig.^2)))./gauss_sig;
         leakage_count_dsig=sqrt(pi/2).*gauss_amp.*gauss_sig.*( (sqrt(2/pi).*(gauss_mean-NR_mean).*exp( (-(gauss_mean-NR_mean).^2)./(2.*gauss_sig.^2)))./(gauss_sig.^2) ...
             - (sqrt(2/pi).*(gauss_mean).*exp( (-(gauss_mean).^2)./(2.*gauss_sig.^2)))./(gauss_sig.^2)) + ...
             sqrt(pi/2).*gauss_amp.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) - erf((gauss_mean-NR_mean)./(sqrt(2).*gauss_sig)));
         
         total_count_err=sqrt( (total_count_damp.*gauss_amp_err).^2 + (total_count_dmean.*gauss_mean_err).^2 + (total_count_dsig.*gauss_sig_err).^2);
         leakage_count_err=sqrt( (leakage_count_damp.*gauss_amp_err).^2 + (leakage_count_dmean.*gauss_mean_err).^2 + (leakage_count_dsig.*gauss_sig_err).^2);
         
         leakage_frac(i)=leakage_count./total_count;
         leakage_frac_err(i)=sqrt( (leakage_count_err./total_count).^2 + (leakage_count.*total_count_err./(total_count.^2)).^2);

         end


         gaussian_leakage_frac_tracker{bin_counter}=leakage_frac;
         gaussian_leakage_frac_err_tracker{bin_counter}=leakage_frac_err;         
              

        %Counting Data Approach          
         clear gauss_mean gauss_sig gauss_amp gauss_amp_err gauss_mean_err gauss_sig_err leakage_frac_err leakage_frac NR_mean
         for i=1:length(mean_ER_band)

         %Number of events in S1 bin
         total_count=ER_bin_count(i);
         leakage_count=ER_bin_belowNR_count(i);
  
         total_count_err=sqrt(total_count);
         leakage_count_err=sqrt(leakage_count);
         
         leakage_frac(i)=leakage_count./total_count;
         leakage_frac_err(i)=sqrt( (leakage_count_err./total_count).^2 + (leakage_count.*total_count_err./(total_count.^2)).^2);

         end


         counting_leakage_frac_tracker{bin_counter}=leakage_frac;
         counting_leakage_frac_err_tracker{bin_counter}=leakage_frac_err;            
         %Based on ER Fit Approach
         i=1;
         clear gauss_mean gauss_sig gauss_amp gauss_amp_err gauss_mean_err gauss_sig_err leakage_frac_err leakage_frac NR_mean
         for S1_val=1:1:49
             
                NR_mean=NR_mean_refit.a.*S1_val.^(NR_mean_refit.b);
             
             gauss_mean= ER_mean_power_fit.a.*S1_val.^(ER_mean_power_fit.b);
             gauss_sig= abs((ER_upper_power_fit.a.*S1_val.^(ER_upper_power_fit.b)) - gauss_mean)/1.28;
             gauss_amp=interp1(ER_bins,gauss_amp_list,S1_val); %estimate gaussian amplitude from the individual bins
             gauss_amp_err=interp1(ER_bins,gauss_amp_list_err,S1_val); %estimate gaussian amplitude error from the individual bins
             gauss_mean_err=interp1(ER_bins,mean_ER_band_err,S1_val); %estimate mean err from the individual bins
             gauss_sig_err=interp1(ER_bins,gauss_sig_list_err,S1_val); %estimate sigma error from the individual bins

         %Using solution of definite integral so that I can propogate errors with analytic approach
         total_count=sqrt(pi./2).*gauss_amp.*gauss_sig.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) + 1 );
         leakage_count=sqrt(pi./2).*gauss_amp.*gauss_sig.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) - erf((gauss_mean-NR_mean)./(sqrt(2).*gauss_sig)));
         
         %derivative w.r.t amplitude
         total_count_damp = sqrt(pi/2).*gauss_sig.*(erf(gauss_mean/(sqrt(2)*gauss_sig)) +1);
         leakage_count_damp= sqrt(pi/2).*gauss_sig.*(erf(gauss_mean/(sqrt(2)*gauss_sig)) - erf((gauss_mean - NR_mean)/(sqrt(2)*gauss_sig)));
         
         %derivative w.r.t mean
         total_count_dmean = gauss_amp.*exp(-(gauss_mean.^2)./(2.*gauss_sig.^2));
         leakage_count_dmean = sqrt(pi/2).*gauss_amp.*gauss_sig.*( ( sqrt(2/pi).*exp(-(gauss_mean.^2)./(2.*gauss_sig.^2)))./gauss_sig - sqrt(2/pi).*exp(-((gauss_mean-NR_mean).^2)./(2.*gauss_sig.^2))./gauss_sig);
         
         %derivative w.r.t sigma
         total_count_dsig=sqrt(pi/2).*gauss_amp.* (erf ( gauss_mean./(sqrt(2).*gauss_sig)) + 1) - (gauss_amp .* gauss_mean .* exp( - (gauss_mean.^2)./(2.*gauss_sig.^2)))./gauss_sig;
         leakage_count_dsig=sqrt(pi/2).*gauss_amp.*gauss_sig.*( (sqrt(2/pi).*(gauss_mean-NR_mean).*exp( (-(gauss_mean-NR_mean).^2)./(2.*gauss_sig.^2)))./(gauss_sig.^2) ...
             - (sqrt(2/pi).*(gauss_mean).*exp( (-(gauss_mean).^2)./(2.*gauss_sig.^2)))./(gauss_sig.^2)) + ...
             sqrt(pi/2).*gauss_amp.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) - erf((gauss_mean-NR_mean)./(sqrt(2).*gauss_sig)));
         
         total_count_err=sqrt( (total_count_damp.*gauss_amp_err).^2 + (total_count_dmean.*gauss_mean_err).^2 + (total_count_dsig.*gauss_sig_err).^2);
         leakage_count_err=sqrt( (leakage_count_damp.*gauss_amp_err).^2 + (leakage_count_dmean.*gauss_mean_err).^2 + (leakage_count_dsig.*gauss_sig_err).^2);
         
         leakage_frac(i)=leakage_count./total_count;
         leakage_frac_err(i)=sqrt( (leakage_count_err./total_count).^2 + (leakage_count.*total_count_err./(total_count.^2)).^2);
         i=i+1;
        
         end
         
         Fit_leakage_frac_tracker{bin_counter}=leakage_frac;
         Fit_leakage_frac_err_tracker{bin_counter}=leakage_frac_err; 
        
         bin_counter=bin_counter+1;
         
          %% Move onto next Z slice
max_r=18;
s1_min=0;
s1_max=50;


pulse_area_cut= (log10(s1_phe_both)+0.7*log10(s2_phe_both)<4.85) & (log10(s1_phe_both)+0.7*log10(s2_phe_both)>2) & ...
    (log10(s1_phe_both)-2*log10(s2_phe_both_xyz)<-4.7) & (log10(s1_phe_both)-2*log10(s2_phe_both)>-6.7);
fiducial_cut=inrange(drift_time,[125 175]) & inrange(s2radius_del,[0,max_r]);

tritium_cut= clean_cut & pulse_area_cut  & fiducial_cut;

  


%% Make Band Plot
   band_sigma=1.28; %90% confidence bounds
    
             
            s2_s1_bin= 1:0.05:5;
           ER_bin_size=3; %3 works
           index_i=1;
           ER_bins=ER_bin_size/2+s1_min:ER_bin_size:s1_max;
               lower_bound=ones(size(ER_bins));
               mean_ER_band=ones(size(ER_bins));
               upper_bound=ones(size(ER_bins));
                counts_ER_hist=ones(size(ER_bins));
               ER_sigma=ones(size(ER_bins));
               ER_sigma_sigma=ones(size(ER_bins));
        
%             temporary, for calculating correct E-lifetime   
            clear mean_ER_band_err lower_bound_err upper_bound_err lower_bound mean_ER_band upper_bound  ...
                gauss_sig_list_err gauss_amp_list gauss_amp_list_err
            
                     %Refit the NR band - ER binning is same as NR binning
             NR_mean_data=NR_mean_fit_points;
             NR_mean_data_err=NR_mean_fit_points_err;
                  
            
        for ER_bin = (ER_bin_size/2)+s1_min:ER_bin_size:(s1_max); %HIST uses bin centers
   
        
            ER_fit_cut=tritium_cut & inrange(s1_phe_both_xyz,[ER_bin-ER_bin_size/2,ER_bin+ER_bin_size/2]);
             
            ER_band_hist=hist(log10(s2_phe_both_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            counts_ER_hist(index_i)=sum(ER_band_hist);
           
            Fit_ER=fit(s2_s1_bin',ER_band_hist','gauss1');
             
            rms_fit=fit_ML_normal(log10(s2_phe_both_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            
            %band=confint(Fit_ER,.1);
%             ER_sigma(index_i)= Fit_ER.c1/sqrt(2);
            ER_sigma(index_i)= sqrt(rms_fit.sig2);
            ER_sigma_sigma(index_i)=ER_sigma(index_i)/sqrt(2*counts_ER_hist(index_i));
            
            lower_bound(index_i)=Fit_ER.b1-ER_sigma(index_i)*band_sigma;
            mean_ER_band(index_i)=Fit_ER.b1;
            upper_bound(index_i)=Fit_ER.b1+ER_sigma(index_i)*band_sigma;


            ER_bin_count(index_i)=length(s1_phe_both_xyz(ER_fit_cut));
            ER_bin_belowNR_count(index_i)=length(s1_phe_both_xyz(ER_fit_cut & log10(s2_phe_both_xyz./s1_phe_both_xyz)<NR_mean_data(index_i)));
            
            temp_err=confint(Fit_ER,0.68);            
            mean_ER_band_err(index_i)=abs(temp_err(1,2)-Fit_ER.b1);
            lower_bound_err(index_i)=sqrt(mean_ER_band_err(index_i).^2 + (band_sigma.*abs(temp_err(1,3)-Fit_ER.c1)/sqrt(2)).^2);
            upper_bound_err(index_i)=sqrt(mean_ER_band_err(index_i).^2 + (band_sigma.*abs(temp_err(1,3)-Fit_ER.c1)/sqrt(2)).^2);

            amp(index_i)=Fit_ER.a1;
            amp_err(index_i)=abs(temp_err(1,1)-Fit_ER.a1);
            
            gauss_sig_list_err(index_i)=abs(temp_err(1,3) - Fit_ER.c1)./sqrt(2);
            gauss_amp_list(index_i)=Fit_ER.a1;
            gauss_amp_list_err(index_i)=abs(temp_err(1,1)-Fit_ER.a1);
            %mean_ER_band_sim(index_i)=Fit_ER_sim.b1;
            
            index_i=index_i+1;
         end
        
         
         % Fit power law to the mean
         g=fittype('a*x^b');
         ER_mean_power_fit=fit(ER_bins(ER_bins>1)',mean_ER_band(ER_bins>1)',g,'startpoint',[ 2.5 -.1]);
         ER_mean_fit{bin_counter}=ER_mean_power_fit.a.*ER_bins.^(ER_mean_power_fit.b);
            
         ER_lower_power_fit=fit(ER_bins(ER_bins>1)',lower_bound(ER_bins>1)',g,'startpoint',[ 2.5 -.1]);
         ER_lower_fit{bin_counter}=ER_lower_power_fit.a.*ER_bins.^(ER_lower_power_fit.b);

         ER_upper_power_fit=fit(ER_bins(ER_bins>1)',upper_bound(ER_bins>1)',g,'startpoint',[ 2.5 -.1]);
         ER_upper_fit{bin_counter}=ER_upper_power_fit.a.*ER_bins.^(ER_upper_power_fit.b);
         
         ER_mean_fit_points{bin_counter}=mean_ER_band;
         ER_lower_fit_points{bin_counter}=lower_bound;
         ER_upper_fit_points{bin_counter}=upper_bound;
                 
         ER_mean_fit_points_err{bin_counter}=mean_ER_band_err;
         ER_lower_fit_points_err{bin_counter}=lower_bound_err;
         ER_upper_fit_points_err{bin_counter}=upper_bound_err;

         ER_bins_tracked{bin_counter}=ER_bins;
         
         
         

             NR_bin_data=NR_bins_tracked;
             NR_mean_refit=fit(NR_bin_data(NR_bin_data>1)',NR_mean_data(NR_bin_data>1)',g,'startpoint',[ 2.5 -.1]);
         
             
         %Gaussian Assumption            
         clear gauss_mean gauss_sig gauss_amp gauss_amp_err gauss_mean_err gauss_sig_err leakage_frac_err leakage_frac NR_mean
         for i=1:length(mean_ER_band)
             
             NR_mean=NR_mean_data(i);            
             gauss_mean= mean_ER_band(i);
             gauss_sig= ER_sigma(i);
             gauss_amp=amp(i); %estimate gaussian amplitude from the individual bins
             gauss_amp_err=amp_err(i); %estimate gaussian amplitude error from the individual bins
             gauss_mean_err=mean_ER_band_err(i); %estimate mean err from the individual bins
             gauss_sig_err=ER_sigma_sigma(i); %estimate sigma error from the individual bins

         %Using solution of definite integral so that I can propogate errors with analytic approach
         total_count=sqrt(pi./2).*gauss_amp.*gauss_sig.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) + 1 );
         leakage_count=sqrt(pi./2).*gauss_amp.*gauss_sig.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) - erf((gauss_mean-NR_mean)./(sqrt(2).*gauss_sig)));
         
         %derivative w.r.t amplitude
         total_count_damp = sqrt(pi/2).*gauss_sig.*(erf(gauss_mean/(sqrt(2)*gauss_sig)) +1);
         leakage_count_damp= sqrt(pi/2).*gauss_sig.*(erf(gauss_mean/(sqrt(2)*gauss_sig)) - erf((gauss_mean - NR_mean)/(sqrt(2)*gauss_sig)));
         
         %derivative w.r.t mean
         total_count_dmean = gauss_amp.*exp(-(gauss_mean.^2)./(2.*gauss_sig.^2));
         leakage_count_dmean = sqrt(pi/2).*gauss_amp.*gauss_sig.*( ( sqrt(2/pi).*exp(-(gauss_mean.^2)./(2.*gauss_sig.^2)))./gauss_sig - sqrt(2/pi).*exp(-((gauss_mean-NR_mean).^2)./(2.*gauss_sig.^2))./gauss_sig);
         
         %derivative w.r.t sigma
         total_count_dsig=sqrt(pi/2).*gauss_amp.* (erf ( gauss_mean./(sqrt(2).*gauss_sig)) + 1) - (gauss_amp .* gauss_mean .* exp( - (gauss_mean.^2)./(2.*gauss_sig.^2)))./gauss_sig;
         leakage_count_dsig=sqrt(pi/2).*gauss_amp.*gauss_sig.*( (sqrt(2/pi).*(gauss_mean-NR_mean).*exp( (-(gauss_mean-NR_mean).^2)./(2.*gauss_sig.^2)))./(gauss_sig.^2) ...
             - (sqrt(2/pi).*(gauss_mean).*exp( (-(gauss_mean).^2)./(2.*gauss_sig.^2)))./(gauss_sig.^2)) + ...
             sqrt(pi/2).*gauss_amp.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) - erf((gauss_mean-NR_mean)./(sqrt(2).*gauss_sig)));
         
         total_count_err=sqrt( (total_count_damp.*gauss_amp_err).^2 + (total_count_dmean.*gauss_mean_err).^2 + (total_count_dsig.*gauss_sig_err).^2);
         leakage_count_err=sqrt( (leakage_count_damp.*gauss_amp_err).^2 + (leakage_count_dmean.*gauss_mean_err).^2 + (leakage_count_dsig.*gauss_sig_err).^2);
         
         leakage_frac(i)=leakage_count./total_count;
         leakage_frac_err(i)=sqrt( (leakage_count_err./total_count).^2 + (leakage_count.*total_count_err./(total_count.^2)).^2);

         end


         gaussian_leakage_frac_tracker{bin_counter}=leakage_frac;
         gaussian_leakage_frac_err_tracker{bin_counter}=leakage_frac_err;         
              

        %Counting Data Approach          
         clear gauss_mean gauss_sig gauss_amp gauss_amp_err gauss_mean_err gauss_sig_err leakage_frac_err leakage_frac NR_mean
         for i=1:length(mean_ER_band)

         %Number of events in S1 bin
         total_count=ER_bin_count(i);
         leakage_count=ER_bin_belowNR_count(i);
  
         total_count_err=sqrt(total_count);
         leakage_count_err=sqrt(leakage_count);
         
         leakage_frac(i)=leakage_count./total_count;
         leakage_frac_err(i)=sqrt( (leakage_count_err./total_count).^2 + (leakage_count.*total_count_err./(total_count.^2)).^2);

         end


         counting_leakage_frac_tracker{bin_counter}=leakage_frac;
         counting_leakage_frac_err_tracker{bin_counter}=leakage_frac_err;            
         %Based on ER Fit Approach
         i=1;
         clear gauss_mean gauss_sig gauss_amp gauss_amp_err gauss_mean_err gauss_sig_err leakage_frac_err leakage_frac NR_mean
         for S1_val=1:1:49
             
                NR_mean=NR_mean_refit.a.*S1_val.^(NR_mean_refit.b);
             
             gauss_mean= ER_mean_power_fit.a.*S1_val.^(ER_mean_power_fit.b);
             gauss_sig= abs((ER_upper_power_fit.a.*S1_val.^(ER_upper_power_fit.b)) - gauss_mean)/1.28;
             gauss_amp=interp1(ER_bins,gauss_amp_list,S1_val); %estimate gaussian amplitude from the individual bins
             gauss_amp_err=interp1(ER_bins,gauss_amp_list_err,S1_val); %estimate gaussian amplitude error from the individual bins
             gauss_mean_err=interp1(ER_bins,mean_ER_band_err,S1_val); %estimate mean err from the individual bins
             gauss_sig_err=interp1(ER_bins,gauss_sig_list_err,S1_val); %estimate sigma error from the individual bins

         %Using solution of definite integral so that I can propogate errors with analytic approach
         total_count=sqrt(pi./2).*gauss_amp.*gauss_sig.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) + 1 );
         leakage_count=sqrt(pi./2).*gauss_amp.*gauss_sig.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) - erf((gauss_mean-NR_mean)./(sqrt(2).*gauss_sig)));
         
         %derivative w.r.t amplitude
         total_count_damp = sqrt(pi/2).*gauss_sig.*(erf(gauss_mean/(sqrt(2)*gauss_sig)) +1);
         leakage_count_damp= sqrt(pi/2).*gauss_sig.*(erf(gauss_mean/(sqrt(2)*gauss_sig)) - erf((gauss_mean - NR_mean)/(sqrt(2)*gauss_sig)));
         
         %derivative w.r.t mean
         total_count_dmean = gauss_amp.*exp(-(gauss_mean.^2)./(2.*gauss_sig.^2));
         leakage_count_dmean = sqrt(pi/2).*gauss_amp.*gauss_sig.*( ( sqrt(2/pi).*exp(-(gauss_mean.^2)./(2.*gauss_sig.^2)))./gauss_sig - sqrt(2/pi).*exp(-((gauss_mean-NR_mean).^2)./(2.*gauss_sig.^2))./gauss_sig);
         
         %derivative w.r.t sigma
         total_count_dsig=sqrt(pi/2).*gauss_amp.* (erf ( gauss_mean./(sqrt(2).*gauss_sig)) + 1) - (gauss_amp .* gauss_mean .* exp( - (gauss_mean.^2)./(2.*gauss_sig.^2)))./gauss_sig;
         leakage_count_dsig=sqrt(pi/2).*gauss_amp.*gauss_sig.*( (sqrt(2/pi).*(gauss_mean-NR_mean).*exp( (-(gauss_mean-NR_mean).^2)./(2.*gauss_sig.^2)))./(gauss_sig.^2) ...
             - (sqrt(2/pi).*(gauss_mean).*exp( (-(gauss_mean).^2)./(2.*gauss_sig.^2)))./(gauss_sig.^2)) + ...
             sqrt(pi/2).*gauss_amp.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) - erf((gauss_mean-NR_mean)./(sqrt(2).*gauss_sig)));
         
         total_count_err=sqrt( (total_count_damp.*gauss_amp_err).^2 + (total_count_dmean.*gauss_mean_err).^2 + (total_count_dsig.*gauss_sig_err).^2);
         leakage_count_err=sqrt( (leakage_count_damp.*gauss_amp_err).^2 + (leakage_count_dmean.*gauss_mean_err).^2 + (leakage_count_dsig.*gauss_sig_err).^2);
         
         leakage_frac(i)=leakage_count./total_count;
         leakage_frac_err(i)=sqrt( (leakage_count_err./total_count).^2 + (leakage_count.*total_count_err./(total_count.^2)).^2);
         i=i+1;
        
         end
         
         Fit_leakage_frac_tracker{bin_counter}=leakage_frac;
         Fit_leakage_frac_err_tracker{bin_counter}=leakage_frac_err;          
        
         bin_counter=bin_counter+1;       
         
          %% Move onto next Z slice
max_r=18;
s1_min=0;
s1_max=50;


pulse_area_cut= (log10(s1_phe_both)+0.7*log10(s2_phe_both)<4.85) & (log10(s1_phe_both)+0.7*log10(s2_phe_both)>2) & ...
    (log10(s1_phe_both)-2*log10(s2_phe_both_xyz)<-4.7) & (log10(s1_phe_both)-2*log10(s2_phe_both)>-6.7);
fiducial_cut=inrange(drift_time,[180 225]) & inrange(s2radius_del,[0,max_r]);

tritium_cut= clean_cut & pulse_area_cut  & fiducial_cut;

  


%% Make Band Plot
   band_sigma=1.28; %90% confidence bounds
    
             
             
            s2_s1_bin= 1:0.05:5;
           ER_bin_size=3; %3 works
           index_i=1;
           ER_bins=ER_bin_size/2+s1_min:ER_bin_size:s1_max;
               lower_bound=ones(size(ER_bins));
               mean_ER_band=ones(size(ER_bins));
               upper_bound=ones(size(ER_bins));
                counts_ER_hist=ones(size(ER_bins));
               ER_sigma=ones(size(ER_bins));
               ER_sigma_sigma=ones(size(ER_bins));
        
%             temporary, for calculating correct E-lifetime   
            clear mean_ER_band_err lower_bound_err upper_bound_err lower_bound mean_ER_band upper_bound  ...
                gauss_sig_list_err gauss_amp_list gauss_amp_list_err
            
                     %Refit the NR band - ER binning is same as NR binning
             NR_mean_data=NR_mean_fit_points;
             NR_mean_data_err=NR_mean_fit_points_err;
                  
            
        for ER_bin = (ER_bin_size/2)+s1_min:ER_bin_size:(s1_max); %HIST uses bin centers
   
        
            ER_fit_cut=tritium_cut & inrange(s1_phe_both_xyz,[ER_bin-ER_bin_size/2,ER_bin+ER_bin_size/2]);
             
            ER_band_hist=hist(log10(s2_phe_both_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            counts_ER_hist(index_i)=sum(ER_band_hist);
           
            Fit_ER=fit(s2_s1_bin',ER_band_hist','gauss1');
             
            rms_fit=fit_ML_normal(log10(s2_phe_both_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            
            %band=confint(Fit_ER,.1);
%             ER_sigma(index_i)= Fit_ER.c1/sqrt(2);
            ER_sigma(index_i)= sqrt(rms_fit.sig2);
            ER_sigma_sigma(index_i)=ER_sigma(index_i)/sqrt(2*counts_ER_hist(index_i));
            
            lower_bound(index_i)=Fit_ER.b1-ER_sigma(index_i)*band_sigma;
            mean_ER_band(index_i)=Fit_ER.b1;
            upper_bound(index_i)=Fit_ER.b1+ER_sigma(index_i)*band_sigma;


            ER_bin_count(index_i)=length(s1_phe_both_xyz(ER_fit_cut));
            ER_bin_belowNR_count(index_i)=length(s1_phe_both_xyz(ER_fit_cut & log10(s2_phe_both_xyz./s1_phe_both_xyz)<NR_mean_data(index_i)));
            
            temp_err=confint(Fit_ER,0.68);            
            mean_ER_band_err(index_i)=abs(temp_err(1,2)-Fit_ER.b1);
            lower_bound_err(index_i)=sqrt(mean_ER_band_err(index_i).^2 + (band_sigma.*abs(temp_err(1,3)-Fit_ER.c1)/sqrt(2)).^2);
            upper_bound_err(index_i)=sqrt(mean_ER_band_err(index_i).^2 + (band_sigma.*abs(temp_err(1,3)-Fit_ER.c1)/sqrt(2)).^2);

            amp(index_i)=Fit_ER.a1;
            amp_err(index_i)=abs(temp_err(1,1)-Fit_ER.a1);
            
            gauss_sig_list_err(index_i)=abs(temp_err(1,3) - Fit_ER.c1)./sqrt(2);
            gauss_amp_list(index_i)=Fit_ER.a1;
            gauss_amp_list_err(index_i)=abs(temp_err(1,1)-Fit_ER.a1);
            %mean_ER_band_sim(index_i)=Fit_ER_sim.b1;
            
            index_i=index_i+1;
         end
        
         
         % Fit power law to the mean
         g=fittype('a*x^b');
         ER_mean_power_fit=fit(ER_bins(ER_bins>1)',mean_ER_band(ER_bins>1)',g,'startpoint',[ 2.5 -.1]);
         ER_mean_fit{bin_counter}=ER_mean_power_fit.a.*ER_bins.^(ER_mean_power_fit.b);
            
         ER_lower_power_fit=fit(ER_bins(ER_bins>1)',lower_bound(ER_bins>1)',g,'startpoint',[ 2.5 -.1]);
         ER_lower_fit{bin_counter}=ER_lower_power_fit.a.*ER_bins.^(ER_lower_power_fit.b);

         ER_upper_power_fit=fit(ER_bins(ER_bins>1)',upper_bound(ER_bins>1)',g,'startpoint',[ 2.5 -.1]);
         ER_upper_fit{bin_counter}=ER_upper_power_fit.a.*ER_bins.^(ER_upper_power_fit.b);
         
         ER_mean_fit_points{bin_counter}=mean_ER_band;
         ER_lower_fit_points{bin_counter}=lower_bound;
         ER_upper_fit_points{bin_counter}=upper_bound;
                 
         ER_mean_fit_points_err{bin_counter}=mean_ER_band_err;
         ER_lower_fit_points_err{bin_counter}=lower_bound_err;
         ER_upper_fit_points_err{bin_counter}=upper_bound_err;

         ER_bins_tracked{bin_counter}=ER_bins;
         
         
         

             NR_bin_data=NR_bins_tracked;
             NR_mean_refit=fit(NR_bin_data(NR_bin_data>1)',NR_mean_data(NR_bin_data>1)',g,'startpoint',[ 2.5 -.1]);
         
             
         %Gaussian Assumption            
         clear gauss_mean gauss_sig gauss_amp gauss_amp_err gauss_mean_err gauss_sig_err leakage_frac_err leakage_frac NR_mean
         for i=1:length(mean_ER_band)
             
             NR_mean=NR_mean_data(i);            
             gauss_mean= mean_ER_band(i);
             gauss_sig= ER_sigma(i);
             gauss_amp=amp(i); %estimate gaussian amplitude from the individual bins
             gauss_amp_err=amp_err(i); %estimate gaussian amplitude error from the individual bins
             gauss_mean_err=mean_ER_band_err(i); %estimate mean err from the individual bins
             gauss_sig_err=ER_sigma_sigma(i); %estimate sigma error from the individual bins

         %Using solution of definite integral so that I can propogate errors with analytic approach
         total_count=sqrt(pi./2).*gauss_amp.*gauss_sig.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) + 1 );
         leakage_count=sqrt(pi./2).*gauss_amp.*gauss_sig.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) - erf((gauss_mean-NR_mean)./(sqrt(2).*gauss_sig)));
         
         %derivative w.r.t amplitude
         total_count_damp = sqrt(pi/2).*gauss_sig.*(erf(gauss_mean/(sqrt(2)*gauss_sig)) +1);
         leakage_count_damp= sqrt(pi/2).*gauss_sig.*(erf(gauss_mean/(sqrt(2)*gauss_sig)) - erf((gauss_mean - NR_mean)/(sqrt(2)*gauss_sig)));
         
         %derivative w.r.t mean
         total_count_dmean = gauss_amp.*exp(-(gauss_mean.^2)./(2.*gauss_sig.^2));
         leakage_count_dmean = sqrt(pi/2).*gauss_amp.*gauss_sig.*( ( sqrt(2/pi).*exp(-(gauss_mean.^2)./(2.*gauss_sig.^2)))./gauss_sig - sqrt(2/pi).*exp(-((gauss_mean-NR_mean).^2)./(2.*gauss_sig.^2))./gauss_sig);
         
         %derivative w.r.t sigma
         total_count_dsig=sqrt(pi/2).*gauss_amp.* (erf ( gauss_mean./(sqrt(2).*gauss_sig)) + 1) - (gauss_amp .* gauss_mean .* exp( - (gauss_mean.^2)./(2.*gauss_sig.^2)))./gauss_sig;
         leakage_count_dsig=sqrt(pi/2).*gauss_amp.*gauss_sig.*( (sqrt(2/pi).*(gauss_mean-NR_mean).*exp( (-(gauss_mean-NR_mean).^2)./(2.*gauss_sig.^2)))./(gauss_sig.^2) ...
             - (sqrt(2/pi).*(gauss_mean).*exp( (-(gauss_mean).^2)./(2.*gauss_sig.^2)))./(gauss_sig.^2)) + ...
             sqrt(pi/2).*gauss_amp.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) - erf((gauss_mean-NR_mean)./(sqrt(2).*gauss_sig)));
         
         total_count_err=sqrt( (total_count_damp.*gauss_amp_err).^2 + (total_count_dmean.*gauss_mean_err).^2 + (total_count_dsig.*gauss_sig_err).^2);
         leakage_count_err=sqrt( (leakage_count_damp.*gauss_amp_err).^2 + (leakage_count_dmean.*gauss_mean_err).^2 + (leakage_count_dsig.*gauss_sig_err).^2);
         
         leakage_frac(i)=leakage_count./total_count;
         leakage_frac_err(i)=sqrt( (leakage_count_err./total_count).^2 + (leakage_count.*total_count_err./(total_count.^2)).^2);

         end


         gaussian_leakage_frac_tracker{bin_counter}=leakage_frac;
         gaussian_leakage_frac_err_tracker{bin_counter}=leakage_frac_err;         
              

        %Counting Data Approach          
         clear gauss_mean gauss_sig gauss_amp gauss_amp_err gauss_mean_err gauss_sig_err leakage_frac_err leakage_frac NR_mean
         for i=1:length(mean_ER_band)

         %Number of events in S1 bin
         total_count=ER_bin_count(i);
         leakage_count=ER_bin_belowNR_count(i);
  
         total_count_err=sqrt(total_count);
         leakage_count_err=sqrt(leakage_count);
         
         leakage_frac(i)=leakage_count./total_count;
         leakage_frac_err(i)=sqrt( (leakage_count_err./total_count).^2 + (leakage_count.*total_count_err./(total_count.^2)).^2);

         end


         counting_leakage_frac_tracker{bin_counter}=leakage_frac;
         counting_leakage_frac_err_tracker{bin_counter}=leakage_frac_err;            
         %Based on ER Fit Approach
         i=1;
         clear gauss_mean gauss_sig gauss_amp gauss_amp_err gauss_mean_err gauss_sig_err leakage_frac_err leakage_frac NR_mean
         for S1_val=1:1:49
             
                NR_mean=NR_mean_refit.a.*S1_val.^(NR_mean_refit.b);
             
             gauss_mean= ER_mean_power_fit.a.*S1_val.^(ER_mean_power_fit.b);
             gauss_sig= abs((ER_upper_power_fit.a.*S1_val.^(ER_upper_power_fit.b)) - gauss_mean)/1.28;
             gauss_amp=interp1(ER_bins,gauss_amp_list,S1_val); %estimate gaussian amplitude from the individual bins
             gauss_amp_err=interp1(ER_bins,gauss_amp_list_err,S1_val); %estimate gaussian amplitude error from the individual bins
             gauss_mean_err=interp1(ER_bins,mean_ER_band_err,S1_val); %estimate mean err from the individual bins
             gauss_sig_err=interp1(ER_bins,gauss_sig_list_err,S1_val); %estimate sigma error from the individual bins

         %Using solution of definite integral so that I can propogate errors with analytic approach
         total_count=sqrt(pi./2).*gauss_amp.*gauss_sig.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) + 1 );
         leakage_count=sqrt(pi./2).*gauss_amp.*gauss_sig.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) - erf((gauss_mean-NR_mean)./(sqrt(2).*gauss_sig)));
         
         %derivative w.r.t amplitude
         total_count_damp = sqrt(pi/2).*gauss_sig.*(erf(gauss_mean/(sqrt(2)*gauss_sig)) +1);
         leakage_count_damp= sqrt(pi/2).*gauss_sig.*(erf(gauss_mean/(sqrt(2)*gauss_sig)) - erf((gauss_mean - NR_mean)/(sqrt(2)*gauss_sig)));
         
         %derivative w.r.t mean
         total_count_dmean = gauss_amp.*exp(-(gauss_mean.^2)./(2.*gauss_sig.^2));
         leakage_count_dmean = sqrt(pi/2).*gauss_amp.*gauss_sig.*( ( sqrt(2/pi).*exp(-(gauss_mean.^2)./(2.*gauss_sig.^2)))./gauss_sig - sqrt(2/pi).*exp(-((gauss_mean-NR_mean).^2)./(2.*gauss_sig.^2))./gauss_sig);
         
         %derivative w.r.t sigma
         total_count_dsig=sqrt(pi/2).*gauss_amp.* (erf ( gauss_mean./(sqrt(2).*gauss_sig)) + 1) - (gauss_amp .* gauss_mean .* exp( - (gauss_mean.^2)./(2.*gauss_sig.^2)))./gauss_sig;
         leakage_count_dsig=sqrt(pi/2).*gauss_amp.*gauss_sig.*( (sqrt(2/pi).*(gauss_mean-NR_mean).*exp( (-(gauss_mean-NR_mean).^2)./(2.*gauss_sig.^2)))./(gauss_sig.^2) ...
             - (sqrt(2/pi).*(gauss_mean).*exp( (-(gauss_mean).^2)./(2.*gauss_sig.^2)))./(gauss_sig.^2)) + ...
             sqrt(pi/2).*gauss_amp.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) - erf((gauss_mean-NR_mean)./(sqrt(2).*gauss_sig)));
         
         total_count_err=sqrt( (total_count_damp.*gauss_amp_err).^2 + (total_count_dmean.*gauss_mean_err).^2 + (total_count_dsig.*gauss_sig_err).^2);
         leakage_count_err=sqrt( (leakage_count_damp.*gauss_amp_err).^2 + (leakage_count_dmean.*gauss_mean_err).^2 + (leakage_count_dsig.*gauss_sig_err).^2);
         
         leakage_frac(i)=leakage_count./total_count;
         leakage_frac_err(i)=sqrt( (leakage_count_err./total_count).^2 + (leakage_count.*total_count_err./(total_count.^2)).^2);
         i=i+1;
        
         end
         
         Fit_leakage_frac_tracker{bin_counter}=leakage_frac;
         Fit_leakage_frac_err_tracker{bin_counter}=leakage_frac_err; 
        
         bin_counter=bin_counter+1;        
 
                   %% Move onto next Z slice
max_r=18;
s1_min=0;
s1_max=50;


pulse_area_cut= (log10(s1_phe_both)+0.7*log10(s2_phe_both)<4.85) & (log10(s1_phe_both)+0.7*log10(s2_phe_both)>2) & ...
    (log10(s1_phe_both)-2*log10(s2_phe_both_xyz)<-4.7) & (log10(s1_phe_both)-2*log10(s2_phe_both)>-6.7);
fiducial_cut=inrange(drift_time,[235 280]) & inrange(s2radius_del,[0,max_r]);

tritium_cut= clean_cut & pulse_area_cut  & fiducial_cut;





%% Make Band Plot
   band_sigma=1.28; %90% confidence bounds
    
             
             
            s2_s1_bin= 1:0.05:5;
           ER_bin_size=3; %3 works
           index_i=1;
           ER_bins=ER_bin_size/2+s1_min:ER_bin_size:s1_max;
               lower_bound=ones(size(ER_bins));
               mean_ER_band=ones(size(ER_bins));
               upper_bound=ones(size(ER_bins));
                counts_ER_hist=ones(size(ER_bins));
               ER_sigma=ones(size(ER_bins));
               ER_sigma_sigma=ones(size(ER_bins));
        
%             temporary, for calculating correct E-lifetime   
            clear mean_ER_band_err lower_bound_err upper_bound_err lower_bound mean_ER_band upper_bound  ...
                gauss_sig_list_err gauss_amp_list gauss_amp_list_err
            
                     %Refit the NR band - ER binning is same as NR binning
             NR_mean_data=NR_mean_fit_points;
             NR_mean_data_err=NR_mean_fit_points_err;
                  
            
        for ER_bin = (ER_bin_size/2)+s1_min:ER_bin_size:(s1_max); %HIST uses bin centers
   
        
            ER_fit_cut=tritium_cut & inrange(s1_phe_both_xyz,[ER_bin-ER_bin_size/2,ER_bin+ER_bin_size/2]);
             
            ER_band_hist=hist(log10(s2_phe_both_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            counts_ER_hist(index_i)=sum(ER_band_hist);
           
            Fit_ER=fit(s2_s1_bin',ER_band_hist','gauss1');
             
            rms_fit=fit_ML_normal(log10(s2_phe_both_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            
            %band=confint(Fit_ER,.1);
%             ER_sigma(index_i)= Fit_ER.c1/sqrt(2);
            ER_sigma(index_i)= sqrt(rms_fit.sig2);
            ER_sigma_sigma(index_i)=ER_sigma(index_i)/sqrt(2*counts_ER_hist(index_i));
            
            lower_bound(index_i)=Fit_ER.b1-ER_sigma(index_i)*band_sigma;
            mean_ER_band(index_i)=Fit_ER.b1;
            upper_bound(index_i)=Fit_ER.b1+ER_sigma(index_i)*band_sigma;


            ER_bin_count(index_i)=length(s1_phe_both_xyz(ER_fit_cut));
            ER_bin_belowNR_count(index_i)=length(s1_phe_both_xyz(ER_fit_cut & log10(s2_phe_both_xyz./s1_phe_both_xyz)<NR_mean_data(index_i)));
            
            temp_err=confint(Fit_ER,0.68);            
            mean_ER_band_err(index_i)=abs(temp_err(1,2)-Fit_ER.b1);
            lower_bound_err(index_i)=sqrt(mean_ER_band_err(index_i).^2 + (band_sigma.*abs(temp_err(1,3)-Fit_ER.c1)/sqrt(2)).^2);
            upper_bound_err(index_i)=sqrt(mean_ER_band_err(index_i).^2 + (band_sigma.*abs(temp_err(1,3)-Fit_ER.c1)/sqrt(2)).^2);

            amp(index_i)=Fit_ER.a1;
            amp_err(index_i)=abs(temp_err(1,1)-Fit_ER.a1);
            
            gauss_sig_list_err(index_i)=abs(temp_err(1,3) - Fit_ER.c1)./sqrt(2);
            gauss_amp_list(index_i)=Fit_ER.a1;
            gauss_amp_list_err(index_i)=abs(temp_err(1,1)-Fit_ER.a1);
            %mean_ER_band_sim(index_i)=Fit_ER_sim.b1;
            
            index_i=index_i+1;
         end
        
         
         % Fit power law to the mean
         g=fittype('a*x^b');
         ER_mean_power_fit=fit(ER_bins(ER_bins>1)',mean_ER_band(ER_bins>1)',g,'startpoint',[ 2.5 -.1]);
         ER_mean_fit{bin_counter}=ER_mean_power_fit.a.*ER_bins.^(ER_mean_power_fit.b);
            
         ER_lower_power_fit=fit(ER_bins(ER_bins>1)',lower_bound(ER_bins>1)',g,'startpoint',[ 2.5 -.1]);
         ER_lower_fit{bin_counter}=ER_lower_power_fit.a.*ER_bins.^(ER_lower_power_fit.b);

         ER_upper_power_fit=fit(ER_bins(ER_bins>1)',upper_bound(ER_bins>1)',g,'startpoint',[ 2.5 -.1]);
         ER_upper_fit{bin_counter}=ER_upper_power_fit.a.*ER_bins.^(ER_upper_power_fit.b);
         
         ER_mean_fit_points{bin_counter}=mean_ER_band;
         ER_lower_fit_points{bin_counter}=lower_bound;
         ER_upper_fit_points{bin_counter}=upper_bound;
                 
         ER_mean_fit_points_err{bin_counter}=mean_ER_band_err;
         ER_lower_fit_points_err{bin_counter}=lower_bound_err;
         ER_upper_fit_points_err{bin_counter}=upper_bound_err;

         ER_bins_tracked{bin_counter}=ER_bins;
         
         
         

             NR_bin_data=NR_bins_tracked;
             NR_mean_refit=fit(NR_bin_data(NR_bin_data>1)',NR_mean_data(NR_bin_data>1)',g,'startpoint',[ 2.5 -.1]);
         
             
         %Gaussian Assumption            
         clear gauss_mean gauss_sig gauss_amp gauss_amp_err gauss_mean_err gauss_sig_err leakage_frac_err leakage_frac NR_mean
         for i=1:length(mean_ER_band)
             
             NR_mean=NR_mean_data(i);            
             gauss_mean= mean_ER_band(i);
             gauss_sig= ER_sigma(i);
             gauss_amp=amp(i); %estimate gaussian amplitude from the individual bins
             gauss_amp_err=amp_err(i); %estimate gaussian amplitude error from the individual bins
             gauss_mean_err=mean_ER_band_err(i); %estimate mean err from the individual bins
             gauss_sig_err=ER_sigma_sigma(i); %estimate sigma error from the individual bins

         %Using solution of definite integral so that I can propogate errors with analytic approach
         total_count=sqrt(pi./2).*gauss_amp.*gauss_sig.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) + 1 );
         leakage_count=sqrt(pi./2).*gauss_amp.*gauss_sig.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) - erf((gauss_mean-NR_mean)./(sqrt(2).*gauss_sig)));
         
         %derivative w.r.t amplitude
         total_count_damp = sqrt(pi/2).*gauss_sig.*(erf(gauss_mean/(sqrt(2)*gauss_sig)) +1);
         leakage_count_damp= sqrt(pi/2).*gauss_sig.*(erf(gauss_mean/(sqrt(2)*gauss_sig)) - erf((gauss_mean - NR_mean)/(sqrt(2)*gauss_sig)));
         
         %derivative w.r.t mean
         total_count_dmean = gauss_amp.*exp(-(gauss_mean.^2)./(2.*gauss_sig.^2));
         leakage_count_dmean = sqrt(pi/2).*gauss_amp.*gauss_sig.*( ( sqrt(2/pi).*exp(-(gauss_mean.^2)./(2.*gauss_sig.^2)))./gauss_sig - sqrt(2/pi).*exp(-((gauss_mean-NR_mean).^2)./(2.*gauss_sig.^2))./gauss_sig);
         
         %derivative w.r.t sigma
         total_count_dsig=sqrt(pi/2).*gauss_amp.* (erf ( gauss_mean./(sqrt(2).*gauss_sig)) + 1) - (gauss_amp .* gauss_mean .* exp( - (gauss_mean.^2)./(2.*gauss_sig.^2)))./gauss_sig;
         leakage_count_dsig=sqrt(pi/2).*gauss_amp.*gauss_sig.*( (sqrt(2/pi).*(gauss_mean-NR_mean).*exp( (-(gauss_mean-NR_mean).^2)./(2.*gauss_sig.^2)))./(gauss_sig.^2) ...
             - (sqrt(2/pi).*(gauss_mean).*exp( (-(gauss_mean).^2)./(2.*gauss_sig.^2)))./(gauss_sig.^2)) + ...
             sqrt(pi/2).*gauss_amp.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) - erf((gauss_mean-NR_mean)./(sqrt(2).*gauss_sig)));
         
         total_count_err=sqrt( (total_count_damp.*gauss_amp_err).^2 + (total_count_dmean.*gauss_mean_err).^2 + (total_count_dsig.*gauss_sig_err).^2);
         leakage_count_err=sqrt( (leakage_count_damp.*gauss_amp_err).^2 + (leakage_count_dmean.*gauss_mean_err).^2 + (leakage_count_dsig.*gauss_sig_err).^2);
         
         leakage_frac(i)=leakage_count./total_count;
         leakage_frac_err(i)=sqrt( (leakage_count_err./total_count).^2 + (leakage_count.*total_count_err./(total_count.^2)).^2);

         end


         gaussian_leakage_frac_tracker{bin_counter}=leakage_frac;
         gaussian_leakage_frac_err_tracker{bin_counter}=leakage_frac_err;         
              

        %Counting Data Approach          
         clear gauss_mean gauss_sig gauss_amp gauss_amp_err gauss_mean_err gauss_sig_err leakage_frac_err leakage_frac NR_mean
         for i=1:length(mean_ER_band)

         %Number of events in S1 bin
         total_count=ER_bin_count(i);
         leakage_count=ER_bin_belowNR_count(i);
  
         total_count_err=sqrt(total_count);
         leakage_count_err=sqrt(leakage_count);
         
         leakage_frac(i)=leakage_count./total_count;
         leakage_frac_err(i)=sqrt( (leakage_count_err./total_count).^2 + (leakage_count.*total_count_err./(total_count.^2)).^2);

         end


         counting_leakage_frac_tracker{bin_counter}=leakage_frac;
         counting_leakage_frac_err_tracker{bin_counter}=leakage_frac_err;            

          %Based on ER Fit Approach
         i=1;
         clear gauss_mean gauss_sig gauss_amp gauss_amp_err gauss_mean_err gauss_sig_err leakage_frac_err leakage_frac NR_mean
         for S1_val=1:1:49
             
                NR_mean=NR_mean_refit.a.*S1_val.^(NR_mean_refit.b);
             
             gauss_mean= ER_mean_power_fit.a.*S1_val.^(ER_mean_power_fit.b);
             gauss_sig= abs((ER_upper_power_fit.a.*S1_val.^(ER_upper_power_fit.b)) - gauss_mean)/1.28;
             gauss_amp=interp1(ER_bins,gauss_amp_list,S1_val); %estimate gaussian amplitude from the individual bins
             gauss_amp_err=interp1(ER_bins,gauss_amp_list_err,S1_val); %estimate gaussian amplitude error from the individual bins
             gauss_mean_err=interp1(ER_bins,mean_ER_band_err,S1_val); %estimate mean err from the individual bins
             gauss_sig_err=interp1(ER_bins,gauss_sig_list_err,S1_val); %estimate sigma error from the individual bins

         %Using solution of definite integral so that I can propogate errors with analytic approach
         total_count=sqrt(pi./2).*gauss_amp.*gauss_sig.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) + 1 );
         leakage_count=sqrt(pi./2).*gauss_amp.*gauss_sig.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) - erf((gauss_mean-NR_mean)./(sqrt(2).*gauss_sig)));
         
         %derivative w.r.t amplitude
         total_count_damp = sqrt(pi/2).*gauss_sig.*(erf(gauss_mean/(sqrt(2)*gauss_sig)) +1);
         leakage_count_damp= sqrt(pi/2).*gauss_sig.*(erf(gauss_mean/(sqrt(2)*gauss_sig)) - erf((gauss_mean - NR_mean)/(sqrt(2)*gauss_sig)));
         
         %derivative w.r.t mean
         total_count_dmean = gauss_amp.*exp(-(gauss_mean.^2)./(2.*gauss_sig.^2));
         leakage_count_dmean = sqrt(pi/2).*gauss_amp.*gauss_sig.*( ( sqrt(2/pi).*exp(-(gauss_mean.^2)./(2.*gauss_sig.^2)))./gauss_sig - sqrt(2/pi).*exp(-((gauss_mean-NR_mean).^2)./(2.*gauss_sig.^2))./gauss_sig);
         
         %derivative w.r.t sigma
         total_count_dsig=sqrt(pi/2).*gauss_amp.* (erf ( gauss_mean./(sqrt(2).*gauss_sig)) + 1) - (gauss_amp .* gauss_mean .* exp( - (gauss_mean.^2)./(2.*gauss_sig.^2)))./gauss_sig;
         leakage_count_dsig=sqrt(pi/2).*gauss_amp.*gauss_sig.*( (sqrt(2/pi).*(gauss_mean-NR_mean).*exp( (-(gauss_mean-NR_mean).^2)./(2.*gauss_sig.^2)))./(gauss_sig.^2) ...
             - (sqrt(2/pi).*(gauss_mean).*exp( (-(gauss_mean).^2)./(2.*gauss_sig.^2)))./(gauss_sig.^2)) + ...
             sqrt(pi/2).*gauss_amp.*( erf(gauss_mean./(sqrt(2).*gauss_sig)) - erf((gauss_mean-NR_mean)./(sqrt(2).*gauss_sig)));
         
         total_count_err=sqrt( (total_count_damp.*gauss_amp_err).^2 + (total_count_dmean.*gauss_mean_err).^2 + (total_count_dsig.*gauss_sig_err).^2);
         leakage_count_err=sqrt( (leakage_count_damp.*gauss_amp_err).^2 + (leakage_count_dmean.*gauss_mean_err).^2 + (leakage_count_dsig.*gauss_sig_err).^2);
         
         leakage_frac(i)=leakage_count./total_count;
         leakage_frac_err(i)=sqrt( (leakage_count_err./total_count).^2 + (leakage_count.*total_count_err./(total_count.^2)).^2);
         i=i+1;
        
         end
         
         Fit_leakage_frac_tracker{bin_counter}=leakage_frac;
         Fit_leakage_frac_err_tracker{bin_counter}=leakage_frac_err;                 

         %% Plot the final result
          colors={[1 0 0],[0 1 0],[0 0 1],[0 0 0],[0.8 0.5 0]};

 figure;
 hold on
 ylabel('Counting Leakage Fraction'); xlabel('S1 corrected (Pulse Area Phe)'); 
 myfigview(16);
%  ylim([0 3.5]);xlim([0 s1_max]);

        for i=1:length(counting_leakage_frac_tracker)
         errorbar(ER_bins_tracked{i},counting_leakage_frac_tracker{i},counting_leakage_frac_err_tracker{i},'Color',colors{i})
        end
  legend('35-75 uSec','75-125 uSec','125-175 uSec','180-225 uSec','235-280 uSec')       

 figure;
 hold on
 ylabel('Gaussian Leakage Fraction'); xlabel('S1 corrected (Pulse Area Phe)'); 
 myfigview(16);
%  ylim([0 3.5]);xlim([0 s1_max]);

        for i=1:length(counting_leakage_frac_tracker)
         errorbar(ER_bins_tracked{i},gaussian_leakage_frac_tracker{i},gaussian_leakage_frac_err_tracker{i},'Color',colors{i})
        end   

        legend('35-75 uSec','75-125 uSec','125-175 uSec','180-225 uSec','235-280 uSec')
        
 figure;
 hold on
 ylabel('Fit Leakage Fraction'); xlabel('S1 corrected (Pulse Area Phe)'); 
 myfigview(16);
%  ylim([0 3.5]);xlim([0 s1_max]);

        for i=1:length(counting_leakage_frac_tracker)
         errorbar([1:1:49],Fit_leakage_frac_tracker{i},Fit_leakage_frac_err_tracker{i},'Color',colors{i})
        end   

        legend('35-75 uSec','75-125 uSec','125-175 uSec','180-225 uSec','235-280 uSec')       