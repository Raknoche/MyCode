%The desired outputs are the "field" corrections,s2z_normalization_fit, norm_S2_xy (and s2_xy_xbins s2_xy_ybins),
%mean_fit_s1z and norm_s1_xy (and s1_xy_xbins s1_xy_ybins)

%Resrict s1a s1b to (inrange(s1a_phe_both,[90 250]) & inrange(s1b_phe_both,[35,90]))


%Apply corrections with
% s2_phe_bottom_xyz_z=s2_phe_bottom_xyz.*polyval(s2z_normalization_fit,0)./polyval(s2z_normalization_fit,drift_time)
% s2_phe_bottom_xyz_xyz=s2_phe_bottom_xyz_z.*interp2(s2_xy_xbins,s2_xy_ybins,norm_S2_xy,s2x,s2y,'spline');
% s1_phe_both_xyz_z=s1_phe_both_xyz.*polyval(mean_fit_s1z,(det_edge-4)/2)./polyval(mean_fit_s1z,drift_time);
% s1_phe_both_xyz_xyz=s1_phe_both_xyz_z.*interp2(s1_xy_xbins,s1_xy_ybins,norm_S1_both,s2x,s2y,'spline');


%%%%%%%%%%%%%%%%%
%First measure the "true" lifetime from CH3T, after removing field dependence based on NEST
%Note that earlier versions of this code assumed no field dependence in CH3T
%%%%%%%%%%%%%%%%%

%Load in the uncorrected CH3T data, using cp 170** (Kr 2p13 IQs and DP2.1)
load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p21\TritiumMethod\libNEST_CH3TModel\CH3T_Sep2015_NoCorrections.mat');

%Construct data vectors
s2radius=(s2x.^2+s2y.^2).^(1/2);
s2_phe_bottom=s2_phe_bottom(s2radius<25);
s2_phe_both=s2_phe_both(s2radius<25);
s1_phe_both=s1_phe_both(s2radius<25);
s1_phe_bottom=s1_phe_bottom(s2radius<25);
s2x=s2x(s2radius<25);
s2y=s2y(s2radius<25);
drift_time=drift_time(s2radius<25);
s2radius=s2radius(s2radius<25);

%Cut down to only CH3T data
h3_area_cut=(inrange(s2_phe_both,[630 14000]) & inrange(s1_phe_both,[1 125]));
h3_s2_phe_bottom=s2_phe_bottom(h3_area_cut);
h3_s2_phe_both=s2_phe_both(h3_area_cut);
h3_drift_time=drift_time(h3_area_cut);
h3_radius=s2radius(h3_area_cut);
h3_x=s2x(h3_area_cut);
h3_y=s2y(h3_area_cut);

%Sticking with Scott's 2D map instead of Lucie's 3D map for 4 reasons
% 1) It is easier to interpolate and extrapolate a 2D map
% 2) The maps are very similar until you get to the extrapolation regions
% 3) Lucie's map claims very strong jump in field right before the grids... and Scott's is more well behaved
% 4) Lucie's map has discrete jumps in the field over the entire lenght of z, due to how she places her charges
% %load in 3D field map from Lucie
% load('C:\Program Files\MATLAB\R2012a\bin\Run04_ERBand\Kr2p20\Sep2015_3Dfield_map.mat');
% %inpaint it with spring method for interpolation
% FieldMap2015=inpaint_nans3(FieldMap2015,1);

%Estimate field at CH3T data points
h3_field=EstimateField_Sep2015(h3_radius,h3_drift_time);

%Plot the estimate of the field
x=h3_radius(inrange(h3_radius,[0 25]) & inrange(h3_drift_time,[0 300])); y=h3_drift_time(inrange(h3_radius,[0 25]) & inrange(h3_drift_time,[0 300]));
p=h3_field(inrange(h3_radius,[0 25]) & inrange(h3_drift_time,[0 300]));
x2=round(x-min(x)+1);
y2=round(y-min(y)+1);
A=[];
for k=1:length(x),
   A(x2(k),y2(k))=p(k);
end
figure
imagesc(A.');
c=colormap(jet);
c(1,:)=[1 1 1];
colormap(c);
h=colorbar;
myfigview(16);
ylabel(h,'Field Strength (V/cm)');
xlabel('Radius (cm)');
ylabel('Drift Time (uSec)');

%Use NEST to remove the field effect in the CH3T data, assuming E=2.5 keV
Energy_estimate=2.5; %Setting to 2.5 for now, maybe use better estimate for each data set in the future... it does have a large effect on r.  Have to change r later in the code if I do so
TI_param=0.046.*h3_field.^-0.14; %based on PandaX work, which is consistent with Pixie work
N_q=Energy_estimate./0.0137; %number of quanta
alpha=0.11; %exciton to ion ratio
N_ion=N_q./(1+alpha);
h3_recomb=1-log(1+(TI_param./4).*N_ion)./((TI_param./4).*N_ion); %From NEST, based on energy (fixed for now) and where the event is in space

%Plot estimate of recombination for CH3T
clear x y p
x=h3_radius(inrange(h3_radius,[0 25]) & inrange(h3_drift_time,[0 300])); y=h3_drift_time(inrange(h3_radius,[0 25]) & inrange(h3_drift_time,[0 300]));
p=h3_recomb(inrange(h3_radius,[0 25]) & inrange(h3_drift_time,[0 300]));
x2=round(x-min(x)+1);
y2=round(y-min(y)+1);
A=[];
for k=1:length(x),
   A(x2(k),y2(k))=p(k);
end
figure
imagesc(A.');
c=colormap(jet);
c(1,:)=[1 1 1];
colormap(c);
h=colorbar;
myfigview(16);
caxis([0.25, 0.32]);
ylabel(h,'Recombination Probability for 2.5 keV ER');
xlabel('Radius (cm)');
ylabel('Drift Time (uSec)');

%%%% Detector Center Detection %%%%%%%
h3x_center=mean(h3_x(h3_drift_time>4)); %exclude extraction field region at uSec<4
h3y_center=mean(h3_y(h3_drift_time>4));
h3z_center=mean(h3_drift_time(h3_drift_time>4));

%Determine the field dependence from CH3T data
%This is doing event by event basis, maybe it'd be better to do a 3D map and interpolate instead.
Nphot=N_q.*(alpha./(1+alpha))+N_ion.*h3_recomb;
Nelec=N_q-Nphot;
Nphot_center=mean(Nphot(inrange(h3_drift_time,[h3z_center-30 h3z_center+30]) & inrange(h3_x,[h3x_center-5 h3x_center+5])  & inrange(h3_y,[h3y_center-5 h3y_center+5])));
Nelec_center=mean(Nelec(inrange(h3_drift_time,[h3z_center-30 h3z_center+30]) & inrange(h3_x,[h3x_center-5 h3x_center+5])  & inrange(h3_y,[h3y_center-5 h3y_center+5])));
h3_r_center=mean(h3_recomb(inrange(h3_drift_time,[h3z_center-30 h3z_center+30]) & inrange(h3_x,[h3x_center-5 h3x_center+5])  & inrange(h3_y,[h3y_center-5 h3y_center+5])));

h3_S1fieldremoval=Nphot_center./Nphot;
h3_S2fieldremoval=Nelec_center./Nelec;

%Removing field from S2.  Don't care about S1 since it isn't used in this code
h3_s2_phe_bottom=h3_s2_phe_bottom.*h3_S2fieldremoval;
h3_s2_phe_both=h3_s2_phe_both.*h3_S2fieldremoval;



%% Get lifetime for s2_bottom -- there's a ~4% systematic based on what I pick s2_limit to be, 2% systematic based on what I pick binning to be.  Keeping dT spacing large helps negate these effects
s2_limit=2000; %For TAGCB=1,7,1,8.5,2 s2_limit=2500, FOR TAGCB=1,7,2,10,2 s2_limit=1200
s2_min=400; %For TAGCB=1,7,1,8.5,2 = s2_min=100, FOR TAGCB=1,7,2,10,2 s2_min=s2_min

dT_step=10;
i=1;

for dT_max=10+dT_step:dT_step:330;
h3_mean(i)=histfitlandau(h3_s2_phe_bottom(inrange(h3_radius,[0 25]) & inrange(h3_drift_time,[dT_max-dT_step dT_max])), 150, 400, (s2_limit+(285-dT_max)*2.0408));
dT_loc(i)=dT_max-dT_step/2;
i=i+1;
end

h3_mean_fit=fit(dT_loc(2:end-1).',h3_mean(2:end-1).','exp1');
fit_points=0:2:310;
h3_fit_line=h3_mean_fit.a.*exp(h3_mean_fit.b.*fit_points);

figure
scatter(dT_loc,h3_mean,'.k');
hold on;
plot(fit_points,h3_fit_line,'-r','LineWidth',2);
xlabel('Drift Time (uSec)');ylabel('Uncorrected CH3T S2Bot Peak (phe)'); myfigview(16);
title('Electron Lifetime Extracted from Landau Fits to CH3T');
e_lifetime_landau=-1/h3_mean_fit.b;
b=confint(h3_mean_fit, 0.683);%1 sigma matrix
sigma_e_lifetime=(1/b(1,2)-1/b(2,2))/2;
legend('S2\_bottom Mean', strcat( '\lambda= ', num2str(e_lifetime_landau,4), ' \pm ', num2str(sigma_e_lifetime,2), ' \mus' )...
        ,'location','northeast');

    
    
s2_bot_norm_z_index=dT_loc;
s2_bot_norm_z_means=h3_mean;
    s2_at_top=RK_1DCubicInterp_LinearExtrap(s2_bot_norm_z_index,s2_bot_norm_z_means,4); %normalize to right below the gate, spline doesnt work well here
    
%Checking interp
figure
hold on;
plot([0:1:350],RK_1DCubicInterp_LinearExtrap(s2_bot_norm_z_index,s2_bot_norm_z_means,[0:1:350]),'.r') %spline doesnt work well here
plot(s2_bot_norm_z_index,s2_bot_norm_z_means,'.k')
xlabel('Drift Time');
ylabel('S2 bot interp');
myfigview(16);
    
    s2_bot_norm_z=RK_1DCubicInterp_LinearExtrap(s2_bot_norm_z_index,s2_bot_norm_z_means,4)./s2_bot_norm_z_means;
    %%Tracker for pseudo-lifetime
    s2_at_bottom=RK_1DCubicInterp_LinearExtrap(s2_bot_norm_z_index,s2_bot_norm_z_means,320);
    pseudo_lifetime_2014_bot=-320/log(s2_at_bottom/s2_at_top);
    
%Show some slices for Landau
dT_step=60;
figure
hold on;
colors={'r','b','g'};
colors2={[1 0.6 0.6],[0.6 0.6 1],[0.6 1 0.6]};
i=1;
    while i <4
for dT_max=30+dT_step:dT_step:310-dT_step;
val=h3_s2_phe_bottom(inrange(h3_drift_time,[dT_max-dT_step dT_max]));
passo=75;
inizio=s2_min;
fine=(s2_limit+(285-dT_max)*2.0408);%s2_limit;
[y,x]=hist_(val,passo,min(val),fine,fine,1);
step(val,[inizio:passo:fine],colors2{i})
x=x+passo/2;
inz=find(x<inizio);
    if isempty(inz)
        inz=1;
    end;
y(1:inz(end))=0;
a0=max(y);
mpv0=min(x(find(y==a0)));
sigma0=std(val);
ftype=fittype('a*landau(x,mpv,sigma,0)');
fitres=fit(x,y,ftype,'StartPoint',[a0 mpv0 sigma0]);
bound=confint(fitres);
xfit=x(1):(passo/10):x(end);
yfit=fitres(xfit);
temp=line(xfit,yfit,'LineWidth',2,'LineStyle','-','Color',colors{i});
hAnnotation = get(temp,'Annotation');
hLegendEntry = get(hAnnotation','LegendInformation');
set(hLegendEntry,'IconDisplayStyle','off')
picco=max(yfit);
mpv=xfit(find(yfit==picco));
mpv=mpv(end);
mv=mean(val);
temp=line([mpv mpv],[0 picco],'LineWidth',2,'LineStyle','--','Color',colors{i});
hAnnotation = get(temp,'Annotation');
hLegendEntry = get(hAnnotation','LegendInformation');
set(hLegendEntry,'IconDisplayStyle','off')
i=i+1;
    end
    end
legend('30-90 uSec','90-150 uSec','150-210 uSec');
xlabel('Uncorrected S2Bot (phe)');ylabel('Counts'); title('Landau Fits to S2Bot Spectra');
myfigview(16);

electron_lifetime_bot=e_lifetime_landau;
% h3_s2_phe_bottom_z=h3_s2_phe_bottom.*exp(h3_drift_time./electron_lifetime_bot);
 h3_s2_phe_bottom_z=h3_s2_phe_bottom.*RK_1DCubicInterp_LinearExtrap(s2_bot_norm_z_index,s2_bot_norm_z, h3_drift_time);

 %Checking interp
 figure
 plot(h3_drift_time,RK_1DCubicInterp_LinearExtrap(s2_bot_norm_z_index,s2_bot_norm_z, h3_drift_time),'.r') %again spline doesnt work well here
 hold on;
 plot(s2_bot_norm_z_index,s2_bot_norm_z,'.k')
 
figure %so landau fit doesnt override previous plot

%Showing that the correction worked
i=1;
dT_step=10;
for dT_max=10+dT_step:dT_step:330;
h3_mean_test(i)=histfitlandau(h3_s2_phe_bottom_z(inrange(h3_radius,[0 25]) & inrange(h3_drift_time,[dT_max-dT_step dT_max])), 30, 400, (s2_limit+(285-dT_max)*2.0408));
dT_loc_test(i)=dT_max-dT_step/2;
i=i+1;
end
figure
scatter(dT_loc_test,h3_mean_test,'.k');
hold on;
scatter(dT_loc_test,h3_mean,'.r');
xlabel('Drift Time (uSec)'); ylabel('CH3T S2 (phe)');
title('Showing the Lifetime Correction Worked - S2 bottom');
legend('Lifetime Corrected','Uncorrected')
myfigview(16);

% clearvars -except h3_s2_phe_bottom_z h3_radius h3_drift_time s2_limit s2_min h3_s2_phe_bottom e_lifetime_landau h3_energy h3_x h3_y

%% Get lifetime from S2_both
clear h3_mean dT_loc
dT_step=10;
i=1;
s2_limit=5000;

for dT_max=10+dT_step:dT_step:330;
h3_mean(i)=histfitlandau(h3_s2_phe_both(inrange(h3_radius,[0 25]) & inrange(h3_drift_time,[dT_max-dT_step dT_max])), 75, 400, (s2_limit+(285-dT_max)*2.5));
dT_loc(i)=dT_max-dT_step/2;
i=i+1;
end

h3_mean_fit=fit(dT_loc(2:end-1).',h3_mean(2:end-1).','exp1');
fit_points=30:2:310;
h3_fit_line=h3_mean_fit.a.*exp(h3_mean_fit.b.*fit_points);

figure
scatter(dT_loc,h3_mean,'.k');
hold on;
plot(fit_points,h3_fit_line,'-r','LineWidth',2);
xlabel('Drift Time (uSec)');ylabel('Uncorrected CH3T S2Both Peak (phe)'); myfigview(16);
title('Electron Lifetime Extracted from Landau Fits to CH3T');
e_lifetime_landau=-1/h3_mean_fit.b;
b=confint(h3_mean_fit, 0.683);%1 sigma matrix
sigma_e_lifetime=(1/b(1,2)-1/b(2,2))/2;
legend('S2\_both Mean', strcat( '\lambda= ', num2str(e_lifetime_landau,4), ' \pm ', num2str(sigma_e_lifetime,2), ' \mus' )...
        ,'location','northeast');
    
s2_both_norm_z_index=dT_loc;
s2_both_norm_z_means=h3_mean;
    s2_at_top=RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,4); %normalize to right below the gate
    
       
%Checking interp
figure
hold on;
plot([0:1:350],RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,[0:1:350]),'.r') %spline doesnt work well here
plot(s2_both_norm_z_index,s2_both_norm_z_means,'.k')
xlabel('Drift Time');
ylabel('S2 both interp');
myfigview(16);
    
    s2_both_norm_z=RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,4)./s2_both_norm_z_means;
    %%Tracker for pseudo-lifetime
    s2_at_bottom=RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,320);
    pseudo_lifetime_2014_both=-320/log(s2_at_bottom/s2_at_top);
    
%Show some slices for Landau
dT_step=60;
figure
hold on;
colors={'r','b','g'};
colors2={[1 0.6 0.6],[0.6 0.6 1],[0.6 1 0.6]};
i=1;
    while i <3
for dT_max=30+dT_step:dT_step:310-dT_step;
val=h3_s2_phe_both(inrange(h3_drift_time,[dT_max-dT_step dT_max]));
passo=30;
inizio=s2_min;
fine=(5000+(285-dT_max)*2.0408);%s2_limit;
[y,x]=hist_(val,passo,min(val),fine,fine,1);
step(val,[inizio:passo:fine],colors2{i})
x=x+passo/2;
inz=find(x<inizio);
    if isempty(inz)
        inz=1;
    end;
y(1:inz(end))=0;
a0=max(y);
mpv0=min(x(find(y==a0)));
sigma0=std(val);
ftype=fittype('a*landau(x,mpv,sigma,0)');
fitres=fit(x,y,ftype,'StartPoint',[a0 mpv0 sigma0]);
bound=confint(fitres);
xfit=x(1):(passo/10):x(end);
yfit=fitres(xfit);
temp=line(xfit,yfit,'LineWidth',2,'LineStyle','-','Color',colors{i});
hAnnotation = get(temp,'Annotation');
hLegendEntry = get(hAnnotation','LegendInformation');
set(hLegendEntry,'IconDisplayStyle','off')
picco=max(yfit);
mpv=xfit(find(yfit==picco));
mpv=mpv(end);
mv=mean(val);
temp=line([mpv mpv],[0 picco],'LineWidth',2,'LineStyle','--','Color',colors{i});
hAnnotation = get(temp,'Annotation');
hLegendEntry = get(hAnnotation','LegendInformation');
set(hLegendEntry,'IconDisplayStyle','off')
i=i+1;
    end
    end
legend('30-90 uSec','90-150 uSec','150-210 uSec');
xlabel('Uncorrected S2Both (phe)');ylabel('Counts'); title('Landau Fits to S2Both Spectra');
myfigview(16);

electron_lifetime_both=e_lifetime_landau;
 h3_s2_phe_both_z=h3_s2_phe_both.*RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z, h3_drift_time);


  %Checking interp
 figure
 plot(h3_drift_time,RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z, h3_drift_time),'.r') %again spline doesnt work well here
 hold on;
 plot(s2_both_norm_z_index,s2_both_norm_z,'.k')

 
figure %so previous isn't overwritten by landau fit

%Showing that the correction worked
dT_step=10;
% dT_max=20+dT_step;
i=1;

clear h3_mean_test dT_loc_test
for dT_max=10+dT_step:dT_step:330;
h3_mean_test(i)=histfitlandau(h3_s2_phe_both_z(inrange(h3_radius,[0 25]) & inrange(h3_drift_time,[dT_max-dT_step dT_max])), 30, 400, (5000+(285-dT_max)*2.0408));
dT_loc_test(i)=dT_max-dT_step/2;
i=i+1;
end
figure
scatter(dT_loc_test,h3_mean_test,'.k');
hold on;
scatter(dT_loc_test,h3_mean,'.r');
xlabel('Drift Time (uSec)'); ylabel('CH3T S2 (phe)');
title('Showing the Lifetime Correction Worked - S2 both');
legend('Lifetime Corrected','Uncorrected')
myfigview(16);


%%
%%%%%%%%%%%%%%%%%%%
% Get S2 XY geomtric corrections based on CH3T
%%%%%%%%%%%%%%%%%%%
clear x2 hist_bin

xx_step=3;
yy_step=3;
%set up matricies to be filled
h3_s2_bottom_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(h3_s2_phe_bottom_z(inrange(h3_drift_time,[10 330]) & inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step]))) > 600;
            h3_s2_bottom_xymeans(j,i)=histfitlandau(h3_s2_phe_bottom_z(inrange(h3_drift_time,[10 330]) & inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step])), 100, 400, 2500);
             
           else 
          
        h3_s2_bottom_xymeans(j,i)=0;     
        
      end
      j=j+1;
    end
    i=i+1;
   end
   
    h3_s2_bottom_xymeans(isnan(h3_s2_bottom_xymeans)) = 0;
     
   
%%Plot the correction S2_XY bottom %%%%%%%%%%%%%%%

s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;

old_s2xbins=s2x_bins;
old_s2ybins=s2y_bins;
    
    color_range_max=max(max(h3_s2_bottom_xymeans));
    color_range_min=min(min(h3_s2_bottom_xymeans(h3_s2_bottom_xymeans>0)));
    vc_step=(color_range_max-color_range_min)/50;
    
    vc=color_range_min:vc_step:color_range_max;

    s2_xy_bottom_fig = figure;
    contourf(s2x_bins,s2y_bins,h3_s2_bottom_xymeans,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('S2 Mean (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);

    %Remove zeros before interpolating
%     h3_s2_bottom_xymeans(h3_s2_bottom_xymeans==0)=NaN;
%     h3_s2_bottom_xymeans=inpaint_nans(h3_s2_bottom_xymeans,3);

    %Calculate corrections matrix and 1 sigma corrections matrix
center=interp2(s2x_bins,s2y_bins,h3_s2_bottom_xymeans,h3x_center,h3y_center,'spline');%Normalize to the center (x=y=0)
norm_S2_bot=center./h3_s2_bottom_xymeans; 
norm_S2_bot(isinf(norm_S2_bot))=1;%no infinity! Do not correct events outside r=25cm
norm_S2_bot(isnan(norm_S2_bot))=1;
% 
% sigma_center=interp2(s2xbins,s2ybins,Sigma_S2_bottom,0,0,'spline',1);%Normalize to the center (x=y=0)
% sigma_norm_S2_bot=sqrt((sigma_center./mean_S2_bottom_xy).^2+(Sigma_S2_bottom.*center./mean_S2_bottom_xy.^2).^2);
% sigma_norm_S2_bot(isinf(sigma_norm_S2_bot))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
% sigma_norm_S2_bot(isnan(sigma_norm_S2_bot))=1;%no nan. Set Sigma=1, which is 100% uncertainty
% 
h3_s2_phe_bottom_xyz=h3_s2_phe_bottom_z.*interp2(s2x_bins,s2y_bins,norm_S2_bot,h3_x,h3_y,'spline',1);

figure
%Showing that the correction worked
xx_step=3;
yy_step=3;
h3_s2_bottom_xymeans_test = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(h3_s2_phe_bottom_z(inrange(h3_drift_time,[10 330]) & inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step]))) > 600;
            h3_s2_bottom_xymeans_test(j,i)=histfitlandau(h3_s2_phe_bottom_xyz(inrange(h3_drift_time,[10 330]) & inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step])), 100, 400, 2500);
             
           else 
          
        h3_s2_bottom_xymeans_test(j,i)=0;     
        
      end
      j=j+1;
    end
    i=i+1;
   end
   
    h3_s2_bottom_xymeans_test(isnan(h3_s2_bottom_xymeans_test)) = 0;
     
   
%%Plot the correction S2_XY bottom %%%%%%%%%%%%%%%

s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;
    
    color_range_max=max(max(h3_s2_bottom_xymeans_test));
    color_range_min=min(min(h3_s2_bottom_xymeans_test(h3_s2_bottom_xymeans_test>0)));
    vc_step=(color_range_max-color_range_min)/50;
    
    vc=color_range_min:vc_step:color_range_max;

    s2_xy_bottom_fig = figure;
    contourf(s2x_bins,s2y_bins,h3_s2_bottom_xymeans_test,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('XYZ Corrected S2 Mean (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);
%% Getting S2 XY for using S2 both

%%%%%%%%%%%%%%%%%%%
% Get S2 XY geomtric corrections based on CH3T
%%%%%%%%%%%%%%%%%%%
clear x2 hist_bin

xx_step=3;
yy_step=3;
%set up matricies to be filled
h3_s2_both_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(h3_s2_phe_both_z(inrange(h3_drift_time,[10 330]) & inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step]))) > 600;
            h3_s2_both_xymeans(j,i)=histfitlandau(h3_s2_phe_both_z(inrange(h3_drift_time,[10 330]) & inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step])), 100, 400, 5000);
             
           else 
          
        h3_s2_both_xymeans(j,i)=0;     
        
      end
      j=j+1;
    end
    i=i+1;
   end
   
    h3_s2_both_xymeans(isnan(h3_s2_both_xymeans)) = 0;
     
   
%%Plot the correction S2_XY bottom %%%%%%%%%%%%%%%

s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;

old_s2xbins_both=s2x_bins;
old_s2ybins_both=s2y_bins;
    
    color_range_max=max(max(h3_s2_both_xymeans));
    color_range_min=min(min(h3_s2_both_xymeans(h3_s2_both_xymeans>0)));
    vc_step=(color_range_max-color_range_min)/50;
    
    vc=color_range_min:vc_step:color_range_max;

    s2_xy_both_fig = figure;
    contourf(s2x_bins,s2y_bins,h3_s2_both_xymeans,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('S2 Mean (Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);

    
        %Remove zeros before interpolating
%     h3_s2_both_xymeans(h3_s2_both_xymeans==0)=NaN;
%     h3_s2_both_xymeans=inpaint_nans(h3_s2_both_xymeans,3);
    
    %Calculate corrections matrix and 1 sigma corrections matrix
center=interp2(s2x_bins,s2y_bins,h3_s2_both_xymeans,h3x_center,h3y_center,'spline');%Normalize to the center (x=y=0)
norm_S2_both=center./h3_s2_both_xymeans; 
norm_S2_both(isinf(norm_S2_both))=1;%no infinity! Do not correct events outside r=25cm
norm_S2_both(isnan(norm_S2_both))=1;

h3_s2_phe_both_xyz=h3_s2_phe_both_z.*interp2(s2x_bins,s2y_bins,norm_S2_both,h3_x,h3_y,'spline',1);

figure
%Showing that the correction worked
xx_step=3;
yy_step=3;
h3_s2_both_xymeans_test = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(h3_s2_phe_both_z(inrange(h3_drift_time,[10 330]) & inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step]))) > 600;
            h3_s2_both_xymeans_test(j,i)=histfitlandau(h3_s2_phe_both_xyz(inrange(h3_drift_time,[10 330]) & inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step])), 100, 400, 5000);
             
           else 
          
        h3_s2_both_xymeans_test(j,i)=0;     
        
      end
      j=j+1;
    end
    i=i+1;
   end
   
    h3_s2_both_xymeans_test(isnan(h3_s2_both_xymeans_test)) = 0;
     
   
%%Plot the correction S2_XY bottom %%%%%%%%%%%%%%%

s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;
    
    color_range_max=max(max(h3_s2_both_xymeans_test));
    color_range_min=min(min(h3_s2_both_xymeans_test(h3_s2_both_xymeans_test>0)));
    vc_step=(color_range_max-color_range_min)/50;
    
    vc=color_range_min:vc_step:color_range_max;

    s2_xy_both_fig = figure;
    contourf(s2x_bins,s2y_bins,h3_s2_both_xymeans_test,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('XYZ Corrected S2 Mean (Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);




%% Loading in Krypton data

%%%%%%%%%%%%%%%%%%%%%%
%Correct the Kr data with the lifetime found above
%%%%%%%%%%%%%%%%%%%%%%
dir_path='C:\Users\Richard\Desktop\LUX Analysis\';
user_name='RichardKnoche';
dp_version='2.0';
algorithm_name='LUXkrypCal';
use_ccv_flag=1;

% folders_to_load{1}='Run04_CH3T_Kr2p13 IQs\lux10_20140903T1918_cp17059'; %tritium is on 20140905
folders_to_load{1}='lux10_20150929T1905_cp18197';

rqs_to_load = {'pulse_area_phe','event_timestamp_samples'...
   ,'pulse_classification' ...
   ,'z_drift_samples' , 's1s2_pairing'...
   ,'z_corrected_pulse_area_all_phe','xyz_corrected_pulse_area_all_phe'...
   ,'xyz_corrected_pulse_area_bot_phe'...
   ,'top_bottom_ratio','x_cm'...
   ,'y_cm','z_corrected_pulse_area_bot_phe','xyz_corrected_pulse_area_bot_phe'...
   ,'full_evt_area_phe','Kr83fit_s1a_area_phe','Kr83fit_s1b_area_phe','Kr83fit_dt_samples'};
    
path=strcat(dir_path,'/',folders_to_load{1});         
d = LUXLoadMultipleRQMs_framework(path,rqs_to_load);                
                  
%%  Defining cuts and variables

submit=0;
delay_min=0;
grid_size=2; %cm XY plane

rcut_min = 0;
rcut_max = 25;%cm
s1area_bound_min = 20;%100
s1area_bound_max = 600;%600

s2area_bound_min = 100;%200 %change to cut out CH3T if needed
s2area_bound_max = 35000;%30000

min_evts_per_bin = 200;
max_num_bins = 65;

S1xybinsize = grid_size;%cm
S1_xbin_min = -25;
S1_xbin_max = 25;
S1_ybin_min = -25;
S1_ybin_max = 25;

S2xybinsize = grid_size;%cm
S2_xbin_min = -25;
S2_xbin_max = 25;
S2_ybin_min = -25;
S2_ybin_max = 25;


% Pulse Classification
 Pulse_Class=1;
  
    d.z_drift_samples(isnan(d.z_drift_samples)) = 0.0; % get rid of NaN   
   
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


%  Old cuts
%     s1_area_cut= inrange(d.pulse_area_phe,[s1area_bound_min,s1area_bound_max]);
%     s2_area_cut= inrange(d.pulse_area_phe,[s2area_bound_min,s2area_bound_max]);
%     s1_class=(d.pulse_classification==1 | d.pulse_classification==9) & s1_area_cut ;
%     s2_class=(d.pulse_classification==2) & s2_area_cut ;  
%     s1_single_cut =logical( (s1_class & d.golden & d.selected_s1_s2).*repmat(sum(s2_class & d.golden & d.selected_s1_s2)==1,10,1));
%     s2_single_cut =logical( (s2_class & d.golden & d.selected_s1_s2).*repmat(sum(s1_class & d.golden & d.selected_s1_s2)==1,10,1));
%     
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
    event_number=d.event_number(evt_cut);
    event_timestamp_samples=d.event_timestamp_samples(evt_cut);
    
    time_wait_cut=event_timestamp_samples/1e8/60 > delay_min; %allow x min for Kr mixing

    
    s2radius = (s2x.^2+s2y.^2).^(0.5);
    clean_cut = (log10(s1_phe_both)+0.6*log10(s2_phe_both)>4.5) & s2radius<25;

         s1_phe_bottom=s1_phe_bottom(clean_cut);
         s1_phe_both=s1_phe_both(clean_cut);
         s2_phe_bottom=s2_phe_bottom(clean_cut);
         s2_phe_both=s2_phe_both(clean_cut);
         drift_time=drift_time(clean_cut);
         s2x=s2x(clean_cut);
         s2y=s2y(clean_cut);
         s2radius = (s2x.^2+s2y.^2).^(0.5);
         event_timestamp_samples=event_timestamp_samples(clean_cut);

         %%%% Detector Center Detection %%%%%%%
x_center=mean(s2x(drift_time>4)); %exclude extraction field region at uSec<4
y_center=mean(s2y(drift_time>4));
z_center=mean(drift_time(drift_time>4));
det_edge=330;
                  
%%%%%%%%%%%%%%%%%
%Correcting S2 Z and finds remaining "field" dependence
%%%%%%%%%%%%%%%%%

         %%%%%%%% CORRECTING S2 SIGNAL BASED ON CH3T %%%%%%%%%%%%%
         clear hist_s2_bottom Fit_s2_bottom means means_error mean_index mean_ratios s2_phe_bottom_zslice
%          electron_lifetime=e_lifetime_landau;

         s2_phe_bottom_z=s2_phe_bottom.*RK_1DCubicInterp_LinearExtrap(s2_bot_norm_z_index,s2_bot_norm_z, drift_time);
         s2_phe_bottom_xyz=s2_phe_bottom_z.*interp2(old_s2xbins,old_s2ybins,norm_S2_bot,s2x,s2y,'spline',1);
         s2_phe_both_z=s2_phe_both.*RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z, drift_time);
         s2_phe_both_xyz=s2_phe_both_z.*interp2(old_s2xbins_both,old_s2ybins_both,norm_S2_both,s2x,s2y,'spline',1);
 
         %Plotting the corrections
         figure
        hold on;
        plot(drift_time,RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z, drift_time).*interp2(old_s2xbins,old_s2ybins,norm_S2_both,s2x,s2y,'spline',1),'.k')
        plot(drift_time,RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z, drift_time),'.r')
        xlabel('Drift Time'); ylabel('S2 corrections'); legend('Z','XYZ');myfigview(16);
         
%Calculating remaining S2 Z dependence using S2 bot
clear dT_loc s2z_mean s2z_fit
dT_step=10;
i=1;
for dT_max=10+dT_step:dT_step:320;
s2z_fit=fit([0:100:16000]',hist(s2_phe_bottom_xyz(inrange(s2radius,[0 25]) & inrange(drift_time,[dT_max-dT_step dT_max])),[0:100:16000])','gauss1');
s2z_mean_bottom(i)=s2z_fit.b1;
znorm_dT_loc(i)=dT_max-dT_step/2;
i=i+1;
end


% % Making a plot of the S2 Z dependence
% figure
% hold on;
% i=1;
% for dT_max=90:60:210;
% s2z_fit_temp=fit([2000:200:16000]',hist(s2_phe_bottom_xyz(inrange(s2radius,[0 25]) & inrange(drift_time,[dT_max-60 dT_max])),[2000:200:16000])','gauss1');
% step(s2_phe_bottom_xyz(inrange(s2radius,[0 25]) & inrange(drift_time,[dT_max-60 dT_max])),[2000:200:16000],colors2{i});
% x=[2000:200:16000]; y=s2z_fit_temp.a1.*exp(-((x-s2z_fit_temp.b1)./s2z_fit_temp.c1).^2); 
% i=i+1;
% end
% legend('30-90 uSec','90-150 uSec','150-210 uSec');
% myfigview(16);
% xlabel('S2_E (phe)'); ylabel('Counts');
% i=1;
% for dT_max=90:60:210;
% s2z_fit_temp=fit([2000:200:16000]',hist(s2_phe_bottom_xyz(inrange(s2radius,[0 25]) & inrange(drift_time,[dT_max-60 dT_max])),[2000:200:16000])','gauss1');
% % step(s2_phe_bottom_xyz(inrange(s2radius,[0 25]) & inrange(drift_time,[dT_max-60 dT_max])),[2000:200:16000],colors{i});
% x=[2000:200:16000]; y=s2z_fit_temp.a1.*exp(-((x-s2z_fit_temp.b1)./s2z_fit_temp.c1).^2); plot(x,y,'LineWidth',2,'Color',colors{i});
% i=i+1;
% end

s2_field_dep_z_bins_bottom=znorm_dT_loc;
s2_field_dep_z_bottom=s2z_mean_bottom;

% s2z_normalization_fit_bot=polyfit(znorm_dT_loc,s2z_mean_bottom,6);
fit_points=0:2:330;
s2z_fit_line=RK_1DCubicInterp_LinearExtrap(s2_field_dep_z_bins_bottom,s2_field_dep_z_bottom,fit_points);

figure
scatter(znorm_dT_loc,s2z_mean_bottom,'.k');
hold on;
plot(fit_points,s2z_fit_line,'-r','LineWidth',2);
xlabel('Drift Time (uSec)');ylabel('CH3T Corrected Kr S2 Peak (phe)'); myfigview(16);
title('Field Effect from CH3T Corrected Kr S2 Peaks');


s2_phe_bottom_xyz_z=s2_phe_bottom_xyz.*RK_1DCubicInterp_LinearExtrap(s2_field_dep_z_bins_bottom,s2_field_dep_z_bottom,z_center)./RK_1DCubicInterp_LinearExtrap(s2_field_dep_z_bins_bottom,s2_field_dep_z_bottom,drift_time);

% %Checking interp
% figure
% plot(drift_time,s2_phe_bottom_xyz_z,'.k')

%%Showing Z correction worked
i=1;
for dT_max=10+dT_step:dT_step:320;
s2z_fit=fit([0:100:16000]',hist(s2_phe_bottom_xyz_z(inrange(s2radius,[0 25]) & inrange(drift_time,[dT_max-dT_step dT_max])),[0:100:16000])','gauss1');
s2z_mean_bottom_z(i)=s2z_fit.b1;
i=i+1;
end

figure
scatter(znorm_dT_loc,s2z_mean_bottom_z,'.k')
hold on;
scatter(znorm_dT_loc,s2z_mean_bottom,'.r');
xlabel('Drift Time (uSec)');ylabel('Kr S2 (phe)');
legend('CH3T and Field Corrected - S2 bot','CH3T Corrected');
title('Showing Field Correction works');
myfigview(16);

%%
   
%Calculating remaining S2 Z dependence using S2 both
clear dT_loc s2z_fit
dT_step=10;
i=1;
for dT_max=10+dT_step:dT_step:320;
s2z_fit_both=fit([0:200:40000]',hist(s2_phe_both_xyz(inrange(s2radius,[0 25]) & inrange(drift_time,[dT_max-dT_step dT_max])),[0:200:40000])','gauss1');
s2z_mean_both(i)=s2z_fit_both.b1;
znorm_dT_loc(i)=dT_max-dT_step/2;
i=i+1;
end

s2_field_dep_z_bins_both=znorm_dT_loc;
s2_field_dep_z_both=s2z_mean_both;

% s2z_normalization_fit_both=polyfit(znorm_dT_loc,s2z_mean_both,6);
fit_points=0:2:330;
s2z_fit_line=RK_1DCubicInterp_LinearExtrap(s2_field_dep_z_bins_both,s2_field_dep_z_both,fit_points);

figure
scatter(znorm_dT_loc,s2z_mean_both,'.k');
hold on;
plot(fit_points,s2z_fit_line,'-r','LineWidth',2);
xlabel('Drift Time (uSec)');ylabel('CH3T Corrected Kr S2Both Peak (phe)'); myfigview(16);
title('Field Effect from CH3T Corrected Kr S2Both Peaks');

% Making a plot of the S2 Z dependence
% figure
% hold on;
% i=1;
% for dT_max=90:60:210;
% s2z_fit_temp=fit([3000:200:35000]',hist(s2_phe_both_xyz(inrange(s2radius,[0 25]) & inrange(drift_time,[dT_max-60 dT_max])),[3000:200:35000])','gauss1');
% step(s2_phe_both_xyz(inrange(s2radius,[0 25]) & inrange(drift_time,[dT_max-60 dT_max])),[3000:200:35000],colors2{i});
% % x=[2000:200:16000]; y=s2z_fit_temp.a1.*exp(-((x-s2z_fit_temp.b1)./s2z_fit_temp.c1).^2); 
% i=i+1;
% end
% legend('30-90 uSec','90-150 uSec','150-210 uSec');
% myfigview(16);
% xlabel('S2_E (phe)'); ylabel('Counts');
% i=1;
% for dT_max=90:60:210;
% s2z_fit_temp=fit([3000:200:35000]',hist(s2_phe_both_xyz(inrange(s2radius,[0 25]) & inrange(drift_time,[dT_max-60 dT_max])),[3000:200:35000])','gauss1');
% % step(s2_phe_bottom_xyz(inrange(s2radius,[0 25]) & inrange(drift_time,[dT_max-60 dT_max])),[2000:200:16000],colors{i});
% x=[3000:200:35000]; y=s2z_fit_temp.a1.*exp(-((x-s2z_fit_temp.b1)./s2z_fit_temp.c1).^2); plot(x,y,'LineWidth',2,'Color',colors{i});
% i=i+1;
% end



s2_phe_both_xyz_z=s2_phe_both_xyz.*RK_1DCubicInterp_LinearExtrap(s2_field_dep_z_bins_both,s2_field_dep_z_both,z_center)./RK_1DCubicInterp_LinearExtrap(s2_field_dep_z_bins_both,s2_field_dep_z_both,drift_time);

%%Showing Z correction worked
i=1;
for dT_max=10+dT_step:dT_step:320;
s2z_fit=fit([0:200:40000]',hist(s2_phe_both_xyz_z(inrange(s2radius,[0 25]) & inrange(drift_time,[dT_max-dT_step dT_max])),[0:200:40000])','gauss1');
s2z_mean_z(i)=s2z_fit.b1;
i=i+1;
end

figure
scatter(znorm_dT_loc,s2z_mean_z,'.k')
hold on;
scatter(znorm_dT_loc,s2z_mean_both,'.r');
xlabel('Drift Time (uSec)');ylabel('Kr S2 (phe)');
legend('CH3T and Field Corrected - S2 both','CH3T Corrected');
title('Showing Field Correction works - S2 both');
myfigview(16);

%%
%Corrects S2 XY and finds remaining "field" dependence

clear s2_phe_bottom_xyz_z_xymeans s2_phe_bottom_xyz_z_xymeans_err
%Calculating remaining S2 XY dependence
xx_step=2;
yy_step=2;
s2_phe_bottom_xyz_z_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));
s2_phe_bottom_xyz_z_xymeans_err = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(s2_phe_bottom_xyz_z(inrange(drift_time,[10 330]) & inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step]))) > 300;
            s2_phe_bottom_xyz_z_xymeans_fit=fit([0:100:40000]',hist(s2_phe_bottom_xyz_z(inrange(drift_time,[10 330]) & inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step])),[0:100:40000])','gauss1');
            s2_phe_bottom_xyz_z_xymeans(j,i)=s2_phe_bottom_xyz_z_xymeans_fit.b1;
            clear temp
            temp=confint(s2_phe_bottom_xyz_z_xymeans_fit,0.683);
            s2_phe_bottom_xyz_z_xymeans_err(j,i)=abs(s2_phe_bottom_xyz_z_xymeans(j,i)-temp(1,2));
      else 
            s2_phe_bottom_xyz_z_xymeans(j,i)=0;  
            s2_phe_bottom_xyz_z_xymeans_err(j,i)=0;
      end
      j=j+1;
    end
    i=i+1;
   end
   
   
       %Remove zeros before interpolating
%     s2_phe_bottom_xyz_z_xymeans(s2_phe_bottom_xyz_z_xymeans==0)=NaN;
%     s2_phe_bottom_xyz_z_xymeans=inpaint_nans(s2_phe_bottom_xyz_z_xymeans,3);
%     s2_phe_bottom_xyz_z_xymeans_err(s2_phe_bottom_xyz_z_xymeans_err==0)=NaN;
%     s2_phe_bottom_xyz_z_xymeans_err=inpaint_nans(s2_phe_bottom_xyz_z_xymeans_err,3);    
   
s2xbins=-25+xx_step/2:xx_step:25-xx_step/2;
s2ybins=-25+yy_step/2:yy_step:25-yy_step/2;
xy_center=interp2(s2xbins,s2ybins,s2_phe_bottom_xyz_z_xymeans,x_center,y_center,'spline');%Normalize to the center (x=y=0)
norm_s2_xy_bot=xy_center./s2_phe_bottom_xyz_z_xymeans; 
norm_s2_xy_bot(isinf(norm_s2_xy_bot))=1;%no infinity! Do not correct events outside r=25cm
norm_s2_xy_bot(isnan(norm_s2_xy_bot))=1;

s2_field_dep_xy_bottom=s2_phe_bottom_xyz_z_xymeans;
s2_field_dep_xy_bottom_err=s2_phe_bottom_xyz_z_xymeans_err;
s2_field_dep_xy_xbins_bottom=s2xbins;
s2_field_dep_xy_ybins_bottom=s2ybins;

s2_phe_bottom_xyz_xyz=s2_phe_bottom_xyz_z.*interp2(s2xbins,s2ybins,norm_s2_xy_bot,s2x,s2y,'spline',1);
s2_xy_xbins_bot=s2xbins;
s2_xy_ybins_bot=s2ybins;
   
%Showing that XY correction worked
xx_step=2;
yy_step=2;
s2_phe_bottom_xyz_xyz_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(s2_phe_bottom_xyz_z(inrange(drift_time,[10 330]) & inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step]))) > 300;
            s2_phe_bottom_xyz_xyz_xymeans_fit=fit([0:100:40000]',hist(s2_phe_bottom_xyz_xyz(inrange(drift_time,[10 330]) & inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step])),[0:100:40000])','gauss1');
            s2_phe_bottom_xyz_xyz_xymeans(j,i)=s2_phe_bottom_xyz_xyz_xymeans_fit.b1;
      else 
            s2_phe_bottom_xyz_xyz_xymeans(j,i)=0;      
      end
      j=j+1;
    end
    i=i+1;
   end
s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;
    
    color_range_max=max(max(s2_phe_bottom_xyz_xyz_xymeans));
    color_range_min=min(min(s2_phe_bottom_xyz_xyz_xymeans(s2_phe_bottom_xyz_xyz_xymeans>0)));
    vc_step=(color_range_max-color_range_min)/50;
    
    vc=color_range_min:vc_step:color_range_max;

    s2_xy_bottom_fig = figure;
    contourf(s2x_bins,s2y_bins,s2_phe_bottom_xyz_xyz_xymeans,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('S2 Mean (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);
%%
%Corrects S2 Both XY and finds remaining "field" dependence

clear s2_phe_both_xyz_z_xymeans
%Calculating remaining S2 XY dependence
xx_step=2;
yy_step=2;
s2_phe_both_xyz_z_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));
s2_phe_both_xyz_z_xymeans_err = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(s2_phe_both_xyz_z(inrange(drift_time,[10 330]) & inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step]))) > 300;
            s2_phe_both_xyz_z_xymeans_fit=fit([0:200:60000]',hist(s2_phe_both_xyz_z(inrange(drift_time,[10 330]) & inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step])),[0:200:60000])','gauss1');
            s2_phe_both_xyz_z_xymeans(j,i)=s2_phe_both_xyz_z_xymeans_fit.b1;
            clear temp
            temp=confint(s2_phe_both_xyz_z_xymeans_fit,0.683);
            s2_phe_both_xyz_z_xymeans_err(j,i)=abs(s2_phe_both_xyz_z_xymeans(j,i)-temp(1,2));
      else 
            s2_phe_both_xyz_z_xymeans(j,i)=0;
            s2_phe_both_xyz_z_xymeans_err(j,i)=0;
      end
      j=j+1;
    end
    i=i+1;
   end
   
      
       %Remove zeros before interpolating
%     s2_phe_both_xyz_z_xymeans(s2_phe_both_xyz_z_xymeans==0)=NaN;
%     s2_phe_both_xyz_z_xymeans=inpaint_nans(s2_phe_both_xyz_z_xymeans,3);
%     s2_phe_both_xyz_z_xymeans_err(s2_phe_both_xyz_z_xymeans_err==0)=NaN;
%     s2_phe_both_xyz_z_xymeans_err=inpaint_nans(s2_phe_both_xyz_z_xymeans_err,3);    
   
   
s2xbins=-25+xx_step/2:xx_step:25-xx_step/2;
s2ybins=-25+yy_step/2:yy_step:25-yy_step/2;
xy_center=interp2(s2xbins,s2ybins,s2_phe_both_xyz_z_xymeans,x_center,y_center,'spline');%Normalize to the center (x=y=0)
norm_S2_xy_both=xy_center./s2_phe_both_xyz_z_xymeans; 
norm_S2_xy_both(isinf(norm_S2_xy_both))=1;%no infinity! Do not correct events outside r=25cm
norm_S2_xy_both(isnan(norm_S2_xy_both))=1;


s2_field_dep_xy_both=s2_phe_both_xyz_z_xymeans;
s2_field_dep_xy_both_err=s2_phe_both_xyz_z_xymeans_err;
s2_field_dep_xy_xbins_both=s2xbins;
s2_field_dep_xy_ybins_both=s2ybins;

s2_phe_both_xyz_xyz=s2_phe_both_xyz_z.*interp2(s2xbins,s2ybins,norm_S2_xy_both,s2x,s2y,'spline',1);
s2_xy_xbins_both=s2xbins;
s2_xy_ybins_both=s2ybins;
   
%Showing that XY correction worked
xx_step=2;
yy_step=2;
s2_phe_both_xyz_xyz_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(s2_phe_both_xyz_z(inrange(drift_time,[10 330]) & inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step]))) > 300;
            s2_phe_both_xyz_xyz_xymeans_fit=fit([0:200:60000]',hist(s2_phe_both_xyz_xyz(inrange(drift_time,[10 330]) & inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step])),[0:200:60000])','gauss1');
            s2_phe_both_xyz_xyz_xymeans(j,i)=s2_phe_both_xyz_xyz_xymeans_fit.b1;
      else 
            s2_phe_both_xyz_xyz_xymeans(j,i)=0;      
      end
      j=j+1;
    end
    i=i+1;
   end
s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;
    
    color_range_max=max(max(s2_phe_both_xyz_xyz_xymeans));
    color_range_min=min(min(s2_phe_both_xyz_xyz_xymeans(s2_phe_both_xyz_xyz_xymeans>0)));
    vc_step=(color_range_max-color_range_min)/50;
    
    vc=color_range_min:vc_step:color_range_max;

    s2_xy_both_fig = figure;
    contourf(s2x_bins,s2y_bins,s2_phe_both_xyz_xyz_xymeans,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('S2 Mean (both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);
    

%Removing zeros from the matrix
% s2_field_dep_xy_both(s2_field_dep_xy_both==0)=NaN;
% s2_field_dep_xy_both=inpaint_nans(s2_field_dep_xy_both,3);
% s2_field_dep_xy_both_err(s2_field_dep_xy_both_err==0)=NaN;
% s2_field_dep_xy_both_err=inpaint_nans(s2_field_dep_xy_both_err,3);
% s2_field_dep_xy_bottom(s2_field_dep_xy_bottom==0)=NaN;
% s2_field_dep_xy_bottom=inpaint_nans(s2_field_dep_xy_bottom,3);
% s2_field_dep_xy_bottom_err(s2_field_dep_xy_bottom_err==0)=NaN;
% s2_field_dep_xy_bottom_err=inpaint_nans(s2_field_dep_xy_bottom_err,3);

 %% Finding 3D field dependence   

det_edge=330;
z_step=30;
s2zbins=floor(0.05*det_edge):z_step:floor(0.05*det_edge)+z_step*15;
r_max=25;
xy_step=3.5;

s2xbins=-r_max+xy_step/2:xy_step:r_max-xy_step/2;
s2ybins=-r_max+xy_step/2:xy_step:r_max-xy_step/2;

Mean_s2_3D_both=zeros(length(s2ybins),length(s2xbins),length(s2zbins));
Mean_s2_3D_bottom=zeros(length(s2ybins),length(s2xbins),length(s2zbins));
Mean_s2_3D_both_err=zeros(length(s2ybins),length(s2xbins),length(s2zbins));
Mean_s2_3D_bottom_err=zeros(length(s2ybins),length(s2xbins),length(s2zbins));
Sigma_s2_3D_both=zeros(length(s2ybins),length(s2xbins),length(s2zbins));
Sigma_s2_3D_bottom=zeros(length(s2ybins),length(s2xbins),length(s2zbins));

for k = s2zbins; 
    for i = (-r_max+xy_step):xy_step:r_max %map to rows going down (y_cm)
        for j = (-r_max+xy_step):xy_step:r_max %map columns across (x_cm)
                       
            l=int8(i/xy_step+(r_max)/xy_step); %real y from 1:25
            m=int8(j/xy_step+(r_max)/xy_step); % real x from 1:25
            n=(k-floor(0.05*det_edge))/z_step + 1; %z from 1:16
            
           %sort, and make the cut. using the variable q
           spatial_cut = s2x<j & s2x>(j-xy_step) & s2y<i & s2y>(i-xy_step) & inrange(drift_time,[(k-z_step/2) , (k+z_step/2)]); %no 0th element!
                    %swap x,y due to matrix x-row, y-column definition                     
                                  
            if (length(s2_phe_both_xyz(spatial_cut)) >= 100) % at least 200 counts before fitting. 

                fit_both= fit([0:500:60000].',hist(s2_phe_both_xyz(spatial_cut),[0:500:60000]).','gauss1');
                fit_bottom= fit([0:500:60000].',hist(s2_phe_bottom_xyz(spatial_cut),[0:500:60000]).','gauss1');
                
                %Peak location (mean)
                Mean_s2_3D_both(l,m,n)=fit_both.b1;
                Mean_s2_3D_bottom(l,m,n)=fit_bottom.b1;
                
                %1-sigma of peak position
                clear temp_both temp_bottom
                temp_both=confint(fit_both,0.683);
                temp_bottom=confint(fit_both,0.683);
                Mean_s2_3D_both_err(l,m,n)=abs(temp_both(1,2)-Mean_s2_3D_both(l,m,n));
                Mean_s2_3D_bottom_err(l,m,n)=abs(temp_bottom(1,2)-Mean_s2_3D_both(l,m,n));
                Sigma_s2_3D_both(l,m,n)=fit_both.c1/sqrt(2);
                Sigma_s2_3D_bottom(l,m,n)=fit_bottom.c1/sqrt(2);   

             else %not enough stats to do the fit
 
                Mean_s2_3D_both(l,m,n)=0;
                Mean_s2_3D_both_err(l,m,n)=0;
                Mean_s2_3D_bottom(l,m,n)=0;
                Mean_s2_3D_bottom_err(l,m,n)=0;     
                Sigma_s2_3D_both(l,m,n)=0;
                Sigma_s2_3D_bottom(l,m,n)=0;                

            end
              

        end
    end
    k
end

s2_3D_xbins=s2xbins;
s2_3D_ybins=s2ybins;
s2_3D_zbins=s2zbins;
s2_field_dep_3D_both=Mean_s2_3D_both;
s2_field_dep_3D_both_err=Mean_s2_3D_both_err;
s2_field_dep_3D_bottom=Mean_s2_3D_bottom;
s2_3D_field_dep_bot_err=Mean_s2_3D_bottom_err;

%Removing zeros from the matrix
% s2_field_dep_3D_both(s2_field_dep_3D_both==0)=NaN;
% s2_field_dep_3D_both=inpaint_nans3(s2_field_dep_3D_both,0);
% s2_field_dep_3D_both_err(s2_field_dep_3D_both_err==0)=NaN;
% s2_field_dep_3D_both_err=inpaint_nans3(s2_field_dep_3D_both_err,0);
% s2_field_dep_3D_bottom(s2_field_dep_3D_bottom==0)=NaN;
% s2_field_dep_3D_bottom=inpaint_nans3(s2_field_dep_3D_bottom,0);
% s2_3D_field_dep_bot_err(s2_3D_field_dep_bot_err==0)=NaN;
% s2_3D_field_dep_bot_err=inpaint_nans3(s2_3D_field_dep_bot_err,0);

%% Convert S2 measurements into S1 measurements
% remember to norm to center instead of top in KrypCal code
% S2 Z Kr/CH3T:  s2_field_dep_z_bins_bottom and s2_field_dep_z_bottom
% S2 XY Kr/CH3T: s2_field_dep_xy_both, s2_field_dep_xy_xbins_both s2_field_dep_xy_ybins_both
% S2 3D Kr/CH3T: s2_3D_xbins,s2_3D_ybins,s2_3D_zbins,s2_field_dep_3D_both (and bottom)

%Calculate field dependence of S1 by Converting S2 field dependence into S1 field dependence
%assume r=0.7 and a=0.2
r_h3=h3_r_center;
a=0.11;
% S1 Z - bins are same as s2_field_dep_z_bins_both (or bottom)
% [s2_fieldpoly_z_both s2_fieldpoly_z_both_err]=polyfit(s2_field_dep_z_bins_both,s2_field_dep_z_both,6);
% [s2_fieldpoly_z_bottom s2_fieldpoly_z_bottom_err]=polyfit(s2_field_dep_z_bins_bottom,s2_field_dep_z_bottom,6);

s2p_s2_both_z=RK_1DCubicInterp_LinearExtrap(s2_field_dep_z_bins_both,s2_field_dep_z_both,[0:1:330])./RK_1DCubicInterp_LinearExtrap(s2_field_dep_z_bins_both,s2_field_dep_z_both,4); %Dividing by this removes the s2 field dependence, normalizing to z=0 (liquid surface); 
s2p_s2_both_z_centernorm=RK_1DCubicInterp_LinearExtrap(s2_field_dep_z_bins_both,s2_field_dep_z_both,[0:1:330])./RK_1DCubicInterp_LinearExtrap(s2_field_dep_z_bins_both,s2_field_dep_z_both,z_center); %Dividing by this removes the s2 field dependence, normalizing to z=center; 
NionsH3=164.4;
NionsKr=2732;
S2_kr_c=RK_1DCubicInterp_LinearExtrap(s2_field_dep_z_bins_both,s2_field_dep_z_both,z_center);
S2_h3_c=histfitlandau(h3_s2_phe_both_xyz(inrange(h3_drift_time,[h3z_center-30 h3z_center+30]) & inrange(h3_x,[h3x_center-5 h3x_center+5])  & inrange(h3_y,[h3y_center-5 h3y_center+5])), 100, 400, 6000); 
r_c=1-(S2_kr_c/S2_h3_c)*(NionsH3/NionsKr)*(1-r_h3);
s1p_s1_both_z=(a+1-s2p_s2_both_z_centernorm.*(1-r_c))/(a+r_c);


%Old incorrect way: [a+1-(s2p_s2_both_z_centernorm).*(1-r)]/(a+r); %Dividing by this removes the s1 field dependence, normalizing to z=center; 

s2p_s2_bot_z=RK_1DCubicInterp_LinearExtrap(s2_field_dep_z_bins_bottom,s2_field_dep_z_bottom,[0:1:330])./RK_1DCubicInterp_LinearExtrap(s2_field_dep_z_bins_bottom,s2_field_dep_z_bottom,4); 
s2p_s2_bot_z_centernorm=RK_1DCubicInterp_LinearExtrap(s2_field_dep_z_bins_bottom,s2_field_dep_z_bottom,[0:1:330])./RK_1DCubicInterp_LinearExtrap(s2_field_dep_z_bins_bottom,s2_field_dep_z_bottom,z_center); 
%Old incorrect way: s1p_s1_bot_z=[a+1-(s2p_s2_bot_z_centernorm).*(1-r)]/(a+r);
S2_kr_c_bot=RK_1DCubicInterp_LinearExtrap(s2_field_dep_z_bins_bottom,s2_field_dep_z_bottom,z_center);
S2_h3_c_bot=histfitlandau(h3_s2_phe_bottom_xyz(inrange(h3_drift_time,[h3z_center-30 h3z_center+30]) & inrange(h3_x,[h3x_center-5 h3x_center+5])  & inrange(h3_y,[h3y_center-5 h3y_center+5])), 100, 400, 6000); %Just average corrects CH3T over detector, since we claim there is no field dependence
r_c=1-(S2_kr_c_bot/S2_h3_c_bot)*(NionsH3/NionsKr)*(1-r_h3);
s1p_s1_bot_z=(a+1-s2p_s2_bot_z_centernorm.*(1-r_c))/(a+r_c);

%Plot the field measurements
figure
hold on;
plot([0:1:330],s2p_s2_both_z_centernorm,'-b','LineWidth',2)
plot([0:1:330],s2p_s2_bot_z_centernorm,'-g','LineWidth',2)
plot([0:1:330],s1p_s1_both_z,'-r','LineWidth',2)
plot([0:1:330],s1p_s1_bot_z,'-m','LineWidth',2)
xlabel('Drift time (uSec)'); ylabel('Field effect at Z/Field effect at center');
legend('S2 Both','S2 Bot','S1 Both','S1 Bot');
myfigview(16);

%CONVERT THIS TO INTERP TO THE CENTER
% S1 XY - bins are the same as s2_field_dep_xy_xbins_both and s2_field_dep_xy_ybins_both
s2p_s2_both_xy=interp2(s2_field_dep_xy_xbins_both,s2_field_dep_xy_ybins_both,s2_field_dep_xy_both,[-25:1:25].',[-25:1:25])./interp2(s2_field_dep_xy_xbins_both,s2_field_dep_xy_ybins_both,s2_field_dep_xy_both,x_center,y_center); %Field normalized to center, divide by this to remove
S2_kr_c=interp2(s2_field_dep_xy_xbins_both,s2_field_dep_xy_ybins_both,s2_field_dep_xy_both,x_center,y_center);
S2_h3_c=histfitlandau(h3_s2_phe_both_xyz(inrange(h3_drift_time,[h3z_center-30 h3z_center+30]) & inrange(h3_x,[h3x_center-5 h3x_center+5])  & inrange(h3_y,[h3y_center-5 h3y_center+5])), 100, 400, 5000); %Just average corrects CH3T over detector, since we claim there is no field dependence
r_c=1-(S2_kr_c/S2_h3_c)*(NionsH3/NionsKr)*(1-r_h3);
s1p_s1_both_xy=(a+1-s2p_s2_both_xy.*(1-r_c))/(a+r_c);
%Old incorrect way: s1p_s1_both_xy=[a+1-(s2p_s2_both_xy).*(1-r)]/(a+r);

s2p_s2_bot_xy=interp2(s2_field_dep_xy_xbins_bottom,s2_field_dep_xy_ybins_bottom,s2_field_dep_xy_bottom,[-25:1:25].',[-25:1:25])./interp2(s2_field_dep_xy_xbins_bottom,s2_field_dep_xy_ybins_bottom,s2_field_dep_xy_bottom,x_center,y_center); 
S2_kr_c=interp2(s2_field_dep_xy_xbins_bottom,s2_field_dep_xy_ybins_bottom,s2_field_dep_xy_bottom,x_center,y_center);
S2_h3_c=histfitlandau(h3_s2_phe_bottom_xyz(inrange(h3_drift_time,[h3z_center-30 h3z_center+30]) & inrange(h3_x,[h3x_center-5 h3x_center+5])  & inrange(h3_y,[h3y_center-5 h3y_center+5])), 100, 400, 5000); %Just average corrects CH3T over detector, since we claim there is no field dependence
r_c=1-(S2_kr_c/S2_h3_c)*(NionsH3/NionsKr)*(1-r_h3);
s1p_s1_bot_xy=(a+1-s2p_s2_bot_xy.*(1-r_c))/(a+r_c);
%Old incorrect way:s1p_s1_bot_xy=[a+1-(s2p_s2_bot_xy).*(1-r)]/(a+r);

%Plot XY Field
    color_range_max=max(max(s2p_s2_both_xy));
    color_range_min=min(min(s2p_s2_both_xy(s2p_s2_both_xy>0.9)));
    vc_step=(color_range_max-color_range_min)/50;    
    vc=color_range_min:vc_step:color_range_max;
    figure;
    contourf([-25:1:25],[-25:1:25],s2p_s2_both_xy,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('S2 Field Effect (Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Field Effect at XY/Field Effect at Center','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);

    %repeat for bottom    
    color_range_max=max(max(s2p_s2_bot_xy));
    color_range_min=min(min(s2p_s2_bot_xy(s2p_s2_bot_xy>0.95)));
    vc_step=(color_range_max-color_range_min)/50;    
    vc=color_range_min:vc_step:color_range_max;
    figure;
    contourf([-25:1:25],[-25:1:25],s2p_s2_bot_xy,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('S2 Field Effect (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Field Effect at XY/Field Effect at Center','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);
    
%Plot S1 XY Field
    color_range_max=1.03;
    color_range_min=0.97;
    vc_step=(color_range_max-color_range_min)/50;    
    vc=color_range_min:vc_step:color_range_max;
    figure;
    contourf([-25:1:25],[-25:1:25],s1p_s1_both_xy,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('S1 Field Effect (Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    c=colormap;
    c(end,1)=1;
    c(end,2)=1;
    c(end,3)=1;
    colormap(c);
    ylabel(g,'Field Effect at XY/Field Effect at Center','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);    
    
    %Repeat for bottom
    color_range_max=1.03;
    color_range_min=0.96;
    vc_step=(color_range_max-color_range_min)/50;    
    vc=color_range_min:vc_step:color_range_max;
    figure;
    contourf([-25:1:25],[-25:1:25],s1p_s1_bot_xy,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('S1 Field Effect (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    clear c
    c=colormap;
    c(end,1)=1;
    c(end,2)=1;
    c(end,3)=1;
    colormap(c);
    ylabel(g,'Field Effect at XY/Field Effect at Center','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);    
    
% S1 3D - bins are the same as s2_3D_xbins, s2_3D_ybins, and s2_3D_zbins
s2p_s2_both_3D=interp3(s2_3D_xbins,s2_3D_ybins,s2_3D_zbins,s2_field_dep_3D_both,[-25:1:25].',[-25:1:25],[0:1:330])./interp3(s2_3D_xbins,s2_3D_ybins,s2_3D_zbins,s2_field_dep_3D_both,x_center,y_center,z_center);
S2_kr_c=interp3(s2_3D_xbins,s2_3D_ybins,s2_3D_zbins,s2_field_dep_3D_both,x_center,y_center,z_center);
S2_h3_c=histfitlandau(h3_s2_phe_both_xyz(inrange(h3_drift_time,[h3z_center-30 h3z_center+30]) & inrange(h3_x,[h3x_center-5 h3x_center+5])  & inrange(h3_y,[h3y_center-5 h3y_center+5])), 100, 400, 5000); %Just average corrects CH3T over detector, since we claim there is no field dependence
r_c=1-(S2_kr_c/S2_h3_c)*(NionsH3/NionsKr)*(1-r_h3);
s1p_s1_both_3D=(a+1-s2p_s2_both_3D.*(1-r_c))/(a+r_c);
%incorrect old way: s1p_s1_both_3D=[a+1-(s2p_s2_both_xy).*(1-r)]/(a+r);

s2p_s2_bot_3D=interp3(s2_3D_xbins,s2_3D_ybins,s2_3D_zbins,s2_field_dep_3D_bottom,[-25:1:25].',[-25:1:25],[0:1:330])./interp3(s2_3D_xbins,s2_3D_ybins,s2_3D_zbins,s2_field_dep_3D_bottom,x_center,y_center,z_center);
S2_kr_c=interp3(s2_3D_xbins,s2_3D_ybins,s2_3D_zbins,s2_field_dep_3D_bottom,x_center,y_center,z_center);
S2_h3_c=histfitlandau(h3_s2_phe_bottom_xyz(inrange(h3_drift_time,[h3z_center-30 h3z_center+30]) & inrange(h3_x,[h3x_center-5 h3x_center+5])  & inrange(h3_y,[h3y_center-5 h3y_center+5])), 100, 400, 5000); %Just average corrects CH3T over detector, since we claim there is no field dependence
r_c=1-(S2_kr_c/S2_h3_c)*(NionsH3/NionsKr)*(1-r_h3);
s1p_s1_bot_3D=(a+1-s2p_s2_bot_3D.*(1-r_c))/(a+r_c);

%incorrect old way: s1p_s1_bot_3D=[a+1-(s2p_s2_bot_3D).*(1-r)]/(a+r);


%% Compare to Kr S1a/S1b maps

% Get s1a s1b info

    s1ab_cut=logical(sum(s1_single_cut(:,:),1));
    s1a_phe_both=d.Kr83fit_s1a_area_phe(s1ab_cut);
    s1b_phe_both=d.Kr83fit_s1b_area_phe(s1ab_cut);
    s1ab_timing=d.Kr83fit_dt_samples(s1ab_cut);
    s1ab_x=d.x_cm(s2_single_cut);
    s1ab_y=d.y_cm(s2_single_cut);
    s1ab_z=d.z_drift_samples(s2_single_cut)/100;
    s1ab_radius=sqrt(s1ab_x.^2+s1ab_y.^2);
  
    s1ab_timing_cut=s1ab_timing>13 & s1b_phe_both>30 & s1ab_radius.'<=25;
       %plotting x=s1a v y=s1a/s1b cuts: y<0.035x-1, y>-0.035.*x+5; y>0.01.*x; %y<4.2; x<255
%     s1ab_areacut= (s1a_phe_both./s1b_phe_both-0.035.*s1a_phe_both<-1) & (s1a_phe_both./s1b_phe_both+0.035.*s1a_phe_both>5)...
%         & (s1a_phe_both./s1b_phe_both > 0.01.* s1a_phe_both) & (s1a_phe_both./s1b_phe_both<4.2) & (s1a_phe_both<255);
    
    s1ab_x=s1ab_x(s1ab_timing_cut);
    s1ab_y=s1ab_y(s1ab_timing_cut);
    s1ab_z=s1ab_z(s1ab_timing_cut);
    s1a_phe_both=s1a_phe_both(s1ab_timing_cut);
    s1b_phe_both=s1b_phe_both(s1ab_timing_cut);
    
    dT_step=det_edge/(length(s1ab_z)/1000); %2000 events per z bin
    
    %% making plots
    %Z Dependence of S1a
%     figure
%     hold on;
% i=1;
%     for z_max=90:60:210;
%         s1a_fit=fit([0:2:400].',hist(s1a_phe_both(s1a_phe_both>0 & inrange(s1ab_z,[z_max-60,z_max]).'),[0:2:400]).','gauss1');
%         step(s1a_phe_both(s1a_phe_both>0 & inrange(s1ab_z,[z_max-60,z_max]).'),[0:2:400],colors2{i});
%         i=i+1;
%     end
%   xlabel('S1a (phe)'); ylabel('Counts'); myfigview(16);
%   legend('30-90 uSec','90-150 uSec','150-210 uSec');
%   i=1;
%     for z_max=90:60:210;
%         s1a_fit=fit([0:2:400].',hist(s1a_phe_both(s1a_phe_both>0 & inrange(s1ab_z,[z_max-60,z_max]).'),[0:2:400]).','gauss1');
%         x=[0:2:400]; y=s1a_fit.a1.*exp(-((x-s1a_fit.b1)./s1a_fit.c1).^2); plot(x,y,'LineWidth',2,'Color',colors{i});
%         i=i+1;
%     end
    
 
    
    clear s1a_z_means s1a_z_means_err s1b_z_means s1b_z_means_err s1ab_bincenters s1ab_z_means s1ab_z_means_err
%Z Dependence of S1a/S1b
i=1;
    for z_max=10+dT_step:dT_step:det_edge;
        s1a_fit=fit([0:2:400].',hist(s1a_phe_both(s1a_phe_both>0 & inrange(s1ab_z,[z_max-dT_step,z_max]).'),[0:2:400]).','gauss1');
        s1b_fit=fit([0:2:300].',hist(s1b_phe_both(s1b_phe_both>0 & inrange(s1ab_z,[z_max-dT_step,z_max]).'),[0:2:300]).','gauss1');
        s1a_z_means(i)=s1a_fit.b1;
        s1a_z_means_err(i)=s1a_fit.c1/sqrt(2)/sqrt(length(s1a_phe_both(s1a_phe_both>0 & inrange(s1ab_z,[z_max-dT_step,z_max]).')));
        s1b_z_means(i)=s1b_fit.b1;
        s1b_z_means_err(i)=s1b_fit.c1/sqrt(2)/sqrt(length(s1b_phe_both(s1b_phe_both>0 & inrange(s1ab_z,[z_max-dT_step,z_max]).')));
        s1ab_bincenters(i)=z_max-dT_step/2;
        i=i+1;
    end

s1ab_z_means=s1a_z_means./s1b_z_means;
s1ab_z_means_err=sqrt( (s1b_z_means_err.*s1a_z_means./(s1b_z_means.^2)).^2 + (s1a_z_means_err./s1b_z_means).^2);
    
%fit polynomial to s1az and s1bz, to remove z dependence later in xy fit
[s1a_P, s1a_S]=polyfit(s1ab_bincenters,s1a_z_means,3);
[s1b_P, s1b_S]=polyfit(s1ab_bincenters,s1b_z_means,3);
[s1ab_P, s1ab_S]=polyfit(s1ab_bincenters,s1ab_z_means,2);

%S1a z dependence plot    
s1az_plot=figure;
errorbar(s1ab_bincenters,s1a_z_means,s1a_z_means_err,'.k')   
xlabel('Drift Time (uSec)');ylabel('S1a Mean (phe)'); title('S1a Z Dependence'); myfigview(16);
hold on;
plot([50:1:320],polyval(s1a_P,[50:1:320]),'-r','LineWidth',2)

%S1b z dependence plot    
s1bz_plot=figure;
errorbar(s1ab_bincenters,s1b_z_means,s1b_z_means_err,'.k')   
xlabel('Drift Time (uSec)');ylabel('S1b Mean (phe)'); title('S1b Z Dependence'); myfigview(16);
hold on;
plot([20:1:320],polyval(s1b_P,[20:1:320]),'-r','LineWidth',2)
    
%S1a/b z dependence plot
s1a_over_bz_plot=figure;
errorbar(s1ab_bincenters,s1ab_z_means,s1ab_z_means_err,'.k')   
xlabel('Drift Time (uSec)');ylabel('S1a/b'); title('S1 a/b Z Dependence'); myfigview(16);
hold on;
plot([20:1:320],polyval(s1ab_P,[20:1:320]),'-r','LineWidth',2)

%Saving variables for plot at the end

s1ab_bincenters_0307=s1ab_bincenters;
s1a_z_means_0307=s1a_z_means;
s1a_z_means_err_0307=s1a_z_means_err;
s1b_z_means_0307=s1b_z_means;
s1b_z_means_err_0307=s1b_z_means_err;
s1ab_z_means_0307=s1ab_z_means;
s1ab_z_means_err_0307=s1ab_z_means_err;
s1a_P_0307=s1a_P;
s1b_P_0307=s1b_P;
s1ab_P_0307=s1ab_P;

%% S1a/b XY map (after z correction)    

%Correcting z dependence to get XY dependence in same manner that Kr/CH3T
%map does
    s1a_phe_both_z=s1a_phe_both.'.*polyval(s1a_P,z_center)./polyval(s1a_P,s1ab_z);
    s1b_phe_both_z=s1b_phe_both.'.*polyval(s1b_P,z_center)./polyval(s1b_P,s1ab_z);
    
s1ab_xbin_min=-25;
s1ab_ybin_min=-25;
s1ab_xbin_max=25;
s1ab_ybin_max=25;
s1ab_xybinsize=sqrt((50*50)/(length(s1ab_z)/500));
clear s1a_xy_mean s1a_xy_mean_err s1b_xy_mean s1b_xy_mean_err

s1ab_xbins=s1ab_xbin_min+s1ab_xybinsize/2:s1ab_xybinsize:s1ab_xbin_max;
s1ab_ybins=s1ab_ybin_min+s1ab_xybinsize/2:s1ab_xybinsize:s1ab_ybin_max; 

s1a_xy_mean=zeros(length(s1ab_ybins),length(s1ab_xbins));
s1b_xy_mean=zeros(length(s1ab_ybins),length(s1ab_xbins));

 for x_bin=s1ab_xbin_min:s1ab_xybinsize:(s1ab_xbin_max-s1ab_xybinsize/2);
        for y_bin=s1ab_ybin_min:s1ab_xybinsize:(s1ab_ybin_max-s1ab_xybinsize/2);          
          x_min = x_bin; x_max = x_bin+s1ab_xybinsize;
          y_min = y_bin; y_max = y_bin+s1ab_xybinsize;      

            x_count=int32(1+x_bin/s1ab_xybinsize+s1ab_xbin_max/s1ab_xybinsize);
            y_count=int32(1+y_bin/s1ab_xybinsize+s1ab_ybin_max/s1ab_xybinsize); 

          bin_cut = inrange(s1ab_z,[10 330]).' & s1a_phe_both>0 & s1b_phe_both>0 & inrange(s1ab_x,[x_min,x_max]).' & inrange(s1ab_y,[y_min,y_max]).'; 
            if length(s1a_phe_both_z(bin_cut))>200;
                s1a_z_fit=fit([0:2:400].',hist(s1a_phe_both_z(bin_cut),[0:2:400]).','gauss1');
                s1a_xy_mean(y_count,x_count)=s1a_z_fit.b1;
                s1a_xy_mean_err(y_count,x_count)=s1a_z_fit.c1/sqrt(2)/length(s1a_phe_both_z(bin_cut));

                
                s1b_z_fit=fit([0:2:300].',hist(s1b_phe_both_z(bin_cut),[0:2:300]).','gauss1');
                s1b_xy_mean(y_count,x_count)=s1b_z_fit.b1;
                s1b_xy_mean_err(y_count,x_count)=s1b_z_fit.c1/sqrt(2)/length(s1b_phe_both_z(bin_cut));
            else
                s1a_xy_mean(y_count,x_count)=0;
                s1a_xy_mean_err(y_count,x_count)=0;
                s1b_xy_mean(y_count,x_count)=0;
                s1b_xy_mean_err(y_count,x_count)=0;   
            end
                
        end
 end
     s1ab_xy_mean=s1a_xy_mean./s1b_xy_mean;
     s1ab_xy_mean_err=sqrt( (s1b_xy_mean_err.*s1a_xy_mean./(s1b_xy_mean.^2)).^2 + (s1a_xy_mean_err./s1b_xy_mean).^2);
   

    %Plot s1a XY means
    color_range_max=max(max(s1a_xy_mean));
    color_range_min=min(min(s1a_xy_mean(s1a_xy_mean>0)));
    vc_step=(color_range_max-color_range_min)/50;
    vc=color_range_min:vc_step:color_range_max;    
    s1a_xy_fig = figure;
    contourf(s1ab_xbins,s1ab_ybins,s1a_xy_mean,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('0307 S1a Mean vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'S1a Mean (phe)','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);
    
    %Plot s1b XY means
    color_range_max=max(max(s1b_xy_mean));
    color_range_min=min(min(s1b_xy_mean(s1b_xy_mean>0)));
    vc_step=(color_range_max-color_range_min)/50;
    vc=color_range_min:vc_step:color_range_max;    
    s1b_xy_fig = figure;
    contourf(s1ab_xbins,s1ab_ybins,s1b_xy_mean,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('0307 S1b Mean vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'S1b Mean (phe)','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);

    %s1a/b XY Ratio
    color_range_max=max(max(s1ab_xy_mean));
    color_range_min=min(min(s1ab_xy_mean(s1ab_xy_mean>0)));
    vc_step=(color_range_max-color_range_min)/50;
    vc=color_range_min:vc_step:color_range_max;   
    s1ab_xy_fig = figure;
    contourf(s1ab_xbins,s1ab_ybins,s1ab_xy_mean,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('0307 S1a/b Ratio vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'S1a/b Ratio','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);

    %Saving variables    
    s1a_xy_mean_0307=s1a_xy_mean;
    s1a_xy_mean_err_0307=s1a_xy_mean_err;
    s1b_xy_mean_0307=s1b_xy_mean;
    s1b_xy_mean_err_0307=s1b_xy_mean_err;
    s1ab_xy_mean_0307=s1ab_xy_mean;
    s1ab_xy_mean_err_0307=s1ab_xy_mean_err;
    s1ab_xbins_0307=s1ab_xbins;
    s1ab_ybins_0307=s1ab_ybins;
    
    
%% 3D Kr s1a/s1b map

s1ab_xyz_numbins=10; %gives ~100 events per bin
s1ab_xyz_zstep=det_edge/s1ab_xyz_numbins;
s1ab_xyz_xstep=50/s1ab_xyz_numbins;
s1ab_xyz_ystep=50/s1ab_xyz_numbins;
s1ab_xyz_xmax=25;
r_max=25;

s1ab_xyz_zbins=50:s1ab_xyz_zstep:50+s1ab_xyz_zstep*s1ab_xyz_numbins; %20 us
s1ab_xyz_xbins=(-r_max+s1ab_xyz_xstep/2):s1ab_xyz_xstep:r_max;
s1ab_xyz_ybins=(-r_max+s1ab_xyz_ystep/2):s1ab_xyz_ystep:r_max;

s1a_xyz_mean=zeros(length(s1ab_xyz_ybins),length(s1ab_xyz_xbins),length(s1ab_xyz_zbins));
s1a_xyz_mean_err=zeros(length(s1ab_xyz_ybins),length(s1ab_xyz_xbins),length(s1ab_xyz_zbins));
s1b_xyz_mean=zeros(length(s1ab_xyz_ybins),length(s1ab_xyz_xbins),length(s1ab_xyz_zbins));
s1b_xyz_mean_err=zeros(length(s1ab_xyz_ybins),length(s1ab_xyz_xbins),length(s1ab_xyz_zbins));

for k = s1ab_xyz_zbins; %s1zbins are the center of the bin    
    for i = (-r_max+s1ab_xyz_ystep):s1ab_xyz_ystep:r_max %map to rows going down (y_cm) -- s1x bins are the max edge of the bin
        for j = (-r_max+s1ab_xyz_xstep):s1ab_xyz_xstep:r_max %map columns across (x_cm) -- s1y bins are the max edge of the bin
                       
            l=int8(i/s1ab_xyz_ystep+s1ab_xyz_xmax/s1ab_xyz_ystep);
            m=int8(j/s1ab_xyz_xstep+s1ab_xyz_xmax/s1ab_xyz_xstep);
            n=int8((k-50)/s1ab_xyz_zstep + 1);
            
           %sort, and make the cut. using the variable q
           q = s1ab_x<j & s1ab_x>(j-s1ab_xyz_xstep) & s1ab_y<i & s1ab_y>(i-s1ab_xyz_ystep) & inrange(s1ab_z,[(k-s1ab_xyz_zstep/2) , (k+s1ab_xyz_zstep/2)]); %no 0th element!
                    %swap x,y due to matrix x-row, y-column definition

            %Count the number of events per bin
            Count_S1ab_3D(l,m,n)=length(s1a_phe_both(q));
            
            if (Count_S1ab_3D(l,m,n) >= 60) % at least 100 counts before fitting. 
                s1a_z_fit=fit([0:2:400].',hist(s1a_phe_both(q.' & s1a_phe_both>0),[0:2:400]).','gauss1');
                s1a_xyz_mean(l,m,n)=s1a_z_fit.b1;
                s1a_xyz_mean_err(l,m,n)=s1a_z_fit.c1/sqrt(2)/length(s1a_phe_both(q.' & s1a_phe_both>0));

                
                s1b_z_fit=fit([0:2:300].',hist(s1b_phe_both(q.' & s1b_phe_both>0 ),[0:2:300]).','gauss1');
                s1b_xyz_mean(l,m,n)=s1b_z_fit.b1;
                s1b_xyz_mean_err(l,m,n)=s1b_z_fit.c1/sqrt(2)/length(s1b_phe_both(q.' & s1b_phe_both>0));
             else %not enough stats to do the fit
                s1a_xyz_mean(l,m,n)=0;
                s1a_xyz_mean_err(l,m,n)=0;
                s1b_xyz_mean(l,m,n)=0;
                s1b_xyz_mean_err(l,m,n)=0;
            end
                                      
        end
    end
k
end

     s1ab_xyz_mean=s1a_xyz_mean./s1b_xyz_mean;
     s1ab_xyz_mean_err=sqrt( (s1b_xyz_mean_err.*s1a_xyz_mean./(s1b_xyz_mean.^2)).^2 + (s1a_xyz_mean_err./s1b_xyz_mean).^2);


%Saving Variables

s1a_xyz_mean_0307=s1a_xyz_mean;
s1a_xyz_mean_err_0307=s1a_xyz_mean_err;
s1b_xyz_mean_0307=s1b_xyz_mean;
s1b_xyz_mean_err_0307=s1b_xyz_mean_err;
s1ab_xyz_mean_0307=s1ab_xyz_mean;
s1ab_xyz_mean_err_0307=s1ab_xyz_mean_err;
s1ab_xyz_zbins_0307=s1ab_xyz_zbins;
s1ab_xyz_xbins_0307=s1ab_xyz_xbins;
s1ab_xyz_ybins_0307=s1ab_xyz_ybins;


%% Compare S1aS1b to Field Measurement -- Z Dep

r=h3_r_center;
a=0.11;
% S1 Z - bins are same as s2_field_dep_z_bins_both (or bottom)
% [s2_fieldpoly_z_both s2_fieldpoly_z_both_err]=polyfit(s2_field_dep_z_bins_both,s2_field_dep_z_both,6);
% [s2_fieldpoly_z_bottom s2_fieldpoly_z_bottom_err]=polyfit(s2_field_dep_z_bins_bottom,s2_field_dep_z_bottom,6);


% s2p_s2_z_s1as1b_bins=polyval(s2_fieldpoly_z_both,s1ab_bincenters_0307)./polyval(s2_fieldpoly_z_both,z_center);
s2p_s2_z_s1as1b_bins=RK_1DCubicInterp_LinearExtrap(s2_field_dep_z_bins_both,s2_field_dep_z_both,s1ab_bincenters_0307)./RK_1DCubicInterp_LinearExtrap(s2_field_dep_z_bins_both,s2_field_dep_z_both,z_center);

S2_kr_c=RK_1DCubicInterp_LinearExtrap(s2_field_dep_z_bins_both,s2_field_dep_z_both,z_center);
[S2_h3_c S2_h3_c_err]=histfitlandau(h3_s2_phe_both_xyz(inrange(h3_drift_time,[h3z_center-30 h3z_center+30]) & inrange(h3_x,[h3x_center-5 h3x_center+5])  & inrange(h3_y,[h3y_center-5 h3y_center+5])), 100, 400, 6000); %Just average corrects CH3T over detector, since we claim there is no field dependence
S2_h3_c_err=10/S2_h3_c; %replace error measurement based on systematic errors involved with changing the bin_max value
r_c=1-(S2_kr_c/S2_h3_c)*(NionsH3/NionsKr)*(1-r_h3);
s1p_s1_z_s1as1b_bins=(a+1-s2p_s2_z_s1as1b_bins.*(1-r_c))/(a+r_c);

%temp variables for propagating errors.
% [junkone temp_err_one]=polyval(s2_fieldpoly_z_both,s1ab_bincenters_0307,s2_fieldpoly_z_both_err); %numerator
% [junktwo temp_err_two]=polyval(s2_fieldpoly_z_both,z_center,s2_fieldpoly_z_both_err); %denominator 
% r_c_err=sqrt( (S2_kr_c_err*(NionsH3/NionsKr)*(1-r_h3)/S2_h3_c).^2 + (S2_h3_c_err*(S2_kr_c/S2_h3_c.^2)*(NionsH3/NionsKr)*(1-r_h3)).^2 );
% s2p_s2_z_s1as1b_bins_err=sqrt( (temp_err_one./polyval(s2_fieldpoly_z_both,z_center)).^2 + (temp_err_two.*polyval(s2_fieldpoly_z_both,s1ab_bincenters_0307)./(polyval(s2_fieldpoly_z_both,z_center).^2)).^2 );
% s1p_s1_z_s1as1b_bins_err=sqrt( ((s2p_s2_z_s1as1b_bins_err.*(1-r_c))/(a+r_c)).^2 + (r_c_err.*s2p_s2_z_s1as1b_bins.*(1+a)/(a+r_c).^2).^2 ); %doesn't include errors on a or r

% %S2 comparison
% figure
% rkploterr(s1ab_z_means,s2p_s2_z_s1as1b_bins,s1ab_z_means_err,[],'.k');
% xlabel('S1a/S1b'); ylabel('CH3T Corrected Kr S2 at Z / CH3T Corrected Kr S2 at Center');
% xlim([2.5 3]);ylim([0.7 1.2]);
% myfigview(16);
% 
% %S1 comparison
% figure
% ploterr(s1ab_z_means,s1p_s1_z_s1as1b_bins,s1ab_z_means_err,[],'.k');
% xlabel('S1a/S1b'); ylabel('CH3T Corrected Kr S1 at Z / CH3T Corrected Kr S1 at Center');
% xlim([2.5 3]);ylim([0.8 1.2]);
% myfigview(16);

%Both at same time
figure
rkploterr(s1ab_z_means,s2p_s2_z_s1as1b_bins,s1ab_z_means_err,[],[1 0 0],'.',100,1);
hold on;
rkploterr(s1ab_z_means,s1p_s1_z_s1as1b_bins,s1ab_z_means_err,[],[0 0 0],'.',100,1);
xlabel('S1a/S1b'); ylabel('CH3T Corrected Kr at Z / CH3T Corrected Kr at Center');
xlim([2.2 3.1]);ylim([0.6 1.5]);
legend('S2 Mapping','S1 Mapping');
% legend('S2','S1'); Ploterr is screwing up the legend colors
myfigview(16);

%% Compare S1aS1b to Field Measurement -- XY Dep

s2p_s2_xy_s1as1b_bins=interp2(s2_field_dep_xy_xbins_both,s2_field_dep_xy_ybins_both,s2_field_dep_xy_both,s1ab_xbins_0307.',s1ab_ybins_0307)./interp2(s2_field_dep_xy_xbins_both,s2_field_dep_xy_ybins_both,s2_field_dep_xy_both,x_center,y_center); 
temp_err_1=interp2(s2_field_dep_xy_xbins_both,s2_field_dep_xy_ybins_both,s2_field_dep_xy_both_err,s1ab_xbins_0307.',s1ab_ybins_0307);
temp_err_2=interp2(s2_field_dep_xy_xbins_both,s2_field_dep_xy_ybins_both,s2_phe_both_xyz_z_xymeans_err,x_center,y_center);
s2p_s2_xy_s1as1b_bins_err=sqrt( (temp_err_1./interp2(s2_field_dep_xy_xbins_both,s2_field_dep_xy_ybins_both,s2_field_dep_xy_both,x_center,y_center)).^2 + (temp_err_2.*interp2(s2_field_dep_xy_xbins_both,s2_field_dep_xy_ybins_both,s2_field_dep_xy_both,s1ab_xbins_0307.',s1ab_ybins_0307)./(interp2(s2_field_dep_xy_xbins_both,s2_field_dep_xy_ybins_both,s2_field_dep_xy_both,x_center,y_center).^2)).^2 );

S2_kr_c=interp2(s2_field_dep_xy_xbins_both,s2_field_dep_xy_ybins_both,s2_field_dep_xy_both,x_center,y_center); 
S2_kr_c_err=interp2(s2_field_dep_xy_xbins_both,s2_field_dep_xy_ybins_both,s2_phe_both_xyz_z_xymeans_err,x_center,y_center);
[S2_h3_c S2_h3_c_err]=histfitlandau(h3_s2_phe_both_xyz(inrange(h3_drift_time,[h3z_center-30 h3z_center+30]) & inrange(h3_x,[h3x_center-5 h3x_center+5])  & inrange(h3_y,[h3y_center-5 h3y_center+5])), 100, 400, 6000); %Just average corrects CH3T over detector, since we claim there is no field dependence
S2_h3_c_err=10/S2_h3_c; %replace error measurement based on systematic errors involved with changing the bin_max value
r_c=1-(S2_kr_c/S2_h3_c)*(NionsH3/NionsKr)*(1-r_h3);
s1p_s1_xy_s1as1b_bins=(a+1-s2p_s2_xy_s1as1b_bins.*(1-r_c))/(a+r_c);
r_c_err=sqrt( (S2_kr_c_err*(NionsH3/NionsKr)*(1-r_h3)/S2_h3_c).^2 + (S2_h3_c_err*(S2_kr_c/S2_h3_c.^2)*(NionsH3/NionsKr)*(1-r_h3)).^2 );
s1p_s1_xy_s1as1b_bins_err=sqrt( ((s2p_s2_xy_s1as1b_bins_err.*(1-r_c))/(a+r_c)).^2 + (r_c_err.*s2p_s2_xy_s1as1b_bins.*(1+a)/(a+r_c).^2).^2 );

%S2 comparison
figure
hold on;
rkploterr(s1ab_xy_mean,s2p_s2_xy_s1as1b_bins,s1ab_xy_mean_err,[],[0 0 0],'.',100,1);
xlabel('S1a/S1b'); ylabel('CH3T Corrected Kr S2 at XY / CH3T Corrected Kr S2 at Center');
xlim([2.5 2.85]);ylim([0.96 1.09]);
myfigview(16);

%S1 comparison
figure
hold on;
rkploterr(s1ab_xy_mean,s1p_s1_xy_s1as1b_bins,s1ab_xy_mean_err,[],[0 0 0],'.',100,1);
xlabel('S1a/S1b'); ylabel('CH3T Corrected Kr S1 at XY / CH3T Corrected Kr S2 at Center');
xlim([2.5 2.85]);ylim([0.96 1.03]);
myfigview(16);




%% Compare S1aS1b to Field Measurement -- 3D Dep
s2_field_dep_3D_both_old=s2_field_dep_3D_both;
s2_field_dep_3D_both(s2_field_dep_3D_both==0)=nan;

%Can't have zeros when interpolating the map. Note: S2 and S1 XYZ maps with zeros are fine, since you interpolate the norm map, where the zeros get turned into NaN. 
%Here zeros are bad for buisness, since they aren't ignored like NaN are
temp_map=s2_field_dep_3D_both;
q= ~isnan(temp_map);
r=nearestpoint(find(~q),find(q));
s2_field_dep_3D_both=temp_map; temp_mapQ=temp_map(q); s2_field_dep_3D_both(~q)=temp_mapQ(r);s2_field_dep_3D_both=reshape(s2_field_dep_3D_both,size(s2_field_dep_3D_both_old,1),size(s2_field_dep_3D_both_old,2),size(s2_field_dep_3D_both_old,3));

s2p_s2_3D_s1as1b_bins=interp3(s2_3D_xbins,s2_3D_ybins,s2_3D_zbins,s2_field_dep_3D_both,s1ab_xyz_xbins.',s1ab_xyz_ybins,s1ab_xyz_zbins,'cubic')./interp3(s2_3D_xbins,s2_3D_ybins,s2_3D_zbins,s2_field_dep_3D_both,x_center,y_center,z_center,'cubic');
clear temp_err_1 temp_err_2
temp_err_1=interp3(s2_3D_xbins,s2_3D_ybins,s2_3D_zbins,s2_field_dep_3D_both_err,s1ab_xyz_xbins.',s1ab_xyz_ybins,s1ab_xyz_zbins,'cubic');
temp_err_2=interp3(s2_3D_xbins,s2_3D_ybins,s2_3D_zbins,s2_field_dep_3D_both_err,x_center,y_center,z_center,'cubic');
s2p_s2_3D_s1as1b_bins_err=sqrt( (temp_err_1./interp3(s2_3D_xbins,s2_3D_ybins,s2_3D_zbins,s2_field_dep_3D_both,x_center,y_center,z_center)).^2 + (temp_err_2.*interp3(s2_3D_xbins,s2_3D_ybins,s2_3D_zbins,s2_field_dep_3D_both,s1ab_xyz_xbins.',s1ab_xyz_ybins,s1ab_xyz_zbins)./(interp3(s2_3D_xbins,s2_3D_ybins,s2_3D_zbins,s2_field_dep_3D_both,x_center,y_center,z_center).^2)).^2 );

S2_kr_c=interp3(s2_3D_xbins,s2_3D_ybins,s2_3D_zbins,s2_field_dep_3D_both,x_center,y_center,z_center,'cubic');
S2_kr_c_err=interp3(s2_3D_xbins,s2_3D_ybins,s2_3D_zbins,s2_field_dep_3D_both_err,x_center,y_center,z_center,'cubic');
[S2_h3_c S2_h3_c_err]=histfitlandau(h3_s2_phe_both_xyz(inrange(h3_drift_time,[h3z_center-30 h3z_center+30]) & inrange(h3_x,[h3x_center-5 h3x_center+5])  & inrange(h3_y,[h3y_center-5 h3y_center+5])), 100, 400, 6000); %Just average corrects CH3T over detector, since we claim there is no field dependence
S2_h3_c_err=10/S2_h3_c; %replace error measurement based on systematic errors involved with changing the bin_max value
r_c=1-(S2_kr_c/S2_h3_c)*(NionsH3/NionsKr)*(1-r_h3);
s1p_s1_3D_s1as1b_bins=(a+1-s2p_s2_3D_s1as1b_bins.*(1-r_c))/(a+r_c);
r_c_err=sqrt( (S2_kr_c_err*(NionsH3/NionsKr)*(1-r_h3)/S2_h3_c).^2 + (S2_h3_c_err*(S2_kr_c/S2_h3_c.^2)*(NionsH3/NionsKr)*(1-r_h3)).^2 );
s1p_s1_3D_s1as1b_bins_err=sqrt( ((s2p_s2_3D_s1as1b_bins_err.*(1-r_c))/(a+r_c)).^2 + (r_c_err.*s2p_s2_3D_s1as1b_bins.*(1+a)/(a+r_c).^2).^2 );


%% Compare 3D to Z only


%Reshaping to fit polynomials
s1ab_xyz_mean_vector=reshape(s1ab_xyz_mean,[size(s1ab_xyz_mean,1)*size(s1ab_xyz_mean,2)*size(s1ab_xyz_mean,3) 1]);
s1p_s1_3D_s1as1b_bins_vector=reshape(s1p_s1_3D_s1as1b_bins,[size(s1p_s1_3D_s1as1b_bins,1)*size(s1p_s1_3D_s1as1b_bins,2)*size(s1p_s1_3D_s1as1b_bins,3) 1]);
s2p_s2_3D_s1as1b_bins_vector=reshape(s2p_s2_3D_s1as1b_bins,[size(s2p_s2_3D_s1as1b_bins,1)*size(s2p_s2_3D_s1as1b_bins,2)*size(s2p_s2_3D_s1as1b_bins,3) 1]);

s1ab_z_means_vector=reshape(s1ab_z_means,[size(s1ab_z_means,1)*size(s1ab_z_means,2)*size(s1ab_z_means,3) 1]);
s1p_s1_z_s1as1b_bins_vector=reshape(s1p_s1_z_s1as1b_bins,[size(s1p_s1_z_s1as1b_bins,1)*size(s1p_s1_z_s1as1b_bins,2)*size(s1p_s1_z_s1as1b_bins,3) 1]);
s2p_s2_z_s1as1b_bins_vector=reshape(s2p_s2_z_s1as1b_bins,[size(s2p_s2_z_s1as1b_bins,1)*size(s2p_s2_z_s1as1b_bins,2)*size(s2p_s2_z_s1as1b_bins,3) 1]);


%removing NaN
s1p_s1_3D_s1as1b_bins_vector(isnan(s1ab_xyz_mean_vector))=[];
s2p_s2_3D_s1as1b_bins_vector(isnan(s1ab_xyz_mean_vector))=[];
s1ab_xyz_mean_vector(isnan(s1ab_xyz_mean_vector))=[];

s2p_s2_z_s1as1b_bins_vector(isnan(s1ab_z_means_vector))=[];
s1p_s1_z_s1as1b_bins_vector(isnan(s1ab_z_means_vector))=[];
s1ab_z_means_vector(isnan(s1ab_z_means_vector))=[];

%fitting polynomials
s1ab_s1_xyz_fit=fit(s1ab_xyz_mean_vector,s1p_s1_3D_s1as1b_bins_vector,'poly2');
s1ab_s2_xyz_fit=fit(s1ab_xyz_mean_vector,s2p_s2_3D_s1as1b_bins_vector,'poly2');
s1ab_s1_z_fit=fit(s1ab_z_means_vector,s1p_s1_z_s1as1b_bins_vector,'poly2');
s1ab_s2_z_fit=fit(s1ab_z_means_vector,s2p_s2_z_s1as1b_bins_vector,'poly2');

%fitting lines
s1ab_s1_xyz_fit_linear=fit(s1ab_xyz_mean_vector,s1p_s1_3D_s1as1b_bins_vector,'poly1');
s1ab_s2_xyz_fit_linear=fit(s1ab_xyz_mean_vector,s2p_s2_3D_s1as1b_bins_vector,'poly1');
s1ab_s1_z_fit_linear=fit(s1ab_z_means_vector,s1p_s1_z_s1as1b_bins_vector,'poly1');
s1ab_s2_z_fit_linear=fit(s1ab_z_means_vector,s2p_s2_z_s1as1b_bins_vector,'poly1');

%Both at same time
figure
hold on;
rkploterr(s1ab_xyz_mean,s1p_s1_3D_s1as1b_bins,s1ab_xyz_mean_err,s1p_s1_3D_s1as1b_bins_err,[0.7 0.7 0.7],'.',100, 2);
rkploterr(s1ab_xyz_mean,s2p_s2_3D_s1as1b_bins,s1ab_xyz_mean_err,s2p_s2_3D_s1as1b_bins_err,[1 0.7 0.7],'.',100, 2);
rkploterr(s1ab_z_means,s1p_s1_z_s1as1b_bins,s1ab_z_means_err,[],[0 0 0],'.',100, 2);
rkploterr(s1ab_z_means,s2p_s2_z_s1as1b_bins,s1ab_z_means_err,[],[1 0 0],'.',100, 2);
x=0:.01:3.1;
y=s1ab_s1_xyz_fit.p1.*x.^2+s1ab_s1_xyz_fit.p2.*x+s1ab_s1_xyz_fit.p3;
plot(x,y,'Color',[0.7 0.7 0.7],'LineWidth',2);
y=s1ab_s2_xyz_fit.p1.*x.^2+s1ab_s2_xyz_fit.p2.*x+s1ab_s2_xyz_fit.p3;
plot(x,y,'Color',[1 0.7 0.7],'LineWidth',2);
y=s1ab_s1_z_fit.p1.*x.^2+s1ab_s1_z_fit.p2.*x+s1ab_s1_z_fit.p3;
plot(x,y,'Color',[0 0 0],'LineWidth',2);
y=s1ab_s2_z_fit.p1.*x.^2+s1ab_s2_z_fit.p2.*x+s1ab_s2_z_fit.p3;
plot(x,y,'Color',[1 0 0],'LineWidth',2);
xlabel('S1a/S1b'); ylabel('CH3T Corrected Kr at XYZ / CH3T Corrected Kr at Center');
myfigview(16);
legend('S1 - 3D Dependence','S2 - 3D Dependence','S1 - Z Dependence Only','S2 - Z Dependence Only');
xlim([2.3 3.2]);ylim([0.5 1.45]);

%Both at same time - S2 only
figure
hold on;
% rkploterr(s1ab_xyz_mean,s1p_s1_3D_s1as1b_bins,s1ab_xyz_mean_err,s1p_s1_3D_s1as1b_bins_err,[0.7 0.7 0.7],'.',100, 2);
rkploterr(s1ab_xyz_mean,s2p_s2_3D_s1as1b_bins,s1ab_xyz_mean_err,s2p_s2_3D_s1as1b_bins_err,[1 0.7 0.7],'.',100, 2);
% rkploterr(s1ab_z_means,s1p_s1_z_s1as1b_bins,s1ab_z_means_err,[],[0 0 0],'.',100, 2);
rkploterr(s1ab_z_means,s2p_s2_z_s1as1b_bins,s1ab_z_means_err,[],[1 0 0],'.',100, 2);
x=0:.01:3.1;
% y=s1ab_s1_xyz_fit.p1.*x.^2+s1ab_s1_xyz_fit.p2.*x+s1ab_s1_xyz_fit.p3;
% plot(x,y,'Color',[0.7 0.7 0.7],'LineWidth',2);
y=s1ab_s2_xyz_fit.p1.*x.^2+s1ab_s2_xyz_fit.p2.*x+s1ab_s2_xyz_fit.p3;
plot(x,y,'Color',[1 0.7 0.7],'LineWidth',2);
% y=s1ab_s1_z_fit.p1.*x.^2+s1ab_s1_z_fit.p2.*x+s1ab_s1_z_fit.p3;
% plot(x,y,'Color',[0 0 0],'LineWidth',2);
y=s1ab_s2_z_fit.p1.*x.^2+s1ab_s2_z_fit.p2.*x+s1ab_s2_z_fit.p3;
plot(x,y,'Color',[1 0 0],'LineWidth',2);
xlabel('S1a/S1b'); ylabel('S2_E (xyz) /S2_E (Center)');
myfigview(16);
legend('S2 - 3D Dependence','S2 - Z Dependence Only');
xlim([2.3 3.2]);ylim([0.5 1.45]);

% For plotting Chi2 result
% a1=0.05; b1=0;y=a1.*x.^2+b1.*x+(1-a1*(2.755^2)-b1*2.755);
% plot(x,y,'-b')

% For comparing saved .mat files
% figure
% hold on;
% rkploterr(s1ab_xyz_mean_map,s1p_s1_3D_s1as1b_bins_map,s1ab_xyz_mean_err_map,s1p_s1_3D_s1as1b_bins_err_map,[0.7 0.7 0.7],'.',100, 2);
% rkploterr(s1ab_xyz_mean_map,s2p_s2_3D_s1as1b_bins_map,s1ab_xyz_mean_err_map,s2p_s2_3D_s1as1b_bins_err_map,[1 0.7 0.7],'.',100, 2);
% rkploterr(s1ab_z_means_map,s1p_s1_z_s1as1b_bins_map,s1ab_z_means_err_map,s1p_s1_z_s1as1b_bins_err_map,[0 0 0],'.',100, 2);
% rkploterr(s1ab_z_means_map,s2p_s2_z_s1as1b_bins_map,s1ab_z_means_err_map,s2p_s2_z_s1as1b_bins_err_map,[1 0 0],'.',100, 2);
% x=0:.01:3.1;
% y=s1ab_s1_xyz_fit_map.p1.*x.^2+s1ab_s1_xyz_fit_map.p2.*x+s1ab_s1_xyz_fit_map.p3;
% plot(x,y,'Color',[0.7 0.7 0.7],'LineWidth',2);
% y=s1ab_s2_xyz_fit_map.p1.*x.^2+s1ab_s2_xyz_fit_map.p2.*x+s1ab_s2_xyz_fit_map.p3;
% plot(x,y,'Color',[1 0.7 0.7],'LineWidth',2);
% y=s1ab_s1_z_fit_map.p1.*x.^2+s1ab_s1_z_fit_map.p2.*x+s1ab_s1_z_fit_map.p3;
% plot(x,y,'Color',[0 0 0],'LineWidth',2);
% y=s1ab_s2_z_fit_map.p1.*x.^2+s1ab_s2_z_fit_map.p2.*x+s1ab_s2_z_fit_map.p3;
% plot(x,y,'Color',[1 0 0],'LineWidth',2);
% xlabel('S1a/S1b'); ylabel('CH3T Corrected Kr at XYZ / CH3T Corrected Kr at Center');
% myfigview(16);
% legend('S1 - 3D Dependence','S2 - 3D Dependence','S1 - Z Dependence Only','S2 - Z Dependence Only');
% 

%fitting polynomials
% s1ab_s1_xyz_fit=fit(s1ab_xyz_mean_vector,s1p_s1_3D_s1as1b_bins_vector,'poly1');
% s1ab_s2_xyz_fit=fit(s1ab_xyz_mean_vector,s2p_s2_3D_s1as1b_bins_vector,'poly1');
% s1ab_s1_z_fit=fit(s1ab_z_means_vector,s1p_s1_z_s1as1b_bins_vector,'poly1');
% s1ab_s2_z_fit=fit(s1ab_z_means_vector,s2p_s2_z_s1as1b_bins_vector,'poly1');

%giving variables unique names from KrypCal code
s1ab_xyz_mean_map=s1ab_xyz_mean;
s1p_s1_3D_s1as1b_bins_map=s1p_s1_3D_s1as1b_bins;
s1ab_xyz_mean_err_map=s1ab_xyz_mean_err;
s1p_s1_3D_s1as1b_bins_err_map=s1p_s1_3D_s1as1b_bins_err;
s2p_s2_3D_s1as1b_bins_map=s2p_s2_3D_s1as1b_bins;
s2p_s2_3D_s1as1b_bins_err_map=s2p_s2_3D_s1as1b_bins_err;
s1ab_z_means_map=s1ab_z_means;
s1p_s1_z_s1as1b_bins_map=s1p_s1_z_s1as1b_bins;
s1ab_z_means_err_map=s1ab_z_means_err;
% s1p_s1_z_s1as1b_bins_err_map=s1p_s1_z_s1as1b_bins_err;
s2p_s2_z_s1as1b_bins_map=s2p_s2_z_s1as1b_bins;
% s2p_s2_z_s1as1b_bins_err_map=s2p_s2_z_s1as1b_bins_err;
s1ab_s1_xyz_fit_map=s1ab_s1_xyz_fit;
s1ab_s2_xyz_fit_map=s1ab_s2_xyz_fit;
s1ab_s1_z_fit_map=s1ab_s1_z_fit;
s1ab_s2_z_fit_map=s1ab_s2_z_fit;

save('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p21\TritiumMethod\libNEST_CH3TModel\S1aS1bfieldmapping_2p21_2015polyfit','s1ab_xyz_mean_map','s1p_s1_3D_s1as1b_bins_map','s1ab_xyz_mean_err_map','s1p_s1_3D_s1as1b_bins_err_map',...
    's2p_s2_3D_s1as1b_bins_map','s2p_s2_3D_s1as1b_bins_err_map','s1ab_z_means_map','s1p_s1_z_s1as1b_bins_map','s1ab_z_means_err_map', ...
    's2p_s2_z_s1as1b_bins_map','s1ab_s1_xyz_fit_map','s1ab_s2_xyz_fit_map','s1ab_s1_z_fit_map','s1ab_s2_z_fit_map','s1ab_s2_z_fit_linear','s1ab_s1_z_fit_linear','s1ab_s2_xyz_fit_linear','s1ab_s1_xyz_fit_linear');

% save('CH3T_field_normalization_DP21_ConvertS2toS1','s2z_normalization_fit_bot','s2z_normalization_fit_both','norm_S1_both','norm_S1_bottom','norm_S2_xy_bot','norm_S2_xy_both','s2_xy_xbins_bot','s2_xy_xbins_both','s2_xy_ybins_bot',...
%     's2_xy_ybins_both','mean_fit_s1z','mean_fit_s1z_bottom','s1_xy_xbins','s1_xy_xbins_bottom','s1_xy_ybins','s1_xy_ybins_bottom');
