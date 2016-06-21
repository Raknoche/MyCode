load('Run03_Dec_Golden_CH3T.mat')

radius=(x_corr.^2 + y_corr.^2).^(1/2);
s1_min=0;
s1_max=50;
    
%     liquid_cut= inrange(drift_time,[30 300]) & inrange(s2radius_c,[0,18]) & inrange(s1area,[s1_min s1_max]) ...
%     & (log10(s1area)+0.5*log10(s2area)<3.8) & s2_phe_both>100 ; %used to be drift time 100 to 250

   
    band_sigma=4; %4 sigma = 1 in 15,000 chance of being outside this
    
             
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
   
        
            ER_fit_cut=inrange(spikyS1,[ER_bin-ER_bin_size/2,ER_bin+ER_bin_size/2]);
             
            ER_band_hist=hist(log10(s2area(ER_fit_cut)./spikyS1(ER_fit_cut)),s2_s1_bin);
            counts_ER_hist(index_i)=sum(ER_band_hist);
           
            Fit_ER=fit(s2_s1_bin',ER_band_hist','gauss1');
             
            rms_fit=fit_ML_normal(log10(s2area(ER_fit_cut)./spikyS1(ER_fit_cut)),s2_s1_bin);
            
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

 scatter(spikyS1,log10(s2area./spikyS1),'.k');
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
%   %text(5,.7,strcat('NEST\_SE_{bot} Fit=',num2str(ER_mean_NEST_fit.SE_ER_fit*73/64*27/24.5,3),'\pm ', num2str(sigma_SE_ER_fit*73/64*27/24.5,2),'[Phe]'),'fontsize',16,'color','m','interpreter','tex');
%   %text(5,.8,'NEST Platinum V4-b', 'fontsize',16,'color','m','interpreter','tex' );
%   text(5,.55,strcat('+',num2str(band_sigma),'\sigma = (',num2str(ER_upper_power_fit.a,3),'\pm', num2str(sigma_a_ER_upper,1) ,')\times S1_c ^{', num2str(ER_upper_power_fit.b,3),'\pm',num2str(sigma_b_ER_upper,1),'}'),'fontsize',16,'interpreter','tex');
%   text(5,.3,strcat('Mean=(',num2str(ER_mean_power_fit.a,3),'\pm', num2str(sigma_a_ER_mean,1),')\times S1_c ^{', num2str(ER_mean_power_fit.b,3),'\pm', num2str(sigma_b_ER_mean,1),'}'),'fontsize',16,'interpreter','tex');
%   text(5,.05,strcat('-',num2str(band_sigma),'\sigma = (',num2str(ER_lower_power_fit.a,3),'\pm', num2str(sigma_a_ER_lower,1),')\times S1_c ^{', num2str(ER_lower_power_fit.b,3),'\pm', num2str(sigma_b_ER_lower,1),'}'),'fontsize',16,'interpreter','tex');
%   

  
%% Cutting on populations above and below the band width

below_band_cut=log10(s2area./spikyS1)< (ER_lower_power_fit.a.*(spikyS1.^ER_lower_power_fit.b));
above_band_cut=log10(s2area./spikyS1)> (ER_upper_power_fit.a.*(spikyS1.^ER_upper_power_fit.b));
in_band_cut=log10(s2area./spikyS1)< (ER_upper_power_fit.a.*(spikyS1.^ER_upper_power_fit.b)) & log10(s2area./spikyS1)> (ER_lower_power_fit.a.*(spikyS1.^ER_lower_power_fit.b));

%Categories of events
hold on;
plot(spikyS1(below_band_cut),log10(s2area(below_band_cut)./spikyS1(below_band_cut)),'+m','MarkerSize',10,'LineWidth',2);
plot(spikyS1(above_band_cut),log10(s2area(above_band_cut)./spikyS1(above_band_cut)),'+b','MarkerSize',10,'LineWidth',2);


%% New version that classifies based on hand scan
%             
%   %Make Plot    
% 
%  no_fid_new=figure;
%  hold on
% 
%  scatter(spikyS1,log10(s2area./spikyS1),'.k');
%  ylabel('Corrected log10(S2/S1)'); xlabel('S1 corrected (Pulse Area Phe)'); 
%  title('Corrected CH3T','fontsize',16,'Interpreter','none');
%  myfigview(16);
%  ylim([0 3]);xlim([0 s1_max]);
% %  plot(ER_bins,lower_bound,'ok','markersize',6,'linewidth',2)
%  plot(ER_bins,ER_lower_fit,'--r','markersize',14,'linewidth',2)
% %  plot(ER_bins,mean_ER_band,'ok','markersize',6,'linewidth',2)%plot mean of ER band
%  plot(ER_bins,ER_mean_fit,'-r','markersize',20,'linewidth',2)
% %  plot(ER_bins,upper_bound,'ok','markersize',6,'linewidth',2)
%  plot(ER_bins,ER_upper_fit,'--r','markersize',14,'linewidth',2)
%  bin_step=2;
% 
%             total_count=length(log10(s2area));
% 
% plot(spikyS1(not_in_vislux),log10(s2area(not_in_vislux)./spikyS1(not_in_vislux)),'+m','MarkerSize',10,'LineWidth',2);
% plot(spikyS1(seem_fine),log10(s2area(seem_fine)./spikyS1(seem_fine)),'+g','MarkerSize',10,'LineWidth',2);
% plot(spikyS1(spheinbetween),log10(s2area(spheinbetween)./spikyS1(spheinbetween)),'+b','MarkerSize',10,'LineWidth',2);
% plot(spikyS1(cutoffs2),log10(s2area(cutoffs2)./spikyS1(cutoffs2)),'+r','MarkerSize',10,'LineWidth',2);


%% Printing Luxstamps
% clear below_band_indices below_band_events
% counter=1;
%     for ii=1:length(below_band_cut);
%         if below_band_cut(1,ii)==1;
%             below_band_indices(counter)=(ii);
%             counter=counter+1;
%         end
%     end
%     
%     for ii=1:length(below_band_indices);
%         below_band_events{ii}=stamps(below_band_indices(ii),:);
%     end
%     %sprintf('%s, \n',below_band_events{:})
% 
%     clear above_band_indices
% counter=1;
%     for ii=1:length(above_band_cut);
%         if above_band_cut(1,ii)~=0;
%             above_band_indices(counter)=(ii);
%             counter=counter+1;
%         end
%     end
%     
%     for ii=1:length(above_band_indices);
%         above_band_events{ii}=stamps(above_band_indices(ii),:);
%     end
    %sprintf('%s, \n',above_band_events{:})
    
belowband_events=luxstamp(below_band_cut);
aboveband_events=luxstamp(above_band_cut);
withinband_events=luxstamp(in_band_cut);
% sprintf('%15d, \n',luxstamp(below_band_cut)')

%% Plotting ER Band with below band events labeled

%The <= >= parts are to account for rounding error
s2_cutoff_events=( (luxstamp >= 9345476575688025 & luxstamp <= 9345476575688026 ) | luxstamp == 9345636686627672 | luxstamp ==9346286182335170 ...
| (luxstamp >=9346338646205858 & luxstamp <=9346338646205859) | luxstamp == 9353211711501746 | (luxstamp >= 9353850500120662 & luxstamp <= 9353850500120663) ...
| (luxstamp >=9354500564535586 & luxstamp <=9354500564535587) ...
| luxstamp == 9355207286744840| luxstamp ==9357718732440018 );

s1_sphe_events= ( (luxstamp >= 9347665873639937 & luxstamp <= 9347665873639939) | (luxstamp >= 9353409153802245 & luxstamp <= 9353409153802247) | (luxstamp >= 9353545137965615 & luxstamp <= 9353545137965617) ...
 | (luxstamp >= 9353836499369423 & luxstamp <= 9353836499369425) | (luxstamp >= 9355338271581951 & luxstamp <= 9355338271581953) | (luxstamp >= 9356366230146049 & luxstamp <= 9356366230146051) ...
 | (luxstamp >= 9357490076460473 & luxstamp <= 9357490076460475 ) ...
 | (luxstamp >= 9358350130729310 & luxstamp <= 9358350130729312 ) | luxstamp == 9344432895860114 | (luxstamp >= 9345470424648252 & luxstamp <= 9345470424648254) ...
 | (luxstamp >= 9346148377561980 & luxstamp <= 9346148377561982 ) ...
 | (luxstamp >= 9353404051387488 & luxstamp <= 9353404051387490 ) | (luxstamp >= 9356324167913316 & luxstamp <= 9356324167913318) ...
 | (luxstamp >= 9357464869795738 & luxstamp <= 9357464869795740));

look_fine_events = ( (luxstamp >= 9344054807376460 & luxstamp <= 9344054807376462) | (luxstamp >= 9344586269804138 & luxstamp <= 9344586269804140) | luxstamp == 9345132068186502 ...
| luxstamp == 9353562076837210 | (luxstamp >= 9354685871936240 & luxstamp <= 9354685871936242) | luxstamp == 9357003913822558 | ( luxstamp >= 9358257329458430 & luxstamp <= 9358257329458432) );

figure

 no_fid_new=figure;
 hold on
 scatter(spikyS1,log10(s2area./spikyS1),'.k');
 plot(spikyS1(s2_cutoff_events),log10(s2area(s2_cutoff_events)./spikyS1(s2_cutoff_events)),'+r','MarkerSize',10,'LineWidth',2);
 plot(spikyS1(s1_sphe_events),log10(s2area(s1_sphe_events)./spikyS1(s1_sphe_events)),'+b','MarkerSize',10,'LineWidth',2);
 plot(spikyS1(look_fine_events),log10(s2area(look_fine_events)./spikyS1(look_fine_events)),'+g','MarkerSize',10,'LineWidth',2);
 ylabel('Corrected log10(S2/S1)'); xlabel('S1 corrected (Pulse Area Phe)'); 
 title('Corrected CH3T','fontsize',16,'Interpreter','none');
 myfigview(16);
 ylim([0 3]);xlim([0 s1_max]);
 plot(ER_bins,ER_lower_fit,'--r','markersize',14,'linewidth',2)
 plot(ER_bins,ER_mean_fit,'-r','markersize',20,'linewidth',2)
 plot(ER_bins,ER_upper_fit,'--r','markersize',14,'linewidth',2)
 bin_step=2;

            total_count=length(log10(s2area));
 
    
 box;
  text(10,2.8,strcat('Total Count= ',num2str(total_count)),'fontsize',16);


%% Plotting XY dist of below band population

figure
plot(x_corr,y_corr,'.k');
hold on;
plot(x_corr(above_band_cut),y_corr(above_band_cut),'+b','MarkerSize',10,'LineWidth',2);
plot(x_corr(below_band_cut),y_corr(below_band_cut),'+r','MarkerSize',10,'LineWidth',2);
circle(0,0,25);
legend('All Events','Above Band','Below Band');
xlabel('Corrected x (cm)');ylabel('Corrected y (cm)');
myfigview(16);

%% Plotting RvZ dist of below and above band events

figure
hold on;
plot(radius.^2,drift,'.k','MarkerSize',10,'LineWidth',2);
plot(radius(above_band_cut).^2,drift(above_band_cut),'+b','MarkerSize',10,'LineWidth',2);
plot(radius(below_band_cut).^2,drift(below_band_cut),'+r','MarkerSize',10,'LineWidth',2);
legend('All Events','Above Band','Below Band');
xlabel('Corrected Radius (cm^2)');ylabel('Drift Time(uSec)');
myfigview(16);

%% Plot the one unaccounted "background event" location

mycut=(inrange(spikyS1,[15.5,16.5]) & inrange(log10(s2area./spikyS1),[1.25,1.3]));

figure
plot(radius.^2,drift,'.k','MarkerSize',10,'LineWidth',2);
hold on;
plot(radius(mycut).^2,drift(mycut),'+r','MarkerSize',10,'LineWidth',2);
xlabel('Corrected Radius (cm^2)');ylabel('Drift Time(uSec)');
myfigview(16);

figure
plot(x_corr,y_corr,'.k');
hold on;
plot(x_corr(mycut),y_corr(mycut),'+r','MarkerSize',10,'LineWidth',2);
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


% sprintf('%d, \n',belowband_events)