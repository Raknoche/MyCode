load('Run03_WS.mat');
load('DecERBand');

s1_min=0;
s1_max=50;
  
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

% Draw 4 sigma ER Band with WS data
 figure
 scatter(spikyS1,log10(s2area./spikyS1),'.k');
 hold on;
 ylabel('Corrected log10(S2/S1)'); xlabel('S1 corrected (Pulse Area Phe)'); 
 title('Corrected WS Data','fontsize',16,'Interpreter','none');
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
  
  
below_band_cut=log10(s2area./spikyS1)< (ER_lower_power_fit.a.*(spikyS1.^ER_lower_power_fit.b)) & spikyS1<50;
above_band_cut=log10(s2area./spikyS1)> (ER_upper_power_fit.a.*(spikyS1.^ER_upper_power_fit.b)) & spikyS1<50;
in_band_cut=log10(s2area./spikyS1)< (ER_upper_power_fit.a.*(spikyS1.^ER_upper_power_fit.b)) & log10(s2area./spikyS1)> (ER_lower_power_fit.a.*(spikyS1.^ER_lower_power_fit.b));
hold on;
plot(spikyS1(below_band_cut),log10(s2area(below_band_cut)./spikyS1(below_band_cut)),'+m','MarkerSize',10,'LineWidth',2);
plot(spikyS1(above_band_cut),log10(s2area(above_band_cut)./spikyS1(above_band_cut)),'+b','MarkerSize',10,'LineWidth',2);

%% Printing Luxstamps
belowband_events=luxstamp(below_band_cut);
aboveband_events=luxstamp(above_band_cut);
withinband_events=luxstamp(in_band_cut);

% sprintf('%d, \n',belowband_events)
