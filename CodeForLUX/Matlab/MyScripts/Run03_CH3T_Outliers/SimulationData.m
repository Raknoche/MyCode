load('Tritium_SIM_cp14477.mat')
spikyS1=to_save(:,1);
s2area=to_save(:,2);
drift=to_save(:,4);
s1_min=0;
s1_max=50;
  
  
             
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
seem_fine=(luxstamp == 9344432895860114 | luxstamp == 9345132068186502 | luxstamp == 9353562076837210 | luxstamp == 9357003913822558 | luxstamp == 9357490076460474);
not_in_vislux=(luxstamp==9344054807376460 | luxstamp == 9344586269804140 | luxstamp ==9345470424648252 | luxstamp == 9345476575688024 | luxstamp == 9346148377561980 | luxstamp == 9346338646205860 | luxstamp == 9353404051387488 | luxstamp ==9353850500120664 | luxstamp == 9354500564535588 | luxstamp == 9354685871936240 | luxstamp == 9356324167913316 | luxstamp == 9356366230146050 | luxstamp == 9357464869795740 | luxstamp == 9358257329458432 | luxstamp == 9358350130729312 );
cutoffs2=(luxstamp==9345636686627672 | luxstamp == 9346286182335170 | luxstamp == 9353211711501746 | luxstamp==9355207286744840 | luxstamp == 9357718732440018);
spheinbetween= (luxstamp == 9347665873639938 | luxstamp == 9353409153802246 | luxstamp ==9353545137965616 | luxstamp == 9355338271581952 | luxstamp == 9353836499369424 | luxstamp == 9356324167913316 );
hold on;
plot(spikyS1(below_band_cut),log10(s2area(below_band_cut)./spikyS1(below_band_cut)),'+m','MarkerSize',10,'LineWidth',2);
plot(spikyS1(above_band_cut),log10(s2area(above_band_cut)./spikyS1(above_band_cut)),'+b','MarkerSize',10,'LineWidth',2);
