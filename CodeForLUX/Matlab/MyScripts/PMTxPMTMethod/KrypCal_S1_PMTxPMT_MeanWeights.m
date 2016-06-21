%Have to use inpainting to avoid zeros in the interpolation
%Can not use cubic interpolation and set NaN -> sum of S1, since there are too many NaN
%The typical S1 correction using a spline that works a little better, so it handles corrections outside the fiducial better
%This improves resoultion inside the fiducial, where we interpolate instead of extrapolate
%This is also much more computationally intensive than the other correction method

%% Get 3D S1 correction 


rqs_to_load = {'pulse_area_phe','full_evt_area_phe','peak_area_phe','event_timestamp_samples','aft_t1_samples','event_number','pulse_classification','top_bottom_ratio','livetime_latch_samples','livetime_end_samples','pulse_start_samples','pulse_end_samples','z_drift_samples','x_cm','y_cm','aft_t0_samples','aft_t2_samples','s1s2_pairing','s1_before_s2','x_corrected','y_corrected'};

path='C:\\Users\\Richard\\Desktop\\LUX Analysis\\TestData';

d = LUXLoadMultipleRQMs_framework(path,rqs_to_load);

d.s1s2_pairing(isnan(d.s1s2_pairing))=0; %by default it returns nan when there is no pair
d.z_drift_samples(isnan(d.z_drift_samples)) = 0.0; % get rid of NaN   

s1s2_pairing = d.s1s2_pairing;
num_paired_pulses = sum(s1s2_pairing);
num_S1_before_S2 = sum(d.s1_before_s2);
      
s1area_bound_min = 150;%150
s1area_bound_max = 500;%500

s2area_bound_min = 1000;%1000
s2area_bound_max = 20000;%20000
    
events=sort(size(d.pulse_area_phe)); %The first element is the number of pulses. sometimes it's 5, sometimes it's 10

s1_area_cut= inrange(d.pulse_area_phe,[s1area_bound_min,s1area_bound_max]);
s2_area_cut= inrange(d.pulse_area_phe,[s2area_bound_min,s2area_bound_max]);

s1_class=d.pulse_classification==1 & s1_area_cut ; %area cut for Kr events
s2_class=d.pulse_classification==2 & s2_area_cut ;      

s1_before_s2_cut = repmat(num_paired_pulses ==2 ,events(1),1) & repmat(num_S1_before_S2 == 1,events(1),1) ...
                  & repmat(sum(s1_class,1)>=1,events(1),1) & repmat(sum(s2_class,1)==1,events(1),1)...
                  & cumsum(s1_class+s2_class)<=2 ; % ignore all S1 that happen after the s2 !

s1_single_cut = s1_class & s1_before_s2_cut; %
s2_single_cut = s2_class & s1_before_s2_cut;


drift_time = d.z_drift_samples(s2_single_cut)/100;  
original_drift_time=drift_time;


%%%% Detector Edge Detection %%%%%%%

%Finding R_cut
edge_r_max=sqrt(4000*(25^2)/length(drift_time));

if edge_r_max <2
    edge_r_max=2;
end

% edge_r_cut=inrange(drift_radius,[0,edge_r_max]);
bin_size=2;

drift_time_bins=0:bin_size:1000;
drift_time_hist=hist(drift_time,drift_time_bins);



%Finds first 10 consecutive zeros... replace this with first 10 less than
%10% of average at start
zero_count=zeros(length(drift_time_hist));
mean_drift_time=mean(drift_time_hist(30:80));


for i=2:length(drift_time_hist)
    if drift_time_hist(i-1)< 0.50*mean_drift_time;
        zero_count(i)=1+zero_count(i-1);
        if zero_count(i)==10;
            det_edge=(i-11)*bin_size;
        end
    end
end


num_bins=(length(drift_time)/900)^(1/3);
r_max=25; 
xy_step=r_max*2/(num_bins-1);% Change to depend on size of pulse_area_phe
s1xbins=-r_max+xy_step/2:xy_step:r_max-xy_step/2;
s1ybins=-r_max+xy_step/2:xy_step:r_max-xy_step/2;
z_step=ceil( (0.95*det_edge-0.05*det_edge)/num_bins );  %xy_step size about = 30 mm, to create 16 z bins
s1zbins=floor(0.05*det_edge):z_step:floor(0.05*det_edge)+z_step*num_bins; %20 us


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculating light maps for each PMT Group ');

tic

s2=d.pulse_area_phe(s2_single_cut);
s2x = d.x_cm(s2_single_cut);
s2y = d.y_cm(s2_single_cut);
radius=sqrt(s2x.^2+s2y.^2);


s2=s2(radius<25);
s2_cut=inrange(s2,[200 30000]);%cut on appropriate S2 size


s2x = s2x(radius<25);
s2y = s2y(radius<25);

drift_time=original_drift_time(radius<25);
drift_time=drift_time(s2_cut);
dT=drift_time;


xc=s2x(s2_cut);
yc=s2y(s2_cut);

%%%%%%%%%%%%%%%%%%%%%%
s1_top=d.pulse_area_phe(s1_single_cut);
s1_top=s1_top(radius<25);
s1_top=s1_top(s2_cut);

hold_s1_top=s1_top;

clear q temp_norm_s1_both_xyz center_phe_both_xyz Sigma_S1_3D_bottom Sigma_S1_3D_both Sig_b Sig Mean_S1_3D_bottom Mean_S1_3D_both;

i_min=-r_max+xy_step; %#ok<*NASGU>
j_min=-r_max+xy_step;
Mean_S1_3D=zeros(length(s1xbins),length(s1ybins),length(s1zbins));

for k = s1zbins; % out to 485 mm % start at 5 mm. steps of 30 mm
    for i = (-r_max+xy_step):xy_step:r_max %map to rows going down (y_cm)
        for j = (-r_max+xy_step):xy_step:r_max %map columns across (x_cm)

            l=int8(i/xy_step+(r_max)/xy_step); %real y from 1:25
            m=int8(j/xy_step+(r_max)/xy_step); % real x from 1:25
            n=(k-floor(0.05*det_edge))/z_step + 1; %z from 1:16
            
           
           clear q;
           q = xc<j & xc>(j-xy_step) & yc<i & yc>(i-xy_step) & inrange(dT,[(k-z_step/2) , (k+z_step/2)]); %no 0th element
           
           if (length(s1_top(q)) >= 30)
            Mean_S1_3D(l,m,n)=mean(s1_top(q));
            clear temp_s1
           else

            Mean_S1_3D(l,m,n)=0;
           end
            
                   if(strcmp(num2str(Mean_S1_3D(l,m,n)),'NaN'))
                      Mean_S1_3D(l,m,n)=0;
                   end
                  
                   
        end
    end
end

Mean_S1_3D(Mean_S1_3D==0)=nan;  
Mean_S1_3D=inpaint_nans3(Mean_S1_3D);

summed_center=interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D,0,0,(det_edge-4)/2,'spline');
clear Mean_S1_3D
%%%%%%%%%%%%%%%%%%%%%%%%

counter=1;
energy_estimate=zeros(length(xc(dT>(floor(0.05*det_edge)-z_step/2) & dT<=(floor(0.05*det_edge)+z_step*10+z_step/2))),1);
variance=zeros(length(xc(dT>(floor(0.05*det_edge)-z_step/2) & dT<=(floor(0.05*det_edge)+z_step*10+z_step/2))),1);
means=zeros(length(xc(dT>(floor(0.05*det_edge)-z_step/2) & dT<=(floor(0.05*det_edge)+z_step*10+z_step/2))),1);

for jj=1:122;
 
    if jj ~= 5 && jj ~= 32 && jj~= 93; 
        clear peak_areas s1 temp_energy_est temp_var normalized_s1
        peak_areas=d.peak_area_phe(:,jj,:);

        drift_time=original_drift_time(radius<25);
        drift_time=drift_time(s2_cut);
        dT=drift_time;

        xc=s2x(s2_cut);
        yc=s2y(s2_cut);
        s1=peak_areas(s1_single_cut);
        s1=s1(radius<25);
        s1=s1(s2_cut);
        s1(s1<=0)=0;


        clear q temp_norm_s1_both_xyz center_phe_both_xyz Sigma_S1_3D_bottom Sigma_S1_3D_both Sig_b Sig Mean_S1_3D_bottom Mean_S1_3D;

        Mean_S1_3D=zeros(length(s1xbins),length(s1ybins),length(s1zbins));
        stand_dev=zeros(length(s1xbins),length(s1ybins),length(s1zbins));

        for k = s1zbins; % out to 485 mm % start at 5 mm. steps of 30 mm
            for i = (-r_max+xy_step):xy_step:r_max %map to rows going down (y_cm)
                for j = (-r_max+xy_step):xy_step:r_max %map columns across (x_cm)

                     l=int8(i/xy_step+(r_max)/xy_step); %real y from 1:25
                     m=int8(j/xy_step+(r_max)/xy_step); % real x from 1:25
                     n=int8((k-floor(0.05*det_edge))/z_step + 1); %z from 1:16

                   clear q temp_s1;
                   q = xc<=j & xc>(j-xy_step) & yc<=i & yc>(i-xy_step) & dT <=(k+z_step/2) & dT>(k-z_step/2); %no 0th element
                   temp_s1=s1(q);

                    if (length(temp_s1(temp_s1>0)) >= 30)
                       Mean_S1_3D(l,m,n)=mean(temp_s1(temp_s1>0));
                       stand_dev(l,m,n)=std(s1(q));
                    else
                       Mean_S1_3D(l,m,n)=0;
                       stand_dev(l,m,n)=0;   
                    end           

                end
            end
        end

        Mean_S1_3D(Mean_S1_3D==0)=nan;  
        Mean_S1_3D=inpaint_nans3(Mean_S1_3D);

        stand_dev(stand_dev==0)=nan;  
        stand_dev=inpaint_nans3(stand_dev);

        temp_energy_est=summed_center*s1./interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D,xc,yc,dT,'spline'); 
        temp_var=interp3(s1xbins,s1ybins,s1zbins,stand_dev,xc,yc,dT,'spline').^2;
        temp_means=interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D,xc,yc,dT,'spline').^2;

        energy_estimate(:,counter)=temp_energy_est; %Only events where the bin had >30 nonzero events for the PMT have nonzero entries.  If the PMT doesn't see >30 nonzero events, we say it didn't see anything at all
        variance(:,counter)=temp_var;
        means(:,counter)=temp_means;

        counter=counter+1;
    end
disp(['Finished PMT ' num2str(jj)]);
end


disp('Calculating weighted arethmetic means')
s1_corrected=zeros(length(s1),1);
for i=1:length(s1);
    clear temp_energy_estimate temp_variance
    
    temp_energy_estimate=energy_estimate(i,:);
    temp_variance=variance(i,:);
    temp_means=means(i,:);
    
    temp_variance(temp_energy_estimate==0)=[];
    temp_means(temp_energy_estimate==0)=[];
    temp_energy_estimate(temp_energy_estimate==0)=[];

s1_corrected(i)=sum(temp_energy_estimate.*((temp_means.^2)./temp_variance))./sum((temp_means.^2)./temp_variance); %should peak at summed_center keV
end


%% KrypCal Method -- added 5/19/2016

xc=s2x;
yc=s2y;
dT=drift_time;

clear q temp_norm_s1_both_xyz center_phe_both_xyz Sigma_S1_3D_bottom Sigma_S1_3D_both Sig_b Sig Mean_S1_3D_bottom Mean_S1_3D_both;
bin=0.1; %was 5
x=100:bin:500;
x=x';
cut_fit=x>100 & x<500; 


i_min=-r_max+xy_step;
j_min=-r_max+xy_step;
s1_phe_both=d.pulse_area_phe(s1_single_cut);
s1_phe_both=s1_phe_both(radius<25);
s1_phe_both=s1_phe_both(s2_cut);


for k = s1zbins; % out to 485 mm % start at 5 mm. steps of 30 mm
    for i = i_min:xy_step:r_max %map to rows going down (y_cm)
        for j = j_min:xy_step:r_max %map columns across (x_cm)
            l=int8(i/xy_step+(r_max)/xy_step); %real y from 1:25
            m=int8(j/xy_step+(r_max)/xy_step); % real x from 1:25
            n=(k-floor(0.05*det_edge))/z_step + 1; %z from 1:16
            
           %sort, and make the cut. using the variable q
           q = xc<j & xc>(j-xy_step) & yc<i & yc>(i-xy_step) & inrange(dT,[(k-z_step/2) , (k+z_step/2)]); %no 0th element!
                    %swap x,y due to matrix x-row, y-column definition 
                    
           % Hist 
            hist_q = hist(s1_phe_both(q) ,x)'/bin;
            
            %Count the number of events per bin
            Count_S1_3D(l,m,n)=sum(hist_q)*bin; %#ok<*SAGROW>
            
             if (Count_S1_3D(l,m,n) >= 30) % at least 30 counts before fitting. 

                yqfit=hist_q;
                    amp_start=max(yqfit(cut_fit));
                    mean_start=sum(yqfit(cut_fit).*x(cut_fit))/sum(yqfit(cut_fit));
                    sigma_start=std(yqfit(cut_fit).*x(cut_fit));    
                Fit_q= fit(x(cut_fit),yqfit(cut_fit),'gauss1','start',[amp_start mean_start sigma_start]);

                %Peak location (mean)
                Mean_S1_3D_both(l,m,n)=Fit_q.b1;

                %1-sigma of peak position
                Sig_b=Fit_q.c1/sqrt(2)/sqrt(Count_S1_3D(l,m,n));


                    %error checking
                if(strcmp(num2str(Sig_b),'NaN'))
                  Mean_S1_3D_both(l,m,n)=0;
                  Sig_b=0;
                end

                %uncertainty in mean
                Sigma_S1_3D_both(l,m,n)=Sig_b;

             else %not enough stats to do the fit

                Fit_q= 0;
                %Find Peak location
                Mean_S1_3D_both(l,m,n)=0;
                %1-sigma
                Sigma_S1_3D_both(l,m,n)=0;
             end
                   
        end
    end
end

s1xbins=-r_max+xy_step/2:xy_step:r_max-xy_step/2;
s1ybins=-r_max+xy_step/2:xy_step:r_max-xy_step/2;

center_phe_both_xyz=interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D_both,0,0,(det_edge-4)/2,'spline');%Normalize to the center (x=y=0. dT=160us)
norm_s1_both_xyz=center_phe_both_xyz./Mean_S1_3D_both; %Normalize to center (x=y=0. dT=160)
norm_s1_both_xyz(isinf(norm_s1_both_xyz))=1;%remove infinity. no correction outside 25cm
norm_s1_both_xyz(isnan(norm_s1_both_xyz))=1;%remove nan. no correction outside 25cm

s1_phe_both_xyz=s1_phe_both.*interp3(s1xbins,s1ybins,s1zbins,norm_s1_both_xyz,s2x,s2y,dT,'spline');


%Compare the two results
figure
step(s1_corrected(radius(s2_cut)<18 & inrange(dT,[40 260])),[150:1:600],'-r') %#ok<*NBRAK>
hold on;
step(s1_phe_both_xyz(radius(s2_cut)<18 & inrange(dT,[40 260])),[150:1:600],'-k')
xlabel('Corrected S1 (phe)'); ylabel('Count');
myfigview(20);
legend('PMTxPMT Map','KrypCal');


%Compare the two results with No Fiducial
figure
step(s1_corrected,[150:1:600],'-r')
hold on;
step(s1_phe_both_xyz,[150:1:600],'-k')
xlabel('Corrected S1 (phe)'); ylabel('Count');
myfigview(20);
legend('PMTxPMT Map','KrypCal');
