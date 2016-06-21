%Need to run this with Gaussian and Poissons Data
%Mainly interested in the result of KrypCal, norm_finalweights, and nonorms
%Everything should be essentially the same under poisson stats, and should be better with PMTxPMT method without poisson stats

PMT_Radius=2.5; %5cm diameter PMTs

disp('Setting Up PMT Coords')
%Set up hexagonal pattern of PMTs
PMT_X_Spacing=(50-PMT_Radius*18)/10;
PMT_Y_Spacing=(50-PMT_Radius*18)/10;
PMT_row=[1,2,3,4,2,3,4,3,4,4,5,6,7,8,5,6,7,5,6,5,9,9,9,9,8,8,8,7,7,6,9,8,7,6,8,7,6,7,6,6,5,4,3,2,5,4,3,5,4,5,1,1,1,1,2,2,2,3,3,4,1,2,3,4,2,3,4,3,4,4,5,6,7,8,5,6,7,5,6,5,9,9,9,9,8,8,8,7,7,6,9,8,7,6,8,7,6,7,6,6,5,4,3,2,5,4,3,5,4,5,1,1,1,1,2,2,2,3,3,4,5,5];
PMT_column=[5,4,3,2,6,5,4,7,6,8,1,2,3,4,3,4,5,5,6,7,5,7,9,11,6,8,10,7,9,8,13,14,15,16,12,13,14,11,12,10,17,16,15,14,15,14,13,13,12,11,13,11,9,7,12,10,8,11,9,10,5,4,3,2,6,5,4,7,6,8,1,2,3,4,3,4,5,5,6,7,5,7,9,11,6,8,10,7,9,8,13,14,15,16,12,13,14,11,12,10,17,16,15,14,15,14,13,13,12,11,13,11,9,7,12,10,8,11,9,10,9,9];
PMT_Y_Coord=(5-PMT_row)*PMT_Radius*sqrt(3)+((5-PMT_row)/2)*PMT_Y_Spacing;
PMT_X_Coord=(PMT_column-9)*PMT_Radius+((PMT_column-9)/2)*PMT_X_Spacing;

 %Detector height is 360 microseconds here
PMT_Z_Coord([1:60,121])=360*.151; %360 uSec converted to cm
PMT_Z_Coord([61:120,122])=0;

%Assign random QE to each PMT from 0.3 to 0.4
PMT_QE=0.3 + (0.4-0.3).*rand(122,1);

%Generate an S1 event from a Gaussian distibution center at 10,000 with a
%Sigma of 667 (correspond to a Sigma of 20 and Mean of 300 in real data)
%and assign a x,y,z coord


Gaussian_Noise_Sigma=0.5;
% num_events=1000000;
num_events=300000;

disp('Setting Up Events')
Event_Signal=10000;%normrnd(10000,670,1,num_events); used to be 500 
Event_X_Coord=-25 + 50.*rand(num_events,1);
Event_Y_Coord=-25 + 50.*rand(num_events,1);
Event_Z_Coord=10+220.*rand(num_events,1).*0.151; %360 uSec Z Coverted into cm, events cant be within 10 cm of PMTs

%Distribute the signal among all PMTs, giving each PMT a portion of the
%gaussian center divided by ((delta theta + delta phi)/360*2) which is the
%fraction of radians the PMT face covers over a sphere from the location of
%the event

Delta_Theta=zeros(num_events,122);
Delta_Phi=zeros(num_events,122);
PMT_Signal=zeros(num_events,122);
PMT_Signal2=zeros(num_events,122);
PMT_Distance=zeros(num_events,122);
Poisson_PMT_Signal=zeros(num_events,122);
Gaussian_Noise_PMT_Signal=zeros(num_events,122);

disp('Setting Up Unattenuated PMT Signals')
% Setting up PMT Signals 

for jj=1:122

 v_center = [Event_X_Coord(:)-PMT_X_Coord(jj),Event_Y_Coord(:)-PMT_Y_Coord(jj),Event_Z_Coord(:)-PMT_Z_Coord(jj)]; % Vector from B to A
 
 if jj<=60 || jj==121
     v_norm=[0,0,-1];
     v_norm=repmat(v_norm,num_events,1);
 else
     v_norm=[0,0,1];
     v_norm=repmat(v_norm,num_events,1);
 end

  Delta_Phi(:,jj) = (180/pi).*atan2( sqrt(sum(cross(v_center(:,:),v_norm(:,:),2).^2,2)),dot(v_norm,v_center,2));
  
  PMT_Distance(:,jj) = sqrt((Event_X_Coord(:)-PMT_X_Coord(jj)).^2+(Event_Y_Coord(:)-PMT_Y_Coord(jj)).^2+(Event_Z_Coord(:)-PMT_Z_Coord(jj)).^2);

PMT_Signal(:,jj)=cosd(Delta_Phi(:,jj)).*((PMT_Radius^2)./(4.*PMT_Distance(:,jj).^2)).*Event_Signal(:).*PMT_QE(jj); %Divides signal based on fraction angle over 4*pi steradians
Poisson_PMT_Signal(:,jj)=poissrnd(PMT_Signal(:,jj));
Gaussian_Noise_PMT_Signal(:,jj)=normrnd(Poisson_PMT_Signal(:,jj),Gaussian_Noise_Sigma);

end

disp('Calculating Summed PMT Signal')
Summed_PMT_Signal=zeros(num_events,1);
Summed_Poisson_PMT_Signal=zeros(num_events,1);
Summed_Gaussian_Noise_PMT_Signal=zeros(num_events,1);
for i=1:num_events;
Summed_PMT_Signal(i)=sum(PMT_Signal(i,:));
Summed_Poisson_PMT_Signal(i)=sum(Poisson_PMT_Signal(i,:));
Summed_Gaussian_Noise_PMT_Signal(i)=sum(Gaussian_Noise_PMT_Signal(i,:));
end


summed_s1=Summed_Poisson_PMT_Signal;

num_bins=(num_events/900)^(1/3);
r_max=25; 
xy_step=r_max*2/(num_bins-1);% Change to depend on size of pulse_area_phe
det_edge=360*.151;
z_step=ceil( (0.95*det_edge-0.05*det_edge)/num_bins );  %xy_step size about = 30 mm, to create 16 z bins
s1xbins=-r_max+xy_step/2:xy_step:r_max-xy_step/2;
s1ybins=-r_max+xy_step/2:xy_step:r_max-xy_step/2;
s1zbins=floor(0.05*det_edge):z_step:floor(0.05*det_edge)+z_step*num_bins; %20 us


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %Calculating PMTxPMT light map with Normalization
% disp('Calculating PMTxPMT Light Map with Normalization')
% 
% %Calculating detector center
disp('Calculating Summed Center')
Mean_S1_3D=zeros(length(s1xbins),length(s1ybins),length(s1zbins));

for k = s1zbins; % out to 485 mm % start at 5 mm. steps of 30 mm
    for i = (-r_max+xy_step):xy_step:r_max %map to rows going down (y_cm)
        for j = (-r_max+xy_step):xy_step:r_max %map columns across (x_cm)

            clear l m n
            l=int8(i/xy_step+(r_max)/xy_step); %real y from 1:25
            m=int8(j/xy_step+(r_max)/xy_step); % real x from 1:25
            n=(k-floor(0.05*det_edge))/z_step + 1; %z from 1:16
            
           clear q;
           q = Event_X_Coord<j & Event_X_Coord>(j-xy_step) & Event_Y_Coord<i & Event_Y_Coord>(i-xy_step) & inrange(Event_Z_Coord,[(k-z_step/2) , (k+z_step/2)]); %no 0th element
           
           if (length(summed_s1(q)) >= 30)
            Mean_S1_3D(l,m,n)=mean(summed_s1(q));
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

counter=1;
energy_estimate=zeros(length(Event_X_Coord),1);
variance=zeros(length(Event_X_Coord),1);
means=zeros(length(Event_X_Coord),1);

for jj=1:122;
    
clear Mean_S1_3D;
Mean_S1_3D=zeros(length(s1xbins),length(s1ybins),length(s1zbins));

s1=Poisson_PMT_Signal(:,jj);

for k = s1zbins; % out to 485 mm % start at 5 mm. steps of 30 mm
    for i = (-r_max+xy_step):xy_step:r_max %map to rows going down (y_cm)
        for j = (-r_max+xy_step):xy_step:r_max %map columns across (x_cm)
            
                        clear l m n q
             l=int8(i/xy_step+(r_max)/xy_step); %real y from 1:25
             m=int8(j/xy_step+(r_max)/xy_step); % real x from 1:25
             n=int8((k-floor(0.05*det_edge))/z_step + 1); %z from 1:16

           q = Event_X_Coord<=j & Event_X_Coord>(j-xy_step) & Event_Y_Coord<=i & Event_Y_Coord>(i-xy_step) & Event_Z_Coord <=(k+z_step/2) & Event_Z_Coord>(k-z_step/2); %no 0th element
                           
            if (length(s1(q)) >= 30)
           Mean_S1_3D(l,m,n)=mean(s1(q));
            else
           Mean_S1_3D(l,m,n)=0;  
            end           
                   
        end
    end
end

Mean_S1_3D(Mean_S1_3D==0)=nan;  
Mean_S1_3D=inpaint_nans3(Mean_S1_3D);

norm_s1=summed_center./Mean_S1_3D; %Normalize to center (x=y=0. Event_Z_Coord=160)
norm_s1(isinf(norm_s1))=1;%remove infinity. no correction outside 25cm
norm_s1(isnan(norm_s1))=1;%remove nan. no correction outside 25cm

normalized_s1=s1.*interp3(s1xbins,s1ybins,s1zbins,norm_s1,Event_X_Coord,Event_Y_Coord,Event_Z_Coord,'spline');

disp('Calculating Summed Normalized Signal')
Summed_Norm_S1=zeros(num_events,1);
for i=1:num_events;
Summed_Norm_S1(i)=sum(normalized_s1(i,:));
end

clear Mean_S1_3D stand_dev;

Mean_S1_3D=zeros(length(s1xbins),length(s1ybins),length(s1zbins));
stand_dev=zeros(length(s1xbins),length(s1ybins),length(s1zbins));

for k = s1zbins; % out to 485 mm % start at 5 mm. steps of 30 mm
    for i = (-r_max+xy_step):xy_step:r_max %map to rows going down (y_cm)
        for j = (-r_max+xy_step):xy_step:r_max %map columns across (x_cm)
            
            clear l m n q
             l=int8(i/xy_step+(r_max)/xy_step); %real y from 1:25
             m=int8(j/xy_step+(r_max)/xy_step); % real x from 1:25
             n=int8((k-floor(0.05*det_edge))/z_step + 1); %z from 1:16
             
           q = Event_X_Coord<=j & Event_X_Coord>(j-xy_step) & Event_Y_Coord<=i & Event_Y_Coord>(i-xy_step) & Event_Z_Coord <=(k+z_step/2) & Event_Z_Coord>(k-z_step/2); %no 0th element
                           
            if (length(normalized_s1(q)) >= 30)
           Mean_S1_3D(l,m,n)=mean(normalized_s1(q));
           stand_dev(l,m,n)=std(normalized_s1(q));
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

energy_estimate(:,counter)=summed_center*normalized_s1./interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D,Event_X_Coord,Event_Y_Coord,Event_Z_Coord,'spline'); %Only events where the bin had >30 nonzero events for the PMT have nonzero entries.  If the PMT doesn't see >30 nonzero events, we say it didn't see anything at all
variance(:,counter)=interp3(s1xbins,s1ybins,s1zbins,stand_dev,Event_X_Coord,Event_Y_Coord,Event_Z_Coord,'spline').^2;
means(:,counter)=interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D,Event_X_Coord,Event_Y_Coord,Event_Z_Coord,'spline'); 

counter=counter+1;

disp(['Finished PMT ', num2str(jj)]);
end


disp('Calculating weighted arethmetic means for Previous Light Map')
s1_corrected_norms=zeros(length(normalized_s1),1);
s1_corrected_norms_meanweights=zeros(length(normalized_s1),1);
s1_corrected_norms_finalweights=zeros(length(normalized_s1),1);
for i=1:length(normalized_s1);
    clear temp_energy_estimate temp_variance
    
    temp_energy_estimate=energy_estimate(i,:);
    temp_variance=variance(i,:);
    temp_means=means(i,:);
    
s1_corrected_norms(i)=sum(temp_energy_estimate.*((temp_means.^2)./temp_variance))./sum((temp_means.^2)./temp_variance); %should peak at summed_center keV
s1_corrected_norms_meanweights(i)=sum(temp_energy_estimate.*(temp_means))./sum(temp_means); %should peak at summed_center keV
s1_corrected_norms_finalweights(i)=sum(temp_energy_estimate./temp_variance)/sum(1./temp_variance); %should peak at summed_center keV
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculating PMTxPMT light map without Normalization
disp('Calculating PMTxPMT Light Map without Normalization')

counter=1;
energy_estimate=zeros(length(Event_X_Coord),1);
variance=zeros(length(Event_X_Coord),1);
means=zeros(length(Event_X_Coord),1);


for jj=1:122;
    
clear Mean_S1_3D;

s1=Poisson_PMT_Signal(:,jj);

Mean_S1_3D=zeros(length(s1xbins),length(s1ybins),length(s1zbins));
stand_dev=zeros(length(s1xbins),length(s1ybins),length(s1zbins));

for k = s1zbins; % out to 485 mm % start at 5 mm. steps of 30 mm
    for i = (-r_max+xy_step):xy_step:r_max %map to rows going down (y_cm)
        for j = (-r_max+xy_step):xy_step:r_max %map columns across (x_cm)
            
            clear l m n q
             l=int8(i/xy_step+(r_max)/xy_step); %real y from 1:25
             m=int8(j/xy_step+(r_max)/xy_step); % real x from 1:25
             n=int8((k-floor(0.05*det_edge))/z_step + 1); %z from 1:16
             
           q = Event_X_Coord<=j & Event_X_Coord>(j-xy_step) & Event_Y_Coord<=i & Event_Y_Coord>(i-xy_step) & Event_Z_Coord <=(k+z_step/2) & Event_Z_Coord>(k-z_step/2); %no 0th element
                           
            if (length(s1(q)) >= 30)
           Mean_S1_3D(l,m,n)=mean(s1(q));
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
     
energy_estimate(:,counter)=summed_center*s1./interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D,Event_X_Coord,Event_Y_Coord,Event_Z_Coord,'spline'); %Only events where the bin had >30 nonzero events for the PMT have nonzero entries.  If the PMT doesn't see >30 nonzero events, we say it didn't see anything at all
variance(:,counter)=interp3(s1xbins,s1ybins,s1zbins,stand_dev,Event_X_Coord,Event_Y_Coord,Event_Z_Coord,'spline').^2;
means(:,counter)=interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D,Event_X_Coord,Event_Y_Coord,Event_Z_Coord,'spline'); 

counter=counter+1;

disp(['Finished PMT ', num2str(jj)]);
end


disp('Calculating weighted arethmetic means for Previous Light Map')
s1_corrected_nonorms=zeros(length(s1),1);
s1_corrected_nonorms_meanweights=zeros(length(s1),1);
s1_corrected_nonorms_finalweights=zeros(length(s1),1);
for i=1:length(s1);
    clear temp_energy_estimate temp_variance
    
    temp_energy_estimate=energy_estimate(i,:);
    temp_variance=variance(i,:);
    temp_means=means(i,:);
    
s1_corrected_nonorms(i)=sum(temp_energy_estimate.*((temp_means.^2)./temp_variance))./sum((temp_means.^2)./temp_variance); %should peak at summed_center keV
s1_corrected_nonorms_meanweights(i)=sum(temp_energy_estimate.*(temp_means))./sum(temp_means); %should peak at summed_center keV
s1_corrected_nonorms_finalweights(i)=sum(temp_energy_estimate./temp_variance)/sum(1./temp_variance); %should peak at summed_center keV

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculating KrypCal Light Map
disp('Calculating Krypcal Light Map')

%Calculating detector center
disp('Calculating Summed Center')
i_min=-r_max+xy_step;
j_min=-r_max+xy_step;
Mean_S1_3D=zeros(length(s1xbins),length(s1ybins),length(s1zbins));

for k = s1zbins; % out to 485 mm % start at 5 mm. steps of 30 mm
    for i = (-r_max+xy_step):xy_step:r_max %map to rows going down (y_cm)
        for j = (-r_max+xy_step):xy_step:r_max %map columns across (x_cm)

            clear l m n
            l=int8(i/xy_step+(r_max)/xy_step); %real y from 1:25
            m=int8(j/xy_step+(r_max)/xy_step); % real x from 1:25
            n=(k-floor(0.05*det_edge))/z_step + 1; %z from 1:16
            
           clear q;
           q = Event_X_Coord<j & Event_X_Coord>(j-xy_step) & Event_Y_Coord<i & Event_Y_Coord>(i-xy_step) & inrange(Event_Z_Coord,[(k-z_step/2) , (k+z_step/2)]); %no 0th element
           
           if (length(summed_s1(q)) >= 30)
            Mean_S1_3D(l,m,n)=mean(summed_s1(q));
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

norm_s1=summed_center./Mean_S1_3D; %Normalize to center (x=y=0. Event_Z_Coord=160)
norm_s1(isinf(norm_s1))=1;%remove infinity. no correction outside 25cm
norm_s1(isnan(norm_s1))=1;%remove nan. no correction outside 25cm

krypcal_corrected_s1=summed_s1.*interp3(s1xbins,s1ybins,s1zbins,norm_s1,Event_X_Coord,Event_Y_Coord,Event_Z_Coord,'spline');

%PLOTS
% 
% figure
% step(s1_corrected_norms,[0:5:1000],'c')
% hold on;
% step(s1_corrected_nonorms,[0:5:1000],'r')
% step(s1_corrected_nonorms_meanweights,[0:5:1000],'g')
% % step(s1_corrected_norms_meanweights,[0:5:1000],'m')
% step(krypcal_corrected_s1,[0:5:1000],'b')
% logy

% save('SimS1_V2_Gauss00')