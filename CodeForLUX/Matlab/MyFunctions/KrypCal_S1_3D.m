%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This function produces XYZ dependence maps of the S1 signal and directly submits 
%       the result to the LUG.
%
%  Inputs:
%       - s1_phe_both          - S1 pulse areas (with both arrays) after golden selection
%       - s1_phe_bottom        - S1 pulse areas (with bottom arrays) after golden selection
%       - s2_phe_both          - S2 pulse areas (with both arrays) after golden selection
%       - s2x                  - X coordinate of the pulses
%       - s2y                  - Y coordinate of the pulses
%       - drift_time           - Z coordinate of the pulses
%       - x_center             - X coordinate of the detector center
%       - y_center             - Y coordinate of the detector center
%       - z_center             - Z coordinate of the detector center
%       - det_edge             - Edge of the drift time histogram
%       - user_name            - User name for the LUG submission
%       - file_id_cp           - Dateset name with cp number
%       - dp_version           - Version of KrypCal for LUG submission
%       - algorithm_name       - Algorithm name from KrypCal
%       - submit               - Option to submit to LUG (1 to submit, 0 not to)
%
%  Outputs:
%       - lug_val              - Cell with information for LUG entries
% 
%  Author:
%       - Attila Dobi, cleaned up by Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [lug_val] = KrypCal_S1_3D(s1_phe_both,s1_phe_bottom,s2_phe_both,s2x,s2y,drift_time,...
    x_center,y_center,z_center,det_edge, user_name,file_id_cp,dp_version, algorithm_name, submit)

    %Set Up Variables
    s2_cut             = inrange(s2_phe_both,[200 70000]);     %cut on appropriate S2 size
    s1                 = s1_phe_both(s2_cut);                  %No Corrections. Use Bottom, Top or Both PMT.
    s1_bottom          = s1_phe_bottom(s2_cut);
    xc                 = s2x(s2_cut);
    yc                 = s2y(s2_cut);
    dT                 = drift_time(s2_cut);
    Kr_events          = length(s1);
    bin                = 5;
    x                  = 40:bin:500;                           %Set up S1 energy bins
    x                  = x';
    cut_fit            = x>50 & x<490;                         %remove bin edges
    xy_step            = 3;                                    % 3 cm bins
    r_max              = 25;                                   % max radius
    z_step             = ceil( (0.95*det_edge-10)/15 );        %xy_step size about = 30 mm, to create 16 z bins
    s1zbins            = floor(10):z_step:floor(10)+z_step*15; %20 us

    %Preallocate memory
    Count_S1_3D        = zeros(16,16,16);
    Count_S1_3D_bottom = zeros(16,16,16);
    Mean_S1_3D_both    = zeros(16,16,16);
    Mean_S1_3D_bottom  = zeros(16,16,16);
    Sigma_S1_3D_both   = zeros(16,16,16);
    Sigma_S1_3D_bottom = zeros(16,16,16);

    %% Calculate the 3D map for both and bottom PMT arrays
    for k = s1zbins; % out to 485 mm % start at about 20 mm bin center. xy_steps of 30 mm
        tic;
        for i = (-r_max+xy_step):xy_step:r_max %map to rows going down (y_cm)
            for j = (-r_max+xy_step):xy_step:r_max %map columns across (x_cm)

                % Calculate map for Both PMT arrays 

                %Matrix indicies
                l             = int8(i/xy_step+(r_max)/xy_step); %real y
                m             = int8(j/xy_step+(r_max)/xy_step); % real x
                n             =(k-floor(10))/z_step + 1; 

                %Select events within the voxel
                q             = xc<j & xc>(j-xy_step) & yc<i & yc>(i-xy_step) & inrange(dT,[(k-z_step/2) , (k+z_step/2)]); %no 0th element!

                %Histogram of data in the voxel 
                hist_q        = hist(s1(q) ,x)'/bin;
                hist_q_bottom = hist(s1_bottom(q) ,x)'/bin;

                %memory clean up. Improve speed by cutting down the size of S1      
                s1            = s1(~q); 
                s1_bottom     = s1_bottom(~q);
                xc            = xc(~q);
                yc            = yc(~q);
                dT            = dT(~q);
                clear 'q'; 

                %Count the number of events per bin
                Count_S1_3D(l,m,n)        = sum(hist_q)*bin;
                Count_S1_3D_bottom(l,m,n) = sum(hist_q_bottom)*bin;

                %If there are enough events, fill in the map entry for the voxel
                %Note, X,Y indicies intentionally swapped due to how interp3D works
                if (Count_S1_3D(l,m,n) >= 50) % at least 50 counts before fitting. 

                    yqfit                    = hist_q;
                    amp_start                = max(yqfit(cut_fit));
                    mean_start               = sum(yqfit(cut_fit).*x(cut_fit))/sum(yqfit(cut_fit));
                    sigma_start              = std(yqfit(cut_fit).*x(cut_fit));    
                    Fit_q                    = fit(x(cut_fit),yqfit(cut_fit),'gauss1','start',[amp_start mean_start sigma_start]);

                    %Peak location (mean)
                    Mean_S1_3D_both(l,m,n)   = Fit_q.b1;

                    %1-sigma of peak position
                    Sig_b                    = Fit_q.c1/sqrt(2)/sqrt(Count_S1_3D(l,m,n));


                    %error checking
                    if(strcmp(num2str(Sig_b),'NaN'))
                      Mean_S1_3D_both(l,m,n) = 0;
                      Sig_b                  = 0;
                    end

                    %uncertainty in mean
                    Sigma_S1_3D_both(l,m,n)  = Sig_b;

                 else %not enough stats to do the fit
                    Mean_S1_3D_both(l,m,n)   = 0;
                    Sigma_S1_3D_both(l,m,n)  = 0;
                end

                %% Calculate map for bottom PMT array

                clear yqfit Sig_b

                if (Count_S1_3D_bottom(l,m,n) >= 50) % at least 10 counts before fitting. 

                    yqfit_bottom             = hist_q_bottom;
                    amp_start                = max(yqfit_bottom(cut_fit));
                    mean_start               = sum( yqfit_bottom(cut_fit).*x(cut_fit))/sum( yqfit_bottom(cut_fit));
                    sigma_start              = std( yqfit_bottom(cut_fit).*x(cut_fit));    
                    Fit_q_bottom             = fit(x(cut_fit),yqfit_bottom(cut_fit),'gauss1','start',[amp_start mean_start sigma_start]);

                    %Peak location (mean)
                    Mean_S1_3D_bottom(l,m,n) = Fit_q_bottom.b1;

                    %1-sigma of peak position
                    Sig_b                    = Fit_q_bottom.c1/sqrt(2)/sqrt(Count_S1_3D_bottom(l,m,n));


                    %error checking
                    if(strcmp(num2str(Sig_b),'NaN'))
                     Mean_S1_3D_bottom(l,m,n)= 0;
                     Sig_b                   = 0;
                    end

                    %uncertainty in mean
                    Sigma_S1_3D_bottom(l,m,n)= Sig_b;

                 else %not enough stats to do the fit               
                    Mean_S1_3D_bottom(l,m,n) = 0;
                    Sigma_S1_3D_bottom(l,m,n)= 0;

                end          

            end
        end
    toc;
    k %#ok<NOPRT>
    end

    s1xbins                                                   = -r_max+xy_step/2:xy_step:r_max-xy_step/2;
    s1ybins                                                   = -r_max+xy_step/2:xy_step:r_max-xy_step/2;

    %get the correction for the top PMT array by taking Both-Bottom
    Mean_S1_3D_top                                            = Mean_S1_3D_both-Mean_S1_3D_bottom;
    Sigma_S1_3D_top                                           = sqrt(Sigma_S1_3D_both.^2+Sigma_S1_3D_bottom.^2);

    %% S1 3D correction. Both PMT 
    center_phe_both_xyz                                       = interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D_both,x_center,y_center,z_center,'cubic');%Normalize to the center 
    norm_s1_both_xyz                                          = center_phe_both_xyz./Mean_S1_3D_both; %Normalize to center 
    norm_s1_both_xyz(isinf(norm_s1_both_xyz))                 = 1;
    norm_s1_both_xyz(isnan(norm_s1_both_xyz))                 = 1;

    %Get 1 sigma of the corrections matrix.
    sigma_center_phe_both_xyz                                 = interp3(s1xbins,s1ybins,s1zbins,Sigma_S1_3D_both,x_center,y_center,z_center,'cubic');%Normalize to the center 
    sigma_norm_s1_both_xyz                                    = sqrt((sigma_center_phe_both_xyz./Mean_S1_3D_both).^2+(Sigma_S1_3D_both.*center_phe_both_xyz./Mean_S1_3D_both.^2).^2);
    sigma_norm_s1_both_xyz(isinf(sigma_norm_s1_both_xyz))     = 1;%no infinity. Set Sigma=1, which is 100% uncertainty
    sigma_norm_s1_both_xyz(isnan(sigma_norm_s1_both_xyz))     = 1;

    mean_center_3D_both                                       = interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D_both,x_center,y_center,z_center,'cubic');
    sigma_mean_center_3D_both                                 = interp3(s1xbins,s1ybins,s1zbins,Sigma_S1_3D_both,x_center,y_center,z_center,'cubic');

    %to Apply the correction use the following.
    %s1_phe_both_xyz=s1_phe_both.*interp3(xx,yy,zz,norm_s1_both_xyz,s2x,s2y,drift_time,'cubic');

    %% S1 3D correction. Bottom PMT 
    center_phe_bottom_xyz                                     = interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D_bottom,x_center,y_center,z_center,'cubic');%Normalize to the center 
    norm_s1_bottom_xyz                                        = center_phe_bottom_xyz./Mean_S1_3D_bottom; %Normalize to center 
    norm_s1_bottom_xyz(isinf(norm_s1_bottom_xyz))             = 1;
    norm_s1_bottom_xyz(isnan(norm_s1_bottom_xyz))             = 1;

    %Get 1 sigma of the corrections matrix.
    sigma_center_phe_bottom_xyz                               = interp3(s1xbins,s1ybins,s1zbins,Sigma_S1_3D_bottom,x_center,y_center,z_center,'cubic');%Normalize to the center 
    sigma_norm_s1_bottom_xyz                                  = sqrt((sigma_center_phe_bottom_xyz./Mean_S1_3D_bottom).^2+(Sigma_S1_3D_bottom.*center_phe_bottom_xyz./Mean_S1_3D_bottom.^2).^2);
    sigma_norm_s1_bottom_xyz(isinf(sigma_norm_s1_bottom_xyz)) = 1;%no infinity. Set Sigma=1, which is 100% uncertainty
    sigma_norm_s1_bottom_xyz(isnan(sigma_norm_s1_bottom_xyz)) = 1;

    mean_center_3D_bottom                                     = interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D_bottom,x_center,y_center,z_center,'cubic');
    sigma_mean_center_3D_bottom                               = interp3(s1xbins,s1ybins,s1zbins,Sigma_S1_3D_bottom,x_center,y_center,z_center,'cubic');

    %to Apply the correction use the following.
    %s1_phe_bottom_xyz=s1_phe_bottom.*interp3(xx,yy,zz,norm_s1_bottom_xyz,s2x,s2y,drift_time,'cubic');

    %% S1 3D correction. Top PMT 
    center_phe_top_xyz                                        = interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D_top,x_center,y_center,z_center,'cubic');%Normalize to the center 
    norm_s1_top_xyz                                           = center_phe_top_xyz./Mean_S1_3D_top; %Normalize to center 
    norm_s1_top_xyz(isinf(norm_s1_top_xyz))                   = 1;
    norm_s1_top_xyz(isnan(norm_s1_top_xyz))                   = 1;

    %Get 1 sigma of the corrections matrix.
    sigma_center_phe_top_xyz                                  = interp3(s1xbins,s1ybins,s1zbins,Sigma_S1_3D_top,x_center,y_center,z_center,'cubic');%Normalize to the center 
    sigma_norm_s1_top_xyz                                     = sqrt((sigma_center_phe_top_xyz./Mean_S1_3D_top).^2+(Sigma_S1_3D_top.*center_phe_top_xyz./Mean_S1_3D_top.^2).^2);
    sigma_norm_s1_top_xyz(isinf(sigma_norm_s1_top_xyz))       = 1;%no infinity. Set Sigma=1, which is 100% uncertainty
    sigma_norm_s1_top_xyz(isnan(sigma_norm_s1_top_xyz))       = 1;

    mean_center_3D_top                                        = interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D_top,x_center,y_center,z_center,'cubic');
    sigma_mean_center_3D_top                                  = interp3(s1xbins,s1ybins,s1zbins,Sigma_S1_3D_top,x_center,y_center,z_center,'cubic');

    %to Apply the correction use the following.
    %s1_phe_top_xyz=s1_phe_top.*interp3(xx,yy,zz,norm_s1_top_xyz,s2x,s2y,drift_time,'cubic');

    %% s1 3d xyz gauss fit

    s1_phe_both_xyz                                           = s1_phe_both.*interp3(s1xbins,s1ybins,s1zbins,norm_s1_both_xyz,s2x,s2y,drift_time,'cubic',1);
    s1_phe_both_xyz3D_gauss_fit_binning                       = (0:10:1000);
    s1_phe_both_xyz3D_gauss_fit_cut                           = inrange(s1_phe_both_xyz,[0 1000]);
    s1_phe_both_xyz3D_gauss_fit                               = fit(s1_phe_both_xyz3D_gauss_fit_binning.',hist(s1_phe_both_xyz(s1_phe_both_xyz3D_gauss_fit_cut),s1_phe_both_xyz3D_gauss_fit_binning).','gauss1');
    s1_phe_both_xyz3D_gauss_fit_mean                          = s1_phe_both_xyz3D_gauss_fit.b1;
    s1_phe_both_xyz3D_gauss_fit_sigma                         = s1_phe_both_xyz3D_gauss_fit.c1/sqrt(2);
    s1_phe_both_xyz3D_gauss_fit_confint                       = confint(s1_phe_both_xyz3D_gauss_fit,0.68);
    s1_phe_both_xyz3D_gauss_fit_mean_error                    = abs(s1_phe_both_xyz3D_gauss_fit_confint(1,2)-s1_phe_both_xyz3D_gauss_fit_mean);
    s1_phe_both_xyz3D_gauss_fit_sigma_error                   = abs(s1_phe_both_xyz3D_gauss_fit_confint(1,3)/sqrt(2) - s1_phe_both_xyz3D_gauss_fit_sigma);

    %% Save The data
    save(strcat('LUX_corrections/',file_id_cp,'/',file_id_cp,'_3D_IQs'),'norm_s1_both_xyz','norm_s1_bottom_xyz','norm_s1_top_xyz', 'sigma_norm_s1_both_xyz'...
    ,'sigma_norm_s1_bottom_xyz','sigma_norm_s1_top_xyz','s1xbins','s1ybins','s1zbins'...
    ,'Count_S1_3D', 'Count_S1_3D_bottom', 'Mean_S1_3D_both', 'Mean_S1_3D_bottom', 'Mean_S1_3D_top', 'Sigma_S1_3D_both', 'Sigma_S1_3D_bottom', 'Sigma_S1_3D_top'...
    ,'mean_center_3D_both','sigma_mean_center_3D_both','mean_center_3D_bottom','sigma_mean_center_3D_bottom','mean_center_3D_top','sigma_mean_center_3D_top','det_edge'...
    ,'s1_phe_both_xyz3D_gauss_fit_mean','s1_phe_both_xyz3D_gauss_fit_sigma','s1_phe_both_xyz3D_gauss_fit_mean_error','s1_phe_both_xyz3D_gauss_fit_sigma_error');

    %% Submit IQ to LUG

    if submit == 1

        user       = user_name;
        source     = 'Kr-83';
        file_id    = file_id_cp(1:19);
        cp         = file_id_cp(23:27);


        % Querying the CP uning the cp numbers from the IQ entries to get the gs numbers and so see if the alternate Chain is being used
         query_str = ['select * from lug_complete_process_record where cp="' cp '" ;' ];
         data2     = MySQLQuery_UMD(query_str,'read'); 
         if str2num(cp) == data2.cp  %#ok<ST2NM>
            gs     = data2.gs;
         end


        version=dp_version;
        clear lug_val Image_paths xmlinput

        lug_val(1) = {file_id};
        lug_val(2) = {cp}; %Proccessing version number
        lug_val(3) = {gs};
        lug_val(4) = {user};
        lug_val(5) = {algorithm_name};
        lug_val(6) = {version};
        lug_val(7) = {Kr_events};              %total count
        lug_val(8) = {s1xbins};                
        lug_val(9) = {s1ybins};                
        lug_val(10) = {s1zbins};
        lug_val(11) = {norm_s1_both_xyz};        %Normalization Matrix
        lug_val(12) = {sigma_norm_s1_both_xyz}; %1 sigma normalization matric
        lug_val(13) = {norm_s1_bottom_xyz};     %Normalization Matrix
        lug_val(14) = {sigma_norm_s1_bottom_xyz};%1 sigma of the normalization matrix
        lug_val(15) = {norm_s1_top_xyz};        %Normalization Matrix
        lug_val(16) = {sigma_norm_s1_top_xyz};  %1 sigma of the normalization matrix
        lug_val(17) = {mean_center_3D_both};        %Phe Mean matrix
        lug_val(18) = {sigma_mean_center_3D_both};
        lug_val(19) = {mean_center_3D_bottom};      
        lug_val(20) = {sigma_mean_center_3D_bottom};
        lug_val(21) = {mean_center_3D_top};      
        lug_val(22) = {sigma_mean_center_3D_top}; 
        lug_val(23) = {Count_S1_3D};    %Count in each bin
        lug_val(24) = {det_edge};
        lug_val(25) = {x_center};
        lug_val(26) = {y_center};
        lug_val(27) = {z_center};
        lug_val(28) = {s1_phe_both_xyz3D_gauss_fit_mean};
        lug_val(29) = {s1_phe_both_xyz3D_gauss_fit_sigma};
        lug_val(30) = {s1_phe_both_xyz3D_gauss_fit_mean_error};
        lug_val(31) = {s1_phe_both_xyz3D_gauss_fit_sigma_error};

        xmlinput.iq.global.filename_prefix = lug_val(1);
        xmlinput.iq.global.cp_number = lug_val(2);
        xmlinput.iq.global.gs_number = lug_val(3);
        xmlinput.iq.global.computed_by = lug_val(4);
        xmlinput.iq.global.computed_date = datestr(now,'yyyymmddTHHMM');
        xmlinput.iq.global.source = source;
        xmlinput.iq.global.source_id = 'unknown';
        xmlinput.iq.global.notes = ' ';
        xmlinput.iq.global.algorithm_name = lug_val(5);
        xmlinput.iq.global.algorithm_version = lug_val(6);
        xmlinput.iq.correction.fit.number_of_kr_events= lug_val(7);
        xmlinput.iq.correction.fit.x_bin_center = lug_val(8); 
        xmlinput.iq.correction.fit.y_bin_center = lug_val(9);
        xmlinput.iq.correction.fit.z_bin_center = lug_val(10);
        xmlinput.iq.correction.fit.norm_s1_both_xyz = lug_val(11);
        xmlinput.iq.correction.fit.sigma_norm_s1_both_xyz = lug_val(12);
        xmlinput.iq.correction.fit.norm_s1_bottom_xyz = lug_val(13);
        xmlinput.iq.correction.fit.sigma_norm_s1_bottom_xyz = lug_val(14);
        xmlinput.iq.correction.fit.norm_s1_top_xyz = lug_val(15);
        xmlinput.iq.correction.fit.sigma_norm_s1_top_xyz = lug_val(16);
        xmlinput.iq.correction.fit.mean_center_both_xyz = lug_val(17);
        xmlinput.iq.correction.fit.sigma_mean_center_both_xyz = lug_val(18);
        xmlinput.iq.correction.fit.mean_center_bottom_xyz = lug_val(19);
        xmlinput.iq.correction.fit.sigma_mean_center_bottom_xyz = lug_val(20);
        xmlinput.iq.correction.fit.mean_center_top_xyz = lug_val(21);
        xmlinput.iq.correction.fit.sigma_mean_center_top_xyz = lug_val(22);
        xmlinput.iq.correction.fit.count_s1_both_xyz = lug_val(23);
        xmlinput.iq.correction.fit.detector_edge = lug_val(24);
        xmlinput.iq.correction.fit.x_center = lug_val(25);
        xmlinput.iq.correction.fit.y_center = lug_val(26);
        xmlinput.iq.correction.fit.z_center = lug_val(27);
        xmlinput.iq.correction.fit.s1_phe_both_xyz3D_gauss_fit_mean = lug_val(28);
        xmlinput.iq.correction.fit.s1_phe_both_xyz3D_gauss_fit_sigma = lug_val(29);
        xmlinput.iq.correction.fit.s1_phe_both_xyz3D_gauss_fit_mean_error = lug_val(30);
        xmlinput.iq.correction.fit.s1_phe_both_xyz3D_gauss_fit_sigma_error = lug_val(31);

        LUXSubmitIQ(user,'s1_xyz_correction',algorithm_name,version,file_id,cp,xmlinput,'Kr83 Calibration','','');

        clear lug_val Image_paths xmlinput

    end   

    fprintf('KrypCal XYZ Finished \n');
    
end