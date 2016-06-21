%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This function produces XY dependence maps of the single electron size
%
%  Inputs:
%       - SE_pulse_area        - Single electron pulse areas after golden selection
%       - SE_x                 - Single electron x position
%       - SE_y                 - Single electron y position
%       - SE_hist_binning      - Pulse area binning to be used in the Guassian fits to each XY bin
%
%  Outputs:
%       - se_xbins             - 1xN vector defining x binning in XY map
%       - se_ybins             - 1xM vector defining y binning in XY map
%       - mean_SE_xy           - MxN matrix of SE means
%       - skew_SE_xy           - MxN matrix of SE skew
%       - sigma_SE_xy          - MxN matrix of SE width (1-sigma)
%       - kurt_SE_xy           - MxN matrix of SE kurtosis
%       - mean_err_SE_xy       - MxN matrix of the error on SE mean measurement
%       - sigma_err_SE_xy      - MxN matrix of the error on SE width measurement
% 
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [se_xbins, se_ybins, mean_SE_xy, skew_SE_xy, sigma_SE_xy, kurt_SE_xy, mean_err_SE_xy,...
    sigma_err_SE_xy] = MeasureSEXYDep(SE_pulse_area, SE_x, SE_y, SE_hist_binning)
 
    %Define map edges
    SE_xbin_min     = -25;
    SE_xbin_max     = 25;
    SE_ybin_min     = -25;
    SE_ybin_max     = 25;

    %Define bin size
    SExybinsize     = sqrt((50*50)/(length(SE_pulse_area)/100));

    %Define bins
    se_xbins        = SE_xbin_min+SExybinsize/2:SExybinsize:SE_xbin_max;
    se_ybins        = SE_ybin_min+SExybinsize/2:SExybinsize:SE_ybin_max;

    %Preallocate memory
    mean_SE_xy      = zeros(length(se_xbins),length(se_ybins));
    mean_err_SE_xy  = zeros(length(se_xbins),length(se_ybins));
    skew_SE_xy      = zeros(length(se_xbins),length(se_ybins));
    sigma_SE_xy     = zeros(length(se_xbins),length(se_ybins));
    kurt_SE_xy      = zeros(length(se_xbins),length(se_ybins));
    sigma_err_SE_xy = zeros(length(se_xbins),length(se_ybins));
   
    %Make SE XY map
    for x_bin = SE_xbin_min:SExybinsize:(SE_xbin_max-SExybinsize);
        for y_bin = SE_ybin_min:SExybinsize:(SE_ybin_max-SExybinsize);          
            x_min   = x_bin; x_max = x_bin+SExybinsize;
            y_min   = y_bin; y_max = y_bin+SExybinsize;      
      
            %Matrix indicies
            x_count = int32(1+x_bin/SExybinsize+SE_xbin_max/SExybinsize);
            y_count = int32(1+y_bin/SExybinsize+SE_ybin_max/SExybinsize); 
      
            %Cut to select events in the XY bin
            bin_cut = (SE_pulse_area>0) & inrange(SE_x,[x_min,x_max]) & inrange(SE_y,[y_min,y_max]);            

            %If there are >60 events, calculate the map entry. Otherwise fill in with zero
            if length(SE_pulse_area(bin_cut)) > 60; 
                
                %Initial fit with a Guassian distribution
                clear temp_SE
                temp_SE            = SE_pulse_area(bin_cut);
                SE_fit_cut         = inrange(temp_SE,[min(SE_hist_binning) max(SE_hist_binning)]);
                SE_fit_xy          = fit(SE_hist_binning',hist(temp_SE(SE_fit_cut),SE_hist_binning)','gauss1');
                SE_Size_xy         = SE_fit_xy.b1;
                SE_sig_xy          = SE_fit_xy.c1/sqrt(2);
                SE_conf_xy         = confint(SE_fit_xy,0.683);
                SE_sig_mu_xy       = abs(SE_Size_xy-SE_conf_xy(1,2)); % one sigma
                SE_sig_sig_xy      = abs(SE_sig_xy-SE_conf_xy(1,3)/sqrt(2)); % one sigma

                %Refit with a skew Gaussian distribution
                n_sig_fit          = 2;
                fit_start          = SE_Size_xy-n_sig_fit*SE_sig_xy;
                fit_end            = SE_Size_xy+n_sig_fit*SE_sig_xy;
                xfit               = fit_start:0.1:fit_end;
                [xxo, yyo, xo, no] = step(temp_SE(inrange(temp_SE,[fit_start fit_end])),SE_hist_binning); %#ok<ASGLU>
                s                  = fitoptions('Method','NonlinearLeastSquares','Robust','On','Algorithm','Trust-Region','Lower',...
                                        [-5,0,0,0],'Upper',[5,Inf,Inf,Inf],'Startpoint',[1,max(no),10,5]);
                skewGauss1         = fittype('a1/c1*exp(-((x-b1)/c1)^2)*(1+erf(a*(x-b1)/c1))','options',s);
                [f_xy,g_xy]        = fit(xo',no',skewGauss1); %#ok<NASGU>
                delta_xy           = f_xy.a/sqrt(1+f_xy.a^2);
                SE_mu_xy           = f_xy.b1 + f_xy.c1/sqrt(2)*delta_xy*sqrt(2/pi);
                SE_sig_xy          = sqrt((f_xy.c1/sqrt(2))^2 * (1 - 2 * (delta_xy)^2/pi));
                skew_xy            = (4-pi)/2 * (delta_xy * sqrt(2/pi))^3 / (1 - 2 * delta_xy^2/pi)^(3/2);
                kurt_xy            = 2*(pi-3) * (delta_xy * sqrt(2/pi))^3 / (1 - 2 * delta_xy^2/pi)^2; 

                %Fill in map (xy indicies swapped intentional due to how interp2D works)
                mean_SE_xy(y_count,x_count)      = SE_mu_xy;  
                skew_SE_xy(y_count,x_count)      = skew_xy; 
                sigma_SE_xy(y_count,x_count)     = SE_sig_xy;
                kurt_SE_xy(y_count,x_count)      = kurt_xy;
                mean_err_SE_xy(y_count,x_count)  = SE_sig_mu_xy;
                sigma_err_SE_xy(y_count,x_count) = SE_sig_sig_xy;    


              else

                mean_SE_xy(y_count,x_count)      = 0;  
                skew_SE_xy(y_count,x_count)      = 0; 
                sigma_SE_xy(y_count,x_count)     = 0;
                kurt_SE_xy(y_count,x_count)      = 0;
                mean_err_SE_xy(y_count,x_count)  = 0;
                sigma_err_SE_xy(y_count,x_count) = 0;    

            end
        end
    end

     
end

