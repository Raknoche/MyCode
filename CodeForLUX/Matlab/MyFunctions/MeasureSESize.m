%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This function fits a skew Gaussian to a single electron pulse area histogram
%
%  Inputs:
%       - SE_hist_binning   - Binning to use in the histogram
%       - SE_pulse_area     - Single electron pulse areas after golden event selection
%       - SE_radius         - Single electron radius after golden event selection
%
%  Outputs:
%       - SE_mu             - Single electron mean
%       - SE_sig_mu         - Error on single electron mean
%       - SE_sig            - Single electron 1-sigma width
%       - SE_sig_sig        - Error on 1-sigma width
%       - skew              - Signle electron skew
%       - SE_cut            - SE pulse areas after applying radial cut
%       - xfit              - X data for skew gaussian fit
%       - yfit              - Y data for skew gaussian fit
%
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ SE_mu, SE_sig_mu, SE_sig, SE_sig_sig, skew, SE_cut, xfit, yfit ] = MeasureSESize( SE_hist_binning, SE_pulse_area, SE_radius)
 
    %Intial fit with Gaussian distribution
    SE_max_radius      = 17; %Radial cut from study on SE size change v radius
    SE_fit_cut         = SE_radius<SE_max_radius & inrange(SE_pulse_area,[min(SE_hist_binning) max(SE_hist_binning)]);
    SE_fit             = fit(SE_hist_binning',hist(SE_pulse_area(SE_fit_cut),SE_hist_binning)','gauss1');
    SE_Size            = SE_fit.b1;
    SE_sig             = SE_fit.c1/sqrt(2);
    SE_conf            = confint(SE_fit,0.683);
    SE_sig_mu          = abs(SE_Size-SE_conf(1,2)); % one sigma
    SE_sig_sig         = abs(SE_sig-SE_conf(1,3)/sqrt(2)); % one sigma

    %Reffiting with a skew Gaussian distribution
    n_sig_fit          = 2;
    fit_start          = SE_Size-n_sig_fit*SE_sig;
    fit_end            = SE_Size+n_sig_fit*SE_sig;
    xfit               = fit_start:0.1:fit_end;
    SE_cut             = SE_pulse_area(SE_radius<SE_max_radius);
    [xxo, yyo, xo, no] = step(SE_cut(inrange(SE_cut,[fit_start fit_end])),SE_hist_binning);

    s                  = fitoptions('Method','NonlinearLeastSquares','Robust','On','Algorithm','Trust-Region','Lower', ...
                            [-5,0,0,0],'Upper',[5,Inf,Inf,Inf],'Startpoint',[1,max(no),10,5]);
    skewGauss1         = fittype('a1/c1*exp(-((x-b1)/c1)^2)*(1+erf(a*(x-b1)/c1))','options',s);

    [f,g]              = fit(xo',no',skewGauss1);
    delta              = f.a/sqrt(1+f.a^2);
    SE_mu              = f.b1 + f.c1/sqrt(2)*delta*sqrt(2/pi);
    SE_sig             = sqrt((f.c1/sqrt(2))^2 * (1 - 2 * (delta)^2/pi));
    skew               = (4-pi)/2 * (delta * sqrt(2/pi))^3 / (1 - 2 * delta^2/pi)^(3/2);
    kurt               = 2*(pi-3) * (delta * sqrt(2/pi))^3 / (1 - 2 * delta^2/pi)^2;
    ci                 = confint(f,0.683);
    yfit               = f(xfit);   

end

