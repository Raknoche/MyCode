%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This function produces golden event selection cuts for Kr S1a and S1b pulses 
%
%  Inputs:
%       - d                  - the data structure to be analysed.
%       - s1_single_cut      - Golden event cut for merged Kr S1 pulses
%       - s2_single_cut      - Golden event cut for merged Kr S2 pulses
%
%  Outputs:
%       - s1ab_cut           - Golden event cut for S1a and S1b pulses
%
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ s1ab_cut ] = GetS1aS1bCuts(d, s1_single_cut, s2_single_cut)

%Make timing variables
s1ab_timing                        = d.Kr83fit_dt_samples;

%Set up cuts
s1ab_golden_cut                    = logical(sum(s1_single_cut(:,:),1));
s1ab_pulsearea_cut                 = d.Kr83fit_s1b_area_phe > 30;
s1ab_timing_cut                    = s1ab_timing > 13;  
s1ab_radius_cut_temp               = (sqrt(d.x_cm.^2 + d.y_cm.^2) < 25) & s2_single_cut;
s1ab_radius_cut                    = logical(sum(s1ab_radius_cut_temp(:,:),1));
s1ab_cut                           = s1ab_golden_cut & s1ab_timing_cut & s1ab_pulsearea_cut & s1ab_radius_cut;

end

