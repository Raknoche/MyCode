%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This function produces single electron selection cuts from golden S1 and S2 events
%
%  Inputs:
%       - d                  - the data structure to be analysed.
%       - s1_single_cut      - 10xN Golden S1 selection cut
%       - s2_single_cut      - 10xN Golden S2 selection cut
%
%  Outputs:
%       - SE_good_cut      - Golden event cut for S1 pulses
%
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ SE_good_cut ] = GetSECuts( d, s1_single_cut, s2_single_cut)
 
    %Pulse classification cut
    s4_class=(d.pulse_classification==4); 
    
    %Cut to remove S1 tails that improperly classified as SE
    d = se_wrt_first_s1_func(d);
    cut_se_s1_z = d.t_btw_se_first_s1 > 100*10; %Cut is based on timing after S1

    %Final cut must be between the Kr S1 and S2, and meet the cuts above
    SE_good_cut =logical( cut_se_s1_z & s4_class & ( cumsum(s1_single_cut)-cumsum(s2_single_cut) ) );

end

