%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This function produces golden event selection cuts between S1 and S2 area bounds
%
%  Inputs:
%       - d                  - the data structure to be analysed.
%       - s1area_bound_min   - Minimum S1 size
%       - s1area_bound_max   - Maximum S1 size
%       - s2area_bound_min   - Minimum S2 size
%       - s2area_bound_max   - Maximum S2 size
%
%  Outputs:
%       - s1_single_cut      -Golden event cut for S1 pulses
%       - s2_single_cut      - Golden event cut for S2 pulses
%
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ s1_single_cut, s2_single_cut ] = GetGoldenCuts( d, s1area_bound_min, s1area_bound_max, s2area_bound_min, s2area_bound_max )

%Pulse area cuts
cut_s2_area = inrange(d.pulse_area_phe, [s2area_bound_min, s2area_bound_max]);
cut_s1_area = inrange(d.pulse_area_phe, [s1area_bound_min, s1area_bound_max]);

%Pulse classification cuts (S2 cut is defined above 100 phe, and paired with an s1)
cut_pulse_s1 = d.pulse_classification == 1 | d.pulse_classification == 9;
cut_pulse_s2 = d.pulse_classification == 2;
cut_s2_with_threshold = d.pulse_area_phe.*cut_pulse_s2 > 100; % subset of cut_pulse_s2
cut_legit_s2_in_legit_event = d.s1s2_pairing.*cut_s2_with_threshold; % this should be used as s2 classification cuts

%Golden cut defines golden events to be events which have one and only one paired S2 above the threshold of 100 phe
%Note there can be multiple S1s in the golden event
cut_golden_event = sum(cut_legit_s2_in_legit_event) == 1; 
cut_s2_in_golden_events = logical(repmat(cut_golden_event,[10,1]).*cut_legit_s2_in_legit_event); %Selects S2 that is in a golden event
cut_s1_in_golden_events = logical(repmat(cut_golden_event,[10,1]).*cut_pulse_s1.*d.s1s2_pairing); %Selects first S1 that is in a golden event

%Combine area cuts and golden cuts
cut_s2_for = cut_s2_in_golden_events.*cut_s2_area; %Selects S2 that is in a golden event and in Kr area bounds
cut_s1_for = cut_s1_in_golden_events.*cut_s1_area; %Selects first S1 that is in a golden event and in Kr area bounds

%Require that "good" golden events have only one S1, that the S1 be within area bounds, and the S2 be within area bounds
%Note sum(cut_s1_for) == 1 above only requires that the first of the S1 in an event be within area bounds, since the S1S2pairing part of cut_s1_in_golden_events is 0 for all subsequent S1s in the events
cut_selected_events = sum(cut_s2_for) == 1 & sum(cut_s1_for) == 1 & sum(d.pulse_classification==1)==1; 

%Final event selection cuts
s1_single_cut = logical(repmat(cut_selected_events,[10,1]).*cut_s1_in_golden_events);
s2_single_cut = logical(repmat(cut_selected_events,[10,1]).*cut_s2_in_golden_events);

end

