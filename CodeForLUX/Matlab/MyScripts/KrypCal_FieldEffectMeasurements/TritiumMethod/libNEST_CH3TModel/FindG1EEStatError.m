load('S1aS1bParamScan_EnergyChi2_S1FloatPoly_S2Poly2015_ImprovedFineRes.mat')


%% Find EE array with everything else held constant
ee_cut=total_chi2_oldsim_list>0 & a2_list==besta2_oldsim & b2_list==bestb2_oldsim &...
    g1_list==bestg1_oldsim;


% Find 1.4 chi2 diff = 0.68 interval = 1 sigma interval
EErange=ee_list(ee_cut);
Chi2range=total_chi2_oldsim_list(ee_cut);

EEupperbound=max(EErange(Chi2range-min(Chi2range) < 1.4));
EE_one_sig=EEupperbound-bestEE_oldsim


%% Repeat for g1
g1_cut=total_chi2_oldsim_list>0 & a2_list==besta2_oldsim & b2_list==bestb2_oldsim &...
    ee_list==bestEE_oldsim;

% Find 1.4 chi2 diff = 0.68 interval = 1 sigma interval
g1range=g1_list(g1_cut);
Chi2range=total_chi2_oldsim_list(g1_cut);

g1upperbound=max(g1range(Chi2range-min(Chi2range) < 1.4));
g1_one_sig=g1upperbound-bestg1_oldsim
