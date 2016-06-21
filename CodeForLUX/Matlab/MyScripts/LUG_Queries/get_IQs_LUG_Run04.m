%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This script queries the LUG and produces a number of plots which track S1, S2, and SE
%       information over the course of Run04. MUST BE LOGGED INTO SANFORD VPN TO RUN THIS
% 
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get electron lifetime and 1-sigma from lug

%get all the IQs that contain the term 'electron_lifetime' and store them
out = MySQLQuery('select values_xml from control.lug_iqs where iq_type = "electron_lifetime" and iq > 8934 and strikeme = 0 and algorithm_version = 2.22;');

for j=1:size(out.values_xml,1); % 22(first entry) to last entry

    %read in electron lifetime and 1 sigma, using the XMLParser function
    value1 = XMLParser(out.values_xml{j});
    lifetime(j)               = value1.iq.correction.fit.electron_attenuation_us; %#ok<*SAGROW>
    sigma_lifetime(j)         = value1.iq.correction.fit.sigma_electron_attenuation_us;
    num_events(j)             = value1.iq.correction.fit.number_of_kr_events;
    pseduo_lifetime(j)        = value1.iq.correction.fit.pseudo_lifetime_both;
    corrected_s2_at_top(j)    = value1.iq.correction.fit.s2_both_mean_at_top_z;
    corrected_s2_at_top_err(j)= value1.iq.correction.fit.s2_both_mean_at_top_err_z;

    %get the date from file name prexix (as Matlab time)
    file_date       = value1.iq.global.filename_prefix;
    Date_MatTime(j) = datenum(file_date(7:19),'yyyymmddTHHMM');

end

%% Make a  plot of the electron lifetime, including circulation outages
%Manually putting in low purity data sets which KrypCal didn't process
manual_lifetime       = [92.16,99.5, 113.15, 124.3, 131.95, 136.2, 48.75, 59.21, 75.88, 91.67, 102.8,...
                        135.3,177.1, 228.5, 232];
manual_sigma_lifetime = [2.3,1.4, 0.7, 0.8, 1.1, 0.95, 1.5,1.6, 1.41, 1.5, 3.9,...
                        1.2,2.28, 2.84, 1.19];
manual_Date_MatTime   = [datenum('20150406T1558','yyyymmddTHHMM'),datenum('20150415T1902','yyyymmddTHHMM'),...
                        datenum('20150420T1904','yyyymmddTHHMM'),datenum('20150423T1723','yyyymmddTHHMM'),...
                        datenum('20150426T1216','yyyymmddTHHMM'),datenum('20150429T1606','yyyymmddTHHMM'),datenum('20150503T1043','yyyymmddTHHMM'),...
                        datenum('20150505T1934','yyyymmddTHHMM'),datenum('20150511T2036','yyyymmddTHHMM'),...
                        datenum('20150514T2141','yyyymmddTHHMM'),datenum('20150517T1053','yyyymmddTHHMM'),...
                        datenum('20150520T1940','yyyymmddTHHMM'),datenum('20150523T1945','yyyymmddTHHMM'),...
                        datenum('20150526T1747','yyyymmddTHHMM'),datenum('20150527T1032','yyyymmddTHHMM')];

Date_MatTime_patched    = [Date_MatTime,manual_Date_MatTime];
lifetime_patched        = [lifetime,manual_lifetime];
sigma_lifetime_patched  = [sigma_lifetime,manual_sigma_lifetime];
pseudo_lifetime_patched = [pseduo_lifetime,manual_lifetime];    


% Get outage data from SC
if ~exist('xml_settings','var')    
    xmlsettings = XMLReader('defaultSCSettings.xml'); %using lux.sanfordlab.org as default
end


%Get MFC1 Settings -- have to do this in steps, otherwise we run out of memory
%Could do this in a loop.. but copy/paste was easier =) (This code won't change in the future, so not a big deal)
disp('Starting MFC1  Data Collecting')
start_time   = (datenum(2014,09,01,0,0,0)-datenum(1970,1,1))*86400;
end_time     =(datenum(2014,10,01,0,0,0)-datenum(1970,1,1))*86400;
querystr     = sprintf('select * from control.sc_sens_ACRS_Flow where time > %d and time < %d order by time asc',start_time,end_time);
query_result = MySQLQuery(querystr);
MFC1_time    = query_result.time;
MFC1_value   = query_result.value;
sprintf('Finished Month of %d',09)
clear query_result
for i=10:11;
start_time   = (datenum(2014,i,01,0,0,0)-datenum(1970,1,1))*86400;
end_time     = (datenum(2014,i+1,01,0,0,0)-datenum(1970,1,1))*86400;
querystr     = sprintf('select * from control.sc_sens_ACRS_Flow where time > %d and time < %d order by time asc',start_time,end_time);
query_result = MySQLQuery(querystr);
MFC1_time    = [MFC1_time;query_result.time]; 
MFC1_value   = [MFC1_value;query_result.value]; %#ok<*AGROW>
sprintf('Finished Month of %d',i)
clear query_result
end
start_time   = (datenum(2014,12,01,0,0,0)-datenum(1970,1,1))*86400;
end_time     = (datenum(2015,1,01,0,0,0)-datenum(1970,1,1))*86400;
querystr     = sprintf('select * from control.sc_sens_ACRS_Flow where time > %d and time < %d order by time asc',start_time,end_time);
query_result = MySQLQuery(querystr);
MFC1_time    = [MFC1_time;query_result.time];
MFC1_value   = [MFC1_value;query_result.value];
sprintf('Finished Month of %d',12)
clear query_result
for i=1:12;
start_time   = (datenum(2015,i,01,0,0,0)-datenum(1970,1,1))*86400;
end_time     = (datenum(2015,i+1,01,0,0,0)-datenum(1970,1,1))*86400;
querystr     = sprintf('select * from control.sc_sens_ACRS_Flow where time > %d and time < %d order by time asc',start_time,end_time);
query_result = MySQLQuery(querystr);
MFC1_time    = [MFC1_time;query_result.time];
MFC1_value   = [MFC1_value;query_result.value];
sprintf('Finished Month of %d',i)
clear query_result
end
clear query_result
for i=1:6;
start_time   = (datenum(2016,i,01,0,0,0)-datenum(1970,1,1))*86400;
end_time     = (datenum(2016,i+1,01,0,0,0)-datenum(1970,1,1))*86400;
querystr     = sprintf('select * from control.sc_sens_ACRS_Flow where time > %d and time < %d order by time asc',start_time,end_time);
query_result = MySQLQuery(querystr);
MFC1_time    = [MFC1_time;query_result.time];
MFC1_value   = [MFC1_value;query_result.value];
sprintf('Finished Month of %d',i)
clear query_result
end

MFC1_date    = (MFC1_time./86400 + datenum(1970,1,1));
noflow_times = (MFC1_value<15);
 
 
%Plot lifetime with flow outages shaded in
figure
hold on;

%Loop to turn fill on/off during circulation outages
noflow=1;
first_noflow=1;
last_noflow=0;
for scan=1:length(noflow_times);
    if noflow_times(scan)==1; %if we hit a one, record the start date
        if first_noflow==1;
            start_noflow_date= MFC1_date(scan);
            first_noflow     = 0; %got the start date, scan through rest without carring about date
            last_noflow      = 1; %the next zero will be recorded as the end date
        end
    else
        if last_noflow==1;
            end_noflow_date  = MFC1_date(scan-1); %zero-1 is the end date
            last_noflow      = 0;
            first_noflow     = 1; %next one will be recored as the start date
            h                = fill([start_noflow_date start_noflow_date end_noflow_date end_noflow_date], [0 4000 4000 0],[0.8 0.5 0.5]);%FILL IN THE COLOR FROM START TO END;
            set(h,'EdgeColor','None');
        end
    end
end
errorbar(Date_MatTime_patched([1:32,34:end]),lifetime_patched([1:32,34:end]),sigma_lifetime_patched([1:32,34:end]),'.k')
errorbar(Date_MatTime_patched([1:32,34:end]),pseudo_lifetime_patched([1:32,34:end]),sigma_lifetime_patched([1:32,34:end]),'.b')
datetick %formats the x axis to be a date label
ylim([0 3000]);
myfigview(16); xlabel('Date'); ylabel('Electron Lifetime (\mus)');
 

%Make a plot of Z and Field Corrected S2 over time
figure
hold on;
h=fill([735781 735781 736500 736500],[mean(corrected_s2_at_top)-std(corrected_s2_at_top) mean(corrected_s2_at_top)+std(corrected_s2_at_top) mean(corrected_s2_at_top)+std(corrected_s2_at_top) mean(corrected_s2_at_top)-std(corrected_s2_at_top)],[0.8 0.8 0.8]);
set(h,'EdgeColor','None');
plot([735781 736500],[mean(corrected_s2_at_top) mean(corrected_s2_at_top)],'--b')
errorbar(Date_MatTime([1:32,34:end]),corrected_s2_at_top([1:32,34:end]),corrected_s2_at_top_err([1:32,34:end]),'.k');
xlabel('Date'); ylabel('Z and Field Corrected Kr S2 Peak (phe)');
datetick('x',3);  
xlim([735831 736500]);ylim([0 25000]);
myfigview(16);

%% Make a plot of S2 XYZ corrected peaks

%get all the IQs that contain the term 'electron_lifetime' and store them
clear Date_MatTime
out = MySQLQuery('select values_xml from control.lug_iqs where iq_type = "s2_xy_correction" and iq > 8934 and strikeme = 0 and algorithm_version = 2.22;');

for j=1:size(out.values_xml,1); % 22(first entry) to last entry

%read in electron lifetime and 1 sigma, using the XMLParser function
value1                         = XMLParser(out.values_xml{j});

xyz_corrected_s2_at_top(j)     = value1.iq.correction.fit.s2_both_mean_at_top_xyz;
xyz_corrected_s2_at_top_err(j) = value1.iq.correction.fit.s2_both_mean_at_top_err_xyz;

%get the date from file name prexix (as Matlab time)
file_date                      = value1.iq.global.filename_prefix;
Date_MatTime(j)                = datenum(file_date(7:19),'yyyymmddTHHMM');

end

figure
hold on;
h=fill([735781 735781 736500 736500],[mean(xyz_corrected_s2_at_top)-std(xyz_corrected_s2_at_top) mean(xyz_corrected_s2_at_top)+std(xyz_corrected_s2_at_top) mean(xyz_corrected_s2_at_top)+std(xyz_corrected_s2_at_top) mean(xyz_corrected_s2_at_top)-std(xyz_corrected_s2_at_top)],[0.8 0.8 0.8]);
set(h,'EdgeColor','None');
plot([735781 736500],[mean(xyz_corrected_s2_at_top) mean(xyz_corrected_s2_at_top)],'--b')
errorbar(Date_MatTime([1:32,34:end]),xyz_corrected_s2_at_top([1:32,34:end]),xyz_corrected_s2_at_top_err([1:32,34:end]),'.k');
xlabel('Date');ylabel('Field and KrypCal XYZ Corrected S2 (phe)');
datetick('x',3);  
xlim([735831 736500]);ylim([0 25000]);
myfigview(16);

%% Make a plot of S1 Z correction IQs
out = MySQLQuery('select values_xml from control.lug_iqs where iq_type = "z_dep_s1_correction" and iq > 8934 and strikeme = 0 and algorithm_version = 2.22;');
clear Date_MatTime j value1

for j=1:size(out.values_xml,1); % 22(first entry) to last entry

%read in electron lifetime and 1 sigma, using the XMLParser function
value1                        = XMLParser(out.values_xml{j});

%technically this is not a fit to a gaussian, it is the value of S1 at the center of the detector
%from the Z dependence plot. S2 above is found from a gaussian fit though
corrected_s1_at_center(j)     = value1.iq.correction.fit.center_s1_both_z;
corrected_s1_at_center_err(j) = value1.iq.correction.fit.center_s1_both_err_z;

%get the date from file name prexix (as Matlab time)
file_date                     = value1.iq.global.filename_prefix;
Date_MatTime(j)               = datenum(file_date(7:19),'yyyymmddTHHMM');

end


%Make a plot of Z and Field Corrected S1 over time
figure
hold on;
h=fill([735781 735781 736500 736500],[mean(corrected_s1_at_center)-std(corrected_s1_at_center) mean(corrected_s1_at_center)+std(corrected_s1_at_center) mean(corrected_s1_at_center)+std(corrected_s1_at_center) mean(corrected_s1_at_center)-std(corrected_s1_at_center)],[0.8 0.8 0.8]);
set(h,'EdgeColor','None');
plot([735781 736500],[mean(corrected_s1_at_center) mean(corrected_s1_at_center)],'--b')
errorbar(Date_MatTime([1:32,34:end]),corrected_s1_at_center([1:32,34:end]),corrected_s1_at_center_err([1:32,34:end]),'.k');
xlabel('Date'); ylabel('Z and Field Corrected Kr S1 Peak (phe)');
datetick('x',3);  
xlim([735831 736500]);ylim([150 300]);
myfigview(16);

%% Make a plot of S1 XYZ correction IQs
clear Date_MatTime
out = MySQLQuery('select values_xml from control.lug_iqs where iq_type = "s1_xy_correction" and iq > 8934 and strikeme = 0 and algorithm_version = 2.22;');

for j=1:size(out.values_xml,1); % 22(first entry) to last entry

%read in electron lifetime and 1 sigma, using the XMLParser function
value1                            = XMLParser(out.values_xml{j});

xyz_corrected_s1_at_center(j)     = value1.iq.correction.fit.center_s1_both_z; %really this is xyz, bad labeling in KrypCal
xyz_corrected_s1_at_center_err(j) = value1.iq.correction.fit.center_s1_both_err_z;

%get the date from file name prexix (as Matlab time)
file_date                         = value1.iq.global.filename_prefix;
Date_MatTime(j)                   = datenum(file_date(7:19),'yyyymmddTHHMM');

end

figure
hold on;
h=fill([735781 735781 736500 736500],[mean(xyz_corrected_s1_at_center)-std(xyz_corrected_s1_at_center) mean(xyz_corrected_s1_at_center)+std(xyz_corrected_s1_at_center) mean(xyz_corrected_s1_at_center)+std(xyz_corrected_s1_at_center) mean(xyz_corrected_s1_at_center)-std(xyz_corrected_s1_at_center)],[0.8 0.8 0.8]);
set(h,'EdgeColor','None');
plot([735781 736500],[mean(xyz_corrected_s1_at_center) mean(xyz_corrected_s1_at_center)],'--b')
errorbar(Date_MatTime([1:32,34:end]),xyz_corrected_s1_at_center([1:32,34:end]),xyz_corrected_s1_at_center_err([1:32,34:end]),'.k');
xlabel('Date');ylabel('Field and KrypCal XYZ Corrected S1 (phe)');
datetick('x',3);  
xlim([735831 736500]);ylim([150 300]);
myfigview(16);

%% Make a plot of SE Size IQs

out = MySQLQuery('select values_xml from control.lug_iqs where iq_type = "single_e_kr" and strikeme = 0 and algorithm_version = 2.22;');

for j=1:size(out.values_xml,1); % 22(first entry) to last entry

    %read in electron lifetime and 1 sigma, using the XMLParser function
    value1 = XMLParser(out.values_xml{j});

    single_e_both(j)               = value1.iq.correction.fit.e_mean_both;
    error_single_e_both(j)         = value1.iq.correction.fit.e_mean_err_both;
    single_e_bottom(j)             = value1.iq.correction.fit.e_mean_bottom;
    error_single_e_bottom(j)       = value1.iq.correction.fit.e_mean_err_bottom;
    single_e_top(j)                = value1.iq.correction.fit.e_mean_top;
    error_single_e_top(j)          = value1.iq.correction.fit.e_mean_err_top;

    sigma_single_e_both(j)         = value1.iq.correction.fit.e_sig_both;
    error_sigma_single_e_both(j)   = value1.iq.correction.fit.e_sig_err_both;
    sigma_single_e_bottom(j)       = value1.iq.correction.fit.e_sig_bottom;
    error_sigma_single_e_bottom(j) = value1.iq.correction.fit.e_sig_err_bottom;
    sigma_single_e_top(j)          = value1.iq.correction.fit.e_sig_top;
    error_sigma_single_e_top(j)    = value1.iq.correction.fit.e_sig_err_top;

    efit_skew_both(j)              = value1.iq.correction.fit.e_fit_skew_both;
    efit_skew_bottom(j)            = value1.iq.correction.fit.e_fit_skew_bottom;
    efit_skew_top(j)               = value1.iq.correction.fit.e_fit_skew_top;


    %get the date from file name prexix (as Matlab time)
    file_date                      = value1.iq.global.filename_prefix;
    readable_date{j}               = file_date;
    Date_MatTime(j)                = datenum(file_date(7:19),'yyyymmddTHHMM');

end


figure
hold on;
h=fill([735781 735781 736500 736500],[mean(single_e_both)-std(single_e_both) mean(single_e_both)+std(single_e_both) mean(single_e_both)+std(single_e_both) mean(single_e_both)-std(single_e_both)],[0.8 0.8 0.8]);
set(h,'EdgeColor','None');
plot([735781 736500],[mean(single_e_both) mean(single_e_both)],'--b')
errorbar(Date_MatTime([1:end]),single_e_both([1:end]),error_single_e_both([1:end]),'.k'); %#ok<*NBRAK>
xlabel('Date');ylabel('SE Size (phe)');
datetick('x',3);  
xlim([735831 736500]);ylim([23 29]);
myfigview(16);
