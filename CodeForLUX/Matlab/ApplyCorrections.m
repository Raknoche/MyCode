%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This function queries the LUG to find and apply the appropriate XY, and Z corrections 
%       for the dataset being processed.  This runs as part of the DP framework, and produces
%       a number of correction RQs that are added to the dataset.
%
%  Inputs:
%       - filename_evt               - String containing name of the current .evt dataset
%       - data_path_evt              - String containing path to evt file
%       - filename_rq                - String containing name of the current .rq dataset
%       - data_path_rq               - String containing path to the dataset
%       - data_processing_xml_path   - String containing path to DP xml
%       - iq_xml_path                - String containing path to IQ xml
%
%  Outputs:
%       - status                     - Flag indicating successful processing
%
%  Example usage: status = ApplyCorrections_Basic(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
% 
%  Author:
%       - Attila Dobi, Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function status = ApplyCorrections(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path) %#ok<*INUSL>

%% Load .rq file

if strcmp(filename_rq(1),'.') % check for '.' infornt of the file name
    filename_prefix = filename_rq(2:20);
else
    filename_prefix = filename_rq(1:19);
end

dp                       = LUXLoadRQ1s_framework(filename_rq, data_path_rq);

settings.evt_settings    = dp.admin.evt_settings;
settings.daq_settings    = dp.admin.daq_settings;
settings.filename_prefix = dp.admin.filename_prefix;

%% Bookkeeping

myname            = 'ApplyCorrections';

fprintf('\n\n *** Starting module %s\n',myname);
fprintf('\n Data processing xml file = %s\n',data_processing_xml_path);

dp_settings_xml   = XMLReader_framework(data_processing_xml_path);
lug_iqs_xml       = XMLReader_framework(iq_xml_path);
module_names      = {dp_settings_xml.data_processing_settings.module.module_name};
index_temp        = strfind(module_names,myname);
index_module      = find(not(cellfun('isempty', index_temp)));
mymodule_settings = dp_settings_xml.data_processing_settings.module(index_module).parameters; %#ok<*FNDSB>
max_num_pulses    = dp_settings_xml.data_processing_settings.global.max_num_pulses;

%% Initialize variables

pmt_chs          = 1:122; % must take from daq settings
N                = length(dp.event_number);
pulse_event_size = [max_num_pulses N];
livetime         = dp.admin.livetime;

% Correcting issue with top_bottom_ratio
[a1 b1]          = size(dp.pulse_area_phe);  
[a2 b2]          = size(dp.top_bottom_ratio);  adif = a1-a2;

blanks           = [];
for i=1:b1
   blanks = [blanks 0]; 
end

if adif>0
   for i=1:adif
     dp.top_bottom_ratio = vertcat(dp.top_bottom_ratio,blanks);    
   end
end

% Initializing corrected values of area
dp.z_corrected_pulse_area_all_phe   = dp.pulse_area_phe;
dp.xyz_corrected_pulse_area_all_phe = dp.pulse_area_phe;
dp.z_corrected_pulse_area_bot_phe   = dp.pulse_area_phe./(1+dp.top_bottom_ratio);
dp.xyz_corrected_pulse_area_bot_phe = dp.pulse_area_phe./(1+dp.top_bottom_ratio);
dp.correction_s2_z_dependence       = zeros(max_num_pulses,N, 'single');
dp.correction_s1_z_dependence       = zeros(max_num_pulses,N, 'single');
dp.correction_s1_xyz_dependence     = zeros(max_num_pulses,N, 'single');
dp.correction_s2_xy_dependence      = zeros(max_num_pulses,N, 'single');

% Initialize variable with negative value so that we're sure it goes through the code 
dp.s1_correction_is_3d              = zeros(1,N,'int8');
dp.s1_correction_is_3d(:)           = -99;

%% Finding the S2 Z Correction IQ

[a num_iqs]      = size(lug_iqs_xml.iq); %#ok<*ASGLU>

%Intializing variables
s2_norm_z_iqs    = {};
s2_norm_iqs      = {};
dataset_times    = [];
filename_prefixs = {};
cp_number        = {};

for i=1:num_iqs
    if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   
        if (isfield(lug_iqs_xml.iq(i).correction.fit,'s2_both_norm_z_index').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1         %#ok<*BDSCA>
           s2_norm_z_iqs    = [s2_norm_z_iqs lug_iqs_xml.iq(i).correction.fit.s2_both_norm_z_index]; %#ok<*AGROW>
           s2_norm_iqs      = [s2_norm_iqs lug_iqs_xml.iq(i).correction.fit.s2_both_norm_z];
           dataset_times    = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix) ];
           filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
           cp_number        = [cp_number lug_iqs_xml.iq(i).global.cp_number];
        end
    end  
end
 
current_data_set_time = filename2epoch_framework(filename_prefix ); %#ok<*NASGU>
[s2_norm_z s2_norm]   = InterpS2IQ_Run04(filename_prefix, dataset_times, s2_norm_z_iqs, s2_norm_iqs);
 
%% Finding the S1 Z-dep values


z_dep_both_values    = zeros(0,3);
z_dep_bottom_values  = zeros(0,3);
dataset_times        = [];
filename_prefixs     = {};
cp_number            = {};

for i=1:num_iqs
    if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   
        if (isfield(lug_iqs_xml.iq(i).correction.fit,'s1_both_zdep_quad_fit').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1        
           z_dep_both_values   = vertcat(z_dep_both_values,lug_iqs_xml.iq(i).correction.fit.s1_both_zdep_quad_fit);
           z_dep_bottom_values = vertcat(z_dep_bottom_values,lug_iqs_xml.iq(i).correction.fit.s1_bottom_zdep_quad_fit);
           dataset_times       = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix )];
           filename_prefixs    = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
           cp_number           = [cp_number lug_iqs_xml.iq(i).global.cp_number];

        end
    end  
end
 
[index iq_time_zDep] = NearestIQ(filename_prefix,dataset_times);
z_dep_par_all        = z_dep_both_values(index,:);
z_dep_par_bot        = z_dep_bottom_values(index,:);

%% Finding the S2 xy correction map values


s2_xy_index      = [];
dataset_times    = [];
filename_prefixs = {};
cp_number        = {};

for i=1:num_iqs
    if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   
        if (isfield(lug_iqs_xml.iq(i).correction.fit,'norm_s2_both').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1        
           s2_xy_index      = [s2_xy_index i];
           dataset_times    = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix ) ];
           filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
           cp_number        = [cp_number lug_iqs_xml.iq(i).global.cp_number];             
        end
    end  
end
 
   [index iq_time_s2xyDep] = NearestIQ(filename_prefix,dataset_times);
      
   s2_x_bins     = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.x_bin_center;
   s2_y_bins     = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.y_bin_center;
   s2_map_all    = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.norm_s2_both;
   s2_map_bottom = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.norm_s2_bottom;
   
%% Finding the S1 xy correction map values
   
s1_xy_index      = [];
dataset_times    = [];
filename_prefixs = {};
cp_number        = {};

for i=1:num_iqs
    if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   
        if (isfield(lug_iqs_xml.iq(i).correction.fit,'norm_s1_both').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1       
           s1_xy_index      = [s1_xy_index i];                      
           dataset_times    = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix ) ];
           filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
           cp_number        = [cp_number lug_iqs_xml.iq(i).global.cp_number];              
        end
    end  
end
 
[index iq_time]  = NearestIQ(filename_prefix,dataset_times);
s1_x_bins        = lug_iqs_xml.iq(s1_xy_index(index)).correction.fit.x_bin_center;
s1_y_bins        = lug_iqs_xml.iq(s1_xy_index(index)).correction.fit.y_bin_center;
s1_map_all       = lug_iqs_xml.iq(s1_xy_index(index)).correction.fit.norm_s1_both;
s1_map_bottom    = lug_iqs_xml.iq(s1_xy_index(index)).correction.fit.norm_s1_bottom;

%% Finding the S1 xyz correction map values
   
s1_xyz_index     = [];
dataset_times    = [];
filename_prefixs = {};
cp_number        = {};

for i=1:num_iqs
    if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   
        if (isfield(lug_iqs_xml.iq(i).correction.fit,'norm_s1_both_xyz').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1       
           s1_xyz_index     = [s1_xyz_index i];                      
           dataset_times    = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix )];
           filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
           cp_number        = [cp_number lug_iqs_xml.iq(i).global.cp_number];             
        end
    end  
end
 
[index iq_time_s1xyzDep] = NearestIQ(filename_prefix,dataset_times);
s1_xyz_x_bins            = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.x_bin_center;
s1_xyz_y_bins            = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.y_bin_center;
s1_xyz_z_bins            = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.z_bin_center;
xx                       = size(s1_xyz_x_bins);
yy                       = size(s1_xyz_y_bins);
zz                       = size(s1_xyz_z_bins);
s1_xyz_map_all           = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.norm_s1_both_xyz; s1_xyz_map_all = reshape(s1_xyz_map_all,xx(2),yy(2),zz(2));
s1_xyz_map_bottom        = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.norm_s1_bottom_xyz;   s1_xyz_map_bottom = reshape(s1_xyz_map_bottom,xx(2),yy(2),zz(2));
    
%% Finding the detector center according to KrypCal
   
det_center_index = [];
dataset_times    = [];
filename_prefixs = {};
cp_number        = {};

for i=1:num_iqs
    if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   
        if (isfield(lug_iqs_xml.iq(i).correction.fit,'x_center').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1       
           det_center_index = [det_center_index i];                      
           dataset_times    = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix )];
           filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
           cp_number        = [cp_number lug_iqs_xml.iq(i).global.cp_number];
        end
    end  
end
 
[index iq_time_det_center] = NearestIQ(filename_prefix,dataset_times);
det_x_center               = lug_iqs_xml.iq(det_center_index(index)).correction.fit.x_center;
det_y_center               = lug_iqs_xml.iq(det_center_index(index)).correction.fit.y_center;
det_z_center               = lug_iqs_xml.iq(det_center_index(index)).correction.fit.z_center;
s1_ref_z_ns                = det_z_center.*1000; %conver uSec to ns
    
%% Finding the drift time, x and y associated with S1 from the largest S2 pulse in the pairing

[a b]                       = size(dp.z_drift_samples); % b is number of events in file 
drift                       = dp.z_drift_samples;
drift(find(isnan(drift)))   = 0.0;
s1_drift_ns                 = +(dp.pulse_classification==1);
s1_x_cm                     = +(dp.pulse_classification==1);
s1_y_cm                     = +(dp.pulse_classification==1);
s2_phe                      =  dp.pulse_area_phe.*(dp.pulse_classification==2);
s2_phe(find(isnan(s2_phe))) = 0.0;
s1s                         = (sum(dp.s1s2_pairing==1)>1); % tells you if its an event with at least one s1 and one s2.

for i=1:length(s1s)  % for each event in the file
    if s1s(i)>0 % if the event had an S1
        if s1s(i)==1
            [v r]             = max(s2_phe(:,i));% value and index (r for row) of Max S2
            [c1 r1 v1]        = find(dp.pulse_classification(:,i)==1,1,'first');
            s1_drift_ns(c1,i) = 10.*drift(r,i); %row c1,r and column i
            s1_x_cm(c1,i)     = dp.x_cm(r,i);
            s1_y_cm(c1,i)     = dp.y_cm(r,i);  
        else
            s1_drift_ns(:,i)  = 0;
            s1_x_cm(:,i)      = 0;
            s1_y_cm(:,i)      = 0;
        end  
    end
end

 
fprintf('Done finding S1 depths :) \n');
s1_drift_ns(find((s1_drift_ns==0))) = nan;
s1_x_cm(find((s1_x_cm==0)))         = nan;
s1_y_cm(find((s1_y_cm==0)))         = nan;

%% Calculate Z corrections

% Calculate S1 Z-correction
s1_z_correction                                       = polyval(z_dep_par_all,s1_ref_z_ns./1000)./polyval(z_dep_par_all,s1_drift_ns./1000);
s1_z_correction(find(isnan(s1_z_correction)))         = 1.0;
s1_z_correction_bot                                   = polyval(z_dep_par_bot,s1_ref_z_ns./1000)./polyval(z_dep_par_bot,s1_drift_ns./1000);
s1_z_correction_bot(find(isnan(s1_z_correction_bot))) = 1.0;
s1_z_correction(dp.pulse_classification~=1)           = 1.0;

% Calculate electron lifetime correction (S2 Z-correction)
% Reading the values of electron lifetime from the LUG and interpolating
% between the two either side
s2_z_correction                                       = interp1(s2_norm_z,s2_norm, double(dp.z_drift_samples)/100,'cubic');
s2_z_correction(find(isnan(s2_z_correction)))         = 1.0;
s2_z_correction(dp.pulse_classification~=2)           = 1.0;

%% Calculate XY corrections

%Calculate S1 XY-only corrections
s2x                                                   = dp.x_cm.*(+(dp.pulse_classification==2)); 
s2x(find(s2x==0))                                     = nan;
s2y                                                   = dp.y_cm.*(+(dp.pulse_classification==2)); 
s2y(find(s2y==0))                                     = nan;

s1xy_correction                                       = interp2(s1_x_bins,s1_y_bins,s1_map_all,s1_x_cm,s1_y_cm,'cubic');
s1xy_correction                                       = s1xy_correction.*(+dp.pulse_classification==1);
s1xy_correction(find(s1xy_correction==0))             = 1.0;  
s1xy_correction(find(isnan(s1xy_correction)))         = 1.0;

s1xy_correction_bot                                   = interp2(s1_x_bins,s1_y_bins,s1_map_bottom,s1_x_cm,s1_y_cm,'cubic');
s1xy_correction_bot                                   = s1xy_correction_bot.*(+dp.pulse_classification==1);
s1xy_correction_bot(find(s1xy_correction_bot==0))     = 1.0;  
s1xy_correction_bot(find(isnan(s1xy_correction_bot))) = 1.0;

%Calculate S2 XY corrections
s2xy_correction                                       = interp2(s2_x_bins,s2_y_bins,s2_map_all,dp.x_cm,dp.y_cm,'cubic');
s2xy_correction                                       = s2xy_correction.*(+dp.pulse_classification==2);
s2xy_correction(find(s2xy_correction==0))             = 1.0; 
s2xy_correction(find(isnan(s2xy_correction)))         = 1.0;

s2xy_correction_bot                                   = interp2(s2_x_bins,s2_y_bins,s2_map_bottom,dp.x_cm,dp.y_cm,'cubic');
s2xy_correction_bot                                   = s2xy_correction_bot.*(+dp.pulse_classification==2);
s2xy_correction_bot(find(s2xy_correction_bot==0))     = 1.0;  
s2xy_correction_bot(find(isnan(s2xy_correction_bot))) = 1.0;

%% Calculate XYZ corrections

% ALWAYS DOING 1D+2D NOW. (3/21/2016)
s1xyz_correction          = s1xy_correction.*s1_z_correction;
s1xyz_correction_bot      = s1xy_correction_bot.*s1_z_correction_bot;
dp.s1_correction_is_3d(:) = 0;

%% add RQs of correction factors

dp.correction_s2_z_dependence                 = single(s2_z_correction);
dp.correction_s1_z_dependence                 = single(s1_z_correction);
dp.correction_s1_z_dependence_bot             = single(s1_z_correction_bot);
dp.correction_s1_xyz_dependence               = single(s1xyz_correction); %this is either 3D or 1D+2D
dp.correction_s1_xyz_dependence_bot           = single(s1xyz_correction_bot); %this is either 3D or 1D+2D
dp.correction_s1_xy_dependence                = single(s1xy_correction);
dp.correction_s1_xy_dependence_bot            = single(s1xy_correction_bot);
dp.correction_s2_xy_dependence                = single(s2xy_correction);
dp.correction_s2_xy_dependence_bot            = single(s2xy_correction_bot);
dp.admin.corrections.s1_z_dependence.iq_time  = iq_time_zDep;
dp.admin.corrections.s1_z_dependence.all      = z_dep_par_all;
dp.admin.corrections.s1_z_dependence.bottom   = z_dep_par_bot;
dp.admin.corrections.s1_xy_dependence.iq_time = iq_time;
dp.admin.corrections.s2_xy_dependence.iq_time = iq_time_s2xyDep;

%% Applying corrections 

%Apply Z corrections
dp.z_corrected_pulse_area_all_phe = dp.z_corrected_pulse_area_all_phe.*s2_z_correction.*s1_z_correction; 
dp.z_corrected_pulse_area_bot_phe = dp.z_corrected_pulse_area_bot_phe.*s2_z_correction.*s1_z_correction_bot; 


%Apply XYZ corrections
dp.xyz_corrected_pulse_area_all_phe = dp.xyz_corrected_pulse_area_all_phe.*s2_z_correction.*s2xy_correction.*s1xyz_correction;
dp.xyz_corrected_pulse_area_bot_phe = dp.xyz_corrected_pulse_area_bot_phe.*s2_z_correction.*s2xy_correction_bot.*s1xyz_correction_bot;

%% Write Output File
status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, dp.event_number, livetime); % Should add livetime input at the end, remove event_number as settings field

